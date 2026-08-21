"""The core structure representation: a struct-of-arrays molecular structure.

Two conventions, stated once and enforced everywhere
----------------------------------------------------
1. **Indices are 0-based positional.** ``residue_index`` counts residues from the start
   of the structure. PDB numbering (``residue_number``) is *data*: it is preserved,
   written back out, and never used to index or key anything. Crossing these two spaces
   was the single largest bug source in the pre-rewrite code -- it made region selection
   silently match nothing on any structure not numbered from 1, which is to say on every
   crystal and EM structure.

2. **Spans are half-open ``[start, stop)``.** Always. Adjacent regions do not share a
   boundary residue. The pre-rewrite code used inclusive bounds in which an IDR's first
   and last residues *were* the terminal residues of the flanking folded domains,
   because those residues serve as geometric anchors -- and then resolved the resulting
   overlap with three different conventions inside a single function. Here, anchors are
   named explicitly (:attr:`Span.n_anchor`, :attr:`Span.c_anchor`) and are not part of
   the span they anchor.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field, replace
from enum import Enum
from itertools import pairwise

import numpy as np
from scipy.spatial import cKDTree

from .constants import (
    ATOMIC_MASSES,
    CA_CLASH_DISTANCE,
    CLASH_EXCLUDE_WITHIN_RESIDUES,
    THREE_TO_ONE,
)
from .exceptions import EmptyStructureError, GeometryError, InvalidRegionError

__all__ = [
    "Chain",
    "Domain",
    "DomainKind",
    "Span",
    "Structure",
]


class DomainKind(str, Enum):
    """What kind of region a :class:`Domain` is.

    Deliberately only two values. The pre-rewrite code had a third,
    "folded-with-loops", which duplicated every code path that handled folded domains.
    A folded domain with rebuildable loops is just a :attr:`FOLDED` domain whose
    :attr:`Domain.loops` is non-empty.
    """

    FOLDED = "folded"
    IDR = "idr"


@dataclass(frozen=True, slots=True)
class Span:
    """A half-open range of residue indices, ``[start, stop)``.

    Attributes
    ----------
    start
        First residue index in the span, 0-based positional, inclusive.
    stop
        One past the last residue index. Exclusive.
    n_anchor
        Residue index of the fixed residue immediately N-terminal to this span, or
        ``None`` if the span starts at a chain terminus. This residue is *not* part of
        the span; it is the geometric anchor a rebuilt region must connect to.
    c_anchor
        Residue index of the fixed residue immediately C-terminal to this span, or
        ``None`` if the span ends at a chain terminus.
    """

    start: int
    stop: int
    n_anchor: int | None = None
    c_anchor: int | None = None

    def __post_init__(self) -> None:
        if self.stop < self.start:
            raise InvalidRegionError(f"Span stop ({self.stop}) precedes start ({self.start}).")
        if self.start < 0:
            raise InvalidRegionError(f"Span start must be non-negative, got {self.start}.")

    def __len__(self) -> int:
        return self.stop - self.start

    def __iter__(self) -> Iterator[int]:
        return iter(range(self.start, self.stop))

    def __contains__(self, residue_index: int) -> bool:
        return self.start <= residue_index < self.stop

    @property
    def slice(self) -> slice:
        """Return this span as a :class:`slice`, for indexing residue-level arrays."""
        return slice(self.start, self.stop)

    @property
    def is_terminal(self) -> bool:
        """True if this span is missing an anchor on either side.

        A terminal span (a tail) is free at one end, so it is built outward toward a
        sampled target rather than closed onto two fixed points.
        """
        return self.n_anchor is None or self.c_anchor is None

    def overlaps(self, other: Span) -> bool:
        """Return True if this span shares any residue with ``other``."""
        return self.start < other.stop and other.start < self.stop


@dataclass(slots=True)
class Domain:
    """A view onto a contiguous region of a :class:`Structure`.

    A ``Domain`` owns no coordinates. Its :attr:`xyz` is a zero-copy view into the
    parent structure's array, so it is always current by construction.
    """

    structure: Structure
    span: Span
    kind: DomainKind
    #: Rebuildable loops inside a folded domain, as spans of residue indices. Always
    #: empty for an IDR (an IDR is rebuilt in its entirety, so a "loop" within it is
    #: not a distinct concept).
    loops: tuple[Span, ...] = ()
    #: Whether DODO **generated** this domain's coordinates. Set only for an IDR whose build
    #: succeeded. Never set for a folded domain, which is translated and rotated but never
    #: regenerated, and never set for a region whose build failed.
    #:
    #: This is the provenance question, and it is NOT the same as :attr:`placed`. Conflating
    #: them made the validators attribute AlphaFold's own broken geometry to DODO -- see
    #: :meth:`generated_spans`.
    rebuilt: bool = False
    #: Whether this domain's coordinates are **final**, so it must be avoided as an obstacle.
    #:
    #: Deliberately a different question from :attr:`rebuilt`, and all three combinations are
    #: real:
    #:
    #: * a folded domain is placed once repositioned, and never generated;
    #: * an IDR that built successfully is both placed and generated;
    #: * an IDR whose build FAILED is placed but not generated -- it keeps its input
    #:   coordinates, and nothing will move them again.
    #:
    #: That last case is why this field exists. While ``rebuilt`` served both purposes, a failed
    #: region was excluded from the obstacle set for the rest of the run, so a later region
    #: could be built straight through it.
    placed: bool = False
    #: Indices into :attr:`loops` of the loops that were **successfully** rebuilt. A loop whose
    #: build failed keeps its input coordinates, so it is absent here.
    rebuilt_loops: set[int] = field(default_factory=set)
    #: Optional label carried through to output, e.g. for beta-factor annotation.
    label: str | None = None

    def __repr__(self) -> str:
        loops = f", {len(self.loops)} loops" if self.loops else ""
        return (
            f"Domain({self.kind.value}, residues "
            f"{self.span.start}-{self.span.stop}{loops}"
            f"{', rebuilt' if self.rebuilt else ''})"
        )

    def __len__(self) -> int:
        return len(self.span)

    def generated_spans(self) -> tuple[Span, ...]:
        """Spans whose coordinates DODO actually generated, as opposed to merely moved.

        The distinction matters and getting it wrong is not cosmetic. Two things are *not*
        DODO-generated even though they live inside a domain marked :attr:`rebuilt`:

        * **A folded domain.** It is translated and rotated as a rigid body and its atoms are
          never regenerated, so every bond inside it is exactly the input's.
        * **A region whose build failed.** Its input coordinates are deliberately left in
          place, so it is input geometry too.

        Both used to be reported as DODO's work. 
        """
        spans: list[Span] = []
        if self.kind is DomainKind.IDR and self.rebuilt:
            spans.append(self.span)
        spans.extend(self.loops[index] for index in sorted(self.rebuilt_loops))
        return tuple(spans)

    @property
    def atom_slice(self) -> slice:
        """Slice selecting this domain's atoms in the parent's atom-level arrays."""
        return self.structure.atom_slice_for_residues(self.span.start, self.span.stop)

    @property
    def xyz(self) -> np.ndarray:
        """Return atomic coordinates as a zero-copy ``(n_atoms, 3)`` view.

        Writing through this view mutates the parent structure, which is intended --
        it is how transforms are applied without copying.
        """
        return self.structure.xyz[self.atom_slice]

    @property
    def ca_xyz(self) -> np.ndarray:
        """Alpha-carbon coordinates, one row per residue, as an ``(n_residues, 3)`` array.

        A copy, not a view: alpha carbons are not contiguous in the atom array.
        """
        return self.structure.ca_xyz[self.span.slice]

    @property
    def sequence(self) -> str:
        """One-letter sequence of this domain."""
        return self.structure.sequence[self.span.slice]

    @property
    def n_terminal_ca(self) -> np.ndarray | None:
        """Coordinates of the alpha carbon anchoring this domain on its N-terminal side.

        ``None`` at a chain start. This reads the *anchor* residue, which lies outside
        the span. The pre-rewrite code derived the equivalent value from "the first
        atom of the domain", which broke as soon as rebuilt residues were appended out
        of order -- it returned an interior loop residue and grafted the next region
        onto the wrong atom.
        """
        if self.span.n_anchor is None:
            return None
        anchor: np.ndarray = self.structure.ca_xyz[self.span.n_anchor]
        return anchor

    @property
    def c_terminal_ca(self) -> np.ndarray | None:
        """Coordinates of the alpha carbon anchoring this domain C-terminally."""
        if self.span.c_anchor is None:
            return None
        anchor: np.ndarray = self.structure.ca_xyz[self.span.c_anchor]
        return anchor

    def centroid(self, mass_weighted: bool = False) -> np.ndarray:
        """Geometric centroid, or mass-weighted centre of mass.

        The pre-rewrite code had two methods both named as if they computed a centre of
        mass, one mass-weighted and one not, and rotated domains about whichever it
        happened to call. The distinction is explicit here.
        """
        coords = self.xyz
        if coords.shape[0] == 0:
            raise GeometryError(f"Cannot take the centroid of an empty domain: {self!r}")
        if not mass_weighted:
            centroid: np.ndarray = coords.mean(axis=0)
            return centroid
        masses = self.structure.atom_masses[self.atom_slice]
        weighted: np.ndarray = (masses[:, None] * coords).sum(axis=0) / masses.sum()
        return weighted

    def translate(self, vector: Sequence[float] | np.ndarray) -> None:
        """Translate this domain in place."""
        self.structure.xyz[self.atom_slice] += np.asarray(vector, dtype=np.float64)

    def rotate(self, rotation: np.ndarray, about: np.ndarray | None = None) -> None:
        """Rotate this domain in place.

        Parameters
        ----------
        rotation
            A ``(3, 3)`` rotation matrix.
        about
            Point to rotate about. Defaults to this domain's centroid, which is what
            "rotate the domain" almost always means. The pre-rewrite ``Chain.rotate``
            and ``Complex.rotate`` rotated about the world origin instead, silently
            translating as well as rotating.
        """
        rotation = np.asarray(rotation, dtype=np.float64)
        if rotation.shape != (3, 3):
            raise GeometryError(f"Rotation must be a (3, 3) matrix, got {rotation.shape}.")
        centre = self.centroid() if about is None else np.asarray(about, dtype=np.float64)
        sl = self.atom_slice
        self.structure.xyz[sl] = (self.structure.xyz[sl] - centre) @ rotation.T + centre


@dataclass(slots=True)
class Chain:
    """A view onto one polymer chain of a :class:`Structure`."""

    structure: Structure
    span: Span
    chain_id: str
    domains: list[Domain] = field(default_factory=list)
    #: UniProt accession, when known from a DBREF record or an mmCIF struct_ref.
    uniprot_id: str | None = None
    #: The full construct sequence as deposited (SEQRES / entity_poly), when known.
    #: Compared against the observed sequence to locate unmodelled residues. This is
    #: deliberately the *deposited* sequence rather than the UniProt canonical isoform:
    #: using the canonical isoform invents phantom terminal residues for any tagged or
    #: truncated construct.
    full_sequence: str | None = None

    def __repr__(self) -> str:
        return (
            f"Chain({self.chain_id!r}, residues {self.span.start}-{self.span.stop}, "
            f"{len(self.domains)} domains)"
        )

    def __len__(self) -> int:
        return len(self.span)

    @property
    def sequence(self) -> str:
        """Observed one-letter sequence of this chain."""
        return self.structure.sequence[self.span.slice]

    @property
    def ca_xyz(self) -> np.ndarray:
        """Alpha-carbon coordinates for this chain, ``(n_residues, 3)``."""
        return self.structure.ca_xyz[self.span.slice]

    def domain_at(self, residue_index: int) -> Domain | None:
        """Return the domain containing ``residue_index``, or ``None`` if unassigned."""
        for domain in self.domains:
            if residue_index in domain.span:
                return domain
        return None

    def validate_domains(self) -> None:
        """Check that this chain's domains tile its residues without gaps or overlaps.

        Raises
        ------
        InvalidRegionError
            If domains overlap, extend outside the chain, or leave residues
            unassigned. The pre-rewrite manual-override path accepted all three
            silently and failed much later with an unrelated error.
        """
        if not self.domains:
            raise InvalidRegionError(f"Chain {self.chain_id!r} has no domains assigned.")

        # Report residues the way the caller wrote them and the way `dodo regions` prints
        # them. These messages used to name raw positional indices, which appear neither in
        # the user's input nor anywhere in their file.
        #
        # Every labeller here is bounds-checked, because several of these branches exist
        # precisely to report a span that runs off the end of the structure -- an error path
        # that raises IndexError while formatting its own message is worse than the bug it
        # was trying to describe.
        n_residues = self.structure.n_residues

        def label(index: int) -> str:
            if 0 <= index < n_residues:
                return self.structure.residue_id(index)
            return f"index {index}"

        def span_label(span: Span) -> str:
            if 0 <= span.start < n_residues and 0 < span.stop <= n_residues:
                return f"{label(span.start)}-{label(span.stop - 1)}"
            return (
                f"positional {span.start}-{span.stop}, outside this structure's "
                f"{n_residues} residues"
            )

        ordered = sorted(self.domains, key=lambda d: d.span.start)
        if ordered[0].span.start != self.span.start:
            raise InvalidRegionError(
                f"Chain {self.chain_id!r} domains start at residue "
                f"{label(ordered[0].span.start)} but the chain starts at "
                f"{label(self.span.start)}."
            )
        if ordered[-1].span.stop != self.span.stop:
            raise InvalidRegionError(
                f"Chain {self.chain_id!r} domains end at residue "
                f"{label(ordered[-1].span.stop - 1)} but the chain ends at "
                f"{label(self.span.stop - 1)}."
            )
        for previous, nxt in pairwise(ordered):
            if previous.span.stop < nxt.span.start:
                raise InvalidRegionError(
                    f"Chain {self.chain_id!r} has unassigned residues "
                    f"{label(previous.span.stop)}-{label(nxt.span.start - 1)}."
                )
            if previous.span.stop > nxt.span.start:
                raise InvalidRegionError(
                    f"Chain {self.chain_id!r} has overlapping domains: "
                    f"{span_label(previous.span)} and {span_label(nxt.span)}."
                )
        for domain in self.domains:
            for loop in domain.loops:
                if loop.start < domain.span.start or loop.stop > domain.span.stop:
                    raise InvalidRegionError(
                        f"Loop {span_label(loop)} lies outside its domain "
                        f"{span_label(domain.span)} of chain {self.chain_id!r}."
                    )


@dataclass(slots=True)
class Structure:
    """A molecular structure held as parallel numpy arrays.

    One instance owns exactly one coordinate array. All views into it -- chains,
    domains, loops -- are index ranges, so no coordinate is ever duplicated and no
    view can go stale.

    Invariants, checked by :meth:`validate`
    ---------------------------------------
    * Atoms are ordered by residue: all atoms of residue ``i`` precede all atoms of
      residue ``i + 1``. This is what makes a residue range map to an atom *slice*
      rather than a fancy index, and therefore what makes views zero-copy.
    * ``residue_atom_offsets`` has length ``n_residues + 1`` and is non-decreasing;
      residue ``i``'s atoms are ``[offsets[i], offsets[i + 1])``.
    * Every atom-level array has length ``n_atoms``; every residue-level array has
      length ``n_residues``.
    """

    # --- atom-level arrays, all length n_atoms -----------------------------------
    #: The coordinates. ``(n_atoms, 3)`` float64. The single source of truth.
    xyz: np.ndarray
    #: PDB atom names, e.g. ``"CA"``, ``"OD1"``. ``(n_atoms,)`` unicode.
    atom_name: np.ndarray
    #: Element symbols, uppercase. ``(n_atoms,)`` unicode.
    element: np.ndarray
    #: 0-based positional residue index for each atom. ``(n_atoms,)`` int64.
    residue_index: np.ndarray

    # --- residue-level arrays, all length n_residues ------------------------------
    #: Three-letter residue names. ``(n_residues,)`` unicode.
    residue_name: np.ndarray
    #: PDB residue numbers. Data, never an index. ``(n_residues,)`` int64.
    residue_number: np.ndarray
    #: PDB insertion codes, ``""`` when absent. Preserving these is what keeps
    #: residues 10 and 10A distinct; the pre-rewrite reader merged them.
    insertion_code: np.ndarray
    #: B-factor per residue, taken from the alpha carbon. For AlphaFold models this is
    #: pLDDT, which is a directly usable disorder signal.
    b_factor: np.ndarray
    #: Occupancy per residue, from the alpha carbon.
    occupancy: np.ndarray
    #: 0-based positional chain index for each residue. ``(n_residues,)`` int64.
    chain_index: np.ndarray
    #: CSR-style offsets into the atom arrays. ``(n_residues + 1,)`` int64.
    residue_atom_offsets: np.ndarray

    # --- structure-level metadata --------------------------------------------------
    chains: list[Chain] = field(default_factory=list)
    #: Where this structure came from, for error messages and provenance.
    source: str | None = None
    #: Free-form provenance notes, e.g. which records were skipped while parsing.
    notes: list[str] = field(default_factory=list)

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __post_init__(self) -> None:
        self.xyz = np.ascontiguousarray(self.xyz, dtype=np.float64)
        if self.xyz.ndim != 2 or self.xyz.shape[1] != 3:
            raise GeometryError(f"xyz must have shape (n_atoms, 3), got {self.xyz.shape}.")

    def __repr__(self) -> str:
        return (
            f"Structure({self.n_atoms} atoms, {self.n_residues} residues, "
            f"{len(self.chains)} chains" + (f", from {self.source!r}" if self.source else "") + ")"
        )

    # ------------------------------------------------------------------
    # Shape
    # ------------------------------------------------------------------

    @property
    def n_atoms(self) -> int:
        """Number of atoms in this structure."""
        return int(self.xyz.shape[0])

    @property
    def n_residues(self) -> int:
        """Number of residues in this structure."""
        return int(self.residue_name.shape[0])

    # ------------------------------------------------------------------
    # Derived views and lookups
    # ------------------------------------------------------------------

    def atom_slice_for_residues(self, start: int, stop: int) -> slice:
        """Atom-array slice covering residues ``[start, stop)``.

        O(1), and returns a real :class:`slice` so downstream indexing produces views
        rather than copies. This is the payoff of keeping atoms residue-ordered.
        """
        if not (0 <= start <= stop <= self.n_residues):
            raise InvalidRegionError(
                f"Residue range [{start}, {stop}) is outside this structure, which has "
                f"{self.n_residues} residues."
            )
        return slice(int(self.residue_atom_offsets[start]), int(self.residue_atom_offsets[stop]))

    @property
    def ca_indices(self) -> np.ndarray:
        """Atom index of each residue's alpha carbon, ``(n_residues,)`` int64.

        Raises
        ------
        EmptyStructureError
            If any residue has no alpha carbon. Every downstream geometric operation
            assumes one exists, so failing here with a clear message beats failing
            later with an index error.
        """
        is_ca = self.atom_name == "CA"
        idx = np.full(self.n_residues, -1, dtype=np.int64)
        idx[self.residue_index[is_ca]] = np.flatnonzero(is_ca)
        missing = np.flatnonzero(idx < 0)
        if missing.size:
            preview = ", ".join(
                f"{self.residue_name[i]}{self.residue_number[i]}" for i in missing[:5]
            )
            more = f" (and {missing.size - 5} more)" if missing.size > 5 else ""
            raise EmptyStructureError(
                f"{missing.size} residue(s) have no CA atom: {preview}{more}. "
                f"DODO needs an alpha carbon for every residue."
            )
        return idx

    @property
    def ca_xyz(self) -> np.ndarray:
        """Alpha-carbon coordinates, ``(n_residues, 3)``.

        A copy, since alpha carbons are not contiguous. Cheap relative to the geometry
        that consumes it, and callers that need to *write* coordinates go through
        :meth:`set_ca_xyz` so the intent is explicit.
        """
        coords: np.ndarray = self.xyz[self.ca_indices]
        return coords

    def set_ca_xyz(self, residue_indices: np.ndarray | Sequence[int], coords: np.ndarray) -> None:
        """Write alpha-carbon coordinates for the given residues."""
        residue_indices = np.asarray(residue_indices, dtype=np.int64)
        coords = np.asarray(coords, dtype=np.float64)
        if coords.shape != (residue_indices.size, 3):
            raise GeometryError(
                f"Expected coords of shape ({residue_indices.size}, 3), got {coords.shape}."
            )
        self.xyz[self.ca_indices[residue_indices]] = coords

    @property
    def sequence(self) -> str:
        """One-letter sequence of the whole structure, in residue order."""
        return "".join(THREE_TO_ONE.get(str(name), "X") for name in self.residue_name)

    @property
    def atom_masses(self) -> np.ndarray:
        """Atomic mass of each atom, ``(n_atoms,)`` float64.

        Unrecognized elements fall back to carbon rather than raising: an exotic
        heteroatom should not stop a visualization, and the error in a centre of mass
        from one mislabelled atom is negligible.
        """
        return np.array(
            [ATOMIC_MASSES.get(str(e).upper(), ATOMIC_MASSES["C"]) for e in self.element],
            dtype=np.float64,
        )

    def residue_id(self, residue_index: int) -> str:
        """Residue number with its insertion code, e.g. ``"142"`` or ``"142A"``.

        The file's own numbering, for output a user will compare against the input
        structure or type into a viewer. Distinct from :meth:`residue_label`, which
        prepends chain and residue name and is meant for error messages.
        """
        return f"{int(self.residue_number[residue_index])}{self.insertion_code[residue_index]!s}"

    def residue_label(self, residue_index: int) -> str:
        """Human-readable residue label, e.g. ``"A:GLU142"``. For error messages."""
        chain = self.chains[int(self.chain_index[residue_index])].chain_id if self.chains else "?"
        icode = str(self.insertion_code[residue_index])
        return (
            f"{chain}:{self.residue_name[residue_index]}"
            f"{int(self.residue_number[residue_index])}{icode}"
        )

    # ------------------------------------------------------------------
    # Domains
    # ------------------------------------------------------------------

    @property
    def domains(self) -> list[Domain]:
        """All domains across all chains, in residue order."""
        return [d for chain in self.chains for d in chain.domains]

    def idrs(self) -> list[Domain]:
        """All IDR domains, in residue order."""
        return [d for d in self.domains if d.kind is DomainKind.IDR]

    def folded_domains(self) -> list[Domain]:
        """All folded domains, in residue order."""
        return [d for d in self.domains if d.kind is DomainKind.FOLDED]

    # ------------------------------------------------------------------
    # Geometry
    # ------------------------------------------------------------------

    def kdtree(self, atom_mask: np.ndarray | None = None) -> cKDTree:
        """Build a KD-tree over this structure's atoms.

        Deliberately *not* cached. The pre-rewrite code cached spatial indexes inside
        the structure objects and got the invalidation wrong in three separate ways --
        including returning a tree built from a different subset than the caller asked
        for. Callers that need a tree across many queries should hold it themselves;
        they know when they mutate coordinates and the structure does not.

        Parameters
        ----------
        atom_mask
            Boolean mask over atoms. Defaults to all atoms.
        """
        coords = self.xyz if atom_mask is None else self.xyz[atom_mask]
        if coords.shape[0] == 0:
            raise GeometryError(
                "Cannot build a spatial index over zero atoms. Callers should check "
                "for an empty obstacle set instead of relying on this."
            )
        return cKDTree(coords)

    def placed_atom_mask(self) -> np.ndarray:
        """Boolean atom mask selecting atoms whose coordinates are final.

        The obstacle set for clash checking during a build. Geometry that has yet to be placed
        should not veto a candidate conformation -- but geometry that is final must, and that
        includes a region whose build FAILED. Its input coordinates are staying, so a later
        region must avoid it. Keying this on :attr:`Domain.rebuilt` instead left failed regions
        invisible for the rest of the run.
        """
        mask = np.zeros(self.n_atoms, dtype=bool)
        for domain in self.domains:
            if domain.placed:
                mask[domain.atom_slice] = True
        return mask

    def clash_mask(
        self,
        query_xyz: np.ndarray,
        *,
        query_residue_index: np.ndarray | Sequence[int] | None = None,
        cutoff: float = CA_CLASH_DISTANCE,
        obstacle_atom_mask: np.ndarray | None = None,
        exclude_within: int = CLASH_EXCLUDE_WITHIN_RESIDUES,
    ) -> np.ndarray:
        """Which query points clash with this structure.

        Parameters
        ----------
        query_xyz
            ``(n_query, 3)`` candidate positions.
        query_residue_index
            Residue index each query point belongs to. When given together with
            ``exclude_within``, contacts to residues within that separation are not
            counted as clashes -- adjacent residues are covalently bonded at 3.8 A and
            next-nearest are angle-constrained to 5-7.5 A, so counting either as a
            clash makes every valid backbone look like a collision. The pre-rewrite
            whole-structure clash check had no such exclusion and duly reported every
            peptide bond.
        cutoff
            Clash distance in Angstroms.
        obstacle_atom_mask
            Which of this structure's atoms count as obstacles. Defaults to all.
        exclude_within
            Residue separation exempt from clash checking.

        Returns
        -------
        np.ndarray
            Boolean mask over the query points, ``True`` where the point clashes.
        """
        query_xyz = np.atleast_2d(np.asarray(query_xyz, dtype=np.float64))
        if query_xyz.shape[1] != 3:
            raise GeometryError(f"query_xyz must be (n, 3), got {query_xyz.shape}.")

        obstacle_indices = (
            np.arange(self.n_atoms)
            if obstacle_atom_mask is None
            else np.flatnonzero(obstacle_atom_mask)
        )
        # An empty obstacle set means nothing can clash. This is the normal state at
        # the very start of a rebuild, and the pre-rewrite equivalent raised a
        # ValueError from scipy here instead of returning the obvious answer.
        if obstacle_indices.size == 0:
            return np.zeros(query_xyz.shape[0], dtype=bool)

        # Through kdtree(), not a bare cKDTree, so the one spatial-index entry point in this
        # class actually has a caller. It was public API with zero call sites and zero coverage,
        # which is how a method quietly stops matching what everything else does.
        mask = np.zeros(self.n_atoms, dtype=bool)
        mask[obstacle_indices] = True
        neighbours = self.kdtree(mask).query_ball_point(query_xyz, cutoff)

        clashing = np.zeros(query_xyz.shape[0], dtype=bool)
        if query_residue_index is None or exclude_within <= 0:
            for i, hits in enumerate(neighbours):
                clashing[i] = len(hits) > 0
            return clashing

        query_residue_index = np.asarray(query_residue_index, dtype=np.int64)
        if query_residue_index.size != query_xyz.shape[0]:
            raise GeometryError(
                f"query_residue_index has {query_residue_index.size} entries but there "
                f"are {query_xyz.shape[0]} query points."
            )
        obstacle_residues = self.residue_index[obstacle_indices]
        for i, hits in enumerate(neighbours):
            if not hits:
                continue
            separation = np.abs(obstacle_residues[hits] - query_residue_index[i])
            clashing[i] = bool(np.any(separation > exclude_within))
        return clashing

    def radius_of_gyration(self, span: Span | None = None, mass_weighted: bool = False) -> float:
        """Radius of gyration in Angstroms, over the CA trace of ``span``.

        Computed on alpha carbons rather than all atoms, because that is the quantity
        ALBATROSS predicts and therefore the one worth comparing against.
        """
        coords = self.ca_xyz if span is None else self.ca_xyz[span.slice]
        if coords.shape[0] < 2:
            raise GeometryError("Radius of gyration needs at least 2 residues.")
        if mass_weighted:
            raise NotImplementedError(
                "Mass-weighted Rg over a CA trace is not meaningful; use all-atom "
                "coordinates if you need it."
            )
        centred = coords - coords.mean(axis=0)
        return float(np.sqrt((centred**2).sum() / coords.shape[0]))

    def end_to_end_distance(self, span: Span | None = None) -> float:
        """CA-to-CA end-to-end distance in Angstroms over ``span``."""
        coords = self.ca_xyz if span is None else self.ca_xyz[span.slice]
        if coords.shape[0] < 2:
            raise GeometryError("End-to-end distance needs at least 2 residues.")
        return float(np.linalg.norm(coords[-1] - coords[0]))

    # ------------------------------------------------------------------
    # Whole-structure transforms
    # ------------------------------------------------------------------

    def translate(self, vector: Sequence[float] | np.ndarray) -> None:
        """Translate every atom in place."""
        self.xyz += np.asarray(vector, dtype=np.float64)

    def rotate(self, rotation: np.ndarray, about: np.ndarray | None = None) -> None:
        """Rotate every atom in place, about ``about`` (default: the centroid)."""
        rotation = np.asarray(rotation, dtype=np.float64)
        if rotation.shape != (3, 3):
            raise GeometryError(f"Rotation must be (3, 3), got {rotation.shape}.")
        centre = self.xyz.mean(axis=0) if about is None else np.asarray(about, dtype=np.float64)
        self.xyz[:] = (self.xyz - centre) @ rotation.T + centre

    # ------------------------------------------------------------------
    # Copying and validation
    # ------------------------------------------------------------------

    def select_atoms(self, keep: np.ndarray) -> Structure:
        """Return a copy containing only the atoms selected by a boolean mask.

        Residues and chains are preserved -- this drops ATOMS, not residues -- so every
        :class:`Domain` and :class:`Chain` span stays valid and is rebound to the new structure.
        A residue that loses every atom would break that invariant, so it is rejected rather
        than silently collapsing the residue numbering.


        Parameters
        ----------
        keep
            Boolean mask over atoms, length ``n_atoms``.

        Returns
        -------
        Structure
            A new structure. The original is unchanged.

        Raises
        ------
        GeometryError
            If the mask is the wrong length, or would leave a residue with no atoms.
        """
        keep = np.asarray(keep, dtype=bool)
        if keep.shape != (self.n_atoms,):
            raise GeometryError(
                f"keep must be a boolean mask of length n_atoms ({self.n_atoms}), got "
                f"shape {keep.shape}."
            )
        if not keep.any():
            raise GeometryError("keep would remove every atom.")

        surviving_per_residue = np.bincount(self.residue_index[keep], minlength=self.n_residues)
        empty = np.flatnonzero(surviving_per_residue == 0)
        if empty.size:
            preview = ", ".join(self.residue_label(int(i)) for i in empty[:5])
            more = f" (and {empty.size - 5} more)" if empty.size > 5 else ""
            raise GeometryError(
                f"keep would leave {empty.size} residue(s) with no atoms at all: {preview}"
                f"{more}. Dropping a residue entirely would invalidate every Domain and Chain "
                f"span, so select_atoms only removes atoms from residues that keep at least one."
            )

        offsets = np.empty(self.n_residues + 1, dtype=np.int64)
        offsets[0] = 0
        np.cumsum(surviving_per_residue, out=offsets[1:])

        new = Structure(
            xyz=self.xyz[keep].copy(),
            atom_name=self.atom_name[keep].copy(),
            element=self.element[keep].copy(),
            residue_index=self.residue_index[keep].copy(),
            residue_name=self.residue_name.copy(),
            residue_number=self.residue_number.copy(),
            insertion_code=self.insertion_code.copy(),
            b_factor=self.b_factor.copy(),
            occupancy=self.occupancy.copy(),
            chain_index=self.chain_index.copy(),
            residue_atom_offsets=offsets,
            source=self.source,
            notes=list(self.notes),
        )
        for chain in self.chains:
            new_chain = Chain(
                structure=new,
                span=chain.span,
                chain_id=chain.chain_id,
                uniprot_id=chain.uniprot_id,
                full_sequence=chain.full_sequence,
            )
            new_chain.domains = [replace(domain, structure=new) for domain in chain.domains]
            new.chains.append(new_chain)

        new.validate()
        return new

    def copy(self) -> Structure:
        """Deep copy, with domain and chain views rebound to the new structure.

        Used to generate independent conformers from one input without re-parsing.
        """
        new = Structure(
            xyz=self.xyz.copy(),
            atom_name=self.atom_name.copy(),
            element=self.element.copy(),
            residue_index=self.residue_index.copy(),
            residue_name=self.residue_name.copy(),
            residue_number=self.residue_number.copy(),
            insertion_code=self.insertion_code.copy(),
            b_factor=self.b_factor.copy(),
            occupancy=self.occupancy.copy(),
            chain_index=self.chain_index.copy(),
            residue_atom_offsets=self.residue_atom_offsets.copy(),
            source=self.source,
            notes=list(self.notes),
        )
        for chain in self.chains:
            new_chain = Chain(
                structure=new,
                span=chain.span,
                chain_id=chain.chain_id,
                uniprot_id=chain.uniprot_id,
                full_sequence=chain.full_sequence,
            )
            new_chain.domains = [replace(domain, structure=new) for domain in chain.domains]
            new.chains.append(new_chain)
        return new

    def validate(self) -> None:
        """Assert every documented invariant. Cheap enough to call after parsing.

        Raises
        ------
        EmptyStructureError
            If the structure has no atoms or residues.
        GeometryError
            If any array length, ordering or offset invariant is violated.
        """
        if self.n_atoms == 0:
            raise EmptyStructureError(f"Structure has no atoms (source: {self.source!r}).")
        if self.n_residues == 0:
            raise EmptyStructureError(f"Structure has no residues (source: {self.source!r}).")

        atom_arrays = {
            "atom_name": self.atom_name,
            "element": self.element,
            "residue_index": self.residue_index,
        }
        for name, array in atom_arrays.items():
            if array.shape[0] != self.n_atoms:
                raise GeometryError(
                    f"{name} has {array.shape[0]} entries but there are {self.n_atoms} atoms."
                )

        residue_arrays = {
            "residue_name": self.residue_name,
            "residue_number": self.residue_number,
            "insertion_code": self.insertion_code,
            "b_factor": self.b_factor,
            "occupancy": self.occupancy,
            "chain_index": self.chain_index,
        }
        for name, array in residue_arrays.items():
            if array.shape[0] != self.n_residues:
                raise GeometryError(
                    f"{name} has {array.shape[0]} entries but there are {self.n_residues} residues."
                )

        if self.residue_atom_offsets.shape[0] != self.n_residues + 1:
            raise GeometryError(
                f"residue_atom_offsets has {self.residue_atom_offsets.shape[0]} entries; "
                f"expected n_residues + 1 = {self.n_residues + 1}."
            )
        if self.residue_atom_offsets[0] != 0 or self.residue_atom_offsets[-1] != self.n_atoms:
            raise GeometryError(
                f"residue_atom_offsets must run from 0 to n_atoms ({self.n_atoms}), got "
                f"{self.residue_atom_offsets[0]} to {self.residue_atom_offsets[-1]}."
            )
        if np.any(np.diff(self.residue_atom_offsets) < 0):
            raise GeometryError("residue_atom_offsets must be non-decreasing.")

        # The ordering invariant that makes views zero-copy.
        if np.any(np.diff(self.residue_index) < 0):
            raise GeometryError(
                "Atoms are not ordered by residue. A residue range must map to a "
                "contiguous atom slice for domain views to be zero-copy."
            )

        # Offsets must agree with residue_index, or slicing silently returns the
        # wrong atoms -- the kind of mismatch that produces plausible wrong output.
        expected = np.searchsorted(self.residue_index, np.arange(self.n_residues + 1))
        if not np.array_equal(expected, self.residue_atom_offsets):
            raise GeometryError(
                "residue_atom_offsets disagrees with residue_index; the two must "
                "describe the same atom-to-residue mapping."
            )

        if self.chains:
            for chain in self.chains:
                if chain.structure is not self:
                    raise GeometryError(
                        f"Chain {chain.chain_id!r} is bound to a different Structure."
                    )
                for domain in chain.domains:
                    if domain.structure is not self:
                        raise GeometryError(f"{domain!r} is bound to a different Structure.")

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------

    @classmethod
    def from_atom_records(
        cls,
        *,
        xyz: np.ndarray,
        atom_name: Sequence[str],
        element: Sequence[str],
        residue_name: Sequence[str],
        residue_number: Sequence[int],
        chain_id: Sequence[str],
        insertion_code: Sequence[str] | None = None,
        b_factor: Sequence[float] | None = None,
        occupancy: Sequence[float] | None = None,
        source: str | None = None,
    ) -> Structure:
        """Build a :class:`Structure` from flat per-atom records.

        This is the single entry point readers use. It groups atoms into residues and
        chains, assigns the positional indices, and computes the CSR offsets, so that
        no reader has to get those right independently.

        A new residue starts whenever ``(chain_id, residue_number, insertion_code)``
        changes from the previous atom. Grouping on the triple rather than on
        ``residue_number`` alone is what keeps residues 10 and 10A distinct; grouping
        on change rather than on value is what keeps two separate runs of the same
        chain id from being merged, which is how the pre-rewrite reader silently lost
        half the atoms in an interleaved file.

        All input sequences must have the same length as ``xyz``.
        """
        xyz = np.ascontiguousarray(xyz, dtype=np.float64)
        n = xyz.shape[0]
        if n == 0:
            raise EmptyStructureError(f"No atoms to build a structure from (source: {source!r}).")

        atom_name_arr = np.asarray(atom_name, dtype="<U4")
        element_arr = np.asarray([str(e).upper() for e in element], dtype="<U2")
        residue_name_arr = np.asarray(residue_name, dtype="<U3")
        residue_number_arr = np.asarray(residue_number, dtype=np.int64)
        chain_id_arr = np.asarray(chain_id, dtype="<U4")
        insertion_arr = np.asarray(
            [""] * n if insertion_code is None else insertion_code, dtype="<U1"
        )
        b_factor_arr = np.asarray(np.zeros(n) if b_factor is None else b_factor, dtype=np.float64)
        occupancy_arr = np.asarray(np.ones(n) if occupancy is None else occupancy, dtype=np.float64)

        for label, arr in (
            ("atom_name", atom_name_arr),
            ("element", element_arr),
            ("residue_name", residue_name_arr),
            ("residue_number", residue_number_arr),
            ("chain_id", chain_id_arr),
            ("insertion_code", insertion_arr),
            ("b_factor", b_factor_arr),
            ("occupancy", occupancy_arr),
        ):
            if arr.shape[0] != n:
                raise GeometryError(f"{label} has {arr.shape[0]} entries but xyz has {n} atoms.")

        # A new residue begins wherever the identifying triple changes.
        changed = np.ones(n, dtype=bool)
        changed[1:] = (
            (chain_id_arr[1:] != chain_id_arr[:-1])
            | (residue_number_arr[1:] != residue_number_arr[:-1])
            | (insertion_arr[1:] != insertion_arr[:-1])
        )
        residue_index_arr = np.cumsum(changed) - 1
        n_residues = int(residue_index_arr[-1]) + 1

        first_atom_of_residue = np.flatnonzero(changed)
        offsets = np.empty(n_residues + 1, dtype=np.int64)
        offsets[:-1] = first_atom_of_residue
        offsets[-1] = n

        # Residue-level annotations are taken from the alpha carbon where there is one,
        # else the residue's first atom. B-factor in particular is per-atom in the file
        # but per-residue in meaning for AlphaFold pLDDT.
        representative = first_atom_of_residue.copy()
        is_ca = atom_name_arr == "CA"
        ca_atom_indices = np.flatnonzero(is_ca)
        representative[residue_index_arr[ca_atom_indices]] = ca_atom_indices

        chain_ids_per_residue = chain_id_arr[first_atom_of_residue]
        # Chain boundaries: a new chain wherever the chain id changes between
        # consecutive residues. Same change-based logic, same reason.
        chain_changed = np.ones(n_residues, dtype=bool)
        chain_changed[1:] = chain_ids_per_residue[1:] != chain_ids_per_residue[:-1]
        chain_index_arr = np.cumsum(chain_changed) - 1

        structure = cls(
            xyz=xyz,
            atom_name=atom_name_arr,
            element=element_arr,
            residue_index=residue_index_arr.astype(np.int64),
            residue_name=residue_name_arr[first_atom_of_residue],
            residue_number=residue_number_arr[first_atom_of_residue],
            insertion_code=insertion_arr[first_atom_of_residue],
            b_factor=b_factor_arr[representative],
            occupancy=occupancy_arr[representative],
            chain_index=chain_index_arr.astype(np.int64),
            residue_atom_offsets=offsets,
            source=source,
        )

        chain_starts = np.flatnonzero(chain_changed)
        chain_stops = np.append(chain_starts[1:], n_residues)
        for start, stop in zip(chain_starts, chain_stops, strict=True):
            structure.chains.append(
                Chain(
                    structure=structure,
                    span=Span(int(start), int(stop)),
                    chain_id=str(chain_ids_per_residue[start]),
                )
            )

        # Grouping on *change* rather than on value means a file whose records are
        # interleaved (chain A, chain B, chain A again; or residue 1, 2, 1) yields
        # repeated residue keys rather than losing atoms. That is the right tradeoff --
        # the pre-rewrite reader silently discarded roughly half the atoms of an
        # interleaved file -- but a repeated key still signals a malformed input, so
        # say so instead of leaving the caller to wonder why the residue count is high.
        identifiers = [
            (str(chain), int(number), str(icode))
            for chain, number, icode in zip(
                chain_ids_per_residue,
                structure.residue_number,
                structure.insertion_code,
                strict=True,
            )
        ]
        seen: set[tuple[str, int, str]] = set()
        repeated_ids: list[tuple[str, int, str]] = []
        for identifier in identifiers:
            if identifier in seen and identifier not in repeated_ids:
                repeated_ids.append(identifier)
            seen.add(identifier)
        if repeated_ids:
            examples = ", ".join(
                f"{chain}:{number}{icode}" for chain, number, icode in repeated_ids[:5]
            )
            more = f" and {len(repeated_ids) - 5} more" if len(repeated_ids) > 5 else ""
            structure.notes.append(
                f"{len(repeated_ids)} residue identifier(s) appear in more than one run "
                f"of atom records ({examples}{more}). The atoms were kept and treated as "
                f"separate residues rather than merged or dropped, but this usually "
                f"means the input file's records are out of order."
            )

        structure.validate()
        return structure
