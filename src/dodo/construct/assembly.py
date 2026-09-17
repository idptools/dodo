"""Rigid units: which folded domains must move together, and why.

DODO's step 3 moves folded domains so that each linker can adopt its predicted end-to-end
distance. On a single chain that is exactly right. On a **complex** it is a disaster: the domains
of two chains that pack against each other have a relative position that came from the
prediction or the experiment, and pulling them apart to satisfy a linker prediction dismantles
the thing the user handed us. Measured on two copies of dnmt3a placed in contact, step 3 took an
interface of 5,540 atom-atom contacts to zero and the rebuild reported ``ok``.

The fix is one concept, applied everywhere:

    A **rigid unit** is a set of folded domains -- possibly from different chains -- whose
    relative positions and orientations DODO will not change. A unit moves as one rigid body,
    or it does not move at all.

Today's single-chain behaviour is the special case "every folded domain is its own unit". A
complex is the case "all of them are in one unit". Both fall out of the same code, which is the
point: there is no separate complex code path to keep in sync.

Why domains and not chains
--------------------------
A chain in a complex routinely has some domains inside the assembly and some dangling off a long
linker. Locking whole chains would forfeit step 3 -- the step that makes DODO work -- for every
chain that touches another. Locking *domains* keeps it for exactly the domains that are free to
move.

Three ways a pair of domains ends up in one unit
------------------------------------------------
1. **Covalent continuity.** Two folded domains of one chain with no residues between them are
   bonded across the boundary. Separating them breaks a peptide bond, and it used to: through
   ``assign_regions_from_spec`` this took a C-N bond from 1.334 A to 93.479 A while reporting
   "their relative position is not ours to change".
2. **A link that will not be rebuilt.** Residues between two folded domains that DODO is going
   to skip -- anything under ``min_length`` -- keep their input coordinates, so they are a rigid
   link and their neighbours are one unit. Without this, ``rebuild(min_length=50)`` on dnmt3a
   stretched a peptide bond to 33.605 A.
3. **An interface.** Two domains in contact. This is the complex case, and it is the only one of
   the three that is a judgement call rather than a fact; see
   :data:`~dodo.constants.INTERFACE_MIN_RESIDUE_PAIRS` for why the threshold errs toward locking.

Rules 1 and 2 are properties of the chain, so they hold for a single-chain input too. Rule 3 is
restricted to pairs of domains in *different* chains by default: two folded domains of one chain
that pack together are exactly the arrangement DODO exists to re-sample, and locking them by
default would change behaviour validated over 23,587 single-chain structures.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field, replace
from itertools import pairwise
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from ..constants import (
    INTERFACE_CONTACT_RADIUS,
    INTERFACE_MIN_RESIDUE_PAIRS,
    MIN_IDR_LENGTH,
)
from ..exceptions import GeometryError, InvalidRegionError
from ..geometry.transforms import superpose
from ..structure import Domain, DomainKind, Structure
from .unmodelled import skipped_for_length

__all__ = [
    "Assembly",
    "Interface",
    "RigidUnit",
    "find_rigid_units",
    "units_from_spec",
]


# ----------------------------------------------------------------------------------
# Why a pair was locked. Strings rather than an enum because they are report text.
# ----------------------------------------------------------------------------------

#: Locked because the two domains are sequence-adjacent: there is a peptide bond across the join.
REASON_ADJACENT = "covalently adjacent"
#: Locked because everything between the two domains keeps its input coordinates.
REASON_UNBUILT_LINK = "joined by residues that will not be rebuilt"
#: Locked because the input never modelled the linker between them.
REASON_UNMODELLED_LINK = "joined by a linker the input did not model"
#: Locked because the two domains are in physical contact.
REASON_INTERFACE = "interface"
#: Locked because the caller said so.
REASON_SPEC = "specified by the caller"
#: Locked because the caller asked for whole chains to be held rigid.
REASON_WHOLE_CHAIN = "whole chain held rigid"


@dataclass(frozen=True, slots=True)
class Interface:
    """A contact between two folded domains, and whether it locks them together.

    Reported whether or not it locked, because a contact DODO decided to break is exactly the
    thing a user needs to see. ``residues_a`` and ``residues_b`` are 1-based inclusive residue
    numbers as the input file numbers them -- the same numbering ``dodo regions`` prints.
    """

    chain_a: str
    residues_a: tuple[int, int]
    chain_b: str
    residues_b: tuple[int, int]
    #: Distinct residue-residue contacts within the contact radius. The quantity thresholded.
    residue_pairs: int
    #: Atom-atom contacts within the contact radius. Reported for scale, never thresholded --
    #: it varies with how many side chains the input models.
    atom_pairs: int
    #: Whether the two domains are in different chains.
    inter_chain: bool
    locked: bool
    #: Whether the two domains belong to the same final rigid unit. This differs from
    #: ``locked`` when another chain of joins connects them transitively.
    preserved: bool = False

    def __str__(self) -> str:
        where = (
            f"{self.chain_a}:{self.residues_a[0]}-{self.residues_a[1]} .. "
            f"{self.chain_b}:{self.residues_b[0]}-{self.residues_b[1]}"
        )
        verdict = (
            "locked"
            if self.locked
            else "preserved through another rigid-unit join"
            if self.preserved
            else "NOT locked"
        )
        return (
            f"{where}: {self.residue_pairs} residue contact(s), "
            f"{self.atom_pairs} atom contact(s) -- {verdict}"
        )


@dataclass(slots=True)
class RigidUnit:
    """A set of folded domains that move together, or not at all.

    Owns no coordinates. Every transform is applied to each member domain about a **single
    common point**, and that is what makes the motion rigid across the whole unit rather than
    merely within each domain. Rotating each domain about its own centroid would preserve every
    domain's internal geometry and still take the unit apart, and that is precisely the bug
    this class exists to make impossible to write.
    """

    index: int
    domains: list[Domain]
    #: Why each merge happened, in the order the merges were made. Report text.
    reasons: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.domains:
            raise InvalidRegionError("A rigid unit must contain at least one folded domain.")
        self.domains.sort(key=lambda d: d.span.start)

    def __len__(self) -> int:
        return len(self.domains)

    def __iter__(self) -> Iterator[Domain]:
        return iter(self.domains)

    def __repr__(self) -> str:
        return f"RigidUnit({self.index}, {len(self.domains)} domain(s), {self.n_atoms} atoms)"

    @property
    def structure(self) -> Structure:
        """The structure these domains view."""
        return self.domains[0].structure

    @property
    def n_atoms(self) -> int:
        """Total atoms across every domain of this unit."""
        return sum(len(d.xyz) for d in self.domains)

    @property
    def atom_mask(self) -> np.ndarray:
        """Boolean mask over the parent structure's atoms selecting this unit."""
        mask = np.zeros(self.structure.n_atoms, dtype=bool)
        for domain in self.domains:
            mask[domain.atom_slice] = True
        return mask

    @property
    def xyz(self) -> np.ndarray:
        """Return the unit's atoms as one ``(n_atoms, 3)`` array.

        A copy, not a view: a unit's domains occupy disjoint slices of the parent array.
        """
        return np.concatenate([d.xyz for d in self.domains], axis=0)

    def chain_ids(self) -> tuple[str, ...]:
        """Chain ids this unit spans, in structure order, without repeats."""
        structure = self.structure
        seen: list[str] = []
        for domain in self.domains:
            chain_id = structure.chains[int(structure.chain_index[domain.span.start])].chain_id
            if chain_id not in seen:
                seen.append(chain_id)
        return tuple(seen)

    def centroid(self) -> np.ndarray:
        """Geometric centroid over every atom of the unit."""
        centroid: np.ndarray = self.xyz.mean(axis=0)
        return centroid

    def bounding_sphere(self) -> tuple[np.ndarray, float]:
        """Centre and radius of a sphere containing every atom. For cheap clash prefiltering."""
        coords = self.xyz
        centre = coords.mean(axis=0)
        return centre, float(np.linalg.norm(coords - centre, axis=1).max())

    def translate(self, vector: Sequence[float] | np.ndarray) -> None:
        """Translate every domain of the unit."""
        for domain in self.domains:
            domain.translate(vector)

    def rotate(self, rotation: np.ndarray, about: np.ndarray) -> None:
        """Rotate every domain of the unit about one common point.

        ``about`` is required, deliberately. :meth:`Domain.rotate` defaults it to the domain's own
        centroid, which for a multi-domain unit is the one thing it must never be.
        """
        about = np.asarray(about, dtype=np.float64)
        for domain in self.domains:
            domain.rotate(rotation, about=about)

    def describe(self, structure: Structure | None = None) -> str:
        """One-line human-readable description, in input residue numbering."""
        structure = structure if structure is not None else self.structure
        parts = []
        for domain in self.domains:
            chain_id = structure.chains[int(structure.chain_index[domain.span.start])].chain_id
            first = int(structure.residue_number[domain.span.start])
            last = int(structure.residue_number[domain.span.stop - 1])
            parts.append(f"{chain_id}:{first}-{last}")
        reasons = f" [{'; '.join(self.reasons)}]" if self.reasons else ""
        return f"unit {self.index}: " + ", ".join(parts) + reasons


@dataclass(slots=True)
class Assembly:
    """Every rigid unit in a structure, and the interfaces that produced them."""

    units: list[RigidUnit] = field(default_factory=list)
    interfaces: list[Interface] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    #: Maps a domain's ``span.start`` to its unit index. Spans tile chains without overlap, so
    #: the start index identifies a domain uniquely -- and unlike ``id(domain)`` it survives
    #: :meth:`Structure.copy`, which rebinds every Domain view to a new structure.
    _unit_of_start: dict[int, int] = field(default_factory=dict)

    def unit_of(self, domain: Domain) -> RigidUnit:
        """Return the unit containing ``domain``."""
        try:
            return self.units[self._unit_of_start[domain.span.start]]
        except KeyError:
            raise InvalidRegionError(
                f"{domain!r} is not part of this assembly. Rigid units cover folded domains "
                f"only, and they are computed for one structure at a time."
            ) from None

    def same_unit(self, first: Domain, second: Domain) -> bool:
        """Report whether both domains move together."""
        return (
            self._unit_of_start[first.span.start] == self._unit_of_start[second.span.start]
        )

    @property
    def locked_interfaces(self) -> list[Interface]:
        """Interfaces that put two domains in one unit."""
        return [i for i in self.interfaces if i.locked]

    @property
    def broken_interfaces(self) -> list[Interface]:
        """Contacts DODO is prepared to break: below threshold, or intra-chain when not locking."""
        return [i for i in self.interfaces if not i.preserved]

    @property
    def n_locked_units(self) -> int:
        """Units holding more than one folded domain."""
        return sum(1 for u in self.units if len(u) > 1)

    def summary(self) -> str:
        """Multi-line human-readable summary."""
        lines = [
            f"{len(self.units)} rigid unit(s) over "
            f"{sum(len(u) for u in self.units)} folded domain(s); "
            f"{self.n_locked_units} hold more than one"
        ]
        lines += [f"  {unit.describe()}" for unit in self.units]
        if self.interfaces:
            lines.append(
                f"  {len(self.locked_interfaces)}/{len(self.interfaces)} contact(s) locked"
            )
            lines += [f"    {i}" for i in self.interfaces]
        lines += [f"  note: {n}" for n in self.notes]
        return "\n".join(lines)


# ----------------------------------------------------------------------------------
# Union-find. Small enough to inline, and a dependency-free structure is worth more
# here than a clever one: the merge ORDER is what makes unit membership deterministic.
# ----------------------------------------------------------------------------------


class _DisjointSet:
    def __init__(self, n: int) -> None:
        self._parent = list(range(n))

    def find(self, item: int) -> int:
        root = item
        while self._parent[root] != root:
            root = self._parent[root]
        while self._parent[item] != root:  # path compression
            self._parent[item], item = root, self._parent[item]
        return root

    def union(self, a: int, b: int) -> bool:
        """Merge, keeping the LOWER root so unit order follows residue order. True if merged."""
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        low, high = (ra, rb) if ra < rb else (rb, ra)
        self._parent[high] = low
        return True


def _folded_domains_by_chain(structure: Structure) -> list[list[Domain]]:
    """Folded domains of each chain, in residue order. Chains with none give an empty list."""
    return [
        sorted(
            (d for d in chain.domains if d.kind is DomainKind.FOLDED),
            key=lambda d: d.span.start,
        )
        for chain in structure.chains
    ]


def _residue_label(structure: Structure, domain: Domain) -> tuple[int, int]:
    return (
        int(structure.residue_number[domain.span.start]),
        int(structure.residue_number[domain.span.stop - 1]),
    )


def _contacts_between_domains(
    structure: Structure,
    domains: Sequence[Domain],
    *,
    radius: float,
    on_progress: Callable[[int], None] | None = None,
) -> dict[tuple[int, int], tuple[int, int]]:
    """Count contacts for every candidate domain pair.

    Returns ``{(i, j): (residue_pairs, atom_pairs)}`` for pairs with at least one contact, where
    ``i < j`` index into ``domains``.

    Two-stage by design. A bounding-sphere test rejects most pairs for the cost of one vector
    norm, which is what keeps this tractable on an assembly: on a large complex the great
    majority of domain pairs are nowhere near each other, and building a spatial index for those
    is pure waste. Only survivors get a KD-tree query.

    Every pair is measured, including the intra-chain pairs that will not be locked: a contact
    DODO is prepared to break is exactly what a user needs to see in the report, and the
    bounding-sphere prefilter makes measuring it nearly free.
    """
    n = len(domains)
    contacts: dict[tuple[int, int], tuple[int, int]] = {}
    if n < 2:
        return contacts

    centres = np.empty((n, 3), dtype=np.float64)
    radii = np.empty(n, dtype=np.float64)
    for k, domain in enumerate(domains):
        coords = domain.xyz
        centres[k] = coords.mean(axis=0)
        radii[k] = np.linalg.norm(coords - centres[k], axis=1).max()

    separation = np.linalg.norm(centres[:, None, :] - centres[None, :, :], axis=-1)
    reachable = separation <= (radii[:, None] + radii[None, :] + radius)
    candidates = np.triu(reachable, k=1)

    trees: dict[int, cKDTree] = {}

    def tree_for(k: int) -> cKDTree:
        if k not in trees:
            trees[k] = cKDTree(domains[k].xyz)
        return trees[k]

    for i, j in zip(*np.nonzero(candidates), strict=True):
        i, j = int(i), int(j)
        if on_progress is not None:
            on_progress(1)
        hits = tree_for(i).sparse_distance_matrix(tree_for(j), radius, output_type="ndarray")
        if hits.size == 0:
            continue
        # Local atom offsets -> structure residue indices, so a contact is counted once per
        # residue pair however many atoms of those residues happen to be in range. Atom counts
        # scale with how many side chains the input models; residue counts do not.
        atoms_i = structure.residue_index[domains[i].atom_slice][hits["i"]]
        atoms_j = structure.residue_index[domains[j].atom_slice][hits["j"]]
        residue_pairs = np.unique(np.stack([atoms_i, atoms_j], axis=1), axis=0).shape[0]
        contacts[(i, j)] = (int(residue_pairs), int(hits.size))
    return contacts


def find_rigid_units(
    structure: Structure,
    *,
    min_length: int = MIN_IDR_LENGTH,
    contact_radius: float = INTERFACE_CONTACT_RADIUS,
    min_residue_pairs: int = INTERFACE_MIN_RESIDUE_PAIRS,
    lock_intra_chain_interfaces: bool = False,
    lock_interfaces: bool = True,
    lock_chains: bool = False,
    lock_unmodelled_links: bool = True,
    on_progress: Callable[[int], None] | None = None,
) -> Assembly:
    """Group a structure's folded domains into rigid units.

    Evaluated on the **input** coordinates, before anything moves, and it must stay that way:
    once step 3 has run, contacts are a description of what DODO did rather than of what the
    input said.

    Parameters
    ----------
    structure
        A structure whose regions have already been assigned.
    min_length
        The rebuild's minimum region length. Regions shorter than this keep their input
        coordinates, which makes them a rigid link between the domains they join. This is why
        units depend on a *rebuild* parameter and are therefore computed inside the pipeline
        rather than during region identification.
    contact_radius
        Heavy-atom distance at which two domains count as touching.
    min_residue_pairs
        Residue-residue contacts at which a contact locks. See
        :data:`~dodo.constants.INTERFACE_MIN_RESIDUE_PAIRS`.
    lock_intra_chain_interfaces
        Also lock two folded domains of the *same* chain that are in contact. Off by default:
        that arrangement is what step 3 exists to re-sample.
    lock_interfaces
        Master switch. ``False`` keeps rules 1 and 2 -- which are facts about covalent geometry,
        not judgement calls -- and drops rule 3, giving one unit per folded domain except where
        the chain is continuous across the join.
    on_progress
        Called with the number of domain pairs measured since the last call. On an assembly the
        contact scan is the only part of this that takes real time, and it is quadratic in the
        number of folded domains, so it is the part worth reporting.
    lock_unmodelled_links
        Hold together two folded domains whose connecting linker was **inserted** rather than
        observed -- see :mod:`dodo.construct.unmodelled`. On by default, and only ever active on
        a structure whose unmodelled residues were filled in, since nothing else sets
        :attr:`Structure.inserted`.
    lock_chains
        Treat every chain as one rigid body: all of a chain's folded domains go in one unit
        whether or not they touch. Interfaces still lock on top, so a complex of rigid chains
        comes out as a complex. For an experimental structure whose linkers were never modelled
        this is often what the user means -- the deposited arrangement is the measurement, and
        DODO's job is only to fill in the disordered regions.

    Returns
    -------
    Assembly
        The units, every contact found (locked or not), and notes.
    """
    domains: list[Domain] = []
    chain_of: list[int] = []
    for chain_index, chain_domains in enumerate(_folded_domains_by_chain(structure)):
        for domain in chain_domains:
            domains.append(domain)
            chain_of.append(chain_index)

    assembly = Assembly()
    if not domains:
        assembly.notes.append("no folded domains, so there is nothing to hold rigid")
        return assembly

    order = np.argsort([d.span.start for d in domains], kind="stable")
    domains = [domains[k] for k in order]
    chain_of = [chain_of[k] for k in order]

    joined = _DisjointSet(len(domains))
    index_of_start = {d.span.start: k for k, d in enumerate(domains)}
    reasons: dict[int, list[str]] = {}

    def merge(a: int, b: int, reason: str) -> None:
        root_a, root_b = joined.find(a), joined.find(b)
        if joined.union(a, b):
            root = joined.find(a)
            merged = [*reasons.pop(root_a, ()), *reasons.pop(root_b, ()), reason]
            reasons[root] = list(dict.fromkeys(merged))

    # --- rule 0, opt-in: the whole chain is one rigid body ----------------------------
    if lock_chains:
        for chain_domains in _folded_domains_by_chain(structure):
            for previous, current in pairwise(chain_domains):
                merge(
                    index_of_start[previous.span.start],
                    index_of_start[current.span.start],
                    REASON_WHOLE_CHAIN,
                )

    # --- rules 1 and 2: the chain is continuous across the join -----------------------
    for chain_domains in _folded_domains_by_chain(structure):
        for previous, current in pairwise(chain_domains):
            a = index_of_start[previous.span.start]
            b = index_of_start[current.span.start]
            between = current.span.start - previous.span.stop
            if between == 0:
                merge(a, b, REASON_ADJACENT)
                continue
            # Everything between them keeps its input coordinates iff no region in the gap will
            # be rebuilt. In a valid tiling the gap is exactly one IDR, but the check is written
            # over whatever is actually there so a caller-supplied assignment cannot slip past.
            gap = [
                d
                for d in structure.chains[
                    int(structure.chain_index[previous.span.start])
                ].domains
                if previous.span.stop <= d.span.start and d.span.stop <= current.span.start
            ]
            if gap and all(
                d.kind is DomainKind.IDR and skipped_for_length(structure, d.span, min_length)
                for d in gap
            ):
                merge(a, b, REASON_UNBUILT_LINK)
                continue
            # A linker the input never modelled says nothing about where these two domains sit
            # relative to each other -- but their own coordinates do, because an experiment or a
            # prediction put them there and DODO has no better information. Moving them apart to
            # satisfy a prediction for residues nobody ever saw would replace a measurement with
            # a guess. So an unobserved linker locks its neighbours; an observed one does not.
            if lock_unmodelled_links and bool(
                structure.inserted[previous.span.stop : current.span.start].any()
            ):
                merge(a, b, REASON_UNMODELLED_LINK)

    # --- rule 3: interfaces ------------------------------------------------------------
    if lock_interfaces:
        contacts = _contacts_between_domains(
            structure, domains, radius=contact_radius, on_progress=on_progress
        )
        interface_members: list[tuple[int, int]] = []
        for (i, j), (residue_pairs, atom_pairs) in sorted(contacts.items()):
            inter_chain = chain_of[i] != chain_of[j]
            locked = residue_pairs >= min_residue_pairs and (
                inter_chain or lock_intra_chain_interfaces
            )
            assembly.interfaces.append(
                Interface(
                    chain_a=structure.chains[chain_of[i]].chain_id,
                    residues_a=_residue_label(structure, domains[i]),
                    chain_b=structure.chains[chain_of[j]].chain_id,
                    residues_b=_residue_label(structure, domains[j]),
                    residue_pairs=residue_pairs,
                    atom_pairs=atom_pairs,
                    inter_chain=inter_chain,
                    preserved=False,
                    locked=locked,
                )
            )
            interface_members.append((i, j))
            if locked:
                merge(i, j, f"{REASON_INTERFACE} ({residue_pairs} residue contacts)")
        assembly.interfaces = [
            replace(interface, preserved=joined.find(i) == joined.find(j))
            for interface, (i, j) in zip(assembly.interfaces, interface_members, strict=True)
        ]
        if not lock_intra_chain_interfaces:
            assembly.notes.append(
                "intra-chain contacts between folded domains were not locked; that arrangement "
                "is what folded-domain repositioning exists to re-sample. Pass "
                "lock_intra_chain_interfaces=True to hold them too."
            )
    else:
        assembly.notes.append(
            "interface locking is off, so folded domains are held together only where the "
            "chain is continuous across the join"
        )

    _build_units(assembly, domains, joined, reasons)
    return assembly


def _build_units(
    assembly: Assembly,
    domains: Sequence[Domain],
    joined: _DisjointSet,
    reasons: Mapping[int, list[str]],
) -> None:
    """Turn the disjoint set into ordered units. Ordering is by earliest domain, always."""
    groups: dict[int, list[int]] = {}
    for k in range(len(domains)):
        groups.setdefault(joined.find(k), []).append(k)

    for unit_index, root in enumerate(sorted(groups, key=lambda r: min(groups[r]))):
        members = sorted(groups[root])
        unit = RigidUnit(
            index=unit_index,
            domains=[domains[k] for k in members],
            reasons=tuple(dict.fromkeys(reasons.get(root, ()))),
        )
        assembly.units.append(unit)
        for k in members:
            assembly._unit_of_start[domains[k].span.start] = unit_index


def units_from_spec(
    structure: Structure,
    spec: Sequence[Sequence[Any]],
    *,
    min_length: int = MIN_IDR_LENGTH,
) -> Assembly:
    """Build units from an explicit caller specification.

    Parameters
    ----------
    structure
        A structure whose regions have already been assigned.
    spec
        A sequence of units, each a sequence of ``(chain_id, residue_number)`` pairs naming any
        residue inside a folded domain that belongs to the unit::

            [[("A", 300), ("B", 1200)], [("A", 700)]]

        Residue numbers are as the input file numbers them, matching
        :func:`~dodo.regions.identify.assign_regions_from_spec` and ``dodo regions``.
        Folded domains not named anywhere become units of their own, so a caller only has to
        write down the domains they want held together.

    Notes
    -----
    Rules 1 and 2 of :func:`find_rigid_units` are still applied on top of the specification. They
    are statements about covalent geometry -- domains bonded to each other, or joined by residues
    that will not be rebuilt -- and honouring a specification that contradicts them would mean
    knowingly writing a broken chain.
    """
    base = find_rigid_units(structure, min_length=min_length, lock_interfaces=False)
    domains = [d for unit in base.units for d in unit.domains]
    domains.sort(key=lambda d: d.span.start)
    index_of_start = {d.span.start: k for k, d in enumerate(domains)}

    joined = _DisjointSet(len(domains))
    reasons: dict[int, list[str]] = {}
    # Re-apply the covalent rules, which base.units already encodes.
    for unit in base.units:
        for previous, current in pairwise(unit.domains):
            a, b = index_of_start[previous.span.start], index_of_start[current.span.start]
            if joined.union(a, b):
                reasons.setdefault(joined.find(a), []).extend(unit.reasons)

    by_chain = {chain.chain_id: chain for chain in structure.chains}
    for entry in spec:
        members: list[int] = []
        for item in entry:
            chain_id, residue_number = item
            chain = by_chain.get(str(chain_id))
            if chain is None:
                raise InvalidRegionError(
                    f"Chain {chain_id!r} is not in this structure. "
                    f"Available: {sorted(by_chain)}."
                )
            chain_index = next(
                k for k, other in enumerate(structure.chains) if other.chain_id == chain.chain_id
            )
            matches = np.flatnonzero(
                (structure.residue_number == int(residue_number))
                & (structure.chain_index == chain_index)
            )
            if matches.size == 0:
                raise InvalidRegionError(
                    f"Chain {chain_id!r} has no residue numbered {residue_number}."
                )
            domain = chain.domain_at(int(matches[0]))
            if domain is None or domain.kind is not DomainKind.FOLDED:
                raise InvalidRegionError(
                    f"Residue {chain_id}:{residue_number} is not inside a folded domain, so it "
                    f"cannot name a rigid unit. Rigid units group folded domains only."
                )
            members.append(index_of_start[domain.span.start])
        for a, b in pairwise(members):
            if joined.union(a, b):
                reasons.setdefault(joined.find(a), []).append(REASON_SPEC)

    assembly = Assembly(notes=["rigid units supplied by the caller; no contacts were measured"])
    _build_units(assembly, domains, joined, reasons)
    return assembly


def verify_units_rigid(
    before: Mapping[int, np.ndarray], assembly: Assembly, *, tolerance: float = 1e-6
) -> None:
    """Assert that every unit moved as one rigid body.

    ``before`` maps a unit index to that unit's atom coordinates before the move, in the order
    :attr:`RigidUnit.xyz` returns them.

    This is the invariant of the whole mechanism, and it is deliberately checked across the
    *unit*, not the domain. Rotating each domain about its own centroid preserves every domain's
    internal geometry perfectly and still takes the assembly apart, so a per-domain check would
    pass on exactly the bug that matters. Fitting the best proper rigid transform and checking
    its per-atom residual is O(n) and detects rearrangements that preserve centroid radii.
    """
    for unit in assembly.units:
        original = before.get(unit.index)
        if original is None:
            continue
        current = unit.xyz
        if original.shape != current.shape:
            raise GeometryError(
                f"Cannot verify unit {unit.index}: it had {original.shape[0]} atoms before and "
                f"has {current.shape[0]} now."
            )
        if original.shape[0] < 2:
            continue
        drift = _rigid_drift(original, current)
        if drift > tolerance:
            raise GeometryError(
                f"Rigid unit {unit.index} was deformed rather than moved: the worst atom's "
                f"distance to the unit centroid changed by {drift:.3e} A. Its {len(unit)} "
                f"folded domain(s) must keep their relative positions exactly."
            )


def _rigid_drift(before: np.ndarray, after: np.ndarray) -> float:
    """Maximum residual after the best proper rigid superposition of two point sets."""
    rotation, translation = superpose(before, after)
    aligned = before @ rotation.T + translation
    return float(np.linalg.norm(aligned - after, axis=1).max())
