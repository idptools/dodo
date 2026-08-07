"""Bond-length validation for every bond in a structure, not just the backbone.

What this checks
----------------
Three classes of bond, each against the measured reference in
:mod:`dodo.validate.reference`:

1. **Intra-residue** bonds, against :data:`~dodo.validate.reference.RESIDUE_BONDS` -- every
   measured bond of all twenty standard amino acids, side chains included. This is the part v1
   never had: its hand-written side-chain template library was wrong by up to 2.8 A and nothing
   ever measured it.
2. **Inter-residue peptide** bonds, C(i)-N(i+1), against
   :data:`~dodo.validate.reference.PEPTIDE_BOND_LENGTH`.
3. **CA-CA virtual** bonds between consecutive residues, against
   :data:`dodo.constants.CA_CA_BOND_LENGTH`. This is the only bond a rebuilt region has, so it
   is the load-bearing check for everything DODO generates.

What DODO's output actually looks like, and why it drives the design
--------------------------------------------------------------------
A DODO model is deliberately NOT a uniform all-atom structure. Folded domains keep every atom of
the input and are moved only as rigid bodies; rebuilt regions (IDRs, and loops inside folded
domains) are **alpha carbon only**, at 3.81 A spacing, with nothing else. The two kinds of
stretch sit in the same chain and the boundary between them is a real chain connection.

So a validator that demanded N, CA, C and O everywhere would emit one useless finding per rebuilt
residue -- 334 on a rebuilt dnmt3a, 1,526 on a rebuilt p300. Every check here therefore asks
first whether a residue is CA-only and, if so, expects nothing from it beyond its virtual CA-CA
bonds.

Real input contains real oddities, and preserving one is not a bug
------------------------------------------------------------------
The distinction that matters is between a defect DODO *introduced* and one it faithfully
*carried through*. Each exclusion below is on a stated, measured basis:

* **cis and strained peptide bonds.** A cis peptide puts CA(i)-CA(i+1) near 3.19 A rather than
  3.83 A (:data:`~dodo.validate.reference.CA_CA_SHORT_LENGTH`, 2.44% of the reference
  population), and a non-planar omega can push it to 4.4 A -- while the actual C-N bond is
  perfect. Measured on this repository's fixtures, 26% of p300's consecutive pairs deviate from
  3.81 A by more than 0.15 A while *every one* of its 2,413 peptide bonds measures between 1.317
  and 1.407 A. The CA-CA virtual bond is a surrogate for connectivity, so where a real peptide
  bond can be measured it is checked instead and the surrogate is skipped. Pass
  ``strict_ca_ca=True`` to check both.
* **preserved CA-only geometry.** Where a CA-CA bond cannot be cross-checked against a peptide
  bond, what counts as acceptable depends on who built it. Geometry DODO generated is held to
  the 3.81 A it builds to; geometry DODO merely preserved is flagged only when it falls outside
  *both* measured reference populations. Measured on the AlphaFold fixtures' CA traces, that is
  the difference between flagging 6.8-22.8% and flagging 0-0.7%.
* **chain breaks at unmodelled residues.** Where the deposited numbering jumps, the two residues
  are not chemically adjacent and neither their peptide bond nor their CA-CA distance means
  anything. Those pairs are counted and noted, never flagged.
* **modified residues** (MSE, SEP, ...). DODO keeps them on purpose. They are absent from the
  reference table, so they are noted with their names and counts rather than flagged.
* **partially built residues.** A residue with a backbone and no side chain is what
  ``rebuild(backbone=True)`` produces, and also how deposited structures model side-chain
  disorder. See :data:`_TIER_BACKBONE`.

Reporting, not raising
----------------------
:func:`validate_bonds` returns a :class:`BondReport` that collects every finding, following the
:class:`~dodo.geometry.metrics.TraceReport` pattern. A caller wants the whole list, not the first
item. Exceptions are reserved for input that cannot be validated at all (a nonsensical tolerance,
something that is not a :class:`~dodo.structure.Structure`); :meth:`BondReport.raise_if_invalid`
is there for callers who want the strict form.

Performance
-----------
Everything is vectorized over numpy arrays or delegated to a KD-tree; nothing loops over atom
pairs in Python. Measured on the 61,511-atom, 7,781-residue 6kn7 EM assembly: **27 ms** for the
whole structure, 62,548 bonds checked -- of which the KD-tree pass over all 62,548 pairs inside
1.90 A is 15 ms. The intra-residue pass runs 20 vectorized passes, one per residue type, rather
than one per residue.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Final, Literal

import numpy as np
from scipy.spatial import cKDTree

from ..constants import CA_CA_BOND_LENGTH
from ..exceptions import GeometryError
from ..structure import Structure
from .reference import (
    BOND_CUTOFF,
    CA_CA_SHORT_LENGTH,
    CA_CA_TRANS_LENGTH,
    PEPTIDE_BOND_LENGTH,
    PEPTIDE_BOND_SD,
    RESIDUE_BONDS,
)

__all__ = [
    "DEFAULT_ABSOLUTE_TOLERANCE",
    "DEFAULT_CA_CA_TOLERANCE",
    "DEFAULT_TOLERANCE_SIGMA",
    "BondClass",
    "BondReport",
    "BondViolation",
    "BondViolationKind",
    "Provenance",
    "validate_bonds",
]

# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------

#: Default slack in standard deviations added outside the measured 0.1-99.9 percentile range.
#:
#: The accepted window for an intra-residue bond is
#: ``[p0.1 - slack, p99.9 + slack]`` with ``slack = max(tolerance_sigma * sd,
#: absolute_tolerance)``. It is anchored on the measured percentiles rather than on ``mean +/-
#: slack`` because several reference distributions are strongly asymmetric -- ALA C-CA is
#: 1.5422 A with an sd of 0.0095 but a 99.9th percentile of 1.5994, which is +6.0 sd on one side
#: and -2.5 sd on the other. A symmetric window either flags the long tail or waves through the
#: short one. Measured, the asymmetric window is strictly better: at a *narrower* minimum
#: half-width it flags 0.52% of 6kn7's bonds where the symmetric form flags 0.91%.
#:
#: MEASURED to this value over the 89,456 intra-residue bonds of the five all-atom fixtures in
#: this repository. Overall flag rate 0.35%, broken down by input type:
#:
#:     arf19 (AlphaFold)   0.000%      p300 (AlphaFold)    0.109%
#:     dnmt3a (AlphaFold)  0.047%      6kn7 (experimental) 0.524%
#:
#: The split is expected and is documented in :mod:`dodo.validate.reference`: the reference is
#: measured from AlphaFold models, and experimental refinement targets differ systematically --
#: aromatic C-C and carboxylate C-O by up to 0.032 A. So an experimental structure legitimately
#: flags an order of magnitude more than a predicted one, and DODO's input is predicted
#: structures. Callers validating crystal or EM input should expect the higher rate or raise
#: ``absolute_tolerance``.
#:
#: Sensitivity is not the constraint. The narrowest window this produces is +/-0.06 A, while the
#: defects this suite exists to catch are an order of magnitude larger: v1's all-atom module put
#: CA-C at exactly 1.000 A instead of 1.535, and dnmt3a's own A:HIS613 has an imidazole ring
#: modelled with ND1-CE1 at 2.14 A instead of 1.33 A. Both are flagged by more than 10x.
DEFAULT_TOLERANCE_SIGMA: Final[float] = 6.0

#: Absolute floor on the slack, in Angstroms.
#:
#: MEASURED. Reference standard deviations are small -- the median is about 0.005 A, because
#: AlphaFold geometry is more uniform than experimental geometry -- so 6 sd is often 0.03 A,
#: which is tighter than the coordinate precision of a 3.5 A cryo-EM model and tight enough that
#: rounding in a PDB file starts to matter. Sweeping the floor over the fixtures:
#:
#:     0.05 A -> 0.67% flagged      0.06 A -> 0.35% flagged      0.07 A -> 0.20% flagged
#:
#: 0.06 is the knee, and going looser buys little. Nothing at this scale affects the module's
#: ability to catch a real defect.
DEFAULT_ABSOLUTE_TOLERANCE: Final[float] = 0.06

#: Default tolerance on the CA-CA virtual bond, in Angstroms.
#:
#: MEASURED. :data:`dodo.constants.CA_CA_BOND_TOLERANCE` is 0.10, which is the right window for
#: *generation* -- DODO builds at 3.81 A and closes onto anchors within 0.10 A, and
#: :func:`dodo.geometry.metrics.validate_ca_trace` is the gate that enforces it. It is too tight
#: for validating geometry DODO did not generate: the deposited 6kn7 trace spans 3.660-3.970 A
#: over 7,750 trans peptides, so a 0.10 A window flags 1.1% of a real experimental structure and
#: 0.15 A still flags 3 bonds of it. At 0.20 A, 6kn7 is clean.
#:
#: It stays far tighter than any defect worth catching: v1's own output
#: (testing_translation.pdb) contains CA-CA distances of 0.000 A and 56.2 A.
DEFAULT_CA_CA_TOLERANCE: Final[float] = 0.20

#: Elements skipped by the unexpected-contact check.
#:
#: The reference table contains no hydrogen bonds, because it was measured from structures that
#: have none. So every X-H bond -- at ~1.0 A, well inside BOND_CUTOFF -- would be reported as an
#: unexpected contact, one per hydrogen. DODO's readers drop hydrogens, but a Structure can be
#: built directly from records that have them.
_HYDROGEN_ELEMENTS: Final[frozenset[str]] = frozenset({"H", "D"})

#: Atoms never reported as missing, because their absence is normal rather than a defect.
#:
#: OXT exists only on a C-terminal residue -- one per chain -- so demanding it everywhere
#: produces a finding for every other residue in the structure. It is still checked for length
#: wherever it is present.
_OPTIONAL_ATOMS: Final[frozenset[str]] = frozenset({"OXT"})

#: Separation below which two CA positions are called the same point, in Angstroms. Cosmetic:
#: it only decides whether a message says so.
_COINCIDENT: Final[float] = 0.5

#: Kinds of finding :func:`validate_bonds` can report.
BondViolationKind = Literal[
    "bond_length",
    "chain_break",
    "seam",
    "missing_atom",
    "unexpected_contact",
    "non_finite",
]

#: Which class of bond a finding concerns. ``"none"`` for findings that are not about a specific
#: bond, i.e. a non-finite coordinate or an absent atom.
BondClass = Literal["intra_residue", "peptide", "ca_ca", "none"]

#: Where a finding's geometry came from, when the structure carries region assignments.
#:
#: ``"rebuilt"`` means the residue lies in a region DODO regenerates (an IDR, or a loop inside a
#: folded domain), so a defect there is DODO's. ``"input"`` means it lies in a folded stretch,
#: which DODO only ever moves as a rigid body -- and a rigid-body move preserves every
#: intra-residue distance exactly, so an intra-residue defect there provably came from the input
#: file. ``"unknown"`` means the structure carries no domain assignments to consult, which is the
#: case for anything freshly read from disk.
#:
#: ``"seam"`` means the finding is the peptide bond at a rebuilt-region/folded-domain junction,
#: whose strain is inherited from the rigidly-repositioned folded atom rather than caused by DODO's
#: rebuild -- approximate by construction, and attributed to neither side alone.
#:
#: Note that :attr:`dodo.structure.Domain.rebuilt` is deliberately *not* the signal: the pipeline
#: sets it on folded domains too, because there it means "already placed" for clash checking.
Provenance = Literal["rebuilt", "input", "seam", "unknown"]


# ---------------------------------------------------------------------------
# Findings
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class BondViolation:
    """One specific bond-geometry defect, with everything needed to act on it.

    Attributes
    ----------
    kind
        Which check failed. ``"bond_length"`` is a measured distance outside the accepted window;
        ``"chain_break"`` is a C-N pair between consecutively numbered residues too far apart to
        be a bond at all; ``"missing_atom"`` is an atom the residue type requires that is absent
        from a residue meant to have it; ``"unexpected_contact"`` is a pair of atoms in one
        residue closer than :data:`dodo.validate.reference.BOND_CUTOFF` that is not a bond of
        that residue type; ``"non_finite"`` is a NaN or infinite coordinate.
    bond_class
        Which of the three bond classes this concerns, or ``"none"``.
    residue_indices
        0-based positional residue indices involved: one for an intra-residue finding, two for an
        inter-residue one.
    residue_labels
        The same residues as :meth:`dodo.structure.Structure.residue_label` strings, e.g.
        ``"A:GLU142"``.
    atoms
        Atom names involved, in the same order as ``residue_indices`` for an inter-residue bond.
    measured
        The measured distance in Angstroms, or NaN when there was nothing to measure.
    expected, expected_sd
        Reference mean and standard deviation in Angstroms, or ``None`` where the check is not a
        comparison against a distribution (the CA-CA virtual bond carries no sd; a missing atom
        carries neither).
    accepted_low, accepted_high
        The window the measurement had to fall inside, in Angstroms. Asymmetric about
        ``expected`` for intra-residue bonds, deliberately -- see
        :data:`DEFAULT_TOLERANCE_SIGMA`.
    provenance
        See :data:`Provenance`. For the CA-CA check this selects the window; elsewhere it is
        advisory and never suppresses a finding.
    message
        A complete sentence naming the residues, the measured value and the expected range.
    """

    kind: BondViolationKind
    bond_class: BondClass
    residue_indices: tuple[int, ...]
    residue_labels: tuple[str, ...]
    atoms: tuple[str, ...]
    measured: float
    expected: float | None
    expected_sd: float | None
    accepted_low: float | None
    accepted_high: float | None
    provenance: Provenance
    message: str

    def __str__(self) -> str:
        return self.message

    @property
    def deviation(self) -> float:
        """Signed ``measured - expected`` in Angstroms, or NaN when either is undefined."""
        if self.expected is None or not math.isfinite(self.measured):
            return float("nan")
        return self.measured - self.expected

    @property
    def deviation_sigma(self) -> float:
        """Deviation in standard deviations of the reference distribution, or NaN."""
        if self.expected_sd is None or self.expected_sd <= 0.0:
            return float("nan")
        deviation = self.deviation
        return float("nan") if math.isnan(deviation) else abs(deviation) / self.expected_sd

    @property
    def overshoot(self) -> float:
        """How far outside the accepted window the measurement is, in Angstroms.

        The one severity measure comparable across every kind of finding, and the order
        :meth:`BondReport.by_severity` uses. Ranking by standard deviations instead would put a
        0.07 A error on a very tight bond above a 50 A chain break.

        One-sided checks contribute their one side: an unexpected contact must be at least
        ``BOND_CUTOFF`` apart and a chain break at most, so both still get a real number. Infinite
        only when there was nothing to measure at all -- an absent atom or a NaN coordinate -- so
        that those sort to the top of the list rather than to the bottom.
        """
        if not math.isfinite(self.measured):
            return float("inf")
        below = self.accepted_low - self.measured if self.accepted_low is not None else 0.0
        above = self.measured - self.accepted_high if self.accepted_high is not None else 0.0
        if self.accepted_low is None and self.accepted_high is None:
            return float("inf")
        return float(max(below, above, 0.0))


@dataclass(frozen=True, slots=True, eq=False)
class BondReport:
    """The result of :func:`validate_bonds`: counts, distributions and every finding.

    ``eq=False`` for the same reason as :class:`~dodo.geometry.metrics.TraceReport`: the
    generated ``__eq__`` would compare numpy arrays with ``==`` and raise "truth value of an
    array is ambiguous" the first time two reports were compared or put in a set.

    Attributes
    ----------
    violations
        Every finding, ordered by residue index then kind. Empty when the structure is clean.
    n_bonds_checked
        Total bonds whose length was compared against a reference.
    n_intra_checked, n_peptide_checked, n_ca_ca_checked
        Per-class breakdown of the above.
    intra_deviations, peptide_deviations, ca_ca_deviations
        Signed ``measured - reference mean`` in Angstroms for every checked bond of each class.
        Signed rather than absolute so a systematic bias is visible, which a mean absolute
        deviation would hide.
    n_residues, n_ca_only_residues
        Residues seen, and how many were CA-only and so exempt from the all-atom checks.
    n_partial_residues
        Residues with a backbone but no side chain past CB, exempt from the side-chain checks.
    n_unknown_residues
        Residues whose type is absent from the reference table (modified residues). Noted, not
        flagged.
    n_sequence_gaps
        Consecutive residue pairs that are not consecutive in the deposited numbering, whose
        inter-residue bonds were therefore not checked.
    n_ca_ca_skipped_peptide
        Pairs whose CA-CA virtual bond was skipped because a real peptide bond was measured
        instead. See the module docstring.
    notes
        Observations that are deliberately not findings: modified residues kept, chain breaks
        skipped, records out of order.
    tolerance_sigma, absolute_tolerance, ca_ca_tolerance
        The tolerances used, echoed so a report is self-describing.
    """

    violations: tuple[BondViolation, ...]
    n_bonds_checked: int
    n_intra_checked: int
    n_peptide_checked: int
    n_ca_ca_checked: int
    intra_deviations: np.ndarray
    peptide_deviations: np.ndarray
    ca_ca_deviations: np.ndarray
    n_residues: int
    n_ca_only_residues: int
    n_partial_residues: int
    n_unknown_residues: int
    n_sequence_gaps: int
    n_ca_ca_skipped_peptide: int
    notes: tuple[str, ...]
    tolerance_sigma: float
    absolute_tolerance: float
    ca_ca_tolerance: float

    # -- state ---------------------------------------------------------------

    @property
    def ok(self) -> bool:
        """True when nothing was violated. Notes do not make a report invalid.

        A ``"seam"`` finding does not count against ``ok``: it is the peptide bond at a
        rebuilt-region/folded-domain junction, strained because the folded neighbour was
        repositioned as a rigid body, not a defect DODO's rebuild introduced. It is still reported
        in :attr:`violations` (and surfaced on the run summary), just not treated as invalidating.
        """
        return not any(v.kind != "seam" for v in self.violations)

    def __bool__(self) -> bool:
        return self.ok

    def __repr__(self) -> str:
        state = "valid" if self.ok else f"{len(self.violations)} violations"
        return f"BondReport({self.n_bonds_checked} bonds checked, {state})"

    # -- slicing -------------------------------------------------------------

    def of_kind(self, kind: BondViolationKind) -> tuple[BondViolation, ...]:
        """All violations of one kind, in residue order."""
        return tuple(v for v in self.violations if v.kind == kind)

    def of_class(self, bond_class: BondClass) -> tuple[BondViolation, ...]:
        """All violations concerning one bond class, in residue order."""
        return tuple(v for v in self.violations if v.bond_class == bond_class)

    def of_provenance(self, provenance: Provenance) -> tuple[BondViolation, ...]:
        """All violations attributed to one provenance, in residue order.

        ``of_provenance("rebuilt")`` is the subset DODO is responsible for, and the one a
        regression test on rebuild output should assert is empty.
        """
        return tuple(v for v in self.violations if v.provenance == provenance)

    def by_severity(self) -> tuple[BondViolation, ...]:
        """Every violation, worst first, ranked by :attr:`BondViolation.overshoot`."""
        return tuple(
            sorted(self.violations, key=lambda v: (-v.overshoot, v.residue_indices[0], v.kind))
        )

    @property
    def worst(self) -> BondViolation | None:
        """The single worst violation by :attr:`BondViolation.overshoot`, or ``None`` if clean."""
        return self.by_severity()[0] if self.violations else None

    @property
    def violating_residues(self) -> tuple[int, ...]:
        """Residue indices involved in at least one violation, ascending and unique."""
        return tuple(sorted({i for v in self.violations for i in v.residue_indices}))

    # -- statistics ----------------------------------------------------------

    def distribution(self) -> dict[str, dict[str, float]]:
        """Absolute-deviation percentiles per bond class, in Angstroms.

        Returns
        -------
        dict
            ``{bond_class: {"n", "mean", "p50", "p99", "max"}}``, omitting any class with no
            measurements. This is the "is the whole structure slightly off, or is one residue
            broken?" view: a rebuild that stretched every bond a little shows up as a shifted
            signed mean, not as a long violation list.
        """
        out: dict[str, dict[str, float]] = {}
        for name, values in (
            ("intra_residue", self.intra_deviations),
            ("peptide", self.peptide_deviations),
            ("ca_ca", self.ca_ca_deviations),
        ):
            finite = values[np.isfinite(values)]
            if finite.size == 0:
                continue
            absolute = np.abs(finite)
            out[name] = {
                "n": float(finite.size),
                "mean": float(finite.mean()),
                "p50": float(np.percentile(absolute, 50)),
                "p99": float(np.percentile(absolute, 99)),
                "max": float(absolute.max()),
            }
        return out

    def summary(self) -> str:
        """One-line summary of what was checked and what came of it."""
        parts = [
            f"{self.n_bonds_checked} bonds checked "
            f"({self.n_intra_checked} intra-residue, {self.n_peptide_checked} peptide, "
            f"{self.n_ca_ca_checked} CA-CA)",
            f"{self.n_residues} residues ({self.n_ca_only_residues} CA-only)",
        ]
        worst = self.worst
        if worst is not None:
            parts.append(
                f"worst {worst.kind} at {'/'.join(worst.residue_labels)}"
                + (
                    f" ({worst.measured:.3f} A, {worst.overshoot:.3f} A outside range)"
                    if math.isfinite(worst.overshoot)
                    else ""
                )
            )
        parts.append("valid" if self.ok else f"{len(self.violations)} violations")
        return ", ".join(parts)

    def describe(self, max_violations: int = 10) -> str:
        """Full multi-line description: summary, deviation distribution, then the findings.

        Parameters
        ----------
        max_violations
            Truncate the list after this many, worst first, with a count of the remainder. A
            badly broken 8,000-residue model produces thousands of findings and an unreadable
            message; the worst few are what identifies the failure.
        """
        lines = [self.summary()]
        for name, stats in self.distribution().items():
            lines.append(
                f"  {name}: |deviation| median {stats['p50']:.3f} A, "
                f"p99 {stats['p99']:.3f} A, max {stats['max']:.3f} A "
                f"(mean signed {stats['mean']:+.4f} A over {int(stats['n'])} bonds)"
            )
        counts = Counter(v.kind for v in self.violations)
        if counts:
            lines.append(
                "  violations by kind: "
                + ", ".join(f"{kind} {n}" for kind, n in sorted(counts.items()))
            )
        for violation in self.by_severity()[:max_violations]:
            lines.append(f"  - {violation.message}")
        remaining = len(self.violations) - max_violations
        if remaining > 0:
            lines.append(f"  - ... and {remaining} more violation(s)")
        for note in self.notes:
            lines.append(f"  note: {note}")
        return "\n".join(lines)

    def raise_if_invalid(self) -> None:
        """Raise :class:`~dodo.exceptions.GeometryError` if anything was violated."""
        if not self.ok:
            raise GeometryError(f"Invalid bond geometry: {self.describe()}")


# ---------------------------------------------------------------------------
# Reference tables, flattened once at import
# ---------------------------------------------------------------------------

# Completeness tiers, used only to decide whether an absent atom is a defect or a deliberate
# partial build. A residue is assigned the lowest tier that contains every atom it actually has,
# and an atom is reported missing only if its own tier is at or below the residue's.
#
# This exists because "rebuilt regions are CA-only" is only the DEFAULT state, not the only one:
# ``rebuild(backbone=True)`` places real N, C and O in rebuilt regions, and side-chain placement is
# the next priority after that. A residue with N/CA/C/O and no side chain, or one with a CB and
# nothing beyond it, is therefore a deliberate output state.
# Measured: without this rule a backbone-only version of a rebuilt dnmt3a produces 2,323
# missing-atom findings, every one of them noise. Truncated-at-CB side chains in deposited
# structures -- the usual way side-chain disorder is modelled -- are exempted for the same reason
# and counted in a note.
_TIER_BACKBONE: Final[int] = 0
_TIER_CB: Final[int] = 1
_TIER_SIDECHAIN: Final[int] = 2
_BACKBONE_ATOM_NAMES: Final[frozenset[str]] = frozenset({"N", "CA", "C", "O", "OXT"})


@dataclass(frozen=True, slots=True)
class _TypeBonds:
    """One residue type's reference bonds, as parallel arrays over vocabulary columns."""

    col_a: np.ndarray
    col_b: np.ndarray
    name_a: tuple[str, ...]
    name_b: tuple[str, ...]
    mean: np.ndarray
    sd: np.ndarray
    #: Measured 0.1 and 99.9 percentiles, which the accepted window is anchored on.
    p_low: np.ndarray
    p_high: np.ndarray
    #: Vocabulary columns of every atom this residue type must have, and their tiers.
    expected_cols: np.ndarray
    expected_names: tuple[str, ...]
    expected_tier: np.ndarray


def _build_vocabulary() -> tuple[np.ndarray, dict[str, int]]:
    """Sorted array of every atom name the reference knows, plus a name-to-column map."""
    names: set[str] = set()
    for table in RESIDUE_BONDS.values():
        for pair in table:
            names.update(pair)
    ordered = sorted(names)
    return np.array(ordered, dtype="<U4"), {name: i for i, name in enumerate(ordered)}


_ATOM_VOCAB, _ATOM_COL = _build_vocabulary()
_RESIDUE_TYPES: Final[np.ndarray] = np.array(sorted(RESIDUE_BONDS), dtype="<U3")
_V: Final[int] = len(_ATOM_VOCAB)
_CA_COL: Final[int] = _ATOM_COL["CA"]
_N_COL: Final[int] = _ATOM_COL["N"]
_C_COL: Final[int] = _ATOM_COL["C"]
_BACKBONE_COLS: Final[np.ndarray] = np.array(
    sorted(_ATOM_COL[a] for a in _BACKBONE_ATOM_NAMES if a in _ATOM_COL), dtype=np.int64
)
_BACKBONE_CB_COLS: Final[np.ndarray] = np.array(
    sorted([*_BACKBONE_COLS.tolist(), _ATOM_COL["CB"]]), dtype=np.int64
)


def _atom_tier(atom: str) -> int:
    """Completeness tier of one atom name. See :data:`_TIER_BACKBONE`."""
    if atom in _BACKBONE_ATOM_NAMES:
        return _TIER_BACKBONE
    return _TIER_CB if atom == "CB" else _TIER_SIDECHAIN


def _flatten(
    table: Mapping[tuple[str, str], tuple[float, float, int, float, float]],
) -> _TypeBonds:
    """Turn one residue type's bond dict into arrays keyed by vocabulary column."""
    pairs = sorted(table)
    required = sorted({atom for pair in pairs for atom in pair} - _OPTIONAL_ATOMS)
    return _TypeBonds(
        col_a=np.array([_ATOM_COL[a] for a, _ in pairs], dtype=np.int64),
        col_b=np.array([_ATOM_COL[b] for _, b in pairs], dtype=np.int64),
        name_a=tuple(a for a, _ in pairs),
        name_b=tuple(b for _, b in pairs),
        mean=np.array([table[p][0] for p in pairs], dtype=np.float64),
        sd=np.array([table[p][1] for p in pairs], dtype=np.float64),
        p_low=np.array([table[p][3] for p in pairs], dtype=np.float64),
        p_high=np.array([table[p][4] for p in pairs], dtype=np.float64),
        expected_cols=np.array([_ATOM_COL[a] for a in required], dtype=np.int64),
        expected_names=tuple(required),
        expected_tier=np.array([_atom_tier(a) for a in required], dtype=np.int64),
    )


_TYPE_BONDS: Final[dict[str, _TypeBonds]] = {
    name: _flatten(table) for name, table in RESIDUE_BONDS.items()
}

#: Highest completeness tier each residue type actually reaches, in the same order as
#: :data:`_RESIDUE_TYPES`. A residue is only "partially built" if it sits below this: GLY is
#: complete with a bare backbone and ALA is complete at CB, so counting either as partial would
#: report 1,137 of 6kn7's 7,781 residues as missing side chains they never had.
_TYPE_REQUIRED_TIER: Final[np.ndarray] = np.array(
    [int(_TYPE_BONDS[str(name)].expected_tier.max()) for name in _RESIDUE_TYPES], dtype=np.int64
)


def _allowed_contact_keys() -> np.ndarray:
    """Sorted keys of every (residue type, atom pair) the reference accepts as bonded.

    Encoded as ``type_code * V * V + low_col * V + high_col`` so the unexpected-contact check is
    a single vectorized :func:`numpy.isin` rather than a dict lookup per close pair.
    """
    keys: list[int] = []
    for type_code, residue in enumerate(_RESIDUE_TYPES):
        for atom_a, atom_b in RESIDUE_BONDS[str(residue)]:
            low, high = sorted((_ATOM_COL[atom_a], _ATOM_COL[atom_b]))
            keys.append(type_code * _V * _V + low * _V + high)
    return np.array(sorted(keys), dtype=np.int64)


_ALLOWED_CONTACT_KEYS: Final[np.ndarray] = _allowed_contact_keys()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _encode(values: np.ndarray, vocabulary: np.ndarray) -> np.ndarray:
    """Map each string in ``values`` to its index in sorted ``vocabulary``, or -1.

    Vectorized via searchsorted. The alternative -- a dict lookup per atom -- is a Python loop
    over the 61,511 atoms this has to stay fast on.
    """
    position = np.searchsorted(vocabulary, values)
    clipped = np.clip(position, 0, vocabulary.size - 1)
    hit = vocabulary[clipped] == values
    return np.where(hit, clipped, -1).astype(np.int64)


_PROVENANCE_NAMES: Final[tuple[Provenance, Provenance, Provenance]] = (
    "unknown",
    "input",
    "rebuilt",
)


#: Per-residue provenance codes produced by :func:`_residue_provenance`.
_UNKNOWN_CODE, _INPUT_CODE, _GENERATED_CODE = 0, 1, 2


def _residue_provenance(structure: Structure) -> np.ndarray:
    """Per-residue provenance codes: 0 unknown, 1 input, 2 generated by DODO.

    Code 2 means DODO produced these coordinates and is answerable for their geometry, so it is
    keyed on :meth:`~dodo.structure.Domain.generated_spans` rather than on region kind.

    Keying it on kind, as this once did, mislabelled two whole categories as DODO's work: folded
    domains are only translated and rotated, and a region whose build *failed* keeps its input
    coordinates. On AF-Q9BTC0-F1 that attributed six of AlphaFold's own short CA-CA bonds and a
    0.634 A atom pair to DODO, all of them present in the input file at identical values.
    """
    codes = np.zeros(structure.n_residues, dtype=np.int8)
    for domain in structure.domains:
        # Everything inside a recognised region is at least known to be input geometry.
        codes[domain.span.start : domain.span.stop] = _INPUT_CODE
    for domain in structure.domains:
        for span in domain.generated_spans():
            codes[span.start : span.stop] = _GENERATED_CODE
    return codes


def _validated_tolerance(name: str, value: float, *, allow_zero: bool = True) -> float:
    """Coerce a tolerance to float and reject values that would silently disable a check."""
    number = float(value)
    # A NaN threshold makes every comparison below False, which turns a check off instead of
    # failing -- the exact failure mode geometry/metrics.py guards against for the same reason.
    if not math.isfinite(number):
        raise GeometryError(f"{name} must be finite, got {number}.")
    if number < 0.0 or (number == 0.0 and not allow_zero):
        raise GeometryError(f"{name} must be positive, got {number}.")
    return number


# ---------------------------------------------------------------------------
# Shared state for the checks
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class _Context:
    """Everything the individual checks need, computed once."""

    structure: Structure
    xyz: np.ndarray
    finite_atom: np.ndarray
    #: ``(n_residues, n_vocabulary)`` atom index per residue and reference atom name, or -1.
    atom_slot: np.ndarray
    type_code: np.ndarray
    ca_only: np.ndarray
    residue_tier: np.ndarray
    provenance: np.ndarray
    tolerance_sigma: float
    absolute_tolerance: float
    violations: list[BondViolation]
    notes: list[str]

    def label(self, residue: int) -> str:
        return self.structure.residue_label(residue)

    def provenance_of(self, *residues: int) -> Provenance:
        """Provenance of a finding: rebuilt wins over input, which wins over unknown."""
        return _PROVENANCE_NAMES[max(int(self.provenance[r]) for r in residues)]


def _prepare(
    structure: Structure, *, tolerance_sigma: float, absolute_tolerance: float
) -> _Context:
    """Build the per-atom and per-residue index tables every check reads."""
    n_residues = structure.n_residues
    xyz = structure.xyz
    finite_atom = np.all(np.isfinite(xyz), axis=1)

    atom_col = _encode(structure.atom_name, _ATOM_VOCAB)
    # A non-finite atom is deliberately left out of the slot table, so no check downstream
    # measures a distance to NaN and reports a meaningless number. The coordinate itself is
    # reported once, by _check_finite.
    usable = (atom_col >= 0) & finite_atom
    atom_slot = np.full((n_residues, _V), -1, dtype=np.int64)
    atom_slot[structure.residue_index[usable], atom_col[usable]] = np.flatnonzero(usable)

    atoms_per_residue = np.diff(structure.residue_atom_offsets)
    # "CA-only" is the discriminator between a rebuilt region and an incomplete residue, and it
    # has to be exactly this: one atom, and it is the alpha carbon. A residue with one atom that
    # is not a CA is an incomplete residue, not a rebuilt one, and is reported as such.
    ca_only = (atoms_per_residue == 1) & (atom_slot[:, _CA_COL] >= 0)

    # Completeness tier: the lowest tier containing every atom the residue actually has. Atoms
    # outside the reference vocabulary count as side-chain atoms, so an exotic atom pushes the
    # residue to the strictest tier rather than relaxing it.
    n_outside_backbone = np.bincount(
        structure.residue_index[~np.isin(atom_col, _BACKBONE_COLS)], minlength=n_residues
    )
    n_outside_cb = np.bincount(
        structure.residue_index[~np.isin(atom_col, _BACKBONE_CB_COLS)], minlength=n_residues
    )
    residue_tier = np.where(
        n_outside_backbone == 0,
        _TIER_BACKBONE,
        np.where(n_outside_cb == 0, _TIER_CB, _TIER_SIDECHAIN),
    ).astype(np.int64)

    return _Context(
        structure=structure,
        xyz=xyz,
        finite_atom=finite_atom,
        atom_slot=atom_slot,
        type_code=_encode(structure.residue_name, _RESIDUE_TYPES),
        ca_only=ca_only,
        residue_tier=residue_tier,
        provenance=_residue_provenance(structure),
        tolerance_sigma=tolerance_sigma,
        absolute_tolerance=absolute_tolerance,
        violations=[],
        notes=[],
    )


# ---------------------------------------------------------------------------
# The checks
# ---------------------------------------------------------------------------


def _check_finite(context: _Context) -> None:
    """Report every atom whose coordinate is not finite."""
    structure = context.structure
    for atom in np.flatnonzero(~context.finite_atom).tolist():
        residue = int(structure.residue_index[atom])
        name = str(structure.atom_name[atom])
        label = context.label(residue)
        context.violations.append(
            BondViolation(
                kind="non_finite",
                bond_class="none",
                residue_indices=(residue,),
                residue_labels=(label,),
                atoms=(name,),
                measured=float("nan"),
                expected=None,
                expected_sd=None,
                accepted_low=None,
                accepted_high=None,
                provenance=context.provenance_of(residue),
                message=(
                    f"Atom {name} of residue {label} has a non-finite coordinate "
                    f"{structure.xyz[atom].tolist()}, so no bond involving it could be measured; "
                    f"a NaN coordinate means an upstream sampler returned its failure instead of "
                    f"raising."
                ),
            )
        )


def _add_length_violation(
    context: _Context,
    *,
    bond_class: BondClass,
    residues: tuple[int, ...],
    atoms: tuple[str, ...],
    measured: float,
    mean: float,
    sd: float | None,
    low: float,
    high: float,
) -> None:
    """Append a ``bond_length`` finding, phrased as a complete sentence."""
    labels = tuple(context.label(r) for r in residues)
    if len(residues) == 1:
        where = f"in residue {labels[0]}"
        bond = f"{atoms[0]}-{atoms[1]}"
    else:
        where = f"between residues {labels[0]} and {labels[1]}"
        bond = f"{atoms[0]}({labels[0]})-{atoms[1]}({labels[1]})"
    spread = (
        f", which is {abs(measured - mean) / sd:.1f} standard deviations from the reference mean"
        if sd
        else ""
    )
    context.violations.append(
        BondViolation(
            kind="bond_length",
            bond_class=bond_class,
            residue_indices=residues,
            residue_labels=labels,
            atoms=atoms,
            measured=measured,
            expected=mean,
            expected_sd=sd,
            accepted_low=low,
            accepted_high=high,
            provenance=context.provenance_of(*residues),
            message=(
                f"The {bond} bond {where} is {measured:.3f} A, outside the accepted "
                f"{low:.3f}-{high:.3f} A around a reference mean of {mean:.3f} A{spread}."
            ),
        )
    )


def _check_intra_residue(context: _Context) -> tuple[int, np.ndarray]:
    """Check every intra-residue bond, and report atoms that should be there and are not.

    Returns the number of bonds checked and their signed deviations from the reference mean.
    """
    deviations: list[np.ndarray] = []
    n_checked = 0
    sigma = context.tolerance_sigma
    floor = context.absolute_tolerance

    for type_code, residue_name in enumerate(_RESIDUE_TYPES):
        selected = np.flatnonzero(context.type_code == type_code)
        if selected.size == 0:
            continue
        bonds = _TYPE_BONDS[str(residue_name)]
        # One vectorized pass per residue TYPE -- 20 for any structure, rather than one per
        # residue. On 6kn7 that is 20 passes instead of 7,781.
        index_a = context.atom_slot[selected[:, None], bonds.col_a[None, :]]
        index_b = context.atom_slot[selected[:, None], bonds.col_b[None, :]]
        rows, columns = np.nonzero((index_a >= 0) & (index_b >= 0))
        if rows.size:
            measured = np.linalg.norm(
                context.xyz[index_a[rows, columns]] - context.xyz[index_b[rows, columns]], axis=1
            )
            mean = bonds.mean[columns]
            sd = bonds.sd[columns]
            slack = np.maximum(sigma * sd, floor)
            low = bonds.p_low[columns] - slack
            high = bonds.p_high[columns] + slack
            n_checked += int(measured.size)
            deviations.append(measured - mean)
            for hit in np.flatnonzero((measured < low) | (measured > high)).tolist():
                column = int(columns[hit])
                _add_length_violation(
                    context,
                    bond_class="intra_residue",
                    residues=(int(selected[rows[hit]]),),
                    atoms=(bonds.name_a[column], bonds.name_b[column]),
                    measured=float(measured[hit]),
                    mean=float(mean[hit]),
                    sd=float(sd[hit]),
                    low=float(low[hit]),
                    high=float(high[hit]),
                )

        _report_missing_atoms(context, selected, bonds)

    stacked = np.concatenate(deviations) if deviations else np.zeros(0, dtype=np.float64)
    return n_checked, stacked


def _report_missing_atoms(context: _Context, selected: np.ndarray, bonds: _TypeBonds) -> None:
    """Report atoms a residue type requires that are absent, respecting completeness tiers.

    Skipping CA-only residues is the single most important exclusion in this module. A rebuilt
    region is alpha carbons and nothing else, by design, so demanding N/C/O there would produce
    one finding per rebuilt residue -- 334 on a rebuilt dnmt3a -- and every one would be noise.
    The tier rule generalizes that to partially built residues; see :data:`_TIER_BACKBONE`. What
    remains reportable is a residue that has side-chain atoms and is still missing one, which is
    an incomplete residue however it arose.
    """
    all_atom = selected[~context.ca_only[selected]]
    if all_atom.size == 0:
        return
    slots = context.atom_slot[all_atom[:, None], bonds.expected_cols[None, :]]
    in_scope = bonds.expected_tier[None, :] <= context.residue_tier[all_atom][:, None]
    rows, columns = np.nonzero((slots < 0) & in_scope)
    for row, column in zip(rows.tolist(), columns.tolist(), strict=True):
        residue = int(all_atom[row])
        atom = bonds.expected_names[column]
        label = context.label(residue)
        n_bonds = sum(1 for a, b in zip(bonds.name_a, bonds.name_b, strict=True) if atom in (a, b))
        residue_name = str(context.structure.residue_name[residue])
        context.violations.append(
            BondViolation(
                kind="missing_atom",
                bond_class="none",
                residue_indices=(residue,),
                residue_labels=(label,),
                atoms=(atom,),
                measured=float("nan"),
                expected=None,
                expected_sd=None,
                accepted_low=None,
                accepted_high=None,
                provenance=context.provenance_of(residue),
                message=(
                    f"Residue {label} has atoms other than CA but no {atom}, so the {n_bonds} "
                    f"reference bond(s) of {residue_name} involving {atom} could not be checked; "
                    f"this residue is incomplete rather than deliberately CA-only."
                ),
            )
        )


def _consecutive_pairs(structure: Structure) -> tuple[np.ndarray, np.ndarray]:
    """Residue indices ``i`` whose successor ``i + 1`` is the next residue of the same chain.

    Returns the contiguously numbered pairs and the ones separated by a numbering gap. Splitting
    them is what keeps an unmodelled-residue chain break from being reported as a broken bond:
    two residues numbered 380 and 393 are not chemically adjacent, so no bond between them is
    expected and none is checked.
    """
    if structure.n_residues < 2:
        empty = np.zeros(0, dtype=np.int64)
        return empty, empty
    first = np.arange(structure.n_residues - 1, dtype=np.int64)
    same_chain = structure.chain_index[:-1] == structure.chain_index[1:]
    delta = structure.residue_number[1:] - structure.residue_number[:-1]
    # Insertion codes make 10 and 10A consecutive in the polymer at the same residue number, so
    # a zero delta with different codes is contiguous too. Without that clause every insertion
    # would be reported as a chain break.
    different_icode = structure.insertion_code[1:] != structure.insertion_code[:-1]
    contiguous = (delta == 1) | ((delta == 0) & different_icode)
    return first[same_chain & contiguous], first[same_chain & ~contiguous]


def _check_peptide_bonds(context: _Context, pairs: np.ndarray) -> tuple[int, np.ndarray]:
    """Check C(i)-N(i+1) for every consecutively numbered pair that has both atoms."""
    tolerance = max(context.tolerance_sigma * PEPTIDE_BOND_SD, context.absolute_tolerance)
    low = PEPTIDE_BOND_LENGTH - tolerance
    high = PEPTIDE_BOND_LENGTH + tolerance

    c_index = context.atom_slot[pairs, _C_COL]
    n_index = context.atom_slot[pairs + 1, _N_COL]
    rows = np.flatnonzero((c_index >= 0) & (n_index >= 0))
    if rows.size == 0:
        return 0, np.zeros(0, dtype=np.float64)

    measured = np.linalg.norm(context.xyz[c_index[rows]] - context.xyz[n_index[rows]], axis=1)
    for hit in np.flatnonzero((measured < low) | (measured > high)).tolist():
        first = int(pairs[rows[hit]])
        distance = float(measured[hit])
        # A C-N pair beyond the bonding cutoff is not a stretched bond, it is a chain that is not
        # connected at all -- a different defect with a different fix, so it gets its own kind.
        if distance > BOND_CUTOFF:
            _add_chain_break(context, first, distance)
        else:
            _add_length_violation(
                context,
                bond_class="peptide",
                residues=(first, first + 1),
                atoms=("C", "N"),
                measured=distance,
                mean=PEPTIDE_BOND_LENGTH,
                sd=PEPTIDE_BOND_SD,
                low=low,
                high=high,
            )
    return int(measured.size), measured - PEPTIDE_BOND_LENGTH


def _add_chain_break(context: _Context, first: int, distance: float) -> None:
    """Report a C-N pair too far apart to be a bond, between consecutively numbered residues.

    A pair that straddles a generated<->input boundary is not a chain break but a SEAM: the rebuilt
    region abuts a folded domain whose boundary nitrogen was moved as a rigid body in step 3 and
    now points where AlphaFold's original chain ran. That junction is a deliberately-placed,
    CONECT-drawn continuous bond whose strain is inherited from the untouched folded atom rather
    than caused by the rebuild, so it is classified ``kind="seam"`` / ``provenance="seam"``:
    reported, but not counted against :attr:`BondReport.ok` and not blamed on ``"rebuilt"``. Keyed
    strictly on the generated<->input boundary, so a genuine chain break in non-DODO input (which
    has no generated residues) is never masked.
    """
    labels = (context.label(first), context.label(first + 1))
    codes = {int(context.provenance[first]), int(context.provenance[first + 1])}
    if codes == {_INPUT_CODE, _GENERATED_CODE}:
        context.violations.append(
            BondViolation(
                kind="seam",
                bond_class="peptide",
                residue_indices=(first, first + 1),
                residue_labels=labels,
                atoms=("C", "N"),
                measured=distance,
                expected=PEPTIDE_BOND_LENGTH,
                expected_sd=PEPTIDE_BOND_SD,
                accepted_low=None,
                accepted_high=BOND_CUTOFF,
                provenance="seam",
                message=(
                    f"Residues {labels[0]} and {labels[1]} span a rebuilt-region/folded-domain "
                    f"seam; their C-N distance is {distance:.3f} A (ideal "
                    f"{PEPTIDE_BOND_LENGTH:.3f}) because the folded neighbour was repositioned as "
                    f"a rigid body and its nitrogen cannot be reached exactly. The bond is drawn "
                    f"and reported, not a defect DODO introduced."
                ),
            )
        )
        return
    context.violations.append(
        BondViolation(
            kind="chain_break",
            bond_class="peptide",
            residue_indices=(first, first + 1),
            residue_labels=labels,
            atoms=("C", "N"),
            measured=distance,
            expected=PEPTIDE_BOND_LENGTH,
            expected_sd=PEPTIDE_BOND_SD,
            accepted_low=None,
            accepted_high=BOND_CUTOFF,
            provenance=context.provenance_of(first, first + 1),
            message=(
                f"Residues {labels[0]} and {labels[1]} are numbered consecutively but their C-N "
                f"distance is {distance:.3f} A, beyond the {BOND_CUTOFF:.2f} A bonding cutoff, so "
                f"the chain is broken there rather than merely strained (a peptide bond is "
                f"{PEPTIDE_BOND_LENGTH:.3f} A)."
            ),
        )
    )


def _check_ca_ca_bonds(
    context: _Context, pairs: np.ndarray, *, tolerance: float, strict: bool
) -> tuple[int, np.ndarray, int]:
    """Check the CA-CA virtual bond for consecutive residues.

    Skipped wherever a real peptide bond was measured, unless ``strict``. See the module
    docstring: the virtual bond is a surrogate for connectivity, and where the covalent bond
    itself is available the surrogate only disagrees with reality.
    """
    ca_index = context.atom_slot[:, _CA_COL]
    has_peptide = (context.atom_slot[pairs, _C_COL] >= 0) & (
        context.atom_slot[pairs + 1, _N_COL] >= 0
    )
    checkable = (ca_index[pairs] >= 0) & (ca_index[pairs + 1] >= 0)
    skipped = int(np.count_nonzero(checkable & has_peptide))
    selected = pairs[checkable] if strict else pairs[checkable & ~has_peptide]
    if selected.size == 0:
        return 0, np.zeros(0, dtype=np.float64), skipped

    measured = np.linalg.norm(
        context.xyz[ca_index[selected]] - context.xyz[ca_index[selected + 1]], axis=1
    )

    # Two windows, chosen per pair by provenance. Geometry DODO GENERATED is held to the value it
    # builds to. Geometry DODO merely PRESERVED is held to the union of both measured reference
    # populations -- trans at 3.8317 and the cis/low-confidence population at 3.1912 -- because a
    # cis peptide in the input is not a defect and, with no C or N atom present, nothing else can
    # prove it is one. Measured on the AlphaFold fixtures' CA traces, the strict window flags
    # 6.8-22.8% of preserved bonds and this one flags 0-0.7%.
    rebuilt = (context.provenance[selected] == 2) | (context.provenance[selected + 1] == 2)
    low = np.where(rebuilt, CA_CA_BOND_LENGTH - tolerance, CA_CA_SHORT_LENGTH - tolerance)
    high = np.where(rebuilt, CA_CA_BOND_LENGTH + tolerance, CA_CA_TRANS_LENGTH + tolerance)

    for hit in np.flatnonzero((measured < low) | (measured > high)).tolist():
        first = int(selected[hit])
        distance = float(measured[hit])
        labels = (context.label(first), context.label(first + 1))
        detail = ""
        if distance < _COINCIDENT:
            detail = (
                " -- they are effectively the same point, which is what a region grafted onto "
                "the wrong anchor looks like"
            )
        elif not rebuilt[hit]:
            detail = (
                " -- neither residue was rebuilt, so this geometry came from the input and the "
                "window used spans both reference populations"
            )
        context.violations.append(
            BondViolation(
                kind="bond_length",
                bond_class="ca_ca",
                residue_indices=(first, first + 1),
                residue_labels=labels,
                atoms=("CA", "CA"),
                measured=distance,
                expected=CA_CA_BOND_LENGTH,
                expected_sd=None,
                accepted_low=float(low[hit]),
                accepted_high=float(high[hit]),
                provenance=context.provenance_of(first, first + 1),
                message=(
                    f"The virtual CA-CA bond between residues {labels[0]} and {labels[1]} is "
                    f"{distance:.3f} A, outside the accepted {low[hit]:.3f}-{high[hit]:.3f} A "
                    f"for a bond DODO builds at {CA_CA_BOND_LENGTH:.2f} A{detail}."
                ),
            )
        )
    return int(measured.size), measured - CA_CA_BOND_LENGTH, skipped


def _check_unexpected_contacts(context: _Context) -> None:
    """Report intra-residue atom pairs closer than BOND_CUTOFF that are not bonds.

    Within a residue the separation is unambiguous: measured over 3,066,489 intra-residue pairs,
    covalent bonds end at 1.85 A and 1-3 contacts begin at 2.15 A. So an unexpected pair below
    1.90 A is not a tight angle, it is two atoms placed on top of each other.

    Restricted to intra-residue pairs on purpose. Cross-residue coincidence is a steric clash and
    belongs to the clash validator; the only cross-residue pair inside 1.90 A in a real structure
    is the peptide bond itself (measured on 6kn7: 7,750 of them and nothing else).
    """
    structure = context.structure
    finite = np.flatnonzero(context.finite_atom)
    if finite.size < 2:
        return
    # A KD-tree, not a nested loop: on the 61,511-atom assembly this returns all 62,548 pairs
    # inside 1.90 A in 15 ms, where the pairwise alternative is 1.9 billion comparisons.
    tree = cKDTree(context.xyz[finite])
    pairs = tree.query_pairs(BOND_CUTOFF, output_type="ndarray")
    if pairs.size == 0:
        return
    atom_a = finite[pairs[:, 0]]
    atom_b = finite[pairs[:, 1]]

    residue = structure.residue_index[atom_a]
    same_residue = residue == structure.residue_index[atom_b]
    hydrogens = list(_HYDROGEN_ELEMENTS)
    heavy = ~np.isin(structure.element[atom_a], hydrogens) & ~np.isin(
        structure.element[atom_b], hydrogens
    )
    # A residue type absent from the reference carries no expectation, so nothing about it can be
    # "unexpected". Those are noted instead, by _note_unknown_residues.
    known_type = context.type_code[residue] >= 0
    candidate = np.flatnonzero(same_residue & heavy & known_type)
    if candidate.size == 0:
        return

    # An atom name outside the vocabulary is mapped to column _V, which cannot match any allowed
    # key, so an unrecognized atom sitting on top of a known one is still reported.
    encoded_a = _encode(structure.atom_name[atom_a[candidate]], _ATOM_VOCAB)
    encoded_b = _encode(structure.atom_name[atom_b[candidate]], _ATOM_VOCAB)
    col_a = np.where(encoded_a < 0, _V, encoded_a)
    col_b = np.where(encoded_b < 0, _V, encoded_b)
    keys = (
        context.type_code[residue[candidate]] * _V * _V
        + np.minimum(col_a, col_b) * _V
        + np.maximum(col_a, col_b)
    )
    for index in candidate[~np.isin(keys, _ALLOWED_CONTACT_KEYS)].tolist():
        i = int(atom_a[index])
        j = int(atom_b[index])
        residue_index = int(structure.residue_index[i])
        label = context.label(residue_index)
        name_i = str(structure.atom_name[i])
        name_j = str(structure.atom_name[j])
        distance = float(np.linalg.norm(context.xyz[i] - context.xyz[j]))
        context.violations.append(
            BondViolation(
                kind="unexpected_contact",
                bond_class="intra_residue",
                residue_indices=(residue_index,),
                residue_labels=(label,),
                atoms=(name_i, name_j),
                measured=distance,
                expected=None,
                expected_sd=None,
                accepted_low=BOND_CUTOFF,
                accepted_high=None,
                provenance=context.provenance_of(residue_index),
                message=(
                    f"Atoms {name_i} and {name_j} of residue {label} are {distance:.3f} A apart, "
                    f"inside the {BOND_CUTOFF:.2f} A bonding cutoff, but "
                    f"{structure.residue_name[residue_index]} has no {name_i}-{name_j} bond; two "
                    f"unbonded atoms this close have been placed on top of each other."
                ),
            )
        )


def _note_unknown_residues(context: _Context) -> int:
    """Note residue types absent from the reference table, and return how many there are."""
    unknown = np.flatnonzero(context.type_code < 0)
    if unknown.size == 0:
        return 0
    counts = Counter(str(name) for name in context.structure.residue_name[unknown])
    listing = ", ".join(f"{name} x{n}" for name, n in sorted(counts.items()))
    examples = ", ".join(context.label(int(i)) for i in unknown[:3])
    context.notes.append(
        f"{unknown.size} residue(s) of {len(counts)} type(s) are not in the reference bond table, "
        f"so their intra-residue geometry was not checked: {listing} (e.g. {examples}). DODO "
        f"keeps modified residues on purpose -- dropping MSE fabricates a phantom chain break -- "
        f"so this is noted, not flagged."
    )
    return int(unknown.size)


def _note_sequence_gaps(context: _Context, gaps: np.ndarray) -> None:
    """Note consecutive residues that are not consecutive in the deposited numbering."""
    if gaps.size == 0:
        return
    structure = context.structure
    backwards = int(
        np.count_nonzero(structure.residue_number[gaps + 1] < structure.residue_number[gaps])
    )
    examples = ", ".join(
        f"{context.label(int(i))} to {context.label(int(i) + 1)}" for i in gaps[:3]
    )
    context.notes.append(
        f"{gaps.size} consecutive residue pair(s) are not consecutive in the deposited numbering, "
        f"so no peptide or CA-CA bond was expected or checked between them ({examples}). "
        f"Unmodelled residues are the usual cause and are not a defect."
    )
    if backwards:
        context.notes.append(
            f"{backwards} of those pair(s) run backwards in residue number, which means the atom "
            f"records are out of chain order rather than merely interrupted."
        )


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def validate_bonds(
    structure: Structure,
    *,
    tolerance_sigma: float = DEFAULT_TOLERANCE_SIGMA,
    absolute_tolerance: float = DEFAULT_ABSOLUTE_TOLERANCE,
    ca_ca_tolerance: float = DEFAULT_CA_CA_TOLERANCE,
    check_intra_residue: bool = True,
    check_peptide: bool = True,
    check_ca_ca: bool = True,
    check_unexpected_contacts: bool = True,
    strict_ca_ca: bool = False,
) -> BondReport:
    """Check every bond in ``structure`` against measured reference geometry.

    Parameters
    ----------
    structure
        The structure to validate. Mixed all-atom and CA-only regions are expected and handled;
        see the module docstring.
    tolerance_sigma
        Slack in standard deviations of that specific bond's measured distribution, added outside
        its measured 0.1-99.9 percentile range. Defaults to :data:`DEFAULT_TOLERANCE_SIGMA`, which
        flags 0.35% of the bonds in this repository's real fixtures.
    absolute_tolerance
        Floor on the slack in Angstroms, so a bond with a very tight measured sd does not become
        unreasonably strict. Defaults to :data:`DEFAULT_ABSOLUTE_TOLERANCE`.
    ca_ca_tolerance
        Allowed deviation of the CA-CA virtual bond in Angstroms, from
        :data:`dodo.constants.CA_CA_BOND_LENGTH` for rebuilt geometry and from the reference
        populations for preserved geometry. Defaults to :data:`DEFAULT_CA_CA_TOLERANCE`. Absolute
        because the reference carries no per-bond sd for a virtual bond.
    check_intra_residue, check_peptide, check_ca_ca, check_unexpected_contacts
        Switch individual checks off. All on by default; they exist so a caller mid-build can run
        just the CA-CA gate cheaply.
    strict_ca_ca
        Also check the CA-CA virtual bond where a real peptide bond was measured. Off by default:
        cis and strained peptides in real input legitimately put CA-CA between 2.9 and 4.4 A, and
        checking it anyway reports 26% of p300's residue pairs while every one of its peptide
        bonds is within tolerance.

    Returns
    -------
    BondReport
        Counts, deviation distributions and every finding. Truthy when nothing was violated.

    Raises
    ------
    GeometryError
        Only if the input cannot be validated at all: not a :class:`~dodo.structure.Structure`,
        empty, or a tolerance that is negative or not finite. Bad *geometry* is reported, never
        raised -- that is this function's job. Call :meth:`BondReport.raise_if_invalid` for the
        strict form.

    Examples
    --------
    >>> import dodo
    >>> from dodo.validate.bonds import validate_bonds
    >>> structure = dodo.read_structure("6kn7.pdb")          # doctest: +SKIP
    >>> report = validate_bonds(structure)                   # doctest: +SKIP
    >>> print(report.describe())                             # doctest: +SKIP
    """
    if not isinstance(structure, Structure):
        raise GeometryError(
            f"validate_bonds needs a Structure, got {type(structure).__name__}. Read a file with "
            f"dodo.read_structure first."
        )
    if structure.n_residues == 0 or structure.n_atoms == 0:
        raise GeometryError("Cannot validate the bonds of a structure with no atoms.")

    sigma = _validated_tolerance("tolerance_sigma", tolerance_sigma)
    floor = _validated_tolerance("absolute_tolerance", absolute_tolerance)
    ca_ca = _validated_tolerance("ca_ca_tolerance", ca_ca_tolerance, allow_zero=False)

    context = _prepare(structure, tolerance_sigma=sigma, absolute_tolerance=floor)
    _check_finite(context)

    n_intra = 0
    intra_deviations = np.zeros(0, dtype=np.float64)
    if check_intra_residue:
        n_intra, intra_deviations = _check_intra_residue(context)
    if check_unexpected_contacts:
        _check_unexpected_contacts(context)

    pairs, gaps = _consecutive_pairs(structure)
    n_peptide = 0
    peptide_deviations = np.zeros(0, dtype=np.float64)
    n_ca_ca = 0
    ca_ca_deviations = np.zeros(0, dtype=np.float64)
    n_skipped = 0
    if pairs.size:
        if check_peptide:
            n_peptide, peptide_deviations = _check_peptide_bonds(context, pairs)
        if check_ca_ca:
            n_ca_ca, ca_ca_deviations, n_skipped = _check_ca_ca_bonds(
                context, pairs, tolerance=ca_ca, strict=strict_ca_ca
            )

    n_unknown = _note_unknown_residues(context)
    _note_sequence_gaps(context, gaps)
    n_ca_only = int(np.count_nonzero(context.ca_only))
    if n_ca_only:
        context.notes.append(
            f"{n_ca_only} of {structure.n_residues} residue(s) have an alpha carbon and nothing "
            f"else. That is what a rebuilt region is, so those residues were checked on their "
            f"CA-CA virtual bonds only and no atom was reported missing from them."
        )
    # Partially built residues are a deliberate state, not a defect, but the reader should know
    # how much of the structure was exempted from the side-chain checks. Only residues below the
    # tier their own type requires count -- a GLY with a bare backbone is complete, not partial.
    known = context.type_code >= 0
    required = np.where(known, _TYPE_REQUIRED_TIER[np.where(known, context.type_code, 0)], 0)
    n_partial = int(np.count_nonzero((context.residue_tier < required) & ~context.ca_only))
    if n_partial:
        context.notes.append(
            f"{n_partial} residue(s) have a backbone but no side chain beyond CB, so no "
            f"side-chain atom was reported missing from them. That is what DODO's partial "
            f"all-atom output looks like, and it is also how deposited structures model "
            f"side-chain disorder."
        )

    # Ordered by residue index, matching TraceReport, so a diff between two runs reads cleanly.
    # by_severity() is the worst-first view.
    context.violations.sort(key=lambda v: (v.residue_indices[0], v.kind, v.atoms))
    return BondReport(
        violations=tuple(context.violations),
        n_bonds_checked=n_intra + n_peptide + n_ca_ca,
        n_intra_checked=n_intra,
        n_peptide_checked=n_peptide,
        n_ca_ca_checked=n_ca_ca,
        intra_deviations=intra_deviations,
        peptide_deviations=peptide_deviations,
        ca_ca_deviations=ca_ca_deviations,
        n_residues=structure.n_residues,
        n_ca_only_residues=n_ca_only,
        n_partial_residues=n_partial,
        n_unknown_residues=n_unknown,
        n_sequence_gaps=int(gaps.size),
        n_ca_ca_skipped_peptide=n_skipped,
        notes=tuple(context.notes),
        tolerance_sigma=sigma,
        absolute_tolerance=floor,
        ca_ca_tolerance=ca_ca,
    )
