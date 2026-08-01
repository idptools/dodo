"""Structural validation: bond lengths, steric clashes, and CONECT sanity.

Three checks, one entry point. :func:`validate_structure` runs all of them and returns a combined
report, which is what most callers want; the individual validators are available for when only
one matters.

    from dodo.validate import validate_structure

    report = validate_structure(structure)
    if not report.ok:
        print(report.describe())

Every reference value these use was measured, not transcribed -- see
:mod:`dodo.validate.reference`, derived from 105,299,848 bonds across 23,587 AlphaFold structures.

The invariant this package enforces
-----------------------------------
Every validator needs exclusions, and every exclusion is a place a defect can hide. All three of
these independently acquired the same blind spot: a pair of atoms at exactly 0.000 A was invisible
to each one whenever the pair fell inside an exclusion -- non-consecutive residue numbering,
within-three-bonds, sequence-adjacent, or beneath a flat distance floor.

That is not an academic gap. Coincident atoms are the defect that actually shipped: DODO rebuilt a
region's alpha carbons and left the rest of each residue behind, and the result was visible in a
viewer as long spurious lines through the structure. Nothing in the codebase caught it.

So :func:`validate_structure` always runs :func:`~dodo.validate.impossible.find_impossible_pairs`
and merges its findings unconditionally, before and independently of anything the three validators
choose to exclude. **No exclusion may suppress a physically impossible distance.** Keeping that
check in its own module, with no ``exclude`` parameter to pass, is what stops the invariant being
quietly re-broken.

What counts as a defect
-----------------------
DODO moves folded-domain atoms as rigid bodies and never rebuilds them, so geometry it faithfully
preserved from a flawed input is not DODO's defect. The bond validator attributes each finding to
``rebuilt`` or ``input`` geometry for exactly that reason, and it matters in practice: AlphaFold's
own geometry degrades where AlphaFold is unsure, with 2.48% of its CA-CA bonds falling below 3.3 A
at a mean pLDDT of 38.8. Validating an unmodified AlphaFold model strictly will flag those. They
are not bugs here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..structure import Structure
from .bonds import BondReport, BondViolation, validate_bonds
from .clashes import ClashReport, ClashViolation, validate_clashes
from .conect import ConectReport, ConectViolation, validate_conect_file, validate_conect_lines
from .impossible import IMPOSSIBLE_SEPARATION, ImpossiblePair, find_impossible_pairs

__all__ = [
    "IMPOSSIBLE_SEPARATION",
    "BondReport",
    "BondViolation",
    "ClashReport",
    "ClashViolation",
    "ConectReport",
    "ConectViolation",
    "ImpossiblePair",
    "StructureReport",
    "find_impossible_pairs",
    "validate_bonds",
    "validate_clashes",
    "validate_conect_file",
    "validate_conect_lines",
    "validate_structure",
]


@dataclass(slots=True)
class StructureReport:
    """Combined result of every structural check.

    Attributes
    ----------
    impossible
        Atom pairs closer than any real bond. These are reported unconditionally, ahead of every
        exclusion the other validators apply, and their presence alone makes a structure invalid.
    bonds
        Bond-length findings.
    clashes
        Steric-clash findings.
    conect
        CONECT findings, present only when a written file was validated.
    """

    impossible: tuple[ImpossiblePair, ...] = ()
    bonds: BondReport | None = None
    clashes: ClashReport | None = None
    conect: ConectReport | None = None
    notes: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """True if every check passed."""
        return (
            not self.impossible
            and (self.bonds is None or self.bonds.ok)
            and (self.clashes is None or self.clashes.ok)
            and (self.conect is None or self.conect.ok)
        )

    def __bool__(self) -> bool:
        return self.ok

    @property
    def n_findings(self) -> int:
        """Total findings across every check."""
        return (
            len(self.impossible)
            + (len(self.bonds.violations) if self.bonds else 0)
            + (len(self.clashes.violations) if self.clashes else 0)
            + (len(self.conect.violations) if self.conect else 0)
        )

    def summary(self) -> str:
        """One-line summary of every check."""
        parts = []
        if self.impossible:
            coincident = sum(1 for p in self.impossible if p.coincident)
            detail = f" ({coincident} coincident)" if coincident else ""
            parts.append(f"{len(self.impossible)} IMPOSSIBLE{detail}")
        if self.bonds is not None:
            parts.append(f"{len(self.bonds.violations)} bond")
        if self.clashes is not None:
            parts.append(f"{len(self.clashes.violations)} clash")
        if self.conect is not None:
            parts.append(f"{len(self.conect.violations)} CONECT")
        status = "valid" if self.ok else "INVALID"
        return f"{status}: " + ", ".join(parts) if parts else status

    def describe(self, max_per_check: int = 10) -> str:
        """Multi-line report, worst findings first.

        Impossible separations come first regardless of count, because they are the only class
        that cannot be explained away by an unusual-but-real conformation.
        """
        lines = [self.summary()]
        if self.impossible:
            lines.append("")
            lines.append("Physically impossible separations (no exclusion applies):")
            lines += [f"  {p.message}" for p in self.impossible[:max_per_check]]
            if len(self.impossible) > max_per_check:
                lines.append(f"  ... and {len(self.impossible) - max_per_check} more")
        for label, report in (
            ("Bond lengths", self.bonds),
            ("Steric clashes", self.clashes),
            ("CONECT records", self.conect),
        ):
            if report is None or report.ok:
                continue
            lines.append("")
            lines.append(f"{label}:")
            lines.append("  " + report.describe(max_per_check).replace("\n", "\n  "))
        lines += [f"note: {n}" for n in self.notes]
        return "\n".join(lines)

    def raise_if_invalid(self) -> None:
        """Raise :class:`~dodo.exceptions.GeometryError` if anything failed."""
        if self.ok:
            return
        from ..exceptions import GeometryError

        raise GeometryError(self.describe())


def validate_structure(
    structure: Structure,
    *,
    path: str | Path | None = None,
    check_bonds: bool = True,
    check_clashes: bool = True,
) -> StructureReport:
    """Run every structural check and return a combined report.

    Parameters
    ----------
    structure
        The structure to validate.
    path
        A written PDB file for this structure. When given, its CONECT records are validated too;
        that check needs the file because CONECT is a file-level construct that a
        :class:`~dodo.structure.Structure` does not carry.
    check_bonds, check_clashes
        Turn individual checks off. The impossible-separation check cannot be turned off, by
        design -- see the module docstring.

    Returns
    -------
    StructureReport
        Findings from every check that ran.
    """
    report = StructureReport()

    # Always, first, and with no exclusions. See the module docstring for why this is not
    # optional and why it does not live inside the other validators.
    report.impossible = find_impossible_pairs(structure)

    if check_bonds:
        report.bonds = validate_bonds(structure)
    if check_clashes:
        report.clashes = validate_clashes(structure)
    if path is not None:
        report.conect = validate_conect_file(path)

    return report
