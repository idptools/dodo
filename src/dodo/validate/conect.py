"""CONECT validation: are the bonds a written PDB file declares actually bonds.

Why this module exists
----------------------
Both halves of this check correspond to a defect that shipped.

*Spurious bonds.* DODO rebuilt a region's alpha carbons and left that residue's other
atoms at their original positions, so :mod:`dodo.io.write` emitted a CONECT record
bonding ``N`` to ``CA`` across 93 A. A viewer draws that as one long straight line
through the middle of the structure. Nothing in the suite noticed, because every test
looked at *which serials* were bonded and none looked at how far apart they were.

*Missing bonds.* The CA-CA virtual bond is 3.81 A, past the automatic bond-detection
cutoff in VMD and PyMOL, so a rebuilt CA-only region whose CONECT records are absent
renders as a cloud of disconnected dots -- the entire deliverable, lost at the last step.
That is why :func:`validate_conect_lines` also reports the bonds that *should* be there
and are not, and why the report carries a completeness fraction rather than only a list
of complaints.

What DODO's output actually looks like
--------------------------------------
A DODO file legitimately mixes all-atom stretches (folded domains, every atom of the
input) with CA-only stretches (rebuilt IDRs and loops), and the boundary between them is
a real chain connection. So this module never assumes a residue has a full backbone: a
junction is bonded ``C``-``N`` where both atoms exist and ``CA``-``CA`` where they do
not, exactly as the writer decides it, and a residue with no ``N`` is not a finding.

Does DODO's writer emit reciprocal records?
-------------------------------------------
**No.** :func:`dodo.io.write.structure_to_pdb_lines` records each bond exactly once,
from the lower serial, because every viewer DODO targets accepts single-direction CONECT.
Every bond in a DODO file is therefore non-reciprocal by construction, which is why
``require_reciprocal`` defaults to ``False``: turning it on for DODO output reports one
violation per bond and drowns everything else. The count is always available as
:attr:`ConectReport.n_non_reciprocal`, and the flag exists for files from other writers
(wwPDB depositions are reciprocal -- ``tests/data/structures/6kn7.pdb`` is, for its
ligands) where a one-directional record may indicate a truncated file.

Multi-model files
-----------------
DODO writes CONECT once after the last ``ENDMDL``, so one set of serials has to describe
every frame. Distances are therefore measured in *every* frame and the worst one is
reported, and each serial's atom identity is compared across frames. That second check is
not hypothetical: TER records consume an atom serial, so two frames that split their
chains differently are numbered differently, and frame 1's CONECT records then bond the
wrong atoms in frame 2 -- including bonding to a serial that is a TER record while a real
atom is left with no bond at all. :func:`dodo.io.write.write_pdb` guards against writing
such a file; this module catches one that already exists.

Everything here is a report, not an exception. Only unusable input -- an unreadable path,
a file with no atom records, a garbled coordinate field -- raises.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field
from itertools import pairwise
from pathlib import Path
from typing import Final, Literal

import numpy as np

from ..constants import (
    BACKBONE_ATOMS,
    C_N_PEPTIDE_BOND_LENGTH,
    CA_CA_BOND_LENGTH,
    THREE_TO_ONE,
)
from ..exceptions import EmptyStructureError, MalformedRecordError, StructureFileError
from ..io.pdb import decode_hybrid36

__all__ = [
    "MAX_BOND_LENGTH",
    "MIN_BOND_LENGTH",
    "ConectReport",
    "ConectViolation",
    "ConectViolationKind",
    "validate_conect_file",
    "validate_conect_lines",
]


# ---------------------------------------------------------------------------
# Thresholds
#
# Derived from dodo.constants, never typed in. The writer's connectivity decisions come
# from the same two constants, so a change to CA_CA_BOND_LENGTH moves both together.
# ---------------------------------------------------------------------------

#: Slack on every distance comparison, in Angstroms.
#:
#: PDB coordinates are written ``%8.3f``, so each of the six coordinates behind a
#: distance is rounded by up to 0.0005 A and the distance itself can differ from the
#: in-memory value by up to ``sqrt(3) * 0.001`` = 0.0017 A. This is that bound rounded up.
#: Without it, a bond written at exactly the threshold is a coin flip.
COORDINATE_PRECISION_MARGIN: Final[float] = 0.005

#: Multiplier on :data:`~dodo.constants.CA_CA_BOND_LENGTH` giving the longest CA-CA
#: separation that still counts as bonded.
#:
#: This is the same 30% slack ``dodo.io.write._MAX_BONDED_CA_CA_DISTANCE`` uses to decide
#: whether to emit a record, and the two are asserted equal in
#: ``tests/unit/test_validate_conect.py`` so they cannot drift apart. Matching the writer
#: is the point: a bond the writer would have refused to draw, present in a file the
#: writer produced, is a defect, and anything looser lets a real one through. The slack
#: is not decoration -- measured on ``ca_only=True`` output of dnmt3a and arf19, real
#: folded-domain CA-CA bonds run 2.978 to 3.914 A (cis-proline at the low end, strained
#: crystal geometry at the high end), so a threshold at
#: ``CA_CA_BOND_LENGTH + CA_CA_BOND_TOLERANCE`` = 3.91 A would flag deposited geometry
#: DODO merely passed through.
MAX_BONDED_CA_CA_FACTOR: Final[float] = 1.3

#: Longest distance a CONECT record may span, in Angstroms. DERIVED.
MAX_BOND_LENGTH: Final[float] = (
    MAX_BONDED_CA_CA_FACTOR * CA_CA_BOND_LENGTH + COORDINATE_PRECISION_MARGIN
)

#: Longest C(i)-N(i+1) separation that counts as a peptide bond, in Angstroms. DERIVED,
#: with the same 50% slack the writer applies, and used for the same purpose: deciding
#: whether a junction with no alpha carbon is a bond or an unmodelled-residue gap.
MAX_PEPTIDE_BOND_LENGTH: Final[float] = 1.5 * C_N_PEPTIDE_BOND_LENGTH

#: Shortest distance a CONECT record may span, in Angstroms.
#:
#: Below this the two atoms are coincident rather than bonded. The shortest bond in
#: :mod:`dodo.validate.reference` is a carbonyl C=O at 1.222 A, and the shortest bond in
#: any structure is X-H at about 0.96 A, so 0.9 A sits under everything real while still
#: catching the two-atoms-at-one-point case that a failed rebuild produces. It is not
#: raised to 1.1 A precisely so that a file carrying hydrogens -- which DODO passes
#: through from its input untouched -- is not reported as broken.
MIN_BOND_LENGTH: Final[float] = 0.9

#: Width of the atom serial field, in columns. Fixed by the PDB format.
_SERIAL_WIDTH: Final[int] = 5

#: 1-based column of the first serial field on a CONECT record (the origin atom).
_CONECT_FIRST_FIELD: Final[int] = 7

#: Bonded-partner fields the specification defines: columns 12-16, 17-21, 22-26 and
#: 27-31. DODO's writer emits at most this many per record and continues on the next.
_CONECT_DOCUMENTED_PARTNERS: Final[int] = 4


ConectViolationKind = Literal[
    "bond_too_long",
    "bond_too_short",
    "unknown_serial",
    "self_bond",
    "duplicate_bond",
    "non_reciprocal",
    "interchain_bond",
    "missing_bond",
    "ambiguous_serial",
    "frame_mismatch",
    "malformed_record",
]


# ---------------------------------------------------------------------------
# Fixed-column parsing
#
# By column, never by split(). A PDB line is a fixed-width record and split() is wrong on
# real files in both directions: a blank altLoc or insertion code makes fields vanish, and
# a negative coordinate abutting its neighbour ("-12.345-6.789") makes two fields merge.
# The 93 A bug this module exists to catch lived in a file that split() reads happily.
# ---------------------------------------------------------------------------


def _columns(line: str, first: int, last: int) -> str:
    """Return 1-based inclusive columns ``first`` to ``last`` of ``line``.

    Short lines are legal -- most writers strip trailing blanks -- so a slice past the
    end yields the empty string rather than raising.
    """
    return line[first - 1 : last]


@dataclass(frozen=True, slots=True)
class _Atom:
    """One ATOM/HETATM record, reduced to what a connectivity check needs."""

    serial: int
    name: str
    residue_label: str
    chain_id: str
    residue_key: tuple[str, str, str, str]
    residue_name: str
    line_number: int
    #: True when a TER record or a model boundary separates this atom from the previous
    #: one, i.e. this atom starts a new polymer chain whatever the chain-id column says.
    break_before: bool

    @property
    def identity(self) -> str:
        """``"A:GLU142 CA"`` -- the residue label plus the atom name."""
        return f"{self.residue_label} {self.name}"


@dataclass(frozen=True, slots=True)
class _Residue:
    """A run of atoms sharing a residue key, in file order."""

    label: str
    chain_id: str
    residue_name: str
    is_polymer: bool
    break_before: bool
    #: Backbone atom name to position within the frame's atom list. Where a residue
    #: carries the same backbone name twice (unresolved alternate conformers, which
    #: DODO's reader preserves), the last one wins -- the same choice
    #: ``dodo.io.write._backbone_positions`` makes, so the expected-bond set names the
    #: atom the writer actually bonded instead of reporting its sibling as unbonded.
    backbone: dict[str, int]


@dataclass(slots=True)
class _Frame:
    """One model's worth of atoms."""

    model_num: int | None
    atoms: list[_Atom] = field(default_factory=list)
    coords: list[tuple[float, float, float]] = field(default_factory=list)
    by_serial: dict[int, int] = field(default_factory=dict)
    #: Serials used by more than one atom in this frame, with their use count. A CONECT
    #: record naming one of these is ambiguous, so it is reported rather than resolved.
    repeated_serials: dict[int, int] = field(default_factory=dict)
    #: Serials consumed by TER records. Referencing one is not a bond to an atom.
    ter_serials: dict[int, int] = field(default_factory=dict)

    @property
    def label(self) -> str:
        """``"model 2"``, or ``"the file"`` when there are no MODEL records."""
        return "the file" if self.model_num is None else f"model {self.model_num}"

    def xyz(self) -> np.ndarray:
        """``(n, 3)`` coordinates in file order."""
        if not self.coords:
            return np.zeros((0, 3), dtype=np.float64)
        return np.asarray(self.coords, dtype=np.float64)


@dataclass(frozen=True, slots=True)
class _Record:
    """One parsed CONECT record: an origin serial and the serials it bonds to."""

    origin: int
    partners: tuple[int, ...]
    line_number: int


def _decode_serial(text: str, line_number: int, line: str, what: str) -> int:
    """Decode a serial field, raising :class:`MalformedRecordError` if it cannot be."""
    try:
        return decode_hybrid36(text, _SERIAL_WIDTH)
    except ValueError as error:
        raise MalformedRecordError(
            f"{what} {text!r} is not a valid atom serial ({error})", line_number, line.rstrip()
        ) from error


def _residue_label(line: str) -> tuple[str, str, str, tuple[str, str, str, str]]:
    """Extract ``(label, chain_id, residue_name, residue_key)`` from an atom record.

    The label format matches :meth:`dodo.structure.Structure.residue_label` exactly
    (``"A:GLU142"``), so a finding from this module and a finding from the geometry
    checks name the same residue the same way.

    The residue *name* is part of the key as well as the number, so the microheterogeneity
    real depositions contain -- two different residues sharing one sequence number -- is
    two residues here rather than one residue with two conflicting backbones.
    """
    residue_name = _columns(line, 18, 20).strip().upper()
    chain_id = _columns(line, 22, 22).strip()
    # Columns 23-26 are the sequence number and 27 the insertion code. The number is kept
    # as written rather than decoded: it is only ever used for a label and for grouping,
    # and a hybrid-36 or otherwise exotic field must not fail a connectivity check.
    resseq = _columns(line, 23, 26).strip()
    icode = _columns(line, 27, 27).strip()
    label = f"{chain_id or '?'}:{residue_name}{resseq}{icode}"
    return label, chain_id, residue_name, (chain_id, resseq, icode, residue_name)


@dataclass(slots=True)
class _ParsedFile:
    """Everything :func:`validate_conect_lines` needs from the text of a PDB file."""

    frames: list[_Frame] = field(default_factory=list)
    records: list[_Record] = field(default_factory=list)
    malformed: list[tuple[int, str, str]] = field(default_factory=list)
    n_conect_lines: int = 0
    #: Partner serials read from columns beyond the four the specification defines.
    n_extended_partners: int = 0
    #: CONECT records that appear before the last ENDMDL rather than after it.
    n_records_inside_models: int = 0


def _parse(lines: Iterable[str]) -> _ParsedFile:
    """Parse ATOM/HETATM/TER/MODEL/CONECT records by column."""
    parsed = _ParsedFile()
    current: _Frame | None = None
    break_pending = False
    inside_model = False

    for line_number, raw in enumerate(lines, start=1):
        line = raw.rstrip("\r\n")
        # Trailing blanks are stripped from the record name because real files do not
        # always pad it: a bare "TER" is common, and matching only the padded "TER   "
        # would miss the chain break it marks -- which would then make the residues on
        # either side look like bonded neighbours.
        record = _columns(line, 1, 6).rstrip()

        if record in ("ATOM", "HETATM"):
            if current is None:
                current = _Frame(model_num=None)
                parsed.frames.append(current)
            serial = _decode_serial(
                _columns(line, 7, 11), line_number, line, "the atom serial field"
            )
            try:
                x = float(_columns(line, 31, 38))
                y = float(_columns(line, 39, 46))
                z = float(_columns(line, 47, 54))
            except ValueError as error:
                raise MalformedRecordError(
                    f"the coordinate columns are not three numbers ({error})",
                    line_number,
                    line.rstrip(),
                ) from error
            label, chain_id, residue_name, key = _residue_label(line)
            position = len(current.atoms)
            if serial in current.by_serial:
                current.repeated_serials[serial] = current.repeated_serials.get(serial, 1) + 1
            else:
                current.by_serial[serial] = position
            current.atoms.append(
                _Atom(
                    serial=serial,
                    name=_columns(line, 13, 16).strip(),
                    residue_label=label,
                    chain_id=chain_id,
                    residue_key=key,
                    residue_name=residue_name,
                    line_number=line_number,
                    break_before=break_pending,
                )
            )
            current.coords.append((x, y, z))
            break_pending = False

        elif record == "TER":
            if current is not None:
                field_text = _columns(line, 7, 11).strip()
                if field_text:
                    current.ter_serials[
                        _decode_serial(
                            _columns(line, 7, 11), line_number, line, "the TER serial field"
                        )
                    ] = line_number
            break_pending = True

        elif record == "MODEL":
            model_field = _columns(line, 11, 14).strip()
            try:
                model_num = int(model_field)
            except ValueError:
                # A MODEL record whose number is unreadable still starts a frame; the
                # frames are ordered, so a placeholder number costs nothing and refusing
                # the file over a cosmetic field would be worse.
                model_num = len(parsed.frames) + 1
            current = _Frame(model_num=model_num)
            parsed.frames.append(current)
            break_pending = True
            inside_model = True

        elif record == "ENDMDL":
            break_pending = True
            inside_model = False

        elif record == "CONECT":
            parsed.n_conect_lines += 1
            if inside_model:
                parsed.n_records_inside_models += 1
            _parse_conect(line, line_number, parsed)

    return parsed


def _parse_conect(line: str, line_number: int, parsed: _ParsedFile) -> None:
    """Append one CONECT record to ``parsed``, or record it as malformed.

    Fields are 5 columns wide from column 7: the origin, then the four bonded partners
    the specification defines (columns 12-31), then anything further. Extra fields are
    read as additional partners, which is what a writer that overflows a record means by
    them, and counted so that the one case where that reading is wrong stays visible: PDB
    format 2.2 and earlier put hydrogen-bonded and salt-bridged partners in columns 32-56,
    and those are not covalent bonds. Nothing has emitted them since 1996 and DODO's
    writer emits continuation records instead, so the count is a footnote, not a check.
    """
    # rstrip first: a record padded to 80 columns would otherwise yield blank fields.
    text = line.rstrip()
    fields: list[str] = []
    start = _CONECT_FIRST_FIELD
    while start <= len(text):
        fields.append(_columns(text, start, start + _SERIAL_WIDTH - 1))
        start += _SERIAL_WIDTH
    if not fields or not fields[0].strip():
        parsed.malformed.append((line_number, text, "it names no origin atom"))
        return
    try:
        origin = decode_hybrid36(fields[0], _SERIAL_WIDTH)
    except ValueError as error:
        parsed.malformed.append(
            (line_number, text, f"its origin serial {fields[0]!r} cannot be decoded ({error})")
        )
        return
    partners: list[int] = []
    undecodable = 0
    for index, raw_field in enumerate(fields[1:]):
        if not raw_field.strip():
            continue
        try:
            partners.append(decode_hybrid36(raw_field, _SERIAL_WIDTH))
        except ValueError as error:
            undecodable += 1
            parsed.malformed.append(
                (
                    line_number,
                    text,
                    f"its partner serial {raw_field!r} cannot be decoded ({error})",
                )
            )
            continue
        if index >= _CONECT_DOCUMENTED_PARTNERS:
            parsed.n_extended_partners += 1
    if not partners:
        # Only when there was nothing to decode in the first place. A record whose single
        # partner field is garbage is one defect, and saying so twice reads as two.
        if not undecodable:
            parsed.malformed.append((line_number, text, "it names an origin atom but no partners"))
        return
    parsed.records.append(_Record(origin=origin, partners=tuple(partners), line_number=line_number))


def _residues(frame: _Frame) -> list[_Residue]:
    """Group a frame's atoms into residues, in file order."""
    residues: list[_Residue] = []
    previous_key: tuple[str, str, str, str] | None = None
    for position, atom in enumerate(frame.atoms):
        if previous_key != atom.residue_key or atom.break_before:
            residues.append(
                _Residue(
                    label=atom.residue_label,
                    chain_id=atom.chain_id,
                    residue_name=atom.residue_name,
                    # THREE_TO_ONE is the same test dodo.io.write uses to choose ATOM
                    # over HETATM, so "polymer" means the same thing in both modules.
                    is_polymer=atom.residue_name in THREE_TO_ONE,
                    break_before=atom.break_before,
                    backbone={},
                )
            )
            previous_key = atom.residue_key
        if atom.name in BACKBONE_ATOMS:
            residues[-1].backbone[atom.name] = position
    return residues


# ---------------------------------------------------------------------------
# Findings
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ConectViolation:
    """One specific defect in a file's CONECT records.

    Attributes
    ----------
    kind
        Which check failed. See :data:`ConectViolationKind`.
    serials
        The atom serial numbers involved, as written in the file. One for
        ``ambiguous_serial``, two for every bond-level finding, empty for
        ``malformed_record``.
    atoms
        The atom identities behind ``serials``, as ``"A:GLU142 CA"``. Shorter than
        ``serials`` when a serial does not resolve to an atom, which is the point of the
        ``unknown_serial`` finding.
    distance
        The measured separation in Angstroms, or ``None`` where the finding is not about
        a distance.
    low, high
        The acceptable range ``distance`` fell outside of, or ``None``.
    line_number
        1-based line in the file the finding is anchored to: the CONECT record for a
        declared bond, the first atom's record for a missing one.
    model
        Which frame the measurement comes from, e.g. ``"model 2"``. ``"the file"`` when
        the file has no MODEL records.
    message
        A complete, self-contained sentence naming the residues, the measured value and
        the expected range.
    """

    kind: ConectViolationKind
    serials: tuple[int, ...]
    atoms: tuple[str, ...]
    distance: float | None
    low: float | None
    high: float | None
    line_number: int
    model: str
    message: str

    def __str__(self) -> str:
        return self.message


@dataclass(frozen=True, slots=True, eq=False)
class ConectReport:
    """The result of :func:`validate_conect_lines`: counts, completeness and violations.

    ``eq=False`` for the same reason :class:`dodo.geometry.metrics.TraceReport` uses it:
    a generated ``__eq__`` over the numeric fields is meaningless here and two reports
    are never usefully compared.

    Attributes
    ----------
    violations
        Every finding, ordered by serial then kind.
    n_records
        CONECT *lines* in the file, including continuation records.
    n_bonds
        Distinct undirected bonds those records declare. Lower than ``n_records`` when
        records carry several partners, and lower than the number of directed pairs when
        the file is reciprocal.
    n_expected_bonds
        Backbone bonds the file's own coordinates imply: intra-residue N-CA, CA-C and C-O
        wherever both atoms are present, plus one bond per bonded residue junction.
    n_expected_present
        How many of those are actually declared. See :attr:`completeness`.
    n_non_reciprocal
        Declared bonds recorded from one atom but not from the other. Expected to equal
        ``n_bonds`` for a DODO file -- see the module docstring -- and reported as
        violations only when ``require_reciprocal`` is set.
    n_atoms, n_models
        Atoms per frame and number of frames.
    notes
        Anomalies that are not violations, and checks that were skipped and why.
    source
        The file the report describes, for messages.
    """

    violations: tuple[ConectViolation, ...]
    n_records: int
    n_bonds: int
    n_expected_bonds: int
    n_expected_present: int
    n_non_reciprocal: int
    n_atoms: int
    n_models: int
    notes: tuple[str, ...]
    source: str | None

    @property
    def ok(self) -> bool:
        """True when nothing was violated."""
        return not self.violations

    def __bool__(self) -> bool:
        return self.ok

    def __repr__(self) -> str:
        return (
            f"ConectReport({self.n_records} records, {self.n_bonds} bonds, "
            f"{'valid' if self.ok else f'{len(self.violations)} violations'})"
        )

    @property
    def completeness(self) -> float:
        """Fraction of the bonds the coordinates imply that the file actually declares.

        1.0 when nothing was expected, so an empty or ligand-only file is not reported as
        0% complete. Anything below 1.0 on a DODO file means some of the structure will
        render as disconnected dots.
        """
        if self.n_expected_bonds == 0:
            return 1.0
        return self.n_expected_present / self.n_expected_bonds

    def of_kind(self, kind: ConectViolationKind) -> tuple[ConectViolation, ...]:
        """All violations of one kind, in report order."""
        return tuple(v for v in self.violations if v.kind == kind)

    @property
    def kinds(self) -> tuple[ConectViolationKind, ...]:
        """The distinct kinds present, in report order."""
        seen: dict[ConectViolationKind, None] = {}
        for violation in self.violations:
            seen[violation.kind] = None
        return tuple(seen)

    def summary(self) -> str:
        """One-line summary of what was checked, valid or not."""
        where = f"{self.source}: " if self.source else ""
        frames = "1 model" if self.n_models == 1 else f"{self.n_models} models"
        parts = [
            f"{where}{self.n_atoms} atoms in {frames}",
            f"{self.n_records} CONECT record(s) declaring {self.n_bonds} bond(s)",
            f"{self.n_expected_present}/{self.n_expected_bonds} expected bond(s) present "
            f"({100.0 * self.completeness:.1f}%)",
        ]
        if self.ok:
            parts.append("valid")
        else:
            counts = Counter(v.kind for v in self.violations)
            detail = ", ".join(f"{count} {kind}" for kind, count in sorted(counts.items()))
            parts.append(f"{len(self.violations)} violation(s): {detail}")
        return "; ".join(parts)

    def describe(self, max_violations: int = 10) -> str:
        """Full multi-line description: the summary, the notes, then the violations.

        Parameters
        ----------
        max_violations
            Truncate the violation list after this many, with a count of the remainder.
            A file with no CONECT records at all produces one violation per bond, and an
            unreadable message is no better than no message.
        """
        lines = [self.summary()]
        lines.extend(f"  note: {note}" for note in self.notes)
        for violation in self.violations[:max_violations]:
            lines.append(f"  - {violation.message}")
        remaining = len(self.violations) - max_violations
        if remaining > 0:
            lines.append(f"  - ... and {remaining} more violation(s)")
        return "\n".join(lines)

    def raise_if_invalid(self) -> None:
        """Raise :class:`~dodo.exceptions.StructureFileError` if anything was violated.

        For a caller that has just written a file and must not hand it to a user if its
        connectivity is wrong.
        """
        if not self.ok:
            raise StructureFileError(f"Invalid CONECT records: {self.describe()}")


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------


def _identity_of(frame: _Frame, serial: int) -> str:
    """Describe the atom a serial names in one frame, for a message."""
    position = frame.by_serial.get(serial)
    if position is None:
        if serial in frame.ter_serials:
            return f"a TER record in {frame.label}"
        return f"nothing in {frame.label}"
    return frame.atoms[position].identity


def _check_serials(
    frames: list[_Frame],
    directed: dict[tuple[int, int], list[int]],
    violations: list[ConectViolation],
) -> set[int]:
    """Report serials that no frame resolves, or that frames disagree about.

    Returns the serials that resolve to the same atom in every frame, which are the only
    ones whose geometry can be measured.
    """
    # One pass to find where each serial is first mentioned. Scanning `directed` per
    # serial instead would be quadratic, and this runs over ~62,000 serials on 6kn7.
    first_line: dict[int, int] = {}
    for (origin, partner), line_numbers in directed.items():
        for serial in (origin, partner):
            earliest = first_line.get(serial)
            if earliest is None or line_numbers[0] < earliest:
                first_line[serial] = line_numbers[0]
    usable: set[int] = set()
    for serial in sorted(first_line):
        line_number = first_line[serial]
        present = [frame for frame in frames if serial in frame.by_serial]
        if not present:
            ter_frames = [frame.label for frame in frames if serial in frame.ter_serials]
            explanation = (
                f"it is a TER record in {', '.join(ter_frames)}, which consumes an atom "
                f"serial but is not an atom"
                if ter_frames
                else "no ATOM or HETATM record in the file has that serial"
            )
            violations.append(
                ConectViolation(
                    kind="unknown_serial",
                    serials=(serial,),
                    atoms=(),
                    distance=None,
                    low=None,
                    high=None,
                    line_number=line_number,
                    model=frames[0].label,
                    message=(
                        f"The CONECT record on line {line_number} refers to atom serial "
                        f"{serial}, which does not exist: {explanation}. A viewer either "
                        f"drops the bond or bonds the wrong atom."
                    ),
                )
            )
            continue
        repeated = [frame for frame in frames if serial in frame.repeated_serials]
        if repeated:
            frame = repeated[0]
            violations.append(
                ConectViolation(
                    kind="ambiguous_serial",
                    serials=(serial,),
                    atoms=(_identity_of(frame, serial),),
                    distance=None,
                    low=None,
                    high=None,
                    line_number=line_number,
                    model=frame.label,
                    message=(
                        f"Atom serial {serial} is used by "
                        f"{frame.repeated_serials[serial]} different atoms in "
                        f"{frame.label}, so the CONECT record on line {line_number} that "
                        f"names it does not identify one bond. Serials must be unique "
                        f"within a model."
                    ),
                )
            )
            continue
        if len(present) != len(frames):
            missing = [frame.label for frame in frames if serial not in frame.by_serial]
            violations.append(
                ConectViolation(
                    kind="frame_mismatch",
                    serials=(serial,),
                    atoms=(_identity_of(present[0], serial),),
                    distance=None,
                    low=None,
                    high=None,
                    line_number=line_number,
                    model=missing[0],
                    message=(
                        f"Atom serial {serial} is {_identity_of(present[0], serial)} in "
                        f"{present[0].label} but is not an atom in {', '.join(missing)}, "
                        f"and the file's single set of CONECT records has to describe "
                        f"every frame. The usual cause is frames whose TER records fall "
                        f"in different places, since a TER record consumes a serial."
                    ),
                )
            )
            continue
        identities = {_identity_of(frame, serial) for frame in frames}
        if len(identities) > 1:
            violations.append(
                ConectViolation(
                    kind="frame_mismatch",
                    serials=(serial,),
                    atoms=tuple(sorted(identities)),
                    distance=None,
                    low=None,
                    high=None,
                    line_number=line_number,
                    model=frames[1].label,
                    message=(
                        f"Atom serial {serial} names different atoms in different frames "
                        f"({', '.join(sorted(identities))}), so the CONECT record on line "
                        f"{line_number} bonds the right atoms in one frame and the wrong "
                        f"atoms in another. The frames must share one topology and one "
                        f"TER placement."
                    ),
                )
            )
            continue
        usable.add(serial)
    return usable


def _measure(
    frames: list[_Frame], bonds: list[tuple[int, int]]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Measure every declared bond in every frame.

    Returns the longest and shortest distance found for each bond, and the index of the
    frame each extreme came from. Every bond is measured in every frame because a
    multi-model file carries one set of CONECT records for all of them, so a bond that is
    fine in frame 1 and 93 A long in frame 3 has to be caught.

    Vectorized per frame rather than per bond. There is no neighbour search here -- the
    pairs are given, not discovered, so a KD-tree would have nothing to do -- but the
    61,511-atom 6kn7 assembly written with CONECT declares ~62,000 bonds, and a Python
    loop over those costs several times what the whole parse does.
    """
    n = len(bonds)
    longest = np.full(n, -np.inf, dtype=np.float64)
    shortest = np.full(n, np.inf, dtype=np.float64)
    longest_frame = np.zeros(n, dtype=np.int64)
    shortest_frame = np.zeros(n, dtype=np.int64)
    if n == 0:
        return longest, shortest, longest_frame, shortest_frame
    for index, frame in enumerate(frames):
        lookup = frame.by_serial
        left = np.fromiter((lookup[a] for a, _ in bonds), dtype=np.int64, count=n)
        right = np.fromiter((lookup[b] for _, b in bonds), dtype=np.int64, count=n)
        xyz = frame.xyz()
        distances: np.ndarray = np.linalg.norm(xyz[left] - xyz[right], axis=1)
        further = distances > longest
        longest[further] = distances[further]
        longest_frame[further] = index
        nearer = distances < shortest
        shortest[nearer] = distances[nearer]
        shortest_frame[nearer] = index
    return longest, shortest, longest_frame, shortest_frame


def _expected_bonds(frame: _Frame, max_bond_length: float) -> list[tuple[int, int, str, float]]:
    """Return the backbone bonds a frame's coordinates imply: ``(serial, serial, why, A)``.

    Mirrors :func:`dodo.io.write._bonds` deliberately, decision for decision, because
    the question being asked is "did the writer emit what its own rules say it should".
    Intra-residue N-CA, CA-C and C-O wherever both atoms exist -- so a CA-only rebuilt
    region contributes none, which is correct and is why this module does not report
    "residue has no N" thousands of times. Then one bond per junction: C-N when both
    peptide atoms exist, CA-CA when they do not.

    A junction is only expected where the two residues are close enough to be bonded.
    That is what keeps a genuine chain break -- unmodelled residues, of which
    ``tests/data/structures/6kn7.pdb`` has several -- from being reported as a missing
    bond, and it is decided on distance rather than on residue numbering for the reason
    the writer decides it that way: numbering is bookkeeping and can be renumbered or
    carry insertion codes, while distance is the truth.

    Only polymer residues (:data:`~dodo.constants.THREE_TO_ONE`) take part. A ligand or a
    nucleotide whose atoms happen to be named ``C``, ``N`` or ``O`` -- ADP in 6kn7 has all
    three -- would otherwise be reported as missing backbone bonds it never had.
    """
    residues = _residues(frame)
    atoms = frame.atoms
    xyz = frame.xyz()
    expected: list[tuple[int, int, str, float]] = []

    def serial(position: int) -> int:
        return atoms[position].serial

    def distance(left: int, right: int) -> float:
        return float(np.linalg.norm(xyz[left] - xyz[right]))

    for residue in residues:
        if not residue.is_polymer:
            continue
        for first, second in (("N", "CA"), ("CA", "C"), ("C", "O")):
            left = residue.backbone.get(first)
            right = residue.backbone.get(second)
            if left is None or right is None:
                continue
            expected.append(
                (serial(left), serial(right), f"the {first}-{second} bond", distance(left, right))
            )

    for residue, following in pairwise(residues):
        if not (residue.is_polymer and following.is_polymer):
            continue
        # A TER record ends a chain even when the next residue reuses the chain id, which
        # deposited files do routinely: 6kn7 parks its ADP ligands in chain A after TER.
        if following.break_before or residue.chain_id != following.chain_id:
            continue
        # Held as pairs rather than as four Optionals so that the two "both atoms exist"
        # tests narrow once, here, instead of being re-asserted at every use.
        alpha_carbons = _pair(residue.backbone.get("CA"), following.backbone.get("CA"))
        peptide = _pair(residue.backbone.get("C"), following.backbone.get("N"))
        if alpha_carbons is not None:
            separation = distance(*alpha_carbons)
            bonded = separation <= max_bond_length
        elif peptide is not None:
            separation = distance(*peptide)
            bonded = separation <= MAX_PEPTIDE_BOND_LENGTH
        else:
            continue
        if not bonded:
            continue
        if peptide is not None:
            expected.append(
                (
                    serial(peptide[0]),
                    serial(peptide[1]),
                    f"the peptide bond to {following.label}",
                    distance(*peptide),
                )
            )
        else:
            assert alpha_carbons is not None  # only reachable via the CA branch above
            expected.append(
                (
                    serial(alpha_carbons[0]),
                    serial(alpha_carbons[1]),
                    f"the CA-CA virtual bond to {following.label}",
                    separation,
                )
            )
    return expected


def _pair(left: int | None, right: int | None) -> tuple[int, int] | None:
    """``(left, right)`` when both are present, else ``None``."""
    if left is None or right is None:
        return None
    return left, right


def _declares_polymer_bonds(frame: _Frame, bonds: Iterable[tuple[int, int]]) -> bool:
    """Whether the file declares any polymer backbone connectivity at all.

    Two legitimate kinds of file declare none, and neither should be buried under one
    missing-bond violation per residue:

    * ``write_pdb(..., conect=False)`` output, which has no CONECT records whatsoever.
    * A wwPDB deposition. Depositions never list backbone bonds -- readers infer those
      from the residue names -- and only record CONECT for ligands, disulfides and other
      links. ``tests/data/structures/6kn7.pdb`` has 405 such records and not one
      backbone bond among them.

    In both cases the completeness fraction is still reported; only the violations are
    suppressed, and a note says so.
    """
    polymer_backbone: set[int] = {
        atom.serial
        for atom in frame.atoms
        if atom.name in BACKBONE_ATOMS and atom.residue_name in THREE_TO_ONE
    }
    return any(left in polymer_backbone and right in polymer_backbone for left, right in bonds)


def validate_conect_lines(
    lines: Iterable[str],
    *,
    source: str | None = None,
    max_bond_length: float = MAX_BOND_LENGTH,
    min_bond_length: float = MIN_BOND_LENGTH,
    require_reciprocal: bool = False,
    require_backbone_conect: bool | None = None,
) -> ConectReport:
    """Validate the CONECT records in PDB text.

    Parameters
    ----------
    lines
        The file's lines, with or without line endings. Consumed once, so a file handle
        is fine.
    source
        Name to use in messages, typically the path.
    max_bond_length
        Longest distance a CONECT record may span, in Angstroms. Defaults to
        :data:`MAX_BOND_LENGTH`, derived from
        :data:`~dodo.constants.CA_CA_BOND_LENGTH`.
    min_bond_length
        Shortest distance a CONECT record may span. Defaults to :data:`MIN_BOND_LENGTH`.
    require_reciprocal
        Report a violation for each bond recorded from only one of its two atoms. Off by
        default because DODO's writer records every bond once, from the lower serial --
        see the module docstring. :attr:`ConectReport.n_non_reciprocal` is counted either
        way.
    require_backbone_conect
        Whether missing backbone bonds are violations. ``None`` (the default) decides per
        file: a file that declares no polymer backbone connectivity at all is taken to be
        ``conect=False`` output or a wwPDB deposition, and its missing bonds are counted
        but not reported. ``True`` reports them regardless, which is what a test asserting
        that ``conect=True`` output is complete wants. ``False`` never reports them.

    Returns
    -------
    ConectReport
        Counts, completeness and every violation. Truthy when the file is valid.

    Raises
    ------
    EmptyStructureError
        If there is not a single ATOM or HETATM record to check against. Nothing can be
        said about connectivity without coordinates.
    MalformedRecordError
        If an atom record's serial or coordinate columns cannot be parsed. A garbled
        CONECT record is *reported*, because that is the thing under test; a garbled
        coordinate makes the file unusable as a reference.
    ValueError
        If the thresholds are not finite and positive with ``min < max``.

    Examples
    --------
    >>> lines = [
    ...     "ATOM      1  CA  MET A   1       0.000   0.000   0.000  1.00  0.00           C",
    ...     "ATOM      2  CA  PRO A   2       3.810   0.000   0.000  1.00  0.00           C",
    ...     "CONECT    1    2",
    ...     "END",
    ... ]
    >>> report = validate_conect_lines(lines)
    >>> report.ok, report.n_bonds, report.completeness
    (True, 1, 1.0)

    Move the second alpha carbon 93 A away and the bond becomes the defect this module
    was written for:

    >>> lines[1] = lines[1].replace("   3.810", "  93.810")
    >>> [v.kind for v in validate_conect_lines(lines).violations]
    ['bond_too_long']

    No ``missing_bond`` accompanies it: at 93 A those two residues are no longer close
    enough to count as bonded, which is the same test that keeps a genuine chain break
    from being reported as a missing record.
    """
    for name, value in (("max_bond_length", max_bond_length), ("min_bond_length", min_bond_length)):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive, got {value}.")
    if min_bond_length >= max_bond_length:
        raise ValueError(
            f"min_bond_length ({min_bond_length}) must be below max_bond_length "
            f"({max_bond_length})."
        )

    parsed = _parse(lines)
    frames = [frame for frame in parsed.frames if frame.atoms]
    if not frames:
        raise EmptyStructureError(
            f"{source or 'This PDB text'} contains no ATOM or HETATM records, so there "
            f"is nothing for its CONECT records to refer to."
        )

    violations: list[ConectViolation] = []
    notes: list[str] = []
    reference = frames[0]

    for line_number, text, why in parsed.malformed:
        violations.append(
            ConectViolation(
                kind="malformed_record",
                serials=(),
                atoms=(),
                distance=None,
                low=None,
                high=None,
                line_number=line_number,
                model=reference.label,
                message=(
                    f"The CONECT record on line {line_number} cannot be read because "
                    f"{why}: {text!r}. CONECT fields are five columns wide from column 7."
                ),
            )
        )

    # Directed pairs, with the line numbers each was recorded on. Duplicate detection
    # needs the count per direction, and reciprocity needs the direction itself.
    directed: dict[tuple[int, int], list[int]] = {}
    for record in parsed.records:
        for partner in record.partners:
            directed.setdefault((record.origin, partner), []).append(record.line_number)

    for (origin, partner), line_numbers in sorted(directed.items()):
        if origin != partner:
            continue
        identity = _identity_of(reference, origin)
        violations.append(
            ConectViolation(
                kind="self_bond",
                serials=(origin, partner),
                atoms=(identity,),
                # No distance: an atom's separation from itself is not a measurement, and
                # reporting 0.000 A would read as a coincident-atom finding.
                distance=None,
                low=None,
                high=None,
                line_number=line_numbers[0],
                model=reference.label,
                message=(
                    f"The CONECT record on line {line_numbers[0]} bonds atom serial "
                    f"{origin} ({identity}) to itself. No atom is bonded to itself; a "
                    f"viewer draws a zero-length bond or rejects the record."
                ),
            )
        )

    for (origin, partner), line_numbers in sorted(directed.items()):
        if len(line_numbers) < 2:
            continue
        # A record may also list the same partner twice in its own fields, in which case
        # every occurrence carries the same line number and "on lines 5, 5" reads as a
        # transcription error in the message rather than in the file.
        unique_lines = sorted(set(line_numbers))
        where = (
            f"line {unique_lines[0]}"
            if len(unique_lines) == 1
            else "lines " + ", ".join(str(number) for number in unique_lines)
        )
        violations.append(
            ConectViolation(
                kind="duplicate_bond",
                serials=(origin, partner),
                atoms=(_identity_of(reference, origin), _identity_of(reference, partner)),
                distance=None,
                low=None,
                high=None,
                line_number=line_numbers[1],
                model=reference.label,
                message=(
                    f"The bond from atom serial {origin} "
                    f"({_identity_of(reference, origin)}) to {partner} "
                    f"({_identity_of(reference, partner)}) is recorded "
                    f"{len(line_numbers)} times from the same atom, on {where}. "
                    f"Recording it once from each atom is the PDB convention and is "
                    f"correct; recording it twice from the same atom is not, and doubles "
                    f"the bond order some viewers infer."
                ),
            )
        )

    # Reciprocity is counted over undirected bonds so that one missing back-reference is
    # one finding, not two.
    undirected: set[tuple[int, int]] = {
        (origin, partner) if origin < partner else (partner, origin)
        for origin, partner in directed
        if origin != partner
    }
    non_reciprocal = sorted(
        pair
        for pair in undirected
        if (pair[0], pair[1]) not in directed or (pair[1], pair[0]) not in directed
    )
    if require_reciprocal:
        for low_serial, high_serial in non_reciprocal:
            listed, unlisted = (
                (low_serial, high_serial)
                if (low_serial, high_serial) in directed
                else (high_serial, low_serial)
            )
            line_number = directed[listed, unlisted][0]
            violations.append(
                ConectViolation(
                    kind="non_reciprocal",
                    serials=(listed, unlisted),
                    atoms=(_identity_of(reference, listed), _identity_of(reference, unlisted)),
                    distance=None,
                    low=None,
                    high=None,
                    line_number=line_number,
                    model=reference.label,
                    message=(
                        f"Atom serial {listed} ({_identity_of(reference, listed)}) lists "
                        f"{unlisted} ({_identity_of(reference, unlisted)}) as bonded on "
                        f"line {line_number}, but {unlisted} does not list {listed}. Some "
                        f"viewers draw such a bond and some do not."
                    ),
                )
            )
    elif non_reciprocal:
        notes.append(
            f"{len(non_reciprocal)} of {len(undirected)} bond(s) are recorded from only "
            f"one of their two atoms. That is how DODO's writer emits every bond, so it "
            f"is not reported; pass require_reciprocal=True to check a file from another "
            f"writer."
        )

    usable = _check_serials(frames, directed, violations)
    bonds = sorted(pair for pair in undirected if pair[0] in usable and pair[1] in usable)

    for low_serial, high_serial in bonds:
        left = reference.atoms[reference.by_serial[low_serial]]
        right = reference.atoms[reference.by_serial[high_serial]]
        if left.chain_id == right.chain_id:
            continue
        line_number = directed.get(
            (low_serial, high_serial), directed.get((high_serial, low_serial), [0])
        )[0]
        violations.append(
            ConectViolation(
                kind="interchain_bond",
                serials=(low_serial, high_serial),
                atoms=(left.identity, right.identity),
                distance=None,
                low=None,
                high=None,
                line_number=line_number,
                model=reference.label,
                message=(
                    f"The CONECT record on line {line_number} bonds {left.identity} "
                    f"(serial {low_serial}, chain {left.chain_id or '?'}) to "
                    f"{right.identity} (serial {high_serial}, chain "
                    f"{right.chain_id or '?'}), which are in different chains. DODO never "
                    f"emits a bond across a chain break; if this file came from elsewhere "
                    f"it may be a genuine inter-chain link such as a disulfide."
                ),
            )
        )

    longest, shortest, longest_frame, shortest_frame = _measure(frames, bonds)
    for index in np.flatnonzero(longest > max_bond_length):
        position = int(index)
        low_serial, high_serial = bonds[position]
        frame = frames[longest_frame[position]]
        measured = float(longest[position])
        left = frame.atoms[frame.by_serial[low_serial]]
        right = frame.atoms[frame.by_serial[high_serial]]
        line_number = directed.get(
            (low_serial, high_serial), directed.get((high_serial, low_serial), [0])
        )[0]
        in_frame = "" if len(frames) == 1 else f" in {frame.label}"
        violations.append(
            ConectViolation(
                kind="bond_too_long",
                serials=(low_serial, high_serial),
                atoms=(left.identity, right.identity),
                distance=measured,
                low=None,
                high=max_bond_length,
                line_number=line_number,
                model=frame.label,
                message=(
                    f"The CONECT record on line {line_number} bonds {left.identity} "
                    f"(serial {low_serial}) to {right.identity} (serial {high_serial}), "
                    f"but they are {measured:.3f} A apart{in_frame}, beyond the "
                    f"{max_bond_length:.3f} A maximum -- the longest bond DODO writes is "
                    f"the {CA_CA_BOND_LENGTH:.2f} A CA-CA virtual bond. A viewer draws "
                    f"this as a straight line across the structure."
                ),
            )
        )
    for index in np.flatnonzero(shortest < min_bond_length):
        position = int(index)
        low_serial, high_serial = bonds[position]
        frame = frames[shortest_frame[position]]
        measured = float(shortest[position])
        left = frame.atoms[frame.by_serial[low_serial]]
        right = frame.atoms[frame.by_serial[high_serial]]
        line_number = directed.get(
            (low_serial, high_serial), directed.get((high_serial, low_serial), [0])
        )[0]
        in_frame = "" if len(frames) == 1 else f" in {frame.label}"
        coincident = " -- they are effectively the same point" if measured < 0.1 else ""
        violations.append(
            ConectViolation(
                kind="bond_too_short",
                serials=(low_serial, high_serial),
                atoms=(left.identity, right.identity),
                distance=measured,
                low=min_bond_length,
                high=None,
                line_number=line_number,
                model=frame.label,
                message=(
                    f"The CONECT record on line {line_number} bonds {left.identity} "
                    f"(serial {low_serial}) to {right.identity} (serial {high_serial}), "
                    f"but they are only {measured:.3f} A apart{in_frame}, below the "
                    f"{min_bond_length:.2f} A minimum for any covalent bond{coincident}."
                ),
            )
        )

    expected = _expected_bonds(reference, max_bond_length)
    present: set[tuple[int, int]] = {
        (origin, partner) if origin < partner else (partner, origin) for origin, partner in directed
    }
    missing = [
        (min(left, right), max(left, right), why, separation)
        for left, right, why, separation in expected
        if (min(left, right), max(left, right)) not in present
    ]
    n_expected_present = len(expected) - len(missing)

    enforce = (
        _declares_polymer_bonds(reference, undirected)
        if require_backbone_conect is None
        else require_backbone_conect
    )
    if missing and not enforce:
        notes.append(
            f"{len(missing)} backbone bond(s) implied by the coordinates have no CONECT "
            f"record, but this file declares no polymer backbone connectivity at all -- "
            f"it is conect=False output, or a deposition, which never lists backbone "
            f"bonds. They are counted in the completeness fraction and not reported. "
            f"Pass require_backbone_conect=True to report them."
        )
    if enforce:
        for low_serial, high_serial, why, separation in missing:
            left = reference.atoms[reference.by_serial[low_serial]]
            right = reference.atoms[reference.by_serial[high_serial]]
            virtual = "CA" in why
            violations.append(
                ConectViolation(
                    kind="missing_bond",
                    serials=(low_serial, high_serial),
                    atoms=(left.identity, right.identity),
                    distance=separation,
                    low=None,
                    high=max_bond_length,
                    line_number=left.line_number,
                    model=reference.label,
                    message=(
                        f"{left.identity} (serial {low_serial}) and {right.identity} "
                        f"(serial {high_serial}) are {separation:.3f} A apart, which is "
                        f"{why}, but no CONECT record bonds them."
                        + (
                            f" A {CA_CA_BOND_LENGTH:.2f} A CA-CA bond is past the "
                            f"automatic bond-detection cutoff in VMD and PyMOL, so this "
                            f"pair renders as disconnected dots."
                            if virtual
                            else " The atom will render detached from its own residue."
                        )
                    ),
                )
            )

    if parsed.n_extended_partners:
        notes.append(
            f"{parsed.n_extended_partners} partner serial(s) were read from columns "
            f"beyond the four the CONECT specification defines (columns 12-31) and "
            f"treated as bonds. In PDB format 2.2 and earlier columns 32-56 held "
            f"hydrogen-bonded and salt-bridged partners, which are not covalent bonds."
        )
    if parsed.n_records_inside_models:
        notes.append(
            f"{parsed.n_records_inside_models} CONECT record(s) appear inside a "
            f"MODEL/ENDMDL block. The specification puts them after the last ENDMDL, "
            f"which is where DODO writes them; the serials were checked regardless."
        )
    if parsed.n_conect_lines == 0:
        notes.append(
            "This file has no CONECT records. That is valid output "
            "(write_pdb(conect=False)), but a CA-only region in it will render as "
            "disconnected dots, because a 3.81 A CA-CA bond is past the automatic "
            "bond-detection cutoff in VMD and PyMOL."
        )

    violations.sort(key=lambda v: (v.serials[:1] or (0,), v.kind, v.serials[1:], v.line_number))
    return ConectReport(
        violations=tuple(violations),
        n_records=parsed.n_conect_lines,
        n_bonds=len(undirected),
        n_expected_bonds=len(expected),
        n_expected_present=n_expected_present,
        n_non_reciprocal=len(non_reciprocal),
        n_atoms=len(reference.atoms),
        n_models=len(frames),
        notes=tuple(notes),
        source=source,
    )


def validate_conect_file(
    path: str | Path,
    *,
    max_bond_length: float = MAX_BOND_LENGTH,
    min_bond_length: float = MIN_BOND_LENGTH,
    require_reciprocal: bool = False,
    require_backbone_conect: bool | None = None,
) -> ConectReport:
    """Validate the CONECT records of a PDB file on disk.

    Parameters
    ----------
    path
        The file to read. Decoded as latin-1, which cannot fail: the PDB format is ASCII,
        and a stray high byte in a header must not stop a connectivity check.
    max_bond_length, min_bond_length, require_reciprocal, require_backbone_conect
        See :func:`validate_conect_lines`.

    Returns
    -------
    ConectReport
        Truthy when the file is valid.

    Raises
    ------
    StructureFileError
        If the file cannot be read.
    EmptyStructureError
        If it has no ATOM or HETATM records.
    MalformedRecordError
        If an atom record's serial or coordinate columns cannot be parsed.

    Examples
    --------
    >>> report = validate_conect_file("out.pdb")   # doctest: +SKIP
    >>> print(report.describe())                   # doctest: +SKIP
    """
    target = Path(path)
    try:
        text = target.read_bytes().decode("latin-1")
    except OSError as error:
        raise StructureFileError(f"Cannot read {str(target)!r}: {error}.") from error
    return validate_conect_lines(
        text.splitlines(),
        source=target.name,
        max_bond_length=max_bond_length,
        min_bond_length=min_bond_length,
        require_reciprocal=require_reciprocal,
        require_backbone_conect=require_backbone_conect,
    )
