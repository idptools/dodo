"""Reader for the legacy PDB coordinate format.

Why fixed columns
-----------------
Every field below is read by column offset, never by ``str.split()``. That is not
pedantry: ``ATOM`` names are right-justified in a four-column field so that " CA " and
"CA  " mean different elements, residue names for nucleic acids are written "  A", and
coordinate fields run together with no separator as soon as one of them needs six
significant figures and a sign (``-100.123-200.456 -12.345``). A whitespace-splitting
reader gets a different number of fields per line on exactly the files this package
targets.

How models are framed
---------------------
A *frame* starts at a ``MODEL`` record, or at the first coordinate record that follows
an ``ENDMDL`` (or the start of the file) with no ``MODEL`` in between; it ends at the
next ``MODEL`` or ``ENDMDL``. Frames are numbered positionally from 1 and the serial in
columns 10-14 of the ``MODEL`` record is never read, because files written by modelling
scripts -- DODO's own v1 output among them -- put it in the wrong columns or number it
inconsistently. One consequence is worth stating plainly: coordinates appearing *before*
a file's first ``MODEL`` record are model 1, and that file's ``MODEL 1`` header then
opens model *2*. Such a file is malformed either way, and of the two readings the other
one merges records from different frames into a single residue -- which is the class of
silently impossible structure this rewrite exists to eliminate.

Every coordinate record belonging to a frame other than the requested one is counted and
reported in ``structure.notes``. So is every ``MODEL`` left unclosed by an ``ENDMDL`` and
every ``ENDMDL`` that closes nothing.

What this reader deliberately does not carry
--------------------------------------------
:class:`~dodo.structure.Structure` has no column for atom serial numbers or formal
charges, so neither is stored. The serial field is still *decoded* (see
:func:`decode_hybrid36`), because a file whose serials DODO cannot read is a file worth
mentioning in ``structure.notes`` -- but an unreadable serial is not fatal, since
nothing downstream uses it.
"""

from __future__ import annotations

import gzip
import re
from collections import Counter
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Final

import numpy as np

from ..constants import ATOMIC_MASSES, THREE_TO_ONE
from ..exceptions import (
    EmptyStructureError,
    MalformedRecordError,
    StructureFileError,
    UnsupportedFormatError,
)
from ..structure import Structure

__all__ = [
    "decode_hybrid36",
    "parse_pdb_lines",
    "read_pdb",
]

# ---------------------------------------------------------------------------
# Format constants
# ---------------------------------------------------------------------------

#: Field widths of the two PDB columns that overflow into hybrid-36 encoding.
_SERIAL_WIDTH: Final[int] = 5
_RESSEQ_WIDTH: Final[int] = 4

#: Minimum line length that can still hold a complete coordinate triple. Everything
#: past column 54 (occupancy, temperature factor, element, charge) is optional, and
#: real files -- especially anything hand-edited or written by a modelling script --
#: routinely stop at 54, 66 or 78 rather than padding to 80.
_MIN_COORDINATE_COLUMNS: Final[int] = 54

#: Elements dropped unless ``keep_hydrogens`` is set. DODO's geometry is heavy-atom
#: based throughout, and deuterium is a hydrogen for that purpose.
_HYDROGEN_ELEMENTS: Final[frozenset[str]] = frozenset({"H", "D"})

#: Residue names that are solvent no matter which record type they arrive on. Some
#: modelling tools write waters as ATOM rather than HETATM, and ``TIP``/``SOL`` come from
#: MD-derived structures. Checked ahead of the polymer whitelist so that a water written
#: as an ATOM record cannot become a residue in the middle of a chain. Three characters
#: only, because that is the width of the resName field.
_SOLVENT_NAMES: Final[frozenset[str]] = frozenset({"HOH", "DOD", "WAT", "TIP", "SOL"})

#: Marker appended to a chain label to keep two runs of the same chain id apart.
#:
_SEGMENT_MARKER: Final[str] = "\x01"

_GZIP_MAGIC: Final[bytes] = b"\x1f\x8b"

#: Accepted names, after any ``.gz`` has been peeled off. ``.pdb1``/``.pdb2`` are
#: biological-assembly files as served by the RCSB, and ``.ent`` is the name the PDB's
#: own FTP mirrors use.
_ACCEPTED_SUFFIX: Final[re.Pattern[str]] = re.compile(r"\.(pdb|ent)\d*\Z", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Hybrid-36
# ---------------------------------------------------------------------------


def decode_hybrid36(field: str, width: int) -> int:
    """Decode a PDB hybrid-36 integer field.

    Hybrid-36 is what real files use once a fixed-width integer column overflows, and
    the EM assemblies DODO targets overflow the five-column atom serial routinely: an
    ordinary ribosome or virus capsid entry passes 99,999 atoms. The encoding is plain
    decimal while the value fits, then unsigned base 36 with an uppercase leading
    digit, then base 36 with a lowercase leading digit. For a five-column field that
    means ``99999`` is followed by ``A0000`` = 100000 and, much later, ``a0000``; for a
    four-column field ``9999`` is followed by ``A000`` = 10000.

    Parameters
    ----------
    field
        The raw column text. Surrounding whitespace is ignored.
    width
        Width of the field in columns: 5 for an atom serial, 4 for a residue sequence
        number. The decimal/base-36 boundary depends on it.

    Returns
    -------
    int
        The decoded value. Decimal fields may be negative, which is legal (and
        common) for residue sequence numbers; base-36 fields never are.

    Raises
    ------
    ValueError
        If ``width`` is not positive, or the field is empty, or it is neither a decimal
        integer nor a well-formed base-36 value of exactly ``width`` characters.

    Examples
    --------
    >>> decode_hybrid36("99999", 5), decode_hybrid36("A0000", 5)
    (99999, 100000)
    >>> decode_hybrid36("9999", 4), decode_hybrid36("A000", 4)
    (9999, 10000)
    """
    if width < 1:
        raise ValueError(f"Field width must be positive, got {width}.")
    text = field.strip()
    if not text:
        raise ValueError("Cannot decode an empty hybrid-36 field.")

    first = text[0]
    if first == "-" or first.isdigit():
        try:
            return int(text)
        except ValueError:
            raise ValueError(f"{text!r} is not a decimal integer.") from None

    if len(text) != width:
        raise ValueError(
            f"A base-36 hybrid-36 value must fill all {width} columns, but {text!r} "
            f"has {len(text)}."
        )
    # The powers are wrapped in int() because typeshed types int.__pow__ as returning
    # Any (a negative exponent yields a float), which would leak Any into the result.
    digit_span = int(36 ** (width - 1))
    # The case of the *leading* character selects the range, so a mixed-case field is
    # ambiguous rather than merely unusual: int(x, 36) is case-insensitive and would
    # silently decode "Aa000" as if it were "AA000".
    if first.isupper():
        if any(character.islower() for character in text):
            raise ValueError(f"Mixed-case hybrid-36 value {text!r} is ambiguous.")
        case_offset = 0
    elif first.islower():
        if any(character.isupper() for character in text):
            raise ValueError(f"Mixed-case hybrid-36 value {text!r} is ambiguous.")
        case_offset = 26 * digit_span
    else:
        raise ValueError(f"{text!r} is not a valid hybrid-36 value.")

    try:
        raw = int(text, 36)
    except ValueError:
        raise ValueError(f"{text!r} is not a valid base-36 value.") from None
    return raw - 10 * digit_span + int(10**width) + case_offset


# ---------------------------------------------------------------------------
# Element inference
# ---------------------------------------------------------------------------


def _infer_element(name_field: str) -> str:
    """Infer an element symbol from an atom-name field that has no element column.

    Older depositions and anything hand-edited or written by a modelling script often
    stop the line at column 66 or 78, leaving columns 76-78 empty. Guessing here beats
    handing the structure on with a blank element, because the element drives both
    hydrogen filtering and the atomic-mass lookup.

    Parameters
    ----------
    name_field
        Columns 12-16 of the record, sliced but *not* stripped: the leading space is
        the only thing distinguishing an alpha carbon from a calcium ion.

    Returns
    -------
    str
        Uppercase element symbol, or ``""`` if the field contains no letters.
    """
    field = name_field.ljust(4)
    stripped = field.strip()
    if not stripped:
        return ""
    if field[0] == " ":
        # Element symbols are right-justified in columns 12-13, so a name that leaves
        # column 12 blank has a one-letter element in column 13: " CA " is C, "CA  "
        # is calcium.
        if field[1].isalpha():
            return field[1].upper()
    else:
        # A fully occupied name field is either a two-letter element -- selenium in a
        # selenomethionine is written "SE  " -- or a four-character hydrogen name such
        # as "HG12" or the old-style digit-first "1HG2". Take the two-letter reading
        # only when it names an element DODO has data for; that keeps MSE's selenium
        # right without inventing elements for files that left-justify their names.
        candidate = field[:2].strip().upper()
        if candidate in ATOMIC_MASSES:
            return candidate
    for character in stripped:
        if character.isalpha():
            return character.upper()
    return ""


# ---------------------------------------------------------------------------
# Alternate conformers
# ---------------------------------------------------------------------------


def _alternate_location_mask(
    segment_labels: Sequence[str],
    residue_numbers: Sequence[int],
    insertion_codes: Sequence[str],
    alt_locs: Sequence[str],
) -> np.ndarray:
    """Select one conformer per residue from records carrying altLoc indicators.

    Parameters
    ----------
    segment_labels
        Per-atom chain label, TER-splitting already applied.
    residue_numbers
        Per-atom decoded residue sequence number.
    insertion_codes
        Per-atom insertion code, ``""`` when absent.
    alt_locs
        Per-atom altLoc character, ``""`` when absent.

    Returns
    -------
    np.ndarray
        Boolean mask over the records, ``True`` for the atoms to keep.
    """
    n = len(alt_locs)
    labels = np.asarray(segment_labels, dtype="<U4")
    numbers = np.asarray(residue_numbers, dtype=np.int64)
    icodes = np.asarray(insertion_codes, dtype="<U1")

    changed = np.ones(n, dtype=bool)
    changed[1:] = (
        (labels[1:] != labels[:-1]) | (numbers[1:] != numbers[:-1]) | (icodes[1:] != icodes[:-1])
    )
    group = np.cumsum(changed) - 1
    n_groups = int(group[-1]) + 1

    codes = np.fromiter((ord(alt) if alt else 0 for alt in alt_locs), dtype=np.int64, count=n)
    blank = codes == 0
    # Alternates are written in descending occupancy, so 'A' is the one to keep;
    # ranking it below every other code makes a per-residue minimum pick it whenever
    # it is present. When it is *not* -- a residue modelled only as 'B' happens in
    # real depositions -- the minimum falls back to the first code that is there.
    # Dropping such a residue outright would fabricate exactly the phantom chain break
    # that the MSE bug used to produce.
    sentinel = np.iinfo(np.int64).max
    rank = np.where(codes == ord("A"), 1, codes)
    preferred = np.full(n_groups, sentinel, dtype=np.int64)
    np.minimum.at(preferred, group, np.where(blank, sentinel, rank))
    mask: np.ndarray = blank | (rank == preferred[group])
    return mask


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def read_pdb(
    path: str | Path,
    *,
    model: int | None = None,
    keep_hydrogens: bool = False,
) -> Structure:
    """Read a PDB-format file into a :class:`~dodo.structure.Structure`.

    Parameters
    ----------
    path
        Path to a ``.pdb``, ``.PDB``, ``.ent`` or biological-assembly ``.pdb1`` file,
        optionally gzipped. Gzip is detected from the file's magic bytes as well as
        from the ``.gz`` suffix, so a compressed file that was renamed still reads.
    model
        Which model to read, 1-based. ``None`` (the default) reads the first.
    keep_hydrogens
        Keep hydrogen and deuterium atoms. Off by default because every geometric
        operation in DODO is heavy-atom based. The count is recorded in
        ``structure.notes`` either way.

    Returns
    -------
    Structure
        The parsed structure, with ``notes`` describing anything skipped, inferred or
        deduplicated along the way.

    Raises
    ------
    UnsupportedFormatError
        If the extension is not one of the PDB-format extensions.
    StructureFileError
        If the file cannot be read or decompressed, or ``model`` does not exist.
    MalformedRecordError
        If an atom record is too short to hold coordinates, or its coordinates or
        residue sequence number cannot be parsed.
    EmptyStructureError
        If the file holds no polymer atoms at all.
    """
    file_path = Path(path)
    name = file_path.name
    if name.lower().endswith(".gz"):
        name = name[: -len(".gz")]
    if not _ACCEPTED_SUFFIX.search(name):
        raise UnsupportedFormatError(
            f"{file_path.name!r} is not a PDB-format file name. Accepted: .pdb, .PDB, "
            f".ent, .pdb1 (and other assembly numbers), each optionally .gz. mmCIF "
            f"files need DODO's mmCIF reader instead."
        )

    try:
        raw = file_path.read_bytes()
    except OSError as error:
        raise StructureFileError(f"Could not read {file_path}: {error}") from error

    if raw[:2] == _GZIP_MAGIC:
        try:
            raw = gzip.decompress(raw)
        except (OSError, EOFError) as error:
            raise StructureFileError(
                f"{file_path} looks gzipped but could not be decompressed: {error}"
            ) from error

    # latin-1, not utf-8: this is a fixed-column format, and latin-1 is the only codec
    # that guarantees one character per byte. A single non-ASCII byte in an author name
    # in a TITLE or JRNL record would otherwise shift every subsequent column of the
    # decoded text, silently corrupting coordinates.
    text = raw.decode("latin-1")
    return parse_pdb_lines(
        text.splitlines(),
        source=str(file_path),
        model=model,
        keep_hydrogens=keep_hydrogens,
    )


def parse_pdb_lines(
    lines: Iterable[str],
    *,
    source: str | None = None,
    model: int | None = None,
    keep_hydrogens: bool = False,
) -> Structure:
    """Parse PDB-format records into a :class:`~dodo.structure.Structure`.

    Parameters
    ----------
    lines
        The records, one per string. Line terminators are tolerated. Consumed once,
        so a file handle or generator is fine.
    source
        Provenance string stored on the structure, for error messages.
    model
        Which model to read, 1-based. ``None`` reads the first. Models are counted
        positionally, in file order; see the module docstring for how frames are
        delimited when ``MODEL``/``ENDMDL`` records are missing or unbalanced.
    keep_hydrogens
        Keep hydrogen and deuterium atoms rather than dropping them.

    Returns
    -------
    Structure
        The parsed structure. ``Chain.uniprot_id`` is filled from ``DBREF`` records
        naming UniProt, and ``Chain.full_sequence`` from ``SEQRES``.

    Raises
    ------
    StructureFileError
        If ``model`` is not positive, or exceeds the number of models present.
    MalformedRecordError
        If an atom record is too short to hold coordinates, or its coordinates or
        residue sequence number cannot be parsed.
    EmptyStructureError
        If no polymer atoms were found.
    """
    target_model = 1 if model is None else int(model)
    if target_model < 1:
        raise StructureFileError(f"Models are numbered from 1; got model={model!r}.")

    xyz_rows: list[tuple[float, float, float]] = []
    atom_names: list[str] = []
    elements: list[str] = []
    residue_names: list[str] = []
    residue_numbers: list[int] = []
    insertion_codes: list[str] = []
    occupancies: list[float] = []
    b_factors: list[float] = []
    segment_labels: list[str] = []
    alt_locs: list[str] = []

    seqres_codes: dict[str, list[str]] = {}
    seqres_declared: dict[str, int] = {}
    dbref_hits: list[tuple[str, str]] = []
    dbref1_uniprot_chains: set[str] = set()

    skipped_names: Counter[str] = Counter()
    unrecognized_polymer_names: Counter[str] = Counter()
    n_hydrogens = 0
    n_unreadable_serials = 0
    n_unreadable_optional_fields = 0
    n_inferred_elements = 0
    n_ter_splits = 0

    # Model framing. ``n_frames`` counts every frame seen so far, whether it was opened
    # by a MODEL record or implicitly by a coordinate record outside one; ``open_frame``
    # is the 1-based number of the frame currently open, or 0 when none is. A file with
    # no MODEL record at all therefore has exactly one, implicit, model.
    n_frames = 0
    open_frame = 0
    n_other_model_atoms = 0
    n_unclosed_models = 0
    n_stray_endmdl = 0

    current_label: str | None = None
    current_chain_id: str | None = None
    break_pending = False

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.rstrip("\r\n")
        record = line[:6]

        if record.startswith(("ATOM", "HETATM")):
            if open_frame == 0:
                # A coordinate record with no model open opens one: the whole of a file
                # that has no MODEL records, a frame separated from the previous one by
                # a bare ENDMDL, or coordinates preceding the first MODEL record.
                n_frames += 1
                open_frame = n_frames
            if open_frame != target_model:
                # Counted, never merely dropped: an atom that leaves this reader
                # unaccounted for is the failure mode the rewrite exists to remove.
                n_other_model_atoms += 1
                continue
            if len(line) < _MIN_COORDINATE_COLUMNS:
                raise MalformedRecordError(
                    f"Atom record is only {len(line)} columns long; coordinates need at "
                    f"least {_MIN_COORDINATE_COLUMNS}",
                    line_number,
                    line,
                )

            residue_name = line[17:20].strip().upper()
            is_hetatm = record.startswith("HETATM")
            # THREE_TO_ONE is the polymer whitelist. Accepting HETATM residues that
            # are in it is what stops a mid-chain selenomethionine from vanishing and
            # fabricating a chain break; rejecting the ones that are not is what keeps
            # waters and ligands out of the polymer.
            if residue_name in _SOLVENT_NAMES or (is_hetatm and residue_name not in THREE_TO_ONE):
                skipped_names[residue_name or "(blank)"] += 1
                continue
            if not is_hetatm and residue_name not in THREE_TO_ONE:
                # Kept, not skipped: an ATOM record is a polymer record by definition,
                # and dropping it would break the chain. It will read as 'X' in the
                # sequence, which is honest.
                unrecognized_polymer_names[residue_name or "(blank)"] += 1

            try:
                decode_hybrid36(line[6:11], _SERIAL_WIDTH)
            except ValueError:
                # Nothing downstream uses the serial, so an unreadable one is worth a
                # note rather than a refusal to read the file.
                n_unreadable_serials += 1

            resseq_field = line[22:26]
            try:
                residue_number = int(resseq_field)
            except ValueError:
                try:
                    residue_number = decode_hybrid36(resseq_field, _RESSEQ_WIDTH)
                except ValueError as error:
                    raise MalformedRecordError(
                        f"Unreadable residue sequence number {resseq_field!r}: {error}",
                        line_number,
                        line,
                    ) from error

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError as error:
                raise MalformedRecordError(
                    f"Unreadable coordinates: {error}", line_number, line
                ) from error

            element = line[76:78].strip().upper()
            if not element:
                element = _infer_element(line[12:16])
                n_inferred_elements += 1
            if element in _HYDROGEN_ELEMENTS:
                n_hydrogens += 1
                if not keep_hydrogens:
                    continue

            occupancy = 1.0
            occupancy_field = line[54:60].strip()
            if occupancy_field:
                try:
                    occupancy = float(occupancy_field)
                except ValueError:
                    n_unreadable_optional_fields += 1
            b_factor = 0.0
            b_factor_field = line[60:66].strip()
            if b_factor_field:
                try:
                    b_factor = float(b_factor_field)
                except ValueError:
                    n_unreadable_optional_fields += 1

            chain_id = line[21:22].strip()
            if current_label is None or chain_id != current_chain_id or break_pending:
                if break_pending and current_chain_id == chain_id:
                    # A TER separated two runs of one chain id. Alternate the marker so
                    # the new run cannot merge into the run that just ended -- either
                    # direction of the toggle is a distinct label, and counting only the
                    # unmarked-to-marked direction is what used to halve this count.
                    marked = chain_id + _SEGMENT_MARKER
                    label = chain_id if current_label == marked else marked
                    n_ter_splits += 1
                else:
                    label = chain_id
                current_label = label
                current_chain_id = chain_id
                break_pending = False

            xyz_rows.append((x, y, z))
            atom_names.append(line[12:16].strip())
            elements.append(element)
            residue_names.append(residue_name)
            residue_numbers.append(residue_number)
            insertion_codes.append(line[26:27].strip())
            occupancies.append(occupancy)
            b_factors.append(b_factor)
            segment_labels.append(current_label)
            alt_locs.append(line[16:17].strip())

        elif record.startswith("TER"):
            # A TER ends a chain even when the next one reuses the id. Ignoring it is
            # how two chains both labelled 'A' used to be merged into one.
            if open_frame == target_model:
                break_pending = True

        elif record.startswith("MODEL"):
            # Counted positionally rather than read from the serial in columns 10-14:
            # files written by modelling scripts contain things like "MODEL 1" with the
            # number left-justified, and DODO's own multi-model outputs are among them.
            if open_frame != 0:
                # No ENDMDL closed the previous frame. Letting this MODEL end it is an
                # inference -- the alternative, treating both blocks as one model, is
                # the merge this reader must never perform -- so it is recorded.
                n_unclosed_models += 1
            n_frames += 1
            open_frame = n_frames
            break_pending = False

        elif record.startswith("ENDMDL"):
            if open_frame == 0:
                n_stray_endmdl += 1
            open_frame = 0
            break_pending = False

        elif record.startswith("SEQRES"):
            seqres_chain = line[11:12].strip()
            cells = line[19:70]
            names = [cells[i : i + 3].strip().upper() for i in range(0, len(cells), 4)]
            seqres_codes.setdefault(seqres_chain, []).extend(
                THREE_TO_ONE.get(name, "X") for name in names if name
            )
            declared = line[13:17].strip()
            if declared.isdigit():
                seqres_declared[seqres_chain] = int(declared)

        elif record.startswith("DBREF"):
            tag = record.strip()
            dbref_chain = line[12:13].strip()
            if tag == "DBREF":
                if line[26:32].strip().upper() == "UNP":
                    dbref_hits.append((dbref_chain, line[33:41].strip()))
            elif tag == "DBREF1":
                # DBREF1/DBREF2 pairs carry accessions too long for the DBREF layout;
                # the database name is on the DBREF1 line and the accession on DBREF2.
                if line[26:32].strip().upper() == "UNP":
                    dbref1_uniprot_chains.add(dbref_chain)
            elif tag == "DBREF2" and dbref_chain in dbref1_uniprot_chains:
                dbref_hits.append((dbref_chain, line[18:40].strip()))

    total_models = max(n_frames, 1)
    if target_model > total_models:
        where = f" in {source}" if source else ""
        raise StructureFileError(
            f"Model {target_model} was requested but there are only {total_models} model(s){where}."
        )

    notes: list[str] = []
    if not xyz_rows:
        detail = (
            f" {sum(skipped_names.values())} record(s) were skipped as solvent or "
            f"ligand ({_summarize(skipped_names)})."
            if skipped_names
            else ""
        )
        dropped_hydrogens = n_hydrogens and not keep_hydrogens
        hydrogen_detail = (
            f" {n_hydrogens} hydrogen atom(s) were dropped." if dropped_hydrogens else ""
        )
        # An empty requested model in a file that does have coordinates elsewhere is
        # much more often a wrong ``model=`` argument than an empty file, so say so.
        model_detail = (
            f" {n_other_model_atoms} atom record(s) belong to the other "
            f"{total_models - 1} model(s) of this file."
            if n_other_model_atoms
            else ""
        )
        raise EmptyStructureError(
            f"No polymer atoms found in model {target_model} of {source or 'the input'}."
            f"{detail}{hydrogen_detail}{model_detail}"
        )

    n_alternates_dropped = 0
    if any(alt_locs):
        keep = _alternate_location_mask(segment_labels, residue_numbers, insertion_codes, alt_locs)
        n_alternates_dropped = int(np.count_nonzero(~keep))
        if n_alternates_dropped:
            kept: list[int] = np.flatnonzero(keep).tolist()
            xyz_rows = [xyz_rows[i] for i in kept]
            atom_names = [atom_names[i] for i in kept]
            elements = [elements[i] for i in kept]
            residue_names = [residue_names[i] for i in kept]
            residue_numbers = [residue_numbers[i] for i in kept]
            insertion_codes = [insertion_codes[i] for i in kept]
            occupancies = [occupancies[i] for i in kept]
            b_factors = [b_factors[i] for i in kept]
            segment_labels = [segment_labels[i] for i in kept]

    structure = Structure.from_atom_records(
        xyz=np.asarray(xyz_rows, dtype=np.float64),
        atom_name=atom_names,
        element=elements,
        residue_name=residue_names,
        residue_number=residue_numbers,
        chain_id=segment_labels,
        insertion_code=insertion_codes,
        b_factor=b_factors,
        occupancy=occupancies,
        source=source,
    )

    for chain in structure.chains:
        if chain.chain_id.endswith(_SEGMENT_MARKER):
            chain.chain_id = chain.chain_id[: -len(_SEGMENT_MARKER)]

    notes.extend(
        _annotate_chains(
            structure,
            seqres_codes=seqres_codes,
            seqres_declared=seqres_declared,
            dbref_hits=dbref_hits,
        )
    )

    # Reported whenever more than one frame exists *or* any record was passed over,
    # never gated on having seen a MODEL record: a file whose frames are separated by
    # bare ENDMDL records has no MODEL records at all, and gating on them is how every
    # atom after the first ENDMDL used to vanish with an empty notes list.
    if total_models > 1 or n_other_model_atoms:
        ignored = (
            f"{n_other_model_atoms} atom record(s) belonging to the other model(s) were ignored"
            if n_other_model_atoms
            else "the other model(s) hold no atom records"
        )
        notes.append(f"Read model {target_model} of {total_models}; {ignored}.")
    if n_unclosed_models:
        notes.append(
            f"{n_unclosed_models} MODEL record(s) were not closed by an ENDMDL; each was "
            f"taken to end the model that preceded it rather than continue it."
        )
    if n_stray_endmdl:
        notes.append(f"{n_stray_endmdl} ENDMDL record(s) closed no open model and were ignored.")
    if skipped_names:
        notes.append(
            f"Skipped {sum(skipped_names.values())} record(s) whose residue name is not "
            f"a recognized polymer residue: {_summarize(skipped_names)}."
        )
    if unrecognized_polymer_names:
        notes.append(
            f"Kept {sum(unrecognized_polymer_names.values())} ATOM record(s) with an "
            f"unrecognized residue name ({_summarize(unrecognized_polymer_names)}); they "
            f"read as 'X' in the sequence."
        )
    if n_hydrogens:
        notes.append(
            f"{n_hydrogens} hydrogen/deuterium atom(s) "
            f"{'kept' if keep_hydrogens else 'dropped (DODO geometry is heavy-atom based)'}."
        )
    if n_alternates_dropped:
        notes.append(
            f"Dropped {n_alternates_dropped} alternate-conformer atom(s); one conformer "
            f"per residue was kept, preferring altLoc 'A'."
        )
    if n_ter_splits:
        notes.append(
            f"{n_ter_splits} TER record(s) separated two runs of the same chain id; they "
            f"were kept as separate chains rather than merged."
        )
    if n_inferred_elements:
        notes.append(
            f"Inferred the element of {n_inferred_elements} atom(s) from the atom name "
            f"because columns 76-78 were empty."
        )
    if n_unreadable_serials:
        notes.append(
            f"{n_unreadable_serials} atom serial(s) could not be decoded as decimal or "
            f"hybrid-36. Serials are not used by DODO, so the atoms were kept."
        )
    if n_unreadable_optional_fields:
        notes.append(
            f"{n_unreadable_optional_fields} occupancy/temperature-factor field(s) were "
            f"unreadable; defaults (1.0 / 0.0) were used."
        )
    structure.notes.extend(notes)
    return structure


def _annotate_chains(
    structure: Structure,
    *,
    seqres_codes: dict[str, list[str]],
    seqres_declared: dict[str, int],
    dbref_hits: list[tuple[str, str]],
) -> list[str]:
    """Attach SEQRES sequences and UniProt accessions to a structure's chains.

    Parameters
    ----------
    structure
        The structure whose chains are annotated, in place.
    seqres_codes
        One-letter deposited sequence per chain id, harvested from SEQRES.
    seqres_declared
        Residue count each chain's SEQRES header claims, for cross-checking.
    dbref_hits
        ``(chain_id, accession)`` pairs from DBREF records naming UniProt, in file
        order.

    Returns
    -------
    list[str]
        Notes describing any inconsistency found while annotating.
    """
    notes: list[str] = []

    uniprot_by_chain: dict[str, str] = {}
    conflicting: list[str] = []
    for chain_id, accession in dbref_hits:
        if not accession:
            continue
        existing = uniprot_by_chain.setdefault(chain_id, accession)
        if existing != accession and chain_id not in conflicting:
            conflicting.append(chain_id)
    if conflicting:
        notes.append(
            f"Chain(s) {', '.join(sorted(conflicting))} have DBREF records naming more "
            f"than one UniProt accession (a chimeric construct); the first was used."
        )

    full_sequences = {chain_id: "".join(codes) for chain_id, codes in seqres_codes.items()}
    for chain_id, sequence in full_sequences.items():
        declared = seqres_declared.get(chain_id)
        if declared is not None and declared != len(sequence):
            notes.append(
                f"SEQRES for chain {chain_id!r} declares {declared} residues but lists "
                f"{len(sequence)}; the listed residues were used."
            )

    present = {chain.chain_id for chain in structure.chains}
    for chain in structure.chains:
        # The *deposited* sequence, deliberately not the UniProt canonical isoform:
        # measuring missing residues against the isoform invents phantom termini for
        # every tagged or truncated construct.
        if chain.chain_id in full_sequences:
            chain.full_sequence = full_sequences[chain.chain_id]
        if chain.chain_id in uniprot_by_chain:
            chain.uniprot_id = uniprot_by_chain[chain.chain_id]

    unmodelled = sorted(set(full_sequences) - present)
    if unmodelled:
        notes.append(
            f"SEQRES describes chain(s) {', '.join(repr(c) for c in unmodelled)} that "
            f"have no coordinates in this model."
        )
    return notes


def _summarize(counts: Counter[str], limit: int = 5) -> str:
    """Render a residue-name histogram compactly for a provenance note.

    Parameters
    ----------
    counts
        Residue name to record count.
    limit
        How many of the most common names to name individually.

    Returns
    -------
    str
        Something like ``"HOH x412, ADP x2 and 3 more"``.
    """
    head = ", ".join(f"{name} x{count}" for name, count in counts.most_common(limit))
    remainder = len(counts) - limit
    return f"{head} and {remainder} more" if remainder > 0 else head
