"""PDB writing.

Why this module is treated as regression-critical
-------------------------------------------------
DODO's product is a picture. Everything it computes reaches the user through a written
file, so a defect here erases the value of everything upstream. Both previous writers
had defects of exactly that kind:

* Neither emitted ``CONECT`` records. The CA-CA virtual bond is 3.80 A
  (:data:`~dodo.constants.CA_CA_BOND_LENGTH`), which is beyond the automatic
  bond-detection cutoff in VMD and PyMOL, so a rebuilt CA-only IDR loaded as a cloud of
  disconnected dots. That is the entire deliverable, lost at the last step.
* Both got the atom-name columns wrong, in opposite directions. v1 right-justified the
  name from column 13 and wrote no element column, so ``CA`` landed in columns 13-14 and
  MDTraj read a 912-residue CA trace as 912 *calcium* atoms.
  ``tests/data/structures/testing_translation.pdb`` is v1 output and still shows it.
* The multi-model writer's entire body was ``pass``.
* ``os.path.dirname("out.pdb")`` is ``""``, which the v1 directory check rejected, so a
  bare relative filename could not be written at all.
* ``float(monomer.occupancy)`` crashed when occupancy was absent.

Each of those has a named test in ``tests/unit/test_io_write.py``.
"""

from __future__ import annotations

import string
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..constants import (
    BACKBONE_ATOMS,
    BETA_DISORDERED,
    BETA_FOLDED,
    C_N_PEPTIDE_BOND_LENGTH,
    CA_CA_BOND_LENGTH,
    MAX_PDB_ATOM_SERIAL,
    MAX_PDB_RESIDUE_NUMBER,
    RESIDUES_PER_SEQRES_LINE,
    THREE_TO_ONE,
)
from ..exceptions import (
    EmptyStructureError,
    GeometryError,
    InvalidRegionError,
    StructureFileError,
)
from ..structure import DomainKind, Structure

__all__ = [
    "structure_to_pdb_lines",
    "write_pdb",
]

# ---------------------------------------------------------------------------
# Format limits and local tuning knobs
#
# Nothing here duplicates constants.py: these are properties of the PDB record
# format itself, or values derived from constants.py.
# ---------------------------------------------------------------------------

#: Field widths, derived from the format limits declared in constants.py so the two
#: cannot drift apart.
_SERIAL_WIDTH: int = len(str(MAX_PDB_ATOM_SERIAL))
_RESSEQ_WIDTH: int = len(str(MAX_PDB_RESIDUE_NUMBER))

#: Range of the ``%8.3f`` coordinate fields. Outside this the field overflows into its
#: neighbour and every column after it shifts, which is how a writer produces a file
#: that loads without error and is wrong.
_COORD_MIN: float = -999.999
_COORD_MAX: float = 9999.999

#: Range of the ``%6.2f`` occupancy and B-factor fields.
_REAL_MIN: float = -99.99
_REAL_MAX: float = 999.99

#: Defaults for missing occupancy / B-factor, matching what every PDB-writing tool
#: emits for "no information".
_DEFAULT_OCCUPANCY: float = 1.00
_DEFAULT_B_FACTOR: float = 0.00

#: Largest CA(i)-CA(i+1) separation that still gets a CONECT record, in Angstroms.
#:
_MAX_BONDED_CA_CA_DISTANCE: float = 1.3 * CA_CA_BOND_LENGTH

#: Largest C(i)-N(i+1) separation that counts as a peptide bond, in Angstroms. Used only
#: at junctions where a CA is missing, so the CA test cannot be applied.
_MAX_BONDED_C_N_DISTANCE: float = 1.5 * C_N_PEPTIDE_BOND_LENGTH

#: Bonded partners per CONECT record, per the PDB specification.
_CONECT_PARTNERS_PER_LINE: int = 4

#: Chain identifiers available for the single character the PDB format allows, used when
#: a structure carries multi-character mmCIF-style ids that will not fit.
_CHAIN_ID_POOL: str = string.ascii_uppercase + string.ascii_lowercase + string.digits

_HYBRID36_DIGITS_UPPER: str = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
_HYBRID36_DIGITS_LOWER: str = _HYBRID36_DIGITS_UPPER.lower()


# ---------------------------------------------------------------------------
# Field encoding
# ---------------------------------------------------------------------------


def _encode_pure(digits: str, value: int) -> str:
    """Encode a non-negative integer in the base implied by ``digits``."""
    if value == 0:
        return digits[0]
    base = len(digits)
    out: list[str] = []
    while value:
        value, remainder = divmod(value, base)
        out.append(digits[remainder])
    return "".join(reversed(out))


def _hybrid36(value: int, width: int, max_decimal: int) -> str:
    """Encode ``value`` in ``width`` columns, using hybrid-36 above ``max_decimal``.

    Hybrid-36 is the encoding real PDB files use once the decimal field is exhausted:
    the value stays decimal while it fits, then switches to upper-case base-36 and
    finally to lower-case base-36. It is what wwPDB tools, cctbx and MDTraj all read
    back, so it is the only honest way to write a structure larger than the format's
    nominal limit.

    Parameters
    ----------
    value
        The number to encode.
    width
        Field width in columns.
    max_decimal
        Largest value written in plain decimal, from :mod:`dodo.constants`.

    Returns
    -------
    str
        Exactly ``width`` characters.

    Raises
    ------
    StructureFileError
        If the value cannot be represented at all. Raising beats writing a field that
        overflows into its neighbour, which shifts every column after it.
    """
    if 1 - 10 ** (width - 1) <= value <= max_decimal:
        return f"{value:{width}d}"
    shifted = value - (max_decimal + 1)
    span = 26 * 36 ** (width - 1)
    offset = 10 * 36 ** (width - 1)
    if 0 <= shifted < span:
        return _encode_pure(_HYBRID36_DIGITS_UPPER, shifted + offset)
    shifted -= span
    if 0 <= shifted < span:
        return _encode_pure(_HYBRID36_DIGITS_LOWER, shifted + offset)
    raise StructureFileError(
        f"{value} cannot be written in a {width}-column PDB field, even with "
        f"hybrid-36 encoding. Refusing to write a corrupt file."
    )


def _format_atom_name(name: str, element: str) -> str:
    """Format an atom name into columns 13-16 (4 characters, no padding removed).

    The rule, and it is the one both previous writers got wrong in opposite directions:
    the *element symbol* is right-justified in columns 13-14 and the rest of the name
    follows. In practice that means a 1-3 character name on a one-letter element starts
    in column 14 (``" CA "``, ``" N  "``, ``" OD1"``), a 4-character name starts in
    column 13 (``"HD11"``), and a two-letter element starts in column 13 too
    (``"SE  "`` for selenomethionine, ``"ZN  "`` for a zinc ion -- both appear in the
    EM assemblies DODO passes through). Viewers infer the element from these columns,
    which is why v1's right-justified ``"CA"`` loaded as calcium.
    """
    name = name.strip()
    element = element.strip().upper()
    if len(name) >= 4:
        return name[:4]
    if len(element) == 2 and name.upper().startswith(element):
        return name.ljust(4)
    return " " + name.ljust(3)


def _format_real(value: float, default: float) -> tuple[str, bool]:
    """Format an occupancy or B-factor into a 6-column ``%6.2f`` field.

    Returns the formatted field and whether the value had to be substituted or clamped,
    so the caller can record that rather than let it pass unremarked.
    """
    substituted = False
    if not np.isfinite(value):
        # A Structure built from a file with no occupancy column carries NaN here.
        # v1 called float(None) and crashed; writing the conventional default is right,
        # but it is still a substitution and the caller says so.
        value = default
        substituted = True
    if value < _REAL_MIN or value > _REAL_MAX:
        value = min(max(value, _REAL_MIN), _REAL_MAX)
        substituted = True
    return f"{value:6.2f}", substituted


# ---------------------------------------------------------------------------
# Per-structure preparation
# ---------------------------------------------------------------------------


def _note(structure: Structure, message: str) -> None:
    """Record an anomaly on the structure, without duplicating it across frames."""
    if message not in structure.notes:
        structure.notes.append(message)


def _fits_the_chain_id_column(chain_id: str) -> bool:
    """Whether ``chain_id`` can be written in the PDB chain-id column as it stands.

    One printable ASCII character. A blank id is not one character and so does not fit,
    which is why an unset id gets a generated one.
    """
    return len(chain_id) == 1 and chain_id.isascii() and chain_id.isprintable()


def _chain_ids(structure: Structure) -> list[str]:
    """Return one single-character PDB chain id per chain index.

    mmCIF chain ids can be several characters (``"AAA"``); PDB has exactly one column.
    Rather than truncate and silently merge two chains, remap the whole set onto
    single characters and record the mapping.

    The column also has to hold one *printable ASCII* character: PDB is an ASCII
    format, and an mmCIF ``label_asym_id`` is not guaranteed to be. A one-character
    non-ASCII id therefore gets remapped too, rather than travelling all the way down
    to the encoder and surfacing as a ``UnicodeEncodeError`` from a module whose
    contract is that everything it cannot represent raises
    :class:`~dodo.exceptions.StructureFileError`.
    """
    n_chains = int(structure.chain_index.max()) + 1 if structure.n_residues else 0
    raw = [
        structure.chains[i].chain_id if i < len(structure.chains) else "" for i in range(n_chains)
    ]
    if len(structure.chains) < n_chains:
        _note(
            structure,
            f"chain_index refers to {n_chains} chains but only {len(structure.chains)} "
            f"Chain views exist; the remaining chains were given generated ids.",
        )
    if all(_fits_the_chain_id_column(chain_id) for chain_id in raw):
        return raw

    if n_chains > len(_CHAIN_ID_POOL):
        raise StructureFileError(
            f"This structure has {n_chains} chains with ids that do not fit the PDB "
            f"format's single chain-id column, and only {len(_CHAIN_ID_POOL)} "
            f"single-character ids exist. Write mmCIF instead."
        )
    remapped = list(_CHAIN_ID_POOL[:n_chains])
    changed = ", ".join(
        f"{old or '(unset)'}->{new}" for old, new in zip(raw, remapped, strict=True) if old != new
    )
    _note(
        structure,
        f"Chain ids were remapped to the single printable ASCII character the PDB "
        f"chain-id column holds: {changed}.",
    )
    return remapped


def _residue_b_factors(structure: Structure, *, annotate_regions: bool) -> np.ndarray:
    """Return the per-residue B-factor column to write.

    With ``annotate_regions``, the column is replaced by the region annotation. Note
    the polarity, which is the one thing about this feature people get wrong: folded is
    :data:`~dodo.constants.BETA_FOLDED` (100) and disordered is
    :data:`~dodo.constants.BETA_DISORDERED` (0), so colouring by B-factor lights up the
    folded core. v1's docstrings claimed the reverse in three places; its code did this.

    Rebuildable loops inside a folded domain are annotated as disordered. They are
    rebuilt geometry, not folded core, and the point of the annotation is to show which
    coordinates DODO generated.
    """
    if not annotate_regions:
        return structure.b_factor

    domains = structure.domains
    if not domains:
        raise InvalidRegionError(
            "annotate_regions=True but this structure has no domains assigned, so "
            "there is nothing to annotate. Identify regions first."
        )
    values = np.full(structure.n_residues, np.nan, dtype=np.float64)
    for domain in domains:
        values[domain.span.slice] = (
            BETA_FOLDED if domain.kind is DomainKind.FOLDED else BETA_DISORDERED
        )
        for loop in domain.loops:
            values[loop.slice] = BETA_DISORDERED
    unassigned = np.flatnonzero(~np.isfinite(values))
    if unassigned.size:
        values[unassigned] = BETA_DISORDERED
        _note(
            structure,
            f"annotate_regions: {unassigned.size} residue(s) belong to no domain "
            f"(first is {structure.residue_label(int(unassigned[0]))}); they were "
            f"written as disordered.",
        )
    return values


def _elements(structure: Structure) -> np.ndarray:
    """Return the element symbol per atom, inferring any that are missing."""
    element = structure.element
    missing = np.flatnonzero(element == "")
    if missing.size == 0:
        return element
    filled: np.ndarray = element.copy()
    for index in missing:
        name = str(structure.atom_name[index]).strip()
        guess = next((char for char in name if char.isalpha()), "")
        filled[index] = guess.upper()
    _note(
        structure,
        f"{missing.size} atom(s) had no element symbol; it was inferred from the atom "
        f"name so viewers do not have to guess. Check the input file.",
    )
    return filled


def _selected_atom_indices(structure: Structure, *, ca_only: bool) -> np.ndarray:
    """Return the indices of the atoms to write, in file order.

    ``ca_only`` means *exactly one alpha carbon per residue*, and specifically the same
    atom :attr:`~dodo.structure.Structure.ca_indices` resolves to, so the file agrees
    with every geometric calculation DODO runs on the same structure. Two things follow,
    and both are recorded rather than passed over:

    * A residue with no CA -- a ligand, a water, an incomplete residue -- is omitted.
    * A residue carrying *more than one* CA contributes only one. This is reachable
      through DODO's own reader, which keeps a blank-altLoc CA alongside an altLoc-'A'
      CA. Writing both would emit an atom that no CONECT record can reach, i.e. exactly
      the disconnected dot this module exists to prevent.

    The accounting used to compute ``n_residues - n_selected`` and call the result
    "residues with no CA", which reported "omitted -1 residue(s)" for the duplicate case
    and never mentioned the duplicate at all.
    """
    if not ca_only:
        return np.arange(structure.n_atoms, dtype=np.int64)
    ca_atoms = np.flatnonzero(structure.atom_name == "CA").astype(np.int64)
    if ca_atoms.size == 0:
        raise EmptyStructureError(
            f"ca_only=True but this structure has no CA atoms (source: {structure.source!r})."
        )
    residue_of_ca = structure.residue_index[ca_atoms]
    # Last CA wins, which is what Structure.ca_indices does with the same assignment.
    chosen = np.full(structure.n_residues, -1, dtype=np.int64)
    chosen[residue_of_ca] = ca_atoms
    selected: np.ndarray = chosen[chosen >= 0]

    per_residue = np.bincount(residue_of_ca, minlength=structure.n_residues)
    missing = np.flatnonzero(per_residue == 0)
    if missing.size:
        _note(
            structure,
            f"ca_only=True omitted {missing.size} residue(s) that have no CA atom "
            f"(ligands, waters or incomplete residues); the first is "
            f"{structure.residue_label(int(missing[0]))}.",
        )
    crowded = np.flatnonzero(per_residue > 1)
    if crowded.size:
        extra = int(per_residue[crowded].sum()) - int(crowded.size)
        _note(
            structure,
            f"ca_only=True writes exactly one CA per residue, the one "
            f"Structure.ca_indices resolves to: {crowded.size} residue(s) carry more "
            f"than one CA atom (the first is "
            f"{structure.residue_label(int(crowded[0]))}), so {extra} CA atom(s) were "
            f"not written. This usually means unresolved alternate conformers.",
        )
    return selected


def _check_coordinates(structure: Structure, atom_indices: np.ndarray) -> None:
    """Reject coordinates the ``%8.3f`` fields cannot represent.

    Both failure modes here produce a file that loads without complaint and is wrong,
    so both raise. Non-finite coordinates in particular were a real v1 output: its
    samplers returned NaN on failure and the writer wrote it.
    """
    coords = structure.xyz[atom_indices]
    bad = np.flatnonzero(~np.isfinite(coords).all(axis=1))
    if bad.size:
        first = int(atom_indices[bad[0]])
        raise GeometryError(
            f"{bad.size} atom(s) have non-finite coordinates and cannot be written; "
            f"the first is atom {first} of residue "
            f"{structure.residue_label(int(structure.residue_index[first]))}."
        )
    out_of_range = np.flatnonzero(((coords < _COORD_MIN) | (coords > _COORD_MAX)).any(axis=1))
    if out_of_range.size:
        first = int(atom_indices[out_of_range[0]])
        raise StructureFileError(
            f"{out_of_range.size} atom(s) lie outside the range the PDB coordinate "
            f"fields can represent ({_COORD_MIN} to {_COORD_MAX} A); the first is atom "
            f"{first} of residue "
            f"{structure.residue_label(int(structure.residue_index[first]))}. "
            f"Translate the structure toward the origin, or write mmCIF."
        )


def _ter_line(serial: int, residue_field: str) -> str:
    """Render a TER record.

    TER consumes an atom serial number, per the specification, which is why CONECT
    reads serials out of :attr:`_Layout.serials` instead of recomputing them.
    """
    ter = f"TER   {_hybrid36(serial, _SERIAL_WIDTH, MAX_PDB_ATOM_SERIAL)}      {residue_field}"
    # The insertion-code column is the last one and is usually blank; do not leave it as
    # trailing whitespace.
    return ter.rstrip()


@dataclass(frozen=True, slots=True)
class _Layout:
    """The written atom block, plus what CONECT and SEQRES generation needs from it."""

    #: ATOM/HETATM/TER lines, in output order.
    atom_lines: list[str]
    #: Indices into the structure's atom arrays, one per written atom.
    atom_indices: np.ndarray
    #: PDB serial number actually written for each entry of ``atom_indices``. Not simply
    #: ``i + 1``: TER records consume a serial, so CONECT must read them from here.
    serials: np.ndarray
    #: Single-character chain id per chain index.
    chain_ids: list[str]
    #: Boolean mask over residues, ``True`` where the residue contributed at least one
    #: atom to ``atom_lines``. Everything derived from the atom block -- SEQRES, the
    #: range checks, the substitution accounting -- is restricted to these, so a residue
    #: that is not in the file cannot fail the write or be reported as written.
    written_residues: np.ndarray


def _layout(structure: Structure, *, ca_only: bool, annotate_regions: bool) -> _Layout:
    """Build the ATOM/TER block for one structure."""
    atom_indices = _selected_atom_indices(structure, ca_only=ca_only)
    _check_coordinates(structure, atom_indices)

    chain_ids = _chain_ids(structure)
    b_factor = _residue_b_factors(structure, annotate_regions=annotate_regions)
    element = _elements(structure)

    xyz = structure.xyz
    atom_name = structure.atom_name
    residue_index = structure.residue_index
    residue_name = structure.residue_name
    insertion_code = structure.insertion_code
    occupancy = structure.occupancy
    chain_index = structure.chain_index

    # Only the residues that actually reach the file are formatted or range-checked.
    # Doing this for every residue meant a ligand with an unrepresentable residue number
    # aborted a perfectly valid CA-trace write, and occupancy substitutions were counted
    # and reported for residues that never appeared in the output.
    written_residues = np.zeros(structure.n_residues, dtype=bool)
    written_residues[residue_index[atom_indices]] = True

    # Residue-level fields are formatted once per residue rather than once per atom;
    # on a 61,511-atom EM assembly that is the difference between comfortable and not.
    residue_fields: list[str] = []
    for i in range(structure.n_residues):
        if not written_residues[i]:
            residue_fields.append("")
            continue
        name = str(residue_name[i])
        chain_char = chain_ids[int(chain_index[i])]
        resseq = _hybrid36(int(structure.residue_number[i]), _RESSEQ_WIDTH, MAX_PDB_RESIDUE_NUMBER)
        icode = (str(insertion_code[i]) or " ")[:1]
        residue_fields.append(f"{name:>3} {chain_char}{resseq}{icode}")
    record_names = [
        "ATOM  " if str(name) in THREE_TO_ONE else "HETATM" for name in structure.residue_name
    ]

    occupancy_fields: list[str] = []
    b_factor_fields: list[str] = []
    substitutions = 0
    for i in range(structure.n_residues):
        if not written_residues[i]:
            occupancy_fields.append("")
            b_factor_fields.append("")
            continue
        occ_field, occ_sub = _format_real(float(occupancy[i]), _DEFAULT_OCCUPANCY)
        b_field, b_sub = _format_real(float(b_factor[i]), _DEFAULT_B_FACTOR)
        occupancy_fields.append(occ_field)
        b_factor_fields.append(b_field)
        substitutions += int(occ_sub) + int(b_sub)
    if substitutions:
        _note(
            structure,
            f"{substitutions} occupancy/B-factor value(s) were missing or outside the "
            f"range the PDB {_REAL_MIN} to {_REAL_MAX} field can hold; defaults "
            f"({_DEFAULT_OCCUPANCY:.2f} / {_DEFAULT_B_FACTOR:.2f}) or clamped values "
            f"were written.",
        )

    lines: list[str] = []
    serials = np.zeros(atom_indices.size, dtype=np.int64)
    serial = 1
    previous_chain = -1
    previous_residue = -1
    for position, atom in enumerate(atom_indices):
        atom = int(atom)
        residue = int(residue_index[atom])
        chain = int(chain_index[residue])
        if previous_chain != -1 and chain != previous_chain:
            lines.append(_ter_line(serial, residue_fields[previous_residue]))
            serial += 1
        x, y, z = xyz[atom]
        lines.append(
            f"{record_names[residue]}"
            f"{_hybrid36(serial, _SERIAL_WIDTH, MAX_PDB_ATOM_SERIAL)} "
            f"{_format_atom_name(str(atom_name[atom]), str(element[atom]))} "
            f"{residue_fields[residue]}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
            f"{occupancy_fields[residue]}{b_factor_fields[residue]}"
            f"          {str(element[atom])[:2]:>2}"
        )
        serials[position] = serial
        serial += 1
        previous_chain = chain
        previous_residue = residue
    if previous_residue != -1:
        lines.append(_ter_line(serial, residue_fields[previous_residue]))

    return _Layout(
        atom_lines=lines,
        atom_indices=atom_indices,
        serials=serials,
        chain_ids=chain_ids,
        written_residues=written_residues,
    )


# ---------------------------------------------------------------------------
# CONECT
# ---------------------------------------------------------------------------


def _backbone_positions(structure: Structure, atom_indices: np.ndarray) -> dict[str, np.ndarray]:
    """Locate each residue's backbone atoms as positions within the written atoms.

    ``atom_indices`` are the structure atoms actually written, in output order (a
    :class:`_Layout`'s ``atom_indices``); a returned position indexes into *that* array.
    Taking the indices directly rather than a whole ``_Layout`` is what lets the mmCIF
    writer share this connectivity logic without first building a PDB atom block.

    ``-1`` means the residue does not have that atom, which is the normal case for a
    CA-only structure and for the incomplete terminal residues real files contain.

    A residue can also carry the *same* backbone name twice -- unresolved alternate
    conformers, which DODO's own reader preserves. Only one of them can hold the
    residue's backbone connectivity, so the others are written with no connectivity
    record at all and render as isolated dots. That is recorded, because an unbonded atom
    in the output is precisely the failure this module exists to prevent.
    """
    names = structure.atom_name[atom_indices]
    residues = structure.residue_index[atom_indices]
    positions: dict[str, np.ndarray] = {}
    repeated = 0
    first_repeat = -1
    for name in BACKBONE_ATOMS:
        found = np.full(structure.n_residues, -1, dtype=np.int64)
        selected = np.flatnonzero(names == name)
        found[residues[selected]] = selected
        positions[name] = found
        counts = np.bincount(residues[selected], minlength=structure.n_residues)
        crowded = np.flatnonzero(counts > 1)
        if crowded.size:
            repeated += int(counts[crowded].sum()) - int(crowded.size)
            first = int(crowded[0])
            first_repeat = first if first_repeat < 0 else min(first_repeat, first)
    if repeated:
        _note(
            structure,
            f"{repeated} atom(s) repeat a backbone atom name within their own residue "
            f"(the first such residue is {structure.residue_label(first_repeat)}). Only "
            f"one atom of each backbone name per residue can carry the backbone "
            f"connectivity, so the repeats were written with no CONECT record and will "
            f"render as isolated dots. Resolve the alternate conformers, or write with "
            f"ca_only=True.",
        )
    return positions


def _bonds(structure: Structure, atom_indices: np.ndarray) -> list[tuple[int, int]]:
    """List bonds as pairs of positions within the written atoms.

    ``atom_indices`` are the structure atoms actually written, in output order; each
    returned pair indexes into that array. The mmCIF writer reuses this so its
    connectivity records describe exactly the bonds the PDB writer's CONECT records do.

    Intra-residue N-CA, CA-C and C-O wherever both atoms are present, plus one bond per
    residue junction. Whether a junction is bonded is decided from the CA-CA separation
    where both alpha carbons exist, and from the C-N separation otherwise: the CA-CA
    virtual bond is the one distance every DODO structure has and the one the geometry
    engines constrain, so it is the more reliable test. What gets *written* is the
    peptide C-N bond when both atoms are present, and the CA-CA virtual bond when they
    are not -- the CA-only case, which is DODO's main output and the one where a missing
    CONECT record turns a rebuilt IDR into disconnected dots.

    A CA-CA record is deliberately *not* added on top of a real peptide bond: it would
    draw a spurious 3.8 A bond across every peptide unit in an all-atom model.
    """
    positions = _backbone_positions(structure, atom_indices)
    xyz = structure.xyz[atom_indices]
    bonds: list[tuple[int, int]] = []

    for first, second in (("N", "CA"), ("CA", "C"), ("C", "O")):
        left, right = positions[first], positions[second]
        both = np.flatnonzero((left >= 0) & (right >= 0))
        bonds.extend((int(left[i]), int(right[i])) for i in both)

    chain_index = structure.chain_index
    ca, carbon, nitrogen = positions["CA"], positions["C"], positions["N"]
    for residue in range(structure.n_residues - 1):
        nxt = residue + 1
        if chain_index[residue] != chain_index[nxt]:
            continue  # never bond across a chain break
        has_peptide = carbon[residue] >= 0 and nitrogen[nxt] >= 0
        if ca[residue] >= 0 and ca[nxt] >= 0:
            separation = float(np.linalg.norm(xyz[ca[residue]] - xyz[ca[nxt]]))
            bonded = separation <= _MAX_BONDED_CA_CA_DISTANCE
        elif has_peptide:
            separation = float(np.linalg.norm(xyz[carbon[residue]] - xyz[nitrogen[nxt]]))
            bonded = separation <= _MAX_BONDED_C_N_DISTANCE
        else:
            continue
        if not bonded:
            continue
        if has_peptide:
            bonds.append((int(carbon[residue]), int(nitrogen[nxt])))
        else:
            bonds.append((int(ca[residue]), int(ca[nxt])))
    return bonds


def _conect_lines(structure: Structure, layout: _Layout) -> list[str]:
    """Render CONECT records for one structure's written atoms.

    Each bond is listed once, from the lower serial, up to
    :data:`_CONECT_PARTNERS_PER_LINE` partners per record with continuation records
    beyond that. Every viewer DODO targets accepts single-direction CONECT.
    """
    serials = layout.serials
    partners: dict[int, list[int]] = {}
    for left, right in _bonds(structure, layout.atom_indices):
        low, high = sorted((int(serials[left]), int(serials[right])))
        partners.setdefault(low, []).append(high)

    lines: list[str] = []
    for serial in sorted(partners):
        bonded = sorted(partners[serial])
        origin = _hybrid36(serial, _SERIAL_WIDTH, MAX_PDB_ATOM_SERIAL)
        for start in range(0, len(bonded), _CONECT_PARTNERS_PER_LINE):
            chunk = bonded[start : start + _CONECT_PARTNERS_PER_LINE]
            fields = "".join(_hybrid36(p, _SERIAL_WIDTH, MAX_PDB_ATOM_SERIAL) for p in chunk)
            lines.append(f"CONECT{origin}{fields}")
    return lines


# ---------------------------------------------------------------------------
# Header records
# ---------------------------------------------------------------------------


def _cryst1_line(box: tuple[float, float, float]) -> str:
    """Render a CRYST1 record for an orthorhombic P 1 cell.

    DODO does not do periodic boundaries; this exists because some viewers want the
    record present, which is why the default is to write no CRYST1 at all. Callers who
    want the record but have no cell in mind pass
    :data:`~dodo.constants.DEFAULT_BOX_DIMENSIONS`.
    """
    a, b, c = box
    for length in (a, b, c):
        if not np.isfinite(length) or length <= 0:
            raise GeometryError(f"CRYST1 box dimensions must be positive and finite, got {box!r}.")
    return f"CRYST1{a:9.3f}{b:9.3f}{c:9.3f}{90.0:7.2f}{90.0:7.2f}{90.0:7.2f} {'P 1':<11}{1:>4}"


def _seqres_lines(structure: Structure, layout: _Layout) -> list[str]:
    """Render SEQRES records, :data:`~dodo.constants.RESIDUES_PER_SEQRES_LINE` per line.

    These describe the residues actually written to this file, not
    :attr:`~dodo.structure.Chain.full_sequence` and not the residues a ``ca_only`` pass
    left out. A SEQRES block listing residues the ATOM records do not contain would
    disagree with the coordinates in its own file, and readers that trust SEQRES would
    then report residues no reader of the coordinates can find.

    Grouping is by *written chain id*, not by run of ``chain_index``. Both DODO readers
    deliberately preserve interleaved records as separate chains, so one chain id can
    cover several runs of residues; a block per run emitted two ``serNum``-1 blocks with
    the same chainID, each declaring only its own run's ``numRes``, which is the same
    contradiction stated a different way. The PDB format cannot express one chain id
    twice, so the runs are concatenated in the order they are written and the merge is
    recorded.
    """
    chain_index = structure.chain_index
    order: list[str] = []
    names_by_chain: dict[str, list[str]] = {}
    runs_by_chain: dict[str, int] = {}
    previous: str | None = None
    for residue in range(structure.n_residues):
        if not layout.written_residues[residue]:
            continue
        chain = layout.chain_ids[int(chain_index[residue])]
        if chain not in names_by_chain:
            order.append(chain)
            names_by_chain[chain] = []
            runs_by_chain[chain] = 0
        if chain != previous:
            runs_by_chain[chain] += 1
        previous = chain
        names_by_chain[chain].append(str(structure.residue_name[residue]))

    interleaved = [chain for chain in order if runs_by_chain[chain] > 1]
    if interleaved:
        _note(
            structure,
            f"SEQRES: chain id(s) {', '.join(interleaved)} cover more than one run of "
            f"residues, which the single SEQRES chainID column cannot express. Each "
            f"chain's residues were listed in one block, in the order they are written. "
            f"The ATOM records are unchanged.",
        )

    lines: list[str] = []
    for chain in order:
        names = names_by_chain[chain]
        for serial, offset in enumerate(range(0, len(names), RESIDUES_PER_SEQRES_LINE), start=1):
            chunk = names[offset : offset + RESIDUES_PER_SEQRES_LINE]
            residues = " ".join(f"{name:>3}" for name in chunk)
            lines.append(f"SEQRES {serial:>3} {chain} {len(names):>4}  {residues}")
    return lines


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def structure_to_pdb_lines(
    structure: Structure,
    *,
    conect: bool = True,
    annotate_regions: bool = False,
    ca_only: bool = False,
    box: tuple[float, float, float] | None = None,
    seqres: bool = False,
    model_num: int | None = None,
) -> list[str]:
    """Render one structure as PDB records.

    Parameters
    ----------
    structure
        The structure to render.
    conect
        Emit CONECT records for backbone connectivity. On by default, and it matters:
        the CA-CA virtual bond is 3.80 A, past the automatic bond-detection cutoff in
        VMD and PyMOL, so a CA-only model without CONECT renders as disconnected dots.
    annotate_regions
        Replace the B-factor column with the region annotation: folded residues get
        :data:`~dodo.constants.BETA_FOLDED` (100) and disordered residues
        :data:`~dodo.constants.BETA_DISORDERED` (0), so colouring by B-factor
        highlights the folded core. Rebuildable loops count as disordered.
    ca_only
        Write exactly one alpha carbon per residue and nothing else -- the same atom
        :attr:`~dodo.structure.Structure.ca_indices` resolves to, so the file agrees with
        every geometric calculation DODO runs on the structure. Residues with no CA are
        omitted; residues carrying more than one contribute one. Both are recorded in
        ``structure.notes``.
    box
        Cell dimensions for a CRYST1 record, in Angstroms. ``None`` means no CRYST1.
        Pass :data:`~dodo.constants.DEFAULT_BOX_DIMENSIONS` for the generic box.
    seqres
        Emit SEQRES records for the residues written to this file.
    model_num
        When given, wrap the atom records in ``MODEL``/``ENDMDL``.

    Returns
    -------
    list[str]
        Records without trailing newlines. No ``END``: that terminates a *file*, and
        this function renders one structure, which may be one frame of many.

    Raises
    ------
    GeometryError
        If any coordinate is non-finite, or ``box`` is not positive and finite.
    StructureFileError
        If a coordinate, atom serial or residue number cannot be represented in the
        PDB format even with hybrid-36 encoding.
    InvalidRegionError
        If ``annotate_regions`` is requested for a structure with no domains.

    Notes
    -----
    Anomalies that do not prevent writing -- a substituted occupancy, a remapped chain
    id, a residue dropped by ``ca_only`` -- are appended to ``structure.notes`` rather
    than passed over in silence.

    Record order follows the specification: SEQRES belongs to the primary-structure
    section and CRYST1 to the crystallographic section that follows it, so SEQRES comes
    first. The same reasoning puts CONECT after ``ENDMDL``.
    """
    layout = _layout(structure, ca_only=ca_only, annotate_regions=annotate_regions)
    # Built before anything is emitted so an invalid box still fails early, emitted after
    # SEQRES so the records land in specification order.
    cryst1 = None if box is None else _cryst1_line(box)
    lines: list[str] = []
    if seqres:
        lines.extend(_seqres_lines(structure, layout))
    if cryst1 is not None:
        lines.append(cryst1)
    if model_num is not None:
        lines.append(f"MODEL     {model_num:>4}")
    lines.extend(layout.atom_lines)
    if model_num is not None:
        lines.append("ENDMDL")
    if conect:
        lines.extend(_conect_lines(structure, layout))
    return lines


def write_pdb(
    structure: Structure | Sequence[Structure],
    path: str | Path,
    *,
    conect: bool = True,
    models_as_frames: bool = True,
    annotate_regions: bool = False,
    ca_only: bool = False,
    box: tuple[float, float, float] | None = None,
    seqres: bool = False,
) -> None:
    """Write one or more structures to a PDB file.

    Parameters
    ----------
    structure
        A single :class:`~dodo.structure.Structure`, or a sequence of them.
    path
        Destination file. A bare relative name such as ``"out.pdb"`` works; v1's
        writer rejected it, because ``os.path.dirname("out.pdb")`` is ``""`` and it
        passed that straight to a directory-exists check.
    conect
        Emit CONECT records. See :func:`structure_to_pdb_lines`.
    models_as_frames
        With a sequence of structures: ``True`` writes them as successive
        ``MODEL``/``ENDMDL`` frames of one file, which is what makes the
        pseudo-trajectory output work in a viewer, and requires every frame to have
        identical atom count and ordering. ``False`` writes one file per structure,
        named ``<stem>_1<suffix>``, ``<stem>_2<suffix>``, and so on, with no constraint
        between them.

        The behaviour depends on whether a sequence was passed, not on how long it is:
        a one-element sequence is still an ensemble, so it still gets a ``MODEL`` record
        or the ``_1`` suffix. A script that writes N conformers should not have its
        output silently renamed when N happens to be 1. A bare
        :class:`~dodo.structure.Structure` is written plainly, with no ``MODEL``
        records, whatever this flag says.
    annotate_regions, ca_only, box, seqres
        See :func:`structure_to_pdb_lines`.

    Raises
    ------
    EmptyStructureError
        If an empty sequence of structures is given.
    StructureFileError
        If ``path`` is not writable as a file -- its directory does not exist, or it is
        itself an existing directory -- if frames disagree on atom count, ordering or
        chain division while ``models_as_frames`` is set, or if some value cannot be
        represented in the PDB format, ASCII included.

    Notes
    -----
    In a multi-frame file the CONECT records are written once, after the last
    ``ENDMDL``, as the specification requires. That is sound precisely because the
    frames are required to share a topology -- atom order *and* chain division, since
    TER records consume atom serials -- so one set of serial numbers describes them all.
    """
    one_structure = isinstance(structure, Structure)
    frames = [structure] if isinstance(structure, Structure) else list(structure)
    if not frames:
        raise EmptyStructureError("No structures given to write.")
    target = Path(path)

    if one_structure:
        _write_lines(
            target,
            [
                *structure_to_pdb_lines(
                    frames[0],
                    conect=conect,
                    annotate_regions=annotate_regions,
                    ca_only=ca_only,
                    box=box,
                    seqres=seqres,
                ),
                "END",
            ],
        )
        return

    if not models_as_frames:
        for index, frame in enumerate(frames, start=1):
            _write_lines(
                target.with_name(f"{target.stem}_{index}{target.suffix}"),
                [
                    *structure_to_pdb_lines(
                        frame,
                        conect=conect,
                        annotate_regions=annotate_regions,
                        ca_only=ca_only,
                        box=box,
                        seqres=seqres,
                    ),
                    "END",
                ],
            )
        return

    _require_matching_topology(frames)
    cryst1 = None if box is None else _cryst1_line(box)
    header: list[str] = []
    body: list[str] = []
    conect_lines: list[str] = []
    for index, frame in enumerate(frames, start=1):
        layout = _layout(frame, ca_only=ca_only, annotate_regions=annotate_regions)
        if index == 1:
            if seqres:
                header.extend(_seqres_lines(frame, layout))
            if conect:
                conect_lines = _conect_lines(frame, layout)
        body.append(f"MODEL     {index:>4}")
        body.extend(layout.atom_lines)
        body.append("ENDMDL")
    if cryst1 is not None:
        header.append(cryst1)
    _write_lines(target, [*header, *body, *conect_lines, "END"])


def _require_matching_topology(frames: Sequence[Structure]) -> None:
    """Reject frames that cannot share one set of atom records.

    A MODEL/ENDMDL series is one topology sampled repeatedly. Frames that disagree on
    atom count or ordering would produce a file whose CONECT records and residue
    identities apply to the first frame only -- readable, and wrong for every frame
    after it.

    ``chain_index`` is compared for a reason that is easy to miss: it is what decides
    where the TER records go, and a TER record consumes an atom serial. Two frames that
    agree on every other field but split their chains differently are written with
    *different serial numbering*, so frame 1's CONECT records -- the only ones in the
    file -- bond the wrong atoms in frame 2, including bonding to a serial that is a TER
    record, while a real atom is left with no bond at all. MDTraj loads such a file
    without complaint, which makes it the worst possible output.
    """
    reference = frames[0]
    for index, frame in enumerate(frames[1:], start=2):
        if frame.n_atoms != reference.n_atoms:
            raise StructureFileError(
                f"Cannot write these structures as frames of one file: frame {index} has "
                f"{frame.n_atoms} atoms but frame 1 has {reference.n_atoms}. Write them "
                f"separately with models_as_frames=False."
            )
        for label, left, right in _topology_fields(reference, frame):
            if not np.array_equal(left, right):
                raise StructureFileError(
                    f"Cannot write these structures as frames of one file: frame "
                    f"{index} disagrees with frame 1 on {label}. Write them separately "
                    f"with models_as_frames=False."
                )
        if [chain.chain_id for chain in frame.chains] != [
            chain.chain_id for chain in reference.chains
        ]:
            raise StructureFileError(
                f"Cannot write these structures as frames of one file: frame {index} "
                f"disagrees with frame 1 on chain ids."
            )


def _topology_fields(
    reference: Structure, frame: Structure
) -> tuple[tuple[str, np.ndarray, np.ndarray], ...]:
    """List the per-atom record arrays two frames of one file must agree on, with labels."""
    return (
        ("atom names", reference.atom_name, frame.atom_name),
        ("residue names", reference.residue_name, frame.residue_name),
        ("residue numbers", reference.residue_number, frame.residue_number),
        ("insertion codes", reference.insertion_code, frame.insertion_code),
        ("residue indices", reference.residue_index, frame.residue_index),
        ("chain indices", reference.chain_index, frame.chain_index),
    )


def matching_topology(first: Structure, second: Structure) -> bool:
    """Return True if the two structures can be written as frames of one multi-model file.

    Exactly the comparison :func:`_require_matching_topology` enforces at write time -- one
    shared implementation (:func:`_topology_fields`), so a caller deciding WHAT to write
    cannot drift from what the writers accept.
    """
    if first.n_atoms != second.n_atoms:
        return False
    if any(not np.array_equal(left, right) for _, left, right in _topology_fields(first, second)):
        return False
    return [chain.chain_id for chain in first.chains] == [chain.chain_id for chain in second.chains]


def _write_lines(path: Path, lines: list[str]) -> None:
    """Write records to ``path``, one per line.

    The parent directory is checked but never created: silently creating directories
    turns a mistyped path into a file the user cannot find.

    Every way this can fail raises :class:`~dodo.exceptions.StructureFileError`. That is
    the module's contract and it used to have two holes: an existing directory passed as
    ``path`` surfaced as ``IsADirectoryError``, because only the *parent* was checked,
    and a field holding a non-ASCII character surfaced as ``UnicodeEncodeError`` from
    inside :mod:`pathlib` -- while every other "cannot be represented in PDB" case
    raised.
    """
    parent = path.parent
    # Path("out.pdb").parent is Path("."), which exists. The v1 writer used
    # os.path.dirname, got "", and raised on every bare relative filename.
    if not parent.is_dir():
        raise StructureFileError(
            f"Cannot write {str(path)!r}: the directory {str(parent)!r} does not exist."
        )
    if path.is_dir():
        raise StructureFileError(
            f"Cannot write {str(path)!r}: it is an existing directory, not a file."
        )
    text = "\n".join(lines) + "\n"
    try:
        text.encode("ascii")
    except UnicodeEncodeError as error:
        offending = next((line for line in lines if not line.isascii()), "")
        raise StructureFileError(
            f"Cannot write {str(path)!r}: the PDB format is ASCII and this structure "
            f"contains the character {text[error.start : error.end]!r}, which has no "
            f"ASCII representation. The record is {offending!r}."
        ) from error
    try:
        path.write_text(text, encoding="ascii")
    except OSError as error:
        raise StructureFileError(f"Cannot write {str(path)!r}: {error}.") from error
