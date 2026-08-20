"""mmCIF (PDBx) writing.

Why this module exists
----------------------
DODO's only structure writer used to be the PDB writer in :mod:`dodo.io.write`. Asking
for ``.cif`` output ran that writer and put PDB records -- ``MODEL``/``ATOM``/``ENDMDL``
-- into a file merely *named* ``.cif``. A viewer that trusts the extension then hands
those bytes to its mmCIF parser, which finds no ``_atom_site`` loop and refuses the file.
That is exactly what "the multi-model .cif does not open in VMD, but the .pdb does" was:
the ``.cif`` was never mmCIF at all.

This module writes real mmCIF, and in particular expresses multiple conformers the way
mmCIF does: one ``_atom_site`` loop whose rows carry a ``pdbx_PDB_model_num`` column,
*not* PDB's ``MODEL``/``ENDMDL`` bracketing. That single difference is the whole point of
writing mmCIF for an ensemble.

What it reuses from the PDB writer, and why
-------------------------------------------
Atom selection (``ca_only``), region annotation, element inference and -- most
importantly -- the backbone/CA-CA bond list are format-neutral facts about a structure,
so they are computed by the same functions the PDB writer uses (:mod:`dodo.io.write`).
Sharing :func:`dodo.io.write._bonds` in particular means the connectivity this module
emits (as ``_struct_conn``) describes exactly the bonds the PDB writer emits as
``CONECT`` -- there is one definition of what is bonded, not two that can drift.

A note on connectivity and viewers
-----------------------------------
The CA-CA virtual bond is 3.81 A, past the automatic distance-based bond cutoff every
common viewer uses, so a CA-only trace needs its connectivity stated explicitly or it
renders as disconnected dots. In PDB that is ``CONECT``; the mmCIF equivalent is
``_struct_conn``. PyMOL, ChimeraX and the wwPDB tools read ``_struct_conn`` and draw the
bonds. Some builds of VMD's mmCIF plugin read coordinates only and re-derive bonds by
distance regardless, so a CA-only mmCIF may still show gaps there even though the file is
correct; the honest fix for a CA-only model in that specific viewer is the PDB output,
which is why ``_struct_conn`` is written but not claimed to be universal.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Final

import numpy as np

from ..constants import THREE_TO_ONE
from ..exceptions import EmptyStructureError, GeometryError, StructureFileError
from ..structure import Structure
from .write import (
    _bonds,
    _elements,
    _note,
    _require_matching_topology,
    _residue_b_factors,
    _selected_atom_indices,
)

__all__ = [
    "structure_to_cif_text",
    "write_cif",
]

#: The ``data_`` block name. mmCIF requires exactly one block header per file.
_BLOCK_NAME: Final = "dodo"

#: Default occupancy / B-factor written when the stored value is not finite, matching the
#: PDB writer's substitutions so the two formats say the same thing about missing data.
_DEFAULT_OCCUPANCY: Final = 1.0
_DEFAULT_B_FACTOR: Final = 0.0

#: The ``.``/``?`` placeholders. ``.`` means "inapplicable", ``?`` means "unknown"; mmCIF
#: uses ``.`` for a residue that carries no alternate conformer and ``?`` for an absent
#: insertion code, and that is the convention RCSB writes, so it is the one round-tripped.
_INAPPLICABLE: Final = "."
_UNKNOWN: Final = "?"

#: ``_atom_site`` columns, in the order RCSB and AlphaFold write them. Both the ``label_``
#: and ``auth_`` identifiers are emitted with the same values: DODO's own reader and every
#: common viewer prefer ``auth_`` for display, while ``label_`` is what the mmCIF
#: dictionary marks mandatory, so writing both leaves nothing to fall back on or infer.
_ATOM_SITE_COLUMNS: Final = (
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "pdbx_PDB_ins_code",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "occupancy",
    "B_iso_or_equiv",
    "auth_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_atom_id",
    "pdbx_PDB_model_num",
)

#: ``_struct_conn`` columns. Bonds are stated by label identifier (asym/comp/seq/atom plus
#: insertion code) rather than by atom serial, because a serial is per-model and this
#: connectivity describes every model at once -- exactly why it is written a single time.
_STRUCT_CONN_COLUMNS: Final = (
    "id",
    "conn_type_id",
    "ptnr1_label_asym_id",
    "ptnr1_label_comp_id",
    "ptnr1_label_seq_id",
    "ptnr1_label_atom_id",
    "pdbx_ptnr1_PDB_ins_code",
    "ptnr2_label_asym_id",
    "ptnr2_label_comp_id",
    "ptnr2_label_seq_id",
    "ptnr2_label_atom_id",
    "pdbx_ptnr2_PDB_ins_code",
)


# ---------------------------------------------------------------------------
# STAR value encoding
# ---------------------------------------------------------------------------

_WHITESPACE: Final = frozenset(" \t\n\r")
#: Characters that begin a STAR construct, so a value starting with one must be quoted.
_LEADING_SPECIALS: Final = frozenset("_#$'\";[]")
_RESERVED_PREFIXES: Final = ("data_", "save_")
_RESERVED_WORDS: Final = frozenset({"loop_", "stop_", "global_"})


def _needs_quoting(value: str) -> bool:
    """Whether a STAR value has to be quoted to survive a round trip through the reader.

    CIF 1.1 permits a bare value to contain an internal quote character (``O5'`` is
    legal, and DODO's own reader reads it), but strict parsers vary and RCSB quotes such
    a value anyway, so this errs the safe way and quotes anything holding a ``'`` or
    ``"``. A value is otherwise bare only if it contains no whitespace, does not begin
    with a character that starts another STAR construct, and is not a reserved word.
    """
    if not value:
        return True
    if any(ch in _WHITESPACE for ch in value):
        return True
    if "'" in value or '"' in value:
        return True
    if value[0] in _LEADING_SPECIALS:
        return True
    lowered = value.lower()
    return lowered in _RESERVED_WORDS or lowered.startswith(_RESERVED_PREFIXES)


def _encode_value(value: str) -> str:
    """Encode one STAR value: bare when it can be, quoted when it must be.

    Chooses the quote character the value does not itself contain, which is how an atom
    name like ``O5'`` becomes ``"O5'"``. A value carrying *both* quote characters cannot
    be single-line quoted at all; rather than emit something a reader would mis-split,
    this raises -- consistent with the PDB writer, which also refuses to write a field it
    cannot represent instead of writing a corrupt one. No atom, residue, chain or element
    name in a structure file contains both, so this is a guard, not a real limitation.
    """
    if not _needs_quoting(value):
        return value
    if "'" not in value:
        return f"'{value}'"
    if '"' not in value:
        return f'"{value}"'
    raise StructureFileError(
        f"Cannot write the mmCIF value {value!r}: it contains both a single and a double "
        f"quote, which STAR cannot represent on one line. Refusing to write a file a "
        f"reader would mis-split."
    )


def _loop(columns: Sequence[str], category: str, rows: Sequence[Sequence[str]]) -> list[str]:
    """Render a ``loop_`` for ``category`` with ``columns``, padding cells for legibility.

    Padding to a common width per column is what RCSB does and is purely cosmetic: the
    reader splits on whitespace, so any number of spaces between values is equivalent. It
    is worth the one extra pass because a column-aligned file is one a human can debug.
    """
    lines = ["loop_", *(f"_{category}.{name}" for name in columns)]
    if not rows:
        return lines
    widths = [max(len(row[i]) for row in rows) for i in range(len(columns))]
    for row in rows:
        lines.append(" ".join(cell.ljust(widths[i]) for i, cell in enumerate(row)).rstrip())
    return lines


# ---------------------------------------------------------------------------
# Per-model preparation
# ---------------------------------------------------------------------------


def _reject_non_finite(structure: Structure, atom_indices: np.ndarray) -> None:
    """Reject non-finite coordinates before they reach the file.

    mmCIF has no coordinate range limit -- that freedom is the reason DODO prefers it for
    large structures -- so unlike the PDB writer there is no over-range check here. A
    ``NaN`` still must not be written: DODO's samplers return ``NaN`` on failure, and a
    file full of ``NaN`` positions loads as a heap of atoms at nowhere.
    """
    coords = structure.xyz[atom_indices]
    bad = np.flatnonzero(~np.isfinite(coords).all(axis=1))
    if bad.size:
        first = int(atom_indices[bad[0]])
        raise GeometryError(
            f"{bad.size} atom(s) have non-finite coordinates and cannot be written; the "
            f"first is atom {first} of residue "
            f"{structure.residue_label(int(structure.residue_index[first]))}."
        )


def _real(value: float, default: float) -> tuple[str, bool]:
    """Format an occupancy or B-factor, substituting ``default`` for a non-finite value.

    Returns the text and whether a substitution happened, so the caller can record it the
    way the PDB writer records the same substitution rather than let it pass silently.
    """
    if not np.isfinite(value):
        return f"{default:.2f}", True
    return f"{value:.2f}", False


def _atom_rows(
    structure: Structure,
    atom_indices: np.ndarray,
    *,
    model_num: int,
    first_serial: int,
    annotate_regions: bool,
) -> tuple[list[list[str]], int, int]:
    """Render one model's ``_atom_site`` rows and return them with the next free serial.

    Serials continue across models (RCSB numbers them globally), so ``first_serial`` is
    threaded in and the next free value is returned. The third return value counts
    occupancy/B-factor substitutions, aggregated by the caller into a single note.
    """
    element = _elements(structure)
    occupancy = structure.occupancy
    # Occupancy and B-factor are both per residue on a Structure (taken from the alpha
    # carbon), so every atom of a residue is written with that residue's value -- the same
    # thing the PDB writer does. With annotate_regions the B-factor column is replaced by
    # the region annotation; reusing the PDB writer's function keeps the two formats
    # identical, BETA_FOLDED/BETA_DISORDERED polarity and "residue in no domain" note alike.
    b_per_residue = _residue_b_factors(structure, annotate_regions=True) if annotate_regions else None

    residue_index = structure.residue_index
    residue_name = structure.residue_name
    residue_number = structure.residue_number
    insertion_code = structure.insertion_code
    chain_index = structure.chain_index
    atom_name = structure.atom_name
    xyz = structure.xyz

    model_field = str(model_num)
    rows: list[list[str]] = []
    substitutions = 0
    serial = first_serial
    for atom in atom_indices.tolist():
        residue = int(residue_index[atom])
        name = str(residue_name[residue])
        group = "ATOM" if name in THREE_TO_ONE else "HETATM"
        chain_id = _encode_value(_chain_id_of(structure, int(chain_index[residue])))
        entity = str(int(chain_index[residue]) + 1)
        seq = str(int(residue_number[residue]))
        icode_raw = str(insertion_code[residue])
        icode = _encode_value(icode_raw) if icode_raw else _UNKNOWN
        atom_field = _encode_value(str(atom_name[atom]))
        comp_field = _encode_value(name)
        symbol = str(element[atom]) or _element_fallback(str(atom_name[atom]))

        occ_field, occ_sub = _real(float(occupancy[residue]), _DEFAULT_OCCUPANCY)
        b_source = b_per_residue if b_per_residue is not None else structure.b_factor
        b_field, b_sub = _real(float(b_source[residue]), _DEFAULT_B_FACTOR)
        substitutions += int(occ_sub) + int(b_sub)

        x, y, z = (float(v) for v in xyz[atom])
        rows.append(
            [
                group,
                str(serial),
                _encode_value(symbol.upper()),
                atom_field,
                _INAPPLICABLE,
                comp_field,
                chain_id,
                entity,
                seq,
                icode,
                f"{x:.3f}",
                f"{y:.3f}",
                f"{z:.3f}",
                occ_field,
                b_field,
                seq,
                comp_field,
                chain_id,
                atom_field,
                model_field,
            ]
        )
        serial += 1
    return rows, serial, substitutions


def _element_fallback(atom_name: str) -> str:
    """First alphabetic character of an atom name, for the rare atom with no element."""
    return next((ch for ch in atom_name if ch.isalpha()), "").upper()


def _chain_id_of(structure: Structure, chain: int) -> str:
    """The chain id string for a chain index, generating one only if the view is missing.

    mmCIF keeps chain ids verbatim -- ``label_asym_id`` may be several characters, so
    unlike the PDB writer there is no remap onto single characters and no chain is ever
    merged to fit a column. A chain index with no :class:`~dodo.structure.Chain` view is
    the one case a value has to be invented; it is given a numeric id and recorded.
    """
    if chain < len(structure.chains):
        return structure.chains[chain].chain_id or str(chain + 1)
    _note(
        structure,
        f"chain_index refers to chain {chain} but no Chain view exists for it; it was "
        f"given the generated mmCIF chain id {chain + 1}.",
    )
    return str(chain + 1)


def _struct_conn_rows(structure: Structure, atom_indices: np.ndarray) -> list[list[str]]:
    """Render ``_struct_conn`` rows for the backbone/CA-CA bonds of one topology.

    The bonds are :func:`dodo.io.write._bonds`, so this connectivity is identical to the
    PDB writer's ``CONECT`` set. Written once for all models, keyed by label identifier.
    """
    residue_index = structure.residue_index
    residue_name = structure.residue_name
    residue_number = structure.residue_number
    insertion_code = structure.insertion_code
    chain_index = structure.chain_index
    atom_name = structure.atom_name

    def partner(position: int) -> list[str]:
        atom = int(atom_indices[position])
        residue = int(residue_index[atom])
        icode_raw = str(insertion_code[residue])
        return [
            _encode_value(_chain_id_of(structure, int(chain_index[residue]))),
            _encode_value(str(residue_name[residue])),
            str(int(residue_number[residue])),
            _encode_value(str(atom_name[atom])),
            _encode_value(icode_raw) if icode_raw else _UNKNOWN,
        ]

    rows: list[list[str]] = []
    for index, (left, right) in enumerate(_bonds(structure, atom_indices), start=1):
        rows.append([f"covale{index}", "covale", *partner(left), *partner(right)])
    return rows


def _entity_poly_seq_rows(structure: Structure, written: np.ndarray) -> list[list[str]]:
    """Render ``_entity_poly_seq`` rows -- the mmCIF equivalent of SEQRES.

    Describes the residues actually written, grouped by chain (its entity id) and
    numbered from 1 within each, so the sequence block cannot claim residues the
    coordinates do not contain. Interleaved runs of one chain id are concatenated, as the
    PDB writer's SEQRES does, because one entity cannot be declared twice.
    """
    chain_index = structure.chain_index
    residue_name = structure.residue_name
    counts: dict[int, int] = {}
    rows: list[list[str]] = []
    for residue in range(structure.n_residues):
        if not written[residue]:
            continue
        entity = int(chain_index[residue]) + 1
        counts[entity] = counts.get(entity, 0) + 1
        rows.append(
            [str(entity), str(counts[entity]), _encode_value(str(residue_name[residue])), "n"]
        )
    return rows


def _cell_lines(box: tuple[float, float, float]) -> list[str]:
    """Render ``_cell`` and ``_symmetry`` for an orthorhombic P 1 cell.

    The mmCIF counterpart of the PDB writer's CRYST1: written only when a caller asks for
    a box, since DODO does no periodic boundaries and most viewers do not want a cell.
    """
    a, b, c = box
    for length in (a, b, c):
        if not np.isfinite(length) or length <= 0:
            raise GeometryError(f"cell dimensions must be positive and finite, got {box!r}.")
    return [
        f"_cell.length_a {a:.3f}",
        f"_cell.length_b {b:.3f}",
        f"_cell.length_c {c:.3f}",
        "_cell.angle_alpha 90.00",
        "_cell.angle_beta 90.00",
        "_cell.angle_gamma 90.00",
        "_symmetry.space_group_name_H-M 'P 1'",
    ]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def structure_to_cif_text(
    structure: Structure | Sequence[Structure],
    *,
    conect: bool = True,
    annotate_regions: bool = False,
    ca_only: bool = False,
    box: tuple[float, float, float] | None = None,
    seqres: bool = False,
) -> str:
    """Render one structure, or an ensemble, as a single mmCIF document.

    A sequence of structures becomes one ``_atom_site`` loop whose rows carry
    ``pdbx_PDB_model_num`` -- the mmCIF way to hold several conformers of one topology.
    The frames must therefore share that topology, checked exactly as the PDB writer
    checks it. A single :class:`~dodo.structure.Structure` is written as model 1.

    Parameters
    ----------
    structure
        One structure or a sequence of same-topology structures.
    conect
        Emit ``_struct_conn`` connectivity. On by default; see the module docstring for
        why a CA-only trace needs it and which viewers honour it.
    annotate_regions
        Replace the B-factor column with the region annotation (folded
        :data:`~dodo.constants.BETA_FOLDED`, disordered
        :data:`~dodo.constants.BETA_DISORDERED`), identically to the PDB writer.
    ca_only
        Write exactly one alpha carbon per residue -- the atom
        :attr:`~dodo.structure.Structure.ca_indices` resolves to.
    box
        Cell dimensions for a ``_cell``/``_symmetry`` block, or ``None`` for neither.
    seqres
        Emit an ``_entity_poly_seq`` block for the residues written.

    Returns
    -------
    str
        Complete mmCIF text, ending in a newline.

    Raises
    ------
    EmptyStructureError
        If an empty sequence is given.
    GeometryError
        If any coordinate is non-finite, or ``box`` is not positive and finite.
    StructureFileError
        If a value cannot be represented as a STAR token (see :func:`_encode_value`).
    """
    frames = [structure] if isinstance(structure, Structure) else list(structure)
    if not frames:
        raise EmptyStructureError("No structures given to write.")
    if len(frames) > 1:
        _require_matching_topology(frames)

    reference = frames[0]
    atom_indices = _selected_atom_indices(reference, ca_only=ca_only)
    for frame in frames:
        _reject_non_finite(frame, atom_indices)

    written = np.zeros(reference.n_residues, dtype=bool)
    written[reference.residue_index[atom_indices]] = True

    atom_rows: list[list[str]] = []
    substitutions = 0
    serial = 1
    for model_num, frame in enumerate(frames, start=1):
        # ca_only selection depends only on atom names and residue grouping, which the
        # topology check has already proved identical across frames, so the reference
        # selection applies to every frame and the serial keeps climbing.
        rows, serial, subs = _atom_rows(
            frame,
            atom_indices,
            model_num=model_num,
            first_serial=serial,
            annotate_regions=annotate_regions,
        )
        atom_rows.extend(rows)
        substitutions += subs
    if substitutions:
        _note(
            reference,
            f"{substitutions} occupancy/B-factor value(s) were missing and the "
            f"conventional defaults ({_DEFAULT_OCCUPANCY:.2f} / {_DEFAULT_B_FACTOR:.2f}) "
            f"were written.",
        )

    lines = [f"data_{_BLOCK_NAME}", "#"]
    if box is not None:
        lines.extend(_cell_lines(box))
        lines.append("#")
    if seqres:
        lines.extend(_loop(("entity_id", "num", "mon_id", "hetero"), "entity_poly_seq",
                           _entity_poly_seq_rows(reference, written)))
        lines.append("#")
    lines.extend(_loop(_ATOM_SITE_COLUMNS, "atom_site", atom_rows))
    lines.append("#")
    if conect:
        conn_rows = _struct_conn_rows(reference, atom_indices)
        if conn_rows:
            lines.extend(_loop(_STRUCT_CONN_COLUMNS, "struct_conn", conn_rows))
            lines.append("#")
    return "\n".join(lines) + "\n"


def write_cif(
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
    """Write one or more structures to an mmCIF file.

    The parameters mirror :func:`dodo.io.write.write_pdb`. The one behavioural difference
    is how an ensemble is stored: ``models_as_frames=True`` (the default) writes the
    frames as one file whose ``_atom_site`` rows carry ``pdbx_PDB_model_num``, which is
    the mmCIF equivalent of PDB's ``MODEL``/``ENDMDL`` frames; ``False`` writes one file
    per structure, named ``<stem>_1<suffix>``, ``<stem>_2<suffix>``, and so on.

    As with the PDB writer, whether an ensemble is written depends on being *given* a
    sequence, not on its length: a one-element sequence is still an ensemble and still
    gets the ``_1`` suffix under ``models_as_frames=False``, while a bare
    :class:`~dodo.structure.Structure` is always written as a single model.

    Parameters
    ----------
    structure
        A single structure or a sequence of same-topology structures.
    path
        Destination file. A bare relative name works.
    conect, annotate_regions, ca_only, box, seqres
        See :func:`structure_to_cif_text`.
    models_as_frames
        See above.

    Raises
    ------
    EmptyStructureError
        If an empty sequence is given.
    StructureFileError
        If ``path`` is not writable, or a value cannot be represented; plus everything
        :func:`structure_to_cif_text` raises.
    """
    one_structure = isinstance(structure, Structure)
    frames = [structure] if one_structure else list(structure)
    if not frames:
        raise EmptyStructureError("No structures given to write.")
    target = Path(path)

    if one_structure or models_as_frames:
        _write_text(
            target,
            structure_to_cif_text(
                structure,
                conect=conect,
                annotate_regions=annotate_regions,
                ca_only=ca_only,
                box=box,
                seqres=seqres,
            ),
        )
        return

    for index, frame in enumerate(frames, start=1):
        _write_text(
            target.with_name(f"{target.stem}_{index}{target.suffix}"),
            structure_to_cif_text(
                frame,
                conect=conect,
                annotate_regions=annotate_regions,
                ca_only=ca_only,
                box=box,
                seqres=seqres,
            ),
        )


def _write_text(path: Path, text: str) -> None:
    """Write mmCIF text to ``path``, raising :class:`StructureFileError` on any failure.

    The parent directory is checked but never created -- silently creating directories
    turns a mistyped path into a file the user cannot find -- and an existing directory
    passed as ``path`` is reported as such rather than surfacing as ``IsADirectoryError``.
    This mirrors the PDB writer's contract, minus its ASCII check: mmCIF is UTF-8, so a
    multi-character or non-ASCII chain id (which the PDB writer must remap or reject) is
    written verbatim.
    """
    parent = path.parent
    if not parent.is_dir():
        raise StructureFileError(
            f"Cannot write {str(path)!r}: the directory {str(parent)!r} does not exist."
        )
    if path.is_dir():
        raise StructureFileError(
            f"Cannot write {str(path)!r}: it is an existing directory, not a file."
        )
    try:
        path.write_text(text, encoding="utf-8")
    except OSError as error:
        raise StructureFileError(f"Cannot write {str(path)!r}: {error}.") from error
