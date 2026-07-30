"""mmCIF / PDBx reading, on top of a real STAR tokenizer.

Why a tokenizer at all
----------------------
The reader this replaces found the ``_atom_site`` loop with line-prefix heuristics: it
looked for lines beginning ``ATOM``/``HETATM``, assumed ``_atom_site.group_PDB`` was the
first tag in the loop (and raised ``ValueError`` on any file where it was not), split
rows on whitespace, and stopped at the first blank or ``#`` line. Every one of those
assumptions is false for some legal file, and each failure was silent or misattributed:

* tags in a loop may appear in any order, so columns must be resolved by tag *name*;
* a row may be spread over several lines, so "one line is one row" loses atoms;
* comment and blank lines are permitted inside a loop, so stopping at one truncates it;
* values may be quoted and contain whitespace (6kn7's ADP atom names are ``"O5'"``), so
  naive splitting shifts every subsequent column;
* semicolon-delimited text fields span lines, and that is how a deposited sequence is
  written, which is why the old code could not read one and fetched it over HTTP instead.

So this module tokenizes properly. That costs a few hundred lines, and it buys back both
the heuristics and the HTTP module that existed only to fetch deposited sequences.

Author versus label identifiers
-------------------------------
Chain and residue identity come from ``auth_asym_id`` / ``auth_seq_id``, falling back to
``label_*`` only when the author fields are absent. Author numbering is what is printed
in papers, what the equivalent PDB file contains, and therefore what DODO must round
trip. The previous reader took chains from ``label_asym_id`` and residue numbers from
``auth_seq_id``, which produced chain ids disagreeing with the same entry's PDB file
(in 6kn7, the ADP bound to author chain ``O`` sits in label chain ``RA``).
"""

from __future__ import annotations

import gzip
import re
from collections import Counter
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Final, NamedTuple

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
    "parse_cif_text",
    "read_cif",
    "read_unobserved_residues",
]

# ---------------------------------------------------------------------------
# Tokenizer
# ---------------------------------------------------------------------------

_TAG: Final = "tag"
_VALUE: Final = "value"
_NULL: Final = "null"  # unquoted '.' (inapplicable) or '?' (unknown)
_LOOP: Final = "loop"
_DATA: Final = "data"
_CONTROL: Final = "control"  # save_ / stop_ / global_

_WHITESPACE: Final = " \t"


class _Token(NamedTuple):
    """One STAR token, with the 1-based line it started on."""

    kind: str
    value: str
    line: int


def _classify(raw: str, line: int) -> _Token:
    """Classify an unquoted whitespace-delimited word."""
    first = raw[0]
    if first == "_":
        return _Token(_TAG, raw, line)
    if len(raw) == 1 and (first == "." or first == "?"):
        return _Token(_NULL, raw, line)
    lowered = raw.lower()
    if lowered == "loop_":
        return _Token(_LOOP, raw, line)
    if lowered.startswith("data_"):
        return _Token(_DATA, raw[5:], line)
    if lowered.startswith("save_") or lowered in ("stop_", "global_"):
        return _Token(_CONTROL, raw, line)
    return _Token(_VALUE, raw, line)


def _scan_line(line: str, lineno: int) -> Iterator[_Token]:
    """Tokenize one line that contains a quote or a possible comment."""
    i = 0
    end = len(line)
    while i < end:
        ch = line[i]
        if ch in _WHITESPACE:
            i += 1
            continue
        # A '#' only starts a comment at the start of a token, i.e. at line start or
        # after whitespace. Inside an unquoted value it is an ordinary character.
        if ch == "#":
            return
        if ch == "'" or ch == '"':
            j = i + 1
            while True:
                j = line.find(ch, j)
                if j < 0:
                    raise MalformedRecordError(f"Unterminated {ch}-quoted value", lineno, line)
                # Per CIF 1.1 the closing quote must be followed by whitespace or the
                # end of the line. That rule is what keeps 6kn7's ADP atom name
                # "O5'" -- and unquoted names like O5' -- from being cut in half.
                if j + 1 >= end or line[j + 1] in _WHITESPACE:
                    break
                j += 1
            yield _Token(_VALUE, line[i + 1 : j], lineno)
            i = j + 1
            continue
        j = i
        while j < end and line[j] not in _WHITESPACE:
            j += 1
        yield _classify(line[i:j], lineno)
        i = j


def _tokenize(text: str) -> Iterator[_Token]:
    """Yield STAR tokens from mmCIF text.

    Handles single- and double-quoted values, semicolon-delimited multi-line text
    fields, comments, blank lines, and the ``.``/``?`` placeholders.

    Parameters
    ----------
    text
        Complete mmCIF text.

    Yields
    ------
    _Token
        Tokens in file order.
    """
    lines = text.splitlines()
    total = len(lines)
    index = 0
    while index < total:
        line = lines[index]
        lineno = index + 1
        index += 1
        if not line or line.isspace():
            continue
        if line[0] == ";":
            # A semicolon text field: the value starts after the ';' and runs until a
            # line whose first character is ';'. This is how deposited sequences and
            # long author lists are written.
            chunks = [line[1:]]
            closed = False
            trailer = ""
            trailer_line = lineno
            while index < total:
                nxt = lines[index]
                index += 1
                if nxt[:1] == ";":
                    closed = True
                    # Per CIF 1.1 the closing ';' need only be followed by whitespace or
                    # the end of the line: any further values on it are legal tokens, and
                    # a loop row may well continue there. Dropping the rest of the line
                    # loses them with no error and no note.
                    trailer = nxt[1:]
                    trailer_line = index
                    break
                chunks.append(nxt)
            if not closed:
                raise MalformedRecordError(
                    "Unterminated semicolon-delimited text field", lineno, line
                )
            yield _Token(_VALUE, "\n".join(chunks), lineno)
            if trailer.strip():
                # Scanned character by character rather than split, so a quoted value or
                # a comment after the closing ';' behaves as it would anywhere else.
                yield from _scan_line(trailer, trailer_line)
            continue
        stripped = line.lstrip()
        if stripped[0] == "#":
            continue
        # Fast path for the overwhelming majority of lines, including nearly every
        # _atom_site row: no quoting and no comment means plain whitespace splitting is
        # correct, and it is several times quicker than scanning character by character.
        if "'" not in line and '"' not in line and "#" not in line:
            for raw in line.split():
                yield _classify(raw, lineno)
            continue
        yield from _scan_line(line, lineno)


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class _Loop:
    """A parsed CIF loop (or a single-row category promoted to loop shape).

    ``names`` are lowercased tag names with the category prefix removed, so a column is
    always addressed by name and never by position.
    """

    category: str
    names: list[str]
    index: dict[str, int]
    rows: list[list[str | None]]
    #: 1-based line each row started on, or 0 when unknown (a promoted single row).
    row_lines: list[int]

    def column(self, *names: str) -> int:
        """Return the position of the first present column, or -1 if none are."""
        for name in names:
            position = self.index.get(name)
            if position is not None:
                return position
        return -1

    def dict_rows(self) -> Iterator[tuple[dict[str, str | None], int]]:
        """Yield ``(row_by_name, line_number)`` for small categories."""
        for row, line in zip(self.rows, self.row_lines, strict=True):
            yield dict(zip(self.names, row, strict=True)), line


@dataclass(slots=True)
class _Block:
    """One ``data_`` block: its non-loop items and its loops, keyed by category."""

    name: str
    items: dict[str, str | None]
    loops: dict[str, _Loop]


def _parse_blocks(text: str, *, categories: frozenset[str] | None = None) -> list[_Block]:
    """Parse mmCIF text into data blocks.

    Parameters
    ----------
    text
        Complete mmCIF text.
    categories
        Lowercased category names to retain, e.g. ``{"atom_site"}``. Values of other
        categories are still tokenized and validated but not stored, which is what keeps
        a 6.8 MB entry from costing several hundred megabytes of Python strings. Pass
        ``None`` to keep everything.

    Returns
    -------
    list of _Block
        Blocks in file order.
    """
    blocks: list[_Block] = []
    current: _Block | None = None
    tokens = _tokenize(text)
    pending: _Token | None = None

    while True:
        token = pending if pending is not None else next(tokens, None)
        pending = None
        if token is None:
            break

        if token.kind == _DATA:
            current = _Block(token.value, {}, {})
            blocks.append(current)
            continue

        if current is None:
            raise MalformedRecordError(
                f"Found {token.value!r} before any data_ block; this is not mmCIF",
                token.line,
            )

        if token.kind == _CONTROL:
            # save_ frames occur in dictionaries, not in coordinate files. Ignoring the
            # marker (rather than its contents) keeps their contents readable as
            # ordinary items, which is harmless here since we only read named
            # categories.
            continue

        if token.kind == _LOOP:
            pending = _parse_loop(token, tokens, current, categories)
            continue

        if token.kind == _TAG:
            value_token = next(tokens, None)
            if value_token is None or value_token.kind not in (_VALUE, _NULL):
                raise MalformedRecordError(
                    f"Tag {token.value} has no value", token.line, token.value
                )
            key = token.value[1:].lower()
            if categories is not None and key.partition(".")[0] not in categories:
                continue
            if key in current.items:
                raise MalformedRecordError(
                    f"Tag {token.value} appears twice in data block {current.name!r}",
                    token.line,
                )
            current.items[key] = None if value_token.kind == _NULL else value_token.value
            continue

        raise MalformedRecordError(
            f"Unexpected value {token.value!r} where a tag or loop_ was expected",
            token.line,
        )

    return blocks


def _parse_loop(
    loop_token: _Token,
    tokens: Iterator[_Token],
    block: _Block,
    categories: frozenset[str] | None,
) -> _Token | None:
    """Consume one ``loop_`` and store it on ``block``; return the first token after it."""
    names: list[str] = []
    seen_categories: set[str] = set()
    token = next(tokens, None)
    while token is not None and token.kind == _TAG:
        category, dot, name = token.value[1:].partition(".")
        if not dot:
            # Plain CIF (no category prefix). mmCIF always has one; keep the whole tag
            # as the name so nothing is lost if we meet such a file.
            category, name = "", token.value[1:]
        names.append(name.lower())
        seen_categories.add(category.lower())
        token = next(tokens, None)

    if not names:
        raise MalformedRecordError("loop_ with no tags", loop_token.line)
    if len(seen_categories) > 1:
        raise MalformedRecordError(
            f"loop_ mixes categories {sorted(seen_categories)}; each loop describes one",
            loop_token.line,
        )
    if len(set(names)) != len(names):
        duplicated = sorted({n for n in names if names.count(n) > 1})
        raise MalformedRecordError(f"loop_ repeats tag(s) {duplicated}", loop_token.line)

    category = seen_categories.pop()
    keep = categories is None or category in categories
    width = len(names)
    rows: list[list[str | None]] = []
    row_lines: list[int] = []
    buffer: list[str | None] = []
    buffer_line = loop_token.line

    while token is not None and (token.kind == _VALUE or token.kind == _NULL):
        if not buffer:
            buffer_line = token.line
        buffer.append(None if token.kind == _NULL else token.value)
        if len(buffer) == width:
            if keep:
                rows.append(buffer)
                row_lines.append(buffer_line)
            buffer = []
        token = next(tokens, None)

    if buffer:
        raise MalformedRecordError(
            f"_{category} loop ended with {len(buffer)} value(s) for {width} tags; "
            f"the last row is incomplete",
            buffer_line,
        )

    if keep:
        if category in block.loops:
            raise MalformedRecordError(
                f"Category _{category} appears in two loops in data block {block.name!r}",
                loop_token.line,
            )
        block.loops[category] = _Loop(
            category=category,
            names=names,
            index={name: position for position, name in enumerate(names)},
            rows=rows,
            row_lines=row_lines,
        )
    return token


def _category_loop(block: _Block, category: str) -> _Loop | None:
    """Return ``category`` as a loop, promoting a single-row (non-loop) category.

    A category with one row may legally be written as plain ``_cat.tag value`` items --
    which is exactly how AlphaFold writes ``_entity_poly`` for a monomer. Promoting it
    means every caller sees one shape.
    """
    loop = block.loops.get(category)
    if loop is not None:
        return loop
    prefix = f"{category}."
    names = [key[len(prefix) :] for key in block.items if key.startswith(prefix)]
    if not names:
        return None
    row: list[str | None] = [block.items[prefix + name] for name in names]
    return _Loop(
        category=category,
        names=names,
        index={name: position for position, name in enumerate(names)},
        rows=[row],
        row_lines=[0],
    )


# ---------------------------------------------------------------------------
# Value conversion
# ---------------------------------------------------------------------------


def _line_or_none(line: int) -> int | None:
    return line or None


def _to_float(value: str | None, *, default: float, field: str, line: int) -> float:
    if value is None:
        return default
    try:
        return float(value)
    except ValueError:
        raise MalformedRecordError(
            f"_atom_site.{field} is not a number: {value!r}", _line_or_none(line)
        ) from None


def _to_int(value: str | None, *, field: str, line: int) -> int:
    if value is None:
        raise MalformedRecordError(
            f"_atom_site.{field} is missing (. or ?), so this atom has no residue number",
            _line_or_none(line),
        )
    try:
        return int(value)
    except ValueError:
        raise MalformedRecordError(
            f"_atom_site.{field} is not an integer: {value!r}", _line_or_none(line)
        ) from None


def _coordinate(value: str | None, *, field: str, line: int) -> float:
    if value is None:
        raise MalformedRecordError(
            f"_atom_site.{field} is a placeholder ({field} unknown), so this atom has "
            f"no position; DODO will not invent one",
            _line_or_none(line),
        )
    try:
        return float(value)
    except ValueError:
        raise MalformedRecordError(
            f"_atom_site.{field} is not a number: {value!r}", _line_or_none(line)
        ) from None


def _element_from_atom_name(atom_name: str) -> str:
    """Guess an element symbol from an atom name, for files with no ``type_symbol``.

    Two-letter guesses are only accepted when they name an element DODO knows a mass
    for, so ``SE`` of selenomethionine survives while ``CA`` (alpha carbon) does not
    become calcium.
    """
    letters = "".join(ch for ch in atom_name if ch.isalpha()).upper()
    if not letters:
        return "C"
    if len(letters) >= 2 and letters[:2] in ATOMIC_MASSES:
        return letters[:2]
    return letters[0]


def _summarize(counter: Counter[str], limit: int = 5) -> str:
    common = counter.most_common(limit)
    text = ", ".join(f"{name} x{count}" for name, count in common)
    if len(counter) > limit:
        text += f", and {len(counter) - limit} more kind(s)"
    return text


# ---------------------------------------------------------------------------
# _atom_site
# ---------------------------------------------------------------------------


@dataclass(slots=True)
class _AtomTable:
    """Flat per-atom records, in file order, ready for ``from_atom_records``."""

    xyz: list[tuple[float, float, float]]
    atom_name: list[str]
    element: list[str]
    residue_name: list[str]
    residue_number: list[int]
    chain_id: list[str]
    insertion_code: list[str]
    b_factor: list[float]
    occupancy: list[float]


def _extract_atoms(
    loop: _Loop,
    *,
    model: int | None,
    keep_hydrogens: bool,
    source: str | None,
) -> tuple[_AtomTable, list[str]]:
    """Turn an ``_atom_site`` loop into flat per-atom records plus provenance notes."""
    notes: list[str] = []

    x_col = loop.column("cartn_x")
    y_col = loop.column("cartn_y")
    z_col = loop.column("cartn_z")
    missing = [
        tag for tag, col in (("Cartn_x", x_col), ("Cartn_y", y_col), ("Cartn_z", z_col)) if col < 0
    ]
    if missing:
        raise StructureFileError(
            f"_atom_site has no {', '.join(missing)} column(s), so it carries no "
            f"coordinates (source: {source!r})."
        )

    name_col = loop.column("label_atom_id", "auth_atom_id")
    if name_col < 0:
        raise StructureFileError(
            f"_atom_site has neither label_atom_id nor auth_atom_id, so atoms cannot be "
            f"named (source: {source!r})."
        )
    comp_col = loop.column("label_comp_id", "auth_comp_id")
    if comp_col < 0:
        raise StructureFileError(
            f"_atom_site has neither label_comp_id nor auth_comp_id, so residues cannot "
            f"be typed (source: {source!r})."
        )

    # Author identifiers first: they are what the equivalent PDB file uses and what DODO
    # writes back out. label_* is a fallback, not a peer.
    chain_col = loop.column("auth_asym_id")
    chain_tag = "auth_asym_id"
    if chain_col < 0:
        chain_col = loop.column("label_asym_id")
        chain_tag = "label_asym_id"
        if chain_col >= 0:
            notes.append(
                "_atom_site has no auth_asym_id; chain ids were taken from "
                "label_asym_id and may disagree with the equivalent PDB file."
            )
    if chain_col < 0:
        raise StructureFileError(
            f"_atom_site has neither auth_asym_id nor label_asym_id, so atoms cannot be "
            f"assigned to chains (source: {source!r})."
        )

    seq_col = loop.column("auth_seq_id")
    seq_tag = "auth_seq_id"
    if seq_col < 0:
        seq_col = loop.column("label_seq_id")
        seq_tag = "label_seq_id"
        if seq_col >= 0:
            notes.append(
                "_atom_site has no auth_seq_id; residue numbers were taken from "
                "label_seq_id, which counts from 1 within each entity and therefore "
                "does not match the numbering used in publications."
            )
    if seq_col < 0:
        raise StructureFileError(
            f"_atom_site has neither auth_seq_id nor label_seq_id, so residues cannot be "
            f"numbered (source: {source!r})."
        )

    group_col = loop.column("group_pdb")
    symbol_col = loop.column("type_symbol")
    alt_col = loop.column("label_alt_id")
    icode_col = loop.column("pdbx_pdb_ins_code")
    occupancy_col = loop.column("occupancy")
    b_col = loop.column("b_iso_or_equiv")
    model_col = loop.column("pdbx_pdb_model_num")
    serial_col = loop.column("id")

    if symbol_col < 0:
        notes.append(
            "_atom_site has no type_symbol; element symbols were inferred from atom names."
        )
    if group_col < 0:
        notes.append(
            "_atom_site has no group_PDB; every record was treated as ATOM, so any "
            "ligand or solvent whose component is a standard residue code was kept."
        )

    selected_model: str | None = None
    if model_col >= 0:
        if model is None:
            for row in loop.rows:
                if row[model_col] is not None:
                    selected_model = row[model_col]
                    break
        else:
            selected_model = str(model)
    elif model is not None and model != 1:
        raise StructureFileError(
            f"model={model} was requested but _atom_site has no pdbx_PDB_model_num "
            f"column, so the file holds a single unnumbered model (source: {source!r})."
        )

    table = _AtomTable([], [], [], [], [], [], [], [], [])
    models_seen: set[str] = set()
    skipped_model = 0
    skipped_hydrogen = 0
    skipped_altloc = 0
    skipped_component: Counter[str] = Counter()
    nonstandard_atom_records: Counter[str] = Counter()
    inferred_elements = 0

    # Which alternate was kept for each residue *identity*, not for each contiguous run
    # of records. A file ordered by conformer (every altloc A atom, then every altloc B
    # atom) is legal, and run-based tracking keeps every alternate in one -- doubling the
    # residue count, giving each residue two CA atoms, and reporting nothing because no
    # alternate was ever seen twice in a row. Only residues that actually carry an
    # alternate get an entry, so this costs nothing on the overwhelming majority of files.
    kept_altloc: dict[tuple[str, int, str], str] = {}

    for row, line in zip(loop.rows, loop.row_lines, strict=True):
        if model_col >= 0:
            row_model = row[model_col]
            if row_model is not None:
                models_seen.add(row_model)
            if row_model != selected_model:
                skipped_model += 1
                continue

        residue_name = row[comp_col]
        if residue_name is None:
            raise MalformedRecordError(
                "_atom_site row has no component id, so its residue cannot be typed",
                _line_or_none(line),
            )
        residue_name = residue_name.upper()

        group = (row[group_col] or "ATOM").upper() if group_col >= 0 else "ATOM"
        if group == "HETATM" and residue_name not in THREE_TO_ONE:
            # Solvent, ions, nucleotides, ligands. Dropping these is the point of the
            # filter; MSE and friends are in THREE_TO_ONE and survive.
            skipped_component[residue_name] += 1
            continue
        if group not in ("ATOM", "HETATM"):
            raise MalformedRecordError(
                f"_atom_site.group_PDB is {group!r}; expected ATOM or HETATM",
                _line_or_none(line),
            )
        if group == "ATOM" and residue_name not in THREE_TO_ONE:
            nonstandard_atom_records[residue_name] += 1

        atom_name = row[name_col]
        if atom_name is None:
            serial = row[serial_col] if serial_col >= 0 else None
            raise MalformedRecordError(
                f"_atom_site row (id {serial}) has no atom name",
                _line_or_none(line),
            )

        symbol = row[symbol_col] if symbol_col >= 0 else None
        if symbol:
            element = symbol.upper()
        else:
            element = _element_from_atom_name(atom_name)
            inferred_elements += 1

        if not keep_hydrogens and element in ("H", "D"):
            skipped_hydrogen += 1
            continue

        icode_raw = row[icode_col] if icode_col >= 0 else None
        insertion_code = icode_raw or ""

        residue_number = _to_int(row[seq_col], field=seq_tag, line=line)
        chain_raw = row[chain_col]
        if chain_raw is None:
            raise MalformedRecordError(
                f"_atom_site.{chain_tag} is a placeholder, so this atom belongs to no chain",
                _line_or_none(line),
            )

        altloc = (row[alt_col] if alt_col >= 0 else None) or ""
        if altloc:
            # Structure has no altloc axis, and two conformers of one residue would give
            # it two CA atoms -- which silently makes ca_xyz pick whichever came last. So
            # keep the first alternate that appears for each residue identity and count
            # the rest, wherever in the file the other conformer's atoms turn up.
            key = (chain_raw, residue_number, insertion_code)
            if kept_altloc.setdefault(key, altloc) != altloc:
                skipped_altloc += 1
                continue

        table.xyz.append(
            (
                _coordinate(row[x_col], field="Cartn_x", line=line),
                _coordinate(row[y_col], field="Cartn_y", line=line),
                _coordinate(row[z_col], field="Cartn_z", line=line),
            )
        )
        table.atom_name.append(atom_name)
        table.element.append(element)
        table.residue_name.append(residue_name)
        table.residue_number.append(residue_number)
        table.chain_id.append(chain_raw)
        table.insertion_code.append(insertion_code)
        table.occupancy.append(
            _to_float(
                row[occupancy_col] if occupancy_col >= 0 else None,
                default=1.0,
                field="occupancy",
                line=line,
            )
        )
        table.b_factor.append(
            _to_float(
                row[b_col] if b_col >= 0 else None, default=0.0, field="B_iso_or_equiv", line=line
            )
        )

    notes.extend(
        _atom_notes(
            n_rows=len(loop.rows),
            n_kept=len(table.xyz),
            chain_tag=chain_tag,
            seq_tag=seq_tag,
            selected_model=selected_model,
            models_seen=models_seen,
            skipped_model=skipped_model,
            skipped_hydrogen=skipped_hydrogen,
            skipped_altloc=skipped_altloc,
            skipped_component=skipped_component,
            nonstandard_atom_records=nonstandard_atom_records,
            inferred_elements=inferred_elements,
        )
    )

    # Checked before emptiness: "model 7 does not exist" is the useful message, and
    # "nothing survived filtering" is only its consequence.
    if model is not None and model_col >= 0 and str(model) not in models_seen:
        raise StructureFileError(
            f"model={model} is not present; _atom_site holds model(s) "
            f"{', '.join(sorted(models_seen))} (source: {source!r})."
        )

    if not table.xyz:
        raise EmptyStructureError(
            f"No polymer atoms survived filtering of {len(loop.rows)} _atom_site "
            f"record(s) (source: {source!r}). "
            + " ".join(notes[-3:] if notes else ["The loop may be empty."])
        )

    return table, notes


def _atom_notes(
    *,
    n_rows: int,
    n_kept: int,
    chain_tag: str,
    seq_tag: str,
    selected_model: str | None,
    models_seen: set[str],
    skipped_model: int,
    skipped_hydrogen: int,
    skipped_altloc: int,
    skipped_component: Counter[str],
    nonstandard_atom_records: Counter[str],
    inferred_elements: int,
) -> list[str]:
    """Compose the provenance notes for one ``_atom_site`` pass.

    Every atom that was read and every atom that was not is accounted for here. The
    reader this replaces dropped records with no record of having done so, which is how
    a structure could lose half its atoms and still look plausible.
    """
    notes = [
        f"Read {n_kept} atom(s) from {n_rows} _atom_site record(s), using {chain_tag} "
        f"for chain ids and {seq_tag} for residue numbers."
    ]
    # Reported whenever anything was skipped, not only when several models were found: a
    # row whose pdbx_PDB_model_num is a placeholder is also skipped, and that would
    # otherwise be an unexplained gap between records read and atoms kept.
    if len(models_seen) > 1 or skipped_model:
        notes.append(
            f"_atom_site contains {len(models_seen)} model(s) "
            f"({', '.join(sorted(models_seen))}); model {selected_model} was read and "
            f"{skipped_model} atom(s) from other or unnumbered models were skipped. "
            f"Pass model=N to choose another."
        )
    if skipped_component:
        notes.append(
            f"Skipped {sum(skipped_component.values())} HETATM atom(s) of "
            f"{len(skipped_component)} non-polymer component(s): "
            f"{_summarize(skipped_component)}."
        )
    if skipped_hydrogen:
        notes.append(
            f"Skipped {skipped_hydrogen} hydrogen/deuterium atom(s); pass "
            f"keep_hydrogens=True to retain them."
        )
    if skipped_altloc:
        notes.append(
            f"Skipped {skipped_altloc} atom(s) in alternate locations; the first "
            f"alternate listed for each residue was kept."
        )
    if nonstandard_atom_records:
        notes.append(
            f"Kept {sum(nonstandard_atom_records.values())} ATOM atom(s) whose "
            f"component is not a standard amino acid "
            f"({_summarize(nonstandard_atom_records)}); these are most likely nucleic "
            f"acid or sugar residues. They were kept because they are ATOM records, but "
            f"DODO's geometry requires an alpha carbon per residue and will reject them."
        )
    if inferred_elements:
        notes.append(
            f"Inferred the element of {inferred_elements} atom(s) from their names "
            f"because type_symbol was absent or blank."
        )
    return notes


# ---------------------------------------------------------------------------
# Polymer metadata: deposited sequence and UniProt accession
# ---------------------------------------------------------------------------

_PARENTHESIZED: Final = re.compile(r"\(([^)]*)\)")


def _expand_one_letter_code(raw: str) -> str:
    """Collapse ``pdbx_seq_one_letter_code`` to plain one-letter codes.

    The non-canonical form writes modified monomers as ``(MSE)``. Only used as a
    fallback when ``pdbx_seq_one_letter_code_can`` is absent.
    """

    def substitute(match: re.Match[str]) -> str:
        return THREE_TO_ONE.get(match.group(1).strip().upper(), "X")

    return _PARENTHESIZED.sub(substitute, "".join(raw.split()))


def _split_strand_ids(raw: str | None) -> list[str]:
    if raw is None:
        return []
    return [part.strip() for part in raw.split(",") if part.strip()]


def _entity_poly_sequences(
    block: _Block,
) -> tuple[dict[str, str], dict[str, list[str]], list[str]]:
    """Read deposited sequences from ``_entity_poly``.

    Returns
    -------
    dict
        Chain id (author, per ``pdbx_strand_id``) to deposited one-letter sequence.
    dict
        Entity id to its strand ids, reused to map ``_struct_ref`` by entity.
    list of str
        Provenance notes.
    """
    notes: list[str] = []
    sequences: dict[str, str] = {}
    entity_strands: dict[str, list[str]] = {}

    loop = _category_loop(block, "entity_poly")
    if loop is None:
        return sequences, entity_strands, notes

    used_fallback = 0
    for row, line in loop.dict_rows():
        strands = _split_strand_ids(row.get("pdbx_strand_id"))
        entity_id = row.get("entity_id")
        if entity_id is not None:
            entity_strands[entity_id] = strands
        canonical = row.get("pdbx_seq_one_letter_code_can")
        if canonical is not None:
            sequence = "".join(canonical.split())
        else:
            raw = row.get("pdbx_seq_one_letter_code")
            if raw is None:
                notes.append(
                    f"_entity_poly entity {entity_id!r} has no sequence; the deposited "
                    f"sequence for chain(s) {strands or '(unnamed)'} is unavailable."
                )
                continue
            sequence = _expand_one_letter_code(raw)
            used_fallback += 1
        if not strands:
            notes.append(
                f"_entity_poly entity {entity_id!r} has no pdbx_strand_id, so its "
                f"deposited sequence could not be attached to a chain."
            )
            continue
        for chain_id in strands:
            if chain_id in sequences and sequences[chain_id] != sequence:
                notes.append(
                    f"_entity_poly maps chain {chain_id!r} to more than one sequence; "
                    f"the first was kept."
                )
                continue
            sequences[chain_id] = sequence
        del line  # line numbers are only useful for atom records

    if used_fallback:
        notes.append(
            f"{used_fallback} _entity_poly entity/entities had no "
            f"pdbx_seq_one_letter_code_can; the sequence was derived from "
            f"pdbx_seq_one_letter_code, with unrecognized modified monomers as X."
        )
    return sequences, entity_strands, notes


def _struct_ref_uniprot(
    block: _Block, entity_strands: dict[str, list[str]]
) -> tuple[dict[str, str], list[str]]:
    """Map chain ids to UniProt accessions via ``_struct_ref`` / ``_struct_ref_seq``."""
    notes: list[str] = []
    mapping: dict[str, str] = {}

    ref_loop = _category_loop(block, "struct_ref")
    if ref_loop is None:
        return mapping, notes

    accession_by_ref: dict[str, str] = {}
    entity_by_ref: dict[str, str] = {}
    non_uniprot: Counter[str] = Counter()
    for row, _line in ref_loop.dict_rows():
        db_name = (row.get("db_name") or "").upper()
        ref_id = row.get("id")
        accession = row.get("pdbx_db_accession")
        if db_name != "UNP":
            if db_name:
                non_uniprot[db_name] += 1
            continue
        if accession is None:
            notes.append(
                f"_struct_ref {ref_id!r} names UNP but has no pdbx_db_accession, so no "
                f"UniProt id could be recorded."
            )
            continue
        if ref_id is not None:
            accession_by_ref[ref_id] = accession
            entity_id = row.get("entity_id")
            if entity_id is not None:
                entity_by_ref[ref_id] = entity_id

    if non_uniprot:
        notes.append(
            f"_struct_ref cites non-UniProt database(s) ({_summarize(non_uniprot)}); "
            f"only UNP references are recorded."
        )
    if not accession_by_ref:
        return mapping, notes

    seq_loop = _category_loop(block, "struct_ref_seq")
    if seq_loop is not None:
        for row, _line in seq_loop.dict_rows():
            ref_id = row.get("ref_id")
            if ref_id is None or ref_id not in accession_by_ref:
                continue
            for chain_id in _split_strand_ids(row.get("pdbx_strand_id")):
                mapping.setdefault(chain_id, accession_by_ref[ref_id])

    if not mapping:
        # Some depositions carry _struct_ref without _struct_ref_seq. The entity is
        # still enough to reach the chains, via _entity_poly.pdbx_strand_id.
        for ref_id, accession in accession_by_ref.items():
            for chain_id in entity_strands.get(entity_by_ref.get(ref_id, ""), []):
                mapping.setdefault(chain_id, accession)
        if mapping:
            notes.append(
                "_struct_ref_seq was absent or unusable; UniProt accessions were mapped "
                "to chains through _entity_poly.pdbx_strand_id instead."
            )
        else:
            notes.append(
                "_struct_ref names a UniProt accession but it could not be attached to "
                "any chain (no _struct_ref_seq.pdbx_strand_id and no matching entity)."
            )
    return mapping, notes


# ---------------------------------------------------------------------------
# Unobserved residues
# ---------------------------------------------------------------------------

_UNOBS_CATEGORY: Final = "pdbx_unobs_or_zero_occ_residues"


def read_unobserved_residues(text: str) -> dict[str, list[tuple[int, str]]]:
    """Read the depositor's list of residues present in the construct but not modelled.

    ``_pdbx_unobs_or_zero_occ_residues`` is the authoritative statement of what is
    missing from a deposited model: it is written by the depositor, who knows the
    construct, and it needs no alignment to recover. It is the correct starting point
    for finding regions DODO must rebuild.

    Parameters
    ----------
    text
        Complete mmCIF text.

    Returns
    -------
    dict
        Author chain id to ``(residue_number, insertion_code)`` pairs, in file order,
        deduplicated. Empty when the file does not carry the category -- which most
        AlphaFold models and many older entries do not.

        The insertion code is part of the key rather than dropped, because 100 and 100A
        are two residues of the construct and reporting one where two are missing
        understates the gap. Kabat CDR numbering makes unmodelled ``100A``-``100K`` loops
        routine in antibody structures, which are exactly the entries with loops to
        rebuild.

    Notes
    -----
    Rows are merged across models: the category carries ``PDB_model_num``, but a residue
    unmodelled in one model is unmodelled in the deposited construct either way.

    This should eventually live on :class:`~dodo.structure.Structure` (as a per-chain
    array of unmodelled residue numbers) rather than being re-read from text; it is a
    free function only because this module may not modify ``structure.py``.

    Raises
    ------
    MalformedRecordError
        If the category is present but a row lacks a chain id or residue number, or the
        residue number is not an integer. Silently skipping such a row would understate
        what is missing from the model, which is the opposite of this function's job.
    """
    result: dict[str, list[tuple[int, str]]] = {}
    for block in _parse_blocks(text, categories=frozenset({_UNOBS_CATEGORY})):
        loop = _category_loop(block, _UNOBS_CATEGORY)
        if loop is None:
            continue
        for row, line in loop.dict_rows():
            chain_id = row.get("auth_asym_id") or row.get("label_asym_id")
            number = row.get("auth_seq_id") or row.get("label_seq_id")
            if chain_id is None or number is None:
                raise MalformedRecordError(
                    f"_{_UNOBS_CATEGORY} row has no usable chain id and residue number "
                    f"(needs auth_asym_id/auth_seq_id or the label_ equivalents)",
                    _line_or_none(line),
                )
            try:
                residue_number = int(number)
            except ValueError:
                raise MalformedRecordError(
                    f"_{_UNOBS_CATEGORY} residue number is not an integer: {number!r}",
                    _line_or_none(line),
                ) from None
            # A '.' or '?' placeholder arrives as None and means "no insertion code".
            insertion_code = row.get("pdb_ins_code") or ""
            residues = result.setdefault(chain_id, [])
            identity = (residue_number, insertion_code)
            if identity not in residues:
                residues.append(identity)
    return result


# ---------------------------------------------------------------------------
# Public reader
# ---------------------------------------------------------------------------

_READ_CATEGORIES: Final = frozenset(
    {
        "atom_site",
        "entity_poly",
        "struct_ref",
        "struct_ref_seq",
        _UNOBS_CATEGORY,
    }
)


#: A ``data_`` block header. Mandatory in mmCIF: a file without one is not mmCIF,
#: whatever else it contains.
_DATA_BLOCK: Final = re.compile(r"(?mi)^[ \t]*data_\S")
#: Any other STAR construct -- ``loop_``, or a ``_category.tag`` item. Present in every
#: mmCIF and in no PDB file, since PDB lines begin with a fixed-width record name.
_STAR_CONSTRUCT: Final = re.compile(r"(?mi)^[ \t]*(?:loop_[ \t]*$|_[A-Za-z0-9]+\.[A-Za-z0-9_]+)")
#: How many lines to inspect when looking for PDB record names.
_SNIFF_LINES: Final = 200


def _looks_like_pdb(text: str) -> bool:
    """Report whether ``text`` is plainly a PDB-format file rather than mmCIF.

    Keyed on PDB record names, which a column-aligned mmCIF ``_atom_site`` row shares: a
    real RCSB row begins with the literal text ``"ATOM  "``, so this cannot be used to
    tell the formats apart on its own. It is only ever consulted once :data:`_DATA_BLOCK`
    and :data:`_STAR_CONSTRUCT` have established that the text is not mmCIF, and then only
    to make the error message name the right reader. Consulting it *after* a failed parse
    instead -- which is what this module used to do -- turns every diagnostic about a
    legal-but-broken mmCIF into "this is PDB format", discarding the real cause and the
    line number that located it.
    """
    pdb_records = ("HEADER", "ATOM  ", "HETATM", "REMARK", "SEQRES", "CRYST1", "MODEL ")
    return any(line[:6] in pdb_records for line in text.splitlines()[:_SNIFF_LINES])


def _select_block(blocks: list[_Block], source: str | None) -> tuple[_Block, _Loop, list[str]]:
    """Choose the data block holding coordinates, and note any others."""
    notes: list[str] = []
    with_atoms = [
        (block, loop)
        for block, loop in ((b, _category_loop(b, "atom_site")) for b in blocks)
        if loop is not None
    ]
    if not with_atoms:
        raise StructureFileError(
            f"No _atom_site records found in {source or 'the mmCIF text'}; "
            f"data block(s) present: {', '.join(block.name for block in blocks) or 'none'}."
        )
    block, loop = with_atoms[0]
    if len(with_atoms) > 1:
        others = ", ".join(other.name for other, _ in with_atoms)
        notes.append(
            f"{len(with_atoms)} data blocks carry coordinates ({others}); block "
            f"{block.name!r} was read."
        )
    elif len(blocks) > 1:
        notes.append(
            f"Read data block {block.name!r}; {len(blocks) - 1} other block(s) carry no "
            f"coordinates and were ignored."
        )
    return block, loop, notes


#: Length of the sentinel used to measure Structure's fixed-width string arrays. Longer
#: than any identifier a structure file can plausibly carry, so whatever survives the
#: round trip is the stored width.
_PROBE_LENGTH: Final = 16


class _StorageWidths(NamedTuple):
    """Characters of each identifier that ``Structure`` actually keeps.

    ``None`` for a field whose sentinel came back whole, i.e. one that is not clipped at
    :data:`_PROBE_LENGTH` and therefore needs no checking.
    """

    atom_name: int | None
    element: int | None
    residue_name: int | None
    chain_id: int | None
    insertion_code: int | None


@lru_cache(maxsize=1)
def _storage_widths() -> _StorageWidths:
    """Measure ``Structure``'s fixed-width string storage by round-tripping a sentinel.

    ``Structure.from_atom_records`` chooses the numpy dtypes, and this module may not
    change them. Restating the widths here would be a second copy of the same fact, free
    to drift out of step with the first; asking ``Structure`` what survives cannot drift.
    Computed once and cached, since it costs one one-atom structure.
    """
    sentinel = "X" * _PROBE_LENGTH
    probe = Structure.from_atom_records(
        xyz=np.zeros((1, 3), dtype=np.float64),
        atom_name=[sentinel],
        element=[sentinel],
        residue_name=[sentinel],
        residue_number=[1],
        chain_id=[sentinel],
        insertion_code=[sentinel],
    )

    def measured(stored: str) -> int | None:
        return None if len(stored) >= _PROBE_LENGTH else len(stored)

    return _StorageWidths(
        atom_name=measured(str(probe.atom_name[0])),
        element=measured(str(probe.element[0])),
        residue_name=measured(str(probe.residue_name[0])),
        # Chain ids reach a caller through Chain.chain_id rather than an array, so they
        # are measured where a caller would read them.
        chain_id=measured(probe.chains[0].chain_id),
        insertion_code=measured(str(probe.insertion_code[0])),
    )


def _over_wide(values: Sequence[str], width: int | None) -> list[str]:
    """Distinct ``values`` longer than ``width``, sorted; empty when ``width`` is None."""
    if width is None:
        return []
    return sorted({value for value in values if len(value) > width})


def _check_storage_widths(table: _AtomTable, source: str | None) -> list[str]:
    """Check every identifier against the width ``Structure`` stores for it.

    Called *before* a ``Structure`` is built, because clipping part of a residue's
    identity is not a cosmetic loss: two residues differing only in insertion code (10
    and 10A) become one residue holding two alpha carbons, and ``ca_indices`` then keeps
    whichever came last without a word. Chain ids merge whole chains the same way. Both
    therefore raise, and both must be caught here -- by the time ``from_atom_records`` has
    grouped the atoms, the evidence that two residues went in has been destroyed.

    Atom names, element symbols and residue names are labels that nothing keys on, so an
    over-wide one is recorded in ``notes`` rather than raising.
    """
    widths = _storage_widths()
    for field_name, values, width, consequence in (
        (
            "chain_id",
            table.chain_id,
            widths.chain_id,
            "keeping them would merge distinct chains",
        ),
        (
            "insertion_code",
            table.insertion_code,
            widths.insertion_code,
            "keeping them would merge residues that differ only by insertion code, leaving "
            "one residue with two alpha carbons",
        ),
    ):
        offenders = _over_wide(values, width)
        if offenders:
            raise StructureFileError(
                f"{len(offenders)} {field_name} value(s) do not survive storage in "
                f"Structure, whose {field_name} array holds {width} character(s): "
                f"{offenders[:5]}. {consequence.capitalize()}. Widen Structure's "
                f"{field_name} dtype to read this file (source: {source!r})."
            )

    notes: list[str] = []
    for label, field_name, values, width in (
        ("atom name", "atom_name", table.atom_name, widths.atom_name),
        ("element symbol", "element", table.element, widths.element),
        ("residue name", "residue_name", table.residue_name, widths.residue_name),
    ):
        offenders = _over_wide(values, width)
        if offenders:
            notes.append(
                f"{len(offenders)} {label}(s) are longer than the {width} character(s) "
                f"Structure's {field_name} array holds and were truncated, e.g. "
                f"{offenders[:5]}."
            )
    return notes


def _component_id_conflicts(structure: Structure, table: _AtomTable) -> list[str]:
    """Report residues whose atoms disagree about the component id.

    ``from_atom_records`` stores one residue name per residue, taken from that residue's
    first atom, so a residue whose atoms carry two different component ids -- point
    microheterogeneity, or a malformed file -- loses the others. That loss is real and
    belongs in the notes. It is emphatically *not* a fixed-width dtype problem, which is
    what comparing a per-atom set of names against a per-residue array made it look like:
    that comparison reported a three-character component id as "truncated" while never
    reporting the component id that was actually dropped.
    """
    per_atom = np.asarray(table.residue_name, dtype=np.str_)
    first_of_residue = per_atom[structure.residue_atom_offsets[:-1]]
    differing = np.flatnonzero(per_atom != first_of_residue[structure.residue_index])
    if differing.size == 0:
        return []
    examples = "; ".join(
        f"{structure.residue_label(int(structure.residue_index[i]))} has atom "
        f"{table.atom_name[i]} labelled {per_atom[i]}"
        for i in differing[:3]
    )
    return [
        f"{differing.size} atom(s) carry a component id that differs from the one stored "
        f"for their residue ({examples}). Structure keeps one component id per residue, "
        f"taken from its first atom, so the others were dropped."
    ]


def parse_cif_text(
    text: str,
    *,
    source: str | None = None,
    model: int | None = None,
    keep_hydrogens: bool = False,
) -> Structure:
    """Parse mmCIF text into a :class:`~dodo.structure.Structure`.

    Atoms are kept in file order, which is already residue-grouped in any file written
    by a deposition or prediction pipeline. Nothing is sorted and nothing is merged:
    residue grouping and chain construction are left entirely to
    :meth:`~dodo.structure.Structure.from_atom_records`, so an out-of-order file is
    reported through ``structure.notes`` rather than quietly repaired.

    Parameters
    ----------
    text
        Complete mmCIF text.
    source
        Provenance string recorded on the structure, normally a path.
    model
        ``pdbx_PDB_model_num`` to read. ``None`` reads the first model in the file;
        merging models (which the previous reader did) puts every model's atoms in one
        coordinate array and doubles every residue.
    keep_hydrogens
        Keep hydrogen and deuterium atoms. Off by default: DODO rebuilds CA traces and
        then adds backbone atoms, so input hydrogens are dead weight.

    Returns
    -------
    Structure
        With ``notes`` describing everything read and everything skipped, and with
        ``Chain.full_sequence`` and ``Chain.uniprot_id`` filled in where ``_entity_poly``
        and ``_struct_ref`` provide them.

    Raises
    ------
    UnsupportedFormatError
        If the text has no ``data_`` block.
    StructureFileError
        If there is no ``_atom_site`` category, if it lacks a column DODO cannot do
        without, if the requested model is absent, or if a chain id or insertion code is
        too long for ``Structure`` to store without merging distinct chains or residues.
    MalformedRecordError
        If a row cannot be parsed, carrying the line number.
    EmptyStructureError
        If no polymer atoms survive filtering.
    """
    # Format is decided before parsing, and on STAR constructs rather than on what the
    # atom rows look like. Deciding it afterwards, from the shape of the lines, mistakes
    # any column-aligned mmCIF -- which is to say any file RCSB wrote -- for PDB the
    # moment anything else about it is wrong, and throws away the real error with its
    # line number.
    if _DATA_BLOCK.search(text) is None:
        if _STAR_CONSTRUCT.search(text) is None and _looks_like_pdb(text):
            raise UnsupportedFormatError(
                f"{source or 'This text'} is PDB format, not mmCIF; use read_pdb."
            )
        raise UnsupportedFormatError(
            f"{source or 'This text'} contains no data_ block, so it is not mmCIF."
        )

    blocks = _parse_blocks(text, categories=_READ_CATEGORIES)
    if not blocks:
        # A data_ header matched but the tokenizer found none, so the only ``data_`` in
        # the file is inside a text field or a comment.
        raise UnsupportedFormatError(
            f"{source or 'This text'} contains no data_ block, so it is not mmCIF."
        )

    block, atom_site, notes = _select_block(blocks, source)

    table, atom_notes = _extract_atoms(
        atom_site, model=model, keep_hydrogens=keep_hydrogens, source=source
    )
    notes.extend(atom_notes)
    notes.extend(_check_storage_widths(table, source))

    structure = Structure.from_atom_records(
        xyz=np.array(table.xyz, dtype=np.float64),
        atom_name=table.atom_name,
        element=table.element,
        residue_name=table.residue_name,
        residue_number=table.residue_number,
        chain_id=table.chain_id,
        insertion_code=table.insertion_code,
        b_factor=table.b_factor,
        occupancy=table.occupancy,
        source=source,
    )
    notes.extend(_component_id_conflicts(structure, table))

    sequences, entity_strands, sequence_notes = _entity_poly_sequences(block)
    notes.extend(sequence_notes)
    uniprot, uniprot_notes = _struct_ref_uniprot(block, entity_strands)
    notes.extend(uniprot_notes)

    observed_ids = {chain.chain_id for chain in structure.chains}
    for chain in structure.chains:
        chain.full_sequence = sequences.get(chain.chain_id)
        chain.uniprot_id = uniprot.get(chain.chain_id)
        if chain.full_sequence is not None and len(chain.sequence) > len(chain.full_sequence):
            notes.append(
                f"Chain {chain.chain_id!r} has {len(chain.sequence)} observed residues "
                f"but its deposited sequence is only {len(chain.full_sequence)} long; "
                f"the _entity_poly strand mapping looks wrong, so treat full_sequence "
                f"for this chain with suspicion."
            )
    unmapped = sorted(set(sequences) - observed_ids)
    if unmapped:
        notes.append(
            f"_entity_poly names chain(s) {unmapped[:5]} that have no atoms in the "
            f"model read; their deposited sequences were not attached to anything."
        )

    # Parse provenance first, so a reader of structure.notes sees what was read before
    # what from_atom_records inferred from it.
    structure.notes[:0] = notes
    return structure


def read_cif(
    path: str | Path,
    *,
    model: int | None = None,
    keep_hydrogens: bool = False,
) -> Structure:
    """Read an mmCIF file into a :class:`~dodo.structure.Structure`.

    Parameters
    ----------
    path
        Path to a ``.cif``/``.mmcif`` file, optionally gzipped (detected by magic
        bytes, since that is how RCSB and the AlphaFold database serve them).
    model
        ``pdbx_PDB_model_num`` to read; ``None`` reads the first model.
    keep_hydrogens
        Keep hydrogen and deuterium atoms.

    Returns
    -------
    Structure
        See :func:`parse_cif_text`.

    Raises
    ------
    StructureFileError
        If the file cannot be read or decompressed, plus everything
        :func:`parse_cif_text` raises.
    """
    file_path = Path(path)
    try:
        raw = file_path.read_bytes()
    except OSError as exc:
        raise StructureFileError(f"Could not read {file_path}: {exc}") from exc

    if raw[:2] == b"\x1f\x8b":
        try:
            raw = gzip.decompress(raw)
        except OSError as exc:
            raise StructureFileError(
                f"{file_path} starts with the gzip magic bytes but could not be decompressed: {exc}"
            ) from exc

    decode_note: str | None = None
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        # mmCIF is meant to be ASCII, but depositor-authored fields (author names,
        # titles) occasionally carry latin-1 bytes. Those fields are not coordinates,
        # so decode permissively rather than refusing the file.
        text = raw.decode("latin-1")
        decode_note = f"{file_path} is not valid UTF-8; it was decoded as latin-1."

    structure = parse_cif_text(
        text, source=str(file_path), model=model, keep_hydrogens=keep_hydrogens
    )
    if decode_note is not None:
        structure.notes.insert(0, decode_note)
    return structure
