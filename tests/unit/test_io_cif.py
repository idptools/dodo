"""Tests for the mmCIF reader.

Most of these are regression tests for the reader being replaced, which located the
``_atom_site`` loop with line-prefix heuristics. Each inline fixture is the smallest
legal file that broke it: tags out of order, quoted values, multi-line text fields,
comments inside a loop, rows wrapped across lines, and several models in one file.
"""

from __future__ import annotations

import gzip
import time
from pathlib import Path

import numpy as np
import pytest

from dodo.constants import THREE_TO_ONE
from dodo.exceptions import (
    EmptyStructureError,
    MalformedRecordError,
    StructureFileError,
    UnsupportedFormatError,
)
from dodo.io.cif import parse_cif_text, read_cif, read_unobserved_residues
from dodo.structure import Structure

DATA = Path(__file__).resolve().parents[1] / "data" / "structures"

# The canonical tag order, as AlphaFold and RCSB write it.
DEFAULT_TAGS = (
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_seq_id",
    "pdbx_PDB_ins_code",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "occupancy",
    "B_iso_or_equiv",
    "auth_seq_id",
    "auth_asym_id",
    "pdbx_PDB_model_num",
)


def make_cif(
    rows: str,
    *,
    tags: tuple[str, ...] = DEFAULT_TAGS,
    preamble: str = "",
    block: str = "test",
) -> str:
    """Assemble a minimal one-block mmCIF file around ``rows``."""
    tag_lines = "\n".join(f"_atom_site.{tag}" for tag in tags)
    return f"data_{block}\n{preamble}\nloop_\n{tag_lines}\n{rows}\n#\n"


def row(
    *,
    group: str = "ATOM",
    serial: int = 1,
    element: str = "C",
    atom: str = "CA",
    alt: str = ".",
    comp: str = "ALA",
    asym: str = "A",
    label_seq: int | str = 1,
    icode: str = "?",
    x: float = 0.0,
    y: float = 0.0,
    z: float = 0.0,
    occupancy: str = "1.00",
    b: str = "50.00",
    auth_seq: int | str = 1,
    auth_asym: str = "A",
    model: int | str = 1,
) -> str:
    """One ``_atom_site`` row in DEFAULT_TAGS order."""
    fields = [
        group,
        serial,
        element,
        atom,
        alt,
        comp,
        asym,
        label_seq,
        icode,
        x,
        y,
        z,
        occupancy,
        b,
        auth_seq,
        auth_asym,
        model,
    ]
    return " ".join(str(f) for f in fields)


def column_aligned(rows: str) -> str:
    """Re-render ``rows`` with the first field padded to six characters.

    This is how RCSB and gemmi write ``_atom_site``: the row literally begins with the
    text ``"ATOM  "``, so its first six characters are identical to a PDB ATOM record's.
    Every fixture built with :func:`row` joins fields with a single space, which is *not*
    how real files look -- and a format sniffer keyed on line prefixes passes on the
    single-spaced fixture while failing on every real file.
    """
    return "\n".join(
        f"{fields[0]:<6} {' '.join(fields[1:])}"
        for fields in (line.split() for line in rows.splitlines())
    )


def backbone(residue_index: int, *, comp: str = "ALA", auth_asym: str = "A", model: int = 1) -> str:
    """Four backbone atoms of one residue, so the structure has a usable CA."""
    lines = []
    for offset, (atom, element) in enumerate((("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"))):
        lines.append(
            row(
                serial=residue_index * 4 + offset,
                element=element,
                atom=atom,
                comp=comp,
                label_seq=residue_index,
                auth_seq=residue_index,
                auth_asym=auth_asym,
                asym=auth_asym,
                model=model,
                x=float(residue_index) * 3.8,
                y=0.0 if atom == "CA" else 0.4 * (offset + 1),
            )
        )
    return "\n".join(lines)


def small_structure(n_residues: int = 3, **kwargs: object) -> str:
    """Build a tiny all-backbone single-chain file."""
    return make_cif(
        "\n".join(backbone(i + 1, **kwargs) for i in range(n_residues))  # type: ignore[arg-type]
    )


# ---------------------------------------------------------------------------
# Tokenizer and loop parsing
# ---------------------------------------------------------------------------


class TestTokenizer:
    def test_minimal_file_reads(self) -> None:
        s = parse_cif_text(small_structure(3))
        assert isinstance(s, Structure)
        assert s.n_residues == 3
        assert s.n_atoms == 12
        assert s.sequence == "AAA"
        assert [c.chain_id for c in s.chains] == ["A"]

    def test_tags_in_any_order(self) -> None:
        """Columns resolve by tag name, never by position.

        The old reader assumed ``group_PDB`` came first and raised
        ``ValueError('Invalid mmCIF format')`` on any legal file that ordered its tags
        differently.
        """
        shuffled = tuple(reversed(DEFAULT_TAGS))
        rows = "\n".join(" ".join(reversed(line.split())) for line in backbone(1).splitlines())
        s = parse_cif_text(make_cif(rows, tags=shuffled))
        assert s.n_atoms == 4
        assert s.residue_name.tolist() == ["ALA"]
        assert s.residue_number.tolist() == [1]
        assert s.chains[0].chain_id == "A"

    def test_quoted_values_are_not_split(self) -> None:
        """A quoted value holds whitespace, and 6kn7's ``"O5'"`` must survive intact."""
        rows = "\n".join(
            [
                row(serial=1, atom='"O5\'"', comp="ALA", element="O"),
                row(serial=2, atom="CA", comp="ALA"),
                row(serial=3, atom="'C A'", comp="ALA", element="C"),
            ]
        )
        s = parse_cif_text(make_cif(rows))
        assert s.n_atoms == 3
        assert s.atom_name.tolist() == ["O5'", "CA", "C A"]

    def test_apostrophe_inside_double_quotes(self) -> None:
        # A closing quote is only a closing quote when followed by whitespace, so the
        # apostrophe in O5' does not terminate the double-quoted token.
        s = parse_cif_text(make_cif(row(atom='"O5\'"', element="O") + "\n" + row(serial=2)))
        assert "O5'" in s.atom_name.tolist()

    def test_unquoted_prime_in_atom_name(self) -> None:
        s = parse_cif_text(make_cif(row(atom="O5'", element="O") + "\n" + row(serial=2)))
        assert "O5'" in s.atom_name.tolist()

    def test_hash_inside_quotes_is_not_a_comment(self) -> None:
        s = parse_cif_text(make_cif(row(atom='"C#1"') + "\n" + row(serial=2)))
        assert "C#1" in s.atom_name.tolist()

    def test_comment_and_blank_lines_inside_a_loop(self) -> None:
        """The old reader stopped at the first blank or ``#`` line, truncating the loop."""
        rows = "\n".join(
            [
                backbone(1),
                "# a comment in the middle of the loop",
                "",
                "   ",
                backbone(2),
                "  # an indented comment",
                backbone(3),
            ]
        )
        s = parse_cif_text(make_cif(rows))
        assert s.n_residues == 3, "no residue may be lost to a comment line"
        assert s.n_atoms == 12

    def test_row_wrapped_across_lines(self) -> None:
        """A row may be split over several lines; values simply continue."""
        one = backbone(1).splitlines()
        wrapped = "\n".join(
            f"{' '.join(line.split()[:6])}\n{' '.join(line.split()[6:])}" for line in one
        )
        s = parse_cif_text(make_cif(wrapped))
        assert s.n_atoms == 4
        assert s.residue_number.tolist() == [1]

    def test_multiline_semicolon_field(self) -> None:
        preamble = (
            "_entity_poly.entity_id 1\n"
            "_entity_poly.pdbx_strand_id A\n"
            "_entity_poly.pdbx_seq_one_letter_code_can\n"
            ";AAAKKK\n"
            "EEE\n"
            ";\n"
        )
        s = parse_cif_text(make_cif(backbone(1), preamble=preamble))
        assert s.chains[0].full_sequence == "AAAKKKEEE"

    def test_multiline_semicolon_field_as_a_loop_value(self) -> None:
        """A wrapped text field is one value in the middle of a row.

        This is how 6kn7 writes ``_entity_poly``: two of its six rows carry a wrapped
        sequence and four do not, so a reader that counts values per line loses the
        alignment between tags and values for every subsequent row.
        """
        preamble = (
            "loop_\n"
            "_entity_poly.entity_id\n"
            "_entity_poly.pdbx_strand_id\n"
            "_entity_poly.pdbx_seq_one_letter_code_can\n"
            "1 A\n;AAA\nKKK\n;\n"
            "2 B GGG\n"
            "#\n"
        )
        rows = "\n".join(
            row(serial=i, atom="CA", auth_asym=c, asym=c) for i, c in enumerate("AB", start=1)
        )
        s = parse_cif_text(make_cif(rows, preamble=preamble))
        assert [(c.chain_id, c.full_sequence) for c in s.chains] == [("A", "AAAKKK"), ("B", "GGG")]

    def test_values_after_a_closing_semicolon_are_not_lost(self) -> None:
        """The closing ``;`` need only be followed by whitespace; more values may follow.

        Discarding the rest of that line drops a whole loop row with no error and no
        note, which is silent data loss in the one tokenizer path that justifies this
        module reading deposited sequences at all.
        """
        preamble = (
            "loop_\n"
            "_entity_poly.entity_id\n"
            "_entity_poly.pdbx_strand_id\n"
            "_entity_poly.pdbx_seq_one_letter_code_can\n"
            "1 A\n;AAAA\n; 2 B GGGG\n"
            "#\n"
        )
        rows = "\n".join(
            row(serial=i, atom="CA", auth_asym=c, asym=c) for i, c in enumerate("AB", start=1)
        )
        s = parse_cif_text(make_cif(rows, preamble=preamble))
        assert [(c.chain_id, c.full_sequence) for c in s.chains] == [("A", "AAAA"), ("B", "GGGG")]

    def test_comment_after_a_closing_semicolon_is_still_a_comment(self) -> None:
        preamble = (
            "_entity_poly.entity_id 1\n"
            "_entity_poly.pdbx_strand_id A\n"
            "_entity_poly.pdbx_seq_one_letter_code_can\n;AAAA\n; # and a comment\n"
        )
        s = parse_cif_text(make_cif(backbone(1), preamble=preamble))
        assert s.chains[0].full_sequence == "AAAA"

    def test_quoted_value_after_a_closing_semicolon_survives(self) -> None:
        preamble = (
            "loop_\n"
            "_entity_poly.entity_id\n"
            "_entity_poly.pdbx_strand_id\n"
            "_entity_poly.pdbx_seq_one_letter_code_can\n"
            "1 A\n;AAAA\n; 2 'B' GGGG\n"
            "#\n"
        )
        rows = "\n".join(
            row(serial=i, atom="CA", auth_asym=c, asym=c) for i, c in enumerate("AB", start=1)
        )
        s = parse_cif_text(make_cif(rows, preamble=preamble))
        assert s.chains[1].full_sequence == "GGGG"

    def test_only_the_indispensable_columns(self) -> None:
        """A loop with no group_PDB, occupancy, B-factor or model column still reads."""
        tags = (
            "label_atom_id",
            "label_comp_id",
            "auth_asym_id",
            "auth_seq_id",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
        )
        s = parse_cif_text(make_cif("CA ALA A 1 0 0 0\nCA GLY A 2 3.8 0 0", tags=tags))
        assert s.sequence == "AG"
        assert s.occupancy.tolist() == [1.0, 1.0]
        assert s.b_factor.tolist() == [0.0, 0.0]
        assert any("no group_PDB" in note for note in s.notes)

    def test_placeholders_become_defaults(self) -> None:
        rows = row(icode="?", occupancy=".", b="?")
        s = parse_cif_text(make_cif(rows))
        assert s.insertion_code.tolist() == [""]
        assert s.occupancy.tolist() == [1.0]
        assert s.b_factor.tolist() == [0.0]

    def test_quoted_dot_is_a_literal_value(self) -> None:
        # Quoting makes '.' data rather than the inapplicable placeholder.
        s = parse_cif_text(make_cif(row(icode="'A'") + "\n" + row(serial=2, icode="?")))
        assert s.insertion_code.tolist() == ["A", ""]

    def test_incomplete_final_row_raises(self) -> None:
        text = make_cif(backbone(1) + "\nATOM 5 C CA . ALA A 2")
        with pytest.raises(MalformedRecordError, match="incomplete"):
            parse_cif_text(text)

    def test_incomplete_row_reports_its_line(self) -> None:
        text = make_cif(backbone(1) + "\nATOM 5 C CA . ALA A 2")
        with pytest.raises(MalformedRecordError) as info:
            parse_cif_text(text)
        assert info.value.line_number is not None

    def test_unterminated_quote_raises(self) -> None:
        with pytest.raises(MalformedRecordError, match="Unterminated"):
            parse_cif_text(make_cif(row(atom="'CA")))

    def test_unterminated_text_field_raises(self) -> None:
        text = "data_x\n_entity_poly.pdbx_seq_one_letter_code_can\n;AAA\nBBB\n"
        with pytest.raises(MalformedRecordError, match="Unterminated semicolon"):
            parse_cif_text(text)

    def test_repeated_tag_in_loop_raises(self) -> None:
        with pytest.raises(MalformedRecordError, match="repeats tag"):
            parse_cif_text(make_cif(backbone(1), tags=(*DEFAULT_TAGS, "Cartn_x")))

    def test_loop_mixing_categories_raises(self) -> None:
        text = "data_x\nloop_\n_atom_site.id\n_other_thing.id\n1 2\n#\n"
        with pytest.raises(MalformedRecordError, match="mixes categories"):
            parse_cif_text(text)

    def test_tag_without_a_value_raises(self) -> None:
        with pytest.raises(MalformedRecordError, match="no value"):
            parse_cif_text(make_cif(backbone(1)) + "\n_entry.id\n")

    def test_duplicate_item_tag_raises(self) -> None:
        text = make_cif(
            backbone(1), preamble="_entity_poly.entity_id 1\n_entity_poly.entity_id 2\n"
        )
        with pytest.raises(MalformedRecordError, match="appears twice"):
            parse_cif_text(text)

    def test_multiple_data_blocks_are_noted(self) -> None:
        text = "data_metadata\n_entry.id XXXX\n#\n" + small_structure(1)
        s = parse_cif_text(text)
        assert any("other block" in note for note in s.notes)


# ---------------------------------------------------------------------------
# Format rejection
# ---------------------------------------------------------------------------


class TestFormatRejection:
    def test_pdb_text_is_rejected_by_name(self) -> None:
        pdb = (
            "HEADER    TEST\n"
            "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N\n"
            "END\n"
        )
        with pytest.raises(UnsupportedFormatError, match="read_pdb"):
            parse_cif_text(pdb)

    def test_text_without_a_data_block_is_rejected(self) -> None:
        with pytest.raises(UnsupportedFormatError, match="no data_ block"):
            parse_cif_text("# just a comment\n")

    def test_column_aligned_mmcif_reads(self) -> None:
        """A real RCSB-style file, whose atom rows begin with the text ``"ATOM  "``."""
        s = parse_cif_text(make_cif(column_aligned(backbone(1))))
        assert s.n_atoms == 4
        assert s.residue_number.tolist() == [1]

    def test_malformed_column_aligned_mmcif_reports_the_real_error(self) -> None:
        """An interrupted download of a real mmCIF must report *why* it is malformed.

        Column alignment must not turn a genuine ``MalformedRecordError`` -- with the line
        number that locates the trouble -- into "this is a PDB file, use read_pdb".
        """
        text = make_cif(column_aligned(backbone(1) + "\nATOM 5 C CA . ALA A 2"))
        with pytest.raises(MalformedRecordError, match="incomplete") as info:
            parse_cif_text(text)
        assert info.value.line_number is not None

    def test_unterminated_quote_in_a_column_aligned_file_is_reported(self) -> None:
        text = make_cif(column_aligned(backbone(1) + "\n" + row(serial=5, atom="'CA")))
        with pytest.raises(MalformedRecordError, match="Unterminated"):
            parse_cif_text(text)

    def test_repeated_tag_in_a_column_aligned_file_is_reported(self) -> None:
        text = make_cif(column_aligned(backbone(1)), tags=(*DEFAULT_TAGS, "Cartn_x"))
        with pytest.raises(MalformedRecordError, match="repeats tag"):
            parse_cif_text(text)

    def test_column_aligned_file_missing_a_value_is_reported(self) -> None:
        text = make_cif(column_aligned(backbone(1))) + "\n_entry.id\n"
        with pytest.raises(MalformedRecordError, match="no value"):
            parse_cif_text(text)

    def test_coordinate_only_mmcif_with_no_metadata_still_reads(self) -> None:
        """The stripped-download case: a data_ block, a loop_, and nothing else."""
        s = parse_cif_text(make_cif(column_aligned("\n".join(backbone(i + 1) for i in range(3)))))
        assert s.n_residues == 3

    def test_no_atom_site_raises(self) -> None:
        with pytest.raises(StructureFileError, match="No _atom_site"):
            parse_cif_text("data_x\n_entry.id X\n")

    def test_missing_coordinate_column_raises(self) -> None:
        tags = tuple(t for t in DEFAULT_TAGS if t != "Cartn_z")
        rows = "\n".join(
            " ".join(f for i, f in enumerate(line.split()) if i != DEFAULT_TAGS.index("Cartn_z"))
            for line in backbone(1).splitlines()
        )
        with pytest.raises(StructureFileError, match="Cartn_z"):
            parse_cif_text(make_cif(rows, tags=tags))

    def test_placeholder_coordinate_raises_rather_than_inventing_one(self) -> None:
        with pytest.raises(MalformedRecordError, match="Cartn_x"):
            parse_cif_text(make_cif(row(x=".")))  # type: ignore[arg-type]

    def test_non_numeric_coordinate_raises(self) -> None:
        with pytest.raises(MalformedRecordError, match="not a number"):
            parse_cif_text(make_cif(row(x="abc")))  # type: ignore[arg-type]

    def test_non_integer_residue_number_raises(self) -> None:
        with pytest.raises(MalformedRecordError, match="not an integer"):
            parse_cif_text(make_cif(row(auth_seq="1A")))

    def test_unknown_group_pdb_raises(self) -> None:
        with pytest.raises(MalformedRecordError, match="expected ATOM or HETATM"):
            parse_cif_text(make_cif(row(group="ANISOU")))


# ---------------------------------------------------------------------------
# auth vs label identifiers
# ---------------------------------------------------------------------------


class TestAuthorIdentifiers:
    def test_auth_wins_for_both_chain_and_residue(self) -> None:
        """Chains and residue numbers both come from ``auth_*``.

        The old reader mixed ``label_asym_id`` for chains with ``auth_seq_id`` for
        residues, so its chain ids disagreed with the same entry's PDB file while its
        residue numbers agreed. 6kn7 is a real example: author chain ``a`` is label
        chain ``AA``.
        """
        rows = "\n".join(
            row(serial=i, atom=a, element=e, asym="AA", label_seq=i, auth_seq=101, auth_asym="a")
            for i, (a, e) in enumerate((("N", "N"), ("CA", "C"), ("C", "C")), start=1)
        )
        s = parse_cif_text(make_cif(rows))
        assert [c.chain_id for c in s.chains] == ["a"]
        assert s.residue_number.tolist() == [101]

    def test_label_fallback_is_used_and_noted(self) -> None:
        tags = tuple(t for t in DEFAULT_TAGS if not t.startswith("auth_"))
        drop = {DEFAULT_TAGS.index("auth_seq_id"), DEFAULT_TAGS.index("auth_asym_id")}
        rows = "\n".join(
            " ".join(f for i, f in enumerate(line.split()) if i not in drop)
            for line in backbone(7).splitlines()
        )
        s = parse_cif_text(make_cif(rows, tags=tags))
        assert s.chains[0].chain_id == "A"
        assert s.residue_number.tolist() == [7]
        assert any("no auth_asym_id" in note for note in s.notes)
        assert any("no auth_seq_id" in note for note in s.notes)


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------


class TestModels:
    def _two_models(self) -> str:
        first = "\n".join(backbone(i + 1, model=1) for i in range(2))
        second = "\n".join(backbone(i + 1, model=2) for i in range(2))
        return make_cif(f"{first}\n{second}")

    def test_first_model_is_read_by_default(self) -> None:
        """The old reader merged models, doubling every residue."""
        s = parse_cif_text(self._two_models())
        assert s.n_residues == 2
        assert s.n_atoms == 8

    def test_model_selection(self) -> None:
        s = parse_cif_text(self._two_models(), model=2)
        assert s.n_residues == 2

    def test_multiple_models_are_noted(self) -> None:
        s = parse_cif_text(self._two_models())
        note = next(n for n in s.notes if "model(s)" in n)
        assert "2 model(s)" in note
        assert "8 atom(s) from other" in note

    def test_placeholder_model_number_is_skipped_and_accounted_for(self) -> None:
        """A skipped row is always explained, even when only one model was found.

        Otherwise the difference between records read and atoms kept has no stated
        cause, which is the failure mode this whole reader exists to avoid.
        """
        rows = backbone(1) + "\n" + row(serial=9, atom="CB", auth_seq=1, model="?")
        s = parse_cif_text(make_cif(rows))
        assert s.n_atoms == 4
        assert any("unnumbered models were skipped" in note for note in s.notes)

    def test_absent_model_raises_and_lists_what_is_there(self) -> None:
        with pytest.raises(StructureFileError, match=r"model\(s\) 1, 2"):
            parse_cif_text(self._two_models(), model=7)

    def test_model_request_without_a_model_column(self) -> None:
        tags = tuple(t for t in DEFAULT_TAGS if t != "pdbx_PDB_model_num")
        rows = "\n".join(line.rsplit(" ", 1)[0] for line in backbone(1).splitlines())
        text = make_cif(rows, tags=tags)
        assert parse_cif_text(text).n_atoms == 4
        with pytest.raises(StructureFileError, match="no pdbx_PDB_model_num"):
            parse_cif_text(text, model=3)


# ---------------------------------------------------------------------------
# Polymer filtering
# ---------------------------------------------------------------------------


class TestPolymerFiltering:
    def test_mse_hetatm_survives(self) -> None:
        """Dropping MSE fabricates a chain break where the polymer is continuous."""
        rows = "\n".join(
            [
                backbone(1),
                "\n".join(
                    row(
                        group="HETATM",
                        serial=10 + i,
                        atom=a,
                        element=e,
                        comp="MSE",
                        label_seq=2,
                        auth_seq=2,
                    )
                    for i, (a, e) in enumerate((("N", "N"), ("CA", "C"), ("C", "C")))
                ),
                backbone(3),
            ]
        )
        s = parse_cif_text(make_cif(rows))
        assert s.sequence == "AMA"
        assert s.n_residues == 3

    def test_solvent_and_ligands_are_dropped_and_counted(self) -> None:
        rows = "\n".join(
            [
                backbone(1),
                row(group="HETATM", serial=20, atom="O", element="O", comp="HOH", auth_seq=500),
                row(group="HETATM", serial=21, atom="O", element="O", comp="HOH", auth_seq=501),
                row(group="HETATM", serial=22, atom="PB", element="P", comp="ADP", auth_seq=600),
            ]
        )
        s = parse_cif_text(make_cif(rows))
        assert s.n_atoms == 4
        note = next(n for n in s.notes if "non-polymer" in n)
        assert "HOH x2" in note
        assert "ADP x1" in note

    def test_atom_records_that_are_not_amino_acids_are_kept_but_flagged(self) -> None:
        """A nucleic acid chain is not silently discarded, nor silently accepted."""
        rows = "\n".join(
            [
                backbone(1),
                row(
                    group="ATOM",
                    serial=30,
                    atom="P",
                    element="P",
                    comp="DA",
                    auth_seq=1,
                    auth_asym="B",
                    asym="B",
                ),
            ]
        )
        s = parse_cif_text(make_cif(rows))
        assert s.n_atoms == 5, "ATOM records are never dropped"
        assert any("not a standard amino acid" in note for note in s.notes)

    def test_everything_filtered_out_raises_with_an_explanation(self) -> None:
        rows = row(group="HETATM", atom="O", element="O", comp="HOH", auth_seq=1)
        with pytest.raises(EmptyStructureError, match="HOH"):
            parse_cif_text(make_cif(rows))

    def test_hydrogens_dropped_by_default(self) -> None:
        rows = backbone(1) + "\n" + row(serial=5, atom="HA", element="H")
        s = parse_cif_text(make_cif(rows))
        assert s.n_atoms == 4
        assert any("hydrogen" in note for note in s.notes)

    def test_hydrogens_kept_on_request(self) -> None:
        rows = backbone(1) + "\n" + row(serial=5, atom="HA", element="H")
        s = parse_cif_text(make_cif(rows), keep_hydrogens=True)
        assert s.n_atoms == 5

    def test_deuterium_counts_as_hydrogen(self) -> None:
        rows = backbone(1) + "\n" + row(serial=5, atom="DA", element="D")
        assert parse_cif_text(make_cif(rows)).n_atoms == 4

    def test_element_inferred_when_type_symbol_missing(self) -> None:
        tags = tuple(t for t in DEFAULT_TAGS if t != "type_symbol")
        drop = DEFAULT_TAGS.index("type_symbol")
        rows = "\n".join(
            " ".join(f for i, f in enumerate(line.split()) if i != drop)
            for line in (backbone(1) + "\n" + row(serial=5, atom="SE", comp="MSE")).splitlines()
        )
        s = parse_cif_text(make_cif(rows, tags=tags))
        # Two-letter guesses are only taken when they name a known element, so the alpha
        # carbon stays carbon instead of becoming calcium.
        assert s.element.tolist() == ["N", "C", "C", "O", "SE"]
        assert any("Inferred the element" in note for note in s.notes)


# ---------------------------------------------------------------------------
# Alternate locations and insertion codes
# ---------------------------------------------------------------------------


class TestAlternateLocations:
    def test_only_the_first_alternate_is_kept(self) -> None:
        """Two conformers of one residue would give it two CA atoms.

        ``Structure`` has no alternate-location axis, and ``ca_indices`` silently keeps
        whichever CA came last, so the reader must resolve this and say that it did.
        """
        rows = "\n".join(
            [
                row(serial=1, atom="N", element="N"),
                row(serial=2, atom="CA", alt="A"),
                row(serial=3, atom="CA", alt="B", x=9.0),
                row(serial=4, atom="C"),
            ]
        )
        s = parse_cif_text(make_cif(rows))
        assert s.n_atoms == 3
        assert s.ca_xyz.shape == (1, 3)
        assert s.ca_xyz[0][0] == pytest.approx(0.0)
        assert any("alternate locations" in note for note in s.notes)

    def test_alternate_tracking_resets_between_residues(self) -> None:
        rows = "\n".join(
            [
                row(serial=1, atom="CA", alt="A", auth_seq=1, label_seq=1),
                row(serial=2, atom="CA", alt="B", auth_seq=2, label_seq=2),
            ]
        )
        s = parse_cif_text(make_cif(rows))
        assert s.n_atoms == 2, "a per-residue choice must not leak into the next residue"

    def test_conformer_ordered_alternates_are_resolved(self) -> None:
        """Conformers need not be contiguous: a file may list all A, then all B.

        Tracking the chosen alternate per *run* of records rather than per residue
        identity keeps every alternate in such a file, which doubles the residue count
        and gives each residue two CA atoms -- the very defect the altloc filter exists
        to prevent -- while ``skipped_altloc`` stays 0 so nothing is even reported.
        """
        rows = []
        serial = 0
        for alt, shift in (("A", 0.0), ("B", 9.0)):
            for resnum in (1, 2):
                for atom, element in (("N", "N"), ("CA", "C"), ("C", "C")):
                    serial += 1
                    rows.append(
                        row(
                            serial=serial,
                            atom=atom,
                            element=element,
                            alt=alt,
                            auth_seq=resnum,
                            label_seq=resnum,
                            x=resnum * 3.8 + shift,
                        )
                    )
        s = parse_cif_text(make_cif("\n".join(rows)))
        assert s.n_residues == 2, "conformer-ordered alternates must not double the residues"
        assert s.n_atoms == 6
        assert s.ca_xyz.shape == (2, 3)
        assert any("alternate locations" in note for note in s.notes)

    def test_insertion_codes_keep_residues_distinct(self) -> None:
        rows = "\n".join(
            [
                row(serial=1, atom="CA", auth_seq=10, icode="?"),
                row(serial=2, atom="CA", auth_seq=10, icode="A", comp="GLY"),
            ]
        )
        s = parse_cif_text(make_cif(rows))
        assert s.n_residues == 2
        assert s.insertion_code.tolist() == ["", "A"]


# ---------------------------------------------------------------------------
# Fixed-width storage in Structure
# ---------------------------------------------------------------------------


class TestFixedWidthStorage:
    """Every identifier is checked against the width ``Structure`` actually stores.

    A clipped identity field merges two distinct residues or chains, which is the exact
    silent data loss this reader exists to prevent, so those raise. A clipped label is
    only lossy, so it is noted.
    """

    def test_over_wide_insertion_code_raises_rather_than_merging_residues(self) -> None:
        rows = "\n".join(
            [
                row(serial=1, atom="N", element="N", auth_seq=10, icode="A"),
                row(serial=2, atom="CA", auth_seq=10, icode="A"),
                row(serial=3, atom="N", element="N", auth_seq=10, icode="AB", x=10.0),
                row(serial=4, atom="CA", auth_seq=10, icode="AB", x=10.0),
            ]
        )
        with pytest.raises(StructureFileError, match="insertion_code"):
            parse_cif_text(make_cif(rows))

    def test_single_character_insertion_code_is_fine(self) -> None:
        rows = "\n".join(
            [
                row(serial=1, atom="CA", auth_seq=10, icode="A"),
                row(serial=2, atom="CA", auth_seq=10, icode="B", x=10.0),
            ]
        )
        s = parse_cif_text(make_cif(rows))
        assert s.n_residues == 2
        assert s.insertion_code.tolist() == ["A", "B"]

    def test_over_wide_chain_id_raises(self) -> None:
        rows = "\n".join(
            [
                row(serial=1, atom="CA", asym="AAAAA", auth_asym="AAAAA"),
                row(serial=2, atom="CA", asym="AAAAB", auth_asym="AAAAB", x=10.0),
            ]
        )
        with pytest.raises(StructureFileError, match="chain_id"):
            parse_cif_text(make_cif(rows))

    def test_over_wide_element_is_noted(self) -> None:
        s = parse_cif_text(make_cif(row(atom="CA", element="XYZQ")))
        assert s.element.tolist() == ["XY"]
        assert any("element" in note and "XYZQ" in note for note in s.notes)

    def test_over_wide_residue_name_is_noted(self) -> None:
        s = parse_cif_text(make_cif(row(atom="CA", comp="ABCD")))
        assert s.residue_name.tolist() == ["ABC"]
        assert any("residue name" in note and "ABCD" in note for note in s.notes)

    def test_over_wide_atom_name_is_noted(self) -> None:
        s = parse_cif_text(make_cif(row(atom="CA") + "\n" + row(serial=2, atom="HHHHH")))
        assert any("atom name" in note and "HHHHH" in note for note in s.notes)

    def test_two_component_ids_for_one_residue_are_reported_as_such(self) -> None:
        """One residue identity carrying two component ids is not a dtype problem.

        ``Structure`` keeps one component id per residue, taken from its first atom, so
        the second is genuinely lost -- but calling a three-character name "truncated by
        a fixed-width array" misdiagnoses it and never reports the actual loss.
        """
        rows = "\n".join(
            [
                row(serial=1, atom="N", element="N", auth_seq=5, comp="SER"),
                row(serial=2, atom="CA", auth_seq=5, comp="ALA"),
            ]
        )
        s = parse_cif_text(make_cif(rows))
        assert s.n_residues == 1
        assert s.residue_name.tolist() == ["SER"]
        assert any("component id" in note for note in s.notes)
        assert not any("truncated" in note for note in s.notes), (
            "'ALA' fits the residue_name dtype exactly; nothing was truncated"
        )


# ---------------------------------------------------------------------------
# Ordering
# ---------------------------------------------------------------------------


class TestOrdering:
    def test_file_order_is_preserved_and_anomalies_are_reported(self) -> None:
        """Records out of order are reported, not silently regrouped.

        ``from_atom_records`` groups on *change* of the identifying triple, so an
        interleaved file yields repeated residue identifiers and a note. Sorting here
        instead would hide a malformed file, and merging would lose atoms.
        """
        rows = "\n".join(
            [
                row(serial=1, atom="CA", auth_seq=1),
                row(serial=2, atom="CA", auth_seq=2),
                row(serial=3, atom="CA", auth_seq=1),
            ]
        )
        s = parse_cif_text(make_cif(rows))
        assert s.n_atoms == 3, "no atom may be dropped"
        assert s.n_residues == 3
        assert any("more than one run" in note for note in s.notes)

    def test_interleaved_chains_keep_every_atom(self) -> None:
        rows = "\n".join(
            [
                row(serial=1, atom="CA", auth_asym="A", auth_seq=1),
                row(serial=2, atom="CA", auth_asym="B", auth_seq=1),
                row(serial=3, atom="CA", auth_asym="A", auth_seq=2),
            ]
        )
        s = parse_cif_text(make_cif(rows))
        assert s.n_atoms == 3
        assert [c.chain_id for c in s.chains] == ["A", "B", "A"]


# ---------------------------------------------------------------------------
# Deposited sequence and UniProt accession
# ---------------------------------------------------------------------------


class TestPolymerMetadata:
    def test_entity_poly_as_plain_items(self) -> None:
        """AlphaFold writes a monomer's ``_entity_poly`` as items, not a loop."""
        preamble = (
            "_entity_poly.entity_id 1\n"
            "_entity_poly.pdbx_strand_id A\n"
            "_entity_poly.pdbx_seq_one_letter_code_can AAAKAAA\n"
        )
        s = parse_cif_text(make_cif(backbone(1), preamble=preamble))
        assert s.chains[0].full_sequence == "AAAKAAA"

    def test_entity_poly_loop_maps_several_strands(self) -> None:
        preamble = (
            "loop_\n"
            "_entity_poly.entity_id\n"
            "_entity_poly.pdbx_strand_id\n"
            "_entity_poly.pdbx_seq_one_letter_code_can\n"
            "1 A,B AAAA\n"
            "2 C KKKK\n"
            "#\n"
        )
        rows = "\n".join(
            row(serial=i, atom="CA", auth_asym=c, asym=c, auth_seq=1)
            for i, c in enumerate("ABC", start=1)
        )
        s = parse_cif_text(make_cif(rows, preamble=preamble))
        assert [c.full_sequence for c in s.chains] == ["AAAA", "AAAA", "KKKK"]

    def test_deposited_sequence_is_longer_than_the_observed_one(self) -> None:
        # This is the whole point of reading it: the gap is what must be rebuilt.
        preamble = (
            "_entity_poly.entity_id 1\n"
            "_entity_poly.pdbx_strand_id A\n"
            "_entity_poly.pdbx_seq_one_letter_code_can\n;AAAAAAAAAA\n;\n"
        )
        s = parse_cif_text(make_cif(backbone(1), preamble=preamble))
        assert s.chains[0].full_sequence is not None
        assert len(s.chains[0].full_sequence) > len(s.chains[0].sequence)

    def test_non_canonical_sequence_fallback_expands_modified_monomers(self) -> None:
        preamble = (
            "_entity_poly.entity_id 1\n"
            "_entity_poly.pdbx_strand_id A\n"
            "_entity_poly.pdbx_seq_one_letter_code AA(MSE)K(XYZ)\n"
        )
        s = parse_cif_text(make_cif(backbone(1), preamble=preamble))
        assert s.chains[0].full_sequence == "AAMKX"
        assert any("pdbx_seq_one_letter_code_can" in note for note in s.notes)

    def test_sequence_for_an_absent_chain_is_noted(self) -> None:
        preamble = (
            "loop_\n"
            "_entity_poly.entity_id\n"
            "_entity_poly.pdbx_strand_id\n"
            "_entity_poly.pdbx_seq_one_letter_code_can\n"
            "1 A AAAA\n"
            "2 Z KKKK\n"
            "#\n"
        )
        s = parse_cif_text(make_cif(backbone(1), preamble=preamble))
        assert any("no atoms in the model" in note for note in s.notes)

    def test_shorter_deposited_sequence_is_flagged(self) -> None:
        preamble = (
            "_entity_poly.entity_id 1\n"
            "_entity_poly.pdbx_strand_id A\n"
            "_entity_poly.pdbx_seq_one_letter_code_can A\n"
        )
        s = parse_cif_text(
            make_cif("\n".join(backbone(i + 1) for i in range(3)), preamble=preamble)
        )
        assert any("looks wrong" in note for note in s.notes)

    def test_uniprot_via_struct_ref_seq(self) -> None:
        preamble = (
            "_struct_ref.id 1\n"
            "_struct_ref.db_name UNP\n"
            "_struct_ref.pdbx_db_accession Q09472\n"
            "_struct_ref.entity_id 1\n"
            "_struct_ref_seq.align_id 1\n"
            "_struct_ref_seq.ref_id 1\n"
            "_struct_ref_seq.pdbx_strand_id A\n"
        )
        s = parse_cif_text(make_cif(backbone(1), preamble=preamble))
        assert s.chains[0].uniprot_id == "Q09472"

    def test_non_uniprot_reference_is_ignored_and_noted(self) -> None:
        preamble = (
            "_struct_ref.id 1\n"
            "_struct_ref.db_name GB\n"
            "_struct_ref.pdbx_db_accession AF12345\n"
            "_struct_ref_seq.ref_id 1\n"
            "_struct_ref_seq.pdbx_strand_id A\n"
        )
        s = parse_cif_text(make_cif(backbone(1), preamble=preamble))
        assert s.chains[0].uniprot_id is None
        assert any("non-UniProt" in note for note in s.notes)

    def test_uniprot_falls_back_to_the_entity_mapping(self) -> None:
        # Some depositions carry _struct_ref with no _struct_ref_seq at all.
        preamble = (
            "_entity_poly.entity_id 1\n"
            "_entity_poly.pdbx_strand_id A\n"
            "_entity_poly.pdbx_seq_one_letter_code_can AAAA\n"
            "_struct_ref.id 1\n"
            "_struct_ref.db_name UNP\n"
            "_struct_ref.pdbx_db_accession P12345\n"
            "_struct_ref.entity_id 1\n"
        )
        s = parse_cif_text(make_cif(backbone(1), preamble=preamble))
        assert s.chains[0].uniprot_id == "P12345"
        assert any("_entity_poly.pdbx_strand_id instead" in note for note in s.notes)

    def test_no_metadata_leaves_the_fields_none(self) -> None:
        s = parse_cif_text(small_structure(2))
        assert s.chains[0].full_sequence is None
        assert s.chains[0].uniprot_id is None


# ---------------------------------------------------------------------------
# Unobserved residues
# ---------------------------------------------------------------------------

UNOBS = """
loop_
_pdbx_unobs_or_zero_occ_residues.id
_pdbx_unobs_or_zero_occ_residues.PDB_model_num
_pdbx_unobs_or_zero_occ_residues.auth_asym_id
_pdbx_unobs_or_zero_occ_residues.auth_comp_id
_pdbx_unobs_or_zero_occ_residues.auth_seq_id
_pdbx_unobs_or_zero_occ_residues.label_asym_id
_pdbx_unobs_or_zero_occ_residues.label_seq_id
1 1 A MET 1 A 1
2 1 A GLY 2 A 2
3 1 B LYS 47 B 47
4 2 A MET 1 A 1
#
"""


#: Kabat-numbered CDR gaps: 100 and 100A are different residues of the construct.
UNOBS_WITH_INSERTION_CODES = """
loop_
_pdbx_unobs_or_zero_occ_residues.id
_pdbx_unobs_or_zero_occ_residues.auth_asym_id
_pdbx_unobs_or_zero_occ_residues.auth_comp_id
_pdbx_unobs_or_zero_occ_residues.auth_seq_id
_pdbx_unobs_or_zero_occ_residues.PDB_ins_code
1 H MET 100 ?
2 H GLY 100 A
3 H SER 101 ?
4 H ALA 0 ?
#
"""


class TestUnobservedResidues:
    def test_reads_per_chain_residue_numbers(self) -> None:
        assert read_unobserved_residues(make_cif(backbone(1), preamble=UNOBS)) == {
            "A": [(1, ""), (2, "")],
            "B": [(47, "")],
        }

    def test_insertion_coded_gaps_stay_distinct(self) -> None:
        """100 and 100A are two unmodelled residues, and both must be reported.

        Deduplicating on the bare integer collapses them, which understates what is
        missing from the model -- the one thing this function exists not to do. Kabat CDR
        numbering makes insertion-coded gaps routine in antibody structures, which are
        exactly the entries with loops to rebuild.
        """
        assert read_unobserved_residues(
            make_cif(backbone(1), preamble=UNOBS_WITH_INSERTION_CODES)
        ) == {"H": [(100, ""), (100, "A"), (101, ""), (0, "")]}

    def test_absent_category_gives_an_empty_mapping(self) -> None:
        assert read_unobserved_residues(small_structure(1)) == {}

    def test_label_fallback(self) -> None:
        text = (
            "data_x\n"
            "loop_\n"
            "_pdbx_unobs_or_zero_occ_residues.label_asym_id\n"
            "_pdbx_unobs_or_zero_occ_residues.label_seq_id\n"
            "C 12\n"
            "#\n"
        )
        assert read_unobserved_residues(text) == {"C": [(12, "")]}

    def test_missing_identifier_raises(self) -> None:
        text = (
            "data_x\n"
            "loop_\n"
            "_pdbx_unobs_or_zero_occ_residues.auth_asym_id\n"
            "_pdbx_unobs_or_zero_occ_residues.auth_seq_id\n"
            "A ?\n"
            "#\n"
        )
        with pytest.raises(MalformedRecordError, match="no usable chain id"):
            read_unobserved_residues(text)

    def test_non_integer_residue_number_raises(self) -> None:
        text = (
            "data_x\n"
            "loop_\n"
            "_pdbx_unobs_or_zero_occ_residues.auth_asym_id\n"
            "_pdbx_unobs_or_zero_occ_residues.auth_seq_id\n"
            "A 12B\n"
            "#\n"
        )
        with pytest.raises(MalformedRecordError, match="not an integer"):
            read_unobserved_residues(text)


# ---------------------------------------------------------------------------
# read_cif: files, gzip, errors
# ---------------------------------------------------------------------------


class TestReadCif:
    def test_reads_a_plain_file(self, tmp_path: Path) -> None:
        path = tmp_path / "tiny.cif"
        path.write_text(small_structure(2))
        s = read_cif(path)
        assert s.n_residues == 2
        assert s.source == str(path)

    def test_reads_a_gzipped_file(self, tmp_path: Path) -> None:
        # RCSB and the AlphaFold database both serve .cif.gz.
        path = tmp_path / "tiny.cif.gz"
        path.write_bytes(gzip.compress(small_structure(2).encode()))
        assert read_cif(path).n_residues == 2

    def test_missing_file_raises_a_dodo_error(self, tmp_path: Path) -> None:
        with pytest.raises(StructureFileError, match="Could not read"):
            read_cif(tmp_path / "nope.cif")

    def test_latin1_bytes_are_tolerated_and_noted(self, tmp_path: Path) -> None:
        path = tmp_path / "accented.cif"
        text = small_structure(1).replace("data_test", "data_test\n_audit_author.name 'Zídek'")
        path.write_bytes(text.encode("latin-1"))
        s = read_cif(path)
        assert any("latin-1" in note for note in s.notes)


# ---------------------------------------------------------------------------
# Integration against real files
# ---------------------------------------------------------------------------


class TestRealFiles:
    def test_p300(self) -> None:
        s = read_cif(DATA / "p300.cif")
        assert [c.chain_id for c in s.chains] == ["A"]
        assert s.n_residues == 2414
        assert s.n_atoms == 18457
        assert s.chains[0].uniprot_id == "Q09472"
        assert s.chains[0].full_sequence == s.sequence
        assert s.sequence.startswith("MAENVVEPGPPSAKRPKLSSPALSASASDGTDFGSLFDLE")
        # An AlphaFold model carries pLDDT in the B-factor column.
        assert float(s.b_factor.min()) >= 0.0
        assert float(s.b_factor.max()) <= 100.0

    def test_arf19(self) -> None:
        s = read_cif(DATA / "arf19.cif")
        assert [c.chain_id for c in s.chains] == ["A"]
        assert s.n_residues == 1086
        assert s.chains[0].uniprot_id == "Q8RYC8"
        assert s.chains[0].full_sequence == s.sequence

    @pytest.mark.slow
    def test_6kn7_large_assembly(self) -> None:
        start = time.perf_counter()
        s = read_cif(DATA / "6kn7.cif")
        elapsed = time.perf_counter() - start
        print(
            f"\n6kn7.cif ({(DATA / '6kn7.cif').stat().st_size / 1e6:.1f} MB) "
            f"read in {elapsed:.2f} s: {s!r}"
        )
        assert len(s.chains) == 29
        assert s.n_atoms == 61511
        # 61916 _atom_site rows minus 405 ADP atoms; every one of them accounted for.
        note = next(n for n in s.notes if "non-polymer" in n)
        assert "ADP x405" in note
        assert elapsed < 30.0, "a 6.8 MB entry must not take half a minute"

    @pytest.mark.slow
    def test_6kn7_uses_author_chain_ids(self) -> None:
        """Author chain ids, which are the ones the PDB file and the paper use.

        In 6kn7 the label ids differ from the author ids for 18 of the 29 chains
        (author ``a`` is label ``AA``), so a label-based reader disagrees with the
        equivalent PDB file about what a chain is called.
        """
        s = read_cif(DATA / "6kn7.cif")
        ids = [c.chain_id for c in s.chains]
        assert {"a", "b", "c"} <= set(ids)
        assert not {"AA", "BA", "CA"} & set(ids)
        assert ids == sorted(set(ids), key=ids.index), "chains must not repeat"

    @pytest.mark.slow
    def test_6kn7_deposited_sequences_cover_every_chain(self) -> None:
        s = read_cif(DATA / "6kn7.cif")
        assert all(c.full_sequence for c in s.chains)
        # Every chain is at least as long deposited as observed; the difference is what
        # the depositor did not model.
        for chain in s.chains:
            assert chain.full_sequence is not None
            assert len(chain.full_sequence) >= len(chain.sequence)

    @pytest.mark.slow
    def test_p300_cif_and_pdb_agree(self) -> None:
        """The same entry read as mmCIF and as PDB must give the same chains and sequence.

        Parsed here with a deliberately independent five-line PDB reader rather than
        DODO's, so this test checks the mmCIF reader against the file format itself and
        not against a sibling module that might share a mistake.
        """
        residues: list[tuple[str, int, str, str]] = []
        with open(DATA / "p300.pdb") as handle:
            for line in handle:
                if line.startswith(("ATOM  ", "HETATM")):
                    key = (line[21], int(line[22:26]), line[26].strip(), line[17:20].strip())
                    if not residues or residues[-1] != key:
                        residues.append(key)

        s = read_cif(DATA / "p300.cif")
        assert [k[0] for k in residues][:1] == [c.chain_id for c in s.chains]
        assert [k[1] for k in residues] == s.residue_number.tolist()
        assert [k[2] for k in residues] == [str(i) for i in s.insertion_code.tolist()]
        pdb_sequence = "".join(THREE_TO_ONE.get(k[3], "X") for k in residues)
        assert pdb_sequence == s.sequence

    @pytest.mark.slow
    def test_real_files_have_no_unobserved_residue_records(self) -> None:
        # None of the fixtures carry the category; assert the documented empty result
        # rather than pretending the feature is exercised by them.
        for name in ("p300.cif", "arf19.cif", "6kn7.cif"):
            assert read_unobserved_residues((DATA / name).read_text()) == {}

    @pytest.mark.slow
    @pytest.mark.parametrize("name", ["p300.cif", "arf19.cif", "6kn7.cif"])
    def test_real_files_satisfy_every_structure_invariant(self, name: str) -> None:
        s = read_cif(DATA / name)
        s.validate()
        # A CA for every residue, which every downstream geometric operation assumes.
        assert s.ca_xyz.shape == (s.n_residues, 3)
        assert np.isfinite(s.xyz).all()
