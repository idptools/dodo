"""Tests for the fixed-column PDB reader.

Most of these are regression tests for defects confirmed in the pre-rewrite reader.
Each such test is named for the failure it prevents, and the docstring says what the
old behaviour was, because "the reader silently lost atoms" is the single property this
rewrite exists to eliminate.
"""

from __future__ import annotations

import gzip
import time
from pathlib import Path

import numpy as np
import pytest

from dodo.exceptions import (
    EmptyStructureError,
    MalformedRecordError,
    StructureFileError,
    UnsupportedFormatError,
)
from dodo.io.pdb import (
    _alternate_location_mask,
    _infer_element,
    decode_hybrid36,
    parse_pdb_lines,
    read_pdb,
)
from dodo.structure import Structure

DATA = Path(__file__).resolve().parents[1] / "data" / "structures"

# ---------------------------------------------------------------------------
# Fixture construction
# ---------------------------------------------------------------------------


def atom_line(
    *,
    record: str = "ATOM",
    serial: str = "    1",
    name: str = " CA ",
    alt_loc: str = " ",
    res_name: str = "ALA",
    chain_id: str = "A",
    res_seq: str = "   1",
    icode: str = " ",
    x: float = 0.0,
    y: float = 0.0,
    z: float = 0.0,
    occupancy: str = "  1.00",
    b_factor: str = " 50.00",
    element: str = " C",
    charge: str = "  ",
) -> str:
    """Assemble one ATOM/HETATM record by column, the way the spec defines it.

    Fields are passed as raw column text wherever justification is part of what is
    being tested (atom name, serial, resSeq), so a test can put " CA " and "CA  " in
    the same file and assert they are read differently.
    """
    return (
        f"{record:<6}{serial:>5} {name:<4}{alt_loc}{res_name:>3} "
        f"{chain_id}{res_seq:>4}{icode}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occupancy}{b_factor}          {element}{charge}"
    )


def ca_trace(n: int, *, chain_id: str = "A", first: int = 1) -> list[str]:
    """Build a straight CA-only chain of ``n`` alanines, as ATOM records."""
    return [
        atom_line(
            serial=f"{i + 1:>5}",
            res_seq=f"{first + i:>4}",
            chain_id=chain_id,
            x=float(i) * 3.8,
        )
        for i in range(n)
    ]


# ---------------------------------------------------------------------------
# The column map itself
# ---------------------------------------------------------------------------


class TestFixedColumnParsing:
    def test_real_deposited_line_is_read_field_by_field(self) -> None:
        """A verbatim line from 6kn7.pdb, to pin the column offsets to reality."""
        line = "ATOM      1  N   ASP A   1     263.884 252.708  23.742  1.00 57.76           N  "
        s = parse_pdb_lines([line])
        assert s.atom_name.tolist() == ["N"]
        assert s.element.tolist() == ["N"]
        assert s.residue_name.tolist() == ["ASP"]
        assert s.residue_number.tolist() == [1]
        assert s.chains[0].chain_id == "A"
        assert s.xyz[0].tolist() == [263.884, 252.708, 23.742]
        assert s.occupancy[0] == pytest.approx(1.00)
        assert s.b_factor[0] == pytest.approx(57.76)

    def test_running_together_coordinates_are_read(self) -> None:
        """Wide coordinates leave no separator, so splitting on whitespace loses fields.

        This is why every field is read by column offset. A four-digit negative
        coordinate abuts its neighbour, and ``str.split()`` then returns eleven fields
        instead of twelve -- shifting the element, occupancy and B-factor of every atom
        in a large translated assembly.
        """
        line = atom_line(x=-100.123, y=-200.456, z=-12.345)
        assert "-100.123-200.456" in line, "fixture must actually have abutting fields"
        assert len(line.split()) != 12, "the naive reading really does lose a field here"
        s = parse_pdb_lines([line])
        assert s.xyz[0].tolist() == [-100.123, -200.456, -12.345]

    def test_four_character_atom_name_fills_the_field_and_is_not_truncated(self) -> None:
        s = parse_pdb_lines(
            [atom_line(name="HG12", element="  ", res_name="LEU")], keep_hydrogens=True
        )
        assert s.atom_name.tolist() == ["HG12"]
        assert s.element.tolist() == ["H"], "a filled name field is not a two-letter element"

    def test_records_other_than_atoms_are_ignored(self) -> None:
        lines = [
            "HEADER                                            01-JUN-22",
            "ANISOU    1  N   ASP A   1     7986   7986   7986      0      0      0",
            "CONECT61541614426154361544",
            *ca_trace(2),
            "END",
        ]
        assert parse_pdb_lines(lines).n_atoms == 2


# ---------------------------------------------------------------------------
# Defect 1: HETATM polymer residues were dropped
# ---------------------------------------------------------------------------


class TestPolymerHetatm:
    def test_mid_chain_selenomethionine_is_kept(self) -> None:
        """A HETATM MSE in the middle of a chain must not vanish.

        The pre-rewrite reader selected only lines starting with 'ATOM', so a
        selenomethionine -- which is deposited as HETATM -- disappeared and fabricated a
        phantom chain break where the polymer is in fact continuous.
        """
        lines = [
            *ca_trace(2),
            atom_line(record="HETATM", serial="    3", res_name="MSE", res_seq="   3", x=7.6),
            *[atom_line(serial=f"{i:>5}", res_seq=f"{i:>4}", x=3.8 * (i - 1)) for i in (4, 5)],
        ]
        s = parse_pdb_lines(lines)
        assert s.n_residues == 5
        assert s.sequence == "AAMAA", "MSE reads as methionine, and the chain is unbroken"
        assert len(s.chains) == 1

    def test_solvent_and_ligand_hetatm_are_skipped_and_recorded(self) -> None:
        lines = [
            *ca_trace(2),
            atom_line(record="HETATM", serial="    3", res_name="ADP", res_seq=" 401", name=" PB "),
            atom_line(record="HETATM", serial="    4", res_name="HOH", res_seq=" 501", name=" O  "),
        ]
        s = parse_pdb_lines(lines)
        assert s.n_atoms == 2
        assert any("ADP x1" in note and "HOH x1" in note for note in s.notes)

    def test_water_written_as_an_atom_record_is_still_solvent(self) -> None:
        """Some modelling tools emit waters as ATOM; they are not polymer residues."""
        water = atom_line(serial="    3", res_name="HOH", res_seq=" 501", name=" O  ")
        s = parse_pdb_lines([*ca_trace(2), water])
        assert s.n_atoms == 2

    def test_unrecognized_atom_residue_is_kept_and_noted(self) -> None:
        """An ATOM record is polymer by definition; dropping it would break the chain."""
        lines = [
            *ca_trace(2),
            atom_line(serial="    3", res_name="XYZ", res_seq="   3", x=7.6),
        ]
        s = parse_pdb_lines(lines)
        assert s.n_residues == 3
        assert s.sequence == "AAX"
        assert any("unrecognized residue name" in note for note in s.notes)


# ---------------------------------------------------------------------------
# Defect 2: insertion codes were ignored
# ---------------------------------------------------------------------------


class TestInsertionCodes:
    def test_residue_10_and_10a_stay_distinct(self) -> None:
        """The pre-rewrite reader keyed on int(resSeq) and merged 10 with 10A."""
        lines = [
            atom_line(serial="    1", res_seq="  10", x=0.0),
            atom_line(serial="    2", res_seq="  10", icode="A", res_name="GLY", x=3.8),
            atom_line(serial="    3", res_seq="  11", x=7.6),
        ]
        s = parse_pdb_lines(lines)
        assert s.n_residues == 3
        assert s.residue_number.tolist() == [10, 10, 11]
        assert s.insertion_code.tolist() == ["", "A", ""]
        assert s.residue_label(1) == "A:GLY10A"


# ---------------------------------------------------------------------------
# Defect 3: altLoc was ignored
# ---------------------------------------------------------------------------


class TestAlternateConformers:
    def test_alternates_do_not_duplicate_atoms(self) -> None:
        """The pre-rewrite reader produced a residue holding two CAs.

        Two conformers of the same atom are one atom, and the 'A' conformer is the one
        the deposition convention lists first / at higher occupancy.
        """
        lines = [
            atom_line(serial="    1", name=" N  ", element=" N", x=0.0),
            atom_line(serial="    2", alt_loc="A", occupancy="  0.65", x=1.0),
            atom_line(serial="    3", alt_loc="B", occupancy="  0.35", x=2.0),
            atom_line(serial="    4", res_seq="   2", x=3.8),
        ]
        s = parse_pdb_lines(lines)
        assert s.n_atoms == 3
        assert s.n_residues == 2
        assert s.xyz[1][0] == pytest.approx(1.0), "altLoc A is the conformer kept"
        assert any("alternate-conformer" in note for note in s.notes)

    def test_residue_modelled_only_as_alt_b_is_not_lost(self) -> None:
        """Keeping strictly 'A' would delete such a residue and break the chain.

        Real depositions contain residues whose only conformer is labelled 'B'. That is
        exactly the situation where a naive filter reintroduces the phantom chain break
        of defect 1, so the preferred conformer is per residue, not global.
        """
        lines = [
            atom_line(serial="    1", res_seq="   1", x=0.0),
            atom_line(serial="    2", res_seq="   2", alt_loc="B", x=3.8),
            atom_line(serial="    3", res_seq="   2", alt_loc="C", x=3.9),
            atom_line(serial="    4", res_seq="   3", x=7.6),
        ]
        s = parse_pdb_lines(lines)
        assert s.n_residues == 3
        assert s.n_atoms == 3
        assert s.xyz[1][0] == pytest.approx(3.8), "the first-listed conformer wins"

    def test_mask_helper_prefers_a_over_a_lower_letter(self) -> None:
        mask = _alternate_location_mask(
            segment_labels=["A", "A", "A"],
            residue_numbers=[1, 1, 1],
            insertion_codes=["", "", ""],
            alt_locs=["B", "A", ""],
        )
        assert mask.tolist() == [False, True, True]


# ---------------------------------------------------------------------------
# Defect 4: MODEL/ENDMDL was mishandled
# ---------------------------------------------------------------------------


def two_model_lines() -> list[str]:
    """Two models with distinguishable coordinates, plus a left-justified MODEL serial.

    "MODEL 1" with the serial not in columns 10-14 is what DODO's own v1 output looked
    like, so models are counted positionally rather than read from that field.
    """
    return [
        "MODEL 1",
        *ca_trace(3),
        "ENDMDL",
        "MODEL        2",
        *[atom_line(serial=f"{i + 1:>5}", res_seq=f"{i + 1:>4}", x=100.0 + i) for i in range(3)],
        "ENDMDL",
        "END",
    ]


def frame(x: float, *, n: int = 3) -> list[str]:
    """One frame of ``n`` alanines, tagged by a shared x so frames are distinguishable.

    Every frame carries the same residue numbers, which is what makes a merge visible:
    concatenating two frames merges their residues as well as their atoms.
    """
    return [
        atom_line(serial=f"{i + 1:>5}", res_seq=f"{i + 1:>4}", x=x, y=float(i) * 3.8)
        for i in range(n)
    ]


def implicit_then_headed_frames() -> list[str]:
    """Three frames where the first has no ``MODEL`` header but is closed by ``ENDMDL``.

    Written by hand and by scripts that emit a first frame before deciding the output
    is multi-model. The header serials are therefore off by one from the positional
    frame numbers, which is exactly why frames are counted positionally.
    """
    return [
        *frame(1.0),
        "ENDMDL",
        "MODEL        2",
        *frame(2.0),
        "ENDMDL",
        "MODEL        3",
        *frame(3.0),
        "ENDMDL",
        "END",
    ]


class TestModels:
    def test_first_model_only_by_default(self) -> None:
        """Merging models produced an impossible structure with every residue twice."""
        s = parse_pdb_lines(two_model_lines())
        assert s.n_residues == 3
        assert s.xyz[0][0] == pytest.approx(0.0)
        assert any("Read model 1 of 2" in note for note in s.notes)

    def test_second_model_can_be_selected(self) -> None:
        s = parse_pdb_lines(two_model_lines(), model=2)
        assert s.n_residues == 3
        assert s.xyz[0][0] == pytest.approx(100.0)

    def test_missing_model_reports_how_many_exist(self) -> None:
        with pytest.raises(StructureFileError, match="only 2 model"):
            parse_pdb_lines(two_model_lines(), model=5)

    def test_model_numbers_start_at_one(self) -> None:
        with pytest.raises(StructureFileError, match="numbered from 1"):
            parse_pdb_lines(two_model_lines(), model=0)

    def test_single_model_file_needs_no_model_record(self) -> None:
        s = parse_pdb_lines(ca_trace(3))
        assert s.n_residues == 3
        assert not any("Read model" in note for note in s.notes)

    def test_records_before_the_first_model_are_a_frame_of_their_own(self) -> None:
        """A ``MODEL`` record starts a frame even mid-file, so nothing is merged.

        This test previously asserted the opposite -- that the pre-``MODEL`` block and
        the block headed by ``MODEL 1`` are one model of three atoms. That assertion
        was the defect: it made the reader's frame counter off by one for every file
        with coordinates before its first ``MODEL``, so two frames were concatenated
        into a single impossible residue set, ``model=n`` returned frame ``n + 1``, and
        the last frame was unreachable. Merging distinct frames is the failure this
        whole reader exists to prevent, so the frame boundary wins.
        """
        lines = [*ca_trace(1), "MODEL        1", *ca_trace(2), "ENDMDL"]
        s = parse_pdb_lines(lines)
        assert s.n_atoms == 1
        assert any("Read model 1 of 2" in note for note in s.notes)
        assert parse_pdb_lines(lines, model=2).n_atoms == 2

    def test_implicit_first_frame_is_counted_as_a_model(self) -> None:
        """Atoms, ENDMDL, MODEL 2..., MODEL 3... is three frames, not two.

        The implicit frame was collected as model 1 but never counted, so the first
        explicit ``MODEL`` re-opened collection and appended frame 2 to it.
        """
        s = parse_pdb_lines(implicit_then_headed_frames())
        assert s.n_atoms == 3
        assert s.n_residues == 3
        assert s.xyz[:, 0].tolist() == [1.0, 1.0, 1.0]
        assert any("Read model 1 of 3" in note for note in s.notes)

    def test_model_selection_after_an_implicit_frame_is_not_off_by_one(self) -> None:
        lines = implicit_then_headed_frames()
        assert parse_pdb_lines(lines, model=2).xyz[:, 0].tolist() == [2.0, 2.0, 2.0]
        assert parse_pdb_lines(lines, model=3).xyz[:, 0].tolist() == [3.0, 3.0, 3.0]
        with pytest.raises(StructureFileError, match="only 3 model"):
            parse_pdb_lines(lines, model=4)

    def test_bare_endmdl_records_separate_frames(self) -> None:
        """Frames separated by ``ENDMDL`` alone, with no ``MODEL`` headers at all.

        Every atom after the first ``ENDMDL`` used to be dropped with an empty notes
        list, because the "Read model" note was gated on having seen a MODEL record.
        """
        lines = [*frame(1.0), "ENDMDL", *frame(2.0), "ENDMDL", *frame(3.0), "END"]
        s = parse_pdb_lines(lines)
        assert s.n_atoms == 3
        assert s.xyz[:, 0].tolist() == [1.0, 1.0, 1.0]
        note = next(note for note in s.notes if "Read model" in note)
        assert "Read model 1 of 3" in note
        assert "6 atom record(s)" in note
        assert parse_pdb_lines(lines, model=3).xyz[:, 0].tolist() == [3.0, 3.0, 3.0]

    def test_atoms_after_the_only_endmdl_are_not_lost_silently(self) -> None:
        """One MODEL/ENDMDL pair followed by more atoms: half the file used to vanish."""
        lines = ["MODEL        1", *frame(1.0), "ENDMDL", *frame(2.0), "END"]
        s = parse_pdb_lines(lines)
        assert s.n_atoms == 3
        assert s.xyz[:, 0].tolist() == [1.0, 1.0, 1.0]
        note = next(note for note in s.notes if "Read model" in note)
        assert "Read model 1 of 2" in note
        assert "3 atom record(s)" in note
        assert parse_pdb_lines(lines, model=2).xyz[:, 0].tolist() == [2.0, 2.0, 2.0]

    def test_ignored_atom_records_are_counted_in_the_note(self) -> None:
        s = parse_pdb_lines(two_model_lines())
        assert any("3 atom record(s)" in note for note in s.notes)

    def test_model_without_a_closing_endmdl_still_starts_a_new_frame(self) -> None:
        lines = ["MODEL        1", *frame(1.0), "MODEL        2", *frame(2.0), "END"]
        s = parse_pdb_lines(lines)
        assert s.n_atoms == 3
        assert s.xyz[:, 0].tolist() == [1.0, 1.0, 1.0]
        assert any("Read model 1 of 2" in note for note in s.notes)
        assert any("ENDMDL" in note for note in s.notes), "the inference must be recorded"

    def test_repeated_endmdl_is_recorded_rather_than_ignored(self) -> None:
        lines = ["MODEL        1", *frame(1.0), "ENDMDL", "ENDMDL", "END"]
        s = parse_pdb_lines(lines)
        assert s.n_atoms == 3
        assert not any("Read model" in note for note in s.notes), "there is only one frame"
        assert any("ENDMDL" in note for note in s.notes)

    def test_every_atom_record_is_either_kept_or_counted_in_a_note(self) -> None:
        """The accounting invariant: nothing leaves this reader unmentioned."""
        lines = [
            *frame(1.0),
            "ENDMDL",
            "MODEL        2",
            *frame(2.0, n=5),
            atom_line(record="HETATM", res_name="HOH", name=" O  ", element=" O"),
            "ENDMDL",
        ]
        n_records = sum(1 for line in lines if line.startswith(("ATOM", "HETATM")))
        s = parse_pdb_lines(lines)
        note = next(note for note in s.notes if "Read model" in note)
        assert f"{n_records - s.n_atoms} atom record(s)" in note


# ---------------------------------------------------------------------------
# Defect 5: hybrid-36
# ---------------------------------------------------------------------------


class TestHybrid36:
    @pytest.mark.parametrize(
        ("field", "expected"),
        [
            ("    0", 0),
            ("   42", 42),
            ("99999", 99999),
            ("A0000", 100000),
            ("A0001", 100001),
            ("A0009", 100009),
            ("A000A", 100010),
            ("A000Z", 100035),
            ("A001Z", 100071),
            ("ZZZZZ", 43770015),
            ("a0000", 43770016),
            ("zzzzz", 87440031),
        ],
    )
    def test_five_column_serial_field(self, field: str, expected: int) -> None:
        assert decode_hybrid36(field, 5) == expected

    @pytest.mark.parametrize(
        ("field", "expected"),
        [
            ("   1", 1),
            (" -12", -12),
            ("9999", 9999),
            ("A000", 10000),
            ("A001", 10001),
            ("A00A", 10010),
            ("ZZZZ", 1223055),
            ("a000", 1223056),
            ("zzzz", 2436111),
        ],
    )
    def test_four_column_resseq_field(self, field: str, expected: int) -> None:
        assert decode_hybrid36(field, 4) == expected

    def test_the_two_widths_disagree_and_must_not_be_conflated(self) -> None:
        """Same text, different value: the boundary depends on the field width."""
        assert decode_hybrid36("A000", 4) == 10000
        assert decode_hybrid36("A0000", 5) == 100000

    @pytest.mark.parametrize("field", ["", "   ", "A00", "A00000", "Aa000", "**", "12x4"])
    def test_unreadable_fields_raise(self, field: str) -> None:
        with pytest.raises(ValueError, match=r".+"):
            decode_hybrid36(field, 4)

    def test_serial_beyond_99999_does_not_crash_the_reader(self) -> None:
        """``int(line[6:11])`` raised ValueError here, which is every EM assembly.

        DODO's target inputs are cryo-EM assemblies, and anything over 99,999 atoms
        writes its serials in hybrid-36.
        """
        lines = [
            atom_line(serial="A0000", res_seq="   1", x=0.0),
            atom_line(serial="A0001", res_seq="   2", x=3.8),
        ]
        s = parse_pdb_lines(lines)
        assert s.n_atoms == 2
        assert not any("serial" in note for note in s.notes)

    def test_residue_number_beyond_9999_is_decoded(self) -> None:
        lines = [
            atom_line(serial="    1", res_seq="9999", x=0.0),
            atom_line(serial="    2", res_seq="A000", x=3.8),
            atom_line(serial="    3", res_seq="A001", x=7.6),
        ]
        s = parse_pdb_lines(lines)
        assert s.residue_number.tolist() == [9999, 10000, 10001]

    def test_negative_residue_numbers_survive(self) -> None:
        """Expression-tag residues are routinely numbered from a negative index."""
        s = parse_pdb_lines(
            [
                atom_line(serial="    1", res_seq="  -2", x=0.0),
                atom_line(serial="    2", res_seq="  -1", x=3.8),
            ]
        )
        assert s.residue_number.tolist() == [-2, -1]

    def test_unreadable_serial_is_noted_not_fatal(self) -> None:
        """Nothing downstream uses the serial, so a bad one must not lose the atom."""
        s = parse_pdb_lines([atom_line(serial="  ?? ", res_seq="   1")])
        assert s.n_atoms == 1
        assert any("serial" in note for note in s.notes)

    def test_residue_number_that_cannot_be_decoded_is_fatal(self) -> None:
        with pytest.raises(MalformedRecordError, match="residue sequence number"):
            parse_pdb_lines([atom_line(res_seq=" ?? ")])


# ---------------------------------------------------------------------------
# Defect 6: only .pdb was accepted
# ---------------------------------------------------------------------------


class TestFileNames:
    @pytest.mark.parametrize("name", ["s.pdb", "s.PDB", "s.ent", "s.pdb1", "s.pdb2", "s.Ent"])
    def test_accepted_extensions(self, tmp_path: Path, name: str) -> None:
        path = tmp_path / name
        path.write_text("\n".join(ca_trace(3)) + "\n")
        assert read_pdb(path).n_residues == 3

    @pytest.mark.parametrize("name", ["s.pdb.gz", "s.ent.gz", "s.pdb1.gz"])
    def test_gzip_is_transparent(self, tmp_path: Path, name: str) -> None:
        path = tmp_path / name
        path.write_bytes(gzip.compress(("\n".join(ca_trace(3)) + "\n").encode()))
        assert read_pdb(path).n_residues == 3

    def test_gzip_detected_from_magic_bytes_when_renamed(self, tmp_path: Path) -> None:
        """A compressed file that lost its .gz suffix still reads."""
        path = tmp_path / "renamed.pdb"
        path.write_bytes(gzip.compress(("\n".join(ca_trace(3)) + "\n").encode()))
        assert read_pdb(path).n_residues == 3

    def test_mmcif_is_rejected_with_a_pointer_to_the_right_reader(self, tmp_path: Path) -> None:
        path = tmp_path / "s.cif"
        path.write_text("data_XXXX\n")
        with pytest.raises(UnsupportedFormatError, match="mmCIF"):
            read_pdb(path)

    def test_missing_file_raises_a_dodo_error(self, tmp_path: Path) -> None:
        with pytest.raises(StructureFileError, match="Could not read"):
            read_pdb(tmp_path / "absent.pdb")

    def test_corrupt_gzip_raises_a_dodo_error(self, tmp_path: Path) -> None:
        path = tmp_path / "broken.pdb.gz"
        path.write_bytes(b"\x1f\x8bnot really gzip")
        with pytest.raises(StructureFileError, match="could not be decompressed"):
            read_pdb(path)

    def test_non_ascii_bytes_do_not_shift_columns(self, tmp_path: Path) -> None:
        """A latin-1 author name in a TITLE record must not move the coordinate columns.

        Decoding as utf-8 would collapse two bytes into one character on that line and
        (with errors="replace") could do so on any line; latin-1 keeps the byte-for-
        character mapping a fixed-column format depends on.
        """
        path = tmp_path / "s.pdb"
        body = "TITLE     STRUCTURE OF THE J\xc3\x98RGENSEN COMPLEX\n" + "\n".join(ca_trace(3))
        path.write_bytes(body.encode("latin-1"))
        assert read_pdb(path).n_residues == 3


# ---------------------------------------------------------------------------
# Defect 7: truncated lines crashed
# ---------------------------------------------------------------------------


class TestTruncatedLines:
    def test_line_ending_at_column_54_uses_default_occupancy_and_b_factor(self) -> None:
        line = atom_line(x=1.0, y=2.0, z=3.0)[:54]
        assert len(line) == 54
        s = parse_pdb_lines([line])
        assert s.xyz[0].tolist() == [1.0, 2.0, 3.0]
        assert s.occupancy[0] == pytest.approx(1.0)
        assert s.b_factor[0] == pytest.approx(0.0)

    def test_line_ending_at_column_66_gets_an_inferred_element(self) -> None:
        s = parse_pdb_lines([atom_line()[:66]])
        assert s.element.tolist() == ["C"]
        assert any("Inferred the element" in note for note in s.notes)

    def test_line_too_short_for_coordinates_names_the_line(self) -> None:
        short = atom_line()[:40]
        with pytest.raises(MalformedRecordError) as info:
            parse_pdb_lines([short])
        message = str(info.value)
        assert "40 columns" in message
        assert short in message
        assert info.value.line_number == 1

    def test_blank_coordinate_field_is_a_malformed_record(self) -> None:
        line = atom_line()
        blanked = line[:30] + " " * 8 + line[38:]
        with pytest.raises(MalformedRecordError, match="Unreadable coordinates"):
            parse_pdb_lines([blanked])

    def test_unreadable_occupancy_falls_back_and_is_noted(self) -> None:
        line = atom_line(occupancy="  N/A ")
        s = parse_pdb_lines([line])
        assert s.occupancy[0] == pytest.approx(1.0)
        assert any("unreadable" in note for note in s.notes)


# ---------------------------------------------------------------------------
# Defect 8: missing element column
# ---------------------------------------------------------------------------


class TestElementInference:
    @pytest.mark.parametrize(
        ("field", "expected"),
        [
            (" CA ", "C"),
            (" N  ", "N"),
            (" OD1", "O"),
            ("SE  ", "SE"),  # selenomethionine
            ("HG12", "H"),  # four-character hydrogen name
            ("1HG2", "H"),  # old-style, digit first
            ("NE2 ", "N"),  # hand-edited, left-justified
            ("    ", ""),
        ],
    )
    def test_inference_rule(self, field: str, expected: str) -> None:
        assert _infer_element(field) == expected

    def test_selenium_of_a_selenomethionine_is_recovered(self) -> None:
        """MSE's selenium must not become carbon; ATOMIC_MASSES differs by 6.5x."""
        line = atom_line(record="HETATM", name="SE  ", res_name="MSE", element="  ")[:66]
        s = parse_pdb_lines([line])
        assert s.element.tolist() == ["SE"]

    def test_declared_element_column_wins_over_inference(self) -> None:
        s = parse_pdb_lines([atom_line(name="CA  ", res_name="ALA", element=" C")])
        assert s.element.tolist() == ["C"]


# ---------------------------------------------------------------------------
# Defect 9: TER was not used to separate chains
# ---------------------------------------------------------------------------


class TestTerRecords:
    def test_two_chains_sharing_an_id_are_not_merged(self) -> None:
        """The pre-rewrite reader ignored TER and merged them into one chain.

        Both runs are labelled 'A' and both are numbered from 1, so without TER the
        residues merge as well and half the atoms end up inside the other half's
        residues.
        """
        lines = [*ca_trace(3), "TER    4      ALA A   3", *ca_trace(3), "END"]
        s = parse_pdb_lines(lines)
        assert s.n_atoms == 6
        assert s.n_residues == 6
        assert [c.chain_id for c in s.chains] == ["A", "A"]
        assert [len(c) for c in s.chains] == [3, 3]
        assert any("TER record" in note for note in s.notes)

    def test_three_consecutive_runs_of_one_id_stay_separate(self) -> None:
        lines = [
            *ca_trace(1),
            "TER",
            *ca_trace(1),
            "TER",
            *ca_trace(1),
        ]
        s = parse_pdb_lines(lines)
        assert [c.chain_id for c in s.chains] == ["A", "A", "A"]

    @pytest.mark.parametrize("runs", [2, 3, 4, 6])
    def test_every_ter_split_is_counted_in_the_note(self, runs: int) -> None:
        """The note must state how many splits happened, not half of them.

        The internal marker that keeps two runs of one id apart alternates on and off,
        and the count was incremented only on the "on" transitions, so every second
        split went unreported by a module whose contract is that every inferential act
        is recorded accurately.
        """
        lines: list[str] = []
        for _ in range(runs):
            lines.extend(ca_trace(1))
            lines.append("TER")
        s = parse_pdb_lines(lines)
        assert [c.chain_id for c in s.chains] == ["A"] * runs
        note = next(note for note in s.notes if "TER record" in note)
        assert f"{runs - 1} TER record(s)" in note

    def test_no_internal_marker_leaks_into_the_chain_id(self) -> None:
        s = parse_pdb_lines([*ca_trace(2), "TER", *ca_trace(2)])
        for chain in s.chains:
            assert chain.chain_id == "A"
            assert chain.chain_id.isalnum()

    def test_ter_between_different_chain_ids_changes_nothing(self) -> None:
        lines = [*ca_trace(2, chain_id="A"), "TER", *ca_trace(2, chain_id="B")]
        s = parse_pdb_lines(lines)
        assert [c.chain_id for c in s.chains] == ["A", "B"]
        assert not any("TER record" in note for note in s.notes)


# ---------------------------------------------------------------------------
# Hydrogens
# ---------------------------------------------------------------------------


class TestHydrogens:
    def hydrogen_lines(self) -> list[str]:
        return [
            atom_line(serial="    1", name=" N  ", element=" N"),
            atom_line(serial="    2", name=" CA ", element=" C", x=1.0),
            atom_line(serial="    3", name=" H  ", element=" H", x=1.5),
            atom_line(serial="    4", name=" D  ", element=" D", x=1.6),
        ]

    def test_hydrogens_dropped_by_default_and_counted(self) -> None:
        s = parse_pdb_lines(self.hydrogen_lines())
        assert s.n_atoms == 2
        assert set(s.element.tolist()) == {"N", "C"}
        assert any("2 hydrogen/deuterium atom(s) dropped" in note for note in s.notes)

    def test_hydrogens_kept_on_request_and_counted(self) -> None:
        s = parse_pdb_lines(self.hydrogen_lines(), keep_hydrogens=True)
        assert s.n_atoms == 4
        assert any("2 hydrogen/deuterium atom(s) kept" in note for note in s.notes)

    def test_inferred_hydrogen_is_also_dropped(self) -> None:
        """A file with no element column still must not smuggle hydrogens in."""
        lines = [
            atom_line(serial="    1", name=" CA ")[:66],
            atom_line(serial="    2", name="HG12", res_name="LEU")[:66],
        ]
        s = parse_pdb_lines(lines)
        assert s.n_atoms == 1


# ---------------------------------------------------------------------------
# DBREF and SEQRES harvesting
# ---------------------------------------------------------------------------


class TestMetadataRecords:
    def test_dbref_populates_uniprot_id(self) -> None:
        lines = [
            "DBREF  XXXX A    1  2414  UNP    Q09472   EP300_HUMAN      1   2414",
            *ca_trace(2),
        ]
        assert parse_pdb_lines(lines).chains[0].uniprot_id == "Q09472"

    def test_non_uniprot_dbref_is_ignored(self) -> None:
        lines = [
            "DBREF  1ABC A    1   100  GB     M12345   AAA12345         1    100",
            *ca_trace(2),
        ]
        assert parse_pdb_lines(lines).chains[0].uniprot_id is None

    def test_dbref1_dbref2_pair_is_supported(self) -> None:
        """Accessions too long for the DBREF layout are split across two records."""
        lines = [
            "DBREF1 1ABC A    1   246  UNP                  A0A123B456_ORYSI",
            "DBREF2 1ABC A     A0A123B456                        1         246",
            *ca_trace(2),
        ]
        assert parse_pdb_lines(lines).chains[0].uniprot_id == "A0A123B456"

    def test_seqres_populates_the_deposited_construct_sequence(self) -> None:
        """full_sequence is the deposited construct, not the UniProt isoform.

        Missing residues must be measured against what was in the crystallised
        construct; using the canonical isoform invents phantom terminal residues for
        every tagged or truncated construct, which is what makes this distinction load
        bearing rather than pedantic.
        """
        lines = [
            "DBREF  1ABC A    1     5  UNP    P12345   TEST_HUMAN      10     14",
            "SEQRES   1 A    5  MET ALA GLY SER MSE",
            *ca_trace(2),
        ]
        chain = parse_pdb_lines(lines).chains[0]
        assert chain.full_sequence == "MAGSM"
        assert chain.uniprot_id == "P12345"
        assert chain.sequence == "AA", "observed sequence is shorter than the construct"

    def test_seqres_spanning_several_lines_is_concatenated(self) -> None:
        lines = [
            "SEQRES   1 A   15  MET ALA GLY SER GLU LYS LEU ILE VAL PHE TRP TYR THR",
            "SEQRES   2 A   15  PRO CYS",
            *ca_trace(2),
        ]
        assert parse_pdb_lines(lines).chains[0].full_sequence == "MAGSEKLIVFWYTPC"

    def test_seqres_count_mismatch_is_noted(self) -> None:
        lines = ["SEQRES   1 A   99  MET ALA GLY", *ca_trace(2)]
        s = parse_pdb_lines(lines)
        assert any("declares 99 residues but lists 3" in note for note in s.notes)

    def test_seqres_for_an_unmodelled_chain_is_noted(self) -> None:
        lines = ["SEQRES   1 B    3  MET ALA GLY", *ca_trace(2)]
        s = parse_pdb_lines(lines)
        assert any("no coordinates" in note for note in s.notes)

    def test_metadata_is_shared_by_ter_split_chains_of_one_id(self) -> None:
        lines = [
            "DBREF  1ABC A    1     3  UNP    P12345   TEST_HUMAN       1      3",
            *ca_trace(2),
            "TER",
            *ca_trace(2),
        ]
        s = parse_pdb_lines(lines)
        assert [c.uniprot_id for c in s.chains] == ["P12345", "P12345"]


# ---------------------------------------------------------------------------
# Failure modes
# ---------------------------------------------------------------------------


class TestEmptyInput:
    def test_no_atoms_at_all(self) -> None:
        with pytest.raises(EmptyStructureError, match="No polymer atoms"):
            parse_pdb_lines(["HEADER", "END"])

    def test_only_solvent_says_what_was_skipped(self) -> None:
        lines = [atom_line(record="HETATM", res_name="HOH", name=" O  ", element=" O")]
        with pytest.raises(EmptyStructureError, match="HOH x1"):
            parse_pdb_lines(lines)

    def test_only_hydrogens_says_so(self) -> None:
        lines = [atom_line(name=" H1 ", element=" H")]
        with pytest.raises(EmptyStructureError, match="hydrogen"):
            parse_pdb_lines(lines)


# ---------------------------------------------------------------------------
# Integration: real files
# ---------------------------------------------------------------------------


class TestRealFiles:
    def test_p300_alphafold_model(self) -> None:
        s = read_pdb(DATA / "p300.pdb")
        assert s.n_residues == 2414
        assert s.n_atoms == 18457
        assert len(s.chains) == 1
        chain = s.chains[0]
        assert chain.chain_id == "A"
        assert chain.uniprot_id == "Q09472"
        assert chain.full_sequence is not None
        assert len(chain.full_sequence) == 2414
        assert chain.sequence == chain.full_sequence, "an AF2 model has no missing residues"
        # AlphaFold puts pLDDT in the B-factor column; DODO reads it as a disorder signal.
        assert 0.0 < float(s.b_factor.max()) <= 100.0
        s.validate()

    def test_dnmt3a_alphafold_model(self) -> None:
        s = read_pdb(DATA / "dnmt3a.pdb")
        assert s.n_residues == 912
        assert s.n_atoms == 7146
        assert s.chains[0].uniprot_id == "Q9Y6K1"
        assert s.chains[0].full_sequence is not None
        assert len(s.chains[0].full_sequence) == 912
        assert s.ca_xyz.shape == (912, 3)
        s.validate()

    def test_6kn7_em_assembly_loses_no_atoms(self) -> None:
        """A 29-chain cryo-EM assembly: the file this reader exists for."""
        path = DATA / "6kn7.pdb"
        expected_atoms = sum(1 for line in path.read_text().splitlines() if line.startswith("ATOM"))
        started = time.perf_counter()
        s = read_pdb(path)
        elapsed = time.perf_counter() - started

        assert s.n_atoms == expected_atoms == 61511, "every ATOM record must survive"
        assert len(s.chains) == 29
        assert [c.chain_id for c in s.chains] == list("ABCDEFGHIJKLMNOPQRSTUVWXYZabc")
        assert all(c.uniprot_id is not None for c in s.chains)
        assert all(c.full_sequence is not None for c in s.chains)
        # The 405 ADP atoms are ligand, correctly excluded -- and said so out loud.
        assert any("ADP x405" in note for note in s.notes)
        s.validate()
        assert elapsed < 10.0, f"parse took {elapsed:.2f}s"

    def test_multi_model_file_is_not_merged(self) -> None:
        path = DATA / "test.pdb"
        total = sum(1 for line in path.read_text().splitlines() if line.startswith("ATOM"))
        first = read_pdb(path)
        assert first.n_atoms == total // 3
        assert any("Read model 1 of 3" in note for note in first.notes)

        third = read_pdb(path, model=3)
        assert third.n_atoms == first.n_atoms
        with pytest.raises(StructureFileError, match="only 3 model"):
            read_pdb(path, model=4)

    def test_every_fixture_round_trips_through_validate(self) -> None:
        for name in ("p300.pdb", "dnmt3a.pdb", "arf19.pdb", "test.pdb", "testing_translation.pdb"):
            s = read_pdb(DATA / name)
            assert isinstance(s, Structure)
            s.validate()
            assert s.ca_indices.shape == (s.n_residues,)
            assert not np.isnan(s.xyz).any()
