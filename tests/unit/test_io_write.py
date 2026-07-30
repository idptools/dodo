"""Tests for the PDB writer.

Most of these are named for a specific defect in one of the two writers this module
replaces. The writer is the last step before the user sees anything, so a defect here
discards all the work upstream of it -- that is why the column-level assertions are
this pedantic.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from dodo.constants import (
    BETA_DISORDERED,
    BETA_FOLDED,
    DEFAULT_BOX_DIMENSIONS,
    MAX_PDB_ATOM_SERIAL,
    MAX_PDB_RESIDUE_NUMBER,
    RESIDUES_PER_SEQRES_LINE,
)
from dodo.exceptions import (
    EmptyStructureError,
    GeometryError,
    InvalidRegionError,
    StructureFileError,
)
from dodo.io.write import _hybrid36, structure_to_pdb_lines, write_pdb
from dodo.structure import Domain, DomainKind, Span, Structure

FIXTURES = Path(__file__).parent.parent / "data" / "structures"

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def make_structure(
    n_residues: int = 5,
    *,
    all_atom: bool = True,
    chain_id: str = "A",
    first_residue_number: int = 1,
    b_factor: float | list[float] = 50.0,
    occupancy: float | list[float] = 1.0,
    residue_name: str = "ALA",
    x_offset: float = 0.0,
) -> Structure:
    """Build a straight-line poly-alanine chain with realistic backbone geometry.

    The intra-residue offsets are chosen so that C(i)-N(i+1) is 1.30 A and CA(i)-CA(i+1)
    is 3.80 A, which is what the CONECT tests need: a fixture with arbitrary jitter would
    make the peptide-bond tests pass or fail for reasons unrelated to the writer.
    """
    names = ("N", "CA", "C", "O") if all_atom else ("CA",)
    elements = ("N", "C", "C", "O") if all_atom else ("C",)
    offsets = ((-1.0, 0.3), (0.0, 0.0), (1.5, 0.3), (1.6, 1.5)) if all_atom else ((0.0, 0.0),)

    xyz: list[list[float]] = []
    atom_name: list[str] = []
    element: list[str] = []
    resname: list[str] = []
    resnum: list[int] = []
    chain: list[str] = []
    for i in range(n_residues):
        base = x_offset + i * 3.8
        for name, symbol, (dx, dy) in zip(names, elements, offsets, strict=True):
            xyz.append([base + dx, dy, 0.0])
            atom_name.append(name)
            element.append(symbol)
            resname.append(residue_name)
            resnum.append(first_residue_number + i)
            chain.append(chain_id)

    n_atoms = len(atom_name)
    b_values = [b_factor] * n_atoms if isinstance(b_factor, float) else b_factor
    occ_values = [occupancy] * n_atoms if isinstance(occupancy, float) else occupancy
    return Structure.from_atom_records(
        xyz=np.array(xyz),
        atom_name=atom_name,
        element=element,
        residue_name=resname,
        residue_number=resnum,
        chain_id=chain,
        b_factor=b_values,
        occupancy=occ_values,
        source="test-fixture",
    )


def make_duplicate_ca_structure(*, all_atom: bool = False) -> Structure:
    """Build a structure whose first residue carries two CA atoms.

    Reachable through DODO's own reader: ``read_pdb`` keeps ``blank | (rank ==
    preferred[group])``, so a residue with a blank-altLoc CA *and* an altLoc-``A`` CA
    survives with both. Every writer path has to have an answer for it.
    """
    if all_atom:
        rows = [
            ("N", "N", -1.0, 0.3),
            ("CA", "C", 0.0, 0.0),
            ("CA", "C", 0.1, 0.0),  # the second conformer's alpha carbon
            ("C", "C", 1.5, 0.3),
            ("O", "O", 1.6, 1.5),
        ]
        second = [("N", "N", 2.8, 0.3), ("CA", "C", 3.8, 0.0), ("C", "C", 5.3, 0.3)]
    else:
        rows = [("CA", "C", 0.0, 0.0), ("CA", "C", 0.1, 0.0)]
        second = [("CA", "C", 3.8, 0.0)]

    xyz: list[list[float]] = []
    atom_name: list[str] = []
    element: list[str] = []
    resnum: list[int] = []
    for index, group in enumerate((rows, second), start=1):
        for name, symbol, x, y in group:
            xyz.append([x, y, 0.0])
            atom_name.append(name)
            element.append(symbol)
            resnum.append(index)
    return Structure.from_atom_records(
        xyz=np.array(xyz),
        atom_name=atom_name,
        element=element,
        residue_name=["ALA"] * len(atom_name),
        residue_number=resnum,
        chain_id=["A"] * len(atom_name),
        source="duplicate-ca-fixture",
    )


def atom_lines(lines: list[str]) -> list[str]:
    """Select the coordinate records from a rendered structure."""
    return [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]


def conect_pairs(lines: list[str]) -> set[tuple[int, int]]:
    """Parse CONECT records into a set of serial pairs."""
    pairs: set[tuple[int, int]] = set()
    for line in lines:
        if not line.startswith("CONECT"):
            continue
        origin = int(line[6:11])
        for start in range(11, len(line.rstrip()), 5):
            field = line[start : start + 5].strip()
            if field:
                pairs.add((origin, int(field)))
    return pairs


# ---------------------------------------------------------------------------
# The atom-name column rule
# ---------------------------------------------------------------------------


class TestAtomNameColumns:
    """The regression that made v1's CA-only output load as 912 calcium atoms.

    Per the PDB specification, the element symbol is right-justified in columns 13-14,
    so a one-letter element with a 1-3 character name starts in column 14. v1
    right-justified the whole name from column 13 and wrote no element column; MDTraj
    then read every ``CA`` as calcium. ``tests/data/structures/testing_translation.pdb``
    is v1 output and still shows the defect, which is what this test locks out.
    """

    def test_ca_and_element_columns(self) -> None:
        line = atom_lines(structure_to_pdb_lines(make_structure(2, all_atom=False)))[0]
        assert line[12:16] == " CA "
        assert line[76:78] == " C"

    def test_v1_output_shows_the_defect_this_test_prevents(self) -> None:
        # Guards the claim in the docstring above: if this fixture is ever regenerated
        # with the fixed writer, the comparison below stops being meaningful.
        with open(FIXTURES / "testing_translation.pdb") as handle:
            first = next(line for line in handle if line.startswith("ATOM"))
        assert first[12:16] == "CA  ", "v1 wrote the name left-justified from column 13"
        assert first[76:78] == " C", "but did eventually write an element column"

    def test_single_character_name_starts_in_column_14(self) -> None:
        lines = atom_lines(structure_to_pdb_lines(make_structure(1)))
        assert lines[0][12:16] == " N  "
        assert lines[0][76:78] == " N"

    @pytest.mark.parametrize(
        ("name", "element", "expected"),
        [
            ("CA", "C", " CA "),
            ("N", "N", " N  "),
            ("OD1", "O", " OD1"),
            ("HD11", "H", "HD11"),  # a 4-character name is the only one starting in 13
            ("SE", "SE", "SE  "),  # selenomethionine: a real two-letter element
            ("ZN", "ZN", "ZN  "),  # and a metal ion, as in tests/data/structures/6kn7
        ],
    )
    def test_name_justification(self, name: str, element: str, expected: str) -> None:
        s = Structure.from_atom_records(
            xyz=np.zeros((1, 3)),
            atom_name=[name],
            element=[element],
            residue_name=["ALA"],
            residue_number=[1],
            chain_id=["A"],
        )
        assert atom_lines(structure_to_pdb_lines(s))[0][12:16] == expected

    def test_full_column_layout(self) -> None:
        s = make_structure(1, all_atom=False, first_residue_number=42, b_factor=77.5)
        line = atom_lines(structure_to_pdb_lines(s))[0]
        assert line[0:6] == "ATOM  "
        assert line[6:11] == "    1"
        assert line[16] == " ", "altLoc column must be blank"
        assert line[17:20] == "ALA"
        assert line[21] == "A"
        assert line[22:26] == "  42"
        assert line[26] == " ", "insertion-code column blank when there is none"
        assert line[30:38] == "   0.000"
        assert line[54:60] == "  1.00"
        assert line[60:66] == " 77.50"
        assert len(line) == 78

    def test_insertion_code_is_written(self) -> None:
        s = Structure.from_atom_records(
            xyz=np.zeros((2, 3)),
            atom_name=["CA", "CA"],
            element=["C", "C"],
            residue_name=["GLU", "GLU"],
            residue_number=[142, 142],
            insertion_code=["", "B"],
            chain_id=["A", "A"],
        )
        lines = atom_lines(structure_to_pdb_lines(s))
        assert lines[0][22:27] == " 142 "
        assert lines[1][22:27] == " 142B"

    def test_nonstandard_residue_is_a_hetatm(self) -> None:
        s = Structure.from_atom_records(
            xyz=np.zeros((1, 3)),
            atom_name=["O"],
            element=["O"],
            residue_name=["HOH"],
            residue_number=[1],
            chain_id=["A"],
        )
        assert structure_to_pdb_lines(s)[0].startswith("HETATM")


# ---------------------------------------------------------------------------
# CONECT
# ---------------------------------------------------------------------------


class TestConect:
    """Without these records a rebuilt CA-only IDR renders as disconnected dots.

    CA-CA is 3.8 A, past the automatic bond-detection cutoff in both VMD and PyMOL, so
    the bonds have to be stated explicitly. Neither previous writer emitted a single
    CONECT record.
    """

    def test_ca_only_chain_is_fully_connected(self) -> None:
        lines = structure_to_pdb_lines(make_structure(5, all_atom=False))
        assert conect_pairs(lines) == {(1, 2), (2, 3), (3, 4), (4, 5)}

    def test_all_atom_backbone_bonds(self) -> None:
        lines = structure_to_pdb_lines(make_structure(2))
        # Serials 1-4 are N, CA, C, O of residue 1; 5-8 the same for residue 2.
        assert conect_pairs(lines) == {
            (1, 2),  # N-CA
            (2, 3),  # CA-C
            (3, 4),  # C-O
            (3, 5),  # the peptide bond, C(i)-N(i+1)
            (5, 6),
            (6, 7),
            (7, 8),
        }

    def test_no_ca_ca_bond_on_top_of_a_peptide_bond(self) -> None:
        # A CA(i)-CA(i+1) record in an all-atom model would draw a spurious 3.8 A bond
        # across every peptide unit.
        assert (2, 6) not in conect_pairs(structure_to_pdb_lines(make_structure(2)))

    def test_conect_can_be_disabled(self) -> None:
        lines = structure_to_pdb_lines(make_structure(3, all_atom=False), conect=False)
        assert not any(line.startswith("CONECT") for line in lines)

    def test_never_bonds_across_a_chain_break(self) -> None:
        """Two chains whose termini happen to sit 3.8 A apart are still two chains."""
        left = make_structure(3, all_atom=False, chain_id="A")
        right = make_structure(3, all_atom=False, chain_id="B", x_offset=3 * 3.8)
        joined = Structure.from_atom_records(
            xyz=np.vstack([left.xyz, right.xyz]),
            atom_name=[*left.atom_name, *right.atom_name],
            element=[*left.element, *right.element],
            residue_name=["ALA"] * 6,
            residue_number=[1, 2, 3, 1, 2, 3],
            chain_id=["A"] * 3 + ["B"] * 3,
        )
        # The A/B junction really is at bonding distance, so only the chain test can
        # exclude it.
        assert np.linalg.norm(joined.xyz[2] - joined.xyz[3]) == pytest.approx(3.8)
        pairs = conect_pairs(structure_to_pdb_lines(joined))
        assert (1, 2) in pairs and (2, 3) in pairs
        assert not any(low <= 3 < high for low, high in pairs), "bonded across the TER"

    def test_never_bonds_across_a_gap_in_the_coordinates(self) -> None:
        """An unmodelled region leaves the flanking CAs far apart; do not draw a bond."""
        s = make_structure(4, all_atom=False)
        s.xyz[2:] += np.array([30.0, 0.0, 0.0])
        pairs = conect_pairs(structure_to_pdb_lines(s))
        assert pairs == {(1, 2), (3, 4)}

    def test_serials_account_for_ter_consuming_one(self) -> None:
        """TER takes a serial number, so CONECT cannot assume serial == index + 1."""
        left = make_structure(2, all_atom=False, chain_id="A")
        right = make_structure(2, all_atom=False, chain_id="B", x_offset=100.0)
        joined = Structure.from_atom_records(
            xyz=np.vstack([left.xyz, right.xyz]),
            atom_name=["CA"] * 4,
            element=["C"] * 4,
            residue_name=["ALA"] * 4,
            residue_number=[1, 2, 1, 2],
            chain_id=["A", "A", "B", "B"],
        )
        lines = structure_to_pdb_lines(joined)
        coordinates = atom_lines(lines)
        assert [line[6:11].strip() for line in coordinates] == ["1", "2", "4", "5"]
        # Chain B's bond is 4-5, not 3-4: getting this wrong bonds the TER record.
        assert conect_pairs(lines) == {(1, 2), (4, 5)}

    def test_atom_order_within_a_residue_does_not_change_the_bonds(self) -> None:
        """Bonds are looked up by atom name, not by position within the residue."""
        s = Structure.from_atom_records(
            xyz=np.array(
                [
                    [0.0, 0.0, 0.0],  # CA written first, then N, C, O
                    [-1.0, 0.3, 0.0],
                    [1.5, 0.3, 0.0],
                    [1.6, 1.5, 0.0],
                ]
            ),
            atom_name=["CA", "N", "C", "O"],
            element=["C", "N", "C", "O"],
            residue_name=["ALA"] * 4,
            residue_number=[1] * 4,
            chain_id=["A"] * 4,
        )
        # Serial 1 is the CA, so both N-CA and CA-C are listed from it.
        assert conect_pairs(structure_to_pdb_lines(s)) == {(1, 2), (1, 3), (3, 4)}

    def test_repeated_backbone_atom_name_is_reported_as_unbonded(self) -> None:
        """An all-atom residue with two CAs writes an atom no CONECT record reaches.

        The coordinates are kept -- an all-atom write is lossless by definition -- but
        one atom ends up with no bond at all, so the writer has to say so instead of
        emitting a silent disconnected dot.
        """
        s = make_duplicate_ca_structure(all_atom=True)
        lines = structure_to_pdb_lines(s)
        written = {int(line[6:11]) for line in atom_lines(lines)}
        bonded = {serial for pair in conect_pairs(lines) for serial in pair}
        assert written - bonded, "fixture should leave the repeated CA unbonded"
        assert any("repeat a backbone atom name" in note for note in s.notes)

    def test_conect_records_stay_inside_their_columns(self) -> None:
        lines = structure_to_pdb_lines(make_structure(6))
        records = [line for line in lines if line.startswith("CONECT")]
        assert records
        # "CONECT" plus at most four 5-column serial fields.
        assert all(len(line) <= 6 + 5 * (1 + 4) for line in records)


# ---------------------------------------------------------------------------
# TER, END, MODEL
# ---------------------------------------------------------------------------


class TestRecordStructure:
    def test_ter_per_chain_and_single_end(self, tmp_path: Path) -> None:
        s = Structure.from_atom_records(
            xyz=np.arange(12, dtype=float).reshape(4, 3),
            atom_name=["CA"] * 4,
            element=["C"] * 4,
            residue_name=["ALA"] * 4,
            residue_number=[1, 2, 1, 2],
            chain_id=["A", "A", "B", "B"],
        )
        path = tmp_path / "two_chains.pdb"
        write_pdb(s, path)
        lines = path.read_text().splitlines()
        assert sum(line.startswith("TER") for line in lines) == 2
        assert lines[-1] == "END"
        assert sum(line == "END" for line in lines) == 1

    def test_ter_names_the_last_residue_of_the_chain(self) -> None:
        lines = structure_to_pdb_lines(make_structure(3, all_atom=False, first_residue_number=7))
        ter = next(line for line in lines if line.startswith("TER"))
        assert ter[17:20] == "ALA"
        assert ter[21] == "A"
        assert ter[22:26] == "   9"

    def test_single_structure_has_no_model_records(self, tmp_path: Path) -> None:
        path = tmp_path / "one.pdb"
        write_pdb(make_structure(3), path)
        text = path.read_text()
        assert "MODEL" not in text and "ENDMDL" not in text

    def test_multi_model_frames(self, tmp_path: Path) -> None:
        """v1 supported this; the v2 attempt's function body was literally ``pass``."""
        frames = [make_structure(3, all_atom=False) for _ in range(3)]
        for i, frame in enumerate(frames):
            frame.xyz += np.array([float(i), 0.0, 0.0])
        path = tmp_path / "traj.pdb"
        write_pdb(frames, path)
        lines = path.read_text().splitlines()
        assert [line for line in lines if line.startswith("MODEL")] == [
            "MODEL        1",
            "MODEL        2",
            "MODEL        3",
        ]
        assert sum(line == "ENDMDL" for line in lines) == 3
        assert lines[-1] == "END"
        # Each frame carries its own coordinates, and they differ.
        first_x = [line[30:38] for line in atom_lines(lines)]
        assert first_x[0] != first_x[3]

    def test_conect_written_once_after_the_last_endmdl(self, tmp_path: Path) -> None:
        path = tmp_path / "traj.pdb"
        write_pdb([make_structure(4, all_atom=False) for _ in range(2)], path)
        lines = path.read_text().splitlines()
        conect_indices = [i for i, line in enumerate(lines) if line.startswith("CONECT")]
        last_endmdl = max(i for i, line in enumerate(lines) if line == "ENDMDL")
        assert conect_indices, "a pseudo-trajectory still needs its bonds"
        assert min(conect_indices) > last_endmdl
        assert conect_pairs(lines) == {(1, 2), (2, 3), (3, 4)}

    def test_frames_must_share_a_topology(self, tmp_path: Path) -> None:
        frames = [make_structure(3, all_atom=False), make_structure(4, all_atom=False)]
        with pytest.raises(StructureFileError, match="frames of one file"):
            write_pdb(frames, tmp_path / "bad.pdb")

    def test_frames_must_share_residue_identities(self, tmp_path: Path) -> None:
        frames = [
            make_structure(3, all_atom=False),
            make_structure(3, all_atom=False, first_residue_number=100),
        ]
        with pytest.raises(StructureFileError, match="residue numbers"):
            write_pdb(frames, tmp_path / "bad.pdb")

    def test_frames_must_share_chain_indices(self, tmp_path: Path) -> None:
        """chain_index decides TER placement, and TER consumes an atom serial.

        Two frames agreeing on atom names, residue identities and the chain-id list can
        still put the chain break in a different place. The CONECT records are written
        once, from frame 1, so frame 2's bonds would then point at the wrong atoms --
        including at a TER serial -- in a file that loads without complaint.
        """

        def frame(chain_ids: list[str]) -> Structure:
            xyz = np.zeros((6, 3))
            xyz[:, 0] = np.arange(6) * 3.8
            return Structure.from_atom_records(
                xyz=xyz,
                atom_name=["CA"] * 6,
                element=["C"] * 6,
                residue_name=["ALA"] * 6,
                residue_number=[1, 2, 3, 4, 5, 6],
                chain_id=chain_ids,
            )

        first = frame(["A", "A", "A", "B", "B", "B"])
        second = frame(["A", "A", "B", "B", "B", "B"])
        # Everything the old check looked at agrees, which is why this got through.
        assert np.array_equal(first.residue_index, second.residue_index)
        assert np.array_equal(first.residue_number, second.residue_number)
        assert [c.chain_id for c in first.chains] == [c.chain_id for c in second.chains]
        with pytest.raises(StructureFileError, match="chain indices"):
            write_pdb([first, second], tmp_path / "bad.pdb")

    def test_models_as_frames_false_writes_one_file_each(self, tmp_path: Path) -> None:
        frames = [make_structure(3, all_atom=False), make_structure(4, all_atom=False)]
        write_pdb(frames, tmp_path / "conformer.pdb", models_as_frames=False)
        first = tmp_path / "conformer_1.pdb"
        second = tmp_path / "conformer_2.pdb"
        assert first.is_file() and second.is_file()
        assert len(atom_lines(first.read_text().splitlines())) == 3
        assert len(atom_lines(second.read_text().splitlines())) == 4

    def test_a_one_element_ensemble_is_still_an_ensemble(self, tmp_path: Path) -> None:
        """N conformers must not be renamed or reshaped just because N happens to be 1."""
        s = make_structure(3, all_atom=False)
        write_pdb([s], tmp_path / "frames.pdb")
        assert "MODEL        1" in (tmp_path / "frames.pdb").read_text()
        write_pdb([s], tmp_path / "conformer.pdb", models_as_frames=False)
        assert (tmp_path / "conformer_1.pdb").is_file()
        assert not (tmp_path / "conformer.pdb").exists()

    def test_empty_sequence_rejected(self, tmp_path: Path) -> None:
        with pytest.raises(EmptyStructureError):
            write_pdb([], tmp_path / "nothing.pdb")


# ---------------------------------------------------------------------------
# CRYST1 and SEQRES
# ---------------------------------------------------------------------------


class TestOptionalRecords:
    def test_no_cryst1_by_default(self) -> None:
        assert not any(
            line.startswith("CRYST1") for line in structure_to_pdb_lines(make_structure(2))
        )

    def test_cryst1_when_box_given(self) -> None:
        line = structure_to_pdb_lines(make_structure(2), box=(120.0, 130.5, 140.0))[0]
        assert line.startswith("CRYST1")
        assert line[6:15] == "  120.000"
        assert line[15:24] == "  130.500"
        assert line[24:33] == "  140.000"
        assert line[33:40] == "  90.00"
        assert line[55:66].strip() == "P 1"

    def test_default_box_is_available_explicitly(self) -> None:
        line = structure_to_pdb_lines(make_structure(2), box=DEFAULT_BOX_DIMENSIONS)[0]
        assert line[6:15] == f"{DEFAULT_BOX_DIMENSIONS[0]:9.3f}"

    def test_degenerate_box_rejected(self) -> None:
        with pytest.raises(GeometryError, match="positive and finite"):
            structure_to_pdb_lines(make_structure(2), box=(0.0, 10.0, 10.0))

    def test_no_seqres_by_default(self) -> None:
        assert not any(
            line.startswith("SEQRES") for line in structure_to_pdb_lines(make_structure(2))
        )

    def test_seqres_wraps_at_the_constant(self) -> None:
        n = RESIDUES_PER_SEQRES_LINE * 2 + 1
        lines = [
            line
            for line in structure_to_pdb_lines(make_structure(n, all_atom=False), seqres=True)
            if line.startswith("SEQRES")
        ]
        assert len(lines) == 3
        assert lines[0][19:].split() == ["ALA"] * RESIDUES_PER_SEQRES_LINE
        assert lines[2][19:].split() == ["ALA"]
        assert lines[0][7:10] == "  1"
        assert lines[2][7:10] == "  3"
        assert lines[0][11] == "A"
        assert int(lines[0][13:17]) == n

    def test_seqres_follows_the_pdb_record_order(self) -> None:
        """SEQRES is primary structure; CRYST1 is crystallographic and comes after it."""
        lines = structure_to_pdb_lines(
            make_structure(2, all_atom=False), seqres=True, box=DEFAULT_BOX_DIMENSIONS
        )
        order = [line[:6] for line in lines]
        assert order.index("SEQRES") < order.index("CRYST1")

    def test_seqres_follows_the_pdb_record_order_in_multi_frame_output(
        self, tmp_path: Path
    ) -> None:
        path = tmp_path / "traj.pdb"
        write_pdb(
            [make_structure(2, all_atom=False) for _ in range(2)],
            path,
            seqres=True,
            box=DEFAULT_BOX_DIMENSIONS,
        )
        order = [line[:6] for line in path.read_text().splitlines()]
        assert order.index("SEQRES") < order.index("CRYST1")
        assert order.index("CRYST1") < order.index("MODEL ")

    def test_seqres_never_repeats_a_chain_id(self) -> None:
        """Interleaved chain runs must not produce two SEQRES blocks for one chain.

        Both DODO readers deliberately preserve interleaved runs as separate chains, so
        a structure can hold chain A twice. Emitting a serNum-1 block per run declares
        chain A as 3 residues while the ATOM records give it 6 -- a file whose SEQRES
        contradicts its own coordinates.
        """
        xyz = np.zeros((9, 3))
        xyz[:, 0] = np.arange(9) * 3.8
        s = Structure.from_atom_records(
            xyz=xyz,
            atom_name=["CA"] * 9,
            element=["C"] * 9,
            residue_name=["ALA"] * 9,
            residue_number=[1, 2, 3, 1, 2, 3, 10, 11, 12],
            chain_id=["A"] * 3 + ["B"] * 3 + ["A"] * 3,
        )
        assert [(c.chain_id, c.span.start, c.span.stop) for c in s.chains] == [
            ("A", 0, 3),
            ("B", 3, 6),
            ("A", 6, 9),
        ]
        lines = [line for line in structure_to_pdb_lines(s, seqres=True) if line[:6] == "SEQRES"]
        chain_ids = [line[11] for line in lines]
        assert len(set(chain_ids)) == len(chain_ids), "one SEQRES block per chain id"
        by_chain = dict(zip(chain_ids, lines, strict=True))
        assert int(by_chain["A"][13:17]) == 6
        assert int(by_chain["B"][13:17]) == 3
        assert by_chain["A"][19:].split() == ["ALA"] * 6
        assert any("SEQRES" in note for note in s.notes)

    def test_seqres_lists_only_the_residues_that_were_written(self) -> None:
        """With ``ca_only`` a ligand is not written, so SEQRES must not declare it."""
        s = Structure.from_atom_records(
            xyz=np.array([[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [20.0, 0.0, 0.0]]),
            atom_name=["CA", "CA", "O"],
            element=["C", "C", "O"],
            residue_name=["ALA", "ALA", "HOH"],
            residue_number=[1, 2, 100],
            chain_id=["A", "A", "A"],
        )
        lines = [
            line
            for line in structure_to_pdb_lines(s, ca_only=True, seqres=True)
            if line[:6] == "SEQRES"
        ]
        assert lines
        assert "HOH" not in " ".join(lines)
        assert int(lines[0][13:17]) == 2

    def test_seqres_per_chain(self) -> None:
        s = Structure.from_atom_records(
            xyz=np.arange(9, dtype=float).reshape(3, 3),
            atom_name=["CA"] * 3,
            element=["C"] * 3,
            residue_name=["GLY", "SER", "TRP"],
            residue_number=[1, 2, 1],
            chain_id=["A", "A", "B"],
        )
        lines = [line for line in structure_to_pdb_lines(s, seqres=True) if line[:6] == "SEQRES"]
        assert lines[0][11] == "A"
        assert lines[0][19:].split() == ["GLY", "SER"]
        assert lines[1][11] == "B"
        assert lines[1][19:].split() == ["TRP"]


# ---------------------------------------------------------------------------
# Region annotation and CA-only output
# ---------------------------------------------------------------------------


class TestAnnotateRegions:
    """Folded is 100 and disordered is 0.

    That polarity means colouring by B-factor lights up the folded core. v1's
    docstrings claimed the reverse in three places; its code, README and CLI help did
    it this way, and the code was right.
    """

    def _annotated(self) -> list[str]:
        s = make_structure(10, all_atom=False)
        chain = s.chains[0]
        chain.domains = [
            Domain(s, Span(0, 4), DomainKind.FOLDED),
            Domain(s, Span(4, 7), DomainKind.IDR),
            Domain(s, Span(7, 10), DomainKind.FOLDED),
        ]
        return atom_lines(structure_to_pdb_lines(s, annotate_regions=True))

    def test_folded_is_beta_folded(self) -> None:
        lines = self._annotated()
        assert float(lines[0][60:66]) == pytest.approx(BETA_FOLDED)
        assert float(lines[9][60:66]) == pytest.approx(BETA_FOLDED)

    def test_idr_is_beta_disordered(self) -> None:
        lines = self._annotated()
        assert [float(line[60:66]) for line in lines[4:7]] == [BETA_DISORDERED] * 3

    def test_folded_is_one_hundred_and_disordered_is_zero(self) -> None:
        # Locks the polarity itself, not just the constants.
        assert BETA_FOLDED == 100.0
        assert BETA_DISORDERED == 0.0

    def test_loops_inside_a_folded_domain_count_as_rebuilt(self) -> None:
        s = make_structure(20, all_atom=False)
        s.chains[0].domains = [
            Domain(s, Span(0, 20), DomainKind.FOLDED, loops=(Span(5, 8),)),
        ]
        lines = atom_lines(structure_to_pdb_lines(s, annotate_regions=True))
        assert [float(line[60:66]) for line in lines[5:8]] == [BETA_DISORDERED] * 3
        assert float(lines[4][60:66]) == pytest.approx(BETA_FOLDED)

    def test_without_the_flag_the_input_b_factors_survive(self) -> None:
        s = make_structure(3, all_atom=False, b_factor=42.25)
        s.chains[0].domains = [Domain(s, Span(0, 3), DomainKind.FOLDED)]
        lines = atom_lines(structure_to_pdb_lines(s))
        assert float(lines[0][60:66]) == pytest.approx(42.25)

    def test_no_domains_is_an_error_not_a_silent_no_op(self) -> None:
        with pytest.raises(InvalidRegionError, match="no domains"):
            structure_to_pdb_lines(make_structure(3), annotate_regions=True)

    def test_unassigned_residues_are_recorded(self) -> None:
        s = make_structure(6, all_atom=False)
        s.chains[0].domains = [Domain(s, Span(0, 3), DomainKind.FOLDED)]
        structure_to_pdb_lines(s, annotate_regions=True)
        assert any("belong to no domain" in note for note in s.notes)


class TestCaOnly:
    def test_writes_only_alpha_carbons(self) -> None:
        lines = atom_lines(structure_to_pdb_lines(make_structure(5), ca_only=True))
        assert len(lines) == 5
        assert all(line[12:16] == " CA " for line in lines)

    def test_ca_only_output_is_still_bonded(self) -> None:
        lines = structure_to_pdb_lines(make_structure(4), ca_only=True)
        assert conect_pairs(lines) == {(1, 2), (2, 3), (3, 4)}

    def test_residues_without_a_ca_are_reported(self) -> None:
        s = Structure.from_atom_records(
            xyz=np.array([[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [20.0, 0.0, 0.0]]),
            atom_name=["CA", "CA", "O"],
            element=["C", "C", "O"],
            residue_name=["ALA", "ALA", "HOH"],
            residue_number=[1, 2, 100],
            chain_id=["A", "A", "A"],
        )
        assert len(atom_lines(structure_to_pdb_lines(s, ca_only=True))) == 2
        assert any("no CA atom" in note for note in s.notes)

    def test_duplicate_ca_writes_one_per_residue_and_reports_it(self) -> None:
        """``ca_only`` means exactly one CA per residue, matching ``Structure.ca_indices``.

        Writing both alpha carbons of a two-conformer residue produces an atom that no
        CONECT record can reach -- the disconnected dot this module exists to prevent --
        and the residue accounting used to invert, reporting "omitted -1 residue(s)".
        """
        s = make_duplicate_ca_structure()
        lines = atom_lines(structure_to_pdb_lines(s, ca_only=True))
        assert len(lines) == s.n_residues == 2
        # The atom written is the one ca_indices resolves to, so the file agrees with
        # every geometric calculation DODO performs on the same structure.
        assert [float(line[30:38]) for line in lines] == pytest.approx(
            s.ca_xyz[:, 0].tolist(), abs=1e-3
        )
        assert any("more than one CA atom" in note for note in s.notes)
        assert not any("-1 residue" in note for note in s.notes)
        assert not any("no CA atom" in note for note in s.notes)

    def test_duplicate_ca_output_is_fully_bonded(self) -> None:
        s = make_duplicate_ca_structure()
        lines = structure_to_pdb_lines(s, ca_only=True)
        assert conect_pairs(lines) == {(1, 2)}

    def test_ca_only_ignores_a_ligand_it_will_not_write(self) -> None:
        """A residue that is excluded from the output must not abort the write.

        Formatting residue-level fields for every residue meant an unwritable ligand
        residue number raised, so the writer refused to produce a perfectly valid CA
        trace because of a ligand it was never going to write.
        """
        s = Structure.from_atom_records(
            xyz=np.array([[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], [20.0, 0.0, 0.0]]),
            atom_name=["CA", "CA", "PB"],
            element=["C", "C", "P"],
            residue_name=["ALA", "ALA", "ADP"],
            residue_number=[1, 2, 9_000_000],
            chain_id=["A", "A", "A"],
        )
        lines = atom_lines(structure_to_pdb_lines(s, ca_only=True))
        assert len(lines) == 2
        assert any("no CA atom" in note for note in s.notes)

    def test_ca_only_reports_no_substitution_for_unwritten_residues(self) -> None:
        """Occupancy/B-factor substitutions are only real if the residue is written."""
        s = Structure.from_atom_records(
            xyz=np.array([[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]]),
            atom_name=["CA", "O"],
            element=["C", "O"],
            residue_name=["ALA", "HOH"],
            residue_number=[1, 100],
            chain_id=["A", "A"],
            occupancy=[1.0, math.nan],
        )
        assert len(atom_lines(structure_to_pdb_lines(s, ca_only=True))) == 1
        assert not any("occupancy/B-factor" in note for note in s.notes)

    def test_no_ca_at_all_raises(self) -> None:
        s = Structure.from_atom_records(
            xyz=np.zeros((1, 3)),
            atom_name=["O"],
            element=["O"],
            residue_name=["HOH"],
            residue_number=[1],
            chain_id=["A"],
        )
        with pytest.raises(EmptyStructureError, match="no CA atoms"):
            structure_to_pdb_lines(s, ca_only=True)


# ---------------------------------------------------------------------------
# Overflow
# ---------------------------------------------------------------------------


def hy36decode(field: str) -> int:
    """Decode a hybrid-36 field, for use as an independent check in these tests."""
    text = field.strip()
    if text.lstrip("-").isdigit():
        return int(text)
    width = len(text)
    if text[0].isupper():
        value = int(text, 36)
        return value - 10 * 36 ** (width - 1) + 10**width
    value = int(text.upper(), 36)
    return value - 10 * 36 ** (width - 1) + 10**width + 26 * 36 ** (width - 1)


class TestOverflow:
    """PDB has 5 columns for the atom serial and 4 for the residue number.

    Above those, real files use hybrid-36 rather than wrapping or letting the field
    overflow into its neighbour -- which would shift every column after it and produce
    a file that loads without error and is wrong.
    """

    @pytest.mark.parametrize(
        ("value", "width", "limit", "expected"),
        [
            (1, 5, MAX_PDB_ATOM_SERIAL, "    1"),
            (99999, 5, MAX_PDB_ATOM_SERIAL, "99999"),
            (100000, 5, MAX_PDB_ATOM_SERIAL, "A0000"),
            (100001, 5, MAX_PDB_ATOM_SERIAL, "A0001"),
            (87440031, 5, MAX_PDB_ATOM_SERIAL, "zzzzz"),
            (9999, 4, MAX_PDB_RESIDUE_NUMBER, "9999"),
            (10000, 4, MAX_PDB_RESIDUE_NUMBER, "A000"),
            (2436111, 4, MAX_PDB_RESIDUE_NUMBER, "zzzz"),
            (-999, 4, MAX_PDB_RESIDUE_NUMBER, "-999"),  # negative resSeq is legal
        ],
    )
    def test_hybrid36_reference_values(
        self, value: int, width: int, limit: int, expected: str
    ) -> None:
        encoded = _hybrid36(value, width, limit)
        assert encoded == expected
        assert len(encoded) == width
        assert hy36decode(encoded) == value

    @pytest.mark.parametrize("value", [87440032, -1000])
    def test_unrepresentable_value_raises(self, value: int) -> None:
        with pytest.raises(StructureFileError, match="cannot be written"):
            _hybrid36(value, 5 if value > 0 else 4, MAX_PDB_ATOM_SERIAL)

    def test_large_residue_number_uses_hybrid36(self) -> None:
        s = Structure.from_atom_records(
            xyz=np.array([[0.0, 0.0, 0.0], [3.8, 0.0, 0.0]]),
            atom_name=["CA", "CA"],
            element=["C", "C"],
            residue_name=["ALA", "ALA"],
            residue_number=[MAX_PDB_RESIDUE_NUMBER, MAX_PDB_RESIDUE_NUMBER + 1],
            chain_id=["A", "A"],
        )
        lines = atom_lines(structure_to_pdb_lines(s))
        assert lines[0][22:26] == "9999"
        assert lines[1][22:26] == "A000"
        assert hy36decode(lines[1][22:26]) == MAX_PDB_RESIDUE_NUMBER + 1

    @pytest.mark.slow
    def test_large_atom_serial_uses_hybrid36(self) -> None:
        n = MAX_PDB_ATOM_SERIAL + 2
        # A lattice keeps every coordinate inside the %8.3f fields.
        side = 101
        grid = np.stack(np.unravel_index(np.arange(n), (side, side, side)), axis=-1).astype(float)
        s = Structure.from_atom_records(
            xyz=grid * 4.0,
            atom_name=["CA"] * n,
            element=["C"] * n,
            residue_name=["ALA"] * n,
            residue_number=list(range(1, n + 1)),
            chain_id=["A"] * n,
        )
        lines = atom_lines(structure_to_pdb_lines(s, conect=False))
        assert lines[MAX_PDB_ATOM_SERIAL - 1][6:11] == "99999"
        assert lines[MAX_PDB_ATOM_SERIAL][6:11] == "A0000"
        assert hy36decode(lines[MAX_PDB_ATOM_SERIAL + 1][6:11]) == MAX_PDB_ATOM_SERIAL + 2

    def test_encoding_matches_what_the_reader_accepts(self) -> None:
        """The writer's hybrid-36 must be the reader's hybrid-36, over the whole range.

        Two independent implementations of an encoding is exactly the situation where a
        file writes cleanly and reads back as different numbers.
        """
        try:
            from dodo.io.pdb import decode_hybrid36
        except ImportError:  # pragma: no cover - depends on which modules have landed
            pytest.skip("dodo.io.pdb (the reader) is not available yet")
        rng = np.random.default_rng(0)
        for width, limit, ceiling in (
            (5, MAX_PDB_ATOM_SERIAL, 87440031),
            (4, MAX_PDB_RESIDUE_NUMBER, 2436111),
        ):
            values = [1 - 10 ** (width - 1), limit, limit + 1, ceiling]
            values += rng.integers(1 - 10 ** (width - 1), ceiling, size=200).tolist()
            for value in values:
                encoded = _hybrid36(int(value), width, limit)
                assert decode_hybrid36(encoded, width) == int(value)

    def test_coordinate_out_of_range_raises(self) -> None:
        s = make_structure(2, all_atom=False)
        s.xyz[1] = [12345.0, 0.0, 0.0]
        with pytest.raises(StructureFileError, match="coordinate"):
            structure_to_pdb_lines(s)

    def test_non_finite_coordinate_raises(self) -> None:
        """v1's samplers returned NaN on failure and the writer wrote it out."""
        s = make_structure(2, all_atom=False)
        s.xyz[1] = [math.nan, 0.0, 0.0]
        with pytest.raises(GeometryError, match="non-finite"):
            structure_to_pdb_lines(s)


# ---------------------------------------------------------------------------
# Missing values and paths
# ---------------------------------------------------------------------------


class TestMissingValuesAndPaths:
    def test_nan_occupancy_and_b_factor_get_defaults(self) -> None:
        """The previous writer called ``float(None)`` here and crashed."""
        s = make_structure(2, all_atom=False)
        s.occupancy[:] = np.nan
        s.b_factor[:] = np.nan
        line = atom_lines(structure_to_pdb_lines(s))[0]
        assert line[54:60] == "  1.00"
        assert line[60:66] == "  0.00"
        assert any("occupancy/B-factor" in note for note in s.notes)

    def test_out_of_range_b_factor_is_clamped_and_reported(self) -> None:
        s = make_structure(1, all_atom=False, b_factor=123456.0)
        line = atom_lines(structure_to_pdb_lines(s))[0]
        assert line[60:66] == "999.99"
        assert len(line) == 78, "a clamped value must not widen the field"
        assert any("occupancy/B-factor" in note for note in s.notes)

    def test_bare_relative_filename(self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
        """``os.path.dirname("out.pdb")`` is ``""``; the previous writer rejected it."""
        monkeypatch.chdir(tmp_path)
        write_pdb(make_structure(2, all_atom=False), "out.pdb")
        assert (tmp_path / "out.pdb").is_file()

    def test_str_and_path_are_both_accepted(self, tmp_path: Path) -> None:
        s = make_structure(2, all_atom=False)
        write_pdb(s, str(tmp_path / "a.pdb"))
        write_pdb(s, tmp_path / "b.pdb")
        assert (tmp_path / "a.pdb").read_text() == (tmp_path / "b.pdb").read_text()

    def test_missing_directory_is_a_clear_error(self, tmp_path: Path) -> None:
        with pytest.raises(StructureFileError, match="does not exist"):
            write_pdb(make_structure(2), tmp_path / "nope" / "out.pdb")

    def test_existing_directory_as_path_is_a_clear_error(self, tmp_path: Path) -> None:
        """Every "cannot write this" failure is a StructureFileError, not an OSError."""
        with pytest.raises(StructureFileError, match="directory"):
            write_pdb(make_structure(2, all_atom=False), tmp_path)

    def test_non_ascii_chain_id_is_remapped_not_a_unicode_error(self, tmp_path: Path) -> None:
        """A single-character chain id still has to be a character PDB can hold."""
        s = Structure.from_atom_records(
            xyz=np.array([[0.0, 0.0, 0.0], [3.8, 0.0, 0.0]]),
            atom_name=["CA", "CA"],
            element=["C", "C"],
            residue_name=["ALA", "ALA"],
            residue_number=[1, 2],
            chain_id=["Å", "Å"],
        )
        path = tmp_path / "out.pdb"
        write_pdb(s, path)
        lines = atom_lines(path.read_text().splitlines())
        assert [line[21] for line in lines] == ["A", "A"]
        assert any("remapped" in note for note in s.notes)

    def test_non_ascii_field_is_a_structure_file_error(self, tmp_path: Path) -> None:
        """The backstop: no record may leave this module as a UnicodeEncodeError."""
        s = Structure.from_atom_records(
            xyz=np.zeros((1, 3)),
            atom_name=["CA"],
            element=["C"],
            residue_name=["ÅLA"],
            residue_number=[1],
            chain_id=["A"],
        )
        with pytest.raises(StructureFileError, match="ASCII"):
            write_pdb(s, tmp_path / "out.pdb")

    def test_multi_character_chain_id_is_remapped_and_reported(self) -> None:
        """An mmCIF id may be ``AAA``; PDB has one column, and truncating merges chains."""
        s = Structure.from_atom_records(
            xyz=np.arange(12, dtype=float).reshape(4, 3),
            atom_name=["CA"] * 4,
            element=["C"] * 4,
            residue_name=["ALA"] * 4,
            residue_number=[1, 2, 1, 2],
            chain_id=["AAA", "AAA", "AAB", "AAB"],
        )
        lines = atom_lines(structure_to_pdb_lines(s))
        assert [line[21] for line in lines] == ["A", "A", "B", "B"]
        assert any("remapped" in note for note in s.notes)

    def test_missing_element_is_inferred_and_reported(self) -> None:
        s = Structure.from_atom_records(
            xyz=np.zeros((1, 3)),
            atom_name=["CA"],
            element=[""],
            residue_name=["ALA"],
            residue_number=[1],
            chain_id=["A"],
        )
        assert atom_lines(structure_to_pdb_lines(s))[0][76:78] == " C"
        assert any("no element symbol" in note for note in s.notes)

    def test_notes_are_not_duplicated_across_frames(self, tmp_path: Path) -> None:
        s = make_structure(2, all_atom=False)
        s.occupancy[:] = np.nan
        write_pdb([s, s.copy()], tmp_path / "traj.pdb")
        write_pdb(s, tmp_path / "again.pdb")
        assert sum("occupancy/B-factor" in note for note in s.notes) == 1


# ---------------------------------------------------------------------------
# Round trip
# ---------------------------------------------------------------------------


class TestRoundTrip:
    """Coordinates, sequence and residue identity must survive a write/read cycle.

    PDB stores coordinates to 0.001 A, so that is the tolerance. The reader is a
    separate module; if it is not present yet these skip rather than fail.
    """

    @staticmethod
    def _read_pdb(path: Path) -> Structure:
        try:
            from dodo.io import pdb as pdb_reader
        except ImportError:  # pragma: no cover - depends on which modules have landed
            pytest.skip("dodo.io.pdb (the reader) is not available yet")
        for name in ("read_pdb", "read_pdb_file", "read_structure", "parse_pdb"):
            reader = getattr(pdb_reader, name, None)
            if callable(reader):
                return reader(path)  # type: ignore[no-any-return]
        pytest.skip("dodo.io.pdb exposes no recognisable read function")

    def test_dodo_reader_round_trip(self, tmp_path: Path) -> None:
        original = Structure.from_atom_records(
            xyz=np.array(
                [
                    [1.234, -5.678, 9.012],
                    [4.321, -2.100, 8.765],
                    [7.777, 0.001, -3.500],
                    [11.111, 2.222, -3.333],
                ]
            ),
            atom_name=["CA", "CA", "CA", "CA"],
            element=["C", "C", "C", "C"],
            residue_name=["MET", "GLU", "TRP", "GLY"],
            residue_number=[7, 8, 8, 9],
            insertion_code=["", "", "A", ""],
            chain_id=["A", "A", "A", "B"],
            b_factor=[91.5, 72.25, 30.0, 12.75],
        )
        path = tmp_path / "round_trip.pdb"
        write_pdb(original, path)
        reloaded = self._read_pdb(path)

        assert reloaded.n_atoms == original.n_atoms
        assert np.allclose(reloaded.xyz, original.xyz, atol=1e-3)
        assert reloaded.sequence == original.sequence
        assert reloaded.residue_number.tolist() == original.residue_number.tolist()
        assert reloaded.insertion_code.tolist() == original.insertion_code.tolist()
        assert [c.chain_id for c in reloaded.chains] == [c.chain_id for c in original.chains]

    def test_mdtraj_round_trip(self, tmp_path: Path) -> None:
        """An independent tool must agree, including on the element assignment.

        This is the check that catches the v1 column bug: with the name right-justified
        from column 13 and no element column, MDTraj read a CA trace as calcium.
        """
        md = pytest.importorskip("mdtraj", reason="no independent PDB reader installed")
        original = make_structure(6, all_atom=False, first_residue_number=17, b_factor=88.0)
        path = tmp_path / "mdtraj.pdb"
        write_pdb(original, path)

        traj = md.load(str(path))
        assert traj.n_atoms == original.n_atoms
        assert np.allclose(traj.xyz[0] * 10.0, original.xyz, atol=1e-3)
        symbols = {atom.element.symbol for atom in traj.topology.atoms}
        assert symbols == {"C"}, "an alpha carbon must not be read as calcium"
        assert [r.resSeq for r in traj.topology.residues] == original.residue_number.tolist()
        assert [r.name for r in traj.topology.residues] == [
            str(name) for name in original.residue_name
        ]
        # The bonds only exist in the file if CONECT does: 3.8 A is past MDTraj's
        # distance-based bond guessing for a CA trace.
        assert traj.topology.n_bonds == original.n_residues - 1

    def test_mdtraj_reads_multi_model_output_as_frames(self, tmp_path: Path) -> None:
        md = pytest.importorskip("mdtraj", reason="no independent PDB reader installed")
        frames = []
        for i in range(4):
            frame = make_structure(5, all_atom=False)
            frame.xyz += np.array([0.0, float(i), 0.0])
            frames.append(frame)
        path = tmp_path / "frames.pdb"
        write_pdb(frames, path)
        traj = md.load(str(path))
        assert traj.n_frames == 4
        assert np.allclose(traj.xyz[3] * 10.0, frames[3].xyz, atol=1e-3)


# ---------------------------------------------------------------------------
# Integration with real files
# ---------------------------------------------------------------------------


def structure_from_mdtraj(path: Path) -> Structure:
    """Build a Structure from a real fixture using MDTraj as an independent reader.

    DODO's own reader lives in another module and may not exist yet; using MDTraj here
    keeps the integration tests honest about what is being tested, which is the writer.
    """
    md = pytest.importorskip("mdtraj", reason="no independent PDB reader installed")
    traj = md.load(str(path))
    topology = traj.topology
    atoms = list(topology.atoms)
    return Structure.from_atom_records(
        xyz=traj.xyz[0] * 10.0,
        atom_name=[atom.name for atom in atoms],
        element=[atom.element.symbol for atom in atoms],
        residue_name=[atom.residue.name for atom in atoms],
        residue_number=[atom.residue.resSeq for atom in atoms],
        chain_id=[atom.residue.chain.chain_id for atom in atoms],
        source=str(path),
    )


@pytest.mark.slow
class TestRealStructures:
    @pytest.mark.parametrize("name", ["dnmt3a.pdb", "p300.pdb", "6kn7.pdb"])
    def test_real_file_survives_a_dodo_read_write_read_cycle(
        self, name: str, tmp_path: Path
    ) -> None:
        """Read a real structure, write it, read it back: nothing may change.

        6kn7 in particular is the 61,511-atom, 29-chain EM assembly whose interleaved
        records the pre-rewrite reader lost half of.
        """
        try:
            from dodo.io.pdb import read_pdb
        except ImportError:  # pragma: no cover - depends on which modules have landed
            pytest.skip("dodo.io.pdb (the reader) is not available yet")

        original = read_pdb(FIXTURES / name)
        path = tmp_path / f"round_trip_{name}"
        write_pdb(original, path)
        reloaded = read_pdb(path)

        assert reloaded.n_atoms == original.n_atoms
        assert reloaded.n_residues == original.n_residues
        assert np.allclose(reloaded.xyz, original.xyz, atol=1e-3)
        assert reloaded.sequence == original.sequence
        assert reloaded.atom_name.tolist() == original.atom_name.tolist()
        assert reloaded.residue_number.tolist() == original.residue_number.tolist()
        assert reloaded.insertion_code.tolist() == original.insertion_code.tolist()
        assert [c.chain_id for c in reloaded.chains] == [c.chain_id for c in original.chains]
        # 0.005 is half of the B-factor column's last digit.
        assert np.allclose(reloaded.b_factor, original.b_factor, atol=0.005)
        assert np.allclose(reloaded.occupancy, original.occupancy, atol=0.005)

    def test_alphafold_model_survives_a_write_read_cycle(self, tmp_path: Path) -> None:
        md = pytest.importorskip("mdtraj", reason="no independent PDB reader installed")
        original = structure_from_mdtraj(FIXTURES / "dnmt3a.pdb")
        path = tmp_path / "dnmt3a_out.pdb"
        write_pdb(original, path)

        reloaded = md.load(str(path))
        assert reloaded.n_atoms == original.n_atoms
        assert np.allclose(reloaded.xyz[0] * 10.0, original.xyz, atol=1e-3)
        assert [atom.name for atom in reloaded.topology.atoms] == [
            str(name) for name in original.atom_name
        ]
        assert [atom.element.symbol for atom in reloaded.topology.atoms] == [
            str(symbol) for symbol in original.element
        ]
        # A single-chain all-atom model: N-CA, CA-C and C-O in every residue, plus one
        # peptide bond per junction. Anything less means a real backbone went unbonded,
        # which is the failure this module exists to prevent.
        n = original.n_residues
        assert len(conect_pairs(path.read_text().splitlines())) == 3 * n + (n - 1)
        assert original.notes == [], "a clean AlphaFold model should raise no anomalies"

    def test_em_assembly_keeps_every_chain_and_atom(self, tmp_path: Path) -> None:
        """6kn7 is a 29-chain, 7796-residue EM assembly plus 15 ADP ligands.

        This is the class of file v1's reader lost half of, and the ligands are what
        makes ``ca_only`` lossy: they have no alpha carbon, so they cannot be written,
        which the writer records in ``notes`` rather than passing over.
        """
        md = pytest.importorskip("mdtraj", reason="no independent PDB reader installed")
        original = structure_from_mdtraj(FIXTURES / "6kn7.pdb")
        ca_atoms = np.flatnonzero(original.atom_name == "CA")
        assert 0 < ca_atoms.size < original.n_residues, "fixture should contain ligands"

        path = tmp_path / "6kn7_out.pdb"
        write_pdb(original, path, ca_only=True)
        assert any("no CA atom" in note for note in original.notes)

        lines = path.read_text().splitlines()
        assert len(atom_lines(lines)) == ca_atoms.size
        chains_written = {int(original.chain_index[i]) for i in original.residue_index[ca_atoms]}
        assert sum(line.startswith("TER") for line in lines) == len(chains_written)

        reloaded = md.load(str(path))
        assert reloaded.n_atoms == ca_atoms.size
        assert np.allclose(reloaded.xyz[0] * 10.0, original.xyz[ca_atoms], atol=1e-3)
        assert {a.element.symbol for a in reloaded.topology.atoms} == {"C"}
        # Every chain's CA trace is bonded, and no bond crosses a chain boundary.
        assert reloaded.topology.n_bonds > 0
        for bond in reloaded.topology.bonds:
            assert bond[0].residue.chain.index == bond[1].residue.chain.index

    def test_em_assembly_all_atom_write_loses_nothing(self, tmp_path: Path) -> None:
        md = pytest.importorskip("mdtraj", reason="no independent PDB reader installed")
        original = structure_from_mdtraj(FIXTURES / "6kn7.pdb")
        path = tmp_path / "6kn7_all_atom.pdb"
        write_pdb(original, path)

        reloaded = md.load(str(path))
        assert reloaded.n_atoms == original.n_atoms, "not one atom may be dropped"
        assert np.allclose(reloaded.xyz[0] * 10.0, original.xyz, atol=1e-3)
        assert reloaded.topology.n_residues == original.n_residues
