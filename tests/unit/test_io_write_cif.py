"""Tests for the mmCIF writer.

The writer exists because DODO's ``.cif`` output was never mmCIF: it was PDB records in a
file named ``.cif``, which a viewer that trusts the extension hands to its mmCIF parser
and which then fails to open. So the load-bearing assertions here are the ones proving the
output is *real* mmCIF -- parseable by DODO's own reader (which explicitly rejects PDB
text) and, where the parser is installed, by biotite -- and that several conformers are
expressed the mmCIF way, as one ``_atom_site`` loop keyed by ``pdbx_PDB_model_num`` rather
than PDB's ``MODEL``/``ENDMDL`` frames.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dodo.constants import BETA_DISORDERED, BETA_FOLDED, DEFAULT_BOX_DIMENSIONS
from dodo.exceptions import (
    EmptyStructureError,
    GeometryError,
    StructureFileError,
    UnsupportedFormatError,
)
from dodo.io import read_cif, write_cif, write_structure
from dodo.io.cif import parse_cif_text
from dodo.io.cif_write import structure_to_cif_text
from dodo.structure import Domain, DomainKind, Span, Structure

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def make_structure(
    n_residues: int = 5,
    *,
    all_atom: bool = True,
    chain_id: str = "A",
    first_residue_number: int = 1,
    residue_name: str = "ALA",
    base_offset: float = 0.0,
) -> Structure:
    """Build a straight-line poly-alanine chain with CA-CA 3.80 A and C-N 1.30 A.

    Same geometry as the PDB writer's fixture, so the connectivity a reader infers from
    the coordinates is the connectivity the writer is meant to record.
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
        base = base_offset + i * 3.8
        for name, symbol, (dx, dy) in zip(names, elements, offsets, strict=True):
            xyz.append([base + dx, dy, 0.0])
            atom_name.append(name)
            element.append(symbol)
            resname.append(residue_name)
            resnum.append(first_residue_number + i)
            chain.append(chain_id)
    return Structure.from_atom_records(
        xyz=np.array(xyz),
        atom_name=atom_name,
        element=element,
        residue_name=resname,
        residue_number=resnum,
        chain_id=chain,
        source="cif-write-fixture",
    )


def atom_site_rows(text: str) -> list[str]:
    """Return the ``_atom_site`` data rows of a rendered mmCIF document."""
    return [line for line in text.splitlines() if line.startswith(("ATOM", "HETATM"))]


def model_numbers(text: str) -> list[str]:
    """Return the ``pdbx_PDB_model_num`` value (last column) of every atom row."""
    return [line.split()[-1] for line in atom_site_rows(text)]


# ---------------------------------------------------------------------------
# The output is real mmCIF -- the regression that motivated the writer
# ---------------------------------------------------------------------------


class TestOutputIsRealMmcif:
    def test_read_cif_accepts_it_rather_than_calling_it_pdb(self, tmp_path: Path) -> None:
        # The exact failure that started this: PDB-in-a-.cif is rejected by the mmCIF
        # reader with "is PDB format". Real mmCIF must parse.
        path = tmp_path / "out.cif"
        write_cif([make_structure(4, all_atom=False), make_structure(4, all_atom=False)], path)
        structure = read_cif(path)
        assert structure.n_residues == 4

    def test_it_has_the_star_constructs_a_pdb_file_never_has(self) -> None:
        text = structure_to_cif_text(make_structure(3))
        assert text.startswith("data_")
        assert "loop_" in text
        assert "_atom_site.Cartn_x" in text

    def test_parse_cif_text_round_trips_the_coordinates(self) -> None:
        original = make_structure(5)
        parsed = parse_cif_text(structure_to_cif_text(original))
        assert np.allclose(parsed.xyz, original.xyz)


# ---------------------------------------------------------------------------
# Several conformers as pdbx_PDB_model_num, not MODEL/ENDMDL
# ---------------------------------------------------------------------------


class TestMultiModel:
    def test_all_models_share_one_atom_site_loop(self) -> None:
        text = structure_to_cif_text([make_structure(3, all_atom=False) for _ in range(3)])
        assert text.count("_atom_site.Cartn_x") == 1, "one loop, not one per model"

    def test_model_column_numbers_every_frame(self) -> None:
        frames = [make_structure(3, all_atom=False) for _ in range(4)]
        nums = model_numbers(structure_to_cif_text(frames))
        # 3 atoms per frame, four frames, numbered 1..4.
        assert nums == ["1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4"]

    def test_serials_are_globally_unique_across_models(self) -> None:
        frames = [make_structure(3, all_atom=False) for _ in range(3)]
        serials = [int(line.split()[1]) for line in atom_site_rows(structure_to_cif_text(frames))]
        assert serials == list(range(1, 10))

    def test_no_model_or_endmdl_records(self) -> None:
        text = structure_to_cif_text([make_structure(2, all_atom=False) for _ in range(2)])
        assert "MODEL" not in text and "ENDMDL" not in text

    def test_each_model_reads_back_at_its_own_number(self, tmp_path: Path) -> None:
        frames = [
            make_structure(3, all_atom=False, base_offset=offset) for offset in (0.0, 10.0, 20.0)
        ]
        path = tmp_path / "ensemble.cif"
        write_cif(frames, path)
        assert np.allclose(read_cif(path).xyz[0], [0.0, 0.0, 0.0])
        assert np.allclose(read_cif(path, model=2).xyz[0], [10.0, 0.0, 0.0])
        assert np.allclose(read_cif(path, model=3).xyz[0], [20.0, 0.0, 0.0])

    def test_a_single_structure_is_model_one(self) -> None:
        assert set(model_numbers(structure_to_cif_text(make_structure(3)))) == {"1"}

    def test_frames_that_disagree_on_topology_are_refused(self) -> None:
        with pytest.raises(StructureFileError, match="frames of one file"):
            structure_to_cif_text([make_structure(3), make_structure(4)])

    def test_models_as_frames_false_writes_one_file_each(self, tmp_path: Path) -> None:
        frames = [make_structure(2, all_atom=False, base_offset=o) for o in (0.0, 5.0)]
        write_cif(frames, tmp_path / "conf.cif", models_as_frames=False)
        assert (tmp_path / "conf_1.cif").is_file()
        assert (tmp_path / "conf_2.cif").is_file()
        assert not (tmp_path / "conf.cif").exists()
        assert np.allclose(read_cif(tmp_path / "conf_2.cif").xyz[0], [5.0, 0.0, 0.0])


# ---------------------------------------------------------------------------
# Connectivity: _struct_conn, mirroring the PDB writer's CONECT
# ---------------------------------------------------------------------------


class TestConnectivity:
    def test_ca_only_writes_one_virtual_bond_per_junction(self) -> None:
        text = structure_to_cif_text(make_structure(6, all_atom=False))
        covales = [line for line in text.splitlines() if line.startswith("covale")]
        assert len(covales) == 5, "n_residues - 1 CA-CA bonds"

    def test_all_atom_writes_backbone_and_peptide_bonds(self) -> None:
        text = structure_to_cif_text(make_structure(3))
        covales = [line for line in text.splitlines() if line.startswith("covale")]
        # 3 intra-residue bonds (N-CA, CA-C, C-O) x 3 residues + 2 peptide C-N junctions.
        assert len(covales) == 11

    def test_conect_false_writes_no_struct_conn(self) -> None:
        text = structure_to_cif_text(make_structure(4, all_atom=False), conect=False)
        assert "_struct_conn" not in text

    def test_connectivity_is_written_once_for_all_models(self) -> None:
        text = structure_to_cif_text([make_structure(4, all_atom=False) for _ in range(3)])
        assert text.count("_struct_conn.id") == 1

    def test_a_single_residue_has_no_bonds_and_no_empty_loop(self) -> None:
        text = structure_to_cif_text(make_structure(1, all_atom=False))
        assert "_struct_conn" not in text


# ---------------------------------------------------------------------------
# Feature parity with the PDB writer
# ---------------------------------------------------------------------------


class TestFeatureParity:
    def test_ca_only_writes_exactly_one_ca_per_residue(self, tmp_path: Path) -> None:
        path = tmp_path / "ca.cif"
        write_cif(make_structure(5, all_atom=True), path, ca_only=True)
        structure = read_cif(path)
        assert structure.n_atoms == 5
        assert {str(a) for a in structure.atom_name} == {"CA"}

    def test_annotate_regions_encodes_folded_and_disordered(self, tmp_path: Path) -> None:
        structure = make_structure(6)
        structure.chains[0].domains = [
            Domain(structure, Span(0, 3), DomainKind.FOLDED),
            Domain(structure, Span(3, 6), DomainKind.IDR),
        ]
        path = tmp_path / "annot.cif"
        write_cif(structure, path, annotate_regions=True)
        back = read_cif(path)
        assert back.b_factor[0] == BETA_FOLDED
        assert back.b_factor[-1] == BETA_DISORDERED

    def test_box_writes_a_cell_and_symmetry(self) -> None:
        text = structure_to_cif_text(make_structure(2), box=DEFAULT_BOX_DIMENSIONS)
        assert "_cell.length_a" in text
        assert "_symmetry.space_group_name_H-M" in text

    def test_box_rejects_non_positive_dimensions(self) -> None:
        with pytest.raises(GeometryError, match="positive and finite"):
            structure_to_cif_text(make_structure(2), box=(0.0, 10.0, 10.0))

    def test_seqres_writes_entity_poly_seq_for_written_residues(self) -> None:
        text = structure_to_cif_text(make_structure(4, all_atom=False), seqres=True)
        assert "_entity_poly_seq.mon_id" in text
        mon_rows = [
            line for line in text.splitlines() if line.split()[:2] == ["1", "1"] and "ALA" in line
        ]
        assert mon_rows, "entity 1, residue 1 should be present"


# ---------------------------------------------------------------------------
# What mmCIF can represent that PDB cannot
# ---------------------------------------------------------------------------


class TestMmcifOnlyCapabilities:
    def test_multi_character_chain_ids_survive_without_remapping(self, tmp_path: Path) -> None:
        # The PDB writer must remap "AAA" onto a single-character column; mmCIF keeps it.
        path = tmp_path / "wide_chain.cif"
        write_cif(make_structure(3, all_atom=False, chain_id="AAA"), path)
        assert [c.chain_id for c in read_cif(path).chains] == ["AAA"]

    def test_coordinates_beyond_the_pdb_range_are_written(self, tmp_path: Path) -> None:
        # 12345.678 overflows the PDB %8.3f field; mmCIF has no such limit.
        structure = Structure.from_atom_records(
            xyz=np.array([[12345.678, 0.0, 0.0], [12349.478, 0.0, 0.0]]),
            atom_name=["CA", "CA"],
            element=["C", "C"],
            residue_name=["ALA", "GLY"],
            residue_number=[1, 2],
            chain_id=["A", "A"],
        )
        path = tmp_path / "far.cif"
        write_cif(structure, path, ca_only=True)
        assert np.allclose(read_cif(path).xyz[0], [12345.678, 0.0, 0.0])

    def test_atom_names_with_a_prime_are_quoted_and_round_trip(self, tmp_path: Path) -> None:
        structure = Structure.from_atom_records(
            xyz=np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]),
            atom_name=["O5'", "CA"],
            element=["O", "C"],
            residue_name=["ALA", "ALA"],
            residue_number=[1, 1],
            chain_id=["A", "A"],
        )
        path = tmp_path / "prime.cif"
        write_cif(structure, path)
        assert '"O5\'"' in path.read_text()
        assert "O5'" in [str(a) for a in read_cif(path).atom_name]


# ---------------------------------------------------------------------------
# Failure modes
# ---------------------------------------------------------------------------


class TestFailureModes:
    def test_non_finite_coordinates_are_refused(self, tmp_path: Path) -> None:
        structure = Structure.from_atom_records(
            xyz=np.array([[0.0, 0.0, 0.0], [np.nan, 0.0, 0.0]]),
            atom_name=["CA", "CA"],
            element=["C", "C"],
            residue_name=["ALA", "GLY"],
            residue_number=[1, 2],
            chain_id=["A", "A"],
        )
        with pytest.raises(GeometryError, match="non-finite"):
            write_cif(structure, tmp_path / "nan.cif", ca_only=True)

    def test_empty_sequence_is_refused(self, tmp_path: Path) -> None:
        with pytest.raises(EmptyStructureError):
            write_cif([], tmp_path / "empty.cif")

    def test_a_missing_directory_is_reported(self, tmp_path: Path) -> None:
        with pytest.raises(StructureFileError, match="does not exist"):
            write_cif(make_structure(2), tmp_path / "nope" / "out.cif")

    def test_a_directory_path_is_reported(self, tmp_path: Path) -> None:
        with pytest.raises(StructureFileError, match="existing directory"):
            write_cif(make_structure(2), tmp_path)

    def test_a_bare_relative_filename_works(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.chdir(tmp_path)
        write_cif(make_structure(2, all_atom=False), "bare.cif")
        assert (tmp_path / "bare.cif").is_file()


# ---------------------------------------------------------------------------
# write_structure dispatch
# ---------------------------------------------------------------------------


class TestWriteStructureDispatch:
    def test_cif_extension_writes_mmcif(self, tmp_path: Path) -> None:
        path = tmp_path / "out.cif"
        write_structure(make_structure(3, all_atom=False), path)
        assert path.read_text().startswith("data_")

    def test_mmcif_extension_writes_mmcif(self, tmp_path: Path) -> None:
        path = tmp_path / "out.mmcif"
        write_structure(make_structure(3, all_atom=False), path)
        assert path.read_text().startswith("data_")

    def test_pdb_extension_writes_pdb(self, tmp_path: Path) -> None:
        path = tmp_path / "out.pdb"
        write_structure(make_structure(3, all_atom=False), path)
        first = path.read_text().splitlines()[0]
        assert first.startswith(("ATOM", "MODEL", "SEQRES", "CRYST"))

    def test_dispatch_forwards_ca_only(self, tmp_path: Path) -> None:
        path = tmp_path / "out.cif"
        write_structure(make_structure(4, all_atom=True), path, ca_only=True)
        assert read_cif(path).n_atoms == 4


# ---------------------------------------------------------------------------
# Independent validation with a third-party mmCIF parser, when installed
# ---------------------------------------------------------------------------


class TestExternalParser:
    def test_biotite_reads_every_model_and_atom(self, tmp_path: Path) -> None:
        pdbx = pytest.importorskip("biotite.structure.io.pdbx")
        frames = [
            make_structure(4, all_atom=False, base_offset=offset) for offset in (0.0, 7.0, 14.0)
        ]
        path = tmp_path / "biotite.cif"
        write_cif(frames, path, ca_only=True)
        cif = pdbx.CIFFile.read(str(path))
        stack = pdbx.get_structure(cif)
        assert stack.stack_depth() == 3
        assert stack.array_length() == 4
        assert np.allclose(stack.coord[1][0], [7.0, 0.0, 0.0])

    def test_a_pdb_named_cif_is_not_mmcif(self, tmp_path: Path) -> None:
        # Guards the premise: the old behaviour -- PDB records in a .cif -- really is not
        # readable as mmCIF, which is why this writer had to exist.
        from dodo.io import write_pdb

        path = tmp_path / "fake.cif"
        write_pdb(make_structure(3, all_atom=False), path)
        with pytest.raises(UnsupportedFormatError, match="PDB format, not mmCIF"):
            read_cif(path)
