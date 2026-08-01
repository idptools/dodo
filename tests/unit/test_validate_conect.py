"""Tests for CONECT record validation.

Two rules shape this file.

*Calibrate on real output first.* Every check here is run against files produced by
``dodo.rebuild`` and ``dodo.io.write_pdb`` -- single model and multi-model, ``conect=True``
and ``conect=False``, all-atom and ``ca_only`` -- and against the deposited fixtures in
``tests/data/structures/``. A validator that reports violations on valid output is worse
than none, because it trains people to ignore it.

*Then break it deliberately.* Each defect class has a test that constructs the defect by
mutating real output and asserts the exact kind is reported. The headline case,
``test_the_93_angstrom_bug_is_caught``, reproduces the bug that shipped: a rebuilt alpha
carbon with the residue's other atoms left where they were, written by the real writer.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pytest

import dodo
from dodo.constants import CA_CA_BOND_LENGTH
from dodo.exceptions import EmptyStructureError, MalformedRecordError, StructureFileError
from dodo.io.pdb import decode_hybrid36
from dodo.io.write import _MAX_BONDED_CA_CA_DISTANCE, structure_to_pdb_lines, write_pdb
from dodo.structure import Structure
from dodo.validate.conect import (
    COORDINATE_PRECISION_MARGIN,
    MAX_BOND_LENGTH,
    MIN_BOND_LENGTH,
    ConectReport,
    validate_conect_file,
    validate_conect_lines,
)

FIXTURES = Path(__file__).parent.parent / "data" / "structures"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def make_structure(
    n_residues: int = 5,
    *,
    all_atom: bool = True,
    chain_id: str = "A",
    first_residue_number: int = 1,
    residue_name: str = "ALA",
    x_offset: float = 0.0,
) -> Structure:
    """Build a straight poly-alanine chain with C(i)-N(i+1) at 1.30 A and CA-CA at 3.81 A.

    The same fixture shape ``tests/unit/test_io_write.py`` uses, for the same reason: a
    chain with arbitrary jitter would make the bond-length assertions pass or fail for
    reasons unrelated to the code under test.
    """
    names = ("N", "CA", "C", "O") if all_atom else ("CA",)
    elements = ("N", "C", "C", "O") if all_atom else ("C",)
    offsets = ((-1.0, 0.3), (0.0, 0.0), (1.5, 0.3), (1.6, 1.5)) if all_atom else ((0.0, 0.0),)

    xyz: list[list[float]] = []
    atom_name: list[str] = []
    element: list[str] = []
    for i in range(n_residues):
        base = x_offset + i * CA_CA_BOND_LENGTH
        for name, symbol, (dx, dy) in zip(names, elements, offsets, strict=True):
            xyz.append([base + dx, dy, 0.0])
            atom_name.append(name)
            element.append(symbol)
    n_atoms = len(atom_name)
    return Structure.from_atom_records(
        xyz=np.array(xyz),
        atom_name=atom_name,
        element=element,
        residue_name=[residue_name] * n_atoms,
        residue_number=[
            first_residue_number + i for i in range(n_residues) for _ in range(len(names))
        ],
        chain_id=[chain_id] * n_atoms,
        source="conect-test-fixture",
    )


def atom_record(
    serial: int,
    name: str,
    residue_name: str,
    chain_id: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    element: str = "C",
    *,
    serial_field: str | None = None,
) -> str:
    """Format one ATOM record by column, so hand-built fixtures are really PDB records."""
    serial_text = serial_field if serial_field is not None else f"{serial:5d}"
    # Columns: 7-11 serial, 13-16 name, 17 altLoc, 18-20 resName, 22 chainID, 23-26
    # resSeq, 27 iCode, 31-38 x. The blank altLoc and iCode columns are the ones a
    # hand-written fixture forgets, and forgetting one shifts every field after it.
    return (
        f"ATOM  {serial_text} {name:<4} {residue_name:>3} {chain_id}{residue_number:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}"
    )


def conect_record(origin: int, *partners: int) -> str:
    """Format one CONECT record by column."""
    return "CONECT" + f"{origin:5d}" + "".join(f"{p:5d}" for p in partners)


def move_atom(lines: list[str], serial: int, shift: tuple[float, float, float]) -> list[str]:
    """Return ``lines`` with the atom of the given serial displaced by ``shift``.

    Rewrites the coordinate columns in place, so the result is still a well-formed record
    -- the point is to corrupt the *geometry*, not the formatting.
    """
    out = list(lines)
    for index, line in enumerate(out):
        if not line.startswith(("ATOM  ", "HETATM")) or decode_hybrid36(line[6:11], 5) != serial:
            continue
        coords = [float(line[30:38]) + shift[0], float(line[38:46]) + shift[1]]
        coords.append(float(line[46:54]) + shift[2])
        out[index] = f"{line[:30]}{coords[0]:8.3f}{coords[1]:8.3f}{coords[2]:8.3f}{line[54:]}"
        return out
    raise AssertionError(f"no atom with serial {serial} in these lines")


def set_coordinates(lines: list[str], serial: int, xyz: tuple[float, float, float]) -> list[str]:
    """Return ``lines`` with the atom of the given serial placed exactly at ``xyz``."""
    out = list(lines)
    for index, line in enumerate(out):
        if line.startswith(("ATOM  ", "HETATM")) and decode_hybrid36(line[6:11], 5) == serial:
            out[index] = f"{line[:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"
            return out
    raise AssertionError(f"no atom with serial {serial} in these lines")


def drop_first_conect(lines: list[str]) -> tuple[list[str], str]:
    """Remove the first CONECT record, returning the remaining lines and the record."""
    for index, line in enumerate(lines):
        if line.startswith("CONECT"):
            return lines[:index] + lines[index + 1 :], line
    raise AssertionError("these lines contain no CONECT record")


@pytest.fixture(scope="session")
def rebuilt(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """Real ``dodo.rebuild`` output written every way ``write_pdb`` can write it.

    Session-scoped because the rebuild is the expensive part and every calibration test
    wants the same files. A fixed seed keeps the files identical between runs, so a
    failure here is reproducible.
    """
    directory = tmp_path_factory.mktemp("rebuilt")
    report = dodo.rebuild(FIXTURES / "dnmt3a.pdb", seed=11, n_models=3)
    assert report.ok, report.summary()
    paths = {
        "single": directory / "single.pdb",
        "single_no_conect": directory / "single_no_conect.pdb",
        "multi": directory / "multi.pdb",
        "multi_annotated": directory / "multi_annotated.pdb",
        "ca_only": directory / "ca_only.pdb",
    }
    write_pdb(report.models[0], paths["single"])
    write_pdb(report.models[0], paths["single_no_conect"], conect=False)
    write_pdb(report.models, paths["multi"])
    write_pdb(report.models, paths["multi_annotated"], seqres=True, annotate_regions=True)
    write_pdb(report.models[0], paths["ca_only"], ca_only=True)
    return paths


# ---------------------------------------------------------------------------
# Calibration: real output must be reported clean
# ---------------------------------------------------------------------------


class TestRealOutputIsValid:
    @pytest.mark.parametrize(
        "which", ["single", "single_no_conect", "multi", "multi_annotated", "ca_only"]
    )
    def test_rebuild_output_passes(self, rebuilt: dict[str, Path], which: str) -> None:
        report = validate_conect_file(rebuilt[which])
        assert report.ok, report.describe()

    @pytest.mark.parametrize("which", ["single", "multi", "multi_annotated", "ca_only"])
    def test_rebuild_output_is_complete(self, rebuilt: dict[str, Path], which: str) -> None:
        """Every bond the coordinates imply is declared -- nothing renders as loose dots."""
        report = validate_conect_file(rebuilt[which], require_backbone_conect=True)
        assert report.ok, report.describe()
        assert report.completeness == 1.0
        assert report.n_expected_bonds > 0

    def test_ca_only_output_is_all_virtual_bonds(self, rebuilt: dict[str, Path]) -> None:
        """A CA-only file's bonds are every consecutive CA pair and nothing else."""
        report = validate_conect_file(rebuilt["ca_only"])
        assert report.ok, report.describe()
        assert report.n_bonds == report.n_expected_bonds
        # One CA per residue, so bonds = residues - chains. dnmt3a is a single chain.
        assert report.n_bonds == report.n_atoms - 1

    def test_conect_false_output_is_valid_but_says_what_was_lost(
        self, rebuilt: dict[str, Path]
    ) -> None:
        """``conect=False`` is legitimate output, so it must not produce violations.

        It must also not pass silently: a CA-only region in such a file renders as
        disconnected dots, and the note is the only place a user finds that out.
        """
        report = validate_conect_file(rebuilt["single_no_conect"])
        assert report.ok, report.describe()
        assert report.n_records == 0
        assert report.completeness == 0.0
        assert any("no CONECT records" in note for note in report.notes)

    def test_conect_false_output_can_be_held_to_the_strict_standard(
        self, rebuilt: dict[str, Path]
    ) -> None:
        strict = validate_conect_file(rebuilt["single_no_conect"], require_backbone_conect=True)
        assert not strict.ok
        assert {v.kind for v in strict.violations} == {"missing_bond"}
        assert len(strict.violations) == strict.n_expected_bonds

    @pytest.mark.parametrize("name", sorted(path.name for path in FIXTURES.glob("*.pdb")))
    def test_deposited_and_legacy_fixtures_pass(self, name: str) -> None:
        """Real input files, including v1 output, must be reported clean.

        These carry every oddity the module has to tolerate: unmodelled-residue chain
        breaks, selenomethionine, ligands, multi-model frames, and (in 6kn7) 405
        depositor-written CONECT records for ADP that list four partners on one line.
        """
        report = validate_conect_file(FIXTURES / name)
        assert report.ok, report.describe()

    def test_deposited_ligand_records_are_reciprocal(self) -> None:
        """The claim in the module docstring, measured rather than asserted from memory."""
        report = validate_conect_file(FIXTURES / "6kn7.pdb", require_reciprocal=True)
        assert report.ok, report.describe()
        assert report.n_bonds == 435
        assert report.n_non_reciprocal == 0

    def test_dodo_writer_output_is_not_reciprocal(self, rebuilt: dict[str, Path]) -> None:
        """The other half of that claim: DODO records each bond once, from the lower serial.

        This is why ``require_reciprocal`` defaults to off. If the writer ever starts
        emitting both directions, this test fails and the default should be revisited.
        """
        report = validate_conect_file(rebuilt["single"])
        assert report.n_non_reciprocal == report.n_bonds
        strict = validate_conect_file(rebuilt["single"], require_reciprocal=True)
        assert len(strict.of_kind("non_reciprocal")) == report.n_bonds

    @pytest.mark.slow
    def test_all_atom_rebuild_output_passes(self, tmp_path: Path) -> None:
        """``all_atom=True`` places N, C and O in rebuilt regions, so it changes the bonds.

        The intra-residue records and the C-N junctions it enables are exactly where a
        connectivity defect would appear next, since placing real backbone atoms in
        rebuilt regions is the priority after this suite.
        """
        report = dodo.rebuild(FIXTURES / "dnmt3a.pdb", seed=5, all_atom=True, sidechains=True)
        assert report.ok, report.summary()
        path = tmp_path / "all_atom.pdb"
        write_pdb(report.models[0], path)
        checked = validate_conect_file(path, require_backbone_conect=True)
        assert checked.ok, checked.describe()
        assert checked.completeness == 1.0

    def test_insertion_codes_do_not_confuse_the_junctions(self) -> None:
        """Residues 10 and 10A share a sequence number and are still bonded neighbours.

        Connectivity is decided on distance for this reason among others: a numbering test
        both misses this bond and would merge the two residues into one.
        """
        chain = make_structure(3, all_atom=False)
        structure = Structure.from_atom_records(
            xyz=chain.xyz,
            atom_name=["CA"] * 3,
            element=["C"] * 3,
            residue_name=["ALA"] * 3,
            residue_number=[10, 10, 11],
            insertion_code=["", "A", ""],
            chain_id=["A"] * 3,
        )
        report = validate_conect_lines(structure_to_pdb_lines(structure))
        assert report.ok, report.describe()
        assert report.n_expected_bonds == 2
        assert report.completeness == 1.0
        broken, _ = drop_first_conect(structure_to_pdb_lines(structure))
        missing = validate_conect_lines(broken).of_kind("missing_bond")
        assert len(missing) == 1
        assert "A:ALA10 CA" in missing[0].message
        assert "A:ALA10A CA" in missing[0].message

    def test_a_full_em_assembly_through_dodos_writer(self, tmp_path: Path) -> None:
        """61,511 atoms and 31,093 bonds, which is the real performance target."""
        structure = dodo.read_structure(FIXTURES / "6kn7.pdb")
        path = tmp_path / "6kn7_dodo.pdb"
        write_pdb(structure, path)
        start = time.perf_counter()
        report = validate_conect_file(path, require_backbone_conect=True)
        elapsed = time.perf_counter() - start
        assert report.ok, report.describe()
        assert report.n_atoms == 61511
        assert report.completeness == 1.0
        # Generous: measured at 0.33 s. The bound exists to catch a quadratic regression,
        # not to police a few milliseconds on a loaded CI box.
        assert elapsed < 10.0, f"validation took {elapsed:.1f} s"


# ---------------------------------------------------------------------------
# The defect this module was written for
# ---------------------------------------------------------------------------


class TestSpuriousBonds:
    def test_the_93_angstrom_bug_is_caught(self, tmp_path: Path) -> None:
        """A rebuilt CA with the residue's other atoms left behind.

        The bug that shipped, reproduced through the real writer rather than by editing
        text: ``dodo.io.write._bonds`` emits intra-residue N-CA, CA-C and C-O with no
        distance test at all, so moving one alpha carbon is enough to make it write a
        CONECT record spanning 93 A. In a viewer that is a straight line through the
        structure; nothing in the suite noticed it before this check existed.
        """
        structure = make_structure(6)
        moved = int(np.flatnonzero(structure.atom_name == "CA")[2])
        structure.xyz[moved] += (93.0, 0.0, 0.0)
        path = tmp_path / "left_behind.pdb"
        write_pdb(structure, path)

        report = validate_conect_file(path)
        assert not report.ok
        long_bonds = report.of_kind("bond_too_long")
        # Exactly the two intra-residue records that span the gap, N-CA and CA-C. The
        # junctions to the neighbouring residues produce nothing, in either direction:
        # the writer decides a junction from its CA-CA separation, which is now 97 A, so
        # it emits no record there and none is expected. That the chain is broken is a
        # geometry finding (validate_ca_trace), not a connectivity one.
        assert len(long_bonds) == 2
        assert {v.atoms for v in long_bonds} == {
            ("A:ALA3 N", "A:ALA3 CA"),
            ("A:ALA3 CA", "A:ALA3 C"),
        }
        assert not report.of_kind("missing_bond")
        worst = max(long_bonds, key=lambda v: v.distance or 0.0)
        assert worst.distance is not None
        assert 92.0 < worst.distance < 95.0
        assert f"{worst.distance:.3f}" in worst.message
        assert f"{MAX_BOND_LENGTH:.3f}" in worst.message

    def test_a_bond_at_the_threshold_is_not_flagged(self) -> None:
        """The writer bonds CA-CA out to 1.3x the ideal, so the validator must accept it.

        Measured on real ``ca_only=True`` output, folded-domain CA-CA bonds run 2.978 to
        3.914 A: cis-proline at the bottom, strained crystal geometry at the top. A
        threshold at the ideal bond length plus its tolerance would flag deposited
        geometry DODO merely passed through.
        """
        for separation in (2.900, 3.810, 3.914, MAX_BOND_LENGTH - 0.001):
            lines = [
                atom_record(1, "CA", "PRO", "A", 1, 0.0, 0.0, 0.0),
                atom_record(2, "CA", "PRO", "A", 2, separation, 0.0, 0.0),
                conect_record(1, 2),
            ]
            report = validate_conect_lines(lines)
            assert report.ok, f"{separation} A: {report.describe()}"

    def test_just_past_the_threshold_is_flagged(self) -> None:
        lines = [
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 2, MAX_BOND_LENGTH + 0.01, 0.0, 0.0),
            conect_record(1, 2),
        ]
        assert len(validate_conect_lines(lines).of_kind("bond_too_long")) == 1

    def test_coincident_atoms_are_caught(self) -> None:
        """Two atoms at one point, which is what a failed placement leaves behind."""
        lines = [
            atom_record(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            conect_record(1, 2),
        ]
        report = validate_conect_lines(lines)
        short = report.of_kind("bond_too_short")
        assert len(short) == 1
        assert short[0].distance == 0.0
        assert "same point" in short[0].message
        assert "A:ALA1 N" in short[0].message and "A:ALA1 CA" in short[0].message

    def test_a_hydrogen_length_bond_is_not_called_too_short(self) -> None:
        """DODO passes input hydrogens through untouched; 0.97 A is a real bond."""
        lines = [
            atom_record(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, element="N"),
            atom_record(2, "H", "ALA", "A", 1, 0.97, 0.0, 0.0, element="H"),
            conect_record(1, 2),
        ]
        assert validate_conect_lines(lines).ok

    def test_a_short_bond_from_real_output_is_caught(self, rebuilt: dict[str, Path]) -> None:
        """The same defect, injected into a real file rather than a two-atom fixture."""
        lines = rebuilt["single"].read_text().splitlines()
        first_bond = next(line for line in lines if line.startswith("CONECT"))
        origin, partner = int(first_bond[6:11]), int(first_bond[11:16])
        origin_line = next(
            line for line in lines if line.startswith("ATOM  ") and int(line[6:11]) == origin
        )
        here = (float(origin_line[30:38]), float(origin_line[38:46]), float(origin_line[46:54]))
        report = validate_conect_lines(set_coordinates(lines, partner, here))
        assert len(report.of_kind("bond_too_short")) == 1


# ---------------------------------------------------------------------------
# Missing bonds -- the disconnected-dots failure
# ---------------------------------------------------------------------------


class TestMissingBonds:
    def test_a_deleted_record_is_reported(self, rebuilt: dict[str, Path]) -> None:
        lines = rebuilt["ca_only"].read_text().splitlines()
        remaining, removed = drop_first_conect(lines)
        report = validate_conect_lines(remaining, source="ca_only.pdb")
        missing = report.of_kind("missing_bond")
        assert len(missing) == 1
        assert missing[0].serials == (int(removed[6:11]), int(removed[11:16]))
        assert report.n_expected_present == report.n_expected_bonds - 1
        assert report.completeness < 1.0
        assert "disconnected dots" in missing[0].message

    def test_the_message_names_both_residues_and_the_distance(
        self, rebuilt: dict[str, Path]
    ) -> None:
        lines, _ = drop_first_conect(rebuilt["ca_only"].read_text().splitlines())
        message = validate_conect_lines(lines).of_kind("missing_bond")[0].message
        assert "A:MET1 CA" in message
        assert "A:PRO2 CA" in message
        assert "3.8" in message
        assert message.endswith(".")

    def test_a_genuine_chain_break_is_not_a_missing_bond(self) -> None:
        """Unmodelled residues leave a gap that must not be bonded or reported.

        This is the case that separates a useful check from a noise generator: 6kn7 has
        several such gaps, and an implementation keyed on residue numbering rather than
        on distance reports every one of them.
        """
        first = make_structure(3, all_atom=False)
        second = make_structure(3, all_atom=False, first_residue_number=40, x_offset=25.0)
        joined = Structure.from_atom_records(
            xyz=np.vstack([first.xyz, second.xyz]),
            atom_name=["CA"] * 6,
            element=["C"] * 6,
            residue_name=["ALA"] * 6,
            residue_number=[1, 2, 3, 40, 41, 42],
            chain_id=["A"] * 6,
        )
        report = validate_conect_lines(structure_to_pdb_lines(joined))
        assert report.ok, report.describe()
        # Two bonds per fragment, and nothing across the gap.
        assert report.n_bonds == 4
        assert report.n_expected_bonds == 4

    def test_a_ca_only_region_next_to_an_all_atom_region_is_clean(self) -> None:
        """The mixed file DODO actually produces: no "residue has no N" findings.

        A rebuilt IDR is alpha carbons only and a folded domain is all-atom, in one file,
        with a real chain connection between them. A validator that expects a full
        backbone everywhere reports thousands of useless findings here.
        """
        folded = make_structure(4, all_atom=True)
        rebuilt_region = make_structure(
            4, all_atom=False, first_residue_number=5, x_offset=4 * CA_CA_BOND_LENGTH
        )
        mixed = Structure.from_atom_records(
            xyz=np.vstack([folded.xyz, rebuilt_region.xyz]),
            atom_name=[*folded.atom_name.tolist(), *rebuilt_region.atom_name.tolist()],
            element=[*folded.element.tolist(), *rebuilt_region.element.tolist()],
            residue_name=["ALA"] * (folded.n_atoms + rebuilt_region.n_atoms),
            residue_number=[*np.repeat([1, 2, 3, 4], 4).tolist(), 5, 6, 7, 8],
            chain_id=["A"] * (folded.n_atoms + rebuilt_region.n_atoms),
        )
        report = validate_conect_lines(structure_to_pdb_lines(mixed))
        assert report.ok, report.describe()
        assert report.completeness == 1.0
        # 3 intra-residue bonds x 4 all-atom residues, 3 peptide bonds inside the folded
        # stretch, then the CA-CA bond across the boundary and 3 inside the CA-only run.
        assert report.n_expected_bonds == 12 + 3 + 1 + 3

    def test_a_ter_record_ends_a_chain_even_when_the_id_repeats(self) -> None:
        """Two runs of chain A separated by TER are two chains, so nothing bridges them.

        Deposited files do this routinely -- 6kn7 parks its ADP ligands in chain A after
        the TER -- and ``dodo.io.pdb`` splits on TER for the same reason. The record name
        is matched with trailing blanks stripped, because a bare ``TER`` with no serial is
        common and missing it would make the two residues look like bonded neighbours.
        """
        lines = [
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0),
            "TER",
            atom_record(3, "CA", "ALA", "A", 3, 2 * CA_CA_BOND_LENGTH, 0.0, 0.0),
            conect_record(1, 2),
        ]
        report = validate_conect_lines(lines)
        assert report.ok, report.describe()
        # Only the bond inside the first run. Residues 2 and 3 are 3.81 A apart and would
        # be a junction but for the TER between them.
        assert report.n_expected_bonds == 1
        assert report.completeness == 1.0

    def test_a_ligand_is_not_expected_to_have_backbone_bonds(self) -> None:
        """ADP has atoms named C, N and O; it is not a residue with a broken backbone."""
        protein = make_structure(3, all_atom=True)
        ligand_names = ["N1", "C2", "N3", "C4", "O5"]
        combined = Structure.from_atom_records(
            xyz=np.vstack([protein.xyz, protein.xyz[:5] + np.array([40.0, 0.0, 0.0])]),
            atom_name=[*protein.atom_name.tolist(), *ligand_names],
            element=[*protein.element.tolist(), "N", "C", "N", "C", "O"],
            residue_name=["ALA"] * protein.n_atoms + ["ADP"] * 5,
            residue_number=[*np.repeat([1, 2, 3], 4).tolist(), *[401] * 5],
            chain_id=["A"] * (protein.n_atoms + 5),
        )
        report = validate_conect_lines(structure_to_pdb_lines(combined))
        assert report.ok, report.describe()
        # 3 x 3 intra-residue plus 2 peptide bonds. The ligand contributes none.
        assert report.n_expected_bonds == 11

    def test_a_deposition_is_not_buried_in_missing_bond_findings(self) -> None:
        """6kn7 declares 435 ligand bonds and no backbone bonds, which is normal.

        wwPDB depositions never list backbone connectivity -- readers infer it from the
        residue names -- so reporting 31,093 missing bonds would be reporting the format.
        """
        report = validate_conect_file(FIXTURES / "6kn7.pdb")
        assert report.ok
        assert report.n_expected_bonds > 30000
        assert report.n_expected_present == 0
        assert any("deposition" in note for note in report.notes)
        strict = validate_conect_file(FIXTURES / "6kn7.pdb", require_backbone_conect=True)
        assert len(strict.of_kind("missing_bond")) == report.n_expected_bonds


# ---------------------------------------------------------------------------
# Illogical records
# ---------------------------------------------------------------------------


class TestIllogicalRecords:
    def test_unknown_serial(self, rebuilt: dict[str, Path]) -> None:
        lines = [*rebuilt["single"].read_text().splitlines(), conect_record(1, 99998)]
        report = validate_conect_lines(lines)
        unknown = report.of_kind("unknown_serial")
        assert len(unknown) == 1
        assert unknown[0].serials == (99998,)
        assert "does not exist" in unknown[0].message

    def test_a_reference_to_a_ter_serial_says_so(self) -> None:
        """TER consumes an atom serial, so this is the mistake a writer actually makes."""
        lines = [
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0),
            "TER       3      ALA A   2",
            atom_record(4, "CA", "ALA", "B", 1, 40.0, 0.0, 0.0),
            conect_record(1, 2),
            conect_record(2, 3),
        ]
        report = validate_conect_lines(lines)
        unknown = report.of_kind("unknown_serial")
        assert len(unknown) == 1
        assert unknown[0].serials == (3,)
        assert "TER record" in unknown[0].message

    def test_self_bond(self) -> None:
        lines = [
            atom_record(1, "CA", "GLU", "A", 142, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 143, CA_CA_BOND_LENGTH, 0.0, 0.0),
            conect_record(1, 2),
            conect_record(2, 2),
        ]
        report = validate_conect_lines(lines)
        self_bonds = report.of_kind("self_bond")
        assert len(self_bonds) == 1
        assert self_bonds[0].serials == (2, 2)
        assert "to itself" in self_bonds[0].message

    def test_duplicate_bond_from_the_same_origin(self, rebuilt: dict[str, Path]) -> None:
        lines = rebuilt["ca_only"].read_text().splitlines()
        first = next(line for line in lines if line.startswith("CONECT"))
        report = validate_conect_lines([*lines, first])
        duplicates = report.of_kind("duplicate_bond")
        assert len(duplicates) == 1
        assert duplicates[0].serials == (int(first[6:11]), int(first[11:16]))
        assert "2 times from the same atom" in duplicates[0].message
        assert "on lines " in duplicates[0].message

    def test_the_same_partner_twice_on_one_record(self) -> None:
        """One record listing a partner twice is one duplicate, reported against one line."""
        lines = [
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0),
            conect_record(1, 2, 2),
        ]
        duplicates = validate_conect_lines(lines).of_kind("duplicate_bond")
        assert len(duplicates) == 1
        assert "on line 3." in duplicates[0].message

    def test_a_reciprocal_pair_is_not_a_duplicate(self) -> None:
        """The distinction that matters: both directions is the PDB convention, not a bug."""
        lines = [
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0),
            conect_record(1, 2),
            conect_record(2, 1),
        ]
        report = validate_conect_lines(lines, require_reciprocal=True)
        assert report.ok, report.describe()
        assert report.n_bonds == 1
        assert report.n_non_reciprocal == 0

    def test_a_continuation_record_is_not_a_duplicate(self) -> None:
        """An atom with more than four partners takes a second record. That is legal."""
        centre = atom_record(1, "C", "ALA", "A", 1, 0.0, 0.0, 0.0)
        neighbours = [
            atom_record(serial, "C", "ALA", "A", 1, *offset)
            for serial, offset in enumerate(
                [
                    (1.5, 0.0, 0.0),
                    (-1.5, 0.0, 0.0),
                    (0.0, 1.5, 0.0),
                    (0.0, -1.5, 0.0),
                    (0.0, 0.0, 1.5),
                ],
                start=2,
            )
        ]
        lines = [centre, *neighbours, conect_record(1, 2, 3, 4, 5), conect_record(1, 6)]
        report = validate_conect_lines(lines)
        assert report.ok, report.describe()
        assert report.n_bonds == 5
        assert not report.of_kind("duplicate_bond")

    def test_non_reciprocal_is_reported_only_when_asked(self) -> None:
        lines = [
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0),
            atom_record(3, "CA", "ALA", "A", 3, 2 * CA_CA_BOND_LENGTH, 0.0, 0.0),
            conect_record(1, 2),
            conect_record(2, 1, 3),
            # 3 does not list 2 back.
        ]
        lenient = validate_conect_lines(lines)
        assert lenient.ok
        assert lenient.n_non_reciprocal == 1
        assert any("one of their two atoms" in note for note in lenient.notes)

        strict = validate_conect_lines(lines, require_reciprocal=True)
        violations = strict.of_kind("non_reciprocal")
        assert len(violations) == 1
        assert violations[0].serials == (2, 3)
        assert "does not list" in violations[0].message

    def test_interchain_bond(self) -> None:
        """DODO never bonds across a chain break, so a record that does is a defect."""
        first = make_structure(3, all_atom=False, chain_id="A")
        second = make_structure(3, all_atom=False, chain_id="B", x_offset=25.0)
        joined = Structure.from_atom_records(
            xyz=np.vstack([first.xyz, second.xyz]),
            atom_name=["CA"] * 6,
            element=["C"] * 6,
            residue_name=["ALA"] * 6,
            residue_number=[1, 2, 3, 1, 2, 3],
            chain_id=["A", "A", "A", "B", "B", "B"],
        )
        lines = structure_to_pdb_lines(joined)
        assert validate_conect_lines(lines).ok
        # Serial 4 is the TER record between the chains, so chain B starts at 5.
        report = validate_conect_lines([*lines, conect_record(3, 5)])
        interchain = report.of_kind("interchain_bond")
        assert len(interchain) == 1
        assert "different chains" in interchain[0].message
        assert "A:ALA3 CA" in interchain[0].message and "B:ALA1 CA" in interchain[0].message

    def test_ambiguous_serial(self) -> None:
        """Two atoms sharing a serial make every record that names it meaningless."""
        lines = [
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(1, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0),
            conect_record(1, 1),
        ]
        report = validate_conect_lines(lines)
        assert {v.kind for v in report.violations} == {"self_bond", "ambiguous_serial"}
        assert "2 different atoms" in report.of_kind("ambiguous_serial")[0].message

    @pytest.mark.parametrize(
        ("record", "expected"),
        [
            ("CONECT", "names no origin atom"),
            ("CONECT    1", "no partners"),
            ("CONECT   ?!    2", "cannot be decoded"),
            ("CONECT    1   ?!", "cannot be decoded"),
        ],
    )
    def test_malformed_records_are_reported_not_raised(self, record: str, expected: str) -> None:
        """A garbled CONECT record is the thing under test, so it is a finding."""
        lines = [
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0),
            conect_record(1, 2),
            record,
        ]
        report = validate_conect_lines(lines)
        malformed = report.of_kind("malformed_record")
        assert len(malformed) == 1
        assert expected in malformed[0].message
        assert malformed[0].line_number == 4


# ---------------------------------------------------------------------------
# Multi-model files
# ---------------------------------------------------------------------------


class TestMultiModel:
    def test_a_bond_broken_only_in_a_later_frame_is_caught(self, rebuilt: dict[str, Path]) -> None:
        """One set of CONECT records has to describe every frame.

        Frame 1 is untouched; an atom in frame 3 is displaced. Checking only the first
        frame -- the obvious implementation -- reports this file as valid.
        """
        lines = rebuilt["multi"].read_text().splitlines()
        third_frame = lines.index("MODEL        3")
        head, tail = lines[:third_frame], lines[third_frame:]
        target = int(next(line for line in tail if line.startswith("ATOM  "))[6:11])
        report = validate_conect_lines([*head, *move_atom(tail, target, (60.0, 0.0, 0.0))])
        long_bonds = report.of_kind("bond_too_long")
        assert long_bonds, report.describe()
        assert all(v.model == "model 3" for v in long_bonds)
        assert "in model 3" in long_bonds[0].message

    def test_frames_that_disagree_on_ter_placement_are_caught(self) -> None:
        """The failure an earlier review found: serials that mean different atoms.

        A TER record consumes an atom serial, so a frame with one fewer TER numbers every
        atom after it differently. The file loads without complaint in MDTraj and its
        CONECT records bond the wrong atoms in the second frame, which is the worst
        possible output. ``write_pdb`` refuses to write this; the validator has to be able
        to recognise it in a file that already exists.
        """
        model_1 = [
            "MODEL        1",
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0),
            "TER       3      ALA A   2",
            atom_record(4, "CA", "ALA", "B", 1, 40.0, 0.0, 0.0),
            atom_record(5, "CA", "ALA", "B", 2, 40.0 + CA_CA_BOND_LENGTH, 0.0, 0.0),
            "ENDMDL",
        ]
        model_2 = [
            "MODEL        2",
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0),
            # No TER, so chain B's alpha carbons become serials 3 and 4.
            atom_record(3, "CA", "ALA", "B", 1, 40.0, 0.0, 0.0),
            atom_record(4, "CA", "ALA", "B", 2, 40.0 + CA_CA_BOND_LENGTH, 0.0, 0.0),
            "ENDMDL",
        ]
        lines = [*model_1, *model_2, conect_record(1, 2), conect_record(4, 5), "END"]
        report = validate_conect_lines(lines)
        mismatched = report.of_kind("frame_mismatch")
        assert {v.serials[0] for v in mismatched} == {4, 5}
        assert any("TER" in v.message for v in mismatched)
        # Serial 4 is chain B residue 1 in model 1 and chain B residue 2 in model 2.
        assert any("different atoms in different frames" in v.message for v in mismatched)
        # Serial 5 does not exist at all in model 2.
        assert any("is not an atom in model 2" in v.message for v in mismatched)

    def test_frames_are_counted_and_named(self, rebuilt: dict[str, Path]) -> None:
        report = validate_conect_file(rebuilt["multi"])
        assert report.n_models == 3
        single = validate_conect_file(rebuilt["single"])
        assert single.n_models == 1
        assert report.n_atoms == single.n_atoms
        assert report.n_bonds == single.n_bonds


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


class TestParsing:
    def test_columns_not_split(self) -> None:
        """Fields are read by column, so merged and missing whitespace cannot shift them.

        ``split()`` gets both of these wrong: the negative coordinates run together into
        one token, and the blank altLoc and insertion-code columns remove tokens that a
        positional reader would then take from the wrong place.
        """
        lines = [
            "ATOM      1  CA  GLU A 142    -112.345-106.789 -10.500  1.00  0.00           C",
            "ATOM      2  CA  ALA A 143    -100.000-100.000  -8.000  1.00  0.00           C",
            conect_record(1, 2),
        ]
        report = validate_conect_lines(lines)
        assert report.n_atoms == 2
        assert len(report.of_kind("bond_too_long")) == 1
        assert "A:GLU142 CA" in report.violations[0].message

    def test_hybrid36_serials(self) -> None:
        """Files over 99,999 atoms switch the serial column to base 36. Both ends decode.

        ``A0000`` is 100,000 and ``A0001`` is 100,001, and a CONECT record written by the
        same encoding has to resolve to those atoms rather than be reported as unknown.
        """
        lines = [
            atom_record(0, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, serial_field="A0000"),
            atom_record(0, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0, serial_field="A0001"),
            "CONECTA0000A0001",
        ]
        report = validate_conect_lines(lines)
        assert report.ok, report.describe()
        assert report.n_bonds == 1
        assert report.completeness == 1.0
        broken = validate_conect_lines(
            [lines[0], move_atom([lines[1]], 100001, (50.0, 0.0, 0.0))[0], lines[2]]
        )
        assert broken.of_kind("bond_too_long")[0].serials == (100000, 100001)

    def test_records_padded_to_eighty_columns(self) -> None:
        """Deposited files pad; the padding must not read as extra partner fields."""
        lines = [
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0),
            conect_record(1, 2).ljust(80),
        ]
        report = validate_conect_lines(lines)
        assert report.ok, report.describe()
        assert report.n_bonds == 1

    def test_extra_partner_fields_are_read_and_noted(self) -> None:
        """A fifth partner sits in columns 32-36, which format 2.2 used for hydrogen bonds."""
        lines = [
            atom_record(1, "C", "ALA", "A", 1, 0.0, 0.0, 0.0),
            *[
                atom_record(serial, "C", "ALA", "A", 1, *offset)
                for serial, offset in enumerate(
                    [
                        (1.5, 0.0, 0.0),
                        (-1.5, 0.0, 0.0),
                        (0.0, 1.5, 0.0),
                        (0.0, -1.5, 0.0),
                        (0.0, 0.0, 1.5),
                    ],
                    start=2,
                )
            ],
            conect_record(1, 2, 3, 4, 5, 6),
        ]
        report = validate_conect_lines(lines)
        assert report.ok, report.describe()
        assert report.n_bonds == 5
        assert any("beyond the four" in note for note in report.notes)

    def test_line_endings_and_whitespace(self) -> None:
        lines = [
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0) + "\n",
            atom_record(2, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0) + "\r\n",
            conect_record(1, 2) + "\n",
        ]
        assert validate_conect_lines(lines).ok

    def test_conect_records_inside_a_model_block_are_noted(self) -> None:
        lines = [
            "MODEL        1",
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0),
            conect_record(1, 2),
            "ENDMDL",
        ]
        report = validate_conect_lines(lines)
        assert report.ok, report.describe()
        assert any("inside a" in note for note in report.notes)

    def test_an_empty_file_raises(self) -> None:
        with pytest.raises(EmptyStructureError, match="no ATOM or HETATM"):
            validate_conect_lines(["HEADER    NOTHING HERE", "END"])

    def test_a_garbled_coordinate_raises(self) -> None:
        """Unusable reference data raises; a bad CONECT record does not."""
        with pytest.raises(MalformedRecordError, match="coordinate columns"):
            validate_conect_lines(
                ["ATOM      1  CA  ALA A   1        n/a   0.000   0.000  1.00  0.00           C"]
            )

    def test_a_garbled_atom_serial_raises(self) -> None:
        with pytest.raises(MalformedRecordError, match="atom serial"):
            validate_conect_lines(
                ["ATOM     ?!  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C"]
            )

    def test_a_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(StructureFileError, match="Cannot read"):
            validate_conect_file(tmp_path / "absent.pdb")

    @pytest.mark.parametrize(
        ("kwargs", "match"),
        [
            ({"max_bond_length": 0.0}, "must be finite and positive"),
            ({"min_bond_length": float("nan")}, "must be finite and positive"),
            ({"min_bond_length": 5.0, "max_bond_length": 4.0}, "must be below"),
        ],
    )
    def test_nonsense_thresholds_raise(self, kwargs: dict[str, float], match: str) -> None:
        lines = [
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0),
        ]
        with pytest.raises(ValueError, match=match):
            validate_conect_lines(lines, **kwargs)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Report shape and thresholds
# ---------------------------------------------------------------------------


class TestReport:
    def test_the_threshold_tracks_the_writer(self) -> None:
        """The validator's ceiling is the writer's own connectivity cutoff, plus rounding.

        If someone retunes the writer's slack, this fails and both numbers move together.
        Without it the two drift and the validator starts flagging records the writer
        deliberately emitted, or stops flagging ones it never would.
        """
        writers_cutoff = _MAX_BONDED_CA_CA_DISTANCE + COORDINATE_PRECISION_MARGIN
        assert abs(MAX_BOND_LENGTH - writers_cutoff) < 1e-12
        assert MAX_BOND_LENGTH > CA_CA_BOND_LENGTH
        assert MIN_BOND_LENGTH < 1.0

    def test_the_margin_covers_pdb_coordinate_rounding(self) -> None:
        """0.001 A per coordinate, so at most sqrt(3) * 0.001 A on a distance."""
        assert np.sqrt(3.0) * 0.001 < COORDINATE_PRECISION_MARGIN

    def test_ok_and_bool_agree(self, rebuilt: dict[str, Path]) -> None:
        report = validate_conect_file(rebuilt["single"])
        assert report.ok and bool(report)
        broken = validate_conect_lines(
            [*rebuilt["single"].read_text().splitlines(), conect_record(1, 1)]
        )
        assert not broken.ok and not bool(broken)

    def test_raise_if_invalid(self, rebuilt: dict[str, Path]) -> None:
        validate_conect_file(rebuilt["single"]).raise_if_invalid()
        broken = validate_conect_lines(
            [*rebuilt["single"].read_text().splitlines(), conect_record(1, 1)]
        )
        with pytest.raises(StructureFileError, match="Invalid CONECT records"):
            broken.raise_if_invalid()

    def test_describe_truncates(self, rebuilt: dict[str, Path]) -> None:
        report = validate_conect_file(rebuilt["single_no_conect"], require_backbone_conect=True)
        assert len(report.violations) > 20
        described = report.describe(max_violations=5)
        assert described.count("\n  - ") == 6  # five violations plus the "and N more" line
        assert "more violation(s)" in described

    def test_summary_reports_counts_and_completeness(self, rebuilt: dict[str, Path]) -> None:
        report = validate_conect_file(rebuilt["multi"])
        summary = report.summary()
        assert "multi.pdb" in summary
        assert "3 models" in summary
        assert f"{report.n_bonds} bond(s)" in summary
        assert "100.0%" in summary
        assert "valid" in summary

    def test_empty_expectations_are_complete_not_zero(self) -> None:
        """A ligand-only file expects no backbone bonds; that is 100%, not 0%."""
        lines = [
            atom_record(1, "O", "HOH", "A", 1, 0.0, 0.0, 0.0, element="O"),
            atom_record(2, "O", "HOH", "A", 2, 5.0, 0.0, 0.0, element="O"),
        ]
        report = validate_conect_lines(lines)
        assert report.n_expected_bonds == 0
        assert report.completeness == 1.0
        assert report.ok

    def test_every_violation_is_actionable(self, rebuilt: dict[str, Path]) -> None:
        """Rule for this suite: a finding a user cannot act on is noise.

        Every message must be a complete sentence, name the residues involved in
        ``chain:RESNUM`` form, and quote the measured value where there is one. Checked
        over one file corrupted in every way at once, so a new kind cannot be added
        without a message that meets the standard.
        """
        lines = rebuilt["single"].read_text().splitlines()
        # Drop the bond between serials 1 and 2, then displace serial 5 -- whose own two
        # records survive and become 75 A long. Displacing an atom whose record was just
        # dropped would produce neither finding, since nothing then declares that bond
        # and nothing expects it either.
        lines, _ = drop_first_conect(lines)
        lines = move_atom(lines, 5, (75.0, 0.0, 0.0))
        lines = [
            *lines,
            conect_record(9, 9),
            conect_record(9, 99998),
            "CONECT",
            *[line for line in lines if line.startswith("CONECT")][:1],
        ]
        report = validate_conect_lines(lines, source="corrupted.pdb")
        kinds = set(report.kinds)
        assert kinds == {
            "bond_too_long",
            "duplicate_bond",
            "malformed_record",
            "missing_bond",
            "self_bond",
            "unknown_serial",
        }, report.describe(40)
        for violation in report.violations:
            assert violation.message.endswith("."), violation.message
            assert violation.message[0].isupper() or violation.message.startswith(("A:", "B:"))
            if violation.kind != "malformed_record":
                assert ":" in violation.message, violation.message
            if violation.distance is not None:
                assert f"{violation.distance:.3f}"[:4] in violation.message, violation.message

    def test_report_repr(self) -> None:
        lines = [
            atom_record(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            atom_record(2, "CA", "ALA", "A", 2, CA_CA_BOND_LENGTH, 0.0, 0.0),
            conect_record(1, 2),
        ]
        report = validate_conect_lines(lines)
        assert isinstance(report, ConectReport)
        assert repr(report) == "ConectReport(1 records, 1 bonds, valid)"
