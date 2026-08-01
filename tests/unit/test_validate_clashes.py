"""Tests for steric clash validation.

Two halves. The first is synthetic and exhaustive: one hand-placed structure per defect
class and per exclusion, so that a regression names the thing it broke. The second is
calibration against the six real fixtures in ``tests/data/structures/`` and against fresh
``dodo.rebuild()`` output, which is the only way to know whether the limits are set
somewhere useful. A clash validator that has never been pointed at a real structure is a
guess, and the specific way it fails -- thousands of findings on perfectly good input --
is invisible in synthetic tests because synthetic structures have no hydrogen bonds, no
disulfides, no cis-prolines and no unresolved side chains.

Every calibration number asserted here was measured, and the measurement is in the
assertion message so a future reader can tell drift from breakage.
"""

from __future__ import annotations

import time
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np
import pytest

from dodo.constants import CA_CLASH_DISTANCE, CLASH_EXCLUDE_WITHIN_RESIDUES
from dodo.exceptions import EmptyStructureError, GeometryError
from dodo.io import read_structure
from dodo.structure import Structure
from dodo.validate.clashes import (
    ALLOWED_OVERLAP,
    DISULFIDE_MAX_DISTANCE,
    HBOND_ALLOWED_OVERLAP,
    VDW_RADII,
    ca_only_residue_mask,
    contact_limit,
    covalent_bond_pairs,
    validate_clashes,
)
from dodo.validate.reference import RESIDUE_BONDS

DATA = Path(__file__).resolve().parents[1] / "data" / "structures"

# ---------------------------------------------------------------------------
# Synthetic structure builder
# ---------------------------------------------------------------------------

#: One residue's worth of atoms: name -> coordinate.
AtomSpec = Mapping[str, Sequence[float]]

#: One residue: (residue name, atoms, chain id, residue number).
ResidueSpec = tuple[str, AtomSpec, str, int]


def build(residues: Sequence[ResidueSpec]) -> Structure:
    """Assemble a structure from explicitly placed atoms.

    Every coordinate in these tests is written out by hand rather than generated, because
    the whole point is to control exactly which pairs are close and by how much.
    """
    xyz: list[Sequence[float]] = []
    names: list[str] = []
    elements: list[str] = []
    residue_names: list[str] = []
    residue_numbers: list[int] = []
    chain_ids: list[str] = []
    for residue_name, atoms, chain_id, residue_number in residues:
        for atom_name, coordinate in atoms.items():
            xyz.append(coordinate)
            names.append(atom_name)
            # Element from the first character of the atom name, which is how the PDB
            # format works for every atom these tests use.
            elements.append("SE" if atom_name == "SE" else atom_name[0])
            residue_names.append(residue_name)
            residue_numbers.append(residue_number)
            chain_ids.append(chain_id)
    return Structure.from_atom_records(
        xyz=np.array(xyz, dtype=np.float64),
        atom_name=names,
        element=elements,
        residue_name=residue_names,
        residue_number=residue_numbers,
        chain_id=chain_ids,
        source="synthetic",
    )


def alanine(origin: Sequence[float]) -> dict[str, list[float]]:
    """Build an ALA residue whose bonds are exactly the measured lengths.

    Bond *angles* are not realistic and deliberately so: every intra-ALA atom pair is
    within three bonds of every other, so no intra-residue pair is ever contact-checked
    and the angles cannot influence any result. What matters is that each bonded pair sits
    at its :data:`~dodo.validate.reference.RESIDUE_BONDS` length, so the bond graph finds
    it.
    """
    table = RESIDUE_BONDS["ALA"]
    x, y, z = origin
    ca_n = table[("CA", "N")][0]
    ca_c = table[("C", "CA")][0]
    ca_cb = table[("CA", "CB")][0]
    c_o = table[("C", "O")][0]
    return {
        "N": [x + ca_n, y, z],
        "CA": [x, y, z],
        "C": [x - ca_c, y, z],
        "O": [x - ca_c, y + c_o, z],
        "CB": [x, y + ca_cb, z],
    }


def ca_only(origin: Sequence[float]) -> dict[str, list[float]]:
    """Build a rebuilt residue: an alpha carbon and nothing else."""
    return {"CA": list(origin)}


#: Spectator residues inserted between each pair by :func:`spaced`. Three, so that the
#: residues either side are four indices apart -- clear of ``CLASH_EXCLUDE_WITHIN_RESIDUES``,
#: which is 2 and is the wider of the two sequence exclusions.
SPACERS_PER_GAP = 3


def spaced(*residues: ResidueSpec) -> list[ResidueSpec]:
    """Interleave far-away spectator residues between the ones given.

    Residue *indices*, not residue numbers, are what the sequence exclusions work on, and
    two residues built one after another are indices 0 and 1 however they are numbered. A
    test that means "these two residues are far apart in sequence" has to say so by putting
    something between them, or the pair is excluded as a neighbour and the test passes for
    the wrong reason.

    The spacers are spread out rather than stacked, because three residues at the same
    point would themselves be the loudest clash in the structure.
    """
    assert SPACERS_PER_GAP > CLASH_EXCLUDE_WITHIN_RESIDUES
    out: list[ResidueSpec] = []
    spacer = 0
    for i, residue in enumerate(residues):
        if i:
            for _ in range(SPACERS_PER_GAP):
                spacer += 1
                out.append(("UNK", {"C": [1000.0 + 20.0 * spacer, 1000.0, 0.0]}, "A", 900 + spacer))
        out.append(residue)
    return out


def isolated_pair(
    element_a: str,
    element_b: str,
    separation: float,
    *,
    chain_b: str = "A",
) -> Structure:
    """Build two single-atom residues at a controlled separation, with nothing else nearby.

    Each residue is named ``UNK`` so it has no bond table, and carries one atom, so the two
    atoms cannot be bonded, 1-3 or 1-4 related.

    Spectator residues sit between them, 1000 A away; see :func:`spaced` for why.
    """
    return build(
        spaced(
            ("UNK", {element_a: [0.0, 0.0, 0.0]}, "A", 1),
            ("UNK", {element_b: [separation, 0.0, 0.0]}, chain_b, 3),
        )
    )


def sulfur_pair(separation: float) -> Structure:
    """Build two cysteine gamma sulfurs at a controlled separation."""
    return build(
        spaced(
            ("CYS", {"SG": [0.0, 0.0, 0.0]}, "A", 1),
            ("CYS", {"SG": [separation, 0.0, 0.0]}, "A", 3),
        )
    )


# ---------------------------------------------------------------------------
# Limits
# ---------------------------------------------------------------------------


class TestContactLimit:
    def test_carbon_carbon_limit_is_radius_sum_less_allowance(self) -> None:
        assert contact_limit("C", "C") == pytest.approx(2.0 * VDW_RADII["C"] - ALLOWED_OVERLAP)
        assert contact_limit("C", "C") == pytest.approx(2.80)

    def test_polar_pairs_get_the_hydrogen_bond_allowance(self) -> None:
        expected = VDW_RADII["N"] + VDW_RADII["O"] - HBOND_ALLOWED_OVERLAP
        assert contact_limit("N", "O") == pytest.approx(expected)
        assert contact_limit("N", "O") == pytest.approx(2.17)

    def test_polar_limit_is_looser_than_ordinary_limit(self) -> None:
        """Otherwise every hydrogen bond in the input is a finding."""
        assert contact_limit("N", "O") < contact_limit("C", "O")

    def test_sulfur_does_not_get_the_polar_allowance(self) -> None:
        """Measured: sulfur contacts never need it. See POLAR_ELEMENTS."""
        assert contact_limit("S", "O") == pytest.approx(
            VDW_RADII["S"] + VDW_RADII["O"] - ALLOWED_OVERLAP
        )

    def test_case_insensitive(self) -> None:
        assert contact_limit("se", "c") == contact_limit("SE", "C")

    def test_unknown_element_falls_back_to_carbon(self) -> None:
        assert contact_limit("XX", "C") == pytest.approx(contact_limit("C", "C"))


# ---------------------------------------------------------------------------
# CA-only detection
# ---------------------------------------------------------------------------


class TestCaOnlyResidueMask:
    def test_single_ca_residue_is_ca_only(self) -> None:
        structure = build([("ALA", ca_only([0.0, 0.0, 0.0]), "A", 1)])
        assert ca_only_residue_mask(structure).tolist() == [True]

    def test_all_atom_residue_is_not_ca_only(self) -> None:
        structure = build([("ALA", alanine([0.0, 0.0, 0.0]), "A", 1)])
        assert ca_only_residue_mask(structure).tolist() == [False]

    def test_backbone_only_residue_is_not_ca_only(self) -> None:
        """``rebuild(all_atom=True)`` places N, C and O but no side chain.

        That residue has real backbone atoms, so element-aware limits apply to it and are
        strictly more informative than one CA exclusion sphere. Treating it as CA-only
        would silently coarsen the check on exactly the output DODO is moving towards.
        """
        atoms = alanine([0.0, 0.0, 0.0])
        del atoms["CB"]
        structure = build([("ALA", atoms, "A", 1)])
        assert ca_only_residue_mask(structure).tolist() == [False]

    def test_single_atom_residue_that_is_not_a_ca_is_not_ca_only(self) -> None:
        """A lone heteroatom is not a rebuilt residue and must not borrow its limit."""
        structure = build([("UNK", {"ZN": [0.0, 0.0, 0.0]}, "A", 1)])
        assert ca_only_residue_mask(structure).tolist() == [False]

    def test_mixed_structure_reports_per_residue(self) -> None:
        structure = build(
            [
                ("ALA", alanine([0.0, 0.0, 0.0]), "A", 1),
                ("ALA", ca_only([20.0, 0.0, 0.0]), "A", 2),
                ("ALA", alanine([40.0, 0.0, 0.0]), "A", 3),
            ]
        )
        assert ca_only_residue_mask(structure).tolist() == [False, True, False]


# ---------------------------------------------------------------------------
# The bond graph, which is what every exclusion rests on
# ---------------------------------------------------------------------------


class TestCovalentBondPairs:
    def test_intra_residue_bonds_come_from_the_measured_table(self) -> None:
        structure = build([("ALA", alanine([0.0, 0.0, 0.0]), "A", 1)])
        bonds = covalent_bond_pairs(structure)
        found = {
            tuple(sorted((str(structure.atom_name[a]), str(structure.atom_name[b]))))
            for a, b in bonds
        }
        # C-OXT is in the table too, but this residue has no OXT.
        assert found == {pair for pair in RESIDUE_BONDS["ALA"] if "OXT" not in pair}

    def test_ca_only_structure_has_no_covalent_bonds(self) -> None:
        """A CA-CA virtual bond is not a covalent bond.

        A CA-CA virtual bond is not a covalent bond, and pretending otherwise would
        exclude everything within three virtual bonds -- six residues of real trace.
        """
        structure = build([("ALA", ca_only([3.81 * i, 0.0, 0.0]), "A", i + 1) for i in range(5)])
        bonds = covalent_bond_pairs(structure)
        assert bonds.shape == (0, 2)

    def test_peptide_bond_is_found(self) -> None:
        first = alanine([0.0, 0.0, 0.0])
        # Place the next residue's N 1.34 A from this residue's C, along -x.
        second = alanine([-6.0, 0.0, 0.0])
        second["N"] = [first["C"][0] - 1.3394, 0.0, 0.0]
        structure = build([("ALA", first, "A", 1), ("ALA", second, "A", 2)])
        bonds = covalent_bond_pairs(structure)
        pairs = {
            tuple(sorted((int(structure.residue_index[a]), int(structure.residue_index[b]))))
            for a, b in bonds
        }
        assert (0, 1) in pairs

    def test_no_peptide_bond_across_a_chain_break(self) -> None:
        """Residues adjacent in index because the ones between are unmodelled.

        Inventing a bond there would exempt everything within three bonds of it at the far
        end of a gap tens of Angstroms wide, which is a silent hole in the check.
        """
        structure = build(
            [
                ("ALA", alanine([0.0, 0.0, 0.0]), "A", 1),
                ("ALA", alanine([40.0, 0.0, 0.0]), "A", 25),
            ]
        )
        bonds = covalent_bond_pairs(structure)
        crossing = [
            (a, b) for a, b in bonds if structure.residue_index[a] != structure.residue_index[b]
        ]
        assert crossing == []

    def test_no_peptide_bond_across_a_chain_boundary(self) -> None:
        first = alanine([0.0, 0.0, 0.0])
        second = alanine([-6.0, 0.0, 0.0])
        second["N"] = [first["C"][0] - 1.3394, 0.0, 0.0]
        structure = build([("ALA", first, "A", 1), ("ALA", second, "B", 1)])
        bonds = covalent_bond_pairs(structure)
        crossing = [
            (a, b) for a, b in bonds if structure.residue_index[a] != structure.residue_index[b]
        ]
        assert crossing == [], "two chains that happen to touch are not covalently joined"

    def test_disulfide_is_a_bond(self) -> None:
        structure = build(
            [
                (
                    "CYS",
                    {"CA": [0.0, 0.0, 0.0], "CB": [1.53, 0.0, 0.0], "SG": [3.33, 0.0, 0.0]},
                    "A",
                    1,
                ),
                (
                    "CYS",
                    {"CA": [10.0, 0.0, 0.0], "CB": [8.47, 0.0, 0.0], "SG": [5.38, 0.0, 0.0]},
                    "A",
                    20,
                ),
            ]
        )
        bonds = covalent_bond_pairs(structure)
        names = {
            tuple(sorted((str(structure.atom_name[a]), str(structure.atom_name[b]))))
            for a, b in bonds
        }
        assert ("SG", "SG") in names

    def test_terminal_oxt_is_bonded_by_the_distance_fallback(self) -> None:
        """A terminal OXT must be bonded to its own carbonyl carbon.

        OXT is in no residue's bond table. Left unbonded it sits 1.23 A from its own C
        and is the closest contact in the structure -- one false finding per chain, 32
        across the fixtures.
        """
        atoms = alanine([0.0, 0.0, 0.0])
        atoms["OXT"] = [atoms["C"][0], atoms["C"][1] - 1.25, atoms["C"][2]]
        structure = build([("ALA", atoms, "A", 1)])
        bonds = covalent_bond_pairs(structure)
        names = {
            tuple(sorted((str(structure.atom_name[a]), str(structure.atom_name[b]))))
            for a, b in bonds
        }
        assert ("C", "OXT") in names
        assert validate_clashes(structure).ok

    def test_unknown_residue_is_bonded_by_distance(self) -> None:
        """A residue with no bond table still gets its bonds, by distance.

        A modified residue with no table -- here a selenomethionine -- still needs its
        bonds, and the C-SE bond at 1.95 A is longer than the light-atom cutoff allows.
        """
        structure = build(
            [
                (
                    "MSE",
                    {"CG": [0.0, 0.0, 0.0], "SE": [1.95, 0.0, 0.0], "CE": [3.90, 0.0, 0.0]},
                    "A",
                    1,
                )
            ]
        )
        bonds = covalent_bond_pairs(structure)
        assert bonds.shape[0] == 2
        assert validate_clashes(structure).ok

    def test_bonds_are_sorted_deduplicated_and_low_index_first(self) -> None:
        structure = build([("ALA", alanine([0.0, 0.0, 0.0]), "A", 1)])
        bonds = covalent_bond_pairs(structure)
        assert np.all(bonds[:, 0] < bonds[:, 1])
        assert len(np.unique(bonds, axis=0)) == len(bonds)


# ---------------------------------------------------------------------------
# What must not be flagged
# ---------------------------------------------------------------------------


class TestExclusions:
    def test_bonded_atoms_are_not_clashes(self) -> None:
        """Every bond in an ALA is 1.2-1.5 A, far inside every limit."""
        structure = build([("ALA", alanine([0.0, 0.0, 0.0]), "A", 1)])
        assert validate_clashes(structure).ok

    def test_one_three_pair_at_the_measured_minimum_is_not_a_clash(self) -> None:
        """A 1-3 pair at the measured minimum separation is legitimate geometry.

        Measured, 1-3 pairs start at 2.096 A -- inside every limit and fixed by bond
        angles, so nothing can be done about them and nothing should be said.
        """
        atoms = {
            "CA": [0.0, 0.0, 0.0],
            "CB": [1.53, 0.0, 0.0],
            # CG is bonded to CB and 2.10 A from CA: a 1-3 pair at the observed minimum.
            "CG": [1.53 + 1.42, 0.0, 0.0],
        }
        atoms["CG"] = [2.10, 0.0, 0.0]
        structure = build([("LEU", atoms, "A", 1)])
        report = validate_clashes(structure)
        assert report.ok, report.describe()

    def test_one_four_pair_in_a_cis_conformation_is_not_a_clash(self) -> None:
        """A 1-4 pair reaching 2.40 A means its dihedral is near cis, which is legitimate.

        Measured, 1-4 pairs in the fixtures start at 2.404 A, below the 2.80 A C-C limit,
        so without this exclusion the four fixtures report 2,016 of them.
        """
        atoms = {
            "CA": [0.0, 0.0, 0.0],
            "CB": [1.53, 0.0, 0.0],
            "CG": [2.10, 1.42, 0.0],
            # CD is bonded to CG and 2.40 A from CA.
            "CD": [0.9, 2.2, 0.0],
        }
        structure = build([("LYS", atoms, "A", 1)])
        bonds = covalent_bond_pairs(structure)
        assert bonds.shape[0] == 3, "CA-CB, CB-CG, CG-CD"
        distance = float(np.linalg.norm(structure.xyz[0] - structure.xyz[3]))
        assert 2.3 < distance < 2.5, f"the test's own geometry: CA-CD is {distance:.3f} A"
        assert validate_clashes(structure).ok

    def test_disulfide_at_two_angstroms_is_not_a_clash(self) -> None:
        report = validate_clashes(sulfur_pair(2.05))
        assert report.ok, report.describe()

    def test_two_sulfurs_too_far_apart_to_be_a_disulfide_are_a_clash(self) -> None:
        """The exemption is a distance window, not a blanket pass for sulfur."""
        report = validate_clashes(sulfur_pair(DISULFIDE_MAX_DISTANCE + 0.2))
        assert not report.ok
        assert report.violations[0].kind == "atom_contact"

    def test_hydrogen_bond_at_two_point_eight_is_not_a_clash(self) -> None:
        """The single commonest false positive. There are 12,063 of these in the fixtures."""
        report = validate_clashes(isolated_pair("N", "O", 2.80))
        assert report.ok, report.describe()

    def test_short_hydrogen_bond_at_the_measured_minimum_is_not_a_clash(self) -> None:
        """The tightest hydrogen bond in the fixture set is not a clash.

        2.293 A, the ASN611 ND2 to HIS613 NE2 bond in dnmt3a. Real, and the tightest
        polar contact in the fixture set.
        """
        assert validate_clashes(isolated_pair("N", "N", 2.293)).ok

    def test_polar_pair_below_the_polar_limit_is_still_a_clash(self) -> None:
        """The allowance is a relaxation, not an exemption."""
        report = validate_clashes(isolated_pair("N", "O", 2.00))
        assert not report.ok

    def test_consecutive_residues_are_not_checked_against_each_other(self) -> None:
        """Sequence neighbours are exempt from the all-atom rule.

        Neighbours are held close by the backbone whether or not the specific pair is
        within three bonds. Measured, the closest such contact in the fixtures is 2.484 A.
        """
        structure = build(
            [
                ("ALA", {"CB": [0.0, 0.0, 0.0]}, "A", 1),
                ("ALA", {"CB": [2.0, 0.0, 0.0]}, "A", 2),
            ]
        )
        assert validate_clashes(structure).ok

    def test_residues_two_apart_are_checked(self) -> None:
        """Two apart is where genuine packing errors start, so the exclusion stops at one."""
        structure = build(
            [
                ("ALA", {"CB": [0.0, 0.0, 0.0]}, "A", 1),
                ("ALA", {"CB": [50.0, 0.0, 0.0]}, "A", 2),
                ("ALA", {"CB": [2.0, 0.0, 0.0]}, "A", 3),
            ]
        )
        report = validate_clashes(structure)
        assert not report.ok
        assert {report.violations[0].residue_a, report.violations[0].residue_b} == {0, 2}

    def test_sequence_exclusion_does_not_cross_chains(self) -> None:
        """The sequence exclusion must not leak across a chain boundary.

        Index adjacency between two different chains means nothing; a collision there is
        a collision.
        """
        structure = build(
            [
                ("ALA", {"CB": [0.0, 0.0, 0.0]}, "A", 1),
                ("ALA", {"CB": [2.0, 0.0, 0.0]}, "B", 1),
            ]
        )
        assert not validate_clashes(structure).ok

    def test_hydrogens_are_skipped(self) -> None:
        """Hydrogens take no part in the check.

        An explicit-H input names hydrogens H, HB2 and so on, none of which are in the
        measured bond tables, so each would be an unbonded atom 1 A from its parent.
        """
        atoms = alanine([0.0, 0.0, 0.0])
        atoms["H"] = [atoms["N"][0] + 1.0, 0.0, 0.0]
        structure = build([("ALA", atoms, "A", 1)])
        report = validate_clashes(structure)
        assert report.ok
        assert report.n_atoms_skipped == 1


# ---------------------------------------------------------------------------
# What must be flagged
# ---------------------------------------------------------------------------


class TestDetection:
    def test_two_carbons_on_top_of_each_other(self) -> None:
        report = validate_clashes(isolated_pair("C", "C", 0.4))
        assert not report.ok
        assert report.violations[0].separation == pytest.approx(0.4)
        assert report.violations[0].limit == pytest.approx(2.80)

    def test_coincident_atoms(self) -> None:
        """0.000 A. The defect the pre-rewrite output fixture is full of."""
        report = validate_clashes(isolated_pair("C", "O", 0.0))
        assert len(report.violations) == 1
        assert report.violations[0].separation == pytest.approx(0.0)

    def test_contact_just_inside_the_limit_is_flagged(self) -> None:
        assert not validate_clashes(isolated_pair("C", "C", 2.79)).ok

    def test_contact_just_outside_the_limit_is_not(self) -> None:
        assert validate_clashes(isolated_pair("C", "C", 2.81)).ok

    def test_element_awareness_changes_the_verdict(self) -> None:
        """2.5 A is a clash between two carbons and a hydrogen bond between an N and an O.

        A single-number cutoff cannot express that, which is the entire argument for
        element-aware limits.
        """
        assert not validate_clashes(isolated_pair("C", "C", 2.50)).ok
        assert validate_clashes(isolated_pair("N", "O", 2.50)).ok

    def test_rebuilt_ca_driven_through_a_folded_side_chain(self) -> None:
        """The CA rule reaches into all-atom geometry, not only CA-to-CA.

        The defect measured in real ``rebuild()`` output: a rebuilt CA 0.121 A from a
        folded domain's glutamate. The CA rule has to reach into all-atom geometry, not
        only CA-to-CA, or this is invisible.
        """
        structure = build(
            spaced(
                ("ALA", ca_only([0.0, 0.0, 0.0]), "A", 1),
                ("GLU", alanine([30.0, 0.0, 0.0]), "A", 12),
            )
        )
        # Put the folded residue's CB right on top of the rebuilt CA.
        structure.xyz[structure.atom_name == "CB"] = [0.12, 0.0, 0.0]
        report = validate_clashes(structure)
        assert len(report.violations) == 1
        violation = report.violations[0]
        assert violation.kind == "ca_contact"
        assert violation.limit == pytest.approx(CA_CLASH_DISTANCE)
        assert "alpha-carbon only" in violation.message


class TestCaOnlyRule:
    def test_two_rebuilt_cas_closer_than_the_clash_distance(self) -> None:
        structure = build(
            [
                ("ALA", ca_only([0.0, 0.0, 0.0]), "A", 1),
                ("ALA", ca_only([100.0, 0.0, 0.0]), "A", 2),
                ("ALA", ca_only([200.0, 0.0, 0.0]), "A", 3),
                ("ALA", ca_only([3.0, 0.0, 0.0]), "A", 4),
            ]
        )
        report = validate_clashes(structure)
        assert len(report.violations) == 1
        assert report.violations[0].kind == "ca_contact"
        assert report.violations[0].limit == pytest.approx(CA_CLASH_DISTANCE)

    def test_cas_at_the_clash_distance_are_fine(self) -> None:
        structure = build(
            [
                ("ALA", ca_only([0.0, 0.0, 0.0]), "A", 1),
                ("ALA", ca_only([100.0, 0.0, 0.0]), "A", 2),
                ("ALA", ca_only([200.0, 0.0, 0.0]), "A", 3),
                ("ALA", ca_only([CA_CLASH_DISTANCE + 0.01, 0.0, 0.0]), "A", 4),
            ]
        )
        assert validate_clashes(structure).ok

    def test_sequence_neighbours_within_the_exclusion_are_not_checked(self) -> None:
        """A CA trace's own near neighbours are exempt.

        A CA trace's own bonded neighbours sit at 3.81 A and its next-nearest at
        5.0-7.5 A; a tight turn brings i and i+2 closer still. Checking them reports every
        peptide bond, which is what the pre-rewrite whole-structure check did.
        """
        structure = build(
            [
                ("ALA", ca_only([0.0, 0.0, 0.0]), "A", 1),
                ("ALA", ca_only([1.0, 0.0, 0.0]), "A", 2),
                ("ALA", ca_only([2.0, 0.0, 0.0]), "A", 3),
            ]
        )
        assert CLASH_EXCLUDE_WITHIN_RESIDUES == 2
        assert validate_clashes(structure).ok

    def test_rebuilt_to_folded_boundary_is_not_a_false_positive(self) -> None:
        """THE regression test for this module.

        At the join between a rebuilt region and a folded domain, the rebuilt CA(i) sits
        about 2.43 A from N(i+1) of the folded residue it connects to -- correct backbone
        geometry, and closer than the 3.20 A CA rule. Residue i is CA-only so there is no
        C(i) and therefore no bond path to exclude the pair by. Only the sequence
        exclusion saves it, and without that every rebuilt region in every DODO structure
        reports a clash at each of its two ends.
        """
        folded = alanine([2.43 + 1.4561, 0.0, 0.0])
        folded["N"] = [2.43, 0.0, 0.0]
        structure = build([("ALA", ca_only([0.0, 0.0, 0.0]), "A", 41), ("ALA", folded, "A", 42)])
        distance = float(
            np.linalg.norm(structure.xyz[0] - structure.xyz[structure.atom_name == "N"][0])
        )
        assert distance == pytest.approx(2.43), "the test's own geometry"
        report = validate_clashes(structure)
        assert report.ok, report.describe()

    def test_ca_rule_applies_even_when_the_other_residue_is_all_atom(self) -> None:
        structure = build(
            spaced(
                ("ALA", ca_only([0.0, 0.0, 0.0]), "A", 1),
                # Two atoms, so this residue is all-atom rather than CA-only.
                ("ALA", {"CA": [2.9, 0.0, 0.0], "CB": [2.9, 1.5337, 0.0]}, "A", 20),
            )
        )
        report = validate_clashes(structure)
        assert not report.ok
        assert all(v.kind == "ca_contact" for v in report.violations)
        assert all(v.limit == pytest.approx(CA_CLASH_DISTANCE) for v in report.violations)

    def test_all_atom_pair_at_the_same_distance_is_not_flagged(self) -> None:
        """The all-atom rule is more permissive than the CA rule, deliberately.

        The counterpart of the previous test: 2.9 A between two carbons of two all-atom
        residues is tight packing, not a clash. The CA rule is stricter on purpose, because
        one atom is standing in for a whole residue.
        """
        structure = build(
            spaced(
                ("ALA", {"CA": [0.0, 0.0, 0.0], "CB": [0.0, 1.5337, 0.0]}, "A", 1),
                ("ALA", {"CA": [2.9, 0.0, 0.0], "CB": [2.9, 1.5337, 0.0]}, "A", 20),
            )
        )
        assert validate_clashes(structure).ok


# ---------------------------------------------------------------------------
# Report surface
# ---------------------------------------------------------------------------


class TestClashReport:
    def test_clean_report_is_truthy_and_ok(self) -> None:
        report = validate_clashes(isolated_pair("C", "C", 5.0))
        assert report.ok
        assert bool(report) is True
        assert report.violations == ()
        assert report.worst is None

    def test_dirty_report_is_falsy(self) -> None:
        report = validate_clashes(isolated_pair("C", "C", 1.0))
        assert not report.ok
        assert bool(report) is False

    def test_message_names_both_residues_the_value_and_the_limit(self) -> None:
        structure = build(
            spaced(
                ("GLU", {"OE1": [0.0, 0.0, 0.0]}, "A", 142),
                ("LYS", {"NZ": [1.5, 0.0, 0.0]}, "B", 300),
            )
        )
        message = validate_clashes(structure).violations[0].message
        assert "A:GLU142" in message
        assert "B:LYS300" in message
        assert "1.500" in message
        assert "2.17" in message
        assert message.endswith(".")

    def test_messages_are_the_structures_own_residue_labels(self) -> None:
        structure = build(
            spaced(
                ("ALA", ca_only([0.0, 0.0, 0.0]), "A", 7),
                ("ALA", ca_only([1.0, 0.0, 0.0]), "A", 99),
            )
        )
        violation = validate_clashes(structure).violations[0]
        assert violation.label_a == structure.residue_label(violation.residue_a)
        assert violation.label_b == structure.residue_label(violation.residue_b or 0)

    def test_worst_ranks_by_overlap_not_by_raw_distance(self) -> None:
        """``worst`` ranks by overlap depth rather than by raw separation.

        A 2.30 A N-O contact is 0.13 A inside its limit; a 2.30 A CA-CA contact is 0.90 A
        inside its. Ranking by distance alone would call them equally bad.
        """
        structure = build(
            spaced(
                ("ALA", ca_only([0.0, 0.0, 0.0]), "A", 1),
                ("ALA", ca_only([2.30, 0.0, 0.0]), "A", 10),
                ("UNK", {"N": [100.0, 0.0, 0.0]}, "A", 20),
                ("UNK", {"O": [102.15, 0.0, 0.0]}, "A", 30),
            )
        )
        report = validate_clashes(structure)
        assert len(report.violations) == 2
        worst = report.worst
        assert worst is not None
        assert worst.kind == "ca_contact"
        assert worst.overlap > 0.8

    def test_of_kind_partitions_the_violations(self) -> None:
        structure = build(
            spaced(
                ("ALA", ca_only([0.0, 0.0, 0.0]), "A", 1),
                ("ALA", ca_only([2.0, 0.0, 0.0]), "A", 10),
                ("UNK", {"C": [100.0, 0.0, 0.0]}, "A", 20),
                ("UNK", {"C": [101.0, 0.0, 0.0]}, "A", 30),
            )
        )
        report = validate_clashes(structure)
        assert len(report.of_kind("ca_contact")) == 1
        assert len(report.of_kind("atom_contact")) == 1
        assert len(report.of_kind("ca_contact")) + len(report.of_kind("atom_contact")) == len(
            report.violations
        )

    def test_violating_residues_covers_both_ends(self) -> None:
        structure = build(
            spaced(
                ("UNK", {"C": [0.0, 0.0, 0.0]}, "A", 1),
                ("UNK", {"C": [1.0, 0.0, 0.0]}, "A", 20),
            )
        )
        assert validate_clashes(structure).violating_residues == (0, 1 + SPACERS_PER_GAP)

    def test_describe_truncates_and_says_how_many_remain(self) -> None:
        structure = build(
            spaced(*[("UNK", {"C": [0.4 * i, 0.0, 0.0]}, "A", 10 * (i + 1)) for i in range(8)])
        )
        report = validate_clashes(structure)
        assert len(report.violations) > 3
        described = report.describe(max_violations=3)
        assert described.count("  - ") == 4  # three violations plus the "and N more" line
        assert "more clash(es)" in described

    def test_describe_lists_the_worst_first(self) -> None:
        structure = build(
            spaced(
                ("UNK", {"C": [0.0, 0.0, 0.0]}, "A", 1),
                ("UNK", {"C": [2.70, 0.0, 0.0]}, "A", 20),
                ("UNK", {"C": [0.0, 2.75, 0.0]}, "A", 40),
            )
        )
        report = validate_clashes(structure)
        assert len(report.violations) == 2
        first_line = report.describe().splitlines()[1]
        assert "2.700" in first_line, "the 2.70 A pair is the deeper overlap"

    def test_summary_reports_what_was_checked(self) -> None:
        structure = isolated_pair("C", "C", 5.0)
        summary = validate_clashes(structure).summary()
        assert f"{structure.n_atoms} atoms" in summary
        assert "no clashes" in summary

    def test_raise_if_invalid_is_silent_when_clean(self) -> None:
        validate_clashes(isolated_pair("C", "C", 5.0)).raise_if_invalid()

    def test_raise_if_invalid_raises_geometry_error(self) -> None:
        report = validate_clashes(isolated_pair("C", "C", 0.5))
        with pytest.raises(GeometryError, match="Steric clashes detected"):
            report.raise_if_invalid()

    def test_report_is_hashable_and_comparable(self) -> None:
        """Two reports compare and hash without raising.

        ``eq=False`` exists so that comparing two reports does not raise "truth value of
        an array is ambiguous".
        """
        a = validate_clashes(isolated_pair("C", "C", 5.0))
        b = validate_clashes(isolated_pair("C", "C", 5.0))
        assert a != b
        assert len({a, b}) == 2

    def test_pair_accounting_shows_the_exclusions_are_load_bearing(self) -> None:
        structure = build([("ALA", alanine([0.0, 0.0, 0.0]), "A", 1)])
        report = validate_clashes(structure)
        assert report.n_pairs_excluded == 8, (
            "every intra-ALA pair inside the search radius is within three bonds; the two "
            "that are not counted here are further apart than the largest limit"
        )
        assert report.n_pairs_checked == 0


# ---------------------------------------------------------------------------
# Unusable input
# ---------------------------------------------------------------------------


class TestUnusableInput:
    def test_empty_structure_raises(self) -> None:
        structure = build([("ALA", ca_only([0.0, 0.0, 0.0]), "A", 1)])
        structure.xyz = np.zeros((0, 3))
        with pytest.raises(EmptyStructureError):
            validate_clashes(structure)

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"allowed_overlap": -0.1},
            {"hbond_allowed_overlap": float("nan")},
            {"ca_clash_distance": float("inf")},
            {"ca_exclude_within": -1},
            {"atom_exclude_within": -1},
        ],
    )
    def test_bad_arguments_raise(self, kwargs: dict[str, float]) -> None:
        with pytest.raises(GeometryError):
            validate_clashes(isolated_pair("C", "C", 5.0), **kwargs)

    def test_single_atom_structure_is_a_clean_report_not_an_error(self) -> None:
        report = validate_clashes(build([("ALA", ca_only([0.0, 0.0, 0.0]), "A", 1)]))
        assert report.ok
        assert report.n_pairs_checked == 0

    def test_non_finite_coordinate_is_reported_not_raised(self) -> None:
        """A non-finite coordinate is reported rather than raised.

        A NaN coordinate makes cKDTree raise something opaque. Reporting it against its
        residue is the whole reason this returns a report.
        """
        structure = isolated_pair("C", "C", 5.0)
        structure.xyz[1] = [np.nan, 0.0, 0.0]
        report = validate_clashes(structure)
        assert len(report.of_kind("non_finite")) == 1
        violation = report.of_kind("non_finite")[0]
        assert violation.atom_b is None
        assert violation.label_a in violation.message
        assert not report.ok

    def test_finite_atoms_are_still_checked_around_a_non_finite_one(self) -> None:
        structure = build(
            [
                ("UNK", {"C": [0.0, 0.0, 0.0]}, "A", 1),
                ("UNK", {"C": [np.nan, 0.0, 0.0]}, "A", 20),
                ("UNK", {"C": [1.0, 0.0, 0.0]}, "A", 40),
            ]
        )
        report = validate_clashes(structure)
        assert len(report.of_kind("non_finite")) == 1
        assert len(report.of_kind("atom_contact")) == 1

    def test_worst_falls_back_to_a_non_finite_violation(self) -> None:
        structure = isolated_pair("C", "C", 5.0)
        structure.xyz[1] = [np.nan, 0.0, 0.0]
        worst = validate_clashes(structure).worst
        assert worst is not None
        assert worst.kind == "non_finite"


# ---------------------------------------------------------------------------
# Calibration against the real fixtures
# ---------------------------------------------------------------------------

#: Read each fixture once; 6kn7 is 5 MB of PDB.
_CACHE: dict[str, Structure] = {}


def fixture_structure(name: str) -> Structure:
    if name not in _CACHE:
        _CACHE[name] = read_structure(DATA / name)
    return _CACHE[name]


class TestFixtureCalibration:
    """The six unmodified fixtures.

    Four are real deposited or predicted structures and should be essentially clean. Two
    -- ``test.pdb`` and ``testing_translation.pdb`` -- are pre-rewrite DODO output and are
    the positive controls: ``test.pdb``'s CA-CA bonds run 1.23 to 562 A, and
    ``testing_translation.pdb`` still has the orphaned non-CA atoms fixed in commit
    026c4f1, some of them coincident to 0.000 A. If the validator ever reports those two
    as clean it has stopped working.
    """

    def test_alphafold_dnmt3a_is_clean(self) -> None:
        report = validate_clashes(fixture_structure("dnmt3a.pdb"))
        assert report.violations == (), report.describe()

    def test_alphafold_arf19_is_clean(self) -> None:
        report = validate_clashes(fixture_structure("arf19.pdb"))
        assert report.violations == (), report.describe()

    def test_alphafold_p300_reports_only_its_own_broken_patch(self) -> None:
        """Measured: 15 violations, 13 of which are MET121 superimposed on PRO659.

        That is genuinely in the input -- both residues have pLDDT near 30, and AlphaFold
        piled two low-confidence stretches 538 residues apart on top of each other, CA to
        CA at 2.336 A. It is not DODO's doing and it is not this module's error; it is a
        defect in the file, and reporting it is correct.
        """
        report = validate_clashes(fixture_structure("p300.pdb"))
        assert len(report.violations) <= 20, report.describe(20)
        labels = {(v.label_a, v.label_b) for v in report.violations}
        assert ("A:MET121", "A:PRO659") in labels
        broken_patch = [v for v in report.violations if v.label_a == "A:MET121"]
        assert len(broken_patch) >= 10

    def test_em_assembly_findings_are_sparse_and_real(self) -> None:
        """Findings on the EM assembly are sparse, and real.

        6kn7 is a 6.60 A electron microscopy model: its side chains are not resolved by the
        data at all, and the deposited coordinates contain real overlaps.

        Measured, 131 violations over 61,511 atoms -- one per 470 atoms -- at 77 distinct
        sequence positions, each repeated across up to 9 of the 29 NCS copies. That
        repetition is the evidence they are features of the model rather than artefacts of
        this check: a validator bug would not land on the same residue pair in nine
        independent chains.
        """
        structure = fixture_structure("6kn7.pdb")
        report = validate_clashes(structure)
        assert len(report.violations) <= 160, report.describe(20)
        assert len(report.violations) < structure.n_atoms / 300
        positions = {
            (v.atom_name_a, v.atom_name_b, int(structure.residue_number[v.residue_a]))
            for v in report.violations
        }
        assert len(positions) < len(report.violations), "findings recur across the NCS copies"

    @pytest.mark.parametrize(
        ("pdb", "cif"),
        [("arf19.pdb", "arf19.cif"), ("p300.pdb", "p300.cif"), ("6kn7.pdb", "6kn7.cif")],
    )
    def test_the_same_structure_read_two_ways_gives_the_same_answer(
        self, pdb: str, cif: str
    ) -> None:
        """The findings must be properties of the coordinates, not of the reader.

        Three fixtures ship in both formats. Measured, the counts agree exactly (0, 15 and
        131), which is worth asserting: an exclusion that depended on atom ordering, on
        chain-id spelling or on how altLocs were resolved would diverge here and nowhere
        else in this file.
        """
        from_pdb = validate_clashes(fixture_structure(pdb))
        from_cif = validate_clashes(fixture_structure(cif))
        assert len(from_pdb.violations) == len(from_cif.violations)
        assert {round(v.separation, 3) for v in from_pdb.violations} == {
            round(v.separation, 3) for v in from_cif.violations
        }

    def test_pre_rewrite_ca_trace_fixture_is_caught(self) -> None:
        report = validate_clashes(fixture_structure("test.pdb"))
        assert len(report.violations) > 1000
        assert report.n_ca_only_residues == report.n_residues

    def test_pre_rewrite_output_with_orphaned_atoms_is_caught(self) -> None:
        """The v1 output fixture with orphaned atoms is caught.

        Measured: 1,301 violations, including atoms at 0.000 A. This fixture is v1
        output from before the orphaned-atom fix, so it contains both regimes at once and
        exercises the CA rule and the all-atom rule in the same structure.
        """
        report = validate_clashes(fixture_structure("testing_translation.pdb"))
        assert len(report.violations) > 500
        assert report.of_kind("ca_contact")
        assert report.of_kind("atom_contact")
        worst = report.worst
        assert worst is not None
        assert worst.separation < 1.0

    @pytest.mark.parametrize(
        "name", ["dnmt3a.pdb", "arf19.pdb", "p300.pdb", "6kn7.pdb", "testing_translation.pdb"]
    )
    def test_every_finding_is_actionable(self, name: str) -> None:
        """Every finding names the residues, the value and the limit.

        A finding a user cannot act on is noise: it must name both residues, both atoms,
        the measured distance and the limit, in one sentence.
        """
        structure = fixture_structure(name)
        for violation in validate_clashes(structure).violations:
            assert violation.label_a == structure.residue_label(violation.residue_a)
            assert violation.residue_b is not None
            assert violation.label_b == structure.residue_label(violation.residue_b)
            assert violation.label_a in violation.message
            assert violation.label_b in violation.message
            assert f"{violation.separation:.3f}" in violation.message
            assert violation.message.endswith(".")
            assert violation.separation < violation.limit

    def test_no_finding_is_between_bonded_or_neighbouring_residues(self) -> None:
        """The exclusions have to hold on real structures, not only on synthetic ones."""
        structure = fixture_structure("6kn7.pdb")
        bonds = {tuple(pair) for pair in covalent_bond_pairs(structure)}
        for violation in validate_clashes(structure).violations:
            assert violation.atom_b is not None
            assert (violation.atom_a, violation.atom_b) not in bonds
            assert violation.residue_b is not None
            same_chain = (
                structure.chain_index[violation.residue_a]
                == structure.chain_index[violation.residue_b]
            )
            if same_chain:
                assert abs(violation.residue_a - violation.residue_b) > 1


class TestMeasuredClaims:
    """Locks in the numbers the module's own docstrings quote as justification.

    If someone widens or narrows an allowance, these fail and point at the sentence that
    is now wrong.
    """

    def test_disabling_the_hydrogen_bond_allowance_costs_findings(self) -> None:
        structure = fixture_structure("dnmt3a.pdb")
        assert validate_clashes(structure).ok
        without = validate_clashes(structure, hbond_allowed_overlap=ALLOWED_OVERLAP)
        assert len(without.violations) == 1, (
            "measured: treating polar pairs like any other adds 1 finding on dnmt3a, "
            "which is a real hydrogen bond"
        )

    def test_molprobity_strictness_reports_ordinary_packing(self) -> None:
        """A crystallographer's strictness reports ordinary packing here.

        0.40 A of allowed overlap is MolProbity's clash criterion. It is right for
        re-refined crystal structures and wrong here: measured, it reports 21 contacts in
        dnmt3a and 16 in arf19, both of which are clean AlphaFold models.
        """
        strict = validate_clashes(
            fixture_structure("dnmt3a.pdb"), allowed_overlap=0.40, hbond_allowed_overlap=0.40
        )
        assert 10 <= len(strict.violations) <= 40

    def test_fixtures_contain_the_disulfides_the_exemption_is_for(self) -> None:
        for name in ("dnmt3a.pdb", "arf19.pdb", "p300.pdb"):
            structure = fixture_structure(name)
            bonds = covalent_bond_pairs(structure)
            sulfur_bridges = [
                (a, b)
                for a, b in bonds
                if structure.atom_name[a] == "SG" and structure.atom_name[b] == "SG"
            ]
            assert len(sulfur_bridges) == 1, f"{name} has one disulfide"
            distance = float(
                np.linalg.norm(
                    structure.xyz[sulfur_bridges[0][0]] - structure.xyz[sulfur_bridges[0][1]]
                )
            )
            assert 2.0 <= distance <= 2.1


class TestPerformance:
    def test_scales_to_a_sixty_thousand_atom_assembly(self) -> None:
        """The check scales to a 61,511-atom assembly.

        6kn7 has 61,511 atoms. Measured on this machine: 0.178 s median over 5 runs,
        of which 0.067 s is building the covalent bond graph. The generous bound here is
        for slower CI hardware; anything that regresses to nested loops takes minutes.
        """
        structure = fixture_structure("6kn7.pdb")
        start = time.perf_counter()
        report = validate_clashes(structure)
        elapsed = time.perf_counter() - start
        assert structure.n_atoms == 61511
        assert elapsed < 10.0, f"took {elapsed:.2f} s"
        assert report.n_pairs_checked + report.n_pairs_excluded > 100_000


# ---------------------------------------------------------------------------
# Against freshly rebuilt output
# ---------------------------------------------------------------------------


class TestRebuiltOutput:
    """Validation of what ``dodo.rebuild()`` actually produces today.

    Not asserted clean, because it is not: measured on dnmt3a with ``seed=7`` the output
    contains 13 clashes, all of them a rebuilt CA driven into the side chain of the folded
    residue that anchors the region (A:ALA274 CA is 0.121 A from A:GLU283 CD, and GLU283
    is the IDR's own C-terminal anchor). That is a real DODO defect and this module's job
    is to surface it, so these tests assert the *shape* of the report rather than a clean
    bill of health.
    """

    @staticmethod
    def rebuilt(name: str, **kwargs: object) -> Structure:
        import dodo

        report = dodo.rebuild(DATA / name, seed=7, **kwargs)
        model: Structure = report.models[0]
        return model

    @pytest.mark.slow
    def test_rebuilt_output_is_a_mix_of_both_regimes(self) -> None:
        model = self.rebuilt("dnmt3a.pdb")
        report = validate_clashes(model)
        assert 0 < report.n_ca_only_residues < report.n_residues, (
            "rebuild() must produce CA-only regions and all-atom folded domains in one "
            "structure; if it does not, this module's two-regime handling is untested"
        )
        assert report.n_pairs_excluded > report.n_pairs_checked

    @pytest.mark.slow
    def test_findings_on_rebuilt_output_are_few_and_local_to_region_boundaries(self) -> None:
        model = self.rebuilt("dnmt3a.pdb")
        report = validate_clashes(model)
        # Measured: 13. A regression that makes every rebuilt residue clash would blow
        # past this bound immediately, and so would a broken exclusion.
        assert len(report.violations) < 0.05 * model.n_residues, report.describe(10)
        for violation in report.violations:
            assert violation.label_a in violation.message
            assert violation.separation < violation.limit

    @pytest.mark.slow
    def test_no_clash_is_reported_inside_an_untouched_folded_domain(self) -> None:
        """Nothing inside an untouched folded domain is ever reported.

        Folded domains move as rigid bodies, so their internal geometry is exactly the
        input's. dnmt3a's input is clean, therefore every finding on its rebuilt output
        must involve at least one rebuilt residue. A finding between two folded residues
        would mean the rigid-body move was not rigid.
        """
        model = self.rebuilt("dnmt3a.pdb")
        ca_only_residues = ca_only_residue_mask(model)
        for violation in validate_clashes(model).violations:
            assert violation.residue_b is not None
            touched = ca_only_residues[violation.residue_a] or ca_only_residues[violation.residue_b]
            assert touched, f"two folded residues clash after a rigid-body move: {violation}"

    @pytest.mark.slow
    @pytest.mark.parametrize("mode", ["compact", "predicted", "expanded"])
    def test_every_mode_produces_a_usable_report(self, mode: str) -> None:
        model = self.rebuilt("dnmt3a.pdb", mode=mode)
        report = validate_clashes(model)
        assert report.n_pairs_checked > 0
        assert not report.of_kind("non_finite"), "a rebuild must never emit NaN coordinates"
