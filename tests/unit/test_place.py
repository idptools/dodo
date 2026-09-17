"""Tests for step 3: moving folded domains so linkers can reach their predicted dimensions.

There was no test module for this before rigid units were added, which is how the two covalent
defects in :mod:`test_assembly` survived. The assertions here are mostly *differential* -- what
the structure looked like going in against what it looks like coming out -- because that is the
only form in which the invariant can be stated: a unit's internal geometry is unchanged, and
everything else is allowed to move.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.spatial import cKDTree

from dodo.constants import CA_CLASH_DISTANCE
from dodo.construct.assembly import find_rigid_units
from dodo.construct.place import (
    PlacementReport,
    reposition_folded_domains,
    verify_rigid,
)
from dodo.exceptions import BuildError, GeometryError
from dodo.io import read_structure
from dodo.regions.identify import assign_regions
from dodo.structure import DomainKind

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "structures"
COMPLEX_FIXTURES = Path(__file__).resolve().parents[1] / "data" / "complexes_full_sequences_modeled"


def _assigned(name: str):
    structure = read_structure(FIXTURES / name)
    assign_regions(structure)
    return structure


def _unit_contacts(structure, assembly, cutoff: float = 5.0):
    """Distances between atoms of one rigid unit that are in contact. The thing that must hold."""
    out = []
    for unit in assembly.units:
        if len(unit) < 2:
            continue
        coords = unit.xyz
        pairs = cKDTree(coords).query_pairs(cutoff, output_type="ndarray")
        if pairs.size:
            out.append(np.linalg.norm(coords[pairs[:, 0]] - coords[pairs[:, 1]], axis=1))
    return np.concatenate(out) if out else np.zeros(0)


class TestSingleChain:
    def test_the_first_folded_domain_holds_the_frame(self) -> None:
        structure = _assigned("dnmt3a.pdb")
        first = min(
            (d for d in structure.domains if d.kind is DomainKind.FOLDED),
            key=lambda d: d.span.start,
        )
        before = first.xyz.copy()
        reposition_folded_domains(structure, rng=np.random.default_rng(0))
        assert np.array_equal(first.xyz, before)

    def test_a_linker_lands_on_its_predicted_separation(self) -> None:
        structure = _assigned("dnmt3a.pdb")
        report = reposition_folded_domains(structure, rng=np.random.default_rng(0))
        moved = [p for p in report.placements if p.moved]
        assert moved
        for placement in moved:
            assert placement.target_separation is not None
            assert placement.achieved_separation == pytest.approx(
                placement.target_separation, abs=0.5
            ), str(placement)

    def test_every_domain_moves_rigidly(self) -> None:
        structure = _assigned("p300.pdb")
        before = {d.span.start: d.xyz.copy() for d in structure.domains}
        reposition_folded_domains(structure, rng=np.random.default_rng(1))
        for domain in structure.domains:
            verify_rigid(before[domain.span.start], domain.xyz)

    def test_the_same_seed_gives_the_same_answer(self) -> None:
        first, second = _assigned("dnmt3a.pdb"), _assigned("dnmt3a.pdb")
        reposition_folded_domains(first, rng=np.random.default_rng(7))
        reposition_folded_domains(second, rng=np.random.default_rng(7))
        assert np.array_equal(first.xyz, second.xyz)

    def test_unassigned_regions_are_refused(self) -> None:
        structure = read_structure(FIXTURES / "dnmt3a.pdb")
        with pytest.raises(BuildError, match="no assigned regions"):
            reposition_folded_domains(structure, rng=np.random.default_rng(0))


class TestComplexes:
    """The interface must survive, exactly, and the report must say so."""

    def test_a_locked_interface_is_preserved_to_machine_precision(self) -> None:
        structure = _assigned("dnmt3a_dimer.pdb")
        assembly = find_rigid_units(structure)
        before = _unit_contacts(structure, assembly)
        assert before.size > 1000, "the fixture has no interface to preserve"
        reposition_folded_domains(
            structure, rng=np.random.default_rng(0), assembly=assembly
        )
        after = _unit_contacts(structure, assembly)
        assert np.abs(after - before).max() < 1e-9

    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_the_interface_survives_at_every_seed(self, seed: int) -> None:
        structure = _assigned("dnmt3a_dimer.pdb")
        assembly = find_rigid_units(structure)
        before = _unit_contacts(structure, assembly)
        reposition_folded_domains(
            structure, rng=np.random.default_rng(seed), assembly=assembly
        )
        assert np.abs(_unit_contacts(structure, assembly) - before).max() < 1e-9

    def test_without_locking_the_interface_is_destroyed(self) -> None:
        """The control, and it has to measure the SAME atom pairs.

        Without it the test above could pass for the wrong reason -- an assertion that
        distances are preserved proves nothing unless something is capable of changing them.
        Measured here: the interface goes from 5,540 atom contacts within 5 A to none.
        """
        structure = _assigned("dnmt3a_dimer.pdb")
        chain_of = np.repeat(
            structure.chain_index, np.diff(structure.residue_atom_offsets)
        )
        pairs = np.array(
            sorted(cKDTree(structure.xyz).query_pairs(5.0)), dtype=np.int64
        )
        across = pairs[chain_of[pairs[:, 0]] != chain_of[pairs[:, 1]]]
        assert across.shape[0] > 1000, "the fixture has no inter-chain interface"
        reposition_folded_domains(
            structure,
            rng=np.random.default_rng(0),
            assembly=find_rigid_units(structure, lock_interfaces=False),
        )
        after = np.linalg.norm(
            structure.xyz[across[:, 0]] - structure.xyz[across[:, 1]], axis=1
        )
        surviving = int((after <= 5.0).sum())
        assert surviving < 0.5 * across.shape[0], (
            f"{surviving} of {across.shape[0]} inter-chain contacts survived without locking, "
            f"so the locked test above is not proving anything"
        )

    def test_a_unit_reached_by_two_linkers_satisfies_both(self) -> None:
        """The multi-constraint path. Each chain's middle domain is held on both sides."""
        structure = _assigned("dnmt3a_dimer.pdb")
        report = reposition_folded_domains(structure, rng=np.random.default_rng(0))
        multi = [u for u in report.units if u.moved and u.n_constraints > 1]
        assert multi, "the fixture no longer exercises the multi-constraint placer"
        for unit in multi:
            assert unit.worst_residual is not None
            assert unit.worst_residual < 2.0, str(unit)

    def test_a_fully_locked_assembly_moves_nothing(self) -> None:
        structure = _assigned("6kn7.pdb")
        before = structure.xyz.copy()
        report = reposition_folded_domains(structure, rng=np.random.default_rng(0))
        assert np.array_equal(structure.xyz, before)
        assert not report.moved

    def test_other_chains_are_obstacles(self) -> None:
        """Every moved unit must clear every already-final atom of every chain.

        Step 3 used to reset its obstacle set per chain, so two chains could be placed straight
        through each other and nothing would notice.
        """
        structure = _assigned("dnmt3a_dimer.pdb")
        assembly = find_rigid_units(structure)
        reposition_folded_domains(structure, rng=np.random.default_rng(0), assembly=assembly)
        folded = np.zeros(structure.n_atoms, dtype=bool)
        owner = np.full(structure.n_atoms, -1)
        for unit in assembly.units:
            mask = unit.atom_mask
            folded |= mask
            owner[mask] = unit.index
        index = np.flatnonzero(folded)
        pairs = np.array(sorted(cKDTree(structure.xyz[index]).query_pairs(1.0)))
        if pairs.size:
            cross = owner[index[pairs[:, 0]]] != owner[index[pairs[:, 1]]]
            assert not cross.any(), "a unit was placed through another unit"

    def test_cycle_closing_cannot_introduce_inter_unit_clashes(self) -> None:
        """A cycle move is subject to the same clash contract as its initial placement."""
        structure = read_structure(COMPLEX_FIXTURES / "fold_yeast_med2_med3_med15_model_0.cif")
        assign_regions(structure)
        assembly = find_rigid_units(structure)
        report = reposition_folded_domains(
            structure,
            rng=np.random.default_rng(5),
            assembly=assembly,
        )

        for unit in assembly.units:
            obstacles = np.concatenate(
                [other.xyz for other in assembly.units if other.index != unit.index], axis=0
            )
            if obstacles.size:
                assert not cKDTree(obstacles).query_ball_point(
                    unit.xyz, CA_CLASH_DISTANCE, return_length=True
                ).any()
        assert report.ok, report.summary()
        assert not [unit for unit in report.units if unit.clashing]


class TestReporting:
    def test_a_dictated_linker_says_so(self) -> None:
        """A linker inside one rigid unit has a span DODO did not choose; that has to be said."""
        structure = _assigned("dnmt3a_dimer.pdb")
        report = reposition_folded_domains(
            structure,
            rng=np.random.default_rng(0),
            assembly=find_rigid_units(structure, lock_chains=True),
        )
        dictated = [linker for linker in report.linkers if linker.dictated]
        assert dictated
        assert "dictated by the complex" in str(dictated[0])

    def test_interfaces_are_reported_whether_or_not_they_locked(self) -> None:
        structure = _assigned("dnmt3a_dimer.pdb")
        report = reposition_folded_domains(structure, rng=np.random.default_rng(0))
        assert any(i.locked for i in report.interfaces)
        assert any(not i.locked for i in report.interfaces)
        assert report.broken_interfaces == [i for i in report.interfaces if not i.preserved]

    def test_summary_names_the_multi_domain_units(self) -> None:
        structure = _assigned("dnmt3a_dimer.pdb")
        report = reposition_folded_domains(structure, rng=np.random.default_rng(0))
        assert "rigid unit(s) hold more than one folded domain" in report.summary()

    def test_ok_is_true_on_a_clean_run(self) -> None:
        structure = _assigned("dnmt3a.pdb")
        assert reposition_folded_domains(structure, rng=np.random.default_rng(0)).ok

    def test_an_empty_report_is_ok(self) -> None:
        assert PlacementReport().ok


class TestVerifyRigid:
    def test_a_scaled_copy_is_caught(self) -> None:
        before = np.random.default_rng(0).normal(size=(50, 3))
        with pytest.raises(GeometryError, match="deformed rather than moved"):
            verify_rigid(before, before * 1.01)

    def test_a_translation_and_rotation_pass(self) -> None:
        from dodo.geometry.transforms import rotation_from_axis_angle

        before = np.random.default_rng(0).normal(size=(50, 3))
        rotation = rotation_from_axis_angle(np.array([1.0, 2.0, 3.0]), 0.7)
        verify_rigid(before, before @ rotation.T + np.array([10.0, 0.0, -4.0]))

    def test_mismatched_shapes_are_refused(self) -> None:
        with pytest.raises(GeometryError, match="shapes differ"):
            verify_rigid(np.zeros((3, 3)), np.zeros((4, 3)))

    def test_equal_centroid_radii_do_not_hide_a_rearrangement(self) -> None:
        before_angles = np.deg2rad([0.0, 90.0, 180.0, 270.0])
        after_angles = np.deg2rad([30.0, -30.0, 150.0, 210.0])
        before = np.column_stack(
            [np.cos(before_angles), np.sin(before_angles), np.zeros(before_angles.size)]
        )
        after = np.column_stack(
            [np.cos(after_angles), np.sin(after_angles), np.zeros(after_angles.size)]
        )
        assert np.allclose(
            np.linalg.norm(before - before.mean(axis=0), axis=1),
            np.linalg.norm(after - after.mean(axis=0), axis=1),
        )
        with pytest.raises(GeometryError, match="deformed rather than moved"):
            verify_rigid(before, after)
