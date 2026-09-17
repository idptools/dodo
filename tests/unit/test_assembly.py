"""Tests for rigid units: which folded domains must move together.

The invariant these exist to protect is one sentence: **two atoms in the same rigid unit are the
same distance apart after a rebuild as before it.** Everything else here is a route to breaking
that, and each route has a named test, because each of the three was a real defect measured on
real data before rigid units existed:

* a complex was dismantled -- two copies of dnmt3a in contact lost an interface of 5,540
  atom-atom contacts, and the rebuild reported ``ok``;
* two folded domains with no residues between them were separated, taking a peptide bond from
  1.334 A to 93.479 A while the report said "their relative position is not ours to change";
* a linker too short to rebuild was stretched, taking a peptide bond to 33.605 A.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dodo.constants import INTERFACE_MIN_RESIDUE_PAIRS
from dodo.construct.assembly import (
    REASON_ADJACENT,
    REASON_UNBUILT_LINK,
    RigidUnit,
    find_rigid_units,
    units_from_spec,
    verify_units_rigid,
)
from dodo.exceptions import GeometryError, InvalidRegionError
from dodo.io import read_structure
from dodo.regions.identify import assign_regions, assign_regions_from_spec
from dodo.structure import DomainKind, Structure

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "structures"


def _assigned(name: str):
    structure = read_structure(FIXTURES / name)
    assign_regions(structure)
    return structure


class TestSingleChainIsUnchanged:
    """A single chain must come out exactly as it did before rigid units existed.

    The behaviour being protected was validated over 23,587 structures, so "one unit per folded
    domain on a single chain" is not a nicety -- it is what makes the whole mechanism safe to
    have added.
    """

    @pytest.mark.parametrize("name", ["dnmt3a.pdb", "p300.pdb", "arf19.pdb"])
    def test_one_unit_per_folded_domain(self, name: str) -> None:
        structure = _assigned(name)
        assembly = find_rigid_units(structure)
        folded = [d for d in structure.domains if d.kind is DomainKind.FOLDED]
        assert len(assembly.units) == len(folded)
        assert all(len(unit) == 1 for unit in assembly.units)

    def test_intra_chain_contacts_are_reported_but_not_locked(self) -> None:
        """Measuring a contact and acting on it are different decisions.

        A user needs to see a contact DODO is prepared to break -- on a single chain breaking it
        is the entire point of step 3 -- so the contact is reported either way.
        """
        structure = _assigned("dnmt3a.pdb")
        assembly = find_rigid_units(structure)
        assert all(not i.inter_chain for i in assembly.interfaces)
        assert all(not i.locked for i in assembly.interfaces)

    def test_locking_intra_chain_contacts_is_available(self) -> None:
        structure = _assigned("dnmt3a.pdb")
        opted_in = find_rigid_units(structure, lock_intra_chain_interfaces=True)
        touching = [
            i for i in opted_in.interfaces if i.residue_pairs >= INTERFACE_MIN_RESIDUE_PAIRS
        ]
        if touching:
            assert len(opted_in.units) < len(find_rigid_units(structure).units)


class TestComplexes:
    def test_an_interface_across_chains_locks(self) -> None:
        """The headline case: the two chains of a dimer end up in one unit."""
        structure = _assigned("dnmt3a_dimer.pdb")
        assembly = find_rigid_units(structure)
        multi = [unit for unit in assembly.units if len(unit) > 1]
        assert multi, "the dimer's interface did not lock anything together"
        spanning = [unit for unit in multi if len(unit.chain_ids()) > 1]
        assert spanning, "no unit spans both chains, so the interface can still be broken"

    def test_the_big_interface_is_locked_and_a_brush_is_not(self) -> None:
        """The threshold has to discriminate, or it is not doing anything.

        Measured on this fixture: the interface is 344 residue contacts and the incidental
        touches are 2 and 3.
        """
        structure = _assigned("dnmt3a_dimer.pdb")
        assembly = find_rigid_units(structure)
        locked = [i for i in assembly.interfaces if i.locked]
        broken = [i for i in assembly.interfaces if not i.locked]
        assert locked and broken
        assert min(i.residue_pairs for i in locked) >= INTERFACE_MIN_RESIDUE_PAIRS
        assert all(
            i.residue_pairs < INTERFACE_MIN_RESIDUE_PAIRS or not i.inter_chain for i in broken
        )

    def test_an_assembly_of_one_domain_per_chain_is_one_unit(self) -> None:
        """6kn7 is an actin filament: everything touches everything, so nothing may move."""
        structure = _assigned("6kn7.pdb")
        assembly = find_rigid_units(structure)
        assert len(assembly.units) == 1
        assert len(assembly.units[0]) > 20

    def test_transitively_joined_interface_is_not_reported_broken(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A below-threshold edge survives when another path puts both domains in one unit."""
        import dodo.construct.assembly as assembly_module

        structure = Structure.from_atom_records(
            xyz=np.array(
                [[0.0, 0.0, 0.0], [20.0, 0.0, 0.0], [40.0, 0.0, 0.0], [43.81, 0.0, 0.0]]
            ),
            atom_name=["CA"] * 4,
            element=["C"] * 4,
            residue_name=["ALA"] * 4,
            residue_number=[1, 1, 1, 2],
            chain_id=["A", "B", "C", "C"],
            source="transitive-unit probe",
        )
        assign_regions_from_spec(
            structure,
            {
                "A": [("folded", 1, 1)],
                "B": [("folded", 1, 1)],
                "C": [("folded", 1, 1), ("folded", 2, 2)],
            },
        )

        def contacts(*_args: object, **_kwargs: object) -> dict[tuple[int, int], tuple[int, int]]:
            return {(0, 2): (10, 10), (0, 3): (1, 1)}

        monkeypatch.setattr(assembly_module, "_contacts_between_domains", contacts)
        assembly = find_rigid_units(structure)
        below_threshold = next(i for i in assembly.interfaces if i.residue_pairs == 1)
        assert not below_threshold.locked
        assert below_threshold.preserved
        assert below_threshold not in assembly.broken_interfaces
        unit = assembly.unit_of(structure.chains[0].domains[0])
        assert any(REASON_ADJACENT in reason for reason in unit.reasons)
        assert any("interface" in reason for reason in unit.reasons)


class TestCovalentRules:
    """Two folded domains DODO must never separate, whatever the contacts say."""

    def test_adjacent_folded_domains_are_one_unit(self) -> None:
        """The 93.479 A peptide bond.

        Reachable from the documented ``assign_regions_from_spec`` path, which accepts two
        folded domains with no residues between them.
        """
        structure = read_structure(FIXTURES / "dnmt3a.pdb")
        assign_regions_from_spec(
            structure,
            {
                "A": [
                    ("idr", 1, 282),
                    ("folded", 283, 432),
                    ("idr", 433, 473),
                    ("folded", 474, 700),
                    ("folded", 701, 912),
                ]
            },
        )
        assembly = find_rigid_units(structure)
        second = next(d for d in structure.domains if d.span.start == 473)
        third = next(d for d in structure.domains if d.span.start == 700)
        assert assembly.same_unit(second, third)
        assert any(REASON_ADJACENT in reason for reason in assembly.unit_of(second).reasons)

    def test_a_link_that_will_not_be_rebuilt_is_one_unit(self) -> None:
        """The 33.605 A peptide bond.

        ``min_length`` is a public ``rebuild`` argument, and a region under it keeps its input
        coordinates -- which makes it a rigid link between the domains it joins.
        """
        structure = _assigned("dnmt3a.pdb")
        linker = next(
            d for d in structure.domains if d.kind is DomainKind.IDR and not d.span.is_terminal
        )
        assembly = find_rigid_units(structure, min_length=len(linker.span) + 1)
        flanking = [
            d
            for d in structure.domains
            if (d.kind is DomainKind.FOLDED
            and d.span.stop == linker.span.start)
            or d.span.start == linker.span.stop
        ]
        assert len(flanking) == 2
        assert assembly.same_unit(flanking[0], flanking[1])
        assert any(REASON_UNBUILT_LINK in r for r in assembly.unit_of(flanking[0]).reasons)

    def test_a_buildable_link_does_not_lock(self) -> None:
        """The mirror of the test above: this is exactly the case step 3 exists for."""
        structure = _assigned("dnmt3a.pdb")
        assembly = find_rigid_units(structure, min_length=4)
        folded = [d for d in structure.domains if d.kind is DomainKind.FOLDED]
        assert not assembly.same_unit(folded[0], folded[1])


class TestOptions:
    def test_lock_chains_holds_every_chain_rigid(self) -> None:
        structure = _assigned("dnmt3a_dimer.pdb")
        assembly = find_rigid_units(structure, lock_chains=True)
        # Both chains touch, so holding each chain rigid and then locking the interface leaves
        # exactly one unit.
        assert len(assembly.units) == 1

    def test_lock_interfaces_false_keeps_only_the_covalent_rules(self) -> None:
        structure = _assigned("dnmt3a_dimer.pdb")
        assembly = find_rigid_units(structure, lock_interfaces=False)
        assert all(len(unit) == 1 for unit in assembly.units)
        assert assembly.interfaces == []

    def test_units_from_spec_groups_what_the_caller_names(self) -> None:
        structure = _assigned("dnmt3a_dimer.pdb")
        assembly = units_from_spec(structure, [[("A", 300), ("B", 300)]])
        a = next(d for d in structure.chains[0].domains if 300 - 1 in d.span)
        b = next(
            d
            for d in structure.chains[1].domains
            if structure.residue_number[d.span.start]
            <= 300
            <= structure.residue_number[d.span.stop - 1]
        )
        assert assembly.same_unit(a, b)

    def test_units_from_spec_rejects_a_residue_outside_a_folded_domain(self) -> None:
        structure = _assigned("dnmt3a.pdb")
        with pytest.raises(InvalidRegionError, match="not inside a folded domain"):
            units_from_spec(structure, [[("A", 5)]])

    def test_units_from_spec_rejects_an_unknown_chain(self) -> None:
        structure = _assigned("dnmt3a.pdb")
        with pytest.raises(InvalidRegionError, match="not in this structure"):
            units_from_spec(structure, [[("Z", 300)]])


class TestUnitTransforms:
    """A unit moves as one body, or the whole mechanism is decoration."""

    def test_rotation_preserves_every_inter_domain_distance(self) -> None:
        from dodo.geometry.transforms import rotation_from_axis_angle

        structure = _assigned("dnmt3a_dimer.pdb")
        assembly = find_rigid_units(structure)
        unit = next(u for u in assembly.units if len(u) > 1)
        before = unit.xyz.copy()
        unit.rotate(rotation_from_axis_angle(np.array([0.3, 0.4, 0.5]), 1.1), about=unit.centroid())
        unit.translate([12.0, -3.0, 7.5])
        after = unit.xyz
        # Compare a sample of pairs spanning different domains, which is the point.
        rng = np.random.default_rng(0)
        pick = rng.choice(before.shape[0], size=2000, replace=False)
        d0 = np.linalg.norm(before[pick][:, None] - before[pick][None, :100], axis=-1)
        d1 = np.linalg.norm(after[pick][:, None] - after[pick][None, :100], axis=-1)
        assert np.abs(d1 - d0).max() < 1e-9

    def test_rotate_requires_an_explicit_centre(self) -> None:
        """Rotating each domain about its own centroid is the bug this class exists to prevent."""
        structure = _assigned("dnmt3a_dimer.pdb")
        unit = next(u for u in find_rigid_units(structure).units if len(u) > 1)
        with pytest.raises(TypeError):
            unit.rotate(np.eye(3))  # type: ignore[call-arg]

    def test_verify_units_rigid_catches_a_per_domain_rotation(self) -> None:
        """The check has to be on the unit; a per-domain check passes on the real defect."""
        from dodo.construct.place import verify_rigid
        from dodo.geometry.transforms import rotation_from_axis_angle

        structure = _assigned("dnmt3a_dimer.pdb")
        assembly = find_rigid_units(structure)
        unit = next(u for u in assembly.units if len(u) > 1)
        before = {unit.index: unit.xyz.copy()}
        rotation = rotation_from_axis_angle(np.array([0.0, 0.0, 1.0]), 0.4)
        for domain in unit.domains:
            domain.rotate(rotation, about=domain.centroid())  # each about its OWN centroid

        # Every domain is individually intact...
        for domain in unit.domains:
            verify_rigid(domain.xyz, domain.xyz)
        # ...and the unit is not.
        with pytest.raises(GeometryError, match="deformed rather than moved"):
            verify_units_rigid(before, assembly)


class TestDeterminism:
    def test_units_are_ordered_by_their_earliest_domain(self) -> None:
        structure = _assigned("dnmt3a_dimer.pdb")
        assembly = find_rigid_units(structure)
        starts = [unit.domains[0].span.start for unit in assembly.units]
        assert starts == sorted(starts)

    def test_repeated_calls_agree(self) -> None:
        structure = _assigned("dnmt3a_dimer.pdb")
        first = find_rigid_units(structure)
        second = find_rigid_units(structure)
        assert [tuple(d.span.start for d in u.domains) for u in first.units] == [
            tuple(d.span.start for d in u.domains) for u in second.units
        ]


class TestEdges:
    def test_no_folded_domains_gives_no_units(self) -> None:
        from dodo.structure import Domain

        structure = read_structure(FIXTURES / "dnmt3a.pdb")
        chain = structure.chains[0]
        chain.domains = [Domain(structure=structure, span=chain.span, kind=DomainKind.IDR)]
        assembly = find_rigid_units(structure)
        assert assembly.units == []
        assert "nothing to hold rigid" in " ".join(assembly.notes)

    def test_unit_of_rejects_a_domain_it_does_not_own(self) -> None:
        structure = _assigned("dnmt3a.pdb")
        assembly = find_rigid_units(structure)
        idr = next(d for d in structure.domains if d.kind is DomainKind.IDR)
        with pytest.raises(InvalidRegionError, match="not part of this assembly"):
            assembly.unit_of(idr)

    def test_an_empty_unit_is_rejected(self) -> None:
        with pytest.raises(InvalidRegionError, match="at least one folded domain"):
            RigidUnit(index=0, domains=[])

    def test_summary_mentions_every_unit(self) -> None:
        structure = _assigned("dnmt3a_dimer.pdb")
        assembly = find_rigid_units(structure)
        text = assembly.summary()
        assert all(f"unit {unit.index}" in text for unit in assembly.units)
