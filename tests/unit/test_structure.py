"""Tests for the struct-of-arrays structure core.

A number of these are explicit regression tests for bugs confirmed in the pre-rewrite
object-graph implementation. Each is named for the failure it prevents, because the
whole justification for the SoA rewrite is that it makes these impossible by
construction -- and that claim is worth holding to.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from dodo.constants import CA_CA_BOND_LENGTH
from dodo.exceptions import EmptyStructureError, GeometryError, InvalidRegionError
from dodo.structure import Domain, DomainKind, Span, Structure

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def make_structure(
    n_residues: int = 10,
    *,
    chain_id: str = "A",
    first_residue_number: int = 1,
    all_atom: bool = True,
    spacing: float = CA_CA_BOND_LENGTH,
) -> Structure:
    """Build a straight-line poly-alanine structure for testing.

    Deliberately parameterised on ``first_residue_number``: several confirmed bugs only
    appeared on structures not numbered from 1, which is to say on every crystal and EM
    structure, so tests need to exercise both.
    """
    atom_names = ["N", "CA", "C", "O"] if all_atom else ["CA"]
    elements = ["N", "C", "C", "O"] if all_atom else ["C"]

    xyz, names, els, resnames, resnums, chains = [], [], [], [], [], []
    for i in range(n_residues):
        for offset, (name, element) in enumerate(zip(atom_names, elements, strict=True)):
            # CA sits exactly on the residue's lattice point; other atoms are nudged so
            # they are distinguishable but do not collide.
            jitter = 0.0 if name == "CA" else 0.4 * (offset + 1)
            xyz.append([i * spacing, jitter, 0.0])
            names.append(name)
            els.append(element)
            resnames.append("ALA")
            resnums.append(first_residue_number + i)
            chains.append(chain_id)

    return Structure.from_atom_records(
        xyz=np.array(xyz),
        atom_name=names,
        element=els,
        residue_name=resnames,
        residue_number=resnums,
        chain_id=chains,
        b_factor=[85.0] * len(names),
        source="test-fixture",
    )


@pytest.fixture
def structure() -> Structure:
    return make_structure(10)


# ---------------------------------------------------------------------------
# Span
# ---------------------------------------------------------------------------


class TestSpan:
    def test_length_is_half_open(self) -> None:
        assert len(Span(3, 7)) == 4

    def test_membership_excludes_stop(self) -> None:
        span = Span(3, 7)
        assert 3 in span
        assert 6 in span
        assert 7 not in span

    def test_iterates_residue_indices(self) -> None:
        assert list(Span(2, 5)) == [2, 3, 4]

    def test_empty_span_is_allowed(self) -> None:
        # A zero-length span is meaningful: a folded domain with no loops has none.
        assert len(Span(4, 4)) == 0

    def test_inverted_span_rejected(self) -> None:
        with pytest.raises(InvalidRegionError, match="precedes start"):
            Span(7, 3)

    def test_negative_start_rejected(self) -> None:
        with pytest.raises(InvalidRegionError, match="non-negative"):
            Span(-1, 5)

    def test_terminal_when_missing_an_anchor(self) -> None:
        assert Span(0, 5, c_anchor=5).is_terminal
        assert Span(5, 9, n_anchor=4).is_terminal
        assert not Span(5, 9, n_anchor=4, c_anchor=9).is_terminal

    def test_overlap_detection(self) -> None:
        assert Span(0, 5).overlaps(Span(4, 8))
        # Adjacent half-open spans do not overlap. This is the property that lets
        # regions tile a chain without the shared-boundary-residue convention the
        # pre-rewrite code used, which needed three different fixups to resolve.
        assert not Span(0, 5).overlaps(Span(5, 8))


# ---------------------------------------------------------------------------
# Structure construction and invariants
# ---------------------------------------------------------------------------


class TestConstruction:
    def test_shape(self, structure: Structure) -> None:
        assert structure.n_residues == 10
        assert structure.n_atoms == 40

    def test_sequence(self, structure: Structure) -> None:
        assert structure.sequence == "A" * 10

    def test_offsets_are_consistent(self, structure: Structure) -> None:
        assert structure.residue_atom_offsets[0] == 0
        assert structure.residue_atom_offsets[-1] == structure.n_atoms
        assert np.all(np.diff(structure.residue_atom_offsets) == 4)

    def test_validate_passes_on_well_formed(self, structure: Structure) -> None:
        structure.validate()  # must not raise

    def test_insertion_codes_keep_residues_distinct(self) -> None:
        """Residues 10 and 10A are two residues, not one.

        The pre-rewrite reader keyed residues on ``int(resSeq)`` alone, so an insertion
        code merged two residues into one object holding two N atoms and two CAs, and
        silently discarded 10A's identity.
        """
        s = Structure.from_atom_records(
            xyz=np.arange(24, dtype=float).reshape(8, 3),
            atom_name=["N", "CA", "C", "O"] * 2,
            element=["N", "C", "C", "O"] * 2,
            residue_name=["GLY"] * 4 + ["ALA"] * 4,
            residue_number=[10] * 8,
            insertion_code=[""] * 4 + ["A"] * 4,
            chain_id=["A"] * 8,
        )
        assert s.n_residues == 2
        assert s.residue_number.tolist() == [10, 10]
        assert s.insertion_code.tolist() == ["", "A"]

    def test_modified_residues_map_to_parent_amino_acid(self) -> None:
        """MSE is a methionine, not an unknown.

        The pre-rewrite reader selected only records starting with 'ATAM'/'ATOM', so a
        mid-chain selenomethionine (which is a HETATM) vanished entirely and fabricated
        a phantom chain break where the polymer is in fact continuous.
        """
        s = Structure.from_atom_records(
            xyz=np.zeros((2, 3)),
            atom_name=["CA", "CA"],
            element=["C", "C"],
            residue_name=["MSE", "SEC"],
            residue_number=[1, 2],
            chain_id=["A", "A"],
        )
        assert s.sequence == "MU"

    def test_multiple_chains_are_split(self) -> None:
        s = Structure.from_atom_records(
            xyz=np.zeros((4, 3)),
            atom_name=["CA"] * 4,
            element=["C"] * 4,
            residue_name=["ALA"] * 4,
            residue_number=[1, 2, 1, 2],
            chain_id=["A", "A", "B", "B"],
        )
        assert [c.chain_id for c in s.chains] == ["A", "B"]
        assert s.chain_index.tolist() == [0, 0, 1, 1]

    def test_interleaved_chain_records_keep_all_atoms(self) -> None:
        """An A,B,A,B interleaved file must not lose atoms.

        The pre-rewrite reader assigned ``chain_dict[id] = lines`` rather than
        extending, so the second run of each chain overwrote the first -- roughly half
        the atoms gone, with no warning. Here every atom survives, the repeated
        identifiers are treated as distinct residues, and the anomaly is recorded in
        ``notes`` so it is visible rather than silent.
        """
        s = Structure.from_atom_records(
            xyz=np.zeros((4, 3)),
            atom_name=["CA"] * 4,
            element=["C"] * 4,
            residue_name=["ALA"] * 4,
            residue_number=[1, 1, 2, 2],
            chain_id=["A", "B", "A", "B"],
        )
        assert s.n_atoms == 4, "no atoms may be dropped"
        assert len(s.chains) == 4, "each run becomes its own chain span"

    def test_repeated_residue_identifier_is_noted(self) -> None:
        s = Structure.from_atom_records(
            xyz=np.zeros((3, 3)),
            atom_name=["CA"] * 3,
            element=["C"] * 3,
            residue_name=["ALA"] * 3,
            residue_number=[1, 2, 1],
            chain_id=["A"] * 3,
        )
        assert any("more than one run" in note for note in s.notes)

    def test_empty_input_rejected(self) -> None:
        with pytest.raises(EmptyStructureError):
            Structure.from_atom_records(
                xyz=np.zeros((0, 3)),
                atom_name=[],
                element=[],
                residue_name=[],
                residue_number=[],
                chain_id=[],
            )

    def test_mismatched_array_lengths_rejected(self) -> None:
        with pytest.raises(GeometryError, match="entries but xyz has"):
            Structure.from_atom_records(
                xyz=np.zeros((3, 3)),
                atom_name=["CA", "CA"],  # one short
                element=["C"] * 3,
                residue_name=["ALA"] * 3,
                residue_number=[1, 2, 3],
                chain_id=["A"] * 3,
            )

    def test_validate_catches_corrupted_offsets(self, structure: Structure) -> None:
        structure.residue_atom_offsets[3] = 999
        with pytest.raises(GeometryError):
            structure.validate()

    def test_missing_ca_reports_which_residues(self) -> None:
        s = Structure.from_atom_records(
            xyz=np.zeros((2, 3)),
            atom_name=["N", "C"],  # no CA
            element=["N", "C"],
            residue_name=["ALA", "GLY"],
            residue_number=[1, 2],
            chain_id=["A", "A"],
        )
        with pytest.raises(EmptyStructureError, match="no CA atom"):
            _ = s.ca_xyz


# ---------------------------------------------------------------------------
# The bug class the SoA design is meant to eliminate
# ---------------------------------------------------------------------------


class TestViewsCannotGoStale:
    """Regression tests for the pre-rewrite cache-invalidation bugs.

    Six confirmed bugs existed only because coordinates were duplicated across four
    cache levels and exactly one of the four containment links participated in
    invalidation. These tests assert the property that replaces the whole protocol.
    """

    def test_domain_xyz_is_a_view(self, structure: Structure) -> None:
        domain = Domain(structure, Span(0, 4), DomainKind.FOLDED)
        assert np.shares_memory(domain.xyz, structure.xyz)
        assert domain.xyz.flags["OWNDATA"] is False

    def test_parent_write_visible_through_domain(self, structure: Structure) -> None:
        domain = Domain(structure, Span(0, 4), DomainKind.FOLDED)
        structure.xyz[0] = [1.0, 2.0, 3.0]
        # No refresh, no invalidation call. There is one copy of the coordinate.
        assert domain.xyz[0].tolist() == [1.0, 2.0, 3.0]

    def test_domain_write_visible_through_parent(self, structure: Structure) -> None:
        domain = Domain(structure, Span(0, 4), DomainKind.FOLDED)
        domain.translate([10.0, 0.0, 0.0])
        assert structure.xyz[0][0] == pytest.approx(10.0)

    def test_two_overlapping_views_agree(self, structure: Structure) -> None:
        """Two views covering the same residue see the same coordinates.

        In the object-graph version, ``Domain.translate`` updated the domain's cache
        while the ``Monomer`` and ``Complex`` caches kept reporting the old position.
        """
        wide = Domain(structure, Span(0, 8), DomainKind.FOLDED)
        narrow = Domain(structure, Span(2, 4), DomainKind.FOLDED)
        wide.translate([5.0, 5.0, 5.0])
        assert np.array_equal(wide.ca_xyz[2:4], narrow.ca_xyz)

    def test_anchor_is_read_from_the_anchor_residue(self, structure: Structure) -> None:
        """Anchors come from the named anchor residue, not "the domain's last atom".

        The pre-rewrite ``get_end_coord()`` returned the last atom of the monomer dict,
        which broke once rebuilt loop residues were appended out of order: it returned
        an interior loop residue, and the next region was anchored 36 A from where it
        belonged, producing two atoms 0.00 A apart.
        """
        idr = Domain(structure, Span(4, 7), DomainKind.IDR)
        idr.span = Span(4, 7, n_anchor=3, c_anchor=7)
        assert idr.n_terminal_ca is not None
        assert np.array_equal(idr.n_terminal_ca, structure.ca_xyz[3])
        assert idr.c_terminal_ca is not None
        assert np.array_equal(idr.c_terminal_ca, structure.ca_xyz[7])

    def test_terminal_span_has_no_anchor_on_the_free_side(self, structure: Structure) -> None:
        tail = Domain(structure, Span(0, 3), DomainKind.IDR)
        tail.span = Span(0, 3, c_anchor=3)
        assert tail.n_terminal_ca is None
        assert tail.c_terminal_ca is not None

    def test_copy_is_independent(self, structure: Structure) -> None:
        chain = structure.chains[0]
        chain.domains.append(Domain(structure, Span(0, 10), DomainKind.FOLDED))
        clone = structure.copy()
        clone.xyz[0] = [99.0, 99.0, 99.0]
        assert structure.xyz[0][0] != pytest.approx(99.0)
        # Views in the clone must be rebound, or mutating the clone silently mutates
        # the original through a dangling reference.
        assert clone.chains[0].structure is clone
        assert clone.chains[0].domains[0].structure is clone
        clone.validate()


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


class TestGeometry:
    def test_translate_preserves_internal_distances(self, structure: Structure) -> None:
        before = np.linalg.norm(structure.ca_xyz[0] - structure.ca_xyz[5])
        structure.translate([13.0, -4.0, 7.5])
        after = np.linalg.norm(structure.ca_xyz[0] - structure.ca_xyz[5])
        assert after == pytest.approx(before)

    def test_rotation_preserves_internal_distances(self, structure: Structure) -> None:
        rotation = Rotation.from_euler("xyz", [17.0, 43.0, 88.0], degrees=True).as_matrix()
        before = np.linalg.norm(structure.ca_xyz[0] - structure.ca_xyz[9])
        structure.rotate(rotation)
        after = np.linalg.norm(structure.ca_xyz[0] - structure.ca_xyz[9])
        assert after == pytest.approx(before)

    def test_rotation_about_centroid_leaves_centroid_fixed(self, structure: Structure) -> None:
        rotation = Rotation.from_euler("z", 90.0, degrees=True).as_matrix()
        before = structure.xyz.mean(axis=0)
        structure.rotate(rotation)
        assert structure.xyz.mean(axis=0) == pytest.approx(before)

    def test_domain_rotation_defaults_to_its_own_centroid(self, structure: Structure) -> None:
        """A domain rotates about itself, not about the world origin.

        The pre-rewrite ``Chain.rotate`` and ``Complex.rotate`` applied
        ``coords @ R.T`` with no centring, so they translated as well as rotated -- and
        were inconsistent with ``Domain.rotate``, which did centre.
        """
        domain = Domain(structure, Span(0, 5), DomainKind.FOLDED)
        before = domain.centroid()
        domain.rotate(Rotation.from_euler("y", 45.0, degrees=True).as_matrix())
        assert domain.centroid() == pytest.approx(before)

    def test_non_3x3_rotation_rejected(self, structure: Structure) -> None:
        with pytest.raises(GeometryError, match=r"\(3, 3\)"):
            structure.rotate(np.eye(4))

    def test_mass_weighted_centroid_differs_from_geometric(self, structure: Structure) -> None:
        """The two are distinct quantities and must not be conflated.

        The pre-rewrite code had one method computing a mass-weighted centre of mass and
        another computing an unweighted centroid, both named as if they did the same
        thing, and rotated domains about whichever it happened to call.
        """
        domain = Domain(structure, Span(0, 10), DomainKind.FOLDED)
        assert domain.centroid(mass_weighted=True) != pytest.approx(
            domain.centroid(mass_weighted=False)
        )

    def test_empty_domain_centroid_raises(self, structure: Structure) -> None:
        domain = Domain(structure, Span(4, 4), DomainKind.FOLDED)
        with pytest.raises(GeometryError, match="empty domain"):
            domain.centroid()

    def test_end_to_end_distance_of_a_straight_chain(self) -> None:
        s = make_structure(11, all_atom=False)
        assert s.end_to_end_distance() == pytest.approx(10 * CA_CA_BOND_LENGTH)

    def test_radius_of_gyration_is_positive(self, structure: Structure) -> None:
        assert structure.radius_of_gyration() > 0

    def test_metrics_need_two_residues(self) -> None:
        s = make_structure(1, all_atom=False)
        with pytest.raises(GeometryError, match="at least 2"):
            s.end_to_end_distance()


# ---------------------------------------------------------------------------
# Clash detection
# ---------------------------------------------------------------------------


class TestClashDetection:
    def test_empty_obstacle_set_means_no_clash(self, structure: Structure) -> None:
        """The state at the start of every rebuild must not raise.

        The pre-rewrite equivalent built a KD-tree over an empty coordinate list and
        raised ``ValueError: data must be of shape (n, m)`` -- which is precisely the
        state before anything has been placed.
        """
        empty = np.zeros(structure.n_atoms, dtype=bool)
        result = structure.clash_mask(np.zeros((3, 3)), obstacle_atom_mask=empty)
        assert result.tolist() == [False, False, False]

    def test_detects_a_real_clash(self, structure: Structure) -> None:
        query = structure.ca_xyz[5] + np.array([0.1, 0.0, 0.0])
        assert structure.clash_mask(query[None, :])[0]

    def test_distant_point_does_not_clash(self, structure: Structure) -> None:
        assert not structure.clash_mask(np.array([[1000.0, 1000.0, 1000.0]]))[0]

    def test_bonded_neighbours_are_not_clashes(self, structure: Structure) -> None:
        """Adjacent residues sit at 3.8 A by construction and are not a collision.

        Without this exclusion every peptide bond registers as a clash, which is exactly
        what the pre-rewrite whole-structure clash check reported.
        """
        query = structure.ca_xyz[[5]]
        clashing = structure.clash_mask(
            query, query_residue_index=[5], cutoff=4.5, exclude_within=2
        )
        assert not clashing[0]

    def test_without_exclusion_neighbours_do_register(self, structure: Structure) -> None:
        query = structure.ca_xyz[[5]]
        assert structure.clash_mask(query, cutoff=4.5)[0]

    def test_exclusion_still_catches_distant_residues(self) -> None:
        """A residue folding back onto one far away in sequence is a real clash."""
        s = make_structure(20, all_atom=False)
        # Move residue 19 on top of residue 0.
        s.set_ca_xyz([19], s.ca_xyz[[0]] + 0.2)
        clashing = s.clash_mask(
            s.ca_xyz[[19]], query_residue_index=[19], cutoff=3.2, exclude_within=2
        )
        assert clashing[0]

    def test_rebuilt_mask_selects_only_rebuilt_domains(self, structure: Structure) -> None:
        chain = structure.chains[0]
        placed = Domain(structure, Span(0, 5), DomainKind.FOLDED, rebuilt=True)
        pending = Domain(structure, Span(5, 10), DomainKind.IDR, rebuilt=False)
        chain.domains = [placed, pending]
        mask = structure.rebuilt_atom_mask()
        assert mask[:20].all()
        assert not mask[20:].any()

    def test_query_index_length_must_match(self, structure: Structure) -> None:
        with pytest.raises(GeometryError, match="query points"):
            structure.clash_mask(np.zeros((3, 3)), query_residue_index=[1])

    def test_malformed_query_shape_rejected(self, structure: Structure) -> None:
        with pytest.raises(GeometryError, match=r"\(n, 3\)"):
            structure.clash_mask(np.zeros((3, 2)))


# ---------------------------------------------------------------------------
# Domain bookkeeping
# ---------------------------------------------------------------------------


class TestChainDomainValidation:
    def test_tiling_chain_validates(self, structure: Structure) -> None:
        chain = structure.chains[0]
        chain.domains = [
            Domain(structure, Span(0, 4), DomainKind.FOLDED),
            Domain(structure, Span(4, 7), DomainKind.IDR),
            Domain(structure, Span(7, 10), DomainKind.FOLDED),
        ]
        chain.validate_domains()

    def test_gap_rejected(self, structure: Structure) -> None:
        chain = structure.chains[0]
        chain.domains = [
            Domain(structure, Span(0, 4), DomainKind.FOLDED),
            Domain(structure, Span(6, 10), DomainKind.FOLDED),
        ]
        with pytest.raises(InvalidRegionError, match="unassigned residues"):
            chain.validate_domains()

    def test_overlap_rejected(self, structure: Structure) -> None:
        chain = structure.chains[0]
        chain.domains = [
            Domain(structure, Span(0, 6), DomainKind.FOLDED),
            Domain(structure, Span(4, 10), DomainKind.FOLDED),
        ]
        with pytest.raises(InvalidRegionError, match="overlapping"):
            chain.validate_domains()

    def test_no_domains_rejected(self, structure: Structure) -> None:
        with pytest.raises(InvalidRegionError, match="no domains"):
            structure.chains[0].validate_domains()

    def test_loop_outside_its_domain_rejected(self, structure: Structure) -> None:
        chain = structure.chains[0]
        chain.domains = [Domain(structure, Span(0, 10), DomainKind.FOLDED, loops=(Span(8, 12),))]
        with pytest.raises(InvalidRegionError, match="outside its domain"):
            chain.validate_domains()

    def test_domain_lookup_by_residue(self, structure: Structure) -> None:
        chain = structure.chains[0]
        folded = Domain(structure, Span(0, 5), DomainKind.FOLDED)
        idr = Domain(structure, Span(5, 10), DomainKind.IDR)
        chain.domains = [folded, idr]
        assert chain.domain_at(2) is folded
        assert chain.domain_at(7) is idr
        assert chain.domain_at(99) is None


# ---------------------------------------------------------------------------
# The index convention
# ---------------------------------------------------------------------------


class TestIndexConvention:
    """PDB numbering is data; positional indices are indices. They never cross.

    Crossing these two spaces was the single largest bug generator in the pre-rewrite
    code: region bounds were positional while the residue container was keyed on PDB
    ``resSeq``, so selection was off by one on an ordinary AF2 file and matched nothing
    at all on any structure not numbered from 1.
    """

    @pytest.mark.parametrize("first", [1, 0, 42, -5, 1000])
    def test_positional_indexing_is_independent_of_pdb_numbering(self, first: int) -> None:
        s = make_structure(10, first_residue_number=first, all_atom=False)
        domain = Domain(s, Span(2, 5), DomainKind.IDR)
        # Three residues selected regardless of how the file numbers them.
        assert len(domain) == 3
        assert domain.ca_xyz.shape == (3, 3)
        # And PDB numbering is preserved for output.
        assert s.residue_number[2] == first + 2

    def test_numbering_gaps_do_not_shift_selection(self) -> None:
        """A structure with missing residues still selects by position correctly."""
        resnums = [1, 2, 3, 50, 51, 52]  # a 46-residue gap
        s = Structure.from_atom_records(
            xyz=np.arange(18, dtype=float).reshape(6, 3),
            atom_name=["CA"] * 6,
            element=["C"] * 6,
            residue_name=["ALA"] * 6,
            residue_number=resnums,
            chain_id=["A"] * 6,
        )
        domain = Domain(s, Span(3, 6), DomainKind.FOLDED)
        assert len(domain) == 3
        assert s.residue_number[domain.span.slice].tolist() == [50, 51, 52]

    def test_out_of_range_span_rejected(self, structure: Structure) -> None:
        with pytest.raises(InvalidRegionError, match="outside this structure"):
            structure.atom_slice_for_residues(0, 999)

    def test_residue_label_includes_insertion_code(self) -> None:
        s = Structure.from_atom_records(
            xyz=np.zeros((2, 3)),
            atom_name=["CA", "CA"],
            element=["C", "C"],
            residue_name=["GLU", "GLU"],
            residue_number=[142, 142],
            insertion_code=["", "B"],
            chain_id=["A", "A"],
        )
        assert s.residue_label(0) == "A:GLU142"
        assert s.residue_label(1) == "A:GLU142B"
