"""Tests for the constants module.

Mostly guards on internal consistency. The pre-rewrite code carried four divergent
copies of the build-mode table (with three different values for ``super_compact``),
three different CA-CA bond lengths, and two disagreeing atomic-mass tables. These tests
assert that there is now one of each and that the derived helpers agree with them.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from dodo import constants as C


class TestModes:
    def test_all_v1_mode_names_survive(self) -> None:
        """The user-facing vocabulary is preserved even though the numbers changed.

        Users and the README are trained on these names; breaking them buys nothing.
        """
        for name in (
            "super_compact",
            "compact",
            "normal",
            "expanded",
            "super_expanded",
            "max_expansion",
            "predicted",
        ):
            assert name in C.MODES

    def test_modes_are_monotonic(self) -> None:
        ordered = [
            "super_compact",
            "compact",
            "normal",
            "expanded",
            "super_expanded",
            "max_expansion",
        ]
        factors = [C.MODES[m] for m in ordered]
        assert factors == sorted(factors), "compactness ordering must be monotonic"

    def test_predicted_is_unity(self) -> None:
        """'predicted' means 'exactly what ALBATROSS says'."""
        assert C.MODES["predicted"] == 1.0

    def test_normal_and_predicted_coincide(self) -> None:
        """Documented behaviour change from v1, asserted so it stays deliberate.

        In v1 these were different: 'normal' was 0.8 A/residue and 'predicted' called
        ALBATROSS. Rebasing the modes onto the prediction makes them the same thing.
        """
        assert C.MODES["normal"] == C.MODES["predicted"]

    def test_resolve_mode(self) -> None:
        assert C.resolve_mode("compact") == C.MODES["compact"]

    def test_unknown_mode_lists_the_valid_ones(self) -> None:
        with pytest.raises(ValueError, match="Valid modes are") as excinfo:
            C.resolve_mode("very_squished")
        # The message must be actionable, since this is the most common user error.
        assert "super_compact" in str(excinfo.value)
        assert "very_squished" in str(excinfo.value)

    def test_default_mode_is_valid(self) -> None:
        assert C.DEFAULT_MODE in C.MODES


class TestDimensionEstimates:
    def test_flory_matches_the_documented_formula(self) -> None:
        assert C.flory_end_to_end(100) == pytest.approx(
            C.FLORY_RE_PREFACTOR * 100**C.FLORY_RE_EXPONENT
        )

    def test_flory_is_monotonic_in_length(self) -> None:
        values = [C.flory_end_to_end(n) for n in (10, 50, 100, 200, 400)]
        assert values == sorted(values)

    def test_flory_is_sublinear(self) -> None:
        """The whole point of the change: Re must scale sublinearly in N.

        v1 used a constant A/residue multiplier, i.e. an effective Flory exponent of
        1.0, so doubling the length doubled the target. Real chains grow as ~N^0.52, so
        doubling the length multiplies Re by ~1.44.
        """
        ratio = C.flory_end_to_end(400) / C.flory_end_to_end(200)
        assert 1.3 < ratio < 1.5
        assert C.FLORY_RE_EXPONENT < 1.0

    def test_flory_stays_under_contour_length(self) -> None:
        """A predicted Re longer than the chain itself would be unsatisfiable."""
        for n in (5, 10, 50, 100, 500, 1000):
            assert C.flory_end_to_end(n) < C.contour_length(n), f"failed at N={n}"

    def test_contour_length_is_bond_length_times_bonds(self) -> None:
        assert C.contour_length(11) == pytest.approx(10 * C.CA_CA_BOND_LENGTH)

    @pytest.mark.parametrize("bad", [0, -1, -100])
    def test_non_positive_lengths_rejected(self, bad: int) -> None:
        with pytest.raises(ValueError, match="must be positive"):
            C.flory_end_to_end(bad)
        with pytest.raises(ValueError, match="must be positive"):
            C.contour_length(bad)


class TestBackboneAngles:
    def test_sampling_window_sits_inside_the_observed_range(self) -> None:
        assert C.BACKBONE_ANGLE_OBSERVED_MIN <= C.BACKBONE_ANGLE_MIN
        assert C.BACKBONE_ANGLE_MAX <= C.BACKBONE_ANGLE_OBSERVED_MAX

    def test_window_brackets_the_mean_on_both_sides(self) -> None:
        """The window spans the bulk of the distribution but not its tails.

        Both bounds have the SAME provenance -- the author's window, tuned against the
        measured AF2 angle distribution -- so both are asserted the same way, against that
        distribution's mean and sd. The upper bound was briefly re-derived from what all-atom
        reconstruction can represent (capped at 150.0) and that was reverted, so a symmetric
        sd-based assertion is the right property for it again; see
        :meth:`test_window_is_measured_rather_than_capped_for_all_atom_reconstruction`.

        Neither bound is a round number of sd and the window is not symmetric about the
        mean, so the band is deliberately loose: tidying it to mean +/- 1 sd would narrow it
        to 99-151 and drop helical geometry.
        """
        low_sd = (C.BACKBONE_ANGLE_MEAN - C.BACKBONE_ANGLE_MIN) / C.BACKBONE_ANGLE_SD
        high_sd = (C.BACKBONE_ANGLE_MAX - C.BACKBONE_ANGLE_MEAN) / C.BACKBONE_ANGLE_SD
        assert 1.0 <= low_sd <= 2.0
        assert 1.0 <= high_sd <= 2.0
        assert C.BACKBONE_ANGLE_MIN < C.BACKBONE_ANGLE_MEAN < C.BACKBONE_ANGLE_MAX

    def test_window_is_measured_rather_than_capped_for_all_atom_reconstruction(self) -> None:
        """The window is the measured distribution, NOT the all-atom-reconstructable range.

        This replaces an assertion that the upper bound sits *under* the all-atom ceiling.
        That policy is reversed: the bound was narrowed to 150.0 so that every generated
        angle would be reconstructable, and it is back at the author's tuned 161.0. The
        priority order is explicit -- (1) the CA-only pipeline, (2) backbones, (3) all-atom,
        (4) STARLING -- and narrowing the *generation* window to serve (3) degrades (1), by
        making chains more compact than the measured distribution supports. The all-atom path
        has to grow its own constraint instead, and a magnitude cap was never that constraint:
        the real one is on the *change* in pseudo-angle between consecutive residues (measured
        allowance ~35 deg at a base angle of 95, collapsing to ~4 deg at 150).

        So the property asserted here is the consequence, recorded rather than forgotten:
        part of the generation window is NOT all-atom reconstructable. Both ceilings are
        derived from :mod:`dodo.construct.backbone` rather than written down, so this test
        follows any change to the peptide geometry or to the tau tolerance.
        """
        from dodo.construct.backbone import _N_CA_C_TOLERANCE, max_reconstructable_ca_angle

        grid = C.backbone_angle_grid()

        # At the tolerance place_backbone actually builds with, the top of the window is over
        # the ceiling -- so the window admits angles no trans backbone can realize.
        tolerated = max_reconstructable_ca_angle(n_ca_c_tolerance=_N_CA_C_TOLERANCE)
        assert tolerated < C.BACKBONE_ANGLE_MAX
        beyond_tolerated = int(np.count_nonzero(grid > tolerated))
        assert beyond_tolerated > 0
        # ... but only the very top of it, which is why the CA-only cost of capping the
        # window is not worth paying for the all-atom benefit.
        assert beyond_tolerated < grid.size

        # At an *ideal* N-CA-C the unreconstructable share is substantial, and that is the
        # honest statement of the coupling: it is not a rounding-error overlap.
        ideal = max_reconstructable_ca_angle()
        assert ideal < C.BACKBONE_ANGLE_MAX
        assert int(np.count_nonzero(grid > ideal)) > beyond_tolerated
        # The coupling bites at the wide end only; helical geometry is unaffected.
        assert ideal > C.BACKBONE_ANGLE_MIN

    def test_ideal_angle_is_inside_the_window(self) -> None:
        assert C.BACKBONE_ANGLE_MIN <= C.BACKBONE_ANGLE_IDEAL <= C.BACKBONE_ANGLE_MAX

    def test_angle_grid_spans_the_window(self) -> None:
        grid = C.backbone_angle_grid()
        assert grid.min() == C.BACKBONE_ANGLE_MIN
        assert grid.max() == C.BACKBONE_ANGLE_MAX
        assert grid.size == C.BACKBONE_ANGLE_MAX - C.BACKBONE_ANGLE_MIN + 1

    def test_angle_grid_is_ordered_ideal_first(self) -> None:
        """Ideal-first ordering is what makes 'first non-clashing' prefer real geometry.

        The pre-rewrite vectorized generator lost this, emitting 91, 92, 93..., which
        silently discarded the bias toward physically likely angles.
        """
        grid = C.backbone_angle_grid()
        assert grid[0] == C.BACKBONE_ANGLE_IDEAL
        deviations = np.abs(grid - C.BACKBONE_ANGLE_IDEAL)
        assert np.all(np.diff(deviations) >= 0), "must be sorted by distance from ideal"

    def test_angle_grid_has_no_duplicates(self) -> None:
        grid = C.backbone_angle_grid()
        assert grid.size == np.unique(grid).size

    def test_minimum_angle_excludes_physically_impossible_kinks(self) -> None:
        """The window must forbid the fold-backs the unconstrained walk produced.

        Measured on the pre-rewrite IDR builder: CA-CA-CA angles as sharp as 47 deg,
        where the only constraint was a 3.0 A clash cutoff. That permits
        2*asin(3.0/7.6) = 46.5 deg, and no real trace goes below ~75.
        """
        clash_permitted = math.degrees(2 * math.asin(3.0 / (2 * C.CA_CA_BOND_LENGTH)))
        assert clash_permitted < C.BACKBONE_ANGLE_MIN
        assert C.BACKBONE_ANGLE_MIN >= C.BACKBONE_ANGLE_OBSERVED_MIN


class TestClashDistances:
    def test_relaxation_ladder_is_descending(self) -> None:
        ladder = C.CLASH_RELAXATION_LADDER
        assert list(ladder) == sorted(ladder, reverse=True)

    def test_ladder_starts_at_the_nominal_clash_distance(self) -> None:
        assert C.CLASH_RELAXATION_LADDER[0] == C.CA_CLASH_DISTANCE

    def test_clash_distance_is_below_bond_length(self) -> None:
        """Otherwise every bonded pair registers as a clash regardless of exclusion."""
        assert C.CA_CLASH_DISTANCE < C.CA_CA_BOND_LENGTH


class TestElementData:
    def test_single_mass_table(self) -> None:
        for element in ("C", "N", "O", "S", "P", "H"):
            assert element in C.ATOMIC_MASSES

    def test_masses_are_plausible(self) -> None:
        assert C.ATOMIC_MASSES["C"] == pytest.approx(12.011, abs=0.01)
        assert C.ATOMIC_MASSES["O"] == pytest.approx(15.999, abs=0.01)

    def test_three_to_one_covers_all_standard_residues(self) -> None:
        for one in "ACDEFGHIKLMNPQRSTVWY":
            three = C.ONE_TO_THREE[one]
            assert C.THREE_TO_ONE[three] == one

    def test_modified_residues_map_to_their_parent(self) -> None:
        """Dropping these fabricates phantom chain breaks."""
        assert C.THREE_TO_ONE["MSE"] == "M"
        assert C.THREE_TO_ONE["SEP"] == "S"
        assert C.THREE_TO_ONE["PTR"] == "Y"

    def test_heavy_atom_counts_cover_standard_residues(self) -> None:
        for one in "ACDEFGHIKLMNPQRSTVWY":
            assert C.ONE_TO_THREE[one] in C.HEAVY_ATOM_COUNTS

    def test_glycine_has_fewest_heavy_atoms(self) -> None:
        """Underpins why the raw contact score was composition-biased."""
        assert C.HEAVY_ATOM_COUNTS["GLY"] == min(C.HEAVY_ATOM_COUNTS.values())
        assert C.HEAVY_ATOM_COUNTS["TRP"] == max(C.HEAVY_ATOM_COUNTS.values())


class TestRegionThresholds:
    def test_gap_and_domain_length_are_independent(self) -> None:
        """These were one knob in v1 and are two distinct quantities.

        They may hold the same value, but they must be separately settable.
        """
        assert C.MAX_INTERNAL_GAP is not None
        assert C.MIN_FOLDED_DOMAIN_LENGTH is not None

    def test_loop_radius_is_tighter_than_contact_radius(self) -> None:
        assert C.LOOP_CONTACT_RADIUS < C.CONTACT_RADIUS

    def test_min_idr_length_permits_geometry(self) -> None:
        """An IDR must have enough residues for a bond angle to be defined."""
        assert C.MIN_IDR_LENGTH >= 3

    def test_plddt_threshold_is_in_range(self) -> None:
        assert 0.0 < C.PLDDT_DISORDER_THRESHOLD < 100.0


class TestBackboneBondGeometry:
    def test_bond_lengths_are_physical(self) -> None:
        """v1 declared these correctly and then never imported them.

        Its all-atom module added unit vectors without scaling, producing CA-C at
        exactly 1.000 A and a 2.87 A peptide bond.
        """
        assert 1.4 < C.N_CA_BOND_LENGTH < 1.5
        assert 1.5 < C.CA_C_BOND_LENGTH < 1.6
        assert 1.2 < C.C_O_BOND_LENGTH < 1.3
        assert 1.3 < C.C_N_PEPTIDE_BOND_LENGTH < 1.4

    def test_backbone_angles_are_tetrahedral_ish(self) -> None:
        for angle in (C.N_CA_C_ANGLE, C.CA_C_N_ANGLE, C.C_N_CA_ANGLE, C.CA_C_O_ANGLE):
            assert 105.0 < angle < 125.0

    def test_virtual_ca_distance_agrees_with_backbone_geometry(self) -> None:
        """The CA-CA virtual bond must be consistent with the atoms bridging it.

        CA-C-N-CA with trans omega gives a CA-CA distance near 3.8 A. If these two sets
        of constants disagreed, all-atom reconstruction would not close onto the CA
        trace it was built from.
        """
        # Build CA1-C-N-CA2 in a plane, which is what omega = 180 (trans) means.
        ca1 = np.array([0.0, 0.0, 0.0])
        c = np.array([C.CA_C_BOND_LENGTH, 0.0, 0.0])

        # C->CA1 points along 180 deg, so C->N sits CA_C_N_ANGLE away from it.
        theta = math.radians(180.0 - C.CA_C_N_ANGLE)
        n = c + C.C_N_PEPTIDE_BOND_LENGTH * np.array([math.cos(theta), math.sin(theta), 0.0])

        # From N->C, turn by C_N_CA_ANGLE to continue the extended zig-zag.
        direction = (c - n) / np.linalg.norm(c - n)
        phi = math.atan2(direction[1], direction[0]) + math.radians(C.C_N_CA_ANGLE)
        ca2 = n + C.N_CA_BOND_LENGTH * np.array([math.cos(phi), math.sin(phi), 0.0])

        assert np.linalg.norm(ca2 - ca1) == pytest.approx(C.CA_CA_BOND_LENGTH, abs=0.05)


class TestStarlingLimits:
    def test_max_length_is_the_observed_cap(self) -> None:
        assert C.STARLING_MAX_LENGTH == 380

    def test_splice_overlap_is_smaller_than_a_segment(self) -> None:
        assert 0 < C.SEGMENT_SPLICE_OVERLAP < C.STARLING_MAX_LENGTH
