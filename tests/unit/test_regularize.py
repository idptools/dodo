"""Tests for CA bond-length and pseudo-angle regularization.

This module exists because STARLING is a diffusion model: it denoises coordinates toward a
learned distribution and nothing in that process enforces a hard geometric constraint, so
neither its CA-CA distances nor its CA-CA-CA pseudo-angles are a backbone's. The correction has
to happen after generation.

The tests that matter most are the conformation-preservation ones. Making bonds exact is easy
-- a naive walk-and-replace does it in one pass. Making them exact *without changing the
conformation* is the whole problem, because the conformation is the thing a generative model
was used for in the first place.

The angle tests exist because bond repair alone measured as *not enough*: on real STARLING
output at 201 residues, after every bond was projected onto 3.81 A exactly, 2.34% of vertices
were still outside the observed 75-179 deg range and the median conformer had 5 of them, so the
engine's screen accepted 0 of 20 conformers and the engine produced nothing. The angle tests
that matter most are the ones about the two constraints *fighting*: they are constraints on the
same atoms, so anything that satisfies one by abandoning the other is a regression this file has
to catch.
"""

from __future__ import annotations

import numpy as np
import pytest

from dodo.constants import (
    BACKBONE_ANGLE_OBSERVED_MAX,
    BACKBONE_ANGLE_OBSERVED_MIN,
    CA_CA_BOND_LENGTH,
)
from dodo.exceptions import GeometryError
from dodo.geometry.metrics import (
    ca_bond_lengths,
    ca_pseudo_angles,
    end_to_end,
    radius_of_gyration,
)
from dodo.geometry.regularize import needs_regularization, regularize_ca_trace

#: The window DODO's STARLING screen enforces, and therefore the one worth repairing into.
OBSERVED = (BACKBONE_ANGLE_OBSERVED_MIN, BACKBONE_ANGLE_OBSERVED_MAX)


def diffusion_like_trace(n: int, seed: int = 0, bond_sigma: float = 0.15) -> np.ndarray:
    """Build a plausible coil whose bond lengths scatter, as a diffusion model's output does.

    Directionally correlated so the result is a coil rather than white noise, with bond lengths
    drawn around the target rather than fixed to it. ``bond_sigma`` is the scatter; 0.15 A is a
    modest, realistic amount and 0.4 A is deliberately harsh.
    """
    rng = np.random.default_rng(seed)
    points = [np.zeros(3), np.array([CA_CA_BOND_LENGTH, 0.0, 0.0])]
    for _ in range(n - 2):
        random_direction = rng.normal(size=3)
        random_direction /= np.linalg.norm(random_direction)
        previous = points[-1] - points[-2]
        previous /= np.linalg.norm(previous)
        direction = 0.6 * previous + 0.4 * random_direction
        direction /= np.linalg.norm(direction)
        points.append(points[-1] + direction * rng.normal(CA_CA_BOND_LENGTH, bond_sigma))
    return np.array(points)


def mds_like_trace(n: int, seed: int = 0) -> np.ndarray:
    """Build a trace with STARLING's *two* defects at once: short bonds and sharp vertices.

    :func:`diffusion_like_trace` is directionally correlated, so it almost never doubles back
    and its angles are mostly fine -- which makes it the wrong fixture for angle repair. What
    STARLING actually produces comes out of an MDS fit to a predicted distance map, which
    constrains neither the bond nor the vertex: MEASURED, a 201-residue ensemble had bonds with
    a median of 3.21 A (16% short) and pseudo-angles reaching 10.5 deg.

    So this reproduces both. A weak directional correlation lets the walk reverse on itself and
    the bond lengths are drawn around 0.85x the target, giving a trace with real out-of-window
    vertices at both ends of the range rather than a tidy coil.
    """
    rng = np.random.default_rng(seed)
    points = [np.zeros(3), np.array([CA_CA_BOND_LENGTH, 0.0, 0.0])]
    for _ in range(n - 2):
        random_direction = rng.normal(size=3)
        random_direction /= np.linalg.norm(random_direction)
        previous = points[-1] - points[-2]
        previous /= np.linalg.norm(previous)
        direction = 0.15 * previous + 0.85 * random_direction
        direction /= np.linalg.norm(direction)
        points.append(points[-1] + direction * rng.normal(0.85 * CA_CA_BOND_LENGTH, 0.35))
    return np.array(points)


def chain_with_angles(angles_degrees: list[float], torsion_degrees: float = 100.0) -> np.ndarray:
    """Build a trace with bonds exactly ``CA_CA_BOND_LENGTH`` and exactly these pseudo-angles.

    Internal-coordinate construction, so the fixture states the geometry under test directly
    instead of hoping a random walk produces it. That matters for the angle tests: the bond pass
    has *nothing* to do on this input, so any bond error or atom movement observed afterwards is
    attributable to the angle projection alone.

    ``angles_degrees`` has one entry per vertex, i.e. two fewer than the residue count.
    """
    bond = CA_CA_BOND_LENGTH
    first = np.radians(angles_degrees[0])
    points = [
        np.zeros(3),
        np.array([bond, 0.0, 0.0]),
        np.array([bond, 0.0, 0.0]) + bond * np.array([-np.cos(first), np.sin(first), 0.0]),
    ]
    torsion = np.radians(torsion_degrees)
    for theta_degrees in angles_degrees[1:]:
        theta = np.radians(theta_degrees)
        a, b, c = points[-3], points[-2], points[-1]
        bc = (c - b) / np.linalg.norm(c - b)
        normal = np.cross(b - a, bc)
        normal /= np.linalg.norm(normal)
        frame = np.column_stack([bc, np.cross(normal, bc), normal])
        offset = bond * np.array(
            [-np.cos(theta), np.sin(theta) * np.cos(torsion), np.sin(theta) * np.sin(torsion)]
        )
        points.append(c + frame @ offset)
    return np.array(points)


def pinned_infeasible_vertex() -> np.ndarray:
    """Three residues whose only pseudo-angle is forced outside the observed window.

    With both termini pinned and both bonds at 3.81 A, the 1-3 distance determines the angle
    outright: ``d = 2 b sin(theta / 2)``. Endpoints 4.0 A apart therefore force
    ``theta = 2 asin(4.0 / 7.62) = 63.3 deg``, which is 11.7 deg below the observed floor and
    which nothing can move. The middle atom starts on the perpendicular bisector so both bonds
    are already exact and the bond constraint has no complaint of its own.
    """
    half = 2.0
    height = float(np.sqrt(CA_CA_BOND_LENGTH**2 - half**2))
    return np.array([[0.0, 0.0, 0.0], [half, height, 0.0], [2.0 * half, 0.0, 0.0]])


class TestExactness:
    @pytest.mark.parametrize("n", [3, 10, 50, 150, 380])
    def test_bonds_become_exact(self, n: int) -> None:
        result = regularize_ca_trace(diffusion_like_trace(n, seed=n))
        assert result.converged
        bonds = ca_bond_lengths(result.ca_coords)
        assert np.allclose(bonds, CA_CA_BOND_LENGTH, atol=1e-6)

    def test_reports_the_error_it_started_with(self) -> None:
        result = regularize_ca_trace(diffusion_like_trace(100, seed=1))
        assert result.max_bond_error_before > 0.1
        assert result.max_bond_error_after < 1e-6

    def test_an_already_exact_trace_is_left_alone(self) -> None:
        exact = np.zeros((20, 3))
        exact[:, 0] = CA_CA_BOND_LENGTH * np.arange(20)
        result = regularize_ca_trace(exact)
        assert result.converged
        assert result.rmsd == pytest.approx(0.0, abs=1e-9)

    def test_heavy_noise_still_converges(self) -> None:
        result = regularize_ca_trace(diffusion_like_trace(200, seed=7, bond_sigma=0.4))
        assert result.converged
        assert result.max_bond_error_before > 1.0
        assert result.max_bond_error_after < 1e-6

    def test_custom_bond_length_is_honoured(self) -> None:
        result = regularize_ca_trace(diffusion_like_trace(30, seed=2), bond_length=3.81)
        assert np.allclose(ca_bond_lengths(result.ca_coords), 3.81, atol=1e-6)


class TestConformationPreservation:
    """The tests that distinguish a projection from a rebuild."""

    @pytest.mark.parametrize("n", [50, 150, 380])
    def test_radius_of_gyration_is_preserved(self, n: int) -> None:
        """Rg must barely move. This is the headline requirement.

        A naive rebuild -- walk from residue 0 placing each CA at exactly the bond length
        along the direction to the next -- fixes every bond but accumulates error down the
        chain, so the far end drifts and the overall dimensions change. That would throw away
        the conformation, which is the only reason to run a generative model at all.
        """
        result = regularize_ca_trace(diffusion_like_trace(n, seed=n))
        assert abs(result.rg_change_fraction) < 0.01, result.summary()

    def test_end_to_end_is_preserved(self) -> None:
        result = regularize_ca_trace(diffusion_like_trace(150, seed=5))
        relative = abs(result.end_to_end_after - result.end_to_end_before)
        assert relative / result.end_to_end_before < 0.01, result.summary()

    def test_displacement_is_small(self) -> None:
        """Each atom moves a fraction of a bond length, not a whole one."""
        result = regularize_ca_trace(diffusion_like_trace(150, seed=6))
        assert result.rmsd < 0.5, result.summary()

    def test_beats_a_naive_rebuild_on_shape_preservation(self) -> None:
        """Explicit comparison against the obvious wrong approach.

        The naive rebuild is what this module exists instead of, so the advantage should be
        demonstrated rather than asserted in prose.
        """
        trace = diffusion_like_trace(200, seed=11, bond_sigma=0.3)
        rg_original = radius_of_gyration(trace)

        # Naive: keep each bond DIRECTION, force each bond LENGTH, accumulate from residue 0.
        naive = np.zeros_like(trace)
        naive[0] = trace[0]
        for i in range(1, trace.shape[0]):
            direction = trace[i] - trace[i - 1]
            direction /= np.linalg.norm(direction)
            naive[i] = naive[i - 1] + direction * CA_CA_BOND_LENGTH
        naive_error = abs(radius_of_gyration(naive) - rg_original) / rg_original

        projected = regularize_ca_trace(trace)
        assert abs(projected.rg_change_fraction) < naive_error
        # Both fix the bonds; only one keeps the shape.
        assert np.allclose(ca_bond_lengths(naive), CA_CA_BOND_LENGTH, atol=1e-9)

    def test_angles_are_not_wrecked(self) -> None:
        """Correcting bonds must move angles a little, but not reshape the chain.

        Some angle change is unavoidable and correct: moving atoms to fix bond lengths
        necessarily perturbs the angles between them. What matters is that the perturbation is
        bounded.
        """
        trace = diffusion_like_trace(200, seed=3)
        before = ca_pseudo_angles(trace)
        after = ca_pseudo_angles(regularize_ca_trace(trace).ca_coords)
        assert float(np.abs(after - before).max()) < 15.0


class TestAngleRepairPremise:
    """The fixtures have the defect the angle tests claim to repair.

    Without these, an angle test that trivially passes because its input was already in range
    looks exactly like a working repair.
    """

    @pytest.mark.parametrize("n", [40, 120, 201])
    def test_the_mds_fixture_really_has_out_of_window_vertices(self, n: int) -> None:
        angles = ca_pseudo_angles(mds_like_trace(n, seed=n))
        outside = np.count_nonzero((angles < OBSERVED[0]) | (angles > OBSERVED[1]))
        assert outside > 0, f"fixture has no bad vertex to repair; angles {angles.min():.1f}+"

    def test_bond_repair_alone_leaves_bad_vertices_behind(self) -> None:
        """The measurement that makes angle repair necessary rather than nice.

        On real STARLING output this was the whole problem: bonds came out exact and 2.34% of
        vertices were still unreconstructable, so a screen that rejects a conformer for one bad
        vertex rejected every conformer. Reproduced here on a fixture so the reason this code
        exists is pinned in the suite and not only in a commit message.
        """
        bonds_only = regularize_ca_trace(mds_like_trace(201, seed=201))
        assert np.allclose(ca_bond_lengths(bonds_only.ca_coords), CA_CA_BOND_LENGTH, atol=1e-6)
        angles = ca_pseudo_angles(bonds_only.ca_coords)
        outside = np.count_nonzero((angles < OBSERVED[0]) | (angles > OBSERVED[1]))
        assert outside > 0
        # And it is not reported, because nothing asked about angles.
        assert bonds_only.angle_window is None
        assert bonds_only.max_angle_violation_after == 0.0

    def test_the_exact_builder_produces_the_angles_it_was_given(self) -> None:
        requested = [120.0, 40.0, 150.0, 95.0, 178.5, 60.0]
        trace = chain_with_angles(requested)
        assert ca_pseudo_angles(trace) == pytest.approx(requested, abs=1e-9)
        assert np.allclose(ca_bond_lengths(trace), CA_CA_BOND_LENGTH, atol=1e-12)


class TestAngleRepair:
    """Every vertex ends up inside the window, and the bonds stay exact while it happens."""

    @pytest.mark.parametrize("n", [3, 10, 40, 120, 201])
    def test_every_vertex_ends_inside_the_window(self, n: int) -> None:
        result = regularize_ca_trace(mds_like_trace(n, seed=n), angle_window=OBSERVED)
        angles = ca_pseudo_angles(result.ca_coords)
        assert float(angles.min()) >= OBSERVED[0], result.summary()
        assert float(angles.max()) <= OBSERVED[1], result.summary()
        assert result.max_angle_violation_after == 0.0
        assert result.converged, result.summary()

    @pytest.mark.parametrize("n", [10, 40, 120, 201])
    def test_bonds_are_still_exact_after_angle_repair(self, n: int) -> None:
        """The constraint the angle repair is most likely to have traded away.

        Bending a vertex is done by moving its neighbours, which are exactly the atoms the bond
        constraint holds at 3.81 A. A repair that satisfied the angles by leaving bonds at
        3.79 would pass every angle test in this class and be useless, because the writer and
        the validator both check bonds.
        """
        result = regularize_ca_trace(mds_like_trace(n, seed=n), angle_window=OBSERVED)
        assert np.allclose(ca_bond_lengths(result.ca_coords), CA_CA_BOND_LENGTH, atol=1e-6)
        assert result.max_bond_error_after < 1e-6, result.summary()

    @pytest.mark.parametrize("sharp", [5.0, 10.5, 30.0, 60.0, 74.0])
    def test_a_single_sharp_vertex_is_opened_out(self, sharp: float) -> None:
        """A chain that doubles back on itself within one step, which is what 10.5 deg means.

        Parameterized down to 5 deg because the measured worst real vertex was 10.5 and a fixed
        floor invites a projection that only works on mild cases.
        """
        trace = chain_with_angles([120.0, sharp, 120.0])
        result = regularize_ca_trace(trace, angle_window=OBSERVED)
        assert result.converged, result.summary()
        assert float(ca_pseudo_angles(result.ca_coords).min()) >= OBSERVED[0]
        assert result.max_angle_violation_before == pytest.approx(OBSERVED[0] - sharp)

    @pytest.mark.parametrize("wide", [179.5, 179.9, 180.0])
    def test_a_near_straight_vertex_is_bent_back(self, wide: float) -> None:
        """The case a 1-3 distance projection cannot do, and 180.0 is the fixed point itself.

        With the angle written as ``d = 2 b sin(theta / 2)``, ``dd/dtheta`` vanishes at 180 deg,
        so the correction shrinks to nothing exactly where it is needed and at 180 deg the two
        atoms are collinear and moving them along their separation leaves them collinear.
        MEASURED on real conformers, that parameterization needed 1272 sweeps to close the last
        0.3 deg and left 3 of 20 unrepaired inside a 500-sweep budget. Projecting on the angle
        removes the degeneracy; 180.0 exercises the deterministic-perpendicular fallback.
        """
        trace = chain_with_angles([120.0, wide, 120.0])
        result = regularize_ca_trace(trace, angle_window=OBSERVED)
        assert result.converged, result.summary()
        assert float(ca_pseudo_angles(result.ca_coords).max()) <= OBSERVED[1]
        assert np.allclose(ca_bond_lengths(result.ca_coords), CA_CA_BOND_LENGTH, atol=1e-6)

    def test_both_bounds_violated_in_one_chain(self) -> None:
        """Sharp and straight vertices adjacent to each other, pulling in opposite directions."""
        trace = chain_with_angles([120.0, 20.0, 179.8, 25.0, 179.9, 130.0, 15.0])
        result = regularize_ca_trace(trace, angle_window=OBSERVED)
        assert result.converged, result.summary()
        angles = ca_pseudo_angles(result.ca_coords)
        assert float(angles.min()) >= OBSERVED[0]
        assert float(angles.max()) <= OBSERVED[1]

    def test_a_trace_already_in_range_is_left_alone(self) -> None:
        """An inequality constraint, so a satisfied vertex must not be moved at all."""
        trace = chain_with_angles([100.0, 120.0, 140.0, 110.0, 130.0])
        before = trace.copy()
        result = regularize_ca_trace(trace, angle_window=OBSERVED)
        assert result.rmsd == pytest.approx(0.0, abs=1e-9)
        assert result.ca_coords == pytest.approx(before, abs=1e-9)

    def test_repair_lands_inside_the_window_with_margin_to_spare(self) -> None:
        """The margin is what makes the result robust to a downstream re-measurement.

        Sitting exactly on 75.000 deg would leave the engine's screen deciding acceptance on the
        last bit of a float. The repair aims a margin inside, so a repaired vertex clears the
        bound by a real amount.
        """
        trace = chain_with_angles([120.0, 10.0, 120.0, 179.9, 120.0])
        result = regularize_ca_trace(trace, angle_window=OBSERVED, angle_margin=1.0)
        angles = ca_pseudo_angles(result.ca_coords)
        assert float(angles.min()) > OBSERVED[0] + 0.5
        assert float(angles.max()) < OBSERVED[1] - 0.5

    def test_a_custom_window_is_honoured(self) -> None:
        result = regularize_ca_trace(
            mds_like_trace(80, seed=4), angle_window=(100.0, 150.0), angle_margin=0.5
        )
        angles = ca_pseudo_angles(result.ca_coords)
        assert float(angles.min()) >= 100.0
        assert float(angles.max()) <= 150.0

    def test_a_two_residue_trace_has_no_vertex_to_constrain(self) -> None:
        """No vertex exists, so the window is trivially satisfied rather than an error."""
        trace = np.array([[0.0, 0.0, 0.0], [4.5, 0.0, 0.0]])
        result = regularize_ca_trace(trace, angle_window=OBSERVED)
        assert result.converged
        assert result.max_angle_violation_after == 0.0
        assert result.angle_window == OBSERVED


class TestAngleRepairPreservesConformation:
    """The requirement that makes this a projection and not a reshaping.

    A repair that opened every sharp vertex by inflating the chain would satisfy the window and
    destroy the reason STARLING was called.

    The thresholds are set from the *fixture*, not from real STARLING output, and are looser than
    they look because :func:`mds_like_trace` is deliberately far more damaged than the real thing:
    it is only weakly direction-correlated, so it has many more out-of-window vertices per
    residue than a STARLING conformer does. On real output the added cost of angle repair over
    bond repair MEASURED at 0.01 percentage points of Rg and 0.04 of end-to-end distance; on this
    fixture it reaches 1.10 and 1.58.
    """

    @pytest.mark.parametrize("n", [40, 120, 201])
    def test_angle_repair_costs_almost_nothing_over_bond_repair(self, n: int) -> None:
        trace = mds_like_trace(n, seed=n)
        bonds_only = regularize_ca_trace(trace)
        both = regularize_ca_trace(trace, angle_window=OBSERVED)
        assert abs(both.rg_change_fraction - bonds_only.rg_change_fraction) < 0.02
        assert abs(both.end_to_end_change_fraction - bonds_only.end_to_end_change_fraction) < 0.03

    @pytest.mark.parametrize("n", [40, 120, 201])
    def test_dimensions_survive_the_repair_outright(self, n: int) -> None:
        result = regularize_ca_trace(mds_like_trace(n, seed=n), angle_window=OBSERVED)
        assert abs(result.rg_change_fraction) < 0.05, result.summary()
        assert abs(result.end_to_end_change_fraction) < 0.10, result.summary()

    def test_the_chain_is_not_inflated_into_a_rod(self) -> None:
        """The specific failure mode: satisfying every angle by straightening everything.

        A rod of 201 residues at 3.81 A would be ~762 A end to end. Opening sharp vertices does
        push slightly toward extension, so what is checked is that the result stays a coil.
        """
        trace = mds_like_trace(201, seed=7)
        result = regularize_ca_trace(trace, angle_window=OBSERVED)
        contour = CA_CA_BOND_LENGTH * (trace.shape[0] - 1)
        assert end_to_end(result.ca_coords) < 0.3 * contour, result.summary()
        assert np.median(ca_pseudo_angles(result.ca_coords)) < 160.0

    @pytest.mark.parametrize("n", [40, 120, 201])
    def test_atom_displacement_stays_local(self, n: int) -> None:
        """A fraction of a bond length per atom, as the bond-only projection also manages."""
        trace = mds_like_trace(n, seed=n)
        result = regularize_ca_trace(trace, angle_window=OBSERVED)
        displacement = np.linalg.norm(result.ca_coords - trace, axis=1)
        assert float(np.median(displacement)) < CA_CA_BOND_LENGTH
        assert result.rmsd < CA_CA_BOND_LENGTH, result.summary()


class TestAngleRepairConvergence:
    """How the two competing constraints actually settle, and what is reported when they do not."""

    @pytest.mark.parametrize("n", [40, 120, 201])
    def test_convergence_is_quick_not_merely_eventual(self, n: int) -> None:
        """A budget guard, because the first working version of this was 40x slower.

        The 1-3 distance parameterization converged too -- eventually, at 1272 sweeps for one
        real conformer. That is not a usable engine at 100 conformers a call, and "it converges"
        alone would not have caught it. MEASURED on real STARLING output, the shipped projection
        needed at most 29 joint sweeps and 38 polish sweeps at 201 residues.
        """
        result = regularize_ca_trace(mds_like_trace(n, seed=n), angle_window=OBSERVED)
        assert result.converged
        assert result.sweeps <= 120, result.summary()
        assert result.polish_sweeps <= 120, result.summary()

    def test_running_out_of_sweeps_is_reported_not_hidden(self) -> None:
        with pytest.raises(GeometryError, match="did not converge"):
            regularize_ca_trace(mds_like_trace(201, seed=1), angle_window=OBSERVED, max_sweeps=2)

    def test_the_failure_message_names_the_angle_when_the_angle_is_what_failed(self) -> None:
        """Not merely "did not converge": which of the two constraints lost.

        Driven with the provably infeasible geometry below rather than a starved sweep budget,
        because a starved budget fails on the bond tolerance first and would assert nothing about
        the angle half of the message.
        """
        with pytest.raises(GeometryError, match="pseudo-angle is still"):
            regularize_ca_trace(
                pinned_infeasible_vertex(), angle_window=OBSERVED, preserve_endpoints=True
            )

    def test_an_infeasible_pinned_vertex_is_rejected_rather_than_forced(self) -> None:
        """The one case where the two constraints provably cannot both hold.

        Three residues with both termini pinned: the 1-3 distance is fixed by the anchors and
        both bonds are fixed at 3.81 A, so the single pseudo-angle is fully determined by
        ``d = 2 b sin(theta / 2)`` and no correction exists. This must come back as a refusal
        with the residual measured, not as coordinates that have been forced.
        """
        trace = pinned_infeasible_vertex()
        # Premise: the geometry really is out of the window before anything is attempted.
        assert float(ca_pseudo_angles(trace)[0]) < OBSERVED[0]
        result = regularize_ca_trace(
            trace, angle_window=OBSERVED, preserve_endpoints=True, raise_on_failure=False
        )
        assert not result.converged
        assert result.max_angle_violation_after > 0.0
        # Unmoved, not half-corrected: there was nowhere for it to go.
        assert result.max_angle_violation_after == pytest.approx(
            result.max_angle_violation_before, abs=1e-3
        )
        assert np.array_equal(result.ca_coords[0], trace[0])
        assert np.array_equal(result.ca_coords[-1], trace[-1])

    def test_pinned_endpoints_still_get_bonds_and_angles_right_when_feasible(self) -> None:
        """The interior case DODO actually relies on: closure preserved, geometry repaired."""
        trace = mds_like_trace(120, seed=12)
        result = regularize_ca_trace(trace, angle_window=OBSERVED, preserve_endpoints=True)
        assert result.converged, result.summary()
        assert np.array_equal(result.ca_coords[0], trace[0])
        assert np.array_equal(result.ca_coords[-1], trace[-1])
        assert np.allclose(ca_bond_lengths(result.ca_coords), CA_CA_BOND_LENGTH, atol=1e-6)
        angles = ca_pseudo_angles(result.ca_coords)
        assert float(angles.min()) >= OBSERVED[0]
        assert float(angles.max()) <= OBSERVED[1]

    def test_determinism(self) -> None:
        """No RNG anywhere, including the collinear fallback direction."""
        trace = chain_with_angles([120.0, 180.0, 120.0, 0.5, 120.0])
        first = regularize_ca_trace(trace, angle_window=OBSERVED, raise_on_failure=False)
        second = regularize_ca_trace(trace, angle_window=OBSERVED, raise_on_failure=False)
        assert np.array_equal(first.ca_coords, second.ca_coords)

    def test_input_is_not_mutated(self) -> None:
        trace = mds_like_trace(60, seed=3)
        original = trace.copy()
        regularize_ca_trace(trace, angle_window=OBSERVED)
        assert np.array_equal(trace, original)


class TestAngleWindowValidation:
    """A malformed window must be refused, not silently reinterpreted."""

    @pytest.mark.parametrize(
        "window",
        [(150.0, 100.0), (100.0, 100.0), (-5.0, 100.0), (100.0, 190.0)],
    )
    def test_an_impossible_window_raises(self, window: tuple[float, float]) -> None:
        with pytest.raises(GeometryError, match="0 <= minimum < maximum <= 180"):
            regularize_ca_trace(diffusion_like_trace(20), angle_window=window)

    @pytest.mark.parametrize("window", [(float("nan"), 179.0), (75.0, float("inf"))])
    def test_a_non_finite_window_raises(self, window: tuple[float, float]) -> None:
        with pytest.raises(GeometryError, match="must be finite"):
            regularize_ca_trace(diffusion_like_trace(20), angle_window=window)

    @pytest.mark.parametrize("margin", [-1.0, float("nan")])
    def test_a_bad_margin_raises(self, margin: float) -> None:
        with pytest.raises(GeometryError, match="finite and non-negative"):
            regularize_ca_trace(
                diffusion_like_trace(20), angle_window=OBSERVED, angle_margin=margin
            )

    def test_a_margin_wider_than_the_window_raises(self) -> None:
        """Reject rather than clamp a margin that consumes the whole window.

        Half the width on each side leaves nothing to aim at, and silently clamping would turn a
        caller's mistake into an invisible change of target.
        """
        with pytest.raises(GeometryError, match="leaves nothing"):
            regularize_ca_trace(
                diffusion_like_trace(20), angle_window=(75.0, 179.0), angle_margin=60.0
            )


class TestAngleReporting:
    def test_the_window_is_recorded_so_silence_is_not_ambiguous(self) -> None:
        """Distinguish "angles were fine" from "angles were never checked".

        Both give max_angle_violation_after == 0.0, so the window itself has to be carried.
        """
        checked = regularize_ca_trace(mds_like_trace(60, seed=5), angle_window=OBSERVED)
        unchecked = regularize_ca_trace(mds_like_trace(60, seed=5))
        assert checked.max_angle_violation_after == unchecked.max_angle_violation_after == 0.0
        assert checked.angle_window == OBSERVED
        assert unchecked.angle_window is None

    def test_the_violation_it_started_with_is_reported(self) -> None:
        result = regularize_ca_trace(chain_with_angles([120.0, 10.0, 120.0]), angle_window=OBSERVED)
        assert result.max_angle_violation_before == pytest.approx(65.0)
        assert result.max_angle_violation_after == 0.0

    def test_summary_names_the_angle_numbers_when_angles_were_repaired(self) -> None:
        summary = regularize_ca_trace(mds_like_trace(60, seed=5), angle_window=OBSERVED).summary()
        for token in ("converged", "bond error", "angle excursion", "polish", "Rg", "Re"):
            assert token in summary, summary

    def test_summary_says_nothing_about_angles_when_none_were_requested(self) -> None:
        summary = regularize_ca_trace(mds_like_trace(60, seed=5)).summary()
        assert "angle excursion" not in summary
        assert "polish" not in summary


class TestEndpointPreservation:
    """Pinning the termini, which an interior IDR requires."""

    @pytest.mark.parametrize("n", [10, 50, 380])
    def test_endpoints_do_not_move(self, n: int) -> None:
        """An interior IDR has been placed against two fixed anchors.

        Letting its termini float during regularization would break the closure that placement
        achieved, which is the specific reason this option exists.
        """
        trace = diffusion_like_trace(n, seed=n)
        result = regularize_ca_trace(trace, preserve_endpoints=True)
        assert np.array_equal(result.ca_coords[0], trace[0])
        assert np.array_equal(result.ca_coords[-1], trace[-1])

    def test_end_to_end_is_exactly_unchanged_when_pinned(self) -> None:
        """Both endpoints fixed means the span is fixed by construction."""
        trace = diffusion_like_trace(120, seed=9)
        result = regularize_ca_trace(trace, preserve_endpoints=True)
        assert result.end_to_end_after == pytest.approx(result.end_to_end_before, abs=1e-12)

    def test_pinned_still_reaches_exact_bonds(self) -> None:
        result = regularize_ca_trace(diffusion_like_trace(150, seed=4), preserve_endpoints=True)
        assert result.converged
        assert np.allclose(ca_bond_lengths(result.ca_coords), CA_CA_BOND_LENGTH, atol=1e-6)

    def test_pinned_preserves_shape_too(self) -> None:
        result = regularize_ca_trace(diffusion_like_trace(200, seed=8), preserve_endpoints=True)
        assert abs(result.rg_change_fraction) < 0.01, result.summary()

    def test_two_residues_cannot_be_pinned(self) -> None:
        """With both ends fixed and no interior, there is nothing to absorb the correction."""
        with pytest.raises(GeometryError, match="at least 3 residues"):
            two_residues = np.array([[0.0, 0.0, 0.0], [4.5, 0.0, 0.0]])
            regularize_ca_trace(two_residues, preserve_endpoints=True)


class TestFailureModes:
    def test_non_finite_input_raises_rather_than_being_projected(self) -> None:
        """A NaN means the generator failed; regularizing it would hide that."""
        trace = diffusion_like_trace(20, seed=1)
        trace[5, 1] = np.nan
        with pytest.raises(GeometryError, match="non-finite"):
            regularize_ca_trace(trace)

    def test_coincident_consecutive_residues_raise(self) -> None:
        """There is no bond direction to correct along, and it signals a generator failure."""
        trace = diffusion_like_trace(20, seed=1)
        trace[6] = trace[5]
        with pytest.raises(GeometryError, match="coincident"):
            regularize_ca_trace(trace)

    @pytest.mark.parametrize("bad_shape", [(5, 2), (5,), (2, 3, 3)])
    def test_bad_shape_raises(self, bad_shape: tuple[int, ...]) -> None:
        with pytest.raises(GeometryError, match=r"\(n, 3\)"):
            regularize_ca_trace(np.zeros(bad_shape))

    def test_single_residue_raises(self) -> None:
        with pytest.raises(GeometryError, match="at least 2 residues"):
            regularize_ca_trace(np.zeros((1, 3)))

    @pytest.mark.parametrize("bad", [0.0, -1.0, float("nan"), float("inf")])
    def test_bad_bond_length_raises(self, bad: float) -> None:
        with pytest.raises(GeometryError, match="finite and positive"):
            regularize_ca_trace(diffusion_like_trace(10), bond_length=bad)

    @pytest.mark.parametrize("bad", [0.0, -1e-6, float("nan")])
    def test_bad_tolerance_raises(self, bad: float) -> None:
        with pytest.raises(GeometryError, match="finite and positive"):
            regularize_ca_trace(diffusion_like_trace(10), tolerance=bad)

    def test_non_convergence_raises_by_default(self) -> None:
        """Silence is the failure mode this whole rewrite exists to remove."""
        with pytest.raises(GeometryError, match="did not converge"):
            regularize_ca_trace(diffusion_like_trace(200, seed=1), max_sweeps=1)

    def test_non_convergence_can_be_reported_instead(self) -> None:
        result = regularize_ca_trace(
            diffusion_like_trace(200, seed=1), max_sweeps=1, raise_on_failure=False
        )
        assert not result.converged
        assert result.max_bond_error_after > 1e-6

    def test_max_sweeps_must_be_positive(self) -> None:
        with pytest.raises(GeometryError, match="at least 1"):
            regularize_ca_trace(diffusion_like_trace(10), max_sweeps=0)


class TestPreCheck:
    def test_detects_a_noisy_trace(self) -> None:
        assert needs_regularization(diffusion_like_trace(50, seed=1))

    def test_accepts_a_corrected_trace(self) -> None:
        corrected = regularize_ca_trace(diffusion_like_trace(50, seed=1)).ca_coords
        assert not needs_regularization(corrected)

    def test_uses_the_package_tolerance_not_the_convergence_tolerance(self) -> None:
        """The pre-check asks "is this good enough for DODO", not "is this exact"."""
        exact = np.zeros((10, 3))
        exact[:, 0] = (CA_CA_BOND_LENGTH + 0.01) * np.arange(10)
        # 0.01 A off is inside CA_CA_BOND_TOLERANCE, so no correction is needed.
        assert not needs_regularization(exact)


class TestReporting:
    def test_summary_names_the_numbers_that_matter(self) -> None:
        summary = regularize_ca_trace(diffusion_like_trace(50, seed=1)).summary()
        for token in ("converged", "bond error", "RMSD", "Rg", "Re"):
            assert token in summary

    def test_determinism(self) -> None:
        """No randomness anywhere: the same input gives bit-identical output."""
        trace = diffusion_like_trace(100, seed=1)
        first = regularize_ca_trace(trace).ca_coords
        second = regularize_ca_trace(trace).ca_coords
        assert np.array_equal(first, second)

    def test_input_is_not_mutated(self) -> None:
        trace = diffusion_like_trace(50, seed=1)
        original = trace.copy()
        regularize_ca_trace(trace)
        assert np.array_equal(trace, original)
