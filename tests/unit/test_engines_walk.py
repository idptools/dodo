"""Tests for the engine interface and the self-avoiding walk.

Physical validity is the acceptance criterion here, not "it ran". The geometric
measurements below -- bond lengths, pseudo-angles, non-bonded approach distances -- are
computed with plain numpy rather than with DODO's own metrics, so that a bug shared
between the generator and the metric cannot hide from its own test.

Several tests are named for a specific confirmed defect in the builder being replaced:

* ``test_no_angle_below_the_window`` etc. -- the old walk sampled uniformly on the 3.8 A
  sphere with no angle constraint and produced measured pseudo-angles from 47 to 178
  degrees.
* ``TestClosureScheduleFollowsTheAnchors`` -- the old radius ramp was seeded from the
  residue count, so it fought its own anchors.
* ``TestJunctions`` -- the old builder wrote the first and last residues on top of the
  anchor CAs, giving 0.00 A coincident atoms at every junction.
* ``test_many_conformations`` -- the old code raised ValueError at a reshape for any
  ``n_conformations`` above 1.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from scipy.spatial import cKDTree

from dodo.constants import (
    BACKBONE_ANGLE_MAX,
    BACKBONE_ANGLE_MEAN,
    BACKBONE_ANGLE_MIN,
    BACKBONE_ANGLE_SD,
    CA_CA_BOND_LENGTH,
    CA_CA_BOND_TOLERANCE,
    CA_CLASH_DISTANCE,
    CLASH_EXCLUDE_WITHIN_RESIDUES,
    CLASH_RELAXATION_LADDER,
    MAX_ATTEMPTS_PER_REGION,
    flory_end_to_end,
)
from dodo.engines.base import ConformationEngine, IDRRequest, IDRResult
from dodo.engines.walk import (
    SelfAvoidingWalk,
    _target_steering_width,
    _WalkPlan,
    max_reach,
    min_reach,
    reachability_tail,
    reachable_envelope,
)
from dodo.exceptions import (
    EngineError,
    ExhaustedAttemptsError,
    GeometryError,
    UnsatisfiableTargetError,
)

# The CA(i)-CA(i+2) distance the old builder's 3.0 A clash cutoff permitted, expressed as
# the pseudo-angle it corresponds to: 2 * asin(3.0 / (2 * 3.8)) = 46.5 degrees.
OLD_BUILDER_WORST_ANGLE = 46.5

# The literal reachability table carried by the code being replaced, for n = 6..1 steps,
# measured against a 3.89 A bond length.
OLD_REACHABILITY_TABLE = (15.8868, 13.6407, 11.0914, 8.6688, 6.441, 3.89)
OLD_TABLE_BOND_LENGTH = 3.89


# ---------------------------------------------------------------------------
# Independent geometric measurements
# ---------------------------------------------------------------------------


def bond_lengths(trace: np.ndarray) -> np.ndarray:
    """Consecutive CA-CA distances along a trace."""
    return np.linalg.norm(np.diff(trace, axis=0), axis=1)


def pseudo_angles(trace: np.ndarray) -> np.ndarray:
    """CA(i-1)-CA(i)-CA(i+1) angles in degrees, one per interior point of the trace."""
    back = trace[:-2] - trace[1:-1]
    forward = trace[2:] - trace[1:-1]
    cosine = np.einsum("ij,ij->i", back, forward) / (
        np.linalg.norm(back, axis=1) * np.linalg.norm(forward, axis=1)
    )
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def closest_non_bonded(trace: np.ndarray) -> float:
    """Smallest distance between two trace points more than two residues apart."""
    index = np.arange(trace.shape[0])
    separation = np.abs(np.subtract.outer(index, index))
    distances = np.linalg.norm(trace[:, None, :] - trace[None, :, :], axis=2)
    mask = separation > CLASH_EXCLUDE_WITHIN_RESIDUES
    if not np.any(mask):
        return float("inf")
    return float(np.min(distances[mask]))


def full_trace(
    coords: np.ndarray, n_anchor: np.ndarray | None, c_anchor: np.ndarray | None
) -> np.ndarray:
    """Stack anchors around a generated region, giving the trace a viewer would draw."""
    pieces = []
    if n_anchor is not None:
        pieces.append(np.asarray(n_anchor, dtype=np.float64)[None, :])
    pieces.append(coords)
    if c_anchor is not None:
        pieces.append(np.asarray(c_anchor, dtype=np.float64)[None, :])
    return np.concatenate(pieces, axis=0)


def build(
    n_residues: int,
    *,
    seed: int = 0,
    n_anchor: np.ndarray | None = None,
    c_anchor: np.ndarray | None = None,
    target: float | None = None,
    obstacles: np.ndarray | None = None,
    n_conformations: int = 1,
    engine: SelfAvoidingWalk | None = None,
) -> tuple[IDRRequest, IDRResult]:
    """Build a region and return the request beside the result."""
    request = IDRRequest(
        sequence="G" * n_residues,
        n_residues=n_residues,
        target_end_to_end=flory_end_to_end(n_residues) if target is None else target,
        n_anchor_xyz=n_anchor,
        c_anchor_xyz=c_anchor,
        n_conformations=n_conformations,
    )
    walk = SelfAvoidingWalk() if engine is None else engine
    return request, walk.generate(request, obstacles, np.random.default_rng(seed))


ORIGIN = np.zeros(3)


# ---------------------------------------------------------------------------
# The request and result contracts
# ---------------------------------------------------------------------------


class TestIDRRequest:
    def test_normalizes_anchors_to_float_arrays(self) -> None:
        request = IDRRequest("GG", 2, 5.0, n_anchor_xyz=np.array([1, 2, 3]))
        assert request.n_anchor_xyz is not None
        assert request.n_anchor_xyz.dtype == np.float64

    def test_sequence_length_must_match_residue_count(self) -> None:
        with pytest.raises(EngineError, match="sequence has 3 residues"):
            IDRRequest("GGG", 5, 10.0)

    def test_rejects_non_positive_target(self) -> None:
        with pytest.raises(EngineError, match="positive and finite"):
            IDRRequest("GGG", 3, 0.0)

    def test_rejects_non_finite_anchor(self) -> None:
        with pytest.raises(GeometryError, match="non-finite"):
            IDRRequest("GGG", 3, 10.0, n_anchor_xyz=np.array([np.nan, 0.0, 0.0]))

    def test_rejects_wrong_shaped_anchor(self) -> None:
        with pytest.raises(GeometryError, match="3-vector"):
            IDRRequest("GGG", 3, 10.0, c_anchor_xyz=np.zeros(4))

    def test_rejects_zero_conformations(self) -> None:
        with pytest.raises(EngineError, match="n_conformations"):
            IDRRequest("GGG", 3, 10.0, n_conformations=0)

    def test_interior_when_both_anchors_present(self) -> None:
        interior = IDRRequest(
            "GGG", 3, 10.0, n_anchor_xyz=ORIGIN, c_anchor_xyz=np.array([9.0, 0, 0])
        )
        assert interior.is_interior
        assert interior.anchor_separation == pytest.approx(9.0)

    def test_anchor_separation_is_none_for_a_tail(self) -> None:
        tail = IDRRequest("GGG", 3, 10.0, n_anchor_xyz=ORIGIN)
        assert not tail.is_interior
        assert tail.anchor_separation is None


class TestIDRResult:
    def test_rejects_mismatched_success_length(self) -> None:
        with pytest.raises(EngineError, match="one entry per conformation"):
            IDRResult(np.zeros((2, 4, 3)), np.array([True]), "test", 1)

    def test_rejects_failed_conformer_holding_coordinates(self) -> None:
        # The old builder's exact failure mode: a full array of (0, 0, 0) rows with no
        # indication that anything went wrong. It cannot be expressed through this type.
        with pytest.raises(GeometryError, match="marked failed but hold finite coordinates"):
            IDRResult(np.zeros((2, 4, 3)), np.array([False, False]), "test", 1)

    def test_rejects_successful_conformer_holding_nan(self) -> None:
        coords = np.zeros((1, 4, 3))
        coords[0, 2, 1] = np.nan
        with pytest.raises(GeometryError, match="marked successful but contain non-finite"):
            IDRResult(coords, np.array([True]), "test", 1)

    def test_from_batch_blanks_failed_rows(self) -> None:
        coords = np.zeros((3, 4, 3))
        result = IDRResult.from_batch(coords, np.array([True, False, True]), "test", 2)
        assert np.all(np.isnan(result.ca_coords[1]))
        assert np.all(np.isfinite(result.ca_coords[0]))
        assert result.n_successful == 2
        assert not result.all_successful
        assert result.successful_coords.shape == (2, 4, 3)

    def test_derived_metrics_are_nan_for_failed_rows(self) -> None:
        coords = np.zeros((2, 4, 3))
        coords[:, :, 0] = np.arange(4)
        result = IDRResult.from_batch(coords, np.array([True, False]), "test", 1)
        assert np.isnan(result.end_to_end_distances[1])
        assert np.isnan(result.radii_of_gyration[1])
        assert result.end_to_end_distances[0] == pytest.approx(3.0)


class TestEngineProtocol:
    def test_walk_satisfies_the_protocol(self) -> None:
        assert isinstance(SelfAvoidingWalk(), ConformationEngine)

    def test_walk_is_always_available(self) -> None:
        # The dependency-light guarantee: numpy and scipy only.
        assert SelfAvoidingWalk().available() is True

    def test_name_is_recorded_on_the_result(self) -> None:
        engine = SelfAvoidingWalk()
        _, result = build(12, n_anchor=ORIGIN, engine=engine)
        assert result.engine == engine.name

    def test_rejects_a_bare_seed_instead_of_a_generator(self) -> None:
        request = IDRRequest("G" * 8, 8, 20.0, n_anchor_xyz=ORIGIN)
        with pytest.raises(TypeError, match=r"numpy\.random\.Generator"):
            SelfAvoidingWalk().generate(request, None, 7)  # type: ignore[arg-type]

    def test_engine_configuration_is_validated(self) -> None:
        with pytest.raises(EngineError, match="batch_size"):
            SelfAvoidingWalk(batch_size=0)
        with pytest.raises(EngineError, match="max_attempts"):
            SelfAvoidingWalk(max_attempts=0)
        with pytest.raises(EngineError, match="tolerance_fraction"):
            SelfAvoidingWalk(tolerance_fraction=0.0)


# ---------------------------------------------------------------------------
# Regression (a): the angle window
# ---------------------------------------------------------------------------

REGION_CASES = [
    pytest.param({"n_anchor": ORIGIN, "c_anchor": np.array([35.0, 0.0, 0.0])}, id="interior"),
    pytest.param({"n_anchor": ORIGIN}, id="c-terminal-tail"),
    pytest.param({"c_anchor": np.array([4.0, -7.0, 2.0])}, id="n-terminal-tail"),
    pytest.param({}, id="free"),
]


class TestAngleWindow:
    @pytest.mark.parametrize("case", REGION_CASES)
    def test_every_generated_angle_is_inside_the_window(self, case: dict) -> None:
        request, result = build(40, seed=3, **case)
        trace = full_trace(result.ca_coords[0], request.n_anchor_xyz, request.c_anchor_xyz)
        angles = pseudo_angles(trace)
        # Every vertex of this trace is a *generated* residue: an anchor sits at one end of
        # the trace and so is never a vertex. The pseudo-angle centred on an anchor is the
        # one angle the walk cannot constrain -- it needs the residue before the anchor,
        # which the request does not carry -- and it is not measured here because it does
        # not exist here.
        n_anchors = sum(x is not None for x in (request.n_anchor_xyz, request.c_anchor_xyz))
        assert angles.size == request.n_residues + n_anchors - 2
        assert angles.min() >= BACKBONE_ANGLE_MIN - 1e-6
        assert angles.max() <= BACKBONE_ANGLE_MAX + 1e-6

    def test_no_angle_below_the_window_across_many_conformers(self) -> None:
        # The old builder's uniform-sphere sampling with a 3.0 A clash cutoff permitted
        # 46.5 degrees and measured 47 in practice. Nothing here comes close: the window's
        # floor is 91.
        _, result = build(30, seed=11, n_anchor=ORIGIN, n_conformations=8)
        angles = np.concatenate(
            [pseudo_angles(full_trace(c, ORIGIN, None)) for c in result.ca_coords]
        )
        assert angles.min() > OLD_BUILDER_WORST_ANGLE
        assert angles.min() >= BACKBONE_ANGLE_MIN - 1e-6
        assert angles.max() <= BACKBONE_ANGLE_MAX + 1e-6

    def test_angle_distribution_is_not_uniform_over_the_window(self) -> None:
        # Candidates are weighted by the measured angle density (mean BACKBONE_ANGLE_MEAN, sd
        # BACKBONE_ANGLE_SD), so the generated angles must not be spread flat across the
        # window. Both bounds below are derived from the constants, not written down.
        _, result = build(60, seed=5, n_anchor=ORIGIN, n_conformations=6)
        angles = np.concatenate(
            [pseudo_angles(full_trace(c, ORIGIN, None)) for c in result.ca_coords]
        )
        # A uniform draw over the window has sd (high - low) / sqrt(12): 20.2 at the tuned
        # 91-161 window, 17.0 when the window was capped at 150. Measured here: 18.7.
        uniform_sd = (BACKBONE_ANGLE_MAX - BACKBONE_ANGLE_MIN) / np.sqrt(12.0)
        assert angles.std() < uniform_sd
        # The other half of the claim, which the sd alone does not make. It used to be asserted
        # as "the sample mean is closer to BACKBONE_ANGLE_MEAN than to the window's midpoint",
        # which only distinguishes anything while the window is lopsided: at 91-150 the
        # midpoint was 120.5 against a mean of 125, but at the tuned 91-161 the midpoint is
        # 126.0 and the two hypotheses are a degree apart -- indistinguishable, and the
        # measured mean (129.8) is above both, because the walk also steers for extension
        # toward its end-to-end target.
        #
        # Asserted instead where the density weighting and every other force in the engine
        # point the same way: the SHARP end of the window. A sharp turn is suppressed both by
        # the Gaussian weight and by self-avoidance (it folds residue i+1 back onto the chain),
        # so its share of the sample must fall short of the width share a flat draw would give
        # it. Measured over eight seeds: 0.04-0.07 of the sample against a flat 0.114.
        sharp = float(BACKBONE_ANGLE_MEAN - BACKBONE_ANGLE_SD)
        flat_share = (sharp - BACKBONE_ANGLE_MIN) / (BACKBONE_ANGLE_MAX - BACKBONE_ANGLE_MIN)
        assert float(np.mean(angles <= sharp)) < flat_share
        # And the window really is populated out to both ends, so the line above is describing
        # a shaped distribution rather than a truncated one.
        assert angles.min() < sharp < angles.max()


# ---------------------------------------------------------------------------
# Regression (c): the junctions
# ---------------------------------------------------------------------------


class TestJunctions:
    def test_first_residue_sits_one_bond_from_the_n_anchor(self) -> None:
        anchor = np.array([3.0, -4.0, 12.0])
        _, result = build(25, seed=1, n_anchor=anchor, c_anchor=anchor + np.array([20.0, 0, 0]))
        offset = float(np.linalg.norm(result.ca_coords[0, 0] - anchor))
        assert offset == pytest.approx(CA_CA_BOND_LENGTH, abs=1e-9)
        assert offset > 0.0

    def test_last_residue_sits_one_bond_from_the_c_anchor(self) -> None:
        c_anchor = np.array([20.0, 0.0, 0.0])
        _, result = build(25, seed=1, n_anchor=ORIGIN, c_anchor=c_anchor)
        offset = float(np.linalg.norm(result.ca_coords[0, -1] - c_anchor))
        # Solved analytically as a two-sphere intersection, so this is exact rather than
        # merely inside the bond tolerance.
        assert offset == pytest.approx(CA_CA_BOND_LENGTH, abs=1e-9)

    def test_n_terminal_tail_attaches_at_its_c_end(self) -> None:
        c_anchor = np.array([-6.0, 2.0, 1.0])
        _, result = build(20, seed=2, c_anchor=c_anchor)
        coords = result.ca_coords[0]
        # Growth runs outward from the anchor, but the output must be in sequence order:
        # the *last* residue is the one bonded to the C-terminal anchor.
        assert float(np.linalg.norm(coords[-1] - c_anchor)) == pytest.approx(
            CA_CA_BOND_LENGTH, abs=1e-9
        )
        assert float(np.linalg.norm(coords[0] - c_anchor)) > 2 * CA_CA_BOND_LENGTH

    @pytest.mark.parametrize("case", REGION_CASES)
    def test_nothing_is_coincident_anywhere_including_junctions(self, case: dict) -> None:
        request, result = build(35, seed=7, **case)
        trace = full_trace(result.ca_coords[0], request.n_anchor_xyz, request.c_anchor_xyz)
        assert closest_non_bonded(trace) >= CA_CLASH_DISTANCE
        # And no two atoms at all, bonded or not, are on top of each other.
        distances = np.linalg.norm(trace[:, None, :] - trace[None, :, :], axis=2)
        np.fill_diagonal(distances, np.inf)
        assert distances.min() > 1.0


class TestBondGeometry:
    @pytest.mark.parametrize("case", REGION_CASES)
    def test_every_bond_is_within_tolerance(self, case: dict) -> None:
        request, result = build(40, seed=4, **case)
        trace = full_trace(result.ca_coords[0], request.n_anchor_xyz, request.c_anchor_xyz)
        deviation = np.abs(bond_lengths(trace) - CA_CA_BOND_LENGTH)
        assert deviation.max() <= CA_CA_BOND_TOLERANCE

    def test_bonds_are_exact_to_machine_precision(self) -> None:
        # Nothing in the walk consumes the bond tolerance: every candidate is generated at
        # exactly one bond length and the closure is solved analytically. The tolerance
        # exists for other engines, not for this one.
        request, result = build(30, seed=8, n_anchor=ORIGIN, c_anchor=np.array([28.0, 0, 0]))
        trace = full_trace(result.ca_coords[0], request.n_anchor_xyz, request.c_anchor_xyz)
        assert np.abs(bond_lengths(trace) - CA_CA_BOND_LENGTH).max() < 1e-9


# ---------------------------------------------------------------------------
# Regression (b): the closure schedule follows the anchors, not the residue count
# ---------------------------------------------------------------------------


class TestClosureScheduleFollowsTheAnchors:
    def test_extension_fraction_is_seeded_from_the_anchor_separation(self) -> None:
        compact = _WalkPlan.build(
            IDRRequest(
                "G" * 40, 40, 40.0, n_anchor_xyz=ORIGIN, c_anchor_xyz=np.array([20.0, 0, 0])
            ),
            0.1,
            np.random.default_rng(0),
        )
        taut = _WalkPlan.build(
            IDRRequest(
                "G" * 40, 40, 40.0, n_anchor_xyz=ORIGIN, c_anchor_xyz=np.array([140.0, 0, 0])
            ),
            0.1,
            np.random.default_rng(0),
        )
        # Same residue count, same target, seven-fold different anchor separation.
        assert taut.extension_fraction > 6 * compact.extension_fraction
        goal_compact = compact.goal_for(10)
        goal_taut = taut.goal_for(10)
        assert goal_compact is not None and goal_taut is not None
        assert goal_taut.want.item() > goal_compact.want.item()

    def test_schedule_starts_at_the_anchor_separation(self) -> None:
        separation = 60.0
        plan = _WalkPlan.build(
            IDRRequest(
                "G" * 30, 30, 40.0, n_anchor_xyz=ORIGIN, c_anchor_xyz=np.array([separation, 0, 0])
            ),
            0.1,
            np.random.default_rng(0),
        )
        goal = plan.goal_for(0)
        assert goal is not None
        # One bond in from the N-anchor, the schedule still wants very nearly the full
        # anchor separation: it is not a function of chain length.
        assert goal.want.item() == pytest.approx(
            separation * max_reach(30) / max_reach(31), rel=1e-9
        )

    @pytest.mark.parametrize("separation", [20.0, 150.0])
    def test_every_residue_stays_inside_the_funnel(self, separation: float) -> None:
        # The old ramp, asked to bridge 20 A, looped out to 32 A before coming back,
        # because its schedule came from the residue count and nothing enforced it. Here
        # the funnel is a hard constraint at every single step, and it is a function of the
        # anchor separation.
        c_anchor = np.array([separation, 0.0, 0.0])
        request, result = build(40, seed=13, n_anchor=ORIGIN, c_anchor=c_anchor)
        plan = _WalkPlan.build(request, 0.1, np.random.default_rng(0))
        distances = np.linalg.norm(result.ca_coords[0] - c_anchor, axis=1)
        for index, distance in enumerate(distances):
            goal = plan.goal_for(index)
            assert goal is not None
            assert goal.lo.item() - 1e-9 <= distance <= goal.hi.item() + 1e-9

    def test_closure_ignores_a_target_that_disagrees_with_the_anchors(self) -> None:
        # Asked to bridge 20 A while carrying a 200 A target, the anchors win and the
        # region still closes exactly. The old builder took the residue-count-derived
        # dimension as gospel and missed its anchor.
        c_anchor = np.array([20.0, 0.0, 0.0])
        _, result = build(40, seed=13, n_anchor=ORIGIN, c_anchor=c_anchor, target=200.0)
        assert result.all_successful
        assert float(np.linalg.norm(result.ca_coords[0, -1] - c_anchor)) == pytest.approx(
            CA_CA_BOND_LENGTH, abs=1e-9
        )

    def test_a_taut_bridge_is_extended_and_a_slack_one_is_not(self) -> None:
        compact = build(50, seed=14, n_anchor=ORIGIN, c_anchor=np.array([20.0, 0, 0]))[1]
        taut = build(50, seed=14, n_anchor=ORIGIN, c_anchor=np.array([150.0, 0, 0]))[1]
        assert taut.radii_of_gyration[0] > 2.0 * compact.radii_of_gyration[0]

    @pytest.mark.parametrize("separation", [10.0, 60.0, 150.0, 185.0])
    def test_bridges_close_across_the_whole_feasible_range(self, separation: float) -> None:
        c_anchor = np.array([separation, 0.0, 0.0])
        request, result = build(50, seed=15, n_anchor=ORIGIN, c_anchor=c_anchor)
        assert result.all_successful
        trace = full_trace(result.ca_coords[0], request.n_anchor_xyz, request.c_anchor_xyz)
        assert np.abs(bond_lengths(trace) - CA_CA_BOND_LENGTH).max() < 1e-9
        angles = pseudo_angles(trace)
        assert angles.min() >= BACKBONE_ANGLE_MIN - 1e-6
        assert angles.max() <= BACKBONE_ANGLE_MAX + 1e-6


# ---------------------------------------------------------------------------
# Reachability
# ---------------------------------------------------------------------------


class TestReachability:
    def test_one_bond_reaches_exactly_one_bond_length(self) -> None:
        assert max_reach(1) == pytest.approx(CA_CA_BOND_LENGTH)
        assert min_reach(1) == pytest.approx(CA_CA_BOND_LENGTH)

    def test_zero_bonds_reach_nowhere(self) -> None:
        assert max_reach(0) == 0.0
        assert min_reach(0) == 0.0

    def test_two_bond_bounds_come_from_the_angle_window(self) -> None:
        assert max_reach(2) == pytest.approx(
            2 * CA_CA_BOND_LENGTH * np.sin(np.radians(BACKBONE_ANGLE_MAX) / 2)
        )
        assert min_reach(2) == pytest.approx(
            2 * CA_CA_BOND_LENGTH * np.sin(np.radians(BACKBONE_ANGLE_MIN) / 2)
        )

    def test_reach_is_monotone_and_below_the_contour_length(self) -> None:
        values = np.array([max_reach(k) for k in range(1, 30)])
        assert np.all(np.diff(values) > 0)
        contour = np.arange(1, 30) * CA_CA_BOND_LENGTH
        assert np.all(values <= contour + 1e-12)

    def test_negative_bond_counts_are_rejected(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            max_reach(-1)
        with pytest.raises(ValueError, match="non-negative"):
            min_reach(-1)

    def test_tail_is_ordered_longest_first(self) -> None:
        tail = reachability_tail(6)
        assert tail.shape == (6,)
        assert np.all(np.diff(tail) < 0)
        assert tail[-1] == pytest.approx(CA_CA_BOND_LENGTH)

    def test_derived_bound_is_looser_than_the_hardcoded_table_it_replaces(self) -> None:
        # The old literals were an empirical envelope, not a reachability bound: about
        # 0.69x the geometric maximum at six steps. A funnel built on them vetoes closures
        # that are perfectly achievable.
        derived = reachability_tail(6, OLD_TABLE_BOND_LENGTH)
        old = np.array(OLD_REACHABILITY_TABLE)
        assert np.all(derived >= old - 1e-9)
        assert derived[0] / old[0] > 1.4
        # Both agree on the one step where geometry leaves no freedom.
        assert derived[-1] == pytest.approx(old[-1])


# ---------------------------------------------------------------------------
# Skipping clash work that cannot change the answer
# ---------------------------------------------------------------------------


class TestReachableEnvelope:
    def test_no_anchor_leaves_the_region_unconstrained(self) -> None:
        assert reachable_envelope(20, None, None) is None

    def test_one_anchor_gives_the_walk_its_full_reach(self) -> None:
        assert reachable_envelope(20, ORIGIN, None) == pytest.approx(max_reach(20))
        assert reachable_envelope(20, None, ORIGIN) == pytest.approx(max_reach(20))

    def test_two_anchors_bound_the_sum_of_the_two_distances(self) -> None:
        # n + 1 bonds are shared between the two halves of the walk, so the sum of the two
        # anchor distances cannot exceed what n + 1 bonds span, plus the transverse slack
        # max_reach carries for odd bond counts.
        limit = reachable_envelope(20, ORIGIN, np.array([30.0, 0.0, 0.0]))
        assert limit is not None
        assert limit >= max_reach(21)
        assert limit < max_reach(21) + 2.0 * CA_CA_BOND_LENGTH

    def test_envelope_actually_contains_every_residue_it_builds(self) -> None:
        # The point of the bound. If a built region ever left the envelope, an obstacle
        # pruned by it could have been hit.
        n_anchor = ORIGIN
        c_anchor = np.array([28.0, 0.0, 0.0])
        limit = reachable_envelope(24, n_anchor, c_anchor)
        assert limit is not None
        for seed in range(6):
            _, result = build(24, seed=seed, n_anchor=n_anchor, c_anchor=c_anchor)
            trace = result.ca_coords[0]
            reach = np.linalg.norm(trace - n_anchor, axis=1) + np.linalg.norm(
                trace - c_anchor, axis=1
            )
            assert np.max(reach) <= limit + 1e-9

    def test_tail_envelope_contains_every_residue_it_builds(self) -> None:
        limit = reachable_envelope(40, ORIGIN, None)
        assert limit is not None
        for seed in range(6):
            _, result = build(40, seed=seed, n_anchor=ORIGIN)
            assert np.max(np.linalg.norm(result.ca_coords[0], axis=1)) <= limit + 1e-9

    def test_rejects_an_empty_region(self) -> None:
        with pytest.raises(ValueError, match="at least 1"):
            reachable_envelope(0, ORIGIN, None)


class TestObstaclePruning:
    @staticmethod
    def _request(n_residues: int, **anchors: np.ndarray | None) -> IDRRequest:
        return IDRRequest(
            sequence="G" * n_residues,
            n_residues=n_residues,
            target_end_to_end=flory_end_to_end(n_residues),
            **anchors,
        )

    def test_none_stays_none(self) -> None:
        request = self._request(10, n_anchor_xyz=ORIGIN)
        assert SelfAvoidingWalk._prune_unreachable(request, None) is None

    def test_an_obstacle_beyond_reach_is_dropped(self) -> None:
        request = self._request(4, n_anchor_xyz=ORIGIN)
        far = np.array([[1000.0, 0.0, 0.0]])
        pruned = SelfAvoidingWalk._prune_unreachable(request, far)
        assert pruned is not None
        assert pruned.shape[0] == 0

    def test_an_obstacle_just_inside_the_dilated_envelope_is_kept(self) -> None:
        request = self._request(4, n_anchor_xyz=ORIGIN)
        limit = reachable_envelope(4, ORIGIN, None)
        assert limit is not None
        edge = np.array([[limit + CA_CLASH_DISTANCE - 1e-6, 0.0, 0.0]])
        pruned = SelfAvoidingWalk._prune_unreachable(request, edge)
        assert pruned is not None
        assert pruned.shape[0] == 1

    def test_nothing_is_dropped_for_a_region_with_no_anchors(self) -> None:
        request = self._request(10)
        far = np.array([[500.0, 0.0, 0.0]])
        assert SelfAvoidingWalk._prune_unreachable(request, far) is far

    def test_a_long_region_keeps_everything_within_a_structure(self) -> None:
        # 150 residues reach far further than any AlphaFold model is wide, which is why this
        # pruning earns close to nothing end to end: it never fires on the regions that cost.
        request = self._request(150, n_anchor_xyz=ORIGIN)
        rng = np.random.default_rng(0)
        obstacles = rng.uniform(-100.0, 100.0, size=(500, 3))
        assert SelfAvoidingWalk._prune_unreachable(request, obstacles) is obstacles

    def test_a_non_finite_obstacle_still_reaches_the_guard_that_rejects_it(self) -> None:
        # NaN makes every comparison False, so a careless implementation drops the row and
        # silently disarms _obstacle_tree's non-finite check.
        request = self._request(12, n_anchor_xyz=ORIGIN)
        obstacles = np.array([[0.0, 0.0, 8.0], [np.nan, 0.0, 0.0]])
        with pytest.raises(GeometryError, match="non-finite"):
            SelfAvoidingWalk().generate(request, obstacles, np.random.default_rng(0))

    def test_pruning_does_not_change_what_gets_built(self) -> None:
        # An obstacle field spanning far beyond reach must give the identical trace to the
        # same field with the unreachable part removed by hand.
        rng = np.random.default_rng(3)
        obstacles = rng.uniform(-400.0, 400.0, size=(4000, 3))
        near = obstacles[np.linalg.norm(obstacles, axis=1) < 120.0]
        assert near.shape[0] < obstacles.shape[0]
        for seed in range(4):
            _, wide = build(20, seed=seed, n_anchor=ORIGIN, obstacles=obstacles)
            _, tight = build(20, seed=seed, n_anchor=ORIGIN, obstacles=near)
            assert np.array_equal(wide.ca_coords, tight.ca_coords)


class TestCandidateCullIsExact:
    """The cull replaces a 710-point clash query with a 1-point one when it cannot matter.

    It is only sound because every candidate for a residue sits *exactly* one bond length from a
    single centre. These tests pin that invariant and the identity of the output, since a wrong
    centre would silently let clashing candidates through rather than fail loudly.
    """

    def test_every_candidate_sits_one_bond_from_the_reported_centre(self) -> None:
        walk = SelfAvoidingWalk()
        plan = _WalkPlan.build(
            IDRRequest(
                sequence="G" * 12,
                n_residues=12,
                target_end_to_end=flory_end_to_end(12),
                n_anchor_xyz=ORIGIN,
                c_anchor_xyz=np.array([20.0, 0.0, 0.0]),
            ),
            0.2,
            np.random.default_rng(0),
        )
        rng = np.random.default_rng(1)
        coords = np.full((1, 12, 3), np.nan)
        live = np.array([0])
        for index in range(plan.n_grown):
            candidates, _, centres = walk._candidates_for(plan, coords, live, index, rng)
            distances = np.linalg.norm(candidates[0] - centres[0], axis=1)
            assert np.allclose(distances, CA_CA_BOND_LENGTH, atol=1e-9)
            coords[0, index] = candidates[0, 0]

    def test_a_clear_row_and_a_queried_row_agree_on_the_clash_predicate(self) -> None:
        rng = np.random.default_rng(0)
        obstacles = rng.uniform(-30.0, 30.0, size=(800, 3))
        tree = cKDTree(obstacles)
        centres = rng.uniform(-30.0, 30.0, size=(40, 3))
        directions = rng.normal(size=(40, 200, 3))
        directions /= np.linalg.norm(directions, axis=2, keepdims=True)
        candidates = centres[:, None, :] + CA_CA_BOND_LENGTH * directions

        culled = SelfAvoidingWalk._nearest_obstacle_distance(candidates, tree, centres)
        exact = SelfAvoidingWalk._nearest_obstacle_distance(candidates, tree, None)
        # Distances may differ (a culled row reports inf), the decision may not.
        assert np.array_equal(culled >= CA_CLASH_DISTANCE, exact >= CA_CLASH_DISTANCE)

    def test_a_culled_row_is_only_ever_one_that_could_not_clash(self) -> None:
        rng = np.random.default_rng(5)
        obstacles = rng.uniform(-30.0, 30.0, size=(600, 3))
        tree = cKDTree(obstacles)
        centres = rng.uniform(-30.0, 30.0, size=(60, 3))
        directions = rng.normal(size=(60, 150, 3))
        directions /= np.linalg.norm(directions, axis=2, keepdims=True)
        candidates = centres[:, None, :] + CA_CA_BOND_LENGTH * directions

        culled = SelfAvoidingWalk._nearest_obstacle_distance(candidates, tree, centres)
        exact = SelfAvoidingWalk._nearest_obstacle_distance(candidates, tree, None)
        skipped = np.all(np.isinf(culled), axis=1) & ~np.all(np.isinf(exact), axis=1)
        assert skipped.any(), "the fixture should exercise the cull at least once"
        assert np.all(exact[skipped] >= CA_CLASH_DISTANCE)

    def test_a_non_finite_centre_falls_back_to_querying_the_row(self) -> None:
        rng = np.random.default_rng(7)
        obstacles = rng.uniform(-10.0, 10.0, size=(200, 3))
        tree = cKDTree(obstacles)
        candidates = rng.uniform(-10.0, 10.0, size=(2, 30, 3))
        centres = np.array([[np.nan, np.nan, np.nan], [0.0, 0.0, 0.0]])
        culled = SelfAvoidingWalk._nearest_obstacle_distance(candidates, tree, centres)
        exact = SelfAvoidingWalk._nearest_obstacle_distance(candidates, tree, None)
        assert np.allclose(culled[0], exact[0])

    @pytest.mark.parametrize("n_residues", [12, 40])
    def test_output_is_unchanged_by_how_far_the_obstacles_sit(self, n_residues: int) -> None:
        # Obstacles hugging the region exercise the queried path; the same obstacles pushed
        # away exercise the culled path. Neither may differ from building with none at all.
        for seed in range(3):
            _, bare = build(n_residues, seed=seed, n_anchor=ORIGIN)
            _, far = build(
                n_residues,
                seed=seed,
                n_anchor=ORIGIN,
                obstacles=np.array([[0.0, 0.0, 900.0]]),
            )
            assert np.array_equal(bare.ca_coords, far.ca_coords)


# ---------------------------------------------------------------------------
# Hitting the target dimension
# ---------------------------------------------------------------------------


class TestTargetDimension:
    @pytest.mark.parametrize("n_residues", [10, 25, 60, 120])
    def test_tail_lands_within_tolerance_of_the_target(self, n_residues: int) -> None:
        request, result = build(n_residues, seed=n_residues, n_anchor=ORIGIN)
        achieved = float(result.end_to_end_distances[0])
        tolerance = max(0.1 * request.target_end_to_end, CA_CA_BOND_LENGTH)
        assert abs(achieved - request.target_end_to_end) <= tolerance

    def test_free_region_lands_within_tolerance_of_the_target(self) -> None:
        request, result = build(50, seed=21)
        tolerance = max(0.1 * request.target_end_to_end, CA_CA_BOND_LENGTH)
        assert abs(float(result.end_to_end_distances[0]) - request.target_end_to_end) <= tolerance

    def test_free_region_starts_at_the_origin(self) -> None:
        # Documented behaviour: with no anchor there is no frame, so the region is built at
        # the origin for the orchestrator to place.
        _, result = build(20, seed=22)
        assert np.allclose(result.ca_coords[0, 0], 0.0)

    def test_radius_of_gyration_tracks_the_ideal_chain_relation(self) -> None:
        # Re / Rg is sqrt(6) for an ideal chain and slightly larger for a self-avoiding
        # one, so Rg should sit just below sqrt(mean(Re ** 2) / 6) -- not far below it,
        # which would mean a collapsed blob, and not above, which would mean a rod.
        _, result = build(80, seed=23, n_anchor=ORIGIN, n_conformations=8)
        ideal = float(np.sqrt(np.mean(result.end_to_end_distances**2) / 6.0))
        measured = float(np.mean(result.radii_of_gyration))
        assert 0.80 * ideal <= measured <= 1.05 * ideal

    def test_interior_region_reports_what_the_anchors_dictate(self) -> None:
        # An interior region's end-to-end distance is set by its anchors, not by the
        # target. The result reports the achieved value so a caller can say so.
        separation = 100.0
        _, result = build(
            40, seed=24, n_anchor=ORIGIN, c_anchor=np.array([separation, 0, 0]), target=30.0
        )
        achieved = float(result.end_to_end_distances[0])
        assert abs(achieved - separation) < 3 * CA_CA_BOND_LENGTH


# ---------------------------------------------------------------------------
# Multiple conformations
# ---------------------------------------------------------------------------


class TestConformerBatches:
    def test_many_conformations(self) -> None:
        # The old code raised ValueError at a reshape for any value above 1, which killed
        # the multi-model pseudo-trajectory feature outright.
        _, result = build(20, seed=31, n_anchor=ORIGIN, n_conformations=9)
        assert result.ca_coords.shape == (9, 20, 3)
        assert result.all_successful
        assert result.success.dtype == np.bool_

    def test_conformers_are_distinct(self) -> None:
        _, result = build(20, seed=32, n_anchor=ORIGIN, n_conformations=5)
        for first in range(5):
            for second in range(first + 1, 5):
                assert not np.allclose(result.ca_coords[first], result.ca_coords[second])

    def test_batching_does_not_change_the_contract(self) -> None:
        # Exercise the chunking path: seven conformers in batches of two.
        _, result = build(
            18,
            seed=33,
            n_anchor=ORIGIN,
            c_anchor=np.array([25.0, 0, 0]),
            n_conformations=7,
            engine=SelfAvoidingWalk(batch_size=2),
        )
        assert result.n_successful == 7
        for coords in result.ca_coords:
            trace = full_trace(coords, ORIGIN, np.array([25.0, 0, 0]))
            assert np.abs(bond_lengths(trace) - CA_CA_BOND_LENGTH).max() < 1e-9

    def test_repeated_seeds_reproduce_exactly(self) -> None:
        first = build(25, seed=34, n_anchor=ORIGIN, n_conformations=3)[1]
        second = build(25, seed=34, n_anchor=ORIGIN, n_conformations=3)[1]
        assert np.array_equal(first.ca_coords, second.ca_coords)

    def test_different_seeds_differ(self) -> None:
        first = build(25, seed=35, n_anchor=ORIGIN)[1]
        second = build(25, seed=36, n_anchor=ORIGIN)[1]
        assert not np.allclose(first.ca_coords, second.ca_coords)


# ---------------------------------------------------------------------------
# Obstacles
# ---------------------------------------------------------------------------


def sphere_of_obstacles(centre: np.ndarray, radius: float, spacing: float = 1.5) -> np.ndarray:
    """Return a dense ball of points, standing in for an already-placed folded domain."""
    axis = np.arange(-radius, radius + spacing / 2, spacing)
    grid = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    grid = grid[np.linalg.norm(grid, axis=1) <= radius]
    return grid + centre


class TestObstacles:
    def test_walk_routes_around_an_obstacle(self) -> None:
        obstacles = sphere_of_obstacles(np.array([20.0, 0.0, 0.0]), 8.0)
        _, result = build(
            45, seed=41, n_anchor=ORIGIN, c_anchor=np.array([40.0, 0, 0]), obstacles=obstacles
        )
        closest = float(
            np.min(np.linalg.norm(result.ca_coords[0][:, None, :] - obstacles[None, :, :], axis=2))
        )
        assert closest >= CA_CLASH_DISTANCE
        assert result.relaxed_to is None

    def test_relaxation_is_reported_not_hidden(self) -> None:
        # A gap too narrow for the strict threshold: the walk may squeeze through, but only
        # if it says so. Either outcome is acceptable; a silent squeeze is not.
        wall = np.stack(
            np.meshgrid(
                np.array([18.0]), np.arange(-30, 30, 2.6), np.arange(-30, 30, 2.6), indexing="ij"
            ),
            axis=-1,
        ).reshape(-1, 3)
        try:
            _, result = build(
                40, seed=42, n_anchor=ORIGIN, c_anchor=np.array([36.0, 0, 0]), obstacles=wall
            )
        except ExhaustedAttemptsError:
            return
        if result.relaxed_to is None:
            closest = float(
                np.min(np.linalg.norm(result.ca_coords[0][:, None, :] - wall[None, :, :], axis=2))
            )
            assert closest >= CA_CLASH_DISTANCE
        else:
            assert result.relaxed_to in CLASH_RELAXATION_LADDER
            assert result.relaxed_to < CA_CLASH_DISTANCE

    def test_non_finite_obstacles_are_rejected(self) -> None:
        with pytest.raises(GeometryError, match="non-finite"):
            build(10, n_anchor=ORIGIN, obstacles=np.array([[np.nan, 0.0, 0.0]]))

    def test_wrongly_shaped_obstacles_are_rejected(self) -> None:
        with pytest.raises(GeometryError, match=r"shape \(n, 3\)"):
            build(10, n_anchor=ORIGIN, obstacles=np.zeros((4, 2)))

    def test_empty_obstacle_array_is_the_same_as_none(self) -> None:
        with_empty = build(20, seed=43, n_anchor=ORIGIN, obstacles=np.zeros((0, 3)))[1]
        with_none = build(20, seed=43, n_anchor=ORIGIN)[1]
        assert np.array_equal(with_empty.ca_coords, with_none.ca_coords)

    def test_free_region_refuses_an_occupied_origin(self) -> None:
        with pytest.raises(GeometryError, match="own frame starting at the origin"):
            build(10, obstacles=np.zeros((1, 3)))


# ---------------------------------------------------------------------------
# Failure is loud
# ---------------------------------------------------------------------------


class TestFailureIsLoud:
    def test_anchors_further_apart_than_the_chain_can_span(self) -> None:
        with pytest.raises(GeometryError, match="Cannot bridge anchors"):
            build(5, n_anchor=ORIGIN, c_anchor=np.array([100.0, 0, 0]))

    def test_anchors_too_close_for_the_angle_window(self) -> None:
        with pytest.raises(GeometryError, match="cannot double back"):
            build(1, n_anchor=ORIGIN, c_anchor=np.array([3.0, 0, 0]), target=3.8)

    def test_target_beyond_the_contour_length(self) -> None:
        with pytest.raises(UnsatisfiableTargetError) as info:
            build(10, n_anchor=ORIGIN, target=500.0)
        assert info.value.target == pytest.approx(500.0)
        assert info.value.achievable is not None

    def test_target_below_what_the_angle_window_allows(self) -> None:
        with pytest.raises(UnsatisfiableTargetError, match="cannot be that compact"):
            build(3, n_anchor=ORIGIN, target=0.5)

    def test_a_region_with_no_room_raises_rather_than_returning_zeros(self) -> None:
        # Bury the anchor in a solid block of obstacles finer than the loosest rung of the
        # relaxation ladder, so no candidate can be placed at any threshold.
        block = sphere_of_obstacles(ORIGIN, 9.0, spacing=1.4)
        block = block[np.linalg.norm(block - ORIGIN, axis=1) > 1.0]
        with pytest.raises(ExhaustedAttemptsError) as info:
            build(8, n_anchor=ORIGIN, obstacles=block)
        assert info.value.attempts == MAX_ATTEMPTS_PER_REGION

    def test_partial_batches_are_refused_by_default(self) -> None:
        # require_all is the default because a caller asking for ten models and silently
        # receiving seven is exactly the class of surprise this rewrite removes.
        assert SelfAvoidingWalk().require_all is True

    def test_no_result_ever_contains_the_origin_placeholder_rows(self) -> None:
        _, result = build(30, seed=51, n_anchor=ORIGIN, n_conformations=4)
        # The old failure signature: rows of exact (0, 0, 0). The N-anchor here *is* the
        # origin, so a generated residue landing on it would be indistinguishable from the
        # old placeholder -- and is impossible, because it would be a coincident atom.
        assert not np.any(np.all(result.ca_coords == 0.0, axis=2))
        assert np.all(np.isfinite(result.ca_coords))


# ---------------------------------------------------------------------------
# Small regions and edge cases
# ---------------------------------------------------------------------------


class TestSmallRegions:
    def test_single_residue_between_two_anchors(self) -> None:
        c_anchor = np.array([6.5, 0.0, 0.0])
        _, result = build(1, seed=61, n_anchor=ORIGIN, c_anchor=c_anchor, target=6.5)
        position = result.ca_coords[0, 0]
        assert float(np.linalg.norm(position - ORIGIN)) == pytest.approx(
            CA_CA_BOND_LENGTH, abs=1e-9
        )
        assert float(np.linalg.norm(position - c_anchor)) == pytest.approx(
            CA_CA_BOND_LENGTH, abs=1e-9
        )

    def test_single_residue_tail(self) -> None:
        _, result = build(1, seed=62, n_anchor=ORIGIN, target=CA_CA_BOND_LENGTH)
        assert result.all_successful
        assert np.isnan(result.end_to_end_distances[0])

    def test_two_residue_interior_region(self) -> None:
        c_anchor = np.array([9.0, 1.0, 0.0])
        request, result = build(2, seed=63, n_anchor=ORIGIN, c_anchor=c_anchor, target=4.0)
        trace = full_trace(result.ca_coords[0], request.n_anchor_xyz, request.c_anchor_xyz)
        assert np.abs(bond_lengths(trace) - CA_CA_BOND_LENGTH).max() < 1e-9
        angles = pseudo_angles(trace)
        assert angles.min() >= BACKBONE_ANGLE_MIN - 1e-6
        assert angles.max() <= BACKBONE_ANGLE_MAX + 1e-6

    def test_minimum_idr_length_builds(self) -> None:
        request, result = build(4, seed=64, n_anchor=ORIGIN, c_anchor=np.array([12.0, 0, 0]))
        assert result.all_successful
        trace = full_trace(result.ca_coords[0], request.n_anchor_xyz, request.c_anchor_xyz)
        assert np.abs(bond_lengths(trace) - CA_CA_BOND_LENGTH).max() < 1e-9


# ---------------------------------------------------------------------------
# Regression (d): the junction pseudo-angle, centred on the anchor itself
#
# Confirmed defect this section exists for: the CA-CA-CA angle centred on an anchor was
# neither constrained during generation nor measured afterwards. Measured over 480
# junctions of interior regions (n = 20/40/50, anchor separations 30-90 A, 20 seeds):
# range 38.62-176.76 deg, 41.2% below BACKBONE_ANGLE_MIN and 3.3% above
# BACKBONE_ANGLE_MAX, and IDRResult.success was True for every one. That is the same
# defect class (47-178 deg) the module docstring claims to fix by construction, merely
# relocated to the seams -- and validate_ca_trace on the spliced trace failed for 171 of
# those 480 conformers while the engine reported success.
# ---------------------------------------------------------------------------


def spliced_trace(
    coords: np.ndarray,
    *,
    n_prev: np.ndarray | None = None,
    n_anchor: np.ndarray | None = None,
    c_anchor: np.ndarray | None = None,
    c_next: np.ndarray | None = None,
) -> np.ndarray:
    """Return the trace a viewer draws once the region is written back into its chain."""
    pieces = [np.asarray(p, dtype=np.float64)[None, :] for p in (n_prev, n_anchor) if p is not None]
    pieces.append(coords)
    pieces.extend(
        np.asarray(p, dtype=np.float64)[None, :] for p in (c_anchor, c_next) if p is not None
    )
    return np.concatenate(pieces, axis=0)


def flanked_interior(
    n_residues: int,
    separation: float,
    *,
    seed: int,
    n_conformations: int = 4,
    with_outer: bool = True,
) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Build an interior region inside a chain that has residues beyond both anchors.

    Returns the fixed context (the four flanking CA positions) and the conformers. The
    outer residues are placed one bond length from their anchor at a 125 degree
    pseudo-angle, i.e. exactly the geometry a real backbone has there.
    """
    n_anchor = np.zeros(3)
    c_anchor = np.array([separation, 0.0, 0.0])
    lean = CA_CA_BOND_LENGTH * np.array([-np.cos(np.radians(55.0)), np.sin(np.radians(55.0)), 0.0])
    context = {
        "n_anchor": n_anchor,
        "c_anchor": c_anchor,
        "n_prev": n_anchor + lean,
        "c_next": c_anchor - np.array([lean[0], -lean[1], 0.0]),
    }
    request = IDRRequest(
        sequence="G" * n_residues,
        n_residues=n_residues,
        target_end_to_end=flory_end_to_end(n_residues),
        n_anchor_xyz=n_anchor,
        c_anchor_xyz=c_anchor,
        n_anchor_prev_xyz=context["n_prev"] if with_outer else None,
        c_anchor_next_xyz=context["c_next"] if with_outer else None,
        n_conformations=n_conformations,
    )
    result = SelfAvoidingWalk().generate(request, None, np.random.default_rng(seed))
    assert result.all_successful
    return context, result.ca_coords


class TestJunctionAngles:
    @pytest.mark.parametrize(
        ("n_residues", "separation"), [(20, 30.0), (50, 30.0), (40, 90.0), (20, 60.0)]
    )
    def test_junction_angles_are_inside_the_window(
        self, n_residues: int, separation: float
    ) -> None:
        # The angle centred on the N-anchor is (residue before the anchor, anchor, first
        # generated residue); on the C-anchor it is (last generated, anchor, residue
        # after). Both are angles of the output structure, so both must be in the window.
        for seed in range(6):
            context, conformers = flanked_interior(n_residues, separation, seed=seed)
            for coords in conformers:
                angles = pseudo_angles(spliced_trace(coords, **context))
                junctions = np.array([angles[0], angles[-1]])
                assert junctions.min() >= BACKBONE_ANGLE_MIN - 1e-6
                assert junctions.max() <= BACKBONE_ANGLE_MAX + 1e-6

    def test_spliced_trace_passes_the_package_wide_gate(self) -> None:
        # The region is written into a real chain, not consumed in isolation, so the gate
        # that matters is the one applied to the assembled trace -- including its
        # non-bonded steric check.
        from dodo.geometry.metrics import validate_ca_trace

        for seed in range(10):
            context, conformers = flanked_interior(30, 40.0, seed=seed, n_conformations=2)
            for coords in conformers:
                report = validate_ca_trace(spliced_trace(coords, **context))
                assert report.ok, report.describe()

    def test_terminal_tail_junction_is_constrained_too(self) -> None:
        anchor = np.array([2.0, 1.0, -3.0])
        prev = anchor + np.array([-2.2, 3.1, 0.0])
        request = IDRRequest(
            "G" * 25,
            25,
            flory_end_to_end(25),
            n_anchor_xyz=anchor,
            n_anchor_prev_xyz=prev,
            n_conformations=5,
        )
        result = SelfAvoidingWalk().generate(request, None, np.random.default_rng(2))
        for coords in result.ca_coords:
            angles = pseudo_angles(spliced_trace(coords, n_prev=prev, n_anchor=anchor))
            assert angles.min() >= BACKBONE_ANGLE_MIN - 1e-6
            assert angles.max() <= BACKBONE_ANGLE_MAX + 1e-6

    def test_n_terminal_tail_junction_is_constrained_at_its_c_end(self) -> None:
        anchor = np.array([-4.0, 0.0, 1.0])
        following = anchor + np.array([3.0, 2.33, 0.0])
        request = IDRRequest(
            "G" * 18,
            18,
            flory_end_to_end(18),
            c_anchor_xyz=anchor,
            c_anchor_next_xyz=following,
            n_conformations=4,
        )
        result = SelfAvoidingWalk().generate(request, None, np.random.default_rng(4))
        for coords in result.ca_coords:
            angles = pseudo_angles(spliced_trace(coords, c_anchor=anchor, c_next=following))
            assert angles.min() >= BACKBONE_ANGLE_MIN - 1e-6
            assert angles.max() <= BACKBONE_ANGLE_MAX + 1e-6

    def test_single_residue_junction_between_two_flanked_anchors(self) -> None:
        # Both junction angles are centred on anchors and the one generated residue is
        # placed by the analytic closure step, so this is the case where every angle in
        # the trace depends on the closure filter.
        context, conformers = flanked_interior(1, 6.6, seed=8, n_conformations=3)
        for coords in conformers:
            angles = pseudo_angles(spliced_trace(coords, **context))
            assert angles.min() >= BACKBONE_ANGLE_MIN - 1e-6
            assert angles.max() <= BACKBONE_ANGLE_MAX + 1e-6

    def test_missing_flanking_residues_are_reported_not_ignored(self) -> None:
        from dodo.engines.walk import UnconstrainedJunctionWarning

        request = IDRRequest(
            "G" * 12,
            12,
            flory_end_to_end(12),
            n_anchor_xyz=ORIGIN,
            c_anchor_xyz=np.array([20.0, 0.0, 0.0]),
        )
        with pytest.warns(UnconstrainedJunctionWarning, match="junction pseudo-angle"):
            result = SelfAvoidingWalk().generate(request, None, np.random.default_rng(0))
        assert result.unconstrained_junctions == ("n_anchor", "c_anchor")

    def test_no_warning_when_both_flanking_residues_are_supplied(self) -> None:
        context, _ = flanked_interior(15, 25.0, seed=1, n_conformations=1)
        request = IDRRequest(
            "G" * 15,
            15,
            flory_end_to_end(15),
            n_anchor_xyz=context["n_anchor"],
            c_anchor_xyz=context["c_anchor"],
            n_anchor_prev_xyz=context["n_prev"],
            c_anchor_next_xyz=context["c_next"],
        )
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            result = SelfAvoidingWalk().generate(request, None, np.random.default_rng(1))
        assert result.unconstrained_junctions == ()

    def test_outer_residue_without_its_anchor_is_refused(self) -> None:
        with pytest.raises(EngineError, match="n_anchor_prev_xyz"):
            IDRRequest("GGG", 3, 10.0, n_anchor_prev_xyz=ORIGIN)


# ---------------------------------------------------------------------------
# Hitting -- or refusing -- a compact target
#
# Confirmed defect: compact targets were silently overshot by up to +658% with
# success=True (n=10 target 0.50 A -> Re 3.28-4.21 A; n=50 target 1.00 A -> 3.97-4.78 A),
# because the tolerance was floored at CA_CA_BOND_LENGTH and min_reach(k) returned 0.0
# for k >= 4, making both the pre-flight gate and the post-hoc audit unreachable.
# ---------------------------------------------------------------------------


class TestCompactTargets:
    @pytest.mark.parametrize(
        ("n_residues", "target"), [(10, 0.5), (50, 1.0), (20, 2.0), (200, 2.0), (10, 2.5)]
    )
    def test_a_target_below_the_steric_floor_raises(self, n_residues: int, target: float) -> None:
        # Two CA atoms three or more residues apart cannot be closer than
        # CA_CLASH_DISTANCE, so no residue count can produce these end-to-end distances.
        # The engine's own docstring promises UnsatisfiableTargetError here.
        with pytest.raises(UnsatisfiableTargetError):
            build(n_residues, target=target, n_anchor=ORIGIN, n_conformations=5)

    @pytest.mark.parametrize(
        ("n_residues", "target"), [(200, 5.0), (10, 3.5), (50, 6.0), (20, 5.0), (100, 10.0)]
    )
    def test_a_target_just_above_the_floor_is_hit_or_refused(
        self, n_residues: int, target: float
    ) -> None:
        # These are physically possible -- a long chain can bring its two ends together --
        # so the engine must either build them accurately or refuse. What it must not do is
        # what it used to: return 4 A for a 0.5 A request with success=True.
        try:
            _, result = build(n_residues, seed=5, target=target, n_anchor=ORIGIN, n_conformations=6)
        except UnsatisfiableTargetError:
            return
        achieved = result.end_to_end_distances[result.success]
        # Per conformer, against the ensemble mean the request named: an individual conformer
        # is steered to its own draw, so the ensemble is what carries the requested value.
        assert abs(float(achieved.mean()) - target) <= 0.25 * target
        assert float(achieved.min()) >= CA_CLASH_DISTANCE - 1e-9

    def test_min_reach_never_claims_two_atoms_can_coincide(self) -> None:
        for k in range(3, 40):
            assert min_reach(k) >= CA_CLASH_DISTANCE

    @pytest.mark.parametrize("n_residues", [8, 20, 60])
    def test_a_feasible_compact_target_is_actually_achieved(self, n_residues: int) -> None:
        # Half the Flory estimate: well above the steric floor, so it must be hit rather
        # than overshot. Previously n=10 at half the estimate came out +26% with no
        # warning, because the 3.8 A tolerance floor swallowed the miss.
        target = 0.5 * flory_end_to_end(n_residues)
        _, result = build(
            n_residues, seed=n_residues, n_anchor=ORIGIN, target=target, n_conformations=12
        )
        achieved = result.end_to_end_distances
        assert result.all_successful
        assert abs(float(achieved.mean()) - target) <= 0.1 * target

    def test_tolerance_is_not_floored_at_a_bond_length(self) -> None:
        plan = _WalkPlan.build(
            IDRRequest("G" * 12, 12, 12.0, n_anchor_xyz=ORIGIN), 0.1, np.random.default_rng(0)
        )
        assert plan.tolerance == pytest.approx(1.2)


# ---------------------------------------------------------------------------
# The ensemble: a predicted end-to-end distance is a mean, not a per-conformer pin
#
# Confirmed defect: CV(Re) across a requested batch measured 0.006-0.045 where a matched
# freely-rotating chain gives 0.35-0.48, so 60 models of a 200-residue IDR spanned 1.9 A
# of extension -- one conformation sampled sixty times. Separately, extended targets
# landed deterministically on the corridor floor: mean(Re)/target 0.962-0.966 for
# max_expansion at n = 50/100/200.
# ---------------------------------------------------------------------------

ENSEMBLE_CV_MIN = 0.20
ENSEMBLE_CV_MAX = 0.60


class TestEnsembleStatistics:
    @pytest.mark.parametrize("n_residues", [50, 100])
    @pytest.mark.parametrize("multiplier", [0.7, 1.0, 1.3])
    def test_ensemble_mean_matches_the_target(self, n_residues: int, multiplier: float) -> None:
        target = multiplier * flory_end_to_end(n_residues)
        _, result = build(n_residues, seed=7, n_anchor=ORIGIN, target=target, n_conformations=30)
        assert result.all_successful
        mean = float(result.end_to_end_distances.mean())
        assert abs(mean - target) <= 0.05 * target

    @pytest.mark.parametrize("n_residues", [50, 100, 200])
    def test_ensemble_is_as_wide_as_a_polymer_ensemble(self, n_residues: int) -> None:
        target = flory_end_to_end(n_residues)
        _, result = build(n_residues, seed=11, n_anchor=ORIGIN, target=target, n_conformations=30)
        achieved = result.end_to_end_distances
        cv = float(achieved.std() / achieved.mean())
        assert ENSEMBLE_CV_MIN <= cv <= ENSEMBLE_CV_MAX

    def test_extended_targets_are_not_pinned_to_the_corridor_floor(self) -> None:
        # mean(Re)/target was 0.9627-0.9656 for max_expansion, with every conformer within
        # 1.9 A of lo = target - 0.5 * tolerance. The bias was one-sided and systematic.
        for n_residues in (50, 100):
            target = 2.0 * flory_end_to_end(n_residues)
            _, result = build(
                n_residues, seed=13, n_anchor=ORIGIN, target=target, n_conformations=20
            )
            ratio = float(result.end_to_end_distances.mean()) / target
            assert 0.98 <= ratio <= 1.02

    def test_one_conformation_is_steered_straight_at_the_target(self) -> None:
        # A batch of one has no ensemble to spread, so its target is the requested value
        # exactly: the spread is around the mean, not a licence to miss.
        target = flory_end_to_end(40)
        plan = _WalkPlan.build(
            IDRRequest("G" * 40, 40, target, n_anchor_xyz=ORIGIN), 0.1, np.random.default_rng(0)
        )
        assert plan.targets.shape == (1,)
        assert float(plan.targets[0]) == pytest.approx(target)

    def test_per_conformer_targets_are_drawn_around_the_requested_mean(self) -> None:
        target = flory_end_to_end(100)
        plan = _WalkPlan.build(
            IDRRequest("G" * 100, 100, target, n_anchor_xyz=ORIGIN, n_conformations=200),
            0.1,
            np.random.default_rng(0),
        )
        assert float(plan.targets.mean()) == pytest.approx(target, rel=1e-6)
        cv = float(plan.targets.std() / plan.targets.mean())
        assert ENSEMBLE_CV_MIN <= cv <= ENSEMBLE_CV_MAX
        # Maxwell-Boltzmann-like: no mass at zero, a tail on the high side.
        assert float(plan.targets.min()) > 0.0
        assert float(np.median(plan.targets)) < target

    def test_targets_are_reproducible_from_the_seed(self) -> None:
        first = build(40, seed=99, n_anchor=ORIGIN, n_conformations=8)[1]
        second = build(40, seed=99, n_anchor=ORIGIN, n_conformations=8)[1]
        assert np.array_equal(first.ca_coords, second.ca_coords)


# ---------------------------------------------------------------------------
# Local shape
#
# Confirmed defect: with a flat TARGET_STEERING_WIDTH the radial preference and the
# pseudo-angle preference were sharp at the same time, and a bond direction satisfying both
# has only two solutions -- mirror images across a plane that turns slowly as the chain
# grows. Steered regions therefore came out as flat zig-zags: measured on p300's 583-residue
# terminal region, CA pseudo-dihedrals piled up at 0 and 180 degrees (planar order 0.618
# against 0.05-0.08 for a real coil) and the pseudo-angle mean was dragged to 131 degrees.
#
# Nothing in the suite caught it, and that is the point of these tests. Every existing check
# on a steered region measures a *global* observable -- the end-to-end distance, the ensemble
# CV, the scaling exponent -- and all of them stayed healthy, because each conformer picked a
# different plane. The defect lives in the local geometry, so it needs a local measurement.
# ---------------------------------------------------------------------------


def _pseudo_dihedrals(trace: np.ndarray) -> np.ndarray:
    """CA pseudo-dihedrals in degrees, ``(n - 3,)``."""
    b0 = trace[1:-2] - trace[:-3]
    b1 = trace[2:-1] - trace[1:-2]
    b2 = trace[3:] - trace[2:-1]
    n1 = np.cross(b0, b1)
    n2 = np.cross(b1, b2)
    m1 = np.cross(n1, b1 / np.linalg.norm(b1, axis=1, keepdims=True))
    return np.degrees(np.arctan2(np.sum(m1 * n2, axis=1), np.sum(n1 * n2, axis=1)))


def _planar_order(trace: np.ndarray) -> float:
    """``|<exp(2i * dihedral)>|``: 1.0 for a planar zig-zag, ~0 for free dihedrals.

    Doubling the angle before averaging is what makes this a *planarity* measure rather than
    a trans-content measure: 0 and 180 degrees are both in-plane and both map to the same
    point, so a chain alternating between them scores 1.0 just as a pure all-trans one does.
    That alternation is exactly the artifact, and an unsigned dihedral histogram misses it.
    """
    radians = np.radians(_pseudo_dihedrals(trace))
    return float(np.hypot(np.cos(2 * radians).mean(), np.sin(2 * radians).mean()))


#: Planar order of a freely-rotating chain carrying DODO's own angle distribution, measured
#: over 20 chains at each of n = 50, 120, 300 and 583: 0.127, 0.083, 0.047, 0.035. It falls
#: with length because the average runs over more dihedrals.
#:
#: Two ceilings rather than one, because a single conformer's planar order scatters: measured
#: at n = 120 over eight conformers the mean is 0.105 and the highest single value 0.258, so a
#: bound on the maximum tight enough to be meaningful would be flaky. The mean is the robust
#: statistic and carries the tight bound; the per-conformer bound only has to be tight enough
#: to catch the regression, whose *mean* was 0.46-0.62.
MAX_MEAN_PLANAR_ORDER = 0.20
MAX_ANY_PLANAR_ORDER = 0.40


class TestSteeredRegionsAreNotPlanar:
    @pytest.mark.parametrize("n_residues", [120, 300])
    def test_dihedrals_are_not_pinned_to_the_plane(self, n_residues: int) -> None:
        # p300's real terminal regions sit at 0.085-0.126 of the geometric reach, which is
        # where the artifact was worst; the predicted dimension lands in that band.
        target = flory_end_to_end(n_residues)
        _, result = build(n_residues, seed=5, n_anchor=ORIGIN, target=target, n_conformations=8)
        assert result.all_successful
        orders = np.array([_planar_order(trace) for trace in result.ca_coords])
        assert orders.mean() < MAX_MEAN_PLANAR_ORDER, (
            f"mean planar order {orders.mean():.3f} at n={n_residues}: the steered walk is "
            f"laying the chain out in a plane again. A real coil scores under 0.13 at this "
            f"length and the flat-width regression scored 0.46-0.60."
        )
        assert orders.max() < MAX_ANY_PLANAR_ORDER, (
            f"one conformer reached planar order {orders.max():.3f} at n={n_residues}, past "
            f"anything per-conformer scatter explains."
        )

    @pytest.mark.parametrize("n_residues", [120, 300])
    def test_dihedrals_avoid_the_perpendicular_hole(self, n_residues: int) -> None:
        # The sharper signature, and the one a summary statistic can hide: the two-solution
        # collapse does not merely favour 0 and 180, it *empties* +-90, because those are the
        # dihedrals that carry the chain out of the plane. Measured on the regression, the
        # +-60-to-120 band held 6% of dihedrals against 33% for a uniform distribution.
        target = flory_end_to_end(n_residues)
        _, result = build(n_residues, seed=6, n_anchor=ORIGIN, target=target, n_conformations=8)
        assert result.all_successful
        dihedrals = np.concatenate([_pseudo_dihedrals(t) for t in result.ca_coords])
        magnitude = np.abs(dihedrals)
        occupancy = float(((magnitude >= 60.0) & (magnitude <= 120.0)).mean())
        assert occupancy > 0.20, (
            f"only {occupancy:.1%} of pseudo-dihedrals are within 30 degrees of "
            f"perpendicular (uniform would be 33%); the out-of-plane directions are being "
            f"starved."
        )

    def test_the_angle_distribution_is_not_dragged_extended(self) -> None:
        # Same cause, separate symptom: rings with a larger cone radius sweep a wider range of
        # distances, so a sharp radial filter over-selects them and the realized pseudo-angle
        # mean rises. Measured at 131 deg on the regression against the 125 deg DODO samples.
        target = flory_end_to_end(200)
        _, result = build(200, seed=7, n_anchor=ORIGIN, target=target, n_conformations=8)
        assert result.all_successful
        angles = np.concatenate([np.degrees(_pseudo_angles(trace)) for trace in result.ca_coords])
        assert abs(float(angles.mean()) - BACKBONE_ANGLE_MEAN) < 3.0

    def test_closure_regions_were_never_affected(self) -> None:
        # The control. A region pinned between two anchors steers with _natural_fluctuation,
        # which is tens of Angstroms wide, so it never had the defect and must not acquire
        # one from any future change to the steered path.
        far = np.array([120.0, 0.0, 0.0])
        _, result = build(
            120, seed=8, n_anchor=ORIGIN, c_anchor=far, target=120.0, n_conformations=6
        )
        assert result.all_successful
        assert np.mean([_planar_order(trace) for trace in result.ca_coords]) < MAX_MEAN_PLANAR_ORDER


class TestCompactDrawsInAMultiModelEnsemble:
    """The width cap keys on the region, not on one conformer's draw.

    Confirmed defect, and it only appeared at ``n_models > 1``. The pipeline draws a
    per-model end-to-end target for every terminal region and then calls the engine once per
    model with ``n_conformations=1``, so from inside the engine a compact draw is
    indistinguishable from a genuinely small region. The cap -- a fraction of the target,
    there to stop short regions paying accuracy for a benefit they cannot get -- therefore
    handed compact draws a near-floor width and they came out flat while their siblings did
    not. Measured on p300's 583-residue tail at twenty models: the two most compact draws
    (27.9 and 72.1 A against a mean of 186) scored 0.474 and 0.132, the other eighteen
    0.011-0.082, and width against planar order correlated at -0.96.

    Isolating width from target settles which one is responsible. Holding the target at 27.9 A
    and forcing the width by hand, planar order runs 0.476 / 0.144 / 0.055 / 0.034 at widths
    0.56 / 1.20 / 2.23 / 4.00 while the achieved ratio stays at 0.9935-1.0448. It is the
    width, and widening it costs the compact conformer essentially nothing.
    """

    def test_the_cap_reads_the_region_mean_not_the_conformer_draw(self) -> None:
        span, mean = 582, 186.0
        compact = np.array([[27.9]])
        from_mean = float(_target_steering_width(compact, span, mean)[0, 0])
        from_draw = float(_target_steering_width(compact, span, 27.9)[0, 0])
        assert from_draw == pytest.approx(0.02 * 27.9, rel=1e-6), (
            "capping on the draw is what produced the flat conformer; this documents it"
        )
        assert from_mean > 2.0, (
            f"a compact draw from a 186 A ensemble got width {from_mean:.2f}; the region is "
            f"583 residues long and its width must not collapse because one draw was compact"
        )

    def test_a_compact_draw_from_a_long_ensemble_is_not_planar(self) -> None:
        n = 400
        mean = flory_end_to_end(n)
        # 0.15 * mean is exactly the low clip dodo.construct.pipeline applies to its own draws,
        # so this is the most compact conformer the multi-model path can actually ask for.
        request = IDRRequest(
            sequence="G" * n,
            n_residues=n,
            target_end_to_end=0.15 * mean,
            n_anchor_xyz=ORIGIN,
            ensemble_mean_end_to_end=mean,
            n_conformations=6,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = SelfAvoidingWalk().generate(request, None, np.random.default_rng(4))
        assert result.all_successful
        orders = np.array([_planar_order(trace) for trace in result.ca_coords])
        assert orders.mean() < MAX_MEAN_PLANAR_ORDER, (
            f"mean planar order {orders.mean():.3f} for the most compact draw a multi-model "
            f"run can produce; the regression scored 0.474 here."
        )
        assert orders.max() < MAX_ANY_PLANAR_ORDER

    def test_the_default_is_the_target_itself(self) -> None:
        # Left unset, the two are the same number, which is correct whenever the engine does
        # its own spreading -- then target_end_to_end already is the mean.
        target = flory_end_to_end(80)
        plan = _WalkPlan.build(
            IDRRequest("G" * 80, 80, target, n_anchor_xyz=ORIGIN),
            0.1,
            np.random.default_rng(0),
        )
        assert plan.mean_target == pytest.approx(target)


def _pseudo_angles(trace: np.ndarray) -> np.ndarray:
    """CA-CA-CA pseudo-angles in radians, ``(n - 2,)``."""
    first = trace[:-2] - trace[1:-1]
    second = trace[2:] - trace[1:-1]
    cosine = np.sum(first * second, axis=1) / (
        np.linalg.norm(first, axis=1) * np.linalg.norm(second, axis=1)
    )
    return np.arccos(np.clip(cosine, -1.0, 1.0))


# ---------------------------------------------------------------------------
# The relaxation ladder
#
# Confirmed defect: _validate_conformer compared the measured closest non-bonded contact
# against the loosest rung the conformer had used rather than against CA_CLASH_DISTANCE,
# so the clash check could never reject what the ladder allowed. Interior regions with
# two flanking CAs in the obstacle set relaxed to 2.5-2.8 A and were returned successful.
# ---------------------------------------------------------------------------


class TestClashRelaxation:
    def test_no_obstacles_means_no_relaxation(self) -> None:
        for target_multiplier in (0.4, 1.0):
            _, result = build(
                60,
                seed=17,
                n_anchor=ORIGIN,
                target=target_multiplier * flory_end_to_end(60),
                n_conformations=10,
            )
            assert result.relaxed_to is None

    def test_the_regions_own_trace_is_never_relaxed(self) -> None:
        # Even where relaxation is legitimate -- squeezing past a placed domain -- the
        # region's own CA trace, anchors and flanking residues included, is held to
        # CA_CLASH_DISTANCE. Relaxing that is not a concession, it is an overlap.
        wall = np.stack(
            np.meshgrid(
                np.array([18.0]), np.arange(-30, 30, 2.6), np.arange(-30, 30, 2.6), indexing="ij"
            ),
            axis=-1,
        ).reshape(-1, 3)
        try:
            request, result = build(
                40, seed=42, n_anchor=ORIGIN, c_anchor=np.array([36.0, 0, 0]), obstacles=wall
            )
        except ExhaustedAttemptsError:
            return
        for coords in result.successful_coords:
            trace = full_trace(coords, request.n_anchor_xyz, request.c_anchor_xyz)
            assert closest_non_bonded(trace) >= CA_CLASH_DISTANCE

    def test_a_relaxed_conformer_is_never_reported_as_clean(self) -> None:
        # The flanking CAs of a real chain, handed in as obstacles as well as through the
        # request. Whatever the walk does with them, a result whose measured clearance is
        # below CA_CLASH_DISTANCE must carry relaxed_to.
        context, _ = flanked_interior(20, 30.0, seed=0, n_conformations=1)
        obstacles = np.vstack([context["n_prev"], context["c_next"]])
        request = IDRRequest(
            "G" * 20,
            20,
            flory_end_to_end(20),
            n_anchor_xyz=context["n_anchor"],
            c_anchor_xyz=context["c_anchor"],
            n_anchor_prev_xyz=context["n_prev"],
            c_anchor_next_xyz=context["c_next"],
            n_conformations=6,
        )
        result = SelfAvoidingWalk().generate(request, obstacles, np.random.default_rng(0))
        for coords in result.successful_coords:
            clearance = float(
                np.min(np.linalg.norm(coords[:, None, :] - obstacles[None, :, :], axis=2))
            )
            if clearance < CA_CLASH_DISTANCE:
                assert result.relaxed_to is not None
                assert clearance >= result.relaxed_to


# ---------------------------------------------------------------------------
# Impossible anchors are diagnosed at plan time
#
# Confirmed defect: min_reach(k) == 0.0 for k >= 4 made the "anchors too tight"
# GeometryError unreachable, so two anchors at the same coordinate -- two distinct
# residues on top of each other, an impossible input -- burned all 40 attempt rounds and
# then raised ExhaustedAttemptsError, blaming the attempt budget.
# ---------------------------------------------------------------------------


class TestImpossibleAnchors:
    @pytest.mark.parametrize("separation", [0.0, 1.0, 2.5])
    @pytest.mark.parametrize("n_residues", [5, 20])
    def test_anchors_too_close_are_refused_immediately(
        self, separation: float, n_residues: int
    ) -> None:
        with pytest.raises(GeometryError, match=r"double back|apart"):
            build(
                n_residues,
                n_anchor=ORIGIN,
                c_anchor=np.array([separation, 0.0, 0.0]),
                target=20.0,
                n_conformations=2,
            )

    def test_flanking_residues_that_clash_with_each_other_are_input_error(self) -> None:
        # Two fixed CAs the caller supplied, closer than two carbons can sit. The engine
        # cannot fix that by sampling, so it must not spend the attempt budget pretending.
        with pytest.raises(GeometryError, match="fixed"):
            request = IDRRequest(
                "G" * 6,
                6,
                20.0,
                n_anchor_xyz=ORIGIN,
                c_anchor_xyz=np.array([21.0, 0.0, 0.0]),
                n_anchor_prev_xyz=np.array([-3.8, 0.0, 0.0]),
                c_anchor_next_xyz=np.array([-3.8, 1.0, 0.0]),
            )
            SelfAvoidingWalk().generate(request, None, np.random.default_rng(0))


# ---------------------------------------------------------------------------
# The final audit is a real criterion, not an arithmetic ritual
#
# Confirmed defects: the end-to-end audit reused a tolerance floored at CA_CA_BOND_LENGTH,
# the same floored number the generator steered by, so it could never reject what the
# generator's corridor had allowed; and the clash audit compared against the loosest
# relaxation rung the conformer had used rather than against CA_CLASH_DISTANCE.
# ---------------------------------------------------------------------------


def plan_for(n_residues: int, target: float, *, n_anchor: np.ndarray | None = ORIGIN) -> _WalkPlan:
    """Return a resolved plan for a tail region, for testing the audit directly."""
    return _WalkPlan.build(
        IDRRequest("G" * n_residues, n_residues, target, n_anchor_xyz=n_anchor),
        0.1,
        np.random.default_rng(0),
    )


def zigzag(n_residues: int, angle: float = 125.0) -> np.ndarray:
    """Return a planar zig-zag CA trace with exact bonds and a constant pseudo-angle."""
    half = np.radians(angle) / 2.0
    coords = np.zeros((n_residues, 3))
    coords[:, 0] = CA_CA_BOND_LENGTH * np.sin(half) * np.arange(n_residues)
    coords[1::2, 1] = CA_CA_BOND_LENGTH * np.cos(half)
    return coords


class TestFinalAudit:
    def test_a_miss_smaller_than_a_bond_length_is_still_a_miss(self) -> None:
        from dodo.engines.walk import _validate_conformer

        # A valid 5-residue zig-zag spanning 10.84 A, offered against an 8.00 A target. The
        # miss is 2.84 A: larger than the 0.80 A tolerance the request implies, smaller than
        # the CA_CA_BOND_LENGTH floor the old audit applied -- so the old audit accepted it,
        # and that is exactly how a 10-residue region asked for 8.43 A returned 12.14 A with
        # success=True.
        coords = zigzag(5, angle=91.0)
        achieved = float(np.linalg.norm(coords[-1] - coords[0]))
        assert abs(achieved - 8.0) < CA_CA_BOND_LENGTH  # the old floor would have swallowed it
        plan = plan_for(5, 8.0, n_anchor=None)
        problem = _validate_conformer(coords, plan, 8.0, float("inf"), CA_CLASH_DISTANCE)
        assert problem is not None
        assert "misses this conformer's" in problem
        # The same conformer against a target it actually meets is accepted.
        assert _validate_conformer(coords, plan, achieved, float("inf"), CA_CLASH_DISTANCE) is None

    def test_an_internal_clash_is_rejected_at_every_rung_of_the_ladder(self) -> None:
        from dodo.engines.walk import _validate_conformer

        # Found by random search over cone-generated traces: every bond is exactly one bond
        # length and every pseudo-angle is inside the window (measured 97-131 degrees), yet
        # residues 3 and 8 sit 0.66 bond lengths apart. The old audit compared that against
        # the loosest rung the conformer had used, so at 2.00 A it reported this as clean --
        # the clash check could not reject what the relaxation ladder had allowed.
        #
        # The search ran at whatever CA_CA_BOND_LENGTH was then, so the trace is rescaled
        # onto whatever it is now rather than carrying a fixed bond length in its literals: a
        # uniform scale leaves every pseudo-angle untouched, and the one contact that makes
        # this fixture a fixture scales with it.
        found = np.array(
            [
                [0.0000, 0.0000, 0.0000],
                [-0.5437, -2.9466, -2.3371],
                [1.9568, -5.6236, -3.3476],
                [-0.1202, -7.1716, -6.1278],
                [0.8624, -5.9305, -9.5823],
                [2.2451, -8.1633, -12.3287],
                [0.6214, -11.5207, -11.5993],
                [-0.9614, -11.1228, -8.1676],
                [0.2904, -9.3792, -5.0318],
            ]
        )
        coords = found * (CA_CA_BOND_LENGTH / float(np.median(bond_lengths(found))))
        assert np.allclose(bond_lengths(coords), CA_CA_BOND_LENGTH, atol=1e-3)
        angles = pseudo_angles(coords)
        assert angles.min() >= BACKBONE_ANGLE_MIN and angles.max() <= BACKBONE_ANGLE_MAX
        # The premise, stated against the ladder instead of as a literal: the contact is
        # closer than CA_CLASH_DISTANCE, so it is a clash, but it is *further* apart than the
        # loosest rung -- which is what makes the two behaviours distinguishable at all. A
        # contact below the last rung would be rejected even by the audit this test replaced.
        contact = closest_non_bonded(coords)
        assert min(CLASH_RELAXATION_LADDER) < contact < CA_CLASH_DISTANCE
        target = float(np.linalg.norm(coords[-1] - coords[0]))
        plan = plan_for(9, target, n_anchor=None)
        for rung in CLASH_RELAXATION_LADDER:
            problem = _validate_conformer(coords, plan, target, float("inf"), rung)
            assert problem is not None, f"accepted a {contact:.3f} A contact at rung {rung}"
            assert "clash distance" in problem

    def test_obstacle_clearance_is_audited_against_the_rung_that_was_granted(self) -> None:
        from dodo.engines.walk import _validate_conformer

        coords = zigzag(8)
        plan = plan_for(8, float(np.linalg.norm(coords[-1] - coords[0])), n_anchor=None)
        target = float(np.linalg.norm(coords[-1] - coords[0]))
        # Clean at the strict threshold.
        assert _validate_conformer(coords, plan, target, 3.5, CA_CLASH_DISTANCE) is None
        # A measured clearance below the threshold in force is a failure, not a footnote.
        problem = _validate_conformer(coords, plan, target, 2.9, CA_CLASH_DISTANCE)
        assert problem is not None and "already-placed atom" in problem


# ---------------------------------------------------------------------------
# The obstacle set and the anchor residues
#
# base.py used to warn that leaving the anchor residues' own atoms in the obstacle set
# would "make every physically valid junction register as a clash", driving regions into
# the relaxation ladder or into ExhaustedAttemptsError. Measured, it costs no builds at
# all -- but it does make a chemically unremarkable 2.9 A CA-to-carbonyl approach consume a
# rung and be reported as relaxed geometry. This test is that measurement, so neither the
# docstring nor the reader has to guess.
# ---------------------------------------------------------------------------


class TestAnchorAtomsInTheObstacleSet:
    # Plausible N, CA, C, O for both anchor residues of an interior region spanning 40 A.
    ANCHOR_ATOMS = np.array(
        [
            [-1.2, 0.9, 0.3],
            [0.0, 0.0, 0.0],
            [1.4, 0.6, -0.2],
            [1.6, 1.8, -0.1],
            [38.8, 0.9, 0.3],
            [40.0, 0.0, 0.0],
            [41.4, 0.6, -0.2],
            [41.6, 1.8, -0.1],
        ]
    )

    @pytest.mark.parametrize("seed", [2, 5, 71])
    def test_including_them_costs_no_builds(self, seed: int) -> None:
        c_anchor = np.array([40.0, 0.0, 0.0])
        kwargs = {"n_anchor": ORIGIN, "c_anchor": c_anchor, "target": 60.0}
        with_atoms = build(20, seed=seed, obstacles=self.ANCHOR_ATOMS, n_conformations=5, **kwargs)[
            1
        ]
        without = build(20, seed=seed, n_conformations=5, **kwargs)[1]
        # The claim in the name: not one conformer is lost. The anchor's own N, C and O sit
        # ~1.5 A from its CA and the first generated residue is one 3.81 A bond from that CA,
        # so those atoms are inside the clash distance by construction -- and the cost of that
        # is a rung of the ladder (next test), never a build.
        assert with_atoms.n_successful == without.n_successful == 5
        # It is a FALSE positive, and this is what makes that word mean something: the trace
        # the relaxed build produced has no CA-CA contact inside the strict clash distance
        # anywhere, anchors included. A caller reading relaxed_to can tell this apart from a
        # real squeeze because there is no squeeze in the CA trace at all.
        for coords in with_atoms.ca_coords:
            trace = full_trace(coords, ORIGIN, c_anchor)
            assert closest_non_bonded(trace) >= CA_CLASH_DISTANCE
        # Attempt counts are NOT equal between the two runs, and asserting that they were was
        # over-strong -- it held only for the seeds it was written against. Relaxing a rung
        # admits candidates the strict pass rejected, so the two searches diverge and the
        # restart count moves in either direction: measured over 200 seeds the difference spans
        # -2 to +4 restarts, with no build lost on either side at any seed. What must hold is
        # that neither search spends its restart budget on the difference -- over those 200
        # seeds the worst either side spent was 5 of the MAX_ATTEMPTS_PER_REGION it may spend,
        # so require a quarter of the budget to still be untouched.
        assert 4 * with_atoms.attempts <= MAX_ATTEMPTS_PER_REGION
        assert 4 * without.attempts <= MAX_ATTEMPTS_PER_REGION

    def test_including_them_biases_the_junction_away_from_real_geometry(self) -> None:
        """The real cost of including them, now that there is no ladder to absorb it.

        This test used to assert that including the anchor's own atoms consumed a rung of
        CLASH_RELAXATION_LADDER. That ladder is gone -- it was responsible for 69 of 79 steric
        clashes across the fixture sweep while preventing zero region failures -- so the cost
        moved. It is no longer a spurious ``relaxed_to``; it is a distorted junction.

        A residue bonded to the anchor legitimately comes closer to the anchor's backbone than
        the clash distance. Measured over 649,658 sequence-neighbour pairs from the human
        proteome, a CA sits at 2.379 A from the next residue's N at the 0.1st percentile and
        3.280 A at the median -- so the median real junction is INSIDE the 3.20 A clash
        distance. Feeding those atoms to the engine as ordinary obstacles therefore does not
        prevent a defect, it forbids the commonest real geometry.
        """
        kwargs = {"n_anchor": ORIGIN, "c_anchor": np.array([40.0, 0.0, 0.0]), "target": 60.0}
        anchor_backbone = self.ANCHOR_ATOMS[[0, 2, 3]]  # the N-terminal anchor's N, C, O

        def closest_first_ca_to_anchor_backbone(*, obstacles: np.ndarray | None) -> float:
            worst = np.inf
            for seed in (2, 5, 71):
                result = build(20, seed=seed, obstacles=obstacles, n_conformations=5, **kwargs)[1]
                for coords in result.ca_coords:
                    d = np.linalg.norm(anchor_backbone - coords[0], axis=1).min()
                    worst = min(worst, float(d))
            return worst

        excluded = closest_first_ca_to_anchor_backbone(obstacles=None)
        included = closest_first_ca_to_anchor_backbone(obstacles=self.ANCHOR_ATOMS)

        # Included: the engine is forced to hold the whole region off the anchor's backbone.
        assert included >= CA_CLASH_DISTANCE
        # Excluded: it is free to adopt the natural approach, and does.
        assert excluded < CA_CLASH_DISTANCE, (
            f"expected the first CA to approach the anchor backbone naturally, got {excluded:.2f}"
        )


class TestChainClearMask:
    """Lock the vectorized clash cull to the per-conformer KD-tree mask.

    The cull must reproduce, bit for bit, the mask the per-conformer
    ``cKDTree(chain_points).query(...) >= CA_CLASH_DISTANCE`` loop produced -- that is the
    invariant that lets it replace the walk's single hottest inner loop.
    """

    @staticmethod
    def _tree_mask(candidates: np.ndarray, chain_points: np.ndarray) -> np.ndarray:
        """Return the old per-conformer KD-tree clash mask, the reference oracle here."""
        n_live, count, _ = candidates.shape
        nearest = np.full((n_live, count), np.inf)
        if chain_points.shape[1] == 0:
            return nearest >= CA_CLASH_DISTANCE
        for row in range(n_live):
            valid = np.all(np.isfinite(candidates[row]), axis=1)
            if not np.any(valid):
                continue
            found, _ = cKDTree(chain_points[row]).query(candidates[row][valid])
            partial = np.full(count, np.inf)
            partial[valid] = found
            nearest[row] = partial
        return nearest >= CA_CLASH_DISTANCE

    def test_matches_the_kdtree_mask_including_the_clash_boundary(self) -> None:
        rng = np.random.default_rng(7)
        mismatches = 0
        for _ in range(120):
            n_live = int(rng.integers(1, 10))
            count = 355
            n_points = int(rng.integers(1, 50))
            centres = rng.normal(size=(n_live, 3)) * 10.0
            # Candidates exactly one bond length from the apex, as the real caller supplies.
            directions = rng.normal(size=(n_live, count, 3))
            directions /= np.linalg.norm(directions, axis=2, keepdims=True)
            candidates = centres[:, None, :] + CA_CA_BOND_LENGTH * directions
            chain = rng.normal(size=(n_live, n_points, 3)) * 12.0 + centres[:, None, :]
            # Plant a point a hair either side of the clash radius from a real candidate, so
            # the two code paths are forced to agree exactly on the >= boundary.
            for row in range(n_live):
                col = int(rng.integers(0, count))
                unit = rng.normal(size=3)
                unit /= np.linalg.norm(unit)
                offset = CA_CLASH_DISTANCE + rng.normal() * 1e-9
                chain[row, int(rng.integers(0, n_points))] = candidates[row, col] + unit * offset
            expected = self._tree_mask(candidates, chain)
            got = SelfAvoidingWalk._chain_clear_mask(candidates, chain, centres)
            mismatches += int(np.count_nonzero(expected != got))
        assert mismatches == 0

    def test_no_chain_points_is_all_clear(self) -> None:
        candidates = np.zeros((3, 355, 3))
        empty = np.zeros((3, 0, 3))
        assert np.all(SelfAvoidingWalk._chain_clear_mask(candidates, empty, np.zeros((3, 3))))

    def test_missing_centres_is_refused(self) -> None:
        candidates = np.zeros((2, 355, 3))
        chain = np.ones((2, 4, 3))
        with pytest.raises(EngineError, match="candidate centres"):
            SelfAvoidingWalk._chain_clear_mask(candidates, chain, None)
