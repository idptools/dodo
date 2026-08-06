"""Tests for the geometry kernel: transforms, sampling and metrics.

Every builder in DODO stands on these three modules, so this file is deliberately
paranoid. Several tests are named for confirmed defects in the pre-rewrite code, which
had 8 of its 20 geometry functions with no caller at all -- and that untested half
contained the worst of the bugs, including a "rotation" with determinant -1 and a
sampler that returned all-NaN.
"""

from __future__ import annotations

import math
import time

import numpy as np
import pytest
from scipy import stats
from scipy.spatial.transform import Rotation

from dodo.constants import (
    BACKBONE_ANGLE_IDEAL,
    BACKBONE_ANGLE_MAX,
    BACKBONE_ANGLE_MIN,
    CA_CA_BOND_LENGTH,
    CA_CLASH_DISTANCE,
    CANDIDATES_PER_ANGLE,
    backbone_angle_grid,
    flory_end_to_end,
)
from dodo.exceptions import GeometryError
from dodo.geometry.metrics import (
    TraceReport,
    ca_bond_lengths,
    ca_pseudo_angles,
    end_to_end,
    radius_of_gyration,
    validate_ca_trace,
)
from dodo.geometry.sampling import (
    cone_candidates,
    cone_candidates_batch,
    cone_template_cache_info,
    random_points_on_sphere,
    random_unit_vectors,
)
from dodo.geometry.transforms import (
    align_frame,
    apply,
    rotation_between_vectors,
    rotation_between_vectors_batch,
    rotation_from_axis_angle,
)
from dodo.structure import Structure

# Directions that between them cover every branch of rotation_between_vectors: the
# world axes (where a naive "cross with x" perpendicular construction fails), a generic
# direction, and a near-degenerate one.
DIRECTIONS = [
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (-1.0, 0.0, 0.0),
    (0.0, 0.0, -1.0),
    (1.0, 1.0, 1.0),
    (0.31, -2.7, 0.004),
    (1e-6, 0.0, 1.0),
]


def is_proper_rotation(matrix: np.ndarray) -> bool:
    """Return True if ``matrix`` is orthonormal with determinant +1."""
    return bool(
        matrix.shape == (3, 3)
        and np.isclose(np.linalg.det(matrix), 1.0, atol=1e-12)
        and np.allclose(matrix @ matrix.T, np.eye(3), atol=1e-12)
    )


# ---------------------------------------------------------------------------
# transforms: rotation_from_axis_angle
# ---------------------------------------------------------------------------


class TestRotationFromAxisAngle:
    def test_quarter_turn_about_z(self) -> None:
        matrix = rotation_from_axis_angle([0, 0, 1], math.pi / 2)
        assert np.allclose(matrix @ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
        assert np.allclose(matrix @ [0.0, 0.0, 1.0], [0.0, 0.0, 1.0])

    def test_matches_scipy(self) -> None:
        axis = np.array([0.3, -0.7, 0.2])
        angle = 1.234
        expected = Rotation.from_rotvec(axis / np.linalg.norm(axis) * angle).as_matrix()
        assert np.allclose(rotation_from_axis_angle(axis, angle), expected)

    def test_axis_need_not_be_normalized(self) -> None:
        a = rotation_from_axis_angle([0, 0, 7.5], 0.4)
        b = rotation_from_axis_angle([0, 0, 1.0], 0.4)
        assert np.allclose(a, b)

    def test_zero_angle_is_identity(self) -> None:
        assert np.allclose(rotation_from_axis_angle([1, 2, 3], 0.0), np.eye(3))

    def test_negative_angle_is_the_inverse(self) -> None:
        forward = rotation_from_axis_angle([1, 2, 3], 0.9)
        backward = rotation_from_axis_angle([1, 2, 3], -0.9)
        assert np.allclose(forward @ backward, np.eye(3))

    @pytest.mark.parametrize("angle", [0.0, 1e-12, 0.5, math.pi / 2, math.pi, 2 * math.pi, -3.0])
    @pytest.mark.parametrize("axis", DIRECTIONS)
    def test_always_a_proper_rotation(self, axis: tuple[float, ...], angle: float) -> None:
        assert is_proper_rotation(rotation_from_axis_angle(axis, angle))

    def test_zero_axis_raises_rather_than_returning_nan(self) -> None:
        # The pre-rewrite code divided by the norm unguarded and produced a NaN matrix.
        with pytest.raises(GeometryError, match="zero length"):
            rotation_from_axis_angle([0.0, 0.0, 0.0], 1.0)

    def test_nan_axis_raises(self) -> None:
        with pytest.raises(GeometryError, match="non-finite"):
            rotation_from_axis_angle([np.nan, 0.0, 1.0], 1.0)

    def test_nan_angle_raises(self) -> None:
        with pytest.raises(GeometryError, match="finite"):
            rotation_from_axis_angle([0.0, 0.0, 1.0], np.nan)

    def test_wrong_shape_raises(self) -> None:
        with pytest.raises(GeometryError, match="3-vector"):
            rotation_from_axis_angle([1.0, 0.0], 1.0)


# ---------------------------------------------------------------------------
# transforms: rotation_between_vectors
# ---------------------------------------------------------------------------


class TestRotationBetweenVectors:
    @pytest.mark.parametrize("b", DIRECTIONS)
    @pytest.mark.parametrize("a", DIRECTIONS)
    def test_maps_a_onto_b_and_is_proper(self, a: tuple[float, ...], b: tuple[float, ...]) -> None:
        matrix = rotation_between_vectors(a, b)
        assert is_proper_rotation(matrix), f"det = {np.linalg.det(matrix)}"
        unit_a = np.asarray(a) / np.linalg.norm(a)
        unit_b = np.asarray(b) / np.linalg.norm(b)
        assert np.allclose(matrix @ unit_a, unit_b, atol=1e-9)

    def test_identical_vectors_give_identity(self) -> None:
        a = np.array([0.3, 0.4, -1.2])
        assert np.allclose(rotation_between_vectors(a, a), np.eye(3))
        assert is_proper_rotation(rotation_between_vectors(a, a))

    def test_parallel_but_different_magnitudes_give_identity(self) -> None:
        assert np.allclose(rotation_between_vectors([0, 0, 1], [0, 0, 42.0]), np.eye(3))

    @pytest.mark.parametrize("a", DIRECTIONS)
    def test_antiparallel_is_a_rotation_not_a_point_inversion(self, a: tuple[float, ...]) -> None:
        """The headline regression test.

        ``find_transform`` returned ``-np.eye(3)`` for the antiparallel case. That maps
        ``a`` to ``-a`` correctly, so it looks right, but its determinant is -1: it is a
        point inversion that mirrors every component perpendicular to ``a`` as well.
        """
        matrix = rotation_between_vectors(a, tuple(-x for x in a))
        unit_a = np.asarray(a) / np.linalg.norm(a)

        assert is_proper_rotation(matrix), f"det = {np.linalg.det(matrix)}"
        assert np.allclose(matrix @ unit_a, -unit_a, atol=1e-9)
        assert not np.allclose(matrix, -np.eye(3)), "returned the old point-inversion matrix"

    def test_antiparallel_rotation_axis_is_perpendicular_to_a(self) -> None:
        a = np.array([0.0, 0.0, 1.0])
        matrix = rotation_between_vectors(a, -a)
        rotvec = Rotation.from_matrix(matrix).as_rotvec()
        assert np.isclose(np.linalg.norm(rotvec), math.pi)
        assert np.isclose(rotvec @ a, 0.0, atol=1e-12)

    def test_antiparallel_preserves_perpendicular_lengths(self) -> None:
        """A mirror would flip a perpendicular component; a rotation only moves it."""
        a = np.array([0.0, 0.0, 1.0])
        matrix = rotation_between_vectors(a, -a)
        perpendicular = np.array([1.0, 0.0, 0.0])
        image = matrix @ perpendicular
        assert np.isclose(np.linalg.norm(image), 1.0)
        assert np.isclose(image @ a, 0.0, atol=1e-12)

    def test_near_antiparallel_is_still_accurate(self) -> None:
        a = np.array([0.0, 0.0, 1.0])
        b = np.array([1e-9, 0.0, -1.0])
        matrix = rotation_between_vectors(a, b)
        assert is_proper_rotation(matrix)
        assert np.allclose(matrix @ a, b / np.linalg.norm(b), atol=1e-8)

    @pytest.mark.parametrize(
        "bad", [([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]), ([1.0, 0.0, 0.0], [0, 0, 0])]
    )
    def test_zero_vector_raises(self, bad: tuple[list[float], list[float]]) -> None:
        with pytest.raises(GeometryError, match="zero length"):
            rotation_between_vectors(*bad)

    def test_nan_raises(self) -> None:
        with pytest.raises(GeometryError, match="non-finite"):
            rotation_between_vectors([1.0, 0.0, 0.0], [np.nan, 1.0, 0.0])

    def test_random_pairs_are_all_proper(self) -> None:
        rng = np.random.default_rng(20260730)
        for _ in range(500):
            a, b = rng.normal(size=(2, 3))
            matrix = rotation_between_vectors(a, b)
            assert np.isclose(np.linalg.det(matrix), 1.0)
            assert np.allclose(matrix @ (a / np.linalg.norm(a)), b / np.linalg.norm(b), atol=1e-9)


class TestRotationBetweenVectorsBatch:
    """Lock the batched rotation to the scalar one, row by row.

    The batched form must equal the scalar one bit for bit; that equality is the whole
    reason it is allowed to replace the walk engine's per-conformer rotation loop.
    """

    def _targets(self) -> np.ndarray:
        rng = np.random.default_rng(4242)
        generic = rng.normal(size=(2000, 3))
        # Inject the branch-boundary rows: exact +z / -z parallel and antiparallel, tiny
        # perturbations either side of antiparallel, and a longer antiparallel vector.
        edge = np.array(
            [
                [0.0, 0.0, 1.0],
                [0.0, 0.0, -1.0],
                [1e-10, 0.0, 1.0],
                [1e-10, 0.0, -1.0],
                [0.0, 0.0, -5.3],
            ]
        )
        return np.vstack([edge, generic])

    def test_fixed_source_matches_scalar_bit_for_bit(self) -> None:
        source = np.array([0.0, 0.0, 1.0])
        targets = self._targets()
        batch = rotation_between_vectors_batch(source, targets)
        for i, target in enumerate(targets):
            assert np.array_equal(batch[i], rotation_between_vectors(source, target))

    def test_per_row_source_matches_scalar_bit_for_bit(self) -> None:
        rng = np.random.default_rng(99)
        sources = rng.normal(size=(500, 3))
        targets = rng.normal(size=(500, 3))
        batch = rotation_between_vectors_batch(sources, targets)
        for i in range(sources.shape[0]):
            assert np.array_equal(batch[i], rotation_between_vectors(sources[i], targets[i]))

    def test_zero_length_row_raises(self) -> None:
        targets = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        with pytest.raises(GeometryError, match="zero-length"):
            rotation_between_vectors_batch([0.0, 0.0, 1.0], targets)

    def test_wrong_shape_raises(self) -> None:
        with pytest.raises(GeometryError, match=r"shape \(n, 3\)"):
            rotation_between_vectors_batch([0.0, 0.0, 1.0], np.zeros((4, 2)))


# ---------------------------------------------------------------------------
# transforms: align_frame
# ---------------------------------------------------------------------------


class TestAlignFrame:
    @pytest.mark.parametrize("direction", DIRECTIONS)
    def test_z_column_is_the_direction(self, direction: tuple[float, ...]) -> None:
        frame = align_frame((0, 0, 0), direction)
        expected = np.asarray(direction) / np.linalg.norm(direction)
        assert np.allclose(frame[:, 2], expected)
        assert np.allclose(frame @ [0.0, 0.0, 1.0], expected)

    @pytest.mark.parametrize("direction", DIRECTIONS)
    def test_is_a_proper_orthonormal_frame(self, direction: tuple[float, ...]) -> None:
        assert is_proper_rotation(align_frame((1.0, -2.0, 3.0), direction))

    def test_direction_is_relative_to_origin(self) -> None:
        origin = np.array([5.0, 5.0, 5.0])
        target = np.array([5.0, 5.0, 9.0])
        frame = align_frame(origin, target)
        assert np.allclose(frame[:, 2], [0.0, 0.0, 1.0])

    def test_right_handed(self) -> None:
        frame = align_frame((0, 0, 0), (0.2, 0.5, -0.9))
        assert np.allclose(np.cross(frame[:, 0], frame[:, 1]), frame[:, 2], atol=1e-12)

    def test_deterministic(self) -> None:
        a = align_frame((0, 0, 0), (0.2, 0.5, -0.9))
        b = align_frame((0, 0, 0), (0.2, 0.5, -0.9))
        assert np.array_equal(a, b)

    def test_maps_local_geometry_into_world(self) -> None:
        origin = np.array([1.0, 2.0, 3.0])
        frame = align_frame(origin, origin + np.array([0.0, 1.0, 0.0]))
        # A point one Angstrom "ahead" in the frame lands one Angstrom along +y.
        assert np.allclose(apply([0.0, 0.0, 1.0], frame) + origin, [1.0, 3.0, 3.0])

    def test_direction_is_a_point_not_a_vector(self) -> None:
        """The documented, deliberate reading: z runs along ``direction - origin``."""
        origin = np.array([0.0, 0.0, 5.0])
        as_point = align_frame(origin, (0.0, 0.0, 6.0))
        assert np.allclose(as_point[:, 2], [0.0, 0.0, 1.0])
        # Read as a vector it would have been (0, 0, 6) -> also +z, so use a case where
        # the two readings genuinely differ.
        assert np.allclose(align_frame(origin, (0.0, 0.0, 4.0))[:, 2], [0.0, 0.0, -1.0])

    def test_coincident_origin_and_direction_raises(self) -> None:
        with pytest.raises(GeometryError, match="coincides"):
            align_frame((1.0, 1.0, 1.0), (1.0, 1.0, 1.0))

    def test_nan_raises(self) -> None:
        with pytest.raises(GeometryError, match="non-finite"):
            align_frame((0, 0, 0), (np.inf, 0.0, 0.0))


# ---------------------------------------------------------------------------
# transforms: apply
# ---------------------------------------------------------------------------


class TestApply:
    def test_rotates_row_vectors(self) -> None:
        coords = np.array([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        matrix = rotation_from_axis_angle([0, 0, 1], math.pi / 2)
        assert np.allclose(apply(coords, matrix), [[0.0, 1.0, 0.0], [-2.0, 0.0, 0.0]])

    def test_single_point_keeps_its_shape(self) -> None:
        matrix = rotation_from_axis_angle([0, 0, 1], math.pi / 2)
        result = apply([1.0, 0.0, 0.0], matrix)
        assert result.shape == (3,)
        assert np.allclose(result, [0.0, 1.0, 0.0])

    def test_about_point_is_fixed(self) -> None:
        centre = np.array([4.0, -1.0, 2.0])
        coords = np.array([centre, centre + np.array([1.0, 0.0, 0.0])])
        matrix = rotation_from_axis_angle([0, 0, 1], 0.77)
        rotated = apply(coords, matrix, about=centre)
        assert np.allclose(rotated[0], centre)
        assert np.isclose(np.linalg.norm(rotated[1] - centre), 1.0)

    def test_default_is_the_world_origin_not_the_centroid(self) -> None:
        coords = np.array([[10.0, 0.0, 0.0], [11.0, 0.0, 0.0]])
        matrix = rotation_from_axis_angle([0, 0, 1], math.pi)
        assert np.allclose(apply(coords, matrix), [[-10.0, 0.0, 0.0], [-11.0, 0.0, 0.0]])

    def test_does_not_modify_input(self) -> None:
        coords = np.array([[1.0, 2.0, 3.0]])
        original = coords.copy()
        apply(coords, rotation_from_axis_angle([0, 1, 0], 1.0))
        assert np.array_equal(coords, original)

    def test_round_trip(self) -> None:
        rng = np.random.default_rng(7)
        coords = rng.normal(size=(20, 3))
        matrix = rotation_between_vectors([1.0, 0.0, 0.0], [0.3, 0.5, -0.8])
        assert np.allclose(apply(apply(coords, matrix), matrix.T), coords)

    def test_accepts_a_read_only_input_and_returns_a_writable_array(self) -> None:
        coords = np.zeros((4, 3))
        coords.setflags(write=False)
        out = apply(coords, np.eye(3))
        out[0, 0] = 1.0  # must not raise
        assert out.flags.writeable

    def test_reflection_raises(self) -> None:
        """A determinant of -1 silently mirrors the structure, so refuse it."""
        mirror = np.diag([1.0, 1.0, -1.0])
        with pytest.raises(GeometryError, match="determinant"):
            apply(np.zeros((2, 3)), mirror)

    def test_point_inversion_raises(self) -> None:
        with pytest.raises(GeometryError, match="determinant"):
            apply(np.zeros((2, 3)), -np.eye(3))

    def test_scaling_matrix_raises(self) -> None:
        with pytest.raises(GeometryError, match="determinant"):
            apply(np.zeros((2, 3)), 2.0 * np.eye(3))

    def test_bad_rotation_shape_raises(self) -> None:
        with pytest.raises(GeometryError, match=r"\(3, 3\)"):
            apply(np.zeros((2, 3)), np.eye(4))

    def test_bad_coords_shape_raises(self) -> None:
        with pytest.raises(GeometryError, match=r"\(n, 3\)"):
            apply(np.zeros((2, 4)), np.eye(3))

    def test_nan_coords_raise(self) -> None:
        coords = np.zeros((3, 3))
        coords[1, 0] = np.nan
        with pytest.raises(GeometryError, match="non-finite"):
            apply(coords, np.eye(3))

    def test_volume_preserving_shear_raises(self) -> None:
        """Det == +1 does not mean rigid: a shear passes that test and stretches bonds.

        Measured before this was checked: this shear turned a trace with uniform 4.294 A
        bonds into bonds of 5.385 / 3.280 / 3.280 / 5.385 A, and apply() returned it with
        no complaint. A caller-supplied matrix is exactly the case apply() exists to
        accept, so the determinant alone is not a sufficient gate.
        """
        shear = np.array([[1.0, 0.6, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        assert np.isclose(np.linalg.det(shear), 1.0)
        with pytest.raises(GeometryError, match="orthonormal"):
            apply(np.array([[1.0, 2.0, 3.0], [4.0, -1.0, 0.5]]), shear)

    def test_shear_does_not_reach_the_coordinates(self) -> None:
        """The shear is refused before any bond length changes."""
        trace = straight_trace(5)
        shear = np.array([[1.0, 0.0, 0.0], [0.4, 1.0, 0.0], [0.0, 0.0, 1.0]])
        with pytest.raises(GeometryError):
            apply(trace, shear)
        assert np.allclose(ca_bond_lengths(trace), CA_CA_BOND_LENGTH)

    def test_a_real_rotation_still_passes_the_orthonormality_check(self) -> None:
        """The stricter gate must not reject the matrices this module itself produces."""
        rng = np.random.default_rng(3)
        coords = rng.normal(size=(8, 3))
        for _ in range(50):
            matrix = rotation_between_vectors(rng.normal(size=3), rng.normal(size=3))
            assert np.allclose(
                np.linalg.norm(apply(coords, matrix), axis=1), np.linalg.norm(coords, axis=1)
            )


# ---------------------------------------------------------------------------
# sampling: random_unit_vectors
# ---------------------------------------------------------------------------


def naive_polar_unit_vectors(n: int, rng: np.random.Generator) -> np.ndarray:
    """Sample a sphere the WRONG way, with uniform theta and phi.

    Present so the uniformity tests can be shown to have teeth -- a statistical check
    that does not reject this is not testing anything.
    """
    theta = rng.uniform(0.0, np.pi, n)
    phi = rng.uniform(0.0, 2.0 * np.pi, n)
    return np.column_stack(
        (np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta))
    )


class TestRandomUnitVectors:
    def test_shape_and_norms(self) -> None:
        vectors = random_unit_vectors(1000, np.random.default_rng(0))
        assert vectors.shape == (1000, 3)
        assert np.allclose(np.linalg.norm(vectors, axis=1), 1.0)

    def test_n_one_is_a_normal_case(self) -> None:
        vectors = random_unit_vectors(1, np.random.default_rng(0))
        assert vectors.shape == (1, 3)

    def test_n_zero_gives_an_empty_array(self) -> None:
        assert random_unit_vectors(0, np.random.default_rng(0)).shape == (0, 3)

    def test_negative_n_raises(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            random_unit_vectors(-1, np.random.default_rng(0))

    @pytest.mark.parametrize("n", [2.9, 0.5, np.float64(3.5)])
    def test_non_integer_n_raises_rather_than_truncating(self, n: float) -> None:
        """2.9 vectors is not a request that can be honoured, so do not guess.

        Before this check ``random_unit_vectors(2.9)`` silently returned 2 vectors. A
        caller computing a count from a length or a ratio gets one fewer sample than it
        asked for, in a stochastic function where the shortfall is invisible.
        """
        with pytest.raises(TypeError, match="integer"):
            random_unit_vectors(n, np.random.default_rng(0))  # type: ignore[arg-type]

    @pytest.mark.parametrize("n", [np.int64(4), np.uint8(4), True])
    def test_integer_like_counts_are_accepted(self, n: int) -> None:
        """Numpy integers are integers; only fractional counts are refused."""
        assert random_unit_vectors(n, np.random.default_rng(0)).shape == (int(n), 3)

    def test_reproducible_from_a_seed(self) -> None:
        a = random_unit_vectors(50, np.random.default_rng(12345))
        b = random_unit_vectors(50, np.random.default_rng(12345))
        assert np.array_equal(a, b)

    def test_different_seeds_differ(self) -> None:
        a = random_unit_vectors(50, np.random.default_rng(1))
        b = random_unit_vectors(50, np.random.default_rng(2))
        assert not np.allclose(a, b)

    def test_rejects_a_bare_seed(self) -> None:
        with pytest.raises(TypeError, match="default_rng"):
            random_unit_vectors(10, 1234)  # type: ignore[arg-type]

    def test_rejects_legacy_randomstate(self) -> None:
        with pytest.raises(TypeError, match="Generator"):
            random_unit_vectors(10, np.random.RandomState(0))  # type: ignore[arg-type]

    def test_does_not_touch_the_global_numpy_state(self) -> None:
        """Rule 3: nothing in DODO may consume the global random stream."""
        # The legacy global API is used deliberately here: this test exists to prove
        # DODO never touches it, and there is no other way to observe that.
        np.random.seed(0)  # noqa: NPY002
        control = np.random.random(5)  # noqa: NPY002
        np.random.seed(0)  # noqa: NPY002
        random_unit_vectors(10_000, np.random.default_rng(0))
        assert np.array_equal(np.random.random(5), control)  # noqa: NPY002

    def test_z_component_is_uniform_on_minus_one_to_one(self) -> None:
        """Archimedes: the z component of a uniform sphere point is Uniform(-1, 1).

        This is the check that discriminates. Uniform-theta sampling piles points at the
        poles, so its z distribution is arcsine-shaped, not uniform.
        """
        vectors = random_unit_vectors(20_000, np.random.default_rng(99))
        result = stats.kstest(vectors[:, 2], stats.uniform(loc=-1.0, scale=2.0).cdf)
        assert result.pvalue > 0.01, f"z component is not Uniform(-1, 1): p = {result.pvalue}"

    def test_the_same_check_rejects_the_naive_sampler(self) -> None:
        naive = naive_polar_unit_vectors(20_000, np.random.default_rng(99))
        result = stats.kstest(naive[:, 2], stats.uniform(loc=-1.0, scale=2.0).cdf)
        assert result.pvalue < 1e-10, "the uniformity test has no teeth"

    def test_equal_area_bands_are_equally_populated(self) -> None:
        vectors = random_unit_vectors(20_000, np.random.default_rng(7))
        # Bands of equal z-extent have equal area on a sphere, so counts should be flat.
        edges = np.linspace(-1.0, 1.0, 11)
        counts, _ = np.histogram(vectors[:, 2], bins=edges)
        result = stats.chisquare(counts)
        assert result.pvalue > 0.01, f"equal-area bands unevenly filled: {counts}"

        naive_counts, _ = np.histogram(
            naive_polar_unit_vectors(20_000, np.random.default_rng(7))[:, 2], bins=edges
        )
        assert stats.chisquare(naive_counts).pvalue < 1e-10

    def test_mean_direction_is_near_zero(self) -> None:
        vectors = random_unit_vectors(50_000, np.random.default_rng(3))
        assert np.allclose(vectors.mean(axis=0), 0.0, atol=0.02)


# ---------------------------------------------------------------------------
# sampling: random_points_on_sphere
# ---------------------------------------------------------------------------


class TestRandomPointsOnSphere:
    def test_points_lie_on_the_sphere(self) -> None:
        centre = np.array([1.0, -2.0, 30.0])
        points = random_points_on_sphere(centre, 12.5, 500, np.random.default_rng(0))
        assert points.shape == (500, 3)
        assert np.allclose(np.linalg.norm(points - centre, axis=1), 12.5)

    def test_n_zero(self) -> None:
        assert random_points_on_sphere((0, 0, 0), 1.0, 0, np.random.default_rng(0)).shape == (0, 3)

    def test_reproducible(self) -> None:
        a = random_points_on_sphere((0, 0, 0), 3.8, 20, np.random.default_rng(5))
        b = random_points_on_sphere((0, 0, 0), 3.8, 20, np.random.default_rng(5))
        assert np.array_equal(a, b)

    @pytest.mark.parametrize("radius", [0.0, -1.0, np.nan, np.inf])
    def test_bad_radius_raises(self, radius: float) -> None:
        with pytest.raises(GeometryError, match="radius"):
            random_points_on_sphere((0, 0, 0), radius, 5, np.random.default_rng(0))

    def test_bad_centre_raises(self) -> None:
        with pytest.raises(GeometryError, match="3-vector"):
            random_points_on_sphere((0, 0), 1.0, 5, np.random.default_rng(0))

    def test_nan_centre_raises(self) -> None:
        with pytest.raises(GeometryError, match="non-finite"):
            random_points_on_sphere((np.nan, 0, 0), 1.0, 5, np.random.default_rng(0))


# ---------------------------------------------------------------------------
# sampling: cone_candidates
# ---------------------------------------------------------------------------


def measured_angle(previous: np.ndarray, current: np.ndarray, candidate: np.ndarray) -> float:
    """CA-CA-CA angle in degrees, computed independently of the metrics module."""
    a = np.asarray(previous) - np.asarray(current)
    b = np.asarray(candidate) - np.asarray(current)
    cosine = (a @ b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def grow_with_selection_policy(n: int, seed: int, policy: str) -> np.ndarray:
    """Grow an ``n``-residue trace from cone candidates under one selection policy.

    ``policy="first"`` takes the first non-clashing candidate in the emitted order, which
    is the policy :func:`cone_candidates` used to recommend in its docstring.
    ``policy="sampled"`` draws uniformly among the non-clashing candidates. Both see the
    identical candidate set and the identical clash rule at every step, so any difference
    between the chains they produce is attributable to the selection policy alone.
    """
    rng = np.random.default_rng(seed)
    coords = np.zeros((n, 3))
    coords[1] = CA_CA_BOND_LENGTH * random_unit_vectors(1, rng)[0]
    for i in range(2, n):
        candidates = cone_candidates(coords[i - 2], coords[i - 1])
        if i > 2:
            # i-1 is bonded and i-2 is angle-constrained, so neither counts as a clash.
            separation = np.linalg.norm(
                candidates[:, None, :] - coords[None, : i - 2, :], axis=2
            ).min(axis=1)
            free = np.flatnonzero(separation >= CA_CLASH_DISTANCE)
        else:
            free = np.arange(candidates.shape[0])
        if free.size == 0:
            raise AssertionError(f"boxed in at residue {i} with policy {policy!r}")
        index = free[0] if policy == "first" else free[int(rng.integers(free.size))]
        coords[i] = candidates[index]
    return coords


class TestConeCandidates:
    @pytest.mark.parametrize("direction", DIRECTIONS)
    def test_every_candidate_is_one_bond_from_the_apex(self, direction: tuple[float, ...]) -> None:
        unit = np.asarray(direction) / np.linalg.norm(direction)
        previous = np.array([10.0, -5.0, 2.0])
        current = previous + CA_CA_BOND_LENGTH * unit
        candidates = cone_candidates(previous, current)
        distances = np.linalg.norm(candidates - current, axis=1)
        assert np.allclose(distances, CA_CA_BOND_LENGTH, atol=1e-9)

    @pytest.mark.parametrize("direction", DIRECTIONS)
    def test_every_candidate_has_the_requested_pseudo_angle(
        self, direction: tuple[float, ...]
    ) -> None:
        unit = np.asarray(direction) / np.linalg.norm(direction)
        previous = np.array([0.0, 0.0, 0.0])
        current = CA_CA_BOND_LENGTH * unit
        grid = backbone_angle_grid()
        candidates = cone_candidates(previous, current)

        expected = np.repeat(grid, CANDIDATES_PER_ANGLE)
        measured = np.array([measured_angle(previous, current, c) for c in candidates])
        assert np.allclose(measured, expected, atol=1e-8)

    def test_antiparallel_axis_does_not_invert_the_angles(self) -> None:
        """Regression: the -I "rotation" turned every angle theta into 180 - theta.

        The template is built around +z, so a chain segment running along -z is the exact
        antiparallel case. With ``-np.eye(3)`` the axial component of every candidate is
        negated, and 125 degrees comes back as 55.
        """
        previous = np.array([0.0, 0.0, 0.0])
        current = np.array([0.0, 0.0, -CA_CA_BOND_LENGTH])
        candidates = cone_candidates(previous, current, angles=[125.0], per_angle=6)
        measured = np.array([measured_angle(previous, current, c) for c in candidates])
        assert np.allclose(measured, 125.0, atol=1e-9)
        assert not np.any(np.isclose(measured, 55.0, atol=1.0))

    def test_ideal_first_ordering_is_preserved(self) -> None:
        previous = np.array([0.0, 0.0, 0.0])
        current = np.array([CA_CA_BOND_LENGTH, 0.0, 0.0])
        candidates = cone_candidates(previous, current)

        first_ring = candidates[:CANDIDATES_PER_ANGLE]
        measured = [measured_angle(previous, current, c) for c in first_ring]
        assert np.allclose(measured, BACKBONE_ANGLE_IDEAL, atol=1e-8)

        # And the whole sequence is ordered by distance from ideal, which is what makes
        # "take the first non-clashing candidate" prefer realistic geometry. The
        # pre-rewrite vectorized generator emitted 91, 92, 93... and lost this.
        all_measured = np.array([measured_angle(previous, current, c) for c in candidates])
        per_ring = all_measured.reshape(-1, CANDIDATES_PER_ANGLE)[:, 0]
        deviation = np.abs(per_ring - BACKBONE_ANGLE_IDEAL)
        assert np.all(np.diff(deviation) >= -1e-9), f"not ideal-first: {per_ring[:6]}"

    def test_custom_angle_order_is_respected(self) -> None:
        previous = np.array([0.0, 0.0, 0.0])
        current = np.array([0.0, CA_CA_BOND_LENGTH, 0.0])
        requested = [161.0, 91.0, 125.0]
        candidates = cone_candidates(previous, current, angles=requested, per_angle=2)
        measured = np.array([measured_angle(previous, current, c) for c in candidates])
        assert np.allclose(measured, np.repeat(requested, 2), atol=1e-8)

    def test_candidate_count(self) -> None:
        candidates = cone_candidates([0, 0, 0], [3.8, 0, 0])
        assert candidates.shape == (backbone_angle_grid().size * CANDIDATES_PER_ANGLE, 3)
        assert cone_candidates([0, 0, 0], [3.8, 0, 0], angles=[100.0], per_angle=1).shape == (1, 3)

    def test_rings_are_azimuthally_offset_from_each_other(self) -> None:
        """Without a per-ring offset all candidates lie in `per_angle` meridional planes."""
        previous = np.array([0.0, 0.0, 0.0])
        current = np.array([0.0, 0.0, CA_CA_BOND_LENGTH])
        candidates = cone_candidates(previous, current, angles=[110.0, 120.0], per_angle=8)
        azimuth = np.arctan2(candidates[:, 1], candidates[:, 0]) % (2 * np.pi)
        first, second = azimuth[:8], azimuth[8:]
        assert not np.any(np.isclose(np.sort(first), np.sort(second), atol=1e-6))

    def test_only_the_direction_of_the_previous_bond_matters(self) -> None:
        near = cone_candidates([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], angles=[125.0], per_angle=4)
        far = cone_candidates([-99.0, 0.0, 0.0], [1.0, 0.0, 0.0], angles=[125.0], per_angle=4)
        assert np.allclose(near, far)

    def test_translation_equivariance(self) -> None:
        shift = np.array([100.0, -50.0, 7.0])
        base = cone_candidates([0.0, 0.0, 0.0], [3.8, 0.0, 0.0])
        moved = cone_candidates(shift, shift + np.array([3.8, 0.0, 0.0]))
        assert np.allclose(moved, base + shift)

    def test_output_is_writable_and_independent_of_the_cache(self) -> None:
        """A cached template that leaked into the output would be poisoned by a caller."""
        first = cone_candidates([0, 0, 0], [3.8, 0, 0], angles=[125.0], per_angle=3)
        first[:] = 999.0
        second = cone_candidates([0, 0, 0], [3.8, 0, 0], angles=[125.0], per_angle=3)
        assert not np.any(second == 999.0)

    def test_template_cache_is_used_and_bit_identical(self) -> None:
        before = cone_template_cache_info()
        a = cone_candidates([0, 0, 0], [3.8, 0, 0])
        b = cone_candidates([0, 0, 0], [3.8, 0, 0])
        after = cone_template_cache_info()
        assert np.array_equal(a, b), "identical input must give bit-identical output"
        assert after.hits > before.hits, "the cone template cache is not being hit"

    def test_defaults_share_a_cache_entry_with_the_explicit_equivalent(self) -> None:
        cone_candidates([0, 0, 0], [3.8, 0, 0])
        before = cone_template_cache_info()
        cone_candidates(
            [0, 0, 0],
            [3.8, 0, 0],
            angles=backbone_angle_grid(),
            per_angle=CANDIDATES_PER_ANGLE,
            bond_length=CA_CA_BOND_LENGTH,
        )
        assert cone_template_cache_info().misses == before.misses

    def test_custom_bond_length_is_honoured(self) -> None:
        candidates = cone_candidates([0, 0, 0], [1.0, 0, 0], bond_length=5.0)
        assert np.allclose(np.linalg.norm(candidates - [1.0, 0, 0], axis=1), 5.0)

    def test_coincident_input_raises_rather_than_returning_nan(self) -> None:
        with pytest.raises(GeometryError, match="same point"):
            cone_candidates([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])

    @pytest.mark.parametrize("angles", [[0.0], [181.0], [-5.0], [90.0, 200.0]])
    def test_out_of_range_angles_raise(self, angles: list[float]) -> None:
        with pytest.raises(GeometryError, match="0, 180"):
            cone_candidates([0, 0, 0], [3.8, 0, 0], angles=angles)

    def test_empty_angles_raise(self) -> None:
        with pytest.raises(GeometryError, match="non-empty"):
            cone_candidates([0, 0, 0], [3.8, 0, 0], angles=[])

    def test_nan_angles_raise(self) -> None:
        with pytest.raises(GeometryError, match="non-finite"):
            cone_candidates([0, 0, 0], [3.8, 0, 0], angles=[np.nan])

    def test_bad_per_angle_raises(self) -> None:
        with pytest.raises(GeometryError, match="per_angle"):
            cone_candidates([0, 0, 0], [3.8, 0, 0], per_angle=0)

    @pytest.mark.parametrize("bond_length", [0.0, -3.8, np.nan])
    def test_bad_bond_length_raises(self, bond_length: float) -> None:
        with pytest.raises(GeometryError, match="bond_length"):
            cone_candidates([0, 0, 0], [3.8, 0, 0], bond_length=bond_length)

    def test_nan_input_coordinates_raise(self) -> None:
        with pytest.raises(GeometryError, match="non-finite"):
            cone_candidates([np.nan, 0, 0], [3.8, 0, 0])

    def test_scalar_angle_is_accepted(self) -> None:
        assert cone_candidates([0, 0, 0], [3.8, 0, 0], angles=125.0, per_angle=2).shape == (2, 3)

    def test_taking_the_first_candidate_is_not_a_selection_policy(self) -> None:
        """The ideal-first ORDERING is a bias on angles, not a rule for picking a step.

        Measured here at n=200, with an identical candidate set and clash rule for both
        policies: always taking the first non-clashing candidate gives Re 334-641 A, 3.4x
        to 6.5x the 98.8 A that ``flory_end_to_end(200)`` predicts, because the azimuth
        order within each ring is fixed and correlates successive dihedrals into a near
        rod. Sampling among the non-clashing candidates gives Re 41-162 A, mean ~100 A.

        Both policies pass ``validate_ca_trace``: bonds, angles and clashes are all fine,
        and no local gate can see a globally over-extended chain. The docstring is the
        only place this can be prevented, so it is asserted here.
        """
        n = 200
        predicted = flory_end_to_end(n)
        greedy = [grow_with_selection_policy(n, seed, "first") for seed in range(5)]
        sampled = [grow_with_selection_policy(n, seed, "sampled") for seed in range(5)]

        greedy_re = np.array([end_to_end(trace) for trace in greedy])
        sampled_re = np.array([end_to_end(trace) for trace in sampled])
        assert all(validate_ca_trace(trace).ok for trace in greedy + sampled)
        assert greedy_re.min() > 3.0 * predicted, f"greedy Re {greedy_re} vs {predicted:.1f}"
        assert 0.4 * predicted < sampled_re.mean() < 2.0 * predicted, (
            f"sampled Re {sampled_re} vs {predicted:.1f}"
        )

        doc = cone_candidates.__doc__ or ""
        assert "automatically prefers realistic geometry" not in doc, (
            "the docstring recommends the greedy policy again; measured above, it yields "
            f"Re {greedy_re.min():.0f}-{greedy_re.max():.0f} A against a prediction of "
            f"{predicted:.0f} A"
        )
        assert "over-extended" in doc, "the docstring must warn what greedy selection costs"
        assert "sample" in doc.lower(), "the docstring must say to sample among the candidates"

    def test_always_taking_candidate_zero_collapses_the_chain(self) -> None:
        """The degenerate limit of greedy selection: a 200-residue planar hexagon.

        From a coplanar start, candidate[0] is the same in-plane choice every step, every
        dihedral comes out at exactly 0 degrees and the chain closes into a hexagon of
        Rg 4.11 A with thousands of exactly coincident CAs. Bonds and pseudo-angles are
        ideal throughout, which is precisely why the ordering must not be sold as a
        selection rule.
        """
        n = 200
        coords = np.zeros((n, 3))
        coords[1] = [CA_CA_BOND_LENGTH, 0.0, 0.0]
        for i in range(2, n):
            coords[i] = cone_candidates(coords[i - 2], coords[i - 1])[0]

        assert radius_of_gyration(coords) < 10.0
        report = validate_ca_trace(coords)
        assert not report.ok
        assert "steric_clash" in {violation.kind for violation in report.violations}

    def test_template_cache_cannot_be_poisoned_through_the_base_array(self) -> None:
        """setflags(write=False) on a view leaves ``view.base`` writable.

        Confirmed before the fix: writing zeros through ``tmpl.base`` made every later
        ``cone_candidates`` call for that cache key return all candidates sitting on the
        apex -- zero-length bonds, from a function whose contract is that every candidate
        is exactly one bond away.
        """
        # Whitebox on purpose: the promise being tested is about the cached object, which
        # is not reachable through the public API by design.
        from dodo.geometry.sampling import _cone_template

        template = _cone_template((125.0,), 4, CA_CA_BOND_LENGTH)
        assert not template.flags.writeable
        base = template.base
        assert base is None or not base.flags.writeable, "the cache is mutable through .base"

        apex = np.array([CA_CA_BOND_LENGTH, 0.0, 0.0])
        candidates = cone_candidates([0.0, 0.0, 0.0], apex, angles=[125.0], per_angle=4)
        assert np.allclose(np.linalg.norm(candidates - apex, axis=1), CA_CA_BOND_LENGTH)

    @pytest.mark.parametrize("per_angle", [2.9, 0.5, np.float64(3.5)])
    def test_non_integer_per_angle_raises_rather_than_truncating(self, per_angle: float) -> None:
        """``per_angle=2.9`` used to silently return 2 candidates.

        A caller that computes its candidate budget arithmetically then searches a
        narrower cone than it believes it asked for, and a build that fails to place a
        residue gives no hint why.
        """
        with pytest.raises(GeometryError, match="integer"):
            cone_candidates(
                [0, 0, 0],
                [3.8, 0, 0],
                angles=[125.0],
                per_angle=per_angle,  # type: ignore[arg-type]
            )

    def test_integer_like_per_angle_is_accepted(self) -> None:
        shape = cone_candidates([0, 0, 0], [3.8, 0, 0], angles=[125.0], per_angle=np.int64(3)).shape
        assert shape == (3, 3)

    def test_a_180_degree_ring_is_a_single_repeated_point(self) -> None:
        """Documented degeneracy: at 180 degrees the cone has zero radius.

        The ring collapses onto the straight continuation of the previous bond, so all
        ``per_angle`` rows coincide. The shape contract is kept (the caller gets the rows
        it asked for) and the docstring says so, rather than the caller discovering that
        10 candidates offer 1 choice.
        """
        apex = np.array([CA_CA_BOND_LENGTH, 0.0, 0.0])
        candidates = cone_candidates([0.0, 0.0, 0.0], apex, angles=[180.0], per_angle=10)
        assert candidates.shape == (10, 3)
        assert np.allclose(candidates, candidates[0])
        assert np.allclose(candidates[0], [2.0 * CA_CA_BOND_LENGTH, 0.0, 0.0])
        assert "180" in (cone_candidates.__doc__ or "")

    def test_is_fast_enough_to_call_per_residue(self) -> None:
        """Measured at 44 us per warm call, against 2.82 ms for the pre-rewrite version.

        The bound is loose because this runs on shared CI hardware; it is here to catch a
        regression of the kind the pre-rewrite code had (a Python loop over 71 angles,
        each transforming its points one at a time), not to police microseconds.
        """
        previous, current = np.zeros(3), np.array([3.8, 0.0, 0.0])
        cone_candidates(previous, current)  # warm the cache
        start = time.perf_counter()
        for _ in range(200):
            cone_candidates(previous, current)
        per_call_ms = (time.perf_counter() - start) / 200 * 1e3
        assert per_call_ms < 1.0, f"{per_call_ms:.3f} ms per call"


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------


def straight_trace(n: int, spacing: float = CA_CA_BOND_LENGTH) -> np.ndarray:
    """Build a perfectly straight CA trace: ideal bonds, 180 degree angles."""
    coords = np.zeros((n, 3))
    coords[:, 0] = spacing * np.arange(n)
    return coords


def realistic_trace(n: int, seed: int = 0) -> np.ndarray:
    """Grow a SELF-AVOIDING trace with cone_candidates.

    Self-avoidance is not optional here, and getting that wrong is instructive: this helper
    originally chose uniformly among all cone candidates with no clash check, and the traces
    it produced were asserted to be valid. They were not. Cone candidates guarantee correct
    bond lengths and pseudo-angles, but those are purely LOCAL properties -- a chain can
    satisfy both perfectly while folding back through itself. Measured on the original
    helper: non-sequential CA pairs down to 0.753 A at n=150 and 0.450 A at n=600.

    The tests only stopped agreeing once ``validate_ca_trace`` learned to check non-bonded
    contacts. Any generator that omits this step produces geometry no viewer should be asked
    to render.
    """
    rng = np.random.default_rng(seed)
    coords = np.zeros((n, 3))
    coords[1] = CA_CA_BOND_LENGTH * random_unit_vectors(1, rng)[0]
    for i in range(2, n):
        candidates = cone_candidates(coords[i - 2], coords[i - 1])
        order = rng.permutation(candidates.shape[0])
        placed = False
        for index in order:
            candidate = candidates[index]
            # Compare against every residue up to i-2; i-1 is the bonded neighbour and i-2
            # is angle-constrained, so neither is a clash by definition.
            if i > 2:
                separations = np.linalg.norm(coords[: i - 2] - candidate, axis=1)
                if float(separations.min()) < CA_CLASH_DISTANCE:
                    continue
            coords[i] = candidate
            placed = True
            break
        if not placed:
            # Boxed in. A real builder backtracks; for a test fixture, restart with a
            # different seed rather than silently emitting a clashing trace.
            return realistic_trace(n, seed=seed + 10_000)
    return coords


class TestConeCandidatesBatch:
    """Lock the batched cone builder to the scalar cone_candidates, row by row.

    One cone per row, equal to the scalar ``cone_candidates`` on each row bit for bit. This
    is what lets the walk engine build every live conformer's candidates in one call.
    """

    def test_matches_scalar_per_row_bit_for_bit(self) -> None:
        rng = np.random.default_rng(2024)
        grid = backbone_angle_grid()
        previous = rng.normal(size=(600, 3)) * 10.0
        # Apex one generic step from `previous`, plus two rows whose axis lies exactly along
        # +z / -z -- the degenerate directions the shared rotation has to get right.
        current = previous + rng.normal(size=(600, 3))
        current[0] = previous[0] + np.array([0.0, 0.0, CA_CA_BOND_LENGTH])
        current[1] = previous[1] + np.array([0.0, 0.0, -CA_CA_BOND_LENGTH])
        batch = cone_candidates_batch(previous, current, angles=grid, per_angle=5)
        assert batch.shape == (600, grid.size * 5, 3)
        for i in range(previous.shape[0]):
            scalar = cone_candidates(previous[i], current[i], angles=grid, per_angle=5)
            assert np.array_equal(batch[i], scalar)

    def test_shares_the_template_cache_with_the_scalar_path(self) -> None:
        grid = backbone_angle_grid()
        previous = np.zeros((3, 3))
        current = np.tile([CA_CA_BOND_LENGTH, 0.0, 0.0], (3, 1))
        cone_candidates(previous[0], current[0], angles=grid, per_angle=5)  # warm the cache
        before = cone_template_cache_info().misses
        cone_candidates_batch(previous, current, angles=grid, per_angle=5)
        # A batched call for an already-seen (angles, per_angle, bond) resolves to the same
        # cached template rather than rebuilding it.
        assert cone_template_cache_info().misses == before

    def test_mismatched_shapes_raise(self) -> None:
        with pytest.raises(GeometryError, match=r"shape \(n, 3\)"):
            cone_candidates_batch(np.zeros((3, 3)), np.zeros((2, 3)), angles=[120.0], per_angle=2)


class TestCaBondLengths:
    def test_straight_chain(self) -> None:
        lengths = ca_bond_lengths(straight_trace(5))
        assert lengths.shape == (4,)
        assert np.allclose(lengths, CA_CA_BOND_LENGTH)

    def test_two_residues_is_enough(self) -> None:
        assert ca_bond_lengths(np.array([[0.0, 0, 0], [1.0, 0, 0]])) == pytest.approx([1.0])

    def test_one_residue_raises(self) -> None:
        with pytest.raises(GeometryError, match="at least 2"):
            ca_bond_lengths(np.zeros((1, 3)))

    def test_bad_shape_raises(self) -> None:
        with pytest.raises(GeometryError, match=r"\(n, 3\)"):
            ca_bond_lengths(np.zeros((5, 2)))

    def test_nan_raises(self) -> None:
        coords = straight_trace(4)
        coords[2, 1] = np.nan
        with pytest.raises(GeometryError, match="non-finite"):
            ca_bond_lengths(coords)


class TestCaPseudoAngles:
    def test_straight_chain_is_exactly_180_not_nan(self) -> None:
        """arccos(-1 - 1e-16) is NaN; the clip in the implementation is load-bearing."""
        angles = ca_pseudo_angles(straight_trace(6))
        assert angles.shape == (4,)
        assert not np.any(np.isnan(angles))
        assert np.allclose(angles, 180.0)

    def test_right_angle(self) -> None:
        coords = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        assert ca_pseudo_angles(coords) == pytest.approx([90.0])

    def test_known_geometry(self) -> None:
        coords = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [-0.5, math.sqrt(3) / 2, 0.0]])
        assert ca_pseudo_angles(coords) == pytest.approx([120.0])

    def test_vertex_indexing(self) -> None:
        """Entry i is the angle at residue i+1, the middle of the triple."""
        coords = straight_trace(5)
        coords[3] = coords[2] + [0.0, CA_CA_BOND_LENGTH, 0.0]
        angles = ca_pseudo_angles(coords)
        assert angles[0] == pytest.approx(180.0)
        assert angles[1] == pytest.approx(90.0)  # vertex is residue 2

    def test_too_few_residues_raises(self) -> None:
        with pytest.raises(GeometryError, match="at least 3"):
            ca_pseudo_angles(straight_trace(2))

    def test_coincident_neighbour_raises_instead_of_returning_nan(self) -> None:
        coords = straight_trace(4)
        coords[2] = coords[1]
        with pytest.raises(GeometryError, match="undefined"):
            ca_pseudo_angles(coords)

    def test_matches_an_independent_computation(self) -> None:
        coords = realistic_trace(30, seed=4)
        angles = ca_pseudo_angles(coords)
        expected = [measured_angle(coords[i], coords[i + 1], coords[i + 2]) for i in range(28)]
        assert np.allclose(angles, expected)


class TestRgAndEndToEnd:
    def test_end_to_end_of_a_straight_chain_is_the_contour_length(self) -> None:
        assert end_to_end(straight_trace(11)) == pytest.approx(10 * CA_CA_BOND_LENGTH)

    def test_end_to_end_ignores_the_path(self) -> None:
        coords = np.array([[0.0, 0, 0], [5.0, 5.0, 0], [3.0, 4.0, 0.0]])
        assert end_to_end(coords) == pytest.approx(5.0)

    def test_rg_of_two_points(self) -> None:
        coords = np.array([[-1.0, 0, 0], [1.0, 0, 0]])
        assert radius_of_gyration(coords) == pytest.approx(1.0)

    def test_rg_is_translation_and_rotation_invariant(self) -> None:
        coords = realistic_trace(40, seed=11)
        matrix = rotation_between_vectors([1.0, 0, 0], [0.2, -0.5, 0.9])
        moved = apply(coords, matrix) + np.array([100.0, 100.0, 100.0])
        assert radius_of_gyration(moved) == pytest.approx(radius_of_gyration(coords))

    def test_agrees_with_the_structure_core(self) -> None:
        """Phase 1 already computes these; the two must not drift apart."""
        coords = realistic_trace(20, seed=2)
        structure = Structure.from_atom_records(
            xyz=coords,
            atom_name=["CA"] * 20,
            element=["C"] * 20,
            residue_name=["ALA"] * 20,
            residue_number=list(range(1, 21)),
            chain_id=["A"] * 20,
        )
        assert radius_of_gyration(coords) == pytest.approx(structure.radius_of_gyration())
        assert end_to_end(coords) == pytest.approx(structure.end_to_end_distance())

    @pytest.mark.parametrize("function", [radius_of_gyration, end_to_end])
    def test_single_point_raises(self, function: object) -> None:
        with pytest.raises(GeometryError, match="at least 2"):
            function(np.zeros((1, 3)))  # type: ignore[operator]

    @pytest.mark.parametrize("function", [radius_of_gyration, end_to_end])
    def test_nan_raises(self, function: object) -> None:
        with pytest.raises(GeometryError, match="non-finite"):
            function(np.full((3, 3), np.nan))  # type: ignore[operator]


class TestValidateCaTrace:
    def test_a_cone_grown_trace_is_valid(self) -> None:
        report = validate_ca_trace(realistic_trace(60, seed=17))
        assert report.ok, report.describe()
        assert bool(report) is True

    @pytest.mark.parametrize("seed", range(5))
    def test_cone_candidates_land_exactly_on_the_window_edges_without_flagging(
        self, seed: int
    ) -> None:
        """91 and 161 degrees are in the grid, so a strict window check would flag them."""
        previous, current = np.zeros(3), np.array([0.0, 0.0, CA_CA_BOND_LENGTH])
        for angle in (BACKBONE_ANGLE_MIN, BACKBONE_ANGLE_MAX):
            candidate = cone_candidates(previous, current, angles=[angle], per_angle=4)[seed % 4]
            report = validate_ca_trace(np.vstack([previous, current, candidate]))
            assert report.ok, report.describe()

    def test_straight_trace_is_rejected_for_its_angles(self) -> None:
        report = validate_ca_trace(straight_trace(6))
        assert not report.ok
        assert not report.of_kind("bond_length")
        assert len(report.of_kind("pseudo_angle")) == 4
        assert all(v.value == pytest.approx(180.0) for v in report.of_kind("pseudo_angle"))

    def test_stretched_bond_is_flagged_with_its_residue_index(self) -> None:
        coords = realistic_trace(10, seed=1)
        coords[5:] += [12.0, 0.0, 0.0]  # blows out the 4-5 bond only
        report = validate_ca_trace(coords)
        bonds = report.of_kind("bond_length")
        assert [v.residue_index for v in bonds] == [4]
        assert bonds[0].value > CA_CA_BOND_LENGTH
        assert "residues 4-5" in bonds[0].message
        assert bonds[0].low is not None and bonds[0].high is not None

    def test_sharp_angle_is_flagged_at_the_vertex(self) -> None:
        # A hairpin: residue 2 folds straight back on residue 0.
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [CA_CA_BOND_LENGTH, 0.0, 0.0],
                [CA_CA_BOND_LENGTH / 2, 1.0, 0.0],
            ]
        )
        coords[2] = coords[1] + CA_CA_BOND_LENGTH * np.array([-0.9, np.sqrt(1 - 0.81), 0.0])
        report = validate_ca_trace(coords)
        angles = report.of_kind("pseudo_angle")
        assert [v.residue_index for v in angles] == [1]
        assert angles[0].value < BACKBONE_ANGLE_MIN

    def test_flags_the_47_degree_angle_the_old_builder_produced(self) -> None:
        """The pre-rewrite IDR builder had no angle constraint and measured 47 to 178."""
        coords = np.array([[CA_CA_BOND_LENGTH, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        coords[2] = CA_CA_BOND_LENGTH * np.array(
            [math.cos(math.radians(47.0)), math.sin(math.radians(47.0)), 0.0]
        )
        report = validate_ca_trace(coords)
        assert not report.ok
        assert report.of_kind("pseudo_angle")[0].value == pytest.approx(47.0)

    def test_flags_a_trace_of_zero_rows(self) -> None:
        """The old builder returned arrays of exact (0, 0, 0) rows on total failure."""
        coords = realistic_trace(8, seed=5)
        coords[4:] = 0.0
        report = validate_ca_trace(coords)
        assert not report.ok
        assert report.of_kind("bond_length")
        assert report.of_kind("undefined_angle")

    def test_residue_offset_shifts_reported_indices(self) -> None:
        report = validate_ca_trace(straight_trace(4), residue_offset=100)
        assert [v.residue_index for v in report.violations] == [101, 102]
        assert "residue 101" in report.violations[0].message

    def test_nan_coordinate_is_reported_not_raised(self) -> None:
        coords = realistic_trace(8, seed=3)
        coords[4] = np.nan
        report = validate_ca_trace(coords)
        assert not report.ok
        non_finite = report.of_kind("non_finite")
        assert [v.residue_index for v in non_finite] == [4]
        # The NaN must not leak out as spurious bond or angle violations with NaN values.
        assert not any(np.isnan(v.value) for v in report.of_kind("bond_length"))
        assert not any(np.isnan(v.value) for v in report.of_kind("pseudo_angle"))

    def test_coincident_residues_report_a_bond_and_an_undefined_angle(self) -> None:
        coords = realistic_trace(6, seed=6)
        coords[3] = coords[2]
        report = validate_ca_trace(coords)
        # The zero-length 2-3 bond makes both bonds around it out of range, and leaves the
        # angles at the two vertices that use that bond undefined -- but not the one at
        # residue 4, whose own two bonds are fine.
        assert [v.residue_index for v in report.of_kind("bond_length")] == [2, 3]
        assert [v.residue_index for v in report.of_kind("undefined_angle")] == [2, 3]
        assert np.isnan(report.pseudo_angles[[1, 2]]).all()
        assert not np.isnan(report.pseudo_angles[3])

    def test_two_residues_reports_bonds_and_no_angles(self) -> None:
        report = validate_ca_trace(np.array([[0.0, 0, 0], [CA_CA_BOND_LENGTH, 0, 0]]))
        assert report.ok
        assert report.bond_lengths.shape == (1,)
        assert report.pseudo_angles.shape == (0,)

    def test_one_residue_raises(self) -> None:
        with pytest.raises(GeometryError, match="at least 2"):
            validate_ca_trace(np.zeros((1, 3)))

    def test_bad_shape_raises(self) -> None:
        with pytest.raises(GeometryError, match=r"\(n, 3\)"):
            validate_ca_trace(np.zeros(3))

    def test_custom_angle_window(self) -> None:
        straight = straight_trace(5)
        assert not validate_ca_trace(straight).ok
        assert validate_ca_trace(straight, angle_window=(90.0, 180.0)).ok

    def test_custom_bond_tolerance(self) -> None:
        stretched = straight_trace(5, spacing=4.5)
        wide_angles = (90.0, 180.0)
        assert validate_ca_trace(stretched, angle_window=wide_angles).of_kind("bond_length")
        assert not validate_ca_trace(
            stretched, bond_tolerance=1.0, angle_window=wide_angles
        ).of_kind("bond_length")

    def test_observed_range_window_accepts_a_real_structure_angle(self) -> None:
        """An 80 degree angle is outside the generation window but is observed in nature."""
        coords = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        coords *= CA_CA_BOND_LENGTH
        coords[2] = CA_CA_BOND_LENGTH * np.array(
            [math.cos(math.radians(80.0)), math.sin(math.radians(80.0)), 0.0]
        )
        assert validate_ca_trace(coords).of_kind("pseudo_angle")
        assert not validate_ca_trace(coords, angle_window=(75.0, 179.0)).of_kind("pseudo_angle")

    def test_inverted_window_raises(self) -> None:
        with pytest.raises(GeometryError, match="exceeds"):
            validate_ca_trace(straight_trace(4), angle_window=(160.0, 90.0))

    def test_negative_tolerance_raises(self) -> None:
        with pytest.raises(GeometryError, match="non-negative"):
            validate_ca_trace(straight_trace(4), bond_tolerance=-0.1)

    def test_report_exposes_the_measurements(self) -> None:
        coords = realistic_trace(12, seed=8)
        report = validate_ca_trace(coords)
        assert report.n_residues == 12
        assert np.allclose(report.bond_lengths, ca_bond_lengths(coords))
        assert np.allclose(report.pseudo_angles, ca_pseudo_angles(coords))
        assert report.bond_range[0] < CA_CA_BOND_LENGTH < report.bond_range[1]
        assert report.angle_window == (BACKBONE_ANGLE_MIN, BACKBONE_ANGLE_MAX)

    def test_describe_names_specifics(self) -> None:
        text = validate_ca_trace(straight_trace(4)).describe()
        assert "180.00 deg" in text
        assert "residue 1" in text
        assert f"{BACKBONE_ANGLE_MIN:.1f}-{BACKBONE_ANGLE_MAX:.1f}" in text

    def test_describe_truncates(self) -> None:
        text = validate_ca_trace(straight_trace(40)).describe(max_violations=3)
        assert text.count("\n  - ") == 4  # 3 violations plus the "and N more" line
        assert "and 35 more" in text

    def test_summary_reports_ranges(self) -> None:
        summary = validate_ca_trace(realistic_trace(20, seed=9)).summary()
        assert "20 residues" in summary
        assert "valid" in summary

    def test_raise_if_invalid(self) -> None:
        validate_ca_trace(realistic_trace(10, seed=10)).raise_if_invalid()
        with pytest.raises(GeometryError, match="Invalid CA trace"):
            validate_ca_trace(straight_trace(5)).raise_if_invalid()

    def test_violating_residues(self) -> None:
        report = validate_ca_trace(straight_trace(5))
        assert report.violating_residues == (1, 2, 3)

    def test_report_can_be_compared_and_hashed_without_ambiguity(self) -> None:
        """A dataclass __eq__ over numpy fields would raise "truth value is ambiguous"."""
        a = validate_ca_trace(straight_trace(4))
        b = validate_ca_trace(straight_trace(4))
        assert a != b  # identity comparison, and crucially does not raise
        assert isinstance(a, TraceReport)
        assert len({a, b}) == 2


# ---------------------------------------------------------------------------
# Cross-module: the kernel must be able to grow a physically valid backbone
# ---------------------------------------------------------------------------


class TestKernelIntegration:
    @pytest.mark.parametrize("seed", range(8))
    def test_random_walk_over_cone_candidates_is_always_valid(self, seed: int) -> None:
        """The acceptance criterion for the whole kernel.

        The pre-rewrite IDR builder had no angle constraint and produced measured
        pseudo-angles from 47 to 178 degrees. A trace grown from these primitives cannot
        do that, whatever the seed.
        """
        report = validate_ca_trace(realistic_trace(150, seed=seed))
        assert report.ok, report.describe()

    def test_the_walk_is_reproducible(self) -> None:
        assert np.array_equal(realistic_trace(50, seed=42), realistic_trace(50, seed=42))

    def test_frames_and_rotations_compose(self) -> None:
        """A candidate expressed in the local frame lands where the world rotation says."""
        previous = np.array([1.0, 2.0, 3.0])
        current = previous + CA_CA_BOND_LENGTH * np.array([0.0, 0.0, -1.0])
        frame = align_frame(previous, current)
        candidate = cone_candidates(previous, current, angles=[125.0], per_angle=1)[0]
        local = apply(candidate - current, frame.T)
        # In the frame whose z is the chain axis, the axial component is bond * cos(55).
        assert local[2] == pytest.approx(CA_CA_BOND_LENGTH * math.cos(math.radians(55.0)))
