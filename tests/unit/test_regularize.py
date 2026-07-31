"""Tests for CA bond-length regularization.

This module exists because STARLING is a diffusion model: it denoises coordinates toward a
learned distribution and nothing in that process enforces a hard geometric constraint, so its
CA-CA distances scatter around the true virtual bond length instead of sitting on it. The
correction has to happen after generation.

The tests that matter most are the conformation-preservation ones. Making bonds exact is easy
-- a naive walk-and-replace does it in one pass. Making them exact *without changing the
conformation* is the whole problem, because the conformation is the thing a generative model
was used for in the first place.
"""

from __future__ import annotations

import numpy as np
import pytest

from dodo.constants import CA_CA_BOND_LENGTH
from dodo.exceptions import GeometryError
from dodo.geometry.metrics import ca_bond_lengths, ca_pseudo_angles, radius_of_gyration
from dodo.geometry.regularize import needs_regularization, regularize_ca_trace


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
