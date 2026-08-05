"""Tests for the batched free-growth generator (Phase 3, slice 1)."""

from __future__ import annotations

import numpy as np
import pytest

from dodo.constants import (
    BACKBONE_ANGLE_MAX,
    BACKBONE_ANGLE_MIN,
    CA_CA_BOND_LENGTH,
    CA_CLASH_DISTANCE,
    flory_end_to_end,
)
from dodo.engines.batch_walk import (
    clears_obstacles_mask,
    close_chain_ends,
    end_to_end_distances,
    generate_free_batch,
    generate_interior_batch,
    grow_free_batch,
    radii_of_gyration,
    self_avoiding_mask,
    steer_to_target,
)
from dodo.exceptions import EngineError


def _pseudo_angles(coords: np.ndarray) -> np.ndarray:
    u = coords[:, :-2, :] - coords[:, 1:-1, :]
    v = coords[:, 2:, :] - coords[:, 1:-1, :]
    u /= np.linalg.norm(u, axis=-1, keepdims=True)
    v /= np.linalg.norm(v, axis=-1, keepdims=True)
    return np.degrees(np.arccos(np.clip((u * v).sum(-1), -1.0, 1.0)))


class TestGrowthGeometry:
    @pytest.mark.parametrize("n", [2, 3, 10, 50, 200])
    def test_bond_lengths_are_exact(self, n: int) -> None:
        coords = grow_free_batch(500, n, np.random.default_rng(n))
        bonds = np.linalg.norm(np.diff(coords, axis=1), axis=-1)
        assert np.max(np.abs(bonds - CA_CA_BOND_LENGTH)) < 1e-9

    @pytest.mark.parametrize("n", [3, 10, 50, 200])
    def test_pseudo_angles_stay_inside_the_window(self, n: int) -> None:
        angles = _pseudo_angles(grow_free_batch(500, n, np.random.default_rng(n)))
        # A hair of floating-point slack on the truncation bounds.
        assert angles.min() >= BACKBONE_ANGLE_MIN - 1e-6
        assert angles.max() <= BACKBONE_ANGLE_MAX + 1e-6

    def test_shape_and_first_residue_at_origin(self) -> None:
        coords = grow_free_batch(7, 12, np.random.default_rng(0))
        assert coords.shape == (7, 12, 3)
        assert np.all(coords[:, 0, :] == 0.0)

    def test_deterministic_under_seed(self) -> None:
        a = grow_free_batch(64, 40, np.random.default_rng(123))
        b = grow_free_batch(64, 40, np.random.default_rng(123))
        assert np.array_equal(a, b)

    def test_a_different_seed_differs(self) -> None:
        a = grow_free_batch(64, 40, np.random.default_rng(1))
        b = grow_free_batch(64, 40, np.random.default_rng(2))
        assert not np.array_equal(a, b)

    def test_two_residue_chain_is_one_bond(self) -> None:
        coords = grow_free_batch(100, 2, np.random.default_rng(0))
        bonds = np.linalg.norm(coords[:, 1, :] - coords[:, 0, :], axis=-1)
        assert np.allclose(bonds, CA_CA_BOND_LENGTH)

    @pytest.mark.parametrize(("chains", "residues"), [(0, 10), (5, 1), (-1, 5)])
    def test_invalid_inputs_are_refused(self, chains: int, residues: int) -> None:
        with pytest.raises(EngineError):
            grow_free_batch(chains, residues, np.random.default_rng(0))


class TestSelfAvoidingMask:
    def _straight(self, n: int) -> np.ndarray:
        return np.stack([[[i * CA_CA_BOND_LENGTH, 0.0, 0.0] for i in range(n)]]).astype(float)

    def test_an_extended_trace_is_kept(self) -> None:
        assert bool(self_avoiding_mask(self._straight(20))[0])

    def test_a_distant_non_local_collision_is_rejected(self) -> None:
        coords = self._straight(20)
        # Residue 8 dropped onto residue 0: 8 residues apart, so it is a scored, real clash.
        coords[0, 8] = coords[0, 0] + np.array([0.4, 0.0, 0.0])
        assert not bool(self_avoiding_mask(coords)[0])

    def test_a_sequence_local_close_pair_is_excluded(self) -> None:
        coords = self._straight(20)
        # Residue 2 placed 1 A from residue 0: within the exclusion window, so NOT a clash.
        coords[0, 2] = coords[0, 0] + np.array([1.0, 0.0, 0.0])
        assert bool(self_avoiding_mask(coords)[0])

    def test_chunking_matches_a_single_block(self) -> None:
        coords = grow_free_batch(300, 60, np.random.default_rng(5))
        whole = self_avoiding_mask(coords, chunk_pairs=10**9)
        chunked = self_avoiding_mask(coords, chunk_pairs=6000)
        assert np.array_equal(whole, chunked)

    def test_kept_traces_really_are_clear(self) -> None:
        coords = grow_free_batch(400, 40, np.random.default_rng(9))
        keep = self_avoiding_mask(coords)
        survivors = coords[keep]
        n = survivors.shape[1]
        separation = np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
        for trace in survivors[:50]:
            d = np.linalg.norm(trace[:, None, :] - trace[None, :, :], axis=-1)
            assert d[separation > 2].min() >= CA_CLASH_DISTANCE - 1e-9


class TestMetrics:
    def test_end_to_end_of_a_straight_trace(self) -> None:
        coords = np.stack([[[i * CA_CA_BOND_LENGTH, 0.0, 0.0] for i in range(11)]]).astype(float)
        assert end_to_end_distances(coords)[0] == pytest.approx(10 * CA_CA_BOND_LENGTH)

    def test_radius_of_gyration_is_positive_and_finite(self) -> None:
        coords = grow_free_batch(50, 30, np.random.default_rng(0))
        rg = radii_of_gyration(coords)
        assert rg.shape == (50,)
        assert np.all(np.isfinite(rg)) and np.all(rg > 0)


class TestBiasedGrowth:
    def test_unbiased_default_is_bit_identical(self) -> None:
        a = grow_free_batch(200, 50, np.random.default_rng(0))
        b = grow_free_batch(200, 50, np.random.default_rng(0), bias_directions=None)
        assert np.array_equal(a, b)

    def test_biased_growth_preserves_bonds_and_angles(self) -> None:
        bias = np.tile([-1.0, 0.0, 0.0], (500, 1))
        coords = grow_free_batch(
            500, 50, np.random.default_rng(1), bias_directions=bias, bias_residues=8
        )
        bonds = np.linalg.norm(np.diff(coords, axis=1), axis=-1)
        assert np.max(np.abs(bonds - CA_CA_BOND_LENGTH)) < 1e-9
        angles = _pseudo_angles(coords)
        assert angles.min() >= BACKBONE_ANGLE_MIN - 1e-6
        assert angles.max() <= BACKBONE_ANGLE_MAX + 1e-6

    def test_bias_skews_the_stub_toward_the_direction(self) -> None:
        bias = np.tile([1.0, 0.0, 0.0], (2000, 1))
        biased = grow_free_batch(
            2000, 40, np.random.default_rng(2), bias_directions=bias, bias_residues=10
        )
        free = grow_free_batch(2000, 40, np.random.default_rng(2))
        assert biased[:, 10, 0].mean() > free[:, 10, 0].mean() + 5.0

    def test_biased_conformers_remain_diverse(self) -> None:
        bias = np.tile([-1.0, 0.0, 0.0], (1000, 1))
        coords = grow_free_batch(
            1000, 40, np.random.default_rng(3), bias_directions=bias, bias_residues=10
        )
        ends = coords[:, -1, :] - coords[:, 0, :]
        ends /= np.linalg.norm(ends, axis=1, keepdims=True)
        cosines = ends @ ends.T
        off_diagonal = cosines[~np.eye(cosines.shape[0], dtype=bool)].mean()
        assert off_diagonal < 0.7  # not collapsed onto one direction


class TestObstacleClearance:
    def _straight(self, n: int) -> np.ndarray:
        return np.stack([[[i * CA_CA_BOND_LENGTH, 0.0, 0.0] for i in range(n)]]).astype(float)

    def test_empty_obstacles_keeps_everything(self) -> None:
        coords = grow_free_batch(10, 20, np.random.default_rng(0))
        assert clears_obstacles_mask(coords, np.empty((0, 3))).all()

    def test_an_obstacle_on_the_chain_is_rejected(self) -> None:
        coords = self._straight(10)
        obstacle = coords[0, 5][None, :] + np.array([[0.3, 0.0, 0.0]])
        assert not bool(clears_obstacles_mask(coords, obstacle)[0])

    def test_a_distant_obstacle_is_cleared(self) -> None:
        coords = self._straight(10)
        assert bool(clears_obstacles_mask(coords, np.array([[500.0, 500.0, 500.0]]))[0])

    def test_chunking_matches_a_single_block(self) -> None:
        coords = grow_free_batch(200, 40, np.random.default_rng(3))
        obstacles = np.random.default_rng(4).normal(scale=30.0, size=(60, 3))
        whole = clears_obstacles_mask(coords, obstacles, chunk_pairs=10**9)
        chunked = clears_obstacles_mask(coords, obstacles, chunk_pairs=5000)
        assert np.array_equal(whole, chunked)


class TestChainClosure:
    def test_restores_bonds_and_pins_ends_when_feasible(self) -> None:
        coords = grow_free_batch(200, 30, np.random.default_rng(0))
        start = coords[:, 0].copy()
        # A reachable target: within the 29-bond contour length.
        end = start + np.array([40.0, 0.0, 0.0])
        closed = close_chain_ends(coords, start, end, iterations=40)
        bonds = np.linalg.norm(np.diff(closed, axis=1), axis=-1)
        assert np.max(np.abs(bonds - CA_CA_BOND_LENGTH)) < 1e-2
        assert np.allclose(closed[:, 0], start)
        assert np.allclose(closed[:, -1], end, atol=1e-2)

    def test_infeasible_span_leaves_the_end_unreached(self) -> None:
        coords = grow_free_batch(20, 10, np.random.default_rng(1))
        start = coords[:, 0].copy()
        # 9 bonds reach at most 9*3.81 ~ 34 A; ask for 200 A.
        end = start + np.array([200.0, 0.0, 0.0])
        closed = close_chain_ends(coords, start, end, iterations=40)
        bonds = np.linalg.norm(np.diff(closed, axis=1), axis=-1)
        # FABRIK keeps every bond exact; infeasibility shows as the end falling far short.
        assert np.max(np.abs(bonds - CA_CA_BOND_LENGTH)) < 1e-6
        assert np.min(np.linalg.norm(closed[:, -1] - end, axis=-1)) > 100.0


class TestInteriorClosure:
    N_ANCHOR = np.array([-21.0, 0.0, 0.0])
    C_ANCHOR = np.array([21.0, 0.0, 0.0])
    N_AWAY = np.array([1.0, 0.0, 0.0])
    C_AWAY = np.array([-1.0, 0.0, 0.0])

    def _generate(self, seed: int = 0, n: int = 40, stub: int = 8) -> np.ndarray:
        return generate_interior_batch(
            300,
            n,
            self.N_ANCHOR,
            self.C_ANCHOR,
            np.random.default_rng(seed),
            n_away=self.N_AWAY,
            c_away=self.C_AWAY,
            stub=stub,
        )

    def test_endpoints_are_one_bond_from_each_anchor(self) -> None:
        coords = self._generate()
        first_gap = np.linalg.norm(coords[:, 0] - self.N_ANCHOR, axis=-1)
        last_gap = np.linalg.norm(coords[:, -1] - self.C_ANCHOR, axis=-1)
        assert np.allclose(first_gap, CA_CA_BOND_LENGTH, atol=1e-6)
        assert np.allclose(last_gap, CA_CA_BOND_LENGTH, atol=1e-6)

    def test_all_bonds_are_close_to_ideal(self) -> None:
        coords = self._generate()
        bonds = np.linalg.norm(np.diff(coords, axis=1), axis=-1)
        assert np.max(np.abs(bonds - CA_CA_BOND_LENGTH)) < 1e-2

    def test_shape_and_determinism(self) -> None:
        a = self._generate(seed=7)
        b = self._generate(seed=7)
        assert a.shape == (300, 40, 3)
        assert np.array_equal(a, b)

    def test_too_short_for_two_stubs_is_refused(self) -> None:
        with pytest.raises(EngineError):
            generate_interior_batch(
                10,
                12,
                self.N_ANCHOR,
                self.C_ANCHOR,
                np.random.default_rng(0),
                n_away=self.N_AWAY,
                c_away=self.C_AWAY,
                stub=8,
            )


class TestSteering:
    def test_hits_the_target_mean_without_collapsing_the_spread(self) -> None:
        n = 50
        target = flory_end_to_end(n)
        res = generate_free_batch(300, n, target, np.random.default_rng(0), oversample=32)
        assert res.coords.shape == (300, n, 3)
        achieved = res.achieved_end_to_end
        assert abs(float(achieved.mean()) - target) / target < 0.10
        # A physical spread, not a spike pinned to the mean.
        assert 0.15 < float(achieved.std()) / float(achieved.mean()) < 0.60

    def test_compact_and_expanded_modes_are_reachable(self) -> None:
        n = 50
        base = flory_end_to_end(n)
        for mult in (0.7, 1.3):
            res = generate_free_batch(200, n, mult * base, np.random.default_rng(1), oversample=48)
            ratio = float(res.achieved_end_to_end.mean()) / (mult * base)
            assert 0.9 < ratio < 1.1

    def test_deterministic_under_seed(self) -> None:
        target = flory_end_to_end(40)
        a = generate_free_batch(50, 40, target, np.random.default_rng(9), oversample=32)
        b = generate_free_batch(50, 40, target, np.random.default_rng(9), oversample=32)
        assert np.array_equal(a.coords, b.coords)

    def test_returned_conformers_are_self_avoiding_with_exact_bonds(self) -> None:
        res = generate_free_batch(
            100, 40, flory_end_to_end(40), np.random.default_rng(5), oversample=32
        )
        assert self_avoiding_mask(res.coords).all()
        bonds = np.linalg.norm(np.diff(res.coords, axis=1), axis=-1)
        assert np.max(np.abs(bonds - CA_CA_BOND_LENGTH)) < 1e-9

    def test_diagnostics_are_in_range(self) -> None:
        res = generate_free_batch(
            100, 50, flory_end_to_end(50), np.random.default_rng(1), oversample=32
        )
        assert res.n_survived <= res.n_grown
        assert 0.0 <= res.reachable_fraction <= 1.0
        assert 0.0 <= res.reuse_fraction <= 1.0

    def test_steering_an_empty_pool_raises(self) -> None:
        with pytest.raises(EngineError):
            steer_to_target(
                np.empty((0, 10, 3)), 30.0, 5, np.random.default_rng(0), low=3.81, high=35.0
            )

    def test_zero_conformers_is_refused(self) -> None:
        with pytest.raises(EngineError):
            generate_free_batch(0, 20, 30.0, np.random.default_rng(0))
