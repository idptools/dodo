"""Tests for the batched conformation engine + terminal generator (Phase 3, slice 4)."""

from __future__ import annotations

import numpy as np

from dodo.constants import CA_CA_BOND_LENGTH
from dodo.engines.base import ConformationEngine, IDRRequest, IDRResult
from dodo.engines.batch_engine import BatchWalkEngine


def _blob(center: np.ndarray, radius: float = 8.0, n: int = 400, seed: int = 0) -> np.ndarray:
    r = np.random.default_rng(seed)
    p = r.normal(size=(n, 3))
    p /= np.linalg.norm(p, axis=1, keepdims=True)
    return center + p * (radius * r.uniform(0, 1, size=n) ** (1 / 3))[:, None]


def _bonds_exact(coords: np.ndarray) -> bool:
    bonds = np.linalg.norm(np.diff(coords, axis=1), axis=-1)
    return bool(np.max(np.abs(bonds - CA_CA_BOND_LENGTH)) < 1e-6)


N_ANCHOR = np.array([-21.0, 0.0, 0.0])
C_ANCHOR = np.array([21.0, 0.0, 0.0])


def _interior_request(n: int = 40, n_conformations: int = 8) -> IDRRequest:
    return IDRRequest(
        "A" * n,
        n,
        42.0,
        n_anchor_xyz=N_ANCHOR,
        c_anchor_xyz=C_ANCHOR,
        n_anchor_prev_xyz=np.array([-24.81, 0.0, 0.0]),
        c_anchor_next_xyz=np.array([24.81, 0.0, 0.0]),
        n_conformations=n_conformations,
    )


def _two_domains() -> np.ndarray:
    return np.vstack(
        [_blob(np.array([-36.0, 0.0, 0.0]), seed=1), _blob(np.array([36.0, 0.0, 0.0]), seed=2)]
    )


class TestProtocol:
    def test_satisfies_the_conformation_engine_protocol(self) -> None:
        engine = BatchWalkEngine()
        assert isinstance(engine, ConformationEngine)
        assert engine.name == "batch"
        assert engine.available()


class TestFreeAndTerminal:
    def test_free_region(self) -> None:
        engine = BatchWalkEngine()
        result = engine.generate(
            IDRRequest("A" * 30, 30, 40.0, n_conformations=6), None, np.random.default_rng(0)
        )
        assert isinstance(result, IDRResult)
        assert result.ca_coords.shape == (6, 30, 3)
        assert result.n_successful == 6
        assert np.all(np.isfinite(result.ca_coords[result.success]))
        assert _bonds_exact(result.ca_coords[result.success])

    def test_terminal_region_anchored_and_clearing(self) -> None:
        engine = BatchWalkEngine()
        obstacles = _blob(np.array([-15.0, 0.0, 0.0]))
        request = IDRRequest(
            "A" * 30,
            30,
            40.0,
            n_anchor_xyz=np.array([0.0, 0.0, 0.0]),
            n_anchor_prev_xyz=np.array([-3.81, 0.0, 0.0]),
            n_conformations=5,
        )
        result = engine.generate(request, obstacles, np.random.default_rng(0))
        assert result.ca_coords.shape == (5, 30, 3)
        assert result.n_successful == 5
        ok = result.ca_coords[result.success]
        first_gap = np.linalg.norm(ok[:, 0] - np.array([0.0, 0.0, 0.0]), axis=1)
        assert np.allclose(first_gap, CA_CA_BOND_LENGTH, atol=1e-3)
        assert _bonds_exact(ok)


class TestInterior:
    def test_builds_closes_and_anchors(self) -> None:
        engine = BatchWalkEngine()
        result = engine.generate(_interior_request(), _two_domains(), np.random.default_rng(0))
        assert result.ca_coords.shape == (8, 40, 3)
        assert result.n_successful == 8
        ok = result.ca_coords[result.success]
        assert _bonds_exact(ok)
        assert np.allclose(
            np.linalg.norm(ok[:, 0] - N_ANCHOR, axis=1), CA_CA_BOND_LENGTH, atol=1e-3
        )
        assert np.allclose(
            np.linalg.norm(ok[:, -1] - C_ANCHOR, axis=1), CA_CA_BOND_LENGTH, atol=1e-3
        )

    def test_deterministic_under_seed(self) -> None:
        engine = BatchWalkEngine()
        a = engine.generate(_interior_request(), _two_domains(), np.random.default_rng(0))
        b = engine.generate(_interior_request(), _two_domains(), np.random.default_rng(0))
        assert np.array_equal(np.nan_to_num(a.ca_coords), np.nan_to_num(b.ca_coords))
        assert np.array_equal(a.success, b.success)


class TestFallback:
    def test_short_interior_falls_back_to_the_walk(self) -> None:
        # Too short for the interior stubs (n < 2*stub+2 = 10): the batch can't build it, so the
        # whole region comes from the careful walk fallback -- and it still succeeds.
        engine = BatchWalkEngine()
        request = IDRRequest(
            "A" * 8,
            8,
            14.0,
            n_anchor_xyz=np.array([-7.0, 0.0, 0.0]),
            c_anchor_xyz=np.array([7.0, 0.0, 0.0]),
            n_anchor_prev_xyz=np.array([-10.8, 0.0, 0.0]),
            c_anchor_next_xyz=np.array([10.8, 0.0, 0.0]),
            n_conformations=3,
        )
        result = engine.generate(request, None, np.random.default_rng(0))
        assert result.ca_coords.shape == (3, 8, 3)
        assert "walk" in result.engine
        assert result.n_successful >= 1

    def test_a_forced_shortfall_is_topped_up_by_the_walk(self) -> None:
        # oversample=1 leaves almost no batch survivors for a low-yield interior, so most conformers
        # must come from the fallback -- the result still has n conformers, name records the mix.
        engine = BatchWalkEngine(oversample=1)
        result = engine.generate(
            _interior_request(n_conformations=12), _two_domains(), np.random.default_rng(0)
        )
        assert result.ca_coords.shape == (12, 40, 3)
        assert result.n_successful == 12
        assert "walk" in result.engine
        assert _bonds_exact(result.ca_coords[result.success])
