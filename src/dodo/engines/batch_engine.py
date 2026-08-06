"""A conformation engine that builds regions in vectorized batches, with a careful fallback.

The hybrid of Ryan's design: the PRIMARY path is the batched grow-then-filter generators in
:mod:`dodo.engines.batch_walk` (free, terminal, interior) -- cheap to over-generate, so a whole
region's conformers come from one vectorized pass and the self-avoiding, obstacle-clearing survivors
are kept. The FALLBACK is the careful self-avoiding walk, used only for the conformers the batch
could not supply: the hard tail (very short, very long, or very crowded regions) where the filter
survival is low, plus any region type the batch does not build.

Nothing here weakens the contract of :class:`~dodo.engines.base.IDRResult`: a returned conformer is
always finite and self-avoiding, a conformer that could not be built is reported through the success
mask (never as degenerate coordinates), and total failure raises.
"""

from __future__ import annotations

from dataclasses import replace
from typing import TYPE_CHECKING

import numpy as np

from ..constants import CA_CA_BOND_LENGTH
from ..exceptions import DodoError
from . import batch_walk as bw
from .base import IDRResult
from .walk import SelfAvoidingWalk

if TYPE_CHECKING:  # pragma: no cover
    from .base import ConformationEngine, IDRRequest

__all__ = ["BatchWalkEngine"]

#: Obstacle atoms within this distance of an anchor define the domain the region emerges from, whose
#: centroid gives the "grow away from here" direction for that anchor's stub.
_NEAR_ANCHOR_RADIUS: float = 25.0


class BatchWalkEngine:
    """Batched primary generation with a careful walk fallback (a ``ConformationEngine``)."""

    name = "batch"

    def __init__(
        self,
        fallback: ConformationEngine | None = None,
        *,
        oversample: int = 16,
    ) -> None:
        self.fallback: ConformationEngine = SelfAvoidingWalk() if fallback is None else fallback
        self.oversample = oversample

    def available(self) -> bool:
        """Return True -- pure geometry, with no optional dependency that could be missing."""
        return True

    def generate(
        self,
        request: IDRRequest,
        obstacles: np.ndarray | None,
        rng: np.random.Generator,
    ) -> IDRResult:
        """Generate one region's conformers: batched survivors first, careful walk for the rest."""
        n = request.n_conformations
        survivors = self._batch_survivors(request, obstacles, rng)

        if survivors.shape[0] >= n:
            selected = self._select(survivors, request, rng, n)
            return IDRResult.from_batch(selected, np.ones(n, dtype=bool), self.name, attempts=1)

        # Shortfall: keep what the batch cleanly produced, ask the careful walk for the remainder.
        clean = survivors
        got = clean.shape[0]
        try:
            fallback = self.fallback.generate(
                replace(request, n_conformations=n - got), obstacles, rng
            )
        except DodoError:
            # The walk could not build any either. A partial batch result is still a result; only a
            # totally empty one is a failure to raise.
            if got == 0:
                raise
            failed = np.full((n - got, request.n_residues, 3), np.nan)
            success = np.concatenate([np.ones(got, dtype=bool), np.zeros(n - got, dtype=bool)])
            return IDRResult.from_batch(
                np.concatenate([clean, failed], axis=0), success, self.name, attempts=1
            )

        coords = fallback.ca_coords if got == 0 else np.concatenate([clean, fallback.ca_coords])
        success = np.concatenate([np.ones(got, dtype=bool), fallback.success])
        engine = self.name if got == n else f"{self.name}+{fallback.engine}"
        return IDRResult.from_batch(
            coords,
            success,
            engine,
            attempts=1 + fallback.attempts,
            relaxed_to=fallback.relaxed_to,
            unconstrained_junctions=fallback.unconstrained_junctions,
        )

    def _batch_survivors(
        self,
        request: IDRRequest,
        obstacles: np.ndarray | None,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """Return a pool of clean (self-avoiding, clearing, closed) conformers, empty on failure.

        Every failure here -- a region too short for the interior stubs, a batch where nothing
        survived, an unsupported shape -- is swallowed to an empty pool, because the caller's answer
        to an empty pool is simply to fall back to the walk.
        """
        pool = max(request.n_conformations * self.oversample, self.oversample)
        try:
            coords, feasible = self._grow_pool(request, obstacles, rng, pool)
        except DodoError:
            return np.empty((0, request.n_residues, 3))
        clean = feasible & bw.self_avoiding_mask(coords)
        if obstacles is not None and np.asarray(obstacles).shape[0]:
            clean = clean & bw.clears_obstacles_mask(coords, np.asarray(obstacles))
        selected: np.ndarray = coords[clean]
        return selected

    def _grow_pool(
        self,
        request: IDRRequest,
        obstacles: np.ndarray | None,
        rng: np.random.Generator,
        pool: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Grow ``pool`` candidates, returning coordinates and a closure-feasible mask."""
        n_res = request.n_residues
        n_anchor = request.n_anchor_xyz
        c_anchor = request.c_anchor_xyz
        if n_anchor is not None and c_anchor is not None:
            n_away = self._away(n_anchor, request.n_anchor_prev_xyz, obstacles)
            c_away = self._away(c_anchor, request.c_anchor_next_xyz, obstacles)
            from ..construct.dimensions import predict_radius_of_gyration

            target_rg, _ = predict_radius_of_gyration(request.sequence)
            result = bw.generate_interior_batch(
                pool,
                n_res,
                n_anchor,
                c_anchor,
                rng,
                n_away=n_away,
                c_away=c_away,
                target_rg=target_rg,
            )
            return result.coords, result.closure_feasible
        if n_anchor is not None or c_anchor is not None:
            if n_anchor is not None:
                anchor, outer = n_anchor, request.n_anchor_prev_xyz
            else:
                assert c_anchor is not None  # the outer branch guarantees one anchor is set
                anchor, outer = c_anchor, request.c_anchor_next_xyz
            away = self._away(anchor, outer, obstacles)
            terminal = bw.generate_terminal_batch(
                pool, n_res, anchor, away, request.target_end_to_end, rng, obstacles=obstacles
            )
            return terminal.coords, np.ones(terminal.coords.shape[0], dtype=bool)
        free = bw.generate_free_batch(pool, n_res, request.target_end_to_end, rng)
        return free.coords, np.ones(free.coords.shape[0], dtype=bool)

    def _select(
        self,
        survivors: np.ndarray,
        request: IDRRequest,
        rng: np.random.Generator,
        n: int,
    ) -> np.ndarray:
        """Pick ``n`` conformers from the survivors, steered to the target where it is free."""
        if request.is_interior:
            # An interior region's end-to-end is fixed by the anchors; there is nothing to steer, so
            # take a random distinct subset and let the paths differ.
            chosen: np.ndarray = survivors[rng.permutation(survivors.shape[0])[:n]]
            return chosen
        low = CA_CA_BOND_LENGTH
        high = (request.n_residues - 1) * CA_CA_BOND_LENGTH * 0.98
        steered: np.ndarray = bw.steer_to_target(
            survivors, request.target_end_to_end, n, rng, low=low, high=high
        )
        return steered

    def _away(
        self,
        anchor: np.ndarray,
        outer: np.ndarray | None,
        obstacles: np.ndarray | None,
    ) -> np.ndarray:
        """Return the unit direction to grow away from the anchor's domain (nearby obstacles)."""
        anchor = np.asarray(anchor, dtype=np.float64)
        if obstacles is not None:
            points = np.asarray(obstacles, dtype=np.float64)
            if points.shape[0]:
                near = points[np.linalg.norm(points - anchor, axis=1) < _NEAR_ANCHOR_RADIUS]
                if near.shape[0] >= 3:
                    direction = anchor - near.mean(axis=0)
                    norm = float(np.linalg.norm(direction))
                    if norm > 1e-6:
                        unit: np.ndarray = direction / norm
                        return unit
        if outer is not None:
            # No usable obstacle cloud: continue the direction the chain entered the anchor from.
            direction = anchor - np.asarray(outer, dtype=np.float64)
            norm = float(np.linalg.norm(direction))
            if norm > 1e-6:
                unit = direction / norm
                return unit
        return np.array([1.0, 0.0, 0.0])
