"""Vectorized batched CA growth: grow many alpha-carbon traces in lockstep, filter after.

This is the foundation of DODO's vectorized IDR engine (Phase 3). It exists to make an arbitrary
number of independent traces -- many candidate IDRs, or many conformers of one IDR -- in a single
pass whose wall time scales with the residue count, not with how many traces are asked for.

The idea, and why it is vectorizable where the self-avoiding walk is not
-----------------------------------------------------------------------
The self-avoiding walk in :mod:`dodo.engines.walk` rejects each step that clashes and retries, so
different chains stall at different points and cannot be advanced in lockstep. Here growth is
**unconditional**: every chain takes a step every round, sampled so the CA-CA bond is exactly
:data:`~dodo.constants.CA_CA_BOND_LENGTH` and the CA-CA-CA pseudo-angle lands inside the measured
window by construction. Nothing is rejected during growth, so a whole batch advances together as
one array operation per residue.

Self-avoidance is then recovered by **filtering**: a batch is over-generated cheaply and the traces
that happen to collide with themselves are dropped afterwards with one vectorized distance check
(:func:`self_avoiding_mask`). Because growth is cheap, over-generating to cover the discards is
cheaper than steering every step around every other atom.

What this module does NOT yet do -- deliberately, and measured before adding
--------------------------------------------------------------------------
It grows FREE chains at their natural dimension. It does not yet steer toward an ALBATROSS-predicted
end-to-end distance, close onto two fixed anchors, or avoid external obstacles. Those are the later
Phase 3 slices; the natural end-to-end distribution and self-avoidance survival that this module
produces are what the design of those slices is measured against, rather than assumed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import truncnorm

from ..constants import (
    BACKBONE_ANGLE_MAX,
    BACKBONE_ANGLE_MEAN,
    BACKBONE_ANGLE_MIN,
    BACKBONE_ANGLE_SD,
    CA_CA_BOND_LENGTH,
    CA_CLASH_DISTANCE,
    CLASH_EXCLUDE_WITHIN_RESIDUES,
)
from ..exceptions import EngineError
from .walk import sample_end_to_end_targets

__all__ = [
    "FreeBatchResult",
    "InteriorBatchResult",
    "clears_obstacles_mask",
    "close_chain_ends",
    "end_to_end_distances",
    "generate_free_batch",
    "generate_interior_batch",
    "generate_terminal_batch",
    "grow_free_batch",
    "radii_of_gyration",
    "self_avoiding_mask",
    "steer_to_target",
]

#: A bond this far from :data:`~dodo.constants.CA_CA_BOND_LENGTH` marks a chain as not cleanly
#: closed -- the funnel's analytic landing could not be satisfied and left a stretched junction.
_CLOSURE_BOND_TOLERANCE: float = 0.1


def _unit_rows(v: np.ndarray) -> np.ndarray:
    """Normalize along the last axis, flooring the divisor so no row can divide by zero."""
    norm = np.linalg.norm(v, axis=-1, keepdims=True)
    unit: np.ndarray = v / np.maximum(norm, 1e-12)
    return unit


def _random_unit_vectors(rng: np.random.Generator, n: int) -> np.ndarray:
    """``(n, 3)`` isotropically-random unit vectors."""
    gaussian = rng.normal(size=(n, 3))
    return _unit_rows(gaussian)


def _perpendicular_frame(u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Two orthonormal vectors perpendicular to each row of the unit array ``u``.

    ``u`` is ``(n, 3)`` and assumed already normalized. The seed axis is z unless a row is nearly
    parallel to z, where y is used instead, so the cross product never degenerates.
    """
    seed = np.where(np.abs(u[..., 2:3]) > 0.9, np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0]))
    e1 = _unit_rows(seed - np.sum(seed * u, axis=-1, keepdims=True) * u)
    e2 = np.cross(u, e1)
    return e1, e2


def grow_free_batch(
    n_chains: int,
    n_residues: int,
    rng: np.random.Generator,
    *,
    bond: float = CA_CA_BOND_LENGTH,
    angle_mean: float = BACKBONE_ANGLE_MEAN,
    angle_sd: float = BACKBONE_ANGLE_SD,
    angle_min: float = BACKBONE_ANGLE_MIN,
    angle_max: float = BACKBONE_ANGLE_MAX,
    bias_directions: np.ndarray | None = None,
    bias_residues: int = 0,
    bias_kappa: float = 2.0,
) -> np.ndarray:
    """Grow ``n_chains`` freely-rotating CA traces of ``n_residues`` in one vectorized pass.

    Every CA-CA bond is exactly ``bond`` and every CA-CA-CA pseudo-angle is drawn from a normal
    distribution ``(angle_mean, angle_sd)`` truncated to ``[angle_min, angle_max]`` -- both by
    construction, so neither needs checking afterwards. Self-avoidance is NOT enforced; pass the
    result through :func:`self_avoiding_mask` and keep the survivors.

    Parameters
    ----------
    n_chains
        Number of independent traces to grow together. This is the batch axis; wall time is
        near-flat in it until memory bandwidth bites.
    n_residues
        Alpha carbons per trace. At least two.
    rng
        Seeded generator. Every random choice is drawn from it, so a fixed seed reproduces the
        batch exactly.
    bond
        CA-CA virtual bond length.
    angle_mean, angle_sd, angle_min, angle_max
        The pseudo-angle distribution, in degrees.
    bias_directions, bias_residues, bias_kappa
        Optional soft skew for the first ``bias_residues`` bonds toward ``bias_directions``
        (``(n_chains, 3)``, one direction per chain -- e.g. away from a folded domain's centre of
        mass), with concentration ``bias_kappa``. Only the azimuth is skewed, so the bond and
        pseudo-angle are untouched, and it is a von Mises DISTRIBUTION, so every conformer's near-
        anchor stub leans toward open space differently rather than being the same fixed shape.
        ``None`` (default) grows isotropically and is bit-identical to the unbiased path.

    Returns
    -------
    np.ndarray
        ``(n_chains, n_residues, 3)`` coordinates. The first residue of every chain is the origin.
    """
    if n_chains < 1:
        raise EngineError(f"n_chains must be at least 1, got {n_chains}.")
    if n_residues < 2:
        raise EngineError(f"n_residues must be at least 2 to have one bond, got {n_residues}.")

    biased = bias_directions is not None and bias_residues >= 1
    bias = None if bias_directions is None else _unit_rows(np.asarray(bias_directions, np.float64))

    coords = np.empty((n_chains, n_residues, 3), dtype=np.float64)
    coords[:, 0, :] = 0.0
    if biased and bias is not None:
        # First bond leans toward the bias direction, with per-chain scatter so conformers differ.
        first = _unit_rows(bias_kappa * bias + rng.normal(size=(n_chains, 3)))
    else:
        first = _random_unit_vectors(rng, n_chains)
    coords[:, 1, :] = coords[:, 0, :] + bond * first
    if n_residues == 2:
        return coords

    # The pseudo-angle theta at CA(i) is CA(i-1)-CA(i)-CA(i+1). With u the continuation direction
    # unit(CA(i) - CA(i-1)), the outgoing direction v satisfies dot(u, v) = -cos(theta), so v sits
    # on a cone of half-angle beta = 180 - theta about u, at an azimuth that is uniform by default
    # and von-Mises-skewed toward the bias direction for the first ``bias_residues`` steps.
    lo = (angle_min - angle_mean) / angle_sd
    hi = (angle_max - angle_mean) / angle_sd
    for i in range(1, n_residues - 1):
        u = _unit_rows(coords[:, i, :] - coords[:, i - 1, :])
        theta = truncnorm.rvs(
            lo, hi, loc=angle_mean, scale=angle_sd, size=n_chains, random_state=rng
        )
        beta = np.radians(180.0 - theta)
        e1, e2 = _perpendicular_frame(u)
        if biased and bias is not None and i < bias_residues:
            phi_star = np.arctan2(np.sum(bias * e2, axis=-1), np.sum(bias * e1, axis=-1))
            azimuth = rng.vonmises(phi_star, bias_kappa)
        else:
            azimuth = rng.uniform(0.0, 2.0 * np.pi, size=n_chains)
        radial = np.cos(azimuth)[:, None] * e1 + np.sin(azimuth)[:, None] * e2
        direction = np.cos(beta)[:, None] * u + np.sin(beta)[:, None] * radial
        coords[:, i + 1, :] = coords[:, i, :] + bond * direction

    return coords


def self_avoiding_mask(
    coords: np.ndarray,
    *,
    clash_distance: float = CA_CLASH_DISTANCE,
    exclude_within: int = CLASH_EXCLUDE_WITHIN_RESIDUES,
    chunk_pairs: int = 4_000_000,
) -> np.ndarray:
    """Boolean mask of the traces in a batch that do not clash with themselves.

    A trace is kept when every pair of alpha carbons more than ``exclude_within`` residues apart is
    at least ``clash_distance`` apart in space. Sequence-local pairs are skipped because their
    separation is set by the fixed bond and angle, not by packing, exactly as
    :meth:`dodo.structure.Structure.clash_mask` excludes them.

    Parameters
    ----------
    coords
        ``(n_chains, n_residues, 3)`` batch.
    clash_distance
        Minimum allowed CA-CA separation for a non-local pair.
    exclude_within
        Sequence separation at or below which a pair is bonded-close and not scored.
    chunk_pairs
        Upper bound on ``chain_chunk * n_residues**2`` held at once, so the pairwise test stays
        within a bounded amount of memory however large the batch is.

    Returns
    -------
    np.ndarray
        ``(n_chains,)`` boolean, True where the trace is self-avoiding.
    """
    coords = np.asarray(coords, dtype=np.float64)
    if coords.ndim != 3 or coords.shape[2] != 3:
        raise EngineError(f"coords must be (n_chains, n_residues, 3), got {coords.shape}.")
    n_chains, n_residues, _ = coords.shape
    keep = np.ones(n_chains, dtype=bool)
    if n_residues <= exclude_within + 1:
        return keep

    separation = np.abs(np.arange(n_residues)[:, None] - np.arange(n_residues)[None, :])
    scored = separation > exclude_within
    threshold_squared = clash_distance * clash_distance

    step = max(1, chunk_pairs // (n_residues * n_residues))
    for start in range(0, n_chains, step):
        block = coords[start : start + step]
        diff = block[:, :, None, :] - block[:, None, :, :]
        squared = np.einsum("cijd,cijd->cij", diff, diff)
        squared[:, ~scored] = np.inf
        keep[start : start + step] = squared.reshape(block.shape[0], -1).min(axis=1) >= (
            threshold_squared
        )
    return keep


def end_to_end_distances(coords: np.ndarray) -> np.ndarray:
    """``(n_chains,)`` CA-to-CA distance from the first residue to the last, per trace."""
    distances: np.ndarray = np.linalg.norm(coords[:, -1, :] - coords[:, 0, :], axis=-1)
    return distances


def radii_of_gyration(coords: np.ndarray) -> np.ndarray:
    """``(n_chains,)`` radius of gyration of each trace's alpha-carbon cloud."""
    centred = coords - coords.mean(axis=1, keepdims=True)
    radii: np.ndarray = np.sqrt((centred**2).sum(axis=(1, 2)) / coords.shape[1])
    return radii


def clears_obstacles_mask(
    coords: np.ndarray,
    obstacles: np.ndarray,
    *,
    clash_distance: float = CA_CLASH_DISTANCE,
    chunk_pairs: int = 4_000_000,
) -> np.ndarray:
    """Boolean mask of the traces whose every alpha carbon clears an external obstacle set.

    The counterpart of :func:`self_avoiding_mask` for atoms the chain did not generate -- folded
    domains especially. A trace is kept when no alpha carbon comes within ``clash_distance`` of any
    obstacle point.
    """
    coords = np.asarray(coords, dtype=np.float64)
    obstacles = np.asarray(obstacles, dtype=np.float64)
    n_chains = coords.shape[0]
    keep = np.ones(n_chains, dtype=bool)
    if obstacles.shape[0] == 0 or coords.shape[1] == 0:
        return keep
    n_residues, n_obstacles = coords.shape[1], obstacles.shape[0]
    threshold_squared = clash_distance * clash_distance
    step = max(1, chunk_pairs // (n_residues * n_obstacles))
    for start in range(0, n_chains, step):
        block = coords[start : start + step]
        diff = block[:, :, None, :] - obstacles[None, None, :, :]
        squared = np.einsum("cikd,cikd->cik", diff, diff)
        keep[start : start + step] = (
            squared.reshape(block.shape[0], -1).min(axis=1) >= threshold_squared
        )
    return keep


def _nearest_indices(sorted_values: np.ndarray, queries: np.ndarray) -> np.ndarray:
    """Index into ascending ``sorted_values`` of the entry nearest each query."""
    top = sorted_values.shape[0] - 1
    pos = np.clip(np.searchsorted(sorted_values, queries), 0, top)
    below = np.clip(pos - 1, 0, top)
    pick_below = np.abs(sorted_values[below] - queries) <= np.abs(sorted_values[pos] - queries)
    chosen: np.ndarray = np.where(pick_below, below, pos)
    return chosen


def _select_from_pool(
    pool_coords: np.ndarray,
    target_end_to_end: float,
    n_conformers: int,
    rng: np.random.Generator,
    low: float,
    high: float,
) -> tuple[np.ndarray, float, float]:
    """Match ``n_conformers`` physical draws to the nearest-Re pool traces. Shared core."""
    re = end_to_end_distances(pool_coords)
    order = np.argsort(re, kind="stable")
    re_sorted = re[order]
    drawn = sample_end_to_end_targets(target_end_to_end, n_conformers, rng, low=low, high=high)
    chosen = order[_nearest_indices(re_sorted, drawn)]
    reachable = float(np.mean((drawn >= re_sorted[0]) & (drawn <= re_sorted[-1])))
    reuse = 1.0 - np.unique(chosen).shape[0] / n_conformers
    return pool_coords[chosen], reachable, float(reuse)


def steer_to_target(
    pool_coords: np.ndarray,
    target_end_to_end: float,
    n_conformers: int,
    rng: np.random.Generator,
    *,
    low: float,
    high: float,
) -> np.ndarray:
    """Select ``n_conformers`` traces from a pool so their end-to-end distances match a target.

    The pool -- typically the self-avoiding survivors of :func:`grow_free_batch` -- has its own
    natural end-to-end distribution. This draws one physical target per conformer around
    ``target_end_to_end`` with :func:`dodo.engines.walk.sample_end_to_end_targets`, the SAME
    distribution the self-avoiding walk steers to, then matches each draw to the pool trace whose
    end-to-end distance is nearest. The result has mean ``target_end_to_end`` and the physical
    spread, rather than a spike at the mean -- selecting a narrow band would collapse the ensemble,
    which is the failure the walk's sampler exists to avoid.

    Matching is with replacement: where the pool is thin near a drawn target the same trace can be
    chosen twice. Grow a larger pool to reduce it; :func:`generate_free_batch` reports how often it
    happened.
    """
    if pool_coords.shape[0] == 0:
        raise EngineError("Cannot steer an empty pool; every grown trace clashed with itself.")
    selected, _, _ = _select_from_pool(pool_coords, target_end_to_end, n_conformers, rng, low, high)
    return selected


@dataclass(frozen=True, slots=True)
class FreeBatchResult:
    """A steered batch of free conformers, with the diagnostics a caller needs to trust it.

    Attributes
    ----------
    coords
        ``(n_conformers, n_residues, 3)`` selected conformers, steered to the target mean.
    achieved_end_to_end
        ``(n_conformers,)`` end-to-end distance of each returned conformer.
    target_end_to_end
        The requested mean.
    n_grown, n_survived
        Traces grown, and how many survived the self-avoidance filter before steering.
    reachable_fraction
        Fraction of the drawn per-conformer targets that fell inside the survivor pool's
        end-to-end range. Below 1 means the target sits partly outside what free growth reaches at
        this length -- a bias or a fallback is needed, not selection.
    reuse_fraction
        Fraction of returned conformers that reused an already-selected survivor. High reuse means
        the pool was too thin: grow more, or fall back.
    """

    coords: np.ndarray
    achieved_end_to_end: np.ndarray
    target_end_to_end: float
    n_grown: int
    n_survived: int
    reachable_fraction: float
    reuse_fraction: float


def generate_free_batch(
    n_conformers: int,
    n_residues: int,
    target_end_to_end: float,
    rng: np.random.Generator,
    *,
    oversample: int = 16,
    low: float | None = None,
    high: float | None = None,
) -> FreeBatchResult:
    """Grow, self-avoidance-filter and steer a batch of free conformers to a target end-to-end.

    The primary path for a free-ended (terminal) region: over-generate cheaply, drop the traces
    that clash with themselves, and select from the survivors so the batch matches the target mean
    with a physical spread. ``oversample`` sets how many traces are grown per requested conformer,
    to leave room after the self-avoidance discards; the returned :class:`FreeBatchResult` reports
    whether that was enough (``n_survived``, ``reuse_fraction``) and whether the target was even in
    range (``reachable_fraction``), which is what a caller uses to decide to fall back.
    """
    if n_conformers < 1:
        raise EngineError(f"n_conformers must be at least 1, got {n_conformers}.")
    n_grown = n_conformers * oversample
    pool = grow_free_batch(n_grown, n_residues, rng)
    survivors = pool[self_avoiding_mask(pool)]
    if survivors.shape[0] == 0:
        raise EngineError(
            f"Every one of {n_grown} grown traces clashed with itself at {n_residues} residues; "
            f"free growth cannot serve this length -- use a fallback instead."
        )
    bond = CA_CA_BOND_LENGTH
    lo = bond if low is None else low
    hi = (n_residues - 1) * bond * 0.98 if high is None else high
    selected, reachable, reuse = _select_from_pool(
        survivors, target_end_to_end, n_conformers, rng, lo, hi
    )
    return FreeBatchResult(
        coords=selected,
        achieved_end_to_end=end_to_end_distances(selected),
        target_end_to_end=target_end_to_end,
        n_grown=n_grown,
        n_survived=int(survivors.shape[0]),
        reachable_fraction=reachable,
        reuse_fraction=reuse,
    )


def generate_terminal_batch(
    n_conformers: int,
    n_residues: int,
    anchor: np.ndarray,
    away: np.ndarray,
    target_end_to_end: float,
    rng: np.random.Generator,
    *,
    obstacles: np.ndarray | None = None,
    stub: int = 4,
    oversample: int = 16,
    low: float | None = None,
    high: float | None = None,
) -> FreeBatchResult:
    """Grow, self-avoidance-filter and steer a free-ended (terminal) IDR anchored at one end.

    The counterpart of :func:`generate_free_batch` for a region with one fixed end. Residue 0 sits
    one bond from ``anchor``; the chain grows away from that anchor's domain -- its first ``stub``
    steps skewed toward ``away`` (the unit direction from the domain's centre of mass to the anchor)
    so it leaves the crowded shell -- and the self-avoiding survivors that also clear ``obstacles``
    are steered to the target end-to-end. ``FreeBatchResult.n_survived`` / ``reuse_fraction`` report
    whether the over-generation covered the discards.
    """
    if n_conformers < 1:
        raise EngineError(f"n_conformers must be at least 1, got {n_conformers}.")
    bond = CA_CA_BOND_LENGTH
    n_grown = n_conformers * oversample
    bias = np.tile(_unit_rows(np.asarray(away, dtype=np.float64)), (n_grown, 1))
    grown = grow_free_batch(n_grown, n_residues + 1, rng, bias_directions=bias, bias_residues=stub)
    grown = grown - grown[:, 0:1] + np.asarray(anchor, dtype=np.float64)
    coords = grown[:, 1:]  # IDR residues 0..n_residues-1; residue 0 one bond from the anchor
    keep = self_avoiding_mask(coords)
    obstacle_points = None if obstacles is None else np.asarray(obstacles, dtype=np.float64)
    if obstacle_points is not None and obstacle_points.shape[0]:
        keep = keep & clears_obstacles_mask(coords, obstacle_points)
    survivors = coords[keep]
    if survivors.shape[0] == 0:
        raise EngineError(
            f"Every one of {n_grown} grown terminal traces clashed at {n_residues} residues; "
            f"grow more or fall back."
        )
    lo = bond if low is None else low
    hi = (n_residues - 1) * bond * 0.98 if high is None else high
    selected, reachable, reuse = _select_from_pool(
        survivors, target_end_to_end, n_conformers, rng, lo, hi
    )
    return FreeBatchResult(
        coords=selected,
        achieved_end_to_end=end_to_end_distances(selected),
        target_end_to_end=target_end_to_end,
        n_grown=n_grown,
        n_survived=int(survivors.shape[0]),
        reachable_fraction=reachable,
        reuse_fraction=reuse,
    )


def close_chain_ends(
    coords: np.ndarray,
    start: np.ndarray,
    end: np.ndarray,
    *,
    bond: float = CA_CA_BOND_LENGTH,
    iterations: int = 16,
) -> np.ndarray:
    """Pin both ends of each chain to ``start``/``end`` and restore every bond to ``bond``.

    A batched FABRIK-style endpoint closure -- the reach-from-both-ends analogue of a
    SHAKE-with-endpoints projection. Each iteration sweeps backward from the fixed ``end`` and then
    forward from the fixed ``start``, driving every consecutive distance to ``bond`` while keeping
    the two ends pinned. It converges when the span ``|start - end|`` is at most the contour
    length ``(m - 1) * bond``; a shorter chain cannot bridge the gap, so it falls short of ``end``
    while keeping every bond exact -- a residual ``|closed[:, -1] - end|`` the caller treats as
    infeasible.

    ``coords`` is ``(n_chains, m, 3)`` and is not modified; a closed copy is returned.
    """
    out = np.array(coords, dtype=np.float64, copy=True)
    m = out.shape[1]
    if m < 2:
        return out
    start = np.asarray(start, dtype=np.float64)
    end = np.asarray(end, dtype=np.float64)
    for _ in range(iterations):
        out[:, -1] = end
        for i in range(m - 2, -1, -1):
            step = out[:, i] - out[:, i + 1]
            step /= np.maximum(np.linalg.norm(step, axis=-1, keepdims=True), 1e-12)
            out[:, i] = out[:, i + 1] + bond * step
        out[:, 0] = start
        for i in range(1, m):
            step = out[:, i] - out[:, i - 1]
            step /= np.maximum(np.linalg.norm(step, axis=-1, keepdims=True), 1e-12)
            out[:, i] = out[:, i - 1] + bond * step
    return out


def _sphere_sphere_point(
    centre_a: np.ndarray, centre_b: np.ndarray, radius: float, toward: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Point at distance ``radius`` from both centres, on that circle, nearest ``toward``.

    Batched. Two equal-radius spheres meet in a circle when their centres are within ``2 * radius``;
    ``valid`` marks the rows where they do. Used to place the last middle residue so both of its
    bonds -- to the previous residue and to the far anchor -- come out exact.
    """
    diff = centre_b - centre_a
    sep = np.linalg.norm(diff, axis=-1)
    axis = diff / np.maximum(sep, 1e-12)[:, None]
    half = sep / 2.0
    height_squared = radius * radius - half * half
    valid = (sep > 1e-9) & (sep <= 2.0 * radius)
    centre = centre_a + half[:, None] * axis
    spoke = toward - centre
    spoke = _unit_rows(spoke - np.sum(spoke * axis, axis=-1, keepdims=True) * axis)
    point: np.ndarray = centre + np.sqrt(np.maximum(height_squared, 0.0))[:, None] * spoke
    return point, valid


def _bridge_funnel(
    p: np.ndarray,
    q: np.ndarray,
    u_start: np.ndarray,
    n_mid: int,
    rng: np.random.Generator,
    *,
    kappa_min: float = 0.5,
    kappa_max: float = 10.0,
    power: float = 3.0,
    reach_factor: float = 0.90,
    excursion: float = 0.0,
    expansion: float = 0.0,
    landing: int = 0,
) -> np.ndarray:
    """Grow the ``n_mid`` middle residues from ``p`` to ``q`` inside a shrinking reachability ball.

    Each step is an in-window cone step (exact bond, pseudo-angle in the measured window), but its
    azimuth is von-Mises-skewed toward a TARGET DIRECTION that superposes three drifts, then the
    funnel front is closed by a single analytic residue. Because closure is enforced during growth
    rather than projected afterwards, the angle window survives, unlike :func:`close_chain_ends`.

    The drifts, each shifting only the mean of the azimuth distribution so every conformer stays
    unique:

    * reachability -- toward the target, concentration rising (``kappa``) as the reachable ball
      shrinks to the remaining contour, which is what actually closes the chain;
    * excursion (``excursion`` > 0) -- a per-conformer RANDOM direction perpendicular to
      ``p``->``q``, faded out as the reach tightens, so a slack chain bulges sideways into a spread
      arc instead of collapsing inward and self-clashing (the random direction also adds diversity);
    * expansion (``expansion`` > 0) -- away from the chain's own running centroid, a soft
      self-repulsion that keeps it from piling up.

    ``landing`` grows that many residues BACKWARD from ``q`` toward ``p`` first and aims the funnel
    at the DEEPEST of them rather than at ``q``, the intent being a near mid-gap soft target. It is
    OFF by default (``landing == 0``, aim straight at ``q``) because it was MEASURED not to help:
    the reverse pad occupies the same mid-gap the funnel front does and, grown blind to it, the two
    collide -- self-avoidance fell from 57% to 46% at n=40 while closure barely moved, since the
    closure bottleneck is the funnel's last-step aiming, not the target distance. Kept as an
    off-by-default knob. ``landing == 0`` with ``excursion == expansion == 0`` is bit-identical to
    the pure reachability funnel.
    """
    bond = CA_CA_BOND_LENGTH
    n_chains = p.shape[0]
    lo = (BACKBONE_ANGLE_MIN - BACKBONE_ANGLE_MEAN) / BACKBONE_ANGLE_SD
    hi = (BACKBONE_ANGLE_MAX - BACKBONE_ANGLE_MEAN) / BACKBONE_ANGLE_SD

    e_perp = np.zeros((n_chains, 3))
    if excursion > 0.0:
        axis_pq = _unit_rows(q - p)
        gaussian = rng.normal(size=(n_chains, 3))
        e_perp = _unit_rows(gaussian - np.sum(gaussian * axis_pq, axis=-1, keepdims=True) * axis_pq)

    # Landing pad grown backward from q toward p: its deepest point is a near, mid-gap soft target.
    n_pad = int(np.clip(landing, 0, max(0, n_mid - 2)))
    landing_res: np.ndarray | None = None
    target = q
    if n_pad >= 1:
        pad = grow_free_batch(
            n_chains, n_pad + 1, rng, bias_directions=_unit_rows(p - q), bias_residues=n_pad
        )
        pad = pad - pad[:, 0:1] + q[:, None, :]  # pad[:,0]=q, pad[:,1..n_pad] head toward p
        target = pad[:, n_pad]  # deepest landing point, ~n_pad bonds into the gap
        landing_res = pad[:, 1 : n_pad + 1][:, ::-1, :]  # sequence order: deepest .. nearest q

    n_front = n_mid - n_pad - 1  # funnel residues before the single closing residue
    r = np.empty((n_chains, n_front + 2, 3), dtype=np.float64)
    r[:, 0] = p
    for i in range(n_front):
        u = u_start if i == 0 else _unit_rows(r[:, i] - r[:, i - 1])
        theta = truncnorm.rvs(
            lo,
            hi,
            loc=BACKBONE_ANGLE_MEAN,
            scale=BACKBONE_ANGLE_SD,
            size=n_chains,
            random_state=rng,
        )
        beta = np.radians(180.0 - theta)
        e1, e2 = _perpendicular_frame(u)
        # Bonds from r_i to the target run through the two closing bonds, hence the + 2.
        reach = (n_front - i + 2) * bond * reach_factor
        tightness = np.clip(np.linalg.norm(r[:, i] - target, axis=-1) / reach, 0.0, 1.0)
        drift = _unit_rows(target - r[:, i])  # reachability: toward the (soft) target
        if excursion > 0.0:
            drift = drift + (excursion * (1.0 - tightness))[:, None] * e_perp
        if expansion > 0.0:
            centroid = r[:, : i + 1].mean(axis=1)
            drift = drift + expansion * _unit_rows(r[:, i] - centroid)
        phi_star = np.arctan2(np.sum(drift * e2, axis=-1), np.sum(drift * e1, axis=-1))
        kappa = kappa_min + (kappa_max - kappa_min) * tightness**power
        azimuth = rng.vonmises(phi_star, kappa)
        radial = np.cos(azimuth)[:, None] * e1 + np.sin(azimuth)[:, None] * e2
        r[:, i + 1] = r[:, i] + bond * (np.cos(beta)[:, None] * u + np.sin(beta)[:, None] * radial)

    # One analytic residue on the intersection circle closes the funnel front onto the target.
    prev = r[:, n_front]
    u_last = _unit_rows(prev - r[:, n_front - 1]) if n_front >= 2 else u_start
    point, _ = _sphere_sphere_point(prev, target, bond, prev + u_last * bond)
    r[:, n_front + 1] = point

    front = r[:, 1 : n_front + 2]  # funnel residues plus the closing residue
    if landing_res is None:
        return front
    joined: np.ndarray = np.concatenate([front, landing_res], axis=1)
    return joined


def _excursion_for_target_rg(
    p: np.ndarray,
    q: np.ndarray,
    u_start: np.ndarray,
    stub_n: np.ndarray,
    stub_c: np.ndarray,
    n_mid: int,
    rng: np.random.Generator,
    target_rg: float,
    *,
    exc_min: float = 2.0,
    exc_max: float = 3.5,
    probe: int = 512,
) -> tuple[float, bool]:
    """Pick the excursion whose spread makes the interior chain's radius of gyration match a target.

    Returns the excursion and whether it had to be CLAMPED because the target fell outside the
    reachable ``[Rg(exc_min), Rg(exc_max)]`` band -- a clamp means the achieved Rg cannot match it.

    The excursion sets how far the slack middle bulges, hence the chain's Rg, monotonically. This
    grows a small probe batch across the band, measures the assembled chain's Rg at each, and
    interpolates for the value that hits ``target_rg``, clamped to ``[exc_min, exc_max]``. That band
    is deliberately NOT ``[0, huge]`` (measured n=40): below ``exc_min`` the chain has too little
    perpendicular spread and self-clashes (self-avoidance < ~45%); above ``exc_max`` the funnel
    wanders and the analytic closure fails (closure < ~65%). Clean yield peaks inside it, and its Rg
    ceiling sits below the FREE-chain predicted Rg -- which a doubly-pinned chain cannot reach at
    healthy closure anyway -- so clamping here targets the pinned-chain Rg instead of chasing an
    unreachable one. A probe rather than a table because the Rg(excursion) curve depends on length,
    gap and the stub geometry, which are only known per call.
    """
    take = min(probe, p.shape[0])
    candidates = np.array([exc_min, 0.5 * (exc_min + exc_max), exc_max])
    achieved = np.empty(candidates.shape[0])
    for k, value in enumerate(candidates):
        mid = _bridge_funnel(p[:take], q[:take], u_start[:take], n_mid, rng, excursion=float(value))
        chain = np.concatenate([stub_n[:take], mid, stub_c[:take]], axis=1)
        achieved[k] = float(radii_of_gyration(chain).mean())
    achieved = np.maximum.accumulate(achieved)  # guard interp against probe noise
    if target_rg <= achieved[0]:
        return exc_min, True
    if target_rg >= achieved[-1]:
        return exc_max, True
    return float(np.interp(target_rg, achieved, candidates)), False


@dataclass(frozen=True, slots=True)
class InteriorBatchResult:
    """An interior IDR batch with the diagnostics a caller needs to filter it and to trust the Rg.

    Attributes
    ----------
    coords
        ``(n_conformers, n_residues, 3)`` traces bridging the two anchors.
    closure_feasible
        ``(n_conformers,)`` bool, True where every bond is within
        :data:`_CLOSURE_BOND_TOLERANCE` of 3.81 A. The funnel closure is over-generate-then-filter:
        a chain whose analytic landing could not be satisfied carries a long junction bond, and this
        flag is how a caller tells a clean chain from a broken one -- a broken chain can otherwise
        carry a perfectly on-target radius of gyration, which is the trap this exists to close.
    radius_of_gyration
        ``(n_conformers,)`` Rg of each trace.
    target_rg
        The requested Rg, or ``None`` if the excursion was set directly.
    rg_clamped
        True if ``target_rg`` fell outside the reachable band and the excursion was clamped, so the
        achieved Rg cannot reach it (a chain pinned at both ends has a spread ceiling and floor).
    excursion
        The perpendicular-excursion magnitude used (probe-chosen when ``target_rg`` was given).
        Not meaningful for ``closure='fabrik'``.
    """

    coords: np.ndarray
    closure_feasible: np.ndarray
    radius_of_gyration: np.ndarray
    target_rg: float | None
    rg_clamped: bool
    excursion: float

    @property
    def feasible_fraction(self) -> float:
        """Fraction of conformers that closed with every bond intact."""
        return float(np.mean(self.closure_feasible))


def generate_interior_batch(
    n_conformers: int,
    n_residues: int,
    n_anchor: np.ndarray,
    c_anchor: np.ndarray,
    rng: np.random.Generator,
    *,
    n_away: np.ndarray,
    c_away: np.ndarray,
    stub: int = 4,
    closure: str = "funnel",
    target_rg: float | None = None,
    funnel_excursion: float = 2.5,
    landing: int = 0,
    closure_iterations: int = 16,
) -> InteriorBatchResult:
    """Generate interior (connecting) IDR traces bridging two fixed anchors, via the two-stub setup.

    Residue 0 sits one bond from ``n_anchor`` and residue ``n_residues - 1`` one bond from
    ``c_anchor``. A short ``stub`` is grown from each anchor with its near-anchor steps skewed AWAY
    from that domain (``n_away`` / ``c_away``, the unit direction from the domain's centre of mass
    toward the anchor) so it leaves the crowded shell into the gap. The middle then bridges the two
    stub ends in the open space between the domains, by one of two closures:

    * ``"funnel"`` (default) grows the middle inside a shrinking reachability sphere
      (:func:`_bridge_funnel`), so it closes onto the far stub while keeping every pseudo-angle in
      the measured window; ``funnel_excursion`` adds a per-conformer perpendicular bulge that
      spreads a slack chain to its natural size so it does not collapse and self-clash, and
      ``target_rg`` (e.g. the ALBATROSS-predicted Rg) sizes that bulge automatically to a target
      radius of gyration instead of the fixed ``funnel_excursion``;
    * ``"fabrik"`` grows the middle free and projects it onto the two ends
      (:func:`close_chain_ends`) -- exact bonds, but the angle window is not preserved.

    The skew is a per-conformer draw, so the batch stays diverse. Returns an
    :class:`InteriorBatchResult`: the coordinates plus a per-conformer ``closure_feasible`` mask
    (a chain whose landing failed carries a long junction bond, which a caller MUST filter on --
    a broken chain can still show an on-target Rg), the achieved ``radius_of_gyration``, and whether
    the Rg target had to be clamped. Self-avoidance and obstacle clearance stay the caller's to
    filter with :func:`self_avoiding_mask` and :func:`clears_obstacles_mask`.
    """
    bond = CA_CA_BOND_LENGTH
    if n_residues < 2 * stub + 2:
        raise EngineError(
            f"n_residues={n_residues} is too short for two stubs of {stub}; "
            f"need at least {2 * stub + 2}."
        )
    if stub < 2:
        raise EngineError(f"stub must be at least 2 to seed the funnel direction, got {stub}.")
    n_anchor = np.asarray(n_anchor, dtype=np.float64)
    c_anchor = np.asarray(c_anchor, dtype=np.float64)
    n_bias = np.tile(_unit_rows(np.asarray(n_away, dtype=np.float64)), (n_conformers, 1))
    c_bias = np.tile(_unit_rows(np.asarray(c_away, dtype=np.float64)), (n_conformers, 1))

    # Each stub is the anchor plus `stub` IDR residues; translate so point 0 is the fixed anchor.
    sa = grow_free_batch(n_conformers, stub + 1, rng, bias_directions=n_bias, bias_residues=stub)
    stub_n = (sa - sa[:, 0:1] + n_anchor)[:, 1:]  # (B, stub, 3): IDR residues 0..stub-1
    sb = grow_free_batch(n_conformers, stub + 1, rng, bias_directions=c_bias, bias_residues=stub)
    # Grown from the C-anchor, so reverse it to end at residue n-1.
    stub_c = (sb - sb[:, 0:1] + c_anchor)[:, 1:][:, ::-1, :]  # (B, stub, 3): residues n-stub..n-1

    # Bridge the two stub ends with the middle residues, closed in open space.
    p = stub_n[:, -1]  # IDR residue stub-1
    q = stub_c[:, 0]  # IDR residue n-stub
    n_mid = n_residues - 2 * stub
    excursion = funnel_excursion
    rg_clamped = False
    if closure == "funnel":
        u_start = _unit_rows(stub_n[:, -1] - stub_n[:, -2])
        if target_rg is not None:
            excursion, rg_clamped = _excursion_for_target_rg(
                p, q, u_start, stub_n, stub_c, n_mid, rng, target_rg
            )
        middle = _bridge_funnel(p, q, u_start, n_mid, rng, excursion=excursion, landing=landing)
    elif closure == "fabrik":
        mid = grow_free_batch(n_conformers, n_mid + 2, rng)
        mid = mid - mid[:, 0:1] + p[:, None, :]
        mid = close_chain_ends(mid, p, q, bond=bond, iterations=closure_iterations)
        middle = mid[:, 1:-1]  # (B, n_mid, 3): IDR residues stub..n-stub-1
    else:
        raise EngineError(f"Unknown closure {closure!r}; use 'funnel' or 'fabrik'.")

    coords = np.concatenate([stub_n, middle, stub_c], axis=1)
    bond_error = np.abs(np.linalg.norm(np.diff(coords, axis=1), axis=-1) - bond).max(axis=1)
    return InteriorBatchResult(
        coords=coords,
        closure_feasible=bond_error < _CLOSURE_BOND_TOLERANCE,
        radius_of_gyration=radii_of_gyration(coords),
        target_rg=target_rg,
        rg_clamped=rg_clamped,
        excursion=float(excursion),
    )
