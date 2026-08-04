"""Random and enumerated candidate positions for a growing CA trace.

Three primitives, one implementation each. That is a deliberate reaction to the code
being replaced, which had a scalar and a "vectorized" sibling for every sampler here:
of three such pairs, two of the batch versions were outright broken (one applied a
forward rotation where the inverse was required, one divided by a zero-length normal and
returned all-NaN) and all three were dead code, because every caller used the scalar
version in a Python loop. The batch form is the only form here, and ``n = 1`` is a
normal case of it.

Reproducibility
---------------
Every stochastic function takes an explicit :class:`numpy.random.Generator`. There is no
module-level generator, no use of the legacy global ``numpy.random`` state and no use of
the standard library's ``random``. The pre-rewrite sampler called ``np.random.normal``
directly, which made a stochastic structure builder impossible to seed and therefore
impossible to regression-test: two runs of the same command produced different
structures with no way to reproduce either.
"""

from __future__ import annotations

from functools import lru_cache
from typing import NamedTuple

import numpy as np

from ..constants import (
    CA_CA_BOND_LENGTH,
    CANDIDATES_PER_ANGLE,
    backbone_angle_grid,
)
from ..exceptions import GeometryError
from .transforms import ZERO_LENGTH_TOLERANCE, _as_vector3, apply, rotation_between_vectors

__all__ = [
    "CacheStatistics",
    "cone_candidates",
    "cone_template_cache_info",
    "random_points_on_sphere",
    "random_unit_vectors",
]

#: Axis the cone template is built around. Any fixed choice works; the template is
#: rotated onto the real chain axis per call.
_TEMPLATE_AXIS = np.array([0.0, 0.0, 1.0])

#: Golden-ratio conjugate, used to offset each angle ring's azimuth.
#:
#: Without an offset every ring uses the same ``per_angle`` azimuths, so all 355 default
#: candidates lie in just 10 meridional half-planes: a single obstacle sitting in one of
#: those planes rules out 71 candidates at once, and the search has far less freedom than
#: its candidate count suggests. Advancing the azimuth by the golden angle per ring
#: spreads the union of the rings evenly around the axis. It is a fixed irrational
#: offset, not a random one, so :func:`cone_candidates` stays deterministic.
_GOLDEN_RATIO_CONJUGATE = 0.5 * (np.sqrt(5.0) - 1.0)

#: Cache size for cone templates. Small on purpose: a build uses one or two distinct
#: parameter combinations (the default grid, plus whatever a closure step narrows it to),
#: so a handful of slots gives a ~100% hit rate.
_TEMPLATE_CACHE_SIZE = 32


def _require_generator(rng: np.random.Generator) -> np.random.Generator:
    """Check that ``rng`` really is a Generator, with an actionable message if not.

    Passing a seed integer or the legacy ``numpy.random.RandomState`` would otherwise
    fail deep inside the function with an ``AttributeError`` about a method nobody asked
    about.
    """
    if not isinstance(rng, np.random.Generator):
        raise TypeError(
            f"rng must be a numpy.random.Generator, got {type(rng).__name__}. "
            f"Construct one with numpy.random.default_rng(seed); DODO takes an explicit "
            f"generator everywhere so that any build can be reproduced exactly."
        )
    return rng


def _require_count(value: object, name: str, error: type[Exception]) -> int:
    """Return ``value`` as an ``int``, refusing anything not already integral.

    ``int(2.9)`` is ``2``, so a float count used to be accepted and silently truncated:
    ``random_unit_vectors(2.9)`` returned 2 vectors and ``cone_candidates(per_angle=2.9)``
    returned 2 candidates. A caller that derives a count arithmetically -- from a length,
    a ratio, a budget divided by a grid size -- then samples a sparser cone or draws fewer
    directions than it asked for, and nothing in the returned array says so. In a
    stochastic builder that shows up as an unexplained failure to place a residue.

    numpy integers are accepted because they are integers, and ``bool`` is accepted
    because it is an ``int`` subclass and ``n=True`` is an unambiguous request for one.
    ``error`` is the exception type to raise, so each caller can keep its own contract
    (``TypeError`` for a plain argument-type mistake, :class:`GeometryError` where the rest
    of the geometry validation raises).
    """
    if isinstance(value, (int, np.integer)):
        return int(value)
    raise error(
        f"{name} must be an integer, got {value!r} of type {type(value).__name__}. "
        f"A fractional count has no meaning here, and truncating it would quietly return "
        f"fewer items than requested."
    )


def random_unit_vectors(n: int, rng: np.random.Generator) -> np.ndarray:
    """Draw ``n`` directions uniformly distributed on the unit sphere.

    Parameters
    ----------
    n
        Number of vectors. Must be an integer; a fractional count is refused rather than
        truncated. ``0`` is allowed and returns an empty ``(0, 3)`` array, which keeps
        vectorized callers from needing a special case.
    rng
        Seeded random generator.

    Returns
    -------
    np.ndarray
        ``(n, 3)`` unit vectors.

    Raises
    ------
    ValueError
        If ``n`` is negative.
    TypeError
        If ``rng`` is not a :class:`numpy.random.Generator`, or ``n`` is not an integer.
    GeometryError
        If a drawn vector cannot be normalized even after resampling, which would
        require repeatedly drawing three Gaussians that all underflow to zero.

    Notes
    -----
    Uses normalized Gaussians. This is the standard correct method: the trivariate
    standard normal is spherically symmetric, so dividing by its norm gives an exactly
    uniform direction.

    The tempting alternative -- sampling the polar angle ``theta`` uniformly on
    ``[0, pi]`` and the azimuth uniformly on ``[0, 2 pi)`` -- is **wrong**, and wrong in
    a way that matters here. The area element on a sphere is ``sin(theta) d(theta)
    d(phi)``, so uniform ``theta`` oversamples the poles by ``1 / sin(theta)``: the
    density diverges at both ends. A random-walk step drawn that way is biased toward
    continuing along, or reversing, the previous axis, which changes the chain's
    persistence length and therefore its radius of gyration -- the one quantity DODO
    exists to get right. Sampling ``cos(theta)`` uniformly on ``[-1, 1]`` is the correct
    inverse-CDF version of that approach; normalized Gaussians are simpler and avoid the
    trigonometry entirely.
    """
    _require_generator(rng)
    count = _require_count(n, "n", TypeError)
    if count < 0:
        raise ValueError(f"n must be non-negative, got {n}.")
    if count == 0:
        return np.zeros((0, 3), dtype=np.float64)

    vectors = rng.standard_normal((count, 3))
    norms = np.linalg.norm(vectors, axis=1, keepdims=True)

    # A Gaussian triple whose norm underflows is astronomically unlikely (it needs all
    # three components below ~1e-160), but "astronomically unlikely" is not "impossible",
    # and the failure mode is a division by zero that yields NaN coordinates. Resample
    # the offending rows instead, and raise rather than return NaN if that somehow fails.
    for _ in range(8):
        degenerate = (norms <= 0.0).ravel()
        if not degenerate.any():
            break
        vectors[degenerate] = rng.standard_normal((int(degenerate.sum()), 3))
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
    else:
        # Reached only if the eighth resample still left a zero-norm row. Checked again
        # rather than assumed, so a successful final resample is not reported as failure.
        if (norms <= 0.0).any():
            raise GeometryError(
                "Could not draw non-degenerate random directions after 8 attempts. This "
                "should be impossible; suspect a broken or exhausted generator."
            )

    result: np.ndarray = vectors / norms
    return result


def random_points_on_sphere(
    centre: np.ndarray | tuple[float, float, float],
    radius: float,
    n: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw ``n`` points uniformly distributed on a sphere's surface.

    Parameters
    ----------
    centre
        Sphere centre, a 3-vector.
    radius
        Sphere radius in Angstroms. Must be strictly positive: a zero radius would
        return ``n`` copies of ``centre``, which is a degenerate answer of the kind that
        makes a caller's later failure hard to trace back to here.
    n
        Number of points. ``0`` returns an empty ``(0, 3)`` array.
    rng
        Seeded random generator.

    Returns
    -------
    np.ndarray
        ``(n, 3)`` points, each exactly ``radius`` from ``centre``.

    Raises
    ------
    GeometryError
        If ``centre`` is not a finite 3-vector, or ``radius`` is not positive and finite.
    """
    array = np.asarray(centre, dtype=np.float64)
    if array.shape != (3,):
        raise GeometryError(f"centre must be a 3-vector, got shape {array.shape}.")
    if not np.all(np.isfinite(array)):
        raise GeometryError(f"centre contains non-finite values ({array.tolist()}).")
    value = float(radius)
    if not np.isfinite(value) or value <= 0.0:
        raise GeometryError(f"radius must be positive and finite, got {radius!r}.")
    return array + value * random_unit_vectors(n, rng)


@lru_cache(maxsize=_TEMPLATE_CACHE_SIZE)
def _cone_template(angles: tuple[float, ...], per_angle: int, bond_length: float) -> np.ndarray:
    """Unit-frame candidate positions for the next CA, with the chain axis along +z.

    The returned array is ``(len(angles) * per_angle, 3)``, ordered angle-block-major:
    rows ``[i * per_angle, (i + 1) * per_angle)`` are the ring for ``angles[i]``. It is
    read-only *and owns its data*, so a caller cannot mutate the cached object and poison
    every later call. Both halves matter: the first version of this returned a read-only
    ``reshape`` view of a writable ``(n_angles, per_angle, 3)`` buffer, and
    ``setflags(write=False)`` on a view says nothing about its base. Writing zeros through
    ``template.base`` made every later :func:`cone_candidates` call for that cache key
    return all candidates sitting on the apex -- zero-length bonds out of the one function
    that guarantees exact bond lengths.

    Cached because the geometry depends only on these three arguments and not at all on
    the coordinates it will be rotated onto -- the pre-rewrite code noted that opportunity
    and never took it, rebuilding every candidate offset from scratch for every residue of every
    attempt, at a measured 2.82 ms per call.

    Measured here, when the default was 71 angles x 10 candidates (it is now x 5, so a call is
    cheaper than these figures): a full
    :func:`cone_candidates` call costs 44 us warm and 58 us with the cache cleared, so the
    cache is worth 1.3x, not the 86x the pre-rewrite note anticipated. The reason is that
    most of the old 2.82 ms was not the trigonometry but the Python loop around it -- 71
    angles, each calling a circle-point helper and then transforming its 10 points one at
    a time. Building the whole template with array operations removes that, and what the
    cache saves on top is the remaining 13 us of arithmetic. It stays because 13 us of 44
    is still 30% of the per-residue budget and the output is bit-identical either way.
    """
    angle_array = np.asarray(angles, dtype=np.float64)

    # The CA(i-1)-CA(i)-CA(i+1) angle is measured at CA(i) between the vector back to
    # CA(i-1) and the vector out to CA(i+1). The template axis +z points *forward* along
    # CA(i-1) -> CA(i), so the polar angle of the new bond measured from +z is the
    # supplement of the pseudo-angle. Getting this supplement backwards yields candidates
    # whose measured angles are 180 - requested: a plausible-looking backbone built from
    # systematically wrong geometry.
    polar = np.deg2rad(180.0 - angle_array)
    axial = bond_length * np.cos(polar)  # (n_angles,) component along +z
    radial = bond_length * np.sin(polar)  # (n_angles,) distance from the +z axis

    ring = 2.0 * np.pi * np.arange(per_angle, dtype=np.float64) / per_angle
    offset = 2.0 * np.pi * _GOLDEN_RATIO_CONJUGATE * np.arange(angle_array.size, dtype=np.float64)
    azimuth = ring[None, :] + offset[:, None]  # (n_angles, per_angle)

    # Allocate the flat array the caller will receive and write through a temporary
    # ring-shaped view of it, rather than allocating the ring shape and returning a
    # reshaped view. The array handed out is then the base, so making it read-only makes
    # the whole cached buffer read-only: there is no writable ``.base`` to reach through.
    flat = np.empty((angle_array.size * per_angle, 3), dtype=np.float64)
    rings = flat.reshape(angle_array.size, per_angle, 3)
    rings[:, :, 0] = radial[:, None] * np.cos(azimuth)
    rings[:, :, 1] = radial[:, None] * np.sin(azimuth)
    rings[:, :, 2] = axial[:, None]
    del rings  # the only writable handle on the cached buffer; drop it before returning

    flat.setflags(write=False)
    return flat


class CacheStatistics(NamedTuple):
    """Hit statistics for the cone template cache.

    The same four fields :meth:`functools.lru_cache.cache_info` reports, restated as a
    typed tuple so callers do not have to reach for a private typeshed name.
    """

    hits: int
    misses: int
    maxsize: int | None
    currsize: int


def cone_template_cache_info() -> CacheStatistics:
    """Return cache statistics for the cone template cache.

    Exposed for tests and profiling: a build showing one miss per residue is rebuilding
    the template every step instead of rotating a cached one.
    """
    info = _cone_template.cache_info()
    return CacheStatistics(info.hits, info.misses, info.maxsize, info.currsize)


def cone_candidates(
    previous_ca: np.ndarray | tuple[float, float, float],
    current_ca: np.ndarray | tuple[float, float, float],
    angles: np.ndarray | tuple[float, ...] | None = None,
    per_angle: int | None = None,
    bond_length: float | None = None,
) -> np.ndarray:
    """Candidate positions for the next CA, given the two most recent ones.

    This is the central primitive of the build engine. Every candidate sits exactly
    ``bond_length`` from ``current_ca`` and forms a CA-CA-CA pseudo-angle at
    ``current_ca`` drawn from ``angles``, so any candidate accepted by a caller is
    physically valid backbone geometry by construction -- as opposed to being checked
    afterwards, or (as in the code this replaces) not at all.

    Parameters
    ----------
    previous_ca
        Position of CA(i-1). Only the *direction* from here to ``current_ca`` is used, so
        this pair need not already be exactly one bond length apart.
    current_ca
        Position of CA(i), the apex of the cone.
    angles
        CA-CA-CA pseudo-angles to generate, in degrees. Defaults to
        :func:`dodo.constants.backbone_angle_grid`. **Order is preserved** in the output:
        the first ``per_angle`` rows are the ring for ``angles[0]``, and so on. The
        default grid is ordered ideal-first, so row index is a usable *preference order*
        over pseudo-angles (see "Choosing among the candidates" below for what to do with
        it, and what not to). The pre-rewrite vectorized generator emitted 91, 92, 93...
        and silently discarded the ordering.

        An angle of exactly 180 degrees is a degenerate ring: the cone radius is zero, so
        all ``per_angle`` rows for that angle coincide on the straight continuation of the
        previous bond. The shape contract is kept (the caller gets the rows it asked for),
        but that ring offers one distinct choice, not ``per_angle``.
    per_angle
        Candidate positions around each angle's ring. Must be an integer; a fractional
        value is refused rather than truncated. Defaults to
        :data:`dodo.constants.CANDIDATES_PER_ANGLE`.
    bond_length
        CA-CA distance in Angstroms. Defaults to
        :data:`dodo.constants.CA_CA_BOND_LENGTH`.

    Returns
    -------
    np.ndarray
        ``(len(angles) * per_angle, 3)`` candidate positions, writable and independent of
        the internal cache.

    Raises
    ------
    GeometryError
        If ``previous_ca`` and ``current_ca`` coincide (the cone axis would be
        undefined), if either is not a finite 3-vector, if ``angles`` is empty or
        contains a value outside ``(0, 180]``, if ``per_angle`` is not an integer or is
        below 1, or if ``bond_length`` is not positive and finite.

    Notes
    -----
    **Choosing among the candidates.** This function does not choose, takes no ``rng``,
    and the row order is deterministic. The caller's selection policy determines the
    chain's global dimensions, and the ordering here is *not* that policy.

    What the ordering is for: rows are grouped ideal-angle-first, so row index encodes a
    preference over pseudo-angles. Use it as a *weight* -- sample uniformly among the
    candidates that pass the caller's clash test, or sample with a probability that decays
    with distance from :data:`dodo.constants.BACKBONE_ANGLE_IDEAL`
    (:func:`dodo.engines.walk` does the latter). Either way the choice must be random.

    What it is not for: taking the first non-clashing candidate. That policy always picks
    the most ideal *angle* still available, and because the azimuths within a ring are a
    fixed sequence it also picks correlated *dihedrals*, so the chain marches off in one
    direction instead of exploring. Measured at n=200 against
    ``constants.flory_end_to_end(200) = 98.8 A``: first-non-clashing selection gives
    Re 334-641 A (3.4x to 6.5x too extended, Re/contour up to 0.85), while uniform
    selection over the identical candidate sets gives Re 41-162 A with a mean near the
    prediction. In the fully degenerate case -- always row 0, from a coplanar start -- the
    chain closes into a planar hexagon: 200 residues at Rg 4.11 A, every dihedral exactly
    0 degrees, thousands of coincident CAs.

    No local check catches this. Every candidate here has an exact bond length and an
    in-window pseudo-angle by construction, and the greedy chains above have no steric
    clashes either, so :func:`dodo.geometry.metrics.validate_ca_trace` passes them. An
    over-extended chain is only visible in a global observable (Re, Rg), which is why
    ``test_taking_the_first_candidate_is_not_a_selection_policy`` measures one.

    Examples
    --------
    >>> import numpy as np
    >>> from dodo.geometry.metrics import ca_pseudo_angles
    >>> candidates = cone_candidates([0.0, 0.0, 0.0], [3.8, 0.0, 0.0], angles=[120.0])
    >>> trace = np.vstack([[0.0, 0.0, 0.0], [3.8, 0.0, 0.0], candidates[0]])
    >>> bool(np.isclose(ca_pseudo_angles(trace)[0], 120.0))
    True
    """
    apex = _as_vector3(current_ca, "current_ca")
    axis = apex - _as_vector3(previous_ca, "previous_ca")
    if float(np.sqrt(axis @ axis)) <= ZERO_LENGTH_TOLERANCE:
        raise GeometryError(
            f"previous_ca and current_ca are the same point ({apex.tolist()}), so the cone "
            f"axis is undefined. Two consecutive CA positions must be distinct; if the "
            f"caller is starting a fresh chain it needs a seed direction, not a repeat."
        )

    resolved_angles = backbone_angle_grid() if angles is None else np.asarray(angles, dtype=float)
    resolved_angles = np.atleast_1d(resolved_angles)
    if resolved_angles.ndim != 1 or resolved_angles.size == 0:
        raise GeometryError(
            f"angles must be a non-empty 1-D sequence of degrees, got shape "
            f"{resolved_angles.shape}."
        )
    if not np.all(np.isfinite(resolved_angles)):
        raise GeometryError("angles contains non-finite values.")
    # 0 and 180 bracket the physically meaningful range for an angle between two bonds.
    # An angle of exactly 0 would place CA(i+1) on top of CA(i-1); outside the range the
    # supplement used below silently aliases onto a different angle, which is worse than
    # refusing.
    out_of_range = resolved_angles[(resolved_angles <= 0.0) | (resolved_angles > 180.0)]
    if out_of_range.size:
        raise GeometryError(
            f"angles must lie in (0, 180] degrees; got {out_of_range[:5].tolist()}. "
            f"DODO's own generation window is "
            f"{backbone_angle_grid().min():.0f}-{backbone_angle_grid().max():.0f} degrees."
        )

    count = (
        CANDIDATES_PER_ANGLE
        if per_angle is None
        else _require_count(per_angle, "per_angle", GeometryError)
    )
    if count < 1:
        raise GeometryError(f"per_angle must be at least 1, got {count}.")
    length = CA_CA_BOND_LENGTH if bond_length is None else float(bond_length)
    if not np.isfinite(length) or length <= 0.0:
        raise GeometryError(f"bond_length must be positive and finite, got {bond_length!r}.")

    # Resolve every default *before* the cache lookup so that cone_candidates(a, b) and
    # cone_candidates(a, b, backbone_angle_grid(), 10, 3.8) share one cache entry. tolist()
    # rather than a float() comprehension: it already yields Python floats, and at 71
    # angles the comprehension costs more than the cache lookup it feeds.
    template = _cone_template(tuple(resolved_angles.tolist()), count, length)

    # rotation_between_vectors handles the exactly-antiparallel case properly: a chain
    # running along -z relative to the template axis used to get -I from the pre-rewrite
    # code, which negated the axial component and turned every requested angle theta into
    # 180 - theta.
    rotation = rotation_between_vectors(_TEMPLATE_AXIS, axis)
    candidates: np.ndarray = apply(template, rotation) + apex
    return candidates
