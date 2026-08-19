"""Rigid-body rotations, with every degenerate case handled explicitly.

Everything in DODO that moves an atom goes through this module.

Conventions
-----------
* A rotation is a ``(3, 3)`` matrix acting on **column** vectors: ``v' = R @ v``. This
  matches :meth:`dodo.structure.Domain.rotate` and scipy's
  :meth:`~scipy.spatial.transform.Rotation.as_matrix`.
* Coordinate arrays are ``(n, 3)`` **row** vectors, so applying a rotation to a whole
  array is ``coords @ R.T``. :func:`apply` does this so no caller has to remember which
  side the transpose goes on -- getting it backwards yields the inverse rotation, which
  is the third confirmed defect in the code this replaces.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial.transform import Rotation

from ..exceptions import GeometryError

__all__ = [
    "align_frame",
    "apply",
    "rotation_between_vectors",
    "rotation_between_vectors_batch",
    "rotation_from_axis_angle",
]

#: Vector norm at or below which a vector is treated as having no direction.
#:
#: Absolute rather than relative, which is safe because DODO's coordinates are in
#: Angstroms with magnitudes of order 1-1000: any vector shorter than this is either
#: exactly zero or the difference of two coordinates that round-tripped through a PDB
#: file as the same number.
ZERO_LENGTH_TOLERANCE = 1e-12

#: ``|sin(angle)|`` between two unit vectors below which they count as parallel or
#: antiparallel and the general axis-angle path is abandoned.
#:
#: The general path normalizes ``cross(a, b)``, whose norm is ``|sin(angle)|``, so as
#: the vectors approach (anti)parallel that division amplifies rounding error without
#: limit. Below this threshold the two special-case branches are used instead. They are
#: accurate where it matters: the parallel branch returns the identity, which differs
#: from the true rotation by at most 1e-8 radians, and for the antiparallel branch the
#: rotation axis is genuinely unconstrained, so any perpendicular axis is as correct as
#: the noisy one the general path would compute.
_PARALLEL_SIN_TOLERANCE = 1e-8

#: Tolerance on the ``det(R) == +1`` and ``R @ R.T == I`` self-checks.
#:
#: Generous relative to the 1e-16 of a well-conditioned 3x3 product, because the point
#: of the check is to catch a *wrong* matrix (a reflection, a scaling, an uninitialized
#: array), not to police floating-point noise.
_ROTATION_TOLERANCE = 1e-9

#: Module-level identity, so the orthonormality self-check does not rebuild one per call.
#: Read-only because it is shared: a caller-visible function that mutated it would break
#: every later check.
_IDENTITY = np.eye(3)
_IDENTITY.setflags(write=False)


def _as_vector3(vector: np.ndarray | tuple[float, float, float], name: str) -> np.ndarray:
    """Coerce ``vector`` to a finite ``(3,)`` float64 array, or raise."""
    array = np.asarray(vector, dtype=np.float64)
    if array.shape != (3,):
        raise GeometryError(f"{name} must be a 3-vector, got shape {array.shape}.")
    if not np.all(np.isfinite(array)):
        raise GeometryError(
            f"{name} contains non-finite values ({array.tolist()}). NaN or infinity here "
            f"means an earlier step already failed; propagating it would write NaN "
            f"coordinates into the output structure."
        )
    return array


def _unit(vector: np.ndarray | tuple[float, float, float], name: str) -> np.ndarray:
    """Return ``vector`` normalized to unit length, raising if it has no direction."""
    array = _as_vector3(vector, name)
    norm = float(np.sqrt(array @ array))
    if norm <= ZERO_LENGTH_TOLERANCE:
        raise GeometryError(
            f"{name} has zero length ({norm:g}) and therefore no direction. This usually "
            f"means two coordinates that should be distinct are identical."
        )
    return array / norm


def _determinant(matrix: np.ndarray) -> float:
    """Compute the determinant of a ``(3, 3)`` matrix by cofactor expansion.

    Hand-expanded rather than :func:`numpy.linalg.det` because this runs on the
    per-residue candidate-generation path, where a LAPACK call and its argument
    marshalling cost several times the eighteen multiplications below.
    """
    m = matrix
    return float(
        m[0, 0] * (m[1, 1] * m[2, 2] - m[1, 2] * m[2, 1])
        - m[0, 1] * (m[1, 0] * m[2, 2] - m[1, 2] * m[2, 0])
        + m[0, 2] * (m[1, 0] * m[2, 1] - m[1, 1] * m[2, 0])
    )


def _require_rotation_shape(rotation: np.ndarray, context: str) -> np.ndarray:
    """Coerce ``rotation`` to a ``(3, 3)`` float64 array, or raise."""
    matrix = np.asarray(rotation, dtype=np.float64)
    if matrix.shape != (3, 3):
        raise GeometryError(f"{context}: rotation must be a (3, 3) matrix, got {matrix.shape}.")
    return matrix


def _require_proper_rotation(rotation: np.ndarray, context: str) -> np.ndarray:
    """Check ``rotation`` is a proper rotation and return it as float64.

    Raises
    ------
    GeometryError
        If the matrix is not ``(3, 3)``, or its determinant is not +1, or it is not
        orthonormal.

    Notes
    -----
    Two independent things are checked, because neither implies the other.

    A determinant of -1 is a reflection: such a matrix maps ``a`` onto ``b`` just as well
    as a rotation does while mirroring everything else, so it produces a structure of the
    correct dimensions and the wrong chirality, and nothing downstream would notice.

    Orthonormality is what makes the map *rigid*, and the determinant does not imply it.
    A shear such as ``[[1, 0.6, 0], [0, 1, 0], [0, 0, 1]]`` has determinant exactly +1
    and is not a rotation: measured, it turns a trace of uniform 4.294 A CA-CA bonds into
    bonds of 5.385 / 3.280 / 3.280 / 5.385 A. This check used to be skipped by
    :func:`apply` as a hot-path saving, on the argument that no code path here can produce
    a shear -- true, but :func:`apply` exists to accept matrices from callers, and a
    stretched backbone is exactly the class of silent corruption this module is for.
    """
    matrix = _require_rotation_shape(rotation, context)
    determinant = _determinant(matrix)
    if abs(determinant - 1.0) > _ROTATION_TOLERANCE:
        kind = "an improper rotation (reflection)" if determinant < 0 else "not a rotation"
        raise GeometryError(
            f"{context}: matrix has determinant {determinant:.12g}, so it is {kind}. "
            f"A proper rotation has determinant +1; anything else changes handedness or "
            f"scale and would silently mirror or resize the structure."
        )
    deviation = float(np.abs(matrix @ matrix.T - _IDENTITY).max())
    if deviation > _ROTATION_TOLERANCE:
        raise GeometryError(
            f"{context}: matrix has determinant +1 but is not orthonormal "
            f"(max |R R^T - I| = {deviation:.3g}), so it is a shear or a "
            f"volume-preserving stretch rather than a rotation. Applying it would change "
            f"bond lengths and angles while leaving the overall volume intact."
        )
    return matrix


def _perpendicular_to(unit_axis: np.ndarray) -> np.ndarray:
    """Return some unit vector perpendicular to the unit vector ``unit_axis``."""
    smallest = int(np.argmin(np.abs(unit_axis)))
    reference = np.zeros(3)
    reference[smallest] = 1.0
    perpendicular = np.cross(unit_axis, reference)
    unit_perpendicular: np.ndarray = perpendicular / np.sqrt(perpendicular @ perpendicular)
    return unit_perpendicular


def rotation_from_axis_angle(
    axis: np.ndarray | tuple[float, float, float], angle_radians: float
) -> np.ndarray:
    """Rotation matrix for a right-handed rotation about ``axis``.

    Parameters
    ----------
    axis
        Rotation axis as a 3-vector. Need not be normalized; it is normalized here.
    angle_radians
        Rotation angle in **radians**, counter-clockwise about ``axis`` looking down it
        toward the origin. Radians, not degrees -- the parameter is named for it because
        the pre-rewrite code mixed the two units in adjacent functions.

    Returns
    -------
    np.ndarray
        A ``(3, 3)`` proper rotation matrix acting on column vectors.

    Raises
    ------
    GeometryError
        If ``axis`` is not a 3-vector, has zero length, or either argument is
        non-finite.

    Examples
    --------
    >>> import numpy as np
    >>> matrix = rotation_from_axis_angle([0.0, 0.0, 1.0], np.pi / 2)
    >>> np.allclose(matrix @ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0])
    True
    """
    unit_axis = _unit(axis, "axis")
    angle = float(angle_radians)
    if not np.isfinite(angle):
        raise GeometryError(f"angle_radians must be finite, got {angle}.")
    # scipy's exponential map is numerically careful about small angles (it uses the
    # series expansion of sin(x)/x rather than dividing by a tiny number). Hand-rolled
    # Rodrigues formulas that do not are one of the reasons this module exists.
    matrix = np.asarray(Rotation.from_rotvec(unit_axis * angle).as_matrix(), dtype=np.float64)
    return _require_proper_rotation(matrix, "rotation_from_axis_angle")


def rotation_between_vectors(
    a: np.ndarray | tuple[float, float, float], b: np.ndarray | tuple[float, float, float]
) -> np.ndarray:
    """Rotation matrix taking the direction of ``a`` onto the direction of ``b``.

    The minimal such rotation: about the axis ``a x b``, by the angle between them. Only
    directions matter, so the magnitudes of ``a`` and ``b`` are irrelevant.

    Parameters
    ----------
    a
        Source 3-vector.
    b
        Target 3-vector.

    Returns
    -------
    np.ndarray
        A ``(3, 3)`` proper rotation ``R`` with ``R @ unit(a) == unit(b)``.

    Raises
    ------
    GeometryError
        If either vector is not a finite 3-vector, or has zero length.

    Notes
    -----
    The antiparallel case (``b`` pointing opposite to ``a``) is the one that matters. A
    rotation by pi about *any* axis perpendicular to ``a`` maps ``a`` to ``-a``; which
    axis is chosen is genuinely arbitrary, so one is chosen deterministically. What must
    not happen is returning ``-I``, which is what the pre-rewrite code did: that is a
    point inversion with determinant -1, and it mirrors the structure.

    Examples
    --------
    >>> import numpy as np
    >>> a = np.array([0.0, 0.0, 1.0])
    >>> matrix = rotation_between_vectors(a, -a)
    >>> np.allclose(matrix @ a, -a), round(float(np.linalg.det(matrix)), 12)
    (True, 1.0)
    """
    unit_a = _unit(a, "a")
    unit_b = _unit(b, "b")

    cross = np.cross(unit_a, unit_b)
    sine = float(np.sqrt(cross @ cross))
    cosine = float(np.clip(unit_a @ unit_b, -1.0, 1.0))

    if sine <= _PARALLEL_SIN_TOLERANCE:
        if cosine > 0.0:
            return np.eye(3)
        return rotation_from_axis_angle(_perpendicular_to(unit_a), np.pi)

    # atan2 rather than arccos(cosine): arccos loses precision as the angle approaches
    # 0 or pi, exactly where this function is most often called from a builder aligning
    # a template that is already nearly in place.
    return rotation_from_axis_angle(cross / sine, float(np.arctan2(sine, cosine)))


def rotation_between_vectors_batch(
    a: np.ndarray | tuple[float, float, float],
    b: np.ndarray,
) -> np.ndarray:
    """Per-row :func:`rotation_between_vectors` for a whole batch, in one vectorized pass.

    Parameters
    ----------
    a
        Source direction: a single 3-vector broadcast to every row, or an ``(n, 3)`` array
        of per-row sources.
    b
        ``(n, 3)`` target directions.

    Returns
    -------
    np.ndarray
        ``(n, 3, 3)`` proper rotations. Row ``i`` equals ``rotation_between_vectors(a_i,
        b[i])`` bit for bit.

    Raises
    ------
    GeometryError
        If ``b`` is not ``(n, 3)``, if ``a`` is neither ``(3,)`` nor ``(n, 3)``, if either
        contains non-finite values, or if any row of either has zero length.

    Notes
    -----
    The same three-branch logic as the scalar function, applied row-wise: the general
    axis-angle path runs the whole batch through one :class:`scipy.spatial.transform.Rotation`
    call, the (near-)parallel rows take the identity, and the antiparallel rows take a pi
    rotation about an arbitrary perpendicular. Only the antiparallel branch loops, because
    its axis depends on the row's own source direction; it is empty for essentially every
    real input (it needs a target within 1e-8 rad of *opposite* the source), so the loop is
    a formality that keeps the degenerate rows identical to the scalar choice rather than a
    cost. Written so that the arithmetic on each row matches the scalar path term for term,
    which is verified against it directly in the geometry tests.
    """
    b_arr = np.asarray(b, dtype=np.float64)
    if b_arr.ndim != 2 or b_arr.shape[1] != 3:
        raise GeometryError(f"b must have shape (n, 3), got {b_arr.shape}.")
    n = b_arr.shape[0]
    a_arr = np.asarray(a, dtype=np.float64)
    if a_arr.shape == (3,):
        a_arr = np.broadcast_to(a_arr, (n, 3))
    elif a_arr.shape != (n, 3):
        raise GeometryError(f"a must have shape (3,) or {(n, 3)}, got {a_arr.shape}.")
    if not (np.all(np.isfinite(a_arr)) and np.all(np.isfinite(b_arr))):
        raise GeometryError("rotation_between_vectors_batch received non-finite input.")

    a_norm = np.sqrt(np.sum(a_arr * a_arr, axis=1))
    b_norm = np.sqrt(np.sum(b_arr * b_arr, axis=1))
    if np.any(a_norm <= ZERO_LENGTH_TOLERANCE) or np.any(b_norm <= ZERO_LENGTH_TOLERANCE):
        raise GeometryError(
            "rotation_between_vectors_batch received a zero-length vector, so its direction "
            "is undefined: two coordinates that should be distinct are identical."
        )
    unit_a = a_arr / a_norm[:, None]
    unit_b = b_arr / b_norm[:, None]

    cross = np.cross(unit_a, unit_b)
    sine = np.sqrt(np.sum(cross * cross, axis=1))
    cosine = np.clip(np.sum(unit_a * unit_b, axis=1), -1.0, 1.0)

    result = np.empty((n, 3, 3), dtype=np.float64)
    general = sine > _PARALLEL_SIN_TOLERANCE
    if np.any(general):
        # cross / sine, then re-normalized: the scalar path divides by sine and then hands
        # the result to rotation_from_axis_angle, which normalizes again. Both steps are
        # reproduced so the rotation vector matches to the last bit.
        raw = cross[general] / sine[general][:, None]
        raw_norm = np.sqrt(np.sum(raw * raw, axis=1))
        unit_axis = raw / raw_norm[:, None]
        angle = np.arctan2(sine[general], cosine[general])
        result[general] = Rotation.from_rotvec(unit_axis * angle[:, None]).as_matrix()

    degenerate = ~general
    if np.any(degenerate):
        result[degenerate & (cosine > 0.0)] = _IDENTITY
        for idx in np.flatnonzero(degenerate & (cosine <= 0.0)):
            result[idx] = rotation_from_axis_angle(_perpendicular_to(unit_a[idx]), np.pi)
    return result


def align_frame(
    origin: np.ndarray | tuple[float, float, float],
    direction: np.ndarray | tuple[float, float, float],
) -> np.ndarray:
    """Orthonormal right-handed frame anchored at ``origin`` with z pointing at ``direction``.

    Parameters
    ----------
    origin
        The point the frame is anchored at.
    direction
        The point the frame's z axis points **toward**, not a direction vector: the z axis
        is the unit vector along ``direction - origin``. Passing ``origin=(0, 0, 0)`` is
        therefore how you use a plain direction vector, and the two readings coincide
        there. This is the interpretation that makes both arguments meaningful --
        ``align_frame(ca[i], ca[i + 1])`` builds the frame of residue ``i`` pointing along
        the chain -- but it does mean that pairing a non-zero ``origin`` with a
        ``direction`` you meant as a vector silently gives the wrong frame, and nothing can
        detect that. Pass two points, or pass a vector from the origin.

    Returns
    -------
    np.ndarray
        A ``(3, 3)`` proper rotation whose **columns** are the frame's x, y and z axes
        in world coordinates. So ``frame @ (0, 0, 1) == unit(direction - origin)``, and
        ``apply(local_coords, frame) + origin`` maps coordinates expressed in the frame
        into world coordinates, while ``frame.T`` maps world directions into the frame.

    Raises
    ------
    GeometryError
        If either argument is not a finite 3-vector, or if ``direction`` coincides with
        ``origin``, which leaves the z axis undefined.

    Notes
    -----
    The x and y axes are not determined by the input -- only z is -- so they are chosen
    deterministically from the world axis least aligned with z. Two calls with the same
    arguments therefore return the same frame, which is what makes any downstream
    sampling reproducible. They are *not* continuous in ``direction``: a small change in
    direction can swing x and y around the z axis. Do not use this frame to interpolate
    between orientations; use :func:`rotation_between_vectors` or scipy's ``Slerp``.
    """
    start = _as_vector3(origin, "origin")
    end = _as_vector3(direction, "direction")
    offset = end - start
    norm = float(np.sqrt(offset @ offset))
    if norm <= ZERO_LENGTH_TOLERANCE:
        raise GeometryError(
            f"direction {end.tolist()} coincides with origin {start.tolist()} (separation "
            f"{norm:g}), so the frame's z axis is undefined."
        )
    z_axis = offset / norm

    smallest = int(np.argmin(np.abs(z_axis)))
    reference = np.zeros(3)
    reference[smallest] = 1.0
    x_axis = np.cross(reference, z_axis)
    x_axis /= np.sqrt(x_axis @ x_axis)
    # y = z x x rather than x x z: the former gives a right-handed frame (det +1), the
    # latter a left-handed one (det -1) that would mirror anything built in it.
    y_axis = np.cross(z_axis, x_axis)

    frame = np.column_stack((x_axis, y_axis, z_axis))
    return _require_proper_rotation(frame, "align_frame")


def apply(
    coords: np.ndarray,
    rotation: np.ndarray,
    about: np.ndarray | tuple[float, float, float] | None = None,
) -> np.ndarray:
    """Rotate coordinates, optionally about a fixed point.

    Parameters
    ----------
    coords
        ``(n, 3)`` coordinates, or a single ``(3,)`` point. The result has the same
        shape as the input.
    rotation
        A ``(3, 3)`` proper rotation matrix, acting on column vectors.
    about
        Point held fixed by the rotation. Defaults to the world origin.


    Returns
    -------
    np.ndarray
        A new array of rotated coordinates. The input is never modified.

    Raises
    ------
    GeometryError
        If ``coords`` is not ``(n, 3)`` or ``(3,)``, if it contains non-finite values,
        or if ``rotation`` is not a proper rotation matrix.

    Notes
    -----
    Non-finite input coordinates raise rather than passing through. A NaN that reaches
    here came from a failed sampler upstream, and a rotation is the last place it can be
    caught before it becomes an ``ATOM`` record reading ``nan nan nan``, which is what
    the pre-rewrite pipeline emitted.

    """
    array = np.asarray(coords, dtype=np.float64)
    single_point = array.ndim == 1
    if single_point:
        array = array.reshape(1, 3) if array.shape == (3,) else array
    if array.ndim != 2 or array.shape[1] != 3:
        raise GeometryError(f"coords must have shape (n, 3) or (3,), got {np.shape(coords)}.")
    if not np.all(np.isfinite(array)):
        bad = int(np.count_nonzero(~np.all(np.isfinite(array), axis=1)))
        raise GeometryError(
            f"coords contains non-finite values in {bad} of {array.shape[0]} rows. "
            f"Rotating NaN produces NaN; the sampler that generated these should have "
            f"raised instead of returning them."
        )
    # Full rigidity check, determinant and orthonormality. Without the second half a
    # volume-preserving shear passes and the "rotated" coordinates come back with
    # different bond lengths and angles than they went in with.
    matrix = _require_proper_rotation(rotation, "apply")

    if about is None:
        rotated: np.ndarray = array @ matrix.T
    else:
        centre = _as_vector3(about, "about")
        rotated = (array - centre) @ matrix.T + centre
    return rotated.reshape(3) if single_point else rotated
