"""Place backbone N, C and O on a CA-only trace, using four consecutive alpha carbons.

The advantage DODO has over general backbone prediction is that it already knows *every* alpha
carbon. That turns backbone placement from a modelling problem into a lookup, because of one fact
about protein geometry:

**The peptide unit is rigid.** Measured over 19,302 peptide units from 100 frames of all-atom IDR
simulation, its internal bonds do not move:

===============  =================  =============================
bond             measured           :mod:`dodo.constants` value
===============  =================  =============================
CA(i)-C(i)       1.525 +/- 0.004 A  ``CA_C_BOND_LENGTH`` 1.525
C(i)-N(i+1)      1.329 +/- 0.004 A  ``C_N_PEPTIDE_BOND_LENGTH`` 1.329
N(i+1)-CA(i+1)   1.458 +/- 0.004 A  ``N_CA_BOND_LENGTH`` 1.458
CA(i)-CA(i+1)    3.813 +/- 0.033 A  ``CA_CA_BOND_LENGTH`` 3.81
===============  =================  =============================

So given both flanking alpha carbons, the unit between them has exactly **one** degree of freedom:
a rotation about the CA-CA axis. Everything else is fixed. The whole problem is predicting that one
angle, and the fourth alpha carbon is what predicts it.

Why four and not three
----------------------
Three consecutive alpha carbons give the CA-CA-CA pseudo-*angle*. The peptide plane's rotation is
set by the pseudo-*dihedral*, which needs a fourth. Measured, held out on alternating frames:

==========  =========  =========  =========
predictor   C(i)       N(i+1)     O(i)
==========  =========  =========  =========
3 CAs       0.471 A    0.311 A    1.495 A
4 CAs       0.274 A    0.192 A    0.843 A
==========  =========  =========  =========

The fourth alpha carbon removes 53% of the peptide plane's angular uncertainty, from 64.4 to
30.5 degrees. It roughly halves the error on C and N.

Three things this module does *not* do, each for a measured reason
------------------------------------------------------------------
**It does not average positions.** Binning the local-frame coordinates and using the bin mean --
which is the obvious implementation -- produces C-N bonds of 1.254 A against an ideal 1.329, because
averaging points scattered on a sphere lands inside it. Every placed atom is instead rescaled onto
its exact ideal bond length, which costs nothing in accuracy (C 0.279 -> 0.281 A) and makes the
output valid by construction.

**It does not predict the carbonyl O.** O is fully determined by CA(i), C(i) and N(i+1): placed by
ideal geometry from the *true* three, it lands within **0.013 A**. Its 0.84 A error above is
entirely inherited from C and N, amplified because O sits 1.769 A from the CA-CA axis where C sits
0.556 A -- so the same angular uncertainty displaces it three times as far. Predicting it
separately is both redundant and worse.

**It does not invent alpha carbons at the termini.** Synthesizing a missing CA by reflection, or by
extrapolating at the mean pseudo-angle, was measured at 3.9-5.3 A error -- larger than the 3.81 A
bond itself, so a synthesized alpha carbon carries no information about the real one. Nor does the
adjacent peptide unit help: copying its plane rotation leaves a 99.0 degree residual, *worse* than
the 64.4 degree marginal, because successive peptide planes alternate rather than track. So a
terminal unit falls back to the marginal for its frame -- the honest 3-CA answer -- and the two
atoms no peptide unit covers are placed by ideal internal geometry.

Provenance
----------
The table is MEASURED from 100 frames of all-atom IDR simulation supplied by the author
(``subset_frames``), 19,302 peptide units, binned at 20 degrees. Every bin holds at least 146
observations. It is IDR data on purpose: this module exists to rebuild disordered regions.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass, replace
from typing import TYPE_CHECKING, Any, Final

import numpy as np

from ..constants import (
    C_N_PEPTIDE_BOND_LENGTH,
    C_O_BOND_LENGTH,
    CA_C_BOND_LENGTH,
    CA_C_O_ANGLE,
    N_CA_BOND_LENGTH,
    N_CA_C_ANGLE,
)
from ..exceptions import GeometryError

if TYPE_CHECKING:  # pragma: no cover
    from ..structure import Structure

__all__ = [
    "BackboneResult",
    "CABackboneResult",
    "SeamStrain",
    "add_backbone_to_rebuilt",
    "backbone_from_ca",
]

#: Width of a table bin, in degrees of CA pseudo-dihedral.
#:
#: MEASURED choice. Narrowing below 20 degrees stops helping -- held-out C error is 0.283 A at 30
#: degrees, 0.281 at 20, 0.281 at 10 and 0.281 at 5 -- because the residual is the peptide plane's
#: intrinsic spread, not binning coarseness. 20 keeps every bin well populated.
_BIN_WIDTH_DEGREES: Final[int] = 20

_C_BY_BIN: Final[tuple[tuple[float, float, float], ...]] = (
    (+1.4223, -0.2170, +0.4056),  # n=423
    (+1.4230, -0.2540, +0.3861),  # n=613
    (+1.4226, -0.2919, +0.3709),  # n=787
    (+1.4217, -0.3135, +0.2853),  # n=855
    (+1.4207, -0.2905, +0.0443),  # n=656
    (+1.4192, -0.1793, -0.2858),  # n=558
    (+1.4169, -0.0321, -0.4767),  # n=948
    (+1.4173, +0.0406, -0.5059),  # n=1478
    (+1.4172, +0.1045, -0.5091),  # n=2144
    (+1.4183, +0.1581, -0.4886),  # n=2614
    (+1.4189, +0.1972, -0.4665),  # n=2824
    (+1.4195, +0.2447, -0.4339),  # n=2277
    (+1.4209, +0.3013, -0.3824),  # n=1487
    (+1.4200, +0.3691, -0.2863),  # n=676
    (+1.4216, +0.3166, -0.0690),  # n=291
    (+1.4210, +0.0592, +0.3254),  # n=146
    (+1.4180, -0.1227, +0.4088),  # n=201
    (+1.4208, -0.2036, +0.4170),  # n=324
)
_N_BY_BIN: Final[tuple[tuple[float, float, float], ...]] = (
    (+2.4047, +0.1532, -0.2535),  # n=423
    (+2.4065, +0.1816, -0.2439),  # n=613
    (+2.4035, +0.2039, -0.2350),  # n=787
    (+2.4051, +0.2202, -0.1704),  # n=855
    (+2.4055, +0.2061, -0.0169),  # n=656
    (+2.4024, +0.1234, +0.1837),  # n=558
    (+2.3988, +0.0380, +0.2837),  # n=948
    (+2.3980, -0.0154, +0.3049),  # n=1478
    (+2.3965, -0.0550, +0.3108),  # n=2144
    (+2.3979, -0.0995, +0.2998),  # n=2614
    (+2.3979, -0.1298, +0.2879),  # n=2824
    (+2.3982, -0.1737, +0.2651),  # n=2277
    (+2.4003, -0.2111, +0.2318),  # n=1487
    (+2.4006, -0.2538, +0.1724),  # n=676
    (+2.4019, -0.2071, +0.0242),  # n=291
    (+2.4023, -0.0430, -0.2006),  # n=146
    (+2.4043, +0.0830, -0.2624),  # n=201
    (+2.4062, +0.1408, -0.2570),  # n=324
)

#: Used when the fourth alpha carbon is missing but the preceding one is present -- the last
#: peptide unit of a chain. This is the 3-CA answer: C 0.471 A, N 0.311 A.
_C_MARGINAL: Final[tuple[float, float, float]] = (+1.4193, +0.0752, -0.2886)
_N_MARGINAL: Final[tuple[float, float, float]] = (+2.3998, -0.0459, +0.1756)

#: Used for the FIRST peptide unit, whose frame must be built forward from CA(i), CA(i+1), CA(i+2)
#: because there is no CA(i-1). Measured in that convention over 19,402 units: C 0.343 A, N 0.225 A.
_C_FORWARD_MARGINAL: Final[tuple[float, float, float]] = (+1.4193, -0.0140, +0.4396)
_N_FORWARD_MARGINAL: Final[tuple[float, float, float]] = (+2.3997, +0.0130, -0.2806)

#: The per-bin tables as arrays, for the vectorized placement path. Identical values to the tuples
#: above; the NumPy form lets a whole chain's bins be gathered in one indexing operation.
_C_TABLE: Final[np.ndarray] = np.array(_C_BY_BIN, dtype=np.float64)
_N_TABLE: Final[np.ndarray] = np.array(_N_BY_BIN, dtype=np.float64)


@dataclass(frozen=True, slots=True)
class CABackboneResult:
    """Backbone atoms placed for a CA-only trace.

    Attributes
    ----------
    n_xyz, c_xyz, o_xyz
        ``(n_residues, 3)`` coordinates. Every row is finite: the two atoms no peptide unit covers
        (N of the first residue, C and O of the last) are placed by ideal internal geometry.
    source
        Per peptide unit, which predictor was used: ``"table"``, ``"marginal"`` or ``"forward"``.
        Length ``n_residues - 1``.
    notes
        What happened, for a caller that wants to report it.
    """

    n_xyz: np.ndarray
    c_xyz: np.ndarray
    o_xyz: np.ndarray
    source: tuple[str, ...]
    notes: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class SeamStrain:
    """A seam peptide bond that could not be satisfied, so was left long and reported.

    Where a rebuilt region meets a residue DODO did not rebuild, the peptide bond crosses from an
    atom this pass placed to one that was already there. When the fixed neighbour -- part of a
    folded domain moved rigidly in step 3 -- sits farther than the ``2.854 A`` a peptide unit can
    bridge, no placement satisfies both bonds, so the atom is aimed as close as the residue's own
    N-CA-C angle allows and the bond is left long. This records that it happened, on the residue it
    happened to, so a caller can surface it rather than leaving the validator to rediscover it.

    Attributes
    ----------
    residue
        0-based index of the rebuilt boundary residue whose placed atom was strained.
    side
        ``"C"`` -- this residue's carbon reaching the next residue's nitrogen -- or ``"N"`` -- this
        residue's nitrogen reaching the previous residue's carbon.
    bond_length
        The resulting peptide-bond length, in Angstroms; ideal is
        :data:`~dodo.constants.C_N_PEPTIDE_BOND_LENGTH` (1.329).
    """

    residue: int
    side: str
    bond_length: float


@dataclass(frozen=True, slots=True)
class BackboneResult:
    """What :func:`add_backbone_to_rebuilt` produced: the structure and what it could not close.

    Attributes
    ----------
    structure
        The rebuilt structure with N, C and O placed on every rebuilt region; folded domains are
        returned exactly as they arrived.
    strained_seams
        Every seam left with a long peptide bond, in residue order. Empty when every seam closed --
        always the case for a from-sequence build, which has no folded domains to seam against.
    """

    structure: Structure
    strained_seams: tuple[SeamStrain, ...] = ()


def _unit(v: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(v))
    if norm < 1e-9:
        raise GeometryError("Cannot normalize a zero-length vector while placing a backbone.")
    return v / norm


def _frame(along: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, bool]:
    """Right-handed frame: x along ``along``, z perpendicular to the along/reference plane.

    Returns the frame and whether it had to be invented. The second element matters: three
    collinear alpha carbons define no plane, and a caller that cannot tell the difference between
    a measured peptide plane and an arbitrary one will report the arbitrary one as though it were
    determined by the trace.
    """
    x = _unit(along)
    normal = np.cross(reference, along)
    degenerate = float(np.linalg.norm(normal)) < 1e-6
    if degenerate:
        # Any perpendicular will do; the peptide plane is genuinely undetermined here, and the
        # caller passes that on through `notes`.
        fallback = np.array([0.0, 0.0, 1.0])
        if abs(float(np.dot(fallback, x))) > 0.9:
            fallback = np.array([0.0, 1.0, 0.0])
        normal = np.cross(fallback, along)
    z = _unit(normal)
    return np.stack([x, np.cross(z, x), z]), degenerate


def _pseudo_dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    axis = _unit(b1)
    v = b0 - float(np.dot(b0, axis)) * axis
    w = b2 - float(np.dot(b2, axis)) * axis
    return float(np.degrees(np.arctan2(float(np.dot(np.cross(axis, v), w)), float(np.dot(v, w)))))


def _rotate(v: np.ndarray, axis: np.ndarray, radians: float) -> np.ndarray:
    axis = _unit(axis)
    cos, sin = float(np.cos(radians)), float(np.sin(radians))
    rotated: np.ndarray = (
        v * cos + np.cross(axis, v) * sin + axis * float(np.dot(axis, v)) * (1.0 - cos)
    )
    return rotated


def _on_two_spheres(
    centre_a: np.ndarray,
    radius_a: float,
    centre_b: np.ndarray,
    radius_b: float,
    toward: np.ndarray,
) -> np.ndarray | None:
    """Point at exact distances from two centres, on that circle, nearest ``toward``.

    Two spheres meet in a circle, so demanding both bonds leaves one free parameter -- which the
    table's prediction then resolves by picking the nearest point on it. Returns ``None`` when the
    spheres do not intersect, which cannot happen for a physically plausible CA-CA separation but
    is checked rather than assumed.
    """
    separation = float(np.linalg.norm(centre_b - centre_a))
    if (
        separation < 1e-9
        or separation > radius_a + radius_b
        or separation < abs(radius_a - radius_b)
    ):
        return None
    axis = (centre_b - centre_a) / separation
    along = (radius_a**2 - radius_b**2 + separation**2) / (2.0 * separation)
    height_squared = radius_a**2 - along**2
    if height_squared < 0.0:
        return None
    centre = centre_a + along * axis
    spoke = toward - centre
    spoke = spoke - float(np.dot(spoke, axis)) * axis
    if float(np.linalg.norm(spoke)) < 1e-9:
        arbitrary = np.array([1.0, 0.0, 0.0])
        spoke = arbitrary - float(np.dot(arbitrary, axis)) * axis
    resolved: np.ndarray = centre + np.sqrt(height_squared) * _unit(spoke)
    return resolved


def _place_carbonyl_oxygen(ca: np.ndarray, c: np.ndarray, n_next: np.ndarray) -> np.ndarray:
    """O is fully determined by CA(i), C(i) and N(i+1); from the true three it lands to 0.013 A."""
    to_ca = _unit(ca - c)
    to_n = _unit(n_next - c)
    normal = np.cross(to_ca, to_n)
    if float(np.linalg.norm(normal)) < 1e-6:
        normal = np.array([0.0, 0.0, 1.0])
    placed: np.ndarray = c + C_O_BOND_LENGTH * _rotate(to_ca, normal, -np.radians(CA_C_O_ANGLE))
    return placed


_Z_AXIS: Final[np.ndarray] = np.array([0.0, 0.0, 1.0])
_Y_AXIS: Final[np.ndarray] = np.array([0.0, 1.0, 0.0])
_X_AXIS: Final[np.ndarray] = np.array([1.0, 0.0, 0.0])


def _unit_rows(v: np.ndarray) -> np.ndarray:
    """Normalize along the last axis, flooring the divisor so no row can divide by zero.

    The scalar :func:`_unit` raises on a zero-length vector; here the only genuinely zero row would
    be a zero virtual bond, which :func:`_backbone_atoms` rejects up front, so the floor never bites
    on valid input and simply keeps the discarded branch of a :func:`numpy.where` finite.
    """
    norm = np.linalg.norm(v, axis=-1, keepdims=True)
    unit: np.ndarray = v / np.maximum(norm, 1e-12)
    return unit


def _rotate_rows(v: np.ndarray, axis: np.ndarray, radians: float) -> np.ndarray:
    """Rodrigues rotation of every row of ``v`` about the matching row of ``axis``."""
    axis = _unit_rows(axis)
    cos, sin = float(np.cos(radians)), float(np.sin(radians))
    dot = np.sum(axis * v, axis=-1, keepdims=True)
    rotated: np.ndarray = v * cos + np.cross(axis, v) * sin + axis * dot * (1.0 - cos)
    return rotated


def _frames_rows(
    along: np.ndarray, reference: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Batched :func:`_frame`: x along ``along``, z perpendicular to the along/reference plane.

    Returns the three basis vectors and a boolean mask of the units whose alpha carbons were
    collinear -- the same ``degenerate`` signal the scalar path reports.
    """
    x = _unit_rows(along)
    normal = np.cross(reference, along)
    degenerate = np.linalg.norm(normal, axis=-1) < 1e-6
    use_y = np.abs(x[..., 2]) > 0.9
    fallback = np.where(use_y[..., None], _Y_AXIS, _Z_AXIS)
    normal = np.where(degenerate[..., None], np.cross(fallback, along), normal)
    z = _unit_rows(normal)
    y = np.cross(z, x)
    return x, y, z, degenerate


def _dihedral_rows(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> np.ndarray:
    """Batched form of :func:`_pseudo_dihedral`, in degrees."""
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    axis = _unit_rows(b1)
    v = b0 - np.sum(b0 * axis, axis=-1, keepdims=True) * axis
    w = b2 - np.sum(b2 * axis, axis=-1, keepdims=True) * axis
    y = np.sum(np.cross(axis, v) * w, axis=-1)
    x = np.sum(v * w, axis=-1)
    degrees: np.ndarray = np.degrees(np.arctan2(y, x))
    return degrees


def _two_spheres_rows(
    centre_a: np.ndarray,
    radius_a: float,
    centre_b: np.ndarray,
    radius_b: float,
    toward: np.ndarray,
) -> np.ndarray:
    """Batched form of :func:`_on_two_spheres`.

    Where the spheres do not intersect -- the ``None`` case in the scalar version -- this returns a
    point one ``radius_a`` bond from ``centre_a`` toward ``toward``, exactly the fallback the scalar
    caller applies.
    """
    diff = centre_b - centre_a
    sep = np.linalg.norm(diff, axis=-1)
    sep_safe = np.maximum(sep, 1e-12)
    axis = diff / sep_safe[..., None]
    along = (radius_a**2 - radius_b**2 + sep**2) / (2.0 * sep_safe)
    height_squared = radius_a**2 - along**2
    valid = (
        (sep >= 1e-9)
        & (sep <= radius_a + radius_b)
        & (sep >= abs(radius_a - radius_b))
        & (height_squared >= 0.0)
    )
    centre = centre_a + along[..., None] * axis
    spoke = toward - centre
    spoke = spoke - np.sum(spoke * axis, axis=-1, keepdims=True) * axis
    arbitrary = _X_AXIS - np.sum(_X_AXIS * axis, axis=-1, keepdims=True) * axis
    spoke = np.where(np.linalg.norm(spoke, axis=-1, keepdims=True) < 1e-9, arbitrary, spoke)
    resolved = centre + np.sqrt(np.maximum(height_squared, 0.0))[..., None] * _unit_rows(spoke)
    fallback = centre_a + radius_a * _unit_rows(toward - centre_a)
    return np.where(valid[..., None], resolved, fallback)


def _carbonyl_rows(ca: np.ndarray, c: np.ndarray, n_next: np.ndarray) -> np.ndarray:
    """Batched form of :func:`_place_carbonyl_oxygen`."""
    to_ca = _unit_rows(ca - c)
    to_n = _unit_rows(n_next - c)
    normal = np.cross(to_ca, to_n)
    normal = np.where(np.linalg.norm(normal, axis=-1, keepdims=True) < 1e-6, _Z_AXIS, normal)
    return c + C_O_BOND_LENGTH * _rotate_rows(to_ca, normal, -np.radians(CA_C_O_ANGLE))


def _terminal_nitrogen_rows(ca: np.ndarray, c: np.ndarray, ca_next: np.ndarray) -> np.ndarray:
    """Batched form of :func:`_terminal_nitrogen`."""
    to_c = _unit_rows(c - ca)
    normal = np.cross(to_c, ca_next - ca)
    normal = np.where(np.linalg.norm(normal, axis=-1, keepdims=True) < 1e-6, _Z_AXIS, normal)
    return ca + N_CA_BOND_LENGTH * _rotate_rows(to_c, normal, np.radians(N_CA_C_ANGLE))


def _terminal_carbon_rows(ca: np.ndarray, n: np.ndarray, ca_prev: np.ndarray) -> np.ndarray:
    """Batched form of :func:`_terminal_carbon`."""
    to_n = _unit_rows(n - ca)
    normal = np.cross(to_n, ca_prev - ca)
    normal = np.where(np.linalg.norm(normal, axis=-1, keepdims=True) < 1e-6, _Z_AXIS, normal)
    return ca + CA_C_BOND_LENGTH * _rotate_rows(to_n, normal, -np.radians(N_CA_C_ANGLE))


def _terminal_oxygen_rows(ca: np.ndarray, c: np.ndarray, ca_prev: np.ndarray) -> np.ndarray:
    """Batched form of :func:`_terminal_oxygen`."""
    to_ca = _unit_rows(ca - c)
    normal = np.cross(to_ca, ca_prev - c)
    normal = np.where(np.linalg.norm(normal, axis=-1, keepdims=True) < 1e-6, _Z_AXIS, normal)
    return c + C_O_BOND_LENGTH * _rotate_rows(to_ca, normal, -np.radians(CA_C_O_ANGLE))


def _backbone_atoms(
    ca: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, tuple[str, ...]]:
    """Place N, C and O for every peptide unit of one or many CA traces, with no per-residue loop.

    ``ca`` is ``(..., n_residues, 3)`` -- a single trace as ``(n, 3)`` or a batch of conformers as
    ``(B, n, 3)``. Every step runs over the residue axis (and any leading batch axes) at once: a
    peptide unit's placement depends only on four alpha carbons that are all known up front, so
    nothing here is sequential and there is nothing to loop over in Python.

    Returns the three atom arrays (same leading shape as ``ca``), the per-unit ``degenerate`` mask,
    and the per-unit predictor names. A unit's predictor is fixed by its position in the chain, so
    ``source`` is one tuple shared by every conformer in a batch.
    """
    n = ca.shape[-2]
    units = n - 1
    along = ca[..., 1:, :] - ca[..., :-1, :]
    if np.any(np.linalg.norm(along, axis=-1) < 1e-9):
        raise GeometryError(
            "Two consecutive alpha carbons coincide, giving a zero-length virtual bond."
        )

    # Reference axis fixing each frame's peptide plane. Past the first unit it is the incoming CA-CA
    # bond; the first unit has no predecessor, so it looks forward (or, for a two-residue chain with
    # nothing to look at, to a fixed axis).
    reference = np.empty_like(along)
    if units >= 2:
        reference[..., 1:, :] = ca[..., 1:units, :] - ca[..., : units - 1, :]
    # ca[1] - ca[2], NOT ca[2] - ca[1]. The forward marginal (_C/_N_FORWARD_MARGINAL) was MEASURED
    # in the frame this sign builds -- verified by re-deriving it in scripts/derive_peptide_table.py
    # -- and the frame's z is cross(reference, along), so the opposite sign flips z (and y) and
    # mis-places the first peptide unit of every region: measured C error 1.06 A against the 0.31 A
    # the marginal delivers in its own frame. Do not "simplify" this to ca[2] - ca[1].
    reference[..., 0, :] = ca[..., 1, :] - ca[..., 2, :] if n >= 3 else _Z_AXIS

    # Predictor per unit, fixed by position: interior units have all four CAs and use the table; the
    # last uses the trailing marginal; the first uses the leading one.
    index = np.arange(units)
    is_table = (index >= 1) & (index <= n - 3)
    is_marginal = (index == n - 2) & (n >= 3)

    local_c = np.broadcast_to(np.array(_C_FORWARD_MARGINAL), along.shape).copy()
    local_n = np.broadcast_to(np.array(_N_FORWARD_MARGINAL), along.shape).copy()
    if n >= 3:
        local_c[..., n - 2, :] = _C_MARGINAL
        local_n[..., n - 2, :] = _N_MARGINAL
    table_units = np.nonzero(is_table)[0]
    if table_units.size:
        tau = _dihedral_rows(
            ca[..., table_units - 1, :],
            ca[..., table_units, :],
            ca[..., table_units + 1, :],
            ca[..., table_units + 2, :],
        )
        bins = np.minimum(((tau + 180.0) // _BIN_WIDTH_DEGREES).astype(np.int64), len(_C_TABLE) - 1)
        local_c[..., table_units, :] = _C_TABLE[bins]
        local_n[..., table_units, :] = _N_TABLE[bins]

    x, y, z, degenerate = _frames_rows(along, reference)
    apex = ca[..., :units, :]
    placed_c = apex + local_c[..., 0:1] * x + local_c[..., 1:2] * y + local_c[..., 2:3] * z
    placed_n = apex + local_n[..., 0:1] * x + local_n[..., 1:2] * y + local_n[..., 2:3] * z

    n_xyz = np.full(ca.shape, np.nan)
    c_xyz = np.full(ca.shape, np.nan)
    o_xyz = np.full(ca.shape, np.nan)

    # N belongs to residue i+1; anchor it to that residue's own alpha carbon so N-CA is exact.
    n_owner = ca[..., 1 : units + 1, :]
    n_xyz[..., 1 : units + 1, :] = n_owner + N_CA_BOND_LENGTH * _unit_rows(placed_n - n_owner)
    # C satisfies both its bonds; the table's prediction chooses where on the intersection circle.
    c_xyz[..., :units, :] = _two_spheres_rows(
        apex, CA_C_BOND_LENGTH, n_xyz[..., 1 : units + 1, :], C_N_PEPTIDE_BOND_LENGTH, placed_c
    )
    o_xyz[..., :units, :] = _carbonyl_rows(
        apex, c_xyz[..., :units, :], n_xyz[..., 1 : units + 1, :]
    )

    # The two atoms no peptide unit covers: N of the first residue, C and O of the last.
    n_xyz[..., 0, :] = _terminal_nitrogen_rows(ca[..., 0, :], c_xyz[..., 0, :], ca[..., 1, :])
    c_xyz[..., -1, :] = _terminal_carbon_rows(ca[..., -1, :], n_xyz[..., -1, :], ca[..., -2, :])
    o_xyz[..., -1, :] = _terminal_oxygen_rows(ca[..., -1, :], c_xyz[..., -1, :], ca[..., -2, :])

    sources = tuple(
        "table" if is_table[i] else "marginal" if is_marginal[i] else "forward"
        for i in range(units)
    )
    return n_xyz, c_xyz, o_xyz, degenerate, sources


def backbone_from_ca(ca_coords: np.ndarray) -> CABackboneResult:
    """Place N, C and O for every residue of a CA-only trace.

    Parameters
    ----------
    ca_coords
        ``(n_residues, 3)`` alpha-carbon coordinates of one continuous chain segment, in order.
        At least two residues.

    Returns
    -------
    CABackboneResult
        Backbone coordinates, all finite, with ideal bond lengths by construction.

    Raises
    ------
    ~dodo.exceptions.GeometryError
        If the input is not ``(n, 3)`` with ``n >= 2``, or contains a non-finite coordinate.
    """
    ca = np.asarray(ca_coords, dtype=np.float64)
    if ca.ndim != 2 or ca.shape[1] != 3:
        raise GeometryError(f"Expected CA coordinates shaped (n_residues, 3), got {ca.shape}.")
    if ca.shape[0] < 2:
        raise GeometryError(
            f"A backbone needs at least two alpha carbons to define one peptide unit, got "
            f"{ca.shape[0]}."
        )
    if not np.all(np.isfinite(ca)):
        raise GeometryError("CA coordinates contain a non-finite value.")

    n_res = ca.shape[0]
    n_xyz, c_xyz, o_xyz, degenerate, sources = _backbone_atoms(ca)
    collinear = int(np.count_nonzero(degenerate))

    notes = [
        f"Placed {n_res} residue(s) of backbone from alpha carbons: "
        f"{sources.count('table')} peptide unit(s) from the four-CA table, "
        f"{sources.count('marginal')} from the trailing marginal, "
        f"{sources.count('forward')} from the leading marginal."
    ]
    if collinear:
        notes.append(f"{collinear} peptide unit(s) had collinear alpha carbons.")
    return CABackboneResult(
        n_xyz=n_xyz, c_xyz=c_xyz, o_xyz=o_xyz, source=tuple(sources), notes=tuple(notes)
    )


def _terminal_nitrogen(ca: np.ndarray, c: np.ndarray, ca_next: np.ndarray) -> np.ndarray:
    """N of the first residue, which no peptide unit covers.

    Fixed by the N-CA-C angle off the already-placed C. The remaining rotation about the CA-C axis
    is this residue's phi, which nothing in a CA trace constrains at a chain start, so it is
    resolved deterministically into the CA(0)/C(0)/CA(1) plane.
    """
    to_c = _unit(c - ca)
    reference = ca_next - ca
    normal = np.cross(to_c, reference)
    if float(np.linalg.norm(normal)) < 1e-6:
        normal = np.array([0.0, 0.0, 1.0])
    return ca + N_CA_BOND_LENGTH * _rotate(to_c, normal, np.radians(N_CA_C_ANGLE))


def _terminal_carbon(ca: np.ndarray, n: np.ndarray, ca_prev: np.ndarray) -> np.ndarray:
    """C of the last residue, the mirror case of :func:`_terminal_nitrogen`."""
    to_n = _unit(n - ca)
    reference = ca_prev - ca
    normal = np.cross(to_n, reference)
    if float(np.linalg.norm(normal)) < 1e-6:
        normal = np.array([0.0, 0.0, 1.0])
    return ca + CA_C_BOND_LENGTH * _rotate(to_n, normal, -np.radians(N_CA_C_ANGLE))


def _terminal_oxygen(ca: np.ndarray, c: np.ndarray, ca_prev: np.ndarray) -> np.ndarray:
    """O of the last residue, which has no following N to orient the carbonyl against.

    The obvious shortcut is to invent that nitrogen by continuing the chain straight out along
    CA->C, and this function exists because that is wrong in a way that is easy to miss: a notional
    N on the CA-C line is *collinear* with it, so the cross product defining the carbonyl plane
    degenerates to zero and the normal falls back to an arbitrary axis. On a 22-residue test
    sequence that put the final O 0.981 A from its own CA -- two unbonded atoms on top of each
    other, in a chain with no seams to blame.

    So the plane is taken from the previous alpha carbon instead, which is always available (this is
    only reached with at least two residues) and never collinear with CA-C in practice. Which side
    of that plane O falls on is this residue's psi, and nothing in a CA trace constrains psi at a
    free C-terminus, so the choice is arbitrary but consistent.
    """
    to_ca = _unit(ca - c)
    normal = np.cross(to_ca, ca_prev - c)
    if float(np.linalg.norm(normal)) < 1e-6:
        normal = np.array([0.0, 0.0, 1.0])
    placed: np.ndarray = c + C_O_BOND_LENGTH * _rotate(to_ca, normal, -np.radians(CA_C_O_ANGLE))
    return placed


def _existing_atom(structure: Structure, residue: int, name: str) -> np.ndarray | None:
    """Coordinates of a named atom already present in a residue, or None if it has none."""
    atoms = structure.atom_slice_for_residues(residue, residue + 1)
    for index in range(atoms.start, atoms.stop):
        if str(structure.atom_name[index]) == name:
            found: np.ndarray = structure.xyz[index]
            return found
    return None


#: Separation a seam placement must keep from the neighbouring residue's atoms, in Angstroms.
#:
#: MEASURED as the smallest value that stops the recurring seam defect. Real non-bonded backbone
#: contacts across a peptide junction run from about 2.7 A (a carbonyl O to the next residue's CB)
#: upward, so a 2.6 A floor rejects a collision without rejecting legitimate geometry. It is
#: checked rather than assumed -- see :func:`_place_on_cone` for the three successive
#: by-construction fixes this replaced.
SEAM_CLASH_DISTANCE: Final[float] = 2.6


def _place_on_cone(
    apex: np.ndarray,
    reference: np.ndarray,
    angle_degrees: float,
    length: float,
    *,
    prefer: np.ndarray | None = None,
    depart: np.ndarray | None = None,
    avoid: np.ndarray | None = None,
    clash: float = SEAM_CLASH_DISTANCE,
    samples: int = 48,
) -> np.ndarray:
    """Place an atom at ``length`` from ``apex`` and ``angle_degrees`` off ``apex``->``reference``.

    That leaves exactly one free parameter -- the azimuth about the reference axis -- and this picks
    it. Collisions are ruled out first, then ``prefer`` is approached and ``depart`` is escaped.

    This exists because the strained-seam fallback kept producing the same class of defect in a
    new place each time it was fixed by construction alone. Aiming the carbon straight at the
    neighbour's nitrogen made the three collinear and sent the carbonyl to an arbitrary axis
    (measured: O 0.6 A from its own alpha carbon on p300). Holding the N-CA-C angle instead fixed
    that, and let the oxygen land 0.975 A from the neighbour's nitrogen. Placing the oxygen trans
    to that nitrogen fixed *that*, and left it 0.96 A from the neighbour's CB. Every fix was
    locally right and blind to one more neighbour; the way to stop is to check instead of
    construct, which is what this does.

    The angle to ``reference`` is never traded away: it is the residue's own internal geometry, and
    damaging that is worse than a long bond to the next residue.
    """
    axis = _unit(reference - apex)
    fallback = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(fallback, axis))) > 0.9:
        fallback = np.array([0.0, 1.0, 0.0])
    first = _unit(fallback - float(np.dot(fallback, axis)) * axis)
    second = np.cross(axis, first)

    angle = np.radians(angle_degrees)
    azimuths = np.linspace(0.0, 2.0 * np.pi, samples, endpoint=False)
    radial = np.cos(azimuths)[:, None] * first + np.sin(azimuths)[:, None] * second
    candidates = apex + length * (np.cos(angle) * axis + np.sin(angle) * radial)

    score = np.zeros(samples)
    if avoid is not None and avoid.shape[0]:
        gaps = np.linalg.norm(candidates[:, None, :] - avoid[None, :, :], axis=2)
        # A hard penalty, not a weighted term: a candidate that collides is not a candidate.
        score += 1e6 * np.square(np.clip(clash - gaps, 0.0, None)).sum(axis=1)
    if prefer is not None:
        score += np.linalg.norm(candidates - prefer, axis=1)
    if depart is not None:
        score -= np.linalg.norm(candidates - depart, axis=1)
    chosen: np.ndarray = candidates[int(np.argmin(score))]
    return chosen


def _seam_obstacles(structure: Structure, residues: Iterable[int]) -> np.ndarray:
    """Atoms of ``residues`` that a seam placement must not collide with.

    Deliberately narrow: the atoms that actually get hit are the neighbouring residue's own, and its
    side chain in particular -- an anchor's CB is the thing a carbonyl oxygen runs into.
    """
    collected: list[np.ndarray] = []
    for residue in residues:
        if not 0 <= residue < structure.n_residues:
            continue
        atoms = structure.atom_slice_for_residues(residue, residue + 1)
        for index in range(atoms.start, atoms.stop):
            collected.append(structure.xyz[index])
    if not collected:
        return np.empty((0, 3))
    return np.asarray(collected, dtype=np.float64)


def _polish_coupled_clashes(
    placed: dict[int, dict[str, np.ndarray]],
    structure: Structure,
    generated: list[tuple[int, int]],
    *,
    rounds: int = 8,
    grid: int = 25,
    angle_window: tuple[float, float] = (80.0, 160.0),
) -> int:
    """Clear steric clashes coupled between two movable peptide units. Modifies ``placed``.

    The per-unit refinement (:func:`~dodo.construct.backbone_refine.refine_backbone`) moves one
    azimuth at a time, so a clash coupled between two movable units -- where clearing it by rotating
    either one alone makes another contact worse -- sits at a local minimum that more sweeps cannot
    escape (measured: 200 sweeps leaves the same clashes as 30, a wide window barely helps). For
    each residual clash this jointly searches the azimuths of the one or two *interior* units (both
    flanking alpha carbons rebuilt, so only rebuilt atoms move -- never a folded-domain or seam
    atom) that place the clashing atoms, keeping any combination that lowers the total van der Waals
    overlap without pushing a moved residue's N-CA-C angle out of ``angle_window``.

    The overlap and its element limits are the validator's own
    (:func:`~dodo.validate.clashes.contact_limit` constants), and only pairs more than one residue
    apart are scored -- consecutive-residue contacts are backbone bonds or their 1-3/1-4 neighbours,
    and seam contacts belong to a different problem. Returns the number of joint moves applied.
    """
    from scipy.spatial import cKDTree

    from ..validate.clashes import (
        ALLOWED_OVERLAP,
        DEFAULT_VDW_RADIUS,
        HBOND_ALLOWED_OVERLAP,
        POLAR_ELEMENTS,
        VDW_RADII,
    )
    from .backbone_refine import _angle, _azimuth_frame, _place_oxygen, _place_unit

    rebuilt: set[int] = set()
    for start, stop in generated:
        rebuilt.update(range(start, stop))
    ca = structure.ca_xyz
    n_res = structure.n_residues
    units = [
        u
        for u in range(n_res - 1)
        if u in rebuilt and (u + 1) in rebuilt and u in placed and (u + 1) in placed
    ]
    if not units:
        return 0
    unit_set = set(units)

    azimuth: dict[int, float] = {}
    for u in units:
        axis, first, second = _azimuth_frame(ca[u], ca[u + 1])
        offset = placed[u]["C"] - ca[u]
        radial = offset - float(np.dot(offset, axis)) * axis
        azimuth[u] = float(
            np.degrees(np.arctan2(float(np.dot(radial, second)), float(np.dot(radial, first))))
        )

    def limit(ea: str, eb: str) -> float:
        ra = VDW_RADII.get(ea, DEFAULT_VDW_RADIUS)
        rb = VDW_RADII.get(eb, DEFAULT_VDW_RADIUS)
        both_polar = ea in POLAR_ELEMENTS and eb in POLAR_ELEMENTS
        return ra + rb - (HBOND_ALLOWED_OVERLAP if both_polar else ALLOWED_OVERLAP)

    query_r = 2.0 * max(VDW_RADII.get(e, DEFAULT_VDW_RADIUS) for e in ("C", "N", "O", "S"))

    # The whole atom cloud: every fixed atom (folded residues all-atom + rebuilt alpha carbons --
    # exactly what structure.xyz holds while rebuilt regions are CA-only) plus the placed N/C/O.
    fixed_xyz = np.asarray(structure.xyz, dtype=np.float64)
    fixed_res = np.asarray(structure.residue_index)
    fixed_elem = np.array([str(e).upper() for e in structure.element], dtype=object)
    mov_keys = [
        (r, nm)
        for r in sorted(rebuilt)
        if r in placed
        for nm in ("N", "C", "O")
        if nm in placed[r]  # a seam atom may have been omitted as unsatisfiable
    ]
    row = {key: len(fixed_xyz) + i for i, key in enumerate(mov_keys)}
    all_xyz = np.vstack([fixed_xyz, np.array([placed[r][nm] for r, nm in mov_keys])])
    all_res = np.concatenate([fixed_res, np.array([r for r, _ in mov_keys], dtype=fixed_res.dtype)])
    all_elem = np.concatenate([fixed_elem, np.array([nm for _, nm in mov_keys], dtype=object)])

    def unit_of(index: int) -> int | None:
        if index < len(fixed_xyz):
            return None  # a fixed atom -- folded, or a rebuilt CA -- never moves
        r, nm = mov_keys[index - len(fixed_xyz)]
        u = r - 1 if nm == "N" else r
        return u if u in unit_set else None

    def place_into(u: int, azimuth_deg: float, target: np.ndarray) -> None:
        c, n_next = _place_unit(ca[u], ca[u + 1], azimuth_deg)
        target[row[(u, "C")]] = c
        target[row[(u + 1, "N")]] = n_next
        target[row[(u, "O")]] = _place_oxygen(ca[u], c, n_next)

    lo, hi = angle_window

    def angle_ok(u: int, positions: np.ndarray) -> bool:
        # Moving unit u sets residue u's C and residue u+1's N; check the N-CA-C at both, using the
        # N and C the move does not touch (placed by the neighbouring units, held fixed here). A
        # residue whose N or C was omitted at a seam has no angle to check, so it is skipped.
        for res in (u, u + 1):
            if (res, "N") in row and (res, "C") in row:
                angle = _angle(positions[row[(res, "N")]], ca[res], positions[row[(res, "C")]])
                if not lo <= angle <= hi:
                    return False
        return True

    def overlap(
        moved_rows: list[int],
        positions: np.ndarray,
        tree: Any,
        txyz: np.ndarray,
        tres: np.ndarray,
        telem: np.ndarray,
    ) -> float:
        total = 0.0
        for r_idx in moved_rows:
            point = positions[r_idx]
            residue = int(all_res[r_idx])
            element = str(all_elem[r_idx])
            for j in tree.query_ball_point(point, query_r):
                if abs(int(tres[j]) - residue) <= 1:
                    continue
                distance = float(np.linalg.norm(point - txyz[j]))
                headroom = limit(element, str(telem[j])) - distance
                if headroom > 0.0:
                    total += headroom
        # Overlap AMONG the moved atoms themselves -- excluded from the static tree, but it is
        # exactly the clash a two-unit joint search is trying to clear (an atom of one moved unit
        # against an atom of the other), so it has to be counted or the search is blind to it.
        for a_pos in range(len(moved_rows)):
            for b_pos in range(a_pos + 1, len(moved_rows)):
                ia, ib = moved_rows[a_pos], moved_rows[b_pos]
                if abs(int(all_res[ia]) - int(all_res[ib])) <= 1:
                    continue
                distance = float(np.linalg.norm(positions[ia] - positions[ib]))
                headroom = limit(str(all_elem[ia]), str(all_elem[ib])) - distance
                if headroom > 0.0:
                    total += headroom
        return total

    def try_polish(unit_group: tuple[int, ...]) -> bool:
        moved = [r for u in unit_group for r in (row[(u, "C")], row[(u + 1, "N")], row[(u, "O")])]
        static = np.ones(len(all_xyz), dtype=bool)
        static[moved] = False
        txyz, tres, telem = all_xyz[static], all_res[static], all_elem[static]
        tree = cKDTree(txyz)
        base = overlap(moved, all_xyz, tree, txyz, tres, telem)
        if base <= 1e-9:
            return False
        trial = all_xyz.copy()
        grid_offsets = np.linspace(-180.0, 180.0, grid)
        best = base
        best_choice: tuple[float, ...] | None = None

        def descend(depth: int, choice: tuple[float, ...]) -> None:
            nonlocal best, best_choice
            u = unit_group[depth]
            for delta in grid_offsets:
                z = azimuth[u] + float(delta)
                place_into(u, z, trial)
                if not angle_ok(u, trial):
                    continue
                if depth + 1 < len(unit_group):
                    descend(depth + 1, (*choice, z))
                else:
                    value = overlap(moved, trial, tree, txyz, tres, telem)
                    if value < best - 1e-9:
                        best = value
                        best_choice = (*choice, z)

        descend(0, ())
        if best_choice is None:
            return False
        for u, z in zip(unit_group, best_choice, strict=True):
            place_into(u, z, all_xyz)
            azimuth[u] = ((z + 180.0) % 360.0) - 180.0
        return True

    moves = 0
    for _ in range(rounds):
        tree = cKDTree(all_xyz)
        groups: set[tuple[int, ...]] = set()
        for i, j in tree.query_pairs(query_r, output_type="ndarray"):
            if abs(int(all_res[i]) - int(all_res[j])) <= 1:
                continue
            distance = float(np.linalg.norm(all_xyz[i] - all_xyz[j]))
            if distance >= limit(str(all_elem[i]), str(all_elem[j])):
                continue
            involved = tuple(sorted({u for u in (unit_of(i), unit_of(j)) if u is not None}))
            if involved:
                groups.add(involved)
        if not groups:
            break
        applied = False
        for group in groups:
            if try_polish(group):
                applied = True
                moves += 1
        if not applied:
            break

    for u in units:
        placed[u]["C"] = all_xyz[row[(u, "C")]].copy()
        placed[u + 1]["N"] = all_xyz[row[(u + 1, "N")]].copy()
        placed[u]["O"] = all_xyz[row[(u, "O")]].copy()
    return moves


def add_backbone_to_rebuilt(
    structure: Structure,
    *,
    refine: bool = True,
    on_region_done: Callable[[int], None] | None = None,
) -> BackboneResult:
    """Return a copy with N, C and O placed on the regions DODO rebuilt.

    Only the rebuilt regions gain atoms. Folded domains are returned exactly as they arrived,
    with every atom they already had -- DODO moves them rigidly and never regenerates them, and
    that is not negotiable here either.

    Parameters
    ----------
    structure
        A rebuilt structure. Regions DODO generated are identified through
        :meth:`~dodo.structure.Domain.generated_spans`, so a region whose build failed is left
        alone rather than given a backbone on top of untouched input coordinates.
    refine
        Run :func:`~dodo.construct.backbone_refine.refine_backbone` afterwards. On by default:
        measured over 3,643 held-out residues it improves every criterion at once -- N 0.188 to
        0.164 A, C 0.283 to 0.210, O 0.816 to 0.614, and the N-CA-C spread 11.70 to 3.47 degrees
        against a real 2.94.

    on_region_done
        Called with each region's residue count as its backbone is finished. This pass is slow
        enough on a large structure to look like a hang without it -- on p300 it is most of the
        wall time -- so a caller needs to report progress through it, not merely before it.

    Returns
    -------
    BackboneResult
        The new structure (the input is not modified) plus every seam left with a strained
        peptide bond, so the caller can report what could not be closed rather than leaving the
        validator to rediscover it.
    """
    from ..structure import Structure

    generated: list[tuple[int, int]] = []
    for domain in structure.domains:
        for span in domain.generated_spans():
            generated.append((span.start, span.stop))
    if not generated:
        return BackboneResult(structure=structure)

    # Everything DODO did not generate: folded domains, and any region whose build failed and so
    # still holds its input coordinates. Also every alpha carbon, including the generated ones,
    # since a new backbone atom must not be placed on top of a neighbouring region's CA.
    inherited = structure.xyz[
        ~np.isin(
            structure.residue_index,
            np.concatenate([np.arange(a, b) for a, b in generated]),
        )
    ]
    obstacle_blocks = [inherited, structure.ca_xyz]

    placed: dict[int, dict[str, np.ndarray]] = {}
    strained: list[SeamStrain] = []
    for start, stop in generated:
        ca = structure.ca_xyz[start:stop]
        if ca.shape[0] < 2 or not np.all(np.isfinite(ca)):
            continue
        built = backbone_from_ca(ca)
        n_xyz, c_xyz, o_xyz = built.n_xyz, built.c_xyz, built.o_xyz
        if refine and ca.shape[0] >= 3:
            from .backbone_refine import refine_backbone

            # Regions are refined one at a time, so each one is given the backbone already placed
            # for the earlier ones. Without this they are mutually blind: refinement avoids the
            # folded domains and, within a region, itself, but two separate regions passing close
            # together would each place atoms into the other's space. Measured on p300, which has
            # six generated regions, that accounted for 4 of 19 remaining steric clashes.
            refined = refine_backbone(ca, n_xyz, c_xyz, obstacles=np.vstack(obstacle_blocks))
            n_xyz, c_xyz, o_xyz = refined.n_xyz, refined.c_xyz, refined.o_xyz
        obstacle_blocks.append(np.vstack([n_xyz, c_xyz, o_xyz]))
        if on_region_done is not None:
            on_region_done(stop - start)
        for offset, residue in enumerate(range(start, stop)):
            placed[residue] = {"N": n_xyz[offset], "C": c_xyz[offset], "O": o_xyz[offset]}

        # SEAMS. A rebuilt region is placed in isolation, so where it meets a residue DODO did not
        # rebuild, the peptide bond across the join is between an atom this function placed and one
        # that was already there -- and nothing has reconciled them. Measured on dnmt3a before this
        # was handled: C(A:PRO282)-N(A:GLU283) came out at 1.222 A against an ideal 1.329, and the
        # loop boundary at A:LYS382/A:LEU383 at 2.325 A.
        #
        # Both ends are fixed by construction against the REAL neighbour rather than a predicted
        # one: the outgoing C onto the next residue's existing N, and the incoming N onto the
        # previous residue's existing C.
        if stop < structure.n_residues:
            neighbour_n = _existing_atom(structure, stop, "N")
            if neighbour_n is not None:
                anchor_ca = structure.ca_xyz[stop - 1]
                seam = _on_two_spheres(
                    anchor_ca,
                    CA_C_BOND_LENGTH,
                    neighbour_n,
                    C_N_PEPTIDE_BOND_LENGTH,
                    placed[stop - 1]["C"],
                )
                if seam is None:
                    # Unreachable, and for a real reason rather than a numerical one. The existing
                    # N belongs to a residue DODO did not rebuild, so it still points toward where
                    # this region ran in the ORIGINAL model -- the domain was moved rigidly and the
                    # region re-drawn beneath it. Measured on dnmt3a's loop, CA(392) sits 3.741 A
                    # from N(393) where a peptide unit needs about 2.43, so no C satisfies both
                    # bonds. Aim C straight at the neighbour instead: the seam bond stays too long
                    # and is reported, but nothing impossible is written.
                    seam = _place_on_cone(
                        anchor_ca,
                        placed[stop - 1]["N"],
                        N_CA_C_ANGLE,
                        CA_C_BOND_LENGTH,
                        prefer=neighbour_n,
                        avoid=_seam_obstacles(structure, (stop, stop + 1)),
                    )
                    strained.append(
                        SeamStrain(
                            residue=stop - 1,
                            side="C",
                            bond_length=float(np.linalg.norm(seam - neighbour_n)),
                        )
                    )
                placed[stop - 1]["C"] = seam
                # The carbonyl needs the CA-C-N plane, which the two-spheres C (or the cone
                # fallback) makes well-defined -- it sits off the CA->N axis -- so O lands *trans*
                # to the neighbour's nitrogen rather than on an arbitrary axis.
                placed[stop - 1]["O"] = _place_carbonyl_oxygen(anchor_ca, seam, neighbour_n)
        if start > 0:
            neighbour_c = _existing_atom(structure, start - 1, "C")
            if neighbour_c is not None:
                first_ca = structure.ca_xyz[start]
                seam = _on_two_spheres(
                    first_ca,
                    N_CA_BOND_LENGTH,
                    neighbour_c,
                    C_N_PEPTIDE_BOND_LENGTH,
                    placed[start]["N"],
                )
                if seam is None:
                    # Mirror of the C-side fallback: hold this residue's own N-CA-C angle and spend
                    # the remaining freedom reaching toward the neighbour's carbon.
                    seam = _place_on_cone(
                        first_ca,
                        placed[start]["C"],
                        N_CA_C_ANGLE,
                        N_CA_BOND_LENGTH,
                        prefer=neighbour_c,
                        avoid=_seam_obstacles(structure, (start - 1, start - 2)),
                    )
                    strained.append(
                        SeamStrain(
                            residue=start,
                            side="N",
                            bond_length=float(np.linalg.norm(seam - neighbour_c)),
                        )
                    )
                placed[start]["N"] = seam

    # With every region placed and the seams reconciled, clear the steric clashes the per-region
    # refinement leaves behind -- the ones coupled between two movable units, which single-azimuth
    # descent cannot escape. Done here, over the whole assembled structure, so a clash between two
    # different rebuilt regions is visible where a per-region pass could not see it.
    if refine:
        _polish_coupled_clashes(placed, structure, generated)

    # Rebuild the flat atom records, inserting N before and C, O after each rebuilt CA so the
    # output keeps the conventional N-CA-C-O order within a residue.
    xyz: list[np.ndarray] = []
    names: list[str] = []
    elements: list[str] = []
    residue_names: list[str] = []
    residue_numbers: list[int] = []
    chain_ids: list[str] = []
    insertion_codes: list[str] = []
    b_factors: list[float] = []
    occupancies: list[float] = []

    for residue in range(structure.n_residues):
        atoms = structure.atom_slice_for_residues(residue, residue + 1)
        chain = structure.chains[int(structure.chain_index[residue])]
        existing = {str(name) for name in structure.atom_name[atoms]}
        extra = placed.get(residue, {})

        # Collect this residue's atoms in conventional N-CA-C-O order, then append once. Built as
        # a list rather than through a closure so nothing captures the loop variable.
        # b_factor and occupancy are per-RESIDUE on a Structure, not per-atom, so a placed atom
        # inherits its residue's values rather than carrying its own.
        b_value = float(structure.b_factor[residue])
        occupancy = float(structure.occupancy[residue])

        records: list[tuple[str, str, np.ndarray]] = []
        if not extra:
            # Nothing is being added to this residue, so it is passed through byte-for-byte -- same
            # atoms, same coordinates, SAME ORDER. The order matters: reordering into the
            # conventional N-CA-C-O below would rewrite 864 of dnmt3a's 912 residues, including
            # every folded-domain residue, and "folded domains are returned exactly as they
            # arrived" has to mean exactly.
            for index in range(atoms.start, atoms.stop):
                records.append(
                    (
                        str(structure.atom_name[index]),
                        str(structure.element[index]),
                        structure.xyz[index],
                    )
                )
            for name, element, point in records:
                xyz.append(np.asarray(point, dtype=np.float64))
                names.append(name)
                elements.append(element)
                residue_names.append(str(structure.residue_name[residue]))
                residue_numbers.append(int(structure.residue_number[residue]))
                chain_ids.append(chain.chain_id)
                insertion_codes.append(str(structure.insertion_code[residue]))
                b_factors.append(b_value)
                occupancies.append(occupancy)
            continue

        # A residue that gains atoms is emitted in conventional N-CA-C-O order, then everything else
        # in the order it arrived. A placed atom is only ever offered for a residue that lacks it,
        # so an existing atom of the same name wins -- which is what keeps this additive.
        emitted: set[str] = set()
        for name in ("N", "CA", "C", "O"):
            if name in extra:
                records.append((name, name[0], extra[name]))
                emitted.add(name)
            elif name in existing:
                for index in range(atoms.start, atoms.stop):
                    if str(structure.atom_name[index]) == name:
                        records.append((name, str(structure.element[index]), structure.xyz[index]))
                        emitted.add(name)
                        break
        # Side chains, and anything else the residue carries, in the order it arrived.
        for index in range(atoms.start, atoms.stop):
            name = str(structure.atom_name[index])
            if name in emitted:
                continue
            records.append((name, str(structure.element[index]), structure.xyz[index]))

        for name, element, point in records:
            xyz.append(np.asarray(point, dtype=np.float64))
            names.append(name)
            elements.append(element)
            residue_names.append(str(structure.residue_name[residue]))
            residue_numbers.append(int(structure.residue_number[residue]))
            chain_ids.append(chain.chain_id)
            insertion_codes.append(str(structure.insertion_code[residue]))
            b_factors.append(b_value)
            occupancies.append(occupancy)

    rebuilt = Structure.from_atom_records(
        xyz=np.asarray(xyz, dtype=np.float64),
        atom_name=names,
        element=elements,
        residue_name=residue_names,
        residue_number=residue_numbers,
        chain_id=chain_ids,
        insertion_code=insertion_codes,
        b_factor=b_factors,
        occupancy=occupancies,
        source=f"{structure.source} + placed backbone",
    )

    # Carry the region assignment across. from_atom_records builds a bare structure, so without
    # this the result has no domains at all -- every downstream consumer that asks what was
    # rebuilt, including the validator's provenance and the region annotation in the writer,
    # would silently see nothing. Residues map one-to-one, so the spans transfer unchanged.
    if rebuilt.n_residues != structure.n_residues:
        raise GeometryError(
            f"Placing a backbone changed the residue count from {structure.n_residues} to "
            f"{rebuilt.n_residues}, which it must never do."
        )
    for source_chain, target_chain in zip(structure.chains, rebuilt.chains, strict=True):
        target_chain.domains = [
            replace(domain, structure=rebuilt, rebuilt_loops=set(domain.rebuilt_loops))
            for domain in source_chain.domains
        ]
    return BackboneResult(structure=rebuilt, strained_seams=tuple(strained))
