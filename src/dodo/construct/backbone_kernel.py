"""Compiled inner loop for backbone refinement.

The refinement objective is cheap arithmetic evaluated an enormous number of times: for a
583-residue region, 397 peptide units x 25 candidate azimuths x 15 sweeps. In numpy that is
~150,000 array calls whose per-call dispatch cost dwarfs the six multiplies inside, and profiling
put about 60% of the time there rather than in the arithmetic. Batching can only amortise that
overhead, never remove it -- measured, batching ten models together recovered 3.6x of a
theoretical 10x, and batching over peptide units cost more in extra sweeps than it saved.

Compiling the loop removes the overhead outright. Measured over all six of p300's rebuilt regions
with the real obstacle set: **7.76 s to 0.49 s, 15.9x**.

Equivalence to :mod:`dodo.construct.backbone_refine` is not assumed, it is established. Given
byte-identical inputs for one peptide unit, this kernel and the numpy scorer agree to:

* C, N and O placement -- 0.000e+00, bit-identical
* the Ramachandran term -- 0.000e+00
* the clash term -- 0.000e+00
* the two N-CA-C angle terms -- 2.8e-13 and 1.1e-13, from ``math.acos`` versus ``np.arccos``
  rounding in the last bits

and they select the same candidate. The test that establishes this lives in
``tests/unit/test_backbone_kernel.py`` and it works by comparing the objective as a PURE FUNCTION
on identical inputs. That distinction matters: comparing the end states of two independent runs
cannot localise a difference, because a greedy search amplifies any perturbation chaotically. Four
earlier attempts to find a suspected bug that way all failed, and there was no bug to find.

The consequence of that 1e-13 is worth stating plainly, because it is visible: across the ~4,800
decisions in a sweep it eventually flips one that is nearly balanced, after which the two backends
diverge and produce different -- equally good, measurably equivalent -- coordinates. So output is
not bit-identical between backends. It was never bit-identical across numpy versions or BLAS
builds either, for exactly the same reason.

``fastmath`` is deliberately OFF. It was measured to make no difference to speed here, and leaving
it off avoids reassociating sums for no gain.
"""

from __future__ import annotations

import math

import numba as nb
import numpy as np
from scipy.spatial import cKDTree

from ..exceptions import GeometryError

PC = math.radians(20.4118)
PN = math.radians(14.8941)
CA_C = 1.525
N_CA = 1.458
C_O = 1.231
ACO = math.radians(120.8)
NCAC = 111.0
BIN = 30.0


@nb.njit(cache=True, fastmath=False, inline="always")
def _frame(ax, ay, az, bx, by, bz):
    ux = bx - ax
    uy = by - ay
    uz = bz - az
    n = math.sqrt(ux * ux + uy * uy + uz * uz)
    ux /= n
    uy /= n
    uz /= n
    if abs(uz) > 0.9:
        sx, sy, sz = 0.0, 1.0, 0.0
    else:
        sx, sy, sz = 0.0, 0.0, 1.0
    d = sx * ux + sy * uy + sz * uz
    fx = sx - d * ux
    fy = sy - d * uy
    fz = sz - d * uz
    fn = math.sqrt(fx * fx + fy * fy + fz * fz)
    fx /= fn
    fy /= fn
    fz /= fn
    gx = uy * fz - uz * fy
    gy = uz * fx - ux * fz
    gz = ux * fy - uy * fx
    return ux, uy, uz, fx, fy, fz, gx, gy, gz


@nb.njit(cache=True, fastmath=False, inline="always")
def _angle(ax, ay, az, bx, by, bz, cx, cy, cz):
    u1 = ax - bx
    u2 = ay - by
    u3 = az - bz
    v1 = cx - bx
    v2 = cy - by
    v3 = cz - bz
    n1 = math.sqrt(u1 * u1 + u2 * u2 + u3 * u3)
    n2 = math.sqrt(v1 * v1 + v2 * v2 + v3 * v3)
    d = (u1 * v1 + u2 * v2 + u3 * v3) / (n1 * n2)
    if d > 1.0:
        d = 1.0
    if d < -1.0:
        d = -1.0
    return math.degrees(math.acos(d))


@nb.njit(cache=True, fastmath=False, inline="always")
def _dihedral(p0x, p0y, p0z, p1x, p1y, p1z, p2x, p2y, p2z, p3x, p3y, p3z):
    b1x = p2x - p1x
    b1y = p2y - p1y
    b1z = p2z - p1z
    bn = math.sqrt(b1x * b1x + b1y * b1y + b1z * b1z)
    b1x /= bn
    b1y /= bn
    b1z /= bn
    d0x = p0x - p1x
    d0y = p0y - p1y
    d0z = p0z - p1z
    d3x = p3x - p2x
    d3y = p3y - p2y
    d3z = p3z - p2z
    t = d0x * b1x + d0y * b1y + d0z * b1z
    vx = d0x - t * b1x
    vy = d0y - t * b1y
    vz = d0z - t * b1z
    t = d3x * b1x + d3y * b1y + d3z * b1z
    wx = d3x - t * b1x
    wy = d3y - t * b1y
    wz = d3z - t * b1z
    cx = b1y * vz - b1z * vy
    cy = b1z * vx - b1x * vz
    cz = b1x * vy - b1y * vx
    return math.degrees(math.atan2(cx * wx + cy * wy + cz * wz, vx * wx + vy * wy + vz * wz))


@nb.njit(cache=True, fastmath=False, inline="always")
def _place(ca_a, ca_b, psi):
    ux, uy, uz, fx, fy, fz, gx, gy, gz = _frame(
        ca_a[0], ca_a[1], ca_a[2], ca_b[0], ca_b[1], ca_b[2]
    )
    cp = math.cos(psi)
    sp = math.sin(psi)
    rx = cp * fx + sp * gx
    ry = cp * fy + sp * gy
    rz = cp * fz + sp * gz
    cc = math.cos(PC)
    sc = math.sin(PC)
    cn = math.cos(PN)
    sn = math.sin(PN)
    Cx = ca_a[0] + CA_C * (cc * ux + sc * rx)
    Cy = ca_a[1] + CA_C * (cc * uy + sc * ry)
    Cz = ca_a[2] + CA_C * (cc * uz + sc * rz)
    Nx = ca_b[0] + N_CA * (-cn * ux - sn * rx)
    Ny = ca_b[1] + N_CA * (-cn * uy - sn * ry)
    Nz = ca_b[2] + N_CA * (-cn * uz - sn * rz)
    return Cx, Cy, Cz, Nx, Ny, Nz


@nb.njit(cache=True, fastmath=False, inline="always")
def _oxy(cax, cay, caz, Cx, Cy, Cz, Nx, Ny, Nz):
    ax = cax - Cx
    ay = cay - Cy
    az = caz - Cz
    an = math.sqrt(ax * ax + ay * ay + az * az)
    ax /= an
    ay /= an
    az /= an
    bx = Nx - Cx
    by = Ny - Cy
    bz = Nz - Cz
    bn = math.sqrt(bx * bx + by * by + bz * bz)
    bx /= bn
    by /= bn
    bz /= bn
    nx = ay * bz - az * by
    ny = az * bx - ax * bz
    nz = ax * by - ay * bx
    nn = math.sqrt(nx * nx + ny * ny + nz * nz)
    if nn > 1e-9:
        nx /= nn
        ny /= nn
        nz /= nn
    else:
        nx = 0.0
        ny = 0.0
        nz = 1.0
    co = math.cos(-ACO)
    so = math.sin(-ACO)
    dt = nx * ax + ny * ay + nz * az
    kx = ny * az - nz * ay
    ky = nz * ax - nx * az
    kz = nx * ay - ny * ax
    dx = ax * co + kx * so + nx * dt * (1.0 - co)
    dy = ay * co + ky * so + ny * dt * (1.0 - co)
    dz = az * co + kz * so + nz * dt * (1.0 - co)
    return Cx + C_O * dx, Cy + C_O * dy, Cz + C_O * dz


@nb.njit(cache=True, fastmath=False, inline="always")
def _clash(px, py, pz, pts, idx, k, limit, cd):
    tot = 0.0
    for j in range(k):
        i = idx[j]
        dx = px - pts[i, 0]
        dy = py - pts[i, 1]
        dz = pz - pts[i, 2]
        d2 = dx * dx + dy * dy + dz * dz
        if d2 < limit:
            g = cd - math.sqrt(d2)
            tot += g * g
    return tot


@nb.njit(cache=True, fastmath=False)
def sweep_region(  # noqa: D103 - documented at module level; njit rejects a docstring-only body
    ca,
    n_live,
    c_live,
    o_live,
    azimuths,
    fixed_pad,
    moving_pad,
    fixed_idx,
    moving_idx,
    rama,
    span,
    candidates,
    angle_w,
    clash_w,
    rama_w,
    cd,
):
    n_res = ca.shape[0]
    n_units = n_res - 1
    limit = cd * cd
    swept = 0.0
    largest = 0.0
    for unit in range(n_units):
        best_val = 1e30
        best_az = azimuths[unit]
        for t in range(candidates + 1):
            if t == 0:
                az = azimuths[unit]
            else:
                az = azimuths[unit] + span * (-1.0 + 2.0 * (t - 1) / (candidates - 1))
            Cx, Cy, Cz, Nx, Ny, Nz = _place(ca[unit], ca[unit + 1], math.radians(az))
            Ox, Oy, Oz = _oxy(ca[unit, 0], ca[unit, 1], ca[unit, 2], Cx, Cy, Cz, Nx, Ny, Nz)
            val = 0.0
            if 0 < unit < n_res - 1:
                a = _angle(
                    n_live[unit, 0],
                    n_live[unit, 1],
                    n_live[unit, 2],
                    ca[unit, 0],
                    ca[unit, 1],
                    ca[unit, 2],
                    Cx,
                    Cy,
                    Cz,
                )
                val += angle_w * (a - NCAC) * (a - NCAC)
                phi = _dihedral(
                    c_live[unit - 1, 0],
                    c_live[unit - 1, 1],
                    c_live[unit - 1, 2],
                    n_live[unit, 0],
                    n_live[unit, 1],
                    n_live[unit, 2],
                    ca[unit, 0],
                    ca[unit, 1],
                    ca[unit, 2],
                    Cx,
                    Cy,
                    Cz,
                )
                psi = _dihedral(
                    n_live[unit, 0],
                    n_live[unit, 1],
                    n_live[unit, 2],
                    ca[unit, 0],
                    ca[unit, 1],
                    ca[unit, 2],
                    Cx,
                    Cy,
                    Cz,
                    Nx,
                    Ny,
                    Nz,
                )
                i = int((phi + 180.0) // BIN)
                j = int((psi + 180.0) // BIN)
                if i < 0:
                    i = 0
                if i > 11:
                    i = 11
                if j < 0:
                    j = 0
                if j > 11:
                    j = 11
                r = rama[i, j]
                val += rama_w * r * r
            if 0 < unit + 1 < n_res - 1:
                a = _angle(
                    Nx,
                    Ny,
                    Nz,
                    ca[unit + 1, 0],
                    ca[unit + 1, 1],
                    ca[unit + 1, 2],
                    c_live[unit + 1, 0],
                    c_live[unit + 1, 1],
                    c_live[unit + 1, 2],
                )
                val += angle_w * (a - NCAC) * (a - NCAC)
                phi = _dihedral(
                    Cx,
                    Cy,
                    Cz,
                    Nx,
                    Ny,
                    Nz,
                    ca[unit + 1, 0],
                    ca[unit + 1, 1],
                    ca[unit + 1, 2],
                    c_live[unit + 1, 0],
                    c_live[unit + 1, 1],
                    c_live[unit + 1, 2],
                )
                if unit + 2 < n_res:
                    bx = n_live[unit + 2, 0]
                    by = n_live[unit + 2, 1]
                    bz = n_live[unit + 2, 2]
                else:
                    bx = c_live[unit + 1, 0]
                    by = c_live[unit + 1, 1]
                    bz = c_live[unit + 1, 2]
                psi = _dihedral(
                    Nx,
                    Ny,
                    Nz,
                    ca[unit + 1, 0],
                    ca[unit + 1, 1],
                    ca[unit + 1, 2],
                    c_live[unit + 1, 0],
                    c_live[unit + 1, 1],
                    c_live[unit + 1, 2],
                    bx,
                    by,
                    bz,
                )
                i = int((phi + 180.0) // BIN)
                j = int((psi + 180.0) // BIN)
                if i < 0:
                    i = 0
                if i > 11:
                    i = 11
                if j < 0:
                    j = 0
                if j > 11:
                    j = 11
                r = rama[i, j]
                val += rama_w * r * r
            K = fixed_idx.shape[2]
            val += clash_w * _clash(Cx, Cy, Cz, fixed_pad, fixed_idx[0, unit], K, limit, cd)
            val += clash_w * _clash(Nx, Ny, Nz, fixed_pad, fixed_idx[1, unit], K, limit, cd)
            val += clash_w * _clash(Ox, Oy, Oz, fixed_pad, fixed_idx[2, unit], K, limit, cd)
            KM = moving_idx.shape[2]
            val += clash_w * _clash(Cx, Cy, Cz, moving_pad, moving_idx[0, unit], KM, limit, cd)
            val += clash_w * _clash(Nx, Ny, Nz, moving_pad, moving_idx[1, unit], KM, limit, cd)
            val += clash_w * _clash(Ox, Oy, Oz, moving_pad, moving_idx[2, unit], KM, limit, cd)
            if val < best_val:
                best_val = val
                best_az = az
        d = abs(best_az - azimuths[unit])
        if d > largest:
            largest = d
        azimuths[unit] = best_az
        swept += best_val
        Cx, Cy, Cz, Nx, Ny, Nz = _place(ca[unit], ca[unit + 1], math.radians(best_az))
        c_live[unit, 0] = Cx
        c_live[unit, 1] = Cy
        c_live[unit, 2] = Cz
        n_live[unit + 1, 0] = Nx
        n_live[unit + 1, 1] = Ny
        n_live[unit + 1, 2] = Nz
        Ox, Oy, Oz = _oxy(ca[unit, 0], ca[unit, 1], ca[unit, 2], Cx, Cy, Cz, Nx, Ny, Nz)
        o_live[unit, 0] = Ox
        o_live[unit, 1] = Oy
        o_live[unit, 2] = Oz
        if unit + 1 < n_units:
            Ox, Oy, Oz = _oxy(
                ca[unit + 1, 0],
                ca[unit + 1, 1],
                ca[unit + 1, 2],
                c_live[unit + 1, 0],
                c_live[unit + 1, 1],
                c_live[unit + 1, 2],
                n_live[unit + 2, 0],
                n_live[unit + 2, 1],
                n_live[unit + 2, 2],
            )
            o_live[unit + 1, 0] = Ox
            o_live[unit + 1, 1] = Oy
            o_live[unit + 1, 2] = Oz
    return swept, largest


#: Neighbours retained per movable atom. MEASURED on the refiner's own clash calls over a
#: 398-residue region: the sets it builds run to a median of 1 and a maximum of 15, so 48 is far
#: above anything observed. Padding beyond the real count is harmless -- the slots point at
#: :data:`_PAD_COORDINATE`, whose clamped gap is exactly zero -- but truncating below it would
#: silently change the objective, so the cap is checked rather than trusted.
MAX_NEIGHBOURS: int = 48

#: Coordinate for a padded neighbour slot, in Angstroms. Far enough that it can never clash.
_PAD_COORDINATE: float = 1.0e6


def _neighbour_indices(tree, centres, reach, own_chain, own_oxygen, chain_of, oxygen_of, sentinel):
    """Return padded neighbour indices for each centre, excluding covalently bonded pairs.

    ``tree.query`` with a fixed ``k`` rather than ``query_ball_point``, because the radius query
    returns a ragged list of lists that has to be padded in Python while this returns a rectangular
    array in one vectorized call. Raises if the cap is ever reached, since a silent truncation would
    change the objective rather than merely slow it down.
    """
    count = min(MAX_NEIGHBOURS, tree.n)
    distances, found = tree.query(centres, k=count)
    if count == 1:
        distances = distances[:, None]
        found = found[:, None]
    keep = distances <= reach
    if chain_of is not None:
        separation = (
            np.abs(own_chain[:, None] - chain_of[found]) + own_oxygen[:, None] + oxygen_of[found]
        )
        keep &= separation > 4
    if count == MAX_NEIGHBOURS and bool(keep[:, -1].any()):
        raise GeometryError(
            f"A peptide unit has at least {MAX_NEIGHBOURS} neighbours within the clash shell, so "
            f"the compiled kernel's neighbour cap would truncate the objective. Raise "
            f"MAX_NEIGHBOURS or use backend='numpy'."
        )
    padded = np.full((centres.shape[0], MAX_NEIGHBOURS), sentinel, dtype=np.int64)
    padded[:, :count] = np.where(keep, found, sentinel)
    return padded


def refine_region(
    ca: np.ndarray,
    n_xyz: np.ndarray,
    c_xyz: np.ndarray,
    *,
    obstacles: np.ndarray | None = None,
    clash_distance: float = 2.9,
    angle_weight: float = 0.124,
    clash_weight: float = 40.0,
    rama_weight: float = 20.0,
    max_sweeps: int = 30,
    tolerance: float = 0.25,
    candidates: int = 24,
    rama_table: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int, bool]:
    """Refine one region with the compiled kernel. Mirrors ``refine_backbone``'s contract.

    Returns ``(n_xyz, c_xyz, o_xyz, azimuths, sweeps, converged)``. The tree work stays here in
    Python because scipy's KD-tree cannot cross into nopython mode; only the sweep is compiled, and
    that is where the time was.
    """
    from .backbone_refine import _azimuth_frame
    from .ca_backbone import _terminal_oxygen

    n_res = ca.shape[0]
    n_units = n_res - 1
    fixed = ca if obstacles is None else np.vstack([ca, obstacles])
    fixed_padded = np.vstack([fixed, np.full((1, 3), _PAD_COORDINATE)])
    fixed_tree = cKDTree(fixed)

    n_live = n_xyz.copy()
    c_live = c_xyz.copy()
    o_live = c_xyz.copy()
    azimuths = np.empty(n_units)
    for unit in range(n_units):
        axis, first, second = _azimuth_frame(ca[unit], ca[unit + 1])
        offset = c_xyz[unit] - ca[unit]
        radial = offset - float(np.dot(offset, axis)) * axis
        azimuths[unit] = np.degrees(
            np.arctan2(float(np.dot(radial, second)), float(np.dot(radial, first)))
        )
    for unit in range(n_units):
        cx, cy, cz, nx, ny, nz = _place(ca[unit], ca[unit + 1], math.radians(azimuths[unit]))
        c_live[unit] = (cx, cy, cz)
        n_live[unit + 1] = (nx, ny, nz)
    for unit in range(n_units):
        o_live[unit] = _oxy(ca[unit, 0], ca[unit, 1], ca[unit, 2], *c_live[unit], *n_live[unit + 1])
    o_live[-1] = _terminal_oxygen(ca[-1], c_live[-1], ca[-2])

    reach = clash_distance + 3.0
    midpoints = 0.5 * (ca[:-1] + ca[1:])
    ca_chain = np.concatenate([3 * np.arange(n_res) + 1, np.full(fixed.shape[0] - n_res, 10**9)])
    ca_oxygen = np.zeros(fixed.shape[0], dtype=np.int64)
    movable_chain = 3 * np.tile(np.arange(n_res), 3) + np.concatenate(
        [np.zeros(n_res, dtype=np.int64), np.full(n_res, 2), np.full(n_res, 2)]
    )
    movable_oxygen = np.concatenate(
        [np.zeros(2 * n_res, dtype=np.int64), np.ones(n_res, dtype=np.int64)]
    )
    kinds = (("C", 0), ("N", 1), ("O", 0))

    def owner(kind, offset):
        units = np.arange(n_units)
        position = 2 if kind in ("C", "O") else 0
        return 3 * (units + offset) + position, np.full(n_units, 1 if kind == "O" else 0)

    fixed_idx = np.empty((3, n_units, MAX_NEIGHBOURS), dtype=np.int64)
    for which, (kind, offset) in enumerate(kinds):
        chain, oxygen = owner(kind, offset)
        fixed_idx[which] = _neighbour_indices(
            fixed_tree, midpoints, reach, chain, oxygen, ca_chain, ca_oxygen, fixed.shape[0]
        )
    moving_idx = np.empty((3, n_units, MAX_NEIGHBOURS), dtype=np.int64)

    if rama_table is None:
        from .backbone_refine import _RAMA_TABLE

        rama_table = _RAMA_TABLE

    previous = None
    sweeps = 0
    converged = False
    for sweep in range(max_sweeps):
        live = np.vstack([n_live, c_live, o_live])
        live_padded = np.vstack([live, np.full((1, 3), _PAD_COORDINATE)])
        live_tree = cKDTree(live)
        for which, (kind, offset) in enumerate(kinds):
            chain, oxygen = owner(kind, offset)
            moving_idx[which] = _neighbour_indices(
                live_tree, midpoints, reach, chain, oxygen, movable_chain, movable_oxygen, 3 * n_res
            )
        sweeps = sweep + 1
        swept, largest = sweep_region(
            ca,
            n_live,
            c_live,
            o_live,
            azimuths,
            fixed_padded,
            live_padded,
            fixed_idx,
            moving_idx,
            rama_table,
            180.0 / (1.0 + sweep),
            candidates,
            angle_weight,
            clash_weight,
            rama_weight,
            clash_distance,
        )
        improvement = float("inf") if previous is None else previous - swept
        previous = swept
        if largest == 0.0 or improvement <= tolerance:
            converged = True
            break
    return n_live, c_live, o_live, azimuths, sweeps, converged
