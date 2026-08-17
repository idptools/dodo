"""Compiled inner loop for backbone refinement.

The refinement objective is cheap arithmetic evaluated an enormous number of times: for a
583-residue region, 397 peptide units x 25 candidate azimuths x 15 sweeps. In numpy that is
~150,000 array calls whose per-call dispatch cost dwarfs the six multiplies inside, and profiling
put about 60% of the time there rather than in the arithmetic. Batching can only amortise that
overhead, never remove it -- measured, batching ten models together recovered 3.6x of a
theoretical 10x, and batching over peptide units cost more in extra sweeps than it saved.

Compiling the loop removes the overhead outright. Measured over all six of p300's rebuilt regions
with the real obstacle set: **7.76 s to 0.49 s, 15.9x**.

The objective is compiled twice over, and the second time was not optional
-----------------------------------------------------------------------------
:func:`region_energy` exists because once the sweep was compiled, the *reporting* was the expensive
part. ``energy_before``/``energy_after`` are read once each per region, and the numpy scorer that
produced them evaluated every peptide unit through ``score_candidates``: profiled over p300's six
regions, 12 calls cost 0.562 s against the compiled sweep's 0.204 s of self time. Two readings cost
2.8x the optimisation they described.

Compiled, those same 12 calls cost 0.026 s -- and because nothing then reads the numpy neighbour
tables on this path, they are no longer built either. Measured over the six regions, 30 interleaved
reps each: refinement **971.9 ms to 507.4 ms, 1.92x** -- and 1478.4 to 767.9 ms, 1.93x, on a repeat
of the same experiment while the machine was busier, which is why the ratio is quoted rather than
the absolute times. The whole ``rebuild --backbone`` command gains 1.19x, and its 1.19 MB output
file is byte-identical before and after.

:func:`_unit_value` is why the two entry points cannot disagree: the sweep and the energy are the
same inlined body, so a term can only be in both or neither. Extracting it was verified to leave
the sweep bit-identical -- all six regions, every coordinate, azimuth, sweep count and convergence
flag, max abs difference 0.000e+00.

Equivalence to :mod:`dodo.construct.backbone_refine` is not assumed, it is established. Given
byte-identical inputs for one peptide unit, this kernel and the numpy scorer agree to:

* C, N and O placement -- 0.000e+00, bit-identical
* the Ramachandran term -- 0.000e+00
* the clash term -- 0.000e+00
* the two N-CA-C angle terms -- 2.8e-13 and 1.1e-13, from ``math.acos`` versus ``np.arccos``
  rounding in the last bits

and they select the same candidate. Summed over a whole region, :func:`region_energy` against the
numpy ``total_energy`` on frozen states -- the six real p300 regions, each both canonicalized and
refined, twelve states -- agrees to a worst absolute difference of 1.5e-11 on the angle term,
9.1e-13 on the Ramachandran term, 1.1e-13 on the clash term and 2.9e-11 on all three together, the
last being 2.8e-16 relative. On the 60-residue test fixture each isolated term is bit-identical and
only their sum differs, by 4.5e-13, because numpy reduces ``20 * square(r)`` where the kernel
evaluates ``20 * r * r``.

The tests that establish this live in ``tests/unit/test_backbone_kernel.py`` and they work by
comparing the objective as a PURE FUNCTION on identical inputs. That distinction matters: comparing
the end states of two independent runs cannot localise a difference, because a greedy search
amplifies any perturbation chaotically. Four earlier attempts to find a suspected bug that way all
failed, and there was no bug to find.

The consequence of that 1e-13 is worth stating plainly, because it is visible: across the ~4,800
decisions in a sweep it eventually flips one that is nearly balanced, after which the two backends
diverge and produce different -- equally good, measurably equivalent -- coordinates. So output is
not bit-identical between backends. It was never bit-identical across numpy versions or BLAS
builds either, for exactly the same reason.

The neighbour search is compiled too, and scipy is off the hot path
-------------------------------------------------------------------
:func:`_neighbour_indices_grid` replaced a ``cKDTree`` that sat in the middle of an otherwise
compiled path: built in Python, queried in Python, its index arrays handed to the sweep. Profiled at
0.199 s cumulative against the compiled sweep's 0.204 s of self time -- the search cost as much as
the optimisation it fed. A uniform cell list in nopython mode costs 0.011 s, **18x**, and takes
``refine_region`` 1.80x overall with ``sweep_region`` unchanged (68 calls, 0.219 vs 0.221 s),
which is what confirms only the search moved.

Two things the compiled version does better, not merely faster. It is a true RADIUS search, where
``cKDTree.query(k=48)`` can only ever examine the 48 nearest -- so the cap is now a checkable
property rather than a silent limit (see :func:`_check_neighbour_cap`, which found that the scipy
guard misses a case it should catch). And all three atom kinds share one scan, because the three
calls it replaced differed only in the covalent-separation filter and were repeating one identical
geometric search three times.

Levers that were measured and rejected
--------------------------------------
``fastmath`` is deliberately OFF, and the reason is not that it is slow. Measured, it is worth 1.14x
on the sweep loop -- but the loop is a fifth of ``refine_backbone``, so that dilutes to 1.2% on
a whole rebuild. Against that, its ``nnan`` licence lets the ``acos`` domain clamp swallow a NaN:
:func:`_angle` returns ``nan`` with it off and ``180.0`` with it on. That path is currently
unreachable, because ``refine_backbone`` rejects non-finite input, but it makes that guard
load-bearing for error behaviour -- a NaN that used to surface visibly would become a plausible
angle. 1.2% does not buy that in a geometry kernel.

Intel SVML is unavailable rather than unhelpful: no ``intel-cmplr-lib-rt`` wheel exists for arm64 on
any platform, and numba declines to use it here regardless. A faster threading layer is also a dead
end -- ``omp`` measured 0.86x against the default ``workqueue``, and neither matters while this code
is serial. ``prange`` over the MODEL axis does scale (16.4x at 18 models), but it needs the pipeline
to collect every model's CA trace before one batched pass, so it is not wired up here.
"""

from __future__ import annotations

import math

import numba as nb
import numpy as np

from ..exceptions import GeometryError

PC = math.radians(20.4118)
PN = math.radians(14.8941)
CA_C = 1.525
N_CA = 1.458
C_O = 1.231
ACO = math.radians(120.8)
NCAC = 111.0
BIN = 30.0
# Hard N-CA-C window (dodo.constants.N_CA_C_WINDOW_MIN/MAX; a test pins the copies equal). The
# penalty is a step, not a slope: it only has to dominate any clash saving a candidate could buy
# by collapsing the angle -- the worst plausible clash sum is a few thousand -- so that the
# argmin never picks an angle the bond validator would flag as two atoms on top of each other.
ANGLE_LO = 80.0
ANGLE_HI = 160.0
ANGLE_PENALTY = 1.0e5


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


@nb.njit(cache=True, fastmath=False, inline="always")
def _unit_value(
    unit,
    az,
    ca,
    n_live,
    c_live,
    fixed_pad,
    moving_pad,
    fixed_idx,
    moving_idx,
    rama,
    angle_w,
    clash_w,
    rama_w,
    cd,
    limit,
):
    """Everything one azimuth can change, for one peptide unit -- the objective's whole body.

    Extracted so :func:`sweep_region` and :func:`region_energy` cannot drift apart. They are the
    optimiser and the number it reports, so a term present in one and not the other would make
    ``energy_before``/``energy_after`` describe something other than what was minimised. Inlined,
    so the sweep pays nothing for the indirection -- verified bit-identical to the version that had
    this body written out inside the candidate loop.
    """
    n_res = ca.shape[0]
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
        if a < ANGLE_LO or a > ANGLE_HI:
            val += ANGLE_PENALTY
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
        if a < ANGLE_LO or a > ANGLE_HI:
            val += ANGLE_PENALTY
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
    return val


@nb.njit(cache=True, fastmath=False)
def region_energy(  # noqa: D103 - documented at module level; njit rejects a docstring-only body
    ca,
    n_live,
    c_live,
    azimuths,
    fixed_pad,
    moving_pad,
    fixed_idx,
    moving_idx,
    rama,
    angle_w,
    clash_w,
    rama_w,
    cd,
):
    n_units = ca.shape[0] - 1
    limit = cd * cd
    total = 0.0
    for unit in range(n_units):
        total += _unit_value(
            unit,
            azimuths[unit],
            ca,
            n_live,
            c_live,
            fixed_pad,
            moving_pad,
            fixed_idx,
            moving_idx,
            rama,
            angle_w,
            clash_w,
            rama_w,
            cd,
            limit,
        )
    # Halved for the same reason the numpy scorer halves it: each unit's score covers the residues
    # on both of its sides, so summing over units double-counts every shared term.
    return total / 2.0


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
            val = _unit_value(
                unit,
                az,
                ca,
                n_live,
                c_live,
                fixed_pad,
                moving_pad,
                fixed_idx,
                moving_idx,
                rama,
                angle_w,
                clash_w,
                rama_w,
                cd,
                limit,
            )
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


#: Neighbours retained per movable atom. MEASURED, twice. Over EVERY clash call of a full p300
#: rebuild (222 calls, 67,791 centres, all six regions, against the real obstacle set) the kept
#: sets run to a median of 3 and a maximum of 41, with the in-shell count before the covalent
#: separation filter peaking at 43 -- which made 48 look like 5 points of headroom. Then the
#: 117-structure corpus ran with the backbone ON (2026-08-13) and a crowded region of Q9C000 put
#: **54** fixed atoms in one unit's shell, so 48 was an overflow crash on real input, not a margin.
#: 96 is ~1.8x that new worst. The cost of the size is small and bounded: pad slots point at
#: :data:`_PAD_COORDINATE`, whose clamped gap is exactly zero, so extra rows cost one distance
#: computation each in the scan -- and the refinement is a single-digit share of ``--backbone``
#: wall time. Truncating below the real count would silently change the objective, so the cap is
#: checked, never trusted (:func:`_check_neighbour_cap`); a region that beats even this cap falls
#: back to the capless numpy path rather than failing the rebuild (see ``refine_backbone``).
MAX_NEIGHBOURS: int = 96

#: Coordinate for a padded neighbour slot, in Angstroms. Far enough that it can never clash.
_PAD_COORDINATE: float = 1.0e6

#: Ceiling on cells in the neighbour grid, which is dense over the bounding box and so costs
#: ``cells * 8`` bytes. MEASURED worst case over a full p300 rebuild: 47x55x67 = 173,195 cells,
#: 1.39 MB, for 10,022 points -- sparse (0.06 points per cell), because a chain fills little of its
#: own box, but cheap. The cost grows with bounding-box VOLUME, and DODO generates extended IDRs, so
#: a conformation spanning several times p300's ~390 A would grow this cubically. Rather than let
#: it, :func:`_cell_list` enlarges the cell until the count fits. 4,000,000 cells is 32 MB, about
#: 23x p300's worst case, so nothing real is expected to reach it.
_MAX_GRID_CELLS: int = 4_000_000


def _check_neighbour_cap(worst: int) -> None:
    """Raise if any centre had more qualifying neighbours than a row can hold.

    :func:`_neighbour_indices_grid` counts qualifying neighbours WITHOUT a ceiling, so ``worst`` is
    the true requirement even for a row that overflowed. Truncating would change the objective
    rather than merely slow it down, so this is checked, never trusted.

    This is strictly stronger than the guard :func:`_neighbour_indices` can manage. That one raises
    only when the 48th *nearest* point survives the separation filter, so a case where the 48th
    nearest is covalently adjacent passes while qualifying points are silently dropped. Constructed
    and confirmed: 60 points inside reach with the 48th-nearest made adjacent kept 47 of 59 and did
    not raise. Consequence worth knowing: a crowded region that previously completed with a quietly
    truncated objective now aborts instead.
    """
    if worst > MAX_NEIGHBOURS:
        raise GeometryError(
            f"A peptide unit has {worst} neighbours within the clash shell but the compiled "
            f"kernel's rows hold {MAX_NEIGHBOURS}, so the objective would be truncated. Raise "
            f"MAX_NEIGHBOURS or use backend='numpy'."
        )


@nb.njit(cache=True, fastmath=False)
def _cell_list(pts, reach, max_cells):
    """Bin ``pts`` into a uniform grid, CSR-style. Returns ``(lo, dims, starts, order, cell)``.

    ``order`` lists point indices grouped by cell and ``starts[c]:starts[c + 1]`` slices the group
    for cell ``c``. Two passes, no per-point allocation, nothing ragged -- numba has no dict or
    list-of-lists worth using here.

    The cell side starts at ``reach`` and is enlarged until the dense ``starts`` array fits inside
    ``max_cells`` (see :data:`_MAX_GRID_CELLS`). Enlarging is always SAFE, never merely tolerable:
    the caller scans the 3x3x3 block around a query's cell, and with a side of at least ``reach``
    any point within ``reach`` of the query necessarily falls in that block. A bigger cell scans
    more points per cell, so it costs time, but it cannot change which neighbours are found.
    """
    n = pts.shape[0]
    lo = np.empty(3, np.float64)
    hi = np.empty(3, np.float64)
    for d in range(3):
        lo[d] = pts[0, d]
        hi[d] = pts[0, d]
    for i in range(1, n):
        for d in range(3):
            v = pts[i, d]
            if v < lo[d]:
                lo[d] = v
            elif v > hi[d]:
                hi[d] = v
    cell = reach
    dims = np.empty(3, np.int64)
    while True:
        for d in range(3):
            dims[d] = int((hi[d] - lo[d]) / cell) + 1
        ncell = dims[0] * dims[1] * dims[2]
        if ncell <= max_cells:
            break
        # Cube-root scaling reaches the target in one or two rounds rather than creeping up on it.
        cell *= (ncell / max_cells) ** (1.0 / 3.0)
    starts = np.zeros(ncell + 1, np.int64)
    cid = np.empty(n, np.int64)
    for i in range(n):
        ix = int((pts[i, 0] - lo[0]) / cell)
        iy = int((pts[i, 1] - lo[1]) / cell)
        iz = int((pts[i, 2] - lo[2]) / cell)
        # Guards the top edge only, where (hi - lo) / cell can round up to dims.
        if ix >= dims[0]:
            ix = dims[0] - 1
        if iy >= dims[1]:
            iy = dims[1] - 1
        if iz >= dims[2]:
            iz = dims[2] - 1
        c = (ix * dims[1] + iy) * dims[2] + iz
        cid[i] = c
        starts[c + 1] += 1
    for c in range(ncell):
        starts[c + 1] += starts[c]
    order = np.empty(n, np.int64)
    fill = starts[:ncell].copy()
    for i in range(n):
        c = cid[i]
        order[fill[c]] = i
        fill[c] += 1
    return lo, dims, starts, order, cell


@nb.njit(cache=True, fastmath=False)
def _neighbour_indices_grid(
    pts, centres, reach, own_chain, own_oxygen, chain_of, oxygen_of, sentinel, cap, max_cells
):
    """Compiled replacement for :func:`_neighbour_indices`, for all three atom kinds at once.

    A true radius search: the cell side is at least ``reach``, so every point within ``reach`` of a
    centre lies in the 3x3x3 block of cells around it and one scan finds them all. ``cKDTree.query``
    with a fixed ``k`` could only ever examine the ``k`` nearest.

    All three kinds (rows of ``own_chain``) share one scan. The three calls this replaces shared
    ``centres`` and ``reach`` and differed only in the covalent-separation filter, so they were
    repeating one identical geometric search three times.

    Returns ``(out, worst)``. ``worst`` is the largest TRUE kept count over all centres, counted
    without a ceiling even once a row is full, so :func:`_check_neighbour_cap` can raise on overflow
    instead of the caller silently truncating.
    """
    nk = own_chain.shape[0]
    nc = centres.shape[0]
    out = np.full((nk, nc, cap), sentinel, np.int64)
    worst = 0
    n = pts.shape[0]
    if n == 0 or nc == 0:
        return out, worst
    lo, dims, starts, order, cell = _cell_list(pts, reach, max_cells)
    reach2 = reach * reach
    # Kept neighbours per kind for the current centre, ordered by ascending squared distance.
    d2buf = np.empty((nk, cap), np.float64)
    idbuf = np.empty((nk, cap), np.int64)
    held = np.zeros(nk, np.int64)  # entries actually stored (<= cap)
    total = np.zeros(nk, np.int64)  # entries that qualified (uncapped)
    for q in range(nc):
        qx = centres[q, 0]
        qy = centres[q, 1]
        qz = centres[q, 2]
        for w in range(nk):
            held[w] = 0
            total[w] = 0
        cx0 = int((qx - lo[0]) / cell)
        cy0 = int((qy - lo[1]) / cell)
        cz0 = int((qz - lo[2]) / cell)
        # A centre outside the point cloud clamps to the rim cell; since the side is at least
        # `reach`, the +-1 block around the rim still covers everything that could be within reach.
        if cx0 < 0:
            cx0 = 0
        elif cx0 >= dims[0]:
            cx0 = dims[0] - 1
        if cy0 < 0:
            cy0 = 0
        elif cy0 >= dims[1]:
            cy0 = dims[1] - 1
        if cz0 < 0:
            cz0 = 0
        elif cz0 >= dims[2]:
            cz0 = dims[2] - 1
        xlo = cx0 - 1 if cx0 > 0 else 0
        xhi = cx0 + 2 if cx0 + 2 <= dims[0] else dims[0]
        ylo = cy0 - 1 if cy0 > 0 else 0
        yhi = cy0 + 2 if cy0 + 2 <= dims[1] else dims[1]
        zlo = cz0 - 1 if cz0 > 0 else 0
        zhi = cz0 + 2 if cz0 + 2 <= dims[2] else dims[2]
        for cx in range(xlo, xhi):
            for cy in range(ylo, yhi):
                base = (cx * dims[1] + cy) * dims[2]
                for cz in range(zlo, zhi):
                    c = base + cz
                    for s in range(starts[c], starts[c + 1]):
                        i = order[s]
                        dx = qx - pts[i, 0]
                        dy = qy - pts[i, 1]
                        dz = qz - pts[i, 2]
                        d2 = dx * dx + dy * dy + dz * dz
                        if d2 > reach2:
                            continue
                        for w in range(nk):
                            sep = (
                                abs(own_chain[w, q] - chain_of[i]) + own_oxygen[w, q] + oxygen_of[i]
                            )
                            if sep <= 4:
                                continue
                            total[w] += 1
                            k = held[w]
                            if k == cap:
                                continue
                            # Insertion sort, ascending d2, stable on ties. This reproduces
                            # cKDTree's by-distance ordering, which keeps the clash sum's
                            # accumulation order -- and so its floating-point result -- identical.
                            j = k
                            while j > 0 and d2buf[w, j - 1] > d2:
                                d2buf[w, j] = d2buf[w, j - 1]
                                idbuf[w, j] = idbuf[w, j - 1]
                                j -= 1
                            d2buf[w, j] = d2
                            idbuf[w, j] = i
                            held[w] = k + 1
        for w in range(nk):
            if total[w] > worst:
                worst = total[w]
            for j in range(held[w]):
                out[w, q, j] = idbuf[w, j]
    return out, worst


def _neighbour_indices(tree, centres, reach, own_chain, own_oxygen, chain_of, oxygen_of, sentinel):
    """Return padded neighbour indices using ``cKDTree``. NOT on the hot path any more.

    Retained deliberately, as the oracle :func:`_neighbour_indices_grid` is checked against -- see
    ``tests/unit/test_backbone_kernel.py::test_the_compiled_neighbour_search_matches_the_kdtree``.
    Note the fixed ``k``: this can only ever examine the 48 nearest points, so it agrees with a true
    radius search exactly when at most 48 lie inside the shell. The measured peak over a full p300
    rebuild is 43, so it agreed there; the compiled version is not limited that way.
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


#: The three atoms one azimuth moves, as ``(kind, residue offset from the unit index)``. Order fixes
#: the meaning of the leading axis of the ``fixed_idx``/``moving_idx`` tables: 0 is the unit's C, 1
#: its N, 2 its O.
_MOVABLE_KINDS: tuple[tuple[str, int], ...] = (("C", 0), ("N", 1), ("O", 0))


class RegionKernel:
    """Compiled refinement for ONE region, with its neighbour tables built once.

    Both compiled entry points -- the sweep and the energy -- read the same tables from here, which
    is the point of the class rather than two functions. The energy is what ``energy_before`` and
    ``energy_after`` report, so evaluating it against different neighbour sets than the sweep
    optimised would make those numbers describe a different objective. Sharing the tables makes
    that impossible instead of merely unlikely, and saves rebuilding a KD-tree over every obstacle
    in the structure -- 12,517 of them on p300's last region -- once per energy reading.

    The tree work stays in Python because scipy's KD-tree cannot cross into nopython mode.
    """

    def __init__(
        self,
        ca: np.ndarray,
        *,
        obstacles: np.ndarray | None = None,
        clash_distance: float = 2.9,
        angle_weight: float = 0.124,
        clash_weight: float = 40.0,
        rama_weight: float = 20.0,
        rama_table: np.ndarray | None = None,
    ) -> None:
        if rama_table is None:
            from .backbone_refine import _RAMA_TABLE

            rama_table = _RAMA_TABLE
        self.ca = np.ascontiguousarray(ca, dtype=np.float64)
        self.clash_distance = float(clash_distance)
        self.angle_weight = float(angle_weight)
        self.clash_weight = float(clash_weight)
        self.rama_weight = float(rama_weight)
        self.rama_table = rama_table

        n_res = self.ca.shape[0]
        n_units = n_res - 1
        fixed = self.ca if obstacles is None else np.vstack([self.ca, obstacles])
        self._fixed_padded = np.vstack([fixed, np.full((1, 3), _PAD_COORDINATE)])
        self._reach = clash_distance + 3.0
        self._midpoints = 0.5 * (self.ca[:-1] + self.ca[1:])
        self._n_res = n_res

        ca_chain = np.concatenate(
            [3 * np.arange(n_res) + 1, np.full(fixed.shape[0] - n_res, 10**9)]
        )
        ca_oxygen = np.zeros(fixed.shape[0], dtype=np.int64)
        self._movable_chain = 3 * np.tile(np.arange(n_res), 3) + np.concatenate(
            [np.zeros(n_res, dtype=np.int64), np.full(n_res, 2), np.full(n_res, 2)]
        )
        self._movable_oxygen = np.concatenate(
            [np.zeros(2 * n_res, dtype=np.int64), np.ones(n_res, dtype=np.int64)]
        )
        units = np.arange(n_units)
        self._owners = tuple(
            (
                3 * (units + offset) + (2 if kind in ("C", "O") else 0),
                np.full(n_units, 1 if kind == "O" else 0),
            )
            for kind, offset in _MOVABLE_KINDS
        )

        # One scan for all three kinds; they share `centres` and `reach` and differ only in the
        # covalent-separation filter, so three separate searches were repeating identical geometry.
        self._own_chain = np.ascontiguousarray([chain for chain, _ in self._owners])
        self._own_oxygen = np.ascontiguousarray([oxygen for _, oxygen in self._owners])
        self._fixed_idx, worst = _neighbour_indices_grid(
            np.ascontiguousarray(fixed, dtype=np.float64),
            self._midpoints,
            self._reach,
            self._own_chain,
            self._own_oxygen,
            ca_chain,
            ca_oxygen,
            fixed.shape[0],
            MAX_NEIGHBOURS,
            _MAX_GRID_CELLS,
        )
        _check_neighbour_cap(worst)

    def _moving_tables(
        self, n_live: np.ndarray, c_live: np.ndarray, o_live: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Rebuild the self-avoidance tables for the CURRENT generated atoms.

        Unlike the fixed ones these cannot be precomputed: the atoms a unit could hit move as
        refinement proceeds.
        """
        live = np.ascontiguousarray(np.vstack([n_live, c_live, o_live]), dtype=np.float64)
        live_padded = np.vstack([live, np.full((1, 3), _PAD_COORDINATE)])
        moving_idx, worst = _neighbour_indices_grid(
            live,
            self._midpoints,
            self._reach,
            self._own_chain,
            self._own_oxygen,
            self._movable_chain,
            self._movable_oxygen,
            3 * self._n_res,
            MAX_NEIGHBOURS,
            _MAX_GRID_CELLS,
        )
        _check_neighbour_cap(worst)
        return live_padded, moving_idx

    def energy(
        self,
        n_live: np.ndarray,
        c_live: np.ndarray,
        o_live: np.ndarray,
        azimuths: np.ndarray,
    ) -> float:
        """Return the whole-region objective at the given state: the compiled ``total_energy``.

        The numpy scorer this replaces evaluates every peptide unit through ``score_candidates``,
        one array call per term per unit, and it is called twice per region for reporting only. That
        made the *reporting* cost 2.8x the compiled sweep it was reporting on.

        Note what it does and does not read. C, N and O are re-placed from ``azimuths`` exactly as a
        candidate would be, not taken from ``c_live``/``n_live``/``o_live``; those are read for the
        neighbouring units' atoms and for the self-avoidance tables. That is what the numpy version
        does too, and it is why passing a state whose azimuths disagree with its coordinates gives
        the energy of the azimuths.
        """
        live_padded, moving_idx = self._moving_tables(n_live, c_live, o_live)
        return float(
            region_energy(
                self.ca,
                np.ascontiguousarray(n_live),
                np.ascontiguousarray(c_live),
                np.ascontiguousarray(azimuths),
                self._fixed_padded,
                live_padded,
                self._fixed_idx,
                moving_idx,
                self.rama_table,
                self.angle_weight,
                self.clash_weight,
                self.rama_weight,
                self.clash_distance,
            )
        )

    def sweep_to_convergence(
        self,
        n_live: np.ndarray,
        c_live: np.ndarray,
        o_live: np.ndarray,
        azimuths: np.ndarray,
        *,
        max_sweeps: int = 30,
        tolerance: float = 0.25,
        candidates: int = 24,
    ) -> tuple[int, bool]:
        """Run coordinate descent IN PLACE on the given state. Returns ``(sweeps, converged)``."""
        previous = None
        sweeps = 0
        converged = False
        for sweep in range(max_sweeps):
            live_padded, moving_idx = self._moving_tables(n_live, c_live, o_live)
            sweeps = sweep + 1
            swept, largest = sweep_region(
                self.ca,
                n_live,
                c_live,
                o_live,
                azimuths,
                self._fixed_padded,
                live_padded,
                self._fixed_idx,
                moving_idx,
                self.rama_table,
                180.0 / (1.0 + sweep),
                candidates,
                self.angle_weight,
                self.clash_weight,
                self.rama_weight,
                self.clash_distance,
            )
            improvement = float("inf") if previous is None else previous - swept
            previous = swept
            if largest == 0.0 or improvement <= tolerance:
                converged = True
                break
        return sweeps, converged


def seed_region(
    ca: np.ndarray, c_xyz: np.ndarray, n_xyz: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Canonicalize a starting backbone into ``(n_live, c_live, o_live, azimuths)``.

    Bit-identical to the seeding :func:`~dodo.construct.backbone_refine.refine_backbone` does with
    ``_place_unit``/``_place_oxygen``; the equivalence test pins that down.
    """
    from .backbone_refine import _azimuth_frame
    from .ca_backbone import _terminal_oxygen

    n_units = ca.shape[0] - 1
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
    return n_live, c_live, o_live, azimuths


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

    Returns ``(n_xyz, c_xyz, o_xyz, azimuths, sweeps, converged)``. Seeds, then sweeps; a caller
    that already holds a seeded state should build a :class:`RegionKernel` and drive it directly, so
    that the seeding is not repeated.
    """
    refiner = RegionKernel(
        ca,
        obstacles=obstacles,
        clash_distance=clash_distance,
        angle_weight=angle_weight,
        clash_weight=clash_weight,
        rama_weight=rama_weight,
        rama_table=rama_table,
    )
    n_live, c_live, o_live, azimuths = seed_region(ca, c_xyz, n_xyz)
    sweeps, converged = refiner.sweep_to_convergence(
        n_live,
        c_live,
        o_live,
        azimuths,
        max_sweeps=max_sweeps,
        tolerance=tolerance,
        candidates=candidates,
    )
    return n_live, c_live, o_live, azimuths, sweeps, converged
