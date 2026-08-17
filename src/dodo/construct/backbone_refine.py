"""Refine a placed backbone: bond angles, dihedrals and steric clashes, alpha carbons held fixed.

Why the alpha carbons do not move
---------------------------------
They are the part DODO is confident about. The walk engine placed them at exactly 3.81 A with the
pseudo-angle inside its measured window, self-avoiding, and steered to the ALBATROSS-predicted
end-to-end distance -- which is the scientific claim the whole package makes. Letting a refinement
nudge them would trade that claim for backbone cosmetics. So they are fixed, and only N, C and O
move.

The parameterization, and why it makes this small
-------------------------------------------------
With both flanking alpha carbons fixed and ideal internal geometry, the peptide unit between them
has exactly **one** degree of freedom: an azimuth about the CA-CA axis. Derived from DODO's ideal
bond lengths and angles:

* C sits :data:`~dodo.constants.CA_C_BOND_LENGTH` from CA(i), at :data:`_POLAR_C` degrees off the
  CA(i)->CA(i+1) direction, at azimuth ``psi``;
* N sits :data:`~dodo.constants.N_CA_BOND_LENGTH` from CA(i+1), at :data:`_POLAR_N` degrees off the
  CA(i+1)->CA(i) direction, at azimuth ``psi + 180``.

So the whole backbone of an ``n``-residue region is described by ``n - 1`` angles, and all three
bonds are exact for *any* azimuth. Bond lengths therefore need no refinement at all, which is what
makes this a well-posed problem in one parameter per unit rather than a 3N-dimensional minimization
that can drift anywhere.

The two polar angles are solved *for DODO's* :data:`~dodo.constants.CA_CA_BOND_LENGTH`, not taken
from a table of ideal internal coordinates. That distinction cost 0.004 A: textbook peptide geometry
determines a CA-CA distance of 3.8040 A, so a unit built from it and then stretched onto a 3.81 A
trace put the discrepancy in the peptide bond, giving 1.3334 A instead of 1.329. Solving the other
way round -- fixing all three bonds and the 3.81 A span, and letting the two *bond angles* absorb
the difference -- costs 0.37 and 0.39 degrees on CA-C-N and C-N-CA, which is nothing next to the 2-3
degrees those angles vary by in real backbone, and buys exact bonds. 3.81 A is itself well founded:
19,500 consecutive alpha carbons in all-atom IDR simulation give a mean of 3.8129 and a median of
3.8112 A.

**That exactness is inherited from the input spacing, and it is worth being precise about why.**
C is placed off CA(i) and N off CA(i+1), so the gap between them -- the peptide bond -- absorbs any
deviation of ``|CA(i) - CA(i+1)|`` from the canonical value, one for one. Measured on a real
simulation frame whose alpha carbons span 3.747-3.877 A, the resulting C-N bond spans 1.288-1.383 A,
correlating with the input spacing at r = 1.0000. This is not a defect in the refinement: a rigid
peptide unit *determines* the CA-CA distance, so a trace that disagrees with it cannot be fitted by
one without something giving, and the peptide bond is the least bad place to put the discrepancy.

In DODO it does not arise. The walk engine emits exactly 3.81 A by construction, and STARLING
output is projected onto exactly 3.81 A before use. It arises for a caller passing an experimental
or unregularized trace straight in, so :func:`refine_backbone` measures the input spacing and says
so in its notes rather than quietly returning stretched bonds. Regularize first if you see it.

What is actually refined
------------------------
Three things the per-unit placement cannot see, because each involves more than one unit or more
than the backbone:

**The N-CA-C angle.** Unit ``i-1`` places N(i) and unit ``i`` places C(i), so the angle between them
at CA(i) is set by two independent azimuths and constrained by neither. Measured against 100 frames
of all-atom simulation, the unrefined backbone gets its centre right -- 110.38 degrees against a
real 109.90 -- but its spread is **3.3x too wide**, 9.29 degrees against 2.84. That is real
information the per-unit table never used.

**Phi and psi.** Backbone dihedrals are a property of consecutive units. A CA trace does not
determine them, but it does not licence every value either: a placement can land in a sterically
forbidden region of the Ramachandran plot, and nothing in the per-unit placement notices.

**Steric clashes.** Newly placed N, C and O atoms can collide with each other and with atoms DODO
did not build -- folded-domain side chains in particular, which the alpha-carbon stage never had to
avoid because it was only placing alpha carbons.

Refinement is a bounded coordinate descent over the azimuths: deterministic, no random restarts, and
no azimuth can change a bond length -- the three bonds inside a unit are fixed by the unit's own
geometry, whatever it is rotated to. The one thing that does move them is the canonicalization on
entry, and only when the input alpha carbons are not uniformly spaced, as described above.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Final

import numpy as np
from scipy.spatial.distance import cdist

from ..constants import (
    C_N_PEPTIDE_BOND_LENGTH,
    C_O_BOND_LENGTH,
    CA_C_BOND_LENGTH,
    CA_C_O_ANGLE,
    CA_CA_BOND_LENGTH,
    N_CA_BOND_LENGTH,
    N_CA_C_ANGLE,
    N_CA_C_WINDOW_MAX,
    N_CA_C_WINDOW_MIN,
)
from ..exceptions import GeometryError, InvalidParameterError, MissingDependencyError
from .ca_backbone import _terminal_carbon, _terminal_nitrogen, _terminal_oxygen

__all__ = [
    "RefinementResult",
    "refine_backbone",
]

#: Polar angle of C off the CA(i)->CA(i+1) direction, in degrees.
#:
#: DERIVED from DODO's ideal peptide geometry, not fitted: build CA-C, C-N and N-CA at their ideal
#: lengths with the ideal CA-C-N and C-N-CA angles, planar and trans, and this angle falls out --
#: along with a CA-CA separation of 3.8040 A, which is DODO's 3.81 to within 0.006 A. That agreement
#: is a check on the constants rather than a coincidence.
_POLAR_C: Final[float] = 20.4118

#: Polar angle of N off the CA(i+1)->CA(i) direction, in degrees. DERIVED with :data:`_POLAR_C` by
#: solving the planar trans unit for exact bonds at :data:`~dodo.constants.CA_CA_BOND_LENGTH`.
_POLAR_N: Final[float] = 14.8941

#: Azimuthal separation of C and N about the CA-CA axis, in degrees. Exactly 180 for a planar trans
#: unit. Real units measure 163.1 +/- 12.7, so this is an idealization -- and a deliberate one: it
#: is what reduces the unit to a single free parameter.
_AZIMUTH_SEPARATION: Final[float] = 180.0

#: Step penalty for an N-CA-C angle outside the hard window
#: (:data:`~dodo.constants.N_CA_C_WINDOW_MIN`/``MAX``). The soft angle term (weight 0.124) is
#: deliberately weak -- tightening it was measured to hurt accuracy -- which leaves it able to lose
#: to a strong clash term: measured on the corpus (2026-08-13), two structures traded the angle
#: down to ~79 degrees, putting a residue's own N and C 1.90 A apart, which the bond validator
#: rightly reports as two atoms on top of each other. The step only has to dominate any clash
#: saving (a few thousand at worst), so candidates inside the window are still ranked purely by
#: the tuned objective. The compiled kernel applies the identical step (``ANGLE_PENALTY`` in
#: backbone_kernel.py); a test pins the two equal.
_ANGLE_WINDOW_PENALTY: Final[float] = 1.0e5

#: Ramachandran penalty, MEASURED from the same 100 frames of all-atom IDR simulation as the
#: peptide-unit table, over 19,602 phi/psi pairs at 30 degrees, circularly smoothed 3x3.
#:
#: Indexed ``[phi_bin][psi_bin]``. Zero in the most populated cell and rising as occupancy falls,
#: capped at 3.0 so an unusual dihedral is discouraged rather than forbidden -- IDRs do visit
#: unusual dihedrals, which is the point of them.
#:
#: This replaces four hand-placed circles taken from literature centres for FOLDED proteins, and
#: the difference matters. Those emphasised the alpha basin; measured, these IDRs are dominated by
#: extended/polyproline-II, median phi -106.9 and psi +122.3, and the table's minimum sits exactly
#: there. Glycine is a distinct population (median phi -49.5, psi +53.7) and is not separated out
#: here -- a per-residue-type table is the obvious refinement.
#:
#: 93.4% of the measured pairs fall in cells with a penalty below 1.0.
_RAMA_BIN_WIDTH: Final[int] = 30
_RAMA_PENALTY: Final[tuple[tuple[float, ...], ...]] = (
    (0.97, 1.73, 1.41, 1.02, 0.98, 1.08, 1.48, 1.36, 0.83, 0.47, 0.39, 0.52),
    (0.68, 1.52, 1.23, 0.76, 0.69, 0.76, 1.16, 1.18, 0.57, 0.19, 0.10, 0.22),
    (0.58, 1.48, 1.19, 0.69, 0.60, 0.67, 1.10, 1.18, 0.50, 0.10, 0.00, 0.12),
    (0.73, 1.64, 1.34, 0.82, 0.73, 0.81, 1.28, 1.45, 0.65, 0.23, 0.13, 0.25),
    (1.09, 2.06, 1.69, 1.16, 1.09, 1.20, 1.82, 2.00, 1.02, 0.56, 0.47, 0.59),
    (2.25, 3.00, 2.46, 2.01, 1.99, 2.17, 3.00, 3.00, 1.86, 1.38, 1.32, 1.47),
    (3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 2.79, 2.75, 2.81, 3.00),
    (2.75, 2.81, 2.86, 2.95, 3.00, 3.00, 2.89, 2.50, 2.28, 2.25, 2.36, 2.60),
    (2.43, 2.38, 2.40, 2.50, 2.68, 2.86, 2.63, 2.31, 2.07, 2.00, 2.10, 2.28),
    (2.22, 2.18, 2.18, 2.23, 2.39, 2.60, 2.45, 2.19, 1.98, 1.92, 1.98, 2.09),
    (2.08, 2.10, 2.01, 2.03, 2.13, 2.49, 2.42, 2.20, 1.93, 1.81, 1.81, 1.90),
    (1.52, 1.95, 1.74, 1.47, 1.51, 1.67, 2.14, 1.84, 1.36, 1.05, 0.98, 1.11),
)


@dataclass(frozen=True, slots=True)
class RefinementResult:
    """Outcome of refining a backbone.

    Attributes
    ----------
    n_xyz, c_xyz, o_xyz
        Refined coordinates. Alpha carbons are unchanged and are not returned.
    azimuths
        The refined azimuth of each peptide unit, in degrees.
    sweeps
        How many coordinate-descent sweeps ran.
    converged
        True if a sweep ended with no azimuth moving more than ``tolerance``.
    energy_before, energy_after
        The objective, for reporting. Lower is better; the units are arbitrary.
    notes
        What changed, in words.
    """

    n_xyz: np.ndarray
    c_xyz: np.ndarray
    o_xyz: np.ndarray
    azimuths: np.ndarray
    sweeps: int
    converged: bool
    energy_before: float
    energy_after: float
    notes: tuple[str, ...]


def _clash_penalty(
    points: np.ndarray, others: np.ndarray, clash_distance: float, weight: float
) -> np.ndarray:
    """Return the summed squared overlap of each point against ``others`` that are too close.

    ``cdist`` with ``sqeuclidean``, not a broadcast ``norm``, and the difference is larger than it
    looks: the broadcast form materialises a ``(n_points, n_others, 3)`` array of differences before
    reducing it, while cdist writes straight into the ``(n_points, n_others)`` result. Measured at
    the shapes this is actually called with -- 25 candidates against a median of 2 neighbours and a
    90th percentile of 6 -- it is 1.6x faster at 2 neighbours and 4.3x at 50, for bit-identical
    output.

    Clamping the squared distance before the square root, rather than clipping the gap after it, is
    what keeps a single expression correct: beyond the cutoff the clamp makes the gap exactly zero.

    A KD-tree would be the obvious thing here and is the wrong thing: building one per call costs
    16-25 us against 4, because the neighbour sets are tiny. The tree earns its keep choosing those
    neighbour sets once per sweep, which is where it is used.
    """
    limit = clash_distance * clash_distance
    squared = cdist(points, others, metric="sqeuclidean")
    np.minimum(squared, limit, out=squared)
    gaps = clash_distance - np.sqrt(squared)
    summed: np.ndarray = weight * np.einsum("ij,ij->i", gaps, gaps)
    return summed


def _cross(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Cross product of last-axis-3 arrays, written out by component.

    Not :func:`numpy.cross`, and the reason is measured rather than stylistic. ``np.cross`` supports
    arbitrary axes, and pays for it on every call with ``normalize_axis_tuple`` and ``moveaxis``
    -- pure Python bookkeeping that dwarfs six multiplies when the arrays are small. Profiling one
    583-residue region showed 219,816 calls to ``np.cross`` costing 6.0 s cumulative, with 1.3
    million calls to ``normalize_axis_tuple`` underneath it doing nothing but validate axes that
    were never anything other than -1.
    """
    ax, ay, az = a[..., 0], a[..., 1], a[..., 2]
    bx, by, bz = b[..., 0], b[..., 1], b[..., 2]
    return np.stack((ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx), axis=-1)


def _norm(v: np.ndarray) -> np.ndarray:
    """Euclidean length along the last axis, keeping the dimension.

    ``np.linalg.norm`` carries the same per-call overhead as ``np.cross`` for the same reason.
    """
    lengths: np.ndarray = np.sqrt(np.einsum("...i,...i->...i", v, v).sum(-1, keepdims=True))
    return lengths


def _unit(v: np.ndarray) -> np.ndarray:
    # np.maximum rather than a `where=` mask with a zeros_like buffer: the buffer was allocated on
    # every one of 193,623 calls per region, and a zero-length vector normalises to zero either way
    # because the numerator is already zero.
    norm = _norm(v)
    normalized: np.ndarray = v / np.maximum(norm, 1e-12)
    return normalized


def _azimuth_frame(ca_a: np.ndarray, ca_b: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Axis plus two perpendiculars defining where azimuth zero points, for one peptide unit."""
    axis = _unit(ca_b - ca_a)
    seed = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(seed, axis))) > 0.9:
        seed = np.array([0.0, 1.0, 0.0])
    first = _unit(seed - float(np.dot(seed, axis)) * axis)
    return axis, first, _cross(axis, first)


#: Position of each backbone atom along the main chain, in bonds from that residue's nitrogen.
#: The backbone is a linear graph -- N(i)-CA(i)-C(i)-N(i+1)-CA(i+1)... -- with each carbonyl oxygen
#: pendant on its own C, which is all :func:`_bond_separation` needs to count bonds.
_CHAIN_POSITION: Final[dict[str, int]] = {"N": 0, "CA": 1, "C": 2}


def _bond_separation(kind_a: str, residue_a: int, kind_b: str, residue_b: int) -> int:
    """Covalent bonds between two backbone atoms, counted along the backbone graph.

    Exists because the clash term needs to tell "close because they are bonded" from "close
    because they are colliding", and a residue-window rule cannot. Excluding everything within two
    residues -- which is what this replaced -- throws away real contacts: two carbonyl oxygens on
    ADJACENT residues are five bonds apart and perfectly free to collide, and on p300 one such
    pair ended up 0.6 A apart with nothing in the objective looking at it. It also over-excluded in
    the other direction, hiding genuine 1-4 and 1-5 contacts up to two residues away.
    """
    if kind_a == "O" and kind_b == "O":
        return abs(3 * residue_a - 3 * residue_b) + 2
    if kind_a == "O":
        return abs((3 * residue_a + 2) - (3 * residue_b + _CHAIN_POSITION[kind_b])) + 1
    if kind_b == "O":
        return abs((3 * residue_a + _CHAIN_POSITION[kind_a]) - (3 * residue_b + 2)) + 1
    return abs(
        (3 * residue_a + _CHAIN_POSITION[kind_a]) - (3 * residue_b + _CHAIN_POSITION[kind_b])
    )


#: Bond separation at or below which two backbone atoms are close for covalent reasons and must not
#: be scored as a clash.
#:
#: MEASURED at 4, not the 3 that looks right on paper. A 1-4 pair such as CA(i)-C(i+1) or
#: O(i)-C(i+1) has its distance set mostly by bond geometry, not by packing, so putting it in a
#: 2.9 A clash term adds a penalty no rotation can relieve -- and the optimizer then trades real
#: clashes away to reduce it. Measured on dnmt3a, including 1-4 pairs took steric contacts from 3
#: to 15; excluding them brings it back. This is the same failure as the original 14,859 fake
#: clashes, one shell further out and correspondingly easier to miss.
#:
#: 4 still admits the pair that motivated counting bonds at all: two carbonyl oxygens on adjacent
#: residues are 5 bonds apart.
_BONDED_SEPARATION: Final[int] = 4


def _place_unit(
    ca_a: np.ndarray, ca_b: np.ndarray, azimuth_degrees: float
) -> tuple[np.ndarray, np.ndarray]:
    """Place C(i) and N(i+1) for one peptide unit at the given azimuth."""
    axis, first, second = _azimuth_frame(ca_a, ca_b)
    psi = np.radians(azimuth_degrees)
    radial = np.cos(psi) * first + np.sin(psi) * second
    opposite = -radial

    polar_c = np.radians(_POLAR_C)
    c = ca_a + CA_C_BOND_LENGTH * (np.cos(polar_c) * axis + np.sin(polar_c) * radial)
    polar_n = np.radians(_POLAR_N)
    n = ca_b + N_CA_BOND_LENGTH * (-np.cos(polar_n) * axis + np.sin(polar_n) * opposite)
    return c, n


def _place_oxygen(ca: np.ndarray, c: np.ndarray, n_next: np.ndarray) -> np.ndarray:
    """O is fully determined by CA(i), C(i) and N(i+1) -- exact to 0.013 A from the true three."""
    to_ca = _unit(ca - c)
    to_n = _unit(n_next - c)
    normal = _cross(to_ca, to_n)
    if float(np.linalg.norm(normal)) < 1e-9:
        normal = np.array([0.0, 0.0, 1.0])
    normal = _unit(normal)
    angle = np.radians(CA_C_O_ANGLE)
    cos, sin = np.cos(-angle), np.sin(-angle)
    direction = (
        to_ca * cos
        + _cross(normal, to_ca) * sin
        + normal * float(np.dot(normal, to_ca)) * (1 - cos)
    )
    placed: np.ndarray = c + C_O_BOND_LENGTH * direction
    return placed


def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    v1, v2 = _unit(a - b), _unit(c - b)
    return float(np.degrees(np.arccos(np.clip(float(np.dot(v1, v2)), -1.0, 1.0))))


def _dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b1 = _unit(p2 - p1)
    v = (p0 - p1) - float(np.dot(p0 - p1, b1)) * b1
    w = (p3 - p2) - float(np.dot(p3 - p2, b1)) * b1
    return float(np.degrees(np.arctan2(float(np.dot(_cross(b1, v), w)), float(np.dot(v, w)))))


def _rama_penalty(phi: np.ndarray, psi: np.ndarray) -> np.ndarray:
    """Look up the measured penalty for phi/psi, vectorized over an array of candidates."""
    n_bins = len(_RAMA_PENALTY)
    i = np.clip(((np.asarray(phi) + 180.0) // _RAMA_BIN_WIDTH).astype(np.int64), 0, n_bins - 1)
    j = np.clip(((np.asarray(psi) + 180.0) // _RAMA_BIN_WIDTH).astype(np.int64), 0, n_bins - 1)
    looked_up: np.ndarray = _RAMA_TABLE[i, j]
    return looked_up


#: The penalty table as an array, built once. Kept alongside the tuple form so the source stays
#: readable while lookups stay vectorized.
_RAMA_TABLE: Final[np.ndarray] = np.asarray(_RAMA_PENALTY, dtype=np.float64)


def _angles_batch(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Angle a-b-c in degrees, broadcasting over leading axes."""
    v1 = _unit(a - b)
    v2 = _unit(c - b)
    degrees: np.ndarray = np.degrees(np.arccos(np.clip((v1 * v2).sum(-1), -1.0, 1.0)))
    return degrees


def _dihedrals_batch(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> np.ndarray:
    """Dihedral p0-p1-p2-p3 in degrees, broadcasting over leading axes."""
    b1 = _unit(p2 - p1)
    d0 = p0 - p1
    d3 = p3 - p2
    v = d0 - (d0 * b1).sum(-1, keepdims=True) * b1
    w = d3 - (d3 * b1).sum(-1, keepdims=True) * b1
    degrees: np.ndarray = np.degrees(np.arctan2((_cross(b1, v) * w).sum(-1), (v * w).sum(-1)))
    return degrees


def _place_units_batch(
    ca_a: np.ndarray, ca_b: np.ndarray, azimuths_degrees: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Place C(i) and N(i+1) for one peptide unit at MANY azimuths at once."""
    axis, first, second = _azimuth_frame(ca_a, ca_b)
    psi = np.radians(np.asarray(azimuths_degrees, dtype=np.float64))[:, None]
    radial = np.cos(psi) * first + np.sin(psi) * second
    polar_c = np.radians(_POLAR_C)
    polar_n = np.radians(_POLAR_N)
    c = ca_a + CA_C_BOND_LENGTH * (np.cos(polar_c) * axis + np.sin(polar_c) * radial)
    n = ca_b + N_CA_BOND_LENGTH * (-np.cos(polar_n) * axis - np.sin(polar_n) * radial)
    return c, n


def _place_oxygen_batch(ca: np.ndarray, c: np.ndarray, n_next: np.ndarray) -> np.ndarray:
    """Vectorized :func:`_place_oxygen`."""
    to_ca = _unit(ca - c)
    to_n = _unit(n_next - c)
    normal = _cross(to_ca, to_n)
    norm = np.linalg.norm(normal, axis=-1, keepdims=True)
    safe = np.divide(normal, np.maximum(norm, 1e-12))
    normal = np.where(norm > 1e-9, safe, np.array([0.0, 0.0, 1.0]))
    angle = np.radians(CA_C_O_ANGLE)
    cos, sin = np.cos(-angle), np.sin(-angle)
    dotted = (normal * to_ca).sum(-1, keepdims=True)
    direction = to_ca * cos + _cross(normal, to_ca) * sin + normal * dotted * (1.0 - cos)
    placed: np.ndarray = c + C_O_BOND_LENGTH * direction
    return placed


def _load_kernel() -> Any:
    """Return the compiled kernel module, or None when numba is not importable.

    Cached on the function, so a missing numba costs one failed import per process rather than
    one per region. Deliberately tolerant: numba is a base dependency, but a platform without a
    wheel, or a numpy it has not caught up with, should fall back rather than fail a rebuild.
    """
    cached = getattr(_load_kernel, "_cached", "unset")
    if cached == "unset":
        try:
            from . import backbone_kernel
        except Exception:
            cached = None
        else:
            cached = backbone_kernel
        _load_kernel._cached = cached  # type: ignore[attr-defined]
    return cached


def refine_backbone(
    ca_xyz: np.ndarray,
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
    backend: str = "auto",
) -> RefinementResult:
    """Refine N, C and O of a backbone, holding the alpha carbons fixed.

    Parameters
    ----------
    ca_xyz
        ``(n_residues, 3)`` alpha carbons. Not modified.
    n_xyz, c_xyz
        Starting backbone, as produced by
        :func:`~dodo.construct.ca_backbone.backbone_from_ca`. Used only to seed each azimuth.
    obstacles
        ``(n_obstacles, 3)`` coordinates of atoms DODO did not build -- folded-domain side chains
        especially -- which the refined backbone must avoid.
    clash_distance
        Minimum separation between a refined backbone atom and any other atom, in Angstroms.
    angle_weight, clash_weight, rama_weight
        Relative weights of the N-CA-C term, the clash term and the Ramachandran term. Clashes
        dominate by design: a steric overlap is a hard error, a strained angle a soft one.

        These defaults were TUNED against truth rather than chosen, and the balance matters more
        than it looks. ``angle_weight`` is 1/2.84**2, which expresses the angle deviation in units
        of its real measured spread. Getting it wrong in either direction makes the refinement
        worse than doing nothing:

        ===================  =====  =====  =====  =========
        weights              N      C      O      angle sd
        ===================  =====  =====  =====  =========
        unrefined            0.190  0.291  0.838  11.54
        angle=1.0 rama=1     0.250  0.349  1.075   1.32
        angle=0.124 rama=1   0.205  0.272  0.821   1.84
        **0.124 / 20**       0.170  0.217  0.634   3.37
        angle=0.05 rama=50   0.174  0.225  0.658   5.52
        ===================  =====  =====  =====  =========

        Errors in Angstroms against all-atom simulation, over 1,640 residues; the real N-CA-C
        spread in that data is 2.93 degrees.

        At ``angle=1.0`` the angle term dominates and drives the N-CA-C spread to 1.32 degrees --
        *tighter than the real 2.93* -- by moving every atom further from where it belongs. An
        over-tightened geometric term is not a better structure.
    max_sweeps, tolerance
        Coordinate-descent limits. A sweep in which no azimuth moves more than ``tolerance``
        degrees ends the refinement.
    candidates
        Azimuths tried per unit per sweep. Deterministic grid, refined around the incumbent.

    Returns
    -------
    RefinementResult
        Refined coordinates and what happened.

    Raises
    ------
    ~dodo.exceptions.GeometryError
        If the inputs disagree in length or contain a non-finite coordinate.
    """
    ca = np.asarray(ca_xyz, dtype=np.float64)
    n_start = np.asarray(n_xyz, dtype=np.float64)
    c_start = np.asarray(c_xyz, dtype=np.float64)
    if ca.ndim != 2 or ca.shape[1] != 3:
        raise GeometryError(f"Expected alpha carbons shaped (n_residues, 3), got {ca.shape}.")
    if n_start.shape != ca.shape or c_start.shape != ca.shape:
        raise GeometryError(
            f"N, C and CA arrays must have matching shapes; got {n_start.shape}, "
            f"{c_start.shape} and {ca.shape}."
        )
    if ca.shape[0] < 3:
        raise GeometryError(
            f"Refinement needs at least three residues to have an N-CA-C angle to refine, got "
            f"{ca.shape[0]}."
        )
    for name, array in (("ca_xyz", ca), ("n_xyz", n_start), ("c_xyz", c_start)):
        if not np.all(np.isfinite(array)):
            raise GeometryError(f"{name} contains a non-finite value.")

    n_res = ca.shape[0]
    n_units = n_res - 1
    obstacle_points = None if obstacles is None else np.asarray(obstacles, dtype=np.float64)

    if backend not in ("auto", "numba", "numpy"):
        raise InvalidParameterError(f"Unknown backend {backend!r}. Use 'auto', 'numba' or 'numpy'.")
    kernel = None if backend == "numpy" else _load_kernel()
    if kernel is None and backend == "numba":
        raise MissingDependencyError(
            "numba",
            "compiling the backbone refinement inner loop, which is 16x faster than the pure-numpy "
            "path; pass backend='numpy' to use that instead",
            extra="fast",
        )

    # One object per region, holding the neighbour tables the compiled sweep AND the compiled energy
    # both read. Built here rather than inside either of them so the fixed tables -- a KD-tree over
    # every obstacle in the structure -- are built once for the region instead of once per use.
    #
    # Not built at all with no sweeps to run: there is then nothing for the kernel to do, so the
    # numpy path is the only path, for the energy as well. That keeps `max_sweeps=0` bit-comparable
    # between backends, which is what `test_zero_sweeps_behaves_the_same_either_way` asks for.
    compiled: Any | None = None
    if kernel is not None and max_sweeps > 0:
        try:
            compiled = kernel.RegionKernel(
                ca,
                obstacles=obstacle_points,
                clash_distance=clash_distance,
                angle_weight=angle_weight,
                clash_weight=clash_weight,
                rama_weight=rama_weight,
                rama_table=_RAMA_TABLE,
            )
        except GeometryError:
            # The kernel's neighbour tables hold MAX_NEIGHBOURS rows per movable atom, and the
            # fixed-atom crowding around a region is set by the INPUT -- folded-domain atoms packed
            # however the structure packs them -- so no fixed cap is provably enough (measured:
            # p300's worst shell held 43, then corpus structure Q9C000 produced 54). A region past
            # the cap is not a failed rebuild; it is a region the compiled path cannot score without
            # silently truncating the objective. The numpy path scores the same objective with no
            # cap, so under ``backend="auto"`` degrade to it for this region. An explicit
            # ``backend="numba"`` keeps the loud error: the caller asked for the kernel by name.
            if backend == "numba":
                raise
            compiled = None
    use_kernel = compiled is not None

    from scipy.spatial import cKDTree

    fixed = ca if obstacle_points is None else np.vstack([ca, obstacle_points])

    # Seed each azimuth from the starting placement, so refinement begins where the table put it.
    azimuths = np.empty(n_units)
    for unit in range(n_units):
        axis, first, second = _azimuth_frame(ca[unit], ca[unit + 1])
        offset = c_start[unit] - ca[unit]
        radial = offset - float(np.dot(offset, axis)) * axis
        azimuths[unit] = np.degrees(
            np.arctan2(float(np.dot(radial, second)), float(np.dot(radial, first)))
        )

    # Live state, updated in place. The objective is LOCAL -- changing one azimuth moves exactly
    # three atoms, C(i), N(i+1) and O(i) -- so a candidate is scored by recomputing only the terms
    # that touch them. Scoring the whole chain per candidate is O(sweeps x units x candidates x N)
    # and was measured at over ten minutes for ten chains; this is O(sweeps x units x candidates).
    # ONE buffer, with n/c/o as views into it rather than three separate arrays.
    #
    # This is the difference between a refinement that finishes and one that looks hung. The clash
    # term needs every movable atom as a single (3N, 3) block to index into, and it used to build
    # that block with np.vstack INSIDE the scoring function -- which runs three times per unit per
    # sweep. On p300 that is roughly 137,000 rebuilds of a 1,500-atom array for one structure, and
    # it was 21.5 of the 21.6 seconds `--backbone` added. Writing through views means the block is
    # always current and never rebuilt.
    live_points = np.empty((3 * n_res, 3), dtype=np.float64)
    n_live = live_points[:n_res]
    c_live = live_points[n_res : 2 * n_res]
    o_live = live_points[2 * n_res :]
    n_live[:] = n_start
    c_live[:] = c_start
    o_live[:] = c_start
    for unit in range(n_units):
        c_live[unit], n_live[unit + 1] = _place_unit(ca[unit], ca[unit + 1], azimuths[unit])
    for unit in range(n_units):
        o_live[unit] = _place_oxygen(ca[unit], c_live[unit], n_live[unit + 1])
    # Not _place_oxygen with an invented next N along CA->C: that N is collinear with the bond, so
    # the carbonyl plane is undefined. See :func:`~dodo.construct.ca_backbone._terminal_oxygen`.
    o_live[-1] = _terminal_oxygen(ca[-1], c_live[-1], ca[-2])

    # Per unit, the fixed atoms it could ever clash with. Precomputed ONCE: the movable atoms
    # travel on small circles (C 0.56 A radius, N 0.38, O further), so a generous local shell
    # covers every candidate azimuth and the inner loop needs no tree query at all. Atoms that are
    # close for covalent reasons are excluded by bond count; ignoring that entirely counted 14,859
    # "clashes" across 2,559 residues where the true number was 1.
    reach = clash_distance + 3.0
    # The backbone must also avoid ITSELF, which is a separate problem from avoiding the fixed
    # atoms above and was originally missed. The alpha carbons are guaranteed 3.2 A apart by the
    # engine that placed them, so a CA trace cannot self-intersect -- but the atoms hung off it
    # can, and a carbonyl oxygen points away from the chain by design. Measured on arf19, scoring
    # against fixed atoms alone left 15 steric clashes where the CA-only output had none, every
    # one of them between two generated backbone atoms 6 to 45 residues apart, the worst an N and
    # an O 1.362 A apart. Nothing in the objective was looking at them.
    #
    # These partners move as refinement proceeds, so unlike the fixed ones they cannot be
    # precomputed once. They are rebuilt once per sweep instead: one tree over 3N points plus
    # n_units queries, which is negligible beside the candidate scoring it guards.
    # The three atoms one azimuth moves, and which residue each belongs to. Every exclusion below
    # is computed per MOVABLE ATOM rather than per unit, because C(unit), N(unit+1) and O(unit) sit
    # at different points in the bond graph and so have genuinely different neighbours.
    movable_kinds: tuple[tuple[str, int], ...] = (("C", 0), ("N", 1), ("O", 0))
    movable_residue_of = np.tile(np.arange(n_res), 3)
    movable_kind_of = np.repeat(np.asarray(["N", "C", "O"]), n_res)

    # Fixed-atom neighbours, per unit and per movable atom. Obstacle points carry no residue, so
    # they are never excluded; the alpha carbons, which occupy the first n_res rows, are.
    #
    # Skipped entirely when the kernel is driving, because then nothing reads them: the compiled
    # sweep and the compiled energy both use the equivalent padded tables inside `compiled`.
    # Building them anyway meant a KD-tree over 12,517 obstacles plus 3 x n_units Python
    # `_bond_separation` comprehensions per region, thrown away.
    fixed_neighbours: list[list[np.ndarray]] = []
    if not use_kernel:
        fixed_tree = cKDTree(fixed)
        for unit in range(n_units):
            midpoint = 0.5 * (ca[unit] + ca[unit + 1])
            found = np.asarray(fixed_tree.query_ball_point(midpoint, reach), dtype=np.int64)
            per_atom: list[np.ndarray] = []
            for kind, offset in movable_kinds:
                is_ca = found < n_res
                separation = np.array(
                    [
                        _bond_separation(kind, unit + offset, "CA", int(j))
                        if j < n_res
                        else _BONDED_SEPARATION + 1
                        for j in found
                    ],
                    dtype=np.int64,
                )
                per_atom.append(found[~(is_ca & (separation <= _BONDED_SEPARATION))])
            fixed_neighbours.append(per_atom)

    # The backbone must also avoid ITSELF, which is a separate problem from avoiding the fixed
    # atoms above and was originally missed. The alpha carbons are guaranteed 3.2 A apart by the
    # engine that placed them, so a CA trace cannot self-intersect -- but the atoms hung off it
    # can, and a carbonyl oxygen points away from the chain by design. Measured on arf19, scoring
    # against fixed atoms alone left 15 steric clashes where the CA-only output had none, every
    # one of them between two generated backbone atoms 6 to 45 residues apart, the worst an N and
    # an O 1.362 A apart. Nothing in the objective was looking at them.
    #
    # These partners move as refinement proceeds, so unlike the fixed ones they cannot be
    # precomputed once. They are rebuilt once per sweep instead: one tree over 3N points plus
    # n_units queries, which is negligible beside the candidate scoring it guards.
    movable_neighbours: list[list[np.ndarray]] = [
        [np.empty(0, dtype=np.int64)] * len(movable_kinds) for _ in range(n_units)
    ]

    def refresh_movable_neighbours() -> None:
        tree = cKDTree(live_points)
        for unit in range(n_units):
            midpoint = 0.5 * (ca[unit] + ca[unit + 1])
            found = np.asarray(tree.query_ball_point(midpoint, reach), dtype=np.int64)
            per_atom: list[np.ndarray] = []
            for kind, offset in movable_kinds:
                separation = np.array(
                    [
                        _bond_separation(
                            kind,
                            unit + offset,
                            str(movable_kind_of[j]),
                            int(movable_residue_of[j]),
                        )
                        for j in found
                    ],
                    dtype=np.int64,
                )
                per_atom.append(found[separation > _BONDED_SEPARATION])
            movable_neighbours[unit] = per_atom

    def clash_batch(points: np.ndarray, unit: int, which: int) -> np.ndarray:
        """Clash penalty for one movable atom at MANY candidate positions.

        ``which`` indexes :data:`movable_kinds`: 0 for this unit's C, 1 for its N, 2 for its O.
        """
        penalty = np.zeros(points.shape[0])
        indices = fixed_neighbours[unit][which]
        if indices.size:
            penalty += _clash_penalty(points, fixed[indices], clash_distance, clash_weight)
        moving = movable_neighbours[unit][which]
        if moving.size:
            penalty += _clash_penalty(points, live_points[moving], clash_distance, clash_weight)
        return penalty

    def score_candidates(unit: int, azimuths_deg: np.ndarray) -> np.ndarray:
        """Everything a single azimuth can change, for every candidate at once."""
        c_new, n_new = _place_units_batch(ca[unit], ca[unit + 1], azimuths_deg)
        o_new = _place_oxygen_batch(ca[unit], c_new, n_new)
        total = np.zeros(azimuths_deg.shape[0])

        # N-CA-C at residue `unit` uses N from the previous unit, which is fixed here.
        if 0 < unit < n_res - 1:
            observed = _angles_batch(n_live[unit][None, :], ca[unit][None, :], c_new)
            total += angle_weight * np.square(observed - N_CA_C_ANGLE)
            outside = (observed < N_CA_C_WINDOW_MIN) | (observed > N_CA_C_WINDOW_MAX)
            total += np.where(outside, _ANGLE_WINDOW_PENALTY, 0.0)
            phi = _dihedrals_batch(
                c_live[unit - 1][None, :], n_live[unit][None, :], ca[unit][None, :], c_new
            )
            psi = _dihedrals_batch(n_live[unit][None, :], ca[unit][None, :], c_new, n_new)
            total += rama_weight * np.square(_rama_penalty(phi, psi))

        # N-CA-C at residue `unit + 1` uses C from the next unit, also fixed here.
        if 0 < unit + 1 < n_res - 1:
            observed = _angles_batch(n_new, ca[unit + 1][None, :], c_live[unit + 1][None, :])
            total += angle_weight * np.square(observed - N_CA_C_ANGLE)
            outside = (observed < N_CA_C_WINDOW_MIN) | (observed > N_CA_C_WINDOW_MAX)
            total += np.where(outside, _ANGLE_WINDOW_PENALTY, 0.0)
            phi = _dihedrals_batch(c_new, n_new, ca[unit + 1][None, :], c_live[unit + 1][None, :])
            psi = _dihedrals_batch(
                n_new,
                ca[unit + 1][None, :],
                c_live[unit + 1][None, :],
                n_live[unit + 2][None, :] if unit + 2 < n_res else c_live[unit + 1][None, :],
            )
            total += rama_weight * np.square(_rama_penalty(phi, psi))

        total += clash_batch(c_new, unit, 0)
        total += clash_batch(n_new, unit, 1)
        total += clash_batch(o_new, unit, 2)
        return total

    def apply(unit: int, azimuth: float) -> None:
        c_live[unit], n_live[unit + 1] = _place_unit(ca[unit], ca[unit + 1], azimuth)
        o_live[unit] = _place_oxygen(ca[unit], c_live[unit], n_live[unit + 1])
        if unit + 1 < n_units:
            o_live[unit + 1] = _place_oxygen(ca[unit + 1], c_live[unit + 1], n_live[unit + 2])

    def numpy_total_energy() -> float:
        value = 0.0
        for unit in range(n_units):
            value += float(score_candidates(unit, np.asarray([azimuths[unit]]))[0])
        # Each unit's score covers the residues on both of its sides, so summing over units
        # double-counts the shared terms. Halving keeps the reported number comparable between
        # the before and after readings, which is all it is used for.
        return value / 2.0

    def total_energy() -> float:
        """Return the objective over the whole region, at the live state.

        Called exactly twice -- once for `before`, once for `after` -- and it was measured costing
        2.8x the compiled sweep it reports on, because the numpy form pushes every peptide unit
        through `score_candidates` one array call per term. The compiled form is the same objective,
        established term by term on byte-identical inputs rather than assumed; see
        ``tests/unit/test_backbone_kernel.py``.

        The numpy backend keeps the numpy scorer. Both readings come from the same one, whichever
        that is, so `before` and `after` stay comparable with each other -- the only thing they are
        used for.
        """
        if compiled is not None:
            return float(compiled.energy(n_live, c_live, o_live, azimuths))
        refresh_movable_neighbours()
        return numpy_total_energy()

    before = total_energy()
    sweeps = 0
    converged = False
    previous: float | None = None
    if compiled is not None:
        # ONLY the sweep loop and the objective it minimises are compiled. Everything around them --
        # the seeding, the terminal atoms, the notes, the non-ideal-spacing check, and the fact that
        # `before` and `after` are both real numbers read off the same scorer -- stays shared with
        # the numpy path, so a RefinementResult means the same thing whichever backend produced it.
        # Returning early from the kernel instead left energy_before/after as NaN and dropped the
        # non-ideal-spacing note, which three existing tests caught.
        #
        # The live views and `azimuths` are refined IN PLACE, which is also what removes the second
        # seeding: the kernel used to re-derive the same azimuths and the same canonicalized atoms
        # from `n_start`/`c_start` and hand them back to be copied over these arrays.
        sweeps, converged = compiled.sweep_to_convergence(
            n_live,
            c_live,
            o_live,
            azimuths,
            max_sweeps=max_sweeps,
            tolerance=tolerance,
            candidates=candidates,
        )
    else:
        for sweep in range(max_sweeps):
            sweeps = sweep + 1
            # The generated atoms moved last sweep, so who each unit can hit has changed with them.
            refresh_movable_neighbours()
            largest_move = 0.0
            swept = 0.0
            span = 180.0 / (1.0 + sweep)
            for unit in range(n_units):
                incumbent = azimuths[unit]
                trials = incumbent + np.linspace(-span, span, candidates)
                trials = np.concatenate([[incumbent], trials])
                values = score_candidates(unit, trials)
                best = int(np.argmin(values))
                chosen = float(trials[best])
                swept += float(values[best])
                apply(unit, chosen)
                azimuths[unit] = chosen
                largest_move = max(largest_move, abs(chosen - incumbent))
            # Converge on the OBJECTIVE, not on how far the azimuths moved.
            #
            # The move-size test this replaced could not fire. Candidates are drawn from a window
            # that
            # narrows as 180/(1+sweep) across `candidates` points, so the smallest non-zero move
            # available is the candidate spacing -- 15.65 degrees on the first sweep and still 0.52
            # on
            # the thirtieth, against a 0.25 degree tolerance. `largest_move <= tolerance` therefore
            # meant `largest_move == 0`, so every long region ran the full sweep budget and reported
            # `converged` as False even once the geometry had stopped changing in any way that
            # mattered.
            #
            # `swept` is the sum of the accepted scores over the sweep. It double-counts terms
            # shared
            # between neighbouring units, exactly as total_energy does, which is why it is only ever
            # compared against itself.
            improvement = float("inf") if previous is None else previous - swept
            previous = swept
            if largest_move == 0.0 or improvement <= tolerance:
                converged = True
                break

    after = total_energy()
    # Three atoms belong to no peptide unit -- N of the first residue, C and O of the last -- so no
    # azimuth ever moved them, and they still sit where they were placed relative to neighbours
    # that have since been refined. Re-derive them, or they keep drifting out of geometry as
    # everything around them improves: measured on dnmt3a, N(0) and C(0) ended up 1.788 A apart
    # against a correct 2.459, because N(0) was aimed at where C(0) used to be.
    # Copies, not the views above: the caller must not hold a window onto a shared buffer.
    n_final, c_final, o_final = n_live.copy(), c_live.copy(), o_live.copy()
    if n_res >= 2:
        n_final[0] = _terminal_nitrogen(ca[0], c_final[0], ca[1])
        c_final[-1] = _terminal_carbon(ca[-1], n_final[-1], ca[-2])
        o_final[-1] = _terminal_oxygen(ca[-1], c_final[-1], ca[-2])
    angles = [_angle(n_final[r], ca[r], c_final[r]) for r in range(1, n_res - 1)]
    notes = [
        f"Refined {n_units} peptide unit(s) over {sweeps} sweep(s); objective "
        f"{before:.1f} -> {after:.1f}. N-CA-C angle now {np.mean(angles):.2f} "
        f"+/- {np.std(angles):.2f} deg against an ideal {N_CA_C_ANGLE}."
    ]
    # A peptide unit determines its own CA-CA distance, so input spacing that disagrees with the
    # canonical value lands in the peptide bond. Report it rather than returning stretched bonds
    # with no indication of why -- see the module docstring.
    peptide = np.linalg.norm(n_final[1:] - c_final[:-1], axis=1)
    drift = float(np.max(np.abs(peptide - C_N_PEPTIDE_BOND_LENGTH)))
    if drift > 0.01:
        spacing = np.linalg.norm(np.diff(ca, axis=0), axis=1)
        notes.append(
            f"Input alpha carbons span {spacing.min():.3f}-{spacing.max():.3f} A rather than a "
            f"uniform {CA_CA_BOND_LENGTH} A, so the peptide C-N bond absorbed the difference and "
            f"now spans {peptide.min():.3f}-{peptide.max():.3f} A against an ideal "
            f"{C_N_PEPTIDE_BOND_LENGTH}. Regularize the trace to fix this."
        )
    return RefinementResult(
        n_xyz=n_final,
        c_xyz=c_final,
        o_xyz=o_final,
        azimuths=azimuths,
        sweeps=sweeps,
        converged=converged,
        energy_before=before,
        energy_after=after,
        notes=tuple(notes),
    )
