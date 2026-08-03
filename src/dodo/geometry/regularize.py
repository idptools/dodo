"""Project a CA trace onto exact bond-length and in-range pseudo-angle geometry.

Why this exists
---------------
Generative models do not produce exact bond lengths. STARLING is a diffusion model: it
denoises coordinates toward a learned distribution, and nothing in that process enforces a
hard geometric constraint. Its output CA traces are close to correct but not exact --
consecutive CA-CA distances scatter around the true virtual bond length rather than sitting
on it. That is inherent to the method, not a bug in it, and it means the coordinates need a
post-generation correction step before DODO can use them.

Bond lengths are not the only thing that needs correcting, and were not the hard part.
MEASURED on 20 STARLING conformers of a 201-residue sequence: after every bond was projected
onto 3.81 A exactly, CA-CA-CA pseudo-angles still spanned 10.5-179.3 deg, with 2.34% of
vertices outside the 75-179 deg observed range and a median of 5 such vertices per conformer.
A 10.5 deg vertex means the chain doubles back on itself within one step, which no all-atom
reconstruction can realize, so a screen that rejects a conformer for any bad vertex rejected
20 of 20. Bond repair alone therefore made the STARLING engine produce nothing.
:func:`regularize_ca_trace` consequently repairs *both*, simultaneously -- see
"Two constraints on one set of coordinates" below.

The same is true of any learned or interpolated source, so this module is written generically
rather than as a STARLING detail.

What "correction" has to mean here
----------------------------------
Naively rebuilding the chain -- walk from residue 0 and place each CA at exactly the bond
length along the direction to the next -- fixes every bond but destroys the conformation:
errors accumulate along the chain, so the far end drifts and the overall dimensions change.
That is unacceptable, because the whole point of using a generative model is that its
*conformation* is worth having. The end-to-end distance and radius of gyration are the
quantities DODO is trying to get right.

So this is a **constrained projection**, not a rebuild: find the nearest configuration in
which every bond is exact, moving each atom as little as possible. The method is the same
iterative pairwise correction used by SHAKE in molecular dynamics -- for each bond, move both
of its atoms along the bond axis by half the error each, and sweep until converged. Splitting
the correction symmetrically is what keeps the displacement local instead of letting it
accumulate down the chain.

Endpoints
---------
:func:`regularize_ca_trace` can hold the first and last residues fixed. That matters for
DODO specifically: an interior IDR has been positioned to connect two fixed folded-domain
anchors, and letting the termini float during regularization would break the closure that was
just achieved. With endpoints pinned, the correction is absorbed by the interior of the chain.

Note the tradeoff, because it is a real one: pinning both endpoints means the total end-to-end
distance is fixed, so a chain whose bonds are systematically too short cannot be corrected by
lengthening bonds alone -- it has to become less extended somewhere. That is the right answer
(the anchors are ground truth and the generated bonds are not), but it does mean a large
systematic bond error shows up as a change in local geometry rather than in overall span.

Two constraints on one set of coordinates
-----------------------------------------
Bonds and pseudo-angles are constraints on the *same* atoms, so they fight: bending a sharp
vertex is done by moving its neighbours, which are exactly the atoms the bond constraint is
holding at 3.81 A. Converging both at once is the whole problem, and two things were needed to
do it.

**Project on the angle, not on the 1-3 distance.** With both bonds at ``b``, an angle bound is
exactly a distance bound on the 1-3 pair, ``d = 2 b sin(theta / 2)``, which is tempting because
it turns the angle into another bond and lets one SHAKE loop handle everything. It is also
badly conditioned at the top of the window: ``dd/dtheta = b cos(theta / 2)`` goes to zero as
theta approaches 180, so a vertex 1 deg from straight needs a 0.0015 A change in ``d``, and the
step that produces is both far too small and directed along the chain rather than across it. At
exactly 180 deg it is a fixed point -- moving two collinear atoms along the line joining them
leaves them collinear. MEASURED, that showed up as a crawl: on the conformers that failed, the
worst angle violation fell from 5.2 deg to 0.3 deg in 50 sweeps and then needed 1272 more to
close the remaining 0.3 deg, and 3 of 20 conformers never closed it inside a 500-sweep budget.
Projecting on ``theta`` itself removes the bad parameterization: those same conformers converge
in at most 29 sweeps. The step used is the Lagrange step along the constraint gradient, which
is also the minimum-displacement correction to first order -- the property that keeps the
conformation.

**Aim inside the window, test against it.** A vertex driven to exactly the boundary chatters
across it at the 1e-6 deg level forever, because the bond pass keeps nudging it, so the loop
never satisfies "every angle is inside the window". MEASURED with a zero margin: 4 of 20
conformers converged. Aiming a margin inside the requested window fixes that and also buys the
headroom for the final phase, which is bond-only: bonds are polished to the requested tolerance
with the angle pass switched off, so a sweep cannot end having traded an exact bond away.
MEASURED, that polish moves angles by at most 0.050 deg, against a 1.0 deg margin.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..constants import CA_CA_BOND_LENGTH, CA_CA_BOND_TOLERANCE
from ..exceptions import GeometryError
from .metrics import ca_bond_lengths, end_to_end, radius_of_gyration

__all__ = [
    "RegularizationResult",
    "regularize_ca_trace",
]

#: Default sweep limit. Convergence is geometric and in practice takes tens of sweeps for a
#: trace already close to correct, which is the case this module is built for. The cap exists
#: so a pathological input cannot spin forever.
#:
#: MEASURED headroom with an angle window in play, on real STARLING output: the joint phase
#: needed at most 29 sweeps and the bond polish at most 38, both on 201-residue conformers.
_DEFAULT_MAX_SWEEPS = 500

#: Degrees the angle projection aims *inside* the requested window.
#:
#: CHOICE, sized from a measurement. Two jobs. It stops the loop chattering on the boundary --
#: with a zero margin the joint phase converged for 4 of 20 measured conformers instead of 20 of
#: 20 -- and it covers the angle change caused by the bond-only polish that runs afterwards,
#: MEASURED at no more than 0.050 deg. 1.0 deg is 20x that, and cost nothing: sweep counts and
#: atom displacements were identical at 0.5 and 1.0 deg. It also leaves repaired vertices off
#: the exact edge of the observed distribution, which is where they least belong.
_DEFAULT_ANGLE_MARGIN = 1.0

#: Bond tolerance the *joint* phase drives to, in Angstroms, before the bond-only polish.
#:
#: DERIVED from the margin. The joint phase only has to make bonds accurate enough that the
#: angle geometry it is correcting is the real one; the polish then takes them to the caller's
#: tolerance. The two numbers are coupled, because a looser joint tolerance leaves the polish
#: more work and so moves angles further: MEASURED, 1e-2 A here moved angles by up to 0.494 deg
#: during the polish, which a 0.25 deg margin could not absorb, while 1e-3 A moved them by at
#: most 0.050 deg. 1e-3 A is also 1% of :data:`dodo.constants.CA_CA_BOND_TOLERANCE`, so the
#: geometry the angle pass sees is already better than DODO accepts anywhere else.
_JOINT_BOND_TOLERANCE = 1e-3

#: ``|sin(theta)|`` below which a vertex counts as collinear and the ``(u, v)`` plane -- and
#: with it the direction the angle gradient points -- is numerically undefined. Below this the
#: projection falls back to a deterministic perpendicular; see :func:`_angle_sweep`.
_COLLINEAR_SINE = 1e-9


@dataclass(frozen=True, slots=True)
class RegularizationResult:
    """A regularized trace, and what the correction cost.

    The before/after dimensions are reported because they are the thing that must NOT change
    much. A regularizer that achieved exact bonds by quietly reshaping the conformation would
    have defeated its own purpose, and the only way to know is to measure.

    Attributes
    ----------
    ca_coords
        The corrected coordinates, ``(n, 3)``.
    converged
        Whether every constraint that was *requested* reached its tolerance -- bonds always,
        pseudo-angles as well when ``angle_window`` was given. It is deliberately a single
        answer to "is this trace usable", because a caller that checked only the bonds is the
        failure mode the angle repair exists to close.
    sweeps
        Sweeps performed. With an ``angle_window`` this is the *joint* phase, and
        :attr:`polish_sweeps` is the bond-only phase that follows it.
    polish_sweeps
        Bond-only sweeps run after the joint phase. ``0`` when no ``angle_window`` was given,
        because then every sweep is already bond-only.
    max_bond_error_before, max_bond_error_after
        Largest absolute deviation from the target bond length, in Angstroms.
    angle_window
        The window pseudo-angles were held inside, in degrees, or ``None`` if angles were not
        constrained. Recorded so a caller can tell "angles were fine" from "angles were never
        checked" -- two very different things to report.
    max_angle_violation_before, max_angle_violation_after
        Largest number of degrees any vertex sat outside ``angle_window``, or ``0.0`` when no
        vertex did. Always ``0.0`` when ``angle_window`` is ``None``, which is why
        :attr:`angle_window` has to be read alongside it.
    rmsd
        Root-mean-square displacement of the atoms, in Angstroms. How far the correction had
        to move things.
    end_to_end_before, end_to_end_after
        End-to-end distance in Angstroms. With ``preserve_endpoints=True`` these are equal by
        construction.
    radius_of_gyration_before, radius_of_gyration_after
        Radius of gyration in Angstroms. The headline check on conformation preservation.
    """

    ca_coords: np.ndarray
    converged: bool
    sweeps: int
    max_bond_error_before: float
    max_bond_error_after: float
    rmsd: float
    end_to_end_before: float
    end_to_end_after: float
    radius_of_gyration_before: float
    radius_of_gyration_after: float
    polish_sweeps: int = 0
    angle_window: tuple[float, float] | None = None
    max_angle_violation_before: float = 0.0
    max_angle_violation_after: float = 0.0

    @property
    def rg_change_fraction(self) -> float:
        """Fractional change in radius of gyration. Near zero means shape was preserved."""
        if self.radius_of_gyration_before == 0.0:
            return 0.0
        return (
            self.radius_of_gyration_after - self.radius_of_gyration_before
        ) / self.radius_of_gyration_before

    @property
    def end_to_end_change_fraction(self) -> float:
        """Fractional change in end-to-end distance. The other half of the shape check."""
        if self.end_to_end_before == 0.0:
            return 0.0
        return (self.end_to_end_after - self.end_to_end_before) / self.end_to_end_before

    def summary(self) -> str:
        """One-line human-readable summary of the correction."""
        status = "converged" if self.converged else "DID NOT CONVERGE"
        sweeps = f"{self.sweeps} sweeps"
        if self.angle_window is not None:
            sweeps = f"{self.sweeps} joint + {self.polish_sweeps} bond-polish sweeps"
        angles = ""
        if self.angle_window is not None:
            angles = (
                f", worst angle excursion outside "
                f"{self.angle_window[0]:.0f}-{self.angle_window[1]:.0f} deg "
                f"{self.max_angle_violation_before:.2f} -> "
                f"{self.max_angle_violation_after:.2e} deg"
            )
        return (
            f"{self.ca_coords.shape[0]} residues, {status} in {sweeps}: "
            f"max bond error {self.max_bond_error_before:.4f} -> "
            f"{self.max_bond_error_after:.2e} A{angles}, atom RMSD {self.rmsd:.4f} A, "
            f"Rg {self.radius_of_gyration_before:.2f} -> "
            f"{self.radius_of_gyration_after:.2f} A "
            f"({self.rg_change_fraction * 100:+.2f}%), "
            f"Re {self.end_to_end_before:.2f} -> {self.end_to_end_after:.2f} A"
        )


def regularize_ca_trace(
    ca_coords: np.ndarray,
    *,
    bond_length: float | None = None,
    tolerance: float = 1e-6,
    angle_window: tuple[float, float] | None = None,
    angle_margin: float = _DEFAULT_ANGLE_MARGIN,
    preserve_endpoints: bool = False,
    max_sweeps: int = _DEFAULT_MAX_SWEEPS,
    raise_on_failure: bool = True,
) -> RegularizationResult:
    """Correct a CA trace so bonds are exactly ``bond_length`` and angles are in range.

    Parameters
    ----------
    ca_coords
        ``(n, 3)`` CA coordinates in chain order, ``n >= 2``.
    bond_length
        Target CA-CA distance in Angstroms. Defaults to
        :data:`dodo.constants.CA_CA_BOND_LENGTH`.
    angle_window
        ``(minimum, maximum)`` CA-CA-CA pseudo-angle in degrees, held at *every* vertex.
        ``None``, the default, corrects bond lengths only and leaves angles wherever the input
        had them -- which is the right thing for a trace whose angles are already known to be
        in range, and the wrong thing for generative-model output. Pass
        ``(BACKBONE_ANGLE_OBSERVED_MIN, BACKBONE_ANGLE_OBSERVED_MAX)`` for the latter. Ignored
        for a trace of fewer than 3 residues, which has no vertex.
    angle_margin
        Degrees inside ``angle_window`` the correction aims for. Defaults to
        :data:`_DEFAULT_ANGLE_MARGIN`; see it for why a nonzero value is required rather than
        merely nice, and the module docstring for how it interacts with the bond polish.
    tolerance
        Convergence tolerance on each bond, in Angstroms. The default of 1e-6 is far tighter
        than :data:`dodo.constants.CA_CA_BOND_TOLERANCE` on purpose: this function's job is to
        make the bonds *exact*, and a loose tolerance here would just move the problem
        downstream into the writer and the viewer.
    preserve_endpoints
        Hold the first and last residues fixed. Required when the trace has already been
        placed against fixed anchors, since letting the termini move would break that
        placement. See the module docstring for the tradeoff this imposes.
    max_sweeps
        Maximum correction sweeps before giving up, applied to each phase separately.
    raise_on_failure
        Raise :class:`~dodo.exceptions.GeometryError` if convergence is not reached. When
        False, the caller gets an unconverged result with ``converged=False`` and is
        responsible for checking it -- which is why it defaults to True.

    Returns
    -------
    RegularizationResult
        The corrected trace plus before/after measurements.

    Notes
    -----
    The two constraints are not always jointly satisfiable, and this function does not pretend
    otherwise. The clearest case is ``preserve_endpoints`` on a 3-residue trace: with both
    termini pinned and both bonds fixed, the 1-3 distance -- and therefore the one pseudo-angle
    -- is fully determined by the anchors, so if that angle is outside the window nothing can
    move it. Such a trace comes back with ``converged=False`` and a
    :attr:`~RegularizationResult.max_angle_violation_after` saying by how much, and the caller
    is expected to reject it rather than to receive coordinates that have been forced.

    Raises
    ------
    GeometryError
        If the input shape is wrong, contains non-finite coordinates, has coincident
        consecutive residues (there is no bond direction to correct along), has an invalid
        ``angle_window`` or ``angle_margin``, or fails to converge while ``raise_on_failure``
        is set.

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.default_rng(0)
    >>> trace = np.cumsum(rng.normal(0, 2.2, size=(40, 3)), axis=0)  # wrong bond lengths
    >>> result = regularize_ca_trace(trace)
    >>> bool(result.converged), float(result.max_bond_error_after) < 1e-6
    (True, True)

    With an angle window, sharp vertices are opened out as well:

    >>> from dodo.constants import BACKBONE_ANGLE_OBSERVED_MAX, BACKBONE_ANGLE_OBSERVED_MIN
    >>> window = (BACKBONE_ANGLE_OBSERVED_MIN, BACKBONE_ANGLE_OBSERVED_MAX)
    >>> fixed = regularize_ca_trace(trace, angle_window=window)
    >>> bool(fixed.converged), float(fixed.max_angle_violation_after)
    (True, 0.0)
    """
    array = np.asarray(ca_coords, dtype=np.float64)
    if array.ndim != 2 or array.shape[1] != 3:
        raise GeometryError(f"ca_coords must have shape (n, 3), got {array.shape}.")
    if array.shape[0] < 2:
        raise GeometryError(
            f"Regularizing bond lengths needs at least 2 residues, got {array.shape[0]}."
        )
    if not np.all(np.isfinite(array)):
        bad = int(np.flatnonzero(~np.all(np.isfinite(array), axis=1))[0])
        raise GeometryError(
            f"ca_coords contains non-finite values, first at residue {bad}. A NaN coordinate "
            f"means the generator failed and returned its failure instead of raising; "
            f"regularizing it would hide that."
        )

    target = CA_CA_BOND_LENGTH if bond_length is None else float(bond_length)
    if not np.isfinite(target) or target <= 0.0:
        raise GeometryError(f"bond_length must be finite and positive, got {target}.")
    if tolerance <= 0.0 or not np.isfinite(tolerance):
        raise GeometryError(f"tolerance must be finite and positive, got {tolerance}.")
    if max_sweeps < 1:
        raise GeometryError(f"max_sweeps must be at least 1, got {max_sweeps}.")
    if angle_window is not None:
        _validate_angle_window(angle_window, angle_margin)

    n = array.shape[0]
    if preserve_endpoints and n < 3:
        raise GeometryError(
            f"preserve_endpoints needs at least 3 residues so the interior can absorb the "
            f"correction, got {n}."
        )

    original = array.copy()
    coords = array.copy()

    before_bonds = ca_bond_lengths(coords)

    # Screen for coincident consecutive residues BEFORE correcting anything.
    #
    # This has to be a pre-scan rather than a check inside the sweep. A Gauss-Seidel sweep
    # applies each bond's correction immediately, so fixing bond i-1 moves residue i and
    # dissolves a coincidence at bond i before that bond is ever examined -- the projection
    # would recover silently and the caller would never learn its generator emitted two atoms
    # at the same point. Recovering is not the goal; a coincident pair means the upstream
    # sampler failed, and whatever else it got wrong is still in the coordinates.
    coincident = np.flatnonzero(before_bonds == 0.0)
    if coincident.size:
        first = int(coincident[0])
        raise GeometryError(
            f"Residues {first} and {first + 1} are coincident"
            + (f" ({coincident.size} such pairs)" if coincident.size > 1 else "")
            + ". There is no bond direction to correct along, and coincident consecutive "
            "atoms mean the generator failed rather than being something to project away."
        )

    max_error_before = float(np.max(np.abs(before_bonds - target)))
    rg_before = radius_of_gyration(coords) if n >= 2 else 0.0
    re_before = end_to_end(coords)

    # A window is only meaningful where there is a vertex to apply it at. A 2-residue trace has
    # none, so the request is recorded and the angle phase is skipped rather than erroring on
    # something that is trivially satisfied.
    constrain_angles = angle_window is not None and n >= 3
    violation_before = (
        _max_angle_violation(coords, angle_window) if angle_window is not None else 0.0
    )

    # Per-atom weights. A pinned endpoint gets zero weight, so the whole of a bond's
    # correction is absorbed by its mobile partner. Using weights rather than branching keeps
    # one code path for both modes, which matters because the pinned case is the one DODO
    # actually relies on and a separate branch would be the less-tested one.
    mobility = np.ones(n, dtype=np.float64)
    if preserve_endpoints:
        mobility[0] = 0.0
        mobility[-1] = 0.0

    sweeps = 0
    polish_sweeps = 0
    joint_converged = True

    if angle_window is not None and constrain_angles:
        low_aim = angle_window[0] + angle_margin
        high_aim = angle_window[1] - angle_margin
        joint_converged = False
        for sweep in range(1, max_sweeps + 1):
            sweeps = sweep
            _angle_sweep(coords, mobility, low_aim, high_aim)
            _bond_sweep(coords, mobility, target)
            # Measure the state the sweep actually left behind, and measure the angles against
            # the REQUESTED window rather than the inset one it aimed at. A vertex parked on the
            # inset bound crosses it by ~1e-6 deg on every sweep for as long as the bond pass
            # keeps nudging it, so testing against the aim point never terminates -- MEASURED, it
            # cost 13 of 20 conformers their convergence. The margin is there to be spent.
            if (
                _max_bond_error(coords, target) <= _JOINT_BOND_TOLERANCE
                and _max_angle_violation(coords, angle_window) <= 0.0
            ):
                joint_converged = True
                break

    # Bond-only phase. When there is no angle window this is the whole algorithm and behaves
    # exactly as it always did. When there is one it is the polish that takes bonds from
    # _JOINT_BOND_TOLERANCE to the caller's tolerance with nothing pulling the other way, so
    # the function cannot return having traded an exact bond for an angle. It is safe for the
    # angles because the margin covers the movement -- MEASURED at 0.050 deg against 1.0.
    bond_converged = False
    for sweep in range(1, max_sweeps + 1):
        if constrain_angles:
            polish_sweeps = sweep
        else:
            sweeps = sweep
        if _bond_sweep(coords, mobility, target) <= tolerance:
            bond_converged = True
            break

    after_bonds = ca_bond_lengths(coords)
    max_error_after = float(np.max(np.abs(after_bonds - target)))
    violation_after = (
        _max_angle_violation(coords, angle_window) if angle_window is not None else 0.0
    )
    converged = bond_converged and joint_converged and violation_after <= 0.0

    if not converged and raise_on_failure:
        raise GeometryError(
            _failure_message(
                max_sweeps=max_sweeps,
                bond_converged=bond_converged,
                max_error_after=max_error_after,
                tolerance=tolerance,
                angle_window=angle_window if constrain_angles else None,
                violation_after=violation_after,
            )
        )

    return RegularizationResult(
        ca_coords=coords,
        converged=converged,
        sweeps=sweeps,
        polish_sweeps=polish_sweeps,
        angle_window=angle_window,
        max_angle_violation_before=violation_before,
        max_angle_violation_after=violation_after,
        max_bond_error_before=max_error_before,
        max_bond_error_after=max_error_after,
        rmsd=float(np.sqrt(np.mean(np.sum((coords - original) ** 2, axis=1)))),
        end_to_end_before=re_before,
        end_to_end_after=end_to_end(coords),
        radius_of_gyration_before=rg_before,
        radius_of_gyration_after=radius_of_gyration(coords),
    )


# ---------------------------------------------------------------------------
# The two constraint sweeps, and the measurements that decide when to stop
# ---------------------------------------------------------------------------


def _bond_sweep(coords: np.ndarray, mobility: np.ndarray, target: float) -> float:
    """One Gauss-Seidel pass setting every consecutive CA-CA distance to ``target`` in place.

    Gauss-Seidel rather than Jacobi: each correction is applied immediately so later bonds in
    the same pass see the updated positions. That converges markedly faster than accumulating
    all corrections and applying them at once, because a chain constraint is sequentially
    coupled.

    Returns
    -------
    float
        The largest absolute bond error *encountered*, i.e. measured before each bond's own
        correction. That is the quantity the convergence test wants: it is zero only if the
        pass found nothing left to fix.
    """
    n = coords.shape[0]
    max_error = 0.0
    for i in range(n - 1):
        delta = coords[i + 1] - coords[i]
        length = float(np.linalg.norm(delta))
        if length == 0.0:
            # Defensive only. The pre-scan in regularize_ca_trace rejects coincident input, a
            # bond correction cannot create a new coincidence (it moves atoms apart along the
            # bond axis), and the angle sweep moves atoms across the chain rather than along it.
            raise GeometryError(
                f"Residues {i} and {i + 1} became coincident during correction, which "
                f"should not be reachable. This is a bug in regularize_ca_trace."
            )

        error = length - target
        max_error = max(max_error, abs(error))

        weight_total = mobility[i] + mobility[i + 1]
        if weight_total == 0.0:
            # Both ends pinned: only possible for a 2-residue trace with preserve_endpoints,
            # which regularize_ca_trace rejects. Guard anyway rather than divide by zero.
            continue

        # Move each atom along the bond axis in proportion to its mobility, so the pair's
        # weighted centroid is unchanged and the correction stays local.
        direction = delta / length
        shift = error * direction
        coords[i] += shift * (mobility[i] / weight_total)
        coords[i + 1] -= shift * (mobility[i + 1] / weight_total)
    return max_error


def _angle_sweep(
    coords: np.ndarray, mobility: np.ndarray, low_degrees: float, high_degrees: float
) -> None:
    """One Gauss-Seidel pass pulling every out-of-range pseudo-angle into range, in place.

    Only *violated* vertices are touched. The constraint is an inequality, so a vertex already
    inside the window is left exactly where it is -- which is what keeps the correction's cost
    proportional to the damage rather than to the chain length.

    The step is the Lagrange step along the constraint gradient. With ``u = x[i-1] - x[i]`` and
    ``v = x[i+1] - x[i]`` of lengths ``p`` and ``q``, and ``c = cos(theta)``,
    ``s = sin(theta)``, the unit vectors that rotate ``u`` and ``v`` in the ``(u, v)`` plane in
    the direction of increasing ``theta`` are ``(c*u_hat - v_hat)/s`` and
    ``(c*v_hat - u_hat)/s``, so the gradients are those over ``p`` and ``q``, and the vertex's
    own gradient is minus their sum by translation invariance. Scaling the whole displacement
    by one ``lambda`` chosen to close the angle error makes this the minimum-displacement
    correction to first order, which is why it perturbs the conformation as little as it does.

    Why the gradient and not the equivalent 1-3 distance constraint: see the module docstring.
    Short version, MEASURED -- the distance form crawls near 180 deg and left 3 of 20 real
    conformers unrepaired after 500 sweeps, while this form finished all 20 within 29.
    """
    low = np.radians(low_degrees)
    high = np.radians(high_degrees)
    for i in range(1, coords.shape[0] - 1):
        u = coords[i - 1] - coords[i]
        v = coords[i + 1] - coords[i]
        p = float(np.linalg.norm(u))
        q = float(np.linalg.norm(v))
        if p == 0.0 or q == 0.0:
            # No angle is defined at all. The bond sweep is the one that reports this.
            continue
        u_hat = u / p
        v_hat = v / q
        cosine = float(np.clip(u_hat @ v_hat, -1.0, 1.0))
        theta = float(np.arccos(cosine))
        if theta < low:
            aim = low
        elif theta > high:
            aim = high
        else:
            continue

        sine = float(np.sqrt(max(0.0, 1.0 - cosine * cosine)))
        if sine < _COLLINEAR_SINE:
            # The three atoms are collinear, so the plane they would be bent in -- and with it
            # the gradient direction -- is undefined. Any perpendicular is as good as any
            # other, so one is chosen deterministically; nothing here may depend on an RNG.
            # The two unit gradients are opposed in the theta -> 0 limit and aligned in the
            # theta -> 180 limit, and getting that relative sign wrong makes the step a
            # translation of the whole vertex, which changes no angle at all.
            perpendicular = _any_perpendicular(u_hat)
            grad_previous = perpendicular / p
            grad_next = (perpendicular if cosine < 0.0 else -perpendicular) / q
        else:
            grad_previous = (cosine * u_hat - v_hat) / (sine * p)
            grad_next = (cosine * v_hat - u_hat) / (sine * q)
        grad_vertex = -(grad_previous + grad_next)

        denominator = (
            mobility[i - 1] * float(grad_previous @ grad_previous)
            + mobility[i] * float(grad_vertex @ grad_vertex)
            + mobility[i + 1] * float(grad_next @ grad_next)
        )
        if denominator == 0.0:
            # Every atom of this vertex is pinned, so its angle is fixed by the anchors and no
            # correction exists. regularize_ca_trace reports that as non-convergence.
            continue
        multiplier = (aim - theta) / denominator
        coords[i - 1] += multiplier * mobility[i - 1] * grad_previous
        coords[i] += multiplier * mobility[i] * grad_vertex
        coords[i + 1] += multiplier * mobility[i + 1] * grad_next


def _any_perpendicular(axis: np.ndarray) -> np.ndarray:
    """Return a unit vector perpendicular to ``axis``, chosen deterministically.

    Crossing with whichever cardinal axis ``axis`` is least aligned to keeps the cross product
    away from zero, so the result is well conditioned for every input direction.
    """
    helper = np.zeros(3)
    helper[int(np.argmin(np.abs(axis)))] = 1.0
    perpendicular = np.cross(axis, helper)
    return perpendicular / float(np.linalg.norm(perpendicular))


def _max_bond_error(coords: np.ndarray, target: float) -> float:
    """Largest absolute deviation of a consecutive CA-CA distance from ``target``."""
    return float(np.abs(np.linalg.norm(np.diff(coords, axis=0), axis=1) - target).max())


def _max_angle_violation(coords: np.ndarray, window: tuple[float, float]) -> float:
    """Degrees by which the worst vertex sits outside ``window``. ``0.0`` if none does.

    Vertices whose angle is undefined -- a coincident consecutive pair -- are skipped rather
    than counted, and that is safe rather than lenient: a zero-length bond puts the bond error
    at the full target length, so such a trace cannot be reported as converged whatever this
    returns. Measuring it non-raising matters because the projection loop calls this on
    intermediate coordinates once per sweep and must not abort part-way through a correction.
    """
    if coords.shape[0] < 3:
        return 0.0
    u = coords[:-2] - coords[1:-1]
    v = coords[2:] - coords[1:-1]
    norms = np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1)
    defined = norms > 0.0
    if not bool(defined.any()):
        return 0.0
    cosine = np.einsum("ij,ij->i", u[defined], v[defined]) / norms[defined]
    angles = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    excursion = np.maximum(np.maximum(window[0] - angles, angles - window[1]), 0.0)
    return float(excursion.max())


def _validate_angle_window(window: tuple[float, float], margin: float) -> None:
    """Reject an angle window or margin that cannot describe a real constraint."""
    if len(window) != 2:
        raise GeometryError(f"angle_window must be a (minimum, maximum) pair, got {window!r}.")
    low, high = (float(value) for value in window)
    if not (np.isfinite(low) and np.isfinite(high)):
        raise GeometryError(f"angle_window bounds must be finite, got {window!r}.")
    if not (0.0 <= low < high <= 180.0):
        raise GeometryError(
            f"angle_window must satisfy 0 <= minimum < maximum <= 180 degrees, got {window!r}."
        )
    if not np.isfinite(margin) or margin < 0.0:
        raise GeometryError(f"angle_margin must be finite and non-negative, got {margin!r}.")
    if 2.0 * margin >= high - low:
        raise GeometryError(
            f"An angle_margin of {margin} deg leaves nothing of the {low}-{high} deg window to "
            f"aim at. The margin is taken off both ends, so it must be under half the width."
        )


def _failure_message(
    *,
    max_sweeps: int,
    bond_converged: bool,
    max_error_after: float,
    tolerance: float,
    angle_window: tuple[float, float] | None,
    violation_after: float,
) -> str:
    """Explain which constraint did not converge, naming both numbers.

    Both constraints are reported when both failed. Naming only the first would hide the more
    interesting case: bonds that look fine because the angle repair gave up, or angles that look
    fine because the bond polish never ran.
    """
    parts: list[str] = []
    if not bond_converged:
        parts.append(
            f"the largest remaining bond error is {max_error_after:.3e} A against a tolerance "
            f"of {tolerance:.1e} A"
        )
    if angle_window is not None and violation_after > 0.0:
        parts.append(
            f"the worst pseudo-angle is still {violation_after:.3f} deg outside the requested "
            f"{angle_window[0]:.0f}-{angle_window[1]:.0f} deg window"
        )
    if not parts:
        # Reachable when the joint phase ran out of sweeps but the final geometry happens to
        # satisfy both constraints. Saying so is more useful than inventing a number.
        parts.append(
            "the joint bond-and-angle phase did not reach its own tolerance, even though the "
            "returned geometry does satisfy both constraints"
        )
    detail = "; and ".join(parts)
    hint = (
        "This usually means the input is far from a valid chain -- check that the coordinates "
        "really are a CA trace in chain order."
    )
    if angle_window is not None:
        hint = (
            "Bond and angle constraints are not always jointly satisfiable -- a vertex whose "
            "neighbours are pinned has no freedom left -- so this can be a property of the "
            "input rather than of the projection. Reject the trace rather than forcing it."
        )
    return f"CA trace regularization did not converge in {max_sweeps} sweeps: {detail}. {hint}"


def needs_regularization(
    ca_coords: np.ndarray,
    *,
    bond_length: float | None = None,
    tolerance: float | None = None,
) -> bool:
    """Return True if any consecutive CA-CA distance is outside tolerance.

    A cheap pre-check, so a trace that is already exact is not put through a pointless
    projection. Defaults to :data:`dodo.constants.CA_CA_BOND_TOLERANCE`, i.e. the tolerance
    DODO accepts elsewhere, rather than the much tighter convergence tolerance
    :func:`regularize_ca_trace` drives to.
    """
    target = CA_CA_BOND_LENGTH if bond_length is None else float(bond_length)
    allowed = CA_CA_BOND_TOLERANCE if tolerance is None else float(tolerance)
    bonds = ca_bond_lengths(np.asarray(ca_coords, dtype=np.float64))
    return bool(np.any(np.abs(bonds - target) > allowed))
