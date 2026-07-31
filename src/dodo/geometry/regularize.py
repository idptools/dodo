"""Project a CA trace onto exact bond-length geometry.

Why this exists
---------------
Generative models do not produce exact bond lengths. STARLING is a diffusion model: it
denoises coordinates toward a learned distribution, and nothing in that process enforces a
hard geometric constraint. Its output CA traces are close to correct but not exact --
consecutive CA-CA distances scatter around the true virtual bond length rather than sitting
on it. That is inherent to the method, not a bug in it, and it means the coordinates need a
post-generation correction step before DODO can use them.

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
_DEFAULT_MAX_SWEEPS = 500


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
        Whether every bond reached the requested tolerance.
    sweeps
        Sweeps performed.
    max_bond_error_before, max_bond_error_after
        Largest absolute deviation from the target bond length, in Angstroms.
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

    @property
    def rg_change_fraction(self) -> float:
        """Fractional change in radius of gyration. Near zero means shape was preserved."""
        if self.radius_of_gyration_before == 0.0:
            return 0.0
        return (
            self.radius_of_gyration_after - self.radius_of_gyration_before
        ) / self.radius_of_gyration_before

    def summary(self) -> str:
        """One-line human-readable summary of the correction."""
        status = "converged" if self.converged else "DID NOT CONVERGE"
        return (
            f"{self.ca_coords.shape[0]} residues, {status} in {self.sweeps} sweeps: "
            f"max bond error {self.max_bond_error_before:.4f} -> "
            f"{self.max_bond_error_after:.2e} A, atom RMSD {self.rmsd:.4f} A, "
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
    preserve_endpoints: bool = False,
    max_sweeps: int = _DEFAULT_MAX_SWEEPS,
    raise_on_failure: bool = True,
) -> RegularizationResult:
    """Correct a CA trace so every consecutive CA-CA distance is exactly ``bond_length``.

    Parameters
    ----------
    ca_coords
        ``(n, 3)`` CA coordinates in chain order, ``n >= 2``.
    bond_length
        Target CA-CA distance in Angstroms. Defaults to
        :data:`dodo.constants.CA_CA_BOND_LENGTH`.
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
        Maximum correction sweeps before giving up.
    raise_on_failure
        Raise :class:`~dodo.exceptions.GeometryError` if convergence is not reached. When
        False, the caller gets an unconverged result with ``converged=False`` and is
        responsible for checking it -- which is why it defaults to True.

    Returns
    -------
    RegularizationResult
        The corrected trace plus before/after measurements.

    Raises
    ------
    GeometryError
        If the input shape is wrong, contains non-finite coordinates, has coincident
        consecutive residues (there is no bond direction to correct along), or fails to
        converge while ``raise_on_failure`` is set.

    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.default_rng(0)
    >>> trace = np.cumsum(rng.normal(0, 2.2, size=(40, 3)), axis=0)  # wrong bond lengths
    >>> result = regularize_ca_trace(trace)
    >>> bool(result.converged), float(result.max_bond_error_after) < 1e-6
    (True, True)
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

    # Per-atom weights. A pinned endpoint gets zero weight, so the whole of a bond's
    # correction is absorbed by its mobile partner. Using weights rather than branching keeps
    # one code path for both modes, which matters because the pinned case is the one DODO
    # actually relies on and a separate branch would be the less-tested one.
    mobility = np.ones(n, dtype=np.float64)
    if preserve_endpoints:
        mobility[0] = 0.0
        mobility[-1] = 0.0

    converged = False
    sweeps = 0
    for sweep in range(1, max_sweeps + 1):
        sweeps = sweep
        max_error = 0.0

        # Gauss-Seidel sweep: apply each correction immediately so later bonds in the same
        # pass see the updated positions. That converges markedly faster than accumulating
        # all corrections and applying them at once, because a chain constraint is
        # sequentially coupled.
        for i in range(n - 1):
            delta = coords[i + 1] - coords[i]
            length = float(np.linalg.norm(delta))
            if length == 0.0:
                # Defensive only: the pre-scan above rejects coincident input, and a
                # correction cannot create a new coincidence (it moves atoms apart along the
                # bond axis). Kept so a future change to the sweep cannot divide by zero.
                raise GeometryError(
                    f"Residues {i} and {i + 1} became coincident during correction, which "
                    f"should not be reachable. This is a bug in regularize_ca_trace."
                )

            error = length - target
            max_error = max(max_error, abs(error))

            weight_total = mobility[i] + mobility[i + 1]
            if weight_total == 0.0:
                # Both ends pinned: only possible for a 2-residue trace with
                # preserve_endpoints, which is rejected above. Guard anyway rather than
                # divide by zero.
                continue

            # Move each atom along the bond axis in proportion to its mobility, so the pair's
            # weighted centroid is unchanged and the correction stays local.
            direction = delta / length
            shift = error * direction
            coords[i] += shift * (mobility[i] / weight_total)
            coords[i + 1] -= shift * (mobility[i + 1] / weight_total)

        if max_error <= tolerance:
            converged = True
            break

    after_bonds = ca_bond_lengths(coords)
    max_error_after = float(np.max(np.abs(after_bonds - target)))

    if not converged and raise_on_failure:
        raise GeometryError(
            f"Bond-length regularization did not converge in {max_sweeps} sweeps: largest "
            f"remaining bond error is {max_error_after:.3e} A against a tolerance of "
            f"{tolerance:.1e} A. This usually means the input is far from a valid chain -- "
            f"check that the coordinates really are a CA trace in chain order."
        )

    return RegularizationResult(
        ca_coords=coords,
        converged=converged,
        sweeps=sweeps,
        max_bond_error_before=max_error_before,
        max_bond_error_after=max_error_after,
        rmsd=float(np.sqrt(np.mean(np.sum((coords - original) ** 2, axis=1)))),
        end_to_end_before=re_before,
        end_to_end_after=end_to_end(coords),
        radius_of_gyration_before=rg_before,
        radius_of_gyration_after=radius_of_gyration(coords),
    )


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
