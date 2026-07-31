"""Measurements of a CA trace, and the acceptance gate every builder must pass.

:func:`validate_ca_trace` is the contract. A generated backbone is acceptable when its
CA-CA bond lengths are within tolerance of :data:`dodo.constants.CA_CA_BOND_LENGTH` and
its CA-CA-CA pseudo-angles lie inside the generation window. "It ran" is not the
criterion, because the code being replaced ran fine: it had no angle constraint at all
and produced measured pseudo-angles from 47 to 178 degrees. No real protein trace
exhibits a 47 degree CA-CA-CA angle -- at 3.8 A bonds it puts CA(i-1) and CA(i+1) 3.1 A
apart, closer than two non-bonded carbons can sit -- and such a trace cannot be
reconstructed to all-atom because there is no backbone dihedral pair that produces it.

The report names every violation with its residue index rather than returning a bare
bool, for the same reason the rest of this package raises instead of returning sentinel
values: a builder that is told *which* residue and *by how much* can retry that step,
and a test failure that says "bond 47 is 5.91 A, expected 3.70-3.90" is diagnosable
where "assert report.ok" is not.

Nothing here is stochastic, so nothing here takes a generator.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.spatial import cKDTree

from ..constants import (
    BACKBONE_ANGLE_MAX,
    BACKBONE_ANGLE_MIN,
    CA_CA_BOND_LENGTH,
    CA_CA_BOND_TOLERANCE,
    CA_CLASH_DISTANCE,
    CLASH_EXCLUDE_WITHIN_RESIDUES,
)
from ..exceptions import GeometryError

__all__ = [
    "TraceReport",
    "TraceViolation",
    "ca_bond_lengths",
    "ca_pseudo_angles",
    "end_to_end",
    "radius_of_gyration",
    "validate_ca_trace",
]

#: Kinds of violation :func:`validate_ca_trace` can report.
ViolationKind = Literal[
    "non_finite",
    "bond_length",
    "pseudo_angle",
    "undefined_angle",
    "steric_clash",
]

#: Separation below which two consecutive CA positions count as coincident, in Angstroms.
#:
#: At or below this the pseudo-angles at both adjacent vertices are undefined, because the
#: vectors whose angle they measure have no direction. That is reported explicitly rather
#: than allowed to surface as NaN.
COINCIDENT_TOLERANCE = 1e-9

#: Slack applied to the pseudo-angle window before flagging a violation, in degrees.
#:
#: The generation grid contains the window endpoints exactly (91 and 150 degrees), and a
#: candidate generated at exactly 91 degrees measures back as 91 +/- ~1e-13 after a
#: rotation and an arccos. Comparing strictly against the window would therefore report
#: spurious violations on perfectly valid geometry roughly half the time it lands on an
#: endpoint. This slack is far below any physically meaningful angle difference.
ANGLE_EDGE_TOLERANCE = 1e-6


def _as_coords(coords: np.ndarray, name: str, minimum_rows: int) -> np.ndarray:
    """Coerce ``coords`` to an ``(n, 3)`` float64 array with at least ``minimum_rows`` rows."""
    array = np.asarray(coords, dtype=np.float64)
    if array.ndim != 2 or array.shape[1] != 3:
        raise GeometryError(f"{name} must have shape (n, 3), got {np.shape(coords)}.")
    if array.shape[0] < minimum_rows:
        raise GeometryError(
            f"{name} needs at least {minimum_rows} coordinates, got {array.shape[0]}."
        )
    return array


def _require_finite(coords: np.ndarray, name: str) -> None:
    """Raise if any coordinate is NaN or infinite, naming the offending rows."""
    bad = np.flatnonzero(~np.all(np.isfinite(coords), axis=1))
    if bad.size:
        preview = ", ".join(str(int(i)) for i in bad[:5])
        more = f" (and {bad.size - 5} more)" if bad.size > 5 else ""
        raise GeometryError(
            f"{name} contains non-finite coordinates at row(s) {preview}{more}. Measuring "
            f"NaN yields NaN; use validate_ca_trace if you want this reported rather than "
            f"raised."
        )


def ca_bond_lengths(ca_coords: np.ndarray) -> np.ndarray:
    """Distances between consecutive alpha carbons, in Angstroms.

    Parameters
    ----------
    ca_coords
        ``(n, 3)`` CA coordinates in chain order, ``n >= 2``.

    Returns
    -------
    np.ndarray
        ``(n - 1,)`` distances. Entry ``i`` is the CA(i)-CA(i+1) distance.

    Raises
    ------
    GeometryError
        If the input is not ``(n, 3)`` with at least 2 rows, or contains non-finite
        values.
    """
    array = _as_coords(ca_coords, "ca_coords", 2)
    _require_finite(array, "ca_coords")
    lengths: np.ndarray = np.linalg.norm(np.diff(array, axis=0), axis=1)
    return lengths


def _pseudo_angles_with_gaps(ca_coords: np.ndarray) -> np.ndarray:
    """Pseudo-angles in degrees, with NaN where the angle is genuinely undefined.

    Private because NaN is exactly what the public API must not hand back. Only
    :func:`validate_ca_trace` uses it, so that a coincident or non-finite coordinate can
    be *reported* against its residue index instead of raising or silently producing a
    NaN some caller compares against a threshold (which is always False, so a NaN angle
    would pass a naive window check).
    """
    vertex = ca_coords[1:-1]
    back = ca_coords[:-2] - vertex
    forward = ca_coords[2:] - vertex

    back_norm = np.linalg.norm(back, axis=1)
    forward_norm = np.linalg.norm(forward, axis=1)

    with np.errstate(divide="ignore", invalid="ignore"):
        cosine = np.einsum("ij,ij->i", back, forward) / (back_norm * forward_norm)
        # Clipping is not cosmetic. For a perfectly straight trace the dot product of the
        # two unit vectors evaluates to -1 - 1e-16, and arccos of that is NaN. A straight
        # segment is not an exotic input: it is what an extended AlphaFold IDR looks like
        # locally, and it is what DODO's own test fixtures are built from.
        angles = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))

    undefined = (
        (back_norm <= COINCIDENT_TOLERANCE)
        | (forward_norm <= COINCIDENT_TOLERANCE)
        | ~np.isfinite(back_norm)
        | ~np.isfinite(forward_norm)
    )
    angles[undefined] = np.nan
    result: np.ndarray = angles
    return result


def ca_pseudo_angles(ca_coords: np.ndarray) -> np.ndarray:
    """CA(i-1)-CA(i)-CA(i+1) pseudo-angles, in degrees.

    Parameters
    ----------
    ca_coords
        ``(n, 3)`` CA coordinates in chain order, ``n >= 3``.

    Returns
    -------
    np.ndarray
        ``(n - 2,)`` angles in degrees, in ``[0, 180]``. Entry ``i`` is the angle at
        vertex CA(i+1), i.e. the angle subtended at the *middle* residue of the triple
        ``(i, i + 1, i + 2)``.

    Raises
    ------
    GeometryError
        If the input is not ``(n, 3)`` with at least 3 rows, contains non-finite values,
        or contains two coincident consecutive positions -- for which the angle is
        undefined. Never returns NaN.
    """
    array = _as_coords(ca_coords, "ca_coords", 3)
    _require_finite(array, "ca_coords")
    angles = _pseudo_angles_with_gaps(array)
    undefined = np.flatnonzero(np.isnan(angles))
    if undefined.size:
        vertices = ", ".join(str(int(i) + 1) for i in undefined[:5])
        more = f" (and {undefined.size - 5} more)" if undefined.size > 5 else ""
        raise GeometryError(
            f"The pseudo-angle at residue index {vertices}{more} is undefined because a "
            f"neighbouring CA is coincident with it (separation <= "
            f"{COINCIDENT_TOLERANCE:g} A). Returning NaN here is how the pre-rewrite "
            f"sampler leaked garbage into finished structures."
        )
    return angles


def radius_of_gyration(coords: np.ndarray) -> float:
    """Radius of gyration of a coordinate set, in Angstroms.

    Unweighted, i.e. every point counts once. Over a CA trace that is the quantity
    ALBATROSS predicts and therefore the one worth comparing a build against.

    Parameters
    ----------
    coords
        ``(n, 3)`` coordinates, ``n >= 2``.

    Returns
    -------
    float
        ``sqrt(mean(|r - r_mean|^2))``, matching
        :meth:`dodo.structure.Structure.radius_of_gyration`.

    Raises
    ------
    GeometryError
        If the input is not ``(n, 3)`` with at least 2 rows, or contains non-finite
        values.
    """
    array = _as_coords(coords, "coords", 2)
    _require_finite(array, "coords")
    centred = array - array.mean(axis=0)
    return float(np.sqrt((centred**2).sum() / array.shape[0]))


def end_to_end(coords: np.ndarray) -> float:
    """Distance from the first coordinate to the last, in Angstroms.

    Parameters
    ----------
    coords
        ``(n, 3)`` coordinates in chain order, ``n >= 2``.

    Returns
    -------
    float
        ``|coords[-1] - coords[0]|``. This is the quantity a
        :class:`~dodo.construct.dimensions.DimensionTarget` sets, so it is the one to
        compare a finished region against.

    Raises
    ------
    GeometryError
        If the input is not ``(n, 3)`` with at least 2 rows, or contains non-finite
        values.
    """
    array = _as_coords(coords, "coords", 2)
    _require_finite(array, "coords")
    return float(np.linalg.norm(array[-1] - array[0]))


@dataclass(frozen=True, slots=True)
class TraceViolation:
    """One specific defect in a CA trace, with the residue it belongs to.

    Attributes
    ----------
    kind
        Which check failed. ``"bond_length"`` and ``"pseudo_angle"`` are out-of-window
        measurements; ``"non_finite"`` is a NaN or infinite coordinate;
        ``"undefined_angle"`` is a pseudo-angle that could not be measured at all because
        a neighbouring CA is coincident; and ``"steric_clash"`` is a pair of
        non-sequential residues closer than the clash distance.
    residue_index
        The residue this violation is reported against, in the caller's numbering (see
        ``residue_offset`` on :func:`validate_ca_trace`). For a bond it is the *first*
        residue of the pair; for an angle it is the vertex, i.e. the middle residue of
        the triple.
    value
        The measured value: Angstroms for a bond, degrees for an angle, NaN when there
        was nothing measurable.
    low, high
        The acceptable range the value fell outside of, or ``None`` when the violation is
        not a range check.
    message
        A complete, self-contained sentence naming the residue and the numbers.
    """

    kind: ViolationKind
    residue_index: int
    value: float
    low: float | None
    high: float | None
    message: str

    def __str__(self) -> str:
        return self.message


@dataclass(frozen=True, slots=True, eq=False)
class TraceReport:
    """The result of :func:`validate_ca_trace`: measurements plus every violation.

    ``eq=False`` because the dataclass-generated ``__eq__`` would compare numpy arrays
    with ``==`` and raise "truth value of an array is ambiguous" the first time anyone
    compared two reports or put one in a set.

    Attributes
    ----------
    n_residues
        Number of CA positions validated.
    bond_lengths
        ``(n - 1,)`` measured CA-CA distances in Angstroms.
    pseudo_angles
        ``(n - 2,)`` measured CA-CA-CA angles in degrees, with NaN at any vertex where the
        angle is undefined. Empty when fewer than 3 residues were given.
    violations
        Every violation found, ordered by residue index.
    bond_range
        The ``(low, high)`` bond length window used, in Angstroms.
    angle_window
        The ``(low, high)`` pseudo-angle window used, in degrees.
    residue_offset
        Added to every reported residue index, so a region validated in isolation reports
        indices in the whole structure's numbering.
    """

    n_residues: int
    bond_lengths: np.ndarray
    pseudo_angles: np.ndarray
    violations: tuple[TraceViolation, ...]
    bond_range: tuple[float, float]
    angle_window: tuple[float, float]
    residue_offset: int

    @property
    def ok(self) -> bool:
        """True when nothing was violated."""
        return not self.violations

    def __bool__(self) -> bool:
        return self.ok

    def __repr__(self) -> str:
        return (
            f"TraceReport({self.n_residues} residues, "
            f"{'valid' if self.ok else f'{len(self.violations)} violations'})"
        )

    def of_kind(self, kind: ViolationKind) -> tuple[TraceViolation, ...]:
        """All violations of one kind, in residue order."""
        return tuple(v for v in self.violations if v.kind == kind)

    @property
    def violating_residues(self) -> tuple[int, ...]:
        """Residue indices involved in at least one violation, ascending and unique."""
        return tuple(sorted({v.residue_index for v in self.violations}))

    def summary(self) -> str:
        """One-line summary of the measured geometry, valid or not."""
        # Filter rather than use nanmin: an all-NaN slice makes nanmin warn and return
        # NaN, so a summary of a totally broken trace would be both noisy and useless.
        bonds = self.bond_lengths[np.isfinite(self.bond_lengths)]
        angles = self.pseudo_angles[~np.isnan(self.pseudo_angles)]
        parts = [f"{self.n_residues} residues"]
        if bonds.size:
            parts.append(f"bonds {bonds.min():.3f}-{bonds.max():.3f} A")
        if angles.size:
            parts.append(f"angles {angles.min():.1f}-{angles.max():.1f} deg")
        parts.append("valid" if self.ok else f"{len(self.violations)} violations")
        return ", ".join(parts)

    def describe(self, max_violations: int = 10) -> str:
        """Full multi-line description: the summary, then the violations one per line.

        Parameters
        ----------
        max_violations
            Truncate the list after this many, with a count of the remainder. A badly
            broken 500-residue chain produces hundreds of violations and an unreadable
            message; the first few are what identifies the failure.
        """
        lines = [self.summary()]
        for violation in self.violations[:max_violations]:
            lines.append(f"  - {violation.message}")
        remaining = len(self.violations) - max_violations
        if remaining > 0:
            lines.append(f"  - ... and {remaining} more violation(s)")
        return "\n".join(lines)

    def raise_if_invalid(self) -> None:
        """Raise :class:`~dodo.exceptions.GeometryError` if anything was violated.

        The one-liner for a builder that has just produced coordinates and must not
        return them if they are unphysical.
        """
        if not self.ok:
            raise GeometryError(f"Invalid CA trace: {self.describe()}")


def validate_ca_trace(
    ca_coords: np.ndarray,
    bond_tolerance: float | None = None,
    angle_window: tuple[float, float] | None = None,
    *,
    clash_distance: float | None = None,
    exclude_within: int = CLASH_EXCLUDE_WITHIN_RESIDUES,
    residue_offset: int = 0,
) -> TraceReport:
    """Check a CA trace for physically valid bond lengths and pseudo-angles.

    Parameters
    ----------
    ca_coords
        ``(n, 3)`` CA coordinates in chain order, ``n >= 2``.
    bond_tolerance
        Allowed deviation from :data:`dodo.constants.CA_CA_BOND_LENGTH`, in Angstroms.
        Defaults to :data:`dodo.constants.CA_CA_BOND_TOLERANCE`.
    angle_window
        ``(low, high)`` acceptable pseudo-angle range in degrees. Defaults to
        ``(BACKBONE_ANGLE_MIN, BACKBONE_ANGLE_MAX)``, DODO's *generation* window.

        Pass ``(BACKBONE_ANGLE_OBSERVED_MIN, BACKBONE_ANGLE_OBSERVED_MAX)`` when checking
        a trace DODO did not generate: real structures populate the full 75-179 degree
        range, and flagging an experimentally observed 80 degree angle as invalid would
        be measuring the input against the generator's taste rather than against physics.
    residue_offset
        Added to every residue index in the report. Pass ``span.start`` when validating a
        region in isolation so that the violations name residues in the structure's own
        numbering -- mixing local and global indices was the single largest bug source in
        the pre-rewrite code.

    Returns
    -------
    TraceReport
        Measurements and violations. Truthy when the trace is valid.

    Raises
    ------
    GeometryError
        Only if the *input shape* is unusable: not ``(n, 3)``, or fewer than 2 rows. Bad
        geometry -- including NaN coordinates -- is reported, not raised, because this
        function's job is to describe what is wrong. Call
        :meth:`TraceReport.raise_if_invalid` to convert the report into an exception.

    Examples
    --------
    >>> import numpy as np
    >>> from dodo.constants import CA_CA_BOND_LENGTH
    >>> straight = np.zeros((5, 3))
    >>> straight[:, 0] = CA_CA_BOND_LENGTH * np.arange(5)
    >>> report = validate_ca_trace(straight)
    >>> report.ok, [v.kind for v in report.violations][:1]
    (False, ['pseudo_angle'])

    A straight line has ideal bonds and 180 degree angles, which is outside the
    generation window -- extended AlphaFold spaghetti is precisely what DODO rebuilds.
    """
    array = _as_coords(ca_coords, "ca_coords", 2)
    tolerance = CA_CA_BOND_TOLERANCE if bond_tolerance is None else float(bond_tolerance)
    # Finiteness matters as much as sign: a NaN threshold makes every comparison below
    # False, which silently disables that half of the gate instead of failing loudly.
    if not math.isfinite(tolerance):
        raise GeometryError(f"bond_tolerance must be finite, got {tolerance}.")
    if tolerance < 0.0:
        raise GeometryError(f"bond_tolerance must be non-negative, got {tolerance}.")
    clash = CA_CLASH_DISTANCE if clash_distance is None else float(clash_distance)
    if not math.isfinite(clash):
        raise GeometryError(f"clash_distance must be finite, got {clash}.")
    if clash < 0.0:
        raise GeometryError(f"clash_distance must be non-negative, got {clash}.")
    clash_distance = clash
    if exclude_within < 0:
        raise GeometryError(f"exclude_within must be non-negative, got {exclude_within}.")
    low_angle, high_angle = (
        (BACKBONE_ANGLE_MIN, BACKBONE_ANGLE_MAX) if angle_window is None else angle_window
    )
    low_angle, high_angle = float(low_angle), float(high_angle)
    if not (math.isfinite(low_angle) and math.isfinite(high_angle)):
        raise GeometryError(f"angle_window must be finite, got ({low_angle}, {high_angle}).")
    if low_angle > high_angle:
        raise GeometryError(
            f"angle_window low ({low_angle}) exceeds high ({high_angle}); the window is "
            f"(low, high) in degrees."
        )
    low_bond = CA_CA_BOND_LENGTH - tolerance
    high_bond = CA_CA_BOND_LENGTH + tolerance

    violations: list[TraceViolation] = []

    # Non-finite coordinates first: they are the cause of any NaN measurement below, so
    # reporting them separately keeps the bond and angle lists honest.
    for row in np.flatnonzero(~np.all(np.isfinite(array), axis=1)):
        index = int(row) + residue_offset
        violations.append(
            TraceViolation(
                kind="non_finite",
                residue_index=index,
                value=float("nan"),
                low=None,
                high=None,
                message=(
                    f"residue {index}: coordinate is not finite "
                    f"({array[row].tolist()}). A NaN coordinate means an upstream "
                    f"sampler failed and returned its failure instead of raising."
                ),
            )
        )

    bonds: np.ndarray = np.linalg.norm(np.diff(array, axis=0), axis=1)
    angles = (
        _pseudo_angles_with_gaps(array) if array.shape[0] >= 3 else np.zeros(0, dtype=np.float64)
    )

    # NaN comparisons are False, so bonds and angles derived from a non-finite coordinate
    # are skipped here rather than reported twice with a meaningless magnitude.
    for i in np.flatnonzero((bonds < low_bond) | (bonds > high_bond)):
        index = int(i) + residue_offset
        violations.append(
            TraceViolation(
                kind="bond_length",
                residue_index=index,
                value=float(bonds[i]),
                low=low_bond,
                high=high_bond,
                message=(
                    f"residues {index}-{index + 1}: CA-CA bond is {bonds[i]:.3f} A, "
                    f"outside {low_bond:.3f}-{high_bond:.3f} A"
                ),
            )
        )

    for i in np.flatnonzero(np.isnan(angles)):
        index = int(i) + 1 + residue_offset
        violations.append(
            TraceViolation(
                kind="undefined_angle",
                residue_index=index,
                value=float("nan"),
                low=None,
                high=None,
                message=(
                    f"residue {index}: CA-CA-CA pseudo-angle is undefined because a "
                    f"neighbouring CA is coincident or non-finite"
                ),
            )
        )

    out_of_window = (angles < low_angle - ANGLE_EDGE_TOLERANCE) | (
        angles > high_angle + ANGLE_EDGE_TOLERANCE
    )
    for i in np.flatnonzero(out_of_window):
        index = int(i) + 1 + residue_offset
        violations.append(
            TraceViolation(
                kind="pseudo_angle",
                residue_index=index,
                value=float(angles[i]),
                low=low_angle,
                high=high_angle,
                message=(
                    f"residue {index}: CA-CA-CA pseudo-angle is {angles[i]:.2f} deg, "
                    f"outside {low_angle:.1f}-{high_angle:.1f} deg"
                ),
            )
        )

    # Non-bonded contacts. Without this the gate is fiction: bonds and angles are purely
    # LOCAL measurements, so a trace can satisfy both perfectly while folding back through
    # itself. A planar regular hexagon of side 3.80 A traversed twice has every bond at
    # exactly 3.800 A and every pseudo-angle at exactly 120.0 deg, yet residues i and i+6
    # occupy the SAME POINT -- and this function reported it valid until the check below
    # existed. Coincident atoms at 0.00 A are the specific defect the gate is advertised
    # to prevent, so it has to actually look for them.
    if clash_distance > 0.0 and array.shape[0] > exclude_within + 1:
        finite_rows = np.flatnonzero(np.all(np.isfinite(array), axis=1))
        if finite_rows.size > exclude_within + 1:
            tree = cKDTree(array[finite_rows])
            for local_i, local_j in tree.query_pairs(clash_distance):
                row_a = int(finite_rows[local_i])
                row_b = int(finite_rows[local_j])
                if abs(row_a - row_b) <= exclude_within:
                    continue
                separation = float(np.linalg.norm(array[row_a] - array[row_b]))
                low_index, high_index = sorted((row_a + residue_offset, row_b + residue_offset))
                violations.append(
                    TraceViolation(
                        kind="steric_clash",
                        residue_index=low_index,
                        value=separation,
                        low=clash_distance,
                        high=None,
                        message=(
                            f"residues {low_index} and {high_index} are {separation:.3f} A "
                            f"apart, closer than the {clash_distance:.2f} A clash distance"
                            + (
                                " -- they are effectively the same point"
                                if separation < 0.5
                                else ""
                            )
                        ),
                    )
                )

    violations.sort(key=lambda v: (v.residue_index, v.kind))
    return TraceReport(
        n_residues=array.shape[0],
        bond_lengths=bonds,
        pseudo_angles=angles,
        violations=tuple(violations),
        bond_range=(low_bond, high_bond),
        angle_window=(low_angle, high_angle),
        residue_offset=residue_offset,
    )
