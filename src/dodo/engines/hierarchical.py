"""Long-IDR assembly: build a chain longer than the generator's cap out of blobs.

STARLING will not generate a region longer than :data:`dodo.constants.STARLING_MAX_LENGTH`
(380 residues, observed against STARLING 2.0.2 and re-read from the install at runtime
where it exposes the number). Real IDRs are routinely longer -- the disordered N-terminus
of p300 alone is over a thousand residues -- so "longer than 380" must be a supported case,
not an error and not a silent downgrade to a random walk.

The strategy: blob assembly
---------------------------
Split the region into ``k`` segments no longer than the cap, generate each with the
sub-engine, and arrange the segments in space.

This is correct on polymer-scaling grounds rather than being a workaround. For a polymer
with Flory exponent ``nu``, a segment of ``N / k`` residues has an end-to-end distance

    ``Re_seg ~ (N / k)**nu``

and arranging ``k`` such blobs self-avoidingly, with a step size of order ``Re_seg``, gives
a global end-to-end distance

    ``Re_global ~ Re_seg * k**nu ~ (N / k)**nu * k**nu = N**nu``.

So the assembled chain scales with length exactly as a single long chain would: the blob
decomposition is a renormalization of the same polymer, not an approximation to it. The
sub-chain statistics come from the generator and the super-chain statistics come from the
arrangement, and the two are consistent by construction.

That argument is a reason to believe the *scaling* is right, and it is not enough on its
own -- it fixes an exponent, not a number. So the arrangement does not rely on it: it
targets the requested full-length end-to-end distance directly, taken from
:mod:`dodo.construct.dimensions`, which computed it for the whole sequence. The machinery
for that is the same reachability funnel a single-chain builder uses one level down: at
each junction, the remaining segments can displace the chain end by at most the sum of
their own end-to-end distances, so a junction that puts the target out of reach is
rejected before it is taken. ``k^(1-nu) > 1`` guarantees the segments have more reach than
the target needs, which is why the funnel has room to work.

Junctions are the weak point, and what is done about it
-------------------------------------------------------
Two independently generated segments butt-joined together have no correct local backbone
statistics at the joint: nothing relates the last bond of one to the first bond of the
next, so the CA-CA distance and the two CA-CA-CA angles spanning the joint are whatever
the arbitrary relative orientation happened to give. That is exactly the class of defect
this rewrite exists to remove.

The fix is :data:`dodo.constants.SEGMENT_SPLICE_OVERLAP` residues of overlap. Adjacent
segments are generated over *overlapping* residue ranges, so for a splice at residue ``m``
both segments contain conformations of residues ``m - 2`` and ``m - 1``. The incoming
segment is then placed by a two-point alignment: its copy of residue ``m - 1`` is
translated onto the placed chain's residue ``m - 1``, and its ``m-1 -> m-2`` direction is
rotated onto the placed chain's. After that:

* the bond ``(m-1, m)`` is the incoming segment's own bond, so it is valid;
* the angle at ``m - 1``, spanning ``m-2`` (kept) / ``m-1`` (coincident) / ``m`` (incoming),
  equals the incoming segment's own angle at that vertex, because the ``m-2`` direction was
  aligned -- so it too is valid;
* the angle at ``m - 2`` involves only kept residues.

Every angle and bond at the junction therefore comes from a single continuous
generator-produced trace, which is what "honest local geometry" means here. The alignment
leaves one free degree of freedom -- the spin about the ``m-1 -> m-2`` axis -- and that is
what the arrangement spends on hitting the global target and avoiding clashes. The choice
of splice point *within* the overlap is a second, discrete degree of freedom used the same
way.

Nothing is taken on trust: the finished chain goes through
:func:`dodo.geometry.metrics.validate_ca_trace`, and a bond or angle violation inside a
junction window raises rather than being noted, because the construction above says there
cannot be one. A violation elsewhere belongs to the sub-engine and is reported against it by
name. A non-bonded *contact* is the exception to the attribution rule: the two-point
alignment guarantees the junction's bond and its two angles and says nothing about global
self-avoidance, so a contact is never blamed on the splice.

Where the two kinds of clash are handled
----------------------------------------
Self-avoidance -- the chain against itself -- is enforced *during* assembly, at every
junction, because it is a property of the chain's own shape and does not depend on where the
chain sits in the world.

External obstacles, meaning the folded domains and everything else already placed, are
handled *after* assembly, when the finished chain is moved onto its anchors by
:func:`dodo.engines.starling.place_between_anchors`. Assembly runs in the first segment's own
frame, so an obstacle avoided while splicing would have been avoided in a frame the chain
does not end up occupying -- which is worse than not trying, because the result looks like it
worked. If no orientation clears the obstacles, that conformation is reported as a failure
through the success mask rather than returned.

That covers contacts *between* segments and contacts with the world. It does not cover a
contact that arrives already inside a segment, which no junction-level check can see and
which used to pass through the whole assembly path undetected -- bonds and angles are local
measurements and cannot notice a chain folding back through itself. So the assembled chain is
also measured as a whole, at the clash threshold the result will disclose through
:attr:`~dodo.engines.base.IDRResult.relaxed_to`. A sub-engine that declares a relaxed
threshold has that threshold carried onto the assembled result, because the segment is
spliced in unchanged; a contact tighter than anything declared fails the conformation.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from ..constants import (
    BACKBONE_ANGLE_OBSERVED_MAX,
    BACKBONE_ANGLE_OBSERVED_MIN,
    CA_CA_BOND_LENGTH,
    CA_CLASH_DISTANCE,
    CLASH_EXCLUDE_WITHIN_RESIDUES,
    FLORY_RE_EXPONENT,
    MAX_ATTEMPTS_PER_REGION,
    SEGMENT_SPLICE_OVERLAP,
    contour_length,
)
from ..exceptions import (
    BuildError,
    EngineError,
    ExhaustedAttemptsError,
    GeometryError,
    UnsatisfiableTargetError,
)
from ..geometry.metrics import TraceReport, validate_ca_trace
from ..geometry.transforms import apply, rotation_between_vectors, rotation_from_axis_angle
from .base import ConformationEngine, IDRRequest, IDRResult
from .starling import (
    END_TO_END_ABSOLUTE_TOLERANCE,
    END_TO_END_RELATIVE_TOLERANCE,
    AnchorPlacement,
    LengthCap,
    desired_internal_span,
    end_to_end_tolerance,
    place_between_anchors,
    starling_max_length,
)

__all__ = [
    # Re-exported from dodo.engines.starling so that "the right dimension" has exactly one
    # definition across the two engines that both check it.
    "END_TO_END_ABSOLUTE_TOLERANCE",
    "END_TO_END_RELATIVE_TOLERANCE",
    "AssemblyReport",
    "HierarchicalEngine",
    "Junction",
    "SegmentSpan",
    "default_sub_engine",
    "resolve_segment_cap",
    "segment_spans",
]

#: Spin angles scored per candidate splice point.
#:
#: CHOICE, and cheap for a specific reason: the position of the chain's far end is a
#: sinusoid in the spin angle (Rodrigues' formula has only ``cos`` and ``sin`` terms), so
#: the squared distance to any fixed point is ``a + b cos(phi) + c sin(phi)``. A single
#: period therefore has one minimum and one maximum, and 64 samples locate the optimum to
#: within ~3 degrees -- there is no structure between samples to miss.
SPIN_SAMPLES = 64

#: Best-scoring spins per splice point that are actually materialized and clash-checked.
#:
#: CHOICE. Scoring a spin costs one rotated *point*; checking one costs a rotated segment
#: and a KD-tree query. Scoring all of them and materializing a few is what keeps a
#: thousand-residue assembly to well under a second.
SPINS_MATERIALIZED = 8

#: Arrangement retries before the segment ensemble itself is regenerated.
#:
#: CHOICE. Re-arranging costs a few hundred rotations; regenerating costs ``k`` sub-engine
#: calls. When an arrangement cannot reach the target it is usually the arrangement's
#: fault, so try that first -- but not forever, because a set of segments that are all far
#: too compact cannot be arranged into a long chain at any spin.
ARRANGEMENTS_PER_ENSEMBLE = 4

#: Factor by which the per-segment end-to-end target is inflated after a whole cycle of
#: arrangements failed to close, and the ceiling on the accumulated inflation.
#:
#: CHOICE. Sizing every segment to ``target * (n_seg / N) ** nu`` is the schedule the blob
#: decomposition is derived from, and it is right on average -- but it gives the arrangement
#: no slack, and a compact interior region (900 residues bridging anchors 60 A apart) needs
#: slack: three blobs each spanning exactly a third of the way cannot be interdigitated
#: without a contact. Inflating the segment target and letting the junction search fold the
#: extra length back out recovers those cases. The ceiling stops it from asking for segments
#: near their own contour length, which no sub-engine can deliver.
SEGMENT_STRETCH_STEP = 1.6
MAX_SEGMENT_STRETCH = 4.0

#: Cost, in Angstroms of end-to-end error, at or below which a junction candidate is taken
#: without scoring the rest. CHOICE: a target hit this closely cannot be improved on in any
#: way that matters, and the search is the expensive part of assembly.
GOOD_ENOUGH_COST = 0.5


@dataclass(frozen=True, slots=True)
class SegmentSpan:
    """A half-open range of the *final* chain that one segment covers, ``[start, stop)``.

    Indices are positions within the region being built, 0-based, in the same half-open
    convention as :class:`dodo.structure.Span`. Adjacent segments deliberately overlap, so
    these spans are *not* a tiling -- ``spans[i].stop > spans[i + 1].start`` by exactly the
    splice overlap.
    """

    start: int
    stop: int

    def __len__(self) -> int:
        return self.stop - self.start

    def __str__(self) -> str:
        return f"[{self.start}, {self.stop})"


@dataclass(frozen=True, slots=True)
class Junction:
    """One spliced joint, with the geometry that was measured across it.

    Exists so a caller can state, per junction, that the joint is physically valid --
    rather than trusting that the construction in the module docstring was implemented
    correctly. ``bond_length`` and ``angle_at_splice`` are measured on the finished chain,
    not predicted.
    """

    conformation: int
    index: int
    splice_residue: int
    overlap: tuple[int, int]
    spin_degrees: float
    bond_length: float
    angle_at_splice: float
    clash_cutoff: float

    def __str__(self) -> str:
        return (
            f"junction {self.index} spliced at residue {self.splice_residue} "
            f"(overlap {self.overlap[0]}-{self.overlap[1]}): bond "
            f"{self.bond_length:.3f} A, angle {self.angle_at_splice:.1f} deg"
        )


@dataclass(frozen=True, slots=True)
class _SegmentSet:
    """One conformer per segment, with everything the sub-engine said about them.

    Exists so that ``relaxed_to`` cannot be dropped on the floor. It used to be: only
    ``ca_coords`` and ``success`` were read off each sub-engine result, so a segment the
    walk had built at a relaxed 2.8 A clash threshold -- and correctly said so -- was
    spliced in, and the assembled chain came back claiming ``relaxed_to=None`` with a
    2.934 A CA-CA contact inside it.
    """

    segments: list[np.ndarray]
    #: True when every segment is inside DODO's *generation* window for bonds and angles.
    #: False loosens only the pseudo-angle window, to the MEASURED 75-179 degree observed
    #: range -- never the bond gate, which stays at CA_CA_BOND_TOLERANCE because the
    #: trans-peptide CA-CA distance is rigid and no sub-engine has a legitimate reason to
    #: produce a 4.25 A virtual bond. Clash relaxation deliberately does not clear this flag
    #: either: a contact is not a reason to stop measuring bonds and angles strictly.
    strict: bool
    violations: list[str]
    #: Loosest clash threshold any segment was built at, or ``None`` if all were strict.
    relaxed_to: float | None

    @property
    def clash_floor(self) -> float:
        """The clash distance the assembled chain can honestly be held to."""
        return (
            CA_CLASH_DISTANCE
            if self.relaxed_to is None
            else min(CA_CLASH_DISTANCE, self.relaxed_to)
        )


@dataclass(frozen=True, slots=True)
class _Conformation:
    """One assembled chain and the provenance of the segments that went into it."""

    chain: np.ndarray
    junctions: list[Junction]
    attempts: int
    strict: bool
    violations: list[str]
    relaxed_to: float | None
    clash_floor: float


@dataclass(frozen=True, slots=True)
class AssemblyReport:
    """Everything about a hierarchical build that :class:`IDRResult` has no field for.

    Attributes
    ----------
    sub_engine
        ``name`` of the engine that generated the segments. Recorded because "this chain
        was assembled from walk segments" and "from STARLING segments" are different
        scientific claims and the output must not conflate them.
    delegated
        True when the region fitted under the cap and was passed straight through, so no
        assembly happened at all.
    cap
        The segment length cap used, and whether it came from the installed STARLING, the
        caller, or DODO's constant.
    spans
        The overlapping segment ranges.
    target_end_to_end
        What the arrangement aimed at, in Angstroms. Equals the request's target unless
        the anchors forced it elsewhere (see :func:`~dodo.engines.starling.desired_internal_span`).
    achieved_end_to_end
        Measured value per conformation.
    junctions
        Every junction of every conformation, measured.
    trace_reports
        Whole-chain validation per conformation.
    segment_violations
        Violations found *inside* segments, i.e. attributable to the sub-engine rather than
        to assembly. Reported, not raised: this module cannot fix another engine's geometry,
        but it will not hide it either.
    attempts
        Total assembly attempts across all conformations.
    """

    sub_engine: str
    delegated: bool
    cap: LengthCap
    n_segments: int
    spans: tuple[SegmentSpan, ...]
    target_end_to_end: float
    achieved_end_to_end: tuple[float, ...]
    junctions: tuple[Junction, ...]
    trace_reports: tuple[TraceReport, ...]
    placements: tuple[AnchorPlacement, ...]
    segment_violations: tuple[str, ...]
    attempts: int
    notes: tuple[str, ...]

    def summary(self) -> str:
        """One-line summary for logs and figure captions."""
        if self.delegated:
            return (
                f"hierarchical: no assembly needed ({self.n_segments} segment, cap "
                f"{self.cap}); delegated to {self.sub_engine}"
            )
        achieved = ", ".join(f"{value:.1f}" for value in self.achieved_end_to_end)
        return (
            f"hierarchical: {self.n_segments} {self.sub_engine} segments, cap {self.cap}, "
            f"{len(self.junctions)} junction(s), target "
            f"{self.target_end_to_end:.1f} A, achieved {achieved} A"
        )


# ---------------------------------------------------------------------------
# Splitting
# ---------------------------------------------------------------------------


def segment_spans(n_residues: int, cap: int, overlap: int) -> tuple[SegmentSpan, ...]:
    """Split ``n_residues`` into overlapping segments no longer than ``cap``.

    The count follows from the overlap. ``k`` segments covering ``n`` residues with
    ``overlap`` shared residues at each of the ``k - 1`` joints must generate
    ``n + (k - 1) * overlap`` residues in total, so the longest segment is at least
    ``(n + (k - 1) * overlap) / k``. Requiring that to be at most ``cap`` and solving for
    ``k`` gives

        ``k >= (n - overlap) / (cap - overlap)``

    and ``k`` is the ceiling of that. Note the denominator: each extra segment buys only
    ``cap - overlap`` residues of new coverage, not ``cap``. Sizing on ``cap`` alone --
    the obvious ``ceil(n / cap)`` -- silently produces segments *longer* than the cap once
    the overlap is added back, which for STARLING means a hard error from inside the model.

    Returns
    -------
    tuple[SegmentSpan, ...]
        Spans in chain order. Consecutive spans overlap by exactly ``overlap`` residues;
        the first starts at 0 and the last stops at ``n_residues``. A region at or under
        the cap yields a single span and no overlap.

    Raises
    ------
    EngineError
        If the cap and overlap cannot produce usable segments. Two ways that happens:
        ``cap <= overlap`` (no progress per segment, so no finite ``k`` works), and
        segments too short to carry two junctions -- a splice needs two overlap residues
        before it, so a segment spliced at both ends must be longer than twice the overlap.
    """
    if n_residues <= 0:
        raise EngineError(f"Cannot split a region of {n_residues} residues.")
    if overlap < 2:
        raise EngineError(
            f"splice_overlap must be at least 2, got {overlap}: the two-point junction "
            f"alignment needs two residues present in both segments."
        )
    if cap <= overlap:
        raise EngineError(
            f"A segment cap of {cap} residues with {overlap} residues of overlap makes no "
            f"progress: each extra segment would add {cap - overlap} new residues."
        )
    if n_residues <= cap:
        return (SegmentSpan(0, n_residues),)

    k = -(-(n_residues - overlap) // (cap - overlap))  # ceiling division
    total = n_residues + (k - 1) * overlap
    base, remainder = divmod(total, k)
    lengths = [base + (1 if i < remainder else 0) for i in range(k)]

    minimum = 2 * overlap + 1
    if min(lengths) < minimum:
        raise EngineError(
            f"Splitting {n_residues} residues at a cap of {cap} with {overlap} residues of "
            f"overlap gives segments as short as {min(lengths)}, but a segment spliced at "
            f"both ends must exceed {minimum - 1} residues for its two splice points to "
            f"stay in order. Use a larger cap or a smaller overlap."
        )

    spans: list[SegmentSpan] = []
    start = 0
    for index, length in enumerate(lengths):
        spans.append(SegmentSpan(start, start + length))
        start = start + length - overlap
        if index == len(lengths) - 1 and spans[-1].stop != n_residues:  # pragma: no cover
            raise GeometryError(
                f"Segment split ended at {spans[-1].stop}, not {n_residues}. This is an "
                f"arithmetic bug in segment_spans, not a bad input."
            )
    return tuple(spans)


def resolve_segment_cap(sub_engine: ConformationEngine, override: int | None = None) -> LengthCap:
    """Decide the maximum segment length, and record where the number came from.

    Order of preference, most specific first:

    1. ``override``, when the caller passed one.
    2. The sub-engine's own advertised limit, via a ``max_length`` attribute or method.
       :class:`~dodo.engines.starling.StarlingEngine` has one and reads it from the
       installed STARLING at runtime, so a STARLING release that raises the limit is picked
       up without a DODO change.
    3. :func:`~dodo.engines.starling.starling_max_length`, which probes an installed
       STARLING and otherwise returns :data:`dodo.constants.STARLING_MAX_LENGTH`.

    An engine with no intrinsic limit -- the self-avoiding walk, for instance -- lands on
    (3). That is deliberate rather than "no cap": segmented assembly over an uncapped
    engine is how this module is tested and how a mixed pipeline stays comparable, and a
    caller who wants the walk to build 2,000 residues in one piece should call the walk.
    """
    if override is not None:
        if override < 1:
            raise EngineError(f"max_segment_length must be positive, got {override}.")
        return LengthCap(int(override), "caller")
    probe = getattr(sub_engine, "max_length", None)
    value = probe() if callable(probe) else probe
    if isinstance(value, LengthCap):
        return value
    if isinstance(value, int | np.integer) and not isinstance(value, bool) and int(value) > 0:
        return LengthCap(int(value), f"{getattr(sub_engine, 'name', 'sub-engine')}.max_length")
    return starling_max_length()


def default_sub_engine() -> ConformationEngine:
    """Return STARLING if it can run here, else the self-avoiding walk.

    A convenience for callers assembling long IDRs without caring which generator does the
    segments -- the engine choice is still visible afterwards through
    :attr:`AssemblyReport.sub_engine`.

    Both imports are function-local. :mod:`dodo.engines.walk` is imported here rather than
    at module scope so that this module stays importable and testable on its own, and
    STARLING is never imported at all: :meth:`StarlingEngine.available` answers from a spec
    lookup.
    """
    from .starling import StarlingEngine
    from .walk import SelfAvoidingWalk

    starling = StarlingEngine()
    if starling.available():
        return starling
    return SelfAvoidingWalk()


# ---------------------------------------------------------------------------
# The engine
# ---------------------------------------------------------------------------


class HierarchicalEngine:
    """Assemble a region longer than the sub-engine's cap out of overlapping segments.

    Satisfies the :class:`~dodo.engines.base.ConformationEngine` protocol, and takes one:
    the sub-engine generates the segments. It works with any of them, which is what makes
    it testable where STARLING is not installed -- assembly over walk-generated segments
    exercises every line of the arrangement, the splicing and the validation.

    Parameters
    ----------
    sub_engine
        The engine that generates each segment.
    max_segment_length
        Override the cap. Leave ``None`` to resolve it via :func:`resolve_segment_cap`.
    splice_overlap
        Residues shared by adjacent segments. Defaults to
        :data:`dodo.constants.SEGMENT_SPLICE_OVERLAP`.
    scaling_exponent
        Flory exponent used for the per-segment target and the running-length schedule.
        Defaults to :data:`dodo.constants.FLORY_RE_EXPONENT`.
    spin_samples, spins_materialized
        Search width at each junction. See :data:`SPIN_SAMPLES`.
    max_attempts
        Assembly attempts per conformation before
        :class:`~dodo.exceptions.ExhaustedAttemptsError`. Defaults to
        :data:`dodo.constants.MAX_ATTEMPTS_PER_REGION`.
    end_to_end_tolerance
        Absolute tolerance in Angstroms on the achieved full-length end-to-end distance.
        Defaults to the larger of :data:`END_TO_END_ABSOLUTE_TOLERANCE` and
        :data:`END_TO_END_RELATIVE_TOLERANCE` times the target.
    clash_cutoff
        Minimum CA-CA approach between residues of different segments, in Angstroms.
    """

    name: str = "hierarchical"

    def __init__(
        self,
        sub_engine: ConformationEngine,
        *,
        max_segment_length: int | None = None,
        splice_overlap: int = SEGMENT_SPLICE_OVERLAP,
        scaling_exponent: float = FLORY_RE_EXPONENT,
        spin_samples: int = SPIN_SAMPLES,
        spins_materialized: int = SPINS_MATERIALIZED,
        max_attempts: int = MAX_ATTEMPTS_PER_REGION,
        end_to_end_tolerance: float | None = None,
        clash_cutoff: float = CA_CLASH_DISTANCE,
    ) -> None:
        if not hasattr(sub_engine, "generate"):
            raise EngineError(
                f"sub_engine must be a ConformationEngine with a generate() method, got "
                f"{type(sub_engine).__name__}."
            )
        if spin_samples < 1 or spins_materialized < 1:
            raise EngineError("spin_samples and spins_materialized must both be positive.")
        if max_attempts < 1:
            raise EngineError(f"max_attempts must be positive, got {max_attempts}.")
        self.sub_engine = sub_engine
        self.max_segment_length = max_segment_length
        self.splice_overlap = int(splice_overlap)
        self.scaling_exponent = float(scaling_exponent)
        self.spin_samples = int(spin_samples)
        self.spins_materialized = int(spins_materialized)
        self.max_attempts = int(max_attempts)
        self.end_to_end_tolerance = end_to_end_tolerance
        self.clash_cutoff = float(clash_cutoff)

    def __repr__(self) -> str:
        return (
            f"HierarchicalEngine(sub_engine={getattr(self.sub_engine, 'name', '?')!r}, "
            f"cap={self.cap()}, splice_overlap={self.splice_overlap})"
        )

    def available(self) -> bool:
        """Whether this engine can run, which is whether its sub-engine can."""
        return bool(self.sub_engine.available())

    def cap(self) -> LengthCap:
        """Return the segment length cap in force, with its provenance."""
        return resolve_segment_cap(self.sub_engine, self.max_segment_length)

    # -- generation --------------------------------------------------------

    def generate(
        self,
        request: IDRRequest,
        obstacles: np.ndarray | None,
        rng: np.random.Generator,
    ) -> IDRResult:
        """Assemble ``request.n_conformations`` conformations. See :meth:`generate_detailed`."""
        result, _ = self.generate_detailed(request, obstacles, rng)
        return result

    def generate_detailed(
        self,
        request: IDRRequest,
        obstacles: np.ndarray | None,
        rng: np.random.Generator,
    ) -> tuple[IDRResult, AssemblyReport]:
        """Assemble the region and return the measurements alongside the coordinates.

        Regions at or under the cap are passed straight through to the sub-engine; the
        report says so through :attr:`AssemblyReport.delegated`, so a caller can never
        mistake a delegated build for an assembled one.

        Parameters
        ----------
        request
            The region to build. ``target_end_to_end`` is the *full-length* target and is
            what the arrangement aims at, adjusted only where the anchors demand it.
        obstacles
            ``(m, 3)`` coordinates the finished chain must clear -- typically
            ``structure.xyz[structure.placed_atom_mask()]``, minus the anchor residues'
            own atoms. Applied when the assembled chain is placed on its anchors, not
            during assembly; see the module docstring for why.
        rng
            Seeded generator. Two calls with equally seeded generators produce identical
            coordinates.

        Raises
        ------
        EngineError
            The request is inconsistent, or the cap and overlap cannot produce usable
            segments.
        UnsatisfiableTargetError
            The target end-to-end distance exceeds what the generated segments could span
            even fully extended. Carries both numbers.
        ExhaustedAttemptsError
            No arrangement reached the target within the attempt budget. Deliberately not
            a silently-returned near-miss: the whole point of this module is that the
            full-length dimension is right.
        BuildError
            A finished junction measured as physically invalid, which the construction says
            is impossible. That is a bug in this module and it is raised as one rather than
            noted.
        """
        # IDRRequest validates itself in __post_init__ (region length, sequence agreement,
        # positive target, finite anchors), so none of that is re-checked here.
        if not isinstance(rng, np.random.Generator):
            raise TypeError(
                f"rng must be a numpy.random.Generator, got {type(rng).__name__}. "
                f"Construct one with numpy.random.default_rng(seed)."
            )

        cap = self.cap()
        sub_name = str(getattr(self.sub_engine, "name", type(self.sub_engine).__name__))
        spans = segment_spans(request.n_residues, cap.value, self.splice_overlap)

        if len(spans) == 1:
            result = self.sub_engine.generate(request, obstacles, rng)
            if result.n_successful == 0:
                # The sub-engine's own contract says total failure raises, but this engine is
                # the one whose return value a caller sees, so it does not pass on an
                # all-NaN IDRResult(0/N built) that reads as a result.
                raise ExhaustedAttemptsError(
                    f"{sub_name} built none of the {request.n_conformations} requested "
                    f"conformation(s) of this {request.n_residues}-residue region and "
                    f"returned an all-NaN result instead of raising, so this engine raises "
                    f"on its behalf.",
                    attempts=int(result.attempts),
                )
            report = AssemblyReport(
                sub_engine=sub_name,
                delegated=True,
                cap=cap,
                n_segments=1,
                spans=spans,
                target_end_to_end=float(request.target_end_to_end),
                achieved_end_to_end=_achieved_end_to_end(result.ca_coords),
                junctions=(),
                trace_reports=(),
                placements=(),
                segment_violations=(),
                attempts=int(result.attempts),
                notes=(
                    f"{request.n_residues} residues is within the {cap} cap, so no "
                    f"segmentation was needed; the request went straight to {sub_name}.",
                ),
            )
            return result, report

        # Aim inside the anchor-feasible window with room to spare on both sides, and cap
        # the acceptance tolerance at the room that actually remains. Both halves matter.
        # Clipped onto the bare window edge, the assembly was asked for a value it could
        # only miss outwards, and a 0.163% undershoot of a 292 A target -- against its own
        # 14.6 A tolerance -- put both anchors out of reach and threw the conformation away:
        # 35% of interior builds returned nothing at all. And a tolerance wider than the
        # window is a tolerance that accepts unplaceable arrangements, so it is clamped to
        # the distance from the aim point to the nearer edge.
        target, window = desired_internal_span(
            request.target_end_to_end,
            request.n_anchor_xyz,
            request.c_anchor_xyz,
            slack=self._tolerance(request.target_end_to_end),
        )
        notes: list[str] = [
            f"{request.n_residues} residues exceeds the {cap} cap, so the region was "
            f"assembled from {len(spans)} {sub_name} segments overlapping by "
            f"{self.splice_overlap} residues."
        ]
        if window is not None and not (window[0] <= request.target_end_to_end <= window[1]):
            notes.append(
                f"The requested target ({request.target_end_to_end:.1f} A) lies outside "
                f"what the anchors permit ([{window[0]:.1f}, {window[1]:.1f}] A); the "
                f"assembly aimed at {target:.1f} A."
            )
        tolerance = self._tolerance(target)
        if window is not None:
            # The uncertainty on the dimension belongs to the *requested* target -- a
            # prediction carrying a few percent of error -- not to the smaller number the
            # anchors clipped it to. Scaling the tolerance off the clipped value makes a
            # 900-residue region bridging 60 A-separated anchors answerable only to 3 A when
            # the prediction it came from is uncertain to 11, and every value in the window is
            # placeable anyway. So take the larger scale, then cap it at the room available.
            tolerance = max(tolerance, self._tolerance(request.target_end_to_end))
            room = min(target - window[0], window[1] - target)
            if room < tolerance:
                notes.append(
                    f"The anchor-feasible end-to-end window is [{window[0]:.1f}, "
                    f"{window[1]:.1f}] A, which leaves only {room:.1f} A of room around the "
                    f"{target:.1f} A aim point, so the assembly tolerance was tightened from "
                    f"{tolerance:.1f} A to that. A looser tolerance would accept arrangements "
                    f"the anchors cannot reach."
                )
                tolerance = max(0.0, room)

        coords = np.empty((request.n_conformations, request.n_residues, 3), dtype=np.float64)
        success = np.zeros(request.n_conformations, dtype=bool)
        junctions: list[Junction] = []
        reports: list[TraceReport] = []
        placements: list[AnchorPlacement] = []
        achieved: list[float] = []
        segment_violations: list[str] = []
        relaxations: list[float] = []
        attempts = 0

        for conformation in range(request.n_conformations):
            built = self._assemble_one(
                request=request,
                spans=spans,
                target=target,
                tolerance=tolerance,
                rng=rng,
                conformation=conformation,
                sub_name=sub_name,
            )
            attempts += built.attempts
            junctions.extend(built.junctions)
            for violation in built.violations:
                if violation not in segment_violations:
                    segment_violations.append(violation)
            if built.relaxed_to is not None:
                relaxations.append(built.relaxed_to)
                note = (
                    f"{sub_name} built at least one segment at a relaxed clash threshold of "
                    f"{built.relaxed_to:.2f} A, and that segment is spliced into the chain "
                    f"unchanged, so the assembled chain is a {built.relaxed_to:.2f} A chain. "
                    f"Reported through relaxed_to rather than hidden."
                )
                if note not in notes:
                    notes.append(note)

            # The whole-chain check, at the threshold the result will disclose. Bonds and
            # angles are LOCAL measurements: a chain can satisfy both at every residue while
            # folding back through itself, so without the non-bonded contact check a clash
            # arriving inside a segment passed through undetected and unreported.
            trace_report = self._validate_assembly(
                built.chain,
                built.junctions,
                strict=built.strict,
                clash_distance=built.clash_floor,
            )
            reports.append(trace_report)

            placement = place_between_anchors(
                built.chain,
                n_anchor_xyz=request.n_anchor_xyz,
                c_anchor_xyz=request.c_anchor_xyz,
                rng=rng,
                obstacles=obstacles,
                desired_end_to_end=target,
                clash_cutoff=self.clash_cutoff,
                internal_clash_cutoff=built.clash_floor,
            )
            placements.append(placement)
            coords[conformation] = placement.ca_coords
            # Success needs both: the chain must be placeable on its anchors *and* be
            # physically valid. A conformation whose geometry fails even the observed-range
            # angle window cannot be reconstructed to all-atom, so calling it a success
            # because the placement worked would hand the caller an unusable chain with a
            # True beside it.
            success[conformation] = placement.ok and trace_report.ok
            achieved.append(placement.achieved_end_to_end)
            if placement.relaxed_to is not None:
                relaxations.append(placement.relaxed_to)
            if not placement.ok:
                notes.append(
                    f"Conformation {conformation} assembled to the target but is not usable "
                    f"where it has to sit: worst anchor gap error "
                    f"{placement.anchor_residual:.3f} A, closest internal CA-CA contact "
                    f"{placement.min_internal_ca_distance:.3f} A against a "
                    f"{placement.internal_clash_cutoff:.2f} A floor, obstacles cleared="
                    f"{placement.clash_free}"
                    + (f". {'; '.join(placement.notes)}" if placement.notes else ".")
                )
            if not trace_report.ok:
                notes.append(
                    f"Conformation {conformation} is not a valid CA trace, so it is reported "
                    f"as a failure. The violations are outside every junction window, which "
                    f"means they came in with a {sub_name} segment: {trace_report.describe(3)}"
                )

        if segment_violations:
            notes.append(
                f"{sub_name} segments contained geometry outside DODO's generation window; "
                f"assembly preserved it rather than smoothing it over. First few: "
                f"{'; '.join(segment_violations[:3])}"
            )

        if not bool(success.any()):
            # base.py: partial success goes through the mask, total failure raises. An
            # IDRResult(0/N built) full of NaN is neither, and a caller written against that
            # contract reads it as a result.
            raise ExhaustedAttemptsError(
                f"Every one of the {request.n_conformations} requested conformation(s) of "
                f"this {request.n_residues}-residue region assembled to the target but failed "
                f"validation or placement, so there is nothing to return. {' '.join(notes)}",
                attempts=attempts,
            )

        # from_batch NaN-fills any conformation that could not be placed, which is
        # IDRResult's contract: a failed row must not look like a structure. The real
        # coordinates stay on the corresponding AnchorPlacement for diagnostics.
        result = IDRResult.from_batch(
            ca_coords=coords,
            success=success,
            engine=f"{self.name}({sub_name})",
            attempts=attempts,
            # min, not max: with more than one rung in play the threshold the chain as a
            # whole was actually built at is the loosest of them. max understates the
            # relaxation, which is the opposite of what relaxed_to exists to do.
            relaxed_to=min(relaxations) if relaxations else None,
        )
        report = AssemblyReport(
            sub_engine=sub_name,
            delegated=False,
            cap=cap,
            n_segments=len(spans),
            spans=spans,
            target_end_to_end=target,
            achieved_end_to_end=tuple(achieved),
            junctions=tuple(junctions),
            trace_reports=tuple(reports),
            placements=tuple(placements),
            segment_violations=tuple(segment_violations),
            attempts=attempts,
            notes=tuple(notes),
        )
        return result, report

    # -- one conformation --------------------------------------------------

    def _assemble_one(
        self,
        *,
        request: IDRRequest,
        spans: tuple[SegmentSpan, ...],
        target: float,
        tolerance: float,
        rng: np.random.Generator,
        conformation: int,
        sub_name: str,
    ) -> _Conformation:
        """Generate segments and arrange them until the full-length target is met.

        Returns the chain, its junctions, how many attempts it took, whether the segments
        were strictly valid (which sets the standard the junctions are held to), any
        sub-engine violations seen, and any clash relaxation the sub-engine declared.
        """
        built: _SegmentSet | None = None
        # Kept across the regeneration cycle so the exhaustion message can still name what
        # was wrong with the segments; `built` is cleared on purpose when they are retried.
        last_violations: list[str] = []
        best_error = float("inf")
        best_reach = 0.0
        attempts = 0
        stretch = 1.0

        for attempt in range(1, self.max_attempts + 1):
            attempts = attempt
            if built is None:
                built = self._generate_segments(request, spans, target * stretch, rng, sub_name)
                last_violations = built.violations
                reach = sum(
                    float(np.linalg.norm(segment[-1] - segment[0])) for segment in built.segments
                )
                best_reach = max(best_reach, reach)
                if reach < target:
                    # The arrangement can place segments but cannot stretch them, so a set
                    # whose own end-to-end distances do not sum to the target cannot reach
                    # it at any spin. Segment generation is stochastic, so try another set
                    # rather than declaring the target impossible on one sample.
                    built = None
                    continue
            arranged = self._arrange(built.segments, spans, target, rng, conformation, built.strict)
            if arranged is not None:
                chain, junctions = arranged
                error = abs(float(np.linalg.norm(chain[-1] - chain[0])) - target)
                best_error = min(best_error, error)
                if error <= tolerance:
                    return _Conformation(
                        chain=chain,
                        junctions=junctions,
                        attempts=attempts,
                        strict=built.strict,
                        violations=built.violations,
                        relaxed_to=built.relaxed_to,
                        clash_floor=built.clash_floor,
                    )
            # Re-arranging is cheap and regenerating is not, so alternate: several fresh
            # arrangements of the same segments, then a fresh set of segments.
            if attempt % ARRANGEMENTS_PER_ENSEMBLE == 0:
                built = None
                if best_error == float("inf"):
                    # Not one arrangement has closed: every junction was rejected, which for
                    # a very compact target means the blobs are so tight that no spin puts the
                    # incoming one clear of what is already placed. Blobs sized exactly to the
                    # scaling schedule leave the arrangement nothing to fold, so ask for more
                    # extended ones and let the junction search take the slack back out. The
                    # ceiling in _segment_target keeps this below each segment's contour
                    # length, and the arrangement still has to hit `target` itself.
                    stretch = min(stretch * SEGMENT_STRETCH_STEP, MAX_SEGMENT_STRETCH)

        if best_error == float("inf") and best_reach < target:
            raise UnsatisfiableTargetError(
                f"Every set of segments generated for this {request.n_residues}-residue "
                f"region was too compact to reach a full-length end-to-end distance of "
                f"{target:.1f} A: the best {len(spans)} segments spanned {best_reach:.1f} A "
                f"in total. Either the sub-engine is producing segments far more compact "
                f"than polymer scaling assumes, or this target is not achievable for this "
                f"sequence.",
                target=target,
                achievable=best_reach,
            )
        # Name what was wrong with the segments as well as the distance missed. When a
        # junction is rejected for unbuildable geometry, *every* candidate is rejected and
        # the arrangement dead-ends -- which reads as "could not reach the target" unless the
        # measured violation is carried into the message, and then the number that is
        # actually wrong is invisible.
        blame = (
            f" The {sub_name} segments themselves measured badly, which is the likelier "
            f"reason no junction was accepted: {'; '.join(last_violations[:3])}"
            if last_violations
            else ""
        )
        raise ExhaustedAttemptsError(
            f"Hierarchical assembly of {request.n_residues} residues from "
            f"{len(spans)} {sub_name} segments could not reach an end-to-end distance of "
            f"{target:.1f} +/- {tolerance:.1f} A; the closest arrangement was off by "
            f"{best_error:.1f} A. Returning it anyway would mean returning a chain with the "
            f"wrong dimensions, which is the failure this engine exists to prevent." + blame,
            attempts=attempts,
        )

    def _generate_segments(
        self,
        request: IDRRequest,
        spans: tuple[SegmentSpan, ...],
        target: float,
        rng: np.random.Generator,
        sub_name: str,
    ) -> _SegmentSet:
        """Generate one conformer per segment, and check what came back.

        Each segment gets its own end-to-end target from the same polymer scaling that
        justifies the decomposition: a segment holding ``n`` of ``N`` residues should span
        ``target * (n / N) ** nu``. ``target`` is the value the *arrangement* aims at -- the
        anchor-clipped one, not ``request.target_end_to_end``. Sizing the segments from the
        unclipped prediction while arranging them to the clipped target is asking for parts
        that do not fit the whole: for an interior region whose anchors sit much closer than
        the prediction, segments came back several times too extended and no arrangement of
        them could be folded down, so requests that were comfortably satisfiable raised.

        Segments are generated in *free space*, with no anchors and no obstacles, because
        they are rigidly re-oriented afterwards -- avoiding an obstacle in a frame the
        segment will not end up in is worse than not trying.
        """
        segments: list[np.ndarray] = []
        violations: list[str] = []
        strict = True
        relaxed_to: float | None = None
        for index, span in enumerate(spans):
            length = len(span)
            sub_request = IDRRequest(
                sequence=request.sequence[span.start : span.stop],
                n_residues=length,
                target_end_to_end=self._segment_target(target, length, request.n_residues),
                n_anchor_xyz=None,
                c_anchor_xyz=None,
                n_conformations=1,
            )
            result = self.sub_engine.generate(sub_request, None, rng)
            coords = np.asarray(result.ca_coords, dtype=np.float64)
            if coords.ndim == 3:
                if not bool(np.asarray(result.success).reshape(-1)[0]):
                    raise BuildError(
                        f"{sub_name} reported failure for segment {index} "
                        f"({span}) of a hierarchical build. Its coordinates are documented "
                        f"as meaningless when success is False, so they are not used."
                    )
                coords = coords[0]
            if coords.shape != (length, 3):
                raise BuildError(
                    f"{sub_name} returned coordinates of shape {coords.shape} for segment "
                    f"{index} ({span}), but {length} residues were requested."
                )
            # A declared relaxation is a fact about the segment, and the segment is spliced
            # in unchanged, so it is a fact about the assembled chain too. Carrying it up is
            # the whole of relaxed_to's purpose; dropping it made a 2.934 A contact arrive in
            # a chain reporting relaxed_to=None and trace_report.ok=True.
            declared = None if result.relaxed_to is None else float(result.relaxed_to)
            if declared is not None and declared < CA_CLASH_DISTANCE:
                relaxed_to = declared if relaxed_to is None else min(relaxed_to, declared)
            # Validate at the threshold the sub-engine itself declared, so a legitimately
            # relaxed segment is not double-counted as a defect -- but a contact tighter than
            # anything it declared is a defect and is reported as one.
            report = validate_ca_trace(
                coords,
                clash_distance=CA_CLASH_DISTANCE if declared is None else declared,
                residue_offset=span.start,
            )
            local = [v for v in report.violations if v.kind != "steric_clash"]
            if local:
                # Only bond and angle violations loosen the standard the junctions are held
                # to, and only for angles. A contact says nothing about whether this
                # sub-engine's angles can be trusted, so it must not buy the segment a wider
                # angle window either.
                strict = False
                violations.append(f"segment {index} {span}: {local[0].message}")
            for clash in report.of_kind("steric_clash"):
                violations.append(f"segment {index} {span}: {clash.message}")
            segments.append(coords)
        return _SegmentSet(
            segments=segments, strict=strict, violations=violations, relaxed_to=relaxed_to
        )

    def _segment_target(self, target: float, n_segment: int, n_total: int) -> float:
        """Per-segment end-to-end target, from the scaling in the module docstring."""
        if n_segment < 2:
            return 0.0
        scaled = float(target) * (n_segment / n_total) ** self.scaling_exponent
        # A segment cannot span more than its own contour length, and asking for that
        # would make the sub-engine burn its whole retry budget on an impossible target.
        ceiling = 0.95 * contour_length(n_segment)
        return float(min(scaled, ceiling))

    # -- the arrangement ---------------------------------------------------

    def _arrange(
        self,
        segments: Sequence[np.ndarray],
        spans: Sequence[SegmentSpan],
        target: float,
        rng: np.random.Generator,
        conformation: int,
        strict: bool,
    ) -> tuple[np.ndarray, list[Junction]] | None:
        """Splice the segments together, steering the running end toward the target.

        Returns ``None`` when some junction admitted no valid candidate at all, which the
        caller treats as "try again" rather than as a failure -- a dead end here is a
        property of this particular set of segments and spins, not of the request.
        """
        n_residues = spans[-1].stop
        chain = np.empty((n_residues, 3), dtype=np.float64)
        chain[spans[0].start : spans[0].stop] = segments[0]
        junctions: list[Junction] = []

        # Reach still available after each junction: what the remaining segments could add
        # to the distance between the chain's current end and its final end.
        spans_reach = [float(np.linalg.norm(segment[-1] - segment[0])) for segment in segments]

        for index in range(1, len(segments)):
            overlap = (spans[index].start, spans[index - 1].stop)
            # What the segments after this one could still add to the distance between the
            # chain's current end and its final end: their own end-to-end distances, plus
            # one bond per junction, which the sum over segment spans leaves out.
            remaining = sum(spans_reach[index + 1 :]) + CA_CA_BOND_LENGTH * (
                len(segments) - index - 1
            )
            is_last = index == len(segments) - 1
            desired_running = (
                target
                if is_last
                else target * (spans[index].stop / n_residues) ** self.scaling_exponent
            )
            chosen = self._best_junction(
                chain=chain,
                segment=segments[index],
                span=spans[index],
                overlap=overlap,
                desired_running=desired_running,
                final_target=target,
                remaining_reach=remaining,
                rng=rng,
                strict=strict,
            )
            if chosen is None:
                return None
            splice_residue, spin, tail = chosen
            chain[splice_residue : spans[index].stop] = tail
            junctions.append(
                Junction(
                    conformation=conformation,
                    index=index,
                    splice_residue=splice_residue,
                    overlap=overlap,
                    spin_degrees=float(np.rad2deg(spin)),
                    bond_length=float(
                        np.linalg.norm(chain[splice_residue] - chain[splice_residue - 1])
                    ),
                    angle_at_splice=_angle_at(chain, splice_residue - 1),
                    clash_cutoff=self.clash_cutoff,
                )
            )
        return chain, junctions

    def _best_junction(
        self,
        *,
        chain: np.ndarray,
        segment: np.ndarray,
        span: SegmentSpan,
        overlap: tuple[int, int],
        desired_running: float,
        final_target: float,
        remaining_reach: float,
        rng: np.random.Generator,
        strict: bool,
    ) -> tuple[int, float, np.ndarray] | None:
        """Choose a splice point and spin for one junction.

        Splice points are tried from the middle of the overlap outwards -- the middle
        leaves the most generated context on both sides, and a splice at the very edge of
        the overlap uses residues that only just exist in both segments.

        For each splice point every spin is *scored* by transforming a single point (the
        segment's far end), which is enough to know the resulting running end-to-end
        distance; only the best few are materialized and clash-checked. See
        :data:`SPIN_SAMPLES`.
        """
        overlap_start, overlap_stop = overlap
        # A splice at m needs residues m-1 and m-2 present in both segments, so m-2 must
        # be at or after the start of the overlap; and the outgoing chain must reach m-1,
        # so m is at most the end of the overlap.
        candidates = [m for m in range(overlap_start + 2, overlap_stop + 1) if m < span.stop]
        if not candidates:
            raise GeometryError(
                f"Overlap [{overlap_start}, {overlap_stop}) admits no splice point. The "
                f"overlap must be at least 2 residues wide; segment_spans should have "
                f"caught this."
            )
        middle = (overlap_start + overlap_stop) // 2
        candidates.sort(key=lambda m: (abs(m - middle), m))

        best: tuple[float, int, float, np.ndarray] | None = None
        phase = float(rng.uniform(0.0, 2.0 * np.pi))
        spins = phase + np.linspace(0.0, 2.0 * np.pi, self.spin_samples, endpoint=False)

        for splice in candidates:
            local = splice - span.start
            pivot = chain[splice - 1]
            axis_vector = chain[splice - 2] - pivot
            axis_norm = float(np.linalg.norm(axis_vector))
            if axis_norm < _DEGENERATE_LENGTH:  # pragma: no cover - would be a broken chain
                continue
            axis = axis_vector / axis_norm
            aligned = self._align_segment(segment, local, pivot, axis_vector)
            if aligned is None:
                continue
            tail_local = aligned[local:]

            # Score every spin on one point: where the far end of the tail lands.
            ends = _spin_points(tail_local[-1] - pivot, axis, spins) + pivot
            running = np.linalg.norm(ends - chain[0], axis=1)
            cost = np.abs(running - desired_running)
            # Reject spins from which the FINAL target is out of reach for the segments
            # that remain: the end can move by at most `remaining_reach` from here, so the
            # target must lie within that of the current running distance. This is the
            # walk's reachability funnel, one level up. Not applied at the last junction,
            # where `remaining_reach` is zero and `cost` is already the final error.
            if remaining_reach > 0.0:
                reachable = np.abs(running - final_target) <= remaining_reach
                cost = np.where(reachable, cost, np.inf)
            order = np.argsort(cost, kind="stable")[: self.spins_materialized]

            for position in order:
                if not np.isfinite(cost[position]):
                    break
                spin = float(spins[position])
                rotation = rotation_from_axis_angle(axis, spin)
                tail = apply(tail_local - pivot, rotation) + pivot
                if not self._junction_valid(chain, splice, tail, strict):
                    continue
                if not self._tail_clear(chain, splice, tail):
                    continue
                score = float(cost[position])
                if best is None or score < best[0]:
                    best = (score, splice, spin, tail)
                if score <= GOOD_ENOUGH_COST:
                    return best[1], best[2], best[3]
                break
            if best is not None and best[0] <= GOOD_ENOUGH_COST:
                break

        if best is None:
            return None
        return best[1], best[2], best[3]

    def _align_segment(
        self, segment: np.ndarray, local: int, pivot: np.ndarray, axis_vector: np.ndarray
    ) -> np.ndarray | None:
        """Two-point alignment of an incoming segment onto the placed chain.

        Maps the segment's residue ``m - 1`` onto the chain's, and its ``m-1 -> m-2``
        direction onto the chain's. This is what makes the junction bond and both junction
        angles come from a single generator-produced trace; see the module docstring.
        """
        own_axis = segment[local - 2] - segment[local - 1]
        if float(np.linalg.norm(own_axis)) < _DEGENERATE_LENGTH:  # pragma: no cover
            return None
        rotation = rotation_between_vectors(own_axis, axis_vector)
        aligned: np.ndarray = apply(segment - segment[local - 1], rotation) + pivot
        return aligned

    def _junction_valid(
        self, chain: np.ndarray, splice: int, tail: np.ndarray, strict: bool
    ) -> bool:
        """Measure the spliced joint rather than assume it.

        The construction guarantees a valid bond and valid angles, so this is a cheap
        assertion in the inner loop rather than a filter that is expected to fire. It is
        still checked every time: "guaranteed by construction" is what the pre-rewrite code
        also believed about its junctions, which were 0.00 A apart.
        """
        window = np.vstack((chain[max(0, splice - 3) : splice], tail[: min(2, tail.shape[0])]))
        if window.shape[0] < 2:  # pragma: no cover - splice >= 2 guarantees 3 rows
            return True
        report = validate_ca_trace(
            window,
            angle_window=None
            if strict
            else (BACKBONE_ANGLE_OBSERVED_MIN, BACKBONE_ANGLE_OBSERVED_MAX),
        )
        return report.ok

    def _tail_clear(
        self,
        chain: np.ndarray,
        splice: int,
        tail: np.ndarray,
    ) -> bool:
        """Whether a candidate tail clears the part of the chain already placed.

        Pairs closer than :data:`dodo.constants.CLASH_EXCLUDE_WITHIN_RESIDUES` in sequence
        are exempt: residues ``m`` and ``m-1`` are covalently bonded at 3.8 A and ``m`` and
        ``m-2`` are angle-constrained to 5-7.5 A, so counting either as a clash would
        reject every correct splice.

        External obstacles are deliberately *not* considered here. Assembly happens in the
        first segment's own frame, and the finished chain is then moved onto its anchors as
        a rigid body, so an obstacle avoided during assembly would have been avoided in a
        frame the chain does not end up occupying -- worse than not trying, because it looks
        like it worked. Obstacle avoidance happens once, at that final placement, where the
        coordinates are the ones that will be written out.
        """
        placed = chain[:splice]
        if placed.shape[0]:
            tree = cKDTree(placed)
            for offset, hits in enumerate(tree.query_ball_point(tail, self.clash_cutoff)):
                if any(
                    (splice + offset) - int(hit) > CLASH_EXCLUDE_WITHIN_RESIDUES for hit in hits
                ):
                    return False
        return True

    # -- validation --------------------------------------------------------

    def _validate_assembly(
        self,
        chain: np.ndarray,
        junctions: Sequence[Junction],
        *,
        strict: bool,
        clash_distance: float = CA_CLASH_DISTANCE,
    ) -> TraceReport:
        """Validate the finished chain and attribute every violation.

        Three measurements, not two. Bonds and angles are purely local, so on their own they
        cannot see a chain that folds back through itself: a planar hexagon of side 3.80 A
        traversed twice has every bond at exactly 3.800 A and every angle at exactly 120
        degrees while residues six apart occupy the same point. ``clash_distance`` is
        therefore checked too, at the threshold this conformation will *disclose* through
        :attr:`~dodo.engines.base.IDRResult.relaxed_to`, so the returned chain is valid at
        the number the result reports and not merely at a number nobody stated.

        A bond or angle violation inside a junction window is this module's fault and raises.
        A contact is not: the two-point alignment guarantees the junction bond and both
        junction angles, and says nothing about global self-avoidance, so a clash whose lower
        residue happens to fall in a junction window is not evidence of a splicing bug.
        Cross-segment contacts are already enforced at :attr:`clash_cutoff` by
        :meth:`_tail_clear`, which means a surviving contact came in inside a segment; it is
        left in the returned report for the caller to see, because silently smoothing another
        engine's geometry would hide which engine produced what.
        """
        report = validate_ca_trace(
            chain,
            angle_window=None
            if strict
            else (BACKBONE_ANGLE_OBSERVED_MIN, BACKBONE_ANGLE_OBSERVED_MAX),
            clash_distance=clash_distance,
        )
        if report.ok:
            return report
        junction_residues = {
            residue
            for junction in junctions
            for residue in (
                junction.splice_residue - 2,
                junction.splice_residue - 1,
                junction.splice_residue,
            )
        }
        offenders = [
            v
            for v in report.violations
            if v.residue_index in junction_residues and v.kind != "steric_clash"
        ]
        if offenders:
            raise BuildError(
                f"Splicing produced invalid geometry at "
                f"{len(offenders)} junction residue(s), which the two-point alignment is "
                f"supposed to make impossible. This is a bug in hierarchical assembly, not "
                f"bad input. First: {offenders[0].message}"
            )
        return report

    def _tolerance(self, target: float) -> float:
        """Absolute tolerance on the achieved full-length end-to-end distance."""
        return end_to_end_tolerance(target, self.end_to_end_tolerance)


# ---------------------------------------------------------------------------
# Small numerics
# ---------------------------------------------------------------------------

#: Below this length a vector has no usable direction, in Angstroms.
_DEGENERATE_LENGTH = 1e-9


def _spin_points(vector: np.ndarray, axis: np.ndarray, angles: np.ndarray) -> np.ndarray:
    """Rotate one vector about ``axis`` by many angles at once, via Rodrigues' formula.

    ``(n_angles, 3)``. Used only to *score* candidate spins: knowing where the far end of
    a segment lands does not require rotating the whole segment, and scoring cheaply is
    what makes a dense spin scan affordable. The spin that wins is then applied with
    :func:`dodo.geometry.transforms.rotation_from_axis_angle`, so the coordinates that
    survive come from the one audited rotation implementation. The two agree to machine
    precision, which the tests check.
    """
    cos = np.cos(angles)[:, None]
    sin = np.sin(angles)[:, None]
    parallel = float(np.dot(vector, axis)) * axis
    perpendicular = vector - parallel
    cross = np.cross(axis, perpendicular)
    rotated: np.ndarray = parallel[None, :] + perpendicular[None, :] * cos + cross[None, :] * sin
    return rotated


def _achieved_end_to_end(ca_coords: np.ndarray) -> tuple[float, ...]:
    """End-to-end distance per conformation of an ``(n_conf, n_residues, 3)`` array.

    Accepts a single ``(n_residues, 3)`` conformation too, because an engine is only
    documented to return the batched shape and this reads another engine's output.
    """
    array = np.asarray(ca_coords, dtype=np.float64)
    if array.ndim == 2:
        array = array[None, :, :]
    if array.ndim != 3 or array.shape[2] != 3:
        raise GeometryError(
            f"Expected (n_conformations, n_residues, 3) coordinates, got {array.shape}."
        )
    if array.shape[1] < 2:
        zeros: tuple[float, ...] = (0.0,) * int(array.shape[0])
        return zeros
    distances: np.ndarray = np.linalg.norm(array[:, -1] - array[:, 0], axis=1)
    values: list[float] = [float(value) for value in distances]
    return tuple(values)


def _angle_at(chain: np.ndarray, vertex: int) -> float:
    """CA-CA-CA pseudo-angle at ``vertex``, in degrees, or NaN at a chain terminus."""
    if vertex < 1 or vertex + 1 >= chain.shape[0]:
        return float("nan")
    first = chain[vertex - 1] - chain[vertex]
    second = chain[vertex + 1] - chain[vertex]
    norms = float(np.linalg.norm(first)) * float(np.linalg.norm(second))
    if norms < _DEGENERATE_LENGTH:  # pragma: no cover - coincident CAs
        return float("nan")
    cosine = float(np.clip(np.dot(first, second) / norms, -1.0, 1.0))
    return float(np.rad2deg(np.arccos(cosine)))
