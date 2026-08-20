"""The self-avoiding walk engine: cone-constrained growth with a closure funnel.

This is DODO's dependency-free conformation engine. It grows a CA trace one residue at a
time from a cone of candidates around the previous two alpha carbons, rejects candidates
that clash or that would make the region's endpoint unreachable, and picks from what
survives with a physically weighted random choice.

Failure is loud
---------------
Nothing here returns a partial or degenerate chain. A conformer that dead-ends is marked
failed and retried with fresh randomness up to
:data:`~dodo.constants.MAX_ATTEMPTS_PER_REGION` rounds; if nothing survives, the engine
raises :class:`~dodo.exceptions.ExhaustedAttemptsError`. A request that is geometrically
impossible raises before any sampling happens, because spending 40 attempt rounds on an
unsatisfiable target and then returning zeros is what the old code did.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Final

import numpy as np
from scipy.spatial import cKDTree

from ..constants import (
    BACKBONE_ANGLE_MAX,
    BACKBONE_ANGLE_MEAN,
    BACKBONE_ANGLE_MIN,
    BACKBONE_ANGLE_SD,
    CA_CA_BOND_LENGTH,
    CA_CLASH_DISTANCE,
    CANDIDATES_PER_ANGLE,
    CLASH_EXCLUDE_WITHIN_RESIDUES,
    MAX_ATTEMPTS_PER_REGION,
    MAX_CANDIDATES_PER_RESIDUE,
    WALK_BATCH_SIZE,
    backbone_angle_grid,
)
from ..exceptions import (
    EngineError,
    ExhaustedAttemptsError,
    GeometryError,
    UnsatisfiableTargetError,
)
from ..geometry.metrics import end_to_end, validate_ca_trace
from ..geometry.sampling import cone_candidates_batch, random_unit_vectors
from ..geometry.transforms import align_frame
from ..geometry.transforms import apply as apply_rotation
from .base import IDRRequest, IDRResult

__all__ = [
    "SelfAvoidingWalk",
    "UnconstrainedJunctionWarning",
    "max_reach",
    "min_reach",
    "reachability_tail",
    "reachable_envelope",
    "sample_end_to_end_targets",
]

#: Identifier recorded in :attr:`dodo.engines.base.IDRResult.engine`.
ENGINE_NAME: Final[str] = "self_avoiding_walk"

_PARALLEL_QUERY_MINIMUM: Final[int] = 2000

#: Clearance around a candidate cloud's centre below which the cloud cannot clash, in A.
_CANDIDATE_CLEAR_DISTANCE: Final[float] = CA_CA_BOND_LENGTH + CA_CLASH_DISTANCE


class UnconstrainedJunctionWarning(UserWarning):
    """A junction pseudo-angle could be neither steered nor measured.

    Raised as a warning, not an exception, because a region with no neighbour beyond its
    anchor is still buildable and the pipeline still has to build it. But it is *not*
    silent: 41% of junction angles built without the neighbour coordinate fell below
    :data:`dodo.constants.BACKBONE_ANGLE_MIN`, and an unmeasured angle looks exactly like
    a valid one to every check downstream.
    """


#: Fractional tolerance on the achieved end-to-end distance of a free or terminal region.
#:
#: 10% is deliberately loose: ALBATROSS's own end-to-end predictions are good to a
#: few percent, and the analytical fallback
#: (:func:`dodo.constants.flory_end_to_end`) is measured at 0.62-2.62x depending on
#: composition. Demanding better agreement than that would be precision theatre.
#:
#: It is applied as a pure fraction, with no absolute floor.
END_TO_END_TOLERANCE_FRACTION: Final[float] = 0.10

#: Width of the closure steering weight in the safe direction, as a multiple of
#: ``CA_CA_BOND_LENGTH * sqrt(bonds left)``.
#:
#: CHOICE, and MEASURED to be inert -- the value does not matter, only that it is not tiny.
#: Swept over an interior bridge at extensions 0.60 and 0.88, output is bit-identical from
#: **0.05 to 1e6**: success 1.000, mean pairwise RMSD 14.6163 A and achieved/requested 0.9855 at
#: every value in that range. It only begins to bind below about 0.05, 15x under this default.
#:
#: That is the term working as intended rather than a defect. It is the width in the *safe*
#: direction, and it is applied as a Gaussian over the span between the schedule and the corridor
#: floor -- a span far narrower than ``0.75 * 3.81 * sqrt(k)`` (20 A at k = 50), so the weight is
#: flat across every candidate and expresses no preference, which is exactly the "essentially
#: unbiased" behaviour described below. Do not spend time tuning it; if you need the closure to
#: behave differently, :data:`SCHEDULE_MARGIN_FRACTION` is the lever that does something.
#:
#: Applies to a *closure*, where the far anchor is a fixed point and the schedule is
#: seeded from the real anchor separation. The scaling is the point rather than the constant:
#: the span of a ``k``-bond subchain fluctuates as ``sqrt(k)``, so a width proportional to
#: ``sqrt(k)`` leaves the walk essentially unbiased while it has plenty of chain left --
#: which is what makes the output a coil rather than a funnel-shaped artifact -- and tightens
#: automatically into real steering as the remaining chain runs out.
#:
#: It does *not* apply to a free or terminal region, where the far endpoint is a target
#: distance rather than a fixed point. Steering that with a ``sqrt(k)`` width is what made
#: the achieved end-to-end distance settle onto the corridor floor at 0.95 of the target for
#: every mode and length; see :data:`TARGET_STEERING_WIDTH`.
STEERING_SIGMA_FACTOR: Final[float] = 0.75

#: Fraction of the remaining headroom used as the closure steering width in the
#: unrecoverable direction. CHOICE, and MEASURED to be nearly inert. Small enough that the walk
#: tracks its schedule with margin to spare rather than grinding along the reachability ceiling,
#: where the discrete candidate set runs out of legal moves; large enough that the walk is not
#: pinned to a curve.
#:
#: Swept 0.10 to 1.00 over six interior bridges (extensions 0.30-0.90, n = 20-120, four seeds):
#: success stays 1.000 throughout, mean pairwise RMSD moves only 15.68-16.78 A and
#: achieved/requested 0.9850-0.9905, with no monotone trend -- a 10x change in the constant is
#: worth less than the seed-to-seed scatter. It is partly floored out by construction: the width
#: is ``max(fraction * gap, 0.5 * CA_CA_BOND_LENGTH)``, so wherever the gap is small the floor,
#: not this fraction, sets the width.
HEADROOM_STEERING_FRACTION: Final[float] = 0.35

#: Fraction of the headroom between the closure schedule and the hard reachability bound
#: that the walk is actually allowed to use.
#:
#: CHOICE, and the difference between a funnel that works and one that does not. The
#: reachability bound from :func:`max_reach` is a cliff: a walk sitting exactly on it has
#: no legal move, because closing the remaining distance requires the perfect zig-zag and
#: the candidate set is discrete. A diffusive walk offered the whole headroom will reach
#: that cliff -- measured, a taut bridge (extension fraction 0.78) lagged to the bound by
#: residue 18 of 50 and dead-ended at residue 20, every seed.
#:
#: Restricting it to half the headroom turns the boundary from absorbing into restoring:
#: from the policy boundary the most inward step available advances by ``b * sin(theta /
#: 2)``, which for the top of the angle window is 3.75 A, while the boundary itself
#: recedes by only ``extension_fraction * 3.75`` A. The walk is therefore pushed back
#: toward its schedule instead of being stuck against a wall. Measured over a 50-residue
#: bridge at ten seeds each, this takes extension fractions of 0.05 through 0.97 from
#: "0.78 and above never builds" to "every fraction builds, almost always on the first
#: attempt".
#: SWEPT 2026-08-18, and this is the one closure constant that does anything. Over an interior
#: bridge at extensions 0.80-0.97 (six conformers, five seeds), success rate by value:
#:
#:     margin   x0.80  x0.88  x0.94  x0.96  x0.97
#:     0.25     1.00   1.00   1.00   1.00   1.00
#:     0.50     1.00   1.00   1.00   1.00   1.00
#:     0.75     1.00   1.00   1.00   1.00   1.00
#:     0.90     1.00   1.00   1.00   1.00   0.00
#:     1.00     0.00   0.00   0.00   0.00   0.00
#:
#: So the cliff is real and closer than the note above implies: **1.0 fails at every extension
#: from about 0.70 upward**, not merely at 0.78, and 0.90 is already inside the failure region at
#: 0.97. The usable band is roughly 0.25-0.75 and 0.5 sits in the middle of it.
#:
#: Raising it to 0.75 was considered and REJECTED. On the synthetic bridge it looks like a free
#: win -- diversity up 29-31% (mean pairwise RMSD 14.3 vs 10.9 A at extension 0.80) at unchanged
#: success -- but that does not survive contact with real structures, where closure competes with
#: folded-domain obstacles and most residues sit in terminal regions whose diversity comes from
#: target sampling instead. Measured on dnmt3a/arf19/p300 at five models: zero failures either
#: way, clashes 4/1/7 against 3/0/8, diversity 104/143/136 A against 110/145/126 -- mixed in both
#: directions and inside single-seed noise. No reason to move a stabilized constant for that.
SCHEDULE_MARGIN_FRACTION: Final[float] = 0.5

#: FLOOR on the steering width for a free or terminal region, in A. The width itself is
#: computed per conformer by :func:`_target_steering_width`; this is only its lower bound.
#:
#: **Read this first, because the rest of this note predates the change.** Until 2026-08-19 this
#: constant *was* the width, applied flat at every step of every steered region, and everything
#: below is the record of tuning it in that role. It is kept because the measurements are still
#: true and still constrain the floor. What it did not measure -- and what turned out to matter
#: more -- is the *local* geometry the width produces. A flat 0.5 A pins the distance to residue
#: 0 so hard that the pseudo-angle preference and the radial preference together admit only two
#: bond directions per step, mirror images across a slowly-turning plane, and the chain comes out
#: as a flat zig-zag. See :data:`STEERING_SLACK_FRACTION` for the measurement, the mechanism and
#: the replacement.
#:
#: MEASURED, and the single number that decides whether the achieved end-to-end distance
#: matches the requested one. It has to be compared against the *step*, not against the
#: chain: one bond changes the distance to the endpoint by at most 7.5 A end to end, so a
#: width of 20 A -- which is what ``0.75 * b * sqrt(index)`` came to at mid-chain -- expresses
#: essentially no preference between advancing 1 A and advancing 3.5 A. The walk then fell
#: behind its schedule by about 0.2 A per bond, and over 99 bonds that is the whole 5%
#: one-sided miss that used to make every extended conformer land on the corridor floor.
#:
#: Measured ratio of achieved to requested end-to-end distance at 0.963 of the geometric
#: reach -- the most extended target :mod:`dodo.construct.dimensions` can produce -- for
#: widths of 0.4 / 0.6 / 0.75 / 1.0 / 2.85 A:
#:
#:     n = 50 : 0.993 / 0.985 / 0.979 / 0.966 / 0.905
#:     n = 200: 0.998 / 0.995 / 0.992 / 0.987 / 0.946
#:
#: That sweep bottomed out at 0.4 and its conclusion -- that 0.5 sat in the flat part of the
#: curve and tighter bought nothing -- did not survive being tested below its own floor; see the
#: re-sweep below. What tightening does *not* cost is unchanged: the candidates that achieve a given
#: distance from the endpoint form a ring rather than a point, so pinning the radial
#: coordinate this hard leaves the transverse freedom -- and therefore the conformational
#: diversity -- untouched. Measured over 20 conformers at the predicted dimension, the
#: internal scaling exponent is 0.53-0.64 and the coefficient of variation of the end-to-end
#: distance is 0.48, both squarely in the polymer range.
#: RE-SWEPT 2026-08-18, 12 seeds with standard errors, and the note above needs qualifying: the
#: sweep it cites never tested a width below **0.4**, and never tested 0.5 itself -- 0.5 was
#: interpolated between 0.4 and 0.6. "Tighter buys nothing measurable" was therefore an
#: extrapolation past the edge of the data, and it does not hold. Accuracy keeps improving
#: monotonically below 0.4:
#:
#:     width          0.15     0.20     0.30     0.40     0.50     0.75     1.00     2.00
#:     achieved/req  0.9998   0.9997   0.9988   0.9979   0.9969   0.9930   0.9874   0.9706
#:     std error     0.0001   0.0001   0.0002   0.0003   0.0005   0.0010   0.0018   0.0036
#:
#: An earlier revision of this note claimed tightening costs about 6% of inter-conformer RMSD.
#: That was seed scatter, not signal: with error bars the diversity is flat across the whole
#: range (254.6 +- 24.9 A at 0.15 against 257.8 +- 25.3 A at 0.5), and on real structures it does
#: not move either (dnmt3a/arf19/p300 at six models: 102/139/136 A at 0.5 against 104/139/141 A
#: at 0.2).
#:
#: Nor does tightening distort the polymer statistics, which is the claim that would actually
#: matter. Measured at the *predicted* dimension -- the regime this ships in -- 24 conformers of a
#: 120-residue region give an internal scaling exponent of 0.583 at width 0.5 and 0.586 at 0.2,
#: with CV(Re) 0.488 either way: indistinguishable, and both sitting on the self-avoiding value of
#: ~0.588. (Measure that at an over-extended target instead and the exponent rises to ~0.87 for
#: every width alike -- that is the target straightening the chain, not the steering width.)
#:
#: The accuracy gain is real but concentrated where the target is extended: negligible at the
#: predicted dimension (1.0008 against 1.0000) and 1.02% -> 0.15% at 0.963 of the geometric reach.
#:
#: **0.2 was tried on that basis and REVERTED, because tightening costs steric clashes.** The
#: sweep above scores accuracy and diversity, and those were the wrong metrics to decide on: a
#: more tightly steered walk has less freedom left to dodge an obstacle, so marginal contacts
#: survive into the backbone stage. Measured over dnmt3a/arf19/p300 at four seeds, introduced
#: clashes total **21 at width 0.5 against 25 at 0.2**, and p300 alone goes 3 -> 7 at seed 0,
#: past the frozen ratchet in ``tests/unit/test_pipeline.py``. A 0.2% end-to-end difference sits
#: about a hundredfold inside the 10% tolerance and is invisible in any output; four extra steric
#: contacts are physical defects a user can see. So 0.5 stands, now for a measured reason rather
#: than an interpolated one.
TARGET_STEERING_WIDTH: Final[float] = 0.5

#: ``<Re> / sqrt(<Re^2>)`` for the radial distribution of a Gaussian chain.
#:
#: DERIVED. For ``P(R) proportional to R^2 * exp(-3 R^2 / (2 <R^2>))`` -- the
#: Maxwell-Boltzmann-like distribution of the end-to-end *distance* of an ideal chain,
#: which is what you get when the end-to-end *vector* is an isotropic Gaussian -- the mean
#: is ``sqrt(8 / (3 * pi)) * sqrt(<R^2>)``. Dividing the requested target by this factor
#: gives the scale whose sampled mean is the target, which is the whole point: the request
#: names a mean, so the ensemble must reproduce it as a mean.
MAXWELL_MEAN_FACTOR: Final[float] = math.sqrt(8.0 / (3.0 * math.pi))

#: Largest fraction of the geometric reach a per-conformer target is allowed to ask for.
#:
#: MEASURED. The high tail of the Maxwell distribution reaches past what a discrete
#: cone-constrained walk can build: a target approaching :func:`max_reach` needs the perfect
#: all-trans zig-zag, and with ``require_all`` one unlucky draw would fail the whole batch.
#: Sampled targets are clipped here and the ensemble mean is restored on the ones still free
#: to move.
#:
#: At 0.98 of the reach, tails of 10, 20, 50, 100, 200 and 380 residues each built all ten
#: requested conformers and achieved 0.95-0.99 of the requested distance -- inside the 10%
#: tolerance, so nothing is being waved through. Note that
#: :mod:`dodo.construct.dimensions` clamps its own targets to ``0.95 * contour_length``,
#: which is 0.963 of :func:`max_reach`: this constant has to sit above that or a legitimate
#: clamped target would be refused as unsatisfiable.
MAX_TARGET_EXTENSION_FRACTION: Final[float] = 0.98

#: Steering width as a fraction of the radial slack the target leaves per bond.
#:
#: MEASURED 2026-08-19, and the fix for the flat-zig-zag artifact. Read with
#: :data:`STEERING_WIDTH_TARGET_CAP`, which bounds what this produces.
#:
#: **The defect.** A step's distance from residue 0 depends only on the angle between the new
#: bond and the radial direction, so holding that distance fixed puts the bond on a cone about
#: the radial axis; holding the pseudo-angle fixed puts it on a cone about the previous bond.
#: Two cones about different axes meet in **two points**. With the old flat 0.5 A width both
#: preferences are sharp at once, so each step was a choice between two mirror-image directions
#: either side of a plane that turns only slowly -- and the chain came out flat. It is not a
#: sampling artifact: re-gridding the candidates to resolve the radial shell (a generator laying
#: rings on the constant-distance cones, built and measured) changed nothing, because the target
#: distribution itself is bimodal in azimuth. The width had to change.
#:
#: **The measurement.** CA pseudo-dihedral planar order -- ``|<exp(2i * dihedral)>|``, which is
#: 1.0 for a perfect planar zig-zag and ~0.05-0.08 for a freely-rotating chain with DODO's own
#: angle distribution -- on p300, dnmt3a and arf19 at three seeds and two modes, per region:
#:
#:     region                          old     new     kind
#:     p300 1832-2414 (583 res)      0.618   0.044    terminal
#:     p300 1-334     (334 res)      0.578   0.053    terminal
#:     dnmt3a 1-282   (282 res)      0.604   0.051    terminal
#:     arf19 1042-1086 (45 res)      0.361   0.113    terminal
#:     every closure region          unchanged to 4 dp
#:
#: Closure regions are untouched by construction -- they steer toward a fixed anchor with
#: :func:`_natural_fluctuation`, which is tens of Angstroms wide and never had the defect. That
#: they come out bit-identical is the control on this change.
#:
#: **Why it scales with slack rather than being another constant.** The width is how much radial
#: freedom the walk is allowed, and how much it *can* be allowed is set by how much the target
#: leaves unused: a target near the geometric ceiling needs nearly every bond to advance
#: maximally and genuinely has no freedom to give, while a compact target has most of a bond
#: length spare per step. Swept over six seeds on the realistic length/extension curve, a flat
#: width fails at one end or the other -- 2.0 A flat misses by 17 A at 0.9 of reach, and the
#: ``sqrt(remaining)`` form the closure funnel uses collapses outright there (3 of 24 conformers
#: built). Scaling by the slack holds accuracy at both ends.
#:
#: Swept 0.4 / 0.6 / 0.8 (and 0.8 capped) at 180 conformers per cell. 0.6 is the knee: 0.4 leaves
#: the planar order at 0.12-0.13 on long regions, against 0.047-0.055 at 0.6 and a
#: freely-rotating reference of 0.047; 0.8 buys no further order and costs three times the
#: end-to-end scatter.
#:
#: **VALIDATED on the whole AlphaFold human proteome** -- all 23,587 structures, old and new width
#: over the *same* structures, at two independent build seeds (0 and 1):
#:
#:     metric                          old            new
#:     crashes / timeouts              0 / 0          0 / 0
#:     regions rebuilt                 99.894%        99.894%   (seed 0)
#:                                     99.905%        99.905%   (seed 1)
#:     blocked structures              49 / 43        49 / 43   -- the same structures, both seeds
#:     rebuilt bond defects            0              0
#:     introduced steric clashes       26,823         22,595    (seed 0)
#:     structures with >= 20 of them   12 -> 3 (seed 0), 20 -> 4 (seed 1)
#:     steered planar order            0.497          0.094
#:     closure planar order            0.086          0.086     -- the control, unchanged
#:     end-to-end, fraction > 5% off   0.01%          1.51%
#:     end-to-end, fraction > 10% off  0.00%          0.00%
#:
#: **Introduced impossible contacts are unchanged, and that took a control to establish.** At seed
#: 0 the count went 11 -> 16, which looks like a regression on the one invariant that is never a
#: ratchet. It is not: the old width alone gives 11 at seed 0 and 15 at seed 1, and at seed 1 the
#: two policies give **15 and 15**. The seed-0 gap sits inside the old code's own seed-to-seed
#: swing. Changing this width changes RNG consumption, so every build differs and the *set* of
#: affected structures churns freely; only the rate is meaningful, and the rate did not move.
#:
#: What did move, reproducibly at both seeds, is the *composition*: the median sequence separation
#: of an introduced impossible contact falls from ~92 to ~48 residues, and the fraction closer than
#: 20 residues rises from 15% to 35%. A freer walk doubles back on itself locally more often and
#: slams into distant structure less often. Both classes are equally impossible and the totals are
#: tiny (about 15 in 23,587), but if anyone chases these next, expect short-range ones now.
#:
#: Measure this the way ``known_failures.md`` defines it -- pairs in the output that are **not
#: already in the input**, matched on atom identity. Several AlphaFold inputs ship with impossible
#: pairs of their own (AF-Q9BTC0-F1 has 92), and some come out cleaner than they went in, so
#: counting output pairs conflates DODO's defects with AlphaFold's and inflates the number
#: threefold.
STEERING_SLACK_FRACTION: Final[float] = 0.6

#: Ceiling on the steering width as a fraction of that conformer's own end-to-end target.
#:
#: MEASURED 2026-08-19. The width is an absolute distance and its accuracy cost is a roughly
#: fixed number of Angstroms, so as a fraction of the target it falls hardest on the shortest
#: regions -- which are also the regions with nothing to gain, because a flat zig-zag needs
#: length to become visible. Measured planar order at 15 residues is 0.257 before the change and
#: 0.258 after: no benefit at all, for the largest accuracy bill of any length.
#:
#: This cap spends the width only where it does something. SWEPT over 500 AlphaFold human
#: structures (both policies on the same structures, one seed), by rebuilt region size --
#: mean planar order, and the fraction of steered regions finishing more than
#: :data:`~dodo.construct.pipeline._END_TO_END_NOTABLE_FRACTION` (5%) from their target:
#:
#:     residues     old       0.03      0.02      0.015     <- cap
#:     < 80       0.406      0.118     0.157      0.214
#:     < 120      0.450      0.089     0.112      0.163
#:     < 200      0.490      0.063     0.076      0.104
#:     < 400      0.567      0.050     0.050      0.052
#:     >= 400     0.616      0.047     0.046      0.045
#:     ---- fraction of steered regions more than 5% off target ----
#:     overall     0.00%     3.43%     0.86%      0.29%
#:     order > 0.25 rate
#:     overall    97.40%     1.10%     4.30%     11.10%
#:
#: 0.02 is the choice. Above 200 residues -- where a flat zig-zag is actually visible, and where
#: the complaint that started this came from -- it is indistinguishable from the uncapped 0.03,
#: and it introduces exactly as many steric clashes (518 either way over the 500 structures).
#: What it buys is a four-fold cut in how often the 5% reporting threshold fires. Tightening
#: further to 0.015 keeps giving accuracy back but starts costing real planar order on the
#: short regions, and 11% of regions above 0.25 is no longer a clean fix.
#:
#: **No region at any cap in any sweep exceeded** :data:`END_TO_END_TOLERANCE_FRACTION`, because
#: :func:`_validate_conformer` rejects a conformer that would and the region is retried. The 10%
#: contract is enforced, not merely aimed at; this cap is about how often the 5% *reporting*
#: threshold fires, not about correctness.
STEERING_WIDTH_TARGET_CAP: Final[float] = 0.02

#: Greatest distance one bond can add along a fixed direction, in Angstroms. DERIVED: the
#: axial advance of the all-trans zig-zag at :data:`~dodo.constants.BACKBONE_ANGLE_MAX`, the
#: same quantity :func:`max_reach` accumulates.
_MAX_ADVANCE_PER_BOND: Final[float] = CA_CA_BOND_LENGTH * math.sin(
    math.radians(BACKBONE_ANGLE_MAX) / 2.0
)

#: Numerical slack on the hard distance corridor, in Angstroms.
#:
#: Not a physical tolerance. The corridor walls are clamped to :func:`min_reach` and
#: :func:`max_reach`, which at the first step collapse onto exactly one bond length; without
#: a few ulps of slack the only geometrically legal move can fall outside a wall computed
#: from the same numbers by a different route.
_CORRIDOR_EPSILON: Final[float] = 1e-9

#: Allowed deviation of a caller-supplied outer residue from one CA-CA bond length, in A.
#: CHOICE, deliberately loose; see :func:`_check_fixed_context`.
_OUTER_BOND_SLACK: Final[float] = 1.0

# ---------------------------------------------------------------------------
# Reachability: how far a chain of k bonds can actually span
# ---------------------------------------------------------------------------


def max_reach(n_bonds: int, bond_length: float = CA_CA_BOND_LENGTH) -> float:
    """Greatest CA(i)-to-CA(i+n) distance spannable by ``n_bonds`` valid bonds.

    Parameters
    ----------
    n_bonds
        Number of CA-CA bonds available. ``0`` returns ``0.0``.
    bond_length
        CA-CA virtual bond length in Angstroms.

    Returns
    -------
    float
        Maximum end-to-end distance in Angstroms.

    Notes
    -----
    **Derivation.** The most extended chain compatible with a fixed CA-CA-CA pseudo-angle
    is the planar all-trans zig-zag at the largest permitted angle,
    :data:`~dodo.constants.BACKBONE_ANGLE_MAX`. Writing ``h`` for half that angle, each
    bond advances ``b * sin(h)`` along the chain axis and displaces ``b * cos(h)``
    transversely, with the transverse displacement alternating in sign. So after ``k``
    bonds the axial separation is ``k * b * sin(h)`` and the transverse separation is
    ``b * cos(h)`` for odd ``k`` and zero for even ``k``::

        max_reach(k) = hypot(k * b * sin(h), b * cos(h) if k odd else 0)

    which correctly returns exactly ``b`` at ``k = 1``.
    """
    if n_bonds < 0:
        raise ValueError(f"n_bonds must be non-negative, got {n_bonds}.")
    if n_bonds == 0:
        return 0.0
    half = math.radians(BACKBONE_ANGLE_MAX) / 2.0
    axial = n_bonds * bond_length * math.sin(half)
    transverse = bond_length * math.cos(half) if n_bonds % 2 else 0.0
    return float(math.hypot(axial, transverse))


def min_reach(n_bonds: int, bond_length: float = CA_CA_BOND_LENGTH) -> float:
    """Smallest CA(i)-to-CA(i+n) distance spannable by ``n_bonds`` valid bonds.

    Parameters
    ----------
    n_bonds
        Number of CA-CA bonds available.
    bond_length
        CA-CA virtual bond length in Angstroms.

    Returns
    -------
    float
        Minimum end-to-end distance in Angstroms.

    Notes
    -----
    Exact for the two cases that geometry pins, and the steric floor beyond them.

    ``k = 1`` is one bond, so the separation is the bond length. ``k = 2`` is fixed by the
    pseudo-angle at the intervening residue: ``2 * b * sin(theta / 2)``, minimized at
    :data:`~dodo.constants.BACKBONE_ANGLE_MIN`.

    For ``k >= 3`` the two atoms are more than
    :data:`~dodo.constants.CLASH_EXCLUDE_WITHIN_RESIDUES` apart along the chain, so they are
    a *non-bonded* pair and cannot approach closer than
    :data:`~dodo.constants.CA_CLASH_DISTANCE`. That is the bound returned. It is not a
    refinement of the old one, it replaces a wrong one: the previous version propagated the
    two-bond minimum backwards, ``two_bond_minimum - (k - 2) * b``, and hit exactly 0.0 by
    ``k = 4``. Zero is the claim that two distinct residues can occupy the same point, and
    it made every low-side feasibility check in this module dead code -- a 10-residue region
    asked for an end-to-end distance of 0.5 A built happily at 4 A, and two anchors at
    separation 0.0 A were not diagnosed at all.

    Still deliberately permissive: a chain of ``k`` bonds folding back to
    :data:`~dodo.constants.CA_CLASH_DISTANCE` also has to fit its intervening residues
    somewhere, so the true minimum for large ``k`` is above this. Over-rejecting is the more
    dangerous error, since it turns a buildable region into a reported failure, so the bound
    stays at the one thing that is certain.

    Note that ``bond_length`` therefore has no effect for ``k >= 3``: the steric floor is a
    property of two carbon atoms, not of the chain connecting them.
    """
    if n_bonds < 0:
        raise ValueError(f"n_bonds must be non-negative, got {n_bonds}.")
    if n_bonds == 0:
        return 0.0
    if n_bonds == 1:
        return float(bond_length)
    two_bond_minimum = 2.0 * bond_length * math.sin(math.radians(BACKBONE_ANGLE_MIN) / 2.0)
    if n_bonds == 2:
        return float(two_bond_minimum)
    return float(CA_CLASH_DISTANCE)


def reachability_tail(n_bonds: int, bond_length: float = CA_CA_BOND_LENGTH) -> np.ndarray:
    """Reachability bounds for the last ``n_bonds`` steps, longest first.

    Parameters
    ----------
    n_bonds
        Length of the tail, in bonds.
    bond_length
        CA-CA virtual bond length in Angstroms.

    Returns
    -------
    np.ndarray
        ``(n_bonds,)`` maximum spannable distances for ``k = n_bonds, ..., 1``.

    Notes
    -----
    The derived replacement for the hardcoded table described in :func:`max_reach`. The
    walk consults :func:`max_reach` directly and does not need this array; it is exposed so
    that a caller -- or a test comparing against the old literals -- can see the schedule.
    """
    return np.array([max_reach(k, bond_length) for k in range(n_bonds, 0, -1)], dtype=np.float64)


def reachable_envelope(
    n_residues: int,
    n_anchor_xyz: np.ndarray | None,
    c_anchor_xyz: np.ndarray | None,
) -> float | None:
    """Bound on how far a region's alpha carbons can stray from its anchors.

    Parameters
    ----------
    n_residues
        Residues in the region.
    n_anchor_xyz
        CA of the residue the region attaches to at its N-terminal end, or ``None``.
    c_anchor_xyz
        CA of the residue the region attaches to at its C-terminal end, or ``None``.

    Returns
    -------
    float or None
        With both anchors, the bound ``L`` on ``|x - N| + |x - C|`` -- the prolate spheroid
        with the anchors as foci that contains every position the region can occupy. With one
        anchor, the bound on ``|x - anchor|``, a sphere. ``None`` when there is no anchor at
        all, where the region is built in its own frame and nothing constrains it.

    Notes
    -----
    DERIVED from :func:`max_reach` alone, so it is as conservative as that bound is. Residue
    ``i`` (0-based, growth order) is ``i + 1`` bonds from the start anchor and, in closure mode,
    ``n_residues - i`` bonds from the far anchor. So
    ``|x - N| <= max_reach(i + 1)`` and ``|x - C| <= max_reach(n_residues - i)``, and the sum is
    bounded by the largest such sum over ``i``. Taking the maximum rather than evaluating at one
    ``i`` matters because :func:`max_reach` is not quite linear: it carries a transverse
    ``b * cos(angle / 2)`` term for odd bond counts.

    Being an *ellipsoid* rather than two spheres is what makes the two-anchor case worth having:
    the intersection of the two spheres is larger than the set actually reachable, because a
    residue cannot be far from both anchors at once when the bonds have to be shared between the
    two halves of the walk.
    """
    if n_residues < 1:
        raise ValueError(f"n_residues must be at least 1, got {n_residues}.")
    if n_anchor_xyz is not None and c_anchor_xyz is not None:
        return max(max_reach(i + 1) + max_reach(n_residues - i) for i in range(n_residues))
    if n_anchor_xyz is None and c_anchor_xyz is None:
        return None
    return max_reach(n_residues)


# ---------------------------------------------------------------------------
# Per-step goal and whole-region plan
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class _Goal:
    """A corridor and a preferred value for one step's distance to a reference point.

    ``lo``/``hi`` are hard: a candidate outside them is provably unable to reach the
    region's endpoint, so it is dropped. The rest is soft, weighting the random choice
    among the survivors.

    Every field is a column, ``(rows, 1)``, so that it broadcasts against the
    ``(n_live, n_candidates)`` distance array. ``rows`` is 1 when the whole batch shares a
    goal -- closure mode, where the anchors dictate the schedule -- and ``n_live`` when each
    conformer is steered to its own end-to-end target, which is what makes a batch an
    ensemble rather than a set of copies.

    In closure mode the two widths are deliberately different, and the asymmetry is the
    single most important thing about that funnel. Falling *behind* the schedule is cheap:
    the chain still has slack and can catch up. Running *ahead* of it is not recoverable,
    because a chain can close the distance to its target by at most one bond length per
    residue while the reachability ceiling drops by nearly a whole bond length per residue
    too -- so a walk that drifts up against the ceiling has no legal move left and dies
    several residues later. With a symmetric width, a taut bridge (extension fraction 0.8
    and up) dead-ends around two thirds of the way through essentially every time.

    For a free or terminal region the widths are deliberately *equal*, and that is a fix
    rather than a simplification. A one-sided width plus a hard low wall meant the walk
    lagged its schedule all the way along and was dragged to the endpoint by the wall,
    landing on it: measured mean(Re)/target was 0.957-0.965 across every expansion mode and
    length, with all 60 conformers of a 200-residue region inside 1.9 A of ``lo``. A
    symmetric width around an unbiased schedule lets the achieved distance scatter about the
    target instead of pressing against one side of it.
    """

    lo: np.ndarray
    hi: np.ndarray
    want: np.ndarray
    #: Width applied to candidates closer to the endpoint than the schedule wants.
    sigma_down: np.ndarray
    #: Width applied to candidates further away than the schedule wants.
    sigma_up: np.ndarray

    def log_weight(self, distance: np.ndarray) -> np.ndarray:
        """Log relative preference for candidates at the given distances from the reference.

        Returned as a logarithm, and that is not a micro-optimization. A 200-residue region
        steered toward 700 A starts 690 A from where it needs to end up, and
        ``exp(-0.5 * (690 / 30) ** 2)`` is 1e-115 -- it underflows to exactly 0.0 for *every*
        candidate. The selection then falls back to a uniform draw, so the steering silently
        switches itself off in precisely the situation that needs it most, the walk fails to
        extend, and the only thing left holding it to the target is the hard corridor wall.
        Working in logs and shifting each row by its own maximum keeps the relative
        preferences intact at any distance.
        """
        width = np.where(distance > self.want, self.sigma_up, self.sigma_down)
        scaled = (distance - self.want) / width
        log_weight: np.ndarray = -0.5 * scaled * scaled
        return log_weight


def _column(value: float | np.ndarray) -> np.ndarray:
    """Reshape a scalar or 1-D per-conformer array into a broadcastable ``(rows, 1)`` column."""
    array = np.atleast_1d(np.asarray(value, dtype=np.float64))
    return array.reshape(-1, 1)


def sample_end_to_end_targets(
    target: float,
    n_conformations: int,
    rng: np.random.Generator,
    *,
    low: float,
    high: float,
) -> np.ndarray:
    """Draw one end-to-end target per conformer, centred on the requested mean.

    Parameters
    ----------
    target
        Requested *mean* end-to-end distance in Angstroms, from
        :attr:`dodo.engines.base.IDRRequest.target_end_to_end`.
    n_conformations
        Number of conformers in the batch.
    rng
        Seeded generator. The only source of randomness here.
    low, high
        Feasible range for a single conformer, in Angstroms: the steric floor from
        :func:`min_reach` and :data:`MAX_TARGET_EXTENSION_FRACTION` times the geometric
        reach from :func:`max_reach`.

    Returns
    -------
    np.ndarray
        ``(n_conformations,)`` targets, inside ``[low, high]``, whose mean is ``target``.

    Notes
    -----
    **Why this exists.** ALBATROSS predicts a mean, and the code being fixed treated it as a
    hard per-conformer constraint. That collapses the ensemble: measured CV(Re) over batches
    of 30 was 0.006-0.045 for every mode and length, against 0.35-0.48 for a
    freely-rotating chain with the same bond length and the same CA-CA-CA angle window.
    Sixty models of a 200-residue IDR spanning 1.9 A of extension is one conformation
    written sixty times, which is worse than useless for a pseudo-trajectory: it looks like
    an ensemble and reports like one.

    **The distribution.** For a Gaussian chain the end-to-end *vector* is isotropic
    Gaussian, so the end-to-end *distance* follows
    ``P(R) proportional to R ** 2 * exp(-3 * R ** 2 / (2 * <R ** 2>))``. Sampling it is one
    line -- the norm of three independent Gaussians -- and gives a coefficient of variation
    of 0.42, squarely inside the measured polymer range.

    **The mean is pinned, the spread is not.** The draw is rescaled so that the *sample*
    mean is exactly ``target`` before clipping, for two reasons. It removes sampling noise
    in the one quantity the caller actually asked for, and it makes ``n_conformations=1``
    behave exactly as before this change: a single conformer is steered at the requested
    value rather than at a random draw from around it. Rescaling is multiplicative, so the
    coefficient of variation -- the part that makes it an ensemble -- is untouched.

    Clipping to the feasible range can pull the mean off ``target`` again, so it is followed
    by a few rounds of redistributing the shortfall over the samples that are still free to
    move. Where the request sits hard against the geometric ceiling this cannot fully
    converge, and the ensemble is then legitimately narrow: an ensemble with mean 0.96 of
    full extension and a 42% spread does not exist.
    """
    if not math.isfinite(low) or not math.isfinite(high) or low > high:
        raise EngineError(
            f"Infeasible target range [{low}, {high}] for the conformer ensemble; this is a "
            f"DODO bug, the plan should have raised UnsatisfiableTargetError first."
        )
    scale = target / MAXWELL_MEAN_FACTOR
    # Three independent Gaussians of variance <R^2>/3 each: the norm is Maxwell-Boltzmann
    # distributed with the intended <R^2>. Drawn from the caller's generator, so a seed
    # reproduces the whole ensemble.
    deviates = rng.normal(0.0, scale / math.sqrt(3.0), size=(n_conformations, 3))
    values = np.linalg.norm(deviates, axis=1)
    mean = float(values.mean())
    if mean > 0.0:
        values = values * (target / mean)

    clipped = np.clip(values, low, high)
    for _ in range(16):
        shortfall = (target - float(clipped.mean())) * n_conformations
        if abs(shortfall) <= 1e-12 * max(target, 1.0) * n_conformations:
            break
        movable = (clipped > low) & (clipped < high) if shortfall < 0.0 else clipped < high
        if not np.any(movable):
            break
        share = clipped[movable] / float(clipped[movable].sum())
        clipped[movable] = np.clip(clipped[movable] + shortfall * share, low, high)
    result: np.ndarray = clipped
    return result


@dataclass(frozen=True, slots=True)
class _WalkPlan:
    """Everything the growth loop needs, resolved once and checked for feasibility.

    Growth always runs *outward from a fixed anchor*, so an N-terminal tail (no N-anchor,
    but a C-anchor to attach to) is grown from its C-anchor backwards and the result is
    flipped into N-to-C order at the end. That is why this holds ``start_anchor`` and
    ``end_anchor`` rather than N and C: one code path covers all four cases.
    """

    n_residues: int
    start_anchor: np.ndarray | None
    end_anchor: np.ndarray | None
    reverse: bool
    target: float
    tolerance: float
    #: Fraction of the available geometric reach the chain must span. Seeded from the real
    #: anchor separation in closure mode and from the target end-to-end otherwise.
    extension_fraction: float
    #: Fixed CA one residue beyond ``start_anchor``, away from the region, or ``None``.
    #: Present only to constrain and measure the pseudo-angle centred on that anchor.
    start_outer: np.ndarray | None = None
    #: Fixed CA one residue beyond ``end_anchor``, away from the region, or ``None``.
    end_outer: np.ndarray | None = None
    #: Per-conformer end-to-end target in Angstroms, ``(n_conformations,)``. All equal to
    #: ``target`` in closure mode, where the anchors dictate the span; drawn from
    #: :func:`sample_end_to_end_targets` otherwise, so that a batch is an ensemble.
    targets: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.float64))
    #: Fractional tolerance on the achieved end-to-end distance, applied per conformer to
    #: that conformer's own target.
    tolerance_fraction: float = END_TO_END_TOLERANCE_FRACTION
    #: The end-to-end distance the *region* was asked for, before any per-conformer spreading.
    #: Equals :attr:`target` unless the caller did the spreading itself; see
    #: :attr:`dodo.engines.base.IDRRequest.ensemble_mean_end_to_end`.
    mean_target: float = 0.0

    @property
    def closes(self) -> bool:
        """True if the walk must land one bond length from a fixed far anchor."""
        return self.end_anchor is not None

    @property
    def n_grown(self) -> int:
        """Residues placed by seeded growth; the rest is the analytic closure step."""
        return self.n_residues - 1 if self.closes else self.n_residues

    @classmethod
    def build(
        cls, request: IDRRequest, tolerance_fraction: float, rng: np.random.Generator
    ) -> _WalkPlan:
        """Resolve a request into a plan, raising if it cannot be satisfied at all.

        Parameters
        ----------
        request
            The region to build.
        tolerance_fraction
            Fractional tolerance on the achieved end-to-end distance.
        rng
            Seeded generator, used to draw the per-conformer end-to-end targets. Taken here
            rather than at growth time so that the ensemble is fixed before any sampling and
            a retry cannot quietly re-roll a conformer an easier target.

        Returns
        -------
        _WalkPlan
            The resolved plan.

        Raises
        ------
        GeometryError
            If the two anchors are closer together or further apart than any valid chain
            of this length can span, or if the fixed residues the caller supplied are
            themselves not valid backbone geometry.
        UnsatisfiableTargetError
            If the target end-to-end distance is outside what this many residues can span.
        """
        n = request.n_residues
        # Pure fraction, no floor: see END_TO_END_TOLERANCE_FRACTION for why the old
        # CA_CA_BOND_LENGTH floor made every low-side check in this module unreachable.
        tolerance = tolerance_fraction * request.target_end_to_end
        # Falls back to the per-conformer target, which is correct whenever the engine is the
        # one doing the spreading -- then they are the same number by construction.
        mean_target = (
            request.target_end_to_end
            if request.ensemble_mean_end_to_end is None
            else float(request.ensemble_mean_end_to_end)
        )
        _check_fixed_context(request)

        if request.n_anchor_xyz is not None and request.c_anchor_xyz is not None:
            separation = float(np.linalg.norm(request.c_anchor_xyz - request.n_anchor_xyz))
            n_bonds = n + 1  # anchor -> r0 -> ... -> r(n-1) -> anchor
            ceiling = max_reach(n_bonds)
            floor = min_reach(n_bonds)
            if separation > ceiling:
                raise GeometryError(
                    f"Cannot bridge anchors {separation:.2f} A apart with {n} residue(s): "
                    f"{n_bonds} CA-CA bonds constrained to pseudo-angles of at most "
                    f"{BACKBONE_ANGLE_MAX:.0f} degrees span at most {ceiling:.2f} A. The "
                    f"region needs more residues, or the flanking domains are too far apart."
                )
            if separation < floor:
                raise GeometryError(
                    f"Cannot bridge anchors only {separation:.2f} A apart with {n} "
                    f"residue(s): {n_bonds} bonds cannot double back that tightly without "
                    f"a pseudo-angle below {BACKBONE_ANGLE_MIN:.0f} degrees or two "
                    f"non-bonded CA atoms closer than {CA_CLASH_DISTANCE:.2f} A (minimum "
                    f"span {floor:.2f} A)."
                )
            return cls(
                n_residues=n,
                start_anchor=request.n_anchor_xyz,
                end_anchor=request.c_anchor_xyz,
                reverse=False,
                target=request.target_end_to_end,
                tolerance=tolerance,
                extension_fraction=separation / ceiling,
                start_outer=request.n_anchor_prev_xyz,
                end_outer=request.c_anchor_next_xyz,
                # The anchors dictate the span, so there is no ensemble of targets to draw:
                # the variety between conformers comes from the walk between the anchors.
                targets=np.full(request.n_conformations, request.target_end_to_end),
                tolerance_fraction=tolerance_fraction,
                mean_target=mean_target,
            )

        # No closure: the achieved end-to-end distance is between the first and last
        # generated residue, which is n - 1 bonds apart.
        span_bonds = n - 1
        fraction = 0.0
        targets = np.full(request.n_conformations, request.target_end_to_end)
        if span_bonds >= 1:
            ceiling = max_reach(span_bonds)
            floor = min_reach(span_bonds)
            if request.target_end_to_end > ceiling + tolerance:
                raise UnsatisfiableTargetError(
                    f"A {n}-residue region cannot reach the requested end-to-end distance.",
                    target=request.target_end_to_end,
                    achievable=ceiling,
                )
            if request.target_end_to_end < floor - tolerance:
                raise UnsatisfiableTargetError(
                    f"A {n}-residue region cannot be that compact: with {span_bonds} bond(s) "
                    f"the first and last CA are more than "
                    f"{CLASH_EXCLUDE_WITHIN_RESIDUES} residues apart, so they cannot "
                    f"approach closer than {floor:.2f} A without overlapping.",
                    target=request.target_end_to_end,
                    achievable=floor,
                )
            fraction = min(1.0, request.target_end_to_end / ceiling)
            targets = sample_end_to_end_targets(
                request.target_end_to_end,
                request.n_conformations,
                rng,
                low=floor,
                high=max(floor, MAX_TARGET_EXTENSION_FRACTION * ceiling),
            )

        reverse = request.n_anchor_xyz is None and request.c_anchor_xyz is not None
        start = request.c_anchor_xyz if reverse else request.n_anchor_xyz
        # Growth runs outward from whichever anchor exists, so for an N-terminal tail the
        # "outer" residue of the growth start is the one C-terminal to the C-anchor.
        start_outer = request.c_anchor_next_xyz if reverse else request.n_anchor_prev_xyz
        return cls(
            n_residues=n,
            start_anchor=start,
            end_anchor=None,
            reverse=reverse,
            target=request.target_end_to_end,
            tolerance=tolerance,
            extension_fraction=fraction,
            start_outer=start_outer,
            end_outer=None,
            targets=targets,
            tolerance_fraction=tolerance_fraction,
            mean_target=mean_target,
        )

    def tolerance_for(self, target: float | np.ndarray) -> np.ndarray:
        """Absolute tolerance on the achieved end-to-end distance of ``target``."""
        return _column(self.tolerance_fraction * np.asarray(target, dtype=np.float64))

    def goal_for(
        self,
        index: int,
        targets: float | np.ndarray | None = None,
        current: np.ndarray | None = None,
    ) -> _Goal | None:
        """Return the distance corridor and preference for residue ``index``, if any.

        Parameters
        ----------
        index
            Index of the residue about to be placed, in growth order.
        targets
            End-to-end target for each live conformer, ``(n_live,)``, or ``None`` to use the
            request's single value. Ignored in closure mode, where the anchors decide.
        current
            Distance of the *previously* placed residue from residue 0, ``(n_live,)``, or
            ``None``. Supplying it makes the schedule self-correcting; see the notes.

        Returns
        -------
        _Goal or None
            ``None`` when the step is unconstrained -- the very first residue of a free or
            terminal region, which can go anywhere.

        Notes
        -----
        Why ``current`` matters, given that a precomputed profile would be simpler. Each step
        lands slightly short of what the schedule asked for -- the candidate set is discrete
        and the angle weighting pulls against extension -- and a schedule computed once from
        ``index`` alone lets that shortfall accumulate over every bond. Measured: about 0.2 A
        per step, which over 99 bonds is the entire 5% one-sided miss the achieved end-to-end
        distance used to show. Recomputing the demand from where the walk actually *is*
        redistributes the residual over the bonds that remain, so a shortfall is repaid
        instead of compounding, and the final miss is one step's worth rather than the sum.
        """
        if self.closes:
            # Bonds still needed to get from this residue to the far anchor: this residue
            # to the last generated one, plus the final closure bond.
            remaining = self.n_residues - index
            ceiling = max_reach(remaining)
            floor = min_reach(remaining)
            # Constant extension fraction: the schedule is the same proportion of the
            # available reach at every step, so the relative margin never shrinks. Seeded
            # from the real anchor separation, which is regression (b).
            want = min(max(self.extension_fraction * ceiling, floor), ceiling)
            return _Goal(
                lo=_column(floor),
                hi=_column(want + SCHEDULE_MARGIN_FRACTION * (ceiling - want)),
                want=_column(want),
                sigma_down=_column(_natural_fluctuation(remaining)),
                sigma_up=_column(_headroom_width(ceiling - want)),
            )

        if index < 1 or self.n_residues < 2:
            # Distance is measured from residue 0, so residue 0 has no goal, and a
            # one-residue region has no end-to-end distance to steer toward.
            return None

        span = self.n_residues - 1
        remaining = span - index
        target = _column(self.target if targets is None else targets)
        tolerance = self.tolerance_for(target)
        reach_left = max_reach(remaining)

        # Hard corridor: the endpoint still has to land ``target`` from residue 0 using
        # ``remaining`` bonds, so a candidate at distance d is usable only while
        # ``|target - d| <= reach_left + tolerance``.
        band_lo = min_reach(index) - _CORRIDOR_EPSILON
        band_hi = max_reach(index) + _CORRIDOR_EPSILON
        lo = np.clip(target - reach_left - tolerance, band_lo, band_hi)
        hi = np.clip(target + reach_left + tolerance, band_lo, band_hi)

        # Soft preference: the mean radial profile of a Gaussian bridge. For an ideal chain of
        # ``span`` bonds pinned at the origin and at distance R, residue ``index`` sits
        # ``R * index / span`` along the end-to-end axis with transverse variance
        # ``b ** 2 * index * remaining / span``, so its expected distance from residue 0 is the
        # hypotenuse of the two. That shape matters in both limits: it grows like ``sqrt``
        # early for a slack target, so a compact region is a coil rather than a spoke, and it
        # becomes linear in ``index`` for a taut one, which is the only profile that reaches an
        # extended target at all.
        axial = target * (index / span)
        transverse = CA_CA_BOND_LENGTH * math.sqrt(index * remaining / span)
        schedule = np.hypot(axial, transverse)

        if current is not None:
            # Self-correcting demand: cover the residual that is actually left, spread evenly
            # over the bonds that are actually left. ``remaining + 1`` counts this step.
            here = _column(current)
            catch_up = here + (target - here) / (remaining + 1)
            # The larger of the two, not the average: the bridge profile is what makes a slack
            # target a coil rather than a spoke (mid-chain it deliberately bulges *past* a very
            # compact target and comes back), while the catch-up term is what stops a taut
            # target from being missed.
            schedule = np.maximum(schedule, catch_up)
        schedule = np.clip(schedule, band_lo, band_hi)

        width = _target_steering_width(target, span, self.mean_target)
        return _Goal(
            lo=lo,
            hi=hi,
            want=schedule,
            # Symmetric: see the class docstring.
            sigma_down=width,
            sigma_up=width,
        )


# ---------------------------------------------------------------------------
# The engine
# ---------------------------------------------------------------------------


@dataclass
class SelfAvoidingWalk:
    """Grow a CA trace by cone-constrained, clash-avoiding, closure-aware steps.

    Attributes
    ----------
    name
        Engine identifier, recorded on every result.
    batch_size
        Conformers grown per batch. Bounds peak memory and is the unit over which the
        obstacle clash query is vectorized. Defaults to
        :data:`dodo.constants.WALK_BATCH_SIZE`.
    max_attempts
        Attempt rounds before giving up on the conformers still outstanding. Defaults to
        :data:`dodo.constants.MAX_ATTEMPTS_PER_REGION`.
    tolerance_fraction
        Fractional tolerance on the achieved end-to-end distance of a free or terminal
        region. See :data:`END_TO_END_TOLERANCE_FRACTION`.
    require_all
        Raise unless every requested conformer was built. Default True, because a caller
        asking for ten models and silently receiving seven is the kind of surprise this
        rewrite exists to remove. Set False to accept a partial batch, which is then
        reported through :attr:`~dodo.engines.base.IDRResult.success` with the failed
        rows NaN-filled.

    Notes
    -----
    Reproducibility. Every random decision comes from the ``rng`` passed to
    :meth:`generate`; two calls with equally seeded generators produce bit-identical
    output. Retries do not reseed -- they simply keep drawing from the caller's generator,
    which is both "fresh randomness" and reproducible.
    """

    name: str = ENGINE_NAME
    batch_size: int = WALK_BATCH_SIZE
    max_attempts: int = MAX_ATTEMPTS_PER_REGION
    tolerance_fraction: float = END_TO_END_TOLERANCE_FRACTION
    require_all: bool = True

    def __post_init__(self) -> None:
        """Validate the engine's own configuration."""
        if self.batch_size < 1:
            raise EngineError(f"batch_size must be at least 1, got {self.batch_size}.")
        if self.max_attempts < 1:
            raise EngineError(f"max_attempts must be at least 1, got {self.max_attempts}.")
        if not 0.0 < self.tolerance_fraction <= 1.0:
            raise EngineError(
                f"tolerance_fraction must be in (0, 1], got {self.tolerance_fraction}."
            )

    def available(self) -> bool:
        """Return True: this engine needs nothing beyond numpy and scipy.

        That is the point of it. It is what keeps ``pip install dodo`` -- with no torch,
        no model weights and no network -- able to rebuild a structure.
        """
        return True

    # ------------------------------------------------------------------
    # Public entry point
    # ------------------------------------------------------------------

    def generate(
        self,
        request: IDRRequest,
        obstacles: np.ndarray | None,
        rng: np.random.Generator,
    ) -> IDRResult:
        """Generate conformers for one disordered region.

        Parameters
        ----------
        request
            The region to build. For an interior region -- both anchors fixed -- the
            achieved end-to-end distance is dictated by the anchor separation, and
            :attr:`~dodo.engines.base.IDRRequest.target_end_to_end` is used only for
            reporting: the anchors are already placed and they win. For a free or terminal
            region it is the *mean* of the batch, not a per-conformer constraint; see
            :func:`sample_end_to_end_targets`.

            Supply :attr:`~dodo.engines.base.IDRRequest.n_anchor_prev_xyz` and
            :attr:`~dodo.engines.base.IDRRequest.c_anchor_next_xyz` wherever the anchors
            have chain neighbours. Without them the pseudo-angle centred on each anchor can
            neither be steered nor measured, and the result says so through
            :attr:`~dodo.engines.base.IDRResult.unconstrained_junctions` and an
            :class:`UnconstrainedJunctionWarning`.
        obstacles
            ``(n_obstacles, 3)`` coordinates to avoid, or ``None``. Should exclude the
            atoms of the anchor residues themselves; see
            :class:`~dodo.engines.base.IDRRequest`. These are the only positions the clash
            threshold is ever relaxed against: the region's own trace, its anchors and its
            flanking residues are held to :data:`~dodo.constants.CA_CLASH_DISTANCE`
            unconditionally.
        rng
            Seeded generator. Required, positionally or by keyword.

        Returns
        -------
        IDRResult
            Conformers and the explicit success mask.

        Raises
        ------
        GeometryError
            If the anchors cannot be bridged by this many residues, if the fixed residues
            handed in are not valid geometry, or if the obstacles block the only available
            starting point.
        UnsatisfiableTargetError
            If the target end-to-end distance is beyond what this many residues can span,
            or below what they can contract to without two alpha carbons overlapping.
        ExhaustedAttemptsError
            If no conformer -- or, with ``require_all``, not every conformer -- could be
            built within :attr:`max_attempts` rounds.

        Warns
        -----
        UnconstrainedJunctionWarning
            If an anchor has no neighbour coordinate, so its junction pseudo-angle is
            neither constrained nor checked.
        """
        if not isinstance(rng, np.random.Generator):
            raise TypeError(
                f"rng must be a numpy.random.Generator, got {type(rng).__name__}. Build one "
                f"with numpy.random.default_rng(seed): every stochastic step in DODO is "
                f"seedable so that any structure can be reproduced exactly."
            )
        plan = _WalkPlan.build(request, self.tolerance_fraction, rng)
        tree = self._obstacle_tree(self._prune_unreachable(request, obstacles))
        self._check_start_is_free(plan, tree)
        loose_junctions = request.unconstrained_junctions
        if loose_junctions:
            warnings.warn(
                f"The junction pseudo-angle centred on {' and '.join(loose_junctions)} cannot be "
                f"constrained or measured: the CA of the residue beyond the anchor was not "
                f"supplied. That angle belongs to the assembled chain, and measured over 480 "
                f"such junctions 41% of them fell below {BACKBONE_ANGLE_MIN:.0f} degrees. Pass "
                f"n_anchor_prev_xyz / c_anchor_next_xyz to close it.",
                UnconstrainedJunctionWarning,
                stacklevel=2,
            )

        n_conformations = request.n_conformations
        coords = np.full((n_conformations, plan.n_residues, 3), np.nan, dtype=np.float64)
        success = np.zeros(n_conformations, dtype=bool)
        relaxed_to: float | None = None
        attempts = 0
        last_failure = "no attempt was made"

        pending = np.flatnonzero(~success)
        while pending.size and attempts < self.max_attempts:
            attempts += 1
            for offset in range(0, pending.size, self.batch_size):
                chunk = pending[offset : offset + self.batch_size]
                # Each conformer keeps its own end-to-end target across retries. Re-drawing
                # it would let a conformer that failed a demanding target quietly get an
                # easier one, which is both irreproducible and a silent miss of the request.
                built, ok, rungs, reason = self._grow_batch(plan, plan.targets[chunk], tree, rng)
                if reason is not None:
                    last_failure = reason
                if not np.any(ok):
                    continue
                coords[chunk[ok]] = built[ok]
                success[chunk[ok]] = True
                worst = float(np.min(rungs[ok]))
                if worst < CA_CLASH_DISTANCE:
                    relaxed_to = worst if relaxed_to is None else min(relaxed_to, worst)
            pending = np.flatnonzero(~success)

        if not np.any(success):
            raise ExhaustedAttemptsError(
                f"The {self.name} engine built none of the {n_conformations} requested "
                f"conformer(s) for a {plan.n_residues}-residue region. Last failure: "
                f"{last_failure}.",
                attempts=attempts,
            )
        if self.require_all and not np.all(success):
            raise ExhaustedAttemptsError(
                f"The {self.name} engine built only {int(np.count_nonzero(success))} of "
                f"{n_conformations} requested conformer(s) for a {plan.n_residues}-residue "
                f"region. Last failure: {last_failure}. Pass require_all=False to accept a "
                f"partial batch, which is reported through IDRResult.success.",
                attempts=attempts,
            )

        if plan.reverse:
            # Grown outward from the C-anchor, so growth order is C-to-N. Flip into
            # sequence order; the caller indexes these by residue and must not have to
            # know which end the walk started from.
            coords = coords[:, ::-1, :]
        return IDRResult.from_batch(
            coords, success, self.name, attempts, relaxed_to, loose_junctions
        )

    # ------------------------------------------------------------------
    # Growth
    # ------------------------------------------------------------------

    def _grow_batch(
        self,
        plan: _WalkPlan,
        targets: np.ndarray,
        obstacle_tree: cKDTree | None,
        rng: np.random.Generator,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, str | None]:
        """Grow one conformer per entry of ``targets``, in lockstep.

        Parameters
        ----------
        plan
            The resolved plan.
        targets
            ``(size,)`` end-to-end target for each conformer in this batch. Unused in
            closure mode, where the anchors set the span.
        obstacle_tree
            Spatial index over the obstacle coordinates, or ``None``.
        rng
            Seeded generator.

        Returns
        -------
        tuple
            ``(coords, ok, rungs, reason)``: ``(size, n_residues, 3)`` coordinates,
            ``(size,)`` success mask, ``(size,)`` loosest clash threshold each conformer
            needed, and a description of the first failure seen, or ``None``.

        Notes
        -----
        Conformers advance one residue at a time together, which is what lets the obstacle
        clash test for the whole batch -- up to ``size * _candidates_per_step()`` positions -- be a
        single KD-tree query. A conformer that dead-ends drops out of the batch
        immediately and its row stays NaN.
        """
        n = plan.n_residues
        size = int(targets.size)
        coords = np.full((size, n, 3), np.nan, dtype=np.float64)
        alive = np.ones(size, dtype=bool)
        rungs = np.full(size, CA_CLASH_DISTANCE, dtype=np.float64)
        reason: str | None = None

        first = 0
        if plan.start_anchor is None:
            # A free region is built in its own frame: put residue 0 at the origin. The
            # orchestrator positions the region afterwards.
            coords[:, 0] = 0.0
            first = 1

        for index in range(first, plan.n_grown):
            live = np.flatnonzero(alive)
            if live.size == 0:
                break
            candidates, weights, centres = self._candidates_for(plan, coords, live, index, rng)
            chosen, ok, rung = self._select(
                candidates,
                weights=weights,
                extra_mask=None,
                reference=self._reference_for(plan, coords, live),
                goal=plan.goal_for(
                    index, targets[live], self._distance_from_origin(plan, coords, live, index)
                ),
                chain_points=self._chain_points(plan, coords, live, index),
                obstacle_tree=obstacle_tree,
                candidate_centres=centres,
                rng=rng,
            )
            coords[live[ok], index] = chosen[ok]
            rungs[live[ok]] = np.minimum(rungs[live[ok]], rung[ok])
            alive[live[~ok]] = False
            if not np.all(ok) and reason is None:
                reason = (
                    f"no candidate for residue {index} of {n} survived the clash, angle "
                    f"and reachability filters"
                )

        if plan.closes and np.any(alive):
            ok_closure, closure_reason = self._close(plan, coords, alive, rungs, obstacle_tree, rng)
            alive &= ok_closure
            if closure_reason is not None and reason is None:
                reason = closure_reason

        # Final independent check. Everything above is meant to be correct by
        # construction, so this should never reject anything -- which is exactly why it is
        # worth running: it is the difference between "success means we believe it worked"
        # and "success means the geometry was measured and is valid".
        for row in np.flatnonzero(alive):
            clearance = float("inf")
            if obstacle_tree is not None:
                found, _ = obstacle_tree.query(coords[row])
                clearance = float(np.min(found))
            problem = _validate_conformer(
                coords[row],
                plan,
                float(targets[row]),
                clearance,
                float(rungs[row]),
            )
            if problem is not None:
                alive[row] = False
                if reason is None:
                    reason = f"generated conformer failed validation: {problem}"
        return coords, alive, rungs, reason

    def _candidates_for(
        self,
        plan: _WalkPlan,
        coords: np.ndarray,
        live: np.ndarray,
        index: int,
        rng: np.random.Generator,
    ) -> tuple[np.ndarray, np.ndarray | None, np.ndarray]:
        """Return candidate positions for one residue, their angle weights, and their centre.

        Parameters
        ----------
        plan
            The resolved plan.
        coords
            ``(size, n_residues, 3)`` working coordinates.
        live
            Indices of the conformers still growing.
        index
            Residue index about to be placed, in growth order.
        rng
            Seeded generator.

        Returns
        -------
        tuple
            ``(candidates, weights, centres)`` with candidates ``(n_live, n_candidates, 3)``,
            weights ``(n_candidates,)`` or ``None`` for an unweighted step, and ``(n_live, 3)``
            centres.

        Notes
        -----
        ``centres`` is the point every candidate in a row sits *exactly* one
        :data:`~dodo.constants.CA_CA_BOND_LENGTH` from -- the sphere centre in the
        unconstrained-junction case below, the cone apex otherwise. It is returned rather than
        recomputed by the caller because which point that is depends on the same three-way branch
        as the candidates themselves, and a copy of that branch elsewhere would be free to drift
        out of step with this one. Getting it wrong would silently disarm the clash query cull in
        :meth:`_nearest_obstacle_distance`.
        """
        n_live = int(live.size)
        count = _candidates_per_step()

        # Two prior points are needed to define a cone. The first generated residue has only
        # one -- the anchor -- unless the caller supplied the residue *beyond* the anchor, in
        # which case the pair (outer, anchor) defines the cone and the pseudo-angle centred
        # on the anchor is in-window by construction, exactly as for an interior residue.
        # Without it there is nothing to form an angle with, so the 3.8 A sphere is sampled
        # and the junction angle of the assembled chain is whatever it happens to be: that is
        # the 38.6-176.8 degree measurement in the module docstring, and it is why
        # generate() warns rather than proceeding quietly.
        if (index == 0 and plan.start_outer is None) or (index == 1 and plan.start_anchor is None):
            centre = (
                np.broadcast_to(_require_anchor(plan.start_anchor), (n_live, 3))
                if index == 0
                else coords[live, 0]
            )
            directions = random_unit_vectors(n_live * count, rng).reshape(n_live, count, 3)
            candidates: np.ndarray = centre[:, None, :] + CA_CA_BOND_LENGTH * directions
            return candidates, None, np.asarray(centre, dtype=np.float64)

        grid = backbone_angle_grid()
        per_angle = _per_angle()
        if index == 0:
            # Junction step: the cone apex is the anchor and the axis runs from the fixed
            # residue beyond it, so the angle this constrains is the one centred on the
            # anchor itself.
            previous = np.broadcast_to(_require_anchor(plan.start_anchor), (n_live, 3))
            before = np.broadcast_to(_require_anchor(plan.start_outer), (n_live, 3))
        else:
            previous = coords[live, index - 1]
            before = (
                np.broadcast_to(_require_anchor(plan.start_anchor), (n_live, 3))
                if index == 1
                else coords[live, index - 2]
            )
        # One batched rotation for the whole live set: cone_candidates_batch rotates the
        # shared cached template onto every conformer's axis at once, and its output is
        # bit-identical to the per-conformer cone_candidates loop this replaces (verified
        # against it in the geometry tests). `previous` is the cone apex, `before` the
        # point defining the axis, exactly as in the scalar call's argument order.
        candidates = cone_candidates_batch(before, previous, angles=grid, per_angle=per_angle)
        # The cone apex: cone_candidates puts every row of its output exactly one bond length
        # from `previous`, which is what licenses the cull in _nearest_obstacle_distance.
        return candidates, _angle_weights(grid, per_angle), np.asarray(previous, dtype=np.float64)

    def _close(
        self,
        plan: _WalkPlan,
        coords: np.ndarray,
        alive: np.ndarray,
        rungs: np.ndarray,
        obstacle_tree: cKDTree | None,
        rng: np.random.Generator,
    ) -> tuple[np.ndarray, str | None]:
        """Place the final residue exactly one bond from the far anchor.

        Parameters
        ----------
        plan
            The resolved plan; must have an end anchor.
        coords
            ``(size, n_residues, 3)`` working coordinates, modified in place.
        alive
            ``(size,)`` mask of conformers still growing.
        rungs
            ``(size,)`` loosest clash threshold used so far, updated in place.
        obstacle_tree
            Spatial index over obstacles, or ``None``.
        rng
            Seeded generator.

        Returns
        -------
        tuple
            ``(ok, reason)``: mask over the whole batch, and a failure description or
            ``None``.

        Notes
        -----
        The last residue X must satisfy ``|X - P| = |X - E| = b`` for the last grown
        residue P and the far anchor E. Those two spheres of equal radius intersect in a
        circle lying in the plane that perpendicularly bisects ``PE``, centred at the
        midpoint, of radius ``sqrt(b ** 2 - (d / 2) ** 2)`` for ``d = |PE|``. Solving it
        analytically rather than sampling for it is what makes both junction bonds exact:
        every point on that circle is at *precisely* one bond length from P and from E,
        so the coincident-anchor-atom defect (c) cannot recur and no bond-length tolerance
        is consumed.

        The pseudo-angle at X is ``2 * asin(d / (2 * b))``, fixed by ``d`` alone, which is
        why the previous step's corridor is exactly ``[min_reach(2), max_reach(2)]``. The
        pseudo-angle at P depends on which point of the circle is used, so those are
        filtered and weighted here.
        """
        end = _require_anchor(plan.end_anchor)
        n = plan.n_residues
        index = n - 1
        live = np.flatnonzero(alive)
        count = _candidates_per_step()

        # The point the closure residue attaches to, and the one before it (needed for the
        # pseudo-angle at the attachment point). For a one-residue region the attachment
        # point *is* the N-anchor, so the angle being filtered here is that anchor's junction
        # angle and the residue before it is the anchor's own neighbour -- which is exactly
        # what start_outer is. When the caller did not supply it there is nothing to form the
        # angle with, and the junction is unconstrained (and warned about).
        if n == 1:
            attach = np.broadcast_to(_require_anchor(plan.start_anchor), (live.size, 3))
            before: np.ndarray | None = (
                np.broadcast_to(plan.start_outer, (live.size, 3))
                if plan.start_outer is not None
                else None
            )
        else:
            attach = coords[live, n - 2]
            if n == 2:
                before = np.broadcast_to(_require_anchor(plan.start_anchor), (live.size, 3))
            else:
                before = coords[live, n - 3]

        axis = end[None, :] - attach
        separation = np.linalg.norm(axis, axis=1)
        radius_squared = CA_CA_BOND_LENGTH**2 - (separation / 2.0) ** 2
        solvable = radius_squared > 0.0

        candidates = np.full((live.size, count, 3), np.nan, dtype=np.float64)
        # Stratified azimuths with a random offset per conformer: an even sweep of the
        # circle covers the accessible pseudo-angles at the attachment point far better
        # than the same number of independent uniform draws.
        base = 2.0 * np.pi * (np.arange(count, dtype=np.float64) / count)
        offsets = rng.random(live.size) * 2.0 * np.pi
        for row in range(live.size):
            if not solvable[row]:
                continue
            radius = math.sqrt(float(radius_squared[row]))
            centre = attach[row] + 0.5 * axis[row]
            frame = align_frame(attach[row], end)
            phi = base + offsets[row]
            local = np.stack([radius * np.cos(phi), radius * np.sin(phi), np.zeros(count)], axis=1)
            candidates[row] = apply_rotation(local, frame) + centre

        # Filter and weight by the pseudo-angle at the attachment point.
        finite = np.all(np.isfinite(candidates), axis=2)
        if before is None:
            weights: np.ndarray | None = None
            extra: np.ndarray = np.broadcast_to(solvable[:, None], (live.size, count)).copy()
        else:
            angles = _angle_at(attach, before, candidates)
            extra = (
                (angles >= BACKBONE_ANGLE_MIN)
                & (angles <= BACKBONE_ANGLE_MAX)
                & solvable[:, None]
                & finite
            )
            weights = _gaussian_weight(angles)

        if plan.end_outer is not None:
            # The other junction: the pseudo-angle centred on the C-anchor is
            # (this candidate, C-anchor, the residue beyond it). It is fixed by which point
            # of the closure circle is used, so it is filtered here or not at all -- and
            # unfiltered it ran from 38.6 to 176.8 degrees in the assembled chain.
            junction = _angle_at(
                np.broadcast_to(end, (live.size, 3)),
                np.broadcast_to(plan.end_outer, (live.size, 3)),
                candidates,
            )
            extra = (
                extra & (junction >= BACKBONE_ANGLE_MIN) & (junction <= BACKBONE_ANGLE_MAX) & finite
            )
            junction_weight = _gaussian_weight(junction)
            weights = junction_weight if weights is None else weights * junction_weight

        chosen, ok, rung = self._select(
            candidates,
            weights=weights,
            extra_mask=extra,
            reference=None,
            goal=None,
            chain_points=self._chain_points(plan, coords, live, index),
            obstacle_tree=obstacle_tree,
            # The closure circle is the intersection of the two bond-length spheres around
            # `attach` and `end`, so every point on it is exactly one bond length from `attach`
            # -- the same guarantee the cone gives, so the same cull applies.
            candidate_centres=attach,
            rng=rng,
        )
        coords[live[ok], index] = chosen[ok]
        rungs[live[ok]] = np.minimum(rungs[live[ok]], rung[ok])

        result = np.zeros(alive.shape[0], dtype=bool)
        result[live[ok]] = True
        reason = None
        if not np.all(ok):
            reason = (
                "the closure step found no position one bond length from both the last "
                "grown residue and the C-terminal anchor with a valid pseudo-angle"
            )
        return result, reason

    # ------------------------------------------------------------------
    # Candidate filtering and selection
    # ------------------------------------------------------------------

    def _select(
        self,
        candidates: np.ndarray,
        *,
        weights: np.ndarray | None,
        extra_mask: np.ndarray | None,
        reference: np.ndarray | None,
        goal: _Goal | None,
        chain_points: np.ndarray,
        obstacle_tree: cKDTree | None,
        candidate_centres: np.ndarray | None,
        rng: np.random.Generator,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Choose one candidate per conformer, or report that none is acceptable.

        Parameters
        ----------
        candidates
            ``(n_live, n_candidates, 3)`` positions to choose among.
        weights
            Relative preference per candidate, ``(n_candidates,)`` or
            ``(n_live, n_candidates)``, or ``None`` for uniform.
        extra_mask
            Additional hard mask, ``(n_live, n_candidates)``, or ``None``.
        reference
            ``(n_live, 3)`` point that ``goal`` measures distance from.
        goal
            Distance corridor and preference, or ``None`` for an unconstrained step.
        chain_points
            ``(n_live, n_points, 3)`` already-placed positions this residue must not
            clash with. May be empty along axis 1.
        obstacle_tree
            Spatial index over external obstacles, or ``None``.
        candidate_centres
            ``(n_live, 3)`` point each row's candidates all sit one bond length from, or
            ``None``. Passed straight to :meth:`_nearest_obstacle_distance`, which uses it to
            skip rows that provably cannot clash.
        rng
            Seeded generator.

        Returns
        -------
        tuple
            ``(chosen, ok, rung)``: ``(n_live, 3)`` positions (NaN where none was found),
            ``(n_live,)`` success mask, and ``(n_live,)`` clash threshold each conformer
            had to accept against the *external obstacles*.

        Notes
        -----
        Two clash tests, not one, and the split is the point. Distances to the region's own
        placed residues, its anchors and its flanking residues are held to
        :data:`~dodo.constants.CA_CLASH_DISTANCE` unconditionally: two alpha carbons of the
        same chain at 2.5 A is not a concession, it is an overlap, and it is exactly what the
        old single test produced -- the relaxation ladder was walked for every request, and
        the final acceptance check then compared against whichever rung had been used, so
        relaxed geometry could not be rejected by construction.

        Only distances to the external obstacle set walk the ladder, and only when there is
        an obstacle set at all. Those coordinates are all-atom positions of already-placed
        residues, where a 2.8 A CA-to-heavy-atom approach is tight rather than impossible,
        and every relaxation is reported through
        :attr:`~dodo.engines.base.IDRResult.relaxed_to`.
        """
        n_live, count, _ = candidates.shape
        acceptable = np.ones((n_live, count), dtype=bool)
        if extra_mask is not None:
            acceptable &= extra_mask

        # Preferences are accumulated in logs and only exponentiated once, after the row is
        # shifted by its own maximum. In linear space the distance preference underflows to
        # exactly 0.0 for every candidate whenever the target is many steering widths away,
        # which turns the weighted draw into a uniform one without saying so.
        log_weight = np.zeros((n_live, count), dtype=np.float64)
        if weights is not None:
            with np.errstate(divide="ignore"):
                log_weight = log_weight + np.log(
                    np.maximum(np.asarray(weights, dtype=np.float64), 0.0)
                )

        if goal is not None:
            if reference is None:
                raise EngineError("A distance goal was given with no reference point.")
            distance = np.linalg.norm(candidates - reference[:, None, :], axis=2)
            # Hard reachability funnel: a candidate from which the region's endpoint can no
            # longer be reached is not a slightly worse choice, it is a dead end several
            # steps later. Dropping it here is what turns a 40-attempt flail into a build.
            acceptable &= (distance >= goal.lo) & (distance <= goal.hi)
            log_weight = log_weight + goal.log_weight(distance)

        # The region's own trace: hard threshold, no ladder, ever.
        acceptable &= self._chain_clear_mask(candidates, chain_points, candidate_centres)

        rung_used = np.full(n_live, CA_CLASH_DISTANCE, dtype=np.float64)
        if obstacle_tree is None:
            chosen_mask = acceptable
        else:
            nearest = self._nearest_obstacle_distance(candidates, obstacle_tree, candidate_centres)
            chosen_mask = acceptable & (nearest >= CA_CLASH_DISTANCE)
        resolved = np.any(chosen_mask, axis=1)

        # Shift each row by the best log-weight among the candidates that are actually
        # available to it, so the largest weight in every row is exactly 1 and the ratios
        # between the surviving candidates are preserved however far the row sits from its
        # preferred distance. Rejected candidates get weight 0 and cannot be picked.
        masked = np.where(chosen_mask, log_weight, -np.inf)
        best = np.max(masked, axis=1, keepdims=True)
        usable = chosen_mask & np.isfinite(masked)
        with np.errstate(invalid="ignore"):
            effective = np.where(usable, np.exp(np.minimum(masked - best, 0.0)), 0.0)
        # A row whose every surviving candidate has weight exactly zero (an angle weight of
        # zero, say) still has to make a choice among them rather than silently returning
        # candidate 0: without this the racing key below is +inf everywhere and argmin picks
        # index 0, which is a clashing candidate written out as a successful placement.
        starved = np.any(chosen_mask, axis=1) & ~np.any(effective > 0.0, axis=1)
        effective[starved] = np.where(chosen_mask[starved], 1.0, 0.0)

        # Weighted choice by exponential racing: for independent U ~ Uniform(0, 1) the
        # index minimizing -log(U) / w is drawn with probability proportional to w. One
        # vectorized draw for the whole batch, and no per-conformer normalization.
        draws = rng.random((n_live, count))
        with np.errstate(divide="ignore", invalid="ignore"):
            # np.where evaluates both branches, so the zero weights of rejected candidates
            # divide before they are discarded. Their result is thrown away either way.
            key = np.where(effective > 0.0, -np.log(draws) / effective, np.inf)
        picked = np.argmin(key, axis=1)

        chosen = np.full((n_live, 3), np.nan, dtype=np.float64)
        ok = resolved & np.isfinite(key[np.arange(n_live), picked])
        chosen[ok] = candidates[np.arange(n_live), picked][ok]
        rung_used[~ok] = CA_CLASH_DISTANCE
        return chosen, ok, rung_used

    @staticmethod
    def _nearest_obstacle_distance(
        candidates: np.ndarray,
        obstacle_tree: cKDTree,
        centres: np.ndarray | None = None,
    ) -> np.ndarray:
        """Distance from every candidate to the nearest *external* obstacle atom.

        Parameters
        ----------
        candidates
            ``(n_live, n_candidates, 3)`` positions.
        obstacle_tree
            Spatial index over external obstacles, shared by the whole batch.
        centres
            ``(n_live, 3)`` point each row's candidates all sit exactly
            :data:`~dodo.constants.CA_CA_BOND_LENGTH` from, or ``None`` to skip the cull. Every
            candidate generator in this engine has such a point; see
            :data:`_CANDIDATE_CLEAR_DISTANCE`.

        Returns
        -------
        np.ndarray
            ``(n_live, n_candidates)`` distances, ``inf`` for a non-finite candidate, which
            is rejected on other grounds anyway, and ``inf`` for every candidate of a row whose
            centre is far enough from the obstacles that none of them can clash.

        Notes
        -----
        The ``inf`` returned for a culled row is not the true distance, and deliberately so: the
        caller compares it against :data:`~dodo.constants.CA_CLASH_DISTANCE` and nothing else, so
        substituting ``inf`` for "provably at least the clash distance" cannot change a single
        decision. That is what makes the cull free rather than approximate -- output stays
        bit-identical for a given seed.

        This is the hot loop of the whole package. Culling here replaces up to 355 queries with
        one, and on p300 it removes 94.4% of the query points; see
        :data:`_CANDIDATE_CLEAR_DISTANCE` for the measurement.
        """
        n_live, count, _ = candidates.shape
        rows = np.arange(n_live)
        if centres is not None:
            centre_array = np.asarray(centres, dtype=np.float64)
            usable = np.all(np.isfinite(centre_array), axis=1)
            clearance = np.zeros(n_live, dtype=np.float64)
            if np.any(usable):
                found, _ = obstacle_tree.query(centre_array[usable])
                clearance[usable] = found
            # A non-finite centre is left at clearance 0, so its row is queried the slow way
            # rather than trusted. Its candidates are non-finite too and get inf regardless.
            rows = np.flatnonzero(clearance < _CANDIDATE_CLEAR_DISTANCE)

        nearest = np.full((n_live, count), np.inf, dtype=np.float64)
        if rows.size == 0:
            return nearest

        flat = candidates[rows].reshape(-1, 3)
        finite = np.all(np.isfinite(flat), axis=1)
        distances = np.full(flat.shape[0], np.inf)
        if np.any(finite):
            # workers=-1 only above _PARALLEL_QUERY_MINIMUM.
            workers = -1 if int(np.count_nonzero(finite)) >= _PARALLEL_QUERY_MINIMUM else 1
            found, _ = obstacle_tree.query(flat[finite], workers=workers)
            distances[finite] = found
        nearest[rows] = distances.reshape(rows.size, count)
        return nearest

    @staticmethod
    def _chain_clear_mask(
        candidates: np.ndarray,
        chain_points: np.ndarray,
        centres: np.ndarray | None,
    ) -> np.ndarray:
        """Mask the candidates clear of every already-placed CA of their own chain.

        Parameters
        ----------
        candidates
            ``(n_live, n_candidates, 3)`` positions, each one bond length from its row's
            entry in ``centres``.
        chain_points
            ``(n_live, n_points, 3)`` positions this residue must not clash with -- the
            conformer's own placed residues, its anchors and its flanking residues, with the
            sequence-local neighbours already excluded.
        centres
            ``(n_live, 3)`` apex each row's candidates orbit, one bond length away.

        Returns
        -------
        np.ndarray
            ``(n_live, n_candidates)`` boolean, ``True`` where the candidate is at least
            :data:`~dodo.constants.CA_CLASH_DISTANCE` from every chain point.

        Notes
        -----
        Replaces a per-conformer :class:`~scipy.spatial.cKDTree` built and queried once per
        residue -- the walk's single hottest cost -- with one vectorized pass over the whole
        live set. It returns the mask ``_nearest_chain_distance(...) >= CA_CLASH_DISTANCE``
        used to produce, which is all the selection consumed; the raw distances went nowhere
        else.

        The pass is exact where it counts. Every candidate sits one bond length
        (:data:`~dodo.constants.CA_CA_BOND_LENGTH`) from its apex, so a chain point farther
        than ``CA_CA_BOND_LENGTH + CA_CLASH_DISTANCE`` from that apex is farther than the
        clash distance from *every* candidate in the row and cannot decide the mask.

        """
        n_live, count, _ = candidates.shape
        clear = np.ones((n_live, count), dtype=bool)
        n_points = chain_points.shape[1]
        if n_points == 0:
            return clear
        if centres is None:
            raise EngineError(
                "_chain_clear_mask needs the candidate centres to cull far chain points; "
                "both the cone step and the closure step supply them."
            )
        finite = np.all(np.isfinite(candidates), axis=2)

        # Cull chain points that cannot reach any candidate in their row, then gather the
        # survivors into a dense (n_live, k) block -- k is the most any row keeps.
        reach = CA_CA_BOND_LENGTH + CA_CLASH_DISTANCE
        to_apex_sq = np.sum((chain_points - centres[:, None, :]) ** 2, axis=2)
        near = to_apex_sq <= reach * reach
        k = int(near.sum(axis=1).max())
        if k == 0:
            return clear
        order = np.argsort(~near, axis=1, kind="stable")[:, :k]
        picked = np.take_along_axis(chain_points, order[:, :, None], axis=1)
        picked_near = np.take_along_axis(near, order, axis=1)

        diff = candidates[:, :, None, :] - picked[:, None, :, :]
        squared = np.sum(diff * diff, axis=3)
        squared = np.where(picked_near[:, None, :], squared, np.inf)
        nearest = np.sqrt(np.min(squared, axis=2))
        # A non-finite candidate hit nothing under the old tree query (its distance came back
        # inf, so it passed the threshold); keep the degenerate rows unchanged too.
        cleared: np.ndarray = (nearest >= CA_CLASH_DISTANCE) | ~finite
        return cleared

    # ------------------------------------------------------------------
    # Small helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _prune_unreachable(request: IDRRequest, obstacles: np.ndarray | None) -> np.ndarray | None:
        """Drop obstacles this region provably cannot reach, before the tree is built.

        Parameters
        ----------
        request
            The region about to be built; supplies the anchors and the residue count.
        obstacles
            ``(n_obstacles, 3)`` coordinates, or ``None``.

        Returns
        -------
        np.ndarray or None
            The surviving obstacles, or ``obstacles`` unchanged when nothing can be dropped.

        Notes
        -----
        Conservative by construction: :func:`reachable_envelope` bounds where the region's alpha
        carbons can go, and the envelope is then dilated by
        :data:`~dodo.constants.CA_CLASH_DISTANCE` -- the only threshold ever applied to the
        obstacle set -- before anything is dropped. An obstacle outside the dilated envelope is
        further than the clash distance from every position the region can occupy, so it can
        neither reject a candidate in :meth:`_select` nor fail a conformer in
        :func:`_validate_conformer`. Output is therefore bit-identical, which is how this is
        tested rather than merely asserted.

        """
        if obstacles is None:
            return None
        array = np.asarray(obstacles, dtype=np.float64)
        if array.ndim != 2 or array.shape[1] != 3 or array.shape[0] == 0:
            # Malformed input is _obstacle_tree's business to diagnose, not this method's.
            return obstacles
        if not np.all(np.isfinite(array)):
            # A NaN coordinate makes every comparison below False, which would drop the row and
            # silently disarm the non-finite guard in _obstacle_tree -- the one thing standing
            # between a NaN obstacle and a region grown straight through a domain. Prune nothing
            # and let that guard fire.
            return obstacles
        limit = reachable_envelope(request.n_residues, request.n_anchor_xyz, request.c_anchor_xyz)
        if limit is None:
            return obstacles
        n_anchor, c_anchor = request.n_anchor_xyz, request.c_anchor_xyz
        keep: np.ndarray
        if n_anchor is not None and c_anchor is not None:
            reach = np.linalg.norm(array - n_anchor, axis=1) + np.linalg.norm(
                array - c_anchor, axis=1
            )
            keep = reach <= limit + 2.0 * CA_CLASH_DISTANCE
        else:
            # reachable_envelope only returns a finite limit for a region with at least one
            # anchor, so exactly one of the two is set here.
            anchor = n_anchor if n_anchor is not None else c_anchor
            if anchor is None:
                return obstacles
            keep = np.linalg.norm(array - anchor, axis=1) <= limit + CA_CLASH_DISTANCE
        if bool(np.all(keep)):
            # The common case for a long region. Return the original so the caller's array is
            # not needlessly copied.
            return obstacles
        pruned: np.ndarray = array[keep]
        return pruned

    @staticmethod
    def _obstacle_tree(obstacles: np.ndarray | None) -> cKDTree | None:
        """Build a spatial index over the obstacle coordinates, or return ``None``.

        Parameters
        ----------
        obstacles
            ``(n_obstacles, 3)`` coordinates, or ``None``.

        Returns
        -------
        scipy.spatial.cKDTree or None
            The index, or ``None`` when there is nothing to avoid.

        Raises
        ------
        GeometryError
            If the obstacle array is not ``(n, 3)`` or contains non-finite values. A NaN
            obstacle silently never matches, so a whole domain could vanish from clash
            checking and the walk would grow straight through it.
        """
        if obstacles is None:
            return None
        array = np.asarray(obstacles, dtype=np.float64)
        if array.ndim != 2 or array.shape[1] != 3:
            raise GeometryError(f"obstacles must have shape (n, 3), got {array.shape}.")
        if array.shape[0] == 0:
            return None
        if not np.all(np.isfinite(array)):
            bad = int(np.count_nonzero(~np.all(np.isfinite(array), axis=1)))
            raise GeometryError(
                f"obstacles contains non-finite values in {bad} of {array.shape[0]} rows. "
                f"Such a point cannot clash with anything, so the region it belongs to "
                f"would be invisible to clash checking."
            )
        return cKDTree(array)

    @staticmethod
    def _check_start_is_free(plan: _WalkPlan, obstacle_tree: cKDTree | None) -> None:
        """Raise if the fixed starting point of a free region is already occupied.

        A free region is built at the world origin, so an obstacle set that covers the
        origin makes every attempt fail identically. Say so immediately instead of burning
        the whole attempt budget to discover it.
        """
        if plan.start_anchor is not None or obstacle_tree is None:
            return
        distance, _ = obstacle_tree.query(np.zeros((1, 3)))
        if float(distance[0]) < CA_CLASH_DISTANCE:
            raise GeometryError(
                f"A region with no anchors is built in its own frame starting at the "
                f"origin, but an obstacle sits {float(distance[0]):.2f} A from the origin "
                f"(minimum {CA_CLASH_DISTANCE:.2f} A). Either pass obstacles=None and "
                f"place the region afterwards, or give it an anchor to grow from."
            )

    @staticmethod
    def _reference_for(plan: _WalkPlan, coords: np.ndarray, live: np.ndarray) -> np.ndarray | None:
        """Return the per-conformer point that this step's distance goal measures from.

        The far anchor in closure mode; residue 0 otherwise, since the end-to-end distance
        of a free or terminal region is measured from there.
        """
        if plan.end_anchor is not None:
            return np.broadcast_to(plan.end_anchor, (live.size, 3))
        if plan.n_residues < 2:
            return None
        origin: np.ndarray = coords[live, 0]
        return origin

    @staticmethod
    def _distance_from_origin(
        plan: _WalkPlan, coords: np.ndarray, live: np.ndarray, index: int
    ) -> np.ndarray | None:
        """Distance of the last placed residue from residue 0, for the catch-up schedule.

        ``None`` in closure mode, where the schedule is seeded from the anchors and not from
        an end-to-end target, and for the first two residues, where there is nothing placed
        to measure yet.
        """
        if plan.closes or index < 1:
            return None
        distances: np.ndarray = np.linalg.norm(coords[live, index - 1] - coords[live, 0], axis=-1)
        return distances

    @staticmethod
    def _chain_points(
        plan: _WalkPlan, coords: np.ndarray, live: np.ndarray, index: int
    ) -> np.ndarray:
        """Collect the placed positions that residue ``index`` must not clash with.

        Parameters
        ----------
        plan
            The resolved plan.
        coords
            ``(size, n_residues, 3)`` working coordinates.
        live
            Indices of the conformers still growing.
        index
            Residue index about to be placed, in growth order.

        Returns
        -------
        np.ndarray
            ``(n_live, n_points, 3)`` positions.

        Notes
        -----
        Residues within :data:`dodo.constants.CLASH_EXCLUDE_WITHIN_RESIDUES` of ``index``
        are excluded: ``|i - j| == 1`` are covalently bonded at 3.8 A and ``|i - j| == 2``
        are angle-constrained to 5.0-7.5 A, so counting either as a clash would make every
        valid backbone look like a collision.

        The anchors take part in this exclusion as virtual residues at ``-1`` and
        ``n_residues``, and the flanking residues beyond them as ``-2`` and
        ``n_residues + 1``. That is what stops a mid-region residue from being placed on top
        of an anchor while still allowing the junction bonds themselves -- and it is why the
        flanking residues do not need to be passed in ``obstacles`` as well: they are part of
        the chain, so they are held to the strict threshold rather than the relaxable one.
        """
        exclude = CLASH_EXCLUDE_WITHIN_RESIDUES
        kept = coords[live, : max(0, index - exclude)]
        extras = []
        if plan.start_anchor is not None and index - (-1) > exclude:
            extras.append(np.broadcast_to(plan.start_anchor, (live.size, 1, 3)))
        if plan.start_outer is not None and index - (-2) > exclude:
            extras.append(np.broadcast_to(plan.start_outer, (live.size, 1, 3)))
        if plan.end_anchor is not None and plan.n_residues - index > exclude:
            extras.append(np.broadcast_to(plan.end_anchor, (live.size, 1, 3)))
        if plan.end_outer is not None and plan.n_residues + 1 - index > exclude:
            extras.append(np.broadcast_to(plan.end_outer, (live.size, 1, 3)))
        if not extras:
            return kept
        stacked: np.ndarray = np.concatenate([*extras, kept], axis=1)
        return stacked


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------


def _per_angle() -> int:
    """Candidate positions per cone ring, capped by the per-residue candidate budget."""
    grid_size = int(backbone_angle_grid().size)
    return max(1, min(CANDIDATES_PER_ANGLE, MAX_CANDIDATES_PER_RESIDUE // grid_size))


def _candidates_per_step() -> int:
    """Total candidates evaluated per residue."""
    return int(backbone_angle_grid().size) * _per_angle()


def _natural_fluctuation(n_bonds: int) -> float:
    """Typical fluctuation in the span of an ``n_bonds`` subchain, in Angstroms."""
    return STEERING_SIGMA_FACTOR * CA_CA_BOND_LENGTH * math.sqrt(max(n_bonds, 1))


def _target_steering_width(target: np.ndarray, span: int, mean_target: float) -> np.ndarray:
    """Steering width on the distance to residue 0, for a free or terminal region.

    Parameters
    ----------
    target
        ``(rows, 1)`` end-to-end target of each live conformer, in Angstroms.
    span
        Bonds between the region's first and last residue.
    mean_target
        The end-to-end distance the *region* was asked for, in Angstroms -- the mean the batch
        is spread around, not any one conformer's draw.

    Returns
    -------
    np.ndarray
        ``(rows, 1)`` width in Angstroms, applied symmetrically about the schedule.

    Notes
    -----
    Three terms, each doing one job. :data:`STEERING_SLACK_FRACTION` of the radial room the
    target leaves unused is the width the walk would like -- wide enough that a whole ring of
    bond directions reaches the distance it wants, instead of the two that a sharp radial
    preference leaves. :data:`STEERING_WIDTH_TARGET_CAP` keeps that from costing more accuracy
    than it is worth on a short region. :data:`TARGET_STEERING_WIDTH` floors it, so a target
    close to the geometric ceiling is still steered as tightly as it always was.

    The slack term is per conformer and the cap is per region, and that split is load-bearing.
    How much radial room a conformer has is a property of *its own* draw, so the slack reads
    ``target``. Whether the region is big enough for a flat zig-zag to be worth paying accuracy
    for is a property of the *region*, so the cap reads ``mean_target`` -- the value the caller
    asked for, before :func:`sample_end_to_end_targets` spread the batch around it.

    Keying the cap on the per-conformer draw instead is a bug that only shows up at
    ``n_models > 1``, and it was one. The draws have a coefficient of variation near 0.42, so in
    a batch of twenty a conformer can draw a target several times more compact than the mean;
    capping at a fraction of *that* handed it a width near the floor and it came out planar
    while its nineteen siblings did not. Measured on p300's 583-residue terminal region at
    twenty models: the two most compact draws (27.9 and 72.1 A against a mean of 186 A) received
    widths of 0.56 and 1.44 A and planar orders of 0.474 and 0.132, against 0.011-0.082 for the
    other eighteen -- a correlation between width and planar order of -0.96. A region does not
    become short because one of its conformers drew a compact target.
    """
    values = np.asarray(target, dtype=np.float64)
    reach = max_reach(span) if span >= 1 else 0.0
    # A zero reach means there is no span to steer along; the floor is the whole answer.
    slack = 1.0 - np.clip(values / reach, 0.0, 1.0) if reach > 0.0 else np.zeros_like(values)
    wanted = STEERING_SLACK_FRACTION * _MAX_ADVANCE_PER_BOND * slack
    ceiling = STEERING_WIDTH_TARGET_CAP * float(mean_target)
    width: np.ndarray = np.maximum(np.minimum(wanted, ceiling), TARGET_STEERING_WIDTH)
    return width


def _headroom_width(gap: float) -> float:
    """Steering width for the direction in which overshooting cannot be recovered from.

    Parameters
    ----------
    gap
        Distance from the schedule to the hard bound in that direction, in Angstroms.

    Returns
    -------
    float
        The width to use.

    Notes
    -----
    Floored at half a bond length. Where the corridor itself is only a couple of Angstroms
    wide -- the last step before a closure -- the corridor is doing the work and a width
    much below it would concentrate every conformer on the same pseudo-angle.
    """
    return max(HEADROOM_STEERING_FRACTION * abs(gap), 0.5 * CA_CA_BOND_LENGTH)


def _angle_weights(grid: np.ndarray, per_angle: int) -> np.ndarray:
    """Relative preference for each cone candidate, from the measured angle distribution.

    Parameters
    ----------
    grid
        Pseudo-angles in degrees, in the order :func:`cone_candidates` was given them.
    per_angle
        Candidates generated per angle.

    Returns
    -------
    np.ndarray
        ``(len(grid) * per_angle,)`` weights.

    Notes
    -----
    Weighting by the measured density -- mean :data:`dodo.constants.BACKBONE_ANGLE_MEAN`,
    standard deviation :data:`~dodo.constants.BACKBONE_ANGLE_SD` -- rather than sampling
    the window uniformly. Uniform sampling would put as many 91-degree turns in the output
    as 125-degree ones, giving a backbone that satisfies the window while having the wrong
    angle distribution and therefore the wrong persistence length.

    :func:`dodo.constants.backbone_angle_grid` returns angles ideal-first and
    :func:`~dodo.geometry.sampling.cone_candidates` preserves that order block-wise, so
    ``np.repeat`` here lines each weight up with its own ring.
    """
    weights: np.ndarray = np.repeat(_gaussian_weight(grid), per_angle)
    return weights


def _gaussian_weight(
    values: np.ndarray,
    centre: float = BACKBONE_ANGLE_MEAN,
    width: float = BACKBONE_ANGLE_SD,
) -> np.ndarray:
    """Unnormalized Gaussian weight, used for both angle and distance preferences."""
    scaled = (np.asarray(values, dtype=np.float64) - centre) / width
    weight: np.ndarray = np.exp(-0.5 * scaled * scaled)
    return weight


def _angle_at(vertex: np.ndarray, first: np.ndarray, second: np.ndarray) -> np.ndarray:
    """Angle in degrees at ``vertex``, between the rays to ``first`` and to ``second``.

    Parameters
    ----------
    vertex
        ``(n, 3)`` vertex positions.
    first
        ``(n, 3)`` first neighbour.
    second
        ``(n, m, 3)`` candidate second neighbours.

    Returns
    -------
    np.ndarray
        ``(n, m)`` angles in degrees.

    Notes
    -----
    The dot product of the normalized rays is clipped to ``[-1, 1]`` before ``arccos``.
    Without the clip, a rounding error of 1e-16 above 1 makes ``arccos`` return NaN, and
    a NaN angle compares False against both window bounds -- so a perfectly straight
    junction would be silently rejected rather than accepted.
    """
    back = first - vertex
    out = second - vertex[:, None, :]
    back_norm = np.linalg.norm(back, axis=1)[:, None]
    out_norm = np.linalg.norm(out, axis=2)
    with np.errstate(invalid="ignore", divide="ignore"):
        cosine = np.einsum("ij,imj->im", back, out) / (back_norm * out_norm)
    degrees: np.ndarray = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    return degrees


def _check_fixed_context(request: IDRRequest) -> None:
    """Raise if the fixed residues around the region are not valid backbone geometry.

    Parameters
    ----------
    request
        The region to build.

    Raises
    ------
    GeometryError
        If two of the fixed CAs that are more than
        :data:`~dodo.constants.CLASH_EXCLUDE_WITHIN_RESIDUES` residues apart sit closer
        than :data:`~dodo.constants.CA_CLASH_DISTANCE`, or if an outer residue is not
        roughly one bond length from its anchor.

    """
    fixed: list[tuple[int, str, np.ndarray]] = []
    if request.n_anchor_prev_xyz is not None:
        fixed.append((-2, "n_anchor_prev_xyz", request.n_anchor_prev_xyz))
    if request.n_anchor_xyz is not None:
        fixed.append((-1, "n_anchor_xyz", request.n_anchor_xyz))
    if request.c_anchor_xyz is not None:
        fixed.append((request.n_residues, "c_anchor_xyz", request.c_anchor_xyz))
    if request.c_anchor_next_xyz is not None:
        fixed.append((request.n_residues + 1, "c_anchor_next_xyz", request.c_anchor_next_xyz))

    for outer, anchor in (
        (request.n_anchor_prev_xyz, request.n_anchor_xyz),
        (request.c_anchor_next_xyz, request.c_anchor_xyz),
    ):
        if outer is None or anchor is None:
            continue
        bond = float(np.linalg.norm(outer - anchor))
        if abs(bond - CA_CA_BOND_LENGTH) > _OUTER_BOND_SLACK:
            raise GeometryError(
                f"The chain is broken next to this region's anchor: the anchor and the fixed "
                f"residue beyond it are {bond:.2f} A apart, where consecutive alpha carbons are "
                f"{CA_CA_BOND_LENGTH:.2f} A. That residue exists only to pin the pseudo-angle "
                f"centred on the anchor, and a separation this far off means it is not the "
                f"anchor's chain neighbour. This is a problem with the coordinates handed in, "
                f"not with the region: no conformer can repair a break in the fixed chain, so "
                f"the region is left as it arrived. Real AlphaFold models do contain these -- "
                f"one measured example has a single 5.26 A step in an otherwise sound chain."
            )

    for i, (index_a, name_a, point_a) in enumerate(fixed):
        for index_b, name_b, point_b in fixed[i + 1 :]:
            if abs(index_a - index_b) <= CLASH_EXCLUDE_WITHIN_RESIDUES:
                continue
            separation = float(np.linalg.norm(point_a - point_b))
            if separation < CA_CLASH_DISTANCE:
                raise GeometryError(
                    f"The fixed residues {name_a} and {name_b} are {separation:.2f} A apart, "
                    f"closer than the {CA_CLASH_DISTANCE:.2f} A two non-bonded alpha carbons "
                    f"can be, and they are {abs(index_a - index_b)} residues apart along the "
                    f"chain. That is a problem with the coordinates handed in, not with the "
                    f"region: no conformer of the region can repair it."
                )


def _require_anchor(anchor: np.ndarray | None) -> np.ndarray:
    """Return ``anchor``, raising if it is missing.

    The plan's case analysis guarantees an anchor exists wherever this is called, so a
    failure here is a bug in this module rather than bad input -- but returning ``None``
    into coordinate arithmetic would produce a TypeError three frames away, and silently
    substituting the origin would produce a plausible wrong structure.
    """
    if anchor is None:
        raise EngineError(
            "Internal error: a growth step needed a fixed anchor that the plan does not "
            "have. This is a DODO bug; please report it."
        )
    return anchor


def _validate_conformer(
    coords: np.ndarray,
    plan: _WalkPlan,
    target: float,
    obstacle_clearance: float,
    obstacle_threshold: float,
) -> str | None:
    """Measure a finished conformer and return the first physical violation, or ``None``.

    Parameters
    ----------
    coords
        ``(n_residues, 3)`` generated CA positions, in growth order.
    plan
        The resolved plan, for the anchors and the flanking residues.
    target
        This conformer's own end-to-end target in Angstroms.
    obstacle_clearance
        Measured distance from this conformer to the nearest external obstacle atom, or
        ``inf`` when there are no obstacles.
    obstacle_threshold
        The clash threshold this conformer was granted against the external obstacles: the
        strict :data:`~dodo.constants.CA_CLASH_DISTANCE` unless a rung of the relaxation
        ladder was needed, in which case that rung -- and the caller is told about it.

    Returns
    -------
    str or None
        A description of the violation, or ``None`` if the conformer is valid.

    Notes
    -----
    Independent verification, deliberately duplicating what the generator guarantees.
    ``success=True`` on an :class:`~dodo.engines.base.IDRResult` from this engine means
    the geometry was *measured* and is valid, not that the code believes it should be. The
    bond, angle and steric measurement is :func:`dodo.geometry.metrics.validate_ca_trace`,
    which is not part of this module -- checking generated geometry with the generator's own
    arithmetic would let a shared mistake pass unnoticed.

    Three things changed here after measurement, all of them cases where the check could not
    fail:

    * The trace validated is the *spliced* one, including the residues beyond the anchors
      where the caller supplied them, so each anchor is an angle **vertex**. Previously each
      anchor sat at one end of the trace and its pseudo-angle was structurally invisible: 41%
      of those angles were below the window floor and none of them was ever measured.
    * The region's own steric check uses :data:`~dodo.constants.CA_CLASH_DISTANCE`, not the
      relaxation rung. Comparing against the rung meant a conformer that relaxed to 2.0 A was
      validated against 2.0 A, so the clash check could not reject what the ladder allowed.
    * The end-to-end audit uses this conformer's own target and a tolerance that is a pure
      fraction of it. It used to reuse a tolerance floored at one bond length -- the same
      floored number the generator steered by -- so it was arithmetically independent and
      logically vacuous.
    """
    if not np.all(np.isfinite(coords)):
        return "the trace contains non-finite coordinates"

    trace, offset = _spliced_trace(coords, plan)
    if trace.shape[0] >= 2:
        report = validate_ca_trace(trace, residue_offset=offset)
        if not report.ok:
            return (
                f"{report.describe()} (the region's residues are 0..{plan.n_residues - 1}; "
                f"negative indices are the fixed residues N-terminal to it)"
            )

    closest = _closest_non_bonded(trace)
    if closest < CA_CLASH_DISTANCE:
        return (
            f"two non-bonded CA atoms of the chain are {closest:.3f} A apart, closer than "
            f"the {CA_CLASH_DISTANCE:.2f} A minimum approach"
        )

    if obstacle_clearance < obstacle_threshold:
        return (
            f"the region comes within {obstacle_clearance:.3f} A of an already-placed atom, "
            f"closer than the {obstacle_threshold:.2f} A threshold in force"
        )

    if not plan.closes and plan.n_residues >= 2:
        achieved = end_to_end(coords)
        tolerance = plan.tolerance_fraction * target
        if abs(achieved - target) > tolerance:
            return (
                f"end-to-end distance {achieved:.2f} A misses this conformer's {target:.2f} A "
                f"target by more than the {tolerance:.2f} A tolerance"
            )
    return None


def _spliced_trace(coords: np.ndarray, plan: _WalkPlan) -> tuple[np.ndarray, int]:
    """Return the generated CAs inside their fixed context, and the index of the first row.

    Parameters
    ----------
    coords
        ``(n_residues, 3)`` generated CA positions, in growth order.
    plan
        The resolved plan.

    Returns
    -------
    tuple
        ``(trace, residue_offset)``. The trace is
        ``[start_outer?, start_anchor?, region..., end_anchor?, end_outer?]`` and the offset
        is what :func:`dodo.geometry.metrics.validate_ca_trace` needs so that the region's
        own residues are numbered from 0, the anchors at -1 and ``n_residues``, and the
        flanking residues at -2 and ``n_residues + 1``.

    Notes
    -----
    The flanking residues are placed at *exactly* one bond length from their anchor along the
    true direction to them, rather than at their own coordinates. The pseudo-angle centred on
    an anchor -- the whole reason they are here -- depends only on that direction, so nothing
    being measured changes; what it buys is that a CA-CA distance in the *input* structure
    which is not exactly 3.80 A (cis-proline sits near 2.9 A) cannot be reported as a
    bond-length violation of geometry this engine did not generate. The real coordinates are
    used for the clash check in :func:`_check_fixed_context`, which runs once at plan time.
    """
    pieces = []
    offset = 0
    if plan.start_anchor is not None:
        if plan.start_outer is not None:
            pieces.append(_at_bond_length(plan.start_anchor, plan.start_outer)[None, :])
            offset -= 1
        pieces.append(plan.start_anchor[None, :])
        offset -= 1
    pieces.append(coords)
    if plan.end_anchor is not None:
        pieces.append(plan.end_anchor[None, :])
        if plan.end_outer is not None:
            pieces.append(_at_bond_length(plan.end_anchor, plan.end_outer)[None, :])
    trace: np.ndarray = np.concatenate(pieces, axis=0)
    return trace, offset


def _at_bond_length(anchor: np.ndarray, outer: np.ndarray) -> np.ndarray:
    """Place ``outer`` at exactly one CA-CA bond length from ``anchor``, same direction."""
    offset = outer - anchor
    length = float(np.linalg.norm(offset))
    if length <= 0.0:
        # Guarded rather than assumed: _check_fixed_context rejects this at plan time, and a
        # zero-length direction here would produce NaN and turn a measured angle into an
        # unmeasured one.
        raise EngineError(
            "Internal error: a flanking residue coincides with its anchor. This is a DODO "
            "bug; please report it."
        )
    scaled: np.ndarray = anchor + offset * (CA_CA_BOND_LENGTH / length)
    return scaled


def _closest_non_bonded(trace: np.ndarray) -> float:
    """Smallest distance between two trace points more than the exclusion apart.

    Returns ``inf`` when the trace is too short to contain a non-bonded pair. The
    exclusion is :data:`dodo.constants.CLASH_EXCLUDE_WITHIN_RESIDUES`: bonded neighbours
    sit at 3.8 A and next-nearest at 5.0-7.5 A by construction, so counting either as a
    clash would report every peptide bond as a violation.
    """
    if trace.shape[0] <= CLASH_EXCLUDE_WITHIN_RESIDUES + 1:
        return float("inf")
    index = np.arange(trace.shape[0])
    separations = np.abs(np.subtract.outer(index, index))
    distances = np.linalg.norm(trace[:, None, :] - trace[None, :, :], axis=2)
    return float(np.min(distances[separations > CLASH_EXCLUDE_WITHIN_RESIDUES]))


if TYPE_CHECKING:
    # Checked by mypy, never executed: if SelfAvoidingWalk ever stops satisfying the
    # engine protocol, that is a type error here rather than a runtime surprise in the
    # orchestrator that selects engines by protocol.
    from .base import ConformationEngine

    _protocol_check: ConformationEngine = SelfAvoidingWalk()
