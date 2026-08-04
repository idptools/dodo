"""The STARLING engine: an idptools STARLING ensemble adapted to DODO's anchored chains.

STARLING (idptools, PyPI ``starling-ensemble``) is a generative model of IDR ensembles.
It is the scientifically strongest generator DODO can call: its conformers come from a
model trained on simulated ensembles, not from a random walk with a distance target
bolted on. That is why it is the *primary* engine when installed.

Three things make wrapping it more than a one-liner.

**It is heavy and optional.** STARLING pulls torch and roughly 2.4 GB of weights. So
nothing here imports it at module scope: the single import site is
:func:`_import_starling`, which is also the seam the tests monkeypatch with a fake.
:meth:`StarlingEngine.available` answers without importing anything, by asking
:mod:`importlib.util` whether the module exists -- an import would cost seconds on a
path (``dodo --help``, engine selection) that must stay instant.

**It generates an ensemble; DODO needs one chain between two fixed anchors.** A STARLING
conformer floats in its own frame with no relation to the folded domains it has to
connect. So this module selects the conformer whose own end-to-end distance is closest
to what the anchor separation demands, then places it as a rigid body such that its
first CA sits exactly one CA-CA bond from the N-anchor and its last CA one bond from the
C-anchor (see :func:`place_between_anchors` for the construction and why it is exact).
When no conformer can satisfy both anchors, that is reported through the success mask and
the measured residual -- never forced, and never returned as if it had worked.

**Its geometry is a model output, not a DODO construction.** STARLING reconstructs 3D
coordinates, so its virtual CA-CA bonds are close to but not exactly 3.80 A and its
pseudo-angles follow the model rather than DODO's sampling window. Conformers are
therefore *screened* (:func:`screen_conformers`) against the physically observed ranges
in :mod:`dodo.constants` before use, and the measured distributions are reported. A
conformer that fails screening is dropped with a note; if none survive, this module
raises rather than returning coordinates it has already measured as unphysical.

Screening looks at three things and every one of them is a hard gate, because each names
geometry with no all-atom realization at all: a bond longer than a trans peptide can make, a
pseudo-angle outside the observed range at *any* vertex, and a non-bonded CA-CA pair inside
:data:`dodo.constants.CA_CLASH_DISTANCE`. None of them is a budget to be spent. The chosen
conformer's *dimension* is a fourth gate, applied after placement: nothing in the rigid-body
construction looks at :attr:`~dodo.engines.base.IDRRequest.target_end_to_end` on a free or
terminal region, so without it the engine answers any request with whatever the ensemble
happens to contain.

Licence note
------------
STARLING's PyPI metadata says MIT while the LICENSE file in its repository says LGPL-3.
Because that ambiguity is unresolved, no STARLING code is vendored, copied or adapted
here. This module only calls the installed package's public entry point, which no licence
restricts.
"""

from __future__ import annotations

import importlib.util
import inspect
import os
import sys
import warnings
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from types import ModuleType
from typing import Any, Final

import numpy as np
from scipy.spatial import cKDTree

from ..constants import (
    BACKBONE_ANGLE_MAX,
    BACKBONE_ANGLE_MIN,
    BACKBONE_ANGLE_OBSERVED_MAX,
    BACKBONE_ANGLE_OBSERVED_MIN,
    CA_CA_BOND_LENGTH,
    CA_CA_BOND_TOLERANCE,
    CA_CLASH_DISTANCE,
    CLASH_EXCLUDE_WITHIN_RESIDUES,
    CLASH_RELAXATION_LADDER,
    STARLING_MAX_LENGTH,
)
from ..exceptions import (
    EngineError,
    EngineUnavailableError,
    ExhaustedAttemptsError,
    GeometryError,
    MissingDependencyError,
)
from ..geometry.metrics import ca_bond_lengths, ca_pseudo_angles
from ..geometry.regularize import regularize_ca_trace
from ..geometry.sampling import random_unit_vectors
from ..geometry.transforms import apply, rotation_between_vectors, rotation_from_axis_angle
from .base import IDRRequest, IDRResult

__all__ = [
    "AnchorPlacement",
    "LengthCap",
    "StarlingEngine",
    "StarlingReport",
    "desired_internal_span",
    "end_to_end_tolerance",
    "place_between_anchors",
    "rank_conformers",
    "screen_conformers",
    "starling_installed",
    "starling_max_length",
]

#: Import name of the STARLING package. Distribution name differs from import name, which
#: is why both are spelled out: ``pip install starling-ensemble`` gives ``import starling``.
STARLING_MODULE = "starling"
STARLING_DISTRIBUTION = "starling-ensemble"

#: Conformers requested from STARLING per :meth:`StarlingEngine.generate` call.
#:
#: CHOICE. Selection needs *spread* in end-to-end distance, not just a good mean: an
#: interior IDR must span whatever distance its flanking domains impose, and the chance
#: that some conformer lands near that value grows with the ensemble size. 100 is
#: STARLING's own documented default order of magnitude and costs one batched forward
#: pass.
DEFAULT_ENSEMBLE_SIZE = 100

#: Minimum number of conformers generated per *requested* conformation. CHOICE. Requesting
#: 50 conformations must not leave selection with 50 candidates and no choice.
MIN_OVERSAMPLE = 4

#: Rigid-body orientations tried per conformer while looking for a clash-free placement.
#:
#: CHOICE. The anchor construction leaves two free rotational degrees of freedom (the
#: plane of the tilt and the spin about the anchor axis), and both are sampled. 32 is
#: cheap -- each trial is one 3x3 matrix product and one KD-tree query.
ORIENTATION_ATTEMPTS = 32

#: Allowed *shortfall* of a STARLING virtual CA-CA bond below :data:`CA_CA_BOND_LENGTH`, in A.
#:
#: CHOICE, and deliberately much looser than :data:`CA_CA_BOND_TOLERANCE` (0.10 A), which
#: is DODO's tolerance for geometry DODO *builds*. STARLING reconstructs coordinates from
#: a predicted distance map, so its bonds scatter around 3.8 A by more than a build
#: tolerance while still describing a real chain. 0.5 A is loose enough not to reject the
#: model's normal output and tight enough to catch the two failures that matter: a
#: nanometre/Angstrom unit mixup (which shows up as bonds near 0.38) and a garbled or
#: transposed array (which shows up as anything else).
BOND_SCREEN_TOLERANCE = 0.5

#: Longest virtual CA-CA bond a screened conformer may contain, in Angstroms.
#:
#: DERIVED, and the reason the screen is asymmetric. :data:`dodo.constants.CA_CA_BOND_LENGTH`
#: documents the trans-peptide CA-CA distance as 3.80-3.81 A and "remarkably rigid": with
#: standard bond lengths and angles there is no way to make a *longer* one, so a virtual bond
#: above this ceiling has no all-atom realization at any dihedral. Applying
#: :data:`BOND_SCREEN_TOLERANCE` symmetrically accepted 3.30-4.30 A, so a conformer whose
#: every bond measured 4.250 A -- 11.8% long, 4.5x outside the build tolerance, and
#: unbuildable -- was screened as usable. The shortfall side stays loose because that is
#: where a distance-map reconstruction's noise lands and where cis peptides (~2.9 A) live,
#: and because it is the net that catches a nanometre/Angstrom unit mixup.
BOND_SCREEN_MAX_LENGTH = CA_CA_BOND_LENGTH + CA_CA_BOND_TOLERANCE

#: Pseudo-angle range a screened conformer must lie inside at **every** vertex.
#:
#: ``BACKBONE_ANGLE_OBSERVED_MIN..MAX`` (75-179 deg, MEASURED from real structures) rather
#: than the tighter ``BACKBONE_ANGLE_MIN..MAX`` generation window, because that window is
#: DODO's own sampling taste and STARLING did not agree to it.
#:
#: Not a fraction. This used to allow 2% of a conformer's vertices to sit at any angle
#: whatsoever, which is the wrong shape of gate for a hard physical constraint: a 50 degree
#: CA-CA-CA angle puts the i-1 and i+1 alpha carbons 3.2 A apart with a whole peptide unit
#: to fit between them, so no all-atom reconstruction of that vertex exists. A percentage
#: allowance also *grows with length* -- 7 arbitrarily sharp vertices in a 380-residue
#: region -- so the longer the region, the more impossible geometry it admitted. One
#: unreconstructable vertex is one too many, so the conformer is dropped and another of the
#: hundred in the ensemble is used instead.
SCREEN_ANGLE_WINDOW: tuple[float, float] = (
    BACKBONE_ANGLE_OBSERVED_MIN,
    BACKBONE_ANGLE_OBSERVED_MAX,
)

#: Relative and absolute tolerance on the achieved end-to-end distance of a built region.
#:
#: CHOICE. Deliberately not tighter: the target itself comes from ALBATROSS, whose own error
#: on Re is several percent, and from an analytical fallback that measures 0.6-0.95x of the
#: prediction for real IDR compositions. Demanding better than 5% of an uncertain number is
#: false precision. Deliberately not absent either, which is what it used to be: without a
#: tolerance a 350 A request was answered with a 104 A chain and reported as a success,
#: which silently discards every :data:`dodo.constants.MODES` multiplier.
#:
#: Shared with :mod:`dodo.engines.hierarchical`, which imports both names from here so that
#: the two engines cannot drift apart on what "the right dimension" means.
END_TO_END_RELATIVE_TOLERANCE = 0.05
END_TO_END_ABSOLUTE_TOLERANCE = 2.0

#: Attribute names probed, in order, to get coordinates out of whatever STARLING returns.
#:
#: ASSUMPTION about STARLING's API -- see the module report. Each is tried as a plain
#: attribute and then as a zero-argument method.
COORDINATE_ATTRIBUTES: tuple[str, ...] = ("coordinates", "coords", "ca_coords", "xyz")

#: Attribute names probed for STARLING's own maximum sequence length.
MAX_LENGTH_ATTRIBUTES: tuple[str, ...] = (
    "MAX_SEQUENCE_LENGTH",
    "MAX_LENGTH",
    "MAXIMUM_SEQUENCE_LENGTH",
    "STARLING_MAX_LENGTH",
)

#: Attribute names probed to decide whether model weights are on disk.
WEIGHTS_PATH_ATTRIBUTES: tuple[str, ...] = (
    "DEFAULT_MODEL_WEIGHTS_PATH",
    "MODEL_WEIGHTS_PATH",
    "DEFAULT_WEIGHTS_PATH",
)


# ---------------------------------------------------------------------------
# Optional-dependency plumbing. Three seams, all cheap, all monkeypatchable.
# ---------------------------------------------------------------------------


def starling_installed() -> bool:
    """Return True if the ``starling`` package could be imported, without importing it.

    A spec lookup touches the filesystem and nothing else, so this stays microseconds
    even though importing STARLING costs seconds of torch startup. Never raises:
    :func:`importlib.util.find_spec` can raise for a package whose parent is broken, and
    "it is broken" and "it is absent" are the same answer to this question.
    """
    try:
        return importlib.util.find_spec(STARLING_MODULE) is not None
    except (ImportError, ValueError, AttributeError):  # pragma: no cover - install-specific
        return False


def _loaded_starling() -> ModuleType | None:
    """Return the ``starling`` module only if it is *already* imported, else None.

    Lets :meth:`StarlingEngine.available` do deeper checks for free when something else
    has already paid the import cost, without ever paying it itself.
    """
    return sys.modules.get(STARLING_MODULE)


def _import_starling() -> ModuleType:
    """Import and return the ``starling`` module. **The** import site for STARLING.

    Raises
    ------
    MissingDependencyError
        If STARLING is not installed. The message names the extra, because
        ``ModuleNotFoundError: starling`` does not tell a user that the install command
        spells the distribution differently from the import.
    """
    try:
        module = importlib.import_module(STARLING_MODULE)
    except ImportError as exc:
        raise MissingDependencyError(
            package=STARLING_DISTRIBUTION,
            purpose="Generating IDR conformers with the STARLING engine",
            extra="starling",
        ) from exc
    return module


def _generate_callable(module: ModuleType) -> Callable[..., Any]:
    """Find STARLING's ensemble-generation entry point on an imported module.

    Probes the documented top-level ``starling.generate`` first, then the submodule path
    it is re-exported from. Raises rather than guessing further: silently calling the
    wrong function would produce coordinates nobody could account for.
    """
    candidate = getattr(module, "generate", None)
    if callable(candidate):
        return candidate  # type: ignore[no-any-return]
    frontend = getattr(module, "frontend", None)
    nested = getattr(getattr(frontend, "ensemble_generation", None), "generate", None)
    if callable(nested):
        return nested  # type: ignore[no-any-return]
    raise EngineUnavailableError(
        f"The installed {STARLING_MODULE!r} package exposes neither "
        f"{STARLING_MODULE}.generate nor "
        f"{STARLING_MODULE}.frontend.ensemble_generation.generate, so DODO cannot call "
        f"it. This is an API mismatch, not a missing install; check the STARLING version."
    )


def _weights_present(module: ModuleType) -> bool | None:
    """Whether STARLING's model weights are on disk.

    Returns ``None`` for "cannot tell", which is a genuinely different answer from False
    and must not be collapsed into it: STARLING may download weights on first use, and
    refusing to run because DODO could not find a file it does not own would be worse
    than trying and reporting the failure.
    """
    for holder in (module, getattr(module, "configs", None)):
        if holder is None:
            continue
        for attribute in WEIGHTS_PATH_ATTRIBUTES:
            value = getattr(holder, attribute, None)
            if isinstance(value, str | os.PathLike):
                return os.path.exists(value)
    return None


@dataclass(frozen=True, slots=True)
class LengthCap:
    """A maximum sequence length and where the number came from.

    ``source`` exists so a run can state whether the cap was read from the installed
    STARLING or from DODO's fallback constant. The two can disagree after a STARLING
    release, and a build that silently used a stale 380 when the installed model accepts
    more would quietly split IDRs it did not need to split.
    """

    value: int
    source: str

    def __str__(self) -> str:
        return f"{self.value} residues (from {self.source})"


def starling_max_length(*, probe: bool = True) -> LengthCap:
    """Return STARLING's maximum sequence length, read from STARLING when possible.

    Parameters
    ----------
    probe
        Import STARLING to ask it, when it is installed. Set False to force DODO's
        constant, which is what a caller that must not pay a torch import should do.

    Notes
    -----
    When STARLING is not installed this costs one spec lookup and returns
    :data:`dodo.constants.STARLING_MAX_LENGTH`. It never guesses: an installed STARLING
    that exposes no length attribute also yields the constant, with a ``source`` saying so.
    """
    fallback = LengthCap(STARLING_MAX_LENGTH, "dodo.constants.STARLING_MAX_LENGTH")
    module = _loaded_starling()
    if module is None:
        if not probe or not starling_installed():
            return fallback
        try:
            module = _import_starling()
        except MissingDependencyError:  # pragma: no cover - race with an uninstall
            return fallback
    for holder_name, holder in (("", module), ("configs.", getattr(module, "configs", None))):
        if holder is None:
            continue
        for attribute in MAX_LENGTH_ATTRIBUTES:
            value = getattr(holder, attribute, None)
            if isinstance(value, int | np.integer) and not isinstance(value, bool) and value > 0:
                return LengthCap(int(value), f"{STARLING_MODULE}.{holder_name}{attribute}")
    return fallback


# ---------------------------------------------------------------------------
# Conformer selection
# ---------------------------------------------------------------------------


def end_to_end_tolerance(
    target: float,
    override: float | None = None,
    *,
    relative: float = END_TO_END_RELATIVE_TOLERANCE,
    absolute: float = END_TO_END_ABSOLUTE_TOLERANCE,
) -> float:
    """Absolute tolerance, in Angstroms, on an achieved end-to-end distance.

    One definition shared by :class:`StarlingEngine` and
    :class:`~dodo.engines.hierarchical.HierarchicalEngine` so that "close enough to the
    target" cannot mean two different things one level apart.
    """
    if override is not None:
        return _validated_tolerance(override)
    return max(float(absolute), float(relative) * abs(float(target)))


def _validated_tolerance(value: float) -> float:
    """Coerce a caller-supplied end-to-end tolerance, rejecting NaN and negatives.

    A NaN tolerance compares False against everything, so it would *disable* the dimension
    check rather than loosen it -- the same silent-half-a-gate failure mode
    :func:`dodo.geometry.metrics.validate_ca_trace` now raises on.
    """
    tolerance = float(value)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise EngineError(f"end_to_end_tolerance must be finite and non-negative, got {value!r}.")
    return tolerance


def desired_internal_span(
    target_end_to_end: float,
    n_anchor_xyz: np.ndarray | None,
    c_anchor_xyz: np.ndarray | None,
    *,
    bond_length: float = CA_CA_BOND_LENGTH,
    slack: float = 0.0,
) -> tuple[float, tuple[float, float] | None]:
    """Resolve what end-to-end distance the generated chain itself should have.

    The target from :mod:`dodo.construct.dimensions` describes the *region*, but an
    interior region also has to fit between two fixed anchors, and the anchors win: they
    are real coordinates in the input structure, while the target is a prediction.

    Geometry. With anchors ``A`` and ``B`` at separation ``D``, a chain whose first CA is
    one bond from ``A`` and whose last CA is one bond from ``B`` can have any internal
    end-to-end distance in ``[D - 2b, D + 2b]`` and no other -- the two bond spheres of
    radius ``b`` bound the reachable separation of the two termini. So the desired value
    is the requested target clipped into that window, and the window itself is returned so
    a caller can report how far the anchors forced it from the prediction.

    Parameters
    ----------
    slack
        Angstroms of room to leave on **both** sides of the returned aim point, so that a
        builder which lands within ``slack`` of what it aimed at is still inside the
        anchor-feasible window.

        Why it exists: clipping onto the bare window put the aim point exactly on an edge,
        with zero room on that side. A builder is then asked to hit a value it may only
        approach from the infeasible direction -- a measured 0.163% undershoot of a 292 A
        target was enough to make both anchors unreachable and throw the whole conformation
        away, which is what made 35% of interior builds return nothing. Pass the builder's
        own end-to-end tolerance. Slack wider than half the window collapses the aim point
        onto the window centre, which is the anchor separation itself: the point with the
        most room on both sides and the one the anchors most directly demand.

    Returns
    -------
    tuple[float, tuple[float, float] | None]
        The desired internal end-to-end distance, and the *full* feasible window (not the
        slack-shrunk one, so a caller reports feasibility against real geometry), or
        ``None`` for a region with fewer than two anchors, where no window applies.
    """
    target = float(target_end_to_end)
    if n_anchor_xyz is None or c_anchor_xyz is None:
        return target, None
    separation = float(np.linalg.norm(np.asarray(c_anchor_xyz) - np.asarray(n_anchor_xyz)))
    low = max(0.0, separation - 2.0 * bond_length)
    high = separation + 2.0 * bond_length
    margin = float(np.clip(slack, 0.0, 0.5 * (high - low)))
    return float(np.clip(target, low + margin, high - margin)), (low, high)


def rank_conformers(
    conformers: np.ndarray,
    desired_end_to_end: float,
    *,
    window: tuple[float, float] | None = None,
) -> np.ndarray:
    """Order conformers by how well their end-to-end distance suits the request.

    Conformers inside ``window`` (the anchor-feasible band from
    :func:`desired_internal_span`) always precede conformers outside it, and within each
    group the ordering is by distance from ``desired_end_to_end``. A caller that walks
    this order and stops at the first placement that works therefore tries every feasible
    conformer before any infeasible one, and reports an honest residual if it has to fall
    back to one that cannot reach.
    """
    array = _as_conformers(conformers)
    spans = np.linalg.norm(array[:, -1, :] - array[:, 0, :], axis=1)
    error = np.abs(spans - float(desired_end_to_end))
    if window is None:
        return np.argsort(error, kind="stable")
    low, high = window
    outside = np.maximum(np.maximum(low - spans, spans - high), 0.0)
    # Lexsort with the last key as primary: feasibility first, then closeness.
    return np.lexsort((error, outside > 0.0))


# ---------------------------------------------------------------------------
# Rigid-body placement between anchors
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AnchorPlacement:
    """One conformer, placed, with every measurement needed to judge the placement.

    ``ca_coords`` is always a real conformation -- the same chain the generator produced,
    rotated and translated. It is never a fabricated or zero-filled array, even when
    ``ok`` is False: a caller that wants to show the user *why* a placement failed needs
    the geometry that failed. Whether it is usable is :attr:`ok`, and nothing else.

    Attributes
    ----------
    ca_coords
        ``(n_residues, 3)`` placed CA coordinates.
    conformer_index
        Index of the conformer chosen from the ensemble handed in.
    achieved_end_to_end
        Measured CA(first)-to-CA(last) distance of the placed chain, in Angstroms.
    desired_end_to_end
        What :func:`desired_internal_span` asked for.
    n_anchor_gap, c_anchor_gap
        Distance from the first/last generated CA to its anchor, in Angstroms, or ``None``
        where there is no anchor on that side. These should equal
        :data:`dodo.constants.CA_CA_BOND_LENGTH`; the residual is how far off they are.
    anchor_residual
        Worst ``|gap - bond_length|`` over the anchors present. ``0.0`` when free.
    min_internal_ca_distance
        Closest non-bonded CA-CA approach *within* the conformer. A rigid-body placement
        cannot change a conformer's own internal geometry, so this is measured rather than
        fixed -- but it is a hard gate on :attr:`ok`, not merely reported.
    internal_clash_cutoff
        The floor :attr:`min_internal_ca_distance` is held to, in Angstroms. Defaults to
        :data:`dodo.constants.CA_CLASH_DISTANCE`; a caller that has *declared* a relaxed
        threshold on its result may lower it to that rung.
    clash_free
        Whether the placed chain clears the obstacles at :attr:`relaxed_to`.
    relaxed_to
        Clash threshold actually used, in Angstroms, or ``None`` if no relaxation was
        needed. Relaxation is recorded, never hidden.
    orientations_tried
        How many rigid-body orientations were evaluated.
    notes
        Human-readable provenance, including every compromise made.
    """

    ca_coords: np.ndarray
    conformer_index: int
    achieved_end_to_end: float
    desired_end_to_end: float
    n_anchor_gap: float | None
    c_anchor_gap: float | None
    anchor_residual: float
    min_internal_ca_distance: float
    clash_free: bool
    relaxed_to: float | None
    orientations_tried: int
    notes: tuple[str, ...] = ()
    internal_clash_cutoff: float = CA_CLASH_DISTANCE

    @property
    def ok(self) -> bool:
        """True when both anchors are met and nothing clashes -- including with itself.

        ``min_internal_ca_distance`` is part of this and used not to be. A conformer that
        folds back through itself has perfectly good bonds and angles, because both are
        purely local measurements, so nothing else here notices it: chains with non-bonded
        CA pairs 2.1 A apart -- two carbon atoms inside each other -- were returned as
        usable placements because they cleared the *external* obstacles.
        """
        return (
            self.clash_free
            and self.anchor_residual <= CA_CA_BOND_TOLERANCE
            and self.min_internal_ca_distance >= self.internal_clash_cutoff
        )


def place_between_anchors(
    conformer: np.ndarray,
    *,
    n_anchor_xyz: np.ndarray | None,
    c_anchor_xyz: np.ndarray | None,
    rng: np.random.Generator,
    obstacles: np.ndarray | None = None,
    conformer_index: int = 0,
    desired_end_to_end: float | None = None,
    bond_length: float = CA_CA_BOND_LENGTH,
    clash_cutoff: float = CA_CLASH_DISTANCE,
    internal_clash_cutoff: float = CA_CLASH_DISTANCE,
    relaxation_ladder: Sequence[float] = CLASH_RELAXATION_LADDER,
    orientation_attempts: int = ORIENTATION_ATTEMPTS,
) -> AnchorPlacement:
    """Place a generated chain as a rigid body so its termini meet their anchors.

    The construction, for the two-anchor case, is exact rather than iterative. Let ``A``
    and ``B`` be the anchor CAs, ``D = |B - A|``, ``u`` the unit vector from ``A`` to
    ``B``, ``v`` any unit vector perpendicular to ``u``, and ``b`` the CA-CA bond length.
    Put the first generated CA at ``X = A + b(cos t u + sin t v)`` and the last at
    ``Y = B - b(cos t u - sin t v)``. Then ``|X - A| = |Y - B| = b`` for every ``t``, and

        ``Y - X = (D - 2b cos t) u``

    -- the perpendicular parts cancel by construction. So the separation of the two
    termini is ``D - 2b cos t``, which is exactly the conformer's own end-to-end distance
    ``d`` when ``cos t = (D - d) / 2b``. That has a solution whenever ``|D - d| <= 2b``,
    which is the feasibility window :func:`desired_internal_span` reports. Rotating the
    conformer's ``last - first`` vector onto ``Y - X`` and translating its first CA to
    ``X`` then satisfies both anchors to machine precision, leaving two free rotations --
    the choice of ``v`` and a spin about the ``X``-``Y`` axis -- to spend on avoiding
    clashes.

    Why this matters: the pre-rewrite builder wrote the first and last generated residues
    directly on top of the anchor CAs, producing 0.00 A coincident atoms at every junction
    of every output structure. Here the junction distance is a bond length by construction
    and is measured and reported afterwards regardless.

    For a single anchor there is no closure constraint: the chain is placed one bond from
    its anchor and oriented so it grows away from it. With no anchors it is centred on the
    origin and randomly oriented.

    Parameters
    ----------
    conformer
        ``(n_residues, 3)`` CA coordinates in their own frame.
    n_anchor_xyz, c_anchor_xyz
        Fixed anchor CA coordinates, or ``None`` at a chain terminus.
    rng
        Seeded generator. Every orientation tried is drawn from it, so a placement is
        reproducible.
    obstacles
        ``(m, 3)`` coordinates that the placed chain must clear. Pass an obstacle set that
        already excludes the anchor residues' own atoms: they sit one bond length from the
        chain's termini by design, and counting that as a clash would reject every correct
        placement.
    desired_end_to_end
        What the selection stage asked for, carried through for reporting only.
    clash_cutoff
        Minimum approach to an *obstacle*, in Angstroms.
    internal_clash_cutoff
        Minimum non-bonded approach the conformer must already have *within itself* for the
        placement to count as usable. No rotation can improve it, so this is a property of
        the conformer that the placement merely inherits -- and refuses to certify.
    relaxation_ladder
        Clash thresholds tried in order when the requested one admits no orientation. Only
        rungs at or below ``clash_cutoff`` are used.

    Returns
    -------
    AnchorPlacement
        The best placement found, with :attr:`AnchorPlacement.ok` saying whether it is
        usable. This function does not raise on failure to satisfy the anchors -- an
        interior IDR whose anchors are unreachable for every conformer is a real situation
        the caller must be able to report per conformer.

    Raises
    ------
    GeometryError
        If ``conformer`` is not ``(n, 3)`` with ``n >= 1``, or contains non-finite values.
    """
    chain = _as_single_conformer(conformer)
    n_residues = chain.shape[0]
    ladder = _clash_ladder(clash_cutoff, relaxation_ladder)
    internal_min = _min_internal_ca_distance(chain)
    span = float(np.linalg.norm(chain[-1] - chain[0])) if n_residues > 1 else 0.0
    desired = span if desired_end_to_end is None else float(desired_end_to_end)
    notes: list[str] = []

    anchor_n = None if n_anchor_xyz is None else np.asarray(n_anchor_xyz, dtype=np.float64)
    anchor_c = None if c_anchor_xyz is None else np.asarray(c_anchor_xyz, dtype=np.float64)

    if n_residues == 1:
        # A one-residue region has no orientation and no internal span. It is still a
        # real request (a single missing residue), and it still must not be written on
        # top of an anchor, so place it one bond length away along the anchor axis.
        return _place_single_residue(
            anchor_n,
            anchor_c,
            rng,
            obstacles,
            ladder,
            bond_length,
            desired,
            internal_min,
            internal_clash_cutoff,
        )

    candidates: list[tuple[np.ndarray, float]] = []
    tried = 0
    for _ in range(max(1, int(orientation_attempts))):
        tried += 1
        placed, orientation_notes = _one_orientation(chain, anchor_n, anchor_c, rng, bond_length)
        for note in orientation_notes:
            if note not in notes:
                notes.append(note)
        residual = _anchor_residual(placed, anchor_n, anchor_c, bond_length)
        # Return as soon as an orientation is good on both counts. Generating the whole
        # budget first would cost an unobstructed placement -- the common case -- 32
        # rotations and 32 KD-tree queries for nothing.
        if residual <= CA_CA_BOND_TOLERANCE and _clash_free(placed, obstacles, ladder[0]):
            return _placement(
                placed,
                conformer_index,
                desired,
                anchor_n,
                anchor_c,
                bond_length,
                internal_min,
                clash_free=True,
                relaxed_to=None,
                internal_clash_cutoff=internal_clash_cutoff,
                tried=tried,
                notes=notes,
            )
        candidates.append((placed, residual))

    # Nothing was clean. Prefer orientations that at least satisfy the anchors, then walk
    # the relaxation ladder. Relaxing the clash distance only after exhausting every
    # orientation at the strict threshold is deliberate: it is a real concession, so it is
    # the last resort and it is recorded.
    order = sorted(range(len(candidates)), key=lambda i: candidates[i][1])
    for threshold in ladder:
        for index in order:
            placed, residual = candidates[index]
            if _clash_free(placed, obstacles, threshold):
                relaxed = None if threshold == ladder[0] else threshold
                if relaxed is not None:
                    notes.append(
                        f"No orientation cleared obstacles at {ladder[0]:.2f} A; accepted "
                        f"one at a relaxed {threshold:.2f} A."
                    )
                return _placement(
                    placed,
                    conformer_index,
                    desired,
                    anchor_n,
                    anchor_c,
                    bond_length,
                    internal_min,
                    clash_free=True,
                    relaxed_to=relaxed,
                    internal_clash_cutoff=internal_clash_cutoff,
                    tried=tried,
                    notes=notes,
                )

    best, _ = candidates[order[0]]
    notes.append(
        f"No orientation of this conformer cleared the obstacles at any threshold down to "
        f"{ladder[-1]:.2f} A after {tried} attempts; returning the least-bad placement "
        f"with clash_free=False."
    )
    return _placement(
        best,
        conformer_index,
        desired,
        anchor_n,
        anchor_c,
        bond_length,
        internal_min,
        clash_free=False,
        relaxed_to=None,
        internal_clash_cutoff=internal_clash_cutoff,
        tried=tried,
        notes=notes,
    )


def _one_orientation(
    chain: np.ndarray,
    anchor_n: np.ndarray | None,
    anchor_c: np.ndarray | None,
    rng: np.random.Generator,
    bond_length: float,
) -> tuple[np.ndarray, list[str]]:
    """Produce one candidate rigid-body placement of ``chain``. See the caller's docstring.

    Returns the placed coordinates and every note the placement earned. A list rather than
    one note because two things can be wrong at once -- coincident anchors *and* a conformer
    that cannot bridge them -- and reporting only the last one loses the first.
    """
    span_vector = chain[-1] - chain[0]
    span = float(np.linalg.norm(span_vector))
    notes: list[str] = []

    if anchor_n is not None and anchor_c is not None:
        separation_vector = anchor_c - anchor_n
        separation = float(np.linalg.norm(separation_vector))
        if separation < _DEGENERATE_LENGTH:
            # Coincident anchors: no axis to build the construction on. Any direction is
            # as good as any other, so pick one at random and report the situation.
            axis = random_unit_vectors(1, rng)[0]
            notes.append(
                "The two anchors are coincident, which no real structure should contain; "
                "the placement axis was chosen at random."
            )
        else:
            axis = separation_vector / separation
        cosine = np.clip((separation - span) / (2.0 * bond_length), -1.0, 1.0)
        exact = abs((separation - span) / (2.0 * bond_length)) <= 1.0
        if not exact:
            notes.append(
                f"This conformer's end-to-end distance ({span:.1f} A) cannot bridge an "
                f"anchor separation of {separation:.1f} A: the reachable window is "
                f"[{max(0.0, separation - 2 * bond_length):.1f}, "
                f"{separation + 2 * bond_length:.1f}] A. Placed to satisfy the N-anchor "
                f"exactly; the C-anchor gap is reported as a residual."
            )
        sine = float(np.sqrt(max(0.0, 1.0 - cosine**2)))
        perpendicular = _perpendicular(axis, rng)
        start = anchor_n + bond_length * (cosine * axis + sine * perpendicular)
        end = anchor_c - bond_length * (cosine * axis - sine * perpendicular)
        return _rigid_place(chain, start, end, rng), notes

    if anchor_n is not None:
        direction = random_unit_vectors(1, rng)[0]
        start = anchor_n + bond_length * direction
        return _rigid_place(chain, start, start + span * direction, rng), notes

    if anchor_c is not None:
        direction = random_unit_vectors(1, rng)[0]
        end = anchor_c + bond_length * direction
        # The chain runs backwards from the C-anchor: its LAST residue is the one that
        # must sit a bond length away, so orient first-from-last along the outward
        # direction and place the last residue at `end`.
        return _rigid_place(chain, end + span * direction, end, rng), notes

    centred = chain - chain.mean(axis=0)
    rotation = rotation_from_axis_angle(
        random_unit_vectors(1, rng)[0], float(rng.uniform(0.0, 2.0 * np.pi))
    )
    return apply(centred, rotation), notes


def _rigid_place(
    chain: np.ndarray, start: np.ndarray, end: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    """Rotate and translate ``chain`` so its first CA is at ``start`` and its last at ``end``.

    ``|end - start|`` must equal the chain's own end-to-end distance for both to be hit
    exactly; when it does not (an infeasible anchor pair), the first CA lands on ``start``
    and the chain points at ``end``, which is the best a rigid body can do and is why the
    caller reports a residual.
    """
    span_vector = chain[-1] - chain[0]
    target_vector = end - start
    if (
        float(np.linalg.norm(span_vector)) < _DEGENERATE_LENGTH
        or float(np.linalg.norm(target_vector)) < _DEGENERATE_LENGTH
    ):
        # Either the conformer's termini coincide or the two placement points do. Neither
        # defines a direction, so there is nothing to align: use a random orientation.
        # Without this guard rotation_between_vectors would be handed a zero vector, and
        # the pre-rewrite code's response to that was a matrix of NaN.
        rotation = rotation_from_axis_angle(
            random_unit_vectors(1, rng)[0], float(rng.uniform(0.0, 2.0 * np.pi))
        )
        randomized: np.ndarray = apply(chain - chain[0], rotation) + start
        return randomized
    rotation = rotation_between_vectors(span_vector, target_vector)
    aligned = apply(chain - chain[0], rotation) + start
    # The spin about the start-end axis is a free degree of freedom: it changes nothing
    # about the anchors and everything about which obstacles the chain runs into.
    spin = rotation_from_axis_angle(target_vector, float(rng.uniform(0.0, 2.0 * np.pi)))
    spun: np.ndarray = apply(aligned - start, spin) + start
    return spun


def _place_single_residue(
    anchor_n: np.ndarray | None,
    anchor_c: np.ndarray | None,
    rng: np.random.Generator,
    obstacles: np.ndarray | None,
    ladder: tuple[float, ...],
    bond_length: float,
    desired: float,
    internal_min: float,
    internal_clash_cutoff: float,
) -> AnchorPlacement:
    """Place a one-residue region: no orientation, but still not on top of an anchor."""
    notes: list[str] = []
    if anchor_n is not None and anchor_c is not None:
        midpoint = 0.5 * (anchor_n + anchor_c)
        separation = float(np.linalg.norm(anchor_c - anchor_n))
        if separation < _DEGENERATE_LENGTH:
            offset = bond_length * random_unit_vectors(1, rng)[0]
        else:
            # Sit off the anchor-anchor axis by whatever it takes to be one bond from
            # each: half the separation along the axis is fixed, so the perpendicular
            # offset is the remaining leg of the isoceles triangle.
            half = 0.5 * separation
            perpendicular_length = float(np.sqrt(max(0.0, bond_length**2 - half**2)))
            offset = perpendicular_length * _perpendicular((anchor_c - anchor_n) / separation, rng)
        placed = (midpoint + offset).reshape(1, 3)
    elif anchor_n is not None:
        placed = (anchor_n + bond_length * random_unit_vectors(1, rng)[0]).reshape(1, 3)
    elif anchor_c is not None:
        placed = (anchor_c + bond_length * random_unit_vectors(1, rng)[0]).reshape(1, 3)
    else:
        placed = np.zeros((1, 3), dtype=np.float64)
    threshold = next((rung for rung in ladder if _clash_free(placed, obstacles, rung)), None)
    if threshold is None:
        notes.append("A single-residue region could not be placed clear of the obstacles.")
    return _placement(
        placed,
        0,
        desired,
        anchor_n,
        anchor_c,
        bond_length,
        internal_min,
        clash_free=threshold is not None,
        relaxed_to=None if threshold in (None, ladder[0]) else threshold,
        internal_clash_cutoff=internal_clash_cutoff,
        tried=1,
        notes=notes,
    )


def _placement(
    coords: np.ndarray,
    conformer_index: int,
    desired: float,
    anchor_n: np.ndarray | None,
    anchor_c: np.ndarray | None,
    bond_length: float,
    internal_min: float,
    *,
    clash_free: bool,
    relaxed_to: float | None,
    internal_clash_cutoff: float,
    tried: int,
    notes: Sequence[str],
) -> AnchorPlacement:
    """Measure a finished placement and package it up."""
    n_gap = None if anchor_n is None else float(np.linalg.norm(coords[0] - anchor_n))
    c_gap = None if anchor_c is None else float(np.linalg.norm(coords[-1] - anchor_c))
    achieved = float(np.linalg.norm(coords[-1] - coords[0])) if coords.shape[0] > 1 else 0.0
    return AnchorPlacement(
        ca_coords=coords,
        conformer_index=conformer_index,
        achieved_end_to_end=achieved,
        desired_end_to_end=desired,
        n_anchor_gap=n_gap,
        c_anchor_gap=c_gap,
        anchor_residual=_anchor_residual(coords, anchor_n, anchor_c, bond_length),
        min_internal_ca_distance=internal_min,
        clash_free=clash_free,
        relaxed_to=relaxed_to,
        orientations_tried=tried,
        notes=tuple(notes),
        internal_clash_cutoff=internal_clash_cutoff,
    )


def _dimension_mask(
    achieved: np.ndarray,
    *,
    desired: float,
    tolerance: float,
    window: tuple[float, float] | None,
    n_residues: int,
) -> tuple[np.ndarray, str | None]:
    """Which conformers have the end-to-end distance the request asked for, and why not.

    Three regimes, because three different things are ground truth.

    * **A window** (both anchors fixed): the anchors are real coordinates in the input
      structure and the target is a prediction, so what matters is landing inside the window
      the anchors permit. Also demanding closeness to the clipped prediction would fail
      conformers that satisfy both anchors exactly, which is the one thing that cannot be
      compromised.
    * **One conformation, no window**: the target is the only statement about the dimension
      and this chain is the whole answer, so it must be within ``tolerance`` of it.
    * **Several conformations, no window**: a predicted end-to-end distance is the *mean* of
      a distribution, not a per-conformer constraint (see
      :class:`~dodo.engines.base.IDRRequest`). A batch containing a 20 A and a 50 A chain
      around a 35 A target is a *correct* ensemble, so the gate is on the batch mean. When
      the mean is wrong the batch as a whole has the wrong dimension and every conformer in
      it is reported as a failure -- there is no subset of it that is the right answer.

    A region of fewer than two residues has no end-to-end distance and is exempt.

    Returns
    -------
    tuple[np.ndarray, str | None]
        Per-conformer mask, and a note naming both numbers when something failed.
    """
    n = int(achieved.shape[0])
    if n_residues < 2:
        return np.ones(n, dtype=bool), None
    if window is not None:
        inside = (achieved >= window[0]) & (achieved <= window[1])
        if bool(inside.all()):
            return inside, None
        worst = float(achieved[~inside][0])
        return inside, (
            f"{int(np.count_nonzero(~inside))} of {n} conformer(s) came out with an "
            f"end-to-end distance the anchors cannot accommodate (e.g. {worst:.1f} A against "
            f"a feasible window of [{window[0]:.1f}, {window[1]:.1f}] A), so they are "
            f"reported as failures rather than returned with unreachable anchors."
        )
    if n == 1:
        error = abs(float(achieved[0]) - desired)
        if error <= tolerance:
            return np.ones(1, dtype=bool), None
        return np.zeros(1, dtype=bool), (
            f"The ensemble's best conformer is {float(achieved[0]):.1f} A end-to-end against "
            f"a target of {desired:.1f} +/- {tolerance:.1f} A (off by {error:.1f} A), so no "
            f"conformer of the requested dimension exists in it. Reported as a failure "
            f"rather than returned at the wrong size."
        )
    mean = float(achieved.mean())
    error = abs(mean - desired)
    if error <= tolerance:
        return np.ones(n, dtype=bool), None
    return np.zeros(n, dtype=bool), (
        f"The {n} selected conformers average {mean:.1f} A end-to-end against a target of "
        f"{desired:.1f} +/- {tolerance:.1f} A (off by {error:.1f} A). A target is the mean of "
        f"a distribution, so the batch is judged on its mean -- and this batch is the wrong "
        f"size, so none of it is reported as a success."
    )


def _anchor_residual(
    coords: np.ndarray,
    anchor_n: np.ndarray | None,
    anchor_c: np.ndarray | None,
    bond_length: float,
) -> float:
    """Worst deviation of a terminal-to-anchor distance from one bond length."""
    residual = 0.0
    if anchor_n is not None:
        residual = max(residual, abs(float(np.linalg.norm(coords[0] - anchor_n)) - bond_length))
    if anchor_c is not None:
        residual = max(residual, abs(float(np.linalg.norm(coords[-1] - anchor_c)) - bond_length))
    return residual


# ---------------------------------------------------------------------------
# Screening: what came back from the model, measured
# ---------------------------------------------------------------------------


def regularize_conformers(
    conformers: np.ndarray,
    *,
    bond_length: float = CA_CA_BOND_LENGTH,
    angle_window: tuple[float, float] | None = SCREEN_ANGLE_WINDOW,
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Repair every conformer's virtual bonds *and* pseudo-angles into physical geometry.

    STARLING is a diffusion model that reconstructs coordinates from a predicted distance
    map, so neither its bonds nor its angles are a backbone's. That is normal model output,
    not a defect -- but it is not a protein either: the trans-peptide CA-CA distance is rigid,
    DODO writes structures whose bonds a viewer will draw, and the validator checks them
    against 3.81 +/- 0.10 A.

    So the geometry is *corrected*, not merely screened. Screening alone -- which is all this
    module used to do -- either rejects usable conformers for ordinary reconstruction error or
    passes that error straight through into the output file.

    Bonds were the easy half and repairing only them was not enough. MEASURED over 20
    conformers of a 201-residue sequence, after every bond was projected onto 3.81 A exactly the
    pseudo-angles still spanned 10.5-179.3 deg, 2.34% of vertices were outside
    :data:`SCREEN_ANGLE_WINDOW` and the median conformer had 5 of them, so
    :func:`screen_conformers` -- which rejects a conformer for a single unreconstructable vertex,
    and rightly -- accepted 0 of 20. Both constraints are therefore projected together; see
    :mod:`dodo.geometry.regularize` for why doing so is more than running two loops in sequence.

    Endpoints are deliberately *not* pinned: the conformer has not been placed against its
    anchors yet, so there is nothing to preserve, and letting the whole chain relax gives a
    better-conditioned projection than holding two ends fixed.

    A conformer the projection cannot repair is returned in whatever state it reached, with a
    note, and is *not* forced. Nothing extra is needed to drop it: its angles are still outside
    the window, so :func:`screen_conformers` rejects it downstream on exactly the criterion it
    failed. That is the honest path and the reason this function does not raise.

    Parameters
    ----------
    angle_window
        Pseudo-angle window to repair into, in degrees. Defaults to
        :data:`SCREEN_ANGLE_WINDOW`, i.e. the same window the screen enforces, so that repair
        and acceptance cannot drift apart. ``None`` repairs bonds only, which is what the
        engine did before angle repair existed and is kept for diagnosing the two halves
        separately.

    Returns
    -------
    tuple[numpy.ndarray, tuple[str, ...]]
        The corrected conformers, same shape as the input, and notes describing what moved.
    """
    from ..geometry.regularize import regularize_ca_trace

    if conformers.ndim != 3 or conformers.shape[2] != 3:
        raise EngineError(
            f"Expected conformers shaped (n_conformers, n_residues, 3), got {conformers.shape}."
        )
    if conformers.shape[1] < 2:
        return conformers, ()

    before = np.linalg.norm(np.diff(conformers, axis=1), axis=2)
    corrected = np.empty_like(conformers)
    unconverged = 0
    worst_angle_before = 0.0
    worst_angle_after = 0.0
    re_change: list[float] = []
    rg_change: list[float] = []
    for index, conformer in enumerate(conformers):
        # raise_on_failure=False, and this is the point of requesting angles at all. A conformer
        # whose bonds and angles are not jointly satisfiable is a real thing an ensemble can
        # contain, and the ensemble has ninety-nine others: raising would throw all of them away
        # over one, which is what this call used to do -- it checked `converged` while leaving
        # the default in place, so the branch below was unreachable.
        result = regularize_ca_trace(
            conformer,
            bond_length=bond_length,
            angle_window=angle_window,
            raise_on_failure=False,
        )
        corrected[index] = result.ca_coords
        if not result.converged:
            unconverged += 1
        worst_angle_before = max(worst_angle_before, result.max_angle_violation_before)
        worst_angle_after = max(worst_angle_after, result.max_angle_violation_after)
        re_change.append(result.end_to_end_change_fraction)
        rg_change.append(result.rg_change_fraction)
    after = np.linalg.norm(np.diff(corrected, axis=1), axis=2)
    displacement = np.linalg.norm(corrected - conformers, axis=2)

    notes = [
        f"Regularized {conformers.shape[0]} conformer(s) onto a {bond_length:.2f} A CA-CA "
        f"bond: worst bond deviation {np.abs(before - bond_length).max():.3f} A before, "
        f"{np.abs(after - bond_length).max():.3f} A after; atom displacement "
        f"{np.median(displacement):.3f} A median, {np.percentile(displacement, 95):.3f} A p95, "
        f"{displacement.max():.3f} A max."
    ]
    if angle_window is not None:
        notes.append(
            f"Repaired CA-CA-CA pseudo-angles into {angle_window[0]:.0f}-{angle_window[1]:.0f} "
            f"deg: worst excursion outside it {worst_angle_before:.2f} deg before, "
            f"{worst_angle_after:.2f} deg after."
        )
    # The measurement that says the repair is a projection and not a reshaping. Mean over the
    # ensemble, because a per-conformer worst case is dominated by whichever chain STARLING
    # broke most and says nothing about whether selection still sees the same distribution.
    notes.append(
        f"Ensemble dimensions moved by {float(np.mean(re_change)) * 100:+.2f}% in end-to-end "
        f"distance and {float(np.mean(rg_change)) * 100:+.2f}% in radius of gyration on average."
    )
    if unconverged:
        notes.append(
            f"{unconverged} of {conformers.shape[0]} conformer(s) could not be projected onto "
            f"both constraints at once within the sweep limit. They are left as they are rather "
            f"than forced, and the screen below rejects them on the criterion they still fail."
        )
    return corrected, tuple(notes)


#: Bond-length window a RAW STARLING conformer must fall inside before repair, in Angstroms.
#:
#: MEASURED over 40 conformers of a 38-residue sequence from STARLING 2.x: virtual CA-CA bonds
#: spanned 1.81-4.86 A with a median of 3.36, and 47.6% of all bonds were more than 0.5 A from
#: 3.81. That is not damage, it is simply what reconstructing coordinates from a predicted
#: distance map produces -- the MDS step fits the whole map and never constrains the bond.
#:
#: So a *physical* bond ceiling cannot be applied to raw output: at the trans-peptide limit of
#: 3.91 A it rejected 100 of 100 conformers and made the engine unusable. This window is
#: deliberately generous, and exists only to catch a chain the model genuinely lost -- the 6 A
#: discontinuity case -- before :func:`regularize_conformers` can quietly repair it. The
#: physical criteria are then applied AFTER repair, where they mean something.
RECONSTRUCTION_BOND_WINDOW: Final[tuple[float, float]] = (1.0, 6.0)

#: Window the MEDIAN virtual bond of a raw ensemble must fall in, in Angstroms.
#:
#: Tighter than :data:`RECONSTRUCTION_BOND_WINDOW`, and deliberately so: a median over thousands
#: of bonds is far more stable than any single bond, so it can be held to a narrower range than
#: the per-bond screen without rejecting real output. MEASURED medians from STARLING 2.x are
#: 3.36 A (38 residues) and 3.17 A (401 residues), comfortably inside this.
#:
#: The split matters. With one shared window at a 1.0 A floor, a synthetic "chain" whose bonds
#: were all exactly 1.0 A passed the unit check on the boundary and only failed later inside
#: bond repair, which reports a convergence failure rather than the real problem -- that the
#: coordinates are not a protein backbone in either unit.
RECONSTRUCTION_MEDIAN_WINDOW: Final[tuple[float, float]] = (2.0, 6.0)


#: Worst non-bonded CA-CA contact a raw conformer may have and still be worth repairing, in
#: Angstroms, expressed as a shortfall below :data:`~dodo.constants.CA_CLASH_DISTANCE`.
#:
#: MEASURED, and this is the cheap-elimination threshold that makes the pipeline affordable.
#: Repair is the expensive step -- 3.8 s for 100 conformers of 380 residues -- and it cannot fix an
#: internal clash: projecting bonds and angles rearranges a chain locally, it does not unthread it.
#: So a conformer that already collides badly should be dropped BEFORE paying to repair it.
#:
#: The threshold has to be loose enough that it never discards a conformer repair would have saved.
#: Over 60 conformers at 38 and 200 residues, every conformer that ultimately passed the physical
#: screen had a raw worst contact of at least 3.536 A, while the rejected ones sat at a median of
#: 2.218 A. Repair itself moves the worst contact by only -0.128 A at the median. A cutoff 0.3 A
#: below the 3.20 A clash distance -- so 2.90 A -- therefore has 0.6 A of margin against the
#: closest accepted conformer, and still eliminates the bulk of the hopeless ones for free.
PREFILTER_CLASH_SHORTFALL: Final[float] = 0.3

#: Largest distance any single alpha carbon may be moved by repair, in Angstroms.
#:
#: MEASURED. Repair is meant to keep a conformer as close to STARLING's output as it can while
#: making the geometry physical, and it normally does: displacement is 0.84 A at the median and
#: 2.10 A at the 99th percentile. But a severe angle violation near a free terminus has nothing to
#: anchor it, and one conformer in twelve had an atom travel 21.8 A to open a 9.9 degree vertex at
#: residue 3 -- a local rearrangement, not a repair. With ~100 conformers generated per round,
#: rejecting such a conformer costs nothing and keeps the promise that a kept conformer is
#: recognisably the one STARLING produced.
MAX_REPAIR_DISPLACEMENT: Final[float] = 4.0

#: Extra rounds of generation allowed when the first round does not yield enough usable conformers.
#:
#: CHOICE, bounded by cost: each round is a full diffusion-plus-MDS pass. Long regions are where
#: yield is poor -- at 200 residues only about 6 of 30 conformers survive, because STARLING's raw
#: output already contains internal clashes -- so one retry roughly doubles the pool for one extra
#: generation, and a hard cap keeps a pathological sequence from generating forever.
MAX_GENERATION_ROUNDS: Final[int] = 3


def _repair_one(chain: np.ndarray) -> tuple[np.ndarray, bool]:
    """Project one conformer onto physical bonds and angles. Returns it and whether to keep it.

    Rejects rather than returns a repair that moved any alpha carbon further than
    :data:`MAX_REPAIR_DISPLACEMENT`. Repair is supposed to keep a conformer recognisably the one
    STARLING produced while making its geometry valid, and it normally does -- displacement is
    0.84 A at the median and 2.10 A at the 99th percentile. But a severe angle violation beside a
    free terminus has nothing to anchor it: one conformer in twelve had an atom travel 21.8 A to
    open a 9.9 degree vertex at residue 3. That is a local rearrangement, not a repair, and with
    ~100 conformers generated per round it is cheaper to discard than to accept.
    """
    if chain.shape[0] < 3:
        # Nothing to project: a two-residue chain has one bond and no vertex, a one-residue chain
        # has neither. regularize_ca_trace rightly refuses to be handed either, so short regions
        # are placed exactly as STARLING produced them.
        return chain, True
    result = regularize_ca_trace(chain, angle_window=SCREEN_ANGLE_WINDOW, raise_on_failure=False)
    if not result.converged:
        return chain, False
    moved = float(np.max(np.linalg.norm(result.ca_coords - chain, axis=1)))
    if moved > MAX_REPAIR_DISPLACEMENT:
        return chain, False
    return result.ca_coords, True


def prefilter_conformers(conformers: np.ndarray) -> tuple[np.ndarray, tuple[str, ...]]:
    """Drop conformers not worth repairing, using only checks that cost nothing.

    This is the efficiency step. Repair is the expensive part of the pipeline -- 3.8 s for 100
    conformers of 380 residues -- and two defects survive it: a chain the model lost entirely, and
    a chain that collides with itself. Projecting bonds and angles rearranges a chain locally; it
    does not unthread one. So both are far cheaper to detect before paying to repair.

    The clash threshold has to be loose enough never to discard a conformer repair would have
    saved. MEASURED over 60 conformers at 38 and 200 residues: every conformer that ultimately
    passed :func:`screen_conformers` had a raw worst contact of at least 3.536 A, the rejected ones
    sat at a median of 2.218 A, and repair moves the worst contact by only -0.128 A at the median.
    See :data:`PREFILTER_CLASH_SHORTFALL`.

    Returns
    -------
    tuple[np.ndarray, tuple[str, ...]]
        Indices worth repairing, ordered **best first** -- fewest pseudo-angle violations, then the
        roomiest internal contact -- and notes recording what was dropped and why. The ordering is
        load-bearing: the caller repairs only as many as it needs, so the conformers likeliest to
        succeed have to be tried first or the saving evaporates.
    """
    array = _as_conformers(conformers)
    n = array.shape[0]
    low, high = RECONSTRUCTION_BOND_WINDOW
    bonds = np.linalg.norm(np.diff(array, axis=1), axis=2)
    intact = (bonds >= low).all(axis=1) & (bonds <= high).all(axis=1)

    floor = CA_CLASH_DISTANCE - PREFILTER_CLASH_SHORTFALL
    contacts = np.array([_min_internal_ca_distance(chain) for chain in array])
    roomy = contacts >= floor

    notes: list[str] = []
    lost = int(np.count_nonzero(~intact))
    if lost:
        notes.append(
            f"Dropped {lost} of {n} conformer(s) before repair: a virtual bond outside "
            f"{low}-{high} A is a lost chain, not reconstruction error."
        )
    hopeless = int(np.count_nonzero(intact & ~roomy))
    if hopeless:
        notes.append(
            f"Dropped {hopeless} of {n} conformer(s) before repair: their closest non-bonded "
            f"contact is under {floor:.2f} A, and repair cannot unthread a self-collision."
        )

    keep = np.flatnonzero(intact & roomy)
    if keep.size and array.shape[1] >= 3:
        window_low, window_high = SCREEN_ANGLE_WINDOW

        def out_of_window(index: int) -> int:
            angles = ca_pseudo_angles(array[index])
            return int(np.count_nonzero((angles < window_low) | (angles > window_high)))

        violations = np.array([out_of_window(int(index)) for index in keep])
        keep = keep[np.lexsort((-contacts[keep], violations))]
    return keep, tuple(notes)


def screen_conformers(
    conformers: np.ndarray,
    *,
    bond_tolerance: float = BOND_SCREEN_TOLERANCE,
    max_bond_length: float = BOND_SCREEN_MAX_LENGTH,
    angle_window: tuple[float, float] | None = SCREEN_ANGLE_WINDOW,
    internal_clash_cutoff: float = CA_CLASH_DISTANCE,
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Drop conformers whose geometry is not usable, and report what was measured.

    A generative model's output is data, not a promise. This measures every conformer
    against the ranges in :mod:`dodo.constants` and returns the indices worth keeping,
    plus notes carrying the measured distributions so a caller can print them.

    Rejection criteria, each with its reason:

    * **Bond length** outside ``[CA_CA_BOND_LENGTH - bond_tolerance, max_bond_length]``
      anywhere. A chain with a 6 A virtual bond is not a chain; no all-atom reconstruction
      can close it. The range is *asymmetric* on purpose -- see
      :data:`BOND_SCREEN_MAX_LENGTH`.
    * **Pseudo-angle** outside ``angle_window`` at *any* vertex -- see
      :data:`SCREEN_ANGLE_WINDOW` for why this is not a percentage. Pass ``None`` to skip
      the angle criterion entirely, which only a diagnostic caller isolating one of the
      other two should do.
    * **Internal clash**: any non-bonded CA-CA pair closer than ``internal_clash_cutoff``,
      which defaults to :data:`dodo.constants.CA_CLASH_DISTANCE`. It used to default to the
      *last rung* of :data:`dodo.constants.CLASH_RELAXATION_LADDER` (2.0 A), which let
      conformers through with alpha carbons 2.1 A apart -- closer than two bonded carbons.
      The ladder is for chains DODO is building and has nowhere else to go; a model
      ensemble has a hundred conformers, so an unusable one is dropped and the next is
      used. A rigid-body placement cannot fix an internal clash.

    Returns
    -------
    tuple[np.ndarray, tuple[str, ...]]
        Indices of kept conformers (ascending, possibly empty) and notes. An empty result
        is returned rather than raised: the caller knows how many it asked for and what to
        say about getting none.
    """
    array = _as_conformers(conformers)
    keep: list[int] = []
    notes: list[str] = []
    bond_extremes: list[tuple[float, float]] = []
    angle_extremes: list[tuple[float, float]] = []
    rejected: dict[str, int] = {"bond_length": 0, "pseudo_angle": 0, "internal_clash": 0}
    worst_contact = float("inf")
    sharpest = float("inf")
    longest = 0.0

    for index, chain in enumerate(array):
        if chain.shape[0] >= 2:
            bonds = ca_bond_lengths(chain)
            bond_extremes.append((float(bonds.min()), float(bonds.max())))
            if float(bonds.max()) > max_bond_length:
                longest = max(longest, float(bonds.max()))
                rejected["bond_length"] += 1
                continue
            if float(bonds.min()) < CA_CA_BOND_LENGTH - bond_tolerance:
                rejected["bond_length"] += 1
                continue
        if chain.shape[0] >= 3 and angle_window is not None:
            angles = ca_pseudo_angles(chain)
            finite = angles[np.isfinite(angles)]
            if finite.size:
                angle_extremes.append((float(finite.min()), float(finite.max())))
            # A NaN angle -- coincident or duplicated CAs -- counts as outside the window.
            # Comparing it away would let the most broken vertex there is pass unnoticed.
            outside = np.count_nonzero(~((angles >= angle_window[0]) & (angles <= angle_window[1])))
            if outside:
                if finite.size:
                    sharpest = min(sharpest, float(finite.min()))
                rejected["pseudo_angle"] += 1
                continue
        contact = _min_internal_ca_distance(chain)
        if contact < internal_clash_cutoff:
            worst_contact = min(worst_contact, contact)
            rejected["internal_clash"] += 1
            continue
        keep.append(index)

    if longest > max_bond_length:
        notes.append(
            f"At least one conformer had a CA-CA virtual bond of {longest:.3f} A, longer than "
            f"the {max_bond_length:.2f} A ceiling a trans peptide can make; no dihedral "
            f"produces such a bond, so the chain cannot be reconstructed to all-atom."
        )
    if angle_window is not None and np.isfinite(sharpest):
        notes.append(
            f"At least one conformer had a CA-CA-CA pseudo-angle of {sharpest:.1f} deg, "
            f"outside the observed {angle_window[0]:.0f}-{angle_window[1]:.0f} deg range; "
            f"such a vertex has no all-atom reconstruction."
        )
    if np.isfinite(worst_contact):
        notes.append(
            f"At least one conformer had non-bonded CA atoms {worst_contact:.3f} A apart, "
            f"inside the {internal_clash_cutoff:.2f} A clash distance; no rigid-body "
            f"placement can undo a conformer's own internal geometry."
        )
    if bond_extremes:
        low = min(pair[0] for pair in bond_extremes)
        high = max(pair[1] for pair in bond_extremes)
        notes.append(f"Ensemble CA-CA bonds spanned {low:.3f}-{high:.3f} A.")
    if angle_extremes:
        low = min(pair[0] for pair in angle_extremes)
        high = max(pair[1] for pair in angle_extremes)
        notes.append(
            f"Ensemble CA-CA-CA pseudo-angles spanned {low:.1f}-{high:.1f} deg "
            f"(DODO's generation window is {BACKBONE_ANGLE_MIN:.0f}-{BACKBONE_ANGLE_MAX:.0f})."
        )
    for reason, count in rejected.items():
        if count:
            notes.append(f"Rejected {count} of {array.shape[0]} conformer(s) on {reason}.")
    return np.asarray(keep, dtype=np.int64), tuple(notes)


# ---------------------------------------------------------------------------
# The engine
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class StarlingReport:
    """Everything about one :meth:`StarlingEngine.generate` call that the result cannot carry.

    :class:`~dodo.engines.base.IDRResult` is a fixed five-field contract shared by every
    engine, so the STARLING-specific measurements -- which cap was used, how big the
    ensemble was, what the anchor residuals came to -- live here and are returned by
    :meth:`StarlingEngine.generate_detailed`. Discarding them would leave a caller unable
    to answer "how far off were the anchors?", which is exactly the question the
    pre-rewrite code made unanswerable.
    """

    sequence_length: int
    cap: LengthCap
    ensemble_size: int
    kept: int
    placements: tuple[AnchorPlacement, ...]
    notes: tuple[str, ...]

    @property
    def worst_anchor_residual(self) -> float:
        """Largest anchor mismatch over the conformers returned, in Angstroms."""
        return max((p.anchor_residual for p in self.placements), default=0.0)

    def summary(self) -> str:
        """One-line summary for logs and figure captions."""
        placed = sum(1 for p in self.placements if p.ok)
        return (
            f"STARLING: {self.sequence_length} residues, cap {self.cap}, "
            f"{self.kept}/{self.ensemble_size} conformers usable, {placed}/"
            f"{len(self.placements)} placed, worst anchor residual "
            f"{self.worst_anchor_residual:.3f} A"
        )


class StarlingEngine:
    """Generate IDR conformers with idptools STARLING and place them between anchors.

    Satisfies the :class:`~dodo.engines.base.ConformationEngine` protocol.

    Parameters
    ----------
    ensemble_size
        Conformers to ask STARLING for. See :data:`DEFAULT_ENSEMBLE_SIZE`.
    device
        Passed through to STARLING when its signature accepts it (``"cpu"``, ``"cuda"``,
        ``"mps"``). ``None`` lets STARLING choose.
    oversample
        Minimum conformers generated per requested conformation.
    orientation_attempts
        Rigid-body orientations tried per conformer while avoiding obstacles.
    screen
        Run :func:`screen_conformers` before use. Leave on; the switch exists so a
        diagnostic caller can see what was rejected.
    end_to_end_tolerance
        Absolute tolerance in Angstroms on the achieved end-to-end distance. ``None`` uses
        :func:`end_to_end_tolerance`, i.e. the larger of
        :data:`END_TO_END_ABSOLUTE_TOLERANCE` and :data:`END_TO_END_RELATIVE_TOLERANCE`
        times the target.
    """

    name: str = "starling"

    def __init__(
        self,
        *,
        ensemble_size: int = DEFAULT_ENSEMBLE_SIZE,
        device: str | None = None,
        oversample: int = MIN_OVERSAMPLE,
        orientation_attempts: int = ORIENTATION_ATTEMPTS,
        regularize: bool = True,
        screen: bool = True,
        end_to_end_tolerance: float | None = None,
    ) -> None:
        if ensemble_size < 1:
            raise EngineError(f"ensemble_size must be at least 1, got {ensemble_size}.")
        if oversample < 1:
            raise EngineError(f"oversample must be at least 1, got {oversample}.")
        self.ensemble_size = int(ensemble_size)
        self.device = device
        self.oversample = int(oversample)
        self.orientation_attempts = int(orientation_attempts)
        self.regularize = bool(regularize)
        self.screen = bool(screen)
        # Validate now rather than at the first generate() call: a NaN tolerance would
        # disable the dimension check silently, which is the defect this field fixes.
        self.end_to_end_tolerance = (
            None if end_to_end_tolerance is None else _validated_tolerance(end_to_end_tolerance)
        )
        # Generated ensembles, keyed by (sequence, count). Per INSTANCE rather than global, so a
        # caller holding two engines with different settings cannot be served the other's output,
        # and the memory is released with the engine. See generate_detailed for why reuse is sound.
        self._ensembles: dict[tuple[str, int, int], tuple[np.ndarray, tuple[str, ...]]] = {}

    def __repr__(self) -> str:
        return (
            f"StarlingEngine(ensemble_size={self.ensemble_size}, device={self.device!r}, "
            f"available={self.available()})"
        )

    # -- availability ------------------------------------------------------

    def available(self) -> bool:
        """Whether STARLING can run here. Cheap, and never raises.

        Costs one :func:`importlib.util.find_spec` call when STARLING has not been
        imported yet. It deliberately does **not** import STARLING to check for model
        weights: that import is the multi-second torch startup this whole module is
        arranged to avoid, and engine selection calls this on every run. When STARLING
        *has* already been imported by someone else, the deeper checks (entry point
        present, weights on disk) are free and are done.

        A True answer therefore means "installed, and nothing cheap says it is broken".
        :meth:`generate` does the remaining checks and raises
        :class:`~dodo.exceptions.EngineUnavailableError` with specifics if they fail.
        """
        module = _loaded_starling()
        if module is None:
            return starling_installed()
        try:
            _generate_callable(module)
        except EngineUnavailableError:
            return False
        return _weights_present(module) is not False

    def max_length(self) -> LengthCap:
        """STARLING's maximum sequence length, read from the install when possible."""
        return starling_max_length(probe=self.available())

    # -- generation --------------------------------------------------------

    def generate(
        self,
        request: IDRRequest,
        obstacles: np.ndarray | None,
        rng: np.random.Generator,
    ) -> IDRResult:
        """Generate ``request.n_conformations`` anchored conformers.

        See :meth:`generate_detailed`, which this wraps.
        """
        result, _ = self.generate_detailed(request, obstacles, rng)
        return result

    def usable_conformers(
        self,
        sequence: str,
        *,
        n_conformations: int = 1,
        rng: np.random.Generator,
    ) -> tuple[np.ndarray, list[str]]:
        """Repaired, physically valid conformers for one sequence, with NO dimension selection.

        The pool :meth:`generate_detailed` would rank and choose from, handed over unranked. Nothing
        here knows or cares what end-to-end distance the caller wants, and that is the point: it
        exists for the conformer-driven placement mode, where the folded domains are moved to suit
        the conformer instead of the conformer being chosen to suit the domains.

        Worth being explicit about why that is not merely a convenience. Selecting conformers by
        end-to-end distance discards most of the ensemble and flattens what is left. Measured on
        dnmt3a's regions, STARLING's end-to-end distances scatter with a standard deviation of 41%
        of the mean, while the selection tolerance is 5% of the target -- so roughly 10% of
        conformers survive, and they are a narrow slice at the mean. The width of that distribution
        is the main thing STARLING knows that a single predicted number does not, so a caller that
        can avoid selecting on it should.

        Raises
        ------
        EngineError
            The sequence is longer than STARLING's cap, or no conformer survived repair.
        """
        _require_generator(rng)
        module = _import_starling()
        entry_point = _generate_callable(module)
        cap = starling_max_length()
        if len(sequence) > cap.value:
            raise EngineError(
                f"STARLING cannot generate a {len(sequence)}-residue sequence: its cap is {cap}."
            )
        # A dummy target, never read: _usable_conformers works from the sequence and the counts.
        # IDRRequest validates that the field is positive and finite, so it cannot simply be None.
        request = IDRRequest(
            sequence=sequence,
            n_residues=len(sequence),
            target_end_to_end=1.0,
            n_anchor_xyz=None,
            c_anchor_xyz=None,
            n_anchor_prev_xyz=None,
            c_anchor_next_xyz=None,
            n_conformations=n_conformations,
        )
        wanted = max(self.ensemble_size, n_conformations * self.oversample)
        usable, notes, _ = self._usable_conformers(request, entry_point, wanted, rng)
        if usable.shape[0] == 0:
            raise EngineError(
                f"No usable conformer for this {len(sequence)}-residue sequence. {'; '.join(notes)}"
            )
        return usable, notes

    def generate_detailed(
        self,
        request: IDRRequest,
        obstacles: np.ndarray | None,
        rng: np.random.Generator,
    ) -> tuple[IDRResult, StarlingReport]:
        """Generate conformers and return the measurements alongside them.

        Raises
        ------
        MissingDependencyError
            STARLING is not installed. Names the ``starling`` extra.
        EngineUnavailableError
            STARLING is installed but cannot be used: no recognizable entry point, or
            weights missing from the path it advertises.
        EngineError
            The request is longer than STARLING's cap (use
            :class:`~dodo.engines.hierarchical.HierarchicalEngine`), the request is
            internally inconsistent, or STARLING returned something this adapter cannot
            interpret as coordinates. Every one of these says which.
        """
        # The request validates itself: IDRRequest.__post_init__ rejects an empty region,
        # a sequence that disagrees with n_residues, a non-positive target and a non-finite
        # anchor. Nothing is re-checked here, because two copies of a validation rule drift.
        _require_generator(rng)

        module = _import_starling()
        entry_point = _generate_callable(module)
        if _weights_present(module) is False:
            raise EngineUnavailableError(
                "STARLING is installed but its model weights are not on disk at the path "
                "it advertises. Run STARLING once directly to download them, or point it "
                "at an existing weights file."
            )

        cap = starling_max_length()
        if request.n_residues > cap.value:
            raise EngineError(
                f"STARLING cannot generate a {request.n_residues}-residue region: its cap "
                f"is {cap}. Wrap this engine in "
                f"dodo.engines.hierarchical.HierarchicalEngine, which splits the region "
                f"into cap-sized segments and assembles them to the full-length target."
            )

        wanted = max(self.ensemble_size, request.n_conformations * self.oversample)
        usable, notes, generated = self._usable_conformers(request, entry_point, wanted, rng)
        if usable.shape[0] == 0:
            # Deduplicated, in order: every round reports the same drop reason, so concatenating
            # them verbatim produced a wall of text that repeated itself three times.
            reasons = list(dict.fromkeys(n for n in notes if n.startswith("Dropped")))
            raise EngineError(
                f"STARLING returned {generated} conformer(s) for this {request.n_residues}-residue "
                f"sequence over {MAX_GENERATION_ROUNDS} round(s) and none of them survived "
                f"geometric screening. {' '.join(reasons)}"
            )

        # Leave the aim point at least a tolerance of room inside the anchor-feasible
        # window, so that a conformer selected as "close enough" is still placeable. Clipped
        # onto the bare window edge, a conformer a fraction of an Angstrom on the wrong side
        # makes both anchors unreachable and the whole conformation is discarded.
        tolerance = end_to_end_tolerance(request.target_end_to_end, self.end_to_end_tolerance)
        desired, window = desired_internal_span(
            request.target_end_to_end,
            request.n_anchor_xyz,
            request.c_anchor_xyz,
            slack=tolerance,
        )
        if window is not None and not (window[0] <= request.target_end_to_end <= window[1]):
            notes.append(
                f"The requested end-to-end target ({request.target_end_to_end:.1f} A) lies "
                f"outside what the anchors permit ([{window[0]:.1f}, {window[1]:.1f}] A); "
                f"selection aimed at {desired:.1f} A instead."
            )
        order = rank_conformers(usable, desired, window=window)
        if request.n_conformations > order.size:
            notes.append(
                f"{request.n_conformations} conformations were requested but only "
                f"{order.size} usable conformer(s) came back, so some are reused with a "
                f"different rigid-body placement. They are not independent samples."
            )

        placements: list[AnchorPlacement] = []
        coords = np.empty((request.n_conformations, request.n_residues, 3), dtype=np.float64)
        success = np.zeros(request.n_conformations, dtype=bool)
        relaxations: list[float] = []
        for slot in range(request.n_conformations):
            # Distinct conformers per requested conformation where the ensemble allows it:
            # returning the same chain n times would be a fake ensemble.
            choice = int(order[slot % order.size])
            placement = place_between_anchors(
                usable[choice],
                n_anchor_xyz=request.n_anchor_xyz,
                c_anchor_xyz=request.c_anchor_xyz,
                rng=rng,
                obstacles=obstacles,
                conformer_index=int(choice),
                desired_end_to_end=desired,
                orientation_attempts=self.orientation_attempts,
            )
            placements.append(placement)
            coords[slot] = placement.ca_coords
            success[slot] = placement.ok
            if placement.relaxed_to is not None:
                relaxations.append(placement.relaxed_to)
            if not placement.ok:
                notes.append(
                    f"Conformation {slot} could not be placed (anchor residual "
                    f"{placement.anchor_residual:.3f} A, closest internal CA-CA contact "
                    f"{placement.min_internal_ca_distance:.3f} A)"
                    + (f": {'; '.join(placement.notes)}" if placement.notes else ".")
                )

        # Two independent things have to hold, and only the first used to be checked. A
        # placement can satisfy both anchors and clear every obstacle while having entirely
        # the wrong dimension: nothing in the rigid-body construction looks at the target on
        # a free or terminal region, so a 350 A request was answered with a 104 A chain and
        # reported as a success -- silently discarding the target and every
        # dodo.constants.MODES multiplier derived from it. The dimension is the scientific
        # content of the request, so it is a gate, not a note.
        dimension, dimension_note = _dimension_mask(
            np.array([p.achieved_end_to_end for p in placements], dtype=np.float64),
            desired=desired,
            tolerance=tolerance,
            window=window,
            n_residues=request.n_residues,
        )
        if dimension_note:
            notes.append(dimension_note)
        success &= dimension

        if not bool(success.any()):
            # base.py: partial success goes through the mask, total failure raises. Returning
            # an all-NaN IDRResult(0/N built) instead means a caller written against that
            # contract treats the object as a result and reads NaN coordinates out of it.
            raise ExhaustedAttemptsError(
                f"None of the {request.n_conformations} requested conformation(s) of this "
                f"{request.n_residues}-residue region could be built from the "
                f"{usable.shape[0]} usable conformer(s) STARLING returned. "
                f"{' '.join(notes)}",
                attempts=generated,
            )

        # from_batch, not the constructor: IDRResult requires the rows of failed conformers
        # to be non-finite so that nothing downstream can mistake them for a structure, and
        # from_batch is what enforces that. The real coordinates of a failed placement are
        # not lost -- they stay on its AnchorPlacement, where a caller can inspect why it
        # failed without any risk of writing them to a file.
        result = IDRResult.from_batch(
            ca_coords=coords,
            success=success,
            engine=self.name,
            attempts=generated,
            # min, not max: with two rungs in play the threshold the chain as a whole was
            # actually built at is the loosest one. max understates the relaxation, which
            # is the opposite of what relaxed_to is for.
            relaxed_to=min(relaxations) if relaxations else None,
        )
        report = StarlingReport(
            sequence_length=request.n_residues,
            cap=cap,
            ensemble_size=generated,
            kept=int(usable.shape[0]),
            placements=tuple(placements),
            notes=tuple(notes),
        )
        return result, report

    # -- the actual call ---------------------------------------------------

    def _usable_conformers(
        self,
        request: IDRRequest,
        entry_point: Callable[..., Any],
        wanted: int,
        rng: np.random.Generator,
    ) -> tuple[np.ndarray, list[str], int]:
        """Obtain enough physically valid conformers, repairing only as many as needed.

        The order of operations is the point, and it is chosen so the expensive step runs least:

        1. **Generate** a round of ``wanted`` conformers, cached per sequence and round.
        2. **Eliminate** what is not worth repairing -- a lost chain, or one that already collides
           with itself -- using only checks that cost nothing. See :func:`prefilter_conformers`.
        3. **Repair** every survivor, best first. Note *every*: an earlier version of this stopped
           as soon as it had the requested number, which looked like the obvious saving and was
           wrong. Selection downstream needs **spread** in end-to-end distance, not merely a
           quantity of conformers -- an interior IDR has to span whatever its flanking domains
           impose. Truncating the pool at 2 left both survivors the wrong size and the region
           failed on dimension with a 12.7 A miss, having repaired only 2 of 16 candidates. The
           saving comes from step 2 instead, which is free and does not touch the pool's spread.
        4. **Reject** a repair that moved an atom further than :data:`MAX_REPAIR_DISPLACEMENT`, or
           that could not satisfy both constraints, or that still fails the physical screen.
        5. **Generate another round** if that did not yield enough, up to
           :data:`MAX_GENERATION_ROUNDS`.

        Yield is the reason this is a loop rather than a single pass: at 200 residues only about 6
        conformers in 30 survive, because STARLING's raw output already contains internal clashes
        that repair cannot mend. A single round can therefore come up short through no fault of
        DODO's, and the honest response is to ask for more rather than fail.
        """
        needed = max(request.n_conformations * self.oversample, request.n_conformations)
        accepted: list[np.ndarray] = []
        notes: list[str] = []
        generated = 0

        for round_index in range(MAX_GENERATION_ROUNDS):
            raw, round_notes = self._ensemble_for(request, entry_point, wanted, round_index, rng)
            notes.extend(round_notes)
            generated += int(raw.shape[0])

            if self.screen:
                order, prefilter_notes = prefilter_conformers(raw)
                notes.extend(prefilter_notes)
            else:
                order = np.arange(raw.shape[0], dtype=np.int64)
                notes.append("Screening disabled by caller; conformer geometry is unchecked.")

            repaired_count = 0
            budget_rejects = 0
            screen_rejects = 0
            for index in order:
                chain = raw[int(index)]
                if not self.regularize:
                    accepted.append(chain)
                    continue
                repaired_count += 1
                fixed, ok = _repair_one(chain)
                if not ok:
                    budget_rejects += 1
                    continue
                if self.screen:
                    survived, _ = screen_conformers(fixed[None, :, :])
                    if survived.size == 0:
                        screen_rejects += 1
                        continue
                accepted.append(fixed)

            notes.append(
                f"Round {round_index + 1}: repaired {repaired_count} of {order.size} candidate(s) "
                f"to keep {len(accepted)} (wanted at least {needed}); {budget_rejects} exceeded "
                f"the repair budget, "
                f"{screen_rejects} still failed the physical screen."
            )
            if len(accepted) >= needed:
                break

        if not accepted:
            return np.empty((0, request.n_residues, 3)), notes, generated
        return np.stack(accepted), notes, generated

    def _ensemble_for(
        self,
        request: IDRRequest,
        entry_point: Callable[..., Any],
        wanted: int,
        round_index: int,
        rng: np.random.Generator,
    ) -> tuple[np.ndarray, list[str]]:
        """Raw conformers for one round, generated once per (sequence, size, round).

        Cached because STARLING conditions on sequence alone: the ensemble it returns for a region
        is the same every time it is asked, so what differs between DODO's models is which
        conformer is SELECTED and where it is PLACED, both of which happen downstream of here.
        Without this, a 2-model rebuild of a 3-region structure ran 6 diffusion-plus-MDS
        generations to obtain 3 distinct ensembles.

        ``round_index`` is part of the key so that asking for a second round yields *new*
        conformers rather than the first round's again.
        """
        # Drawn HERE, before the cache is consulted, so the generator advances by the same amount
        # whether or not this ensemble already existed. With the draw on the miss path only, a
        # second request for the same sequence left the rng in a different state and selection and
        # placement diverged -- an existing reproducibility test caught it, two calls with one seed
        # returning conformers 41.9 A apart.
        starling_seed = int(rng.integers(0, 2**31 - 1))
        key = (request.sequence, wanted, round_index)
        cached = self._ensembles.get(key)
        if cached is not None:
            return cached[0].copy(), [
                f"Reused the round-{round_index + 1} ensemble already generated for this "
                f"{request.n_residues}-residue sequence."
            ]
        raw, notes = self._call_starling(entry_point, request.sequence, wanted, starling_seed)
        conformers, unit_notes = _to_angstroms(raw)
        notes.extend(unit_notes)
        if round_index:
            notes.append(
                f"Generated an additional round of {wanted} conformer(s): earlier rounds did not "
                f"yield enough usable ones."
            )
        self._ensembles[key] = (conformers.copy(), tuple(notes))
        return conformers, notes

    def _call_starling(
        self,
        entry_point: Callable[..., Any],
        sequence: str,
        n_conformers: int,
        seed: int,
    ) -> tuple[np.ndarray, list[str]]:
        """Call STARLING's entry point, passing only keywords its signature accepts.

        STARLING's signature has changed across releases, so passing a fixed keyword set
        would break on any version that renamed one. Introspecting instead means an
        unknown keyword is dropped with a note rather than raising a TypeError from inside
        someone else's package.

        Reproducibility: DODO's rule is that every stochastic call is seeded. STARLING
        samples with torch, so the only lever is its own seed argument. When the signature
        exposes one it gets an integer drawn from ``rng``; when it does not, that is
        recorded as a note, because a caller relying on seed-identical output deserves to
        be told this engine cannot promise it.
        """
        notes: list[str] = []
        wanted: dict[str, Any] = {
            "conformations": n_conformers,
            # THE keyword this adapter cannot do without. STARLING's generate() defaults to
            # return_structures=False, which returns an Ensemble carrying only DISTANCE MAPS --
            # no 3D coordinates at all. Without this every request failed at coordinate
            # extraction with "could not find coordinates on the object STARLING returned
            # (Ensemble)", which reads like a version mismatch and is really a missing argument.
            "return_structures": True,
            "return_data": True,
            "verbose": False,
            "show_progress_bar": False,
            # STARLING draws a second, per-step bar. DODO has its own progress reporting.
            "show_per_step_progress_bar": False,
            "seed": seed,
        }
        if self.device is not None:
            wanted["device"] = self.device
        accepted, takes_kwargs = _acceptable_keywords(entry_point)
        if accepted is None or takes_kwargs:
            passed = wanted
        else:
            passed = {key: value for key, value in wanted.items() if key in accepted}
            dropped = sorted(set(wanted) - set(passed))
            if dropped:
                notes.append(
                    f"The installed STARLING's generate() does not accept {dropped}; "
                    f"called without them."
                )
        if "conformations" not in passed and not takes_kwargs:
            raise EngineUnavailableError(
                "The installed STARLING's generate() has no 'conformations' parameter, so "
                "DODO cannot ask it for a specific number of conformers. This adapter "
                "targets STARLING 2.x; check the installed version."
            )
        if "seed" not in passed:
            notes.append(
                "The installed STARLING's generate() takes no seed, so its sampling is not "
                "reproducible from DODO's rng. Conformer *selection* and placement still "
                "are."
            )

        try:
            raw = entry_point(sequence, **passed)
        except Exception as exc:
            # Anything STARLING raises -- a missing weight file, a CUDA error, a shape
            # mismatch -- must surface as a DODO error naming the engine, not as a bare
            # third-party traceback the user cannot act on. The original is chained.
            raise EngineError(
                f"STARLING failed while generating {n_conformers} conformer(s) for a "
                f"{len(sequence)}-residue sequence: {type(exc).__name__}: {exc}"
            ) from exc

        coords = _coordinates_from_result(raw, sequence, len(sequence))
        return coords, notes


# ---------------------------------------------------------------------------
# Interpreting STARLING's return value
# ---------------------------------------------------------------------------


def _acceptable_keywords(function: Callable[..., Any]) -> tuple[set[str] | None, bool]:
    """Keyword names ``function`` accepts, and whether it has a ``**kwargs`` catch-all.

    ``(None, False)`` means the signature could not be read at all, which happens for
    C-implemented callables; the caller then passes everything and lets it fail loudly.
    """
    try:
        signature = inspect.signature(function)
    except (TypeError, ValueError):  # pragma: no cover - C callables
        return None, False
    names: set[str] = set()
    takes_kwargs = False
    for parameter in signature.parameters.values():
        if parameter.kind is inspect.Parameter.VAR_KEYWORD:
            takes_kwargs = True
        elif parameter.kind in (
            inspect.Parameter.POSITIONAL_OR_KEYWORD,
            inspect.Parameter.KEYWORD_ONLY,
        ):
            names.add(parameter.name)
    return names, takes_kwargs


def _coordinates_from_result(result: Any, sequence: str, n_residues: int) -> np.ndarray:
    """Extract ``(n_conformers, n_residues, 3)`` coordinates from STARLING's return value.

    STARLING returns a mapping of input name to an ``Ensemble`` object. Both the mapping
    and the attribute holding coordinates are ASSUMPTIONS about its API (see
    :data:`COORDINATE_ATTRIBUTES`), so every step is checked and a failure says exactly
    what was found instead. Guessing here would mean handing the rest of DODO an array of
    the wrong thing -- a distance map, say -- which would still be numbers.
    """
    payload = result
    if isinstance(payload, dict):
        # MEASURED: STARLING keys this dict "sequence_1", "sequence_2", ... -- positionally,
        # NOT by the sequence string. A single-sequence request therefore only ever worked
        # through the len == 1 branch below.
        if sequence in payload:
            payload = payload[sequence]
        elif len(payload) == 1:
            payload = next(iter(payload.values()))
        else:
            raise EngineError(
                f"STARLING returned a mapping with {len(payload)} entries and none keyed "
                f"by the requested sequence, so DODO cannot tell which ensemble is the "
                f"answer. Keys: {sorted(map(str, payload))[:5]}"
            )

    array = _extract_array(payload)
    if array is None:
        raise EngineError(
            f"DODO could not find coordinates on the object STARLING returned "
            f"({type(payload).__name__}). It probed the attributes "
            f"{list(COORDINATE_ATTRIBUTES)} and a .trajectory.xyz fallback. This is an "
            f"API mismatch; check the installed STARLING version."
        )

    coords = np.asarray(array, dtype=np.float64)
    if coords.ndim == 2:
        coords = coords[None, :, :]
    if coords.ndim != 3 or coords.shape[2] != 3:
        raise EngineError(
            f"STARLING coordinates have shape {coords.shape}; DODO expects "
            f"(n_conformers, n_residues, 3). A 2D array of shape (n, n) would be a "
            f"distance map, which is not what this adapter asked for."
        )
    if coords.shape[1] != n_residues:
        raise EngineError(
            f"STARLING returned {coords.shape[1]} residues per conformer for a "
            f"{n_residues}-residue sequence. Refusing to guess whether the difference is "
            f"a terminal cap, a truncation or a transposed array."
        )
    if coords.shape[0] == 0:
        raise EngineError("STARLING returned an ensemble with zero conformers.")
    finite = np.isfinite(coords).all(axis=(1, 2))
    if not finite.all():
        bad = int(np.count_nonzero(~finite))
        if not finite.any():
            raise EngineError(
                f"All {coords.shape[0]} STARLING conformer(s) contain non-finite coordinates, "
                f"so there is nothing to build from. This is a STARLING generation failure, not "
                f"a DODO one -- check the installed version and the sequence."
            )
        # Dropped, NOT fatal. This used to raise whenever ANY conformer was non-finite, which is
        # the wrong inference from a correct premise: a NaN must never reach a structure file, but
        # the way to ensure that is to discard the conformers carrying it, not to throw away the
        # good ones with them. MEASURED on a p300 rebuild, where it made every region fail: the
        # bad fraction ran 2-21% per region, so 79-98% of each ensemble was usable and discarded.
        # A region failing here keeps its INPUT coordinates while the folded domains are still
        # repositioned, so the visible symptom was a structure whose IDRs connected to domains
        # that had moved out from under them.
        warnings.warn(
            f"Dropped {bad} of {coords.shape[0]} STARLING conformers for non-finite "
            f"coordinates; building from the remaining {int(finite.sum())}. A NaN here comes "
            f"from STARLING's generation, not from DODO.",
            UserWarning,
            stacklevel=2,
        )
        coords = coords[finite]
    return coords


def _extract_array(payload: Any) -> Any:
    """Probe ``payload`` for something array-like. Returns None if nothing matches."""
    if isinstance(payload, np.ndarray):
        return payload
    for attribute in COORDINATE_ATTRIBUTES:
        value = getattr(payload, attribute, None)
        if value is None:
            continue
        if callable(value):
            try:
                value = value()
            except TypeError:
                # A method needing arguments is not the accessor we are looking for.
                continue
        if value is not None and np.ndim(value) >= 2:
            return value
    # MEASURED against STARLING 2.x: `Ensemble.trajectory` is a soursop SSProtein, which has
    # no `.xyz` of its own -- the array lives one level further down on the mdtraj object it
    # wraps, at `.trajectory.traj.xyz`, shaped (n_conformers, n_residues, 3) and in NANOMETRES.
    # The bare `.trajectory.xyz` probe below it is kept for older layouts that did expose it.
    trajectory = getattr(payload, "trajectory", None)
    for path in ("traj.xyz", "xyz"):
        probed: Any = trajectory
        for part in path.split("."):
            probed = getattr(probed, part, None)
            if probed is None:
                break
        if probed is not None and np.ndim(probed) >= 2:
            return probed
    return None


def _to_angstroms(coords: np.ndarray) -> tuple[np.ndarray, list[str]]:
    """Return coordinates in Angstroms, converting from nanometres if that is what they are.

    Why this exists: the mdtraj trajectory STARLING hands out is in NANOMETRES, while some routes
    give Angstroms. A silent factor of ten would produce a structure with 0.38 A bonds, which every
    viewer draws -- wrongly -- as a solid blob, so the unit is measured rather than assumed.

    The unit and the *quality* of the geometry are deliberately separate questions here, and
    conflating them is how two successive versions of this function got it wrong. The first
    required the median bond to sit within 0.5 A of 3.81 or within 0.05 of 0.381; the second used
    the same 0.5 A on the scaled-up value. Both are too strict, because STARLING's reconstruction
    error is large and grows with length: a 38-residue sequence measured a median of 3.36 A and a
    401-residue one 3.17 A, so the second version rejected a perfectly ordinary 401-mer as having
    units that were "neither Angstroms nor nanometres".

    Units differ by a factor of ten, so that decision is unambiguous by ORDER OF MAGNITUDE and
    needs no tolerance at all. Whether the geometry is usable is then a separate test, made
    against :data:`RECONSTRUCTION_BOND_WINDOW` -- wide enough for real reconstruction error, and
    still narrow enough that a genuinely garbled array fails it.
    """
    notes: list[str] = []
    if coords.shape[1] < 2:
        return coords, notes
    median_bond = float(np.median(np.linalg.norm(np.diff(coords, axis=1), axis=2)))

    # Nanometres or Angstroms, whichever leaves the median nearer a real CA-CA bond. The two
    # candidates are a factor of ten apart, so this is never a close call.
    if abs(median_bond * 10.0 - CA_CA_BOND_LENGTH) < abs(median_bond - CA_CA_BOND_LENGTH):
        coords = coords * 10.0
        median_bond *= 10.0
        notes.append(
            f"Coordinates were in nanometres (median virtual bond {median_bond / 10.0:.4f}) and "
            f"have been multiplied by 10."
        )

    low, high = RECONSTRUCTION_MEDIAN_WINDOW
    if not low <= median_bond <= high:
        raise EngineError(
            f"STARLING coordinates have a median CA-CA distance of {median_bond:.3f} A, outside "
            f"the {low}-{high} A window any reconstruction of a real chain falls in, so neither "
            f"unit interpretation makes this a protein backbone. DODO will not write it to a file."
        )
    return coords, notes


# ---------------------------------------------------------------------------
# Small shared numerics
# ---------------------------------------------------------------------------

#: Below this length a vector is treated as having no direction, in Angstroms. Matches the
#: tolerance in :mod:`dodo.geometry.transforms`, which raises rather than returning NaN.
_DEGENERATE_LENGTH = 1e-9


def _require_generator(rng: np.random.Generator) -> np.random.Generator:
    """Reject a seed integer or a legacy RandomState with an actionable message."""
    if not isinstance(rng, np.random.Generator):
        raise TypeError(
            f"rng must be a numpy.random.Generator, got {type(rng).__name__}. Construct "
            f"one with numpy.random.default_rng(seed); every stochastic path in DODO takes "
            f"an explicit generator so that a build can be reproduced exactly."
        )
    return rng


def _as_conformers(conformers: np.ndarray) -> np.ndarray:
    """Coerce input to ``(n_conformers, n_residues, 3)``, raising on anything else."""
    array = np.asarray(conformers, dtype=np.float64)
    if array.ndim == 2:
        array = array[None, :, :]
    if array.ndim != 3 or array.shape[2] != 3 or array.shape[0] == 0 or array.shape[1] == 0:
        raise GeometryError(
            f"Expected conformers of shape (n_conformers, n_residues, 3) with both counts "
            f"non-zero, got {np.shape(conformers)}."
        )
    if not np.all(np.isfinite(array)):
        raise GeometryError(
            "Conformer coordinates contain non-finite values. A NaN here came from a "
            "failed sampler upstream, which should have raised instead of returning it."
        )
    return array


def _as_single_conformer(conformer: np.ndarray) -> np.ndarray:
    """Coerce input to ``(n_residues, 3)``, raising on anything else."""
    array = np.asarray(conformer, dtype=np.float64)
    if array.ndim != 2 or array.shape[1] != 3 or array.shape[0] == 0:
        raise GeometryError(
            f"Expected a conformer of shape (n_residues, 3), got {np.shape(conformer)}."
        )
    if not np.all(np.isfinite(array)):
        raise GeometryError("Conformer coordinates contain non-finite values.")
    return array


def _perpendicular(unit_axis: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Draw a random unit vector perpendicular to ``unit_axis``.

    Drawn by projecting a uniform direction off the axis, so the result is uniform on the
    circle. Retried on the measure-zero event that the draw was parallel to the axis --
    normalizing the near-zero residual there would amplify float error into a vector that
    is not perpendicular at all.
    """
    for _ in range(8):
        candidate = random_unit_vectors(1, rng)[0]
        residual = candidate - np.dot(candidate, unit_axis) * unit_axis
        norm = float(np.linalg.norm(residual))
        if norm > 1e-6:
            unit: np.ndarray = residual / norm
            return unit
    # Deterministic fallback: cross with whichever basis vector is least aligned.
    basis = np.eye(3)[int(np.argmin(np.abs(unit_axis)))]
    residual = np.cross(unit_axis, basis)
    fallback: np.ndarray = residual / float(np.linalg.norm(residual))
    return fallback


def _clash_ladder(cutoff: float, ladder: Sequence[float]) -> tuple[float, ...]:
    """Return the clash thresholds to try, strictest first, starting from ``cutoff``."""
    rungs = [float(cutoff)] + [float(rung) for rung in ladder if float(rung) < float(cutoff)]
    return tuple(sorted(set(rungs), reverse=True))


def _clash_free(coords: np.ndarray, obstacles: np.ndarray | None, cutoff: float) -> bool:
    """Whether no point of ``coords`` comes within ``cutoff`` of an obstacle."""
    if obstacles is None:
        return True
    obstacle_array = np.asarray(obstacles, dtype=np.float64)
    if obstacle_array.ndim != 2 or obstacle_array.shape[1] != 3:
        raise GeometryError(
            f"obstacles must have shape (n_obstacles, 3), got {obstacle_array.shape}."
        )
    if obstacle_array.shape[0] == 0:
        return True
    distances, _ = cKDTree(obstacle_array).query(coords, k=1)
    return bool(np.min(distances) >= cutoff)


def _min_internal_ca_distance(chain: np.ndarray) -> float:
    """Closest approach between CA pairs separated by more than the bonded exclusion.

    Pairs within :data:`dodo.constants.CLASH_EXCLUDE_WITHIN_RESIDUES` in sequence are
    excluded: ``i, i+1`` are covalently bonded at 3.8 A and ``i, i+2`` are held at
    5.0-7.5 A by the backbone angle, so counting either as a clash reports every peptide
    bond as a collision -- which is precisely what the pre-rewrite whole-structure check
    did.

    Returns ``inf`` for a chain too short to have any non-bonded pair, which is the honest
    answer: there is no closest non-bonded approach.
    """
    n = chain.shape[0]
    if n <= CLASH_EXCLUDE_WITHIN_RESIDUES + 1:
        return float("inf")
    best = float("inf")
    # One vectorized pass per sequence-separation band, rather than an n^2 distance
    # matrix: bands keep peak memory at O(n) for a chain of any length DODO builds.
    for offset in range(CLASH_EXCLUDE_WITHIN_RESIDUES + 1, n):
        gaps = np.linalg.norm(chain[offset:] - chain[:-offset], axis=1)
        best = min(best, float(gaps.min()))
    return best
