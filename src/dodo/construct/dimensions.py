"""Target dimensions for a disordered region.

This module is the scientific premise of DODO. Everything else -- parsing, region
identification, coordinate generation -- exists to serve one question: *how big should
this disordered region be?* An IDR rebuilt to the wrong dimensions is just differently
wrong from AlphaFold's extended spaghetti.

The answer comes from ALBATROSS (via sparrow), which predicts mean end-to-end distance
from sequence. When sparrow is not installed, an analytical polymer scaling law stands in;
it is sequence-blind and documented as an approximation, not a substitute.

What changed from v1, and why
-----------------------------
v1 expressed its named build modes as **Angstroms per residue** -- ``normal`` was 0.8
A/residue, ``max_expansion`` 1.65, and so on -- with ``predicted`` as a separate mode that
called ALBATROSS. That is linear in chain length, but real IDR end-to-end distance scales
as roughly :math:`N^{0.55}`, so a fixed multiplier can only agree with the prediction at
one length. ``normal`` gives 80 A at N=100, which is about right, and 400 A at N=500,
where the prediction is nearer 190 A. The error grows without bound.

Worse, the first v2 attempt dropped ALBATROSS entirely and substituted two *disagreeing*
placeholders -- ``1.0 * n_residues`` in the IDR builder and ``1.4 * n_residues`` for
folded-domain spacing -- so the walk was actively fighting its own anchors.

Here the modes are multipliers **on the predicted end-to-end distance**. ``expanded``
means 1.3x whatever this sequence's predicted dimension is, at any length. The vocabulary
users know is preserved; the numbers underneath are length-independent and physically
meaningful. One consequence is documented and deliberate: ``normal`` and ``predicted`` are
now synonyms, which they were not in v1.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from functools import lru_cache
from typing import Literal

from ..constants import (
    ALBATROSS_MIN_LENGTH,
    CA_CA_BOND_LENGTH,
    DEFAULT_MODE,
    contour_length,
    flory_end_to_end,
    resolve_mode,
)
from ..exceptions import MissingDependencyError, UnsatisfiableTargetError

__all__ = [
    "DimensionTarget",
    "albatross_available",
    "predict_end_to_end",
    "predict_radius_of_gyration",
    "predict_scaling_exponent",
    "target_dimensions",
]

#: Where a predicted dimension came from. Carried on every result so a caller -- or a
#: figure caption -- can state which, rather than guessing from install state.
PredictionSource = Literal["albatross", "analytical"]

#: Fraction of the contour length a target is clamped to when it would otherwise exceed
#: what the chain can physically span.
#:
#: Not 1.0: a chain stretched to its full contour length is a straight rod with exactly
#: one conformation, so the sampler would have no freedom and every clash would be fatal.
#: 0.95 leaves enough slack to find a conformation while still honouring an extreme
#: request as far as geometry allows.
_MAX_CONTOUR_FRACTION = 0.95


@dataclass(frozen=True, slots=True)
class DimensionTarget:
    """The dimension a region should be built to, with its full provenance.

    Every field exists so that a build is reproducible and reportable. The pre-rewrite
    code passed a bare float around, which meant that when output dimensions looked wrong
    there was no way to tell whether the prediction, the mode multiplier or the clamp was
    responsible.

    Attributes
    ----------
    end_to_end
        The target CA-to-CA end-to-end distance in Angstroms. This is what the builder
        aims for.
    n_residues
        Length of the region, in residues.
    mode
        The named build mode requested.
    factor
        The multiplier that ``mode`` resolved to.
    predicted
        The unscaled predicted end-to-end distance, before ``factor`` was applied.
    source
        Which predictor produced ``predicted``.
    clamped
        True if ``end_to_end`` had to be reduced to stay physically achievable. When this
        is set, ``end_to_end`` is *not* ``predicted * factor``.
    max_achievable
        The ceiling used when clamping, i.e. a fraction of the contour length.
    """

    end_to_end: float
    n_residues: int
    mode: str
    factor: float
    predicted: float
    source: PredictionSource
    clamped: bool
    max_achievable: float

    def __str__(self) -> str:
        detail = f"{self.end_to_end:.1f} A"
        if self.factor != 1.0:
            detail += f" ({self.mode}, {self.factor}x of {self.predicted:.1f} A)"
        else:
            detail += f" ({self.mode})"
        if self.clamped:
            detail += f" [clamped from {self.predicted * self.factor:.1f} A]"
        return f"{detail} via {self.source} over {self.n_residues} residues"

    @property
    def per_residue(self) -> float:
        """Target end-to-end distance divided by residue count, in Angstroms.

        Provided only for comparison against v1's A/residue modes when auditing a change
        in output. Do not build anything on it: it is length-dependent, which is the
        property this module exists to remove.
        """
        return self.end_to_end / self.n_residues


def albatross_available() -> bool:
    """Return True if sparrow is importable, so ALBATROSS predictions can be made.

    Feature-detected by import and attribute presence rather than by version comparison.
    sparrow's own README advertises a 1.0.x version while an installed copy reports
    0.2.2, and the project publishes no releases, so a version string is not a usable
    signal about what the API offers.
    """
    return _sparrow_predictor_factory() is not None


@lru_cache(maxsize=1)
def _sparrow_predictor_factory() -> type | None:
    """Import and cache sparrow's ``Protein`` class, or return None if unavailable.

    Imported lazily and exactly once. sparrow pulls in torch, so a module-scope import
    would make even ``dodo --help`` pay several seconds of framework startup -- which is
    the reason sparrow is an optional extra in the first place.
    """
    try:
        from sparrow import Protein
    except ImportError:
        return None
    # Guard against an installed-but-incompatible sparrow: the predictor attribute and the
    # specific method are what we actually depend on.
    probe = getattr(Protein, "predictor", None)
    if probe is None:
        return None
    factory: type = Protein
    return factory


def _require_sequence(sequence: str) -> str:
    """Validate and normalize a sequence for prediction."""
    cleaned = sequence.strip().upper()
    if not cleaned:
        raise ValueError("Cannot predict dimensions for an empty sequence.")
    if not cleaned.isalpha():
        offenders = sorted({c for c in cleaned if not c.isalpha()})
        raise ValueError(
            f"Sequence contains non-alphabetic characters {offenders}. Strip gaps, "
            f"numbering and whitespace before predicting dimensions."
        )
    return cleaned


def predict_end_to_end(
    sequence: str,
    *,
    prefer_albatross: bool = True,
    warn_on_fallback: bool = True,
) -> tuple[float, PredictionSource]:
    """Predict the mean end-to-end distance of a disordered sequence.

    Parameters
    ----------
    sequence
        One-letter amino acid sequence of the disordered region.
    prefer_albatross
        Use ALBATROSS when sparrow is installed. Set False to force the analytical
        estimate, which is useful for testing and for reproducing a result on a machine
        without sparrow.
    warn_on_fallback
        Emit a :class:`UserWarning` when sparrow is missing and the analytical estimate is
        used instead. On by default because the two are not interchangeable and a silent
        downgrade would make two runs of the "same" command disagree with no visible
        cause.

        Only an *unintentional* downgrade warns. Passing ``prefer_albatross=False`` is a
        deliberate choice and is never warned about, regardless of this setting.

    Returns
    -------
    tuple[float, PredictionSource]
        The predicted end-to-end distance in Angstroms, and which predictor produced it.

    Notes
    -----
    Short sequences need no special handling. sparrow exposes ``use_scaled`` on
    :meth:`end_to_end_distance` and documents ``MIN_LENGTH_ALBATROSS_RE_RG`` = 35 as the
    length below which the scaled network is required, but measurement shows sparrow
    already ignores the flag below that threshold and takes the correct path itself
    (at N=34 scaled and unscaled agree exactly; at N=35 they diverge). So DODO passes
    ``use_scaled=True``, which is also sparrow's own default, and does not second-guess it.
    This matters because most loops DODO rebuilds are shorter than 35 residues.
    """
    cleaned = _require_sequence(sequence)
    n = len(cleaned)

    if prefer_albatross:
        protein_cls = _sparrow_predictor_factory()
        if protein_cls is not None:
            predicted = float(protein_cls(cleaned).predictor.end_to_end_distance(use_scaled=True))
            if predicted <= 0.0:
                raise ValueError(
                    f"ALBATROSS returned a non-positive end-to-end distance "
                    f"({predicted}) for a {n}-residue sequence. This is a sparrow "
                    f"problem, not a DODO one; report it upstream."
                )
            return predicted, "albatross"

        if warn_on_fallback:
            warnings.warn(
                "sparrow is not installed, so DODO is estimating end-to-end distance "
                "from an analytical polymer scaling law instead of predicting it with "
                "ALBATROSS. The estimate is blind to sequence composition and errs "
                "compact for charged or proline-rich sequences. For sequence-specific "
                "predictions install sparrow directly -- it is not on PyPI, so there is no "
                "extra for it: pip install git+https://github.com/idptools/sparrow.git",
                UserWarning,
                stacklevel=2,
            )

    return flory_end_to_end(n), "analytical"


def predict_radius_of_gyration(sequence: str) -> tuple[float, PredictionSource]:
    """Predict the radius of gyration of a disordered sequence, in Angstroms.

    Used to validate built conformations rather than to drive them: DODO builds to an
    end-to-end target, and Rg is the independent check that the resulting coil is not
    merely the right length end to end while being the wrong shape overall.

    Falls back to ``Re / sqrt(6)`` (the ideal-chain relation) without sparrow, which is
    cruder still than the end-to-end fallback since it compounds two approximations.
    """
    cleaned = _require_sequence(sequence)
    protein_cls = _sparrow_predictor_factory()
    if protein_cls is not None:
        value = float(protein_cls(cleaned).predictor.radius_of_gyration(use_scaled=True))
        return value, "albatross"
    return flory_end_to_end(len(cleaned)) / 6.0**0.5, "analytical"


def predict_scaling_exponent(sequence: str) -> float:
    """Predict the apparent Flory scaling exponent for a sequence.

    Diagnostic rather than structural. A value near 0.5-0.6 indicates a genuinely
    disordered, solvent-expanded chain; a value well below that (measurement puts a
    hydrophobic-rich sequence at ~0.20) means the sequence is predicted to collapse into a
    globule and is not really an IDR -- in which case rebuilding it as an expanded coil is
    the wrong thing to do, and DODO should say so rather than do it.

    Raises
    ------
    MissingDependencyError
        There is no meaningful analytical stand-in: a length-only law has a *fixed*
        exponent by construction, so returning one would answer a different question than
        the caller asked.
    """
    cleaned = _require_sequence(sequence)
    protein_cls = _sparrow_predictor_factory()
    if protein_cls is None:
        raise MissingDependencyError(
            package="sparrow",
            purpose="Predicting a sequence-specific scaling exponent",
            extra="albatross",
        )
    return float(protein_cls(cleaned).predictor.scaling_exponent())


def target_dimensions(
    sequence: str,
    *,
    mode: str = DEFAULT_MODE,
    prefer_albatross: bool = True,
    warn_on_fallback: bool = True,
    clamp: bool = True,
    warn_on_clamp: bool = True,
) -> DimensionTarget:
    """Resolve a sequence and a build mode into a concrete dimension target.

    This is the function the builders call.

    Parameters
    ----------
    sequence
        One-letter sequence of the disordered region.
    mode
        A named build mode: ``super_compact``, ``compact``, ``normal``/``predicted``,
        ``expanded``, ``super_expanded`` or ``max_expansion``. Each is a multiplier on the
        predicted end-to-end distance.
    prefer_albatross
        Use ALBATROSS if available.
    warn_on_fallback
        Warn when the analytical estimate is used instead.
    clamp
        Reduce a target that exceeds what the chain can physically span. When False, an
        unachievable target raises instead, which is what a caller wanting to fail loudly
        should choose.
    warn_on_clamp
        Warn when a target is reduced. On by default: the predicted distance is what a user
        believes they are getting, so substituting a different one silently would change the
        science without saying so. Turn it off only when the caller surfaces the
        :attr:`DimensionTarget.clamped` flag itself.

    Returns
    -------
    DimensionTarget
        The target and its provenance.

    Raises
    ------
    ValueError
        If ``mode`` is not recognized, or the sequence is empty or non-alphabetic.
    UnsatisfiableTargetError
        If ``clamp`` is False and the requested target exceeds the contour length.

    Examples
    --------
    >>> target = target_dimensions("GSGSGSGSGSGSGSGSGSGS", mode="compact")  # doctest: +SKIP
    >>> target.end_to_end < target.predicted  # compact means smaller than predicted
    True
    """
    cleaned = _require_sequence(sequence)
    n = len(cleaned)
    factor = resolve_mode(mode)

    predicted, source = predict_end_to_end(
        cleaned, prefer_albatross=prefer_albatross, warn_on_fallback=warn_on_fallback
    )
    requested = predicted * factor

    # A chain of n residues has n-1 bonds, so it cannot span more than (n-1) * 3.8 A no
    # matter what is asked of it. max_expansion on a long charged sequence gets close
    # enough to matter, and the pre-rewrite code would spin its full retry budget against
    # an impossible target and then return degenerate coordinates.
    ceiling = _MAX_CONTOUR_FRACTION * contour_length(n) if n > 1 else 0.0
    if n > 1 and requested > ceiling:
        if not clamp:
            raise UnsatisfiableTargetError(
                f"Mode {mode!r} asks for an end-to-end distance no {n}-residue chain can reach.",
                target=requested,
                achievable=ceiling,
            )
        if warn_on_clamp:
            # Say so. The predicted distance is the default and the thing a user believes they
            # are getting, so silently substituting a different one changes the science without
            # telling them. Clamping is still the right behaviour -- the alternative is refusing
            # to build a region at all -- but it has to be visible.
            #
            # Two distinct cases, and conflating them would make the message wrong. A request
            # above the contour length is unreachable at any conformation. A request between the
            # ceiling and the contour length IS reachable, but only as an almost-straight rod
            # with essentially one conformation, which is not a plausible chain and leaves the
            # sampler no freedom to avoid a clash.
            full = contour_length(n)
            why = (
                f"no chain of {n} residues can reach it: {n - 1} bonds of "
                f"{CA_CA_BOND_LENGTH:.2f} A span at most {full:.1f} A fully extended"
                if requested > full
                else (
                    f"a chain of {n} residues can only reach it by straightening almost "
                    f"completely, since its fully extended span is {full:.1f} A -- leaving no "
                    f"conformational freedom and no way to avoid a clash"
                )
            )
            warnings.warn(
                f"The requested end-to-end distance for this {n}-residue region is "
                f"{requested:.1f} A, and {why}. Building to {ceiling:.1f} A instead, which is "
                f"{_MAX_CONTOUR_FRACTION:.0%} of the contour length and the largest distance "
                f"that still admits a physically plausible conformation. The region will come "
                f"out more extended than a real chain of this composition would be, because the "
                f"request itself was not physically achievable.",
                UserWarning,
                stacklevel=2,
            )
        return DimensionTarget(
            end_to_end=ceiling,
            n_residues=n,
            mode=mode,
            factor=factor,
            predicted=predicted,
            source=source,
            clamped=True,
            max_achievable=ceiling,
        )

    return DimensionTarget(
        end_to_end=requested,
        n_residues=n,
        mode=mode,
        factor=factor,
        predicted=predicted,
        source=source,
        clamped=False,
        max_achievable=ceiling,
    )


def albatross_short_sequence_threshold() -> int:
    """Return the length below which ALBATROSS needs its scaled network.

    Read from sparrow when installed so the value tracks upstream, with DODO's own
    constant as the fallback. Exposed mainly so a caller can report it; DODO does not
    branch on it, for the reason given in :func:`predict_end_to_end`.
    """
    try:
        from sparrow.data.configs import MIN_LENGTH_ALBATROSS_RE_RG
    except ImportError:
        return ALBATROSS_MIN_LENGTH
    return int(MIN_LENGTH_ALBATROSS_RE_RG)
