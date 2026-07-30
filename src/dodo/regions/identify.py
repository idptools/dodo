"""Assign folded domains, IDRs and loops to a structure.

Given a structure, decide which residues belong to a folded domain, which form an
intrinsically disordered region, and which are rebuildable loops sitting inside a folded
domain. Everything downstream depends on getting this right: a folded domain misread as an
IDR gets replaced by a random walk, and an IDR swallowed into a folded domain never gets
rebuilt at all -- so the output silently presents an AlphaFold artifact as structure.

Strategies
----------
Three ways to make the call, behind one enum:

``CONTACT``
    Geometric burial from the coordinates (:mod:`dodo.regions.contact`). Works on any
    structure, including experimental ones with no confidence metric.
``PLDDT``
    AlphaFold's own per-residue confidence, which is already sitting in the B-factor
    column. Free, and it is the model's own assessment of where it is guessing. The
    pre-rewrite code parsed it into a field and then never looked at it, choosing instead to
    re-derive confidence geometrically from the coordinates AlphaFold produced *because* it
    was uncertain.
``METAPREDICT``
    Sequence-only disorder prediction. Needs no structure at all, so it is the strategy for
    building from sequence, and it is far faster than the geometric route.

:data:`Strategy.AUTO` picks ``PLDDT`` when the structure looks like an AlphaFold model and
``CONTACT`` otherwise.

The two merge bugs
------------------
v1 and the first v2 attempt shared a block-merging routine with two reproducible defects,
both caused by index-range arithmetic over the block list:

1. ``for fd_ind in range(0, len(bounds) - 1)`` never executes for a *single* candidate
   block, so a protein with one clean folded domain and no internal sub-threshold dip came
   out with zero folded domains -- the entire chain classified as one IDR.
2. The guard ``if fd_ind < len(bounds) - 2`` meant the gap before the *last* block was never
   tested, so an IDR lying between the last two blocks was absorbed into the final folded
   domain. Bounds ``[[40,100],[103,163],[283,433]]`` collapsed to a single 393-residue
   "folded domain" containing a 120-residue IDR.

:func:`merge_blocks` here is a fold over consecutive pairs with no index arithmetic and no
cardinality special cases, so both are gone by construction rather than by patch.

Critically, these had to be fixed *together with* the score normalization in
:mod:`dodo.regions.contact`. The old raw score was noisy enough to fragment a domain into
dozens of blocks, which is what kept execution away from the single-block and last-gap paths
where the bugs lived. Improving the metric is what exposes them.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np

from ..constants import (
    CONTACT_SCORE_THRESHOLD,
    MAX_INTERNAL_GAP,
    MIN_FOLDED_DOMAIN_LENGTH,
    MIN_FOLDED_SEED_RUN,
    MIN_IDR_LENGTH,
    MIN_LOOP_LENGTH,
    PLDDT_DISORDER_THRESHOLD,
)
from ..exceptions import InvalidRegionError, MissingDependencyError
from ..structure import Chain, Domain, DomainKind, Span, Structure
from .contact import contact_profile, is_loop_like, loop_contact_counts

__all__ = [
    "RegionAssignment",
    "Strategy",
    "assign_regions",
    "assign_regions_from_spec",
    "find_runs",
    "merge_blocks",
]


class Strategy(str, Enum):
    """How to decide which residues are folded."""

    #: Geometric burial from coordinates.
    CONTACT = "contact"
    #: AlphaFold pLDDT, read from the B-factor column.
    PLDDT = "plddt"
    #: Sequence-only disorder prediction via metapredict.
    METAPREDICT = "metapredict"
    #: pLDDT for AlphaFold models, contact scoring otherwise.
    AUTO = "auto"


@dataclass(frozen=True, slots=True)
class RegionAssignment:
    """The outcome of region identification for one chain, with its evidence.

    The score profile and threshold are retained deliberately. Region identification is a
    judgement call with tuned constants, and a user who disagrees with a boundary needs to
    be able to audit *why* it was drawn there rather than re-run with different numbers and
    hope. The pre-rewrite code returned bare index pairs.

    Attributes
    ----------
    chain_id
        The chain these regions belong to.
    domains
        The assigned domains, in residue order, tiling the chain with no gaps or overlaps.
    strategy
        Which strategy produced the call.
    score
        Per-residue score that was thresholded, ``(n_chain_residues,)``.
    threshold
        The cutoff applied to ``score``.
    folded_mask
        Which residues the threshold called folded, before merging and length filtering.
    notes
        Anything a caller should know: domains rejected for length, a fully disordered
        chain, a strategy that had to be substituted.
    """

    chain_id: str
    domains: tuple[Domain, ...]
    strategy: Strategy
    score: np.ndarray
    threshold: float
    folded_mask: np.ndarray
    notes: tuple[str, ...] = ()

    @property
    def n_folded(self) -> int:
        """Number of folded domains assigned."""
        return sum(1 for d in self.domains if d.kind is DomainKind.FOLDED)

    @property
    def n_idrs(self) -> int:
        """Number of IDRs assigned."""
        return sum(1 for d in self.domains if d.kind is DomainKind.IDR)

    @property
    def fully_disordered(self) -> bool:
        """True if no folded domain was found anywhere in the chain."""
        return self.n_folded == 0

    def describe(self) -> str:
        """Return a one-line human-readable summary, 1-based for display.

        Output is 1-based because that is what users see in PDB files and papers, while
        everything internal is 0-based positional. That conversion happens here and
        nowhere else.
        """
        parts = []
        for domain in self.domains:
            label = "FD" if domain.kind is DomainKind.FOLDED else "IDR"
            start = domain.span.start + 1
            stop = domain.span.stop
            piece = f"{label} {start}-{stop}"
            if domain.loops:
                loops = ", ".join(f"{loop.start + 1}-{loop.stop}" for loop in domain.loops)
                piece += f" (loops: {loops})"
            parts.append(piece)
        return f"chain {self.chain_id}: " + "; ".join(parts)


def find_runs(mask: np.ndarray, *, min_length: int = 1) -> list[tuple[int, int]]:
    """Find maximal runs of True in a boolean mask, as half-open ``(start, stop)`` pairs.

    Parameters
    ----------
    mask
        Boolean array.
    min_length
        Discard runs shorter than this.

    Returns
    -------
    list[tuple[int, int]]
        Runs in order. Empty if ``mask`` is empty or all False.

    Examples
    --------
    >>> find_runs(np.array([False, True, True, False, True]))
    [(1, 3), (4, 5)]
    >>> find_runs(np.array([True, True, False, True]), min_length=2)
    [(0, 2)]
    """
    if mask.size == 0:
        return []
    # Pad with False at both ends so every run has an explicit rising and falling edge;
    # this removes the "run touches the array boundary" special cases that index-based
    # implementations get wrong.
    padded = np.concatenate([[False], mask.astype(bool), [False]])
    edges = np.diff(padded.astype(np.int8))
    starts = np.flatnonzero(edges == 1)
    stops = np.flatnonzero(edges == -1)
    return [(int(a), int(b)) for a, b in zip(starts, stops, strict=True) if b - a >= min_length]


def merge_blocks(blocks: list[tuple[int, int]], *, max_gap: int) -> list[tuple[int, int]]:
    """Merge blocks separated by a gap of at most ``max_gap`` residues.

    A fold over consecutive pairs. There is no index arithmetic over the block list and no
    special case for its length, which is what makes this correct for zero, one, and n
    blocks alike -- and which is what the two pre-rewrite merge bugs both came from.

    Parameters
    ----------
    blocks
        Half-open ``(start, stop)`` pairs in ascending order.
    max_gap
        Largest gap that still counts as internal to one domain. A gap of exactly this
        many residues merges; one residue more does not.

    Returns
    -------
    list[tuple[int, int]]
        Merged blocks, in order.

    Examples
    --------
    A single block survives. The v1 loop produced nothing here, classifying a whole
    single-domain protein as disordered:

    >>> merge_blocks([(40, 100)], max_gap=25)
    [(40, 100)]

    The gap before the last block is tested like any other. v1 never tested it, so this
    collapsed into one 393-residue "folded domain" containing a 120-residue IDR:

    >>> merge_blocks([(40, 100), (103, 163), (283, 433)], max_gap=25)
    [(40, 163), (283, 433)]

    >>> merge_blocks([], max_gap=25)
    []
    """
    if not blocks:
        return []
    if max_gap < 0:
        raise ValueError(f"max_gap must be non-negative, got {max_gap}")

    merged: list[tuple[int, int]] = [blocks[0]]
    for start, stop in blocks[1:]:
        previous_start, previous_stop = merged[-1]
        if start < previous_stop:
            raise InvalidRegionError(
                f"Blocks must be sorted and non-overlapping; ({start}, {stop}) overlaps "
                f"({previous_start}, {previous_stop})."
            )
        if start - previous_stop <= max_gap:
            merged[-1] = (previous_start, stop)
        else:
            merged.append((start, stop))
    return merged


def _folded_mask_from_contact(
    structure: Structure, chain: Chain, threshold: float
) -> tuple[np.ndarray, np.ndarray]:
    """Score by geometric burial. Returns (score, folded_mask) for the chain's residues."""
    profile = contact_profile(structure)
    score = profile.smoothed[chain.span.slice]
    return score, score >= threshold


def _folded_mask_from_plddt(
    structure: Structure, chain: Chain, threshold: float
) -> tuple[np.ndarray, np.ndarray]:
    """Score by AlphaFold pLDDT, read from the B-factor column."""
    score = structure.b_factor[chain.span.slice].astype(np.float64)
    return score, score >= threshold


def _folded_mask_from_metapredict(chain: Chain, threshold: float) -> tuple[np.ndarray, np.ndarray]:
    """Score by sequence-only disorder prediction.

    Returns *order* (1 - disorder) so that higher means more folded, matching the sign
    convention of the other two strategies. Keeping every strategy's score oriented the
    same way is what lets one thresholding path serve all three.
    """
    try:
        import metapredict
    except ImportError as exc:
        raise MissingDependencyError(
            package="metapredict",
            purpose="Sequence-based region identification",
            extra="predictors",
        ) from exc

    disorder = np.asarray(metapredict.predict_disorder(chain.sequence), dtype=np.float64)
    if disorder.shape[0] != len(chain):
        raise InvalidRegionError(
            f"metapredict returned {disorder.shape[0]} scores for a {len(chain)}-residue "
            f"chain. This is a version incompatibility; DODO expects one score per residue."
        )
    order = 1.0 - disorder
    return order, order >= threshold


def _looks_like_alphafold(structure: Structure, chain: Chain) -> bool:
    """Guess whether a chain's B-factors are pLDDT rather than crystallographic B-factors.

    pLDDT is a percentage in [0, 100] and, in a real model, spans a wide range with most
    residues confident. Crystallographic B-factors are unbounded above and typically sit in
    the 5-60 range with a different shape. The discriminating check is the upper bound
    together with some spread: a B-factor column that never exceeds 100 and reaches above 70
    somewhere is almost certainly pLDDT.

    Deliberately conservative -- a wrong guess here silently changes which residues get
    rebuilt -- so :data:`Strategy.AUTO` reports which strategy it chose.
    """
    values = structure.b_factor[chain.span.slice]
    if values.size == 0:
        return False
    if np.all(values == values[0]):
        # A constant column carries no information either way; do not claim pLDDT.
        return False
    return bool(values.max() <= 100.0 and values.max() >= 70.0 and values.min() >= 0.0)


def _detect_loops(
    structure: Structure,
    span: Span,
    loop_counts: np.ndarray,
    *,
    min_loop_length: int,
) -> tuple[Span, ...]:
    """Find rebuildable loops inside a folded domain.

    A loop is a run of low-CA-contact residues at least ``min_loop_length`` long, lying
    strictly inside the domain. Runs touching either boundary are excluded: those are the
    domain's own termini, and a "loop" that is anchored on only one side is a tail, which
    the IDR path already handles.

    Note the comparison is ``>=`` against ``min_loop_length``. The pre-rewrite code used
    strict ``>``, making the effective minimum 11 while the documentation said 10.
    """
    if len(span) == 0:
        return ()
    local = loop_counts[span.slice]
    candidates = find_runs(is_loop_like(local), min_length=min_loop_length)

    loops: list[Span] = []
    for start, stop in candidates:
        absolute_start = span.start + start
        absolute_stop = span.start + stop
        # Strictly interior: a loop needs a fixed residue on each side to anchor to.
        if absolute_start <= span.start or absolute_stop >= span.stop:
            continue
        loops.append(
            Span(
                absolute_start,
                absolute_stop,
                n_anchor=absolute_start - 1,
                c_anchor=absolute_stop,
            )
        )
    return tuple(loops)


def _absorb_short_gaps(
    blocks: list[tuple[int, int]],
    *,
    chain_start: int,
    chain_stop: int,
    min_idr_length: int,
) -> tuple[list[tuple[int, int]], list[str]]:
    """Extend folded blocks over disordered gaps too short to rebuild.

    A 2-residue disordered stretch has no meaningful polymer statistics to impose -- there is
    nothing to be gained by replacing two residues with a random walk, and doing so throws
    away real coordinates. Such a gap is absorbed into the folded block beside it, so those
    residues simply keep their input coordinates.

    Done *before* tiling rather than during it, so that the tiling stays a strict alternation
    of folded and disordered regions. Absorbing during the tiling walk produced two adjacent
    folded domains, which is a valid tiling but a confusing one for the builder to consume.

    Returns the adjusted blocks and a note for each absorbed gap, since silently reclassifying
    residues is exactly the kind of invisible decision this package tries not to make.
    """
    if not blocks:
        return [], []

    notes: list[str] = []

    def note(start: int, stop: int, where: str) -> None:
        notes.append(
            f"residues {start + 1}-{stop} ({stop - start} residue(s)) scored as disordered "
            f"but are shorter than the {min_idr_length}-residue minimum, so they were "
            f"absorbed into the adjacent folded domain ({where}) and keep their original "
            f"coordinates"
        )

    adjusted = [blocks[0]]
    for start, stop in blocks[1:]:
        previous_start, previous_stop = adjusted[-1]
        gap = start - previous_stop
        if 0 < gap < min_idr_length:
            note(previous_stop, start, "interior gap")
            adjusted[-1] = (previous_start, stop)
        else:
            adjusted.append((start, stop))

    # Terminal tails: a short N- or C-terminal disordered stretch is absorbed the same way.
    first_start, first_stop = adjusted[0]
    leading = first_start - chain_start
    if 0 < leading < min_idr_length:
        note(chain_start, first_start, "N-terminal tail")
        adjusted[0] = (chain_start, first_stop)

    last_start, last_stop = adjusted[-1]
    trailing = chain_stop - last_stop
    if 0 < trailing < min_idr_length:
        note(last_stop, chain_stop, "C-terminal tail")
        adjusted[-1] = (last_start, chain_stop)

    return adjusted, notes


def _tile_chain(
    structure: Structure,
    chain: Chain,
    folded_blocks: list[tuple[int, int]],
    loop_counts: np.ndarray,
    *,
    min_loop_length: int,
) -> list[Domain]:
    """Turn folded blocks into a complete, gapless tiling of the chain.

    Everything not in a folded block is an IDR. Spans are half-open and adjacent regions do
    not share a residue; instead each IDR records the flanking fixed residues as explicit
    anchors. The pre-rewrite convention had adjacent regions *overlap* by one residue --
    because that residue serves as the geometric anchor -- and then needed three different
    conventions inside one function to undo the overlap.

    Assumes short gaps have already been absorbed by :func:`_absorb_short_gaps`, so every
    gap here becomes a real IDR.
    """
    domains: list[Domain] = []
    cursor = chain.span.start

    def add_idr(start: int, stop: int) -> None:
        if stop <= start:
            return
        domains.append(
            Domain(
                structure=structure,
                span=Span(
                    start,
                    stop,
                    n_anchor=start - 1 if start > chain.span.start else None,
                    c_anchor=stop if stop < chain.span.stop else None,
                ),
                kind=DomainKind.IDR,
            )
        )

    for start, stop in folded_blocks:
        add_idr(cursor, start)
        loops = _detect_loops(
            structure, Span(start, stop), loop_counts, min_loop_length=min_loop_length
        )
        domains.append(
            Domain(
                structure=structure,
                span=Span(start, stop),
                kind=DomainKind.FOLDED,
                loops=loops,
            )
        )
        cursor = stop

    add_idr(cursor, chain.span.stop)
    return domains


def assign_regions(
    structure: Structure,
    *,
    strategy: Strategy | str = Strategy.AUTO,
    threshold: float | None = None,
    max_internal_gap: int = MAX_INTERNAL_GAP,
    min_folded_length: int = MIN_FOLDED_DOMAIN_LENGTH,
    min_loop_length: int = MIN_LOOP_LENGTH,
    min_idr_length: int = MIN_IDR_LENGTH,
    min_seed_run: int = MIN_FOLDED_SEED_RUN,
) -> list[RegionAssignment]:
    """Assign folded domains, IDRs and loops to every chain of a structure.

    Mutates ``structure``: each chain's ``domains`` list is replaced with the assignment.
    The returned objects carry the evidence behind each call.

    Parameters
    ----------
    structure
        The structure to annotate.
    strategy
        Which signal to threshold. See :class:`Strategy`.
    threshold
        Cutoff on the strategy's score. Defaults to the strategy's own tuned value:
        :data:`~dodo.constants.CONTACT_SCORE_THRESHOLD` for contacts,
        :data:`~dodo.constants.PLDDT_DISORDER_THRESHOLD` for pLDDT, and 0.5 for metapredict.
    max_internal_gap
        Longest run of non-folded residues that stays inside one folded domain. Split from
        ``min_folded_length``: the pre-rewrite code used a single ``gap_thresh`` knob for
        both, which are unrelated quantities.
    min_folded_length
        Shortest run of folded residues that counts as a domain.
    min_loop_length
        Shortest rebuildable loop inside a folded domain.
    min_idr_length
        Shortest IDR worth rebuilding. Shorter disordered stretches keep their input
        coordinates, and that is recorded in the assignment's notes.
    min_seed_run
        Consecutive above-threshold residues needed to seed a candidate block.

    Returns
    -------
    list[RegionAssignment]
        One per chain, in chain order.

    Raises
    ------
    MissingDependencyError
        If ``strategy`` is ``METAPREDICT`` and metapredict is not installed.
    InvalidRegionError
        If the resulting tiling is inconsistent, which would be a bug here rather than bad
        input -- the check exists so such a bug surfaces immediately instead of as a
        confusing failure during the build.
    """
    strategy = Strategy(strategy)
    assignments: list[RegionAssignment] = []

    # One structure-wide contact pass, reused across chains. Burial is inherently a
    # whole-structure property: a residue at a chain-chain interface is buried by the
    # partner chain, and scoring each chain in isolation would call those interfaces
    # disordered.
    loop_counts = loop_contact_counts(structure)

    for chain in structure.chains:
        notes: list[str] = []
        resolved = strategy
        if resolved is Strategy.AUTO:
            resolved = (
                Strategy.PLDDT if _looks_like_alphafold(structure, chain) else Strategy.CONTACT
            )
            reason = (
                "B-factors look like pLDDT"
                if resolved is Strategy.PLDDT
                else "B-factors do not look like pLDDT"
            )
            notes.append(f"strategy auto-selected as {resolved.value} ({reason})")

        if resolved is Strategy.PLDDT:
            cutoff = PLDDT_DISORDER_THRESHOLD if threshold is None else threshold
            score, folded_mask = _folded_mask_from_plddt(structure, chain, cutoff)
        elif resolved is Strategy.METAPREDICT:
            cutoff = 0.5 if threshold is None else threshold
            score, folded_mask = _folded_mask_from_metapredict(chain, cutoff)
        else:
            cutoff = CONTACT_SCORE_THRESHOLD if threshold is None else threshold
            score, folded_mask = _folded_mask_from_contact(structure, chain, cutoff)

        # Seed -> merge -> length filter. Offsets are chain-local until here.
        seeds = find_runs(folded_mask, min_length=min_seed_run)
        merged = merge_blocks(seeds, max_gap=max_internal_gap)
        kept = [(a, b) for a, b in merged if b - a >= min_folded_length]
        rejected = [(a, b) for a, b in merged if b - a < min_folded_length]
        for a, b in rejected:
            notes.append(
                f"candidate folded block at residues {chain.span.start + a + 1}-"
                f"{chain.span.start + b} ({b - a} residues) was below the "
                f"{min_folded_length}-residue minimum and treated as disordered"
            )

        # Lift chain-local offsets to structure-wide residue indices.
        absolute = [(chain.span.start + a, chain.span.start + b) for a, b in kept]

        absolute, absorb_notes = _absorb_short_gaps(
            absolute,
            chain_start=chain.span.start,
            chain_stop=chain.span.stop,
            min_idr_length=min_idr_length,
        )
        notes.extend(absorb_notes)

        domains = _tile_chain(
            structure,
            chain,
            absolute,
            loop_counts,
            min_loop_length=min_loop_length,
        )

        if not absolute:
            notes.append("no folded domain found; the whole chain is treated as disordered")

        chain.domains = domains
        # Fail fast on an inconsistent tiling: this is our bug if it fires, and catching it
        # here beats a confusing geometry error three phases later.
        chain.validate_domains()

        assignments.append(
            RegionAssignment(
                chain_id=chain.chain_id,
                domains=tuple(domains),
                strategy=resolved,
                score=score,
                threshold=cutoff,
                folded_mask=folded_mask,
                notes=tuple(notes),
            )
        )

    return assignments


def assign_regions_from_spec(
    structure: Structure,
    spec: dict[str, list[tuple[str, int, int]]],
) -> list[RegionAssignment]:
    """Assign regions from an explicit caller-supplied specification.

    The manual override, for when a user knows the domain boundaries and does not want them
    inferred. The pre-rewrite equivalent was undocumented, called two methods that did not
    exist, and had demonstrably never been executed.

    Parameters
    ----------
    structure
        The structure to annotate.
    spec
        Chain id to a list of ``(kind, start, stop)`` tuples, where ``kind`` is ``"folded"``
        or ``"idr"`` and the bounds are **1-based inclusive** residue positions, matching
        what a user reads off a PDB file or a paper. They are converted to 0-based half-open
        internally.

    Returns
    -------
    list[RegionAssignment]
        One per specified chain.

    Raises
    ------
    InvalidRegionError
        If a chain id is unknown, a kind is unrecognized, bounds are out of range, or the
        regions do not tile the chain without gaps or overlaps. All four were accepted
        silently by the pre-rewrite code and surfaced much later as unrelated errors.
    """
    by_id = {chain.chain_id: chain for chain in structure.chains}
    assignments: list[RegionAssignment] = []

    for chain_id, entries in spec.items():
        chain = by_id.get(chain_id)
        if chain is None:
            raise InvalidRegionError(
                f"Chain {chain_id!r} is not in this structure. Available: {sorted(by_id)}."
            )

        domains: list[Domain] = []
        for kind_name, start_1based, stop_1based in sorted(entries, key=lambda e: e[1]):
            try:
                kind = DomainKind(kind_name.lower())
            except ValueError:
                raise InvalidRegionError(
                    f"Unknown region kind {kind_name!r}. Use 'folded' or 'idr'."
                ) from None
            # 1-based inclusive -> 0-based half-open, offset into the structure.
            start = chain.span.start + start_1based - 1
            stop = chain.span.start + stop_1based
            if start < chain.span.start or stop > chain.span.stop:
                raise InvalidRegionError(
                    f"Region {kind_name} {start_1based}-{stop_1based} lies outside chain "
                    f"{chain_id!r}, which has {len(chain)} residues."
                )
            domains.append(
                Domain(
                    structure=structure,
                    span=Span(
                        start,
                        stop,
                        n_anchor=start - 1 if start > chain.span.start else None,
                        c_anchor=stop if stop < chain.span.stop else None,
                    ),
                    kind=kind,
                )
            )

        chain.domains = domains
        chain.validate_domains()
        assignments.append(
            RegionAssignment(
                chain_id=chain_id,
                domains=tuple(domains),
                strategy=Strategy.CONTACT,
                score=np.zeros(len(chain), dtype=np.float64),
                threshold=float("nan"),
                folded_mask=np.zeros(len(chain), dtype=bool),
                notes=("regions supplied explicitly by the caller; nothing was inferred",),
            )
        )

    return assignments
