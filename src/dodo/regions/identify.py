"""Assign folded domains, IDRs and loops to a structure.

Given a structure, decide which residues belong to a folded domain, which form an
intrinsically disordered region, and which are rebuildable loops sitting inside a folded
domain. Everything downstream depends on getting this right: a folded domain misread as an
IDR gets replaced by a random walk, and an IDR swallowed into a folded domain never gets
rebuilt at all -- so the output silently presents an AlphaFold artifact as structure.

Strategies
----------
``DENSITY`` (the default)
    DODO's original all-atom density score: all-atom pairs within 8 A per residue, thresholded
    at 480. This is the method the package was built and validated on, and the author reports
    it draws better boundaries than sequence-based disorder predictors. It is reimplemented
    over a KD-tree rather than changed -- same numbers, 10.1 s down to 7 ms on a 1086-residue
    model.
``CONTACT``
    An alternative CA-only burial score. Composition-free (every residue has exactly one CA)
    and invariant to whether side chains are modelled, which the density score is not. Useful
    for comparison and for CA-only input, but it is not the validated method.
``PLDDT``
    AlphaFold's own per-residue confidence, from the B-factor column. Cheap, and it is the
    model's own account of where it was guessing -- but the density method reportedly beats
    disorder-based calls, so this is an explicit opt-in rather than a default.

:data:`Strategy.AUTO` resolves to ``DENSITY`` for all-atom input and ``CONTACT`` for CA-only
input, where a pair count cannot be compared against the tuned threshold.

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
    CA_CONTACT_SCORE_THRESHOLD,
    CONTACT_SCORE_THRESHOLD,
    MAX_INTERNAL_GAP,
    MIN_FOLDED_DOMAIN_LENGTH,
    MIN_FOLDED_SEED_RUN,
    MIN_IDR_LENGTH,
    MIN_LOOP_LENGTH,
    PLDDT_DISORDER_THRESHOLD,
)
from ..exceptions import InvalidRegionError
from ..structure import Chain, Domain, DomainKind, Span, Structure
from .contact import contact_profile, density_profile, is_loop_like, loop_contact_counts

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

    #: DODO's original all-atom density score. The default, and the method the package was
    #: built and validated on -- the author reports it outperforms sequence-based disorder
    #: predictors at drawing region boundaries.
    DENSITY = "density"
    #: Alternative CA-only burial score: composition-free and invariant to whether side chains
    #: are modelled, but not the validated method. Available for comparison.
    CONTACT = "contact"
    #: AlphaFold pLDDT, read from the B-factor column.
    PLDDT = "plddt"
    #: Resolves to :attr:`DENSITY` when the structure has side chains to measure, and to
    #: :attr:`CONTACT` when it is CA-only (where a pair count is not comparable to the tuned
    #: threshold). Deliberately does NOT prefer pLDDT: the density method was found to work
    #: better, and pLDDT is available as an explicit choice.
    AUTO = "auto"
    #: Do not identify anything: use the regions already attached to the structure.
    #:
    #: This is the granular-control path. Assign regions however you like -- run
    #: :func:`assign_regions` and edit the result, call :func:`assign_regions_from_spec`, or
    #: construct :class:`~dodo.structure.Domain` objects directly -- then hand the structure to
    #: :func:`~dodo.rebuild` with this strategy and it will build exactly those regions.
    #:
    #: It replaces v1's ``regions_dict=`` parameter, which took a parallel stringly-typed
    #: description of the structure that had to be validated and kept in step with the real
    #: object. Passing the actual objects means there is one representation, it is already
    #: validated, and it carries the score profile and threshold that produced it.
    PRESET = "preset"


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


def _folded_mask_from_density(
    structure: Structure, chain: Chain, threshold: float
) -> tuple[np.ndarray, np.ndarray]:
    """Score with DODO's original all-atom density metric. Returns (score, folded_mask)."""
    profile = density_profile(structure)
    score = profile.smoothed[chain.span.slice]
    return score, score >= threshold


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


def _detect_loops(
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
        loops = _detect_loops(Span(start, stop), loop_counts, min_loop_length=min_loop_length)
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
        :data:`~dodo.constants.PLDDT_DISORDER_THRESHOLD` for pLDDT.
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
    InvalidRegionError
        If the resulting tiling is inconsistent, which would be a bug here rather than bad
        input -- the check exists so such a bug surfaces immediately instead of as a
        confusing failure during the build.
    """
    strategy = Strategy(strategy)
    assignments: list[RegionAssignment] = []

    if strategy is Strategy.PRESET:
        # Nothing to identify: the caller has already decided. Report what is there and leave
        # every Domain object exactly as given, so a hand-tuned boundary survives verbatim.
        for chain in structure.chains:
            if not chain.domains:
                raise InvalidRegionError(
                    f"strategy={Strategy.PRESET.value!r} builds the regions already attached to "
                    f"the structure, but chain {chain.chain_id} has none. Assign them first -- "
                    f"with assign_regions(), assign_regions_from_spec(), or by constructing "
                    f"Domain objects -- or choose a strategy that identifies them."
                )
            folded = np.zeros(len(chain.span), dtype=bool)
            for domain in chain.domains:
                if domain.kind is DomainKind.FOLDED:
                    start = domain.span.start - chain.span.start
                    folded[start : start + len(domain.span)] = True
            assignments.append(
                RegionAssignment(
                    chain_id=chain.chain_id,
                    domains=tuple(chain.domains),
                    strategy=Strategy.PRESET,
                    # No score was computed and no threshold was applied. NaN says that
                    # honestly; a zero would read as a real measurement.
                    score=np.full(len(chain.span), np.nan),
                    threshold=float("nan"),
                    folded_mask=folded,
                    notes=("regions supplied by the caller; none were identified",),
                )
            )
        return assignments

    # One structure-wide contact pass, reused across chains. Burial is inherently a
    # whole-structure property: a residue at a chain-chain interface is buried by the
    # partner chain, and scoring each chain in isolation would call those interfaces
    # disordered.
    loop_counts = loop_contact_counts(structure)

    for chain in structure.chains:
        notes: list[str] = []
        resolved = strategy
        if resolved is Strategy.AUTO:
            # The density metric needs side chains to count. On a CA-only structure a pair
            # count is not comparable to the tuned 480 threshold, so fall back to the CA-only
            # score there rather than silently mis-thresholding.
            has_side_chains = bool(np.any(~np.isin(structure.atom_name, ("N", "CA", "C", "O"))))
            resolved = Strategy.DENSITY if has_side_chains else Strategy.CONTACT
            # Only worth saying when auto did something you would not have guessed. Density on
            # all-atom input is the expected, validated path and is what almost every run takes,
            # so announcing it every time is noise that trains people to skip the notes.
            if resolved is not Strategy.DENSITY:
                notes.append(
                    f"strategy auto-selected as {resolved.value} -- the input has no side chains, "
                    "so the density threshold does not apply"
                )

        if resolved is Strategy.PLDDT:
            cutoff = PLDDT_DISORDER_THRESHOLD if threshold is None else threshold
            score, folded_mask = _folded_mask_from_plddt(structure, chain, cutoff)
        elif resolved is Strategy.CONTACT:
            cutoff = CA_CONTACT_SCORE_THRESHOLD if threshold is None else threshold
            score, folded_mask = _folded_mask_from_contact(structure, chain, cutoff)
        else:
            cutoff = CONTACT_SCORE_THRESHOLD if threshold is None else threshold
            score, folded_mask = _folded_mask_from_density(structure, chain, cutoff)

        # Seed -> merge -> length filter. Offsets are chain-local until here.
        seeds = find_runs(folded_mask, min_length=min_seed_run)
        merged = merge_blocks(seeds, max_gap=max_internal_gap)
        kept = [(a, b) for a, b in merged if b - a >= min_folded_length]
        rejected = [(a, b) for a, b in merged if b - a < min_folded_length]
        if rejected:
            # One summary line, not one per block. A structure with a dozen short folded-looking
            # patches used to print a dozen near-identical notes, which buries anything that
            # actually needs reading.
            lengths = ", ".join(str(b - a) for a, b in rejected)
            one = len(rejected) == 1
            notes.append(
                f"{len(rejected)} short {'stretch' if one else 'stretches'} ({lengths} residues) "
                f"scored as folded but came in under the {min_folded_length}-residue minimum for a "
                f"folded domain, so {'it is' if one else 'they are'} treated as disordered"
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
