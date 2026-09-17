"""Putting back the residues a structure did not model.

A crystal or cryo-EM structure contains what could be resolved. What could not be resolved is,
overwhelmingly, disordered -- which is to say it is exactly what DODO exists to rebuild, and it
is not in the file. Handing DODO such a structure and a reference sequence lets it answer the
question the file cannot: **every residue the reference says is there and the structure does not
show is treated as disordered, and rebuilt.**

The reference comes from whichever of these is available, in this order:

* a FASTA the caller supplies (see :mod:`dodo.io.fasta`), which is the general case;
* ``Chain.full_sequence``, already parsed from ``SEQRES`` or mmCIF ``_entity_poly``.

Where the inserted residues go
------------------------------
Insertion happens **after** region identification, and that ordering is load-bearing. An
inserted residue has no coordinates; it gets a placeholder alpha carbon so the arrays stay
well-formed, and that placeholder is fiction. Burial is scored from coordinates, so letting
region identification see the placeholders would let fiction decide which residues are folded --
a straight line drawn through a protein core scores as buried, and the region DODO most needs to
rebuild would come out classified as structure.

So the regions are decided on the observed geometry, and the inserted residues are then spliced
into that decision by where they fall:

* strictly inside a folded domain -- an unresolved surface loop -- becomes a **loop** of that
  domain. It is rebuilt between two fixed anchors, and the domain stays one rigid body, which is
  what keeps the two halves of a split domain from being pulled apart;
* inside or against a disordered region extends that region;
* between two folded domains with nothing between them becomes a new connecting IDR;
* past either end of the chain extends, or creates, a terminal IDR.

What is not guessed
-------------------
A gap whose flanking anchors are further apart than the missing residues can physically span is
detected before anything is inserted and reported, not attempted: no conformation exists, and
inserting the residues would only turn a statement about the input into a build failure. And a
region whose build fails afterwards has its placeholder residues **removed from the output**
rather than written -- there is no input geometry to fall back on, so writing them would be
writing a straight line and calling it a structure.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from itertools import pairwise

import numpy as np

from ..constants import CA_CA_BOND_LENGTH
from ..engines.walk import max_reach, min_reach
from ..exceptions import InvalidRegionError
from ..structure import Chain, Domain, DomainKind, Span, Structure

__all__ = [
    "ChainInsertion",
    "InsertionReport",
    "SequenceMapping",
    "hold_observed_residues",
    "insert_unmodelled_residues",
    "map_to_reference",
    "skipped_for_length",
]


def hold_observed_residues(structure: Structure) -> list[str]:
    """Re-tile every chain so that anything observed is fixed and only insertions are rebuilt.

    The rule for an experimental structure, in Ryan's words: a residue resolved in the map is
    "sufficiently static to be resolved", so it stays exactly where the experiment put it --
    **however few of them there are**. It is not a claim that the residues form a folded domain.
    What is missing is what is dynamic, and that is what gets built.

    Region identification cannot express this, because it answers a different question. It scores
    burial and needs :data:`~dodo.constants.MIN_FOLDED_DOMAIN_LENGTH` residues before it will call
    anything folded -- so a chain with 19 residues resolved gets *no* folded domain, the whole
    chain becomes one anchor-free region, and DODO regenerates it and lands it on a centroid
    computed from placeholder coordinates. Measured on 7R5J: 48 of the 56 Nup98 copies model only
    19 residues each, and every one of them was rebuilt from scratch and flung up to 855 A out of
    the pore, leaving 38 detached islands in the output.

    So this replaces the burial-based tiling outright: each run of observed residues becomes a
    folded domain, each run of inserted residues an IDR anchored to whatever flanks it. Length
    does not enter into it.

    Returns notes describing what it did, for the report.
    """
    notes: list[str] = []
    held = rebuilt = 0
    for chain in structure.chains:
        inserted = structure.inserted[chain.span.slice]
        domains: list[Domain] = []
        for start, stop, is_inserted in _alternating_runs(inserted, chain.span.start):
            if is_inserted:
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
                rebuilt += stop - start
            else:
                domains.append(
                    Domain(structure=structure, span=Span(start, stop), kind=DomainKind.FOLDED)
                )
                held += stop - start
        chain.domains = domains
        chain.validate_domains()

    notes.append(
        f"experimental input: {held} observed residue(s) are held exactly as deposited and "
        f"{rebuilt} inserted residue(s) are rebuilt around them. A residue resolved in the "
        f"experiment is treated as static however short the stretch it belongs to, so region "
        f"identification's folded/disordered call is not used here."
    )
    if not rebuilt:
        notes.append(
            "nothing was rebuilt: every residue in this structure was observed. Supply the "
            "full-length sequence with fasta= (or --fasta) to build the residues the file does "
            "not model, or pass units='predicted' to re-sample its disordered regions."
        )
    return notes


def _alternating_runs(mask: np.ndarray, offset: int) -> list[tuple[int, int, bool]]:
    """Split a boolean mask into maximal runs, as ``(start, stop, value)`` absolute indices."""
    if mask.size == 0:
        return []
    edges = np.flatnonzero(np.diff(mask)) + 1
    bounds = np.concatenate([[0], edges, [mask.size]])
    return [
        (offset + int(a), offset + int(b), bool(mask[a])) for a, b in pairwise(bounds)
    ]


def skipped_for_length(structure: Structure, span: Span, min_length: int) -> bool:
    """Whether a region is left alone for being too short to be worth rebuilding.

    The one predicate every caller must agree on -- the builder that skips the region, the
    pre-marker that puts it into the obstacle set before anything is built, and the rigid-unit
    rule that treats an unrebuilt link as part of the rigid body. When they disagreed, a region
    was skipped by one and treated as buildable by another.

    **An inserted region is never skipped, however short.** The skip exists because a short
    region's input coordinates are good enough to leave alone; an inserted region has no input
    coordinates at all, only a placeholder, so skipping it would leave a straight line in the
    output. Two inserted residues between two folded domains have to be built or removed, and
    building them is trivially possible.
    """
    if bool(structure.inserted[span.slice].any()):
        return False
    return (span.stop - span.start) < min_length

#: One-letter code to the three-letter name written for an inserted residue.
#:
#: Inverted from :data:`~dodo.constants.THREE_TO_ONE`, preferring the standard name where
#: several map to the same letter -- an inserted residue is one the structure never showed, so
#: there is no evidence it is selenomethionine rather than methionine, and writing MSE would be
#: inventing a modification.
_ONE_TO_THREE: dict[str, str] = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE", "G": "GLY", "H": "HIS",
    "I": "ILE", "K": "LYS", "L": "LEU", "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN",
    "R": "ARG", "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
    "O": "PYL", "U": "SEC", "X": "UNK",
}

#: Fraction of a chain's observed residues that may disagree with the reference before the
#: reference is rejected for that chain.
#:
#: CHOICE. A few disagreements are ordinary -- engineered point mutations, a modified residue
#: read as X, a cloning artefact. A tenth of the chain disagreeing means the reference is a
#: different protein, and rebuilding a chain against the wrong sequence produces a confident,
#: wrong answer, which is worse than declining.
MAX_REFERENCE_MISMATCH_FRACTION = 0.1

#: Shifts searched when trying to explain author numbering as a constant offset into the
#: reference. CHOICE: wide enough for a cleaved tag or a construct numbered from its own start,
#: narrow enough that the scan stays trivial.
_MAX_NUMBERING_SHIFT = 5000

#: Tie-break weight, per residue of disagreement with author numbering, in the alignment.
#: Small enough that it never outweighs a single matched residue; large enough to choose
#: between two otherwise equal embeddings, which is the only thing it is for.
_POSITION_TIE_BREAK = 1e-3


@dataclass(frozen=True, slots=True)
class SequenceMapping:
    """Where each observed residue of a chain sits in its reference sequence."""

    #: Index into the reference for each observed residue, strictly increasing.
    positions: np.ndarray
    #: How the mapping was found: ``"author numbering"`` or ``"alignment"``.
    method: str
    #: Observed residues whose identity disagrees with the reference at their mapped position.
    mismatches: int
    #: Constant offset such that ``reference_index == residue_number - 1 - shift``, when the
    #: author numbering explains the mapping. ``None`` when it does not.
    shift: int | None

    def __len__(self) -> int:
        return int(self.positions.shape[0])


@dataclass(frozen=True, slots=True)
class ChainInsertion:
    """What was inserted into one chain."""

    chain_id: str
    observed: int
    reference: int
    inserted: int
    leading: int
    trailing: int
    internal_gaps: int
    method: str
    mismatches: int
    #: Gaps that were NOT filled because no chain of that length could bridge their anchors,
    #: as ``(first_residue_number, last_residue_number, separation, violated_bound)``.
    impossible_gaps: tuple[tuple[int, int, float, float], ...] = ()

    def __str__(self) -> str:
        if not self.inserted:
            return f"chain {self.chain_id}: nothing missing ({self.observed} residues)"
        parts = []
        if self.leading:
            parts.append(f"{self.leading} at the N terminus")
        if self.internal_gaps:
            parts.append(f"{self.internal_gaps} internal gap(s)")
        if self.trailing:
            parts.append(f"{self.trailing} at the C terminus")
        detail = f" ({', '.join(parts)})" if parts else ""
        mismatch = (
            f", {self.mismatches} mismatch(es) against the reference" if self.mismatches else ""
        )
        skipped = (
            f"; {len(self.impossible_gaps)} gap(s) left empty as unbridgeable"
            if self.impossible_gaps
            else ""
        )
        return (
            f"chain {self.chain_id}: {self.observed} observed of {self.reference} in the "
            f"reference, {self.inserted} inserted{detail} by {self.method}{mismatch}{skipped}"
        )


@dataclass(slots=True)
class InsertionReport:
    """Outcome of filling in a structure's unmodelled residues."""

    chains: list[ChainInsertion] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def n_inserted(self) -> int:
        """Total residues inserted across every chain."""
        return sum(c.inserted for c in self.chains)

    @property
    def impossible_gaps(self) -> list[tuple[str, tuple[int, int, float, float]]]:
        """Gaps left empty because no chain of that length could bridge them."""
        return [(c.chain_id, gap) for c in self.chains for gap in c.impossible_gaps]

    def summary(self) -> str:
        """Multi-line human-readable summary."""
        lines = [
            f"{self.n_inserted} unmodelled residue(s) inserted across "
            f"{sum(1 for c in self.chains if c.inserted)} chain(s)"
        ]
        lines += [f"  {c}" for c in self.chains if c.inserted or c.mismatches]
        for chain_id, (first, last, separation, bound) in self.impossible_gaps:
            limit = "at most" if separation > bound else "at least"
            lines.append(
                f"  chain {chain_id} residues {first}-{last}: NOT inserted -- the flanking "
                f"residues are {separation:.1f} A apart and this many residues span {limit} "
                f"{bound:.1f} A, so no conformation exists"
            )
        lines += [f"  note: {n}" for n in self.notes]
        return "\n".join(lines)


def _numbering_mapping(
    observed: str, residue_numbers: np.ndarray, reference: str
) -> SequenceMapping | None:
    """Try to explain the mapping as ``reference_index = residue_number - 1 - shift``.

    The common case by a wide margin: author numbering usually *is* the position in the
    deposited sequence, sometimes offset by a constant when a tag was cleaved. When it works it
    gives an exact answer with no alignment ambiguity at all.

    The shift is not searched exhaustively. Only a shift that lands the first observed residue
    on a matching reference residue can possibly work, so the candidates are read straight off
    the positions where that residue occurs -- a few dozen, against the ten thousand a blind
    scan would test. Each candidate is then rejected on a handful of sampled residues before
    anything compares the whole chain.
    """
    numbers = np.asarray(residue_numbers, dtype=np.int64)
    if numbers.size == 0 or np.any(np.diff(numbers) <= 0):
        # Non-increasing author numbering cannot be an index into anything.
        return None
    reference_array = np.frombuffer(reference.encode("ascii", "replace"), dtype="S1")
    observed_array = np.frombuffer(observed.encode("ascii", "replace"), dtype="S1")

    lowest, highest = int(numbers[0]), int(numbers[-1])
    low_shift = max(highest - len(reference), -_MAX_NUMBERING_SHIFT)
    high_shift = min(lowest - 1, _MAX_NUMBERING_SHIFT)
    if high_shift < low_shift:
        return None

    occurrences = np.flatnonzero(reference_array == observed_array[0])
    candidates = np.unique(np.concatenate([[0], lowest - 1 - occurrences]))
    candidates = candidates[(candidates >= low_shift) & (candidates <= high_shift)]

    probe = np.unique(np.linspace(0, numbers.size - 1, num=min(8, numbers.size)).astype(np.int64))
    for shift in candidates:
        positions = numbers - 1 - int(shift)
        if positions[0] < 0 or positions[-1] >= len(reference):
            continue
        if not np.array_equal(reference_array[positions[probe]], observed_array[probe]):
            continue
        if np.array_equal(reference_array[positions], observed_array):
            return SequenceMapping(
                positions=positions.astype(np.int64),
                method="author numbering",
                mismatches=0,
                shift=int(shift),
            )
    return None


def _alignment_mapping(
    observed: str, residue_numbers: np.ndarray, reference: str
) -> SequenceMapping:
    """Best embedding of the observed sequence into the reference, allowing mismatches.

    Dynamic programming over ``(observed residue, reference position)``: every observed residue
    must be placed, positions must strictly increase, and reference positions may be skipped
    freely -- a skipped position is precisely a residue the structure did not model. Cost is one
    per mismatched residue, plus a thousandth per residue of disagreement with the author
    numbering, which decides between embeddings that are otherwise exactly as good. Without that
    tie-break a repetitive sequence -- and disordered regions are full of them -- would put a run
    of glycines at the leftmost place it fits rather than where the numbering says it belongs.

    Vectorised over the reference, so the Python loop runs once per observed residue: the prefix
    minimum along the reference is what makes "any earlier position" an O(1) lookup instead of a
    second loop.
    """
    n, m = len(observed), len(reference)
    if n == 0 or n > m:
        raise InvalidRegionError(
            f"Cannot map {n} observed residues into a reference of {m}: the reference must be "
            f"at least as long as what was modelled."
        )
    reference_array = np.frombuffer(reference.encode("ascii", "replace"), dtype="S1")
    observed_array = np.frombuffer(observed.encode("ascii", "replace"), dtype="S1")

    numbers = np.asarray(residue_numbers, dtype=np.int64)
    guess = _numbering_shift_guess(numbers, n, m)
    # Normalised so the tie-break summed over the WHOLE chain stays below one mismatch. It has
    # to be able to choose between two embeddings that match equally well, and it must never be
    # able to buy one by accepting a mismatch.
    tie_weight = _POSITION_TIE_BREAK / max(1, n * m)

    infinity = np.float64(np.inf)
    best = np.full(m, infinity)
    choice = np.zeros((n, m), dtype=np.int32)

    for i in range(n):
        cost = np.where(reference_array == observed_array[i], 0.0, 1.0)
        cost += tie_weight * np.abs(np.arange(m) - (numbers[i] - 1 - guess))
        if i == 0:
            current = cost.copy()
            choice[0] = -1
        else:
            running = np.minimum.accumulate(best)
            argmin = _running_argmin(best)
            current = np.full(m, infinity)
            current[1:] = cost[1:] + running[:-1]
            choice[i, 1:] = argmin[:-1]
            choice[i, 0] = -1
        best = current

    end = int(np.argmin(best))
    positions = np.empty(n, dtype=np.int64)
    for i in range(n - 1, -1, -1):
        positions[i] = end
        end = int(choice[i, end])
    mismatches = int(np.sum(reference_array[positions] != observed_array))
    return SequenceMapping(
        positions=positions, method="alignment", mismatches=mismatches, shift=None
    )


def _running_argmin(values: np.ndarray) -> np.ndarray:
    """Index of the minimum of ``values[:k + 1]`` for each ``k``. The argmin of a prefix min."""
    indices = np.arange(values.shape[0])
    improved = values < np.minimum.accumulate(
        np.concatenate([[np.inf], values[:-1]])
    )
    return np.maximum.accumulate(np.where(improved, indices, -1))


def _numbering_shift_guess(numbers: np.ndarray, n: int, m: int) -> int:
    """Guess a constant offset between author numbering and reference position.

    Only a tie-break hint, so being approximately right is the whole requirement. Clamped so the
    implied window fits inside the reference; otherwise every residue is pulled against the same
    wall and the hint stops discriminating between embeddings, which is all it is for.
    """
    if numbers.size == 0:
        return 0
    return int(max(int(numbers[0]) - 1, int(numbers[-1]) - m))


def map_to_reference(
    observed: str, residue_numbers: Sequence[int] | np.ndarray, reference: str
) -> SequenceMapping:
    """Map a chain's observed residues onto its reference sequence.

    Tries the author numbering first, because when it works it is exact and free, and falls back
    to alignment when it does not.
    """
    numbers = np.asarray(residue_numbers, dtype=np.int64)
    exact = _numbering_mapping(observed, numbers, reference)
    if exact is not None:
        return exact
    return _alignment_mapping(observed, numbers, reference)


# ----------------------------------------------------------------------------------
# Building the filled-in structure
# ----------------------------------------------------------------------------------


@dataclass(slots=True)
class _ChainPlan:
    """Per-chain bookkeeping between deciding what to insert and building the arrays."""

    chain: Chain
    reference: str
    mapping: SequenceMapping
    #: Reference positions to insert, sorted.
    positions: list[int]
    #: Author residue numbers for the inserted residues, parallel to ``positions``.
    numbers: list[int]
    renumbered: bool
    impossible: list[tuple[int, int, float, float]]


def _gaps(mapping: SequenceMapping, reference_length: int) -> list[tuple[int, int, int, int]]:
    """Find runs of unmodelled reference positions, as ``(before, after, start, stop)``.

    ``before`` and ``after`` are observed residue indices flanking the run (``-1`` and ``n`` at
    the termini), and ``[start, stop)`` are reference positions.
    """
    positions = mapping.positions
    runs: list[tuple[int, int, int, int]] = []
    if positions[0] > 0:
        runs.append((-1, 0, 0, int(positions[0])))
    for k in range(len(positions) - 1):
        lo, hi = int(positions[k]), int(positions[k + 1])
        if hi > lo + 1:
            runs.append((k, k + 1, lo + 1, hi))
    if positions[-1] < reference_length - 1:
        runs.append((len(positions) - 1, len(positions), int(positions[-1]) + 1, reference_length))
    return runs


def _plan_chain(
    structure: Structure, chain: Chain, reference: str, report: InsertionReport
) -> _ChainPlan | None:
    """Decide what to insert into one chain, or return None to leave it alone."""
    observed = chain.sequence
    numbers = structure.residue_number[chain.span.slice]
    if not observed:
        return None
    if len(observed) > len(reference):
        report.notes.append(
            f"chain {chain.chain_id}: the reference is {len(reference)} residues but "
            f"{len(observed)} are modelled, so it is not this chain's sequence; left as-is."
        )
        return None

    mapping = map_to_reference(observed, numbers, reference)
    if mapping.mismatches > MAX_REFERENCE_MISMATCH_FRACTION * len(observed):
        report.notes.append(
            f"chain {chain.chain_id}: {mapping.mismatches} of {len(observed)} modelled residues "
            f"disagree with the reference, which is more than "
            f"{MAX_REFERENCE_MISMATCH_FRACTION:.0%}; the reference is a different sequence, so "
            f"the chain was left as-is."
        )
        return None

    ca = structure.ca_xyz[chain.span.slice]
    positions: list[int] = []
    impossible: list[tuple[int, int, float, float]] = []
    n_reference = len(reference)

    for before, after, start, stop in _gaps(mapping, n_reference):
        length = stop - start
        if before >= 0 and after < len(observed):
            # An internal gap has two fixed anchors, so it is bridgeable or it is not, and that
            # is decided by geometry alone. Saying so here beats inserting residues that no
            # conformation can place and calling it a build failure later.
            separation = float(np.linalg.norm(ca[after] - ca[before]))
            ceiling = max_reach(length + 1)
            floor = min_reach(length + 1)
            if separation > ceiling or separation < floor:
                impossible.append(
                    (
                        int(numbers[before]),
                        int(numbers[after]),
                        separation,
                        ceiling if separation > ceiling else floor,
                    )
                )
                continue
        positions.extend(range(start, stop))

    if not positions:
        return _empty_plan(chain, reference, mapping, impossible)

    numbers_for_inserted, renumbered = _number_inserted(mapping, numbers, positions, n_reference)
    return _ChainPlan(
        chain=chain,
        reference=reference,
        mapping=mapping,
        positions=positions,
        numbers=numbers_for_inserted,
        renumbered=renumbered,
        impossible=impossible,
    )


def _empty_plan(
    chain: Chain,
    reference: str,
    mapping: SequenceMapping,
    impossible: list[tuple[int, int, float, float]],
) -> _ChainPlan:
    return _ChainPlan(
        chain=chain,
        reference=reference,
        mapping=mapping,
        positions=[],
        numbers=[],
        renumbered=False,
        impossible=impossible,
    )


def _number_inserted(
    mapping: SequenceMapping,
    observed_numbers: np.ndarray,
    positions: Sequence[int],
    reference_length: int,
) -> tuple[list[int], bool]:
    """Choose author residue numbers for the inserted residues.

    Two residues of one chain must not share a ``(number, insertion code)`` pair -- readers group
    residues on exactly that, so a collision silently merges two residues into one. When the
    author numbering has room for the missing residues, which is the usual case because a gap in
    the numbering is *how* an unmodelled stretch is recorded, they take the numbers the file left
    free. When it does not, the whole chain is renumbered from the reference and the report says
    so, rather than producing a file whose numbering quietly disagrees with itself.
    """
    if mapping.shift is not None:
        # The numbering IS the reference position, so there is by construction room for
        # everything and the inserted residues take the numbers the gaps left free.
        return [p + 1 + mapping.shift for p in positions], False

    taken = {int(n) for n in observed_numbers}
    by_position = {int(p): int(o) for o, p in enumerate(mapping.positions)}
    first_position = int(mapping.positions[0])
    first_number = int(observed_numbers[0])
    numbers: list[int] = []
    fits = True
    for position in positions:
        # Walk back to the nearest observed residue and count forward from its number. Before
        # the first observed residue there is nothing to count forward from, so count BACK from
        # it instead -- counting forward from it would number the leading tail in reverse.
        anchor = position - 1
        while anchor >= 0 and anchor not in by_position:
            anchor -= 1
        if anchor >= 0:
            candidate = int(observed_numbers[by_position[anchor]]) + (position - anchor)
        else:
            candidate = first_number - (first_position - position)
        if candidate in taken:
            fits = False
            break
        taken.add(candidate)
        numbers.append(candidate)

    if fits:
        return numbers, False
    return [p + 1 for p in positions], True


def _placeholder_coordinates(
    structure: Structure, chain: Chain, mapping: SequenceMapping, positions: Sequence[int]
) -> np.ndarray:
    """Provisional alpha carbons for the inserted residues.

    Fiction, and treated as such everywhere downstream: an internal run is spread evenly along
    the straight line between its anchors, and a terminal run walks outward from the last
    observed residue at the CA-CA bond length. It exists so the arrays are well-formed and the
    rebuilt region has somewhere to start; nothing reads it as data, and a region that fails to
    build has these residues removed from the output rather than written.
    """
    ca = structure.ca_xyz[chain.span.slice]
    by_position = {int(p): int(o) for o, p in enumerate(mapping.positions)}
    highest = int(mapping.positions[-1])
    coordinates = np.empty((len(positions), 3), dtype=np.float64)

    for k, position in enumerate(positions):
        before = position - 1
        while before >= 0 and before not in by_position:
            before -= 1
        after = position + 1
        while after <= highest and after not in by_position:
            after += 1

        if before >= 0 and after <= highest:
            start, stop = ca[by_position[before]], ca[by_position[after]]
            fraction = (position - before) / (after - before)
            coordinates[k] = start + fraction * (stop - start)
        elif after <= highest:
            # Leading tail: walk outward from the first observed residue, away from the chain.
            first = by_position[after]
            direction = ca[first] - ca[min(first + 1, len(ca) - 1)]
            coordinates[k] = ca[first] + _unit(direction) * CA_CA_BOND_LENGTH * (after - position)
        elif before >= 0:
            last = by_position[before]
            direction = ca[last] - ca[max(last - 1, 0)]
            coordinates[k] = ca[last] + _unit(direction) * CA_CA_BOND_LENGTH * (position - before)
        else:  # pragma: no cover - a chain with no observed residue has no atoms to read
            coordinates[k] = 0.0
    return coordinates


def _unit(vector: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if norm < 1e-9:
        return np.array([1.0, 0.0, 0.0])
    return vector / norm


def insert_unmodelled_residues(
    structure: Structure,
    sequences: Mapping[str, str] | None = None,
    *,
    use_full_sequence: bool = True,
    on_chain_done: Callable[[int], None] | None = None,
) -> tuple[Structure, InsertionReport]:
    """Insert the residues a reference sequence says are missing from ``structure``.

    ``structure`` must already have regions assigned -- see the module docstring for why the
    order matters. The returned structure has the same regions, with the inserted residues
    spliced into them, and :attr:`Structure.inserted` marking every residue that was added.

    Parameters
    ----------
    structure
        A structure with regions assigned. Not modified.
    sequences
        Chain id to reference sequence. Chains absent from the mapping fall back to
        ``Chain.full_sequence`` when ``use_full_sequence`` is set, and are otherwise untouched.
    use_full_sequence
        Fall back to the deposited sequence parsed from ``SEQRES`` or mmCIF ``_entity_poly``.
    on_chain_done
        Called with ``1`` after each chain is planned, for a caller showing progress. The
        alignment fallback is the expensive path and it runs per chain, so on an assembly this
        is the difference between a moving count and a frozen one.

    Returns
    -------
    tuple
        The new structure and a report of what was inserted, per chain.
    """
    sequences = dict(sequences or {})
    report = InsertionReport()

    plans: list[_ChainPlan | None] = []
    for chain in structure.chains:
        reference = sequences.get(chain.chain_id)
        if reference is None and use_full_sequence:
            reference = chain.full_sequence
        if not reference:
            plans.append(None)
            if on_chain_done is not None:
                on_chain_done(1)
            continue
        if not chain.domains:
            raise InvalidRegionError(
                f"Chain {chain.chain_id!r} has no assigned regions. Regions must be identified "
                f"on the OBSERVED geometry before unmodelled residues are inserted, so that "
                f"placeholder coordinates cannot influence which residues are called folded."
            )
        plans.append(_plan_chain(structure, chain, reference, report))
        if on_chain_done is not None:
            on_chain_done(1)

    for chain, plan in zip(structure.chains, plans, strict=True):
        if plan is None:
            continue
        report.chains.append(
            ChainInsertion(
                chain_id=chain.chain_id,
                observed=len(chain.span),
                reference=len(plan.reference),
                inserted=len(plan.positions),
                leading=sum(1 for p in plan.positions if p < int(plan.mapping.positions[0])),
                trailing=sum(1 for p in plan.positions if p > int(plan.mapping.positions[-1])),
                internal_gaps=sum(
                    1
                    for before, after, _, _ in _gaps(plan.mapping, len(plan.reference))
                    if before >= 0 and after < len(plan.mapping)
                ),
                method=plan.mapping.method,
                mismatches=plan.mapping.mismatches,
                impossible_gaps=tuple(plan.impossible),
            )
        )
        if plan.renumbered:
            report.notes.append(
                f"chain {chain.chain_id}: the author numbering had no room for the missing "
                f"residues, so the whole chain was renumbered 1..{len(plan.reference)} from the "
                f"reference. Residue numbers in the output will not match the input for this "
                f"chain."
            )

    if not any(plan is not None and plan.positions for plan in plans):
        report.notes.append("no unmodelled residues to insert")
        return structure.copy(), report

    return _build_filled_structure(structure, plans), report


def _build_filled_structure(
    structure: Structure, plans: Sequence[_ChainPlan | None]
) -> Structure:
    """Assemble the new arrays, then splice the regions onto them."""
    # A per-new-residue source: an index into the old residue arrays, or -1 for an insertion.
    residue_source: list[np.ndarray] = []
    inserted_flags: list[np.ndarray] = []
    new_numbers: list[np.ndarray] = []
    new_names: list[np.ndarray] = []
    placeholder_xyz: list[np.ndarray] = []
    # old residue index -> new residue index, for remapping every span.
    remap = np.full(structure.n_residues, -1, dtype=np.int64)
    cursor = 0

    for chain, plan in zip(structure.chains, plans, strict=True):
        old_indices = np.arange(chain.span.start, chain.span.stop)
        if plan is None or not plan.positions:
            residue_source.append(old_indices)
            inserted_flags.append(np.zeros(old_indices.size, dtype=bool))
            new_numbers.append(structure.residue_number[chain.span.slice])
            new_names.append(structure.residue_name[chain.span.slice])
            remap[old_indices] = np.arange(cursor, cursor + old_indices.size)
            cursor += old_indices.size
            continue

        order = np.argsort(
            np.concatenate([plan.mapping.positions, np.asarray(plan.positions, dtype=np.int64)]),
            kind="stable",
        )
        source = np.concatenate(
            [old_indices, np.full(len(plan.positions), -1, dtype=np.int64)]
        )[order]
        flags = np.concatenate(
            [np.zeros(old_indices.size, dtype=bool), np.ones(len(plan.positions), dtype=bool)]
        )[order]
        numbers = np.concatenate(
            [
                structure.residue_number[chain.span.slice],
                np.asarray(plan.numbers, dtype=np.int64),
            ]
        )[order]
        if plan.renumbered:
            all_positions = np.concatenate(
                [plan.mapping.positions, np.asarray(plan.positions, dtype=np.int64)]
            )[order]
            numbers = all_positions + 1
        names = np.concatenate(
            [
                structure.residue_name[chain.span.slice],
                np.asarray(
                    [_ONE_TO_THREE.get(plan.reference[p], "UNK") for p in plan.positions],
                    dtype="<U3",
                ),
            ]
        )[order]

        residue_source.append(source)
        inserted_flags.append(flags)
        new_numbers.append(numbers)
        new_names.append(names)
        placeholder_xyz.append(
            _placeholder_coordinates(structure, chain, plan.mapping, plan.positions)
        )
        observed_new = cursor + np.flatnonzero(~flags)
        remap[old_indices] = observed_new
        cursor += source.size

    source = np.concatenate(residue_source)
    inserted = np.concatenate(inserted_flags)
    n_residues = source.size

    old_counts = np.diff(structure.residue_atom_offsets)
    counts = np.where(inserted, 1, old_counts[np.maximum(source, 0)])
    offsets = np.empty(n_residues + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(counts, out=offsets[1:])
    n_atoms = int(offsets[-1])

    atom_source = np.full(n_atoms, -1, dtype=np.int64)
    kept = np.flatnonzero(~inserted)
    if kept.size:
        starts_new = offsets[kept]
        starts_old = structure.residue_atom_offsets[source[kept]]
        lengths = counts[kept]
        # One flat gather: for every kept atom, its index in the old array.
        atom_source[_ranges(starts_new, lengths)] = _ranges(starts_old, lengths)

    xyz = np.empty((n_atoms, 3), dtype=np.float64)
    atom_name = np.empty(n_atoms, dtype=structure.atom_name.dtype)
    element = np.empty(n_atoms, dtype=structure.element.dtype)
    real = atom_source >= 0
    xyz[real] = structure.xyz[atom_source[real]]
    atom_name[real] = structure.atom_name[atom_source[real]]
    element[real] = structure.element[atom_source[real]]
    xyz[~real] = np.concatenate(placeholder_xyz) if placeholder_xyz else np.zeros((0, 3))
    atom_name[~real] = "CA"
    element[~real] = "C"

    filled = Structure(
        xyz=xyz,
        atom_name=atom_name,
        element=element,
        residue_index=np.repeat(np.arange(n_residues), counts),
        residue_name=np.concatenate(new_names).astype(structure.residue_name.dtype),
        residue_number=np.concatenate(new_numbers),
        insertion_code=np.where(
            inserted, "", structure.insertion_code[np.maximum(source, 0)]
        ).astype(structure.insertion_code.dtype),
        b_factor=np.where(inserted, 0.0, structure.b_factor[np.maximum(source, 0)]),
        # Occupancy zero is the crystallographic convention for a position that was not
        # observed, and it is the one field a viewer will show without being asked.
        occupancy=np.where(inserted, 0.0, structure.occupancy[np.maximum(source, 0)]),
        chain_index=np.repeat(
            np.arange(len(structure.chains)), [len(part) for part in residue_source]
        ),
        residue_atom_offsets=offsets,
        inserted=inserted,
        source=structure.source,
        experimental_method=structure.experimental_method,
        notes=list(structure.notes),
    )

    _splice_regions(structure, filled, plans, remap)
    filled.validate()
    for chain in filled.chains:
        chain.validate_domains()
    return filled


def _ranges(starts: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Concatenate ``range(s, s + n)`` for each ``(s, n)``, without a Python loop."""
    total = int(lengths.sum())
    if total == 0:
        return np.zeros(0, dtype=np.int64)
    offsets = np.zeros(total, dtype=np.int64)
    boundaries = np.cumsum(lengths)[:-1]
    offsets[boundaries] = 1
    group = np.cumsum(offsets)
    within = np.arange(total) - np.repeat(np.concatenate([[0], boundaries]), lengths)
    flat: np.ndarray = starts[group] + within
    return flat


def _splice_regions(
    original: Structure,
    filled: Structure,
    plans: Sequence[_ChainPlan | None],
    remap: np.ndarray,
) -> None:
    """Carry the region assignment across, absorbing each inserted run into the right region."""
    for chain_index, (chain, plan) in enumerate(zip(original.chains, plans, strict=True)):
        span_start = int(np.flatnonzero(filled.chain_index == chain_index)[0])
        span_stop = int(np.flatnonzero(filled.chain_index == chain_index)[-1]) + 1
        new_chain = Chain(
            structure=filled,
            span=Span(span_start, span_stop),
            chain_id=chain.chain_id,
            uniprot_id=chain.uniprot_id,
            full_sequence=plan.reference if plan is not None else chain.full_sequence,
        )
        filled.chains.append(new_chain)

        ordered = sorted(chain.domains, key=lambda d: d.span.start)
        # Every domain keeps its first and last OBSERVED residue, and swallows whatever was
        # inserted strictly inside it. A folded domain therefore stays one domain when an
        # unresolved loop is filled in -- which is what keeps its two halves from being pulled
        # apart, since a rigid unit is a set of whole domains.
        domains: list[Domain] = []
        for domain in ordered:
            start = int(remap[domain.span.start])
            stop = int(remap[domain.span.stop - 1]) + 1
            loops = [
                Span(int(remap[loop.start]), int(remap[loop.stop - 1]) + 1)
                for loop in domain.loops
            ]
            if domain.kind is DomainKind.FOLDED:
                # An unresolved stretch strictly inside a folded domain is a loop: fixed
                # geometry on both sides, so its span is dictated rather than predicted. A run
                # touching either boundary is not -- that is a tail, and _absorb_between hands
                # it to the neighbouring disordered region instead.
                for run_start, run_stop in _runs(filled.inserted, start, stop):
                    if run_start > start and run_stop < stop:
                        loops.append(Span(run_start, run_stop))
            domains.append(
                Domain(
                    structure=filled,
                    span=Span(start, stop),
                    kind=domain.kind,
                    # Anchors, not just bounds. A loop is rebuilt between two FIXED residues,
                    # and a Span without them reads as a region with one free end -- which the
                    # builder rejects outright as "not a loop". Merging first, then anchoring,
                    # so a merged pair anchors on the residues that actually flank the union.
                    loops=tuple(
                        Span(loop.start, loop.stop, n_anchor=loop.start - 1, c_anchor=loop.stop)
                        for loop in _merge_spans(loops)
                    ),
                    label=domain.label,
                )
            )

        domains = _absorb_between(filled, new_chain, domains)
        for index, domain in enumerate(domains):
            if domain.kind is not DomainKind.IDR:
                continue
            domains[index] = Domain(
                structure=filled,
                span=Span(
                    domain.span.start,
                    domain.span.stop,
                    n_anchor=domain.span.start - 1 if domain.span.start > span_start else None,
                    c_anchor=domain.span.stop if domain.span.stop < span_stop else None,
                ),
                kind=DomainKind.IDR,
                label=domain.label,
            )
        new_chain.domains = domains


def _runs(mask: np.ndarray, start: int, stop: int) -> list[tuple[int, int]]:
    """Contiguous True runs of ``mask`` within ``[start, stop)``, as absolute indices."""
    window = mask[start:stop]
    if not window.any():
        return []
    padded = np.concatenate([[False], window, [False]])
    edges = np.flatnonzero(padded[1:] != padded[:-1])
    return [
        (start + int(a), start + int(b)) for a, b in zip(edges[::2], edges[1::2], strict=True)
    ]


def _merge_spans(spans: Sequence[Span]) -> list[Span]:
    """Sort and merge overlapping or abutting spans, so loops stay disjoint and ordered."""
    if not spans:
        return []
    ordered = sorted(spans, key=lambda s: s.start)
    merged = [ordered[0]]
    for span in ordered[1:]:
        last = merged[-1]
        if span.start <= last.stop:
            merged[-1] = Span(last.start, max(last.stop, span.stop))
        else:
            merged.append(span)
    return merged


def _absorb_between(filled: Structure, chain: Chain, domains: list[Domain]) -> list[Domain]:
    """Give every residue between (and around) the domains to a disordered region.

    After the domains have taken their observed residues, whatever is left is inserted geometry
    lying between two regions or past a chain end. It extends an adjacent IDR when there is one,
    and becomes a new IDR when there is not -- which is the case of two folded domains whose
    connecting linker was never modelled at all.
    """
    result: list[Domain] = []
    cursor = chain.span.start

    def new_idr(start: int, stop: int) -> Domain:
        return Domain(structure=filled, span=Span(start, stop), kind=DomainKind.IDR)

    for domain in domains:
        if domain.span.start > cursor:
            if domain.kind is DomainKind.IDR:
                domain = Domain(
                    structure=filled,
                    span=Span(cursor, domain.span.stop),
                    kind=DomainKind.IDR,
                    label=domain.label,
                )
            elif result and result[-1].kind is DomainKind.IDR:
                previous = result[-1]
                result[-1] = Domain(
                    structure=filled,
                    span=Span(previous.span.start, domain.span.start),
                    kind=DomainKind.IDR,
                    label=previous.label,
                )
            else:
                result.append(new_idr(cursor, domain.span.start))
        result.append(domain)
        cursor = domain.span.stop

    if cursor < chain.span.stop:
        if result and result[-1].kind is DomainKind.IDR:
            previous = result[-1]
            result[-1] = Domain(
                structure=filled,
                span=Span(previous.span.start, chain.span.stop),
                kind=DomainKind.IDR,
                label=previous.label,
            )
        else:
            result.append(new_idr(cursor, chain.span.stop))
    return result
