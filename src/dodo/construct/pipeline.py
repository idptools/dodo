"""End-to-end structure rebuilding: the layer that wires everything together.

This is the one call most users want. Everything below it -- parsing, region identification,
dimension prediction, conformation generation, backbone placement, writing -- is separately
usable and separately tested, but assembling them by hand means knowing which anchors to pass
where, and getting that wrong is how the pre-rewrite code ended up grafting regions onto the
wrong residues.

Design notes worth knowing
--------------------------
**Anchors, including the outer ones.** A rebuilt region is closed onto the fixed residues on
either side of it. But constraining the *pseudo-angle* centred on an anchor also needs the
residue one step beyond it, because that angle is
``CA(anchor-1) - CA(anchor) - CA(first rebuilt)``. This function has the whole structure, so it
supplies those outer residues; a caller wiring the engine directly usually forgets, and the
engine then warns that 41% of such junctions fall outside the valid angle window.

**Models are independent conformers, and the folded domains do not move between them.** v1
placed folded domains once, outside the model loop, so every model in a multi-model run shared
one identical domain arrangement *and* essentially one end-to-end distance -- a
"pseudo-trajectory" that was an ensemble at fixed dimensions, which is the one thing an IDR
ensemble should not be. Here each model re-samples its own target from the physical
distribution, so the models genuinely differ.

**Failure is per-region and explicit.** A region that cannot be built is reported in
:attr:`RebuildReport.failures` with the reason, and its input coordinates are left in place.
Nothing is silently replaced with degenerate coordinates.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..constants import CA_CA_BOND_LENGTH, DEFAULT_MODE, MIN_IDR_LENGTH
from ..exceptions import BuildError, DodoError
from ..regions.identify import RegionAssignment, Strategy, assign_regions
from ..structure import Domain, DomainKind, Span, Structure
from .backbone import place_backbone_for_domain
from .dimensions import DimensionTarget, target_dimensions
from .place import DomainPlacement, reposition_folded_domains

__all__ = [
    "RebuildReport",
    "RegionOutcome",
    "build_from_sequence",
    "rebuild",
]


@dataclass(frozen=True, slots=True)
class RegionOutcome:
    """What happened to one region in one model."""

    model: int
    chain_id: str
    #: Residue span, in the 1-based inclusive numbering a user reads off a PDB file.
    residues: tuple[int, int]
    n_residues: int
    built: bool
    target: DimensionTarget | None = None
    achieved_end_to_end: float | None = None
    #: The target this model actually steered to. Equals ``target.end_to_end`` for a single
    #: model; for a multi-model run of a free-ended region it is this model's own draw from
    #: the physical distribution, which is why it can differ substantially from the mean.
    requested_end_to_end: float | None = None
    reason: str | None = None

    def __str__(self) -> str:
        where = (
            f"model {self.model} chain {self.chain_id} "
            f"residues {self.residues[0]}-{self.residues[1]}"
        )
        if not self.built:
            return f"{where}: NOT BUILT ({self.reason})"
        detail = ""
        if self.achieved_end_to_end is not None:
            asked = self.requested_end_to_end
            if asked is None and self.target is not None:
                asked = self.target.end_to_end
            if asked is not None:
                detail = f", end-to-end {self.achieved_end_to_end:.1f} A against {asked:.1f} A"
                if self.target is not None and abs(asked - self.target.end_to_end) > 0.5:
                    detail += f" (this model's draw; ensemble mean {self.target.end_to_end:.1f} A)"
        suffix = f" -- {self.reason}" if self.reason else ""
        return f"{where}: built{detail}{suffix}"


@dataclass(slots=True)
class RebuildReport:
    """Everything that happened during a rebuild.

    Returned alongside the models rather than printed, so a caller can act on it. The
    pre-rewrite code communicated exclusively through ``print()``, which meant a script could
    not tell whether a region had actually been rebuilt.
    """

    models: list[Structure] = field(default_factory=list)
    assignments: list[RegionAssignment] = field(default_factory=list)
    outcomes: list[RegionOutcome] = field(default_factory=list)
    #: Where every folded domain ended up. Folded domains are moved as rigid bodies, never
    #: rebuilt, so this records position and orientation changes only.
    placements: list[DomainPlacement] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def failures(self) -> list[RegionOutcome]:
        """Regions that could not be built. Their input coordinates were left in place."""
        return [o for o in self.outcomes if not o.built]

    @property
    def n_built(self) -> int:
        """How many region-model pairs were successfully rebuilt."""
        return sum(1 for o in self.outcomes if o.built)

    @property
    def ok(self) -> bool:
        """True if every region in every model was rebuilt."""
        return not self.failures

    def summary(self) -> str:
        """Multi-line human-readable summary."""
        lines = [
            f"{len(self.models)} model(s), {self.n_built}/{len(self.outcomes)} "
            f"region-model pair(s) rebuilt"
        ]
        lines += [f"  {a.describe()}" for a in self.assignments]
        moved = [p for p in self.placements if p.moved]
        if moved:
            lines.append(f"  {len(moved)} folded domain(s) repositioned:")
            lines += [f"    {p}" for p in moved]
        if self.failures:
            lines.append(f"  {len(self.failures)} failure(s):")
            lines += [f"    {f}" for f in self.failures]
        lines += [f"  note: {n}" for n in self.notes]
        return "\n".join(lines)


def _outer_anchor(structure: Structure, domain: Domain, *, side: str) -> np.ndarray | None:
    """CA of the residue one step beyond an anchor, or None if there is no such residue.

    Needed to constrain the pseudo-angle centred *on* the anchor. Without it that angle is
    unconstrained and unmeasured, which is where the pre-rewrite builder's 38-177 degree
    junction angles came from.
    """
    span = domain.span
    anchor = span.n_anchor if side == "n" else span.c_anchor
    if anchor is None:
        return None
    outer = anchor - 1 if side == "n" else anchor + 1
    chain = structure.chains[int(structure.chain_index[anchor])]
    if not (chain.span.start <= outer < chain.span.stop):
        return None
    coords: np.ndarray = structure.ca_xyz[outer]
    return coords


def _obstacles_for(structure: Structure, domain: Domain) -> np.ndarray | None:
    """Coordinates the region must avoid: every already-placed atom outside it.

    The engine wants a plain ``(n, 3)`` array, not a Structure, so this flattens the
    already-rebuilt geometry.

    Two exclusions matter. The region's own residues are excluded because they are what is
    being replaced. Its immediate anchor residues are excluded because the engine already
    treats them as bonded neighbours and holds the region to a clash distance against them --
    passing them as obstacles too would make the engine reject the only positions that can
    actually connect to them.
    """
    mask = structure.rebuilt_atom_mask()
    # Drop the region itself.
    mask[domain.atom_slice] = False
    # Drop the flanking anchors, which the engine handles as chain neighbours.
    for anchor in (domain.span.n_anchor, domain.span.c_anchor):
        if anchor is not None:
            mask[structure.atom_slice_for_residues(anchor, anchor + 1)] = False
    if not mask.any():
        return None
    coords: np.ndarray = structure.xyz[mask]
    return coords


def _drop_non_ca_from_rebuilt(structure: Structure) -> Structure:
    """Return a copy holding only alpha carbons in rebuilt regions, all atoms elsewhere.

    Rebuilt regions are CA-only. Their other atoms cannot be left in place: the region's alpha
    carbons have moved, and the rest have not, so every residue ends up split between two
    locations. Measured on p300, that split was ~93 A per residue -- the CA at its new position
    and N/C/O/side chain still on AlphaFold's original path. The writer then emits a CONECT
    record bonding N to CA across that gap, which renders as a long spurious straight line, and
    the orphaned atoms trail along the old path as disconnected dots. Both artefacts are visible
    in a viewer and neither is a rendering problem.

    This is also simply what DODO produces and always has: all atoms in folded domains, alpha
    carbons only in rebuilt regions. Placing real backbone atoms for those residues is the next
    priority, not this one.

    Loops inside folded domains count as rebuilt too, so they are stripped as well.
    """
    keep = np.ones(structure.n_atoms, dtype=bool)
    is_ca = structure.atom_name == "CA"

    def strip(start: int, stop: int) -> None:
        atoms = structure.atom_slice_for_residues(start, stop)
        keep[atoms] = is_ca[atoms]

    for domain in structure.domains:
        if domain.kind is DomainKind.IDR and domain.rebuilt:
            strip(domain.span.start, domain.span.stop)
        for loop in domain.loops:
            strip(loop.start, loop.stop)

    if keep.all():
        return structure
    return structure.select_atoms(keep)


def _build_region(
    structure: Structure,
    *,
    span: Span,
    sequence: str,
    target: DimensionTarget | None,
    model: int,
    rng: np.random.Generator,
    engine: object,
    min_length: int,
    label: str,
    requested_override: float | None = None,
    domain: Domain | None = None,
) -> RegionOutcome:
    """Build one region -- a loop, a connecting IDR or a terminal IDR -- and report the outcome.

    One function for all three because they differ only in what sets the target: a loop's span
    is dictated by the folded-domain geometry it bridges, while an IDR's comes from its predicted
    dimensions. Sharing the code means the anchor handling cannot drift between them, which is
    where the pre-rewrite version accumulated its junction defects.
    """
    from ..engines.base import IDRRequest

    chain = structure.chains[int(structure.chain_index[span.start])]
    residues = (
        int(structure.residue_number[span.start]),
        int(structure.residue_number[span.stop - 1]),
    )

    def outcome(**kwargs: object) -> RegionOutcome:
        return RegionOutcome(
            model=model,
            chain_id=chain.chain_id,
            residues=residues,
            n_residues=len(sequence),
            **kwargs,  # type: ignore[arg-type]
        )

    if len(sequence) < min_length:
        return outcome(
            built=False,
            reason=f"{label}: shorter than the {min_length}-residue minimum; left as-is",
        )

    n_anchor_xyz = structure.ca_xyz[span.n_anchor] if span.n_anchor is not None else None
    c_anchor_xyz = structure.ca_xyz[span.c_anchor] if span.c_anchor is not None else None

    # A loop's target is simply the distance it has to bridge. Anything else would be a
    # constraint we invented rather than one the structure imposes.
    if target is None:
        if n_anchor_xyz is None or c_anchor_xyz is None:
            return outcome(
                built=False,
                reason=f"{label}: needs two anchors but only has one; not a loop",
            )
        requested = float(np.linalg.norm(c_anchor_xyz - n_anchor_xyz))
    else:
        requested = requested_override if requested_override is not None else target.end_to_end

    try:
        result = engine.generate(  # type: ignore[attr-defined]
            IDRRequest(
                sequence=sequence,
                n_residues=len(sequence),
                target_end_to_end=requested,
                n_anchor_xyz=n_anchor_xyz,
                c_anchor_xyz=c_anchor_xyz,
                n_anchor_prev_xyz=_outer_ca(structure, span.n_anchor, step=-1),
                c_anchor_next_xyz=_outer_ca(structure, span.c_anchor, step=+1),
                n_conformations=1,
            ),
            obstacles=_obstacles_for_span(structure, span),
            rng=rng,
        )
    except DodoError as exc:
        return outcome(built=False, reason=f"{label}: {type(exc).__name__}: {exc}")

    if not bool(result.success[0]):
        return outcome(
            built=False, reason=f"{label}: the engine reported failure for this conformer"
        )

    structure.set_ca_xyz(np.arange(span.start, span.stop), result.ca_coords[0])
    if domain is not None:
        domain.rebuilt = True

    achieved = float(np.linalg.norm(result.ca_coords[0][-1] - result.ca_coords[0][0]))
    return outcome(
        built=True,
        target=target,
        achieved_end_to_end=achieved,
        requested_end_to_end=requested,
    )


def _outer_ca(structure: Structure, anchor: int | None, *, step: int) -> np.ndarray | None:
    """CA one step beyond an anchor, for constraining the junction angle centred on it."""
    if anchor is None:
        return None
    outer = anchor + step
    chain = structure.chains[int(structure.chain_index[anchor])]
    if not (chain.span.start <= outer < chain.span.stop):
        return None
    coords: np.ndarray = structure.ca_xyz[outer]
    return coords


def _obstacles_for_span(structure: Structure, span: Span) -> np.ndarray | None:
    """Already-placed atoms the region must avoid, excluding itself and its own anchors."""
    mask = structure.rebuilt_atom_mask()
    mask[structure.atom_slice_for_residues(span.start, span.stop)] = False
    for anchor in (span.n_anchor, span.c_anchor):
        if anchor is not None:
            mask[structure.atom_slice_for_residues(anchor, anchor + 1)] = False
    if not mask.any():
        return None
    coords: np.ndarray = structure.xyz[mask]
    return coords


def _rebuild_one_model(
    structure: Structure,
    *,
    model: int,
    mode: str,
    engine_name: str,
    rng: np.random.Generator,
    min_length: int,
    sidechains: bool | None,
    all_atom: bool,
    model_targets: dict[tuple[str, int, int], np.ndarray],
) -> list[RegionOutcome]:
    """Rebuild every loop and IDR of one model, in place. Returns per-region outcomes."""
    from ..engines.walk import SelfAvoidingWalk

    engine: object
    if engine_name == "walk":
        engine = SelfAvoidingWalk()
    elif engine_name == "starling":
        from ..engines.starling import StarlingEngine

        engine = StarlingEngine()
        if not engine.available():
            raise BuildError(
                "The starling engine was requested but is not available. Install it with "
                "pip install 'dodo[starling]', or pass engine='walk'."
            )
    else:
        raise ValueError(f"Unknown engine {engine_name!r}. Use 'walk' or 'starling'.")

    outcomes: list[RegionOutcome] = []
    # Folded domains are already positioned (step 3 ran before this), so they are obstacles
    # from the outset. Their atoms are never rebuilt.
    for domain in structure.folded_domains():
        domain.rebuilt = True

    # BUILD ORDER, and it is deliberate: loops, then connecting IDRs, then terminal IDRs.
    #
    # The ordering follows how constrained each region is. A loop is pinned at both ends inside
    # a single folded domain, so it has the least freedom and must be satisfied first. A
    # connecting IDR is pinned at both ends but between two domains that have already been
    # positioned to accommodate it. A terminal IDR is free at one end and can go almost
    # anywhere, so it is built last, into whatever space remains.
    #
    # Building in residue order instead -- which an earlier version of this pipeline did --
    # lets a floppy terminal tail occupy space that a tightly constrained loop then cannot
    # avoid.
    loops: list[tuple[Domain, Span]] = [
        (domain, loop) for domain in structure.folded_domains() for loop in domain.loops
    ]
    idrs = structure.idrs()
    connecting = [d for d in idrs if not d.span.is_terminal]
    terminal = [d for d in idrs if d.span.is_terminal]

    for parent, loop in loops:
        outcomes.append(
            _build_region(
                structure,
                span=loop,
                sequence=structure.sequence[loop.slice],
                # A loop gets NO dimension prediction. Its span is dictated by the folded
                # domain it bridges: the two anchors are fixed atoms of a domain that is not
                # ours to move, so the distance is already decided. Predicting an end-to-end
                # distance here and then failing to achieve it would be inventing a constraint
                # the geometry has already settled.
                target=None,
                model=model,
                rng=rng,
                engine=engine,
                min_length=min_length,
                label=f"loop in FD {parent.span.start + 1}-{parent.span.stop}",
            )
        )

    for domain in connecting + terminal:
        key = (
            structure.chains[int(structure.chain_index[domain.span.start])].chain_id,
            domain.span.start,
            domain.span.stop,
        )
        drawn = model_targets.get(key)
        outcomes.append(
            _build_region(
                structure,
                span=domain.span,
                sequence=domain.sequence,
                target=target_dimensions(domain.sequence, mode=mode)
                if len(domain.sequence) >= min_length
                else None,
                requested_override=float(drawn[model - 1]) if drawn is not None else None,
                model=model,
                rng=rng,
                engine=engine,
                min_length=min_length,
                label="connecting IDR" if not domain.span.is_terminal else "terminal IDR",
                domain=domain,
            )
        )

    return outcomes


def rebuild(
    source: str | Path | Structure,
    *,
    mode: str = DEFAULT_MODE,
    n_models: int = 1,
    strategy: Strategy | str = Strategy.AUTO,
    engine: str = "walk",
    all_atom: bool = False,
    sidechains: bool = False,
    min_length: int = MIN_IDR_LENGTH,
    seed: int | None = None,
) -> RebuildReport:
    """Rebuild the disordered regions of a structure.

    Parameters
    ----------
    source
        Path to a PDB or mmCIF file, or an already-parsed :class:`~dodo.structure.Structure`.
    mode
        Named build mode. A multiplier on the predicted end-to-end distance:
        ``super_compact``, ``compact``, ``normal``/``predicted``, ``expanded``,
        ``super_expanded``, ``max_expansion``.
    n_models
        Number of independent conformers to generate. Folded domains stay put; each model
        re-samples its own IDR dimensions from the physical distribution, so the models
        genuinely differ.
    strategy
        How to identify regions. See :class:`~dodo.regions.identify.Strategy`.
    engine
        ``"walk"`` (always available) or ``"starling"`` (needs ``pip install 'dodo[starling]'``).
    all_atom
        Place N, C and O for rebuilt regions. Off by default because the analytic placement
        currently refuses a fraction of generated traces; see the README.
    sidechains
        Place CB as well. Only meaningful with ``all_atom=True``.
    min_length
        Shortest region worth rebuilding. Shorter ones keep their input coordinates.
    seed
        Seed for reproducibility. With a fixed seed the output is bit-identical.

    Returns
    -------
    RebuildReport
        The models, the region assignments, and a per-region outcome for each. Check
        :attr:`RebuildReport.ok` or :attr:`RebuildReport.failures`.

    Examples
    --------
    >>> report = rebuild("AF-P04637-F1-model_v6.pdb", mode="expanded")  # doctest: +SKIP
    >>> print(report.summary())                                        # doctest: +SKIP
    >>> write_pdb(report.models, "p53_dodo.pdb")                       # doctest: +SKIP
    """
    from ..io import read_structure

    if n_models < 1:
        raise ValueError(f"n_models must be at least 1, got {n_models}.")

    original = read_structure(source) if not isinstance(source, Structure) else source
    report = RebuildReport()
    report.notes.extend(original.notes)

    rng = np.random.default_rng(seed)

    # Draw every model's end-to-end target for every region BEFORE building anything.
    #
    # A predicted end-to-end distance is the MEAN of a distribution, so an ensemble of models
    # must scatter around it -- otherwise the "pseudo-trajectory" is one conformation repeated,
    # which is what v1 shipped. The draws have to happen in a single call per region because
    # the sampler deliberately returns the mean exactly when asked for one conformer (a
    # single-model run should hit the target, not a random draw from around it).
    #
    # Only regions with a free end are scattered. An interior region with both anchors fixed
    # has its end-to-end distance *determined* by the anchor separation; there is nothing to
    # sample, and its models differ in path rather than in span.
    model_targets: dict[tuple[str, int, int], np.ndarray] = {}
    if n_models > 1:
        from ..engines.walk import sample_end_to_end_targets

        probe = original.copy()
        for assignment in assign_regions(probe, strategy=strategy):
            for domain in assignment.domains:
                if domain.kind is not DomainKind.IDR or not domain.span.is_terminal:
                    continue
                sequence = domain.sequence
                if len(sequence) < min_length:
                    continue
                mean = target_dimensions(sequence, mode=mode).end_to_end
                ceiling = 0.95 * (len(sequence) - 1) * CA_CA_BOND_LENGTH
                model_targets[(assignment.chain_id, domain.span.start, domain.span.stop)] = (
                    sample_end_to_end_targets(
                        mean,
                        n_models,
                        rng,
                        low=max(CA_CA_BOND_LENGTH, 0.15 * mean),
                        high=min(ceiling, 3.0 * mean),
                    )
                )

    # STEP 3, done ONCE for the whole run: move the folded domains so each linker's flanking
    # anchors sit at its predicted end-to-end distance.
    #
    # Without this the linker is built into whatever gap AlphaFold happened to leave, which is
    # essentially arbitrary. Measured on real models, AlphaFold packs domains 2-3.6x CLOSER
    # than the linkers between them predict (p300's 151-residue linker: a 26.1 A gap against a
    # 94.9 A prediction), because it has no way to know where to put domains joined by a long
    # disordered region. Building into that gap gives a compact blob wedged between domains.
    #
    # ONCE, not per model, and that is deliberate. The multi-model output exists so a viewer can
    # flick between conformers and see the disordered regions move; if the folded domains were
    # re-placed for every model they would jump around between frames, which looks wrong and
    # destroys the illusion. So the domain arrangement is shared and only the disordered regions
    # differ -- which also means a connecting IDR's end-to-end distance is fixed by its anchors
    # across all models (its path varies, its span cannot), while a terminal IDR is free and
    # does scatter.
    base = original.copy()
    base_assignments = assign_regions(base, strategy=strategy)
    report.assignments = base_assignments
    for assignment in base_assignments:
        report.notes.extend(assignment.notes)

    placement = reposition_folded_domains(base, mode=mode, rng=rng)
    report.placements.extend(placement.placements)
    report.notes.extend(placement.notes)
    for clashing in placement.clashing:
        report.notes.append(str(clashing))

    for model_number in range(1, n_models + 1):
        # Each model starts from the REPOSITIONED coordinates, so every model shares one domain
        # arrangement and a failed region in one cannot contaminate the next.
        # Regions are NOT re-assigned here. They were assigned once on the input geometry and
        # carried through by copy(), which rebinds every Domain view to the new structure.
        #
        # Re-assigning would be wrong, not merely wasteful: the density metric reads the
        # coordinates, and repositioning has just moved the folded domains apart, which changes
        # their contact density and shifts the boundaries. Measured on p300, re-assigning after
        # repositioning moved a folded domain's bounds from 569-650 to 570-644 -- so the anchors
        # that drove the placement would no longer be the anchors used to build against it.
        working = base.copy()

        with warnings.catch_warnings():
            # The pipeline supplies outer anchors, so the engine's warning about missing ones
            # should not fire. Promote it to an error so a regression here is loud rather than
            # buried in a user's terminal.
            warnings.simplefilter("error", category=_junction_warning())
            outcomes = _rebuild_one_model(
                working,
                model=model_number,
                mode=mode,
                engine_name=engine,
                rng=rng,
                min_length=min_length,
                sidechains=sidechains,
                all_atom=all_atom,
                model_targets=model_targets,
            )
        report.outcomes.extend(outcomes)
        # Rebuilt regions are CA-only, so drop the atoms that did not move with them.
        report.models.append(_drop_non_ca_from_rebuilt(working))

    return report


def _junction_warning() -> type[Warning]:
    """Return the engine's unconstrained-junction warning class, imported lazily."""
    from ..engines.walk import UnconstrainedJunctionWarning

    return UnconstrainedJunctionWarning


def build_from_sequence(
    sequence: str,
    *,
    mode: str = DEFAULT_MODE,
    n_models: int = 1,
    engine: str = "walk",
    all_atom: bool = False,
    sidechains: bool = False,
    seed: int | None = None,
) -> RebuildReport:
    """Build coordinates for a disordered sequence with no input structure.

    The sequence-only path: there are no folded domains and no anchors, so the whole chain is
    generated as one free region steered toward its predicted dimensions.

    Parameters
    ----------
    sequence
        One-letter amino acid sequence.
    mode, n_models, engine, all_atom, sidechains, seed
        As for :func:`rebuild`.

    Returns
    -------
    RebuildReport
        Models and outcomes, exactly as :func:`rebuild` returns.

    Examples
    --------
    >>> report = build_from_sequence("GSGSGSGSGSGSGSGSGSGS")  # doctest: +SKIP
    >>> write_pdb(report.models, "idr.pdb")                   # doctest: +SKIP
    """
    from ..constants import ONE_TO_THREE
    from ..engines.base import IDRRequest
    from ..engines.walk import SelfAvoidingWalk
    from ..structure import Span

    cleaned = sequence.strip().upper()
    if not cleaned.isalpha() or not cleaned:
        raise ValueError(
            f"sequence must be a non-empty string of one-letter amino acid codes, got {sequence!r}."
        )
    if n_models < 1:
        raise ValueError(f"n_models must be at least 1, got {n_models}.")
    if engine != "walk":
        raise ValueError(
            f"build_from_sequence currently supports engine='walk' only, got {engine!r}."
        )

    n = len(cleaned)
    rng = np.random.default_rng(seed)
    walk = SelfAvoidingWalk()
    report = RebuildReport()

    target = target_dimensions(cleaned, mode=mode)
    result = walk.generate(
        IDRRequest(
            sequence=cleaned,
            n_residues=n,
            target_end_to_end=target.end_to_end,
            n_anchor_xyz=None,
            c_anchor_xyz=None,
            n_conformations=n_models,
        ),
        obstacles=None,
        rng=rng,
    )

    for model_number in range(n_models):
        built = bool(result.success[model_number])
        coords = result.ca_coords[model_number]
        if built:
            structure = Structure.from_atom_records(
                xyz=coords,
                atom_name=["CA"] * n,
                element=["C"] * n,
                residue_name=[ONE_TO_THREE.get(c, "UNK") for c in cleaned],
                residue_number=list(range(1, n + 1)),
                chain_id=["A"] * n,
                source=f"sequence ({n} residues)",
            )
            structure.chains[0].domains = [
                Domain(structure=structure, span=Span(0, n), kind=DomainKind.IDR, rebuilt=True)
            ]
            if all_atom:
                place_backbone_for_domain(
                    structure, structure.chains[0].domains[0], rng=rng, sidechains=sidechains
                )
            report.models.append(structure)

        report.outcomes.append(
            RegionOutcome(
                model=model_number + 1,
                chain_id="A",
                residues=(1, n),
                n_residues=n,
                built=built,
                target=target if built else None,
                achieved_end_to_end=(
                    float(np.linalg.norm(coords[-1] - coords[0])) if built else None
                ),
                reason=None if built else "the engine reported failure for this conformer",
            )
        )

    return report
