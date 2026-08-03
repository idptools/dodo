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

import sys
import warnings
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..constants import (
    ANCHOR_ALWAYS_EXEMPT_ATOMS,
    ANCHOR_EXEMPT_ATOMS,
    ANCHOR_EXEMPT_ATOMS_BY_RESIDUE,
    CA_CA_BOND_LENGTH,
    DEFAULT_MODE,
    MIN_IDR_LENGTH,
    SHORT_REGION_TOLERANCE,
)
from ..exceptions import BuildError, DodoError, InvalidParameterError
from ..regions.identify import RegionAssignment, Strategy, assign_regions
from ..structure import Domain, DomainKind, Span, Structure
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

    @property
    def tolerated(self) -> bool:
        """True if this region was not rebuilt but that is acceptable, because it is short.

        A short disordered stretch left as AlphaFold drew it does not look wrong in a figure.
        It is the long regions that trail across the image as extended spaghetti, and those are
        what DODO exists to fix -- so under
        :data:`~dodo.constants.SHORT_REGION_TOLERANCE` residues a failure is reported and the run
        still succeeds. See :attr:`RebuildReport.ok`.
        """
        return not self.built and self.n_residues < SHORT_REGION_TOLERANCE

    def __str__(self) -> str:
        where = (
            f"model {self.model} chain {self.chain_id} "
            f"residues {self.residues[0]}-{self.residues[1]}"
        )
        if not self.built:
            label = "not rebuilt, left as-is" if self.tolerated else "NOT BUILT"
            return f"{where}: {label} ({self.reason})"
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
        """Every region that was not rebuilt. Their input coordinates were left in place.

        Includes the short ones that :attr:`ok` tolerates -- use :attr:`blocking_failures` for
        the subset that makes a run unsuccessful.
        """
        return [o for o in self.outcomes if not o.built]

    @property
    def blocking_failures(self) -> list[RegionOutcome]:
        """Regions long enough that failing to rebuild them is a real failure.

        See :attr:`RegionOutcome.tolerated` for where the line is drawn and why.
        """
        return [o for o in self.failures if not o.tolerated]

    @property
    def tolerated_failures(self) -> list[RegionOutcome]:
        """Regions not rebuilt, but short enough that it does not matter for a figure."""
        return [o for o in self.failures if o.tolerated]

    @property
    def n_built(self) -> int:
        """How many region-model pairs were successfully rebuilt."""
        return sum(1 for o in self.outcomes if o.built)

    @property
    def ok(self) -> bool:
        """True if every region that matters was rebuilt.

        A region shorter than :data:`~dodo.constants.SHORT_REGION_TOLERANCE` residues that could
        not be rebuilt does NOT make this False. Such a region is reported, and left with its
        input coordinates, but a few residues of AlphaFold geometry do not spoil a figure --
        whereas a long extended region does, which is the whole reason DODO exists.

        This is what the CLI's exit status reflects.
        """
        return not self.blocking_failures

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
        if self.blocking_failures:
            lines.append(f"  {len(self.blocking_failures)} failure(s):")
            lines += [f"    {f}" for f in self.blocking_failures]
        if self.tolerated_failures:
            lines.append(
                f"  {len(self.tolerated_failures)} short region(s) left as-is, under the "
                f"{SHORT_REGION_TOLERANCE}-residue threshold where this is not treated as a "
                f"failure:"
            )
            lines += [f"    {f}" for f in self.tolerated_failures]
        # Surface relaxed builds here too, not only on the individual outcome. A region built
        # against the relaxed anchor exemption may sit closer to its anchors' backbone than the
        # clash distance; that is a deliberate trade against leaving the region unbuilt, but it
        # has to be visible in the summary a CLI user actually reads.
        relaxed = [o for o in self.outcomes if o.built and o.reason]
        if relaxed:
            lines.append(
                f"  {len(relaxed)} region(s) needed a relaxed anchor exemption to be buildable:"
            )
            lines += [f"    {o}" for o in relaxed]
        lines += [f"  note: {n}" for n in self.notes]
        return "\n".join(lines)


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
    doomed = _doomed_atom_mask(structure)
    if not doomed.any():
        return structure
    return structure.select_atoms(~doomed)


def _doomed_atom_mask(structure: Structure) -> np.ndarray:
    """Atoms that :func:`_drop_non_ca_from_rebuilt` will delete: non-CA atoms of rebuilt regions.

    Factored out so that the two callers cannot drift apart, because they must agree. One
    deletes these atoms from the finished model; the other must not offer them to the engine as
    obstacles. When they disagreed, the engine spent its effort avoiding atoms that would not
    exist in the output -- 42-52% of the obstacle set on p300, measured per region.
    """
    doomed = np.zeros(structure.n_atoms, dtype=bool)
    not_ca = structure.atom_name != "CA"

    # Only regions DODO actually generated. A region whose build FAILED keeps its input
    # coordinates, so stripping it would destroy real side-chain data to no purpose and then
    # make the result look like DODO's work. Loops used to be stripped unconditionally, which
    # did exactly that.
    for domain in structure.domains:
        for span in domain.generated_spans():
            atoms = structure.atom_slice_for_residues(span.start, span.stop)
            doomed[atoms] = not_ca[atoms]
    return doomed


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

    request = IDRRequest(
        sequence=sequence,
        n_residues=len(sequence),
        target_end_to_end=requested,
        n_anchor_xyz=n_anchor_xyz,
        c_anchor_xyz=c_anchor_xyz,
        n_anchor_prev_xyz=_outer_ca(structure, span.n_anchor, step=-1),
        c_anchor_next_xyz=_outer_ca(structure, span.c_anchor, step=+1),
        n_conformations=1,
    )

    # TWO PASSES over the anchor exemption, strictest first.
    #
    # A region is bonded to its anchors' alpha carbons, so those must always be exempt from
    # clash checking -- without that there is no way to attach the region at all. The anchors'
    # BACKBONE atoms are a different matter. The residue actually bonded to an anchor does come
    # closer to them than the clash distance (measured over 649,658 sequence-neighbour pairs, an
    # alpha carbon sits 2.379 A from the next residue's N at the 0.1st percentile and 3.280 A at
    # the median), but a residue further along the region has no such licence -- and exempting
    # them region-wide let residues 3 to 16 positions away be placed against anchor backbone.
    #
    # So: try with only the alpha carbons exempt, which is the honest constraint, and fall back
    # to exempting the backbone only for a region that would otherwise not be built at all.
    # Measured over 36 structures, the strict pass alone takes contacts from 19 to 1 but takes
    # unbuilt regions from 2 to 8, and an unbuilt region is far more visible in a figure than a
    # contact a fraction of an Angstrom inside the limit. This gets both: strict where it is
    # free, relaxed only where the alternative is leaving the region as AlphaFold spaghetti.
    result = None
    failure: str | None = None
    relaxed = False
    for exempt_backbone in (False, True):
        try:
            attempt = engine.generate(  # type: ignore[attr-defined]
                request,
                obstacles=_obstacles_for_span(
                    structure, span, exempt_anchor_backbone=exempt_backbone
                ),
                rng=rng,
            )
        except DodoError as exc:
            failure = f"{label}: {type(exc).__name__}: {exc}"
            continue
        if bool(attempt.success[0]):
            result = attempt
            relaxed = exempt_backbone
            break
        failure = f"{label}: the engine reported failure for this conformer"

    if result is None:
        return outcome(built=False, reason=failure)

    structure.set_ca_xyz(np.arange(span.start, span.stop), result.ca_coords[0])
    if domain is not None:
        domain.rebuilt = True
        domain.placed = True

    achieved = float(np.linalg.norm(result.ca_coords[0][-1] - result.ca_coords[0][0]))
    return outcome(
        built=True,
        target=target,
        achieved_end_to_end=achieved,
        requested_end_to_end=requested,
        reason=(
            "built against a relaxed anchor exemption: the anchors' backbone atoms had to be "
            "exempted from clash checking for this region to be buildable at all, so it may "
            "sit closer to them than the clash distance"
            if relaxed
            else None
        ),
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


def _obstacles_for_span(
    structure: Structure, span: Span, *, exempt_anchor_backbone: bool = True
) -> np.ndarray | None:
    """Already-placed atoms the region must avoid.

    The region's own residues are excluded, since the engine handles self-avoidance internally.
    Its anchors are only *partly* excluded, and how far is the caller's choice -- see
    :func:`_build_region`, which tries the strict setting first.

    Parameters
    ----------
    structure
        The structure being built.
    span
        The region about to be rebuilt.
    exempt_anchor_backbone
        Also exempt the anchors' N, C and O (and proline's CD), not just their alpha carbons.

        The alpha carbons are *always* exempt: the region is bonded to them, and treating them
        as obstacles would leave no way to attach the region at all.

        The backbone is a judgement call. The residue actually bonded to an anchor does come
        closer to that backbone than the clash distance -- measured over 649,658
        sequence-neighbour pairs, an alpha carbon sits 2.379 A from the next residue's N at the
        0.1st percentile and 3.280 A at the median, so the median real junction is already inside
        the 3.20 A limit. But this exemption is not per-residue, so granting it also lets a
        residue 3 to 16 positions further along be placed against anchor backbone, which has no
        such justification. Hence: off for the first attempt, on only as a fallback.

    Notes
    -----
    Exempting the whole anchor *residue*, as this once did, is a defect rather than a
    conservative simplification: the region is bonded to the anchor's CA but has no bonded
    relationship to its *side chain*, so a blanket exemption let the walk place an alpha carbon
    straight through it, producing overlaps at 0.871-0.944 A -- below the shortest bond in any
    protein -- in output the pipeline reported as clean.
    """
    mask = structure.placed_atom_mask()
    mask[structure.atom_slice_for_residues(span.start, span.stop)] = False
    # Do not ask the engine to avoid atoms that will not survive into the output. Rebuilt
    # regions are reduced to alpha carbons by _drop_non_ca_from_rebuilt, so their N/C/O and side
    # chains are phantoms: avoiding them makes the build strictly harder than the real problem
    # without making the result any better. Measured on p300, they were 42-52% of the obstacle
    # set for every region.
    mask &= ~_doomed_atom_mask(structure)
    for anchor in (span.n_anchor, span.c_anchor):
        if anchor is None:
            continue
        atoms = structure.atom_slice_for_residues(anchor, anchor + 1)
        exempt = set(ANCHOR_ALWAYS_EXEMPT_ATOMS)
        if exempt_anchor_backbone:
            exempt |= ANCHOR_EXEMPT_ATOMS
            exempt |= ANCHOR_EXEMPT_ATOMS_BY_RESIDUE.get(
                str(structure.residue_name[anchor]), frozenset()
            )
        names = structure.atom_name[atoms]
        mask[atoms] = mask[atoms] & ~np.isin(names, list(exempt))
    if not mask.any():
        return None
    coords: np.ndarray = structure.xyz[mask]
    return coords


def _warn_starling_isolation() -> None:
    """Warn that STARLING never sees the folded domains. Short on purpose.

    The previous version of this ran to eleven lines. It was accurate and nobody read it, which
    makes it worse than a short one: the single fact a user has to take away is that the ensemble
    is not conditioned on the rest of the structure. The full explanation lives in the docs.

    Called from :func:`_make_engine`, which runs once per rebuild -- so a ten-model run warns once
    rather than ten times, and no module-level "have I said this yet" flag is needed to arrange it.
    """
    warnings.warn(
        "STARLING generates each region from its SEQUENCE ALONE -- it never sees the folded "
        "domains, so a conformer cannot avoid them and its statistics do not account for them. "
        "DODO picks the conformer that best fits the anchor separation and places it rigidly. "
        "Read a STARLING region as a realistic IDR conformation that has been positioned, not "
        "one sampled in the presence of the domains around it.",
        UserWarning,
        stacklevel=3,
    )


class _Progress:
    """A residue-weighted progress bar over the regions a rebuild will attempt.

    Wraps tqdm rather than exposing it, for two reasons. A rebuild is slow enough that silence
    reads as a hang -- a 2,400-residue structure spends most of a minute inside one region -- and
    the bar has to be genuinely optional: writing carriage returns into a piped log or a notebook
    cell is worse than showing nothing.
    """

    __slots__ = ("_bar",)

    def __init__(self, total: int) -> None:
        from tqdm.auto import tqdm

        self._bar = tqdm(
            total=max(total, 1),
            unit="res",
            desc="rebuilding",
            # Residues, not iterations: "1.2k/2.4k res" is meaningful where "3/6" is not.
            bar_format="{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} {unit} [{elapsed}]",
            leave=False,
        )

    def advance(self, residues: int) -> None:
        self._bar.update(residues)

    def next_model(self, done: int, total: int) -> None:
        if total > 1:
            # set_description_str, not set_description: the latter appends its own colon, and the
            # bar_format already has one, which reads as "rebuilding (model 2/3): :".
            self._bar.set_description_str(f"rebuilding (model {min(done + 1, total)}/{total})")

    def close(self) -> None:
        self._bar.close()


class _NoProgress:
    """The do-nothing tracker, so the pipeline never branches on whether a bar exists."""

    __slots__ = ()

    def advance(self, residues: int) -> None:
        return

    def next_model(self, done: int, total: int) -> None:
        return

    def close(self) -> None:
        return


def _progress_bar(requested: bool | None, total: int) -> _Progress | _NoProgress:
    """Build a progress tracker, or a no-op stand-in.

    ``None`` means decide from the environment: a bar on an interactive terminal, silence when
    stderr is redirected. tqdm is a hard dependency, but an ImportError is still tolerated rather
    than allowed to take down a rebuild -- failing a structure build over a progress bar would be
    absurd.
    """
    if requested is False:
        return _NoProgress()
    if requested is None and not sys.stderr.isatty():
        return _NoProgress()
    try:
        return _Progress(total)
    except ImportError:
        return _NoProgress()


def _make_engine(engine_name: str) -> object:
    """Build the conformation engine for a rebuild. Called ONCE, not once per model.

    The lifetime matters. This used to be constructed inside :func:`_rebuild_one_model`, so a
    multi-model run got a fresh engine per model -- and with it a fresh, empty ensemble cache.
    STARLING conditions on sequence alone, so every model then re-ran the diffusion and MDS for
    regions whose ensembles were already in hand: a 2-model, 3-region dnmt3a rebuild performed 6
    generations to obtain 3 distinct ensembles. Hoisting construction here is what lets
    :class:`~dodo.engines.starling.StarlingEngine` reuse them.
    """
    from ..engines.walk import SelfAvoidingWalk

    if engine_name == "walk":
        return SelfAvoidingWalk()
    if engine_name == "starling":
        from ..engines.hierarchical import HierarchicalEngine
        from ..engines.starling import StarlingEngine

        # Wrapped in HierarchicalEngine unconditionally, not only when a region turns out to be
        # too long. STARLING caps sequences at starling.configs.MAX_SEQUENCE_LENGTH (380), and
        # real IDRs routinely exceed it -- p300 has a 401-residue linker and a 583-residue tail.
        # Leaving the bare engine here meant those regions hard-failed with an error telling the
        # USER to wrap it, which is not something a user of `dodo rebuild --engine starling` can
        # act on. The wrapper splits an over-cap region into cap-sized segments and assembles
        # them; below the cap it delegates straight through, so wrapping costs nothing.
        engine = HierarchicalEngine(StarlingEngine())
        if not engine.available():
            raise BuildError(
                "The starling engine was requested but is not available. Install it with "
                "pip install 'idptools-dodo[starling]', or pass engine='walk'."
            )
        _warn_starling_isolation()
        return engine
    raise InvalidParameterError(f"Unknown engine {engine_name!r}. Use 'walk' or 'starling'.")


def _rebuild_one_model(
    structure: Structure,
    *,
    model: int,
    mode: str,
    engine: object,
    rng: np.random.Generator,
    min_length: int,
    model_targets: dict[tuple[str, int, int], np.ndarray],
    on_region_done: Callable[[int], None] | None = None,
) -> list[RegionOutcome]:
    """Rebuild every loop and IDR of one model, in place. Returns per-region outcomes.

    Takes a built ``engine`` rather than its name so that one instance spans every model; see
    :func:`_make_engine`.

    ``on_region_done`` is called with each region's residue count as it finishes, built or not,
    which is what drives the progress bar. Weighted by residues rather than counting regions
    because region lengths differ by two orders of magnitude -- on p300 they run from 10 to 583 --
    so a bar that advanced once per region would sit still through the only part that takes time.
    """
    outcomes: list[RegionOutcome] = []
    # Folded domains are already positioned (step 3 ran before this), so they are obstacles
    # from the outset. `placed`, not `rebuilt`: their atoms are moved rigidly, never generated.
    for domain in structure.folded_domains():
        domain.placed = True

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
    loops: list[tuple[Domain, int, Span]] = [
        (domain, index, loop)
        for domain in structure.folded_domains()
        for index, loop in enumerate(domain.loops)
    ]
    idrs = structure.idrs()
    connecting = [d for d in idrs if not d.span.is_terminal]
    terminal = [d for d in idrs if d.span.is_terminal]

    for parent, loop_index, loop in loops:
        loop_outcome = _build_region(
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
        outcomes.append(loop_outcome)
        if on_region_done is not None:
            on_region_done(len(loop))
        # Record success per loop. Without this, a loop that failed to build would still have
        # its side chains stripped and its inherited geometry attributed to DODO.
        if loop_outcome.built:
            parent.rebuilt_loops.add(loop_index)

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
        # Placed whatever happened. On success _build_region has already set it; on failure or
        # a skip the region keeps its input coordinates, and nothing will move them again -- so
        # it is final, and every region built after this one must avoid it. Leaving a failed
        # region out of the obstacle set is how AF-O14683-F1 ended up with a rebuilt alpha
        # carbon 1.27 A from an atom of the region that had just failed beside it.
        domain.placed = True
        if on_region_done is not None:
            on_region_done(len(domain.span))

    return outcomes


def rebuild(
    source: str | Path | Structure,
    *,
    mode: str = DEFAULT_MODE,
    n_models: int = 1,
    strategy: Strategy | str = Strategy.AUTO,
    engine: str = "walk",
    backbone: bool = False,
    min_length: int = MIN_IDR_LENGTH,
    seed: int | None = None,
    progress: bool | None = None,
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
        ``"walk"`` (always available) or ``"starling"``, which needs
        ``pip install 'idptools-dodo[starling]'``.
    backbone
        Place N, C and O on the rebuilt regions, from four consecutive alpha carbons, then refine.
        **Opt-in and off by default.** Folded domains are untouched either way: they keep every
        atom they arrived with, and only regions DODO generated gain a backbone.

        Held out against all-atom simulation: N 0.16 A, C 0.22 A, O 0.63 A, with every bond
        length exact by construction.

        What keeps it opt-in is the seams. Where a rebuilt region meets a folded domain, an exact
        peptide bond is not merely hard but geometrically impossible: a peptide unit reaches
        2.854 A from an alpha carbon to the nitrogen it bonds to, and a rebuilt alpha carbon sits
        3.2-4.5 A from that fixed nitrogen. DODO aims the carbon at it, leaving the bond long --
        about 2.2 A against an ideal 1.329 -- rather than writing two atoms on top of each other,
        and reports every seam it does this at.
    min_length
        Shortest region worth rebuilding. Shorter ones keep their input coordinates.
    seed
        Seed for reproducibility. With a fixed seed the output is bit-identical.
    progress
        Show a progress bar on stderr. ``None`` (the default) shows one when stderr is a terminal
        and stays silent otherwise, so piped output and log files do not fill with carriage
        returns. ``False`` disables it outright.

        Weighted by residues rather than by regions: region lengths span two orders of magnitude
        -- on p300 they run from 10 to 583 -- so a bar advancing once per region would appear
        stuck through the only part that takes real time.

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
        raise InvalidParameterError(f"n_models must be at least 1, got {n_models}.")

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

    # Total work for the progress bar: every region every model will attempt, weighted by length.
    # Counted from the assignment rather than guessed, so the bar cannot overrun or stop short.
    total_work = n_models * sum(
        len(span)
        for assignment in base_assignments
        for domain in assignment.domains
        for span in ([domain.span] if domain.kind is DomainKind.IDR else list(domain.loops))
    )
    tracker = _progress_bar(progress, total_work)
    # One engine for the whole run, so a per-sequence cache inside it survives across models.
    engine_instance = _make_engine(engine)

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
                engine=engine_instance,
                rng=rng,
                min_length=min_length,
                model_targets=model_targets,
                on_region_done=tracker.advance,
            )
        report.outcomes.extend(outcomes)
        # Rebuilt regions are CA-only, so drop the atoms that did not move with them. This has to
        # happen BEFORE any backbone placement: it deletes every non-CA atom in a generated region,
        # which would otherwise take the N, C and O placed just below straight back out again.
        final = _drop_non_ca_from_rebuilt(working)
        if backbone:
            from .ca_backbone import add_backbone_to_rebuilt

            final = add_backbone_to_rebuilt(final)
        report.models.append(final)
        tracker.next_model(model_number, n_models)

    tracker.close()
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
    backbone: bool = False,
    seed: int | None = None,
) -> RebuildReport:
    """Build coordinates for a disordered sequence with no input structure.

    The sequence-only path: there are no folded domains and no anchors, so the whole chain is
    generated as one free region steered toward its predicted dimensions.

    Parameters
    ----------
    sequence
        One-letter amino acid sequence.
    mode, n_models, engine, backbone, seed
        As for :func:`rebuild`. Note that ``backbone`` is cleaner here than on
        :func:`rebuild`: with no input structure there are no folded domains and so no seams,
        which is where the strained bonds in a rebuild come from.

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
        raise InvalidParameterError(
            f"sequence must be a non-empty string of one-letter amino acid codes, got {sequence!r}."
        )
    if n_models < 1:
        raise InvalidParameterError(f"n_models must be at least 1, got {n_models}.")
    if engine != "walk":
        raise InvalidParameterError(
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
            if backbone:
                from .ca_backbone import add_backbone_to_rebuilt

                structure = add_backbone_to_rebuilt(structure)
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
