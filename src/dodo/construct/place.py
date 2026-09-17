"""Reposition folded domains so linker IDRs can adopt their predicted dimensions.

This is step 3 of DODO's algorithm, and it is the step that makes the whole thing work.

AlphaFold places folded domains at whatever separation its prediction happened to produce,
which for a long disordered linker is essentially arbitrary. If a linker's predicted end-to-end
distance is 100 A and AlphaFold left its flanking domains 300 A apart, then rebuilding the
linker *in place* can only produce a stretched rod -- the chain is forced to span a distance its
sequence says it should not. The fix is not a better linker builder. It is to **move the
domains** so their separation matches the prediction, and only then build the linker.

What moves, and what moves together
-----------------------------------
Folded-domain atoms are **never rebuilt**. They come from AlphaFold (or a crystal structure) and
are trusted. A domain only ever moves as part of a **rigid unit** -- see
:mod:`dodo.construct.assembly` -- and a unit moves as one rigid body or not at all.

That indirection is the whole complex story. On a single chain every folded domain is its own
unit and this module does exactly what it always did. In a complex the domains that pack against
each other are one unit, so satisfying a linker prediction can never take an interface apart:
there is no code path that moves one domain of a unit without the others, because every
transform here goes through :meth:`~dodo.construct.assembly.RigidUnit.rotate`, which requires a
common centre of rotation.

Which units can move
--------------------
Only a chain is covalent, so the only thing that relates two units geometrically is a connecting
IDR inside one chain. Units and those linkers form a graph; each connected component is placed
independently, anchored on its earliest unit, and units in different components keep whatever
relative position the input gave them, because nothing in the input says what else it should be.

A consequence worth stating plainly: in a complex where every folded domain packs into one
assembly there is a single unit, no linker crosses between units, and **this module correctly
does nothing**.

Placing one unit
----------------
A unit reached by a single linker is placed the way DODO has always placed a domain: sample a
direction, put the attachment alpha carbon at the predicted separation along it, turn the unit so
its body extends away from the domain it connects back to, perturb, and reject on clash.

A unit reached by two or more linkers cannot be placed by choosing one direction, and in a
complex that is common rather than exotic. Those are placed by projecting each attachment point
onto the sphere its constraint defines and superposing onto the result, iterated -- which
degenerates to a pure translation for one constraint, hits the exact answer for two when one
exists, and returns the least-squares compromise when the constraints cannot all be met. The
same step, swept over a whole component, is what closes a cycle.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Callable
from dataclasses import dataclass, field, replace
from itertools import pairwise

import numpy as np
from scipy.spatial import cKDTree

from ..constants import (
    CA_CLASH_DISTANCE,
    DEFAULT_MODE,
    MAX_FD_PLACEMENT_ATTEMPTS,
    MAX_UNIT_RELAX_ITERATIONS,
    MAX_UNIT_RELAX_SWEEPS,
    MIN_IDR_LENGTH,
    UNIT_PLACEMENT_TOLERANCE,
)
from ..engines.walk import max_reach, min_reach
from ..exceptions import BuildError, GeometryError
from ..geometry.transforms import rotation_between_vectors, superpose
from ..structure import Chain, Domain, DomainKind, Structure
from .assembly import Assembly, Interface, RigidUnit, _rigid_drift, find_rigid_units
from .dimensions import DimensionTarget, target_dimensions

__all__ = [
    "DomainPlacement",
    "LinkerOutcome",
    "PlacementReport",
    "UnitPlacement",
    "reposition_folded_domains",
    "verify_rigid",
]

#: Maximum random perturbation applied to the "face the neighbour" orientation, in degrees.
#:
#: CHOICE. Large enough that repeated runs explore genuinely different arrangements, small
#: enough that the perturbation does not undo the facing constraint it perturbs. v1 used fully
#: random directions and the author noted that worked better than being clever about the
#: system's centre -- this keeps most of that freedom while still preferring sane linker routes.
PERTURBATION_DEGREES = 60.0


@dataclass(frozen=True, slots=True)
class DomainPlacement:
    """Where one folded domain ended up, and against what target."""

    chain_id: str
    #: 1-based inclusive residue range, as a user reads it off a PDB file.
    residues: tuple[int, int]
    moved: bool
    #: Target anchor separation in Angstroms, from the linker's predicted dimensions. ``None``
    #: when this domain carries no linker of its own -- which happens when it moved as part of a
    #: multi-domain rigid unit that some *other* domain's linker positioned.
    target_separation: float | None = None
    achieved_separation: float | None = None
    attempts: int = 0
    clashing: bool = False
    reason: str | None = None
    #: Index of the rigid unit this domain belongs to.
    unit: int = 0

    def __str__(self) -> str:
        where = f"chain {self.chain_id} FD {self.residues[0]}-{self.residues[1]}"
        if not self.moved:
            return f"{where}: held in place ({self.reason})"
        detail = ""
        if self.target_separation is not None and self.achieved_separation is not None:
            detail = (
                f", anchor separation {self.achieved_separation:.1f} A "
                f"against a {self.target_separation:.1f} A target"
            )
        elif self.reason:
            detail = f", {self.reason}"
        flag = " [CLASHING]" if self.clashing else ""
        return f"{where}: moved in {self.attempts} attempt(s){detail}{flag}"


@dataclass(frozen=True, slots=True)
class UnitPlacement:
    """What happened to one rigid unit."""

    index: int
    chains: tuple[str, ...]
    domains: tuple[str, ...]
    n_atoms: int
    moved: bool
    reason: str | None = None
    attempts: int = 0
    clashing: bool = False
    clashing_atoms: int = 0
    #: Worst ``|achieved - target|`` over the linkers that positioned this unit, in Angstroms.
    worst_residual: float | None = None
    n_constraints: int = 0

    def __str__(self) -> str:
        what = (
            f"unit {self.index} ({len(self.domains)} domain(s) over "
            f"chain(s) {', '.join(self.chains)}, {self.n_atoms} atoms)"
        )
        if not self.moved:
            return f"{what}: held in place ({self.reason})"
        detail = f" against {self.n_constraints} linker constraint(s)"
        if self.worst_residual is not None:
            detail += f", worst error {self.worst_residual:.2f} A"
        flag = f" [CLASHING: {self.clashing_atoms} atoms]" if self.clashing else ""
        return f"{what}: moved in {self.attempts} attempt(s){detail}{flag}"


@dataclass(frozen=True, slots=True)
class LinkerOutcome:
    """One connecting IDR's span: what was predicted, and what the geometry actually allows."""

    chain_id: str
    #: 1-based inclusive residue range of the linker itself.
    residues: tuple[int, int]
    n_residues: int
    target: float
    achieved: float
    #: True when both flanking domains are in the same rigid unit, so the span was not ours to
    #: set. The prediction is reported anyway: a large disagreement is real information about the
    #: input, and hiding it would misrepresent what DODO did.
    dictated: bool
    #: True when the span exceeds what this many residues can physically bridge, so the region
    #: cannot be built at all.
    unbridgeable: bool
    ceiling: float

    def __str__(self) -> str:
        where = f"chain {self.chain_id} linker {self.residues[0]}-{self.residues[1]}"
        if self.unbridgeable:
            return (
                f"{where}: UNBRIDGEABLE -- its flanking domains are {self.achieved:.1f} A apart "
                f"but {self.n_residues} residue(s) span at most {self.ceiling:.1f} A"
            )
        if self.dictated:
            return (
                f"{where}: span {self.achieved:.1f} A, dictated by the complex "
                f"(both flanking domains are in one rigid unit); predicted {self.target:.1f} A"
            )
        return f"{where}: span {self.achieved:.1f} A against a {self.target:.1f} A target"


@dataclass(slots=True)
class PlacementReport:
    """Outcome of repositioning every folded domain in a structure."""

    placements: list[DomainPlacement] = field(default_factory=list)
    units: list[UnitPlacement] = field(default_factory=list)
    linkers: list[LinkerOutcome] = field(default_factory=list)
    interfaces: list[Interface] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def moved(self) -> list[DomainPlacement]:
        """Domains that were repositioned."""
        return [p for p in self.placements if p.moved]

    @property
    def clashing(self) -> list[DomainPlacement]:
        """Domains that could not be placed without a clash."""
        return [p for p in self.placements if p.clashing]

    @property
    def rigid_units(self) -> list[UnitPlacement]:
        """Units holding more than one folded domain: the ones locking a complex together."""
        return [u for u in self.units if len(u.domains) > 1]

    @property
    def unbridgeable_linkers(self) -> list[LinkerOutcome]:
        """Linkers whose flanking domains are too far apart to connect."""
        return [linker for linker in self.linkers if linker.unbridgeable]

    @property
    def broken_interfaces(self) -> list[Interface]:
        """Contacts between folded domains that were not locked, and so may have been broken."""
        return [i for i in self.interfaces if not i.preserved]

    @property
    def ok(self) -> bool:
        """True if every moved unit was placed without a clash and every linker can be built."""
        return not self.clashing and not self.unbridgeable_linkers

    def summary(self) -> str:
        """Multi-line human-readable summary."""
        lines = [f"{len(self.moved)}/{len(self.placements)} folded domain(s) repositioned"]
        if self.rigid_units:
            lines.append(
                f"  {len(self.rigid_units)} rigid unit(s) hold more than one folded domain "
                f"and were moved, if at all, as one body:"
            )
            lines += [f"    {u}" for u in self.rigid_units]
        lines += [f"  {p}" for p in self.placements]
        if self.unbridgeable_linkers:
            lines.append(f"  {len(self.unbridgeable_linkers)} unbridgeable linker(s):")
            lines += [f"    {linker}" for linker in self.unbridgeable_linkers]
        lines += [f"  note: {n}" for n in self.notes]
        return "\n".join(lines)


def verify_rigid(before: np.ndarray, after: np.ndarray, *, tolerance: float = 1e-6) -> None:
    """Assert that a set of atoms survived a transform unchanged.

    Folded-domain atoms are never rebuilt, so any transform applied to one must be a rigid
    motion: every internal pairwise distance is preserved. This checks that directly rather than
    trusting the arithmetic, because a subtly non-orthonormal rotation would deform a real
    domain while still looking plausible in a viewer.

    Fits the best proper rigid transform and checks every aligned atom's residual. This is O(n)
    and, unlike centroid radii, cannot accept a non-rigid rearrangement that happens to leave
    every atom the same distance from the centroid.

    **Called on a whole rigid unit, not on one domain.** Rotating each domain of a unit about
    its own centroid preserves every domain's internal geometry perfectly and still takes the
    assembly apart, so a per-domain check would pass on exactly the bug that matters.

    Raises
    ------
    GeometryError
        If the transform was not rigid.
    """
    if before.shape != after.shape:
        raise GeometryError(
            f"Cannot compare geometry: shapes differ, {before.shape} vs {after.shape}."
        )
    if before.shape[0] < 2:
        return
    drift = _rigid_drift(before, after)
    if drift > tolerance:
        raise GeometryError(
            f"A folded domain was deformed rather than moved rigidly: the worst atom's "
            f"distance to the centroid changed by {drift:.3e} A. Folded-domain atoms "
            f"must never be rebuilt or rescaled, only rotated and translated."
        )


def _exit_ca(structure: Structure, domain: Domain) -> np.ndarray:
    """CA of the domain's last residue: where a C-terminal linker departs from."""
    coords: np.ndarray = structure.ca_xyz[domain.span.stop - 1]
    return coords


def _entry_ca(structure: Structure, domain: Domain) -> np.ndarray:
    """CA of the domain's first residue: where an N-terminal linker arrives at."""
    coords: np.ndarray = structure.ca_xyz[domain.span.start]
    return coords


def _linker_between(chain: Chain, first: Domain, second: Domain) -> Domain | None:
    """Return the IDR between two folded domains, if they are consecutive with one between."""
    ordered = sorted(chain.domains, key=lambda d: d.span.start)
    starts = [d.span.start for d in ordered]
    try:
        index = starts.index(first.span.start)
    except ValueError:  # pragma: no cover - callers pass domains from this chain
        return None
    if index + 2 >= len(ordered) or ordered[index + 2].span.start != second.span.start:
        return None
    middle = ordered[index + 1]
    return middle if middle.kind is DomainKind.IDR else None


# ----------------------------------------------------------------------------------
# Constraints: one per connecting IDR
# ----------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class _Constraint:
    """A connecting IDR, as a distance between one atom of each of two rigid units.

    ``lo`` is the domain the linker leaves (so the attachment point is its *exit* alpha carbon)
    and ``hi`` the one it arrives at (its *entry* alpha carbon). Which of the two gets moved
    depends on the traversal order, so both directions have to be expressible.
    """

    chain_id: str
    residues: tuple[int, int]
    n_residues: int
    lo: Domain
    hi: Domain
    lo_unit: int
    hi_unit: int
    target: float
    ceiling: float
    floor: float

    @property
    def inter_unit(self) -> bool:
        return self.lo_unit != self.hi_unit

    def point_lo(self, structure: Structure) -> np.ndarray:
        return _exit_ca(structure, self.lo)

    def point_hi(self, structure: Structure) -> np.ndarray:
        return _entry_ca(structure, self.hi)

    def separation(self, structure: Structure) -> float:
        return float(np.linalg.norm(self.point_hi(structure) - self.point_lo(structure)))

    def other_unit(self, unit_index: int) -> int:
        return self.hi_unit if unit_index == self.lo_unit else self.lo_unit

    def attachment(self, structure: Structure, unit_index: int) -> np.ndarray:
        """Where the linker attaches to the unit being moved."""
        return (
            self.point_lo(structure) if unit_index == self.lo_unit else self.point_hi(structure)
        )

    def anchor(self, structure: Structure, unit_index: int) -> np.ndarray:
        """Where the linker attaches to the *other* unit, which is not being moved."""
        return (
            self.point_hi(structure) if unit_index == self.lo_unit else self.point_lo(structure)
        )

    def attached_domain(self, unit_index: int) -> Domain:
        return self.lo if unit_index == self.lo_unit else self.hi


def _constraints(
    structure: Structure,
    assembly: Assembly,
    *,
    mode: str,
    on_linker_done: Callable[[int], None] | None = None,
) -> list[_Constraint]:
    """One constraint per connecting IDR between two folded domains, in structure order.

    Each one costs a dimension prediction, and on a large assembly that is the slowest thing
    this module does -- measured on the 808-chain nuclear pore, 68 s, against under a second for
    the placement itself. Hence the callback: without it the stage shows "0/1 units" and sits
    there, which is the frozen line a progress indicator exists to prevent.
    """
    constraints: list[_Constraint] = []
    for chain in structure.chains:
        folded = [
            d
            for d in sorted(chain.domains, key=lambda d: d.span.start)
            if d.kind is DomainKind.FOLDED
        ]
        for previous, current in pairwise(folded):
            linker = _linker_between(chain, previous, current)
            if linker is None:
                # No IDR between them: either they are covalently adjacent or the tiling put
                # something else there. Either way there is no linker to satisfy, and
                # find_rigid_units has already put covalently adjacent domains in one unit.
                continue
            n = len(linker.span)
            target: DimensionTarget = target_dimensions(linker.sequence, mode=mode)
            ceiling = max_reach(n + 1)
            constraints.append(
                _Constraint(
                    chain_id=chain.chain_id,
                    residues=(
                        int(structure.residue_number[linker.span.start]),
                        int(structure.residue_number[linker.span.stop - 1]),
                    ),
                    n_residues=n,
                    lo=previous,
                    hi=current,
                    lo_unit=assembly.unit_of(previous).index,
                    hi_unit=assembly.unit_of(current).index,
                    # A target above the reach of this many bonds is unreachable at any
                    # conformation, so aiming at it would place the domains where the linker
                    # cannot follow. dimensions clamps to 0.95 * contour_length, which is a
                    # slightly different (and looser) bound than the pseudo-angle reach.
                    target=min(target.end_to_end, ceiling),
                    ceiling=ceiling,
                    floor=min_reach(n + 1),
                )
            )
            if on_linker_done is not None:
                on_linker_done(1)
    return constraints


# ----------------------------------------------------------------------------------
# Moving a unit
# ----------------------------------------------------------------------------------


def _snapshot(unit: RigidUnit) -> list[np.ndarray]:
    return [domain.xyz.copy() for domain in unit.domains]


def _restore(unit: RigidUnit, snapshot: list[np.ndarray]) -> None:
    for domain, coords in zip(unit.domains, snapshot, strict=True):
        domain.xyz[:] = coords


def _apply_transform(unit: RigidUnit, rotation: np.ndarray, translation: np.ndarray) -> None:
    """Apply ``coords @ R.T + t`` to a whole unit, rotating about its centroid.

    Mathematically identical to rotating about the origin and translating, but done about the
    centroid so the arithmetic stays conditioned: an assembly can sit a thousand Angstroms from
    the origin, and rotating there throws away significant digits for no reason.
    """
    centre = unit.centroid()
    unit.rotate(rotation, about=centre)
    unit.translate(translation + rotation @ centre - centre)


def _orient_away_from(
    structure: Structure,
    unit: RigidUnit,
    *,
    attachment: np.ndarray,
    from_point: np.ndarray,
    rng: np.random.Generator,
    perturb: bool,
) -> None:
    """Rotate ``unit`` in place so its body extends away from ``from_point``.

    The construction: the vector from the attachment alpha carbon to the unit's centroid should
    point in the same direction as the vector from ``from_point`` to that alpha carbon. Then the
    attachment point is the closest part of the unit to the one it connects back to, and the
    linker arrives at the near face rather than having to wrap around to a far one.

    Rotation is about the attachment point, so it stays put and a subsequent translation can
    place it exactly.
    """
    centroid = unit.centroid()
    body = centroid - attachment
    if float(np.linalg.norm(body)) < 1e-9:
        # A single-residue unit has no body direction to align. Nothing to orient.
        return
    approach = attachment - np.asarray(from_point, dtype=np.float64)
    if float(np.linalg.norm(approach)) < 1e-9:
        # The attachment sits on the other unit's own attachment point; any orientation is as
        # good as another.
        return

    unit.rotate(rotation_between_vectors(body, approach), about=attachment)

    if perturb:
        # A small random rotation about the attachment point, so repeated runs differ without
        # abandoning the facing constraint just established.
        axis = rng.normal(size=3)
        axis /= np.linalg.norm(axis)
        angle = np.deg2rad(rng.uniform(-PERTURBATION_DEGREES, PERTURBATION_DEGREES))
        from ..geometry.transforms import rotation_from_axis_angle

        unit.rotate(rotation_from_axis_angle(axis, float(angle)), about=attachment)


def _place_at_separation(
    unit: RigidUnit,
    attachment: np.ndarray,
    *,
    from_point: np.ndarray,
    separation: float,
    direction: np.ndarray,
) -> None:
    """Translate ``unit`` so its attachment CA sits ``separation`` A from ``from_point``."""
    unit_vector = np.asarray(direction, dtype=np.float64)
    norm = float(np.linalg.norm(unit_vector))
    if norm < 1e-9:
        raise GeometryError("Cannot place a unit along a zero-length direction.")
    wanted = np.asarray(from_point, dtype=np.float64) + (unit_vector / norm) * separation
    unit.translate(wanted - attachment)


def _clash_count(unit: RigidUnit, obstacles: cKDTree | None, cutoff: float) -> int:
    """Atoms of ``unit`` within ``cutoff`` of an already-placed atom.

    Queries a tree built once per placement rather than going through
    :meth:`Structure.clash_mask`, which rebuilds its index on every call. Measured on a
    61,511-atom assembly with 50,000 obstacle atoms and a 3,000-atom query: 7.6 ms per call
    against 0.69 ms with the tree already built, so at the 500-attempt budget this is the
    difference between 3.8 s and 0.35 s for a single unit. Identical semantics -- both count a
    query atom as clashing when any obstacle lies at or within the cutoff.
    """
    if obstacles is None:
        return 0
    lengths = obstacles.query_ball_point(unit.xyz, cutoff, return_length=True)
    return int(np.count_nonzero(lengths))


def _place_one_constraint(
    structure: Structure,
    unit: RigidUnit,
    constraint: _Constraint,
    *,
    rng: np.random.Generator,
    perturb: bool,
    obstacles: cKDTree | None,
    max_attempts: int,
    clash_distance: float,
) -> tuple[int, int]:
    """Place a unit held by exactly one linker. Returns ``(attempts, clashing_atoms)``.

    This is DODO's original folded-domain placement, generalised from a domain to a unit, and it
    is kept in that exact shape on purpose: for a single chain every unit holds one domain, so
    this path draws the same random numbers in the same order as before rigid units existed and
    the output is byte-identical. The regression suite asserts that.
    """
    from_point = constraint.anchor(structure, unit.index)
    separation = constraint.target
    before = _snapshot(unit)

    best_direction: np.ndarray | None = None
    best_clash_count = np.inf
    attempts = 0

    for attempt in range(1, max_attempts + 1):
        attempts = attempt
        # Restore, then re-place: each attempt is independent, so a rejected orientation cannot
        # bias the next one.
        _restore(unit, before)

        direction = rng.normal(size=3)
        if float(np.linalg.norm(direction)) < 1e-9:
            continue
        attachment = constraint.attachment(structure, unit.index)
        _place_at_separation(
            unit, attachment, from_point=from_point, separation=separation, direction=direction
        )
        _orient_away_from(
            structure,
            unit,
            attachment=constraint.attachment(structure, unit.index),
            from_point=from_point,
            rng=rng,
            perturb=perturb,
        )
        # Orientation rotates about the attachment CA, which the translation had already
        # positioned, so the separation still holds. Re-place anyway: the rotation is about a
        # point that itself came from the pre-rotation geometry, and re-placing costs nothing
        # while removing any doubt.
        attachment = constraint.attachment(structure, unit.index)
        _place_at_separation(
            unit,
            attachment,
            from_point=from_point,
            separation=separation,
            direction=attachment - from_point,
        )

        clash_count = _clash_count(unit, obstacles, clash_distance)
        if clash_count == 0:
            return attempts, 0
        if clash_count < best_clash_count:
            best_clash_count = clash_count
            best_direction = direction

    # Exhausted attempts. Re-apply the least-bad arrangement rather than leaving whatever the
    # final attempt happened to produce, and report the clash instead of hiding it.
    if best_direction is not None:
        _restore(unit, before)
        attachment = constraint.attachment(structure, unit.index)
        _place_at_separation(
            unit,
            attachment,
            from_point=from_point,
            separation=separation,
            direction=best_direction,
        )
        _orient_away_from(
            structure,
            unit,
            attachment=constraint.attachment(structure, unit.index),
            from_point=from_point,
            rng=rng,
            perturb=False,
        )
        attachment = constraint.attachment(structure, unit.index)
        _place_at_separation(
            unit,
            attachment,
            from_point=from_point,
            separation=separation,
            direction=attachment - from_point,
        )
    return attempts, int(best_clash_count) if best_direction is not None else 0


def _relax_unit(
    structure: Structure,
    unit: RigidUnit,
    constraints: list[_Constraint],
    *,
    tolerance: float = UNIT_PLACEMENT_TOLERANCE,
    max_iterations: int = MAX_UNIT_RELAX_ITERATIONS,
) -> float:
    """Move one unit to satisfy several distance constraints at once. Returns the worst error.

    Alternating projection. Each constraint says the attachment point must lie on a sphere; the
    projection step moves each point to the nearest point of its own sphere, and the superposition
    step finds the single rigid motion that comes closest to hitting all of those at once. Repeat.

    Why this rather than an optimiser: it needs no gradients or step sizes, it is rigid by
    construction (the only thing ever applied is a proper rotation plus a translation), and it
    degenerates correctly at every constraint count -- a pure translation onto the sphere for one,
    the exact answer for two whenever one exists, and the least-squares compromise for three or
    more, which is the honest answer to an over-determined system.
    """
    if not constraints:
        return 0.0
    anchors = np.empty((len(constraints), 3), dtype=np.float64)
    targets = np.empty(len(constraints), dtype=np.float64)
    for k, constraint in enumerate(constraints):
        anchors[k] = constraint.anchor(structure, unit.index)
        targets[k] = constraint.target

    worst = np.inf
    for _ in range(max_iterations):
        points = np.array(
            [c.attachment(structure, unit.index) for c in constraints], dtype=np.float64
        )
        offsets = points - anchors
        norms = np.linalg.norm(offsets, axis=1)
        worst = float(np.abs(norms - targets).max())
        if worst <= tolerance:
            return worst
        # An attachment sitting exactly on its anchor gives no direction to project along. Use
        # the direction from the anchor to the unit's centroid, which is always defined for a
        # unit with more than one atom and is the sensible way to push the unit out.
        fallback = unit.centroid() - anchors
        fallback_norms = np.linalg.norm(fallback, axis=1)
        directions = np.where(
            norms[:, None] > 1e-9,
            offsets / np.maximum(norms, 1e-12)[:, None],
            np.where(
                fallback_norms[:, None] > 1e-9,
                fallback / np.maximum(fallback_norms, 1e-12)[:, None],
                np.array([1.0, 0.0, 0.0]),
            ),
        )
        wanted = anchors + directions * targets[:, None]
        rotation, translation = superpose(points, wanted)
        _apply_transform(unit, rotation, translation)
    return worst


def _place_multi_constraint(
    structure: Structure,
    unit: RigidUnit,
    constraints: list[_Constraint],
    *,
    rng: np.random.Generator,
    perturb: bool,
    obstacles: cKDTree | None,
    max_attempts: int,
    clash_distance: float,
    tolerance: float,
) -> tuple[int, int, float]:
    """Place a unit held by two or more linkers. Returns ``(attempts, clashing_atoms, residual)``.

    Each attempt seeds the relaxation from a fresh single-constraint placement, so the restarts
    explore genuinely different basins rather than re-converging on the one the first draw found.
    The budget is spent on restarts because the relaxation itself is deterministic: once started,
    it goes where the constraints send it.
    """
    before = _snapshot(unit)
    best: tuple[float, float, list[np.ndarray]] | None = None
    attempts = 0
    # Fewer restarts than the single-constraint budget: each one runs a full relaxation, and a
    # unit pinned by two or more linkers has far less freedom to explore in the first place.
    restarts = max(1, max_attempts // 25)

    for attempt in range(1, restarts + 1):
        attempts = attempt
        _restore(unit, before)
        _place_one_constraint(
            structure,
            unit,
            constraints[0],
            rng=rng,
            perturb=perturb,
            obstacles=None,  # seeding only; the placement is about to move again
            max_attempts=1,
            clash_distance=clash_distance,
        )
        residual = _relax_unit(structure, unit, constraints, tolerance=tolerance)
        clashes = _clash_count(unit, obstacles, clash_distance)
        if clashes == 0 and residual <= tolerance:
            return attempts, 0, residual
        # Rank on clashes first: a placement that overlaps another chain is wrong in a way that
        # a linker half an Angstrom off its predicted mean is not.
        score = (float(clashes), residual)
        if best is None or score < (best[0], best[1]):
            best = (float(clashes), residual, _snapshot(unit))

    if best is not None:
        _restore(unit, best[2])
        return attempts, int(best[0]), best[1]
    return attempts, 0, 0.0


# ----------------------------------------------------------------------------------
# The traversal
# ----------------------------------------------------------------------------------


def _components(n_units: int, constraints: list[_Constraint]) -> list[list[int]]:
    """Group the units into connected components, each sorted, ordered by lowest unit."""
    adjacency: dict[int, set[int]] = {k: set() for k in range(n_units)}
    for constraint in constraints:
        if constraint.inter_unit:
            adjacency[constraint.lo_unit].add(constraint.hi_unit)
            adjacency[constraint.hi_unit].add(constraint.lo_unit)

    seen: set[int] = set()
    components: list[list[int]] = []
    for start in range(n_units):
        if start in seen:
            continue
        stack, group = [start], []
        seen.add(start)
        while stack:
            node = stack.pop()
            group.append(node)
            for neighbour in sorted(adjacency[node]):
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
        components.append(sorted(group))
    return components


def _obstacle_tree(structure: Structure, atom_mask: np.ndarray) -> cKDTree | None:
    if not atom_mask.any():
        return None
    return cKDTree(structure.xyz[atom_mask])


def _unit_clash_count(
    unit: RigidUnit,
    assembly: Assembly,
    cutoff: float,
    *,
    obstacle_mask: np.ndarray | None = None,
) -> int:
    """Count atoms of ``unit`` clashing with atoms in every other folded rigid unit."""
    if obstacle_mask is None:
        obstacles = np.zeros(unit.structure.n_atoms, dtype=bool)
        for other in assembly.units:
            obstacles |= other.atom_mask
    else:
        obstacles = obstacle_mask.copy()
    obstacles[unit.atom_mask] = False
    return _clash_count(unit, _obstacle_tree(unit.structure, obstacles), cutoff)


def reposition_folded_domains(
    structure: Structure,
    *,
    mode: str = DEFAULT_MODE,
    rng: np.random.Generator,
    max_attempts: int = MAX_FD_PLACEMENT_ATTEMPTS,
    clash_distance: float = CA_CLASH_DISTANCE,
    perturb: bool = True,
    assembly: Assembly | None = None,
    min_length: int = MIN_IDR_LENGTH,
    lock_interfaces: bool = True,
    lock_intra_chain_interfaces: bool = False,
    tolerance: float = UNIT_PLACEMENT_TOLERANCE,
    on_unit_done: Callable[[int], None] | None = None,
    on_linker_done: Callable[[int], None] | None = None,
) -> PlacementReport:
    """Move folded domains so each linker IDR can reach its predicted dimensions.

    Mutates ``structure`` in place. Domains move only as members of a rigid unit, and only as
    rigid bodies; the transform is verified unchanged across the whole unit.

    Parameters
    ----------
    structure
        A structure whose regions have already been assigned.
    mode
        Build mode, passed through to :func:`~dodo.construct.dimensions.target_dimensions` to
        get each linker's target end-to-end distance.
    rng
        Random generator, for the orientation perturbation and the search directions.
    max_attempts
        Placement attempts per unit before accepting the least-bad option.
    clash_distance
        Minimum acceptable approach between atoms of different units.
    perturb
        Apply the random orientation perturbation. Turn off for a deterministic,
        purely geometric arrangement.
    assembly
        Pre-computed rigid units. Built here when omitted, which is the normal case; the
        pipeline passes one in so the same grouping can be reported and reused.
    min_length, lock_interfaces, lock_intra_chain_interfaces
        Passed to :func:`~dodo.construct.assembly.find_rigid_units` when building the assembly.
    tolerance
        Worst per-linker error at which a multi-constraint placement is accepted.
    on_unit_done
        Called with ``1`` after each rigid unit is settled, for a caller showing progress.
        Placing one unit can spend the whole attempt budget, so on an assembly this is a real
        wait rather than a formality.
    on_linker_done
        Called with ``1`` after each connecting IDR's target is worked out. That is where the
        time actually goes on a large assembly -- one dimension prediction per linker -- so a
        caller showing progress wants both.

    Returns
    -------
    PlacementReport
        Where every domain and unit ended up, every linker's target against what it got, and
        every interface, locked or not.

    Raises
    ------
    BuildError
        If a chain's domains have not been assigned.
    GeometryError
        If any transform turns out not to be rigid, which would be a bug here.
    """
    for chain in structure.chains:
        if not chain.domains:
            raise BuildError(
                f"Chain {chain.chain_id!r} has no assigned regions. Call assign_regions() "
                f"before repositioning folded domains."
            )

    if assembly is None:
        assembly = find_rigid_units(
            structure,
            min_length=min_length,
            lock_interfaces=lock_interfaces,
            lock_intra_chain_interfaces=lock_intra_chain_interfaces,
        )

    report = PlacementReport(interfaces=list(assembly.interfaces))
    report.notes.extend(assembly.notes)
    if not assembly.units:
        return report

    constraints = _constraints(
        structure, assembly, mode=mode, on_linker_done=on_linker_done
    )
    inter_unit = [c for c in constraints if c.inter_unit]
    components = _components(len(assembly.units), constraints)

    before_unit: dict[int, np.ndarray] = {u.index: u.xyz for u in assembly.units}
    outcome: dict[int, UnitPlacement] = {}

    # Everything that will never move is an obstacle from the outset: each component's anchor
    # unit, and every unit that has no linker to anywhere. Waiting until a unit's turn to add it
    # would let a unit placed early be put straight through geometry that was never going to
    # move -- the same forward-blindness the region builder had to fix separately.
    placed_mask = np.zeros(structure.n_atoms, dtype=bool)
    anchors: dict[int, int] = {}
    for group in components:
        anchor = group[0]
        anchors[anchor] = anchor
        placed_mask |= assembly.units[anchor].atom_mask
        reason = (
            "the earliest rigid unit of its connected group; it defines the frame of reference"
            if len(group) > 1
            else "no linker connects this unit to another, so nothing says where else to put it"
        )
        outcome[anchor] = _unit_outcome(assembly.units[anchor], moved=False, reason=reason)
        if on_unit_done is not None:
            on_unit_done(1)

    for group in components:
        if len(group) == 1:
            continue
        anchor = group[0]
        adjacency: dict[int, list[_Constraint]] = {k: [] for k in group}
        for constraint in inter_unit:
            if constraint.lo_unit in adjacency:
                adjacency[constraint.lo_unit].append(constraint)
                adjacency[constraint.hi_unit].append(constraint)

        settled = {anchor}
        queue: deque[int] = deque([anchor])
        while queue:
            current = queue.popleft()
            for constraint in adjacency[current]:
                nxt = constraint.other_unit(current)
                if nxt in settled:
                    continue
                unit = assembly.units[nxt]
                held = [c for c in adjacency[nxt] if c.other_unit(nxt) in settled]
                obstacles = _obstacle_tree(structure, placed_mask)
                if len(held) == 1:
                    attempts, clashes = _place_one_constraint(
                        structure,
                        unit,
                        held[0],
                        rng=rng,
                        perturb=perturb,
                        obstacles=obstacles,
                        max_attempts=max_attempts,
                        clash_distance=clash_distance,
                    )
                    residual = abs(held[0].separation(structure) - held[0].target)
                else:
                    attempts, clashes, residual = _place_multi_constraint(
                        structure,
                        unit,
                        held,
                        rng=rng,
                        perturb=perturb,
                        obstacles=obstacles,
                        max_attempts=max_attempts,
                        clash_distance=clash_distance,
                        tolerance=tolerance,
                    )
                outcome[nxt] = _unit_outcome(
                    unit,
                    moved=True,
                    attempts=attempts,
                    clashing=clashes > 0,
                    clashing_atoms=clashes,
                    worst_residual=residual,
                    n_constraints=len(held),
                    reason=(
                        f"no clash-free arrangement found in {attempts} attempt(s); kept the "
                        f"one with fewest clashing atoms ({clashes})"
                        if clashes
                        else None
                    ),
                )
                settled.add(nxt)
                placed_mask |= unit.atom_mask
                queue.append(nxt)
                if on_unit_done is not None:
                    on_unit_done(1)

        _close_cycles(
            structure,
            assembly,
            group=group,
            anchor=anchor,
            constraints=[c for c in inter_unit if c.lo_unit in adjacency],
            tolerance=tolerance,
            clash_distance=clash_distance,
            obstacle_mask=placed_mask,
            report=report,
        )

    # Cycle relaxation runs after the tree-placement outcomes were created. Recompute from the
    # final coordinates so neither a cycle move nor a clash against a moved neighbour can leave
    # the report describing stale geometry.
    folded_mask = np.zeros(structure.n_atoms, dtype=bool)
    for unit in assembly.units:
        folded_mask |= unit.atom_mask
    for unit in assembly.units:
        if unit.index not in outcome:
            continue
        clashes = _unit_clash_count(
            unit, assembly, clash_distance, obstacle_mask=folded_mask
        )
        final_outcome = outcome[unit.index]
        final_reason: str | None = final_outcome.reason
        if clashes and not final_outcome.clashing:
            final_reason = (
                f"{final_reason}; final rigid-unit arrangement contains {clashes} clashing atom(s)"
                if final_reason
                else f"final rigid-unit arrangement contains {clashes} clashing atom(s)"
            )
        outcome[unit.index] = replace(
            final_outcome,
            clashing=clashes > 0,
            clashing_atoms=clashes,
            reason=final_reason,
        )

    # Rigidity is the invariant of the whole mechanism, so it is checked over the unit -- a
    # per-domain check passes on exactly the bug that matters.
    for unit in assembly.units:
        verify_rigid(before_unit[unit.index], unit.xyz)

    _fill_report(report, structure, assembly, constraints, outcome)
    return report


def _unit_outcome(
    unit: RigidUnit,
    *,
    moved: bool,
    reason: str | None = None,
    attempts: int = 0,
    clashing: bool = False,
    clashing_atoms: int = 0,
    worst_residual: float | None = None,
    n_constraints: int = 0,
) -> UnitPlacement:
    structure = unit.structure
    return UnitPlacement(
        index=unit.index,
        chains=unit.chain_ids(),
        domains=tuple(
            f"{structure.chains[int(structure.chain_index[d.span.start])].chain_id}:"
            f"{int(structure.residue_number[d.span.start])}-"
            f"{int(structure.residue_number[d.span.stop - 1])}"
            for d in unit.domains
        ),
        n_atoms=unit.n_atoms,
        moved=moved,
        reason=reason,
        attempts=attempts,
        clashing=clashing,
        clashing_atoms=clashing_atoms,
        worst_residual=worst_residual,
        n_constraints=n_constraints,
    )


def _close_cycles(
    structure: Structure,
    assembly: Assembly,
    *,
    group: list[int],
    anchor: int,
    constraints: list[_Constraint],
    tolerance: float,
    clash_distance: float,
    obstacle_mask: np.ndarray,
    report: PlacementReport,
) -> None:
    """Sweep a component until its cycle-closing constraints are satisfied, or stop trying.

    A spanning tree can satisfy every edge it contains and nothing else, so a component with a
    cycle -- two units joined by two independent linker paths, which a homodimer of a two-domain
    protein produces immediately -- arrives here with the closing edge unsatisfied. Block
    coordinate descent fixes it: relax each non-anchor unit against *all* of its constraints in
    turn, and repeat.

    **Skipped entirely when the tree placement already satisfies every constraint**, which is
    every acyclic component and therefore every single-chain input. That is what keeps this
    whole mechanism invisible to the behaviour DODO already had.
    """

    def worst() -> float:
        errors = [abs(c.separation(structure) - c.target) for c in constraints]
        return max(errors) if errors else 0.0

    start = worst()
    if start <= tolerance:
        return

    movable = [k for k in group if k != anchor]
    per_unit = {
        k: [c for c in constraints if k in (c.lo_unit, c.hi_unit)] for k in movable
    }
    previous = start
    rejected = 0
    for sweep in range(MAX_UNIT_RELAX_SWEEPS):
        for k in movable:
            unit = assembly.units[k]
            before = _snapshot(unit)
            clashes_before = _unit_clash_count(
                unit, assembly, clash_distance, obstacle_mask=obstacle_mask
            )
            _relax_unit(structure, unit, per_unit[k], tolerance=tolerance)
            clashes_after = _unit_clash_count(
                unit, assembly, clash_distance, obstacle_mask=obstacle_mask
            )
            if clashes_after > clashes_before:
                _restore(unit, before)
                rejected += 1
        current = worst()
        if current <= tolerance:
            report.notes.append(
                f"units {group} form a cycle; {sweep + 1} relaxation sweep(s) closed it "
                f"(worst linker error {start:.1f} A -> {current:.2f} A)"
            )
            if rejected:
                report.notes.append(
                    f"cycle relaxation rejected {rejected} move(s) that would have introduced "
                    "additional inter-unit clashes"
                )
            return
        # Converged short of the tolerance: the constraints are mutually unsatisfiable, and more
        # sweeps only burn time. Say so with the number rather than looping to the cap.
        if previous - current < 1e-3:
            report.notes.append(
                f"units {group} form a cycle whose linker predictions cannot all be met; "
                f"relaxation settled at a worst linker error of {current:.1f} A after "
                f"{sweep + 1} sweep(s). The spans DODO could not set are reported per linker."
            )
            if rejected:
                report.notes.append(
                    f"cycle relaxation rejected {rejected} move(s) that would have introduced "
                    "additional inter-unit clashes"
                )
            return
        previous = current
    report.notes.append(
        f"units {group} form a cycle that had not converged after {MAX_UNIT_RELAX_SWEEPS} "
        f"sweeps; worst linker error {worst():.1f} A."
    )
    if rejected:
        report.notes.append(
            f"cycle relaxation rejected {rejected} move(s) that would have introduced "
            "additional inter-unit clashes"
        )


def _fill_report(
    report: PlacementReport,
    structure: Structure,
    assembly: Assembly,
    constraints: list[_Constraint],
    outcome: dict[int, UnitPlacement],
) -> None:
    """Turn unit-level outcomes into the per-domain and per-linker reports users read."""
    report.units = [outcome[u.index] for u in assembly.units if u.index in outcome]

    for constraint in constraints:
        achieved = constraint.separation(structure)
        report.linkers.append(
            LinkerOutcome(
                chain_id=constraint.chain_id,
                residues=constraint.residues,
                n_residues=constraint.n_residues,
                target=constraint.target,
                achieved=achieved,
                dictated=not constraint.inter_unit,
                unbridgeable=achieved > constraint.ceiling,
                ceiling=constraint.ceiling,
            )
        )

    # A domain's own linker, if it has one, so the per-domain line can still quote a target and
    # an achieved separation the way it always did. A domain that moved only because its unit
    # moved has none, and says so rather than borrowing another domain's numbers.
    by_domain: dict[int, _Constraint] = {}
    for constraint in constraints:
        if constraint.inter_unit:
            by_domain.setdefault(constraint.hi.span.start, constraint)
            by_domain.setdefault(constraint.lo.span.start, constraint)

    for unit in assembly.units:
        unit_outcome = outcome.get(unit.index)
        if unit_outcome is None:  # pragma: no cover - every unit gets an outcome
            continue
        for domain in unit.domains:
            chain_id = structure.chains[int(structure.chain_index[domain.span.start])].chain_id
            residues = (
                int(structure.residue_number[domain.span.start]),
                int(structure.residue_number[domain.span.stop - 1]),
            )
            own = by_domain.get(domain.span.start) if unit_outcome.moved else None
            report.placements.append(
                DomainPlacement(
                    chain_id=chain_id,
                    residues=residues,
                    moved=unit_outcome.moved,
                    target_separation=own.target if own is not None else None,
                    achieved_separation=(own.separation(structure) if own is not None else None),
                    attempts=unit_outcome.attempts,
                    clashing=unit_outcome.clashing,
                    reason=(
                        unit_outcome.reason
                        if not unit_outcome.moved or unit_outcome.reason
                        else (
                            None if own is not None else f"moved as part of rigid unit {unit.index}"
                        )
                    ),
                    unit=unit.index,
                )
            )
    report.placements.sort(key=lambda p: (p.chain_id, p.residues))
