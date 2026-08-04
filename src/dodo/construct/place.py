"""Reposition folded domains so linker IDRs can adopt their predicted dimensions.

This is step 3 of DODO's algorithm, and it is the step that makes the whole thing work.

AlphaFold places folded domains at whatever separation its prediction happened to produce,
which for a long disordered linker is essentially arbitrary. If a linker's predicted end-to-end
distance is 100 A and AlphaFold left its flanking domains 300 A apart, then rebuilding the
linker *in place* can only produce a stretched rod -- the chain is forced to span a distance its
sequence says it should not. The fix is not a better linker builder. It is to **move the
domains** so their separation matches the prediction, and only then build the linker.

Getting this backwards is exactly what an earlier version of this rewrite did, and the output
looked wrong for precisely this reason.

What moves and what does not
----------------------------
Folded-domain atoms are **never rebuilt**. They come from AlphaFold (or AlphaFold3, or a
crystal structure) and are trusted. A domain only ever moves as a **rigid body**: its internal
geometry must come out bit-identical, while its position and orientation change. Every
transform here is a rotation plus a translation applied to whole domains, and
:func:`verify_rigid` exists to assert that property rather than assume it.

The first domain of a chain is left exactly where it is, so the output stays in the input's
frame of reference. Subsequent domains are placed relative to their predecessor.

Orientation
-----------
v1 translated domains without rotating them, which can leave a linker having to wrap around a
domain to reach an entry point facing the wrong way. Here each domain is first rotated so that
its body extends *away* from the domain it connects back to -- its entry alpha carbon becomes
the closest point to the previous domain's exit -- and then given a random perturbation, with
clash rejection. That keeps linkers on the outside of the domains they join.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import pairwise
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from ..constants import (
    CA_CA_BOND_LENGTH,
    CA_CLASH_DISTANCE,
    DEFAULT_MODE,
    MAX_FD_PLACEMENT_ATTEMPTS,
)
from ..exceptions import BuildError, GeometryError
from ..geometry.transforms import rotation_between_vectors
from ..structure import Chain, Domain, DomainKind, Structure
from .dimensions import DimensionTarget, target_dimensions

if TYPE_CHECKING:
    from .pipeline import RegionOutcome

__all__ = [
    "DomainPlacement",
    "PlacementReport",
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
    #: Target anchor separation in Angstroms, from the linker's predicted dimensions.
    target_separation: float | None = None
    achieved_separation: float | None = None
    attempts: int = 0
    clashing: bool = False
    reason: str | None = None

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
        flag = " [CLASHING]" if self.clashing else ""
        return f"{where}: moved in {self.attempts} attempt(s){detail}{flag}"


@dataclass(slots=True)
class PlacementReport:
    """Outcome of repositioning every folded domain in a structure."""

    placements: list[DomainPlacement] = field(default_factory=list)
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
    def ok(self) -> bool:
        """True if every moved domain was placed without a clash."""
        return not self.clashing

    def summary(self) -> str:
        """Multi-line human-readable summary."""
        lines = [f"{len(self.moved)}/{len(self.placements)} folded domain(s) repositioned"]
        lines += [f"  {p}" for p in self.placements]
        lines += [f"  note: {n}" for n in self.notes]
        return "\n".join(lines)


def verify_rigid(before: np.ndarray, after: np.ndarray, *, tolerance: float = 1e-6) -> None:
    """Assert that a domain's internal geometry survived a transform unchanged.

    Folded-domain atoms are never rebuilt, so any transform applied to one must be a rigid
    motion: every internal pairwise distance is preserved. This checks that directly rather than
    trusting the arithmetic, because a subtly non-orthonormal rotation would deform a real
    domain while still looking plausible in a viewer.

    Compares the distance from every atom to the domain's own centroid, which is O(n) and
    detects any scaling, shearing or reflection.

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
    radii_before = np.linalg.norm(before - before.mean(axis=0), axis=1)
    radii_after = np.linalg.norm(after - after.mean(axis=0), axis=1)
    drift = float(np.abs(radii_before - radii_after).max())
    if drift > tolerance:
        raise GeometryError(
            f"A folded domain was deformed rather than moved rigidly: the worst atom's "
            f"distance to the domain centroid changed by {drift:.3e} A. Folded-domain atoms "
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
    try:
        index = ordered.index(first)
    except ValueError:  # pragma: no cover - callers pass domains from this chain
        return None
    if index + 2 >= len(ordered) or ordered[index + 2] is not second:
        return None
    middle = ordered[index + 1]
    return middle if middle.kind is DomainKind.IDR else None


def _orient_away_from(
    structure: Structure,
    domain: Domain,
    *,
    from_point: np.ndarray,
    rng: np.random.Generator,
    perturb: bool,
) -> None:
    """Rotate ``domain`` in place so its body extends away from ``from_point``.

    The construction: the vector from the domain's entry CA to its centroid should point in the
    same direction as the vector from ``from_point`` to that entry CA. Then the entry CA is the
    closest part of the domain to the previous one, and the linker arrives at the near face
    rather than having to wrap around to a far one.

    Rotation is about the entry CA, so that point stays put and the subsequent translation can
    place it exactly.
    """
    entry = _entry_ca(structure, domain)
    centroid = domain.centroid()

    body = centroid - entry
    if float(np.linalg.norm(body)) < 1e-9:
        # A single-residue domain has no body direction to align. Nothing to orient.
        return
    approach = entry - np.asarray(from_point, dtype=np.float64)
    if float(np.linalg.norm(approach)) < 1e-9:
        # Entry sits on the previous exit; any orientation is as good as another.
        return

    domain.rotate(rotation_between_vectors(body, approach), about=entry)

    if perturb:
        # A small random rotation about the entry point, so repeated runs differ without
        # abandoning the facing constraint just established.
        axis = rng.normal(size=3)
        axis /= np.linalg.norm(axis)
        angle = np.deg2rad(rng.uniform(-PERTURBATION_DEGREES, PERTURBATION_DEGREES))
        from ..geometry.transforms import rotation_from_axis_angle

        domain.rotate(rotation_from_axis_angle(axis, float(angle)), about=entry)


def _place_at_separation(
    structure: Structure,
    domain: Domain,
    *,
    from_point: np.ndarray,
    separation: float,
    direction: np.ndarray,
) -> None:
    """Translate ``domain`` so its entry CA sits ``separation`` A from ``from_point``."""
    unit = np.asarray(direction, dtype=np.float64)
    norm = float(np.linalg.norm(unit))
    if norm < 1e-9:
        raise GeometryError("Cannot place a domain along a zero-length direction.")
    unit = unit / norm
    wanted = np.asarray(from_point, dtype=np.float64) + unit * separation
    domain.translate(wanted - _entry_ca(structure, domain))


def _placed_atom_mask(structure: Structure, placed: list[Domain]) -> np.ndarray:
    """Atom mask over the domains already positioned, for clash checking."""
    mask = np.zeros(structure.n_atoms, dtype=bool)
    for domain in placed:
        mask[domain.atom_slice] = True
    return mask


def place_domain_after(
    structure: Structure,
    *,
    region: Domain,
    downstream: Domain | None,
    engine: object,
    rng: np.random.Generator,
    model: int,
    mode: str,
    clash_distance: float = CA_CLASH_DISTANCE,
) -> RegionOutcome:
    """Build ``region`` from a generated conformer, then move ``downstream`` to meet its far end.

    The conformer-driven counterpart to :func:`reposition_folded_domains`. That function moves a
    domain to a *predicted* separation and leaves a region to be built into the gap; this one takes
    the region's own conformer as given and moves the domain to wherever it ends.

    Each candidate conformer is attached to the upstream anchor -- which is already final, since the
    chain is walked N to C -- and the downstream domain is then translated so its entry alpha carbon
    sits one virtual bond from the conformer's last one. A candidate is rejected if that leaves the
    domain colliding with anything already placed, and the next is tried. Dimension never enters it.

    ``downstream=None`` covers a terminal IDR, where there is no domain beyond the region to move.
    The conformer is still taken as generated rather than selected for its size, which is the
    consistent thing to do: under this mode a tail's length is STARLING's statement about that
    sequence, and filtering it would reintroduce exactly the selection the mode exists to avoid --
    on one region it narrowed the achieved spread from 44.4 A to 11.3 A.
    """
    from ..construct.pipeline import RegionOutcome
    from ..engines.starling import place_between_anchors

    chain = structure.chains[int(structure.chain_index[region.span.start])]
    residues = (region.span.start + 1, region.span.stop)
    sequence = structure.sequence[region.span.slice]

    anchor = region.span.n_anchor
    if anchor is None:
        raise GeometryError("place_domain_after needs a region with an N-terminal anchor.")
    anchor_xyz = structure.ca_xyz[anchor]

    try:
        pool, _pool_notes = engine.usable_conformers(  # type: ignore[attr-defined]
            sequence, n_conformations=1, rng=rng
        )
    except Exception as exc:
        return RegionOutcome(
            model=model,
            chain_id=chain.chain_id,
            residues=residues,
            n_residues=len(region.span),
            built=False,
            reason=f"{type(exc).__name__}: {exc}",
        )

    obstacles = structure.xyz[structure.placed_atom_mask()]
    original = None if downstream is None else downstream.xyz.copy()
    for index in range(pool.shape[0]):
        placed = place_between_anchors(
            pool[index],
            n_anchor_xyz=anchor_xyz,
            c_anchor_xyz=None,
            rng=rng,
            obstacles=obstacles,
            conformer_index=index,
        )
        if not placed.ok:
            continue
        coords = placed.ca_coords
        # One virtual bond beyond the conformer's last alpha carbon, continuing its own direction.
        if downstream is not None:
            direction = coords[-1] - coords[-2] if coords.shape[0] >= 2 else coords[-1] - anchor_xyz
            downstream.xyz[:] = original
            _place_at_separation(
                structure,
                downstream,
                from_point=coords[-1],
                separation=CA_CA_BOND_LENGTH,
                direction=direction,
            )
            if _collides(structure, downstream, obstacles, coords, clash_distance):
                continue
        structure.set_ca_xyz(np.arange(region.span.start, region.span.stop), coords)
        region.rebuilt = True
        region.placed = True
        return RegionOutcome(
            model=model,
            chain_id=chain.chain_id,
            residues=residues,
            n_residues=len(region.span),
            built=True,
            achieved_end_to_end=float(np.linalg.norm(coords[-1] - coords[0])),
            reason=(
                "conformer taken as generated"
                if downstream is None
                else "domain positioned to match this conformer"
            ),
        )

    if downstream is not None and original is not None:
        downstream.xyz[:] = original
    return RegionOutcome(
        model=model,
        chain_id=chain.chain_id,
        residues=residues,
        n_residues=len(region.span),
        built=False,
        reason=(
            f"none of {pool.shape[0]} conformer(s) could be placed clear of what was already there"
        ),
    )


def _collides(
    structure: Structure,
    domain: Domain,
    obstacles: np.ndarray,
    region_coords: np.ndarray,
    clash_distance: float,
) -> bool:
    """Return True if a freshly moved domain is too close to the obstacles or its own linker."""
    moved = domain.xyz
    if obstacles.shape[0]:
        tree = cKDTree(obstacles)
        if tree.query(moved)[0].min() < clash_distance:
            return True
    # The linker that positioned it is excluded near its own junction, where proximity is bonding
    # rather than collision -- the entry CA is deliberately one bond from the conformer's last.
    if region_coords.shape[0] > 2:
        tree = cKDTree(region_coords[:-2])
        if tree.query(moved)[0].min() < clash_distance:
            return True
    return False


def reposition_folded_domains(
    structure: Structure,
    *,
    mode: str = DEFAULT_MODE,
    rng: np.random.Generator,
    max_attempts: int = MAX_FD_PLACEMENT_ATTEMPTS,
    clash_distance: float = CA_CLASH_DISTANCE,
    perturb: bool = True,
) -> PlacementReport:
    """Move folded domains so each linker IDR can reach its predicted dimensions.

    Mutates ``structure`` in place. Each domain moves as a rigid body only; internal geometry is
    verified unchanged.

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
        Placement attempts per domain before accepting the least-bad option.
    clash_distance
        Minimum acceptable approach between atoms of different domains.
    perturb
        Apply the random orientation perturbation. Turn off for a deterministic,
        purely geometric arrangement.

    Returns
    -------
    PlacementReport
        Where every domain ended up, with target and achieved separations.

    Raises
    ------
    BuildError
        If a chain's domains have not been assigned.
    GeometryError
        If any transform turns out not to be rigid, which would be a bug here.
    """
    report = PlacementReport()

    for chain in structure.chains:
        if not chain.domains:
            raise BuildError(
                f"Chain {chain.chain_id!r} has no assigned regions. Call assign_regions() "
                f"before repositioning folded domains."
            )

        folded = [
            d
            for d in sorted(chain.domains, key=lambda d: d.span.start)
            if d.kind is DomainKind.FOLDED
        ]
        if len(folded) < 2:
            if folded:
                report.notes.append(
                    f"chain {chain.chain_id}: only one folded domain, so there is no linker "
                    f"separation to set and nothing was moved"
                )
            continue

        placed: list[Domain] = [folded[0]]
        report.placements.append(
            DomainPlacement(
                chain_id=chain.chain_id,
                residues=(
                    int(structure.residue_number[folded[0].span.start]),
                    int(structure.residue_number[folded[0].span.stop - 1]),
                ),
                moved=False,
                reason="first folded domain of the chain; it defines the frame of reference",
            )
        )

        for previous, current in pairwise(folded):
            residues = (
                int(structure.residue_number[current.span.start]),
                int(structure.residue_number[current.span.stop - 1]),
            )
            linker = _linker_between(chain, previous, current)
            if linker is None:
                report.placements.append(
                    DomainPlacement(
                        chain_id=chain.chain_id,
                        residues=residues,
                        moved=False,
                        reason="no IDR between this domain and the previous one, so their "
                        "relative position is not ours to change",
                    )
                )
                placed.append(current)
                continue

            # The target is the linker's predicted end-to-end distance, measured between the
            # exit CA of the previous domain and the entry CA of this one -- i.e. the span the
            # linker has to bridge.
            target: DimensionTarget = target_dimensions(linker.sequence, mode=mode)
            separation = target.end_to_end

            before = current.xyz.copy()
            exit_point = _exit_ca(structure, previous)
            obstacles = _placed_atom_mask(structure, placed)

            best_direction: np.ndarray | None = None
            best_clash_count = np.inf
            attempts = 0

            for attempt in range(1, max_attempts + 1):
                attempts = attempt
                # Restore, then re-place: each attempt is independent, so a rejected
                # orientation cannot bias the next one.
                current.xyz[:] = before

                direction = rng.normal(size=3)
                if float(np.linalg.norm(direction)) < 1e-9:
                    continue
                _place_at_separation(
                    structure,
                    current,
                    from_point=exit_point,
                    separation=separation,
                    direction=direction,
                )
                _orient_away_from(
                    structure, current, from_point=exit_point, rng=rng, perturb=perturb
                )
                # Orientation rotates about the entry CA, which the translation had already
                # positioned, so the separation still holds. Re-place anyway: the rotation is
                # about a point that itself came from the pre-rotation geometry, and re-placing
                # costs nothing while removing any doubt.
                _place_at_separation(
                    structure,
                    current,
                    from_point=exit_point,
                    separation=separation,
                    direction=_entry_ca(structure, current) - exit_point,
                )

                clashes = structure.clash_mask(
                    current.xyz,
                    cutoff=clash_distance,
                    obstacle_atom_mask=obstacles,
                    exclude_within=0,
                )
                clash_count = int(clashes.sum())
                if clash_count == 0:
                    best_direction = None
                    break
                if clash_count < best_clash_count:
                    best_clash_count = clash_count
                    best_direction = direction
            else:
                # Exhausted attempts. Re-apply the least-bad arrangement rather than leaving
                # whatever the final attempt happened to produce, and report the clash instead
                # of hiding it.
                if best_direction is not None:
                    current.xyz[:] = before
                    _place_at_separation(
                        structure,
                        current,
                        from_point=exit_point,
                        separation=separation,
                        direction=best_direction,
                    )
                    _orient_away_from(
                        structure, current, from_point=exit_point, rng=rng, perturb=False
                    )
                    _place_at_separation(
                        structure,
                        current,
                        from_point=exit_point,
                        separation=separation,
                        direction=_entry_ca(structure, current) - exit_point,
                    )

            verify_rigid(before, current.xyz)
            achieved = float(np.linalg.norm(_entry_ca(structure, current) - exit_point))
            still_clashing = best_direction is not None

            report.placements.append(
                DomainPlacement(
                    chain_id=chain.chain_id,
                    residues=residues,
                    moved=True,
                    target_separation=separation,
                    achieved_separation=achieved,
                    attempts=attempts,
                    clashing=still_clashing,
                    reason=(
                        f"no clash-free arrangement found in {max_attempts} attempts; kept the "
                        f"one with fewest clashing atoms ({int(best_clash_count)})"
                        if still_clashing
                        else None
                    ),
                )
            )
            placed.append(current)

    return report
