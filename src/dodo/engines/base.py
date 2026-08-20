"""The conformation engine interface: what a builder asks for, and what it gets back.

One request type, one result type, one protocol. Every engine speaks exactly this, so the
orchestrator can choose one at runtime without knowing anything about how it works. One ships
today -- the self-avoiding walk -- and the protocol is what lets another be added later without
the orchestrator changing.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

from ..exceptions import EngineError, GeometryError

__all__ = [
    "ConformationEngine",
    "IDRRequest",
    "IDRResult",
]


def _as_anchor(value: np.ndarray | None, name: str) -> np.ndarray | None:
    """Coerce an optional anchor to a finite ``(3,)`` float64 array.

    Parameters
    ----------
    value
        Candidate anchor coordinate, or ``None`` for "no anchor on this side".
    name
        Field name, for the error message.

    Returns
    -------
    np.ndarray or None
        The anchor as a fresh float64 array, or ``None``.

    Raises
    ------
    GeometryError
        If the anchor is not a finite 3-vector. A NaN anchor would propagate into every
        candidate position generated from it, and the resulting all-NaN region would be
        reported as a *successful* build.
    """
    if value is None:
        return None
    array = np.array(value, dtype=np.float64).reshape(-1)
    if array.shape != (3,):
        raise GeometryError(
            f"{name} must be a 3-vector of CA coordinates, got shape {array.shape}."
        )
    if not np.all(np.isfinite(array)):
        raise GeometryError(
            f"{name} contains non-finite values ({array.tolist()}). An anchor comes from a "
            f"fixed residue that is already placed, so a NaN here means the caller passed "
            f"coordinates from a region that was never built."
        )
    return array


# ``eq=False`` on both dataclasses below: the default ``__eq__`` compares field tuples,
# and comparing two numpy arrays returns an array, whose truth value is ambiguous. A
# generated ``__eq__`` would therefore raise ValueError instead of answering, so identity
# comparison -- which is what callers actually want for a request object -- is used.


@dataclass(frozen=True, slots=True, eq=False)
class IDRRequest:
    """One region to build, fully specified.

    Attributes
    ----------
    sequence
        One-letter sequence of the region. Carried because sequence-aware engines
        (STARLING, ALBATROSS-derived targets) need it; the self-avoiding walk uses only
        its length.
    n_residues
        Number of residues to generate. Must equal ``len(sequence)``.
    target_end_to_end
        Target CA-to-CA distance between the region's first and last residue, in
        Angstroms. This comes from
        :attr:`dodo.construct.dimensions.DimensionTarget.end_to_end`; engines must not
        invent their own.

        Note what it is *not*: for a region with both anchors fixed, the achievable
        end-to-end distance is dictated by the anchors, not by this number. See
        :attr:`n_anchor_xyz`.
    n_anchor_xyz
        CA coordinate of the fixed residue immediately N-terminal to the region, or
        ``None`` if the region starts at a chain terminus. The anchor lies **outside**
        the region (:attr:`dodo.structure.Span.n_anchor`), so the first generated
        residue must sit one CA-CA bond length away from it -- not on top of it.
    c_anchor_xyz
        CA coordinate of the fixed residue immediately C-terminal to the region, or
        ``None``. The last generated residue must sit one bond length from it.
    n_anchor_prev_xyz
        CA coordinate of the fixed residue one *further* N-terminal, i.e. the neighbour of
        :attr:`n_anchor_xyz` on the far side from the region. Optional, and the reason it
        exists is not decoration: the CA-CA-CA pseudo-angle centred on the N-anchor is
        ``(this residue, n_anchor, first generated residue)``, so without it that angle
        cannot be constrained during generation and cannot be measured afterwards. It is
        an angle of the *output structure* -- the anchor is a residue of the file DODO
        writes -- and leaving it free is how the pre-rewrite builder ended up with 38-177
        degree seams while reporting success. Supply it whenever the anchor has a
        neighbour; only the direction from the anchor to it is used.
    c_anchor_next_xyz
        CA coordinate of the fixed residue one further C-terminal, the mirror image of
        :attr:`n_anchor_prev_xyz`. Constrains and exposes the pseudo-angle centred on the
        C-anchor.
    ensemble_mean_end_to_end
        The end-to-end distance the *region* was asked for, in Angstroms, when
        :attr:`target_end_to_end` is one draw from an ensemble spread around it rather than
        the mean itself. ``None`` means they are the same thing, which is the case whenever
        the engine does its own spreading.

        This exists because a caller may do the spreading instead. :mod:`dodo.construct.pipeline`
        draws a per-model target for every terminal region up front and then calls the engine once
        per model with ``n_conformations=1``, so from inside the engine that draw is
        indistinguishable from a region that was simply asked for a compact target. The
        distinction matters to :func:`dodo.engines.walk._target_steering_width`, which decides
        whether a region is long enough to be worth spending accuracy on -- a property of the
        region, not of one draw. Leave it unset and a compact draw is treated as a small region.
    n_conformations
        How many independent conformers to generate. Used for multi-model output; the
        engine reports one success flag per conformer.

        Note what a batch means physically. A predicted end-to-end distance is the *mean*
        of a distribution, not a per-conformer constraint, so an engine that pins every
        conformer to :attr:`target_end_to_end` returns one conformation sampled
        ``n_conformations`` times. Engines are expected to spread the batch around the
        target and to match it in the mean; see
        :func:`dodo.engines.walk.sample_end_to_end_targets`.

    Notes
    -----
    Anchor coordinates are the CA positions *only*. An engine checks clashes against an
    obstacle set supplied separately, and that obstacle set should **exclude the atoms of
    the two anchor residues themselves**: a correctly placed junction CA sits about 2.4 A
    from the anchor residue's carbonyl carbon by construction, and an engine's clash
    threshold is a CA-to-CA threshold, so those atoms are not obstacles in any useful
    sense.
    """

    sequence: str
    n_residues: int
    target_end_to_end: float
    n_anchor_xyz: np.ndarray | None = None
    c_anchor_xyz: np.ndarray | None = None
    n_anchor_prev_xyz: np.ndarray | None = None
    c_anchor_next_xyz: np.ndarray | None = None
    ensemble_mean_end_to_end: float | None = None
    n_conformations: int = 1

    def __post_init__(self) -> None:
        """Validate the request and normalize the anchor arrays."""
        if self.n_residues < 1:
            raise EngineError(f"n_residues must be at least 1, got {self.n_residues}.")
        if len(self.sequence) != self.n_residues:
            raise EngineError(
                f"sequence has {len(self.sequence)} residues but n_residues is "
                f"{self.n_residues}. These are two statements of the same fact and the "
                f"engine cannot tell which one the caller meant."
            )
        if self.n_conformations < 1:
            raise EngineError(f"n_conformations must be at least 1, got {self.n_conformations}.")
        target = float(self.target_end_to_end)
        if not np.isfinite(target) or target <= 0.0:
            raise EngineError(
                f"target_end_to_end must be positive and finite, got {self.target_end_to_end!r}. "
                f"Take it from dodo.construct.target_dimensions(...).end_to_end."
            )
        # frozen dataclass: normalize through object.__setattr__ so that every engine sees
        # float64 (3,) arrays and none of them has to re-validate.
        object.__setattr__(self, "target_end_to_end", target)
        object.__setattr__(self, "n_anchor_xyz", _as_anchor(self.n_anchor_xyz, "n_anchor_xyz"))
        object.__setattr__(self, "c_anchor_xyz", _as_anchor(self.c_anchor_xyz, "c_anchor_xyz"))
        object.__setattr__(
            self, "n_anchor_prev_xyz", _as_anchor(self.n_anchor_prev_xyz, "n_anchor_prev_xyz")
        )
        object.__setattr__(
            self, "c_anchor_next_xyz", _as_anchor(self.c_anchor_next_xyz, "c_anchor_next_xyz")
        )
        # An outer neighbour with no anchor is not a partially specified junction, it is a
        # contradiction: the only thing the coordinate is for is the pseudo-angle centred
        # on the anchor, and without the anchor there is no such angle. Accepting it
        # silently would leave the junction unconstrained while looking constrained.
        for outer, anchor in (
            ("n_anchor_prev_xyz", "n_anchor_xyz"),
            ("c_anchor_next_xyz", "c_anchor_xyz"),
        ):
            if getattr(self, outer) is not None and getattr(self, anchor) is None:
                raise EngineError(
                    f"{outer} was given without {anchor}. It names the residue beyond the "
                    f"anchor, and its only use is the pseudo-angle centred on that anchor, "
                    f"so there is nothing for it to constrain when the anchor is absent."
                )

    def __repr__(self) -> str:
        anchors = {
            (True, True): "interior",
            (True, False): "C-terminal tail",
            (False, True): "N-terminal tail",
            (False, False): "free",
        }[(self.n_anchor_xyz is not None, self.c_anchor_xyz is not None)]
        return (
            f"IDRRequest({self.n_residues} residues, {anchors}, "
            f"target Re {self.target_end_to_end:.1f} A, "
            f"{self.n_conformations} conformation(s))"
        )

    @property
    def is_interior(self) -> bool:
        """True if both anchors are fixed, so the region is a closure problem."""
        return self.n_anchor_xyz is not None and self.c_anchor_xyz is not None

    @property
    def unconstrained_junctions(self) -> tuple[str, ...]:
        """Names of the anchors whose pseudo-angle this request cannot constrain.

        An anchor with no outer neighbour supplied is a vertex of the assembled chain whose
        CA-CA-CA angle nothing in the pipeline can either steer or measure. Engines report
        this through :attr:`IDRResult.unconstrained_junctions` and warn, because a junction
        angle that is merely unmeasured looks exactly like a junction angle that is fine.
        """
        missing = []
        if self.n_anchor_xyz is not None and self.n_anchor_prev_xyz is None:
            missing.append("n_anchor")
        if self.c_anchor_xyz is not None and self.c_anchor_next_xyz is None:
            missing.append("c_anchor")
        return tuple(missing)

    @property
    def anchor_separation(self) -> float | None:
        """Distance between the two anchors in Angstroms, or ``None`` unless interior.

        This -- not :attr:`n_residues`, and not :attr:`target_end_to_end` -- is what a
        closure schedule must be seeded from. The builder being replaced seeded its
        radius ramp from the residue count, so it fought its own anchors: asked to bridge
        150 A it produced a near-perfect rod, and asked to bridge 20 A it looped out to
        32 A before coming back.
        """
        if self.n_anchor_xyz is None or self.c_anchor_xyz is None:
            return None
        return float(np.linalg.norm(self.c_anchor_xyz - self.n_anchor_xyz))


@dataclass(frozen=True, slots=True, eq=False)
class IDRResult:
    """Generated conformers, with an explicit record of which ones worked.

    Attributes
    ----------
    ca_coords
        ``(n_conformations, n_residues, 3)`` CA coordinates. **Rows for conformers whose
        :attr:`success` entry is False are meaningless** and are required to be
        non-finite; see :meth:`from_batch`.
    success
        ``(n_conformations,)`` boolean. True where that conformer is a complete,
        physically valid region. Never implied, never inferred from the coordinates.
    engine
        Name of the engine that produced this result, for provenance and error messages.
    attempts
        Total attempt rounds spent, including successful ones. A number close to the
        engine's attempt budget is a warning that the request is near-infeasible even
        when everything succeeded.
    relaxed_to
        The clash threshold actually used, in Angstroms, if it had to be relaxed below
        :data:`dodo.constants.CA_CLASH_DISTANCE` at any point; ``None`` if the strict
        threshold held throughout. Reported rather than hidden: a caller that cares about
        contacts needs to know a rung of
        :data:`dodo.constants.CLASH_RELAXATION_LADDER` was used.
    unconstrained_junctions
        Names of the anchors (``"n_anchor"``, ``"c_anchor"``) whose CA-CA-CA pseudo-angle
        was neither steered nor measured, because
        :attr:`IDRRequest.n_anchor_prev_xyz` / :attr:`IDRRequest.c_anchor_next_xyz` were
        not supplied. Empty when every junction angle in the assembled chain was checked.
        This is reported rather than assumed benign: measured over 480 junctions built
        without those coordinates, 41% of the junction angles fell below
        :data:`dodo.constants.BACKBONE_ANGLE_MIN` and 3% above
        :data:`~dodo.constants.BACKBONE_ANGLE_MAX`, with every conformer marked successful.
    """

    ca_coords: np.ndarray
    success: np.ndarray
    engine: str
    attempts: int
    relaxed_to: float | None = None
    unconstrained_junctions: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Validate shapes and enforce the failed-rows-are-not-coordinates rule."""
        coords = np.asarray(self.ca_coords, dtype=np.float64)
        if coords.ndim != 3 or coords.shape[2] != 3:
            raise GeometryError(
                f"ca_coords must have shape (n_conformations, n_residues, 3), got {coords.shape}."
            )
        flags = np.asarray(self.success, dtype=bool)
        if flags.shape != (coords.shape[0],):
            raise EngineError(
                f"success must have one entry per conformation: expected shape "
                f"({coords.shape[0]},), got {flags.shape}."
            )
        finite_rows = np.all(np.isfinite(coords), axis=(1, 2))
        if np.any(flags & ~finite_rows):
            bad = int(np.count_nonzero(flags & ~finite_rows))
            raise GeometryError(
                f"{bad} conformer(s) are marked successful but contain non-finite "
                f"coordinates. A successful conformer is a complete set of finite CA "
                f"positions; anything else is a failure and must be reported as one."
            )
        if np.any(~flags & finite_rows):
            bad = int(np.count_nonzero(~flags & finite_rows))
            raise GeometryError(
                f"{bad} conformer(s) are marked failed but hold finite coordinates. "
                f"Failed conformers must not carry anything that can be mistaken for a "
                f"structure -- writing an unbuilt region as a block of (0, 0, 0) atoms is "
                f"precisely the failure this type exists to prevent. Use "
                f"IDRResult.from_batch, which NaN-fills them for you."
            )
        if self.attempts < 0:
            raise EngineError(f"attempts must be non-negative, got {self.attempts}.")
        object.__setattr__(self, "ca_coords", coords)
        object.__setattr__(self, "success", flags)

    def __repr__(self) -> str:
        relaxed = f", relaxed to {self.relaxed_to:.2f} A" if self.relaxed_to is not None else ""
        loose = (
            f", unchecked junctions: {'+'.join(self.unconstrained_junctions)}"
            if self.unconstrained_junctions
            else ""
        )
        return (
            f"IDRResult({self.n_successful}/{self.n_conformations} built, "
            f"{self.n_residues} residues, engine={self.engine!r}, "
            f"attempts={self.attempts}{relaxed}{loose})"
        )

    @classmethod
    def from_batch(
        cls,
        ca_coords: np.ndarray,
        success: np.ndarray,
        engine: str,
        attempts: int,
        relaxed_to: float | None = None,
        unconstrained_junctions: tuple[str, ...] = (),
    ) -> IDRResult:
        """Build a result, NaN-filling the rows of conformers that failed.

        The one blessed way to construct a partially successful result. NaN rather than
        zeros, and NaN rather than "whatever the partial walk had reached": NaN cannot be
        mistaken for a coordinate, propagates through any arithmetic that touches it, and
        makes :func:`dodo.io.write_pdb` produce an obviously broken file instead of a
        plausible wrong one.

        Parameters
        ----------
        ca_coords
            ``(n_conformations, n_residues, 3)`` coordinates. Rows where ``success`` is
            False are overwritten with NaN and need not be meaningful.
        success
            ``(n_conformations,)`` boolean mask.
        engine
            Engine name.
        attempts
            Attempt rounds spent.
        relaxed_to
            Clash threshold actually used, if relaxed.
        unconstrained_junctions
            Anchors whose junction pseudo-angle could not be checked.

        Returns
        -------
        IDRResult
            The result, with failed rows blanked.
        """
        coords = np.array(ca_coords, dtype=np.float64, copy=True)
        flags = np.asarray(success, dtype=bool)
        if coords.ndim != 3 or flags.shape != (coords.shape[0],):
            # Defer to __post_init__ for the full message rather than duplicating it.
            return cls(coords, flags, engine, attempts, relaxed_to, unconstrained_junctions)
        coords[~flags] = np.nan
        return cls(coords, flags, engine, attempts, relaxed_to, unconstrained_junctions)

    @property
    def n_conformations(self) -> int:
        """Number of conformers in this result, successful or not."""
        return int(self.ca_coords.shape[0])

    @property
    def n_residues(self) -> int:
        """Number of residues per conformer."""
        return int(self.ca_coords.shape[1])

    @property
    def n_successful(self) -> int:
        """How many conformers were built successfully."""
        return int(np.count_nonzero(self.success))

    @property
    def all_successful(self) -> bool:
        """True if every requested conformer was built."""
        return bool(np.all(self.success))

    @property
    def successful_coords(self) -> np.ndarray:
        """``(n_successful, n_residues, 3)`` coordinates of the conformers that worked.

        The safe way to consume a partial result: it cannot silently include a failed
        conformer, because the mask does the selecting.
        """
        selected: np.ndarray = self.ca_coords[self.success]
        return selected

    @property
    def end_to_end_distances(self) -> np.ndarray:
        """Achieved CA-to-CA end-to-end distance of each conformer, in Angstroms.

        ``(n_conformations,)``, NaN for failed conformers and for a single-residue
        region, where the quantity is undefined. Compare against
        :attr:`IDRRequest.target_end_to_end` to report achieved-versus-target error.
        """
        if self.n_residues < 2:
            return np.full(self.n_conformations, np.nan)
        distances: np.ndarray = np.linalg.norm(
            self.ca_coords[:, -1, :] - self.ca_coords[:, 0, :], axis=-1
        )
        return distances

    @property
    def radii_of_gyration(self) -> np.ndarray:
        """Radius of gyration of each conformer's CA trace, in Angstroms.

        ``(n_conformations,)``, NaN for failed conformers. The independent check on shape:
        a conformer can have the right end-to-end distance and still be the wrong shape,
        and for an ideal chain ``Rg == sqrt(mean(Re ** 2) / 6)``.
        """
        if self.n_residues < 2:
            return np.full(self.n_conformations, np.nan)
        centred = self.ca_coords - self.ca_coords.mean(axis=1, keepdims=True)
        radii: np.ndarray = np.sqrt((centred**2).sum(axis=(1, 2)) / self.n_residues)
        return radii


@runtime_checkable
class ConformationEngine(Protocol):
    """What every conformation engine must provide.

    A protocol rather than a base class: the walk is pure geometry and a generative engine
    would wrap a neural network, so implementations need share no code -- only this surface.
    """

    #: Short, stable identifier, recorded in :attr:`IDRResult.engine`.
    name: str

    def available(self) -> bool:
        """Return True if this engine can actually run here.

        Separate from importability on purpose: an engine can be installed and still
        unusable -- STARLING with no model weights downloaded, for instance. Callers
        select engines with this; they do not catch ImportError.
        """
        ...

    def generate(
        self,
        request: IDRRequest,
        obstacles: np.ndarray | None,
        rng: np.random.Generator,
    ) -> IDRResult:
        """Generate conformers for one region.

        Parameters
        ----------
        request
            The region to build.
        obstacles
            ``(n_obstacles, 3)`` coordinates of already-placed atoms that generated
            residues must not collide with, or ``None`` for free space. Obtain it from
            ``structure.xyz[structure.placed_atom_mask()]``, minus the atoms of the two
            anchor residues -- see :class:`IDRRequest`.
        rng
            Seeded generator. Every stochastic decision an engine makes must come from
            here: there is no randomness in DODO that a caller cannot reproduce.

        Returns
        -------
        IDRResult
            Conformers plus the explicit success mask. A partially successful build is
            reported through the mask; a build that produced nothing raises.

        Raises
        ------
        dodo.exceptions.BuildError
            Subclasses thereof, when the request cannot be satisfied at all:
            :class:`~dodo.exceptions.UnsatisfiableTargetError` for a geometrically
            impossible request, :class:`~dodo.exceptions.ExhaustedAttemptsError` when the
            attempt budget ran out. An engine never signals failure by returning
            coordinates.
        """
        ...
