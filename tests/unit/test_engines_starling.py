"""Tests for the STARLING engine and hierarchical long-IDR assembly.

STARLING is not installed in this environment and is not installable in CI (2.4 GB of
weights plus torch), so every STARLING code path here runs against a **fake** injected
through the module's own import seam, ``starling._import_starling``. The fake is not a
mock of the model -- it returns real, geometrically valid CA traces built by the same cone
construction the walk engine uses -- so the adapter is exercised on realistic data. What is
*not* covered is whether the real package's API matches the assumptions the fake encodes;
that is listed explicitly in the module report.

Hierarchical assembly is tested for real, over a locally defined cone-walk engine (and, when
the real ``SelfAvoidingWalk`` is importable, over that too). Everything the brief asks to be
verified about an assembled chain -- bond lengths, pseudo-angles across every junction,
clashes, and the achieved full-length end-to-end distance -- is measured, not asserted by
proxy.
"""

from __future__ import annotations

import contextlib
from dataclasses import dataclass
from itertools import pairwise
from types import SimpleNamespace
from typing import Any

import numpy as np
import pytest

from dodo import constants as C
from dodo.engines import hierarchical as H
from dodo.engines import starling as S
from dodo.engines.base import IDRRequest, IDRResult
from dodo.exceptions import (
    BuildError,
    EngineError,
    EngineUnavailableError,
    ExhaustedAttemptsError,
    GeometryError,
    MissingDependencyError,
    UnsatisfiableTargetError,
)
from dodo.geometry.metrics import ca_bond_lengths, ca_pseudo_angles, validate_ca_trace
from dodo.geometry.transforms import apply, rotation_from_axis_angle

# ---------------------------------------------------------------------------
# Helpers: valid CA traces without needing any engine
# ---------------------------------------------------------------------------


def cone_chain(
    n_residues: int,
    rng: np.random.Generator,
    *,
    target: float | None = None,
    bond: float = C.CA_CA_BOND_LENGTH,
    angle_range: tuple[float, float] | None = None,
    clash_floor: float = C.CA_CLASH_DISTANCE,
) -> np.ndarray:
    """Grow one physically valid CA trace: exact bonds, angles inside the window.

    Deliberately independent of any DODO engine, so a failure in this file points at the
    module under test rather than at a shared helper. ``angle_range`` overrides DODO's
    generation window, which is how a "sub-engine with its own geometry" is simulated;
    ``clash_floor`` lowers the self-avoidance threshold, which is how a sub-engine that lets
    its own CA atoms come too close is simulated.
    """
    if n_residues == 1:
        return np.zeros((1, 3), dtype=np.float64)
    angles = np.deg2rad(
        C.backbone_angle_grid()
        if angle_range is None
        else np.arange(angle_range[0], angle_range[1] + 1.0)
    )
    for _ in range(40):
        coords = np.zeros((n_residues, 3), dtype=np.float64)
        coords[1] = np.array([bond, 0.0, 0.0])
        for i in range(2, n_residues):
            axis = coords[i - 1] - coords[i - 2]
            axis /= np.linalg.norm(axis)
            remaining = n_residues - 1 - i
            # Where the running span should be at this residue if the chain is to finish at
            # `target`: a self-avoiding walk's span grows as (residues) ** nu. Steering
            # toward this rather than merely refusing to overshoot it is what makes the
            # finished chain actually *have* the requested dimension. Filtering alone caps
            # the span from above and never drives it up, so an ensemble built that way sits
            # systematically below its own targets -- which left the fake STARLING unable to
            # answer any request near its nominal dimension.
            scheduled = (
                None if target is None else target * (i / (n_residues - 1)) ** C.FLORY_RE_EXPONENT
            )
            best: np.ndarray | None = None
            best_cost = np.inf
            for _ in range(200 if scheduled is None else 24):
                theta = float(rng.choice(angles))
                perpendicular = np.cross(axis, rng.normal(size=3))
                norm = float(np.linalg.norm(perpendicular))
                if norm < 1e-9:
                    continue
                perpendicular /= norm
                direction = axis * np.cos(np.pi - theta) + perpendicular * np.sin(np.pi - theta)
                candidate = coords[i - 1] + bond * direction
                if (
                    i > 2
                    and float(np.linalg.norm(coords[: i - 2] - candidate, axis=1).min())
                    < clash_floor
                ):
                    continue
                span = float(np.linalg.norm(candidate - coords[0]))
                if scheduled is None:
                    best = candidate
                    break
                # Reachability funnel: reject a step from which the target can no longer be
                # reached even by folding all the way back.
                if span - remaining * bond > target:
                    continue
                cost = abs(span - scheduled)
                if cost < best_cost:
                    best, best_cost = candidate, cost
                if cost <= 0.5:
                    break
            if best is None:
                # A dead end: the funnel and the self-avoidance filter left nothing. Start
                # the chain again rather than forcing a step that violates one of them.
                break
            coords[i] = best
        else:
            return coords
    raise AssertionError(  # pragma: no cover - 40 restarts effectively never all fail
        f"cone_chain could not grow {n_residues} residues to a target of {target}"
    )


class ConeWalk:
    """A minimal ConformationEngine that grows valid free chains toward a target Re.

    Stands in for ``SelfAvoidingWalk`` so that hierarchical assembly is testable in this
    tree regardless of what the walk agent has landed. It satisfies the protocol DODO's
    engines are written against and nothing more.
    """

    name: str = "cone-walk"

    def __init__(self, *, tolerance: float = 6.0, attempts: int = 60) -> None:
        self.tolerance = tolerance
        self.attempts = attempts
        self.calls = 0

    def available(self) -> bool:
        return True

    def generate(
        self, request: IDRRequest, obstacles: np.ndarray | None, rng: np.random.Generator
    ) -> IDRResult:
        self.calls += 1
        coords = np.empty((request.n_conformations, request.n_residues, 3), dtype=np.float64)
        used = 0
        for slot in range(request.n_conformations):
            best: np.ndarray | None = None
            best_error = np.inf
            for _ in range(self.attempts):
                used += 1
                chain = cone_chain(request.n_residues, rng, target=request.target_end_to_end)
                error = abs(float(np.linalg.norm(chain[-1] - chain[0])) - request.target_end_to_end)
                if error < best_error:
                    best, best_error = chain, error
                if best_error <= self.tolerance:
                    break
            assert best is not None
            coords[slot] = best
        return IDRResult.from_batch(
            ca_coords=coords,
            success=np.ones(request.n_conformations, dtype=bool),
            engine=self.name,
            attempts=used,
        )


def span_of(chain: np.ndarray) -> float:
    """End-to-end distance of one ``(n, 3)`` CA trace, in Angstroms."""
    return float(np.linalg.norm(chain[-1] - chain[0]))


def request_for(
    n_residues: int,
    *,
    target: float,
    n_anchor: np.ndarray | None = None,
    c_anchor: np.ndarray | None = None,
    n_conformations: int = 1,
    sequence: str | None = None,
) -> IDRRequest:
    """Build an IDRRequest with every field spelled out."""
    return IDRRequest(
        sequence=("G" * n_residues) if sequence is None else sequence,
        n_residues=n_residues,
        target_end_to_end=target,
        n_anchor_xyz=n_anchor,
        c_anchor_xyz=c_anchor,
        n_conformations=n_conformations,
    )


# ---------------------------------------------------------------------------
# The fake STARLING
# ---------------------------------------------------------------------------


@dataclass
class FakeEnsemble:
    """Stands in for STARLING's Ensemble object: coordinates behind an attribute."""

    coordinates: np.ndarray


def fake_starling(
    *,
    n_conformers: int = 24,
    seed: int = 0,
    attribute: str = "coordinates",
    accept_seed: bool = True,
    scale: float = 1.0,
    as_dict: bool = True,
    fail_with: Exception | None = None,
    max_length: int | None = None,
    weights_path: str | None = None,
    record: dict[str, Any] | None = None,
) -> SimpleNamespace:
    """Build a module-like fake whose ``generate`` returns valid CA traces.

    Every knob corresponds to one assumption the adapter makes about STARLING's real API,
    so the tests can vary them independently: the name of the coordinate attribute, whether
    the signature takes a seed, the units, whether a mapping or a bare ensemble comes back.
    """

    def generate(
        sequence: str,
        conformations: int = 200,
        return_data: bool = True,
        verbose: bool = False,
        show_progress_bar: bool = True,
        **kwargs: Any,
    ) -> Any:
        if record is not None:
            record.update(sequence=sequence, conformations=conformations, kwargs=dict(kwargs))
        if fail_with is not None:
            raise fail_with
        rng = np.random.default_rng(seed + int(kwargs.get("seed", 0)))
        # A spread of dimensions around the polymer-scaling estimate, which is what a real
        # ensemble looks like and what conformer selection needs to have any choice. The
        # spread is deliberately narrower than it once was: 0.6-1.6x of the prediction from
        # a handful of conformers leaves gaps of several Angstroms in the ensemble's own
        # dimension distribution, so no conformer lands near the request and the engine has
        # to report failure. A real 100-conformer STARLING ensemble is far denser than that.
        predicted = C.flory_end_to_end(len(sequence))
        stack = np.stack(
            [
                cone_chain(
                    len(sequence),
                    rng,
                    target=float(
                        np.clip(
                            predicted * rng.uniform(0.75, 1.35),
                            0.0,
                            0.9 * C.contour_length(len(sequence)),
                        )
                    ),
                )
                for _ in range(min(conformations, n_conformers))
            ]
        )
        payload: Any = FakeEnsemble(coordinates=stack * scale)
        if attribute != "coordinates":
            payload = SimpleNamespace(**{attribute: stack * scale})
        return {sequence: payload} if as_dict else payload

    if not accept_seed:

        def generate_no_seed(
            sequence: str,
            conformations: int = 200,
            return_data: bool = True,
        ) -> Any:
            return generate(sequence, conformations, return_data)

        entry = generate_no_seed
    else:
        entry = generate

    module = SimpleNamespace(generate=entry)
    if max_length is not None:
        module.configs = SimpleNamespace(MAX_SEQUENCE_LENGTH=max_length)
    if weights_path is not None:
        module.DEFAULT_MODEL_WEIGHTS_PATH = weights_path
    return module


@pytest.fixture
def installed_starling(monkeypatch: pytest.MonkeyPatch):
    """Install a fake STARLING through the module's seams and return a factory."""

    def install(**kwargs: Any) -> SimpleNamespace:
        module = fake_starling(**kwargs)
        monkeypatch.setattr(S, "_import_starling", lambda: module)
        monkeypatch.setattr(S, "starling_installed", lambda: True)
        monkeypatch.setattr(S, "_loaded_starling", lambda: module)
        return module

    return install


# ---------------------------------------------------------------------------
# Availability and the optional-dependency contract
# ---------------------------------------------------------------------------


def test_starling_is_not_installed_here() -> None:
    # The premise of the whole fake-injection approach. If this ever fails, the STARLING
    # paths should be tested for real instead.
    assert S.starling_installed() is False


def test_module_does_not_import_starling_or_torch_at_module_scope() -> None:
    import subprocess
    import sys

    code = (
        "import sys; import dodo.engines.starling, dodo.engines.hierarchical; "
        "print([m for m in ('starling', 'torch', 'sparrow') if m in sys.modules])"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "[]"


def test_available_is_false_and_quiet_without_starling() -> None:
    assert StarlingEngineFactory().available() is False


def StarlingEngineFactory(**kwargs: Any) -> S.StarlingEngine:  # noqa: N802
    """Construct the engine with test-sized defaults."""
    kwargs.setdefault("ensemble_size", 16)
    kwargs.setdefault("oversample", 2)
    return S.StarlingEngine(**kwargs)


def test_generate_raises_missing_dependency_naming_the_extra() -> None:
    engine = StarlingEngineFactory()
    request = request_for(20, target=30.0)
    with pytest.raises(MissingDependencyError) as excinfo:
        engine.generate(request, None, np.random.default_rng(0))
    assert excinfo.value.extra == "starling"
    assert "dodo[starling]" in str(excinfo.value)


def test_available_true_when_module_is_importable(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: None)
    assert StarlingEngineFactory().available() is True


def test_available_false_when_loaded_module_has_no_entry_point(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: SimpleNamespace())
    assert StarlingEngineFactory().available() is False


def test_available_false_when_advertised_weights_are_absent(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    module = fake_starling(weights_path=str(tmp_path / "definitely-absent.pt"))
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    assert StarlingEngineFactory().available() is False


def test_generate_raises_engine_unavailable_when_weights_are_absent(
    monkeypatch: pytest.MonkeyPatch, tmp_path
) -> None:
    module = fake_starling(weights_path=str(tmp_path / "absent.pt"))
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    with pytest.raises(EngineUnavailableError, match="weights"):
        StarlingEngineFactory().generate(
            request_for(20, target=30.0), None, np.random.default_rng(0)
        )


def test_unrecognized_api_raises_rather_than_guessing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(S, "_import_starling", lambda: SimpleNamespace())
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    with pytest.raises(EngineUnavailableError, match="generate"):
        StarlingEngineFactory().generate(
            request_for(20, target=30.0), None, np.random.default_rng(0)
        )


def test_max_length_falls_back_to_the_constant_when_absent() -> None:
    cap = S.starling_max_length()
    assert cap.value == C.STARLING_MAX_LENGTH
    assert "dodo.constants" in cap.source


def test_max_length_is_read_from_starling_when_it_exposes_one(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = fake_starling(max_length=512)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    cap = S.starling_max_length()
    assert cap.value == 512
    assert cap.source == "starling.configs.MAX_SEQUENCE_LENGTH"


# ---------------------------------------------------------------------------
# Request validation: no silent nonsense
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"n_residues": 0, "target": 10.0}, "at least 1"),
        ({"n_residues": 5, "target": -1.0}, "positive and finite"),
        ({"n_residues": 5, "target": 10.0, "n_conformations": 0}, "at least 1"),
        ({"n_residues": 5, "target": 10.0, "sequence": "GGG"}, "n_residues is"),
    ],
)
def test_bad_requests_raise(kwargs: dict[str, Any], match: str) -> None:
    # Validation lives on IDRRequest itself, so no engine can be handed one of these. This
    # test documents that the engines under test rely on that and do not re-check it.
    with pytest.raises(EngineError, match=match):
        request_for(**kwargs)


def test_non_finite_anchor_is_rejected() -> None:
    with pytest.raises(GeometryError):
        request_for(5, target=10.0, n_anchor=np.array([0.0, np.nan, 0.0]))


def test_generate_rejects_a_seed_integer_instead_of_a_generator(
    installed_starling,
) -> None:
    installed_starling()
    with pytest.raises(TypeError, match="default_rng"):
        StarlingEngineFactory().generate(request_for(20, target=30.0), None, 12345)  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# The anchor construction
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("offset_bonds", [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0])
def test_placement_puts_termini_exactly_one_bond_from_each_anchor(offset_bonds: float) -> None:
    # Sweep the whole feasible window rather than arbitrary separations: with the
    # conformer's own span as the reference, an offset of +/-2 bond lengths is exactly the
    # boundary of what the two bond spheres can reach (cos t = -1 and +1), and 0 is the
    # collinear middle. Every value in between must be hit to machine precision.
    rng = np.random.default_rng(7)
    chain = cone_chain(30, rng, target=40.0)
    span = float(np.linalg.norm(chain[-1] - chain[0]))
    separation = span + offset_bonds * C.CA_CA_BOND_LENGTH
    anchor_n = np.zeros(3)
    anchor_c = np.array([separation, 0.0, 0.0])
    placement = S.place_between_anchors(
        chain, n_anchor_xyz=anchor_n, c_anchor_xyz=anchor_c, rng=rng
    )
    assert placement.ok
    assert placement.n_anchor_gap == pytest.approx(C.CA_CA_BOND_LENGTH, abs=1e-9)
    assert placement.c_anchor_gap == pytest.approx(C.CA_CA_BOND_LENGTH, abs=1e-9)
    assert placement.anchor_residual < 1e-9
    assert placement.achieved_end_to_end == pytest.approx(span, abs=1e-9)


def test_selection_finds_a_conformer_for_a_given_anchor_separation() -> None:
    rng = np.random.default_rng(21)
    ensemble = np.stack([cone_chain(30, rng, target=float(t)) for t in np.linspace(8.0, 90.0, 14)])
    for separation in (10.0, 25.0, 50.0):
        anchor_c = np.array([separation, 0.0, 0.0])
        desired, window = S.desired_internal_span(separation, np.zeros(3), anchor_c)
        order = S.rank_conformers(ensemble, desired, window=window)
        placement = S.place_between_anchors(
            ensemble[order[0]],
            n_anchor_xyz=np.zeros(3),
            c_anchor_xyz=anchor_c,
            rng=rng,
            desired_end_to_end=desired,
        )
        assert placement.ok, f"no conformer could bridge {separation} A"


def test_placement_preserves_internal_geometry_exactly() -> None:
    rng = np.random.default_rng(11)
    chain = cone_chain(40, rng, target=50.0)
    span = float(np.linalg.norm(chain[-1] - chain[0]))
    placement = S.place_between_anchors(
        chain,
        n_anchor_xyz=np.zeros(3),
        c_anchor_xyz=np.array([span, 0.0, 0.0]),
        rng=rng,
    )
    before = validate_ca_trace(chain)
    after = validate_ca_trace(placement.ca_coords)
    assert before.ok and after.ok
    np.testing.assert_allclose(after.bond_lengths, before.bond_lengths, atol=1e-9)
    np.testing.assert_allclose(after.pseudo_angles, before.pseudo_angles, atol=1e-6)


def test_placement_never_puts_a_residue_on_top_of_an_anchor() -> None:
    # The specific pre-rewrite defect: first/last generated residue written onto the
    # anchor CA, giving 0.00 A coincident atoms at every junction.
    rng = np.random.default_rng(3)
    for separation in (5.0, 15.0, 35.0, 60.0, 120.0):
        chain = cone_chain(25, rng, target=separation)
        placement = S.place_between_anchors(
            chain,
            n_anchor_xyz=np.zeros(3),
            c_anchor_xyz=np.array([separation, 0.0, 0.0]),
            rng=rng,
        )
        assert placement.n_anchor_gap is not None and placement.c_anchor_gap is not None
        assert placement.n_anchor_gap > C.CA_CLASH_DISTANCE
        assert placement.c_anchor_gap > C.CA_CLASH_DISTANCE


def test_unreachable_anchors_are_reported_not_forced() -> None:
    rng = np.random.default_rng(5)
    chain = cone_chain(20, rng, target=25.0)
    span = float(np.linalg.norm(chain[-1] - chain[0]))
    placement = S.place_between_anchors(
        chain,
        n_anchor_xyz=np.zeros(3),
        c_anchor_xyz=np.array([span + 50.0, 0.0, 0.0]),
        rng=rng,
    )
    assert not placement.ok
    assert placement.anchor_residual > C.CA_CA_BOND_TOLERANCE
    assert placement.n_anchor_gap == pytest.approx(C.CA_CA_BOND_LENGTH, abs=1e-9)
    assert any("cannot bridge" in note for note in placement.notes)
    # And the coordinates are still a real conformation, not zeros or NaN.
    assert np.all(np.isfinite(placement.ca_coords))
    assert validate_ca_trace(placement.ca_coords).ok


def test_single_anchor_tail_is_placed_one_bond_out() -> None:
    rng = np.random.default_rng(13)
    chain = cone_chain(15, rng, target=25.0)
    for n_anchor, c_anchor in ((np.zeros(3), None), (None, np.zeros(3))):
        placement = S.place_between_anchors(
            chain, n_anchor_xyz=n_anchor, c_anchor_xyz=c_anchor, rng=rng
        )
        assert placement.ok
        assert placement.anchor_residual < 1e-9


def test_single_residue_region_is_placed_between_its_anchors() -> None:
    rng = np.random.default_rng(17)
    placement = S.place_between_anchors(
        np.zeros((1, 3)),
        n_anchor_xyz=np.zeros(3),
        c_anchor_xyz=np.array([4.0, 0.0, 0.0]),
        rng=rng,
    )
    assert placement.ok
    assert placement.n_anchor_gap == pytest.approx(C.CA_CA_BOND_LENGTH, abs=1e-9)
    assert placement.c_anchor_gap == pytest.approx(C.CA_CA_BOND_LENGTH, abs=1e-9)


def test_coincident_anchors_are_flagged_rather_than_dividing_by_zero() -> None:
    rng = np.random.default_rng(19)
    chain = cone_chain(12, rng, target=15.0)
    placement = S.place_between_anchors(
        chain, n_anchor_xyz=np.zeros(3), c_anchor_xyz=np.zeros(3), rng=rng
    )
    assert np.all(np.isfinite(placement.ca_coords))
    assert any("coincident" in note for note in placement.notes)


def test_placement_avoids_obstacles_and_reports_relaxation() -> None:
    rng = np.random.default_rng(23)
    chain = cone_chain(30, rng, target=40.0)
    span = float(np.linalg.norm(chain[-1] - chain[0]))
    anchor_c = np.array([span, 0.0, 0.0])
    # A wall of obstacles on one side of the anchor axis: some orientations clash, some do
    # not, so a clash-free placement exists and must be found.
    grid = np.stack(
        np.meshgrid(np.linspace(0.0, span, 12), np.linspace(-20.0, 20.0, 12), np.array([12.0])),
        axis=-1,
    ).reshape(-1, 3)
    placement = S.place_between_anchors(
        chain,
        n_anchor_xyz=np.zeros(3),
        c_anchor_xyz=anchor_c,
        rng=rng,
        obstacles=grid,
    )
    assert placement.ok
    from scipy.spatial import cKDTree

    assert float(cKDTree(grid).query(placement.ca_coords, k=1)[0].min()) >= C.CA_CLASH_DISTANCE


def test_placement_reports_failure_when_boxed_in_on_all_sides() -> None:
    rng = np.random.default_rng(29)
    chain = cone_chain(20, rng, target=30.0)
    span = float(np.linalg.norm(chain[-1] - chain[0]))
    # Obstacles filling the whole neighbourhood: no orientation can be clash-free.
    axis = np.linspace(-30.0, span + 30.0, 25)
    box = np.stack(
        np.meshgrid(axis, np.linspace(-30.0, 30.0, 25), np.linspace(-30.0, 30.0, 25)),
        axis=-1,
    ).reshape(-1, 3)
    placement = S.place_between_anchors(
        chain,
        n_anchor_xyz=np.zeros(3),
        c_anchor_xyz=np.array([span, 0.0, 0.0]),
        rng=rng,
        obstacles=box,
        orientation_attempts=4,
    )
    assert not placement.ok
    assert not placement.clash_free
    assert any("cleared the obstacles" in note for note in placement.notes)
    assert np.all(np.isfinite(placement.ca_coords))


def test_desired_internal_span_clips_into_the_anchor_window() -> None:
    anchor_n = np.zeros(3)
    anchor_c = np.array([100.0, 0.0, 0.0])
    low = 100.0 - 2 * C.CA_CA_BOND_LENGTH
    high = 100.0 + 2 * C.CA_CA_BOND_LENGTH
    assert S.desired_internal_span(10.0, anchor_n, anchor_c) == (low, (low, high))
    assert S.desired_internal_span(1000.0, anchor_n, anchor_c) == (high, (low, high))
    assert S.desired_internal_span(95.0, anchor_n, anchor_c) == (95.0, (low, high))
    # No window without two anchors: a tail is free.
    assert S.desired_internal_span(55.0, anchor_n, None) == (55.0, None)


def test_rank_conformers_prefers_feasible_then_closest() -> None:
    rng = np.random.default_rng(31)
    chains = [cone_chain(20, rng, target=t) for t in (10.0, 30.0, 50.0)]
    spans = [float(np.linalg.norm(c[-1] - c[0])) for c in chains]
    stack = np.stack(chains)
    order = S.rank_conformers(stack, spans[1])
    assert order[0] == 1
    # With a window excluding the middle conformer, it must come last despite being the
    # closest to the desired value.
    window = (spans[0] - 1.0, spans[0] + 1.0)
    order = S.rank_conformers(stack, spans[1], window=window)
    assert order[0] == 0


# ---------------------------------------------------------------------------
# The STARLING adapter, against the fake
# ---------------------------------------------------------------------------


def test_generate_places_conformers_between_anchors(installed_starling) -> None:
    installed_starling(n_conformers=32)
    engine = StarlingEngineFactory()
    anchor_n = np.zeros(3)
    anchor_c = np.array([45.0, 0.0, 0.0])
    request = request_for(60, target=45.0, n_anchor=anchor_n, c_anchor=anchor_c)
    result, report = engine.generate_detailed(request, None, np.random.default_rng(0))
    assert result.ca_coords.shape == (1, 60, 3)
    assert bool(result.success[0])
    assert result.engine == "starling"
    assert report.worst_anchor_residual < 1e-9
    assert report.kept > 0
    gap_n = float(np.linalg.norm(result.ca_coords[0, 0] - anchor_n))
    gap_c = float(np.linalg.norm(result.ca_coords[0, -1] - anchor_c))
    assert gap_n == pytest.approx(C.CA_CA_BOND_LENGTH, abs=1e-9)
    assert gap_c == pytest.approx(C.CA_CA_BOND_LENGTH, abs=1e-9)


def test_generate_is_reproducible_for_one_seed(installed_starling) -> None:
    installed_starling(n_conformers=16)
    engine = StarlingEngineFactory()
    request = request_for(40, target=35.0, n_anchor=np.zeros(3), c_anchor=np.array([35.0, 0, 0]))
    first = engine.generate(request, None, np.random.default_rng(99))
    second = engine.generate(request, None, np.random.default_rng(99))
    np.testing.assert_array_equal(first.ca_coords, second.ca_coords)
    third = engine.generate(request, None, np.random.default_rng(100))
    assert not np.array_equal(first.ca_coords, third.ca_coords)


def test_multiple_conformations_are_distinct(installed_starling) -> None:
    installed_starling(n_conformers=64)
    engine = StarlingEngineFactory(ensemble_size=64)
    request = request_for(40, target=C.flory_end_to_end(40), n_conformations=5)
    result = engine.generate(request, None, np.random.default_rng(1))
    assert result.ca_coords.shape == (5, 40, 3)
    assert result.success.all()
    spans = np.linalg.norm(result.ca_coords[:, -1] - result.ca_coords[:, 0], axis=1)
    assert len(np.unique(np.round(spans, 6))) > 1


def test_seed_is_passed_through_when_the_signature_accepts_it(
    installed_starling,
) -> None:
    record: dict[str, Any] = {}
    installed_starling(record=record)
    engine = StarlingEngineFactory()
    engine.generate(request_for(30, target=30.0), None, np.random.default_rng(4))
    assert "seed" in record["kwargs"]


def test_missing_seed_parameter_is_reported_not_ignored(installed_starling) -> None:
    installed_starling(accept_seed=False)
    engine = StarlingEngineFactory()
    _, report = engine.generate_detailed(
        request_for(30, target=30.0), None, np.random.default_rng(4)
    )
    assert any("no seed" in note for note in report.notes)
    assert any("show_progress_bar" in note or "verbose" in note for note in report.notes)


def test_nanometre_coordinates_are_detected_and_converted(installed_starling) -> None:
    installed_starling(scale=0.1)
    engine = StarlingEngineFactory()
    _, report = engine.generate_detailed(
        request_for(30, target=30.0), None, np.random.default_rng(6)
    )
    assert any("nanometres" in note for note in report.notes)


def test_coordinates_in_no_recognizable_unit_raise(installed_starling) -> None:
    installed_starling(scale=3.0)
    with pytest.raises(EngineError, match="neither"):
        StarlingEngineFactory().generate(
            request_for(30, target=30.0), None, np.random.default_rng(6)
        )


def test_alternate_coordinate_attribute_is_found(installed_starling) -> None:
    installed_starling(attribute="xyz")
    result = StarlingEngineFactory().generate(
        request_for(30, target=30.0), None, np.random.default_rng(8)
    )
    assert result.ca_coords.shape == (1, 30, 3)


def test_bare_ensemble_without_a_mapping_is_accepted(installed_starling) -> None:
    installed_starling(as_dict=False)
    result = StarlingEngineFactory().generate(
        request_for(30, target=30.0), None, np.random.default_rng(8)
    )
    assert result.ca_coords.shape == (1, 30, 3)


def test_object_with_no_coordinates_raises_an_api_mismatch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = SimpleNamespace(generate=lambda sequence, conformations=1, **kw: object())
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    with pytest.raises(EngineError, match="could not find coordinates"):
        StarlingEngineFactory().generate(
            request_for(20, target=20.0), None, np.random.default_rng(0)
        )


def test_wrong_residue_count_raises_rather_than_being_padded(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rng = np.random.default_rng(0)
    wrong = np.stack([cone_chain(19, rng)])
    module = SimpleNamespace(
        generate=lambda sequence, conformations=1, **kw: FakeEnsemble(coordinates=wrong)
    )
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    with pytest.raises(EngineError, match="19 residues"):
        StarlingEngineFactory().generate(
            request_for(20, target=20.0), None, np.random.default_rng(0)
        )


def test_nan_from_starling_is_fatal(monkeypatch: pytest.MonkeyPatch) -> None:
    rng = np.random.default_rng(0)
    coords = np.stack([cone_chain(20, rng), cone_chain(20, rng)])
    coords[1, 5] = np.nan
    module = SimpleNamespace(
        generate=lambda sequence, conformations=1, **kw: FakeEnsemble(coordinates=coords)
    )
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    with pytest.raises(EngineError, match="non-finite"):
        StarlingEngineFactory().generate(
            request_for(20, target=20.0), None, np.random.default_rng(0)
        )


def test_starling_exceptions_are_wrapped_with_context(installed_starling) -> None:
    installed_starling(fail_with=RuntimeError("CUDA out of memory"))
    with pytest.raises(EngineError, match="CUDA out of memory") as excinfo:
        StarlingEngineFactory().generate(
            request_for(20, target=20.0), None, np.random.default_rng(0)
        )
    assert isinstance(excinfo.value.__cause__, RuntimeError)


def test_over_cap_request_points_at_the_hierarchical_engine(installed_starling) -> None:
    installed_starling()
    request = request_for(C.STARLING_MAX_LENGTH + 1, target=100.0)
    with pytest.raises(EngineError, match="HierarchicalEngine"):
        StarlingEngineFactory().generate(request, None, np.random.default_rng(0))


def test_ensemble_with_no_usable_conformer_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    # A "chain" whose bonds are 1 A: measured, rejected, and reported -- not returned.
    broken = np.zeros((3, 20, 3))
    broken[:, :, 0] = np.arange(20) * 1.0
    module = SimpleNamespace(
        generate=lambda sequence, conformations=1, **kw: FakeEnsemble(coordinates=broken)
    )
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    with pytest.raises(EngineError, match="median CA-CA distance"):
        StarlingEngineFactory().generate(
            request_for(20, target=20.0), None, np.random.default_rng(0)
        )


def test_screening_rejects_bad_geometry_and_says_why() -> None:
    rng = np.random.default_rng(41)
    good = cone_chain(30, rng, target=40.0)
    stretched = good.copy()
    stretched[15:] += np.array([6.0, 0.0, 0.0])  # one 9-ish A bond
    folded = good.copy()
    folded[10] = folded[8] + 0.5 * (folded[9] - folded[8])  # a clash and a sharp angle
    kept, notes = S.screen_conformers(np.stack([good, stretched, folded]))
    assert kept.tolist() == [0]
    assert any("bond_length" in note for note in notes)
    assert any("pseudo-angles spanned" in note for note in notes)


def test_report_summary_mentions_the_cap_and_the_residual(installed_starling) -> None:
    installed_starling()
    _, report = StarlingEngineFactory().generate_detailed(
        request_for(30, target=30.0), None, np.random.default_rng(0)
    )
    text = report.summary()
    assert "STARLING" in text and "cap" in text and "anchor residual" in text


# ---------------------------------------------------------------------------
# Segmentation arithmetic
# ---------------------------------------------------------------------------


def test_single_span_under_the_cap() -> None:
    assert H.segment_spans(100, 380, 10) == (H.SegmentSpan(0, 100),)
    assert H.segment_spans(380, 380, 10) == (H.SegmentSpan(0, 380),)


@pytest.mark.parametrize(
    ("n", "cap", "overlap"),
    [(400, 380, 10), (760, 380, 10), (1000, 380, 10), (2000, 380, 10), (121, 40, 6)],
)
def test_segment_spans_cover_the_region_and_respect_the_cap(n: int, cap: int, overlap: int) -> None:
    spans = H.segment_spans(n, cap, overlap)
    assert spans[0].start == 0
    assert spans[-1].stop == n
    assert all(len(span) <= cap for span in spans)
    for previous, following in pairwise(spans):
        assert previous.stop - following.start == overlap
    # Every residue is covered by at least one segment.
    covered = set()
    for span in spans:
        covered.update(range(span.start, span.stop))
    assert covered == set(range(n))


def test_naive_ceiling_division_would_exceed_the_cap() -> None:
    # 400 residues at a 380 cap: ceil(400/380) = 2 segments, which with 10 residues of
    # overlap must generate 410 residues -- 205 each, fine. But 770 residues would give
    # ceil(770/380) = 3 and 790/3 = 264, while the overlap-aware count is what keeps the
    # LONGEST segment under the cap for every n. Check the property directly.
    for n in range(381, 4000, 37):
        spans = H.segment_spans(n, 380, 10)
        assert max(len(span) for span in spans) <= 380


def test_segment_spans_rejects_impossible_parameters() -> None:
    with pytest.raises(EngineError, match="no progress"):
        H.segment_spans(1000, 10, 10)
    with pytest.raises(EngineError, match="at least 2"):
        H.segment_spans(1000, 100, 1)
    with pytest.raises(EngineError, match="as short as"):
        H.segment_spans(100, 25, 12)


def test_resolve_segment_cap_prefers_override_then_engine_then_constant() -> None:
    walk = ConeWalk()
    assert H.resolve_segment_cap(walk, 50) == H.LengthCap(50, "caller")
    assert H.resolve_segment_cap(walk).value == C.STARLING_MAX_LENGTH

    class Capped(ConeWalk):
        def max_length(self) -> int:
            return 123

    assert H.resolve_segment_cap(Capped()).value == 123
    assert H.resolve_segment_cap(Capped()).source.endswith("max_length")


# ---------------------------------------------------------------------------
# Hierarchical assembly over a real sub-engine
# ---------------------------------------------------------------------------


def assembled(
    n_residues: int = 121,
    *,
    cap: int = 40,
    overlap: int = 6,
    target: float | None = None,
    seed: int = 0,
    n_anchor: np.ndarray | None = None,
    c_anchor: np.ndarray | None = None,
    obstacles: np.ndarray | None = None,
    n_conformations: int = 1,
) -> tuple[IDRResult, H.AssemblyReport, H.HierarchicalEngine]:
    """Assemble a long region from cone-walk segments, at test-sized caps."""
    engine = H.HierarchicalEngine(
        ConeWalk(), max_segment_length=cap, splice_overlap=overlap, max_attempts=16
    )
    if target is None:
        target = C.flory_end_to_end(n_residues)
    request = request_for(
        n_residues,
        target=target,
        n_anchor=n_anchor,
        c_anchor=c_anchor,
        n_conformations=n_conformations,
    )
    result, report = engine.generate_detailed(request, obstacles, np.random.default_rng(seed))
    return result, report, engine


def test_assembly_splits_and_reports_which_cap_it_used() -> None:
    result, report, _ = assembled()
    assert not report.delegated
    assert report.n_segments == len(report.spans) >= 3
    assert report.cap.source == "caller"
    assert report.sub_engine == "cone-walk"
    assert result.engine == "hierarchical(cone-walk)"
    assert "assembled from" in report.notes[0]


def test_delegates_without_assembling_when_under_the_cap() -> None:
    engine = H.HierarchicalEngine(ConeWalk(), max_segment_length=200)
    request = request_for(50, target=C.flory_end_to_end(50))
    result, report = engine.generate_detailed(request, None, np.random.default_rng(0))
    assert report.delegated
    assert report.junctions == ()
    assert result.engine == "cone-walk"
    assert any("no segmentation" in note for note in report.notes)


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_assembled_chain_is_physically_valid_everywhere(seed: int) -> None:
    result, report, _ = assembled(seed=seed)
    chain = result.ca_coords[0]
    trace = validate_ca_trace(chain)
    assert trace.ok, trace.describe()
    bonds = ca_bond_lengths(chain)
    angles = ca_pseudo_angles(chain)
    assert abs(bonds - C.CA_CA_BOND_LENGTH).max() <= C.CA_CA_BOND_TOLERANCE
    assert angles.min() >= C.BACKBONE_ANGLE_MIN - 1e-6
    assert angles.max() <= C.BACKBONE_ANGLE_MAX + 1e-6
    # Junctions in particular, measured against the same window.
    assert report.junctions
    for junction in report.junctions:
        assert abs(junction.bond_length - C.CA_CA_BOND_LENGTH) <= C.CA_CA_BOND_TOLERANCE
        assert junction.angle_at_splice >= C.BACKBONE_ANGLE_MIN - 1e-6
        assert junction.angle_at_splice <= C.BACKBONE_ANGLE_MAX + 1e-6


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_assembled_chain_has_no_clashes(seed: int) -> None:
    from scipy.spatial import cKDTree

    chain = assembled(seed=seed)[0].ca_coords[0]
    tree = cKDTree(chain)
    pairs = tree.query_pairs(r=C.CA_CLASH_DISTANCE, output_type="ndarray")
    non_bonded = [
        (int(i), int(j)) for i, j in pairs if abs(int(i) - int(j)) > C.CLASH_EXCLUDE_WITHIN_RESIDUES
    ]
    assert non_bonded == []


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4, 5])
def test_achieved_end_to_end_matches_the_full_length_target(seed: int) -> None:
    result, report, _ = assembled(seed=seed)
    achieved = float(np.linalg.norm(result.ca_coords[0, -1] - result.ca_coords[0, 0]))
    tolerance = max(
        H.END_TO_END_ABSOLUTE_TOLERANCE,
        H.END_TO_END_RELATIVE_TOLERANCE * report.target_end_to_end,
    )
    assert abs(achieved - report.target_end_to_end) <= tolerance
    # The report must agree with the coordinates it describes.
    assert report.achieved_end_to_end[0] == pytest.approx(achieved, abs=1e-9)


@pytest.mark.parametrize("factor", [0.6, 1.0, 1.4])
def test_a_range_of_targets_is_hit(factor: float) -> None:
    n = 121
    target = C.flory_end_to_end(n) * factor
    result, _, _ = assembled(n, target=target, seed=7)
    achieved = float(np.linalg.norm(result.ca_coords[0, -1] - result.ca_coords[0, 0]))
    tolerance = max(H.END_TO_END_ABSOLUTE_TOLERANCE, H.END_TO_END_RELATIVE_TOLERANCE * target)
    assert abs(achieved - target) <= tolerance


def test_assembly_between_anchors_meets_both_of_them() -> None:
    anchor_n = np.zeros(3)
    anchor_c = np.array([70.0, 0.0, 0.0])
    result, report, _ = assembled(
        121, target=C.flory_end_to_end(121), n_anchor=anchor_n, c_anchor=anchor_c, seed=2
    )
    assert bool(result.success[0])
    chain = result.ca_coords[0]
    assert float(np.linalg.norm(chain[0] - anchor_n)) == pytest.approx(
        C.CA_CA_BOND_LENGTH, abs=1e-9
    )
    assert float(np.linalg.norm(chain[-1] - anchor_c)) == pytest.approx(
        C.CA_CA_BOND_LENGTH, abs=1e-9
    )
    assert validate_ca_trace(chain).ok
    # 70 A apart, the anchors permit an internal span of 62.4-77.6 A. The prediction (76.0 A)
    # is inside that, but only 1.6 A from the upper edge -- less than the engine's own 3.8 A
    # tolerance -- so the aim point is pulled in to leave that much room on both sides.
    # Aiming closer to the edge than the tolerance would accept arrangements that overshoot
    # the window and cannot be placed on the anchors at all.
    window_high = 70.0 + 2 * C.CA_CA_BOND_LENGTH
    slack = max(H.END_TO_END_ABSOLUTE_TOLERANCE, H.END_TO_END_RELATIVE_TOLERANCE * 76.0)
    assert report.target_end_to_end == pytest.approx(window_high - slack, abs=0.2)
    assert report.target_end_to_end < C.flory_end_to_end(121)


def test_anchors_override_a_target_they_cannot_accommodate() -> None:
    # 40 A apart, the anchors cap the internal span at 47.6 A, well under the 76 A
    # prediction. The anchors are real coordinates and the prediction is not, so the anchors
    # win -- and the report says so instead of quietly missing the target.
    result, report, _ = assembled(
        121,
        target=C.flory_end_to_end(121),
        n_anchor=np.zeros(3),
        c_anchor=np.array([40.0, 0.0, 0.0]),
        seed=4,
    )
    window_high = 40.0 + 2 * C.CA_CA_BOND_LENGTH
    slack = max(
        H.END_TO_END_ABSOLUTE_TOLERANCE,
        H.END_TO_END_RELATIVE_TOLERANCE * C.flory_end_to_end(121),
    )
    assert report.target_end_to_end == pytest.approx(window_high - slack, abs=1e-6)
    assert any("anchors permit" in note for note in report.notes)
    if bool(result.success[0]):
        chain = result.ca_coords[0]
        assert float(np.linalg.norm(chain[0] - np.zeros(3))) == pytest.approx(
            C.CA_CA_BOND_LENGTH, abs=1e-9
        )


def test_assembly_avoids_obstacles() -> None:
    from scipy.spatial import cKDTree

    slab = np.stack(
        np.meshgrid(np.linspace(-40.0, 40.0, 20), np.linspace(-40.0, 40.0, 20), np.array([25.0])),
        axis=-1,
    ).reshape(-1, 3)
    result, _, _ = assembled(121, obstacles=slab, seed=3)
    chain = result.ca_coords[0]
    assert bool(result.success[0])
    assert float(cKDTree(slab).query(chain, k=1)[0].min()) >= C.CA_CLASH_DISTANCE


def test_assembly_is_reproducible_for_one_seed() -> None:
    first = assembled(seed=12)[0]
    second = assembled(seed=12)[0]
    np.testing.assert_array_equal(first.ca_coords, second.ca_coords)
    third = assembled(seed=13)[0]
    assert not np.array_equal(first.ca_coords, third.ca_coords)


def test_multiple_conformations_assemble_independently() -> None:
    result, report, _ = assembled(n_conformations=3, seed=5)
    assert result.ca_coords.shape == (3, 121, 3)
    assert result.success.all()
    assert len({j.conformation for j in report.junctions}) == 3
    spans = np.linalg.norm(result.ca_coords[:, -1] - result.ca_coords[:, 0], axis=1)
    assert len(np.unique(np.round(spans, 6))) == 3
    for chain in result.ca_coords:
        assert validate_ca_trace(chain).ok


def test_junction_count_is_one_fewer_than_the_segment_count() -> None:
    _, report, _ = assembled()
    assert len(report.junctions) == report.n_segments - 1
    splices = [j.splice_residue for j in report.junctions]
    assert splices == sorted(splices)
    # Every splice lies strictly inside its overlap, with two residues of context.
    for junction in report.junctions:
        low, high = junction.overlap
        assert low + 2 <= junction.splice_residue <= high


def test_scoring_rotation_agrees_with_the_audited_one() -> None:
    # _spin_points is only used to score candidates cheaply; if it disagreed with
    # transforms.rotation_from_axis_angle the chosen spin would not be the scored one.
    from dodo.geometry.transforms import rotation_from_axis_angle

    rng = np.random.default_rng(0)
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    vector = rng.normal(size=3) * 10.0
    angles = np.linspace(0.0, 2.0 * np.pi, 17)
    scored = H._spin_points(vector, axis, angles)
    for angle, point in zip(angles, scored, strict=True):
        expected = rotation_from_axis_angle(axis, float(angle)) @ vector
        np.testing.assert_allclose(point, expected, atol=1e-9)


def test_unreachable_target_raises_with_both_numbers() -> None:
    # Ask for an end-to-end distance beyond the contour length of the whole region: no
    # arrangement of any segments can reach it.
    engine = H.HierarchicalEngine(
        ConeWalk(), max_segment_length=40, splice_overlap=6, max_attempts=4
    )
    n = 121
    request = request_for(n, target=C.contour_length(n) * 1.5)
    with pytest.raises(UnsatisfiableTargetError) as excinfo:
        engine.generate(request, None, np.random.default_rng(0))
    assert excinfo.value.target is not None and excinfo.value.achievable is not None
    assert excinfo.value.achievable < excinfo.value.target


def test_exhausted_attempts_rather_than_a_wrong_dimension() -> None:
    # A tolerance of zero cannot be met by a discrete spin search, so the engine must give
    # up loudly instead of returning a chain of the wrong size.
    engine = H.HierarchicalEngine(
        ConeWalk(),
        max_segment_length=40,
        splice_overlap=6,
        max_attempts=2,
        end_to_end_tolerance=1e-9,
    )
    request = request_for(121, target=C.flory_end_to_end(121))
    with pytest.raises(ExhaustedAttemptsError) as excinfo:
        engine.generate(request, None, np.random.default_rng(0))
    assert excinfo.value.attempts == 2
    assert "wrong dimensions" in str(excinfo.value)


def test_sub_engine_failure_is_not_silently_used() -> None:
    class FailingWalk(ConeWalk):
        name = "failing"

        def generate(
            self, request: IDRRequest, obstacles: np.ndarray | None, rng: np.random.Generator
        ) -> IDRResult:
            result = super().generate(request, obstacles, rng)
            return IDRResult.from_batch(
                ca_coords=result.ca_coords,
                success=np.zeros(request.n_conformations, dtype=bool),
                engine=self.name,
                attempts=result.attempts,
            )

    engine = H.HierarchicalEngine(FailingWalk(), max_segment_length=40, splice_overlap=6)
    with pytest.raises(BuildError, match="reported failure"):
        engine.generate(request_for(121, target=60.0), None, np.random.default_rng(0))


def test_sub_engine_returning_the_wrong_shape_is_caught() -> None:
    class ShortWalk(ConeWalk):
        name = "short"

        def generate(
            self, request: IDRRequest, obstacles: np.ndarray | None, rng: np.random.Generator
        ) -> IDRResult:
            return IDRResult.from_batch(
                ca_coords=np.zeros((request.n_conformations, request.n_residues - 1, 3)),
                success=np.ones(request.n_conformations, dtype=bool),
                engine=self.name,
                attempts=1,
            )

    engine = H.HierarchicalEngine(ShortWalk(), max_segment_length=40, splice_overlap=6)
    with pytest.raises(BuildError, match="shape"):
        engine.generate(request_for(121, target=60.0), None, np.random.default_rng(0))


def test_loose_sub_engine_geometry_is_attributed_not_hidden() -> None:
    # A sub-engine whose segments have exact bonds but pseudo-angles below DODO's 91-150
    # generation window (80-88 degrees, inside the MEASURED observed range of 75-179).
    # Assembly must proceed, hold the junctions to the sub-engine's own standard, and report
    # the violations against it by name rather than smoothing them over.
    class LooseWalk(ConeWalk):
        name = "loose"

        def generate(
            self, request: IDRRequest, obstacles: np.ndarray | None, rng: np.random.Generator
        ) -> IDRResult:
            coords = np.stack(
                [
                    cone_chain(request.n_residues, rng, angle_range=(80.0, 88.0))
                    for _ in range(request.n_conformations)
                ]
            )
            return IDRResult.from_batch(
                ca_coords=coords,
                success=np.ones(request.n_conformations, dtype=bool),
                engine=self.name,
                attempts=1,
            )

    engine = H.HierarchicalEngine(
        LooseWalk(), max_segment_length=40, splice_overlap=6, max_attempts=12
    )
    n = 121
    # Measured reach of four such segments is around 110 A, so 60 A is comfortably inside
    # what an arrangement of them can make.
    request = request_for(n, target=60.0)
    result, report = engine.generate_detailed(request, None, np.random.default_rng(0))
    assert result.success.all()
    assert report.segment_violations
    assert any("outside DODO's generation window" in note for note in report.notes)
    chain = result.ca_coords[0]
    # The junctions themselves are still measured, and against the sub-engine's own standard
    # -- exact bonds and angles inside the observed range -- they hold.
    for junction in report.junctions:
        assert abs(junction.bond_length - C.CA_CA_BOND_LENGTH) <= C.CA_CA_BOND_TOLERANCE
        assert junction.angle_at_splice >= C.BACKBONE_ANGLE_OBSERVED_MIN
        assert junction.angle_at_splice <= C.BACKBONE_ANGLE_OBSERVED_MAX
    # And the loose geometry is preserved, not quietly corrected into the window.
    assert ca_pseudo_angles(chain).min() < C.BACKBONE_ANGLE_MIN


def test_anchors_beyond_reach_raise_with_both_numbers() -> None:
    # Anchors 1000 A apart with 121 residues: the anchors force a target no arrangement of
    # segments could span. Must name both numbers rather than returning a stretched chain.
    engine = H.HierarchicalEngine(
        ConeWalk(), max_segment_length=40, splice_overlap=6, max_attempts=6
    )
    request = request_for(
        121, target=60.0, n_anchor=np.zeros(3), c_anchor=np.array([1000.0, 0.0, 0.0])
    )
    with pytest.raises(UnsatisfiableTargetError) as excinfo:
        engine.generate(request, None, np.random.default_rng(0))
    assert excinfo.value.target is not None and excinfo.value.achievable is not None
    assert excinfo.value.achievable < excinfo.value.target


def test_tiny_regions_produce_real_coordinates() -> None:
    # A 1-3 residue region is under any cap, so it delegates -- but it must still come back
    # as finite coordinates with an honest success mask, not an empty array.
    engine = H.HierarchicalEngine(ConeWalk(), max_segment_length=40, splice_overlap=6)
    for n_residues in (1, 2, 3):
        # 6.5 A is inside what 3 residues can span (2 * 3.8 * sin(theta / 2) for theta in
        # the angle window is 5.4-7.5 A); asking for less would be geometrically impossible.
        result = engine.generate(
            request_for(n_residues, target=6.5), None, np.random.default_rng(0)
        )
        assert result.ca_coords.shape == (1, n_residues, 3)
        assert result.success.all()
        assert np.all(np.isfinite(result.ca_coords))


def test_production_configuration_assembles_a_long_idr() -> None:
    # The real thing: a 900-residue IDR at STARLING's actual 380-residue cap and the real
    # splice overlap, over the real walk engine. This is the case the module exists for.
    walk = pytest.importorskip("dodo.engines.walk")
    engine = H.HierarchicalEngine(
        walk.SelfAvoidingWalk(),
        max_segment_length=C.STARLING_MAX_LENGTH,
        splice_overlap=C.SEGMENT_SPLICE_OVERLAP,
        max_attempts=12,
    )
    n = 900
    target = C.flory_end_to_end(n)
    result, report = engine.generate_detailed(
        request_for(n, target=target), None, np.random.default_rng(0)
    )
    assert report.n_segments == 3
    assert all(len(span) <= C.STARLING_MAX_LENGTH for span in report.spans)
    chain = result.ca_coords[0]
    trace = validate_ca_trace(chain)
    assert trace.ok, trace.describe()
    achieved = float(np.linalg.norm(chain[-1] - chain[0]))
    tolerance = max(H.END_TO_END_ABSOLUTE_TOLERANCE, H.END_TO_END_RELATIVE_TOLERANCE * target)
    assert abs(achieved - target) <= tolerance
    # No clash anywhere in 900 residues, junctions included.
    from scipy.spatial import cKDTree

    pairs = cKDTree(chain).query_pairs(r=C.CA_CLASH_DISTANCE, output_type="ndarray")
    assert [
        (int(i), int(j)) for i, j in pairs if abs(int(i) - int(j)) > C.CLASH_EXCLUDE_WITHIN_RESIDUES
    ] == []


def test_bad_sub_engine_object_is_rejected_at_construction() -> None:
    with pytest.raises(EngineError, match="generate"):
        H.HierarchicalEngine(object())  # type: ignore[arg-type]


def test_rng_must_be_a_generator() -> None:
    engine = H.HierarchicalEngine(ConeWalk(), max_segment_length=40, splice_overlap=6)
    with pytest.raises(TypeError, match="default_rng"):
        engine.generate(request_for(121, target=60.0), None, 7)  # type: ignore[arg-type]


def test_no_placeable_conformation_raises_rather_than_returning_nan() -> None:
    # base.py's contract: total failure raises. Box the region in so placement cannot
    # succeed anywhere, and check the engine says so instead of handing back a result object
    # whose every row is NaN.
    axis = np.linspace(-60.0, 60.0, 26)
    box = np.stack(np.meshgrid(axis, axis, axis), axis=-1).reshape(-1, 3)
    engine = H.HierarchicalEngine(
        ConeWalk(), max_segment_length=40, splice_overlap=6, max_attempts=8
    )
    request = request_for(
        121,
        target=C.flory_end_to_end(121),
        n_anchor=np.zeros(3),
        c_anchor=np.array([40.0, 0.0, 0.0]),
    )
    with pytest.raises(ExhaustedAttemptsError, match="obstacles cleared=False"):
        engine.generate_detailed(request, box, np.random.default_rng(0))


def test_a_partially_successful_batch_nan_fills_only_the_failures(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # IDRResult's contract for the in-between case: the failed conformation carries nothing
    # that reads as coordinates, and the successful one is untouched. Forced through the
    # placement seam because a partial failure is not something obstacles can be aimed at.
    real = H.place_between_anchors
    calls = {"n": 0}

    def flaky(*args: Any, **kwargs: Any) -> S.AnchorPlacement:
        placement = real(*args, **kwargs)
        calls["n"] += 1
        if calls["n"] == 1:
            # Same coordinates, declared unplaceable: only the verdict changes.
            return S.AnchorPlacement(
                ca_coords=placement.ca_coords,
                conformer_index=placement.conformer_index,
                achieved_end_to_end=placement.achieved_end_to_end,
                desired_end_to_end=placement.desired_end_to_end,
                n_anchor_gap=placement.n_anchor_gap,
                c_anchor_gap=placement.c_anchor_gap,
                anchor_residual=10.0,
                min_internal_ca_distance=placement.min_internal_ca_distance,
                clash_free=False,
                relaxed_to=None,
                orientations_tried=placement.orientations_tried,
                notes=("forced failure",),
            )
        return placement

    monkeypatch.setattr(H, "place_between_anchors", flaky)
    result, report = assembled(121, n_conformations=2, seed=7)[:2]
    assert result.success.tolist() == [False, True]
    assert np.all(np.isnan(result.ca_coords[0]))
    assert np.all(np.isfinite(result.ca_coords[1]))
    # The real geometry survives on the placement, for diagnosis rather than for writing.
    assert np.all(np.isfinite(report.placements[0].ca_coords))
    assert validate_ca_trace(report.placements[0].ca_coords).ok


def test_hierarchical_over_starling_segments(installed_starling) -> None:
    # The configuration that matters in production: STARLING under hierarchical assembly.
    # Still a fake STARLING, but the whole path -- ensemble, screening, selection,
    # splicing, validation -- runs.
    installed_starling(n_conformers=48)
    sub = StarlingEngineFactory(ensemble_size=48, oversample=1)
    engine = H.HierarchicalEngine(sub, max_segment_length=40, splice_overlap=6, max_attempts=16)
    n = 121
    result, report = engine.generate_detailed(
        request_for(n, target=C.flory_end_to_end(n)), None, np.random.default_rng(0)
    )
    assert report.sub_engine == "starling"
    assert report.n_segments >= 3
    assert result.ca_coords.shape == (1, n, 3)
    assert validate_ca_trace(result.ca_coords[0]).ok
    for junction in report.junctions:
        assert abs(junction.bond_length - C.CA_CA_BOND_LENGTH) <= C.CA_CA_BOND_TOLERANCE


def test_cap_comes_from_the_sub_engine_when_not_overridden(
    monkeypatch: pytest.MonkeyPatch, installed_starling
) -> None:
    installed_starling(max_length=256)
    engine = H.HierarchicalEngine(StarlingEngineFactory())
    cap = engine.cap()
    assert cap.value == 256
    assert "starling" in cap.source


def test_default_sub_engine_falls_back_to_the_walk() -> None:
    walk = pytest.importorskip("dodo.engines.walk")
    engine = H.default_sub_engine()
    assert isinstance(engine, walk.SelfAvoidingWalk)


def test_assembly_over_the_real_walk_engine_if_present() -> None:
    walk = pytest.importorskip("dodo.engines.walk")
    engine = H.HierarchicalEngine(
        walk.SelfAvoidingWalk(), max_segment_length=60, splice_overlap=8, max_attempts=16
    )
    n = 200
    result, report = engine.generate_detailed(
        request_for(n, target=C.flory_end_to_end(n)), None, np.random.default_rng(0)
    )
    chain = result.ca_coords[0]
    assert validate_ca_trace(chain).ok, validate_ca_trace(chain).describe()
    tolerance = max(
        H.END_TO_END_ABSOLUTE_TOLERANCE,
        H.END_TO_END_RELATIVE_TOLERANCE * report.target_end_to_end,
    )
    achieved = float(np.linalg.norm(chain[-1] - chain[0]))
    assert abs(achieved - report.target_end_to_end) <= tolerance


def test_report_summaries_are_informative() -> None:
    _, report, _ = assembled()
    text = report.summary()
    assert "segments" in text and "junction" in text and "target" in text
    assert str(report.junctions[0]).startswith("junction ")
    assert str(report.spans[0]).startswith("[")


def test_placement_rejects_bad_obstacle_shapes() -> None:
    rng = np.random.default_rng(0)
    chain = cone_chain(10, rng)
    with pytest.raises(GeometryError, match="obstacles"):
        S.place_between_anchors(
            chain,
            n_anchor_xyz=None,
            c_anchor_xyz=None,
            rng=rng,
            obstacles=np.zeros((4, 2)),
        )


def test_placement_rejects_degenerate_conformer_arrays() -> None:
    rng = np.random.default_rng(0)
    with pytest.raises(GeometryError):
        S.place_between_anchors(np.zeros((0, 3)), n_anchor_xyz=None, c_anchor_xyz=None, rng=rng)
    with pytest.raises(GeometryError):
        S.rank_conformers(np.zeros((2, 4)), 10.0)


# ---------------------------------------------------------------------------
# The remaining guarded assumptions about STARLING's API, one test each
# ---------------------------------------------------------------------------


def test_entry_point_is_found_on_the_submodule_path(monkeypatch: pytest.MonkeyPatch) -> None:
    # STARLING re-exports generate() from starling.frontend.ensemble_generation; a release
    # that stops re-exporting it at top level must still work.
    inner = fake_starling()
    module = SimpleNamespace(
        frontend=SimpleNamespace(ensemble_generation=SimpleNamespace(generate=inner.generate))
    )
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    result = StarlingEngineFactory().generate(
        request_for(24, target=C.flory_end_to_end(24)), None, np.random.default_rng(0)
    )
    assert result.ca_coords.shape == (1, 24, 3)


def test_max_length_imports_starling_when_not_already_loaded(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = fake_starling(max_length=411)
    monkeypatch.setattr(S, "_loaded_starling", lambda: None)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    assert S.starling_max_length().value == 411
    # probe=False must not import anything, so it falls back to the constant.
    assert S.starling_max_length(probe=False).value == C.STARLING_MAX_LENGTH


def test_trajectory_xyz_in_nanometres_is_understood(monkeypatch: pytest.MonkeyPatch) -> None:
    # The other route to coordinates: an mdtraj trajectory, whose units are nanometres.
    rng = np.random.default_rng(0)
    chains = [cone_chain(20, rng) for _ in range(3)]
    stack = np.stack(chains) / 10.0
    payload = SimpleNamespace(trajectory=SimpleNamespace(xyz=stack))
    module = SimpleNamespace(generate=lambda sequence, conformations=1, **kw: payload)
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    # These three chains are whatever dimension a free cone walk gave them, so the request
    # asks for the dimension one of them has, in Angstroms: this test is about interpreting
    # the units of STARLING's return value, not about hitting a target. A literal target
    # cannot do that -- 20.0 A was reachable for a 20-residue walk sampling angles up to 161
    # degrees and is not reachable now that the window stops at 150 and the same walk comes
    # back more compact.
    target = span_of(chains[0])
    result, report = StarlingEngineFactory().generate_detailed(
        request_for(20, target=target), None, np.random.default_rng(0)
    )
    assert any("nanometres" in note for note in report.notes)
    assert abs(ca_bond_lengths(result.ca_coords[0]) - C.CA_CA_BOND_LENGTH).max() < 1e-9


def test_bare_array_and_two_dimensional_return_are_accepted(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rng = np.random.default_rng(0)
    single = cone_chain(20, rng)
    module = SimpleNamespace(generate=lambda sequence, conformations=1, **kw: single)
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    # The fixture is one specific chain, so the request must ask for the dimension it has:
    # this test is about interpreting STARLING's return value, not about hitting a target.
    result = StarlingEngineFactory().generate(
        request_for(20, target=span_of(single)), None, np.random.default_rng(0)
    )
    assert result.ca_coords.shape == (1, 20, 3)


def test_ambiguous_mapping_raises_rather_than_picking_one(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rng = np.random.default_rng(0)
    stack = np.stack([cone_chain(20, rng)])
    module = SimpleNamespace(
        generate=lambda sequence, conformations=1, **kw: {
            "other": FakeEnsemble(stack),
            "another": FakeEnsemble(stack),
        }
    )
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    with pytest.raises(EngineError, match="which ensemble"):
        StarlingEngineFactory().generate(
            request_for(20, target=20.0), None, np.random.default_rng(0)
        )


def test_mapping_with_one_unrecognized_key_is_accepted(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rng = np.random.default_rng(0)
    stack = np.stack([cone_chain(20, rng)])
    module = SimpleNamespace(
        generate=lambda sequence, conformations=1, **kw: {"protein_1": FakeEnsemble(stack)}
    )
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    result = StarlingEngineFactory().generate(
        request_for(20, target=span_of(stack[0])), None, np.random.default_rng(0)
    )
    assert result.ca_coords.shape == (1, 20, 3)


def test_empty_ensemble_and_distance_map_shapes_raise(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    for payload, match in (
        (np.zeros((0, 20, 3)), "zero conformers"),
        (np.zeros((20, 20)), "distance map"),
        (np.zeros((2, 20, 4)), "expects"),
    ):
        module = SimpleNamespace(
            generate=lambda sequence, conformations=1, _p=payload, **kw: FakeEnsemble(_p)
        )
        monkeypatch.setattr(S, "_import_starling", lambda _m=module: _m)
        monkeypatch.setattr(S, "starling_installed", lambda: True)
        monkeypatch.setattr(S, "_loaded_starling", lambda _m=module: _m)
        with pytest.raises(EngineError, match=match):
            StarlingEngineFactory().generate(
                request_for(20, target=20.0), None, np.random.default_rng(0)
            )


def test_a_coordinate_accessor_that_needs_arguments_is_skipped(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rng = np.random.default_rng(0)
    stack = np.stack([cone_chain(20, rng)])

    class Awkward:
        def coordinates(self, frame: int) -> np.ndarray:  # needs an argument
            return stack[frame]

        xyz = stack

    module = SimpleNamespace(generate=lambda sequence, conformations=1, **kw: Awkward())
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    result = StarlingEngineFactory().generate(
        request_for(20, target=span_of(stack[0])), None, np.random.default_rng(0)
    )
    assert result.ca_coords.shape == (1, 20, 3)


def test_generate_without_a_conformations_parameter_is_an_api_mismatch(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = SimpleNamespace(generate=lambda sequence: None)
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    with pytest.raises(EngineUnavailableError, match="conformations"):
        StarlingEngineFactory().generate(
            request_for(20, target=20.0), None, np.random.default_rng(0)
        )


def test_device_is_passed_through(installed_starling) -> None:
    record: dict[str, Any] = {}
    installed_starling(record=record)
    StarlingEngineFactory(device="cpu").generate(
        request_for(20, target=C.flory_end_to_end(20)), None, np.random.default_rng(0)
    )
    assert record["kwargs"].get("device") == "cpu"


def test_screening_can_be_disabled_and_says_so(installed_starling) -> None:
    installed_starling()
    _, report = StarlingEngineFactory(screen=False).generate_detailed(
        request_for(20, target=C.flory_end_to_end(20)), None, np.random.default_rng(0)
    )
    assert any("Screening disabled" in note for note in report.notes)


def test_an_ensemble_of_stretched_conformers_is_rejected_wholesale(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Median bond stays 3.8 so the unit check passes, but every conformer has one 6 A bond.
    rng = np.random.default_rng(0)
    stack = np.stack([cone_chain(30, rng) for _ in range(4)])
    stack[:, 20:] += np.array([6.0, 0.0, 0.0])
    module = SimpleNamespace(generate=lambda sequence, conformations=1, **kw: FakeEnsemble(stack))
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    with pytest.raises(EngineError, match="none of them survived"):
        StarlingEngineFactory().generate(
            request_for(30, target=30.0), None, np.random.default_rng(0)
        )


def test_screening_rejects_an_internal_clash() -> None:
    # A hairpin with exact bonds whose two strands are 1.5 A apart: no rigid-body placement
    # can fix that, so the conformer is unusable rather than merely imperfect. The angle
    # criterion is switched off here so that the clash criterion is what fires.
    n = 20
    chain = np.zeros((n, 3))
    for i in range(n):
        if i < n // 2:
            chain[i] = [3.8 * i, 0.0, 0.0]
        else:
            chain[i] = [3.8 * (n - 1 - i), 1.5, 0.0]
    kept, notes = S.screen_conformers(chain[None, :, :], angle_window=None, bond_tolerance=4.0)
    assert kept.size == 0
    assert any("internal_clash" in note for note in notes)


def test_engine_configuration_is_validated() -> None:
    with pytest.raises(EngineError, match="ensemble_size"):
        S.StarlingEngine(ensemble_size=0)
    with pytest.raises(EngineError, match="oversample"):
        S.StarlingEngine(oversample=0)
    with pytest.raises(EngineError, match="spin_samples"):
        H.HierarchicalEngine(ConeWalk(), spin_samples=0)
    with pytest.raises(EngineError, match="max_attempts"):
        H.HierarchicalEngine(ConeWalk(), max_attempts=0)
    with pytest.raises(EngineError, match="max_segment_length"):
        H.HierarchicalEngine(ConeWalk(), max_segment_length=0).cap()
    with pytest.raises(EngineError, match="Cannot split"):
        H.segment_spans(0, 40, 6)


def test_repr_and_availability_pass_through() -> None:
    engine = H.HierarchicalEngine(ConeWalk(), max_segment_length=40)
    assert engine.available() is True
    assert "cone-walk" in repr(engine)
    assert "available=False" in repr(S.StarlingEngine())


def test_default_sub_engine_prefers_starling_when_available(installed_starling) -> None:
    installed_starling()
    assert H.default_sub_engine().name == "starling"


def test_delegated_report_summary_says_no_assembly_was_needed() -> None:
    engine = H.HierarchicalEngine(ConeWalk(), max_segment_length=200)
    _, report = engine.generate_detailed(
        request_for(50, target=C.flory_end_to_end(50)), None, np.random.default_rng(0)
    )
    assert "no assembly needed" in report.summary()


def test_coincident_termini_do_not_produce_nan() -> None:
    # A conformer whose first and last CA coincide has no end-to-end direction. The old code
    # divided by that norm and emitted NaN coordinates; the guard must return real geometry.
    rng = np.random.default_rng(0)
    chain = cone_chain(20, rng)
    chain[-1] = chain[0]
    placement = S.place_between_anchors(
        chain,
        n_anchor_xyz=np.zeros(3),
        c_anchor_xyz=np.array([5.0, 0.0, 0.0]),
        rng=rng,
    )
    assert np.all(np.isfinite(placement.ca_coords))


def test_single_residue_placement_handles_every_anchor_configuration() -> None:
    rng = np.random.default_rng(0)
    single = np.zeros((1, 3))
    for n_anchor, c_anchor in (
        (np.zeros(3), None),
        (None, np.zeros(3)),
        (None, None),
        (np.zeros(3), np.zeros(3)),
    ):
        placement = S.place_between_anchors(
            single, n_anchor_xyz=n_anchor, c_anchor_xyz=c_anchor, rng=rng
        )
        assert np.all(np.isfinite(placement.ca_coords))
        if n_anchor is not None and c_anchor is not None and np.allclose(n_anchor, c_anchor):
            # Coincident anchors cannot both be one bond away from one residue; the
            # construction puts it one bond from the midpoint and the residual says so.
            assert placement.n_anchor_gap == pytest.approx(C.CA_CA_BOND_LENGTH, abs=1e-9)


def test_empty_obstacle_array_is_not_a_clash() -> None:
    rng = np.random.default_rng(0)
    chain = cone_chain(10, rng)
    placement = S.place_between_anchors(
        chain,
        n_anchor_xyz=None,
        c_anchor_xyz=None,
        rng=rng,
        obstacles=np.zeros((0, 3)),
    )
    assert placement.clash_free


def test_junction_validation_raises_on_a_fabricated_violation() -> None:
    # The guard that says "splicing cannot produce invalid geometry": prove it fires by
    # handing it a chain that has been broken at a junction residue.
    engine = H.HierarchicalEngine(ConeWalk(), max_segment_length=40, splice_overlap=6)
    rng = np.random.default_rng(0)
    chain = cone_chain(30, rng)
    chain[15:] += np.array([9.0, 0.0, 0.0])  # a 12 A bond at residue 14-15
    junction = H.Junction(
        conformation=0,
        index=1,
        splice_residue=15,
        overlap=(10, 16),
        spin_degrees=0.0,
        bond_length=12.0,
        angle_at_splice=120.0,
        clash_cutoff=C.CA_CLASH_DISTANCE,
    )
    with pytest.raises(BuildError, match="bug in hierarchical assembly"):
        engine._validate_assembly(chain, [junction], strict=True)


def test_segment_target_and_helpers_handle_edges() -> None:
    engine = H.HierarchicalEngine(ConeWalk(), max_segment_length=40)
    assert engine._segment_target(100.0, 1, 100) == 0.0
    # A segment target is capped by the segment's own contour length.
    assert engine._segment_target(10_000.0, 10, 100) <= 0.95 * C.contour_length(10)
    assert H._achieved_end_to_end(np.zeros((5, 3)))[0] == 0.0
    assert H._achieved_end_to_end(np.zeros((2, 1, 3))) == (0.0, 0.0)
    with pytest.raises(GeometryError, match="Expected"):
        H._achieved_end_to_end(np.zeros((4, 4)))
    assert np.isnan(H._angle_at(np.zeros((3, 3)), 0))
    assert np.isnan(H._angle_at(np.zeros((3, 3)), 2))


def test_segments_that_cannot_fold_report_exhaustion() -> None:
    # Near-rod segments asked for a compact full-length target: every junction leaves the
    # target out of reach, so the arrangement dead-ends and the engine says so instead of
    # returning a rod and calling it compact.
    class RodWalk(ConeWalk):
        name = "rod"

        def generate(
            self, request: IDRRequest, obstacles: np.ndarray | None, rng: np.random.Generator
        ) -> IDRResult:
            coords = np.stack(
                [
                    cone_chain(request.n_residues, rng, angle_range=(178.0, 179.0))
                    for _ in range(request.n_conformations)
                ]
            )
            return IDRResult.from_batch(
                ca_coords=coords,
                success=np.ones(request.n_conformations, dtype=bool),
                engine=self.name,
                attempts=1,
            )

    engine = H.HierarchicalEngine(
        RodWalk(), max_segment_length=40, splice_overlap=6, max_attempts=5
    )
    with pytest.raises(ExhaustedAttemptsError, match="wrong dimensions"):
        engine.generate(request_for(121, target=30.0), None, np.random.default_rng(0))


def test_a_segment_that_is_invalid_even_loosely_fails_the_conformation() -> None:
    # A sub-engine whose segments contain a 9 A virtual bond: outside even the loosened
    # standard, so the assembled chain cannot be reconstructed to all-atom. It must come back
    # as a failure with the violations named, not as a success with a broken bond in it.
    class BrokenWalk(ConeWalk):
        name = "broken"

        def generate(
            self, request: IDRRequest, obstacles: np.ndarray | None, rng: np.random.Generator
        ) -> IDRResult:
            chains = []
            for _ in range(request.n_conformations):
                chain = cone_chain(request.n_residues, rng)
                chain[request.n_residues // 2 :] += np.array([9.0, 0.0, 0.0])
                chains.append(chain)
            return IDRResult.from_batch(
                ca_coords=np.stack(chains),
                success=np.ones(request.n_conformations, dtype=bool),
                engine=self.name,
                attempts=1,
            )

    engine = H.HierarchicalEngine(
        BrokenWalk(), max_segment_length=40, splice_overlap=6, max_attempts=6
    )
    with pytest.raises(ExhaustedAttemptsError, match="not a valid CA trace") as excinfo:
        engine.generate_detailed(request_for(121, target=90.0), None, np.random.default_rng(0))
    assert "CA-CA bond" in str(excinfo.value)


def test_a_one_residue_starling_request_is_placed(installed_starling) -> None:
    # Degenerate but real: a single missing residue. No bonds and no angles exist to measure,
    # so every measurement path has to cope with an empty array rather than divide by zero.
    installed_starling(n_conformers=4)
    result, report = StarlingEngineFactory().generate_detailed(
        request_for(1, target=1.0, n_anchor=np.zeros(3)), None, np.random.default_rng(0)
    )
    assert result.ca_coords.shape == (1, 1, 3)
    assert bool(result.success[0])
    assert report.placements[0].n_anchor_gap == pytest.approx(C.CA_CA_BOND_LENGTH, abs=1e-9)


def test_non_finite_conformers_are_rejected_by_the_array_guards() -> None:
    rng = np.random.default_rng(0)
    chain = cone_chain(10, rng)
    broken = chain.copy()
    broken[4] = np.nan
    with pytest.raises(GeometryError, match="non-finite"):
        S.place_between_anchors(broken, n_anchor_xyz=None, c_anchor_xyz=None, rng=rng)
    with pytest.raises(GeometryError, match="non-finite"):
        S.rank_conformers(broken[None, :, :], 10.0)


# ---------------------------------------------------------------------------
# Physical-validity gates
#
# Every test below reproduces a defect that was measured on this module: geometry that
# cannot be reconstructed to an all-atom backbone, returned with success=True. The gate
# they all check against is dodo.geometry.metrics.validate_ca_trace, which measures
# non-bonded contacts as well as bonds and angles -- bonds and angles are purely local, so
# a chain can satisfy both perfectly while folding back through itself.
# ---------------------------------------------------------------------------


#: A CA-CA-CA pseudo-angle, in degrees, that is comfortably *inside* DODO's generation
#: window and near its blunt end -- so a chain built from it is nearly extended and
#: self-avoids without effort.
#:
#: DERIVED from :data:`~dodo.constants.BACKBONE_ANGLE_MAX` rather than fixed at a literal.
#: The helpers below use it for the geometry that is supposed to be *valid*, and a literal
#: cannot do that job: 160.0 was inside the window when the upper bound was 161.0, and when
#: the bound came down to 150.0 the "valid except for two angles" chain started reporting a
#: violation at all 116 of its blunt vertices.
IN_WINDOW_BLUNT_ANGLE = C.BACKBONE_ANGLE_MAX - 2.0


def sharp_angle_chain(
    n_residues: int = 120,
    *,
    sharp_vertices: tuple[int, ...] = (40, 80),
    sharp_angle: float = 50.0,
    blunt_angle: float = IN_WINDOW_BLUNT_ANGLE,
    bond: float = C.CA_CA_BOND_LENGTH,
    seed: int = 7,
) -> np.ndarray:
    """Grow a self-avoiding chain with exact bonds and a few impossible angles.

    Everything about this trace is valid except ``len(sharp_vertices)`` pseudo-angles, which
    sit far below even the observed 75-179 degree range. It exists to check that a *fraction*
    of impossible angles is not treated as acceptable noise.
    """
    rng = np.random.default_rng(seed)
    coords = np.zeros((n_residues, 3), dtype=np.float64)
    coords[1] = np.array([bond, 0.0, 0.0])
    for i in range(2, n_residues):
        axis = coords[i - 1] - coords[i - 2]
        axis /= np.linalg.norm(axis)
        theta = np.deg2rad(sharp_angle if (i - 1) in sharp_vertices else blunt_angle)
        for _ in range(500):
            perpendicular = np.cross(axis, rng.normal(size=3))
            norm = float(np.linalg.norm(perpendicular))
            if norm < 1e-9:  # pragma: no cover - measure-zero draw
                continue
            perpendicular /= norm
            direction = axis * np.cos(np.pi - theta) + perpendicular * np.sin(np.pi - theta)
            candidate = coords[i - 1] + bond * direction
            if i > 2 and float(np.linalg.norm(coords[: i - 2] - candidate, axis=1).min()) < (
                C.CA_CLASH_DISTANCE + 0.1
            ):
                continue
            coords[i] = candidate
            break
        else:  # pragma: no cover - 500 draws per step effectively never all fail
            raise AssertionError(f"sharp_angle_chain dead-ended at residue {i}")
    return coords


def helix_residues_per_turn(margin: float = 2.0) -> int:
    """Residues per turn whose helical pseudo-angle stays inside the generation window.

    DERIVED. A planar regular n-gon has interior angle ``180 - 360 / n``, and adding a rise
    turns it into a helix whose pseudo-angle sits a fraction of a degree *above* that, so
    :data:`~dodo.constants.BACKBONE_ANGLE_MAX` caps ``n`` and ``margin`` degrees of headroom
    absorbs the rise.

    Not a literal, for the same reason as :data:`IN_WINDOW_BLUNT_ANGLE`: 12 residues per turn
    is exactly 150.000 degrees before the rise and 150.033 after it, which was well inside
    the 161-degree window this helper was written against and is outside the 150-degree one.
    At the current window this is 11 residues per turn, or 147.3 degrees.
    """
    return int(np.floor(360.0 / (180.0 - (C.BACKBONE_ANGLE_MAX - margin))))


def touching_helix(
    n_residues: int = 40,
    *,
    residues_per_turn: int | None = None,
    contact: float = 2.11,
    bond: float = C.CA_CA_BOND_LENGTH,
) -> np.ndarray:
    """Build a helix whose turns touch: exact bonds, in-window angles, CAs ``contact`` apart.

    Closed-form so the numbers are exact rather than sampled. Residues one turn apart sit at
    the same angular position, so their separation is the rise per turn -- set directly to
    ``contact``. This is the shape of conformer that a purely local bond-and-angle screen
    cannot see: every bond is exactly ``CA_CA_BOND_LENGTH`` and every pseudo-angle is inside
    DODO's window, by :func:`helix_residues_per_turn`.
    """
    if residues_per_turn is None:
        residues_per_turn = helix_residues_per_turn()
    rise = contact / residues_per_turn
    chord = float(np.sqrt(bond**2 - rise**2))
    radius = chord / (2.0 * np.sin(np.pi / residues_per_turn))
    angles = 2.0 * np.pi * np.arange(n_residues) / residues_per_turn
    return np.stack(
        [radius * np.cos(angles), radius * np.sin(angles), rise * np.arange(n_residues)], axis=1
    )


def install_fixed_ensemble(monkeypatch: pytest.MonkeyPatch, stack: np.ndarray) -> SimpleNamespace:
    """Install a fake STARLING that always returns exactly ``stack``."""
    module = SimpleNamespace(
        generate=lambda sequence, conformations=1, **kw: FakeEnsemble(coordinates=stack)
    )
    monkeypatch.setattr(S, "_import_starling", lambda: module)
    monkeypatch.setattr(S, "starling_installed", lambda: True)
    monkeypatch.setattr(S, "_loaded_starling", lambda: module)
    return module


def test_the_sharp_angle_helper_is_valid_except_for_its_angles() -> None:
    # Premise of the next two tests: nothing but the two angles is wrong with this chain.
    chain = sharp_angle_chain()
    report = validate_ca_trace(chain)
    assert {v.kind for v in report.violations} == {"pseudo_angle"}
    assert len(report.violations) == 2
    assert ca_pseudo_angles(chain).min() == pytest.approx(50.0, abs=1e-6)


def test_screen_rejects_a_conformer_with_even_one_impossible_angle() -> None:
    # 2 of 118 vertices is 1.7%, which the old 2%-of-vertices allowance waved through. A
    # 50 degree CA-CA-CA angle puts the i-1 and i+1 CAs 3.2 A apart and leaves no room for
    # the intervening backbone at all, so one is one too many.
    chain = sharp_angle_chain()
    kept, notes = S.screen_conformers(np.stack([chain] * 4))
    assert kept.size == 0
    assert any("pseudo_angle" in note for note in notes)


def test_engine_refuses_an_ensemble_whose_conformers_all_have_sharp_angles(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    install_fixed_ensemble(monkeypatch, np.stack([sharp_angle_chain()] * 4))
    with pytest.raises(EngineError, match="none of them survived"):
        StarlingEngineFactory().generate(
            request_for(120, target=100.0), None, np.random.default_rng(0)
        )


def test_the_touching_helix_helper_is_valid_except_for_its_contacts() -> None:
    # Premise of the next three tests.
    chain = touching_helix()
    report = validate_ca_trace(chain)
    assert {v.kind for v in report.violations} == {"steric_clash"}
    assert S._min_internal_ca_distance(chain) == pytest.approx(2.11, abs=1e-6)
    assert float(np.abs(ca_bond_lengths(chain) - C.CA_CA_BOND_LENGTH).max()) < 1e-9
    angles = ca_pseudo_angles(chain)
    assert angles.min() >= C.BACKBONE_ANGLE_MIN and angles.max() <= C.BACKBONE_ANGLE_MAX


def test_screen_rejects_a_conformer_whose_own_ca_atoms_touch() -> None:
    # The internal-clash floor must be CA_CLASH_DISTANCE, not the last rung of the
    # relaxation ladder: two carbon atoms 2.11 A apart are inside each other.
    kept, notes = S.screen_conformers(np.stack([touching_helix()] * 3))
    assert kept.size == 0
    assert any("internal_clash" in note for note in notes)


def test_placement_is_not_ok_when_the_conformer_clashes_with_itself() -> None:
    # A rigid-body placement cannot change a conformer's internal geometry, so ok must
    # consult min_internal_ca_distance. Without this, a self-intersecting chain is reported
    # as a usable placement purely because it clears the external obstacles.
    placement = S.place_between_anchors(
        touching_helix(),
        n_anchor_xyz=None,
        c_anchor_xyz=None,
        rng=np.random.default_rng(0),
        obstacles=None,
    )
    assert placement.min_internal_ca_distance == pytest.approx(2.11, abs=1e-6)
    assert not placement.ok
    assert not validate_ca_trace(placement.ca_coords).ok


def test_engine_refuses_an_ensemble_of_self_touching_conformers(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    install_fixed_ensemble(monkeypatch, np.stack([touching_helix()] * 3))
    with pytest.raises(EngineError, match="none of them survived"):
        StarlingEngineFactory().generate(
            request_for(40, target=16.0), None, np.random.default_rng(0)
        )


def test_engine_reports_failure_when_the_ensemble_cannot_reach_the_target(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # A free region: no anchor constrains the dimension, so request.target_end_to_end is the
    # only thing that does. An ensemble of compact conformers answered a 350 A request with a
    # ~40 A chain and called it a success, silently discarding every constants.MODES
    # multiplier along with the target.
    rng = np.random.default_rng(3)
    stack = np.stack([cone_chain(100, rng, target=40.0) for _ in range(6)])
    spans = np.linalg.norm(stack[:, -1] - stack[:, 0], axis=1)
    assert spans.max() < 200.0  # premise
    install_fixed_ensemble(monkeypatch, stack)
    with pytest.raises(ExhaustedAttemptsError, match="end-to-end"):
        StarlingEngineFactory().generate(
            request_for(100, target=350.0), None, np.random.default_rng(0)
        )


def test_engine_accepts_a_target_the_ensemble_can_reach(monkeypatch: pytest.MonkeyPatch) -> None:
    # Positive control for the test above: the tolerance must not reject a good answer.
    rng = np.random.default_rng(4)
    stack = np.stack([cone_chain(100, rng, target=t) for t in (30.0, 60.0, 90.0, 120.0)])
    spans = np.linalg.norm(stack[:, -1] - stack[:, 0], axis=1)
    install_fixed_ensemble(monkeypatch, stack)
    target = float(spans[2])
    result, report = StarlingEngineFactory().generate_detailed(
        request_for(100, target=target), None, np.random.default_rng(0)
    )
    assert bool(result.success[0])
    achieved = float(result.end_to_end_distances[0])
    assert abs(achieved - target) <= S.end_to_end_tolerance(target)
    assert validate_ca_trace(result.ca_coords[0]).ok
    assert report.placements[0].ok


def test_engine_raises_rather_than_returning_all_nan(monkeypatch: pytest.MonkeyPatch) -> None:
    # Zero successful conformations is total failure, and base.py's contract is that total
    # failure raises. Returning IDRResult(0/2 built) with NaN rows looks like a result.
    rng = np.random.default_rng(5)
    stack = np.stack([cone_chain(30, rng, target=30.0) for _ in range(4)])
    install_fixed_ensemble(monkeypatch, stack)
    # A dense obstacle cloud around the origin: no orientation of any conformer clears it.
    grid = np.linspace(-60.0, 60.0, 25)
    obstacles = np.stack(np.meshgrid(grid, grid, grid), axis=-1).reshape(-1, 3)
    with pytest.raises(ExhaustedAttemptsError):
        StarlingEngineFactory().generate(
            request_for(30, target=30.0, n_conformations=2), obstacles, np.random.default_rng(0)
        )


# -- hierarchical assembly --------------------------------------------------


class RelaxedWalk(ConeWalk):
    """A sub-engine that builds valid segments but declares a relaxed clash threshold.

    The walk engine does exactly this when no strict-threshold step is available: it returns
    coordinates and reports the rung it had to use. Assembly must carry that disclosure onto
    the assembled result -- a segment built at 2.8 A is spliced into the chain unchanged, so
    the assembled chain is a 2.8 A chain.
    """

    name = "relaxed"

    def __init__(self, relaxed_to: float = 2.8) -> None:
        super().__init__()
        self.relaxed_to = relaxed_to

    def generate(
        self, request: IDRRequest, obstacles: np.ndarray | None, rng: np.random.Generator
    ) -> IDRResult:
        result = super().generate(request, obstacles, rng)
        return IDRResult.from_batch(
            ca_coords=result.ca_coords,
            success=result.success,
            engine=self.name,
            attempts=result.attempts,
            relaxed_to=self.relaxed_to,
        )


def test_hierarchical_propagates_a_sub_engine_relaxed_clash_threshold() -> None:
    engine = H.HierarchicalEngine(
        RelaxedWalk(2.8), max_segment_length=40, splice_overlap=6, max_attempts=16
    )
    result, report = engine.generate_detailed(
        request_for(121, target=C.flory_end_to_end(121)), None, np.random.default_rng(0)
    )
    assert result.success.all()
    # The disclosure, not None: a caller that cares about contacts must be told.
    assert result.relaxed_to == pytest.approx(2.8)
    assert any("2.8" in note for note in report.notes)


def test_hierarchical_reports_the_loosest_rung_a_segment_used() -> None:
    # Two different rungs across the segments: the threshold "actually used" for the chain as
    # a whole is the loosest one. Reporting the strictest would understate the relaxation.
    class MixedWalk(ConeWalk):
        name = "mixed"

        def __init__(self) -> None:
            super().__init__()
            self.rungs = [3.2, 2.8, 2.5, 2.8, 3.2, 2.5]

        def generate(
            self, request: IDRRequest, obstacles: np.ndarray | None, rng: np.random.Generator
        ) -> IDRResult:
            result = super().generate(request, obstacles, rng)
            rung = self.rungs[(self.calls - 1) % len(self.rungs)]
            return IDRResult.from_batch(
                ca_coords=result.ca_coords,
                success=result.success,
                engine=self.name,
                attempts=result.attempts,
                relaxed_to=None if rung >= C.CA_CLASH_DISTANCE else rung,
            )

    engine = H.HierarchicalEngine(
        MixedWalk(), max_segment_length=40, splice_overlap=6, max_attempts=16
    )
    result, _ = engine.generate_detailed(
        request_for(121, target=C.flory_end_to_end(121)), None, np.random.default_rng(1)
    )
    assert result.success.all()
    assert result.relaxed_to == pytest.approx(2.5)


def test_hierarchical_fails_a_conformation_for_an_undeclared_segment_clash() -> None:
    # A sub-engine that returns a self-touching segment and reports relaxed_to=None. The
    # contact is spliced into the assembled chain, so the chain is not physically valid and
    # the assembled result must say so rather than reporting success with a 2.11 A pair in it.
    class TouchingWalk(ConeWalk):
        name = "touching"

        def generate(
            self, request: IDRRequest, obstacles: np.ndarray | None, rng: np.random.Generator
        ) -> IDRResult:
            coords = np.stack(
                [
                    folded_chain(request.n_residues, rng, target=request.target_end_to_end)
                    for _ in range(request.n_conformations)
                ]
            )
            return IDRResult.from_batch(
                ca_coords=coords,
                success=np.ones(request.n_conformations, dtype=bool),
                engine=self.name,
                attempts=1,
            )

    engine = H.HierarchicalEngine(
        TouchingWalk(), max_segment_length=40, splice_overlap=6, max_attempts=4
    )
    with pytest.raises(BuildError, match="clash distance") as excinfo:
        engine.generate(request_for(121, target=60.0), None, np.random.default_rng(0))
    assert "not a valid CA trace" in str(excinfo.value)


def test_assembled_chain_is_checked_for_whole_chain_clashes() -> None:
    # _validate_assembly must measure non-bonded contacts, and must attribute them to the
    # sub-engine rather than raising as if the splice had produced them.
    engine = H.HierarchicalEngine(ConeWalk(), max_segment_length=40, splice_overlap=6)
    chain = touching_helix(60)
    report = engine._validate_assembly(chain, (), strict=True, clash_distance=C.CA_CLASH_DISTANCE)
    assert not report.ok
    assert report.of_kind("steric_clash")


def test_hierarchical_raises_rather_than_returning_all_nan() -> None:
    grid = np.linspace(-80.0, 80.0, 25)
    obstacles = np.stack(np.meshgrid(grid, grid, grid), axis=-1).reshape(-1, 3)
    engine = H.HierarchicalEngine(
        ConeWalk(), max_segment_length=40, splice_overlap=6, max_attempts=6
    )
    with pytest.raises(ExhaustedAttemptsError):
        engine.generate(
            request_for(121, target=C.flory_end_to_end(121)),
            obstacles,
            np.random.default_rng(0),
        )


def test_desired_internal_span_leaves_slack_inside_the_anchor_window() -> None:
    anchor_n = np.zeros(3)
    anchor_c = np.array([100.0, 0.0, 0.0])
    low = 100.0 - 2 * C.CA_CA_BOND_LENGTH
    high = 100.0 + 2 * C.CA_CA_BOND_LENGTH
    # Without slack the aim point sits exactly on the window edge, so any undershoot at all
    # puts both anchors out of reach and the conformation is thrown away.
    assert S.desired_internal_span(10.0, anchor_n, anchor_c) == (low, (low, high))
    aim, window = S.desired_internal_span(10.0, anchor_n, anchor_c, slack=3.0)
    assert window == (low, high)
    assert aim == pytest.approx(low + 3.0)
    aim, _ = S.desired_internal_span(1000.0, anchor_n, anchor_c, slack=3.0)
    assert aim == pytest.approx(high - 3.0)
    # A target already inside the slack-shrunk window is untouched.
    assert S.desired_internal_span(99.0, anchor_n, anchor_c, slack=3.0)[0] == pytest.approx(99.0)
    # Slack wider than half the window collapses onto the centre, which is the anchor
    # separation -- the point with the most room on both sides.
    aim, _ = S.desired_internal_span(10.0, anchor_n, anchor_c, slack=100.0)
    assert aim == pytest.approx(100.0)


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4, 5])
def test_interior_assembly_is_usable_for_a_prediction_far_from_the_anchors(seed: int) -> None:
    # The module's primary case: a long IDR between two folded domains whose separation
    # disagrees with the predicted Re. The clipped target used to land exactly on the window
    # edge, so a fraction of an Angstrom of undershoot made both anchors unreachable and the
    # conformation came back as an all-NaN failure -- 35% of interior builds.
    separation = 150.0
    result, report, _engine = assembled(
        121,
        target=C.flory_end_to_end(121),
        n_anchor=np.zeros(3),
        c_anchor=np.array([separation, 0.0, 0.0]),
        seed=seed,
    )
    assert bool(result.success[0]), report.notes
    chain = result.ca_coords[0]
    assert validate_ca_trace(chain).ok
    assert float(np.linalg.norm(chain[0])) == pytest.approx(C.CA_CA_BOND_LENGTH, abs=1e-6)
    assert float(np.linalg.norm(chain[-1] - np.array([separation, 0.0, 0.0]))) == pytest.approx(
        C.CA_CA_BOND_LENGTH, abs=1e-6
    )


def test_segment_targets_follow_the_clipped_target_not_the_prediction() -> None:
    # The arrangement aims at the anchor-clipped target, so the segments must be sized from
    # the same number. Sized from the unclipped prediction they come back many times too
    # extended -- or, for a prediction far above the anchors, too compact to bridge them --
    # and the arrangement cannot fold them to fit.
    seen: list[float] = []

    class RecordingWalk(ConeWalk):
        name = "recording"

        def generate(
            self, request: IDRRequest, obstacles: np.ndarray | None, rng: np.random.Generator
        ) -> IDRResult:
            seen.append(request.target_end_to_end)
            return super().generate(request, obstacles, rng)

    engine = H.HierarchicalEngine(
        RecordingWalk(), max_segment_length=40, splice_overlap=6, max_attempts=8
    )
    n = 121
    request = request_for(
        n,
        target=20.0,  # far below what 20 A-separated anchors... see below
        n_anchor=np.zeros(3),
        c_anchor=np.array([150.0, 0.0, 0.0]),
    )
    # The per-segment schedule is what is under test, not whether this particular build
    # closes, so a BuildError from the arrangement is not a failure of this test.
    with contextlib.suppress(BuildError):
        engine.generate_detailed(request, None, np.random.default_rng(0))
    assert seen, "the sub-engine was never called"
    clipped = 150.0 - 2 * C.CA_CA_BOND_LENGTH
    # Every per-segment target must be consistent with a full-length target inside the
    # anchor window, not with the 20 A prediction.
    smallest_segment = min(len(span) for span in H.segment_spans(n, 40, 6))
    floor = clipped * (smallest_segment / n) ** C.FLORY_RE_EXPONENT * 0.5
    assert min(seen) >= floor, (min(seen), floor)


def folded_chain(
    n_residues: int,
    rng: np.random.Generator,
    *,
    target: float,
    contact: tuple[float, float] = (2.2, 3.15),
) -> np.ndarray:
    """Fold a valid chain back on itself until two non-bonded CAs are ``contact`` apart.

    Built by taking a clean cone chain and spinning its second half about the axis through
    residues ``m-2`` and ``m-1``. That rotation is an isometry of the tail and leaves the
    pivot bond and the angle at ``m-1`` untouched, so **every bond and every pseudo-angle in
    the result is exactly what the generator produced** -- the only thing wrong with it is a
    non-bonded contact. That is the shape of defect a local bond-and-angle check cannot see,
    and the reason validate_ca_trace has to measure contacts.
    """
    for _ in range(20):
        chain = cone_chain(n_residues, rng, target=target)
        splice = n_residues // 2
        pivot = chain[splice - 1]
        axis_vector = chain[splice - 2] - pivot
        axis = axis_vector / float(np.linalg.norm(axis_vector))
        best: np.ndarray | None = None
        best_error = np.inf
        for angle in np.linspace(0.0, 2.0 * np.pi, 721, endpoint=False):
            rotation = rotation_from_axis_angle(axis, float(angle))
            tail = apply(chain[splice:] - pivot, rotation) + pivot
            candidate = np.vstack((chain[:splice], tail))
            closest = S._min_internal_ca_distance(candidate)
            if not contact[0] <= closest <= contact[1]:
                continue
            # The spin also moves the far end, so keep the one that leaves the chain closest
            # to the dimension it was asked for: a segment folded into a ball is rejected by
            # the arrangement for being unable to reach anything, which would stop the test
            # before it got to the check it is about.
            error = abs(span_of(candidate) - target)
            if error < best_error:
                best, best_error = candidate, error
        if best is not None and best_error <= 0.2 * target:
            return best
    raise AssertionError(  # pragma: no cover - 20 chains x 721 spins always hit one
        f"folded_chain could not bring two CAs of a {n_residues}-residue chain to within "
        f"{contact} A while spanning about {target:.1f} A"
    )


def test_a_delegated_all_nan_result_is_not_passed_on() -> None:
    # A sub-engine that breaks its own contract by returning IDRResult(0/1 built) instead of
    # raising. HierarchicalEngine is the object the caller holds, so it does not forward an
    # all-NaN result that reads as a success.
    class GivesUp(ConeWalk):
        name = "gives-up"

        def generate(
            self, request: IDRRequest, obstacles: np.ndarray | None, rng: np.random.Generator
        ) -> IDRResult:
            return IDRResult.from_batch(
                ca_coords=np.zeros((request.n_conformations, request.n_residues, 3)),
                success=np.zeros(request.n_conformations, dtype=bool),
                engine=self.name,
                attempts=3,
            )

    engine = H.HierarchicalEngine(GivesUp(), max_segment_length=200, splice_overlap=6)
    with pytest.raises(ExhaustedAttemptsError, match="all-NaN"):
        engine.generate(
            request_for(50, target=C.flory_end_to_end(50)), None, np.random.default_rng(0)
        )


def test_the_new_gates_are_bit_identical_across_equal_seeds() -> None:
    # Every gate added here reads the rng (screening order, spins, orientations), so identical
    # seeds must still give identical coordinates -- and different seeds must not.
    first, _, _ = assembled(
        121,
        target=C.flory_end_to_end(121),
        n_anchor=np.zeros(3),
        c_anchor=np.array([150.0, 0.0, 0.0]),
        seed=4,
    )
    second, _, _ = assembled(
        121,
        target=C.flory_end_to_end(121),
        n_anchor=np.zeros(3),
        c_anchor=np.array([150.0, 0.0, 0.0]),
        seed=4,
    )
    np.testing.assert_array_equal(first.ca_coords, second.ca_coords)
    assert first.relaxed_to == second.relaxed_to


def test_a_declared_relaxation_is_honoured_and_a_tighter_contact_is_not() -> None:
    """Pin the one policy decision: a success below CA_CLASH_DISTANCE must be *declared*.

    ``relaxed_to`` exists in base.py precisely so that a chain built at a rung of
    CLASH_RELAXATION_LADDER can be returned with that fact attached, and the walk engine
    already does exactly this. So an assembled chain is validated at the threshold the
    result will disclose -- and only at that threshold. A contact tighter than anything the
    sub-engine declared is not covered by the disclosure and fails the conformation.
    """

    def walk_declaring(rung: float, contact: tuple[float, float]) -> ConeWalk:
        class Declaring(ConeWalk):
            name = "declaring"

            def generate(
                self, request: IDRRequest, obstacles: np.ndarray | None, rng: np.random.Generator
            ) -> IDRResult:
                coords = np.stack(
                    [
                        folded_chain(
                            request.n_residues,
                            rng,
                            target=request.target_end_to_end,
                            contact=contact,
                        )
                        for _ in range(request.n_conformations)
                    ]
                )
                return IDRResult.from_batch(
                    ca_coords=coords,
                    success=np.ones(request.n_conformations, dtype=bool),
                    engine=self.name,
                    attempts=1,
                    relaxed_to=rung,
                )

        return Declaring()

    # Contact inside the declared rung: accepted, and the disclosure is on the result.
    engine = H.HierarchicalEngine(
        walk_declaring(2.8, (2.95, 3.15)),
        max_segment_length=40,
        splice_overlap=6,
        max_attempts=6,
    )
    result, _ = engine.generate_detailed(
        request_for(121, target=60.0), None, np.random.default_rng(0)
    )
    chain = result.ca_coords[0]
    assert result.success.all()
    assert result.relaxed_to == pytest.approx(2.8)
    closest = S._min_internal_ca_distance(chain)
    assert 2.8 <= closest < C.CA_CLASH_DISTANCE
    # Valid at the number the result reports, which is the whole point of reporting it.
    assert validate_ca_trace(chain, clash_distance=result.relaxed_to).ok
    # And knowingly not valid at the strict default, which relaxed_to is what discloses.
    assert not validate_ca_trace(chain).ok

    # Contact tighter than the declared rung: not covered by the disclosure, so it fails.
    engine = H.HierarchicalEngine(
        walk_declaring(2.8, (2.2, 2.5)),
        max_segment_length=40,
        splice_overlap=6,
        max_attempts=4,
    )
    with pytest.raises(BuildError, match="clash distance"):
        engine.generate(request_for(121, target=60.0), None, np.random.default_rng(0))


def stretched_chain(n_residues: int = 40, *, bond: float = 4.25, seed: int = 11) -> np.ndarray:
    """Scale a chain up so every CA-CA bond is ``bond`` A, leaving the angles in window.

    Scaling is a similarity transform, so the pseudo-angles are untouched and the non-bonded
    distances scale up -- there is nothing wrong with this trace except that every one of its
    virtual bonds is longer than a trans peptide can make.
    """
    chain = cone_chain(n_residues, np.random.default_rng(seed))
    return chain * (bond / C.CA_CA_BOND_LENGTH)


def test_screen_rejects_bonds_longer_than_a_peptide_can_make() -> None:
    # constants.py: the trans-peptide CA-CA virtual bond is 3.80-3.81 A and "remarkably
    # rigid". A 4.25 A virtual bond therefore has no all-atom realization at all, and a
    # symmetric 0.5 A screen accepted every one of them (3.30-4.30 A).
    chain = stretched_chain()
    assert float(np.abs(ca_bond_lengths(chain) - 4.25).max()) < 1e-9  # premise
    angles = ca_pseudo_angles(chain)
    assert angles.min() >= C.BACKBONE_ANGLE_MIN and angles.max() <= C.BACKBONE_ANGLE_MAX
    assert S._min_internal_ca_distance(chain) >= C.CA_CLASH_DISTANCE
    kept, notes = S.screen_conformers(np.stack([chain] * 3))
    assert kept.size == 0
    assert any("bond_length" in note for note in notes)
    assert any("longer than" in note for note in notes)


def test_assembly_does_not_loosen_the_bond_gate_for_a_sub_engine() -> None:
    # The loose standard exists for pseudo-angles, which real structures populate over
    # 75-179 degrees. It must not extend to bonds: a 4.25 A CA-CA bond is unreconstructable
    # whoever produced it, so an assembled chain full of them is not a success.
    class StretchedWalk(ConeWalk):
        name = "stretched"

        def generate(
            self, request: IDRRequest, obstacles: np.ndarray | None, rng: np.random.Generator
        ) -> IDRResult:
            coords = np.stack(
                [
                    cone_chain(request.n_residues, rng) * (4.25 / C.CA_CA_BOND_LENGTH)
                    for _ in range(request.n_conformations)
                ]
            )
            return IDRResult.from_batch(
                ca_coords=coords,
                success=np.ones(request.n_conformations, dtype=bool),
                engine=self.name,
                attempts=1,
            )

    engine = H.HierarchicalEngine(
        StretchedWalk(), max_segment_length=40, splice_overlap=6, max_attempts=4
    )
    with pytest.raises(BuildError, match="CA-CA bond"):
        engine.generate(request_for(121, target=90.0), None, np.random.default_rng(0))
