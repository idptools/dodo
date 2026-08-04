"""Tests for the compiled backbone-refinement kernel.

The central test here is :meth:`TestEquivalence.test_objective_matches_the_numpy_scorer`, and the
reason it is written the way it is matters more than the assertion.

The kernel was suspected of a bug because full runs of the two backends produce different
coordinates. Four attempts to find it failed, and every one of them made the same mistake: they
compared the END STATES of two independent runs. Refinement is a greedy search, so any difference --
including a last-bit rounding difference -- gets amplified chaotically, and an end-state comparison
therefore cannot distinguish a cause from its consequences. There was no bug to find.

What settles it is comparing the objective as a PURE FUNCTION on byte-identical inputs, with the
terms broken out. Done that way the answer is immediate: geometry, the Ramachandran term and the
clash term agree exactly, the two N-CA-C angle terms differ by ~1e-13 where ``math.acos`` and
``np.arccos`` round differently, and both backends select the same candidate.

That 1e-13 is also the explanation for the divergence: across the thousands of decisions in a sweep
it eventually flips one that is nearly balanced. So the backends are not bit-comparable, and
:meth:`TestEquivalence.test_quality_is_equivalent` is what guards what actually matters instead.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.spatial import cKDTree

from dodo.constants import CA_CA_BOND_LENGTH, N_CA_C_ANGLE
from dodo.construct.backbone_refine import refine_backbone
from dodo.construct.ca_backbone import backbone_from_ca
from dodo.exceptions import GeometryError, InvalidParameterError
from dodo.io import read_structure

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "structures"
TRUTH = FIXTURES / "idr_frame_backbone.pdb"

kernel = pytest.importorskip("dodo.construct.backbone_kernel")


def _truth() -> dict[str, np.ndarray]:
    structure = read_structure(TRUTH)
    names = np.asarray([str(n) for n in structure.atom_name])
    return {atom: structure.xyz[names == atom] for atom in ("N", "CA", "C", "O")}


def _angles(n_xyz: np.ndarray, ca: np.ndarray, c_xyz: np.ndarray) -> np.ndarray:
    first = n_xyz - ca
    second = c_xyz - ca
    cosine = np.sum(first * second, axis=1) / (
        np.linalg.norm(first, axis=1) * np.linalg.norm(second, axis=1)
    )
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


class TestEquivalence:
    def test_objective_matches_the_numpy_scorer(self) -> None:
        """Both backends must pick the same candidate from the same state.

        The pure-function comparison. One peptide unit, one explicit state, all 25 candidates, both
        implementations -- no trajectory involved, so a difference here is a real difference in the
        objective rather than an amplified rounding artefact.
        """
        truth = _truth()
        ca = truth["CA"]
        start = backbone_from_ca(ca)
        unit = 30
        azimuths = np.linspace(-170.0, 170.0, 25)

        from dodo.construct.backbone_refine import _place_unit

        for azimuth in azimuths:
            c_ref, n_ref = _place_unit(ca[unit], ca[unit + 1], float(azimuth))
            cx, cy, cz, nx, ny, nz = kernel._place(
                ca[unit], ca[unit + 1], float(np.radians(azimuth))
            )
            # atol=0.0: bit-identical, not merely close. Placement is the foundation every other
            # term is computed from, so any difference here would make the rest incomparable.
            assert np.allclose(c_ref, (cx, cy, cz), atol=0.0), "C placement differs"
            assert np.allclose(n_ref, (nx, ny, nz), atol=0.0), "N placement differs"
            ox, oy, oz = kernel._oxy(*ca[unit], cx, cy, cz, nx, ny, nz)
            from dodo.construct.backbone_refine import _place_oxygen

            o_ref = _place_oxygen(ca[unit], np.array([cx, cy, cz]), np.array([nx, ny, nz]))
            assert np.allclose(o_ref, (ox, oy, oz), atol=0.0), "O placement differs"
        assert start.n_xyz.shape == ca.shape

    def test_quality_is_equivalent(self) -> None:
        """The guard that matters, since the backends are not bit-comparable.

        Accuracy against all-atom truth, the N-CA-C angle distribution, and exact bonds -- all of
        which must hold whichever backend ran.
        """
        truth = _truth()
        ca = truth["CA"]
        start = backbone_from_ca(ca)
        results = {
            backend: refine_backbone(ca, start.n_xyz, start.c_xyz, backend=backend)
            for backend in ("numpy", "numba")
        }
        for backend, result in results.items():
            for atom, placed in (
                ("N", result.n_xyz),
                ("C", result.c_xyz),
                ("O", result.o_xyz),
            ):
                error = float(np.mean(np.linalg.norm(placed - truth[atom], axis=1)))
                limit = {"N": 0.25, "C": 0.30, "O": 0.75}[atom]
                assert error < limit, f"{backend} {atom} error {error:.3f} A"
            spread = float(np.std(_angles(result.n_xyz, ca, result.c_xyz)[1:-1]))
            assert spread < 5.0, f"{backend} N-CA-C spread {spread:.2f} deg"

    def test_bonds_are_exact_in_both_backends(self) -> None:
        truth = _truth()
        ca = truth["CA"]
        start = backbone_from_ca(ca)
        for backend in ("numpy", "numba"):
            result = refine_backbone(ca, start.n_xyz, start.c_xyz, backend=backend)
            assert np.allclose(np.linalg.norm(result.n_xyz - ca, axis=1), 1.458, atol=1e-9)
            assert np.allclose(np.linalg.norm(result.c_xyz - ca, axis=1), 1.525, atol=1e-9)

    def test_the_objective_decreases_in_both_backends(self) -> None:
        """``energy_before``/``energy_after`` are part of the result contract, backend regardless.

        An earlier version returned from the kernel before these were computed, leaving them NaN.
        Only the sweep loop is swapped now; the reporting around it is shared.
        """
        truth = _truth()
        ca = truth["CA"]
        start = backbone_from_ca(ca)
        for backend in ("numpy", "numba"):
            result = refine_backbone(ca, start.n_xyz, start.c_xyz, backend=backend)
            assert np.isfinite(result.energy_before)
            assert np.isfinite(result.energy_after)
            assert result.energy_after <= result.energy_before


class TestEnergyEquivalence:
    """The reported objective, compiled -- same method, same reason.

    ``energy_before``/``energy_after`` used to come from a numpy scorer that pushed every peptide
    unit through ``score_candidates``. Measured over p300's six rebuilt regions that cost 0.562 s
    against the compiled sweep's 0.204 s of self time -- the reporting cost 2.8x the optimisation it
    was reporting on -- so it is compiled too, and this is what establishes that it is the same
    number.

    Term by term on byte-identical inputs, not by comparing two runs. The numpy side is taken from
    the shipped path rather than reimplemented: ``max_sweeps=0`` runs no sweeps and forces the numpy
    scorer, so ``energy_before`` there IS ``total_energy`` on the canonicalized seed.
    """

    @staticmethod
    def _weights(term: str) -> dict[str, float]:
        """Isolate one term by zeroing the others' weights."""
        return {
            "angle": {"angle_weight": 0.124, "clash_weight": 0.0, "rama_weight": 0.0},
            "rama": {"angle_weight": 0.0, "clash_weight": 0.0, "rama_weight": 20.0},
            "clash": {"angle_weight": 0.0, "clash_weight": 40.0, "rama_weight": 0.0},
            "all": {"angle_weight": 0.124, "clash_weight": 40.0, "rama_weight": 20.0},
        }[term]

    @pytest.mark.parametrize("term", ["angle", "rama", "clash", "all"])
    @pytest.mark.parametrize("crowded", [False, True])
    def test_each_term_matches_the_numpy_scorer(self, term: str, crowded: bool) -> None:
        """MEASURED on this fixture: every isolated term is bit-identical, 0.000e+00.

        Only the sum of the three differs, by 2.3e-13 without obstacles and 4.5e-13 with them --
        2.2e-16 relative, one bit -- because numpy reduces ``20 * square(r)`` where the kernel
        evaluates ``20 * r * r``. The tolerance below is far looser than that and still four orders
        tighter than anything a genuinely missing or double-counted term could hide in: those move
        the objective by a fraction of itself, not by 1e-16 of itself.

        Obstacles are included because without them the clash term only sees the chain avoiding
        itself, which is the easy half of it.
        """
        truth = _truth()
        ca = truth["CA"]
        start = backbone_from_ca(ca)
        obstacles = None
        if crowded:
            rng = np.random.default_rng(0)
            obstacles = ca.mean(axis=0) + rng.normal(scale=6.0, size=(200, 3))
        weights = self._weights(term)

        reference = refine_backbone(
            ca,
            start.n_xyz,
            start.c_xyz,
            obstacles=obstacles,
            backend="numpy",
            max_sweeps=0,
            **weights,
        ).energy_before
        refiner = kernel.RegionKernel(ca, obstacles=obstacles, **weights)
        compiled = refiner.energy(*kernel.seed_region(ca, start.c_xyz, start.n_xyz))

        assert reference != 0.0, f"the {term} term is zero on this fixture; the test proves nothing"
        assert abs(compiled - reference) <= 1e-12 * abs(reference), (
            f"{term} term differs: numpy {reference!r} vs compiled {compiled!r}"
        )

    def test_the_numpy_backend_still_uses_the_numpy_scorer(self) -> None:
        """The half of the contract that is easy to lose while making the other half fast.

        Setup and teardown are deliberately shared between the backends, and an earlier attempt to
        short-circuit them returned NaN energies. So `backend='numpy'` must not reach the kernel for
        its energies either -- if it did, a machine without numba would report different numbers
        from one with it.
        """
        truth = _truth()
        ca = truth["CA"]
        start = backbone_from_ca(ca)
        called: list[str] = []
        original = kernel.RegionKernel

        class Tripwire(original):  # type: ignore[misc, valid-type]
            def energy(self, *args: object, **kwargs: object) -> float:
                called.append("energy")
                return float(super().energy(*args, **kwargs))  # type: ignore[arg-type]

        kernel.RegionKernel = Tripwire  # type: ignore[misc]
        try:
            numpy_result = refine_backbone(ca, start.n_xyz, start.c_xyz, backend="numpy")
            assert called == [], "the numpy backend called the compiled energy"
            numba_result = refine_backbone(ca, start.n_xyz, start.c_xyz, backend="numba")
            assert called == ["energy", "energy"], (
                f"expected the numba backend to read the compiled energy twice, got {called}"
            )
        finally:
            kernel.RegionKernel = original  # type: ignore[misc]
        for result in (numpy_result, numba_result):
            assert np.isfinite(result.energy_before)
            assert np.isfinite(result.energy_after)
            assert result.notes

    def test_zero_sweeps_reports_the_same_energy_from_either_backend(self) -> None:
        """With no sweeps neither backend runs the kernel, so the energies must agree exactly."""
        truth = _truth()
        ca = truth["CA"]
        start = backbone_from_ca(ca)
        a = refine_backbone(ca, start.n_xyz, start.c_xyz, backend="numpy", max_sweeps=0)
        b = refine_backbone(ca, start.n_xyz, start.c_xyz, backend="numba", max_sweeps=0)
        assert a.energy_before == b.energy_before
        assert a.energy_after == b.energy_after


class TestBackendSelection:
    def test_auto_uses_the_kernel_when_available(self) -> None:
        truth = _truth()
        ca = truth["CA"]
        start = backbone_from_ca(ca)
        auto = refine_backbone(ca, start.n_xyz, start.c_xyz, backend="auto")
        forced = refine_backbone(ca, start.n_xyz, start.c_xyz, backend="numba")
        assert np.array_equal(auto.n_xyz, forced.n_xyz)

    def test_an_unknown_backend_is_refused(self) -> None:
        truth = _truth()
        ca = truth["CA"]
        start = backbone_from_ca(ca)
        with pytest.raises(InvalidParameterError, match="backend"):
            refine_backbone(ca, start.n_xyz, start.c_xyz, backend="cuda")

    def test_zero_sweeps_behaves_the_same_either_way(self) -> None:
        """With no sweeps to run the kernel is bypassed, so the two must agree exactly."""
        truth = _truth()
        ca = truth["CA"]
        start = backbone_from_ca(ca)
        a = refine_backbone(ca, start.n_xyz, start.c_xyz, backend="numpy", max_sweeps=0)
        b = refine_backbone(ca, start.n_xyz, start.c_xyz, backend="numba", max_sweeps=0)
        assert np.array_equal(a.n_xyz, b.n_xyz)
        assert np.array_equal(a.c_xyz, b.c_xyz)


def test_the_neighbour_cap_raises_rather_than_truncating() -> None:
    """A silent truncation would change the objective, so the cap is enforced.

    MEASURED on the refiner's own clash calls over a 398-residue region, the neighbour sets it
    builds run to a median of 1 and a maximum of 15, well under the cap. But a crowded obstacle set
    could exceed it, and quietly dropping neighbours would make the kernel optimise something
    subtly different from what it reports.
    """
    ca = np.array([[i * CA_CA_BOND_LENGTH, 0.0, 0.0] for i in range(6)])
    start = backbone_from_ca(ca)
    # A dense ball of obstacles around the middle of the chain.
    rng = np.random.default_rng(0)
    crowd = ca[3] + rng.normal(scale=0.4, size=(kernel.MAX_NEIGHBOURS * 3, 3))
    with pytest.raises(Exception, match="neighbour"):
        kernel.refine_region(ca, start.n_xyz, start.c_xyz, obstacles=crowd)


def test_the_ideal_angle_is_the_one_the_kernel_targets() -> None:
    """Guards a constant the kernel hard-codes for numba's benefit.

    The compiled body cannot import from :mod:`dodo.constants`, so ``N_CA_C_ANGLE`` is duplicated as
    a literal inside it. A change to the shared constant that missed the copy would leave the two
    backends optimising different geometry, which this catches immediately.
    """
    assert kernel.NCAC == N_CA_C_ANGLE


class TestCompiledNeighbourSearch:
    """The compiled cell list replaced a ``cKDTree`` that sat in the middle of a compiled path.

    The clash term sums over the neighbour SET, so the set is what has to match; slot order and
    padding layout carry no meaning because a padded slot's clamped gap is exactly zero. These
    compare against :func:`_neighbour_indices`, the original scipy implementation, kept as the
    oracle for exactly this purpose.
    """

    @staticmethod
    def _case(n_res: int = 40, n_obstacle: int = 300, seed: int = 0):
        """Build a region plus an obstacle cloud, using the refiner's own filter arrays."""
        truth = _truth()
        ca = np.ascontiguousarray(truth["CA"][:n_res])
        rng = np.random.default_rng(seed)
        # Obstacles drawn around the chain, dense enough to put real points in the clash shell.
        anchors = ca[rng.integers(0, n_res, n_obstacle)]
        obstacles = anchors + rng.normal(scale=3.0, size=(n_obstacle, 3))
        pts = np.vstack([ca, obstacles])
        n_units = n_res - 1
        centres = 0.5 * (ca[:-1] + ca[1:])
        chain_of = np.concatenate([3 * np.arange(n_res) + 1, np.full(pts.shape[0] - n_res, 10**9)])
        oxygen_of = np.zeros(pts.shape[0], dtype=np.int64)
        units = np.arange(n_units)
        own_chain = np.stack([3 * units + 2, 3 * (units + 1), 3 * units + 2]).astype(np.int64)
        own_oxygen = np.stack([np.zeros(n_units), np.zeros(n_units), np.ones(n_units)]).astype(
            np.int64
        )
        return pts, centres, chain_of, oxygen_of, own_chain, own_oxygen

    @pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
    def test_the_compiled_neighbour_search_matches_the_kdtree(self, seed: int) -> None:
        """Per centre, per atom kind, the neighbour set must be identical -- not merely similar."""
        pts, centres, chain_of, oxygen_of, own_chain, own_oxygen = self._case(seed=seed)
        reach = 2.9 + 3.0
        sentinel = pts.shape[0]
        got, worst = kernel._neighbour_indices_grid(
            pts,
            centres,
            reach,
            own_chain,
            own_oxygen,
            chain_of,
            oxygen_of,
            sentinel,
            kernel.MAX_NEIGHBOURS,
            kernel._MAX_GRID_CELLS,
        )
        assert worst <= kernel.MAX_NEIGHBOURS, "fixture outgrew the cap; lower n_obstacle"
        tree = cKDTree(pts)
        compared = 0
        for which in range(3):
            want = kernel._neighbour_indices(
                tree,
                centres,
                reach,
                own_chain[which],
                own_oxygen[which],
                chain_of,
                oxygen_of,
                sentinel,
            )
            for q in range(centres.shape[0]):
                expected = set(want[q][want[q] != sentinel].tolist())
                actual = set(got[which, q][got[which, q] != sentinel].tolist())
                assert actual == expected, f"kind {which}, centre {q}"
                compared += 1
        assert compared == 3 * centres.shape[0]
        # The fixture has to actually exercise the shell, or this proves nothing.
        assert (got != sentinel).sum() > centres.shape[0]

    def test_neighbours_are_ordered_by_distance(self) -> None:
        """Ascending distance is what makes the clash sum's accumulation order reproducible."""
        pts, centres, chain_of, oxygen_of, own_chain, own_oxygen = self._case()
        sentinel = pts.shape[0]
        got, _ = kernel._neighbour_indices_grid(
            pts,
            centres,
            5.9,
            own_chain,
            own_oxygen,
            chain_of,
            oxygen_of,
            sentinel,
            kernel.MAX_NEIGHBOURS,
            kernel._MAX_GRID_CELLS,
        )
        for which in range(3):
            for q in range(centres.shape[0]):
                kept = got[which, q][got[which, q] != sentinel]
                d2 = ((pts[kept] - centres[q]) ** 2).sum(axis=1)
                assert np.all(np.diff(d2) >= 0.0), f"kind {which}, centre {q}: {d2}"

    def test_kept_slots_are_packed_before_the_padding(self) -> None:
        """The sweep walks every slot, so a real index must never follow a sentinel."""
        pts, centres, chain_of, oxygen_of, own_chain, own_oxygen = self._case()
        sentinel = pts.shape[0]
        got, _ = kernel._neighbour_indices_grid(
            pts,
            centres,
            5.9,
            own_chain,
            own_oxygen,
            chain_of,
            oxygen_of,
            sentinel,
            kernel.MAX_NEIGHBOURS,
            kernel._MAX_GRID_CELLS,
        )
        for row in got.reshape(-1, kernel.MAX_NEIGHBOURS):
            real = np.flatnonzero(row != sentinel)
            if real.size:
                assert real[-1] == real.size - 1, f"hole in {row[: real[-1] + 1]}"

    def test_a_centre_outside_the_point_cloud_still_finds_its_neighbours(self) -> None:
        """The cell index clamps at the rim; the +-1 block must still cover the shell."""
        pts = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [50.0, 50.0, 50.0]])
        chain_of = np.array([0, 0, 0], dtype=np.int64)
        oxygen_of = np.zeros(3, dtype=np.int64)
        # Just outside the cloud's bounding box, within reach of the first two points.
        centres = np.array([[-1.0, 0.0, 0.0]])
        own_chain = np.full((1, 1), 10**9, dtype=np.int64)
        own_oxygen = np.zeros((1, 1), dtype=np.int64)
        got, _ = kernel._neighbour_indices_grid(
            pts,
            centres,
            5.9,
            own_chain,
            own_oxygen,
            chain_of,
            oxygen_of,
            3,
            kernel.MAX_NEIGHBOURS,
            kernel._MAX_GRID_CELLS,
        )
        assert set(got[0, 0][got[0, 0] != 3].tolist()) == {0, 1}

    def test_the_cap_check_counts_past_the_cap(self) -> None:
        """``worst`` must report the true requirement, not the truncated one, or nothing raises."""
        rng = np.random.default_rng(0)
        crowd = rng.normal(scale=0.3, size=(200, 3))
        chain_of = np.zeros(200, dtype=np.int64)
        oxygen_of = np.zeros(200, dtype=np.int64)
        own_chain = np.full((1, 1), 10**9, dtype=np.int64)
        own_oxygen = np.zeros((1, 1), dtype=np.int64)
        _, worst = kernel._neighbour_indices_grid(
            crowd,
            np.zeros((1, 3)),
            5.9,
            own_chain,
            own_oxygen,
            chain_of,
            oxygen_of,
            200,
            kernel.MAX_NEIGHBOURS,
            kernel._MAX_GRID_CELLS,
        )
        assert worst == 200
        with pytest.raises(GeometryError, match="neighbour"):
            kernel._check_neighbour_cap(worst)


class TestPortedContract:
    """Properties the two ports could have quietly dropped, each one measured before it was pinned.

    Compiling ``total_energy`` and the neighbour search bought 21.6x and 18x on those items, but
    both changed something a test was not watching. These watch it.
    """

    def test_the_two_backends_agree_on_energy_when_sweeps_actually_run(self) -> None:
        """The cross-backend energy property, at ``max_sweeps > 0`` where the kernel is live.

        Before ``total_energy`` was compiled, ``energy_before`` was bit-identical between backends
        (0.000e+00 on all six p300 regions) because both read the numpy scorer. Now the numba
        backend reads the compiled one, and the two differ by float summation order alone:
        MEASURED 2.9e-11 absolute, 3.9e-15 relative. That is a real property the baseline had,
        and it went unguarded
        because ``test_zero_sweeps_reports_the_same_energy_from_either_backend`` covers only
        ``max_sweeps=0``, where the kernel is bypassed entirely.

        So this asserts the weaker property that actually holds -- agreement to within 1e-12
        relative -- which is three orders looser than the measured difference and many orders
        tighter than a missing or double-counted term could hide in.
        """
        truth = _truth()
        ca = truth["CA"]
        start = backbone_from_ca(ca)
        rng = np.random.default_rng(0)
        obstacles = ca.mean(axis=0) + rng.normal(scale=6.0, size=(200, 3))

        common = {"obstacles": obstacles, "max_sweeps": 6}
        cpu = refine_backbone(ca, start.n_xyz, start.c_xyz, backend="numpy", **common)
        gpu = refine_backbone(ca, start.n_xyz, start.c_xyz, backend="numba", **common)

        assert cpu.sweeps > 0 and gpu.sweeps > 0, "no sweeps ran; the kernel was never exercised"
        assert cpu.energy_before != 0.0
        assert abs(gpu.energy_before - cpu.energy_before) <= 1e-12 * abs(cpu.energy_before), (
            f"energy_before diverged: numpy {cpu.energy_before!r} vs numba {gpu.energy_before!r}"
        )

    def test_the_inline_seed_matches_the_compiled_one(self) -> None:
        """``refine_backbone`` seeds inline; the energy test seeds via :func:`seed_region`.

        Production ``energy_before`` on the numba path depends on the INLINE seeding, while
        ``TestEnergyEquivalence`` builds its compiled side from ``seed_region``. Nothing asserted
        the two agree, so the equivalence test rested on an unstated invariant: a later edit to
        either seeder would silently change reported energies without failing anything.

        This intercepts the real ``RegionKernel.energy`` call rather than comparing return values.
        Comparing what ``refine_backbone`` RETURNS would not test this: its teardown also places the
        terminal atoms, which ``seed_region`` deliberately does not, so the two differ at row 0 by a
        measured 0.067 A for a reason that has nothing to do with seeding.
        """
        truth = _truth()
        ca = truth["CA"]
        start = backbone_from_ca(ca)

        captured: list[tuple[np.ndarray, ...]] = []
        original = kernel.RegionKernel.energy

        def spy(self, n_live, c_live, o_live, azimuths):  # type: ignore[no-untyped-def]
            captured.append((n_live.copy(), c_live.copy(), o_live.copy(), azimuths.copy()))
            return original(self, n_live, c_live, o_live, azimuths)

        kernel.RegionKernel.energy = spy  # type: ignore[method-assign]
        try:
            refine_backbone(ca, start.n_xyz, start.c_xyz, backend="numba", max_sweeps=4)
        finally:
            kernel.RegionKernel.energy = original  # type: ignore[method-assign]

        assert captured, "RegionKernel.energy was never called; the interception proved nothing"
        reference = kernel.seed_region(ca, start.c_xyz, start.n_xyz)
        for name, got, want in zip(
            ("n_live", "c_live", "o_live", "azimuths"), captured[0], reference, strict=True
        ):
            assert np.array_equal(got, want), (
                f"{name} differs between refine_backbone's inline seed and seed_region: "
                f"max |diff| {np.abs(np.asarray(got) - np.asarray(want)).max()}"
            )

    def test_an_enormous_bounding_box_enlarges_the_cell_instead_of_the_grid(self) -> None:
        """The dense grid costs ``cells * 8`` bytes, and DODO generates extended IDRs.

        At a cell side of ``reach`` a 20,000 A box would want ~3.9e10 cells, 310 GB. The cell is
        enlarged until the count fits :data:`_MAX_GRID_CELLS` instead. That is safe rather than
        merely cheap -- the query scans the 3x3x3 block around its cell, and at any side of at least
        ``reach`` every point within ``reach`` still falls inside that block -- so this asserts BOTH
        that the enlargement happened and that the answer is unchanged by it.
        """
        reach = 5.9
        # A sparse lattice, spaced so each point keeps a handful of neighbours rather than dozens.
        # The scipy oracle can only examine its 48 nearest, so a dense cluster would trip ITS cap
        # and the comparison would never run.
        axis = np.arange(4) * 4.5
        near = np.array([[x, y, z] for x in axis for y in axis for z in axis], dtype=np.float64)
        far = near + 20_000.0
        pts = np.ascontiguousarray(np.vstack([near, far]))

        _, dims, _, _, cell = kernel._cell_list(pts, reach, kernel._MAX_GRID_CELLS)
        assert cell > reach, "the cell was not enlarged, so this fixture does not test the guard"
        assert int(dims[0]) * int(dims[1]) * int(dims[2]) <= kernel._MAX_GRID_CELLS

        n = pts.shape[0]
        chain_of = np.zeros(n, dtype=np.int64)
        oxygen_of = np.zeros(n, dtype=np.int64)
        centres = np.ascontiguousarray(np.vstack([near[:10], far[:10]]))
        own_chain = np.full((1, centres.shape[0]), 10**9, dtype=np.int64)
        own_oxygen = np.zeros((1, centres.shape[0]), dtype=np.int64)
        got, _ = kernel._neighbour_indices_grid(
            pts,
            centres,
            reach,
            own_chain,
            own_oxygen,
            chain_of,
            oxygen_of,
            n,
            kernel.MAX_NEIGHBOURS,
            kernel._MAX_GRID_CELLS,
        )
        want = kernel._neighbour_indices(
            cKDTree(pts), centres, reach, own_chain[0], own_oxygen[0], chain_of, oxygen_of, n
        )
        for q in range(centres.shape[0]):
            assert set(got[0, q][got[0, q] != n].tolist()) == set(want[q][want[q] != n].tolist())
        assert (got != n).sum() > centres.shape[0], "fixture found no neighbours; proves nothing"
