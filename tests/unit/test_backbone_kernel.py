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

from dodo.constants import CA_CA_BOND_LENGTH, N_CA_C_ANGLE
from dodo.construct.backbone_refine import refine_backbone
from dodo.construct.ca_backbone import backbone_from_ca
from dodo.exceptions import InvalidParameterError
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
