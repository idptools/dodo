"""Tests for the backbone refinement pass.

Refinement holds every alpha carbon fixed and rotates each peptide unit about the CA-CA axis --
the single degree of freedom a rigid peptide unit has once both alpha carbons are pinned. It is
scored on four things at once: N-CA-C angle, steric clashes, phi/psi plausibility, and not
disturbing what was already good.

That "at once" is the whole difficulty, and the tests here exist because two separate attempts at
the weighting made things worse in ways that a single-criterion check would have called a success:

* the first weighting improved bond angles while pushing dihedrals *out* of their basins, from
  87.9% down to 68.4%;
* the clash term initially counted covalent bonds as clashes -- 14,859 fake contacts against one
  real one, because an N-CA bond at 1.458 A trivially fails a 1.0 A separation test. Any weighting
  tuned against that signal is tuned against noise.

So the central test is :meth:`TestRefinementHelps.test_every_criterion_improves_or_holds`: nothing
is allowed to regress in exchange for something else improving.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dodo.constants import (
    C_N_PEPTIDE_BOND_LENGTH,
    C_O_BOND_LENGTH,
    CA_C_BOND_LENGTH,
    N_CA_BOND_LENGTH,
    N_CA_C_ANGLE,
)
from dodo.construct.backbone_refine import refine_backbone
from dodo.construct.ca_backbone import backbone_from_ca
from dodo.io import read_structure

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "structures"
TRUTH = FIXTURES / "idr_frame_backbone.pdb"


def _truth() -> dict[str, np.ndarray]:
    structure = read_structure(TRUTH)
    names = np.asarray([str(n) for n in structure.atom_name])
    return {atom: structure.xyz[names == atom] for atom in ("N", "CA", "C", "O")}


def _angles(n_xyz: np.ndarray, ca: np.ndarray, c_xyz: np.ndarray) -> np.ndarray:
    """N-CA-C bond angle at every residue, in degrees."""
    a = n_xyz - ca
    b = c_xyz - ca
    cosine = np.sum(a * b, axis=1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def _dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> np.ndarray:
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    axis = b1 / np.linalg.norm(b1, axis=1)[:, None]
    v = b0 - np.sum(b0 * axis, axis=1)[:, None] * axis
    w = b2 - np.sum(b2 * axis, axis=1)[:, None] * axis
    return np.degrees(np.arctan2(np.sum(np.cross(axis, v) * w, axis=1), np.sum(v * w, axis=1)))


def _phi_psi(n_xyz: np.ndarray, ca: np.ndarray, c_xyz: np.ndarray) -> tuple[np.ndarray, ...]:
    phi = _dihedral(c_xyz[:-1], n_xyz[1:], ca[1:], c_xyz[1:])
    psi = _dihedral(n_xyz[:-1], ca[:-1], c_xyz[:-1], n_xyz[1:])
    return phi, psi


def _in_populated_region(phi: np.ndarray, psi: np.ndarray) -> float:
    """Fraction of residues in the part of the map real IDRs actually occupy.

    Measured on 19,602 phi/psi pairs from all-atom IDR simulation: the distribution is dominated by
    extended and PPII geometry with a median around phi -107, psi +122 -- NOT the alpha-helical
    basin a folded-protein intuition reaches for. A generous box around that, plus the left-handed
    region glycine uses.
    """
    extended = (phi > -180) & (phi < -30) & ((psi > 60) | (psi < -150))
    helical = (phi > -160) & (phi < -30) & (psi > -80) & (psi < 10)
    left = (phi > 20) & (phi < 100) & (psi > 0) & (psi < 100)
    return float(np.mean(extended | helical | left))


@pytest.fixture(scope="module")
def initial() -> tuple[np.ndarray, ...]:
    """Alpha carbons plus the unrefined placement, which refinement is asked to improve."""
    truth = _truth()
    ca = truth["CA"]
    start = backbone_from_ca(ca)
    return ca, start.n_xyz, start.c_xyz, start.o_xyz


class TestInvariants:
    """What refinement is structurally forbidden to change."""

    def test_alpha_carbons_are_untouched(self, initial) -> None:
        """Alpha carbons are DODO's answer, and this pass is not allowed to revise them."""
        ca, n_xyz, c_xyz, _ = initial
        before = ca.copy()
        refine_backbone(ca, n_xyz, c_xyz)
        assert np.array_equal(ca, before), "refinement mutated the alpha carbons it was given"

    def test_inputs_are_not_mutated(self, initial) -> None:
        """Callers keep the unrefined placement to compare against, so it must survive the call."""
        ca, n_xyz, c_xyz, _ = initial
        n_before, c_before = n_xyz.copy(), c_xyz.copy()
        refine_backbone(ca, n_xyz, c_xyz)
        assert np.array_equal(n_xyz, n_before)
        assert np.array_equal(c_xyz, c_before)

    @pytest.mark.parametrize(
        ("bond", "ideal"),
        [
            ("N-CA", N_CA_BOND_LENGTH),
            ("CA-C", CA_C_BOND_LENGTH),
            ("C-O", C_O_BOND_LENGTH),
        ],
    )
    def test_bonds_to_alpha_carbons_stay_exact(self, initial, bond: str, ideal: float) -> None:
        """A rigid unit rotated about the CA-CA axis cannot change a bond length.

        So this is really a test that the rotation is a rotation. Trading bond length for a better
        angle would be the easiest way for this pass to "improve" its score while producing a worse
        structure.

        The peptide C-N bond is deliberately not in this list; see
        :meth:`test_peptide_bond_absorbs_non_ideal_input_spacing`.
        """
        ca, n_xyz, c_xyz, _ = initial
        result = refine_backbone(ca, n_xyz, c_xyz)
        observed = {
            "N-CA": np.linalg.norm(result.n_xyz - ca, axis=1),
            "CA-C": np.linalg.norm(result.c_xyz - ca, axis=1),
            "C-O": np.linalg.norm(result.o_xyz - result.c_xyz, axis=1),
        }[bond]
        assert np.allclose(observed, ideal, atol=1e-6), (
            f"{bond} drifted to {observed.min():.4f}-{observed.max():.4f} A"
        )

    def test_peptide_bond_is_exact_on_an_ideally_spaced_trace(self) -> None:
        """On the traces DODO actually produces, C-N comes out constant.

        The walk engine emits exactly 3.81 A and STARLING output is projected onto it, so this is
        the case that matters in practice. It is 1.3334 A rather than 1.329 because the canonical
        peptide unit determines a CA-CA distance of 3.8040 and DODO uses 3.81; the 0.004 A is a
        known, documented consequence of that choice, not drift.
        """
        import dodo

        ca = dodo.build_from_sequence("GRNQNGGGYQNYNNQGYQGHGGQHQNNYNQY", seed=0).models[0].ca_xyz
        assert np.allclose(np.linalg.norm(np.diff(ca, axis=0), axis=1), 3.81, atol=1e-6)
        start = backbone_from_ca(ca)
        result = refine_backbone(ca, start.n_xyz, start.c_xyz)
        peptide = np.linalg.norm(result.n_xyz[1:] - result.c_xyz[:-1], axis=1)
        assert peptide.std() < 1e-6, "C-N varies on a uniformly spaced trace"
        assert abs(float(peptide.mean()) - C_N_PEPTIDE_BOND_LENGTH) < 0.005

    def test_peptide_bond_absorbs_non_ideal_input_spacing(self, initial) -> None:
        """And the caller is told, rather than getting stretched bonds with no explanation.

        A rigid peptide unit *determines* its CA-CA distance, so a trace that disagrees cannot be
        fitted by one without something giving. C is placed off CA(i) and N off CA(i+1), so the
        discrepancy lands in the peptide bond between them -- one for one, at r = 1.0 against the
        input spacing. The ground-truth fixture is a real simulation frame spanning 3.747-3.877 A,
        so it exercises exactly this.

        Refinement is not free to fix it (the alpha carbons are fixed input), so the requirement is
        that it says so. Silently returning a 1.288 A peptide bond is the failure mode here.
        """
        ca, n_xyz, c_xyz, _ = initial
        spacing = np.linalg.norm(np.diff(ca, axis=0), axis=1)
        assert spacing.std() > 1e-3, "fixture no longer has non-uniform spacing; test is void"

        result = refine_backbone(ca, n_xyz, c_xyz)
        peptide = np.linalg.norm(result.n_xyz[1:] - result.c_xyz[:-1], axis=1)
        assert np.corrcoef(spacing, peptide)[0, 1] > 0.99
        assert any("Regularize" in note for note in result.notes), (
            f"non-ideal spacing was not reported; notes were {result.notes}"
        )

    def test_output_is_finite(self, initial) -> None:
        ca, n_xyz, c_xyz, _ = initial
        result = refine_backbone(ca, n_xyz, c_xyz)
        for name, array in (("N", result.n_xyz), ("C", result.c_xyz), ("O", result.o_xyz)):
            assert np.isfinite(array).all(), f"{name} contains non-finite coordinates"


class TestRefinementHelps:
    def test_the_objective_goes_down(self, initial) -> None:
        ca, n_xyz, c_xyz, _ = initial
        result = refine_backbone(ca, n_xyz, c_xyz)
        assert result.energy_after <= result.energy_before

    def test_every_criterion_improves_or_holds(self, initial) -> None:
        """The test the first two weightings would have failed.

        Refinement optimises four things at once, and a weighting that buys one with another is
        worse than not refining. Each criterion is checked separately against the unrefined
        placement, with a small tolerance so noise in an already-good measure is not a failure.
        """
        ca, n_xyz, c_xyz, _ = initial
        result = refine_backbone(ca, n_xyz, c_xyz)

        angle_before = float(np.std(_angles(n_xyz, ca, c_xyz)[1:-1]))
        angle_after = float(np.std(_angles(result.n_xyz, ca, result.c_xyz)[1:-1]))
        assert angle_after <= angle_before + 0.1, (
            f"N-CA-C spread worsened: {angle_before:.2f} -> {angle_after:.2f} deg"
        )

        before = _in_populated_region(*_phi_psi(n_xyz, ca, c_xyz))
        after = _in_populated_region(*_phi_psi(result.n_xyz, ca, result.c_xyz))
        assert after >= before - 0.02, (
            f"fraction of residues in populated phi/psi regions fell: {before:.3f} -> {after:.3f}"
        )

    def test_bond_angles_approach_ideal(self, initial) -> None:
        ca, n_xyz, c_xyz, _ = initial
        result = refine_backbone(ca, n_xyz, c_xyz)
        after = _angles(result.n_xyz, ca, result.c_xyz)[1:-1]
        assert abs(float(np.mean(after)) - N_CA_C_ANGLE) < 2.0
        # Real IDR backbone holds this angle to about 2.8 degrees; allow a little more.
        assert float(np.std(after)) < 5.0

    @pytest.mark.parametrize("atom", ["N", "C", "O"])
    def test_accuracy_against_truth_improves(self, initial, atom: str) -> None:
        """Refinement is scored on physics, but it must not walk *away* from the real answer.

        Nothing in the objective mentions the true coordinates, so this could in principle come out
        either way. It is the check that the physical terms are proxies for something real.
        """
        ca, n_xyz, c_xyz, o_xyz = initial
        truth = _truth()
        result = refine_backbone(ca, n_xyz, c_xyz)
        start = {"N": n_xyz, "C": c_xyz, "O": o_xyz}[atom]
        end = {"N": result.n_xyz, "C": result.c_xyz, "O": result.o_xyz}[atom]
        before = float(np.mean(np.linalg.norm(start - truth[atom], axis=1)))
        after = float(np.mean(np.linalg.norm(end - truth[atom], axis=1)))
        assert after <= before + 0.02, f"{atom} error grew: {before:.3f} -> {after:.3f} A"


class TestClashes:
    def test_covalent_bonds_are_not_counted_as_clashes(self, initial) -> None:
        """The 14,859-fake-contacts bug.

        An N-CA bond is 1.458 A, so any clash test that does not exclude bonded and near-bonded
        atoms fires on every residue in the chain. The symptom was not a crash but a useless
        gradient: the clash term dominated the objective with signal that no rotation could reduce.

        A clean chain refined against itself should report essentially no clashes.
        """
        ca, n_xyz, c_xyz, _ = initial
        result = refine_backbone(ca, n_xyz, c_xyz)
        placed = np.vstack([result.n_xyz, ca, result.c_xyz, result.o_xyz])
        residue = np.tile(np.arange(ca.shape[0]), 4)
        gaps = np.linalg.norm(placed[:, None, :] - placed[None, :, :], axis=2)
        far_apart = np.abs(residue[:, None] - residue[None, :]) > 2
        assert int(np.sum((gaps < 2.5) & far_apart)) // 2 == 0

    def test_obstacles_are_avoided(self, initial) -> None:
        """An obstacle placed on an atom must move it.

        Obstacles are how a rebuilt region learns about the folded domains around it, so this
        failing means backbone placement inside a real structure ignores the protein it sits in.
        """
        ca, n_xyz, c_xyz, _ = initial
        target = 20
        obstacle = c_xyz[target][None, :].copy()
        result = refine_backbone(ca, n_xyz, c_xyz, obstacles=obstacle)
        moved = float(np.linalg.norm(result.c_xyz[target] - c_xyz[target]))
        assert moved > 0.5, f"atom sitting on an obstacle moved only {moved:.3f} A"
        assert float(np.linalg.norm(result.c_xyz[target] - obstacle[0])) > moved / 2

    def test_empty_obstacles_are_accepted(self, initial) -> None:
        """A region with nothing around it is normal, not an edge case to crash on."""
        ca, n_xyz, c_xyz, _ = initial
        result = refine_backbone(ca, n_xyz, c_xyz, obstacles=np.zeros((0, 3)))
        assert np.isfinite(result.n_xyz).all()


class TestControl:
    def test_zero_sweeps_still_canonicalizes(self, initial) -> None:
        """``max_sweeps=0`` runs no refinement but is not an identity.

        Entry re-places every unit from the canonical peptide geometry at the azimuth measured off
        the input, which is what makes the search one-parameter-per-unit. So zero sweeps gives
        "snapped to ideal geometry, not otherwise adjusted" rather than the input back. Worth
        pinning down, because it is the natural way to ask for an unrefined baseline and the
        distinction is invisible in the output file.
        """
        ca, n_xyz, c_xyz, _ = initial
        result = refine_backbone(ca, n_xyz, c_xyz, max_sweeps=0)
        assert result.sweeps == 0
        assert result.energy_after == result.energy_before
        # Snapped, so not identical -- but close, and with the alpha-carbon bonds now exact.
        assert not np.allclose(result.n_xyz, n_xyz)
        assert float(np.max(np.linalg.norm(result.n_xyz - n_xyz, axis=1))) < 0.5
        assert np.allclose(np.linalg.norm(result.n_xyz - ca, axis=1), N_CA_BOND_LENGTH, atol=1e-6)

    def test_convergence_is_reported(self, initial) -> None:
        ca, n_xyz, c_xyz, _ = initial
        result = refine_backbone(ca, n_xyz, c_xyz, max_sweeps=100)
        assert result.converged
        assert result.sweeps < 100

    def test_stopping_early_is_reported_as_unconverged(self, initial) -> None:
        """A caller must be able to tell a finished refinement from one that ran out of sweeps."""
        ca, n_xyz, c_xyz, _ = initial
        result = refine_backbone(ca, n_xyz, c_xyz, max_sweeps=1)
        assert not result.converged

    def test_is_deterministic(self, initial) -> None:
        """Coordinate descent over a fixed candidate set; two runs must agree exactly."""
        ca, n_xyz, c_xyz, _ = initial
        first = refine_backbone(ca, n_xyz, c_xyz)
        second = refine_backbone(ca, n_xyz, c_xyz)
        assert np.array_equal(first.n_xyz, second.n_xyz)
        assert np.array_equal(first.c_xyz, second.c_xyz)

    def test_notes_report_the_resulting_angle(self, initial) -> None:
        ca, n_xyz, c_xyz, _ = initial
        result = refine_backbone(ca, n_xyz, c_xyz)
        assert result.notes
        assert "N-CA-C" in result.notes[0]


def test_runs_fast_enough_to_be_default_one_day() -> None:
    """Speed is a feature here, so it gets an assertion rather than a comment.

    An early version scored the whole chain for every candidate azimuth and took over ten minutes
    for ten chains; localising the objective and vectorising over candidates brought a 900-residue
    structure to a few seconds. This bound is loose enough for a slow machine but would catch a
    return to whole-chain scoring, which is orders of magnitude slower rather than marginally.
    """
    import time

    truth = _truth()
    ca = truth["CA"]
    start = backbone_from_ca(ca)
    began = time.perf_counter()
    refine_backbone(ca, start.n_xyz, start.c_xyz)
    elapsed = time.perf_counter() - began
    assert elapsed < 10.0, f"refining 60 residues took {elapsed:.1f} s"
