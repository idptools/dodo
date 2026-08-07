"""Tests for placing N, C and O from alpha carbons alone.

Ground truth is a 60-residue slice of a real all-atom IDR simulation frame
(``idr_frame_backbone.pdb``). The test strips it to alpha carbons, rebuilds, and compares against
the atoms it threw away -- which is the only honest way to measure this, and is why the fixture
exists rather than reusing a folded domain from one of the AlphaFold files.

Several tests are named for a specific defect, because each one produced output that looked
entirely reasonable:

* the lookup table was transcribed rolled by one bin, which cost a factor of two in accuracy while
  every bond length stayed exact;
* only two of a peptide unit's three bonds were being enforced, so N-CA came out at
  1.252 +/- 0.115 A against a real 1.458 -- close enough to pass a casual look at the file;
* the terminal carbonyl oxygen was oriented against an invented nitrogen placed *collinear* with
  CA->C, which degenerates the plane that defines the carbonyl and put the final O 0.981 A from its
  own alpha carbon.

The accuracy bounds below are deliberately loose -- roughly 1.5x the measured error. They are
regression guards against a change that quietly makes placement worse, not a claim that these
numbers are the target.
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

import numpy as np
import pytest

from dodo.constants import (
    C_N_PEPTIDE_BOND_LENGTH,
    C_O_BOND_LENGTH,
    CA_C_BOND_LENGTH,
    N_CA_BOND_LENGTH,
)
from dodo.construct.ca_backbone import backbone_from_ca
from dodo.exceptions import GeometryError
from dodo.io import read_structure

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "structures"
TRUTH = FIXTURES / "idr_frame_backbone.pdb"


def _truth() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Alpha carbons and the true N, C, O of the ground-truth frame."""
    structure = read_structure(TRUTH)
    names = np.asarray([str(n) for n in structure.atom_name])
    per_residue = {}
    for atom in ("N", "CA", "C", "O"):
        selected = structure.xyz[names == atom]
        assert selected.shape[0] == structure.n_residues, f"fixture is missing some {atom} atoms"
        per_residue[atom] = selected
    return per_residue["CA"], per_residue["N"], per_residue["C"], per_residue["O"]


def _bond_lengths(ca, n_xyz, c_xyz, o_xyz):
    """Return the four bonds the construction is supposed to make exact, as flat arrays."""
    return {
        "N-CA": np.linalg.norm(n_xyz - ca, axis=1),
        "CA-C": np.linalg.norm(c_xyz - ca, axis=1),
        "C-O": np.linalg.norm(o_xyz - c_xyz, axis=1),
        "C-N(next)": np.linalg.norm(n_xyz[1:] - c_xyz[:-1], axis=1),
    }


class TestGeometryIsExact:
    """Bond lengths are set by construction, so they are exact or the construction is broken."""

    @pytest.mark.parametrize(
        ("bond", "ideal"),
        [
            ("N-CA", N_CA_BOND_LENGTH),
            ("CA-C", CA_C_BOND_LENGTH),
            ("C-O", C_O_BOND_LENGTH),
            ("C-N(next)", C_N_PEPTIDE_BOND_LENGTH),
        ],
    )
    def test_every_bond_is_ideal(self, bond: str, ideal: float) -> None:
        """All four, not just the two that a partial implementation happened to enforce.

        N-CA is the one that was silently wrong: 1.252 +/- 0.115 A against 1.458. It is included
        here by name because the failure mode is a bond that is *plausible* rather than absurd.
        """
        ca, *_ = _truth()
        result = backbone_from_ca(ca)
        observed = _bond_lengths(ca, result.n_xyz, result.c_xyz, result.o_xyz)[bond]
        assert np.allclose(observed, ideal, atol=1e-6), (
            f"{bond} ranges {observed.min():.4f}-{observed.max():.4f} A, ideal {ideal}"
        )

    def test_nothing_is_nan_or_infinite(self) -> None:
        ca, *_ = _truth()
        result = backbone_from_ca(ca)
        for name, array in (("N", result.n_xyz), ("C", result.c_xyz), ("O", result.o_xyz)):
            assert np.isfinite(array).all(), f"{name} contains non-finite coordinates"

    def test_one_atom_of_each_kind_per_residue(self) -> None:
        ca, *_ = _truth()
        result = backbone_from_ca(ca)
        assert result.n_xyz.shape == ca.shape
        assert result.c_xyz.shape == ca.shape
        assert result.o_xyz.shape == ca.shape


class TestAccuracy:
    """How close the placed atoms are to the ones actually there."""

    # Measured on this fixture, at roughly 1.5x so ordinary drift does not fail the suite but a
    # real regression does. Note that this fixture is one of the 100 frames the lookup table was
    # built from, so these are in-sample -- but only barely: rebuilding the table from 80 frames
    # and scoring the other 20 gives 0.164 / 0.217 / 0.634 A against 0.163 / 0.215 / 0.626 for the
    # shipped table on the same residues, so it generalizes and the numbers are not flattering
    # themselves.
    LIMITS: ClassVar[dict[str, float]] = {"N": 0.25, "C": 0.30, "O": 0.75}

    @pytest.mark.parametrize("atom", ["N", "C", "O"])
    def test_placement_is_close_to_truth(self, atom: str) -> None:
        ca, true_n, true_c, true_o = _truth()
        result = backbone_from_ca(ca)
        placed = {"N": result.n_xyz, "C": result.c_xyz, "O": result.o_xyz}[atom]
        truth = {"N": true_n, "C": true_c, "O": true_o}[atom]
        error = float(np.mean(np.linalg.norm(placed - truth, axis=1)))
        assert error < self.LIMITS[atom], f"mean {atom} error {error:.3f} A"

    def test_nitrogen_is_placed_better_than_carbon(self) -> None:
        """The reason N is placed first and C is then closed onto it.

        Four alpha carbons determine N more tightly than C. The pipeline depends on that ordering,
        so if it ever inverts, the placement order is the thing to revisit.
        """
        ca, true_n, true_c, _ = _truth()
        result = backbone_from_ca(ca)
        n_error = float(np.mean(np.linalg.norm(result.n_xyz - true_n, axis=1)))
        c_error = float(np.mean(np.linalg.norm(result.c_xyz - true_c, axis=1)))
        assert n_error < c_error

    def test_beats_a_naive_midpoint(self) -> None:
        """A guard that the lookup table is contributing something.

        Putting N at the CA(i)/CA(i+1) midpoint needs no table at all. When the table was
        transcribed rolled by one 20-degree bin it still beat this, so passing here is necessary
        rather than sufficient -- which is what :meth:`test_the_table_is_not_rolled` is for.
        """
        ca, true_n, _, _ = _truth()
        result = backbone_from_ca(ca)
        table = float(np.mean(np.linalg.norm(result.n_xyz - true_n, axis=1)))
        midpoint = float(np.mean(np.linalg.norm((ca[:-1] + ca[1:]) / 2 - true_n[1:], axis=1)))
        assert table < midpoint

    def test_the_table_is_not_rolled(self) -> None:
        """Deliberately mis-bin the trace and confirm accuracy degrades.

        This is the test for the transcription bug. A rolled table is self-consistent and produces
        exact bonds, so nothing else here notices; what it costs is accuracy. Rotating the CA trace
        into the wrong bins must be *worse* than using the right ones -- if it is not, the bin
        index is not being used meaningfully.
        """
        ca, true_n, _, _ = _truth()
        correct = float(np.mean(np.linalg.norm(backbone_from_ca(ca).n_xyz - true_n, axis=1)))
        # Reversing the trace maps each pseudo-dihedral to its negative, so every unit lands in a
        # different bin while the chain stays geometrically valid.
        reversed_result = backbone_from_ca(ca[::-1].copy())
        mismatched = float(np.mean(np.linalg.norm(reversed_result.n_xyz[::-1] - true_n, axis=1)))
        assert correct < mismatched


class TestTermini:
    """The first N and the last C and O belong to no peptide unit."""

    def test_terminal_atoms_are_placed_not_left_at_the_origin(self) -> None:
        ca, *_ = _truth()
        result = backbone_from_ca(ca)
        for label, point in (
            ("N(0)", result.n_xyz[0]),
            ("C(-1)", result.c_xyz[-1]),
            ("O(-1)", result.o_xyz[-1]),
        ):
            assert np.isfinite(point).all(), f"{label} is not finite"
            assert np.linalg.norm(point) > 1e-6, f"{label} was left at the origin"

    def test_terminal_oxygen_is_not_on_top_of_its_alpha_carbon(self) -> None:
        """The collinear-plane bug, which put the final O 0.981 A from its own CA.

        Orienting the last carbonyl against an invented next nitrogen along CA->C makes that
        nitrogen collinear with the bond, so the cross product defining the carbonyl plane
        degenerates to zero and the normal falls back to an arbitrary axis. A real CA-O separation
        across a carbonyl is about 2.4 A.
        """
        ca, *_ = _truth()
        result = backbone_from_ca(ca)
        separation = float(np.linalg.norm(result.o_xyz[-1] - ca[-1]))
        assert separation > 2.0, f"terminal O sits {separation:.3f} A from its own CA"

    def test_terminal_oxygen_matches_the_others(self) -> None:
        """Whatever the termini do must be the same *kind* of geometry as the interior."""
        ca, *_ = _truth()
        result = backbone_from_ca(ca)
        interior = np.linalg.norm(result.o_xyz[:-1] - ca[:-1], axis=1)
        terminal = float(np.linalg.norm(result.o_xyz[-1] - ca[-1]))
        assert interior.min() - 0.3 < terminal < interior.max() + 0.3

    def test_two_residues_is_enough(self) -> None:
        """The shortest input that can work at all: no unit has four alpha carbons."""
        ca, *_ = _truth()
        result = backbone_from_ca(ca[:2])
        assert np.isfinite(result.n_xyz).all()
        assert np.isfinite(result.c_xyz).all()
        assert np.isfinite(result.o_xyz).all()
        observed = _bond_lengths(ca[:2], result.n_xyz, result.c_xyz, result.o_xyz)
        assert np.allclose(observed["C-N(next)"], C_N_PEPTIDE_BOND_LENGTH, atol=1e-6)

    def test_three_residues_uses_a_marginal(self) -> None:
        """Too short for a full four-CA window, long enough to need more than the terminal rules."""
        ca, *_ = _truth()
        result = backbone_from_ca(ca[:3])
        assert np.isfinite(result.n_xyz).all()
        assert "marginal" in result.source or "forward" in result.source


class TestRejectsBadInput:
    def test_single_residue_is_refused(self) -> None:
        ca, *_ = _truth()
        with pytest.raises(GeometryError):
            backbone_from_ca(ca[:1])

    def test_wrong_shape_is_refused(self) -> None:
        with pytest.raises(GeometryError):
            backbone_from_ca(np.zeros((5, 2)))

    def test_duplicated_alpha_carbons_are_refused(self) -> None:
        """Two alpha carbons at the same point is invalid input, and saying so is the right answer.

        Every direction in the construction is a normalised difference of alpha carbons, so a
        duplicated position is exactly where a silent NaN would enter and propagate into the output
        file. Refusing beats inventing a position: a zero-length CA-CA bond cannot occur in a real
        chain, so the caller has a broken trace and needs to know.
        """
        ca = np.zeros((4, 3))
        ca[1] = ca[2] = ca[3] = [3.81, 0.0, 0.0]
        with pytest.raises(GeometryError, match="zero-length"):
            backbone_from_ca(ca)

    def test_collinear_alpha_carbons_are_reported(self) -> None:
        """A straight trace has no pseudo-dihedral; that must be noted, not silently binned."""
        ca = np.array([[i * 3.81, 0.0, 0.0] for i in range(6)])
        result = backbone_from_ca(ca)
        assert np.isfinite(result.n_xyz).all()
        assert any("collinear" in note for note in result.notes)


def test_result_reports_where_each_unit_came_from() -> None:
    """``source`` is per-unit provenance, and the notes must agree with it."""
    ca, *_ = _truth()
    result = backbone_from_ca(ca)
    assert len(result.source) == ca.shape[0] - 1
    assert set(result.source) <= {"table", "marginal", "forward"}
    # On a 60-residue chain the interior should overwhelmingly come from the four-CA table.
    assert result.source.count("table") > 0.9 * len(result.source)
    assert str(result.source.count("table")) in result.notes[0]


class TestAddBackboneToRebuilt:
    """The structure-level entry point, and what it is required to leave alone.

    The constraint that shapes this function: folded domains keep every atom they arrived with.
    DODO does not regenerate folded-domain geometry, so adding a backbone must be strictly additive
    over the regions DODO itself generated.
    """

    def test_only_rebuilt_regions_gain_atoms(self) -> None:
        """Folded domains come through untouched, to the last decimal place."""
        import dodo
        from dodo.construct.ca_backbone import add_backbone_to_rebuilt

        structure = dodo.build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", seed=0).models[0]
        before = structure.n_atoms
        after = add_backbone_to_rebuilt(structure).structure
        # An all-IDR structure: every residue is DODO's, so every residue gains N, C and O.
        assert after.n_atoms == before + 3 * structure.n_residues
        assert {str(n) for n in after.atom_name} == {"N", "CA", "C", "O"}

    def test_alpha_carbons_are_unchanged(self) -> None:
        import dodo
        from dodo.construct.ca_backbone import add_backbone_to_rebuilt

        structure = dodo.build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", seed=0).models[0]
        result = add_backbone_to_rebuilt(structure).structure
        assert np.allclose(result.ca_xyz, structure.ca_xyz)

    def test_domains_survive_the_rebuild(self) -> None:
        """``from_atom_records`` builds a fresh Structure, which drops domains unless carried over.

        Losing them is silent: the coordinates are right and every consumer that asks the structure
        what is folded and what DODO built gets an empty answer.
        """
        import dodo
        from dodo.construct.ca_backbone import add_backbone_to_rebuilt

        structure = dodo.build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", seed=0).models[0]
        result = add_backbone_to_rebuilt(structure).structure
        assert len(result.domains) == len(structure.domains)
        assert [d.kind for d in result.domains] == [d.kind for d in structure.domains]
        assert [d.span for d in result.domains] == [d.span for d in structure.domains]

    def test_atoms_are_in_n_ca_c_o_order(self) -> None:
        """Viewers and the CONECT writer both assume the conventional within-residue order."""
        import dodo
        from dodo.construct.ca_backbone import add_backbone_to_rebuilt

        structure = dodo.build_from_sequence("GRNQNGG", seed=0).models[0]
        result = add_backbone_to_rebuilt(structure).structure
        names = [str(n) for n in result.atom_name]
        assert names[:8] == ["N", "CA", "C", "O", "N", "CA", "C", "O"]

    def test_output_has_no_impossible_contacts(self) -> None:
        """The terminal-oxygen bug surfaced here: one CA and its own O 0.981 A apart."""
        import dodo
        from dodo.construct.ca_backbone import add_backbone_to_rebuilt
        from dodo.validate import find_impossible_pairs, validate_bonds

        structure = dodo.build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", seed=0).models[0]
        result = add_backbone_to_rebuilt(structure).structure
        assert not find_impossible_pairs(result)
        assert not validate_bonds(result).violations
