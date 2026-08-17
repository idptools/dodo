"""Tests for the combined validation entry point, :func:`dodo.validate.validate_structure`.

The three individual validators are tested exhaustively elsewhere. What is tested here is the
layer that combines them, and specifically the two things that layer alone can get wrong:

* the impossible-separation check must run unconditionally, because a validator you can talk out
  of reporting a 0.000 A contact is not doing its job;
* a finding must not be blamed on DODO when DODO did not produce the geometry.

That second point is not pedantry about wording. DODO moves folded domains as rigid bodies and
never regenerates their atoms, so a bad bond in the input is still there afterwards -- preserved
exactly, which is correct behaviour. Before this, ``dodo validate`` on DODO's own freshly rebuilt
dnmt3a printed ``INVALID: 3 bond`` and exited 2, and all three findings were AlphaFold's own
distorted HIS613 imidazole ring: 2.547 A in the input file, 2.548 A in the output. A user reading
that concludes DODO broke their structure.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

from dodo.io import read_structure, write_pdb
from dodo.structure import Structure
from dodo.validate import validate_structure

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "structures"
DNMT3A = FIXTURES / "dnmt3a.pdb"

#: The residue AlphaFold ships with a torn imidazole ring in dnmt3a. Its CE1-ND1 and CD2-NE2
#: bonds measure 2.140 A against a 1.342 A reference -- 169.7 standard deviations out. Nothing
#: DODO does creates or repairs it; it is the natural test case for inherited geometry.
BROKEN_INPUT_RESIDUE = "A:HIS613"


class TestInheritedFindings:
    def test_input_defect_is_present_in_the_output_and_identified_as_inherited(self) -> None:
        from dodo.construct.pipeline import rebuild

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(DNMT3A, seed=0, backbone=False)
        result = validate_structure(report.models[0])

        assert result.bonds is not None
        assert result.bonds.violations, "expected the input's HIS613 defect to survive"
        # Every one of them is on geometry DODO did not build.
        assert result.n_inherited_bond_violations == len(result.bonds.violations)
        assert not result.bonds.of_provenance("rebuilt")
        assert all(BROKEN_INPUT_RESIDUE in v.residue_labels for v in result.bonds.violations), [
            v.residue_labels for v in result.bonds.violations
        ]
        # And the summary says so, rather than reading as DODO's fault.
        assert "inherited from the input" in result.summary()

    def test_the_defect_really_is_in_the_input_at_the_same_distance(self) -> None:
        """Guards the premise. If this drifts, the test above is measuring something else."""
        from dodo.construct.pipeline import rebuild

        original = read_structure(DNMT3A)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            rebuilt = rebuild(DNMT3A, seed=0).models[0]

        def ring_bond(structure: Structure) -> float:
            names = np.asarray(structure.atom_name).astype(str)
            residues = np.asarray(structure.residue_index)
            wanted = {}
            for index in range(structure.n_atoms):
                if structure.residue_label(int(residues[index])) != BROKEN_INPUT_RESIDUE:
                    continue
                if names[index] in ("CE1", "ND1"):
                    wanted[names[index]] = structure.xyz[index]
            assert set(wanted) == {"CE1", "ND1"}
            return float(np.linalg.norm(wanted["CE1"] - wanted["ND1"]))

        before, after = ring_bond(original), ring_bond(rebuilt)
        assert before > 2.0, f"expected a broken ring in the input, measured {before:.3f} A"
        # Rigid-body motion preserves internal geometry, so it must come through unchanged.
        assert abs(after - before) < 1e-3, f"input {before:.4f} A vs output {after:.4f} A"

    def test_a_plain_deposited_file_reports_nothing_as_inherited(self) -> None:
        """Saying "DODO did not build this" is vacuous for a file DODO never touched.

        Saying it anyway would be noise on every structure a user validates, so the inference
        is gated on the structure showing CA-only residues -- DODO's signature.
        """
        result = validate_structure(read_structure(DNMT3A))
        assert result.bonds is not None
        assert result.bonds.violations, "the fixture is expected to have its own defect"
        assert result.n_inherited_bond_violations == 0
        assert "inherited" not in result.summary()

    def test_provenance_survives_a_write_read_cycle_via_the_ca_only_signature(
        self, tmp_path: Path
    ) -> None:
        """A PDB file carries no region metadata, so exact provenance is gone after writing.

        It is still recoverable as an inference rather than a guess: in 2.0 DODO only ever builds
        alpha carbons, so a CA-CA virtual bond is its only possible bond contribution and any
        other bond violation is inherited by construction. This is the path ``dodo validate``
        takes, and it is the one that used to blame DODO for its input.
        """
        from dodo.construct.pipeline import rebuild

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(DNMT3A, seed=0, backbone=False)
        path = tmp_path / "out.pdb"
        write_pdb(report.models[0], path)

        from_file = validate_structure(read_structure(path))
        assert from_file.bonds is not None
        assert from_file.bonds.n_ca_only_residues > 0, "the CA-only signature must be present"
        # Metadata is gone, so nothing is labelled exactly...
        assert not from_file.bonds.of_provenance("input")
        # ...and the inference still attributes every finding correctly.
        assert from_file.n_inherited_bond_violations == len(from_file.bonds.violations)


class TestTheImpossibleCheckCannotBeDisabled:
    def _coincident_pair(self) -> Structure:
        """Two alpha carbons at exactly the same point, inside every normal exclusion."""
        n = 4
        xyz = np.array([[0.0, 0.0, 0.0], [3.81, 0.0, 0.0], [3.81, 0.0, 0.0], [7.62, 0.0, 0.0]])
        return Structure.from_atom_records(
            xyz=xyz,
            atom_name=["CA"] * n,
            element=["C"] * n,
            residue_name=["ALA"] * n,
            residue_number=list(range(1, n + 1)),
            chain_id=["A"] * n,
        )

    @pytest.mark.parametrize(
        ("check_bonds", "check_clashes"),
        [(True, True), (True, False), (False, True), (False, False)],
    )
    def test_reported_whatever_else_is_switched_off(
        self, *, check_bonds: bool, check_clashes: bool
    ) -> None:
        result = validate_structure(
            self._coincident_pair(), check_bonds=check_bonds, check_clashes=check_clashes
        )
        assert result.impossible, "a 0.000 A pair must be reported with every check disabled"
        assert not result.ok
        assert "IMPOSSIBLE" in result.summary()
        assert any(pair.coincident for pair in result.impossible)

    def test_a_sound_structure_reports_nothing(self) -> None:
        result = validate_structure(read_structure(DNMT3A), check_bonds=False, check_clashes=False)
        assert not result.impossible
        assert result.ok
