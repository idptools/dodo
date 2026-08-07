"""The peptide-plane table is auditable: re-derivable from committed frames, with a pinned accuracy.

``dodo.construct.ca_backbone`` bakes a measured 4-CA lookup table into source. Before this, the
table's provenance lived outside the repo -- a set of magic numbers nobody could check. These tests
close that: :mod:`scripts.derive_peptide_table` regenerates the shipped constants exactly from the
committed frames in ``tests/data/backbone/frames`` (100 all-atom IDR frames, stripped to backbone
atoms and gzipped), and the placement's accuracy on those frames is pinned so a regression in either
the table or the placement math is caught.

They also record a confirmed defect: the first peptide unit of every region is placed with the
forward marginal applied at the wrong reference sign, ~3x worse than it should be. See the xfail at
the bottom and ``backbone_plan.md`` (BB-1).
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

from dodo.construct.ca_backbone import backbone_from_ca

REPO_ROOT = Path(__file__).resolve().parents[2]
FRAMES_DIR = REPO_ROOT / "tests" / "data" / "backbone" / "frames"


def _load_deriver() -> object:
    """Load the (non-package) derivation script under ``scripts/`` as a module."""
    spec = importlib.util.spec_from_file_location(
        "derive_peptide_table", REPO_ROOT / "scripts" / "derive_peptide_table.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_DERIVE = _load_deriver()


def _frames() -> list[Path]:
    return sorted(FRAMES_DIR.glob("*.pdb.gz"))


class TestPeptideTableIsReproducible:
    def test_committed_frames_regenerate_the_shipped_table(self) -> None:
        """The derivation, run on the committed frames, reproduces every shipped constant.

        This is what makes the table auditable rather than magic numbers: change the derivation
        method and this fails. Covers ``_C_BY_BIN`` / ``_N_BY_BIN`` and both marginals.
        """
        frames = _frames()
        assert len(frames) == 100, "the committed frame corpus changed size"
        result = _DERIVE.derive(frames)  # type: ignore[attr-defined]
        problems = _DERIVE.verify(result, tolerance=1e-3)  # type: ignore[attr-defined]
        assert not problems, "; ".join(problems)

    def test_the_corpus_is_the_documented_size(self) -> None:
        """19,302 interior peptide units -- the figure ca_backbone is written against."""
        result = _DERIVE.derive(_frames())  # type: ignore[attr-defined]
        assert sum(result["counts"]) == 19302


class TestShippedTableAccuracy:
    """Pin the placement's accuracy on the committed frames.

    Place backbones from CA and measure the error, so a regression in the table or the placement
    math becomes visible.
    """

    @staticmethod
    def _errors() -> tuple[float, float, float]:
        c_err: list[float] = []
        n_err: list[float] = []
        first_c: list[float] = []
        for frame in _frames():
            for residues in _DERIVE.read_backbone(frame).values():  # type: ignore[attr-defined]
                cas = np.array([r["CA"] for r in residues])
                m = len(cas)
                if m < 5:
                    continue
                built = backbone_from_ca(cas)
                # Interior table units are placed correctly; the terminal atoms and the first unit
                # are not, so they are excluded from the interior pins and measured separately.
                for i in range(1, m - 2):
                    c_err.append(float(np.linalg.norm(built.c_xyz[i] - residues[i]["C"])))
                for i in range(2, m - 1):
                    n_err.append(float(np.linalg.norm(built.n_xyz[i] - residues[i]["N"])))
                first_c.append(float(np.linalg.norm(built.c_xyz[0] - residues[0]["C"])))
        return float(np.mean(c_err)), float(np.mean(n_err)), float(np.mean(first_c))

    def test_interior_accuracy_matches_the_documentation(self) -> None:
        """Interior units reproduce the cited held-out errors (C 0.210, N 0.164 A) in-sample."""
        c_err, n_err, _ = self._errors()
        assert c_err < 0.24, f"interior C error {c_err:.3f} A regressed (documented ~0.21)"
        assert n_err < 0.18, f"interior N error {n_err:.3f} A regressed (documented ~0.16)"

    def test_first_peptide_unit_is_placed_accurately(self) -> None:
        """The first unit hits the accuracy its forward marginal delivers (BB-1 sign fix).

        The forward marginal is measured in the frame ``reference = CA[i+1]-CA[i+2]`` builds, and
        ``_backbone_atoms`` now applies it with that same sign. Before the fix the placement used
        the opposite sign and the first unit landed ~1.06 A out; it is now ~0.34 A. This guards
        against a regression that flips the sign back.
        """
        _, _, first_c = self._errors()
        assert first_c < 0.40, (
            f"first-unit C error {first_c:.3f} A -- forward-marginal sign regressed"
        )
