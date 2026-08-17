#!/usr/bin/env python
"""Re-derive DODO's peptide-plane lookup table from all-atom IDR simulation frames.

The table baked into :mod:`dodo.construct.ca_backbone` (``_C_BY_BIN`` / ``_N_BY_BIN`` and the two
marginals) turns backbone placement into a lookup: given four alpha carbons, the peptide unit
between the middle two has one degree of freedom, and the table predicts it. This script is how
that table was measured, committed so the numbers are reproducible rather than a one-off nobody can
check.

Run it over the committed fixture frames (``tests/data/backbone/frames``) and it reproduces the
shipped constants exactly; run it over a larger frame set to re-derive them. ``tests/unit/
test_backbone_table.py`` runs the ``--verify`` path in CI.

Method (must match the placement in ``ca_backbone._backbone_atoms`` term for term)
----------------------------------------------------------------------------------
For a peptide unit ``i`` between residues ``i`` and ``i+1``:

* the frame is anchored at ``CA[i]``: x along ``CA[i+1]-CA[i]``, z along ``cross(reference,
  along)``, y = ``cross(z, x)``;
* the true carbonyl ``C(i)`` and amide ``N(i+1)`` are projected into that frame; their per-bin
  means are the table.

Three predictors, by a unit's position in the chain, each with its own ``reference``:

* **interior** units (both flanking CAs plus a fourth ahead) use ``reference = CA[i]-CA[i-1]`` (the
  incoming virtual bond) and are binned by the pseudo-dihedral ``CA[i-1..i+2]`` -> ``_C_BY_BIN``;
* the **trailing marginal** is the same incoming-reference projection pooled over every interior
  unit, unconditional on the dihedral -> ``_C_MARGINAL`` (used for the last unit, which has no
  fourth CA);
* the **forward marginal** is for the first unit, which has no ``CA[i-1]``. It is measured in the
  forward-looking frame over every unit that has a ``CA[i+2]`` -> ``_C_FORWARD_MARGINAL``.

.. note::

   The forward marginal is measured with ``reference = CA[i+1]-CA[i+2]`` -- the only sign that
   reproduces the shipped ``_C_FORWARD_MARGINAL`` / ``_N_FORWARD_MARGINAL``, and (since the BB-1 fix
   of 2026-08) the sign ``_backbone_atoms`` places the first unit with. An earlier placement used
   the opposite sign and mis-placed the first peptide unit of every region -- ~1.06 A C error
   against the 0.31 A the marginal delivers. Keep the two in step: this derivation and the unit-0
   ``reference`` in ``ca_backbone`` must use the same sign, guarded by
   ``tests/unit/test_backbone_table.py``.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

BIN_WIDTH_DEGREES = 20
N_BINS = 360 // BIN_WIDTH_DEGREES
_BACKBONE_ATOMS = ("N", "CA", "C")


def _open(path: Path):
    """Open a PDB file, transparently decompressing a ``.gz``."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def read_backbone(path: Path) -> dict[str, list[dict[str, np.ndarray]]]:
    """Return, per chain, the residues (in file order) that carry all of N, CA and C.

    Only backbone atoms are read; a stripped N/CA/C(/O) frame is all this needs.
    """
    per_residue: dict[tuple[str, int], dict[str, np.ndarray]] = {}
    order: dict[str, list[int]] = {}
    with _open(path) as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            name = line[12:16].strip()
            if name not in _BACKBONE_ATOMS:
                continue
            chain = line[21]
            resseq = int(line[22:26])
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            per_residue.setdefault((chain, resseq), {})[name] = xyz
            seqs = order.setdefault(chain, [])
            if not seqs or seqs[-1] != resseq:
                seqs.append(resseq)
    chains: dict[str, list[dict[str, np.ndarray]]] = {}
    for chain, seqs in order.items():
        residues = [per_residue[(chain, s)] for s in seqs]
        chains[chain] = [r for r in residues if all(a in r for a in _BACKBONE_ATOMS)]
    return chains


def _local(
    apex: np.ndarray, reference: np.ndarray, along: np.ndarray, point: np.ndarray
) -> np.ndarray:
    """Project ``point`` into the peptide-unit frame anchored at ``apex``.

    Mirrors ``ca_backbone._frames_rows``: x along ``along``, z along ``cross(reference, along)``,
    y = ``cross(z, x)``.
    """
    x = along / np.linalg.norm(along)
    normal = np.cross(reference, along)
    z = normal / np.linalg.norm(normal)
    y = np.cross(z, x)
    return (point - apex) @ np.column_stack((x, y, z))


def _pseudo_dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    """Signed CA pseudo-dihedral in degrees; mirrors ``ca_backbone._dihedral_rows``."""
    axis = (p2 - p1) / np.linalg.norm(p2 - p1)
    v = (p1 - p0) - np.dot(p1 - p0, axis) * axis
    w = (p3 - p2) - np.dot(p3 - p2, axis) * axis
    return float(np.degrees(np.arctan2(np.dot(np.cross(axis, v), w), np.dot(v, w))))


class _Accumulator:
    """Sums of local C and N coordinates, for one predictor."""

    def __init__(self) -> None:
        self.c: list[np.ndarray] = []
        self.n: list[np.ndarray] = []

    def add(self, c: np.ndarray, n: np.ndarray) -> None:
        self.c.append(c)
        self.n.append(n)

    def mean(self) -> tuple[np.ndarray, np.ndarray, int]:
        return np.mean(self.c, axis=0), np.mean(self.n, axis=0), len(self.c)


def derive(frames: list[Path]) -> dict[str, object]:
    """Derive the 1D and 2D binned tables and both marginals from ``frames``.

    Returns a dict with the 1D ``c_by_bin`` / ``n_by_bin`` ``(N_BINS, 3)`` arrays and per-bin
    counts, the 2D ``c_by_bin_2d`` / ``n_by_bin_2d`` ``(N_BINS, N_BINS, 3)`` arrays with their
    counts, and the trailing and forward marginals.

    The 2D table keys a unit on the pseudo-dihedral of its own virtual bond (``tau_i``, over
    ``CA[i-1..i+2]``) *and* the next one (``tau_{i+1}``, over ``CA[i..i+3]``) -- effectively a 5-CA
    predictor for the peptide plane. A unit needs a fifth alpha carbon to have ``tau_{i+1}``, so the
    last interior unit of a chain contributes to the 1D table only. Every 2D cell no unit populated
    is backfilled from the 1D row for ``tau_i``, so the committed table is fully defined and the
    placement path never has to branch on an empty cell.
    """
    bins = [_Accumulator() for _ in range(N_BINS)]
    bins2d: dict[tuple[int, int], _Accumulator] = defaultdict(_Accumulator)
    trailing = _Accumulator()
    forward = _Accumulator()
    for path in frames:
        for residues in read_backbone(path).values():
            cas = [r["CA"] for r in residues]
            m = len(cas)
            for i in range(m - 2):
                along = cas[i + 1] - cas[i]
                c_true, n_true = residues[i]["C"], residues[i + 1]["N"]
                # Forward marginal: every unit with a CA ahead, forward reference (see the
                # module note on its sign).
                fwd_ref = cas[i + 1] - cas[i + 2]
                forward.add(
                    _local(cas[i], fwd_ref, along, c_true),
                    _local(cas[i], fwd_ref, along, n_true),
                )
                if i == 0:
                    continue
                # Interior / trailing: incoming reference.
                back_ref = cas[i] - cas[i - 1]
                local_c = _local(cas[i], back_ref, along, c_true)
                local_n = _local(cas[i], back_ref, along, n_true)
                trailing.add(local_c, local_n)
                if i <= m - 3:  # has a fourth CA -> 1D-binnable on tau_i
                    tau_i = _pseudo_dihedral(cas[i - 1], cas[i], cas[i + 1], cas[i + 2])
                    a = min(int((tau_i + 180.0) // BIN_WIDTH_DEGREES), N_BINS - 1)
                    bins[a].add(local_c, local_n)
                    if i <= m - 4:  # has a fifth CA too -> 2D-binnable on (tau_i, tau_{i+1})
                        tau_j = _pseudo_dihedral(cas[i], cas[i + 1], cas[i + 2], cas[i + 3])
                        b = min(int((tau_j + 180.0) // BIN_WIDTH_DEGREES), N_BINS - 1)
                        bins2d[(a, b)].add(local_c, local_n)

    c_by_bin = np.array([b.mean()[0] for b in bins])
    n_by_bin = np.array([b.mean()[1] for b in bins])
    counts = [b.mean()[2] for b in bins]

    c_by_bin_2d = np.empty((N_BINS, N_BINS, 3))
    n_by_bin_2d = np.empty((N_BINS, N_BINS, 3))
    counts_2d = np.zeros((N_BINS, N_BINS), dtype=np.int64)
    for a in range(N_BINS):
        for b in range(N_BINS):
            acc = bins2d.get((a, b))
            if acc is not None and acc.c:
                mc, mn, cnt = acc.mean()
                c_by_bin_2d[a, b], n_by_bin_2d[a, b], counts_2d[a, b] = mc, mn, cnt
            else:  # unpopulated cell: fall back to the 1D row for tau_i
                c_by_bin_2d[a, b], n_by_bin_2d[a, b] = c_by_bin[a], n_by_bin[a]

    t_c, t_n, _ = trailing.mean()
    f_c, f_n, _ = forward.mean()
    return {
        "c_by_bin": c_by_bin,
        "n_by_bin": n_by_bin,
        "counts": counts,
        "c_by_bin_2d": c_by_bin_2d,
        "n_by_bin_2d": n_by_bin_2d,
        "counts_2d": counts_2d,
        "trailing_c": t_c,
        "trailing_n": t_n,
        "forward_c": f_c,
        "forward_n": f_n,
    }


def _fmt_row(row: np.ndarray, count: int | None = None) -> str:
    body = "(" + ", ".join(f"{v:+.4f}" for v in row) + ")"
    return f"    {body},  # n={count}" if count is not None else body


def _fmt_triple(row: np.ndarray) -> str:
    return "(" + ", ".join(f"{v:+.4f}" for v in row) + ")"


def _fmt_2d_table(name: str, table: np.ndarray, counts: np.ndarray) -> list[str]:
    """Format a ``(N_BINS, N_BINS, 3)`` table as a nested tuple, 3 triples per line.

    One outer entry per ``tau_i`` bin. The ``# n=`` comment sums that row's populated units, so a
    reader can see which rows lean on the 1D backfill.
    """
    lines = [f"{name} = ("]
    for a in range(N_BINS):
        lo = -180 + a * BIN_WIDTH_DEGREES
        lines.append(
            f"    (  # tau_i bin {a} [{lo:+d}, {lo + BIN_WIDTH_DEGREES:+d}) deg, "
            f"n={int(counts[a].sum())}"
        )
        triples = [_fmt_triple(table[a, b]) for b in range(N_BINS)]
        for start in range(0, N_BINS, 3):
            lines.append("        " + ", ".join(triples[start : start + 3]) + ",")
        lines.append("    ),")
    lines.append(")")
    return lines


def emit(result: dict[str, object]) -> str:
    """Format the derived constants as the Python literals in ``ca_backbone``."""
    lines = ["_C_BY_BIN = ("]
    for row, count in zip(result["c_by_bin"], result["counts"], strict=True):  # type: ignore[arg-type]
        lines.append(_fmt_row(row, count))
    lines.append(")")
    lines.append("_N_BY_BIN = (")
    for row, count in zip(result["n_by_bin"], result["counts"], strict=True):  # type: ignore[arg-type]
        lines.append(_fmt_row(row, count))
    lines.append(")")
    lines.append(f"_C_MARGINAL = {_fmt_row(result['trailing_c'])}  # type: ignore[arg-type]")
    lines.append(f"_N_MARGINAL = {_fmt_row(result['trailing_n'])}")
    lines.append(f"_C_FORWARD_MARGINAL = {_fmt_row(result['forward_c'])}")
    lines.append(f"_N_FORWARD_MARGINAL = {_fmt_row(result['forward_n'])}")
    counts_2d = result["counts_2d"]
    lines += _fmt_2d_table("_C_BY_BIN_2D", result["c_by_bin_2d"], counts_2d)  # type: ignore[arg-type]
    lines += _fmt_2d_table("_N_BY_BIN_2D", result["n_by_bin_2d"], counts_2d)  # type: ignore[arg-type]
    return "\n".join(lines)


def verify(result: dict[str, object], tolerance: float) -> list[str]:
    """Compare a derived table to the shipped constants; return a list of mismatches."""
    from dodo.construct import ca_backbone as cb

    problems: list[str] = []

    def check(name: str, got: np.ndarray, want: np.ndarray) -> None:
        """Record a mismatch if ``got`` deviates from ``want`` beyond the tolerance."""
        deviation = float(np.max(np.abs(np.asarray(got) - np.asarray(want))))
        if deviation > tolerance:
            problems.append(f"{name}: max|deviation| {deviation:.5f} A > tolerance {tolerance}")

    check("_C_BY_BIN", result["c_by_bin"], cb._C_TABLE)  # type: ignore[arg-type]
    check("_N_BY_BIN", result["n_by_bin"], cb._N_TABLE)  # type: ignore[arg-type]
    check("_C_MARGINAL", result["trailing_c"], cb._C_MARGINAL)  # type: ignore[arg-type]
    check("_N_MARGINAL", result["trailing_n"], cb._N_MARGINAL)  # type: ignore[arg-type]
    check("_C_FORWARD_MARGINAL", result["forward_c"], cb._C_FORWARD_MARGINAL)  # type: ignore[arg-type]
    check("_N_FORWARD_MARGINAL", result["forward_n"], cb._N_FORWARD_MARGINAL)  # type: ignore[arg-type]
    check("_C_BY_BIN_2D", result["c_by_bin_2d"], cb._C_TABLE_2D)  # type: ignore[arg-type]
    check("_N_BY_BIN_2D", result["n_by_bin_2d"], cb._N_TABLE_2D)  # type: ignore[arg-type]
    return problems


def main(argv: list[str] | None = None) -> int:
    """Derive the table from a frames directory; optionally emit literals and/or verify."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "frames",
        type=Path,
        help="Directory of all-atom (or N/CA/C/O) IDR PDB frames, plain or .gz.",
    )
    parser.add_argument(
        "--emit", action="store_true", help="Print the constants as Python literals."
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Compare the derived table to the shipped ca_backbone constants.",
    )
    parser.add_argument(
        "--tolerance", type=float, default=1e-3, help="Verify tolerance in Angstroms."
    )
    args = parser.parse_args(argv)

    paths = sorted([*args.frames.glob("*.pdb"), *args.frames.glob("*.pdb.gz")])
    if not paths:
        parser.error(f"no .pdb or .pdb.gz frames found in {args.frames}")
    result = derive(paths)
    total = sum(result["counts"])  # type: ignore[arg-type]
    print(f"derived from {len(paths)} frame(s): {total} interior peptide units", file=sys.stderr)

    if args.emit:
        print(emit(result))
    if args.verify:
        problems = verify(result, args.tolerance)
        if problems:
            print("VERIFY FAILED:", file=sys.stderr)
            for problem in problems:
                print(f"  {problem}", file=sys.stderr)
            return 1
        print("VERIFY OK: derived table matches the shipped constants", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
