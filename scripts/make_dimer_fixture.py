"""Regenerate ``tests/data/structures/dnmt3a_dimer.pdb``, the minimal complex fixture.

Two copies of ``dnmt3a.pdb`` as chains A and B, the second rotated 180 degrees and pushed along
+y until the two C-terminal folded domains touch at exactly 4.00 A. That gives the smallest
structure that exercises everything the complex code has to get right:

* two chains, two folded domains each, so folded-domain repositioning has something to move;
* one large inter-chain interface (5,540 atom-atom contacts within 5 A) that must survive;
* an interface-induced folded domain -- a 30-residue helix that reads as disordered in the
  monomer and as folded in the dimer, because burial is scored over the whole structure.

Synthetic on purpose. The real AlphaFold 3 complexes in
``tests/data/complexes_full_sequences_modeled`` are what keep the thresholds honest; this one
isolates a single property per assertion and is small enough to run in the unit suite.

Deterministic: no random numbers, so re-running reproduces the file byte for byte.

    python scripts/make_dimer_fixture.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from dodo.geometry.transforms import rotation_from_axis_angle
from dodo.io import read_structure
from dodo.regions.identify import assign_regions

#: Closest heavy-atom approach between the two chains' C-terminal folded domains, in Angstroms.
#: Inside the 5 A contact radius, outside the 3.2 A clash distance -- a real interface that is
#: not also a clash the pipeline would have to resolve.
TARGET_APPROACH = 4.0

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tests" / "data" / "structures" / "dnmt3a.pdb"
DESTINATION = ROOT / "tests" / "data" / "structures" / "dnmt3a_dimer.pdb"


def main() -> None:
    lines = [line for line in SOURCE.read_text().splitlines(keepends=True) if line.startswith("ATOM")]
    structure = read_structure(SOURCE)
    assign_regions(structure)

    folded = [d for d in structure.chains[0].domains if d.kind.value == "folded"]
    interface_domain = folded[-1]
    atoms = structure.atom_slice_for_residues(
        interface_domain.span.start, interface_domain.span.stop
    )
    centre = structure.xyz[atoms].mean(axis=0)

    flipped = (structure.xyz - centre) @ rotation_from_axis_angle(
        np.array([1.0, 0.0, 0.0]), np.pi
    ).T + centre
    direction = np.array([0.0, 1.0, 0.0])
    tree = cKDTree(structure.xyz[atoms])

    # Bisect the push distance rather than solving it: the closest approach between two rigid
    # atom clouds along a line is monotone once they have separated, and 60 halvings is exact to
    # far below the file's own 0.001 A precision.
    low, high = 0.0, 500.0
    for _ in range(60):
        middle = 0.5 * (low + high)
        if tree.query(flipped[atoms] + direction * middle)[0].min() < TARGET_APPROACH:
            low = middle
        else:
            high = middle
    approach = float(tree.query(flipped[atoms] + direction * high)[0].min())

    def rewrite(coords: np.ndarray, chain_id: str) -> list[str]:
        return [
            line[:21] + chain_id + line[22:30] + "{:8.3f}{:8.3f}{:8.3f}".format(*coords[i]) + line[54:]
            for i, line in enumerate(lines)
        ]

    with DESTINATION.open("w") as handle:
        handle.write(
            f"REMARK   DODO test fixture: two copies of dnmt3a.pdb, chain B rotated 180 deg\n"
            f"REMARK   about x and translated {high:.3f} A along +y so the two C-terminal\n"
            f"REMARK   folded domains approach to {approach:.3f} A. See scripts/{Path(__file__).name}.\n"
        )
        handle.writelines(rewrite(structure.xyz, "A"))
        handle.write("TER\n")
        handle.writelines(rewrite(flipped + direction * high, "B"))
        handle.write("TER\nEND\n")

    print(f"wrote {DESTINATION} (translation {high:.3f} A, approach {approach:.3f} A)")
    print(read_structure(DESTINATION))


if __name__ == "__main__":
    main()
