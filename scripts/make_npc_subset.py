"""Extract a small, committable slice of the 7R5J nuclear-pore assembly.

``tests/data/complexes_missing_residues/7R5J-assembly1.cif`` is 597 MB and 4.9 million atoms --
useful for a scale check, far too large to commit or to iterate on. This pulls five chains out
of it into a minimal mmCIF that still exercises everything the missing-residue path has to get
right, and that the FASTA in the same directory covers:

* ``00`` RanBP2 -- 756 residues modelled of 3,224, so 2,468 are missing, almost all of them one
  enormous C-terminal disordered region;
* ``A0`` Nup93 -- one residue missing, the trivial case;
* ``E0`` Ndc1 -- an internal gap as well as a terminal one;
* ``F0`` Nup35 -- a long N-terminal tail;
* ``W0`` Nup88 -- several folded domains with observed linkers between them.

Only the ``_atom_site`` loop is kept: everything DODO needs for this path is in the coordinates
and the FASTA, and copying the 900,000 lines of header would defeat the point.

    python scripts/make_npc_subset.py
"""

from __future__ import annotations

from pathlib import Path

#: Author chain ids to keep. One copy of each; the assembly repeats every chain eight times.
KEEP = frozenset({"00", "A0", "E0", "F0", "W0"})

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data" / "complexes_missing_residues"
SOURCE = DATA / "7R5J-assembly1.cif"
DESTINATION = DATA / "7R5J_subset.cif"


def main() -> None:
    if not SOURCE.exists():
        raise SystemExit(
            f"{SOURCE} is not present. It is a 597 MB download from the RCSB "
            f"(https://www.rcsb.org/structure/7R5J, biological assembly 1) and is deliberately "
            f"not committed."
        )

    header: list[str] = []
    atoms: list[str] = []
    in_atom_site = False
    with SOURCE.open() as handle:
        for line in handle:
            if line.startswith("_atom_site."):
                in_atom_site = True
                header.append(line)
            elif in_atom_site and line.startswith(("ATOM", "HETATM")):
                fields = line.split()
                # auth_asym_id is the third field from the end: ..., auth_comp_id, auth_asym_id,
                # auth_atom_id, model_num. Reading it positionally is safe because the header
                # above pins the column order for this file.
                if len(fields) > 20 and fields[-3] in KEEP:
                    atoms.append(line)

    if not atoms:
        raise SystemExit("no atoms matched; has the column order changed?")

    with DESTINATION.open("w") as out:
        out.write("data_7R5J_subset\n#\nloop_\n")
        out.writelines(header)
        out.writelines(atoms)
        out.write("#\n")
    print(f"wrote {DESTINATION} ({len(atoms)} atoms from chains {sorted(KEEP)})")


if __name__ == "__main__":
    main()
