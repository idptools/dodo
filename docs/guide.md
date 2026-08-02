# Using DODO

## The command line

One command, with subcommands.

```bash
# Rebuild a local structure
dodo rebuild AF-P04637-F1-model_v6.pdb -o p53_dodo.pdb

# Download from the AlphaFold database and rebuild
dodo fetch P04637 -o p53_dodo.pdb
dodo fetch "human p53" -o p53_dodo.pdb        # name resolution needs [lookup]

# Build a disordered region from sequence alone
dodo sequence GRNQNGGGYQNYNNQGYQGHGGQHQNNYNQYPCNYFGPGYNN -o my_idr.pdb

# Report what DODO thinks your structure looks like, without rebuilding
dodo regions AF-P04637-F1-model_v6.pdb

# Check any structure's geometry
dodo validate p53_dodo.pdb
```

Exit status is `0` on success, `2` if some regions could not be built, `1` on error.

### Shared flags

| Flag | Default | Meaning |
|---|---|---|
| `-o`, `--out` | *required* | Output PDB path |
| `-m`, `--mode` | `predicted` | Multiplier on the predicted end-to-end distance |
| `-n`, `--models` | `1` | Number of conformers |
| `-e`, `--engine` | `walk` | `walk` or `starling` |
| `-s`, `--strategy` | `auto` | Region identification: `auto`, `density`, `contact`, `plddt`, `metapredict` |
| `--seed` | none | Makes output bit-identical |
| `--ca-only` | off | Alpha carbons only, folded domains included |
| `-b`, `--annotate-regions` | off | Encode region type in the B-factor column, for colouring |
| `--no-conect` | off | Omit CONECT records — **not recommended**, see below |
| `-q`, `--quiet` | off | Suppress the per-region report |

### Build modes

Each mode is a multiplier on the ALBATROSS-predicted end-to-end distance, so it means the same
thing at any chain length.

| Mode | Multiplier |
|---|---|
| `super_compact` | 0.4× |
| `compact` | 0.7× |
| `normal`, `predicted` | 1.0× |
| `expanded` | 1.3× |
| `super_expanded` | 1.6× |
| `max_expansion` | 2.0× |

If a mode asks for more than the chain can physically span, DODO builds to 95% of the contour
length instead and warns, saying which of the two reasons applies — the request exceeds the fully
extended length, or it is reachable only by straightening so completely that no conformation
remains.

:::{warning}
These modes meant something different in 1.x, where they were Ångströms per residue. The same mode
name gives a different structure, and short regions come out **larger** than they did. See the
changelog.
:::

### Why CONECT records matter

CA–CA spacing is 3.81 Å, past the automatic bond-detection cutoff in both VMD and PyMOL. Without
CONECT records a rebuilt region renders as a cloud of disconnected dots. This is not cosmetic
polish; its absence defeats the tool. They are written by default — leave them on.

## The Python API

```python
import dodo

report = dodo.rebuild("model.pdb", mode="expanded", n_models=10, seed=0)
print(report.summary())
dodo.write_pdb(report.models, "out.pdb")
```

`rebuild()` returns a report rather than writing a file, so you can inspect what happened before
committing to output:

```python
report.ok            # True if every region in every model was rebuilt
report.models        # list[Structure] -- the conformers
report.failures      # regions that could not be built, each with a reason
report.assignments   # what DODO decided was folded vs disordered, with the evidence
report.outcomes      # per region per model: target, achieved dimension, or why it failed
report.placements    # where every folded domain ended up
```

From sequence alone:

```python
report = dodo.build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", n_models=5, seed=0)
```

### Working with the stages directly

Every stage is separately usable, which is the point of the rewrite.

```python
from dodo.io import read_structure, write_pdb
from dodo.regions import assign_regions, Strategy
from dodo.construct import target_dimensions

structure = read_structure("model.cif")          # PDB or mmCIF, gzip fine
assignment = assign_regions(structure, strategy=Strategy.PLDDT)[0]
print(assignment.describe())                     # chain A: IDR 1-31; FD 32-290; ...
print(assignment.score, assignment.threshold)    # the evidence behind every boundary

target = target_dimensions("GSGSGSGS" * 18, mode="compact")
print(target)   # 62.6 A (compact, 0.7x of 89.4 A) via albatross over 144 residues
```

### Overriding the regions

If you disagree with a boundary, say so, and DODO will build exactly what you asked for:

```python
import dodo

structure = dodo.read_structure("model.pdb")
dodo.assign_regions_from_spec(structure, {"A": [("idr", 1, 40), ("folded", 41, 290)]})
report = dodo.rebuild(structure, strategy="preset")
```

`strategy="preset"` means *identify nothing, build what is already there*. Bounds are **1-based
inclusive**, matching what you read off a PDB file. Overlaps, gaps and out-of-range bounds are
rejected with an explanation.

You can equally start from DODO's own answer and adjust it:

```python
structure = dodo.read_structure("model.pdb")
dodo.assign_regions(structure)                       # DODO's call, with its evidence
structure.chains[0].domains[1].span                  # inspect, then edit as you like
report = dodo.rebuild(structure, strategy="preset")   # build your version
```

This replaces 1.x's `regions_dict=`, which took a separate description of the structure alongside
the real one, so the two could disagree.

### The structure object

`Structure` is a struct-of-arrays type: one coordinate array, with `Domain` and `Chain` as
zero-copy views into it. Indices are 0-based positional throughout; PDB numbering is carried as
data and never used as an index. Spans are half-open `[start, stop)`.

```python
structure.xyz            # (n_atoms, 3)
structure.ca_xyz         # (n_residues, 3)
structure.sequence       # one-letter string
structure.domains        # every Domain across every chain
domain.xyz               # a zero-copy view; writing through it moves the domain
```

## Reproducibility

Pass `seed=` or `--seed` and the output is bit-identical across runs. Omit it and each run differs.
Both matter: the first for a figure you need to regenerate, the second for sampling.
