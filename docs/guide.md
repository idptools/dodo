# Using DODO

## The command line

One command, with subcommands.

```bash
# Rebuild a local structure
dodo rebuild AF-P04637-F1-model_v6.pdb -o p53_dodo.pdb

# Download from the AlphaFold database and rebuild
dodo fetch P04637 -o p53_dodo.pdb
dodo fetch "human p53" -o p53_dodo.pdb        # names work out of the box

# Build a disordered region from sequence alone
dodo sequence GRNQNGGGYQNYNNQGYQGHGGQHQNNYNQYPCNYFGPGYNN -o my_idr.pdb

# Report what DODO thinks your structure looks like, without rebuilding
dodo regions AF-P04637-F1-model_v6.pdb

# Check any structure's geometry
dodo validate p53_dodo.pdb
```

Exit status is `0` on success, `2` if a region of 10 residues or more could not be rebuilt, `1` on
error. A shorter region that could not be rebuilt is reported and left as it arrived, and does not
fail the run — a few residues of AlphaFold geometry do not spoil a figure, whereas a long extended
region does.

### Shared flags

| Flag | Default | Meaning |
|---|---|---|
| `-o`, `--out` | *required* | Output PDB path |
| `-m`, `--mode` | `predicted` | Multiplier on the predicted end-to-end distance |
| `-n`, `--models` | `1` | Number of conformers |
| `-s`, `--strategy` | `auto` | Region identification: `auto`, `density`, `contact`, `plddt` |
| `--seed` | none | Makes output bit-identical |
| `--backbone` | off | Also place N, C and O on the rebuilt regions — see below |
| `--ca-only` | off | Alpha carbons only, folded domains included |
| `-b`, `--annotate-regions` | off | Encode region type in the B-factor column, for colouring |
| `--no-conect` | off | Omit CONECT records — **not recommended**, see below |
| `-q`, `--quiet` | off | Suppress the per-region report and the progress bar |

### Caching

DODO caches two things per user, both on by default, and `--no-cache` (or `DODO_NO_CACHE=1`) opts
out of both:

- **ALBATROSS predictions**, keyed by sequence hash and sparrow version. About 116 bytes each. This
  one is worth more than it looks: the first prediction in a process imports sparrow, parrot and
  torch, which measures **1.57 s against 0.17 s of actual rebuilding** for a 912-residue structure.
  A cache hit returns before sparrow is imported at all, so a repeat run never loads torch.
- **Downloaded structures.** Mean 0.25 MB and largest 1.51 MB, measured over 259 cached AlphaFold
  models. The path is printed on every fetch.

Neither cache can change a result: predictions are deterministic per sequence, and the key includes
the sparrow version so an upgrade that changes the network invalidates them rather than serving
stale values.

A progress bar is shown on stderr when stderr is a terminal, weighted by residues rather than by
region, because region lengths span two orders of magnitude and a per-region bar sits still through
the slowest part. It is suppressed automatically when output is piped, and by `-q`. From the Python
API, pass `progress=True` or `progress=False` to override.

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

A mode is a multiplier, but **no mode can exceed what the chain can physically span.** A region of
`N` residues has `N - 1` virtual CA–CA bonds of 3.81 Å, so its end-to-end distance can never exceed
`(N - 1) × 3.81 Å` — the fully extended contour length. Asking for more would mean breaking bonds.

DODO caps every target at **95% of the contour length** and warns when it has to. The 5% margin is
not timidity: a chain at exactly its contour length is a straight rod with one conformation, leaving
the sampler no freedom to avoid a clash.

This bites `max_expansion` on short regions, because the prediction it multiplies grows as roughly
`N^0.55` while the ceiling grows as `N`. Measured, with `*` marking a capped request:

| Region length | 95% ceiling | generic IDR | poly-E | poly-P | poly-G |
|---|---|---|---|---|---|
| 8 residues | 25.3 Å | 34 Å `*` | 39 Å `*` | 48 Å `*` | 30 Å `*` |
| 10 residues | 32.6 Å | 38 Å `*` | 46 Å `*` | 55 Å `*` | 34 Å `*` |
| 15 residues | 50.7 Å | 50 Å | 64 Å `*` | 72 Å `*` | 42 Å |
| 20 residues | 68.8 Å | 57 Å | 80 Å `*` | 86 Å `*` | 49 Å |
| 30 residues | 105.0 Å | 71 Å | 109 Å `*` | 108 Å `*` | 61 Å |
| 50 residues | 177.4 Å | 97 Å | 155 Å | 143 Å | 78 Å |

Where the cap bites depends on composition, because the prediction does: an expanded sequence like
poly-proline hits the ceiling out to about 30 residues, while poly-glycine never does. So on a short
region `max_expansion`, `super_expanded` and `expanded` can all converge on the same capped answer —
there is no more chain to give. The warning names which of the two reasons applied: the request
exceeded the contour length outright, or it was reachable only by straightening so completely that
no conformation remained.

:::{warning}
These modes meant something different in 1.x, where they were Ångströms per residue. The same mode
name gives a different structure, and short regions come out **larger** than they did. See the
changelog.
:::

### Backbone reconstruction (`--backbone`)

DODO places alpha carbons. `--backbone` additionally infers **N, C and O** for the regions it
rebuilt, from the alpha carbons alone. It is **opt-in** while it earns confidence; every bond it
writes inside a rebuilt region is exact; the limits are at the seams and in a handful of marginal
steric contacts, both described below.

The idea is that a CA-only trace is far more informative than it looks. Four consecutive alpha
carbons define a pseudo-dihedral, and that angle largely determines where the intervening peptide
unit's N and C sit. DODO bins that angle at 20° and looks the offsets up in a table measured from
100 frames of all-atom IDR simulation, then closes each carbon onto the nitrogen already placed so
that every bond inside a peptide unit is exact by construction. A second pass rotates each unit
about its CA–CA axis — the single degree of freedom it has left — to improve bond angles, steric
clashes and φ/ψ plausibility together.

Accuracy, measured properly: the lookup table was rebuilt from 80 simulation frames and tested on
the 20 it had never seen, over 3,643 residues whose alpha carbons were the only input.

| Atom | Mean error | Median | 90th percentile |
|---|---|---|---|
| N | 0.16 Å | 0.12 Å | 0.30 Å |
| C | 0.22 Å | 0.14 Å | 0.41 Å |
| O | 0.63 Å | 0.40 Å | 1.27 Å |

Rebuilding the table on 80 frames instead of 100 changed the result by 0.001 Å, so this is a
generalization figure rather than a fit: the shipped table scores 0.163 / 0.215 / 0.626 on the same
residues.

Oxygen is the loose one, and unavoidably: its position turns on ψ, which a CA trace constrains only
weakly. The N–CA–C bond angle comes out at 3.4° spread against 2.8° in real backbone, and all four
bond lengths are exact by construction.

**Folded domains are never affected.** They keep every atom they arrived with, `--backbone` or not.
Only the regions DODO generated gain atoms, which is what makes this additive rather than a
rewrite — and `--ca-only` remains the way to strip folded domains down.

**The seams are strained, and cannot be made exact.** Where a rebuilt region meets a folded domain,
the domain's existing nitrogen still points toward where the region ran in AlphaFold's original
model — DODO moved the domain rigidly and redrew the region beneath it, and folded-domain atoms are
not DODO's to move. A peptide unit reaches at most 2.854 Å from an alpha carbon to the nitrogen it
bonds to; measured across 17 seams in three structures, a rebuilt alpha carbon sits **3.2–4.5 Å**
away. The bond is therefore not merely hard to get right, it is geometrically unsatisfiable.

DODO aims the carbon at the nitrogen and leaves the bond long — around 2.2 Å against an ideal
1.33 — rather than writing two atoms into the same space, and reports every seam where it does.
Measured over three structures at five seeds each, that is 4 such bonds on dnmt3a (of 911 residues),
6 on arf19 and 10 on p300, identically at every seed.

Where the carbon and its oxygen go is **checked, not merely constructed**, and that distinction was
earned the hard way. Three successive versions placed them by geometry alone and each was blind to
one more neighbour: aiming the carbon at the anchor's nitrogen made the two collinear and threw the
carbonyl onto an arbitrary axis (O 0.6 Å from its own alpha carbon); holding the residue's own
N–CA–C angle instead left the oxygen 0.975 Å from that nitrogen; putting the oxygen trans to the
nitrogen left it 0.96 Å from the anchor's CB. The placement now sweeps the one free azimuth and
rejects any position that collides. Over the same 15 runs, `--backbone` introduces **zero** atom
pairs closer than the 1.00 Å floor below which no real bond exists, and **zero** damage to any
residue's own internal geometry.

:::{note}
The obvious fix is to leave the seam residue un-rebuilt so its alpha carbon is input geometry, and
it very nearly works — that distance becomes 2.45–2.52 Å at all 17 seams, comfortably reachable.
What defeats it is the side chain. Closing an exact bond onto that residue means re-placing its
nitrogen, and its side chain was built around where that nitrogen used to be, so the new one is
driven into the residue's own CB (measured 1.405 Å against a correct 2.45) and, for proline, snaps
the ring bond to CD (3.444 Å against 1.47). Trading a strained backbone bond for a broken proline
ring is not a trade worth making, so DODO does not.

Doing this properly means constraining the walk so its closing alpha carbon lands within peptide
reach of the anchor's nitrogen. Measured on ten seams, 5 have 20–43% of the closure circle in
reach and 5 have none, so it needs constraints reaching further back than the closing step. That is
future work, and it is what would make this the default.
:::

**It can introduce marginal steric contacts.** Alpha carbons are placed 3.20 Å apart and the trace
cannot self-intersect, but atoms hung off that trace can — a carbonyl oxygen points away from the
chain by design. Refinement scores each unit against the folded domains, against the rest of its own
region, and against the regions already placed before it, and the azimuth it can turn is a single
degree of freedom, so it cannot always find a clean answer. Measured against the CA-only output of
the same structures:

| Structure | CA-only | `--backbone` | Worst new overlap |
|---|---|---|---|
| arf19 | 0 | 1–15 | 0.48 Å |
| dnmt3a | 0 | 1–3 | 0.11 Å |
| p300 | 2 | 13–31 | 0.47 Å |

Ranges are three seeds, and the spread is the honest answer: the count turns on whether a conformer
happens to fold a hairpin back on itself, not on anything the placement does differently. Every
contact is marginal — under half an Ångström into the van der Waals limit, and nothing near
the 1.00 Å floor below which no real bond exists. This is the same trade DODO already makes for
anchor exemption: minimise, report, do not pretend to eliminate.

**φ/ψ are IDR-like, not folded-protein-like, and that is correct.** Measured over 19,602 pairs of
real simulated IDR backbone, the distribution is dominated by extended and PPII geometry with a
median near φ −107°, ψ +122° — not the α-helical basin a folded-structure intuition expects. A
Ramachandran plot of a rebuilt IDR that looked like a folded protein's would be evidence of a bug.

Building from sequence is the clean case — no folded domains, so no seams at all, and the output has
no bond violations and no impossible contacts whatsoever:

```bash
dodo sequence GRNQNGGGYQNYNNQGYQGHGGQHQNNYNQYPCNYFGPGYNN --backbone -o my_idr.pdb
```

### Why CONECT records matter, and what to do in VMD

CA–CA spacing is 3.81 Å, well past the distance cutoff either VMD or PyMOL uses to infer bonds. So a
CA-only region has no bonds a viewer can find for itself, and without CONECT records it renders as a
cloud of disconnected dots. They are written by default — leave them on.

DODO's CONECT output is complete: measured on p300, all 1,520 consecutive CA–CA pairs inside rebuilt
regions carry a record, in the fixed columns the format specifies, and `dodo validate` checks this.
**PyMOL honours them and draws the trace.**

**VMD frequently does not.** VMD infers bonds by distance when it loads a PDB and its handling of
CONECT is unreliable, so a rebuilt region can still come up as dots however correct the file is. That
is not something DODO can fix from this side. The fix is to use a representation that does not
consult the bond list at all:

```
# In the VMD representation dialog, set Drawing Method to:
Trace        # connects consecutive alpha carbons directly
Tube         # same, drawn thicker
```

`Trace` and `Tube` follow residue order rather than bonds, which is exactly right for a CA-only
model. If you need the bonds themselves — for a selection or an analysis rather than a picture —
`topo readvarxyz` / the `topotools` plugin can add them explicitly, or convert to a format that
carries topology, such as PSF.

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
report.ok            # True if every region that matters was rebuilt (see below)
report.models        # list[Structure] -- the conformers
report.failures      # regions that could not be built, each with a reason
report.assignments   # what DODO decided was folded vs disordered, with the evidence
report.outcomes      # per region per model: target, achieved dimension, or why it failed
report.placements    # where every folded domain ended up
```

`report.ok` tolerates short regions. A region under 10 residues that could not be rebuilt is
reported but does not make the run unsuccessful, because a few residues left as AlphaFold drew them
do not look wrong — it is the long extended regions DODO exists to fix. The three relevant fields:

```python
report.failures            # every region not rebuilt, short ones included
report.blocking_failures   # the subset long enough to matter; report.ok is `not this`
report.tolerated_failures  # the short ones, reported and left alone
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
print(target)   # 55.9 A (compact, 0.7x of 79.8 A) via albatross over 144 residues
```

### Overriding the regions

If you disagree with a boundary, say so, and DODO will build exactly what you asked for:

```python
import dodo

structure = dodo.read_structure("model.pdb")
dodo.assign_regions_from_spec(structure, {"A": [("idr", 1, 40), ("folded", 41, 912)]})
report = dodo.rebuild(structure, strategy="preset")
```

`strategy="preset"` means *identify nothing, build what is already there*. Bounds are **1-based
inclusive**, matching what you read off a PDB file, and the regions must **tile the whole chain** —
every residue belongs to exactly one. Overlaps, gaps and out-of-range bounds are rejected with an
explanation naming the offending residue.

You can equally start from DODO's own answer and adjust it:

```python
import dodo

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
