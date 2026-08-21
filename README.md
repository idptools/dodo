DODO: re<ins>D</ins>esign c<ins>O</ins>mputationally generate<ins>D pr<ins>O</ins>teins
==============================

## What is DODO?

Protien structure predictors have revolutionized biology. However, many predicted 
structures have low-confidence regions that do not adopt any predicted secondary 
structure. DODO takes these regions and rebuilds them under the assumption that they
are disordered regions. At this time, the rebuilt regions include the backbone of the 
amino acids but do not include the side chains. This is an improvement over the first
version of DODO because we now have the full backbone instead of just alpha carbons. I'm 
hoping to implement all atom IDR rebuilds in the future. 

I want to be extremely clear that **the DODO-generated IDRs are not simulations** and should not be
treated as such. The main constraints are simply bond lengths, angles, and to avoid clashing. 
That's about it. 

By default, DODO identifies each disordered region and predicts its end-to-end distance from sequence
— with ALBATROSS ([sparrow](https://github.com/idptools/sparrow)).
You can also ask for regions more compact or more expanded than predicted. And you
can generate several conformers into one multi-model PDB, which in VMD looks *like* a
simulation trajectory. Once again, to be very clear, it is **not** equivalent to a simulation, but it's
nice for visualization.


![DODO_EXAMPLE](https://github.com/idptools/dodo/blob/main/images/DODOv2_example.png)

> **DODO 2.0 is a rewrite and it breaks the 1.x API.** The `build.pdb_from_name()` /
> `pdb_from_pdb()` / `pdb_from_sequence()` functions and the three `pdb-from-*` console
> scripts are gone, replaced by `dodo.rebuild()` / `dodo.build_from_sequence()` and a single
> `dodo` command. See [What changed, and migrating from 1.x](#what-changed-and-migrating-from-1x).

## Installation

Requires Python 3.10 or newer.

```bash
pip install git+https://github.com/idptools/dodo.git
```

The base install depends on **numpy, scipy, numba, tqdm, getSequence, and sparrow**.


## Command line

One command with subcommands.

```bash
# Rebuild a local structure
dodo rebuild AF-P04637-F1-model_v6.pdb -o p53_dodo.pdb

# Download from the AlphaFold database and rebuild
dodo fetch P04637 -o p53_dodo.pdb
dodo fetch "human p53" -o p53_dodo.pdb        # names work out of the box

# Build a disordered region from sequence alone
dodo sequence GRNQNGGGYQNYNNQGYQGHGGQHQNNYNQYPCNYFGPGYNN -o my_idr.pdb

# Just tell me what you think my structure looks like
dodo regions AF-P04637-F1-model_v6.pdb
```

### Shared flags

Shared by the commands that build — `rebuild`, `fetch`, `sequence`. `regions` takes only
`-s/--strategy` and `--scores`; `validate` takes only `--max-findings`, `--no-clashes` and
`--no-bonds`. Neither accepts `-o` or any of the build flags below.

| Flag | Default | Meaning |
|---|---|---|
| `-o`, `--out` | *required* | Output path. The format follows the extension: `.cif`/`.mmcif` write mmCIF, anything else writes PDB. See [multi-model output](#multi-model-output-is-now-a-spread-of-conformers) for which to use in which viewer |
| `-m`, `--mode` | `predicted` | Target dimension as a multiplier on the predicted end-to-end distance |
| `-n`, `--models` | `1` | Number of conformers. Folded domains are positioned once and held fixed across all models; only the disordered regions differ |
| `-s`, `--strategy` | `auto` | How to identify regions: `auto`, `density`, `contact`, `plddt` |
| `--seed` | none | Makes output reproducible |
| `--backbone` / `--no-backbone` | **on** | Place N, C and O on the rebuilt regions, inferred from the alpha carbons; `--no-backbone` writes alpha carbons only |
| `--ca-only` | off | Alpha carbons only, folded domains included |
| `-b`, `--annotate-regions` | off | Encode region type in the B-factor column, for colouring |
| `--no-conect` | off | Omit CONECT records — **not recommended**, see [below](#why-conect-records-matter) |
| `-q`, `--quiet` | off | Suppress the per-region report and the progress bar |
| `--cache-structures` (`fetch` only) | off | Keep the downloaded AlphaFold model on disk; otherwise it is deleted when the command exits |

### Caching

**Downloaded structures are not kept.** `dodo fetch` downloads to a temporary directory and
deletes it on exit. Pass `--cache-structures` if you want repeat fetches of the same accession to
skip the download. A model averages 0.25 MB and runs to 1.51 MB — small individually, but a
directory you never chose to grow, so keeping them is your decision rather than DODO's.

**ALBATROSS predictions are cached, and this one is on by default. However, we won't destroy your disk.** 
The reason is speed, not convenience: a cache miss has to import sparrow, which pulls in parrot and torch. 
Measured on a 912-residue structure, that import is **1.6 s against 0.17 s of actual rebuilding** — it 
dominates the run. Every CLI invocation is a fresh process, so without the cache you pay it every single
time:

| | full `dodo rebuild` |
|---|---|
| prediction cache warm | **0.36 s** |
| prediction cache disabled | **1.96 s** |

The disk cost is negligible. We measured this. An entry is a sequence hash
and a float — **116 bytes**. Rebuilding *the entire AlphaFold human proteome*, all 23,587
structures and every disordered region in them **totalled 1.07 MB**. The whole human proteome 
costs about one megabyte.

The cache lives at `~/.cache/dodo/<generation>/` (or the platform equivalent) and is keyed by
sequence hash *and* sparrow version, so upgrading sparrow invalidates old values rather than
silently serving predictions from a different network.

**It is capped at 10 MiB**, roughly ten times the whole human proteome, with the oldest entries
dropped first once it is reached. Nothing you do in normal use will approach that; the cap is
there so no workload, however unusual, can grow the file without bound.

**To turn it off:** `--no-cache` on the commands that predict — `rebuild`, `fetch`, `sequence` — or
`DODO_NO_CACHE=1` in the environment. Both disable prediction caching and structure caching
together. DODO still works exactly the same — it is only slower. `regions` and `validate` never call
the predictor, so they reject `--no-cache` rather than accept a flag that would do nothing.

Exit status is `0` on success, `2` if a region of 10 residues or more could not be rebuilt, `1` on
error. A shorter region that could not be rebuilt is reported and left as it arrived, and does not
fail the run — a few residues of AlphaFold geometry do not spoil a rebuild, whereas a long extended
region does.

## Python API

```python
import dodo

report = dodo.rebuild("AF-P04637-F1-model_v6.pdb", mode="expanded", n_models=10, seed=0)
print(report.summary())
dodo.write_pdb(report.models, "p53_dodo.pdb")
```

`rebuild()` returns a `RebuildReport` rather than writing a file, so you can inspect what
happened before committing to output:

```python
report.ok            # True if every region that matters was rebuilt (see below)
report.models        # list[Structure] -- the conformers
report.failures      # regions that could not be built, each with a reason
report.assignments   # what DODO decided was folded vs disordered, with the evidence
report.outcomes      # per region per model: target, achieved dimension, or why it failed
```

From sequence alone:

```python
report = dodo.build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", n_models=5, seed=0)
```

### Working with the pieces directly

Every stage is separately usable, which is the point of the rewrite:

```python
from dodo.io import read_structure, write_pdb
from dodo.regions import assign_regions, Strategy
from dodo.construct import target_dimensions

structure = read_structure("model.cif")          # PDB or mmCIF, gzip fine
assignment = assign_regions(structure, strategy=Strategy.PLDDT)[0]
print(assignment.describe())                     # one line per chain, as `dodo regions` prints
print(assignment.score, assignment.threshold)    # the evidence behind every boundary

target = target_dimensions("GSGSGSGS" * 18, mode="compact")
print(target)   # 55.9 A (compact, 0.7x of 79.8 A) via albatross over 144 residues
```

The core `Structure` is a struct-of-arrays type: one coordinate array, with `Domain` and
`Chain` as zero-copy views into it. Indices are 0-based positional throughout; PDB numbering is
carried as data and never used as an index. Spans are half-open `[start, stop)`.

## How it works

The order of these steps is not incidental — it's the algorithm.

1. **Identify** folded domains, IDRs, and loops (IDRs tethered at both ends inside a *single*
   folded domain), using the all-atom density metric.
2. **Predict** each IDR's end-to-end distance with ALBATROSS. Loops get no prediction — their
   span is already dictated by the folded-domain geometry they bridge.
3. **Reposition the folded domains.** Translate and rotate each one as a rigid body so that
   consecutive domains sit at the predicted end-to-end distance of the IDR between them.
4. **Rebuild loops.**
5. **Rebuild connecting IDRs** between folded domains.
6. **Rebuild terminal IDRs.**

Steps 4–6 run in that order because that is decreasing order of constraint: a loop is pinned at
both ends inside one domain, a connecting IDR is pinned between two domains that have already
been positioned for it, and a terminal IDR is free at one end and can go almost anywhere. 

**Folded-domain atoms are never rebuilt.** They come from AlphaFold (or AlphaFold3, or a
crystal structure, or wherever else) and are trusted. A domain only ever moves as a rigid body, and DODO checks
that: after each domain is repositioned.

## Build modes

Modes are **multipliers on the predicted end-to-end distance**:

| Mode | Factor |
|---|---|
| `super_compact` | 0.4× |
| `compact` | 0.7× |
| `normal` / `predicted` | 1.0× |
| `expanded` | 1.3× |
| `super_expanded` | 1.6× |
| `max_expansion` | 2.0× |

**This changed from 1.x.** v1 expressed modes as Ångströms *per residue* —
`normal` was 0.8 Å/residue. That's linear in chain length, but real IDR end-to-end distance
scales as roughly $N^{0.52}$, so a fixed per-residue multiplier can only agree with the
prediction at one length. `normal` gave 80 Å at N=100 (about right) and 400 Å at N=500, where
the prediction is nearer 190 Å — and the error grew without bound.
One consequence: `normal` and `predicted` are now synonyms, which they weren't in 1.x.

A mode is a multiplier, but **no mode can exceed what the chain can physically span.** A region of
`N` residues has `N - 1` virtual CA–CA bonds of 3.81 Å, so its end-to-end distance can never exceed
`(N - 1) × 3.81 Å` — the fully extended contour length. Asking for more would mean breaking bonds.

DODO caps every target at **95% of the contour length** and warns when it has to. The 5% margin is
not timidity: a chain at exactly its contour length is a straight rod with one conformation, leaving
the sampler no freedom to avoid a clash.

## Region identification

Four strategies behind one flag:

- **`density`** — DODO's original all-atom density metric: all-atom pairs within 8 Å per
  residue, thresholded at 480. **This is the method DODO was built and validated on**.
- **`contact`** — a CA-only alternative. Composition-free (every residue has exactly one CA)
  and invariant to whether side chains are modelled, which the density score is not. Useful for
  comparison and for input without side chains, but it is not the validated method.
- **`plddt`** — AlphaFold's own per-residue confidence, from the B-factor column. Cheap, but
  density is DODO's validated default, so this is an explicit opt-in.
- **`auto`** (default) — `density` for all-atom input, `contact` when the input has no side
  chains at all, CA-only or full-backbone alike (a pair count can't be compared against a
  threshold tuned on packed cores). It never picks pLDDT on its own, and it says so whenever it
  picks anything other than `density`; on all-atom input, the expected case, it stays quiet.

Every assignment keeps its score profile and threshold, so a boundary you disagree with can be
audited rather than just re-run with different numbers.

### Checking the call before you rebuild

`dodo regions` runs region identification and stops: it reads the structure, prints what it found,
and writes nothing. Nothing is predicted, so sparrow — and torch behind it — is never imported, and
the whole command takes about a third of a second on a 912-residue model. Speed is not really the
point; seeing what DODO will treat as disordered *before* it rebuilds anything is.

The files below are the regression fixtures in `tests/data/structures/`, so reproducing these lines
needs a git clone rather than a `pip install`. Output is shown as `#` comments.

```bash
dodo regions tests/data/structures/dnmt3a.pdb
# chain A: IDR 1-282; FD 283-432 (loops: 383-393); IDR 433-473; FD 474-912
```

One line per chain, tiling it with nothing left over. `383-393` is reported *inside* `FD 283-432`
rather than as a region of its own because it is a loop — an IDR tethered at both ends within that
one domain. Ranges are the input file's own residue numbers, so they can be read straight off
against the structure or typed into a viewer.

Region lines go to **stdout** and `note:` lines to **stderr**, which is what makes the command
scriptable: `2>/dev/null` leaves pasteable region lines, `2>&1` keeps the explanation with them.

```bash
dodo regions tests/data/structures/p300.pdb
# chain A: IDR 1-334; FD 335-417; IDR 418-568; FD 569-650; IDR 651-1048; FD 1049-1527; IDR 1528-1577; FD 1578-1831 (loops: 1710-1719); IDR 1832-2414
#   note: 1 short stretch (19 residues) scored as folded but came in under the 25-residue minimum for a folded domain, so it is treated as disordered
```

A note is where a tuned minimum changed what will be rebuilt: 25 residues is the shortest run that
may become a folded domain, so p300's 19-residue stretch is called disordered instead — and is
therefore rebuilt. The [guide](docs/guide.md#reading-the-notes) lists what each note means.

**The same structure under all four strategies** — the comparison the command exists for. The
`chain A: ` prefix is dropped from each row:

| `-s` | Cutoff | Regions in `dnmt3a.pdb` |
|---|---|---|
| `auto`, resolving to `density` here | 480 | `IDR 1-282; FD 283-432 (loops: 383-393); IDR 433-473; FD 474-912` |
| `contact` | 10 | `IDR 1-174; FD 175-221; IDR 222-279; FD 280-434 (loops: 383-393); IDR 435-467; FD 468-912` |
| `plddt` | 70 | `IDR 1-282; FD 283-430 (loops: 383-393); IDR 431-473; FD 474-912` |

Three different scores, so the cutoffs do not compare across strategies — only the boundaries do.
pLDDT lands within two residues of density here; `contact` does not, promoting a 47-residue stretch
(`175-221`) that density calls disordered. When they disagree, `--scores` prints the per-residue
profile behind whichever one you ran, as `residue<TAB>score` under the region line and on the same
numbering, which is the audit trail for a boundary you doubt.

Once you have picked one, the same `-s` goes to `rebuild`:

```bash
dodo rebuild tests/data/structures/dnmt3a.pdb -s plddt -o dnmt3a_dodo.pdb
```

`dodo regions` exits 0, or 1 if the structure cannot be read.

### Specifying the regions yourself

The metrics above are a convenience, not a requirement. If you already know where your domains are
— from a paper, a UniProt annotation, or your own judgement — you can state them and DODO will
build exactly that, running no metric at all.

Give each chain a list of `(kind, start, stop)` regions, then rebuild with `strategy="preset"`,
which means *identify nothing, build what is already there*:

```python
import dodo

structure = dodo.read_structure("model.pdb")
dodo.assign_regions_from_spec(structure, {
    "A": [
        ("idr",     1, 282),      # rebuilt
        ("folded", 283, 912),     # moved as a rigid body, never regenerated
    ],
})
report = dodo.rebuild(structure, strategy="preset", seed=0)
dodo.write_pdb(report.models, "my_regions.pdb")
```

Three rules, and the first is the one people trip over:

1. **The regions must tile the chain.** Every residue belongs to exactly one region — no gaps, no
   overlaps, and the first and last must reach the chain's ends. A gap is rejected, not silently
   filled.
2. **Bounds are inclusive residue numbers, exactly as the input file numbers them** — the same
   numbering `dodo regions` prints, so you can copy a boundary from one straight into the other.
   A residue distinguished by an insertion code is named as a string, `"10A"`.
3. **`folded` is a promise, not a request.** A folded region's atoms are never regenerated; it is
   only translated and rotated as a rigid body. Everything you mark `idr` is rebuilt from scratch.

A folded domain can carry **rebuildable loops** — stretches inside it that you do want regenerated
— as an optional fourth element. This is how you say *"283 to 912 is one domain, and the piece in
the middle is a loop within it, not a linker between two domains"*:

```python
import dodo

structure = dodo.read_structure("model.pdb")
assignment = dodo.assign_regions_from_spec(structure, {
    "A": [("idr", 1, 282), ("folded", 283, 912, [(433, 473)])],
})[0]
print(assignment.describe())
# chain A: IDR 1-282; FD 283-912 (loops: 433-473)
```

A loop must lie **strictly inside** its folded domain, since it is rebuilt between two fixed anchor
residues and needs one on each side. A stretch that reaches a domain boundary is a tail — make it
its own `idr` region instead.

Every rejection names the residue you typed and what was wrong with it, so a bad spec fails
immediately rather than surfacing later as unrelated geometry:

```python
import dodo
from dodo.exceptions import InvalidRegionError

structure = dodo.read_structure("model.pdb")
try:
    dodo.assign_regions_from_spec(structure, {"A": [("idr", 1, 282), ("folded", 300, 912)]})
except InvalidRegionError as error:
    print(error)
# Chain 'A' has unassigned residues 283-299.
```

You can equally start from DODO's own answer and adjust the part you disagree with, which keeps its
score profile and threshold as evidence for the boundaries you did not touch:

```python
import dodo
from dodo.structure import Span

structure = dodo.read_structure("model.pdb")
dodo.assign_regions(structure)                       # DODO's call, with its evidence
domain = structure.chains[0].domains[1]
print(domain.span)                                   # inspect
domain.span = Span(300, 400)                         # a Span is frozen: replace, don't mutate
report = dodo.rebuild(structure, strategy="preset")  # build your version
```

A spec, by contrast, carries no evidence — nothing was measured, so its score and threshold are
`nan` rather than a zero that would read as a measurement.

This is a Python-only path: `strategy="preset"` needs region objects, and there is no command-line
syntax for them, so `dodo rebuild -s preset` is not accepted. It replaces v1's `regions_dict=`,
which took a parallel stringly-typed description of the structure that could disagree with the real
one. The [guide](docs/guide.md#overriding-the-regions) covers multi-chain specs, the full error
list, and what a preset assignment does and does not carry.

## Multi-model output is now a spread of conformers

`-n 10` writes ten conformers as MODEL/ENDMDL frames. The folded domains are positioned once
and held fixed across every model, so a viewer flicking between frames sees only the disordered
regions move — if the domains were re-placed per model they would jump around and destroy the
illusion. One consequence worth knowing: a *connecting* IDR's end-to-end distance is fixed by
its anchors, so across models its path varies but its span cannot. A *terminal* IDR is free at
one end and does scatter.

**In 1.x they effectively didn't.** v1 placed folded domains once outside the model loop and
targeted only the *mean* predicted end-to-end distance, so every model shared one arrangement
and essentially one dimension. Measured on an early 2.0 engine, the coefficient of variation of
end-to-end distance across conformers was 0.006–0.045, where a freely-rotating chain with the same
bond length and CA–CA–CA angle window gives 0.35–0.48. Sixty models of a 200-residue IDR spanning
1.9 Å of extension is one conformation sampled sixty times.

A predicted end-to-end distance is the **mean of a distribution**, so each model now draws its
own target from the ideal-chain radial distribution $P(R) \propto R^2 e^{-3R^2/2\langle R^2\rangle}$,
whose own coefficient of variation is 0.42. The suite holds the achieved spread inside 0.20–0.60,
and the ensemble mean still matches the prediction.

Two details worth knowing:

- A **single** model still hits the target exactly. One conformer should be the predicted
  dimension, not a random draw around it.
- Only regions with a **free end** scatter. An interior region pinned between two folded
  domains has its end-to-end distance determined by the anchor separation — there's nothing to
  sample, and its models differ in path rather than in span.

### File format, and a VMD caveat

The output format follows the `-o` extension. `.pdb` writes the conformers as `MODEL`/`ENDMDL`
frames; `.cif`/`.mmcif` writes a single `_atom_site` loop keyed by `pdbx_PDB_model_num`, which is
the mmCIF way to carry several models. Both are correct, and PyMOL, ChimeraX and most tools read
the models from either.

**VMD is the exception for mmCIF.** Its mmCIF reader loads every model into a *single* frame, so a
multi-model `.cif` opens as all conformers overlaid at once — a hairball — instead of a browsable
trajectory. This is a VMD limitation, not a problem with the file: VMD does the same to the PDB's
own deposited NMR ensembles (e.g. 1L2Y's 38 models load as one frame), while other readers see the
models correctly. So:

- **Opening a multi-model ensemble in VMD → write PDB** (`-o out.pdb`). VMD reads its
  `MODEL`/`ENDMDL` frames as a proper trajectory you can step through.
- **PyMOL or ChimeraX → either format works.** mmCIF is the better choice for very large single
  models, where PDB's fixed-width columns run out of room.

A *single*-model `.cif` opens fine in VMD; only the multi-model case is affected. `dodo` prints a
reminder to this effect whenever it writes a multi-model mmCIF.

NOTE: VMD 2.0 is in the works at this time. We do not yet know if it will support multiple models
in a .cif file. However, the VMD limitation is confirmed for VMD < 2.0.

## Reproducibility

Everything stochastic takes a seed. Same seed, bit-identical output:

```bash
dodo rebuild model.pdb -o out.pdb --seed 42
```

v1 had no seed anywhere, so a stochastic structure builder wasn't reproducible and couldn't be
regression-tested at all.

## Why CONECT records matter

DODO writes CONECT records by default, and you should leave them on. They matter most under
`--no-backbone`: CA–CA spacing is 3.81 Å, past the automatic bond-detection cutoff in both VMD and
PyMOL, so without CONECT an alpha-carbon-only region renders as a cloud of disconnected dots. With
the backbone on (the default) a viewer can find the real N–CA and CA–C bonds for itself, but the
records still declare the chain explicitly rather than leaving it to each viewer's heuristics.

## Current limitations

1. **~~Conformations aren't scientifically useful.~~**
   They are improved over V1, which used a random walk. 
   Now, the walk is now a self-avoiding, angle-constrained growth walk that hits
   the predicted dimensions and produces a spread of distinct conformers across models — but it is
   a geometric sampler, not a force field, so that spread is not a thermodynamic ensemble and is
   not a substitute for a simulation.

2. **~~Rebuilt IDRs do not contain side chains.~~** 
   However, rebuilt IDRs **do now contain the entire backbone!** This is a major improvement over V1.
   N, CA, C, O backbone by default, and `--no-backbone` returns the alpha-carbon-only output if you
   want it. Note the backbone-carbon-only limitation does NOT apply to
   the regions DODO leaves alone, which keep every atom.

3. **Assembly rebuilding is not implemented.** Multi-chain structures are read and written
   correctly, and regions are identified per chain, but rebuilding unmodelled regions of an EM
   assembly against the deposited sequence isn't wired up yet.

4. **Cis-peptide bonds are not modelled.** DODO builds every virtual CA–CA bond at 3.81 Å, 
   the trans value. A cis peptide — in practice almost always X–Pro — sits near 2.9 Å, and 
   DODO does not produce one.

5. **Anchor alpha carbons are exempt from clashing.** 
   Rebuilding a region means attaching it to the fixed residues on either side, and that requires
   exempting those *anchors* from clash checking to some degree. This is unavoidable, so it is worth
   being plain about it. The anchors' **alpha carbons** are always exempt. The region is bonded to 
   them; treating them as obstacles would make every valid attachment register as a clash and there 
   would be nothing to build. There is no version of the algorithm without this.

## Improvements over DODO V1

1. **~~Unusual-bond warnings in VMD.~~** Addressed: correct CONECT records, correct atom-name
   columns, and the element column written. (v1 right-justified atom names from column 13 and
   omitted the element, with the result that MDTraj read its CA-only output as 912 *calcium*
   atoms.)

2. **~~Some visualization modes don't work in VMD; tube and trace fail.~~** Should be addressed,
   but please report if you still see it.



### Backbone reconstruction (`--backbone` / `--no-backbone`)

By default DODO places **N, C and O** on the rebuilt regions, inferred from the alpha carbons alone
— consecutive alpha carbons largely determine where the peptide unit between them sits, and DODO
looks that up in a table measured from 100 frames of all-atom IDR simulation, then settles each
unit's one remaining degree of freedom against bond angles, clashes and φ/ψ together.

Every bond length inside a rebuilt region is exact by construction. Side chains are still not built.

The one place it cannot be exact is the **seams**. Where a rebuilt region meets a folded domain,
that domain's nitrogen still points toward where the region ran in AlphaFold's model, and
folded-domain atoms are not DODO's to move. A peptide unit reaches at most 2.854 Å from an alpha
carbon to the nitrogen it bonds to; a rebuilt alpha carbon measures well beyond that from it. The
bond is unsatisfiable, so DODO aims the atom as close as the residue's own N–CA–C angle allows,
leaves the bond long — measured 2.6–3.7 Å against an ideal 1.33 (mean ~3.0) — and labels and
reports it (on `RebuildReport.backbone_seams` and in the run summary). Nothing impossible is
written: measured over three structures at three seeds, backbone placement introduces zero atom
pairs closer than the 1.00 Å floor below which no real bond exists.

Building from sequence has no seams and comes out completely clean. See the
[user guide](https://dodo.readthedocs.io) for the full accounting, including why leaving the seam
residue un-rebuilt does not fix it.

## What changed, and migrating from 1.x

**2.0 breaks the 1.x API deliberately.** This section is enough to port a script. The exhaustive
list — every behaviour change with its measurements — is in [CHANGELOG.md](CHANGELOG.md).

### The three that bite hardest

1. **Build modes mean something different, so your figures will change.** 1.x modes were
   Ångströms *per residue*; 2.0 modes are *multipliers on the predicted end-to-end distance*.
   Short regions come out **larger** than 1.x made them, long ones smaller, crossing over near 80
   residues. Regenerating a 1.x figure will not reproduce it — and 2.0's version is the correct
   one, because real IDR dimensions scale as roughly N^0.55, not linearly.

   | Mode | 1.x (Å per residue) | 2.0 (× predicted) |
   |---|---|---|
   | `super_compact` | 0.3 | 0.4× |
   | `compact` | 0.55 | 0.7× |
   | `normal` | 0.8 | 1.0× |
   | `predicted` | *(used the prediction)* | 1.0× |
   | `expanded` | 1.05 | 1.3× |
   | `super_expanded` | 1.3 | 1.6× |
   | `max_expansion` | 1.65 | 2.0× |

   Note `normal` and `predicted` were **different** in 1.x and are exact synonyms in 2.0.

2. **The entry points are gone**, replaced by one command and two functions.
3. **Functions return a report** instead of returning `None` and printing, so a script can finally
   tell whether a region was actually rebuilt.

### Commands

| 1.x | 2.0 |
|---|---|
| `pdb-from-pdb in.pdb -o out.pdb` | `dodo rebuild in.pdb -o out.pdb` |
| `pdb-from-name "human p53" -o out.pdb` | `dodo fetch "human p53" -o out.pdb` |
| `pdb-from-sequence SEQ -o out.pdb` | `dodo sequence SEQ -o out.pdb` |
| *(none)* | `dodo regions in.pdb` — what DODO thinks your structure is, without rebuilding |
| *(none)* | `dodo validate out.pdb` — bond lengths, clashes and CONECT records |

Multi-word protein names now need quoting: `dodo fetch "human p53"`.

### Flags

| 1.x | 2.0 |
|---|---|
| `-o`, `--out_path` | `-o`, `--out` |
| `-m`, `--mode` | `-m`, `--mode` (same names, **different meaning** — see above) |
| `-n`, `--num_models` | `-n`, `--models` |
| `-c`, `--no_CONECT_lines` | `--no-conect` (no short form) |
| `-f`, `--no_FD_atoms` | `--ca-only` |
| `-b`, `--beta_for_FD_IDR` | `-b`, `--annotate-regions` |
| `-s`, `--silent` | `-q`, `--quiet` — **`-s` now means `--strategy`** |
| `-u`, `--use_metapredict` | removed; `density` is the default and now takes 7 ms |
| `-apr`, `--attempts_per_region` | removed; fixed at 40 internally |
| `-apc`, `--attempts_per_coord` | removed |
| `-apr`/`-api` on `pdb-from-sequence` | removed (note `-apr` meant *per residue* here and *per region* elsewhere) |
| `-l`, `--linear_placement` | removed |
| `-j`, `--just_fds` | removed |
| *(none)* | `--seed`, `--backbone`/`--no-backbone`, `--no-cache`, `--cache-structures` |

### Python API

| 1.x | 2.0 |
|---|---|
| `build.pdb_from_pdb(path, out_path=...)` | `dodo.rebuild(path)` then `dodo.write_pdb(...)` |
| `build.pdb_from_name(name, out_path=...)` | `dodo.fetch_alphafold(acc)` then `dodo.rebuild(...)` |
| `build.pdb_from_sequence(seq, out_path=...)` | `dodo.build_from_sequence(seq)` |
| `num_models=N` | `n_models=N` |
| `CONECT_lines=False` | `write_pdb(..., conect=False)` |
| `include_FD_atoms=False` | `write_pdb(..., ca_only=True)` |
| `beta_for_FD_IDR=True` | `write_pdb(..., annotate_regions=True)` |
| `verbose=False` | `progress=False` (and the report is returned, not printed) |
| `regions_dict={...}` | `assign_regions_from_spec(...)` with `strategy='preset'` |
| `use_metapredict=True` | removed; see above |
| `graph=True` | removed — it plotted with matplotlib, which is no longer a dependency |
| `end_coord=(x,y,z)` on `pdb_from_sequence` | removed |
| `just_fds=True` | removed |
| functions returned `None` and printed | functions return a `RebuildReport` |
| `except dodoException` | `except DodoError` |

```python
# 1.x
import dodo
dodo.build.pdb_from_pdb("in.pdb", out_path="out.pdb", mode="expanded", num_models=10)

# 2.0
import dodo
report = dodo.rebuild("in.pdb", mode="expanded", n_models=10, seed=0)
print(report.summary())          # what happened, per region
dodo.write_pdb(report.models, "out.pdb")
```

### Other behaviour worth knowing

- **A multi-model run is now an ensemble.** 1.x placed folded domains once outside the model
  loop, so every model shared one arrangement *and* essentially one end-to-end distance. Each
  model now draws its own dimensions; folded domains are still positioned once and held fixed, so
  only the disordered regions move between frames.
- **`--seed` exists.** 1.x had no seed anywhere, so a stochastic builder was not reproducible.
- **Rebuilt regions get a backbone by default.** N, C and O are placed on them; `--no-backbone`
  returns the alpha-carbon-only output 1.x produced. Side chains are still not built.
- **Exit status means something:** `0` success, `2` a region of 10+ residues could not be rebuilt,
  `1` error.
- **Installation is lighter.** metapredict, matplotlib and cython are gone from the dependencies.

## Copyright

Copyright (c) 2023-2026, Ryan Emenecker — Holehouse Lab

#### Acknowledgements

Originally based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms)
version 1.1.
