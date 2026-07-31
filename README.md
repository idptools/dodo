DODO: re<ins>D</ins>esign AlphaF<ins>O</ins>ld <ins>D</ins>isordered regi<ins>O</ins>ns
==============================

## What is DODO?

DODO takes a predicted protein structure, works out which parts are folded domains and which
are intrinsically disordered regions, and rebuilds the disordered parts so they adopt realistic
polymer dimensions instead of AlphaFold's characteristic extended spaghetti.

To be clear: the work DeepMind did on AlphaFold is **amazing**, and none of this takes away
from it in *any way*. AlphaFold simply isn't trying to represent a disordered region as an
ensemble — an IDR has no single structure to predict. But for figures and presentations it's
useful to have those regions *look* like what they are. DODO does that.

It identifies each disordered region, predicts its end-to-end distance from sequence with
ALBATROSS ([sparrow](https://github.com/idptools/sparrow)), and rebuilds it to those
dimensions. You can also ask for regions more compact or more expanded than predicted. And you
can generate several conformers into one multi-model PDB, which in VMD looks *like* a
simulation trajectory — to be very clear, it is **not** equivalent to a simulation, but it's
nice for visualization.

![DODO_EXAMPLE](https://github.com/idptools/dodo/blob/main/images/DODO_example.png)

> **DODO 2.0 is a rewrite and it breaks the 1.x API.** The `build.pdb_from_name()` /
> `pdb_from_pdb()` / `pdb_from_sequence()` functions and the three `pdb-from-*` console
> scripts are gone, replaced by `dodo.rebuild()` / `dodo.build_from_sequence()` and a single
> `dodo` command. See [Migrating from 1.x](#migrating-from-1x). The scientific behaviour also
> changed in ways worth reading about: build modes are now length-independent, and a
> multi-model run produces a genuine ensemble rather than one conformation repeated.

## Installation

Requires Python 3.10 or newer.

```bash
pip install git+https://github.com/idptools/dodo.git
```

The base install depends only on **numpy and scipy**, so it's fast and light. It gets you
structure reading and writing, region identification, and IDR rebuilding with the random-walk
engine. Without ALBATROSS it falls back to an analytical polymer scaling law for target
dimensions — see [Dimension prediction](#dimension-prediction) for how good that actually is.

For sequence-specific predictions, install the extra you need:

```bash
pip install "dodo[albatross]"    # ALBATROSS dimension prediction via sparrow
pip install "dodo[predictors]"   # metapredict, for sequence-only region identification
pip install "dodo[lookup]"       # resolve protein names to UniProt accessions
pip install "dodo[starling]"     # STARLING generative IDR ensembles (large: ~2.4 GB of weights)
pip install "dodo[viz]"          # matplotlib debug plotting
pip install "dodo[all]"          # everything above
```

`dodo[albatross]` is the one most users want. It's a separate extra because sparrow pulls in
torch, and a visualization tool shouldn't force that on someone who just wants to read a PDB.

## Command line

One command with subcommands.

```bash
# Rebuild a local structure
dodo rebuild AF-P04637-F1-model_v6.pdb -o p53_dodo.pdb

# Download from the AlphaFold database and rebuild
dodo fetch P04637 -o p53_dodo.pdb
dodo fetch "human p53" -o p53_dodo.pdb        # needs dodo[lookup]

# Build a disordered region from sequence alone
dodo sequence GRNQNGGGYQNYNNQGYQGHGGQHQNNYNQYPCNYFGPGYNN -o my_idr.pdb

# Just tell me what you think my structure looks like
dodo regions AF-P04637-F1-model_v6.pdb
```

### Shared flags

| Flag | Default | Meaning |
|---|---|---|
| `-o`, `--out` | *required* | Output PDB path |
| `-m`, `--mode` | `predicted` | Target dimension as a multiplier on the predicted end-to-end distance |
| `-n`, `--models` | `1` | Number of conformers. Folded domains are positioned once and held fixed across all models; only the disordered regions differ |
| `-e`, `--engine` | `walk` | `walk` or `starling` |
| `-s`, `--strategy` | `auto` | How to identify regions: `auto`, `density`, `contact`, `plddt`, `metapredict` |
| `--all-atom` | off | Place N, C and O for rebuilt regions ([caveat](#all-atom-output)) |
| `--sidechains` | off | Also place CB; only with `--all-atom` |
| `--seed` | none | Makes output reproducible |
| `--no-conect` | off | Omit CONECT records — **not recommended**, see [below](#why-conect-records-matter) |
| `-q`, `--quiet` | off | Suppress the per-region report |

Exit status is `0` on success, `2` if some regions could not be built, `1` on error.

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
report.ok            # True if every region in every model was rebuilt
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
print(assignment.describe())                     # chain A: IDR 1-31; FD 32-290; ...
print(assignment.score, assignment.threshold)    # the evidence behind every boundary

target = target_dimensions("GSGSGSGS...", mode="compact")
print(target)   # 62.6 A (compact, 0.7x of 89.4 A) via albatross over 144 residues
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
been positioned for it, and a terminal IDR is free at one end and can go almost anywhere. Build
the loose regions first and they occupy space the tight ones then cannot avoid.

**Step 3 is the one that makes the rest work, and it is easy to miss.** AlphaFold has no way to
know where to put two domains joined by a long disordered linker, so it packs them together —
measured on real models, **2–3.6× closer than the linker between them predicts** (p300's
151-residue linker: a 26.1 Å gap against a 94.9 Å prediction). Rebuilding the linker into that
gap can only produce a compact blob wedged between domains. The fix isn't a better linker
builder; it's moving the domains.

**Folded-domain atoms are never rebuilt.** They come from AlphaFold (or AlphaFold3, or a
crystal structure) and are trusted. A domain only ever moves as a rigid body, and DODO verifies
that: internal geometry is asserted unchanged to ~10⁻¹³ Å after every transform.

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

**This changed from 1.x, and it matters.** v1 expressed modes as Ångströms *per residue* —
`normal` was 0.8 Å/residue. That's linear in chain length, but real IDR end-to-end distance
scales as roughly $N^{0.52}$, so a fixed per-residue multiplier can only agree with the
prediction at one length. `normal` gave 80 Å at N=100 (about right) and 400 Å at N=500, where
the prediction is nearer 190 Å — and the error grew without bound.

Now a mode means the same thing at every length. One consequence: `normal` and `predicted` are
now synonyms, which they weren't in 1.x.

A target that exceeds what the chain can physically span is clamped to 95% of contour length
rather than being chased fruitlessly, and the clamping is reported.

## Dimension prediction

With `dodo[albatross]` installed, targets come from ALBATROSS. Without it, DODO falls back to
an analytical scaling law, $R_e = 6.22\,N^{0.522}$ (from Kohn *et al.* 2004), and **warns you**
that it did — a silent downgrade would make two runs of the same command disagree with no
visible cause.

How good is the fallback? Benchmarked against ALBATROSS over 72 sequences across six
compositional classes, its ratio to the prediction is 0.97 on average. Per class at N=500:

| | ratio | | ratio |
|---|---|---|---|
| polar | 0.95× | polyampholyte | 0.92× |
| proline-rich | 0.74× | polyanionic | 0.62× |
| polycationic | 0.64× | hydrophobic | 2.62× |

So for genuine IDR compositions it lands within roughly 0.6–0.95×, erring compact. The charged
classes are worst because polyelectrolyte expansion pushes their scaling exponent to 0.60–0.64.

The 2.62× outlier is worth understanding: a sequence drawn uniformly from all 20 amino acids is
hydrophobic-rich, and ALBATROSS correctly predicts a **collapsed globule** (measured scaling
exponent 0.20). Such a sequence isn't disordered at all, so no length-only law can describe it.
The error here is *compositional*, not a matter of tuning — true $R_e$ spans 3.4× across
compositions at fixed length — so refitting the exponent doesn't help, and was tried.

Use the fallback to keep a light install working, not to avoid installing sparrow for real work.

## Region identification

Four strategies behind one flag:

- **`density`** — DODO's original all-atom density metric: all-atom pairs within 8 Å per
  residue, thresholded at 480. **This is the method DODO was built and validated on, and it
  draws better boundaries than sequence-based disorder predictors.** Reimplemented over a
  KD-tree rather than changed — same numbers, 10.1 s down to 7 ms on a 1,086-residue model.
- **`contact`** — a CA-only alternative. Composition-free (every residue has exactly one CA)
  and invariant to whether side chains are modelled, which the density score is not. Useful for
  comparison and for CA-only input, but it is not the validated method.
- **`plddt`** — AlphaFold's own per-residue confidence, from the B-factor column. Cheap, but
  the density method beats disorder-based calls, so this is an explicit opt-in.
- **`metapredict`** — sequence-only. The backup, and the only option with no structure at all.
- **`auto`** (default) — `density` for all-atom input, `contact` for CA-only input (where a
  pair count can't be compared against the tuned threshold). It never picks pLDDT on its own,
  and it tells you what it chose.

Every assignment keeps its score profile and threshold, so a boundary you disagree with can be
audited rather than just re-run with different numbers.

You can also supply regions explicitly:

```python
from dodo.regions import assign_regions_from_spec

assign_regions_from_spec(structure, {"A": [("idr", 1, 40), ("folded", 41, 290)]})
```

Bounds here are **1-based inclusive**, matching what you read off a PDB file. Overlaps, gaps
and out-of-range bounds are rejected with an explanation — v1 accepted all of them silently and
failed later with something unrelated.

## Multi-model output is now a real ensemble

`-n 10` writes ten conformers as MODEL/ENDMDL frames. The folded domains are positioned once
and held fixed across every model, so a viewer flicking between frames sees only the disordered
regions move — if the domains were re-placed per model they would jump around and destroy the
illusion. One consequence worth knowing: a *connecting* IDR's end-to-end distance is fixed by
its anchors, so across models its path varies but its span cannot. A *terminal* IDR is free at
one end and does scatter.

**In 1.x they effectively didn't.** v1 placed folded domains once outside the model loop and
targeted only the *mean* predicted end-to-end distance, so every model shared one arrangement
and essentially one dimension. Measured on an early 2.0 engine, the coefficient of variation of
end-to-end distance across conformers was 0.006–0.045, where a matched physical reference gives
0.35–0.48. Sixty models of a 200-residue IDR spanning 1.9 Å of extension is one conformation
sampled sixty times.

A predicted end-to-end distance is the **mean of a distribution**, so each model now draws its
own target from the physical radial distribution $P(R) \propto R^2 e^{-3R^2/2\langle R^2\rangle}$.
Measured CV is now 0.36–0.53, and the ensemble mean still matches the prediction.

Two details worth knowing:

- A **single** model still hits the target exactly. One conformer should be the predicted
  dimension, not a random draw around it.
- Only regions with a **free end** scatter. An interior region pinned between two folded
  domains has its end-to-end distance determined by the anchor separation — there's nothing to
  sample, and its models differ in path rather than in span.

## Reproducibility

Everything stochastic takes a seed. Same seed, bit-identical output:

```bash
dodo rebuild model.pdb -o out.pdb --seed 42
```

v1 had no seed anywhere, so a stochastic structure builder wasn't reproducible and couldn't be
regression-tested at all.

## Why CONECT records matter

DODO writes CONECT records by default, and you should leave them on. CA–CA spacing is 3.81 Å,
past the automatic bond-detection cutoff in both VMD and PyMOL — so without CONECT a rebuilt
region renders as a cloud of disconnected dots. This isn't cosmetic polish; its absence defeats
the tool.

## Current limitations

Honestly stated, with what's fixed since 1.x marked.

1. **~~Rebuilding uses a simple random walk, so conformations aren't scientifically useful.~~**
   Partly addressed. The walk is now a self-avoiding, angle-constrained growth walk that hits
   the predicted dimensions and produces a genuine ensemble across models. It is still a
   geometric sampler, not a force field — for ensemble-grade conformations, use
   `--engine starling`.

2. **~~Rebuilt IDRs contain only alpha carbons.~~** Addressed, with a caveat — see
   [All-atom output](#all-atom-output) below.

3. **~~Unusual-bond warnings in VMD.~~** Addressed: correct CONECT records, correct atom-name
   columns, and the element column written. (v1 right-justified atom names from column 13 and
   omitted the element, with the result that MDTraj read its CA-only output as 912 *calcium*
   atoms.)

4. **~~Some visualization modes don't work in VMD; tube and trace fail.~~** Should be addressed
   by 2 and 3 together, but please report if you still see it.

5. **Assembly rebuilding is not implemented.** Multi-chain structures are read and written
   correctly, and regions are identified per chain, but rebuilding unmodelled regions of an EM
   assembly against the deposited sequence isn't wired up yet.

### All-atom output

`--all-atom` places N, C and O for rebuilt regions. It is **off by default** because it
currently refuses a fraction of generated traces, and that fraction grows with chain length:
roughly 12/12 accepted at 20 residues, 6/12 at 100, 2/12 at 200.

The cause is understood and measured. Reconstructability constrains the **change** in CA–CA–CA
pseudo-angle between consecutive residues, not its magnitude: the allowed change is about 35°
at a base angle of 95° and collapses to 4° at 150°. The generator doesn't yet enforce that joint
constraint, so a long chain almost always contains one unreconstructable junction. The fix is a
constraint on consecutive angles during candidate selection.

Side chains are CB only. A full rotamer library is deliberately not implemented rather than
fabricated.

## Development

```bash
git clone https://github.com/idptools/dodo.git
cd dodo
pip install -e ".[dev]"
pytest                      # 1400 tests
pytest -m "not slow"        # the fast subset
ruff check src tests && mypy
pre-commit install
```

Tests marked `network` hit the AlphaFold database, RCSB and UniProt; CI deselects them from the
main matrix and runs them separately so an upstream outage doesn't redden a pull request.

## Changes

### 2.0 — in development

A full rewrite. The 1.x API is gone; see [Migrating from 1.x](#migrating-from-1x).

**Correctness fixes that changed output**

- **AlphaFold downloads work again.** v1 hardcoded `AF-{id}-F1-model_v4.pdb`; EBI has since
  published `model_v6` and retired the older URLs, so v4 and v5 both 404 and *every*
  `pdb_from_name` call in 1.x fails today. The URL is now resolved from the AFDB API, which
  survives future version bumps.
- **Build modes are length-independent** (see [Build modes](#build-modes)).
- **Multi-model output is a real ensemble** (see [above](#multi-model-output-is-now-a-real-ensemble)).
- **Region identification: two domain-merging bugs fixed.** A single candidate folded block
  yielded *zero* folded domains, so a clean single-domain protein came out entirely disordered
  and its real domain was replaced by a random walk. And the gap before the last block was
  never tested, so an IDR between the last two blocks was absorbed into the final domain and
  never rebuilt.
- **The folded/disordered score no longer depends on composition.** v1 thresholded a raw
  atom-*pair* count, which scales with a residue's own heavy-atom count: measured within one
  folded domain, every glycine fell below the cutoff while 94% of Trp/Phe/Tyr sat above it. The
  score now counts neighbouring residues by their alpha carbons, so composition can't bias it —
  and it's invariant to whether side chains are present, which matters because DODO must handle
  full models, structures with unmodelled side chains, and its own output.
- **Reader data loss fixed.** Mid-chain selenomethionine no longer vanishes (which fabricated a
  phantom chain break); insertion codes keep residues 10 and 10A distinct; alternate conformers
  no longer duplicate atoms; multi-model files are no longer merged into one impossible
  structure; and hybrid-36 serials are decoded, so files over 99,999 atoms — i.e. the EM
  assemblies this tool targets — no longer crash outright.
- **CA–CA bond length is 3.81 Å everywhere.** v1 had 3.8, 3.856 and 3.89 live simultaneously.
- **Generated angles are restricted to what can be reconstructed to all-atom.** A CA
  pseudo-angle is coupled to N–CA–C for a trans peptide, so generating a wide angle can make the
  output un-reconstructable by construction.
- **Failure is never silent.** v1 builders returned coordinate arrays full of exact `(0,0,0)`
  rows on total failure, samplers returned NaN, and the domain placer returned positions it had
  already determined to be clashing. Everything now raises or reports an explicit success mask.

**New**

- Single `dodo` CLI with `rebuild` / `fetch` / `sequence` / `regions` subcommands.
- mmCIF reading, including `entity_poly`, `struct_ref` and unobserved-residue records.
- pLDDT-based region identification.
- Bond-length regularization for generative-model output (`dodo.geometry.regularize`). STARLING
  is a diffusion model, so its CA–CA distances scatter rather than sitting on the bond length.
  This is a constrained projection, not a rebuild: measured bond error 0.43–1.13 Å → 5×10⁻⁷ Å
  while preserving radius of gyration to within 0.03%.
- Seeds throughout, so output is reproducible.
- Type annotations, `mypy --strict` clean, `py.typed` shipped.

**Packaging and infrastructure**

- Python 3.10+; 3.8 and 3.9 are end-of-life.
- Base install needs only numpy and scipy. Everything heavy is an extra.
- **The distribution now contains the package.** 1.x shipped a wheel that excluded the entire
  backend — 12.7 KB with none of the engine — while `pip install .` still exited 0.
- **CI now runs.** The workflow's trigger key was `off:`, which YAML parses as boolean false, so
  GitHub never ran anything on a 3,800-line numerics package.

### 1.x

<details>
<summary>Earlier changelog</summary>

- **0.15** — 2024-07-09. Fixed building from local PDBs with metapredict when the structure was
  predicted fully disordered; fixed multiple models for fully disordered local PDBs.
- **0.14** — 2024-04-09. Fixed manual region assignment in `pdb_from_pdb`.
- **0.13** — 2023-12-11. Fixed fully-disordered proteins failing to build. Thanks to GitHub user
  alexpmagalhaes for a well-reported bug.
- **0.12** — 2023-11-15. Save predicted folded domains as individual PDBs.
- **0.11** — 2023-10-30. Small fixes, docs, typos.
- **0.10** — 2023-10-24. Multiple IDR models in one PDB for simulation-like visualization;
  command-line functionality.
- **0.06** — 2023-10-17. Approximate linear arrangement of folded domains.
- **0.05** — 2023-10-05. Backend overhaul; IDR-from-sequence.
- **0.04** — 2023-09-29. Keep folded-domain atoms.
- **0.03** — 2023-09-28. Save sequence-only IDRs.
- **0.02** — 2023-09-26. IDR coordinates from sequence.
- **0.01** — 2023-09-25. Initial release.

</details>

## Migrating from 1.x

| 1.x | 2.0 |
|---|---|
| `build.pdb_from_pdb(path, out_path=...)` | `dodo.rebuild(path)` → `dodo.write_pdb(report.models, out)` |
| `build.pdb_from_name(name, out_path=...)` | `dodo.fetch_alphafold(accession)` then `dodo.rebuild(...)`, or `dodo fetch` |
| `build.pdb_from_sequence(seq, out_path=...)` | `dodo.build_from_sequence(seq)` |
| `pdb-from-pdb file -o out` | `dodo rebuild file -o out` |
| `pdb-from-name "human p53" -o out` | `dodo fetch "human p53" -o out` |
| `pdb-from-sequence SEQ -o out` | `dodo sequence SEQ -o out` |
| `num_models=N` | `n_models=N` / `-n N` |
| `mode='normal'` (0.8 Å/residue) | `mode='predicted'` (1.0× predicted) — see [Build modes](#build-modes) |
| `use_metapredict=True` | `strategy='metapredict'` |
| `regions_dict={...}` | `assign_regions_from_spec(...)`, 1-based inclusive |
| `attempts_per_region`, `attempts_per_coord` | removed; retry budgets are internal and reported on failure |
| `graph=True` | removed; use `dodo.regions` and plot the models yourself |
| `beta_for_FD_IDR=True` | `write_pdb(..., annotate_regions=True)` |
| `just_fds=True` | not yet reimplemented |
| functions returned `None` and printed | functions return a `RebuildReport` |

## Copyright

Copyright (c) 2023-2026, Ryan Emenecker — Holehouse Lab

#### Acknowledgements

Originally based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms)
version 1.1.
