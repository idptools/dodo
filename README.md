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
> `dodo` command. See [What changed, and migrating from 1.x](#what-changed-and-migrating-from-1x). The scientific behaviour also
> changed in ways worth reading about: build modes are now length-independent, and a
> multi-model run produces a genuine ensemble rather than one conformation repeated.

## Installation

Requires Python 3.10 or newer.

```bash
pip install git+https://github.com/idptools/dodo.git
```

Two installs, and that is the whole story:

```bash
pip install idptools-dodo                 # everything you need to rebuild a structure
pip install "idptools-dodo[starling]"     # adds STARLING ensembles (large: ~2.4 GB of weights)
```

The base install depends on **numpy, scipy and getSequence** — small, fast, and no torch. It gets
you structure reading and writing, region identification, IDR rebuilding, protein-name lookup for
`dodo fetch`, and the validator.

STARLING is the only extra, because it is the only dependency heavy enough to be worth opting into.

**Most users should also install sparrow**, which provides ALBATROSS:

```bash
pip install git+https://github.com/idptools/sparrow.git
```

Without it DODO falls back to an analytical polymer scaling law that is blind to sequence
composition — for a 100-residue poly-glutamate region ALBATROSS predicts 122.2 Å and the fallback
estimates 68.8 Å. DODO warns when it has fallen back, but install sparrow. It is a separate step
rather than an extra because sparrow is not published on PyPI, and PyPI does not permit an extra to
reference a git URL.

The distribution is named `idptools-dodo` because PyPI's `dodo` belongs to an unrelated 2014
project. The import name is unaffected — it is always `import dodo`.

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

| Flag | Default | Meaning |
|---|---|---|
| `-o`, `--out` | *required* | Output PDB path |
| `-m`, `--mode` | `predicted` | Target dimension as a multiplier on the predicted end-to-end distance |
| `-n`, `--models` | `1` | Number of conformers. Folded domains are positioned once and held fixed across all models; only the disordered regions differ |
| `-e`, `--engine` | `walk` | `walk` or `starling` |
| `-s`, `--strategy` | `auto` | How to identify regions: `auto`, `density`, `contact`, `plddt` |
| `--seed` | none | Makes output reproducible |
| `--no-conect` | off | Omit CONECT records — **not recommended**, see [below](#why-conect-records-matter) |
| `-q`, `--quiet` | off | Suppress the per-region report |

Exit status is `0` on success, `2` if a region of 10 residues or more could not be rebuilt, `1` on
error. A shorter region that could not be rebuilt is reported and left as it arrived, and does not
fail the run — a few residues of AlphaFold geometry do not spoil a figure, whereas a long extended
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
print(assignment.describe())                     # chain A: IDR 1-31; FD 32-290; ...
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

With sparrow installed, targets come from ALBATROSS. Without it, DODO falls back to
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
- **`auto`** (default) — `density` for all-atom input, `contact` for CA-only input (where a
  pair count can't be compared against the tuned threshold). It never picks pLDDT on its own,
  and it tells you what it chose.

Every assignment keeps its score profile and threshold, so a boundary you disagree with can be
audited rather than just re-run with different numbers.

### Overriding the regions yourself

If you disagree with a boundary, say so and DODO will build exactly what you asked for:

```python
import dodo

structure = dodo.read_structure("model.pdb")
dodo.assign_regions_from_spec(structure, {"A": [("idr", 1, 40), ("folded", 41, 912)]})
report = dodo.rebuild(structure, strategy="preset")
```

`strategy="preset"` means *identify nothing, build what is already there*. Bounds are **1-based
inclusive**, matching what you read off a PDB file, and the regions must **tile the whole chain** —
every residue belongs to exactly one. Overlaps, gaps and out-of-range bounds are rejected with an
explanation naming the offending residue, where v1 accepted all of them silently and failed later
with something unrelated.

This replaces v1's `regions_dict=` parameter. That took a separate, stringly-typed description of
the structure alongside the real one, so the two could disagree; here there is one representation,
it is already validated, and it carries the score profile that produced it. You can equally start
from DODO's own answer and adjust it:

```python
import dodo

structure = dodo.read_structure("model.pdb")
dodo.assign_regions(structure)                       # DODO's call, with its evidence
structure.chains[0].domains[1].span                  # inspect, then edit as you like
report = dodo.rebuild(structure, strategy="preset")  # build your version
```

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

2. **Rebuilt IDRs contain only alpha carbons.** Still true in 2.0, and deliberately so — getting
   the alpha-carbon trace right comes first. Note this has never applied to the regions DODO
   leaves alone, which keep every atom; see [Atoms in the output](#atoms-in-the-output) below.

3. **~~Unusual-bond warnings in VMD.~~** Addressed: correct CONECT records, correct atom-name
   columns, and the element column written. (v1 right-justified atom names from column 13 and
   omitted the element, with the result that MDTraj read its CA-only output as 912 *calcium*
   atoms.)

4. **~~Some visualization modes don't work in VMD; tube and trace fail.~~** Should be addressed
   by 2 and 3 together, but please report if you still see it.

5. **Assembly rebuilding is not implemented.** Multi-chain structures are read and written
   correctly, and regions are identified per chain, but rebuilding unmodelled regions of an EM
   assembly against the deposited sequence isn't wired up yet.

6. **A single marginal steric contact survives, across the whole test corpus.** Measured over
   117 structures: **one** contact, at 3.02 Å against a 3.20 Å limit. None is below the 1.00 Å
   impossible floor, and none appears in an unmodified input. See
   [Anchor exemptions](#anchor-exemptions) for why one remains and why that is the right place to
   stop.

7. **Three regions out of 117 structures are left unbuilt, and two are not DODO's doing.**
   Measured:

   - One input has two *fixed* residues 3.04 Å apart — already closer than the clash distance
     DODO is required to satisfy — so it is asked to thread a 16-residue loop between anchors
     that clash with each other.
   - One input contains a chain break: consecutive alpha carbons 5.26 Å apart where a real chain
     is 3.81 Å. DODO needs that residue to constrain a junction angle and cannot.
   - One is a genuine build failure: a 7-residue tail at the very end of a 189-residue chain,
     where the walk exhausted its attempts.

   In every case the region keeps its input coordinates and the reason is reported. DODO never
   substitutes degenerate output for a region it could not build.

### Anchor exemptions

Rebuilding a region means attaching it to the fixed residues on either side, and that requires
exempting those *anchors* from clash checking to some degree. This is unavoidable, so it is worth
being plain about it.

The anchors' **alpha carbons** are always exempt. The region is bonded to them; treating them as
obstacles would make every valid attachment register as a clash and there would be nothing to
build. There is no version of the algorithm without this.

The anchors' **backbone** — N, C, O — is a judgement call, and DODO makes it per region rather than
globally. A residue bonded to an anchor genuinely does come closer to its backbone than the clash
distance: measured over 649,658 sequence-neighbour pairs from the human proteome, an alpha carbon
sits 2.379 Å from the next residue's N at the 0.1st percentile and 3.280 Å at the median, so the
*median* real junction is already inside the 3.20 Å limit. But that exemption is not per-residue,
so granting it also licenses a residue 3 to 16 positions further along, which has no such claim.

So DODO tries twice. First with only the alpha carbons exempt, which is the honest constraint. Only
if the region cannot be built that way does it retry with the backbone exempt as well — and when it
does, it says so, both on the region's outcome and in the run summary.

The measured effect of doing it this way, over the 117-structure corpus:

| | contacts DODO introduces | regions left unbuilt |
|---|---|---|
| backbone always exempt | 21 | 4 |
| backbone never exempt | 1 | 8 |
| **two passes, strict first** | **1** | **3** |

Both halves improve, which is why this is the default. A region left unbuilt stays as AlphaFold
spaghetti and is glaringly visible in a figure; a contact a fraction of an Ångström inside the limit
is not visible at all. Given DODO is a visualization tool, that is the trade worth making — and the
remaining contact is one, not a class.

### Atoms in the output

Regions DODO does **not** rebuild keep every atom they arrived with. Folded domains are moved as
rigid bodies — translated and rotated, never regenerated — so each one retains its full atomic
detail exactly as AlphaFold produced it, down to the last decimal place. Only the regions DODO
actually rebuilds are reduced to one alpha carbon per residue.

On dnmt3a, for example: 4,636 atoms across the 578 residues DODO leaves alone, preserved
bit-for-bit, and 334 alpha carbons for the 334 residues it rebuilds.

Building backbone (N, C, O) or side-chain atoms **for rebuilt regions** is not part of 2.0. The
machinery exists behind the `all_atom=` and `sidechains=` keyword arguments of `dodo.rebuild`, but
it is not exposed on the command line, because it currently refuses a fraction of generated traces
and that fraction grows with chain length: roughly 12/12 accepted at 20 residues, 6/12 at 100,
2/12 at 200.

The cause is understood and measured. Reconstructability constrains the **change** in CA–CA–CA
pseudo-angle between consecutive residues, not its magnitude: the allowed change is about 35°
at a base angle of 95° and collapses to 4° at 150°. The generator doesn't yet enforce that joint
constraint, so a long chain almost always contains one unreconstructable junction. The fix is a
constraint on consecutive angles during candidate selection, and it is the next feature after 2.0.

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

## What changed, and migrating from 1.x

**2.0 breaks the 1.x API deliberately.** The full list — every behaviour change, the complete
old-to-new translation table, and the known limitations with their measurements — lives in
[CHANGELOG.md](CHANGELOG.md), which is the single copy so the two cannot drift.

The three that will bite hardest:

1. **Build modes mean something different.** 1.x's were Ångströms per residue; 2.0's are
   multipliers on the predicted end-to-end distance. Short regions come out **larger** than 1.x
   made them, long ones smaller, crossing over near 80 residues. Regenerating a 1.x figure will
   not reproduce it, and 2.0's version is the correct one.
2. **The entry points are gone.** `build.pdb_from_name()` / `pdb_from_pdb()` /
   `pdb_from_sequence()` and the three `pdb-from-*` commands are replaced by `dodo.rebuild()` /
   `dodo.build_from_sequence()` and a single `dodo` command.
3. **Functions return a report** instead of returning `None` and printing, so a script can finally
   tell whether a region was actually rebuilt.

## Copyright

Copyright (c) 2023-2026, Ryan Emenecker — Holehouse Lab

#### Acknowledgements

Originally based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms)
version 1.1.
