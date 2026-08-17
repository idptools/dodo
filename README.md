DODO: re<ins>D</ins>esign c<ins>O</ins>mputationally generate<ins>D pr<ins>O</ins>teins
==============================

## What is DODO?

Protien structure predictors have revolutionized biology. However, many predicted 
structures have low-confidence regions that do not adopt any predicted secondary 
structure (typically IDRs). DODO takes these regions and rebuilds them under the assumption that they
are disordered regions.

To be clear: the work that the various groups do to make structure predictors is **amazing**, 
and none of this takes away from it in *any way*. But for figures and presentations it's
useful for disordered regions *look* like what they are.

By default, DODO identifies each disordered region and predicts its end-to-end distance from sequence
— with ALBATROSS ([sparrow](https://github.com/idptools/sparrow)) when sparrow is installed, and an
analytical scaling law otherwise — then repositions the folded domains as needed and rebuilds the IDR
to the predicted dimensions. 
You can also ask for regions more compact or more expanded than predicted. And you
can generate several conformers into one multi-model PDB, which in VMD looks *like* a
simulation trajectory — to be very clear, it is **not** equivalent to a simulation, but it's
nice for visualization.

![DODO_EXAMPLE](https://github.com/idptools/dodo/blob/main/images/DODO_example.png)

> **DODO 2.0 is a rewrite and it breaks the 1.x API.** The `build.pdb_from_name()` /
> `pdb_from_pdb()` / `pdb_from_sequence()` functions and the three `pdb-from-*` console
> scripts are gone, replaced by `dodo.rebuild()` / `dodo.build_from_sequence()` and a single
> `dodo` command. See [What changed, and migrating from 1.x](#what-changed-and-migrating-from-1x). The scientific behaviour also
> changed in ways worth reading about: build modes are now length-independent, and a
> multi-model run produces a spread of distinct conformers rather than one conformation repeated.

## Installation

Requires Python 3.10 or newer.

```bash
pip install git+https://github.com/idptools/dodo.git
```

It is not yet on PyPI, so the git URL is the install for now. When the wheel is published the command
becomes `pip install idptools-dodo` — one step, no extra sources.

The base install depends on **numpy, scipy, numba, tqdm and getSequence** — small, fast, and no torch. It
gets you structure reading and writing, region identification, IDR rebuilding, protein-name lookup
for `dodo fetch`, and the validator.

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
| `-s`, `--strategy` | `auto` | How to identify regions: `auto`, `density`, `contact`, `plddt` |
| `--seed` | none | Makes output reproducible |
| `--backbone` / `--no-backbone` | **on** | Place N, C and O on the rebuilt regions, inferred from the alpha carbons; `--no-backbone` writes alpha carbons only |
| `--ca-only` | off | Alpha carbons only, folded domains included |
| `-b`, `--annotate-regions` | off | Encode region type in the B-factor column, for colouring |
| `--no-conect` | off | Omit CONECT records — **not recommended**, see [below](#why-conect-records-matter) |
| `-q`, `--quiet` | off | Suppress the per-region report and the progress bar |

DODO caches ALBATROSS predictions (~116 bytes each; a hit skips importing torch, worth 1.5 s) and
downloaded structures (mean 0.25 MB, largest 1.51 MB over 259 measured files). Both are on by
default; `--no-cache` or `DODO_NO_CACHE=1` opts out of both.

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
crystal structure) and are trusted. A domain only ever moves as a rigid body, and DODO checks
that: after each domain is repositioned, `verify_rigid` asserts its internal geometry is unchanged
to within **10⁻⁶ Å** — a bound set to catch a non-rigid transform, not floating-point noise. The
residual of a float64 rigid move comes in far under it: measured across these structures, ~10⁻¹³ Å.

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
  residue, thresholded at 480. **This is the method DODO was built and validated on**; the author
  reports it draws better region boundaries than sequence-based disorder predictors, though that
  comparison predates this repository and is not benchmarked here. Reimplemented over a
  KD-tree rather than changed — same numbers, 10.1 s down to 7 ms on a 1,086-residue model.
- **`contact`** — a CA-only alternative. Composition-free (every residue has exactly one CA)
  and invariant to whether side chains are modelled, which the density score is not. Useful for
  comparison and for CA-only input, but it is not the validated method.
- **`plddt`** — AlphaFold's own per-residue confidence, from the B-factor column. Cheap, but
  density is DODO's validated default, so this is an explicit opt-in.
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
   the predicted dimensions and produces a spread of distinct conformers across models — but it is
   a geometric sampler, not a force field, so that spread is not a thermodynamic ensemble and is
   not a substitute for a simulation.

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

### Backbone reconstruction (`--backbone` / `--no-backbone`)

By default DODO places **N, C and O** on the rebuilt regions, inferred from the alpha carbons alone
— consecutive alpha carbons largely determine where the peptide unit between them sits, and DODO
looks that up in a table measured from 100 frames of all-atom IDR simulation, then settles each
unit's one remaining degree of freedom against bond angles, clashes and φ/ψ together. The lookup is
keyed on **five** alpha carbons (both pseudo-dihedrals flanking the unit), which is measurably
better than the four-carbon form it replaced: held out, C improves 5.1% and N 3.8%, each with a
paired 95% confidence interval excluding zero. Pass `--no-backbone` (or `backbone=False` in the
API) for alpha-carbon-only output.

Held out properly — table rebuilt from 80 frames, scored on the 20 it had never seen, 3,643
residues:

| Atom | Mean error | Median |
|---|---|---|
| N | 0.16 Å | 0.12 Å |
| C | 0.22 Å | 0.14 Å |
| O | 0.63 Å | 0.40 Å |

The table is reproducible, not a set of magic numbers: `scripts/derive_peptide_table.py` re-derives
it exactly from the committed frames in `tests/data/backbone/`, and `tests/unit/test_backbone_table.py`
pins that regeneration and the placement accuracy in CI.

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

## Development

```bash
git clone https://github.com/idptools/dodo.git
cd dodo
pip install -e ".[dev]"
pytest                      # ~2000 tests
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
