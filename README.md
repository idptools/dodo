DODO: re<ins>D</ins>esign c<ins>O</ins>mputationally generate<ins>D pr<ins>O</ins>teins
==============================

## What is DODO?

Protien structure predictors have revolutionized biology. However, many predicted 
structures have low-confidence regions that do not adopt any predicted secondary 
structure. DODO takes these regions and rebuilds them under the assumption that they
are disordered regions. At this time, the rebuilt regions include the backbone of the 
amino acids but do not include the side chains. This is an improvement over the first
version of DODO because we now have the backbone instead of just alpha carbons. I'm 
hoping to implement all atom IDR rebuilds in the future. 

To be clear: the work that the various groups do to make structure predictors is **amazing**, 
and none of this takes away from it in *any way*. But for figures and presentations it's
useful for disordered regions *look* like what they are.

By default, DODO identifies each disordered region and predicts its end-to-end distance from sequence
— with ALBATROSS ([sparrow](https://github.com/idptools/sparrow)).
You can also ask for regions more compact or more expanded than predicted. And you
can generate several conformers into one multi-model PDB, which in VMD looks *like* a
simulation trajectory — to be very clear, it is **not** equivalent to a simulation, but it's
nice for visualization.


![DODO_EXAMPLE](https://github.com/idptools/dodo/blob/main/images/DODO_ensemble.gif)

> **DODO 2.0 is a rewrite and it breaks the 1.x API.** The `build.pdb_from_name()` /
> `pdb_from_pdb()` / `pdb_from_sequence()` functions and the three `pdb-from-*` console
> scripts are gone, replaced by `dodo.rebuild()` / `dodo.build_from_sequence()` and a single
> `dodo` command. See [What changed, and migrating from 1.x](#what-changed-and-migrating-from-1x).

## Installation

Requires Python 3.10 or newer.

```bash
pip install git+https://github.com/idptools/dodo.git
```

The base install depends on **numpy, scipy, numba, tqdm, getSequence, and SPARROW**.


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
| `--cache-structures` (`fetch` only) | off | Keep the downloaded AlphaFold model on disk; otherwise it is deleted when the command exits |

### Caching

**Downloaded structures are not kept.** `dodo fetch` downloads to a temporary directory and
deletes it on exit. Pass `--cache-structures` if you want repeat fetches of the same accession to
skip the download. A model averages 0.25 MB and runs to 1.51 MB — small individually, but a
directory you never chose to grow, so keeping them is your decision rather than DODO's.

**ALBATROSS predictions are cached, and this one is on by default.** The reason is speed, not
convenience: a cache miss has to import sparrow, which pulls in parrot and torch. Measured on a
912-residue structure, that import is **1.6 s against 0.17 s of actual rebuilding** — it dominates
the run. Every CLI invocation is a fresh process, so without the cache you pay it every single
time:

| | full `dodo rebuild` |
|---|---|
| prediction cache warm | **0.36 s** |
| prediction cache disabled | **1.96 s** |

The disk cost is negligible, and that is measured rather than assumed. An entry is a sequence hash
and a float — **116 bytes**. Rebuilding *the entire AlphaFold human proteome*, all 23,587
structures and every disordered region in them, produced **9,172 entries totalling 1.07 MB**. The
whole human proteome costs about one megabyte.

The cache lives at `~/.cache/dodo/<generation>/` (or the platform equivalent) and is keyed by
sequence hash *and* sparrow version, so upgrading sparrow invalidates old values rather than
silently serving predictions from a different network.

**It is capped at 10 MiB**, roughly ten times the whole human proteome, with the oldest entries
dropped first once it is reached. Nothing you do in normal use will approach that; the cap is
there so no workload, however unusual, can grow the file without bound.

**To turn it off:** `--no-cache` on any command, or `DODO_NO_CACHE=1` in the environment. Both
disable prediction caching and structure caching together. DODO still works exactly the same — it
is only slower.

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

DODO writes CONECT records by default, and you should leave them on. They matter most under
`--no-backbone`: CA–CA spacing is 3.81 Å, past the automatic bond-detection cutoff in both VMD and
PyMOL, so without CONECT an alpha-carbon-only region renders as a cloud of disconnected dots. With
the backbone on (the default) a viewer can find the real N–CA and CA–C bonds for itself, but the
records still declare the chain explicitly rather than leaving it to each viewer's heuristics.

## Current limitations

Honestly stated, with what's fixed since 1.x marked.

1. **~~Rebuilding uses a simple random walk, so conformations aren't scientifically useful.~~**
   Partly addressed. The walk is now a self-avoiding, angle-constrained growth walk that hits
   the predicted dimensions and produces a spread of distinct conformers across models — but it is
   a geometric sampler, not a force field, so that spread is not a thermodynamic ensemble and is
   not a substitute for a simulation.

2. **~~Rebuilt IDRs contain only alpha carbons.~~** Addressed: rebuilt regions now get a full
   N, CA, C, O backbone by default, and `--no-backbone` returns the alpha-carbon-only output if you
   want it. Side chains are still not built. Note the alpha-carbon-only limitation never applied to
   the regions DODO leaves alone, which keep every atom; see
   [Atoms in the output](#atoms-in-the-output) below.

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

7. **Cis-peptide bonds are not modelled, so a cis-proline inside a rebuilt region comes out
   trans.** DODO builds every virtual CA–CA bond at 3.81 Å, the trans value. A cis peptide — in
   practice almost always X–Pro — sits near 2.9 Å, and DODO cannot produce one.

   This only affects regions DODO *rebuilds*. Folded domains are moved as rigid bodies and never
   regenerated, so a cis-proline in one survives untouched, to the last decimal place.

   Measured on 1,200 AlphaFold human structures, 17,406 CA–CA bonds are short (< 3.30 Å) and 4,482
   of those are X–Pro. Most are not real cis-prolines: short bonds have a mean pLDDT of 38.8
   against 72.2 elsewhere, so the bulk of that population is AlphaFold producing compressed
   geometry where it is unsure. Filtering to confident ones (pLDDT ≥ 70) leaves 831, and **32 of
   them — 3.9% — fall in a region DODO rebuilds.** That is about one affected bond per 38
   structures, or roughly 600 across the entire human proteome. Where the short bond is *not*
   confident, rebuilding it to clean geometry is an improvement rather than a loss.

   It is not planned. Supporting cis would make the bond length per-bond rather than a constant,
   and that constant is load-bearing in 65 places across 12 modules — reachability schedules,
   closure geometry, the cone sampler, the clash exclusions and the peptide-plane table, which was
   measured on trans units only and would need a cis counterpart. The deeper problem is that it is
   not well-posed for *de novo* generation: DODO invents a new conformation and has no way to know
   which prolines should be cis, and in a real IDR the two states interconvert rather than one
   being correct. If you need cis-prolines preserved in a specific region, exclude that region from
   rebuilding rather than expecting DODO to reproduce it.

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
detail exactly as AlphaFold produced it, down to the last decimal place. The regions DODO actually
rebuilds get a backbone — N, CA, C and O per residue — and no side chain.

On dnmt3a, for example: 4,636 atoms across the 578 residues DODO leaves alone, preserved
bit-for-bit, plus 1,336 backbone atoms for the 334 residues it rebuilds — 5,972 in total. Under
`--no-backbone` those 334 residues are one alpha carbon each instead, for 4,970.

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

- **A multi-model run is now a real ensemble.** 1.x placed folded domains once outside the model
  loop, so every model shared one arrangement *and* essentially one end-to-end distance. Each
  model now draws its own dimensions; folded domains are still positioned once and held fixed, so
  only the disordered regions move between frames.
- **`--seed` exists.** 1.x had no seed anywhere, so a stochastic builder was not reproducible.
- **Rebuilt regions get a backbone by default.** N, C and O are placed on them; `--no-backbone`
  returns the alpha-carbon-only output 1.x produced. Side chains are still not built.
- **Exit status means something:** `0` success, `2` a region of 10+ residues could not be rebuilt,
  `1` error.
- **Installation is lighter.** metapredict, matplotlib and cython are gone from the dependencies;
  sparrow is now an explicit separate install rather than a git dependency resolved for you.

## Copyright

Copyright (c) 2023-2026, Ryan Emenecker — Holehouse Lab

#### Acknowledgements

Originally based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms)
version 1.1.
