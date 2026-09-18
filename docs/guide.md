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
| `--backbone` / `--no-backbone` | **on** | Place N, C and O on the rebuilt regions; `--no-backbone` writes alpha carbons only — see below |
| `--ca-only` | off | Alpha carbons only, folded domains included |
| `-b`, `--annotate-regions` | off | Encode region type in the B-factor column, for colouring |
| `--no-conect` | off | Omit CONECT records — **not recommended**, see below |
| `-q`, `--quiet` | off | Suppress the per-region report and the progress bar |
| `--cache-structures` | off | `fetch` only: keep the downloaded model on disk — see [Caching](#caching) |

### Caching

DODO caches two things per user, and only the first is on by default. `--no-cache` (or
`DODO_NO_CACHE=1`) turns off both:

- **ALBATROSS predictions — on by default**, keyed by sequence hash and sparrow version. About
  116 bytes each. This one is worth more than it looks: the first prediction in a process imports
  sparrow, parrot and torch, which measures **1.57 s against 0.17 s of actual rebuilding** for a
  912-residue structure. A cache hit returns before sparrow is imported at all, so a repeat run
  never loads torch — measured end to end, 0.36 s against 1.96 s.

  The disk cost is genuinely negligible: rebuilding the entire AlphaFold human proteome (23,587
  structures) produced 9,172 entries totalling **1.07 MB**. It is nonetheless **capped at 10 MiB**,
  about ten times that, with the oldest entries dropped first — a backstop against a pathological
  workload rather than something normal use will meet.
- **Downloaded structures — off by default, opt in with `dodo fetch --cache-structures`.** Mean
  0.25 MB and largest 1.51 MB, measured over 259 AlphaFold models. Individually small, but they
  accumulate on disk unnoticed, so DODO does not keep them unless you say so: without the flag the
  model is downloaded to a temporary directory and deleted when the command exits. Either way the
  path it used is printed on every fetch, so you can see which happened.

Neither cache can change a result: predictions are deterministic per sequence, and the key includes
the sparrow version so an upgrade that changes the network invalidates them rather than serving
stale values.

### The progress indicator

One line on stderr, relabelled as the rebuild moves through its stages:

```text
reading 7R5J-assembly1.cif: 4.21M records [01:43]
identifying regions:  62%|██████    | 503/808 chains [00:55]
filling in unmodelled residues:  88%|████████▊ | 711/808 chains [01:24]
finding rigid units: 9.4k domain pairs [00:02]
positioning folded domains:  41%|████      | 760/1841 steps [00:28]
rebuilding (model 2/5):  17%|█▋        | 41.0k/241k res [03:12]
```

Every stage is covered, including reading the file, which on a large assembly is the longest
single wait — measured on the 597 MB, 4.9-million-atom nuclear pore, 125 s of the run happens
before region identification even starts. Every stage also *advances*, not just announces
itself: a bar only redraws when something moves, so a stage that sets a label and then works
silently leaves a frozen line, which is worse than no line at all.

Within the rebuild stage the count is weighted by residues rather than by region, because region
lengths span two orders of magnitude and a per-region bar sits still through the slowest part.

It is shown when stderr is a terminal, suppressed automatically when output is piped, and
suppressed by `-q`. From the Python API, pass `progress=True` or `progress=False` to override.
`dodo regions` and `dodo units` show it too — they are the commands you point at a large
structure first.

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

**The target is hit to within 10%, not exactly.** That allowance is deliberate rather than
sloppy: ALBATROSS's own predictions are good to a few percent, and the analytical fallback used
when it is unavailable is looser still, so insisting on tighter agreement would be precision
theatre. In practice the walk does much better — measured across the test corpus, a region steered
to a target lands within 0.27% of it at the median and 1.32% at the worst. The run summary names
any region that uses more than half the allowance, so you never have to take the fit on trust.

This applies to regions with a **free end**, which are the ones a target is aimed at. A region
pinned between two folded domains has no target to hit: its span is whatever the anchor positions
dictate, so DODO reports it as a span rather than scoring it against a number it was never steered
to. That is also why the domains are moved first — step 3 puts the anchors at the predicted
distance, and the linker is then built to fit them.

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

### Backbone reconstruction (`--backbone` / `--no-backbone`)

DODO places alpha carbons, and then — **by default** — infers **N, C and O** for the regions it
rebuilt, from those alpha carbons alone. Pass `--no-backbone` (or `backbone=False`) for
alpha-carbon-only output. Every bond it writes inside a rebuilt region is exact; the limits are at
the seams and in a handful of marginal steric contacts, both described below.

The idea is that a CA-only trace is far more informative than it looks. Four consecutive alpha
carbons define a pseudo-dihedral, and that angle largely determines where the intervening peptide
unit's N and C sit. DODO bins that angle at 20° and looks the offsets up in a table measured from
100 frames of all-atom IDR simulation, then closes each carbon onto the nitrogen already placed so
that every bond inside a peptide unit is exact by construction. A second pass rotates each unit
about its CA–CA axis — the single degree of freedom it has left — to improve bond angles, steric
clashes and φ/ψ plausibility together.

One angle is not the most a CA trace will tell you, though. Conditioning a unit on **both** the
pseudo-dihedrals that flank it — a five-carbon window rather than four — is measurably better, so
that is what ships: held out over 19,302 units, placed-atom error falls from 0.314 to 0.298 Å for C
(−5.1%) and 0.194 to 0.187 Å for N (−3.8%), each with a paired 95% confidence interval excluding
zero. Units at the end of a region have no fifth carbon and fall back to the four-carbon table.

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

The table is committed with its provenance: `scripts/derive_peptide_table.py` re-derives every
constant exactly from the frames in `tests/data/backbone/`, and `tests/unit/test_backbone_table.py`
runs that check in CI along with the placement accuracy — so the numbers above are reproducible
rather than a one-off measurement.

Oxygen is the loose one, and unavoidably: its position turns on ψ, which a CA trace constrains only
weakly. The N–CA–C bond angle comes out at 3.4° spread against 2.8° in real backbone, and all four
bond lengths are exact by construction.

**Folded domains are never affected.** They keep every atom they arrived with, `--backbone` or not.
Only the regions DODO generated gain atoms, which is what makes this additive rather than a
rewrite — and `--ca-only` remains the way to strip folded domains down.

**The seams are exact.** The old order of operations generated a CA-only region first and asked the
backbone pass to repair its folded-domain junction afterwards. That was too late: a peptide unit
can reach only 2.787 Å from a folded carbon to the next alpha carbon, or 2.854 Å from a generated
alpha carbon to the next folded nitrogen, while the frozen trace often left 3–5 Å. The visible gap
was therefore real geometry, not a viewer problem.

DODO now supplies the folded boundary C/N/O coordinates to the walk. The first and last generated
alpha carbons are accepted only when an exact peptide unit is reachable. At an N-terminal seam it
also rejects a CA if *every* exact nitrogen on the two-bond intersection circle would overlap the
folded carbonyl oxygen. Once the trace is fixed, backbone placement searches the full exact-bond
circle rather than testing only the point nearest its lookup-table prediction. Folded-domain atoms
still never move.

The paired A/B over dnmt3a, arf19 and p300 at seeds 0–2 changed **60 strained seams to zero**. Both
paths built every region and introduced zero impossible contacts or rebuilt-provenance bond
defects. Introduced steric-clash findings fell from 11 to 1 across the nine runs. Total measured
backbone runtime was 8.32 s before and 6.33 s after on the test machine; the constraint did not add
a slow minimization step and in difficult closures reduced wasted retries.

The conservative fallback still exists for inputs whose fixed folded geometry offers no valid
exact point: it leaves and reports a long bond on `RebuildReport.backbone_seams` instead of writing
an overlap. Normal committed fixtures do not exercise it.

:::{note}
Leaving the boundary residue un-rebuilt is still the wrong repair. Although its original alpha
carbon makes the peptide bond reachable, replacing that residue's nitrogen can drive it into its
existing side chain and can break a proline ring. Boundary-aware generation moves only generated
alpha carbons and preserves every folded-domain atom.
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
dodo sequence GRNQNGGGYQNYNNQGYQGHGGQHQNNYNQYPCNYFGPGYNN -o my_idr.pdb
```

### Why CONECT records matter, and what to do in VMD

This matters most under `--no-backbone`. CA–CA spacing is 3.81 Å, well past the distance cutoff
either VMD or PyMOL uses to infer bonds, so an alpha-carbon-only region has no bonds a viewer can
find for itself and renders as a cloud of disconnected dots without CONECT records. With the
backbone on — the default — the real N–CA and CA–C bonds are short enough for a viewer to find
unaided, but the records still declare the chain rather than leaving it to each viewer's
heuristics. They are written by default either way: leave them on.

DODO's CONECT output is complete, and `dodo validate` checks it: every bond it writes gets a record,
in the fixed columns the format specifies — the backbone bonds within and between rebuilt residues,
and, under `--no-backbone`, the consecutive CA–CA pairs instead (measured on p300, all 1,520 of
them). **PyMOL honours them and draws the trace.**

**VMD frequently does not.** VMD infers bonds by distance when it loads a PDB and its handling of
CONECT is unreliable, so a rebuilt region can still come up as dots however correct the file is. That
is not something DODO can fix from this side. The fix is to use a representation that does not
consult the bond list at all:

```
# In the VMD representation dialog, set Drawing Method to:
Trace        # connects consecutive alpha carbons directly
Tube         # same, drawn thicker
```

`Trace` and `Tube` follow residue order rather than bonds, which sidesteps the problem entirely and
is exactly right for an alpha-carbon-only model. If you need the bonds themselves — for a selection
or an analysis rather than a picture —
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
print(assignment.describe())                     # one line per chain, as `dodo regions` prints
print(assignment.score, assignment.threshold)    # the evidence behind every boundary

target = target_dimensions("GSGSGSGS" * 18, mode="compact")
print(target)   # 55.9 A (compact, 0.7x of 79.8 A) via albatross over 144 residues
```

### Reading the notes

`assign_regions` emits a note whenever it made a decision the region line alone cannot show. On the
command line these go to **stderr**, indented under the chain they belong to, so
`dodo regions x.pdb 2>/dev/null` gives clean region lines and `2>&1` keeps the reasoning. In Python
they are `assignment.notes`. There are four:

| Note | What happened | Consequence |
|---|---|---|
| `N short stretch(es) (… residues) scored as folded but came in under the 25-residue minimum for a folded domain` | A run scored folded but was too short to be one | It is treated as disordered, so it **is** rebuilt |
| `residues A-B (n residue(s)) scored as disordered but are shorter than the 4-residue minimum, so they were absorbed into the adjacent folded domain` | A disordered stretch too short to be worth a random walk | It keeps its input coordinates and is **not** rebuilt |
| `strategy auto-selected as contact -- the input has no side chains, so the density threshold does not apply` | `auto` declined `density` because the input has no side chains at all | The call is made on a metric DODO was not tuned on |
| `no folded domain found; the whole chain is treated as disordered` | Nothing cleared the threshold anywhere | The entire chain is rebuilt |

The two minimums are the tuned constants `min_folded_length` (25) and `min_idr_length` (4). Neither
is a command-line flag, but `assign_regions` takes both, alongside `threshold`, `min_loop_length`,
`max_internal_gap` and `min_seed_run`:

```python
import dodo

structure = dodo.read_structure("model.pdb")
print(dodo.assign_regions(structure)[0].describe())
# chain A: IDR 1-282; FD 283-432 (loops: 383-393); IDR 433-473; FD 474-912

structure = dodo.read_structure("model.pdb")
print(dodo.assign_regions(structure, min_folded_length=200)[0].describe())
# chain A: IDR 1-473; FD 474-912
```

DNMT3A's 150-residue folded domain does not clear a 200-residue minimum, so it is reclassified as
disordered and folded into the leading IDR — and will therefore be rebuilt. That is the mechanism
behind the first note in the table, seen from the other side. Moving a minimum steps outside the
defaults DODO was validated on, which is a different thing from fixing a boundary you disagree with;
for that, [override the region](#overriding-the-regions) directly.

The last two notes are worth treating as warnings rather than information. A whole chain called
disordered is a fair answer for a backbone-only fragment, but it is also what a misapplied strategy
looks like:

```bash
dodo regions tests/data/structures/idr_frame_backbone.pdb
# chain A: IDR 1-60
#   note: strategy auto-selected as contact -- the input has no side chains, so the density threshold does not apply
#   note: no folded domain found; the whole chain is treated as disordered
```

Asking for `-s density` on that same file prints the same region line and *only* the second note.
The auto-selection note is the one thing `auto` tells you that an explicit `-s` cannot, which is a
reason to leave the default alone.

### Residue numbering: three surfaces, two axes

DODO is 0-based positional internally and never uses PDB numbering as an index. Two of the three
places that number residues for a human convert to the **file's own numbering**; the third does not,
and the difference is invisible on the structures where they agree.

| Surface | Axis |
|---|---|
| `RegionAssignment.describe()`, and so `dodo regions` | the file's `residue_number`, with insertion code |
| `dodo regions --scores` labels, and the absorption note | the same |
| `assign_regions_from_spec` bounds | 1-based counted from the start of that chain |

They coincide exactly when a chain is numbered contiguously from 1, which is most AlphaFold models
and almost no crystal structures. `Structure.residue_id` is the converter in either direction:

```python
import dodo

structure = dodo.read_structure("model.pdb")
domain = dodo.assign_regions(structure)[0].domains[1]
print(structure.residue_id(domain.span.start))   # what describe() prints for this boundary
print(domain.span.start - structure.chains[0].span.start + 1)  # what a spec would take
```

Three cases make this more than pedantry, all of them real in `tests/data/structures/6kn7.cif`: a
chain numbered from `-1` (an expression tag), a chain whose numbering is not contiguous (138
residues spanning 87 to 272), and residues distinguished only by insertion code. A range is
therefore **endpoints, not a count** — use `len(domain)` for the number of residues — and on a
negative start the range reads as `-1-29`, which is the honest rendering of an unhelpful numbering.

### What an assignment carries

`assign_regions` returns one `RegionAssignment` per chain. On `model.pdb` (DNMT3A, 912 residues):

```python
import dodo

structure = dodo.read_structure("model.pdb")
assignment = dodo.assign_regions(structure)[0]

print(assignment.chain_id)          # 'A'
print(assignment.strategy.value)    # 'density' -- what auto actually resolved to
print(assignment.threshold)         # 480.0, the cutoff applied to score
print(assignment.n_folded, assignment.n_idrs)   # 2 2
print(assignment.fully_disordered)  # False
print(assignment.score.shape)       # (912,) -- per residue, chain-length
print(int(assignment.folded_mask.sum()))        # 485
print(assignment.notes)             # ()
```

Two of those mislead if taken at face value. `threshold` is only comparable within one strategy —
density thresholds at 480, contact at 10, pLDDT at 70, because they are three different scores.
And `folded_mask` is the **raw** thresholding result, before runs are merged and short ones
dropped: 485 residues here, against 589 actually covered by the final folded domains. For the
regions themselves, read `assignment.domains`, not the mask.

### Overriding the regions

DODO's metrics are a convenience, not a requirement. When you already know where the domains are,
state them and DODO will build exactly that, running no metric at all. Give each chain a list of
`(kind, start, stop)` regions and rebuild with `strategy="preset"` — *identify nothing, build what
is already there*:

```python
import dodo

structure = dodo.read_structure("model.pdb")
dodo.assign_regions_from_spec(structure, {
    "A": [
        ("idr",     1, 282),      # rebuilt from scratch
        ("folded", 283, 912),     # moved as a rigid body, never regenerated
    ],
})
report = dodo.rebuild(structure, strategy="preset", seed=0)
dodo.write_pdb(report.models, "my_regions.pdb")
```

#### The bounds

Inclusive residue numbers, **exactly as the input file numbers them**, which is the same numbering
`dodo regions` and `describe()` print — so a boundary copies straight from one into the other. See
[residue numbering](#residue-numbering-three-surfaces-two-axes) for why that is worth stating.
A residue distinguished only by an insertion code is named as a string:

```python
import dodo

structure = dodo.read_structure("model.pdb")
print(structure.residue_id(0), structure.residue_id(structure.n_residues - 1))
```

Bounds used to be a count from the start of the chain. On a chain numbered from 20, or on any chain
after the first, that selected the wrong residues silently — no error, and nothing in the output to
reveal it. If you have a script written against the old behaviour, a chain numbered contiguously
from 1 is unaffected; any other chain needs its bounds rewritten as real residue numbers.

#### Loops inside a folded domain

A folded domain's atoms are never regenerated, so a stretch inside it that you *do* want rebuilt has
to be declared as a loop. That is an optional fourth element:

```python
import dodo

structure = dodo.read_structure("model.pdb")
assignment = dodo.assign_regions_from_spec(structure, {
    "A": [("idr", 1, 282), ("folded", 283, 912, [(433, 473)])],
})[0]
print(assignment.describe())
# chain A: IDR 1-282; FD 283-912 (loops: 433-473)
```

That spec says *283-912 is one domain, and 433-473 is a loop within it* — as against DODO's own
reading of the same model, which is two folded domains with a linker between them. The distinction
matters to the builder: a linker's length sets how far apart two domains are placed, while a loop is
rebuilt in place between two fixed anchors and moves nothing.

A loop must lie **strictly inside** its domain, because it is rebuilt between an anchor residue on
each side. A stretch touching either boundary is a tail; declare it as its own `idr` region. Loops
may not overlap each other, and only a `folded` region may have them — an IDR is rebuilt in its
entirety, so a loop inside one is not a distinct region.

Earlier versions had no syntax for this, so a `describe()` → spec → rebuild round trip quietly
turned a loop into ordinary folded interior and never rebuilt it.

#### Several chains, and overriding only some of them

A spec takes one entry per chain, each in **that chain's own numbering** — two identical chains of a
homodimer both read `1` to `375`, not `1-375` and `376-750`. Chains you leave out keep whatever
regions they already have, loops included, so the usual pattern is to let the metric do the bulk and
replace only the chain you disagree with:

```python
import dodo

structure = dodo.read_structure("model.pdb")
dodo.assign_regions(structure)                    # every chain, by metric
dodo.assign_regions_from_spec(structure, {        # then replace just this one
    "A": [("idr", 1, 100), ("folded", 101, 912)],
})
report = dodo.rebuild(structure, strategy="preset", seed=0)
```

`strategy="preset"` then builds whatever is attached to each chain, however it got there — so a
run can mix DODO's own boundaries on most chains with your own on one.

#### What gets rejected

Every failure raises `InvalidRegionError` naming the residues you typed, before any building starts:

| Spec | Message |
|---|---|
| a gap between regions | `Chain 'A' has unassigned residues 283-299.` |
| regions that overlap | `Chain 'A' has overlapping domains: 1-300 and 250-912.` |
| not reaching a chain end | `Chain 'A' domains end at residue 900 but the chain ends at 912.` |
| a residue number the chain lacks | `Chain 'A' has no residue '9999' (start bound). …` |
| `stop` before `start` | `Region folded 432-283 of chain 'A' ends before it starts.` |
| a loop touching its domain edge | `Loop 283-393 must lie strictly inside its folded domain 283-432 …` |
| a loop on an `idr` | `Region idr 1-282 of chain 'A' declares loops, but only a folded domain can have them …` |
| a 2- or 5-element entry | `Region entry ('folded', 1) for chain 'A' has 2 elements. Use (kind, start, stop), or …` |

The pre-rewrite override path accepted the first four of these silently and failed much later with
something unrelated.

#### What a preset assignment does and does not carry

```python
import dodo
import numpy as np

structure = dodo.read_structure("model.pdb")
assignment = dodo.assign_regions_from_spec(structure, {
    "A": [("idr", 1, 282), ("folded", 283, 912)],
})[0]

print(assignment.strategy.value)                 # 'preset' -- no metric was run
print(np.isnan(assignment.threshold))            # True: nothing was thresholded
print(bool(np.isnan(assignment.score).all()))    # True: nothing was scored
print(int(assignment.folded_mask.sum()))         # 630 -- what you declared folded
```

`score` and `threshold` are `nan` rather than zero on purpose: a zero would read as a real
measurement. `folded_mask`, by contrast, is a fact about your own spec and is filled in.

If you want DODO's evidence for the boundaries you *didn't* dispute, start from its answer and edit
the one you did. A `Span` is frozen, so replace it rather than mutating it:

```python
import dodo
from dodo.structure import Span

structure = dodo.read_structure("model.pdb")
dodo.assign_regions(structure)                       # DODO's call, with its evidence
domain = structure.chains[0].domains[1]
print(domain.span)                                   # inspect
domain.span = Span(300, 400)                         # frozen: replace, don't mutate
report = dodo.rebuild(structure, strategy="preset")  # build your version
```

#### Python only

`strategy="preset"` builds region objects attached to the structure, and there is no command-line
syntax for expressing them, so `dodo rebuild -s preset` is rejected — the CLI offers `auto`,
`density`, `contact` and `plddt`. Manual specification is a Python API.

This replaces 1.x's `regions_dict=`, which took a separate stringly-typed description of the
structure alongside the real one, so the two could disagree.

## Complexes

### Rigid units: what may move, and what may not

Step 3 of DODO's algorithm moves folded domains so each linker IDR can reach its predicted
end-to-end distance. On one chain that is the point of the whole package. On a complex, applied
naively, it is a disaster: the domains of two chains that pack against each other were put there
by a prediction or an experiment, and pulling them apart to satisfy a linker prediction dismantles
the thing you handed DODO.

So folded domains are grouped into **rigid units**. A unit is a set of folded domains — possibly
from different chains — whose relative positions and orientations DODO will not change. A unit
moves as one rigid body or it does not move at all, and the invariant is exact: two atoms in the
same unit are the same distance apart in the output as in the input, to 1e-6 Å.

Domains rather than chains, deliberately. A chain in a complex routinely has some domains inside
the assembly and some dangling off a long linker; locking whole chains would forfeit step 3 for
every chain that touches another, while locking domains keeps it for exactly the domains that are
free to move.

Three things put two domains in one unit:

1. **Covalent continuity.** Two folded domains of one chain with no residues between them are
   bonded across the boundary. This is not a judgement call.
2. **A link that will not be rebuilt.** Residues between two folded domains that DODO will skip —
   anything under `min_length` — keep their input coordinates, so they are a rigid link and their
   neighbours are one unit.
3. **An interface.** Two domains in contact, by default only across chains. This is the only one
   of the three that is a judgement call; see the threshold discussion below.

### Looking before you build

```bash
dodo units complex.cif
```

```
3 rigid unit(s) over 6 folded domain(s); 1 hold more than one
  unit 0: A:193-222, A:473-912, B:193-222, B:473-912 [interface (344 residue contacts)]
  unit 1: A:281-433
  unit 2: B:281-433
  3/7 contact(s) locked
    A:193-222 .. A:473-912: 2 residue contact(s), 8 atom contact(s) -- NOT locked
    A:193-222 .. B:473-912: 344 residue contact(s), 5540 atom contact(s) -- locked
```

Every measured contact is listed, locked or not, because a contact DODO is prepared to break is
exactly what you need to see. Adjust with `--min-residue-pairs` and `--contact-radius`, or name the
units yourself in Python.

### Predicted or experimental

The choice is not really about units, it is about the input: is the arrangement of the folded
domains a guess DODO may improve on, or a measurement it must not touch? DODO's whole premise —
that a predicted inter-domain arrangement across a long linker is arbitrary — is a statement about
*predictions*. It is not true of a cryo-EM map.

| Setting | What it holds together |
|---|---|
| `--units predicted` | covalent continuity, plus contacts **across chains**. Domains joined by a linker are repositioned. |
| `--units experimental` | everything that arrived with coordinates. Nothing moves; only the disordered regions are rebuilt. |
| `--units auto` (default) | resolves to one of the two, from the file, and says which |
| `--units none` | covalent continuity only — the behaviour before rigid units existed |
| `--lock-intra-chain-interfaces` | under `predicted`, also hold two touching folded domains of the *same* chain |

`auto` decides on two facts about the file rather than on a guess about the science:

1. the file declares an experimental method — `EXPDTA` in a PDB, `_exptl.method` in an mmCIF.
   Measured across every fixture class present: AlphaFold DB models and AlphaFold 3 server output
   declare nothing at all, while 6kn7 and 7R5J declare `ELECTRON MICROSCOPY` in both formats.
   `THEORETICAL MODEL` is the PDB's own marker for a computed structure and does *not* count;
2. residues were filled in from a reference sequence. A prediction models every residue it was
   given, so having unmodelled residues is itself a property of an experiment.

Either one means experimental. The resolution is written into `report.notes` every time, because
a default that silently decides whether someone's measured domains may move is one they need to be
able to see and override.

```python
report = dodo.rebuild("model.pdb", seed=0, progress=False)
print([n for n in report.notes if "units=auto" in n][0])
```

`--units chains` is accepted as a literal-minded alias for `experimental`: holding every chain
rigid is exactly what "do not move anything that was in the input" comes to, because a linker only
ever joins two domains of the same chain.

The default interface threshold is three residue-residue contacts within 5 Å, and it is
deliberately low. The two errors are not symmetric. Locking a pair that did not need it costs only
that the pair stays where the input put it; *not* locking a pair that needed it takes a real
complex apart. Measured over five AlphaFold 3 complexes, 95 inter-chain folded-domain pairs are in
contact at all and the distribution is strongly bimodal — 7 pairs touch at exactly one residue pair
and 3 more at two, then the next value is 3 and the median is 37.

In Python, name the units outright when you know better than the heuristic:

```python
import dodo

report = dodo.rebuild(
    "complex.pdb",
    units=[[("A", 300), ("B", 300)]],   # any residue inside each folded domain to be held
    seed=0,
    progress=False,
)
for unit in report.rigid_units:
    print(unit)
```

### What the report tells you

`RebuildReport` gains four things on a complex:

```python
report.units                 # every rigid unit: what it holds, whether it moved, how far off
report.rigid_units           # the subset holding more than one folded domain
report.linkers               # every connecting IDR: predicted span against what it got
report.interfaces            # every measured contact, locked or not
report.unbridgeable_linkers  # flanking domains further apart than the residues can span
```

A linker whose two flanking domains are in one unit did not have its span set by DODO, and says so
rather than quietly reporting a prediction it never used. An unbridgeable linker counts against
`report.ok`: it is a property of the input and of the locking decision, knowable before any region
is attempted, and reporting it as a generic build failure would hide that.

### Placing a unit held by more than one linker

A unit reached by a single linker is placed the way DODO has always placed a domain: sample a
direction, put the attachment alpha carbon at the predicted separation along it, turn the unit so
its body extends away from the domain it connects back to, perturb, reject on clash.

A unit reached by two or more cannot be placed by choosing one direction, and in a complex that is
common rather than exotic — a heterodimer whose two chains each contribute a domain to each of two
units produces it immediately. Those are placed by projecting each attachment point onto the sphere
its own constraint defines, superposing onto the result, and iterating. The same step swept over a
whole component is what closes a cycle, and it is skipped entirely when the spanning-tree placement
already satisfies every constraint — which is every acyclic component, and therefore every
single-chain input.

## Structures with unresolved residues

A crystal or cryo-EM structure contains what could be resolved. What could not be resolved is,
overwhelmingly, disordered — which is to say it is exactly what DODO exists to rebuild, and it is
not in the file. Give DODO the full-length sequence and it fills it in.

```bash
dodo rebuild 7R5J.cif --fasta 7R5J_chain_sequences.fasta -o filled.cif
```

```python
import dodo

report = dodo.rebuild("model.pdb", sequences={"A": "MSKGEELFT..."}, seed=0, progress=False)
print(report.insertions.summary())
```

### Matching the FASTA to the chains

The format is trivial; the mapping is not. A FASTA from the RCSB names each chain twice — the
mmCIF `label_asym_id` and the author's `auth_asym_id` — and a biological assembly renames every
chain again for each symmetry copy:

```text
>7R5J_11|Chains TA[auth I0], UA[auth I1], VA[auth I2], WA[auth I3]|Nucleoporin p58/p45|Homo sapiens
```

DODO keys chains on the author id, so `I0` is the name that wins. Matching runs four passes,
cheapest and most certain first: the id verbatim; the id with a biological assembly's `-2`/`-3`
copy suffix stripped; the chain's observed sequence being a subsequence of exactly one record (or
of several records that all carry the *same* sequence, since identical entities are
interchangeable); and finally, nothing — which is reported rather than guessed at. A hand-written
`>A` or `>A some description` works too.

`--fill-missing` falls back to the sequence the file itself declares — `SEQRES` in a PDB,
`_entity_poly` in an mmCIF — for chains no FASTA record covers. It is off by default because it
changes what gets built for any file carrying a `SEQRES` record, which is most experimental
structures.

### Where the inserted residues go

Insertion happens **after** region identification, and that ordering is load-bearing. An inserted
residue has no coordinates; it gets a placeholder alpha carbon so the arrays stay well-formed, and
that placeholder is fiction. Burial is scored from coordinates, so letting region identification
see the placeholders would let fiction decide which residues are folded — a straight line drawn
through a protein core scores as buried, and the region you most wanted rebuilt would come out
classified as structure.

Regions are therefore decided on the observed geometry, and the inserted residues spliced into that
decision by where they fall:

- strictly inside a folded domain — an unresolved surface loop — becomes a **loop** of that domain,
  rebuilt between two fixed anchors. The domain stays one domain, which is what stops its two
  halves being pulled apart;
- inside or against a disordered region extends that region;
- between two folded domains with nothing between them becomes a new connecting IDR;
- past either end of the chain extends, or creates, a terminal IDR.

Two folded domains whose connecting linker was **inserted** rather than observed are held rigid.
A linker nobody ever saw says nothing about where those domains sit relative to each other, but
their own coordinates do, and moving them apart to satisfy a prediction for residues nobody
observed would replace a measurement with a guess.

### What is refused rather than guessed

- A reference that disagrees with more than 10% of a chain's modelled residues is a different
  protein; the chain is left alone and the report says so.
- A gap whose flanking residues are further apart than the missing residues can physically span
  has no solution. It is reported before anything is inserted, not attempted.
- A region that fails to build has its placeholder residues **removed from the output**. Every
  other region DODO fails to build keeps its input coordinates, which is honest because they came
  from the input; here there is nothing to keep, so writing the placeholder would mean writing a
  straight line and calling it a structure.

Residue numbering is preserved: an inserted residue takes the number the file's own numbering gap
left free. When the author numbering has no room — which means the numbering disagrees with the
sequence — the whole chain is renumbered from the reference and the report says so, rather than
producing a file whose numbering quietly disagrees with itself. Inserted residues are marked on the
model, so you can always tell which is which:

```python
model = report.models[0]
model.inserted            # (n_residues,) bool: True where DODO added the residue
```

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
