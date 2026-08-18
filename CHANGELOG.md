# Changelog

This changelog starts at 2.0. Earlier releases predate it.

It records what changes **for you as a user** — the commands and functions you call, the
structures that come out, and the things you have to do differently. Internal restructuring is
not listed, however large, unless you can see it from outside.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] — 2.0.0

DODO 2.0 is a rewrite. **It breaks the 1.x API deliberately**, so read the migration notes below
before upgrading a script.

The scientific behaviour changed too, in ways you can see in your output. The two that matter
most: build modes now mean something different, and a multi-model run is now a real ensemble.

### Installing is simpler: one package, no extras

```bash
pip install idptools-dodo                 # everything you need
```

`import dodo` is unchanged and always will be. Only the distribution name moved, because PyPI's
`dodo` belongs to an unrelated project from 2014.

**metapredict is gone**, and with it the `[predictors]` extra and the `metapredict` region-
identification strategy. Its only purpose in 1.x was faster region identification, and that reason
no longer exists: the all-atom density metric it was a workaround for now runs in **7 ms** on a
1,086-residue model, down from 10.1 s. metapredict requires torch, pytorch-lightning, cython and
matplotlib, so dropping it is most of what makes the base install light.

**`[lookup]` is gone** — getSequence is now installed by default, so `dodo fetch "human p53"` works
out of the box. It is three small packages (urllib3, requests, protfasta) with no torch in the
chain; putting it behind an extra bought nothing and cost a confusing failure for anyone who had
not read the README.

**`[viz]` is gone.** It installed matplotlib for a debug plotter that did not exist.

**`[all]` is gone.** With no extras left there is nothing to union.

There is no `[albatross]` extra either, and there never should be. It pointed at a PyPI package
called `sparrow` that is **not** the sparrow DODO needs — PyPI's `sparrow` is an RDF/SPARQL library,
and installing that extra quietly put the wrong package on your import path. sparrow is not
published on PyPI and an extra may not reference a git URL, so it is an explicit step:

```bash
pip install git+https://github.com/idptools/sparrow.git
```

Install it. Without it DODO falls back to a sequence-blind approximation, and says so.

### Short regions that cannot be rebuilt no longer fail the run

A region under **10 residues** that DODO cannot rebuild is now reported and left with its input
coordinates, and the run still succeeds — `report.ok` stays `True` and the CLI exits `0`. A few
residues of AlphaFold geometry do not look wrong in a figure; it is the long extended regions that
do, and those still fail loudly.

`report.blocking_failures` and `report.tolerated_failures` split the two, and `report.failures`
remains everything that was not rebuilt.

### Build modes mean something different — check your figures

In 1.x a mode was **Ångströms per residue**: `normal` was 0.8 Å/residue, so a 100-residue IDR got
an 80 Å target and a 500-residue IDR got 400 Å. Real IDR end-to-end distance scales as roughly
N^0.55, not linearly, so a fixed per-residue value can only be right at one length.

In 2.0 a mode is a **multiplier on the ALBATROSS-predicted end-to-end distance**, so it is
length-independent. `expanded` means 1.3× whatever this particular sequence is predicted to be.

**The same mode name now gives a different structure, and the difference runs in both
directions.** Measured on a generic disordered composition, across all six modes:

| Region length | What 2.0 gives, relative to 1.x |
|---|---|
| ~11 residues | **2.0–2.4× larger** |
| ~20 residues | 1.7–1.9× larger |
| ~40 residues | 1.3–1.4× larger |
| ~75–95 residues | about the same (the crossover) |
| ~150 residues | 0.74–0.81× |
| ~280 residues | 0.56–0.62× |
| ~1200 residues | 0.31–0.34× |

Short regions — most loops, and many linkers — come out **bigger** than 1.x made them. Long tails
come out considerably smaller. If you are regenerating a figure from a 1.x script, expect it to
look different, and treat 2.0's version as the correct one.

The ratio also depends on composition now, which is the point: at 20 residues the same
`super_compact` request gives 1.65× for poly-glycine and 2.85× for poly-proline, because 2.0 asks
ALBATROSS what *this* sequence does instead of applying one number to every sequence.

Also: `mode='normal'` and `mode='predicted'` are now exact synonyms. In 1.x they differed, by
1.70× on a 282-residue IDR.

### Multi-model output is a real ensemble now

`-n 10` writes ten conformers. In 1.x every model shared one arrangement *and* essentially one
end-to-end distance, so it was an ensemble at fixed dimensions — the one thing an IDR ensemble
should not be. Each model now draws its own dimensions from the physical distribution.

Folded domains are positioned once and held fixed across every model, so flicking between frames
in a viewer shows only the disordered regions moving.

### Added

- **A single `dodo` command** with subcommands, replacing three separate scripts. `dodo rebuild`,
  `dodo fetch`, `dodo sequence`, `dodo regions`, and `dodo validate`.
- **`dodo validate`** checks bond lengths, steric clashes and CONECT records on any PDB or mmCIF
  file, not just DODO's output. Its reference values were measured from 105,299,848 bonds across
  23,587 AlphaFold structures. Findings on atoms DODO did not build are labelled as inherited, so
  a defect that arrived in your input is not reported as DODO's.
- **`dodo regions`** prints what DODO thinks your structure looks like, with the score profile and
  threshold behind every boundary, without rebuilding anything.
- **N, C and O are placed on rebuilt regions by default**, from the alpha carbons alone — new in
  2.0. Pass **`--no-backbone`** (or `backbone=False`) for alpha-carbon-only output; rebuilt alpha
  carbons are bit-identical either way, so the backbone is purely additive. Consecutive alpha
  carbons determine pseudo-dihedrals, which largely determine where the intervening peptide unit
  sits; DODO bins those at 20° against a table measured from 100 frames of all-atom IDR simulation,
  then rotates each unit about its CA–CA axis to settle bond angles, clashes and φ/ψ together.
  Rebuilding the table on 80 frames and scoring the 20 it had never seen, over 3,643 residues:
  **N 0.16 Å, C 0.22 Å, O 0.63 Å** mean error, with all four bond lengths exact by construction.
  The lookup is keyed on **five** consecutive alpha carbons — both pseudo-dihedrals flanking a unit
  — which beats the four-carbon form by 5.1% on C and 3.8% on N (held out over 19,302 units, paired
  95% CI excluding zero); units with no fifth carbon fall back to the four-carbon table. Folded
  domains are unaffected and keep every atom either way. Side chains are not built.
- **`dodo validate` still tells you which findings are the input's**, on backbone output as well as
  alpha-carbon-only output. Validating a *file* has no region metadata to read, so provenance is
  inferred from the structural signature DODO leaves — a residue carrying fewer atoms than its own
  type requires is one DODO built, and a finding touching none of those is inherited. Validating
  rebuilt dnmt3a reports `7 bond (3 inherited from the input)`: the three are AlphaFold's own
  distorted HIS613 ring, present in the input at identical values, and the other four are the
  seams, which DODO owns rather than blames on the input.
- **Reports instead of print statements.** `dodo.rebuild()` returns a `RebuildReport` with
  `.ok`, `.models`, `.failures`, `.assignments` and `.outcomes`. In 1.x these functions returned
  `None` and printed, so a script could not tell whether a region had actually been rebuilt.
- **Region identification strategies**, selectable with `-s`: `density` (DODO's original all-atom
  metric, and the default), `contact` (a CA-only alternative), and `plddt`. `auto` picks between
  the first two based on whether your input has side chains, and tells you.
- **Granular control from Python.** Assign regions however you like, then pass
  `strategy='preset'` and DODO builds exactly those:

  ```python
  import dodo
  structure = dodo.read_structure("model.pdb")
  dodo.assign_regions_from_spec(structure, {"A": [("idr", 1, 60), ("folded", 61, 912)]})
  report = dodo.rebuild(structure, strategy="preset")
  ```

  This replaces 1.x's `regions_dict=`, which took a separate description of the structure that
  could disagree with the real one.
- **mmCIF reading**, alongside PDB. Both accept `.gz`.
- **Reproducibility.** `seed=` / `--seed` makes output bit-identical.
- **`--ca-only`** and **`-b`/`--annotate-regions`**, which were 1.x's `-f/--no_FD_atoms` and
  `-b/--beta_for_FD_IDR`.

### Fixed

- **Backbone placement could put two atoms on top of each other at a region/domain seam.** The exact
  seam placement satisfies both bond lengths and was blind to everything else nearby, so it could
  land the placed atom on a folded one: measured on PTBP2 (Q9UKA9), a rebuilt N **0.797 Å** from a
  folded carbonyl oxygen, and on Q15642 an induced carbonyl oxygen on a folded atom 56 residues away
  in sequence — too far apart to be caught by a neighbouring-residue check. The exact placement is
  now verified before it is accepted (against every atom nearby, against the oxygen it induces, and
  against the residue's own N–CA–C angle) and falls back to the obstacle-avoiding placement when it
  fails. Across all 117 corpus structures, backbone placement now introduces **zero** separations
  below the 1.00 Å floor.
- **Refinement could collapse a residue's own N–CA–C angle under clash pressure**, putting that
  residue's N and C 1.90 Å apart — which the validator correctly reports as two atoms on top of each
  other. Measured on Q8N8A8 and O60271 at ~79°. The soft angle term is deliberately weak (tightening
  it was measured to make placement worse), so it could lose to the much stronger clash term. Both
  the compiled and NumPy refinement paths now hold a hard 80–160° window on that angle.
- **A crowded region could fail a rebuild outright** with `GeometryError` from the compiled
  refinement kernel's neighbour cap. See the `MAX_NEIGHBOURS` note under performance below.
- **`--backbone` was more than twice as slow as it needed to be, for two separate reasons.** On p300
  it went from **23.6 s to 11.0 s** (the backbone pass itself from 21.6 s to 9.0 s), with placement
  accuracy unchanged — N 0.147 to 0.148 Å, C 0.150 to 0.152, O 0.448 to 0.451, N–CA–C spread 3.06 to
  3.09 degrees.

  The first was `numpy.cross`. It supports arbitrary axes and pays for that on every call with
  `normalize_axis_tuple` and `moveaxis` — Python-level bookkeeping that dwarfs six multiplies on
  small arrays. Profiling one 583-residue region found 219,816 calls costing 6.0 s cumulative, with
  1.3 million calls to `normalize_axis_tuple` beneath them validating an axis that was never anything
  but -1. Written out by component instead. `numpy.linalg.norm` had the same problem.

  The second was that **refinement never converged.** Candidates are drawn from a window narrowing as
  `180/(1+sweep)`, so the smallest non-zero move available is the candidate spacing — 15.65 degrees
  on the first sweep and still 0.52 on the thirtieth, against a 0.25 degree tolerance. The test
  `largest_move <= tolerance` therefore meant `largest_move == 0`: every long region burned the full
  30-sweep budget and reported `converged=False` even after the geometry had stopped changing.
  Convergence is now judged on the objective, which finishes a 583-residue region in 15 sweeps and a
  60-residue one in 7.

  A third change came from evaluating scipy alternatives at the shapes the code actually uses: the
  clash term now uses `scipy.spatial.distance.cdist` with `sqeuclidean` rather than a broadcast
  `numpy.linalg.norm`, which never materialises the `(n_points, n_others, 3)` intermediate. Measured
  1.6x faster at 2 neighbours and 4.3x at 50, bit-identical output, and 11.0 s to 9.9 s on p300.

  Two scipy options measured and **rejected**, both for the same reason the `np.cross` fix worked:
  the arrays here are tiny, so per-call overhead decides. Building a `cKDTree` per clash test costs
  16-25 us against 4 (it earns its keep choosing neighbour sets once per sweep, which is where it is
  used), and `scipy.spatial.transform.Rotation` on a 24-candidate batch is 11.8 us against 9.0 for
  the hand-rolled equivalent. A third, expressing refinement as a `scipy.optimize` problem, does not
  apply at all: the Ramachandran term is a binned table lookup, so the objective is piecewise
  constant and a gradient-based solver would see no gradient from it.

  One thing that looked like the obvious culprit and was **not**: the clash term rebuilt a `(3N, 3)`
  array with `np.vstack` on every call, roughly 137,000 times per structure. Fixing it — the live
  arrays are now views into one buffer — saved 0.2 s of 21.5. Worth keeping, but the profile, not the
  guess, found the real cost.

- **The progress bar stopped at the end of the alpha-carbon rebuild** and left the backbone pass
  unreported, which on p300 was most of the wait. It now accounts for both passes and names the
  region it is working on, so a single long region is distinguishable from a stall.



- **`--backbone` no longer damages geometry at region seams.** The seam bond itself stays strained —
  that is unavoidable and documented — but three successive attempts to place the seam carbon and its
  oxygen *by construction* each introduced a different collision, because each was blind to one more
  neighbour. Aiming the carbon at the anchor's nitrogen made the two collinear and threw the carbonyl
  onto an arbitrary axis (O 0.6 Å from its own alpha carbon). Holding the residue's own N–CA–C angle
  fixed that and left the oxygen 0.975 Å from that nitrogen. Placing the oxygen trans to the nitrogen
  fixed that and left it 0.96 Å from the anchor's CB.

  The placement now sweeps its one free azimuth and rejects any position that collides. Measured over
  three structures at five seeds each: **zero** atom pairs closer than the 1.00 Å floor, and **zero**
  damage to any residue's own internal geometry.

  Also fixed: `--backbone` no longer rewrites the atom order of residues it does not touch. It had
  been normalising every residue to N–CA–C–O order, which rewrote 864 of dnmt3a's 912 residues
  including every folded-domain one. Untouched residues now come through byte-identical, which is
  what "folded domains are returned exactly as they arrived" has to mean.

### Changed

- **A region whose build fails keeps its input coordinates, while the folded domains are
  repositioned regardless.** So when regions fail, the output is the original IDR geometry connected
  to domains that have moved out from under it. The run is correctly reported as failed and the CLI
  exits 2, but the file is still written.
  `tests/data/structures/testing_translation.pdb` reproduces this -- two domains moved, a
  280-residue terminal IDR left unbuilt, a written model with a 31.6 A gap. So check the exit
  status, not just the file. Whether a declared-failed rebuild should write at all is unresolved.
- **Backbone refinement now runs a compiled kernel by default**, `dodo.construct.backbone_kernel`,
  with `refine_backbone(backend="numpy")` keeping the pure-numpy path reachable. This adds `numba` as
  a base dependency; an import failure falls back to numpy silently rather than failing a rebuild.

  The refinement objective is cheap arithmetic evaluated an enormous number of times -- for a
  583-residue region, 397 units x 25 candidates x 15 sweeps -- and profiling put about 60% of the
  time in numpy's per-call dispatch rather than the arithmetic. Batching can only amortise that, never
  remove it: batching ten models recovered 3.6x of a theoretical 10x, and batching over peptide units
  cost more in extra sweeps than it saved. Compiling removes it outright. Measured on a 398-residue
  region, **2.47 s to 0.35 s (7.1x)**.

  Three things are compiled, because each one became the bottleneck as soon as the previous was fixed:

  - The **sweep loop**, 16x on its own.
  - The **objective**, once the sweep was fast enough that *reporting* it dominated. `energy_before`
    and `energy_after` are read once each per region, and the numpy scorer pushed every peptide unit
    through `score_candidates`: over p300's six regions, 12 calls cost 0.562 s against the compiled
    sweep's 0.204 s of self time -- two readings cost 2.8x the optimisation they described. Compiled
    they cost 0.026 s, taking refinement **1.92x**.
  - The **neighbour search**, which was a `cKDTree` built and queried in Python in the middle of an
    otherwise compiled path, at 0.199 s against the sweep's 0.204 s. A uniform cell list in nopython
    mode costs 0.011 s, **18x**, taking `refine_region` a further **1.80x**. `sweep_region` is
    unchanged across that change (68 calls, 0.219 vs 0.221 s), which is what shows only the search
    moved.

  End to end, a backbone rebuild now takes about **0.3 s on dnmt3a and 1.7 s on p300** (4 models of
  p300: 5.2 s), against 23.6 s for p300 before any of this work. Seam counts, intra-residue damage
  and impossible pairs are unchanged across three seeds on all three structures.

  Most of what remained after the compiled kernel landed was not the kernel at all. Profiling put
  **78%** of backbone wall time in the coupled-clash polish, which was still pure Python: a KD-tree
  query per atom per candidate, and a per-candidate Python placement-and-score loop. Precomputing
  each moved atom's neighbours once per group (the static cloud is fixed during a group's search,
  and each atom rides a bounded circle, so the shell is a provable superset) and evaluating a whole
  azimuth grid as array operations took that polish from **6.16 s to 0.40 s, 15x**, for
  **bit-identical** output — worst coordinate difference 0.0 over the corpus at four seeds.

  The compiled neighbour search is also *better*, not merely faster, in two ways. It is a true radius
  search, where `cKDTree.query(k=MAX_NEIGHBOURS)` can only ever examine that many nearest points; and
  it counts qualifying neighbours without a ceiling, so the cap is checked rather than trusted. That
  found a real latent bug in the guard it replaced: the old one raised only if the *last* nearest
  point survived the covalent-separation filter, so a constructed case kept 47 of 59 qualifying
  neighbours and **silently dropped 12** without raising. A crowded region that previously completed
  with a quietly truncated objective is now detected rather than ignored.

  `MAX_NEIGHBOURS` is **96**, and the cap is a checked property rather than an assumption. It was
  48, chosen against a worst measured count of 43 over a full p300 rebuild — which looked like
  headroom until the 117-structure corpus ran with the backbone on and a crowded region of Q9C000
  put **54** atoms in one unit's shell. Rather than trade one arbitrary number for another, a region
  that exceeds even 96 now falls back to the uncapped NumPy path instead of failing the rebuild;
  an explicit `backend="numba"` still raises, because the caller asked for the kernel by name.

  Two levers were measured and rejected, recorded so they are not re-litigated. `fastmath=True` is
  worth 1.14x on the sweep loop but only **1.2%** on a whole rebuild, and its `nnan` licence lets the
  `acos` domain clamp swallow a NaN -- `_angle` returns `nan` with it off and `180.0` with it on --
  which is not worth 1.2% in a geometry kernel. Intel **SVML is unavailable** rather than unhelpful:
  no `intel-cmplr-lib-rt` wheel exists for arm64 on any platform. A faster **threading layer** is a
  dead end too, with `omp` measuring 0.86x against the default `workqueue`; neither matters while this
  code is serial. `prange` over the *model* axis does scale (16.4x at 18 models) but needs the pipeline
  to collect every model's CA trace before one batched pass, so it is not wired up.

  Equivalence is established rather than assumed. On byte-identical inputs the two agree bit-for-bit
  on C, N and O placement, on the Ramachandran term and on the clash term, and to 1e-13 on the two
  N-CA-C angle terms where `math.acos` and `np.arccos` round differently -- and both select the same
  candidate. Accuracy against all-atom truth is equivalent (N 0.148 vs 0.139 A, O 0.451 vs 0.445,
  identical angle statistics, zero clashes, bonds exact).

  **Output is not bit-identical between backends.** That 1e-13 eventually flips a nearly balanced
  decision, after which a greedy search diverges to different -- measurably equivalent -- coordinates.
  It was never bit-identical across numpy versions or BLAS builds either, for the same reason.

  Worth recording for anyone who revisits this: a suspected bug here took four failed attempts to
  chase, and every one compared the END STATES of two independent runs. A greedy search amplifies any
  perturbation chaotically, so that comparison cannot separate a cause from its consequences. The
  pure-function comparison on identical inputs settled it in one shot, and there was no bug. The test
  that does it is in `tests/unit/test_backbone_kernel.py`.



- **ALBATROSS predictions are cached on disk**, keyed by sequence hash and sparrow version. The
  first prediction in a process imports sparrow, parrot and torch, which measured **1.57 s against
  0.17 s of actual rebuilding** for a 912-residue structure — so on a repeat run the import, not
  the geometry, was what you waited for. A cache hit returns before sparrow is imported at all:
  1.47 s to 0.18 s, with torch never loaded. Entries are about 116 bytes; set `DODO_NO_CACHE=1` to
  opt out.

- **Obstacle clash-testing is pruned.** Every candidate position for one residue sits exactly one
  bond length from a single centre, so when that centre clears every obstacle by more than
  `CA_CA_BOND_LENGTH + CA_CLASH_DISTANCE` a single one-point query replaces 355. Measured over 4
  seeds per structure: `_nearest_obstacle_distance` on p300 fell from 1.683 s to 0.088 s (19x),
  and whole rebuilds from 5.68 s to 3.46 s on arf19 and 19.6 s to 10.4 s for p300 at 10 models.
  Output is **byte-identical**, verified across 55 outputs spanning PDB, mmCIF, a 55-chain
  assembly, `--backbone`, `--ca-only` and multi-model runs.

  Two prunings that did **not** pay, measured rather than assumed: spatial reachability pruning is
  inert on the long regions where the time actually goes (a 150-residue region already reaches
  564 A, wider than any test structure), and burial pruning qualified only 1.2-1.4% of atoms while
  costing 240-530 ms — several times the entire remaining obstacle budget.


- **Rebuilds are 1.6–4× faster**, with no change to what is built. Measured end to end:

  | | before | after |
  |---|---|---|
  | arf19 | 9.7 s | 2.4 s |
  | p300 | 9.8 s | 4.8 s |
  | arf19, 10 models | 63.8 s | 21.7 s |
  | p300, 10 models | 34.2 s | 21.1 s |

  Around 90% of a rebuild is region building, and almost all of that is one KD-tree clash query
  per residue. Two changes: candidate positions per step dropped from 710 to 355, and those
  queries now run single-threaded. Handing scipy an 18-thread pool for a 710-point query cost more
  than it saved — measured, thread setup is a roughly fixed 0.3 ms and parallelism only pays above
  about 2,000 points.

  The candidate count buys retries rather than per-step correctness: 99.9% of steps have every
  candidate geometrically valid before the clash test. Halving it was validated on the full
  117-structure corpus (identical results) and over 5 structures × 4 seeds (same failure count).
  Beware single-seed timings when revisiting this — retry counts are wildly stochastic, and one
  structure measured 3.4× on one seed and 1.3× averaged.

  Two things that looked like wins and were not, both measured: replacing the per-step KD-tree with
  brute-force numpy is **10× slower**, and widening the conformer batch does not help because
  per-conformer cost is flat in batch width (0.35–0.46 s per conformer from 1 row to 32).

- **A progress bar**, on stderr, weighted by residues. Shown when stderr is a terminal, suppressed
  when piped or under `-q`. `progress=True`/`False` overrides from the API. This adds `tqdm` as a
  dependency — pure Python, no dependencies of its own, ~60 kB.

- **Quieter notes.** The `strategy auto-selected as density` note is gone: density on all-atom input
  is the expected path that nearly every run takes, and printing it every time trained people to
  skip the notes. It still fires for the CA-only fallback to `contact`, where the choice changes
  where boundaries land. Short folded blocks below the length minimum are now one summary line
  instead of one line per block.

- **`dodo fetch --no-cache`** downloads to a temporary file instead of the per-user cache. Caching
  stays on by default: measured over 259 cached AlphaFold models, the mean file is 0.25 MB and the
  largest 1.51 MB.


- **Downloads no longer break.** 1.x built AlphaFold URLs with a hardcoded `model_v4`, which now
  returns HTTP 404 for every accession. The model version is resolved from the AFDB API instead.
- **Rebuilt regions land where the sequence says.** Folded domains are now translated and rotated
  so that the gap between them matches the predicted end-to-end distance of the linker joining
  them. AlphaFold packs domains 2–3.6× closer than the linker predicts, so without this the
  linker is built into a gap unrelated to its sequence.
- **Structures render correctly in VMD and PyMOL.** CONECT records are written by default, atom
  names occupy the correct columns, and the element column is populated. 1.x's CA-only output was
  read by MDTraj as *calcium* atoms.
- **Failures are explicit.** A region that cannot be built is reported with a reason and keeps its
  input coordinates. 1.x could return coordinate arrays of exact `(0, 0, 0)` rows, or NaN.
- **Requires Python 3.10+, numpy 2.0+, scipy 1.13+ and numba.** numba is a base dependency, not an
  extra: it compiles the backbone refinement kernel, which is on by default.
- Exceptions all derive from `DodoError`. Bad arguments raise `InvalidParameterError`, which is
  both a `DodoError` and a `ValueError`, so either `except` clause works.

### Removed

- `build.pdb_from_name()`, `build.pdb_from_pdb()`, `build.pdb_from_sequence()`, and the
  `pdb-from-name` / `pdb-from-pdb` / `pdb-from-sequence` commands. See the migration table.
- `regions_dict=`. Use `strategy='preset'` as shown above.
- `attempts_per_region` / `attempts_per_coord` / `attempts_per_residue` / `attempts_per_IDR`.
  Retry budgets are internal and reported when a region fails.
- `graph=True` and `plot_structure()`. No plotting is included.
- `linear_placement`. Folded domains are placed by facing the next domain and perturbing, then
  rejecting arrangements that clash.
- `just_fds=True`. Not reimplemented.
- `--all-atom` and `--sidechains`. **Side-chain** building for rebuilt regions is not in 2.0.
  Backbone building is, under the new `--backbone` flag — see Added.

  Note that neither flag did what its name suggests. 1.x's "all atom" meant keeping the atoms of
  the regions DODO does **not** rebuild, and 2.0 does that unconditionally, so folded domains keep
  full atomic detail regardless. And in the 2.0 development tree `rebuild(all_atom=True)` was a
  silent no-op — byte-identical output — while `build_from_sequence(all_atom=True)` raised. Both
  are gone rather than fixed, because a flag that lies is worse than no flag.

### Known limitations

Stated with measurements, over a 117-structure corpus stratified by region topology.

- **One marginal steric contact across 117 structures**, at 3.02 Å against a 3.20 Å limit. Nothing
  below the 1.00 Å floor below which no real bond exists, and nothing that was not already in the
  input. Attaching a rebuilt region to the structure requires exempting its *anchors* from clash
  checking; DODO exempts only their alpha carbons by default and falls back to their backbone only
  for a region that would otherwise not be built at all — reporting it when it does. Measured over
  the corpus, that two-pass approach gives 1 contact and 3 unbuilt regions, against 21 and 4 for
  always exempting the backbone, and 1 and 8 for never exempting it.
- **Three regions out of 117 structures are left unbuilt, and two of those are input defects**: one
  file has two fixed residues 3.04 Å apart, already closer than the clash distance DODO must
  satisfy, and one contains a chain break with consecutive alpha carbons 5.26 Å apart. The third is
  a genuine failure on a 7-residue terminal tail. Every one keeps its input coordinates and is
  reported with a reason.
- **Side chains are not built.** Rebuilt regions get a complete N/CA/C/O backbone and nothing more.
  Regions DODO does *not* rebuild keep every atom they arrived with, side chains included.
- **Backbone placement introduces a small number of marginal steric contacts.** The alpha carbons
  cannot self-intersect, but the atoms hung off them can, and one azimuth per peptide unit is not
  always enough freedom to avoid it — the carbonyl oxygen especially, which is fully determined by
  its own carbon and so cannot be moved independently. Refinement scores each unit against the
  folded domains, the rest of its own region, and the regions placed before it, and a joint
  two-unit polish then clears the coupled cases single-azimuth descent cannot escape. Measured
  against CA-only output at seed 0: dnmt3a 2 contacts, arf19 0, p300 3; over nine AlphaFold
  structures, six are clean and the worst is 4. Overlaps run to 0.8 Å at worst, and **nothing is
  near the 1.00 Å floor** below which no real bond exists — that is asserted absolutely, across all
  117 corpus structures.
- **A strained peptide bond is left at every seam, and this cannot be fixed by
  construction.** A peptide unit reaches at most 2.854 Å from an alpha carbon to the nitrogen it
  bonds to. Where a rebuilt region meets a folded domain, that domain's nitrogen still points toward
  where the region ran in AlphaFold's model, because DODO moved the domain rigidly and redrew the
  region beneath it. Measured across 17 seams in three structures, the rebuilt alpha carbon sits
  **3.2–4.5 Å** from that nitrogen — so the bond is not merely hard to get right, it is
  unsatisfiable.

  DODO aims the carbon at the nitrogen and leaves the bond long (about 2.2 Å against an ideal 1.33)
  rather than writing two atoms into the same space, and **labels and reports** every one: they
  appear on `RebuildReport.backbone_seams`, in the run summary, and as `kind="seam"` /
  `provenance="seam"` findings from the bond validator — reported, but not blamed on the rebuild,
  whose own bond geometry is clean. On dnmt3a that is 4–6 bonds out of 911. **Zero impossible
  contacts are introduced.** Building from sequence has no seams and is completely clean.

  Two ways to close it were tried and both rejected on measurement. Leaving the seam residue
  un-rebuilt does make the bond reachable (2.45–2.52 Å at all 17 seams), but closing onto it means
  re-placing its nitrogen, and that residue's side chain was built around where its nitrogen used to
  be: the new one is driven into the residue's own CB (1.405 Å against a correct 2.45) and, for
  proline, snaps the ring bond to CD (3.444 Å against 1.47). Constraining the walk so its closing
  alpha carbon lands within peptide reach does not work either: over 83 closures across three
  structures and three seeds, only 22 have an in-reach point that also satisfies the CA–CA–CA
  pseudo-angle window — and those already close — while 40 never bring the closure circle within
  reach at all and 21 could only reach by kinking that angle. The strained seam is what independent
  rigid-body repositioning and chain rebuilding cost.
- **Reported failures are honest, not silent.** A region DODO cannot build keeps its input
  coordinates and appears in `report.failures` with a reason. It is never replaced with degenerate
  output — 1.x could return coordinate arrays of exact `(0, 0, 0)` rows, or NaN.

### Migrating from 1.x

| 1.x | 2.0 |
|---|---|
| `build.pdb_from_pdb(path, out_path=...)` | `dodo.rebuild(path)` then `dodo.write_pdb(...)` |
| `build.pdb_from_name(name, out_path=...)` | `dodo.fetch_alphafold(accession)` then `dodo.rebuild(...)` |
| `build.pdb_from_sequence(seq, out_path=...)` | `dodo.build_from_sequence(seq)` |
| `pdb-from-pdb in.pdb -o out.pdb` | `dodo rebuild in.pdb -o out.pdb` |
| `pdb-from-name "human p53" -o out.pdb` | `dodo fetch "human p53" -o out.pdb` |
| `pdb-from-sequence SEQ -o out.pdb` | `dodo sequence SEQ -o out.pdb` |
| `-o`, `--out_path` | `-o`, `--out` |
| `-c`, `--no_CONECT_lines` | `--no-conect` (no short form) |
| `-s`, `--silent` | `-q`, `--quiet` (`-s` now means `--strategy`) |
| `-u`, `--use_metapredict` | removed; `density` is the default and now takes 7 ms |
| `-f`, `--no_FD_atoms` | `--ca-only` |
| `-b`, `--beta_for_FD_IDR` | `-b`, `--annotate-regions` |
| `num_models=N` | `n_models=N`, or `-n N` |
| `mode='normal'` (0.8 Å/residue) | `mode='predicted'` — **different result**, see above |
| `use_metapredict=True` | removed; see above |
| `regions_dict={...}` | `assign_regions_from_spec(...)` + `strategy='preset'` |
| `beta_for_FD_IDR=True` | `write_pdb(..., annotate_regions=True)` |
| functions returned `None` and printed | functions return a `RebuildReport` |
| `except dodoException` | `except DodoError` |

Multi-word protein names now need quoting: `dodo fetch "human p53"`.

Exit status: `0` on success, `2` if some regions could not be built, `1` on error.
