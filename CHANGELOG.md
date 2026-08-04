# Changelog

This changelog starts at 2.0. Earlier releases predate it.

It records what changes **for you as a user** — the commands and functions you call, the
structures that come out, and the things you have to do differently. Internal restructuring is
not listed, however large, unless you can see it from outside.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[semantic versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] — unreleased

DODO 2.0 is a rewrite. **It breaks the 1.x API deliberately**, so read the migration notes below
before upgrading a script.

The scientific behaviour changed too, in ways you can see in your output. The two that matter
most: build modes now mean something different, and a multi-model run is now a real ensemble.

### Installing is simpler: one package, one extra

```bash
pip install idptools-dodo                 # everything you need
pip install "idptools-dodo[starling]"     # adds STARLING ensembles (~2.4 GB of weights)
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

**`[all]` is gone.** With one extra left there is nothing to union.

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
- **`--backbone` / `backbone=True` places N, C and O on rebuilt regions**, from the alpha carbons
  alone. New in 2.0 and **opt-in**, because of the seams (see Known limitations). Four consecutive
  alpha carbons determine a pseudo-dihedral, which largely determines where the intervening peptide
  unit sits; DODO bins that angle at 20° against a table measured from 100 frames of all-atom IDR
  simulation, then rotates each unit about its CA–CA axis to settle bond angles, clashes and φ/ψ
  together. Rebuilding the table on 80 frames and scoring the 20 it had never seen, over 3,643
  residues: **N 0.16 Å, C 0.22 Å, O 0.63 Å** mean error, with all four bond lengths exact by
  construction. Folded domains are unaffected and keep every atom either way. Side chains are not
  built.
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

- **The STARLING engine did not work at all**, and the error blamed the wrong thing. Every region
  failed with "could not find coordinates on the object STARLING returned (Ensemble) ... this is an
  API mismatch". The real cause was a missing argument: `starling.generate()` defaults to
  `return_structures=False`, which returns an ensemble carrying only **distance maps**. Four
  separate problems, each of which would have been worse had it silently "worked":

  - `return_structures=True` is now passed.
  - Coordinates are read from `ensemble.trajectory.traj.xyz`. `trajectory` is a soursop `SSProtein`
    and has no `.xyz` of its own, which is what the old probe looked for.
  - **They are in nanometres.** Had extraction succeeded without conversion, DODO would have
    written structures ten times too small.
  - STARLING keys its return dict `sequence_1`, not by the sequence string.

- **Regions longer than STARLING's 380-residue cap no longer fail.** `--engine starling` now wraps
  `HierarchicalEngine` automatically. Previously a 401-residue linker errored with advice to wrap
  the engine yourself, which is not something a user of `dodo rebuild` can act on.

- **STARLING's output is now repaired before it is judged.** Its coordinates come from MDS on a
  predicted distance map and are not a physically valid backbone: measured, virtual CA-CA bonds
  span 1.81-4.86 A with a median of 3.36, **47.6% of them more than 0.5 A off ideal**, and
  pseudo-angles reach 5.2 degrees — a vertex where the chain doubles back within one step. Applying
  the physical screen to raw output rejected **100 of 100** conformers.

  Bond and pseudo-angle repair now run before the physical screen, and the reordering is safe
  because repair barely moves the ensemble: end-to-end distance changes by +0.1-0.25% and radius of
  gyration by +0.25-1.14%. Gross breakage is still caught on raw coordinates first, so a chain the
  model genuinely lost cannot be quietly repaired into looking fine. Measured acceptance went from
  0/20 to **20/20** at 38 residues, and a 200-residue region went from raising outright to building
  every requested conformer.

  Remaining honest limitation: at ~200 residues only about 5 of 20 conformers survive, and the
  cause is **internal clashes already present in STARLING's raw output** — 18 of 28 raw conformers
  clash before DODO touches them. Repair does not introduce them; no conformer went from clash-free
  to clashing.

- **The STARLING isolation warning is emitted once per run, not once per model**, and is 392
  characters instead of about 800. A ten-model run printed it ten times, which is how a warning
  worth reading becomes one nobody reads.

### Added

- **`--domain-placement conformer`** (`domain_placement="conformer"`), opt-in and requiring
  `--engine starling`. Positions the folded domains to match each generated conformer instead of
  generating conformers to match a predicted domain separation.

  The default path predicts a linker's end-to-end distance, moves the domains to it, and selects
  conformers that match — which is right for the walk engine, because it builds *to* a dimension.
  STARLING does not: it samples a distribution. Measured on dnmt3a, its end-to-end distances scatter
  with a standard deviation of 41% of the mean against a 5% selection tolerance, so about a tenth of
  conformers survive and they form a narrow band at the mean. The two models agree on the mean
  (133.0 Å against 136.6 Å on one region); what selection throws away is the spread.

  Measured on dnmt3a over three models: regions built went from 6/9 to **9/9**, and the achieved
  end-to-end spread from **sd 10.8 Å to sd 44.9 Å**. A conformer is now rejected only if the domain
  placement it implies collides with something already placed.

  Two consequences: the domains differ between models rather than sharing one arrangement, and
  `--mode` has nothing to multiply for the affected regions. Loops are unaffected — they are pinned
  inside a single domain, so nothing can be repositioned for them.

### Changed

- **STARLING conformers are now filtered before they are repaired, and generated in rounds.** Repair
  is the expensive step — 3.8 s for 100 conformers of 380 residues — and two defects survive it: a
  chain the model lost, and one that already collides with itself. Both are now detected by checks
  that cost nothing, so hopeless conformers are dropped before any repair is paid for. The clash
  threshold is calibrated never to discard a conformer repair would have saved: measured over 60
  conformers, every one that ultimately passed the physical screen had a raw worst contact of at
  least 3.536 Å, against a 2.90 Å cutoff.

  Survivors are repaired best-first, and a repair that moves any alpha carbon more than 4.0 Å is
  rejected rather than accepted — repair is meant to keep a conformer recognisably STARLING's
  (displacement is 0.84 Å median, 2.10 Å at the 99th percentile), but a severe angle violation beside
  a free terminus has nothing to anchor it and one conformer in twelve had an atom travel 21.8 Å.
  If a round does not yield enough usable conformers, another is generated, up to three.

  One thing that looked like the obvious saving and was **wrong**: stopping repair as soon as enough
  conformers were in hand. Selection downstream needs *spread* in end-to-end distance, not a count —
  truncating the pool at 2 left both survivors the wrong size and a region failed on dimension by
  12.7 Å having repaired only 2 of 16 candidates. All survivors are repaired; the saving comes from
  the free elimination step instead.



- **STARLING ensembles are generated once per sequence and reused across models.** STARLING
  conditions on sequence alone, so the ensemble it returns for a region is identical every time;
  what differs per model is which conformer is selected and where it is placed. A 2-model dnmt3a
  rebuild ran 6 diffusion-plus-MDS generations for 3 regions before this.

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
- **Requires Python 3.10+, numpy 2.0+ and scipy 1.13+.**
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
- **Rebuilt regions are alpha-carbon only by default.** Deliberate. Regions DODO does *not* rebuild
  keep every atom. `--backbone` adds N, C and O to the rebuilt regions; side chains are still not
  built.
- **`--backbone` introduces a small number of marginal steric contacts.** The alpha carbons cannot
  self-intersect, but the atoms hung off them can, and one azimuth per peptide unit is not always
  enough freedom to avoid it. Refinement scores each unit against the folded domains, the rest of
  its own region, and the regions placed before it. Measured against CA-only output over three
  seeds: dnmt3a 0 → 1–3 contacts, arf19 0 → 1–15, p300 2 → 13–31. The spread is the honest answer,
  turning on whether a conformer folds a hairpin back on itself. Worst overlap across all of them is
  0.48 Å; none is near the 1.00 Å floor below which no real bond exists.
- **`--backbone` leaves a strained peptide bond at every seam, and this cannot be fixed by
  construction.** A peptide unit reaches at most 2.854 Å from an alpha carbon to the nitrogen it
  bonds to. Where a rebuilt region meets a folded domain, that domain's nitrogen still points toward
  where the region ran in AlphaFold's model, because DODO moved the domain rigidly and redrew the
  region beneath it. Measured across 17 seams in three structures, the rebuilt alpha carbon sits
  **3.2–4.5 Å** from that nitrogen — so the bond is not merely hard to get right, it is
  unsatisfiable.

  DODO aims the carbon at the nitrogen and leaves the bond long (about 2.2 Å against an ideal 1.33)
  rather than writing two atoms into the same space, and reports every seam. On dnmt3a that is 4–6
  bonds out of 911. **Zero impossible contacts are introduced**, measured over three structures at
  three seeds each. Building from sequence has no seams and is completely clean.

  Leaving the seam residue un-rebuilt was tried and rejected. Its alpha carbon being input geometry
  does make the bond reachable — 2.45–2.52 Å at all 17 seams — but closing onto it requires
  re-placing its nitrogen, and that residue's side chain was built around where its nitrogen used to
  be. The new one is driven into the residue's own CB (1.405 Å against a correct 2.45) and, for
  proline, snaps the ring bond to CD (3.444 Å against 1.47). The real fix is to constrain the walk so
  its closing alpha carbon lands within peptide reach, which needs constraints reaching further back
  than the closing step, and is future work.
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
