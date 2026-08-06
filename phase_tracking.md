# DODO v2 refactor — phase tracking

Living document tracking objectives and verified success for the v2 refactor.
Decision record: **refactor, not restart** (v2 is well-tested and scientifically correct;
the work is a lean-out plus a bounded vectorization refactor). See the plan below.

Verification bar for every phase (all must hold before a phase is "done"):
- `python -m pytest` green (0 failed, 0 errored) in the `hh1` conda env.
- `import dodo` works; `dodo rebuild` on the dnmt3a fixture still produces valid output.
- `mypy --strict` and `ruff` clean (the branch's standing invariant).
- The default CA path (`engine="walk"`, folded-domain repositioning on) is behaviorally unchanged.

Baseline before Phase 0: 2202 passed, 116 skipped, 0 failed (full suite, ~5m41s).

---

## Phase 0 — Delete dead code + excise the STARLING web  ✅ DONE (2026-08-05)

**Objective:** reclaim ~25% of the tree by removing code that never runs in the default
CA path, with zero behavior change to that path. STARLING returns later as a point update;
it stays recoverable from git history at HEAD (`51b2ad6`). No branches created (Ryan owns git).

Decisions:
- KEEP the engine seam (`--engine`, `engine=` param, `_make_engine`) but walk-only, so STARLING
  re-adds as a small diff.
- REMOVE the `domain_placement="conformer"` mode entirely (STARLING-only inverse placement).

### 0a. Delete dead `backbone.py` (superseded by `ca_backbone.py`; only tests import it)
- [x] Delete `src/dodo/construct/backbone.py` (2874 lines)
- [x] Delete `tests/unit/test_backbone.py` (imports from the dead module)
- [x] Fix `tests/unit/test_constants.py` (imports `_N_CA_C_TOLERANCE`, `max_reconstructable_ca_angle` from `construct.backbone`)
- [x] Confirm no `src/` importers remain

### 0b. Excise the STARLING web
Delete files:
- [x] `src/dodo/engines/starling.py` (2262)
- [x] `src/dodo/engines/hierarchical.py` (1399)
- [x] `src/dodo/geometry/regularize.py` (722; imported only by starling.py)
Delete/trim tests:
- [x] `tests/unit/test_engines_starling.py`
- [x] `tests/unit/test_regularize.py`
- [x] hierarchical test if present
- [x] STARLING/conformer tests in `tests/unit/test_pipeline.py`, `tests/unit/test_documentation.py`
Rewire src (remove STARLING/conformer paths):
- [x] `pipeline.py`: drop `_rebuild_from_conformers`, `_warn_starling_isolation`, the `starling` branch of `_make_engine`, and all `conformer_placed`/`domain_placement` branching
- [x] `place.py`: drop `place_domain_after` (imports `place_between_anchors` from starling)
- [x] `cli.py`: `_ENGINE_CHOICES = ("walk",)`, simplify `--engine` help, remove `--domain-placement` flag + its call-site kwargs
- [x] `constants.py`: remove `STARLING_MAX_LENGTH`, segment-overlap, and the regularize reference comment
- [x] `exceptions.py`: drop STARLING docstring examples
- [x] `pyproject.toml`: remove the `[starling]` extra

### 0c. Verify
- [x] Full pytest green
- [x] `import dodo` + `dodo rebuild` on dnmt3a fixture valid
- [x] mypy --strict + ruff clean
- [x] Record new line count and test totals here

**Result (2026-08-05):**
- Deleted 7 files: `backbone.py`, `starling.py`, `hierarchical.py`, `regularize.py` (src) +
  `test_backbone.py`, `test_engines_starling.py`, `test_regularize.py`.
- Rewired `pipeline.py`, `place.py`, `cli.py`, `constants.py`, `pyproject.toml` and trimmed
  `test_constants.py` / `test_pipeline.py` / `test_documentation.py`.
- Kept the `engine` seam (`--engine`, walk-only) for STARLING's later return; removed the
  `domain_placement="conformer"` mode entirely.
- **src: 29,901 → 22,223 lines (−7,678, −25.7%).** Net across src+tests: −12,951 lines (20 insertions).
- Verification: ruff clean, `ruff format` clean, `mypy --strict` clean (34 files), `import dodo` OK.
  Full suite **1820 passed, 116 skipped, 0 failed** (was 2202; the −382 are the deleted tests).
  Real `dodo rebuild dnmt3a` → exit 0, 4970 atoms. `--engine starling` rejected by argparse;
  `--domain-placement` gone from help.
- STARLING remains recoverable from git history at HEAD (`51b2ad6`). A few passing docstring
  mentions of STARLING remain in `base.py` / `backbone_refine.py` / `exceptions.py` — cosmetic,
  deferred to the Phase 4 docs pass. No branch was created (Ryan owns git).

---

## Phase 1 — Wire up batching that already exists  ⚠️ REASSESSED (2026-08-05)
Finding (verified in code): the two items are NOT standalone wire-ups.
- **Rebuild conformer batch is BLOCKED by per-model inter-region obstacles.** `placed_atom_mask()`
  (structure.py:610) includes each region as it is built (pipeline.py:365/:674), so within a model
  region B avoids region A's *this-model* coords. `generate(request, obstacles, rng)` takes ONE
  shared obstacle array per call, so batching conformers across models can't give conformer m the
  right per-model obstacle. Correct batching needs per-conformer obstacles → **this is Phase 3**
  (grow-then-filter with a global clash pass), not a wire-up.
- **Backbone refinement batch (16.4×)** needs a batched numba kernel + batched placement entry
  point → that IS the **Phase 2** backbone-vectorization work.
- `build_from_sequence` / `dodo sequence -n N` already batches (`n_conformations=n_models` in one
  call); nothing to wire up there, though the walk's two per-conformer inner loops still slow it.

Resolution: fold item 1 into Phase 3 and item 2 into Phase 2. Next clean win = Phase 2.

## Phase 2 — Vectorize deterministic backbone placement  ✅ DONE (2026-08-05)
Rewrote `ca_backbone.backbone_from_ca`'s per-residue Python loop as leading-dim-agnostic batched
numpy (`_backbone_atoms`, handles `(n,3)` and `(B,n,3)`), plus batched helpers (`_frames_rows`,
`_dihedral_rows`, `_two_spheres_rows`, `_carbonyl_rows`, terminals). Public API + `CABackboneResult`
unchanged; scalar seam-helpers kept for `add_backbone_to_rebuilt`.

**Verification (evidence):**
- Equivalence to the old per-residue impl across 16 traces (edge/random/compact, n=2..200):
  **max abs coord diff 1.78e-14** (machine epsilon). `source`/`notes` identical.
- Batched `(B=50,n=100)` vs looped-single: **exactly 0.0** — batch axis correct.
- Speed: n=583 single trace **32.20 ms → 0.37 ms (87.7×)**; 50×100 conformers loop 10.79 ms →
  one batched call 1.69 ms (**6.4×**).
- ruff + `ruff format` clean; `mypy --strict` clean (34 files); full suite **1820 passed, 0 failed**
  (one error-message contract test adjusted to keep the "zero-length" vocabulary).

**Deferred (the "backbone batch", = old Phase 1 item 2) — DECISION 2026-08-05: not doing now.**
Investigated in depth: placement is already batched (Phase 2) and is now ~0.37 ms/region, so the
only remaining win is the refinement's measured 16.4×, which is numba `prange` thread-parallelism
over the model axis (kernel header lines 6–8, 91–92). Delivering it needs the whole convergence loop
compiled as `njit(parallel=True)` (16-arg `sweep_region` + per-sweep grid rebuilds, stacked across
models) plus restructuring `rebuild()` to run backbone as a post-loop batched pass. Payoff is narrow
— opt-in `--backbone` × multi-model × multi-core only; nothing for the default path. Ryan chose to
prioritize Phase 3. Left as a documented low-priority follow-up.

## Phase 3 — Hybrid vectorized IDR generator  ⏳ IN PROGRESS (2026-08-05)
Staged so each slice is verifiable and the tested `walk.py` pipeline stays stable until integration.
- **Slice 1 ✅ DONE (2026-08-05): batched free growth + self-avoidance filter.** New module
  `src/dodo/engines/batch_walk.py` (`grow_free_batch`, `self_avoiding_mask`, `end_to_end_distances`,
  `radii_of_gyration`) + `tests/unit/test_batch_walk.py` (23 tests). Additive; nothing imports it yet.
  ruff/mypy clean; full suite 1843 passed / 0 failed.
  MEASURED (the data that drives the rest):
  - Geometry exact by construction: bond dev ≤1.5e-14 Å; every pseudo-angle in [91,161].
  - Natural Re vs analytical Flory target — ratio 0.91 (n=10), 0.99 (20), 1.05 (50), 1.12 (100),
    1.17 (200); CV ~0.28–0.38. So Slice 2 target-steering is mostly SELECTION near the target, not
    heavy biasing (for normal modes).
  - Self-avoidance survival collapses with length: 97% (n=10), 85% (20), 48% (50), 17% (100),
    2% (200), 0.1% (380). → grow-then-filter is great for short/medium, impractical for long free
    chains. Confirms the hybrid: filter-primary + careful fallback for long/crowded.
  - Speed: 21 µs/chain (1000×100-mer in 21 ms) vs the current walk's 12,649 ms for 1000 conformers.
- **Slice 2 ✅ DONE (2026-08-05): target-steering + crowded-environment survival.** Added
  `steer_to_target`, `generate_free_batch` (+`FreeBatchResult`), and `clears_obstacles_mask` to
  `batch_walk.py` (+11 tests, 34 total). Steering = grow-and-filter a pool, draw physical targets
  via the walk's own `sample_end_to_end_targets`, match each to the nearest-Re survivor. ruff/mypy
  clean; full suite 1854 passed / 0 failed.
  MEASURED:
  - Steering accuracy across compact/normal/expanded × n=20–200: achieved-mean ratio **0.975–1.003**,
    CV preserved at **0.35–0.44** (no ensemble collapse), reachability 0.97–1.00. → hitting mode
    targets is a solved SELECTION problem, not biasing.
  - Attrition remains the only real constraint: at n=200 survivors ≈223, reuse 0.44–0.62 (pool thin)
    → long regions need big oversample or fallback.
  - Crowding by ONE adjacent domain is moderate: a r=15 Å blob 5 Å from the anchor halves survival
    (39%→19%), recovering toward free by 30 Å. Manageable via oversample + first-step biasing; the
    harder interior/between-two-domains case is Slice 3.
- **Slice 3 ⏳ IN PROGRESS: crowd-aware growth + interior closure (Ryan's "clever setup").**
  Strategy: bias the few near-anchor steps AWAY from each domain's COM so the chain leaves the
  crowded shell cheaply, then finish in open space; keep it a soft SKEW (per-conformer scatter),
  not a fixed stub. Interior = two biased stubs (one per anchor) + a middle segment closed in open
  space. Domain "face-each-other" orientation lives in step 3 / `place.py`.
  - **Keystone DONE + validated:** `grow_free_batch` gained an optional von-Mises azimuth skew
    (`bias_directions`/`bias_residues`/`bias_kappa`) that leans near-anchor steps toward a
    direction WITHOUT touching the exact bond/angle. Unbiased path bit-identical (max diff 0.0);
    +4 tests (38 total); ruff/mypy clean.
    MEASURED (terminal, r=15 Å domain blob): naive self-avoid+clear 19.3/25.8/29.7% (gap 5/10/15)
    → **biased 41.6/42.2/42.5%** — recovers ABOVE the free no-obstacle baseline (39%), because
    leaning outward also cuts early self-clashes. Conformers stay diverse (mean pairwise
    end-direction cosine ~0.28–0.30, i.e. NOT collapsed) — the skew works as intended.
  - **Two-stub + FABRIK closure DONE + measured:** added `close_chain_ends` (batched FABRIK-style
    endpoint closure) and `generate_interior_batch` (two biased stubs + free middle closed onto the
    stub ends). +8 tests (52 total in `test_batch_walk.py`); ruff/mypy clean.
    MEASURED (n=40, gap 42 Å, two r=15 blobs, B=20000): endpoints exactly one bond from each anchor;
    **bond error 5.8e-15, 100% feasible closure**; clears both domains 78%. BUT the FABRIK middle
    ignores the angle window → **only 82% of pseudo-angles in [91,161]**, self-avoid+clears 19%
    (≈ naive 24%). Diagnosis: the stubs bring P,Q close, so the middle is very slack and FABRIK
    CRUMPLES it (bad angles, self-clash). FABRIK is the wrong middle tool.
  - **Shrinking-sphere funnel middle DONE + measured:** added `_bridge_funnel` (in-window cone
    steps with a von-Mises azimuth skewed toward the far anchor, concentration rising as the
    reachability ball shrinks; last residue placed on the two-bond intersection circle) +
    `_sphere_sphere_point`; `generate_interior_batch` now takes `closure="funnel"` (default) or
    `"fabrik"`. 47 tests; ruff/mypy clean. Tuned stub→4 (stub=8 overshot the 42 Å gap and CROSSED)
    and a gentle kappa schedule (kmin 0.5 / kmax 10 / reach 0.90).
    MEASURED (n=40, gap 42 Å, two blobs): **funnel angle-in-window 77%→95.5%** and clears 65→75.5%
    — the shrinking-sphere works. But a genuine TRADEOFF: funnel feasible 73.5% (endgame misses the
    2-bond landing 26% of the time) and clean yield 8.6% vs FABRIK's 100% feasible / 30.2% clean —
    because keeping angles realistic over a very slack span (contour ~148 vs gap 42) forces coiling
    that self-clashes more than FABRIK's kinky-but-avoidant chains. Interior yield is inherently low
    (slack coil, double anchor, two domains) — this is the fallback-to-walk regime.
  - **Excursion drift (Ryan's "skew in a different direction") DONE + measured — big win.** Added a
    per-conformer RANDOM perpendicular excursion drift (faded as the reach tightens) + an optional
    expansion drift to `_bridge_funnel`, exposed as `funnel_excursion` (default 2.5). The excursion
    spreads a slack interior chain so it stops collapsing/self-clashing; expansion (away-from-centroid)
    was measured to KILL closure (radial push fights the funnel) — so excursion only.
    MEASURED (full interior, n=40, gap 42 Å): funnel exc=0 → exc=2.5 takes **CLEAN yield 8.6% → 36.3%**
    (4.2×), Rg 12.6 → 15.1 (toward the free-chain 16.5 — anti-collapse), angle 95.5→95.9%, clears
    75.5→83%. **The funnel now beats FABRIK on every quality metric** (fabrik: 30.2% clean, 77% angle,
    65% clears). Conformers stay unique (excursion is a distribution, not a fixed direction).
    → the interior generator is now genuinely usable (~3× oversample gives enough); realistic angles,
    natural spread, closes, clears domains. 48 tests; ruff/mypy clean.
  - **Excursion tied to predicted Rg DONE + measured.** `generate_interior_batch` now takes
    `target_rg` (e.g. ALBATROSS's Rg). Because the Rg(excursion) curve depends on length, gap and
    stub geometry (and OVERSHOOTS at high excursion for long chains), a fixed map is fragile — so it
    probe-calibrates: grow a 512-conformer probe at 4 excursion values, measure assembled-chain Rg,
    interpolate for the target, clamp to [0, 5]. `_excursion_for_target_rg`. 49 tests; ruff/mypy clean.
    MEASURED: n=80 target 25.0 → achieved 25.8; n=40 target 17.4 → 16.6; n=20 target 12.1 → 9.3
    (clamped — a short pinned chain physically can't spread that far). Achieved Rg is monotone in the
    target and tracks within ~1 A inside the reachable band.
    **Adversarial verification (6-agent workflow) done.** The CALIBRATION passes: Rg tracks target
    within ~1.5 A across a 6.7× length range, strictly monotone (compact target → 4× the self-contacts
    of expanded — a real size change), adapts to gap (excursion 5.0→3.85→0.08 compact→expanded), NO
    overshoot / crash / NaN / false-success. Determinism bit-exact; diversity rich (35 A pairwise RMSD,
    all unique); the probe is UNBIASED (matches a manually-set excursion to 4 sig figs). But it surfaced
    real, mostly PRE-EXISTING issues to fix before pipeline use:
    1. **Endgame closure** (`_sphere_sphere_point`) is the recurring weak point: closes 22–40% at
       short-n / high-excursion, ~74% mid-range. Pushing excursion high for Rg makes it worse.
    2. **Free-chain predicted Rg is ~0.8 A ABOVE the reachable ceiling** for a chain pinned at gap=Re
       (physics: doubly-pinned spreads less than free) → targeting it drives excursion to max. Target
       the pinned-chain Rg, or accept the ceiling.
    3. **Min-excursion floor needed:** at wide gaps the probe zeros excursion (gap alone meets Rg) →
       self-avoidance craters 70%→22%. Keep a spread floor even when Rg is already met.
    4. **Surface diagnostics:** silent clamp + "Rg looks right while chain is broken" (stub-too-large →
       Rg exact, closes 0%). Return achieved-Rg / closure-feasibility / reachable flag so a caller
       can't trust a broken chain. 5. self-avoid collapses at very long n (100→21%, 200→0.1%): fallback
       regime.
  - **FIX #1 DONE (surface diagnostics):** `generate_interior_batch` now returns an
    `InteriorBatchResult` (like `FreeBatchResult`): `coords`, per-conformer `closure_feasible` mask
    (bond within 0.1 A everywhere — the broken-chain filter that closes the "Rg-right-while-broken"
    trap), `radius_of_gyration`, `target_rg`, `rg_clamped`, `excursion`, `.feasible_fraction`. Tests
    updated to `.coords`; +3 tests (51 total); full suite 1871 passed / 0 failed; ruff/mypy clean.
  - **FIX #2 (endgame) — landing-pad trick MEASURED, does NOT help; reverted.** Implemented Ryan's
    idea: grow `landing` residues backward from Q, aim the funnel at the deepest (mid-gap soft
    target); kept as an off-by-default `landing` knob on `_bridge_funnel`/`generate_interior_batch`
    (default 0 = old single sphere-sphere close, bit-identical). MEASURED (n=40): landing 0→3 drops
    clean 35.9%→28.3% and self-avoid 57%→46% while closure barely moves (75.2%→74.2%). WHY: the
    reverse pad lands in the same mid-gap the funnel front occupies and, grown blind to it, they
    self-collide; and it doesn't fix closure because the bottleneck is the funnel's LAST-STEP AIMING
    (cone can't turn fast enough at the very end), not the target distance. Endgame left as-is:
    ~75% closure at normal n is acceptable via over-generate-then-filter + fix-#1's `closure_feasible`
    flag. A possible future attempt: a multi-bond analytic close that RESERVES funnel residues (no
    separate pad → no self-collision). 51 tests; ruff/mypy clean.
  - **FIX #3 DONE + measured (pinned-chain Rg + min-excursion floor).** The full-interior excursion
    tradeoff (measured n=40): clean yield PEAKS at excursion ~3.0 (closes 74/self-avoid 65), poor
    below ~2 (self-avoid <45) and above ~3.5 (closes <65). So `_excursion_for_target_rg` now clamps
    the probe result to the closure-healthy, self-avoiding band **[exc_min=2.0, exc_max=3.5]** — one
    move that is both the min-excursion floor AND "target the pinned-chain Rg" (the band's Rg ceiling
    ~15.9 sits below the unreachable free-chain 17.4). MEASURED wins vs the verification's flags:
    expanded-gap self-avoid **22% → 66%** (clean 15→43%); free-target closure **~40% → 68%** (no
    longer driven to exc 5); normal n=50 Rg 19.4 / closes 71 / clean 43. Rg tracks EXACTLY in-band
    (target 15.0→15.0, not clamped) and clamps honestly outside (free 17.4→15.8, `rg_clamped=True`).
    Tests updated (in-band + below-floor-clamp cases); 52 tests; ruff/mypy clean.
  - **All three post-verification fixes done.** Interior generator now: surfaces diagnostics
    (`InteriorBatchResult` + `closure_feasible`), endgame left as-is (landing pad measured not to
    help; filterable via the flag), Rg driven to the pinned-chain value with a self-avoidance floor.
  - Remaining for Phase 3: Slice 4 (careful `walk.py` fallback for the hard tail + wire the whole
    thing behind a `ConformationEngine` into the pipeline), Slice 5 (failures.txt corpus benchmark).
  - Remaining for Slice 3/4: robustify endgame feasibility (75% → higher via a multi-bond landing);
    then Slice 4 (careful `walk.py` fallback + wire into a `ConformationEngine`), Slice 5 (corpus).
- **Slice 4 ⏳ IN PROGRESS — increment 1 DONE (2026-08-05): the engine adapter + terminal generator.**
  - `generate_terminal_batch` (new, `batch_walk.py`): the missing one-anchor generator — grow from the
    anchor biased away from its domain, filter self-avoid + clears, steer to target Re. Completes the
    free/terminal/interior trio.
  - `BatchWalkEngine` (new module `engines/batch_engine.py`): a `ConformationEngine` (verified via
    `isinstance`) that dispatches an `IDRRequest` by shape (free/terminal/interior), derives `n_away`/
    `c_away` from the obstacle cloud near each anchor and `target_rg` from the sequence, over-generates
    a pool, filters clean (self-avoid + clears + `closure_feasible`), and **falls back to the wrapped
    `SelfAvoidingWalk` for the shortfall** — returning a protocol-correct `IDRResult` (finite iff
    success, NaN-filled failures, total failure raises, engine name records the mix).
  - +11 tests (7 engine, 4 terminal); full suite **1883 passed / 0 failed**; ruff/mypy clean; ADDITIVE
    (nothing in the pipeline imports it yet). NOTE (from Phase 1): the rebuild path calls generate with
    n_conformations=1 per region, so in a rebuild the engine is fast per-region grow-then-filter + walk
    fallback; its many-at-once speedup lands on the ensemble / `build_from_sequence` path.
  - **n=1 speed — MEASURED and DEFERRED by decision (2026-08-05).** At n=1 the engine is currently
    SLOWER than the walk (free 0.9×, terminal 0.2×, interior 0.4×) from fixed overheads: double
    over-generation (engine pool 16 × generator's own oversample 16 = 256 grown for 1), the Rg probe
    (1,536 grows ≈ 9 ms of interior's 19 ms), and double filtering. These AMORTIZE at ensemble scale
    (1000 conformers 21 ms batch vs 12,649 ms walk ≈ 600×). Ryan's call: leave it — sub-second at n=1
    is the bar and it's met; the dream is fast at N=1 AND N>1, satisfied. Fixes (oversample=1 from the
    engine, skip/cache the probe for small pools, filter once) are noted for later, not now.
  - **Increment 2 DONE (2026-08-06): wired into the pipeline + end-to-end verified.**
    - `_make_engine` (`pipeline.py`) grows a third branch: `engine == "batch"` → lazily imports and
      returns `BatchWalkEngine()`; unknown names still raise `InvalidParameterError`. `cli.py`
      `_ENGINE_CHOICES = ("walk", "batch")`. This two-branch + one-choice edit is the ONLY pipeline
      touch — everything else stays additive.
    - End-to-end `rebuild(dnmt3a, engine="batch", seed=0)`: ok=True, builds all 3 regions, **no
      blocking failures, no NEW bond violations** (same 3 inherited from the source as the walk).
      ~6.1 s vs the walk's 0.1 s — the deferred n=1 overhead above, accepted (sub-second bar is about
      per-conformer cost; this run is dominated by the one-time ALBATROSS Rg predictor + double pool).
    - **The anchor-pseudo-angle gap was measured, not a problem.** Rebuilt-region angles on dnmt3a:
      **batch 99.7 % in [91,161]° vs walk 99.1 %** — the stub grows in-window cone steps and the
      away-bias keeps the anchor junction sensible, so batch junctions are as good as (slightly better
      than) the walk's WITHOUT an explicit n_anchor_prev constraint. Concern retired.
    - +2 committed regression tests (`test_pipeline.py::TestBatchEngine`): rebuild-with-batch builds
      every region with exact bonds; unknown engine refused. Full touched suite green; ruff/mypy clean.
    - Not run: the failures.txt corpus needs network/AFDB (offline here) → folded into Slice 5.
- Slice 5: benchmark survival/closure feasibility against the failures.txt corpus (don't assume; needs
  network to fetch the AFDB structures).

### Walk.py hot-loop vectorization  ✅ DONE (2026-08-06)
Killed both per-conformer inner loops in `SelfAvoidingWalk`, acceptance bar = **bit-identical output**
(not just close) + measured speedup. Both met.
- **Loop 1 — cone candidates (`_candidates_for`).** Was one `cone_candidates` call per live conformer
  per residue. New `cone_candidates_batch` (sampling.py) rotates the *shared cached template* onto every
  conformer's axis in one pass via new `rotation_between_vectors_batch` (transforms.py). No RNG on this
  path, so it cannot perturb draw order. Verified bit-identical to the scalar per-row call (max|diff| 0.0
  over 2,000 rotations incl. ±z antiparallel edge cases + 1,500 cone builds); scipy's batched
  `from_rotvec` matches scalar to the bit.
- **Loop 2 — chain-clash KD-tree (`_nearest_chain_distance` → `_chain_clear_mask`).** Was a `cKDTree`
  built + queried per live conformer per residue — the profiler's #1 cost (~65 % of `generate`). Its
  output fed exactly one consumer: `>= CA_CLASH_DISTANCE`. Replaced with a vectorized cull: every
  candidate sits one bond length from its apex, so any chain point beyond `bond + clash` (7.01 Å) of the
  apex cannot clash with *any* candidate → prune to the near points (triangle inequality proves the
  offending point is always kept), gather to a dense block, one brute-force min. `cKDTree` distances are
  NOT bit-identical to `sqrt(Σsq)` (~3.5e-15), but the RNG draw is drawn unconditionally in `_select`
  (independent of this mask), so draw order is provably unchanged and only a distance landing within
  ~4e-15 of exactly 3.2 Å could flip the mask — measure-zero, and checked: 0 mismatches over 853 k
  comparisons with points planted at 3.2 ± 1e-9 Å.
- **Bit-identical end to end:** 42 walk before/after cases (free/terminal/interior-closure/obstacles ×6
  seeds) all max|diff| 0.0; full `rebuild(dnmt3a, engine="walk", seed=0)` max|diff| 0.0.
- **Speed (batch fills the pool):** free n=150 ×30 **3.1×** (549→178 ms), free n=60 ×30 **3.2×**,
  interior n=40 ×30 **2.9×**, free n=150 ×1 **1.5×** (even n=1 wins — loop 2's per-residue tree is gone).
- +10 committed regression tests (rotation/cone batch = scalar bit-for-bit incl. edges + cache sharing;
  `_chain_clear_mask` = KD-tree mask incl. the clash boundary). Full unit suite **1600 passed / 1 skip**;
  ruff + mypy clean. Files touched: `geometry/transforms.py`, `geometry/sampling.py`, `engines/walk.py`
  (+ tests). Additive apart from the two in-place loop swaps.

## Phase 4 — Honesty + CI gate  ✅ DONE (2026-08-06), one item deferred by decision

### README/docs honesty  ✅ DONE (2026-08-06)
Every claim verified against the code first (6-agent workflow + my own re-measurement — do not bake in
an agent's number). All five flagged overstatements fixed, plus three clear contradictions the audit
surfaced. Each fix quotes the code that grounds it.
- **"by default ALBATROSS"** (README ~L15): OVERSTATED. `predict_end_to_end`/`target_dimensions` default
  `prefer_albatross=True`, but with no sparrow (`from sparrow import Protein` → ImportError,
  dimensions.py:156-166) they warn and `return flory_end_to_end(n)` (6.22·N^0.522, sequence-blind). So
  ALBATROSS is the default *only if sparrow is installed*. Reworded to "with ALBATROSS when sparrow is
  installed, and an analytical scaling law otherwise".
- **"~10⁻¹³ Å asserted"** (README ~L185): MISLEADING word "asserted". The enforced tolerance is
  `verify_rigid(..., tolerance=1e-6)` (place.py:127), run once per repositioned domain (not "every
  transform"). 10⁻¹³ is the *observed residual* — I re-measured on real rebuilds: dnmt3a 1.4e-13, arf19
  0.9e-13, p300 max 2.5e-13, so 10⁻¹³ is right as the residual (the audit agent's 1.95e-12 was from
  random rotations, not real transforms). Reworded: asserts to within **1e-6 Å**, residual ~10⁻¹³ Å.
- **"genuine/real ensemble"** (README ~L31, L319, L371): mechanism is TRUE (sample_end_to_end_targets
  draws each conformer's target from the Maxwell radial dist, CV 0.42; walk.py:574-658), but "genuine/
  real ensemble" reads as thermodynamic and it is a geometric sampler. Softened to "a spread of distinct
  conformers", kept the "geometric sampler, not a force field / not a thermodynamic ensemble" caveat
  adjacent. The one unbacked number — "Measured CV is now 0.36–0.53" — replaced with what is backed
  (target CV 0.42; suite asserts achieved spread in [0.20,0.60], test_engines_walk.py:1209-1210). The
  "0.35–0.48" reference is real — it is the freely-rotating-chain CV the code docstring cites
  (walk.py:607-608) — so kept and attributed precisely.
- **region-ID "beats disorder predictors"** (README ~L271,L278 + docs/algorithm.md:74,89): UNSUPPORTED.
  No benchmark anywhere in the repo; metapredict removed (TestMetapredictIsGone); every code docstring
  already hedges it as "the author reports" (constants.py:188, identify.py, contact.py). README/docs
  stated it flatly. Reworded to attribute-and-hedge ("the author reports … not benchmarked here") and
  dropped the flat "beats disorder-based calls".
- **PyPI/install** (README ~L37-63): FALSE — `idptools-dodo` is 404 on PyPI (verified: unrelated `dodo`
  and `numpy` both 200), and the section gave two contradictory install commands. Removed the
  `pip install idptools-dodo` "whole story" block; kept the real git install with a note that PyPI is
  future. **Test count** "1400 tests" → real count is 2011 (`pytest --collect-only`), reworded "~2000".
  **Base deps** listed 4, pyproject has 5 — added the missing `tqdm`.
- Verified-accurate, left unchanged: dnmt3a "4,636 atoms / 334 CA" (reproduced exactly), the
  117-structure corpus (real: tests/data/corpus.json), Python 3.10 floor, cache-size numbers.

### Docs CI job was red — FIXED (2026-08-06, found while validating the doc edits)
Building `sphinx -b html -W docs` (exactly the CI `docs` job) failed with 20 warnings-as-errors — all
pre-existing, none from my edits. Cause: Phase 0 deleted `construct.backbone`, `engines.starling`,
`engines.hierarchical`, `geometry.regularize`, but `docs/api/index.rst` still autosummary-referenced
all four (verified DEAD by import). Removed exactly those four lines (the other 25 modules import
fine). Clean rebuild from a wiped `docs/api/generated/` (gitignored artifacts) → **build succeeded, 0
warnings**. So the docs gate is green again. Unrelated to the historical-failures work below.

### CI offline-gate  ⏸️ DEFERRED by decision (2026-08-06)
Ryan's call: skip the historical-failures offline gate for now — no AFDB downloads, no fixture commits.
The 9 invariants keep running in the `network` job (main/schedule, non-gating), exactly as today. To be
picked up as a separate follow-up when the cross-platform shakeout can be watched. The `continue-on-error`
line stays (it is correct as-is). Full analysis of what the gate would take is preserved below.

Analysis retained for the follow-up:
The plan line "turn off `continue-on-error`" is **stale/wrong**: continue-on-error sits only on the
`network` job (CI.yaml:237), deliberately, so upstream AFDB/RCSB/UniProt outages don't redden main —
and that job still runs genuinely-live tests (TestLive, TestCorpusRebuild). Do NOT flip it. The real
gap: the 9 historical-failure rebuild INVARIANTS are `@network @slow`, so `-m "not network"` skips them
and they never gate PRs. Closing that (workflow CI analysis) needs, and is blocked on:
  1. **Fresh AFDB v6 downloads for all 9** (~3.3 MB raw / 0.74 MB gz). The committed `p300.pdb` is a
     pre-v6 release with *different coordinates* (verified live) — cannot be reused; `CLASH_CEILING=0`
     is a ratchet measured on the v6 bytes. Downloading + committing needs the user's OK.
  2. **Re-validate the 0-clash / 0-impossible-separation ratchet** against the committed bytes before
     gating on them.
  3. **Cross-platform determinism risk**: the invariants use absolute thresholds and `atol=0`; moving
     them off ubuntu-only into the macOS/Windows matrix could flake on BLAS differences. Safer variant:
     gate offline but ubuntu-only.
  Mechanism (once unblocked): commit the 9 v6 PDBs under tests/data/structures/afdb/, make `_fetched`
  resolve to the committed file offline (fetch only when absent), drop `@network @slow` from
  TestHistoricalFailures, share one rebuild across the class's methods to bound runtime.
