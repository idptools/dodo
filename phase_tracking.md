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

**Deferred (the "backbone batch", = old Phase 1 item 2):** `_backbone_atoms` is batch-ready, but
`add_backbone_to_rebuilt` still calls it per region/model. Wiring the `(B,n,3)` path through it —
and adding a batch axis to the numba **refinement** kernel (the measured 16.4×, still unchanged) —
is the follow-up, gated on restructuring the per-region seam/refine loop.

## Phase 3 — Hybrid vectorized IDR generator  ⬜ NOT STARTED
- PRIMARY: batched grow-freely-then-filter across regions & conformers (terminal + connecting IDRs); local clash check only.
- Closure: batched SHAKE-with-endpoints for interior IDRs.
- FALLBACK: current `walk.py` reject-during-growth, per region, only when a batch fully fails.
- Loops stay on the careful path (most constrained).
- Also kill the two per-conformer inner loops in the walk.
- Measure survival rate / closure feasibility against the failures.txt corpus (don't assume).

## Phase 4 — Honesty + CI gate  ⬜ NOT STARTED
- Fix the five README overstatements (region-ID claim, ~10⁻¹³Å→1e-6Å, "physical/genuine ensemble", "by default ALBATROSS", PyPI state).
- Make the 9 historical failures gate CI offline from cached fixtures; turn off `continue-on-error`.
