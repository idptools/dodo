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

## Phase 1 — Wire up batching that already exists  ⬜ NOT STARTED
- Pipeline: one batched `generate(n_conformations=n_models)` call instead of the external model loop with `n_conformations=1`.
- Backbone: gather all models' CA traces → one batched refinement pass (authors measured 16.4×, left unwired).

## Phase 2 — Vectorize deterministic backbone placement  ⬜ NOT STARTED
- Rewrite `ca_backbone.backbone_from_ca` (per-residue loop) as batched numpy over `(batch × residues)`; keep numba refine kernel, add batch axis. Guarded by existing RMSD-ceiling tests.

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
