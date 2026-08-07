# DODO backbone (N, C, O) redo — phased plan

Living plan for making DODO's **analytic** CA→backbone path first-in-class, in the same
spirit as the CA-engine work: measure → correctness → accuracy → vectorize → promote.

**Scope decision (2026-08-06, Ryan):** invest in the analytic path itself as a long-term
first-class method. The external CA→all-atom model is **explicitly out of scope** here — do
not plan around it, do not add a pluggable backend for it. Side chains remain out of scope.

**Verification bar (every phase):** `pytest` green in `hh1`; `ruff` + `mypy --strict` clean;
the CA-only path stays bit-identical (backbone work must never perturb alpha carbons — see the
core invariant below); each claimed improvement is *measured* on the corpus, not asserted.

**Core invariant that constrains everything below:** alpha carbons are DODO's scientific
claim and are held fixed; folded-domain atoms are never moved except as a rigid body. Any
backbone fix that would move a CA or a folded-domain atom trades the thing DODO is trusted for
and must be justified explicitly.

---

## Measured baseline (2026-08-06, grounding workflow; seeds 0/7, dnmt3a·arf19·p300)

The CA-only path is valid (0 rebuilt-provenance bond violations, 0 clashes). `--backbone` is
what introduces every defect below (CA coords are bit-identical backbone on/off, verified).

| metric | dnmt3a | arf19 | p300 | notes |
|---|---|---|---|---|
| introduced steric clashes | 0→9 | 0→1 | 2→12 (10 new) | all 2.1–2.8 Å, **0 impossible (<1 Å)** |
| strained/long seam peptide bonds | 4/5 | 6/6 | 10/10 | 20–21 total, 2.6–4.0 Å vs ideal 1.329 |
| rebuilt-region internal bonds | exact | exact | exact | N–CA/CA–C ~1e-14 Å; interior C–N 4.3e-7 Å |
| wall-clock added by `--backbone` | +0.072 s | +0.182 s | +0.317 s | refinement is ~90 % of it (refine=False adds ~0.01 s) |

**Two things make the output validate as INVALID**, and they are the whole job:

1. **Introduced clashes** — corpus total splits **9 interior · 5 rebuilt↔folded · 6
   inter-region**; mostly recoverable (a refinement-quality problem).
2. **Seam chain-breaks** (21/21, every one geometrically unsatisfiable) — the deep blocker.

**Seam root cause (measured):** the folded domain is moved rigidly in step 3 to hit the IDR's
predicted CA end-to-end distance, while the walk redraws the IDR CA trace *without ever seeing
the folded N/C*. So the fixed neighbour N/C still points where AlphaFold ran the region — 3.0–
5.2 Å from the rebuilt boundary CA (mean 3.85), past the **2.854 Å** (=1.525+1.329) a peptide
unit can bridge. The written seam C–N bonds land at 2.6–4.0 Å → all flagged `chain_break`.

**Accuracy is not reproducible in-repo:** the 4-CA peptide-plane table is a baked-in 18-row
literal (`ca_backbone.py:101-155`); its raw provenance (`subset_frames`, 100 IDR frames) and
the derivation script live *outside* the repo. The cited held-out errors (N 0.164 / C 0.210 /
O 0.614 Å) are one-off prose, not reproducible from anything committed.

---

## Phase BB-0 — Reproducible baseline & honest accounting  ✅ DONE (2026-08-06)
Foundation: you cannot claim "first-class" against numbers you cannot reproduce, and the
current output lies by omission.
- ✅ **Seam accounting surfaced on `RebuildReport`.** The dead `seams_strained` local is now a
  returned value: `add_backbone_to_rebuilt` returns `BackboneResult(structure, strained_seams)`;
  `RebuildReport.backbone_seams` carries `SeamStrain(residue, side, bond_length)`, one per model
  per seam; `summary()` reports the distinct count and worst bond. (Introduced *clashes* are left
  to the regression test below rather than run in the hot path — the validator isn't cheap and
  the report already carries the seam story; revisit if a per-run clash count is wanted.)
- ✅ **Baseline frozen as a regression floor.** `test_pipeline.py::TestBackboneBaseline` ratchets
  introduced clashes ≤ {dnmt3a 9, arf19 1, p300 10} and strained seams ≤ {4, 6, 10} (seed 0,
  committed fixtures), asserts 0 impossible contacts (hard invariant), that interior rebuilt
  bonds stay exact (every rebuilt-provenance violation is a seam `chain_break`), and that the
  report's seam count matches the validator's independent chain_break count. p300 leg is `@slow`.
- ✅ **Doc/code drift fixed.** README & guide described the seam as "~2.2 Å from aiming C at the
  N"; the code holds the residue's N-CA-C angle (`_place_on_cone`) and the bonds measure **2.6–3.7
  Å, mean ~3.0** (reproduced on the committed corpus). Also fixed the guide's internal 17-vs-20
  seam-count inconsistency. Docs build clean under `-W`.
- ✅ **Reproducible derivation pipeline committed** (Ryan's call: commit the pipeline). The 100 IDR
  frames, stripped to backbone atoms and gzipped (**1.4 MB**), are committed to
  `tests/data/backbone/frames/`; `scripts/derive_peptide_table.py` re-derives every shipped
  constant (`_C_BY_BIN`/`_N_BY_BIN` + both marginals) **exactly** (max dev 5e-5 Å, 19,302 units);
  `tests/unit/test_backbone_table.py` runs `--verify` + pins the placement accuracy in CI; README &
  guide now point at the harness. The exact frame/dihedral/bin convention was reverse-engineered and
  proven bit-for-bit against the committed table. This unblocks BB-3 (the 2D-table re-derivation).
- 🐞 **The harness immediately caught a real bug** (now a BB-1 item, see below): the forward
  marginal is *derived* with `reference = CA[i+1]-CA[i+2]` but *applied* by `_backbone_atoms` with
  the opposite sign, so the **first peptide unit of every region is mis-placed — C error 1.06 Å vs
  the 0.31 Å the marginal actually delivers (~3×)**. Recorded as an xfail in `test_backbone_table.py`.
- Effort: **S–M**, risk: **low**. Verification: ruff/mypy clean, docs `-W` clean, 113 tests pass +
  1 documented xfail; CA-only path untouched (never calls the backbone pass).

## Phase BB-1 — Correctness: first-unit sign bug + introduced clashes  ⏳ IN PROGRESS (2026-08-06)
- ✅ **Forward-marginal sign bug FIXED (2026-08-06).** `_backbone_atoms` set `reference[..., 0, :] =
  ca[2] - ca[1]` for the first unit, but the committed `_C/_N_FORWARD_MARGINAL` were measured with
  `reference = ca[i+1] - ca[i+2]` (proven in BB-0), so the first peptide unit of every region was
  placed in a sign-flipped frame. Flipped the unit-0 reference to `ca[1] - ca[2]` (one line +
  a "do not simplify" comment). Measured: first-unit C error **1.055 → 0.345 Å (~3×)**, now
  delivering the forward marginal's own accuracy. **Baseline unchanged** — corpus clashes/seams
  stay {9/4, 1/6, 10/10}, 0 impossible (the first-unit move creates no clash and doesn't alter seam
  reachability), so the BB-0 ratchets still hold. Confirmed isolated: the fix touches only unit 0;
  interior units and the `_terminal_*` helpers (own references, already correct) are untouched.
  `test_backbone_table.py`'s xfail is now a passing regression guard; the derivation script's note
  updated. **Full suite: 1546 passed, 2 skipped, 0 failed**, ruff/mypy clean. One side-effect worth
  recording: the fix *removed a self-clash* on the `_truth()` fixture that the buggy first unit was
  creating, which had been the only thing making `test_backbone_kernel`'s no-obstacle clash-term
  equivalence case non-zero — that case is now skipped (the corrected chain is self-avoiding; the
  crowded/obstacle case still proves the numpy↔kernel clash equivalence). i.e. a test had been
  quietly depending on the bug.
- ✅ **Introduced clashes: joint 2-unit polish DONE (2026-08-06).** Diagnosed first: every introduced
  clash is already inside the refinement objective's threshold, so it is not a blind-objective
  problem — it is a *coupled local minimum* single-azimuth descent cannot escape (proven: 200
  sweeps, a constant-wide window, and a denser grid all leave dnmt3a at 9→7-8; joint moves reach 1
  in a prototype). The fix is a new structure-level **joint 2-unit azimuth polish**
  (`_polish_coupled_clashes` in ca_backbone.py, wired into `add_backbone_to_rebuilt` after seam
  stitching): for each residual clash it jointly searches the azimuths of the one or two *interior*
  units (both flanking CAs rebuilt, so no folded or seam atom moves) that place the clashing atoms,
  keeping any combination that lowers the validator's own van der Waals overlap without pushing a
  moved residue's N-CA-C angle out of [80,160]°. Deterministic; runs over the assembled structure
  so inter-region clashes are visible.
  - **Result: corpus introduced clashes 20 → 6** — dnmt3a 9→2, arf19 1→0, p300 10→4. Geometry stays
    valid (0 rebuilt bad bonds, 0 impossible), overall N-CA-C distribution unchanged, +~0.1 s.
  - **Residuals are the fixable floor / near-fundamental:** dnmt3a's 2 both involve a *fixed* rebuilt
    CA (CA232) the loop's atoms cannot escape; p300's include 2 inherited AlphaFold clashes, a
    fixed-CA one, and a couple of tight both-movable clusters that would need 3-way moves
    (diminishing returns → a documented follow-up). Fixed-CA cases are CA-engine/BB-2 territory.
  - Tightened the BB-0 ratchet to {dnmt3a 2, arf19 0, p300 4}; +1 test proving the polish earns its
    place (stubbing it out makes the output worse). Full suite **1547 passed, 2 skipped**; ruff/mypy
    clean.
  - The grounding's other refinement ideas (element-vdW objective, clash-aware kernel termination)
    proved **unnecessary**: the objective already sees these clashes, and the joint polish clears
    them without touching the numpy/numba refinement or its equivalence tests.
  - Follow-up (open): 3-way / cluster moves for the few remaining both-movable p300 clusters.

## Phase BB-2 — Solve the seam (the first-class blocker)  ⏳ SCOPED — needs a focused effort (2026-08-07)
The crux of "first-class", and the two easy roads were tried and rejected on evidence:

- **BB-2a — honest omission: TRIED, REVERTED.** Omitting the unsatisfiable seam atom (C+O / N)
  does remove the impossible long bond (chain_break → 0), but it (1) does NOT reach
  `validate_bonds.ok=True` — it trades chain_break for `missing_atom` (an incomplete residue),
  both leave ok=False; and (2) breaks the core invariant *every rebuilt residue has a complete
  N/CA/C/O backbone* that 5+ tests enforce, producing gappy output. Aligned with DODO's "nothing
  impossible written" ethos, but a bigger, more disruptive change than a "safety net" — reverted to
  the reported long-bond (BB-0 already surfaces it). Kept the harmless polish guards for absent atoms.

- **BB-2b — seam-aware CA closing: FEASIBILITY MEASURED, only a partial fix.** Steering the terminal
  rebuilt CA on its closure circle toward the folded neighbour's N reaches within the 2.854 Å a
  peptide unit spans for only **~73% of C-side seams** (measured, seed 0: dnmt3a **1/3**, arf19
  **3/3**, p300 **4/5**). For the rest, the *entire* closure circle — every point 3.81 Å from both
  CA(stop-2) and the anchor CA(stop) — is >2.854 Å from N(stop) (dnmt3a res392 circle-min 3.07,
  res472 3.08; p300 res567 3.76), because the rigidly-repositioned anchor/domain geometry puts the
  fixed N out of reach of the whole circle. So single-CA steering cannot close those; a **complete**
  fix needs one of:
  - multi-CA steering (move CA(stop-2) too, changing the circle) — cascades deeper into the walk; or
  - **seam-aware domain repositioning** (step 3): orient/place the folded domain so its boundary N/C
    faces the IDR's approach, which is the true root (the domain was placed for the CA end-to-end
    target, blind to where its own boundary N points). This is option (f), previously dismissed as
    "fights the end-to-end target" — but it may be the only route to 100%.
- **Design for the partial (73%) BB-2b, if pursued:** add seam-target fields to `IDRRequest`
  (folded boundary N on the C side, C on the N side), populate them from the structure in the
  pipeline, and in the walk's closure prefer the circle point nearest the seam target that still
  holds the CA-CA-CA pseudo-angle in-window; keep it seed-deterministic and gated behind
  measurement (CA-CA exact, angle window, end-to-end on target, no walk-test regression). Residual
  unreachable seams keep the reported long-bond fallback.
- Rejected: (b) relaxing the boundary CA needs a median 0.71 Å / max 2.32 Å CA move → breaks the
  exact-CA invariant; (c) folding the boundary residue in is documented to fail (guide:171-184).
- **Recommendation:** BB-2 is genuinely hard — neither easy road reaches first-class, and the real
  fix touches either the closure or step 3, both risky changes to the just-stabilized CA engine, for
  a ~73%-or-needs-more payoff. It warrants its own careful, measured session rather than being
  crammed in. BB-0 + BB-1 stand on their own and are committable now.
- Effort: **L**, risk: **high**.

## Phase BB-3 — Accuracy: the 2D peptide-plane table  ⬜
The single measured win over the shipped table.
- **2D table indexed on (τ_i, τ_{i+1})** — effectively a 5-CA predictor. Measured (5-fold CV,
  19,302 units): C 0.2286→0.2168 Å, N 0.1565→0.1499 Å; paired improvement CI excludes 0; a
  continuous-kNN ceiling (0.215) shows a plain 2D bin table captures nearly all the signal, so
  no fancy regression is needed. Keep the current 1D row as the fallback for the last interior
  unit (no CA[i+3]). Bonds stay exact by construction; O improves ~5 % by inheritance.
- **Do NOT** pursue finer bins, bin interpolation, or a separate O predictor — all measured to
  give nothing (the 1-dihedral residual is the plane's intrinsic spread; O is geometrically
  determined by C/N).
- **Dependency:** re-deriving the table needs the external sim data + a committed derivation
  script → do this *with* BB-0's reproducibility work, not before it.
- Effort: **M**, risk: **low** (accuracy-only, no invariant touched).

## Phase BB-4 — Vectorize the refinement  ⬜
Bring the backbone path up to the CA path's batching philosophy.
- Wire the deferred **prange-over-models** refinement (measured **20.9×** on dnmt3a's 282-res
  region; ~16.4× end-to-end on multi-model `--backbone`). Restructure `rebuild()` to run
  backbone as a **post-model batched pass**, stacking every model's regions and compiling the
  sweep with `njit(parallel=True)`. Placement is already batched (Phase 2).
- Payoff is opt-in `--backbone` × multi-model × multi-core only — schedule after BB-1/BB-2
  make the output actually worth producing at scale.
- Effort: **L**, risk: **medium** (numba parallel kernel + rebuild restructure).

## Phase BB-5 — Promote to first-class  ⬜
- **Corpus regression tests:** after BB-1/BB-2, assert 0 introduced clashes and 0
  rebuilt-provenance seam chain-breaks across the fixtures (and add to the offline gate if that
  lands).
- **Honest, updated docs;** consider softening `--backbone`'s caveats and, eventually, whether
  it can be trusted enough to change its default posture.
- Effort: **M**, risk: **low**.

---

## Suggested order & rationale
BB-0 (foundation) → **BB-1 + BB-2a** (fastest route to *valid* output: clashes down, seams
omitted-not-lied-about) → **BB-2b** (the real seam fix, the hard/ambitious part) → BB-3
(accuracy) and BB-4 (speed) in parallel once validity holds → BB-5 (promote). BB-2b is the one
item that reaches back into the CA engine; everything else is contained to the backbone modules.
