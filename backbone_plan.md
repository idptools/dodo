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
  - ✅ **Follow-up RESOLVED (2026-08-07) — and the "3-way" premise was refuted.** Measured the
    coupling graph of the residual clashes: there is **no three-unit cluster** to move. p300's two
    residual coupled components are exactly two units each ({925,933}, {2190,2195}); the rest are a
    movable unit against a *fixed* CA (dnmt3a res386 vs a folded PRO233 CA; p300 res2191 vs a fixed
    ALA2197 CA), which no azimuth can move. The 2-unit joint search already forms and scores those
    groups — what the borderline ones needed was azimuth **resolution**. Raising the grid from
    15 deg (`grid=25`) to 5 deg (`grid=73`) clears one more p300 contact (introduced **4 → 3**) with
    the N-CA-C angle distribution byte-identical and ~+1 s; ratchet tightened to **p300 (3, 10)**.
    Rejected: widening `angle_window` to (70,170) clears one further contact but only by placing a
    165.7 deg N-CA-C angle (real ~111 deg) — trading a 0.05 A borderline clash for a distorted
    backbone angle. The two remaining p300 contacts sit against fixed CAs (unclearable by the 1-DOF
    azimuth model); genuinely clearing those is CA-engine territory, not a backbone-polish job.

## Phase BB-2 — Solve the seam (the first-class blocker)  ✅ DONE (2026-08-07) — solution is C

**Outcome: the seam cannot be forced closed without violating a hard invariant. The principled
solution is C (honest validator reframe), already committed. D (seam-aware closure) was implemented,
measured to close zero additional seams, and fully reverted.** Decision confirmed with Ryan: finalize
C, stop trying to force closure, move to BB-3/4/5.

- ✅ **C — validator reframe DONE (2026-08-07).** A long C-N across a generated↔input boundary is no
  longer a `chain_break` blamed on `"rebuilt"`; it is `kind="seam"`/`provenance="seam"` (validate/
  bonds.py) — reported, not counted against `ok`, not attributed to the rebuild (its strain is
  inherited from the rigidly-repositioned folded N). Result: `--backbone` introduces **0
  rebuilt-provenance bond violations** on the corpus (seams: dnmt3a 4, arf19 6, p300 10, all now
  "seam"; arf19 reaches `ok=True`, dnmt3a/p300 `ok=False` only from inherited *input* defects — the
  same differential-clean bar the CA path meets). +safety test: the identical geometry read back
  without region info is a `chain_break`, so a real break in non-DODO input is never masked. Zero
  engine risk; full suite 1548 passed.
- ❌ **D — seam-aware closure: IMPLEMENTED, MEASURED USELESS, REVERTED (2026-08-07).** Threaded the
  folded N through `IDRRequest.c_anchor_n_xyz` → `_WalkPlan.end_seam_target` and biased `_close`'s
  candidate selection toward it under the existing hard angle filter. **Closes 0 additional seams.**
  The earlier "73% reachable" estimate was measured over the full closure circle (`circle_min`) at a
  single seed; the *definitive* measurement (corpus × 3 seeds, 83 C-side closures, accounting for the
  angle filter) shows why single-step biasing can't work:
  - **CLOSE (22/83):** the near-N arc is *also* angle-valid → these already close in the seam-blind
    walk. Biasing them is a no-op.
  - **FAR (40/83):** the whole arrived circle is >2.854 Å from N (`circle_min > 2.854`). The last CA
    sits on a circle fixed by where the walk *arrived* (CA(n-2)); no last-step choice reaches N.
  - **ANGLE (21/83):** the circle reaches N but the [91,161]° pseudo-angle window forbids the near-N
    arc (`circle_min ≤ 2.854 < valid_min`, `n_valid_reach=0`). Using it would kink the invariant.
  So 59/83 are unfixable by closure-step selection and the fixable 22 already close. Empirically: a
  soft-Gaussian weight left rebuilt CA coords **bit-identical** (the racing draw `argmin(-log U/eff)`
  swamps modest weight ratios — near candidate wins only `1/(1+Σ eff_far)`); a hard reach-threshold
  weight moved one p300 CA but closed 0 seams and made one 0.026 Å *worse*. Reverted to the committed
  C-only state (`git diff` empty, 241 tests green, mypy clean).

**Why no bounded fix exists:** the folded N points where AlphaFold's *original* chain ran; step-3
rigid repositioning places the domain for the region's anchors/end-to-end, blind to its own boundary
N. Closing the seam short needs either kinking the CA-CA-CA pseudo-angle (approach B / CCD: 100% but
41-79° vs the 91° floor) or moving a folded domain non-rigidly — both hard-invariant violations. The
strained seam is an honest, irreducible artifact of independent rigid-repositioning + chain rebuild;
C surfaces it truthfully. Seam counts are frozen ratchets in `_BACKBONE_BASELINE` (down-only).

Evidence trail for the rejected/earlier roads:

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
- **Resolution (2026-08-07):** BB-2 is closed. The measured conclusion is that no closure-side fix
  can force the seam short without kinking the pseudo-angle (approach B/D) or moving a folded domain
  non-rigidly (domain repositioning) — both hard-invariant violations. C is the principled answer and
  is committed. If a *complete* seam fix is ever revisited, the only remaining route is seam-aware
  step-3 domain placement (orient the folded domain so its boundary N/C faces the IDR's approach),
  which fights the end-to-end target and is out of scope for v2.0.
- Effort spent: **L**; outcome: C (low-risk, committed), D (reverted).

## Phase BB-3 — Accuracy: the 2D peptide-plane table  ⏸️ BUILT & MEASURED, HELD (2026-08-07)
The one measured accuracy win over the shipped table — implemented in full, then held because the
win is too modest to justify its cost. **All BB-3 code reverted; the 1D table stands.**

- **The win is real but modest and could not be reproduced at the plan's headline size.** Built the
  2D table indexed on (τ_i, τ_{i+1}) — a 5-CA predictor — deriving it from the committed 19,302
  units, baking it in, and wiring placement (2D where a fifth CA exists, 1D row for the last
  interior unit and for backfilled cells). Independent 5-fold CV on **placed-atom** error (the
  number that matters): **C 0.3139→0.2978 Å (−5.1 %), N 0.1941→0.1868 Å (−3.8 %), paired 95 % CI
  excludes 0 for both.** The direction/significance match the plan, but the plan's absolute figures
  (0.2286→0.2168) are a different projection; the honest placed-atom gain is ~0.016 Å on C. A
  min-count threshold on sparse 2D cells only *hurt* held-out accuracy (monotonically), so the full
  table is the best 2D variant — the sparse cells are signal, not overfit.
- **Blocker: it trips the p300 clash ratchet 4→5.** The (more accurate) 2D placement shifts the
  2191-2197 backbone cluster so one more *borderline* contact crosses the vdW cutoff — the new one
  is `GLN2192 CA / ALA2197 N` at **2.648 Å vs a 2.65 Å minimum, 0.002 Å over**. All five introduced
  clashes are 0.002-0.086 Å borderline; four are in the same coupled cluster that is the documented
  **3-way-clash-move follow-up** (the 2-unit polish can't escape a 3+-body cluster). Noise-level, but
  a regression by the down-only gate's definition.
- **Decision (Ryan, 2026-08-07): HOLD.** The ~0.016 Å placed-C win does not justify (a) a 324-entry
  (~292-line) table during the v2 de-bloat refactor [[dodo-v2-refactor-not-restart]], nor (b) a
  clash-ratchet regression that would itself require the separate 3-way-polish work to avoid —
  especially as the real path to backbone/all-atom quality is Ryan's external CA→all-atom model
  [[dodo-external-all-atom-model]], which makes DODO's analytical backbone a stopgap. If revisited,
  the clean order is: 3-way clash polish first (clears the cluster, ratchet stays at 4), then BB-3.
- **Do NOT** pursue finer bins, bin interpolation, or a separate O predictor — all measured to
  give nothing (the 1-dihedral residual is the plane's intrinsic spread; O is geometrically
  determined by C/N). The derivation script + `--verify`/CV harness for the 2D table are in this
  session's scratchpad if the work is picked up again.
- Effort spent: **M**; outcome: reverted, 1D table retained.

## Phase BB-4 — Speed  ✅ DONE (2026-08-13) — but the plan's target was wrong
**The plan said "vectorize the refinement (prange over models)." Profiling refuted that:** the
compiled refinement kernel (`sweep_region`) is only **~2.5 %** of `--backbone` wall time — already
fast. The real cost was `_polish_coupled_clashes`, the coupled-clash polish (BB-1/BB-2.5), which was
**pure Python and 78 %** of `--backbone` time (p300 ×2 models: 6.16 s of 7.92 s), amplified further
by the BB-2.5 grid 25→73 change (8.5× more candidates per 2-unit group).

- ✅ **Vectorized the clash polish (`ca_backbone._polish_coupled_clashes`).** Two changes, both
  behaviour-preserving:
  1. **Precompute each moved atom's static neighbours once per group** instead of a
     `cKDTree.query_ball_point` per atom per candidate. The static cloud is fixed during a group's
     search and each moved atom rides a bounded circle, so a `query_r + reach` shell is a provable
     superset (far neighbours add zero headroom) — exact same overlap sum.
  2. **Batch the deepest unit's whole azimuth grid** (`_place_units_batch` / `_place_oxygen_batch` /
     `_angles_batch`) with a `const + var` split, so `argmin` over the valid candidates is the same
     global minimum the scalar loop found — replacing the per-candidate Python placement/angle/score.
- **Result: polish 6.16 s → 0.40 s (15×); p300 ×4 `--backbone` 13.8 s → 4.59 s (3×)**; single-model
  p300 2.5 → 1.4 s, dnmt3a 1.1 → 0.27 s. **Output is bit-identical** to the scalar polish (worst
  |Δ| = 0.0 over corpus × 4 seeds), clash ratchet unchanged (2/0/3), determinism preserved.
- **prange-over-models NOT wired** — with the polish at 0.40 s the refinement kernel is no longer the
  bottleneck, so parallelising it would chase ~2.5 %. Revisit only if profiling shifts.
- Effort spent: **M**; risk turned out low (pure-numpy, bit-identical, no numba parallel kernel).

## Phase BB-5 — Promote to first-class  ✅ DONE (2026-08-13)
The release is backbone-only and aims to be first-in-class, so `--backbone` is now the default and
is held to corpus-wide invariants. **Promoting it surfaced three real defects that opt-in status had
been hiding** — the corpus harnesses had only ever run the CA-only path.

- ✅ **Default flipped.** `backbone=True` on `rebuild()` and `build_from_sequence()`;
  `[--backbone | --no-backbone]` on the CLI (argparse `BooleanOptionalAction`). Rebuilt alpha
  carbons are bit-identical either way, so the backbone stays purely additive. 19 test call sites
  that used a plain `rebuild()` as a *CA-only* reference now pass `backbone=False` explicitly.
- ✅ **New offline regression gate:** `TestBackboneIsFirstClass` (tests/unit/test_pipeline.py) over
  every committed fixture × seeds 0-2 — the guarantees that must hold at *every* seed, where
  `TestBackboneBaseline` pins only seed-0 ceilings: nothing impossible introduced, zero
  rebuilt-provenance bond defects, **complete N/CA/C/O on every rebuilt residue with N-CA / CA-C /
  C-O exact to 1e-6**, seam count matching the report, and byte-identical output on re-run. Plus
  `test_clashes_within_ceiling` in the historical suite now asserts BOTH regimes (CA-only holds its
  absolute 0; backbone gets a measured `BACKBONE_CLASH_CEILING = 5`).
- ✅ **Three defects found and fixed** (all pre-existing, all corpus-only):
  1. **Impossible contacts at seams.** The exact two-sphere seam placement satisfies both bond
     lengths and was blind to everything else — PTBP2 put a rebuilt N 0.797 Å from a folded carbonyl
     O; Q15642's exact C induced an O onto a folded atom **56 residues away in sequence** (invisible
     to the narrow neighbouring-residue obstacle set). The exact point is now *checked before it is
     accepted* — against a 4.0 Å shell of every structure atom, against the carbonyl O it induces
     (O is determined by C, so only a different C can move it), and against the residue's own
     N-CA-C angle — and falls back to the obstacle-avoiding cone, which now scores the induced O too.
  2. **Refinement collapsing N-CA-C under clash pressure.** Q8N8A8 and O60271 reached ~79°, putting
     a residue's own N and C 1.90 Å apart — exactly what the bond validator flags as two atoms on
     top of each other. The soft angle term (weight 0.124) is deliberately weak and loses to the
     clash term (40). Both backends now apply a hard 80-160° window
     (`constants.N_CA_C_WINDOW_MIN/MAX`, derived from the validator's 1.90 Å floor) as a step
     penalty; a test pins the kernel's literal copies to the shared constants.
  3. **A crowded region failing the rebuild outright.** `MAX_NEIGHBOURS = 48` overflowed on Q9C000
     (54 neighbours) and raised `GeometryError`. Raised to 96 (~1.8× the new worst), and
     `backend="auto"` now degrades a still-over-cap region to the uncapped numpy path instead of
     failing; explicit `backend="numba"` still raises.
- ✅ **Docs de-hedged and corrected** across README, `docs/guide.md`, `docs/algorithm.md`, CHANGELOG,
  the CLI help and both API docstrings: no more "opt-in for now" / "while it earns confidence" /
  "what would make this the default" / backbone listed under "what comes next". The stale
  Known-limitations clash figures (dnmt3a 0→1-3, arf19 0→1-15, p300 2→13-31) are replaced with the
  post-BB-1/BB-2.5 measurements (2 / 0 / 3), and the seam entry now records both rejected closure
  approaches with their evidence.
- **Result: 1914 passed, 117 skipped, 0 failed** — including the full 117-structure corpus and
  9-structure historical suites running with the backbone ON. `--backbone` introduces **zero**
  separations below the 1.00 Å floor across all 117 structures.
- Effort spent: **M**; risk: the default flip is a deliberate breaking change to default output.

---

## Suggested order & rationale
BB-0 (foundation ✅) → **BB-1** (clashes down ✅) → **BB-2** (seam: solved by C, honest labeling;
closure-forcing proven infeasible ✅) → **BB-3** (accuracy: built, measured, HELD ⏸️ — modest win not
worth the table bloat + clash-ratchet friction) → **BB-4** (speed ✅ — vectorised the clash polish,
the real 78 % bottleneck, 3× on multi-model `--backbone`, bit-identical) → **BB-5** (promoted to
first-class ✅ — default-on, corpus-wide gate, three latent defects fixed). **All phases are now
resolved.** The release is backbone-only and first-in-class (all-atom deferred), so the open
question worth revisiting is **BB-3**: with the clash ratchet now defended corpus-wide and the
polish 15× faster, the 2D table's accuracy win may be worth its cost after all — re-measure before
deciding.
