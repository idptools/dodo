# Multi-chain complexes: implementation plan

Covers ROADMAP item 1: rebuild the IDRs of a multi-protein complex without perturbing the
complex. The relative positions of the folded subunits must survive the rebuild exactly.

Everything in section 0 was measured on this branch (`main`, f3216a5) before the plan was
written. Every number below is reproducible with the scripts named beside it.

**Decisions taken 2026-08-21, before writing any code:**

1. Locking crosses chains only. Two folded domains of the same chain that pack against each other
   are *not* locked by default -- re-sampling that arrangement is what step 3 is for, and the
   single-chain path stays byte-identical. Intra-chain locking is an opt-in flag (2.2, C-5).
2. Protein only for this release. Nucleic acids get a clear refusal naming why; dropped ligands get
   an explicit note. Rigid passengers are a later phase (12).
3. Missing-residue complexes (ROADMAP item 1, case 2) are a separate effort, sketched as C-6 and not
   blocking this one.

---

## 0. What DODO does to a complex today

### 0.1 The headline

A two-chain complex built from two copies of `tests/data/structures/dnmt3a.pdb`, positioned so
the chains' folded domains touch (4.00 A closest heavy-atom approach, 5,540 atom-atom contacts
within 5 A across the interface):

```
dimer rebuild: 2.1 s
report.ok: True
1 model(s), 8/8 region-model pair(s) rebuilt
  4 folded domain(s) repositioned:
    chain A FD 281-433: moved in 1 attempt(s), anchor separation 58.7 A against a 58.7 A target
    chain A FD 473-912: moved in 1 attempt(s), anchor separation 46.4 A against a 46.4 A target
    chain B FD 281-433: ...
    chain B FD 473-912: ...
  15 backbone seam(s) left with a strained peptide bond (worst 141.83 A vs 1.33 A ideal)
```

The interface is gone. Contacting folded-domain pairs: **7 before, 0 after**; the two largest
interfaces go from 5,540 atom contacts to 0. Domain centroids move 62-121 A. `report.ok` is
**True**.

### 0.2 Root cause, from the code

`reposition_folded_domains` (`src/dodo/construct/place.py`) iterates **one chain at a time**, and
`placed: list[Domain] = [folded[0]]` is re-initialised inside that per-chain loop. So the obstacle
set for a placement contains only the current chain. Instrumented on the dimer:

```
obstacle set for one placement:    261 atoms, chains ['A']
obstacle set for one placement:   1485 atoms, chains ['A']
obstacle set for one placement:    261 atoms, chains ['B']
obstacle set for one placement:   1485 atoms, chains ['B']
total atoms in the structure: 14292
```

Two independent defects follow: (a) no cross-chain geometry is preserved, because nothing knows
the chains are related; (b) no cross-chain geometry is *avoided* either, because the clash check
cannot see the other chain.

### 0.3 Two more defects of the same kind, both reachable from the public API

These are not multi-chain bugs. They are the same missing abstraction -- "some folded domains
must move together" -- showing up on a single chain, and the fix in this plan removes them too.

**Adjacent folded domains are torn apart.** `assign_regions_from_spec` (shipped in ca1e057)
accepts two folded domains with no IDR between them. `_linker_between` returns `None`, and the
code leaves the second domain at its input coordinates *after having moved the first*:

```
chain A FD 474-700: moved in 1 attempt(s), anchor separation 47.6 A against a 47.6 A target
chain A FD 701-912: held in place (no IDR between this domain and the previous one, so their
                    relative position is not ours to change)

peptide bond C(700)-N(701): 1.334 A before -> 93.479 A after
```

The reason string states the opposite of what happened.

**A region that will not be rebuilt is still stretched.** `rebuild(min_length=...)` skips regions
shorter than the minimum -- they keep their input coordinates -- but step 3 still moves the
domains on either side to the linker's predicted separation:

```
rebuild("dnmt3a.pdb", min_length=50, seed=0)
  chain A FD 474-912: moved in 1 attempt(s), anchor separation 47.6 A against a 47.6 A target
  model 1 chain A residues 433-473: NOT BUILT (shorter than the 50-residue minimum; left as-is)

peptide bond C(473)-N(474): 1.342 A input -> 33.605 A output
```

`dodo validate` catches this (`the chain is broken there rather than merely strained`). The
rebuild report does not mention it.

### 0.4 What already works, and must not regress

* **Region identification is already complex-aware.** `contact_profile` scores burial over the
  whole structure, so a residue buried by a partner chain reads as folded. Measured: monomer
  dnmt3a gives `folded 282-432, folded 473-912`; the dimer gives an *extra* folded domain at
  192-222 -- the helix that packs against the partner. That is the behaviour we want (a bound
  segment is not a free IDR), and section 5 depends on it.
* **IDR building is already complex-aware.** `_obstacles_for_span` uses
  `Structure.placed_atom_mask()`, which spans every chain, and `_rebuild_one_model` marks every
  folded domain `placed` before building anything. Rebuilt IDRs do avoid other chains.
* **A complex with one folded domain per chain rebuilds cleanly today.**
  `tests/data/structures/6kn7.pdb` (29 chains, 61,511 atoms, 7,781 residues) rebuilds in **6.2 s**
  with `report.ok == True`, 25/25 regions built, because no chain has two folded domains so step 3
  is a no-op. This is the performance and correctness baseline.
* **117-structure corpus and 23,587-structure proteome run.** Every one is single-chain. Any
  change to the single-chain code path is a regression risk against work that is already
  validated; section 3.6 makes preserving it an explicit, testable requirement.

---

## 1. The model: rigid units

Introduce one concept, and route everything through it.

> A **rigid unit** is a set of folded domains -- possibly from different chains -- whose relative
> positions and orientations DODO will not change. A unit moves as one rigid body or not at all.

Today's behaviour is the special case "every folded domain is its own unit". The complex case is
"every folded domain in the assembly is in one unit". Both fall out of the same code.

Why a set of domains rather than a set of chains: a chain in a complex routinely has some domains
in the assembly and some dangling off a long linker. Locking whole chains would forfeit step 3 --
the step that makes DODO work -- for every chain that touches another. Locking domains keeps it
for exactly the domains that are free to move.

**The invariant, stated so it can be tested:** for any two atoms in the same unit, the distance
between them after a rebuild equals the distance before it, to 1e-6 A. This subsumes "the
interface is preserved", "internal geometry is preserved", and "the peptide bond between adjacent
folded domains survives".

---

## 2. Building the units

Three rules, applied in this order, then a union-find. All three are evaluated on the **input**
coordinates, before anything moves.

### 2.1 Covalent continuity (hard, non-negotiable)

Two folded domains in the same chain with **no residues between them** are in the same unit.
There is a peptide bond across the boundary; separating them breaks it. This is defect 0.3a.

Extends to: two folded domains separated only by residues that **will not be rebuilt** -- a region
under `min_length`, or one the caller has pinned. Those residues keep their input coordinates, so
they are a rigid link between their neighbours, and the flanking domains are one unit. This is
defect 0.3b.

The two cases are the same rule: *if the geometry between two folded domains is not going to be
regenerated, that geometry is part of the rigid body.* Note this makes unit construction depend on
`min_length`, which is a `rebuild()` argument -- so units are built inside the pipeline where
`min_length` is known, not in region identification.

### 2.2 Interfaces (the complex case)

Two folded domains are in the same unit if they are in physical contact. Proposed measure, to be
calibrated in phase C-0:

* count residue pairs (i in domain u, j in domain v) with any heavy-atom pair within
  `INTERFACE_CONTACT_RADIUS = 5.0` A -- the standard interface definition, and the one that makes
  the number comparable to published buried-surface estimates;
* lock when the count reaches `INTERFACE_MIN_RESIDUE_PAIRS`.

**The threshold should be low, and the asymmetry is the argument.** Over-locking costs accuracy of
a kind DODO can afford: a unit that could have been repositioned stays where the input put it, and
the output is closer to the input than it needed to be. Under-locking costs the thing this feature
exists to prevent: a real complex is dismantled and reported clean. So the default errs toward
locking, and loosening it is the user's call, not ours. Starting proposal: 3 residue pairs.

Restricted **by default to pairs of domains in different chains** (decision 1). Two folded domains
of the same chain that pack against each other are exactly the case DODO exists to re-sample --
AlphaFold's inter-domain placement across a long linker is what step 3 replaces -- and locking them
by default would change the behaviour validated over 23,587 single-chain structures. Available as
`lock_intra_chain_interfaces=True` for callers who want it, and measured in C-5 before any proposal
to change the default.

By contact, not by chain: this is what handles a bound SLiM. A short segment of chain A that
threads across chain B's surface is already classified folded by the burial score (measured in
0.4), and this rule locks it to B. It is then not rebuilt, the interaction survives, and the
linker back to chain A's own domains becomes a placement constraint -- which is exactly right.

### 2.3 Caller override

`units=<spec>` names units explicitly, in the same style as `assign_regions_from_spec`: chain id
plus inclusive input residue numbers. Also `units="none"` for today's behaviour, and
`units="chains"` (one unit per chain's folded domains) for people who know what they have.

### 2.4 Transitivity is intended

Union-find over the locked pairs. On a filament or a large assembly this merges everything into
one unit and step 3 becomes a no-op -- which is the correct answer for 6kn7. The failure mode to
watch is a chain of weak contacts merging two things that should have been free; the report
(section 7) prints the contact count for every locked pair so the user can see which one did it,
and `dodo units` (section 8) shows the grouping before a build is run.

### 2.5 Determinism

Units are ordered by their earliest domain (chain index, then residue index), and the list of
domains inside a unit likewise. Every downstream loop iterates that order, so the RNG sequence is
a function of the input and the seed only.

---

## 3. Placing the units

### 3.1 The graph

* **Nodes:** units.
* **Edges:** one per connecting IDR. For consecutive folded domains `u`, `v` in a chain with a
  single IDR between them, the edge carries: attachment point `a` = exit CA of `u`, `b` = entry CA
  of `v`, target `d = target_dimensions(linker.sequence, mode).end_to_end`, and hard bounds
  `min_reach(n+1)`, `max_reach(n+1)` from `engines.walk`.
* An edge whose endpoints are in the **same** unit is not a constraint -- the separation is already
  decided. Record it, compare the achieved span to the prediction, and disclose (section 7). If the
  achieved span exceeds `max_reach`, the linker is unbuildable; say so up front rather than letting
  the engine raise mid-build.
* **Terminal IDRs contribute no edges.** A free end constrains nothing.
* **Edges only ever exist within a chain**, because only a chain is covalent. Two units that share
  no chain have no geometric relationship DODO is entitled to invent.

Consequence worth stating plainly: in a typical AF3 complex where every folded domain packs into a
single assembly, there is one unit, no edges, and **step 3 does nothing**. The whole feature, for
the most common case, is "detect that and skip".

### 3.2 Components

Place each connected component independently. Within a component, the **reference unit** is the one
containing the earliest folded domain (chain order, then residue index) -- which for a single chain
of singleton units is exactly today's "first folded domain of the chain; it defines the frame of
reference". Units in other components do not move: nothing constrains them. They are obstacles, and
if a moved unit ends up overlapping one, that is disclosed rather than silently accepted.

("Anchor the largest unit instead" would move less mass. It also changes the output frame for every
existing single-chain input. Worth measuring in C-5; not worth the regression now.)

### 3.3 One unit, one constraint (k = 1)

Today's algorithm, generalised from a domain to a unit:

1. draw a random direction, translate the unit so `b` sits `d` A from `a`;
2. rotate the unit about `b` so the unit's centroid points away from `a`, then perturb by up to
   `PERTURBATION_DEGREES`;
3. re-place, clash-check against everything already final, retry to `MAX_FD_PLACEMENT_ATTEMPTS`,
   keep the least-bad and report it.

Two changes, both required: the centroid is the **unit's**, not the entry domain's; and rotations
are applied to every domain of the unit about a **single common point**, never per domain.

For a single chain whose units are all singletons, this must consume the RNG in the same order and
produce the same coordinates as today. That is a test, not an aspiration (3.6).

### 3.4 One unit, several constraints (k >= 2)

A unit reached by two or more linkers from already-placed units cannot be placed by choosing one
direction. Two heterodimer chains that both contribute a domain to each of two units is the common
case, and it is not exotic.

**Sphere projection + Kabsch, iterated.** Each constraint says "attachment point `b_i` must sit on
the sphere of radius `d_i` about `q_i`". Given a current placement:

1. for each constraint, project the current position of `b_i` onto its sphere:
   `t_i = q_i + d_i * unit(T b_i - q_i)`;
2. find the rigid transform taking `{b_i}` to `{t_i}` in least squares -- Kabsch;
3. apply it, and repeat until the worst residual `| ||T b_i - q_i|| - d_i |` stops improving.

It converges in a handful of iterations, needs no gradients, and degenerates correctly: with k=1 it
is a pure translation onto the sphere (which is why k=1 keeps the dedicated placer, for the
orientation heuristic); with k=2 it lands on the exact solution when one exists; with k>=3 it
returns the least-squares compromise, which is the honest answer for an over-determined system.

Initialise from the k=1 placer on the first constraint, so the search starts somewhere sane; on a
clash, re-randomise that initial direction and re-converge. Budget restarts the same way
`MAX_FD_PLACEMENT_ATTEMPTS` is budgeted today, keep the best by (fewest clashing atoms, then
smallest residual), and report both.

Deferred until measured: a closed-form k=2 solver. With two constraints the feasible set is exactly
a 3-parameter family plus a free rotation about the `t_1`-`t_2` axis, so an exact solver would give
zero residual and reduce the clash search to a 1-D scan. Implement the general iteration first;
add the special case only if the measurements show the iteration leaving residuals or the clash
search struggling on real dimers.

### 3.5 Cycles

Two units joined by two independent linker paths -- a homodimer of a two-domain protein where both
domains dimerise, or a chain that folds back so its first and third domains share a unit -- give a
cycle. Greedy placement over a spanning tree cannot satisfy the closing edge.

Handle it with the same primitive: after the tree placement, run **block-coordinate descent** over
the component. Repeatedly pick each non-reference unit and apply 3.4 using *all* of its edges, not
just the tree edge. Sweep until the worst residual converges or the budget runs out.

**Skip the sweep entirely when the tree placement already satisfies every edge within tolerance.**
Every acyclic component satisfies this, which is every single-chain input, which is how the
existing behaviour survives untouched.

Feasibility first: a cycle whose target distances violate the triangle inequality has no solution,
and iterating is a waste. Check the cycle basis up front, report an infeasible cycle by name, and
fall back to weighted least squares with the residual disclosed.

**A softer constraint model, worth measuring in C-3.** The predicted end-to-end distance is the
*mean* of a distribution, and a linker can be shorter than its mean far more easily than it can be
longer than its contour length. Treating each edge as "soft target `d`, hard ceiling
`max_reach(n+1)`" rather than as an equality makes over-determined and cyclic systems far more
often feasible, at no cost in physical realism. Measure the residual distribution both ways before
choosing.

### 3.6 The regression requirement

Everything above must reduce, exactly, to what ships today when the input is a single chain whose
folded domains form no units.

* **Byte-identical output** for `dnmt3a.pdb`, `p300.pdb`, `arf19.pdb`, `testing_translation.pdb`
  and `test.pdb` at seeds 0-2, with and without `--backbone`, before and after the change.
* **The 117-structure corpus** and the 9 historical-failure structures, unchanged.
* This is achievable because the only single-chain code path that changes is "iterate over units of
  one domain each" instead of "iterate over domains", with the same placement function, the same
  RNG draws in the same order, and the relaxation sweep skipped. If it comes out non-identical,
  that is a bug in this work, not an acceptable cost.

---

## 4. Clashes

### 4.1 Correctness

The obstacle set for placing a unit is **every atom already final, in every chain**: units already
placed, units in other components that will never move, and folded domains held as reference. It is
the whole-structure `placed_atom_mask` idea that `_obstacles_for_span` already uses, applied to step
3 where today it is per-chain. This alone fixes 0.2b.

### 4.2 Cost, and what to do about it

`Structure.clash_mask` rebuilds its KD-tree on every call. Measured on 6kn7, 50,000 obstacle atoms
against a 3,000-atom query:

```
clash_mask (rebuilds the KD-tree every call): 7.6 ms/attempt
prebuilt tree: build 5.8 ms once, query 0.69 ms/attempt
```

11x per attempt. At 500 attempts that is 3.8 s versus 0.35 s, for one unit. Build the obstacle tree
once per placement and query the moving unit's atoms against it.

Three more, in order of expected value, each to be measured rather than assumed:

* **Bounding spheres.** Precompute a centroid and radius per unit; skip the tree query entirely when
  the centres are further apart than the sum of radii plus the cutoff. In a large assembly most
  pairs never touch.
* **Coarse-then-exact.** Search on alpha carbons with an inflated cutoff, verify the accepted
  placement on all atoms. Roughly 8x fewer query points on a full-atom structure. Must be verified
  not to change accepted placements, or it is not worth having.
* **Query the smaller set.** Whichever of (moving unit, obstacles) has fewer atoms should be the
  query, not the tree.

### 4.3 Big units are harder to place

A unit of several domains sweeps a much larger volume than a single domain, so a random direction
clears less often. Mitigations, in order: bias the sampled directions away from the centroid of the
already-placed geometry rather than sampling the sphere uniformly; spend the budget on the *free*
degrees of freedom (the rotation about the constraint axis for k=2 is one cheap scan) rather than on
fresh random draws; keep the existing least-bad fallback and disclosure. This is the same lesson as
the anchor-free region placement work -- the budget belongs in orientations, not in draws.

---

## 5. Region identification in a complex

No change to the algorithm; three things to document and one to check.

* Burial is scored whole-structure, so interface residues read as folded (measured in 0.4). Keep it.
* **The same chain gives different regions alone and in a complex.** All indices here are 0-based
  positional, not residue numbers. dnmt3a's first folded domain runs 282-432 as a monomer and
  280-433 in the dimer, and a 30-residue interface helix at 192-222 appears that the monomer does
  not have. This is correct and it will surprise people. Document it in the guide, and print it in
  the report.
* **`--strategy plddt` is interface-blind.** Per-residue pLDDT says nothing about a partner chain, so
  a bound SLiM stays an IDR and gets rebuilt off the interface. Warn when `plddt` is combined with a
  multi-chain input.
* **To check:** whether the contact radii merge genuinely separate chains that merely pass close --
  8 A all-atom for the default `density` strategy, 12 A CA-CA for `contact`. Both were tuned on
  single chains, where the only thing at that distance is the same molecule. If they over-merge,
  region boundaries in dense assemblies are wrong before units are ever built.

---

## 6. Build order and obstacles during the rebuild

Mostly already right; two things to measure.

* **Terminal-tail starvation.** Tails are built last, all of them, in residue order across chains --
  so in a homo-oligomer the last chain's tail is always built into whatever space the others left.
  Measure the per-chain failure rate on a 10+ chain complex. If chain index predicts failure,
  interleave by chain or build shortest-first, and measure again. This is the "bowl of spaghetti"
  concern from ROADMAP item 2, arriving early.
* **Obstacle set size.** Every IDR build in a 60,000-atom assembly builds a tree over most of it.
  Profile it on 6kn7 against the 6.2 s baseline; the same prebuilt-tree fix from 4.2 may apply.

---

## 7. Reporting

The codebase's standing rule is that a thing DODO could not do is labelled, not hidden. Additions:

* `UnitReport`: for each unit -- id, the chains and domains in it, atom count, whether it moved, and
  the reason it is a unit (covalent continuity / interface with N residue-pair contacts / caller
  spec). Enough for a user to see why something did or did not move.
* `InterfaceReport`: every locked pair, its contact count, and confirmation that the contact map is
  unchanged after the build.
* Per-edge linker reporting: target separation, achieved separation, residual, and for an intra-unit
  linker an explicit `separation dictated by the complex; predicted X A, actual Y A`.
* **Unbridgeable linker**, named as such when the span a locked unit forces exceeds `max_reach`.
* **`report.ok` must go False** when a unit-level invariant fails: a broken covalent link, an
  introduced impossible contact between units, or an unresolvable clash. The dimer measured in 0.1
  must not report `ok: True`.
* Related, and worth a decision while this code is open: a `--backbone` seam of 141.83 A is
  disclosed but does not count against `ok`, and it is not the same object as a 2.2 A seam. Consider
  a magnitude above which a seam is a failure rather than a strain.

---

## 8. API and CLI

```python
rebuild(
    source,
    units="auto",                     # "auto" | "none" | "chains" | {spec}
    interface_contact_radius=5.0,
    interface_min_residue_pairs=3,
    lock_intra_chain_interfaces=False,
    ...
)
```

```
dodo rebuild complex.cif --units auto|none|chains
dodo rebuild complex.cif --units-from units.json
dodo units complex.cif          # print the detected units and interfaces, build nothing
```

`dodo units` matters more than it looks. Unit detection is a judgement call on real data, and the
user needs to see and correct the grouping before spending a build on it -- the same role `dodo
regions` plays for region assignment.

New module `src/dodo/construct/assembly.py` (`RigidUnit`, `find_rigid_units`, `units_from_spec`,
`InterfaceReport`); `src/dodo/construct/place.py` reworked to place units, keeping
`reposition_folded_domains` as the entry point; `superpose` (Kabsch) added to
`src/dodo/geometry/transforms.py`, which has no least-squares superposition today; new constants in
`src/dodo/constants.py`, each documented with its measurement the way the existing ones are.

---

## 9. Phases

Each phase ends shippable and is judged on measurements, not on being finished.

### C-0 -- Fixtures and baseline

Build the multi-chain corpus (section 10). Record, for each: interface contact counts, `report.ok`,
runtime, per-chain region assignment, and the damage the current code does. Calibrate
`INTERFACE_MIN_RESIDUE_PAIRS` against the deposited biological assemblies of the experimental
fixtures -- print the distribution of pairwise contact counts and pick the threshold at the observed
gap, rather than choosing a number and defending it afterwards.

*Done when:* the harness reproduces every number in section 0 and fails on today's code.

### C-1 -- Rigid units, freeze only

Build the unit graph (section 2), expose `dodo units`, and make step 3 refuse to break a unit: any
unit of more than one domain is not moved at all. Singleton units are placed exactly as today, but
against a whole-structure obstacle set.

This is small, it is low-risk, and it delivers the headline requirement plus both defects in 0.3.
The only thing it gives up is repositioning a multi-domain unit that legitimately should move.

*Done when:* interface contact maps are bit-preserved on every fixture; the 93.479 A and 33.605 A
tears are gone; single-chain output is byte-identical at seeds 0-2; the dimer no longer reports
`ok: True`.

### C-2 -- Units move

Generalise the placer from domain to unit (3.3), with a common rotation centre and unit-level
`verify_rigid`. Add the prebuilt-tree and bounding-sphere work from 4.2 and report the speed
change.

*Done when:* a multi-domain unit hanging off a long linker is placed at its predicted separation
with its interfaces intact; single-chain output still byte-identical; 6kn7 no slower than 6.2 s.

### C-3 -- Multiple constraints and cycles

Kabsch relaxation (3.4), component sweeps (3.5), feasibility checks, soft-target/hard-ceiling
measurement.

*Done when:* the two-interface heterodimer fixture places with both linker residuals inside
tolerance; an infeasible cycle is named rather than thrashed; single-chain output still
byte-identical.

### C-4 -- Reporting, validation, docs

Section 7 in full; `dodo validate` gains a unit-integrity check; guide and README sections; the
build-order measurements from section 6.

*Done when:* every invariant in section 11 is asserted by a test, and the guide documents the
monomer-versus-complex region difference.

### C-5 -- Policy questions, measured

Intra-chain interface locking; largest-unit-as-reference; optional AF PAE input as the locking
signal for predicted complexes (inter-domain PAE is the *right* answer to "is this relative
placement confident", and it is what the interface heuristic is approximating). Each one lands only
with numbers attached.

### C-6 -- Missing residues (ROADMAP item 1, case 2)

Independent of everything above and shippable separately; deliberately not blocking this work
(decision 3).

An experimental complex is missing its IDRs entirely, so there is nothing to rebuild until the
residues exist. Needs: `Chain.full_sequence` (already parsed from SEQRES and `_entity_poly`) or a
caller-supplied FASTA chain map; an alignment of observed against full sequence; and a new
`insert_unmodelled_residues` that rebuilds the atom arrays through
`Structure.from_atom_records`, marking the inserted residues so they are always rebuilt, never
obstacles until built, and excluded from burial scoring and interface detection. Note the mmCIF
reader does not currently parse `_pdbx_poly_seq_scheme`, which is the only exact mapping from
deposited sequence position to author numbering; without it, internal gaps have to be inferred from
gaps in `residue_number`, which is workable but not equivalent.

The interaction with this plan is small and favourable: in case 2 everything present is one rigid
unit, so DODO freezes the complex and fills the gaps.

---

## 10. Fixtures

Committed fixtures today are single-chain except 6kn7, which has one folded domain per chain and so
exercises none of this. Needed:

| Fixture | What it is for |
|---|---|
| `dnmt3a_dimer.pdb` (synthetic, script committed) | Two chains, two domains each, one interface. The minimal case; reproduces every number in 0.1. |
| A real AF3 two-chain prediction with tails and linkers | The primary target use case. |
| A two-interface heterodimer | The cycle case (3.5). Can be synthesised the same way as the dimer. |
| A complex with one chain in the assembly and one free | Components (3.2); the free chain must not move and must not be built through. |
| A bound-SLiM complex | A short segment of one chain across another's surface; must be locked, not rebuilt. |
| `6kn7.pdb` | Regression: one unit, no edges, step 3 a no-op, 6.2 s. |
| A 10+ chain homo-oligomer with tails | Terminal-tail starvation (section 6). |

Synthetic fixtures are legitimate here and preferable for the geometric cases: they isolate one
property, and the generator script is the documentation of what is being tested. The AF3 and
experimental fixtures are what keep the thresholds honest.

---

## 11. Invariants to assert

1. **Unit rigidity.** `verify_rigid` over each unit's concatenated atoms, 1e-6 A.
2. **Interface preservation.** Every atom pair within 5 A in the input and in the same unit is at
   the same distance, to 1e-6 A, in the output.
3. **Covalent continuity.** Every consecutive-residue C-N bond at or under 1.90 A in the input is at
   or under 1.90 A in the output, unless the residues between were rebuilt.
4. **No new impossible contacts** between units (< 1.00 A), measured differentially against the
   input, per the corpus convention.
5. **Linker targets.** Inter-unit linkers land within tolerance of their target; intra-unit linkers
   are reported as dictated, with the prediction alongside.
6. **No single-chain regression.** Byte-identical output on every committed single-chain fixture at
   seeds 0-2, both backbone modes, plus the corpus and historical suites.
7. **Determinism.** Byte-identical re-run at a fixed seed, on multi-chain input.
8. **Honest `ok`.** A destroyed interface, a broken bond, or an unresolved inter-unit clash makes
   `report.ok` False.

---

## 12. Out of scope, stated so it is not discovered later

* **Nucleic acids -- protein-only this release (decision 2).** A protein-DNA or protein-RNA complex
  cannot be read at all: nucleotides come in as polymer residues with no CA, and
  `Structure.ca_indices` raises `EmptyStructureError`. Common in real complexes. For now, refuse
  with a message that names nucleic acids as the reason rather than the generic missing-CA error.
  The later answer is a "rigid passenger" concept -- non-protein entities that are read, never
  rebuilt, and move with the nearest unit. The rigid-unit abstraction makes that natural, but it
  touches the reader and both writers.
* **Ligands and metals are silently dropped.** The reader skips HETATM residues outside the polymer
  whitelist; 6kn7 loses its 15 ADP (7,796 residues deposited, 7,781 read). A zinc or nucleotide at
  an interface disappears from the output. This release makes the drop explicit in the report rather
  than a parser note; carrying them as rigid passengers is the same later phase as nucleic acids.
* ROADMAP items 2 (multiple copies per frame) and 3 (dynamic folded domains) are separate, though
  both will want rigid units once they exist.

---

## 13. One unrelated thing found while reading

`src/dodo/construct/pipeline.py` around line 1005 contains a block labelled
`# MUTATION EXPERIMENT (temporary)` which builds loops outside the `attempts` list, so
`_repair_forward_contacts` never sees a loop. It is committed (ca1e057) and the comment says it is
the revert of a reviewed fix. Not part of this work, but it should not stay.
