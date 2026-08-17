# How DODO works

Seven steps. The order is not incidental — it *is* the algorithm.

1. **Identify** folded domains, IDRs, and loops, using the all-atom density metric. A *loop* is a
   disordered stretch tethered at both ends inside a **single** connected folded domain, which
   makes it a different problem from an IDR joining two domains.
2. **Predict** each IDR's end-to-end distance from its sequence with ALBATROSS. Loops get no
   prediction — see below.
3. **Reposition the folded domains.** Translate and rotate each one as a rigid body so that
   consecutive domains sit at the predicted end-to-end distance of the IDR between them.
4. **Rebuild loops.**
5. **Rebuild connecting IDRs** between folded domains.
6. **Rebuild terminal IDRs.**
7. **Place the backbone.** N, C and O are added to the rebuilt regions by default;
   `--no-backbone` leaves them alpha-carbon only. Side chains are not built. See {doc}`guide`.

## Why folded domains are never rebuilt

DODO does not regenerate a single folded-domain atom. It only translates and rotates them, as
rigid bodies. Whatever AlphaFold produced — side-chain positions, bond lengths, the lot — comes
through exactly, to the last decimal place.

This has a consequence worth stating plainly: **DODO faithfully preserves defects in its input.**
AlphaFold ships genuinely broken residues; one of the structures in DODO's own test corpus contains
92 pairs of atoms closer together than the shortest bond in any protein, one of them at 0.046 Å.
DODO does not repair those and does not claim to. The validation suite therefore distinguishes
geometry DODO *generated* from geometry it *inherited* — see {doc}`validation`.

## Why step 3 is the one that matters

This is the step that is easy to miss, and an earlier version of the v2 rewrite missed it entirely.

AlphaFold has no way to know how far apart two domains should sit when a disordered linker joins
them, so it packs them far closer than the linker's sequence implies. Measured on real models,
**2–3.6× closer**. On p300, a 151-residue linker connects domains sitting 26.1 Å apart, where the
predicted end-to-end distance for that sequence is 94.9 Å.

Skip step 3 and the linker gets built into a 26 Å gap. The result is a compact blob wedged between
two domains — which is not what the sequence says, and is arguably worse than the spaghetti it
replaced, because it looks deliberate.

Domains are placed by orienting each toward the next, perturbing the orientation, and rejecting
arrangements that clash. Internal geometry is verified to survive the motion: `verify_rigid` asserts
each atom's distance to the domain centroid is unchanged to within 1e-6 Å, and the measured drift is
around 1e-13 Å — floating-point noise rather than deformation.

## Why loops get no dimension prediction

A loop's ends are pinned to two anchors inside one folded domain, and that domain is not DODO's to
move. The distance is therefore already decided by geometry. Predicting an end-to-end distance for
it and then failing to achieve it would be inventing a constraint that the structure has already
settled.

Connecting IDRs are different: their anchors belong to two *different* domains, and step 3 has
already moved those domains to the predicted separation.

## Why the build order is loops, then connecting, then terminal

Decreasing order of constraint.

- A **loop** is pinned at both ends inside one domain that cannot move. Least freedom.
- A **connecting IDR** is pinned at both ends, but between two domains that were positioned *for
  it* in step 3.
- A **terminal IDR** is free at one end and can go almost anywhere. Most freedom.

Build the loose regions first and they occupy space the tight ones then cannot avoid. Building in
residue order instead — which an earlier version did — lets a floppy terminal tail wander through
the space a loop needs.

## Region identification

The default is DODO's original **all-atom density** metric: all-atom pairs within 8 Å, per residue,
thresholded at 480. This is the method the package was built and validated on; the author reports it
draws better boundaries than sequence-based disorder predictors, though that comparison predates this
repository and is not benchmarked here. It is reimplemented over a KD-tree rather than
changed — same numbers, 10.1 s down to 7 ms on a 1,086-residue model.

Alternatives, selectable with `-s` / `strategy=`:

`density`
: The default and the validated method. Needs side chains to count.

`contact`
: A CA-only burial score. Composition-free and invariant to whether side chains are modelled,
  which the density score is not. Not the validated method; useful for comparison and for CA-only
  input.

`plddt`
: AlphaFold's own per-residue confidence, from the B-factor column. An explicit opt-in — density
  is DODO's validated default.

`auto`
: `density` for all-atom input, `contact` for CA-only input, where a pair count cannot be compared
  against the tuned threshold. It reports which it chose. It deliberately never picks `plddt`.

`preset`
: Identify nothing; build the regions already attached to the structure. See {doc}`guide`.

Every assignment keeps its score profile and threshold, so a boundary you disagree with can be
audited rather than just re-run with different numbers.

## What the disordered regions are built with

A self-avoiding growth walk, angle-constrained. Consecutive alpha carbons are placed 3.81 Å apart,
and the CA–CA–CA pseudo-angle is held inside a 91–161° window measured from real structures.

Candidate positions are ordered ideal-first and rejected if they come within 3.20 Å of anything
already placed. There is no relaxation of that threshold: an earlier version dropped it as far as
2.00 Å when no candidate fitted, on the theory that a tight contact beats a failed region. Measured
across three structures at four seeds, that ladder caused 69 of 79 steric clashes and prevented
zero region failures, so it was deleted. If a region genuinely cannot be built at 3.20 Å, DODO
reports the failure and leaves the input coordinates alone.

Each model in a multi-model run draws its own target from the physical distribution
P(R) ∝ R²exp(−3R²/2⟨R²⟩), so an ensemble genuinely varies in size. Folded domains are positioned
once and held fixed across every model, so a viewer flicking between frames sees only the
disordered regions move.

## The STARLING engine, and what it does not know

`--engine starling` replaces the self-avoiding walk with conformers from
[STARLING](https://github.com/idptools/starling), a generative model of disordered ensembles. It is
the only optional dependency, because it is large — roughly 2.4 GB of weights.

### Loops always use the walk engine

`--engine starling` applies to terminal and connecting IDRs only. Loops are built by the
self-avoiding walk regardless, and this is structural rather than a gap to be closed.

A loop is pinned at both ends inside a **single** folded domain, and that domain is not DODO's to
move, so the loop's end-to-end distance is decided exactly by geometry before anything is generated.
STARLING conditions on sequence alone and returns a *distribution* — measured on dnmt3a's regions,
its end-to-end distances scatter with a standard deviation of 41% of the mean — so the chance a
conformer lands on the one distance the anchors permit is small, and unlike a connecting IDR there
is no second lever, because the domain cannot be moved to suit the conformer.

Asking STARLING for loops anyway is what produced `residues 1712-1726: NOT BUILT (loop in FD
1578-1830)` on p300: every loop failed while the IDRs around them built.

### STARLING models IDRs *alone*

This is the most important thing to understand about the engine, and you cannot see it in the
output, so DODO warns about it on every run.

**STARLING is given a sequence and nothing else.** Not the folded domains, not their positions, not
the space they occupy. It was trained on isolated disordered regions and that is what it models. So
a STARLING conformer:

- cannot avoid the folded domains, because it was never shown them;
- has conformational statistics that are **not** conditioned on them either.

What DODO does with that conformer is pick the one whose own end-to-end distance best matches what
the anchors demand, then place it as a rigid body so its first alpha carbon sits one bond from the
N-anchor and its last one bond from the C-anchor. Regions it cannot fit are reported, never forced.

So the region's **internal geometry is STARLING's** and its **placement is DODO's**. Read a STARLING
region as a realistic IDR conformation that has been positioned — not as one sampled in the presence
of the domains it sits between. If that distinction matters for what you are doing, it matters a
lot.

The walk engine is the opposite trade: its conformations are geometric rather than learned, but it
is aware of every already-placed atom and avoids all of them.

### Regions longer than 380 residues

STARLING will not generate a region longer than 380 residues, and real IDRs routinely are — p300's
disordered N-terminus alone is over a thousand. DODO handles this rather than erroring or silently
downgrading to a random walk.

The region is split into segments within the cap, each generated separately, and the segments are
then arranged in space. That is sound on polymer-scaling grounds rather than being a workaround: for
a chain with Flory exponent `ν`, a segment of `N/k` residues has end-to-end distance `~(N/k)^ν`, and
arranging `k` of them self-avoidingly with a step of that order gives `~(N/k)^ν × k^ν = N^ν` overall.
The assembled chain scales with length exactly as one long chain would.

Two things the assembly has to get right, because independently generated segments know nothing
about each other:

- **Segments overlap where they join.** Adjacent segments share residues, and the splice is chosen
  from candidate rotations about the join, so the chain is continuous rather than merely adjacent.
- **Clashes are checked between segments as well as within them.** A segment generated in isolation
  has no idea another segment occupies the same space, so the arrangement measures and rejects.

### Bond lengths are corrected, not just checked

STARLING is a diffusion model that reconstructs coordinates from a predicted distance map, so its
virtual CA–CA bonds scatter around 3.8 Å rather than sitting on it. That is normal model output, but
it is not a protein — the trans-peptide CA–CA distance is rigid.

So DODO projects every conformer onto an exact 3.81 Å bond before using it, with a SHAKE-style
iterative correction. Screening alone, which is what this used to do, either rejects usable
conformers for ordinary diffusion noise or passes that noise into the output file. On synthetic
diffusion-like traces the worst bond deviation goes from 0.19–1.02 Å to 0.0000 Å, moving atoms by
0.12–0.72 Å depending on how noisy the input was.

## What comes next

In priority order:

1. Perfecting the alpha-carbon approach — this, and it is where the effort has gone.
2. Backbone building (N, C, O) for rebuilt regions. **Shipped and on by default** — see
   {doc}`guide`. Held out against all-atom simulation: N 0.16 Å, C 0.22 Å, O 0.63 Å, with every
   bond length inside a region exact. Its one irreducible limit is the seams: an exact peptide bond
   onto an untouched folded domain is geometrically unsatisfiable from a rebuilt alpha carbon, so
   that bond is drawn as close as possible, left long, and reported. `--no-backbone` returns
   alpha-carbon-only output.
3. Full all-atom reconstruction.
4. STARLING-generated ensembles as the default engine.
