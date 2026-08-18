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

## What comes next

Shipped, and where the effort has gone: the alpha-carbon approach itself, and backbone building
(N, C, O) for rebuilt regions — on by default, held out against all-atom simulation at N 0.16 Å,
C 0.22 Å, O 0.63 Å, with every bond length inside a region exact. The peptide-plane lookup is keyed
on five consecutive alpha carbons, worth 5.1% on C and 3.8% on N over the four-carbon form. Its one
irreducible limit is the seams: an exact peptide bond onto an untouched folded domain is
geometrically unsatisfiable from a rebuilt alpha carbon, so that bond is drawn as close as possible,
left long, and reported. `--no-backbone` returns alpha-carbon-only output.

Still ahead, in priority order:

1. Full all-atom reconstruction — side chains, which DODO does not build today.
2. STARLING-generated ensembles as an engine option.
