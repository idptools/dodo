r"""All-atom backbone reconstruction from a CA-only trace.

What this module does
---------------------
Given alpha-carbon coordinates and a sequence, place ``N``, ``C`` and ``O`` for every
residue -- and optionally ``CB`` -- so that the result is a physically valid protein
backbone rather than a CA trace with extra dots in it. This is what closes three of the
four limitations DODO's own README lists: the ``unusual bond`` warnings VMD emits for a
3.8 A "bond", the broken tube and cartoon representations, and CA-only IDRs.

Why it was written from nothing
------------------------------
v1's all-atom module is dead code that its own docstring calls "hot garbage", and the
description is accurate. It added *unit* vectors without scaling them by a bond length,
so CA-C came out at exactly 1.000 A against a real 1.53 and the peptide bond at 2.87 A
against 1.33 -- while the same file declared the correct constants and never imported
them. Its per-residue side-chain templates were worse: every template had N-CA = 1.220 A
against the 1.46 the file itself declared, the arginine template carried a C-terminal
``OXT`` and had ``NH1`` coordinates copy-pasted out of glutamate's ``OE2`` slot, and the
only consumer applied templates by pure translation with no rotation, so every side chain
pointed the same absolute lab direction. None of it is transcribed here. Everything below
is built from internal coordinates and the constants in :mod:`dodo.constants`.

The geometry, derived
---------------------
**Step 1: the peptide unit is rigid, and it fixes the CA-CA virtual bond.**

Five atoms make up one peptide unit: ``CA(i)``, ``C(i)``, ``O(i)``, ``N(i+1)``,
``CA(i+1)``. With omega (the ``CA(i)-C(i)-N(i+1)-CA(i+1)`` dihedral) trans, the unit is
planar, and once the three bond lengths ``CA-C``, ``C-N``, ``N-CA`` and the two bond
angles ``CA-C-N``, ``C-N-CA`` are fixed the whole unit is a rigid body. Building it from
the declared constants gives a ``CA(i)-CA(i+1)`` distance of **3.8040 A**, against a
declared :data:`~dodo.constants.CA_CA_BOND_LENGTH` of 3.81 -- 0.006 A apart, which is inside
the peptide-angle strain this module already applies to close onto a real trace (see
:data:`_MAX_PEPTIDE_ANGLE_STRAIN`, whose window is 3.659-3.932 A). That agreement is not a
coincidence, it is the consistency check this module rests on: a correct reconstruction
*closes onto the CA trace it was handed*, so the input alpha carbons come out unchanged
rather than merely approximately reproduced.

Placing the O needs one derived quantity. The carbonyl carbon is sp2, so its three
substituents ``CA``, ``N`` and ``O`` are coplanar with it and their three angles sum to
360 degrees. Hence ``N-C-O`` = 360 - :data:`~dodo.constants.CA_C_N_ANGLE` -
:data:`~dodo.constants.CA_C_O_ANGLE` = 123.0 degrees, and the O is placed *cis* to
``CA(i+1)`` about the ``C-N`` bond, i.e. on the opposite side of the peptide plane from
``CA(i)``. That reproduces :data:`~dodo.constants.CA_C_O_ANGLE` exactly, which is
asserted in the tests.

**Step 2: one free rotation per peptide unit.**

Both alpha carbons lie in the unit's plane, so knowing them fixes the unit up to a
rotation about the ``CA(i)-CA(i+1)`` axis. Write that angle ``theta_j`` for unit ``j``.
In an orthonormal frame ``(xhat_j, e1_j, e2_j)`` with ``xhat_j`` along
``CA(j+1) - CA(j)``, every atom of the unit sits at ``CA(j) + px * xhat_j + py *
yhat_j(theta_j)`` with ``yhat_j(theta) = cos(theta) * e1_j + sin(theta) * e2_j`` and
``(px, py)`` the atom's fixed in-plane coordinates. So the whole reconstruction is
``n - 1`` scalar unknowns.

**Step 3: what constrains those unknowns.**

``N(i)`` comes from unit ``i-1`` and ``C(i)`` from unit ``i``, so the intra-residue
``N-CA-C`` angle (tau) at residue ``i`` couples exactly two consecutive unknowns. There
are ``n - 2`` such couplings for ``n - 1`` unknowns, which leaves a one-parameter family
-- and worse, each coupling has two roots (the classic peptide-plane flip), so tau alone
admits on the order of ``2**n`` distinct solutions that all satisfy it. Tau is therefore
necessary and *not* sufficient, which is why the literature methods (PULCHRA, BBQ,
SABBAC) key a statistical library on the CA pseudo-dihedral. DODO ships no such library
and this module does not invent one.

What it uses instead is the other thing that is genuinely determined by geometry:
**hard-sphere sterics**, which is how Ramachandran derived the allowed phi/psi regions in
the first place. Both phi and psi at residue ``i`` are functions of ``(theta_{i-1},
theta_i)`` and nothing else -- phi needs ``C(i-1)``, ``N(i)``, ``CA(i)``, ``C(i)`` and
psi needs ``N(i)``, ``CA(i)``, ``C(i)``, ``N(i+1)`` -- so a steric term stays exactly as
pairwise as the tau term. Residue ``i``'s ``CB`` is likewise a function of the same two
unknowns, which matters because ``CB`` is the atom that makes the disallowed regions
disallowed.

**Step 4: solve it exactly on a grid.**

Because every term couples only consecutive unknowns, the problem is a linear chain and
Viterbi dynamic programming finds the *global* optimum over a discretized ``theta``
without any local-minimum risk and without an iterative sampler. Cost is
``O(n * grid**2)``, vectorized as one small matrix product per residue.

**Step 5: the part step 3 is structurally blind to.**

Steps 3 and 4 are *local*, and that is not a detail. A pairwise chain can only see atom
pairs drawn from residues ``i-1``, ``i``, ``i+1``, so a contact at larger sequence
separation is outside the objective entirely -- which means the global optimum of that
objective can be, and was, a structure with physically impossible contacts. Rebuilding
dnmt3a's own deposited CA trace for residues 282-414 put the carbonyl oxygens of two
residues nine apart at **1.282 A**, shorter than an O=O double bond, against 4.827 A in the
deposited model, at 25 out of 25 rng seeds. Every bond length, every bond angle and every
local steric term was satisfied.

So the solve is iterated against a long-range field (:func:`_clash_field`): every unit's
atoms are swept over the whole rotation grid and scored against every already-placed atom
the hard-sphere limits apply to, and the chain is re-solved. That is a mean-field
approximation -- the partners are held where they currently are -- so exact single-unit
descent (:func:`_polish_clashes`) follows it, and the whole thing is measured
(:func:`_steric_survey`) rather than assumed: a reconstruction with a contact inside its
limit is refused, never returned. Which pairs the limits apply to is itself measured; see
:data:`_NONBONDED_RESIDUE_SEPARATION`.

Terminal residues
-----------------
Residue 0 has no ``CA(-1)`` so no unit supplies its ``N``, and residue ``n-1`` has no
``CA(n)`` so no unit supplies its ``C`` or ``O``. Neither is dropped. Both are placed
from internal coordinates with the one remaining dihedral taken from
:data:`_TERMINAL_PHI` / :data:`_TERMINAL_PSI`; see those for why those values -- except
when that dihedral would put the atom inside another atom, in which case it is rotated to
the nearest value that does not (:func:`_choose_terminal`). With ``allow_chain_breaks`` the
two segments either side of a break each have such a terminus pointing at the other, and
placing them from fixed dihedrals with no knowledge of each other put ``C(i)`` and
``N(i+1)`` **1.124 A** apart across dnmt3a's cis-prolines -- shorter than a C-N triple bond
-- which the PDB writer then emitted a CONECT record for.

What is *not* done
------------------
* No cis peptide bonds. omega is always trans, per
  :data:`~dodo.constants.OMEGA_TRANS`, so a CA-CA virtual bond near 2.9 A (cis-proline)
  is rejected rather than modelled. :mod:`dodo.constants` already states that DODO does
  not model cis-proline.
* No ``OXT`` on the C-terminus, and no amide or alpha hydrogens.
* With ``sidechains=True``, ``CB`` **and nothing beyond it**. That is a deliberate choice
  over shipping a fabricated rotamer library: CB alone already fixes rendering, and see
  :func:`place_sidechain_cb` for the three-angle derivation and the measurement that
  validates its chirality.

Two hard limits worth knowing about
-----------------------------------
**Angles.** ``alpha_C + alpha_N`` -- the angles ``C(i)`` and ``N(i+1)`` subtend at their own
alpha carbons, measured from the CA-CA axis -- is 35.50 degrees and is *invariant* to how
the peptide bond angles are strained (see :func:`ca_angle_budget`). Since those two angles
plus tau bracket the ``CA(i-1)-CA(i)-CA(i+1)`` pseudo-angle, no all-atom-reconstructable
trans backbone can have a CA pseudo-angle above ``N_CA_C_ANGLE + 35.50`` = 146.5 degrees, or
above 160.5 at the default ``N-CA-C`` tolerance
(:func:`max_reconstructable_ca_angle`). :data:`~dodo.constants.BACKBONE_ANGLE_MAX is 161.0,
which is ABOVE that ceiling: the window is the author's measured AF2 distribution, deliberately
not capped to what all-atom reconstruction can represent (CA-only correctness is the higher
priority).
:func:`unreconstructable_ca_angles` answers the per-residue question before anything is
built, which is what a generator needs; an empty answer is necessary, not sufficient.

The per-residue ceiling is not a per-*trace* one, because adjacent residues share a peptide
unit and cannot both draw most of its budget. Measured on self-avoiding random walks of 60
and 150 residues, six seeds each: with the CA-CA-CA window capped at 140 degrees, 12 of 12
reconstruct with no failure of any kind; at 145.5 -- a degree *under* the ideal ceiling --
still no ``N-CA-C`` failure, and the one refusal is steric; at 150 degrees, two of six
150-residue traces fail on ``N-CA-C``. The sampling window's actual upper bound is 161.0, well
above all of these, so a full-window trace fails more often still. That threshold moved out with
the tau tolerance: at the earlier :data:`_N_CA_C_TOLERANCE` of 8.0 it was 145.5 that failed
three of four.

**Density.** A CA trace also has to leave room for the atoms hanging off it. dnmt3a's
closest non-bonded CA-CA pair is 3.78 A, essentially the bonded distance, and it manages
that because its local backbone geometry cooperates; a random walk's does not.
Measured over walks of 60 residues, acceptance was 4 of 4 at a 5.0 A self-avoidance radius,
3 of 4 at 4.5 A and 2 of 4 at 3.8 A -- against
:data:`~dodo.constants.CA_CLASH_DISTANCE` of 3.20, which admits traces no backbone can be
built on. Refusals of this kind are steric and say so; they are not ``N-CA-C`` failures.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from ..constants import (
    C_N_CA_ANGLE,
    C_N_PEPTIDE_BOND_LENGTH,
    C_O_BOND_LENGTH,
    CA_C_BOND_LENGTH,
    CA_C_N_ANGLE,
    CA_C_O_ANGLE,
    N_CA_BOND_LENGTH,
    N_CA_C_ANGLE,
    OMEGA_TRANS,
    ONE_TO_THREE,
)
from ..exceptions import GeometryError
from ..structure import Domain, Structure

__all__ = [
    "BackboneResult",
    "StericClash",
    "ca_angle_budget",
    "max_reconstructable_ca_angle",
    "place_backbone",
    "place_backbone_for_domain",
    "place_sidechain_cb",
    "reconstructable_ca_bond_range",
    "unreconstructable_ca_angles",
]

# ---------------------------------------------------------------------------
# Local constants
#
# Nothing here duplicates constants.py. Each entry is either DERIVED from a constant
# declared there, or a value that constants.py does not yet carry -- in which case its
# provenance is recorded in full, because the whole point of constants.py is that a
# future reader can tell a measured value from a tuned knob.
# ---------------------------------------------------------------------------

#: The ``N-C-O`` angle at the carbonyl carbon for an unstrained peptide unit, in degrees.
#:
#: DERIVED. The carbonyl carbon is sp2, so ``CA``, ``N`` and ``O`` are coplanar with it
#: and the three angles at it sum to 360 degrees.
#:
#: Note which of the three is the derived one. When a unit has to be strained to close
#: onto a CA-CA distance, the strain lands on ``CA-C-N`` and something has to give, since
#: the three cannot all keep their declared values and still sum to 360. It gives here,
#: in the derived angle, so that both *declared* constants -- CA_C_N_ANGLE plus its
#: reported strain, and CA_C_O_ANGLE exactly -- are honoured. Letting ``CA-C-O`` absorb it
#: instead would silently put a declared constant off by the strain, which is the sort of
#: quiet drift constants.py exists to prevent.
_N_C_O_ANGLE: float = 360.0 - CA_C_N_ANGLE - CA_C_O_ANGLE

#: Largest peptide-angle strain, in degrees, used to close onto a CA-CA distance that is
#: not exactly 3.804 A.
#:
#: CHOICE, with a measured justification. A real CA trace does not have a uniform CA-CA
#: virtual bond, so something has to absorb the difference. Bond *lengths* are the
#: stiffest thing in a protein (Engh & Huber spreads are ~0.02 A) and bond *angles* are
#: an order of magnitude softer, so the discrepancy is absorbed entirely into
#: CA_C_N_ANGLE and C_N_CA_ANGLE, along the direction in that two-dimensional space that
#: changes the CA-CA distance fastest -- i.e. the minimum-norm perturbation for a given
#: change in distance. This constant caps the Euclidean norm of that perturbation.
#:
#: 12.0 degrees admits a CA-CA range of 3.659-3.932 A. Measured on the folded domains of
#: tests/data/structures/dnmt3a.pdb (an AlphaFold model), the trans CA-CA virtual bond
#: spans 3.769-3.913 A, so this window covers real input with margin while still
#: rejecting cis-peptides at ~2.9-3.0 A. Beyond the cap, place_backbone raises rather
#: than silently emitting a peptide bond nobody would recognise.
_MAX_PEPTIDE_ANGLE_STRAIN: float = 12.0

#: Number of samples in the peptide-plane rotation grid, i.e. 3-degree steps.
#:
#: CHOICE, and the honest justification is weaker than the one this comment used to
#: carry. It claimed grids of 90, 120, 180 and 360 agree on the reconstructed N/C/O
#: positions to within 0.005 A RMSD. Measured on dnmt3a 699-912, seed 20260730, the
#: all-atom RMSD against grid 120 is 0.111 A (grid 60), 0.069 A (90), 0.069 A (180) and
#: 0.053 A (360) -- an order of magnitude more than the claimed figure. So the objective
#: is *not* smooth at 3-degree resolution: it has near-degenerate minima that a shift of
#: the grid moves between, which is the same residual freedom the rng phase exposes on
#: purpose. What the grid does buy is bounded error -- the discretization can misplace a
#: peptide plane by at most half a step -- and cost is quadratic in it. 120 is where the
#: marginal accuracy stops paying for itself.
_PLANE_GRID: int = 120

#: Largest all-atom RMSD, in Angstroms, between reconstructions of the same trace at grids
#: of 90, 120, 180 and 360. MEASURED at 0.069 (see :data:`_PLANE_GRID`) and pinned just
#: above, so the documented number cannot drift away from the code again.
_PLANE_GRID_AGREEMENT_RMSD: float = 0.08

#: Scale, in degrees, on the ``N-CA-C`` deviation term of the objective.
#:
#: MEASURED, then rounded. The per-residue spread of ``N-CA-C`` in a real structure is
#: 1.7 degrees sd (dnmt3a folded domains, mean 112.2). This is *not* a hard bound: it is
#: the scale that makes a one-sigma tau deviation cost the same as a mild steric
#: overlap, so the two terms are commensurable.
_TAU_SD: float = 2.0

#: Weight on the hard-sphere steric term, per squared Angstrom of overlap.
#:
#: TUNED against tests/data/structures/dnmt3a.pdb, and the measurements are recorded here
#: because the tuning set is otherwise unrecoverable. Mean N/C/O RMSD over four folded
#: segments (700-800, 282-414, 600-698, 480-578):
#:
#:     weight 0   -> N 0.307  C 0.343  O 1.063
#:     weight 1   -> N 0.307  C 0.343  O 1.063
#:     weight 3   -> N 0.295  C 0.317  O 0.973
#:     weight 10  -> N 0.290  C 0.304  O 0.928
#:     weight 30  -> N 0.292  C 0.297  O 0.914
#:
#: 10 is where the curve flattens.
_STERIC_WEIGHT: float = 10.0

#: Weight on the positive-phi term of the objective.
#:
#: TUNED, and the single most valuable term after tau. L-amino acids overwhelmingly
#: disfavour phi > 0; the hard-sphere term captures part of that but not the whole of it,
#: because the clash that forbids positive phi involves the amide hydrogen this module
#: does not place. Measured on dnmt3a residues 700-800: without this term 21% of rebuilt
#: residues had phi > 0 against 6% in the reference, and O RMSD was 1.59 A; with it, 3%
#: and 1.12 A. The term is normalised so phi = +30 degrees costs 1.0, which is why the
#: weight itself is 1.
_POSITIVE_PHI_WEIGHT: float = 1.0

#: Normalisation for the positive-phi term, in squared degrees. DERIVED: 30**2, so a
#: 30-degree excursion into positive phi costs exactly one unit.
_POSITIVE_PHI_SCALE: float = 900.0

#: Tolerance on the reconstructed ``N-CA-C`` angle (tau), in degrees.
#:
#: MEASURED, and widened from an earlier 8.0 after re-measuring on a broader sample.
#:
#: The 8.0 value came from dnmt3a's FOLDED DOMAINS, where tau spans 108.6-116.1 degrees.
#: That is a real measurement but an unrepresentative one: folded cores hold tau tightly,
#: and the regions DODO rebuilds are loops and IDRs, which do not.
#:
#: Re-measured over 12,193 residues across four whole deposited structures (dnmt3a, p300,
#: arf19, 6kn7 -- folded domains, loops and termini alike): tau has mean 111.87, sd 3.36,
#: and spans 98.6-129.7 degrees. Fraction of those real residues a given gate would
#: REJECT as unreconstructable:
#:
#:     +/- 8.0 deg -> 4.19%      +/-12.0 deg -> 0.97%
#:     +/-10.0 deg -> 2.19%      +/-14.0 deg -> 0.25%
#:
#: Rejecting one real residue in 24 is not a defensible gate for a tool whose input is
#: real structures. 14.0 admits 99.75% of deposited geometry, and it is where the value
#: stops: 16.0 would admit tau below 95 degrees, which nothing in the sample exhibits.
#:
#: This is a correction of an over-strict gate, not a loosening to make output pass. The
#: effect is large -- generated-trace acceptance goes from 6/24 to 20/24 -- but the
#: justification is the deposited distribution, not the acceptance rate.
_N_CA_C_TOLERANCE: float = 14.0

#: Dihedral used to place the C-terminal residue's own ``C``, in degrees, as the
#: ``C(n-2)-N(n-1)-CA(n-1)-C(n-1)`` angle -- that is, phi of the last residue.
#:
#: CHOICE. The last residue's phi is unconstrained by the CA trace, so it has to be
#: chosen. -120 degrees is the centre of the broad beta / extended region, which is the
#: most populated part of the Ramachandran plot for a residue with no secondary
#: structure, and it is on the correct side for an L-amino acid. Anything in that region
#: would do; what matters is that the choice is stated rather than left to fall out of an
#: arbitrary cross product.
_TERMINAL_PHI: float = -120.0

#: Dihedral used to place the N-terminal residue's own ``N``, in degrees, as the
#: ``N(0)-CA(0)-C(0)-N(1)`` angle -- that is, psi of the first residue. Also used, offset
#: by 180 degrees, to place the last residue's carbonyl oxygen. CHOICE, for the same
#: reason as :data:`_TERMINAL_PHI`; +140 is the extended-beta partner of phi = -120.
_TERMINAL_PSI: float = 140.0

#: ``CA-CB`` bond length in Angstroms. MEASURED (Engh & Huber 1991).
_CA_CB_BOND_LENGTH: float = 1.530

#: ``N-CA-CB`` and ``C-CA-CB`` angles in degrees. MEASURED (Engh & Huber 1991). Together
#: with N_CA_C_ANGLE these three angles determine CB up to a reflection, and the
#: reflection is fixed by L-amino-acid chirality; see :func:`place_sidechain_cb`.
_N_CA_CB_ANGLE: float = 110.5
_C_CA_CB_ANGLE: float = 110.1

#: Minimum hard-sphere contact distances between backbone heavy atoms, in Angstroms,
#: keyed on the element pair.
#:
#: MEASURED, and citable: these are the classic "normally allowed" contact limits
#: Ramachandran, Ramakrishnan & Sasisekharan used in 1963 to derive the allowed phi/psi
#: regions from sterics alone. They are used here for exactly that purpose, which is why
#: this module needs no statistical Ramachandran library.
#:
#: These belong in constants.py; this module does not own that file, so they live here
#: with their provenance and the need is reported.
_MIN_CONTACT: dict[frozenset[str], float] = {
    frozenset(("C", "C")): 3.20,
    frozenset(("C", "N")): 2.90,
    frozenset(("C", "O")): 2.80,
    frozenset(("N", "N")): 2.70,
    frozenset(("N", "O")): 2.70,
    frozenset(("O", "O")): 2.80,
}

#: Residue separation at which a heavy-atom pair is treated as a hard non-bonded
#: contact, i.e. one that :data:`_MIN_CONTACT` applies to as a *limit* rather than as a
#: soft preference.
#:
#: MEASURED, and the measurement is the whole justification. Over all 912 residues of
#: tests/data/structures/dnmt3a.pdb, using N/CA/C/O/CB:
#:
#:     residue separation >= 2   0 of 507 pairs below their limit (closest margin +0.05 A)
#:     residue separation  = 1   463 of 1490 below (C(i)..C(i+1) reaches 2.92 A vs 3.20)
#:     residue separation  = 0   35 of 901 below (N(i)..O(i) reaches 2.60 A vs 2.70)
#:
#: So the 1963 "normally allowed" limits are a real, hard constraint two residues apart
#: and are *routinely* violated within one residue, where four-bond pairs across a planar
#: peptide unit have nowhere else to go. Applying them as hard limits locally would
#: reject every alpha helix in the PDB; applying them as hard limits at separation two or
#: more is what the deposited structure actually satisfies. Local pairs are still scored,
#: softly, by :func:`_pair_cost` -- that is what breaks the peptide-plane flip.
#:
#: Matches :data:`~dodo.constants.CLASH_EXCLUDE_WITHIN_RESIDUES` (2), which
#: :func:`~dodo.geometry.metrics.validate_ca_trace` uses for the same reason at the CA
#: level.
_NONBONDED_RESIDUE_SEPARATION: int = 2

#: Cost added per violated hard constraint on the ``N-CA-C`` angle, and per violated hard
#: non-bonded contact, in the Viterbi objective.
#:
#: DERIVED from what these numbers have to do, not tuned. Both are constraints, not
#: preferences, so they have to dominate every soft term (largest soft cost is on the
#: order of 200) rather than be traded against them. Before they existed the solver
#: routinely bought a 0.5 A steric improvement with a 10-degree ``N-CA-C`` error and then
#: the call raised: on random-walk traces whose pseudo-angles were all inside the
#: geometric ceiling, 23 of 62 residues came out beyond the tolerance when only 8 of them
#: were impossible.
#:
#: There are three tiers and the ordering between them is load-bearing, not cosmetic:
#:
#:     tau feasibility        1e9   exact, and decided by the trace: unrepairable
#:     local hard contact     1e6   exact, a function of both rotations in the pair term
#:     long-range contact     1e4   estimated against a *stale* neighbour position
#:
#: An estimate must never override an exact constraint, so the accumulated long-range
#: field is capped below the local hard cost (:data:`_MAX_CLASH_FIELD`). Without that cap
#: the observed failure was exactly the inversion it prevents: enough accumulated
#: long-range penalty on one unit made the solver buy relief by violating an
#: exactly-scored two-residue contact instead, turning ``O(22)..N(24)`` into 2.279 A
#: against a 2.70 A limit while chasing a CB-CB overlap thirty residues away.
_INFEASIBLE_TAU_COST: float = 1.0e9
_INFEASIBLE_CLASH_COST: float = 1.0e6

#: Cost added per long-range contact event by :func:`_clash_field`, for one rotation of one
#: peptide unit.
_CLASH_FIELD_COST: float = 1.0e4

#: Where an atom that does not exist is parked so a radius query cannot find it. Any
#: coordinate far outside a protein does; this one is obviously not a real position.
_FAR_AWAY: float = 1.0e6

#: Cap on any single accumulated long-range clash penalty. Keeps the mean-field term
#: dominant over every soft cost (largest is on the order of 200) and subordinate to both
#: exactly-evaluated constraints.
_MAX_CLASH_FIELD: float = 1.0e5

#: How many times the solver may re-solve with clash penalties added before giving up and
#: raising. CHOICE. Measured on dnmt3a's four folded domains and on self-avoiding
#: random-walk traces, every repairable clash is gone within three passes; the extra
#: passes cost nothing when they are not needed because the loop exits as soon as the
#: trace is clean.
_MAX_CLASH_REPAIR_PASSES: int = 8

#: How many exact descent sweeps :func:`_polish_clashes` may take per repair pass.
#: CHOICE: each sweep either strictly improves the model or ends the loop, so this only
#: bounds the work, not the correctness.
_MAX_POLISH_SWEEPS: int = 4


#: Number of samples used when choosing a terminal residue's one free dihedral, i.e.
#: 5-degree steps. CHOICE: the terminal atoms have no neighbour to close onto, so the
#: only thing this resolution has to be fine enough for is clash avoidance.
_TERMINAL_DIHEDRAL_GRID: int = 72

#: Non-bonded pairs scored by the steric term, as ``(name, name)`` keys into the atom
#: table built per residue in :func:`_pair_cost`.
#:
#: Every pair listed is separated by at least three bonds, so none of them is fixed by a
#: bond length or bond angle -- they vary with phi or psi, which is what makes them
#: informative. ``Cp``/``Op`` are residue ``i-1``'s carbon and oxygen, ``Nn``/``CAn`` are
#: residue ``i+1``'s nitrogen and alpha carbon, unprefixed names belong to residue ``i``.
#: 1-2 and 1-3 pairs are deliberately absent: they are constants of the geometry and
#: scoring them would add a constant to every candidate's cost, which is at best wasted
#: work and at worst a bias.
#:
#: The list is complete with respect to :func:`_steric_violations`: every pair of atoms
#: drawn from residues ``i-1``, ``i``, ``i+1`` that this cost function can see and whose
#: separation is not fixed by the geometry appears here. That completeness is what lets
#: :func:`_repair_clashes` skip these pairs when it builds its long-range penalties --
#: they are already handled here, exactly, as a function of both rotations rather than
#: against a stale neighbour.
_STERIC_PAIRS: tuple[tuple[str, str], ...] = (
    # Determined by phi: what rotating about N(i)-CA(i) moves past residue i-1.
    ("CAp", "C"),
    ("CAp", "O"),
    ("CAp", "CB"),
    ("Cp", "C"),
    ("Cp", "O"),
    ("Cp", "CB"),
    ("Op", "C"),
    ("Op", "O"),
    ("Op", "CB"),
    # Determined by psi: what rotating about CA(i)-C(i) moves past residue i+1.
    ("N", "O"),
    ("N", "Nn"),
    ("N", "CAn"),
    ("CB", "O"),
    ("CB", "Nn"),
    ("CB", "CAn"),
    # Cross terms, sensitive to both. These are the two-residue-separation pairs, so they
    # are the ones _MIN_CONTACT bounds absolutely rather than softly; see
    # _NONBONDED_RESIDUE_SEPARATION. CAp-CAn is absent because both atoms are input
    # alpha carbons: their separation is a property of the trace, not of any rotation.
    ("CAp", "Nn"),
    ("Cp", "Nn"),
    ("Op", "Nn"),
    ("Cp", "CAn"),
    ("Op", "CAn"),
)

#: Element of each name used in :data:`_STERIC_PAIRS`.
_STERIC_ELEMENTS: dict[str, str] = {
    "CAp": "C",
    "Cp": "C",
    "Op": "O",
    "N": "N",
    "CB": "C",
    "C": "C",
    "O": "O",
    "Nn": "N",
    "CAn": "C",
}

#: Residue offset, relative to residue ``i``, of each name used in :data:`_STERIC_PAIRS`.
_STERIC_OFFSETS: dict[str, int] = {
    "CAp": -1,
    "Cp": -1,
    "Op": -1,
    "N": 0,
    "CB": 0,
    "C": 0,
    "O": 0,
    "Nn": 1,
    "CAn": 1,
}

#: Residues for which no ``CB`` is placed. Glycine has none; ``UNK`` stands for a residue
#: whose identity is not known, and inventing a side chain for it would be inventing
#: information.
_NO_CB_RESIDUES: frozenset[str] = frozenset(("GLY", "UNK"))

#: Atom names emitted per residue, in the canonical within-residue order. CB is appended
#: when side chains are requested and the residue has one.
_BACKBONE_ORDER: tuple[str, ...] = ("N", "CA", "C", "O")

#: Element of each atom this module places.
_ELEMENT_OF: dict[str, str] = {"N": "N", "CA": "C", "C": "C", "O": "O", "CB": "C"}

#: Which entries of :data:`_STERIC_PAIRS` are hard limits. DERIVED: exactly those whose
#: two atoms are :data:`_NONBONDED_RESIDUE_SEPARATION` residues apart or more.
_STERIC_PAIRS_HARD: tuple[bool, ...] = tuple(
    abs(_STERIC_OFFSETS[right] - _STERIC_OFFSETS[left]) >= _NONBONDED_RESIDUE_SEPARATION
    for left, right in _STERIC_PAIRS
)

#: The pairs :data:`_STERIC_PAIRS` already scores, keyed as
#: ``(residue separation, name of earlier residue's atom, name of later residue's atom)``.
#: DERIVED from :data:`_STERIC_PAIRS` and :data:`_STERIC_OFFSETS`.
#:
#: :func:`_repair_clashes` skips these when it penalizes long-range contacts, because
#: :func:`_pair_cost` already scores them as an exact function of both peptide-plane
#: rotations. Penalizing them again from a stale neighbour position would let a stale
#: measurement veto a rotation that is fine in the joint solution.
#: Real atom name behind each :data:`_STERIC_PAIRS` name.
_STERIC_BASE_NAME: dict[str, str] = {
    "CAp": "CA",
    "Cp": "C",
    "Op": "O",
    "N": "N",
    "CB": "CB",
    "C": "C",
    "O": "O",
    "Nn": "N",
    "CAn": "CA",
}


def _locally_scored_keys() -> frozenset[tuple[int, str, str]]:
    """Build :data:`_LOCALLY_SCORED`. Order within a pair follows residue order."""
    keys: set[tuple[int, str, str]] = set()
    for left, right in _STERIC_PAIRS:
        low, high = left, right
        if _STERIC_OFFSETS[left] > _STERIC_OFFSETS[right] or (
            _STERIC_OFFSETS[left] == _STERIC_OFFSETS[right]
            and _STERIC_BASE_NAME[left] > _STERIC_BASE_NAME[right]
        ):
            low, high = right, left
        separation = abs(_STERIC_OFFSETS[right] - _STERIC_OFFSETS[left])
        keys.add((separation, _STERIC_BASE_NAME[low], _STERIC_BASE_NAME[high]))
    return frozenset(keys)


_LOCALLY_SCORED: frozenset[tuple[int, str, str]] = _locally_scored_keys()

#: Largest entry in :data:`_MIN_CONTACT`; the search radius for any contact query.
_MAX_CONTACT: float = max(_MIN_CONTACT.values())


# ---------------------------------------------------------------------------
# Internal-coordinate placement
# ---------------------------------------------------------------------------


def _place_atom(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    bond: float,
    angle: float,
    dihedral: float,
) -> np.ndarray:
    """Place atom ``d`` from internal coordinates relative to ``a``, ``b``, ``c``.

    The standard "natural extension reference frame" construction: ``d`` is placed at
    distance ``bond`` from ``c``, making angle ``angle`` at ``c`` with ``b``, and
    dihedral ``dihedral`` about the ``b-c`` bond measured from ``a``.

    Broadcasts over leading axes, so ``a``, ``b`` and ``c`` may be ``(3,)`` or any
    ``(..., 3)`` that broadcast together.

    Parameters
    ----------
    a, b, c
        The three reference positions, ``(..., 3)``.
    bond
        ``c-d`` distance in Angstroms.
    angle
        ``b-c-d`` angle in degrees.
    dihedral
        ``a-b-c-d`` dihedral in degrees.

    Returns
    -------
    np.ndarray
        Position of ``d``, broadcast to the common shape.

    Raises
    ------
    GeometryError
        If ``a``, ``b`` and ``c`` are collinear, in which case the dihedral reference
        plane is undefined and any answer would be arbitrary. The pre-rewrite code took
        a cross product here without checking and propagated the resulting NaN into
        written coordinates.
    """
    theta = np.radians(angle)
    phi = np.radians(dihedral)

    bc = c - b
    bc = bc / np.linalg.norm(bc, axis=-1, keepdims=True)
    normal = np.cross(b - a, bc)
    normal_length = np.linalg.norm(normal, axis=-1, keepdims=True)
    if bool(np.any(normal_length < 1e-8)):
        raise GeometryError(
            "Cannot place an atom from internal coordinates: the three reference atoms "
            "are collinear, so the dihedral has no reference plane. This means the CA "
            "trace contains three consecutive alpha carbons in a straight line."
        )
    normal = normal / normal_length
    in_plane = np.cross(normal, bc)

    placed: np.ndarray = (
        c
        + (-bond * np.cos(theta)) * bc
        + (bond * np.sin(theta) * np.cos(phi)) * in_plane
        + (bond * np.sin(theta) * np.sin(phi)) * normal
    )
    return placed


def _place_atom_over_dihedrals(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    bond: float,
    angle: float,
    dihedrals: np.ndarray,
) -> np.ndarray:
    """:func:`_place_atom` for many dihedral values at once, returning ``(k, 3)``.

    Same construction, one vectorized pass. Terminal placement is searched over a whole
    circle of dihedrals on every repair pass, and doing that as a Python loop over
    :func:`_place_atom` dominated the runtime of a rebuild.
    """
    theta = np.radians(angle)
    phi = np.radians(np.asarray(dihedrals, dtype=np.float64))

    bc = c - b
    bc = bc / np.linalg.norm(bc)
    normal = np.cross(b - a, bc)
    length = float(np.linalg.norm(normal))
    if length < 1e-8:
        raise GeometryError(
            "Cannot place an atom from internal coordinates: the three reference atoms "
            "are collinear, so the dihedral has no reference plane. This means the CA "
            "trace contains three consecutive alpha carbons in a straight line."
        )
    normal = normal / length
    in_plane = np.cross(normal, bc)
    placed: np.ndarray = (
        c[None, :]
        + (-bond * np.cos(theta)) * bc[None, :]
        + (bond * np.sin(theta) * np.cos(phi))[:, None] * in_plane[None, :]
        + (bond * np.sin(theta) * np.sin(phi))[:, None] * normal[None, :]
    )
    return placed


def _dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> np.ndarray:
    """Dihedral angle ``p0-p1-p2-p3`` in degrees, broadcasting over leading axes."""
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    axis = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * axis).sum(-1, keepdims=True) * axis
    w = b2 - (b2 * axis).sum(-1, keepdims=True) * axis
    result: np.ndarray = np.degrees(np.arctan2((np.cross(axis, v) * w).sum(-1), (v * w).sum(-1)))
    return result


def _angle_between(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Angle in degrees between two vector fields, ``(..., 3)``."""
    cosine = (u * v).sum(-1) / (np.linalg.norm(u, axis=-1) * np.linalg.norm(v, axis=-1))
    result: np.ndarray = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    return result


# ---------------------------------------------------------------------------
# The rigid peptide unit
# ---------------------------------------------------------------------------


@lru_cache(maxsize=1)
def _strain_direction() -> tuple[float, float]:
    """Return the direction in bond-angle space that changes CA-CA fastest.

    This is the normalised gradient of the CA-CA virtual bond length with respect to the
    two peptide bond angles, in degrees. Strain is applied along this direction because
    it is, to first order, the *minimum-norm* angular perturbation that achieves a given
    change in the CA-CA distance -- i.e. it distorts the peptide unit as little as
    possible for the closure it has to achieve. Measured value is (0.817, 0.577): the
    ``CA-C-N`` angle does about 40% more of the work than ``C-N-CA``.
    """
    step = 1e-5
    base = _unit_ca_distance(CA_C_N_ANGLE, C_N_CA_ANGLE)
    gradient = np.array(
        [
            (_unit_ca_distance(CA_C_N_ANGLE + step, C_N_CA_ANGLE) - base) / step,
            (_unit_ca_distance(CA_C_N_ANGLE, C_N_CA_ANGLE + step) - base) / step,
        ]
    )
    gradient /= np.linalg.norm(gradient)
    return float(gradient[0]), float(gradient[1])


def _unit_ca_distance(ca_c_n: float, c_n_ca: float) -> float:
    """CA(i)-CA(i+1) distance for one trans peptide unit with the given bond angles."""
    coords = _peptide_unit(np.array([ca_c_n]), np.array([c_n_ca]))
    return float(np.linalg.norm(coords["CA2"][0]))


def _peptide_unit(ca_c_n: np.ndarray, c_n_ca: np.ndarray) -> dict[str, np.ndarray]:
    """Build ``n`` trans peptide units in their own plane, one per pair of bond angles.

    ``CA1`` is at the origin and ``C`` on the +x axis; everything lies in ``z = 0``,
    which is what trans omega means for these five atoms. Returns the raw positions,
    still in the construction frame -- :func:`_peptide_unit_local` reframes them onto the
    CA-CA axis.
    """
    n = ca_c_n.shape[0]
    ca1 = np.zeros((n, 3))
    carbon = np.zeros((n, 3))
    carbon[:, 0] = CA_C_BOND_LENGTH

    radians = np.radians(ca_c_n)
    nitrogen = carbon.copy()
    nitrogen[:, 0] += C_N_PEPTIDE_BOND_LENGTH * -np.cos(radians)
    nitrogen[:, 1] += C_N_PEPTIDE_BOND_LENGTH * np.sin(radians)

    # omega trans puts CA(i) and CA(i+1) on opposite sides of the C-N bond, which for a
    # planar unit is exactly a dihedral of 180 degrees.
    ca2 = _place_atom_varying(ca1, carbon, nitrogen, N_CA_BOND_LENGTH, c_n_ca, OMEGA_TRANS)
    # The O is cis to CA(i+1) about C-N, i.e. on the same side as CA(i+1) and opposite
    # CA(i). Placing it this way -- rather than declaring an angle from CA -- is what
    # makes the resulting CA-C-O angle come out at exactly CA_C_O_ANGLE, since the three
    # angles at the sp2 carbon then close to 360 by construction. Note the N-C-O angle is
    # computed from this unit's *actual* CA-C-N angle rather than from the unstrained
    # constant: with the strained value it is CA_C_O_ANGLE that stays exact, which is the
    # constant that was declared.
    oxygen = _place_atom_varying(
        ca2, nitrogen, carbon, C_O_BOND_LENGTH, 360.0 - ca_c_n - CA_C_O_ANGLE, 0.0
    )
    return {"CA1": ca1, "C": carbon, "N": nitrogen, "CA2": ca2, "O": oxygen}


def _place_atom_varying(
    a: np.ndarray,
    b: np.ndarray,
    c: np.ndarray,
    bond: float,
    angle: np.ndarray,
    dihedral: float,
) -> np.ndarray:
    """:func:`_place_atom` with a per-row bond angle. Same construction, vectorized."""
    theta = np.radians(angle)[:, None]
    phi = np.radians(dihedral)

    bc = c - b
    bc = bc / np.linalg.norm(bc, axis=-1, keepdims=True)
    normal = np.cross(b - a, bc)
    normal = normal / np.linalg.norm(normal, axis=-1, keepdims=True)
    in_plane = np.cross(normal, bc)
    placed: np.ndarray = (
        c
        + (-bond * np.cos(theta)) * bc
        + (bond * np.sin(theta) * np.cos(phi)) * in_plane
        + (bond * np.sin(theta) * np.sin(phi)) * normal
    )
    return placed


def _peptide_unit_local(strain: np.ndarray) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """In-plane coordinates of ``C``, ``N`` and ``O`` on the CA-CA axis.

    Parameters
    ----------
    strain
        Per-unit angular strain in degrees, applied along :func:`_strain_direction`.

    Returns
    -------
    tuple
        The CA-CA distance of each unit, ``(n,)``, and a mapping from atom name to
        ``(n, 2)`` coordinates in the frame whose +x axis points from ``CA(i)`` to
        ``CA(i+1)`` and whose +y axis is the in-plane perpendicular. Two coordinates,
        not three: every atom of a trans peptide unit is in that plane, and carrying a
        third that is identically zero invites someone to put something in it.
    """
    g1, g2 = _strain_direction()
    coords = _peptide_unit(CA_C_N_ANGLE + strain * g1, C_N_CA_ANGLE + strain * g2)

    distance = np.linalg.norm(coords["CA2"], axis=-1)
    x_axis = coords["CA2"] / distance[:, None]
    # The unit is planar in z = 0, so the in-plane perpendicular is just x rotated by 90
    # degrees about z. Deriving it this way keeps the y sign consistent across units,
    # which matters only for readability -- the rotation angle spans the full circle.
    y_axis = np.stack([-x_axis[:, 1], x_axis[:, 0], np.zeros_like(x_axis[:, 0])], axis=-1)

    local: dict[str, np.ndarray] = {}
    for name in ("C", "N", "O"):
        point = coords[name]
        local[name] = np.stack([(point * x_axis).sum(-1), (point * y_axis).sum(-1)], axis=-1)
    return distance, local


@lru_cache(maxsize=1)
def _strain_table() -> tuple[np.ndarray, np.ndarray]:
    """Monotone table mapping angular strain to CA-CA distance, for inversion.

    Built lazily rather than at import: :mod:`dodo.constants` is imported by the CLI on
    every invocation and a few thousand small numpy calls is not something ``dodo --help``
    should pay for.
    """
    strain = np.linspace(-_MAX_PEPTIDE_ANGLE_STRAIN, _MAX_PEPTIDE_ANGLE_STRAIN, 4001)
    distance, _ = _peptide_unit_local(strain)
    return strain, distance


def reconstructable_ca_bond_range() -> tuple[float, float]:
    """CA-CA virtual bond lengths this module can close onto, in Angstroms.

    A trans peptide unit built from the declared bond lengths and *ideal* bond angles has
    a fixed CA-CA distance of 3.804 A. Anything else has to be absorbed by straining the
    two peptide bond angles, and :data:`_MAX_PEPTIDE_ANGLE_STRAIN` caps how far that
    goes. Measured window: 3.659 to 3.932 A.

    A cis peptide bond puts CA-CA near 2.9 A and is outside this window by design.
    """
    _, distance = _strain_table()
    return float(distance[0]), float(distance[-1])


def ca_angle_budget() -> float:
    """Sum of the angles ``C(i)`` and ``N(i+1)`` subtend at their own alpha carbons.

    Precisely: ``angle(CA(i+1), CA(i), C(i)) + angle(CA(i), CA(i+1), N(i+1))``, in
    degrees, for one trans peptide unit at zero strain.

    This number is the reason :func:`max_reconstructable_ca_angle` exists. Measured value
    is 35.50 degrees, and -- this is the useful part -- it is essentially *invariant* to
    strain. Along the one-parameter family of ``(CA-C-N, C-N-CA)`` pairs that hold the
    CA-CA distance fixed at 3.804 A, the two angles trade off against each other while
    their sum barely moves::

        CA-C-N   C-N-CA   alpha_C   alpha_N   sum
        116.20   121.70     20.50     15.00   35.50
        120.20   117.22     16.38     19.36   35.74
        124.20   114.28     12.88     22.80   35.68
        128.20   112.28      9.75     25.67   35.43
        132.20   110.94      6.88     28.14   35.02

    Sixteen degrees of trade-off moves alpha_C by 13.6 degrees and the sum by 0.5. So a
    residue cannot buy angular headroom on both of its sides at once, and two adjacent
    residues cannot both draw on the same peptide unit's budget -- which is why a trace
    with two consecutive wide CA pseudo-angles is unreconstructable even when either one
    alone would be fine.

    For comparison, the same two angles measured on dnmt3a's folded domains average
    20.62 and 14.79 degrees, reaching 24.1 and 16.5 where the pseudo-angle demands it.
    """
    distance, local = _peptide_unit_local(np.array([0.0]))
    alpha_c = float(np.degrees(np.arctan2(-local["C"][0, 1], local["C"][0, 0])))
    alpha_n = float(np.degrees(np.arctan2(local["N"][0, 1], distance[0] - local["N"][0, 0])))
    return alpha_c + alpha_n


def max_reconstructable_ca_angle(*, n_ca_c_tolerance: float = 0.0) -> float:
    """Largest ``CA(i-1)-CA(i)-CA(i+1)`` pseudo-angle that admits a valid backbone.

    With ``n_ca_c_tolerance`` at its default of zero this is the ceiling for an *ideal*
    ``N-CA-C``. A caller who accepts a deviation of ``t`` degrees gets ``t`` more degrees
    of ceiling, because the bound below is on ``|pseudo-angle - N-CA-C|``, so the number a
    trace generator needs is this function evaluated at the tolerance it will build with:
    146.5 degrees at zero, and wider as the tolerance grows (see _N_CA_C_TOLERANCE).

    ``N(i)`` lies on a cone of half-angle ``alpha_N`` about the ``CA(i) -> CA(i-1)``
    axis and ``C(i)`` on a cone of half-angle ``alpha_C`` about ``CA(i) -> CA(i+1)``,
    so the angle between them -- which is ``N-CA-C`` -- can be no further from the angle
    between those two axes than ``alpha_N + alpha_C``. The angle between the axes *is*
    the CA pseudo-angle. Hence the pseudo-angle cannot exceed
    ``N_CA_C_ANGLE + ca_angle_budget()``, and :func:`ca_angle_budget` shows that budget
    cannot be enlarged by straining the peptide angles.

    Measured value: 146.5 degrees, and 160.5 at :data:`_N_CA_C_TOLERANCE`, against a
    :data:`~dodo.constants.BACKBONE_ANGLE_MAX` of 161.0. So the sampling window extends PAST
    this ceiling at its wide end -- deliberately, because the window is the author's measured
    AF2 distribution and getting CA-only geometry right outranks all-atom reconstructability.
    Two consequences follow, and both are real: a few angles from the top of the window are
    individually unreconstructable, and even angles inside the ceiling compound, because
    adjacent residues that each need most of the budget share one peptide unit and cannot both
    have it. The second is the dominant effect, which is why capping the window's magnitude
    never fixed acceptance -- the constraint that matters is on the CHANGE in pseudo-angle
    between consecutive residues.

    Raises
    ------
    GeometryError
        If ``n_ca_c_tolerance`` is negative or not finite. An infinite tolerance would
        return an infinite ceiling, which is the same silent pass
        :func:`_validate_tolerance` refuses.
    """
    tolerance = float(n_ca_c_tolerance)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise GeometryError(
            f"n_ca_c_tolerance must be finite and non-negative, got {n_ca_c_tolerance}."
        )
    return N_CA_C_ANGLE + ca_angle_budget() + tolerance


def unreconstructable_ca_angles(
    ca_coords: np.ndarray, *, n_ca_c_tolerance: float = _N_CA_C_TOLERANCE
) -> np.ndarray:
    """Residue indices whose CA pseudo-angle no trans backbone can realize.

    A pre-flight check, so that a generator or a caller can find out *before* building
    which residues are impossible and why -- rather than being told afterwards that the
    whole trace was refused, or being handed a structure with a 130-degree ``N-CA-C``
    angle in it.

    The answer is exact, not a heuristic: a pseudo-angle above
    ``max_reconstructable_ca_angle(n_ca_c_tolerance=...)`` forces ``N-CA-C`` outside the
    tolerance no matter what the rest of the trace does (see
    :func:`max_reconstructable_ca_angle`). The converse is not true residue by residue --
    two adjacent residues can each be inside the ceiling and still compete for one peptide
    unit's angular budget -- so an empty return is a necessary rather than a sufficient
    condition for :func:`place_backbone` to succeed.

    Parameters
    ----------
    ca_coords
        ``(n_residues, 3)`` alpha-carbon coordinates.
    n_ca_c_tolerance
        The tolerance the caller intends to build with, in degrees.

    Returns
    -------
    np.ndarray
        0-based residue indices, ascending. Both termini are never included: they have no
        pseudo-angle, and their ``N-CA-C`` is placed exactly.
    """
    ca = np.ascontiguousarray(np.asarray(ca_coords, dtype=np.float64))
    if ca.ndim != 2 or ca.shape[1] != 3:
        raise GeometryError(f"ca_coords must have shape (n_residues, 3), got {ca.shape}.")
    ceiling = max_reconstructable_ca_angle(n_ca_c_tolerance=n_ca_c_tolerance)
    pseudo = _ca_pseudo_angles(ca)
    # NaN at the termini compares False, which is the intended answer rather than an
    # accident: nothing constrains a terminal residue's pseudo-angle because it has none.
    return np.flatnonzero(pseudo > ceiling)


# ---------------------------------------------------------------------------
# CB
# ---------------------------------------------------------------------------


def place_sidechain_cb(n_xyz: np.ndarray, ca_xyz: np.ndarray, c_xyz: np.ndarray) -> np.ndarray:
    """Place ``CB`` from a residue's own ``N``, ``CA`` and ``C``.

    Derivation, rather than the hardcoded magic vector that circulates for this. Let
    ``nhat`` and ``chat`` be the unit vectors from ``CA`` to ``N`` and to ``C``. The unit
    vector ``w`` from ``CA`` to ``CB`` satisfies two linear equations,
    ``w . nhat = cos(110.5)`` and ``w . chat = cos(110.1)``, from the two measured bond
    angles. Those fix ``w``'s components in the ``N-CA-C`` plane; ``|w| = 1`` then fixes
    the magnitude of its out-of-plane component, and only the *sign* is left. That sign
    is the residue's chirality, and L-amino acids take the positive side of
    ``nhat x chat``.

    That chirality choice is not asserted, it is measured: applying this construction to
    the reference ``N``, ``CA`` and ``C`` of tests/data/structures/dnmt3a.pdb reproduces
    the deposited ``CB`` to **0.078 A RMSD** with the positive sign and 2.392 A with the
    negative one. v1's only side-chain consumer applied templates by pure translation
    with no rotation at all, so every side chain pointed the same lab direction; nothing
    about chirality was even expressible.

    Parameters
    ----------
    n_xyz, ca_xyz, c_xyz
        Backbone positions, ``(..., 3)``, broadcasting together.

    Returns
    -------
    np.ndarray
        ``CB`` positions, broadcast to the common shape.

    Raises
    ------
    GeometryError
        If ``N``, ``CA`` and ``C`` are collinear, so the plane is undefined.
    """
    n_dir = n_xyz - ca_xyz
    n_dir = n_dir / np.linalg.norm(n_dir, axis=-1, keepdims=True)
    c_dir = c_xyz - ca_xyz
    c_dir = c_dir / np.linalg.norm(c_dir, axis=-1, keepdims=True)

    overlap = (n_dir * c_dir).sum(-1, keepdims=True)
    determinant = 1.0 - overlap**2
    if bool(np.any(determinant < 1e-8)):
        raise GeometryError(
            "Cannot place CB: this residue's N, CA and C are collinear, so the "
            "backbone plane that fixes CB's position does not exist."
        )
    cos_n = np.cos(np.radians(_N_CA_CB_ANGLE))
    cos_c = np.cos(np.radians(_C_CA_CB_ANGLE))
    coeff_n = (cos_n - cos_c * overlap) / determinant
    coeff_c = (cos_c - cos_n * overlap) / determinant
    in_plane = coeff_n * n_dir + coeff_c * c_dir

    normal = np.cross(n_dir, c_dir)
    normal = normal / np.linalg.norm(normal, axis=-1, keepdims=True)
    # Clipped at zero because the three angles are mutually consistent only to rounding;
    # a tiny negative radicand here would otherwise become NaN in the output file.
    out_of_plane = np.sqrt(np.clip(1.0 - (in_plane * in_plane).sum(-1, keepdims=True), 0.0, None))

    placed: np.ndarray = ca_xyz + _CA_CB_BOND_LENGTH * (in_plane + out_of_plane * normal)
    return placed


# ---------------------------------------------------------------------------
# The solver
# ---------------------------------------------------------------------------


def _local_frames(ca: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-peptide-unit orthonormal frame on the CA-CA axis.

    Returns ``(distance, x_axis, e1, e2)``, each with one row per unit. ``e1`` is built
    from whichever world axis is least aligned with ``x_axis``, which keeps the cross
    product well-conditioned; the specific choice is immaterial because the rotation
    parameter it seeds spans the whole circle.
    """
    delta = np.diff(ca, axis=0)
    distance = np.linalg.norm(delta, axis=1)
    if bool(np.any(distance < 1e-8)):
        duplicated = int(np.argmin(distance))
        raise GeometryError(
            f"Alpha carbons {duplicated} and {duplicated + 1} are at the same position, "
            f"so there is no CA-CA axis to build a peptide unit on."
        )
    x_axis = delta / distance[:, None]

    reference = np.zeros_like(x_axis)
    reference[np.arange(x_axis.shape[0]), np.argmin(np.abs(x_axis), axis=1)] = 1.0
    e1 = reference - np.einsum("ij,ij->i", reference, x_axis)[:, None] * x_axis
    e1 /= np.linalg.norm(e1, axis=1)[:, None]
    e2 = np.cross(x_axis, e1)
    return distance, x_axis, e1, e2


def _contact_limit(left: str, right: str) -> float:
    """Hard-sphere minimum separation for two atom names, in Angstroms."""
    return _MIN_CONTACT[frozenset((_ELEMENT_OF[left], _ELEMENT_OF[right]))]


def _pair_key(separation: int, earlier: str, later: str) -> tuple[int, str, str]:
    """Canonical key into :data:`_LOCALLY_SCORED` for one atom pair."""
    if separation == 0 and earlier > later:
        return separation, later, earlier
    return separation, earlier, later


def _pair_cost(
    ca: np.ndarray,
    residue: int,
    c_vec: np.ndarray,
    n_vec: np.ndarray,
    o_vec: np.ndarray,
    has_cb: bool,
    n_ca_c_tolerance: float,
) -> np.ndarray:
    """Cost of every ``(theta_{i-1}, theta_i)`` pair at one interior residue.

    Four terms, and each is here for a stated reason:

    1. ``N-CA-C`` deviation, as a hard constraint plus a soft cost. The tolerance is the
       caller's acceptance criterion, so the optimizer has to treat it as a *constraint*:
       a purely soft version lets the search buy a fraction of an Angstrom of steric
       relief with a 10-degree tau error and then :func:`place_backbone` refuses the
       whole trace. Measured on 62-residue random-walk traces, that trade cost 23 of 62
       residues their tolerance when only 8 were geometrically impossible.
    2. Hard-sphere overlap between backbone atoms three or more bonds apart. This is the
       Ramachandran criterion in its original form and it is what breaks the
       peptide-plane flip degeneracy that term 1 alone leaves wide open. Hard for the
       pairs that are two residues apart, where the 1963 limits are a real bound, soft
       within one residue, where real structures sit below them; see
       :data:`_NONBONDED_RESIDUE_SEPARATION`.
    3. A penalty on positive phi. L-amino acids disfavour it, and the specific clash that
       forbids it involves the amide hydrogen, which this module does not place -- so
       term 2 cannot see it and it has to be stated separately.

    Returns an ``(grid, grid)`` array indexed by the previous unit's rotation then this
    unit's.
    """
    grid = c_vec.shape[1]
    previous = residue - 1

    # Both operands are unit-length by construction, so the dot product is the cosine
    # directly -- but clip anyway: an unclipped value of 1 + 1e-16 makes arccos NaN, and
    # a NaN here would win the argmin and be written to a file.
    n_unit = n_vec[previous] / N_CA_BOND_LENGTH
    c_unit = c_vec[residue] / CA_C_BOND_LENGTH
    tau = np.degrees(np.arccos(np.clip(n_unit @ c_unit.T, -1.0, 1.0)))
    deviation = np.abs(tau - N_CA_C_ANGLE)
    cost: np.ndarray = (deviation / _TAU_SD) ** 2 + _INFEASIBLE_TAU_COST * (
        deviation > n_ca_c_tolerance
    )

    # Shaped so that broadcasting does the outer product: atoms fixed by the previous
    # unit vary along axis 0, atoms fixed by this unit along axis 1.
    atoms: dict[str, np.ndarray] = {
        "CAp": ca[previous].reshape(1, 1, 3),
        "CAn": ca[residue + 1].reshape(1, 1, 3),
        "Cp": (ca[previous] + c_vec[previous])[:, None, :],
        "Op": (ca[previous] + o_vec[previous])[:, None, :],
        "N": (ca[residue] + n_vec[previous])[:, None, :],
        "C": (ca[residue] + c_vec[residue])[None, :, :],
        "O": (ca[residue] + o_vec[residue])[None, :, :],
        "Nn": (ca[residue + 1] + n_vec[residue])[None, :, :],
    }
    centre = ca[residue].reshape(1, 1, 3)
    if has_cb:
        atoms["CB"] = place_sidechain_cb(atoms["N"], centre, atoms["C"])

    overlap = np.zeros((grid, grid))
    breaches = np.zeros((grid, grid))
    for (left, right), hard in zip(_STERIC_PAIRS, _STERIC_PAIRS_HARD, strict=True):
        if not has_cb and "CB" in (left, right):
            continue
        limit = _MIN_CONTACT[frozenset((_STERIC_ELEMENTS[left], _STERIC_ELEMENTS[right]))]
        separation = np.linalg.norm(atoms[left] - atoms[right], axis=-1)
        excess = np.clip(limit - separation, 0.0, None)
        overlap += excess**2
        if hard:
            breaches += excess > 0.0
    cost = cost + _STERIC_WEIGHT * overlap + _INFEASIBLE_CLASH_COST * breaches

    phi = _dihedral(atoms["Cp"], atoms["N"], centre, atoms["C"])
    cost = cost + _POSITIVE_PHI_WEIGHT * np.clip(phi, 0.0, None) ** 2 / _POSITIVE_PHI_SCALE
    return cost


def _fixed_arm_tau_cost(
    fixed_arm: np.ndarray, moving_arm: np.ndarray, n_ca_c_tolerance: float
) -> np.ndarray:
    """``N-CA-C`` cost over the grid when one of the two arms is a fixed atom.

    Used at a segment end whose ``N`` or ``C`` is supplied by the caller -- an anchored
    domain rebuild, where that atom belongs to a peptide unit outside the span and must
    not move. Without this term the anchor residue's own tau is optimized against a
    placeholder atom and then measured against the real one.
    """
    fixed_unit = fixed_arm / np.linalg.norm(fixed_arm)
    moving_unit = moving_arm / np.linalg.norm(moving_arm, axis=-1, keepdims=True)
    tau = np.degrees(np.arccos(np.clip(moving_unit @ fixed_unit, -1.0, 1.0)))
    deviation = np.abs(tau - N_CA_C_ANGLE)
    cost: np.ndarray = (deviation / _TAU_SD) ** 2 + _INFEASIBLE_TAU_COST * (
        deviation > n_ca_c_tolerance
    )
    return cost


@dataclass(frozen=True, slots=True)
class _Segment:
    """One intact stretch of trace, with its peptide units precomputed on the grid.

    ``c_vec``, ``o_vec`` and ``n_vec`` are ``(n_units, grid, 3)`` *offsets*: ``c_vec[j,
    k]`` and ``o_vec[j, k]`` from ``ca[j]``, ``n_vec[j, k]`` from ``ca[j+1]``. Keeping
    them as offsets rather than positions is what lets a repair pass re-evaluate one
    unit's whole rotation range with two array operations.

    ``fixed_first_n`` and ``fixed_last_c`` / ``fixed_last_o``, when present, are atoms the
    caller owns: they are used by the objective and copied into the output rather than
    placed. See :func:`place_backbone_for_domain`.
    """

    start: int
    ca: np.ndarray
    sequence: str
    strain: np.ndarray
    c_vec: np.ndarray
    o_vec: np.ndarray
    n_vec: np.ndarray
    fixed_first_n: np.ndarray | None
    fixed_last_c: np.ndarray | None
    fixed_last_o: np.ndarray | None

    @property
    def n_residues(self) -> int:
        return int(self.ca.shape[0])

    @property
    def n_units(self) -> int:
        return int(self.ca.shape[0]) - 1


def _segment_geometry(
    ca: np.ndarray,
    sequence: str,
    start: int,
    phase: float,
    *,
    fixed_first_n: np.ndarray | None = None,
    fixed_last_c: np.ndarray | None = None,
    fixed_last_o: np.ndarray | None = None,
) -> _Segment:
    """Precompute one segment's peptide units at every grid rotation, or raise."""
    if ca.shape[0] < 2:
        raise GeometryError(
            "A chain break left a single isolated residue, whose C, O and N are not "
            "determined by anything. Trim the trace instead of asking for a backbone "
            "that cannot exist."
        )
    distance, x_axis, e1, e2 = _local_frames(ca)

    low, high = reconstructable_ca_bond_range()
    offending = np.flatnonzero((distance < low) | (distance > high))
    if offending.size:
        preview = ", ".join(
            f"{start + int(i)}-{start + int(i) + 1} at {distance[i]:.3f} A" for i in offending[:5]
        )
        more = f" (and {offending.size - 5} more)" if offending.size > 5 else ""
        raise GeometryError(
            f"{offending.size} CA-CA virtual bond(s) are outside the range a trans "
            f"peptide unit can span, {low:.3f}-{high:.3f} A: {preview}{more}. A "
            f"separation near 2.9 A is a cis peptide bond, which DODO does not model; a "
            f"much larger one is a chain break. Pass allow_chain_breaks=True to build "
            f"each intact stretch separately, or split the trace yourself."
        )

    strain = np.interp(distance, _strain_table()[1], _strain_table()[0])
    _, local = _peptide_unit_local(strain)

    step = np.radians(360.0 / _PLANE_GRID)
    theta = np.radians(phase) + np.arange(_PLANE_GRID) * step
    # y_hat[j, k] is unit j's in-plane perpendicular at grid rotation k.
    y_hat = (
        np.cos(theta)[None, :, None] * e1[:, None, :]
        + np.sin(theta)[None, :, None] * e2[:, None, :]
    )
    c_vec = local["C"][:, 0, None, None] * x_axis[:, None, :] + local["C"][:, 1, None, None] * y_hat
    o_vec = local["O"][:, 0, None, None] * x_axis[:, None, :] + local["O"][:, 1, None, None] * y_hat
    # N belongs to the unit's *second* residue, so express it relative to that CA.
    n_vec = (local["N"][:, 0] - distance)[:, None, None] * x_axis[:, None, :] + local["N"][
        :, 1, None, None
    ] * y_hat
    return _Segment(
        start=start,
        ca=ca,
        sequence=sequence,
        strain=strain,
        c_vec=c_vec,
        o_vec=o_vec,
        n_vec=n_vec,
        fixed_first_n=fixed_first_n,
        fixed_last_c=fixed_last_c,
        fixed_last_o=fixed_last_o,
    )


def _solve_planes(segment: _Segment, n_ca_c_tolerance: float, field: np.ndarray) -> np.ndarray:
    """Choose one peptide-plane rotation per unit; return the grid index of each.

    Exact global minimisation of the objective in :func:`_pair_cost`, plus the per-unit
    ``field``, over the discretized rotation grid, by Viterbi dynamic programming: every
    term couples only consecutive units, so the problem is a linear chain and the optimum
    is reachable without any iterative sampling and with no local-minimum risk. That
    matters more than it sounds, because the alternative -- greedy forward placement --
    has to commit to one of two roots at every residue and there are ``2 ** n`` root
    sequences that all satisfy ``N-CA-C`` exactly.

    ``field[j, k]`` is a unary cost on unit ``j`` taking rotation ``k``. It is how
    :func:`_repair_clashes` tells the chain about contacts that are *not* local, which
    the pairwise terms are structurally unable to see: before it existed, rebuilding
    dnmt3a 282-414 from its own deposited alpha carbons put two carbonyl oxygens nine
    residues apart at 1.282 A, shorter than an O=O double bond, and the objective's
    global optimum was exactly that structure.
    """
    grid = _PLANE_GRID
    n_units = segment.n_units
    accumulated = field[0].copy()
    if segment.fixed_first_n is not None:
        accumulated = accumulated + _fixed_arm_tau_cost(
            segment.fixed_first_n - segment.ca[0], segment.c_vec[0], n_ca_c_tolerance
        )

    backpointer = np.zeros((n_units, grid), dtype=np.int32)
    for residue in range(1, n_units):
        has_cb = ONE_TO_THREE.get(segment.sequence[residue], "UNK") not in _NO_CB_RESIDUES
        total = (
            accumulated[:, None]
            + _pair_cost(
                segment.ca,
                residue,
                segment.c_vec,
                segment.n_vec,
                segment.o_vec,
                has_cb,
                n_ca_c_tolerance,
            )
            + field[residue][None, :]
        )
        backpointer[residue] = np.argmin(total, axis=0)
        accumulated = total[backpointer[residue], np.arange(grid)]

    if segment.fixed_last_c is not None:
        accumulated = accumulated + _fixed_arm_tau_cost(
            segment.fixed_last_c - segment.ca[-1], segment.n_vec[-1], n_ca_c_tolerance
        )

    chosen = np.zeros(n_units, dtype=np.int64)
    chosen[n_units - 1] = int(np.argmin(accumulated))
    for residue in range(n_units - 1, 0, -1):
        chosen[residue - 1] = backpointer[residue][chosen[residue]]
    return chosen


# ---------------------------------------------------------------------------
# Non-bonded contacts
#
# Bond lengths, bond angles and the pairwise objective above are all LOCAL. A backbone
# can satisfy every one of them and still fold two carbonyl oxygens through each other,
# and one did: the section below is the part of this module that can tell.
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class StericClash:
    """One pair of non-bonded heavy atoms closer than their hard-sphere limit.

    Attributes
    ----------
    residue_i, residue_j
        0-based positional residue indices, ``residue_i < residue_j``.
    atom_i, atom_j
        PDB atom names.
    distance
        Measured separation in Angstroms.
    minimum
        The :data:`_MIN_CONTACT` limit for this element pair.
    """

    residue_i: int
    atom_i: str
    residue_j: int
    atom_j: str
    distance: float
    minimum: float

    @property
    def overlap(self) -> float:
        """How far inside the limit the pair sits, in Angstroms."""
        return self.minimum - self.distance

    def __str__(self) -> str:
        return (
            f"{self.atom_i}({self.residue_i})..{self.atom_j}({self.residue_j}) "
            f"at {self.distance:.3f} A, limit {self.minimum:.2f} A"
        )


@dataclass(frozen=True, slots=True)
class _AtomTable:
    """Flat atom list for contact queries: coordinates, residue index, name, element."""

    xyz: np.ndarray
    residue: np.ndarray
    name: np.ndarray
    code: np.ndarray


#: Element ordering behind :data:`_LIMIT_MATRIX`.
_ELEMENT_CODE: dict[str, int] = {"N": 0, "C": 1, "O": 2}

#: :data:`_MIN_CONTACT` as a lookup matrix over :data:`_ELEMENT_CODE`. DERIVED.
_LIMIT_MATRIX: np.ndarray = np.array(
    [
        [_MIN_CONTACT[frozenset((left, right))] for right in ("N", "C", "O")]
        for left in ("N", "C", "O")
    ]
)


def _atom_table(
    ca: np.ndarray,
    nitrogen: np.ndarray,
    carbon: np.ndarray,
    oxygen: np.ndarray,
    beta_carbon: np.ndarray,
) -> _AtomTable:
    """Flatten a model into one array per field, ``CB`` rows included where finite.

    ``CB`` is carried here even when the caller did not ask for side chains, because a
    backbone that can only avoid a clash by leaving out an atom every residue except
    glycine has is not a valid backbone. It is the atom that makes the disallowed regions
    of the Ramachandran plot disallowed.
    """
    n = int(ca.shape[0])
    with_cb = np.flatnonzero(np.isfinite(beta_carbon).all(axis=1))
    xyz = np.concatenate([nitrogen, ca, carbon, oxygen, beta_carbon[with_cb]])
    residue = np.concatenate([np.arange(n), np.arange(n), np.arange(n), np.arange(n), with_cb])
    name = np.array(
        ["N"] * n + ["CA"] * n + ["C"] * n + ["O"] * n + ["CB"] * int(with_cb.size), dtype="<U2"
    )
    code = np.concatenate(
        [
            np.full(n, _ELEMENT_CODE["N"]),
            np.full(n, _ELEMENT_CODE["C"]),
            np.full(n, _ELEMENT_CODE["C"]),
            np.full(n, _ELEMENT_CODE["O"]),
            np.full(with_cb.size, _ELEMENT_CODE["C"]),
        ]
    )
    return _AtomTable(xyz=xyz, residue=residue, name=name, code=code)


def _eligible_pairs(
    residue_low: np.ndarray,
    residue_high: np.ndarray,
    name_low: np.ndarray,
    name_high: np.ndarray,
    breaks: tuple[int, ...],
) -> np.ndarray:
    """Mask of atom pairs whose separation :data:`_MIN_CONTACT` actually bounds.

    Two residues apart or more, or one apart across a chain break -- where nothing joins
    the two residues at all, so the full non-bonded limit applies and ``C(i)..N(i+1)``
    is a contact rather than a peptide bond. See :data:`_NONBONDED_RESIDUE_SEPARATION`
    for the measurement that sets the cutoff.

    One exception, across a break only: the two alpha carbons. They are input, not output,
    and a break is *defined* by their separation being outside what a trans peptide unit
    spans -- which for a cis peptide bond means about 2.95 to 3.02 A, below the 3.20 A C-C
    limit. That distance is already reported, by :func:`reconstructable_ca_bond_range`
    refusing it and by :attr:`BackboneResult.chain_breaks` listing it, and no placement of
    any other atom can change it. Flagging it here would turn "DODO does not model cis
    peptide bonds" into "DODO cannot open a structure containing one", while saying nothing
    the caller does not already know.
    """
    separation = residue_high - residue_low
    eligible: np.ndarray = separation >= _NONBONDED_RESIDUE_SEPARATION
    if breaks:
        across = (separation == 1) & np.isin(residue_low, np.asarray(breaks))
        eligible = eligible | (across & ~((name_low == "CA") & (name_high == "CA")))
    return eligible


def _steric_survey(
    table: _AtomTable, breaks: tuple[int, ...], radius: float
) -> tuple[tuple[StericClash, ...], float]:
    """Every non-bonded pair inside its limit, plus the closest eligible separation.

    ``radius`` bounds the neighbour search. Any radius at least :data:`_MAX_CONTACT`
    finds every violation; a larger one also makes the reported closest contact
    meaningful on a small or extended trace where nothing is within a limit at all.
    """
    pairs = cKDTree(table.xyz).query_pairs(radius, output_type="ndarray")
    if pairs.shape[0] == 0:
        return (), float("inf")
    residue_a, residue_b = table.residue[pairs[:, 0]], table.residue[pairs[:, 1]]
    first_is_low = residue_a <= residue_b
    low = np.where(first_is_low, pairs[:, 0], pairs[:, 1])
    high = np.where(first_is_low, pairs[:, 1], pairs[:, 0])
    eligible = _eligible_pairs(
        table.residue[low], table.residue[high], table.name[low], table.name[high], breaks
    )
    if not bool(eligible.any()):
        return (), float("inf")
    low, high = low[eligible], high[eligible]
    limits = _LIMIT_MATRIX[table.code[low], table.code[high]]
    distance = np.linalg.norm(table.xyz[low] - table.xyz[high], axis=1)
    closest = float(distance.min())
    inside = np.flatnonzero(distance < limits)
    clashes = [
        StericClash(
            residue_i=int(table.residue[low[index]]),
            atom_i=str(table.name[low[index]]),
            residue_j=int(table.residue[high[index]]),
            atom_j=str(table.name[high[index]]),
            distance=float(distance[index]),
            minimum=float(limits[index]),
        )
        for index in inside
    ]
    clashes.sort(key=lambda clash: (-clash.overlap, clash.residue_i, clash.residue_j))
    return tuple(clashes), closest


# ---------------------------------------------------------------------------
# Assembling a model, and repairing what the pairwise objective cannot see
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class _Model:
    """One complete reconstruction of a whole trace, and what is wrong with it."""

    ca: np.ndarray
    nitrogen: np.ndarray
    carbon: np.ndarray
    oxygen: np.ndarray
    beta_carbon: np.ndarray
    clashes: tuple[StericClash, ...]

    @property
    def rank(self) -> tuple[int, float]:
        """Sort key: fewest clashes first, then least total overlap."""
        return len(self.clashes), float(sum(clash.overlap for clash in self.clashes))


def _terminal_dihedral_candidates(default: float) -> np.ndarray:
    """Dihedral values to try for a terminal atom, nearest the default first.

    Ordering matters: the search takes the first minimum, so a terminal atom that has no
    clash to avoid keeps exactly the documented :data:`_TERMINAL_PHI` /
    :data:`_TERMINAL_PSI` value rather than drifting to an equally good neighbour.
    """
    step = 360.0 / _TERMINAL_DIHEDRAL_GRID
    offsets = [0.0]
    for index in range(1, _TERMINAL_DIHEDRAL_GRID // 2 + 1):
        offsets.append(index * step)
        offsets.append(-index * step)
    return default + np.array(offsets[:_TERMINAL_DIHEDRAL_GRID])


def _choose_terminal(
    positions: np.ndarray,
    name: str,
    residue: int,
    table: _AtomTable,
    breaks: tuple[int, ...],
) -> int:
    """Index of the candidate position with the least steric overlap.

    A terminal atom is placed from an arbitrary dihedral because the CA trace does not
    determine it. That is fine in isolation and not fine next to another segment: with
    ``allow_chain_breaks``, two segments' terminal atoms used to be placed from fixed
    dihedrals with no knowledge of each other, and across dnmt3a's cis-prolines that put
    ``C(i)`` and ``N(i+1)`` at 1.124 A -- shorter than a C-N triple bond -- which the
    writer then emitted a CONECT record for.
    """
    eligible = _eligible_against(residue, name, table, breaks)
    partners = table.xyz[eligible]
    if partners.shape[0] == 0:
        return 0
    limits = _LIMIT_MATRIX[_ELEMENT_CODE[_ELEMENT_OF[name]], table.code[eligible]]
    distance = np.linalg.norm(positions[:, None, :] - partners[None, :, :], axis=-1)
    excess = np.clip(limits[None, :] - distance, 0.0, None)
    cost = _INFEASIBLE_CLASH_COST * (excess > 0.0).sum(axis=1) + _STERIC_WEIGHT * (excess**2).sum(
        axis=1
    )
    return int(np.argmin(cost))


def _assemble(
    ca: np.ndarray,
    residue_names: np.ndarray,
    segments: list[_Segment],
    chosen: list[np.ndarray],
    breaks: tuple[int, ...],
) -> _Model:
    """Place every atom of a whole trace, refine the terminals, and measure the result.

    Terminal refinement runs twice. Once is not enough when two free terminal atoms face
    each other across a chain break: the second sweep lets each see where the other
    actually ended up.
    """
    n_residues = int(ca.shape[0])
    nitrogen = np.empty((n_residues, 3))
    carbon = np.empty((n_residues, 3))
    oxygen = np.empty((n_residues, 3))

    for index, segment in enumerate(segments):
        start = segment.start
        stop = start + segment.n_residues
        units = np.arange(segment.n_units)
        picked = chosen[index]
        carbon[start : stop - 1] = segment.ca[:-1] + segment.c_vec[units, picked]
        oxygen[start : stop - 1] = segment.ca[:-1] + segment.o_vec[units, picked]
        nitrogen[start + 1 : stop] = segment.ca[1:] + segment.n_vec[units, picked]
        if segment.fixed_first_n is not None:
            nitrogen[start] = segment.fixed_first_n
        else:
            nitrogen[start] = _place_atom(
                nitrogen[start + 1],
                carbon[start],
                ca[start],
                N_CA_BOND_LENGTH,
                N_CA_C_ANGLE,
                _TERMINAL_PSI,
            )
        if segment.fixed_last_c is not None and segment.fixed_last_o is not None:
            carbon[stop - 1] = segment.fixed_last_c
            oxygen[stop - 1] = segment.fixed_last_o
        else:
            carbon[stop - 1] = _place_atom(
                carbon[stop - 2],
                nitrogen[stop - 1],
                ca[stop - 1],
                CA_C_BOND_LENGTH,
                N_CA_C_ANGLE,
                _TERMINAL_PHI,
            )
            oxygen[stop - 1] = _place_atom(
                nitrogen[stop - 1],
                ca[stop - 1],
                carbon[stop - 1],
                C_O_BOND_LENGTH,
                CA_C_O_ANGLE,
                _TERMINAL_PSI - 180.0,
            )

    beta_carbon = _beta_carbons(ca, nitrogen, carbon, residue_names)
    for _sweep in range(2):
        for segment in segments:
            table = _atom_table(ca, nitrogen, carbon, oxygen, beta_carbon)
            start = segment.start
            stop = start + segment.n_residues
            if segment.fixed_first_n is None:
                candidates = _place_atom_over_dihedrals(
                    nitrogen[start + 1],
                    carbon[start],
                    ca[start],
                    N_CA_BOND_LENGTH,
                    N_CA_C_ANGLE,
                    _terminal_dihedral_candidates(_TERMINAL_PSI),
                )
                nitrogen[start] = candidates[
                    _choose_terminal(candidates, "N", start, table, breaks)
                ]
            if segment.fixed_last_c is None:
                last = stop - 1
                candidates = _place_atom_over_dihedrals(
                    carbon[last - 1],
                    nitrogen[last],
                    ca[last],
                    CA_C_BOND_LENGTH,
                    N_CA_C_ANGLE,
                    _terminal_dihedral_candidates(_TERMINAL_PHI),
                )
                carbon[last] = candidates[_choose_terminal(candidates, "C", last, table, breaks)]
                oxygen_candidates = _place_atom_over_dihedrals(
                    nitrogen[last],
                    ca[last],
                    carbon[last],
                    C_O_BOND_LENGTH,
                    CA_C_O_ANGLE,
                    _terminal_dihedral_candidates(_TERMINAL_PSI) - 180.0,
                )
                oxygen[last] = oxygen_candidates[
                    _choose_terminal(oxygen_candidates, "O", last, table, breaks)
                ]
            beta_carbon = _beta_carbons(ca, nitrogen, carbon, residue_names)

    table = _atom_table(ca, nitrogen, carbon, oxygen, beta_carbon)
    clashes, _ = _steric_survey(table, breaks, _MAX_CONTACT)
    return _Model(
        ca=ca,
        nitrogen=nitrogen,
        carbon=carbon,
        oxygen=oxygen,
        beta_carbon=beta_carbon,
        clashes=clashes,
    )


def _beta_carbons(
    ca: np.ndarray, nitrogen: np.ndarray, carbon: np.ndarray, residue_names: np.ndarray
) -> np.ndarray:
    """``CB`` for every residue that has one, NaN for glycine and ``UNK``."""
    beta_carbon = place_sidechain_cb(nitrogen, ca, carbon)
    absent = np.array([str(name) in _NO_CB_RESIDUES for name in residue_names])
    beta_carbon[absent] = np.nan
    return beta_carbon


def _movable_positions(
    segment: _Segment, residue: int, name: str, model: _Model
) -> list[tuple[int, np.ndarray]]:
    """Rotation-dependent positions of one atom, as ``(unit, (grid, 3))`` entries.

    An atom placed by a peptide unit sweeps a circle as that unit's plane rotates, so a
    contact it is involved in can be scored over the whole grid at once. ``CB`` appears
    twice because it depends on both of its residue's neighbouring units, through ``N``
    and through ``C``; the other one is held at its current position, which is what makes
    this a mean field rather than an exact joint cost.

    Atoms that no unit owns -- alpha carbons, terminal atoms, caller-supplied anchor
    atoms -- return nothing. They are immovable here by construction, which is why
    :func:`_choose_terminal` exists.
    """
    local = residue - segment.start
    centre = segment.ca[local]
    out: list[tuple[int, np.ndarray]] = []
    if name == "C" and local < segment.n_units:
        out.append((local, centre + segment.c_vec[local]))
    elif name == "O" and local < segment.n_units:
        out.append((local, centre + segment.o_vec[local]))
    elif name == "N" and local >= 1:
        out.append((local - 1, centre + segment.n_vec[local - 1]))
    elif name == "CB":
        if local < segment.n_units:
            out.append(
                (
                    local,
                    place_sidechain_cb(
                        model.nitrogen[residue], centre, centre + segment.c_vec[local]
                    ),
                )
            )
        if local >= 1:
            out.append(
                (
                    local - 1,
                    place_sidechain_cb(
                        centre + segment.n_vec[local - 1], centre, model.carbon[residue]
                    ),
                )
            )
    return out


def _owning_segment(segments: list[_Segment], residue: int) -> int | None:
    """Index of the segment containing ``residue``, or None."""
    for index, segment in enumerate(segments):
        if segment.start <= residue < segment.start + segment.n_residues:
            return index
    return None


def _candidate_atoms(segment: _Segment, model: _Model) -> list[tuple[str, np.ndarray, np.ndarray]]:
    """Every rotation-dependent atom of a segment, as ``(name, residues, positions)``.

    ``positions`` is ``(n_units, grid, 3)`` and ``residues`` is ``(n_units,)``: entry
    ``[j, k]`` is where that atom sits if unit ``j`` takes rotation ``k``. ``CB`` appears
    twice, once for each of the two units it depends on, with the other one held at its
    current position -- that is the mean-field approximation, and it is the only one in
    here.
    """
    units = np.arange(segment.n_units)
    first = segment.start + units
    second = first + 1
    out: list[tuple[str, np.ndarray, np.ndarray]] = [
        ("C", first, segment.ca[:-1, None, :] + segment.c_vec),
        ("O", first, segment.ca[:-1, None, :] + segment.o_vec),
        ("N", second, segment.ca[1:, None, :] + segment.n_vec),
    ]
    has_cb = np.isfinite(model.beta_carbon).all(axis=1)
    if bool(has_cb[first].any()):
        cb_first = place_sidechain_cb(
            model.nitrogen[first][:, None, :],
            segment.ca[:-1, None, :],
            out[0][2],
        )
        # A residue with no CB must not contribute a phantom contact. Pushing the point
        # far away is how it is dropped, since the query is by radius.
        cb_first[~has_cb[first]] = _FAR_AWAY
        out.append(("CB", first, cb_first))
    if bool(has_cb[second].any()):
        cb_second = place_sidechain_cb(
            out[2][2],
            segment.ca[1:, None, :],
            model.carbon[second][:, None, :],
        )
        cb_second[~has_cb[second]] = _FAR_AWAY
        out.append(("CB", second, cb_second))
    return out


def _clash_field(
    segment: _Segment, model: _Model, table: _AtomTable, breaks: tuple[int, ...]
) -> np.ndarray:
    """Long-range contact cost of every rotation of every unit, as ``(n_units, grid)``.

    This is the term the pairwise chain structurally cannot express, and its absence is
    what let the objective's *global* optimum be a structure with two carbonyl oxygens
    nine residues apart at 1.282 A -- shorter than an O=O double bond -- while every local
    term was satisfied. Each unit's atoms are swept over the whole rotation grid and
    scored against every already-placed atom the contact limits apply to, so the Viterbi
    pass minimises a cost that can actually see the rest of the protein.

    The approximation is that the partners are held at their current positions, which is
    stale for anything that also moves. That is why the result is iterated (see
    :func:`_repair_clashes`) and why exact single-unit descent runs afterwards. Pairs
    already scored exactly by :func:`_pair_cost` are skipped unless they span a chain
    break, where no local term exists: an estimate must not overrule an exact term.
    """
    field = np.zeros((segment.n_units, _PLANE_GRID))
    if table.xyz.shape[0] == 0:
        return field
    tree = cKDTree(table.xyz)
    break_array = np.asarray(breaks) if breaks else None
    for name, residues, positions in _candidate_atoms(segment, model):
        flat = positions.reshape(-1, 3)
        pairs = cKDTree(flat).sparse_distance_matrix(tree, _MAX_CONTACT, output_type="ndarray")
        if pairs.shape[0] == 0:
            continue
        point = pairs["i"]
        partner = pairs["j"]
        distance = pairs["v"]
        unit = point // _PLANE_GRID
        rotation = point % _PLANE_GRID
        atom_residue = residues[unit]
        partner_residue = table.residue[partner]
        separation = np.abs(partner_residue - atom_residue)
        eligible = separation >= _NONBONDED_RESIDUE_SEPARATION
        if break_array is not None:
            spans_break = (separation == 1) & np.isin(
                np.minimum(partner_residue, atom_residue), break_array
            )
            eligible |= spans_break
        else:
            spans_break = np.zeros(eligible.shape, dtype=bool)
        # Drop the pairs _pair_cost owns. Doing it by key rather than by residue window
        # keeps the two terms exactly complementary: nothing is scored twice and nothing
        # is missed.
        earlier_is_partner = partner_residue < atom_residue
        for other in ("N", "CA", "C", "O", "CB"):
            same = table.name[partner] == other
            if not bool(same.any()):
                continue
            for offset in (0, 1):
                keyed = same & (separation == offset)
                if not bool(keyed.any()):
                    continue
                low_first = _pair_key(offset, other, name) in _LOCALLY_SCORED
                high_first = _pair_key(offset, name, other) in _LOCALLY_SCORED
                covered = np.where(earlier_is_partner, low_first, high_first)
                eligible &= ~(keyed & covered & ~spans_break)
        if not bool(eligible.any()):
            continue
        unit, rotation, distance = unit[eligible], rotation[eligible], distance[eligible]
        limits = _LIMIT_MATRIX[_ELEMENT_CODE[_ELEMENT_OF[name]], table.code[partner[eligible]]]
        excess = np.clip(limits - distance, 0.0, None)
        contributions = _CLASH_FIELD_COST * (excess > 0.0) + _STERIC_WEIGHT * excess**2
        np.add.at(field, (unit, rotation), contributions)
    return field


def _unit_atoms(segment: _Segment, unit: int, model: _Model) -> list[tuple[int, str, np.ndarray]]:
    """Atoms whose position depends on one unit's rotation, as ``(residue, name, (grid, 3))``.

    ``CB`` of either flanking residue is included when that residue has one: it depends on
    this unit through ``N`` or through ``C``, with the other held at its current position.
    """
    residue = segment.start + unit
    carbon = segment.ca[unit] + segment.c_vec[unit]
    nitrogen = segment.ca[unit + 1] + segment.n_vec[unit]
    moving: list[tuple[int, str, np.ndarray]] = [
        (residue, "C", carbon),
        (residue, "O", segment.ca[unit] + segment.o_vec[unit]),
        (residue + 1, "N", nitrogen),
    ]
    if bool(np.isfinite(model.beta_carbon[residue]).all()):
        moving.append(
            (
                residue,
                "CB",
                place_sidechain_cb(model.nitrogen[residue], segment.ca[unit], carbon),
            )
        )
    if bool(np.isfinite(model.beta_carbon[residue + 1]).all()):
        moving.append(
            (
                residue + 1,
                "CB",
                place_sidechain_cb(nitrogen, segment.ca[unit + 1], model.carbon[residue + 1]),
            )
        )
    return moving


def _eligible_against(
    atom_residue: int, name: str, table: _AtomTable, breaks: tuple[int, ...]
) -> np.ndarray:
    """Mask of table rows whose contact with a named atom at ``atom_residue`` is bounded.

    The array form of the rule :func:`_eligible_pairs` applies to whole pair lists, for the
    case of one moving atom against everything else.
    """
    separation = np.abs(table.residue - atom_residue)
    eligible: np.ndarray = separation >= _NONBONDED_RESIDUE_SEPARATION
    if breaks:
        lower = np.minimum(table.residue, atom_residue)
        across = (separation == 1) & np.isin(lower, np.asarray(breaks))
        if name == "CA":
            across = across & (table.name != "CA")
        eligible = eligible | across
    return eligible


def _contact_cost_over_grid(
    moving: list[tuple[int, str, np.ndarray]],
    table: _AtomTable,
    breaks: tuple[int, ...],
) -> np.ndarray:
    """Exact contact cost, per rotation, of a set of moving atoms against fixed ones."""
    cost = np.zeros(_PLANE_GRID)
    for atom_residue, name, positions in moving:
        eligible = _eligible_against(atom_residue, name, table, breaks)
        partners = table.xyz[eligible]
        if partners.shape[0] == 0:
            continue
        limits = _LIMIT_MATRIX[_ELEMENT_CODE[_ELEMENT_OF[name]], table.code[eligible]]
        distance = np.linalg.norm(positions[:, None, :] - partners[None, :, :], axis=-1)
        excess = np.clip(limits[None, :] - distance, 0.0, None)
        cost += _INFEASIBLE_CLASH_COST * (excess > 0.0).sum(axis=1) + _STERIC_WEIGHT * (
            excess**2
        ).sum(axis=1)
    return cost


def _unit_tau_cost(
    segment: _Segment, unit: int, model: _Model, n_ca_c_tolerance: float
) -> np.ndarray:
    """Both ``N-CA-C`` costs one unit takes part in, per rotation, as hard constraints."""
    residue = segment.start + unit
    cost: np.ndarray = _fixed_arm_tau_cost(
        model.nitrogen[residue] - segment.ca[unit], segment.c_vec[unit], n_ca_c_tolerance
    ) + _fixed_arm_tau_cost(
        model.carbon[residue + 1] - segment.ca[unit + 1], segment.n_vec[unit], n_ca_c_tolerance
    )
    return cost


def _rotation_cost(
    segment: _Segment,
    unit: int,
    model: _Model,
    table: _AtomTable,
    breaks: tuple[int, ...],
    n_ca_c_tolerance: float,
) -> np.ndarray:
    """Exact cost of every rotation of one peptide unit, everything else held fixed.

    "Exact" is the point. :func:`_clash_field` estimates a long-range contact against the
    *other* atom's current position, which is stale as soon as that atom moves too, and two
    partners that both flee the same overlap can chase each other -- measured, a CB-CB pair
    oscillated through 2.17, 2.53, 2.99, 2.62 A over four passes without landing. This
    function has no such blind spot: it scores the real contact set of the real model as a
    function of one rotation. What it cannot do is move two units at once, which is what
    the Viterbi pass against :func:`_clash_field` is for.

    Includes both ``N-CA-C`` angles that this unit takes part in, as hard constraints, so
    the descent cannot trade a clash for an angle the caller will refuse.

    Every atom that moves with this unit lives in the two residues the unit joins, both of
    which are within one residue of each moving atom, so the eligibility rule has already
    excluded them: no atom is ever scored against a stale copy of itself.
    """
    moving = _unit_atoms(segment, unit, model)
    cost: np.ndarray = _contact_cost_over_grid(moving, table, breaks) + _unit_tau_cost(
        segment, unit, model, n_ca_c_tolerance
    )
    return cost


def _polish_clashes(
    ca: np.ndarray,
    residue_names: np.ndarray,
    segments: list[_Segment],
    chosen: list[np.ndarray],
    model: _Model,
    breaks: tuple[int, ...],
    n_ca_c_tolerance: float,
) -> _Model:
    """Exact descent on the units involved in a clash: one unit at a time, then two.

    Only clash-involved units are visited, because they are the only ones whose rotation
    can remove a measured contact. Every accepted move is verified by re-measuring the
    whole model, so the clash count and total overlap are monotone: this can only improve a
    model, never damage one. ``chosen`` is mutated so the caller's rotations stay in step
    with the returned model.

    Coordinated moves of two or more units are *not* attempted here. They were, exhaustively
    and exactly, over both units owning each clashing atom pair: measured over 28 traces
    (four dnmt3a domains and 24 self-avoiding walks of 20 to 150 residues) it changed the
    outcome of exactly none of them while costing 60 per cent more time, because on the
    traces that still fail every joint alternative violates either the N-CA-C tolerance or
    another contact. Multi-unit escape is the Viterbi pass's job, and it has the whole chain
    to work with rather than two units.
    """
    for _sweep in range(_MAX_POLISH_SWEEPS):
        if not model.clashes:
            break
        changed = False
        for index, unit in _clash_units(segments, model):
            table = _atom_table(ca, model.nitrogen, model.carbon, model.oxygen, model.beta_carbon)
            cost = _rotation_cost(segments[index], unit, model, table, breaks, n_ca_c_tolerance)
            current = int(chosen[index][unit])
            best = int(np.argmin(cost))
            # Strictly better only, and verified on the reassembled model rather than
            # trusted from the per-unit cost. The per-unit cost is exact for the pairs it
            # scores, but reassembly also re-chooses the terminal dihedrals, so the only
            # way to be sure a move is an improvement is to measure the whole model. An
            # unverified descent step made a 3.184 A contact into a 2.900 A one.
            if cost[best] >= cost[current] - 1e-9:
                continue
            chosen[index][unit] = best
            candidate_model = _assemble(ca, residue_names, segments, chosen, breaks)
            if candidate_model.rank < model.rank:
                model = candidate_model
                changed = True
            else:
                chosen[index][unit] = current
        if not changed:
            break
    return model


def _clash_units(segments: list[_Segment], model: _Model) -> list[tuple[int, int]]:
    """Every ``(segment, unit)`` whose rotation moves an atom in a measured clash."""
    units: list[tuple[int, int]] = []
    for clash in model.clashes:
        for residue, name in ((clash.residue_i, clash.atom_i), (clash.residue_j, clash.atom_j)):
            index = _owning_segment(segments, residue)
            if index is None:
                continue
            for unit, _ in _movable_positions(segments[index], residue, name, model):
                if (index, unit) not in units:
                    units.append((index, unit))
    return units


def _repair_clashes(
    ca: np.ndarray,
    residue_names: np.ndarray,
    segments: list[_Segment],
    n_ca_c_tolerance: float,
    breaks: tuple[int, ...],
) -> _Model:
    """Solve, measure non-bonded contacts, penalize, re-solve; return the best model.

    Two mechanisms, because the two failure modes are different. A whole-chain Viterbi
    re-solve against the long-range field (:func:`_clash_field`) can flip several peptide
    units together, which is what a deep overlap needs: the 1.282 A O..O contact in
    dnmt3a 282-414 cannot be fixed by rotating either unit alone, because the correct
    rotation only becomes available once the neighbouring flipped units move too. Exact
    single-unit descent (:func:`_polish_clashes`) then removes what is left, without the
    stale-neighbour blind spot the field has.

    The field is *accumulated* rather than recomputed from scratch, which is what stops two
    partners from fleeing the same overlap into each other's old positions over and over:
    a rotation that was bad in an earlier pass stays charged. It is capped
    (:data:`_MAX_CLASH_FIELD`) so that a long enough search can never outweigh an
    exactly-scored constraint.

    The loop exits as soon as a model has no clashes, or as soon as a pass has nothing
    left to penalize -- which happens when every remaining clash is between atoms no
    peptide-plane rotation can move, i.e. alpha carbons of the input trace. The best model
    seen is returned either way, and :func:`place_backbone` refuses it if it is not clean.
    Reporting a rebuilt backbone with atoms at 1.28 A as a success is the defect this whole
    section exists to prevent.
    """
    fields = [np.zeros((segment.n_units, _PLANE_GRID)) for segment in segments]
    chosen = [
        _solve_planes(segment, n_ca_c_tolerance, fields[index])
        for index, segment in enumerate(segments)
    ]
    best = _polish_clashes(
        ca,
        residue_names,
        segments,
        chosen,
        _assemble(ca, residue_names, segments, chosen, breaks),
        breaks,
        n_ca_c_tolerance,
    )
    for _pass in range(_MAX_CLASH_REPAIR_PASSES):
        if not best.clashes:
            break
        # Measure the field on the best model so far rather than on whatever the last pass
        # produced. Otherwise the search drifts: contacts get rediscovered in new places
        # and the field fills up with penalties on rotations no good solution uses.
        table = _atom_table(ca, best.nitrogen, best.carbon, best.oxygen, best.beta_carbon)
        added = False
        for index, segment in enumerate(segments):
            penalty = _clash_field(segment, best, table, breaks)
            if not bool((penalty > 0.0).any()):
                continue
            fields[index] = np.minimum(fields[index] + penalty, _MAX_CLASH_FIELD)
            added = True
        if not added:
            break
        chosen = [
            _solve_planes(segment, n_ca_c_tolerance, fields[index])
            for index, segment in enumerate(segments)
        ]
        candidate = _polish_clashes(
            ca,
            residue_names,
            segments,
            chosen,
            _assemble(ca, residue_names, segments, chosen, breaks),
            breaks,
            n_ca_c_tolerance,
        )
        if candidate.rank < best.rank:
            best = candidate
    return best


# ---------------------------------------------------------------------------
# Result
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class BackboneResult:
    """A reconstructed backbone as flat atom records, plus what it cost geometrically.

    The coordinate arrays are laid out exactly as
    :meth:`~dodo.structure.Structure.from_atom_records` wants them: atoms grouped by
    residue, in ``N, CA, C, O[, CB]`` order within each residue.

    The geometry fields are not decoration. Bond lengths and omega are exact here by
    construction, but ``N-CA-C`` is not always achievable (see
    :func:`max_reconstructable_ca_angle`), and neither is a completely unstrained peptide
    angle. Rather than average that away, the per-residue and per-unit numbers are
    carried out so a caller can report them -- or refuse the result.

    Attributes
    ----------
    xyz
        ``(n_atoms, 3)`` coordinates.
    atom_names
        ``(n_atoms,)`` PDB atom names.
    elements
        ``(n_atoms,)`` element symbols.
    residue_index
        ``(n_atoms,)`` 0-based positional residue index per atom.
    residue_names
        ``(n_residues,)`` three-letter residue names, in residue order.
    n_ca_c_angles
        ``(n_residues,)`` measured ``N-CA-C`` angle in degrees.
    peptide_angle_strain
        ``(n_residues - 1,)`` signed strain in degrees applied to each peptide unit's
        bond angles to close onto the given CA-CA distance. See
        :data:`_MAX_PEPTIDE_ANGLE_STRAIN`. ``NaN`` at an index listed in
        :attr:`chain_breaks`, where there is no peptide unit to strain -- NaN rather than
        0.0, because 0.0 would read as "an unstrained peptide bond" when in fact there is
        no peptide bond at all.
    chain_breaks
        Residue indices ``i`` such that no peptide bond was built between ``i`` and
        ``i + 1``. Empty unless ``allow_chain_breaks`` was set. The ``C(i)`` to
        ``N(i+1)`` distance across one of these is **not** a peptide bond and must not be
        treated as one -- but it is still held to the non-bonded contact limit, because
        two atoms that are not bonded cannot be 1.1 A apart either.
    min_nonbonded_contact
        Closest separation, in Angstroms, between any two heavy atoms that
        :data:`_MIN_CONTACT` bounds -- two residues apart or more, or either side of a
        chain break. ``inf`` if no such pair is within measuring range of a limit. This
        is the number that was missing: a rebuild with two carbonyl oxygens at 1.282 A
        used to be reported with no field, no warning and no exception.
    steric_clashes
        Every pair below its limit, worst overlap first. Always empty on a result
        :func:`place_backbone` returned, because it refuses rather than returning a
        clashing backbone; the field exists so the measurement is on the object rather
        than only in an error message, and because ``CB`` is included in it even when
        ``sidechains`` is False.
    sidechains
        Whether ``CB`` atoms are present.
    """

    xyz: np.ndarray
    atom_names: np.ndarray
    elements: np.ndarray
    residue_index: np.ndarray
    residue_names: np.ndarray
    n_ca_c_angles: np.ndarray
    peptide_angle_strain: np.ndarray
    chain_breaks: tuple[int, ...]
    min_nonbonded_contact: float
    steric_clashes: tuple[StericClash, ...]
    sidechains: bool

    def __repr__(self) -> str:
        return (
            f"BackboneResult({self.n_atoms} atoms, {self.n_residues} residues, "
            f"max |N-CA-C - {N_CA_C_ANGLE}| = {self.max_n_ca_c_deviation:.2f} deg, "
            f"max |strain| = {self.max_peptide_angle_strain:.2f} deg, "
            f"closest non-bonded contact {self.min_nonbonded_contact:.2f} A"
            + (f", {len(self.chain_breaks)} chain breaks" if self.chain_breaks else "")
            + (f", {len(self.steric_clashes)} CLASHES" if self.steric_clashes else "")
            + (", with CB" if self.sidechains else "")
            + ")"
        )

    @property
    def n_atoms(self) -> int:
        """Number of atoms placed."""
        return int(self.xyz.shape[0])

    @property
    def n_residues(self) -> int:
        """Number of residues."""
        return int(self.residue_names.shape[0])

    @property
    def max_n_ca_c_deviation(self) -> float:
        """Largest absolute deviation of ``N-CA-C`` from its constant, in degrees."""
        return float(np.max(np.abs(self.n_ca_c_angles - N_CA_C_ANGLE)))

    @property
    def max_peptide_angle_strain(self) -> float:
        """Largest absolute peptide bond-angle strain used, in degrees.

        NaN entries -- chain-break positions, which have no peptide unit -- are skipped
        rather than propagated, so a single break does not turn the summary of 400 good
        peptide units into ``nan``.
        """
        finite = self.peptide_angle_strain[np.isfinite(self.peptide_angle_strain)]
        if finite.size == 0:
            return 0.0
        return float(np.max(np.abs(finite)))

    def atom_records(self, *, chain_id: str = "A", first_residue_number: int = 1) -> dict[str, Any]:
        """Keyword arguments for :meth:`~dodo.structure.Structure.from_atom_records`.

        Provided because the per-atom expansion of the residue-level fields is fiddly
        and getting it wrong produces a structure that validates and is wrong.
        """
        return {
            "xyz": self.xyz,
            "atom_name": list(self.atom_names),
            "element": list(self.elements),
            "residue_name": [str(self.residue_names[i]) for i in self.residue_index],
            "residue_number": [first_residue_number + int(i) for i in self.residue_index],
            "chain_id": [chain_id] * self.n_atoms,
        }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def _validate_inputs(ca_coords: np.ndarray, sequence: str) -> tuple[np.ndarray, str]:
    """Normalize and check the CA array and sequence, or raise."""
    ca = np.ascontiguousarray(np.asarray(ca_coords, dtype=np.float64))
    if ca.ndim != 2 or ca.shape[1] != 3:
        raise GeometryError(f"ca_coords must have shape (n_residues, 3), got {ca.shape}.")
    if not bool(np.isfinite(ca).all()):
        bad = int(np.flatnonzero(~np.isfinite(ca).all(axis=1))[0])
        raise GeometryError(
            f"ca_coords contains non-finite values, first at residue {bad}. A NaN here "
            f"is a failed upstream sampler, and building a backbone on it would launder "
            f"that failure into a structure file."
        )
    cleaned = sequence.strip().upper()
    if len(cleaned) != ca.shape[0]:
        raise GeometryError(
            f"sequence has {len(cleaned)} residues but ca_coords has {ca.shape[0]}. "
            f"These index the same residues, so a mismatch means one of them is for a "
            f"different region."
        )
    if ca.shape[0] < 2:
        raise GeometryError(
            f"Cannot build a backbone from {ca.shape[0]} residue(s): a peptide unit needs "
            f"two consecutive alpha carbons. Nothing determines C, O or N for an isolated "
            f"residue, so there is no honest answer to return."
        )
    return ca, cleaned


def _segment_bounds(ca: np.ndarray, allow_chain_breaks: bool) -> list[tuple[int, int]]:
    """Split the trace at CA-CA separations no trans peptide unit can span.

    Returns half-open ``[start, stop)`` residue ranges. With ``allow_chain_breaks``
    False this is always a single segment covering everything, and
    :func:`_solve_planes` raises on the offending bond -- which is the right default,
    because a break is usually a caller mistake rather than an intention.
    """
    if not allow_chain_breaks:
        return [(0, ca.shape[0])]
    distance = np.linalg.norm(np.diff(ca, axis=0), axis=1)
    low, high = reconstructable_ca_bond_range()
    breaks = np.flatnonzero((distance < low) | (distance > high))
    bounds: list[tuple[int, int]] = []
    start = 0
    for index in breaks:
        bounds.append((start, int(index) + 1))
        start = int(index) + 1
    bounds.append((start, ca.shape[0]))
    return bounds


def _validate_tolerance(n_ca_c_tolerance: float) -> float:
    """Check the ``N-CA-C`` tolerance, or raise.

    Non-finite is refused, and that is a behaviour change with a reason. The old error
    message recommended ``n_ca_c_tolerance=float("inf")`` as a way to inspect a failing
    trace. Measured, that turned a loud failure into a silent one: ten walk-engine
    conformers rebuilt that way came back looking like any other success while carrying
    ``N-CA-C`` angles from 101.6 to 132.9 degrees, 11.6 per cent of residues more than 8
    degrees out, and 29 heavy-atom pairs under 2.4 A. To inspect a trace before building,
    call :func:`unreconstructable_ca_angles`; it answers the same question without
    producing a structure file that nobody can tell is wrong.
    """
    tolerance = float(n_ca_c_tolerance)
    if not np.isfinite(tolerance):
        raise GeometryError(
            f"n_ca_c_tolerance must be finite, got {n_ca_c_tolerance}. An infinite "
            f"tolerance does not inspect a failing trace, it accepts one: the returned "
            f"result is indistinguishable from a success while carrying N-CA-C angles no "
            f"residue exhibits. Call unreconstructable_ca_angles(ca_coords) to find out "
            f"which residues cannot work, and why."
        )
    if tolerance <= 0.0:
        raise GeometryError(
            f"n_ca_c_tolerance must be positive, got {tolerance}. Zero would demand an "
            f"exact {N_CA_C_ANGLE} deg N-CA-C at every residue, which no real structure "
            f"achieves either."
        )
    return tolerance


def place_backbone(
    ca_coords: np.ndarray,
    sequence: str,
    rng: np.random.Generator,
    sidechains: bool = False,
    *,
    allow_chain_breaks: bool = False,
    n_ca_c_tolerance: float = _N_CA_C_TOLERANCE,
) -> BackboneResult:
    """Reconstruct ``N``, ``C``, ``O`` -- and optionally ``CB`` -- from a CA trace.

    The alpha carbons are **input, not output**: they are copied into the result
    unchanged, bit for bit. That is not a nicety, it is the check that the reconstruction
    is self-consistent, since a trans peptide unit built from the declared constants spans
    exactly the declared CA-CA virtual bond (3.8040 A against 3.80). See the module
    docstring for the derivation.

    Three things can make this refuse, and all three are physical rather than stylistic:
    a CA-CA virtual bond no trans peptide unit can span, an ``N-CA-C`` angle the trace
    makes impossible, or a non-bonded contact inside its hard-sphere limit. The last one
    is checked over the whole trace, not just between neighbours; see
    :data:`_NONBONDED_RESIDUE_SEPARATION`.

    Parameters
    ----------
    ca_coords
        ``(n_residues, 3)`` alpha-carbon coordinates for one contiguous chain segment,
        in residue order. At least two residues.
    sequence
        One-letter sequence of the same residues. Used for residue naming and to decide
        which residues have a ``CB``; glycine and ``UNK`` do not.
    rng
        Seeded generator. Used for exactly one thing: a random phase offset on the
        peptide-plane rotation grid. The solver itself is deterministic dynamic
        programming, so the same ``rng`` state gives bit-identical output, and a
        different seed explores a slightly different discretization -- which is honest
        about the fact that a CA trace does not determine its backbone uniquely.
    sidechains
        Place ``CB`` in addition to the backbone. **Only ``CB``**: nothing beyond it.
        That is a deliberate refusal to ship a fabricated rotamer library, and CB on its
        own is what fixes cartoon and licorice rendering. See
        :func:`place_sidechain_cb`. Note that ``CB`` is *modelled* either way -- it
        constrains the search and is held to the contact limits regardless -- so this
        flag changes what is emitted, not where the backbone goes.
    allow_chain_breaks
        Treat a CA-CA separation outside :func:`reconstructable_ca_bond_range` as a chain
        break and build each intact stretch separately, instead of raising. The breaks
        are listed in :attr:`BackboneResult.chain_breaks`; no peptide bond is built
        across one, so the ``C``-to-``N`` distance there is not a bond -- but it is still
        required to respect the non-bonded C-N limit of
        ``2.90 A``, because atoms that are not bonded cannot overlap either. Off by
        default.
    n_ca_c_tolerance
        Largest permitted deviation of the reconstructed ``N-CA-C`` angle from
        :data:`~dodo.constants.N_CA_C_ANGLE`, in degrees. Must be finite and positive.
        The solver treats it as a constraint, not a preference, so a trace that admits a
        conforming backbone gets one.

    Returns
    -------
    BackboneResult
        Flat atom records plus the measured geometry.

    Raises
    ------
    GeometryError
        If the input is malformed or non-finite; if fewer than two residues are given; if
        a CA-CA separation is outside what a trans peptide unit can span and
        ``allow_chain_breaks`` is not set; if the reconstructed ``N-CA-C`` angle misses by
        more than ``n_ca_c_tolerance``; or if any non-bonded heavy-atom pair ends up
        inside its hard-sphere limit. The ``N-CA-C`` failure is usually not a bug here: a
        CA pseudo-angle above ``max_reconstructable_ca_angle(n_ca_c_tolerance=...)`` makes
        it geometrically impossible, and :func:`unreconstructable_ca_angles` names those
        residues without building anything.

    Examples
    --------
    >>> import numpy as np
    >>> from dodo.constants import CA_CA_BOND_LENGTH
    >>> ca = np.zeros((4, 3))
    >>> ca[:, 0] = np.arange(4) * CA_CA_BOND_LENGTH * 0.8
    >>> ca[1::2, 1] = CA_CA_BOND_LENGTH * 0.6
    >>> result = place_backbone(ca, "AAAA", np.random.default_rng(0))
    >>> bool(np.array_equal(result.xyz[result.atom_names == "CA"], ca))
    True
    """
    ca, cleaned = _validate_inputs(ca_coords, sequence)
    tolerance = _validate_tolerance(n_ca_c_tolerance)
    return _reconstruct(
        ca,
        cleaned,
        rng,
        sidechains,
        allow_chain_breaks=allow_chain_breaks,
        n_ca_c_tolerance=tolerance,
    )


def _reconstruct(
    ca: np.ndarray,
    sequence: str,
    rng: np.random.Generator,
    sidechains: bool,
    *,
    allow_chain_breaks: bool,
    n_ca_c_tolerance: float,
    fixed_first_n: np.ndarray | None = None,
    fixed_last_c: np.ndarray | None = None,
    fixed_last_o: np.ndarray | None = None,
) -> BackboneResult:
    """Shared body of :func:`place_backbone` and :func:`place_backbone_for_domain`.

    The ``fixed_*`` atoms belong to peptide units outside the trace being solved. They
    are used by the objective and copied through, never re-placed. Only an anchored
    in-place rebuild has them, and it is what they exist for: without them the anchor
    residue's ``C`` is solved against a placeholder ``N`` and the seam peptide bond comes
    out at 1.38 A instead of 1.329.
    """
    n_residues = ca.shape[0]

    # One phase offset for the whole trace rather than one per unit: the grid step is the
    # only thing randomness touches here, and a per-unit offset would make the reported
    # geometry vary between adjacent units for no physical reason.
    phase = float(rng.uniform(0.0, 360.0 / _PLANE_GRID))

    residue_names = np.array([ONE_TO_THREE.get(letter, "UNK") for letter in sequence], dtype="<U3")
    bounds = _segment_bounds(ca, allow_chain_breaks)
    breaks = tuple(stop - 1 for _, stop in bounds[:-1])
    segments = [
        _segment_geometry(
            ca[start:stop],
            sequence[start:stop],
            start,
            phase,
            fixed_first_n=fixed_first_n if start == 0 else None,
            fixed_last_c=fixed_last_c if stop == n_residues else None,
            fixed_last_o=fixed_last_o if stop == n_residues else None,
        )
        for start, stop in bounds
    ]

    model = _repair_clashes(ca, residue_names, segments, n_ca_c_tolerance, breaks)

    strain = np.full(n_residues - 1, np.nan)
    for segment in segments:
        strain[segment.start : segment.start + segment.n_units] = segment.strain

    tau = _angle_between(model.nitrogen - ca, model.carbon - ca)
    _raise_on_n_ca_c(tau, ca, n_ca_c_tolerance)

    table = _atom_table(ca, model.nitrogen, model.carbon, model.oxygen, model.beta_carbon)
    clashes, closest = _steric_survey(table, breaks, 2.0 * _MAX_CONTACT)
    if clashes:
        _raise_on_clashes(clashes)

    xyz_rows: list[np.ndarray] = []
    names: list[str] = []
    elements: list[str] = []
    indices: list[int] = []
    for residue in range(n_residues):
        placed = {
            "N": model.nitrogen[residue],
            "CA": ca[residue],
            "C": model.carbon[residue],
            "O": model.oxygen[residue],
        }
        order = list(_BACKBONE_ORDER)
        if sidechains and str(residue_names[residue]) not in _NO_CB_RESIDUES:
            placed["CB"] = model.beta_carbon[residue]
            order.append("CB")
        for name in order:
            xyz_rows.append(placed[name])
            names.append(name)
            elements.append(_ELEMENT_OF[name])
            indices.append(residue)

    return BackboneResult(
        xyz=np.ascontiguousarray(np.array(xyz_rows, dtype=np.float64)),
        atom_names=np.array(names, dtype="<U4"),
        elements=np.array(elements, dtype="<U2"),
        residue_index=np.array(indices, dtype=np.int64),
        residue_names=residue_names,
        n_ca_c_angles=tau,
        peptide_angle_strain=strain,
        chain_breaks=breaks,
        min_nonbonded_contact=closest,
        steric_clashes=clashes,
        sidechains=sidechains,
    )


def _raise_on_n_ca_c(tau: np.ndarray, ca: np.ndarray, n_ca_c_tolerance: float) -> None:
    """Refuse a reconstruction whose ``N-CA-C`` angles miss, and say why."""
    n_residues = int(ca.shape[0])
    deviation = np.abs(tau - N_CA_C_ANGLE)
    over = np.flatnonzero(deviation > n_ca_c_tolerance)
    if not over.size:
        return
    ceiling = max_reconstructable_ca_angle(n_ca_c_tolerance=n_ca_c_tolerance)
    pseudo = _ca_pseudo_angles(ca)
    preview = ", ".join(
        f"residue {int(i)}: N-CA-C {tau[i]:.1f} deg, CA pseudo-angle {pseudo[i]:.1f} deg"
        if 0 < i < n_residues - 1
        else f"residue {int(i)}: N-CA-C {tau[i]:.1f} deg"
        for i in over[:5]
    )
    more = f" (and {over.size - 5} more)" if over.size > 5 else ""
    raise GeometryError(
        f"{over.size} of {n_residues} residue(s) ended up with an N-CA-C angle more "
        f"than {n_ca_c_tolerance:.1f} deg from {N_CA_C_ANGLE} deg: {preview}{more}. No "
        f"trans backbone can realize a CA pseudo-angle above {ceiling:.1f} deg at this "
        f"tolerance (see max_reconstructable_ca_angle), and adjacent residues near that "
        f"ceiling compete for the same peptide unit, so a trace sampled up to "
        f"BACKBONE_ANGLE_MAX is not always reconstructable. Call "
        f"unreconstructable_ca_angles(ca_coords) to list the residues that cannot work, "
        f"and regenerate the trace with a CA-CA-CA angle window inside the ceiling."
    )


def _raise_on_clashes(clashes: tuple[StericClash, ...]) -> None:
    """Refuse a reconstruction with overlapping atoms, and say what can be done.

    Two failure modes get different advice, because they have different causes. If both
    atoms of the worst clash are alpha carbons, the *input trace* intersects itself and no
    placement of N, C and O can help -- that is a defect in whatever generated the trace,
    and :func:`~dodo.geometry.metrics.validate_ca_trace` finds it directly. Otherwise the
    solver could not find a set of peptide-plane rotations that avoids the contact.
    """
    worst = clashes[0]
    preview = "; ".join(str(clash) for clash in clashes[:5])
    more = f" (and {len(clashes) - 5} more)" if len(clashes) > 5 else ""
    ca_only = all(clash.atom_i == "CA" and clash.atom_j == "CA" for clash in clashes)
    cause = (
        "The clashing atoms are alpha carbons of the input trace, so no backbone exists "
        "for this trace at all: validate_ca_trace() in dodo.geometry.metrics reports the "
        "same thing on the CA coordinates alone."
        if ca_only
        else "No assignment of peptide-plane rotations avoids these contacts. This is a "
        "property of the CA trace: two alpha carbons close enough to each other leave no "
        "room for the atoms hanging off them, even though every CA-CA distance is legal."
    )
    raise GeometryError(
        f"{len(clashes)} non-bonded heavy-atom pair(s) are closer than their hard-sphere "
        f"contact limit, worst {worst.overlap:.3f} A inside it: {preview}{more}. "
        f"{cause} Returning this would report a backbone with overlapping atoms as "
        f"valid, which is worse than refusing: a rebuilt backbone whose atoms are "
        f"{worst.distance:.2f} A apart is not usable output."
    )


def _ca_pseudo_angles(ca: np.ndarray) -> np.ndarray:
    """``CA(i-1)-CA(i)-CA(i+1)`` angles in degrees, NaN at both termini.

    Reported in error messages, because the pseudo-angle is what actually explains an
    unreachable ``N-CA-C``, and a message that names the symptom without the cause sends
    the reader looking in this module for a bug that is upstream.
    """
    angles = np.full(ca.shape[0], np.nan)
    if ca.shape[0] >= 3:
        angles[1:-1] = _angle_between(ca[:-2] - ca[1:-1], ca[2:] - ca[1:-1])
    return angles


def _anchor_atom(structure: Structure, residue: int, name: str) -> np.ndarray:
    """One named atom of an anchor residue, or raise.

    An anchor's atoms are inputs to the solve, not outputs of it: the peptide unit joining
    the anchor to the span is built onto them. A missing one is therefore not something to
    work around -- substituting an arbitrary dihedral for it is exactly what left the seam
    peptide bonds at 1.38 and 1.27 A.
    """
    atom_slice = structure.atom_slice_for_residues(residue, residue + 1)
    found = np.flatnonzero(structure.atom_name[atom_slice] == name)
    if found.size == 0:
        raise GeometryError(
            f"Anchor residue {structure.residue_label(residue)} has no {name} atom. The "
            f"anchor's own N, C and O are what the rebuilt span's terminal peptide unit "
            f"closes onto, so an anchor with missing backbone atoms cannot anchor "
            f"anything. Move the anchor, or rebuild the span without one."
        )
    return np.asarray(structure.xyz[atom_slice.start + int(found[0])], dtype=np.float64)


def place_backbone_for_domain(
    structure: Structure,
    domain: Domain,
    rng: np.random.Generator,
    sidechains: bool = False,
    *,
    allow_chain_breaks: bool = False,
    n_ca_c_tolerance: float = _N_CA_C_TOLERANCE,
) -> BackboneResult:
    """Rebuild one domain's ``N``, ``C``, ``O`` -- and optionally ``CB`` -- in place.

    The domain's *anchor* residues are included in the solve. That is what
    :attr:`~dodo.structure.Span.n_anchor` and :attr:`~dodo.structure.Span.c_anchor` are
    for: an anchored domain's first and last residues are interior residues of the chain,
    and treating them as chain termini would substitute the arbitrary
    :data:`_TERMINAL_PHI` / :data:`_TERMINAL_PSI` for geometry the neighbouring residue
    actually determines.

    Which rows get written is not the same question as which rows get solved, and getting
    that wrong is how both seam peptide bonds came out broken. A peptide unit spans *two*
    residues: the unit joining the N-side anchor to the span owns ``C`` and ``O`` of the
    anchor and ``N`` of the span's first residue. Writing only the span's rows pairs a
    freshly solved ``N`` with the anchor's original ``C``, and the bond between them is
    then whatever the two happen to be -- measured on dnmt3a, 1.3767 A at one seam and
    1.2672 A with omega 162.7 deg at the other, while all 198 interior bonds were
    1.3290000 A at omega exactly 180. So the write set is the span plus ``C``/``O`` of the
    N-side anchor and ``N`` of the C-side anchor. Each anchor keeps the atoms that face
    away from the span, which belong to units outside it, and each anchor's own
    ``N-CA-C`` is solved against its real fixed atom rather than a placeholder.

    Only coordinates are written. This function cannot *add* atoms, because a
    :class:`~dodo.structure.Structure` owns one fixed-size coordinate array and this
    module does not own :mod:`dodo.structure`. So a residue that has no ``N``, ``C`` or
    ``O`` slot is an error, not something to skip: skipping would leave the caller with a
    structure that is half rebuilt and looks finished.

    For the CA-only case -- which is DODO's own IDR output -- call :func:`place_backbone`
    and construct a new :class:`~dodo.structure.Structure` from
    :meth:`BackboneResult.atom_records`.

    Parameters
    ----------
    structure
        The structure to write into. Mutated, and only if everything can be written.
    domain
        The domain whose residues are rebuilt. An interior span must carry an anchor on
        each side that has a neighbour; see Raises.
    rng
        Seeded generator; see :func:`place_backbone`.
    sidechains
        Also rewrite ``CB`` where the residue has one and the slot exists. The anchors'
        ``CB`` is rewritten too when present, because ``CB`` is placed from ``N`` and
        ``C`` and one of those has moved.
    allow_chain_breaks, n_ca_c_tolerance
        See :func:`place_backbone`.

    Returns
    -------
    BackboneResult
        The solved geometry, covering the span *and its anchors*, so that the caller can
        inspect what was written -- ``N-CA-C`` angles, peptide strain, the closest
        non-bonded contact -- instead of having to re-measure the mutated structure by
        hand. Returning ``None`` here meant the broken seam bonds could not be detected
        at all without doing that.

    Raises
    ------
    GeometryError
        If the domain is shorter than two residues; if the span has a neighbouring residue
        on a side where it declares no anchor, since the peptide bond to that neighbour
        cannot then be closed; if any residue in the span, or either anchor, lacks an
        ``N``, ``C`` or ``O`` atom; or for any reason :func:`place_backbone` raises.
        Nothing is written unless everything can be.
    """
    span = domain.span
    if len(span) < 2:
        raise GeometryError(
            f"Cannot rebuild the backbone of {domain!r}: it has {len(span)} residue(s) "
            f"and a peptide unit needs two."
        )

    # Extend over the anchors so the domain's own terminal residues are solved as
    # interior residues.
    start = span.n_anchor if span.n_anchor is not None else span.start
    stop = (span.c_anchor + 1) if span.c_anchor is not None else span.stop
    if not (0 <= start < stop <= structure.n_residues):
        raise GeometryError(
            f"{domain!r} with its anchors spans residues [{start}, {stop}), which is "
            f"outside the structure's {structure.n_residues} residues."
        )
    unanchored: list[tuple[str, int]] = []
    if span.n_anchor is None and span.start > 0:
        unanchored.append(("n_anchor", span.start - 1))
    if span.c_anchor is None and span.stop < structure.n_residues:
        unanchored.append(("c_anchor", span.stop))
    if unanchored:
        detail = ", ".join(
            f"{side} should be {structure.residue_label(neighbour)}"
            for side, neighbour in unanchored
        )
        raise GeometryError(
            f"{domain!r} has a neighbouring residue but no anchor on that side "
            f"({detail}), so the peptide bond across the boundary cannot be closed: a "
            f"freshly solved N or C would be paired with an atom this call is not allowed "
            f"to move, leaving a bond of whatever length falls out -- 1.38 A and 1.27 A "
            f"when this was measured. Set the anchor, or rebuild a span that reaches the "
            f"end of the chain."
        )

    fixed_first_n = _anchor_atom(structure, start, "N") if span.n_anchor is not None else None
    fixed_last_c = _anchor_atom(structure, stop - 1, "C") if span.c_anchor is not None else None
    fixed_last_o = _anchor_atom(structure, stop - 1, "O") if span.c_anchor is not None else None

    solved = _reconstruct(
        np.ascontiguousarray(structure.ca_xyz[start:stop], dtype=np.float64),
        structure.sequence[start:stop].strip().upper(),
        rng,
        sidechains,
        allow_chain_breaks=allow_chain_breaks,
        n_ca_c_tolerance=_validate_tolerance(n_ca_c_tolerance),
        fixed_first_n=fixed_first_n,
        fixed_last_c=fixed_last_c,
        fixed_last_o=fixed_last_o,
    )

    # The write set: every atom of the span, plus the anchor atoms that belong to the
    # peptide units joining the anchors to the span. CB follows N and C, so an anchor
    # whose C or N moved needs its CB recomputed too.
    wanted: dict[int, list[str]] = {}
    for residue in range(span.start, span.stop):
        wanted[residue] = ["N", "C", "O", "CB"] if sidechains else ["N", "C", "O"]
    if span.n_anchor is not None:
        wanted[start] = ["C", "O", "CB"] if sidechains else ["C", "O"]
    if span.c_anchor is not None:
        wanted[stop - 1] = ["N", "CB"] if sidechains else ["N"]

    # Resolve every target slot before writing anything, so a missing atom cannot leave
    # the structure half updated.
    targets: dict[tuple[int, str], int] = {}
    missing: list[str] = []
    for residue, names_wanted in sorted(wanted.items()):
        atom_slice = structure.atom_slice_for_residues(residue, residue + 1)
        names = structure.atom_name[atom_slice]
        residue_name = str(structure.residue_name[residue])
        for name in names_wanted:
            if name == "CB" and residue_name in _NO_CB_RESIDUES:
                continue
            found = np.flatnonzero(names == name)
            if found.size == 0:
                missing.append(f"{structure.residue_label(residue)}:{name}")
                continue
            targets[(residue, name)] = atom_slice.start + int(found[0])
    if missing:
        preview = ", ".join(missing[:5])
        more = f" (and {len(missing) - 5} more)" if len(missing) > 5 else ""
        raise GeometryError(
            f"{len(missing)} atom slot(s) needed by the rebuilt backbone do not exist in "
            f"this structure: {preview}{more}. place_backbone_for_domain rewrites "
            f"coordinates and cannot add atoms. For a CA-only region, call "
            f"place_backbone and build a new Structure from "
            f"BackboneResult.atom_records()."
        )

    before = _structure_contacts(structure)
    saved = {index: structure.xyz[index].copy() for index in targets.values()}
    for (residue, name), atom_index in targets.items():
        source = np.flatnonzero(
            (solved.residue_index == residue - start) & (solved.atom_names == name)
        )
        if source.size == 0:
            # The residue has the slot but the solved model has no such atom: the
            # structure's residue_name and its one-letter sequence disagree about whether
            # this residue has a CB. Leaving the existing coordinate alone is the only
            # honest option -- there is nothing to write.
            continue
        structure.xyz[atom_index] = solved.xyz[int(source[0])]

    # The solve only sees the residues it was given, so nothing in it knows about the rest
    # of the structure. That makes this check necessary rather than belt-and-braces: a span
    # can be rebuilt into a perfectly self-consistent backbone that runs through the
    # neighbouring domain. Only contacts this call *created or tightened* count -- a
    # structure's pre-existing overlaps are not this function's to adjudicate -- and if any
    # exist the write is undone, so a refused rebuild leaves the structure as it was.
    introduced = {
        pair: distance
        for pair, distance in _structure_contacts(structure).items()
        if (pair & set(targets.values())) and distance < before.get(pair, np.inf) - 1e-9
    }
    if introduced:
        for atom_index, coordinates in saved.items():
            structure.xyz[atom_index] = coordinates
        worst = min(introduced, key=lambda pair: introduced[pair])
        first, second = sorted(worst)
        preview = ", ".join(
            f"{structure.atom_name[i]}({structure.residue_label(int(structure.residue_index[i]))})"
            f"..{structure.atom_name[j]}"
            f"({structure.residue_label(int(structure.residue_index[j]))}) "
            f"at {introduced[pair]:.3f} A"
            for pair in sorted(introduced, key=lambda pair: introduced[pair])[:5]
            for i, j in [sorted(pair)]
        )
        raise GeometryError(
            f"The rebuilt span would put {len(introduced)} heavy-atom pair(s) inside their "
            f"hard-sphere contact limit against the rest of the structure, closest "
            f"{introduced[worst]:.3f} A ({structure.atom_name[first]}/"
            f"{structure.atom_name[second]}): {preview}. Nothing was written. The solve is "
            f"given only the span and its anchors, so it cannot avoid atoms outside them; "
            f"a different rng seed explores a different set of peptide-plane rotations and "
            f"may fit, but a span that has to thread through occupied space needs the "
            f"domain moved, not the backbone re-solved."
        )
    return solved


def _structure_contacts(structure: Structure) -> dict[frozenset[int], float]:
    """Every non-bonded heavy-atom pair in a structure that is inside its contact limit.

    Keyed by atom index so that two surveys of the same structure can be compared pair by
    pair. Only ``N``, ``CA``, ``C``, ``O`` and ``CB`` are considered, because those are the
    atoms this module places and the ones :data:`_MIN_CONTACT` is calibrated on.
    """
    keep = np.isin(structure.atom_name, list(_ELEMENT_OF))
    indices = np.flatnonzero(keep)
    if indices.size == 0:
        return {}
    codes = np.array(
        [_ELEMENT_CODE[_ELEMENT_OF[str(name)]] for name in structure.atom_name[indices]]
    )
    table = _AtomTable(
        xyz=np.ascontiguousarray(structure.xyz[indices], dtype=np.float64),
        residue=structure.residue_index[indices],
        name=structure.atom_name[indices],
        code=codes,
    )
    clashes, _ = _steric_survey(table, (), _MAX_CONTACT)
    found: dict[frozenset[int], float] = {}
    for clash in clashes:
        first = int(
            indices[
                np.flatnonzero((table.residue == clash.residue_i) & (table.name == clash.atom_i))[0]
            ]
        )
        second = int(
            indices[
                np.flatnonzero((table.residue == clash.residue_j) & (table.name == clash.atom_j))[0]
            ]
        )
        found[frozenset((first, second))] = clash.distance
    return found
