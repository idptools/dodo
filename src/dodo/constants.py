"""Single source of truth for every physical and algorithmic constant in DODO.

Nothing in this package should hardcode a distance, angle, or threshold. v1 and the
first v2 attempt both accumulated several divergent copies of the same numbers -- at
one point three different values of ``super_compact`` and three different CA-CA bond
lengths were live simultaneously -- and reconciling them after the fact was the
single most expensive part of the rewrite. Everything lands here instead, with its
provenance recorded, so a future reader can tell a measured value from a tuned knob
from a guess.

Provenance tags used in comments below:

MEASURED
    Derived from structural data; the source is named.
TUNED
    Fit empirically against AF2 structures by the original author. The tuning dataset was
    not preserved, so these should not be changed without re-tuning.
DERIVED
    Computed from other constants in this file.
CHOICE
    An arbitrary but deliberate engineering decision.
"""

from __future__ import annotations

from typing import Final

import numpy as np

# ---------------------------------------------------------------------------
# Backbone geometry
# ---------------------------------------------------------------------------

#: Virtual CA(i)-CA(i+1) distance in Angstroms.
#:
#: MEASURED. The trans-peptide CA-CA virtual bond length is 3.80-3.81 A and is
#: remarkably rigid (cis-proline, at ~2.9 A, is the only common exception and DODO
#: does not model it).
#:
#: This value was inconsistent across the pre-rewrite code: v1's ``parameters.py``
#: said 3.8, v2's cone generator used 3.856, and one distance table used 3.89. There is
#: now exactly one, on Ryan's instruction: 3.81. It applies everywhere, including as the
#: target that generative-model output is projected onto -- see
#: :mod:`dodo.geometry.regularize`, since a diffusion model does not produce exact bonds.
#:
#: Do not introduce a second value for any specific consumer. Reconciling the three that
#: were live at once was among the most expensive parts of this rewrite.
CA_CA_BOND_LENGTH: Final[float] = 3.81

#: Acceptable spread on the CA-CA bond when closing onto a fixed anchor, in A.
#: CHOICE. Wide enough that closure is achievable, tight enough that a viewer still
#: draws the bond.
CA_CA_BOND_TOLERANCE: Final[float] = 0.10

#: Minimum non-bonded CA-CA approach distance in Angstroms.
#:
#: TUNED. v1's ``parameters.py`` established 3.2 and that is what we keep. It is
#: permissive: two CAs in a real packed core rarely come closer than ~4.5 A. But
#: DODO builds compact coils into cavities left by folded domains, and a 4.5 A
#: exclusion makes many closures geometrically impossible. 3.2 is the compromise.
CA_CLASH_DISTANCE: Final[float] = 3.20

#: Relaxation ladder applied when no candidate satisfies CA_CLASH_DISTANCE.
#:
#: CHOICE. Tried in order; the first value yielding a non-empty candidate set wins.
#: Accepting a mild clash is strictly better than aborting the build, but the caller
#: is told which rung was used so it can be reported rather than hidden.
CLASH_RELAXATION_LADDER: Final[tuple[float, ...]] = (3.20, 2.80, 2.50, 2.00)

#: Residue separation below which a CA pair is exempt from clash checking.
#:
#: DERIVED from chain connectivity: ``|i-j| == 1`` are covalently bonded at 3.8 A and
#: ``|i-j| == 2`` are geometrically constrained by the backbone angle to 5.0-7.5 A.
#: Neither is a clash. The pre-rewrite code had no such exclusion, which is why
#: its whole-structure clash check reported every peptide bond as a violation.
CLASH_EXCLUDE_WITHIN_RESIDUES: Final[int] = 2

# ---------------------------------------------------------------------------
# CA-CA-CA pseudo-angle model
# ---------------------------------------------------------------------------
#
# MEASURED from AlphaFold2 structures by the original author: the CA(i-1)-CA(i)-CA(i+1)
# pseudo-angle has mean ~125 deg and standard deviation ~26 deg, with an observed
# range of 75-179 deg.
#
# Restricting generated angles to a window is what keeps a backbone from folding back
# on itself: the pre-rewrite IDR builder had no angle constraint at all and produced
# measured angles as sharp as 47 deg, which no real trace exhibits and which cannot be
# reconstructed to all-atom.

BACKBONE_ANGLE_MEAN: Final[float] = 125.0  # MEASURED
BACKBONE_ANGLE_SD: Final[float] = 26.0  # MEASURED
BACKBONE_ANGLE_OBSERVED_MIN: Final[float] = 75.0  # MEASURED
BACKBONE_ANGLE_OBSERVED_MAX: Final[float] = 179.0  # MEASURED

#: Sampling window for generated backbones, in degrees.
#:
#: TUNED against AlphaFold2 structures by the package author, and this is the value that
#: stands. Relative to the measured distribution it is mean -1.31 sd to mean +1.38 sd, so it
#: covers roughly the central 80% of observed angles -- tighter than the full observed 75-179
#: range, to keep sampling away from the rare extremes while still spanning helix-like (~91
#: deg) through extended (~161 deg) geometry.
#:
#: It is *not* symmetric about the mean and not a round number of standard deviations, so do
#: not "tidy" it into mean +/- 1 sd: that would narrow it to 99-151 and exclude helical
#: geometry.
#:
#: DO NOT NARROW THIS TO SERVE ALL-ATOM RECONSTRUCTION. It was briefly capped at 150 for that
#: reason and reverted. A CA pseudo-angle is coupled to N-CA-C for a trans peptide, so wide
#: angles are harder to reconstruct -- measured, an isolated angle reconstructs up to about
#: 154.6 deg at the backbone module's tau tolerance. But getting the CA-only geometry right is
#: the first priority, and capping the window degrades it (more compact chains than the
#: measured distribution supports) to benefit a later one. Capping magnitude did not even fix
#: all-atom acceptance, because the real constraint is on the *change* in pseudo-angle between
#: consecutive residues: measured allowance is ~35 deg at a base angle of 95 and collapses to
#: ~4 deg at 150. That belongs in the all-atom path as a joint constraint on consecutive
#: angles, not here as a magnitude cap.
BACKBONE_ANGLE_MIN: Final[float] = 91.0
BACKBONE_ANGLE_MAX: Final[float] = 161.0

#: Preferred angle. Candidates are ordered by ``|angle - ideal|`` so that a
#: first-non-clashing search naturally prefers realistic geometry.
BACKBONE_ANGLE_IDEAL: Final[float] = BACKBONE_ANGLE_MEAN

#: Candidate positions generated per angle in the cone template. CHOICE: 10 x 71
#: angles = 710 candidates per step, which empirically finds a non-clashing option
#: on the first pass in the overwhelming majority of cases.
CANDIDATES_PER_ANGLE: Final[int] = 10

# ---------------------------------------------------------------------------
# All-atom backbone geometry
# ---------------------------------------------------------------------------
#
# MEASURED (standard peptide geometry). v1 declared correct values in parameters.py
# and then never imported them: its all-atom module added unit vectors without
# scaling by bond length, producing CA-C at exactly 1.000 A and peptide bonds at
# 2.87 A. These are the real numbers and they are used.

N_CA_BOND_LENGTH: Final[float] = 1.458  # MEASURED
CA_C_BOND_LENGTH: Final[float] = 1.525  # MEASURED
C_O_BOND_LENGTH: Final[float] = 1.231  # MEASURED
C_N_PEPTIDE_BOND_LENGTH: Final[float] = 1.329  # MEASURED

N_CA_C_ANGLE: Final[float] = 111.0  # MEASURED, degrees
CA_C_N_ANGLE: Final[float] = 116.2  # MEASURED, degrees
C_N_CA_ANGLE: Final[float] = 121.7  # MEASURED, degrees
CA_C_O_ANGLE: Final[float] = 120.8  # MEASURED, degrees

#: Backbone omega dihedral for a trans peptide, in degrees. MEASURED.
OMEGA_TRANS: Final[float] = 180.0

# ---------------------------------------------------------------------------
# Region identification
# ---------------------------------------------------------------------------

#: All-atom contact radius in Angstroms for the primary folded/disordered score. TUNED.
#:
#: This is the author's original value and it belongs to :func:`~dodo.regions.contact
#: .atom_pair_counts`, the density metric DODO was built and validated on. Larger values start
#: merging folded domains that AlphaFold happens to park close together in space, which is the
#: specific failure this was tuned to avoid.
CONTACT_RADIUS: Final[float] = 8.0

#: Folded/disordered cutoff on the primary density score: number of all-atom PAIRS within
#: :data:`CONTACT_RADIUS`. TUNED, and the author reports it outperforms sequence-based disorder
#: predictors at drawing region boundaries.
#:
#: Known caveat, recorded rather than acted on: a pair count scales with the residue's own
#: heavy-atom count, so it is composition-sensitive. Measured within one folded domain, the
#: correlation with heavy-atom count is r = 0.65, every glycine falls below this threshold
#: (mean 292), and 94% of Trp/Phe/Tyr sit above it (mean 943). The smoothing and merging stages
#: downstream absorb a good deal of that noise, and the metric demonstrably works, so it stays
#: the default. :data:`CA_CONTACT_RADIUS` below is the composition-free alternative for
#: comparison, not a replacement.
CONTACT_SCORE_THRESHOLD: Final[float] = 480.0

#: Radius in Angstroms for the ALTERNATIVE CA-only folded/disordered score.
#:
#: MEASURED. Chosen by sweeping 8-16 A against arf19, dnmt3a and p300 and scoring how well the
#: counts separate known folded domains from known IDRs: balanced accuracy 86/99/90% at 8 A,
#: 96/100/95% at 12 A, 99/100/96% at 14 A. Accuracy keeps climbing past 12, but so does the
#: risk of merging nearby domains, and the per-structure optimal threshold starts diverging
#: (at 14 A the three want 15, 6 and 20; at 12 A the two informative ones agree on 10 and 11).
CA_CONTACT_RADIUS: Final[float] = 12.0

#: Cutoff for the ALTERNATIVE CA-only score: number of non-local CA neighbours within
#: :data:`CA_CONTACT_RADIUS`.
#:
#: MEASURED, the lower of the two per-structure optima (10 and 11). The lower value is chosen
#: deliberately because the two errors are not symmetric: calling a folded domain disordered
#: replaces real structure with a random walk, while calling an IDR folded merely leaves
#: AlphaFold's own coordinates in place. So err toward folded.
#:
#: This score counts residues by their alpha carbons, so composition cannot bias it -- every
#: residue has exactly one CA. It is also scale-invariant: measured on arf19, the all-atom
#: density score of the same structure stripped to CA-only came out at 0.26x its full value,
#: whereas this score is identical either way (measured ratio 1.000). That matters for input
#: with unmodelled side chains and for DODO's own CA-only output.
CA_CONTACT_SCORE_THRESHOLD: Final[float] = 10.0

#: CA-only contact radius in Angstroms for LOOP detection. TUNED.
#:
#: Deliberately tighter than either folded/disordered radius: a loop is about *local* backbone
#: packing, not about how much of the domain happens to be nearby. CA-only because a loop is
#: defined by where the backbone goes -- the side chains of a floppy loop can make plenty of
#: contacts without the backbone being packed at all.
LOOP_CONTACT_RADIUS: Final[float] = 7.0

#: Smoothing window (residues) applied to the contact score before thresholding.
#: CHOICE. The raw per-residue score is noisy enough to fragment a single domain
#: into dozens of blocks, which is what masked v1's two domain-merge bugs.
CONTACT_SCORE_SMOOTHING_WINDOW: Final[int] = 7

#: A residue is "loop-like" if fewer than this many other CAs lie within
#: LOOP_CONTACT_RADIUS. TUNED. Physical basis: a packed-core CA has ~6-10 CA
#: neighbours within 7 A; an extended loop has 2-4.
LOOP_CONTACT_CUTOFF: Final[int] = 6

#: Minimum length in residues of a rebuildable loop inside a folded domain.
#:
#: TUNED at 10, but note the pre-rewrite code compared with strict ``>``, making the
#: effective minimum 11 while the docs said 10. We compare with ``>=`` so the
#: constant means what it says. That is a deliberate one-residue behaviour change.
MIN_LOOP_LENGTH: Final[int] = 10

#: Longest run of low-contact residues that stays *inside* one folded domain.
#:
#: TUNED at 25. The pre-rewrite code used a single ``gap_thresh`` knob for this AND
#: for the minimum acceptable folded-domain length, which are unrelated quantities.
#: They are split here.
MAX_INTERNAL_GAP: Final[int] = 25

#: Shortest run of folded residues that counts as a folded domain. Split from
#: MAX_INTERNAL_GAP, above; TUNED to the same starting value.
MIN_FOLDED_DOMAIN_LENGTH: Final[int] = 25

#: Minimum number of consecutive above-threshold residues to seed a folded block.
MIN_FOLDED_SEED_RUN: Final[int] = 2

#: Shortest IDR worth rebuilding. Below this, leave the input coordinates alone --
#: there is no meaningful polymer statistics to impose on 3 residues.
MIN_IDR_LENGTH: Final[int] = 4

#: pLDDT below which an AF2 residue is treated as disordered. CHOICE, following the
#: widely used AF2 confidence bands (<50 very low, 50-70 low, 70-90 confident).
#: pLDDT is already in the B-factor column and the pre-rewrite code never looked at
#: it, re-deriving confidence geometrically from coordinates instead.
PLDDT_DISORDER_THRESHOLD: Final[float] = 70.0

# ---------------------------------------------------------------------------
# Chain dimension targets
# ---------------------------------------------------------------------------

#: Named build modes, as multipliers on the predicted end-to-end distance.
#:
#: This is the substantive scientific change in v2. v1 expressed these as Angstroms
#: *per residue*, i.e. linear in N, but real IDR end-to-end distance scales as
#: roughly N^0.52. A fixed A/residue multiplier can therefore only agree with the
#: prediction at one chain length: v1's "normal" (0.8 A/res) gives 80 A at N=100
#: (about right) and 400 A at N=500, where the prediction is closer to 200 A.
#:
#: Rebasing onto the prediction makes every knob length-independent: "expanded"
#: means 1.3x whatever this sequence's predicted dimension is, at any length.
#:
#: DERIVED from v1's multipliers by dividing through by v1's "normal" (0.8), then
#: rounded. Note this makes "normal" and "predicted" synonyms, which they were
#: not in v1 -- a deliberate, documented break.
MODES: Final[dict[str, float]] = {
    "super_compact": 0.4,
    "compact": 0.7,
    "normal": 1.0,
    "predicted": 1.0,
    "expanded": 1.3,
    "super_expanded": 1.6,
    "max_expansion": 2.0,
}

#: Default mode when the caller does not specify one.
DEFAULT_MODE: Final[str] = "predicted"

#: Analytical fallback for end-to-end distance when ALBATROSS is unavailable
#: (the dependency-light "lite" install), as Re = PREFACTOR * N ** EXPONENT in A.
#:
#: MEASURED, indirectly: Kohn et al. (2004) PNAS 101:12491 report Rg = 2.54 * N^0.522
#: for chemically denatured proteins, and Re = sqrt(6) * Rg for an ideal chain, giving
#: a prefactor of sqrt(6) * 2.54 = 6.22.
#:
#: HOW GOOD IS IT? Benchmarked against ALBATROSS over 72 sequences spanning six
#: compositional classes and N = 20-1000. Ratio of this estimate to the prediction:
#: mean 0.97, and per class at N = 500 --
#:
#:     polar          0.95x       polyampholyte  0.92x
#:     proline-rich   0.74x       polyanionic    0.62x
#:     polycationic   0.64x       hydrophobic    2.62x
#:
#: So for genuine IDR compositions it lands within roughly 0.6-0.95x, erring compact.
#: The charged classes are the worst of those because polyelectrolyte expansion pushes
#: their scaling exponent to 0.60-0.64, well above the 0.52 baked in here.
#:
#: The 2.62x outlier is instructive and is NOT a fitting failure: a sequence drawn
#: uniformly from all 20 amino acids is hydrophobic-rich, and ALBATROSS correctly
#: predicts a collapsed globule (measured scaling exponent 0.20, Re = 61 A at N = 500).
#: Such a sequence is not disordered at all, so no length-only law can describe it.
#:
#: The lesson: the irreducible error here is COMPOSITIONAL, not a matter of tuning. At
#: N = 500 the true Re spans 3.4x across compositions at fixed length, so refitting the
#: exponent cannot help -- it was tried, and a global refit (6.32 * N^0.537) came out
#: worse-centred (mean ratio 1.06). These values stay because they are citable and
#: better-centred.
#:
#: Use this to keep the dependency-light install functional and physically sane, not to
#: avoid installing sparrow for real work.
FLORY_RE_PREFACTOR: Final[float] = 6.22
FLORY_RE_EXPONENT: Final[float] = 0.522

#: Sequence length below which sparrow's ALBATROSS networks need their "scaled"
#: variant. Mirrors ``sparrow.data.configs.MIN_LENGTH_ALBATROSS_RE_RG``. Most loops
#: DODO rebuilds are shorter than this, so it matters more than it looks.
#: Re-read from sparrow at runtime when available rather than trusted blindly.
ALBATROSS_MIN_LENGTH: Final[int] = 35

# ---------------------------------------------------------------------------
# Conformation engines
# ---------------------------------------------------------------------------

#: Hard cap on sequence length for the STARLING engine, in residues.
#:
#: Observed against STARLING 2.0.2. Longer IDRs are handled by hierarchical segment
#: assembly (see engines/hierarchical.py) rather than by erroring out. Queried from
#: STARLING at runtime where it exposes the limit; this is the fallback.
STARLING_MAX_LENGTH: Final[int] = 380

#: Residues of overlap between adjacent STARLING segments in hierarchical assembly.
#: CHOICE. Splicing inside an overlap preserves locally correct backbone statistics
#: across the junction; butt-joining two independently generated segments does not.
SEGMENT_SPLICE_OVERLAP: Final[int] = 10

# ---------------------------------------------------------------------------
# Sampling budgets
# ---------------------------------------------------------------------------

#: Candidate positions evaluated per residue before the step is declared failed.
MAX_CANDIDATES_PER_RESIDUE: Final[int] = 1000

#: Restarts of a whole region before the build is declared failed.
MAX_ATTEMPTS_PER_REGION: Final[int] = 40

#: Candidate placements tried when positioning a folded domain.
MAX_FD_PLACEMENT_ATTEMPTS: Final[int] = 500

#: Conformers generated per batch in the vectorized walk.
WALK_BATCH_SIZE: Final[int] = 50

# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

#: CRYST1 box dimensions in Angstroms. Purely cosmetic -- DODO does not do periodic
#: boundaries -- but some viewers want the record present.
DEFAULT_BOX_DIMENSIONS: Final[tuple[float, float, float]] = (500.0, 500.0, 500.0)

#: Residues per SEQRES line, per the PDB format specification.
RESIDUES_PER_SEQRES_LINE: Final[int] = 13

#: B-factor values written when annotating regions for visualization.
#:
#: Note the polarity: folded = 100, disordered = 0, so that colouring by B-factor in
#: a viewer highlights the folded core. v1's docstrings claimed the opposite in three
#: places while its README and CLI help were correct; the code did this.
BETA_FOLDED: Final[float] = 100.0
BETA_DISORDERED: Final[float] = 0.0

#: Largest residue number representable in the PDB format's 4-column resSeq field.
MAX_PDB_RESIDUE_NUMBER: Final[int] = 9999

#: Largest atom serial representable in the PDB format's 5-column serial field.
#: Beyond this, real files switch to hybrid-36 encoding.
MAX_PDB_ATOM_SERIAL: Final[int] = 99999

# ---------------------------------------------------------------------------
# Element data
# ---------------------------------------------------------------------------

#: Atomic masses in unified atomic mass units, for centre-of-mass calculations.
#:
#: The pre-rewrite code carried two disagreeing copies of this table -- one if/elif
#: ladder and one dict, differing on C, O, S and P and on whether Se existed at all.
#: This is the only copy.
ATOMIC_MASSES: Final[dict[str, float]] = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "S": 32.06,
    "P": 30.974,
    "SE": 78.971,  # selenomethionine
}

#: Three-letter to one-letter residue codes, including the modified residues that
#: appear in real deposited structures as HETATM records.
THREE_TO_ONE: Final[dict[str, str]] = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    # Modified residues that are part of the polymer backbone and must not be
    # dropped. Discarding MSE in particular fabricates a phantom chain break,
    # which the pre-rewrite reader did.
    "MSE": "M",  # selenomethionine
    "SEC": "U",  # selenocysteine
    "PYL": "O",  # pyrrolysine
    "HYP": "P",  # hydroxyproline
    "SEP": "S",  # phosphoserine
    "TPO": "T",  # phosphothreonine
    "PTR": "Y",  # phosphotyrosine
    "MLY": "K",  # methyllysine
    "CSO": "C",  # oxidized cysteine
    "CME": "C",  # modified cysteine
    "UNK": "X",
}

ONE_TO_THREE: Final[dict[str, str]] = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
    "X": "UNK",
}

#: Backbone atom names, in the canonical order they appear within a residue.
BACKBONE_ATOMS: Final[tuple[str, ...]] = ("N", "CA", "C", "O")

#: Atoms of an *anchor* residue that a newly placed neighbouring CA may legitimately come
#: closer to than :data:`CA_CLASH_DISTANCE`, and which are therefore exempt from the obstacle
#: set when building the region that anchor holds.
#:
#: MEASURED over 649,658 sequence-neighbour CA / atom pairs sampled from the human proteome.
#: The 0.1th percentile of the distance from a residue's atoms to an adjacent residue's CA:
#:
#: ===== ====== ==================================================
#: atom  p0.1   why
#: ===== ====== ==================================================
#: C     2.410  1-3 through the peptide bond
#: N     2.379  1-3 through the peptide bond
#: O     2.712  1-4, but pinned nearly eclipsing the peptide plane
#: CA    3.018  the bonded partner itself
#: ===== ====== ==================================================
#:
#: Every *side-chain* atom stays above 3.29 A at the same percentile, with one exception:
#: proline CD reaches 2.778 A (minimum 2.245 A), because it is covalently bonded to proline's
#: own backbone N and so is 1-3 from the preceding residue's C. GLU/GLN/LYS/ARG CD -- the other
#: residues with a CD -- stay above 3.34 A, so the exemption is proline's alone.
#:
#: This exists because exempting the *whole* anchor residue is wrong. Doing so let the walk
#: place a CA on top of an anchor's side chain: measured collisions at 0.871 A (LEU CD1),
#: 0.937 A (ASN ND2) and 0.944 A (LYS CD) in output that every other check called clean.
ANCHOR_EXEMPT_ATOMS: Final[frozenset[str]] = frozenset(BACKBONE_ATOMS)

#: The anchor atoms that are exempt from clash checking **unconditionally**.
#:
#: A rebuilt region is bonded to its anchors' alpha carbons -- that bond is what attaches the
#: region to the structure at all -- so treating them as obstacles would make every valid
#: attachment a clash and leave nothing to build. This exemption is not a trade-off; there is no
#: version of the algorithm without it.
#:
#: :data:`ANCHOR_EXEMPT_ATOMS` above is the *discretionary* part, granted only when a region
#: cannot otherwise be built. See ``_obstacles_for_span``.
ANCHOR_ALWAYS_EXEMPT_ATOMS: Final[frozenset[str]] = frozenset({"CA"})

#: Per-residue additions to :data:`ANCHOR_EXEMPT_ATOMS`. See above for the measurement.
ANCHOR_EXEMPT_ATOMS_BY_RESIDUE: Final[dict[str, frozenset[str]]] = {
    "PRO": frozenset({"CD"}),
}

#: Heavy-atom counts per residue, used to normalize the contact score. DERIVED by
#: counting non-hydrogen atoms in each standard residue.
HEAVY_ATOM_COUNTS: Final[dict[str, int]] = {
    "GLY": 4,
    "ALA": 5,
    "SER": 6,
    "CYS": 6,
    "PRO": 7,
    "VAL": 7,
    "THR": 7,
    "MET": 8,
    "ASN": 8,
    "LEU": 8,
    "ILE": 8,
    "ASP": 8,
    "LYS": 9,
    "GLU": 9,
    "GLN": 9,
    "HIS": 10,
    "PHE": 11,
    "TYR": 12,
    "ARG": 11,
    "TRP": 14,
}


def flory_end_to_end(n_residues: int) -> float:
    """Analytical end-to-end distance estimate in Angstroms.

    The dependency-free fallback used when sparrow/ALBATROSS is not installed. See
    :data:`FLORY_RE_PREFACTOR` for provenance and limitations -- in particular this
    is blind to sequence composition and underestimates long expanded chains.

    Parameters
    ----------
    n_residues
        Number of residues in the disordered region. Must be positive.

    Returns
    -------
    float
        Estimated mean end-to-end distance in Angstroms.
    """
    if n_residues <= 0:
        raise ValueError(f"n_residues must be positive, got {n_residues}")
    return float(FLORY_RE_PREFACTOR * n_residues**FLORY_RE_EXPONENT)


def contour_length(n_residues: int) -> float:
    """Fully extended CA-trace length in Angstroms.

    The physical ceiling on any end-to-end distance: a target above this is
    unsatisfiable and callers should clamp to it rather than spin.
    """
    if n_residues <= 0:
        raise ValueError(f"n_residues must be positive, got {n_residues}")
    return float((n_residues - 1) * CA_CA_BOND_LENGTH)


def resolve_mode(mode: str) -> float:
    """Map a named build mode to its multiplier on predicted end-to-end distance.

    Raises
    ------
    ~dodo.exceptions.InvalidParameterError
        If ``mode`` is not a recognized mode name. The message lists the valid options, since
        this is the most common user-facing input error. That class is both a
        :class:`~dodo.exceptions.DodoError` and a :exc:`ValueError`, so either ``except`` works.
    """
    from .exceptions import InvalidParameterError

    try:
        return MODES[mode]
    except KeyError:
        valid = ", ".join(sorted(MODES))
        raise InvalidParameterError(f"Unknown mode {mode!r}. Valid modes are: {valid}") from None


def backbone_angle_grid() -> np.ndarray:
    """Integer-degree CA-CA-CA angles to sample, ordered by preference.

    Ordered by distance from :data:`BACKBONE_ANGLE_IDEAL` so that a caller taking
    the first non-clashing candidate automatically prefers realistic geometry. The
    pre-rewrite vectorized generator lost this ordering (it emitted 91, 92, 93...),
    silently discarding the bias toward physically likely angles.
    """
    angles = np.arange(BACKBONE_ANGLE_MIN, BACKBONE_ANGLE_MAX + 1.0)
    return angles[np.argsort(np.abs(angles - BACKBONE_ANGLE_IDEAL), kind="stable")]
