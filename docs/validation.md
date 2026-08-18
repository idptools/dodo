# Validation

DODO ships a validator for structural geometry, and points it at its own output.

```bash
dodo validate p53_dodo.pdb
```

```python
from dodo.validate import validate_structure

report = validate_structure(structure)
if not report.ok:
    print(report.describe())
```

It works on any PDB or mmCIF file, not just DODO's output.

## Where the reference numbers come from

Measured, not transcribed. Across **105,299,848 bonds in 23,587 AlphaFold structures** — the human
proteome — giving a mean and standard deviation for each of 173 distinct bonds across the 20
standard residues.

A few of the values that matter:

| Quantity | Measured |
|---|---|
| Peptide C–N bond | 1.3376 ± 0.0086 Å |
| CA–CA, trans peptide | 3.8317 ± 0.0842 Å |
| CA–CA, cis peptide | 3.1912 Å (2.44% of bonds) |
| Shortest real bond of any kind | 1.2026 Å at its 0.1st percentile (ASN CG–OD1) |

That last one sets the floor below which no bonding relationship, no strain and no crystallographic
oddity can put two heavy atoms.

## The invariant: no exclusion may hide an impossible distance

Every structural validator needs exclusions. One that flags covalently bonded atoms, or 1-3 pairs
held at 2.2 Å by a bond angle, is useless — it buries real findings under thousands of correct ones.

But an exclusion is a statement about *why* two atoms are close, and it quietly becomes a statement
that they may be *arbitrarily* close. All three of DODO's validators independently acquired the same
blind spot: a pair of atoms at exactly 0.000 Å was invisible to each of them whenever the pair fell
inside an exclusion.

- Bond validation skipped any geometric check on residues whose deposited numbering was not
  consecutive.
- Clash detection excluded pairs within three bonds unconditionally on distance, so 198 of 300
  sampled coincident pairs reported nothing.
- CONECT validation applied a single flat 0.9 Å floor, so two alpha carbons declared bonded at
  1.0 Å passed.

That is not hypothetical. Coincident atoms are the defect that actually shipped: DODO rebuilt a
region's alpha carbons and left the rest of each residue behind, and the result was visible in a
viewer as long spurious lines through the structure. Nothing in the codebase caught it.

So {func}`~dodo.validate.find_impossible_pairs` runs first, with **no exclusions at all** and no
`exclude` parameter to pass, and its findings are merged unconditionally. Keeping it in its own
module with nothing to configure is what stops the invariant being quietly re-broken.

## Generated geometry versus inherited geometry

DODO moves folded domains as rigid bodies and never regenerates their atoms, so a defect that
arrives in the input is still there afterwards — faithfully. Reporting that without saying so makes
DODO look responsible for its input.

Validating a rebuilt dnmt3a reports seven bond findings. Three of them are AlphaFold's own distorted
HIS613 imidazole ring, whose CE1–NE2 bond measures 2.547 Å against a 1.341 Å reference — 232
standard deviations out, and present at that value in the input file. The other four are the seams,
where a rebuilt region meets a folded domain.

So findings are attributed rather than merely counted:

```
INVALID: 7 bond (3 inherited from the input), 2 clash, 0 CONECT
```

Provenance is derivable exactly when validating a `Structure` in process, where the domains still
know which regions were rebuilt. Validating a *file*, that metadata is gone — but DODO leaves a
structural signature: it rebuilds regions with a backbone and never builds side chains, so a residue
carrying fewer atoms than its own residue type requires is one DODO built. A finding that touches
none of those residues is on geometry DODO did not build, and is reported as inherited. Glycine,
which is complete with a bare backbone, is not mistaken for a rebuilt residue.

The seams are deliberately *not* counted as inherited. They touch a residue DODO built, so they are
DODO's own compromise to report — not something to blame on the input.

This is why DODO's own tests assert a *differential* invariant — that DODO introduces no defect the
input did not already have — rather than that the output is unconditionally clean. Several real
AlphaFold inputs are not clean, and one in DODO's test corpus contains 92 impossible pairs of its
own.

## What is checked

**Bond lengths.** Every intra-residue bond against its measured reference, every peptide C–N, and
every CA–CA virtual bond in rebuilt regions. Findings carry a provenance and a deviation in
standard deviations.

**Steric clashes.** Element-aware van der Waals limits for all-atom contacts, and
`CA_CLASH_DISTANCE` for contacts involving a CA-only residue. Pairs within two residues in sequence
are exempt, being held close by the backbone rather than clashing.

**CONECT records.** That every declared bond is geometrically plausible, that every expected
backbone bond is declared, and that serial numbers resolve. This is a file-level check, so it needs
the written file rather than a `Structure`.

## Calibration

A validator that has never been pointed at a real structure is a guess, and the specific way it
fails — thousands of findings on perfectly good input — is invisible in synthetic tests, because
synthetic structures have no hydrogen bonds, no disulfides, no cis-prolines and no unresolved side
chains.

So it was calibrated against 300 held-out AlphaFold structures: **0.0666% of 1,555,229 bonds
flagged, with 288 of 300 structures completely clean**. Over a 117-structure corpus stratified by
region topology, DODO's own output contains **zero** impossible separations absent from the input
and **zero** bad bonds in geometry it generated.

Runtime is 0.22 s over 61,511 atoms.

## A note on AlphaFold's own geometry

It degrades where AlphaFold is unsure, which is worth knowing before validating an unmodified model
strictly. 2.44% of its CA–CA bonds fall below 3.3 Å, and those occur at a mean pLDDT of 38.8 against
72.2 for normal-length bonds. AlphaFold also ships genuinely broken residues: dnmt3a's HIS613, above,
is 232 standard deviations out on one bond.

None of that is a bug in DODO, and none of it is repaired by DODO.
