"""Steric clash detection for finished DODO structures.

What makes this hard is not finding close contacts -- a KD-tree does that in milliseconds
-- it is knowing which close contacts are *wrong*. A folded domain taken from a crystal
structure is full of atom pairs well inside any plausible clash cutoff: covalent bonds at
1.2-1.8 A, angle-fixed 1-3 pairs from 2.10 A, disulfides at 2.03 A, and every
secondary-structure hydrogen bond at 2.6-3.2 A. Measured on the fixtures in this
repository, a naive "flag everything closer than 3.2 A" check reports 210,335 pairs on the
6kn7 EM assembly and 361,957 across the four deposited and predicted fixtures, of which
146 are defects. A validator like that is worse than no validator, because the user learns
to ignore it.

So this module is mostly a careful account of what is *not* a clash, and every exclusion
below is justified by a measurement on the six fixture structures rather than by taste.

The two regimes DODO output actually contains
---------------------------------------------
DODO does not emit a uniform all-atom structure. Folded domains keep every atom of the
input and are moved only as rigid bodies; rebuilt regions -- IDRs, and loops inside folded
domains -- are ALPHA CARBON ONLY, because placing real backbone atoms there is the next
priority and is not done yet. One output file therefore contains all-atom stretches and
CA-only stretches, joined at a real chain connection.

Both regimes are handled explicitly:

* Between two all-atom residues, contacts are judged against element-aware minima derived
  from van der Waals radii (see :data:`VDW_RADII`) less a measured allowed overlap.
* Any pair involving the alpha carbon of a CA-only residue is judged against
  :data:`dodo.constants.CA_CLASH_DISTANCE` instead. A lone CA stands in for a whole
  residue's excluded volume, so a carbon van der Waals radius is the wrong yardstick for
  it, and 3.20 A is the value DODO's own builder enforces while placing those atoms.

Distinguishing a defect DODO introduced from one it preserved
-------------------------------------------------------------
It cannot be done from the output file alone, and this module does not pretend otherwise:
it reports what it measures. What it does instead is refuse to flag the things real
deposited structures legitimately contain, so that a finding on DODO output is worth
reading. Measured, over the unmodified fixtures:

===========================  =========  ================================================
fixture                      violations what they are
===========================  =========  ================================================
dnmt3a.pdb (AF2, 7146 atoms)         0   -
arf19.pdb (AF2, 8467 atoms)          0   -
p300.pdb (AF2, 18457 atoms)         15   13 from one genuine overlap in the input, where
                                         AlphaFold superimposed MET121 on PRO659 in a
                                         pLDDT-30 region; 2 short HIS-TRP contacts
6kn7.pdb (EM, 61511 atoms)         131   real side-chain overlaps in a 6.60 A model,
                                         where side chains are not resolved at all; 77
                                         distinct sequence positions repeated across the
                                         29 NCS copies
test.pdb (CA-only)              14,441   a degenerate pre-rewrite trace: CA-CA bonds
                                         1.23-562 A, 4,259 pseudo-angles below 75 deg
testing_translation.pdb (v1)     1,301   pre-rewrite output with the orphaned-atom bug
                                         fixed in commit 026c4f1: atoms at 0.000 A
===========================  =========  ================================================

Zero on well-formed input, five figures on known-bad input, and the two intermediate cases
are defects that are genuinely present in the file. Three of these fixtures also ship as
mmCIF, and the counts agree exactly there (0, 15 and 131), which says the findings are
properties of the coordinates rather than of the reader. That is the calibration; the tests
in ``tests/unit/test_validate_clashes.py`` assert it so it cannot silently drift.

Measured cost on the largest fixture: 0.185 s for 61,511 atoms, of which 0.073 s is
building the covalent bond graph, 0.063 s is expanding it to 1-3 and 1-4 pairs, and 0.019 s
is the KD-tree and its 183,225 candidate pairs.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final, Literal

import numpy as np
from scipy.spatial import cKDTree

from ..constants import CA_CLASH_DISTANCE, CLASH_EXCLUDE_WITHIN_RESIDUES
from ..exceptions import EmptyStructureError, GeometryError
from ..structure import Structure
from .reference import BOND_CUTOFF, RESIDUE_BONDS

__all__ = [
    "ALLOWED_OVERLAP",
    "DEFAULT_VDW_RADIUS",
    "DISULFIDE_MAX_DISTANCE",
    "EXCLUDE_WITHIN_BONDS",
    "HBOND_ALLOWED_OVERLAP",
    "HEAVY_BOND_CUTOFF",
    "LARGE_ATOM_BOND_CUTOFF",
    "POLAR_ELEMENTS",
    "SKIPPED_ELEMENTS",
    "VDW_RADII",
    "ClashKind",
    "ClashReport",
    "ClashViolation",
    "ca_only_residue_mask",
    "contact_limit",
    "covalent_bond_pairs",
    "validate_clashes",
]

#: Van der Waals radii in Angstroms, from Bondi (1964) *J. Phys. Chem.* 68:441, which is
#: the table crystallographic clash analysis has used ever since. Only the elements that
#: appear in protein structures are listed; anything else falls back to
#: :data:`DEFAULT_VDW_RADIUS`.
#:
#: These are radii, not limits. The limits this module applies are radius sums less a
#: measured allowed overlap -- see :func:`contact_limit`.
VDW_RADII: Final[dict[str, float]] = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "SE": 1.90,
    "P": 1.80,
    "F": 1.47,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "ZN": 1.39,
    "MG": 1.73,
    "CA": 2.31,
    "FE": 2.00,
    "MN": 2.00,
    "NA": 2.27,
    "K": 2.75,
}

#: Radius used for an element not in :data:`VDW_RADII`, in Angstroms. Carbon's, because it
#: is the commonest protein element and the largest of the four that make up 99.9% of the
#: atoms in the fixtures. Guessing large here is the safe direction only for detection; it
#: is the unsafe direction for false positives, which is why the value is carbon's rather
#: than the largest in the table.
DEFAULT_VDW_RADIUS: Final[float] = 1.70

#: Elements treated as hydrogen-bond donors or acceptors for the purposes of the relaxed
#: limit in :func:`contact_limit`.
#:
#: Sulfur is deliberately NOT included even though it can weakly hydrogen bond. Measured
#: over the four deposited fixtures, sulfur does not need the allowance: the closest
#: non-bonded O-S contact is 2.754 A (0.57 A of overlap, inside the ordinary allowance) and
#: the closest S-S is 3.338 A. The only sulfur contacts that fall below their limits are the
#: deepest overlaps measured anywhere in the fixtures -- a 2.303 A C-S and a 2.239 A N-S,
#: both in the 6.60 A EM model -- and those are exactly what should be reported.
POLAR_ELEMENTS: Final[frozenset[str]] = frozenset({"N", "O"})

#: Allowed van der Waals overlap for an ordinary non-bonded contact, in Angstroms.
#:
#: MEASURED, from the two fixtures that contain no defects. Over all non-bonded,
#: non-neighbour contacts in the AlphaFold models dnmt3a and arf19, the deepest overlap
#: involving at least one non-polar atom is 0.448 A (dnmt3a) and 0.252 A (arf19). 0.60
#: clears the larger of those by 0.15 A, which is the margin for structures outside the
#: fixture set. The resulting limits are C-C 2.80, C-N 2.65, C-O 2.62, C-S 2.90,
#: S-S 3.00 A.
#:
#: It is deliberately well above what a crystallographer would use. MolProbity's clashscore
#: calls a contact bad at 0.40 A of overlap; applied here that reports 1,279 contacts across
#: the four deposited and predicted fixtures, 21 of them in dnmt3a and 16 in arf19, which
#: are ordinary tight packing in an AlphaFold model rather than errors. DODO's job is to
#: catch geometry DODO broke, and measured on fresh ``rebuild()`` output those defects are
#: 0.1-2.5 A contacts, not 2.7 A ones. Being strict here would cost far more in false
#: positives than it buys in sensitivity.
ALLOWED_OVERLAP: Final[float] = 0.60

#: Allowed van der Waals overlap when BOTH atoms are in :data:`POLAR_ELEMENTS`, in
#: Angstroms.
#:
#: MEASURED. A hydrogen bond puts two heavy atoms 2.6-3.2 A apart, and every alpha helix and
#: beta sheet in the input is held together by them: across the four deposited fixtures
#: there are 12,063 non-bonded N/O-to-N/O pairs between 2.6 and 3.2 A (659 in dnmt3a, 476 in
#: arf19, 1,194 in p300, 9,031 in 6kn7). Under a single 3.2 A cutoff every one of them is a
#: false finding, which is the mistake this module exists to avoid.
#:
#: The element-aware limits above already clear the bulk of them -- N-O at 3.07 - 0.60 =
#: 2.47 A sits below the normal hydrogen-bonding range -- so what this constant buys is the
#: short tail. The deepest polar overlap measured is 0.832 A (an N-O pair at 2.238 A in
#: 6kn7) and the deepest in a high-quality model is 0.807 A (ASN611 ND2 to HIS613 NE2 at
#: 2.293 A in dnmt3a, radius sum 3.10 -- a short but perfectly real hydrogen bond). 0.90
#: clears both, giving limits N-N 2.20, N-O 2.17, O-O 2.14 A.
#:
#: Setting this equal to :data:`ALLOWED_OVERLAP` -- that is, treating polar pairs like any
#: other -- adds 99 findings across the four fixtures (1 in dnmt3a, 0 in arf19, 1 in p300,
#: 98 in 6kn7), every one of them a real hydrogen bond.
HBOND_ALLOWED_OVERLAP: Final[float] = 0.90

#: Longest SG-SG separation treated as a disulfide bond rather than a clash, in Angstroms.
#:
#: MEASURED: the fixtures contain 3 disulfides, at 2.029, 2.034 and 2.044 A, and the
#: next-closest SG-SG pair anywhere in them is 3.338 A (two cysteines in a metal site).
#: 2.50 A sits in that gap with room for a strained bridge, and no pair of free thiols can
#: sit that close, so nothing is masked by choosing the upper end. Without this exemption
#: each disulfide is the deepest "clash" in its structure, at 1.0 A inside the S-S limit.
DISULFIDE_MAX_DISTANCE: Final[float] = 2.50

#: Distance below which two heavy atoms of the same residue are taken to be bonded when
#: the residue is not in :data:`~dodo.validate.reference.RESIDUE_BONDS`, in Angstroms.
#: Reuses the measured cutoff from that module, which separates covalent bonds (1.00-1.85)
#: from 1-3 contacts (2.10 and up) with only 4 exceptions in 197,906 measured pairs.
HEAVY_BOND_CUTOFF: Final[float] = BOND_CUTOFF

#: Same, for pairs where either atom is sulfur, selenium or phosphorus, in Angstroms.
#: Those bonds are genuinely longer than :data:`HEAVY_BOND_CUTOFF` allows -- C-SE in
#: selenomethionine is 1.95 A and S-S is 2.03 A -- so a single light-atom cutoff would
#: leave the selenium of an MSE residue unbonded and then report it as clashing with the
#: carbon it is bonded to.
LARGE_ATOM_BOND_CUTOFF: Final[float] = 2.10

#: Elements skipped entirely, because DODO never places them and heavy-atom contact
#: analysis is the standard measure. Keeping them would require a hydrogen bonding model:
#: an explicit-H input names its hydrogens H, HA, HB2 and so on, none of which appear in
#: the measured bond tables, so every hydrogen would be an unbonded atom 1.0 A from its
#: parent and would be reported as a clash.
SKIPPED_ELEMENTS: Final[frozenset[str]] = frozenset({"H", "D"})

#: Number of bonds within which a pair is exempt from contact checking.
#:
#: 3, i.e. 1-2 (bonded), 1-3 and 1-4 pairs are all excluded. Measured on the four deposited
#: fixtures, 1-2 pairs run 1.010-3.246 A, 1-3 pairs start at 2.096 A and 1-4 pairs start at
#: 2.404 A. All three ranges reach below the ordinary limit of 2.80 A, so all three have to
#: go, and the counts say by how much: without these exclusions the four fixtures report
#: 97,457 bonded pairs, 115,389 1-3 pairs and 2,016 1-4 pairs as clashes.
#:
#: They are fixed by bond lengths, bond angles and one rotatable dihedral respectively. A
#: 1-4 pair reaching 2.40 A simply means that dihedral is near cis, which is legitimate
#: geometry and not a collision.
#:
#: The cost of excluding 1-4 pairs is real but small: it means this module cannot see a
#: defect confined to a single dihedral. That is bond-and-angle validation's job, not
#: clash detection's, and DODO does not alter intra-residue geometry -- folded domains move
#: as rigid bodies.
EXCLUDE_WITHIN_BONDS: Final[int] = 3

#: Kinds of violation :func:`validate_clashes` can report.
ClashKind = Literal["atom_contact", "ca_contact", "non_finite"]


def contact_limit(
    element_a: str,
    element_b: str,
    *,
    allowed_overlap: float = ALLOWED_OVERLAP,
    hbond_allowed_overlap: float = HBOND_ALLOWED_OVERLAP,
) -> float:
    """Minimum acceptable separation for a non-bonded contact between two elements.

    ``r(a) + r(b) - overlap``, where the radii are Bondi's (:data:`VDW_RADII`) and the
    overlap is :data:`HBOND_ALLOWED_OVERLAP` when both elements are in
    :data:`POLAR_ELEMENTS` and :data:`ALLOWED_OVERLAP` otherwise.

    Parameters
    ----------
    element_a, element_b
        Element symbols. Case-insensitive; an unrecognized element gets
        :data:`DEFAULT_VDW_RADIUS`.
    allowed_overlap, hbond_allowed_overlap
        Overrides for the two measured overlap allowances, in Angstroms.

    Returns
    -------
    float
        The limit in Angstroms. Contacts strictly below it are clashes.

    Examples
    --------
    >>> round(contact_limit("C", "C"), 3)
    2.8
    >>> round(contact_limit("N", "O"), 3)
    2.17
    """
    a = element_a.upper()
    b = element_b.upper()
    radius_sum = VDW_RADII.get(a, DEFAULT_VDW_RADIUS) + VDW_RADII.get(b, DEFAULT_VDW_RADIUS)
    overlap = (
        hbond_allowed_overlap if a in POLAR_ELEMENTS and b in POLAR_ELEMENTS else allowed_overlap
    )
    return float(radius_sum - overlap)


def ca_only_residue_mask(structure: Structure) -> np.ndarray:
    """Which residues are represented by an alpha carbon and nothing else.

    These are DODO's rebuilt regions. The distinction is drawn from the atoms actually
    present rather than from the region assignment, so it holds for a structure read back
    from a file with no domain annotation -- which is how a user validates DODO output.

    Parameters
    ----------
    structure
        The structure to inspect.

    Returns
    -------
    np.ndarray
        Boolean mask over residues, ``(n_residues,)``. ``True`` where the residue has
        exactly one atom and that atom is named ``CA``.

    Notes
    -----
    A residue with N, CA, C and O but no side chain (what ``rebuild(all_atom=True)``
    produces) is *not* CA-only: it has real backbone atoms, so element-aware contact
    limits apply to it and are more informative than a single CA exclusion sphere.
    """
    counts = np.diff(structure.residue_atom_offsets)
    single = counts == 1
    mask = np.zeros(structure.n_residues, dtype=bool)
    if not single.any():
        return mask
    first_atom = structure.residue_atom_offsets[:-1][single]
    mask[single] = structure.atom_name[first_atom] == "CA"
    return mask


def _intra_residue_bonds(structure: Structure) -> list[tuple[int, int]]:
    """Covalent bonds inside residues, from the measured tables with a distance fallback."""
    bonds: list[tuple[int, int]] = []
    names = structure.atom_name
    offsets = structure.residue_atom_offsets
    elements = structure.element
    for residue in range(structure.n_residues):
        start = int(offsets[residue])
        stop = int(offsets[residue + 1])
        if stop - start < 2:
            continue
        table = RESIDUE_BONDS.get(str(structure.residue_name[residue]))
        local: dict[str, int] = {}
        for atom in range(start, stop):
            # First occurrence wins. A duplicated atom name inside one residue means the
            # reader kept two alternate conformers, and bonding both copies would invent a
            # bond between them; the duplicate is left unbonded and picked up by the
            # fallback below, which bonds it to its real neighbours.
            local.setdefault(str(names[atom]), atom)
        bonded = np.zeros(stop - start, dtype=bool)
        if table is not None:
            for atom_a, atom_b in table:
                index_a = local.get(atom_a)
                index_b = local.get(atom_b)
                if index_a is not None and index_b is not None:
                    bonds.append((min(index_a, index_b), max(index_a, index_b)))
                    bonded[index_a - start] = True
                    bonded[index_b - start] = True
        # Atoms the table said nothing about: the selenium of a selenomethionine deposited
        # as MET, the phosphate of a phosphoserine deposited as SER, or every atom of a
        # residue with no table at all. Bond them by distance. Without this, each one is an
        # unbonded atom sitting 1.2-1.9 A from its own parent and is reported as clashing
        # with it. The terminal OXT used to be the common case here -- 32 of them across the
        # fixtures, one per chain -- and is now covered by the reference table directly; the
        # fallback stays because modified residues are not, and cannot all be enumerated.
        unbonded = np.flatnonzero(~bonded) + start
        if unbonded.size:
            coords = structure.xyz[start:stop]
            distances = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
            large = np.array(
                [str(e).upper() in {"S", "SE", "P"} for e in elements[start:stop]], dtype=bool
            )
            cutoff = np.where(
                large[:, None] | large[None, :], LARGE_ATOM_BOND_CUTOFF, HEAVY_BOND_CUTOFF
            )
            candidate = (distances < cutoff) & (distances > 0.0)
            for orphan in unbonded:
                for other in np.flatnonzero(candidate[orphan - start]) + start:
                    bonds.append((int(min(orphan, other)), int(max(orphan, other))))
    return bonds


def _inter_residue_bonds(structure: Structure) -> list[tuple[int, int]]:
    """Peptide C(i)-N(i+1) bonds and disulfide SG-SG bridges."""
    bonds: list[tuple[int, int]] = []
    names = structure.atom_name
    offsets = structure.residue_atom_offsets

    def atom_of(residue: int, name: str) -> int | None:
        start = int(offsets[residue])
        stop = int(offsets[residue + 1])
        hits = np.flatnonzero(names[start:stop] == name)
        return start + int(hits[0]) if hits.size else None

    chain_index = structure.chain_index
    for residue in range(structure.n_residues - 1):
        # Only within a chain: two residues adjacent in index but on different chains are
        # not bonded, and a real collision between two chain termini must stay visible.
        if chain_index[residue] != chain_index[residue + 1]:
            continue
        carbon = atom_of(residue, "C")
        nitrogen = atom_of(residue + 1, "N")
        if carbon is None or nitrogen is None:
            continue
        # Distance-gated so that a chain break -- residues adjacent in index because the
        # unmodelled ones in between are simply absent -- does not get a phantom bond
        # spanning tens of Angstroms, which would then exempt everything within 3 bonds of
        # it at the far end of the gap.
        if float(np.linalg.norm(structure.xyz[carbon] - structure.xyz[nitrogen])) <= BOND_CUTOFF:
            bonds.append((min(carbon, nitrogen), max(carbon, nitrogen)))

    # Disulfides. A CYS SG-SG pair at ~2.05 A is a bond; there are 20 in the fixtures and
    # every one would otherwise be the worst "clash" in its structure.
    sulfurs = np.flatnonzero(names == "SG")
    if sulfurs.size > 1:
        tree = cKDTree(structure.xyz[sulfurs])
        for local_a, local_b in tree.query_pairs(DISULFIDE_MAX_DISTANCE):
            atom_a = int(sulfurs[local_a])
            atom_b = int(sulfurs[local_b])
            bonds.append((min(atom_a, atom_b), max(atom_a, atom_b)))
    return bonds


def covalent_bond_pairs(structure: Structure) -> np.ndarray:
    """Every covalent bond in a structure, as atom index pairs.

    Three sources, in order of authority: the measured intra-residue tables in
    :mod:`dodo.validate.reference`, a distance fallback for atoms and residues those
    tables do not cover, and the two inter-residue bonds that exist in a protein --
    peptide C(i)-N(i+1) and disulfide SG-SG.

    Parameters
    ----------
    structure
        The structure to bond.

    Returns
    -------
    np.ndarray
        ``(n_bonds, 2)`` int64 array of atom index pairs, each with the lower index first,
        sorted and deduplicated. Empty with shape ``(0, 2)`` for a CA-only structure,
        which has no covalent bonds at all -- a CA-CA virtual bond is not one.

    Examples
    --------
    >>> from dodo.io import read_structure
    >>> bonds = covalent_bond_pairs(read_structure("x.pdb"))   # doctest: +SKIP
    """
    pairs = _intra_residue_bonds(structure) + _inter_residue_bonds(structure)
    if not pairs:
        return np.zeros((0, 2), dtype=np.int64)
    array = np.asarray(pairs, dtype=np.int64)
    array = np.unique(array, axis=0)
    return array


def _pair_keys(atom_a: np.ndarray, atom_b: np.ndarray, n_atoms: int) -> np.ndarray:
    """Encode unordered atom index pairs as sortable int64 keys."""
    low = np.minimum(atom_a, atom_b).astype(np.int64)
    high = np.maximum(atom_a, atom_b).astype(np.int64)
    keys: np.ndarray = low * n_atoms + high
    return keys


def _neighbour_table(bonds: np.ndarray, n_atoms: int) -> tuple[np.ndarray, np.ndarray]:
    """Build a padded adjacency table ``(n_atoms, max_degree)`` and a per-atom degree."""
    degree = np.zeros(n_atoms, dtype=np.int64)
    np.add.at(degree, bonds[:, 0], 1)
    np.add.at(degree, bonds[:, 1], 1)
    max_degree = int(degree.max()) if n_atoms else 0
    table = np.full((n_atoms, max(max_degree, 1)), -1, dtype=np.int64)
    cursor = np.zeros(n_atoms, dtype=np.int64)
    for atom_a, atom_b in bonds:
        table[atom_a, cursor[atom_a]] = atom_b
        cursor[atom_a] += 1
        table[atom_b, cursor[atom_b]] = atom_a
        cursor[atom_b] += 1
    return table, degree


def _bonded_within_keys(bonds: np.ndarray, n_atoms: int) -> np.ndarray:
    """Sorted keys of every pair separated by at most three bonds.

    Built from the adjacency table directly rather than by multiplying a sparse matrix
    three times: the graph has maximum degree 4, so the three relations are cheap to write
    out. 1-2 pairs are the bonds; 1-3 pairs are two distinct neighbours of a common atom;
    1-4 pairs are a neighbour of one bond end crossed with a neighbour of the other.

    "Three" is :data:`EXCLUDE_WITHIN_BONDS`, which is where the justification lives. It is
    written out rather than looped over because each of the three relations has its own
    closed form; a general breadth-first version would be slower and no clearer, and the
    depth is not something a caller should be tuning.
    """
    if bonds.shape[0] == 0:
        return np.zeros(0, dtype=np.int64)
    assert EXCLUDE_WITHIN_BONDS == 3, "the three relations below are written out by hand"
    table, _ = _neighbour_table(bonds, n_atoms)
    width = table.shape[1]

    keys = [_pair_keys(bonds[:, 0], bonds[:, 1], n_atoms)]

    # 1-3: every pair of distinct neighbours of the same atom.
    for i in range(width):
        for j in range(i + 1, width):
            left = table[:, i]
            right = table[:, j]
            valid = (left >= 0) & (right >= 0)
            if valid.any():
                keys.append(_pair_keys(left[valid], right[valid], n_atoms))

    # 1-4: neighbours of one bond end against neighbours of the other, skipping the two
    # ends themselves (those pairs are 1-2 or 1-3 and are already covered).
    ends_a = bonds[:, 0]
    ends_b = bonds[:, 1]
    for i in range(width):
        left = table[ends_a, i]
        for j in range(width):
            right = table[ends_b, j]
            valid = (
                (left >= 0) & (right >= 0) & (left != right) & (left != ends_b) & (right != ends_a)
            )
            if valid.any():
                keys.append(_pair_keys(left[valid], right[valid], n_atoms))

    merged = np.unique(np.concatenate(keys))
    return merged


@dataclass(frozen=True, slots=True)
class ClashViolation:
    """One pair of atoms placed closer together than physics allows.

    Attributes
    ----------
    kind
        ``"atom_contact"`` for a contact judged against element-aware van der Waals
        limits, ``"ca_contact"`` for one involving the alpha carbon of a CA-only residue
        and judged against :data:`dodo.constants.CA_CLASH_DISTANCE`, and ``"non_finite"``
        for an atom whose coordinate is NaN or infinite -- which cannot be tested for
        clashes at all, so it is reported rather than silently dropped.
    atom_a, atom_b
        Atom indices, ``atom_a < atom_b``. ``atom_b`` is ``None`` for ``"non_finite"``,
        which concerns one atom.
    residue_a, residue_b
        Residue indices the atoms belong to, 0-based positional. ``residue_b`` is ``None``
        for ``"non_finite"``.
    label_a, label_b
        The residues as :meth:`dodo.structure.Structure.residue_label` renders them, e.g.
        ``"A:GLU142"``. ``label_b`` is ``None`` for ``"non_finite"``.
    atom_name_a, atom_name_b
        PDB atom names, e.g. ``"CA"``, ``"OD1"``.
    separation
        Measured distance in Angstroms; NaN for ``"non_finite"``.
    limit
        The minimum acceptable separation this pair was judged against, in Angstroms; NaN
        for ``"non_finite"``.
    message
        A complete, self-contained sentence naming both residues, the measured distance and
        the limit.
    """

    kind: ClashKind
    atom_a: int
    atom_b: int | None
    residue_a: int
    residue_b: int | None
    label_a: str
    label_b: str | None
    atom_name_a: str
    atom_name_b: str | None
    separation: float
    limit: float
    message: str

    def __str__(self) -> str:
        return self.message

    @property
    def overlap(self) -> float:
        """How far inside the limit this contact sits, in Angstroms. NaN if unmeasurable."""
        return self.limit - self.separation


@dataclass(frozen=True, slots=True, eq=False)
class ClashReport:
    """The result of :func:`validate_clashes`: what was checked, and what failed.

    ``eq=False`` for the same reason :class:`dodo.geometry.metrics.TraceReport` uses it:
    a generated ``__eq__`` would compare numpy arrays with ``==`` and raise "truth value of
    an array is ambiguous" the first time two reports were compared or put in a set.

    Attributes
    ----------
    n_atoms
        Atoms in the structure.
    n_residues
        Residues in the structure.
    n_ca_only_residues
        How many residues are alpha-carbon only, i.e. how much of the structure is a
        rebuilt region. 0 for an unmodified all-atom input.
    n_pairs_checked
        Atom pairs whose separation was measured against a limit -- that is, pairs inside
        the search radius that survived every exclusion. Not the total number of pairs in
        the structure, which is quadratic and never enumerated.
    n_pairs_excluded
        Pairs inside the search radius that were excluded as bonded, 1-3, 1-4, sequence
        neighbours or disulfide-bridged. Reported because it is the number that shows
        whether the exclusions are doing anything.
    n_atoms_skipped
        Atoms left out of the geometry entirely: hydrogens and non-finite coordinates.
    violations
        Every violation found, ordered by residue then atom index.
    atom_limits
        The ``(low, high)`` range of element-aware limits actually applied, in Angstroms,
        or ``None`` when no all-atom pair was checked.
    ca_clash_distance
        The limit applied to pairs involving a CA-only residue, in Angstroms.
    search_radius
        The KD-tree query radius used, in Angstroms: the largest limit in play, since no
        pair beyond it can violate anything.
    """

    n_atoms: int
    n_residues: int
    n_ca_only_residues: int
    n_pairs_checked: int
    n_pairs_excluded: int
    n_atoms_skipped: int
    violations: tuple[ClashViolation, ...]
    atom_limits: tuple[float, float] | None
    ca_clash_distance: float
    search_radius: float

    @property
    def ok(self) -> bool:
        """True when nothing was violated."""
        return not self.violations

    def __bool__(self) -> bool:
        return self.ok

    def __repr__(self) -> str:
        return (
            f"ClashReport({self.n_atoms} atoms, "
            f"{'no clashes' if self.ok else f'{len(self.violations)} clashes'})"
        )

    def of_kind(self, kind: ClashKind) -> tuple[ClashViolation, ...]:
        """All violations of one kind, in report order."""
        return tuple(v for v in self.violations if v.kind == kind)

    @property
    def worst(self) -> ClashViolation | None:
        """The violation with the deepest overlap, or ``None`` if there are none.

        Deepest overlap rather than smallest separation, because the limits differ between
        pairs: a 2.30 A N-O contact (0.13 A inside its limit) is a much milder defect than
        a 2.30 A CA-CA contact (0.90 A inside its limit), and ranking by raw distance would
        put them the wrong way round.
        """
        measured = [v for v in self.violations if math.isfinite(v.overlap)]
        if not measured:
            return self.violations[0] if self.violations else None
        return max(measured, key=lambda v: v.overlap)

    @property
    def violating_residues(self) -> tuple[int, ...]:
        """Residue indices involved in at least one violation, ascending and unique."""
        found: set[int] = set()
        for violation in self.violations:
            found.add(violation.residue_a)
            if violation.residue_b is not None:
                found.add(violation.residue_b)
        return tuple(sorted(found))

    def summary(self) -> str:
        """One-line summary of what was checked and what was found."""
        parts = [
            f"{self.n_atoms} atoms",
            f"{self.n_residues} residues ({self.n_ca_only_residues} CA-only)",
            f"{self.n_pairs_checked} pairs checked, {self.n_pairs_excluded} excluded",
        ]
        if self.ok:
            parts.append("no clashes")
        else:
            parts.append(f"{len(self.violations)} clashes")
            worst = self.worst
            if worst is not None and math.isfinite(worst.overlap):
                parts.append(f"worst {worst.separation:.3f} A vs {worst.limit:.2f} A limit")
        return ", ".join(parts)

    def describe(self, max_violations: int = 10) -> str:
        """Full multi-line description: the summary, then the violations one per line.

        Parameters
        ----------
        max_violations
            Truncate the list after this many, with a count of the remainder. Two
            interpenetrating domains produce thousands of pairs and an unreadable message;
            the worst few identify the failure. Violations are listed worst-overlap first
            for exactly that reason.
        """
        lines = [self.summary()]
        ordered = sorted(
            self.violations,
            key=lambda v: -v.overlap if math.isfinite(v.overlap) else -math.inf,
        )
        for violation in ordered[:max_violations]:
            lines.append(f"  - {violation.message}")
        remaining = len(self.violations) - max_violations
        if remaining > 0:
            lines.append(f"  - ... and {remaining} more clash(es)")
        return "\n".join(lines)

    def raise_if_invalid(self) -> None:
        """Raise :class:`~dodo.exceptions.GeometryError` if anything clashed.

        The one-liner for a caller that must not hand on a structure containing
        interpenetrating atoms.
        """
        if not self.ok:
            raise GeometryError(f"Steric clashes detected: {self.describe()}")


def _check_float(name: str, value: float) -> float:
    """Coerce to float and reject non-finite or negative values."""
    number = float(value)
    # A NaN limit makes every comparison below False, which disables the check silently
    # instead of failing loudly -- the exact failure mode this package exists to avoid.
    if not math.isfinite(number):
        raise GeometryError(f"{name} must be finite, got {number}.")
    if number < 0.0:
        raise GeometryError(f"{name} must be non-negative, got {number}.")
    return number


def validate_clashes(
    structure: Structure,
    *,
    allowed_overlap: float = ALLOWED_OVERLAP,
    hbond_allowed_overlap: float = HBOND_ALLOWED_OVERLAP,
    ca_clash_distance: float = CA_CLASH_DISTANCE,
    ca_exclude_within: int = CLASH_EXCLUDE_WITHIN_RESIDUES,
    atom_exclude_within: int = 1,
) -> ClashReport:
    """Find atoms placed unphysically close together.

    Every close contact that a real structure legitimately contains is excluded first:
    covalent bonds, 1-3 and 1-4 pairs, disulfide bridges, sequence neighbours, and -- via
    a relaxed limit for N/O-to-N/O pairs -- hydrogen bonds. See the module docstring for
    the measurements behind each.

    Parameters
    ----------
    structure
        The structure to check. Both all-atom and CA-only regions are handled; see the
        module docstring.
    allowed_overlap
        Van der Waals overlap tolerated in an ordinary contact, in Angstroms. Defaults to
        the measured :data:`ALLOWED_OVERLAP`. Lower it to be stricter -- 0.40 reproduces
        MolProbity's clash criterion, which on these fixtures reports ordinary tight
        packing as well as real defects.
    hbond_allowed_overlap
        Overlap tolerated when both atoms are in :data:`POLAR_ELEMENTS`, in Angstroms.
        Defaults to the measured :data:`HBOND_ALLOWED_OVERLAP`. Setting it equal to
        ``allowed_overlap`` turns the hydrogen-bond allowance off, which is only useful for
        demonstrating how many findings it suppresses.
    ca_clash_distance
        Minimum separation for any pair involving a CA-only residue's alpha carbon, in
        Angstroms. Defaults to :data:`dodo.constants.CA_CLASH_DISTANCE`, the value DODO's
        builder enforces while placing those atoms.
    ca_exclude_within
        Sequence separation exempt from the CA rule, within a chain. Defaults to
        :data:`dodo.constants.CLASH_EXCLUDE_WITHIN_RESIDUES`. This exclusion is what keeps
        the boundary between a rebuilt region and a folded domain from being reported: a
        rebuilt CA(i) sits about 2.4 A from N(i+1) of the folded residue it connects to,
        which is correct backbone geometry and not a collision.
    atom_exclude_within
        Sequence separation exempt from the all-atom rule, within a chain. Defaults to 1,
        because consecutive residues are held close by the backbone whether or not the
        pair is within 3 bonds: measured on the deposited fixtures the closest non-bonded
        i/i+1 contact beyond 1-4 is 2.484 A (6kn7) and 2.569 A (p300), and 12 such pairs
        in 6kn7 fall below their element limits. Raising it to 2 would suppress real
        information -- residues two apart are where genuine packing errors start showing
        up, and the closest i/i+2 contact in the fixtures is a 2.293 A hydrogen bond that
        the polar allowance already clears on its own.

    Returns
    -------
    ClashReport
        Everything measured and every violation. Truthy when there are no clashes.

    Raises
    ------
    EmptyStructureError
        If the structure has no atoms or no residues. There is nothing to report about.
    GeometryError
        If a threshold argument is negative or non-finite, or an ``exclude_within`` is
        negative. Bad *geometry* is reported rather than raised -- that is the point of the
        report -- but an unusable argument is a caller bug.

    Examples
    --------
    >>> import dodo
    >>> from dodo.validate.clashes import validate_clashes
    >>> report = validate_clashes(dodo.read_structure("model.pdb"))   # doctest: +SKIP
    >>> print(report.describe())                                     # doctest: +SKIP
    """
    if structure.n_atoms == 0:
        raise EmptyStructureError(f"Structure has no atoms (source: {structure.source!r}).")
    if structure.n_residues == 0:
        raise EmptyStructureError(f"Structure has no residues (source: {structure.source!r}).")
    allowed_overlap = _check_float("allowed_overlap", allowed_overlap)
    hbond_allowed_overlap = _check_float("hbond_allowed_overlap", hbond_allowed_overlap)
    ca_clash_distance = _check_float("ca_clash_distance", ca_clash_distance)
    if ca_exclude_within < 0:
        raise GeometryError(f"ca_exclude_within must be non-negative, got {ca_exclude_within}.")
    if atom_exclude_within < 0:
        raise GeometryError(f"atom_exclude_within must be non-negative, got {atom_exclude_within}.")

    n_atoms = structure.n_atoms
    elements = np.char.upper(structure.element.astype("<U2"))
    violations: list[ClashViolation] = []

    # --- decide which atoms take part -------------------------------------------------
    is_hydrogen = np.isin(elements, list(SKIPPED_ELEMENTS))
    finite = np.all(np.isfinite(structure.xyz), axis=1)
    for atom in np.flatnonzero(~finite):
        residue = int(structure.residue_index[atom])
        label = structure.residue_label(residue)
        violations.append(
            ClashViolation(
                kind="non_finite",
                atom_a=int(atom),
                atom_b=None,
                residue_a=residue,
                residue_b=None,
                label_a=label,
                label_b=None,
                atom_name_a=str(structure.atom_name[atom]),
                atom_name_b=None,
                separation=float("nan"),
                limit=float("nan"),
                message=(
                    f"{label} atom {structure.atom_name[atom]} has a non-finite coordinate "
                    f"({structure.xyz[atom].tolist()}), so it cannot be checked for clashes; "
                    f"a NaN coordinate means an upstream sampler returned its failure "
                    f"instead of raising."
                ),
            )
        )
    included = finite & ~is_hydrogen
    n_skipped = int((~included).sum())

    # --- limits, and therefore the search radius --------------------------------------
    radii = np.array(
        [VDW_RADII.get(str(e), DEFAULT_VDW_RADIUS) for e in elements], dtype=np.float64
    )
    is_polar = np.isin(elements, list(POLAR_ELEMENTS))
    ca_only = ca_only_residue_mask(structure)
    atom_is_ca_only = ca_only[structure.residue_index]

    present = radii[included] if included.any() else np.array([DEFAULT_VDW_RADIUS])
    largest_atom_limit = float(2.0 * present.max() - min(allowed_overlap, hbond_allowed_overlap))
    search_radius = max(largest_atom_limit, ca_clash_distance if ca_only.any() else 0.0)

    if not included.any() or search_radius <= 0.0:
        # Nothing measurable. Still a real report, not an exception: a structure of one
        # atom has no pairs and that is not an error.
        return ClashReport(
            n_atoms=n_atoms,
            n_residues=structure.n_residues,
            n_ca_only_residues=int(ca_only.sum()),
            n_pairs_checked=0,
            n_pairs_excluded=0,
            n_atoms_skipped=n_skipped,
            violations=tuple(violations),
            atom_limits=None,
            ca_clash_distance=ca_clash_distance,
            search_radius=float(search_radius),
        )

    # --- candidate pairs --------------------------------------------------------------
    atom_indices = np.flatnonzero(included)
    tree = cKDTree(structure.xyz[atom_indices])
    local_pairs = tree.query_pairs(search_radius, output_type="ndarray")
    if local_pairs.shape[0] == 0:
        pair_a = np.zeros(0, dtype=np.int64)
        pair_b = np.zeros(0, dtype=np.int64)
    else:
        pair_a = atom_indices[local_pairs[:, 0]]
        pair_b = atom_indices[local_pairs[:, 1]]
    low = np.minimum(pair_a, pair_b)
    high = np.maximum(pair_a, pair_b)
    pair_a, pair_b = low, high

    residue_a = structure.residue_index[pair_a]
    residue_b = structure.residue_index[pair_b]
    separation = np.linalg.norm(structure.xyz[pair_a] - structure.xyz[pair_b], axis=1)

    # --- classify: CA rule or all-atom rule -------------------------------------------
    ca_pair = atom_is_ca_only[pair_a] | atom_is_ca_only[pair_b]
    limits = np.where(
        is_polar[pair_a] & is_polar[pair_b],
        radii[pair_a] + radii[pair_b] - hbond_allowed_overlap,
        radii[pair_a] + radii[pair_b] - allowed_overlap,
    )
    limits = np.where(ca_pair, ca_clash_distance, limits)

    # --- exclusions -------------------------------------------------------------------
    same_chain = structure.chain_index[residue_a] == structure.chain_index[residue_b]
    sequence_gap = np.abs(residue_a - residue_b)
    exclude_gap = np.where(ca_pair, ca_exclude_within, atom_exclude_within)
    excluded = same_chain & (sequence_gap <= exclude_gap) & (residue_a != residue_b)

    bonds = covalent_bond_pairs(structure)
    bonded_keys = _bonded_within_keys(bonds, n_atoms)
    if bonded_keys.size and pair_a.size:
        keys = _pair_keys(pair_a, pair_b, n_atoms)
        position = np.clip(np.searchsorted(bonded_keys, keys), 0, bonded_keys.size - 1)
        excluded |= bonded_keys[position] == keys

    checked = ~excluded
    clashing = np.flatnonzero(checked & (separation < limits))

    # --- report -----------------------------------------------------------------------
    for raw_index in clashing:
        index = int(raw_index)
        label_a = structure.residue_label(int(residue_a[index]))
        label_b = structure.residue_label(int(residue_b[index]))
        name_a = str(structure.atom_name[pair_a[index]])
        name_b = str(structure.atom_name[pair_b[index]])
        if ca_pair[index]:
            which = (
                f"{label_a} is"
                if atom_is_ca_only[pair_a[index]] and not atom_is_ca_only[pair_b[index]]
                else f"{label_b} is"
                if atom_is_ca_only[pair_b[index]] and not atom_is_ca_only[pair_a[index]]
                else f"{label_a} and {label_b} are"
            )
            message = (
                f"{label_a} {name_a} and {label_b} {name_b} are "
                f"{separation[index]:.3f} A apart, closer than the minimum CA-CA approach "
                f"of {limits[index]:.2f} A; {which} alpha-carbon only, so one atom stands "
                f"for the whole residue's excluded volume."
            )
            kind: ClashKind = "ca_contact"
        else:
            element_a = str(elements[pair_a[index]])
            element_b = str(elements[pair_b[index]])
            message = (
                f"{label_a} {name_a} and {label_b} {name_b} are "
                f"{separation[index]:.3f} A apart, closer than the {limits[index]:.2f} A "
                f"minimum for a non-bonded {element_a}-{element_b} contact (van der Waals "
                f"radii {radii[pair_a[index]]:.2f} + {radii[pair_b[index]]:.2f} A less the "
                f"allowed overlap)."
            )
            kind = "atom_contact"
        violations.append(
            ClashViolation(
                kind=kind,
                atom_a=int(pair_a[index]),
                atom_b=int(pair_b[index]),
                residue_a=int(residue_a[index]),
                residue_b=int(residue_b[index]),
                label_a=label_a,
                label_b=label_b,
                atom_name_a=name_a,
                atom_name_b=name_b,
                separation=float(separation[index]),
                limit=float(limits[index]),
                message=message,
            )
        )

    violations.sort(key=lambda v: (v.residue_a, v.residue_b or -1, v.atom_a, v.atom_b or -1))
    atom_checked = checked & ~ca_pair
    atom_limits = (
        (float(limits[atom_checked].min()), float(limits[atom_checked].max()))
        if atom_checked.any()
        else None
    )
    return ClashReport(
        n_atoms=n_atoms,
        n_residues=structure.n_residues,
        n_ca_only_residues=int(ca_only.sum()),
        n_pairs_checked=int(checked.sum()),
        n_pairs_excluded=int(excluded.sum()),
        n_atoms_skipped=n_skipped,
        violations=tuple(violations),
        atom_limits=atom_limits,
        ca_clash_distance=ca_clash_distance,
        search_radius=float(search_radius),
    )
