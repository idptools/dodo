"""The one check no exclusion may suppress: atoms closer than any real bond.

Why this module exists separately
---------------------------------
Every structural validator needs exclusions. A clash detector that flags covalently bonded
atoms, or 1-3 pairs held at 2.2 A by a bond angle, or sequence neighbours held close by the
backbone, is useless -- it buries the real findings under thousands of correct ones.

But an exclusion is a statement about *why* two atoms are close, and it silently becomes a
statement that they *may be arbitrarily close*. That is how all three of DODO's validators
independently acquired the same blind spot: a pair of atoms at exactly 0.000 A was invisible to
each of them whenever the pair happened to fall inside an exclusion --

* bond validation skipped any geometric check on residues whose deposited numbering was not
  consecutive, so two alpha carbons at 0.000 A escaped;
* clash detection excluded pairs within three bonds, or adjacent in sequence, unconditionally on
  distance, so 198 of 300 sampled coincident pairs reported zero violations;
* CONECT validation applied a single flat 0.9 A floor to every declared bond, so two alpha
  carbons declared bonded at 1.0 A -- physically impossible -- passed.

Coincident atoms are not a hypothetical. They are the defect that actually shipped: DODO rebuilt
a region's alpha carbons while leaving the rest of each residue behind, and the result was
visible in a viewer as long spurious lines. A validation suite that can be talked out of
reporting a 0.000 A separation is not doing its job.

The principle
-------------
**No exclusion may suppress a physically impossible distance.** Two distinct heavy atoms have a
floor below which no bonding relationship, no strain and no crystallographic oddity can put them.
Below that floor the question "are these two atoms bonded?" is irrelevant; the geometry is broken
either way.

Every validator calls :func:`find_impossible_pairs` *before* applying any exclusion, and merges
its findings unconditionally.

Where the floor comes from
--------------------------
Measured. Across 105,299,848 bonds in 23,587 AlphaFold structures, the shortest bond of any kind
at its 0.1st percentile is ASN CG-OD1 at 1.2026 A, followed by the other carbonyl and
carboxylate C=O bonds at 1.207-1.213 A. So the shortest real heavy-atom separation in a protein
is a carbonyl double bond just above 1.20 A.

:data:`IMPOSSIBLE_SEPARATION` is set below that with margin, so a genuinely short-but-real bond
is never flagged here, while anything actually broken is.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from ..exceptions import GeometryError
from ..structure import Structure

__all__ = [
    "IMPOSSIBLE_SEPARATION",
    "ImpossiblePair",
    "find_impossible_pairs",
]

#: Separation in Angstroms below which two distinct heavy atoms cannot be, whatever their
#: bonding relationship.
#:
#: MEASURED. The shortest bond observed in 105,299,848 measurements is a carbonyl C=O whose 0.1st
#: percentile is 1.2026 A (ASN CG-OD1). 1.00 sits below that with room to spare, so a real bond
#: at the extreme low end of its distribution is never reported here -- this check is for
#: geometry that is broken, not geometry that is tight.
#:
#: Deliberately NOT set close to 1.20. A validator that fires on unusual-but-real geometry gets
#: ignored, and being ignored is the only way this check can fail.
IMPOSSIBLE_SEPARATION: float = 1.00

#: Separation below which two atoms are reported as effectively the same point rather than merely
#: too close. Purely for message wording; both are reported either way.
COINCIDENT_SEPARATION: float = 0.10


@dataclass(frozen=True, slots=True)
class ImpossiblePair:
    """Two atoms closer together than any real bond could put them.

    Attributes
    ----------
    atom_indices
        The two atom indices, lower first.
    residue_indices
        Their residue indices, in the same order.
    residue_labels
        Human-readable residue labels, e.g. ``("A:GLU142", "A:PRO143")``.
    atom_names
        Their atom names.
    separation
        Measured distance in Angstroms.
    coincident
        True if the two atoms are effectively at the same point.
    message
        A complete, self-contained sentence naming the atoms and the distance.
    """

    atom_indices: tuple[int, int]
    residue_indices: tuple[int, int]
    residue_labels: tuple[str, str]
    atom_names: tuple[str, str]
    separation: float
    coincident: bool
    message: str


def find_impossible_pairs(
    structure: Structure,
    *,
    threshold: float = IMPOSSIBLE_SEPARATION,
    include_hydrogens: bool = False,
) -> tuple[ImpossiblePair, ...]:
    """Find every pair of atoms closer than ``threshold``, with no exclusions whatsoever.

    Deliberately has no ``exclude`` parameter. Bonded pairs, 1-3 pairs, sequence neighbours,
    alternate conformers, atoms of the same residue -- none of them are exempt, because none of
    them can legitimately be this close. Adding an exclusion here would reintroduce exactly the
    blind spot this function exists to close.

    Parameters
    ----------
    structure
        The structure to check.
    threshold
        Separation in Angstroms below which a pair is reported. Defaults to
        :data:`IMPOSSIBLE_SEPARATION`. A caller may lower it, but raising it above the shortest
        real bond (1.20 A) would start flagging genuine carbonyl bonds, so that is rejected.
    include_hydrogens
        Include hydrogens. Off by default, since DODO works with heavy atoms and an X-H bond is
        legitimately ~1.0 A -- which is why the default threshold applies to heavy atoms only.

    Returns
    -------
    tuple[ImpossiblePair, ...]
        Every offending pair, ordered by separation, closest first. Empty if the structure is
        sound.

    Raises
    ------
    GeometryError
        If ``threshold`` is not finite and positive, or is high enough to flag real bonds.
    """
    if not np.isfinite(threshold) or threshold <= 0.0:
        raise GeometryError(f"threshold must be finite and positive, got {threshold}.")
    if threshold > 1.20:
        raise GeometryError(
            f"threshold of {threshold} A would flag real carbonyl bonds, whose 0.1st percentile "
            f"is 1.2026 A. This check is for geometry that is broken, not geometry that is "
            f"tight; use the clash validator for contacts above a bond length."
        )

    selected = (
        np.arange(structure.n_atoms)
        if include_hydrogens
        else np.flatnonzero(~np.isin(structure.element, ("H", "D")))
    )
    if selected.size < 2:
        return ()

    coords = structure.xyz[selected]
    # A non-finite coordinate cannot be closer than anything; excluding it here keeps the KD-tree
    # well-defined and is not a loss, since the bond validator reports non-finite atoms itself.
    finite = np.all(np.isfinite(coords), axis=1)
    selected = selected[finite]
    coords = coords[finite]
    if selected.size < 2:
        return ()

    tree = cKDTree(coords)
    found: list[ImpossiblePair] = []
    for local_a, local_b in tree.query_pairs(threshold):
        atom_a = int(selected[local_a])
        atom_b = int(selected[local_b])
        if atom_b < atom_a:
            atom_a, atom_b = atom_b, atom_a
        separation = float(np.linalg.norm(structure.xyz[atom_a] - structure.xyz[atom_b]))
        residue_a = int(structure.residue_index[atom_a])
        residue_b = int(structure.residue_index[atom_b])
        label_a = structure.residue_label(residue_a)
        label_b = structure.residue_label(residue_b)
        name_a = str(structure.atom_name[atom_a])
        name_b = str(structure.atom_name[atom_b])
        coincident = separation < COINCIDENT_SEPARATION

        detail = "occupy the same point" if coincident else f"are {separation:.3f} A apart"
        message = (
            f"{label_a} {name_a} and {label_b} {name_b} {detail}, closer than the shortest bond "
            f"in any protein ({threshold:.2f} A floor; the shortest measured bond is a carbonyl "
            f"C=O at 1.20 A). No bonding relationship makes this possible, so it is reported "
            f"regardless of any exclusion."
        )
        found.append(
            ImpossiblePair(
                atom_indices=(atom_a, atom_b),
                residue_indices=(residue_a, residue_b),
                residue_labels=(label_a, label_b),
                atom_names=(name_a, name_b),
                separation=separation,
                coincident=coincident,
                message=message,
            )
        )

    found.sort(key=lambda pair: pair.separation)
    return tuple(found)
