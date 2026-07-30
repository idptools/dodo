"""Contact-based scoring of how buried each residue is.

The signal DODO uses to tell a folded domain from a disordered region: a residue packed into
a folded core has many neighbours, a residue dangling in an IDR has few.

Both scores here count **neighbouring residues by their alpha carbons**. That single choice
fixes two independent defects in the pre-rewrite scoring.

Composition bias
----------------
v1 and the first v2 attempt scored a residue by counting *atom pairs* within 8 A -- every
heavy atom of residue *i* against every heavy atom of every other residue -- and thresholded
the raw count at 480. That count scales with residue *i*'s own heavy-atom count, so it
measured composition as much as burial. Measured on a real 1086-residue AlphaFold model,
within a single folded domain:

* correlation between a residue's heavy-atom count and its score: r = 0.65
* every one of the 18 glycines scored below the 480 threshold (mean 292)
* 94% of Trp/Phe/Tyr scored above it (mean 943)

Glycine has 4 heavy atoms and tryptophan 14, so a buried glycine and an exposed tryptophan
could land on the same side of the cutoff for reasons unrelated to folding. Counting
*residues* via their alpha carbons removes the bias entirely: every residue has exactly one
CA, whatever it is made of.

Scale invariance
----------------
An all-atom score is also not comparable between inputs. Measured on arf19, the same
structure stripped to CA-only scored 0.26x its all-atom value, because side chains reach
further than an alpha carbon alone -- so no single threshold can serve both. That matters
here specifically, because DODO must handle full AlphaFold models, experimental structures
with unmodelled side chains, and its *own* output, which is all-atom in folded domains and
partly CA-only in rebuilt regions. A CA-only score is invariant by construction: measured
ratio 1.000 between the two forms of the same structure.

Smoothing, and why it is not cosmetic
-------------------------------------
The per-residue score is noisy enough to fragment a single domain into dozens of
above/below-threshold blocks. That noise had a load-bearing side effect: it *masked* two
domain-merging bugs downstream (see :mod:`dodo.regions.identify`), because a domain arriving
as 45 fragments never reaches the single-block or last-gap code paths where those bugs lived.
Smoothing is therefore inseparable from fixing them -- improving the metric is what triggers
them.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from ..constants import (
    CONTACT_RADIUS,
    CONTACT_SCORE_SMOOTHING_WINDOW,
    LOOP_CONTACT_CUTOFF,
    LOOP_CONTACT_RADIUS,
)
from ..structure import Structure

__all__ = [
    "ContactProfile",
    "ca_neighbour_counts",
    "contact_profile",
    "is_loop_like",
    "loop_contact_counts",
    "smooth",
]


@dataclass(frozen=True, slots=True)
class ContactProfile:
    """Per-residue burial scores for one structure.

    Attributes
    ----------
    raw
        Non-local CA neighbour count per residue, ``(n_residues,)``.
    smoothed
        ``raw`` after a centred moving average. This is what gets thresholded.
    window
        The smoothing window actually used, in residues.
    radius
        The contact radius used, in Angstroms.
    """

    raw: np.ndarray
    smoothed: np.ndarray
    window: int
    radius: float

    def __len__(self) -> int:
        return int(self.raw.shape[0])


def smooth(values: np.ndarray, window: int) -> np.ndarray:
    """Centred moving average, with edges handled by shrinking the window.

    Edge handling matters here: the first and last residues of a chain are exactly the
    terminal-IDR boundaries DODO cares most about getting right, and zero-padding would drag
    their scores toward zero and manufacture disorder that is not in the structure. Shrinking
    the window instead means an edge residue is averaged over only the neighbours it has.

    Parameters
    ----------
    values
        The sequence to smooth, ``(n,)``.
    window
        Window width in residues. Values below 2 return ``values`` unchanged. An even window
        is rounded up to the next odd number so the average stays centred.
    """
    n = values.shape[0]
    if window < 2 or n == 0:
        return values.astype(np.float64, copy=True)
    if window % 2 == 0:
        window += 1
    half = window // 2

    # Cumulative sums give an O(n) moving average with exact edge normalization.
    padded = np.concatenate([[0.0], np.cumsum(values, dtype=np.float64)])
    starts = np.maximum(np.arange(n) - half, 0)
    stops = np.minimum(np.arange(n) + half + 1, n)
    totals = padded[stops] - padded[starts]
    averaged: np.ndarray = totals / (stops - starts)
    return averaged


def ca_neighbour_counts(
    structure: Structure,
    *,
    radius: float,
    exclude_within: int,
) -> np.ndarray:
    """Count non-local CA neighbours within ``radius`` of each residue's alpha carbon.

    The shared primitive behind both scores in this module. A single KD-tree pair query,
    which is what makes it tractable on large structures: v1's equivalent was a pure-Python
    O(N^2 A^2) double loop that took 10.1 s on a 1086-residue model, and the first v2 attempt
    vectorized it with a dense ``pdist`` that allocates an N-squared matrix (about 500 MB at
    1000 residues, 2 GB at 2000), trading a time problem for a memory ceiling on exactly the
    large assemblies this package targets. A KD-tree has neither problem.

    Parameters
    ----------
    structure
        The structure to score.
    radius
        Contact radius in Angstroms.
    exclude_within
        Residues within this sequence separation are not counted. See the callers for why
        the two scores make different choices here.

    Returns
    -------
    np.ndarray
        Neighbour count per residue, ``(n_residues,)`` int64.
    """
    if structure.n_residues == 0:
        return np.zeros(0, dtype=np.int64)

    tree = cKDTree(structure.ca_xyz)
    pairs = tree.query_pairs(radius, output_type="ndarray")

    counts = np.zeros(structure.n_residues, dtype=np.int64)
    if pairs.size:
        keep = np.abs(pairs[:, 0] - pairs[:, 1]) > exclude_within
        if np.any(keep):
            # query_pairs reports each pair once, so credit both residues.
            contributors = np.concatenate([pairs[keep, 0], pairs[keep, 1]])
            counts += np.bincount(contributors, minlength=structure.n_residues)
    return counts


def contact_profile(
    structure: Structure,
    *,
    radius: float = CONTACT_RADIUS,
    window: int = CONTACT_SCORE_SMOOTHING_WINDOW,
    exclude_within: int = 2,
) -> ContactProfile:
    """Compute per-residue burial scores for a structure.

    Parameters
    ----------
    structure
        The structure to score.
    radius
        CA-CA contact radius in Angstroms. See
        :data:`~dodo.constants.CONTACT_RADIUS` for how 12 A was chosen and why larger
        values are not better.
    window
        Smoothing window in residues.
    exclude_within
        Residues within this sequence separation are not counted as contacts. Sequence
        neighbours are in contact by covalent geometry regardless of folding, so counting
        them adds a near-constant offset to every residue and dilutes the signal. The
        pre-rewrite code excluded only the residue itself.

    Returns
    -------
    ContactProfile
        Raw and smoothed scores, with the parameters used.
    """
    raw = ca_neighbour_counts(structure, radius=radius, exclude_within=exclude_within).astype(
        np.float64
    )
    return ContactProfile(
        raw=raw,
        smoothed=smooth(raw, window),
        window=window,
        radius=radius,
    )


def loop_contact_counts(
    structure: Structure,
    *,
    radius: float = LOOP_CONTACT_RADIUS,
    exclude_within: int = 0,
) -> np.ndarray:
    """Count CA neighbours per residue, for detecting loops inside folded domains.

    A separate, tighter measure from :func:`contact_profile`. The physical basis: a CA in a
    packed core has roughly 6-10 other CAs within 7 A, while a CA in an extended surface loop
    has 2-4. That gap is what :data:`~dodo.constants.LOOP_CONTACT_CUTOFF` thresholds. The
    tighter radius is deliberate -- a loop is about *local* backbone packing, not about how
    much of the domain happens to be nearby.

    Why ``exclude_within`` defaults to 0 here, unlike :func:`contact_profile`
    ------------------------------------------------------------------------
    Sequence-adjacent CAs are always within 7 A by covalent geometry, so they contribute a
    near-constant offset of about 4 to every residue's count. Removing that offset does not
    improve discrimination -- measured on two structures, the balanced separation between
    folded-domain and IDR residues is identical (95% and 97%) either way -- it only shifts
    the scale.

    But it shifts it by enough to invalidate the threshold.
    :data:`~dodo.constants.LOOP_CONTACT_CUTOFF` was tuned against counts that *included*
    sequence neighbours, and excluding them drops the folded-domain median from 7-8 to 4, at
    which point a cutoff of 6 flags roughly three quarters of a real folded domain as loop --
    handing genuine secondary structure to the rebuilder. Measured fraction of a folded
    domain flagged loop-like at ``exclude_within=0`` (arf19 / dnmt3a):

        cutoff 4:  2.1% / 0.7%     cutoff 6: 19.1% / 13.9%   <- the tuned value
        cutoff 5: 10.9% / 6.6%     cutoff 7: 33.8% / 30.3%

    Real domains run about 15-30% loop content, so the tuned 6 is in range and 7 begins
    eating secondary structure. Keeping the offset keeps the constant meaningful.

    Returns
    -------
    np.ndarray
        Neighbour count per residue, ``(n_residues,)`` int64.
    """
    counts: np.ndarray = ca_neighbour_counts(
        structure, radius=radius, exclude_within=exclude_within
    )
    return counts


def is_loop_like(counts: np.ndarray, *, cutoff: int = LOOP_CONTACT_CUTOFF) -> np.ndarray:
    """Return a boolean mask of residues whose CA neighbour count is below ``cutoff``."""
    return counts < cutoff
