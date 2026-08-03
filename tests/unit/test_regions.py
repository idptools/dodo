"""Tests for region identification.

The merge tests are the important ones. v1 and the first v2 attempt shared a block-merging
routine with two reproducible defects that between them could dissolve a real folded domain
or swallow a whole IDR, and both were masked by a noisy score -- so improving the metric is
what exposes them. Each has a named test here.
"""

from __future__ import annotations

import numpy as np
import pytest

from dodo.constants import (
    CA_CONTACT_SCORE_THRESHOLD,
    CONTACT_SCORE_THRESHOLD,
    LOOP_CONTACT_CUTOFF,
    MIN_LOOP_LENGTH,
    PLDDT_DISORDER_THRESHOLD,
)
from dodo.exceptions import InvalidRegionError
from dodo.regions.contact import contact_profile, loop_contact_counts, smooth
from dodo.regions.identify import (
    Strategy,
    _absorb_short_gaps,
    assign_regions,
    assign_regions_from_spec,
    find_runs,
    merge_blocks,
)
from dodo.structure import DomainKind, Structure

from .test_structure import make_structure

# ---------------------------------------------------------------------------
# Block merging: the two v1 bugs
# ---------------------------------------------------------------------------


class TestMergeBlocks:
    def test_single_block_survives(self) -> None:
        """A lone candidate block must not vanish.

        v1 iterated ``range(0, len(bounds) - 1)``, whose body never executes for one block,
        so ``final_fds`` stayed empty and a protein with one clean folded domain came out
        classified as entirely disordered -- meaning its real domain was replaced by a
        random walk.
        """
        assert merge_blocks([(40, 100)], max_gap=25) == [(40, 100)]

    def test_gap_before_the_last_block_is_tested(self) -> None:
        """The final gap is not exempt from the merge test.

        v1 guarded with ``if fd_ind < len(bounds) - 2``, so the gap before the last block
        was never examined. These exact bounds collapsed to ``[(40, 433)]``: a 393-residue
        "folded domain" containing a 120-residue IDR that consequently never got rebuilt.
        """
        blocks = [(40, 100), (103, 163), (283, 433)]
        assert merge_blocks(blocks, max_gap=25) == [(40, 163), (283, 433)]

    @pytest.mark.parametrize(
        ("blocks", "expected"),
        [
            ([], []),
            ([(0, 10)], [(0, 10)]),
            ([(0, 10), (12, 20)], [(0, 20)]),
            ([(0, 10), (50, 60)], [(0, 10), (50, 60)]),
            ([(0, 5), (6, 10), (11, 15), (100, 110)], [(0, 15), (100, 110)]),
        ],
    )
    def test_every_cardinality(
        self, blocks: list[tuple[int, int]], expected: list[tuple[int, int]]
    ) -> None:
        """Zero, one and many blocks all behave. Both v1 bugs were cardinality bugs."""
        assert merge_blocks(blocks, max_gap=25) == expected

    def test_gap_exactly_at_the_threshold_merges(self) -> None:
        # 25 residues between blocks, max_gap 25: inclusive, so it merges.
        assert merge_blocks([(0, 10), (35, 45)], max_gap=25) == [(0, 45)]

    def test_gap_one_over_the_threshold_does_not(self) -> None:
        assert merge_blocks([(0, 10), (36, 46)], max_gap=25) == [(0, 10), (36, 46)]

    def test_overlapping_blocks_rejected(self) -> None:
        """Out-of-order input is a caller bug and must not be silently mangled."""
        with pytest.raises(InvalidRegionError, match="sorted and non-overlapping"):
            merge_blocks([(10, 20), (5, 15)], max_gap=25)

    def test_negative_gap_rejected(self) -> None:
        with pytest.raises(ValueError, match="non-negative"):
            merge_blocks([(0, 10)], max_gap=-1)


class TestFindRuns:
    @pytest.mark.parametrize(
        ("mask", "expected"),
        [
            ([], []),
            ([False, False], []),
            ([True, True, True], [(0, 3)]),
            ([False, True, True, False, True], [(1, 3), (4, 5)]),
            ([True, False, True], [(0, 1), (2, 3)]),
            ([True, True, False, True], [(0, 2), (3, 4)]),
        ],
    )
    def test_runs(self, mask: list[bool], expected: list[tuple[int, int]]) -> None:
        """Runs touching either array boundary are found correctly.

        The padding trick removes the boundary special cases that index-based
        implementations habitually get wrong.
        """
        assert find_runs(np.array(mask, dtype=bool)) == expected

    def test_min_length_filter(self) -> None:
        mask = np.array([True, False, True, True, True])
        assert find_runs(mask, min_length=2) == [(2, 5)]


# ---------------------------------------------------------------------------
# Contact scoring
# ---------------------------------------------------------------------------


class TestSmooth:
    def test_window_below_two_is_a_passthrough(self) -> None:
        values = np.array([1.0, 5.0, 2.0])
        assert np.array_equal(smooth(values, 1), values)

    def test_constant_input_is_unchanged(self) -> None:
        values = np.full(20, 3.0)
        assert np.allclose(smooth(values, 7), 3.0)

    def test_edges_are_not_dragged_toward_zero(self) -> None:
        """Edge handling shrinks the window rather than zero-padding.

        The first and last residues of a chain are exactly the terminal-IDR boundaries the
        package cares most about. Zero-padding would depress their scores and manufacture
        disorder that is not in the structure.
        """
        values = np.full(10, 4.0)
        result = smooth(values, 5)
        assert result[0] == pytest.approx(4.0)
        assert result[-1] == pytest.approx(4.0)

    def test_smoothing_reduces_variance(self) -> None:
        rng = np.random.default_rng(0)
        noisy = rng.normal(10.0, 3.0, size=200)
        assert smooth(noisy, 7).std() < noisy.std()

    def test_even_window_is_rounded_up_to_stay_centred(self) -> None:
        values = np.arange(11, dtype=float)
        # A centred average of a linear ramp returns the ramp itself.
        assert smooth(values, 4)[5] == pytest.approx(5.0)


class TestContactProfile:
    def test_counts_neighbouring_residues_not_atom_pairs(self) -> None:
        """The score counts residues via their alpha carbons, so composition cannot bias it.

        A raw atom-pair count scales with the residue's own atom count: measured on a real
        model, within one folded domain, glycine (4 heavy atoms) averaged 292 against a 480
        threshold while Trp/Phe/Tyr (12-14) averaged 943. Every residue has exactly one CA,
        so counting CAs removes the bias by construction.
        """
        structure = make_structure(12)
        profile = contact_profile(structure)
        assert profile.raw.shape == (12,)
        # A straight extended chain has few non-local contacts.
        assert np.all(profile.raw < 5.0)

    def test_score_is_invariant_to_side_chains(self) -> None:
        """The same structure, all-atom or CA-only, must score identically.

        An all-atom score is not comparable across inputs: measured on arf19, the CA-only
        form scored 0.26x the all-atom value, so one threshold could not serve both. DODO has
        to handle full AlphaFold models, experimental structures with unmodelled side chains,
        and its own partly-CA-only output, so invariance here is a correctness requirement,
        not a nicety.
        """
        all_atom = make_structure(30, all_atom=True)
        ca_only = make_structure(30, all_atom=False)
        assert np.allclose(contact_profile(all_atom).raw, contact_profile(ca_only).raw)

    def test_compact_scores_higher_than_extended(self) -> None:
        """The signal must actually distinguish packed from extended."""
        extended = make_structure(30, all_atom=False, spacing=3.8)
        compact = make_structure(30, all_atom=False, spacing=3.8)
        # Collapse the compact one into a tight 3D cluster.
        rng = np.random.default_rng(1)
        compact.xyz[:] = rng.normal(0.0, 4.0, size=compact.xyz.shape)
        assert contact_profile(compact).raw.mean() > contact_profile(extended).raw.mean()

    def test_empty_structure_is_handled(self) -> None:
        structure = make_structure(1, all_atom=False)
        profile = contact_profile(structure)
        assert len(profile) == 1

    def test_sequence_neighbours_are_excluded(self) -> None:
        """Sequence neighbours are in contact by covalent geometry, not by folding."""
        structure = make_structure(10, all_atom=False)
        with_exclusion = contact_profile(structure, exclude_within=2).raw.sum()
        without = contact_profile(structure, exclude_within=0).raw.sum()
        assert with_exclusion < without


class TestLoopContactCounts:
    def test_defaults_to_including_sequence_neighbours(self) -> None:
        """The loop measure keeps the offset the tuned threshold was calibrated against.

        Excluding sequence neighbours drops the folded-domain median from 7-8 to 4, at which
        point the tuned cutoff of 6 flags roughly three quarters of a real folded domain as
        loop -- handing genuine secondary structure to the rebuilder. Discrimination is
        identical either way; only the scale moves. So the default keeps the offset.
        """
        structure = make_structure(20, all_atom=False)
        default = loop_contact_counts(structure)
        excluded = loop_contact_counts(structure, exclude_within=2)
        assert default.sum() > excluded.sum()

    def test_extended_chain_has_only_sequence_neighbours(self) -> None:
        """A straight chain at 3.8 A spacing sees only i+-1 within 7 A.

        i+-2 sits at 7.6 A, just outside the 7 A loop radius -- which is a useful sanity
        check on the radius itself: it is tight enough that an extended backbone registers
        almost no packing.
        """
        structure = make_structure(20, all_atom=False, spacing=3.8)
        counts = loop_contact_counts(structure)
        assert counts[10] == 2, "interior residue sees only i-1 and i+1"
        assert counts[0] == 1, "terminal residue sees only i+1"

    def test_extended_chain_is_loop_like(self) -> None:
        structure = make_structure(20, all_atom=False, spacing=3.8)
        counts = loop_contact_counts(structure)
        assert np.all(counts < LOOP_CONTACT_CUTOFF)


# ---------------------------------------------------------------------------
# Assignment
# ---------------------------------------------------------------------------


#: Phenylalanine's heavy atoms, in PDB order. Used to give the all-atom fixture real side
#: chains: :data:`~dodo.constants.CONTACT_SCORE_THRESHOLD` counts all-atom PAIRS within 8 A, so
#: a backbone-only structure cannot reach it however tightly it is packed -- and a structure
#: with no atom outside N/CA/C/O is CA-only as far as strategy selection is concerned, which is
#: the distinction :meth:`TestAssignRegions.test_auto_picks_density_for_all_atom_input` turns on.
PHE_ATOMS = ("N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ")


def _atom_shell(n_atoms: int, radius: float = 1.5) -> np.ndarray:
    """``n_atoms`` offsets spread evenly over a sphere of ``radius``, deterministically.

    A Fibonacci lattice. The point is only that a residue's non-CA atoms do not sit on top of
    each other or all in one direction, so that pair counts in the fixture's cores behave like
    a packed core's rather than like a line of beads.
    """
    index = np.arange(n_atoms)
    z = 1.0 - 2.0 * (index + 0.5) / n_atoms
    ring = np.sqrt(np.maximum(1.0 - z * z, 0.0))
    phi = np.pi * (3.0 - np.sqrt(5.0)) * index
    return radius * np.stack([ring * np.cos(phi), ring * np.sin(phi), z], axis=1)


def make_two_domain_structure(
    fd1: int = 40, idr: int = 60, fd2: int = 40, *, seed: int = 0, all_atom: bool = False
) -> Structure:
    """Build a structure with two compact blobs joined by an extended linker.

    Built geometrically rather than mocked, so the whole scoring path is exercised.

    With ``all_atom`` every residue carries phenylalanine's eleven heavy atoms instead of a
    lone CA, on the same alpha-carbon layout. That is what the density strategy needs to be
    exercised at all: it counts all-atom pairs against a threshold tuned on packed cores.
    """
    rng = np.random.default_rng(seed)
    coords = []
    # First compact domain: a tight cluster.
    coords.extend(rng.normal(0.0, 5.0, size=(fd1, 3)))
    # Extended linker marching away along +x.
    start = np.array([40.0, 0.0, 0.0])
    coords.extend(start + np.arange(idr)[:, None] * np.array([3.8, 0.0, 0.0]))
    # Second compact domain, far away.
    centre = start + np.array([idr * 3.8 + 40.0, 0.0, 0.0])
    coords.extend(centre + rng.normal(0.0, 5.0, size=(fd2, 3)))

    ca_xyz = np.array(coords)
    n = fd1 + idr + fd2
    per_residue_b = [90.0] * fd1 + [30.0] * idr + [90.0] * fd2

    if not all_atom:
        return Structure.from_atom_records(
            xyz=ca_xyz,
            atom_name=["CA"] * n,
            element=["C"] * n,
            residue_name=["ALA"] * n,
            residue_number=list(range(1, n + 1)),
            chain_id=["A"] * n,
            b_factor=per_residue_b,
            source="synthetic-two-domain",
        )

    offsets = np.zeros((len(PHE_ATOMS), 3))
    offsets[[i for i, name in enumerate(PHE_ATOMS) if name != "CA"]] = _atom_shell(
        len(PHE_ATOMS) - 1
    )
    elements = ["N" if name[0] == "N" else "O" if name[0] == "O" else "C" for name in PHE_ATOMS]
    return Structure.from_atom_records(
        xyz=(ca_xyz[:, None, :] + offsets[None, :, :]).reshape(-1, 3),
        atom_name=list(PHE_ATOMS) * n,
        element=elements * n,
        residue_name=["PHE"] * (n * len(PHE_ATOMS)),
        residue_number=[i for i in range(1, n + 1) for _ in PHE_ATOMS],
        chain_id=["A"] * (n * len(PHE_ATOMS)),
        b_factor=[b for b in per_residue_b for _ in PHE_ATOMS],
        source="synthetic-two-domain-all-atom",
    )


class TestAssignRegions:
    def test_finds_two_domains_and_a_linker(self) -> None:
        structure = make_two_domain_structure()
        assignment = assign_regions(structure, strategy=Strategy.CONTACT)[0]
        assert assignment.n_folded == 2
        assert assignment.n_idrs == 1
        kinds = [d.kind for d in assignment.domains]
        assert kinds == [DomainKind.FOLDED, DomainKind.IDR, DomainKind.FOLDED]

    def test_plddt_strategy_uses_b_factors(self) -> None:
        """Explicitly requested pLDDT still works -- and asking by name is now the only route.

        AUTO never selects it (see the two tests below), so this is the whole of pLDDT's
        remaining route into a build and it has to keep working: the fixture's B-factor column
        is 90 over the domains and 30 over the linker, and thresholding it at
        PLDDT_DISORDER_THRESHOLD must recover both domains.
        """
        structure = make_two_domain_structure()
        assignment = assign_regions(structure, strategy=Strategy.PLDDT)[0]
        assert assignment.strategy is Strategy.PLDDT
        assert assignment.threshold == PLDDT_DISORDER_THRESHOLD
        assert assignment.n_folded == 2

    def test_auto_picks_density_for_all_atom_input(self) -> None:
        """AUTO resolves to the author's density metric when there are side chains to count.

        This replaces an assertion that AUTO prefers pLDDT when the B-factor column looks like
        pLDDT. That preference is reversed: DENSITY -- all-atom pairs within
        CONTACT_RADIUS against CONTACT_SCORE_THRESHOLD -- is the metric DODO was built and
        validated on, and the author reports it outperforms sequence-based disorder predictors,
        so it is the default and pLDDT is an explicit opt-in.

        The fixture's B-factor column still reads exactly like pLDDT (90 over the domains, 30
        over the linker), which is the point: that is no longer evidence for anything.
        """
        structure = make_two_domain_structure(all_atom=True)
        assert structure.b_factor.max() <= 100.0  # premise: still looks like pLDDT
        assignment = assign_regions(structure, strategy=Strategy.AUTO)[0]
        assert assignment.strategy is Strategy.DENSITY
        assert assignment.threshold == CONTACT_SCORE_THRESHOLD
        # Deliberately NOT noted. Density on all-atom input is the expected, validated path that
        # almost every run takes, and announcing it every time trained people to skip the notes --
        # so the note now fires only when auto does something surprising, which is the CA-only
        # fallback to `contact`. Pinned here so it does not creep back.
        assert not any("auto-selected" in note for note in assignment.notes)
        # Resolved to it *and* it works: the density metric draws the two domains and the
        # linker, with the scored profile either side of its own threshold.
        assert [d.kind for d in assignment.domains] == [
            DomainKind.FOLDED,
            DomainKind.IDR,
            DomainKind.FOLDED,
        ]
        idr = next(d for d in assignment.domains if d.kind is DomainKind.IDR)
        assert assignment.score[idr.span.slice].max() < CONTACT_SCORE_THRESHOLD

    def test_auto_picks_contact_for_ca_only_input_and_never_plddt(self) -> None:
        """AUTO falls back to the CA-only score on CA-only input -- not to pLDDT.

        A pair count is not comparable to CONTACT_SCORE_THRESHOLD once the side chains that
        produced those pairs are missing (measured: the same structure stripped to CA-only
        scores 0.26x), so the density metric is not applicable here. The alternative is the
        composition-free CA neighbour count, thresholded at CA_CONTACT_SCORE_THRESHOLD -- not
        pLDDT, however pLDDT-shaped the B-factor column looks.

        "Never pLDDT" is asserted so that it cannot pass by coincidence: the B-factor column is
        INVERTED, confident over the linker and unconfident over the two blobs. Read as pLDDT
        that column calls the linker the folded domain, which is the opposite of the truth and
        of what the geometry says -- so if AUTO were still consulting it, the answer here would
        differ, and it does not.
        """
        structure = make_two_domain_structure()
        structure.b_factor[:] = np.where(structure.b_factor > 50.0, 30.0, 90.0)
        assert structure.b_factor.max() <= 100.0  # premise: still shaped like pLDDT

        misled = assign_regions(structure, strategy=Strategy.PLDDT)[0]
        assert misled.n_folded == 1  # premise: the column, if read, says something else

        assignment = assign_regions(structure, strategy=Strategy.AUTO)[0]
        assert assignment.strategy is Strategy.CONTACT
        assert assignment.threshold == CA_CONTACT_SCORE_THRESHOLD
        # This one IS noted, and it is the only auto case that is. Falling back to `contact`
        # means the validated density threshold could not be applied, which changes where the
        # boundaries land, so a reader needs to know. Contrast
        # :meth:`test_auto_picks_density_for_all_atom_input`, where the note would be noise.
        assert any("auto-selected" in note for note in assignment.notes)
        assert any("no side chains" in note for note in assignment.notes)
        assert assignment.n_folded == 2
        assert [d.kind for d in assignment.domains] == [
            DomainKind.FOLDED,
            DomainKind.IDR,
            DomainKind.FOLDED,
        ]

    def test_auto_picks_contact_when_b_factors_are_crystallographic(self) -> None:
        structure = make_two_domain_structure()
        # Crystallographic B-factors: a plausible range that never approaches 100.
        structure.b_factor[:] = np.linspace(8.0, 45.0, structure.n_residues)
        assignment = assign_regions(structure, strategy=Strategy.AUTO)[0]
        assert assignment.strategy is Strategy.CONTACT

    def test_auto_does_not_claim_plddt_for_a_constant_column(self) -> None:
        """A constant B-factor column carries no information and must not be read as pLDDT."""
        structure = make_two_domain_structure()
        structure.b_factor[:] = 0.0
        assert assign_regions(structure, strategy=Strategy.AUTO)[0].strategy is Strategy.CONTACT

    def test_domains_tile_the_chain_without_gaps(self) -> None:
        structure = make_two_domain_structure()
        assign_regions(structure, strategy=Strategy.CONTACT)
        chain = structure.chains[0]
        chain.validate_domains()  # raises on any gap or overlap
        covered = sum(len(d) for d in chain.domains)
        assert covered == structure.n_residues

    def test_idr_anchors_point_at_the_flanking_residues(self) -> None:
        """An IDR's anchors are the neighbouring fixed residues, outside its own span."""
        structure = make_two_domain_structure()
        assignment = assign_regions(structure, strategy=Strategy.CONTACT)[0]
        idr = next(d for d in assignment.domains if d.kind is DomainKind.IDR)
        assert idr.span.n_anchor == idr.span.start - 1
        assert idr.span.c_anchor == idr.span.stop
        assert not idr.span.is_terminal

    def test_terminal_idr_has_one_free_end(self) -> None:
        structure = make_two_domain_structure(fd1=0, idr=40, fd2=40)
        assignment = assign_regions(structure, strategy=Strategy.CONTACT)[0]
        first = assignment.domains[0]
        if first.kind is DomainKind.IDR:
            assert first.span.n_anchor is None
            assert first.span.is_terminal

    def test_fully_disordered_chain_is_reported(self) -> None:
        structure = make_structure(60, all_atom=False, spacing=3.8)
        assignment = assign_regions(structure, strategy=Strategy.CONTACT)[0]
        assert assignment.fully_disordered
        assert any("no folded domain" in note for note in assignment.notes)
        # Still a valid tiling: one IDR covering everything.
        assert sum(len(d) for d in assignment.domains) == 60

    def test_short_interior_gap_never_becomes_an_idr(self) -> None:
        """A 2-residue gap between domains is not rebuilt.

        With the default knobs this is handled upstream, by the ``max_internal_gap`` merge
        (25 >= 2), rather than by short-gap absorption -- so no absorption note is expected
        here. The outcome is what matters: two residues do not become an IDR.
        """
        structure = make_two_domain_structure(fd1=40, idr=2, fd2=40)
        assignment = assign_regions(structure, strategy=Strategy.CONTACT)[0]
        assert assignment.n_idrs == 0
        assert assignment.n_folded == 1, "the two blocks merge across a 2-residue gap"

    def test_absorption_helper_handles_tails_and_interior_gaps(self) -> None:
        """Short-gap absorption, tested directly on the helper.

        Tested at the unit level rather than end to end because the path is genuinely hard
        to reach through real geometry: a 2-residue tail hanging off a domain is *in contact*
        with that domain, so scoring it folded is correct, and interior gaps are already
        merged upstream by ``max_internal_gap`` (25 >= ``min_idr_length`` 4). The helper still
        has to be right for callers who tighten ``max_internal_gap`` below
        ``min_idr_length``, and for genuinely distant terminal residues.
        """
        blocks, notes = _absorb_short_gaps(
            [(3, 40), (60, 100)], chain_start=0, chain_stop=102, min_idr_length=4
        )
        # 3-residue N-terminal tail and 2-residue C-terminal tail both absorbed;
        # the 20-residue interior gap is left alone as a real IDR.
        assert blocks == [(0, 40), (60, 102)]
        assert any("N-terminal tail" in note for note in notes)
        assert any("C-terminal tail" in note for note in notes)
        assert not any("interior gap" in note for note in notes)

    def test_absorption_leaves_rebuildable_gaps_alone(self) -> None:
        blocks, notes = _absorb_short_gaps(
            [(0, 40), (60, 100)], chain_start=0, chain_stop=100, min_idr_length=4
        )
        assert blocks == [(0, 40), (60, 100)]
        assert notes == []

    def test_rejected_short_folded_block_is_recorded(self) -> None:
        structure = make_two_domain_structure(fd1=40, idr=60, fd2=40)
        assignment = assign_regions(structure, strategy=Strategy.CONTACT, min_folded_length=200)[0]
        assert assignment.fully_disordered
        # One summary line for all rejected blocks, not one line each: a structure with a dozen
        # short folded-looking patches used to print a dozen near-identical notes.
        matching = [n for n in assignment.notes if "minimum for a folded domain" in n]
        assert len(matching) == 1, assignment.notes
        assert "treated as disordered" in matching[0]

    def test_describe_is_one_based(self) -> None:
        """Display is 1-based; internals are 0-based positional. Converted in one place."""
        structure = make_two_domain_structure()
        assignment = assign_regions(structure, strategy=Strategy.CONTACT)[0]
        assert assignment.describe().startswith("chain A: FD 1-")

    def test_score_and_threshold_are_retained_for_audit(self) -> None:
        structure = make_two_domain_structure()
        assignment = assign_regions(structure, strategy=Strategy.CONTACT)[0]
        assert assignment.score.shape == (structure.n_residues,)
        assert assignment.folded_mask.shape == (structure.n_residues,)
        assert np.isfinite(assignment.threshold)

    def test_every_chain_of_a_multi_chain_structure_is_assigned(self) -> None:
        rng = np.random.default_rng(3)
        n_per = 50
        coords = np.concatenate(
            [rng.normal(0.0, 5.0, (n_per, 3)), rng.normal(200.0, 5.0, (n_per, 3))]
        )
        structure = Structure.from_atom_records(
            xyz=coords,
            atom_name=["CA"] * 2 * n_per,
            element=["C"] * 2 * n_per,
            residue_name=["ALA"] * 2 * n_per,
            residue_number=list(range(1, n_per + 1)) * 2,
            chain_id=["A"] * n_per + ["B"] * n_per,
        )
        assignments = assign_regions(structure, strategy=Strategy.CONTACT)
        assert [a.chain_id for a in assignments] == ["A", "B"]


class TestLoopDetection:
    def test_loops_are_strictly_interior(self) -> None:
        """A run touching a domain boundary is a tail, not an anchored loop."""
        structure = make_two_domain_structure()
        assign_regions(structure, strategy=Strategy.CONTACT)
        for domain in structure.domains:
            for loop in domain.loops:
                assert loop.start > domain.span.start
                assert loop.stop < domain.span.stop
                assert loop.span_anchors_present if hasattr(loop, "span_anchors_present") else True
                assert loop.n_anchor is not None
                assert loop.c_anchor is not None

    def test_min_loop_length_uses_inclusive_comparison(self) -> None:
        """``>=``, not ``>``.

        v1 compared with strict ``>``, so the effective minimum was 11 while the
        documentation said 10. The constant now means what it says.
        """
        assert MIN_LOOP_LENGTH == 10
        mask = np.array([False] * 5 + [True] * 10 + [False] * 5)
        assert find_runs(mask, min_length=MIN_LOOP_LENGTH) == [(5, 15)]


class TestManualSpec:
    def test_explicit_regions_are_used_verbatim(self) -> None:
        structure = make_structure(100, all_atom=False)
        assignments = assign_regions_from_spec(
            structure, {"A": [("folded", 1, 40), ("idr", 41, 70), ("folded", 71, 100)]}
        )
        assignment = assignments[0]
        assert [d.kind for d in assignment.domains] == [
            DomainKind.FOLDED,
            DomainKind.IDR,
            DomainKind.FOLDED,
        ]
        # 1-based inclusive input -> 0-based half-open internally.
        assert assignment.domains[0].span.start == 0
        assert assignment.domains[0].span.stop == 40
        assert assignment.domains[1].span.start == 40

    def test_unknown_chain_is_rejected(self) -> None:
        structure = make_structure(10, all_atom=False)
        with pytest.raises(InvalidRegionError, match="not in this structure"):
            assign_regions_from_spec(structure, {"Z": [("folded", 1, 10)]})

    def test_unknown_kind_is_rejected(self) -> None:
        structure = make_structure(10, all_atom=False)
        with pytest.raises(InvalidRegionError, match="Unknown region kind"):
            assign_regions_from_spec(structure, {"A": [("squishy", 1, 10)]})

    def test_out_of_range_is_rejected(self) -> None:
        structure = make_structure(10, all_atom=False)
        with pytest.raises(InvalidRegionError, match="outside chain"):
            assign_regions_from_spec(structure, {"A": [("folded", 1, 50)]})

    def test_gap_in_the_spec_is_rejected(self) -> None:
        """All four of these were accepted silently by the pre-rewrite override path."""
        structure = make_structure(100, all_atom=False)
        with pytest.raises(InvalidRegionError, match="unassigned residues"):
            assign_regions_from_spec(structure, {"A": [("folded", 1, 40), ("folded", 61, 100)]})

    def test_overlap_in_the_spec_is_rejected(self) -> None:
        structure = make_structure(100, all_atom=False)
        with pytest.raises(InvalidRegionError, match="overlapping"):
            assign_regions_from_spec(structure, {"A": [("folded", 1, 60), ("folded", 40, 100)]})
