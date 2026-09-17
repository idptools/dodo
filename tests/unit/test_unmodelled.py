"""Tests for filling in the residues a structure did not model.

The strongest test available here is a **round trip**: take a complete structure, delete known
residues, hand back the original sequence, and check that exactly those residues come back, in
the right places, with the right numbering. Anything the mapping gets wrong shows up as a
sequence that does not match, and the fixtures carve out one of each case the splicer has to
handle -- a leading tail, an internal gap inside a disordered region, an internal gap inside a
folded domain, and a trailing tail.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dodo.construct.unmodelled import (
    MAX_REFERENCE_MISMATCH_FRACTION,
    insert_unmodelled_residues,
    map_to_reference,
    skipped_for_length,
)
from dodo.exceptions import InvalidRegionError
from dodo.io import read_structure
from dodo.regions.identify import assign_regions
from dodo.structure import DomainKind, Span

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "structures"

#: One of each case the splicer has to handle, as residue index ranges into dnmt3a.pdb:
#: a leading tail, a gap inside a disordered linker, a gap strictly inside a folded domain,
#: and a trailing tail.
CARVED = ((0, 40), (440, 460), (600, 612), (887, 912))


def _carved():
    """dnmt3a with CARVED removed, plus the sequence it should be restored to."""
    full = read_structure(FIXTURES / "dnmt3a.pdb")
    drop = np.zeros(full.n_residues, dtype=bool)
    for start, stop in CARVED:
        drop[start:stop] = True
    partial = full.select_residues(~drop)
    assign_regions(partial)
    return partial, full.sequence


class TestMapping:
    def test_author_numbering_is_used_when_it_explains_the_mapping(self) -> None:
        partial, reference = _carved()
        mapping = map_to_reference(
            partial.sequence, partial.residue_number, reference
        )
        assert mapping.method == "author numbering"
        assert mapping.mismatches == 0
        assert mapping.shift == 0
        # The mapping IS the numbering, so it can be checked independently of the code above.
        assert np.array_equal(mapping.positions, partial.residue_number - 1)

    def test_a_constant_offset_is_found(self) -> None:
        """A cleaved tag numbers the construct from something other than 1."""
        partial, reference = _carved()
        shifted = partial.residue_number + 1000
        mapping = map_to_reference(partial.sequence, shifted, reference)
        assert mapping.method == "author numbering"
        assert mapping.shift == 1000

    def test_alignment_takes_over_when_the_numbering_is_meaningless(self) -> None:
        partial, reference = _carved()
        nonsense = np.arange(1, partial.n_residues + 1) * 3
        mapping = map_to_reference(partial.sequence, nonsense, reference)
        assert mapping.method == "alignment"
        assert mapping.mismatches == 0
        assert "".join(reference[p] for p in mapping.positions) == partial.sequence

    def test_alignment_reports_a_genuine_mismatch(self) -> None:
        partial, reference = _carved()
        mutated = list(reference)
        target = int(partial.residue_number[100]) - 1
        mutated[target] = "W" if mutated[target] != "W" else "A"
        nonsense = np.arange(1, partial.n_residues + 1) * 3
        mapping = map_to_reference(partial.sequence, nonsense, "".join(mutated))
        assert mapping.mismatches >= 1

    def test_a_reference_shorter_than_the_structure_is_refused(self) -> None:
        partial, _ = _carved()
        with pytest.raises(InvalidRegionError, match="at least as long"):
            map_to_reference(partial.sequence, partial.residue_number, "MKV")


class TestRoundTrip:
    def test_every_carved_residue_comes_back(self) -> None:
        partial, reference = _carved()
        filled, report = insert_unmodelled_residues(partial, {"A": reference})
        assert filled.sequence == reference
        assert int(filled.inserted.sum()) == sum(stop - start for start, stop in CARVED)
        assert report.n_inserted == int(filled.inserted.sum())

    def test_the_original_numbering_is_kept(self) -> None:
        partial, reference = _carved()
        filled, _ = insert_unmodelled_residues(partial, {"A": reference})
        assert np.array_equal(filled.residue_number, np.arange(1, len(reference) + 1))

    def test_observed_atoms_are_untouched(self) -> None:
        """Insertion adds; it must not disturb a single observed coordinate."""
        partial, reference = _carved()
        filled, _ = insert_unmodelled_residues(partial, {"A": reference})
        observed = ~filled.inserted
        counts = np.diff(filled.residue_atom_offsets)
        assert int(counts[observed].sum()) == partial.n_atoms
        atom_observed = np.repeat(observed, counts)
        assert np.array_equal(filled.xyz[atom_observed], partial.xyz)

    def test_structure_metadata_survives_insertion(self) -> None:
        partial, reference = _carved()
        partial.experimental_method = "X-RAY DIFFRACTION"
        partial.notes = ["source note"]
        filled, _ = insert_unmodelled_residues(partial, {"A": reference})
        assert filled.experimental_method == partial.experimental_method
        assert filled.notes == partial.notes
        assert filled.notes is not partial.notes

    def test_an_inserted_residue_carries_one_placeholder_alpha_carbon(self) -> None:
        partial, reference = _carved()
        filled, _ = insert_unmodelled_residues(partial, {"A": reference})
        counts = np.diff(filled.residue_atom_offsets)
        assert set(counts[filled.inserted].tolist()) == {1}
        atom_inserted = np.repeat(filled.inserted, counts)
        assert set(filled.atom_name[atom_inserted].tolist()) == {"CA"}
        # Occupancy zero is the crystallographic way of saying "not observed".
        assert float(filled.occupancy[filled.inserted].max()) == 0.0

    def test_nothing_to_do_is_not_an_error(self) -> None:
        structure = read_structure(FIXTURES / "dnmt3a.pdb")
        assign_regions(structure)
        filled, report = insert_unmodelled_residues(structure, {"A": structure.sequence})
        assert report.n_inserted == 0
        assert filled.n_residues == structure.n_residues


class TestSplicing:
    def test_a_gap_inside_a_folded_domain_becomes_a_loop(self) -> None:
        """It must not split the domain: the two halves would then be free to move apart."""
        partial, reference = _carved()
        filled, _ = insert_unmodelled_residues(partial, {"A": reference})
        folded = [d for d in filled.chains[0].domains if d.kind is DomainKind.FOLDED]
        assert folded
        covered = np.zeros(filled.n_residues, dtype=bool)
        for domain in folded:
            for loop in domain.loops:
                covered[loop.slice] = True
        inside = filled.inserted & _inside_any(folded, filled.n_residues)
        assert inside.any()
        assert bool(np.all(covered[inside])), "an inserted run inside a domain is not a loop"

    def test_every_loop_has_two_anchors(self) -> None:
        """A loop is built between two FIXED residues, so its span must carry both anchors.

        Without them the builder rejects it outright: "needs two anchors but only has one; not
        a loop".
        """
        partial, reference = _carved()
        filled, _ = insert_unmodelled_residues(partial, {"A": reference})
        for domain in filled.chains[0].domains:
            for loop in domain.loops:
                assert loop.n_anchor == loop.start - 1
                assert loop.c_anchor == loop.stop

    def test_terminal_runs_extend_the_disordered_regions(self) -> None:
        partial, reference = _carved()
        filled, _ = insert_unmodelled_residues(partial, {"A": reference})
        domains = sorted(filled.chains[0].domains, key=lambda d: d.span.start)
        assert domains[0].kind is DomainKind.IDR
        assert domains[0].span.start == 0
        assert domains[-1].kind is DomainKind.IDR
        assert domains[-1].span.stop == filled.n_residues

    def test_the_tiling_stays_valid(self) -> None:
        partial, reference = _carved()
        filled, _ = insert_unmodelled_residues(partial, {"A": reference})
        filled.validate()
        for chain in filled.chains:
            chain.validate_domains()

    def test_regions_are_never_identified_from_placeholders(self) -> None:
        """Insertion needs regions already assigned, and refuses to run without them.

        Burial is scored from coordinates, and a placeholder is a straight line -- through a
        protein core it scores as buried, so letting it be scored would let fiction decide which
        residues are folded.
        """
        partial, reference = _carved()
        for chain in partial.chains:
            chain.domains = []
        with pytest.raises(InvalidRegionError, match="assigned regions"):
            insert_unmodelled_residues(partial, {"A": reference})


class TestRefusals:
    def test_the_wrong_protein_is_declined(self) -> None:
        partial, reference = _carved()
        scrambled = "".join("W" if r != "W" else "A" for r in reference)
        _filled, report = insert_unmodelled_residues(partial, {"A": scrambled})
        assert report.n_inserted == 0
        assert any("disagree with the reference" in note for note in report.notes)
        assert MAX_REFERENCE_MISMATCH_FRACTION < 1.0

    def test_a_shorter_reference_is_declined_with_a_reason(self) -> None:
        partial, _ = _carved()
        _filled, report = insert_unmodelled_residues(partial, {"A": "MKVLA"})
        assert report.n_inserted == 0
        assert any("not this chain's sequence" in note for note in report.notes)

    def test_an_unbridgeable_gap_is_reported_not_attempted(self) -> None:
        """An impossible gap is named, not filled.

        One residue cannot span 60 A, so inserting it would only turn a fact about the input
        into a build failure. The gap is left empty and reported instead.
        """
        from dodo.engines.walk import max_reach

        full = read_structure(FIXTURES / "dnmt3a.pdb")
        cut_start, cut_stop = 440, 470
        drop = np.zeros(full.n_residues, dtype=bool)
        drop[cut_start:cut_stop] = True
        partial = full.select_residues(~drop)
        assign_regions(partial)

        separation = float(
            np.linalg.norm(full.ca_xyz[cut_stop] - full.ca_xyz[cut_start - 1])
        )
        # A reference claiming exactly ONE missing residue where thirty were removed: two bonds
        # to cross a distance two bonds cannot cross.
        assert separation > max_reach(2), "the fixture no longer makes an impossible gap"
        reference = full.sequence[:cut_start] + "G" + full.sequence[cut_stop:]

        filled, report = insert_unmodelled_residues(partial, {"A": reference})
        assert report.impossible_gaps, "an impossible gap was filled instead of being reported"
        for _chain_id, (_first, _last, span, ceiling) in report.impossible_gaps:
            assert span > ceiling
        assert int(filled.inserted.sum()) == 0
        assert "no conformation exists" in report.summary()

    def test_a_too_close_gap_reports_the_minimum_reach(self) -> None:
        from dodo.construct.unmodelled import ChainInsertion, InsertionReport

        report = InsertionReport(
            chains=[
                ChainInsertion(
                    chain_id="A",
                    observed=2,
                    reference=3,
                    inserted=0,
                    leading=0,
                    trailing=0,
                    internal_gaps=1,
                    method="probe",
                    mismatches=0,
                    impossible_gaps=((10, 12, 1.0, 5.4),),
                )
            ]
        )
        summary = report.summary()
        assert "span at least 5.4 A" in summary
        assert "span at most 5.4 A" not in summary


class TestSkipPredicate:
    def test_a_short_observed_region_is_skipped(self) -> None:
        structure = read_structure(FIXTURES / "dnmt3a.pdb")
        assert skipped_for_length(structure, Span(10, 12), 4)
        assert not skipped_for_length(structure, Span(10, 20), 4)

    def test_an_inserted_region_is_never_skipped(self) -> None:
        """It has no input coordinates to leave alone, only a placeholder."""
        structure = read_structure(FIXTURES / "dnmt3a.pdb")
        structure.inserted[10:12] = True
        assert not skipped_for_length(structure, Span(10, 12), 4)


def _inside_any(domains, n_residues: int) -> np.ndarray:
    mask = np.zeros(n_residues, dtype=bool)
    for domain in domains:
        mask[domain.span.start + 1 : domain.span.stop - 1] = True
    return mask
