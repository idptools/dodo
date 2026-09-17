"""Tests for reading FASTA and deciding which chain each record belongs to.

The format is trivial; the mapping is not. A FASTA downloaded from the RCSB names each chain
twice -- once as the mmCIF ``label_asym_id`` and once as the author's ``auth_asym_id`` -- and a
biological assembly then renames every chain again for each symmetry copy. DODO keys chains on
the author id, so that is the name that has to win, and a chain no header names has to be
matched on the one thing renaming cannot touch: its sequence.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dodo.exceptions import StructureFileError
from dodo.io.fasta import FastaRecord, match_chains, parse_fasta, read_fasta

DATA = Path(__file__).resolve().parents[1] / "data" / "complexes_missing_residues"
REAL_FASTA = DATA / "7R5J_chain_sequences.fasta"


class TestParsing:
    def test_a_bare_header_names_its_chain(self) -> None:
        (record,) = parse_fasta(">A\nMKVLA\n")
        assert record.chain_ids == ("A",)
        assert record.sequence == "MKVLA"

    def test_a_description_after_the_id_is_not_part_of_it(self) -> None:
        (record,) = parse_fasta(">A my favourite protein\nMKV\n")
        assert record.chain_ids == ("A",)
        assert record.description == "my favourite protein"

    def test_the_author_id_wins_over_the_label_id(self) -> None:
        (record,) = parse_fasta(">7R5J_1|Chain WC[auth W0]|Nup88|Homo sapiens\nMKV\n")
        assert record.chain_ids[0] == "W0"
        assert "WC" in record.chain_ids

    def test_a_multi_chain_header_names_every_chain(self) -> None:
        (record,) = parse_fasta(
            ">7R5J_11|Chains TA[auth I0], UA[auth I1], VA[auth I2]|p58|Homo sapiens\nMKV\n"
        )
        assert record.chain_ids[:3] == ("I0", "I1", "I2")
        assert record.description == "p58"

    def test_chains_without_auth_ids_still_parse(self) -> None:
        (record,) = parse_fasta(">1ABC_1|Chains A, B|Thing|Org\nMKV\n")
        assert record.chain_ids == ("A", "B")

    def test_a_database_header_names_no_chain(self) -> None:
        """Reject a database header.

        ``>sp|P04637|P53_HUMAN`` is a UniProt accession, not a chain id, and guessing that it is
        one would silently assign the wrong reference to a chain called 'sp'.
        """
        (record,) = parse_fasta(">sp|P04637|P53_HUMAN Cellular tumor antigen p53\nMKV\n")
        assert record.chain_ids == ()

    def test_sequences_are_joined_and_cleaned(self) -> None:
        (record,) = parse_fasta(">A\nMKV\nLA-\nGG*\n")
        assert record.sequence == "MKVLAGG"

    def test_a_file_with_no_header_is_refused(self) -> None:
        with pytest.raises(StructureFileError, match="No FASTA records"):
            parse_fasta("MKVLA\n")

    def test_a_header_with_no_sequence_is_refused(self) -> None:
        with pytest.raises(StructureFileError, match="no sequence"):
            parse_fasta(">A\n>B\nMKV\n")

    def test_the_real_file_parses(self) -> None:
        records = read_fasta(REAL_FASTA)
        assert len(records) == 25
        assert all(record.chain_ids for record in records)
        assert all(record.sequence.isalpha() for record in records)


class TestMatching:
    def _records(self) -> list[FastaRecord]:
        return parse_fasta(
            ">1ABC_1|Chains A[auth X], B[auth Y]|First|Org\nMKVLAGGWW\n"
            ">1ABC_2|Chain C[auth Z]|Second|Org\nQQQPPPRRR\n"
        )

    def test_an_exact_author_id_matches(self) -> None:
        sequences, notes = match_chains(self._records(), ["X", "Z"])
        assert sequences == {"X": "MKVLAGGWW", "Z": "QQQPPPRRR"}
        assert notes == []

    def test_a_label_id_matches_too(self) -> None:
        sequences, _ = match_chains(self._records(), ["B"])
        assert sequences == {"B": "MKVLAGGWW"}

    def test_an_assembly_copy_suffix_is_stripped(self) -> None:
        """A biological assembly numbers its symmetry copies ``X-2``, ``X-3``, ..."""
        sequences, notes = match_chains(self._records(), ["X-3"])
        assert sequences == {"X-3": "MKVLAGGWW"}
        assert any("copy suffix" in note for note in notes)

    def test_an_unnamed_chain_matches_on_its_sequence(self) -> None:
        """The fallback that renaming cannot defeat."""
        sequences, notes = match_chains(
            self._records(), ["whatever"], {"whatever": "MKLGW"}
        )
        assert sequences == {"whatever": "MKVLAGGWW"}
        assert any("by sequence" in note for note in notes)

    def test_an_ambiguous_sequence_is_left_unmatched(self) -> None:
        """Between two DIFFERENT references, guessing would rebuild against the wrong one."""
        records = parse_fasta(">1|Chain A|x|o\nMKVLA\n>2|Chain B|y|o\nMKVQQ\n")
        sequences, notes = match_chains(records, ["Q"], {"Q": "MKV"})
        assert sequences == {}
        assert any("two or more DIFFERENT" in note for note in notes)

    def test_identical_references_are_interchangeable(self) -> None:
        """Two chains of the same entity carry the same sequence, so either answer is right."""
        records = parse_fasta(">1|Chain A|x|o\nMKVLA\n>2|Chain B|y|o\nMKVLA\n")
        sequences, _ = match_chains(records, ["Q"], {"Q": "MKV"})
        assert sequences == {"Q": "MKVLA"}

    def test_conflicting_duplicate_chain_ids_are_refused(self) -> None:
        records = parse_fasta(">Chain A\nMAAA\n>Chain A\nMAAT\n")
        with pytest.raises(StructureFileError, match="claimed by two records"):
            match_chains(records, ["A"], {"A": "MAA"})

    def test_duplicate_chain_ids_with_the_same_sequence_are_allowed(self) -> None:
        records = parse_fasta(">1|Chain A|first\nMAAA\n>2|Chain A|second\nMAAA\n")
        sequences, notes = match_chains(records, ["A"], {"A": "MAA"})
        assert sequences == {"A": "MAAA"}
        assert notes == []

    def test_an_unmatched_chain_is_named(self) -> None:
        sequences, notes = match_chains(self._records(), ["nope"], {"nope": "CCCCC"})
        assert sequences == {}
        assert any("no matching FASTA record" in note for note in notes)

    @pytest.mark.skipif(
        not (DATA / "7R5J_subset.cif").exists(),
        reason="7R5J_subset.cif is generated by scripts/make_npc_subset.py from a 597 MB download",
    )
    def test_the_real_file_matches_the_real_structure(self) -> None:
        """Every chain of a real, renamed, biological-assembly slice finds its reference."""
        from dodo.io import read_structure

        structure = read_structure(DATA / "7R5J_subset.cif")
        observed = {chain.chain_id: chain.sequence for chain in structure.chains}
        sequences, _ = match_chains(read_fasta(REAL_FASTA), observed.keys(), observed)
        assert set(sequences) == set(observed)
        for chain_id, reference in sequences.items():
            assert len(reference) >= len(observed[chain_id])
