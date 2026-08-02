"""The broad corpus: 117 AlphaFold structures chosen to span DODO's region topologies.

Why a manifest instead of "grab some structures"
------------------------------------------------
Sampling by file size, which is what an earlier sweep did, gives you whatever the proteome's bulk
happens to be -- mostly single-domain proteins with a terminal tail. It found 50 of 60 structures
completely clean and concluded the pipeline was in good shape.

This corpus was selected by *topology* instead: 4,000 structures were scanned and classified by
whether they contain terminal IDRs, connecting IDRs between folded domains, and loops enclosed
within a single folded domain, then stratified by the length of each. It immediately found three
defects the size-sampled sweep had not, all in structures the random sample had no reason to
include -- a 2,240-residue protein with a very long connecting IDR, one with five folded domains,
and one whose loop cannot be closed because its own anchors are 3.04 A apart.

That is the argument for the manifest: the interesting failures live in the tails.

What it covers, from ``tests/data/corpus.json``
-----------------------------------------------
* terminal IDRs from 4 to 1,729 residues
* connecting IDRs from 26 to 1,612 residues
* loops from 10 to 25 residues (the loop definition itself bounds these)
* 6 structures that are essentially fully folded, where "nothing to rebuild" must be a valid
  answer rather than a crash
* deliberate extremes: five or more folded domains, three or more loops in one domain, very long
  connecting IDRs, and chains under 120 residues

Inherited geometry, and why the assertions are phrased against the input
------------------------------------------------------------------------
Several of these inputs are not clean. AF-Q9BTC0-F1 contains 92 pairs of atoms closer than the
shortest bond in any protein, one at 0.046 A. DODO moves folded domains as rigid bodies and never
regenerates their atoms, so it preserves all of that exactly -- which is correct, and means an
absolute assertion of "no impossible separations in the output" would be a test of AlphaFold, not
of DODO.

So every assertion here is differential: the output must contain no defect that the **input** did
not already have. That is the invariant DODO can actually be held to.

These tests fetch from the AlphaFold database and are marked ``network`` and ``slow``. Fetches are
cached by :func:`dodo.io.fetch_alphafold`, so only the first run downloads anything.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pytest

if TYPE_CHECKING:  # pragma: no cover
    from dodo.structure import Structure

MANIFEST = Path(__file__).resolve().parents[1] / "data" / "corpus.json"


def _load() -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = json.loads(MANIFEST.read_text())
    return records


CORPUS = _load()
IDS = [f"{r['accession']}-{r['stratum'].replace('/', '-')}" for r in CORPUS]


def _clash_keys(structure: Structure) -> set[tuple[str, str, str, str]]:
    from dodo.validate import validate_clashes

    return {
        (v.label_a, v.atom_name_a or "", v.label_b or "", v.atom_name_b or "")
        for v in validate_clashes(structure).violations
    }


def _impossible_keys(structure: Structure) -> set[tuple[str, str, str, str]]:
    from dodo.validate import find_impossible_pairs

    return {
        (p.residue_labels[0], p.atom_names[0], p.residue_labels[1], p.atom_names[1])
        for p in find_impossible_pairs(structure)
    }


class TestCorpusManifest:
    """Guards on the manifest itself. If coverage silently narrows, these fail."""

    def test_size(self) -> None:
        assert len(CORPUS) >= 100, f"the corpus must stay at 100+ structures, got {len(CORPUS)}"

    def test_accessions_are_unique(self) -> None:
        accessions = [r["accession"] for r in CORPUS]
        assert len(set(accessions)) == len(accessions)

    def test_every_region_topology_is_represented(self) -> None:
        assert sum(1 for r in CORPUS if r["terminal"]) >= 20, "too few terminal IDRs"
        assert sum(1 for r in CORPUS if r["connecting"]) >= 20, "too few connecting IDRs"
        assert sum(1 for r in CORPUS if r["loops"]) >= 20, "too few loops"
        fully_folded = [r for r in CORPUS if not (r["terminal"] or r["connecting"] or r["loops"])]
        assert fully_folded, "no fully-folded structure; 'nothing to rebuild' must be covered"

    def test_lengths_span_a_wide_range(self) -> None:
        """Varying lengths were the requirement, so assert the spread rather than mere presence."""
        for key, floor, ceiling in (("terminal", 20, 500), ("connecting", 40, 300)):
            lengths = [n for r in CORPUS for n in r[key]]
            assert min(lengths) <= floor, f"no short {key} region (min {min(lengths)})"
            assert max(lengths) >= ceiling, f"no long {key} region (max {max(lengths)})"
        loops = [n for r in CORPUS for n in r["loops"]]
        assert min(loops) <= 12 and max(loops) >= 20, f"loops span only {min(loops)}-{max(loops)}"


@pytest.mark.network
@pytest.mark.slow
@pytest.mark.parametrize("record", CORPUS, ids=IDS)
class TestCorpusRebuild:
    @pytest.fixture(scope="class")
    def _fetch(self) -> Any:
        from dodo.io import fetch_alphafold

        def fetch(accession: str) -> Path:
            try:
                return fetch_alphafold(accession)
            except Exception as exc:  # a fetch failure must skip, not fail
                pytest.skip(f"could not fetch {accession}: {type(exc).__name__}: {exc}")

        return fetch

    def test_introduces_no_new_defects(self, record: dict[str, Any], _fetch: Any) -> None:
        """The differential invariant: DODO adds no impossible separation and no bad bond.

        Not "the output is clean" -- several of these inputs are not clean, and DODO preserving
        an input defect faithfully is correct behaviour. What DODO is answerable for is anything
        the input did not already contain.
        """
        from dodo.construct.pipeline import rebuild
        from dodo.io import read_structure
        from dodo.validate import validate_bonds

        path = _fetch(record["accession"])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            original = read_structure(path)
            report = rebuild(path, seed=0)
        model = report.models[0]

        new_impossible = _impossible_keys(model) - _impossible_keys(original)
        assert not new_impossible, (
            f"{record['accession']}: DODO introduced {len(new_impossible)} impossible "
            f"separation(s) absent from the input: {sorted(new_impossible)[:3]}"
        )

        generated = [v for v in validate_bonds(model).violations if v.provenance == "rebuilt"]
        assert not generated, (
            f"{record['accession']}: {len(generated)} bad bond(s) in geometry DODO generated: "
            + "; ".join(v.message for v in generated[:3])
        )

    def test_rebuild_completes_and_reports_honestly(
        self, record: dict[str, Any], _fetch: Any
    ) -> None:
        """It builds, and any region it could not build is reported rather than faked."""
        from dodo.construct.pipeline import rebuild

        path = _fetch(record["accession"])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(path, seed=0)

        assert report.models, f"{record['accession']}: produced no model"
        assert report.models[0].n_residues == record["n_residues"], (
            "the manifest's residue count no longer matches the fetched structure; AFDB may have "
            "published a new model version, in which case the manifest needs regenerating"
        )
        for failure in report.failures:
            # A failure is allowed, but it must carry a reason and leave input coordinates.
            assert failure.reason, f"{record['accession']}: a region failed with no reason given"

    def test_failed_regions_keep_every_atom(self, record: dict[str, Any], _fetch: Any) -> None:
        """A region DODO could not build must not lose its side chains.

        Loops used to be stripped to alpha carbons unconditionally, whether or not the build
        succeeded, so a loop DODO declined to touch still came out mutilated -- and then looked
        like DODO's output to the validators.
        """
        from dodo.construct.pipeline import rebuild

        path = _fetch(record["accession"])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(path, seed=0)
        if not report.failures:
            pytest.skip("no region failed for this structure")

        model = report.models[0]
        for failure in report.failures:
            start, stop = failure.residues[0] - 1, failure.residues[1]
            names = set()
            for residue in range(start, stop):
                atoms = model.atom_slice_for_residues(residue, residue + 1)
                names.update(str(n) for n in model.atom_name[atoms])
            assert names > {"CA"}, (
                f"{record['accession']}: region {failure.residues} failed to build but was "
                f"stripped to {sorted(names)} anyway"
            )
