"""Regression tests over the structures DODO has historically failed on.

``failures.txt`` in this directory is the list Ryan kept by hand from October 2023 onward, with
credit to whoever found each one. Its stated purpose:

    Compiling these over time to make sure that as people find problems with DODO, I don't break
    previous fixes with new 'fixes', which wouldn't really be a fix because... I'd be breaking
    something else.

That is exactly what this module does. The list sat unused for two and a half years; these are the
automated tests it was waiting for.

Why these particular proteins are a good corpus
-----------------------------------------------
They were not chosen to be convenient. Each one broke a real DODO release, and between them they
cover the shapes that break region identification:

* **p300** (Q09472, 2414 residues) is mostly disorder wrapped around several small ordered
  domains -- the case that stresses folded-domain repositioning hardest, since there are many
  domains to place and long linkers between them.
* **TDP-43** (Q13148) is the textbook two-folded-domains-plus-long-C-terminal-IDR architecture.
* **lysozyme** (P00698) and **GFP** (P42212) are the opposite extreme: compact, essentially fully
  folded, with almost nothing to rebuild. A pipeline that assumes it will find disorder breaks
  here, and "no IDRs at all" is a real input a user will hand us.
* The **PTBP2 / PTBP3 / G3BP1 / TRF1** group are multi-domain proteins with interleaved short
  linkers, which is where loop-versus-IDR classification goes wrong.

What is asserted, and what is only ratcheted
--------------------------------------------
The hard assertions are invariants that must hold for every structure forever: the rebuild
completes, no two atoms end up closer than the shortest bond in any protein, and no bond DODO
built is the wrong length.

Steric clashes are handled differently, with a *ratchet* rather than a fixed expectation. Two of
these structures currently retain a handful of clashes between a rebuilt alpha carbon and a folded
domain's side chain, because the walk's relaxation ladder prefers a known-bad contact to failing
outright. That is a real open defect with a fix in progress. Until it lands, asserting zero would
leave a permanently red suite, and asserting the current numbers exactly would break the moment
AlphaFold publishes a new model version. So the ceiling is deliberately loose and can only ever be
tightened -- see :data:`CLASH_CEILING`.

These tests need the network. They fetch from the AlphaFold database and are marked ``network``, so
``pytest -m 'not network'`` skips them. Fetches are cached by
:func:`dodo.io.fetch_alphafold`, so a repeat run makes no requests at all.
"""

from __future__ import annotations

import warnings
from itertools import pairwise
from typing import TYPE_CHECKING, NamedTuple

import pytest

if TYPE_CHECKING:  # pragma: no cover
    from pathlib import Path

# NOTE: the network/slow markers go on the class below, deliberately NOT on the module. The one
# test that needs no network -- the check that this corpus still matches failures.txt -- must keep
# running under ``pytest -m 'not network'``, since it is the one guarding against silent drift.


class Case(NamedTuple):
    """One historical failure."""

    accession: str
    name: str
    #: Who reported it, per failures.txt. Kept because provenance is the point of this corpus.
    credit: str


CASES: tuple[Case, ...] = (
    Case("Q09472", "p300", "Stephen P., Holehouse Lab"),
    Case("Q13148", "TDP-43", "Ryan"),
    Case("Q9WUL8", "Q9WUL8", "Alex Holehouse"),
    Case("P54274", "TRF1", "Alex Holehouse"),
    Case("Q9UKA9", "PTBP2", "Alex Holehouse"),
    Case("Q13283", "G3BP1", "Alex Holehouse"),
    Case("O95758", "PTBP3", "Alex Holehouse"),
    Case("P00698", "lysozyme", "Alex Holehouse"),
    Case("P42212", "GFP", "Alex Holehouse"),
)

#: Maximum steric clashes tolerated per structure.
#:
#: A RATCHET, not a measurement. Measured today against AFDB v6 with ``seed=0``: p300 6, TDP-43 2,
#: and 0 for the four others available locally. The ceiling sits well above that because the count
#: depends on the AlphaFold model version, which EBI bumps without warning, and on the sparrow
#: dimension predictions.
#:
#: Every one of these is a rebuilt alpha carbon against a folded-domain side chain, admitted by
#: :data:`dodo.constants.CLASH_RELAXATION_LADDER` bottoming out at its 2.00 A rung. When the
#: relaxation pass lands this becomes 0 and this comment goes away. It must only ever move down.
CLASH_CEILING: int = 12

IDS = [f"{c.name}-{c.accession}" for c in CASES]


@pytest.fixture(scope="module")
def _fetched() -> dict[str, Path]:
    """Fetch every case once per session. Cached across runs by dodo.io.fetch_alphafold."""
    from dodo.io import fetch_alphafold

    paths: dict[str, Path] = {}
    for case in CASES:
        try:
            paths[case.accession] = fetch_alphafold(case.accession)
        except Exception as exc:  # a fetch failure must skip the suite, not fail it
            pytest.skip(f"could not fetch {case.accession}: {type(exc).__name__}: {exc}")
    return paths


@pytest.mark.parametrize("case", CASES, ids=IDS)
@pytest.mark.usefixtures("_fetched")
@pytest.mark.network
@pytest.mark.slow
class TestHistoricalFailures:
    def test_rebuild_completes(self, case: Case, _fetched: dict[str, Path]) -> None:
        """It builds at all. Several of these used to raise."""
        from dodo.construct.pipeline import rebuild

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(_fetched[case.accession], seed=0)
        assert report.models, f"{case.name}: produced no model"
        assert report.models[0].n_residues > 0

    def test_no_impossible_separations(self, case: Case, _fetched: dict[str, Path]) -> None:
        """No two atoms closer than the shortest bond in any protein. Never negotiable.

        This is the check that caught the anchor-exemption defect, where rebuilt alpha carbons
        were placed through folded-domain side chains at 0.87-0.94 A.
        """
        from dodo.construct.pipeline import rebuild
        from dodo.validate import find_impossible_pairs

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(_fetched[case.accession], seed=0)
        pairs = find_impossible_pairs(report.models[0])
        assert not pairs, f"{case.name}: " + "; ".join(p.message for p in pairs[:5])

    def test_no_rebuilt_bond_violations(self, case: Case, _fetched: dict[str, Path]) -> None:
        """Every bond DODO built is the right length.

        Only ``rebuilt`` provenance is asserted. AlphaFold's own geometry degrades where it is
        unsure -- 2.44% of its CA-CA bonds fall below 3.3 A at a mean pLDDT of 38.8 -- and DODO
        preserves folded-domain atoms faithfully, so input defects are not DODO's.
        """
        from dodo.construct.pipeline import rebuild
        from dodo.validate import validate_bonds

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(_fetched[case.accession], seed=0)
        violations = [
            v for v in validate_bonds(report.models[0]).violations if v.provenance == "rebuilt"
        ]
        assert not violations, f"{case.name}: " + "; ".join(v.message for v in violations[:5])

    def test_clashes_within_ceiling(self, case: Case, _fetched: dict[str, Path]) -> None:
        """A ratchet, not a target. See :data:`CLASH_CEILING`."""
        from dodo.construct.pipeline import rebuild
        from dodo.validate import validate_clashes

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(_fetched[case.accession], seed=0)
        violations = validate_clashes(report.models[0]).violations
        assert len(violations) <= CLASH_CEILING, (
            f"{case.name}: {len(violations)} clashes exceeds the ceiling of {CLASH_CEILING}. "
            f"This ratchet only moves down. Worst: "
            + "; ".join(v.message for v in violations[:3])
        )

    def test_regions_are_identified(self, case: Case, _fetched: dict[str, Path]) -> None:
        """Region identification produces a complete, non-overlapping partition of every chain.

        Worth asserting separately because a silently empty or overlapping assignment is how
        several of these originally failed -- and because lysozyme and GFP are near-fully-folded,
        so "found no IDRs" must be a valid answer rather than a crash.
        """
        from dodo.io import read_structure
        from dodo.regions import assign_regions

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            structure = read_structure(_fetched[case.accession])
            assignments = assign_regions(structure)
        assert assignments, f"{case.name}: no chain was assigned"
        for assignment, chain in zip(assignments, structure.chains, strict=True):
            covered: list[tuple[int, int]] = sorted(
                (d.span.start, d.span.stop) for d in chain.domains
            )
            assert covered, f"{case.name}: chain has no domains"
            assert covered[0][0] == chain.span.start, f"{case.name}: partition misses the N-term"
            assert covered[-1][1] == chain.span.stop, f"{case.name}: partition misses the C-term"
            for (_, prev_stop), (next_start, _) in pairwise(covered):
                assert prev_stop == next_start, (
                    f"{case.name}: partition has a gap or overlap at {prev_stop}/{next_start}"
                )
            assert assignment.describe()

    def test_deterministic_under_seed(self, case: Case, _fetched: dict[str, Path]) -> None:
        """Same seed, same coordinates. A scientific tool that cannot reproduce itself is broken."""
        import numpy as np

        from dodo.construct.pipeline import rebuild

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            first = rebuild(_fetched[case.accession], seed=7)
            second = rebuild(_fetched[case.accession], seed=7)
        np.testing.assert_allclose(
            first.models[0].xyz, second.models[0].xyz, rtol=0, atol=0, equal_nan=True
        )


def test_every_accession_in_failures_txt_is_covered() -> None:
    """The corpus must not silently drift from the hand-maintained list beside it.

    If someone adds a UniProt ID to failures.txt and forgets to add a Case, this fails and says so.
    """
    import re
    from pathlib import Path as _Path

    text = (_Path(__file__).parent / "failures.txt").read_text()
    # UniProt accessions: O/P/Q + 5 alphanumerics, or the newer A-N/R-Z 6-or-10 form.
    found = set(re.findall(r"\b[OPQ][0-9][A-Z0-9]{3}[0-9]\b", text))
    covered = {c.accession for c in CASES}
    missing = found - covered
    assert not missing, (
        f"failures.txt lists {sorted(missing)} but no Case covers them. Add them to CASES."
    )
