"""Every code example in the shipped documentation must actually run.

Not a stylistic check. Three examples were broken in ways a reader would hit immediately by
copying them, and nothing noticed because nothing ran them:

* ``target_dimensions("GSGSGSGS...", mode="compact")`` -- the literal ``...`` is not a sequence, so
  it raised ``ValueError`` instead of printing the result the comment claimed;
* the preset-regions example specified regions covering residues 1-290 of a 912-residue chain, so
  it raised ``InvalidRegionError`` -- region specs must tile the whole chain, which the surrounding
  prose did not say either;
* one block used ``dodo.`` with no ``import dodo`` of its own.

A fourth was subtler: the example's own comment claimed an output of ``62.6 A (0.7x of 89.4 A)``
where the real result is ``55.9 A (0.7x of 79.8 A)``. A stale number in a docstring is a small lie
that a reader has no way to detect, so the asserted outputs are checked too where they are exact.
"""

from __future__ import annotations

import contextlib
import io
import re
import shutil
import warnings
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURES = REPO_ROOT / "tests" / "data" / "structures"

#: Files whose fenced python blocks are executed.
DOCUMENTS = ("README.md", "docs/guide.md", "docs/index.md")

#: Blocks that are attribute listings rather than runnable programs. Matched on their first token,
#: so a block that becomes runnable stops being skipped automatically.
_LISTING_PREFIXES = ("report.ok", "report.failures", "structure.xyz")


def _python_blocks(document: str) -> list[str]:
    text = (REPO_ROOT / document).read_text()
    return re.findall(r"```python\n(.*?)```", text, re.S)


@pytest.fixture
def sandbox(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Create a directory containing the filenames the examples refer to."""
    shutil.copy(FIXTURES / "dnmt3a.pdb", tmp_path / "model.pdb")
    shutil.copy(FIXTURES / "arf19.cif", tmp_path / "model.cif")
    shutil.copy(FIXTURES / "dnmt3a.pdb", tmp_path / "AF-P04637-F1-model_v6.pdb")
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.mark.slow
@pytest.mark.parametrize("document", DOCUMENTS)
def test_every_python_example_runs(document: str, sandbox: Path) -> None:
    """Blocks share one namespace, which is how a reader meets them: in order, top to bottom."""
    blocks = _python_blocks(document)
    assert blocks, f"{document} has no python examples; did the fence style change?"

    namespace: dict[str, object] = {}
    for index, block in enumerate(blocks, start=1):
        if block.strip().startswith(_LISTING_PREFIXES):
            continue
        try:
            with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
                warnings.simplefilter("ignore")
                exec(compile(block, f"<{document} block {index}>", "exec"), namespace)
        except Exception as exc:  # the point is to report WHICH block failed, and why
            pytest.fail(f"{document} block {index} raised {type(exc).__name__}: {exc}\n\n{block}")


@pytest.mark.slow
def test_the_documented_target_dimension_is_still_what_dodo_produces() -> None:
    """The one example whose exact output is quoted in a comment.

    It read 62.6 A against a real 55.9 A. Numbers in documentation are load-bearing for trust, so
    this one is pinned rather than left to drift again.
    """
    from dodo.construct import target_dimensions

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        target = target_dimensions("GSGSGSGS" * 18, mode="compact")

    rendered = str(target)
    for document in ("README.md", "docs/guide.md"):
        text = (REPO_ROOT / document).read_text()
        quoted = re.search(r"print\(target\)\s*#\s*(.+)", text)
        assert quoted, f"{document} no longer quotes the target_dimensions output"
        assert quoted.group(1).strip() == rendered, (
            f"{document} claims {quoted.group(1).strip()!r} but DODO produces {rendered!r}"
        )


def test_no_document_advertises_a_removed_extra() -> None:
    """No user-facing extras ship now. Removed ones must not reappear in prose."""
    removed = ("[predictors]", "[lookup]", "[viz]", "[all]", "[albatross]")
    for document in (*DOCUMENTS, "CHANGELOG.md", "docs/algorithm.md", "docs/validation.md"):
        text = (REPO_ROOT / document).read_text()
        for extra in removed:
            install = f'pip install "idptools-dodo{extra}"'
            assert install not in text, f"{document} still tells users to run {install}"
