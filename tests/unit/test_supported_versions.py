"""The supported-version claims must agree with each other, and with what CI tests.

DODO states which Python versions it supports in three separate places, and they drifted: the CI
matrix omitted 3.11, which is the version the package is actually developed on -- so the one leg
whose result was already known was the one nobody ran. Nothing detected that, because nothing
compared the three lists.

Why Python 3.9 gets a test of its own
-------------------------------------
Python 3.9 cannot run this package. Counted by AST over every source module: 51
``dataclass(slots=True)``, 10 ``zip(..., strict=True)`` and 3 ``itertools.pairwise`` imports, all
3.10-or-later, plus a runtime ``X | Y`` in an ``isinstance``. ``requires-python = ">=3.10"``
declares that, and pip enforces it -- a real install into a 3.9 environment fails with
``ERROR: Package 'idptools-dodo' requires a different Python: 3.9.12 not in '>=3.10'``.

The trap, and the reason this file exists rather than a CI smoke test: **``import dodo`` SUCCEEDS
on Python 3.9.** The package resolves everything lazily through PEP 562 ``__getattr__``, so the
import itself touches nothing that 3.9 lacks. ``python3.9 -m compileall src/dodo`` also exits 0,
because none of the 3.10 features are *syntax*. The first real failure is on attribute access:
``dodo.Structure`` raises ``ImportError: cannot import name 'pairwise' from 'itertools'``.

So the obvious test -- "prove the package fails to import on 3.9" -- would pass while proving
nothing, and would keep passing if the floor were accidentally lowered. What is actually worth
asserting is the metadata: that the declared floor exists, excludes 3.9, and matches everything
else. That is checkable from any interpreter.
"""

from __future__ import annotations

import re
from itertools import pairwise
from pathlib import Path

import pytest
import tomllib

REPO_ROOT = Path(__file__).resolve().parents[2]
PYPROJECT = REPO_ROOT / "pyproject.toml"
CI_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "CI.yaml"

#: The floor. Below this the package cannot run, not merely "is not tested".
MINIMUM_PYTHON = (3, 10)


@pytest.fixture(scope="module")
def metadata() -> dict:
    return tomllib.loads(PYPROJECT.read_text())


def _classifier_versions(metadata: dict) -> set[tuple[int, int]]:
    found = set()
    for classifier in metadata["project"]["classifiers"]:
        match = re.fullmatch(r"Programming Language :: Python :: (\d+)\.(\d+)", classifier)
        if match:
            found.add((int(match.group(1)), int(match.group(2))))
    return found


def _ci_matrix_versions() -> set[tuple[int, int]]:
    """Parse the matrix without a YAML dependency, which the test extra does not include."""
    text = CI_WORKFLOW.read_text()
    match = re.search(r"python-version:\s*\[([^\]]+)\]", text)
    assert match, "could not find the CI python-version matrix"
    return {(int(a), int(b)) for a, b in re.findall(r'"(\d+)\.(\d+)"', match.group(1))}


class TestTheThreeListsAgree:
    def test_requires_python_is_the_declared_floor(self, metadata: dict) -> None:
        requires = metadata["project"]["requires-python"]
        expected = f">={MINIMUM_PYTHON[0]}.{MINIMUM_PYTHON[1]}"
        assert requires == expected, f"requires-python is {requires!r}, expected {expected!r}"

    def test_classifiers_start_at_the_floor_and_are_contiguous(self, metadata: dict) -> None:
        versions = sorted(_classifier_versions(metadata))
        assert versions, "no 'Programming Language :: Python :: X.Y' classifiers"
        assert versions[0] == MINIMUM_PYTHON, f"classifiers start at {versions[0]}, not the floor"
        for earlier, later in pairwise(versions):
            assert later == (earlier[0], earlier[1] + 1), f"gap between {earlier} and {later}"

    def test_ci_tests_every_version_the_classifiers_claim(self, metadata: dict) -> None:
        """The drift this file was written for: 3.11 was claimed and never tested."""
        claimed = _classifier_versions(metadata)
        tested = _ci_matrix_versions()
        assert not claimed - tested, (
            f"claimed but never tested in CI: {sorted(claimed - tested)}. Either add the version "
            f"to the CI matrix or stop advertising it."
        )
        assert not tested - claimed, (
            f"tested in CI but not advertised: {sorted(tested - claimed)}. Add the classifier."
        )


class TestPython39IsExcluded:
    def test_the_floor_excludes_3_9(self, metadata: dict) -> None:
        """Not a style preference. 3.9 lacks language features this package uses everywhere."""
        requires = metadata["project"]["requires-python"]
        assert requires == ">=3.10"
        assert "3.9" not in {f"{major}.{minor}" for major, minor in _classifier_versions(metadata)}

    def test_the_features_that_make_3_10_the_floor_are_really_used(self) -> None:
        """Guards the justification, so the floor cannot be lowered without this failing.

        If someone removes the last ``zip(strict=True)`` and the last ``pairwise``, the floor
        becomes a free choice again and this test should be revisited deliberately -- rather than
        the docstring above quietly becoming false.
        """
        sources = list((REPO_ROOT / "src" / "dodo").rglob("*.py"))
        assert sources
        text = "\n".join(path.read_text() for path in sources)
        assert text.count("slots=True") >= 20, "dataclass(slots=True) is a 3.10+ feature"
        assert "strict=True" in text, "zip(..., strict=True) is a 3.10+ feature"
        assert "from itertools import pairwise" in text, "pairwise is a 3.10+ feature"


#: Dependencies declared as a direct reference (``name @ url``) rather than a version range.
#:
#: These are exempt from the floor rules below, because a direct reference names a repository
#: rather than a version and there is no floor to declare or to install. Listing them here
#: rather than pattern-matching them away is deliberate: a new one has to be added on purpose,
#: which is the moment to notice that it also blocks publishing to PyPI.
#:
#: ``sparrow`` is not on PyPI *yet* -- only ``git+https://github.com/idptools/sparrow.git`` --
#: so this is not a convenience, it is the only way to declare it. It requires
#: ``tool.hatch.metadata.allow-direct-references``; without that every build fails outright.
#: DODO therefore installs from GitHub for now. When sparrow is published, this set becomes
#: empty and the hatch flag can go with it.
KNOWN_DIRECT_REFERENCES = frozenset({"sparrow"})


def _requirement_name(requirement: str) -> str:
    """Return the distribution name from a requirement, direct reference or not."""
    return re.split(r"[><=!\[@]", requirement, maxsplit=1)[0].strip()


class TestDependencyFloorsAreDeclaredAsRanges:
    def test_direct_references_are_the_ones_we_know_about(self, metadata: dict) -> None:
        """A direct reference cannot carry a floor, so it must be an acknowledged exception.

        It also cannot be uploaded to PyPI: the index rejects metadata containing one. So a new
        direct reference silently converts the distribution from publishable to
        checkout-only, and that should not happen by accident.
        """
        found = {_requirement_name(r) for r in metadata["project"]["dependencies"] if "@" in r}
        assert found == KNOWN_DIRECT_REFERENCES, (
            f"direct references changed: {sorted(found)} against the expected "
            f"{sorted(KNOWN_DIRECT_REFERENCES)}. Adding one blocks publishing to PyPI; "
            f"removing one is good news and should update this set."
        )

    def test_build_config_permits_those_direct_references(self, metadata: dict) -> None:
        """Without this hatchling flag the direct reference makes every build fail.

        Not a style check. `pip install -e .` dies in "Preparing editable metadata" with
        "Dependency #6 of field `project.dependencies` cannot be a direct reference unless
        field `tool.hatch.metadata.allow-direct-references` is set to `true`", and so does
        every wheel build, including the ones the packaging tests do.
        """
        if not KNOWN_DIRECT_REFERENCES:
            pytest.skip("no direct references to permit")
        allowed = metadata.get("tool", {}).get("hatch", {}).get("metadata", {})
        assert allowed.get("allow-direct-references") is True, (
            "project.dependencies contains a direct reference but "
            "tool.hatch.metadata.allow-direct-references is not true, so no build can succeed"
        )

    def test_runtime_dependencies_use_lower_bounds_not_pins(self, metadata: dict) -> None:
        """Users must be able to co-install DODO with other packages.

        An exact pin would make the package hard to use in any environment but its own, so the
        floors are ranges. The low end of those ranges is verified by the dedicated ``floors``
        job in CI, which installs the declared minimums rather than whatever pip resolves.
        """
        for requirement in metadata["project"]["dependencies"]:
            if _requirement_name(requirement) in KNOWN_DIRECT_REFERENCES:
                continue
            assert ">=" in requirement, f"{requirement!r} should declare a lower bound"
            assert "==" not in requirement, f"{requirement!r} pins an exact version"

    def test_a_ci_job_actually_installs_those_minimums(self, metadata: dict) -> None:
        """A floor nobody installs is a claim, not a guarantee."""
        text = CI_WORKFLOW.read_text()
        assert "floors:" in text, "no CI job tests the declared dependency floors"
        for requirement in metadata["project"]["dependencies"]:
            name = _requirement_name(requirement)
            if name in KNOWN_DIRECT_REFERENCES:
                continue
            floor = requirement.split(">=", 1)[1].strip()
            major_minor = ".".join(floor.split(".")[:2])
            assert f'"{name}=={major_minor}.*"' in text, (
                f"the floors job does not install {name}=={major_minor}.*, so the declared "
                f"floor {requirement!r} is never exercised"
            )
