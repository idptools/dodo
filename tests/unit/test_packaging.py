"""Tests that the distribution actually contains the package.

These exist because DODO shipped a broken distribution for months without anyone
noticing. Two independent causes, both silent, both leaving ``pip install .`` exiting 0:

1. ``dodo/backend/`` had no ``__init__.py``, and setuptools' ``packages.find`` with
   ``namespaces = false`` therefore skipped it -- excluding 3,410 of the package's
   3,848 lines from every wheel and sdist. The published wheel was 12.7 KB.
2. During the rewrite, a subpackage named ``build`` was silently excluded because the
   repository's ``.gitignore`` contains ``build/`` (the standard pattern for the build
   artifact directory) and hatchling honours gitignore. Same outcome, different cause.

The lesson is that "the package imports in my checkout" proves nothing about what
users receive, so the invariant has to be asserted against a built artifact.
"""

from __future__ import annotations

import subprocess
import sys
import zipfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src" / "dodo"


def _source_subpackages() -> set[str]:
    """Every subpackage present in the source tree, as dotted names."""
    return {
        f"dodo.{path.parent.relative_to(SRC_ROOT).as_posix().replace('/', '.')}"
        for path in SRC_ROOT.rglob("__init__.py")
        if path.parent != SRC_ROOT
    }


def test_source_tree_has_subpackages() -> None:
    """Sanity check on the test itself: the source tree must have subpackages to find.

    Without this, a refactor that flattened the layout would make the wheel test
    vacuously pass.
    """
    assert _source_subpackages(), "found no subpackages under src/dodo/"


def test_every_source_module_is_importable() -> None:
    """Every module in the source tree imports without error.

    Catches the v1 failure directly: ``dodo/build.py`` imported six modules that had
    been deleted, so ``import dodo`` raised, and nothing in CI noticed because CI was
    disabled.
    """
    import importlib

    failures: list[str] = []
    for path in sorted(SRC_ROOT.rglob("*.py")):
        rel = path.relative_to(SRC_ROOT)
        parts = list(rel.parts)
        if parts[-1] == "__init__.py":
            parts.pop()
        else:
            parts[-1] = parts[-1].removesuffix(".py")
        name = ".".join(["dodo", *parts]) if parts else "dodo"
        try:
            importlib.import_module(name)
        except Exception as exc:
            failures.append(f"{name}: {type(exc).__name__}: {exc}")

    assert not failures, "modules failed to import:\n  " + "\n  ".join(failures)


@pytest.mark.slow
def test_wheel_contains_every_subpackage(tmp_path: Path) -> None:
    """A built wheel contains every subpackage found in the source tree.

    This is the test whose absence let a 12.7 KB engine-less wheel ship.
    """
    # A real guard, unlike the `except FileNotFoundError` this replaces. That could never fire:
    # the subprocess is `sys.executable -m build`, and sys.executable always exists, so a missing
    # backend exits 1 and raises CalledProcessError -- a FAILURE, not the intended skip. It went
    # unnoticed because developer machines have `build`; CI installs only `.[test]`, which did
    # not, so every matrix leg would have gone red. `build` is now pinned in that extra too.
    pytest.importorskip("build")
    try:
        subprocess.run(
            [sys.executable, "-m", "build", "--wheel", "--outdir", str(tmp_path)],
            cwd=REPO_ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:  # pragma: no cover
        pytest.fail(f"wheel build failed:\n{exc.stdout}\n{exc.stderr}")

    wheels = list(tmp_path.glob("*.whl"))
    assert len(wheels) == 1, f"expected exactly one wheel, got {wheels}"

    with zipfile.ZipFile(wheels[0]) as archive:
        names = archive.namelist()

    packaged = {
        name.rsplit("/", 1)[0].replace("/", ".") for name in names if name.endswith("__init__.py")
    }
    missing = _source_subpackages() - packaged
    assert not missing, (
        f"these subpackages exist in src/ but are absent from the wheel: "
        f"{sorted(missing)}. A likely cause is a .gitignore pattern matching the "
        f"directory name -- hatchling honours gitignore, so a subpackage named e.g. "
        f"'build' or 'dist' disappears silently."
    )

    # The type marker has to ship or downstream type checking silently does nothing.
    assert "dodo/py.typed" in names, "py.typed is missing from the wheel"


@pytest.mark.slow
def test_wheel_installs_and_imports(tmp_path: Path) -> None:
    """A built wheel installs into a clean environment and imports.

    The end-to-end invariant that no amount of source-tree testing can establish: this
    is what a user actually gets.
    """
    pytest.importorskip("build")
    subprocess.run(
        [sys.executable, "-m", "build", "--wheel", "--outdir", str(tmp_path / "dist")],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )

    venv = tmp_path / "venv"
    subprocess.run([sys.executable, "-m", "venv", str(venv)], check=True, capture_output=True)
    python = venv / ("Scripts" if sys.platform == "win32" else "bin") / "python"
    wheel = next((tmp_path / "dist").glob("*.whl"))

    subprocess.run(
        [str(python), "-m", "pip", "install", "--quiet", str(wheel)],
        check=True,
        capture_output=True,
        text=True,
    )
    result = subprocess.run(
        [
            str(python),
            "-c",
            "import dodo; from dodo.structure import Structure; print(dodo.__version__)",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.stdout.strip(), "importing the installed package produced no version"


def test_version_is_pep440_clean() -> None:
    """The reported version carries no local segment.

    versioningit previously produced versions like ``0.0.0+36.g45f3738``. PyPI rejects
    any distribution whose version has a local segment, so this would have blocked
    publication at upload time -- and ``twine check`` does not catch it.
    """
    packaging_version = pytest.importorskip("packaging.version")

    import dodo

    version = packaging_version.Version(dodo.__version__)
    assert version.local is None, (
        f"version {dodo.__version__!r} has local segment {version.local!r}; PyPI will reject it"
    )
