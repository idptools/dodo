"""Sphinx configuration.

Replaces the 2018 Computational Molecular Science cookiecutter config, which declared
``project = 'dodo'`` with an empty version and pointed autosummary at a ``dodo.canvas`` module
that does not exist in 2.0.

Two choices worth explaining, because both are load-bearing:

**The API reference is generated, not written.** ``autosummary_generate`` with ``:recursive:``
walks the real package, so a new module or a new public name appears in the docs without anyone
remembering to add a page. A hand-maintained file per module drifts within days: both
``assign_regions_from_spec`` and ``reposition_folded_domains`` became public during a single
afternoon's work on this branch.

**sparrow and STARLING are mocked.** Neither is a hard dependency -- sparrow is not on PyPI and
STARLING is ~2.4 GB of weights -- and a docs build must not need either. DODO imports them lazily
and inside functions, so autodoc only ever meets them in annotations.
"""

from __future__ import annotations

import importlib.metadata

project = "DODO"
copyright = "2023-2026, Ryan Emenecker - Holehouse Lab"  # noqa: A001
author = "Ryan Emenecker - Holehouse Lab"

# Read the installed version rather than hardcoding one. hatch-vcs derives it from git tags, so
# any literal here would be wrong the moment a tag moved.
try:
    release = importlib.metadata.version("idptools-dodo")
except importlib.metadata.PackageNotFoundError:  # pragma: no cover - built from a bare checkout
    release = "0.0.0.dev0"
version = ".".join(release.split(".")[:2])

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store", "README.md", "_static/README.md"]

# -- Autodoc -----------------------------------------------------------------

autosummary_generate = True
autosummary_imported_members = False

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}
# Types are noise in the signature and useful in the parameter list. This also avoids the
# duplicated type rendering that sphinx-autodoc-typehints produces alongside napoleon, which is
# why that extension is deliberately not enabled here.
autodoc_typehints = "description"
autodoc_typehints_description_target = "documented_params"
autodoc_preserve_defaults = True
autodoc_mock_imports = ["sparrow", "starling"]

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
# Render an "Attributes" section as :ivar: fields rather than as separate object descriptions.
# Without this, every dataclass field is documented twice on the same page -- once from the
# class docstring's Attributes section and once by autodoc as a class member -- which produced
# 159 "duplicate object description" warnings and would block fail_on_warning.
napoleon_use_ivar = True

# -- Cross-references --------------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}

# -- MyST --------------------------------------------------------------------

myst_enable_extensions = ["colon_fence", "deflist"]
myst_heading_anchors = 3

# -- HTML --------------------------------------------------------------------

html_theme = "furo"
html_title = f"DODO {version}"
html_static_path = ["_static"]
