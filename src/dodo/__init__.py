"""DODO: rebuild the disordered regions of predicted protein structures.

DODO takes a predicted structure -- typically from AlphaFold -- identifies its folded
domains, intrinsically disordered regions and loops, and rebuilds the disordered parts
so they adopt realistic polymer dimensions instead of AlphaFold's characteristic
extended spaghetti. The result is a structure that looks like what the protein actually
is, which matters for figures and presentations.

Quick start
-----------
>>> import dodo
>>> models = dodo.rebuild("AF-P04637-F1-model_v6.pdb")   # doctest: +SKIP
>>> dodo.write_pdb(models, "p53_dodo.pdb")               # doctest: +SKIP

Or from a sequence alone:

>>> models = dodo.build_from_sequence("MEEPQSDPSV...")   # doctest: +SKIP

Importing is cheap
------------------
This module resolves its submodules lazily (PEP 562 / Scientific Python SPEC 1). That
is not cosmetic: the optional prediction dependencies pull in torch, and eagerly
importing them made ``dodo --help`` pay several seconds of framework startup before
printing anything. Nothing heavy is imported until you touch an attribute that needs it.
"""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING, Any

try:  # pragma: no cover - trivial, and depends on install state
    from ._version import __version__
except ImportError:  # pragma: no cover
    # Running from a source tree with no build metadata. Report something honest rather
    # than crashing on import, which is what the v1 layout did. Deliberately carries no
    # PEP 440 local version segment ("+unknown" and friends), because a local segment
    # makes a distribution unpublishable and this string can end up in one.
    __version__ = "0.0.0.dev0"

# Names resolved lazily on first access, mapped to the submodule that defines them.
_LAZY_ATTRS: dict[str, str] = {
    # Core data model
    "Structure": "dodo.structure",
    "Chain": "dodo.structure",
    "Domain": "dodo.structure",
    "DomainKind": "dodo.structure",
    "Span": "dodo.structure",
    # Exceptions
    "DodoError": "dodo.exceptions",
    "BuildError": "dodo.exceptions",
    "GeometryError": "dodo.exceptions",
    "RegionIdentificationError": "dodo.exceptions",
    "StructureFileError": "dodo.exceptions",
    "MissingDependencyError": "dodo.exceptions",
}

_LAZY_MODULES: tuple[str, ...] = (
    "constants",
    "exceptions",
    "structure",
)

if TYPE_CHECKING:
    # Give type checkers and IDEs the real thing; this branch never runs.
    from . import constants, exceptions, structure
    from .exceptions import (
        BuildError,
        DodoError,
        GeometryError,
        MissingDependencyError,
        RegionIdentificationError,
        StructureFileError,
    )
    from .structure import Chain, Domain, DomainKind, Span, Structure


def __getattr__(name: str) -> Any:
    """Resolve submodules and re-exported names on first access (PEP 562)."""
    if name in _LAZY_MODULES:
        module = importlib.import_module(f".{name}", __name__)
        globals()[name] = module
        return module
    if name in _LAZY_ATTRS:
        module = importlib.import_module(_LAZY_ATTRS[name])
        value = getattr(module, name)
        globals()[name] = value
        return value
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    """Include lazily-resolved names in ``dir()`` and tab completion."""
    return sorted({*globals(), *_LAZY_MODULES, *_LAZY_ATTRS})


__all__ = [
    "__version__",
    *_LAZY_MODULES,
    *_LAZY_ATTRS,
]
