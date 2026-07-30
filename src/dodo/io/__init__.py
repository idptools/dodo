"""Structure file reading and writing, and remote structure retrieval.

Format dispatch
---------------
:func:`read_structure` is the entry point most callers want: it picks the reader from the
file's extension and content, so nothing upstream has to care whether a structure arrived
as PDB or mmCIF.

Why the parsers are hand-rolled
-------------------------------
DODO parses PDB and mmCIF itself rather than depending on biotite, gemmi or MDAnalysis.
That is a deliberate choice to keep the runtime dependency set at numpy plus scipy, which
is what makes the base install light enough for a visualization tool. The cost is that
every format edge case is ours to handle, so the readers are written defensively and every
inference or omission is recorded in :attr:`~dodo.structure.Structure.notes` rather than
applied silently.
"""

from __future__ import annotations

import gzip
from pathlib import Path

from ..exceptions import StructureFileError, UnsupportedFormatError
from ..structure import Structure
from .cif import parse_cif_text, read_cif, read_unobserved_residues
from .fetch import (
    AlphaFoldModel,
    default_cache_dir,
    fetch_alphafold,
    fetch_alphafold_model,
    fetch_pdb,
    fetch_uniprot_sequence,
    resolve_uniprot_accession,
)
from .pdb import decode_hybrid36, parse_pdb_lines, read_pdb
from .write import structure_to_pdb_lines, write_pdb

__all__ = [
    "AlphaFoldModel",
    "decode_hybrid36",
    "default_cache_dir",
    "fetch_alphafold",
    "fetch_alphafold_model",
    "fetch_pdb",
    "fetch_uniprot_sequence",
    "parse_cif_text",
    "parse_pdb_lines",
    "read_cif",
    "read_pdb",
    "read_structure",
    "read_unobserved_residues",
    "resolve_uniprot_accession",
    "structure_to_pdb_lines",
    "write_pdb",
]

#: Extensions handled by the PDB reader, lowercased, with any ``.gz`` already stripped.
_PDB_SUFFIXES = frozenset({".pdb", ".ent"})
#: Extensions handled by the mmCIF reader.
_CIF_SUFFIXES = frozenset({".cif", ".mmcif"})
#: How many lines to inspect when sniffing a file's format.
_SNIFF_LINES = 200


def _content_looks_like_cif(path: Path) -> bool:
    """Return True if the file appears to be mmCIF, judging by its content.

    Extension dispatch alone is not enough: files arrive named ``.pdb`` that are really
    mmCIF often enough to be worth one cheap read.

    Keyed on STAR constructs (``data_``, ``loop_``, ``_atom_site.``) rather than on what
    the atom records look like, because a column-aligned mmCIF ``_atom_site`` row begins
    with the literal text ``ATOM  `` and so is indistinguishable from a PDB ATOM record by
    prefix alone. Getting that backwards is how a legal mmCIF ends up reported as "this is
    a PDB file".
    """
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    try:
        with opener(path, "rt", encoding="latin-1", errors="replace") as handle:
            for _ in range(_SNIFF_LINES):
                line = handle.readline()
                if not line:
                    break
                stripped = line.lstrip()
                if stripped.startswith(("data_", "loop_", "_atom_site.")):
                    return True
                if stripped.startswith(("ATOM  ", "HETATM", "MODEL ", "CRYST1", "SEQRES")):
                    # Reached PDB coordinate or primary-structure records without seeing a
                    # STAR construct, so this is PDB.
                    return False
    except OSError as exc:
        raise StructureFileError(f"Could not read {path}: {exc}") from exc
    return False


def read_structure(
    path: str | Path,
    *,
    model: int | None = None,
    keep_hydrogens: bool = False,
) -> Structure:
    """Read a structure file, choosing the reader from its extension and content.

    Parameters
    ----------
    path
        Path to a PDB or mmCIF file, optionally gzipped.
    model
        Which model to read from a multi-model file, 1-based. Defaults to the first.
        Reading only one is deliberate: merging frames yields several atoms per position,
        which is geometrically meaningless.
    keep_hydrogens
        Keep hydrogens. DODO's geometry is heavy-atom based, so they are dropped by default.

    Returns
    -------
    Structure
        The parsed structure. Check :attr:`~dodo.structure.Structure.notes` for anything
        the reader inferred, skipped or deduplicated.

    Raises
    ------
    UnsupportedFormatError
        If the file is neither PDB nor mmCIF.
    StructureFileError
        If the file cannot be read.
    """
    path = Path(path)
    suffixes = [s.lower() for s in path.suffixes]
    # Peel a trailing .gz so "x.pdb.gz" dispatches on ".pdb".
    if suffixes and suffixes[-1] == ".gz":
        suffixes = suffixes[:-1]
    suffix = suffixes[-1] if suffixes else ""

    # A numbered biological assembly is ".pdb1", ".pdb2", and so on.
    if suffix.startswith(".pdb") and suffix[4:].isdigit():
        suffix = ".pdb"

    if suffix in _CIF_SUFFIXES:
        return read_cif(path, model=model, keep_hydrogens=keep_hydrogens)

    if suffix in _PDB_SUFFIXES:
        # Trust content over extension: mmCIF misnamed .pdb is common enough that
        # dispatching on the name alone produces a baffling parse error.
        if _content_looks_like_cif(path):
            return read_cif(path, model=model, keep_hydrogens=keep_hydrogens)
        return read_pdb(path, model=model, keep_hydrogens=keep_hydrogens)

    # Unknown extension: sniff rather than refuse, since structures are routinely
    # downloaded with no extension at all.
    if _content_looks_like_cif(path):
        return read_cif(path, model=model, keep_hydrogens=keep_hydrogens)
    try:
        return read_pdb(path, model=model, keep_hydrogens=keep_hydrogens)
    except StructureFileError as exc:
        raise UnsupportedFormatError(
            f"{path} is neither PDB nor mmCIF, as far as DODO can tell. Reading it as PDB "
            f"failed with: {exc}"
        ) from exc
