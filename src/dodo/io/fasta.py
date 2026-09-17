"""Reading FASTA, and working out which chain each record belongs to.

DODO needs a reference sequence when the input structure is missing residues -- a crystal or
cryo-EM structure models what it could see, and what it could not see is very often exactly the
disordered region DODO exists to rebuild. The reference says what should be there.

The hard part is not the format, it is the mapping. A FASTA downloaded from the RCSB names its
chains like this::

    >7R5J_11|Chains TA[auth I0], UA[auth I1], VA[auth I2], WA[auth I3]|Nucleoporin p58/p45|...

so one record covers four chains, each under two different names: the mmCIF ``label_asym_id``
and the author's own ``auth_asym_id``. DODO keys chains on the **author** id, because that is
what the reader uses and what a user sees in a viewer, so ``I0`` is the name that has to match --
but a hand-written FASTA is just as likely to say ``>I0`` or ``>chain I0``, and a biological
assembly renames every chain again (``I0-2``, ``I0-3``, ... for the symmetry copies).

Rather than pick one convention and make everyone else's file an error, :func:`match_chains`
tries the cheap exact answers first and falls back to matching on the sequence itself, which is
the one thing that cannot be renamed.
"""

from __future__ import annotations

import gzip
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from ..exceptions import StructureFileError

__all__ = [
    "FastaRecord",
    "match_chains",
    "parse_fasta",
    "read_fasta",
]

#: ``LABEL[auth AUTH]`` or a bare ``LABEL`` inside an RCSB "Chains ..." header field.
_CHAIN_TOKEN = re.compile(r"\s*([^\s,\[\]]+)\s*(?:\[\s*auth\s+([^\]]+?)\s*\])?\s*$")

#: Trailing ``-2``, ``-3``, ... that a biological-assembly file appends to duplicate a chain.
_ASSEMBLY_COPY_SUFFIX = re.compile(r"^(.*?)-\d+$")

#: Header prefixes that carry no chain information, so the first token is not a chain id.
_DATABASE_PREFIXES = frozenset({"sp", "tr", "gi", "ref", "pdb", "gb", "emb", "dbj"})


@dataclass(frozen=True, slots=True)
class FastaRecord:
    """One FASTA entry, with every chain id its header names."""

    header: str
    sequence: str
    #: Chain ids this record claims, most specific first: author ids before label ids, because
    #: the author id is what DODO keys chains on.
    chain_ids: tuple[str, ...]
    description: str = ""

    def __str__(self) -> str:
        chains = ", ".join(self.chain_ids) if self.chain_ids else "no chain named"
        return f"{len(self.sequence)} residues [{chains}] {self.description}".rstrip()


def _strip_assembly_copy(chain_id: str) -> str:
    """``"I0-3"`` -> ``"I0"``. A biological assembly numbers its symmetry copies this way."""
    match = _ASSEMBLY_COPY_SUFFIX.match(chain_id)
    return match.group(1) if match else chain_id


def _chain_ids_from_header(header: str) -> tuple[tuple[str, ...], str]:
    """Extract chain ids and a description from a FASTA header (without the leading ``>``)."""
    fields = [f.strip() for f in header.split("|")]
    auth_ids: list[str] = []
    label_ids: list[str] = []
    description = ""

    for index, field in enumerate(fields):
        lowered = field.lower()
        if not (lowered.startswith("chains ") or lowered.startswith("chain ")):
            continue
        listing = field.split(None, 1)[1] if " " in field else ""
        for token in listing.split(","):
            match = _CHAIN_TOKEN.match(token)
            if match is None:
                continue
            label, auth = match.group(1), match.group(2)
            if auth:
                auth_ids.append(auth)
            if label:
                label_ids.append(label)
        # The RCSB puts the molecule name in the field after the chain list.
        if index + 1 < len(fields):
            description = fields[index + 1]
        break
    else:
        # No "Chains" field. Either a bare ``>A`` / ``>A some description``, or a database
        # header like ``>sp|P04637|P53_HUMAN``, which names no chain at all.
        first = fields[0]
        token = first.split()[0] if first.split() else ""
        description = first[len(token) :].strip()
        if token and token.lower() not in _DATABASE_PREFIXES and len(fields) == 1:
            auth_ids.append(token)
        elif len(fields) > 1:
            description = " ".join(f for f in fields[1:] if f)

    # Author ids first, then any label id not already claimed: match order is preference order.
    ordered = list(dict.fromkeys([*auth_ids, *label_ids]))
    return tuple(ordered), description


def parse_fasta(text: str, *, source: str | None = None) -> list[FastaRecord]:
    """Parse FASTA text into records.

    Raises
    ------
    StructureFileError
        If the text contains no ``>`` header, or a header with no sequence under it. Both mean
        the file is not what the caller thinks it is, and guessing would produce a rebuild
        against the wrong sequence -- a failure that looks like a modelling problem rather than
        an input problem.
    """
    where = f" in {source}" if source else ""
    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        if header is None:
            return
        sequence = "".join(chunks).replace("*", "").replace("-", "").upper()
        if not sequence:
            raise StructureFileError(f"FASTA record {header!r}{where} has no sequence.")
        chain_ids, description = _chain_ids_from_header(header)
        records.append(
            FastaRecord(
                header=header, sequence=sequence, chain_ids=chain_ids, description=description
            )
        )

    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            chunks = []
        elif header is not None:
            chunks.append(re.sub(r"\s+", "", line))
    flush()

    if not records:
        raise StructureFileError(
            f"No FASTA records found{where}. A FASTA file has at least one line starting "
            f"with '>'."
        )
    return records


def read_fasta(path: str | Path) -> list[FastaRecord]:
    """Read a FASTA file, ``.gz`` accepted."""
    path = Path(path)
    try:
        if path.suffix == ".gz":
            with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
                text = handle.read()
        else:
            text = path.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise StructureFileError(f"Could not read {path}: {exc}") from exc
    return parse_fasta(text, source=str(path))


def _is_subsequence(observed: str, full: str) -> bool:
    """Report whether every residue of ``observed`` appears in ``full`` in order.

    The relationship between a modelled chain and its reference: the structure shows a subset of
    the residues, in order, with gaps where nothing was resolved. Cheap -- one pass -- and it is
    what makes sequence matching usable as a fallback on a file with 808 renamed chains.
    """
    if len(observed) > len(full):
        return False
    position = 0
    for residue in observed:
        position = full.find(residue, position)
        if position < 0:
            return False
        position += 1
    return True


def match_chains(
    records: Sequence[FastaRecord],
    chain_ids: Iterable[str],
    observed: Mapping[str, str] | None = None,
) -> tuple[dict[str, str], list[str]]:
    """Work out which FASTA record belongs to each chain.

    Parameters
    ----------
    records
        Parsed FASTA records.
    chain_ids
        The structure's chain ids, as DODO reads them (author ids).
    observed
        Each chain's observed one-letter sequence. Optional, but without it the only thing
        available is the names in the headers, and a biological assembly renames every chain.

    Returns
    -------
    tuple
        ``(sequences, notes)``: chain id to reference sequence for every chain that matched, and
        human-readable notes describing how the awkward ones were resolved and which were left
        unmatched.

    Notes
    -----
    Four passes, cheapest and most certain first:

    1. the chain id appears verbatim in a header;
    2. it appears after stripping a biological assembly's ``-2``/``-3`` copy suffix;
    3. the chain's observed sequence is a subsequence of exactly one record's sequence, or of
       several records that all carry the *same* sequence -- identical entities are
       interchangeable by definition;
    4. nothing matched, which is reported rather than guessed at.

    Pass 3 is the one that earns its keep. Matching on sequence cannot be defeated by renaming,
    and renaming is exactly what a biological-assembly file does to every chain it duplicates.
    """
    chain_ids = list(chain_ids)
    observed = dict(observed or {})
    by_id: dict[str, FastaRecord] = {}
    for record in records:
        for name in record.chain_ids:
            previous = by_id.get(name)
            if previous is not None and previous.sequence != record.sequence:
                raise StructureFileError(
                    f"FASTA chain id {name!r} is claimed by two records with different "
                    f"sequences ({previous.header!r} and {record.header!r}); DODO cannot choose "
                    "which sequence to insert. Give each record an unambiguous chain id."
                )
            by_id.setdefault(name, record)

    sequences: dict[str, str] = {}
    notes: list[str] = []
    renamed: list[str] = []
    by_sequence: list[str] = []
    unmatched: list[str] = []
    ambiguous: list[str] = []

    for chain_id in chain_ids:
        named = by_id.get(chain_id)
        if named is not None:
            sequences[chain_id] = named.sequence
            continue

        stripped = _strip_assembly_copy(chain_id)
        if stripped != chain_id and stripped in by_id:
            sequences[chain_id] = by_id[stripped].sequence
            renamed.append(chain_id)
            continue

        chain_sequence = observed.get(chain_id)
        if not chain_sequence:
            unmatched.append(chain_id)
            continue
        candidates = [r for r in records if _is_subsequence(chain_sequence, r.sequence)]
        distinct = {r.sequence for r in candidates}
        if len(distinct) == 1:
            sequences[chain_id] = candidates[0].sequence
            by_sequence.append(chain_id)
        elif not candidates:
            unmatched.append(chain_id)
        else:
            ambiguous.append(chain_id)

    if renamed:
        notes.append(
            f"{len(renamed)} chain(s) matched a FASTA record after stripping a biological "
            f"assembly copy suffix (e.g. {renamed[0]!r})."
        )
    if by_sequence:
        notes.append(
            f"{len(by_sequence)} chain(s) matched a FASTA record by sequence rather than by "
            f"name (e.g. {by_sequence[0]!r}); the headers did not name them."
        )
    if ambiguous:
        notes.append(
            f"{len(ambiguous)} chain(s) are a subsequence of two or more DIFFERENT FASTA "
            f"records and were left unmatched rather than guessed at: "
            f"{', '.join(ambiguous[:8])}"
            + (f" and {len(ambiguous) - 8} more" if len(ambiguous) > 8 else "")
            + ". Name them in the headers to disambiguate."
        )
    if unmatched:
        notes.append(
            f"{len(unmatched)} chain(s) had no matching FASTA record and were left as-is: "
            f"{', '.join(unmatched[:8])}"
            + (f" and {len(unmatched) - 8} more" if len(unmatched) > 8 else "")
            + "."
        )
    return sequences, notes
