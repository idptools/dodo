"""Remote retrieval of AlphaFold models, PDB entries and UniProt sequences.

Why this module exists in the form it does
-----------------------------------------
v1 built AlphaFold download URLs by string formatting::

    f'https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb'

EBI has since published ``model_v6`` and retired the older files: as of this writing
``model_v4`` and ``model_v5`` both return HTTP 404 while ``model_v6`` returns 200, so
*every* fetch-by-name call in the released package fails. A hardcoded version number is
a time bomb, and this is it going off.

Nothing here ever constructs a file URL. The current URL is read from the prediction
API::

    https://alphafold.ebi.ac.uk/api/prediction/{ACCESSION}  ->  [...]['pdbUrl']

which reports whatever the live release is and will keep working across future version
bumps. The same response carries the model's mean pLDDT and the fraction of residues in
each confidence band, which is directly useful for region identification, so it is
returned to the caller in :class:`AlphaFoldModel` rather than thrown away.

Three further v1 behaviours are deliberately not reproduced:

* No request had a timeout, so a hung connection hung DODO forever. Every request here
  takes one, and it is a required argument of the private helpers so it cannot be
  forgotten. The one call DODO does not make itself -- the optional ``getSequence``
  name lookup, which uses ``requests`` internally and accepts no timeout -- is bounded
  from the outside instead, so no path through this module can wait indefinitely.
* ``[0]`` of the API response was assumed to be the model of the requested accession.
  It is not: for UBR4 (Q5T4S7, 5183 residues) the canonical sequence has no model at
  all and the response contains only alternative-isoform entries, the first of which is
  212 residues long. Silently modelling a 212-residue isoform as if it were the protein
  the user asked for is exactly the class of plausible-but-wrong output this rewrite
  exists to eliminate, so entries are matched on ``uniprotAccession``.
* ``-F1`` was assumed to be the whole protein, and when its sequence did not match
  UniProt the code raised ``'AF2 sequence does not match sequence retrieved from
  Uniprot!'``. For a long protein split into fragments that message is actively
  misleading -- the mismatch is expected. Fragments are detected and described.

Requests use :mod:`urllib` from the standard library. ``requests`` is not a declared
dependency of this package and does not become one for this.
"""

from __future__ import annotations

import http.client
import importlib
import io
import json
import os
import re
import shutil
import sys
import threading
import urllib.error
import urllib.parse
import urllib.request
import warnings
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from ..exceptions import (
    MissingDependencyError,
    RetrievalError,
    StructureNotFoundError,
    UnsupportedFormatError,
)

__all__ = [
    "AlphaFoldModel",
    "default_cache_dir",
    "fetch_alphafold",
    "fetch_alphafold_model",
    "fetch_pdb",
    "fetch_uniprot_sequence",
    "resolve_uniprot_accession",
]

# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------
#
# These live here rather than in constants.py because they are facts about external
# web services, not physical or algorithmic constants of the method.

#: AlphaFold DB prediction metadata. The *only* supported way to learn a file URL.
AFDB_PREDICTION_API: Final[str] = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

#: RCSB coordinate file download.
RCSB_DOWNLOAD_URL: Final[str] = "https://files.rcsb.org/download/{entry}.{extension}"

#: UniProtKB single-entry FASTA.
UNIPROT_FASTA_URL: Final[str] = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"

#: UniProtKB free-text search, used to resolve a protein name to an accession.
UNIPROT_SEARCH_URL: Final[str] = "https://rest.uniprot.org/uniprotkb/search"

#: Hosts DODO will talk to. Checked on every request, including URLs that came *out of*
#: a remote response: ``pdbUrl`` is remote data, and remote data does not get to point
#: DODO at ``file:///`` or at an arbitrary third-party host.
_ALLOWED_HOSTS: Final[frozenset[str]] = frozenset(
    {"alphafold.ebi.ac.uk", "files.rcsb.org", "rest.uniprot.org"}
)

#: Sequence length above which AlphaFold DB splits a prediction into ``-F1``, ``-F2``,
#: ... fragments, or omits it entirely. Published by EBI alongside the database.
AFDB_FRAGMENT_LENGTH_LIMIT: Final[int] = 2700

#: Rows requested from the UniProt search. Only the first is used; see
#: :func:`_search_uniprot` for why this cannot be 1.
_UNIPROT_SEARCH_PAGE_SIZE: Final[int] = 5

_AFDB_CACHE_SUBDIR: Final[str] = "alphafold"
_PDB_CACHE_SUBDIR: Final[str] = "pdb"

#: Official UniProtKB accession syntax, plus an optional ``-N`` isoform suffix (the
#: AFDB API accepts isoform accessions, e.g. ``P04637-2``, and returns their models).
_UNIPROT_ACCESSION_RE: Final[re.Pattern[str]] = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?$"
)

#: Classic 4-character PDB IDs and the extended ``pdb_XXXXXXXX`` form.
_PDB_ID_RE: Final[re.Pattern[str]] = re.compile(r"^(?:[0-9][A-Z0-9]{3}|PDB_[A-Z0-9]{8})$")

#: ``-F<n>`` suffix of an AFDB entry id.
_FRAGMENT_SUFFIX_RE: Final[re.Pattern[str]] = re.compile(r"-F([0-9]+)$")

_FORMAT_ALIASES: Final[dict[str, str]] = {
    "cif": "cif",
    "mmcif": "cif",
    "pdb": "pdb",
}


def _user_agent() -> str:
    """Return the User-Agent string DODO identifies itself with.

    EBI, RCSB and UniProt all ask automated clients to identify themselves, and a
    named agent is what lets them tell us apart from a scraper if we ever misbehave.
    """
    from .. import __version__

    return f"dodo/{__version__} (https://github.com/idptools/dodo)"


# ---------------------------------------------------------------------------
# Cache location
# ---------------------------------------------------------------------------


#: Cache layout generation, appended to the cache directory.
#:
#: Bumping this abandons every previously cached file rather than migrating it. That is
#: the point: generation 1 was written before downloads were length-verified, so a cache
#: entry from it may be a silently truncated fragment, and nothing recorded its expected
#: size, so there is no way to tell a good entry from a bad one after the fact. The
#: cache-hit gate only checks that a file is non-empty, which a fragment satisfies.
#: Abandoning the generation is cheap (the files re-download on demand) and is the only
#: honest way to guarantee a poisoned entry is never served again.
#:
#: Bump this whenever a change could make previously cached content untrustworthy.
CACHE_GENERATION = "v2"


def default_cache_dir() -> Path:
    r"""Return the per-user cache directory DODO downloads into.

    Follows the same platform conventions as ``platformdirs`` but computed with the
    standard library only, since DODO declares no dependency on it:

    * macOS: ``~/Library/Caches/dodo/<generation>``
    * Windows: ``%LOCALAPPDATA%\dodo\Cache\<generation>``
    * everything else: ``$XDG_CACHE_HOME/dodo/<generation>``, falling back to
      ``~/.cache/dodo/<generation>``

    ``DODO_CACHE_DIR`` overrides the platform base, which is what makes a CI run
    hermetic. The generation segment is appended to the override too, so an override
    cannot accidentally resurrect an untrustworthy earlier layout.

    Returns
    -------
    Path
        The cache directory. Not created; the fetchers create it on demand.
    """
    override = os.environ.get("DODO_CACHE_DIR")
    if override:
        return Path(override).expanduser() / CACHE_GENERATION
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Caches" / "dodo" / CACHE_GENERATION
    if sys.platform == "win32":
        local_app_data = os.environ.get("LOCALAPPDATA")
        base = Path(local_app_data) if local_app_data else Path.home() / "AppData" / "Local"
        return base / "dodo" / "Cache" / CACHE_GENERATION
    xdg = os.environ.get("XDG_CACHE_HOME")
    base_dir = Path(xdg).expanduser() if xdg else Path.home() / ".cache"
    return base_dir / "dodo" / CACHE_GENERATION


# ---------------------------------------------------------------------------
# HTTP plumbing
# ---------------------------------------------------------------------------


#: Transport failures that can come out of a *body read*, after ``urlopen`` has already
#: returned a response.
#:
#: ``OSError`` alone is not enough, and the omission is not academic:
#: ``http.client.IncompleteRead`` -- the error the standard library raises when a
#: connection closes before ``Content-Length`` bytes have arrived, which is the single
#: most likely way a download dies -- is an ``HTTPException``, and ``HTTPException`` is
#: not an ``OSError``. Catching only ``OSError`` let it escape as a raw stdlib exception
#: past the documented ``Raises: RetrievalError`` contract of every public fetcher,
#: taking the "this is transient, retry" advice with it.
_MID_BODY_ERRORS: Final[tuple[type[BaseException], ...]] = (
    OSError,  # covers URLError and TimeoutError, both subclasses
    http.client.HTTPException,
)


class _HttpStatusError(Exception):
    """Internal carrier for a non-2xx response.

    Not part of the public hierarchy: every call site translates it into a
    :class:`~dodo.exceptions.RetrievalError` or
    :class:`~dodo.exceptions.StructureNotFoundError`, because "404" means something
    different for each endpoint and the caller needs the endpoint-specific advice.
    """

    def __init__(self, status: int, url: str, body: bytes = b""):
        self.status = status
        self.url = url
        self.body = body
        super().__init__(f"HTTP {status} from {url}")


def _urlopen(request: urllib.request.Request, timeout: float) -> Any:
    """Open ``request``, the single point at which DODO touches the network.

    Every remote read in DODO funnels through here. That is deliberate: it is the seam
    the offline tests monkeypatch, so the whole retrieval layer is testable without a
    network, and it is the one place a timeout can be verified to be present.
    """
    return urllib.request.urlopen(request, timeout=timeout)


def _check_url(url: str) -> None:
    """Reject any URL that is not https to a known service host.

    Raises
    ------
    RetrievalError
        If the scheme is not https or the host is not one DODO talks to. This guards
        the URLs taken from AFDB's JSON response, which are remote data.
    """
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme != "https" or parsed.hostname not in _ALLOWED_HOSTS:
        raise RetrievalError(
            f"Refusing to fetch {url!r}: DODO only makes https requests to "
            f"{', '.join(sorted(_ALLOWED_HOSTS))}."
        )


def _open(url: str, *, timeout: float, accept: str | None = None) -> Any:
    """Open a URL, translating transport failures into DODO exceptions.

    Parameters
    ----------
    url
        Absolute https URL on a known service host.
    timeout
        Socket timeout in seconds. Required, and must be positive.
    accept
        Optional ``Accept`` header value.

    Returns
    -------
    Any
        An open file-like response, usable as a context manager.

    Raises
    ------
    RetrievalError
        On a transient transport failure (DNS, connection refused, timeout, a response
        the HTTP client cannot parse), or an invalid timeout.
    _HttpStatusError
        On a non-2xx status, for the caller to translate.
    """
    if timeout <= 0:
        raise RetrievalError(f"timeout must be positive, got {timeout!r}.")
    _check_url(url)
    headers = {"User-Agent": _user_agent()}
    if accept is not None:
        headers["Accept"] = accept
    request = urllib.request.Request(url, headers=headers, method="GET")
    try:
        return _urlopen(request, timeout)
    except urllib.error.HTTPError as exc:
        # Read the body now: HTTPError's file handle is closed once we leave this frame,
        # and the body is where UniProt and AFDB put their explanation of a 400.
        try:
            body = exc.read()
        except Exception:  # a diagnostic read must never mask the error itself
            body = b""
        raise _HttpStatusError(exc.code, url, body) from exc
    except (urllib.error.URLError, TimeoutError, http.client.HTTPException) as exc:
        raise RetrievalError(
            f"Could not reach {url}: {exc}. This is a transient network failure "
            f"(timeout was {timeout:g} s), so retrying, raising the timeout, or "
            f"checking your connection is the right response."
        ) from exc


def _get_bytes(url: str, *, timeout: float, accept: str | None = None) -> bytes:
    """Fetch a URL into memory. See :func:`_open` for the failure modes."""
    try:
        with _open(url, timeout=timeout, accept=accept) as response:
            return bytes(response.read())
    except _MID_BODY_ERRORS as exc:
        # A connection dropped mid-body raises out of read(), past _open's handler.
        # See _MID_BODY_ERRORS for why OSError alone is not the right net to cast.
        raise RetrievalError(
            f"Reading {url} failed part-way through: {type(exc).__name__}: {exc}. "
            f"This is transient; retry."
        ) from exc


def _get_json(url: str, *, timeout: float) -> Any:
    """Fetch a URL and parse it as JSON.

    Raises
    ------
    RetrievalError
        If the body is not valid JSON. A service returning an HTML error page with a
        200 status is a real occurrence and must not surface as a JSONDecodeError.
    """
    raw = _get_bytes(url, timeout=timeout, accept="application/json")
    try:
        return json.loads(raw)
    except json.JSONDecodeError as exc:
        preview = raw[:200].decode("utf-8", errors="replace")
        raise RetrievalError(f"{url} did not return JSON. First bytes: {preview!r}") from exc


def _content_length(response: Any) -> int | None:
    """Return the body length a response declared, or ``None`` if there is none to trust.

    ``None`` means "there is nothing to verify against", not "verified".

    Returns ``None`` for a chunked response even when a ``Content-Length`` header is
    present. Per RFC 9112 section 6.1, chunked framing takes precedence and
    ``Content-Length`` must be ignored; ``http.client`` implements exactly that
    (``if length and not self.chunked`` in ``HTTPResponse.begin``, which leaves
    ``response.length`` as ``None``). Reading the header regardless would reject a
    *complete* chunked download as truncated whenever an intermediary adds a
    ``Content-Length`` -- and every endpoint DODO downloads from serves chunked today, so
    a corporate proxy or CDN in the path would break all structure downloads with a
    misleading "bytes are missing" message and no workaround.

    A header a proxy has mangled into something non-numeric is no more usable than an
    absent one, so that also yields ``None``.
    """
    # Mirror http.client's own framing decision rather than second-guessing it.
    if getattr(response, "chunked", False):
        return None

    getter = getattr(getattr(response, "headers", None), "get", None)
    if getter is None:
        return None
    raw = getter("Content-Length")
    if raw is None:
        return None
    try:
        declared = int(str(raw).strip())
    except (TypeError, ValueError):
        return None
    return declared if declared >= 0 else None


def _download_to(url: str, destination: Path, *, timeout: float) -> None:
    """Stream ``url`` to ``destination``, atomically, and verify its length.

    Three protections, because each covers a failure the others do not:

    * The download lands in a sibling ``.part`` file and is renamed into place only once
      every check has passed, so neither a process killed mid-download nor a body
      rejected below can leave a file in the cache.
    * The bytes written are compared against the response's ``Content-Length``. The
      atomic rename does *not* protect against a truncated *response*, and a response
      sent with a fixed length does not raise when it is cut short:
      ``http.client.HTTPResponse.readinto`` deliberately tolerates a short body ("it
      might break compatibility"), so ``read()`` returns ``b""`` and
      ``shutil.copyfileobj`` finishes normally on a body that lost most of its content.
      Measured without this check: a 100-byte fragment of a 4.7 kB model passed the only
      gate there was (non-empty), was committed by ``os.replace``, and was then served
      as a cache hit forever -- one dropped connection turning into permanent silent
      corruption of every structure built from that file.
    * Bodies sent with ``Transfer-Encoding: chunked`` declare no length, and this is the
      common case rather than an exotic one: measured against the live services, both
      AFDB file downloads and RCSB downloads are chunked, and only UniProt sends a
      ``Content-Length``. Chunked framing is self-terminating, so a stream that stops
      early fails the framing check inside ``http.client`` and surfaces as
      ``IncompleteRead`` -- which is an ``HTTPException``, not an ``OSError``, and so is
      only translated into :class:`~dodo.exceptions.RetrievalError` here because
      :data:`_MID_BODY_ERRORS` catches both. Verified end to end against a local server
      that hangs up mid-stream in both framings.

    Raises
    ------
    RetrievalError
        On transport failure, if the response body is empty, or if the number of bytes
        received disagrees with the declared ``Content-Length``. In every case nothing is
        written to ``destination``.
    _HttpStatusError
        On a non-2xx status.
    """
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(f"{destination.name}.{os.getpid()}.part")
    declared: int | None = None
    try:
        try:
            with _open(url, timeout=timeout) as response, partial.open("wb") as handle:
                declared = _content_length(response)
                shutil.copyfileobj(response, handle)
        except _MID_BODY_ERRORS as exc:
            # A connection dropped mid-body raises out of read(), past _open's handler.
            # See _MID_BODY_ERRORS for why OSError alone is not the right net to cast.
            raise RetrievalError(
                f"Download of {url} failed part-way through: {type(exc).__name__}: "
                f"{exc}. Nothing was cached; this is transient, so retry."
            ) from exc
        written = partial.stat().st_size
        if written == 0:
            raise RetrievalError(
                f"{url} returned an empty file. Nothing was cached; retry, and if it "
                f"persists the entry may have been withdrawn."
            )
        if declared is not None and written < declared:
            raise RetrievalError(
                f"Download of {url} is incomplete: the server declared "
                f"Content-Length {declared} bytes but only {written} arrived, so "
                f"{declared - written} bytes of the file are missing. Nothing was "
                f"cached -- a short file would be served from the cache as if it were "
                f"the real structure. This is transient, so retry."
            )
        if declared is not None and written > declared:
            raise RetrievalError(
                f"Download of {url} does not match its own Content-Length: the server "
                f"declared {declared} bytes and {written} arrived. Nothing was cached, "
                f"because a body that does not match its description cannot be trusted "
                f"to be the file that was asked for."
            )
        os.replace(partial, destination)
    finally:
        partial.unlink(missing_ok=True)


def _write_json(path: Path, payload: Any) -> None:
    """Write JSON to ``path`` atomically, for the same reason as :func:`_download_to`."""
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(f"{path.name}.{os.getpid()}.part")
    try:
        partial.write_text(json.dumps(payload), encoding="utf-8")
        os.replace(partial, path)
    finally:
        partial.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Identifier normalization
# ---------------------------------------------------------------------------


def _normalize_accession(accession: str) -> str:
    """Uppercase and syntax-check a UniProt accession.

    Checking locally means an obvious typo produces an immediate, specific error
    instead of an opaque HTTP 400, and does not spend a request on EBI to find out.
    Anything this accepts but a service rejects still surfaces as a 400 downstream.

    Raises
    ------
    RetrievalError
        If ``accession`` is not valid UniProtKB accession syntax.
    """
    candidate = accession.strip().upper()
    if not _UNIPROT_ACCESSION_RE.match(candidate):
        raise RetrievalError(
            f"{accession!r} is not a UniProt accession. Expected something like "
            f"'P04637' or 'A0A123B4C5', optionally with an isoform suffix ('P04637-2'). "
            f"If you have a protein *name*, resolve it first with "
            f"resolve_uniprot_accession()."
        )
    return candidate


def _looks_like_accession(name: str) -> bool:
    """Return True if ``name`` is already UniProtKB accession syntax."""
    return bool(_UNIPROT_ACCESSION_RE.match(name.strip().upper()))


# ---------------------------------------------------------------------------
# AlphaFold DB
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AlphaFoldModel:
    """A downloaded AlphaFold model and the metadata that came with it.

    Attributes
    ----------
    path
        Local path to the downloaded coordinate file.
    accession
        The UniProt accession that was requested, normalized to uppercase.
    entry_id
        AFDB entry id of the model, e.g. ``"AF-P04637-F1"``.
    model_url
        The URL the file actually came from. Recorded because it encodes the model
        version, and because the whole point of this module is that the version is
        discovered rather than assumed.
    model_version
        AFDB release version of the file (``6`` at time of writing), or ``None`` if the
        API did not report one.
    mean_plddt
        Mean pLDDT over the model, from the API's ``globalMetricValue``. ``None`` if
        absent. A whole-model confidence summary: low values indicate a largely
        disordered protein.
    fraction_plddt_very_low, fraction_plddt_low, fraction_plddt_confident,
    fraction_plddt_very_high
        Fraction of residues in each AF2 confidence band (<50, 50-70, 70-90, >90).
        Per-residue pLDDT is in the file's B-factor column; these are the summary.
    sequence
        Sequence of the modelled span, as reported by the API.
    uniprot_length
        Length of the full UniProt sequence, which for a fragment is longer than
        ``sequence``.
    residue_range
        ``(start, stop)`` UniProt residue numbers covered by this model, inclusive.
    fragment_entry_ids
        Entry ids of all fragments of this accession, when there is more than one.
        Empty for the ordinary single-model case.
    from_cache
        True if the coordinate file was already on disk and no download happened.
    notes
        Human-readable anomalies: fragmentation, a cif fallback, a discarded corrupt
        cache entry. Never silently empty when something unusual happened.
    """

    path: Path
    accession: str
    entry_id: str
    model_url: str
    model_version: int | None = None
    mean_plddt: float | None = None
    fraction_plddt_very_low: float | None = None
    fraction_plddt_low: float | None = None
    fraction_plddt_confident: float | None = None
    fraction_plddt_very_high: float | None = None
    sequence: str = ""
    uniprot_length: int = 0
    residue_range: tuple[int, int] = (0, 0)
    fragment_entry_ids: tuple[str, ...] = ()
    from_cache: bool = False
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def n_modelled_residues(self) -> int:
        """Number of UniProt residues this model covers."""
        start, stop = self.residue_range
        return max(0, stop - start + 1)

    @property
    def is_fragment(self) -> bool:
        """True if this model covers only part of the UniProt sequence.

        The case v1 mishandled. When it is True, the model's sequence is *expected* to
        differ from the UniProt sequence, and comparing the two is not a validity check.
        """
        if bool(self.fragment_entry_ids):
            return True
        return 0 < self.n_modelled_residues < self.uniprot_length

    @property
    def fraction_low_confidence(self) -> float | None:
        """Fraction of residues with pLDDT below 70, or ``None`` if unreported.

        The band boundary DODO treats as disordered
        (:data:`dodo.constants.PLDDT_DISORDER_THRESHOLD`), so this is the quickest
        estimate of how much of the protein DODO will want to rebuild.
        """
        if self.fraction_plddt_very_low is None or self.fraction_plddt_low is None:
            return None
        return self.fraction_plddt_very_low + self.fraction_plddt_low


def _afdb_prediction(accession: str, *, timeout: float) -> list[dict[str, Any]]:
    """Fetch the AFDB prediction metadata list for ``accession``.

    Raises
    ------
    StructureNotFoundError
        On HTTP 404, or an empty result list: AFDB genuinely has no model.
    RetrievalError
        On HTTP 400 (malformed accession), any other HTTP status, or a transport
        failure.
    """
    url = AFDB_PREDICTION_API.format(accession=urllib.parse.quote(accession, safe=""))
    try:
        payload = _get_json(url, timeout=timeout)
    except _HttpStatusError as exc:
        if exc.status == 404:
            raise StructureNotFoundError(
                f"AlphaFold DB has no model for {accession}. Retrying will not help. "
                f"Entries are genuinely missing for some UniProt accessions, most often "
                f"proteins longer than {AFDB_FRAGMENT_LENGTH_LIMIT} residues -- titin "
                f"(Q8WZ42, 34350 residues) is the canonical example, and it 404s here "
                f"too. Check the accession at https://alphafold.ebi.ac.uk/, and if the "
                f"protein really is absent you will need a locally predicted structure."
            ) from exc
        if exc.status == 400:
            raise RetrievalError(
                f"AlphaFold DB rejected {accession!r} as a malformed identifier "
                f"(HTTP 400). This is a bad accession, not a missing structure and not "
                f"a network problem: check the spelling."
            ) from exc
        raise RetrievalError(
            f"AlphaFold DB returned HTTP {exc.status} for {accession}. Server-side "
            f"errors here are usually transient; retrying is reasonable."
        ) from exc

    if not isinstance(payload, list) or not payload:
        raise StructureNotFoundError(
            f"AlphaFold DB returned no predictions for {accession} (empty response). "
            f"There is no model to download."
        )
    entries = [entry for entry in payload if isinstance(entry, dict)]
    if not entries:
        raise RetrievalError(
            f"AlphaFold DB returned an unrecognized payload for {accession}: expected a "
            f"list of prediction objects, got {type(payload[0]).__name__} entries."
        )
    return entries


def _cached_metadata_problem(cached: Any) -> str | None:
    """Describe why cached AFDB metadata is unusable, or return None if it is fine.

    :func:`_write_json` only ever writes a non-empty list of prediction objects here, so
    anything else means the file was damaged (a partial write, a disk fault) or edited
    after DODO wrote it. Re-fetching is the right recovery, but doing it without saying
    so is the silent self-heal this package exists to avoid: valid-JSON-wrong-shape is
    the same class of anomaly as unparseable JSON and gets the same note.
    """
    if not isinstance(cached, list):
        return f"was valid JSON but a {type(cached).__name__}, not a list of predictions"
    if not cached:
        return "was an empty list, which is not something DODO writes"
    bad = [entry for entry in cached if not isinstance(entry, dict)]
    if bad:
        kinds = ", ".join(sorted({type(entry).__name__ for entry in bad}))
        return (
            f"held {len(bad)} of {len(cached)} entries that were not prediction objects ({kinds})"
        )
    return None


def _fragment_index(entry_id: str) -> int:
    """Return the ``-F<n>`` fragment number of an AFDB entry id, or 1 if it has none."""
    match = _FRAGMENT_SUFFIX_RE.search(entry_id)
    return int(match.group(1)) if match else 1


def _select_entries(accession: str, entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Pick the entries that really are models of ``accession``, in fragment order.

    Raises
    ------
    StructureNotFoundError
        If none of the returned entries is a model of this accession. That happens for
        long proteins: the request succeeds, but every entry is an alternative isoform
        with a different sequence. v1 took ``[0]`` here and modelled the wrong protein.
    """
    matches = [
        entry
        for entry in entries
        if str(entry.get("uniprotAccession", "")).strip().upper() == accession
    ]
    if not matches:
        available = sorted(
            {str(entry.get("uniprotAccession") or entry.get("entryId") or "?") for entry in entries}
        )
        raise StructureNotFoundError(
            f"AlphaFold DB returned {len(entries)} prediction(s) for the query "
            f"{accession}, but none of them models {accession} itself; they are for "
            f"{', '.join(available)}. Those are different sequences -- normally "
            f"alternative isoforms -- and using one as if it were {accession} would give "
            f"you a structure of the wrong protein. AFDB omits the canonical sequence of "
            f"proteins longer than {AFDB_FRAGMENT_LENGTH_LIMIT} residues, which is the "
            f"usual reason for this. To use one of the isoforms deliberately, pass its "
            f"accession (e.g. {available[0]!r}) instead."
        )
    return sorted(matches, key=lambda entry: _fragment_index(str(entry.get("entryId", ""))))


def _coerce_float(value: Any) -> float | None:
    """Return ``value`` as a float, or ``None`` if it is absent or not a number."""
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _coerce_int(value: Any, default: int = 0) -> int:
    """Return ``value`` as an int, or ``default`` if it is absent or not a number."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def fetch_alphafold_model(
    accession: str,
    *,
    timeout: float = 30.0,
    cache_dir: Path | None = None,
    refresh: bool = False,
) -> AlphaFoldModel:
    """Download the current AlphaFold model for a UniProt accession.

    The URL is resolved from the prediction API rather than constructed, so this keeps
    working when EBI bumps the model version. Both the API response and the coordinate
    file are cached, so a repeat call for the same accession makes no request at all --
    which matters because the regression suite fetches the same handful of accessions
    over and over.

    Parameters
    ----------
    accession
        UniProt accession, e.g. ``"P04637"``. Case-insensitive. An isoform accession
        (``"P04637-2"``) is accepted and fetches that isoform's model. Protein *names*
        are not accepted here; use :func:`resolve_uniprot_accession` first.
    timeout
        Per-request socket timeout in seconds.
    cache_dir
        Cache root. Defaults to :func:`default_cache_dir`.
    refresh
        Ignore cached metadata and coordinates and re-fetch. The cache key includes the
        model version, so a new AFDB release is normally picked up automatically --
        except that the *metadata* telling us the new version is itself cached. This is
        the escape hatch for that; deleting the cache directory has the same effect.

    Returns
    -------
    AlphaFoldModel
        The local path plus the API's confidence metrics and fragment information.

    Raises
    ------
    StructureNotFoundError
        AFDB has no model for this accession (HTTP 404), or has models only for other
        isoforms. Retrying will not help.
    RetrievalError
        The accession is malformed, the service returned an unexpected status, or the
        network failed. Network failures are transient and worth retrying.

    Warns
    -----
    UserWarning
        If the model is only a fragment of the protein. See
        :attr:`AlphaFoldModel.is_fragment`.

    Examples
    --------
    >>> model = fetch_alphafold_model("P04637")      # doctest: +SKIP
    >>> model.path.name                              # doctest: +SKIP
    'AF-P04637-F1-model_v6.pdb'
    >>> model.mean_plddt                             # doctest: +SKIP
    75.06
    """
    normalized = _normalize_accession(accession)
    directory = (default_cache_dir() if cache_dir is None else Path(cache_dir)) / _AFDB_CACHE_SUBDIR
    metadata_path = directory / f"{normalized}.prediction.json"
    notes: list[str] = []

    entries: list[dict[str, Any]] | None = None
    if not refresh and metadata_path.is_file():
        try:
            cached = json.loads(metadata_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            # Self-heal rather than fail, but say so: a truncated cache file is a
            # symptom (interrupted run, full disk) the user may want to know about.
            notes.append(
                f"Cached AFDB metadata at {metadata_path} was unreadable ({exc}); "
                f"it was re-fetched."
            )
        else:
            problem = _cached_metadata_problem(cached)
            if problem is None:
                entries = list(cached)
            else:
                notes.append(
                    f"Cached AFDB metadata at {metadata_path} {problem}; it was "
                    f"discarded and re-fetched."
                )

    if not entries:
        entries = _afdb_prediction(normalized, timeout=timeout)
        _write_json(metadata_path, entries)

    matches = _select_entries(normalized, entries)
    entry = matches[0]
    entry_id = str(entry.get("entryId") or f"AF-{normalized}-F1")

    model_url = str(entry.get("pdbUrl") or "")
    if not model_url:
        model_url = str(entry.get("cifUrl") or "")
        if not model_url:
            raise RetrievalError(
                f"AlphaFold DB reported entry {entry_id} for {normalized} but gave "
                f"neither a pdbUrl nor a cifUrl, so there is nothing to download."
            )
        notes.append(
            f"AFDB entry {entry_id} has no PDB-format file; downloaded the mmCIF "
            f"({model_url}) instead."
        )

    filename = Path(urllib.parse.urlparse(model_url).path).name
    if not filename:
        raise RetrievalError(
            f"AlphaFold DB gave an unusable model URL for {entry_id}: {model_url!r}"
        )

    # The remote filename carries the model version (AF-P04637-F1-model_v6.pdb), so a
    # release bump changes the cache key and a stale v4 file can never be served as a
    # hit for a v6 request.
    path = directory / filename
    from_cache = path.is_file() and path.stat().st_size > 0
    if refresh or not from_cache:
        try:
            _download_to(model_url, path, timeout=timeout)
        except _HttpStatusError as exc:
            if exc.status == 404:
                raise StructureNotFoundError(
                    f"AlphaFold DB advertised {model_url} for {entry_id} but the file is "
                    f"not there (HTTP 404). This usually means a release changed while "
                    f"cached metadata was still in use: retry with refresh=True."
                ) from exc
            raise RetrievalError(
                f"Downloading {model_url} failed with HTTP {exc.status}. Retrying is "
                f"reasonable if this persists."
            ) from exc
        from_cache = False

    uniprot_length = len(str(entry.get("uniprotSequence") or ""))
    start = _coerce_int(entry.get("uniprotStart"), 1)
    stop = _coerce_int(entry.get("uniprotEnd"), 0)
    sequence = str(entry.get("sequence") or "")
    if stop <= 0:
        stop = start + max(len(sequence), 1) - 1
    if uniprot_length <= 0:
        uniprot_length = max(stop, len(sequence))

    fragment_ids: tuple[str, ...] = ()
    if len(matches) > 1:
        fragment_ids = tuple(str(match.get("entryId") or "?") for match in matches)

    fragment_message = _fragment_message(
        accession=normalized,
        entry_id=entry_id,
        residue_range=(start, stop),
        uniprot_length=uniprot_length,
        fragment_ids=fragment_ids,
    )
    if fragment_message is not None:
        notes.append(fragment_message)

    model = AlphaFoldModel(
        path=path,
        accession=normalized,
        entry_id=entry_id,
        model_url=model_url,
        model_version=_coerce_int(entry.get("latestVersion")) or None,
        mean_plddt=_coerce_float(entry.get("globalMetricValue")),
        fraction_plddt_very_low=_coerce_float(entry.get("fractionPlddtVeryLow")),
        fraction_plddt_low=_coerce_float(entry.get("fractionPlddtLow")),
        fraction_plddt_confident=_coerce_float(entry.get("fractionPlddtConfident")),
        fraction_plddt_very_high=_coerce_float(entry.get("fractionPlddtVeryHigh")),
        sequence=sequence,
        uniprot_length=uniprot_length,
        residue_range=(start, stop),
        fragment_entry_ids=fragment_ids,
        from_cache=from_cache,
        notes=tuple(notes),
    )

    if fragment_message is not None:
        warnings.warn(fragment_message, UserWarning, stacklevel=2)

    return model


def _fragment_message(
    *,
    accession: str,
    entry_id: str,
    residue_range: tuple[int, int],
    uniprot_length: int,
    fragment_ids: tuple[str, ...],
) -> str | None:
    """Describe a fragmentary model, or return None if the model is complete.

    v1 detected this state indirectly -- the model sequence did not equal the UniProt
    sequence -- and raised ``'AF2 sequence does not match sequence retrieved from
    Uniprot!'``, which blames the data for something the database documents. The
    mismatch is expected for a fragment, so say what is actually going on.
    """
    start, stop = residue_range
    n_modelled = max(0, stop - start + 1)
    if not fragment_ids and not 0 < n_modelled < uniprot_length:
        return None
    message = (
        f"{entry_id} is a *fragment*: it covers UniProt residues {start}-{stop} of "
        f"{accession}, which is {uniprot_length} residues long. AlphaFold DB splits "
        f"predictions for proteins over {AFDB_FRAGMENT_LENGTH_LIMIT} residues into "
        f"overlapping fragments, so this model's sequence is expected to differ from the "
        f"full UniProt sequence -- that is not an error and not a corrupt file. "
    )
    if len(fragment_ids) > 1:
        message += f"The fragments are: {', '.join(fragment_ids)}. "
    return message + (
        "Region identification on this file describes the fragment only, and the "
        "numbering to interpret it against is the UniProt numbering above."
    )


def fetch_alphafold(
    accession: str,
    *,
    timeout: float = 30.0,
    cache_dir: Path | None = None,
    refresh: bool = False,
) -> Path:
    """Download the current AlphaFold model and return the local file path.

    Thin wrapper over :func:`fetch_alphafold_model` for callers that only want the
    file. Use the other one if you want the pLDDT summary or need to know whether the
    model is a fragment -- this signature cannot report either, and a fragment is only
    surfaced here as a :class:`UserWarning`.

    Parameters
    ----------
    accession
        UniProt accession, e.g. ``"P04637"``.
    timeout
        Per-request socket timeout in seconds.
    cache_dir
        Cache root. Defaults to :func:`default_cache_dir`.
    refresh
        Re-fetch even if cached.

    Returns
    -------
    Path
        Local path to the downloaded model.

    Raises
    ------
    StructureNotFoundError
        AFDB has no model for this accession.
    RetrievalError
        Malformed accession, unexpected status, or network failure.
    """
    return fetch_alphafold_model(
        accession, timeout=timeout, cache_dir=cache_dir, refresh=refresh
    ).path


# ---------------------------------------------------------------------------
# RCSB PDB
# ---------------------------------------------------------------------------


def fetch_pdb(
    pdb_id: str,
    *,
    format: str = "cif",
    timeout: float = 30.0,
    cache_dir: Path | None = None,
    refresh: bool = False,
) -> Path:
    """Download an experimental structure from RCSB.

    Parameters
    ----------
    pdb_id
        Four-character PDB ID (``"6KN7"``) or extended ID (``"pdb_00006kn7"``).
        Case-insensitive.
    format
        ``"cif"`` (default) or ``"pdb"``. ``"mmcif"`` is accepted as a synonym of
        ``"cif"``. mmCIF is the default because it is the only format in which large
        assemblies are distributed at all.
    timeout
        Per-request socket timeout in seconds.
    cache_dir
        Cache root. Defaults to :func:`default_cache_dir`.
    refresh
        Re-download even if cached.

    Returns
    -------
    Path
        Local path to the downloaded file.

    Raises
    ------
    UnsupportedFormatError
        If ``format`` is not one DODO can read.
    StructureNotFoundError
        If RCSB has no such entry, or no file of that format for it (HTTP 404).
    RetrievalError
        Malformed ID, unexpected status, or network failure.
    """
    extension = _FORMAT_ALIASES.get(format.strip().lower().lstrip("."))
    if extension is None:
        raise UnsupportedFormatError(
            f"format must be one of {sorted(set(_FORMAT_ALIASES))}, got {format!r}."
        )

    entry = pdb_id.strip().upper()
    if not _PDB_ID_RE.match(entry):
        raise RetrievalError(
            f"{pdb_id!r} is not a PDB ID. Expected four characters starting with a digit "
            f"(e.g. '6KN7') or an extended id (e.g. 'pdb_00006kn7')."
        )

    directory = (default_cache_dir() if cache_dir is None else Path(cache_dir)) / _PDB_CACHE_SUBDIR
    path = directory / f"{entry}.{extension}"
    if not refresh and path.is_file() and path.stat().st_size > 0:
        return path

    url = RCSB_DOWNLOAD_URL.format(entry=entry, extension=extension)
    try:
        _download_to(url, path, timeout=timeout)
    except _HttpStatusError as exc:
        if exc.status == 404:
            hint = ""
            if extension == "pdb":
                # Real case: entries with more than 99999 atoms or 62 chains cannot be
                # expressed in the fixed-width PDB format and RCSB publishes no .pdb
                # for them at all. Large EM assemblies -- exactly what DODO is pointed
                # at -- are the common victims.
                hint = (
                    " Note that RCSB does not distribute a PDB-format file for entries "
                    "too large for it (over 99999 atoms or 62 chains); try "
                    'format="cif".'
                )
            raise StructureNotFoundError(
                f"RCSB has no {extension} file for entry {entry} (HTTP 404). Retrying "
                f"will not help.{hint}"
            ) from exc
        raise RetrievalError(
            f"Downloading {url} failed with HTTP {exc.status}. This is usually transient."
        ) from exc
    return path


# ---------------------------------------------------------------------------
# UniProt
# ---------------------------------------------------------------------------


def _search_uniprot(name: str, *, timeout: float) -> str | None:
    """Return the top-hit accession for a free-text protein name, or None.

    Uses the UniProtKB REST search directly, so the base install can resolve names
    without the optional ``lookup`` extra (which pulls torch transitively).

    Raises
    ------
    RetrievalError
        On a transport failure, an unexpected status, or a response whose ``results``
        are not the documented shape. All three are distinguished from "no hits", which
        returns ``None``, because they call for different advice -- and because a hit
        DODO cannot read is not evidence that no hit exists.
    """
    query = urllib.parse.urlencode(
        {
            "query": name,
            "format": "json",
            # NOT size=1. Measured against the live service: size=1 does not return the
            # top relevance hit, it returns an apparently arbitrary one, while size>=2
            # returns the properly ranked list. 'human p53' gives Q89743 at size=1 and
            # P04637 at size=2; 'dnmt3a human' gives Q6PJ37 at size=1 and Q9Y6K1 at
            # size=2. Asking for a small page and taking the first row is the only way
            # to get the ranked answer, and getting this wrong means silently resolving
            # to the wrong protein.
            "size": _UNIPROT_SEARCH_PAGE_SIZE,
            # Only the accession is needed; the default field set is ~30 KB per hit.
            "fields": "accession",
        }
    )
    url = f"{UNIPROT_SEARCH_URL}?{query}"
    try:
        payload = _get_json(url, timeout=timeout)
    except _HttpStatusError as exc:
        if exc.status == 400:
            raise RetrievalError(
                f"UniProt rejected the search query {name!r} as malformed (HTTP 400)."
            ) from exc
        raise RetrievalError(
            f"UniProt search for {name!r} returned HTTP {exc.status}. Usually transient."
        ) from exc

    results = payload.get("results") if isinstance(payload, dict) else None
    if results is None or (isinstance(results, list) and not results):
        return None
    # Everything below this line is defensive against a payload that is JSON but not the
    # documented shape -- an unannounced API change, or a proxy's error page. Indexing
    # into it blind used to raise a bare AttributeError out of the public API, past the
    # RetrievalError contract every caller wraps this in.
    if not isinstance(results, list):
        raise RetrievalError(
            f"UniProt search for {name!r} returned an unrecognized payload: 'results' "
            f"is a {type(results).__name__}, not a list of hits."
        )
    top = results[0]
    if not isinstance(top, dict):
        raise RetrievalError(
            f"UniProt search for {name!r} returned an unrecognized payload: the first "
            f"hit is a {type(top).__name__} ({top!r}), not an object with a "
            f"'primaryAccession'."
        )
    accession = str(top.get("primaryAccession") or "").strip().upper()
    if not accession:
        # Not the same as "no hits": there *was* a hit and DODO could not read it.
        # Reporting it as no-match would send the caller off to the getSequence fallback
        # -- or to a MissingDependencyError -- for what is really a payload problem.
        raise RetrievalError(
            f"UniProt search for {name!r} returned an unrecognized payload: the top hit "
            f"has no primaryAccession (keys: {sorted(top)})."
        )
    return accession


class _DeadlineExceededError(Exception):
    """Internal: a bounded call did not return in time. Never leaves this module."""


def _call_with_deadline(work: Callable[[], Any], *, timeout: float, description: str) -> Any:
    """Run ``work()`` and give up waiting for it after ``timeout`` seconds.

    The escape hatch for a third-party call that does its own networking and exposes no
    timeout of its own. ``timeout`` is a required keyword of every helper in this module
    precisely so that an unbounded network call cannot exist; a callable that ignores it
    would be a hole in that, so the wait is bounded from the outside instead.

    The worker runs on a *daemon* thread. Python cannot cancel a thread, so an abandoned
    call keeps running until it finishes on its own -- but a daemon thread cannot delay
    interpreter exit, so it costs the caller nothing but a file descriptor.

    Raises
    ------
    _DeadlineExceededError
        If ``work`` has not returned after ``timeout`` seconds.
    Exception
        Whatever ``work`` raised, re-raised on the calling thread.
    """
    completed: list[Any] = []
    failed: list[BaseException] = []

    def target() -> None:
        try:
            completed.append(work())
        except BaseException as exc:  # re-raised on the calling thread below
            failed.append(exc)

    thread = threading.Thread(target=target, name=f"dodo-{description}", daemon=True)
    thread.start()
    thread.join(timeout)
    if thread.is_alive():
        raise _DeadlineExceededError(f"{description} did not return within {timeout:g} s")
    if failed:
        raise failed[0]
    return completed[0] if completed else None


def _accession_from_getsequence(name: str, *, timeout: float) -> str | None:
    """Resolve a name with the optional ``getSequence`` package, or return None.

    ``getseq`` returns ``[header, sequence]`` with a header like
    ``sp|P04637|P53_HUMAN ...``, but the exact shape has changed across releases, so
    the accession is pattern-matched out of the result rather than indexed by position:
    v1 did ``result[0].split('|')[1]`` and would raise IndexError on any other shape.

    This is the one path in the module that does not go through :func:`_urlopen`.
    ``getSequence`` issues its own HTTP requests with ``requests``, so DODO's timeout,
    its ``User-Agent`` and its host allowlist do not apply to them, and the hosts it
    talks to are whatever that package decides. What DODO *can* do is refuse to wait
    forever, which is what :func:`_call_with_deadline` is for -- the alternative is the
    unbounded network call this rewrite removed everywhere else.

    Parameters
    ----------
    name
        Free-text protein name.
    timeout
        Seconds to wait for ``getseq`` before giving up on it. Bounds DODO's wait, not
        the underlying socket: ``getseq`` takes no timeout argument in any release, so
        the request itself cannot be cut short.

    Raises
    ------
    ImportError
        If ``getSequence`` is not installed. The caller decides what that means.
    RetrievalError
        If ``getseq`` does not return within ``timeout`` seconds.
    """
    module = importlib.import_module("getSequence")
    getseq = getattr(module, "getseq", None)
    if getseq is None:
        return None
    try:
        result = _call_with_deadline(
            lambda: getseq(name), timeout=timeout, description="getSequence lookup"
        )
    except _DeadlineExceededError as exc:
        raise RetrievalError(
            f"getSequence did not resolve {name!r} within the {timeout:g} s timeout. It "
            f"performs its own HTTP requests, which DODO cannot interrupt, so the "
            f"abandoned lookup may keep running in the background until this process "
            f"exits. Retry, raise the timeout, or pass the accession directly."
        ) from exc
    except Exception:  # third-party package, undocumented failure modes
        return None
    for token in re.split(r"[|\s,\[\]'\"]+", str(result)):
        if _looks_like_accession(token):
            return token.strip().upper()
    return None


def resolve_uniprot_accession(name: str, *, timeout: float = 30.0) -> str:
    """Resolve a protein name to a UniProt accession.

    Resolution order, in precedence:

    1. If ``name`` is already accession syntax it is returned uppercased, with no
       network request at all.
    2. The UniProtKB REST search, top hit. Standard library only, so the base install
       can do this.
    3. The optional ``getSequence`` package, only if step 2 found nothing and only if
       it happens to be installed. It is imported lazily and inside this function: it
       depends on torch transitively, and a module-scope import would make even
       ``dodo --help`` pay several seconds of framework startup.

    The search takes the top hit, exactly as v1 did, so an ambiguous name gives an
    arbitrary-but-plausible answer -- pass an accession if you need certainty.

    Parameters
    ----------
    name
        Protein name (``"human p53"``) or an accession (``"P04637"``).
    timeout
        Per-request socket timeout in seconds, for step 2. Step 3 is a third-party
        package that issues its own requests and accepts no timeout argument, so there
        ``timeout`` bounds how long DODO waits for an answer rather than the socket
        itself; see :func:`_accession_from_getsequence`.

    Returns
    -------
    str
        The accession, uppercased.

    Raises
    ------
    RetrievalError
        If the name cannot be resolved, the network failed, or the ``getSequence``
        fallback did not answer within ``timeout``.
    MissingDependencyError
        If the built-in search found nothing and ``getSequence`` -- which might have
        resolved it -- is not installed.

    Warns
    -----
    UserWarning
        If the answer came from step 3. That path resolves free text with a third
        party's own ranking, against hosts outside DODO's allowlist, and returns a
        plausible accession even for a name with no real match, so the caller is told
        which resolver produced the answer rather than left to assume.

    Examples
    --------
    >>> resolve_uniprot_accession("human p53")   # doctest: +SKIP
    'P04637'
    """
    if not name or not name.strip():
        raise RetrievalError("Cannot resolve an empty protein name.")
    if _looks_like_accession(name):
        return name.strip().upper()

    transport_failure: RetrievalError | None = None
    hit: str | None = None
    try:
        hit = _search_uniprot(name, timeout=timeout)
    except RetrievalError as exc:
        transport_failure = exc
    if hit is not None:
        return hit

    try:
        fallback = _accession_from_getsequence(name, timeout=timeout)
    except RetrievalError as exc:
        # The fallback ran out of time. If the built-in resolver had already failed on
        # the network, that is the cause worth reporting; a hung optional package is a
        # symptom of the same dead connection.
        if transport_failure is not None:
            raise transport_failure from exc
        raise
    except ImportError as exc:
        if transport_failure is not None:
            # The built-in resolver never got an answer *because the network failed*.
            # Blaming a missing optional package here would send the user to install
            # torch to fix their wifi.
            raise transport_failure from exc
        raise MissingDependencyError(
            package="getSequence",
            purpose=(
                f"Resolving the protein name {name!r} to a UniProt accession, after "
                f"DODO's built-in UniProt search returned no match for it"
            ),
            extra="lookup",
        ) from exc

    if fallback is not None:
        warnings.warn(
            f"Resolved the name {name!r} to {fallback} with the optional getSequence "
            f"package, not with DODO's own UniProt search, which found no match for it. "
            f"getSequence applies its own ranking, queries services outside DODO's "
            f"allowlist, and returns a plausible-looking accession even for a name that "
            f"has no real match, so confirm that {fallback} is the protein you meant "
            f"before building anything from it.",
            UserWarning,
            stacklevel=2,
        )
        return fallback
    if transport_failure is not None:
        raise transport_failure
    raise RetrievalError(
        f"No UniProt entry matches {name!r}. Try a more specific name (for example "
        f"'human p53' rather than 'p53'), or pass the accession directly."
    )


def fetch_uniprot_sequence(accession: str, *, timeout: float = 30.0) -> str:
    """Fetch the amino acid sequence of a UniProt entry.

    Parameters
    ----------
    accession
        UniProt accession, optionally with an isoform suffix. Case-insensitive.
    timeout
        Socket timeout in seconds.

    Returns
    -------
    str
        The one-letter sequence, uppercased, with no whitespace.

    Raises
    ------
    StructureNotFoundError
        If UniProt has no such entry (HTTP 404). Retrying will not help.
    RetrievalError
        If the accession is malformed, the response is not FASTA, or the network
        failed.

    Notes
    -----
    Deliberately not cached. It is a few hundred bytes, it changes when UniProt
    revises an entry, and it is the sequence DODO checks a downloaded model against --
    caching it would let a stale copy outlive the model it is compared with.
    """
    normalized = _normalize_accession(accession)
    url = UNIPROT_FASTA_URL.format(accession=urllib.parse.quote(normalized, safe="-"))
    try:
        raw = _get_bytes(url, timeout=timeout, accept="text/plain")
    except _HttpStatusError as exc:
        if exc.status == 404:
            raise StructureNotFoundError(
                f"UniProt has no entry {normalized} (HTTP 404). Retrying will not help; "
                f"the accession may be obsolete or mistyped."
            ) from exc
        if exc.status == 400:
            raise RetrievalError(
                f"UniProt rejected {normalized!r} as a malformed accession (HTTP 400)."
            ) from exc
        raise RetrievalError(
            f"UniProt returned HTTP {exc.status} for {normalized}. Usually transient."
        ) from exc

    text = raw.decode("utf-8", errors="replace")
    if not text.lstrip().startswith(">"):
        preview = text[:120]
        raise RetrievalError(
            f"Expected FASTA from {url} but got something else. First characters: {preview!r}"
        )
    lines = [line.strip() for line in io.StringIO(text)]
    sequence = "".join(line for line in lines if line and not line.startswith(">")).upper()
    if not sequence:
        raise RetrievalError(
            f"UniProt returned a FASTA record for {normalized} with no sequence in it."
        )
    if not sequence.isalpha():
        offenders = sorted({character for character in sequence if not character.isalpha()})
        raise RetrievalError(
            f"UniProt sequence for {normalized} contains non-letter characters "
            f"{offenders}, so it is not a usable amino acid sequence."
        )
    return sequence
