"""Tests for remote structure and sequence retrieval.

The bulk of these run offline: every network read in :mod:`dodo.io.fetch` goes through
``fetch._urlopen``, which is monkeypatched here with a recording transport that serves
canned responses. That seam is what makes the interesting cases -- a 404 from AFDB, a
malformed accession, a timeout, a fragmented protein -- testable at all, since none of
them can be provoked on demand against a live service.

The handful of tests that do hit the network are marked ``network`` so CI can deselect
them with ``-m 'not network'``.
"""

from __future__ import annotations

import email.message
import http.client
import io
import json
import sys
import threading
import time
import types
import urllib.error
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pytest

from dodo.exceptions import (
    MissingDependencyError,
    RetrievalError,
    StructureNotFoundError,
    UnsupportedFormatError,
)
from dodo.io import fetch

# ---------------------------------------------------------------------------
# Offline transport
# ---------------------------------------------------------------------------

MODEL_BODY = b"ATOM      1  N   MET A   1      10.000  10.000  10.000  1.00 80.00           N\n"


class FakeResponse:
    """Minimal stand-in for the object ``urlopen`` returns.

    Carries a ``Content-Length`` header, because every service DODO talks to declares
    one for a static file and it is the only thing that distinguishes a complete body
    from a truncated one. A fake without it cannot exercise the truncation check at all.
    """

    def __init__(self, body: bytes, status: int = 200) -> None:
        self._buffer = io.BytesIO(body)
        self.status = status
        self.declared_length: int | None = len(body)
        #: Verbatim header value, for testing a header that is not a number at all.
        self.raw_content_length: str | None = None

    @property
    def headers(self) -> email.message.Message:
        message = email.message.Message()
        value = self.raw_content_length
        if value is None and self.declared_length is not None:
            value = str(self.declared_length)
        if value is not None:
            message["Content-Length"] = value
        return message

    def read(self, size: int | None = None) -> bytes:
        return self._buffer.read() if size is None else self._buffer.read(size)

    def __enter__(self) -> FakeResponse:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self._buffer.close()


class BrokenResponse(FakeResponse):
    """A response that dies part-way through the body, like a dropped connection."""

    def read(self, size: int | None = None) -> bytes:
        raise OSError("connection reset by peer")


class TruncatedResponse(FakeResponse):
    """A response that promises ``Content-Length`` bytes and delivers fewer.

    The failure mode a ``BytesIO``-backed fake structurally cannot produce: ``read()``
    never raises, it just returns ``b""`` early, so ``shutil.copyfileobj`` completes
    "successfully" on a body that lost most of its content. This is what a connection
    closed mid-transfer looks like to the reader when nothing checks the declared size.
    """

    def __init__(self, full_body: bytes, delivered: int) -> None:
        super().__init__(full_body[:delivered])
        self.declared_length = len(full_body)


class HeaderlessResponse(FakeResponse):
    """A response with no ``Content-Length`` at all, as in a chunked transfer."""

    def __init__(self, body: bytes, status: int = 200) -> None:
        super().__init__(body, status)
        self.declared_length = None


class IncompleteReadResponse(FakeResponse):
    """A response whose body read fails with the stdlib's own truncation error.

    ``http.client.IncompleteRead`` is what real urllib raises when the connection closes
    before ``Content-Length`` bytes arrive. It is an ``HTTPException``, *not* an
    ``OSError``, which is exactly why it used to escape the mid-body handlers.
    """

    def __init__(self, partial: bytes = b"MEEP", expected: int = 389) -> None:
        super().__init__(partial)
        self._partial = partial
        self._expected = expected

    def read(self, size: int | None = None) -> bytes:
        raise http.client.IncompleteRead(self._partial, self._expected)


@dataclass(frozen=True)
class Call:
    url: str
    timeout: float


Handler = Callable[[str], FakeResponse]


class Transport:
    """Records every request and answers it from a routing table.

    Routes are matched by substring against the URL, in insertion order. An unmatched
    URL is a test bug, not a pass -- a silently unmatched route would let a test claim
    to exercise a code path it never reached.
    """

    def __init__(self, routes: dict[str, Handler]) -> None:
        self.routes = routes
        self.calls: list[Call] = []

    def __call__(self, request: Any, timeout: float) -> FakeResponse:
        url = request.full_url
        self.calls.append(Call(url, timeout))
        assert request.get_header("User-agent", "").startswith("dodo/")
        for fragment, handler in self.routes.items():
            if fragment in url:
                return handler(url)
        raise AssertionError(f"test made an unrouted request to {url}")

    @property
    def urls(self) -> list[str]:
        return [call.url for call in self.calls]


def json_route(payload: Any) -> Handler:
    """Serve ``payload`` as JSON."""
    return lambda _url: FakeResponse(json.dumps(payload).encode())


def body_route(body: bytes) -> Handler:
    """Serve raw bytes."""
    return lambda _url: FakeResponse(body)


def status_route(status: int, body: bytes = b"{}") -> Handler:
    """Fail with an HTTP status."""

    def handler(url: str) -> FakeResponse:
        raise urllib.error.HTTPError(
            url, status, "canned", email.message.Message(), io.BytesIO(body)
        )

    return handler


def raising_route(exc: BaseException) -> Handler:
    """Fail with a transport-level exception."""

    def handler(_url: str) -> FakeResponse:
        raise exc

    return handler


def broken_route() -> Handler:
    """Succeed, then die mid-body."""
    return lambda _url: BrokenResponse(b"")


def truncated_route(full_body: bytes, delivered: int) -> Handler:
    """Promise ``full_body`` but deliver only its first ``delivered`` bytes."""
    return lambda _url: TruncatedResponse(full_body, delivered)


def headerless_route(body: bytes) -> Handler:
    """Serve raw bytes with no ``Content-Length``."""
    return lambda _url: HeaderlessResponse(body)


def incomplete_read_route() -> Handler:
    """Succeed, then raise the stdlib's ``IncompleteRead`` out of the body read."""
    return lambda _url: IncompleteReadResponse()


@pytest.fixture
def transport(monkeypatch: pytest.MonkeyPatch) -> Callable[[dict[str, Handler]], Transport]:
    def install(routes: dict[str, Handler]) -> Transport:
        fake = Transport(routes)
        monkeypatch.setattr(fetch, "_urlopen", fake)
        return fake

    return install


InstallTransport = Callable[[dict[str, Handler]], Transport]


# ---------------------------------------------------------------------------
# Canned AFDB payloads
# ---------------------------------------------------------------------------


def afdb_entry(
    accession: str = "P04637",
    *,
    entry_id: str | None = None,
    uniprot_start: int = 1,
    uniprot_end: int = 393,
    uniprot_length: int | None = None,
    version: int = 6,
) -> dict[str, Any]:
    """One AFDB prediction record, shaped like the live API's.

    Field names and values copied from a real response for P04637.
    """
    entry = entry_id or f"AF-{accession}-F1"
    length = uniprot_length if uniprot_length is not None else uniprot_end
    return {
        "entryId": entry,
        "uniprotAccession": accession,
        "uniprotId": "P53_HUMAN",
        "uniprotStart": uniprot_start,
        "uniprotEnd": uniprot_end,
        "uniprotSequence": "M" * length,
        "sequence": "M" * (uniprot_end - uniprot_start + 1),
        "pdbUrl": f"https://alphafold.ebi.ac.uk/files/{entry}-model_v{version}.pdb",
        "cifUrl": f"https://alphafold.ebi.ac.uk/files/{entry}-model_v{version}.cif",
        "latestVersion": version,
        "globalMetricValue": 75.06,
        "fractionPlddtVeryLow": 0.298,
        "fractionPlddtLow": 0.104,
        "fractionPlddtConfident": 0.071,
        "fractionPlddtVeryHigh": 0.527,
    }


def afdb_routes(entries: list[dict[str, Any]]) -> dict[str, Handler]:
    return {
        "api/prediction": json_route(entries),
        "alphafold.ebi.ac.uk/files/": body_route(MODEL_BODY),
    }


# ---------------------------------------------------------------------------
# Cache directory
# ---------------------------------------------------------------------------


class TestDefaultCacheDir:
    """Platform conventions, plus the generation segment.

    Every path ends in ``fetch.CACHE_GENERATION``. That segment exists so a cache written
    by a DODO old enough to predate download length-verification is abandoned rather than
    trusted: such an entry may be a silently truncated fragment, nothing recorded its
    expected size, and the cache-hit gate only checks that a file is non-empty -- which a
    fragment satisfies. Re-downloading is cheap; serving corruption forever is not.
    """

    def test_macos(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("DODO_CACHE_DIR", raising=False)
        monkeypatch.setattr(sys, "platform", "darwin")
        expected = Path.home() / "Library" / "Caches" / "dodo" / fetch.CACHE_GENERATION
        assert fetch.default_cache_dir() == expected

    def test_linux_respects_xdg(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("DODO_CACHE_DIR", raising=False)
        monkeypatch.setattr(sys, "platform", "linux")
        monkeypatch.setenv("XDG_CACHE_HOME", "/somewhere/cache")
        assert fetch.default_cache_dir() == Path("/somewhere/cache/dodo") / fetch.CACHE_GENERATION

    def test_linux_without_xdg(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("DODO_CACHE_DIR", raising=False)
        monkeypatch.delenv("XDG_CACHE_HOME", raising=False)
        monkeypatch.setattr(sys, "platform", "linux")
        assert fetch.default_cache_dir() == Path.home() / ".cache" / "dodo" / fetch.CACHE_GENERATION

    def test_windows(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.delenv("DODO_CACHE_DIR", raising=False)
        monkeypatch.setattr(sys, "platform", "win32")
        monkeypatch.setenv("LOCALAPPDATA", r"C:\Users\x\AppData\Local")
        assert fetch.default_cache_dir().parts[-3:] == (
            "dodo",
            "Cache",
            fetch.CACHE_GENERATION,
        )

    def test_env_override_wins(self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
        monkeypatch.setenv("DODO_CACHE_DIR", str(tmp_path / "elsewhere"))
        assert fetch.default_cache_dir() == tmp_path / "elsewhere" / fetch.CACHE_GENERATION

    def test_override_cannot_resurrect_an_older_generation(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """The generation applies to the override too, not just the platform defaults."""
        monkeypatch.setenv("DODO_CACHE_DIR", str(tmp_path))
        assert fetch.default_cache_dir().name == fetch.CACHE_GENERATION


# ---------------------------------------------------------------------------
# AlphaFold: the URL must come from the API
# ---------------------------------------------------------------------------


class TestAlphaFoldUrlResolution:
    def test_downloads_the_url_the_api_reports(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        calls = transport(afdb_routes([afdb_entry()]))
        path = fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert path.read_bytes() == MODEL_BODY
        assert path.name == "AF-P04637-F1-model_v6.pdb"
        assert calls.urls == [
            "https://alphafold.ebi.ac.uk/api/prediction/P04637",
            "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v6.pdb",
        ]

    def test_never_constructs_a_versioned_url(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        """The shipped bug: a hardcoded model_v4 URL, now a 404 for every accession.

        A future release will be model_v7 and this test still holds, because the version
        is only ever read out of the API response.
        """
        calls = transport(afdb_routes([afdb_entry(version=7)]))
        path = fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert path.name.endswith("model_v7.pdb")
        assert not any("model_v4" in url or "model_v5" in url for url in calls.urls)

    def test_accession_is_case_insensitive(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport(afdb_routes([afdb_entry()]))
        model = fetch.fetch_alphafold_model("p04637", cache_dir=tmp_path)
        assert model.accession == "P04637"

    def test_cif_fallback_is_noted_not_silent(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        entry = afdb_entry()
        del entry["pdbUrl"]
        transport(afdb_routes([entry]))
        model = fetch.fetch_alphafold_model("P04637", cache_dir=tmp_path)
        assert model.path.suffix == ".cif"
        assert any("no PDB-format file" in note for note in model.notes)

    def test_no_downloadable_url_raises(self, transport: InstallTransport, tmp_path: Path) -> None:
        entry = afdb_entry()
        del entry["pdbUrl"]
        del entry["cifUrl"]
        transport({"api/prediction": json_route([entry])})
        with pytest.raises(RetrievalError, match="neither a pdbUrl nor a cifUrl"):
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)

    def test_url_from_the_response_must_be_an_expected_host(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        """``pdbUrl`` is remote data and does not get to redirect DODO anywhere it likes."""
        entry = afdb_entry()
        entry["pdbUrl"] = "http://evil.example.com/AF-P04637-F1-model_v6.pdb"
        transport({"api/prediction": json_route([entry])})
        with pytest.raises(RetrievalError, match="only makes https requests"):
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)


class TestAlphaFoldMetadata:
    def test_plddt_metrics_are_returned(self, transport: InstallTransport, tmp_path: Path) -> None:
        transport(afdb_routes([afdb_entry()]))
        model = fetch.fetch_alphafold_model("P04637", cache_dir=tmp_path)
        assert model.mean_plddt == pytest.approx(75.06)
        assert model.fraction_plddt_very_high == pytest.approx(0.527)
        # 0.298 very low + 0.104 low: the fraction below the pLDDT 70 disorder band.
        assert model.fraction_low_confidence == pytest.approx(0.402)
        assert model.model_version == 6
        assert model.entry_id == "AF-P04637-F1"

    def test_missing_metrics_are_none_not_zero(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        entry = afdb_entry()
        for key in ("globalMetricValue", "fractionPlddtVeryLow", "fractionPlddtLow"):
            del entry[key]
        transport(afdb_routes([entry]))
        model = fetch.fetch_alphafold_model("P04637", cache_dir=tmp_path)
        assert model.mean_plddt is None
        assert model.fraction_low_confidence is None

    def test_complete_model_is_not_a_fragment(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport(afdb_routes([afdb_entry()]))
        model = fetch.fetch_alphafold_model("P04637", cache_dir=tmp_path)
        assert not model.is_fragment
        assert model.residue_range == (1, 393)
        assert model.notes == ()


# ---------------------------------------------------------------------------
# AlphaFold: failure modes, each distinct
# ---------------------------------------------------------------------------


class TestAlphaFoldFailures:
    def test_404_means_no_model_and_says_retrying_will_not_help(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        """Verified live: titin (Q8WZ42) 404s, because AFDB has no model for it."""
        transport({"api/prediction": status_route(404)})
        with pytest.raises(StructureNotFoundError) as excinfo:
            fetch.fetch_alphafold("Q8WZ42", cache_dir=tmp_path)
        message = str(excinfo.value)
        assert "Retrying will not help" in message
        assert "longer than 2700 residues" in message

    def test_400_is_a_bad_accession_not_a_missing_structure(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport({"api/prediction": status_route(400, b'{"error":"Invalid identifier"}')})
        with pytest.raises(RetrievalError) as excinfo:
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert not isinstance(excinfo.value, StructureNotFoundError)
        assert "malformed identifier" in str(excinfo.value)

    def test_malformed_accession_is_rejected_without_a_request(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        calls = transport({})
        with pytest.raises(RetrievalError, match="not a UniProt accession"):
            fetch.fetch_alphafold("not-an-accession", cache_dir=tmp_path)
        assert calls.calls == []

    def test_timeout_is_transient_and_says_so(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport({"api/prediction": raising_route(TimeoutError("timed out"))})
        with pytest.raises(RetrievalError) as excinfo:
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path, timeout=5.0)
        message = str(excinfo.value)
        assert "transient network failure" in message
        assert "timeout was 5 s" in message

    def test_connection_failure_is_transient(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport({"api/prediction": raising_route(urllib.error.URLError("no route to host"))})
        with pytest.raises(RetrievalError, match="transient network failure"):
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)

    def test_every_request_carries_the_timeout(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        """v1 passed no timeout at all, so a hung connection hung DODO forever."""
        calls = transport(afdb_routes([afdb_entry()]))
        fetch.fetch_alphafold("P04637", cache_dir=tmp_path, timeout=12.5)
        assert calls.calls
        assert all(call.timeout == 12.5 for call in calls.calls)

    def test_nonpositive_timeout_is_rejected(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport({})
        with pytest.raises(RetrievalError, match="timeout must be positive"):
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path, timeout=0)

    def test_server_error_mentions_retrying(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport({"api/prediction": status_route(503)})
        with pytest.raises(RetrievalError, match="HTTP 503"):
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)

    def test_non_json_response_is_reported_as_such(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport({"api/prediction": body_route(b"<html>maintenance</html>")})
        with pytest.raises(RetrievalError, match="did not return JSON"):
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)

    def test_empty_prediction_list_is_not_found(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport({"api/prediction": json_route([])})
        with pytest.raises(StructureNotFoundError, match="no predictions"):
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)

    def test_empty_body_is_not_cached_as_a_valid_file(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport({"api/prediction": json_route([afdb_entry()]), "files/": body_route(b"")})
        with pytest.raises(RetrievalError, match="empty file"):
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert list((tmp_path / "alphafold").glob("*.pdb")) == []

    def test_interrupted_download_leaves_nothing_behind(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        """A truncated cache file would be reused forever as a cache hit."""
        transport({"api/prediction": json_route([afdb_entry()]), "files/": broken_route()})
        with pytest.raises(RetrievalError, match="failed part-way through"):
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert list((tmp_path / "alphafold").glob("*.pdb")) == []
        assert list((tmp_path / "alphafold").glob("*.part")) == []


# ---------------------------------------------------------------------------
# Truncated bodies. The one failure mode that becomes permanent if it is missed:
# a short file committed to the cache is served as a hit forever.
# ---------------------------------------------------------------------------

FULL_MODEL_BODY = MODEL_BODY * 60


class TestTruncatedDownloads:
    def test_short_body_is_rejected_not_cached(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        """Content-Length says 4980 bytes, the connection delivers 100 of them.

        ``read()`` never fails here -- it just returns ``b""`` early -- so nothing but
        the declared length can reveal the loss. Committing this to the cache turns one
        dropped connection into a permanently wrong structure.
        """
        transport(
            {
                "api/prediction": json_route([afdb_entry()]),
                "files/": truncated_route(FULL_MODEL_BODY, 100),
            }
        )
        with pytest.raises(RetrievalError) as excinfo:
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        message = str(excinfo.value)
        assert "100" in message and str(len(FULL_MODEL_BODY)) in message
        assert list((tmp_path / "alphafold").glob("*.pdb")) == []
        assert list((tmp_path / "alphafold").glob("*.part")) == []

    def test_a_truncated_body_never_becomes_a_cache_hit(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        """The second call must not be served the corpse of the first."""
        calls = transport(
            {
                "api/prediction": json_route([afdb_entry()]),
                "files/": truncated_route(FULL_MODEL_BODY, 100),
            }
        )
        for _ in range(2):
            with pytest.raises(RetrievalError):
                fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert list((tmp_path / "alphafold").glob("*.pdb")) == []
        # Metadata is cached, so the second attempt is the file only: 2 + 1 requests.
        assert len(calls.calls) == 3, "a rejected download must be retried, not reused"

    def test_fetch_pdb_rejects_a_short_body(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        full = b"data_1UBQ\n" + b"ATOM\n" * 500
        transport({"files.rcsb.org": truncated_route(full, 40)})
        with pytest.raises(RetrievalError, match="incomplete"):
            fetch.fetch_pdb("1UBQ", cache_dir=tmp_path)
        assert list((tmp_path / "pdb").glob("*")) == []

    def test_an_overlong_body_is_also_rejected(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        """Bytes beyond Content-Length mean the stream is not what was described."""

        def overlong(_url: str) -> FakeResponse:
            response = FakeResponse(FULL_MODEL_BODY)
            response.declared_length = 10
            return response

        transport({"api/prediction": json_route([afdb_entry()]), "files/": overlong})
        with pytest.raises(RetrievalError, match=r"incomplete|does not match"):
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert list((tmp_path / "alphafold").glob("*.pdb")) == []

    def test_a_complete_body_of_the_declared_length_is_accepted(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport(afdb_routes([afdb_entry()]))
        path = fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert path.read_bytes() == MODEL_BODY

    def test_a_body_with_no_declared_length_is_accepted(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        """Chunked responses declare no length, so there is nothing to check against.

        Not a hole: chunked framing is self-terminating, so a chunked body that stops
        early fails inside ``http.client`` and arrives here as ``IncompleteRead`` --
        covered by :class:`TestMidBodyTruncationErrors`. That matters because chunked is
        what the live download endpoints actually use: measured against the services,
        AFDB and RCSB both send ``Transfer-Encoding: chunked`` and no ``Content-Length``,
        so between them these two tests cover one real endpoint each.
        """
        transport(
            {
                "api/prediction": json_route([afdb_entry()]),
                "files/": headerless_route(MODEL_BODY),
            }
        )
        path = fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert path.read_bytes() == MODEL_BODY

    def test_unparseable_content_length_does_not_break_the_download(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        def garbled(_url: str) -> FakeResponse:
            response = FakeResponse(MODEL_BODY)
            response.raw_content_length = "banana"
            return response

        transport({"api/prediction": json_route([afdb_entry()]), "files/": garbled})
        assert fetch.fetch_alphafold("P04637", cache_dir=tmp_path).read_bytes() == MODEL_BODY


class TestMidBodyTruncationErrors:
    """``IncompleteRead`` is an ``HTTPException``, not an ``OSError``.

    It is the canonical stdlib error for "the connection closed before Content-Length
    bytes arrived", so it is precisely what these handlers exist to translate. Letting
    it escape breaks the documented ``Raises: RetrievalError`` contract of every public
    fetcher, and with it the "this is transient, retry" advice.
    """

    def test_download_incomplete_read_becomes_a_retrieval_error(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport({"api/prediction": json_route([afdb_entry()]), "files/": incomplete_read_route()})
        with pytest.raises(RetrievalError, match="failed part-way through"):
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert list((tmp_path / "alphafold").glob("*.pdb")) == []
        assert list((tmp_path / "alphafold").glob("*.part")) == []

    def test_sequence_incomplete_read_becomes_a_retrieval_error(
        self, transport: InstallTransport
    ) -> None:
        transport({"rest.uniprot.org": incomplete_read_route()})
        with pytest.raises(RetrievalError, match="part-way through"):
            fetch.fetch_uniprot_sequence("P04637")

    def test_metadata_incomplete_read_becomes_a_retrieval_error(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport({"api/prediction": incomplete_read_route()})
        with pytest.raises(RetrievalError, match="part-way through"):
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)

    def test_an_http_exception_from_urlopen_is_a_retrieval_error(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport(
            {"api/prediction": raising_route(http.client.BadStatusLine("garbage from server"))}
        )
        with pytest.raises(RetrievalError):
            fetch.fetch_alphafold("P04637", cache_dir=tmp_path)


# ---------------------------------------------------------------------------
# AlphaFold: fragments and isoforms
# ---------------------------------------------------------------------------


class TestFragments:
    @staticmethod
    def fragment_entries() -> list[dict[str, Any]]:
        # Shaped after a >2700-residue protein: two overlapping fragments of a 5183
        # residue sequence.
        return [
            afdb_entry(
                "Q5T4S7",
                entry_id="AF-Q5T4S7-F1",
                uniprot_start=1,
                uniprot_end=2700,
                uniprot_length=5183,
            ),
            afdb_entry(
                "Q5T4S7",
                entry_id="AF-Q5T4S7-F2",
                uniprot_start=2501,
                uniprot_end=5183,
                uniprot_length=5183,
            ),
        ]

    def test_fragment_is_detected_and_explained(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport(afdb_routes(self.fragment_entries()))
        with pytest.warns(UserWarning, match="is a \\*fragment\\*") as warned:
            model = fetch.fetch_alphafold_model("Q5T4S7", cache_dir=tmp_path)
        assert model.is_fragment
        assert model.entry_id == "AF-Q5T4S7-F1"
        assert model.residue_range == (1, 2700)
        assert model.uniprot_length == 5183
        assert model.fragment_entry_ids == ("AF-Q5T4S7-F1", "AF-Q5T4S7-F2")
        assert any("fragment" in note for note in model.notes)
        message = str(warned[0].message)
        # v1's message for this state blamed the data. The mismatch is expected.
        assert "does not match sequence retrieved from Uniprot" not in message
        assert "expected to differ" in message
        assert "AF-Q5T4S7-F2" in message

    def test_first_fragment_is_chosen_regardless_of_response_order(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport(afdb_routes(list(reversed(self.fragment_entries()))))
        with pytest.warns(UserWarning):
            model = fetch.fetch_alphafold_model("Q5T4S7", cache_dir=tmp_path)
        assert model.entry_id == "AF-Q5T4S7-F1"

    def test_partial_coverage_alone_counts_as_a_fragment(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        entries = [
            afdb_entry("Q5T4S7", uniprot_start=1, uniprot_end=2700, uniprot_length=5183),
        ]
        transport(afdb_routes(entries))
        with pytest.warns(UserWarning, match="fragment"):
            model = fetch.fetch_alphafold_model("Q5T4S7", cache_dir=tmp_path)
        assert model.is_fragment
        assert model.fragment_entry_ids == ()

    def test_isoform_only_response_does_not_silently_model_the_wrong_protein(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        """Verified live for Q5T4S7: the API returns only isoform entries.

        ``[0]`` there is a 212-residue isoform, so v1's indexing would have handed back
        a model of a different sequence entirely.
        """
        entries = [
            afdb_entry("Q5T4S7-6", entry_id="AF-Q5T4S7-6-F1", uniprot_end=212),
            afdb_entry("Q5T4S7-5", entry_id="AF-Q5T4S7-5-F1", uniprot_end=2456),
        ]
        calls = transport(afdb_routes(entries))
        with pytest.raises(StructureNotFoundError) as excinfo:
            fetch.fetch_alphafold("Q5T4S7", cache_dir=tmp_path)
        message = str(excinfo.value)
        assert "Q5T4S7-5" in message and "Q5T4S7-6" in message
        assert "longer than 2700 residues" in message
        assert len(calls.calls) == 1, "must not download the wrong protein's model"

    def test_an_isoform_can_be_requested_explicitly(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        entries = [afdb_entry("Q5T4S7-6", entry_id="AF-Q5T4S7-6-F1", uniprot_end=212)]
        transport(afdb_routes(entries))
        model = fetch.fetch_alphafold_model("Q5T4S7-6", cache_dir=tmp_path)
        assert model.accession == "Q5T4S7-6"
        assert not model.is_fragment


# ---------------------------------------------------------------------------
# Caching
# ---------------------------------------------------------------------------


class TestCache:
    def test_second_call_makes_no_request(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        calls = transport(afdb_routes([afdb_entry()]))
        first = fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert len(calls.calls) == 2  # metadata, then the model file
        second = fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert second == first
        assert len(calls.calls) == 2, "a cache hit must not touch the network"

    def test_cache_hit_is_reported(self, transport: InstallTransport, tmp_path: Path) -> None:
        transport(afdb_routes([afdb_entry()]))
        assert not fetch.fetch_alphafold_model("P04637", cache_dir=tmp_path).from_cache
        assert fetch.fetch_alphafold_model("P04637", cache_dir=tmp_path).from_cache

    def test_refresh_forces_a_refetch(self, transport: InstallTransport, tmp_path: Path) -> None:
        calls = transport(afdb_routes([afdb_entry()]))
        fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        fetch.fetch_alphafold("P04637", cache_dir=tmp_path, refresh=True)
        assert len(calls.calls) == 4

    def test_a_version_bump_is_not_served_from_the_old_cache(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        """The cache key is the remote filename, which carries the version.

        This is the property that makes the whole model_v4 -> v6 class of failure
        impossible: a stale v4 file cannot answer a v6 request.
        """
        transport(afdb_routes([afdb_entry(version=6)]))
        old = fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        transport(afdb_routes([afdb_entry(version=7)]))
        new = fetch.fetch_alphafold("P04637", cache_dir=tmp_path, refresh=True)
        assert old != new
        assert old.exists() and new.exists()

    def test_corrupt_metadata_cache_is_refetched_and_noted(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport(afdb_routes([afdb_entry()]))
        fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        (tmp_path / "alphafold" / "P04637.prediction.json").write_text("{tru", encoding="utf-8")
        model = fetch.fetch_alphafold_model("P04637", cache_dir=tmp_path)
        assert any("unreadable" in note for note in model.notes)

    @pytest.mark.parametrize(
        "payload",
        [
            pytest.param({"entryId": "AF-P04637-F1"}, id="dict-not-list"),
            pytest.param(["AF-P04637-F1"], id="list-of-strings"),
            pytest.param([], id="empty-list"),
            pytest.param("AF-P04637-F1", id="bare-string"),
        ],
    )
    def test_wrongly_shaped_metadata_cache_is_refetched_and_noted(
        self, transport: InstallTransport, tmp_path: Path, payload: Any
    ) -> None:
        """Valid JSON of the wrong shape is still a corrupt cache, and must say so.

        Self-healing is right; doing it silently is not. This file was written by DODO,
        so a shape DODO would never write is evidence of something worth reporting --
        a partial write, a disk fault, or a hand-edited cache.
        """
        transport(afdb_routes([afdb_entry()]))
        fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        metadata = tmp_path / "alphafold" / "P04637.prediction.json"
        metadata.write_text(json.dumps(payload), encoding="utf-8")
        model = fetch.fetch_alphafold_model("P04637", cache_dir=tmp_path)
        assert model.path.is_file()
        assert any("Cached AFDB metadata" in note for note in model.notes), model.notes

    def test_a_partly_unusable_metadata_cache_is_noted(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        """A list whose entries are not all prediction objects: usable, but not intact."""
        transport(afdb_routes([afdb_entry()]))
        fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        metadata = tmp_path / "alphafold" / "P04637.prediction.json"
        entries = json.loads(metadata.read_text(encoding="utf-8"))
        metadata.write_text(json.dumps([*entries, "junk"]), encoding="utf-8")
        model = fetch.fetch_alphafold_model("P04637", cache_dir=tmp_path)
        assert any("Cached AFDB metadata" in note for note in model.notes), model.notes

    def test_a_well_shaped_metadata_cache_produces_no_note(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport(afdb_routes([afdb_entry()]))
        fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert fetch.fetch_alphafold_model("P04637", cache_dir=tmp_path).notes == ()

    def test_missing_model_file_is_redownloaded_from_cached_metadata(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        calls = transport(afdb_routes([afdb_entry()]))
        path = fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        path.unlink()
        fetch.fetch_alphafold("P04637", cache_dir=tmp_path)
        assert path.exists()
        assert len(calls.calls) == 3, "only the file, not the metadata, needed refetching"


# ---------------------------------------------------------------------------
# RCSB
# ---------------------------------------------------------------------------


class TestFetchPdb:
    def test_defaults_to_mmcif(self, transport: InstallTransport, tmp_path: Path) -> None:
        calls = transport({"files.rcsb.org": body_route(b"data_6KN7\n")})
        path = fetch.fetch_pdb("6kn7", cache_dir=tmp_path)
        assert path.name == "6KN7.cif"
        assert calls.urls == ["https://files.rcsb.org/download/6KN7.cif"]

    def test_pdb_format(self, transport: InstallTransport, tmp_path: Path) -> None:
        calls = transport({"files.rcsb.org": body_route(b"HEADER\n")})
        path = fetch.fetch_pdb("1UBQ", format="pdb", cache_dir=tmp_path)
        assert path.name == "1UBQ.pdb"
        assert calls.urls == ["https://files.rcsb.org/download/1UBQ.pdb"]

    def test_mmcif_alias(self, transport: InstallTransport, tmp_path: Path) -> None:
        transport({"files.rcsb.org": body_route(b"data_1UBQ\n")})
        assert fetch.fetch_pdb("1UBQ", format="mmCIF", cache_dir=tmp_path).suffix == ".cif"

    def test_unsupported_format(self, tmp_path: Path) -> None:
        with pytest.raises(UnsupportedFormatError, match="format must be one of"):
            fetch.fetch_pdb("1UBQ", format="xyz", cache_dir=tmp_path)

    def test_malformed_id(self, transport: InstallTransport, tmp_path: Path) -> None:
        calls = transport({})
        with pytest.raises(RetrievalError, match="not a PDB ID"):
            fetch.fetch_pdb("ubiquitin", cache_dir=tmp_path)
        assert calls.calls == []

    def test_404_suggests_cif_when_pdb_was_asked_for(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        """Large EM assemblies -- DODO's usual input -- have no PDB-format file at all."""
        transport({"files.rcsb.org": status_route(404, b"<html>404</html>")})
        with pytest.raises(StructureNotFoundError) as excinfo:
            fetch.fetch_pdb("8J07", format="pdb", cache_dir=tmp_path)
        assert 'format="cif"' in str(excinfo.value)

    def test_404_for_cif_does_not_suggest_cif(
        self, transport: InstallTransport, tmp_path: Path
    ) -> None:
        transport({"files.rcsb.org": status_route(404)})
        with pytest.raises(StructureNotFoundError) as excinfo:
            fetch.fetch_pdb("ZZZZ".replace("Z", "9", 1), cache_dir=tmp_path)
        assert 'format="cif"' not in str(excinfo.value)

    def test_cache_hit(self, transport: InstallTransport, tmp_path: Path) -> None:
        calls = transport({"files.rcsb.org": body_route(b"data_1UBQ\n")})
        first = fetch.fetch_pdb("1UBQ", cache_dir=tmp_path)
        second = fetch.fetch_pdb("1ubq", cache_dir=tmp_path)
        assert first == second
        assert len(calls.calls) == 1

    def test_extended_id(self, transport: InstallTransport, tmp_path: Path) -> None:
        transport({"files.rcsb.org": body_route(b"data_x\n")})
        assert fetch.fetch_pdb("pdb_00006kn7", cache_dir=tmp_path).name == "PDB_00006KN7.cif"


# ---------------------------------------------------------------------------
# UniProt name resolution
# ---------------------------------------------------------------------------


def search_payload(accession: str) -> dict[str, Any]:
    return {"results": [{"primaryAccession": accession}]}


@pytest.fixture
def no_getsequence(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    """Make ``import getSequence`` fail, whether or not it is installed."""
    monkeypatch.setitem(sys.modules, "getSequence", None)
    yield


@pytest.fixture
def fake_getsequence(monkeypatch: pytest.MonkeyPatch) -> Iterator[list[str]]:
    """Install a stand-in ``getSequence`` and record the names it is asked about."""
    seen: list[str] = []
    module = types.ModuleType("getSequence")

    def getseq(name: str) -> list[str]:
        seen.append(name)
        return ["sp|P04637|P53_HUMAN Cellular tumor antigen p53", "MEEPQSDPSV"]

    module.getseq = getseq  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "getSequence", module)
    yield seen


class TestResolveUniprotAccession:
    def test_accession_short_circuits_without_a_request(self, transport: InstallTransport) -> None:
        calls = transport({})
        assert fetch.resolve_uniprot_accession("p04637") == "P04637"
        assert calls.calls == []

    def test_isoform_accession_short_circuits(self, transport: InstallTransport) -> None:
        transport({})
        assert fetch.resolve_uniprot_accession("P04637-2") == "P04637-2"

    def test_rest_search_is_the_preferred_path(
        self, transport: InstallTransport, no_getsequence: None
    ) -> None:
        calls = transport(
            {"rest.uniprot.org/uniprotkb/search": json_route(search_payload("P04637"))}
        )
        assert fetch.resolve_uniprot_accession("human p53") == "P04637"
        assert len(calls.calls) == 1
        assert "query=human+p53" in calls.urls[0]
        assert "fields=accession" in calls.urls[0]

    def test_search_never_asks_for_a_single_row(
        self, transport: InstallTransport, no_getsequence: None
    ) -> None:
        """size=1 is a trap: the live API does not relevance-rank a one-row page.

        Measured: 'human p53' resolves to Q89743 at size=1 and to P04637 at size=2, so
        a one-row request silently resolves to the wrong protein.
        """
        calls = transport({"search": json_route(search_payload("P04637"))})
        fetch.resolve_uniprot_accession("human p53")
        assert "size=1&" not in calls.urls[0] and not calls.urls[0].endswith("size=1")

    def test_top_of_a_multi_row_page_is_taken(
        self, transport: InstallTransport, no_getsequence: None
    ) -> None:
        payload = {"results": [{"primaryAccession": acc} for acc in ("P04637", "P02340")]}
        transport({"search": json_route(payload)})
        assert fetch.resolve_uniprot_accession("human p53") == "P04637"

    def test_getsequence_is_not_used_when_the_search_succeeds(
        self, transport: InstallTransport, fake_getsequence: list[str]
    ) -> None:
        transport({"search": json_route(search_payload("P04637"))})
        assert fetch.resolve_uniprot_accession("human p53") == "P04637"
        assert fake_getsequence == [], "the stdlib path must take precedence"

    def test_falls_back_to_getsequence_when_the_search_finds_nothing(
        self, transport: InstallTransport, fake_getsequence: list[str]
    ) -> None:
        transport({"search": json_route({"results": []})})
        with pytest.warns(UserWarning, match="getSequence"):
            assert fetch.resolve_uniprot_accession("tumour protein 53") == "P04637"
        assert fake_getsequence == ["tumour protein 53"]

    def test_missing_getsequence_names_the_extra(
        self, transport: InstallTransport, no_getsequence: None
    ) -> None:
        transport({"search": json_route({"results": []})})
        with pytest.raises(MissingDependencyError) as excinfo:
            fetch.resolve_uniprot_accession("something obscure")
        message = str(excinfo.value)
        assert "dodo[lookup]" in message
        assert "built-in UniProt search returned no match" in message

    def test_network_failure_is_not_blamed_on_a_missing_package(
        self, transport: InstallTransport, no_getsequence: None
    ) -> None:
        """A dead connection must not send the user off to install torch."""
        transport({"search": raising_route(TimeoutError("timed out"))})
        with pytest.raises(RetrievalError) as excinfo:
            fetch.resolve_uniprot_accession("human p53")
        assert not isinstance(excinfo.value, MissingDependencyError)
        assert "transient network failure" in str(excinfo.value)

    def test_getsequence_covers_a_network_failure(
        self, transport: InstallTransport, fake_getsequence: list[str]
    ) -> None:
        transport({"search": raising_route(TimeoutError("timed out"))})
        with pytest.warns(UserWarning, match="getSequence"):
            assert fetch.resolve_uniprot_accession("human p53") == "P04637"

    def test_no_hits_and_a_useless_getsequence_reports_no_match(
        self, transport: InstallTransport, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        module = types.ModuleType("getSequence")
        module.getseq = lambda name: ["no accession here", ""]  # type: ignore[attr-defined]
        monkeypatch.setitem(sys.modules, "getSequence", module)
        transport({"search": json_route({"results": []})})
        with pytest.raises(RetrievalError, match="No UniProt entry matches"):
            fetch.resolve_uniprot_accession("gibberish protein")

    def test_getsequence_exception_does_not_escape(
        self, transport: InstallTransport, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        def boom(name: str) -> list[str]:
            raise RuntimeError("torch exploded")

        module = types.ModuleType("getSequence")
        module.getseq = boom  # type: ignore[attr-defined]
        monkeypatch.setitem(sys.modules, "getSequence", module)
        transport({"search": json_route({"results": []})})
        with pytest.raises(RetrievalError, match="No UniProt entry matches"):
            fetch.resolve_uniprot_accession("gibberish protein")

    def test_malformed_search_payload_is_a_retrieval_error_not_an_attributeerror(
        self, transport: InstallTransport, no_getsequence: None
    ) -> None:
        """``results`` full of strings must not crash out as a bare AttributeError.

        Every other malformed-payload branch in this module is defended; a caller
        wrapping ``except RetrievalError`` around this gets an opaque crash otherwise.
        """
        transport({"search": json_route({"results": ["P04637"]})})
        with pytest.raises(RetrievalError) as excinfo:
            fetch.resolve_uniprot_accession("human p53")
        assert not isinstance(excinfo.value, AttributeError)
        assert "unrecognized" in str(excinfo.value)

    def test_search_result_without_an_accession_is_reported(
        self, transport: InstallTransport, no_getsequence: None
    ) -> None:
        """A hit DODO cannot read is not the same thing as "no hits"."""
        transport({"search": json_route({"results": [{"uniProtkbId": "P53_HUMAN"}]})})
        with pytest.raises(RetrievalError, match=r"unrecognized|no primaryAccession"):
            fetch.resolve_uniprot_accession("human p53")

    def test_search_results_of_a_non_list_type_is_reported(
        self, transport: InstallTransport, no_getsequence: None
    ) -> None:
        transport({"search": json_route({"results": {"primaryAccession": "P04637"}})})
        with pytest.raises(RetrievalError, match="unrecognized"):
            fetch.resolve_uniprot_accession("human p53")

    def test_empty_name_rejected(self) -> None:
        with pytest.raises(RetrievalError, match="empty protein name"):
            fetch.resolve_uniprot_accession("   ")

    def test_search_400_is_reported(
        self, transport: InstallTransport, no_getsequence: None
    ) -> None:
        transport({"search": status_route(400)})
        with pytest.raises(RetrievalError, match="malformed"):
            fetch.resolve_uniprot_accession("bad ) query")


# ---------------------------------------------------------------------------
# The getSequence fallback: a third party doing its own networking
# ---------------------------------------------------------------------------


@pytest.fixture
def hanging_getsequence(monkeypatch: pytest.MonkeyPatch) -> Iterator[threading.Event]:
    """Install a ``getSequence`` whose ``getseq`` blocks, like a hung connection.

    The real package makes its own ``requests`` calls with a hardcoded timeout that its
    public API does not expose, so from DODO's side it is a call of unbounded duration.
    The event is set on teardown so the stand-in does not outlive the test.
    """
    release = threading.Event()
    module = types.ModuleType("getSequence")

    def getseq(name: str) -> list[str]:
        release.wait(5.0)
        return ["sp|P04637|P53_HUMAN Cellular tumor antigen p53", "MEEPQSDPSV"]

    module.getseq = getseq  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "getSequence", module)
    yield release
    release.set()


class TestGetSequenceFallbackTimeout:
    def test_a_hanging_getsequence_does_not_hang_dodo(
        self, transport: InstallTransport, hanging_getsequence: threading.Event
    ) -> None:
        """``timeout`` is documented as applying to this call, so it must.

        v1's defect was an unbounded network call, and this path is the last one where
        it survived: the fallback made its own requests and honoured nothing DODO asked
        for.
        """
        transport({"search": json_route({"results": []})})
        started = time.monotonic()
        with pytest.raises(RetrievalError) as excinfo:
            fetch.resolve_uniprot_accession("human p53", timeout=0.05)
        elapsed = time.monotonic() - started
        assert elapsed < 2.0, f"waited {elapsed:.2f} s on a 0.05 s timeout"
        assert "getSequence" in str(excinfo.value)
        assert "0.05" in str(excinfo.value)

    def test_a_hanging_getsequence_does_not_mask_a_network_failure(
        self, transport: InstallTransport, hanging_getsequence: threading.Event
    ) -> None:
        """If the stdlib search failed on the network, say that; it is the real cause."""
        transport({"search": raising_route(TimeoutError("timed out"))})
        with pytest.raises(RetrievalError) as excinfo:
            fetch.resolve_uniprot_accession("human p53", timeout=0.05)
        assert "transient network failure" in str(excinfo.value)

    def test_the_fallback_result_is_flagged_as_unverified(
        self, transport: InstallTransport, fake_getsequence: list[str]
    ) -> None:
        """An answer from a third party that DODO cannot audit is not a silent answer.

        getSequence resolves free text with its own ranking, off DODO's allowlisted
        hosts, and will return a plausible accession for a query that has no real match
        (measured: 'some obscure protein' -> Q7S936). The caller gets told.
        """
        transport({"search": json_route({"results": []})})
        with pytest.warns(UserWarning) as warned:
            assert fetch.resolve_uniprot_accession("some obscure protein") == "P04637"
        message = str(warned[0].message)
        assert "getSequence" in message
        assert "P04637" in message

    def test_the_stdlib_path_warns_about_nothing(
        self, transport: InstallTransport, fake_getsequence: list[str]
    ) -> None:
        import warnings as warnings_module

        transport({"search": json_route(search_payload("P04637"))})
        with warnings_module.catch_warnings():
            warnings_module.simplefilter("error")
            assert fetch.resolve_uniprot_accession("human p53") == "P04637"


# ---------------------------------------------------------------------------
# UniProt sequence
# ---------------------------------------------------------------------------

FASTA = (
    b">sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53\n"
    b"MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP\n"
    b"DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGF\n"
)


class TestFetchUniprotSequence:
    def test_parses_fasta(self, transport: InstallTransport) -> None:
        calls = transport({"rest.uniprot.org": body_route(FASTA)})
        sequence = fetch.fetch_uniprot_sequence("p04637")
        assert sequence.startswith("MEEPQSDPSV")
        assert sequence.endswith("YQGSYGF")
        assert "\n" not in sequence and ">" not in sequence
        assert calls.urls == ["https://rest.uniprot.org/uniprotkb/P04637.fasta"]

    def test_404_is_not_found(self, transport: InstallTransport) -> None:
        transport({"rest.uniprot.org": status_route(404)})
        with pytest.raises(StructureNotFoundError, match="no entry"):
            fetch.fetch_uniprot_sequence("P04637")

    def test_400_is_a_retrieval_error(self, transport: InstallTransport) -> None:
        transport({"rest.uniprot.org": status_route(400)})
        with pytest.raises(RetrievalError) as excinfo:
            fetch.fetch_uniprot_sequence("P04637")
        assert not isinstance(excinfo.value, StructureNotFoundError)

    def test_malformed_accession_rejected_without_a_request(
        self, transport: InstallTransport
    ) -> None:
        calls = transport({})
        with pytest.raises(RetrievalError, match="not a UniProt accession"):
            fetch.fetch_uniprot_sequence("human p53")
        assert calls.calls == []

    def test_non_fasta_body_is_reported(self, transport: InstallTransport) -> None:
        transport({"rest.uniprot.org": body_route(b"<html>oops</html>")})
        with pytest.raises(RetrievalError, match="Expected FASTA"):
            fetch.fetch_uniprot_sequence("P04637")

    def test_header_only_body_is_reported(self, transport: InstallTransport) -> None:
        transport({"rest.uniprot.org": body_route(b">sp|P04637|P53_HUMAN\n\n")})
        with pytest.raises(RetrievalError, match="no sequence in it"):
            fetch.fetch_uniprot_sequence("P04637")

    def test_non_letter_characters_rejected(self, transport: InstallTransport) -> None:
        transport({"rest.uniprot.org": body_route(b">sp|P04637|P53\nMEEP*QSD\n")})
        with pytest.raises(RetrievalError, match="non-letter characters"):
            fetch.fetch_uniprot_sequence("P04637")

    def test_timeout_is_passed_through(self, transport: InstallTransport) -> None:
        calls = transport({"rest.uniprot.org": body_route(FASTA)})
        fetch.fetch_uniprot_sequence("P04637", timeout=7.5)
        assert [call.timeout for call in calls.calls] == [7.5]


# ---------------------------------------------------------------------------
# Live services. Marked so CI can deselect them: -m 'not network'
# ---------------------------------------------------------------------------


@pytest.mark.network
@pytest.mark.slow
class TestLive:
    """Live checks against EBI, RCSB and UniProt.

    Marked slow as well as network: individual cases were measured at 0.2-1.2 s on a
    good connection, and any of them can cross a second on a worse one.
    """

    def test_alphafold_p53(self, tmp_path: Path) -> None:
        model = fetch.fetch_alphafold_model("P04637", cache_dir=tmp_path)
        # Asserting on the *shape* of the URL, not the version: the version is exactly
        # what this module refuses to assume, and pinning v6 here would re-create the
        # bug at the next release.
        assert model.model_url.startswith("https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v")
        assert model.model_url.endswith(".pdb")
        assert not model.is_fragment
        assert model.mean_plddt is not None and 0.0 < model.mean_plddt <= 100.0
        assert model.path.is_file() and model.path.stat().st_size > 1000
        assert model.path.read_text(encoding="utf-8").startswith(("HEADER", "ATOM", "REMARK"))

    def test_alphafold_p300(self, tmp_path: Path) -> None:
        """Q09472 is 2414 residues -- just under the fragmentation limit."""
        model = fetch.fetch_alphafold_model("Q09472", cache_dir=tmp_path)
        assert model.entry_id == "AF-Q09472-F1"
        assert model.residue_range == (1, 2414)
        assert not model.is_fragment

    def test_titin_has_no_model(self, tmp_path: Path) -> None:
        with pytest.raises(StructureNotFoundError, match="Retrying will not help"):
            fetch.fetch_alphafold("Q8WZ42", cache_dir=tmp_path)

    def test_malformed_accession_at_the_api(self, tmp_path: Path) -> None:
        # Valid syntax locally, so this reaches EBI and comes back 400 or 404: either
        # way it must not be mistaken for a transient failure.
        with pytest.raises((RetrievalError, StructureNotFoundError)):
            fetch.fetch_alphafold("Q0Q0Q1", cache_dir=tmp_path)

    def test_live_cache_hit(self, tmp_path: Path) -> None:
        first = fetch.fetch_alphafold_model("P04637", cache_dir=tmp_path)
        second = fetch.fetch_alphafold_model("P04637", cache_dir=tmp_path)
        assert second.from_cache and not first.from_cache

    def test_resolve_name(self) -> None:
        assert fetch.resolve_uniprot_accession("human p53") == "P04637"
        assert fetch.resolve_uniprot_accession("p300 human") == "Q09472"

    def test_uniprot_sequence(self) -> None:
        sequence = fetch.fetch_uniprot_sequence("P04637")
        assert len(sequence) == 393
        assert sequence.startswith("MEEPQSDPSV")

    def test_rcsb_cif(self, tmp_path: Path) -> None:
        path = fetch.fetch_pdb("1UBQ", cache_dir=tmp_path)
        assert path.read_text(encoding="utf-8").startswith("data_1UBQ")

    def test_rcsb_missing_entry(self, tmp_path: Path) -> None:
        with pytest.raises(StructureNotFoundError):
            fetch.fetch_pdb("9ZZZ", cache_dir=tmp_path)
