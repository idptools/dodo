"""Transport-level download tests against a real socket and real urllib.

These live apart from ``test_io_fetch.py`` because they cannot be written with a mock. The
defect that motivated them was *introduced by a fix* and was invisible to the existing
test suite for a structural reason: the fakes there model a response as a ``BytesIO`` with
a headers dict, which can express "has a Content-Length" and "has none", but cannot
express ``Transfer-Encoding: chunked``. Chunked framing is a property of ``http.client``'s
parser, not of the object DODO sees, so only a real HTTP exchange exercises it.

The history worth preserving:

1. ``_download_to`` used ``shutil.copyfileobj`` and never compared bytes written against
   the declared length. ``http.client`` deliberately tolerates a short identity-length
   body, so a truncated download completed "successfully" and was committed to the cache
   as permanent silent corruption.
2. The fix added a ``Content-Length`` check -- and read the header unconditionally. But
   RFC 9112 section 6.1 says chunked framing takes precedence and ``Content-Length`` must
   be ignored, which ``http.client`` implements. So a *complete* chunked download was
   rejected as truncated whenever an intermediary added a ``Content-Length``. Every
   endpoint DODO downloads from serves chunked, so a corporate proxy or CDN would have
   broken all structure downloads with a misleading "bytes are missing" error.

Both directions are asserted here, so neither can regress into the other.
"""

from __future__ import annotations

import socket
import threading
import urllib.request
from pathlib import Path
from typing import ClassVar

import pytest

from dodo.exceptions import RetrievalError
from dodo.io import fetch

#: A stand-in body. Contents are irrelevant; only its length is under test.
BODY = b"ATOM  " + b"x" * 5000

#: An allowlisted URL. The probe server is local, but ``_download_to`` checks the URL it
#: is given against DODO's host allowlist before opening anything, and that check should
#: stay exercised rather than be bypassed by passing a localhost URL.
ALLOWED_URL = "https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v6.pdb"


def _chunked_frames(body: bytes) -> list[bytes]:
    """Encode ``body`` as a single well-formed HTTP chunk plus the terminator."""
    return [f"{len(body):x}\r\n".encode(), body, b"\r\n", b"0\r\n\r\n"]


def _serve_once(headers: bytes, frames: list[bytes]) -> int:
    """Serve exactly one HTTP response on an ephemeral port, then close. Returns the port."""
    server = socket.socket()
    server.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    server.bind(("127.0.0.1", 0))
    server.listen(1)
    port = int(server.getsockname()[1])

    def run() -> None:
        try:
            conn, _ = server.accept()
            conn.recv(65535)
            conn.sendall(b"HTTP/1.1 200 OK\r\n" + headers + b"\r\n")
            for frame in frames:
                conn.sendall(frame)
            conn.close()
        except OSError:
            pass
        finally:
            server.close()

    threading.Thread(target=run, daemon=True).start()
    return port


def _download(port: int, destination: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Run ``_download_to`` against the local probe server using real urllib."""
    local = f"http://127.0.0.1:{port}/probe.pdb"
    monkeypatch.setattr(
        fetch,
        "_urlopen",
        lambda url, timeout: urllib.request.urlopen(local, timeout=timeout),
    )
    fetch._download_to(ALLOWED_URL, destination, timeout=10.0)


class TestIdentityLengthVerification:
    """A declared Content-Length must be enforced."""

    def test_truncated_body_is_rejected(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        port = _serve_once(f"Content-Length: {len(BODY)}\r\n".encode(), [BODY[:100]])
        destination = tmp_path / "out.pdb"
        with pytest.raises(RetrievalError, match=r"incomplete|missing"):
            _download(port, destination, monkeypatch)

    def test_truncated_body_is_not_committed_to_the_cache(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The corruption that mattered was permanent, not transient."""
        port = _serve_once(f"Content-Length: {len(BODY)}\r\n".encode(), [BODY[:100]])
        destination = tmp_path / "out.pdb"
        with pytest.raises(RetrievalError):
            _download(port, destination, monkeypatch)
        assert not destination.exists(), "a truncated download must leave nothing behind"
        assert not list(tmp_path.iterdir()), "no partial file may survive either"

    def test_complete_body_is_accepted(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        port = _serve_once(f"Content-Length: {len(BODY)}\r\n".encode(), [BODY])
        destination = tmp_path / "out.pdb"
        _download(port, destination, monkeypatch)
        assert destination.read_bytes() == BODY


class TestChunkedFraming:
    """Chunked framing takes precedence over Content-Length (RFC 9112 section 6.1)."""

    def test_complete_chunked_body_is_accepted(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        port = _serve_once(b"Transfer-Encoding: chunked\r\n", _chunked_frames(BODY))
        destination = tmp_path / "out.pdb"
        _download(port, destination, monkeypatch)
        assert destination.read_bytes() == BODY

    @pytest.mark.parametrize("declared", [99999, 10, 0])
    def test_chunked_body_ignores_a_contradictory_content_length(
        self, declared: int, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A complete chunked download must not be rejected over a stale Content-Length.

        This is the regression. An intermediary that adds a Content-Length to a chunked
        response would otherwise break every structure download -- and all of AFDB's and
        RCSB's endpoints serve chunked, so it would break all of them at once. Both a
        too-large and a too-small declared value are tested, since the naive check had a
        separate failure branch for each.
        """
        port = _serve_once(
            b"Transfer-Encoding: chunked\r\nContent-Length: %d\r\n" % declared,
            _chunked_frames(BODY),
        )
        destination = tmp_path / "out.pdb"
        _download(port, destination, monkeypatch)
        assert destination.read_bytes() == BODY

    def test_truncated_chunked_body_is_still_caught(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Ignoring Content-Length must not mean giving up on chunked truncation.

        A chunk header promising more than arrives is a framing violation, so
        ``http.client`` raises ``IncompleteRead`` -- which is an ``HTTPException``, not an
        ``OSError``, and so needs its own handler to become a ``RetrievalError``.
        """
        frames = [f"{len(BODY):x}\r\n".encode(), BODY[:100]]  # promises more, then hangs up
        port = _serve_once(b"Transfer-Encoding: chunked\r\n", frames)
        destination = tmp_path / "out.pdb"
        with pytest.raises(RetrievalError):
            _download(port, destination, monkeypatch)
        assert not destination.exists()


class TestContentLengthHelper:
    """Unit-level checks on the framing decision itself."""

    def test_chunked_response_declares_nothing_to_verify(self) -> None:
        class Chunked:
            chunked = True
            headers: ClassVar[dict[str, str]] = {"Content-Length": "99999"}

        assert fetch._content_length(Chunked()) is None

    def test_identity_response_declares_its_length(self) -> None:
        class Identity:
            chunked = False
            headers: ClassVar[dict[str, str]] = {"Content-Length": "5006"}

        assert fetch._content_length(Identity()) == 5006

    @pytest.mark.parametrize("raw", ["", "abc", "12.5", "-1", None])
    def test_unusable_header_declares_nothing(self, raw: str | None) -> None:
        class Mangled:
            chunked = False
            headers = {"Content-Length": raw} if raw is not None else {}

        assert fetch._content_length(Mangled()) is None


class TestCacheGeneration:
    """The cache directory carries a generation so untrustworthy entries are abandoned."""

    def test_cache_dir_includes_the_generation(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv("DODO_CACHE_DIR", "/tmp/dodo-test-cache")
        assert fetch.default_cache_dir().name == fetch.CACHE_GENERATION

    def test_generation_applies_to_the_override_too(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """An override must not be able to resurrect a pre-verification cache layout."""
        monkeypatch.setenv("DODO_CACHE_DIR", "/tmp/dodo-test-cache")
        resolved = fetch.default_cache_dir()
        assert resolved == Path("/tmp/dodo-test-cache") / fetch.CACHE_GENERATION
