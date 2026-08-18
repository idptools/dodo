"""Exception hierarchy for DODO.

Every failure path in DODO raises one of these. That is a deliberate change from the
pre-rewrite code, whose dominant failure mode was returning a plausible-looking wrong
answer: builders returned coordinate arrays full of exact ``(0, 0, 0)`` rows on total
failure, samplers returned NaN, and the folded-domain placer returned positions it had
already determined to be clashing. For a tool whose product is a picture, a silent
wrong answer is worse than a crash -- it becomes a figure in a paper.

The rule this package follows: if a routine cannot do what it was asked, it raises. If
partial success is meaningful, it returns an explicit success mask alongside the
results. It never returns degenerate coordinates and hopes.
"""

from __future__ import annotations


class DodoError(Exception):
    """Base class for every DODO-specific exception.

    Catch this to catch anything DODO raises on purpose.
    """


class InvalidParameterError(DodoError, ValueError):
    """An argument was outside its permitted range or set of values.

    Deliberately inherits from **both** :class:`DodoError` and :exc:`ValueError`. A caller
    migrating from 1.x writes ``except DodoError`` as the translation of 1.x's
    ``dodoException``, and the thing they most often wrapped it around was a mistyped mode name
    -- so bad arguments have to be inside the hierarchy or that translation silently fails to
    catch the commonest mistake. But a mistyped keyword argument is also exactly what
    :exc:`ValueError` means in Python, and code already written as ``except ValueError`` should
    not break to serve the first point.

    Inheriting from both costs nothing and makes either instinct correct. It also means the CLI
    catches these through its single ``except DodoError`` and prints one plain line, instead of
    dumping a traceback at a user who mistyped a flag value.
    """


# ---------------------------------------------------------------------------
# Input and parsing
# ---------------------------------------------------------------------------


class StructureFileError(DodoError):
    """A structure file could not be read, or is not the format it claims to be."""


class UnsupportedFormatError(StructureFileError):
    """The file extension or content is a format DODO does not handle."""


class MalformedRecordError(StructureFileError):
    """A record in an otherwise readable file could not be parsed.

    Carries the offending line so the user can find it, since real PDB files in the
    wild routinely contain one bad line in an otherwise fine structure.
    """

    def __init__(self, message: str, line_number: int | None = None, line: str | None = None):
        self.line_number = line_number
        self.line = line
        if line_number is not None:
            message = f"{message} (line {line_number})"
        if line is not None:
            message = f"{message}\n  {line!r}"
        super().__init__(message)


class EmptyStructureError(StructureFileError):
    """The file parsed successfully but contained no usable atoms."""


# ---------------------------------------------------------------------------
# Network retrieval
# ---------------------------------------------------------------------------


class RetrievalError(DodoError):
    """A remote structure or sequence could not be retrieved."""


class StructureNotFoundError(RetrievalError):
    """The requested accession has no structure available.

    Distinct from a transient network failure: retrying will not help. The AlphaFold
    database genuinely has no model for some UniProt entries, very long proteins
    (titin) among them.
    """


# ---------------------------------------------------------------------------
# Region identification
# ---------------------------------------------------------------------------


class RegionIdentificationError(DodoError):
    """Folded domains, IDRs and loops could not be assigned."""


class InvalidRegionError(RegionIdentificationError):
    """A caller-supplied region definition is inconsistent.

    Raised for overlapping regions, out-of-range indices, inverted spans, and gaps
    that leave residues unassigned. The pre-rewrite manual-override path accepted all
    of these silently and crashed much later with an unrelated message.
    """


# ---------------------------------------------------------------------------
# Building
# ---------------------------------------------------------------------------


class BuildError(DodoError):
    """A structure or region could not be built."""


class GeometryError(BuildError):
    """A geometric operation was asked for something impossible.

    For example: closing a chain of ``n`` residues onto anchors further apart than
    ``(n - 1) * CA_CA_BOND_LENGTH``, or building a region of zero length.
    """


class UnsatisfiableTargetError(BuildError):
    """The requested dimension cannot be achieved by any conformation.

    Typically a target end-to-end distance exceeding the contour length of the
    sequence. Carries both numbers so the caller can report or clamp.
    """

    def __init__(self, message: str, target: float | None = None, achievable: float | None = None):
        self.target = target
        self.achievable = achievable
        if target is not None and achievable is not None:
            message = f"{message} (requested {target:.1f} A, maximum achievable {achievable:.1f} A)"
        super().__init__(message)


class ExhaustedAttemptsError(BuildError):
    """A region could not be built within the allotted number of attempts.

    Distinguishes "ran out of budget" from "geometrically impossible": this one is
    worth retrying with a larger budget or a different seed, and
    :class:`UnsatisfiableTargetError` is not.
    """

    def __init__(self, message: str, attempts: int | None = None):
        self.attempts = attempts
        if attempts is not None:
            message = f"{message} (after {attempts} attempts)"
        super().__init__(message)


# ---------------------------------------------------------------------------
# Optional dependencies
# ---------------------------------------------------------------------------


class MissingDependencyError(DodoError):
    """An optional dependency is required for the requested operation.

    The message always names the install command, because a bare "ModuleNotFoundError"
    tells a user nothing about how to resolve it.
    """

    def __init__(self, package: str, purpose: str, extra: str | None = None):
        self.package = package
        self.extra = extra
        install = f"pip install 'idptools-dodo[{extra}]'" if extra else f"pip install {package}"
        super().__init__(
            f"{purpose} requires {package}, which is not installed.\n  Install it with: {install}"
        )


class EngineError(DodoError):
    """A conformation engine failed or was misused."""


class EngineUnavailableError(EngineError):
    """The requested engine cannot run in this environment.

    Separate from :class:`MissingDependencyError` because an engine can be installed
    and still unavailable -- STARLING with no model weights downloaded, for instance.
    """
