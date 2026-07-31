"""Command-line interface.

One command with subcommands, replacing v1's three separate console scripts
(``pdb-from-name``, ``pdb-from-pdb``, ``pdb-from-sequence``). Those were near-duplicate
argparse files whose shared flags had drifted apart -- ``-apr`` meant
``--attempts_per_region`` in two of them and ``--attempts_per_residue`` in the third, and the
README documented several defaults that disagreed with the code.

Startup cost is deliberate: nothing heavy is imported at module scope, so ``dodo --help``
does not pay for torch. The scientific dependencies load only when a subcommand that needs
them actually runs.
"""

from __future__ import annotations

import argparse
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Sequence

#: Modes, duplicated here as strings so ``--help`` needs no heavy import. Validated for real
#: against :data:`dodo.constants.MODES` once a command runs.
_MODE_CHOICES = (
    "super_compact",
    "compact",
    "normal",
    "predicted",
    "expanded",
    "super_expanded",
    "max_expansion",
)
_STRATEGY_CHOICES = ("auto", "contact", "plddt", "metapredict")
_ENGINE_CHOICES = ("walk", "starling")


def _add_common_build_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the flags shared by every building subcommand, defined once."""
    parser.add_argument("-o", "--out", required=True, help="output PDB path")
    parser.add_argument(
        "-m",
        "--mode",
        default="predicted",
        choices=_MODE_CHOICES,
        help=(
            "target dimension as a multiplier on the ALBATROSS-predicted end-to-end "
            "distance (default: predicted, i.e. 1.0x)"
        ),
    )
    parser.add_argument(
        "-n",
        "--models",
        type=int,
        default=1,
        metavar="N",
        help=(
            "number of conformers. Folded domains stay put; each model draws its own IDR "
            "dimensions from the physical distribution (default: 1)"
        ),
    )
    parser.add_argument(
        "-e",
        "--engine",
        default="walk",
        choices=_ENGINE_CHOICES,
        help="conformation engine (default: walk; starling needs pip install 'dodo[starling]')",
    )
    parser.add_argument(
        "--all-atom",
        action="store_true",
        help="place N, C and O for rebuilt regions (see README for current limitations)",
    )
    parser.add_argument(
        "--sidechains",
        action="store_true",
        help="also place CB; only meaningful with --all-atom",
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="random seed; makes output reproducible"
    )
    parser.add_argument(
        "--no-conect",
        action="store_true",
        help=(
            "omit CONECT records. Not recommended: CA-CA spacing exceeds the auto-bond "
            "cutoff in VMD and PyMOL, so without them a rebuilt region renders as "
            "disconnected dots"
        ),
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="suppress the per-region report")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dodo",
        description=(
            "Rebuild the disordered regions of a predicted protein structure so they adopt "
            "realistic polymer dimensions."
        ),
    )
    parser.add_argument("--version", action="store_true", help="print the version and exit")
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND")

    rebuild_parser = subparsers.add_parser(
        "rebuild",
        help="rebuild the IDRs of a local PDB or mmCIF file",
        description="Rebuild the disordered regions of a structure file.",
    )
    rebuild_parser.add_argument("structure", help="path to a PDB or mmCIF file (.gz accepted)")
    rebuild_parser.add_argument(
        "-s",
        "--strategy",
        default="auto",
        choices=_STRATEGY_CHOICES,
        help=(
            "how to identify regions: auto picks pLDDT for AlphaFold models and geometric "
            "contacts otherwise (default: auto)"
        ),
    )
    _add_common_build_arguments(rebuild_parser)

    fetch_parser = subparsers.add_parser(
        "fetch",
        help="download a structure from the AlphaFold database and rebuild it",
        description=(
            "Resolve a UniProt accession or protein name, download its AlphaFold model, and "
            "rebuild its disordered regions."
        ),
    )
    fetch_parser.add_argument(
        "target", help="UniProt accession (e.g. P04637) or a protein name (e.g. 'human p53')"
    )
    fetch_parser.add_argument(
        "-s", "--strategy", default="auto", choices=_STRATEGY_CHOICES, help="(default: auto)"
    )
    _add_common_build_arguments(fetch_parser)

    sequence_parser = subparsers.add_parser(
        "sequence",
        help="build coordinates for a disordered sequence, with no input structure",
        description="Build a disordered region from sequence alone.",
    )
    sequence_parser.add_argument("sequence", help="one-letter amino acid sequence")
    _add_common_build_arguments(sequence_parser)

    regions_parser = subparsers.add_parser(
        "regions",
        help="report the folded domains and IDRs of a structure without rebuilding it",
        description=(
            "Identify and print regions. Useful for checking what DODO thinks your structure "
            "looks like before spending time rebuilding it."
        ),
    )
    regions_parser.add_argument("structure", help="path to a PDB or mmCIF file")
    regions_parser.add_argument(
        "-s", "--strategy", default="auto", choices=_STRATEGY_CHOICES, help="(default: auto)"
    )
    regions_parser.add_argument(
        "--scores", action="store_true", help="also print the per-residue score profile"
    )

    return parser


def _report(report: object, *, quiet: bool) -> int:
    """Print a rebuild report and return the process exit status."""
    from .construct.pipeline import RebuildReport

    assert isinstance(report, RebuildReport)
    if not quiet:
        print(report.summary(), file=sys.stderr)
    if not report.models:
        print("dodo: nothing could be built", file=sys.stderr)
        return 1
    # A partial success is not a failure, but it must not look like a clean run either.
    return 0 if report.ok else 2


def main(argv: Sequence[str] | None = None) -> int:
    """Entry point. Returns a process exit status rather than calling sys.exit."""
    parser = _build_parser()
    args = parser.parse_args(argv)

    if getattr(args, "version", False):
        from . import __version__

        print(__version__)
        return 0

    if args.command is None:
        parser.print_help()
        return 0

    from .exceptions import DodoError

    try:
        if args.command == "regions":
            from .io import read_structure
            from .regions import assign_regions

            structure = read_structure(args.structure)
            for note in structure.notes:
                print(f"note: {note}", file=sys.stderr)
            for assignment in assign_regions(structure, strategy=args.strategy):
                print(assignment.describe())
                for note in assignment.notes:
                    print(f"  note: {note}", file=sys.stderr)
                if args.scores:
                    for index, value in enumerate(assignment.score):
                        print(f"  {index + 1}\t{value:.3f}")
            return 0

        from .construct.pipeline import build_from_sequence, rebuild
        from .io import write_pdb

        common = {
            "mode": args.mode,
            "n_models": args.models,
            "engine": args.engine,
            "all_atom": args.all_atom,
            "sidechains": args.sidechains,
            "seed": args.seed,
        }

        if args.command == "sequence":
            report = build_from_sequence(args.sequence, **common)
        elif args.command == "rebuild":
            report = rebuild(args.structure, strategy=args.strategy, **common)
        elif args.command == "fetch":
            from .io import fetch_alphafold, resolve_uniprot_accession

            accession = args.target
            # A bare accession has no spaces and is short; anything else is a name to resolve.
            if " " in accession or len(accession) > 10:
                accession = resolve_uniprot_accession(args.target)
                print(f"resolved {args.target!r} to {accession}", file=sys.stderr)
            path = fetch_alphafold(accession)
            print(f"using {path}", file=sys.stderr)
            report = rebuild(path, strategy=args.strategy, **common)
        else:  # pragma: no cover - argparse restricts the choices
            parser.error(f"unknown command {args.command!r}")

        status = _report(report, quiet=args.quiet)
        if report.models:
            write_pdb(report.models, args.out, conect=not args.no_conect)
            if not args.quiet:
                print(f"wrote {args.out}", file=sys.stderr)
        return status

    except DodoError as exc:
        # A DODO error is an expected, explained failure. Print it plainly rather than
        # dumping a traceback at a user who cannot act on it.
        print(f"dodo: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:  # pragma: no cover
        print("dodo: interrupted", file=sys.stderr)
        return 130


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
