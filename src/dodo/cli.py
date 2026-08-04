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
import os
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
_STRATEGY_CHOICES = ("auto", "density", "contact", "plddt")
_ENGINE_CHOICES = ("walk", "starling")


def _add_common_build_arguments(parser: argparse.ArgumentParser) -> None:
    """Add the flags shared by every building subcommand, defined once."""
    parser.add_argument("-o", "--out", required=True, help="output PDB path")
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help=(
            "do not read or write DODO's per-user cache. Two things are cached: ALBATROSS "
            "predictions (about 116 bytes each, and a hit avoids importing torch, worth 1.5 s) "
            "and downloaded structures (mean 0.25 MB, largest 1.51 MB over 259 measured files). "
            "Both are on by default; this opts out of both, as does DODO_NO_CACHE=1"
        ),
    )
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
        help=(
            "conformation engine (default: walk). 'starling' is UNDOCUMENTED in 2.0 and not "
            "verified end to end: check the exit status, because a failed region keeps its "
            "input coordinates while the folded domains move anyway"
        ),
    )
    # Backbone placement for REBUILT regions is --backbone, below. Side-chain placement is not in
    # 2.0 at all, and there is no flag for it.
    #
    # v1's ``all_atom`` and ``sidechains`` keyword arguments are gone rather than hidden. In the 2.0
    # development tree ``rebuild(all_atom=True)`` was a silent no-op -- byte-identical output --
    # while ``build_from_sequence(all_atom=True)`` raised. A flag that lies is worse than no flag.
    #
    # Note this is NOT the same feature as v1's "all atom" option, which meant *keeping* every
    # atom of the regions DODO does not rebuild. v2 does that unconditionally: folded domains
    # retain their full atomic detail and only rebuilt regions are reduced to alpha carbons.
    parser.add_argument(
        "--seed", type=int, default=None, help="random seed; makes output reproducible"
    )
    # v1 parity. These were `-f/--no_FD_atoms` and `-b/--beta_for_FD_IDR` on pdb-from-pdb and
    # pdb-from-name. Both behaviours already existed in v2's writer but had no command line, which
    # would have forced anyone with a shell script using them to rewrite it in Python.
    parser.add_argument(
        "--backbone",
        action="store_true",
        help=(
            "also place N, C and O on the rebuilt regions, inferred from four consecutive alpha "
            "carbons. Opt-in for now: the geometry inside a rebuilt region is exact, but the "
            "peptide bond where a region meets a folded domain is unsatisfiable and is left long"
        ),
    )
    parser.add_argument(
        "--ca-only",
        action="store_true",
        help=(
            "write alpha carbons only, for the folded domains too (v1's -f/--no_FD_atoms). "
            "By default the regions DODO does not rebuild keep every atom they arrived with"
        ),
    )
    parser.add_argument(
        "-b",
        "--annotate-regions",
        action="store_true",
        help=(
            "encode each residue's region type in the B-factor column so a viewer can colour "
            "by it (v1's -b/--beta_for_FD_IDR)"
        ),
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
            "how to identify regions. density is DODO's original all-atom density metric and "
            "the default for all-atom input; contact is a CA-only alternative; plddt is an "
            "explicit opt-in (default: auto)"
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
    for sub in (rebuild_parser, fetch_parser):
        sub.add_argument(
            "--domain-placement",
            default="predicted",
            choices=("predicted", "conformer"),
            help=(
                "how the folded domains get their positions. 'predicted' (default) moves "
                "them to ALBATROSS's predicted linker dimension and builds the linker into "
                "that gap. 'conformer' requires --engine starling, which is undocumented in "
                "2.0, so it is undocumented too"
            ),
        )

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

    validate_parser = subparsers.add_parser(
        "validate",
        help="check a structure's bond lengths, clashes and CONECT records",
        description=(
            "Validate structural geometry against reference bond lengths measured across 23,587 "
            "AlphaFold structures. Works on DODO output and on any other PDB or mmCIF file."
        ),
    )
    validate_parser.add_argument("structure", help="path to a PDB or mmCIF file")
    validate_parser.add_argument(
        "--max-findings",
        type=int,
        default=10,
        metavar="N",
        help="findings to print per check (default: 10)",
    )
    validate_parser.add_argument(
        "--no-clashes", action="store_true", help="skip the steric clash check"
    )
    validate_parser.add_argument(
        "--no-bonds", action="store_true", help="skip the bond-length check"
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

        if args.command == "validate":
            from .io import read_structure
            from .validate import validate_structure

            structure = read_structure(args.structure)
            path = (
                args.structure if str(args.structure).lower().endswith((".pdb", ".ent")) else None
            )
            result = validate_structure(
                structure,
                path=path,
                check_bonds=not args.no_bonds,
                check_clashes=not args.no_clashes,
            )
            print(result.describe(args.max_findings))
            if result.n_inherited_bond_violations:
                # Say this out loud. Otherwise a user validating DODO's own output reads
                # "INVALID" and concludes DODO broke their structure, when DODO faithfully
                # preserved a defect that arrived in the file -- folded-domain atoms are moved
                # rigidly and never regenerated.
                print(
                    f"\nnote: {result.n_inherited_bond_violations} of these bond findings are on "
                    f"atoms DODO did not build. If this file is DODO output, that geometry came "
                    f"in with the input and was preserved, not introduced.",
                    file=sys.stderr,
                )
            # 0 clean, 2 findings. Distinct from 1, which means the file could not be read.
            return 0 if result.ok else 2

        from .construct.pipeline import build_from_sequence, rebuild
        from .io import write_pdb

        common = {
            "mode": args.mode,
            "n_models": args.models,
            "engine": args.engine,
            "backbone": args.backbone,
            "seed": args.seed,
        }
        # None lets the pipeline decide from the terminal; -q means silence, bar included.
        progress: bool | None = False if args.quiet else None
        if getattr(args, "no_cache", False):
            # Set in the environment rather than threaded as a parameter: the prediction cache is
            # consulted deep inside dimension prediction, on paths that do not otherwise know a CLI
            # ran, and one switch covering all of them cannot drift out of sync with the flag.
            os.environ["DODO_NO_CACHE"] = "1"

        if args.command == "sequence":
            report = build_from_sequence(args.sequence, **common)
        elif args.command == "rebuild":
            report = rebuild(
                args.structure,
                strategy=args.strategy,
                progress=progress,
                domain_placement=args.domain_placement,
                **common,
            )
        elif args.command == "fetch":
            from .io import fetch_alphafold, resolve_uniprot_accession

            # Always resolve. resolve_uniprot_accession returns accession-shaped input
            # unchanged without touching the network, so there is nothing to save by guessing
            # first -- and the guess this used to make ("no spaces and 10 characters or fewer
            # must be an accession") was wrong in both directions: `dodo fetch p53` and
            # `dodo fetch calmodulin` failed, while `dodo fetch transthyretin` worked.
            accession = resolve_uniprot_accession(args.target)
            if accession != args.target:
                print(f"resolved {args.target!r} to {accession}", file=sys.stderr)
            if args.no_cache:
                # A TemporaryDirectory rather than the cache root, so the file is gone when the
                # command exits. The rebuild has to happen INSIDE the context manager; doing it
                # after would be reading a path that no longer exists.
                import tempfile
                from pathlib import Path

                with tempfile.TemporaryDirectory(prefix="dodo-fetch-") as scratch:
                    path = fetch_alphafold(accession, cache_dir=Path(scratch))
                    print(f"using {path} (not cached)", file=sys.stderr)
                    report = rebuild(
                        path,
                        strategy=args.strategy,
                        progress=progress,
                        domain_placement=args.domain_placement,
                        **common,
                    )
            else:
                path = fetch_alphafold(accession)
                print(f"using {path}", file=sys.stderr)
                report = rebuild(
                    path,
                    strategy=args.strategy,
                    progress=progress,
                    domain_placement=args.domain_placement,
                    **common,
                )
        else:  # pragma: no cover - argparse restricts the choices
            parser.error(f"unknown command {args.command!r}")

        status = _report(report, quiet=args.quiet)
        if report.models:
            write_pdb(
                report.models,
                args.out,
                conect=not args.no_conect,
                ca_only=args.ca_only,
                annotate_regions=args.annotate_regions,
            )
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
