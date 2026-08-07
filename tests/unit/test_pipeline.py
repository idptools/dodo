"""Tests for the end-to-end rebuild pipeline and the CLI.

The pipeline is the layer v1 had and the first v2 attempt did not: everything below it was
separately usable, but nothing wired it together. These tests cover the wiring, and in
particular the two things that are easy to get wrong when assembling the pieces by hand --
supplying the outer anchors so junction angles stay constrainable, and drawing a fresh
dimension target per model so a multi-model run is an ensemble rather than one conformation
repeated.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pytest

from dodo.cli import main
from dodo.constants import CA_CA_BOND_LENGTH
from dodo.construct.pipeline import build_from_sequence, rebuild
from dodo.geometry.metrics import end_to_end, validate_ca_trace
from dodo.io import read_structure
from dodo.structure import DomainKind
from dodo.validate import find_impossible_pairs, validate_bonds, validate_clashes

FIXTURES = Path(__file__).resolve().parents[1] / "data" / "structures"
DNMT3A = FIXTURES / "dnmt3a.pdb"

#: A generic disordered composition for the sequence-only path.
IDR_SEQUENCE = "SGQNTEKDRSGQNTPKAE" * 3


@pytest.mark.slow
class TestRebuild:
    def test_rebuilds_every_idr(self) -> None:
        report = rebuild(DNMT3A, seed=0)
        assert report.ok, report.summary()
        assert report.n_built > 0
        assert len(report.models) == 1

    def test_folded_domains_move_rigidly_and_are_never_rebuilt(self) -> None:
        """Folded-domain atoms are transformed, never regenerated.

        This test previously asserted folded-domain coordinates were bit-identical, which is
        the wrong invariant and hid the fact that step 3 of the algorithm was missing entirely.
        Folded domains DO move -- that is the whole point of repositioning them so a linker can
        reach its predicted dimensions -- they just move as rigid bodies. Their internal
        geometry must survive exactly; their position and orientation must not be assumed to.

        Loop residues are excluded, because a loop inside a folded domain IS rebuilt.
        """
        original = read_structure(DNMT3A)
        report = rebuild(DNMT3A, seed=0)
        rebuilt = report.models[0]

        for domain in rebuilt.domains:
            if domain.kind is not DomainKind.FOLDED:
                continue
            residues = np.arange(domain.span.start, domain.span.stop)
            in_loop = np.zeros(residues.size, dtype=bool)
            for loop in domain.loops:
                in_loop |= (residues >= loop.start) & (residues < loop.stop)

            before = original.ca_xyz[domain.span.slice][~in_loop]
            after = rebuilt.ca_xyz[domain.span.slice][~in_loop]
            if before.shape[0] < 2:
                continue
            # Distance to the domain's own centroid is invariant under any rigid motion, and
            # detects scaling, shearing or reflection.
            radii_before = np.linalg.norm(before - before.mean(axis=0), axis=1)
            radii_after = np.linalg.norm(after - after.mean(axis=0), axis=1)
            drift = float(np.abs(radii_before - radii_after).max())
            assert drift < 1e-6, f"{domain!r} was deformed, not moved: drift {drift:.2e} A"

    def test_repositioning_actually_happens(self) -> None:
        """Step 3 must run. Its absence is what made an earlier version's output wrong.

        AlphaFold packs domains joined by a long linker far closer than the linker predicts --
        measured, 2-3.6x closer on real models -- so if nothing moves, the linker is built into
        a gap that bears no relation to its sequence.
        """
        report = rebuild(DNMT3A, seed=0)
        assert report.placements, "no folded domains were considered for repositioning"
        moved = [p for p in report.placements if p.moved]
        assert moved, "no folded domain moved; step 3 did not run"
        for placement in moved:
            assert placement.target_separation is not None
            assert placement.achieved_separation == pytest.approx(
                placement.target_separation, abs=0.5
            ), str(placement)

    def test_linkers_reach_their_predicted_dimensions(self) -> None:
        """The point of the whole exercise: a connecting IDR ends up the right size."""
        report = rebuild(DNMT3A, seed=0)
        connecting = [
            o
            for o in report.outcomes
            if o.built and o.target is not None and o.requested_end_to_end is not None
        ]
        assert connecting
        for outcome in connecting:
            error = abs(outcome.achieved_end_to_end - outcome.requested_end_to_end)
            relative = error / outcome.requested_end_to_end
            assert relative < 0.15, str(outcome)

    def test_loops_are_rebuilt(self) -> None:
        """Loops are a distinct region type and must actually be built.

        An earlier version identified them and then only ever rebuilt IDRs, so loops were
        reported in the region assignment and then silently left alone.
        """
        report = rebuild(DNMT3A, seed=0)
        has_loops = any(d.loops for d in report.models[0].domains)
        assert has_loops, "fixture needs a folded domain with a loop for this to mean anything"
        assert any("loop" in (o.reason or "") or o.target is None for o in report.outcomes), (
            "no loop appears in the outcomes; loops are not being rebuilt"
        )

    def test_rebuilt_regions_contain_only_alpha_carbons(self) -> None:
        """A rebuilt region is CA-only; a folded domain keeps every atom.

        Regression test for a bug that was visible in a viewer. set_ca_xyz() writes only alpha
        carbons, so the N/C/O and side-chain atoms of a rebuilt residue stayed at their original
        AlphaFold positions while the CA moved. Measured on p300, that left each rebuilt residue
        split across ~93 A: the writer then emitted a CONECT record bonding N to CA over that
        distance, which renders as a long spurious straight line, and the orphaned atoms trailed
        along the region's old path as disconnected dots.
        """
        report = rebuild(DNMT3A, seed=0)
        structure = report.models[0]
        for domain in structure.domains:
            atoms = structure.atom_slice_for_residues(domain.span.start, domain.span.stop)
            names = set(structure.atom_name[atoms].tolist())
            if domain.kind is DomainKind.IDR and domain.rebuilt:
                assert names == {"CA"}, f"{domain!r} kept non-CA atoms: {sorted(names)}"
            elif domain.kind is DomainKind.FOLDED and not domain.loops:
                assert len(names) > 4, f"{domain!r} lost its side chains: {sorted(names)}"

    def test_no_residue_is_split_across_two_locations(self) -> None:
        """Every atom stays near its own alpha carbon.

        The direct assertion of the defect above: a residue's furthest heavy atom is ~6-7 A from
        its CA (arginine's terminal nitrogens). Anything in the tens or hundreds means the
        residue's atoms did not move together.
        """
        report = rebuild(DNMT3A, seed=0)
        structure = report.models[0]
        worst = 0.0
        for residue in range(structure.n_residues):
            atoms = structure.atom_slice_for_residues(residue, residue + 1)
            names = structure.atom_name[atoms]
            if "CA" not in names:
                continue
            coords = structure.xyz[atoms]
            ca = coords[list(names).index("CA")]
            worst = max(worst, float(np.linalg.norm(coords - ca, axis=1).max()))
        assert worst < 10.0, f"a residue's atoms are {worst:.1f} A apart; it was split"

    def test_conect_records_are_all_physical_bonds(self, tmp_path: Path) -> None:
        """No CONECT record may span more than a real bond.

        CA-CA is 3.81 A and every other bond DODO writes is shorter, so anything much above
        that means CONECT is bonding the wrong atoms -- which is what a split residue produced.
        """
        from dodo.io import write_pdb

        report = rebuild(DNMT3A, seed=0)
        out = tmp_path / "conect.pdb"
        write_pdb(report.models, out, conect=True)

        positions: dict[int, np.ndarray] = {}
        bonds: list[tuple[int, int]] = []
        for line in out.read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                positions[int(line[6:11])] = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
            elif line.startswith("CONECT"):
                origin = int(line[6:11])
                for column in range(11, min(len(line), 31), 5):
                    field = line[column : column + 5].strip()
                    if field:
                        bonds.append((origin, int(field)))

        assert bonds, "no CONECT records were written"
        lengths = np.array(
            [
                np.linalg.norm(positions[a] - positions[b])
                for a, b in bonds
                if a in positions and b in positions
            ]
        )
        assert lengths.max() < 4.5, (
            f"longest CONECT bond is {lengths.max():.1f} A; CONECT is bonding distant atoms"
        )

    def test_regions_are_not_reassigned_after_repositioning(self) -> None:
        """The reported assignment must match the model's actual domains.

        Region identification reads the coordinates, and repositioning moves folded domains
        apart -- which changes their contact density. Re-assigning afterwards shifted a p300
        domain's bounds from 569-650 to 570-644, so the anchors that drove the placement were
        no longer the anchors built against.
        """
        report = rebuild(DNMT3A, seed=0)
        reported = [(d.kind, d.span.start, d.span.stop) for d in report.assignments[0].domains]
        actual = [(d.kind, d.span.start, d.span.stop) for d in report.models[0].domains]
        assert reported == actual

    def test_rebuilt_regions_are_physically_valid(self) -> None:
        """Every rebuilt region passes the clash-aware gate, junctions included.

        Validated over the region *plus its flanking anchors*, because the junction angles
        belong to the assembled chain rather than to the region in isolation -- which is
        exactly the distinction the pre-rewrite code missed.
        """
        report = rebuild(DNMT3A, seed=0)
        structure = report.models[0]
        for domain in structure.domains:
            if domain.kind is not DomainKind.IDR or not domain.rebuilt:
                continue
            start = domain.span.n_anchor if domain.span.n_anchor is not None else domain.span.start
            stop = (
                domain.span.c_anchor + 1 if domain.span.c_anchor is not None else domain.span.stop
            )
            trace = structure.ca_xyz[start:stop]
            report_ = validate_ca_trace(trace, residue_offset=start)
            assert report_.ok, f"{domain!r}: {report_.describe()}"

    def test_a_single_model_hits_the_target(self) -> None:
        """One model should reach the predicted dimension, not a random draw around it."""
        report = rebuild(DNMT3A, seed=0)
        built = [o for o in report.outcomes if o.built and o.target is not None]
        assert built
        for outcome in built:
            assert outcome.requested_end_to_end == pytest.approx(
                outcome.target.end_to_end, rel=1e-9
            )

    def test_multiple_models_are_an_ensemble(self) -> None:
        """The headline scientific fix: models must scatter around the predicted mean.

        v1 placed folded domains once outside the model loop and targeted only the mean, so
        every model shared one arrangement and essentially one end-to-end distance. Measured on
        the first v2 engine, CV of Re across conformers was 0.006-0.045 where a matched physical
        reference gives 0.35-0.48 -- sixty models of a 200-residue IDR spanning 1.9 A of
        extension is one conformation sampled sixty times.
        """
        report = rebuild(DNMT3A, n_models=6, seed=1)
        assert len(report.models) == 6

        terminal = [
            d
            for d in report.models[0].domains
            if d.kind is DomainKind.IDR and d.span.is_terminal and len(d) > 20
        ]
        assert terminal, "fixture needs a long free-ended IDR for this test to mean anything"
        span = terminal[0].span

        distances = np.array([end_to_end(model.ca_xyz[span.slice]) for model in report.models])
        cv = float(distances.std() / distances.mean())
        assert cv > 0.15, f"models barely differ (CV {cv:.3f}); this is not an ensemble"

    def test_models_are_independent_of_each_other(self) -> None:
        """A failure or an odd conformation in one model must not affect the next."""
        report = rebuild(DNMT3A, n_models=3, seed=2)
        coords = [m.xyz for m in report.models]
        assert not np.array_equal(coords[0], coords[1])
        assert not np.array_equal(coords[1], coords[2])

    def test_is_reproducible(self) -> None:
        first = rebuild(DNMT3A, n_models=2, seed=7)
        second = rebuild(DNMT3A, n_models=2, seed=7)
        for a, b in zip(first.models, second.models, strict=True):
            assert np.array_equal(a.xyz, b.xyz)

    def test_different_seeds_differ(self) -> None:
        first = rebuild(DNMT3A, seed=1)
        second = rebuild(DNMT3A, seed=2)
        assert not np.array_equal(first.models[0].xyz, second.models[0].xyz)

    @pytest.mark.parametrize("mode", ["compact", "predicted", "expanded"])
    def test_modes_change_the_target(self, mode: str) -> None:
        report = rebuild(DNMT3A, mode=mode, seed=0)
        built = [o for o in report.outcomes if o.built and o.target is not None]
        assert built
        assert all(o.target.mode == mode for o in built)

    def test_accepts_an_already_parsed_structure(self) -> None:
        structure = read_structure(DNMT3A)
        report = rebuild(structure, seed=0)
        assert report.ok

    def test_report_carries_the_region_assignment(self) -> None:
        report = rebuild(DNMT3A, seed=0)
        assert report.assignments
        assert "chain A" in report.assignments[0].describe()

    def test_zero_models_rejected(self) -> None:
        with pytest.raises(ValueError, match="at least 1"):
            rebuild(DNMT3A, n_models=0)

    def test_unknown_engine_rejected(self) -> None:
        with pytest.raises(ValueError, match="Unknown engine"):
            rebuild(DNMT3A, engine="magic")


class TestAnchorObstacles:
    """The anchor of a rebuilt region is only *partly* exempt from clash checking.

    Exempting the whole anchor residue looks conservative and is not. The first rebuilt
    residue is bonded to the anchor's CA but has no bonded relationship to the anchor's side
    chain, so a blanket exemption let the walk place a CA straight through it -- producing
    overlaps at 0.871 A (LEU CD1), 0.937 A (ASN ND2) and 0.944 A (LYS CD) in output that
    every other check called clean. Those are below the shortest bond in any protein.
    """

    def test_anchor_backbone_is_exempt_but_side_chain_is_not(self) -> None:
        from dodo.construct.pipeline import _obstacles_for_span
        from dodo.regions import assign_regions

        structure = read_structure(DNMT3A)
        assign_regions(structure)
        # Mark everything placed, so the obstacle set is limited only by the exemptions.
        # `placed`, not `rebuilt`: the obstacle set is about final coordinates, not provenance.
        for domain in structure.domains:
            domain.placed = True
        span = next(
            d.span
            for d in structure.domains
            if d.kind is DomainKind.IDR and d.span.c_anchor is not None
        )
        anchor = span.c_anchor
        assert anchor is not None

        obstacles = _obstacles_for_span(structure, span)
        assert obstacles is not None
        present = {tuple(np.round(row, 4)) for row in obstacles}

        atoms = structure.atom_slice_for_residues(anchor, anchor + 1)
        names = [str(n) for n in structure.atom_name[atoms]]
        coords = structure.xyz[atoms]
        # A residue with a side chain, or the assertion below proves nothing.
        assert any(n not in ("N", "CA", "C", "O") for n in names)

        for name, xyz in zip(names, coords, strict=True):
            key = tuple(np.round(xyz, 4))
            if name in ("N", "CA", "C", "O"):
                assert key not in present, f"anchor backbone {name} must be exempt"
            else:
                assert key in present, f"anchor side-chain {name} must remain an obstacle"

    def test_the_alpha_carbon_exemption_is_unconditional(self) -> None:
        """Always exempt, and separately from the discretionary backbone exemption.

        A rebuilt region is bonded to its anchors' alpha carbons. Treating those as obstacles
        would make every valid attachment register as a clash, so there is no version of the
        algorithm without this exemption -- it is not a trade-off to be tuned.
        """
        from dodo.constants import ANCHOR_ALWAYS_EXEMPT_ATOMS, ANCHOR_EXEMPT_ATOMS

        assert frozenset({"CA"}) == ANCHOR_ALWAYS_EXEMPT_ATOMS
        # The discretionary set is the anchor BACKBONE. It happens to contain CA as well, which
        # is harmless -- but the unconditional set is what guarantees the alpha carbon stays
        # exempt even when the backbone exemption is withheld.
        assert "CA" in ANCHOR_EXEMPT_ATOMS
        assert "N" not in ANCHOR_ALWAYS_EXEMPT_ATOMS, "N must be discretionary, not unconditional"

    def test_strict_pass_keeps_anchor_backbone_as_an_obstacle(self) -> None:
        """The strict setting exempts only the alpha carbons; the fallback adds the backbone."""
        from dodo.construct.pipeline import _obstacles_for_span
        from dodo.regions import assign_regions

        structure = read_structure(DNMT3A)
        assign_regions(structure)
        for domain in structure.domains:
            domain.placed = True
        span = next(
            d.span
            for d in structure.domains
            if d.kind is DomainKind.IDR and d.span.c_anchor is not None
        )
        anchor = span.c_anchor
        assert anchor is not None
        atoms = structure.atom_slice_for_residues(anchor, anchor + 1)
        names = [str(n) for n in structure.atom_name[atoms]]
        coords = structure.xyz[atoms]

        def present(*, exempt_backbone: bool) -> set[tuple[float, ...]]:
            obstacles = _obstacles_for_span(structure, span, exempt_anchor_backbone=exempt_backbone)
            assert obstacles is not None
            return {tuple(np.round(row, 4)) for row in obstacles}

        strict, relaxed = present(exempt_backbone=False), present(exempt_backbone=True)
        for name, xyz in zip(names, coords, strict=True):
            key = tuple(np.round(xyz, 4))
            if name == "CA":
                assert key not in strict, "the anchor CA must be exempt in BOTH passes"
                assert key not in relaxed
            elif name in ("N", "C", "O"):
                assert key in strict, f"strict pass must keep anchor {name} as an obstacle"
                assert key not in relaxed, f"fallback must exempt anchor {name}"
            else:
                assert key in strict, f"anchor side-chain {name} is never exempt"
                assert key in relaxed, f"anchor side-chain {name} is never exempt"

    @pytest.mark.slow
    def test_a_relaxed_build_says_so(self) -> None:
        """Whenever the fallback is used it must be visible, not silent.

        The relaxed pass lets a region sit closer to its anchors' backbone than the clash
        distance. That is a deliberate trade -- an unbuilt region is far more visible in a figure
        than a marginal contact -- but the user has to be able to tell it happened.
        """
        report = rebuild(DNMT3A, seed=0)
        for outcome in report.outcomes:
            if outcome.built and outcome.reason:
                assert "relaxed anchor exemption" in outcome.reason, outcome.reason
        # And the wording reaches the summary a CLI user actually reads.
        relaxed = [o for o in report.outcomes if o.built and o.reason]
        if relaxed:
            assert "relaxed anchor exemption" in report.summary()

    def test_proline_cd_is_exempt(self) -> None:
        """Proline's CD is bonded to its own backbone N, so it is 1-3 from the preceding C.

        Measured minimum to a neighbouring CA is 2.245 A -- below the clash distance, and
        legitimately so. It is the only side-chain atom with that exemption.
        """
        from dodo.constants import ANCHOR_EXEMPT_ATOMS, ANCHOR_EXEMPT_ATOMS_BY_RESIDUE

        assert "CD" in ANCHOR_EXEMPT_ATOMS_BY_RESIDUE["PRO"]
        assert "CD" not in ANCHOR_EXEMPT_ATOMS
        for residue in ("GLU", "GLN", "LYS", "ARG"):
            assert residue not in ANCHOR_EXEMPT_ATOMS_BY_RESIDUE

    @pytest.mark.slow
    def test_output_has_no_impossible_separations(self) -> None:
        """The end-to-end guard: no rebuilt structure may contain a sub-bond-length contact."""
        from dodo.validate import find_impossible_pairs

        for seed in (0, 1, 2, 3):
            report = rebuild(DNMT3A, seed=seed)
            pairs = find_impossible_pairs(report.models[0])
            assert not pairs, f"seed {seed}: {[p.message for p in pairs]}"


class TestBuildFromSequence:
    def test_builds_a_free_chain(self) -> None:
        report = build_from_sequence(IDR_SEQUENCE, seed=0)
        assert report.n_built == 1
        assert len(report.models) == 1
        assert report.models[0].n_residues == len(IDR_SEQUENCE)

    def test_sequence_round_trips(self) -> None:
        report = build_from_sequence(IDR_SEQUENCE, seed=0)
        assert report.models[0].sequence == IDR_SEQUENCE

    def test_geometry_is_valid(self) -> None:
        report = build_from_sequence(IDR_SEQUENCE, seed=0)
        assert validate_ca_trace(report.models[0].ca_xyz).ok

    def test_multiple_conformers_scatter(self) -> None:
        report = build_from_sequence(IDR_SEQUENCE, n_models=8, seed=3)
        distances = np.array([end_to_end(m.ca_xyz) for m in report.models])
        assert float(distances.std() / distances.mean()) > 0.15

    def test_is_reproducible(self) -> None:
        first = build_from_sequence(IDR_SEQUENCE, seed=5)
        second = build_from_sequence(IDR_SEQUENCE, seed=5)
        assert np.array_equal(first.models[0].xyz, second.models[0].xyz)

    @pytest.mark.parametrize("bad", ["", "   ", "ACDE FGH", "ACDE-1"])
    def test_malformed_sequence_rejected(self, bad: str) -> None:
        with pytest.raises(ValueError):
            build_from_sequence(bad)

    def test_starling_is_refused_for_now(self) -> None:
        with pytest.raises(ValueError, match="engine='walk'"):
            build_from_sequence(IDR_SEQUENCE, engine="starling")


class TestCli:
    def test_help_exits_cleanly(self, capsys: pytest.CaptureFixture[str]) -> None:
        assert main([]) == 0
        assert "rebuild" in capsys.readouterr().out

    def test_version(self, capsys: pytest.CaptureFixture[str]) -> None:
        assert main(["--version"]) == 0
        assert capsys.readouterr().out.strip()

    def test_regions_subcommand(self, capsys: pytest.CaptureFixture[str]) -> None:
        assert main(["regions", str(DNMT3A)]) == 0
        assert "chain A" in capsys.readouterr().out

    @pytest.mark.slow
    def test_rebuild_subcommand_writes_a_file(self, tmp_path: Path) -> None:
        out = tmp_path / "out.pdb"
        assert main(["rebuild", str(DNMT3A), "-o", str(out), "--seed", "0", "-q"]) == 0
        assert out.exists()
        # Must be readable back by DODO's own reader.
        assert read_structure(out).n_residues > 0

    def test_sequence_subcommand_writes_a_file(self, tmp_path: Path) -> None:
        out = tmp_path / "seq.pdb"
        assert main(["sequence", IDR_SEQUENCE, "-o", str(out), "--seed", "0", "-q"]) == 0
        assert out.exists()

    @pytest.mark.slow
    def test_multi_model_output_has_model_records(self, tmp_path: Path) -> None:
        """Without MODEL/ENDMDL the pseudo-trajectory feature has no output format."""
        out = tmp_path / "models.pdb"
        assert main(["rebuild", str(DNMT3A), "-o", str(out), "-n", "3", "--seed", "0", "-q"]) == 0
        lines = out.read_text().splitlines()
        assert sum(1 for line in lines if line.startswith("MODEL")) == 3
        assert sum(1 for line in lines if line.startswith("ENDMDL")) == 3

    @pytest.mark.slow
    def test_conect_records_are_written_by_default(self, tmp_path: Path) -> None:
        """CA-CA spacing exceeds viewer auto-bond cutoffs, so this is not optional polish."""
        out = tmp_path / "c.pdb"
        main(["rebuild", str(DNMT3A), "-o", str(out), "--seed", "0", "-q"])
        assert any(line.startswith("CONECT") for line in out.read_text().splitlines())

    def test_a_missing_file_is_a_clean_error_not_a_traceback(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        status = main(["rebuild", str(tmp_path / "nope.pdb"), "-o", str(tmp_path / "x.pdb")])
        assert status == 1
        assert "dodo:" in capsys.readouterr().err

    def test_unknown_mode_is_rejected_by_argparse(self) -> None:
        with pytest.raises(SystemExit):
            main(["rebuild", str(DNMT3A), "-o", "x.pdb", "-m", "very_squished"])


class TestPresetRegions:
    """The granular-control path, and the replacement for v1's ``regions_dict=``.

    v1 took a parallel, stringly-typed description of the structure and tried to reconcile it
    with the real one. The author's assessment of that design was blunt -- "a very bad idea" --
    and the failure mode bears it out: the two representations could disagree, and v1 accepted
    overlaps, gaps and out-of-range bounds silently before failing much later with something
    unrelated.

    So there is no ``regions`` parameter. Instead the caller assigns regions onto the structure
    however they like and asks :func:`~dodo.rebuild` to build exactly those. One representation,
    already validated, carrying the score profile and threshold that produced it.
    """

    def test_rebuild_honours_caller_supplied_regions_verbatim(self) -> None:
        from dodo.regions import assign_regions_from_spec

        structure = read_structure(DNMT3A)
        spec = [("idr", 1, 60), ("folded", 61, 800), ("idr", 801, 912)]
        assign_regions_from_spec(structure, {"A": spec})

        report = rebuild(structure, strategy="preset", seed=0)

        got = [
            (d.kind.value, d.span.start + 1, d.span.stop)
            for d in sorted(report.models[0].domains, key=lambda d: d.span.start)
        ]
        assert got == spec, f"preset regions were not honoured: {got}"
        assert report.ok, report.summary()
        # Both IDRs were rebuilt, and the folded domain was not.
        assert report.n_built == 2

    def test_preset_differs_from_what_dodo_would_have_chosen(self) -> None:
        """Guards the premise: if the spec matched the automatic call, the test above is vacuous."""
        from dodo.regions import assign_regions, assign_regions_from_spec

        automatic = assign_regions(read_structure(DNMT3A))[0]
        auto_bounds = [(d.kind.value, d.span.start + 1, d.span.stop) for d in automatic.domains]

        structure = read_structure(DNMT3A)
        spec = [("idr", 1, 60), ("folded", 61, 800), ("idr", 801, 912)]
        assign_regions_from_spec(structure, {"A": spec})
        assert auto_bounds != spec, "the spec must differ from the automatic assignment"

    def test_preset_without_any_assignment_explains_itself(self) -> None:
        from dodo.exceptions import InvalidRegionError

        with pytest.raises(InvalidRegionError, match="assign_regions_from_spec"):
            rebuild(DNMT3A, strategy="preset")

    def test_preset_reports_that_it_identified_nothing(self) -> None:
        """A NaN score and threshold, because none was computed. Zero would read as measured."""
        from dodo.regions import assign_regions, assign_regions_from_spec

        structure = read_structure(DNMT3A)
        assign_regions_from_spec(structure, {"A": [("idr", 1, 60), ("folded", 61, 912)]})
        assignment = assign_regions(structure, strategy="preset")[0]

        assert assignment.strategy.value == "preset"
        assert np.isnan(assignment.threshold)
        assert np.all(np.isnan(assignment.score))
        assert any("supplied by the caller" in note for note in assignment.notes)
        # folded_mask is still real, since it is derivable from the domains themselves.
        assert assignment.folded_mask[0] is np.False_ or not assignment.folded_mask[0]
        assert assignment.folded_mask[-1]


class TestShortRegionsAreTolerated:
    """A short region DODO cannot rebuild is reported, not treated as a failed run.

    DODO is a visualization tool first. A handful of residues left as AlphaFold drew them does not
    look wrong in a figure -- it is the long regions, the ones that trail across the image as
    extended spaghetti, that DODO exists to fix. So the threshold is about what a reader would
    actually notice, not about what the builder would prefer.

    Measured on the 117-structure corpus, this changes one outcome: a 7-residue terminal tail on
    AF-O14683-F1 that the walk cannot fit. The 16-residue loop and 71-residue linker that also fail
    stay failures, and both of those are input defects -- one file has two fixed residues 3.04 A
    apart, the other a chain break with consecutive alpha carbons 5.26 A apart.
    """

    def _outcome(self, *, n_residues: int, built: bool) -> object:
        from dodo.construct.pipeline import RegionOutcome

        return RegionOutcome(
            model=1,
            chain_id="A",
            residues=(1, n_residues),
            n_residues=n_residues,
            built=built,
            reason=None if built else "could not be built",
        )

    def test_the_threshold_is_ten_residues(self) -> None:
        from dodo.constants import SHORT_REGION_TOLERANCE

        assert SHORT_REGION_TOLERANCE == 10

    @pytest.mark.parametrize(
        ("n_residues", "tolerated"), [(1, True), (9, True), (10, False), (71, False)]
    )
    def test_tolerance_is_decided_by_length(self, n_residues: int, *, tolerated: bool) -> None:
        outcome = self._outcome(n_residues=n_residues, built=False)
        assert outcome.tolerated is tolerated

    def test_a_built_region_is_never_tolerated(self) -> None:
        """`tolerated` describes a failure, so a success must not report it."""
        assert self._outcome(n_residues=3, built=True).tolerated is False

    def test_ok_ignores_short_failures_but_not_long_ones(self) -> None:
        from dodo.construct.pipeline import RebuildReport

        short = RebuildReport(outcomes=[self._outcome(n_residues=7, built=False)])
        assert short.ok, "a 7-residue region left as-is must not fail the run"
        assert short.failures and not short.blocking_failures
        assert short.tolerated_failures

        long = RebuildReport(outcomes=[self._outcome(n_residues=71, built=False)])
        assert not long.ok, "a 71-residue region left unbuilt is a real failure"
        assert long.blocking_failures and not long.tolerated_failures

    def test_both_kinds_are_named_distinctly_in_the_summary(self) -> None:
        """A tolerated region must not be printed as though it were a failure."""
        from dodo.construct.pipeline import RebuildReport

        report = RebuildReport(
            outcomes=[
                self._outcome(n_residues=7, built=False),
                self._outcome(n_residues=71, built=False),
            ]
        )
        summary = report.summary()
        assert "1 failure(s)" in summary
        assert "left as-is" in summary
        assert "NOT BUILT" in summary


class TestMetapredictIsGone:
    """metapredict was dropped along with its only reason for existing.

    In 1.x it provided faster region identification than the all-atom density metric. That metric
    now runs in 7 ms on a 1,086-residue model, down from 10.1 s, so the tradeoff is gone -- and
    metapredict requires torch, pytorch-lightning, cython and matplotlib, which is most of the
    weight a light install avoids.
    """

    def test_the_strategy_no_longer_exists(self) -> None:
        from dodo.regions import Strategy

        assert not hasattr(Strategy, "METAPREDICT")
        with pytest.raises(ValueError, match="metapredict"):
            Strategy("metapredict")

    def test_the_cli_does_not_offer_it(self) -> None:
        from dodo.cli import _STRATEGY_CHOICES

        assert "metapredict" not in _STRATEGY_CHOICES
        assert set(_STRATEGY_CHOICES) == {"auto", "density", "contact", "plddt"}

    def test_nothing_imports_it(self) -> None:
        import pathlib

        root = pathlib.Path(__file__).resolve().parents[2] / "src" / "dodo"
        offenders = [
            path.relative_to(root)
            for path in root.rglob("*.py")
            if "metapredict" in path.read_text()
        ]
        assert not offenders, f"metapredict still referenced in {offenders}"


class TestBackboneFlag:
    """The opt-in ``backbone=`` flag on the two entry points.

    Opt-in rather than default because of the seams: where a rebuilt region joins a residue DODO did
    not touch, that residue's existing N or C still points toward where the region ran in the
    *original* model, and folded-domain atoms are not DODO's to move. See
    :meth:`test_folded_domains_keep_every_atom`.
    """

    def test_off_by_default(self) -> None:
        report = build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", seed=0)
        assert {str(n) for n in report.models[0].atom_name} == {"CA"}

    def test_on_when_asked(self) -> None:
        report = build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", seed=0, backbone=True)
        assert {str(n) for n in report.models[0].atom_name} == {"N", "CA", "C", "O"}

    def test_alpha_carbons_are_identical_either_way(self) -> None:
        """The flag adds atoms; it must not change DODO's actual answer.

        Same seed, so the alpha carbons are the same coordinates to the bit. If this drifts, the
        backbone pass is perturbing the trace rather than decorating it.
        """
        plain = build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", seed=0).models[0]
        with_backbone = build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", seed=0, backbone=True).models[
            0
        ]
        assert np.array_equal(plain.ca_xyz, with_backbone.ca_xyz)

    def test_folded_domains_keep_every_atom(self) -> None:
        """The constraint that makes this additive: folded domains are untouched.

        DODO never regenerates folded-domain geometry, so a folded domain must come through with
        its side chains intact whether or not the rebuilt regions gained a backbone.
        """
        source = FIXTURES / "dnmt3a.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(source, seed=0, backbone=True)
        model = report.models[0]
        names = {str(n) for n in model.atom_name}
        # Side-chain atoms only a real folded domain has; a CA-plus-backbone output has four names.
        assert len(names) > 10, f"folded domains lost their side chains; only {names}"
        assert {"CB", "CG", "N", "C", "O", "CA"} <= names

    def test_every_rebuilt_residue_gets_exactly_a_backbone(self) -> None:
        """Generated residues end up with N, CA, C, O and nothing else."""
        source = FIXTURES / "dnmt3a.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fancy = rebuild(source, seed=0, backbone=True).models[0]
        generated = {
            residue
            for domain in fancy.domains
            for span in domain.generated_spans()
            for residue in range(span.start, span.stop)
        }
        assert generated, "nothing was rebuilt, so this proves nothing"
        for residue in sorted(generated):
            atoms = fancy.atom_slice_for_residues(residue, residue + 1)
            names = {str(n) for n in fancy.atom_name[atoms]}
            assert names == {"N", "CA", "C", "O"}, f"residue {residue} has {sorted(names)}"

    def test_every_generated_bond_violation_is_at_a_seam(self) -> None:
        """The seams are strained, and nothing else is. That is the measured limit of this flag.

        An exact peptide bond onto a folded domain is not reachable from a rebuilt alpha carbon: a
        peptide unit spans at most 2.854 A to the nitrogen it bonds to, and the measured distance
        across 17 seams in three structures was 3.2-4.5 A. So the requirement is not that these
        bonds be ideal -- they cannot be -- but that every violation sit AT a seam, so the interior
        of a rebuilt region is known to be clean.

        Leaving the seam residue un-rebuilt does not fix it, which is worth recording because it is
        the obvious idea and it very nearly works. Its alpha carbon being input geometry does make
        the bond reachable -- 2.45-2.52 A at all 17 seams. But closing onto it means re-placing its
        nitrogen, and that residue's side chain was built around where its nitrogen used to be, so
        the new one is driven into its own CB (measured 1.405 A against a correct 2.45) and, for
        proline, snaps the ring bond to CD (3.444 A against 1.47).
        """
        source = FIXTURES / "dnmt3a.pdb"
        for seed in (0, 1, 2):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fancy = rebuild(source, seed=seed, backbone=True).models[0]
            seams = {
                residue
                for domain in fancy.domains
                for span in domain.generated_spans()
                for residue in (span.start, span.start - 1, span.stop - 1, span.stop)
            }
            for violation in validate_bonds(fancy).violations:
                if violation.provenance != "rebuilt":
                    continue
                assert set(violation.residue_indices) & seams, (
                    f"seed {seed}: a violation away from any seam -- {violation.message}"
                )

    def test_introduces_no_impossible_contacts(self) -> None:
        """Whatever the seams do, they must not put two atoms on top of each other.

        The seam fallback exists for this: when no carbon can satisfy both the CA-C bond and the
        peptide bond to an untouched neighbour, it aims at the neighbour and leaves the seam bond
        long rather than writing an impossible pair. A strained bond is a visible, reportable
        compromise; a 0.6 A contact is a broken file.
        """
        source = FIXTURES / "dnmt3a.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fancy = rebuild(source, seed=0, backbone=True).models[0]
        inherited = {
            (p.residue_labels, p.atom_names) for p in find_impossible_pairs(read_structure(source))
        }
        introduced = [
            p
            for p in find_impossible_pairs(fancy)
            if (p.residue_labels, p.atom_names) not in inherited
        ]
        assert introduced == [], f"backbone placement introduced {introduced}"


class TestBackboneDoesNotDamageInputGeometry:
    """The guard that catches the seam "fix" that keeps looking correct and is not.

    A rebuilt alpha carbon sits 3.2-4.5 A from the folded domain's fixed nitrogen, where a peptide
    unit reaches only 2.854 A, so the seam bond is geometrically unsatisfiable and DODO leaves it
    long. The obvious repair is to stop rebuilding one residue short, so the seam falls on input
    geometry -- and measured, that does make the distance 2.45-2.52 A at all 17 seams.

    It was implemented on that basis and it is wrong, because closing an exact bond onto that
    residue means re-placing its nitrogen, and the residue's side chain was built around where that
    nitrogen used to be. Measured on dnmt3a: PRO282's ring bond CD-N stretched to 2.263 A against an
    accepted 1.378-1.705, and GLU473's new N landed 1.683 A from its own CB. Trading four strained
    backbone bonds for a snapped proline ring is not a trade worth making.

    What makes this worth a test rather than a comment is how it hid: both defects are reported with
    ``input`` provenance, because a residue that was deliberately not rebuilt is not DODO's work by
    the provenance rules. Any check filtered to ``provenance == "rebuilt"`` -- the natural thing to
    write when asking "did DODO break anything" -- sees a clean run.
    """

    def test_no_new_violation_inside_a_single_residue(self) -> None:
        """The check that catches the seam experiment without flagging the accepted compromise.

        Deliberately scoped to violations WITHIN one residue. The strained seam bonds are between
        two residues and are the known, reported trade; flagging those would only mean deleting this
        test the first time somebody read it. What must never appear is a residue whose own internal
        geometry got damaged -- exactly how the seam experiment failed, on PRO282's ring bond CD-N
        and GLU473's N-to-CB separation.

        Provenance is deliberately not filtered. Both of those defects were reported as ``input``,
        because a residue DODO did not rebuild is not DODO's work by the provenance rules, so the
        natural check -- filtered to ``rebuilt`` -- saw a clean run.
        """
        source = FIXTURES / "dnmt3a.pdb"
        baseline = {v.message for v in validate_bonds(read_structure(source)).violations}
        for seed in (0, 1, 2):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model = rebuild(source, seed=seed, backbone=True).models[0]
            introduced = [
                v
                for v in validate_bonds(model).violations
                if v.message not in baseline and len(set(v.residue_indices)) == 1
            ]
            assert introduced == [], (
                f"seed {seed} damaged a residue's own geometry: "
                f"{[v.message[:90] for v in introduced]}"
            )

    def test_seam_bonds_are_the_only_thing_introduced(self) -> None:
        """And they are between-residue, bounded, and reported rather than silent."""
        source = FIXTURES / "dnmt3a.pdb"
        baseline = {v.message for v in validate_bonds(read_structure(source)).violations}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = rebuild(source, seed=0, backbone=True).models[0]
        introduced = [v for v in validate_bonds(model).violations if v.message not in baseline]
        assert introduced, "expected the known seam compromise to show up"
        assert all(len(set(v.residue_indices)) == 2 for v in introduced)
        assert len(introduced) < 10, f"{len(introduced)} is more seams than dnmt3a has"

    def test_side_chain_geometry_of_untouched_residues_survives(self) -> None:
        """No atom of a residue DODO did not rebuild may move at all.

        Stated as exact equality rather than a tolerance, because there is no legitimate reason for
        one of these coordinates to change: a residue outside every generated span is either part of
        a rigidly-moved folded domain or was left alone entirely, and in both cases the backbone
        pass has no business writing to it.
        """
        source = FIXTURES / "dnmt3a.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plain = rebuild(source, seed=0).models[0]
            fancy = rebuild(source, seed=0, backbone=True).models[0]
        generated = {
            residue
            for domain in fancy.domains
            for span in domain.generated_spans()
            for residue in range(span.start, span.stop)
        }
        untouched = [r for r in range(fancy.n_residues) if r not in generated]
        assert untouched, "nothing was left untouched, so this proves nothing"
        for residue in untouched:
            a = plain.atom_slice_for_residues(residue, residue + 1)
            b = fancy.atom_slice_for_residues(residue, residue + 1)
            assert {str(n) for n in plain.atom_name[a]} == {str(n) for n in fancy.atom_name[b]}, (
                f"residue {residue} gained or lost an atom"
            )
            assert np.array_equal(plain.xyz[a], fancy.xyz[b]), f"residue {residue} moved"

    def test_alpha_carbons_are_identical_with_and_without_backbone(self) -> None:
        """``--backbone`` is purely additive: it decorates the trace, it does not change it.

        This held before the seam experiment, stopped holding during it (shortening every region
        moved the anchor the walk closed onto) and holds again. Worth pinning: a flag whose name
        promises extra atoms should not quietly return a different structure.
        """
        source = FIXTURES / "dnmt3a.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            plain = rebuild(source, seed=0).models[0]
            fancy = rebuild(source, seed=0, backbone=True).models[0]
        assert np.array_equal(plain.ca_xyz, fancy.ca_xyz)


class TestBatchEngine:
    """The vectorized batch engine, wired in as ``engine='batch'``, rebuilds valid structures."""

    def test_rebuild_with_the_batch_engine_builds_every_region_with_exact_bonds(self) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(DNMT3A, engine="batch", seed=0, progress=False)
        assert report.ok
        assert not report.blocking_failures
        # The batch engine (with its walk fallback) builds every region the walk does, and the
        # rebuilt alpha-carbon bonds are exact -- so it introduces no impossible separations.
        model = report.models[0]
        for domain in model.domains:
            for span in domain.generated_spans():
                ca = model.ca_xyz[span.start : span.stop]
                if ca.shape[0] >= 2:
                    bonds = np.linalg.norm(np.diff(ca, axis=0), axis=1)
                    assert np.max(np.abs(bonds - CA_CA_BOND_LENGTH)) < 1e-6

    def test_unknown_engine_is_refused(self) -> None:
        with pytest.raises(ValueError, match="Unknown engine"):
            rebuild(DNMT3A, engine="nonsense", seed=0, progress=False)


#: (introduced-clash ceiling, strained-seam ceiling) per fixture, engine=walk, seed 0. A ratchet
#: for :class:`TestBackboneBaseline` -- these only ever move DOWN as the backbone redo improves.
_BACKBONE_BASELINE: dict[str, tuple[int, int]] = {
    "dnmt3a": (2, 4),
    "arf19": (0, 6),
    "p300": (4, 10),
}


class TestBackboneBaseline:
    """Frozen 2026-08 baseline for ``--backbone`` quality on the committed corpus (BB-0 floor).

    The analytic backbone currently leaves introduced steric clashes and strained seams. The
    ceilings below are the measured baseline and are **ratchets**: they may only ever move DOWN as
    the backbone redo (BB-1 clashes, BB-2 seams) improves them. Raising one to make a change pass
    hides a regression -- fix the change instead. Impossible contacts are a hard invariant, never a
    ratchet. This also guards BB-0's own deliverable: that strained seams are surfaced on the
    report rather than left for the validator to rediscover.
    """

    def _check(self, name: str) -> None:
        from dodo.construct.ca_backbone import SeamStrain

        source = FIXTURES / f"{name}.pdb"
        clash_ceiling, seam_ceiling = _BACKBONE_BASELINE[name]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ca = rebuild(source, seed=0, progress=False)
            bb = rebuild(source, backbone=True, seed=0, progress=False)
        ca_model, bb_model = ca.models[0], bb.models[0]

        # Hard invariant, never a ratchet: the backbone introduces nothing physically impossible.
        assert not find_impossible_pairs(bb_model), (
            f"{name}: --backbone introduced an impossible contact"
        )

        # Introduced steric clashes must not exceed the frozen baseline (ratchet -- down only).
        introduced = len(validate_clashes(bb_model).violations) - len(
            validate_clashes(ca_model).violations
        )
        assert introduced <= clash_ceiling, (
            f"{name}: {introduced} introduced clashes exceeds the baseline ceiling "
            f"{clash_ceiling}. This ratchet only moves down; fix the change, do not raise it."
        )

        # The rebuild introduces ZERO bond defects. The unavoidable seam bonds are reclassified as
        # kind="seam"/provenance="seam" (approximate by construction, inherited from the
        # rigidly-repositioned folded neighbour), so nothing is attributed to "rebuilt" -- the same
        # differential-clean bar the CA-only path already meets.
        bond_report = validate_bonds(bb_model)
        assert not bond_report.of_provenance("rebuilt"), (
            f"{name}: rebuild introduced a bond defect: "
            + "; ".join(v.message for v in bond_report.of_provenance("rebuilt"))
        )
        # The seams the validator independently sees (kind="seam") match the seams the report
        # records; both are within baseline and typed. CA-only placement strains no seam.
        assert len(bond_report.of_kind("seam")) == len(bb.backbone_seams)
        assert all(isinstance(s, SeamStrain) for s in bb.backbone_seams)
        assert len(bb.backbone_seams) <= seam_ceiling
        assert not ca.backbone_seams

    @pytest.mark.parametrize("name", ["dnmt3a", "arf19"])
    def test_backbone_quality_is_within_the_frozen_baseline(self, name: str) -> None:
        self._check(name)

    def test_seam_reclassification_needs_the_generated_boundary(self, tmp_path: Path) -> None:
        """The 'seam' exemption must not mask a real chain break in non-DODO input.

        The rebuilt structure labels its unclosable seams kind='seam' (not counted against ok). The
        IDENTICAL geometry, written out and read back as a plain PDB with no region assignments,
        must instead be chain_breaks -- proving the exemption keys strictly on the
        generated<->input boundary and cannot silently accept a genuine break.
        """
        import dodo

        source = FIXTURES / "dnmt3a.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(source, seed=0, backbone=True)
        with_regions = validate_bonds(report.models[0])
        seams = with_regions.of_kind("seam")
        assert seams, "expected the rebuilt structure to carry labelled seams"

        out = tmp_path / "backbone.pdb"
        dodo.write_pdb(report.models, out)
        bare = read_structure(out)  # a fresh read carries no generated spans
        without = validate_bonds(bare)
        assert not without.of_kind("seam"), "a region-less structure must not get the exemption"
        assert len(without.of_kind("chain_break")) >= len(seams), (
            "the same long C-N bonds must read back as honest chain_breaks without region info"
        )

    def test_the_joint_clash_polish_earns_its_place(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Stubbing out the coupled-clash polish must make the output measurably worse.

        Guards against the polish silently no-op'ing: single-azimuth refinement alone leaves 9
        introduced clashes on dnmt3a; the joint polish takes that to 2.
        """
        import dodo.construct.ca_backbone as ca_backbone

        source = FIXTURES / "dnmt3a.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with monkeypatch.context() as patch:
                patch.setattr(ca_backbone, "_polish_coupled_clashes", lambda *a, **k: 0)
                unpolished = rebuild(source, backbone=True, seed=0, progress=False).models[0]
            polished = rebuild(source, backbone=True, seed=0, progress=False).models[0]
        without = len(validate_clashes(unpolished).violations)
        with_polish = len(validate_clashes(polished).violations)
        assert with_polish < without, f"polish did not reduce clashes: {without} -> {with_polish}"

    @pytest.mark.slow
    def test_backbone_quality_is_within_the_frozen_baseline_p300(self) -> None:
        self._check("p300")
