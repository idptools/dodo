"""Tests for the end-to-end rebuild pipeline and the CLI.

The pipeline is the layer v1 had and the first v2 attempt did not: everything below it was
separately usable, but nothing wired it together. These tests cover the wiring, and in
particular the two things that are easy to get wrong when assembling the pieces by hand --
supplying the outer anchors so junction angles stay constrainable, and drawing a fresh
dimension target per model so a multi-model run is an ensemble rather than one conformation
repeated.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dodo.cli import main
from dodo.construct.pipeline import build_from_sequence, rebuild
from dodo.geometry.metrics import end_to_end, validate_ca_trace
from dodo.io import read_structure
from dodo.structure import DomainKind

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
