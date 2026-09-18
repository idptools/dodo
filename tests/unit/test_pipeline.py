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
from dodo.constants import (
    C_O_BOND_LENGTH,
    CA_C_BOND_LENGTH,
    CA_CLASH_DISTANCE,
    N_CA_BOND_LENGTH,
)
from dodo.construct.pipeline import build_from_sequence, rebuild
from dodo.geometry.metrics import end_to_end, validate_ca_trace
from dodo.io import read_structure
from dodo.structure import DomainKind, Span, Structure
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
        report = rebuild(DNMT3A, seed=0, backbone=False)
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

    def test_an_unknown_engine_is_refused(self) -> None:
        with pytest.raises(ValueError, match="Use 'walk'"):
            build_from_sequence(IDR_SEQUENCE, engine="nonesuch")

    def test_ensemble_shares_a_centre_not_a_pinned_atom(self) -> None:
        """Models are centred on their centroids, not nailed together at one atom.

        The walk builds every free conformer in its own frame with residue 0 at the origin.
        Shipped like that, every model of a multi-model build had its first CA at exactly
        (0, 0, 0) and the ensemble fanned out from that one pinned point.
        """
        report = build_from_sequence(IDR_SEQUENCE, n_models=6, seed=7)
        assert report.n_built == 6
        first_cas = np.array([m.ca_xyz[0] for m in report.models])
        assert not np.allclose(first_cas, first_cas[0]), "first CA is pinned across models"
        for m in report.models:
            assert np.allclose(m.ca_xyz.mean(axis=0), 0.0, atol=1e-9)


class TestFullyDisorderedRebuildStaysPut:
    """A chain with no folded domain is rebuilt with no anchors, in the engine's own frame.

    The pipeline must land each model back on the input chain's centroid. Shipped without
    that, every model came out with its first CA at the engine's origin -- pinning the
    ensemble to (0, 0, 0) and abandoning wherever the input actually placed the chain.
    """

    def test_every_model_lands_on_the_input_centroid(self) -> None:
        source = build_from_sequence(IDR_SEQUENCE, seed=1, backbone=False).models[0]
        # Move the input well away from the origin so "kept where the input was" is
        # distinguishable from "left in the engine's origin frame".
        source.xyz += np.array([100.0, 50.0, -30.0])
        input_centroid = source.ca_xyz.mean(axis=0)

        report = rebuild(source, n_models=3, seed=2, progress=False)
        assert report.n_built == 3
        first_cas = np.array([m.ca_xyz[0] for m in report.models])
        assert not np.allclose(first_cas, first_cas[0]), "first CA is pinned across models"
        for m in report.models:
            assert np.allclose(m.ca_xyz.mean(axis=0), input_centroid, atol=1e-6)


def _two_chain_structure(centre_a: np.ndarray, centre_b: np.ndarray) -> Structure:
    """Build a synthetic two-chain input for the anchor-free placement tests.

    Chain A: a 60-residue helical CA globule at ``centre_a``, preset as folded.
    Chain B: a 54-residue CA ring centred on ``centre_b``, preset as one anchor-free IDR.
    """
    from dodo.regions import assign_regions_from_spec

    def helix(n: int, centre: np.ndarray) -> np.ndarray:
        t = np.arange(n)
        coords = np.stack([4.0 * np.cos(t * 1.75), 4.0 * np.sin(t * 1.75), t * 1.5], axis=1)
        return coords - coords.mean(axis=0) + centre

    def ring(n: int, centre: np.ndarray) -> np.ndarray:
        angles = np.arange(n) * 2 * np.pi / n
        coords = 32.7 * np.stack([np.cos(angles), np.sin(angles), np.zeros(n)], axis=1)
        return coords - coords.mean(axis=0) + centre

    n_a, n_b = 60, 54
    structure = Structure.from_atom_records(
        xyz=np.vstack([helix(n_a, centre_a), ring(n_b, centre_b)]),
        atom_name=["CA"] * (n_a + n_b),
        element=["C"] * (n_a + n_b),
        residue_name=["ALA"] * n_a + ["SER"] * n_b,
        residue_number=list(range(1, n_a + 1)) + list(range(1, n_b + 1)),
        chain_id=["A"] * n_a + ["B"] * n_b,
        source="synthetic two-chain",
    )
    assign_regions_from_spec(structure, {"A": [("folded", 1, n_a)], "B": [("idr", 1, n_b)]})
    return structure


class TestFreeRegionPlacementRespectsObstacles:
    """The clash guarantee for an anchor-free region is enforced at its FINAL position.

    The region is generated in the engine's own frame and landed on its input centroid
    afterwards. First shipped without that, the region was clash-checked at the world origin and
    translated away from everything the check saw: a disordered chain could be written
    straight through a folded partner while the report said built and ok.
    """

    def _min_inter_chain(self, model: Structure) -> float:
        a = model.ca_xyz[model.chains[0].span.slice]
        b = model.ca_xyz[model.chains[1].span.slice]
        return float(np.linalg.norm(a[:, None] - b[None], axis=-1).min())

    def test_built_regions_are_clash_free_where_they_land(self) -> None:
        from dodo.constants import CA_CLASH_DISTANCE

        # Both chains' input centroids coincide, so a landing that ignored the folded chain
        # would overlap it. Every model must either land clear or refuse to build.
        centre = np.array([50.0, 50.0, 50.0])
        structure = _two_chain_structure(centre, centre)
        report = rebuild(
            structure, strategy="preset", n_models=3, seed=3, backbone=False, progress=False
        )
        outcomes = [o for o in report.outcomes if o.chain_id == "B"]
        for model, outcome in zip(report.models, outcomes, strict=True):
            if outcome.built:
                assert self._min_inter_chain(model) >= CA_CLASH_DISTANCE
            else:
                b = model.chains[1].span.slice
                assert np.allclose(model.ca_xyz[b], structure.ca_xyz[b]), (
                    "an unbuilt region must keep its input coordinates"
                )

    def test_an_obstacle_at_the_world_origin_does_not_abort_the_build(self) -> None:
        from dodo.constants import CA_CLASH_DISTANCE

        # Real AlphaFold files have atoms within the clash distance of the world origin. That
        # must not matter to a region whose input position is 120 A away: the engine's
        # origin-frame guard fired here when world-frame obstacles were (wrongly) passed in.
        structure = _two_chain_structure(np.zeros(3), np.array([120.0, 0.0, 0.0]))
        # Put a chain-A atom exactly on the origin, the worst case for the old guard.
        shift = -structure.xyz[np.argmin(np.linalg.norm(structure.xyz[:60], axis=1))]
        structure.xyz[:60] += shift
        assert float(np.linalg.norm(structure.xyz[:60], axis=1).min()) < 1e-9

        report = rebuild(
            structure, strategy="preset", n_models=1, seed=0, backbone=False, progress=False
        )
        outcome = next(o for o in report.outcomes if o.chain_id == "B")
        assert outcome.built, outcome.reason
        assert self._min_inter_chain(report.models[0]) >= CA_CLASH_DISTANCE

    def test_a_buried_centroid_is_refused_not_silently_clashed(self) -> None:
        from dodo.regions import assign_regions_from_spec

        # A disordered chain whose input centroid sits inside a dense folded blob cannot land
        # clear anywhere near it. The honest outcome is NOT BUILT with the input kept, never
        # a model that threads the blob.
        centre = np.array([50.0, 50.0, 50.0])
        grid = np.stack(np.meshgrid(*[np.arange(4) * 3.9] * 3), axis=-1).reshape(-1, 3)
        blob = grid - grid.mean(axis=0) + centre
        angles = np.arange(24) * 2 * np.pi / 24
        ring = 14.6 * np.stack([np.cos(angles), np.sin(angles), np.zeros(24)], axis=1)
        ring = ring - ring.mean(axis=0) + centre
        n_a, n_b = len(blob), len(ring)
        structure = Structure.from_atom_records(
            xyz=np.vstack([blob, ring]),
            atom_name=["CA"] * (n_a + n_b),
            element=["C"] * (n_a + n_b),
            residue_name=["ALA"] * n_a + ["SER"] * n_b,
            residue_number=list(range(1, n_a + 1)) + list(range(1, n_b + 1)),
            chain_id=["A"] * n_a + ["B"] * n_b,
            source="synthetic buried centroid",
        )
        assign_regions_from_spec(structure, {"A": [("folded", 1, n_a)], "B": [("idr", 1, n_b)]})

        report = rebuild(
            structure, strategy="preset", n_models=1, seed=0, backbone=False, progress=False
        )
        outcome = next(o for o in report.outcomes if o.chain_id == "B")
        assert not outcome.built
        assert "clashed" in (outcome.reason or "")
        b = report.models[0].chains[1].span.slice
        assert np.allclose(report.models[0].ca_xyz[b], structure.ca_xyz[b])


def _min_built_to_kept_input(model: Structure, report: object) -> float:
    """Closest approach between DODO-generated alpha carbons and geometry it left as input.

    The guarantee under test: a region DODO built must clear every region that kept its input
    coordinates by at least CA_CLASH_DISTANCE. Hydrogens are excluded, matching the rule
    ``validate_clashes`` applies. Returns ``inf`` when the model has nothing of one kind.
    """
    from scipy.spatial import cKDTree

    generated: list[np.ndarray] = []
    kept: list[np.ndarray] = []
    heavy = ~np.isin(np.char.upper(model.element.astype("<U2")), ["H", "D"])
    for domain in model.domains:
        atoms = model.atom_slice_for_residues(domain.span.start, domain.span.stop)
        if domain.kind is DomainKind.IDR and domain.rebuilt:
            generated.append(model.ca_xyz[domain.span.slice])
        elif domain.kind is DomainKind.IDR:
            keep = np.flatnonzero(heavy[atoms]) + atoms.start
            if keep.size:
                kept.append(model.xyz[keep])
    if not generated or not kept:
        return float("inf")
    tree = cKDTree(np.vstack(kept))
    return float(min(tree.query(block)[0].min() for block in generated))


class TestForwardObstacleVisibility:
    """A region built early must still avoid regions whose coordinates stay as they arrived.

    The obstacle set can only hold geometry that is already final, so a region built early
    cannot see one built later. That is harmless when the later region is itself rebuilt -- it
    avoided this one. It is a defect when the later region keeps its input coordinates, because
    then nothing ever moves out of the way.
    """

    def test_a_short_region_is_an_obstacle_before_anything_is_built(self) -> None:
        # A region under the length minimum is skipped, so its input coordinates are final
        # before the first build starts. It must therefore be in the obstacle set of every
        # region, including ones built ahead of it in the loop.
        import dodo.construct.pipeline as pipeline
        from dodo.regions import assign_regions

        source = read_structure(DNMT3A)
        assign_regions(source)
        idr_spans = [
            d.span for chain in source.chains for d in chain.domains if d.kind is DomainKind.IDR
        ]
        assert idr_spans, "fixture must have IDRs for this test to mean anything"

        original = pipeline._build_region
        seen: list[tuple[Span, np.ndarray]] = []

        def spy(structure: Structure, **kwargs: object) -> object:
            span = kwargs["span"]
            seen.append((span, structure.placed_atom_mask().copy()))  # type: ignore[arg-type]
            return original(structure, **kwargs)  # type: ignore[arg-type]

        monkeypatched = pipeline._build_region
        pipeline._build_region = spy  # type: ignore[assignment]
        try:
            # min_length high enough that the SHORTEST IDR is skipped while a longer one builds.
            lengths = sorted(len(s) for s in idr_spans)
            if len(lengths) < 2 or lengths[0] == lengths[-1]:
                pytest.skip("fixture needs IDRs of differing lengths")
            rebuild(
                DNMT3A,
                n_models=1,
                seed=0,
                backbone=False,
                progress=False,
                min_length=lengths[0] + 1,
            )
        finally:
            pipeline._build_region = monkeypatched  # type: ignore[assignment]

        short = [s for s in idr_spans if len(s) == lengths[0]]
        assert short, "no short region to check"
        # Every region's obstacle mask -- including regions built before the short one's own
        # turn -- must already cover the short region's atoms.
        for span, mask in seen:
            if span in short:
                continue
            for short_span in short:
                atoms = source.atom_slice_for_residues(short_span.start, short_span.stop)
                assert mask[atoms].all(), (
                    f"region {span.start}-{span.stop} was built without the skipped region "
                    f"{short_span.start}-{short_span.stop} in its obstacle set"
                )

    @pytest.mark.slow
    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_a_failed_region_does_not_keep_an_earlier_build_on_top_of_it(self, seed: int) -> None:
        # testing_translation.pdb has a 280-residue terminal IDR whose build fails on a broken
        # input chain, and a connecting IDR that is built BEFORE it. Without the repair pass the
        # connecting IDR landed 1.37-2.79 A from it, in output the report called built.
        report = rebuild(
            FIXTURES / "testing_translation.pdb",
            n_models=1,
            seed=seed,
            backbone=False,
            progress=False,
        )
        assert report.failures, "fixture must have a failed region for this to mean anything"
        closest = _min_built_to_kept_input(report.models[0], report)
        assert closest >= CA_CLASH_DISTANCE, (
            f"a built region sits {closest:.2f} A from a region that kept its input "
            f"coordinates (limit {CA_CLASH_DISTANCE})"
        )

    def test_a_clean_structure_pays_for_no_retries(self) -> None:
        # The repair pass must cost nothing when there is nothing to repair: one build attempt
        # per region, no second pass.
        import dodo.construct.pipeline as pipeline

        original = pipeline._build_region
        calls: list[Span] = []

        def spy(structure: Structure, **kwargs: object) -> object:
            calls.append(kwargs["span"])  # type: ignore[arg-type]
            return original(structure, **kwargs)  # type: ignore[arg-type]

        pipeline._build_region = spy  # type: ignore[assignment]
        try:
            report = rebuild(DNMT3A, n_models=1, seed=0, backbone=False, progress=False)
        finally:
            pipeline._build_region = original  # type: ignore[assignment]

        assert not report.failures, "fixture must build cleanly for this test to mean anything"
        assert len(calls) == len(set(calls)), (
            f"a region was built twice with nothing to repair: "
            f"{[s for s in calls if calls.count(s) > 1]}"
        )


class TestRepairPassMechanism:
    """The repair pass is kind-agnostic: loops and IDRs go through one build-ordered list.

    That matters because loops are built before every IDR, so a rebuilt loop can end up on the
    input coordinates of an IDR that fails later -- the same exposure, across the kind boundary.
    These drive the mechanism directly, since no committed fixture happens to pair a rebuilt loop
    with a region that kept its input.
    """

    def _pair(self, separation: float) -> tuple[Structure, list[object], list[object]]:
        """Two 2-residue regions ``separation`` apart: the first built, the second kept as input."""
        from dodo.construct.pipeline import RegionOutcome, _RegionAttempt

        xyz = np.array(
            [[0.0, 0.0, 0.0], [3.81, 0.0, 0.0], [0.0, separation, 0.0], [3.81, separation, 0.0]]
        )
        structure = Structure.from_atom_records(
            xyz=xyz,
            atom_name=["CA"] * 4,
            element=["C"] * 4,
            residue_name=["GLY"] * 4,
            residue_number=[1, 2, 3, 4],
            chain_id=["A"] * 4,
            source="synthetic repair pair",
        )
        outcomes = [
            RegionOutcome(model=1, chain_id="A", residues=(1, 2), n_residues=2, built=True),
            RegionOutcome(
                model=1,
                chain_id="A",
                residues=(3, 4),
                n_residues=2,
                built=False,
                reason="engine failed",
            ),
        ]
        retry_outcome = RegionOutcome(
            model=1, chain_id="A", residues=(1, 2), n_residues=2, built=True, reason="retried"
        )
        calls: list[str] = []

        def run_built() -> object:
            calls.append("retry")
            return retry_outcome

        attempts = [
            _RegionAttempt(
                span=Span(0, 2),
                run=run_built,  # type: ignore[arg-type]
                record=lambda _outcome: calls.append("record"),
                built=True,
                outcome_index=0,
            ),
            _RegionAttempt(
                span=Span(2, 4),
                run=lambda: outcomes[1],  # type: ignore[arg-type,return-value]
                record=lambda _outcome: None,
                built=False,
                outcome_index=1,
            ),
        ]
        return structure, attempts, [outcomes, calls]  # type: ignore[list-item]

    def test_a_contact_with_a_later_kept_region_triggers_one_retry(self) -> None:
        from dodo.construct.pipeline import _repair_forward_contacts

        structure, attempts, (outcomes, calls) = self._pair(separation=1.0)  # type: ignore[misc]
        _repair_forward_contacts(structure, attempts=attempts, outcomes=outcomes)  # type: ignore[arg-type]
        assert calls == ["retry", "record"], "the region was not rebuilt exactly once"
        assert outcomes[0].reason == "retried", "the retry's outcome did not replace the original"

    def test_a_clear_separation_triggers_nothing(self) -> None:
        from dodo.construct.pipeline import _repair_forward_contacts

        structure, attempts, (outcomes, calls) = self._pair(separation=25.0)  # type: ignore[misc]
        _repair_forward_contacts(structure, attempts=attempts, outcomes=outcomes)  # type: ignore[arg-type]
        assert calls == [], "a region well clear of everything was rebuilt anyway"
        assert outcomes[0].reason is None

    def test_a_failed_retry_keeps_the_conformer_and_discloses(self) -> None:
        from dodo.construct.pipeline import RegionOutcome, _repair_forward_contacts

        structure, attempts, (outcomes, _calls) = self._pair(separation=1.0)  # type: ignore[misc]
        attempts[0].run = lambda: RegionOutcome(  # type: ignore[assignment]
            model=1,
            chain_id="A",
            residues=(1, 2),
            n_residues=2,
            built=False,
            reason="no candidate survived",
        )
        _repair_forward_contacts(structure, attempts=attempts, outcomes=outcomes)  # type: ignore[arg-type]
        kept = outcomes[0]
        assert kept.built, "a failed retry must not mark the region unbuilt"
        assert "could not avoid" in (kept.reason or ""), kept.reason
        assert "3-4" in (kept.reason or ""), "the disclosure must name the region it sits on"

    def test_scheduler_gives_loops_to_the_repair_pass(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        import dodo.construct.pipeline as pipeline
        from dodo.regions import assign_regions

        structure = read_structure(DNMT3A)
        assign_regions(structure)
        loops = [loop for domain in structure.folded_domains() for loop in domain.loops]
        assert loops, "fixture must contain a rebuildable loop"
        seen: list[Span] = []

        def fake_build(
            structure: Structure, *, model: int, span: Span, **_kwargs: object
        ) -> object:
            residue = span.start
            chain = structure.chains[int(structure.chain_index[residue])]
            return pipeline.RegionOutcome(
                model=model,
                chain_id=chain.chain_id,
                residues=(
                    int(structure.residue_number[span.start]),
                    int(structure.residue_number[span.stop - 1]),
                ),
                n_residues=len(span),
                built=False,
                reason="probe",
            )

        def capture(
            _structure: Structure,
            *,
            attempts: list[object],
            outcomes: list[object],
        ) -> None:
            del outcomes
            seen.extend(attempt.span for attempt in attempts)  # type: ignore[attr-defined]

        monkeypatch.setattr(pipeline, "_build_region", fake_build)
        monkeypatch.setattr(pipeline, "_repair_forward_contacts", capture)
        pipeline._rebuild_one_model(
            structure,
            model=1,
            mode="predicted",
            engine=object(),
            rng=np.random.default_rng(0),
            min_length=structure.n_residues + 1,
            model_targets={},
        )
        assert all(loop in seen for loop in loops)


class TestReportDistinguishesItsNotes:
    """A built region can carry a note for two unrelated reasons. The summary must not conflate.

    ``summary()`` used to bucket every built-with-a-reason outcome under "needed a relaxed anchor
    exemption". When the repair pass began attaching a created-contact disclosure to a built
    outcome, that heading started stating something false about the one thing the disclosure
    exists to say.
    """

    def _outcome(self, **kwargs: object) -> object:
        from dodo.construct.pipeline import RegionOutcome

        base: dict[str, object] = {
            "model": 1,
            "chain_id": "A",
            "residues": (1, 40),
            "n_residues": 40,
            "built": True,
        }
        base.update(kwargs)
        return RegionOutcome(**base)  # type: ignore[arg-type]

    def test_a_created_contact_is_not_filed_as_an_anchor_exemption(self) -> None:
        from dodo.construct.pipeline import RebuildReport

        report = RebuildReport(
            outcomes=[
                self._outcome(unresolved_contact=True, reason="sits on kept input coordinates")
            ]
        )
        text = report.summary()
        assert "relaxed anchor exemption" not in text, text
        assert "closer than the clash distance to input coordinates that were kept" in text
        assert report.unresolved_contacts

    def test_a_relaxed_build_is_still_reported_as_one(self) -> None:
        from dodo.construct.pipeline import RebuildReport

        report = RebuildReport(
            outcomes=[self._outcome(relaxed_anchors=True, reason="built against a relaxed ...")]
        )
        text = report.summary()
        assert "relaxed anchor exemption" in text
        assert not report.unresolved_contacts
        assert report.ok, "a relaxed build is a deliberate trade, not a failure"

    def test_an_unresolved_contact_makes_the_run_not_ok(self) -> None:
        # The hole this closes: the neighbouring failed region is short enough to be tolerated,
        # so nothing else made the run unsuccessful, and DODO exited 0 having written a contact
        # it created itself.
        from dodo.construct.pipeline import RebuildReport

        tolerated = self._outcome(
            residues=(41, 47), n_residues=7, built=False, reason="engine failed"
        )
        assert tolerated.tolerated, "fixture must be a tolerated failure for this to mean anything"
        report = RebuildReport(
            outcomes=[
                self._outcome(unresolved_contact=True, reason="sits on kept input coordinates"),
                tolerated,
            ]
        )
        assert not report.blocking_failures, "only the contact should make this run not-ok"
        assert not report.ok

    def test_a_clashing_folded_unit_makes_the_run_not_ok(self) -> None:
        from dodo.construct.pipeline import RebuildReport
        from dodo.construct.place import DomainPlacement, UnitPlacement

        report = RebuildReport(
            models=[object()],  # type: ignore[list-item]
            placements=[
                DomainPlacement(chain_id="A", residues=(1, 10), moved=True, clashing=True)
            ],
            units=[
                UnitPlacement(
                    index=0,
                    chains=("A",),
                    domains=("A:1-10",),
                    n_atoms=10,
                    moved=True,
                    clashing=True,
                    clashing_atoms=2,
                )
            ],
        )
        assert not report.ok
        from dodo.cli import _report

        assert _report(report, quiet=True) == 2


class TestEnsembleTopologyReconciliation:
    """Mixed built/failed regions across models must not make the ensemble unwritable.

    A region that builds in some models and fails in others gives frames with different atom
    counts, which cannot be written into one multi-model file. rebuild() must return a writable
    majority and disclose the dropped models -- not crash the write with nothing saved.
    """

    def _model_with(self, n_atoms: int, first_atom: str = "CA") -> Structure:
        return Structure.from_atom_records(
            xyz=np.arange(n_atoms * 3, dtype=float).reshape(n_atoms, 3),
            atom_name=[first_atom] + ["CA"] * (n_atoms - 1),
            element=["C"] * n_atoms,
            residue_name=["GLY"] * n_atoms,
            residue_number=list(range(1, n_atoms + 1)),
            chain_id=["A"] * n_atoms,
            source="synthetic",
        )

    def test_the_minority_topology_is_dropped_and_disclosed(self) -> None:
        from dodo.construct.pipeline import RebuildReport, _reconcile_ensemble_topology

        report = RebuildReport(
            models=[self._model_with(10), self._model_with(12), self._model_with(10)],
            model_numbers=[1, 2, 3],
        )
        _reconcile_ensemble_topology(report)
        assert [m.n_atoms for m in report.models] == [10, 10]
        assert report.model_numbers == [1, 3]
        assert any("model(s) 2" in note for note in report.notes)
        assert any("MODEL 2 = build model 3" in note for note in report.notes)

    def test_equal_atom_counts_with_different_records_still_diverge(self) -> None:
        # Two identical anchor-free chains can swap which one fails between models: the
        # frames then have EQUAL atom counts but different atom records, which the writers
        # reject just the same. Reconciliation must group on the writers' full rule.
        from dodo.construct.pipeline import RebuildReport, _reconcile_ensemble_topology

        report = RebuildReport(
            models=[
                self._model_with(10),
                self._model_with(10, first_atom="N"),
                self._model_with(10),
            ],
            model_numbers=[1, 2, 3],
        )
        _reconcile_ensemble_topology(report)
        assert report.model_numbers == [1, 3]
        assert any("model(s) 2" in note for note in report.notes)

    def test_the_kept_group_is_what_the_writer_accepts(self) -> None:
        from dodo.construct.pipeline import RebuildReport, _reconcile_ensemble_topology
        from dodo.io.write import _require_matching_topology

        report = RebuildReport(
            models=[
                self._model_with(10),
                self._model_with(12),
                self._model_with(10, first_atom="N"),
                self._model_with(10),
            ],
            model_numbers=[1, 2, 3, 4],
        )
        _reconcile_ensemble_topology(report)
        _require_matching_topology(report.models)  # must not raise
        assert report.model_numbers == [1, 4]

    def test_a_tie_keeps_the_earlier_models(self) -> None:
        from dodo.construct.pipeline import RebuildReport, _reconcile_ensemble_topology

        report = RebuildReport(
            models=[self._model_with(12), self._model_with(10)], model_numbers=[1, 2]
        )
        _reconcile_ensemble_topology(report)
        assert [m.n_atoms for m in report.models] == [12]
        assert report.model_numbers == [1]
        assert any("model(s) 2" in note for note in report.notes)

    def test_a_consistent_ensemble_is_untouched(self) -> None:
        from dodo.construct.pipeline import RebuildReport, _reconcile_ensemble_topology

        report = RebuildReport(
            models=[self._model_with(10), self._model_with(10)], model_numbers=[1, 2]
        )
        _reconcile_ensemble_topology(report)
        assert len(report.models) == 2
        assert report.model_numbers == [1, 2]
        assert not report.notes


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

    @staticmethod
    def _ca_only_pdb(path: Path, *, first_residue_number: int, n: int = 30) -> Path:
        """Write a straight CA-only chain under an author numbering of our choosing."""
        lines = [
            f"ATOM  {i + 1:>5}  CA  ALA A{first_residue_number + i:>4}    "
            f"{i * 3.81:>8.3f}{0.0:>8.3f}{0.0:>8.3f}  1.00 50.00           C"
            for i in range(n)
        ]
        path.write_text("\n".join(lines) + "\nEND\n")
        return path

    def test_regions_reports_the_files_own_numbering(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """A chain numbered from 100 must be reported as 100-129, not 1-30."""
        pdb = self._ca_only_pdb(tmp_path / "from100.pdb", first_residue_number=100)
        assert main(["regions", str(pdb)]) == 0
        assert capsys.readouterr().out.strip() == "chain A: IDR 100-129"

    def test_regions_scores_share_the_axis_of_the_region_line(
        self, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """The score profile is the evidence for the boundaries, so it uses their numbering.

        A profile labelled 1..30 against a region line reading 100-129 would be an audit
        trail a reader cannot line up against either the region line or the input file.
        """
        pdb = self._ca_only_pdb(tmp_path / "from100.pdb", first_residue_number=100)
        assert main(["regions", str(pdb), "--scores"]) == 0
        lines = capsys.readouterr().out.splitlines()

        assert lines[0] == "chain A: IDR 100-129"
        labels = [line.split("\t")[0].strip() for line in lines[1:]]
        assert labels == [str(n) for n in range(100, 130)]

    def test_regions_exits_one_on_an_unreadable_file(self, tmp_path: Path) -> None:
        """`regions` is scriptable, so a read failure must not look like a clean run."""
        assert main(["regions", str(tmp_path / "does-not-exist.pdb")]) == 1

    def test_units_uses_experimental_rebuild_preprocessing(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        path = FIXTURES / "6kn7.pdb"
        assert main(["units", str(path), "--units", "experimental", "-q"]) == 0
        assert "over 29 folded domain(s)" in capsys.readouterr().out

    def test_units_accepts_the_same_fasta_as_rebuild(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        data = FIXTURES.parent / "complexes_missing_residues"
        assert (
            main(
                [
                    "units",
                    str(data / "7R5J_subset.cif"),
                    "--fasta",
                    str(data / "7R5J_chain_sequences.fasta"),
                    "-q",
                ]
            )
            == 0
        )
        captured = capsys.readouterr()
        assert "over 6 folded domain(s)" in captured.out
        assert "2686 inserted residue(s)" in captured.err

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
    """The ``backbone=`` flag on the two entry points, and its default.

    On by default: the backbone is the point of a rebuild for most callers. ``backbone=False`` opts
    back out to alpha carbons only. Both modes use the same seam-compatible CA trace, so the flag
    adds atoms without quietly returning a different conformation.
    """

    def test_on_by_default(self) -> None:
        report = build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", seed=0)
        assert {str(n) for n in report.models[0].atom_name} == {"N", "CA", "C", "O"}

    def test_off_when_opted_out(self) -> None:
        report = build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", seed=0, backbone=False)
        assert {str(n) for n in report.models[0].atom_name} == {"CA"}

    def test_on_when_asked(self) -> None:
        report = build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", seed=0, backbone=True)
        assert {str(n) for n in report.models[0].atom_name} == {"N", "CA", "C", "O"}

    def test_alpha_carbons_are_identical_either_way(self) -> None:
        """The flag adds atoms; it must not change DODO's actual answer.

        Same seed, so the alpha carbons are the same coordinates to the bit. If this drifts, the
        backbone pass is perturbing the trace rather than decorating it.
        """
        plain = build_from_sequence("GRNQNGGGYQNYNNQGYQGHGG", seed=0, backbone=False).models[0]
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

    def test_every_generated_bond_is_exact_including_seams(self) -> None:
        """Boundary-aware CA generation closes the former folded-domain gaps exactly."""
        source = FIXTURES / "dnmt3a.pdb"
        for seed in (0, 1, 2):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                report = rebuild(source, seed=seed, backbone=True)
            bonds = validate_bonds(report.models[0])
            assert not report.backbone_seams
            assert not bonds.of_kind("seam")
            assert not bonds.of_provenance("rebuilt")

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
    """Closing seams must never come at the cost of damaging folded-domain geometry.

    A rejected shortcut kept the input boundary residue and replaced its nitrogen. That drove the
    nitrogen into its existing side chain and broke a proline ring. The production fix instead
    constrains generated alpha carbons and leaves every folded-domain atom untouched.
    """

    def test_no_new_violation_inside_a_single_residue(self) -> None:
        """Catch the rejected repair that damaged PRO282 and GLU473 internally.

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

    def test_no_bond_violation_is_introduced(self) -> None:
        """The exact seam repair leaves neither seam nor internal bond violations."""
        source = FIXTURES / "dnmt3a.pdb"
        baseline = {v.message for v in validate_bonds(read_structure(source)).violations}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            model = rebuild(source, seed=0, backbone=True).models[0]
        introduced = [v for v in validate_bonds(model).violations if v.message not in baseline]
        assert introduced == []

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


_BACKBONE_BASELINE: dict[str, tuple[int, int]] = {
    "dnmt3a": (2, 0),
    "arf19": (0, 0),
    "p300": (3, 0),  # clashes 4 -> 3: finer coupled-clash azimuth grid (5 deg vs 15 deg)
}


class TestEndToEndToleranceIsDisclosed:
    """The 10% end-to-end allowance must be visible where it applies -- and only there.

    It applies to a *steered* region: one with a free end, whose span the walk actually aims at a
    target. It does NOT apply to an interior region, whose span is dictated by its two fixed
    anchors; the engine neither samples a target for one nor checks it afterwards. Comparing the
    two numbers there compares a region's own span against the separation of the anchors outside
    it, which differ by the direction of two terminal bonds -- up to 7.62 A of pure geometry.
    """

    def test_interior_regions_are_never_flagged(self) -> None:
        """The bug this guards: scoring interior regions flagged 7 of 7, none of them real.

        dnmt3a seed 2 is the sharpest case -- residues 433-473 span 52.2 A between anchors that
        are 47.6 A apart, which is a correct closure, not a 9.5% miss.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(DNMT3A, seed=2, progress=False)
        interior = [o for o in report.outcomes if o.built and not o.steered]
        assert interior, "fixture must contain an interior region for this to mean anything"
        wide = [
            o
            for o in interior
            if abs(o.achieved_end_to_end - o.requested_end_to_end) / o.requested_end_to_end > 0.05
        ]
        assert wide, "seed 2 no longer produces a wide interior span; re-pick one"
        assert "from the requested end-to-end distance" not in report.summary()
        # ...and it is described as a span, not as a missed target.
        assert "set by its anchors" in str(wide[0])

    def test_steered_regions_hit_their_target_on_the_corpus(self) -> None:
        """Where steering applies it is accurate, so the summary stays quiet."""
        worst = 0.0
        for name in ("dnmt3a", "arf19"):
            for seed in (0, 1, 2):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    report = rebuild(FIXTURES / f"{name}.pdb", seed=seed, progress=False)
                for o in report.outcomes:
                    if o.built and o.steered and o.achieved_end_to_end and o.requested_end_to_end:
                        rel = abs(o.achieved_end_to_end - o.requested_end_to_end)
                        worst = max(worst, rel / o.requested_end_to_end)
                assert "from the requested end-to-end distance" not in report.summary()
        assert worst < 0.05, f"steered end-to-end accuracy regressed: worst {worst:.1%}"

    def test_a_steered_region_that_leans_on_the_tolerance_is_named(self) -> None:
        """The disclosure path itself, driven directly so it does not depend on a lucky seed."""
        from dodo.construct.pipeline import RebuildReport, RegionOutcome

        report = RebuildReport(
            outcomes=[
                RegionOutcome(
                    model=1,
                    chain_id="A",
                    residues=(1, 40),
                    n_residues=40,
                    built=True,
                    achieved_end_to_end=90.0,
                    requested_end_to_end=100.0,
                    steered=True,
                )
            ]
        )
        summary = report.summary()
        assert "1 region(s) finished more than 5% from the requested end-to-end distance" in summary
        assert "10.0%" in summary
        assert "1-40" in summary


class TestBackboneBaseline:
    """Frozen 2026-08 baseline for ``--backbone`` quality on the committed corpus (BB-0 floor).

    The ceilings below are measured **ratchets**: they may only move down. Seam ceilings are now
    zero; clash ceilings retain the best committed seed-0 result. Raising one to make a change pass
    hides a regression. Impossible contacts are a hard invariant, never a ratchet.
    """

    def _check(self, name: str) -> None:
        from dodo.construct.ca_backbone import SeamStrain

        source = FIXTURES / f"{name}.pdb"
        clash_ceiling, seam_ceiling = _BACKBONE_BASELINE[name]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ca = rebuild(source, seed=0, backbone=False, progress=False)
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

        # The rebuild introduces ZERO bond defects, including at folded-domain seams.
        bond_report = validate_bonds(bb_model)
        assert not bond_report.of_provenance("rebuilt"), (
            f"{name}: rebuild introduced a bond defect: "
            + "; ".join(v.message for v in bond_report.of_provenance("rebuilt"))
        )
        # Validator and report agree, and the zero seam ceiling is exact rather than aspirational.
        assert len(bond_report.of_kind("seam")) == len(bb.backbone_seams)
        assert all(isinstance(s, SeamStrain) for s in bb.backbone_seams)
        assert len(bb.backbone_seams) <= seam_ceiling
        assert not ca.backbone_seams

    @pytest.mark.parametrize("name", ["dnmt3a", "arf19"])
    def test_backbone_quality_is_within_the_frozen_baseline(self, name: str) -> None:
        self._check(name)

    def test_seam_reclassification_needs_the_generated_boundary(self, tmp_path: Path) -> None:
        """The 'seam' exemption must not mask a real chain break in non-DODO input.

        Normal fixture rebuilds now close every seam, so make one deliberately long after the
        rebuild. With region provenance it is a seam; the IDENTICAL geometry, written and read
        back without assignments, must be a chain break. This proves the exemption keys strictly
        on the generated/input boundary and cannot silently accept a genuine break.
        """
        import dodo

        source = FIXTURES / "dnmt3a.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(source, seed=0, backbone=True)
        model = report.models[0].copy()
        span = next(
            span
            for domain in model.domains
            for span in domain.generated_spans()
            if span.stop < model.n_residues
        )
        boundary = span.stop - 1
        anchor = span.stop
        c_atoms = model.atom_slice_for_residues(boundary, boundary + 1)
        n_atoms = model.atom_slice_for_residues(anchor, anchor + 1)
        c_index = next(i for i in range(c_atoms.start, c_atoms.stop) if model.atom_name[i] == "C")
        n_index = next(i for i in range(n_atoms.start, n_atoms.stop) if model.atom_name[i] == "N")
        direction = model.xyz[n_index] - model.xyz[c_index]
        model.xyz[n_index] = model.xyz[c_index] + 4.0 * direction / np.linalg.norm(direction)

        with_regions = validate_bonds(model)
        seams = with_regions.of_kind("seam")
        assert seams, "expected the rebuilt structure to carry labelled seams"

        out = tmp_path / "backbone.pdb"
        dodo.write_pdb([model], out)
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


def _rebuilt_residues(model: Structure) -> list[int]:
    """Every residue index DODO rebuilt -- the union of the generated spans."""
    return [
        residue
        for domain in model.domains
        for span in domain.generated_spans()
        for residue in range(span.start, span.stop)
    ]


def _atoms_by_name(model: Structure, residue: int) -> dict[str, np.ndarray]:
    """Return the atoms of one residue keyed by atom name."""
    atoms = model.atom_slice_for_residues(residue, residue + 1)
    names = model.atom_name[atoms]
    xyz = model.xyz[atoms]
    return {str(name): xyz[i] for i, name in enumerate(names)}


class TestBackboneIsFirstClass:
    """First-class guarantees for ``--backbone``, at EVERY seed on EVERY committed fixture.

    :class:`TestBackboneBaseline` pins the seed-0 quality floor. These guarantees must hold for
    every conformer: no impossible contact, exact rebuilt and seam bonds, a complete N/CA/C/O on
    every rebuilt residue, and reproducible output. A regression that surfaced only at seed 1
    would pass a seed-0 baseline untouched; this closes that gap.
    """

    @pytest.mark.parametrize("seed", [0, 1, 2])
    @pytest.mark.parametrize("name", ["dnmt3a", "arf19"])
    def test_first_class_invariants(self, name: str, seed: int) -> None:
        self._check(name, seed)

    @pytest.mark.slow
    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_first_class_invariants_p300(self, seed: int) -> None:
        self._check("p300", seed)

    def _check(self, name: str, seed: int) -> None:
        from dodo.construct.ca_backbone import SeamStrain

        source = FIXTURES / f"{name}.pdb"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            report = rebuild(source, seed=seed, backbone=True, progress=False)
        model = report.models[0]
        rebuilt = _rebuilt_residues(model)
        assert rebuilt, f"{name}: fixture has no rebuilt regions, so it proves nothing here"

        # 1. Nothing physically impossible beyond what the input already contained -- a hard
        #    invariant at every seed, never a ratchet.
        inherited = {
            (p.residue_labels, p.atom_names) for p in find_impossible_pairs(read_structure(source))
        }
        introduced = [
            p
            for p in find_impossible_pairs(model)
            if (p.residue_labels, p.atom_names) not in inherited
        ]
        assert introduced == [], f"{name} seed {seed}: --backbone introduced {introduced}"

        # 2. The rebuild introduces zero bond defects of its own.
        bonds = validate_bonds(model)
        assert not bonds.of_provenance("rebuilt"), (
            f"{name} seed {seed}: rebuilt-provenance bond defect: "
            + "; ".join(v.message for v in bonds.of_provenance("rebuilt"))
        )

        # 3. Every rebuilt residue carries a COMPLETE N/CA/C/O backbone, and the three bonds one
        #    residue determines are exact by construction. The cross-residue C-N bonds, including
        #    folded-domain seams, are covered by (2) and (4).
        for residue in rebuilt:
            atoms = _atoms_by_name(model, residue)
            missing = {"N", "CA", "C", "O"} - set(atoms)
            assert not missing, (
                f"{name} seed {seed}: rebuilt residue {residue} is missing {missing}"
            )
            for first, second, ideal, label in (
                ("N", "CA", N_CA_BOND_LENGTH, "N-CA"),
                ("CA", "C", CA_C_BOND_LENGTH, "CA-C"),
                ("C", "O", C_O_BOND_LENGTH, "C-O"),
            ):
                bond = float(np.linalg.norm(atoms[first] - atoms[second]))
                assert abs(bond - ideal) < 1e-6, (
                    f"{name} seed {seed} residue {residue}: {label} bond {bond:.6f} A vs {ideal}"
                )

        # 4. Boundary-aware CA generation makes every peptide seam exactly closable. Keep the
        #    validator/report agreement check too, so a future fallback cannot go unreported.
        assert len(bonds.of_kind("seam")) == len(report.backbone_seams)
        assert all(isinstance(s, SeamStrain) for s in report.backbone_seams)
        assert not report.backbone_seams, (
            f"{name} seed {seed}: boundary-aware generation left a strained seam"
        )

        # 5. Reproducible: the same seed yields byte-identical backbone atoms on a second run.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            again = rebuild(source, seed=seed, backbone=True, progress=False).models[0]
        assert np.array_equal(model.xyz, again.xyz), (
            f"{name} seed {seed}: --backbone output is not reproducible"
        )


class TestComplexes:
    """Multi-chain input, through the whole pipeline.

    The invariant is stated on the finished model rather than on the placement step, because
    that is what a user gets: two atoms DODO held rigid together are the same distance apart in
    the output as in the input, and a complex that arrives intact leaves intact.
    """

    DIMER = FIXTURES / "dnmt3a_dimer.pdb"

    @staticmethod
    def _interface_pairs(structure, cutoff: float = 5.0):
        from scipy.spatial import cKDTree

        chain_of = np.repeat(structure.chain_index, np.diff(structure.residue_atom_offsets))
        pairs = np.array(sorted(cKDTree(structure.xyz).query_pairs(cutoff)), dtype=np.int64)
        return pairs[chain_of[pairs[:, 0]] != chain_of[pairs[:, 1]]]

    def test_the_interface_survives_a_full_rebuild(self) -> None:
        from dodo.construct.assembly import find_rigid_units
        from dodo.io import read_structure
        from dodo.regions.identify import assign_regions

        original = read_structure(self.DIMER)
        probe = read_structure(self.DIMER)
        assign_regions(probe)
        assembly = find_rigid_units(probe)

        report = rebuild(self.DIMER, seed=0, n_models=1, progress=False)
        model = report.models[0]

        # Folded atoms are never regenerated, so they can be matched by identity. Loops inside a
        # folded domain ARE rebuilt, so they are excluded: they are not part of the rigid body.
        def index(structure):
            out = {}
            for atom in range(structure.n_atoms):
                residue = structure.residue_index[atom]
                out[
                    (
                        int(structure.chain_index[residue]),
                        int(structure.residue_number[residue]),
                        str(structure.atom_name[atom]),
                    )
                ] = atom
            return out

        before_index, after_index = index(original), index(model)
        worst = 0.0
        counted = 0
        for unit in assembly.units:
            if len(unit) < 2:
                continue
            keys = []
            for domain in unit.domains:
                loop_residues = {r for loop in domain.loops for r in range(loop.start, loop.stop)}
                atoms = probe.atom_slice_for_residues(domain.span.start, domain.span.stop)
                for atom in range(atoms.start, atoms.stop):
                    residue = probe.residue_index[atom]
                    if int(residue) in loop_residues:
                        continue
                    key = (
                        int(probe.chain_index[residue]),
                        int(probe.residue_number[residue]),
                        str(probe.atom_name[atom]),
                    )
                    if key in before_index and key in after_index:
                        keys.append(key)
            if len(keys) < 2:
                continue
            before = original.xyz[[before_index[k] for k in keys]]
            after = model.xyz[[after_index[k] for k in keys]]
            from scipy.spatial import cKDTree

            pairs = cKDTree(before).query_pairs(5.0, output_type="ndarray")
            if pairs.size == 0:
                continue
            d0 = np.linalg.norm(before[pairs[:, 0]] - before[pairs[:, 1]], axis=1)
            d1 = np.linalg.norm(after[pairs[:, 0]] - after[pairs[:, 1]], axis=1)
            worst = max(worst, float(np.abs(d1 - d0).max()))
            counted += len(pairs)
        assert counted > 10_000, "the fixture has no interface to check"
        assert worst < 1e-9, f"a rigid unit's internal geometry moved by {worst:.2e} A"

    def test_units_none_reproduces_the_old_behaviour(self) -> None:
        """Check the escape hatch actually escapes."""
        locked = rebuild(self.DIMER, seed=0, n_models=1, backbone=False, progress=False)
        loose = rebuild(
            self.DIMER, units="none", seed=0, n_models=1, backbone=False, progress=False
        )
        assert len(locked.rigid_units) == 1
        assert loose.rigid_units == []
        assert not np.array_equal(locked.models[0].xyz, loose.models[0].xyz)

    def test_experimental_moves_nothing(self) -> None:
        report = rebuild(
            self.DIMER, units="experimental", seed=0, n_models=1, backbone=False, progress=False
        )
        assert not [p for p in report.placements if p.moved]

    def test_predicted_repositions_linker_connected_domains(self) -> None:
        report = rebuild(
            self.DIMER, units="predicted", seed=0, n_models=1, backbone=False, progress=False
        )
        assert [p for p in report.placements if p.moved]
        # ...but not at the cost of the interface: the units holding it did not move.
        assert all(not u.moved for u in report.rigid_units)

    def test_an_unknown_units_mode_is_refused(self) -> None:
        from dodo.exceptions import InvalidParameterError

        with pytest.raises(InvalidParameterError, match="Unknown units"):
            rebuild(self.DIMER, units="everything", n_models=1, progress=False)

    def test_no_peptide_bond_is_invented_across_a_chain_break(self) -> None:
        """No peptide bond is invented across a chain break.

        A region ending at a chain's C terminus used to have its carbonyl aimed at the first
        nitrogen of the NEXT chain, producing "seams" of 71 and 149 A that were reported as
        strained peptide bonds rather than as the chain breaks they are.
        """
        report = rebuild(self.DIMER, seed=0, n_models=1, progress=False)
        assert report.backbone_seams
        assert max(seam.bond_length for seam in report.backbone_seams) < 10.0


class TestUnmodelledResidues:
    """Rebuilding what a structure never modelled, from a reference sequence."""

    @staticmethod
    def _carved():
        from dodo.io import read_structure

        full = read_structure(FIXTURES / "dnmt3a.pdb")
        drop = np.zeros(full.n_residues, dtype=bool)
        for start, stop in ((0, 40), (440, 460), (600, 612), (887, 912)):
            drop[start:stop] = True
        return full.select_residues(~drop), full.sequence

    def test_the_missing_residues_are_built(self) -> None:
        partial, reference = self._carved()
        report = rebuild(
            partial, sequences={"A": reference}, seed=0, n_models=1, progress=False
        )
        model = report.models[0]
        assert model.sequence == reference
        assert int(model.inserted.sum()) == 97
        assert report.insertions is not None
        assert report.insertions.n_inserted == 97

    def test_inserted_residues_get_real_geometry(self) -> None:
        """An inserted alpha carbon must not still sit on the line it was parked on."""
        partial, reference = self._carved()
        report = rebuild(
            partial,
            sequences={"A": reference},
            seed=0,
            n_models=1,
            backbone=False,
            progress=False,
        )
        model = report.models[0]
        ca = model.ca_xyz
        bonds = np.linalg.norm(np.diff(ca, axis=0), axis=1)
        touching = model.inserted[1:] | model.inserted[:-1]
        assert touching.any()
        assert float(bonds[touching].min()) > 3.5
        assert float(bonds[touching].max()) < 4.1

    def test_the_provenance_flag_survives_to_the_model(self) -> None:
        partial, reference = self._carved()
        report = rebuild(
            partial, sequences={"A": reference}, seed=0, n_models=1, progress=False
        )
        assert int(report.models[0].inserted.sum()) == 97

    def test_nothing_happens_without_a_reference(self) -> None:
        partial, _ = self._carved()
        report = rebuild(partial, seed=0, n_models=1, progress=False)
        assert report.insertions is None
        assert int(report.models[0].inserted.sum()) == 0

    def test_a_fasta_file_works_end_to_end(self, tmp_path: Path) -> None:
        partial, reference = self._carved()
        path = tmp_path / "reference.fasta"
        path.write_text(f">A dnmt3a\n{reference}\n")
        report = rebuild(partial, fasta=path, seed=0, n_models=1, progress=False)
        assert report.models[0].sequence == reference

    def test_fill_missing_uses_the_deposited_sequence(self) -> None:
        """A chain with SEQRES carries its own reference; fill_missing opts into using it."""
        partial, reference = self._carved()
        partial.chains[0].full_sequence = reference
        off = rebuild(partial, seed=0, n_models=1, progress=False)
        on = rebuild(partial, fill_missing=True, seed=0, n_models=1, progress=False)
        assert off.models[0].n_residues < on.models[0].n_residues
        assert on.models[0].sequence == reference


class TestInputKind:
    """``units="auto"`` has to tell a measurement from a prediction, and say which it chose."""

    def test_a_prediction_resolves_to_predicted(self) -> None:
        """AlphaFold DB models and AlphaFold 3 server output declare no experimental method."""
        report = rebuild(FIXTURES / "dnmt3a.pdb", seed=0, n_models=1, progress=False)
        note = next(n for n in report.notes if "units=auto" in n)
        assert "'predicted'" in note
        assert [p for p in report.placements if p.moved]

    def test_a_declared_experimental_method_resolves_to_experimental(self) -> None:
        """6kn7 declares ELECTRON MICROSCOPY in both formats; nothing measured may move."""
        from dodo.io import read_structure

        original = read_structure(FIXTURES / "6kn7.pdb")
        report = rebuild(FIXTURES / "6kn7.pdb", seed=0, n_models=1, backbone=False, progress=False)
        note = next(n for n in report.notes if "units=auto" in n)
        assert "'experimental'" in note and "ELECTRON MICROSCOPY" in note
        assert not [p for p in report.placements if p.moved]
        # Every folded-domain atom is exactly where the file put it. Rebuilt loops inside a
        # folded domain are excluded: those residues are regenerated by design.
        from dodo.regions.identify import assign_regions
        from dodo.structure import DomainKind

        assign_regions(original)
        model = report.models[0]
        rigid = np.zeros(original.n_atoms, dtype=bool)
        for domain in original.domains:
            if domain.kind is not DomainKind.FOLDED:
                continue
            rigid[domain.atom_slice] = True
            for loop in domain.loops:
                rigid[original.atom_slice_for_residues(loop.start, loop.stop)] = False
        kept = {
            (
                int(model.chain_index[model.residue_index[a]]),
                int(model.residue_number[model.residue_index[a]]),
                str(model.atom_name[a]),
            ): a
            for a in range(model.n_atoms)
        }
        checked = 0
        for atom in np.flatnonzero(rigid)[::97]:  # a spread sample; the full set is 61,511
            residue = original.residue_index[atom]
            key = (
                int(original.chain_index[residue]),
                int(original.residue_number[residue]),
                str(original.atom_name[atom]),
            )
            if key in kept:
                assert np.array_equal(model.xyz[kept[key]], original.xyz[atom])
                checked += 1
        assert checked > 500

    def test_filling_in_residues_resolves_to_experimental(self) -> None:
        """Only a structure that did not model everything has residues to fill in."""
        partial, reference = TestUnmodelledResidues._carved()
        report = rebuild(
            partial, sequences={"A": reference}, seed=0, n_models=1, progress=False
        )
        note = next(n for n in report.notes if "units=auto" in n)
        assert "'experimental'" in note
        assert not [p for p in report.placements if p.moved]

    def test_an_explicit_mode_overrides_the_detection(self) -> None:
        partial, reference = TestUnmodelledResidues._carved()
        report = rebuild(
            partial,
            sequences={"A": reference},
            units="predicted",
            seed=0,
            n_models=1,
            progress=False,
        )
        assert not [n for n in report.notes if "units=auto" in n]

    def test_a_theoretical_model_is_not_experimental(self) -> None:
        """``EXPDTA THEORETICAL MODEL`` is the PDB's own marker for a computed structure."""
        from dodo.structure import classify_experimental_method

        assert classify_experimental_method("THEORETICAL MODEL") is None
        assert classify_experimental_method(None) is None
        assert classify_experimental_method("?") is None
        assert classify_experimental_method("X-RAY DIFFRACTION") == "X-RAY DIFFRACTION"
        assert classify_experimental_method("  electron   microscopy ") == "ELECTRON MICROSCOPY"


class TestProgressStages:
    """Progress has to cover the whole rebuild, not just the region loop.

    The bar used to be created after the file was read and sized only for the regions, so on a
    597 MB assembly the first two minutes -- and then region identification, filling in
    unmodelled residues, finding rigid units and positioning them -- were silent. A stage that
    sets a label and then works without advancing is no better: tqdm only redraws when something
    moves, so a frozen line is exactly the thing this is meant to prevent. Hence the assertion
    is that every stage's counter *reaches its total*, not merely that a stage was announced.
    """

    class _Spy:
        """A tracker that records what each stage counted, in place of a real bar."""

        def __init__(self) -> None:
            self.stages: list[tuple[str, int, int | None]] = []
            self._label = "start"
            self._total: int | None = None
            self._n = 0

        def stage(self, label: str, total: int | None = None, unit: str = "") -> None:
            self.stages.append((self._label, self._n, self._total))
            self._label, self._total, self._n = label, total, 0

        def batched(self, batch: int):
            return self.advance

        def advance(self, amount: int = 1) -> None:
            self._n += amount

        def describe(self, text: str) -> None:
            return

        def next_model(self, done: int, total: int) -> None:
            return

        def close(self) -> None:
            self.stages.append((self._label, self._n, self._total))

    def _run(self, monkeypatch: pytest.MonkeyPatch, **kwargs: object):
        from dodo.construct import pipeline as pipeline_module

        spy = self._Spy()
        monkeypatch.setattr(pipeline_module, "_progress_bar", lambda requested: spy)
        report = rebuild(progress=True, seed=0, n_models=1, **kwargs)  # type: ignore[arg-type]
        return spy, report

    def test_every_stage_of_a_complex_rebuild_is_covered(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        spy, _ = self._run(monkeypatch, source=FIXTURES / "dnmt3a_dimer.pdb")
        labels = [label for label, _n, _total in spy.stages]
        assert any(label.startswith("reading") for label in labels)
        for expected in (
            "identifying regions",
            "finding rigid units",
            "positioning folded domains",
            "rebuilding",
        ):
            assert expected in labels, f"no progress stage for {expected!r}; labels were {labels}"

    def test_reading_is_a_stage_of_its_own(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """The single longest silent stretch on a large file, and it used to have no bar at all."""
        spy, _ = self._run(monkeypatch, source=FIXTURES / "dnmt3a_dimer.pdb")
        label, counted, _total = next(
            entry for entry in spy.stages if entry[0].startswith("reading")
        )
        assert "dnmt3a_dimer.pdb" in label
        assert counted == 14292, "the reader reported a different number of records than it read"

    def test_every_measured_stage_reaches_its_total(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        spy, _ = self._run(monkeypatch, source=FIXTURES / "dnmt3a_dimer.pdb")
        for label, counted, total in spy.stages:
            if total is None:
                continue
            assert counted == total, f"stage {label!r} stopped at {counted} of {total}"

    def test_filling_in_residues_is_its_own_stage(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        partial, reference = TestUnmodelledResidues._carved()
        spy, _ = self._run(monkeypatch, source=partial, sequences={"A": reference})
        labels = [label for label, _n, _total in spy.stages]
        assert "filling in unmodelled residues" in labels

    def test_a_structure_passed_in_memory_has_no_reading_stage(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        from dodo.io import read_structure

        structure = read_structure(FIXTURES / "dnmt3a.pdb")
        spy, _ = self._run(monkeypatch, source=structure)
        assert not [label for label, _n, _total in spy.stages if label.startswith("reading")]

    def test_progress_false_costs_nothing(self) -> None:
        """The no-op tracker must answer every call the real one does."""
        from dodo.construct.pipeline import _NoProgress, _Progress, _progress_bar

        tracker = _progress_bar(False)
        assert isinstance(tracker, _NoProgress)
        tracker.stage("x", total=3, unit="y")
        tracker.batched(10)(5)
        tracker.advance(1)
        tracker.describe("z")
        tracker.next_model(1, 2)
        tracker.close()
        assert set(_NoProgress.__slots__ or ()) == set()
        assert {name for name in dir(_Progress) if not name.startswith("_")} <= {
            name for name in dir(_NoProgress) if not name.startswith("_")
        }


class TestObservedResiduesAreStatic:
    """In an experimental structure, a resolved residue stays put -- however few there are.

    The rule, in Ryan's words: a residue resolved in the map is "sufficiently static to be
    resolved", so it is left exactly where the experiment put it. It is not a claim that those
    residues form a folded domain. What is missing is what is dynamic, and that is what gets
    built.

    This class exists because the opposite happened. Region identification needs
    MIN_FOLDED_DOMAIN_LENGTH residues before it calls anything folded, so a chain with 19
    residues resolved got NO folded domain, the whole chain became one anchor-free region, and
    DODO regenerated it and landed it on a centroid computed from placeholder coordinates.
    Measured on the 7R5J nuclear pore: 48 of the 56 Nup98 copies model only 19 residues each,
    and every one was flung up to 855 A out of the pore -- 38 detached islands in the output,
    visible as free-floating IDRs in ChimeraX.
    """

    @staticmethod
    def _sparse_chain():
        """dnmt3a with 19 residues resolved, and the full-length sequence to restore."""
        from dodo.io import read_structure

        full = read_structure(FIXTURES / "dnmt3a.pdb")
        keep = np.zeros(full.n_residues, dtype=bool)
        keep[596:615] = True
        return full.select_residues(keep), full.sequence

    def test_the_observed_residues_do_not_move_at_all(self) -> None:
        partial, reference = self._sparse_chain()
        before = partial.ca_xyz.copy()
        report = rebuild(
            partial, sequences={"A": reference}, seed=0, n_models=1, progress=False
        )
        model = report.models[0]
        observed = np.flatnonzero(~model.inserted)
        assert observed.size == 19
        assert np.abs(model.ca_xyz[observed] - before).max() == 0.0

    def test_too_few_to_be_a_domain_still_anchors_the_chain(self) -> None:
        """19 residues is below MIN_FOLDED_DOMAIN_LENGTH, and that must not matter here."""
        from dodo.constants import MIN_FOLDED_DOMAIN_LENGTH
        from dodo.structure import DomainKind

        partial, reference = self._sparse_chain()
        assert MIN_FOLDED_DOMAIN_LENGTH > 19
        report = rebuild(
            partial, sequences={"A": reference}, seed=0, n_models=1, progress=False
        )
        domains = report.models[0].chains[0].domains
        folded = [d for d in domains if d.kind is DomainKind.FOLDED]
        assert len(folded) == 1 and len(folded[0].span) == 19
        # ...and the missing residues hang off it as anchored regions, not as a free chain.
        assert [d.kind for d in domains].count(DomainKind.IDR) == 2

    def test_nothing_floats(self) -> None:
        """The symptom itself: one spatially connected piece, not two."""
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        from scipy.spatial import cKDTree

        partial, reference = self._sparse_chain()
        report = rebuild(
            partial, sequences={"A": reference}, seed=0, n_models=1, progress=False
        )
        ca = report.models[0].ca_xyz
        pairs = cKDTree(ca).query_pairs(6.0, output_type="ndarray")
        rows = np.concatenate([np.arange(len(ca) - 1), pairs[:, 0]])
        cols = np.concatenate([np.arange(1, len(ca)), pairs[:, 1]])
        count, _ = connected_components(
            coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(ca), len(ca))),
            directed=False,
        )
        assert count == 1

    def test_the_chain_comes_out_continuous(self) -> None:
        partial, reference = self._sparse_chain()
        report = rebuild(
            partial,
            sequences={"A": reference},
            seed=0,
            n_models=1,
            backbone=False,
            progress=False,
        )
        bonds = np.linalg.norm(np.diff(report.models[0].ca_xyz, axis=0), axis=1)
        assert bonds.min() > 3.5 and bonds.max() < 4.1

    def test_predicted_mode_is_untouched_by_the_rule(self) -> None:
        """The rule is about measurements. A prediction models everything and re-samples freely."""
        report = rebuild(FIXTURES / "dnmt3a.pdb", seed=0, n_models=1, progress=False)
        assert not [n for n in report.notes if "experimental input" in n]
        assert report.n_built > 0

    def test_an_experimental_file_with_nothing_missing_says_so(self) -> None:
        """Following the rule to its end: if everything was observed, there is nothing to build."""
        report = rebuild(
            FIXTURES / "6kn7.pdb", seed=0, n_models=1, backbone=False, progress=False
        )
        assert report.n_built == 0
        assert any("nothing was rebuilt" in note for note in report.notes)
        assert any("--fasta" in note for note in report.notes)
