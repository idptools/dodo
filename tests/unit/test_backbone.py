"""Tests for all-atom backbone reconstruction.

Several of these are named for a specific defect in v1's all-atom module, because every
one of those defects would have been caught by a single assertion on a bond length and
none of them were: v1 added unit vectors without scaling by bond length, so CA-C came out
at 1.000 A against a real 1.53 and the peptide bond at 2.87 A against 1.33, in a file
that declared the correct values and never imported them.

The load-bearing test is :meth:`TestClosure.test_peptide_unit_reproduces_ca_ca_bond`:
building CA-C-N-CA with trans omega from the declared constants must reproduce the
declared CA-CA virtual bond, or nothing else in the module can close onto the trace it
was given.
"""

from __future__ import annotations

from pathlib import Path
from typing import ClassVar

import numpy as np
import pytest
from scipy.spatial import cKDTree

from dodo import constants as C
from dodo.construct import backbone as B
from dodo.construct.backbone import (
    BackboneResult,
    ca_angle_budget,
    max_reconstructable_ca_angle,
    place_backbone,
    place_backbone_for_domain,
    place_sidechain_cb,
    reconstructable_ca_bond_range,
    unreconstructable_ca_angles,
)
from dodo.exceptions import GeometryError
from dodo.io.pdb import read_pdb
from dodo.structure import Domain, DomainKind, Span, Structure

DATA = Path(__file__).resolve().parents[1] / "data" / "structures"
DNMT3A = DATA / "dnmt3a.pdb"

#: Bond lengths are exact by construction here, not approximate. This tolerance is float
#: noise from a handful of trig operations, nothing physical.
EXACT = 1e-8

#: Largest CA-CA-CA pseudo-angle :func:`place_backbone` can reconstruct at the ``N-CA-C``
#: tolerance it builds with by default, in degrees.
#:
#: DERIVED, and deliberately not a literal. It is
#: ``N_CA_C_ANGLE + ca_angle_budget() + _N_CA_C_TOLERANCE``, so it moves whenever any of
#: those three moves -- and it has moved: re-measuring tau over four whole deposited
#: structures widened the tolerance from 8.0 to 14.0 degrees, which lifted this ceiling from
#: 154.5 to 160.5. Every test below that has to construct a trace the module must *refuse*
#: derives its angle from here, because a literal picked against the old ceiling quietly
#: became reconstructable and the test went on passing while exercising nothing.
RECONSTRUCTABLE_CEILING = max_reconstructable_ca_angle(n_ca_c_tolerance=B._N_CA_C_TOLERANCE)

#: A CA pseudo-angle that no trans backbone can realize, in degrees. DERIVED from
#: :data:`RECONSTRUCTABLE_CEILING` with two degrees of margin: orders of magnitude more than
#: the float noise in a pseudo-angle measured back through an arccos, and well inside the
#: headroom the MEASURED observed range leaves above the ceiling.
UNRECONSTRUCTABLE_ANGLE = RECONSTRUCTABLE_CEILING + 2.0

#: Upper CA-CA-CA bound, in degrees, for a :func:`random_walk` that must *contain* residues
#: no backbone can be built on.
#:
#: DERIVED: :data:`~dodo.constants.BACKBONE_ANGLE_OBSERVED_MAX`, the measured top of the
#: observed pseudo-angle range, which sits above :data:`RECONSTRUCTABLE_CEILING` -- as
#: constants.py explains, an angle that wide implies a tau no real residue exhibits, so the
#: observed range is not all trans backbone geometry. Sampling up to
#: :data:`~dodo.constants.BACKBONE_ANGLE_MAX` cannot produce an unreconstructable residue at
#: all, because that window is now capped *below* the ceiling on purpose, so the tests that
#: need one cannot use it.
STRAINING_ANGLE_MAX = C.BACKBONE_ANGLE_OBSERVED_MAX


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Angle a-b-c in degrees, broadcasting over leading axes."""
    u, v = a - b, c - b
    cosine = (u * v).sum(-1) / (np.linalg.norm(u, axis=-1) * np.linalg.norm(v, axis=-1))
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> np.ndarray:
    """Dihedral p0-p1-p2-p3 in degrees."""
    return B._dihedral(p0, p1, p2, p3)


def zigzag(n: int, *, pseudo_angle: float = 120.0) -> np.ndarray:
    """Make a planar CA zigzag with exact CA_CA_BOND_LENGTH and a fixed pseudo-angle.

    Planar on purpose: it is the simplest trace with a well-defined pseudo-angle, so a
    failure here is a failure of the placement geometry and not of the search.
    """
    half = np.radians(pseudo_angle) / 2.0
    xyz = np.zeros((n, 3))
    xyz[:, 0] = np.arange(n) * C.CA_CA_BOND_LENGTH * np.sin(half)
    xyz[1::2, 1] = C.CA_CA_BOND_LENGTH * np.cos(half)
    return xyz


#: Self-avoidance radius, in Angstroms, used by :func:`random_walk` unless told otherwise.
#:
#: MEASURED, and it is not :data:`~dodo.constants.CA_CLASH_DISTANCE`. Over all 912 residues
#: of dnmt3a the closest non-bonded CA-CA pair is 3.78 A -- essentially the bonded distance
#: -- so a trace that lets two alpha carbons approach the 3.20 A hard-sphere limit is
#: already outside what real protein backbones do, and the atoms hanging off those two
#: alpha carbons then have nowhere to go. A *free* random walk is worse still: at 60
#: residues it puts non-bonded alpha carbons at 0.2 A of each other, which is not a trace
#: any all-atom reconstruction can rescue. See
#: :meth:`TestReconstructableWindow.test_acceptance_depends_on_how_tightly_the_trace_packs`
#: for the measured acceptance curve against this radius.
SELF_AVOID = 4.5


def random_walk(
    n: int,
    rng: np.random.Generator,
    *,
    angle_max: float,
    self_avoid: float = SELF_AVOID,
) -> np.ndarray:
    """Make a self-avoiding CA trace with exact CA_CA_BOND_LENGTH bonds and bounded angles.

    This is the shape of trace DODO's own IDR engine produces, which is why it is worth
    testing against rather than only using real structures. Self-avoiding by rejection with
    backtracking, and deterministic given ``rng``.
    """
    xyz = np.zeros((n, 3))
    xyz[1, 0] = C.CA_CA_BOND_LENGTH
    index = 2
    restarts = 0
    while index < n:
        placed = False
        for _ in range(400):
            theta = np.radians(rng.uniform(C.BACKBONE_ANGLE_MIN, angle_max))
            previous = xyz[index - 1] - xyz[index - 2]
            previous /= np.linalg.norm(previous)
            perpendicular = rng.normal(size=3)
            perpendicular -= perpendicular @ previous * previous
            perpendicular /= np.linalg.norm(perpendicular)
            step = np.cos(np.pi - theta) * previous + np.sin(np.pi - theta) * perpendicular
            candidate = xyz[index - 1] + C.CA_CA_BOND_LENGTH * step
            if index < 2 or bool(
                np.all(np.linalg.norm(xyz[: index - 1] - candidate, axis=1) >= self_avoid)
            ):
                xyz[index] = candidate
                placed = True
                break
        if not placed:
            index = max(2, index - 3)
            restarts += 1
            if restarts > 400:
                raise AssertionError(
                    f"could not build a {n}-residue walk self-avoiding at {self_avoid} A"
                )
            continue
        index += 1
    return xyz


def split_atoms(result: BackboneResult) -> dict[str, np.ndarray]:
    """Group a result's coordinates by atom name."""
    return {name: result.xyz[result.atom_names == name] for name in np.unique(result.atom_names)}


@pytest.fixture(scope="module")
def dnmt3a() -> Structure:
    return read_pdb(DNMT3A)


def reference_atom(structure: Structure, name: str, start: int, stop: int) -> np.ndarray:
    """Deposited coordinates of one atom name over a residue range, NaN where absent."""
    out = np.full((stop - start, 3), np.nan)
    for residue in range(start, stop):
        atoms = structure.atom_slice_for_residues(residue, residue + 1)
        found = np.flatnonzero(structure.atom_name[atoms] == name)
        if found.size:
            out[residue - start] = structure.xyz[atoms][found[0]]
    return out


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Per-atom RMSD between two coordinate sets, ignoring NaN rows."""
    return float(np.sqrt(np.nanmean(((a - b) ** 2).sum(-1))))


#: Ramachandran 1963 "normally allowed" hard-sphere contact limits, in Angstroms, keyed
#: on the element pair. Transcribed here rather than imported so that the contact
#: measurement below is independent of the module under test: a bug in the module's own
#: table must not be able to hide a clash from these tests.
CONTACT_LIMIT: dict[frozenset[str], float] = {
    frozenset(("C", "C")): 3.20,
    frozenset(("C", "N")): 2.90,
    frozenset(("C", "O")): 2.80,
    frozenset(("N", "N")): 2.70,
    frozenset(("N", "O")): 2.70,
    frozenset(("O", "O")): 2.80,
}


def nonbonded_contacts(
    xyz: np.ndarray,
    residue: np.ndarray,
    names: np.ndarray,
    *,
    breaks: tuple[int, ...] = (),
) -> list[tuple[float, float, str]]:
    """Every non-bonded heavy-atom pair, closest first, with its contact limit.

    Deliberately written from scratch instead of calling the module's own steric code:
    the point of these tests is that atoms at 1.28 A cannot be reported as a valid
    backbone, and reusing the checker under test would make that unfalsifiable.

    Only pairs at least two residues apart are counted, plus every pair spanning a chain
    break, where no bond joins the two residues at all. Closer than that the 1963 limits
    do not apply and real structures say so: over all of dnmt3a, 0 of 507 pairs two or
    more residues apart are below their limit (closest margin +0.05 A), while 463 of 1490
    consecutive-residue pairs are -- ``C(i)..C(i+1)`` reaches 2.92 A against a 3.20 A
    C-C limit throughout the helices, because that pair is four bonds apart across a
    planar peptide unit and has nowhere else to go.
    """
    tree = cKDTree(xyz)
    element = {"N": "N", "O": "O", "C": "C", "CA": "C", "CB": "C"}
    out: list[tuple[float, float, str]] = []
    for i, j in tree.query_pairs(max(CONTACT_LIMIT.values())):
        low, high = (i, j) if residue[i] <= residue[j] else (j, i)
        separation = int(residue[high] - residue[low])
        broken = int(residue[low]) in breaks and separation == 1
        if separation < 2 and not broken:
            continue
        limit = CONTACT_LIMIT[frozenset((element[str(names[low])], element[str(names[high])]))]
        distance = float(np.linalg.norm(xyz[low] - xyz[high]))
        if distance < limit:
            out.append(
                (
                    distance,
                    limit,
                    f"{names[low]}{residue[low]}..{names[high]}{residue[high]}",
                )
            )
    return sorted(out)


def closest_nonbonded(
    xyz: np.ndarray, residue: np.ndarray, names: np.ndarray, *, breaks: tuple[int, ...] = ()
) -> float:
    """Smallest separation among the pairs :func:`nonbonded_contacts` judges."""
    closest = float("inf")
    tree = cKDTree(xyz)
    for i, j in tree.query_pairs(2 * max(CONTACT_LIMIT.values())):
        low, high = (i, j) if residue[i] <= residue[j] else (j, i)
        separation = int(residue[high] - residue[low])
        broken = int(residue[low]) in breaks and separation == 1
        if separation < 2 and not broken:
            continue
        if broken and str(names[low]) == "CA" and str(names[high]) == "CA":
            continue
        closest = min(closest, float(np.linalg.norm(xyz[low] - xyz[high])))
    return closest


def result_contacts(result: BackboneResult, offset: int = 0) -> list[tuple[float, float, str]]:
    """:func:`nonbonded_contacts` on a :class:`BackboneResult`."""
    return nonbonded_contacts(
        result.xyz,
        result.residue_index + offset,
        result.atom_names,
        breaks=result.chain_breaks,
    )


# ---------------------------------------------------------------------------
# The consistency check the whole module rests on
# ---------------------------------------------------------------------------


class TestClosure:
    def test_peptide_unit_reproduces_ca_ca_bond(self) -> None:
        # Build CA-C-N-CA with trans omega from the declared constants and nothing else.
        # If this does not reproduce CA_CA_BOND_LENGTH, a reconstruction cannot close
        # onto the trace it was handed and every other test here is meaningless.
        distance = B._unit_ca_distance(C.CA_C_N_ANGLE, C.C_N_CA_ANGLE)
        assert distance == pytest.approx(3.803955, abs=1e-6)
        assert distance == pytest.approx(C.CA_CA_BOND_LENGTH, abs=0.01)

    def test_reconstructable_window_brackets_the_constant(self) -> None:
        low, high = reconstructable_ca_bond_range()
        assert low < C.CA_CA_BOND_LENGTH < high
        # Wide enough for a real AlphaFold model's trans bonds (measured 3.769-3.913 on
        # dnmt3a's folded domains) and narrow enough to exclude cis-proline at ~2.9 A.
        assert low < 3.769
        assert high > 3.913
        assert low > 3.2

    def test_ca_angle_budget_is_invariant_to_strain(self) -> None:
        # alpha_C + alpha_N is what caps the reconstructable CA pseudo-angle, and the
        # cap is only meaningful if straining the peptide bond angles cannot lift it.
        budget = ca_angle_budget()
        assert budget == pytest.approx(35.50, abs=0.05)
        # The claim is about the family of bond-angle pairs that hold the CA-CA distance
        # *fixed*, since that is the family a reconstruction is free to move along. These
        # pairs are the table in ca_angle_budget's docstring.
        iso_distance = np.array(
            [
                [116.20, 121.70],
                [120.20, 117.22],
                [124.20, 114.28],
                [128.20, 112.28],
                [132.20, 110.94],
            ]
        )
        built = B._peptide_unit(iso_distance[:, 0], iso_distance[:, 1])
        distance = np.linalg.norm(built["CA2"], axis=1)
        assert np.allclose(distance, 3.803955, atol=5e-4)
        x_axis = built["CA2"] / distance[:, None]
        alpha_c = angle(built["C"], built["CA1"], built["CA2"])
        alpha_n = angle(built["N"], built["CA2"], built["CA1"])
        # alpha_C moves by more than 13 degrees across this family; the sum moves by less
        # than one. That is what makes the pseudo-angle ceiling a real ceiling.
        assert np.ptp(alpha_c) > 13.0
        assert np.ptp(alpha_c + alpha_n) < 1.0
        assert np.allclose(alpha_c + alpha_n, budget, atol=0.6)
        assert x_axis.shape == (5, 3)

    def test_max_reconstructable_ca_angle_is_below_the_sampling_window(self) -> None:
        # Not an aspiration: BACKBONE_ANGLE_MAX is genuinely above the geometric
        # ceiling, which is why place_backbone has to be able to refuse a trace.
        assert max_reconstructable_ca_angle() == pytest.approx(146.5, abs=0.1)
        assert max_reconstructable_ca_angle() < C.BACKBONE_ANGLE_MAX


# ---------------------------------------------------------------------------
# Bond lengths and angles -- the v1 regressions
# ---------------------------------------------------------------------------


class TestBondGeometry:
    @pytest.fixture(scope="class")
    def built(self) -> tuple[np.ndarray, BackboneResult]:
        ca = zigzag(24)
        return ca, place_backbone(ca, "A" * 24, np.random.default_rng(7))

    def test_ca_c_bond_is_not_one_angstrom(self, built: tuple[np.ndarray, BackboneResult]) -> None:
        # v1 added a unit vector without scaling it, so CA-C came out at exactly 1.000 A.
        ca, result = built
        atoms = split_atoms(result)
        lengths = np.linalg.norm(atoms["C"] - ca, axis=1)
        assert np.allclose(lengths, C.CA_C_BOND_LENGTH, atol=EXACT)
        assert not np.any(np.isclose(lengths, 1.0, atol=0.01))

    def test_peptide_bond_is_not_two_point_eight_seven(
        self, built: tuple[np.ndarray, BackboneResult]
    ) -> None:
        # v1's peptide bond was 2.87 A against a declared 1.329.
        _, result = built
        atoms = split_atoms(result)
        lengths = np.linalg.norm(atoms["N"][1:] - atoms["C"][:-1], axis=1)
        assert np.allclose(lengths, C.C_N_PEPTIDE_BOND_LENGTH, atol=EXACT)

    def test_n_ca_and_c_o_bonds_are_exact(self, built: tuple[np.ndarray, BackboneResult]) -> None:
        ca, result = built
        atoms = split_atoms(result)
        assert np.allclose(np.linalg.norm(atoms["N"] - ca, axis=1), C.N_CA_BOND_LENGTH, atol=EXACT)
        assert np.allclose(
            np.linalg.norm(atoms["O"] - atoms["C"], axis=1), C.C_O_BOND_LENGTH, atol=EXACT
        )

    def test_omega_is_trans(self, built: tuple[np.ndarray, BackboneResult]) -> None:
        ca, result = built
        atoms = split_atoms(result)
        omega = dihedral(ca[:-1], atoms["C"][:-1], atoms["N"][1:], ca[1:])
        assert np.allclose(np.abs(omega), C.OMEGA_TRANS, atol=1e-6)

    def test_peptide_bond_angles_match_within_declared_strain(
        self, built: tuple[np.ndarray, BackboneResult]
    ) -> None:
        ca, result = built
        atoms = split_atoms(result)
        ca_c_n = angle(ca[:-1], atoms["C"][:-1], atoms["N"][1:])
        c_n_ca = angle(atoms["C"][:-1], atoms["N"][1:], ca[1:])
        # This trace's bonds are exactly CA_CA_BOND_LENGTH, within 0.01 A of the unstrained
        # peptide unit's own 3.804 A, so the strain needed to close onto them is tiny.
        assert np.allclose(ca_c_n, C.CA_C_N_ANGLE, atol=0.5)
        assert np.allclose(c_n_ca, C.C_N_CA_ANGLE, atol=0.5)
        assert result.max_peptide_angle_strain < 1.0

    def test_carbonyl_is_planar_and_ca_c_o_angle_is_exact(
        self, built: tuple[np.ndarray, BackboneResult]
    ) -> None:
        # CA_C_O_ANGLE is never used to place an interior O; it falls out of placing the
        # O in the sp2 plane. Checking it is therefore a real check on the derivation of
        # _N_C_O_ANGLE, not a tautology.
        ca, result = built
        atoms = split_atoms(result)
        assert np.allclose(
            angle(ca[:-1], atoms["C"][:-1], atoms["O"][:-1]), C.CA_C_O_ANGLE, atol=1e-6
        )
        # The derived N-C-O angle closes the sp2 carbon to 360 degrees.
        assert abs(B._N_C_O_ANGLE - 123.0) < 1e-9
        assert np.allclose(
            angle(atoms["N"][1:], atoms["C"][:-1], atoms["O"][:-1])
            + angle(ca[:-1], atoms["C"][:-1], atoms["O"][:-1])
            + angle(ca[:-1], atoms["C"][:-1], atoms["N"][1:]),
            360.0,
            atol=1e-6,
        )
        # The five atoms of the peptide unit are coplanar.
        for i in range(len(ca) - 1):
            unit = np.array([ca[i], atoms["C"][i], atoms["O"][i], atoms["N"][i + 1], ca[i + 1]])
            centred = unit - unit.mean(axis=0)
            singular = np.linalg.svd(centred, compute_uv=False)
            assert singular[2] < 1e-6

    def test_n_ca_c_angle_is_within_tolerance(
        self, built: tuple[np.ndarray, BackboneResult]
    ) -> None:
        _, result = built
        assert result.max_n_ca_c_deviation < B._N_CA_C_TOLERANCE
        # A 120-degree zigzag is well inside the reconstructable window, so tau should be
        # essentially ideal rather than merely inside tolerance.
        assert result.max_n_ca_c_deviation < 1.0


# ---------------------------------------------------------------------------
# The alpha carbons are input, not output
# ---------------------------------------------------------------------------


class TestCaCoordinatesUnchanged:
    def test_ca_is_bit_identical(self) -> None:
        ca = zigzag(30, pseudo_angle=132.0)
        result = place_backbone(ca, "G" * 30, np.random.default_rng(0), sidechains=True)
        assert np.array_equal(result.xyz[result.atom_names == "CA"], ca)

    def test_ca_is_bit_identical_on_a_real_trace(self, dnmt3a: Structure) -> None:
        ca = dnmt3a.ca_xyz[700:780]
        result = place_backbone(ca, dnmt3a.sequence[700:780], np.random.default_rng(0))
        assert np.array_equal(result.xyz[result.atom_names == "CA"], ca)

    def test_input_array_is_not_mutated(self) -> None:
        ca = zigzag(12)
        before = ca.copy()
        place_backbone(ca, "A" * 12, np.random.default_rng(0))
        assert np.array_equal(ca, before)


# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------


class TestSeeding:
    def test_same_seed_is_bit_identical(self) -> None:
        ca = zigzag(20, pseudo_angle=128.0)
        first = place_backbone(ca, "A" * 20, np.random.default_rng(11))
        second = place_backbone(ca, "A" * 20, np.random.default_rng(11))
        assert np.array_equal(first.xyz, second.xyz)

    def test_different_seeds_differ_but_stay_valid(self) -> None:
        # A CA trace does not determine its backbone uniquely, and the rng is how that
        # residual freedom is exposed. Both answers must still be physically valid.
        ca = zigzag(20, pseudo_angle=128.0)
        first = place_backbone(ca, "A" * 20, np.random.default_rng(1))
        second = place_backbone(ca, "A" * 20, np.random.default_rng(9999))
        assert not np.array_equal(first.xyz, second.xyz)
        for result in (first, second):
            atoms = split_atoms(result)
            assert np.allclose(
                np.linalg.norm(atoms["N"][1:] - atoms["C"][:-1], axis=1),
                C.C_N_PEPTIDE_BOND_LENGTH,
                atol=EXACT,
            )
            assert result.max_n_ca_c_deviation < B._N_CA_C_TOLERANCE

    def test_no_global_numpy_state_is_used(self) -> None:
        # Rule: nothing in DODO may consume numpy's global random state, or two runs of
        # the "same" seeded command disagree.
        ca = zigzag(16)
        np.random.seed(0)  # noqa: NPY002
        first = place_backbone(ca, "A" * 16, np.random.default_rng(5))
        np.random.seed(12345)  # noqa: NPY002
        second = place_backbone(ca, "A" * 16, np.random.default_rng(5))
        assert np.array_equal(first.xyz, second.xyz)


# ---------------------------------------------------------------------------
# Terminal residues
# ---------------------------------------------------------------------------


class TestTerminals:
    def test_both_termini_get_a_full_complement(self) -> None:
        # The failure mode being prevented is dropping atoms at the ends, which produces
        # a file whose first residue has no N and whose last has no C or O.
        ca = zigzag(9)
        result = place_backbone(ca, "A" * 9, np.random.default_rng(0))
        atoms = split_atoms(result)
        for name in ("N", "CA", "C", "O"):
            assert atoms[name].shape == (9, 3)
        assert np.isfinite(result.xyz).all()

    def test_terminal_dihedrals_are_the_declared_values(self) -> None:
        ca = zigzag(9)
        result = place_backbone(ca, "A" * 9, np.random.default_rng(0))
        atoms = split_atoms(result)
        psi_first = dihedral(atoms["N"][0], ca[0], atoms["C"][0], atoms["N"][1])
        assert float(psi_first) == pytest.approx(B._TERMINAL_PSI, abs=1e-6)
        phi_last = dihedral(atoms["C"][-2], atoms["N"][-1], ca[-1], atoms["C"][-1])
        assert float(phi_last) == pytest.approx(B._TERMINAL_PHI, abs=1e-6)

    def test_terminal_n_ca_c_is_ideal(self) -> None:
        ca = zigzag(9)
        result = place_backbone(ca, "A" * 9, np.random.default_rng(0))
        # Both termini have a free dihedral, so their tau is placed exactly rather than
        # traded off against anything.
        assert result.n_ca_c_angles[0] == pytest.approx(C.N_CA_C_ANGLE, abs=1e-6)
        assert result.n_ca_c_angles[-1] == pytest.approx(C.N_CA_C_ANGLE, abs=1e-6)

    def test_two_residues_is_the_minimum(self) -> None:
        ca = zigzag(2)
        result = place_backbone(ca, "AA", np.random.default_rng(0))
        assert result.n_residues == 2
        atoms = split_atoms(result)
        assert np.linalg.norm(atoms["N"][1] - atoms["C"][0]) == pytest.approx(
            C.C_N_PEPTIDE_BOND_LENGTH, abs=EXACT
        )

    def test_no_oxt_and_no_hydrogens(self) -> None:
        result = place_backbone(zigzag(6), "AAAAAA", np.random.default_rng(0), sidechains=True)
        assert set(np.unique(result.atom_names)) == {"N", "CA", "C", "O", "CB"}
        assert set(np.unique(result.elements)) == {"N", "C", "O"}


# ---------------------------------------------------------------------------
# CB
# ---------------------------------------------------------------------------


class TestSidechainCb:
    def test_cb_absent_by_default(self) -> None:
        result = place_backbone(zigzag(6), "AAAAAA", np.random.default_rng(0))
        assert "CB" not in set(np.unique(result.atom_names))
        assert result.sidechains is False

    def test_glycine_gets_no_cb(self) -> None:
        result = place_backbone(zigzag(5), "AGAGA", np.random.default_rng(0), sidechains=True)
        assert (result.atom_names == "CB").sum() == 3
        cb_residues = set(result.residue_index[result.atom_names == "CB"].tolist())
        assert cb_residues == {0, 2, 4}

    def test_unknown_residue_gets_no_cb(self) -> None:
        # Inventing a side chain for a residue whose identity is unknown is inventing
        # information, which is the specific sin this module exists to avoid.
        result = place_backbone(zigzag(4), "AXAX", np.random.default_rng(0), sidechains=True)
        assert set(result.residue_index[result.atom_names == "CB"].tolist()) == {0, 2}

    def test_cb_geometry_matches_the_declared_angles(self) -> None:
        ca = zigzag(10)
        result = place_backbone(ca, "A" * 10, np.random.default_rng(0), sidechains=True)
        atoms = split_atoms(result)
        assert np.allclose(
            np.linalg.norm(atoms["CB"] - ca, axis=1), B._CA_CB_BOND_LENGTH, atol=1e-9
        )
        assert np.allclose(angle(atoms["N"], ca, atoms["CB"]), B._N_CA_CB_ANGLE, atol=1e-6)
        assert np.allclose(angle(atoms["C"], ca, atoms["CB"]), B._C_CA_CB_ANGLE, atol=1e-6)

    def test_cb_chirality_matches_a_real_l_amino_acid_structure(self, dnmt3a: Structure) -> None:
        # The chirality sign is the one thing about CB that cannot be derived from bond
        # angles, so it is measured against deposited coordinates rather than asserted.
        start, stop = 700, 800
        ca = dnmt3a.ca_xyz[start:stop]
        reference_n = reference_atom(dnmt3a, "N", start, stop)
        reference_c = reference_atom(dnmt3a, "C", start, stop)
        reference_cb = reference_atom(dnmt3a, "CB", start, stop)
        placed = place_sidechain_cb(reference_n, ca, reference_c)
        has_cb = np.isfinite(reference_cb).all(axis=1)
        assert rmsd(placed[has_cb], reference_cb[has_cb]) < 0.15
        # The mirror image would be about 2.4 A away, so the sign is unambiguous.
        mirrored = 2 * ca - placed
        assert rmsd(mirrored[has_cb], reference_cb[has_cb]) > 1.5

    def test_cb_rejects_a_collinear_backbone(self) -> None:
        straight = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]])
        with pytest.raises(GeometryError, match="collinear"):
            place_sidechain_cb(straight[0], straight[1], straight[2])


# ---------------------------------------------------------------------------
# Failures are loud
# ---------------------------------------------------------------------------


class TestFailuresAreLoud:
    def test_single_residue_raises(self) -> None:
        with pytest.raises(GeometryError, match="peptide unit needs"):
            place_backbone(np.zeros((1, 3)), "A", np.random.default_rng(0))

    def test_sequence_length_mismatch_raises(self) -> None:
        with pytest.raises(GeometryError, match="sequence has 3 residues"):
            place_backbone(zigzag(5), "AAA", np.random.default_rng(0))

    def test_wrong_shape_raises(self) -> None:
        with pytest.raises(GeometryError, match=r"\(n_residues, 3\)"):
            place_backbone(np.zeros((5, 2)), "AAAAA", np.random.default_rng(0))

    def test_nan_input_raises_rather_than_propagating(self) -> None:
        # v1's samplers returned NaN on failure and the writer wrote it. A NaN reaching
        # this module means an upstream failure, and laundering it into a structure file
        # is the worst available outcome.
        ca = zigzag(6)
        ca[3, 1] = np.nan
        with pytest.raises(GeometryError, match="non-finite"):
            place_backbone(ca, "A" * 6, np.random.default_rng(0))

    def test_cis_peptide_bond_raises_and_says_so(self) -> None:
        ca = zigzag(8)
        # Pull two alpha carbons to a cis-proline separation.
        direction = ca[4] - ca[3]
        ca[4:] -= direction * (1.0 - 2.95 / np.linalg.norm(direction))
        with pytest.raises(GeometryError, match="cis peptide bond"):
            place_backbone(ca, "A" * 8, np.random.default_rng(0))

    def test_duplicate_alpha_carbons_raise(self) -> None:
        ca = zigzag(6)
        ca[3] = ca[2]
        with pytest.raises(GeometryError, match="same position"):
            place_backbone(ca, "A" * 6, np.random.default_rng(0))

    def test_unreconstructable_pseudo_angle_raises_and_names_the_ceiling(self) -> None:
        # A trace whose CA pseudo-angle exceeds the ceiling for the tolerance being built
        # with cannot have an N-CA-C angle inside that tolerance at all. Emitting one quietly
        # is what the module must not do; the message has to point upstream. The angle comes
        # from RECONSTRUCTABLE_CEILING rather than from BACKBONE_ANGLE_MAX, which is capped
        # below the ceiling on purpose and therefore cannot express this failure.
        ca = zigzag(12, pseudo_angle=UNRECONSTRUCTABLE_ANGLE)
        with pytest.raises(GeometryError, match="max_reconstructable_ca_angle"):
            place_backbone(ca, "A" * 12, np.random.default_rng(0))

    def test_infinite_tolerance_is_refused_rather_than_honoured(self) -> None:
        # The old error message recommended n_ca_c_tolerance=float("inf") as a way to
        # inspect a failing trace. Measured, it did not inspect anything: it returned a
        # result indistinguishable from a success while carrying N-CA-C angles from 101.6
        # to 132.9 deg and heavy-atom pairs at 1.19 A. A loud failure turned into a silent
        # one, which is the opposite of the intent.
        ca = zigzag(12, pseudo_angle=UNRECONSTRUCTABLE_ANGLE)
        with pytest.raises(GeometryError, match="must be finite"):
            place_backbone(ca, "A" * 12, np.random.default_rng(0), n_ca_c_tolerance=float("inf"))
        with pytest.raises(GeometryError, match="must be positive"):
            place_backbone(ca, "A" * 12, np.random.default_rng(0), n_ca_c_tolerance=0.0)

    def test_a_failing_trace_is_inspectable_without_building_it(self) -> None:
        # What replaces the escape hatch: the same question, answered before anything is
        # placed, and without producing a structure file nobody can tell is wrong.
        ca = zigzag(12, pseudo_angle=UNRECONSTRUCTABLE_ANGLE)
        impossible = unreconstructable_ca_angles(ca)
        assert impossible.size == 10  # every interior residue of this zigzag
        assert list(impossible) == list(range(1, 11))


# ---------------------------------------------------------------------------
# Chain breaks
# ---------------------------------------------------------------------------


class TestChainBreaks:
    def _broken(self) -> tuple[np.ndarray, str]:
        ca = zigzag(14)
        ca[7:] += np.array([12.0, 0.0, 0.0])
        return ca, "A" * 14

    def test_break_raises_by_default(self) -> None:
        ca, sequence = self._broken()
        with pytest.raises(GeometryError, match="allow_chain_breaks"):
            place_backbone(ca, sequence, np.random.default_rng(0))

    def test_break_is_reported_when_allowed(self) -> None:
        ca, sequence = self._broken()
        result = place_backbone(ca, sequence, np.random.default_rng(0), allow_chain_breaks=True)
        assert result.chain_breaks == (6,)
        atoms = split_atoms(result)
        # Every peptide bond except the one across the break is exact.
        lengths = np.linalg.norm(atoms["N"][1:] - atoms["C"][:-1], axis=1)
        intact = np.ones(lengths.shape[0], dtype=bool)
        intact[6] = False
        assert np.allclose(lengths[intact], C.C_N_PEPTIDE_BOND_LENGTH, atol=EXACT)
        assert lengths[6] > 5.0
        # And the break position carries no strain value, rather than a misleading zero.
        assert np.isnan(result.peptide_angle_strain[6])
        assert np.isfinite(result.max_peptide_angle_strain)

    def test_isolated_residue_between_two_breaks_raises(self) -> None:
        ca = zigzag(9)
        ca[4] += np.array([0.0, 0.0, 14.0])
        ca[5:] += np.array([0.0, 0.0, 28.0])
        with pytest.raises(GeometryError, match="isolated residue"):
            place_backbone(ca, "A" * 9, np.random.default_rng(0), allow_chain_breaks=True)


# ---------------------------------------------------------------------------
# The arrays are usable
# ---------------------------------------------------------------------------


class TestAtomRecords:
    def test_records_build_a_valid_structure(self) -> None:
        ca = zigzag(15, pseudo_angle=126.0)
        result = place_backbone(ca, "ACDEFGHIKLMNPQR", np.random.default_rng(3), sidechains=True)
        structure = Structure.from_atom_records(
            **result.atom_records(chain_id="B", first_residue_number=42),
            source="test",
        )
        structure.validate()
        assert structure.n_residues == 15
        assert structure.sequence == "ACDEFGHIKLMNPQR"
        assert structure.residue_number[0] == 42
        assert structure.chains[0].chain_id == "B"
        # And the CA trace survives the round trip untouched.
        assert np.allclose(structure.ca_xyz, ca, atol=0.0)

    def test_atoms_are_residue_ordered(self) -> None:
        result = place_backbone(zigzag(8), "AAGAAGAA", np.random.default_rng(0), sidechains=True)
        assert np.all(np.diff(result.residue_index) >= 0)
        # N, CA, C, O then CB: the canonical within-residue order, which is what makes
        # the writer's CONECT logic find the backbone.
        first = result.atom_names[result.residue_index == 0]
        assert list(first) == ["N", "CA", "C", "O", "CB"]

    def test_repr_reports_the_geometry(self) -> None:
        result = place_backbone(zigzag(8), "A" * 8, np.random.default_rng(0))
        text = repr(result)
        assert "N-CA-C" in text
        assert "strain" in text


# ---------------------------------------------------------------------------
# Behaviour on the traces DODO itself generates
# ---------------------------------------------------------------------------


class TestGeneratedTraces:
    def test_inside_the_ceiling_the_angle_is_never_the_reason_for_a_refusal(self) -> None:
        # With the CA-CA-CA window kept below the geometric ceiling, N-CA-C is achievable at
        # every residue, so it must never be what a refusal is about. A refusal can still
        # happen -- a trace can pack its own carbonyls into each other -- but then the
        # message says so, and the two causes stay distinguishable. This is the property the
        # old solver broke: it refused traces on N-CA-C that had no impossible residue in
        # them at all.
        accepted = 0
        for seed in range(4):
            ca = random_walk(50, np.random.default_rng(seed), angle_max=135.0)
            try:
                result = place_backbone(ca, "A" * 50, np.random.default_rng(seed), sidechains=True)
            except GeometryError as refusal:
                assert "N-CA-C angle more than" not in str(refusal), seed
                assert "hard-sphere" in str(refusal), seed
                continue
            accepted += 1
            assert result.max_n_ca_c_deviation <= B._N_CA_C_TOLERANCE
            assert result_contacts(result) == []
        # MEASURED: 3 of 4 at this length and self-avoidance radius. Pinned as a floor so a
        # regression that refuses everything cannot pass by refusing loudly.
        assert accepted >= 3

    def test_full_sampling_window_is_not_always_reconstructable(self) -> None:
        # This is the finding, asserted so it cannot regress unnoticed: a trace sampled
        # up to BACKBONE_ANGLE_MAX is sometimes impossible to reconstruct, and the module
        # refuses rather than emitting a 20-degree N-CA-C angle.
        failures = 0
        for seed in range(6):
            ca = random_walk(50, np.random.default_rng(seed), angle_max=C.BACKBONE_ANGLE_MAX)
            try:
                place_backbone(ca, "A" * 50, np.random.default_rng(0))
            except GeometryError:
                failures += 1
        assert failures > 0


# ---------------------------------------------------------------------------
# The real-structure regression: strip to CA, rebuild, compare
# ---------------------------------------------------------------------------


class TestAgainstDeposited:
    #: Thresholds are set just above the measured values so a regression trips them.
    #: Measured per segment (dnmt3a, four folded segments, 562 residues, seed 20260730):
    #: N 0.167-0.348, C 0.224-0.381, O 0.734-1.042, CB 0.234-0.469 A RMSD.
    #:
    #: O is the worst by a factor of three for a structural reason and not a bug: it sits
    #: 1.75 A off the CA-CA axis against C's 0.54, so it amplifies whatever residual
    #: peptide-plane rotation error is left, and resolving that last rotation is exactly
    #: what PULCHRA and BBQ use a statistical library keyed on the CA pseudo-dihedral for.
    #: DODO ships no such library, so this is the honest ceiling of a pure-geometry method.
    #: An RMSD is the wrong summary of it on its own, because the distribution is
    #: heavy-tailed: see :class:`TestCarbonylTail`.
    #:
    #: Two things these RMSDs include without distinguishing, both of which are expected to
    #: be worse than the interior and neither of which is a defect. First, the *terminal*
    #: atoms -- the first residue's N and the last residue's C and O -- have a free dihedral
    #: that no CA determines, so they are placed from a stated default (or rotated off it to
    #: avoid a clash). Measured errors: N(first) 0.21 / 0.46 / 2.39 / 2.74 A, C(last) 0.32 /
    #: 0.16 / 2.61 / 0.48 A, O(last) 0.54 / 2.10 / 4.11 / 0.53 A. Second, a segment cut out
    #: of a larger protein has no neighbour to close onto at either end, which is exactly
    #: what :attr:`~dodo.structure.Span.n_anchor` exists to fix for an in-place rebuild.
    LIMITS: ClassVar[dict[str, float]] = {"N": 0.40, "C": 0.45, "O": 1.20, "CB": 0.55}

    SEGMENTS: ClassVar[tuple[tuple[int, int], ...]] = (
        (699, 912),
        (579, 698),
        (282, 414),
        (480, 578),
    )

    @pytest.mark.parametrize(("start", "stop"), SEGMENTS)
    def test_rebuilt_backbone_matches_deposited(
        self, dnmt3a: Structure, start: int, stop: int
    ) -> None:
        ca = dnmt3a.ca_xyz[start:stop]
        result = place_backbone(
            ca, dnmt3a.sequence[start:stop], np.random.default_rng(20260730), sidechains=True
        )
        atoms = split_atoms(result)
        for name in ("N", "C", "O"):
            deposited = reference_atom(dnmt3a, name, start, stop)
            assert rmsd(atoms[name], deposited) < self.LIMITS[name], name
        deposited_cb = reference_atom(dnmt3a, "CB", start, stop)
        has_cb = np.isfinite(deposited_cb).all(axis=1)
        assert rmsd(atoms["CB"], deposited_cb[has_cb]) < self.LIMITS["CB"]

    def test_rebuilt_backbone_is_physically_valid_on_a_real_trace(self, dnmt3a: Structure) -> None:
        start, stop = 699, 912
        ca = dnmt3a.ca_xyz[start:stop]
        result = place_backbone(ca, dnmt3a.sequence[start:stop], np.random.default_rng(20260730))
        atoms = split_atoms(result)
        assert np.allclose(np.linalg.norm(atoms["N"] - ca, axis=1), C.N_CA_BOND_LENGTH, atol=EXACT)
        assert np.allclose(np.linalg.norm(atoms["C"] - ca, axis=1), C.CA_C_BOND_LENGTH, atol=EXACT)
        assert np.allclose(
            np.linalg.norm(atoms["N"][1:] - atoms["C"][:-1], axis=1),
            C.C_N_PEPTIDE_BOND_LENGTH,
            atol=EXACT,
        )
        assert result.max_n_ca_c_deviation < B._N_CA_C_TOLERANCE
        assert result.max_peptide_angle_strain < B._MAX_PEPTIDE_ANGLE_STRAIN

    def test_positive_phi_stays_rare(self, dnmt3a: Structure) -> None:
        # L-amino acids disfavour phi > 0. Measured: 6% in the deposited model, and the
        # rebuild must not be dramatically worse or the chirality term is not working.
        start, stop = 699, 912
        ca = dnmt3a.ca_xyz[start:stop]
        result = place_backbone(ca, dnmt3a.sequence[start:stop], np.random.default_rng(20260730))
        atoms = split_atoms(result)
        phi = dihedral(atoms["C"][:-1], atoms["N"][1:], ca[1:], atoms["C"][1:])
        assert float((phi > 0).mean()) < 0.12


# ---------------------------------------------------------------------------
# In-place domain rebuild
# ---------------------------------------------------------------------------


class TestPlaceBackboneForDomain:
    def _structure(self, n: int = 20) -> Structure:
        ca = zigzag(n, pseudo_angle=126.0)
        result = place_backbone(ca, "A" * n, np.random.default_rng(0))
        return Structure.from_atom_records(
            **result.atom_records(),
            source="fixture",
        )

    def test_rewrites_in_place_without_touching_alpha_carbons(self) -> None:
        structure = self._structure()
        domain = Domain(
            structure=structure, span=Span(4, 16, n_anchor=3, c_anchor=16), kind=DomainKind.IDR
        )
        ca_before = structure.ca_xyz.copy()
        n_before = structure.xyz[structure.atom_name == "N"].copy()
        place_backbone_for_domain(structure, domain, np.random.default_rng(5))
        assert np.array_equal(structure.ca_xyz, ca_before)
        # Something inside the domain actually moved, or the call did nothing.
        n_after = structure.xyz[structure.atom_name == "N"]
        assert not np.array_equal(n_before[4:16], n_after[4:16])
        # Nothing before the N-side anchor moved. The C-side anchor's own N *does* move,
        # and must: it is the second half of the peptide unit that closes the span onto
        # that anchor, and leaving it in place is what left that bond at 1.2672 A.
        assert np.array_equal(n_before[:4], n_after[:4])
        assert not np.array_equal(n_before[16], n_after[16])
        assert np.array_equal(n_before[17:], n_after[17:])

    def test_result_is_still_valid_geometry(self) -> None:
        structure = self._structure()
        domain = Domain(
            structure=structure, span=Span(2, 18, n_anchor=1, c_anchor=18), kind=DomainKind.IDR
        )
        place_backbone_for_domain(structure, domain, np.random.default_rng(6))
        ca = structure.ca_xyz
        nitrogen = structure.xyz[structure.atom_name == "N"]
        carbon = structure.xyz[structure.atom_name == "C"]
        assert np.allclose(
            np.linalg.norm(nitrogen[2:18] - ca[2:18], axis=1), C.N_CA_BOND_LENGTH, atol=EXACT
        )
        assert np.allclose(
            np.linalg.norm(carbon[2:18] - ca[2:18], axis=1), C.CA_C_BOND_LENGTH, atol=EXACT
        )
        assert np.allclose(
            angle(nitrogen[2:18], ca[2:18], carbon[2:18]),
            C.N_CA_C_ANGLE,
            atol=B._N_CA_C_TOLERANCE,
        )

    def test_missing_atom_slots_raise_rather_than_being_skipped(self) -> None:
        # A CA-only structure has nowhere to put an N. Skipping those residues would
        # leave a half-rebuilt structure that looks finished, so it is an error and the
        # message says what to do instead.
        ca = zigzag(10)
        structure = Structure.from_atom_records(
            xyz=ca,
            atom_name=["CA"] * 10,
            element=["C"] * 10,
            residue_name=["ALA"] * 10,
            residue_number=list(range(1, 11)),
            chain_id=["A"] * 10,
        )
        domain = Domain(structure=structure, span=Span(0, 10), kind=DomainKind.IDR)
        with pytest.raises(GeometryError, match="cannot add atoms"):
            place_backbone_for_domain(structure, domain, np.random.default_rng(0))

    def test_too_short_domain_raises(self) -> None:
        structure = self._structure()
        domain = Domain(structure=structure, span=Span(5, 6), kind=DomainKind.IDR)
        with pytest.raises(GeometryError, match="peptide unit needs two"):
            place_backbone_for_domain(structure, domain, np.random.default_rng(0))

    def test_an_interior_span_without_anchors_is_refused(self) -> None:
        # An anchored domain's own first residue is an interior residue of the chain, so it
        # must not be given the arbitrary terminal dihedral -- and the peptide bond to the
        # residue outside the span has to be closed by something. Without an anchor there
        # is nothing to close it onto, so this is refused rather than quietly producing the
        # 1.38 A "peptide bond" it used to.
        structure = self._structure()
        domain = Domain(structure=structure, span=Span(5, 15), kind=DomainKind.IDR)
        before = structure.xyz.copy()
        with pytest.raises(GeometryError, match="no anchor on that side"):
            place_backbone_for_domain(structure, domain, np.random.default_rng(3))
        assert np.array_equal(structure.xyz, before)

    def test_a_span_reaching_the_chain_end_needs_no_anchor_there(self) -> None:
        structure = self._structure()
        domain = Domain(structure=structure, span=Span(0, 20), kind=DomainKind.IDR)
        result = place_backbone_for_domain(structure, domain, np.random.default_rng(3))
        assert result.n_residues == 20
        assert result.steric_clashes == ()

    def test_the_anchors_own_geometry_is_solved_and_reported(self) -> None:
        # The anchor's N is an input to the solve, not an output of it, and the anchor's own
        # N-CA-C is therefore a real measurement rather than an artefact of a placeholder.
        # Both facts are checkable on the mutated structure: the reported angle is the one
        # the structure now has, and it is inside the tolerance. Before, the anchor's C was
        # solved against an N placed from an arbitrary dihedral, and nothing measured the
        # result -- the function returned None.
        structure = self._structure()
        span = Span(5, 15, n_anchor=4, c_anchor=15)
        domain = Domain(structure=structure, span=span, kind=DomainKind.IDR)
        anchor_n_before = reference_atom(structure, "N", 4, 5)[0]
        result = place_backbone_for_domain(structure, domain, np.random.default_rng(3))
        assert result.n_residues == 12  # the span's ten residues plus both anchors
        anchor_n_after = reference_atom(structure, "N", 4, 5)[0]
        assert np.array_equal(anchor_n_before, anchor_n_after)
        measured = float(
            angle(
                anchor_n_after,
                structure.ca_xyz[4],
                reference_atom(structure, "C", 4, 5)[0],
            )
        )
        assert measured == pytest.approx(float(result.n_ca_c_angles[0]), abs=1e-9)
        assert abs(measured - C.N_CA_C_ANGLE) <= B._N_CA_C_TOLERANCE

    def test_a_rebuild_that_would_collide_with_the_rest_is_refused_and_undone(self) -> None:
        # The solve is handed the span and its anchors and nothing else, so it cannot see
        # the rest of the structure. Without the check that follows the write, a span could
        # be rebuilt into a self-consistent backbone that runs straight through a
        # neighbouring domain and the caller would be told nothing at all.
        #
        # The intruder is placed where this exact seeded rebuild is known to put N(10), so
        # the collision is deterministic rather than hoped for.
        probe = self._structure()
        span = Span(8, 15, n_anchor=7, c_anchor=15)
        place_backbone_for_domain(
            probe, Domain(structure=probe, span=span, kind=DomainKind.IDR), np.random.default_rng(4)
        )
        landing = reference_atom(probe, "N", 10, 11)[0]

        structure = self._structure()
        intruder = structure.atom_slice_for_residues(0, 1).start  # residue 0's N
        structure.xyz[intruder] = landing + np.array([0.8, 0.0, 0.0])
        before = structure.xyz.copy()
        with pytest.raises(GeometryError, match="rest of the structure"):
            place_backbone_for_domain(
                structure,
                Domain(structure=structure, span=span, kind=DomainKind.IDR),
                np.random.default_rng(4),
            )
        # Nothing written: a refused rebuild leaves the structure exactly as it was.
        assert np.array_equal(structure.xyz, before)

    def test_real_folded_domain_round_trips(self, dnmt3a: Structure) -> None:
        # The end-to-end case the module exists for: take a real folded domain, keep only
        # its alpha carbons, rebuild N/C/O in place, and check the result against what
        # was there before.
        structure = dnmt3a.copy()
        # 700 onward: the 698-699 junction is a cis-proline, which is a chain break as
        # far as this module is concerned, and using it as an anchor would strand
        # residue 698 on its own.
        span = Span(700, 899, n_anchor=699, c_anchor=899)
        domain = Domain(structure=structure, span=span, kind=DomainKind.FOLDED)
        deposited = {
            name: reference_atom(structure, name, span.start, span.stop) for name in ("N", "C", "O")
        }
        place_backbone_for_domain(structure, domain, np.random.default_rng(20260730))
        for name, before in deposited.items():
            after = reference_atom(structure, name, span.start, span.stop)
            assert rmsd(after, before) < TestAgainstDeposited.LIMITS[name], name
        assert np.array_equal(structure.ca_xyz, dnmt3a.ca_xyz)


# ---------------------------------------------------------------------------
# Non-bonded contacts
#
# Bond lengths and bond angles are purely local: a backbone can satisfy every one of
# them and still fold two carbonyl oxygens through each other. Before these tests
# existed, rebuilding dnmt3a's own deposited CA trace put O(296) and O(305) at 1.282 A,
# shorter than an O=O double bond, and every assertion in this file passed.
# ---------------------------------------------------------------------------


class TestNonBondedContacts:
    @pytest.mark.parametrize(("start", "stop"), TestAgainstDeposited.SEGMENTS)
    def test_deposited_reference_satisfies_the_contact_limits(
        self, dnmt3a: Structure, start: int, stop: int
    ) -> None:
        # The gate is calibrated against reality, not the other way round: the deposited
        # model passes the same limits the rebuild is held to, so a violation in a
        # rebuild is a defect in the rebuild and not an unreasonable threshold.
        xyz, residue, names = [], [], []
        for index in range(start, stop):
            atoms = dnmt3a.atom_slice_for_residues(index, index + 1)
            for atom in range(atoms.start, atoms.stop):
                if str(dnmt3a.atom_name[atom]) in ("N", "CA", "C", "O", "CB"):
                    xyz.append(dnmt3a.xyz[atom])
                    residue.append(index)
                    names.append(str(dnmt3a.atom_name[atom]))
        found = nonbonded_contacts(np.array(xyz), np.array(residue), np.array(names))
        assert found == []

    @pytest.mark.parametrize(("start", "stop"), TestAgainstDeposited.SEGMENTS)
    def test_rebuilt_real_trace_has_no_overlapping_atoms(
        self, dnmt3a: Structure, start: int, stop: int
    ) -> None:
        result = place_backbone(
            dnmt3a.ca_xyz[start:stop],
            dnmt3a.sequence[start:stop],
            np.random.default_rng(20260730),
            sidechains=True,
        )
        assert result_contacts(result, offset=start) == []
        assert result.steric_clashes == ()
        assert result.min_nonbonded_contact >= min(CONTACT_LIMIT.values())

    def test_the_modules_own_measurement_agrees_with_an_independent_one(
        self, dnmt3a: Structure
    ) -> None:
        # BackboneResult.min_nonbonded_contact is the number a caller reports, so it has to
        # be the same number an independent measurement gets. Checked on a rebuild with a
        # chain break in it, because that is where the two could most easily disagree about
        # which pairs count.
        start, stop = 473, 912
        result = place_backbone(
            dnmt3a.ca_xyz[start:stop],
            dnmt3a.sequence[start:stop],
            np.random.default_rng(20260730),
            sidechains=True,
            allow_chain_breaks=True,
        )
        independent = closest_nonbonded(
            result.xyz,
            result.residue_index,
            result.atom_names,
            breaks=result.chain_breaks,
        )
        assert result.min_nonbonded_contact == pytest.approx(independent, abs=1e-9)

    def test_generated_trace_has_no_overlapping_atoms(self) -> None:
        # Self-avoiding at the CA level says nothing about the atoms hanging off it: the
        # measured worst contact on traces like these was 0.234 A, two effectively
        # coincident oxygens, reported as a valid backbone.
        for seed in range(6):
            ca = random_walk(20, np.random.default_rng(seed), angle_max=135.0)
            result = place_backbone(ca, "A" * 20, np.random.default_rng(seed), sidechains=True)
            assert result_contacts(result) == [], seed
            assert result.steric_clashes == ()

    def test_a_clash_is_never_reported_as_a_valid_backbone(self) -> None:
        # Two CA strands run antiparallel 3.0 A apart: the alpha carbons themselves are
        # closer than any two non-bonded carbons can be, so no placement of N, C and O
        # can rescue it and the call must refuse rather than return.
        ca = np.zeros((16, 3))
        ca[:8] = zigzag(8)
        ca[8:] = zigzag(8)[::-1] + np.array([0.0, 3.0, 0.0])
        with pytest.raises(GeometryError, match=r"hard-sphere"):
            place_backbone(ca, "A" * 16, np.random.default_rng(0), allow_chain_breaks=True)

    def test_chain_break_junction_is_not_a_sub_bond_length_contact(self, dnmt3a: Structure) -> None:
        # allow_chain_breaks used to place each segment's terminal atoms from fixed
        # dihedrals with no regard for the next segment: across dnmt3a's two cis-prolines
        # that put C(i) and N(i+1) at 1.124 A, shorter than a C-N triple bond, and the
        # writer then emitted a CONECT record for the pair.
        start, stop = 473, 912
        result = place_backbone(
            dnmt3a.ca_xyz[start:stop],
            dnmt3a.sequence[start:stop],
            np.random.default_rng(20260730),
            sidechains=True,
            allow_chain_breaks=True,
        )
        assert len(result.chain_breaks) == 2
        atoms = split_atoms(result)
        for index in result.chain_breaks:
            junction = float(np.linalg.norm(atoms["N"][index + 1] - atoms["C"][index]))
            assert junction >= C.CA_CLASH_DISTANCE - 0.30, (index, junction)
        assert result_contacts(result, offset=start) == []


# ---------------------------------------------------------------------------
# The seam of an anchored in-place rebuild
# ---------------------------------------------------------------------------


class TestDomainSeam:
    @staticmethod
    def _peptide_bond(structure: Structure, first: int) -> tuple[float, float]:
        """``C(first)-N(first+1)`` length and ``|omega|`` across that bond."""
        carbon = reference_atom(structure, "C", first, first + 1)[0]
        nitrogen = reference_atom(structure, "N", first + 1, first + 2)[0]
        omega = dihedral(structure.ca_xyz[first], carbon, nitrogen, structure.ca_xyz[first + 1])
        return float(np.linalg.norm(nitrogen - carbon)), abs(float(omega))

    def test_seam_peptide_bonds_are_as_good_as_interior_ones(self, dnmt3a: Structure) -> None:
        # The two bonds the anchoring machinery exists to get right were the only two it
        # got wrong: 1.3767 A and 1.2672 A with omega 162.7 deg, while all 198 interior
        # bonds were 1.3290000 A and omega exactly 180.
        structure = dnmt3a.copy()
        span = Span(700, 899, n_anchor=699, c_anchor=899)
        domain = Domain(structure=structure, span=span, kind=DomainKind.FOLDED)
        place_backbone_for_domain(structure, domain, np.random.default_rng(20260730))
        for first in (699, 750, 898):
            length, omega = self._peptide_bond(structure, first)
            assert length == pytest.approx(C.C_N_PEPTIDE_BOND_LENGTH, abs=1e-6), first
            assert omega == pytest.approx(C.OMEGA_TRANS, abs=1e-4), first

    def test_seam_bonds_are_exact_on_a_synthetic_trace(self) -> None:
        ca = zigzag(24, pseudo_angle=126.0)
        built = place_backbone(ca, "A" * 24, np.random.default_rng(0))
        structure = Structure.from_atom_records(**built.atom_records(), source="fixture")
        domain = Domain(
            structure=structure, span=Span(10, 20, n_anchor=9, c_anchor=20), kind=DomainKind.IDR
        )
        place_backbone_for_domain(structure, domain, np.random.default_rng(2))
        for first in (9, 14, 19):
            length, omega = self._peptide_bond(structure, first)
            assert length == pytest.approx(C.C_N_PEPTIDE_BOND_LENGTH, abs=1e-6), first
            assert omega == pytest.approx(C.OMEGA_TRANS, abs=1e-4), first

    def test_anchor_alpha_carbons_and_outside_atoms_are_untouched(self) -> None:
        ca = zigzag(24, pseudo_angle=126.0)
        built = place_backbone(ca, "A" * 24, np.random.default_rng(0))
        structure = Structure.from_atom_records(**built.atom_records(), source="fixture")
        before = structure.xyz.copy()
        domain = Domain(
            structure=structure, span=Span(10, 20, n_anchor=9, c_anchor=20), kind=DomainKind.IDR
        )
        place_backbone_for_domain(structure, domain, np.random.default_rng(2))
        assert np.array_equal(structure.ca_xyz, ca)
        # The anchors keep the atoms that face away from the rebuilt span: the N-side
        # anchor's own N, the C-side anchor's own C and O. Those belong to peptide units
        # outside the span and moving them would break the bonds beyond the anchors.
        outside = [(9, "N"), (20, "C"), (20, "O")]
        for residue, name in outside:
            atoms = structure.atom_slice_for_residues(residue, residue + 1)
            index = atoms.start + int(np.flatnonzero(structure.atom_name[atoms] == name)[0])
            assert np.array_equal(structure.xyz[index], before[index]), (residue, name)
        for residue in range(0, 9):
            atoms = structure.atom_slice_for_residues(residue, residue + 1)
            assert np.array_equal(structure.xyz[atoms], before[atoms]), residue


# ---------------------------------------------------------------------------
# The path from DODO's own CA-trace generators to all-atom output
# ---------------------------------------------------------------------------


#: CA-CA-CA angle ceiling, in degrees, that a generator can sample up to and still expect
#: an all-atom backbone. MEASURED at self-avoidance 5.0 A over walks of 60 and 150 residues,
#: six seeds each: 140 and 143 deg both give 12 of 12 accepted with no failure of any kind,
#: 145.5 -- one degree under the *ideal* ceiling max_reconstructable_ca_angle() = 146.5 --
#: gives 11 of 12 and its one refusal is steric rather than N-CA-C, and BACKBONE_ANGLE_MAX
#: itself gives 9 of 12 with two N-CA-C failures, both at 150 residues.
#:
#: So a trace-level N-CA-C failure appears between 145.5 and 150 deg, not below the ideal
#: ceiling of 146.5: with the tau tolerance re-measured to 14.0 deg, the per-residue ceiling
#: is 160.5 and what refuses these traces is purely the budget competition
#: ca_angle_budget() describes -- two adjacent residues that both want most of a peptide
#: unit's 35.5 deg cannot both have it, so the per-residue ceiling is not a per-*trace* one.
#: At the earlier tolerance of 8.0 deg the same sweep put 3 of 4 N-CA-C failures at 145.5.
#: 140 is kept because it is the value with margin on both counts.
GENERATOR_ANGLE_MAX = 140.0


class TestReconstructableWindow:
    def test_traces_inside_the_generator_window_are_accepted_at_every_length(self) -> None:
        # The working path from a CA-trace generator to all-atom output, which is what did
        # not exist: with default settings the solver refused 100% of the IDR engine's
        # traces, 19 of 62 residues at a time, on N-CA-C angles that were achievable.
        for n in (20, 60, 150):
            for seed in range(2):
                ca = random_walk(
                    n,
                    np.random.default_rng(seed),
                    angle_max=GENERATOR_ANGLE_MAX,
                    self_avoid=5.0,
                )
                result = place_backbone(ca, "A" * n, np.random.default_rng(seed), sidechains=True)
                assert result.max_n_ca_c_deviation <= B._N_CA_C_TOLERANCE, (n, seed)
                assert result_contacts(result) == [], (n, seed)

    def test_acceptance_depends_on_how_tightly_the_trace_packs(self) -> None:
        # The honest limitation, measured and pinned. A trace that lets non-bonded alpha
        # carbons approach each other leaves the atoms hanging off them nowhere to go, and
        # then no assignment of peptide-plane rotations is clash-free. What must never
        # happen is a *silent* pass: every accepted result is clean, and every refusal says
        # which contacts it could not fix.
        #
        # MEASURED at 60 residues, four seeds: 4 of 4 accepted at 5.0 A self-avoidance,
        # 3 of 4 at 4.5, 2 of 4 at 3.8 -- against dnmt3a's own closest non-bonded CA-CA
        # pair of 3.78 A. Density is the reason: these walks put 4.7 to 11.7 non-bonded
        # CA pairs under 4.5 A per 100 residues, where dnmt3a has 3.2.
        rates: dict[float, int] = {}
        for radius in (3.8, 5.0):
            accepted = 0
            for seed in range(4):
                ca = random_walk(
                    60, np.random.default_rng(seed), angle_max=135.0, self_avoid=radius
                )
                try:
                    result = place_backbone(
                        ca, "A" * 60, np.random.default_rng(seed), sidechains=True
                    )
                except GeometryError as refusal:
                    assert "hard-sphere" in str(refusal), (radius, seed)
                    continue
                accepted += 1
                assert result_contacts(result) == [], (radius, seed)
            rates[radius] = accepted
        assert rates[5.0] == 4
        assert 1 <= rates[3.8] <= rates[5.0]

    def test_the_ceiling_accounts_for_the_tolerance(self) -> None:
        # max_reconstructable_ca_angle() is the ceiling for an *ideal* N-CA-C. A caller
        # who accepts a tolerance can go that much further, and the number the caller
        # needs is the one the error message quotes.
        assert max_reconstructable_ca_angle() == pytest.approx(C.N_CA_C_ANGLE + ca_angle_budget())
        assert max_reconstructable_ca_angle(n_ca_c_tolerance=B._N_CA_C_TOLERANCE) == pytest.approx(
            max_reconstructable_ca_angle() + B._N_CA_C_TOLERANCE
        )

    def test_unreconstructable_residues_are_named_before_building(self) -> None:
        # Premise: the walk has to be able to *produce* an impossible residue, which a walk
        # sampled inside the generation window no longer can -- BACKBONE_ANGLE_MAX is capped
        # below the ceiling on purpose. So it samples the measured observed range instead.
        assert STRAINING_ANGLE_MAX > RECONSTRUCTABLE_CEILING
        ca = random_walk(60, np.random.default_rng(0), angle_max=STRAINING_ANGLE_MAX)
        impossible = unreconstructable_ca_angles(ca)
        assert impossible.size > 0
        pseudo = B._ca_pseudo_angles(ca)
        ceiling = max_reconstructable_ca_angle(n_ca_c_tolerance=B._N_CA_C_TOLERANCE)
        assert bool((pseudo[impossible] > ceiling).all())
        # And every residue it does not name is reconstructable, which is what makes the
        # function usable as a pre-flight check rather than a hint.
        keep = np.setdiff1d(np.arange(1, ca.shape[0] - 1), impossible)
        assert bool((pseudo[keep] <= ceiling).all())

    def test_refusal_is_bounded_by_the_residues_that_strain_the_budget(self) -> None:
        # Only a residue that draws on a peptide unit's angular budget beyond what an ideal
        # N-CA-C leaves -- pseudo-angle above max_reconstructable_ca_angle() -- can be
        # refused, either on its own account or by pinning a unit its neighbour also needs
        # (see ca_angle_budget). So the count in the message cannot exceed the number of
        # such residues. It used to: 23 of 62 residues were refused on a trace with 20
        # strained ones and only 8 impossible ones, because the objective treated the
        # tolerance as a preference and traded it away for a fraction of an Angstrom of
        # steric relief.
        #
        # Sampled up to the measured observed maximum rather than to BACKBONE_ANGLE_MAX, for
        # the reason above: the generation window is capped below the ceiling, so a trace
        # drawn from it has no residue that is impossible on its own account and there is no
        # refusal to bound.
        assert STRAINING_ANGLE_MAX > RECONSTRUCTABLE_CEILING
        for seed in range(4):
            ca = random_walk(60, np.random.default_rng(seed), angle_max=STRAINING_ANGLE_MAX)
            pseudo = B._ca_pseudo_angles(ca)
            strained = int((pseudo > max_reconstructable_ca_angle()).sum())
            impossible = unreconstructable_ca_angles(ca).size
            assert impossible > 0, seed
            with pytest.raises(GeometryError) as raised:
                place_backbone(ca, "A" * 60, np.random.default_rng(0))
            refused = int(str(raised.value).split()[0])
            assert impossible <= refused <= strained, (seed, impossible, refused, strained)


# ---------------------------------------------------------------------------
# Reported accuracy, including the tail
# ---------------------------------------------------------------------------


class TestCarbonylTail:
    def test_carbonyl_error_tail_is_reported_and_bounded(self, dnmt3a: Structure) -> None:
        """The carbonyl O error distribution is heavy-tailed, and the tail is the finding.

        MEASURED on dnmt3a, seed 20260730, per segment: O median 0.36-0.44 A, p95
        1.31-2.50 A, max 3.09-4.12 A, and 10 to 22 per cent of carbonyls more than 1 A out.
        Quoting only RMSD (0.73-1.04) and median hides that: an O that is 3 A from where it
        belongs is a carbonyl pointing the wrong way, not a slightly misplaced one, and it
        is the one thing a pure-geometry method cannot pin down without the statistical
        library keyed on the CA pseudo-dihedral that PULCHRA and BBQ use and DODO does not
        ship.

        These are a re-measurement, and what moved is instructive rather than incidental.
        At the earlier ``_N_CA_C_TOLERANCE`` of 8.0 the same four segments gave 11 to 19 per
        cent out and a p95 of 1.71-2.50; re-measuring tau over four whole deposited
        structures widened that gate to 14.0, and the extra six degrees of tau freedom is
        freedom the solver *spends* -- on the steric term, which is the term that competes
        with tau. Segment 480-578 is where it shows: 17.4 per cent of its carbonyls out at
        the old gate, 21.4 per cent at the new one, with the RMSD essentially unchanged
        (0.89 to 1.01 A). That is the trade the wider gate buys, and it belongs in the
        tail statistic rather than hidden behind an RMSD that barely moves.
        """
        worst_p95 = 0.0
        worst_fraction = 0.0
        for start, stop in TestAgainstDeposited.SEGMENTS:
            result = place_backbone(
                dnmt3a.ca_xyz[start:stop],
                dnmt3a.sequence[start:stop],
                np.random.default_rng(20260730),
            )
            atoms = split_atoms(result)
            deposited = reference_atom(dnmt3a, "O", start, stop)
            error = np.linalg.norm(atoms["O"] - deposited, axis=1)
            worst_p95 = max(worst_p95, float(np.nanpercentile(error, 95)))
            worst_fraction = max(worst_fraction, float(np.nanmean(error > 1.0)))
        # The tail is the honest statistic for this distribution. These are pinned just
        # above the measured values so that an improvement in the median cannot mask a
        # regression in the tail.
        assert worst_p95 < 2.60, worst_p95
        assert worst_fraction < 0.22, worst_fraction
        # And it is not a hidden statistic any more: the same tail is quoted in the
        # docstring above, next to the RMSD limits it qualifies.
        assert worst_fraction > 0.05


# ---------------------------------------------------------------------------
# The rotation grid, and the claim made for it
# ---------------------------------------------------------------------------


class TestPlaneGrid:
    def test_grid_refinement_agreement_matches_the_documented_number(
        self, dnmt3a: Structure
    ) -> None:
        # The docstring used to claim 0.005 A RMSD agreement between grids of 90, 120,
        # 180 and 360. Measured, it is an order of magnitude worse, so the claim was
        # replaced by the measurement and pinned here.
        start, stop = 699, 912
        ca = dnmt3a.ca_xyz[start:stop]
        sequence = dnmt3a.sequence[start:stop]
        original = B._PLANE_GRID
        reference: np.ndarray | None = None
        worst = 0.0
        try:
            for grid in (120, 90, 180, 360):
                B._PLANE_GRID = grid
                result = place_backbone(ca, sequence, np.random.default_rng(20260730))
                if reference is None:
                    reference = result.xyz
                    continue
                worst = max(worst, rmsd(result.xyz, reference))
        finally:
            B._PLANE_GRID = original
        assert worst > 0.01, "grids agree better than documented; update the docstring"
        assert worst < B._PLANE_GRID_AGREEMENT_RMSD, worst
