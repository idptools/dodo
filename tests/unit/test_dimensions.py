"""Tests for target dimension prediction.

This is the module whose absence made the first v2 attempt unable to do the thing DODO
exists for. It had no sparrow import anywhere, and substituted two *disagreeing* hardcoded
placeholders (``1.0 * n_residues`` in the IDR builder, ``1.4 * n_residues`` for folded-domain
spacing), so the builder was working against its own anchors.

The tests that matter most here are the length-independence ones. v1's modes were Angstroms
per residue, i.e. linear in N, while real IDRs scale as roughly N^0.55 -- so the whole point
of the rewrite is that a mode now means the same thing at every chain length.
"""

from __future__ import annotations

import json
import warnings

import pytest

from dodo.constants import MODES, contour_length
from dodo.construct.dimensions import (
    DimensionTarget,
    albatross_available,
    albatross_short_sequence_threshold,
    predict_end_to_end,
    predict_radius_of_gyration,
    predict_scaling_exponent,
    target_dimensions,
)
from dodo.exceptions import MissingDependencyError, UnsatisfiableTargetError

#: A generic disordered composition: polar and charged, low hydrophobic.
IDR = "SGQNTEKDRSGQNTPKAE" * 8  # 144 residues

requires_albatross = pytest.mark.skipif(
    not albatross_available(), reason="sparrow/ALBATROSS is not installed"
)


class TestPredictEndToEnd:
    def test_analytical_path_needs_no_sparrow(self) -> None:
        value, source = predict_end_to_end(IDR, prefer_albatross=False)
        assert source == "analytical"
        assert value > 0

    def test_explicit_analytical_choice_does_not_warn(self) -> None:
        """Only an *unintentional* downgrade warns; asking for it is not a downgrade."""
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            predict_end_to_end(IDR, prefer_albatross=False)

    @requires_albatross
    def test_albatross_is_used_when_available(self) -> None:
        _, source = predict_end_to_end(IDR)
        assert source == "albatross"

    @requires_albatross
    def test_albatross_and_analytical_are_the_same_order(self) -> None:
        """A sanity bound on the fallback, not a tight agreement claim.

        Measured across six compositional classes the fallback runs 0.6-0.95x of the
        prediction for genuine IDR sequences, so a factor-of-two bound is the honest
        assertion. See the FLORY_RE_PREFACTOR docstring for the per-class numbers.
        """
        predicted, _ = predict_end_to_end(IDR)
        analytical, _ = predict_end_to_end(IDR, prefer_albatross=False)
        assert 0.5 < analytical / predicted < 2.0

    @requires_albatross
    @pytest.mark.parametrize("n", [5, 10, 20, 34, 35, 36, 50])
    def test_short_sequences_do_not_need_special_handling(self, n: int) -> None:
        """Sparrow handles sequences below its 35-residue threshold itself.

        Most loops DODO rebuilds are shorter than 35 residues, so this path matters. The
        concern was that ``use_scaled`` would have to be toggled below the threshold;
        measurement shows sparrow already ignores the flag there and takes the right path.
        """
        value, source = predict_end_to_end(IDR[:n])
        assert source == "albatross"
        assert value > 0

    @pytest.mark.parametrize("bad", ["", "   ", "ACDE FGH", "ACDE-FGH", "ACDE1"])
    def test_malformed_sequences_are_rejected(self, bad: str) -> None:
        with pytest.raises(ValueError):
            predict_end_to_end(bad)

    def test_sequence_is_case_insensitive(self) -> None:
        upper, _ = predict_end_to_end(IDR.upper(), prefer_albatross=False)
        lower, _ = predict_end_to_end(IDR.lower(), prefer_albatross=False)
        assert upper == pytest.approx(lower)


class TestModes:
    def test_every_mode_resolves(self) -> None:
        for mode in MODES:
            target = target_dimensions(IDR, mode=mode, prefer_albatross=False)
            assert target.end_to_end > 0
            assert target.mode == mode

    def test_modes_are_ordered_by_compactness(self) -> None:
        ordered = [
            "super_compact",
            "compact",
            "normal",
            "expanded",
            "super_expanded",
            "max_expansion",
        ]
        values = [
            target_dimensions(IDR, mode=m, prefer_albatross=False).end_to_end for m in ordered
        ]
        assert values == sorted(values)

    def test_predicted_returns_the_prediction_unscaled(self) -> None:
        target = target_dimensions(IDR, mode="predicted", prefer_albatross=False)
        assert target.factor == 1.0
        assert target.end_to_end == pytest.approx(target.predicted)

    def test_unknown_mode_is_rejected_with_the_valid_list(self) -> None:
        with pytest.raises(ValueError, match="Valid modes are"):
            target_dimensions(IDR, mode="very_squished")

    @pytest.mark.parametrize("mode", ["compact", "normal", "expanded"])
    def test_factor_is_length_independent(self, mode: str) -> None:
        """The point of the rewrite: a mode means the same thing at every length.

        v1 applied a fixed Angstroms-per-residue multiplier, so "expanded" meant something
        different at N=50 than at N=800. Here the *factor* is constant and the per-residue
        value is what varies -- the opposite of v1, and the correct way round, because real
        end-to-end distance scales sublinearly in N.
        """
        factors = set()
        per_residue = []
        for n in (50, 200, 800):
            target = target_dimensions((IDR * 10)[:n], mode=mode, prefer_albatross=False)
            factors.add(target.factor)
            per_residue.append(target.per_residue)
        assert len(factors) == 1, "the multiplier must not depend on length"
        assert per_residue[0] > per_residue[-1], "per-residue value must fall with length"

    def test_per_residue_value_falls_with_length(self) -> None:
        """Direct expression of sublinear scaling, which is what v1 got wrong."""
        short = target_dimensions((IDR * 10)[:50], prefer_albatross=False)
        long = target_dimensions((IDR * 10)[:800], prefer_albatross=False)
        assert short.per_residue > 2 * long.per_residue


class TestClamping:
    def test_impossible_target_is_clamped(self) -> None:
        """A 10-residue chain cannot span what max_expansion asks for."""
        target = target_dimensions("GSGSGSGSGS", mode="max_expansion", prefer_albatross=False)
        assert target.clamped
        assert target.end_to_end <= contour_length(10)
        assert target.end_to_end < target.predicted * target.factor

    def test_clamp_can_be_refused(self) -> None:
        with pytest.raises(UnsatisfiableTargetError) as excinfo:
            target_dimensions(
                "GSGSGSGSGS", mode="max_expansion", prefer_albatross=False, clamp=False
            )
        # Both numbers must be reported so a caller can decide what to do.
        assert excinfo.value.target is not None
        assert excinfo.value.achievable is not None

    def test_reasonable_target_is_not_clamped(self) -> None:
        assert not target_dimensions(IDR, mode="normal", prefer_albatross=False).clamped

    def test_target_never_exceeds_contour_length(self) -> None:
        for n in (5, 10, 50, 200, 1000):
            for mode in MODES:
                target = target_dimensions((IDR * 20)[:n], mode=mode, prefer_albatross=False)
                assert target.end_to_end <= contour_length(n), f"N={n} mode={mode}"


class TestProvenance:
    def test_target_records_how_it_was_derived(self) -> None:
        target = target_dimensions(IDR, mode="expanded", prefer_albatross=False)
        assert isinstance(target, DimensionTarget)
        assert target.n_residues == len(IDR)
        assert target.mode == "expanded"
        assert target.factor == MODES["expanded"]
        assert target.source == "analytical"
        assert target.end_to_end == pytest.approx(target.predicted * target.factor)

    def test_str_is_informative(self) -> None:
        text = str(target_dimensions(IDR, mode="compact", prefer_albatross=False))
        assert "compact" in text
        assert "analytical" in text
        assert "residues" in text

    def test_str_mentions_clamping_when_it_happened(self) -> None:
        text = str(target_dimensions("GSGSGSGSGS", mode="max_expansion", prefer_albatross=False))
        assert "clamped" in text


class TestRadiusOfGyration:
    def test_returns_a_positive_value(self) -> None:
        value, source = predict_radius_of_gyration(IDR)
        assert value > 0
        assert source in {"albatross", "analytical"}

    @requires_albatross
    def test_rg_is_smaller_than_re(self) -> None:
        """A basic polymer sanity check: Rg < Re for any real chain."""
        rg, _ = predict_radius_of_gyration(IDR)
        re, _ = predict_end_to_end(IDR)
        assert rg < re


class TestScalingExponent:
    @requires_albatross
    def test_idr_composition_gives_a_coil_exponent(self) -> None:
        assert 0.45 < predict_scaling_exponent(IDR) < 0.75

    @requires_albatross
    def test_hydrophobic_composition_is_detectably_collapsed(self) -> None:
        """The diagnostic that says "this sequence is not really an IDR".

        A hydrophobic-rich sequence is predicted to collapse into a globule, so rebuilding it
        as an expanded coil would be actively wrong. Measurement puts such a sequence near
        0.20-0.38, well below the 0.5-0.6 of a solvent-expanded chain.
        """
        collapsed = predict_scaling_exponent("ADEFGHIKLMNPQRSTVWY" * 8)
        expanded = predict_scaling_exponent(IDR)
        assert collapsed < expanded

    def test_raises_without_sparrow_rather_than_inventing_an_answer(self) -> None:
        """There is no meaningful fallback: a length-only law has a fixed exponent.

        Returning that fixed value would answer a different question than the caller asked,
        which is the class of silent-wrong-answer this rewrite exists to eliminate.
        """
        if albatross_available():
            pytest.skip("sparrow is installed; the raising path is covered elsewhere")
        with pytest.raises(MissingDependencyError, match="sparrow"):
            predict_scaling_exponent(IDR)


class TestThresholdReporting:
    def test_threshold_is_read_from_sparrow_when_present(self) -> None:
        assert albatross_short_sequence_threshold() == 35


class TestPredictionCacheCap:
    """The prediction cache is persistent, so it needs a ceiling it cannot be argued past.

    A size cap rather than an entry count: disk is the thing being protected, and an entry count
    only bounds bytes through an assumed cost per entry, which drifted once already (documented at
    90 bytes, measured at 116). Measuring the serialized payload makes the guarantee hold whatever
    a key or value costs.
    """

    def test_cache_is_trimmed_to_stay_under_the_cap(self, monkeypatch: pytest.MonkeyPatch) -> None:
        from dodo.construct import dimensions

        monkeypatch.setattr(dimensions, "_PREDICTION_CACHE_MAX_BYTES", 50_000)
        cache = {f"v1|end_to_end|{i:064x}": float(i) for i in range(2000)}
        assert len(json.dumps(cache).encode("utf-8")) > 50_000, "fixture must exceed the cap"

        payload = dimensions._trim_to_cap(cache)

        assert len(payload.encode("utf-8")) <= 50_000
        assert json.loads(payload) == cache, "what is written must match what is kept in memory"

    def test_it_evicts_oldest_first(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """FIFO, so the entries a user is actively hitting are the ones that survive."""
        from dodo.construct import dimensions

        monkeypatch.setattr(dimensions, "_PREDICTION_CACHE_MAX_BYTES", 50_000)
        cache = {f"v1|end_to_end|{i:064x}": float(i) for i in range(2000)}
        dimensions._trim_to_cap(cache)

        assert f"v1|end_to_end|{1999:064x}" in cache, "the newest entry must survive"
        assert f"v1|end_to_end|{0:064x}" not in cache, "the oldest must be the one dropped"

    def test_a_small_cache_is_left_alone(self) -> None:
        """The common case must not pay for the cap: no eviction, no re-serialization loop."""
        from dodo.construct import dimensions

        cache = {"v1|end_to_end|abc": 42.0}
        payload = dimensions._trim_to_cap(cache)

        assert cache == {"v1|end_to_end|abc": 42.0}
        assert json.loads(payload) == cache

    def test_the_cap_is_far_above_a_real_workload(self) -> None:
        """Rebuilding the whole human proteome produced 9,172 entries and 1.07 MB, measured.

        The cap must stay comfortably above that or it stops being a backstop and starts evicting
        during ordinary use.
        """
        from dodo.construct import dimensions

        measured_human_proteome_bytes = 1_066_900
        assert 5 * measured_human_proteome_bytes < dimensions._PREDICTION_CACHE_MAX_BYTES
