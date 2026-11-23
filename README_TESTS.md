"""
Tests for src/primer3_design.py

The Primer3 engine is mocked for most tests so the results are deterministic.
primer3-py must still be installed because primer3_design imports it.
"""

from __future__ import annotations

import pytest

import primer3_design
from primer3_design import (
    Primer3Settings,
    clean_sequence,
    design_primers,
    score_primer_pair,
)


def test_clean_sequence_removes_whitespace_and_uppercases() -> None:
    assert clean_sequence("acgt\nacgt  acgt") == "ACGTACGTACGT"


def test_clean_sequence_rejects_invalid_characters() -> None:
    with pytest.raises(ValueError, match="unsupported character"):
        clean_sequence("ACGT-X")


def test_clean_primer_pair_receives_high_score() -> None:
    score, notes = score_primer_pair(
        pair_penalty=0.2,
        forward_tm=60.0,
        reverse_tm=60.3,
        forward_gc=50.0,
        reverse_gc=50.0,
        forward_self_any=5.0,
        forward_self_end=2.0,
        forward_hairpin=20.0,
        reverse_self_any=5.0,
        reverse_self_end=2.0,
        reverse_hairpin=20.0,
        pair_any=5.0,
        pair_end=2.0,
    )

    assert score >= 90.0
    assert notes == "Clean Primer3 profile"


def test_large_tm_difference_reduces_score() -> None:
    good_score, _ = score_primer_pair(
        pair_penalty=0.2,
        forward_tm=60.0,
        reverse_tm=60.5,
        forward_gc=50.0,
        reverse_gc=50.0,
        forward_self_any=5.0,
        forward_self_end=2.0,
        forward_hairpin=20.0,
        reverse_self_any=5.0,
        reverse_self_end=2.0,
        reverse_hairpin=20.0,
        pair_any=5.0,
        pair_end=2.0,
    )

    poor_score, notes = score_primer_pair(
        pair_penalty=0.2,
        forward_tm=60.0,
        reverse_tm=66.0,
        forward_gc=50.0,
        reverse_gc=50.0,
        forward_self_any=5.0,
        forward_self_end=2.0,
        forward_hairpin=20.0,
        reverse_self_any=5.0,
        reverse_self_end=2.0,
        reverse_hairpin=20.0,
        pair_any=5.0,
        pair_end=2.0,
    )

    assert poor_score < good_score
    assert "Tm difference" in notes


def test_high_hairpin_value_reduces_score() -> None:
    normal_score, _ = score_primer_pair(
        pair_penalty=0.2,
        forward_tm=60.0,
        reverse_tm=60.0,
        forward_gc=50.0,
        reverse_gc=50.0,
        forward_self_any=5.0,
        forward_self_end=2.0,
        forward_hairpin=20.0,
        reverse_self_any=5.0,
        reverse_self_end=2.0,
        reverse_hairpin=20.0,
        pair_any=5.0,
        pair_end=2.0,
    )

    hairpin_score, notes = score_primer_pair(
        pair_penalty=0.2,
        forward_tm=60.0,
        reverse_tm=60.0,
        forward_gc=50.0,
        reverse_gc=50.0,
        forward_self_any=5.0,
        forward_self_end=2.0,
        forward_hairpin=52.0,
        reverse_self_any=5.0,
        reverse_self_end=2.0,
        reverse_hairpin=20.0,
        pair_any=5.0,
        pair_end=2.0,
    )

    assert hairpin_score < normal_score
    assert "forward hairpin" in notes


def test_design_primers_rejects_template_shorter_than_product_minimum() -> None:
    settings = Primer3Settings(product_size_min=380)

    with pytest.raises(ValueError, match="minimum requested product size"):
        design_primers(
            sequence="ACGT" * 50,
            settings=settings,
        )


def test_design_primers_parses_and_ranks_mocked_results(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fake_result = {
        "PRIMER_PAIR_NUM_RETURNED": 2,

        "PRIMER_LEFT_0_SEQUENCE": "AAAACCCCGGGGTTTTAAAA",
        "PRIMER_RIGHT_0_SEQUENCE": "TTTTGGGGCCCCAAAATTTT",
        "PRIMER_LEFT_0": [10, 20],
        "PRIMER_RIGHT_0": [430, 20],
        "PRIMER_LEFT_0_TM": 60.0,
        "PRIMER_RIGHT_0_TM": 60.2,
        "PRIMER_LEFT_0_GC_PERCENT": 50.0,
        "PRIMER_RIGHT_0_GC_PERCENT": 50.0,
        "PRIMER_PAIR_0_PRODUCT_SIZE": 440,
        "PRIMER_PAIR_0_PENALTY": 0.2,
        "PRIMER_LEFT_0_PENALTY": 0.1,
        "PRIMER_RIGHT_0_PENALTY": 0.1,
        "PRIMER_LEFT_0_SELF_ANY_TH": 5.0,
        "PRIMER_LEFT_0_SELF_END_TH": 2.0,
        "PRIMER_LEFT_0_HAIRPIN_TH": 20.0,
        "PRIMER_RIGHT_0_SELF_ANY_TH": 5.0,
        "PRIMER_RIGHT_0_SELF_END_TH": 2.0,
        "PRIMER_RIGHT_0_HAIRPIN_TH": 20.0,
        "PRIMER_PAIR_0_COMPL_ANY_TH": 5.0,
        "PRIMER_PAIR_0_COMPL_END_TH": 2.0,

        "PRIMER_LEFT_1_SEQUENCE": "ACGTACGTACGTACGTACGT",
        "PRIMER_RIGHT_1_SEQUENCE": "TGCATGCATGCATGCATGCA",
        "PRIMER_LEFT_1": [20, 20],
        "PRIMER_RIGHT_1": [420, 20],
        "PRIMER_LEFT_1_TM": 58.0,
        "PRIMER_RIGHT_1_TM": 63.0,
        "PRIMER_LEFT_1_GC_PERCENT": 50.0,
        "PRIMER_RIGHT_1_GC_PERCENT": 50.0,
        "PRIMER_PAIR_1_PRODUCT_SIZE": 420,
        "PRIMER_PAIR_1_PENALTY": 2.0,
        "PRIMER_LEFT_1_PENALTY": 1.0,
        "PRIMER_RIGHT_1_PENALTY": 1.0,
        "PRIMER_LEFT_1_SELF_ANY_TH": 35.0,
        "PRIMER_LEFT_1_SELF_END_TH": 12.0,
        "PRIMER_LEFT_1_HAIRPIN_TH": 46.0,
        "PRIMER_RIGHT_1_SELF_ANY_TH": 35.0,
        "PRIMER_RIGHT_1_SELF_END_TH": 12.0,
        "PRIMER_RIGHT_1_HAIRPIN_TH": 46.0,
        "PRIMER_PAIR_1_COMPL_ANY_TH": 35.0,
        "PRIMER_PAIR_1_COMPL_END_TH": 12.0,
    }

    def fake_design_primers(sequence_args, global_args):
        assert sequence_args["SEQUENCE_ID"] == "test_candidate"
        assert len(sequence_args["SEQUENCE_TEMPLATE"]) == 500
        assert global_args["PRIMER_NUM_RETURN"] == 5
        return fake_result

    monkeypatch.setattr(
        primer3_design.primer3.bindings,
        "design_primers",
        fake_design_primers,
    )

    pairs = design_primers(
        sequence="ACGT" * 125,
        sequence_id="test_candidate",
    )

    assert len(pairs) == 2
    assert pairs[0].rank == 1
    assert pairs[0].forward_sequence == "AAAACCCCGGGGTTTTAAAA"
    assert pairs[0].product_size == 440
    assert pairs[0].quality_score > pairs[1].quality_score
    assert pairs[0].suggested_annealing_temperature == pytest.approx(58.0)
