"""
Tests for src/region_scanner.py
"""

from __future__ import annotations

import pytest

from region_scanner import (
    ScannerSettings,
    clean_sequence,
    evaluate_region,
    gc_percent,
    longest_homopolymer,
    scan_sequence,
    shannon_entropy,
)


def balanced_sequence(length: int = 1000) -> str:
    """Return a predictable, balanced A/C/G/T test sequence."""
    motif = "ACGT"
    repeats = (length // len(motif)) + 1
    return (motif * repeats)[:length]


def test_clean_sequence_removes_fasta_header_and_spaces() -> None:
    raw = """>example_sequence
    ACGT ACGT
    ACGT
    """

    assert clean_sequence(raw) == "ACGTACGTACGT"


def test_clean_sequence_preserves_lowercase_masking() -> None:
    raw = ">masked\nACGTacgt"

    assert clean_sequence(raw) == "ACGTacgt"


def test_clean_sequence_rejects_invalid_letters() -> None:
    with pytest.raises(ValueError, match="Unsupported sequence character"):
        clean_sequence(">bad\nACGTXYZ")


def test_gc_percent_is_calculated_correctly() -> None:
    assert gc_percent("AACCGGTT") == pytest.approx(50.0)


def test_longest_homopolymer_is_detected() -> None:
    assert longest_homopolymer("ACGTTTTTTACG") == 6


def test_balanced_sequence_has_high_entropy() -> None:
    entropy = shannon_entropy(balanced_sequence(500))

    assert entropy == pytest.approx(2.0, abs=0.01)


def test_scan_sequence_returns_500_bp_windows() -> None:
    settings = ScannerSettings(
        window_size=500,
        step_size=100,
        minimum_score=0.0,
        maximum_candidates=20,
    )

    candidates = scan_sequence(
        balanced_sequence(1000),
        settings,
    )

    assert len(candidates) == 6
    assert all(len(candidate.sequence) == 500 for candidate in candidates)


def test_scan_sequence_uses_one_based_coordinates() -> None:
    settings = ScannerSettings(
        window_size=500,
        step_size=100,
        minimum_score=0.0,
        maximum_candidates=20,
    )

    candidates = scan_sequence(
        balanced_sequence(1000),
        settings,
    )

    coordinate_pairs = {
        (candidate.start, candidate.end)
        for candidate in candidates
    }

    assert (1, 500) in coordinate_pairs
    assert (501, 1000) in coordinate_pairs


def test_short_sequence_raises_clear_error() -> None:
    settings = ScannerSettings(window_size=500)

    with pytest.raises(ValueError, match="at least 500 bp"):
        scan_sequence("ACGT" * 100, settings)


def test_region_with_ns_receives_penalty() -> None:
    clean_region = balanced_sequence(500)
    region_with_ns = "N" * 10 + clean_region[10:]

    clean_score = evaluate_region(clean_region)[0]
    n_score = evaluate_region(region_with_ns)[0]

    assert n_score < clean_score


def test_low_complexity_region_scores_below_balanced_region() -> None:
    balanced_score = evaluate_region(balanced_sequence(500))[0]
    low_complexity_score = evaluate_region("A" * 500)[0]

    assert low_complexity_score < balanced_score