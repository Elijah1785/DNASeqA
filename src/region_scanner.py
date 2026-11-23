#!/usr/bin/env python3
"""
500bp_scanner.py

Scan a DNA sequence for high-quality overlapping 500 bp candidate regions.

This module is designed to be imported by main.py, but it can also be run
directly for region scanning only.
"""

from __future__ import annotations

import math
import re
from collections import Counter
from dataclasses import dataclass


@dataclass
class ScannerSettings:
    window_size: int = 500
    step_size: int = 100
    minimum_score: float = 60.0
    maximum_candidates: int = 20


@dataclass
class CandidateRegion:
    rank: int
    start: int
    end: int
    score: float
    rating: str
    gc_percent: float
    n_count: int
    masked_percent: float
    longest_homopolymer: int
    entropy: float
    simple_repeat_percent: float
    repeated_kmer_percent: float
    notes: str
    sequence: str


def clean_sequence(text: str) -> str:
    """
    Clean pasted FASTA or plain DNA text.

    FASTA headers, spaces, numbers, and punctuation are removed. Letter case is
    preserved so lowercase masking can be measured by the scanner.
    """
    lines = []

    for line in text.splitlines():
        stripped = line.strip()

        if not stripped or stripped.startswith(">"):
            continue

        lines.append(re.sub(r"[^A-Za-z]", "", stripped))

    sequence = "".join(lines)

    if not sequence:
        raise ValueError("No DNA sequence was detected.")

    invalid = sorted(set(sequence.upper()) - set("ACGTN"))

    if invalid:
        raise ValueError(
            "Unsupported sequence character(s): " + ", ".join(invalid)
        )

    return sequence


def gc_percent(sequence: str) -> float:
    valid = [base for base in sequence.upper() if base in "ACGT"]

    if not valid:
        return 0.0

    return 100.0 * (valid.count("G") + valid.count("C")) / len(valid)


def shannon_entropy(sequence: str) -> float:
    valid = [base for base in sequence.upper() if base in "ACGT"]

    if not valid:
        return 0.0

    counts = Counter(valid)
    total = len(valid)

    return -sum(
        (count / total) * math.log2(count / total)
        for count in counts.values()
    )


def longest_homopolymer(sequence: str) -> int:
    longest = 0

    for match in re.finditer(r"([ACGTacgt])\1*", sequence):
        longest = max(longest, len(match.group(0)))

    return longest


def simple_repeat_percent(sequence: str) -> float:
    covered: set[int] = set()

    patterns = [
        re.compile(r"((?P<di>[ACGT]{2})(?P=di){3,})", re.IGNORECASE),
        re.compile(r"((?P<tri>[ACGT]{3})(?P=tri){2,})", re.IGNORECASE),
    ]

    for pattern in patterns:
        for match in pattern.finditer(sequence):
            covered.update(range(match.start(), match.end()))

    return 100.0 * len(covered) / len(sequence) if sequence else 0.0


def repeated_kmer_percent(sequence: str, k: int = 20) -> float:
    upper = sequence.upper()

    if len(upper) < k:
        return 100.0

    positions_by_kmer: dict[str, list[int]] = {}

    for index in range(len(upper) - k + 1):
        kmer = upper[index:index + k]

        if any(base not in "ACGT" for base in kmer):
            continue

        positions_by_kmer.setdefault(kmer, []).append(index)

    covered: set[int] = set()

    for positions in positions_by_kmer.values():
        if len(positions) > 1:
            for start in positions:
                covered.update(range(start, start + k))

    return 100.0 * len(covered) / len(sequence)


def evaluate_region(sequence: str) -> tuple[
    float, str, float, int, float, int, float, float, float, str
]:
    upper = sequence.upper()
    length = len(sequence)

    n_count = upper.count("N")
    masked_count = sum(
        1 for base in sequence
        if base.isalpha() and base.islower()
    )

    masked_percent = 100.0 * masked_count / length if length else 0.0
    gc = gc_percent(sequence)
    entropy = shannon_entropy(sequence)
    longest_run = longest_homopolymer(sequence)
    simple_repeat = simple_repeat_percent(sequence)
    repeated_kmers = repeated_kmer_percent(sequence)

    score = 100.0
    reasons: list[str] = []

    if n_count:
        score -= min(50.0, n_count * 5.0)
        reasons.append(f"{n_count} N base(s)")

    if masked_percent:
        score -= min(45.0, masked_percent * 1.5)
        reasons.append(f"{masked_percent:.1f}% lowercase/masked")

    if 45.0 <= gc <= 55.0:
        pass
    elif 40.0 <= gc < 45.0 or 55.0 < gc <= 60.0:
        score -= 5.0
        reasons.append(f"acceptable GC ({gc:.1f}%)")
    elif 35.0 <= gc < 40.0 or 60.0 < gc <= 65.0:
        score -= 18.0
        reasons.append(f"poor GC ({gc:.1f}%)")
    else:
        score -= 35.0
        reasons.append(f"extreme GC ({gc:.1f}%)")

    if longest_run >= 10:
        score -= 30.0
        reasons.append(f"homopolymer length {longest_run}")
    elif longest_run >= 8:
        score -= 18.0
        reasons.append(f"homopolymer length {longest_run}")
    elif longest_run >= 7:
        score -= 10.0
        reasons.append(f"homopolymer length {longest_run}")
    elif longest_run >= 6:
        score -= 4.0
        reasons.append(f"homopolymer length {longest_run}")

    if entropy < 1.60:
        score -= 30.0
        reasons.append(f"very low entropy ({entropy:.3f})")
    elif entropy < 1.75:
        score -= 18.0
        reasons.append(f"low entropy ({entropy:.3f})")
    elif entropy < 1.85:
        score -= 8.0
        reasons.append(f"moderate entropy ({entropy:.3f})")

    if simple_repeat >= 20.0:
        score -= 30.0
        reasons.append(f"{simple_repeat:.1f}% simple repeats")
    elif simple_repeat >= 10.0:
        score -= 18.0
        reasons.append(f"{simple_repeat:.1f}% simple repeats")
    elif simple_repeat >= 4.0:
        score -= 8.0
        reasons.append(f"{simple_repeat:.1f}% simple repeats")

    if repeated_kmers >= 30.0:
        score -= 25.0
        reasons.append(f"{repeated_kmers:.1f}% repeated 20-mers")
    elif repeated_kmers >= 15.0:
        score -= 15.0
        reasons.append(f"{repeated_kmers:.1f}% repeated 20-mers")
    elif repeated_kmers >= 5.0:
        score -= 6.0
        reasons.append(f"{repeated_kmers:.1f}% repeated 20-mers")

    score = max(0.0, min(100.0, score))

    hard_failure = (
        n_count > 0
        or masked_percent > 10.0
        or gc < 35.0
        or gc > 65.0
        or longest_run >= 10
        or entropy < 1.60
    )

    if hard_failure or score < 60:
        rating = "BAD"
    elif score < 80:
        rating = "REVIEW"
    elif score < 90:
        rating = "GOOD"
    else:
        rating = "EXCELLENT"

    notes = "; ".join(reasons) if reasons else "No composition concerns"

    return (
        score,
        rating,
        gc,
        n_count,
        masked_percent,
        longest_run,
        entropy,
        simple_repeat,
        repeated_kmers,
        notes,
    )


def scan_sequence(
    sequence: str,
    settings: ScannerSettings | None = None,
) -> list[CandidateRegion]:
    """Return ranked candidate regions from a DNA sequence."""
    settings = settings or ScannerSettings()
    cleaned = clean_sequence(sequence)

    if len(cleaned) < settings.window_size:
        raise ValueError(
            f"Sequence is {len(cleaned)} bp; at least "
            f"{settings.window_size} bp are required."
        )

    candidates: list[CandidateRegion] = []

    for zero_based_start in range(
        0,
        len(cleaned) - settings.window_size + 1,
        settings.step_size,
    ):
        window = cleaned[
            zero_based_start:zero_based_start + settings.window_size
        ]

        metrics = evaluate_region(window)

        if metrics[0] < settings.minimum_score:
            continue

        candidates.append(
            CandidateRegion(
                rank=0,
                start=zero_based_start + 1,
                end=zero_based_start + settings.window_size,
                score=metrics[0],
                rating=metrics[1],
                gc_percent=metrics[2],
                n_count=metrics[3],
                masked_percent=metrics[4],
                longest_homopolymer=metrics[5],
                entropy=metrics[6],
                simple_repeat_percent=metrics[7],
                repeated_kmer_percent=metrics[8],
                notes=metrics[9],
                sequence=window.upper(),
            )
        )

    candidates.sort(
        key=lambda candidate: (
            -candidate.score,
            candidate.masked_percent,
            candidate.simple_repeat_percent,
            candidate.repeated_kmer_percent,
            abs(candidate.gc_percent - 50.0),
            candidate.start,
        )
    )

    selected = candidates[:settings.maximum_candidates]

    for rank, candidate in enumerate(selected, start=1):
        candidate.rank = rank

    return selected


def print_candidates(candidates: list[CandidateRegion]) -> None:
    if not candidates:
        print("No candidate regions passed the scanner settings.")
        return

    print(
        f"{'Rank':<6}{'Start':>9}{'End':>9}{'Score':>9}"
        f"{'Rating':>13}{'GC%':>9}{'Longest run':>14}"
    )
    print("-" * 69)

    for candidate in candidates:
        print(
            f"{candidate.rank:<6}"
            f"{candidate.start:>9}"
            f"{candidate.end:>9}"
            f"{candidate.score:>9.2f}"
            f"{candidate.rating:>13}"
            f"{candidate.gc_percent:>9.2f}"
            f"{candidate.longest_homopolymer:>14}"
        )


if __name__ == "__main__":
    print("Paste a FASTA or DNA sequence. Press Enter on a blank line:")
    lines: list[str] = []

    while True:
        try:
            line = input()
        except EOFError:
            break

        if not line.strip():
            break

        lines.append(line)

    try:
        results = scan_sequence("\n".join(lines))
        print_candidates(results)
    except ValueError as exc:
        print(f"Error: {exc}")