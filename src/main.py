#!/usr/bin/env python3
"""
main.py

Run the complete candidate-region scanning and Primer3 design pipeline.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from primer3_design import Primer3Settings, PrimerPair, design_primers
from region_scanner import (
    CandidateRegion,
    ScannerSettings,
    print_candidates,
    scan_sequence,
)

@dataclass
class RankedResult:
    rank: int
    combined_score: float
    candidate: CandidateRegion
    primer_pair: PrimerPair


def read_pasted_sequence() -> str:
    print(
        "Paste a FASTA or DNA sequence below.\n"
        "Press Enter on a completely blank line when finished.\n"
    )

    lines: list[str] = []

    while True:
        try:
            line = input()
        except EOFError:
            break

        if not line.strip():
            break

        lines.append(line)

    return "\n".join(lines)


def rank_complete_results(
    candidates: list[CandidateRegion],
    primer_settings: Primer3Settings,
    maximum_results: int = 10,
) -> list[RankedResult]:
    results: list[RankedResult] = []

    for candidate in candidates:
        try:
            pairs = design_primers(
                sequence=candidate.sequence,
                sequence_id=f"candidate_{candidate.rank}",
                settings=primer_settings,
            )
        except (ValueError, RuntimeError) as exc:
            print(
                f"Primer3 skipped candidate {candidate.rank}: {exc}"
            )
            continue

        for pair in pairs:
            combined_score = (
                0.40 * candidate.score
                + 0.60 * pair.quality_score
            )

            results.append(
                RankedResult(
                    rank=0,
                    combined_score=combined_score,
                    candidate=candidate,
                    primer_pair=pair,
                )
            )

    results.sort(
        key=lambda item: (
            -item.combined_score,
            item.primer_pair.pair_penalty,
            -item.candidate.score,
            item.candidate.start,
        )
    )

    selected = results[:maximum_results]

    for rank, result in enumerate(selected, start=1):
        result.rank = rank

    return selected


def print_results(results: list[RankedResult]) -> None:
    if not results:
        print("No complete region + primer results were produced.")
        return

    print("\nFinal ranked primer candidates")
    print("-" * 116)
    print(
        f"{'Rank':<6}{'Region':>8}{'Combined':>11}{'Primer':>9}"
        f"{'Product':>10}{'Left Tm':>10}{'Right Tm':>11}"
        f"{'Pair Pen.':>11}{'Pair End':>10}{'Coordinates':>18}"
    )
    print("-" * 116)

    for result in results:
        pair = result.primer_pair
        candidate = result.candidate

        print(
            f"{result.rank:<6}"
            f"{candidate.rank:>8}"
            f"{result.combined_score:>11.2f}"
            f"{pair.quality_score:>9.2f}"
            f"{pair.product_size:>10}"
            f"{pair.forward_tm:>10.2f}"
            f"{pair.reverse_tm:>11.2f}"
            f"{pair.pair_penalty:>11.3f}"
            f"{pair.pair_end:>10.2f}"
            f"{f'{candidate.start}-{candidate.end}':>18}"
        )


def show_result(result: RankedResult) -> None:
    candidate = result.candidate
    pair = result.primer_pair

    print("\n" + "=" * 78)
    print(f"Final candidate #{result.rank}")
    print(f"Combined score: {result.combined_score:.2f}")
    print(
        f"Region: scanner candidate {candidate.rank}, "
        f"coordinates {candidate.start}-{candidate.end}"
    )
    print(
        f"Region score: {candidate.score:.2f} "
        f"({candidate.rating})"
    )
    print(f"Region GC: {candidate.gc_percent:.2f}%")
    print(f"Region notes: {candidate.notes}")

    print("\nPrimer pair")
    print(f"Forward: {pair.forward_sequence}")
    print(f"Reverse: {pair.reverse_sequence}")
    print(
        f"Lengths: {pair.forward_length} / "
        f"{pair.reverse_length} bp"
    )
    print(
        f"Tm: {pair.forward_tm:.2f} / "
        f"{pair.reverse_tm:.2f}°C"
    )
    print(
        f"GC: {pair.forward_gc:.2f} / "
        f"{pair.reverse_gc:.2f}%"
    )
    print(f"Expected product: {pair.product_size} bp")
    print(f"Pair penalty: {pair.pair_penalty:.3f}")
    print(
        f"Self Any: {pair.forward_self_any:.2f} / "
        f"{pair.reverse_self_any:.2f}"
    )
    print(
        f"Self End: {pair.forward_self_end:.2f} / "
        f"{pair.reverse_self_end:.2f}"
    )
    print(
        f"Hairpin: {pair.forward_hairpin:.2f} / "
        f"{pair.reverse_hairpin:.2f}°C"
    )
    print(
        f"Pair Any: {pair.pair_any:.2f}; "
        f"Pair End: {pair.pair_end:.2f}"
    )
    print(f"Primer quality score: {pair.quality_score:.2f}")
    print(f"Primer notes: {pair.notes}")
    print(
        "Suggested initial annealing temperature: "
        f"{pair.suggested_annealing_temperature:.1f}°C"
    )

    print("\nCandidate sequence:")

    for index in range(0, len(candidate.sequence), 80):
        print(candidate.sequence[index:index + 80])


def save_results(
    candidates: list[CandidateRegion],
    results: list[RankedResult],
) -> None:
    candidates_path = Path("scanner_candidates.csv")
    results_path = Path("ranked_primer_pairs.csv")
    fasta_path = Path("ranked_candidate_regions.fasta")

    with candidates_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "rank",
            "start",
            "end",
            "score",
            "rating",
            "gc_percent",
            "n_count",
            "masked_percent",
            "longest_homopolymer",
            "entropy",
            "simple_repeat_percent",
            "repeated_kmer_percent",
            "notes",
            "sequence",
        ])

        for candidate in candidates:
            writer.writerow([
                candidate.rank,
                candidate.start,
                candidate.end,
                f"{candidate.score:.2f}",
                candidate.rating,
                f"{candidate.gc_percent:.2f}",
                candidate.n_count,
                f"{candidate.masked_percent:.2f}",
                candidate.longest_homopolymer,
                f"{candidate.entropy:.4f}",
                f"{candidate.simple_repeat_percent:.2f}",
                f"{candidate.repeated_kmer_percent:.2f}",
                candidate.notes,
                candidate.sequence,
            ])

    with results_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "final_rank",
            "combined_score",
            "region_rank",
            "region_start",
            "region_end",
            "region_score",
            "forward_primer",
            "reverse_primer",
            "forward_tm",
            "reverse_tm",
            "forward_gc",
            "reverse_gc",
            "product_size",
            "pair_penalty",
            "forward_self_any",
            "forward_self_end",
            "forward_hairpin",
            "reverse_self_any",
            "reverse_self_end",
            "reverse_hairpin",
            "pair_any",
            "pair_end",
            "primer_quality_score",
            "primer_notes",
            "suggested_annealing_temperature",
        ])

        for result in results:
            candidate = result.candidate
            pair = result.primer_pair

            writer.writerow([
                result.rank,
                f"{result.combined_score:.2f}",
                candidate.rank,
                candidate.start,
                candidate.end,
                f"{candidate.score:.2f}",
                pair.forward_sequence,
                pair.reverse_sequence,
                f"{pair.forward_tm:.2f}",
                f"{pair.reverse_tm:.2f}",
                f"{pair.forward_gc:.2f}",
                f"{pair.reverse_gc:.2f}",
                pair.product_size,
                f"{pair.pair_penalty:.3f}",
                f"{pair.forward_self_any:.2f}",
                f"{pair.forward_self_end:.2f}",
                f"{pair.forward_hairpin:.2f}",
                f"{pair.reverse_self_any:.2f}",
                f"{pair.reverse_self_end:.2f}",
                f"{pair.reverse_hairpin:.2f}",
                f"{pair.pair_any:.2f}",
                f"{pair.pair_end:.2f}",
                f"{pair.quality_score:.2f}",
                pair.notes,
                f"{pair.suggested_annealing_temperature:.1f}",
            ])

    with fasta_path.open("w", encoding="utf-8") as handle:
        for result in results:
            candidate = result.candidate
            pair = result.primer_pair

            handle.write(
                f">candidate_{result.rank}|"
                f"region={candidate.start}-{candidate.end}|"
                f"forward={pair.forward_sequence}|"
                f"reverse={pair.reverse_sequence}|"
                f"product={pair.product_size}\n"
            )

            for index in range(0, len(candidate.sequence), 80):
                handle.write(candidate.sequence[index:index + 80] + "\n")

    print(f"\nSaved: {candidates_path.resolve()}")
    print(f"Saved: {results_path.resolve()}")
    print(f"Saved: {fasta_path.resolve()}")


def main() -> None:
    print("=" * 78)
    print("PCR Candidate Scanner + Primer3 Pipeline")
    print("=" * 78)

    sequence_text = read_pasted_sequence()

    scanner_settings = ScannerSettings(
        window_size=500,
        step_size=100,
        minimum_score=60.0,
        maximum_candidates=20,
    )

    primer_settings = Primer3Settings(
        number_to_return=5,
        product_size_min=380,
        product_size_max=480,
    )

    try:
        candidates = scan_sequence(
            sequence_text,
            scanner_settings,
        )
    except ValueError as exc:
        print(f"Scanner error: {exc}")
        return

    print("\nTop candidate regions")
    print_candidates(candidates)

    if not candidates:
        return

    print("\nRunning Primer3...")
    results = rank_complete_results(
        candidates,
        primer_settings,
        maximum_results=10,
    )

    print_results(results)

    if not results:
        return

    print(
        "\nCommands:\n"
        "  show 1       Display full details for result #1\n"
        "  list         Show the final ranking again\n"
        "  regions      Show scanner candidates again\n"
        "  save         Save CSV and FASTA results\n"
        "  quit         Exit\n"
    )

    while True:
        command = input("Command: ").strip().lower()

        if command in {"quit", "q", "exit"}:
            print("Finished.")
            break

        if command == "list":
            print_results(results)
            continue

        if command == "regions":
            print_candidates(candidates)
            continue

        if command == "save":
            save_results(candidates, results)
            continue

        if command.startswith("show"):
            parts = command.split()

            if len(parts) != 2 or not parts[1].isdigit():
                print("Use a command such as: show 1")
                continue

            rank = int(parts[1])

            if not 1 <= rank <= len(results):
                print(
                    f"Choose a result from 1 to {len(results)}."
                )
                continue

            show_result(results[rank - 1])
            continue

        print("Unknown command. Use show 1, list, regions, save, or quit.")


if __name__ == "__main__":
    main()