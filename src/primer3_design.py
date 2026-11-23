#!/usr/bin/env python3
"""
primer3_design.py

Reusable Primer3 design and primer-ranking module.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Iterable

try:
    import primer3  # type: ignore
except ImportError as exc:
    raise ImportError(
        "primer3-py is not installed. Install it with:\n"
        "    py -m pip install primer3-py"
    ) from exc


@dataclass(frozen=True)
class Primer3Settings:
    number_to_return: int = 5

    primer_min_size: int = 18
    primer_opt_size: int = 20
    primer_max_size: int = 25

    primer_min_tm: float = 57.0
    primer_opt_tm: float = 60.0
    primer_max_tm: float = 63.0

    primer_min_gc: float = 40.0
    primer_max_gc: float = 60.0

    product_size_min: int = 380
    product_size_max: int = 480

    max_ns_accepted: int = 0

    thermodynamic_oligo_alignment: bool = True
    thermodynamic_template_alignment: bool = False


@dataclass
class PrimerPair:
    rank: int
    primer3_index: int

    forward_sequence: str
    reverse_sequence: str

    forward_start: int
    forward_length: int
    reverse_start: int
    reverse_length: int

    forward_tm: float
    reverse_tm: float
    forward_gc: float
    reverse_gc: float

    product_size: int

    pair_penalty: float
    forward_penalty: float
    reverse_penalty: float

    forward_self_any: float
    forward_self_end: float
    forward_hairpin: float

    reverse_self_any: float
    reverse_self_end: float
    reverse_hairpin: float

    pair_any: float
    pair_end: float

    quality_score: float
    notes: str

    @property
    def suggested_annealing_temperature(self) -> float:
        return min(self.forward_tm, self.reverse_tm) - 2.0

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["suggested_annealing_temperature"] = (
            self.suggested_annealing_temperature
        )
        return result


def clean_sequence(sequence: str) -> str:
    cleaned = "".join(sequence.split()).upper()

    invalid = sorted(set(cleaned) - set("ACGTN"))

    if invalid:
        raise ValueError(
            "Sequence contains unsupported character(s): "
            + ", ".join(invalid)
        )

    if not cleaned:
        raise ValueError("Sequence is empty.")

    return cleaned


def _number(
    result: dict[str, Any],
    key: str,
    default: float = 0.0,
) -> float:
    value = result.get(key, default)

    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _location(
    result: dict[str, Any],
    key: str,
) -> tuple[int, int]:
    value = result.get(key)

    if isinstance(value, (list, tuple)) and len(value) >= 2:
        return int(value[0]), int(value[1])

    return 0, 0


def score_primer_pair(
    *,
    pair_penalty: float,
    forward_tm: float,
    reverse_tm: float,
    forward_gc: float,
    reverse_gc: float,
    forward_self_any: float,
    forward_self_end: float,
    forward_hairpin: float,
    reverse_self_any: float,
    reverse_self_end: float,
    reverse_hairpin: float,
    pair_any: float,
    pair_end: float,
) -> tuple[float, str]:
    score = 100.0
    notes: list[str] = []

    tm_difference = abs(forward_tm - reverse_tm)

    if tm_difference > 2.0:
        score -= min(20.0, tm_difference * 5.0)
        notes.append(f"Tm difference {tm_difference:.1f}°C")
    else:
        score -= tm_difference * 2.0

    score -= min(25.0, pair_penalty * 4.0)

    for label, gc_value in (
        ("forward GC", forward_gc),
        ("reverse GC", reverse_gc),
    ):
        if gc_value < 35.0 or gc_value > 65.0:
            score -= 18.0
            notes.append(f"{label} {gc_value:.1f}%")
        elif gc_value < 40.0 or gc_value > 60.0:
            score -= 8.0
            notes.append(f"{label} {gc_value:.1f}%")

    for label, value in (
        ("forward hairpin", forward_hairpin),
        ("reverse hairpin", reverse_hairpin),
    ):
        if value >= 50.0:
            score -= 25.0
            notes.append(f"{label} {value:.1f}°C")
        elif value >= 45.0:
            score -= 14.0
            notes.append(f"{label} {value:.1f}°C")
        elif value >= 40.0:
            score -= 7.0
            notes.append(f"{label} {value:.1f}°C")

    for label, value in (
        ("forward self-end", forward_self_end),
        ("reverse self-end", reverse_self_end),
        ("pair end", pair_end),
    ):
        if value >= 20.0:
            score -= 20.0
            notes.append(f"{label} {value:.1f}°C")
        elif value >= 10.0:
            score -= 10.0
            notes.append(f"{label} {value:.1f}°C")

    for label, value in (
        ("forward self-any", forward_self_any),
        ("reverse self-any", reverse_self_any),
        ("pair any", pair_any),
    ):
        if value >= 40.0:
            score -= 10.0
            notes.append(f"{label} {value:.1f}°C")
        elif value >= 30.0:
            score -= 5.0
            notes.append(f"{label} {value:.1f}°C")

    score = max(0.0, min(100.0, score))
    note_text = "; ".join(notes) if notes else "Clean Primer3 profile"

    return score, note_text


def _primer3_arguments(
    settings: Primer3Settings,
) -> dict[str, Any]:
    return {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": 0,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_NUM_RETURN": settings.number_to_return,
        "PRIMER_MIN_SIZE": settings.primer_min_size,
        "PRIMER_OPT_SIZE": settings.primer_opt_size,
        "PRIMER_MAX_SIZE": settings.primer_max_size,
        "PRIMER_MIN_TM": settings.primer_min_tm,
        "PRIMER_OPT_TM": settings.primer_opt_tm,
        "PRIMER_MAX_TM": settings.primer_max_tm,
        "PRIMER_MIN_GC": settings.primer_min_gc,
        "PRIMER_MAX_GC": settings.primer_max_gc,
        "PRIMER_MAX_NS_ACCEPTED": settings.max_ns_accepted,
        "PRIMER_PRODUCT_SIZE_RANGE": [
            [settings.product_size_min, settings.product_size_max]
        ],
        "PRIMER_EXPLAIN_FLAG": 1,
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": int(
            settings.thermodynamic_oligo_alignment
        ),
        "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT": int(
            settings.thermodynamic_template_alignment
        ),
    }


def design_primers(
    sequence: str,
    sequence_id: str = "candidate",
    settings: Primer3Settings | None = None,
) -> list[PrimerPair]:
    settings = settings or Primer3Settings()
    template = clean_sequence(sequence)

    if len(template) < settings.product_size_min:
        raise ValueError(
            f"Sequence length is {len(template)} bp, but the minimum "
            f"requested product size is {settings.product_size_min} bp."
        )

    sequence_args = {
        "SEQUENCE_ID": sequence_id,
        "SEQUENCE_TEMPLATE": template,
    }

    try:
        result = primer3.bindings.design_primers(
            sequence_args,
            _primer3_arguments(settings),
        )
    except Exception as exc:
        raise RuntimeError(
            f"Primer3 failed for {sequence_id}: {exc}"
        ) from exc

    returned = int(result.get("PRIMER_PAIR_NUM_RETURNED", 0))
    pairs: list[PrimerPair] = []

    for index in range(returned):
        forward_sequence = str(
            result.get(f"PRIMER_LEFT_{index}_SEQUENCE", "")
        )
        reverse_sequence = str(
            result.get(f"PRIMER_RIGHT_{index}_SEQUENCE", "")
        )

        if not forward_sequence or not reverse_sequence:
            continue

        forward_start, forward_length = _location(
            result,
            f"PRIMER_LEFT_{index}",
        )
        reverse_start, reverse_length = _location(
            result,
            f"PRIMER_RIGHT_{index}",
        )

        forward_tm = _number(result, f"PRIMER_LEFT_{index}_TM")
        reverse_tm = _number(result, f"PRIMER_RIGHT_{index}_TM")
        forward_gc = _number(
            result,
            f"PRIMER_LEFT_{index}_GC_PERCENT",
        )
        reverse_gc = _number(
            result,
            f"PRIMER_RIGHT_{index}_GC_PERCENT",
        )

        product_size = int(
            _number(result, f"PRIMER_PAIR_{index}_PRODUCT_SIZE")
        )

        pair_penalty = _number(
            result,
            f"PRIMER_PAIR_{index}_PENALTY",
        )
        forward_penalty = _number(
            result,
            f"PRIMER_LEFT_{index}_PENALTY",
        )
        reverse_penalty = _number(
            result,
            f"PRIMER_RIGHT_{index}_PENALTY",
        )

        forward_self_any = _number(
            result,
            f"PRIMER_LEFT_{index}_SELF_ANY_TH",
            _number(result, f"PRIMER_LEFT_{index}_SELF_ANY"),
        )
        forward_self_end = _number(
            result,
            f"PRIMER_LEFT_{index}_SELF_END_TH",
            _number(result, f"PRIMER_LEFT_{index}_SELF_END"),
        )
        forward_hairpin = _number(
            result,
            f"PRIMER_LEFT_{index}_HAIRPIN_TH",
        )

        reverse_self_any = _number(
            result,
            f"PRIMER_RIGHT_{index}_SELF_ANY_TH",
            _number(result, f"PRIMER_RIGHT_{index}_SELF_ANY"),
        )
        reverse_self_end = _number(
            result,
            f"PRIMER_RIGHT_{index}_SELF_END_TH",
            _number(result, f"PRIMER_RIGHT_{index}_SELF_END"),
        )
        reverse_hairpin = _number(
            result,
            f"PRIMER_RIGHT_{index}_HAIRPIN_TH",
        )

        pair_any = _number(
            result,
            f"PRIMER_PAIR_{index}_COMPL_ANY_TH",
            _number(result, f"PRIMER_PAIR_{index}_COMPL_ANY"),
        )
        pair_end = _number(
            result,
            f"PRIMER_PAIR_{index}_COMPL_END_TH",
            _number(result, f"PRIMER_PAIR_{index}_COMPL_END"),
        )

        quality_score, notes = score_primer_pair(
            pair_penalty=pair_penalty,
            forward_tm=forward_tm,
            reverse_tm=reverse_tm,
            forward_gc=forward_gc,
            reverse_gc=reverse_gc,
            forward_self_any=forward_self_any,
            forward_self_end=forward_self_end,
            forward_hairpin=forward_hairpin,
            reverse_self_any=reverse_self_any,
            reverse_self_end=reverse_self_end,
            reverse_hairpin=reverse_hairpin,
            pair_any=pair_any,
            pair_end=pair_end,
        )

        pairs.append(
            PrimerPair(
                rank=0,
                primer3_index=index,
                forward_sequence=forward_sequence,
                reverse_sequence=reverse_sequence,
                forward_start=forward_start,
                forward_length=forward_length,
                reverse_start=reverse_start,
                reverse_length=reverse_length,
                forward_tm=forward_tm,
                reverse_tm=reverse_tm,
                forward_gc=forward_gc,
                reverse_gc=reverse_gc,
                product_size=product_size,
                pair_penalty=pair_penalty,
                forward_penalty=forward_penalty,
                reverse_penalty=reverse_penalty,
                forward_self_any=forward_self_any,
                forward_self_end=forward_self_end,
                forward_hairpin=forward_hairpin,
                reverse_self_any=reverse_self_any,
                reverse_self_end=reverse_self_end,
                reverse_hairpin=reverse_hairpin,
                pair_any=pair_any,
                pair_end=pair_end,
                quality_score=quality_score,
                notes=notes,
            )
        )

    pairs.sort(
        key=lambda pair: (
            -pair.quality_score,
            pair.pair_penalty,
            abs(pair.forward_tm - pair.reverse_tm),
            pair.primer3_index,
        )
    )

    for rank, pair in enumerate(pairs, start=1):
        pair.rank = rank

    return pairs


def design_primers_for_candidates(
    candidates: Iterable[Any],
    *,
    settings: Primer3Settings | None = None,
    sequence_attribute: str = "sequence",
    id_attribute: str = "rank",
) -> list[tuple[Any, list[PrimerPair]]]:
    settings = settings or Primer3Settings()
    completed: list[tuple[Any, list[PrimerPair]]] = []

    for candidate in candidates:
        sequence = getattr(candidate, sequence_attribute)
        candidate_id = getattr(
            candidate,
            id_attribute,
            len(completed) + 1,
        )

        pairs = design_primers(
            sequence=sequence,
            sequence_id=f"candidate_{candidate_id}",
            settings=settings,
        )

        completed.append((candidate, pairs))

    return completed


def print_primer_pairs(pairs: list[PrimerPair]) -> None:
    if not pairs:
        print("Primer3 returned no acceptable primer pairs.")
        return

    print(
        f"{'Rank':<6}"
        f"{'Forward primer':<28}"
        f"{'Reverse primer':<28}"
        f"{'Product':>9}"
        f"{'Tm F/R':>15}"
        f"{'Score':>9}"
    )
    print("-" * 95)

    for pair in pairs:
        tm_text = f"{pair.forward_tm:.1f}/{pair.reverse_tm:.1f}"

        print(
            f"{pair.rank:<6}"
            f"{pair.forward_sequence:<28}"
            f"{pair.reverse_sequence:<28}"
            f"{pair.product_size:>9}"
            f"{tm_text:>15}"
            f"{pair.quality_score:>9.1f}"
        )