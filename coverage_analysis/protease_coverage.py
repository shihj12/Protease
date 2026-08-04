"""Pure in-silico digestion and sequence-coverage calculations.

The implementation deliberately does not depend on the legacy Streamlit app in
``ase_assay``. A peptide is detectable when its length is within the configured
window. Parallel digests are combined by taking the union of residue masks;
they are never treated as a single-tube co-digest.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from statistics import median
from typing import Callable, Iterable, Iterator


MIN_PEPTIDE_LENGTH = 7
MAX_PEPTIDE_LENGTH = 52
MISSED_CLEAVAGES = 2

LABELS: dict[str, frozenset[str]] = {
    "K+R": frozenset("KR"),
    "L": frozenset("L"),
    "T": frozenset("T"),
    "I": frozenset("I"),
    "V": frozenset("V"),
    "H": frozenset("H"),
}

KR_PLUS_ONE_LABELS: dict[str, frozenset[str]] = {
    f"K+R+{residue}": frozenset(f"KR{residue}")
    for residue in ("L", "T", "I", "V", "H")
}

# Arginine is the SILAC label that misbehaves metabolically (Arg-to-Pro
# conversion), so these sets ask what a lysine-anchored scheme recovers with
# arginine dropped entirely. Every set is a union at the peptide level: a
# peptide counts when it carries at least one of the listed residues.
LYSINE_PARTNER_RESIDUES = ("L", "T", "I", "V", "H")

LYSINE_LABELS: dict[str, frozenset[str]] = {
    "K": frozenset("K"),
    **{
        f"K+{residue}": frozenset(f"K{residue}")
        for residue in LYSINE_PARTNER_RESIDUES
    },
    **{
        f"K+{first}+{second}": frozenset(f"K{first}{second}")
        for first, second in combinations(LYSINE_PARTNER_RESIDUES, 2)
    },
}

# Carried through the same pass so every arginine-free figure can show what was
# given up: the incumbent K+R and the best K+R-plus-one scheme.
ARGININE_REFERENCE_LABELS: dict[str, frozenset[str]] = {
    "R": frozenset("R"),
    "K+R": frozenset("KR"),
    "K+R+L": frozenset("KRL"),
}

# One pass covers the arginine-free sets, the single residues, and the
# arginine references, so the cached label matrix answers all of them.
LYSINE_STUDY_LABELS: dict[str, frozenset[str]] = {
    **LABELS,
    **ARGININE_REFERENCE_LABELS,
    **LYSINE_LABELS,
}

NATURE_ENZYMES = (
    "Trypsin",
    "Lys-C",
    "Chymotrypsin",
    "Lys-N",
    "Glu-C",
    "Asp-N",
)

REQUESTED_GROUPS = (
    "Trypsin",
    "Trypsin/Lys-C",
    "Chymotrypsin",
    "Lys-N",
    "Glu-C",
    "Asp-N",
    "Arg-N",
    "ProAlanase",
    "Elastase",
    "Pepsin",
)

ALTERNATIVE_GROUPS = tuple(
    name for name in REQUESTED_GROUPS if name not in {"Trypsin", "Trypsin/Lys-C"}
)


@dataclass(frozen=True)
class Protein:
    identifier: str
    sequence: str


@dataclass(frozen=True)
class ProteaseRule:
    name: str
    cut: Callable[[str, int], bool]
    description: str


def _after(residues: str, blocked_next: str = "") -> Callable[[str, int], bool]:
    residues_set = frozenset(residues)
    blocked_set = frozenset(blocked_next)

    def predicate(sequence: str, left_index: int) -> bool:
        return (
            sequence[left_index] in residues_set
            and sequence[left_index + 1] not in blocked_set
        )

    return predicate


def _before(residues: str) -> Callable[[str, int], bool]:
    residues_set = frozenset(residues)

    def predicate(sequence: str, left_index: int) -> bool:
        return sequence[left_index + 1] in residues_set

    return predicate


def _pepsin_ph13(sequence: str, left_index: int) -> bool:
    """PeptideCutter's deterministic pH 1.3 approximation.

    Pepsin is broad and condition-dependent, so this is explicitly a theoretical
    preference model, not a claim of strict specificity.
    """

    p1 = sequence[left_index]
    p1_prime = sequence[left_index + 1]
    p2 = sequence[left_index - 1] if left_index >= 1 else None
    p3 = sequence[left_index - 2] if left_index >= 2 else None
    p2_prime = sequence[left_index + 2] if left_index + 2 < len(sequence) else None
    common_allowed = p3 not in frozenset("HKR") and p2 != "P" and p2_prime != "P"
    if not common_allowed:
        return False
    return p1 in frozenset("FL") or (p1_prime in frozenset("FL") and p1 != "R")


RULES: dict[str, ProteaseRule] = {
    "Trypsin": ProteaseRule(
        "Trypsin", _after("KR", "P"), "C-terminal to K/R, except before P."
    ),
    "Lys-C": ProteaseRule("Lys-C", _after("K"), "C-terminal to K."),
    "Chymotrypsin": ProteaseRule(
        "Chymotrypsin",
        _after("FYW", "P"),
        "High-specificity rule: C-terminal to F/Y/W, except before P.",
    ),
    "Lys-N": ProteaseRule("Lys-N", _before("K"), "N-terminal to K."),
    "Glu-C": ProteaseRule("Glu-C", _after("E"), "C-terminal to E."),
    "Asp-N": ProteaseRule("Asp-N", _before("D"), "N-terminal to D."),
    "Arg-N": ProteaseRule(
        "Arg-N", _before("R"), "Theoretical N-terminal arginine-specific rule."
    ),
    "ProAlanase": ProteaseRule(
        "ProAlanase",
        _after("PA"),
        "C-terminal to P and A; preferential rather than perfectly specific.",
    ),
    "Elastase": ProteaseRule(
        "Elastase",
        _after("ALIV", "P"),
        "Proteome Discoverer rule: C-terminal to A/L/I/V, except before P.",
    ),
    "Pepsin": ProteaseRule(
        "Pepsin", _pepsin_ph13, "PeptideCutter pH 1.3 preference model."
    ),
}


def read_fasta(path: Path) -> list[Protein]:
    proteins: list[Protein] = []
    identifier: str | None = None
    sequence_parts: list[str] = []

    def flush() -> None:
        if identifier is not None:
            sequence = "".join(sequence_parts).upper()
            if sequence:
                proteins.append(Protein(identifier, sequence))

    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:]
                pipe_parts = header.split("|")
                identifier = pipe_parts[1] if len(pipe_parts) >= 3 else header.split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line)
    flush()
    return proteins


def cleavage_boundaries(sequence: str, rule: ProteaseRule) -> list[int]:
    """Return zero-based peptide boundaries, including both protein termini."""

    return [0] + [
        index + 1
        for index in range(len(sequence) - 1)
        if rule.cut(sequence, index)
    ] + [len(sequence)]


def peptide_intervals(
    sequence: str,
    rule: ProteaseRule,
    min_length: int = MIN_PEPTIDE_LENGTH,
    max_length: int = MAX_PEPTIDE_LENGTH,
    missed_cleavages: int = MISSED_CLEAVAGES,
) -> Iterator[tuple[int, int]]:
    """Yield half-open intervals for all peptides within the search window."""

    boundaries = cleavage_boundaries(sequence, rule)
    fragment_count = len(boundaries) - 1
    for start_fragment in range(fragment_count):
        for missed in range(missed_cleavages + 1):
            end_fragment = start_fragment + missed
            if end_fragment >= fragment_count:
                break
            start = boundaries[start_fragment]
            end = boundaries[end_fragment + 1]
            if min_length <= end - start <= max_length:
                yield start, end


def interval_mask(start: int, end: int) -> int:
    return ((1 << (end - start)) - 1) << start


def digest_mask(
    sequence: str,
    rule: ProteaseRule,
    labels: dict[str, frozenset[str]] | None = None,
    min_length: int = MIN_PEPTIDE_LENGTH,
    max_length: int = MAX_PEPTIDE_LENGTH,
    missed_cleavages: int = MISSED_CLEAVAGES,
) -> tuple[int, dict[str, int], int]:
    covered = 0
    labelled = {name: 0 for name in labels or {}}
    peptide_count = 0
    for start, end in peptide_intervals(
        sequence, rule, min_length, max_length, missed_cleavages
    ):
        mask = interval_mask(start, end)
        covered |= mask
        peptide_count += 1
        if labels:
            peptide_residues = frozenset(sequence[start:end])
            for label_name, residues in labels.items():
                if not peptide_residues.isdisjoint(residues):
                    labelled[label_name] |= mask
    return covered, labelled, peptide_count


def digest_proteome(
    proteins: list[Protein],
    enzyme_names: Iterable[str],
    labels: dict[str, frozenset[str]] | None = None,
    min_length: int = MIN_PEPTIDE_LENGTH,
    max_length: int = MAX_PEPTIDE_LENGTH,
    missed_cleavages: int = MISSED_CLEAVAGES,
    progress: Callable[[str, int, int], None] | None = None,
) -> tuple[dict[str, list[int]], dict[str, dict[str, list[int]]], dict[str, int]]:
    enzyme_masks: dict[str, list[int]] = {}
    label_masks: dict[str, dict[str, list[int]]] = {}
    peptide_counts: dict[str, int] = {}
    total = len(proteins)
    for enzyme_name in enzyme_names:
        rule = RULES[enzyme_name]
        masks: list[int] = []
        per_label = {name: [] for name in labels or {}}
        peptide_count = 0
        for index, protein in enumerate(proteins, start=1):
            mask, labelled, count = digest_mask(
                protein.sequence,
                rule,
                labels,
                min_length,
                max_length,
                missed_cleavages,
            )
            masks.append(mask)
            peptide_count += count
            for label_name, label_mask in labelled.items():
                per_label[label_name].append(label_mask)
            if progress and (index % 2500 == 0 or index == total):
                progress(enzyme_name, index, total)
        enzyme_masks[enzyme_name] = masks
        label_masks[enzyme_name] = per_label
        peptide_counts[enzyme_name] = peptide_count
    return enzyme_masks, label_masks, peptide_counts


def group_to_enzymes(group: str) -> tuple[str, ...]:
    return ("Trypsin", "Lys-C") if group == "Trypsin/Lys-C" else (group,)


def display_combination(groups: Iterable[str]) -> str:
    return " + ".join(groups)


def requested_designs() -> tuple[
    list[tuple[str, tuple[str, ...]]],
    list[tuple[str, tuple[str, ...]]],
    list[tuple[str, tuple[str, ...]]],
]:
    singles = [(group, group_to_enzymes(group)) for group in REQUESTED_GROUPS]
    pairs = [
        (
            display_combination(("Trypsin/Lys-C", alternative)),
            ("Trypsin", "Lys-C", alternative),
        )
        for alternative in ALTERNATIVE_GROUPS
    ]
    triples = [
        (
            display_combination(("Trypsin/Lys-C", first, second)),
            ("Trypsin", "Lys-C", first, second),
        )
        for first, second in combinations(ALTERNATIVE_GROUPS, 2)
    ]
    return singles, pairs, triples


def nature_designs() -> tuple[
    list[tuple[str, tuple[str, ...]]],
    list[tuple[str, tuple[str, ...]]],
    list[tuple[str, tuple[str, ...]]],
]:
    alternatives = ("Chymotrypsin", "Lys-N", "Glu-C", "Asp-N")
    singles = [
        ("Trypsin", ("Trypsin",)),
        ("Trypsin/Lys-C", ("Trypsin", "Lys-C")),
        ("Chymotrypsin", ("Chymotrypsin",)),
        ("Lys-N", ("Lys-N",)),
        ("Glu-C", ("Glu-C",)),
        ("Asp-N", ("Asp-N",)),
    ]
    pairs = [
        (
            display_combination(("Trypsin/Lys-C", alternative)),
            ("Trypsin", "Lys-C", alternative),
        )
        for alternative in alternatives
    ]
    triples = [
        (
            display_combination(("Trypsin/Lys-C", first, second)),
            ("Trypsin", "Lys-C", first, second),
        )
        for first, second in combinations(alternatives, 2)
    ]
    return singles, pairs, triples


def coverage_statistics(
    proteins: list[Protein],
    enzyme_masks: dict[str, list[int]],
    enzymes: Iterable[str],
) -> dict[str, float | int]:
    enzyme_order = tuple(dict.fromkeys(enzymes))
    covered_residues = 0
    total_residues = 0
    protein_percentages: list[float] = []
    for index, protein in enumerate(proteins):
        mask = 0
        for enzyme in enzyme_order:
            mask |= enzyme_masks[enzyme][index]
        covered = mask.bit_count()
        length = len(protein.sequence)
        covered_residues += covered
        total_residues += length
        protein_percentages.append(100.0 * covered / length if length else 0.0)
    return {
        "n_proteins": len(proteins),
        "covered_residues": covered_residues,
        "total_residues": total_residues,
        "residue_weighted_pct": 100.0 * covered_residues / total_residues,
        "median_pct": median(protein_percentages),
    }


def rank_designs(
    proteins: list[Protein],
    enzyme_masks: dict[str, list[int]],
    designs: Iterable[tuple[str, tuple[str, ...]]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for combination, enzymes in designs:
        row: dict[str, object] = {
            "combination": combination,
            "enzymes": ";".join(enzymes),
        }
        row.update(coverage_statistics(proteins, enzyme_masks, enzymes))
        rows.append(row)
    rows.sort(
        key=lambda row: (
            float(row["residue_weighted_pct"]),
            float(row["median_pct"]),
            str(row["combination"]),
        ),
        reverse=True,
    )
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows


def labelled_design_statistics(
    proteins: list[Protein],
    label_masks: dict[str, dict[str, list[int]]],
    designs: Iterable[tuple[str, tuple[str, ...]]],
    labels: Iterable[str] = LABELS,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    total_residues = sum(len(protein.sequence) for protein in proteins)
    for combination, enzymes in designs:
        for label in labels:
            covered_total = 0
            protein_percentages: list[float] = []
            for index, protein in enumerate(proteins):
                mask = 0
                for enzyme in dict.fromkeys(enzymes):
                    mask |= label_masks[enzyme][label][index]
                covered = mask.bit_count()
                covered_total += covered
                protein_percentages.append(100.0 * covered / len(protein.sequence))
            rows.append(
                {
                    "combination": combination,
                    "enzymes": ";".join(enzymes),
                    "label": label,
                    "covered_residues": covered_total,
                    "total_residues": total_residues,
                    "residue_weighted_pct": 100.0 * covered_total / total_residues,
                    "median_pct": median(protein_percentages),
                }
            )
    return rows
