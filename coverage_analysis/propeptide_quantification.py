"""Paired quantification of intact precursor and processed mature forms.

An activation junction is jointly quantifiable when parallel digestion yields:

* a protein-unique, fully-specific peptide spanning the intact junction; and
* a protein-unique, semi-specific peptide beginning or ending at the mature
  product boundary.

If a digestion enzyme cuts at the biological processing bond, a missed-cleavage
spanning peptide can still report intact precursor, but its ordinary terminal
peptide is not mature-state-specific.
"""

from __future__ import annotations

import sys
from dataclasses import asdict, dataclass
from itertools import chain
from pathlib import Path
from typing import Iterable

from coverage_analysis.protease_coverage import (
    MAX_PEPTIDE_LENGTH,
    MIN_PEPTIDE_LENGTH,
    MISSED_CLEAVAGES,
    ProteaseRule,
    RULES,
    cleavage_boundaries,
)


ROOT = Path(__file__).resolve().parents[1]
ASE_ASSAY_DIR = ROOT / "ase_assay"
if str(ASE_ASSAY_DIR) not in sys.path:
    sys.path.insert(0, str(ASE_ASSAY_DIR))

import proteome  # noqa: E402  (legacy module is intentionally loaded by path)
from presets import PRESETS  # noqa: E402
from uniprot import fetch_many  # noqa: E402


@dataclass(frozen=True)
class ProcessingJunction:
    accession: str
    gene: str
    protein: str
    position: int
    side: str
    product_start: int
    product_end: int
    product: str
    junction_type: str
    evidence: str

    @property
    def key(self) -> tuple[str, int, str]:
        return self.accession, self.position, self.side

    @property
    def bond(self) -> str:
        return f"{self.position}|{self.position + 1}"


@dataclass(frozen=True)
class PeptideCandidate:
    enzyme: str
    start: int
    end: int
    sequence: str
    missed_cleavages: int
    unique: bool

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def preset_symbols() -> list[str]:
    """All preset gene symbols, once each and in preset order."""

    return list(dict.fromkeys(chain.from_iterable(PRESETS.values())))


def load_preset_entries():
    entries, errors = fetch_many(preset_symbols())
    if errors:
        detail = "; ".join(f"{symbol}: {message}" for symbol, message in errors)
        raise RuntimeError(f"Unable to resolve all preset proteins: {detail}")
    return entries


def processing_junctions(entry) -> list[ProcessingJunction]:
    """Return terminal pro-peptide junctions plus PGRN granulin boundaries.

    UniProt represents conventional zymogen activation with adjacent ``PROPEP``
    and ``CHAIN`` features. PGRN is different: its processed products are eight
    ``PEPTIDE`` features, so internal granulin boundaries are added explicitly.
    The signal-peptide/paragranulin boundary is excluded because it is secretion
    processing rather than granulin release.
    """

    by_key: dict[tuple[int, str], ProcessingJunction] = {}

    def consider(candidate: ProcessingJunction) -> None:
        key = (candidate.position, candidate.side)
        current = by_key.get(key)
        if current is None:
            by_key[key] = candidate
            return
        candidate_length = candidate.product_end - candidate.product_start
        current_length = current.product_end - current.product_start
        if candidate_length > current_length:
            by_key[key] = candidate

    propeptides = entry.of_kind("PROPEP")
    for mature in entry.of_kind("CHAIN"):
        if mature.start > 1:
            position = mature.start - 1
            adjacent = next((f for f in propeptides if f.end == position), None)
            if adjacent is not None:
                consider(
                    ProcessingJunction(
                        entry.accession,
                        entry.gene,
                        entry.protein,
                        position,
                        "N",
                        mature.start,
                        mature.end,
                        mature.description or entry.protein,
                        "propeptide",
                        adjacent.evidence,
                    )
                )
        if mature.end < entry.length:
            position = mature.end
            adjacent = next((f for f in propeptides if f.start == position + 1), None)
            if adjacent is not None:
                consider(
                    ProcessingJunction(
                        entry.accession,
                        entry.gene,
                        entry.protein,
                        position,
                        "C",
                        mature.start,
                        mature.end,
                        mature.description or entry.protein,
                        "propeptide",
                        adjacent.evidence,
                    )
                )

    if entry.gene.upper() == "GRN" or entry.accession == "P28799":
        signal_ends = {feature.end for feature in entry.of_kind("SIGNAL")}
        for product in entry.of_kind("PEPTIDE"):
            if "granulin" not in product.description.lower():
                continue
            if product.start > 1 and product.start - 1 not in signal_ends:
                consider(
                    ProcessingJunction(
                        entry.accession,
                        entry.gene,
                        entry.protein,
                        product.start - 1,
                        "N",
                        product.start,
                        product.end,
                        product.description,
                        "granulin",
                        product.evidence,
                    )
                )
            if product.end < entry.length:
                consider(
                    ProcessingJunction(
                        entry.accession,
                        entry.gene,
                        entry.protein,
                        product.end,
                        "C",
                        product.start,
                        product.end,
                        product.description,
                        "granulin",
                        product.evidence,
                    )
                )

    return [by_key[key] for key in sorted(by_key)]


def _fully_specific_intervals(
    sequence: str,
    rule: ProteaseRule,
    min_length: int,
    max_length: int,
    missed_cleavages: int,
) -> list[tuple[int, int, int]]:
    boundaries = cleavage_boundaries(sequence, rule)
    intervals: list[tuple[int, int, int]] = []
    for start_fragment in range(len(boundaries) - 1):
        for missed in range(missed_cleavages + 1):
            end_fragment = start_fragment + missed
            if end_fragment >= len(boundaries) - 1:
                break
            start = boundaries[start_fragment]
            end = boundaries[end_fragment + 1]
            if min_length <= end - start <= max_length:
                intervals.append((start, end, missed))
    return intervals


def _candidate(
    sequence: str, enzyme: str, start: int, end: int, missed: int
) -> PeptideCandidate:
    peptide = sequence[start:end]
    return PeptideCandidate(
        enzyme=enzyme,
        start=start + 1,
        end=end,
        sequence=peptide,
        missed_cleavages=missed,
        unique=bool(proteome.is_unique(peptide)),
    )


def junction_candidates(
    entry,
    junction: ProcessingJunction,
    enzyme: str,
    min_length: int = MIN_PEPTIDE_LENGTH,
    max_length: int = MAX_PEPTIDE_LENGTH,
    missed_cleavages: int = MISSED_CLEAVAGES,
) -> tuple[list[PeptideCandidate], list[PeptideCandidate]]:
    """Candidate intact-spanning and mature-terminal peptides for one aliquot."""

    sequence = entry.sequence.upper()
    rule = RULES[enzyme]
    boundaries = cleavage_boundaries(sequence, rule)
    position = junction.position  # also the zero-based boundary index

    cuts_at_junction = position in boundaries[1:-1]

    intact = [
        _candidate(sequence, enzyme, start, end, missed)
        for start, end, missed in _fully_specific_intervals(
            sequence, rule, min_length, max_length, missed_cleavages
        )
        if start < position < end
    ]

    # At an exact digestive cut, a spanning peptide necessarily contains a
    # missed cleavage and remains specific to physically intact precursor. The
    # corresponding ordinary terminal peptide, however, also arises when the
    # precursor was never biologically processed and is not mature-specific.
    if cuts_at_junction:
        return intact, []

    product_start = junction.product_start - 1
    product_end = junction.product_end
    internal = [
        boundary
        for boundary in boundaries[1:-1]
        if product_start < boundary < product_end
    ]
    mature: list[PeptideCandidate] = []
    if junction.side == "N":
        endpoints = [boundary for boundary in internal if boundary > position]
        endpoints.append(product_end)
        for missed, end in enumerate(endpoints[: missed_cleavages + 1]):
            if min_length <= end - position <= max_length:
                mature.append(_candidate(sequence, enzyme, position, end, missed))
    else:
        starts = [product_start]
        starts.extend(boundary for boundary in internal if boundary < position)
        starts = sorted(set(starts), reverse=True)
        for missed, start in enumerate(starts[: missed_cleavages + 1]):
            if min_length <= position - start <= max_length:
                mature.append(_candidate(sequence, enzyme, start, position, missed))
    return intact, mature


def _best(candidates: Iterable[PeptideCandidate]) -> PeptideCandidate | None:
    candidates = list(candidates)
    if not candidates:
        return None
    return min(
        candidates,
        key=lambda peptide: (
            not peptide.unique,
            peptide.missed_cleavages,
            abs(peptide.length - 15),
            peptide.length,
            peptide.enzyme,
        ),
    )


def quantify_designs(entries, designs):
    """Return per-junction results and aggregate rows for parallel designs."""

    entries_by_accession = {entry.accession: entry for entry in entries}
    junctions = [
        junction
        for entry in entries
        for junction in processing_junctions(entry)
    ]
    enzymes = sorted({enzyme for _, _, group in designs for enzyme in group})
    candidates: dict[
        tuple[tuple[str, int, str], str],
        tuple[list[PeptideCandidate], list[PeptideCandidate]],
    ] = {}
    for junction in junctions:
        entry = entries_by_accession[junction.accession]
        for enzyme in enzymes:
            candidates[(junction.key, enzyme)] = junction_candidates(
                entry, junction, enzyme
            )

    detail_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for family, combination, group in designs:
        family_rows: list[dict[str, object]] = []
        for junction in junctions:
            intact_candidates: list[PeptideCandidate] = []
            mature_candidates: list[PeptideCandidate] = []
            for enzyme in dict.fromkeys(group):
                intact, mature = candidates[(junction.key, enzyme)]
                intact_candidates.extend(intact)
                mature_candidates.extend(mature)

            intact_best = _best(intact_candidates)
            mature_best = _best(mature_candidates)
            intact_detectable = bool(intact_candidates)
            mature_detectable = bool(mature_candidates)
            intact_unique = any(peptide.unique for peptide in intact_candidates)
            mature_unique = any(peptide.unique for peptide in mature_candidates)
            intact_unique_zero = any(
                peptide.unique and peptide.missed_cleavages == 0
                for peptide in intact_candidates
            )
            mature_unique_zero = any(
                peptide.unique and peptide.missed_cleavages == 0
                for peptide in mature_candidates
            )
            if intact_unique and mature_unique:
                outcome = "both"
            elif intact_unique:
                outcome = "intact only"
            elif mature_unique:
                outcome = "mature only"
            else:
                outcome = "neither"

            row: dict[str, object] = {
                "family": family,
                "combination": combination,
                "enzymes": ";".join(group),
                **asdict(junction),
                "bond": junction.bond,
                "intact_detectable": intact_detectable,
                "mature_detectable": mature_detectable,
                "both_detectable": intact_detectable and mature_detectable,
                "intact_unique": intact_unique,
                "mature_unique": mature_unique,
                "both_unique": intact_unique and mature_unique,
                "both_unique_zero_missed": (
                    intact_unique_zero and mature_unique_zero
                ),
                "outcome": outcome,
                "intact_peptide": intact_best.sequence if intact_best else "",
                "intact_enzyme": intact_best.enzyme if intact_best else "",
                "intact_start": intact_best.start if intact_best else "",
                "intact_end": intact_best.end if intact_best else "",
                "intact_missed": (
                    intact_best.missed_cleavages if intact_best else ""
                ),
                "mature_peptide": mature_best.sequence if mature_best else "",
                "mature_enzyme": mature_best.enzyme if mature_best else "",
                "mature_start": mature_best.start if mature_best else "",
                "mature_end": mature_best.end if mature_best else "",
                "mature_missed": (
                    mature_best.missed_cleavages if mature_best else ""
                ),
            }
            family_rows.append(row)
            detail_rows.append(row)

        counts = {
            outcome: sum(row["outcome"] == outcome for row in family_rows)
            for outcome in ("both", "intact only", "mature only", "neither")
        }
        summary_rows.append(
            {
                "family": family,
                "combination": combination,
                "enzymes": ";".join(group),
                "total_junctions": len(junctions),
                "both_unique": counts["both"],
                "both_unique_zero_missed": sum(
                    bool(row["both_unique_zero_missed"]) for row in family_rows
                ),
                "intact_only": counts["intact only"],
                "mature_only": counts["mature only"],
                "neither": counts["neither"],
                "both_detectable": sum(
                    bool(row["both_detectable"]) for row in family_rows
                ),
                "proteins_with_both": len(
                    {row["accession"] for row in family_rows if row["both_unique"]}
                ),
                "total_proteins": len(entries),
                "pct_both_unique": (
                    100.0 * counts["both"] / len(junctions) if junctions else 0.0
                ),
            }
        )

    return junctions, detail_rows, summary_rows
