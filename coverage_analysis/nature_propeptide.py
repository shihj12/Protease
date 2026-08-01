"""Match precursor/mature diagnostic candidates to Nature-observed peptides."""

from __future__ import annotations

import csv
import gzip
from dataclasses import asdict
from pathlib import Path

from coverage_analysis.propeptide_quantification import (
    PeptideCandidate,
    _best,
    junction_candidates,
    processing_junctions,
)
from coverage_analysis.protease_coverage import NATURE_ENZYMES


def _observed_bits_for_candidates(
    peptide_cache: Path, candidate_sequences: set[str]
) -> dict[str, int]:
    observed: dict[str, int] = {}
    with gzip.open(peptide_cache, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            sequence = row["sequence"]
            if sequence in candidate_sequences:
                observed[sequence] = observed.get(sequence, 0) | int(row["enzyme_bits"])
    return observed


def _is_observed(
    peptide: PeptideCandidate, observed_bits: dict[str, int], enzyme_bit: dict[str, int]
) -> bool:
    return bool(observed_bits.get(peptide.sequence, 0) & enzyme_bit[peptide.enzyme])


def quantify_nature_propeptides(entries, designs, peptide_cache: Path):
    """Apply the theoretical candidate definitions to PXD024364 observations.

    The Nature search required specific cleavage. Consequently an observed
    semi-specific mature neo-terminal peptide is not expected; this function
    still tests exact sequences so that the limitation is measured, not assumed.
    """

    junctions = [
        junction
        for entry in entries
        for junction in processing_junctions(entry)
    ]
    entries_by_accession = {entry.accession: entry for entry in entries}
    enzymes = sorted({enzyme for _, _, group in designs for enzyme in group})
    enzyme_bit = {enzyme: 1 << i for i, enzyme in enumerate(NATURE_ENZYMES)}

    inventory: dict[
        tuple[tuple[str, int, str], str],
        tuple[list[PeptideCandidate], list[PeptideCandidate]],
    ] = {}
    candidate_sequences: set[str] = set()
    for junction in junctions:
        entry = entries_by_accession[junction.accession]
        for enzyme in enzymes:
            intact, mature = junction_candidates(entry, junction, enzyme)
            inventory[(junction.key, enzyme)] = intact, mature
            candidate_sequences.update(peptide.sequence for peptide in intact)
            candidate_sequences.update(peptide.sequence for peptide in mature)

    observed_bits = _observed_bits_for_candidates(peptide_cache, candidate_sequences)
    detail_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for family, combination, group in designs:
        design_rows: list[dict[str, object]] = []
        for junction in junctions:
            intact_candidates: list[PeptideCandidate] = []
            mature_candidates: list[PeptideCandidate] = []
            for enzyme in dict.fromkeys(group):
                intact, mature = inventory[(junction.key, enzyme)]
                intact_candidates.extend(
                    peptide
                    for peptide in intact
                    if _is_observed(peptide, observed_bits, enzyme_bit)
                )
                mature_candidates.extend(
                    peptide
                    for peptide in mature
                    if _is_observed(peptide, observed_bits, enzyme_bit)
                )

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
            intact_best = _best(intact_candidates)
            mature_best = _best(mature_candidates)
            row: dict[str, object] = {
                "family": family,
                "combination": combination,
                "enzymes": ";".join(group),
                **asdict(junction),
                "bond": junction.bond,
                "intact_observed": bool(intact_candidates),
                "mature_observed": bool(mature_candidates),
                "both_observed": bool(intact_candidates and mature_candidates),
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
                "reference_search_specificity": "specific cleavage",
            }
            design_rows.append(row)
            detail_rows.append(row)

        counts = {
            outcome: sum(row["outcome"] == outcome for row in design_rows)
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
                    bool(row["both_unique_zero_missed"]) for row in design_rows
                ),
                "intact_only": counts["intact only"],
                "mature_only": counts["mature only"],
                "neither": counts["neither"],
                "intact_observed_unique": sum(
                    bool(row["intact_unique"]) for row in design_rows
                ),
                "mature_observed_unique": sum(
                    bool(row["mature_unique"]) for row in design_rows
                ),
                "both_detectable": sum(
                    bool(row["both_observed"]) for row in design_rows
                ),
                "proteins_with_intact": len(
                    {row["accession"] for row in design_rows if row["intact_unique"]}
                ),
                "total_proteins": len(entries),
                "search_specificity": "specific cleavage",
            }
        )
    return junctions, detail_rows, summary_rows, len(observed_bits)
