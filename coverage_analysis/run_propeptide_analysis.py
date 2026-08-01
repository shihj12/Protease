"""Generate paired precursor/mature-state quantification tables and figures."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from coverage_analysis.propeptide_plots import (
    combination_plot,
    outcome_plot,
    pgrn_heatmap,
    protein_heatmap,
)
from coverage_analysis.propeptide_quantification import (
    load_preset_entries,
    preset_symbols,
    processing_junctions,
    quantify_designs,
)
from coverage_analysis.protease_coverage import (
    MAX_PEPTIDE_LENGTH,
    MIN_PEPTIDE_LENGTH,
    MISSED_CLEAVAGES,
    REQUESTED_GROUPS,
    requested_designs,
)


DEFAULT_OUTPUT = Path(__file__).resolve().parent / "output"


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] = []
    for row in rows:
        columns.extend(key for key in row if key not in columns)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def _designs():
    singles, pairs, triples = requested_designs()
    return [
        *[("single", combination, enzymes) for combination, enzymes in singles],
        *[("two", combination, enzymes) for combination, enzymes in pairs],
        *[("three", combination, enzymes) for combination, enzymes in triples],
    ]


def run(output_dir: Path) -> None:
    entries = load_preset_entries()
    designs = _designs()
    junctions, detail_rows, summary_rows = quantify_designs(entries, designs)

    junction_rows = [
        {
            "gene": junction.gene,
            "accession": junction.accession,
            "protein": junction.protein,
            "bond": junction.bond,
            "side": junction.side,
            "product_start": junction.product_start,
            "product_end": junction.product_end,
            "product": junction.product,
            "junction_type": junction.junction_type,
            "evidence": junction.evidence,
        }
        for junction in junctions
    ]
    _write_csv(output_dir / "propeptide_junctions.csv", junction_rows)
    _write_csv(output_dir / "propeptide_quantifiability.csv", detail_rows)
    _write_csv(output_dir / "propeptide_design_summary.csv", summary_rows)

    singles = [row for row in summary_rows if row["family"] == "single"]
    single_order = {group: index for index, group in enumerate(REQUESTED_GROUPS)}
    singles.sort(key=lambda row: single_order[str(row["combination"])])
    pairs = [row for row in summary_rows if row["family"] == "two"]
    triples = [row for row in summary_rows if row["family"] == "three"]

    outcome_plot(singles, output_dir / "graph_8_propeptide_group_outcomes.png")
    combination_plot(pairs, output_dir / "graph_9_propeptide_two_groups.png")
    combination_plot(triples, output_dir / "graph_10_propeptide_three_groups.png")
    protein_heatmap(
        detail_rows,
        list(REQUESTED_GROUPS),
        output_dir / "graph_11_propeptide_protein_heatmap.png",
    )
    pgrn_heatmap(
        detail_rows,
        list(REQUESTED_GROUPS),
        output_dir / "graph_12_pgrn_junction_heatmap.png",
    )

    ranked_pairs = sorted(
        pairs,
        key=lambda row: (
            int(row["both_unique"]),
            int(row["both_unique_zero_missed"]),
            int(row["both_detectable"]),
        ),
        reverse=True,
    )
    ranked_triples = sorted(
        triples,
        key=lambda row: (
            int(row["both_unique"]),
            int(row["both_unique_zero_missed"]),
            int(row["both_detectable"]),
        ),
        reverse=True,
    )
    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "preset_symbols": preset_symbols(),
        "proteins": len(entries),
        "processing_junctions": len(junctions),
        "conventional_propeptide_junctions": sum(
            junction.junction_type == "propeptide" for junction in junctions
        ),
        "pgrn_granulin_junctions": sum(
            junction.junction_type == "granulin" for junction in junctions
        ),
        "minimum_peptide_length": MIN_PEPTIDE_LENGTH,
        "maximum_peptide_length": MAX_PEPTIDE_LENGTH,
        "missed_cleavages": MISSED_CLEAVAGES,
        "joint_quantification": (
            "A protein-unique intact-junction spanning peptide and a protein-unique "
            "mature neo-terminal peptide, potentially from different parallel aliquots."
        ),
        "exact_digest_cut_rule": (
            "A missed-cleavage spanning peptide may report intact precursor, but "
            "the enzyme-generated terminal peptide is not mature-state-specific."
        ),
        "uniqueness": "Human FASTA substring uniqueness with I/L collapsed.",
        "top_two_group_designs": [
            {
                "combination": row["combination"],
                "both_unique": row["both_unique"],
                "both_unique_zero_missed": row["both_unique_zero_missed"],
            }
            for row in ranked_pairs[:3]
        ],
        "top_three_group_designs": [
            {
                "combination": row["combination"],
                "both_unique": row["both_unique"],
                "both_unique_zero_missed": row["both_unique_zero_missed"],
            }
            for row in ranked_triples[:3]
        ],
        "pgrn_source": "https://www.uniprot.org/uniprotkb/P28799/entry",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "propeptide_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(
        f"Finished: {len(entries)} proteins, {len(junctions)} junctions, "
        f"{len(designs)} designs -> {output_dir}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    run(arguments.output)


if __name__ == "__main__":
    main()
