"""Generate Nature-observed SILAC and pro-peptide reference analyses."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from coverage_analysis.nature_propeptide import quantify_nature_propeptides
from coverage_analysis.nature_reference import nature_coverage, nature_label_coverage
from coverage_analysis.plots import silac_panel_plot
from coverage_analysis.propeptide_plots import outcome_plot, pgrn_heatmap, protein_heatmap
from coverage_analysis.propeptide_quantification import load_preset_entries
from coverage_analysis.protease_coverage import (
    LABELS,
    KR_PLUS_ONE_LABELS,
    NATURE_ENZYMES,
    nature_designs,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CACHE = Path(__file__).resolve().parent / "reference_cache"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "nature_reference_output"


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def _family_designs():
    singles, pairs, triples = nature_designs()
    return singles, pairs, triples, [
        *[("single", name, enzymes) for name, enzymes in singles],
        *[("two", name, enzymes) for name, enzymes in pairs],
        *[("three", name, enzymes) for name, enzymes in triples],
    ]


def _rank_reference(rows, designs):
    wanted = {name for name, _ in designs}
    return sorted(
        [row for row in rows if row["combination"] in wanted],
        key=lambda row: (
            float(row["residue_weighted_pct"]), float(row["median_pct"])
        ),
        reverse=True,
    )


def run(output_dir: Path, cache_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    singles, pairs, triples, propeptide_designs = _family_designs()
    coverage_designs = singles + pairs + triples + [("All six", NATURE_ENZYMES)]
    print("Loading Nature-observed coverage rankings")
    coverage_rows, coverage_provenance = nature_coverage(
        cache_dir, coverage_designs, print
    )
    ranked_pairs = _rank_reference(coverage_rows, pairs)
    ranked_triples = _rank_reference(coverage_rows, triples)
    silac_designs = [
        (
            str(row["combination"]),
            tuple(str(row["enzymes"]).split(";")),
        )
        for row in ranked_pairs[:3] + ranked_triples[:3]
    ]

    print("Calculating Nature-observed SILAC label coverage")
    silac_rows, silac_provenance = nature_label_coverage(
        cache_dir, silac_designs, LABELS, print
    )
    _write_csv(output_dir / "nature_silac_quantifiable_coverage.csv", silac_rows)
    silac_panel_plot(
        silac_rows,
        [name for name, _ in silac_designs],
        output_dir / "graph_1_nature_silac_top_combinations.png",
    )
    print("Calculating Nature-observed K+R-plus-one-label coverage")
    kr_plus_one_rows, kr_plus_one_provenance = nature_label_coverage(
        cache_dir, silac_designs, KR_PLUS_ONE_LABELS, print
    )
    kr_baselines = {
        str(row["combination"]): float(row["residue_weighted_pct"])
        for row in silac_rows
        if row["label"] == "K+R"
    }
    for row in kr_plus_one_rows:
        baseline = kr_baselines[str(row["combination"])]
        row["kr_baseline_pct"] = baseline
        row["gain_over_kr_pct"] = float(row["residue_weighted_pct"]) - baseline
    _write_csv(
        output_dir / "nature_silac_kr_plus_one_coverage.csv", kr_plus_one_rows
    )
    silac_panel_plot(
        kr_plus_one_rows,
        [name for name, _ in silac_designs],
        output_dir / "graph_7_nature_silac_kr_plus_one.png",
        label_order=list(KR_PLUS_ONE_LABELS),
        highlight_label="K+R+L",
    )

    print("Matching pro-peptide state candidates to Nature-observed peptides")
    entries = load_preset_entries()
    peptide_cache = cache_dir / "nature_peptides_v1.tsv.gz"
    junctions, detail_rows, summary_rows, matched_sequences = (
        quantify_nature_propeptides(entries, propeptide_designs, peptide_cache)
    )
    _write_csv(output_dir / "nature_propeptide_quantifiability.csv", detail_rows)
    _write_csv(output_dir / "nature_propeptide_design_summary.csv", summary_rows)

    single_rows = [row for row in summary_rows if row["family"] == "single"]
    pair_rows = sorted(
        [row for row in summary_rows if row["family"] == "two"],
        key=lambda row: (
            int(row["both_unique"]), int(row["intact_observed_unique"])
        ),
        reverse=True,
    )
    triple_rows = sorted(
        [row for row in summary_rows if row["family"] == "three"],
        key=lambda row: (
            int(row["both_unique"]), int(row["intact_observed_unique"])
        ),
        reverse=True,
    )
    outcome_plot(
        single_rows,
        output_dir / "graph_2_nature_propeptide_group_outcomes.png",
    )
    outcome_plot(
        pair_rows,
        output_dir / "graph_3_nature_propeptide_two_groups.png",
    )
    outcome_plot(
        triple_rows,
        output_dir / "graph_4_nature_propeptide_three_groups.png",
    )
    single_names = [name for name, _ in singles]
    protein_heatmap(
        detail_rows,
        single_names,
        output_dir / "graph_5_nature_propeptide_protein_heatmap.png",
        metric="intact_unique",
        colorbar_label="Observed intact junctions (%)",
    )
    pgrn_heatmap(
        detail_rows,
        single_names,
        output_dir / "graph_6_nature_pgrn_junction_heatmap.png",
        metric="intact_unique",
    )

    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "output_scope": "Nature reference SILAC and pro-peptide analyses",
        "dataset": "PXD024364 / MSV000086944",
        "nature_enzymes": list(NATURE_ENZYMES),
        "silac_top_designs": [name for name, _ in silac_designs],
        "propeptide_proteins": len(entries),
        "propeptide_junctions": len(junctions),
        "candidate_sequences_observed": matched_sequences,
        "reference_search_specificity": "specific cleavage",
        "mature_state_limitation": (
            "The deposited MaxQuant search required specific cleavage. "
            "Semi-specific mature neo-terminal peptides were not searched, so "
            "absence of mature-state evidence is not evidence of absent processing."
        ),
        "coverage_provenance": coverage_provenance,
        "silac_provenance": silac_provenance,
        "kr_plus_one_provenance": kr_plus_one_provenance,
        "sources": {
            "study": "https://doi.org/10.1038/s41587-023-01714-x",
            "dataset": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD024364",
        },
    }
    (output_dir / "nature_extension_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(f"Finished Nature extension outputs: {output_dir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    arguments = parser.parse_args()
    run(arguments.output, arguments.cache)


if __name__ == "__main__":
    main()
