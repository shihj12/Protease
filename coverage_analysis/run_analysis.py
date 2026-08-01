"""Generate all requested coverage figures and their source tables."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from coverage_analysis.nature_reference import nature_coverage
from coverage_analysis.plots import horizontal_coverage_plot, silac_panel_plot
from coverage_analysis.protease_coverage import (
    LABELS,
    MAX_PEPTIDE_LENGTH,
    MIN_PEPTIDE_LENGTH,
    MISSED_CLEAVAGES,
    NATURE_ENZYMES,
    RULES,
    digest_proteome,
    labelled_design_statistics,
    nature_designs,
    rank_designs,
    read_fasta,
    requested_designs,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FASTA = (
    ROOT
    / "ase_assay"
    / "uniprotkb_human_AND_model_organism_9606_2026_07_21.fasta"
)
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "output"
DEFAULT_CACHE = Path(__file__).resolve().parent / "reference_cache"


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


def _progress(enzyme: str, completed: int, total: int) -> None:
    print(f"  {enzyme}: {completed:,}/{total:,} proteins", flush=True)


def _ranked_subset(
    all_rows: list[dict[str, object]],
    requested: list[tuple[str, tuple[str, ...]]],
) -> list[dict[str, object]]:
    wanted = {name for name, _ in requested}
    return [row for row in all_rows if str(row["combination"]) in wanted]


def run(fasta: Path, output_dir: Path, cache_dir: Path, skip_nature: bool) -> None:
    print(f"Loading theoretical proteome: {fasta}")
    proteins = read_fasta(fasta)
    total_residues = sum(len(protein.sequence) for protein in proteins)
    print(f"Loaded {len(proteins):,} proteins / {total_residues:,} residues")

    requested_singles, requested_pairs, requested_triples = requested_designs()
    enzyme_names = list(RULES)
    print("Calculating theoretical detectable-peptide residue masks")
    enzyme_masks, _, peptide_counts = digest_proteome(
        proteins, enzyme_names, progress=_progress
    )
    single_rows = rank_designs(proteins, enzyme_masks, requested_singles)
    pair_rows = rank_designs(proteins, enzyme_masks, requested_pairs)
    triple_rows = rank_designs(proteins, enzyme_masks, requested_triples)
    for family, rows in (("single", single_rows), ("two", pair_rows), ("three", triple_rows)):
        for row in rows:
            row["family"] = family
    theoretical_rows = single_rows + pair_rows + triple_rows
    _write_csv(output_dir / "theoretical_coverage.csv", theoretical_rows)

    fixed_order = {name: i for i, (name, _) in enumerate(requested_singles)}
    graph1_rows = sorted(single_rows, key=lambda row: fixed_order[str(row["combination"])])
    horizontal_coverage_plot(
        graph1_rows,
        output_dir / "graph_1_theoretical_groups.png",
        preserve_order=True,
    )
    horizontal_coverage_plot(
        pair_rows, output_dir / "graph_2_theoretical_two_groups.png"
    )
    horizontal_coverage_plot(
        triple_rows, output_dir / "graph_3_theoretical_three_groups.png"
    )

    top_pairs = pair_rows[:3]
    top_triples = triple_rows[:3]
    silac_designs = [
        (str(row["combination"]), tuple(str(row["enzymes"]).split(";")))
        for row in top_pairs + top_triples
    ]
    silac_enzymes = sorted({enzyme for _, enzymes in silac_designs for enzyme in enzymes})
    print("Calculating label-bearing theoretical peptide masks for the top designs")
    _, label_masks, _ = digest_proteome(
        proteins, silac_enzymes, labels=LABELS, progress=_progress
    )
    silac_rows = labelled_design_statistics(proteins, label_masks, silac_designs)
    _write_csv(output_dir / "silac_quantifiable_coverage.csv", silac_rows)
    silac_panel_plot(
        silac_rows,
        [combination for combination, _ in silac_designs],
        output_dir / "graph_7_silac_top_combinations.png",
    )

    nature_provenance: dict[str, object] | None = None
    if not skip_nature:
        nature_singles, nature_pairs, nature_triples = nature_designs()
        validation_design = [("All six", NATURE_ENZYMES)]
        all_nature_designs = nature_singles + nature_pairs + nature_triples + validation_design
        print("Calculating empirical coverage from the original Nature MaxQuant output")
        nature_rows, nature_provenance = nature_coverage(
            cache_dir, all_nature_designs, print
        )
        for family, designs in (
            ("single", nature_singles),
            ("two", nature_pairs),
            ("three", nature_triples),
        ):
            names = {name for name, _ in designs}
            subset = [row for row in nature_rows if row["combination"] in names]
            subset.sort(
                key=lambda row: (
                    float(row["residue_weighted_pct"]), float(row["median_pct"])
                ),
                reverse=True,
            )
            for rank, row in enumerate(subset, start=1):
                row["rank"] = rank
                row["family"] = family

        _write_csv(output_dir / "nature_observed_coverage.csv", nature_rows)
        single_order = {name: i for i, (name, _) in enumerate(nature_singles)}
        nature_single_rows = _ranked_subset(nature_rows, nature_singles)
        nature_single_rows.sort(key=lambda row: single_order[str(row["combination"])])
        nature_pair_rows = sorted(
            _ranked_subset(nature_rows, nature_pairs),
            key=lambda row: float(row["residue_weighted_pct"]),
            reverse=True,
        )
        nature_triple_rows = sorted(
            _ranked_subset(nature_rows, nature_triples),
            key=lambda row: float(row["residue_weighted_pct"]),
            reverse=True,
        )
        horizontal_coverage_plot(
            nature_single_rows,
            output_dir / "graph_4_nature_observed_groups.png",
            axis_label="Observed sequence coverage (%)",
            preserve_order=True,
        )
        horizontal_coverage_plot(
            nature_pair_rows,
            output_dir / "graph_5_nature_observed_two_groups.png",
            axis_label="Observed sequence coverage (%)",
        )
        horizontal_coverage_plot(
            nature_triple_rows,
            output_dir / "graph_6_nature_observed_three_groups.png",
            axis_label="Observed sequence coverage (%)",
        )

    metadata = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "theoretical_fasta": str(fasta.resolve()),
        "theoretical_proteins": len(proteins),
        "theoretical_residues": total_residues,
        "minimum_peptide_length": MIN_PEPTIDE_LENGTH,
        "maximum_peptide_length": MAX_PEPTIDE_LENGTH,
        "missed_cleavages": MISSED_CLEAVAGES,
        "coverage_metric_in_figures": "residue-weighted sequence coverage",
        "parallel_digest_rule": "union of independently digested aliquot residue masks",
        "peptide_instances_by_enzyme": peptide_counts,
        "protease_rules": {name: rule.description for name, rule in RULES.items()},
        "top_two_group_designs": [row["combination"] for row in top_pairs],
        "top_three_group_designs": [row["combination"] for row in top_triples],
        "nature": nature_provenance,
        "citations": {
            "nature_study": "https://doi.org/10.1038/s41587-023-01714-x",
            "nature_dataset": "https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD024364",
            "spectronaut_manual": "https://biognosys.com/content/uploads/2024/09/Spectronaut-19-manual-v4.pdf",
            "expasy_rules": "https://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html",
            "proalanase": "https://www.promega.com/products/mass-spectrometry/proteases-and-surfactants/proalanase-mass-spec-grade/",
            "elastase_rule": "https://docs.thermofisher.com/r/Proteome-Discoverer-3.1-User-Guide/en-US1324471691v1",
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(f"Finished. Outputs: {output_dir}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument(
        "--skip-nature",
        action="store_true",
        help="Generate theoretical figures only; do not access PXD024364.",
    )
    arguments = parser.parse_args()
    run(arguments.fasta, arguments.output, arguments.cache, arguments.skip_nature)


if __name__ == "__main__":
    main()
