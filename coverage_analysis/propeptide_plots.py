"""Minimal figures for paired precursor/mature-state quantification."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


TEXT = "#102a43"
COLORS = {
    "both_unique": "#2f6f8f",
    "intact_only": "#7aa6b8",
    "mature_only": "#c5d5dc",
    "neither": "#e8ecef",
}


def _save(figure, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=300, bbox_inches="tight", facecolor="white")
    figure.savefig(destination.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _minimal(axis) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(False)
    axis.tick_params(axis="y", length=0, colors=TEXT, labelsize=8.5)
    axis.tick_params(axis="x", colors=TEXT, labelsize=8)
    axis.xaxis.label.set_color(TEXT)


def outcome_plot(rows: list[dict[str, object]], destination: Path) -> None:
    labels = [str(row["combination"]) for row in rows][::-1]
    figure, axis = plt.subplots(figsize=(9.2, 5.3))
    left = [0] * len(rows)
    keys = ("both_unique", "intact_only", "mature_only", "neither")
    names = ("Both", "Intact only", "Mature only", "Neither")
    for key, name in zip(keys, names):
        values = [int(row[key]) for row in rows][::-1]
        bars = axis.barh(
            labels, values, left=left, height=0.64, color=COLORS[key], label=name
        )
        for bar, value in zip(bars, values):
            if value:
                axis.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_y() + bar.get_height() / 2,
                    str(value),
                    ha="center",
                    va="center",
                    fontsize=7.5,
                    color=TEXT,
                )
        left = [previous + value for previous, value in zip(left, values)]
    _minimal(axis)
    axis.set_xlabel("Processing junctions (n)", fontsize=9)
    axis.legend(frameon=False, fontsize=8, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.1))
    figure.tight_layout()
    _save(figure, destination)


def combination_plot(rows: list[dict[str, object]], destination: Path) -> None:
    rows = sorted(
        rows,
        key=lambda row: (
            int(row["both_unique"]),
            int(row["both_unique_zero_missed"]),
            int(row["both_detectable"]),
        ),
        reverse=True,
    )
    labels = [str(row["combination"]) for row in rows][::-1]
    fully_cleaved = [int(row["both_unique_zero_missed"]) for row in rows][::-1]
    rescued = [
        int(row["both_unique"]) - int(row["both_unique_zero_missed"])
        for row in rows
    ][::-1]
    total = int(rows[0]["total_junctions"]) if rows else 0
    unresolved = [
        total - fully - rescue for fully, rescue in zip(fully_cleaved, rescued)
    ]
    figure, axis = plt.subplots(figsize=(9.2, max(3.3, 0.38 * len(rows) + 1)))
    left = [0] * len(rows)
    for values, color, name in (
        (fully_cleaved, COLORS["both_unique"], "0 missed"),
        (rescued, COLORS["intact_only"], "+ missed"),
        (unresolved, COLORS["neither"], "Not joint"),
    ):
        bars = axis.barh(
            labels, values, left=left, color=color, height=0.64, label=name
        )
        for bar, value in zip(bars, values):
            if value:
                axis.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_y() + bar.get_height() / 2,
                    str(value),
                    ha="center",
                    va="center",
                    fontsize=7.5,
                    color=TEXT,
                )
        left = [previous + value for previous, value in zip(left, values)]
    _minimal(axis)
    axis.set_xlabel("Processing junctions (n)", fontsize=9)
    axis.set_xlim(0, total)
    axis.legend(
        frameon=False, fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.06)
    )
    figure.tight_layout()
    _save(figure, destination)


def protein_heatmap(
    detail_rows: list[dict[str, object]],
    combinations: list[str],
    destination: Path,
    metric: str = "both_unique",
    colorbar_label: str = "Jointly quantifiable (%)",
) -> None:
    selected = [row for row in detail_rows if row["combination"] in combinations]
    genes = list(dict.fromkeys(str(row["gene"]) for row in selected))
    total_by_gene = {
        gene: len(
            {
                (row["position"], row["side"])
                for row in selected
                if row["gene"] == gene and row["combination"] == combinations[0]
            }
        )
        for gene in genes
    }
    matrix: list[list[float]] = []
    annotations: list[list[str]] = []
    for gene in genes:
        totals = total_by_gene[gene]
        values = []
        labels = []
        for combination in combinations:
            count = sum(
                bool(row[metric])
                for row in selected
                if row["gene"] == gene and row["combination"] == combination
            )
            values.append(100.0 * count / totals if totals else 0.0)
            labels.append(f"{count}/{totals}")
        matrix.append(values)
        annotations.append(labels)

    cmap = LinearSegmentedColormap.from_list("coverage", ["#f3f5f7", "#2f6f8f"])
    figure, axis = plt.subplots(figsize=(12.0, max(8.5, 0.34 * len(genes) + 1.8)))
    image = axis.imshow(matrix, cmap=cmap, vmin=0, vmax=100, aspect="auto")
    axis.set_xticks(range(len(combinations)), combinations, rotation=45, ha="right", fontsize=8)
    axis.set_yticks(range(len(genes)), genes, fontsize=8)
    axis.set_xlabel("Protease group", fontsize=9, color=TEXT)
    axis.set_ylabel("Protein", fontsize=9, color=TEXT)
    axis.tick_params(length=0, colors=TEXT)
    for row_index, labels in enumerate(annotations):
        for column_index, label in enumerate(labels):
            value = matrix[row_index][column_index]
            axis.text(
                column_index,
                row_index,
                label,
                ha="center",
                va="center",
                fontsize=6.4,
                color="white" if value >= 55 else TEXT,
            )
    colorbar = figure.colorbar(image, ax=axis, fraction=0.025, pad=0.02)
    colorbar.set_label(colorbar_label, fontsize=8, color=TEXT)
    colorbar.ax.tick_params(labelsize=7, colors=TEXT)
    for spine in axis.spines.values():
        spine.set_visible(False)
    figure.tight_layout()
    _save(figure, destination)


def pgrn_heatmap(
    detail_rows: list[dict[str, object]],
    combinations: list[str],
    destination: Path,
    metric: str = "both_unique",
) -> None:
    selected = [
        row
        for row in detail_rows
        if row["gene"] == "GRN" and row["combination"] in combinations
    ]
    junction_keys = list(
        dict.fromkeys(
            (int(row["position"]), str(row["side"]), str(row["product"]))
            for row in selected
        )
    )
    matrix = [
        [
            int(
                any(
                    row[metric]
                    for row in selected
                    if int(row["position"]) == position
                    and row["side"] == side
                    and row["combination"] == combination
                )
            )
            for combination in combinations
        ]
        for position, side, _ in junction_keys
    ]
    labels = [
        f"{product} {side} {position}|{position + 1}"
        for position, side, product in junction_keys
    ]
    cmap = LinearSegmentedColormap.from_list("binary", ["#edf1f3", "#2f6f8f"])
    figure, axis = plt.subplots(figsize=(11.5, 6.8))
    axis.imshow(matrix, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    axis.set_xticks(range(len(combinations)), combinations, rotation=45, ha="right", fontsize=8)
    axis.set_yticks(range(len(labels)), labels, fontsize=7.5)
    axis.set_xlabel("Protease group", fontsize=9, color=TEXT)
    axis.set_ylabel("PGRN processing junction", fontsize=9, color=TEXT)
    axis.tick_params(length=0, colors=TEXT)
    for row_index, values in enumerate(matrix):
        for column_index, value in enumerate(values):
            axis.text(
                column_index,
                row_index,
                str(value),
                ha="center",
                va="center",
                fontsize=7,
                color="white" if value else TEXT,
            )
    for spine in axis.spines.values():
        spine.set_visible(False)
    figure.tight_layout()
    _save(figure, destination)
