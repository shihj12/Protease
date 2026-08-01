"""Minimal plot functions for the requested figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["svg.hashsalt"] = "protease-coverage"


COLOR = "#486581"
LIGHT_COLOR = "#829ab1"
TEXT_COLOR = "#102a43"


def _save_svg_without_trailing_space(figure, destination: Path) -> None:
    figure.savefig(destination, bbox_inches="tight", facecolor="white")
    text = destination.read_text(encoding="utf-8")
    lines = [line.rstrip() for line in text.splitlines()]
    with destination.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\n".join(lines) + "\n")


def _style_axis(axis) -> None:
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_visible(False)
    axis.tick_params(axis="y", length=0, colors=TEXT_COLOR, labelsize=9)
    axis.tick_params(axis="x", colors=TEXT_COLOR, labelsize=8)
    axis.xaxis.label.set_color(TEXT_COLOR)


def horizontal_coverage_plot(
    rows: list[dict[str, object]],
    destination: Path,
    metric: str = "residue_weighted_pct",
    axis_label: str = "Sequence coverage (%)",
    preserve_order: bool = False,
) -> None:
    if not preserve_order:
        rows = sorted(rows, key=lambda row: float(row[metric]), reverse=True)
    labels = [str(row["combination"]) for row in rows][::-1]
    values = [float(row[metric]) for row in rows][::-1]
    height = max(3.2, 0.37 * len(rows) + 1.0)
    figure, axis = plt.subplots(figsize=(8.6, height))
    bars = axis.barh(labels, values, color=COLOR, height=0.64)
    _style_axis(axis)
    axis.set_xlabel(axis_label, fontsize=9)
    maximum = max(values, default=0.0)
    axis.set_xlim(0, min(100.0, maximum + max(4.0, maximum * 0.08)))
    axis.bar_label(
        bars, labels=[f"{value:.1f}" for value in values], padding=3, fontsize=8
    )
    figure.tight_layout()
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=300, bbox_inches="tight", facecolor="white")
    _save_svg_without_trailing_space(figure, destination.with_suffix(".svg"))
    plt.close(figure)


def silac_panel_plot(
    rows: list[dict[str, object]],
    combinations: list[str],
    destination: Path,
    metric: str = "residue_weighted_pct",
    label_order: list[str] | None = None,
    highlight_label: str | None = "K+R",
) -> None:
    by_combination: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        by_combination.setdefault(str(row["combination"]), []).append(row)
    label_order = label_order or ["K+R", "L", "T", "I", "V", "H"]
    figure, axes = plt.subplots(2, 3, figsize=(14.2, 7.7), sharex=True)
    for axis, combination in zip(axes.flat, combinations):
        lookup = {str(row["label"]): row for row in by_combination[combination]}
        labels = label_order[::-1]
        values = [float(lookup[label][metric]) for label in labels]
        bars = axis.barh(labels, values, color=LIGHT_COLOR, height=0.62)
        if highlight_label in labels:
            bars[labels.index(highlight_label)].set_color(COLOR)
        _style_axis(axis)
        axis.set_title(combination, fontsize=9.5, color=TEXT_COLOR, pad=8)
        axis.bar_label(
            bars,
            labels=[f"{value:.1f}" for value in values],
            padding=2,
            fontsize=7.5,
        )
        axis.set_xlim(0, 100)
    for axis in axes[1, :]:
        axis.set_xlabel("Quantifiable sequence coverage (%)", fontsize=8.5)
    figure.tight_layout(w_pad=2.0, h_pad=2.0)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=300, bbox_inches="tight", facecolor="white")
    _save_svg_without_trailing_space(figure, destination.with_suffix(".svg"))
    plt.close(figure)
