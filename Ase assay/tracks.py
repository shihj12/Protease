"""Genome-browser-style track figure for a single protein.

Top track = molecule-processing features (signal / propeptide / mature chain).
Then one track per protease showing its cut sites, sequence coverage, and — right
at each activation junction — whether it yields a diagnostic pro-form peptide
(green, spans the junction) or cuts ambiguously through it (red).

Kept UI-free so it can be rendered/tested without Streamlit.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D

import analysis

FEATURE_COLORS = {
    "SIGNAL": "#4C78A8",
    "TRANSIT": "#72B7B2",
    "PROPEP": "#F58518",
    "CHAIN": "#9ECAE1",
    "PEPTIDE": "#BAB0AC",
    "INIT_MET": "#D9D9D9",
}
COV_COLOR = "#94A3B8"       # coverage fill
CUT_COLOR = "#334155"       # cut-site ticks
DIAG_COLOR = "#2E7D32"      # diagnostic spanning peptide
AMBIG_COLOR = "#C62828"     # ambiguous junction
TERM_COLOR = "#F58518"      # mature neo-terminal peptide
JUNC_COLOR = "#C62828"      # junction guide line


def _diag_icon(diag_rows):
    """Icon summarising a protease's diagnostic power over propeptide junctions."""
    pj = [d for d in diag_rows if d["is_propeptide_junction"]]
    if not pj:
        return ""
    if any(d["status"].startswith("Diagnostic") for d in pj):
        return "  ✓"   # font-safe (DejaVu Sans)
    return "  ✗"


def _nterm_icon(nterm_rows):
    pj = [r for r in nterm_rows if r["is_propeptide_junction"]]
    if not pj:
        return ""
    return "  ✓" if any(r["detectable"] for r in pj) else "  ✗"


def protein_track_figure(entry, protease_names, missed=2, min_len=7, max_len=35,
                         cov_region="full", propeptide_only=True, mode="diagnostic"):
    """mode='diagnostic' highlights the pro-form spanning peptide (green) or the
    ambiguous junction (red). mode='nterm' highlights the captured neo-N-terminal
    peptide, coloured by whether it is a detectable length."""
    L = entry.length
    junctions = analysis.find_junctions(entry)
    if propeptide_only:
        junctions = [j for j in junctions if j.is_propeptide]
    if mode == "nterm":
        shown_junctions = [j for j in junctions if j.side == "N"]
    else:
        shown_junctions = junctions

    n_rows = 1 + len(protease_names)
    fig_h = 0.9 + 0.5 * n_rows
    fig, ax = plt.subplots(figsize=(12, fig_h))

    row_labels = ["Features"]
    y_of = {}  # row index -> y centre (top row highest)
    for i in range(n_rows):
        y_of[i] = n_rows - i

    band = 0.30  # half-height of a track band

    # --- Junction guide lines across the whole figure ---
    for jn in shown_junctions:
        ax.axvline(jn.pos + 0.5, color=JUNC_COLOR, lw=1.1, ls="--", alpha=0.7,
                   zorder=1)

    # --- Features track (row 0) ---
    yf = y_of[0]
    ax.hlines(yf, 0, L, color="#CBD5E1", lw=1.0, zorder=1)  # backbone
    # Draw chains (light) first, then signal/propep/transit on top.
    order = {"CHAIN": 0, "PEPTIDE": 0, "SIGNAL": 1, "TRANSIT": 1, "PROPEP": 2}
    for f in sorted(entry.features, key=lambda f: order.get(f.kind, 0)):
        color = FEATURE_COLORS.get(f.kind, "#999999")
        h = band if f.kind in ("CHAIN", "PEPTIDE") else band * 1.5
        ax.add_patch(Rectangle((f.start - 1, yf - h), f.end - f.start + 1, 2 * h,
                               facecolor=color, edgecolor="none", zorder=2))
        if f.kind in ("PROPEP", "SIGNAL") and (f.end - f.start + 1) > L * 0.02:
            ax.text((f.start - 1 + f.end) / 2, yf + band * 1.7,
                    "pro" if f.kind == "PROPEP" else "sig",
                    ha="center", va="bottom", fontsize=7, color=color)

    # --- Protease tracks ---
    for k, pname in enumerate(protease_names, start=1):
        yc = y_of[k]
        cov = analysis.coverage(entry, pname, missed, min_len, max_len, cov_region)

        # Coverage fill (contiguous runs).
        flags = cov["covered_flags"]
        start = None
        for i in range(L + 1):
            c = flags[i] if i < L else 0
            if c and start is None:
                start = i
            elif not c and start is not None:
                ax.add_patch(Rectangle((start, yc - band * 0.55), i - start,
                                       band * 1.1, facecolor=COV_COLOR, alpha=0.45,
                                       edgecolor="none", zorder=2))
                start = None

        # Cut-site ticks.
        for s in cov["cut_sites"]:
            ax.vlines(s, yc - band, yc + band, color=CUT_COLOR, lw=0.4,
                      alpha=0.55, zorder=3)

        if mode == "nterm":
            nrows = analysis.nterminomics_capture(entry, pname, missed, min_len, max_len)
            row_labels.append(pname + _nterm_icon(nrows))
            nmap = {r["neo_nterm_pos"]: r for r in nrows}
            for jn in shown_junctions:
                r = nmap.get(jn.mature_start)
                if r is None:
                    continue
                col = DIAG_COLOR if r["detectable"] else AMBIG_COLOR
                ax.add_patch(Rectangle((r["start"] - 1, yc - band),
                                       r["end"] - r["start"] + 1, 2 * band,
                                       facecolor=col, edgecolor="white", lw=0.4,
                                       zorder=4))
        else:
            diag = analysis.junction_diagnostics(entry, pname, missed, min_len, max_len)
            row_labels.append(pname + _diag_icon(diag))
            dmap = {(d["junction"], d["side"]): d for d in diag}
            for jn in shown_junctions:
                d = dmap.get((f"{jn.pos}|{jn.pos + 1}", jn.side))
                if d is None:
                    continue
                if d["status"].startswith("Diagnostic"):
                    ax.add_patch(Rectangle(
                        (d["spanning_start"] - 1, yc - band), d["spanning_len"],
                        2 * band, facecolor=DIAG_COLOR, edgecolor="white", lw=0.4,
                        zorder=4))
                else:  # ambiguous: mark junction + show the ambiguous mature peptide
                    m_len = d["mature_term_len"]
                    m_start = jn.pos + 1 if jn.side == "N" else jn.pos - m_len + 1
                    ax.add_patch(Rectangle((m_start - 1, yc - band), m_len, 2 * band,
                                           facecolor=TERM_COLOR, alpha=0.6,
                                           edgecolor="none", zorder=3))
                    ax.vlines(jn.pos + 0.5, yc - band * 1.2, yc + band * 1.2,
                              color=AMBIG_COLOR, lw=2.2, zorder=5)

    # --- Axes cosmetics ---
    ax.set_xlim(-L * 0.01, L * 1.01)
    ax.set_ylim(0.3, n_rows + 1.2)
    ax.set_yticks([y_of[i] for i in range(n_rows)])
    ax.set_yticklabels(row_labels, fontsize=9)
    ax.set_xlabel("Residue position")
    title = f"{entry.gene or entry.name}  ({entry.accession})  — {L} aa"
    if not shown_junctions:
        title += "  — no propeptide junction"
    ax.set_title(title, fontsize=11, loc="left")
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)

    legend = [
        Patch(facecolor=FEATURE_COLORS["PROPEP"], label="propeptide"),
        Patch(facecolor=FEATURE_COLORS["SIGNAL"], label="signal"),
        Patch(facecolor=FEATURE_COLORS["CHAIN"], label="mature chain"),
        Patch(facecolor=COV_COLOR, alpha=0.45, label="coverage"),
    ]
    if mode == "nterm":
        legend += [
            Patch(facecolor=DIAG_COLOR, label="captured N-term peptide (detectable)"),
            Patch(facecolor=AMBIG_COLOR, label="captured N-term (too short/long)"),
        ]
    else:
        legend += [
            Patch(facecolor=DIAG_COLOR, label="diagnostic peptide"),
            Line2D([0], [0], color=AMBIG_COLOR, lw=2.2, label="ambiguous (cuts junction)"),
        ]
    legend.append(
        Line2D([0], [0], color=JUNC_COLOR, lw=1.1, ls="--", label="activation junction"))
    ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, -0.25),
              ncol=4, fontsize=7.5, frameon=False)
    fig.tight_layout()
    return fig
