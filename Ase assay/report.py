"""Assemble a multi-page PDF report of the current analysis.

Uses matplotlib's PdfPages (no extra dependencies). Pages:
  1. Cover — parameters and the queried proteins.
  2. Protease scorecard across the whole set.
  3. SILAC optimal-label-per-protease + protease x label matrix.
  4. One diagnostic track figure per protein (with a propeptide junction).
  5. One N-terminomics track figure per protein (N-terminal junctions).

`build_pdf(...)` returns the PDF as bytes, ready for st.download_button.
"""

import datetime
import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import analysis
import tracks

LANDSCAPE = (11.0, 8.5)


def _cover_page(pdf, entries, params):
    fig = plt.figure(figsize=LANDSCAPE)
    fig.text(0.07, 0.90, "Propeptide Activation / Protease Explorer",
             fontsize=20, weight="bold")
    fig.text(0.07, 0.855, "Analysis report", fontsize=13, color="#555555")
    today = datetime.date.today().isoformat()
    lines = [
        f"Generated: {today}",
        f"Proteases: {', '.join(params['proteases'])}",
        f"Missed cleavages: {params['missed']}    "
        f"Detectable length: {params['min_len']}-{params['max_len']} aa",
        f"Coverage denominator: {params['cov_region']}    "
        f"Junctions: {'propeptide only' if params['prop_only'] else 'all processing'}",
        f"Candidate SILAC labels: {', '.join(params['candidate_aa'])}",
        f"Proteins queried: {len(entries)}",
    ]
    fig.text(0.07, 0.80, "\n".join(lines), fontsize=11, va="top", linespacing=1.8)

    # Protein list.
    plist = [f"{e.gene or e.name} ({e.accession}) — {e.protein}" for e in entries]
    fig.text(0.07, 0.52, "Proteins:", fontsize=11, weight="bold", va="top")
    # two columns if many
    half = (len(plist) + 1) // 2
    fig.text(0.09, 0.48, "\n".join(plist[:half]), fontsize=8.5, va="top",
             family="DejaVu Sans", linespacing=1.5)
    if len(plist) > half:
        fig.text(0.52, 0.48, "\n".join(plist[half:]), fontsize=8.5, va="top",
                 family="DejaVu Sans", linespacing=1.5)
    fig.text(0.07, 0.04, "In-silico analysis for experiment planning. 'Detectable' "
             "= peptide length only; no m/z or ionisation modelling.",
             fontsize=7.5, color="#888888")
    pdf.savefig(fig)
    plt.close(fig)


def _table_page(pdf, title, columns, rows, note="", max_rows=26, col_widths=None):
    """Render rows (list of lists) as one or more landscape table pages."""
    if not rows:
        return
    for start in range(0, len(rows), max_rows):
        chunk = rows[start:start + max_rows]
        fig, ax = plt.subplots(figsize=LANDSCAPE)
        ax.axis("off")
        suffix = "" if len(rows) <= max_rows else \
            f"  (rows {start + 1}-{start + len(chunk)} of {len(rows)})"
        ax.set_title(title + suffix, fontsize=14, weight="bold", loc="left", pad=18)
        tbl = ax.table(cellText=[[str(c) for c in r] for r in chunk],
                       colLabels=columns, loc="upper center",
                       colWidths=col_widths, cellLoc="left")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(8.5)
        tbl.scale(1, 1.4)
        for (r, _c), cell in tbl.get_celld().items():
            if r == 0:
                cell.set_text_props(weight="bold", color="white")
                cell.set_facecolor("#334155")
            else:
                cell.set_facecolor("#F8FAFC" if r % 2 else "#FFFFFF")
            cell.set_edgecolor("#E2E8F0")
        if note:
            fig.text(0.07, 0.06, note, fontsize=8, color="#666666")
        pdf.savefig(fig)
        plt.close(fig)


def build_pdf(entries, proteases, missed=2, min_len=7, max_len=35,
              cov_region="full", candidate_aa=None, prop_only=True) -> bytes:
    candidate_aa = candidate_aa or (["KR"] + analysis.LABELABLE_AA)
    params = {
        "proteases": proteases, "missed": missed, "min_len": min_len,
        "max_len": max_len, "cov_region": cov_region, "prop_only": prop_only,
        "candidate_aa": candidate_aa,
    }
    buf = io.BytesIO()
    with PdfPages(buf) as pdf:
        _cover_page(pdf, entries, params)

        # Protease scorecard.
        score = analysis.protease_scorecard(entries, proteases, missed, min_len,
                                            max_len, prop_only)
        if score and score[0]["total_junctions"]:
            tot = score[0]["total_junctions"]
            rows = [[s["protease"],
                     f"{s['diagnostic']}/{tot} ({s['pct_diagnostic']:.0f}%)",
                     s["ambiguous"], s["no_spanning_pep"], s["nterm_detectable"]]
                    for s in score]
            _table_page(
                pdf, f"Protease scorecard — {len(entries)} proteins, {tot} junctions",
                ["Protease", "Diagnostic", "Ambiguous", "No span", "N-term detect."],
                rows, note="Diagnostic = a fully-specific peptide spans the junction "
                           "(the protease can confirm activation).")

        # SILAC label optimisation.
        best_rows, matrix = [], {}
        for p in proteases:
            opt = analysis.silac_label_optimization(entries, p, missed, min_len,
                                                    max_len, candidate_aa, prop_only)
            if not opt or not opt[0]["total_peptides"]:
                continue
            matrix[p] = {o["aa"]: o["pct_quantifiable"] for o in opt}
            t = opt[0]
            best_rows.append([p, t["amino_acid"],
                              f"{t['quantifiable']}/{t['total_peptides']} "
                              f"({t['pct_quantifiable']:.0f}%)",
                              t["avg_labels_per_pep"]])
        if best_rows:
            _table_page(pdf, "SILAC — optimal label per protease",
                        ["Protease", "Optimal label", "Quantifiable", "Avg/pep"],
                        best_rows,
                        note="Which heavy amino acid makes the most diagnostic "
                             "peptides quantifiable. Only essential residues tested.")
            # Matrix page.
            aa_keys = ["+".join(sorted(set(t.upper()))) for t in candidate_aa]
            aa_keys = [a for a in dict.fromkeys(aa_keys)]
            mcols = ["Protease"] + aa_keys
            mrows = [[p] + [f"{matrix[p].get(a, 0):.0f}" for a in aa_keys]
                     for p in matrix]
            _table_page(pdf, "SILAC — % diagnostic peptides quantifiable (protease x label)",
                        mcols, mrows,
                        note="Higher = more diagnostic peptides carry >=1 heavy residue.")

        # Diagnostic track figures.
        for e in entries:
            jns = [j for j in analysis.find_junctions(e) if j.is_propeptide]
            if prop_only and not jns:
                continue
            fig = tracks.protein_track_figure(e, proteases, missed, min_len, max_len,
                                              cov_region, prop_only, mode="diagnostic")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

        # N-terminomics track figures (proteins with an N-terminal junction).
        for e in entries:
            jns = [j for j in analysis.find_junctions(e)
                   if j.side == "N" and (j.is_propeptide or not prop_only)]
            if not jns:
                continue
            fig = tracks.protein_track_figure(e, proteases, missed, min_len, max_len,
                                              cov_region, prop_only, mode="nterm")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    return buf.getvalue()
