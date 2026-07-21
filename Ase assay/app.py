"""Propeptide Activation / Protease Explorer — Streamlit app.

Run:  streamlit run app.py

Enter gene symbols (CTSD) or UniProt accessions, pick proteases, and read a
track-per-protease view of which enzymes can actually capture propeptide
activation — plus coverage, SILAC quantifiability, and N-terminomics capture.
"""

import pandas as pd
import streamlit as st

import analysis
import report
import tracks
import uniprot
from presets import PRESETS
from proteases import protease_names, PROTEASES

st.set_page_config(page_title="Propeptide / Protease Explorer", layout="wide")

# --- Sidebar -------------------------------------------------------------------

st.sidebar.title("Inputs")

if "acc_text" not in st.session_state:
    st.session_state.acc_text = " ".join(PRESETS["Lysosomal hydrolases"])


def _set_preset(name):
    st.session_state.acc_text = " ".join(PRESETS[name])


st.sidebar.caption("Load a preset:")
pcols = st.sidebar.columns(len(PRESETS))
for col, name in zip(pcols, PRESETS):
    col.button(name, on_click=_set_preset, args=(name,), width="stretch")

st.sidebar.text_area(
    "Gene symbols or UniProt accessions", key="acc_text", height=110,
    help="Space / comma / newline separated. e.g. CTSD CTSB P07339",
)
use_cache = st.sidebar.checkbox("Use local cache", value=True)
fetch = st.sidebar.button("Fetch / refresh", type="primary")

st.sidebar.markdown("---")
sel_proteases = st.sidebar.multiselect(
    "Proteases", protease_names(),
    default=["Trypsin", "Lys-C", "Glu-C (D/E)", "Asp-N", "Chymotrypsin"],
)
missed = st.sidebar.slider("Missed cleavages", 0, 4, 2)
min_len, max_len = st.sidebar.slider("Detectable peptide length", 4, 60, (7, 35))

st.sidebar.markdown("---")
label_res_text = st.sidebar.text_input(
    "Candidate SILAC labels", value="KR," + ",".join(analysis.LABELABLE_AA),
    help="Labels to test. Use one letter for a single amino acid, or several "
         "together (e.g. KR) for a multi-label scheme. Defaults to the standard "
         "Lys+Arg pair plus each essential amino acid singly.",
)
candidate_aa = [
    r.strip().upper() for r in label_res_text.split(",") if r.strip()
] or (["KR"] + analysis.LABELABLE_AA)
cov_region = st.sidebar.radio(
    "Coverage denominator", ["full", "chain"], horizontal=True,
    help="'chain' restricts % coverage to mature-chain residues only.",
)

# --- Fetch ---------------------------------------------------------------------

if fetch:
    tokens = uniprot.parse_accessions(st.session_state.acc_text)
    prog = st.progress(0.0, text="Fetching from UniProt...")
    entries, errors = [], []
    for i, tok in enumerate(tokens):
        try:
            entries.append(uniprot.fetch_input(tok, use_cache=use_cache))
        except ValueError as e:
            errors.append((tok, str(e)))
        prog.progress((i + 1) / max(len(tokens), 1), text=f"Fetched {tok}")
    prog.empty()
    st.session_state["entries"] = entries
    st.session_state["errors"] = errors

entries = st.session_state.get("entries", [])
errors = st.session_state.get("errors", [])

st.title("Propeptide Activation / Protease Explorer")
st.caption(
    "A protease that cuts *exactly* at the propeptide/mature junction cannot "
    "confirm activation — the mature terminal peptide looks identical with or "
    "without the propeptide. Green = a peptide spans the junction (diagnostic); "
    "red = the enzyme cuts through it (ambiguous)."
)

if errors:
    st.warning("Could not fetch: " + " · ".join(f"**{a}** ({m})" for a, m in errors))
if not entries:
    st.info("Pick a preset or type gene symbols, then press **Fetch / refresh**.")
    st.stop()
if not sel_proteases:
    st.warning("Select at least one protease in the sidebar.")
    st.stop()

# Split entries by whether they have a propeptide junction.
prop_entries = [e for e in entries if any(
    j.is_propeptide for j in analysis.find_junctions(e))]
noprop = [e for e in entries if e not in prop_entries]

tab_tracks, tab_rank, tab_silac, tab_nterm, tab_feat = st.tabs(
    ["Tracks", "Ranking", "SILAC", "N-terminomics", "Features"]
)

# --- Tracks (primary) ----------------------------------------------------------

with tab_tracks:
    c1, c2 = st.columns([3, 1])
    prop_only = c2.toggle("Propeptide junctions only", value=True)

    # Protease scorecard across the whole queried set.
    score = analysis.protease_scorecard(
        entries, sel_proteases, missed, min_len, max_len, prop_only)
    if score and score[0]["total_junctions"]:
        tot = score[0]["total_junctions"]
        st.markdown(f"#### Protease scorecard — {len(entries)} proteins, "
                    f"{tot} activation junctions")
        sdf = pd.DataFrame(score)
        sdf["diagnostic (of total)"] = sdf.apply(
            lambda r: f"{r['diagnostic']} / {r['total_junctions']} "
                      f"({r['pct_diagnostic']:.0f}%)", axis=1)
        st.dataframe(
            sdf[["protease", "diagnostic (of total)", "ambiguous",
                 "no_spanning_pep", "nterm_detectable"]],
            hide_index=True, width="stretch")
        st.bar_chart(sdf.set_index("protease")["diagnostic"], height=200)
        st.markdown("---")

    pool = prop_entries if prop_only else entries
    labels = {f"{e.gene or e.name} ({e.accession})": e for e in pool}
    picked = c1.multiselect("Proteins to display", list(labels.keys()),
                            default=list(labels.keys()))
    if noprop and prop_only:
        st.caption("No propeptide junction (hidden): "
                   + ", ".join(e.gene or e.name for e in noprop))
    for lab in picked:
        e = labels[lab]
        fig = tracks.protein_track_figure(
            e, sel_proteases, missed, min_len, max_len, cov_region, prop_only)
        st.pyplot(fig)
        # Per-protein best-protease line.
        allr, _ = analysis.rank_proteases(e, sel_proteases, missed, min_len, max_len)
        pj = [r for r in allr if r["is_propeptide_junction"]]
        if pj:
            best = max(pj, key=lambda r: r["score"])
            diag = [r["protease"] for r in pj
                    if r["status"].startswith("Diagnostic")]
            msg = (f"**{e.gene or e.name}** — best: **{best['protease']}** "
                   f"({best['status']}).")
            if diag:
                msg += f" Diagnostic proteases: {', '.join(sorted(set(diag)))}."
            else:
                msg += " ⚠️ No selected protease is diagnostic here."
            st.markdown(msg)
        st.markdown("---")

# --- Ranking table -------------------------------------------------------------

with tab_rank:
    only_pp = st.checkbox("Only propeptide-activation junctions", value=True)
    rows = []
    for e in entries:
        allr, _ = analysis.rank_proteases(e, sel_proteases, missed, min_len, max_len)
        rows.extend(allr)
    if only_pp:
        rows = [r for r in rows if r["is_propeptide_junction"]]
    if not rows:
        st.info("No matching junctions.")
    else:
        df = pd.DataFrame(rows)
        df["gene"] = df["name"].str.split("_").str[0]
        st.markdown("#### Best protease per junction")
        best = (df.sort_values("score", ascending=False)
                  .groupby(["accession", "junction", "side"], as_index=False).first())
        st.dataframe(best[["gene", "accession", "junction", "side", "adjacent",
                           "protease", "status", "score"]],
                     hide_index=True, width="stretch")
        st.markdown("#### All protease × junction combinations")
        st.dataframe(
            df[["gene", "junction", "side", "protease", "status", "score",
                "spanning_peptide", "spanning_len", "spanning_missed",
                "mature_term_peptide", "mature_term_len", "mature_term_detectable"]]
            .sort_values(["gene", "junction", "score"], ascending=[True, True, False]),
            hide_index=True, width="stretch")

# --- SILAC ---------------------------------------------------------------------

with tab_silac:
    st.markdown(
        "**Which amino acid is the optimal SILAC label?** For each protease, we "
        "count how many of its *diagnostic* peptides (the pro-form spanning peptide "
        "and the mature neo-terminal peptide) would carry ≥1 heavy residue — i.e. "
        "become quantifiable — for each candidate label. Only essential amino acids "
        "are tested (the rest can't be metabolically labelled to completion)."
    )
    # Matrix: protease x candidate amino acid -> % of diagnostic peptides quantifiable.
    matrix, best_rows = {}, []
    for p in sel_proteases:
        opt = analysis.silac_label_optimization(
            entries, p, missed, min_len, max_len, candidate_aa, prop_only)
        if not opt or not opt[0]["total_peptides"]:
            continue
        matrix[p] = {o["aa"]: o["pct_quantifiable"] for o in opt}
        top = opt[0]
        best_rows.append({
            "protease": p,
            "optimal label": top["amino_acid"],
            "quantifiable": f"{top['quantifiable']} / {top['total_peptides']} "
                            f"({top['pct_quantifiable']:.0f}%)",
            "avg labels/peptide": top["avg_labels_per_pep"],
        })
    if not matrix:
        st.info("No diagnostic peptides to optimise over.")
    else:
        st.markdown("#### Optimal label per protease")
        st.dataframe(pd.DataFrame(best_rows), hide_index=True, width="stretch")
        st.markdown("#### % of diagnostic peptides quantifiable — protease × label")
        mdf = pd.DataFrame(matrix).T  # rows = protease, cols = amino acid
        col_order = ["+".join(sorted(set(t.upper()))) for t in candidate_aa]
        mdf = mdf[[c for c in col_order if c in mdf.columns]]
        st.dataframe(mdf.style.background_gradient(cmap="Greens", axis=None,
                                                   vmin=0, vmax=100).format("{:.0f}"),
                     width="stretch")
        st.caption("Higher = more diagnostic peptides become SILAC-quantifiable "
                   "with that single label. Trypsin/Lys-C favour Lys; non-tryptic "
                   "enzymes often need Leu or another essential residue.")

# --- N-terminomics -------------------------------------------------------------

with tab_nterm:
    st.markdown(
        "N-terminomics (TAILS): the peptide recovered from a processed protein "
        "runs from the **mature N-terminus** to the first enzymatic cut site. "
        "Each protease track shows that captured peptide — **green if a detectable "
        "length, red if too short/long**. (C-terminal propeptides need "
        "C-terminomics and are not shown here.)"
    )
    # Only proteins with an N-terminal propeptide junction are meaningful here.
    npool = [e for e in (prop_entries if prop_only else entries)
             if any(j.side == "N" and (j.is_propeptide or not prop_only)
                    for j in analysis.find_junctions(e))]
    nlabels = {f"{e.gene or e.name} ({e.accession})": e for e in npool}
    if not nlabels:
        st.info("No N-terminal propeptide junctions in the current set.")
    else:
        npicked = st.multiselect("Proteins to display", list(nlabels.keys()),
                                 default=list(nlabels.keys()), key="nterm_pick")
        for lab in npicked:
            e = nlabels[lab]
            fig = tracks.protein_track_figure(
                e, sel_proteases, missed, min_len, max_len, cov_region,
                prop_only, mode="nterm")
            st.pyplot(fig)
        st.markdown("---")
        st.caption("Detail table:")
        rows = []
        for e in npool:
            for p in sel_proteases:
                rows.extend(analysis.nterminomics_capture(
                    e, p, missed, min_len, max_len))
        if rows:
            df = pd.DataFrame(rows)
            df["gene"] = df["name"].str.split("_").str[0]
            if prop_only:
                df = df[df["is_propeptide_junction"]]
            st.dataframe(df[["gene", "protease", "neo_nterm_pos", "adjacent",
                             "captured_peptide", "length", "detectable"]],
                         hide_index=True, width="stretch")

# --- Features overview ---------------------------------------------------------

with tab_feat:
    for e in entries:
        st.subheader(f"{e.gene or e.name} — {e.protein}  ·  {e.accession} "
                     f"({e.length} aa)")
        feat_rows = [
            {"feature": f.kind, "start": f.start, "end": f.end,
             "length": f.end - f.start + 1, "description": f.description}
            for f in e.features
        ]
        if feat_rows:
            st.dataframe(pd.DataFrame(feat_rows), hide_index=True,
                         width="stretch")
        else:
            st.write("No molecule-processing features annotated.")
        st.markdown("---")

# --- PDF export (sidebar) ------------------------------------------------------
# Placed last so `prop_only` (defined in the Tracks tab) is available.

st.sidebar.markdown("---")
if st.sidebar.button("Build PDF report"):
    with st.spinner("Building PDF report..."):
        st.session_state["pdf_bytes"] = report.build_pdf(
            entries, sel_proteases, missed, min_len, max_len, cov_region,
            candidate_aa, prop_only)
    st.session_state["pdf_name"] = (
        f"propeptide_report_{len(entries)}proteins.pdf")

if st.session_state.get("pdf_bytes"):
    st.sidebar.download_button(
        "⬇ Download PDF", data=st.session_state["pdf_bytes"],
        file_name=st.session_state.get("pdf_name", "propeptide_report.pdf"),
        mime="application/pdf")
    st.sidebar.caption("Report reflects the settings used when built. Rebuild "
                       "after changing parameters.")
