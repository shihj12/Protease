"""Propeptide Activation / Protease Explorer — Streamlit app.

Run:  streamlit run app.py

Enter gene symbols (CTSD) or UniProt accessions, pick proteases, and read a
track-per-protease view of which enzymes can actually capture propeptide
activation — plus coverage, SILAC quantifiability, and N-terminomics capture.
"""

import pandas as pd
import streamlit as st

import analysis
import coverage_stats
import peptideatlas
import proteome
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

# Primary SILAC label used for quantifiability / normalizer checks (first
# candidate; the standard KR pair by default).
primary_label = tuple(candidate_aa[0]) if candidate_aa else ("K", "R")

st.sidebar.markdown("---")
st.sidebar.caption("**Confidence filters**")
require_unique = st.sidebar.checkbox(
    "Require proteome-unique peptides", value=True,
    help="Flag diagnostic/normalizer peptides that also occur in another human "
         "protein (paralogues) — they can't specifically quantify this one.")
require_observed = st.sidebar.checkbox(
    "Require PeptideAtlas-observed normalizer", value=False,
    help="Only count a normalizer if it appears in the loaded PeptideAtlas build.")

# --- Data-source status (proteome FASTA + optional PeptideAtlas build) ---------
st.sidebar.markdown("---")
with st.sidebar.expander("Data sources", expanded=False):
    if proteome.available():
        ci = proteome.corpus_info()
        st.caption(f"✓ Proteome uniqueness: {ci['n_proteins']:,} human proteins")
    else:
        st.caption("✗ No proteome FASTA found — uniqueness disabled. Add a "
                   "human `*.fasta` next to app.py.")
    if peptideatlas.available():
        pi = peptideatlas.info()
        st.caption(f"✓ PeptideAtlas: {pi['n_peptides']:,} observed peptides")
    else:
        st.caption("PeptideAtlas not loaded (optional). Drop a build file into a "
                   "`.peptideatlas/` folder to enable observation checks. See "
                   "peptideatlas.org/builds/.")

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

(tab_tracks, tab_rank, tab_conf, tab_silac, tab_nterm, tab_feat,
 tab_proteome) = st.tabs(
    ["Tracks", "Ranking", "Confidence", "SILAC", "N-terminomics", "Features",
     "Proteome"]
)

# --- Tracks (primary) ----------------------------------------------------------

with tab_tracks:
    c1, c2 = st.columns([3, 1])
    prop_only = c2.toggle("Propeptide junctions only", value=True)

    # Protease scorecard across the whole queried set.
    score = analysis.protease_scorecard(
        entries, sel_proteases, missed, min_len, max_len, prop_only, primary_label)
    if score and score[0]["total_junctions"]:
        tot = score[0]["total_junctions"]
        st.markdown(f"#### Protease scorecard — {len(entries)} proteins, "
                    f"{tot} activation junctions")
        sdf = pd.DataFrame(score)
        sdf["diagnostic (of total)"] = sdf.apply(
            lambda r: f"{r['diagnostic']} / {r['total_junctions']} "
                      f"({r['pct_diagnostic']:.0f}%)", axis=1)
        sdf["actionable (of total)"] = sdf.apply(
            lambda r: f"{r['actionable']} / {r['total_junctions']} "
                      f"({r['pct_actionable']:.0f}%)", axis=1)
        st.dataframe(
            sdf[["protease", "diagnostic (of total)", "actionable (of total)",
                 "ambiguous", "no_spanning_pep", "nterm_detectable"]],
            hide_index=True, width="stretch",
            column_config={"actionable (of total)": st.column_config.TextColumn(
                "actionable (of total)",
                help="Diagnostic AND the peptide is proteome-unique AND the protein "
                     "has a constitutive normalizer peptide — i.e. actually usable.")})
        st.bar_chart(sdf.set_index("protease")[["diagnostic", "actionable"]],
                     height=200)
        st.caption("**Actionable** = diagnostic peptide is proteome-unique *and* "
                   "the protein has ≥1 constitutive normalizer peptide to anchor "
                   "total abundance. A diagnostic peptide with no normalizer can't "
                   "separate activation from expression change.")
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

        # --- Double-digest rescue -------------------------------------------
        st.markdown("---")
        st.markdown("#### Double-digest rescue")
        st.caption(
            "Junctions where **no single enzyme** yields a usable peptide but a "
            "**co-digestion (two enzymes together)** does — the extra cut sites "
            "carve an over-long spanning peptide into the detectable window, or "
            "shorten an undetectable mature neo-N-terminus into range.")
        dd_all = st.checkbox("Show all pairs (not just improvements)", value=False,
                             key="dd_all")
        dd_rows = []
        for e in entries:
            dd_rows.extend(analysis.double_digest_search(
                e, sel_proteases, missed, min_len, max_len, primary_label,
                propeptide_only=only_pp, only_improving=not dd_all))
        if not dd_rows:
            st.info("No junction needs a double digest with the current proteases "
                    "and length window — every rescuable junction is already "
                    "covered by a single enzyme (or none can be rescued by a pair).")
        else:
            ddf = pd.DataFrame(dd_rows)
            ddf["combo unique"] = ddf["combo_spanning_unique"].map(
                {True: "✓", False: "✗ shared"}).fillna("—")
            show = ddf[["gene", "junction", "side", "rescue", "single_best",
                        "single_status", "combo", "combo_status",
                        "combo_spanning_peptide", "combo_spanning_len", "combo unique"]]
            st.dataframe(show.rename(columns={
                "single_best": "best single", "single_status": "single result",
                "combo": "enzyme pair", "combo_status": "pair result",
                "combo_spanning_peptide": "pair peptide",
                "combo_spanning_len": "len"}),
                hide_index=True, width="stretch")

# --- Confidence (uniqueness / observed / normalizer / evidence) ----------------

with tab_conf:
    st.markdown(
        "**Is each diagnostic peptide actually usable?** A green bar in *Tracks* "
        "only means a peptide spans the junction. For a real experiment it must "
        "also be **proteome-unique** (not shared with a paralogue), ideally "
        "**observed** in LC-MS (PeptideAtlas), and the protein needs a "
        "**constitutive normalizer** peptide to anchor total abundance. The "
        "junction annotation's **evidence** tells you whether the boundary itself "
        "is experimentally proven or only inferred by similarity."
    )
    only_actionable = st.checkbox("Show only fully-actionable rows", value=False)

    def _fmt_unique(u):
        return "✓ unique" if u is True else ("✗ shared" if u is False else "?")

    def _fmt_obs(o, applicable=True):
        if not peptideatlas.available():
            return "n/a (no build)"
        if not applicable:
            return "n/a (semi-tryptic)"
        return str(o) if o else "not seen"

    conf_rows = []
    for e in entries:
        norm_cache = {p: analysis.normalizer_summary(
            e, p, missed, min_len, max_len, primary_label,
            require_unique, require_observed) for p in sel_proteases}
        for p in sel_proteases:
            norm = norm_cache[p]
            for d in analysis.junction_diagnostics(
                    e, p, missed, min_len, max_len, primary_label):
                if not d["is_propeptide_junction"]:
                    continue
                diag = d["status"].startswith("Diagnostic")
                unique = d["spanning_unique"]
                actionable = (diag and unique is not False
                              and norm["has_normalizer"])
                gaps = []
                if not diag:
                    gaps.append("not diagnostic")
                if unique is False:
                    gaps.append(f"shared with {', '.join(d['spanning_shared_with'][:3])}")
                if not norm["has_normalizer"]:
                    gaps.append("no normalizer")
                if only_actionable and not actionable:
                    continue
                conf_rows.append({
                    "gene": e.gene or e.name,
                    "junction": d["junction"],
                    "side": d["side"],
                    "junction evidence": d["adj_evidence"] or "—",
                    "protease": p,
                    "diagnostic": "✓" if diag else "✗",
                    "spanning peptide": d["spanning_peptide"] or "—",
                    "unique": _fmt_unique(unique),
                    "observed (PA)": _fmt_obs(d["spanning_observed"]),
                    "GRAVY": d["spanning_gravy"],
                    "flyer": d["spanning_flyer"] or "—",
                    "mod-risk": ", ".join(d["spanning_mod_risks"]) or "—",
                    "normalizers": norm["n_normalizers"],
                    "best normalizer": norm["best_normalizer"] or "—",
                    "actionable": "✅" if actionable else "—",
                    "gaps": "; ".join(gaps),
                })
    if not conf_rows:
        st.info("No propeptide-activation junctions in the current set "
                "(or none pass the filter).")
    else:
        cdf = pd.DataFrame(conf_rows).sort_values(
            ["gene", "junction", "actionable"], ascending=[True, True, False])
        st.dataframe(cdf, hide_index=True, width="stretch")
        st.caption("**unique** = maps to exactly one human protein (I/L collapsed). "
                   "**observed (PA)** = PeptideAtlas observation count; semi-tryptic "
                   "mature-terminus peptides are marked n/a because those builds are "
                   "fully-tryptic. **mod-risk** = residues that split/shift the signal "
                   "(Met-ox, N-term pyroGlu, N-glyco sequon, Cys, labile Asp-Pro). "
                   "**junction evidence** = UniProt evidence for the propeptide "
                   "boundary (experimental ▸ curator inference ▸ by similarity ▸ "
                   "automatic).")

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
             "length": f.end - f.start + 1, "evidence": f.evidence or "—",
             "description": f.description}
            for f in e.features
        ]
        if feat_rows:
            st.dataframe(pd.DataFrame(feat_rows), hide_index=True,
                         width="stretch")
        else:
            st.write("No molecule-processing features annotated.")
        st.markdown("---")

# --- Proteome-wide coverage ----------------------------------------------------

with tab_proteome:
    st.markdown(
        "**Whole-proteome coverage for orthogonal digestions.** Each selected "
        "enzyme is a **separate aliquot** (Trypsin in one, Glu-C in another…); a "
        "protein counts if *any* aliquot delivers a qualifying peptide (union at "
        "the protein level — not co-digestion). Numbers are shown for fully-cleaved "
        "peptides (headline) and with missed cleavages (secondary), over the whole "
        "human proteome and — side by side — your **queried proteins** as the "
        "priority set."
    )
    if not coverage_stats.available():
        st.warning("No proteome FASTA found next to app.py — this view needs the "
                   "bundled human `*.fasta`.")
    else:
        aliquots = st.multiselect(
            "Orthogonal aliquots (enzymes, each a separate digest)",
            protease_names(), default=sel_proteases,
            help="Union at the protein level. Order matters only for the "
                 "orthogonality-gain table below.")
        pa_label = st.text_input(
            "SILAC label residues", value="".join(primary_label),
            help="A peptide is 'quantifiable' if it contains any of these residues.")
        plabel = tuple(c for c in pa_label.upper() if c.isalpha()) or ("K", "R")

        atlas_ok = peptideatlas.available()
        observed_gate = st.checkbox(
            "Realistic mode — require PeptideAtlas-observed peptides",
            value=False, disabled=not atlas_ok,
            help="Only count peptides a mass spectrometer has actually seen. This "
                 "collapses the theoretical ceiling toward reality (it folds in "
                 "flyability and practical abundance).")
        if not atlas_ok:
            st.caption("⚠️ Realistic mode needs a PeptideAtlas build — drop one into "
                       "`.peptideatlas/` (see peptideatlas.org/builds). Without it, "
                       "numbers are the **theoretical ceiling** (a labelable peptide "
                       "merely *exists*), which massively overstates real coverage.")
        elif observed_gate:
            tryptic = {"Trypsin", "Trypsin/P", "Lys-C"}
            non_tryptic = [a for a in aliquots if a not in tryptic]
            if non_tryptic:
                st.caption("⚠️ PeptideAtlas builds are essentially fully-tryptic, so "
                           f"observed counts for {', '.join(non_tryptic)} will be "
                           "near-zero — the gate is only meaningful for tryptic "
                           "aliquots (Trypsin / Lys-C).")
        st.caption(f"Using detectable length {min_len}–{max_len} aa and up to "
                   f"{missed} missed cleavages from the sidebar. Priority set = the "
                   f"{len(entries)} queried protein(s).")

        run = st.button("Compute proteome statistics", type="primary")
        if run and aliquots:
            bar = st.progress(0.0, text="Digesting the proteome…")

            def _cb(done, total):
                bar.progress(done / total, text=f"Scanning proteins {done:,}/{total:,}")

            summary = coverage_stats.proteome_summary(
                aliquots, plabel, min_len, max_len, missed,
                priority_accs=[e.accession for e in entries],
                require_observed=observed_gate, progress=_cb)
            bar.empty()
            st.session_state["proteome_summary"] = summary
        elif run:
            st.warning("Select at least one aliquot enzyme.")

        summary = st.session_state.get("proteome_summary")
        if summary:
            mode = ("empirical (PeptideAtlas-observed)" if summary.get("observed_gate")
                    else "theoretical ceiling")
            st.markdown(f"#### {' + '.join(summary['enzymes'])}  ·  label "
                        f"{summary['label']}  ·  {summary['n_proteins']:,} proteins")
            st.caption(f"Mode: **{mode}**. " + (
                "These are peptides actually seen by LC-MS." if summary.get("observed_gate")
                else "A labelable peptide merely *exists* — an upper bound, not a "
                     "prediction. Enable *Realistic mode* (needs a PeptideAtlas build) "
                     "for an empirically grounded estimate."))

            def _metric_rows(block):
                h, m = block["0"], block["M"]
                order = [("Coverage (%)", "coverage_pct"),
                         ("Quantifiable (%)", "quant_pct"),
                         ("Uniquely quantifiable (%)", "unique_pct"),
                         ("Median seq. coverage (%)", "median_cov"),
                         ("Hydrophobic peptides (%)", "hydrophobic_pct")]
                return {lab: {"fully-cleaved": h[k], "+ missed": m[k]}
                        for lab, k in order}

            c1, c2 = st.columns(2)
            with c1:
                st.markdown("**Whole proteome**")
                st.dataframe(pd.DataFrame(_metric_rows(summary["whole"])).T,
                             width="stretch")
            with c2:
                st.markdown(f"**Priority set** ({summary['n_priority']} proteins)")
                if summary["priority"]:
                    st.dataframe(pd.DataFrame(_metric_rows(summary["priority"])).T,
                                 width="stretch")
                else:
                    st.info("No queried proteins matched the proteome FASTA.")

            st.markdown("#### Orthogonality gain")
            st.caption("Uniquely-quantifiable coverage (fully-cleaved) as each "
                       "aliquot is added — how much a second/third enzyme buys you.")
            odf = pd.DataFrame(summary["orthogonality"])
            cols = ["aliquots", "added", "whole_unique_pct"]
            if odf["priority_unique_pct"].notna().any():
                cols.append("priority_unique_pct")
            st.dataframe(odf[cols].rename(columns={
                "aliquots": "aliquot set", "added": "+ enzyme",
                "whole_unique_pct": "whole uniquely-quant (%)",
                "priority_unique_pct": "priority uniquely-quant (%)"}),
                hide_index=True, width="stretch")
            st.caption("Headline = fully-cleaved (0 missed) peptides only; "
                       "**+ missed** allows up to the sidebar's missed-cleavage "
                       "setting. 'Quantifiable' needs a label residue; 'uniquely "
                       "quantifiable' also needs the peptide to be proteome-unique "
                       "(substring, I/L collapsed).")
        else:
            st.info("Choose your aliquots and press **Compute proteome statistics**. "
                    "First run builds a per-enzyme index (~20 s/enzyme) and caches "
                    "it; later runs are instant.")

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
