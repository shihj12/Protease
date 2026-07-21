# Propeptide Activation / Protease Explorer

A lightweight personal tool for planning proteomics experiments that detect
**propeptide activation** of proteins such as lysosomal hydrolases.

Many hydrolases are made as an inactive **pro**-form and activated by proteolytic
removal of a **propeptide**, exposing a mature neo-N-terminus. This app helps you
answer, before you run anything at the bench:

- **Which protease can actually capture the activation event?** If your digestion
  enzyme cuts *exactly* at the propeptide/mature junction, the resulting mature
  N-terminal peptide is identical whether or not the propeptide was ever removed —
  so that protease **cannot confirm activation**. The app flags this.
- **What is my sequence coverage** across many proteins and proteases at once?
- **Are the diagnostic peptides SILAC-quantifiable?** (default Lys+8 / Arg+10).
- **N-terminomics:** which mature neo-N-termini would be captured as detectable
  N-terminal peptides.

## What counts as "diagnostic"

For each activation junction (mature-chain N-terminus preceded by a propeptide):

| Case | Meaning |
| --- | --- |
| **Diagnostic** | The enzyme does *not* cut at the junction, so a fully-specific peptide spans the pro/mature bond. Its presence proves the pro-form; its loss reports activation. |
| **Ambiguous** | The enzyme cuts *exactly* at the junction. The fully-cleaved mature N-terminal peptide looks the same with or without the propeptide. Only a missed-cleavage peptide could span it (unreliable). |
| **Mature N-term peptide** | The semi-specific peptide from the mature N-terminus to the first cut site. Its detection confirms the mature form and is the N-terminomics readout. |

## How to use it

1. In the sidebar, click a **preset** (Lysosomal hydrolases / Proteasome subunits)
   or type your own **gene symbols** (`CTSD CTSB`) and/or UniProt accessions.
2. Choose proteases and parameters, press **Fetch / refresh**.
3. Read the **Tracks** tab: each protein gets one figure with a features track on
   top and **one track per protease** below. Green bar = a peptide spans the
   junction (diagnostic); red bar = the enzyme cuts through it (ambiguous). The
   ✓/✗ next to each protease name summarises it at a glance.

Handles both **N-terminal** propeptides (e.g. cathepsins) and **C-terminal**
propeptides (e.g. legumain, proteasome β-subunits).

### Export

Press **Build PDF report** in the sidebar to compile a multi-page PDF (cover +
parameters, protease scorecard, SILAC label optimisation, and every diagnostic /
N-terminomics track figure), then **Download PDF**. Rebuild after changing
parameters — the report reflects the settings used when built.

## Setup

Streamlit wheels may lag on very new Python. Use a Python **3.11 or 3.12** env:

```bash
conda create -n propep python=3.12
conda activate propep
pip install -r requirements.txt
```

## Run

```bash
streamlit run app.py
```

Then in the sidebar: paste UniProt accessions (defaults are cathepsin D `P07339`,
glucocerebrosidase `P04062`, cathepsin B `P07858`, acid alpha-glucosidase
`P10253`), pick proteases and parameters, and press **Fetch / refresh**.

Sequences **and** the signal/propeptide/chain feature annotations are pulled live
from the UniProt REST API and cached under `.uniprot_cache/` for offline reuse.

## Files

- `proteases.py` — cleavage-rule table + in-silico digestion (extend `PROTEASES`).
- `uniprot.py` — UniProt REST fetch by gene symbol or accession, feature parsing,
  disk cache. Gene resolution prefers the entry whose *primary* gene name matches.
- `analysis.py` — junctions (N- and C-terminal), protease diagnostics, coverage,
  SILAC, N-terminomics.
- `tracks.py` — the per-protein track figure (features + one track per protease),
  `mode='diagnostic'` or `mode='nterm'`.
- `presets.py` — curated gene lists (lysosomal hydrolases, proteasome subunits).
- `report.py` — multi-page PDF report (matplotlib PdfPages, no extra deps).
- `app.py` — Streamlit UI (Tracks / Ranking / SILAC / N-terminomics / Features).

## Notes / assumptions

- Digestion is rule-based (no missed-modification chemistry); "detectable" = peptide
  length within the chosen window. It does not model m/z, charge, or ionisation.
- SILAC quantifiability = peptide contains ≥1 labelled residue (K/R by default).
- Proteases included: Trypsin, Trypsin/P, Lys-C, Arg-C, Glu-C (E and D/E),
  Chymotrypsin, Asp-N, Lys-N. Add more in `proteases.py`.
