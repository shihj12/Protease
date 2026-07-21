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
- **Is the diagnostic peptide actually _usable_?** A peptide that spans the
  junction is only the first hurdle. The app also checks whether it is
  **proteome-unique** (not shared with a paralogue), whether it has been
  **observed by LC-MS** (PeptideAtlas), whether the protein has a **constitutive
  normalizer** peptide to anchor total abundance, and whether the peptide is a
  reasonable MS flyer (hydrophobicity / not membrane-buried). The **Confidence**
  tab and the *Actionable* scorecard column fold these together.
- **How trustworthy is the junction itself?** Each junction carries the UniProt
  **evidence level** for the propeptide boundary (experimental vs. by-similarity),
  and membrane topology (**transmembrane**, glycosylation, disulfides) is overlaid
  so you can spot ectodomain-shedding cases and poorly-detectable regions.
- **Will the peptide behave in the instrument?** Diagnostic peptides are flagged
  for **modification risks** that split or shift the signal — Met oxidation,
  N-terminal Gln/Glu → pyroGlu (common on neo-N-termini), the N-x-S/T
  glycosylation sequon (and annotated glyco sites), Cys, and the labile Asp-Pro
  bond.
- **Can a second enzyme rescue a hopeless junction?** The **double-digest search**
  (in the *Ranking* tab) tries every enzyme *pair*: co-digestion adds cut sites
  that can carve a single enzyme's over-long spanning peptide into the detectable
  window, or shorten an undetectable neo-N-terminus — surfacing junctions that
  only a two-enzyme digest can capture.
- **How good is my digest across the _whole_ proteome?** The **Proteome** tab
  scores a set of **orthogonal aliquots** (e.g. Trypsin in one tube, Glu-C in
  another) over all ~20k human proteins: % covered, % quantifiable with your SILAC
  label, % *uniquely* quantifiable, median sequence coverage, and the hydrophobic
  (poor-flyer) burden — with your queried proteins shown as a **priority set**
  beside the global numbers, and an **orthogonality-gain** table showing how much
  each added aliquot buys you.

## What counts as "diagnostic"

For each activation junction (mature-chain N-terminus preceded by a propeptide):

| Case | Meaning |
| --- | --- |
| **Diagnostic** | The enzyme does *not* cut at the junction, so a fully-specific peptide spans the pro/mature bond. Its presence proves the pro-form; its loss reports activation. |
| **Ambiguous** | The enzyme cuts *exactly* at the junction. The fully-cleaved mature N-terminal peptide looks the same with or without the propeptide. Only a missed-cleavage peptide could span it (unreliable). |
| **Mature N-term peptide** | The semi-specific peptide from the mature N-terminus to the first cut site. Its detection confirms the mature form and is the N-terminomics readout. |
| **Actionable** | Diagnostic **and** the spanning peptide is proteome-unique **and** the protein has ≥1 constitutive normalizer peptide (in the mature chain, unique, quantifiable) to divide by. Only actionable junctions give an experiment that can separate *activation* from a change in *expression*. |

## Confidence checks and data sources

Two local reference sets sharpen "detectable" into "actually works". Both are
optional and the app degrades gracefully if they are missing.

- **Proteome uniqueness** — the bundled Swiss-Prot human FASTA
  (`uniprotkb_human_*.fasta`) is indexed so every diagnostic/normalizer peptide
  can be checked for **proteome-uniqueness** (I and L are collapsed, because they
  are isobaric in MS). Paralogue families — cathepsins, proteasome β-subunits,
  arylsulfatases — routinely share peptides that then cannot be specifically
  quantified. *First launch builds a k-mer index (~20 s) and caches it to
  `.proteome_index.pkl`; later launches load in ~2 s.*
- **PeptideAtlas (observed by LC-MS)** — PeptideAtlas blocks automated web
  queries and asks users to use its bulk downloads, so this tool reads a
  **downloaded build file** instead of hitting the network. Drop a peptide list
  from <https://peptideatlas.org/builds/> (a Human build) into a `.peptideatlas/`
  folder next to `app.py` (any TSV/text file with a peptide-sequence column;
  an observation-count column is used if present). **Caveat:** PeptideAtlas
  builds are essentially fully-tryptic, so a *semi-specific* mature neo-N-terminal
  peptide will usually be absent even when perfectly detectable — the app marks
  those `n/a` rather than "not seen".

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

- `proteases.py` — cleavage-rule table + in-silico digestion, single and
  multi-enzyme co-digestion (`digest_multi`); extend `PROTEASES`.
- `uniprot.py` — UniProt REST fetch by gene symbol or accession, feature parsing
  (molecule processing + membrane topology + PTMs, with evidence codes), disk
  cache. Gene resolution prefers the entry whose *primary* gene name matches.
- `analysis.py` — junctions (N- and C-terminal), protease diagnostics, coverage,
  SILAC, N-terminomics, normalizer peptides, double-digest rescue search, and
  per-peptide confidence (uniqueness, PeptideAtlas, GRAVY/detectability,
  modification risks).
- `proteome.py` — indexes the bundled human FASTA for in-silico peptide
  uniqueness (stdlib only; disk-cached k-mer index, seed-intersection fast path).
- `coverage_stats.py` — whole-proteome coverage/quantifiability statistics for
  orthogonal aliquots (union at the protein level), with a priority subset and
  orthogonality gain; heavy scan cached to disk independent of the priority set.
- `peptideatlas.py` — offline lookup of LC-MS-observed peptides from a downloaded
  PeptideAtlas build; self-disables if no build file is present.
- `tracks.py` — the per-protein track figure (features + membrane topology + one
  track per protease), `mode='diagnostic'` or `mode='nterm'`.
- `presets.py` — curated gene lists (lysosomal hydrolases, proteasome subunits).
- `report.py` — multi-page PDF report (matplotlib PdfPages, no extra deps).
- `app.py` — Streamlit UI (Tracks / Ranking / Confidence / SILAC / N-terminomics
  / Features / Proteome).

## Orthogonal vs. co-digestion

Two different multi-enzyme ideas, kept separate:

- **Co-digestion** (Ranking → *Double-digest rescue*): both enzymes in **one
  tube**, so the cut sites are the *union* and you get new hybrid peptides. Useful
  to rescue a single junction.
- **Orthogonal digestion** (Proteome tab): each enzyme in a **separate aliquot**,
  measured independently; a protein counts if *any* aliquot delivers a qualifying
  peptide (union at the protein level). This is how you buy proteome-wide coverage.

The first Proteome computation scans all ~20k proteins (~20 s per enzyme) and
caches the result to `.coverage_cache/`; changing the priority set or SILAC label
re-aggregates instantly without re-scanning.

## Notes / assumptions

- Digestion is rule-based (no missed-modification chemistry); "detectable" = peptide
  length within the chosen window. It does not model m/z, charge, or ionisation.
  Hydrophobicity (GRAVY) is used only as a coarse flyer flag, not a retention-time
  predictor.
- SILAC quantifiability = peptide contains ≥1 labelled residue (K/R by default).
- Proteome-uniqueness collapses I/L (isobaric) and is only as complete as the
  bundled FASTA (reviewed human by default). A "unique" call means unique *within
  that set*.
- Junction evidence is UniProt's own ECO evidence for the propeptide/signal
  feature; "by similarity" boundaries are predictions, not proven cleavage sites.
- Proteases included: Trypsin, Trypsin/P, Lys-C, Arg-C, Glu-C (E and D/E),
  Chymotrypsin, Asp-N, Lys-N. Add more in `proteases.py`.
