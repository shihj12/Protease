"""Domain analysis: activation junctions, protease diagnostics, coverage,
SILAC labelling, and N-terminomics capture.

All functions are pure over (Entry, parameters) and return plain dicts/lists so
the Streamlit layer can turn them into DataFrames without domain logic leaking
into the UI.

Coordinates everywhere are 1-based inclusive, matching UniProt and proteases.py.
"""

from dataclasses import dataclass

import peptideatlas
import proteome
from proteases import (PROTEASES, Protease, Peptide, cleavage_sites,
                       cleavage_sites_multi, digest, digest_multi)

# Default SILAC scheme: heavy Lys (+8.0142) and Arg (+10.0083). Only the residue
# set matters for "is this peptide quantifiable"; masses are carried for display.
DEFAULT_SILAC = {"K": 8.0142, "R": 10.0083}

# Amino acids that can actually be used as metabolic (SILAC) labels: the essential
# amino acids (cells cannot synthesise them de novo, so incorporation is complete),
# plus Arg. Non-essential residues are excluded because in-cell synthesis dilutes
# the label. These are the only candidates worth "optimising over".
LABELABLE_AA = ["R", "K", "L", "I", "V", "F", "M", "H", "T", "W"]
AA_NAMES = {
    "R": "Arg", "K": "Lys", "L": "Leu", "I": "Ile", "V": "Val", "F": "Phe",
    "M": "Met", "H": "His", "T": "Thr", "W": "Trp", "Y": "Tyr",
}


# --- Peptide biochemistry & detectability confidence ---------------------------

# Kyte-Doolittle hydropathy. GRAVY (the mean over a peptide) predicts LC-MS
# behaviour: strongly positive = hydrophobic, poorly recovered / membrane-buried;
# strongly negative = very hydrophilic, elutes in the void with little retention.
_KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# GRAVY thresholds beyond which a tryptic peptide is a risky MS flyer.
GRAVY_HYDROPHOBIC = 1.0
GRAVY_HYDROPHILIC = -2.0


def gravy(seq: str) -> float:
    """Grand average of hydropathy over the peptide (Kyte-Doolittle)."""
    seq = seq.upper()
    vals = [_KD[a] for a in seq if a in _KD]
    return sum(vals) / len(vals) if vals else 0.0


def flyer_flag(seq: str) -> str:
    """Coarse detectability call from hydrophobicity: 'ok' / 'hydrophobic' /
    'hydrophilic'. Hydrophobic peptides are the ones that go missing."""
    g = gravy(seq)
    if g > GRAVY_HYDROPHOBIC:
        return "hydrophobic"
    if g < GRAVY_HYDROPHILIC:
        return "hydrophilic"
    return "ok"


import re as _re

_SEQUON_RE = _re.compile(r"N[^P][ST]")  # N-glycosylation consensus N-x-S/T (x != P)


def mod_risks(seq: str, start: int = None, end: int = None, entry=None) -> list[str]:
    """Sequence-level modifications that split or shift a peptide's signal and
    undermine clean quantification. Returns short human labels.

    Flags: Met oxidation; N-terminal Gln/Glu -> pyroGlu cyclisation (common on
    neo-N-termini); the N-x-S/T N-glycosylation sequon; Cys (needs alkylation);
    the labile Asp-Pro bond. If `entry`/coords are given, an *annotated*
    glycosylation site inside the peptide is flagged explicitly (stronger than
    the sequon heuristic).
    """
    seq = seq.upper()
    risks: list[str] = []
    if "M" in seq:
        risks.append("Met-ox")
    if seq[:1] in ("Q", "E"):
        risks.append(f"pyroGlu (N-term {seq[0]})")
    if _SEQUON_RE.search(seq):
        risks.append("N-glyco sequon")
    if "C" in seq:
        risks.append("Cys (alkylation)")
    if "DP" in seq:
        risks.append("Asp-Pro (labile)")
    if entry is not None and start is not None and end is not None:
        for f in entry.features:
            if f.kind == "CARBOHYD" and not (end < f.start or start > f.end):
                risks.append("glyco site (annotated)")
                break
    return risks


def overlaps_tm(entry, start: int, end: int) -> bool:
    """Does [start, end] overlap a transmembrane (or intramembrane) span?

    Peptides buried in the membrane are hard to solubilise and typically absent
    from bottom-up datasets — worth flagging for membrane proteins whose
    "activation" is really ectodomain shedding."""
    for f in entry.features:
        if f.kind in ("TRANSMEM", "INTRAMEM") and not (end < f.start or start > f.end):
            return True
    return False


def annotate_peptide(seq, start, end, entry, label_residues=("K", "R"),
                     specific=True) -> dict:
    """Confidence annotations for one diagnostic/normalizer peptide.

    Combines proteome uniqueness (is it specific to this protein?), PeptideAtlas
    observation (has it ever been seen by LC-MS?), SILAC quantifiability, and
    hydrophobicity/membrane detectability. `specific=False` marks a semi-specific
    peptide (a mature neo-terminus), for which a PeptideAtlas miss is uninformative
    because those builds are essentially fully-tryptic.
    """
    hits = proteome.protein_hits(seq)
    unique = None if not hits["available"] else (hits["n_proteins"] == 1)
    other_genes = [g for g in hits["genes"] if g and g] if hits["available"] else []
    obs = peptideatlas.observations(seq) if specific else None
    return {
        "unique": unique,                       # True / False / None (no FASTA)
        "n_proteins": hits["n_proteins"],
        "shared_with": [g for g in other_genes][:5],
        "quantifiable": silac_count(seq, label_residues) >= 1,
        "gravy": round(gravy(seq), 2),
        "flyer": flyer_flag(seq),
        "in_tm": overlaps_tm(entry, start, end),
        "mod_risks": mod_risks(seq, start, end, entry),
        "observed": obs,                        # int n_obs / None (miss or n/a)
        "atlas_applicable": specific,
    }


# --- Activation junctions ------------------------------------------------------


@dataclass
class Junction:
    pos: int           # J: the processing bond is between residue J and J+1
    side: str          # 'N' = mature lies to the right; 'C' = mature to the left
    mature_start: int  # first residue of the mature chain
    mature_end: int    # last residue of the mature chain
    chain_desc: str
    adjacent: str      # feature on the propeptide side: PROPEP / SIGNAL / TRANSIT / none
    adj_desc: str
    adj_evidence: str = ""       # evidence label for the propeptide-side feature
    adj_evidence_rank: int = -1  # 3=experimental .. 0=automatic; -1=unstated

    @property
    def is_propeptide(self) -> bool:
        return self.adjacent == "PROPEP"


def find_junctions(entry) -> list[Junction]:
    """Identify processing junctions where a mature chain abuts a propeptide (or
    signal/transit peptide).

    Two geometries are detected:
      * N-terminal removal  — a chain whose start-1 is the end of a propeptide/
        signal/transit peptide (mature lies to the right; neo-N-terminus at start).
      * C-terminal removal  — a chain whose end+1 is the start of a propeptide
        (mature lies to the left; e.g. legumain, proteasome beta subunits).

    Junctions are flagged so the caller knows which are true propeptide-activation
    events (`is_propeptide`) and which side the mature chain is on.
    """
    features = entry.features
    n = entry.length
    by_key: dict = {}

    def consider(cand: Junction):
        key = (cand.pos, cand.side)
        cur = by_key.get(key)
        if (cur is None
                or (cand.is_propeptide and not cur.is_propeptide)
                or (cand.is_propeptide == cur.is_propeptide
                    and (cand.mature_end - cand.mature_start)
                        > (cur.mature_end - cur.mature_start))):
            by_key[key] = cand

    for chain in entry.of_kind("CHAIN"):
        # N-terminal processing bond (chain begins after another feature).
        if chain.start > 1:
            j = chain.start - 1
            adj, desc, ev, evr = "none", "", "", -1
            for f in features:
                if f.kind in ("PROPEP", "SIGNAL", "TRANSIT") and f.end == j:
                    adj, desc, ev, evr = f.kind, f.description, f.evidence, f.evidence_rank
                    if f.kind == "PROPEP":
                        break
            # Skip bare initiator-methionine removal (chain starts at residue 2 with
            # nothing annotated before it) — not an activation event.
            if not (j == 1 and adj == "none"):
                consider(Junction(j, "N", chain.start, chain.end,
                                  chain.description, adj, desc, ev, evr))
        # C-terminal processing bond (a propeptide begins right after the chain).
        if chain.end < n:
            j = chain.end
            for f in features:
                if f.kind == "PROPEP" and f.start == j + 1:
                    consider(Junction(j, "C", chain.start, chain.end,
                                      chain.description, "PROPEP", f.description,
                                      f.evidence, f.evidence_rank))
                    break

    return [by_key[k] for k in sorted(by_key)]


# --- Protease diagnostics ------------------------------------------------------


def _mature_terminus_peptide(seq: str, sites, jn: "Junction"):
    """The semi-specific peptide at the mature neo-terminus generated by removal
    of the propeptide, running to the first enzymatic cut site inside the chain.

    N-side: mature N-terminus (J+1) -> first cut site downstream.
    C-side: last cut site upstream -> mature C-terminus (J).
    `sites` is a sorted list of cut positions. Returns (start, end, seq).
    """
    j = jn.pos
    if jn.side == "N":
        start = j + 1
        downstream = [s for s in sites if s > j]
        end = downstream[0] if downstream else len(seq)
    else:  # 'C' — mature ends at J
        end = j
        upstream = [s for s in sites if s < j]
        start = (upstream[-1] + 1) if upstream else 1
    return start, end, seq[start - 1:end]


def junction_diagnostics(entry, protease, missed_cleavages=2, min_len=7, max_len=40,
                         label_residues=("K", "R")):
    """For each junction, decide whether the protease can report activation.

    Returns a list of row dicts (one per junction), including a status and a
    numeric score for ranking proteases, plus confidence annotations (uniqueness,
    PeptideAtlas observation, hydrophobicity) on the diagnostic peptides.
    """
    if isinstance(protease, str):
        protease = PROTEASES[protease]
    seq = entry.sequence
    peptides = digest(seq, protease, missed_cleavages, min_len, max_len)
    sites_list = cleavage_sites(seq, protease)
    return _diagnose_rows(entry, protease.name, seq, peptides, sites_list,
                          min_len, max_len, label_residues)


def junction_diagnostics_combo(entry, protease_list, missed_cleavages=2, min_len=7,
                               max_len=40, label_residues=("K", "R")):
    """Junction diagnostics for a co-digestion (union of several enzymes' cuts).

    A combined digest can produce a *shorter* fully-specific peptide that spans a
    junction where a single enzyme's spanning peptide fell outside the detectable
    window — the basis of the double-digest rescue search."""
    names = [p if isinstance(p, str) else p.name for p in protease_list]
    seq = entry.sequence
    peptides = digest_multi(seq, protease_list, missed_cleavages, min_len, max_len)
    sites_list = cleavage_sites_multi(seq, protease_list)
    return _diagnose_rows(entry, " + ".join(names), seq, peptides, sites_list,
                          min_len, max_len, label_residues)


def _diagnose_rows(entry, label, seq, peptides, sites_list, min_len, max_len,
                   label_residues):
    """Core per-junction diagnostic logic shared by the single- and multi-protease
    paths. `label` names the enzyme(s); `sites_list` is the (union of) cut sites."""
    sites = set(sites_list)
    rows = []
    for jn in find_junctions(entry):
        j = jn.pos
        # Peptides spanning the pro/mature bond and within the detectable window.
        spanning = [p for p in peptides if p.spans_bond(j)]
        spanning.sort(key=lambda p: (p.missed, p.length))
        best = spanning[0] if spanning else None
        cuts_at_junction = j in sites

        # Mature neo-terminal (semi-specific) peptide, on the correct side.
        m_start, m_end, m_seq = _mature_terminus_peptide(seq, sites_list, jn)
        m_len = m_end - m_start + 1
        m_detectable = min_len <= m_len <= max_len

        # A spanning peptide is only *reliably* diagnostic when the junction is
        # NOT an enzymatic cut site. If the enzyme cuts exactly at the junction,
        # the fully-cleaved mature terminal peptide is identical whether or not
        # the propeptide was ever present (the ambiguity case) — the only spanning
        # peptide then requires a missed cleavage at the junction and is
        # unreliable, so we down-rank it rather than call it diagnostic.
        best_bonus = 1 if m_detectable else 0
        if cuts_at_junction:
            if best is not None:
                status = "Ambiguous - spans only via a missed-cleavage peptide"
                score = 1 + best_bonus
            else:
                status = "Ambiguous (cuts exactly at junction)"
                score = 0 + best_bonus
        else:
            if best is not None:
                status = "Diagnostic (fully-specific pro-form peptide)"
                score = 3 + best_bonus
            else:
                status = "No detectable spanning peptide"
                score = 1 + best_bonus

        # Confidence annotations. The spanning peptide is fully-specific (PeptideAtlas
        # applies); the mature terminal peptide is semi-specific (it does not).
        span_ann = (annotate_peptide(best.seq, best.start, best.end, entry,
                                     label_residues) if best else None)
        term_ann = (annotate_peptide(m_seq, m_start, m_end, entry, label_residues,
                                     specific=False) if m_seq else None)

        row = {
            "accession": entry.accession,
            "name": entry.name,
            "protease": label,
            "junction": f"{j}|{j + 1}",
            "side": jn.side,
            "mature_start": jn.mature_start,
            "adjacent": jn.adjacent,
            "adj_evidence": jn.adj_evidence,
            "adj_evidence_rank": jn.adj_evidence_rank,
            "is_propeptide_junction": jn.is_propeptide,
            "status": status,
            "score": score,
            "cuts_at_junction": cuts_at_junction,
            "spanning_peptide": best.seq if best else "",
            "spanning_start": best.start if best else None,
            "spanning_end": best.end if best else None,
            "spanning_len": best.length if best else None,
            "spanning_missed": best.missed if best else None,
            "spanning_unique": span_ann["unique"] if span_ann else None,
            "spanning_observed": span_ann["observed"] if span_ann else None,
            "spanning_gravy": span_ann["gravy"] if span_ann else None,
            "spanning_flyer": span_ann["flyer"] if span_ann else None,
            "spanning_mod_risks": span_ann["mod_risks"] if span_ann else [],
            "spanning_shared_with": span_ann["shared_with"] if span_ann else [],
            "mature_term_peptide": m_seq,
            "mature_term_len": m_len,
            "mature_term_detectable": m_detectable,
            "mature_term_unique": term_ann["unique"] if term_ann else None,
            "mature_term_gravy": term_ann["gravy"] if term_ann else None,
            "mature_term_flyer": term_ann["flyer"] if term_ann else None,
            "mature_term_in_tm": term_ann["in_tm"] if term_ann else None,
            "mature_term_mod_risks": term_ann["mod_risks"] if term_ann else [],
        }
        rows.append(row)
    return rows


def rank_proteases(entry, protease_names, missed_cleavages=2, min_len=7, max_len=40):
    """Best score per (junction, protease). Returns rows for all combos, and a
    per-junction 'best protease' summary."""
    all_rows = []
    for name in protease_names:
        all_rows.extend(
            junction_diagnostics(entry, name, missed_cleavages, min_len, max_len)
        )
    # Best protease per junction.
    best_by_junction = {}
    for r in all_rows:
        key = (r["junction"], r["side"])
        cur = best_by_junction.get(key)
        if cur is None or r["score"] > cur["score"]:
            best_by_junction[key] = r
    return all_rows, list(best_by_junction.values())


# --- Double-digest (multi-protease) rescue -------------------------------------


def double_digest_search(entry, protease_names, missed_cleavages=2, min_len=7,
                         max_len=40, label_residues=("K", "R"),
                         propeptide_only=True, only_improving=True):
    """Find enzyme *pairs* (co-digestion) that beat every single enzyme at a
    junction.

    Adding a second protease's cuts can carve a single enzyme's over-long spanning
    peptide down into the detectable window (a "spanning rescue"), or shorten an
    undetectable mature neo-N-terminal peptide into range (an "N-term rescue").
    For each junction this compares the best single protease against the best pair
    and returns the pairs that strictly improve on it.
    """
    from itertools import combinations

    def best_per_junction(row_iter):
        best = {}
        for d in row_iter:
            if propeptide_only and not d["is_propeptide_junction"]:
                continue
            key = (d["junction"], d["side"])
            cur = best.get(key)
            # Higher score wins; ties break toward the shorter spanning peptide.
            if (cur is None or d["score"] > cur["score"]
                    or (d["score"] == cur["score"]
                        and (d["spanning_len"] or 10**9) < (cur["spanning_len"] or 10**9))):
                best[key] = d
        return best

    singles = []
    for name in protease_names:
        singles += junction_diagnostics(entry, name, missed_cleavages, min_len,
                                        max_len, label_residues)
    best_single = best_per_junction(singles)

    combos = []
    for pair in combinations(sorted(protease_names), 2):
        combos += junction_diagnostics_combo(entry, list(pair), missed_cleavages,
                                             min_len, max_len, label_residues)
    best_combo = best_per_junction(combos)

    results = []
    for key, single in best_single.items():
        combo = best_combo.get(key)
        if combo is None:
            continue
        single_diag = single["status"].startswith("Diagnostic")
        combo_diag = combo["status"].startswith("Diagnostic")
        span_rescue = (combo_diag and not single_diag
                       and combo["spanning_unique"] is not False)
        nterm_rescue = (combo["side"] == "N" and combo["mature_term_detectable"]
                        and not single["mature_term_detectable"])
        improves = span_rescue or nterm_rescue
        if only_improving and not improves:
            continue
        results.append({
            "accession": entry.accession, "gene": entry.gene or entry.name,
            "junction": key[0], "side": key[1],
            "single_best": single["protease"], "single_status": single["status"],
            "combo": combo["protease"], "combo_status": combo["status"],
            "combo_spanning_peptide": combo["spanning_peptide"],
            "combo_spanning_len": combo["spanning_len"],
            "combo_spanning_unique": combo["spanning_unique"],
            "combo_nterm_peptide": combo["mature_term_peptide"],
            "combo_nterm_len": combo["mature_term_len"],
            "rescue": "spanning" if span_rescue else "neo-N-term",
            "improves": improves,
        })
    return results


# --- Coverage ------------------------------------------------------------------


def coverage(entry, protease, missed_cleavages=2, min_len=7, max_len=40, region="full"):
    """Sequence coverage by detectable peptides.

    region: 'full' = whole precursor; 'chain' = restrict denominator to mature
    chain residues only. Returns dict with percent and per-residue covered flags.
    """
    if isinstance(protease, str):
        protease = PROTEASES[protease]
    seq = entry.sequence
    n = len(seq)
    peptides = digest(seq, protease, missed_cleavages, min_len, max_len)

    covered = bytearray(n)  # 0/1 per residue (0-based)
    for p in peptides:
        for i in range(p.start - 1, p.end):
            covered[i] = 1

    if region == "chain":
        denom_positions = set()
        for ch in entry.of_kind("CHAIN"):
            denom_positions.update(range(ch.start - 1, ch.end))
        if not denom_positions:
            denom_positions = set(range(n))
    else:
        denom_positions = set(range(n))

    denom = len(denom_positions)
    hit = sum(1 for i in denom_positions if covered[i])
    pct = 100.0 * hit / denom if denom else 0.0
    return {
        "accession": entry.accession,
        "name": entry.name,
        "protease": protease.name,
        "percent": pct,
        "covered_residues": hit,
        "total_residues": denom,
        "covered_flags": bytes(covered),   # length n, 1 where covered
        "cut_sites": cleavage_sites(seq, protease),
        "n_peptides": len(peptides),
    }


# --- SILAC ---------------------------------------------------------------------


def silac_count(pep_seq: str, label_residues=("K", "R")) -> int:
    return sum(pep_seq.upper().count(r) for r in label_residues)


def silac_annotate(peptides, label_residues=("K", "R"), silac_masses=None):
    """Annotate a list of Peptide objects (or (seq,) strings) with label counts."""
    silac_masses = silac_masses or DEFAULT_SILAC
    rows = []
    for p in peptides:
        seq = p.seq if isinstance(p, Peptide) else str(p)
        n_label = silac_count(seq, label_residues)
        heavy_shift = sum(
            seq.upper().count(r) * silac_masses.get(r, 0.0) for r in label_residues
        )
        row = {
            "peptide": seq,
            "length": len(seq),
            "label_sites": n_label,
            "quantifiable": n_label >= 1,
            "heavy_delta_Da": round(heavy_shift, 4),
        }
        if isinstance(p, Peptide):
            row["start"] = p.start
            row["end"] = p.end
            row["missed"] = p.missed
        rows.append(row)
    return rows


def diagnostic_silac(entry, protease, missed_cleavages, min_len, max_len,
                     label_residues=("K", "R"), silac_masses=None):
    """SILAC quantifiability of the diagnostic peptides for each junction:
    the pro-form spanning peptide and the mature neo-N-terminal peptide."""
    diag = junction_diagnostics(entry, protease, missed_cleavages, min_len, max_len)
    rows = []
    for d in diag:
        term_label = "mature N-term" if d["side"] == "N" else "mature C-term"
        for kind, seq in (("pro-form spanning", d["spanning_peptide"]),
                          (term_label, d["mature_term_peptide"])):
            if not seq:
                continue
            n_label = silac_count(seq, label_residues)
            rows.append({
                "accession": d["accession"],
                "name": d["name"],
                "protease": d["protease"],
                "junction": d["junction"],
                "peptide_role": kind,
                "peptide": seq,
                "length": len(seq),
                "label_sites": n_label,
                "quantifiable": n_label >= 1,
            })
    return rows


# --- N-terminomics -------------------------------------------------------------


def nterminomics_capture(entry, protease, missed_cleavages=2, min_len=7, max_len=40):
    """Which mature neo-N-termini (from propeptide/signal removal) would be
    captured as detectable N-terminal peptides.

    In an N-terminomics workflow (e.g. TAILS / reductive dimethylation), protein
    N-termini are blocked and internal tryptic peptides depleted; the peptide you
    recover from a processed protein runs from the mature N-terminus to the first
    enzymatic cut site. Its detection reports the activation event directly.
    """
    if isinstance(protease, str):
        protease = PROTEASES[protease]
    seq = entry.sequence
    sites_list = cleavage_sites(seq, protease)
    rows = []
    for jn in find_junctions(entry):
        if jn.side != "N":
            continue  # N-terminomics captures N-termini only (C-side needs C-TAILS)
        m_start, m_end, m_seq = _mature_terminus_peptide(seq, sites_list, jn)
        m_len = m_end - m_start + 1
        rows.append({
            "accession": entry.accession,
            "name": entry.name,
            "protease": protease.name,
            "neo_nterm_pos": jn.mature_start,
            "adjacent": jn.adjacent,
            "is_propeptide_junction": jn.is_propeptide,
            "captured_peptide": m_seq,
            "start": m_start,
            "end": m_end,
            "length": m_len,
            "detectable": min_len <= m_len <= max_len,
        })
    return rows


# --- Normalizer (reference) peptides -------------------------------------------


def _mature_positions(entry) -> set:
    """1-based residue positions that belong to a mature chain and are NOT part of
    any propeptide/signal/transit peptide — i.e. present in both the pro- and the
    activated form. A peptide lying entirely here is a *constitutive* readout of
    total protein, independent of activation state.
    """
    pos = set()
    for ch in entry.of_kind("CHAIN"):
        pos.update(range(ch.start, ch.end + 1))
    if not pos:  # no chain annotated: fall back to the whole precursor
        pos = set(range(1, entry.length + 1))
    for f in entry.features:
        if f.kind in ("PROPEP", "SIGNAL", "TRANSIT", "INIT_MET"):
            pos.difference_update(range(f.start, f.end + 1))
    return pos


def normalizer_peptides(entry, protease, missed_cleavages=2, min_len=7, max_len=40,
                        label_residues=("K", "R"), require_unique=True,
                        require_observed=False):
    """Constitutive reference peptides for normalising the activation signal.

    Idea #2: a diagnostic peptide is worthless unless *another* quantifiable
    peptide on the same protein anchors total abundance — otherwise you cannot
    tell activation from a change in expression. A valid normalizer is:
      * fully-cleaved (no missed cleavage) and a detectable length,
      * entirely within the mature chain (present in pro- and mature forms),
      * SILAC-quantifiable (carries a label residue),
      * proteome-unique (`require_unique`) so it measures *this* protein, and
      * optionally PeptideAtlas-observed (`require_observed`).

    Returns rows sorted best-first (observed & well-behaved peptides first).
    """
    if isinstance(protease, str):
        protease = PROTEASES[protease]
    seq = entry.sequence
    mature = _mature_positions(entry)
    peptides = digest(seq, protease, missed_cleavages, min_len, max_len)

    out = []
    for p in peptides:
        if p.missed != 0:
            continue  # a reference peptide should be robustly fully-cleaved
        if not (mature.issuperset(range(p.start, p.end + 1))):
            continue  # must lie wholly in the constitutive mature region
        if silac_count(p.seq, label_residues) < 1:
            continue
        ann = annotate_peptide(p.seq, p.start, p.end, entry, label_residues)
        if require_unique and ann["unique"] is False:
            continue
        if require_observed and not ann["observed"]:
            continue
        out.append({
            "accession": entry.accession, "protease": protease.name,
            "peptide": p.seq, "start": p.start, "end": p.end, "length": p.length,
            "unique": ann["unique"], "observed": ann["observed"],
            "gravy": ann["gravy"], "flyer": ann["flyer"], "in_tm": ann["in_tm"],
            "mod_risks": ann["mod_risks"],
        })
    # Prefer observed, non-hydrophobic, non-membrane peptides.
    out.sort(key=lambda r: (
        r["observed"] or 0, r["flyer"] == "ok", not r["in_tm"]), reverse=True)
    return out


def normalizer_summary(entry, protease, missed_cleavages=2, min_len=7, max_len=40,
                       label_residues=("K", "R"), require_unique=True,
                       require_observed=False) -> dict:
    """Whether a protein/protease pair yields at least one valid normalizer."""
    ns = normalizer_peptides(entry, protease, missed_cleavages, min_len, max_len,
                             label_residues, require_unique, require_observed)
    name = protease if isinstance(protease, str) else protease.name
    return {
        "accession": entry.accession, "name": entry.name, "protease": name,
        "n_normalizers": len(ns), "has_normalizer": bool(ns),
        "best_normalizer": ns[0]["peptide"] if ns else "",
    }


# --- Aggregate protease scorecard ---------------------------------------------


def protease_scorecard(entries, protease_names, missed_cleavages=2,
                       min_len=7, max_len=40, propeptide_only=True,
                       label_residues=("K", "R")):
    """Aggregate result across the whole queried set: for each protease, how many
    activation junctions it can diagnose ("Trypsin: X of Y diagnostic").

    Counts a junction once per protein/junction (not per protease-row). Returns a
    list of dict rows sorted best-first.
    """
    # Total junctions in scope (constant across proteases).
    scope = []
    for e in entries:
        for jn in find_junctions(e):
            if propeptide_only and not jn.is_propeptide:
                continue
            scope.append((e.accession, f"{jn.pos}|{jn.pos + 1}", jn.side))
    total = len(scope)

    rows = []
    for name in protease_names:
        diagnostic = ambiguous = no_span = nterm_ok = actionable = 0
        # Normalizer availability is per protein/protease, so compute once per entry.
        has_norm = {e.accession: normalizer_summary(
            e, name, missed_cleavages, min_len, max_len,
            label_residues)["has_normalizer"] for e in entries}
        for e in entries:
            for d in junction_diagnostics(e, name, missed_cleavages, min_len,
                                          max_len, label_residues):
                if propeptide_only and not d["is_propeptide_junction"]:
                    continue
                is_diag = d["status"].startswith("Diagnostic")
                if is_diag:
                    diagnostic += 1
                    # "Actionable" only if the diagnostic peptide is proteome-unique
                    # AND the protein has a constitutive normalizer to divide by.
                    if d["spanning_unique"] is not False and has_norm[e.accession]:
                        actionable += 1
                elif d["status"].startswith("Ambiguous"):
                    ambiguous += 1
                else:
                    no_span += 1
                if d["side"] == "N" and d["mature_term_detectable"]:
                    nterm_ok += 1
        rows.append({
            "protease": name,
            "diagnostic": diagnostic,
            "actionable": actionable,
            "ambiguous": ambiguous,
            "no_spanning_pep": no_span,
            "total_junctions": total,
            "pct_diagnostic": round(100.0 * diagnostic / total, 1) if total else 0.0,
            "pct_actionable": round(100.0 * actionable / total, 1) if total else 0.0,
            "nterm_detectable": nterm_ok,
        })
    rows.sort(key=lambda r: (r["actionable"], r["diagnostic"], r["nterm_detectable"]),
              reverse=True)
    return rows


# --- SILAC label optimisation --------------------------------------------------


def diagnostic_peptides(entry, protease, missed_cleavages=2, min_len=7, max_len=40,
                        propeptide_only=True):
    """The peptides that report activation at each junction: the pro-form spanning
    peptide (present in the pro-form) and the mature neo-terminal peptide (present
    after activation). Returns list of dict rows."""
    out = []
    for d in junction_diagnostics(entry, protease, missed_cleavages, min_len, max_len):
        if propeptide_only and not d["is_propeptide_junction"]:
            continue
        term_role = "mature N-term" if d["side"] == "N" else "mature C-term"
        if d["spanning_peptide"]:
            out.append({"accession": entry.accession, "junction": d["junction"],
                        "role": "pro-form spanning", "peptide": d["spanning_peptide"]})
        if d["mature_term_peptide"]:
            out.append({"accession": entry.accession, "junction": d["junction"],
                        "role": term_role, "peptide": d["mature_term_peptide"]})
    return out


def _label_display(token: str) -> str:
    aas = list(token.upper())
    return "+".join(AA_NAMES.get(a, a) for a in aas)


def silac_label_optimization(entries, protease, missed_cleavages=2, min_len=7,
                             max_len=40, candidates=None, propeptide_only=True):
    """Which label best covers a protease's diagnostic peptides.

    Each candidate is one or more labelable residues (e.g. "K", "L", or "KR" for
    the standard Lys+Arg SILAC pair — a peptide counts if it contains ANY of them).
    Reports how many diagnostic peptides across the queried set would carry >=1
    heavy residue (become SILAC-quantifiable). Ranked best-first; the top row is
    the optimal label. A single internal residue can out-score K or R alone,
    because a tryptic peptide ends in only one of K/R — which is exactly why the
    standard "KR" pair exists.
    """
    candidates = candidates or (["KR"] + LABELABLE_AA)
    peps = []
    for e in entries:
        peps.extend(diagnostic_peptides(e, protease, missed_cleavages, min_len,
                                        max_len, propeptide_only))
    total = len(peps)
    rows = []
    for token in candidates:
        residues = set(token.upper())
        counts = [sum(p["peptide"].upper().count(a) for a in residues) for p in peps]
        quant = sum(1 for c in counts if c >= 1)
        aa_key = "+".join(sorted(residues))
        rows.append({
            "amino_acid": f"{aa_key} ({_label_display(token)})",
            "aa": aa_key,
            "quantifiable": quant,
            "total_peptides": total,
            "pct_quantifiable": round(100.0 * quant / total, 1) if total else 0.0,
            "avg_labels_per_pep": round(sum(counts) / total, 2) if total else 0.0,
            "min_labels": min(counts) if counts else 0,
        })
    rows.sort(key=lambda r: (r["quantifiable"], r["avg_labels_per_pep"]), reverse=True)
    return rows
