"""Whole-proteome digestion statistics for experiment design.

Answers "how good is this enzyme (or set of orthogonal aliquots) across the
*entire* human proteome" alongside the targeted per-protein view:

  * coverage        - % of proteins with >=1 detectable peptide
  * quantifiable    - % with >=1 detectable, SILAC-label-bearing peptide
  * unique quant.   - % with >=1 proteotypic (substring-unique) label peptide
  * median coverage - per-protein % sequence covered
  * hydrophobic     - fraction of quantifiable peptides that are poor flyers

Orthogonal digestion (e.g. Trypsin in one aliquot, Glu-C in another) is modelled
as the **union at the protein level**: a protein counts if *any* aliquot delivers
a qualifying peptide. This is distinct from co-digestion (one tube, union of cut
sites) handled by analysis.double_digest_search.

Every statistic is reported twice: a headline number on fully-cleaved (0-missed)
peptides and a secondary number that also allows missed cleavages.

Performance: uniqueness is the expensive part, so it is resolved *lazily* — a
protein needs only one unique peptide to count, so we stop checking as soon as we
find it (most proteins hit one immediately). The final summary is cached to disk
keyed by all parameters.
"""

import hashlib
import os
import pickle
import re

import proteases
import proteome

_HERE = os.path.dirname(os.path.abspath(__file__))
_CACHE_DIR = os.path.join(_HERE, ".coverage_cache")
_CACHE_VERSION = 2

# Kyte-Doolittle hydropathy as a byte-indexed table (fast GRAVY over bytes).
_KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
       "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
       "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
_KD_TBL = [0.0] * 256
for _a, _v in _KD.items():
    _KD_TBL[ord(_a)] = _v

GRAVY_HYDROPHOBIC = 1.0

_orig = None  # [(acc, gene, seq)] in FASTA order, original (un-collapsed) sequences


def _load_original():
    """Parse the bundled FASTA once, keeping original (un-collapsed) sequences.

    proteome.py stores I/L-collapsed sequences for uniqueness; digestion, GRAVY
    and label counting need the real residues, so we parse our own copy here.
    """
    global _orig
    if _orig is not None:
        return _orig
    path = proteome.find_fasta()
    if not path or not os.path.exists(path):
        _orig = []
        return _orig
    out = []
    acc, gene, buf = None, "", []
    gn_re = re.compile(r"\bGN=(\S+)")

    def flush():
        if acc is not None:
            out.append((acc, gene, "".join(buf).upper()))

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                buf = []
                parts = line[1:].split("|")
                acc = parts[1] if len(parts) > 1 else line[1:].split()[0]
                m = gn_re.search(line)
                gene = m.group(1) if m else ""
            else:
                buf.append(line.strip())
        flush()
    _orig = out
    return _orig


def available() -> bool:
    return bool(_load_original())


def corpus_info() -> dict:
    return {"n_proteins": len(_load_original())}


def _is_hydrophobic(seq: str) -> bool:
    b = seq.encode()
    return (sum(_KD_TBL[c] for c in b) / len(b)) > GRAVY_HYDROPHOBIC if b else False


# --- Summary (cached) ---------------------------------------------------------


def _median(values) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    n = len(s)
    mid = n // 2
    return s[mid] if n % 2 else (s[mid - 1] + s[mid]) / 2.0


def _blank():
    return {"detect": 0, "quant": 0, "unique": 0, "cov": [],
            "quant_pep": 0, "hyd_pep": 0}


def _finalise(acc, n):
    return {
        "n": n,
        "coverage_pct": round(100.0 * acc["detect"] / n, 1) if n else 0.0,
        "quant_pct": round(100.0 * acc["quant"] / n, 1) if n else 0.0,
        "unique_pct": round(100.0 * acc["unique"] / n, 1) if n else 0.0,
        "median_cov": round(_median(acc["cov"]), 1),
        "hydrophobic_pct": (round(100.0 * acc["hyd_pep"] / acc["quant_pep"], 1)
                            if acc["quant_pep"] else 0.0),
    }


def _scan_key(order, label_residues, min_len, max_len, missed):
    payload = repr([
        _CACHE_VERSION, list(order), sorted(r.upper() for r in label_residues),
        min_len, max_len, missed,
    ])
    return hashlib.md5(payload.encode()).hexdigest()[:16]


def _scan(order, label_residues, min_len, max_len, missed, progress=None,
          use_cache=True):
    """Per-protein records for one aliquot set — the heavy pass, cached to disk.

    Independent of any priority set, so changing your queried proteins does NOT
    force a recompute. Each record is:
        (accession, level0_tuple, levelM_tuple, per_enzyme_unique0)
    where each level tuple is (detect, quant, unique, cov_pct, quant_pep, hyd_pep).
    """
    key = _scan_key(order, label_residues, min_len, max_len, missed)
    path = os.path.join(_CACHE_DIR, f"scan_{key}.pkl")
    if use_cache and os.path.exists(path):
        try:
            with open(path, "rb") as fh:
                return pickle.load(fh)
        except (OSError, pickle.PickleError, EOFError):
            pass

    corpus = _load_original()
    proteome.load()  # ensure the uniqueness index is ready
    n = len(corpus)
    label_set = frozenset(r.upper() for r in label_residues)
    records = []
    for i, (acc, _gene, seq) in enumerate(corpus):
        length = len(seq)
        peps_by_enz = [proteases.digest(seq, name, missed, min_len, max_len)
                       for name in order]
        levels = {}
        enz_unique0 = [False] * len(order)
        for level in ("0", "M"):
            covered = bytearray(length)
            detect = quant = unique = False
            qp = hp = 0
            for ei, peps in enumerate(peps_by_enz):
                for p in peps:
                    if level == "0" and p.missed != 0:
                        continue
                    detect = True
                    covered[p.start - 1:p.end] = b"\x01" * (p.end - p.start + 1)
                    if not label_set.isdisjoint(p.seq):
                        quant = True
                        qp += 1
                        if _is_hydrophobic(p.seq):
                            hp += 1
                        need = (not unique) or (level == "0" and not enz_unique0[ei])
                        if need and proteome.is_unique(p.seq):
                            unique = True
                            if level == "0":
                                enz_unique0[ei] = True
            cov_pct = 100.0 * sum(covered) / length if length else 0.0
            levels[level] = (int(detect), int(quant), int(unique), cov_pct, qp, hp)
        records.append((acc, levels["0"], levels["M"], tuple(enz_unique0)))
        if progress is not None and (i % 500 == 0 or i == n - 1):
            progress(i + 1, n)

    os.makedirs(_CACHE_DIR, exist_ok=True)
    try:
        with open(path, "wb") as fh:
            pickle.dump(records, fh, protocol=5)
    except OSError:
        pass
    return records


def _aggregate(records):
    g = {"0": _blank(), "M": _blank()}
    for r in records:
        for li, level in enumerate(("0", "M")):
            detect, quant, unique, cov, qp, hp = r[1 + li]
            a = g[level]
            a["detect"] += detect
            a["quant"] += quant
            a["unique"] += unique
            a["cov"].append(cov)
            a["quant_pep"] += qp
            a["hyd_pep"] += hp
    return g


def _cumulative_unique(records, n_enz):
    """Cumulative count of uniquely-quantifiable proteins (0-missed) as aliquots
    are added in order."""
    cum = [0] * n_enz
    for r in records:
        eu = r[3]
        seen = False
        for k in range(n_enz):
            seen = seen or eu[k]
            if seen:
                cum[k] += 1
    return cum


def proteome_summary(enzyme_names, label_residues=("K", "R"), min_len=7, max_len=40,
                     missed=2, priority_accs=None, progress=None, use_cache=True):
    """Whole-proteome (and priority-subset) statistics for orthogonal aliquots.

    The heavy per-protein scan is cached independently of the priority set, so
    aggregation for the whole proteome and any queried subset is instant once the
    scan exists.
    """
    order = list(enzyme_names)
    records = _scan(order, label_residues, min_len, max_len, missed, progress,
                    use_cache)
    n = len(records)
    prio_set = {a.upper() for a in (priority_accs or [])}
    prio_records = [r for r in records if r[0].upper() in prio_set]
    n_prio = len(prio_records)

    whole = _aggregate(records)
    prio = _aggregate(prio_records) if n_prio else None
    cum_whole = _cumulative_unique(records, len(order))
    cum_prio = _cumulative_unique(prio_records, len(order)) if n_prio else None

    ortho = []
    for k in range(len(order)):
        ortho.append({
            "aliquots": " + ".join(order[:k + 1]),
            "added": order[k],
            "whole_unique_pct": round(100.0 * cum_whole[k] / n, 1) if n else 0.0,
            "priority_unique_pct": (round(100.0 * cum_prio[k] / n_prio, 1)
                                    if n_prio else None),
        })

    return {
        "n_proteins": n,
        "n_priority": n_prio,
        "enzymes": order,
        "label": "+".join(sorted(r.upper() for r in label_residues)),
        "whole": {"0": _finalise(whole["0"], n), "M": _finalise(whole["M"], n)},
        "priority": ({"0": _finalise(prio["0"], n_prio),
                      "M": _finalise(prio["M"], n_prio)} if n_prio else None),
        "orthogonality": ortho,
    }
