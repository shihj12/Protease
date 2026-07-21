"""Offline PeptideAtlas lookup: has a peptide actually been observed by LC-MS?

PeptideAtlas blocks automated per-peptide web queries and asks users to use its
bulk downloads instead, so this module indexes a downloaded build file rather
than hitting the network. Drop a PeptideAtlas peptide list into this folder
(or a `.peptideatlas/` subfolder) and it is picked up automatically.

Where to get one:
    https://peptideatlas.org/builds/  ->  a Human build  ->  the peptide list
    (e.g. "APD_Hs_all.tsv" or the build's `*_peptides.tsv`). Any text/TSV file
    with a peptide-sequence column works; an observation-count column is used if
    present.

Interpretation caveat: PeptideAtlas builds are overwhelmingly **fully-tryptic,
fully-specific** peptides. So a "seen" result strongly supports a pro-form
spanning peptide or a constitutive normalizer, but a semi-specific mature
neo-N-terminal peptide will usually be *absent even if it is perfectly
detectable*. Callers should treat "not found" as informative only for
fully-specific peptides and as "n/a" for semi-specific ones.
"""

import glob
import os
import re
import threading

_HERE = os.path.dirname(os.path.abspath(__file__))

_AA = "ACDEFGHIKLMNPQRSTVWY"
_PEP_RE = re.compile(rf"^[{_AA}]{{5,}}$")
# Header names that likely hold the sequence / observation count.
_SEQ_COLS = ("peptide", "sequence", "peptide_sequence", "pep_seq", "peptideseq")
_OBS_COLS = ("n_observations", "n_obs", "nobs", "norm_psms", "n_psms", "observations",
             "psm_count", "spectra")

_lock = threading.Lock()
_loaded = False
_obs: dict[str, int] = {}     # peptide sequence -> observation count (0 if unknown)
_source: str | None = None


def find_file() -> str | None:
    """Locate a PeptideAtlas peptide list next to this module."""
    patterns = [
        os.path.join(_HERE, ".peptideatlas", "*"),
        os.path.join(_HERE, "*peptideatlas*.tsv"),
        os.path.join(_HERE, "*peptideatlas*.txt"),
        os.path.join(_HERE, "APD_*.tsv"),
        os.path.join(_HERE, "*_peptides.tsv"),
        os.path.join(_HERE, "peptideatlas.*"),
    ]
    for pat in patterns:
        hits = sorted(g for g in glob.glob(pat) if os.path.isfile(g))
        if hits:
            return hits[0]
    return None


def _pick_columns(header: list[str]):
    low = [h.strip().lower() for h in header]
    seq_i = next((i for i, h in enumerate(low) if h in _SEQ_COLS), None)
    obs_i = next((i for i, h in enumerate(low) if h in _OBS_COLS), None)
    return seq_i, obs_i


def _parse(path: str) -> dict[str, int]:
    obs: dict[str, int] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        first = fh.readline()
        delim = "\t" if "\t" in first else ("," if "," in first else None)
        cols = first.rstrip("\n").split(delim) if delim else first.split()
        seq_i, obs_i = _pick_columns(cols)
        header_consumed = seq_i is not None
        if not header_consumed:
            _index_line(first, obs, delim)  # first line was data, not a header
        for line in fh:
            if header_consumed:
                parts = line.rstrip("\n").split(delim) if delim else line.split()
                if seq_i is not None and seq_i < len(parts):
                    seq = parts[seq_i].strip().upper()
                    if _PEP_RE.match(seq):
                        n = _to_int(parts[obs_i]) if (obs_i is not None
                                                      and obs_i < len(parts)) else 0
                        obs[seq] = max(obs.get(seq, 0), n)
            else:
                _index_line(line, obs, delim)
    return obs


def _to_int(tok: str) -> int:
    try:
        return int(float(tok.strip()))
    except (ValueError, AttributeError):
        return 0


def _index_line(line: str, obs: dict, delim):
    """Header-less fallback: take the first peptide-looking token as the sequence
    and the first integer token as the observation count."""
    parts = line.rstrip("\n").split(delim) if delim else line.split()
    seq = next((p.strip().upper() for p in parts
                if _PEP_RE.match(p.strip().upper())), None)
    if seq is None:
        return
    n = next((_to_int(p) for p in parts if p.strip().isdigit()), 0)
    obs[seq] = max(obs.get(seq, 0), n)


def load(path: str | None = None) -> bool:
    """Index a PeptideAtlas build file. Idempotent; True if any peptides loaded."""
    global _loaded, _obs, _source
    with _lock:
        if _loaded:
            return bool(_obs)
        path = path or find_file()
        if not path or not os.path.exists(path):
            _loaded = True
            return False
        try:
            _obs = _parse(path)
            _source = path
        except OSError:
            _obs = {}
        _loaded = True
        return bool(_obs)


def available() -> bool:
    load()
    return bool(_obs)


def info() -> dict:
    load()
    return {"source": _source, "n_peptides": len(_obs)}


def observations(peptide: str):
    """Observation count for an exactly-matching peptide, or None if not present.

    None means "not in the loaded build" (or no build loaded — check available()).
    Matching is exact (case-insensitive); I/L are NOT collapsed because
    PeptideAtlas records the real observed sequence.
    """
    if not load():
        return None
    return _obs.get(peptide.strip().upper())


def seen(peptide: str):
    """True/False if a build is loaded, else None (unknown — no build present)."""
    if not load():
        return None
    return peptide.strip().upper() in _obs
