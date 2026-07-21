"""In-silico peptide uniqueness against the human proteome.

A diagnostic or normalizer peptide is only *specifically* quantifiable if its
sequence maps to exactly one protein. Paralogous families (the cathepsins,
proteasome subunits, arylsulfatases...) share peptides, so a peptide that looks
perfect on paper can be silently non-quantifiable because it is also produced by
a sibling gene.

This module indexes the bundled Swiss-Prot human FASTA and answers, for any
peptide, *how many distinct proteins contain it* (and which). It is pure/offline
and uses the stdlib only.

MS note: isoleucine and leucine are isobaric and generally indistinguishable in
a normal LC-MS/MS experiment, so both are collapsed to a single symbol before
comparison. Uniqueness is therefore evaluated the way the mass spectrometer
actually "sees" the sequence.
"""

import functools
import glob
import os
import pickle
import re
import threading

_HERE = os.path.dirname(os.path.abspath(__file__))

# A 5-mer inverted index turns the ~20k full-length substring scans per query
# into a scan of only the handful of proteins that share a seed with the peptide.
_KMER = 5

# Building the index is ~20s; caching it to disk makes subsequent launches ~4s.
# Bump the version to invalidate old caches after a format change.
_CACHE_VERSION = 2
_INDEX_CACHE = os.path.join(_HERE, ".proteome_index.pkl")

_lock = threading.Lock()
_loaded = False
_accs: list[str] = []
_genes: list[str] = []
_seqs: list[str] = []            # I/L-collapsed sequences, index-aligned with _accs
_kmer_index: dict[str, list[int]] = {}
_fasta_path: str | None = None


def collapse(seq: str) -> str:
    """Collapse the sequence the way MS sees it: I and L are isobaric."""
    return seq.upper().replace("I", "L").replace("J", "L")


def find_fasta() -> str | None:
    """Locate a bundled human proteome FASTA next to this module."""
    hits = sorted(glob.glob(os.path.join(_HERE, "*.fasta")))
    # Prefer a human/9606 file if several are present.
    for h in hits:
        base = os.path.basename(h).lower()
        if "human" in base or "9606" in base:
            return h
    return hits[0] if hits else None


def _parse_fasta(path: str):
    accs, genes, seqs = [], [], []
    acc, gene, buf = None, "", []

    def flush():
        if acc is not None:
            accs.append(acc)
            genes.append(gene)
            seqs.append(collapse("".join(buf)))

    gn_re = re.compile(r"\bGN=(\S+)")
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                buf = []
                # >sp|P07339|CATD_HUMAN ... GN=CTSD ...
                parts = line[1:].split("|")
                acc = parts[1] if len(parts) > 1 else line[1:].split()[0]
                m = gn_re.search(line)
                gene = m.group(1) if m else ""
            else:
                buf.append(line.strip())
        flush()
    return accs, genes, seqs


def _build_index(seqs) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for i, s in enumerate(seqs):
        # Unique 5-mers of this protein -> its index. Using a set avoids adding
        # the same protein many times for a repeated k-mer.
        for kmer in {s[j:j + _KMER] for j in range(len(s) - _KMER + 1)}:
            index.setdefault(kmer, []).append(i)
    return index


def _load_disk_cache(path: str):
    """Return cached (accs, genes, seqs, index) if fresh for `path`, else None."""
    if not os.path.exists(_INDEX_CACHE):
        return None
    try:
        if os.path.getmtime(_INDEX_CACHE) < os.path.getmtime(path):
            return None  # FASTA is newer than the cache
        with open(_INDEX_CACHE, "rb") as fh:
            blob = pickle.load(fh)
        if blob.get("version") != _CACHE_VERSION or blob.get("path") != path:
            return None
        return blob["accs"], blob["genes"], blob["seqs"], blob["index"]
    except (OSError, pickle.PickleError, KeyError, EOFError):
        return None


def _save_disk_cache(path, accs, genes, seqs, index):
    try:
        with open(_INDEX_CACHE, "wb") as fh:
            pickle.dump({"version": _CACHE_VERSION, "path": path, "accs": accs,
                         "genes": genes, "seqs": seqs, "index": index},
                        fh, protocol=5)
    except OSError:
        pass  # caching is best-effort; scanning still works without it


def load(path: str | None = None) -> bool:
    """Load and index the proteome. Idempotent; returns True if a corpus is ready.

    The k-mer index is cached to disk (`.proteome_index.pkl`), so the ~20s build
    happens only the first time (or after the FASTA changes).
    """
    global _loaded, _accs, _genes, _seqs, _kmer_index, _fasta_path
    with _lock:
        if _loaded:
            return bool(_seqs)
        path = path or find_fasta()
        if not path or not os.path.exists(path):
            _loaded = True
            return False
        cached = _load_disk_cache(path)
        if cached is not None:
            _accs, _genes, _seqs, _kmer_index = cached
        else:
            _accs, _genes, _seqs = _parse_fasta(path)
            _kmer_index = _build_index(_seqs)
            _save_disk_cache(path, _accs, _genes, _seqs, _kmer_index)
        _fasta_path = path
        _loaded = True
        return bool(_seqs)


def available() -> bool:
    load()
    return bool(_seqs)


def corpus_info() -> dict:
    load()
    return {"path": _fasta_path, "n_proteins": len(_seqs)}


@functools.lru_cache(maxsize=300_000)
def _seed_set(seed: str) -> frozenset:
    """Proteins containing this 5-mer, as a frozenset. Cached because the same
    seeds recur across millions of peptide queries (converting a large posting to
    a set is the dominant cost of uniqueness checks)."""
    idxs = _kmer_index.get(seed)
    return frozenset(idxs) if idxs else frozenset()


def _candidate_indices(pep: str):
    """Proteins that share every seed with the peptide are the only ones that can
    contain it. Returns None if the peptide is shorter than one k-mer (fall back
    to a full scan). Intersects rarest-seed-first for speed."""
    if len(pep) < _KMER:
        return None
    seeds = {pep[j:j + _KMER] for j in range(0, len(pep) - _KMER + 1, _KMER)}
    seeds.add(pep[len(pep) - _KMER:])  # cover the C-terminus
    sets = []
    for seed in seeds:
        ss = _seed_set(seed)
        if not ss:
            return set()  # a seed absent everywhere => peptide occurs nowhere
        sets.append(ss)
    sets.sort(key=len)          # start from the smallest posting
    cand = sets[0]
    for ss in sets[1:]:
        cand = cand & ss
        if not cand:
            return set()
    return cand


@functools.lru_cache(maxsize=100_000)
def _hits(pep: str) -> tuple:
    """(accessions, genes) of proteins containing the collapsed peptide `pep`."""
    cand = _candidate_indices(pep)
    idxs = range(len(_seqs)) if cand is None else cand
    hit_accs, hit_genes = [], []
    for i in idxs:
        if pep in _seqs[i]:
            hit_accs.append(_accs[i])
            hit_genes.append(_genes[i])
    return tuple(hit_accs), tuple(hit_genes)


def multiplicity(peptide: str, cap: int = 2) -> int:
    """Number of proteins containing `peptide` as a substring, counted up to `cap`.

    Fast path for uniqueness: if the peptide's k-mer seeds intersect to a single
    protein, that protein is necessarily its source and no other protein can
    contain it — so it is unique with NO substring verification. Only when several
    proteins share every seed do we verify substrings (early-exiting at `cap`).
    """
    if not load():
        return 0
    pep = collapse(peptide)
    if not pep:
        return 0
    cand = _candidate_indices(pep)
    if cand is None:  # peptide shorter than one k-mer: fall back to a full scan
        c = 0
        for s in _seqs:
            if pep in s:
                c += 1
                if c >= cap:
                    return c
        return c
    if len(cand) <= 1:
        return len(cand)  # only the source protein shares the seeds -> unique
    c = 0
    for i in cand:
        if pep in _seqs[i]:
            c += 1
            if c >= cap:
                return c
    return c


def protein_hits(peptide: str, limit: int = 6) -> dict:
    """How many distinct human proteins contain `peptide` (I/L-collapsed).

    Returns {n_proteins, accessions, genes, available}. `available` is False when
    no FASTA is bundled (the caller should then treat uniqueness as unknown, not
    as non-unique).
    """
    if not load():
        return {"n_proteins": None, "accessions": [], "genes": [], "available": False}
    pep = collapse(peptide)
    if not pep:
        return {"n_proteins": 0, "accessions": [], "genes": [], "available": True}
    hit_accs, hit_genes = _hits(pep)
    return {
        "n_proteins": len(hit_accs),
        "accessions": list(hit_accs[:limit]),
        "genes": list(hit_genes[:limit]),
        "available": True,
    }


def is_unique(peptide: str) -> bool | None:
    """True/False if the corpus is loaded, else None (unknown — no FASTA).

    Uses the fast seed-intersection path (multiplicity with cap=2)."""
    if not load():
        return None
    return multiplicity(peptide, cap=2) == 1
