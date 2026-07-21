"""Fetch protein sequence + molecule-processing features from UniProt.

Uses the stdlib only (urllib). Results are cached to disk so repeated runs and
re-parameterising the app do not re-hit the network. Molecule-processing
features (signal peptide, transit peptide, propeptide, chain, peptide) carry the
coordinates we need to locate activation junctions.
"""

import json
import os
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, asdict, field

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".uniprot_cache")
BASE_URL = "https://rest.uniprot.org/uniprotkb/{acc}.json"
SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

# UniProt primary-accession pattern (used to tell a gene symbol from an accession).
import re as _re
_ACC_RE = _re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)


def looks_like_accession(token: str) -> bool:
    return bool(_ACC_RE.match(token.strip().upper()))

# Bump when the cached feature schema changes so stale caches are refreshed
# (adds membrane/PTM feature kinds and per-feature evidence).
SCHEMA_VERSION = 2

# UniProt feature "type" strings we care about. Molecule processing drives the
# activation-junction logic; the topology/PTM group provides "alternative
# stability" context (transmembrane spans, disulfides, glycosylation).
FEATURE_TYPES = {
    # molecule processing
    "Signal": "SIGNAL",
    "Transit peptide": "TRANSIT",
    "Propeptide": "PROPEP",
    "Chain": "CHAIN",
    "Peptide": "PEPTIDE",
    "Initiator methionine": "INIT_MET",
    # topology / membrane
    "Transmembrane": "TRANSMEM",
    "Intramembrane": "INTRAMEM",
    "Topological domain": "TOPO_DOM",
    # stability-relevant PTMs / bonds
    "Disulfide bond": "DISULFID",
    "Glycosylation": "CARBOHYD",
}

# ECO evidence codes -> (human label, confidence rank). Higher rank = stronger
# evidence that the annotation (e.g. a propeptide boundary) is real rather than
# predicted by homology. Used to flag which activation junctions are actually
# experimentally supported.
_ECO = {
    "ECO:0000269": ("experimental", 3),
    "ECO:0007744": ("experimental (combinatorial)", 3),
    "ECO:0000305": ("curator inference", 2),
    "ECO:0000303": ("author statement", 2),
    "ECO:0000250": ("by similarity", 1),
    "ECO:0000255": ("sequence model", 1),
    "ECO:0000256": ("automatic annotation", 0),
    "ECO:0000312": ("imported", 0),
    "ECO:0000313": ("imported", 0),
}


def _summarise_evidence(evidences: list):
    """Reduce a feature's evidence list to (label, rank); ('', -1) if none."""
    best = ("", -1)
    for ev in evidences or []:
        label, rank = _ECO.get(ev.get("evidenceCode", ""), ("other", 0))
        if rank > best[1]:
            best = (label, rank)
    return best


@dataclass
class Feature:
    kind: str        # normalised: SIGNAL / TRANSIT / PROPEP / CHAIN / PEPTIDE / ...
    start: int       # 1-based inclusive
    end: int         # 1-based inclusive
    description: str
    evidence: str = ""       # human label for the strongest supporting evidence
    evidence_rank: int = -1  # 3=experimental .. 0=automatic; -1=unstated


@dataclass
class Entry:
    accession: str
    name: str        # gene/entry name for display
    protein: str     # recommended protein name
    sequence: str
    features: list    # list[Feature]
    gene: str = ""   # primary gene symbol
    query: str = ""  # the input token that produced this entry (gene or accession)

    @property
    def length(self) -> int:
        return len(self.sequence)

    def of_kind(self, kind: str) -> list:
        return [f for f in self.features if f.kind == kind]


def _cache_path(accession: str) -> str:
    return os.path.join(CACHE_DIR, f"{accession.upper()}.json")


def _parse(raw: dict) -> Entry:
    accession = raw.get("primaryAccession", "")
    name = raw.get("uniProtkbId", accession)
    protein = ""
    try:
        protein = raw["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except (KeyError, TypeError):
        protein = name
    sequence = raw.get("sequence", {}).get("value", "")
    gene = ""
    try:
        gene = raw["genes"][0]["geneName"]["value"]
    except (KeyError, IndexError, TypeError):
        gene = name.split("_")[0]

    features: list[Feature] = []
    for f in raw.get("features", []):
        kind = FEATURE_TYPES.get(f.get("type"))
        if kind is None:
            continue
        loc = f.get("location", {})
        try:
            start = int(loc["start"]["value"])
            end = int(loc["end"]["value"])
        except (KeyError, TypeError, ValueError):
            continue  # skip features with unknown/fuzzy positions
        ev_label, ev_rank = _summarise_evidence(f.get("evidences", []))
        features.append(Feature(kind, start, end, f.get("description", ""),
                                ev_label, ev_rank))
    features.sort(key=lambda x: (x.start, x.end))
    return Entry(accession, name, protein, sequence, features, gene)


def _entry_from_dict(d: dict) -> Entry:
    # Tolerate cache files that predate the evidence/topology fields.
    feats = []
    for ff in d["features"]:
        feats.append(Feature(
            ff["kind"], ff["start"], ff["end"], ff.get("description", ""),
            ff.get("evidence", ""), ff.get("evidence_rank", -1)))
    # Older cache files predate the `gene` field; derive it from the entry name.
    gene = d.get("gene") or d.get("name", "").split("_")[0]
    return Entry(d["accession"], d["name"], d["protein"], d["sequence"], feats, gene)


def _cache_is_current(d: dict) -> bool:
    """A cache is current only if it was written with the present feature schema."""
    return d.get("schema_version") == SCHEMA_VERSION


def fetch(accession: str, use_cache: bool = True, timeout: int = 20) -> Entry:
    """Fetch a single UniProt entry. Raises ValueError on failure."""
    accession = accession.strip().upper()
    if not accession:
        raise ValueError("empty accession")

    os.makedirs(CACHE_DIR, exist_ok=True)
    path = _cache_path(accession)
    stale = None  # a cache in an older schema, kept as an offline fallback
    if use_cache and os.path.exists(path):
        with open(path, "r", encoding="utf-8") as fh:
            cached = json.load(fh)
        if _cache_is_current(cached):
            return _entry_from_dict(cached)
        stale = cached  # re-fetch to pick up evidence/topology; fall back if offline

    url = BASE_URL.format(acc=accession)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        if stale is not None:
            return _entry_from_dict(stale)
        raise ValueError(f"{accession}: HTTP {e.code} (check the accession)") from e
    except urllib.error.URLError as e:
        if stale is not None:
            return _entry_from_dict(stale)
        raise ValueError(f"{accession}: network error ({e.reason})") from e

    entry = _parse(raw)
    if not entry.sequence:
        raise ValueError(f"{accession}: no sequence returned")

    with open(path, "w", encoding="utf-8") as fh:
        json.dump(
            {
                "schema_version": SCHEMA_VERSION,
                "accession": entry.accession,
                "name": entry.name,
                "protein": entry.protein,
                "sequence": entry.sequence,
                "gene": entry.gene,
                "features": [asdict(f) for f in entry.features],
            },
            fh,
        )
    return entry


_GENE_CACHE = os.path.join(CACHE_DIR, "_gene_map.json")


def _load_gene_map() -> dict:
    if os.path.exists(_GENE_CACHE):
        try:
            with open(_GENE_CACHE, "r", encoding="utf-8") as fh:
                return json.load(fh)
        except (OSError, ValueError):
            return {}
    return {}


def _save_gene_map(m: dict) -> None:
    os.makedirs(CACHE_DIR, exist_ok=True)
    with open(_GENE_CACHE, "w", encoding="utf-8") as fh:
        json.dump(m, fh)


def resolve_gene(symbol: str, organism_id: int = 9606, reviewed: bool = True,
                 use_cache: bool = True, timeout: int = 20) -> str:
    """Resolve a gene symbol (e.g. 'CTSD') to a UniProt accession.

    Prefers reviewed (Swiss-Prot) human entries. Raises ValueError if none found.
    """
    symbol = symbol.strip().upper()
    key = f"{symbol}|{organism_id}|{int(reviewed)}"
    gmap = _load_gene_map()
    if use_cache and key in gmap:
        return gmap[key]

    query = f"gene_exact:{symbol} AND organism_id:{organism_id}"
    if reviewed:
        query += " AND reviewed:true"
    params = urllib.parse.urlencode({
        "query": query,
        "fields": "accession,gene_primary,protein_name,length",
        "format": "json",
        "size": "5",
    })
    url = f"{SEARCH_URL}?{params}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        raise ValueError(f"{symbol}: search HTTP {e.code}") from e
    except urllib.error.URLError as e:
        raise ValueError(f"{symbol}: network error ({e.reason})") from e

    results = data.get("results", [])
    if not results:
        raise ValueError(f"{symbol}: no reviewed {organism_id} entry found")

    # gene_exact also matches gene *synonyms*, so UniProt may rank a different
    # protein first (e.g. ARSA -> GET3). Prefer the entry whose PRIMARY gene name
    # equals the requested symbol; only fall back to the top hit if none match.
    def primary_gene(r):
        try:
            return r["genes"][0]["geneName"]["value"].upper()
        except (KeyError, IndexError, TypeError):
            return ""

    # Prefer an entry whose PRIMARY gene name equals the query (fixes synonym
    # mis-ranking like ARSA->GET3). If none match, the query is likely an old/alt
    # symbol (e.g. GBA -> GBA1); fall back to the top reviewed hit.
    chosen = next((r for r in results if primary_gene(r) == symbol), None)
    if chosen is None:
        chosen = results[0]
    acc = chosen.get("primaryAccession", "")
    if not acc:
        raise ValueError(f"{symbol}: no accession in search result")
    gmap[key] = acc
    _save_gene_map(gmap)
    return acc


def fetch_input(token: str, use_cache: bool = True) -> Entry:
    """Fetch by accession or gene symbol. The originating query is remembered on
    the returned Entry as `.query` for display."""
    token = token.strip()
    if looks_like_accession(token):
        entry = fetch(token, use_cache=use_cache)
        entry.query = token.upper()
        return entry
    acc = resolve_gene(token, use_cache=use_cache)
    entry = fetch(acc, use_cache=use_cache)
    entry.query = token.upper()
    return entry


def fetch_many(tokens, use_cache: bool = True):
    """Fetch several inputs (gene symbols and/or accessions). Returns
    (entries, errors) where errors is a list of (token, message) tuples."""
    entries, errors = [], []
    for tok in tokens:
        tok = tok.strip()
        if not tok:
            continue
        try:
            entries.append(fetch_input(tok, use_cache=use_cache))
        except ValueError as e:
            errors.append((tok, str(e)))
    return entries, errors


def parse_accessions(text: str) -> list[str]:
    """Split a free-text blob (commas, spaces, newlines) into tokens."""
    import re

    tokens = re.split(r"[\s,;]+", text.strip())
    seen, out = set(), []
    for t in tokens:
        t = t.strip().upper()
        if t and t not in seen:
            seen.add(t)
            out.append(t)
    return out
