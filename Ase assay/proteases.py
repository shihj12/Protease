"""In-silico protease digestion.

A small, transparent rules engine — no external dependencies. Each protease is
defined by which residues it cleaves next to and any neighbouring-residue
exceptions. Cleavage positions are computed once, then peptides (including
missed-cleavage variants) are generated with 1-based [start, end] coordinates
matching UniProt conventions.

Extend PROTEASES with new entries as needed; the digest logic reads the table.
"""

from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class Peptide:
    """A peptide with 1-based inclusive coordinates into the parent sequence."""

    start: int  # 1-based position of first residue
    end: int    # 1-based position of last residue (inclusive)
    seq: str
    missed: int  # number of missed cleavage sites contained

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def spans_bond(self, j: int) -> bool:
        """True if this peptide covers the peptide bond between residue j and j+1."""
        return self.start <= j and self.end >= j + 1


# --- Cleavage rule definitions -------------------------------------------------
#
# A rule decides, for the bond between seq[i] and seq[i+1] (0-based), whether the
# protease cleaves there. `residues` is the set the enzyme recognises; `terminus`
# says whether it cuts C-terminal ("C") or N-terminal ("N") to those residues;
# `blocked_by` lists neighbouring residues on the far side that suppress cleavage
# (e.g. trypsin does not cut before proline).


@dataclass(frozen=True)
class Protease:
    name: str
    residues: frozenset
    terminus: str            # "C" = cut after residue, "N" = cut before residue
    blocked_by: frozenset = frozenset()
    description: str = ""


PROTEASES: dict[str, Protease] = {
    "Trypsin": Protease(
        "Trypsin", frozenset("KR"), "C", frozenset("P"),
        "Cuts C-terminal to K/R, but not before P (classic)."
    ),
    "Trypsin/P": Protease(
        "Trypsin/P", frozenset("KR"), "C", frozenset(),
        "Cuts C-terminal to K/R including before P."
    ),
    "Lys-C": Protease(
        "Lys-C", frozenset("K"), "C", frozenset(),
        "Cuts C-terminal to K."
    ),
    "Arg-C": Protease(
        "Arg-C", frozenset("R"), "C", frozenset(),
        "Cuts C-terminal to R."
    ),
    "Glu-C (E)": Protease(
        "Glu-C (E)", frozenset("E"), "C", frozenset(),
        "V8 in ammonium bicarbonate: cuts C-terminal to E."
    ),
    "Glu-C (D/E)": Protease(
        "Glu-C (D/E)", frozenset("DE"), "C", frozenset(),
        "V8 in phosphate buffer: cuts C-terminal to D and E."
    ),
    "Chymotrypsin": Protease(
        "Chymotrypsin", frozenset("FYW"), "C", frozenset("P"),
        "Cuts C-terminal to F/Y/W, not before P."
    ),
    "Asp-N": Protease(
        "Asp-N", frozenset("D"), "N", frozenset(),
        "Cuts N-terminal to D."
    ),
    "Lys-N": Protease(
        "Lys-N", frozenset("K"), "N", frozenset(),
        "Cuts N-terminal to K."
    ),
}


def cleavage_sites(seq: str, protease: Protease) -> list[int]:
    """Return sorted 1-based positions J such that the enzyme cuts the bond J|J+1.

    J is the residue number after which cleavage occurs (so the peptide bond
    broken is between residue J and residue J+1).
    """
    seq = seq.upper()
    sites: list[int] = []
    n = len(seq)
    for i in range(n - 1):  # bond between seq[i] and seq[i+1] -> J = i+1
        left, right = seq[i], seq[i + 1]
        if protease.terminus == "C":
            recognised = left in protease.residues
            far = right  # residue on the far (C) side of the bond
        else:  # "N" — cut before the recognised residue
            recognised = right in protease.residues
            far = left
        if recognised and far not in protease.blocked_by:
            sites.append(i + 1)  # 1-based J
    return sites


def is_cleavage_site(seq: str, protease: Protease, j: int) -> bool:
    """True if the enzyme cuts the bond between residue j and j+1 (1-based j)."""
    return j in set(cleavage_sites(seq, protease))


def cleavage_sites_multi(seq: str, proteases) -> list[int]:
    """Union of cut sites for several proteases used together (a co-digestion).

    Simultaneous/sequential digestion with two enzymes cuts wherever *either*
    would — i.e. the union of their sites. `proteases` is an iterable of Protease
    objects or names.
    """
    s: set[int] = set()
    for p in proteases:
        if isinstance(p, str):
            p = PROTEASES[p]
        s.update(cleavage_sites(seq, p))
    return sorted(s)


def _peptides_from_sites(seq: str, sites: list[int], missed_cleavages: int,
                         min_len: int, max_len: int) -> list[Peptide]:
    """Generate peptides (with missed cleavages) from a sorted list of cut sites.

    Shared by single- and multi-protease digestion. `missed_cleavages` counts
    internal sites left uncut. Coordinates are 1-based inclusive.
    """
    n = len(seq)
    if n == 0:
        return []
    # Boundaries between consecutive peptides, as 1-based end positions.
    # Fragment k spans (bounds[k]+1 .. bounds[k+1]) in 1-based coords.
    bounds = [0] + list(sites) + [n]

    peptides: list[Peptide] = []
    n_frag = len(bounds) - 1
    for i in range(n_frag):
        for mc in range(missed_cleavages + 1):
            k = i + mc
            if k >= n_frag:
                break
            start = bounds[i] + 1
            end = bounds[k + 1]
            length = end - start + 1
            if length < min_len or length > max_len:
                continue
            peptides.append(Peptide(start, end, seq[start - 1:end], mc))
    return peptides


def digest(
    seq: str,
    protease: Protease | str,
    missed_cleavages: int = 2,
    min_len: int = 7,
    max_len: int = 40,
) -> list[Peptide]:
    """Digest `seq`, returning peptides with up to `missed_cleavages` missed sites.

    Peptides are filtered to [min_len, max_len]. Coordinates are 1-based inclusive.
    """
    if isinstance(protease, str):
        protease = PROTEASES[protease]
    seq = seq.upper()
    return _peptides_from_sites(seq, cleavage_sites(seq, protease),
                                missed_cleavages, min_len, max_len)


def digest_multi(
    seq: str,
    proteases,
    missed_cleavages: int = 2,
    min_len: int = 7,
    max_len: int = 40,
) -> list[Peptide]:
    """Co-digestion with several enzymes at once: peptides bounded by the union of
    all enzymes' cut sites. Yields fragments neither enzyme produces alone."""
    seq = seq.upper()
    return _peptides_from_sites(seq, cleavage_sites_multi(seq, proteases),
                                missed_cleavages, min_len, max_len)


def protease_names() -> list[str]:
    return list(PROTEASES.keys())
