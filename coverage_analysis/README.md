# Protease coverage analysis

This directory is a from-scratch, scriptable replacement for the previous
workbook/charts. It answers the requested questions with seven Python-generated
figures and exports every plotted value to CSV.

Run from the repository root:

```powershell
python -m coverage_analysis.run_analysis
```

The first full run streams the required MaxQuant tables and four searched FASTA
files from the original study deposit (PXD024364/MSV000086944). The 19.9 GB ZIP
is **not** downloaded: the script reads only the compressed ZIP members it needs
and keeps reduced caches under `coverage_analysis/reference_cache/`. Later runs
use those caches. To produce only the theoretical figures without network access:

```powershell
python -m coverage_analysis.run_analysis --skip-nature
```

## Definitions

- Figures 1–3 use the bundled 2026-07-21 reviewed human UniProt FASTA.
- Detectable theoretical peptides are 7–52 amino acids with up to two missed
  cleavages, matching the common Spectronaut search window.
- Every displayed coverage value in figures 1–6 is residue-weighted sequence
  coverage (covered residues divided by all residues in the relevant proteome).
  CSV files also include the median per-protein coverage for comparison with the
  paper's reported medians.
- `Trypsin/Lys-C` means the union of two separately digested aliquots and is
  treated as one baseline group. All combinations are parallel digestions.
- Figure 2 adds one alternative enzyme to that baseline. Figure 3 adds two.
- Figures 4–6 use the original Nature study's MaxQuant `peptides.txt` and
  `proteinGroups.txt`, not PeptideAtlas. Peptides are attributed to an enzyme by
  the experiment PSM-count columns, then localized in the first majority protein
  sequence as described by the paper. Only its six enzymes are shown.
- Figure 7 is a 2 × 3 panel: the top three two-group designs on the first row and
  the top three three-group designs on the second. Each subplot reports the
  theoretical residue-weighted sequence coverage obtainable from peptides carrying at
  least one K+R, L, T, I, V, or H residue.

Pepsin, ProAlanase, elastase, and Arg-N are necessarily rule-based projections.
Pepsin in particular has broad, condition-dependent specificity; its graph value
uses the deterministic ExPASy pH 1.3 preference model and should not be read as a
high-confidence experimental prediction.

## Precursor / mature-state paired quantification

The preset-target analysis, including PGRN (`GRN`), is generated separately:

```powershell
python -m coverage_analysis.run_propeptide_analysis
```

It evaluates 34 conventional UniProt pro-peptide boundaries and 15 internal
PGRN/granulin-product boundaries. A junction counts as jointly quantifiable only
when the parallel aliquots provide both a human-proteome-unique peptide spanning
the intact precursor junction and a unique mature neo-terminal peptide. I/L are
collapsed for uniqueness. An analytical enzyme cut at the biological junction
does not make the ordinary terminal peptide mature-specific, although an intact
missed-cleavage spanning peptide can still report precursor.

Figures 8-10 show single-group outcomes and all two-/three-group parallel
designs. The combination figures split robust zero-missed-cleavage pairs from
pairs rescued by up to two missed cleavages. Figure 11 is the protein-by-enzyme
matrix, and figure 12 resolves every PGRN boundary. The detailed CSV gives the
selected peptide sequence, enzyme, coordinates, uniqueness, and missed-cleavage
count for both biological states.

## Nature-reference SILAC and pro-peptide outputs

The matching analyses using only peptides observed in PXD024364 are written to
the separate `coverage_analysis/nature_reference_output/` directory:

```powershell
python -m coverage_analysis.run_nature_extensions
```

The SILAC panel uses the top three observed two-group and three-group designs
from figures 5 and 6, then calculates residue coverage from observed peptides
carrying K+R, L, T, I, V, or H.

The pro-peptide figures match the exact theoretical candidate sequences to the
Nature `peptides.txt` table, retaining the experiment/enzyme assignment. The
study used a specific-cleavage MaxQuant search. It can therefore support
observed intact-junction-spanning peptides, but it did not search for the
semi-specific neo-terminal peptides required to identify the mature cleaved
state. A zero mature-state count in these outputs is a search-design limitation,
not evidence that biological processing did not occur. Testing that state
empirically requires a semi-specific/unspecific re-search or N-terminomics.

## Sources

- Sinitcyn et al., *Nature Biotechnology* 41, 1776–1786 (2023):
  <https://doi.org/10.1038/s41587-023-01714-x>
- Original data deposit: <https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD024364>
- Spectronaut 19 manual: <https://biognosys.com/content/uploads/2024/09/Spectronaut-19-manual-v4.pdf>
- ExPASy cleavage rules: <https://web.expasy.org/peptide_cutter/peptidecutter_enzymes.html>
- ProAlanase specificity: <https://www.promega.com/products/mass-spectrometry/proteases-and-surfactants/proalanase-mass-spec-grade/>
- Proteome Discoverer elastase rule: <https://docs.thermofisher.com/r/Proteome-Discoverer-3.1-User-Guide/en-US1324471691v1>
- Human progranulin (P28799): <https://www.uniprot.org/uniprotkb/P28799/entry>
