"""Curated gene-symbol presets of propeptide-bearing proteins.

Each list holds human gene symbols; the app resolves them to reviewed UniProt
accessions at fetch time. Membership was validated against UniProt molecule-
processing annotations (all listed entries carry a propeptide-activation junction,
N- or C-terminal). Edit freely.
"""

# Lysosomal / granule hydrolases synthesised as inactive proenzymes and activated
# by proteolytic propeptide removal (mostly cathepsins + a few glycosidases).
LYSOSOMAL_HYDROLASES = [
    "CTSB",   # Cathepsin B
    "CTSC",   # Dipeptidyl peptidase 1 (cathepsin C)
    "CTSD",   # Cathepsin D
    "CTSE",   # Cathepsin E
    "CTSF",   # Cathepsin F
    "CTSH",   # Cathepsin H
    "CTSK",   # Cathepsin K
    "CTSL",   # Cathepsin L
    "CTSO",   # Cathepsin O
    "CTSS",   # Cathepsin S
    "CTSV",   # Cathepsin V (CTSL2)
    "CTSW",   # Cathepsin W
    "CTSZ",   # Cathepsin Z/X
    "CTSG",   # Cathepsin G (granule serine protease)
    "LGMN",   # Legumain (C-terminal propeptide)
    "TPP1",   # Tripeptidyl-peptidase 1
    "GAA",    # Lysosomal acid alpha-glucosidase
    "GLB1",   # Beta-galactosidase
    "HEXA",   # Beta-hexosaminidase subunit alpha
    "HEXB",   # Beta-hexosaminidase subunit beta
    "GRN",    # Progranulin (PGRN; precursor of granulin peptides)
]

# 20S proteasome catalytic and structural beta subunits made with propeptides
# (autocatalytically/processively removed on assembly). Includes constitutive,
# immuno- (PSMB8/9/10) and thymo- (PSMB11) subunits.
PROTEASOME_SUBUNITS = [
    "PSMB1",   # beta6
    "PSMB4",   # beta7
    "PSMB5",   # beta5 (catalytic, chymotrypsin-like)
    "PSMB6",   # beta1 (catalytic, caspase-like)
    "PSMB7",   # beta2 (catalytic, trypsin-like)
    "PSMB8",   # beta5i / LMP7 (immunoproteasome)
    "PSMB9",   # beta1i / LMP2 (immunoproteasome)
    "PSMB10",  # beta2i / MECL-1 (immunoproteasome)
    "PSMB11",  # beta5t (thymoproteasome)
]

PRESETS = {
    "Lysosomal hydrolases": LYSOSOMAL_HYDROLASES,
    "Proteasome subunits": PROTEASOME_SUBUNITS,
}
