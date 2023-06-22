complements = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N",
    "t": "t",
    "t": "a",
    "g": "c",
    "c": "g",
    "n": "n",
}

# Transversion, in molecular biology, refers to a point mutation in DNA
# in which a single (two ring) purine (A or G) is changed for a (one ring)
# pyrimidine (T or C), or vice versa.[1] A transversion can be spontaneous,
# or it can be caused by ionizing radiation or alkylating agents. It can
# only be reversed by a spontaneous reversion.
transversions = {
    "A": ["C", "T"],  ## A -> C or T
    "C": ["A", "G"],  ## C -> A or G
    "T": ["A", "G"],  ## T -> A or G
    "G": ["C", "T"],  ## G -> C or T
    "N": ["N", "N"],  ## nothing
}

# Transition, in genetics and molecular biology, refers to a point mutation
# that changes a purine nucleotide to another purine (A ↔ G), or a pyrimidine
# nucleotide to another pyrimidine (C ↔ T). Approximately two out of three
# single nucleotide polymorphisms (SNPs) are transitions.
transitions = {
    "A": "G",  ## A -> G
    "G": "A",  ## G -> A
    "C": "T",  ## C -> T
    "T": "C",  ## T -> C
    "N": "N",  ## none
}

atcgn = ["A", "T", "C", "G", "N"]

atcg = ["A", "T", "C", "G"]
