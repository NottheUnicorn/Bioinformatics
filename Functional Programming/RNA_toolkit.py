nucleotides = ["A", "U", "G", "T"]

#check the sequence to make sure it is a DNA string
def validate_seq(rna_seq):
    tmpseq = rna_seq.upper()
    for nuc in tmpseq:
        if nuc not in nucleotides:
            return False
    return tmpseq