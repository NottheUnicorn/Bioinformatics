import collections
from structures import *

#check the sequence to make sure it is a DNA string
def validate_seq(seq):
    tmpseq = seq.upper()
    for nuc in tmpseq:
        if nuc not in nucleotides:
            return False
    return tmpseq

def count_nuc_freq(seq):
    tmpFreqDict = {"A":0, "C":0, "G":0, "T":0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict
# return dict(collections.Counter(seq))

#DNA -> RNA Transcription
def transcription(seq):
    return seq.replace("T", "U")

def compliment(seq):
    return ''.join([dna_reverse_complement[nuc] for nuc in seq])

#Reverse Compliment Function
def reverse_compliment(seq):
    return ''.join([dna_reverse_complement[nuc] for nuc in seq])[::-1]

#mapping = str.maketrans('ATCG', 'TAGC')
#return seq.translate(mapping)[::-1]

def gc_content(seq):
    """GC Content in a DNA/RNA sequence"""
    return round((seq.count('C') + seq.count('G')) / len(seq) * 100)

def gc_content_subsec(seq, k=20):
    res= []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(gc_content(subseq))
    return res

def translate_seq(seq, init_pos=0):
    """Translates a DNA sequence into an amino acid"""
    return [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) -2,3)]

def codon_usage(seq, aminoacid):
    """Provides the frequency of each codon encoding a given amino acid in a DNA sequence"""
    tmp_list = []
    for i in range(0, len(seq) - 2, 3):
        if DNA_Codons[seq[i:i + 3]] == aminoacid:
            tmp_list.append(seq[i:i + 3])
            
        freq_dict = dict(collections.Counter(tmp_list))
        total_wight = sum(freq_dict.values())
        for seq in freq_dict:
            freq_dict[seq] = round(freq_dict[seq] / total_wight, 2)
        return freq_dict
    
def gen_reading_frames(seq):
    """Generate the six reading frames of a DNA sequence, including reverse compliment"""
    frames = []
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq, 1))
    frames.append(translate_seq(seq, 2))
    frames.append(translate_seq(reverse_compliment(seq), 0))
    frames.append(translate_seq(reverse_compliment(seq), 1))
    frames.append(translate_seq(reverse_compliment(seq), 2))
    return frames

def proteins_from_rf(aa_seq):
    """compute all possible proteins in an amino acid seq and return a list of possible proteins"""
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            #STOP accumalating amino acids if _* STOP codon was found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
            else:
                #START accumalating amino acids if M - START codon was found
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
    return proteins

# Generate all RF
# Extract all prots
# Return a list sorted/unsorted

def all_proteins_from_orfs(seq, startReadPos=0, endReadPos=0, ordered=False):
    """Compute all possible proteins for all open reading frames"""
    """Protein Search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2"""
    """API can be used to pull protein info"""
    if endReadPos > startReadPos:
        rfs = gen_reading_frames(seq[startRead: endRead])
    else:
        rfs = gen_reading_frames(seq)
        
    res = []
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)
            
    if ordered:
        return sorted(res, key=len, reverse=True)
    return res