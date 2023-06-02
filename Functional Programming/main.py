from DNA_toolkit import *
from utilities import *
import random

rndDNAStr = ''.join([random.choice(nucleotides)
                    for nuc in range(50)])

DNAStr = validate_seq(rndDNAStr)


print(f'\nSequence: {DNAStr}\n')
print(f'[1] + Sequence Length: {len(DNAStr)}\n')
print(f'[2] + Nucleotide Frequencey: {count_nuc_freq(DNAStr)}\n')

print(f'[3] + DNA/RNA Transcription: {transcription(DNAStr)}\n')

print(f"[4] + DNA String + Complement:\n5' {DNAStr} 3' ")
print(f"   {''.join(['|' for c in range(len(DNAStr))])}")
print(f"3' {compliment(DNAStr)} 5'\n")
print(f"Reverse Compliment:\n5' {reverse_compliment(DNAStr)} 3'\n")
print(f'[5] + GC Content: {gc_content(DNAStr)}%\n')
print(f'[6] + GC Content in Subsection k=5: {gc_content_subsec(DNAStr, k=5)}\n')

print(f'[7] + Amino Acids Sequence from DNA: {translate_seq(DNAStr,0)}\n')

print(f'[8 + Codon frequency (L): {codon_usage(DNAStr, "L")}]\n')

print(f'[9] + Reading_frames:')
for frame in gen_reading_frames(DNAStr):
    print(frame)
    
print('\n[10] + All proteins in 6 open reading frames:')
for prot in all_proteins_from_orfs(NM_001185097, 0, 0, True):
    print(f'{prot}')


# test_rf_frame = [ 'L', 'M', 'T', 'A', 'L', 'V', 'V', 'L', 'L', 'R', 'R', 'G', 'S', '_', 'G', 'H',]

# print(proteins_from_rf(test_rf_frame))