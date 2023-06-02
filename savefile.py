from bio_seq import bio_seq
from utilities import read_FASTA, readTextFile, writeTextFile


writeTextFile("test.txt", test_dna.seq)
for rf in test_dna.gen_reading_frames():
    writeTextFile("test.txt", test_dna.seq, str(rf), 'a')
    
fasta = read_FASTA("fasta_samples.txt")
print(fasta)
