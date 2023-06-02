from bio_structures import DNA_Codons, NUCLEOTIDE_BASE, RNA_Codons
from collections import Counter 
import random


class bio_seq:
    """DNA sequence class. Default Value: ATCG, DNA, No label"""
    
    def __init__(self, seq="ATCG", seq_type="DNA", label = 'No Label'):
        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be correct {self.seq_type}"
    
    def __validate(self):
        """check the sequence to make sure it is a valid DNA string"""
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)
    
    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type
    
    def get_seq_info(self):
        """Returns 4 strings. Full sequence information"""
        return f"[Label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}\n[Length]: {len(self.seq)}"
        
    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        """Generate a random DNA sequence, provided the length"""
        seq = ''.join([random.choice(NUCLEOTIDE_BASE[seq_type]) for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")
        
    def nucleotide_frequency(self):
        """Counts Nucleotides in a given sequence. Returns a dictionary"""
        return dict(Counter(self.seq))
    
    def transcription(self):
        """DNA -> RNA Transcription.Replacing Thymine with Uracil"""
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "Not a DNA Sequence"
    
    def reverse_compliment(self):
        """Swapping adenine with thymine and guanine with cytosine. Reversing newly generated string"""
        if self.seq_type == "DNA":
            mapping =str.maketrans('ATCG', 'TAGC')
        else:
            mapping = str.maketrans('AUCG', 'UAGC')
        return self.seq.translate(mapping)[::-1]
    
    def gc_content(self):
        """GC content in a DNA/RNA sequence"""
        return round((self.seq.count('C') + self.seq.count('G')) / len(self.seq) * 100)
    
    def gc_content_subsec(self, k=20):
        """GC content in a DNA/RNA sub-sequence length k. k = 20 by default"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(round((subseq.count('C') + subseq.count('G')) / len(subseq) * 100))
        return res
    
    def translate_seq(self, init_pos=0):
        """Translates a DNA sequence into an amino acid sequence"""
        if self.seq_type == "RNA":
            return [RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
        elif self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]


    def codon_usage(self, aminoacid):
        """Provides the frequency of each codon encoding a given amino acid in a DNA sequence"""
        tmp_list = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmp_list.append(self.seq[i:i + 3])
            
        if self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmp_list.append(self.seq[i:i + 3])
                    
        freq_dict = dict(Counter(tmp_list))
        total_wight = sum(freq_dict.values())
        for seq in freq_dict:
            freq_dict[seq] = round(freq_dict[seq] / total_wight, 2)
        return freq_dict
    
    def gen_reading_frames(self):
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_compliment(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames
    
    def proteins_from_rf(self, aa_seq):
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
    
    def all_proteins_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False):
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(self.seq[startReadPos: endReadPos, self.seq_type])
            rfs = tmp_seq.gen_reading_frames
        else:
            rfs = self.gen_reading_frames()
        
        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)
            
        if ordered:
            return sorted(res, key=len, reverse=True)
        return res