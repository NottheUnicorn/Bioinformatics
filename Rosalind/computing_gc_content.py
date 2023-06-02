# read data from the (FATSA formatted file)
def read_file(filepath):
    with open(filepath, 'r') as f:
        return [l.strip() for l in f.readlines()]
    
def gc_content(seq):
    return ((seq.count('C') + seq.count('G')) / len(seq) * 100)
    
# Calculate4 GC content from FATSA formatted text file:

# Read data from the file(FATSA formatted file)
# storing file in a list
FASTA_file = read_file(r"C:\Users\Mozil\Desktop\New folder\Rebel_Bioinformatics\Rosalind\rosalind_gc (2).txt")
# Dictionary for labels and Data
FASTA_dict = {}
# String holding for current label 
FASTA_label = ""

#Converting FATSA/List file data into a dictionary
for line in FASTA_file:
    if '>' in line:
        FASTA_label = line
        FASTA_dict[FASTA_label] = ""
    else:
        FASTA_dict[FASTA_label] += line
        
print(FASTA_dict)

# Clean/prepare the data (format with gc_content function)
# format the data (Store the data in a convient way)
# Run needed operation of Data (pass the data to our gc_content function)

RESULT_dict = {key: gc_content(value) for (key,value) in FASTA_dict.items()}

print(RESULT_dict)


# Collect results

max_gc_key = max(RESULT_dict, key=RESULT_dict.get)

print(f'{max_gc_key}\n{RESULT_dict[max_gc_key]}')

