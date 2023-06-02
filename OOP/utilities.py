def colored(seq):
    bcolors = {
        'A': '\033[92m',
        'C': '\033[94m',
        'G': '\033[93m',
        'T': '\033[91m',
        'U': '\033[91m',
        'reset': '\033[0;0',
    }
    
    tmp_str = ""
    
    for nuc in seq:
        if nuc in bcolors:
            tmp_str += bcolors[nuc] + nuc
        else:
            tmp_str += bcolors['reset'] + nuc
            
        return tmp_str + '\033[0;0m'
    
def readTextFile(filePath):
    with open(filePath, 'r') as f:
        return "".join([l.strip() for l in f.readlines])
    
def writeTextFile(filePath, seq, mode='w'):
    with open(filePath, mode) as f:
        f.write(seq+ '\n')
        
def read_FASTA(filePath):
    with open(filePath, 'r') as f:
        FASTAFile = [l.strip() for l in f.readlines()]

    FASTADict = {}
    FASTALabel = ""

    for line in FASTAFile:
        if '>' in line:
            FASTALabel = line
            FASTADict[FASTALabel] = ""
        else:
            FASTADict[FASTALabel] += line

    return FASTADict
