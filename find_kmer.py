import sys
from Bio import SeqIO
if len(sys.argv) != 3:
    print("usage: python find_kmer.py sequences.fasta kmer")

    sys.exit()

file_name = sys.argv[1]
kmer = sys.argv[2]
with open(file_name, "r") as f:
    i = 0
    for record in SeqIO.parse(f, 'fasta'):
        line = record.seq
        start_index = 0
        keepGoing = True
        while keepGoing:
            loc = line.find(kmer, start_index)
            if loc != -1:
                print('line: {0}, index: {1}'.format(i, loc))
                start_index = loc + 1
            else:
                keepGoing = False
        i += 1
