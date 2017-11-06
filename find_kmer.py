import sys
if len(sys.argv) != 3:
    print("usage: python find_kmer.py sequences.txt kmer")
    print("Note that sequences.txt should have a sequence on each line; not FASTA")
    sys.exit()

file_name = sys.argv[1]
kmer = sys.argv[2]
with open(file_name, "r") as f:
    i = 0
    for line in f:
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
