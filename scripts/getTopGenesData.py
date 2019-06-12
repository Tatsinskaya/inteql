# import needed libraries
import numpy as np
import gzip
import sys

# Read imput line
ind_file = sys.argv[1] 
rna_file = sys.argv[2]
out_file = sys.argv[3]


# For each gene test if it is in the top 5k genes
# If so, it is printed on a file
selected_ind = set(np.load(ind_file)[0].tolist())
c, p = -3, 0
out_file = open(out_file, 'w')
with gzip.open(rna_file) as infile:
    for line in infile:
        if c < -1:
            print(line.strip())
        elif c == -1:
            out_file.write(line)
        if c in selected_ind:
            out_file.write(line)
        c += 1
