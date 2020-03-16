# Import needed libraries
import numpy as np
import gzip
import sys

# Read input line
ind_file = sys.argv[1]  # output01.npy
rna_file = sys.argv[2]  # GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
out_file = sys.argv[3]  # GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top.tab = output02.tab

# For each gene test if it is in the top genes defined in previous script
# If so, it is printed on a file
selected_ind = set(np.load(ind_file)[0].tolist())
c, p = -3, 0
out_file = open(out_file, 'w')
with gzip.open(rna_file) as infile:
    for line in infile:
        if c < -1:
            pass
        elif c == -1:
            out_file.write(line.decode('utf-8'))
        if c in selected_ind:
            out_file.write(line.decode('utf-8'))
        c += 1
