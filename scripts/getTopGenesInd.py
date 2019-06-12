# Import needed libraries
import numpy as np
import gzip
import sys

# Read imput line
in_file  = sys.argv[1] 
out_file = sys.argv[2]
l = int(sys.argv[3])


# Compute the mean expression of each gene and save it on an array
c = 0
means = []
with gzip.open(in_file) as infile:
    for line in infile:
        c += 1
        if c <= 3:  continue
        line = line.strip().split()
        if len(line) < 40: continue
        line = tuple(float(i) for i in line[2:])
        means.append(np.mean(line))

# Select the id of the tissues that are in the top 5k
m = sorted(means, reverse=True)[l]
means = np.array(means)
selected_ind = np.where(means>m)

# Save the results
np.save(out_file, selected_ind)