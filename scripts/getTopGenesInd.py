import numpy as np
import gzip
import sys

in_file  = sys.argv[1] 
out_file = sys.argv[2]
l = int(sys.argv[3])


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

m = sorted(means, reverse=True)[l]
means = np.array(means)
selected_ind = np.where(means>m)

np.save(out_file, selected_ind)
