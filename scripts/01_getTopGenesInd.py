# Import needed libraries
import numpy as np
import gzip
import sys
import pyensembl

# import requests, json

# def fetch_endpoint(server, request, content_type):
#     r = requests.get(server + request, headers={"Content-Type": content_type})
#     if not r.ok:
#         r.raise_for_status()
#         sys.exit()
#     return r.json()

# Read imput line
in_file = sys.argv[1]
out_file = sys.argv[2]
length = int(sys.argv[3])

# server = "http://grch37.rest.ensembl.org/"
# con = "application/json"
ensembl = pyensembl.EnsemblRelease(release=75)
excluding = [["6", 28477797, 33448354], ["17", 44165260, 44784489],
             ["22", 42371706, 43416680]]  # "CHROMOSOME", START POSITION, END POSITION

# Compute the mean expression of each gene and save it on an array
c = 0
means = []
with gzip.open(in_file) as infile:
    for line in infile:
        c += 1
        if c <= 3:  continue
        line = line.strip().split()
        if len(line) < 40: continue
        geneid = line[0].decode('UTF-8').split(".")[0]
        try:
            gene = ensembl.gene_by_id(geneid)
        except:
            print('Gene ', geneid, 'could not be found.')
            continue
        break_loop = False
        if any(gene.contig in sublist for sublist in excluding):
            for i in excluding:
                if gene.contig == i[0]:
                    if i[1] <= ((gene.start + gene.end) / 2) <= i[2]:
                        break_loop = True
                        continue
        if break_loop:
            print('Skipped ', geneid, gene.contig, gene.start)
            means.append(0)
        else:
            line = tuple(float(i) for i in line[2:])
            means.append(np.mean(line))

# Select the id of the tissues that are in the top 5k
m = sorted(means, reverse=True)[length]
means = np.array(means)
selected_ind = np.where(means > m)

# Save the results
np.save(out_file, selected_ind)
