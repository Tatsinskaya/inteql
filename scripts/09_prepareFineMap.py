import sys
import numpy as np
import pandas as pd
import scipy
from pathlib import Path
import subprocess
import time

sys.path.append("/homes/jhidalgo/lib/postgap/lib/")
sys.path.append("/Users/jhidalgo/Desktop/inteql/exec/Dependencies/postgap-master/lib/")

import postgap
import postgap.Globals
import postgap.Integration
from postgap.Utils import *
import requests
import os

#### Argument passing ####
postgap.Globals.DATABASES_DIR = sys.argv[1]   # /nfs/research1/zerbino/jhidalgo/databases
pairs_f = sys.argv[2]     # data/inter_data/data/inter_data/output08.csv.gz
chromosome = sys.argv[3]     # 4
dbsnp_dir = sys.argv[4]       # '/Users/jhidalgo/ebi-cli/nfs/jhidalgo/inteql/data/original-data/dbSNP/chr-header'
output_file = sys.argv[5]       # GTEx_Analysis_v7_eQTL_EVB_linearData_hiC_z = output09.csv.gz
output_dir = sys.argv[6]+'chr_'+chromosome+'/'
genepairs = sys.argv[7]         #'/nfs/research1/zerbino/jhidalgo/inteql/data/original-data/GTEx_Analysis_v7_eQTL/Cells_EBV-transformed_lymphocytes.v7.signif_variant_gene_pairs.txt.gz'

dbsnp_file = dbsnp_dir + 'GTEx_2016_dbSNPs_table_chr_' + chromosome + '.txt'
server = "http://rest.ensembl.org"
postgap.Globals.SPECIES = 'human'
rs_name = 'rs_id_dbSNP147_GRCh37p13'
data = pd.read_csv(pairs_f, compression='gzip')
data = data[data['variant_id'].apply(lambda x: x.split('_')[0] == chromosome)]

rs = pd.read_csv(dbsnp_file, sep='\t')
rs.columns = ['chr', 'variant_pos', 'variant_id', 'ref', 'alt', 'num_alt_per_site', 'rs_id_dbSNP147_GRCh37p13']
data = data.merge(rs, on=['variant_id'], how='inner')
rsIDs = data[rs_name].unique()
# extract SNP objects
snps = postgap.Ensembl_lookup.get_snp_locations(rsIDs)
# get LD matrix
rsIDs, matrix = postgap.LD.get_pairwise_ld(snps, 'EUR')
len(rsIDs), len(set(rsIDs))
data = data.loc[data[rs_name].isin(rsIDs),]
M = {}
for i, a in enumerate(rsIDs):
    if a not in M: M[a] = {}
    for j, b in enumerate(rsIDs):
        M[a][b] = matrix[i, j]

data['rs_name'] = data[rs_name] + '_' + data['promoter_id'] + data['enhancer_id']
data = data.drop_duplicates(subset='rs_name')
c = 0
n = 0
id_beta = {}
decoded = []
for i in data.index:
    gene_stbl_id = data.loc[i,]['gene_id']
    variant_name = data.loc[i,][rs_name]
    ext = "/eqtl/id/homo_sapiens/" + gene_stbl_id + "?statistic=beta;variant_name=" + variant_name + ";tissue=Cells_EBV-transformed_lymphocytes"
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    while r.status_code == 429:
        time.sleep(0.2)
        #print("Slept 0.2 seconds to avoid Rate limit, "+r.headers['X-RateLimit-Remaining']+" remaining. X-RateLimit-Reset: "+r.headers['X-RateLimit-Reset'])
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    if not r.ok:
        try:
            r.raise_for_status()
        except:
            n += 1
    try:
        decoded = r.json()
    except:
        n += 1
    if "error" not in repr(decoded) and decoded != []:
        if (gene_stbl_id, variant_name) not in id_beta:
            id_beta[(gene_stbl_id, variant_name)] = decoded[0]['value']
            # print(decoded)
        else:
            if decoded[0]['value'] != id_beta[(gene_stbl_id, variant_name)]:
                pass
                # print('ERROR!', decoded[0]['value']," != ",id_beta[(gene_stbl_id, variant_name)])
        c += 1
    # else:
        # print(decoded)


ib = {}
keys = tuple(id_beta.keys())
ib['gene_id'], ib[rs_name] = zip(*keys)
ib['beta'] = tuple(id_beta[k] for k in keys)
data = data.merge(pd.DataFrame(ib), right_on=['gene_id', rs_name], left_on=['gene_id', rs_name])
E = []
for a in data[rs_name]:
    e = []
    for b in data[rs_name]:
        e.append(M[a][b])
    E.append(e)

E = np.array(E)
evb = pd.read_csv(genepairs,sep='\t', compression='gzip')
evb['gene_id'] = evb.gene_id.apply(lambda x: x.split('.')[0])
data = data.merge(evb, on=['variant_id', 'gene_id'], how='inner')
data_z = data[['rs_name', 'chr', 'variant_pos', 'ref', 'alt', 'maf', 'beta', 'slope_se']]
data_z.columns = ['rsid', 'chromosome', 'position', 'allele1', 'allele2', 'maf', 'beta', 'se']
#data_z.shape, E.shape
data_z.to_csv(output_dir + 'output09-chromosome' + chromosome + '_to_run_finemap.z', index=False, sep=' ')
np.savetxt(output_dir + 'output09-chromosome' + chromosome + '_to_run_finemap.ld', E, delimiter=' ')
# write the master:
headers = 'z;ld;snp;config;cred;log;n_samples\n'
names = output_dir + 'output09-chromosome'+chromosome+'_to_run_finemap.z;' + output_dir + 'output09-chromosome'+chromosome+'_to_run_finemap.ld;' + output_dir + 'output09-chromosome'+chromosome+'_to_run_finemap.snp;' + output_dir + 'output09-chromosome'+chromosome+'_to_run_finemap.config;' + output_dir + 'output09-chromosome'+chromosome+'_to_run_finemap.cred;' + output_dir + 'output09-chromosome'+chromosome+'_to_run_finemap.log;'+str(data_z.shape[0])
master_file = output_dir + 'output09-master_' + chromosome
with open(master_file, 'w') as f:
    f.write(headers + names)
#print(data_z.shape[0])

#FINEMAP#
finemap_call = 'finemap --sss --in-files ' + master_file + ' --dataset 1 --log'
print('Master: '+master_file)
print('Finemap: '+finemap_call)
subprocess.call(finemap_call, shell=True)
chr_file = output_dir + 'output09-chromosome' + chromosome + '_to_run_finemap.snp'
max_snps = 4
while os.stat(chr_file).st_size == 0:
    print('ERROR: '+chr_file+' is empty, reducing maximum number of allowed causal SNPs to: '+str(max_snps))
    finemap_call = 'finemap --sss --in-files ' + master_file + ' --dataset 1 --log --n-causal-snps '+str(max_snps)
    print('Finemap: ' + finemap_call)
    subprocess.call(finemap_call, shell=True)
    max_snps=max_snps-1
    if max_snps == 2:
        print('ERROR: '+chr_file+' could not be produced, omitting chromosome.')
        sys.exit()


snp = pd.read_csv(chr_file, sep=' ', index_col=0)
data = data.merge(snp, right_on='rsid', left_on='rs_name')
data.to_csv(output_file, index=False, compression='gzip')


