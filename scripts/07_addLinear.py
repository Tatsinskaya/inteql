import sys
import pandas as pd
import pybedtools
from utils import *
import statistics
import gffpandas.gffpandas as gffpd

# Read input line
tpm_top_tab = sys.argv[1]   # '/nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output02.tab'
slope_top = sys.argv[2]     # '/nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output05.csv.gz'
bedfolder = sys.argv[3]     # '/nfs/research1/zerbino/jhidalgo/inteql/data/original-data/TargetFinder/'
outfile = sys.argv[4]       # '/nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output07.csv.gz'  GTEx_Analysis_v7_eQTL_EVB_linearData.csv.gz = output05.csv.gz
annotation_gff_file = sys.argv[5] #'/nfs/research1/zerbino/jhidalgo/inteql/data/original-data/homo_sapiens.GRCh37.GM12878.Regulatory_Build.regulatory_activity.20180925.gff'

genes_top = set()
header = True
with open(tpm_top_tab) as f:
    for line in f:
        if header:
            header = False
            continue
        genes_top.add(geneIdVersion2geneId(line.split()[0]))

data = pd.read_csv(slope_top, compression='gzip')
data['gene_id'] = data['gene_id'].apply(geneIdVersion2geneId)
# Filter by top genes
data = data[data['gene_id'].apply(lambda x: x in genes_top)]
data = data[pd.isna(data['Cells_EBV-transformed_lymphocytes']).apply(lambda x: not (x))]

# For the PCA and tsne I am going to use only the eQTLs presents everywhere dimData = data.dropna() dimData.shapedimData.to_csv('../data/GTEx_Analysis_v7_eQTL_allTissues_slope_top_complete.csv.gz', index=False, compression='gzip')
pairs_vg = tuple(zip(data['variant_id'], data['gene_id'].apply(geneIdVersion2geneId)))
# Variant id to enhancer table
enhancers_f = bedfolder+'enhancers.bed.gz'
# Variant bed table
variants = pybedtools.BedTool.from_dataframe(pd.DataFrame({
    'variant_chr': data['variant_id'].apply(variantId2chr),
    'variant_start': data['variant_id'].apply(variantId2pos),
    'variant_end': data['variant_id'].apply(variantId2end),
    'variant_id': data['variant_id']}) \
    .reindex(
    columns=['variant_chr', 'variant_start', 'variant_end', 'variant_id']))
# Enhancer bed table
enhancers = pybedtools.BedTool(enhancers_f)
# Intersection
variants_enhancers_names = ['variant_chr', 'variant_start', 'variant_end', 'variant_id',
                            'enhancer_chr', 'enhancer_start', 'enhancer_end', 'enhancer_id']
variants_enhancers = variants.intersect(enhancers, wa=True, wb=True, loj=True).to_dataframe(names=variants_enhancers_names)

# Save it as a dictionary
variantID__enhancersID = {}
for i in variants_enhancers.index:
    v = variants_enhancers.loc[i,]['variant_id']
    e = variants_enhancers.loc[i,]['enhancer_id']
    if v not in variantID__enhancersID:
        variantID__enhancersID[v] = [e]
    else:
        variantID__enhancersID[v].append(e)

# Read genes table
genes = pybedtools.BedTool(bedfolder+'genes.bed.gz')
# Read promoters table
promoters = pybedtools.BedTool(bedfolder+'promoters.bed.gz')
# Intersection
genes_promoter_names = ['gene_chr', 'gene_start', 'gene_end', 'gene_id', 'smth1', 'strand', 'annotation', 'type',
                        'smth2',
                        'promoter_chr', 'promoter_start', 'promoter_end', 'promoter_id']
genes_promoters = genes.intersect(promoters, wa=True, wb=True,loj=True).to_dataframe(names=genes_promoter_names)
# Save a dictionary for the future
geneID__promotersID = {}
for i in genes_promoters.index:
    g = genes_promoters.loc[i,]['gene_id']
    p = genes_promoters.loc[i,]['promoter_id']
    if g not in geneID__promotersID:
        geneID__promotersID[g] = [p]
    else:
        geneID__promotersID[g].append(p)

genes_promoters = genes_promoters[['gene_chr', 'gene_start', 'gene_end', 'gene_id', 'promoter_chr', 'promoter_start', 'promoter_end', 'promoter_id']]
genes_promoters = genes_promoters[genes_promoters['gene_id'].apply(lambda x: x in genes_top)]

# Get promoter to enhancer table
variantID = []
geneID = []
enhancerID = []
promoterID = []
# distance = []
var_prom_distance = []
var_enh_distance = []

c = 0
for v, g in pairs_vg:
    if v in variantID__enhancersID and g in geneID__promotersID:
        promoters_list = geneID__promotersID[g]
        enhancers_list = variantID__enhancersID[v]
        pairs_ep = tuple((enhancers_list[i], promoters_list[u]) for i in range(len(enhancers_list)) for u in range(len(promoters_list)))
        for e, p in pairs_ep:
            variantID.append(v)
            var_position = int(v.split("_", 2)[1])
            geneID.append(g)
            if p != '.':
                promoterID.append(p)
                prom_position = int(sum(list(int(x) for x in p.split(":", 1)[1].split("-", 1))) / 2)
                var_prom_distance.append(abs(var_position - prom_position))
            else:
                promoterID.append(0)
                var_prom_distance.append(0)
            if e != '.':
                enhancerID.append(e)
                enh_position = int(sum(list(int(x) for x in e.split(":",1)[1].split("-",1)))/2)
                var_enh_distance.append(abs(var_position - enh_position))
            else:
                enhancerID.append(0)
                var_enh_distance.append(0)
            c += 1

# Set a dataframe
pairs_df = pd.DataFrame({
    'variant_id': variantID,
    'gene_id': geneID,
    'enhancer_id': enhancerID,
    'promoter_id': promoterID,
    # 'distance': distance,
    'var_prom_distance': var_prom_distance,
    'var_enh_distance': var_enh_distance,})

pairs_df = pairs_df.drop_duplicates()

## Add linear data
peaks = pybedtools.BedTool(bedfolder+'peaks.bed.gz')
methylation = pybedtools.BedTool(bedfolder+'methylation.bed.gz')
cage = pybedtools.BedTool(bedfolder+'cage.bed.gz')
peaks = peaks.cat(*[methylation, cage], postmerge=False).sort()
enhancer_names = ['enhancer_chr', 'enhancer_start', 'enhancert_end', 'enhancer_id']
promoter_names = ['promoter_chr', 'promoter_start', 'promoter_end', 'promoter_id']
peak_names = ['peak_chr', 'peak_start', 'peak_end', 'peak_name', 'peak_value']
peaks_df = peaks.to_dataframe(names=peak_names)

enhancer_peaks = enhancers.intersect(peaks, loj=True, wa=True, wb=True) \
    .to_dataframe(names=enhancer_names + peak_names)[['enhancer_id', 'peak_name', 'peak_value']]
enhancer_peaks = enhancer_peaks.groupby(['enhancer_id', 'peak_name']).sum().reset_index()
enhancer_peaks = enhancer_peaks.pivot_table(index='enhancer_id', columns='peak_name', values='peak_value').reset_index().fillna(0).drop(columns='.')
enhancer_peaks.columns = ['enhancer_id'] + ['enhancer_' + i for i in enhancer_peaks.columns[1:]]
promoter_peaks = promoters.intersect(peaks, loj=True, wa=True, wb=True) \
    .to_dataframe(names=promoter_names + peak_names)[['promoter_id', 'peak_name', 'peak_value']]

promoter_peaks = promoter_peaks.groupby(['promoter_id', 'peak_name']).sum().reset_index()
promoter_peaks = promoter_peaks.pivot_table(index='promoter_id', columns='peak_name', values='peak_value') \
    .reset_index().fillna(0)
promoter_peaks.columns = ['promoter_id'] + ['promoter_' + i for i in promoter_peaks.columns[1:]]

####Add regulatory build data
annotation_gff = gffpd.read_gff3(annotation_gff_file).attributes_to_columns()[['seq_id','bound_start','bound_end','regulatory_feature_stable_id','type','activity']]
annotation_gff['seq_id'] = 'chr' + annotation_gff['seq_id'].astype(str)
annotation_gff['name'] = 'GM12878|'+annotation_gff['seq_id'].astype(str)+':'+annotation_gff['bound_start'].astype(str)+'-'+annotation_gff['bound_end'].astype(str)

bed_ann = pybedtools.BedTool.from_dataframe(annotation_gff)
variant_names = ['chrom','start','end','variant_id']
annotation_names = ['seq_id','bound_start','bound_end','regulatory_feature_stable_id','type','activity','annotation_id']
annotation_peaks = variants.intersect(bed_ann, wa=True, wb=True).to_dataframe(names=variant_names+annotation_names)[['variant_id','annotation_id','type','activity']].drop_duplicates()
ann_data = annotation_peaks.join(pd.get_dummies(annotation_peaks[['type','activity']]))
ann_data['type'] = 1
ann_data = ann_data.rename(columns={'type': 'RegElementInfo'}).drop(columns='activity')
pairs_df = pairs_df.merge(ann_data, on='variant_id', how='left').fillna(0)
pairs_df = pairs_df.merge(enhancer_peaks, on='enhancer_id', how='left')
pairs_df = pairs_df.merge(promoter_peaks, on='promoter_id', how='left')
pairs_df = pairs_df.drop_duplicates()
pairs_df.to_csv(outfile, index=False, compression='gzip')