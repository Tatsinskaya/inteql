# import libraries needed
import sys
import pandas as pd
import gffpandas.gffpandas as gffpd
import pybedtools
sys.path.append('../inteql/')
from utils import *
from contactMatrix import *

# Read the input
in_folder = sys.argv[1]  # data/original-data/GM12878_combined/5kb_resolution_intrachromosomal
pairs_file = sys.argv[2]  # output07.csv
out_file = sys.argv[3]  # GTEx_Analysis_v7_eQTL_EVB_linearData_hiC.csvz = output08.csv
bedfolder = sys.argv[4]  # '/nfs/research1/zerbino/jhidalgo/inteql/data/original-data/TargetFinder/'
regfolder = sys.argv[5] # '/nfs/research1/zerbino/jhidalgo/inteql/data/original-data/RegBuild/'

pairs_df = pd.read_csv(pairs_file, compression='gzip')
# pairs_df = pairs_df[pairs_df['variant_id'].apply(lambda x: x.split('_')[0] != '22')]
pairs_df = pairs_df.sort_values(by=['variant_id'])
genes = pybedtools.BedTool(bedfolder+'genes.bed.gz').to_dataframe()
enhancers = gffpd.read_gff3(regfolder+'enhancer.gff').attributes_to_columns()[['seq_id','bound_start','bound_end','regulatory_feature_stable_id','type','activity']]
promoters = gffpd.read_gff3(regfolder+'promoter.gff').attributes_to_columns()[['seq_id','bound_start','bound_end','regulatory_feature_stable_id','type','activity']]

# Empty columns
raw = []
normalized_kr = []
normalized_vc = []
normalized_sq = []
normalizedExpected_kr = []
normalizedExpected_vc = []
normalizedExpected_sq = []

j, c = 0, 0
hi_c = [0] * pairs_df.shape[0]
chromosomes = set()
for i in pairs_df.index:
    # For each eQTL sets the chromosome,
    # the position of the first element, x, 
    # and the position of the second element, y.
    # if pairs_df.loc[i, ]['enhancer_id'] == '0' or pairs_df.loc[i,]['promoter_id'] == '0':
    if pairs_df.loc[i, ]['enhancer_id'] == '0' or pairs_df.loc[i,]['gene_id'] == '0':
        raw.append(0)
        normalizedExpected_kr.append(0)
        normalizedExpected_vc.append(0)
        normalizedExpected_sq.append(0)
        normalized_kr.append(0)
        normalized_vc.append(0)
        normalized_sq.append(0)
        continue
    if pairs_df.loc[i,]['promoter_id'] == '0':
        gene = genes[genes['name'] == pairs_df.loc[i,]['gene_id']]
        if gene['strand'].item() == '-':
            p = position2matrixBin(int(gene['end'].item()))
        if gene['strand'].item() == '+':
            p = position2matrixBin(int(gene['start'].item()))
    else:
        promoter = promoters[promoters['regulatory_feature_stable_id'] == pairs_df.loc[i, ]['promoter_id']]
        p_start=int(promoter['bound_start'].item())
        p_end=int(promoter['bound_end'].item())
        p = position2matrixBin(int(p_start+(p_end-p_start)/2)) #todo function?
    enhancer = enhancers[enhancers['regulatory_feature_stable_id'] == pairs_df.loc[i, ]['enhancer_id']]
    chromosome = str(enhancer['seq_id'].item())
    e_start=int(enhancer['bound_start'].item())
    e_end=int(enhancer['bound_end'].item())
    e = position2matrixBin(int(e_start+(e_end-e_start)/2))
    x, y = sorted([p, e])
    if chromosome not in chromosomes:
        # Note that the data is sorted by chromosome, so, each contact matrix is computed once.
        # And used until next chromosome is found.
        chromosomes.add(chromosome)
        in_folder_chr = in_folder + '/chr' + chromosome + '/MAPQGE30'
        # print(in_folder_chr)
        contactMatrix = ContactMatrix(in_folder_chr, chrm=chromosome, resl=5000)  # todo adjustable resolution (also on position2matrixBin)
        # print(contactMatrix._matrixPath)
    # Check that the lower bin is on the first dimension
    # and the upper in the second dimension
    # If so, add the contact value 
    if x in contactMatrix.rawMatrix and y in contactMatrix.rawMatrix[x]:
        c += 1
        raw.append(contactMatrix.rawMatrix[x][y])
        normalizedExpected_kr.append(contactMatrix.KRnormExpectedMatrix[x][y])
        normalizedExpected_vc.append(contactMatrix.VCnormExpectedMatrix[x][y])
        normalizedExpected_sq.append(contactMatrix.SQRTVCnormExpectedMatrix[x][y])
        normalized_kr.append(contactMatrix.KRnormMatrix[x][y])
        normalized_vc.append(contactMatrix.VCnormMatrix[x][y])
        normalized_sq.append(contactMatrix.SQRTVCnormMatrix[x][y])
    else:
        # Other wise the contact is 0.
        raw.append(0)
        normalizedExpected_kr.append(0)
        normalizedExpected_vc.append(0)
        normalizedExpected_sq.append(0)
        normalized_kr.append(0)
        normalized_vc.append(0)
        normalized_sq.append(0)

# Add the new columns to the data set
pairs_df['hiC_raw'] = raw
pairs_df['hiC_kr'] = normalized_kr
pairs_df['hiC_vc'] = normalized_vc
pairs_df['hiC_sq'] = normalized_sq
pairs_df['hiC_krExpected'] = normalizedExpected_kr
pairs_df['hiC_vcExpected'] = normalizedExpected_vc
pairs_df['hiC_sqExpected'] = normalizedExpected_sq

# Save data.
pairs_df.to_csv(out_file, index=False, compression='gzip')

# Just check it :)
# print(len(normalized_kr), pairs_df.shape)
