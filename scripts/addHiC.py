# import libraries needed
import sys
import pandas as pd

sys.path.append('../inqtl/')
from utils import *


# Read the input
in_folder = sys.argv[1]
pairs_file = sys.argv[2]
out_file = sys.argv[3]

pairs_df = pd.read_csv(pairs_file)
pairs_df = pairs_df[pairs_df['variant_id'].apply(lambda x: x.split('_')[0] != '22')]
pairs_df = pairs_df.sort_values(by=['variant_id'])



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
    # For each eQTL ets the chromosome, 
    # the position of the first element, x, 
    # and the position of the second element, y.
    
    chromosome = str(variantId2chrNum(pairs_df.loc[i,]['variant_id']))
    e = position2matrixBin(region2start(pairs_df.loc[i,]['enhancer_id']))
    p = position2matrixBin(region2start(pairs_df.loc[i,]['promoter_id']))
    x, y = sorted([p, e])
    if chromosome not in chromosomes:
        # Note that the data is sorted by chromosome, so, each contact matrix is computed once.
        # And used until next chromosome is found.
        chromosomes.add(chromosome)
        in_folder_chr = in_folder+'/chr'+chromosome+'/MAPQGE30'
        print(in_folder_chr)

        contactMatrix = ContactMatrix(in_folder_chr, chrm=chromosome, resl=5000)
        print(contactMatrix._matrixPath)


    # Check that the lower bin is on the first dimension
    # and the upper in the second dimension
    # If so, add the contact value 
    if x in contactMatrix.rawMatrix and y in contactMatrix.rawMatrix[x]:
        c+=1
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
pairs_df.to_csv(out_file, index=False)

# Just check it :)
print(len(normalized_kr), pairs_df.shape)