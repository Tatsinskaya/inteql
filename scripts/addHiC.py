import pandas as pd
def region2start(x): return x.split(':')[1].split('-')[0]
import methodName

import sys
in_folder = sys.argv[1]
pairs_file = sys.argv[2]
out_file = sys.argv[3]



raw = []
normalized_kr = []
normalized_vc = []
normalized_sq = []
normalizedExpected_kr = []
normalizedExpected_vc = []
normalizedExpected_sq = []

pairs_df = pd.read_csv(pairs_file)
pairs_df = pairs_df[pairs_df['variant_id'].apply(lambda x: x.split('_')[0] != '22')]
pairs_df = pairs_df.sort_values(by=['variant_id'])


j, c = 0, 0
hi_c = [0] * pairs_df.shape[0]
chromosomes = set()
for i in pairs_df.index:
    chromosome = str(methodName.variantId2chrNum(pairs_df.loc[i,]['variant_id']))
    e = methodName.position2matrixBin(region2start(pairs_df.loc[i,]['enhancer_id']))
    p = methodName.position2matrixBin(region2start(pairs_df.loc[i,]['promoter_id']))
    x, y = sorted([p, e])
    if chromosome not in chromosomes:
        chromosomes.add(chromosome)
        in_folder_chr = in_folder+'/chr'+chromosome+'/MAPQGE30'
        print(in_folder_chr)

        contactMatrix = methodName.ContactMatrix(in_folder_chr, chrm=chromosome, resl=5000)
        print(contactMatrix._matrixPath)


    # check that the lower bin is on the first dimension
    # and the upper in the second dimension
    if x in contactMatrix.rawMatrix and y in contactMatrix.rawMatrix[x]:
        c+=1
        raw.append(contactMatrix.rawMatrix[x][y]) #KRnormMatrix
        normalizedExpected_kr.append(contactMatrix.KRnormExpectedMatrix[x][y])
        normalizedExpected_vc.append(contactMatrix.VCnormExpectedMatrix[x][y])
        normalizedExpected_sq.append(contactMatrix.SQRTVCnormExpectedMatrix[x][y])
        normalized_kr.append(contactMatrix.KRnormMatrix[x][y])
        normalized_vc.append(contactMatrix.VCnormMatrix[x][y])
        normalized_sq.append(contactMatrix.SQRTVCnormMatrix[x][y])
    else:
        raw.append(0)
        normalizedExpected_kr.append(0)
        normalizedExpected_vc.append(0)
        normalizedExpected_sq.append(0)
        normalized_kr.append(0)
        normalized_vc.append(0)
        normalized_sq.append(0)


pairs_df['hiC_raw'] = raw
pairs_df['hiC_kr'] = normalized_kr
pairs_df['hiC_vc'] = normalized_vc
pairs_df['hiC_sq'] = normalized_sq
pairs_df['hiC_krExpected'] = normalizedExpected_kr
pairs_df['hiC_vcExpected'] = normalizedExpected_vc
pairs_df['hiC_sqExpected'] = normalizedExpected_sq

pairs_df.to_csv(out_file, index=False)

print(len(normalized_kr), pairs_df.shape)
