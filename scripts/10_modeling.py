from __future__ import print_function
import pandas as pd
#import numpy as np
#from scipy.stats import linregress
#from sklearn.tree import export_graphviz
#from sklearn.ensemble import ExtraTreesRegressor
import sys
sys.path.append('../inteql/')
from modelFunctions import *

finemap_file = sys.argv[1]      #'/nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output09.csv.gz'
eqtl_file = sys.argv[2]         #'/nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output05.csv.gz'
chr = sys.argv[3]               # Chromosome todo unnecessary?
outfile = sys.argv[4]           # /nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output10_chr.csv.gz

#Finemap data
data_z = pd.read_csv(finemap_file)
eqtl = pd.read_csv(eqtl_file)
data_z = data_z.merge(eqtl, on=['variant_id', 'gene_id']).drop_duplicates()
data_z = data_z.fillna(0)
data_z['z'] = data_z['z'].astype(float)
data_z['Chromosome'] = data_z['variant_id'].str.split("_",1,expand=True)[0]

# Get list of features names for each subset
epigenomicFeatures = list(i for i in data_z.columns if i.startswith('enha') or i.startswith('prom'))[2:]
regbuild_typeFeatures = list(i for i in data_z.columns if i.startswith('type'))
regbuild_activityFeatures = list(i for i in data_z.columns if i.startswith('activity'))
regbuild_dataFeature = ['RegElementInfo']
hiCFeatures = list(i for i in data_z.columns if i.startswith('hi'))
eQTLFeatures = list(eqtl.columns[2:-1])
distanceFeature = ['var_prom_distance','var_enh_distance']
chrFeature = ['Chromosome']

combinations = {}
combinations['Epigenomics'] = [epigenomicFeatures]
#combinations['RegBuild Type'] = [regbuild_typeFeatures]
#combinations['RegBuild Act'] = [regbuild_activityFeatures]
#combinations['RegBuild Data'] = [regbuild_dataFeature]
#combinations['All Reg'] = [regbuild_dataFeature+regbuild_activityFeatures+regbuild_typeFeatures]
#combinations['All Reg + eQTL'] = [regbuild_dataFeature+regbuild_activityFeatures+regbuild_typeFeatures+eQTLFeatures]
#combinations['All Reg + Epi'] = [regbuild_dataFeature+regbuild_activityFeatures+regbuild_typeFeatures+epigenomicFeatures]
combinations['HiC'] = [hiCFeatures]
combinations['eQTL'] = [eQTLFeatures]
combinations['Distance'] = [distanceFeature]
combinations['Distance+Chromosome'] = [distanceFeature+chrFeature]
combinations['HiC + Epigenomic'] = [hiCFeatures+epigenomicFeatures]
combinations['Epi + eQTL'] = [epigenomicFeatures+eQTLFeatures]
combinations['eQTL + Distance'] = [eQTLFeatures+distanceFeature]
#combinations['Epi + eQTL + RegBuild Type'] = [epigenomicFeatures+distanceFeature+regbuild_typeFeatures]
#combinations['Epi + eQTL + RegBuild Act'] = [epigenomicFeatures+distanceFeature+regbuild_activityFeatures]
#combinations['Epi + eQTL + RegBuild Data'] = [epigenomicFeatures+distanceFeature+regbuild_dataFeature]
combinations['All'] = [epigenomicFeatures+hiCFeatures+eQTLFeatures+distanceFeature]

random_state = 42

with open(outfile,'w+') as f:
    print('Chromosome ' + chr, 'RMSE and R Value for Pre-Finemap and Finemap data', sep='\t', file=f)
    print('Data', 'RF RMSE Slope', 'RF RMSE Z', 'R-value Slope', 'R-Value Z-score', sep='\t', file=f)
    for i in combinations:
        X_label = [item for sublist in combinations[i] for item in sublist]
        random_forest = random_forest_regressor(data_z, X_label, 'Cells_EBV-transformed_lymphocytes', random_state=random_state)
        # decision_tree = decision_tree_regressor(data_z, X_label, 'Cells_EBV-transformed_lymphocytes', random_state=random_state)
        random_forest_z = random_forest_regressor(data_z, X_label, 'z', random_state=random_state)
        # decision_tree_z = decision_tree_regressor(data_z, X_label, 'z', random_state=random_state)
        print("{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}".format(i,random_forest['rmse'],random_forest_z['rmse'], random_forest['r_value'], random_forest_z['r_value']), file=f)
    dummy = dummy_regressor(data_z, X_label, 'z', random_state)
    print('Dummy: {:.4f}'.format(dummy['rmse']), file=f)
