from __future__ import print_function
import pandas as pd
import sys
from sklearn.tree import export_graphviz
#import numpy as np
#from scipy.stats import linregress
#from sklearn.ensemble import ExtraTreesRegressor
sys.path.append('../inteql/')
from modelFunctions import *

finemap_file = sys.argv[1]      #'/nfs/research1/zerbino/jhidalgo/inteql/data/output/output09.csv.gz'
eqtl_file = sys.argv[2]         #'/nfs/research1/zerbino/jhidalgo/inteql/data/output/output05.csv.gz'
outfile = sys.argv[3]           # /nfs/research1/zerbino/jhidalgo/inteql/data/output/output10_chr.csv.gz
outfolder = sys.argv[4]

#Finemap data
data_z = pd.read_csv(finemap_file)
eqtl = pd.read_csv(eqtl_file)
data_z = data_z.merge(eqtl, on=['variant_id', 'gene_id']).drop_duplicates()
data_z = data_z.fillna(0)
data_z['z'] = data_z['z'].astype(float)
data_z['Chromosome'] = data_z['variant_id'].str.split("_",1,expand=True)[0]
data_z['position'] = data_z['variant_id'].str.split('_',expand=True)[1]

# Get list of features names for each subset
epigenomicFeatures = list(i for i in data_z.columns if i.startswith('enha') or i.startswith('prom'))[2:]
regbuild_typeFeatures = list(i for i in data_z.columns if i.startswith('type'))
regbuild_activityFeatures = list(i for i in data_z.columns if i.startswith('activity'))
regbuild_dataFeature = ['RegElementInfo']
hiCFeatures = list(i for i in data_z.columns if i.startswith('hi'))
eQTLFeatures = list(eqtl.columns[2:-1])
distanceFeature = ['var_prom_distance','var_enh_distance']
chrFeature = ['Chromosome']
posfeature = ['position']

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
# combinations['Distance'] = [distanceFeature]
# combinations['Distance+Chromosome'] = [distanceFeature+chrFeature]
# combinations['HiC + Epigenomic'] = [hiCFeatures+epigenomicFeatures]
combinations['Epi + eQTL'] = [epigenomicFeatures+eQTLFeatures]
combinations['eQTL + Distance'] = [eQTLFeatures+distanceFeature]
#combinations['Epi + eQTL + RegBuild Type'] = [epigenomicFeatures+distanceFeature+regbuild_typeFeatures]
#combinations['Epi + eQTL + RegBuild Act'] = [epigenomicFeatures+distanceFeature+regbuild_activityFeatures]
#combinations['Epi + eQTL + RegBuild Data'] = [epigenomicFeatures+distanceFeature+regbuild_dataFeature]
combinations['All'] = [epigenomicFeatures+hiCFeatures+eQTLFeatures+distanceFeature+chrFeature]
combinations['All + Pos'] = [epigenomicFeatures+hiCFeatures+eQTLFeatures+distanceFeature+chrFeature+posfeature]

random_state = 42

with open(outfile,'w+') as f:
    print(outfile, 'RMSE and R Value for Pre-Finemap and Finemap data', sep='\t', file=f)
    print('Data', 'RF RMSE Slope', 'RF RMSE Z', 'R-value Slope', 'R-Value Z-score', sep='\t', file=f)
    for i in combinations:
        X_label = [item for sublist in combinations[i] for item in sublist]
        random_forest = random_forest_regressor(data_z, X_label, 'Cells_EBV-transformed_lymphocytes', random_state=random_state)
        # decision_tree = decision_tree_regressor(data_z, X_label, 'Cells_EBV-transformed_lymphocytes', random_state=random_state)
        random_forest_z = random_forest_regressor(data_z, X_label, 'z', random_state=random_state)
        # decision_tree_z = decision_tree_regressor(data_z, X_label, 'z', random_state=random_state)
        print("{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}".format(i,random_forest['rmse'],random_forest_z['rmse'], random_forest['r_value'], random_forest_z['r_value']), file=f)
        pd.DataFrame(random_forest_z['y_test']).reset_index(drop=True).join(pd.DataFrame(random_forest_z['y_pred'])).to_csv(outfolder+i+'_Real_vs_pred.csv')
        data_z.loc[random_forest_z['y_test'].index].to_csv(outfolder+i+'_real.csv')
        pd.DataFrame(['y_pred']).to_csv(outfolder+i+'_predicted.csv')
        # data_z.loc[random_forest_z['y_train'].index].to_csv(outfolder+'ytrain.csv')
        importances = random_forest_z['importances']
        feature_list = random_forest_z['features']
        feature_importances = [(feature, round(importance, 2)) for feature, importance in zip(feature_list, importances)]
        feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)
        with open(outfolder+str(i)+'_top_feature_importance.csv','w+') as impfile:
            for pair in feature_importances[:10]: print('Variable:\t{}\tImportance:\t{}'.format(*pair),file=impfile)
    dummy = dummy_regressor(data_z, X_label, 'z', random_state)
    print('Dummy: {:.4f}'.format(dummy['rmse']), file=f)
