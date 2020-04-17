from __future__ import print_function
import pandas as pd
import sys

sys.path.append('../inteql/')
from modelFunctions import *

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error
from scipy.stats import linregress


def random_forest_regressor_fixed_split(traindata_model, testdata_model, X_label, y_label, max_depth=83, n_estimators=136,
                                        min_samples_split=2, min_samples_leaf=1, max_features='sqrt'):
    # Split data
    X_train = traindata_model[X_label]
    y_train = traindata_model[y_label]
    X_test = testdata_model[X_label]
    y_test = testdata_model[y_label]
    # Create a model and train it
    model = RandomForestRegressor(random_state=42, max_depth=max_depth, n_estimators=n_estimators)
    model = model.fit(X_train, y_train)
    # Predict values for test set and asses error
    y_pred = model.predict(X_test)
    rmse = np.sqrt(mean_squared_error(y_test, y_pred))
    slope, intercept, r_value, p_value, std_err = linregress(y_test, y_pred)
    importances = list(model.feature_importances_)
    feature_list = list(traindata_model[X_label].columns)
    # Return result object
    return {'rmse': rmse, 'model': model, 'y_pred': y_pred, 'r_value': r_value, 'y_test': y_test,
            'importances': importances, 'features': feature_list}


finemap_file = sys.argv[1]  # '/nfs/research1/zerbino/jhidalgo/inteql/data/output/output09.csv.gz'
eqtl_file = sys.argv[2]  # '/nfs/research1/zerbino/jhidalgo/inteql/data/output/output05.csv.gz'
outfile = sys.argv[3]  # /nfs/research1/zerbino/jhidalgo/inteql/data/output/output10_chr.csv.gz
outfolder = sys.argv[4]

# Finemap data
data_z = pd.read_csv(finemap_file)
eqtl = pd.read_csv(eqtl_file)
data_z = data_z.merge(eqtl, on=['variant_id', 'gene_id']).drop_duplicates()
data_z = data_z.fillna(0)
data_z['z'] = data_z['z'].astype(float)

temptestdata = pd.DataFrame()
temptraindata = pd.DataFrame()
testdata = pd.DataFrame()
traindata = pd.DataFrame()
for i in pd.unique(data_z['chr']):
    limit = np.percentile(np.array(data_z[data_z['chr'] == i]['variant_pos']), 33)
    temptestdata = pd.DataFrame(data_z[data_z['chr'] == i][data_z[data_z['chr'] == i]['variant_pos'] >= limit])
    temptraindata = pd.DataFrame(data_z[data_z['chr'] == i][data_z[data_z['chr'] == i]['variant_pos'] < limit])
    frames1 = [temptestdata, testdata]
    frames2 = [temptraindata, traindata]
    testdata = pd.concat(frames1)
    traindata = pd.concat(frames2)

# Get list of features names for each subset
epigenomicFeatures = list(i for i in data_z.columns if i.startswith('enha') or i.startswith('prom'))[2:]
regbuild_typeFeatures = list(i for i in data_z.columns if i.startswith('type'))
regbuild_activityFeatures = list(i for i in data_z.columns if i.startswith('activity'))
regbuild_dataFeature = ['RegElementInfo']
hiCFeatures = list(i for i in data_z.columns if i.startswith('hi'))
eQTLFeatures = list(eqtl.columns[2:-1])
distanceFeature = ['var_prom_distance', 'var_enh_distance', 'tss_distance']
chrFeature = ['chr']
posfeature = ['variant_pos']

combinations = {}
combinations['Epigenomics'] = [epigenomicFeatures]
# combinations['RegBuild Type'] = [regbuild_typeFeatures]
# combinations['RegBuild Act'] = [regbuild_activityFeatures]
# combinations['RegBuild Data'] = [regbuild_dataFeature]
# combinations['All Reg'] = [regbuild_dataFeature+regbuild_activityFeatures+regbuild_typeFeatures]
# combinations['All Reg + eQTL'] = [regbuild_dataFeature+regbuild_activityFeatures+regbuild_typeFeatures+eQTLFeatures]
# combinations['All Reg + Epi'] = [regbuild_dataFeature+regbuild_activityFeatures+regbuild_typeFeatures+epigenomicFeatures]
combinations['HiC'] = [hiCFeatures]
combinations['eQTL'] = [eQTLFeatures]
# combinations['Distance'] = [distanceFeature]
# combinations['Distance+Chromosome'] = [distanceFeature+chrFeature]
# combinations['HiC + Epigenomic'] = [hiCFeatures+epigenomicFeatures]
combinations['Epi + eQTL'] = [epigenomicFeatures + eQTLFeatures]
combinations['eQTL + Distance'] = [eQTLFeatures + distanceFeature]
# combinations['Epi + eQTL + RegBuild Type'] = [epigenomicFeatures+distanceFeature+regbuild_typeFeatures]
# combinations['Epi + eQTL + RegBuild Act'] = [epigenomicFeatures+distanceFeature+regbuild_activityFeatures]
# combinations['Epi + eQTL + RegBuild Data'] = [epigenomicFeatures+distanceFeature+regbuild_dataFeature]
combinations['All'] = [epigenomicFeatures + hiCFeatures + eQTLFeatures + distanceFeature + chrFeature]
combinations['All + Pos'] = [
    epigenomicFeatures + hiCFeatures + eQTLFeatures + distanceFeature + chrFeature + posfeature]

with open(outfile, 'w+') as f:
    print(outfile, 'RMSE and R Value for Pre-Finemap and Finemap data', sep='\t', file=f)
    print('Data', 'RF RMSE Slope', 'RF RMSE Z', 'R-value Slope', 'R-Value Z-score', sep='\t', file=f)
    for i in combinations:
        X_label = [item for sublist in combinations[i] for item in sublist]
        random_forest = random_forest_regressor_fixed_split(traindata, testdata, X_label,
                                                            'Cells_EBV-transformed_lymphocytes')
        random_forest_z = random_forest_regressor_fixed_split(traindata, testdata, X_label, 'z')
        print("{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}".format(i, random_forest['rmse'], random_forest_z['rmse'],
                                                          random_forest['r_value'], random_forest_z['r_value']), file=f)
        # Z output
        pd.DataFrame(random_forest_z['y_test']).reset_index(drop=True).join(
            pd.DataFrame(random_forest_z['y_pred'])).to_csv(outfolder + i + '_Real_vs_pred_Z.csv')
        data_z.loc[random_forest_z['y_test'].index].to_csv(outfolder + i + '_real_Z.csv')
        pd.DataFrame(random_forest_z['y_pred']).to_csv(outfolder + i + '_predicted_Z.csv')

        importances = random_forest_z['importances']
        feature_list = random_forest_z['features']
        feature_importances = [(feature, round(importance, 2)) for feature, importance in
                               zip(feature_list, importances)]
        feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)
        with open(outfolder + str(i) + '_top_feature_importance_Z.csv', 'a') as impfile:
            for pair in feature_importances[:10]:
                print('Variable:\t{}\tImportance:\t{}'.format(*pair), file=impfile)
        # Raw score output
        pd.DataFrame(random_forest['y_test']).reset_index(drop=True).join(
            pd.DataFrame(random_forest['y_pred'])).to_csv(outfolder + i + '_Real_vs_pred_raw.csv')
        data_z.loc[random_forest['y_test'].index].to_csv(outfolder + i + '_real_raw.csv')
        pd.DataFrame(random_forest['y_pred']).to_csv(outfolder + i + '_predicted_raw.csv')
        importances = random_forest_z['importances']
        feature_list = random_forest_z['features']
        feature_importances = [(feature, round(importance, 2)) for feature, importance in
                               zip(feature_list, importances)]
        feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)
        with open(outfolder + str(i) + '_top_feature_importance_raw.csv', 'a') as impfile:
            for pair in feature_importances[:10]:
                print('Variable:\t{}\tImportance:\t{}'.format(*pair), file=impfile)
    dummy = dummy_regressor(data_z, X_label, 'z', 42)
    print('Dummy: {:.4f}'.format(dummy['rmse']), file=f)
