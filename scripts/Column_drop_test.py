from __future__ import print_function
from sklearn.base import clone
import pandas as pd
import sys
import numpy as np
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.dummy import DummyRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from scipy.stats import linregress
#
# sys.path.append('../inteql/')
# from modelFunctions import *


def imp_df(column_names, importances):
    df = pd.DataFrame({'feature': column_names,
                       'feature_importance': importances}) \
        .sort_values('feature_importance', ascending=False) \
        .reset_index(drop=True)
    return df


def drop_col_feat_imp(model, X_train, y_train, random_state=42):
    # clone the model to have the exact same specification as the one initially trained
    model_clone = clone(model)
    # set random_state for comparability
    model_clone.random_state = random_state
    # training and scoring the benchmark model
    model_clone.fit(X_train, y_train)
    benchmark_score = model_clone.score(X_train, y_train)
    print('Benchmark score: '+benchmark_score)
    # list for storing feature importances
    importances = []
    # iterating over all columns and storing feature importance (difference between benchmark and new model)
    for col in X_train.columns:
        model_clone = clone(model)
        model_clone.random_state = random_state
        model_clone.fit(X_train.drop(col, axis=1), y_train)
        drop_col_score = model_clone.score(X_train.drop(col, axis=1), y_train)
        importances.append(benchmark_score - drop_col_score)
    importances_df = imp_df(X_train.columns, importances)
    print(importances_df)
    return importances_df


finemap_file = sys.argv[1]  # '/nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output09.csv.gz'
eqtl_file = sys.argv[2]  # '/nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output05.csv.gz'
outfile = sys.argv[3]  # /nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output10_chr.csv.gz

# Finemap data
data_z = pd.read_csv(finemap_file)
eqtl = pd.read_csv(eqtl_file)
data_z = data_z.merge(eqtl, on=['variant_id', 'gene_id']).drop_duplicates()
data_z = data_z.fillna(0)
data_z['z'] = data_z['z'].astype(float)
data_z['Chromosome'] = data_z['variant_id'].str.split("_", 1, expand=True)[0]

np.random.seed(seed=42)
data_z['random'] = np.random.random(size=len(data_z))

# Get list of features names for each subset
epigenomicFeatures = list(i for i in data_z.columns if i.startswith('enha') or i.startswith('prom'))[2:]
regbuild_typeFeatures = list(i for i in data_z.columns if i.startswith('type'))
regbuild_activityFeatures = list(i for i in data_z.columns if i.startswith('activity'))
regbuild_dataFeature = ['RegElementInfo']
hiCFeatures = list(i for i in data_z.columns if i.startswith('hi'))
eQTLFeatures = list(eqtl.columns[2:-1])
distanceFeature = ['var_prom_distance', 'var_enh_distance']
chrFeature = ['Chromosome']
randomFeature = ['random']

X_label = epigenomicFeatures + hiCFeatures + eQTLFeatures + distanceFeature + randomFeature + chrFeature
X = data_z[X_label]
y = data_z['z']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)
# Create a model and train it
model = RandomForestRegressor(random_state=42, max_depth=None, n_estimators=100)
model = model.fit(X_train, y_train)
# Predict values for test set and asses error
y_pred = model.predict(X_test)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))
slope, intercept, r_value, p_value, std_err = linregress(y_test, y_pred)
feature_list = list(X.columns)
importances = list(model.feature_importances_)
feature_importances = [(feature, round(importance, 2)) for feature, importance in zip(feature_list, importances)]
feature_importances = sorted(feature_importances, key=lambda x: x[1], reverse=True)

# with open(outfile, 'w+') as f:
[print('Variable: {:20} Importance: {}'.format(*pair)) for pair in feature_importances]

drop = pd.DataFrame(drop_col_feat_imp(model,X_train,y_train))
drop.to_csv(outfile)
