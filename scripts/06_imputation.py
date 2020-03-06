import numpy as np
import pandas as pd
from scipy.stats import linregress
from sklearn.metrics import mean_squared_error

eqtl_file = '../data/output/output05.csv.gz' # GTEx_Analysis_v7_eQTL_allTissues_slope_top.csv.gz

# Read eqtl
eqtl = pd.read_csv(eqtl_file)
eqtl = eqtl.drop(columns=['gene_id', 'variant_id'])

# Get labels
y_label = 'Cells_EBV-transformed_lymphocytes'
X_label = list(i for i in eqtl.columns if i != y_label)
# Get the data the eQTL for that tissue
data = eqtl[np.invert(pd.isna(eqtl[y_label]))]
# Get the eQTL that are not specific
data = data[np.invert(pd.isna(data[X_label])).sum(axis=1) != 0]
mean = data[X_label].mean(axis=1)
# mean[:3]

howMany = np.invert(pd.isna(data[X_label])).sum(axis=1)
#print(howMany[:3])

slope, intercept, rvalue, pvalue, stderr = linregress(data[y_label], mean)
#print(pvalue < 10E6)
print('The rvalue is %f and the pvalue %f' % (rvalue, pvalue))

np.sqrt(mean_squared_error(data[y_label], mean))




# Second Part
#import os
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

pandas2ri.activate()
readRDS = robjects.r['readRDS']

from sklearn.metrics import pairwise_distances, mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor

from scipy.stats import linregress

tissues_id = []
with open('/nfs/research1/zerbino/jhidalgo/inteql/data/original-data/Annotations/GTEx_v7_Annotations_TissuesId.txt') as f:
    for line in f:
        tissues_id.append(line.strip())
# Weights by PCA
pc = readRDS('../objects/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top_pcsX.rds')
pc = pandas2ri.ri2py(pc)
pc = pd.DataFrame(pc, columns=['PC1', 'PC2'])
pc['tissue'] = tissues_id
pc.head()
centroids_pc = pc.groupby('tissue')['PC1', 'PC2'].apply(np.mean)
centroids_pc.shape
df_pc = pd.DataFrame(1 / pairwise_distances(centroids_pc), columns=centroids_pc.index,
                     index=centroids_pc.index).replace(np.inf, 0)
# Weights by Tsne
ts = readRDS('../objects/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top_tsne.rds')
ts = pandas2ri.ri2py(ts)
ts = pd.DataFrame(ts, columns=['x', 'y'])
ts['tissue'] = tissues_id
ts.head()
centroids_ts = ts.groupby('tissue')['x', 'y'].apply(np.mean)
centroids_ts.shape
df_ts = pd.DataFrame(1 / pairwise_distances(centroids_ts), columns=centroids_ts.index,
                     index=centroids_ts.index).replace(np.inf, 0)
eqtl_s = eqtl[np.invert(pd.isna(eqtl)).sum(axis=1) != 1]
print(eqtl_s.shape)
eqtl.head()
f_e = open('../data/Plotting/data_RMSE_figure3.tab', 'w')
f_r = open('../data/Plotting/data_rval_figure3.tab', 'w')
print("Tissue\tn_eQTL\tRF\tmean\tpca\ttsne", file=f_e)
print("Tissue\tn_eQTL\tRF\tmean\tpca\ttsne", file=f_r)

for tissue in eqtl_s.columns:

    print(tissue, end='\t')
    y_label = tissue
    W = df_pc[y_label]
    U = df_ts[y_label]

    X_label = list(data.columns)
    X_label.remove(y_label)

    data = eqtl_s[np.invert(pd.isna(eqtl_s[y_label]))]
    print(data.shape)
    data = data[np.invert(pd.isna(data[X_label])).sum(axis=1) != 1]
    print(data.shape)

    # Split data
    X = data[X_label]
    y = data[y_label]
    X_train, X_test, y_train, y_test = train_test_split(X.fillna(0), y.fillna(0), test_size=0.3, random_state=12)

    print('Model...')
    # Model
    model = RandomForestRegressor()
    model = model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    # Regression
    x, y = y_test, y_pred
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    r_model = r_value
    e_model = np.sqrt(mean_squared_error(x, y))

    print('Mean...', end='\t')
    # Normal mean
    mean = X[X_label].mean(axis=1).to_list()

    # Regression
    x, y = data[y_label], mean
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    r_mean = r_value
    e_mean = np.sqrt(mean_squared_error(x, y))

    print('PCA...', end='\t')
    # PCA
    weighted_mean = []
    for i, r in X.iterrows():
        w = W[r.keys()[np.invert(pd.isna(r))]]
        w = w / w.sum()
        m = 0
        for t, v in r.iteritems():
            if not pd.isna(v):
                m += w[t] * v
        weighted_mean.append(m)

    # Regression
    x, y = data[y_label].values, weighted_mean
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    r_pca = r_value
    e_pca = np.sqrt(mean_squared_error(x, y))

    print('tSNE...')
    # tSNE
    weighted_mean = []
    for i, r in X.iterrows():

        u = U[r.keys()[np.invert(pd.isna(r))]]
        u = u / u.sum()
        m = 0
        for t, v in r.iteritems():
            if not pd.isna(v):
                m += u[t] * v
        weighted_mean.append(m)

    # Regression
    x, y = data[y_label].values, weighted_mean
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    r_tsne = r_value
    e_tsne = np.sqrt(mean_squared_error(x, y))

    print("{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}".format(y_label, data.shape[0], e_model, e_mean, e_pca, e_tsne),
          file=f_e)
    print("{}\t{}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}".format(y_label, data.shape[0], r_model, r_mean, r_pca, r_tsne),
          file=f_r)

f_e.close()
f_r.close()
