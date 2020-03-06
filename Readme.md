# INTEQL

This fork is based on Laura Aviñó's Final Bachelor's Degree Project, updated to be automatised on an LSF platform.

The pipeline is comprised of 8 steps:
1) The top N (by default = 5K) top-expressed genes are identified across GTEx RNA-Seq Data (Gene TPMs) and artefact-prone regions are masked.
2) Said RNA-Seq Data is filtered for top-expressed genes.
3) Using GTEx tissue-specific eQTL data, variant-gene associations for the top-expressed genes are filtered.
4) eQTL Slope value for each tissue is matched with each variant-gene pairs.
5) Variant-gene pairs with eQTL data are matched with linear data such as enhancers, promoters, methylation, CAGE (from TargetFinder), the distance between elements and ENSEMBL Regulatory Build elements.
6) HiC data is added for each pair.
7) eQTL data is finemapped to take into account LD.
8) The whole dataset is split and combined into different features to feed a Random Forest model to obtain RMSE and R-value scores to compare each feature or combination of features. GM12878 (EBV transformed lymphocytes) eQTLs are used to measure the performance of the model and the rest of the tissues as a training/test model in a 70/30 split.

eQTL data imputation script was also developed but it is currently not used in the process.

# HOW TO USE

All the data required can be downloaded using data/data_download.sh, which can take a long time due to the size of the files.
Use pipelines/pipeline-inteql.sh as an example of how to run the pipeline on an LSF platform, adjusting any desired path or value.

Increasing the number of top expressed genes used may dramatically affect memory requirements. 
