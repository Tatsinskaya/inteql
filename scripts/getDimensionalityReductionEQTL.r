# Input and output files names:
in_file = '../data/GTEx_Analysis_v7_eQTL_allTissues_slope_top_complete.csv.gz'
pcSd_file = '../objects/GTEx_Analysis_v7_eQTL_allTissues_slope_top_complete_pcSd.rds'
pcsX_file = '../objects/GTEx_Analysis_v7_eQTL_allTissues_slope_top_complete_pcsX.rds'
tsne_file = '../objects/GTEx_Analysis_v7_eQTL_allTissues_slope_top_complete_tsne.rds'

# Read data
raw <- read.csv(in_file, stringsAsFactors = F)

# Extract values
X <- raw[,3:ncol(raw)]

# Perform PCA
pca <- prcomp(t(X))

# Save objects
saveRDS(pca$sdv, pcSd_file)
saveRDS(pca$x[,1:2], pcsX_file)

# TSNE
library(Rtsne)

# Perform tsne
tsne <- Rtsne::Rtsne(t(X), perplexity=15)

# Save object
saveRDS(tsne$Y, tsne_file)                                                                                                                                                                                                                                                                                                                             # Perform tsne                                                                                                                                                              tsne <- Rtsne::Rtsne(t(Y), perplexity=15)                                                                                                                                                                                                                                                                                                               # Save object                                                                                                                                                               saveRDS(tsne$Y, tsne_file)  


