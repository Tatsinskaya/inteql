# Input and output files names:
in_file = 'data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_normalized_top.csv'
pcSd_file = 'objects/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top_pcSd.rds'
pcsX_file = 'objects/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top_pcsX.rds'
cmdX_file = 'objects/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top_cmdX.rds'
tsne_file = 'objects/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top_tsne.rds'


# Read data
Y = read.csv(in_file)

# Perform PCA
pca <- prcomp(t(Y))

# Save objects
saveRDS(pca$sdv, pcSd_file)
saveRDS(pca$x[,1:2], pcsX_file)

# CMD
# Extract the distance matrix
D <- dist(t(Y))
# Perform the cmd
cmd <- cmdscale(D)
# Save object
saveRDS(cmd, cmdX_file)

# TSNE
library(Rtsne)

# Perform tsne
tsne <- Rtsne::Rtsne(t(Y), perplexity=15)

# Save object
saveRDS(tsne$Y, tsne_file)
