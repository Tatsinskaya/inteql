# Input and output files names: todo Fix Input/Output according to argv
in_file <- '/nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output02.tab'
pcSd_file <- '/nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output03-tpm_pcSd.rds'
pcsX_file <- '/nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output03-tpm_pcsX.rds'
cmdX_file <- '/nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output03-tpm_cmdX.rds'
tsne_file <- '/nfs/research1/zerbino/jhidalgo/inteql/data/inter_data/output03-tpm_tsne.rds'

library("limma")

# Read data
raw <- read.csv(in_file, stringsAsFactors = F, sep='\t')

# Extract values
X <- raw[,3:ncol(raw)]
# Normalize
Y <- limma::voom(X)$E

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
