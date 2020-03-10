#!/bin/bash
chromosomes="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X"

cd "original-data" || exit

#########################
### GTEx RNA-Seq TPMs ###
#########################
wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz

##############################
### GTEx Association files ###
##############################
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL.tar.gz
mkdir GTEx_Analysis_v7_eQTL
tar xvzf GTEx_Analysis_v7_eQTL.tar.gz
rm GTEx_Analysis_v7_eQTL.tar.gz

############################
### Rao et al. HiC files ###
############################
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_combined_intrachromosomal_contact_matrices.tar.gz
tar xvzf GSE63525_GM12878_combined_intrachromosomal_contact_matrices.tar.gz


##########################
### GTEx Variant files ###
##########################
mkdir dbSNP
cd dbSNP || exit
mkdir chr-header
wget https://storage.googleapis.com/gtex_analysis_v7/reference/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz
for i in $chromosomes; do
  zcat GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz | head -n1 > "chr-header/GTEx_2016_dbSNPs_table_chr_"$i".txt"
  zcat GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt.gz | awk -F $'\t' -v CHR="$i" '{if($1 == CHR) print $0}' >> "chr-header/GTEx_2016_dbSNPs_table_chr_"$i".txt"
done
cd ..

################################
### Ensembl Regulatory build ###
################################
wget ftp://ftp.ensembl.org/pub/grch37/release-99/regulation/homo_sapiens/RegulatoryFeatureActivity/GM12878/homo_sapiens.GRCh37.GM12878.Regulatory_Build.regulatory_activity.20180925.gff.gz
gzip -d homo_sapiens.GRCh37.GM12878.Regulatory_Build.regulatory_activity.20180925.gff.gz
rm homo_sapiens.GRCh37.GM12878.Regulatory_Build.regulatory_activity.20180925.gff.gz
