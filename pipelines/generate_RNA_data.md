# GTEx RNA preprocessing

For this analysis we use only the top 5000 most expressed genes. That is the 5000 genes that have higher distance.

The raw data we are using is:  ```GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz``` from @Gtex link

Before identifying the top genes we can perform a sanity check.
First, the dimensions of the data:

```{sh}
$ zcat GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz | head -n2
#1.2
56202   11688

$ zcat GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz | wc -l
56205

$ zcat GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz | head -n3 | tail -n1 | cut -f1-5
Name    Description     GTEX-1117F-0226-SM-5GZZ7        GTEX-111CU-1826-SM-5GZYN        GTEX-111FC-0226-SM-5N9B8
```

The data is 56202 transcripts each with data for 11688 samples. It has 3 additional lines at the begining: 
- The version of the file
- The dimension of the data
- The headers of the data

Second, check for duplicate genes:
```{sh}
$ zcat GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz | tail -n56202 | cut -f1 | sort | cut -d"." -f1| uniq -c | sort -rn | head -n1
      1 ENSGR0000266731
```

Now, let's extract the genes with a mean value in the top 5000. 
For the sake of memory this step is split in two. 
1. Go through the data and get the indexes of the transcripts with the mean in the top.
2. Go through the data and extract the information of the transcripts in the top.

```{sh}
$ python scripts/getTopGenesInd.py data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_top_indices 5000
$ python scripts/getTopGenesData.py data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_top_indices.npy data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top.tab
```

Finally, we check that the table looks good:
```{sh}
$ head data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top.tab | cut -f1-5
Name    Description     GTEX-1117F-0226-SM-5GZZ7        GTEX-111CU-1826-SM-5GZYN        GTEX-111FC-0226-SM-5N9B8
ENSG00000225972.1       MTND1P23        10.75   11.69   9.167
ENSG00000225630.1       MTND2P28        356.2   757.9   796.6
ENSG00000237973.1       hsa-mir-6723    32.41   25.3    37.08
ENSG00000229344.1       RP5-857K21.7    13.26   13.98   10.42
ENSG00000248527.1       MTATP6P1        1420    2694    3066
ENSG00000198744.5       RP5-857K21.11   25.01   21.17   14.33
ENSG00000188976.6       NOC2L   90.38   66.14   52
ENSG00000188290.6       HES4    13.66   28.55   30.93
ENSG00000187608.5       ISG15   12.63   39.96   59.52
```
