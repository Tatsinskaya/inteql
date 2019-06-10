To perform a weighted mean we need to get some coeficients that represent the "distance" between tissues. One way to do so is using the RNA expression of tissues. Before using the python scripts @pythonscriptFilter we got the top 5000 most expressed genes and the expression for each sample in GTEx. Then we can run the following:

```{sh}
$ Rscript scripts/getDimensionalityReduction.r
```

In case you have already the data normalized run:
```{sh}
$ Rscript scripts/getDimensionalityReductionY.r
```

This script reads the data and normalize it using log CPM normalization with the function limma::voom().
After that, it performs dimensionality reduction using PCA (`prcomp()`), CMD (`cmdscale()`) and tnse (`Rtsne::Rtsne()`).
Finally, it saves the objects of interest.

The list of files involved:
```{sh}
in_file = '../data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top.tab'
pcSd_file = '../objects/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top_pcSd.rds'
pcsX_file = '../objects/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top_pcsX.rds'
cmdX_file = '../objects/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top_cmdX.rds'
tsne_file = '../objects/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top_tsne.rds'
```

To visualize the results we can use the following R notebook @@.
To extract the distances and coefficients for the we can use the @@ jupyter notebook


To get low dimensionality eQTL space we can run:
```{sh}
$ Rscript getDimensionalityReductionEQTL.r
```

After we can visualize the results using @@.


# REFINE PCA:

Using my notebook I will perform limma::voom().
Then I will run, and PCA() as well.

