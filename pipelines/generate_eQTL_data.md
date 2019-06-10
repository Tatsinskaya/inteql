# eQTL preprocessing

We are going to only use the pairs that its gene is present in the list of top expressed genes, thus we get the list of pairs that fullfield this criteria.
```{sh}
$ python scripts/getPairsIdTop.py data/GTEx_Analysis_v7_eQTL/ data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm_top.tab data/GTEx_Analysis_v7_eQTL_allTissues_pairs_top.csv.gz
```
This script print the total number of pairs after adding a new tissue (for the last tissue _Cells_EBV-transformed_lymphocytes_ we have a total of  _1349165_ pairs). To sanity check:
```{sh}
$ zcat data/GTEx_Analysis_v7_eQTL_allTissues_pairs_top.csv.gz | wc -l
1349166
```

Then, we get the data for each of the tissues:
```{sh}
$ python scripts/getPairsSlopeTop.py data/GTEx_Analysis_v7_eQTL/ data/GTEx_Analysis_v7_eQTL_allTissues_pairs_top.csv.gz data/GTEx_Analysis_v7_eQTL_allTissues_slope_top.csv.gz
```

A bit of cheking:
```{sh}
$ zcat data/GTEx_Analysis_v7_eQTL_allTissues_slope_top.csv.gz | wc -l
1349166

$ zcat data/GTEx_Analysis_v7_eQTL_allTissues_slope_top.csv.gz | head -n1 | sed 's/[^,]//g' | wc -c
50
```