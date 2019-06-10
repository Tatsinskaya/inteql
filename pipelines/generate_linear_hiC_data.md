# Adding tissue specific data

To add the lienar data we use the notebook: ```addLinear.ipynb```

After, to add Hi-C data we just run:
```{sh}
bsub -M 10000 -R "rusage[mem=10000]" -n 8 python scripts/addHiC.py data/GM12878_replicate/5kb_resolution_intrachromosomal/ data/GTEx_Analysis_v7_eQTL_EVB_linearData.csv.gz data/GTEx_Analysis_v7_eQTL_EVB_linearData_hiC.csv.gz
```