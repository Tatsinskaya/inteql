# Import needed libraries
import sys
import gzip
import pandas as pd

# list of GTEx eQTL tissues
tissues_names = ['Adipose_Subcutaneous',  'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial', 'Brain_Amygdala', 
                 'Brain_Anterior_cingulate_cortex_BA24', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum', 
                 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus', 'Brain_Nucleus_accumbens_basal_ganglia',
                 'Brain_Putamen_basal_ganglia', 'Brain_Spinal_cord_cervical_c-1', 'Brain_Substantia_nigra', 'Breast_Mammary_Tissue', 'Cells_Transformed_fibroblasts', 
                 'Colon_Sigmoid', 'Colon_Transverse', 'Esophagus_Gastroesophageal_Junction', 'Esophagus_Mucosa', 'Esophagus_Muscularis', 'Heart_Atrial_Appendage',
                 'Heart_Left_Ventricle', 'Liver', 'Lung', 'Minor_Salivary_Gland', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary', 'Pancreas',
                 'Pituitary', 'Prostate', 'Skin_Not_Sun_Exposed_Suprapubic', 'Skin_Sun_Exposed_Lower_leg', 
                 'Small_Intestine_Terminal_Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus','Vagina', 'Whole_Blood', 'Cells_EBV-transformed_lymphocytes']
     
# Read the input
eqtl_folder = sys.argv[1]   # original-data/GTEx_Analysis_v7_eQTL/
exp_file = sys.argv[2]      # output02.tab
out_file = sys.argv[3]      # GTEx_Analysis_v7_eQTL_allTissues_pairs_top.csv.gz = output04.csv.gz

# Get the set of the top genes from the RNA file
top_genes_set = set()
header = True
with open(exp_file, 'r') as f:
    for line in f:
        if header:
            header = False
            continue
        top_genes_set.add(line.split()[0].split('.')[0])
        

# Get the id of the eQTLs for each tissue that is not been previously seen and the gene is in the top 5k
print('Tissue', 'Total number of pairs')
pairs = set()
for tissue in tissues_names:
    eqtl_file = eqtl_folder+tissue+'.v7.signif_variant_gene_pairs.txt.gz'
    header = True
    with gzip.open(eqtl_file) as f:
        for line in f:
            if header:
                header = False
                continue
            line = line.split()
            gene = line[1].split('.'.encode())[0]
            if gene.decode() in top_genes_set:
                pairs.add((line[0].decode(), line[1].decode()))
    print(tissue, len(pairs))

# Save the data
variant_id, gene_id = tuple(zip(*pairs))
pd.DataFrame({'variant_id': variant_id, 'gene_id':gene_id }).to_csv(out_file, index=False, compression='gzip')
