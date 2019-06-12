# Import needed libraries
import sys
import pandas as pd

# List of GTEx eQTL tissues
tissues_names =  ['Adipose_Subcutaneous',  'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial', 'Brain_Amygdala',
                 'Brain_Anterior_cingulate_cortex_BA24', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum',
                 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus', 'Brain_Nucleus_accumbens_basal_ganglia',
                 'Brain_Putamen_basal_ganglia', 'Brain_Spinal_cord_cervical_c-1', 'Brain_Substantia_nigra', 'Breast_Mammary_Tissue', 'Cells_Transformed_fibroblasts',
                 'Colon_Sigmoid', 'Colon_Transverse', 'Esophagus_Gastroesophageal_Junction', 'Esophagus_Mucosa', 'Esophagus_Muscularis', 'Heart_Atrial_Appendage',
                 'Heart_Left_Ventricle', 'Liver', 'Lung', 'Minor_Salivary_Gland', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary', 'Pancreas',
                 'Pituitary', 'Prostate', 'Skin_Not_Sun_Exposed_Suprapubic', 'Skin_Sun_Exposed_Lower_leg',
                 'Small_Intestine_Terminal_Ileum', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus','Vagina', 'Whole_Blood', 'Cells_EBV-transformed_lymphocytes']

# Read the input line
eqtl_folder = sys.argv[1]
in_file = sys.argv[2]
out_file = sys.argv[3]

# For each tissue merge the eQTL with the list of eQTL we are interested in
pairs_tissues = pd.read_csv(in_file)
pairs_tissues['gene_id'] = pairs_tissues['gene_id'].apply(lambda x: x.split('.')[0])
for tissue in tissues_names:
    eqtl_file = eqtl_folder+tissue+'.v7.signif_variant_gene_pairs.txt.gz'
    tissue_data = pd.read_csv(eqtl_file, sep='\t')
    tissue_data['gene_id'] = tissue_data['gene_id'].apply(lambda x: x.split('.')[0])
    tissue_data = tissue_data[['variant_id', 'gene_id', 'slope']]
    pairs_tissues = pairs_tissues.merge(tissue_data, on=['variant_id', 'gene_id'], how='left')
    pairs_tissues[tissue] = pairs_tissues['slope']
    pairs_tissues = pairs_tissues.drop('slope', axis=1)

    
# Save the data
pairs_tissues.to_csv(out_file, compression='gzip', index=False)
