cd /Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/gentropy

module load igmm/apps/python/3.10.14 
pip install gentropy
export HYDRA_FULL_ERROR=1

mkdir -p variant-index
wget -e robots=off -r -np -nH --cut-dirs=6 -R "index.html*" \
  -P "$WD/variant-index" \
  "https://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/lut/variant-index/"

mkdir -p "$WD/genes-index"
wget -e robots=off -r -np -nH --cut-dirs=6 -R "index.html*" \
  -P "$WD/genes-index" \
  "https://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/lut/genes-index/"

export CS_PATH="$WD/resources/credible_sets.parquet" 

gentropy \
  step=locus_to_gene_feature_matrix \
  step.credible_set_path="$CS_PATH" \
  step.variant_index_path="$WD/variant-index" \
  step.target_index_path="$WD/genes-index" \
  step.feature_matrix_path="$WD/l2g_feature_matrix.parquet"



 