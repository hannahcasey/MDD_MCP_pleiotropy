cd /Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/

chromosomes=(5 13 7 11 2 9 4 3 20 10 18) 
for c in $chromosomes; do
  gcloud storage cp \
    --billing-project=united-button-442820-q5 \
    "gs://gtex-resources/GTEx_Analysis_v10_QTLs/GTEx_Analysis_v10_eQTL_all_associations/Brain_*.v8.EUR.allpairs.chr${c}.parquet" \
    summary_statistics/eQTL/brain/full/
done

awk -F, 'NR>1{
  g=$NF; gsub(/"/,"",g);                 # drop quotes
  n=split(g,a,/;/);                      # split on ;
  for(i=1;i<=n;i++) if(a[i]!="") print a[i];
}' ../output/V2G/cs_with_v2g_genes.csv | sort -u > genes.txt

typeset -a GENES
GENES=("${(@f)$(< genes.txt)}")

PARQ_DIR="summary_statistics/eQTL/brain/full"
OUT_DIR="summary_statistics/eQTL/brain/by_gene"
mkdir -p "$OUT_DIR"

for gene in "${GENES[@]}"; do

  out_file="${OUT_DIR}/${gene}.v8.EUR.brain.allchr.csv"
  
duckdb -c "
    COPY (
      SELECT
        t.*,
        regexp_extract(
          filename,
          '.*/([^/]+)\\.v8\\.EUR\\.allpairs\\.chr[0-9]+\\.parquet',
          1
        ) AS tissue
      FROM read_parquet('${PARQ_DIR}/*.v8.EUR.allpairs.chr*.parquet', filename=true) AS t
      WHERE starts_with(t.phenotype_id, '${gene}.')
    ) TO '${out_file}' (FORMAT CSV, HEADER, DELIMITER ',');
  "
done