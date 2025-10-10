# cd /exports/eddie/scratch/s2093301/V2G/scores

# wget -r -np -nH --cut-dirs=5 \
#      -A '*.parquet' -e robots=off --reject 'index.html*' \
#      -N --show-progress --continue \
#      https://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/v2g_scored/

# cd ../bin
# wget -q https://github.com/duckdb/duckdb/releases/latest/download/duckdb_cli-linux-amd64.zip
# unzip -o duckdb_cli-linux-amd64.zip
# chmod +x duckdb
# export PATH="$HOME/bin:$PATH"

cd /exports/eddie/scratch/s2093301/V2G/

DDB=bin/duckdb
TEMP_DIR=/exports/eddie/scratch/s2093301/V2G/duckdb_tmp
mkdir -p $TEMP_DIR

OUT="v2g_top_gene_per_variant.csv"
rm -f "$OUT"
first=1

for CHR in {1..22} X Y; do
  tmp="chunk_${CHR}.csv"
  $DDB -c "
  PRAGMA threads=1;
  PRAGMA temp_directory='$TEMP_DIR';
  COPY (
    SELECT
      chr_id,
      position,
      ref_allele,
      alt_allele,
      arg_max(gene_id, overall_score) AS gene_id,
      max(overall_score)               AS overall_score
    FROM read_parquet('scores/*.parquet')
    WHERE chr_id='${CHR}'
    GROUP BY chr_id, position, ref_allele, alt_allele
  ) TO '${tmp}' WITH (HEADER, DELIM ',');
  "
  if [ $first -eq 1 ]; then
    cat "$tmp" > "$OUT"; first=0
  else
    tail -n +2 "$tmp" >> "$OUT"
  fi
  rm -f "$tmp"
done