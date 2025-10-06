cd /exports/eddie/scratch/s2093301/V2G/scores

wget -r -np -nH --cut-dirs=5 \
     -A '*.parquet' -e robots=off --reject 'index.html*' \
     -N --show-progress --continue \
     https://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/v2g_scored/

mkdir ../bin 
cd ../bin

wget -q https://github.com/duckdb/duckdb/releases/latest/download/duckdb_cli-linux-amd64.zip
unzip -o duckdb_cli-linux-amd64.zip
chmod +x duckdb
export PATH="$HOME/bin:$PATH"


bin/duckdb -c "
PRAGMA threads=4;
COPY (
  SELECT
    chr_id,
    position,
    gene_id,
    overall_score
  FROM read_parquet('scores/*.parquet')   -- use '**/*.parquet' if files are in nested subfolders
) TO 'v2g_scored_subset.csv' WITH (HEADER, DELIM ',');
"

