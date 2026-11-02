
source activate ldsc

LDSC_DIR=./resources/ldsc
HM3_SNPS=./resources/w_hm3.snplist

## Multisite chronic pain
FILE=./resources/summary_statistics/MCP/ChronicPain.stats
OUT_DIR=./resources/summary_statistics/MCP/

python $LDSC_DIR/munge_sumstats.py \
        --sumstats "$FILE" \
        --snp SNP \
        --a1 ALLELE1 \
        --a2 ALLELE0 \
        --N 387649 \
        --out "$OUT_DIR/MCP_munged" \
        --merge-alleles "$HM3_SNPS"


## Depression
FILE=./resources/summary_statistics/MDD/pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv
OUT_DIR=./resources/summary_statistics/MDD/

python $LDSC_DIR/munge_sumstats.py \
        --sumstats "$FILE" \
        --snp ID \
        --a1 EA \
        --a2 NEA \
        --N 1639572 \
        --p PVAL \
        --out "$OUT_DIR/MDD_munged" \
        --merge-alleles "$HM3_SNPS"      

