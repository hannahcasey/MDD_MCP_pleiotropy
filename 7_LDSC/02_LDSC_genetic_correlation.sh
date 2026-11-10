
source activate ldsc

LDSC_DIR=./resources/ldsc/
LDSCORE_DIR=./resources/ldsc/eur_w_ld_chr/
OUT_DIR=./output/rg/


python "$LDSC_DIR/ldsc.py" \
  --rg ./resources/summary_statistics/MCP/MCP_munged.sumstats.gz,./resources/summary_statistics/MDD/MDD_munged.sumstats.gz\
  --ref-ld-chr "$LDSCORE_DIR" \
  --w-ld-chr "$LDSCORE_DIR" \
  --out ${OUT_DIR}/MCP_vs_MDD_rg

