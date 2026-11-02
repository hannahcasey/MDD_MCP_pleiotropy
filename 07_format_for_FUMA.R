library(data.table)

# Load GWAS
MDD_GWAS <- fread("./resources/summary_statistics/MDD//pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv", sep = "\t")
MCP_GWAS <- read.table("./resources/summary_statistics/MCP/ChronicPain.stats", header = T)

MDD_fuma <- data.frame(
  ID = MDD_GWAS$ID,
  P = MDD_GWAS$PVAL
)

MCP_fuma <- data.frame(
  ID = MCP_GWAS$SNP,
  P = MCP_GWAS$P
)

write.table(MDD_fuma, "./resources/FUMA/MDD_fuma.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(MCP_fuma, "./resources/FUMA/MCP_fuma.txt", sep = "\t", row.names = FALSE, quote = FALSE)

