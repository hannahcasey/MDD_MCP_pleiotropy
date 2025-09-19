library(dplyr)
library(data.table)

# MDD GWAS ----
file <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MDD//pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv"
MDD_GWAS <- fread(file, sep = "\t")
MDD_GWAS$Z <- MDD_GWAS$BETA / MDD_GWAS$SE

MDD_sumstats <- data.frame(
  SNP = MDD_GWAS$ID,
  CHR = as.integer(MDD_GWAS$CHROM),
  BP  = as.integer(MDD_GWAS$POS),
  A1  = toupper(MDD_GWAS$EA),
  A2  = toupper(MDD_GWAS$NEA),
  N   = as.numeric(MDD_GWAS$NEFF),
  Z   = as.numeric(MDD_GWAS$Z),
  P   = MDD_GWAS$PVAL,
  INFO= MDD_GWAS$IMPINFO
)

gz <- gzfile("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MDD/MDD_MiXeR.sumstats.gz", "wt")
write.table(MDD_sumstats, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)

# MCP GWAS ----
MCP_GWAS <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MCP/ChronicPain.stats", header = T)

MCP_GWAS$Z <- MCP_GWAS$BETA / MCP_GWAS$SE

MCP_sumstats <- data.frame(
  SNP = MCP_GWAS$SNP,
  CHR = as.integer(MCP_GWAS$CHR),
  BP  = as.integer(MCP_GWAS$BP),
  A1  = toupper(MCP_GWAS$ALLELE1),
  A2  = toupper(MCP_GWAS$ALLELE0),
  N   = rep(387649, nrow(MCP_GWAS)),
  Z   = as.numeric(MCP_GWAS$Z),
  P   = MCP_GWAS$P,
  INFO= MCP_GWAS$INFO
)

gz <- gzfile("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MCP/MCP_MiXeR.sumstats.gz", "wt")
write.table(MCP_sumstats, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)
