library(ggplot2)
library(locuszoomr)
library(patchwork)
library(data.table)
library(ggforce)
#BiocManager::install("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)
library(cowplot)
library(dplyr)

# Coloc results ----
eQTL_coloc_results <- readRDS("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/eQTL.rds")

# MDD ----
file <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MDD/pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv"
MDD_GWAS <- fread(file, sep = "\t")
MDD_GWAS$EAF <- ((MDD_GWAS$FCAS * MDD_GWAS$NCAS) + (MDD_GWAS$FCON * MDD_GWAS$NCON)) / (MDD_GWAS$NCAS + MDD_GWAS$NCON)
MDD_n <- 1639572
MDD_s <- 357636 /	(357636 + 1281936)

# MCP ----
MCP_GWAS <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MCP/ChronicPain.stats", header = T)
setDT(MCP_GWAS)
MCP_n <- 387649

# eQTL ----
dir <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/eQTL/brain/v8/by_gene"
files <- list.files(
  dir, pattern = "*.csv", full.names = TRUE
)

gene_names <- str_replace(basename(files), "\\.v8\\.EUR\\.brain\\.allchr\\.csv$", "")
eqtl_list <- set_names(
  map(files, ~ read.csv(.x)),
  gene_names
)
eqtl_list <- map(eqtl_list, ~ dplyr::filter(.x, !is.na(slope)))
ch <- import.chain("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/hg38ToHg19.over.chain")

lift_eqtl_df <- function(df, chain = ch, id_col = "variant_id") {
  if (nrow(df) == 0) return(df)
  if (id_col %in% names(df)) {
    m <- str_match(df[[id_col]], "^chr([0-9XYMT]+)_([0-9]+)_[ACGT]+_[ACGT]+_b38$")
    chr38 <- paste0("chr", m[,2])
    pos38 <- as.integer(m[,3])
  } else if (all(c("chr","pos") %in% names(df))) {
    chr38 <- ifelse(grepl("^chr", df$chr), df$chr, paste0("chr", df$chr))
    pos38 <- as.integer(df$pos)
  }
  
  gr38  <- GRanges(chr38, IRanges(pos38, pos38))
  hits  <- liftOver(gr38, chain)
  lens  <- elementNROWS(hits)
  ok    <- lens == 1L
  
  df$chr_hg38 <- chr38
  df$pos_hg38 <- pos38
  df$chr_hg19 <- NA_character_
  df$pos_hg19 <- NA_integer_
  df$liftover_status <- ifelse(lens==0, "unmapped", ifelse(lens>1, "multi", "ok"))
  
  mapped <- unlist(hits[ok])
  df$chr_hg19[ok] <- as.character(seqnames(mapped))
  df$pos_hg19[ok] <- start(mapped)
  
  df
}
eqtl_list <- map(eqtl_list, lift_eqtl_df)

# Locuszooms ----
## PPP6C Brain_Caudate_basal_ganglia ----
### MDD ----
R    <- as.matrix(fread(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MDD.ld.gz"), header = FALSE))
snps <- readLines(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MDD.snplist"))
rownames(R) <- colnames(R) <- snps
R <- (R + t(R))/2
diag(R) <- 1

snps <- eQTL_coloc_results$Brain_Caudate_basal_ganglia_PPP6C_MDD$meta$snps
R <- R[snps, snps]

eQTL_GWAS <- eqtl_list[["ENSG00000119414"]]
eQTL_GWAS <- subset(eQTL_GWAS, tissue == "Brain_Caudate_basal_ganglia")
eQTL_GWAS$chr_hg19 <- as.numeric(gsub("chr", "", eQTL_GWAS$chr_hg19))
eQTL_GWAS$A1 <- do.call(rbind, strsplit(eQTL_GWAS$variant_id, "_", fixed = TRUE))[, 4]
rsid_map <- dplyr::select(MDD_GWAS, CHROM, POS, ID)
eQTL_GWAS <- left_join(eQTL_GWAS, rsid_map, by = c("chr_hg19" = "CHROM", "pos_hg19" = "POS"))
setDT(eQTL_GWAS)
m_eqtl <- eQTL_GWAS[match(snps, eQTL_GWAS$ID),  .(SNP = ID, A1, BETA = slope, SE = slope_se, EAF = maf, P = pval_nominal)]
m_disorder <- MDD_GWAS[match(snps, MDD_GWAS$ID),  .(ID, EA, BETA, SE, PVAL, POS)]

lead_SNPs_disorder <- eQTL_coloc_results$Brain_Caudate_basal_ganglia_PPP6C_MDD$coloc$summary$hit1[eQTL_coloc_results$Brain_Caudate_basal_ganglia_PPP6C_MDD$coloc$summary$PP.H4.abf > 0.8]
lead_SNPs_eQTL <- eQTL_coloc_results$Brain_Caudate_basal_ganglia_PPP6C_MDD$coloc$summary$hit2[eQTL_coloc_results$Brain_Caudate_basal_ganglia_PPP6C_MDD$coloc$summary$PP.H4.abf > 0.8]

r2_vec <- as.numeric(R[which(snps == lead_SNPs_disorder), ]^2)
dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$POS, `–log10(P)` = -log10(m_disorder$PVAL), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$POS),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_disorder, LD = "r2"
)

p1 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("MDD (", lead_SNPs_disorder[1], ")"))

p_genes <- locuszoomr::gg_genetracks(loc)

dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$POS, `–log10(P)` = -log10(m_eqtl$P), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$POS),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_eQTL, LD = "r2"
)

p2 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("PPP6C eQTL Brain_Caudate_basal_ganglia (", lead_SNPs_eQTL[1], ")"))

p <- p1/p2/p_genes

ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/plots/MDD_PPP6C_eQTL_Brain_Caudate_basal_ganglia.jpg", p, width = 8, height = 12, dpi = 600)

### MCP ----
R    <- as.matrix(fread(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MCP.ld.gz"), header = FALSE))
snps <- readLines(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MCP.snplist"))
rownames(R) <- colnames(R) <- snps
R <- (R + t(R))/2
diag(R) <- 1

snps <- eQTL_coloc_results$Brain_Caudate_basal_ganglia_PPP6C_MCP$meta$snps
R <- R[snps, snps]

eQTL_GWAS <- eqtl_list[["ENSG00000119414"]]
eQTL_GWAS <- subset(eQTL_GWAS, tissue == "Brain_Caudate_basal_ganglia")
eQTL_GWAS$chr_hg19 <- as.numeric(gsub("chr", "", eQTL_GWAS$chr_hg19))
eQTL_GWAS$A1 <- do.call(rbind, strsplit(eQTL_GWAS$variant_id, "_", fixed = TRUE))[, 4]
rsid_map <- dplyr::select(MDD_GWAS, CHROM, POS, ID)
eQTL_GWAS <- left_join(eQTL_GWAS, rsid_map, by = c("chr_hg19" = "CHROM", "pos_hg19" = "POS"))
setDT(eQTL_GWAS)
m_eqtl <- eQTL_GWAS[match(snps, eQTL_GWAS$ID),  .(SNP = ID, P = pval_nominal)]
m_disorder <- MCP_GWAS[match(snps, MCP_GWAS$SNP),  .(SNP, P, BP)]

df <- eQTL_coloc_results$Brain_Caudate_basal_ganglia_PPP6C_MCP$coloc$full$results
lead_SNPs_disorder <- df$snp[which.max(df$SNP.PP.H4)]
lead_SNPs_eQTL <- df$snp[which.max(df$SNP.PP.H4)]

r2_vec <- as.numeric(R[which(snps == lead_SNPs_disorder), ]^2)
dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$BP, `–log10(P)` = -log10(m_disorder$P), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$BP),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_disorder, LD = "r2"
)

p1 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("MCP (", lead_SNPs_disorder[1], ")"))

p_genes <- locuszoomr::gg_genetracks(loc)

dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$BP, `–log10(P)` = -log10(m_eqtl$P), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$BP),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_eQTL, LD = "r2"
)

p2 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("PPP6C eQTL Brain_Caudate_basal_ganglia (", lead_SNPs_disorder[1], ")"))

p <- p1/p2/p_genes

ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/plots/MCP_PPP6C_eQTL_Brain_Caudate_basal_ganglia.jpg", p, width = 8, height = 12, dpi = 600)


## PPP6C Brain_Cortex ----
### MDD ----
R    <- as.matrix(fread(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MDD.ld.gz"), header = FALSE))
snps <- readLines(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MDD.snplist"))
rownames(R) <- colnames(R) <- snps
R <- (R + t(R))/2
diag(R) <- 1

snps <- eQTL_coloc_results$Brain_Cortex_PPP6C_MDD$meta$snps
R <- R[snps, snps]

eQTL_GWAS <- eqtl_list[["ENSG00000119414"]]
eQTL_GWAS <- subset(eQTL_GWAS, tissue == "Brain_Cortex")
eQTL_GWAS$chr_hg19 <- as.numeric(gsub("chr", "", eQTL_GWAS$chr_hg19))
eQTL_GWAS$A1 <- do.call(rbind, strsplit(eQTL_GWAS$variant_id, "_", fixed = TRUE))[, 4]
rsid_map <- dplyr::select(MDD_GWAS, CHROM, POS, ID)
eQTL_GWAS <- left_join(eQTL_GWAS, rsid_map, by = c("chr_hg19" = "CHROM", "pos_hg19" = "POS"))
setDT(eQTL_GWAS)
m_eqtl <- eQTL_GWAS[match(snps, eQTL_GWAS$ID),  .(SNP = ID, A1, BETA = slope, SE = slope_se, EAF = maf, P = pval_nominal)]
m_disorder <- MDD_GWAS[match(snps, MDD_GWAS$ID),  .(ID, EA, BETA, SE, PVAL, POS)]

df <- eQTL_coloc_results$Brain_Cortex_PPP6C_MDD$coloc$full$results
lead_SNPs_disorder <- df$snp[which.max(df$SNP.PP.H4)]
lead_SNPs_eQTL <- df$snp[which.max(df$SNP.PP.H4)]

r2_vec <- as.numeric(R[which(snps == lead_SNPs_disorder), ]^2)
dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$POS, `–log10(P)` = -log10(m_disorder$PVAL), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$POS),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_disorder, LD = "r2"
)

p1 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("MDD (", lead_SNPs_disorder[1], ")"))

p_genes <- locuszoomr::gg_genetracks(loc)

dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$POS, `–log10(P)` = -log10(m_eqtl$P), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$POS),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_eQTL, LD = "r2"
)

p2 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("PPP6C eQTL Brain_Cortex (", lead_SNPs_eQTL[1], ")"))

p <- p1/p2/p_genes

ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/plots/MDD_PPP6C_eQTL_Brain_Cortex.jpg", p, width = 8, height = 12, dpi = 600)

### MCP ----
R    <- as.matrix(fread(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MCP.ld.gz"), header = FALSE))
snps <- readLines(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MCP.snplist"))
rownames(R) <- colnames(R) <- snps
R <- (R + t(R))/2
diag(R) <- 1

snps <- eQTL_coloc_results$Brain_Cortex_PPP6C_MCP$meta$snps
R <- R[snps, snps]

eQTL_GWAS <- eqtl_list[["ENSG00000119414"]]
eQTL_GWAS <- subset(eQTL_GWAS, tissue == "Brain_Cortex")
eQTL_GWAS$chr_hg19 <- as.numeric(gsub("chr", "", eQTL_GWAS$chr_hg19))
eQTL_GWAS$A1 <- do.call(rbind, strsplit(eQTL_GWAS$variant_id, "_", fixed = TRUE))[, 4]
rsid_map <- dplyr::select(MDD_GWAS, CHROM, POS, ID)
eQTL_GWAS <- left_join(eQTL_GWAS, rsid_map, by = c("chr_hg19" = "CHROM", "pos_hg19" = "POS"))
setDT(eQTL_GWAS)
m_eqtl <- eQTL_GWAS[match(snps, eQTL_GWAS$ID),  .(SNP = ID, P = pval_nominal)]
m_disorder <- MCP_GWAS[match(snps, MCP_GWAS$SNP),  .(SNP, P, BP)]

df <- eQTL_coloc_results$Brain_Cortex_PPP6C_MCP$coloc$full$results
lead_SNPs_disorder <- df$snp[which.max(df$SNP.PP.H4)]
lead_SNPs_eQTL <- df$snp[which.max(df$SNP.PP.H4)]

r2_vec <- as.numeric(R[which(snps == lead_SNPs_disorder), ]^2)
dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$BP, `–log10(P)` = -log10(m_disorder$P), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$BP),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_disorder, LD = "r2"
)

p1 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("MCP (", lead_SNPs_disorder[1], ")"))

p_genes <- locuszoomr::gg_genetracks(loc)

dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$BP, `–log10(P)` = -log10(m_eqtl$P), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$BP),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_eQTL, LD = "r2"
)

p2 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("PPP6C eQTL Brain_Cortex (", lead_SNPs_disorder[1], ")"))

p <- p1/p2/p_genes

ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/plots/MCP_PPP6C_eQTL_Brain_Cortex.jpg", p, width = 8, height = 12, dpi = 600)

## PPP6C Brain_Frontal_Cortex_BA9 ----
### MDD ----
R    <- as.matrix(fread(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MDD.ld.gz"), header = FALSE))
snps <- readLines(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MDD.snplist"))
rownames(R) <- colnames(R) <- snps
R <- (R + t(R))/2
diag(R) <- 1

snps <- eQTL_coloc_results$Brain_Frontal_Cortex_BA9_PPP6C_MDD$meta$snps
R <- R[snps, snps]

eQTL_GWAS <- eqtl_list[["ENSG00000119414"]]
eQTL_GWAS <- subset(eQTL_GWAS, tissue == "Brain_Frontal_Cortex_BA9")
eQTL_GWAS$chr_hg19 <- as.numeric(gsub("chr", "", eQTL_GWAS$chr_hg19))
eQTL_GWAS$A1 <- do.call(rbind, strsplit(eQTL_GWAS$variant_id, "_", fixed = TRUE))[, 4]
rsid_map <- dplyr::select(MDD_GWAS, CHROM, POS, ID)
eQTL_GWAS <- left_join(eQTL_GWAS, rsid_map, by = c("chr_hg19" = "CHROM", "pos_hg19" = "POS"))
setDT(eQTL_GWAS)
m_eqtl <- eQTL_GWAS[match(snps, eQTL_GWAS$ID),  .(SNP = ID, A1, BETA = slope, SE = slope_se, EAF = maf, P = pval_nominal)]
m_disorder <- MDD_GWAS[match(snps, MDD_GWAS$ID),  .(ID, EA, BETA, SE, PVAL, POS)]

df <- eQTL_coloc_results$Brain_Frontal_Cortex_BA9_PPP6C_MDD$coloc$full$results
lead_SNPs_disorder <- df$snp[which.max(df$SNP.PP.H4)]
lead_SNPs_eQTL <- df$snp[which.max(df$SNP.PP.H4)]

r2_vec <- as.numeric(R[which(snps == lead_SNPs_disorder), ]^2)
dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$POS, `–log10(P)` = -log10(m_disorder$PVAL), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$POS),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_disorder, LD = "r2"
)

p1 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("MDD (", lead_SNPs_disorder[1], ")"))

p_genes <- locuszoomr::gg_genetracks(loc)

dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$POS, `–log10(P)` = -log10(m_eqtl$P), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$POS),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_eQTL, LD = "r2"
)

p2 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("PPP6C eQTL Brain_Frontal_Cortex_BA9 (", lead_SNPs_eQTL[1], ")"))

p <- p1/p2/p_genes

ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/plots/MDD_PPP6C_eQTL_Brain_Frontal_Cortex_BA9.jpg", p, width = 8, height = 12, dpi = 600)

### MCP ----
R    <- as.matrix(fread(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MCP.ld.gz"), header = FALSE))
snps <- readLines(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MCP.snplist"))
rownames(R) <- colnames(R) <- snps
R <- (R + t(R))/2
diag(R) <- 1

snps <- eQTL_coloc_results$Brain_Frontal_Cortex_BA9_PPP6C_MCP$meta$snps
R <- R[snps, snps]

eQTL_GWAS <- eqtl_list[["ENSG00000119414"]]
eQTL_GWAS <- subset(eQTL_GWAS, tissue == "Brain_Frontal_Cortex_BA9")
eQTL_GWAS$chr_hg19 <- as.numeric(gsub("chr", "", eQTL_GWAS$chr_hg19))
eQTL_GWAS$A1 <- do.call(rbind, strsplit(eQTL_GWAS$variant_id, "_", fixed = TRUE))[, 4]
rsid_map <- dplyr::select(MDD_GWAS, CHROM, POS, ID)
eQTL_GWAS <- left_join(eQTL_GWAS, rsid_map, by = c("chr_hg19" = "CHROM", "pos_hg19" = "POS"))
setDT(eQTL_GWAS)
m_eqtl <- eQTL_GWAS[match(snps, eQTL_GWAS$ID),  .(SNP = ID, P = pval_nominal)]
m_disorder <- MCP_GWAS[match(snps, MCP_GWAS$SNP),  .(SNP, P, BP)]

df <- eQTL_coloc_results$Brain_Frontal_Cortex_BA9_PPP6C_MCP$coloc$full$results
lead_SNPs_disorder <- df$snp[which.max(df$SNP.PP.H4)]
lead_SNPs_eQTL <- df$snp[which.max(df$SNP.PP.H4)]

r2_vec <- as.numeric(R[which(snps == lead_SNPs_disorder), ]^2)
dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$BP, `–log10(P)` = -log10(m_disorder$P), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$BP),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_disorder, LD = "r2"
)

p1 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("MCP (", lead_SNPs_disorder[1], ")"))

p_genes <- locuszoomr::gg_genetracks(loc)

dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$BP, `–log10(P)` = -log10(m_eqtl$P), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$BP),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_eQTL, LD = "r2"
)

p2 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("PPP6C eQTL Brain_Frontal_Cortex_BA9 (", lead_SNPs_disorder[1], ")"))

p <- p1/p2/p_genes

ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/plots/MCP_PPP6C_eQTL_Brain_Frontal_Cortex_BA9.jpg", p, width = 8, height = 12, dpi = 600)

## PPP6C Brain_Putamen_basal_ganglia ----
### MDD ----
R    <- as.matrix(fread(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MDD.ld.gz"), header = FALSE))
snps <- readLines(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MDD.snplist"))
rownames(R) <- colnames(R) <- snps
R <- (R + t(R))/2
diag(R) <- 1

snps <- eQTL_coloc_results$Brain_Putamen_basal_ganglia_PPP6C_MDD$meta$snps
R <- R[snps, snps]

eQTL_GWAS <- eqtl_list[["ENSG00000119414"]]
eQTL_GWAS <- subset(eQTL_GWAS, tissue == "Brain_Putamen_basal_ganglia")
eQTL_GWAS$chr_hg19 <- as.numeric(gsub("chr", "", eQTL_GWAS$chr_hg19))
eQTL_GWAS$A1 <- do.call(rbind, strsplit(eQTL_GWAS$variant_id, "_", fixed = TRUE))[, 4]
rsid_map <- dplyr::select(MDD_GWAS, CHROM, POS, ID)
eQTL_GWAS <- left_join(eQTL_GWAS, rsid_map, by = c("chr_hg19" = "CHROM", "pos_hg19" = "POS"))
setDT(eQTL_GWAS)
m_eqtl <- eQTL_GWAS[match(snps, eQTL_GWAS$ID),  .(SNP = ID, A1, BETA = slope, SE = slope_se, EAF = maf, P = pval_nominal)]
m_disorder <- MDD_GWAS[match(snps, MDD_GWAS$ID),  .(ID, EA, BETA, SE, PVAL, POS)]

df <- eQTL_coloc_results$Brain_Putamen_basal_ganglia_PPP6C_MDD$coloc$full$results
lead_SNPs_disorder <- df$snp[which.max(df$SNP.PP.H4)]
lead_SNPs_eQTL <- df$snp[which.max(df$SNP.PP.H4)]

r2_vec <- as.numeric(R[which(snps == lead_SNPs_disorder), ]^2)
dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$POS, `–log10(P)` = -log10(m_disorder$PVAL), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$POS),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_disorder, LD = "r2"
)

p1 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("MDD (", lead_SNPs_disorder[1], ")"))

p_genes <- locuszoomr::gg_genetracks(loc)

dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$POS, `–log10(P)` = -log10(m_eqtl$P), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$POS),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_eQTL, LD = "r2"
)

p2 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("PPP6C eQTL Brain_Putamen_basal_ganglia (", lead_SNPs_eQTL[1], ")"))

p <- p1/p2/p_genes

ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/plots/MDD_PPP6C_eQTL_Brain_Putamen_basal_ganglia.jpg", p, width = 8, height = 12, dpi = 600)

### MCP ----
R    <- as.matrix(fread(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MCP.ld.gz"), header = FALSE))
snps <- readLines(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/PPP6C_MCP.snplist"))
rownames(R) <- colnames(R) <- snps
R <- (R + t(R))/2
diag(R) <- 1

snps <- eQTL_coloc_results$Brain_Putamen_basal_ganglia_PPP6C_MCP$meta$snps
R <- R[snps, snps]

eQTL_GWAS <- eqtl_list[["ENSG00000119414"]]
eQTL_GWAS <- subset(eQTL_GWAS, tissue == "Brain_Putamen_basal_ganglia")
eQTL_GWAS$chr_hg19 <- as.numeric(gsub("chr", "", eQTL_GWAS$chr_hg19))
eQTL_GWAS$A1 <- do.call(rbind, strsplit(eQTL_GWAS$variant_id, "_", fixed = TRUE))[, 4]
rsid_map <- dplyr::select(MDD_GWAS, CHROM, POS, ID)
eQTL_GWAS <- left_join(eQTL_GWAS, rsid_map, by = c("chr_hg19" = "CHROM", "pos_hg19" = "POS"))
setDT(eQTL_GWAS)
m_eqtl <- eQTL_GWAS[match(snps, eQTL_GWAS$ID),  .(SNP = ID, P = pval_nominal)]
m_disorder <- MCP_GWAS[match(snps, MCP_GWAS$SNP),  .(SNP, P, BP)]

df <- eQTL_coloc_results$Brain_Putamen_basal_ganglia_PPP6C_MCP$coloc$full$results
lead_SNPs_disorder <- df$snp[which.max(df$SNP.PP.H4)]
lead_SNPs_eQTL <- df$snp[which.max(df$SNP.PP.H4)]

r2_vec <- as.numeric(R[which(snps == lead_SNPs_disorder), ]^2)
dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$BP, `–log10(P)` = -log10(m_disorder$P), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$BP),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_disorder, LD = "r2"
)

p1 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("MCP (", lead_SNPs_disorder[1], ")"))

p_genes <- locuszoomr::gg_genetracks(loc)

dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$BP, `–log10(P)` = -log10(m_eqtl$P), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$BP),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_eQTL, LD = "r2"
)

p2 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("PPP6C eQTL Brain_Putamen_basal_ganglia (", lead_SNPs_disorder[1], ")"))

p <- p1/p2/p_genes

ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/plots/MCP_PPP6C_eQTL_Brain_Putamen_basal_ganglia.jpg", p, width = 8, height = 12, dpi = 600)


## SP4 Brain_Nucleus_accumbens_basal_ganglia ----
### MDD ----
R    <- as.matrix(fread(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/SP4_MDD.ld.gz"), header = FALSE))
snps <- readLines(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/eqtl/SP4_MDD.snplist"))
rownames(R) <- colnames(R) <- snps
R <- (R + t(R))/2
diag(R) <- 1

snps <- eQTL_coloc_results$Brain_Nucleus_accumbens_basal_ganglia_SP4_MDD$meta$snps
R <- R[snps, snps]

eQTL_GWAS <- eqtl_list[["ENSG00000105866"]]
eQTL_GWAS <- subset(eQTL_GWAS, tissue == "Brain_Nucleus_accumbens_basal_ganglia")
eQTL_GWAS$chr_hg19 <- as.numeric(gsub("chr", "", eQTL_GWAS$chr_hg19))
eQTL_GWAS$A1 <- do.call(rbind, strsplit(eQTL_GWAS$variant_id, "_", fixed = TRUE))[, 4]
rsid_map <- dplyr::select(MDD_GWAS, CHROM, POS, ID)
eQTL_GWAS <- left_join(eQTL_GWAS, rsid_map, by = c("chr_hg19" = "CHROM", "pos_hg19" = "POS"))
setDT(eQTL_GWAS)
m_eqtl <- eQTL_GWAS[match(snps, eQTL_GWAS$ID),  .(SNP = ID, A1, BETA = slope, SE = slope_se, EAF = maf, P = pval_nominal)]
m_disorder <- MDD_GWAS[match(snps, MDD_GWAS$ID),  .(ID, EA, BETA, SE, PVAL, POS)]

df <- eQTL_coloc_results$Brain_Nucleus_accumbens_basal_ganglia_SP4_MDD$coloc$full$results
lead_SNPs_disorder <- df$snp[which.max(df$SNP.PP.H4)]
lead_SNPs_eQTL <- df$snp[which.max(df$SNP.PP.H4)]

r2_vec <- as.numeric(R[which(snps == lead_SNPs_disorder), ]^2)
dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$POS, `–log10(P)` = -log10(m_disorder$PVAL), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$POS),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_disorder, LD = "r2"
)

p1 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("MDD (", lead_SNPs_disorder[1], ")"))

p_genes <- locuszoomr::gg_genetracks(loc)

dt <- data.table(CHR = eQTL_GWAS$chr_hg19[1], BP = m_disorder$POS, `–log10(P)` = -log10(m_eqtl$P), SNP = snps, r2 = r2_vec)
loc <- locuszoomr::locus(
  data = dt, seqname = dt$CHR[1], xrange = range(m_disorder$POS),
  ens_db = EnsDb.Hsapiens.v75, chrom = "CHR", pos = "BP", yvar = "–log10(P)", labs = "SNP",
  index_snp = lead_SNPs_eQTL, LD = "r2"
)

p2 <- locuszoomr::gg_scatter(loc) + 
  ggtitle(paste0("SP4 eQTL Brain_Nucleus_accumbens_basal_ganglia (", lead_SNPs_eQTL[1], ")"))

p <- p1/p2/p_genes

ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/plots/MDD_SP4_eQTL_Brain_Nucleus_accumbens_basal_ganglia.jpg", p, width = 8, height = 12, dpi = 600)
