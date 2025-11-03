library(data.table)
library(susieR)
library(coloc)
library(dplyr)

# Load in data ----
## MDD ----
file <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MDD/pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv"
MDD_GWAS <- fread(file, sep = "\t")
MDD_GWAS$EAF <- ((MDD_GWAS$FCAS * MDD_GWAS$NCAS) + (MDD_GWAS$FCON * MDD_GWAS$NCON)) / (MDD_GWAS$NCAS + MDD_GWAS$NCON)
MDD_n <- 1639572
MDD_s <- 357636 /	(357636 + 1281936)

## MCP ----
MCP_GWAS <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MCP/ChronicPain.stats", header = T)
setDT(MCP_GWAS)
MCP_n <- 387649

## Share Loci ----
FUMA_loci <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.loci.csv", header = T)

## 1000g BIM file ----
BIM <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/1000g/EUR.bim")
setDT(BIM)

# LD matrices -----
LD_matrices <- list()

PLINK <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/plink_mac_20250615/plink"
REF_PREFIX <- "/Users/hannahcasey/Desktop/PhD/projects/UKB_chronic_pain_depressive_symptoms/resources/1000g/EUR"
LD_DIR <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/coloc"

for(loci_idx in seq_len(nrow(FUMA_loci))){
  
  SNP_rsid <- FUMA_loci$LEAD_SNP[loci_idx]
  CHR <- FUMA_loci$CHR[loci_idx]
  start_pos  <- FUMA_loci$MinBP[loci_idx]
  end_pos  <- FUMA_loci$MaxBP[loci_idx]

  region_df <- BIM %>%
    filter(V1 == CHR, between(V4, start_pos, end_pos)) %>%
    distinct(V2, V4) %>%
    rename(rsid = V2, pos = V4)
  
  common_SNPs <- intersect(intersect(region_df$rsid, MDD_GWAS$ID), MCP_GWAS$SNP)
  
  snplist_file <- file.path(LD_DIR, paste0("snplist_", SNP_rsid, ".txt"))
  write.table(common_SNPs, snplist_file, quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  out_prefix <- file.path(LD_DIR, SNP_rsid)
  
  cmd <- paste(
    shQuote(PLINK),
    "--bfile", shQuote(REF_PREFIX),
    "--extract", shQuote(snplist_file),
    "--r square gz",
    "--write-snplist",
    "--a1-allele /Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/1000g/ref_alleles.txt 2 1",
    "--out", shQuote(out_prefix)
  )
  
  system(cmd)
}

# SuSie + coloc ----
coloc_results <- list()

for(loci_idx in seq_len(nrow(FUMA_loci))){
  
  SNP_rsid <- FUMA_loci$LEAD_SNP[loci_idx]
  R    <- as.matrix(fread(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/coloc/", SNP_rsid, ".ld.gz"), header = FALSE))
  snps <- readLines(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/coloc/", SNP_rsid, ".snplist"))
  if(length(snps) < 2) next
  rownames(R) <- colnames(R) <- snps

  BIM_sub <- BIM[match(snps, BIM$V2), .(SNP = V2, REF = V6, ALT = V5)]

  m_mcp <- MCP_GWAS[match(snps, MCP_GWAS$SNP), .(SNP, A1 = ALLELE1, BETA, SE, EAF = A1FREQ, POS = BP, P, INFO)]
  m_mdd  <- MDD_GWAS[match(snps, MDD_GWAS$ID),  .(SNP = ID, A1 = EA, BETA, SE, EAF, P = PVAL, INFO = IMPINFO)]
  
  
  ## flip allele ----
  flip_mcp <- fifelse(m_mcp$A1 == BIM_sub$ALT, -1,
                      fifelse(m_mcp$A1 == BIM_sub$REF,  1, NA_real_))
  flip_mdd <- fifelse(m_mdd$A1 == BIM_sub$ALT, -1,
                      fifelse(m_mdd$A1 == BIM_sub$REF,  1, NA_real_))
  
  ## palindromic ----
  is_pal <- function(a,b)
    ( (a %chin% c("A","T") & b %chin% c("A","T")) |
        (a %chin% c("C","G") & b %chin% c("C","G")) )
  
  pal <- is_pal(BIM_sub$ALT, BIM_sub$REF)
  maf  <- pmin(m_mcp$EAF, 1 - m_mcp$EAF)
  
  keep <- !(pal & maf > 0.40) & !is.na(flip_mcp) & !is.na(flip_mdd)
  idx <- which(keep)
  
  R <- R[idx, idx, drop=FALSE]
  R <- (R + t(R)) / 2
  diag(R) <- 1
  
  snps      <- snps[idx]
  BIM_sub  <- BIM_sub[idx]
  m_mcp    <- m_mcp[idx]
  m_mdd     <- m_mdd[idx]
  flip_mcp <- flip_mcp[idx]
  flip_mdd <- flip_mdd[idx]
  
  m_mcp$BETA <- flip_mcp * m_mcp$BETA
  m_mdd$BETA  <- flip_mdd  * m_mdd$BETA
  
  d_mcp <- list(
    beta    = m_mcp$BETA,
    varbeta = (m_mcp$SE)^2,
    snp     = snps,
    position = m_mcp$POS, 
    type    = "quant",
    N       = MCP_n,
    MAF     = pmin(m_mcp$EAF, 1 - m_mcp$EAF), 
    LD      = R
  )
  
  d_mdd <- list(
    beta    = m_mdd$BETA,
    varbeta = (m_mdd$SE)^2,
    snp     = snps,
    position = m_mcp$POS, 
    type    = "cc",
    N       = MDD_n,
    MAF     = pmin(m_mdd$EAF, 1 - m_mdd$EAF), 
    LD      = R,
    s = MDD_s
  )
  
  # z <- m_mcp$BETA/m_mcp$SE
  # condz = kriging_rss(z, R, n=MCP_n)
  # condz$plot
  # 
  # z <- m_mdd$BETA/m_mdd$SE
  # condz = kriging_rss(z, R, n=MDD_n)
  # condz$plot
  
  ### Run SuSiE ----
  fit_mcp <- runsusie(d_mcp, maxit = 1000, L = 10)
  fit_mdd  <- runsusie(d_mdd, maxit = 1000, L = 10)
  
  ### Run coloc ----
  if(is.null(fit_mcp$sets$cs_index) | is.null(fit_mdd$sets$cs_index)){
    
    method_used <- "abf"
    d1_abf <- d_mcp[c("beta","varbeta","snp","type","N","MAF")]
    d2_abf <- d_mdd[c("beta","varbeta","snp","type","N","MAF","s")]
    cres <- coloc.abf(d1_abf, d2_abf)
    cres$summary <- as.data.frame(t(cres$summary))
    
  } else{
    
    method_used <- "susie"
    cres <- coloc.susie(fit_mcp, fit_mdd)
  }
  
    coloc_results[[SNP_rsid]] <- list(
      fits = list(MCP = fit_mcp, MDD = fit_mdd),
      coloc = list(full = cres, summary = cres$summary, method = method_used),
      meta  = list(snps = snps, n_vars = length(snps))
    )
}

# Save outputs ----
saveRDS(coloc_results, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/cross_trait.rds")

# Format results ----
pph4_from_summary <- function(s) {
  if (is.null(s)) return(NA_real_)
  if (is.data.frame(s) && "PP.H4.abf" %in% names(s)) return(as.numeric(s$PP.H4.abf[1]))
  if (is.atomic(s) && !is.list(s) && !is.null(names(s)) && "PP.H4.abf" %in% names(s)) return(as.numeric(s[["PP.H4.abf"]]))
  if (is.data.frame(s) && "PP.H4.abf" %in% rownames(s)) return(as.numeric(s["PP.H4.abf", 1]))
  NA_real_
}

out <- vector("list", 0L)

for (k in names(coloc_results)) {
  x <- coloc_results[[k]]
  m <- x$coloc$method

  fit_MCP  <- x$fits$MCP
  fit_MDD <- x$fits$MDD
  
  if (m == "susie") {
    ss <- x$coloc$summary
    
    for (i in seq_len(nrow(ss))) {
      cs1 <- ss$idx1[i]; cs2 <- ss$idx2[i]
      cs1_idx <- fit_MCP$sets$cs[[paste0("L", cs1)]]
      cs2_idx <- fit_MDD$sets$cs[[paste0("L", cs2)]]
      
      top1 <- cs1_idx[ which.max(fit_MCP$pip[cs1_idx]) ]
      top2 <- cs2_idx[ which.max(fit_MDD$pip[cs2_idx]) ]
      
      out[[length(out)+1L]] <- data.frame(
        Lead_SNP = k,
        coloc_method       = m,
        PP4                = ss$PP.H4.abf[i],
        MCP_CS_size   = length(cs1_idx),
        MDD_CS_size       = length(cs2_idx),
        top_MCP_snp   = x$meta$snps[top1],
        top_MDD_snp       = x$meta$snps[top2],
        MCP_CS_snps   = paste(x$meta$snps[cs1_idx], collapse=";"),
        MDD_CS_snps       = paste(x$meta$snps[cs2_idx], collapse=";"),
        stringsAsFactors   = FALSE
      )
    }
    
  } else if (m == "abf") {
    
    pp4 <- pph4_from_summary(x$coloc$summary)
    res <- x$coloc$full$results
    ord <- order(-res$SNP.PP.H4)
    top_shared     <- res$snp[ord][1]
    
    out[[length(out)+1L]] <- data.frame(
      Lead_SNP = k,
      coloc_method = m,
      PP4 = pp4,
      MDD_CS_size = NA_integer_,
      MCP_CS_size = NA_integer_,
      top_MCP_snp = top_shared,
      top_MDD_snp = top_shared,
      MCP_CS_snps = NA_character_,
      MDD_CS_snps = NA_character_,
      stringsAsFactors = FALSE
    )
  }
}

out_long <- if (length(out)) do.call(rbind, out) else
  data.frame(
    Lead_SNP = character(),
    coloc_method = character(),
    PP4 = numeric(),
    MCP_CS_size = integer(),
    MDD_CS_size = integer(),
    top_MCP_snp = character(),
    top_MDD_snp = character(),
    MCP_CS_snps = character(),
    MDD_CS_snps = character(),
    stringsAsFactors = FALSE
  )


FUMA_loci <- read.csv("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.loci.csv", sep = "\t")

out_long <- out_long %>%
  left_join(FUMA_loci, by = c("Lead_SNP" = "LEAD_SNP")) %>%
  dplyr::select(locusnum, CHR, MinBP, MaxBP, Lead_SNP, LEAD_BP, 
         coloc_method, 
         MDD_CS_size, MDD_CS_snps, top_MDD_snp,
         MCP_CS_size, MCP_CS_snps, top_MCP_snp,
         PP4) %>%
  dplyr::rename(`Locus number conjFDR` = locusnum,
                Chromosome = CHR,
                `Min pos` = MinBP,
                `Max pos` = MaxBP,
                `Lead SNP` = Lead_SNP,
                `Lead SNP pos` = LEAD_BP,
                `Coloc approach` = coloc_method,
                `MDD credible set (n)` = MDD_CS_size,
                `MCP credible set (n)` = MCP_CS_size,
                `Top MDD SNP` = top_MDD_snp,
                `Top MCP SNP` = top_MCP_snp,
                `MDD credible set SNPs` = MDD_CS_snps,
                `MCP credible set SNPs` = MCP_CS_snps
  )



out_short <- subset(out_long, PP4 > 0.80)
out_short <- out_short %>%
  distinct(across(all_of(c("Top MCP SNP", "Top MDD SNP", "MDD credible set SNPs", "MCP credible set SNPs"))), .keep_all = TRUE)

write.csv(out_long, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/cs_cross_trait.csv", row.names = F)
write.csv(out_short, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/cs_cross_trait_significant.csv", row.names = F)
