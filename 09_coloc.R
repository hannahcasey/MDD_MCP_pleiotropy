library(data.table)
library(susieR)
library(coloc)
library(LDlinkR)
library(dplyr)

# Inputs ----
LDLINK_TOKEN <- "448bee827178"
ld_pop <- "EUR"                        
build <- "grch37"
kb_window <- 250

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

## Lead SNPs ----
FUMA_leads <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.lead.csv", header = T)

## 1000g BIM file ----
BIM <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/1000g/EUR.bim")
setDT(BIM)

# LD matrices -----
LD_matrices <- list()

PLINK <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/plink_mac_20250615/plink"
REF_PREFIX <- "/Users/hannahcasey/Desktop/PhD/projects/UKB_chronic_pain_depressive_symptoms/resources/1000g/EUR"
LD_DIR <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/coloc"

for(lead_idx in seq_len(nrow(FUMA_leads))){
  
  SNP_chr  <- FUMA_leads$CHR[lead_idx]
  SNP_pos  <- FUMA_leads$LEAD_BP[lead_idx]
  SNP_rsid <- FUMA_leads$LEAD_SNP[lead_idx]
  start_pos <- SNP_pos - 250000L
  end_pos   <- SNP_pos + 250000L
  
  region_df <- BIM %>%
    filter(V1 == SNP_chr, between(V4, start_pos, end_pos)) %>%
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

for(lead_idx in seq_len(nrow(FUMA_leads))){
  
  SNP_rsid <- FUMA_leads$LEAD_SNP[lead_idx]
  R    <- as.matrix(fread(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/coloc/", SNP_rsid, ".ld.gz"), header = FALSE))
  snps <- readLines(paste0("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/LD/coloc/", SNP_rsid, ".snplist"))
  rownames(R) <- colnames(R) <- snps
  # R <- (R + t(R))/2
  # diag(R) <- 1

  BIM_sub <- BIM[match(snps, BIM$V2), .(SNP = V2, REF = V6, ALT = V5)]

  m_mcp <- MCP_GWAS[match(snps, MCP_GWAS$SNP), .(SNP, A1 = ALLELE1, BETA, SE, EAF = A1FREQ, POS = BP, P, INFO)]
  m_mdd  <- MDD_GWAS[match(snps, MDD_GWAS$ID),  .(SNP = ID, A1 = EA, BETA, SE, EAF, P = PVAL, INFO = IMPINFO)]
  
  
  ### flip allele ----
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
  
  R         <- R[idx, idx, drop=FALSE]
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

coloc_pub <- data.frame(
  FUMA_lead      = character(),
  coloc_method        = character(),
  PP4           = numeric(),     # PP.H4.abf
  MCP_CS_size   = integer(),
  MDD_CS_size   = integer(),
  MCP_CS_PIPsum = numeric(),
  MDD_CS_PIPsum = numeric(),
  top_MCP_snp   = character(),
  top_MDD_snp   = character(),
  SNP_PP4_top   = numeric(),   
  MCP_CS_snps   = character(),
  MDD_CS_snps   = character(),
  stringsAsFactors = FALSE
)

for(k in names(coloc_results)){
  x <- coloc_results[[k]]
  m <- x$coloc$method
  
  if (m == "susie") {
    ss <- x$coloc$summary
    ss <- ss[ss$PP.H4.abf >= 0.8, , drop = FALSE]
    if (nrow(ss) == 0) next
    
    for (i in seq_len(nrow(ss))) {
      cs1 <- ss$idx1[i]; cs2 <- ss$idx2[i]
      
      cs1_idx <- x$fits$MCP$sets$cs[[cs1]]
      cs2_idx <- x$fits$MDD$sets$cs[[cs2]]
      l1 <- x$fits$MCP$sets$cs_index[cs1]
      l2 <- x$fits$MDD$sets$cs_index[cs2]
      
      # rsID lists per CS
      cs1_rs <- x$meta$snps[cs1_idx]
      cs2_rs <- x$meta$snps[cs2_idx]
      
      # sizes and PIP sums
      cs1_size <- length(cs1_idx)
      cs2_size <- length(cs2_idx)
      cs1_pipsum <- sum(x$fits$MCP$alpha[l1, cs1_idx], na.rm = TRUE)
      cs2_pipsum <- sum(x$fits$MDD$alpha[l2, cs2_idx], na.rm = TRUE)
      
      # top SNP (max PIP) and its PIP
      top1 <- cs1_idx[ which.max(x$fits$MCP$pip[cs1_idx]) ]
      top2 <- cs2_idx[ which.max(x$fits$MDD$pip[cs2_idx]) ]
      
      coloc_pub <- rbind(coloc_pub, data.frame(
        FUMA_lead      = k,
        coloc_method        = m,
        PP4           = ss$PP.H4.abf[i],
        MCP_CS_size   = cs1_size,
        MDD_CS_size   = cs2_size,
        MCP_CS_PIPsum = cs1_pipsum,
        MDD_CS_PIPsum = cs2_pipsum,
        top_MCP_snp   = x$meta$snps[top1],
        top_MDD_snp   = x$meta$snps[top2],
        SNP_PP4_top   = NA, 
        MCP_CS_snps   = paste(cs1_rs, collapse=";"),
        MDD_CS_snps   = paste(cs2_rs, collapse=";"),
        stringsAsFactors = FALSE
      ))
    }
    
  } else if (m == "abf") {
    if (x$coloc$summary$PP.H4.abf >= 0.8) {
      # top shared SNP under H4 = max SNP.PP.H4
      ord <- order(-x$coloc$full$results$SNP.PP.H4)
      top_shared <- x$coloc$full$results$snp[ord][1]
      top_shared_pp4 <- x$coloc$full$results$SNP.PP.H4[ord][1]
      
      coloc_pub <- rbind(coloc_pub, data.frame(
        FUMA_lead      = k,
        coloc_method        = m,
        PP4           = x$coloc$summary$PP.H4.abf,
        MCP_CS_size   = NA_integer_,
        MDD_CS_size   = NA_integer_,
        MCP_CS_PIPsum = NA_real_,
        MDD_CS_PIPsum = NA_real_,
        top_MCP_snp   = top_shared,
        top_MDD_snp   = top_shared,
        SNP_PP4_top   = top_shared_pp4,
        MCP_CS_snps   = NA,
        MDD_CS_snps   = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
}

coloc_pub <- coloc_pub[!duplicated(coloc_pub[, c("top_MCP_snp","top_MDD_snp")]), ]
write.csv(coloc_pub, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/cs.csv", row.names = F)

