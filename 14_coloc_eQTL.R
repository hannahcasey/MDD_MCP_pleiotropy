library(data.table)
library(susieR)
library(coloc)
library(dplyr)
library(biomaRt)
library(stringr)
library(purrr)
library(tidyr)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(rrapply)
library(tibble)

# Load in data ----
## MDD ----
file <- "./resources/summary_statistics/MDD/pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv"
MDD_GWAS <- fread(file, sep = "\t")
MDD_GWAS$EAF <- ((MDD_GWAS$FCAS * MDD_GWAS$NCAS) + (MDD_GWAS$FCON * MDD_GWAS$NCON)) / (MDD_GWAS$NCAS + MDD_GWAS$NCON)
MDD_n <- 1639572
MDD_s <- 357636 /	(357636 + 1281936)

## MCP ----
MCP_GWAS <- read.table("./resources/summary_statistics/MCP/ChronicPain.stats", header = T)
setDT(MCP_GWAS)
MCP_n <- 387649

## 1000g BIM file ----
BIM <- read.table("./resources/1000g/EUR.bim")
setDT(BIM)

## eQTL ----
dir <- "./resources/summary_statistics/eQTL/brain/by_gene"
files <- list.files(
  dir, pattern = "*.csv", full.names = TRUE
)

gene_names <- str_replace(basename(files), "\\.v8\\.EUR\\.brain\\.allchr\\.csv$", "")

eqtl_list <- set_names(
  map(files, ~ read.csv(.x)),
  gene_names
)

eqtl_list <- map(eqtl_list, ~ filter(.x, !is.na(slope)))

ch <- import.chain("./resources/hg38ToHg19.over.chain")

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

## Genes ----
genes <- read.csv("./resources/genes.cvs", header = T, sep = " ")

genes_ensembl <- sort(unique(unlist(strsplit(genes$ensembl_gene_id, ", ", fixed = TRUE))))

mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
gene_coords <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name",
                 "start_position", "end_position", "strand"),
  filters = "ensembl_gene_id",
  values = genes_ensembl,
  mart = mart
)

gene_coords <- left_join(gene_coords, genes, by = "ensembl_gene_id")

# LD matrices -----
LD_matrices <- list()

PLINK <- "./resources/plink_mac_20250615/plink"
REF_PREFIX <- "./resources/1000g/EUR"
LD_DIR <- "./resources/LD/"

for(gene_idx in seq_len(nrow(gene_coords))){
  
  gene_eQTL <- eqtl_list[[gene_coords$ensembl_gene_id[[gene_idx]]]]
  if(nrow(gene_eQTL) == 0) next
  gene_eQTL <- dplyr::select(gene_eQTL, chr_hg19, pos_hg19)
  gene_eQTL$chr_hg19 <- as.numeric(gsub("chr", "", gene_eQTL$chr_hg19))
  
  gene_chr  <- gene_coords$chromosome_name[gene_idx]
  gene_name <- gene_coords$SYMBOL[gene_idx]
  start_pos <- gene_coords$start_position[gene_idx] - 500000L
  end_pos   <- gene_coords$end_position[gene_idx] + 500000L
  
  region_df <- BIM %>%
    dplyr::rename(chr = V1, 
           rsid = V2, 
           dummy = V3, 
           pos = V4, 
           A1 = V5, 
           A2 = V6) %>%
    filter(chr == gene_chr, between(pos, start_pos, end_pos)) %>%
    inner_join(gene_eQTL, by = c("chr" = "chr_hg19", "pos" = "pos_hg19")) %>%
    distinct(rsid, pos)
  
  MDD_common_SNPs <- intersect(region_df$rsid, MDD_GWAS$ID)
  MCP_common_SNPs <- intersect(region_df$rsid, MCP_GWAS$SNP)
  
  MDD_snplist_file <- file.path(LD_DIR, paste0("eqtl/snplist_MDD_", gene_name, ".txt"))
  MCP_snplist_file <- file.path(LD_DIR, paste0("eqtl/snplist_MCP_", gene_name, ".txt"))
  
  write.table(MDD_common_SNPs, MDD_snplist_file, quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(MCP_common_SNPs, MCP_snplist_file, quote=FALSE, row.names=FALSE, col.names=FALSE)
  
  MDD_out_prefix <- file.path(LD_DIR,  "eqtl", paste0(gene_name, "_MDD"))
  cmd <- paste(
    shQuote(PLINK),
    "--bfile", shQuote(REF_PREFIX),
    "--extract", shQuote(MDD_snplist_file),
    "--r square gz",
    "--write-snplist",
    "--a1-allele ./resources/1000g/ref_alleles.txt 2 1",
    "--out", shQuote(MDD_out_prefix)
  )
  system(cmd)
  
  MCP_out_prefix <- file.path(LD_DIR,  "eqtl", paste0(gene_name, "_MCP"))
  cmd <- paste(
    shQuote(PLINK),
    "--bfile", shQuote(REF_PREFIX),
    "--extract", shQuote(MCP_snplist_file),
    "--r square gz",
    "--write-snplist",
    "--a1-allele ./resources/1000g/ref_alleles.txt 2 1",
    "--out", shQuote(MCP_out_prefix)
  )
  system(cmd)
}

# SuSie + coloc ----
coloc_results <- list()

tissue_list <- unique(eqtl_list$ENSG00000078295$tissue)
rsid_map <- dplyr::select(MDD_GWAS, CHROM, POS, ID)

gtex_eur_brain_sample_sizes <- tribble(
  ~tissue,                                 ~N,
  "Brain_Amygdala",                        119,
  "Brain_Anterior_cingulate_cortex_BA24",  135,
  "Brain_Caudate_basal_ganglia",           172,
  "Brain_Cerebellar_Hemisphere",           157,
  "Brain_Cerebellum",                      188,
  "Brain_Cortex",                          184,
  "Brain_Frontal_Cortex_BA9",              158,
  "Brain_Hippocampus",                     150,
  "Brain_Hypothalamus",                    157,
  "Brain_Nucleus_accumbens_basal_ganglia", 181,
  "Brain_Putamen_basal_ganglia",           153,
  "Brain_Spinal_cord_cervical_c-1",        115,
  "Brain_Substantia_nigra",                101
)

tissue_n_lookup <- setNames(gtex_eur_brain_sample_sizes$N,
                            gtex_eur_brain_sample_sizes$tissue)

for(gene_idx in 44:nrow(gene_coords)){
  
  gene_name <- gene_coords$SYMBOL[gene_idx]
  ensembl_gene_id <- gene_coords$ensembl_gene_id[gene_idx]
  
  gene_eQTL <- eqtl_list[[gene_coords$ensembl_gene_id[[gene_idx]]]]
  if(nrow(gene_eQTL) == 0) next
  
  for(disorder in c("MDD", "MCP")){
    
    for(brain_tissue in tissue_list){
      
      R    <- as.matrix(fread(paste0("./resources/LD/eqtl/", gene_name, "_", disorder, ".ld.gz"), header = FALSE))
      snps <- readLines(paste0("./resources/LD/eqtl/", gene_name, "_", disorder, ".snplist"))
      rownames(R) <- colnames(R) <- snps
      BIM_sub <- BIM[match(snps, BIM$V2), .(SNP = V2, REF = V6, ALT = V5)]
      
      if(disorder == "MCP"){
        m_disorder <- MCP_GWAS[match(snps, MCP_GWAS$SNP), .(SNP, A1 = ALLELE1, BETA, SE, EAF = A1FREQ, CHR, POS = BP, P, INFO)]
      } else{
        m_disorder  <- MDD_GWAS[match(snps, MDD_GWAS$ID),  .(SNP = ID, CHR = CHROM, POS, A1 = EA, BETA, SE, EAF, P = PVAL, INFO = IMPINFO)]
      }
      
      eQTL_GWAS <- eqtl_list[[ensembl_gene_id]]
      eQTL_GWAS <- subset(eQTL_GWAS, tissue == brain_tissue)
      if (nrow(eQTL_GWAS) == 0) next ## No eQTL when gene not sufficiently expressed
      eQTL_GWAS$chr_hg19 <- as.numeric(gsub("chr", "", eQTL_GWAS$chr_hg19))
      eQTL_GWAS$A1 <- do.call(rbind, strsplit(eQTL_GWAS$variant_id, "_", fixed = TRUE))[, 4]
      eQTL_GWAS <- left_join(eQTL_GWAS, rsid_map, by = c("chr_hg19" = "CHROM", "pos_hg19" = "POS"))
      eQTL_GWAS <- subset(eQTL_GWAS, !is.infinite(slope))
      setDT(eQTL_GWAS)
      m_eqtl <- eQTL_GWAS[match(snps, eQTL_GWAS$ID),  .(SNP = ID, A1, BETA = slope, SE = slope_se, EAF = maf, P = pval_nominal)]
      
      ## flip allele ----
      flip_disorder <- fifelse(m_disorder$A1 == BIM_sub$ALT, -1,
                          fifelse(m_disorder$A1 == BIM_sub$REF,  1, NA_real_))
      flip_eqtl <- fifelse(m_eqtl$A1 == BIM_sub$ALT, -1,
                          fifelse(m_eqtl$A1 == BIM_sub$REF,  1, NA_real_))
      
      ## palindromic ----
      is_pal <- function(a,b)
        ( (a %chin% c("A","T") & b %chin% c("A","T")) |
            (a %chin% c("C","G") & b %chin% c("C","G")) )
      
      pal <- is_pal(BIM_sub$ALT, BIM_sub$REF)
      maf  <- pmin(m_disorder$EAF, 1 - m_disorder$EAF)
      
      keep <- !(pal & maf > 0.40) & !is.na(flip_disorder) & !is.na(flip_eqtl)
      idx <- which(keep)
      
      R <- R[idx, idx, drop=FALSE]
      R <- (R + t(R)) / 2
      diag(R) <- 1
      
      snps <- snps[idx]
      BIM_sub <- BIM_sub[idx]
      m_eqtl <- m_eqtl[idx]
      m_disorder <- m_disorder[idx]
      flip_disorder <- flip_disorder[idx]
      flip_eqtl <- flip_eqtl[idx]
      
      m_disorder$BETA <- flip_disorder * m_disorder$BETA
      m_eqtl$BETA  <- flip_eqtl  * m_eqtl$BETA
      
      eqtl_n <- tissue_n_lookup[brain_tissue]
      
      if(disorder == "MDD"){
      d_disorder <- list(
        beta    = m_disorder$BETA,
        varbeta = (m_disorder$SE)^2,
        snp     = snps,
        position = m_disorder$POS, 
        type    = "cc",
        N       = MDD_n,
        MAF     = pmin(m_disorder$EAF, 1 - m_disorder$EAF), 
        LD      = R,
        s = MDD_s
      )
      } else{
        d_disorder <- list(
          beta    = m_disorder$BETA,
          varbeta = (m_disorder$SE)^2,
          snp     = snps,
          position = m_disorder$POS, 
          type    = "quant",
          N       = MCP_n,
          MAF     = pmin(m_disorder$EAF, 1 - m_disorder$EAF), 
          LD      = R
        )
      }
      
      d_eqtl <- list(
        beta    = m_eqtl$BETA,
        varbeta = (m_eqtl$SE)^2,
        p = m_eqtl$P,
        snp     = snps,
        position = m_disorder$POS, 
        type    = "quant",
        N       = eqtl_n,
        MAF     = pmin(m_eqtl$EAF, 1 - m_eqtl$EAF), 
        LD      = R
      )

      # z <- m_disorder$BETA/m_disorder$SE
      # condz = kriging_rss(z, R, n=MCP_n)
      # condz$plot
      # 
      # z <- m_eqtl$BETA/m_eqtl$SE
      # condz = kriging_rss(z, R, n=715)
      # condz$plot
      
      ### Run SuSiE ----
      fit_disorder <- runsusie(d_disorder, maxit = 1000, L = 10)
      fit_eqtl  <- runsusie(d_eqtl, maxit = 1000, L = 10)
      
      ### Run coloc ----
      if(is.null(fit_disorder$sets$cs_index) | is.null(fit_eqtl$sets$cs_index)){
        
        method_used <- "abf"
        
        if(disorder == "MDD"){
          d1_abf <- d_disorder[c("beta","varbeta","snp","type","N","MAF", "s")]
        } else{
          d1_abf <- d_disorder[c("beta","varbeta","snp","type","N","MAF")]
        }
        d2_abf <- d_eqtl[c("beta","varbeta","snp","type","N","MAF")]
        cres <- coloc.abf(d1_abf, d2_abf)
        cres$summary <- as.data.frame(t(cres$summary))
        
      } else{
        
        method_used <- "susie"
        cres <- coloc.susie(fit_disorder, fit_eqtl)
      }
      
      coloc_results[[paste0(brain_tissue,"_", gene_name,"_", disorder)]] <- list(
        fits = list(disorder = fit_disorder, eQTL = fit_eqtl),
        coloc = list(full = cres, summary = cres$summary, method = method_used),
        meta  = list(snps = snps, n_vars = length(snps))
      )
    }
  }
}

# Save outputs ----
saveRDS(coloc_results, "./output/coloc/eQTL.rds")

# Format results ----
parse_key <- function(k) {
  parts <- strsplit(k, "_", fixed = TRUE)[[1]]
  n <- length(parts)
  list(
    brain_tissue = if (n > 2) paste(parts[1:(n-2)], collapse = "_") else NA_character_,
    gene         = if (n > 1) parts[n-1] else NA_character_,
    disorder     = parts[n]
  )
}

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
  parts <- parse_key(k)
  
  fit_dis  <- x$fits$disorder
  fit_eqtl <- x$fits$eQTL
  
  if (m == "susie") {
    ss <- x$coloc$summary
    
    for (i in seq_len(nrow(ss))) {
      cs1 <- ss$idx1[i]; cs2 <- ss$idx2[i]
      cs1_idx <- fit_dis$sets$cs[[paste0("L", cs1)]]
      cs2_idx <- fit_eqtl$sets$cs[[paste0("L", cs2)]]
      
      top1 <- cs1_idx[ which.max(fit_dis$pip[cs1_idx]) ]
      top2 <- cs2_idx[ which.max(fit_eqtl$pip[cs2_idx]) ]
      
      out[[length(out)+1L]] <- data.frame(
        Tissue             = parts$brain_tissue,  
        Disorder           = parts$disorder,
        Gene               = parts$gene,
        `Coloc appraoch`       = m,
        PP4                = ss$PP.H4.abf[i],
        
        
        `Disorder credible set (n)`   = length(cs1_idx),
        `Disorder credible set SNPs`   = paste(x$meta$snps[cs1_idx], collapse=";"),
        `Top Disorder SNP`   = x$meta$snps[top1],
        `eQTL credible set (n)`   = length(cs2_idx),
        `eQTL credible set SNPs`   = paste(x$meta$snps[cs2_idx], collapse=";"),
        `Top eQTL SNP`   = x$meta$snps[top2],
        stringsAsFactors   = FALSE,
        check.names = FALSE
      )
    }
    
  } else if (m == "abf") {
    
    pp4 <- pph4_from_summary(x$coloc$summary)
    res <- x$coloc$full$results
    ord <- order(-res$SNP.PP.H4)
    top_shared     <- res$snp[ord][1]
    
    out[[length(out)+1L]] <- data.frame(
      Tissue             = parts$brain_tissue,
      Disorder           = parts$disorder,
      Gene               = parts$gene,
      `Coloc appraoch`       = m,
      PP4                = pp4,
      `Disorder credible set (n)`   = NA_integer_,
      `Disorder credible set SNPs`   = NA_integer_,
      `Top Disorder SNP`   = top_shared,
      `eQTL credible set (n)`   = NA_integer_,
      `eQTL credible set SNPs`   = NA_integer_,
      `Top eQTL SNP`   = top_shared,
      stringsAsFactors   = FALSE,
      check.names = FALSE
    )
  }
}

coloc_full <- if (length(out)) do.call(rbind, out) else
  data.frame(
    Tissue = character(),
    Disorder = character(),
    Gene = character(),
    `Coloc appraoch` = character(),
    PP4 = numeric(),
    `Disorder credible set (n)` = character(),
    `Disorder credible set SNPs` = character(),
    `Top Disorder SNP` = character(),
    `eQTL credible set (n)`= character(),
    `eQTL credible set SNPs`= character(),
    `Top eQTL SNP`= character(),
    stringsAsFactors   = FALSE,
    check.names = FALSE
  )

coloc_full <- arrange(coloc_full, Gene, Tissue, Disorder)
coloc_full$PP4 <- round(coloc_full$PP4, digits =3)
  
coloc_paper <- coloc_full %>%
  transmute(Tissue, Gene, disorder = tolower(Disorder), PP4) %>%
  group_by(Tissue, Gene, disorder) %>%
  summarise(PP4 = max(PP4, na.rm = TRUE), .groups = "drop") %>%
  mutate(disorder = recode(disorder, mdd = "PP4 MDD", mcp = "PP4 MCP")) %>%
  tidyr::pivot_wider(names_from = disorder, values_from = PP4) %>%
  dplyr::filter(dplyr::if_any(dplyr::starts_with("PP4"), ~ .x > 0.80)) %>%
  dplyr::arrange(Gene, Tissue)

write.csv(coloc_full, "./output/coloc/cs_eQTL.csv", row.names = F)
write.csv(coloc_paper, "./output/coloc/cs_eQTL_in_text.csv", row.names = F)

