library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(biomaRt)
library(Rsamtools)
library(GenomicRanges)
library(purrr)
library(haploR)
library(data.table)

# Load conjFDR results
MDD_MCP <- readRDS("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_AnalysisObject.rds")
MDD_MCP_res <- MDD_MCP$get_trait_results()

FUMA_loci <- read.csv("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.loci.csv", sep = "\t")
coloc_loci <- read.csv("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/cs_cross_trait_significant.csv")

coloc_loci <- left_join(coloc_loci, FUMA_loci, by = c("Lead_SNP" = "LEAD_SNP"))

## Load GWAS
MDD_GWAS <- fread("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MDD/pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv", sep = "\t")
MCP_GWAS <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MCP/ChronicPain.stats", header = T)

## Get locus bands
ses <- browserSession("UCSC"); genome(ses) <- "hg19"
cb <- getTable(ucscTableQuery(ses, table = "cytoBand"))
cyto <- GRanges(cb$chrom, IRanges(cb$chromStart + 1, cb$chromEnd), band = cb$name)

get_band <- function(chr, start, end) {
  gr <- GRanges(paste0("chr", chr), IRanges(start, end))
  h  <- findOverlaps(gr, cyto)
  if (length(h) == 0) return(NA_character_)
  b  <- mcols(cyto)$band[subjectHits(h)]
  chr_n <- gsub("^chr","", as.character(seqnames(gr)))
  if (length(b) == 1) paste0(chr_n, b) else paste0(chr_n, b[1], "–", b[length(b)])
}

bands <- coloc_loci %>%
  mutate(Band = mapply(get_band, CHR, MinBP, MaxBP))

## SNP annotation
## Consequence
get_vep_consequence <- function(rsid) {
  url <- paste0("https://rest.ensembl.org/vep/human/id/", rsid, "?content-type=application/json")
  r <- tryCatch(httr::GET(url), error = function(e) NULL); if (is.null(r) || r$status_code != 200) return(NA_character_)
  jj <- jsonlite::fromJSON(rawToChar(r$content))
  msc <- jj$most_severe_consequence
  return(msc)
}

## CADD score
get_vep_cadd <- function(rsid) {
  url <- sprintf(
    "https://rest.ensembl.org/vep/human/id/%s?content-type=application/json&CADD=1",
    rsid
  )
  r <- httr::GET(
    url   = "https://rest.ensembl.org/vep/human/id",
    query = list(id = rsid, CADD = "1")
  )
  
  if (is.null(r) || httr::status_code(r) != 200) {
    return(data.frame(id = rsid, most_severe_consequence = NA, cadd_raw = NA, cadd_phred = NA))
  }
  x <- jsonlite::fromJSON(rawToChar(r$content), simplifyDataFrame = F)
  rec <- x[[1]]
  
  # CADD can appear on various consequence blocks; take the max PHRED across any we see
  take_max <- function(df, col) if (!is.null(df) && col %in% names(df)) suppressWarnings(max(as.numeric(df[[col]]), na.rm = TRUE)) else NA_real_
  phreds <- c(
    suppressWarnings(as.numeric(rec$CADD_PHRED)),
    take_max(rec$transcript_consequences, "CADD_PHRED"),
    take_max(rec$regulatory_feature_consequences, "CADD_PHRED"),
    take_max(rec$motif_feature_consequences, "CADD_PHRED"),
    take_max(rec$intergenic_consequences, "CADD_PHRED")
  )
  raws <- c(
    suppressWarnings(as.numeric(rec$CADD_RAW)),
    take_max(rec$transcript_consequences, "CADD_RAW"),
    take_max(rec$regulatory_feature_consequences, "CADD_RAW"),
    take_max(rec$motif_feature_consequences, "CADD_RAW"),
    take_max(rec$intergenic_consequences, "CADD_RAW")
  )
  data.frame(
    id = rec$id %||% rec$input %||% rsid,
    most_severe_consequence = rec$most_severe_consequence %||% NA_character_,
    cadd_raw = if (all(is.na(raws))) NA_real_ else max(raws, na.rm = TRUE),
    cadd_phred = if (all(is.na(phreds))) NA_real_ else max(phreds, na.rm = TRUE)
  )
}

## Regulomedb rank
get_reg <- function(rsid) {
  reg <- haploR::queryRegulome(query = rsid, genomeAssembly = "GRCh37")
  return(unique(reg[[rsid]]$regulome_score[,2]))
}

top_MCP_SNPs <- coloc_loci$top_MCP_snp
top_MDD_SNPs <- coloc_loci$top_MDD_snp

top_MCP_SNPs_ann <- lapply(top_MCP_SNPs, function(rs) {
  conseq <- get_vep_consequence(rs)
  cr <- round(cadd_scores$CADD_PHRED_max[cadd_scores$rsid == rs], digits = 2)
  rb <- get_reg(rs)
  
  tibble(
    rsid = rs,
    conseq = conseq,
    cr = cr,
    rb = rb
  )
}) |> dplyr::bind_rows()

top_MDD_SNPs_ann <- lapply(top_MDD_SNPs, function(rs) {
  conseq <- get_vep_consequence(rs)
  cr <- round(cadd_scores$CADD_PHRED_max[cadd_scores$rsid == rs], digits = 2)
  rb <- get_reg(rs)
  
  tibble(
    rsid = rs,
    conseq = conseq,
    cr = cr,
    rb = rb
  )
}) |> dplyr::bind_rows()

## Format table
annotated <- bands %>%
  dplyr::select(CHR, LEAD_BP, Band, MinBP, MaxBP, PP4, top_MDD_snp, top_MCP_snp)

for(row in 1:nrow(annotated)){
  if(annotated[row,"top_MDD_snp"] == annotated[row,"top_MCP_snp"]){
    annotated$conseq[row] <- top_MDD_SNPs_ann$conseq[row]
    annotated$cr[row] <-  top_MDD_SNPs_ann$cr[row]
    annotated$rb[row] <- top_MDD_SNPs_ann$rb[row]
  } else{
    annotated$conseq[row] <- paste0(top_MDD_SNPs_ann$conseq[row], ", ", top_MCP_SNPs_ann$conseq[row])
    annotated$cr[row] <-  paste0(top_MDD_SNPs_ann$cr[row], ", ", top_MCP_SNPs_ann$cr[row])
    annotated$rb[row] <- paste0(top_MDD_SNPs_ann$rb[row], ", ", top_MCP_SNPs_ann$rb[row])
  }
}

write.csv(annotated, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/cs_cross_trait_significant_annotated.csv", row.names = F)
