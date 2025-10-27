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

FUMA_loci <- read.csv("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/loci_formatted.csv", check.names = FALSE)
coloc_loci <- read.csv("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/cs_cross_trait_significant.csv", check.names = FALSE)

coloc_loci <- left_join(coloc_loci, FUMA_loci, by = c("Lead SNP", "Max pos", "Min pos", "Chromosome"))

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
  mutate(Band = mapply(get_band, Chromosome, `Min pos`, `Max pos`))

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
  r <- httr::GET(url)
  
  x <- jsonlite::fromJSON(rawToChar(r$content), simplifyDataFrame = F)
  rec <- x[[1]]
  
  tc <- rec$transcript_consequences
  if (is.null(tc) || length(tc) == 0) tc <- rec$intergenic_consequences
    
  ph <- numeric(0)

  for (d in tc) {
      ph <- c(ph, d$cadd_phred)
  }
  
  return(max(ph))
  
}

## Regulomedb rank
get_reg <- function(rsid) {
  reg <- haploR::queryRegulome(query = rsid, genomeAssembly = "GRCh37")
  return(unique(reg[[rsid]]$regulome_score[,2]))
}

top_MDD_SNPs <- coloc_loci$`Top MDD SNP`

top_MCP_SNPs <- coloc_loci$`Top MCP SNP`

top_MDD_SNPs_ann <- lapply(top_MDD_SNPs, function(rs) {
  conseq <- get_vep_consequence(rs)
  cr <- round(get_vep_cadd(rs), digits = 2)
  rb <- get_reg(rs)
  
  tibble(
    rsid = rs,
    conseq = conseq,
    cr = cr,
    rb = rb
  )
}) |> dplyr::bind_rows()

top_MCP_SNPs_ann <- lapply(top_MCP_SNPs, function(rs) {
  conseq <- get_vep_consequence(rs)
  cr <- round(get_vep_cadd(rs), digits = 2)
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
  dplyr::select(Band,"Min pos", "Max pos", PP4, "Top MDD SNP", "Top MCP SNP")

for(row in 1:nrow(annotated)){
  if(annotated[row,"Top MDD SNP"] == annotated[row,"Top MCP SNP"]){
    annotated$`Most severe consequence`[row] <- top_MDD_SNPs_ann$conseq[row]
    annotated$CADD[row] <-  top_MDD_SNPs_ann$cr[row]
    annotated$RegulomeDB[row] <- top_MDD_SNPs_ann$rb[row]
  } else{
    
    if(top_MDD_SNPs_ann$conseq[row] == top_MCP_SNPs_ann$conseq[row]){
      annotated$`Most severe consequence`[row] <- top_MDD_SNPs_ann$conseq[row]
    }else{
      annotated$`Most severe consequence`[row] <- paste0(top_MDD_SNPs_ann$conseq[row], ", ", top_MCP_SNPs_ann$conseq[row])
    }
    
    if(top_MDD_SNPs_ann$cr[row] == top_MCP_SNPs_ann$cr[row]){
      annotated$CADD[row] <- top_MDD_SNPs_ann$cr[row]
    }else{
      annotated$CADD[row] <- paste0(top_MDD_SNPs_ann$cr[row], ", ", top_MCP_SNPs_ann$cr[row])
    }
    
    if(top_MDD_SNPs_ann$rb[row] == top_MCP_SNPs_ann$rb[row]){
      annotated$RegulomeDB[row] <- top_MDD_SNPs_ann$rb[row]
    }else{
      annotated$RegulomeDB[row] <- paste0(top_MDD_SNPs_ann$rb[row], ", ", top_MCP_SNPs_ann$rb[row])
    }
  }
}

write.csv(annotated, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/cs_cross_trait_significant_annotated.csv", row.names = F)

