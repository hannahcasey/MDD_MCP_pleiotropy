library(ghql)
library(jsonlite)
library(dplyr)
library(purrr)
library(tibble)
library(rtracklayer)
library(GenomicRanges)
library(stringr)
library(tidyr)
library(biomaRt)

## Load data ----
### Cross-trait coloc results ----
cs_cross_trait_significant <- read.csv("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/cs_cross_trait_significant.csv", header = T, check.names = F)

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

bands <- cs_cross_trait_significant %>%
  mutate(Band = mapply(get_band, Chromosome, `Min pos`, `Max pos`))

## Get variant ID
rsIDs <- na.omit(c(cs_cross_trait_significant$`Top MDD SNP`, unlist(str_split(cs_cross_trait_significant$`MDD credible set SNPs`, ";")), unlist(str_split(cs_cross_trait_significant$`MCP credible set SNPs`, ";"))))

ensembl_fetch_variant <- function(
    rsids,
    endpoint = "http://rest.ensembl.org/variation/homo_sapiens"
) {
  cli::cli_alert_info(sprintf("Fetching rsIds from %s", endpoint))
  httr::POST(
    endpoint,
    httr::content_type("application/json"),
    httr::accept("application/json"),
    body = jsonlite::toJSON(list(ids = rsids))
  ) |>
    httr::stop_for_status() |>
    httr::content() |>
    jsonlite::toJSON() |>
    jsonlite::fromJSON() |>
    ensembl_parse_response()
}

ensembl_parse_response <- function(resp) {
  cli::cli_alert_info("Parsing response from Ensembl API")
  if (length(resp) == 0) {
    cli::cli_alert_danger("No variants found for the provided rsIDs.")
    return(data.frame())
  }
  purrr::map(resp, function(obj) as.data.frame(obj$mapping)) |>
    dplyr::bind_rows() |>
    tidyr::unnest(cols = dplyr::everything()) |>
    dplyr::mutate(
      rsid = names(resp)
    )
}

SNP_ID_ls <- list()

for (rsID in rsIDs){
  
  res <- ensembl_fetch_variant(rsID)
  variant_ids <- with(res, {
    alleles <- strsplit(allele_string, "/", fixed = TRUE)[[1]]
    ref <- alleles[1]
    alts <- alleles[-1]
    paste(seq_region_name, start, ref, alts, sep = "_")
  })
  SNP_ID_ls[[rsID]]$var_ID <- variant_ids
}

## Get L2G ----
cli <- GraphqlClient$new(url = "https://api.platform.opentargets.org/api/v4/graphql")
qry <- Query$new()
qry$query("QTLCredibleSetsQuery", '
  query QTLCredibleSetsQuery($vid: String!) {
    variant(variantId: $vid) {
      id
      credibleSets(studyTypes: [gwas], page: {index: 0, size: 50}) {
        rows {
          studyId
          studyLocusId
          l2GPredictions(page: {index: 0, size: 50}) {
            rows {
              score
              target { id approvedSymbol }
            }
          }
        }
      }
    }
  }')


L2G_ls <- list()

for (rs_index in seq_along(SNP_ID_ls)){
  
  var_IDs <- SNP_ID_ls[[rs_index]]$var_ID
  rsID <- names(SNP_ID_ls)[rs_index]
  
  for (var_ID in var_IDs){
    
    res_json <- cli$exec(qry$queries$QTLCredibleSetsQuery,
                         variables = list(vid = var_ID))
    dat <- fromJSON(res_json, simplifyVector = TRUE)
    
    lst <- dat$data$variant$credibleSets$rows$l2GPredictions$rows
    
    if (!is.null(lst)){
      
      df <- bind_rows(lapply(seq_along(lst), function(i) {
        x <- lst[[i]]
        if (is.null(x) || NROW(x) == 0) return(NULL)
        dplyr::mutate(tibble::as_tibble(x), cs_idx = i, .before = 1)
      }))
      
      L2G_ls[[rsID]] <- df
      
    } 
  }
}

rsid_genes_df <- do.call(rbind, lapply(names(L2G_ls), function(rsid) {
  x <- L2G_ls[[rsid]]
  if(nrow(x) == 0){
    data.frame(rsID = rsid, genes = NA, stringsAsFactors = FALSE)
  } else{
    gene_symbols <- unique(x$target$approvedSymbol)
    data.frame(rsID = rsid, genes = paste(gene_symbols, collapse = ", "), stringsAsFactors = FALSE)
  }
}))
rownames(rsid_genes_df) <- NULL

## Format tables ----
split_semis <- function(x) str_trim(unlist(strsplit(x, ";", fixed = TRUE)))

rs_lookup <- rsid_genes_df %>%
  mutate(genes = ifelse(is.na(genes), "", genes)) %>%
  { setNames(.$genes, .$rsID) }

rs_lookup <- setNames(rsid_genes_df$genes, rsid_genes_df$rsID)

loci_in_text <- bands %>%
  rowwise() %>%
  mutate(
    .all_rs = list(unique(na.omit(c(
      `Top MDD SNP`,
      `Top MCP SNP`,
      split_semis(`MDD credible set SNPs`),
      split_semis(`MCP credible set SNPs`)
    )))),
    Symbol = {
      gene_strs <- unname(rs_lookup[unlist(.all_rs)])
      gene_strs <- gene_strs[!is.na(gene_strs) & nzchar(gene_strs)]
      syms <- trimws(unlist(strsplit(gene_strs, ",", fixed = TRUE), use.names = FALSE))
      syms <- sort(unique(syms[syms != ""]))
      if (length(syms)) paste(syms, collapse = ", ") else NA_character_
    }
  ) %>%
  ungroup() %>%
  dplyr::select(Band, Chromosome, `Min pos`, `Max pos`, `Top MDD SNP`, `Top MCP SNP`,  PP4, Symbol)


rs_map <- cs_cross_trait_significant %>%
  mutate(
    rsIDs = pmap(
      list(`Lead SNP`, `Top MDD SNP`, `Top MCP SNP`,
           `MDD credible set SNPs`, `MCP credible set SNPs`),
      ~ unique(na.omit(c(..1, ..2, ..3, split_semis(..4), split_semis(..5))))
    )
  ) %>%
  transmute(
    `Locus number conjFDR`,
    Chromosome,
    Position = `Lead SNP pos`,
    rsIDs
  ) %>%
  unnest_longer(rsIDs, values_to = "rsID") %>%   # <- long, not wide
  filter(!is.na(rsID), rsID != "") %>%
  distinct(rsID, .keep_all = TRUE)

L2G_long <- rsid_genes_df %>%
  left_join(rs_map, by = "rsID") %>%
  dplyr::rename(SNP = rsID,
         `L2G Genes` = genes) %>%
  dplyr::select(SNP, `Locus number conjFDR`, Chromosome, Position, `L2G Genes`) %>%
  arrange(`Locus number conjFDR`, Chromosome, Position)

## Get gene list ----
genes_hgnc <- sort(unique(unlist(strsplit(rsid_genes_df$genes, ", ", fixed = TRUE))))

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
            filters = "hgnc_symbol",
            values = genes_hgnc,
            mart = mart)

genes <- tibble(SYMBOL = genes_hgnc) %>%
  left_join(distinct(bm, hgnc_symbol, ensembl_gene_id),
            by = c("SYMBOL" = "hgnc_symbol"))


write.csv(loci_in_text, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/L2G/loci_in_text.csv")
write.csv(L2G_long, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/L2G/L2G_long.csv")
write.csv(genes, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/genes.csv")
write.table(genes$ensembl_gene_id, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/gene_list.txt", row.names = F, col.names = T, quote = F)
