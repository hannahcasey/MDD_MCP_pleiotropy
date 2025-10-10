library(ghql)
library(jsonlite)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(data.table)
library(biomaRt)
library(tidyr)

# Load credible sets ----
cs <- read.csv("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/cs_cross_trait_significant.csv", header = T)

.ot_run_list <- function(query, variables = list()) {
  con <- ghql::GraphqlClient$new("https://api.platform.opentargets.org/api/v4/graphql")
  qry <- ghql::Query$new(); qry$query("q", query)
  fromJSON(con$exec(qry$queries$q, variables), simplifyVector = FALSE)
}

`%||%` <- function(x, y) if (!is.null(x)) x else y

# rsID -> Open Targets variantId ----
rsids_to_variant_ids <- function(rsids, batch_size = 200L) {
  rsids <- unique(na.omit(rsids))
  if (!length(rsids)) return(tibble(rsid = character(), variantId = character()))
  q <- '
  query($terms:[String!]!){
    mapIds(queryTerms:$terms, entityNames:["variant"]) {
      mappings {
        term
        hits { object { __typename ... on Variant { id rsIds } } }
      }
    }
  }'
  batches <- split(rsids, ceiling(seq_along(rsids)/batch_size))
  out <- map(batches, function(bs) {
    res  <- .ot_run_list(q, list(terms = bs))
    maps <- res$data$mapIds$mappings %||% list()
    map_dfr(maps, function(m) {
      term <- m$term %||% NA_character_
      hits <- m$hits %||% list()
      cand <- keep(hits, ~ !is.null(.x$object$rsIds) && term %in% .x$object$rsIds)
      if (!length(cand)) cand <- hits
      variantId <- if (!length(cand) || is.null(cand[[1]]$object$id)) NA_character_ else cand[[1]]$object$id
      tibble(rsid = term, variantId = variantId)
    })
  }) %>% list_rbind()
  full_join(tibble(rsid = rsids), out, by = "rsid")
}

rsIDs <- na.omit(c(cs$top_MCP_snp, unlist(str_split(cs$MCP_CS_snps, ";")), unlist(str_split(cs$MDD_CS_snps, ";"))))
variant_ids <- rsids_to_variant_ids(unique(rsIDs))
temp <- strsplit(variant_ids$variantId, "_")
variant_ids$chr <- as.numeric(vapply(temp, `[`, character(1), 1))
variant_ids$pos <- as.numeric(vapply(temp, `[`, character(1), 2))

# Join with V2G scores - (REF script) ----
v2g_csv   <- normalizePath("~/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/v2g_top_gene_per_variant.csv",
                           winslash = "/", mustWork = TRUE)
pairs_csv <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/v2g_pairs.csv"

pairs <- variant_ids %>%
  mutate(chr = sub("^chr","", as.character(chr), ignore.case = TRUE),
         pos = as.integer(pos)) %>%
  distinct(chr, pos)

fwrite(pairs, pairs_csv, col.names = FALSE)

cmd <- sprintf(
  'awk -F, \'NR==FNR{a[$1","$2]=1; next} FNR==1 || a[$1","$2]\' "%s" "%s"',
  pairs_csv, v2g_csv
)

v2g_scores <- fread(
  cmd = cmd,
  select = c("chr_id","position","ref_allele","alt_allele","gene_id","overall_score")
)

V2G_long <- left_join(variant_ids, v2g_scores, by = c("chr" = "chr_id", "pos" = "position"))

V2G_long <- V2G_long %>%
  rename(V2G_score = overall_score) %>%
  dplyr::select(rsid, variantId, gene_id, V2G_score)

# Format credible sets ----
key_tbl_all <- V2G_long %>%
  group_by(rsid) %>%
  summarise(genes = list(unique(gene_id)), .groups = "drop")

rsid2genes_all <- setNames(key_tbl_all$genes, key_tbl_all$rsid)

split_rs <- function(x) {
  if (is.na(x) || !nzchar(x)) character(0)
  else str_split(x, ";", simplify = FALSE)[[1]] |> str_trim() |> discard(~ .x == "")
}

clump.lead <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.lead.csv", header = T)
clump.lead <- dplyr::select(clump.lead, CHR, LEAD_SNP, LEAD_BP) 

cs_annot <- cs %>%
  mutate(
    MCP_rs = map(MCP_CS_snps, split_rs),
    MDD_rs = map(MDD_CS_snps, split_rs),
    all_rs = pmap(list(MCP_rs, MDD_rs, top_MDD_snp), function(mcp, mdd, topd) {
      vec <- c(mcp, mdd, if (!is.null(topd) && !is.na(topd) && nzchar(topd)) topd else character(0))
      unique(vec)
    }),
    v2g_genes = map(all_rs, function(rs_vec) {
      genes_list <- rsid2genes_all[rs_vec]
      genes <- unique(unlist(genes_list, use.names = FALSE))
      genes <- genes[!is.na(genes)]
      paste(genes, collapse = ";")
    }),
    v2g_genes = as.character(v2g_genes)
  ) %>%
  left_join(clump.lead, by = c("Lead_SNP" = "LEAD_SNP")) %>%
  rename(chr = CHR, pos = LEAD_BP) %>%
  dplyr::select(Lead_SNP, chr, pos, coloc_method, PP4, MCP_CS_size, MDD_CS_size,
         top_MCP_snp, top_MDD_snp, v2g_genes) %>%
  arrange(coloc_method)

## Get gene names ----
gene_list_long <- V2G_long$gene_id
mart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
gene_coords <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = gene_list_long,
  mart = mart
)

V2G_long <- left_join(V2G_long, gene_coords, by = c("gene_id" = "ensembl_gene_id"))

cs_annot <- cs_annot %>%
  mutate(row_id = dplyr::row_number(),
         v2g_genes = na_if(v2g_genes, "")) %>%
  separate_rows(v2g_genes, sep = ";") %>%      # splits "A;B"
  mutate(v2g_genes_clean = sub("\\.\\d+$", "", v2g_genes)) %>%    # drop .version
  left_join(gene_coords, by = c("v2g_genes_clean" = "ensembl_gene_id")) %>%
  group_by(row_id) %>%
  summarise(
    across(c(Lead_SNP, chr, pos, coloc_method, PP4, top_MCP_snp, top_MDD_snp), dplyr::first),
    v2g_genes   = paste(unique(na.omit(v2g_genes)), collapse = ";"),
    v2g_symbols = paste(unique(na.omit(hgnc_symbol)), collapse = ";"),
    .groups = "drop"
  ) %>%
  dplyr::select(-row_id) %>%
  arrange(chr, pos)

# Save outputs ----
write.csv(V2G_long, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/V2G/v2g_long.csv", row.names = FALSE)
write.csv(cs_annot, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/V2G/cs_cross_trait_in_text.csv", row.names = FALSE)
