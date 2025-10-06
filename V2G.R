library(ghql)
library(jsonlite)
library(dplyr)
library(purrr)
library(tibble)
library(stringr)
library(data.table)

# Load credible sets ----
cs <- read.csv("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/coloc/cs.csv", header = T)

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
  select(rsid, variantId, gene_id, V2G_score)

# Format credible sets ----
key_tbl_all <- V2G_long %>%
  group_by(rsid) %>%
  summarise(genes = list(unique(gene_id)), .groups = "drop")

rsid2genes_all <- setNames(key_tbl_all$genes, key_tbl_all$rsid)

split_rs <- function(x) {
  if (is.na(x) || !nzchar(x)) character(0)
  else str_split(x, ";", simplify = FALSE)[[1]] |> str_trim() |> discard(~ .x == "")
}

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
  select(-c("MCP_rs", "MDD_rs", "all_rs")) %>%
  arrange(coloc_method)

# Save outputs ----
write.csv(V2G_long, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/V2G/v2g_long.csv", row.names = FALSE)
write.csv(cs_annot, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/V2G/cs_with_v2g_genes.csv", row.names = FALSE)
