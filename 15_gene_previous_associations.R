library(httr)
library(jsonlite)

# Gene list ----
genes <- read.csv("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/genes.csv", header = T, check.names = F)[2:3]

# Get EFO ID ----
get_efo <- function(term, n = 5) {
  q <- 'query($q:String!, $n:Int!){
    search(queryString:$q, entityNames:["disease"], page:{index:0,size:$n}) {
      hits { id name }
    }
  }'
  r <- POST("https://api.platform.opentargets.org/api/v4/graphql",
            body = list(query = q, variables = list(q = term, n = n)),
            encode = "json", timeout(30))
  stop_for_status(r)
  hits <- fromJSON(content(r, "text", encoding = "UTF-8"), flatten = TRUE)$data$search$hits
  if (is.null(hits) || NROW(hits) == 0) return(data.frame(query = term, name = NA, id = NA))
  
  # require the input term to appear in the disease *name*
  term_clean  <- gsub("\\s+", " ", trimws(term))
  names_clean <- gsub("\\s+", " ", hits$name)
  keep <- grepl(term_clean, names_clean, ignore.case = TRUE, fixed = TRUE)
  
  out <- data.frame(query = term, name = hits$name[keep], id = hits$id[keep], row.names = NULL)
  if (NROW(out) == 0) out <- data.frame(query = term, name = NA, id = NA)
  out
}

MDD_efo  <- unique(c(get_efo("depressive", n = 1000)$id, get_efo("depression", n = 1000)$id))
pain_efo <- unique(c(get_efo("pain", n = 1000)$id, 
                     get_efo("migraine", n = 1000)$id, 
                     get_efo("gout", n = 1000)$id, 
                     get_efo("colitis", n = 1000)$id, 
                     get_efo("arthritis", n = 1000)$id, 
                     get_efo("fibromyalgia", n = 1000)$id))

# Gene query ----
ot <- function(query, variables) {
  r <- POST("https://api.platform.opentargets.org/api/v4/graphql",
            body = list(query = query, variables = variables),
            encode = "json", timeout(30))
  stop_for_status(r)
  fromJSON(content(r, "text", encoding = "UTF-8"), flatten = TRUE)$data
}

fetch_assoc_rows <- function(ensg, page_size = 500, max_pages = 50) {
  q <- 'query($ensg:String!,$i:Int!,$n:Int!){
    target(ensemblId:$ensg){
      associatedDiseases(
        page:{index:$i,size:$n},
        enableIndirect:false,                                    # ← no ontology roll-up
        datasources:[{id:"ot_genetics", weight:1, propagate:false}]  # ← genetics only, no propagation
      ){
        count
        rows {
          score
          disease { id name }
        }
      }
    }
  }'
  acc <- NULL
  for (i in 0:(max_pages - 1)) {
    d <- ot(q, list(ensg = ensg, i = i, n = page_size))
    rows <- d$target$associatedDiseases$rows
    if (is.null(rows) || NROW(rows) == 0) break
    
    raw_n <- if (is.data.frame(rows)) nrow(rows) else length(rows)
    
    # keep only genetics associations with score > 0
    if (is.data.frame(rows)) {
      rows <- subset(rows, !is.na(score) & score > 0)
    } else {
      rows <- Filter(function(x) !is.null(x$score) && x$score > 0, rows)
    }
    
    if (NROW(rows) > 0) acc <- rbind(acc, rows)
    if (raw_n < page_size) break
  }
  if (is.null(acc)) return(NULL)
  acc <- acc[!duplicated(acc$disease.id), , drop = FALSE]
  acc
}

has_gwas_ev <- function(ensg, efo) {
  q <- 'query($efo:String!, $ensg:String!){
    disease(efoId:$efo){
      evidences(ensemblIds:[$ensg], datasourceIds:["ot_genetics_portal"], size:1){ count }
    }
  }'
  d <- ot(q, list(efo = efo, ensg = ensg))
  !is.null(d$disease$evidences$count) && d$disease$evidences$count > 0
}

# Check gene - disease associations ----
gene_disease_associations <- data.frame(
  gene = character(0),
  depression = logical(0),
  pain = logical(0),
  stringsAsFactors = FALSE
)

for (gene_index in 1:nrow(genes)){
  
  ensembl_id <- genes$ensembl_gene_id[gene_index]
  symbol <- genes$SYMBOL[gene_index]
  
  rows <- fetch_assoc_rows(ensg = ensembl_id, page_size = 500,  max_pages = 10)
  assoc_ids <- if (is.null(rows)) {
    character(0)
  } else if (is.data.frame(rows)) {
    rows$disease.id
  } else {
    vapply(rows, function(x) x$disease$id, character(1))
  }
  
  depression_found <- any(assoc_ids %in% MDD_efo)
  pain_found       <- any(assoc_ids %in% pain_efo)
  
  gene_disease_associations <- rbind(
    gene_disease_associations,
    data.frame(gene = symbol,
               depression = depression_found,
               pain = pain_found,
               stringsAsFactors = FALSE)
  )
  
}

