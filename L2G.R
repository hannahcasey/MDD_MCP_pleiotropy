library(ghql)
library(jsonlite)
library(dplyr)
library(purrr)
library(tibble)

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

res <- ensembl_fetch_variant(c("rs2237303"))
variant_id <- with(res, paste(seq_region_name, start, gsub("/", "_", allele_string), sep = "_"))


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

cli <- GraphqlClient$new(url = "https://api.platform.opentargets.org/api/v4/graphql")

res_json <- cli$exec(qry$queries$QTLCredibleSetsQuery,
                     variables = list(vid = variant_id))
dat <- fromJSON(res_json, simplifyVector = TRUE)
dat$data$variant$credibleSets$rows


lst <- dat$data$variant$credibleSets$rows$l2GPredictions$rows

df <- bind_rows(lapply(seq_along(lst), function(i) {
  x <- lst[[i]]
  if (is.null(x) || NROW(x) == 0) return(NULL)
  dplyr::mutate(tibble::as_tibble(x), cs_idx = i, .before = 1)
}))

