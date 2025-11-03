library(dplyr)
library(data.table)
library(purrr)

# Load conjFDR results
MDD_MCP <- readRDS("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_AnalysisObject.rds")
MDD_MCP_res <- MDD_MCP$get_trait_results()
FUMA_loci <- read.csv("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.loci.csv", sep = "\t")
FUMA_lead <- read.csv("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.lead.csv", sep = "\t")

## Load GWAS
MDD_GWAS <- fread("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MDD/pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv", sep = "\t")
MCP_GWAS <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MCP/ChronicPain.stats", header = T)

## Map leads to locus
lead_to_locus <- with(FUMA_lead, setNames(locusnum, LEAD_SNP))
leads_by_locus <- tapply(names(lead_to_locus), lead_to_locus, function(x) sort(unique(x)))

## Make loci table
loci_table <- FUMA_loci %>%
  dplyr::select(locusnum, CHR, MinBP, MaxBP, LEAD_SNP) %>%
  dplyr::rename(`Locus number conjFDR` = locusnum,
                Chromosome = CHR,
                `Min pos` = MinBP,
                `Max pos` = MaxBP,
                `Lead SNP` = LEAD_SNP)

loci_table$`Additional Lead SNP(s)` <- vapply(seq_len(nrow(loci_table)), function(i) {
  locus  <- loci_table$`Locus number conjFDR`[i]
  this_lead <- loci_table[["Lead SNP"]][i]
  vec <- leads_by_locus[[as.character(locus)]]
  if (is.null(vec)) return(NA_character_)
  others <- setdiff(vec, this_lead)
  if (length(others)) paste(others, collapse = ";") else NA_character_
}, character(1))


loci_table <- loci_table %>%
  left_join(dplyr::select(FUMA_lead, LEAD_SNP, LEAD_BP, FDR), by = c("Lead SNP" = "LEAD_SNP")) %>%
  dplyr::rename(pos = LEAD_BP)


loci_table <- loci_table %>%
  left_join(MDD_MCP_res, by = c("Lead SNP" = "SNP")) %>%
  dplyr::mutate(`MDD p` = 10^(-.data$LOG10PVAL1),
         `MCP p` = 10^(-.data$LOG10PVAL2)) %>%
  dplyr::rename(`MDD beta` = BETA1,
         `MCP beta` = BETA2,
         `conj FDR` = FDR)

out <- loci_table %>%
  left_join(MDD_GWAS, by = c("Lead SNP" = "ID")) %>%
  mutate(
    A1 = EA,
    A2 = NEA,
  ) %>%
  mutate(across(c(`MDD beta`, `MCP beta`), ~ round(.x, 3))) %>%
  dplyr::select(
    `Locus number conjFDR`, Chromosome, `Min pos`, `Max pos`, `Lead SNP`, `Additional Lead SNP(s)`,  # locus info
    pos, A1, A2, `MDD beta`, `MDD p`, `MCP beta`, `MCP p`, `conj FDR`,   # lead info
  )

write.csv(out, "~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/loci_formatted.csv", row.names = F)
  