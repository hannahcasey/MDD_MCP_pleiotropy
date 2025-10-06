library(dplyr)

conjFDR.clump.lead <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.lead.csv", header = T)

MDD_MCP <- readRDS("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_AnalysisObject.rds")
MDD_MCP_res <- MDD_MCP$get_trait_results()

fuma_formatted <- conjFDR.clump.lead %>%
  left_join(MDD_MCP_res, by = c("LEAD_SNP" = "SNP")) %>%
  filter(sign(BETA1) == sign(BETA2)) %>%
  select(LEAD_SNP, CHR, LEAD_BP)

nrow(conjFDR.clump.lead)
nrow(fuma_formatted)

write.table(fuma_formatted, "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.lead.fuma.txt", quote = F, row.names = F, col.names = T)
