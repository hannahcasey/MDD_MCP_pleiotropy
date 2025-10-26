library(dplyr)
library(data.table)

# Load conjFDR results
MDD_MCP <- readRDS("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_AnalysisObject.rds")
MDD_MCP_res <- MDD_MCP$get_trait_results()
FUMA_loci <- read.csv("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.loci.csv", sep = "\t")

## Load GWAS
MDD_GWAS <- fread("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MDD/pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv", sep = "\t")
MCP_GWAS <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MCP/ChronicPain.stats", header = T)






Chromosome 
Min
Max
Lead SNP
other SNP

pos
bp
A1
A2

beta_1
P_1

beta_2
P_2

condfdr