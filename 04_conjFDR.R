library(data.table)
# remotes::install_github("alexploner/cfdr.pleio", build_vignettes = TRUE)
library(cfdr.pleio)
library(dplyr)

# MDD GWAS ----
file <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MDD//pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv"
MDD_GWAS <- fread(file, sep = "\t")

MDD_minimal <- data.frame(
  SNP = MDD_GWAS$ID,
  BETA = MDD_GWAS$BETA,
  PVAL   = MDD_GWAS$PVAL,
  BP = MDD_GWAS$POS,
  CHR = MDD_GWAS$CHROM
)

setDT(MDD_minimal)
MDD_minimal <- MDD_minimal[!(CHR == 6 & BP >= 25e6 & BP <= 35e6)]
MDD_minimal <- MDD_minimal[!(CHR == 8 & BP >= 7.2e6 & BP <= 12.5e6)]

# MCP GWAS ----
MCP_GWAS <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MCP/ChronicPain.stats", header = T)

## Harmonize with MDD
merged <- inner_join(MCP_GWAS, MDD_GWAS, by = c("SNP" = "ID"))
flip <- merged$ALLELE1 != merged$EA



MCP_minimal <- data.frame(
  SNP = MCP_GWAS$SNP,
  BETA = MCP_GWAS$BETA,
  PVAL   = MCP_GWAS$P,
  BP = MCP_GWAS$BP,
  CHR = MCP_GWAS$CHR
)

setDT(MCP_minimal)
MCP_minimal <- MCP_minimal[!(CHR == 6 & BP >= 25e6 & BP <= 35e6)]
MCP_minimal <- MCP_minimal[!(CHR == 8 & BP >= 7.2e6 & BP <= 12.5e6)]

# conjFDR ---- 
#download.file("https://zenodo.org/record/5750318/files/genref.zip",
#             destfile = "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/genref.zip", mode = "wb")
#unzip("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/genref.zip", exdir = "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/genref")
REF_DIR <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/genref"

MDD_MCP <- cfdr_pleio$new()

MDD_MCP$init_data(trait1 = MDD_minimal,
                  trait2 = MCP_minimal, 
                  trait_names = c("MDD", "MCP"),
                  refdat = refdata_location(REF_DIR), 
                  local_refdat_path = "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/",
                  verbose = T)

MDD_MCP$initialize_pruning_index(n_iter = 50, seed = 151097, verbose = FALSE)

MDD_MCP$calculate_cond_fdr(fdr_trait = 1, verbose = FALSE)
MDD_MCP$calculate_cond_fdr(fdr_trait = 2, verbose = FALSE)

MDD_MCP_res <- MDD_MCP$get_trait_results()

saveRDS(MDD_MCP, file = paste("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_AnalysisObject.rds", sep="/"))
write.table(MDD_MCP_res, "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_res.txt", row.names = F, quote = F)
