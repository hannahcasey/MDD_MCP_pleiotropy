library(data.table)
#remotes::install_github("alexploner/cfdr.pleio", build_vignettes = TRUE)
library(cfdr.pleio)
library(dplyr)
library(tidyr)

# MDD GWAS ----
## MDD (EUR, UKB 23andMe excluded) ----
MDD_GWAS <- fread("./resources/summary_statistics/MDD//pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv", sep = "\t")
## MDD (EHR, EUR) ----
MDD_EHR_GWAS <- fread("./resources/summary_statistics/MDD/pgc-mdd2025_EHR_eur_v3-49-24-11.tsv", sep = "\t")
## MDD (Questionnaire, EUR, UKB excluded) ----
MDD_quest_GWAS <- fread("./resources/summary_statistics/MDD/pgc-mdd2025_Quest-noUKBB_eur_v3-49-24-11.tsv", sep = "\t")

# MCP GWAS ----
MCP_GWAS <- read.table("./resources/summary_statistics/MCP/ChronicPain.stats", header = T)

# Harmonise ----
mdd <- MDD_GWAS %>%
  transmute(
    SNP = ID,                
    EA_MDD = toupper(EA),   
    NEA_MDD = toupper(NEA),  
    BETA_MDD = BETA,
    P_MDD = PVAL,
    CHR_MDD = CHROM,
    POS_MDD = POS
  )

mdd_EHR <- MDD_EHR_GWAS %>%
  transmute(
    SNP = ID,                
    EA_MDD = toupper(EA),   
    NEA_MDD = toupper(NEA),  
    BETA_MDD = BETA,
    P_MDD = PVAL,
    CHR_MDD = CHROM,
    POS_MDD = POS
  )

mdd_quest <- MDD_quest_GWAS %>%
  transmute(
    SNP = ID,                
    EA_MDD = toupper(EA),   
    NEA_MDD = toupper(NEA),  
    BETA_MDD = BETA,
    P_MDD = PVAL,
    CHR_MDD = CHROM,
    POS_MDD = POS
  )

mcp <- MCP_GWAS %>%
  transmute(
    SNP,
    A1_MCP = toupper(ALLELE1), 
    A2_MCP = toupper(ALLELE0),  
    BETA_MCP = BETA,
    P_MCP = P,
    CHR_MCP = CHR,
    POS_MCP = BP
  )

merged_MCP_MDD <- inner_join(mcp, mdd, by = "SNP")
merged_MCP_MDD_EHR <- inner_join(mcp, mdd_EHR, by = "SNP")
merged_MCP_MDD_quest <- inner_join(mcp, mdd_quest, by = "SNP")

comp <- function(a) chartr("ATCG", "TAGC", a)
is_pal <- function(a, b) (a=="A"&b=="T")|(a=="T"&b=="A")|(a=="C"&b=="G")|(a=="G"&b=="C")

merged_MCP_MDD <- merged_MCP_MDD %>% filter(!is_pal(A1_MCP, A2_MCP))
merged_MCP_MDD_EHR <- merged_MCP_MDD_EHR %>% filter(!is_pal(A1_MCP, A2_MCP))
merged_MCP_MDD_quest <- merged_MCP_MDD_quest %>% filter(!is_pal(A1_MCP, A2_MCP))

harmonise_datasets <- function(merged_df) {
  
  harmonised_df <- merged_df %>%
    mutate(
      cA1 = comp(A1_MCP),
      cA2 = comp(A2_MCP),
      direct_match = (A1_MCP == EA_MDD & A2_MCP == NEA_MDD),
      direct_swap  = (A1_MCP == NEA_MDD & A2_MCP == EA_MDD),
      strand_match = (cA1 == EA_MDD & cA2 == NEA_MDD),
      strand_swap  = (cA1 == NEA_MDD & cA2 == EA_MDD),
      keep = direct_match | direct_swap | strand_match | strand_swap,
      flip = direct_swap | strand_swap
    ) %>%
    filter(keep) %>%
    mutate(
      BETA_MCP_HARM = ifelse(flip, -BETA_MCP, BETA_MCP),
      same_dir = sign(BETA_MCP_HARM) == sign(BETA_MDD)
    )
  return(harmonised_df)
}

harmonised_MCP_MDD <- harmonise_datasets(merged_MCP_MDD)
harmonised_MCP_MDD_EHR <- harmonise_datasets(merged_MCP_MDD_EHR)
harmonised_MCP_MDD_quest <- harmonise_datasets(merged_MCP_MDD_quest)

harmonised_MCP_MDD <- harmonised_MCP_MDD %>%
  filter(!(CHR_MDD == 6 & between(POS_MDD, 25e6, 35e6))) %>%
  filter(!(CHR_MDD == 8 & between(POS_MDD, 7.2e6, 12.5e6)))

harmonised_MCP_MDD_EHR <- harmonised_MCP_MDD_EHR %>%
  filter(!(CHR_MDD == 6 & between(POS_MDD, 25e6, 35e6))) %>%
  filter(!(CHR_MDD == 8 & between(POS_MDD, 7.2e6, 12.5e6)))

harmonised_MCP_MDD_quest <- harmonised_MCP_MDD_quest %>%
  filter(!(CHR_MDD == 6 & between(POS_MDD, 25e6, 35e6))) %>%
  filter(!(CHR_MDD == 8 & between(POS_MDD, 7.2e6, 12.5e6)))

# conjFDR ---- 
# download.file("https://zenodo.org/record/5750318/files/genref.zip",
#              destfile = "./resources/conjFDR/genref.zip", mode = "wb")
# unzip("./resources/conjFDR/genref.zip", exdir = "./resources/conjFDR/genref")
REF_DIR <- "./resources/conjFDR/genref"

## MCP x MDD (EUR, UKB 23andMe excluded) ----
MCP <- data.frame(
  SNP = harmonised_MCP_MDD$SNP,
  BETA = harmonised_MCP_MDD$BETA_MCP_HARM,
  PVAL   = harmonised_MCP_MDD$P_MCP
)

MDD <- data.frame(
  SNP = harmonised_MCP_MDD$SNP,
  BETA = harmonised_MCP_MDD$BETA_MDD,
  PVAL   = harmonised_MCP_MDD$P_MDD
)

setDT(MCP); setDT(MDD)
MDD_MCP <- cfdr_pleio$new()

MDD_MCP$init_data(trait1 = MDD,
                  trait2 = MCP, 
                  trait_names = c("MDD", "MCP"),
                  refdat = refdata_location(REF_DIR), 
                  local_refdat_path = REF_DIR,
                  verbose = T)

MDD_MCP$initialize_pruning_index(n_iter = 50, seed = 151097, verbose = T, force=TRUE)

MDD_MCP$calculate_cond_fdr(fdr_trait = 1, verbose = T)
MDD_MCP$calculate_cond_fdr(fdr_trait = 2, verbose = T)

MDD_MCP_res <- MDD_MCP$get_trait_results()

## Save results
saveRDS(MDD_MCP, file = paste("./output/conjFDR/MCP_MDD_AnalysisObject.rds", sep="/"))
write.table(MDD_MCP_res, "./output/conjFDR/MCP_MDD_res.txt", row.names = F, quote = F)

## Format for sumstats.py
bim_ref <- fread("./resources/python_convert/9279485_ref.bim")
sumstats <- MDD_MCP_res %>%
  inner_join(bim_ref %>% select(SNP, A1, A2, CHR, BP), by = "SNP") %>%
  transmute(SNP, CHR = as.integer(CHR), BP = as.integer(BP),
            A1, A2,
            FDR = conj_fdr) 
fwrite(sumstats, "./resources/python_convert/conjFDR_MCP_MDD_for_clump.csv", sep = "\t", quote = FALSE)

## MCP x MDD (EHR, EUR) ----
MCP <- data.frame(
  SNP = harmonised_MCP_MDD_EHR$SNP,
  BETA = harmonised_MCP_MDD_EHR$BETA_MCP_HARM,
  PVAL   = harmonised_MCP_MDD_EHR$P_MCP
)

MDD <- data.frame(
  SNP = harmonised_MCP_MDD_EHR$SNP,
  BETA = harmonised_MCP_MDD_EHR$BETA_MDD,
  PVAL   = harmonised_MCP_MDD_EHR$P_MDD
)

setDT(MCP); setDT(MDD)
MDD_MCP <- cfdr_pleio$new()

MDD_MCP$init_data(trait1 = MDD,
                  trait2 = MCP, 
                  trait_names = c("MDD", "MCP"),
                  refdat = refdata_location(REF_DIR), 
                  local_refdat_path = REF_DIR,
                  verbose = T)

MDD_MCP$initialize_pruning_index(n_iter = 50, seed = 151097, verbose = T, force=TRUE)

MDD_MCP$calculate_cond_fdr(fdr_trait = 1, verbose = T)
MDD_MCP$calculate_cond_fdr(fdr_trait = 2, verbose = T)

MDD_MCP_res <- MDD_MCP$get_trait_results()

## Save results
saveRDS(MDD_MCP, file = paste("./output/conjFDR/MCP_MDD_EHR_AnalysisObject.rds", sep="/"))
write.table(MDD_MCP_res, "./output/conjFDR/MCP_MDD_EHR_res.txt", row.names = F, quote = F)

## Format for sumstats.py
sumstats <- MDD_MCP_res %>%
  inner_join(bim_ref %>% select(SNP, A1, A2, CHR, BP), by = "SNP") %>%
  transmute(SNP, CHR = as.integer(CHR), BP = as.integer(BP),
            A1, A2,
            FDR = conj_fdr) 
fwrite(sumstats, "./resources/python_convert/conjFDR_MCP_MDD_EHR_for_clump.csv", sep = "\t", quote = FALSE)

 