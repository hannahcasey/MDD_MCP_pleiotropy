library(data.table)
# remotes::install_github("alexploner/cfdr.pleio", build_vignettes = TRUE)
library(cfdr.pleio)
library(dplyr)
library(tidyr)

# MDD GWAS ----
file <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MDD//pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv"
MDD_GWAS <- fread(file, sep = "\t")

# MCP GWAS ----
MCP_GWAS <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MCP/ChronicPain.stats", header = T)

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

merged <- inner_join(mcp, mdd, by = "SNP")

comp <- function(a) chartr("ATCG", "TAGC", a)
is_pal <- function(a, b) (a=="A"&b=="T")|(a=="T"&b=="A")|(a=="C"&b=="G")|(a=="G"&b=="C")

merged <- merged %>% filter(!is_pal(A1_MCP, A2_MCP))

harmonised <- merged %>%
  mutate(
    cA1 =  comp(A1_MCP), cA2 = comp(A2_MCP),
    direct_match = (A1_MCP==EA_MDD & A2_MCP==NEA_MDD),
    direct_swap  = (A1_MCP==NEA_MDD & A2_MCP==EA_MDD),
    strand_match = (cA1==EA_MDD & cA2  ==NEA_MDD),
    strand_swap  = (cA1==NEA_MDD & cA2  ==EA_MDD),
    keep = direct_match | direct_swap | strand_match | strand_swap,
    flip = direct_swap  | strand_swap
  ) %>%
  filter(keep) %>%
  mutate(
    BETA_MCP_HARM = ifelse(flip, -BETA_MCP, BETA_MCP),
    same_dir = sign(BETA_MCP_HARM) == sign(BETA_MDD)  # TRUE = same direction
  )

setDT(harmonised)
harmonised <- harmonised[!(CHR_MDD == 6 & POS_MDD >= 25e6 & POS_MDD <= 35e6)]
harmonised <- harmonised[!(CHR_MDD == 8 & POS_MDD >= 7.2e6 & POS_MDD <= 12.5e6)]

MCP <- data.frame(
  SNP = harmonised$SNP,
  BETA = harmonised$BETA_MCP_HARM,
  PVAL   = harmonised$P_MCP
)

MDD <- data.frame(
  SNP = harmonised$SNP,
  BETA = harmonised$BETA_MDD,
  PVAL   = harmonised$P_MDD
)

setDT(MCP); setDT(MDD)

# conjFDR ---- 
# download.file("https://zenodo.org/record/5750318/files/genref.zip",
#              destfile = "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/genref.zip", mode = "wb")
#unzip("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/genref.zip", exdir = "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/genref")
REF_DIR <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/genref"

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

saveRDS(MDD_MCP, file = paste("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_AnalysisObject.rds", sep="/"))
write.table(MDD_MCP_res, "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_res.txt", row.names = F, quote = F)

## Format for sumstats.py
ref <- fread("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/python_convert/9279485_ref.bim")

sumstats <- MDD_MCP_res %>%
  inner_join(ref %>% select(SNP, A1, A2), by = "SNP") %>%
  transmute(SNP, CHR = as.integer(CHR), BP = as.integer(BP),
            A1, A2,
            FDR = conj_fdr) 

fwrite(sumstats, "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/python_convert/conjFDR_for_clump.csv", sep = "\t", quote = FALSE)


