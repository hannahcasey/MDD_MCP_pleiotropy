library(data.table)
# remotes::install_github("alexploner/cfdr.pleio", build_vignettes = TRUE)
library(cfdr.pleio)
library(dplyr)
library(bigsnpr)
library(ggrepel)
library(ggplot2)

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
#download.file("https://zenodo.org/record/5750318/files/genref.zip",
#             destfile = "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/genref.zip", mode = "wb")
#unzip("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/genref.zip", exdir = "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/genref")
REF_DIR <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/genref"

MDD_MCP <- cfdr_pleio$new()

MDD_MCP$init_data(trait1 = MDD,
                  trait2 = MCP, 
                  trait_names = c("MDD", "MCP"),
                  refdat = refdata_location(REF_DIR), 
                  local_refdat_path = "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/conjFDR/",
                  verbose = T)

MDD_MCP$initialize_pruning_index(n_iter = 50, seed = 151097, verbose = T)

MDD_MCP$calculate_cond_fdr(fdr_trait = 1, verbose = FALSE)
MDD_MCP$calculate_cond_fdr(fdr_trait = 2, verbose = FALSE)

MDD_MCP_res <- MDD_MCP$get_trait_results()

saveRDS(MDD_MCP, file = paste("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_AnalysisObject.rds", sep="/"))
write.table(MDD_MCP_res, "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_res.txt", row.names = F, quote = F)

# Lead SNPs ----
# MDD_MCP <- readRDS("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_AnalysisObject.rds")
# MDD_MCP_res <- MDD_MCP$get_trait_results()

MDD_MCP_res <- left_join(MDD_MCP_res, select(MCP_GWAS, SNP, CHR, BP), by = "SNP")

ld_ref_path <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/1000g/EUR"
plink_path  <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/plink_mac_20250615/plink"
out_dir     <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/"

shared_sig <- MDD_MCP_res %>% filter(conj_fdr < 0.05)
f1 <- paste0(out_dir, "conj_sig.tsv")
fwrite(shared_sig[, .(SNP, P = conj_fdr, CHR, BP)], f1, sep = "\t", quote = FALSE, na = "NA")

system(paste(
  shQuote(plink_path), "--bfile", shQuote(ld_ref_path),
  "--clump", shQuote(f1),
  "--clump-p1 1 --clump-p2 1 --clump-r2 0.6 --clump-kb 1000000",
  "--threads 8 --memory 8000",
  "--out", paste0(out_dir, "clump_r2_0.6")
))

idx <- fread(paste0(out_dir, "clump_r2_0.6.clumped"))$SNP
f2 <- paste0(out_dir, "conj_sig_idx.tsv")
fwrite(shared_sig[SNP %in% idx, .(SNP, P = conj_fdr, CHR, BP)], f2, sep = "\t", quote = FALSE, na = "NA")

system(paste(
  shQuote(plink_path), "--bfile", shQuote(ld_ref_path),
  "--clump", shQuote(f2),
  "--clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 1000000",
  "--threads 8 --memory 8000",
  "--out", paste0(out_dir, "clump_r2_0.1")
))

lead <- fread(paste0(out_dir, "clump_r2_0.1.clumped"))[, .(SNP)]
lead <- lead %>% left_join(shared_sig %>% select(SNP, CHR, BP, conj_fdr), by = "SNP") %>% distinct()

lead_merged <- lead %>%
  arrange(CHR, BP) %>%
  group_by(CHR) %>%
  mutate(gap = BP - lag(BP, default = first(BP)),
         block = cumsum(ifelse(is.na(gap), 0L, gap > 250000))) %>%
  group_by(CHR, block) %>%
  slice_min(conj_fdr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(SNP, CHR, BP, conj_fdr)
 
MDD_MCP_res %>%
  filter(SNP %in% lead_merged$SNP) %>%
  mutate(consistent = sign(BETA1) == sign(BETA2)) %>%
  count(consistent)

#TRUE 389
#FALSE 18

fwrite(lead_merged, paste0(out_dir, "conjFDR_lead_SNPs.tsv"), sep = "\t", quote = FALSE)
fwrite(select(lead_merged, SNP, CHR, BP), paste0(out_dir, "conjFDR_lead_SNPs_fuma.tsv"), sep = "\t", quote = FALSE, col.names = F)

# Manhattan plots ----
df <- MDD_MCP_res %>%
  mutate(
    CHR = suppressWarnings(as.integer(as.character(CHR))),
    BP  = as.numeric(BP),
    log10_conj = -log10(conj_fdr)
  ) %>%
  filter(CHR %in% 1:22, is.finite(BP), is.finite(log10_conj))

chr_info <- df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(
    chr_len = as.numeric(chr_len),
    offset  = dplyr::lag(cumsum(chr_len), default = 0),
    center  = offset + chr_len/2
  )

df <- df %>%
  left_join(select(chr_info, CHR, offset), by = "CHR") %>%
  mutate(
    cumBP     = BP + offset,
    chrom_col = factor(CHR %% 2),
    is_lead   = SNP %in% lead_merged$SNP
  )

cons <- MDD_MCP_res %>%
  transmute(SNP, consistent = sign(BETA1) == sign(BETA2))

lead_df <- df %>%
  filter(is_lead) %>%
  left_join(cons, by = "SNP") %>%
  mutate(highlight_opposite = ifelse(consistent, "same", "opposite"))

axis_df <- chr_info

manhattan_p <- ggplot(df, aes(cumBP, log10_conj)) +
  geom_point(aes(color = chrom_col), size = 0.5, alpha = 0.7, show.legend = FALSE) +
  geom_point(data = lead_df,
             aes(x = cumBP, y = log10_conj, color = chrom_col, shape = highlight_opposite, fill = highlight_opposite),
             size = 3, inherit.aes = FALSE, show.legend = FALSE) +
  scale_shape_manual(values = c(same = 18, opposite = 23)) +
  scale_fill_manual(values = c(same = "white", opposite = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR, expand = c(0.01, 0)) +
  labs(x = "Chromosome", y = expression(-log[10]~"(conjFDR)"),
       title = "Manhattan plot of conjunctional -log10(FDR) for MDD and MCP") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/mahattan.png", plot = manhattan_p, width = 11, height = 5, dpi = 600)

# Q-Q plots ----
dat <- inner_join(MDD %>% select(SNP, p_MDD = PVAL),
                  MCP %>% select(SNP, p_MCP = PVAL),
                  by = "SNP")

condqq_df <- function(p_main, p_cond, thresholds = c(Inf, 0.1, 0.01, 0.001),
                      labels = c("All SNPs", expression(p[MCP] <= 0.1),
                                 expression(p[MCP] <= 0.01), expression(p[MCP] <= 0.001))) {
  stopifnot(length(thresholds) == length(labels))
  out <- vector("list", length(thresholds))
  for (i in seq_along(thresholds)) {
    thr <- thresholds[i]
    keep <- if (is.infinite(thr)) rep(TRUE, length(p_main)) else p_cond <= thr
    n <- sum(keep)
    if (n < 2) next
    p_sorted <- sort(p_main[keep])
    obs <- -log10(p_sorted)
    exp <- -log10((seq_len(n)) / (n + 1))  # expected order stats
    out[[i]] <- data.frame(
      expected = exp,
      observed = obs,
      stratum  = labels[i]
    )
  }
  do.call(rbind, out)
}

lab_left  <- c("All SNPs",
               expression(p[MCP] <= 0.1),
               expression(p[MCP] <= 0.01),
               expression(p[MCP] <= 0.001))

lab_right <- c("All SNPs",
               expression(p[MDD] <= 0.1),
               expression(p[MDD] <= 0.01),
               expression(p[MDD] <= 0.001))

cols <- c("All SNPs" = "#1f77b4", "p[MCP] \u2264 0.1" = "#ff7f0e",
          "p[MCP] \u2264 0.01" = "#2ca02c", "p[MCP] \u2264 0.001" = "#d62728")

p1 <- ggplot(df_left, aes(expected, observed, colour = stratum)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_line(size = 1) +
  scale_color_manual(values = cols, name = NULL,
                     breaks = names(cols), labels = lab_left) +
  labs(title = "MDD | MCP",
       x = expression("Expected  -log"[10]*"(p"[MDD]*")"),
       y = expression("Observed -log"[10]*"(p"[MDD]*")")) +
  coord_cartesian(xlim = c(0, max(df_left$expected)), ylim = c(0, max(df_left$observed))) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.72, 0.25),
    legend.direction  = "vertical",
    legend.background = element_rect(fill = scales::alpha("white", 0.8), colour = NA),
    legend.key        = element_blank(),
    legend.text       = element_text(size = 10)
  ) +  
  coord_cartesian(ylim = c(0, 7), xlim = c(0, 7))

cols2 <- c("All SNPs" = "#1f77b4", "p[MDD] \u2264 0.1" = "#ff7f0e",
           "p[MDD] \u2264 0.01" = "#2ca02c", "p[MDD] \u2264 0.001" = "#d62728")

p2 <- ggplot(df_right, aes(expected, observed, colour = stratum)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_line(size = 1) +
  scale_color_manual(values = cols2, name = NULL,
                     breaks = names(cols2), labels = lab_right) +
  labs(title = "MCP | MDD",
       x = expression("Expected  -log"[10]*"(p"[MCP]*")"),
       y = expression("Observed -log"[10]*"(p"[MCP]*")")) +
  coord_cartesian(xlim = c(0, max(df_right$expected)), ylim = c(0, max(df_right$observed))) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.72, 0.25),
    legend.direction  = "vertical",
    legend.background = element_rect(fill = scales::alpha("white", 0.8), colour = NA),
    legend.key        = element_blank(),
    legend.text       = element_text(size = 10)
  ) +  
  coord_cartesian(ylim = c(0, 7), xlim = c(0, 7))

qq_p <- p1 + p2

ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/qq.png", plot = qq_p, width = 11, height = 6, dpi = 600)


 