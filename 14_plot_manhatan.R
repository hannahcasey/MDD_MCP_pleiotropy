library(ggplot2)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(ggrepel)
library(data.table)
library(tidyr)

# Load conjFDR results
MDD_MCP <- readRDS("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_AnalysisObject.rds")
MDD_MCP_res <- MDD_MCP$get_trait_results()
conjFDR.clump.lead <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.lead.csv", header = T)
conjFDR.clump.loci <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.loci.csv", header = T)

file <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MDD//pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv"
MDD_GWAS <- fread(file, sep = "\t")
MDD_MCP_res <- left_join(MDD_MCP_res, dplyr::select(MDD_GWAS, ID, CHROM, POS), by = c("SNP" = "ID"))

# Load genes
genes <- read.csv("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/L2G/loci_in_text.csv", header = T, check.names = F)

genes <- left_join(genes, conjFDR.clump.loci, by = c("Min pos" = "MinBP",   "Max pos" = "MaxBP"))

names(genes)[is.na(names(genes)) | names(genes) == ""] <- paste0("unnamed_", seq_len(sum(is.na(names(genes)) | names(genes) == "")))

gene_list <- genes %>%
  dplyr::mutate(v2g_list = strsplit(Symbol, ";", fixed = TRUE)) %>%
  pull(v2g_list) %>%
  unlist(use.names = FALSE) %>%
  unique()

# Manhattan plots ----
manhattan_df <- MDD_MCP_res %>%
  mutate(
    CHR        = suppressWarnings(as.integer(as.character(CHROM))), # force numeric chr
    BP         = suppressWarnings(as.numeric(POS)),
    log10_conj = -log10(conj_fdr)
  ) %>%
  left_join(conjFDR.clump.loci, by = c("SNP" = "LEAD_SNP", "CHR")) %>%
  arrange(CHR, BP)

chr_info <- manhattan_df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(
    chr_len = as.numeric(chr_len),
    offset  = dplyr::lag(cumsum(chr_len), default = 0),
    center  = offset + chr_len / 2
  )

manhattan_df <- manhattan_df %>%
  left_join(dplyr::select(chr_info, CHR, offset), by = "CHR") %>%
  dplyr::mutate(
    cumBP     = BP + offset,
    chrom_col = factor(CHR %% 2),
    is_lead   = SNP %in% conjFDR.clump.lead$LEAD_SNP,
    is_colocalised = SNP %in% genes$LEAD_SNP
  )

cons <- MDD_MCP_res %>%
  transmute(SNP, consistent = sign(BETA1) == sign(BETA2))

lead_df_same <- manhattan_df %>% filter(is_lead, sign(BETA1) == sign(BETA2))
lead_df_opp  <- manhattan_df %>% filter(is_lead, sign(BETA1) != sign(BETA2))

lead_df_colocalised <- manhattan_df %>% filter(is_colocalised)

axis_df <- chr_info

df_plot <- manhattan_df
lead_df_same_plot <- lead_df_same

labels_df <- manhattan_df %>%
  filter(is_colocalised) %>%              
  left_join(genes %>% dplyr::select("Min pos", "Max pos", Symbol),
            by = c("MinBP" = "Min pos", "MaxBP" = "Max pos")) %>%
  separate_rows(Symbol, sep = ";") %>%
  mutate(label_gene = Symbol) %>%
  group_by(SNP) %>%
  summarise(
    cumBP = dplyr::first(cumBP),
    log10_conj = dplyr::first(log10_conj),
    label = paste(unique(label_gene), collapse = ", "),
    .groups = "drop"
  )

labels_df$label[labels_df$label == "NA"] <- ""

manhattan_p <- ggplot(df_plot, aes(x = cumBP, y = log10_conj)) +
  geom_point(aes(color = chrom_col), size = 0.5, alpha = 0.6, show.legend = FALSE) +
  geom_point(
    data = lead_df_same_plot,
    aes(x = cumBP, y = log10_conj, color = chrom_col),
    fill = "white", shape = 23, size = 3, inherit.aes = FALSE, 
    show.legend = FALSE,
  ) +
  geom_point(
    data = lead_df_opp,
    aes(x = cumBP, y = log10_conj, color = chrom_col),
    shape = 23, size = 3, inherit.aes = FALSE, show.legend = FALSE,
    fill = "grey70",
  ) +
  geom_point(
    data = lead_df_colocalised,
    aes(x = cumBP, y = log10_conj, color = chrom_col),
    shape = 23, size = 3, inherit.aes = FALSE, show.legend = FALSE,
    fill = "red",
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR,
    expand = c(0.01, 0)
  ) +
  labs(
    x = "Chromosome",
    y = expression(-log[10]~"(conjFDR)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  geom_label_repel(
    data = labels_df,
    aes(x = cumBP, y = log10_conj, label = label),
    size = 3,
    max.overlaps = Inf,      # or set a number (e.g., 50) to limit
    box.padding = 0.35,
    point.padding = 0.25,
    segment.size = 0.2,
    min.segment.length = 0
  )


ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/manhattan.png", plot = manhattan_p, width = 12, height = 6, dpi = 600)
