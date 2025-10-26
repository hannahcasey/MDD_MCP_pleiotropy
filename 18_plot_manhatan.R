library(ggplot2)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
library(ggrepel)

# Load conjFDR results
MDD_MCP <- readRDS("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_AnalysisObject.rds")
MDD_MCP_res <- MDD_MCP$get_trait_results()
conjFDR.clump.lead <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.lead.csv", header = T)
conjFDR.clump.loci <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.loci.csv", header = T)

file <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MDD//pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv"
MDD_GWAS <- fread(file, sep = "\t")
MDD_MCP_res <- left_join(MDD_MCP_res, dplyr::select(MDD_GWAS, ID, CHROM, POS), by = c("SNP" = "ID"))

# Load genes
genes <- read.csv("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/V2G/loci_annotated_in_text.csv", header = T)
gene_list <- genes %>%
  mutate(v2g_list = strsplit(v2g_symbols, ";", fixed = TRUE)) %>%
  pull(v2g_list) %>%
  unlist(use.names = FALSE) %>%
  unique()

# Manhattan plots ----
df <- MDD_MCP_res %>%
  mutate(
    CHR        = suppressWarnings(as.integer(as.character(CHROM))), # force numeric chr
    BP         = suppressWarnings(as.numeric(POS)),
    log10_conj = -log10(conj_fdr)
  ) %>%
  arrange(CHR, BP)

chr_info <- df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(
    chr_len = as.numeric(chr_len),
    offset  = dplyr::lag(cumsum(chr_len), default = 0),
    center  = offset + chr_len / 2
  )

df <- df %>%
  left_join(dplyr::select(chr_info, CHR, offset), by = "CHR") %>%
  dplyr::mutate(
    cumBP     = BP + offset,
    chrom_col = factor(CHR %% 2),
    is_lead   = SNP %in% conjFDR.clump.lead$LEAD_SNP,
    is_colocalised = SNP %in% genes$RSID
  )

cons <- MDD_MCP_res %>%
  transmute(SNP, consistent = sign(BETA1) == sign(BETA2))

lead_df_same <- df %>% filter(is_lead, sign(BETA1) == sign(BETA2))
lead_df_opp  <- df %>% filter(is_lead, sign(BETA1) != sign(BETA2))

lead_df_colocalised <- df %>% filter(is_colocalised)

axis_df <- chr_info

df_plot <- df
lead_df_same_plot <- lead_df_same

labels_df <- df %>%
  filter(is_colocalised) %>%              
  left_join(genes %>% dplyr::select(LEAD_BP, Symbol),
            by = c("SNP" = "RSID")) %>%
  separate_rows(Symbol, sep = ";") %>%
  mutate(label_gene = Symbol) %>%
  group_by(SNP) %>%
  summarise(
    cumBP = dplyr::first(cumBP),
    log10_conj = dplyr::first(log10_conj),
    label = paste(unique(label_gene), collapse = ", "),
    .groups = "drop"
  )

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
