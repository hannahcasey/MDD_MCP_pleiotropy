library(ggvenn)
library(dplyr)
library(stringr)
library(scales)
library(ggplot2)
library(patchwork)
library(data.table)
library(tidyr)

## Load FUMA output ----
FUMA_annov_stats <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/FUMA_job663071//annov.stats.txt", header = T)
FUMA_annot <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/FUMA_job663071/annot.txt", header = T)
FUMA_snps <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/FUMA_job663071/snps.txt", header = T)
FUMA_genes<- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/FUMA_job663071/genes.txt", header = T)

# Load GWAS and conjFDR results
MDD_MCP <- readRDS("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/MDD_MCP_AnalysisObject.rds")
MCP_GWAS <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MCP/ChronicPain.stats", header = T)

MDD_MCP_res <- MDD_MCP$get_trait_results()
conjFDR.clump.lead <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.lead.csv", header = T)
file <- "/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/resources/summary_statistics/MDD//pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv"
MDD_GWAS <- fread(file, sep = "\t")
MDD_MCP_res <- left_join(MDD_MCP_res, dplyr::select(MDD_GWAS, ID, CHROM, POS), by = c("SNP" = "ID"))

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
  mutate(
    cumBP     = BP + offset,
    chrom_col = factor(CHR %% 2),
    is_lead   = SNP %in% conjFDR.clump.lead$LEAD_SNP
  )

cons <- MDD_MCP_res %>%
  transmute(SNP, consistent = sign(BETA1) == sign(BETA2))

lead_df_same <- df %>% filter(is_lead, sign(BETA1) == sign(BETA2))
lead_df_opp  <- df %>% filter(is_lead, sign(BETA1) != sign(BETA2))

axis_df <- chr_info

df_plot <- df
lead_df_same_plot <- lead_df_same

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
    fill = "grey50",
  ) +
  geom_hline(yintercept = -log10(0.01), linetype = 2) +
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
  )

## Functional consequence ----
annot_levels <- c("ncRNA_splicing","splicing","UTR5","exonic","ncRNA_exonic",
                  "UTR3","downstream","upstream","ncRNA_intronic","intergenic","intronic")

df <- FUMA_annov_stats %>%
  mutate(
    annot = factor(annot, levels = annot_levels),
  ) %>%
  filter(!is.na(annot), !is.na(prop)) %>%
  arrange(annot)

m <- nlevels(df$annot)

functional_consequence_p <- ggplot(df, aes(x = prop, y = annot)) +
  geom_col(width = 0.85, fill = "#F8766D", linewidth = 0.3) +
  scale_y_discrete(name = "", expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(name = "Proportion",
                     expand = expansion(mult = c(0, 0))) +
  labs(x = NULL) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x  = element_text(angle = 60, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(colour = "grey40"),
    legend.position    = "inside",
    legend.position.inside = c(0.9, 0.5)
  ) +
  coord_cartesian(clip = "off")

## Minimum chromatin state ----
df <- tibble(minChrState = FUMA_snps$minChrState) %>%
  count(minChrState, name = "n") %>%
  right_join(tibble(minChrState = 1:15), by = "minChrState") %>%
  mutate(n = replace_na(n, 0L)) %>%
  arrange(minChrState) %>%
  mutate(state_f = factor(minChrState, levels = rev(1:15)))  # 15 at top like the figure

min_chromatic_states_p <- ggplot(df, aes(x = n, y = state_f)) +
  geom_col(fill = "#00BFC4") +
  labs(
    x = "Number of SNP",
    y = "Minimum Chromatin State"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(margin = margin(t = 6))
  )

## RegulomeDB Score ----
rdb_vec <- FUMA_snps$RDB

rdb_levels <- c("7","6","5","4","3a","2c","2b","2a","1f","1d","1c","1b","NA")

df <- tibble(RDB = rdb_vec) %>%
  mutate(RDB = ifelse(is.na(RDB), "NA", RDB)) %>%
  count(RDB, name = "n") %>%
  right_join(tibble(RDB = rdb_levels), by = "RDB") %>%   # include zeros
  mutate(n = replace_na(n, 0L),
         RDB_f = factor(RDB, levels = rdb_levels))       # order as panel C

RDB_p <- ggplot(df, aes(x = n, y = RDB_f)) +
  geom_col(fill = "#F8766D") +
  labs(x = "Number of SNP", y = "RegulomeDB Score") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(margin = margin(t = 6))
  )

## CADD ----
cadd_med <- median(FUMA_snps$CADD)
thr_top1 <- 12.37   # canonical "top 1%" CADD cut

CADD_p <- ggplot(FUMA_snps, aes(x = CADD)) +
  geom_density(fill = "#00BFC4", colour = "white",  alpha = 1, adjust = 1) +
  geom_vline(xintercept = cadd_med, linetype = "solid",  color = "grey10", linewidth = 0.6) +
  geom_vline(xintercept = thr_top1, linetype = "solid",  color = "grey10", linewidth = 0.6) +
  labs(x = "CADD Score", y = "Density") +
  annotate("text", x = cadd_med, y = 0, label = sprintf("median = %.2f", cadd_med),
           vjust = -1.2, hjust = 0, size = 3.3, color = "grey10") +
  annotate("text", x = thr_top1, y = 0, label = "CADD≈12.37", vjust = -3,
           hjust = -0.1, size = 3.3, color = "grey30") +
  theme_classic(base_size = 12)

## Gene-mapping ----
gene_id <- FUMA_genes$symbol
genes_pos <- gene_id[FUMA_genes$posMapSNPs > 0] %>% unique()
genes_eqtl <- gene_id[FUMA_genes$eqtlMapSNPs > 0] %>% unique()
genes_ci <- gene_id[FUMA_genes$ciMap == "Yes"] %>% unique()

gene_sets <- list(
  Positional = genes_pos,
  `eQTL (GTEx)` = genes_eqtl,
  `Chromatin Interaction `  = genes_ci
)

sapply(gene_sets, length)

cols <- c("#F8766D", "#00BFC4", "grey")

gene_p <- ggvenn(
  gene_sets,
  fill_color   = alpha(cols, 0.5),
  stroke_color = "white",
  set_name_size = 4,
  text_size     = 4,
  show_percentage = FALSE) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))


small_theme_up <- theme(plot.margin = margin(t = 1, unit = "cm"),
                     legend.position = "none")

small_theme_up_w_legend <- theme(plot.margin = margin(t = 1, unit = "cm"))

small_theme <- theme(plot.margin = margin(t = 0, unit = "cm"),
                     legend.position = "none")

small_row_centered <-
  ((functional_consequence_p + small_theme_up_w_legend) |
   (min_chromatic_states_p + small_theme_up)   |
   (RDB_p + small_theme_up)) /
  (plot_spacer() | (CADD_p + small_theme) | plot_spacer() | (gene_p + small_theme) | plot_spacer()) +
  plot_layout(widths = c(1,1,1,1), heights = c(1,1.5))

FUMA_plot <-
  manhattan_p /
  small_row_centered +
  plot_layout(heights = c(1, 2), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"))

ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/results.png", plot = FUMA_plot, width = 15, height = 13, dpi = 600)






lead_full <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.lead.csv", header = T)
locus_full <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.loci.csv", header = T)


lead <- read.table("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/conjFDR/conjFDR.clump.lead.fuma.txt", header = T)



table(lead_full$is_locus_lead)
max(locus_full$FDR)

## Gene-mapping ----
gene_id <- FUMA_genes$symbol
genes_pos <- gene_id[FUMA_genes$posMapSNPs > 0] %>% unique()
genes_eqtl <- gene_id[FUMA_genes$eqtlMapSNPs > 0] %>% unique()
genes_ci <- gene_id[FUMA_genes$ciMap == "Yes"] %>% unique()

gene_sets <- list(
  Positional = genes_pos,
  `eQTL (GTEx)` = genes_eqtl,
  `Chromatin Interaction `  = genes_ci
)

sapply(gene_sets, length)

cols <- scales::hue_pal()(length(gene_sets))

gene_p <- ggvenn(
  gene_sets,
  fill_color   = alpha(cols, 0.5),
  stroke_color = "grey35",
  set_name_size = 4,
  text_size     = 4,
  show_percentage = FALSE) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

