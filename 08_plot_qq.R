library(ggplot2)
library(data.table)
library(patchwork)
# Load GWAS and conjFDR results
MDD_MCP <- readRDS("./output/conjFDR/MDD_MCP_AnalysisObject.rds")
MCP_GWAS <- read.table("./resources/summary_statistics/MCP/ChronicPain.stats", header = T)

MDD_MCP_res <- MDD_MCP$get_trait_results()
conjFDR.clump.lead <- read.table("./output/conjFDR/conjFDR.clump.lead.csv", header = T)
file <- "./resources/summary_statistics/MDD//pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11.tsv"
MDD_GWAS <- fread(file, sep = "\t")

MDD_MCP_res <- left_join(MDD_MCP_res, select(MDD_GWAS, ID, CHROM, POS), by = c("SNP" = "ID"))

# Manhattan plots ----
df <- MDD_MCP_res %>%
  mutate(
    CHR = suppressWarnings(as.integer(as.character(CHROM))),
    BP  = as.numeric(POS),
    log10_conj = -log10(conj_fdr)
  )

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
    is_lead   = SNP %in% conjFDR.clump.lead$LEAD_SNP
  )

cons <- MDD_MCP_res %>%
  transmute(SNP, consistent = sign(BETA1) == sign(BETA2))

lead_df <- df %>%
  filter(is_lead) %>%
  left_join(cons, by = "SNP") %>%
  mutate(highlight_opposite = ifelse(consistent, "same", "opposite"))

axis_df <- chr_info

manhattan_p <- ggplot(df, aes(cumBP, log10_conj)) +
  geom_point(aes(color = chrom_col), size = 0.5, alpha = 0.6, show.legend = FALSE) +
  geom_point(data = lead_df,
             aes(x = cumBP, y = log10_conj, alpha = 0.4, color = chrom_col, shape = highlight_opposite, fill = highlight_opposite),
             size = 3, inherit.aes = FALSE, show.legend = FALSE) +
  scale_shape_manual(values = c(same = 18, opposite = 23)) +
  scale_fill_manual(values = c(same = "white", opposite = "")) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR, expand = c(0.01, 0)) +
  labs(x = "Chromosome", y = expression(-log[10]~"(conjFDR)"),
       title = "Manhattan plot of conjunctional -log10(FDR) for MDD and MCP") +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

ggsave("./output/conjFDR/manhattan.png", plot = manhattan_p, width = 11, height = 5, dpi = 600)

# Q-Q plots ----
MDD <- data.frame(
  SNP = MDD_GWAS$ID,
  BETA = MDD_GWAS$BETA,
  PVAL   = MDD_GWAS$PVAL
)

MCP <- data.frame(
  SNP = MCP_GWAS$SNP,
  BETA = MCP_GWAS$BETA,
  PVAL   = MCP_GWAS$P
)

dat <- inner_join(MDD %>% select(SNP, p_MDD = PVAL),
                  MCP %>% select(SNP, p_MCP = PVAL),
                  by = "SNP")

condqq_df <- function(p_main, p_cond,
                      thresholds = c(Inf, 0.1, 0.01, 0.001),
                      level_names = c("All SNPs",
                                      "p[MCP] <= 0.1",
                                      "p[MCP] <= 0.01",
                                      "p[MCP] <= 0.001")) {
  
  stopifnot(length(thresholds) == length(level_names))
  
  # Clean p-values: drop NA, cap to (0,1], replace zeros with tiny value
  p_main <- as.numeric(p_main)
  p_cond <- as.numeric(p_cond)
  p_main[p_main <= 0 | is.na(p_main)] <- NA_real_
  p_cond[is.na(p_cond)] <- NA_real_
  
  out <- vector("list", length(thresholds))
  
  for (i in seq_along(thresholds)) {
    thr  <- thresholds[i]
    keep <- if (is.infinite(thr)) rep(TRUE, length(p_main)) else p_cond <= thr
    keep <- keep & !is.na(p_main) & !is.na(p_cond)
    
    n <- sum(keep)
    if (n < 2) next
    
    p_sorted <- sort(p_main[keep])
    # guard against exact zeros (just in case)
    p_sorted[p_sorted == 0] <- .Machine$double.xmin
    
    obs <- -log10(p_sorted)
    exp <- -log10((seq_len(n)) / (n + 1))
    
    out[[i]] <- data.frame(
      expected = exp,
      observed = obs,
      stratum  = level_names[i],
      stringsAsFactors = FALSE
    )
  }
  
  d <- do.call(rbind, out)
  # Make stratum a factor so breaks work predictably
  d$stratum <- factor(d$stratum, levels = level_names)
  d
}

# Left plot: MDD conditioned on MCP
df_left <- condqq_df(
  p_main      = dat$p_MDD,
  p_cond      = dat$p_MCP,
  thresholds  = c(Inf, 0.1, 0.01, 0.001),
  level_names = c("All SNPs",
                  "p[MCP] <= 0.1",
                  "p[MCP] <= 0.01",
                  "p[MCP] <= 0.001")
)

# Right plot: MCP conditioned on MDD
df_right <- condqq_df(
  p_main      = dat$p_MCP,
  p_cond      = dat$p_MDD,
  thresholds  = c(Inf, 0.1, 0.01, 0.001),
  level_names = c("All SNPs",
                  "p[MDD] <= 0.1",
                  "p[MDD] <= 0.01",
                  "p[MDD] <= 0.001")
)

# Expression labels for pretty legends
lab_left  <- c("All SNPs",
               expression(p[MCP] <= 0.1),
               expression(p[MCP] <= 0.01),
               expression(p[MCP] <= 0.001))

lab_right <- c("All SNPs",
               expression(p[MDD] <= 0.1),
               expression(p[MDD] <= 0.01),
               expression(p[MDD] <= 0.001))

cols_left <- c(
  "All SNPs"        = "#1f77b4",
  "p[MCP] <= 0.1"   = "#ff7f0e",
  "p[MCP] <= 0.01"  = "#2ca02c",
  "p[MCP] <= 0.001" = "#d62728"
)

cols_right <- c(
  "All SNPs"        = "#1f77b4",
  "p[MDD] <= 0.1"   = "#ff7f0e",
  "p[MDD] <= 0.01"  = "#2ca02c",
  "p[MDD] <= 0.001" = "#d62728"
)

df_left$stratum  <- factor(as.character(df_left$stratum),  levels = names(cols_left))
df_right$stratum <- factor(as.character(df_right$stratum), levels = names(cols_right))

p1 <- ggplot(df_left, aes(expected, observed, colour = stratum)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = cols_left, name = NULL,
                     breaks = names(cols_left), labels = lab_left) +
  labs(title = "MDD | MCP",
       x = expression("Expected -log"[10]*"(p"[MDD]*")"),
       y = expression("Observed -log"[10]*"(p"[MDD]*")")) +
  coord_cartesian(xlim = c(0, max(df_left$expected)), ylim = c(0, max(df_left$observed))) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.72, 0.25),
    legend.direction = "vertical",
    legend.background = element_rect(fill = scales::alpha("white", 0.8), colour = NA),
    legend.key = element_blank(),
    legend.text = element_text(size = 10)
  ) + 
  coord_cartesian(ylim = c(0, 7), xlim = c(0, 7))

p2 <- ggplot(df_right, aes(expected, observed, colour = stratum)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = cols_right, name = NULL,
                     breaks = names(cols_right), labels = lab_right) +
  labs(title = "MCP | MDD",
       x = expression("Expected -log"[10]*"(p"[MCP]*")"),
       y = expression("Observed -log"[10]*"(p"[MCP]*")")) +
  coord_cartesian(xlim = c(0, max(df_right$expected)), ylim = c(0, max(df_right$observed))) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.72, 0.25),
    legend.direction = "vertical",
    legend.background = element_rect(fill = scales::alpha("white", 0.8), colour = NA),
    legend.key = element_blank(),
    legend.text = element_text(size = 10)
  ) + 
  coord_cartesian(ylim = c(0, 7), xlim = c(0, 7))

qq_p <- p1 + p2

ggsave("./output/conjFDR/qq.png", plot = qq_p, width = 11, height = 6, dpi = 600)
