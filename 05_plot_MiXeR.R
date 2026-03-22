library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)
library(patchwork)
#install.packages("eulerr")
library(eulerr)
library(scales)
library(cowplot)

# Load MiXeR output ----
MDD_and_MCP <- read.csv("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/MiXeR/MDD_and_MCP.test.csv", sep = "\t", header = T)
MDD_vs_MCP <- read.csv("/Users/hannahcasey/Desktop/PhD/projects/MCP_MDD_colocalisation/output/MiXeR/MDD_vs_MCP.csv", sep = "\t", header = T)

# Univartiate plots
MDD_and_MCP_plot <- MDD_and_MCP %>%
  mutate(
    trait = c("Major Depressive Disorder", "Multisite Chronic Pain"),
    poly = `nc.p9..mean.`,
    poly_low = pmax(0, `nc.p9..mean.` - `nc.p9..std.`),
    poly_high =       `nc.p9..mean.` + `nc.p9..std.`,
    disc = `sig2_beta..mean.`,
    disc_low = pmax(1e-12, `sig2_beta..mean.` - `sig2_beta..std.`),
    disc_high =           `sig2_beta..mean.` + `sig2_beta..std.`,
    h2 = `h2..mean.`,
    h2_low = pmax(0, `h2..mean.` - `h2..std.`),
    h2_high =            `h2..mean.` + `h2..std.`
  ) %>%
  mutate(trait = fct_reorder(trait, h2))

MDD_and_MCP_plot$trait

dot_panel <- function(data, x, xmin, xmax, xlab, log10_x = FALSE) {
  p <- ggplot(data, aes(x = !!sym(x), y = trait, color = trait)) +
    geom_point(size = 2.4, show.legend = FALSE) +
    geom_errorbarh(aes(xmin = !!sym(xmin), xmax = !!sym(xmax)), height = 0.15) +
    labs(x = xlab, y = NULL) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      # >>> ADD BORDER AROUND EACH PANEL <<<
      panel.border = element_rect(color = "darkgrey", fill = NA, linewidth = 0.6)
    )
  # axis scaling
  if (log10_x) {
    p <- p + scale_x_continuous(trans = "log10", limits = c(max(min(data[[xmin]]), 1e-12), NA))
  } else {
    p <- p + scale_x_continuous(limits = c(0, NA))
  }
  p
}

p_h2   <- dot_panel(MDD_and_MCP_plot, "h2",  "h2_low",  "h2_high",  "Heritability (h²)")
p_poly <- dot_panel(MDD_and_MCP_plot, "poly","poly_low","poly_high","Polygenicity (nc@p9)")
p_disc <- dot_panel(MDD_and_MCP_plot, "disc","disc_low","disc_high",
                    expression(Discoverability~(sigma[beta]^2)), log10_x = F)

p_univariate <- p_h2 / p_poly / p_disc

# Bivariate
row <- MDD_vs_MCP %>%
  dplyr::filter(grepl("\\.fit\\.json$", fname)) %>%
  dplyr::slice(1)

fmtk <- function(x) sprintf("%.1fk", round(x/100, 1)/10)  # one decimal in k
MDD_n_var     <- paste0(fmtk(row$`nc1.p9..mean.`), "\n(", fmtk(row$`nc1.p9..std.`), ")")
MCP_n_var     <- paste0(fmtk(row$`nc2.p9..mean.`), "\n(", fmtk(row$`nc2.p9..std.`), ")")
overlap_n_var <- paste0(fmtk(row$`nc12.p9..mean.`), "\n(", fmtk(row$`nc12.p9..std.`), ")")

lab1 <- "MDD"
lab2 <- "MCP"

pi1  <- row$`pi1..mean.`          # MDD polygenicity (fraction)
pi2  <- row$`pi2..mean.`          # MCP polygenicity (fraction)
dice <- row$`dice..mean.`         # Dice overlap

shared_frac <- dice * (pi1 + pi2) / 2

M <- 1e6
A  <- pi1  * M         # area for MDD
B  <- pi2  * M         # area for MCP
AB <- shared_frac * M * 3  # overlap area

fit <- euler(c(
  "MDD"        = A,
  "MCP"        = B,
  "MDD&MCP"    = AB
))

gg_cols <- scales::hue_pal()(2)

venn_diagram <- plot(
  fit,
  fills      = list(fill = gg_cols, alpha = 0.6),
  edges      = list(col = "white", lwd = 1),
  labels     = list(labels = c(lab1, lab2), col = "black", font = 1, cex = 1.3),
  quantities = list(
    # supply custom text labels in region order: MDD, MCP, overlap
    labels = c(MDD_n_var, MCP_n_var, overlap_n_var),
    cex    = 1.1,
    font   = 1,
    col    = "black"
  )
)


df_rg <- tibble(
  mean = row$`rg..mean.`,
  sd   = row$`rg..std.`,
  xmin = mean - sd,
  xmax = mean + sd
)

p_rg <- ggplot(df_rg, aes(y = " ", x = mean)) +
  geom_col(width = 0.45, fill = "#92ea8b") +
  geom_errorbarh(aes(xmin = xmin, xmax = xmax), height = 0.18, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
  labs(x = expression("Genetic correlation ("*r[g]*")"), y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )


df_conc <- tibble(
  mean = row$`fraction_concordant_within_shared..mean.`,
  sd   = row$`fraction_concordant_within_shared..std.`,
  xmin = pmax(0, mean - sd),
  xmax = pmin(1, mean + sd)
)

p_conc <- ggplot(df_conc, aes(y = " ", x = mean)) +
  geom_col(width = 0.45, fill = "#92ea8b") +
  geom_errorbarh(aes(xmin = xmin, xmax = xmax), height = 0.18, linewidth = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dotted", linewidth = 0.4) +
  scale_x_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1, suffix = "")) +
  labs(x = "Concordance within shared variants (%)", y = NULL) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

A_panel <- patchwork::wrap_plots(p_univariate)
B_panel <- patchwork::wrap_plots(venn_diagram)
C_panel <- patchwork::wrap_plots(p_rg / p_conc)

A_tagged <- cowplot::ggdraw(A_panel) +
  cowplot::draw_label("A", x = 0.01, y = 0.99,
                      hjust = 0, vjust = 1,
                      fontface = "bold", size = 16)

B_tagged <- cowplot::ggdraw(B_panel) +
  cowplot::draw_label("B", x = 0.01, y = 0.99,
                      hjust = 0, vjust = 1,
                      fontface = "bold", size = 16)

C_tagged <- cowplot::ggdraw(C_panel) +
  cowplot::draw_label("C", x = 0.01, y = 0.99,
                      hjust = 0, vjust = 1,
                      fontface = "bold", size = 16)

pad <- function(p, t=6, r=6, b=6, l=6) p + theme(plot.margin = margin(t, r, b, l))
A_pad <- pad(A_tagged, t = 8, r = 12, b = 8, l = 8)
B_pad <- pad(B_tagged, t = 8, r = 8,  b = 6, l = 12)
C_pad <- pad(C_tagged, t = 8, r = 8,  b = 8, l = 12)

final_plot <- (A_pad | (B_pad / C_pad + plot_layout(heights = c(1,0.7)))) +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(12, 12, 12, 12)))


final_plot

ggsave("~/Desktop/PhD/projects/MCP_MDD_colocalisation/output/MiXeR/results.png", plot = final_plot, width = 12, height = 7, dpi = 600)

