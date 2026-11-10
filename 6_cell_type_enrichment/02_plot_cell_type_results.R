library(dplyr)
library(stringr)
library(readr)
library(purrr)
library(ggplot2)

# Get MCP results
MCP_files <- list.files(
  path = "./output/cell_type_enrichment/FUMA_celltype_MCP_Siletti/",
  pattern = ".gsa\\.out$",
  full.names = TRUE
)

MCP_cell_type_res <- map_dfr(MCP_files, function(file) {
  read.table(file, comment = "#", header = T) %>%
    mutate(source_file = basename(file))
})

# Get MDD results
MDD_files <- list.files(
  path = "./output/cell_type_enrichment/FUMA_celltype_MDD_Siletti//",
  pattern = ".gsa\\.out$",
  full.names = TRUE
)

MDD_cell_type_res <- map_dfr(MDD_files, function(file) {
  read.table(file, comment = "#", header = T) %>%
    mutate(source_file = basename(file))
})


MCP_cell_type_res <- na.omit(MCP_cell_type_res)
MDD_cell_type_res <- na.omit(MDD_cell_type_res)

MCP_cell_type_res$p_adjust_new <- p.adjust(MCP_cell_type_res$P, method = "bonferroni")
MDD_cell_type_res$p_adjust_new <- p.adjust(MDD_cell_type_res$P, method = "bonferroni")

MCP_cell_type_res_sig <- MCP_cell_type_res[MCP_cell_type_res$p_adjust_new < 0.05, ]
MDD_cell_type_res_sig <- MDD_cell_type_res[MDD_cell_type_res$p_adjust_new < 0.05, ]

MCP_cell_type_res_nominal <- MCP_cell_type_res[MCP_cell_type_res$P < 0.05, ]
MDD_cell_type_res_nominal <- MDD_cell_type_res[MDD_cell_type_res$P < 0.05, ]

cell_type_res_nominal <- inner_join(MCP_cell_type_res_nominal, MDD_cell_type_res_nominal, 
                                    by = c("source_file", "FULL_NAME", "VARIABLE", "TYPE"), suffix = c("_MCP", "_MDD")) %>%
  mutate(
    subregion = str_extract(source_file, "(?<=Siletti_)[^_]+"),
    brain_region = str_extract(subregion, "^[^.]+")
  ) %>%
  rename(
    cell_type = FULL_NAME
  ) %>%
  select(brain_region, subregion, cell_type, BETA_MDD, BETA_MCP, SE_MDD, SE_MCP, P_MDD, P_MCP, p_adjust_new_MDD, p_adjust_new_MCP) %>%
  arrange(brain_region, subregion, cell_type)
  
# Save nominal results

results_save <- cell_type_res_nominal %>%
  rename(
         `Brain region` = brain_region,
         `Brain subregion` = subregion,
         `Cell type` = cell_type,
         `MDD beta` = BETA_MDD,
         `MCP beta` = BETA_MCP,
         `MDD SE` = SE_MDD,
         `MCP SE` = SE_MCP,
         `MDD p` = P_MDD,
         `MCP p` = P_MCP,
         `MDD p adjsuted` = p_adjust_new_MDD,
         `MCP p adjsuted` = p_adjust_new_MCP
         ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  )

write.csv(results_save, "./output/cell_type_enrichment/Siletti_nominal_new.csv", row.names = F)

# Plot findings

heatmap_data <- cell_type_res_nominal %>%
  mutate(
    mean_beta = rowMeans(
      cbind(BETA_MDD, BETA_MCP),
      na.rm = TRUE
    )
  ) %>%
  group_by(brain_region, subregion, cell_type) %>%
  summarise(
    mean_beta = mean(mean_beta, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    subregion = factor(
      subregion,
      levels = rev(unique(subregion[order(brain_region, subregion)]))
    )
  )

ggplot(heatmap_data,
       aes(x = cell_type, y = subregion, fill = mean_beta)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_stepsn(
    colours = c(
      "#3B0F70",
      "#6A51A3",
      "#2C7FB8",
      "#41B6C4",
      "#A1DAB4",
      "#FFFFB2",
      "#FE9929",
      "#CC4C02",
      "#7F0000"
    ),
    n.breaks = 9,
    name = "Mean beta\n(MDD and MCP)"
  ) +
  labs(
    title = "",
    x = "Cell type",
    y = "Brain subregion"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
    axis.text.y = element_text(size = 5),
    panel.grid = element_blank()
  ) + coord_flip()

