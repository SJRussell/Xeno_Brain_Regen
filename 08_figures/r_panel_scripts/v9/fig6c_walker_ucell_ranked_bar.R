#!/usr/bin/env Rscript
# Fig 6.c (v9) — Walker etv1+/glutamatergic UCell score per cluster, ranked descending.
# 4:3 wide. Plain 11pt sans title. No subtitle.
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(ggplot2); library(yaml)
})

PROJ  <- "/project/def-srussell/nirosha_xeno"
UCEL  <- file.path(PROJ, "comparators/walker/ucell_signatures_per_cluster_v2.tsv")
FIG_D <- file.path(PROJ, "figures_v9/fig6")
dir.create(FIG_D, showWarnings = FALSE, recursive = TRUE)

TS  <- format(Sys.time(), "%Y%m%dT%H%M%S")
LOG <- file(file.path(FIG_D, "fig6c_log.txt"), "w")
log <- function(...) { m <- paste0(..., collapse=""); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
log("== Fig 6.c (v9) Walker UCell ranked bar — ", TS, " ==")

ucell <- read_tsv(UCEL, show_col_types = FALSE)
score_col   <- "glut_etv1plus_regen_UCellv2"
cluster_col <- "seurat_clusters"

ranked <- ucell %>%
  group_by(.data[[cluster_col]]) %>%
  summarise(mean_score = mean(.data[[score_col]], na.rm = TRUE),
            sem        = sd(.data[[score_col]], na.rm = TRUE) / sqrt(n()),
            n_treatments = n(), .groups = "drop") %>%
  arrange(desc(mean_score)) %>%
  mutate(rank = row_number(),
         cluster_lbl = factor(.data[[cluster_col]], levels = .data[[cluster_col]]))
log("Ranked clusters by mean glut_etv1+_regen UCell across treatments:")
log(capture.output(print(ranked)))

top    <- ranked$mean_score[1];   top_cl <- ranked[[cluster_col]][1]
runner <- ranked$mean_score[2];   run_cl <- ranked[[cluster_col]][2]
fold   <- top / runner
write_tsv(ranked, file.path(FIG_D, "fig6c_cluster_ranking.tsv"))

ranked <- ranked %>% mutate(
  highlight = case_when(
    .data[[cluster_col]] == 4 ~ "Cluster 4 (rank 1)",
    .data[[cluster_col]] == 7 ~ "Cluster 7 (rank 2)",
    TRUE ~ "Other clusters"),
  highlight = factor(highlight, levels = c("Cluster 4 (rank 1)", "Cluster 7 (rank 2)", "Other clusters")))

HL_COL <- c("Cluster 4 (rank 1)" = "#D7263D",
            "Cluster 7 (rank 2)" = "#F4A261",
            "Other clusters"     = "grey75")

panel_theme <- theme_classic(base_size = 11) +
  theme(axis.text       = element_text(size = 9.5, color = "black"),
        axis.title      = element_text(size = 10,  color = "black"),
        plot.title      = element_text(size = 11,  face = "plain", family = "sans", hjust = 0),
        legend.position = "right",
        legend.text     = element_text(size = 8.5),
        legend.title    = element_text(size = 9),
        plot.margin     = margin(8, 14, 6, 8))

panel_c <- ggplot(ranked, aes(x = cluster_lbl, y = mean_score, fill = highlight)) +
  geom_col(width = 0.75, color = "grey20", linewidth = 0.2) +
  geom_errorbar(aes(ymin = pmax(0, mean_score - sem),
                    ymax = mean_score + sem),
                width = 0.2, linewidth = 0.35, color = "grey25") +
  scale_fill_manual(values = HL_COL, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Walker etv1+/glutamatergic UCell score by cluster",
       x = "Cluster (ranked by mean UCell, descending)",
       y = "Mean UCell score") +
  panel_theme

ggsave(file.path(FIG_D, "c_panel.pdf"), panel_c, width = 200, height = 110, units = "mm")
ggsave(file.path(FIG_D, "c_panel.png"), panel_c, width = 200, height = 110, units = "mm", dpi = 300)

key_claim <- sprintf("Cluster %d leads at mean UCell = %.4f vs runner-up cluster %d = %.4f (%.2fx lead) across %d clusters.",
                     top_cl, top, run_cl, runner, fold, nrow(ranked))

write_yaml(list(
  panel              = "Fig 6c",
  title              = "Walker etv1+/glutamatergic UCell score by cluster",
  axis_x_label       = "Cluster (ranked by mean UCell, descending)",
  axis_y_label       = "Mean UCell score",
  data_source_path   = UCEL,
  analysis_layer     = "Walker scRNA UCell mean across 4 treatments per Louvain cluster (v2 re-integration)",
  n_obs              = nrow(ranked),
  key_numerical_claim = key_claim,
  top_cluster        = unname(top_cl),
  top_score          = signif(top, 4),
  runner_up_cluster  = unname(run_cl),
  runner_up_score    = signif(runner, 4),
  fold_lead          = signif(fold, 4),
  generated_utc      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script             = "scripts/v9/fig6c_walker_ucell_ranked_bar.R",
  styling            = list(width_mm = 200, height_mm = 110,
                            highlight_colors = as.list(HL_COL))
), file.path(FIG_D, "panel_meta_c.yaml"))

log("Fig 6.c render COMPLETE.")
close(LOG)
