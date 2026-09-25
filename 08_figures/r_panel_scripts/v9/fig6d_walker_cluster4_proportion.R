#!/usr/bin/env Rscript
# Fig 6.d (v9) — Cluster 4 proportion across 4 Walker treatments with Fisher OR (vs uninjured).
# 1:1 aspect. Bar with error bars (Fisher OR 95% CIs annotated).
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(readr); library(tibble); library(ggplot2); library(yaml)
})

PROJ  <- "/project/def-srussell/nirosha_xeno"
SEUR  <- file.path(PROJ, "comparators/walker/walker_integrated_v2.rds")
FIG_D <- file.path(PROJ, "figures_v9/fig6")
dir.create(FIG_D, showWarnings = FALSE, recursive = TRUE)

TS  <- format(Sys.time(), "%Y%m%dT%H%M%S")
LOG <- file(file.path(FIG_D, "fig6d_log.txt"), "w")
log <- function(...) { m <- paste0(..., collapse=""); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
log("== Fig 6.d (v9) Cluster 4 proportion + Fisher OR — ", TS, " ==")

log("Loading Seurat object: ", SEUR, " ...")
seur <- readRDS(SEUR)

cells_df <- tibble(
  cluster   = as.integer(as.character(seur$seurat_clusters)),
  treatment = as.character(seur$treatment)
)
n_cells <- nrow(cells_df)

TARGET <- 4L
trt_levels <- c("uninjured","tail","limb","muscimol")
trt_levels <- trt_levels[trt_levels %in% cells_df$treatment]

frac_df <- cells_df %>% group_by(treatment) %>%
  summarise(n_cells = n(),
            n_target = sum(cluster == TARGET),
            pct_target = 100 * n_target / n_cells,
            .groups = "drop") %>%
  arrange(match(treatment, trt_levels))
log("Cluster ", TARGET, " fraction per treatment:")
log(capture.output(print(frac_df)))

# Fisher OR vs uninjured
ref <- "uninjured"
ref_row <- frac_df %>% filter(treatment == ref)
or_rows <- list()
for (tr in trt_levels) {
  if (tr == ref) {
    or_rows[[tr]] <- tibble(treatment = tr, OR = 1, CI95_lo = NA_real_, CI95_hi = NA_real_, p = NA_real_)
    next
  }
  this_row <- frac_df %>% filter(treatment == tr)
  m <- matrix(c(this_row$n_target,             ref_row$n_target,
                this_row$n_cells - this_row$n_target,
                ref_row$n_cells  - ref_row$n_target),
              nrow = 2, byrow = FALSE,
              dimnames = list(c(tr, ref), c("target","other")))
  ft <- fisher.test(m)
  or_rows[[tr]] <- tibble(
    treatment = tr,
    OR        = unname(ft$estimate),
    CI95_lo   = unname(ft$conf.int[1]),
    CI95_hi   = unname(ft$conf.int[2]),
    p         = ft$p.value
  )
}
or_tbl <- bind_rows(or_rows)
plot_df <- frac_df %>% left_join(or_tbl, by = "treatment") %>%
  mutate(treatment = factor(treatment, levels = trt_levels))
log("Plot data:")
log(capture.output(print(plot_df)))
write_tsv(plot_df, file.path(FIG_D, "fig6d_cluster4_proportion_fisher.tsv"))

# Annotation labels: OR (95% CI), p
plot_df <- plot_df %>% mutate(
  ann = ifelse(treatment == ref,
               "OR = 1 (ref)",
               sprintf("OR = %.2f\n[%.2f, %.2f]\np = %.2g", OR, CI95_lo, CI95_hi, p))
)

TRT_COL <- c(uninjured = "#7F7F7F", tail = "#3B7DBF",
             limb      = "#D7263D", muscimol = "#9467BD")

panel_theme <- theme_classic(base_size = 11) +
  theme(axis.text       = element_text(size = 9.5, color = "black"),
        axis.title      = element_text(size = 10,  color = "black"),
        plot.title      = element_text(size = 11,  face = "plain", family = "sans", hjust = 0),
        legend.position = "none",
        plot.margin     = margin(8, 12, 6, 8))

y_max <- max(plot_df$pct_target) * 1.45

panel_d <- ggplot(plot_df, aes(x = treatment, y = pct_target, fill = treatment)) +
  geom_col(width = 0.7, color = "grey20", linewidth = 0.2) +
  geom_text(aes(y = pct_target + y_max * 0.04, label = ann),
            size = 2.7, lineheight = 0.95, vjust = 0, color = "black") +
  scale_fill_manual(values = TRT_COL) +
  scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0))) +
  labs(title = "Cluster 4 proportion across treatments",
       x = "Treatment",
       y = "Cluster 4 cells (%)") +
  panel_theme

ggsave(file.path(FIG_D, "d_panel.pdf"), panel_d, width = 130, height = 130, units = "mm")
ggsave(file.path(FIG_D, "d_panel.png"), panel_d, width = 130, height = 130, units = "mm", dpi = 300)

key_claim <- paste(sprintf("%s=%.2f%% (OR=%s)",
                           plot_df$treatment, plot_df$pct_target,
                           ifelse(is.na(plot_df$OR) | plot_df$treatment == ref,
                                  "ref", sprintf("%.2f", plot_df$OR))),
                   collapse = "; ")

write_yaml(list(
  panel              = "Fig 6d",
  title              = "Cluster 4 proportion across treatments",
  axis_x_label       = "Treatment",
  axis_y_label       = "Cluster 4 cells (%)",
  data_source_path   = SEUR,
  analysis_layer     = "scRNA cluster proportion + Fisher exact OR vs uninjured (Walker v2 re-integration)",
  n_obs              = n_cells,
  reference_treatment = ref,
  key_numerical_claim = key_claim,
  fisher_table       = lapply(seq_len(nrow(plot_df)), function(i) as.list(plot_df[i, ])),
  generated_utc      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script             = "scripts/v9/fig6d_walker_cluster4_proportion.R",
  styling            = list(width_mm = 130, height_mm = 130,
                            treatment_colors = as.list(TRT_COL))
), file.path(FIG_D, "panel_meta_d.yaml"))

log("Fig 6.d render COMPLETE.")
close(LOG)
