#!/usr/bin/env Rscript
# Fig 1.D — BRAIN PCA, 12 samples, homeolog-summed voom-normalized, PC1 vs PC2
# Single panel; treatment-colored, timepoint-shaped; convex hulls per (treatment, timepoint).
suppressPackageStartupMessages({
  library(limma); library(dplyr); library(readr); library(tibble); library(ggplot2); library(yaml)
})

PROJ  <- "/project/def-srussell/nirosha_xeno"
VOOM  <- file.path(PROJ, "deg/voom_BRAIN_homeolog_sum_strand_corrected.rds")
META  <- file.path(PROJ, "deg/sample_metadata.tsv")
FIG1  <- file.path(PROJ, "figures/fig1"); dir.create(FIG1, showWarnings = FALSE, recursive = TRUE)

TS  <- format(Sys.time(), "%Y%m%dT%H%M%S")
LOG <- file(file.path(FIG1, "fig1D_log.txt"), "w")
log <- function(...) { m <- paste0(..., collapse=""); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
log("== Fig 1.D BRAIN PCA — ", TS, " ==")

# ---- SO-12 preservation ----
for (p in c("D_panel.pdf","D_panel.png","panel_metadata_D.yaml")) {
  src <- file.path(FIG1, p)
  if (file.exists(src)) {
    ext  <- tools::file_ext(p); base <- tools::file_path_sans_ext(p)
    dst  <- file.path(FIG1, paste0(base, "_pre_legend_polish_", TS, ".", ext))
    file.rename(src, dst); log("SO-12 preserved (pre-polish): ", p, " -> ", basename(dst))
  }
}

# ---- Load voom (homeolog-summed, strand-corrected) ----
obj <- readRDS(VOOM)
log("Top-level class: ", paste(class(obj), collapse=","), "  names: ", paste(names(obj), collapse=", "))
# The strand-corrected homeolog-sum RDS bundles voom + fit + results; the EList is at $v
v <- if (is.list(obj) && "v" %in% names(obj)) obj$v else obj
log("Voom EList class: ", paste(class(v), collapse=","))
emat <- NULL
if (!is.null(v$E)) emat <- v$E
if (is.null(emat) && inherits(v, "EList")) emat <- as.matrix(v$E)
if (is.null(emat)) stop("Could not locate E matrix in voom RDS")
log("E matrix dim (genes x samples): ", paste(dim(emat), collapse=" x "))
log("First 4 sample cols: ", paste(head(colnames(emat),4), collapse=", "))

# ---- Subset to 12 brain samples ----
meta <- read_tsv(META, show_col_types = FALSE) %>%
  filter(tissue == "BRAIN") %>%
  mutate(treatment = factor(treatment, levels = c("ND","CT")),
         time      = factor(time,      levels = c("24HPA","7DPA")))
brain_ids <- meta$sample
stopifnot(length(brain_ids) == 12)
miss <- setdiff(brain_ids, colnames(emat))
if (length(miss)) stop("Missing BRAIN samples in voom E matrix: ", paste(miss, collapse=", "))
emat <- emat[, brain_ids]
log("After subset, E dim: ", paste(dim(emat), collapse=" x "))

# ---- PCA on top variable genes (standard PCA-on-voom convention) ----
# Use top 5000 most-variable genes if available, else all
g_var <- apply(emat, 1, var)
g_top <- names(sort(g_var, decreasing = TRUE))[1:min(5000, length(g_var))]
emat_top <- emat[g_top, ]
pc <- prcomp(t(emat_top), center = TRUE, scale. = FALSE)
varexp <- (pc$sdev^2 / sum(pc$sdev^2)) * 100
log(sprintf("PC1 variance = %.1f%%, PC2 variance = %.1f%%, top-5k genes used", varexp[1], varexp[2]))

pca_df <- as_tibble(pc$x[, 1:4], rownames = "sample") %>%
  left_join(meta, by = "sample") %>%
  mutate(group = paste0(treatment, "_", time))

write_tsv(pca_df, file.path(PROJ, "deg/pca_BRAIN_homeolog_sum_strand_corrected.tsv"))
log("Wrote deg/pca_BRAIN_homeolog_sum_strand_corrected.tsv")

# ---- Convex hulls per (treatment, time) ----
hull_df <- pca_df %>% group_by(group, treatment, time) %>%
  slice(chull(PC1, PC2)) %>% ungroup()

# ---- Panel ----
TRT_COL <- c(ND = "#3B7DBF", CT = "#D9614C")
TIME_SH <- c(`24HPA` = 16, `7DPA` = 17)

panel_theme <- theme_classic(base_size = 10) +
  theme(axis.text  = element_text(size = 9,  color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        plot.title    = element_text(size = 10, face = "plain", hjust = 0),
        plot.subtitle = element_text(size = 9,  color = "grey30"),
        legend.text   = element_text(size = 8),
        legend.title  = element_text(size = 9),
        plot.margin   = margin(8, 14, 6, 8))

ct24_pc2 <- pca_df %>% filter(treatment == "CT", time == "24HPA") %>% arrange(PC2) %>% pull(PC2)
log(sprintf("CT 24HPA PC2 values (sorted): %.1f, %.1f, %.1f", ct24_pc2[1], ct24_pc2[2], ct24_pc2[3]))
outlier_subtitle <- sprintf(
  "Convex hulls per (treatment, timepoint); top 5000 variable genes. NOTE: one CT 24HPA brain sample sits at PC2 = %+.0f vs the other two at PC2 = %+.0f / %+.0f - large within-group PC2 variance (sensitivity-drop candidate).",
  ct24_pc2[1], ct24_pc2[2], ct24_pc2[3])

panel_D <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_polygon(data = hull_df,
               aes(x = PC1, y = PC2, group = group, fill = treatment, linetype = time),
               color = "grey25", alpha = 0.18, linewidth = 0.45,
               show.legend = c(linetype = TRUE, fill = FALSE)) +
  geom_point(aes(color = treatment, shape = time), size = 3.2, stroke = 0.6) +
  scale_color_manual(values = TRT_COL, name = "Treatment") +
  scale_fill_manual(values  = TRT_COL, guide = "none") +
  scale_shape_manual(values = TIME_SH, name = "Timepoint") +
  scale_linetype_manual(values = c(`24HPA` = "dashed", `7DPA` = "solid"),
                        name = "Timepoint") +
  guides(linetype = guide_legend(override.aes = list(color = "grey25", fill = NA)),
         shape    = guide_legend(override.aes = list(color = "grey25"))) +
  labs(title    = "BRAIN PCA, 12 samples, homeolog-summed voom-normalized",
       subtitle = outlier_subtitle,
       x = sprintf("PC1 (%.1f%% variance)", varexp[1]),
       y = sprintf("PC2 (%.1f%% variance)", varexp[2]),
       color = "Treatment", shape = "Timepoint") +
  panel_theme

ggsave(file.path(FIG1, "D_panel.pdf"), panel_D, width = 170, height = 100, units = "mm")
ggsave(file.path(FIG1, "D_panel.png"), panel_D, width = 170, height = 100, units = "mm", dpi = 300)

write_yaml(list(
  panel = "Fig1.D",
  role  = "BRAIN PCA, 12 samples, homeolog-summed voom-normalized; PC1 vs PC2; treatment colored, time shaped; convex hulls per (treatment, time)",
  generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script = "scripts/fig1D_brain_pca.R",
  inputs = list(voom_rds = VOOM, sample_metadata_tsv = META),
  pca = list(n_samples = 12, n_genes_used = nrow(emat_top), gene_selection = "top 5000 by variance",
             pc1_pct = signif(varexp[1], 4), pc2_pct = signif(varexp[2], 4),
             pc3_pct = signif(varexp[3], 4), pc4_pct = signif(varexp[4], 4)),
  hull_groups = sort(unique(pca_df$group)),
  styling = list(width_mm = 170, height_mm = 100,
                 treatment_colors = list(ND = "#3B7DBF", CT = "#D9614C"),
                 time_shapes = list(`24HPA` = "circle", `7DPA` = "triangle"),
                 hull_alpha = 0.18, hull_outline = "grey25 linewidth 0.4",
                 hull_linetype_by_time = list(`24HPA` = "dashed", `7DPA` = "solid"),
                 regression_line = FALSE)
), file.path(FIG1, "panel_metadata_D.yaml"))

log("\nFig 1.D render COMPLETE.")
close(LOG)
