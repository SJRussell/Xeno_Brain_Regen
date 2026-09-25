#!/usr/bin/env Rscript
# Fig 2.C — BLAST PCA, 18 samples, homeolog-summed voom-normalized (strand-corrected), PC1 vs PC2
# Single panel; treatment-colored, timepoint-shaped (3 shapes); convex hulls per (treatment, timepoint).
suppressPackageStartupMessages({
  library(limma); library(dplyr); library(readr); library(tibble); library(ggplot2); library(yaml)
})

PROJ <- "/project/def-srussell/nirosha_xeno"
VOOM <- file.path(PROJ, "deg/voom_BLAST_homeolog_sum_strand_corrected.rds")
META <- file.path(PROJ, "deg/sample_metadata.tsv")
FIG2 <- file.path(PROJ, "figures/fig2"); dir.create(FIG2, showWarnings = FALSE, recursive = TRUE)

TS  <- format(Sys.time(), "%Y%m%dT%H%M%S")
LOG <- file(file.path(FIG2, "fig2C_log.txt"), "w")
log <- function(...) { m <- paste0(..., collapse=""); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
log("== Fig 2.C BLAST PCA — ", TS, " ==")

# ---- SO-12 preservation ----
for (p in c("C_panel.pdf","C_panel.png","panel_metadata_C.yaml")) {
  src <- file.path(FIG2, p)
  if (file.exists(src)) {
    ext <- tools::file_ext(p); base <- tools::file_path_sans_ext(p)
    dst <- file.path(FIG2, paste0(base, "_", TS, ".", ext))
    file.rename(src, dst); log("SO-12 preserved: ", p, " -> ", basename(dst))
  }
}

# ---- Load voom EList from bundled RDS ----
obj <- readRDS(VOOM)
log("Top-level class: ", paste(class(obj), collapse=","), "  names: ", paste(names(obj), collapse=", "))
v <- if (is.list(obj) && "v" %in% names(obj)) obj$v else obj
emat <- v$E
log("E matrix dim (genes x samples): ", paste(dim(emat), collapse=" x "))

# ---- Subset to 18 BLAST samples ----
meta <- read_tsv(META, show_col_types = FALSE) %>%
  filter(tissue == "BLAST") %>%
  mutate(treatment = factor(treatment, levels = c("ND","CT")),
         time      = factor(time,      levels = c("11HPA","24HPA","7DPA")))
blast_ids <- meta$sample
stopifnot(length(blast_ids) == 18)
miss <- setdiff(blast_ids, colnames(emat))
if (length(miss)) stop("Missing BLAST samples in voom E: ", paste(miss, collapse=", "))
emat <- emat[, blast_ids]
log("After subset, E dim: ", paste(dim(emat), collapse=" x "))

# ---- PCA on top 5000 variable genes ----
g_var <- apply(emat, 1, var)
g_top <- names(sort(g_var, decreasing = TRUE))[1:min(5000, length(g_var))]
emat_top <- emat[g_top, ]
pc <- prcomp(t(emat_top), center = TRUE, scale. = FALSE)
varexp <- (pc$sdev^2 / sum(pc$sdev^2)) * 100
log(sprintf("PC1 variance = %.1f%%, PC2 variance = %.1f%%, top-5k genes used", varexp[1], varexp[2]))

pca_df <- as_tibble(pc$x[, 1:4], rownames = "sample") %>%
  left_join(meta, by = "sample") %>%
  mutate(group = paste0(treatment, "_", time))

write_tsv(pca_df, file.path(PROJ, "deg/pca_BLAST_homeolog_sum_strand_corrected.tsv"))
log("Wrote deg/pca_BLAST_homeolog_sum_strand_corrected.tsv")

# ---- Convex hulls per (treatment, time) — 6 groups, n=3 each ----
hull_df <- pca_df %>% group_by(group, treatment, time) %>%
  slice(chull(PC1, PC2)) %>% ungroup()

# ---- Styling ----
TRT_COL <- c(ND = "#3B7DBF", CT = "#D9614C")
TIME_SH <- c(`11HPA` = 15, `24HPA` = 16, `7DPA` = 17)
TIME_LT <- c(`11HPA` = "dotted", `24HPA` = "dashed", `7DPA` = "solid")

panel_theme <- theme_classic(base_size = 10) +
  theme(axis.text  = element_text(size = 9,  color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        plot.title    = element_text(size = 10, face = "plain", hjust = 0),
        plot.subtitle = element_text(size = 9,  color = "grey30"),
        legend.text   = element_text(size = 8),
        legend.title  = element_text(size = 9),
        plot.margin   = margin(8, 14, 6, 8))

panel_C <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_polygon(data = hull_df,
               aes(x = PC1, y = PC2, group = group, fill = treatment, linetype = time),
               color = "grey25", alpha = 0.18, linewidth = 0.4, show.legend = FALSE) +
  geom_point(aes(color = treatment, shape = time), size = 3.2, stroke = 0.6) +
  scale_color_manual(values = TRT_COL) +
  scale_fill_manual(values  = TRT_COL, guide = "none") +
  scale_shape_manual(values = TIME_SH) +
  scale_linetype_manual(values = TIME_LT, guide = "none") +
  labs(title    = "BLAST PCA, 18 samples, homeolog-summed voom-normalized",
       subtitle = "Convex hulls per (treatment, timepoint); top 5000 variable genes",
       x = sprintf("PC1 (%.1f%% variance)", varexp[1]),
       y = sprintf("PC2 (%.1f%% variance)", varexp[2]),
       color = "Treatment", shape = "Timepoint") +
  panel_theme

ggsave(file.path(FIG2, "C_panel.pdf"), panel_C, width = 170, height = 100, units = "mm")
ggsave(file.path(FIG2, "C_panel.png"), panel_C, width = 170, height = 100, units = "mm", dpi = 300)

write_yaml(list(
  panel = "Fig2.C",
  role  = "BLAST PCA, 18 samples, homeolog-summed voom-normalized (strand-corrected); PC1 vs PC2; treatment colored, time shaped; convex hulls per (treatment, time)",
  generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script = "scripts/fig2C_blast_pca.R",
  inputs = list(voom_rds = VOOM, sample_metadata_tsv = META),
  pca = list(n_samples = 18, n_genes_used = nrow(emat_top), gene_selection = "top 5000 by variance",
             pc1_pct = signif(varexp[1], 4), pc2_pct = signif(varexp[2], 4),
             pc3_pct = signif(varexp[3], 4), pc4_pct = signif(varexp[4], 4)),
  hull_groups = sort(unique(pca_df$group)),
  styling = list(width_mm = 170, height_mm = 100,
                 treatment_colors = list(ND = "#3B7DBF", CT = "#D9614C"),
                 time_shapes = list(`11HPA` = "square", `24HPA` = "circle", `7DPA` = "triangle"),
                 hull_alpha = 0.18,
                 hull_linetype_by_time = list(`11HPA` = "dotted", `24HPA` = "dashed", `7DPA` = "solid"),
                 regression_line = FALSE)
), file.path(FIG2, "panel_metadata_C.yaml"))

log("\nFig 2.C render COMPLETE.")
close(LOG)
