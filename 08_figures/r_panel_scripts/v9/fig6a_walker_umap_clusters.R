#!/usr/bin/env Rscript
# Fig 6.a (v9) — Walker axolotl brain scRNA UMAP, all cells, coloured by Louvain cluster.
# 23 categorical colours via ggsci d3 category20 palette extended; cluster 4 enlarged + outlined.
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(readr); library(tibble)
  library(ggplot2); library(yaml); library(ggsci); library(scales)
})

PROJ  <- "/project/def-srussell/nirosha_xeno"
SEUR  <- file.path(PROJ, "comparators/walker/walker_integrated_v2.rds")
FIG_D <- file.path(PROJ, "figures_v9/fig6")
dir.create(FIG_D, showWarnings = FALSE, recursive = TRUE)

TS  <- format(Sys.time(), "%Y%m%dT%H%M%S")
LOG <- file(file.path(FIG_D, "fig6a_log.txt"), "w")
log <- function(...) { m <- paste0(..., collapse=""); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
log("== Fig 6.a (v9) Walker UMAP by cluster — ", TS, " ==")

log("Loading Seurat object: ", SEUR, " (~1.6 GB) ...")
seur <- readRDS(SEUR)
log("Loaded; n_cells = ", ncol(seur), "  n_clusters = ", length(unique(seur$seurat_clusters)))

umap_df <- as.data.frame(Embeddings(seur, "umap")) %>% rownames_to_column("cell")
colnames(umap_df)[2:3] <- c("UMAP_1","UMAP_2")
umap_df$cluster <- as.integer(as.character(seur$seurat_clusters[umap_df$cell]))
n_cells <- nrow(umap_df); n_clusters <- length(unique(umap_df$cluster))
log("n_cells = ", n_cells, "  n_clusters = ", n_clusters)

cluster_levels <- sort(unique(umap_df$cluster))
umap_df$cluster_f <- factor(umap_df$cluster, levels = cluster_levels)

# 23 colours via ggsci d3 category20 extended by colorRampPalette
base_pal  <- pal_d3("category20")(20)
clust_pal <- setNames(colorRampPalette(base_pal)(length(cluster_levels)), cluster_levels)

TARGET <- 4L
bg_df  <- umap_df %>% filter(cluster != TARGET)
hl_df  <- umap_df %>% filter(cluster == TARGET)

# Centroid label for cluster 4 to draw the eye
cent_4 <- hl_df %>% summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

panel_theme <- theme_classic(base_size = 11) +
  theme(axis.text       = element_text(size = 9.5, color = "black"),
        axis.title      = element_text(size = 10,  color = "black"),
        plot.title      = element_text(size = 11,  face = "plain", family = "sans", hjust = 0),
        legend.position = "right",
        legend.text     = element_text(size = 8),
        legend.title    = element_text(size = 9),
        legend.key.size = unit(3.5, "mm"),
        plot.margin     = margin(8, 12, 6, 8))

panel_a <- ggplot() +
  geom_point(data = bg_df, aes(UMAP_1, UMAP_2, colour = cluster_f),
             size = 1.0, alpha = 0.7, stroke = 0) +
  geom_point(data = hl_df, aes(UMAP_1, UMAP_2, colour = cluster_f),
             size = 1.5, alpha = 0.9, stroke = 0.15) +
  annotate("text", x = cent_4$UMAP_1, y = cent_4$UMAP_2, label = "4",
           size = 4.0, fontface = "bold", colour = "black") +
  scale_colour_manual(values = clust_pal, name = "Cluster",
                      guide = guide_legend(override.aes = list(size = 2.5, alpha = 1),
                                           ncol = 2)) +
  coord_fixed() +
  labs(title = "Walker axolotl brain scRNA, 23 clusters",
       x = "UMAP 1", y = "UMAP 2") +
  panel_theme

ggsave(file.path(FIG_D, "a_panel.pdf"), panel_a, width = 150, height = 130, units = "mm")
ggsave(file.path(FIG_D, "a_panel.png"), panel_a, width = 150, height = 130, units = "mm", dpi = 300)

write_yaml(list(
  panel              = "Fig 6a",
  title              = "Walker axolotl brain scRNA, 23 clusters",
  axis_x_label       = "UMAP 1",
  axis_y_label       = "UMAP 2",
  data_source_path   = SEUR,
  analysis_layer     = "scRNA Louvain clusters (Walker v2 4-sample re-integration)",
  n_obs              = n_cells,
  n_clusters         = n_clusters,
  key_numerical_claim = sprintf("UMAP of %d Walker axolotl brain cells across %d Louvain clusters; cluster 4 (etv1+ winner per v2 re-integration) enlarged.",
                                n_cells, n_clusters),
  generated_utc      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script             = "scripts/v9/fig6a_walker_umap_clusters.R",
  styling            = list(width_mm = 150, height_mm = 130,
                            bg_point_size = 1.0, bg_alpha = 0.7,
                            target_point_size = 1.5, target_alpha = 0.9,
                            palette = "ggsci::pal_d3('category20') extended via colorRampPalette to 23")
), file.path(FIG_D, "panel_meta_a.yaml"))

log("Fig 6.a render COMPLETE.")
close(LOG)
