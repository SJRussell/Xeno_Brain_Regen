#!/usr/bin/env Rscript
# Stage C: QC + outlier triage
# - Reads counts/star_gene_counts.tsv
# - Builds sample metadata
# - Filters low-count genes (CPM > 1 in >=3 samples)
# - TMM normalize via edgeR
# - VST + log2-CPM
# - Per-tissue PCA
# - Spearman + Pearson sample-distance heatmaps
# - Outlier detection (confirm CT_T3_7DPA_BRAIN)
# - Save RDS objects + figure PDFs

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(DESeq2)
    library(ggplot2)
    library(pheatmap)
    library(RColorBrewer)
    library(matrixStats)
})

PROJ <- "/project/def-srussell/nirosha_xeno"
COUNTS_DIR <- file.path(PROJ, "counts")
QC_DIR <- file.path(PROJ, "qc")
DEG_DIR <- file.path(PROJ, "deg")
FIG_DIR <- file.path(PROJ, "figures", "stage_C")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DEG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("[", format(Sys.time(), "%FT%T%z"), "] Loading STAR gene counts...\n")
counts_path <- file.path(COUNTS_DIR, "star_gene_counts.tsv")
mat <- read.table(counts_path, header = TRUE, sep = "\t", row.names = 1,
                  check.names = FALSE, stringsAsFactors = FALSE)
cat("  Counts matrix:", nrow(mat), "genes x", ncol(mat), "samples\n")

# ----- Build sample metadata from column names -----
# colnames look like: CT_T1_11HPA_BLAST_S3
parse_meta <- function(s) {
    parts <- strsplit(s, "_")[[1]]
    list(
        treatment = parts[1],     # ND or CT
        rep       = parts[2],     # T1, T2, T3
        time      = parts[3],     # 11HPA, 24HPA, 7DPA
        tissue    = parts[4]      # BLAST, BRAIN
    )
}
meta_list <- lapply(colnames(mat), parse_meta)
meta <- data.frame(
    sample    = colnames(mat),
    treatment = factor(sapply(meta_list, `[[`, "treatment"), levels = c("ND", "CT")),
    rep       = factor(sapply(meta_list, `[[`, "rep"),       levels = c("T1", "T2", "T3")),
    time      = factor(sapply(meta_list, `[[`, "time"),      levels = c("11HPA", "24HPA", "7DPA")),
    tissue    = factor(sapply(meta_list, `[[`, "tissue"),    levels = c("BLAST", "BRAIN")),
    stringsAsFactors = FALSE
)
meta$group <- factor(paste(meta$treatment, meta$time, meta$tissue, sep = "_"))
cat("  Sample groups:\n")
print(table(meta$treatment, meta$time, meta$tissue))
write.table(meta, file.path(DEG_DIR, "sample_metadata.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ----- Filter low-count genes -----
# Per pipeline_plan: CPM > 1 in >= 3 samples (matches BBcore 2019 threshold)
y_full <- DGEList(counts = as.matrix(mat), samples = meta)
cat("\n[", format(Sys.time(), "%FT%T%z"), "] Library sizes (counted reads):\n")
print(summary(y_full$samples$lib.size))

cpm_mat <- cpm(y_full)
keep <- rowSums(cpm_mat > 1) >= 3
cat("\nGene filter: CPM > 1 in >= 3 samples\n")
cat("  Pre-filter:  ", nrow(y_full), "\n")
cat("  Post-filter: ", sum(keep), "\n")

y <- y_full[keep, , keep.lib.sizes = FALSE]

# ----- TMM normalize -----
y <- calcNormFactors(y, method = "TMM")
cat("\nTMM norm factors (range): [",
    round(min(y$samples$norm.factors), 3), ",",
    round(max(y$samples$norm.factors), 3), "]\n")

# ----- VST via DESeq2 -----
cat("\n[", format(Sys.time(), "%FT%T%z"), "] Computing VST via DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(mat[keep, ])),
                              colData = meta,
                              design = ~ tissue + time + treatment)
vsd <- vst(dds, blind = TRUE)
saveRDS(vsd, file.path(DEG_DIR, "vsd_blind.rds"))

# ----- Sample-distance heatmaps -----
cat("\n[", format(Sys.time(), "%FT%T%z"), "] Sample-distance heatmaps...\n")
dist_mat <- as.matrix(dist(t(assay(vsd))))
rownames(dist_mat) <- colnames(dist_mat) <- paste(meta$treatment, meta$rep, meta$time, meta$tissue, sep = "_")
ann <- data.frame(
    Treatment = meta$treatment,
    Time      = meta$time,
    Tissue    = meta$tissue,
    row.names = rownames(dist_mat)
)

pdf(file.path(FIG_DIR, "sample_distance_euclid.pdf"), width = 10, height = 9)
pheatmap(dist_mat,
         clustering_distance_rows = as.dist(dist_mat),
         clustering_distance_cols = as.dist(dist_mat),
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         annotation_col = ann,
         main = "Sample-distance (Euclidean on VST)")
dev.off()

# Spearman correlation distance
spearman_corr <- cor(assay(vsd), method = "spearman")
spearman_dist <- 1 - spearman_corr
# Relabel rownames/colnames to match the annotation_col
rownames(spearman_dist) <- colnames(spearman_dist) <- rownames(dist_mat)
pdf(file.path(FIG_DIR, "sample_distance_spearman.pdf"), width = 10, height = 9)
pheatmap(spearman_dist,
         clustering_distance_rows = as.dist(spearman_dist),
         clustering_distance_cols = as.dist(spearman_dist),
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         annotation_col = ann,
         main = "Sample-distance (1 - Spearman correlation)")
dev.off()

# ----- PCA per tissue -----
cat("\n[", format(Sys.time(), "%FT%T%z"), "] PCA per tissue...\n")

run_pca <- function(vsd_subset, meta_subset, tissue_name) {
    mat_t <- assay(vsd_subset)
    # Use top 2000 variable genes
    vars <- rowVars(mat_t)
    top <- order(vars, decreasing = TRUE)[seq_len(min(2000, nrow(mat_t)))]
    pca <- prcomp(t(mat_t[top, ]), scale. = FALSE)
    pcs <- as.data.frame(pca$x[, 1:4])
    pcs$sample    <- rownames(pcs)
    pcs$treatment <- meta_subset$treatment
    pcs$rep       <- meta_subset$rep
    pcs$time      <- meta_subset$time
    pve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
    list(pcs = pcs, pve = pve)
}

for (tiss in c("BRAIN", "BLAST")) {
    sub_idx <- meta$tissue == tiss
    vsd_sub <- vsd[, sub_idx]
    meta_sub <- meta[sub_idx, ]
    pca_res <- run_pca(vsd_sub, meta_sub, tiss)

    p1 <- ggplot(pca_res$pcs, aes(PC1, PC2, color = treatment, shape = time)) +
        geom_point(size = 4, alpha = 0.85) +
        geom_text(aes(label = paste(rep)), vjust = -1.2, size = 3) +
        labs(
            title = paste0("PCA — ", tiss, " samples (top 2000 variable genes, VST)"),
            x = sprintf("PC1 (%.1f%%)", pca_res$pve[1]),
            y = sprintf("PC2 (%.1f%%)", pca_res$pve[2])
        ) +
        theme_minimal(base_size = 12) +
        scale_color_manual(values = c("ND" = "#1f77b4", "CT" = "#d62728"))

    ggsave(file.path(FIG_DIR, paste0("pca_pc1_pc2_", tiss, ".pdf")), p1,
           width = 8, height = 7)
    write.table(pca_res$pcs, file.path(DEG_DIR, paste0("pca_", tiss, ".tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  %s: PVE [PC1=%.1f%%, PC2=%.1f%%, PC3=%.1f%%, PC4=%.1f%%]\n",
                tiss, pca_res$pve[1], pca_res$pve[2], pca_res$pve[3], pca_res$pve[4]))
}

# ----- Outlier detection -----
# Compute per-sample distance to the centroid of its same-condition peers (group)
cat("\n[", format(Sys.time(), "%FT%T%z"), "] Outlier diagnostics...\n")
cat("\nPer-group within-group mean Euclidean distance:\n")

within_group_dist <- function(mat_t, meta_t) {
    out <- data.frame(
        sample = colnames(mat_t),
        group  = as.character(meta_t$group),
        mean_dist_to_peers = NA_real_,
        median_dist_to_peers = NA_real_,
        stringsAsFactors = FALSE
    )
    dmat <- as.matrix(dist(t(mat_t)))
    for (i in seq_along(out$sample)) {
        peers <- which(meta_t$group == meta_t$group[i])
        peers <- setdiff(peers, i)
        if (length(peers) > 0) {
            out$mean_dist_to_peers[i]   <- mean(dmat[i, peers])
            out$median_dist_to_peers[i] <- median(dmat[i, peers])
        }
    }
    out
}
od <- within_group_dist(assay(vsd), meta)
od$z_within_group <- ave(od$mean_dist_to_peers, od$group,
                         FUN = function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
write.table(od, file.path(DEG_DIR, "outlier_diagnostics.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nSamples with |z| > 1 within-group:\n")
print(od[!is.na(od$z_within_group) & abs(od$z_within_group) > 1, ])

# Within-tissue PCA distance to tissue centroid (catches global outliers)
cat("\nPer-tissue PC1+PC2 distance to tissue centroid:\n")
for (tiss in c("BRAIN", "BLAST")) {
    pcs <- read.table(file.path(DEG_DIR, paste0("pca_", tiss, ".tsv")),
                      header = TRUE, sep = "\t")
    centroid <- colMeans(pcs[, c("PC1", "PC2")])
    pcs$dist_to_centroid <- sqrt((pcs$PC1 - centroid["PC1"])^2 + (pcs$PC2 - centroid["PC2"])^2)
    pcs <- pcs[order(-pcs$dist_to_centroid), ]
    cat("\n", tiss, " samples by distance from PC1+PC2 centroid (largest first):\n")
    print(pcs[, c("sample", "treatment", "rep", "time", "PC1", "PC2", "dist_to_centroid")])
}

# ----- Library-size + counted-rate plot -----
cat("\n[", format(Sys.time(), "%FT%T%z"), "] Library-size summary plot...\n")
star_summary <- read.table(file.path(COUNTS_DIR, "star_gene_summary.tsv"),
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)
star_summary$tissue <- ifelse(grepl("BRAIN", star_summary$sample), "BRAIN", "BLAST")
star_summary$treatment <- ifelse(grepl("^CT", star_summary$sample), "CT", "ND")
star_summary$time <- sub(".*_(11HPA|24HPA|7DPA)_.*", "\\1", star_summary$sample)

p_lib <- ggplot(star_summary, aes(x = sample, y = counted_reads / 1e6,
                                   fill = treatment, color = time)) +
    geom_col() +
    facet_grid(. ~ tissue, scales = "free_x", space = "free_x") +
    coord_flip() +
    labs(title = "Counted reads per sample (M)",
         x = NULL, y = "Counted reads (millions)") +
    theme_minimal(base_size = 9) +
    scale_fill_manual(values = c("ND" = "#1f77b4", "CT" = "#d62728"))
ggsave(file.path(FIG_DIR, "library_size.pdf"), p_lib, width = 9, height = 9)

# ----- Save the full DGEList for Stage D -----
saveRDS(y, file.path(DEG_DIR, "dgelist_filtered_tmm.rds"))
saveRDS(meta, file.path(DEG_DIR, "sample_metadata.rds"))

cat("\n[", format(Sys.time(), "%FT%T%z"), "] Stage C done. Outputs:\n")
cat("  ", FIG_DIR, "/sample_distance_{euclid,spearman}.pdf\n")
cat("  ", FIG_DIR, "/pca_pc1_pc2_{BRAIN,BLAST}.pdf\n")
cat("  ", FIG_DIR, "/library_size.pdf\n")
cat("  ", DEG_DIR, "/sample_metadata.tsv\n")
cat("  ", DEG_DIR, "/pca_{BRAIN,BLAST}.tsv\n")
cat("  ", DEG_DIR, "/outlier_diagnostics.tsv\n")
cat("  ", DEG_DIR, "/dgelist_filtered_tmm.rds\n")
cat("  ", DEG_DIR, "/vsd_blind.rds\n")
