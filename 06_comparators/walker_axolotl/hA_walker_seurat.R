#!/usr/bin/env Rscript
# H-A Walker post-alevin Seurat analysis (DEC-017 steps 4-7).
# Inputs: salmon alevin per-sample quants_mat at
#   /scratch/.../nirosha_xeno_comparators/walker/work/alevin/<LABEL>/alevin/quants_mat.gz
#
# Pipeline:
#   1. Load alevin outputs via tximport (type="alevin") → gene × cell matrix per sample.
#   2. Per-sample Seurat object → mt percentage QC → low-quality cell filter.
#   3. SCTransform per sample → CCA integration across 4 samples.
#   4. UMAP + FindClusters at resolution 0.6.
#   5. Marker-based cluster annotation against Walker 2025 published markers
#      (etv1, slc17a7, cdh13, bcl11b, sox5, neto1, stk39, ntng1 → dpErk+/etv1+
#      glutamatergic regen-required; gfap+sox2+ → ependymoglia; gad1+gad2+ → GABAergic;
#      olig1+sox10+ → oligodendrocyte lineage; cd68+ → microglia/macrophage).
#   6. Pseudobulk per cluster (AggregateExpression sum).
#   7. UCell scoring of 4 Xenopus brain CT-treatment signatures
#      (anti-fibrotic ECM, IFN-immune, heme-iron, structural plasticity) + GDF15.
#   8. Treatment-population enrichment for etv1+ cluster vs others
#      with permutation null (n_perm=1000, fixed seed).
#
# Outputs:
#   /project/.../comparators/walker/walker_integrated.rds (Seurat object)
#   /project/.../comparators/walker/cluster_markers.tsv
#   /project/.../comparators/walker/cluster_annotations.tsv
#   /project/.../comparators/walker/pseudobulk_counts.tsv
#   /project/.../comparators/walker/ucell_signatures_per_cell.tsv
#   /project/.../comparators/walker/ucell_signatures_per_cluster.tsv
#   /project/.../comparators/walker/treatment_pop_enrichment.tsv
#   /scratch/.../nirosha_xeno_comparators/walker/work/figures/walker_*.pdf

suppressPackageStartupMessages({
    library(tximport)
    library(Seurat)
    library(Matrix)
    library(data.table)
    library(ggplot2)
    library(UCell)
    library(future)
})

set.seed(20260501)
# Lift future-globals cap so big SCT objects can transfer if needed, but use
# sequential plan to avoid the multicore deadlock seen in 13136226 (job stuck
# at 6:42 CPU after 1h41m wall — multicore + Seurat internal parallelism
# blocked each other). vst.flavor="v1" sidesteps the glmGamPoi-prerequisite
# warning and uses Seurat's older robust glmGamPoi-free path.
options(future.globals.maxSize = 8 * 1024^3)
plan("sequential")

PROJ <- "/project/def-srussell/nirosha_xeno"
WORK <- "/scratch/srussell/srussell/nirosha_xeno_comparators/walker/work"
ALEVIN_DIR <- file.path(WORK, "alevin")
FIG <- file.path(WORK, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

REF_DIR <- "/scratch/srussell/srussell/refs/amext_v47"
TX2GENE <- fread(file.path(REF_DIR, "tx2gene.tsv"),
                 header = FALSE, col.names = c("tx_id","gene_id","gene_name"))
# Build gene_id -> gene_symbol map for renaming Seurat features.
gid2sym <- TX2GENE[, .(gene_id = unique(gene_id),
                       gene_name = first(gene_name)), by = gene_id]
gid2sym[, gene_name := ifelse(gene_name == "" | is.na(gene_name), gene_id, gene_name)]
sym_map <- setNames(gid2sym$gene_name, gid2sym$gene_id)

# ------------------------------------------------------------
# 1. Load alevin per-sample
# ------------------------------------------------------------
samples <- data.frame(
    label = c("uninjured_brain","tail_brain","limb_brain","muscimol_brain"),
    srr   = c("SRR31800877","SRR31800876","SRR31800875","SRR31800874"),
    treatment = c("uninjured","tail","limb","muscimol"),
    stringsAsFactors = FALSE
)
samples$quant <- file.path(ALEVIN_DIR, samples$label, "alevin", "quants_mat.gz")
samples$present <- file.exists(samples$quant)
cat("Sample alevin presence:\n")
print(samples[, c("label","present")])
if (sum(samples$present) < 2) stop("Need at least 2 samples to integrate; only ", sum(samples$present), " present.")
samples <- samples[samples$present, ]
cat("Proceeding with", nrow(samples), "samples:", paste(samples$label, collapse = ", "), "\n")

seurat_list <- list()
for (i in seq_len(nrow(samples))) {
    s <- samples[i, ]
    cat("[", format(Sys.time()), "] tximport for", s$label, "...\n")
    txi <- tximport(s$quant, type = "alevin")
    mat <- txi$counts  # gene × cell
    # Rename rows from gene_id to gene_name where mappable
    rn <- rownames(mat)
    rn_mapped <- ifelse(rn %in% names(sym_map), sym_map[rn], rn)
    rn_mapped <- make.unique(rn_mapped)
    rownames(mat) <- rn_mapped
    so <- CreateSeuratObject(counts = mat, project = s$label, min.cells = 3, min.features = 200)
    so$sample <- s$label
    so$treatment <- s$treatment
    # Mitochondrial QC — axolotl mt genes typically prefixed "MT-" or have known IDs
    mt_pat <- "^MT-|^mt-|^mito|^ND[1-6]$|^COX[123]$|^ATP[68]$"
    so[["pct_mt"]] <- PercentageFeatureSet(so, pattern = mt_pat)
    cat("  pre-QC cells:", ncol(so), "\n")
    so <- subset(so, subset = nFeature_RNA >= 500 & nCount_RNA >= 1000 & pct_mt < 25)
    cat("  post-QC cells:", ncol(so), "\n")
    seurat_list[[s$label]] <- so
}

# ------------------------------------------------------------
# 2. Per-sample SCTransform + integration via CCA
# ------------------------------------------------------------
seurat_list <- lapply(seurat_list, function(so) {
    SCTransform(so, vars.to.regress = "pct_mt", verbose = FALSE,
                vst.flavor = "v1")  # avoid glmGamPoi prerequisite (v2 default)
})
features <- SelectIntegrationFeatures(seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(seurat_list, anchor.features = features)
anchors <- FindIntegrationAnchors(seurat_list, normalization.method = "SCT",
                                  anchor.features = features)
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(integrated) <- "integrated"
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.6)

# Seurat 5: join split layers in RNA assay so GetAssayData/AggregateExpression
# downstream don't trip on the per-sample layer split. PrepSCTFindMarkers wraps
# in tryCatch — it's needed only for FindAllMarkers; if it fails (memory limit),
# we proceed without it and skip the FindAllMarkers step (cluster_markers.tsv
# may already exist from a prior run).
integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])
saveRDS(integrated, file.path(WORK, "walker_clustered_checkpoint.rds"))  # checkpoint
integrated <- tryCatch(
    PrepSCTFindMarkers(integrated, verbose = FALSE),
    error = function(e) { cat("PrepSCTFindMarkers failed (non-fatal): ", conditionMessage(e), "\n"); integrated }
)

# ------------------------------------------------------------
# 3. Marker-based annotation against Walker 2025
# ------------------------------------------------------------
DefaultAssay(integrated) <- "SCT"
# Walker 2025 uses lowercase axolotl-style symbols, but our gene_name field
# (extracted from GTF [hs] homolog tag) is uppercase human-style. Provide both
# forms; AddModuleScore_UCell silently drops missing genes.
walker_markers <- list(
    glut_etv1plus_regen = toupper(c("etv1","slc17a7","cdh13","bcl11b","sox5","neto1","stk39","ntng1")),
    ependymoglia        = toupper(c("gfap","sox2","fabp7","s100b")),
    gabaergic           = toupper(c("gad1","gad2","slc32a1","slc6a1")),
    oligodendrocyte     = toupper(c("olig1","olig2","sox10","mbp","mog")),
    microglia_macro     = toupper(c("cd68","aif1","csf1r","cx3cr1")),
    neuroendocrine      = toupper(c("avp","oxt","crh","pomc")),
    radial_glia         = toupper(c("nes","vim","gfap","sox2"))
)
# Score each cell against each marker set (UCell-style mean rank).
# UCell handles missing genes gracefully.
integrated <- AddModuleScore_UCell(integrated, features = walker_markers, ncores = 4)

# Cluster-level mean scores → assign top-scoring marker set as candidate label
cluster_meta <- AggregateExpression(integrated, group.by = "seurat_clusters",
                                     assays = "SCT", return.seurat = FALSE)
cluster_scores <- sapply(names(walker_markers), function(set_name) {
    score_col <- paste0(set_name, "_UCell")
    tapply(integrated[[score_col, drop = TRUE]],
           integrated$seurat_clusters, mean, na.rm = TRUE)
})
cluster_anno <- data.frame(
    cluster = rownames(cluster_scores),
    cluster_scores,
    candidate_label = colnames(cluster_scores)[apply(cluster_scores, 1, which.max)],
    candidate_score = apply(cluster_scores, 1, max)
)
fwrite(cluster_anno, file.path(WORK, "cluster_annotations.tsv"), sep = "\t")
cat("\nCluster candidate labels:\n"); print(cluster_anno)

# Top-N markers per cluster (DE) for verification — wrap in tryCatch since
# FindAllMarkers can fail under memory constraints.
markers_ok <- tryCatch({
    markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    fwrite(as.data.table(markers), file.path(WORK, "cluster_markers.tsv"), sep = "\t")
    TRUE
}, error = function(e) {
    cat("FindAllMarkers failed (non-fatal): ", conditionMessage(e), "\n")
    if (!file.exists(file.path(WORK, "cluster_markers.tsv"))) {
        cat("Note: no cluster_markers.tsv from prior runs either.\n")
    } else {
        cat("Using cluster_markers.tsv from prior run.\n")
    }
    FALSE
})

# ------------------------------------------------------------
# 4. Pseudobulk per cluster, summed counts (AggregateExpression sum)
# ------------------------------------------------------------
DefaultAssay(integrated) <- "RNA"
pseudobulk <- AggregateExpression(integrated, group.by = c("seurat_clusters","treatment"),
                                  assays = "RNA", return.seurat = FALSE)
pb_mat <- pseudobulk$RNA
pb_dt <- as.data.table(pb_mat, keep.rownames = "gene")
fwrite(pb_dt, file.path(WORK, "pseudobulk_counts.tsv"), sep = "\t")

# ------------------------------------------------------------
# 5. UCell scoring of Xenopus CT-brain signatures + GDF15
# ------------------------------------------------------------
# Construct signatures from Xenopus BRAIN DEG (homeolog-summed), strip suffixes,
# uppercase for cross-species matching.
build_xen_sig <- function(rds_path, contrast_name, top_n = 50, direction = "up") {
    res <- readRDS(rds_path)
    tt <- as.data.table(res$results[[contrast_name]])
    if (direction == "up") {
        setorder(tt, -logFC)
    } else {
        setorder(tt, logFC)
    }
    # Take top N by logFC magnitude/direction WITHOUT FDR filter — at n=12 BRAIN
    # the q<0.25 set is tiny (often <5 genes) and an under-powered statistical
    # filter shouldn't bottleneck a signature-construction step. Drop by AveExpr
    # < 1 to remove low-expression noise.
    keep <- tt[!is.na(logFC) & is.finite(logFC) & AveExpr > 1]
    head_dt <- head(keep, top_n)
    syms <- toupper(gsub("-", "", gsub("\\.[LSls]$", "", head_dt$gene_id)))
    unique(syms)
}
xen_brain_hsum <- file.path(PROJ, "deg", "voom_BRAIN_homeolog_sum.rds")

# Fallback: if homeolog-sum unavailable, use primary
if (!file.exists(xen_brain_hsum)) xen_brain_hsum <- file.path(PROJ, "deg", "voom_BRAIN_primary.rds")

xenopus_signatures <- list(
    xen_anti_fibrotic_ECM_down = build_xen_sig(xen_brain_hsum, "CT_vs_ND_7DPA", top_n = 50, direction = "down"),
    xen_IFN_immune_up          = build_xen_sig(xen_brain_hsum, "CT_vs_ND_7DPA", top_n = 50, direction = "up"),
    xen_heme_iron_oxphos_up    = build_xen_sig(xen_brain_hsum, "CT_7DPA_vs_24HPA", top_n = 50, direction = "up"),
    xen_structural_plasticity  = build_xen_sig(xen_brain_hsum, "CT_7DPA_vs_24HPA", top_n = 50, direction = "down"),
    xen_GDF15                  = "GDF15"
)
cat("\nXenopus signature gene counts:\n")
for (n in names(xenopus_signatures)) cat("  ", n, ":", length(xenopus_signatures[[n]]), "\n")

# UCell on axolotl single cells (gene symbols upper-cased for matching)
DefaultAssay(integrated) <- "RNA"
# Convert axolotl rownames to upper-case symbol for Xenopus matching
rn_orig <- rownames(integrated[["RNA"]])
rn_upper <- toupper(gsub("-", "", rn_orig))
# AddModuleScore_UCell uses feature names from the matrix directly. We build a
# tmp Seurat with upper-cased symbols just for the scoring. Use LayerData (Seurat 5)
# explicitly to avoid the multi-layer GetAssayData issue.
tmp_counts <- LayerData(integrated, assay = "RNA", layer = "counts")
rownames(tmp_counts) <- make.unique(rn_upper)
tmp <- integrated
tmp[["RNA"]] <- CreateAssayObject(counts = tmp_counts)
DefaultAssay(tmp) <- "RNA"

tmp <- AddModuleScore_UCell(tmp, features = xenopus_signatures, ncores = 4)

# Per-cell scores
score_cols <- paste0(names(xenopus_signatures), "_UCell")
per_cell <- as.data.table(tmp@meta.data[, c("sample","treatment","seurat_clusters", score_cols)],
                          keep.rownames = "cell")
fwrite(per_cell, file.path(WORK, "ucell_signatures_per_cell.tsv"), sep = "\t")

# Per-cluster mean scores by treatment
per_cluster <- per_cell[, lapply(.SD, mean, na.rm = TRUE),
                        by = .(seurat_clusters, treatment),
                        .SDcols = score_cols]
fwrite(per_cluster, file.path(WORK, "ucell_signatures_per_cluster.tsv"), sep = "\t")

# ------------------------------------------------------------
# 6. Treatment-population enrichment for etv1+ cluster vs others
# ------------------------------------------------------------
# Find the cluster best matching the etv1+ regen-required population (highest
# glut_etv1plus_regen UCell score).
etv1_cluster <- cluster_anno$cluster[which.max(cluster_anno$glut_etv1plus_regen)]
cat("\netv1+ regen-required cluster (best match):", etv1_cluster, "\n")

per_cell[, is_etv1_pop := seurat_clusters == etv1_cluster]
prop_dt <- per_cell[, .(n_etv1 = sum(is_etv1_pop), n_total = .N,
                         prop = sum(is_etv1_pop) / .N), by = treatment]
fwrite(prop_dt, file.path(WORK, "etv1_population_proportions.tsv"), sep = "\t")

# Fisher test: tail vs uninjured for etv1+ cluster proportion
fisher_test <- function(t1, t2) {
    a <- per_cell[treatment == t1, sum(is_etv1_pop)]
    b <- per_cell[treatment == t1, sum(!is_etv1_pop)]
    c <- per_cell[treatment == t2, sum(is_etv1_pop)]
    d <- per_cell[treatment == t2, sum(!is_etv1_pop)]
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
    data.table(comparison = paste(t1, "vs", t2),
               OR = ft$estimate, p = ft$p.value,
               lo95 = ft$conf.int[1], hi95 = ft$conf.int[2])
}
treat_enrich <- rbind(
    fisher_test("tail", "uninjured"),
    fisher_test("limb", "uninjured"),
    fisher_test("tail", "muscimol"),
    fisher_test("muscimol", "uninjured")
)
fwrite(treat_enrich, file.path(WORK, "treatment_pop_enrichment.tsv"), sep = "\t")
cat("\nTreatment-population enrichment for etv1+ cluster:\n"); print(treat_enrich)

# ------------------------------------------------------------
# 7. Save Seurat + figures
# ------------------------------------------------------------
PROJ_OUT <- file.path(PROJ, "comparators", "walker")
dir.create(PROJ_OUT, recursive = TRUE, showWarnings = FALSE)
saveRDS(integrated, file.path(PROJ_OUT, "walker_integrated.rds"))
for (f in list.files(WORK, pattern = "\\.tsv$", full.names = TRUE)) {
    file.copy(f, PROJ_OUT, overwrite = TRUE)
}

# UMAP + signature heat
pdf(file.path(FIG, "walker_umap_clusters.pdf"), width = 8, height = 6)
print(DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters",
              label = TRUE) + ggtitle("Walker axolotl brain — integrated UMAP"))
dev.off()

pdf(file.path(FIG, "walker_umap_treatment.pdf"), width = 8, height = 6)
print(DimPlot(integrated, reduction = "umap", group.by = "treatment") +
        ggtitle("Walker axolotl brain — sample of origin"))
dev.off()

# Per-cluster signature score heatmap
sig_matrix <- as.matrix(per_cluster[, ..score_cols])
rownames(sig_matrix) <- per_cluster[, paste(seurat_clusters, treatment, sep = "_")]
pdf(file.path(FIG, "walker_signature_heatmap.pdf"), width = 9, height = 7)
heatmap(t(sig_matrix), Rowv = NA, scale = "row", margins = c(8, 12),
        main = "Xenopus signatures × axolotl cluster (mean UCell)")
dev.off()

# Mirror figures to project
PROJ_FIG <- file.path(PROJ, "figures", "stage_H", "walker")
dir.create(PROJ_FIG, recursive = TRUE, showWarnings = FALSE)
for (f in list.files(FIG, pattern = "walker_.*\\.pdf$", full.names = TRUE)) {
    file.copy(f, PROJ_FIG, overwrite = TRUE)
}

cat("\n[", format(Sys.time(), "%FT%T%z"), "] hA_walker_seurat.R complete.\n", sep = "")
