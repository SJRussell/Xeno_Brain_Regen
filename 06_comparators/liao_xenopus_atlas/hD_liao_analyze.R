#!/usr/bin/env Rscript
# H-D comparator: Liao 2022 XCL adult Xenopus brain atlas (GSE195790, Figshare 19152839).
# 23,906 brain cells across 4 Brain* samples × ~50 celltype labels.
# Goal: score each Liao brain celltype's marker signature against our 12 Xenopus
# brain bulk samples; test treatment × time interaction on per-celltype scores.
# Same-species matching (Xenopus laevis L/S subgenome symbols, no ortholog mapping).
#
# Output dir: /scratch/.../nirosha_xeno_comparators/liao/work/
#   per_celltype_markers.tsv         (top 50 markers per celltype, log-pseudobulk-FC)
#   per_celltype_scores.tsv          (UCell scores per celltype × per Xenopus brain sample)
#   per_celltype_interaction.tsv     (treatment × time interaction lm per celltype)
#   figures/{celltype_score_violin, sig_celltype_heatmap}.pdf
# Mirrored to /project/.../comparators/liao/

suppressPackageStartupMessages({
    library(Matrix)
    library(data.table)
    library(SingleCellExperiment)
    library(scran)
    library(UCell)
    library(DESeq2)
    library(SummarizedExperiment)
    library(ggplot2)
})

PROJ <- "/project/def-srussell/nirosha_xeno"
DEG_DIR <- file.path(PROJ, "deg")
LIAO_RAW <- "/scratch/srussell/srussell/nirosha_xeno_comparators/liao/raw"
WORK <- "/scratch/srussell/srussell/nirosha_xeno_comparators/liao/work"
FIG  <- file.path(WORK, "figures")
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG,  recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Cell metadata + cellID → celltype map
# ------------------------------------------------------------
cat("[", format(Sys.time()), "] Loading cell_info...\n", sep = "")
cell_info <- fread(file.path(LIAO_RAW, "cell_info/dge_cell_info/Brain_cell_info.csv"))
setnames(cell_info, "V1", "rowname")
# cellID format: "BrainN_X.YY.<barcode>"; matches dge column names "<barcode>" prefixed by file name
cell_info[, sample_file := sub("\\.([ACGT]+)$", "", cellID)]  # e.g., Brain1_3.25.CTGT... → Brain1_3.25
cat("Brain cells annotated:", nrow(cell_info), "\n")
cat("Top celltypes:\n"); print(head(cell_info[, .N, by = celltype][order(-N)], 25))

# ------------------------------------------------------------
# 2. Read each Brain*_dge.txt.gz — extract from zip in-memory
# ------------------------------------------------------------
brain_files <- c("Brain1_03.25_dge.txt.gz", "Brain2_03.25_dge.txt.gz",
                 "Brain3_04.25_dge.txt.gz", "Brain4_04.25_dge.txt.gz")
EXTRACT_DIR <- file.path(LIAO_RAW, "brain_extracted")
dir.create(EXTRACT_DIR, showWarnings = FALSE)

read_brain_matrix <- function(fname) {
    out_path <- file.path(EXTRACT_DIR, fname)
    if (!file.exists(out_path)) {
        # Extract from zip
        zip_path <- file.path(LIAO_RAW, "dge_raw_data.zip")
        cmd <- sprintf("cd %s && unzip -p %s dge_raw_data/%s > %s",
                       shQuote(EXTRACT_DIR), shQuote(zip_path),
                       fname, shQuote(out_path))
        system(cmd)
    }
    cat("Reading", fname, "...\n")
    dt <- fread(cmd = paste("zcat", shQuote(out_path)), header = TRUE)
    gene_col <- colnames(dt)[1]  # "GENE"
    genes <- dt[[gene_col]]
    m <- as.matrix(dt[, -1, with = FALSE])
    storage.mode(m) <- "integer"
    rownames(m) <- genes
    Matrix(m, sparse = TRUE)
}

mat_list <- lapply(brain_files, read_brain_matrix)
# Cell column names need to be prefixed by sample-file ID to match cell_info$cellID
sample_ids <- sub("_dge\\.txt\\.gz$", "", brain_files)  # Brain1_03.25 etc.
# But cell_info uses "Brain1_3.25" (no leading 0). Normalize.
sample_ids_norm <- gsub("_0", "_", sample_ids)
for (i in seq_along(mat_list)) {
    colnames(mat_list[[i]]) <- paste0(sample_ids_norm[i], ".", colnames(mat_list[[i]]))
}

# Align all matrices on union of genes
all_genes <- Reduce(union, lapply(mat_list, rownames))
align_mat <- function(m) {
    out <- Matrix(0, nrow = length(all_genes), ncol = ncol(m), sparse = TRUE)
    rownames(out) <- all_genes
    colnames(out) <- colnames(m)
    out[rownames(m), ] <- m
    out
}
mat_list_aligned <- lapply(mat_list, align_mat)
big_mat <- do.call(cbind, mat_list_aligned)
cat("Combined Brain matrix:", nrow(big_mat), "genes ×", ncol(big_mat), "cells\n")

# ------------------------------------------------------------
# 3. Match cells to annotations
# ------------------------------------------------------------
keep <- colnames(big_mat) %in% cell_info$cellID
big_mat <- big_mat[, keep]
cat("After annotation match:", ncol(big_mat), "cells\n")

# Order cell_info by big_mat colnames
cell_info_ord <- cell_info[match(colnames(big_mat), cellID)]
stopifnot(identical(cell_info_ord$cellID, colnames(big_mat)))

# Drop celltypes with <30 cells (signature too noisy)
ct_n <- cell_info_ord[, .N, by = celltype][N >= 30]
keep_celltype <- cell_info_ord$celltype %in% ct_n$celltype
big_mat <- big_mat[, keep_celltype]
cell_info_ord <- cell_info_ord[keep_celltype]
cat("After ≥30-cell filter:", ncol(big_mat), "cells across",
    length(unique(cell_info_ord$celltype)), "celltypes\n")

# ------------------------------------------------------------
# 4. Per-celltype top markers via simple log-fold-change of CPM-normalized pseudobulk
# ------------------------------------------------------------
cat("[", format(Sys.time()), "] Computing per-celltype markers...\n", sep = "")
# Pseudobulk per celltype
ct_levels <- unique(cell_info_ord$celltype)
pb <- sapply(ct_levels, function(ct) {
    sel <- cell_info_ord$celltype == ct
    Matrix::rowSums(big_mat[, sel, drop = FALSE])
})
colnames(pb) <- ct_levels
# Library-size normalize per celltype (CPM)
pb_cpm <- sweep(pb, 2, colSums(pb) / 1e6, "/")
# log2-CPM (with pseudocount)
pb_log <- log2(pb_cpm + 1)

# Markers: log fold change of celltype vs mean of other celltypes
markers_list <- list()
for (ct in ct_levels) {
    other_mean <- rowMeans(pb_log[, setdiff(colnames(pb_log), ct), drop = FALSE])
    fc <- pb_log[, ct] - other_mean
    top <- names(sort(fc, decreasing = TRUE))[1:50]
    markers_list[[ct]] <- top
}
markers_dt <- rbindlist(lapply(names(markers_list), function(ct) {
    data.table(celltype = ct, gene = markers_list[[ct]],
               rank = seq_along(markers_list[[ct]]))
}))
fwrite(markers_dt, file.path(WORK, "per_celltype_markers.tsv"), sep = "\t")

# ------------------------------------------------------------
# 5. Score celltype signatures against Xenopus brain bulk via UCell-on-bulk
# ------------------------------------------------------------
cat("[", format(Sys.time()), "] Loading Xenopus brain bulk + scoring...\n", sep = "")
vsd  <- readRDS(file.path(DEG_DIR, "vsd_blind.rds"))
meta <- readRDS(file.path(DEG_DIR, "sample_metadata.rds"))
brain_idx <- which(meta$tissue == "BRAIN")
brain_meta <- meta[brain_idx, , drop = FALSE]
brain_mat  <- assay(vsd)[, brain_idx, drop = FALSE]
cat("Xenopus BRAIN bulk:", dim(brain_mat)[1], "genes ×", dim(brain_mat)[2], "samples\n")
cat("Liao gene overlap with Xenopus bulk:",
    length(intersect(rownames(brain_mat), markers_dt$gene)), "of",
    length(unique(markers_dt$gene)), "marker genes\n")

# UCell expects features × cells. For bulk, treat each sample as a "cell".
ucell_scores <- ScoreSignatures_UCell(
    brain_mat,
    features = markers_list,
    name = "",
    ncores = 4
)
ucell_scores <- as.data.table(ucell_scores, keep.rownames = "sample")
fwrite(ucell_scores, file.path(WORK, "per_celltype_scores.tsv"), sep = "\t")

# ------------------------------------------------------------
# 6. Per-celltype treatment × time interaction (only 24HPA + 7DPA — paired-tissue design)
# ------------------------------------------------------------
cat("[", format(Sys.time()), "] Per-celltype interaction tests...\n", sep = "")
brain_meta$pair_id <- paste(brain_meta$rep, brain_meta$time, brain_meta$treatment, sep = "_")
# UCell scores aligned by sample name
score_long <- melt(ucell_scores, id.vars = "sample",
                   variable.name = "celltype", value.name = "score")
score_long <- merge(score_long,
                    data.table(sample = brain_meta$sample,
                               rep = brain_meta$rep,
                               time = brain_meta$time,
                               treatment = brain_meta$treatment),
                    by = "sample")

interaction_results <- score_long[, {
    if (length(unique(treatment)) < 2 || length(unique(time)) < 2) {
        list(p_treat = NA_real_, p_time = NA_real_, p_int = NA_real_,
             est_treat = NA_real_, est_time = NA_real_, est_int = NA_real_, n = .N)
    } else {
        df <- data.frame(score = score, rep = rep,
                         treatment = factor(treatment, c("ND", "CT")),
                         time      = factor(time,      c("24HPA", "7DPA")))
        fit <- tryCatch(lm(score ~ treatment * time + rep, data = df),
                        error = function(e) NULL)
        if (is.null(fit)) {
            list(p_treat = NA_real_, p_time = NA_real_, p_int = NA_real_,
                 est_treat = NA_real_, est_time = NA_real_, est_int = NA_real_, n = .N)
        } else {
            cf <- summary(fit)$coefficients
            get <- function(term, col) if (term %in% rownames(cf)) cf[term, col] else NA_real_
            list(p_treat = get("treatmentCT", 4), p_time = get("time7DPA", 4),
                 p_int = get("treatmentCT:time7DPA", 4),
                 est_treat = get("treatmentCT", 1), est_time = get("time7DPA", 1),
                 est_int = get("treatmentCT:time7DPA", 1),
                 n = .N)
        }
    }
}, by = celltype]

interaction_results[, fdr_int := p.adjust(p_int, method = "BH")]
interaction_results[, fdr_treat := p.adjust(p_treat, method = "BH")]
setorder(interaction_results, p_int)
fwrite(interaction_results, file.path(WORK, "per_celltype_interaction.tsv"), sep = "\t")
cat("\n=== Top 10 celltypes by interaction p ===\n")
print(head(interaction_results, 10))
cat("\n=== Significant interaction (FDR<0.10) ===\n")
print(interaction_results[fdr_int < 0.10])

# ------------------------------------------------------------
# 7. Quick figures
# ------------------------------------------------------------
sig_ct <- interaction_results[fdr_int < 0.25, celltype]
if (length(sig_ct) > 0) {
    pdf(file.path(FIG, "celltype_score_violin.pdf"), width = 10, height = 6)
    plot_dt <- score_long[celltype %in% sig_ct]
    plot_dt[, group := paste(treatment, time, sep = "_")]
    print(ggplot(plot_dt, aes(group, score, color = treatment)) +
              geom_violin(alpha = 0.4) + geom_jitter(width = 0.1, alpha = 0.7) +
              facet_wrap(~ celltype, scales = "free_y") +
              theme_minimal() +
              labs(title = "Liao celltype-marker UCell scores in Xenopus brain bulk (FDR_int<0.25)"))
    dev.off()
}

# ------------------------------------------------------------
# 8. Mirror to /project
# ------------------------------------------------------------
PROJ_OUT <- file.path(PROJ, "comparators", "liao")
dir.create(PROJ_OUT, recursive = TRUE, showWarnings = FALSE)
for (f in list.files(WORK, pattern = "\\.tsv$", full.names = TRUE)) {
    file.copy(f, PROJ_OUT, overwrite = TRUE)
}
PROJ_FIG <- file.path(PROJ, "figures", "stage_H", "liao")
dir.create(PROJ_FIG, recursive = TRUE, showWarnings = FALSE)
for (f in list.files(FIG, pattern = "\\.pdf$", full.names = TRUE)) {
    file.copy(f, PROJ_FIG, overwrite = TRUE)
}

cat("[", format(Sys.time()), "] hD_liao_analyze.R complete.\n", sep = "")
