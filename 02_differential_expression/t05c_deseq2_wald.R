#!/usr/bin/env Rscript
# t05c_deseq2_wald.R — DESeq2 Wald PRIMARY differential expression.
#
# This now supersedes voom_*_primary.rds as the brain DE primary (and is run
# alongside it for BLAST for two-engine cross-validation). limma-voom outputs
# are retained as confirmatory; downstream pipeline should preferentially read
# deseq2_*_primary.rds for FDR/logFC.
#
# Outputs:
#   deg/deseq2_BRAIN_primary.rds          full dds + per-contrast results list
#   deg/deseq2_BLAST_primary.rds
#   deg/deseq2_BRAIN_sensitivity.rds      brain w/ CT_T3_7DPA_BRAIN dropped (SO-04)
#   deg/deseq2_*_primary_summary.tsv      sig_q05 / sig_q25 per contrast
#   deg/deseq2_*_primary.xlsx             one sheet per contrast + summary
#   deg/all_deseq2_wald_summary.tsv       combined summary

suppressPackageStartupMessages({
    library(DESeq2); library(edgeR); library(openxlsx)
})

PROJ    <- "/project/def-srussell/nirosha_xeno"
DEG_DIR <- file.path(PROJ, "deg")
COUNTS  <- file.path(PROJ, "counts", "star_gene_counts.tsv")
META    <- file.path(DEG_DIR, "sample_metadata.rds")

cat("[", format(Sys.time(), "%FT%T%z"), "] t05c DESeq2 Wald primary start\n")

mat <- as.matrix(read.table(COUNTS, header = TRUE, sep = "\t", row.names = 1,
                             check.names = FALSE, stringsAsFactors = FALSE))
storage.mode(mat) <- "integer"
meta <- readRDS(META)
stopifnot(all(colnames(mat) == meta$sample))

run_deseq2_wald <- function(mat_t, meta_t, contrasts, label, out_prefix) {
    # Combined-group design lets us use arbitrary pairwise contrasts uniformly
    meta_t$grp <- factor(paste(meta_t$treatment, meta_t$time, sep = "_"))
    cat("\n  --- ", label, " ---\n", sep = "")
    cat("  Group factor levels: ", paste(levels(meta_t$grp), collapse = ", "), "\n")

    # CPM-based prefilter to match the limma/voom filter (CPM > 1 in >= 3 samples)
    keep <- rowSums(edgeR::cpm(mat_t) > 1) >= 3
    mat_t <- mat_t[keep, , drop = FALSE]
    cat("  Genes after CPM>1 in >=3 filter: ", nrow(mat_t), "\n")

    dds <- DESeqDataSetFromMatrix(
        countData = mat_t,
        colData   = meta_t,
        design    = ~ 0 + grp
    )
    dds <- DESeq(dds)
    cat("  DESeq2 fit complete; sizeFactors range: [",
        round(min(sizeFactors(dds)), 3), ", ",
        round(max(sizeFactors(dds)), 3), "]\n", sep = "")

    results_list <- list()
    for (cn in names(contrasts)) {
        lev <- contrasts[[cn]]
        if (!(lev[1] %in% levels(meta_t$grp)) ||
            !(lev[2] %in% levels(meta_t$grp))) {
            cat("  [skip] ", cn, ": missing level\n"); next
        }
        # Use group factor with name="grp", numerator = lev[1], denominator = lev[2]
        res <- results(dds, contrast = c("grp", lev[1], lev[2]),
                       independentFiltering = TRUE, alpha = 0.05)
        df <- as.data.frame(res)
        df$gene_id  <- rownames(df)
        df$contrast <- cn
        df <- df[order(df$pvalue, na.last = TRUE), ]
        # Normalise column names so downstream code can read uniformly
        df$logFC     <- df$log2FoldChange
        df$P.Value   <- df$pvalue
        df$adj.P.Val <- df$padj
        results_list[[cn]] <- df
    }

    sig_summary <- do.call(rbind, lapply(names(results_list), function(cn) {
        df <- results_list[[cn]]
        data.frame(
            contrast = cn,
            n_total  = sum(!is.na(df$padj)),
            sig_q05  = sum(df$padj < 0.05, na.rm = TRUE),
            sig_q25  = sum(df$padj < 0.25, na.rm = TRUE),
            up_q05   = sum(df$padj < 0.05 & df$logFC > 0, na.rm = TRUE),
            down_q05 = sum(df$padj < 0.05 & df$logFC < 0, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    }))
    cat("  Sig summary:\n"); print(sig_summary)

    saveRDS(list(dds = dds, results = results_list, summary = sig_summary,
                 contrasts = contrasts, design_label = "deseq2_wald_grp_combined"),
            file.path(DEG_DIR, paste0(out_prefix, ".rds")))
    write.table(sig_summary,
                file.path(DEG_DIR, paste0(out_prefix, "_summary.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    wb <- createWorkbook()
    for (cn in names(results_list)) {
        addWorksheet(wb, cn)
        writeData(wb, cn, results_list[[cn]])
    }
    addWorksheet(wb, "summary"); writeData(wb, "summary", sig_summary)
    saveWorkbook(wb, file.path(DEG_DIR, paste0(out_prefix, ".xlsx")),
                 overwrite = TRUE)

    sig_summary$run <- out_prefix
    sig_summary
}

# ---- BRAIN: two timepoints (24HPA, 7DPA) ----
brain_idx <- meta$tissue == "BRAIN"
brain_contrasts <- list(
    CT_vs_ND_24HPA      = c("CT_24HPA", "ND_24HPA"),
    CT_vs_ND_7DPA       = c("CT_7DPA",  "ND_7DPA"),
    CT_7DPA_vs_24HPA    = c("CT_7DPA",  "CT_24HPA"),
    ND_7DPA_vs_24HPA    = c("ND_7DPA",  "ND_24HPA")
)
brain_meta <- droplevels(meta[brain_idx, ])
brain_summary <- run_deseq2_wald(
    as.matrix(mat[, brain_idx]), brain_meta, brain_contrasts,
    "BRAIN primary", "deseq2_BRAIN_primary"
)

# Sensitivity: drop CT_T3_7DPA_BRAIN
sens_keep <- brain_meta$sample != "CT_T3_7DPA_BRAIN_S27"
sens_summary <- run_deseq2_wald(
    as.matrix(mat[, brain_idx])[, sens_keep],
    droplevels(brain_meta[sens_keep, ]),
    brain_contrasts, "BRAIN sensitivity (drop CT_T3_7DPA)",
    "deseq2_BRAIN_sensitivity"
)

# ---- BLAST: three timepoints (11HPA, 24HPA, 7DPA) ----
blast_idx <- meta$tissue == "BLAST"
blast_contrasts <- list(
    CT_vs_ND_11HPA      = c("CT_11HPA", "ND_11HPA"),
    CT_vs_ND_24HPA      = c("CT_24HPA", "ND_24HPA"),
    CT_vs_ND_7DPA       = c("CT_7DPA",  "ND_7DPA"),
    CT_24HPA_vs_11HPA   = c("CT_24HPA", "CT_11HPA"),
    ND_24HPA_vs_11HPA   = c("ND_24HPA", "ND_11HPA"),
    CT_7DPA_vs_24HPA    = c("CT_7DPA",  "CT_24HPA"),
    ND_7DPA_vs_24HPA    = c("ND_7DPA",  "ND_24HPA"),
    CT_7DPA_vs_11HPA    = c("CT_7DPA",  "CT_11HPA"),
    ND_7DPA_vs_11HPA    = c("ND_7DPA",  "ND_11HPA")
)
blast_summary <- run_deseq2_wald(
    as.matrix(mat[, blast_idx]), droplevels(meta[blast_idx, ]),
    blast_contrasts, "BLAST primary", "deseq2_BLAST_primary"
)

# ---- combined roll-up ----
combined <- do.call(rbind, list(brain_summary, sens_summary, blast_summary))
write.table(combined, file.path(DEG_DIR, "all_deseq2_wald_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== ALL DESeq2 Wald summary ===\n")
print(combined)

cat("\n[", format(Sys.time(), "%FT%T%z"), "] t05c done\n")
