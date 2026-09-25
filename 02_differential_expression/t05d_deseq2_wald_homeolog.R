#!/usr/bin/env Rscript
# t05d_deseq2_wald_homeolog.R — DESeq2 Wald on homeolog-summed counts.
# Mirrors t05c but on parent-symbol aggregation (matches t10's logic).
#
# Outputs:
#   deg/deseq2_BRAIN_primary_homeolog_sum.{rds,xlsx,_summary.tsv}
#   deg/deseq2_BLAST_primary_homeolog_sum.{rds,xlsx,_summary.tsv}
#   deg/deseq2_BRAIN_sensitivity_homeolog_sum.{rds,xlsx,_summary.tsv}
#   deg/all_deseq2_wald_homeolog_summary.tsv

suppressPackageStartupMessages({
    library(DESeq2); library(edgeR); library(openxlsx); library(data.table)
})

PROJ    <- "/project/def-srussell/nirosha_xeno"
DEG_DIR <- file.path(PROJ, "deg")
COUNTS  <- file.path(PROJ, "counts", "star_gene_counts.tsv")
META    <- file.path(DEG_DIR, "sample_metadata.rds")

cat("[", format(Sys.time(), "%FT%T%z"), "] t05d DESeq2 Wald homeolog-sum start\n")

ct <- fread(COUNTS, header = TRUE)
gene_id <- ct$gene_id
mat_full <- as.matrix(ct[, -"gene_id"])
storage.mode(mat_full) <- "integer"
rownames(mat_full) <- gene_id
meta <- readRDS(META)
stopifnot(all(colnames(mat_full) == meta$sample))

# Homeolog parent collapse — same rule as t10 (strip trailing .L/.S/.l/.s)
parent_symbol <- gsub("\\.[LSls]$", "", gene_id)
mat_h <- rowsum(mat_full, group = parent_symbol, reorder = FALSE)
storage.mode(mat_h) <- "integer"
cat("  Homeolog-summed matrix:", nrow(mat_h), "parents x", ncol(mat_h), "samples\n")

run_deseq2_wald <- function(mat_t, meta_t, contrasts, label, out_prefix) {
    meta_t$grp <- factor(paste(meta_t$treatment, meta_t$time, sep = "_"))
    cat("\n  --- ", label, " ---\n", sep = "")
    keep <- rowSums(edgeR::cpm(mat_t) > 1) >= 3
    mat_t <- mat_t[keep, , drop = FALSE]
    cat("  After CPM>1 in >=3 filter: ", nrow(mat_t), " parents\n")

    dds <- DESeqDataSetFromMatrix(countData = mat_t, colData = meta_t,
                                   design = ~ 0 + grp)
    dds <- DESeq(dds)

    results_list <- list()
    for (cn in names(contrasts)) {
        lev <- contrasts[[cn]]
        if (!(lev[1] %in% levels(meta_t$grp)) ||
            !(lev[2] %in% levels(meta_t$grp))) next
        res <- results(dds, contrast = c("grp", lev[1], lev[2]),
                       independentFiltering = TRUE, alpha = 0.05)
        df <- as.data.frame(res)
        df$gene_id  <- rownames(df)
        df$contrast <- cn
        df$logFC     <- df$log2FoldChange
        df$P.Value   <- df$pvalue
        df$adj.P.Val <- df$padj
        df <- df[order(df$pvalue, na.last = TRUE), ]
        results_list[[cn]] <- df
    }

    sig <- do.call(rbind, lapply(names(results_list), function(cn) {
        df <- results_list[[cn]]
        data.frame(contrast = cn,
                   n_total  = sum(!is.na(df$padj)),
                   sig_q05  = sum(df$padj < 0.05, na.rm = TRUE),
                   sig_q25  = sum(df$padj < 0.25, na.rm = TRUE),
                   up_q05   = sum(df$padj < 0.05 & df$logFC > 0, na.rm = TRUE),
                   down_q05 = sum(df$padj < 0.05 & df$logFC < 0, na.rm = TRUE),
                   stringsAsFactors = FALSE)
    }))
    print(sig)

    saveRDS(list(dds = dds, results = results_list, summary = sig,
                 contrasts = contrasts,
                 design_label = "deseq2_wald_grp_homeolog_sum"),
            file.path(DEG_DIR, paste0(out_prefix, ".rds")))
    write.table(sig, file.path(DEG_DIR, paste0(out_prefix, "_summary.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    wb <- createWorkbook()
    for (cn in names(results_list)) {
        addWorksheet(wb, cn); writeData(wb, cn, results_list[[cn]])
    }
    addWorksheet(wb, "summary"); writeData(wb, "summary", sig)
    saveWorkbook(wb, file.path(DEG_DIR, paste0(out_prefix, ".xlsx")), overwrite = TRUE)
    sig$run <- out_prefix
    sig
}

brain_idx <- meta$tissue == "BRAIN"
brain_contrasts <- list(
    CT_vs_ND_24HPA   = c("CT_24HPA", "ND_24HPA"),
    CT_vs_ND_7DPA    = c("CT_7DPA",  "ND_7DPA"),
    CT_7DPA_vs_24HPA = c("CT_7DPA",  "CT_24HPA"),
    ND_7DPA_vs_24HPA = c("ND_7DPA",  "ND_24HPA")
)
brain_meta <- droplevels(meta[brain_idx, ])
b <- run_deseq2_wald(mat_h[, brain_idx], brain_meta, brain_contrasts,
                     "BRAIN homeolog-sum primary",
                     "deseq2_BRAIN_primary_homeolog_sum")

sens_keep <- brain_meta$sample != "CT_T3_7DPA_BRAIN_S27"
s <- run_deseq2_wald(mat_h[, brain_idx][, sens_keep],
                     droplevels(brain_meta[sens_keep, ]),
                     brain_contrasts, "BRAIN homeolog-sum sensitivity",
                     "deseq2_BRAIN_sensitivity_homeolog_sum")

blast_idx <- meta$tissue == "BLAST"
blast_contrasts <- list(
    CT_vs_ND_11HPA    = c("CT_11HPA", "ND_11HPA"),
    CT_vs_ND_24HPA    = c("CT_24HPA", "ND_24HPA"),
    CT_vs_ND_7DPA     = c("CT_7DPA",  "ND_7DPA"),
    CT_24HPA_vs_11HPA = c("CT_24HPA", "CT_11HPA"),
    ND_24HPA_vs_11HPA = c("ND_24HPA", "ND_11HPA"),
    CT_7DPA_vs_24HPA  = c("CT_7DPA",  "CT_24HPA"),
    ND_7DPA_vs_24HPA  = c("ND_7DPA",  "ND_24HPA"),
    CT_7DPA_vs_11HPA  = c("CT_7DPA",  "CT_11HPA"),
    ND_7DPA_vs_11HPA  = c("ND_7DPA",  "ND_11HPA")
)
l <- run_deseq2_wald(mat_h[, blast_idx], droplevels(meta[blast_idx, ]),
                     blast_contrasts, "BLAST homeolog-sum primary",
                     "deseq2_BLAST_primary_homeolog_sum")

combined <- do.call(rbind, list(b, s, l))
write.table(combined, file.path(DEG_DIR, "all_deseq2_wald_homeolog_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\n=== ALL DESeq2 Wald homeolog summary ===\n"); print(combined)
cat("\n[", format(Sys.time(), "%FT%T%z"), "] t05d done\n")
