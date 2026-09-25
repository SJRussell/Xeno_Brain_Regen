#!/usr/bin/env Rscript
# Stage D: differential expression
# - D1: limma-voom + voomWithQualityWeights (PRIMARY) per tissue
# - D2: DESeq2 LRT for interaction (confirmatory)
# - Sensitivity D1 with CT_T3_7DPA_BRAIN dropped (per SO-04)
# - Contrasts:
#     - per timepoint: CT vs ND
#     - per treatment: 7DPA vs 24HPA  (and 24 vs 11 in BLAST)
#     - INTERACTION (KEY): (CT_late - CT_early) - (ND_late - ND_early)

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(DESeq2)
    library(openxlsx)
})

PROJ <- "/project/def-srussell/nirosha_xeno"
DEG_DIR <- file.path(PROJ, "deg")
FIG_DIR <- file.path(PROJ, "figures", "stage_D")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("[", format(Sys.time(), "%FT%T%z"), "] Loading filtered DGEList from Stage C...\n")
y <- readRDS(file.path(DEG_DIR, "dgelist_filtered_tmm.rds"))
meta <- readRDS(file.path(DEG_DIR, "sample_metadata.rds"))
stopifnot(all(rownames(y$samples) == meta$sample))
cat("  Genes:", nrow(y), "  Samples:", ncol(y), "\n")

# ============================================================
#  Helpers
# ============================================================

run_voom_de <- function(y_sub, meta_sub, contrasts_named, out_prefix, use_qw = TRUE) {
    # Design with no intercept, group factor = treatment_time
    grp <- factor(paste(meta_sub$treatment, meta_sub$time, sep = "_"))
    design <- model.matrix(~ 0 + grp)
    colnames(design) <- levels(grp)

    if (use_qw) {
        v <- voomWithQualityWeights(y_sub, design = design, plot = FALSE)
    } else {
        v <- voom(y_sub, design = design, plot = FALSE)
    }
    fit <- lmFit(v, design)

    # Build the contrast matrix from named string vector
    contr_obj <- makeContrasts(contrasts = contrasts_named, levels = design)
    colnames(contr_obj) <- names(contrasts_named)
    fit2 <- contrasts.fit(fit, contr_obj)
    fit2 <- eBayes(fit2)

    # topTable per contrast
    res_list <- lapply(colnames(contr_obj), function(cn) {
        tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "P")
        tt$gene_id <- rownames(tt)
        tt$contrast <- cn
        tt[, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "contrast")]
    })
    names(res_list) <- colnames(contr_obj)

    # Summary: significant genes per contrast at FDR<0.05 and FDR<0.25
    sig_summary <- data.frame(
        contrast    = names(res_list),
        n_total     = sapply(res_list, nrow),
        sig_q05     = sapply(res_list, function(d) sum(d$adj.P.Val < 0.05)),
        sig_q25     = sapply(res_list, function(d) sum(d$adj.P.Val < 0.25)),
        up_q05      = sapply(res_list, function(d) sum(d$adj.P.Val < 0.05 & d$logFC > 0)),
        down_q05    = sapply(res_list, function(d) sum(d$adj.P.Val < 0.05 & d$logFC < 0)),
        stringsAsFactors = FALSE
    )

    list(v = v, fit = fit2, results = res_list, summary = sig_summary,
         contrasts = contr_obj, design = design)
}

run_deseq2_lrt <- function(y_sub, meta_sub, full_formula, reduced_formula, contrast_label) {
    # DESeq2 with treatment*time design; LRT compares full vs reduced (drops interaction)
    cd <- DESeqDataSetFromMatrix(
        countData = round(as.matrix(y_sub$counts)),
        colData   = meta_sub,
        design    = full_formula
    )
    cd <- DESeq(cd, test = "LRT", reduced = reduced_formula)
    res <- results(cd)
    df <- as.data.frame(res)
    df$gene_id <- rownames(df)
    df$contrast <- contrast_label
    df <- df[order(df$padj, na.last = TRUE), ]
    list(dds = cd, results = df,
         summary = data.frame(
            contrast = contrast_label,
            n_total  = nrow(df),
            sig_q05  = sum(df$padj < 0.05, na.rm = TRUE),
            sig_q25  = sum(df$padj < 0.25, na.rm = TRUE)
         ))
}

# ============================================================
#  PER-TISSUE DE
# ============================================================

all_summaries <- list()
all_results <- list()

for (tiss in c("BRAIN", "BLAST")) {
    cat("\n==========================================\n")
    cat("Tissue: ", tiss, "\n")
    cat("==========================================\n")

    sub_idx <- meta$tissue == tiss
    y_sub <- y[, sub_idx, keep.lib.sizes = FALSE]
    meta_sub <- meta[sub_idx, ]

    # Refilter for this tissue (CPM > 1 in >= 3 samples within tissue)
    cpm_sub <- cpm(y_sub)
    keep_sub <- rowSums(cpm_sub > 1) >= 3
    y_sub <- y_sub[keep_sub, , keep.lib.sizes = FALSE]
    y_sub <- calcNormFactors(y_sub, method = "TMM")
    cat("  Per-tissue gene set: ", nrow(y_sub), "\n")

    # Build contrasts based on which timepoints exist
    times_present <- levels(droplevels(meta_sub$time))
    cat("  Timepoints: ", paste(times_present, collapse = ", "), "\n")

    contrasts_named <- list()
    if ("11HPA" %in% times_present) {
        contrasts_named[["CT_vs_ND_11HPA"]] <- "CT_11HPA - ND_11HPA"
        contrasts_named[["CT_24HPA_vs_11HPA"]] <- "CT_24HPA - CT_11HPA"
        contrasts_named[["ND_24HPA_vs_11HPA"]] <- "ND_24HPA - ND_11HPA"
        contrasts_named[["interaction_24v11"]] <- "(CT_24HPA - CT_11HPA) - (ND_24HPA - ND_11HPA)"
    }
    contrasts_named[["CT_vs_ND_24HPA"]] <- "CT_24HPA - ND_24HPA"
    contrasts_named[["CT_vs_ND_7DPA"]]  <- "CT_7DPA - ND_7DPA"
    contrasts_named[["CT_7DPA_vs_24HPA"]] <- "CT_7DPA - CT_24HPA"
    contrasts_named[["ND_7DPA_vs_24HPA"]] <- "ND_7DPA - ND_24HPA"
    # KEY contrast (per SO-03):
    contrasts_named[["interaction_7Dv24H"]] <- "(CT_7DPA - CT_24HPA) - (ND_7DPA - ND_24HPA)"
    contrasts_str <- unlist(contrasts_named)

    # PRIMARY: voomWithQualityWeights
    cat("\n[", format(Sys.time(), "%FT%T%z"), "] Primary DE: limma-voom + qualityWeights\n")
    primary <- run_voom_de(y_sub, meta_sub, contrasts_str, paste0(tiss, "_primary"), use_qw = TRUE)
    cat("  Quality weights (range): [",
        round(min(primary$v$sample.weights), 3), ",",
        round(max(primary$v$sample.weights), 3), "]\n")
    cat("  Significance summary (q<0.05 / q<0.25):\n")
    print(primary$summary)

    saveRDS(primary, file.path(DEG_DIR, paste0("voom_", tiss, "_primary.rds")))
    write.table(primary$summary,
                file.path(DEG_DIR, paste0("voom_", tiss, "_primary_summary.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # Excel workbook with one sheet per contrast
    wb <- createWorkbook()
    for (cn in names(primary$results)) {
        addWorksheet(wb, cn)
        writeData(wb, cn, primary$results[[cn]])
    }
    addWorksheet(wb, "summary")
    writeData(wb, "summary", primary$summary)
    saveWorkbook(wb, file.path(DEG_DIR, paste0("voom_", tiss, "_primary.xlsx")), overwrite = TRUE)

    all_summaries[[paste0(tiss, "_primary")]] <- primary$summary
    all_results[[paste0(tiss, "_primary")]] <- primary$results

    # SENSITIVITY (BRAIN only): drop CT_T3_7DPA_BRAIN per SO-04
    if (tiss == "BRAIN") {
        cat("\n[", format(Sys.time(), "%FT%T%z"), "] Sensitivity DE: drop CT_T3_7DPA_BRAIN (SO-04)\n")
        keep_samples <- meta_sub$sample != "CT_T3_7DPA_BRAIN_S27"
        y_sens <- y_sub[, keep_samples, keep.lib.sizes = FALSE]
        meta_sens <- meta_sub[keep_samples, ]
        y_sens <- calcNormFactors(y_sens, method = "TMM")
        sens <- run_voom_de(y_sens, meta_sens, contrasts_str, paste0(tiss, "_sensitivity"),
                            use_qw = TRUE)
        cat("  Sensitivity significance summary:\n")
        print(sens$summary)
        saveRDS(sens, file.path(DEG_DIR, paste0("voom_", tiss, "_sensitivity.rds")))
        write.table(sens$summary,
                    file.path(DEG_DIR, paste0("voom_", tiss, "_sensitivity_summary.tsv")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        wb2 <- createWorkbook()
        for (cn in names(sens$results)) {
            addWorksheet(wb2, cn)
            writeData(wb2, cn, sens$results[[cn]])
        }
        addWorksheet(wb2, "summary")
        writeData(wb2, "summary", sens$summary)
        saveWorkbook(wb2, file.path(DEG_DIR, paste0("voom_", tiss, "_sensitivity.xlsx")),
                     overwrite = TRUE)
        all_summaries[[paste0(tiss, "_sensitivity")]] <- sens$summary
        all_results[[paste0(tiss, "_sensitivity")]] <- sens$results
    }

    # CONFIRMATORY DESeq2 LRT for INTERACTION term
    cat("\n[", format(Sys.time(), "%FT%T%z"), "] DESeq2 LRT for interaction (treatment:time)...\n")
    if ("11HPA" %in% times_present) {
        # BLAST: use 24HPA + 7DPA timepoints (since 11HPA pattern is in interaction_24v11)
        # Run two LRT tests: early (11 v 24) and late (24 v 7)
        for (subset_name in c("early_24v11", "late_7v24")) {
            if (subset_name == "early_24v11") {
                keep_t <- meta_sub$time %in% c("11HPA", "24HPA")
            } else {
                keep_t <- meta_sub$time %in% c("24HPA", "7DPA")
            }
            ms <- meta_sub[keep_t, ]
            ms$time <- droplevels(ms$time)
            y_t <- y_sub[, keep_t, keep.lib.sizes = FALSE]
            y_t <- calcNormFactors(y_t, method = "TMM")
            lrt <- run_deseq2_lrt(y_t, ms,
                                  full_formula = ~ time + treatment + treatment:time,
                                  reduced_formula = ~ time + treatment,
                                  contrast_label = paste0(tiss, "_", subset_name, "_interaction_LRT"))
            saveRDS(lrt, file.path(DEG_DIR, paste0("deseq2_lrt_", tiss, "_", subset_name, ".rds")))
            cat("  ", subset_name, " interaction LRT:\n")
            print(lrt$summary)
            all_summaries[[paste0("deseq2_", tiss, "_", subset_name)]] <- lrt$summary
        }
    } else {
        # BRAIN: just 24HPA + 7DPA
        keep_t <- meta_sub$time %in% c("24HPA", "7DPA")
        ms <- meta_sub[keep_t, ]
        ms$time <- droplevels(ms$time)
        y_t <- y_sub[, keep_t, keep.lib.sizes = FALSE]
        y_t <- calcNormFactors(y_t, method = "TMM")
        lrt <- run_deseq2_lrt(y_t, ms,
                              full_formula = ~ time + treatment + treatment:time,
                              reduced_formula = ~ time + treatment,
                              contrast_label = paste0(tiss, "_late_7v24_interaction_LRT"))
        saveRDS(lrt, file.path(DEG_DIR, paste0("deseq2_lrt_", tiss, "_late_7v24.rds")))
        cat("  late 7v24 interaction LRT:\n")
        print(lrt$summary)
        all_summaries[[paste0("deseq2_", tiss, "_late_7v24")]] <- lrt$summary
    }
}

# ============================================================
#  Combined summary
# ============================================================

combined_summary <- do.call(rbind, all_summaries)
combined_summary$run <- rownames(combined_summary)
rownames(combined_summary) <- NULL
write.table(combined_summary,
            file.path(DEG_DIR, "all_de_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n===============================\n")
cat("All DE runs combined summary:\n")
cat("===============================\n")
print(combined_summary)

cat("\n[", format(Sys.time(), "%FT%T%z"), "] Stage D done.\n")
