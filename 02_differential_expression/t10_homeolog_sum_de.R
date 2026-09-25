#!/usr/bin/env Rscript
# Homeolog-summed sensitivity DE.
# Per Stewart's 2026-05-01 spec: SECONDARY analysis, do NOT overwrite primary.
# Map .L/.S (and lowercase .l/.s) variants to a parent symbol; sum raw
# featureCounts BEFORE limma-voom (voom assumes integer counts).
# For X.<num>.<L|S> numbered paralogs, strip only the .<L|S> suffix —
# numbered variants are summed separately under X.<num>, not X.
# LOC* and tRNA/MGC*/symbol-only genes are kept as-is.
#
# Output:
#   deg/voom_BRAIN_homeolog_sum.{xlsx,rds}, voom_BRAIN_homeolog_sum_summary.tsv
#   deg/voom_BLAST_homeolog_sum.{xlsx,rds}, voom_BLAST_homeolog_sum_summary.tsv
#   deg/all_de_summary_homeolog_sum.tsv  — combined summary
#   deg/v1_manuscript_concordance_homeolog_sum.tsv — concordance vs Stewart's named-gene list

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(openxlsx)
    library(data.table)
})

PROJ <- "/project/def-srussell/nirosha_xeno"
DEG_DIR <- file.path(PROJ, "deg")
COUNTS_PATH <- file.path(PROJ, "counts", "star_gene_counts.tsv")
META_PATH <- file.path(DEG_DIR, "sample_metadata.rds")

cat("[", format(Sys.time(), "%FT%T%z"), "] Loading counts and metadata...\n", sep = "")
ct <- fread(COUNTS_PATH, header = TRUE)
meta <- readRDS(META_PATH)

stopifnot("gene_id" %in% colnames(ct))
gene_id <- ct$gene_id
mat_full <- as.matrix(ct[, -"gene_id"])
storage.mode(mat_full) <- "integer"
rownames(mat_full) <- gene_id
cat("Raw matrix:", nrow(mat_full), "genes ×", ncol(mat_full), "samples\n")

# ------------------------------------------------------------
# 1. Derive parent symbol
# ------------------------------------------------------------
# Rule: strip trailing .L, .S, .l, .s (one suffix only). Numbered paralogs
# X.<n>.L become X.<n>; LOC* untouched; tRNA/MGC* untouched.
parent_symbol <- gsub("\\.[LSls]$", "", gene_id)
n_distinct_parents <- length(unique(parent_symbol))
n_collapsed <- nrow(mat_full) - n_distinct_parents
cat("Parent symbols:", n_distinct_parents, " (collapsed", n_collapsed, "rows)\n")

# Verify: SYNPO2L, IL1B, COL17A1 should now exist as parent symbols if the L/S pair was present
key_targets <- c("synpo2l", "il1b", "col17a1", "col1a1", "col1a2", "col3a1", "lamb3",
                 "mmp8", "ptgs2", "tlr5", "ifitm1", "irf1", "irf8", "ccl5", "itga2b",
                 "itgb3", "thbs1", "slc25a37", "cpt1b", "cox7a2l", "synm", "syne3",
                 "coq10b", "cskmt", "gdf15", "synpo2")
for (k in key_targets) {
    hits <- sort(unique(parent_symbol[tolower(parent_symbol) == k]))
    if (length(hits) > 0) {
        n_homeologs <- sum(tolower(gene_id) %in% paste0(k, c(".l", ".s")))
        cat("  ", k, " parent present:", paste(hits, collapse = "/"),
            " (homeologs collapsed:", n_homeologs, ")\n", sep = "")
    } else {
        cat("  ", k, ": ABSENT in XENLA_10.1\n", sep = "")
    }
}

# ------------------------------------------------------------
# 2. Aggregate counts by parent symbol (sum)
# ------------------------------------------------------------
# Use rowsum for fast aggregation. Counts must remain integer for voom.
mat_homeolog <- rowsum(mat_full, group = parent_symbol, reorder = FALSE)
storage.mode(mat_homeolog) <- "integer"
cat("Homeolog-summed matrix:", nrow(mat_homeolog), "parent symbols ×", ncol(mat_homeolog), "samples\n")

# ------------------------------------------------------------
# 3. Build per-tissue DGELists + filter
# ------------------------------------------------------------
# Per t05_stage_D_de.R: filter cpm > 1 in >= 3 samples within tissue, TMM normalize.
run_voom_de_tissue <- function(tiss) {
    sub_idx <- meta$tissue == tiss
    cat("\n=== Tissue:", tiss, "(n=", sum(sub_idx), ") ===\n")
    counts_sub <- mat_homeolog[, meta$sample[sub_idx], drop = FALSE]
    meta_sub <- meta[sub_idx, ]

    y <- DGEList(counts = counts_sub)
    cpm_sub <- cpm(y)
    keep <- rowSums(cpm_sub > 1) >= 3
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y, method = "TMM")
    cat("Per-tissue gene set:", nrow(y), "\n")

    # Build group factor and design (no intercept) — same as primary
    grp <- factor(paste(meta_sub$treatment, meta_sub$time, sep = "_"))
    design <- model.matrix(~ 0 + grp)
    colnames(design) <- levels(grp)

    v <- voomWithQualityWeights(y, design = design, plot = FALSE)
    fit <- lmFit(v, design)

    # Contrasts — same 5 per tissue (BLAST adds early 11→24 set)
    times_present <- levels(droplevels(meta_sub$time))
    contrasts_named <- list()
    if ("11HPA" %in% times_present) {
        contrasts_named[["CT_vs_ND_11HPA"]]    <- "CT_11HPA - ND_11HPA"
        contrasts_named[["CT_24HPA_vs_11HPA"]] <- "CT_24HPA - CT_11HPA"
        contrasts_named[["ND_24HPA_vs_11HPA"]] <- "ND_24HPA - ND_11HPA"
        contrasts_named[["interaction_24v11"]] <- "(CT_24HPA - CT_11HPA) - (ND_24HPA - ND_11HPA)"
    }
    contrasts_named[["CT_vs_ND_24HPA"]]    <- "CT_24HPA - ND_24HPA"
    contrasts_named[["CT_vs_ND_7DPA"]]     <- "CT_7DPA - ND_7DPA"
    contrasts_named[["CT_7DPA_vs_24HPA"]]  <- "CT_7DPA - CT_24HPA"
    contrasts_named[["ND_7DPA_vs_24HPA"]]  <- "ND_7DPA - ND_24HPA"
    contrasts_named[["interaction_7Dv24H"]]<- "(CT_7DPA - CT_24HPA) - (ND_7DPA - ND_24HPA)"

    contrasts_str <- unlist(contrasts_named)
    contr_obj <- makeContrasts(contrasts = contrasts_str, levels = design)
    colnames(contr_obj) <- names(contrasts_named)
    fit2 <- contrasts.fit(fit, contr_obj)
    fit2 <- eBayes(fit2)

    res_list <- lapply(colnames(contr_obj), function(cn) {
        tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "P")
        tt$gene_id <- rownames(tt)
        tt$contrast <- cn
        tt[, c("gene_id", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "contrast")]
    })
    names(res_list) <- colnames(contr_obj)

    sig_summary <- data.frame(
        contrast = names(res_list),
        n_total  = sapply(res_list, nrow),
        sig_q05  = sapply(res_list, function(d) sum(d$adj.P.Val < 0.05)),
        sig_q25  = sapply(res_list, function(d) sum(d$adj.P.Val < 0.25)),
        up_q05   = sapply(res_list, function(d) sum(d$adj.P.Val < 0.05 & d$logFC > 0)),
        down_q05 = sapply(res_list, function(d) sum(d$adj.P.Val < 0.05 & d$logFC < 0)),
        stringsAsFactors = FALSE
    )
    cat("Significance summary:\n"); print(sig_summary)

    list(v = v, fit = fit2, results = res_list, summary = sig_summary,
         contrasts = contr_obj, design = design)
}

results_brain <- run_voom_de_tissue("BRAIN")
results_blast <- run_voom_de_tissue("BLAST")

# ------------------------------------------------------------
# 4. Save outputs
# ------------------------------------------------------------
save_tissue <- function(res, tiss) {
    saveRDS(res, file.path(DEG_DIR, paste0("voom_", tiss, "_homeolog_sum.rds")))
    write.table(res$summary,
                file.path(DEG_DIR, paste0("voom_", tiss, "_homeolog_sum_summary.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    wb <- createWorkbook()
    for (cn in names(res$results)) {
        addWorksheet(wb, cn)
        writeData(wb, cn, res$results[[cn]])
    }
    addWorksheet(wb, "summary")
    writeData(wb, "summary", res$summary)
    saveWorkbook(wb, file.path(DEG_DIR, paste0("voom_", tiss, "_homeolog_sum.xlsx")),
                 overwrite = TRUE)
}
save_tissue(results_brain, "BRAIN")
save_tissue(results_blast, "BLAST")

# Combined summary
combined <- rbind(
    cbind(tissue = "BRAIN", results_brain$summary),
    cbind(tissue = "BLAST", results_blast$summary)
)
write.table(combined, file.path(DEG_DIR, "all_de_summary_homeolog_sum.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ------------------------------------------------------------
# 5. v1 manuscript concordance
# ------------------------------------------------------------
# Stewart's named-gene list (from /loop prompt; presumed BLAST since it's the
# fibrotic ECM/immune signature). 24h and 7d entries.
v1 <- data.table(
    gene  = tolower(c("COL17A1","SYNPO2L","COQ10B","CSKMT",
                       "COL1A1","COL1A2","COL3A1","LAMB3","MMP8","IL1B","PTGS2","TLR5",
                       "IFITM1","IRF1","IRF8","CCL5","ITGA2B","ITGB3","THBS1","SLC25A37",
                       "CPT1B","COX7A2L","SYNM","SYNE3")),
    timepoint = c(rep("24h", 4), rep("7d", 20)),
    v1_logFC  = c(-1.90, +3.65, +0.79, +0.31,
                   -1.83, -2.10, -1.15, -4.38, -4.20, -3.02, -1.16, -3.42,
                   +4.12, +1.76, +0.75, +0.79, +1.54, +1.65, +0.97, +3.29,
                   +1.38, +0.32, +2.96, +0.30)
)

# Look up these genes in BLAST results (homeolog-summed parent symbols are lowercase)
get_logfc <- function(res, gene, contrast_name) {
    tt <- res$results[[contrast_name]]
    # Normalise: lowercase + strip hyphens (Xenbase symbol convention uses hyphens for some
    # canonical names, e.g. "timp-1.L" while v1 manuscript uses "TIMP1"; same gene).
    norm <- function(x) gsub("-", "", tolower(x))
    hit <- tt[norm(tt$gene_id) == norm(gene), ]
    if (nrow(hit) == 0) return(list(logFC = NA_real_, q = NA_real_, present = FALSE))
    list(logFC = as.numeric(hit$logFC[1]), q = as.numeric(hit$adj.P.Val[1]), present = TRUE)
}

conc <- v1[, {
    contrast_name <- if (timepoint == "24h") "CT_vs_ND_24HPA" else "CT_vs_ND_7DPA"
    blast_h <- get_logfc(results_blast, gene, contrast_name)
    brain_h <- get_logfc(results_brain, gene, contrast_name)
    list(
        contrast_name = contrast_name,
        v1_logFC      = as.numeric(v1_logFC),
        BLAST_homeolog_sum_logFC = blast_h$logFC,
        BLAST_homeolog_sum_q     = blast_h$q,
        BRAIN_homeolog_sum_logFC = brain_h$logFC,
        BRAIN_homeolog_sum_q     = brain_h$q,
        BLAST_present = blast_h$present,
        BRAIN_present = brain_h$present,
        BLAST_dir_match = if (blast_h$present) sign(blast_h$logFC) == sign(v1_logFC) else NA,
        BRAIN_dir_match = if (brain_h$present) sign(brain_h$logFC) == sign(v1_logFC) else NA
    )
}, by = .(gene, timepoint)]

write.table(conc, file.path(DEG_DIR, "v1_manuscript_concordance_homeolog_sum.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== v1 manuscript concordance (homeolog-summed) ===\n")
print(conc)
cat("\nBLAST present:", sum(conc$BLAST_present), "/", nrow(conc), "\n")
cat("BLAST direction match:", sum(conc$BLAST_dir_match, na.rm = TRUE), "/", sum(conc$BLAST_present), "\n")
cat("BRAIN present:", sum(conc$BRAIN_present), "/", nrow(conc), "\n")

cat("\n[", format(Sys.time(), "%FT%T%z"), "] Homeolog-sum DE complete.\n", sep = "")
