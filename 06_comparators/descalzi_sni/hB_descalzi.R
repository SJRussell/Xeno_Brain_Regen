#!/usr/bin/env Rscript
# H-B comparator: Descalzi 2017 mouse SNI bulk RNA-seq (GSE91396)
# 36 samples × 3 regions (mPFC, NAc, PAG) × 2 conditions (sham vs SNI) × ~6 reps
# Pre-computed RPKM matrix (GSE91396_rpkm.txt.gz, 7.3 MB) — limma+trend on log2(RPKM+1).
# Cross-species concordance with Xenopus 7DPA CT signature via Ensembl Compara orthologs.
#
# Output dir: /scratch/srussell/srussell/nirosha_xeno_comparators/descalzi/
#   work/de_<region>.tsv
#   work/gsea_<region>.tsv
#   work/xenopus_concordance.tsv
#   work/figures/{volcano_<region>, hallmark_concordance}.pdf
# Final: copied to /project/.../comparators/descalzi/

suppressPackageStartupMessages({
    library(GEOquery)
    library(limma)
    library(edgeR)
    library(biomaRt)
    library(clusterProfiler)
    library(msigdbr)
    library(data.table)
    library(ggplot2)
})

PROJ <- "/project/def-srussell/nirosha_xeno"
DEG_DIR <- file.path(PROJ, "deg")
PATHWAY_DIR <- file.path(PROJ, "pathway")
WORK <- "/scratch/srussell/srussell/nirosha_xeno_comparators/descalzi/work"
RAW  <- "/scratch/srussell/srussell/nirosha_xeno_comparators/descalzi/raw"
FIG  <- file.path(WORK, "figures")
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
dir.create(RAW,  recursive = TRUE, showWarnings = FALSE)
dir.create(FIG,  recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Download GSE91396 expression + metadata
# ------------------------------------------------------------
GEO_DIR <- file.path(RAW, "GSE91396")
dir.create(GEO_DIR, showWarnings = FALSE)
options(timeout = 600)
cat("[", format(Sys.time()), "] Fetching GSE91396 from GEO...\n", sep = "")
gse <- getGEO("GSE91396", destdir = GEO_DIR, GSEMatrix = TRUE, getGPL = FALSE)
eset <- gse[[1]]
pdata <- pData(eset)
cat("Sample table:\n")
print(table(pdata$source_name_ch1))
cat("Sample title head:\n")
print(head(pdata$title))

# Download supplementary RPKM matrix (preferred — eset's matrix may be empty for newer GEO entries)
suppl_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE91396&format=file&file=GSE91396_rpkm.txt.gz"
suppl_path <- file.path(RAW, "GSE91396_rpkm.txt.gz")
if (!file.exists(suppl_path)) {
    download.file(suppl_url, suppl_path, mode = "wb")
}

# Read RPKM matrix
rpkm <- fread(suppl_path, header = TRUE)
cat("RPKM matrix:", dim(rpkm)[1], "genes ×", dim(rpkm)[2], "cols\n")
cat("Cols head:", paste(head(colnames(rpkm), 6), collapse = ", "), "\n")

# GSE91396_rpkm.txt.gz layout (verified):
#   col 1     : gene_id (ENSEMBL ID, e.g. ENSMUSG00000000001)
#   cols 2-37 : 36 sample RPKM columns named "SampleN_<COND>_<REGION>"
#   col 38    : gene_name (e.g. "Gnai3" — what we want for cross-species matching)
#   col 39    : biotype
#   col 40    : length
sample_cols <- grep("^Sample\\d+_", colnames(rpkm), value = TRUE)
gene_id_col <- "gene_id"
gene_name_col <- "gene_name"
stopifnot(length(sample_cols) > 0, gene_id_col %in% colnames(rpkm))
cat("Detected", length(sample_cols), "sample columns\n")

mat <- as.matrix(rpkm[, ..sample_cols])
storage.mode(mat) <- "numeric"  # belt-and-suspenders against any stray non-numeric

# Use gene_name (mouse symbol) as primary identifier for cross-species comparison.
# Fall back to gene_id if gene_name is missing for a row.
gene_label <- rpkm[[gene_name_col]]
gene_label[is.na(gene_label) | gene_label == ""] <- rpkm[[gene_id_col]][is.na(gene_label) | gene_label == ""]
rownames(mat) <- gene_label
mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
mat <- mat[rownames(mat) != "" & !is.na(rownames(mat)), , drop = FALSE]
expr_cols <- colnames(mat)
cat("RPKM matrix post-dedup:", dim(mat)[1], "genes ×", dim(mat)[2], "samples\n")

# ------------------------------------------------------------
# 2. Build sample metadata: region × condition
# ------------------------------------------------------------
# Match pdata rows to RPKM columns. Common GEO pattern: title contains region+condition+rep.
# Extract from pdata if columns align; otherwise parse title strings.
build_meta <- function(pdata, expr_cols) {
    # Try to match expr_cols to pdata rownames or pdata$title
    if (all(expr_cols %in% rownames(pdata))) {
        out <- pdata[expr_cols, , drop = FALSE]
    } else {
        # Try matching by title fragments
        idx <- sapply(expr_cols, function(c) {
            hits <- which(grepl(c, pdata$title, fixed = TRUE) |
                          grepl(c, pdata$description, fixed = TRUE) |
                          grepl(c, pdata$geo_accession, fixed = TRUE))
            if (length(hits) == 1) hits else NA_integer_
        })
        if (any(is.na(idx))) {
            cat("WARN: could not match all expr_cols to pdata. Falling back to ordinal mapping.\n")
            idx <- seq_along(expr_cols)
        }
        out <- pdata[idx, , drop = FALSE]
    }
    out
}
meta <- build_meta(pdata, expr_cols)

# Parse region (mPFC/NAc/PAG) and condition (sham/SNI) from title or characteristics
parse_region <- function(s) {
    s <- tolower(s)
    if (grepl("pfc", s) || grepl("prefrontal", s)) return("mPFC")
    if (grepl("nac", s) || grepl("accumbens", s))  return("NAc")
    if (grepl("pag", s) || grepl("periaqueductal", s)) return("PAG")
    NA_character_
}
parse_condition <- function(s) {
    s <- tolower(s)
    if (grepl("sni", s) || grepl("nerve injury", s) || grepl("spared", s)) return("SNI")
    if (grepl("sham", s)) return("sham")
    NA_character_
}

region_src <- if (!is.null(meta$source_name_ch1)) meta$source_name_ch1 else meta$title
cond_src   <- meta$title
if (!is.null(meta$characteristics_ch1)) cond_src <- paste(cond_src, meta$characteristics_ch1)

meta$region    <- vapply(region_src, parse_region,    character(1))
meta$condition <- vapply(cond_src,   parse_condition, character(1))
cat("\nParsed sample table:\n")
print(table(meta$region, meta$condition, useNA = "ifany"))

# Drop unparseable
keep <- !is.na(meta$region) & !is.na(meta$condition)
meta <- meta[keep, , drop = FALSE]
mat  <- mat[, expr_cols[keep], drop = FALSE]
cat("After parse-filter:", ncol(mat), "samples retained\n")

# ------------------------------------------------------------
# 3. Per-region limma+trend DE: SNI vs sham
# ------------------------------------------------------------
log_rpkm <- log2(mat + 1)
de_results <- list()
gsea_results <- list()
H_msig <- msigdbr(species = "Mus musculus", category = "H")
H_term2gene <- as.data.frame(H_msig[, c("gs_name", "gene_symbol")])

for (region in unique(meta$region)) {
    cat("\n=== Region:", region, "===\n")
    sel <- meta$region == region
    sub_mat <- log_rpkm[, sel, drop = FALSE]
    sub_meta <- meta[sel, , drop = FALSE]
    cond <- factor(sub_meta$condition, levels = c("sham", "SNI"))
    design <- model.matrix(~ cond)
    fit <- lmFit(sub_mat, design)
    fit <- eBayes(fit, trend = TRUE)
    tt  <- topTable(fit, coef = "condSNI", number = Inf, sort.by = "none")
    tt$gene <- rownames(tt)
    de_dt   <- as.data.table(tt)
    setorder(de_dt, P.Value)
    de_dt[, region := region]
    fwrite(de_dt, file.path(WORK, paste0("de_", region, ".tsv")), sep = "\t")
    de_results[[region]] <- de_dt
    cat("  Top 5 by P:\n"); print(head(de_dt, 5))
    cat("  q<0.05:", sum(de_dt$adj.P.Val < 0.05), "\n")

    # GSEA on Hallmark, signed -log10(P) ranked metric
    rank_metric <- with(de_dt, sign(logFC) * -log10(P.Value))
    names(rank_metric) <- de_dt$gene
    rank_metric <- sort(rank_metric, decreasing = TRUE)
    gsea_res <- tryCatch({
        GSEA(rank_metric, TERM2GENE = H_term2gene,
             pvalueCutoff = 0.25, eps = 0, scoreType = "std")
    }, error = function(e) {
        cat("  GSEA error:", conditionMessage(e), "\n")
        NULL
    })
    if (!is.null(gsea_res)) {
        gd <- as.data.table(gsea_res@result)
        gd[, region := region]
        fwrite(gd, file.path(WORK, paste0("gsea_", region, ".tsv")), sep = "\t")
        gsea_results[[region]] <- gd
        cat("  Hallmark sig (q<0.25):", nrow(gd[qvalue < 0.25]), "\n")
    }

    # Volcano
    pdf(file.path(FIG, paste0("volcano_", region, ".pdf")), width = 6, height = 5)
    de_dt[, sig := adj.P.Val < 0.05]
    print(ggplot(de_dt, aes(logFC, -log10(P.Value), color = sig)) +
              geom_point(alpha = 0.4, size = 0.6) +
              scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "red")) +
              theme_minimal() +
              labs(title = paste0("Descalzi 2017 ", region, ": SNI vs sham (log2 RPKM, limma+trend)")))
    dev.off()
}

# ------------------------------------------------------------
# 4. Cross-species concordance: Xenopus 7DPA CT brain signature vs mouse SNI brain
# ------------------------------------------------------------
# Load Xenopus brain DE (CT vs ND at 7DPA from Stage D)
xen_de <- tryCatch(readRDS(file.path(DEG_DIR, "voom_BRAIN_primary.rds")),
                   error = function(e) NULL)
xen_table <- NULL
if (!is.null(xen_de) && is.list(xen_de) && "results" %in% names(xen_de)) {
    nm <- names(xen_de$results)
    target <- grep("CT_vs_ND.*7DPA|7DPA.*CT_vs_ND", nm, value = TRUE, ignore.case = TRUE)
    if (length(target) >= 1) {
        xt <- xen_de$results[[target[1]]]
        xen_table <- data.table(gene = rownames(xt), logFC = xt$logFC,
                                P = xt$P.Value, adj.P = xt$adj.P.Val)
        # Strip .L/.S subgenome suffix to compare with mouse symbols
        xen_table[, gene_clean := gsub("\\.[LS]$", "", gene)]
        xen_table[, gene_clean := toupper(gene_clean)]
    } else {
        cat("WARN: could not find CT_vs_ND_7DPA in voom_BRAIN_primary$results (names: ",
            paste(nm, collapse = ", "), ")\n")
    }
}

if (!is.null(xen_table)) {
    cat("\n=== Cross-species concordance ===\n")
    concordance <- list()
    for (region in names(de_results)) {
        mouse_de <- de_results[[region]]
        mouse_de[, gene_upper := toupper(gene)]
        # Direct symbol overlap (axolotl + Xenopus + mouse share many gene symbols verbatim)
        m <- merge(xen_table[, .(gene_clean, xen_logFC = logFC, xen_p = P)],
                   mouse_de[, .(gene_upper, mouse_logFC = logFC, mouse_p = P.Value)],
                   by.x = "gene_clean", by.y = "gene_upper")
        if (nrow(m) > 100) {
            cat(region, ": overlap n =", nrow(m), "\n")
            ct <- cor.test(m$xen_logFC, m$mouse_logFC, method = "spearman")
            cat("  Spearman ρ =", round(ct$estimate, 3), "p =", format(ct$p.value, digits = 3), "\n")
            # Gene-level concordance: direction agreement among genes with sig in either species
            sig_xen <- m[xen_p < 0.05]
            agree <- sum(sign(sig_xen$xen_logFC) == sign(sig_xen$mouse_logFC))
            cat("  Direction agreement among Xen-sig genes:", agree, "/", nrow(sig_xen), "\n")
            concordance[[region]] <- data.table(
                region = region, n_overlap = nrow(m),
                spearman_rho = ct$estimate, spearman_p = ct$p.value,
                xen_sig_n = nrow(sig_xen), agree = agree)
        }
    }
    if (length(concordance) > 0) {
        conc_dt <- rbindlist(concordance)
        fwrite(conc_dt, file.path(WORK, "xenopus_concordance.tsv"), sep = "\t")
        cat("Concordance summary written.\n")
        print(conc_dt)
    }
}

# ------------------------------------------------------------
# 5. Hallmark NES correlation: mouse Descalzi vs Xenopus brain
# ------------------------------------------------------------
xen_gsea_path <- file.path(PATHWAY_DIR, "all_enrichment_summary.tsv")
if (file.exists(xen_gsea_path) && length(gsea_results) > 0) {
    xen_gsea <- fread(xen_gsea_path)
    # Filter to Xenopus BRAIN, CT_vs_ND 7DPA contrast, Hallmark
    xen_h <- xen_gsea[grepl("BRAIN", contrast, ignore.case = TRUE) &
                      grepl("CT_vs_ND.*7DPA|7DPA.*CT_vs_ND", contrast, ignore.case = TRUE) &
                      grepl("hallmark", db, ignore.case = TRUE)]
    if (nrow(xen_h) > 0) {
        for (region in names(gsea_results)) {
            mouse_h <- gsea_results[[region]]
            joined <- merge(xen_h[, .(ID, xen_NES = NES)],
                            mouse_h[, .(ID, mouse_NES = NES)], by = "ID")
            if (nrow(joined) >= 5) {
                ct <- cor.test(joined$xen_NES, joined$mouse_NES, method = "spearman")
                cat("\nHallmark NES correlation: Xen-BRAIN-7DPA-CT vs Mouse-", region,
                    ": ρ=", round(ct$estimate, 3), "p=", format(ct$p.value, digits = 3),
                    " (n=", nrow(joined), ")\n", sep = "")
                fwrite(joined, file.path(WORK, paste0("hallmark_nes_xen_vs_", region, ".tsv")), sep = "\t")
            }
        }
    } else {
        cat("WARN: no Xenopus BRAIN CT_vs_ND_7DPA Hallmark entries found in pathway/all_enrichment_summary.tsv. Column names: ",
            paste(colnames(xen_gsea), collapse = ", "), "\n")
    }
}

# ------------------------------------------------------------
# 6. Copy small final outputs to /project mirror
# ------------------------------------------------------------
PROJ_OUT <- file.path(PROJ, "comparators", "descalzi")
dir.create(PROJ_OUT, recursive = TRUE, showWarnings = FALSE)
for (f in list.files(WORK, pattern = "\\.tsv$", full.names = TRUE)) {
    file.copy(f, PROJ_OUT, overwrite = TRUE)
}
PROJ_FIG <- file.path(PROJ, "figures", "stage_H", "descalzi")
dir.create(PROJ_FIG, recursive = TRUE, showWarnings = FALSE)
for (f in list.files(FIG, pattern = "\\.pdf$", full.names = TRUE)) {
    file.copy(f, PROJ_FIG, overwrite = TRUE)
}

cat("[", format(Sys.time()), "] hB_descalzi.R complete.\n", sep = "")
