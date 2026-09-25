#!/usr/bin/env Rscript
# GSEA on homeolog-summed primary DE results.
# Mirrors t06_stage_E_enrichment.R but pointed at voom_*_homeolog_sum.rds.
# Strand-correct count matrix, parent-symbol space.
#
# Output:
#   pathway/homeolog_sum_strand_corrected_enrichment.tsv  (combined long table)
#   pathway/homeolog_sum_strand_corrected_results.rds     (per-contrast RDS)
#   pathway/voom_BRAIN_homeolog_sum_*_enrichment.xlsx
#   pathway/voom_BLAST_homeolog_sum_*_enrichment.xlsx

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Xl.eg.db)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(msigdbr)
    library(openxlsx)
})

PROJ <- "/project/def-srussell/nirosha_xeno"
DEG_DIR  <- file.path(PROJ, "deg")
PATH_DIR <- file.path(PROJ, "pathway")
dir.create(PATH_DIR, recursive = TRUE, showWarnings = FALSE)

cat("[", format(Sys.time(), "%FT%T%z"), "] Loading homeolog-summed DE results...\n", sep = "")

# ----- MSigDB Hallmark -----
hallmark_human <- msigdbr(species = "Homo sapiens", collection = "H")
entrez_col <- intersect(c("ncbi_gene", "human_entrez_gene", "entrez_gene"),
                        colnames(hallmark_human))[1]
hallmark_h <- as.data.frame(hallmark_human[, c("gs_name", entrez_col)])
colnames(hallmark_h) <- c("gs_name", "entrez")
hallmark_h <- unique(hallmark_h[!is.na(hallmark_h$entrez), ])
cat("Hallmark gene sets:", length(unique(hallmark_h$gs_name)), "\n")

# ----- Symbol mappings -----
# In homeolog-summed space, gene_ids are PARENT symbols (no .L/.S suffix), so
# the strip step is a no-op for already-summed entries but still safe.
strip_xla_suffix <- function(s) sub("\\.[LSls]$", "", s)
sym_hs <- keys(org.Hs.eg.db, keytype = "SYMBOL")
sym_to_hs_entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                          keys = sym_hs, keytype = "SYMBOL",
                                          columns = c("SYMBOL", "ENTREZID"))
hs_lookup <- setNames(sym_to_hs_entrez$ENTREZID,
                      toupper(sym_to_hs_entrez$SYMBOL))
xla_to_hs_entrez <- function(syms) hs_lookup[toupper(strip_xla_suffix(syms))]

# ----- Helpers -----
build_rank <- function(res_df) {
    gid <- if (!is.null(res_df$gene_id)) res_df$gene_id else rownames(res_df)
    rk <- sign(res_df$logFC) * (-log10(res_df$P.Value))
    names(rk) <- gid
    rk <- rk[!is.na(rk) & is.finite(rk)]
    sort(rk, decreasing = TRUE)
}

run_gsea_hallmark <- function(rank_vec, contrast_name) {
    hs_entrez <- xla_to_hs_entrez(names(rank_vec))
    keep <- !is.na(hs_entrez) & !duplicated(hs_entrez)
    rk_h <- rank_vec[keep]
    names(rk_h) <- hs_entrez[keep]
    rk_h <- sort(rk_h, decreasing = TRUE)
    if (length(rk_h) < 100) {
        cat("  Skip Hallmark (", contrast_name, "): n=", length(rk_h), "\n", sep = ""); return(NULL)
    }
    set.seed(42)
    res <- tryCatch({
        GSEA(rk_h, TERM2GENE = hallmark_h[, c("gs_name", "entrez")],
             pvalueCutoff = 0.25, minGSSize = 15, maxGSSize = 500, eps = 0)
    }, error = function(e) { cat("    Hallmark err:", conditionMessage(e), "\n"); NULL })
    if (is.null(res) || nrow(res@result) == 0) {
        cat("  Hallmark: 0 enriched (", contrast_name, ")\n", sep = ""); return(NULL)
    }
    df <- res@result[, c("ID", "NES", "pvalue", "p.adjust", "qvalue", "core_enrichment")]
    df$contrast <- contrast_name; df$db <- "MSigDB_Hallmark"
    cat("  Hallmark: ", nrow(df), " (", contrast_name, ")\n", sep = "")
    df
}

run_gsea_GO <- function(rank_vec, contrast_name) {
    xla_entrez <- AnnotationDbi::select(org.Xl.eg.db,
                                        keys = names(rank_vec),
                                        keytype = "SYMBOL",
                                        columns = c("SYMBOL", "ENTREZID"))
    xla_entrez <- xla_entrez[!is.na(xla_entrez$ENTREZID) & !duplicated(xla_entrez$SYMBOL), ]
    rk_e <- rank_vec[match(xla_entrez$SYMBOL, names(rank_vec))]
    names(rk_e) <- xla_entrez$ENTREZID
    rk_e <- rk_e[!is.na(rk_e)]
    rk_e <- sort(rk_e, decreasing = TRUE)
    if (length(rk_e) < 100) {
        cat("  Skip GO BP (", contrast_name, "): n=", length(rk_e), "\n", sep = ""); return(NULL)
    }
    set.seed(42)
    res <- tryCatch({
        gseGO(rk_e, OrgDb = org.Xl.eg.db, keyType = "ENTREZID",
              ont = "BP", pvalueCutoff = 0.25, minGSSize = 15, maxGSSize = 500,
              eps = 0, verbose = FALSE)
    }, error = function(e) NULL)
    if (is.null(res) || nrow(res@result) == 0) {
        cat("  GO BP: 0 enriched (", contrast_name, ")\n", sep = ""); return(NULL)
    }
    df <- res@result[, c("ID", "Description", "NES", "pvalue", "p.adjust", "qvalue", "core_enrichment")]
    df$contrast <- contrast_name; df$db <- "GO_BP"
    cat("  GO BP: ", nrow(df), " (", contrast_name, ")\n", sep = "")
    df
}

# ----- Driver -----
all_enrich <- list()
for (run in c("deseq2_BRAIN_primary_homeolog_sum",
              "deseq2_BLAST_primary_homeolog_sum",
              "voom_BRAIN_homeolog_sum",
              "voom_BLAST_homeolog_sum")) {
    cat("\n========== ", run, " ==========\n", sep = "")
    obj <- readRDS(file.path(DEG_DIR, paste0(run, ".rds")))
    for (cn in names(obj$results)) {
        rank_vec <- build_rank(obj$results[[cn]])
        cat("\n--- ", run, "::", cn, "  (n=", length(rank_vec), " ranked)\n", sep = "")
        h <- run_gsea_hallmark(rank_vec, paste0(run, "::", cn))
        if (!is.null(h)) all_enrich[[paste0(run, "_", cn, "_hallmark")]] <- h
        g <- run_gsea_GO(rank_vec, paste0(run, "::", cn))
        if (!is.null(g)) all_enrich[[paste0(run, "_", cn, "_goBP")]] <- g
    }
}

# ----- Combine -----
all_enrich <- lapply(all_enrich, function(d) {
    if (!"Description" %in% names(d)) d$Description <- d$ID; d
})
combined <- do.call(rbind, lapply(all_enrich, function(d) {
    d[, c("contrast", "db", "ID", "Description", "NES", "pvalue", "p.adjust", "qvalue")]
}))
combined$count_provenance <- "STAR ReadsPerGene col4 (-s 2 reverse-stranded)"
combined$strand_audit_2026_05_05 <- "VERIFIED_REVERSE_STRANDED"
combined$homeolog_layer <- "homeolog_summed"

write.table(combined,
            file.path(PATH_DIR, "homeolog_sum_strand_corrected_enrichment.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nWrote ", PATH_DIR, "/homeolog_sum_strand_corrected_enrichment.tsv (",
    nrow(combined), " rows)\n", sep = "")

saveRDS(all_enrich, file.path(PATH_DIR, "homeolog_sum_strand_corrected_results.rds"))

# Per-run xlsx
xl_groups <- split(seq_along(all_enrich), sub("(_hallmark|_goBP)$", "", names(all_enrich)))
for (g in names(xl_groups)) {
    wb <- createWorkbook()
    for (i in xl_groups[[g]]) {
        nm <- names(all_enrich)[i]
        sheet <- substr(sub(paste0("^", g, "_"), "", nm), 1, 31)
        addWorksheet(wb, sheet)
        writeData(wb, sheet, all_enrich[[i]])
    }
    saveWorkbook(wb, file.path(PATH_DIR, paste0(g, "_strand_corrected_enrichment.xlsx")),
                 overwrite = TRUE)
}

cat("\n[", format(Sys.time(), "%FT%T%z"), "] Stage E homeolog-sum done.\n")
