#!/usr/bin/env Rscript
# Stage E (fgsea-based fallback): per-contrast GSEA without clusterProfiler.
# Uses fgsea + GO.db + org.Hs.eg.db. Maps Xenopus symbols → human via
# case-folded symbol match (works for ~80% of orthologs; documented as a
# methodological caveat).
#
# Gene-set sources:
#   1. GO BP via human GO annotation (GO BP terms → human Entrez via org.Hs.eg.db)
#      Mapped through Xla → human symbol case-fold.
#   2. MSigDB Hallmark v2024.1.Hs.entrez.gmt (download to /home/srussell)
#      Auto-downloaded if missing.

suppressPackageStartupMessages({
    library(fgsea)
    library(GO.db)
    library(GSEABase)
    library(org.Hs.eg.db)
    library(AnnotationDbi)
    library(openxlsx)
})

PROJ <- "/project/def-srussell/nirosha_xeno"
DEG_DIR <- file.path(PROJ, "deg")
PATH_DIR <- file.path(PROJ, "pathway")
GMT_CACHE <- file.path(PROJ, "ref", "msigdb")
dir.create(PATH_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GMT_CACHE, recursive = TRUE, showWarnings = FALSE)

cat("[", format(Sys.time(), "%FT%T%z"), "] Stage E — fgsea fallback\n")

# ----- 1. Build GO BP gene sets keyed on human symbols -----
cat("Building GO BP gene sets from org.Hs.eg.db...\n")
all_hs_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
sym_to_go <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys = all_hs_symbols, keytype = "SYMBOL",
                                    columns = c("SYMBOL", "GO", "ONTOLOGY"))
sym_to_go <- sym_to_go[sym_to_go$ONTOLOGY == "BP" & !is.na(sym_to_go$GO), ]
go_to_sym <- split(toupper(sym_to_go$SYMBOL), sym_to_go$GO)
go_to_sym <- lapply(go_to_sym, unique)
sizes <- sapply(go_to_sym, length)
go_to_sym <- go_to_sym[sizes >= 15 & sizes <= 500]
cat("  GO BP gene sets (size 15-500):", length(go_to_sym), "\n")

# Term names from GO.db
go_terms <- AnnotationDbi::select(GO.db, keys = names(go_to_sym), keytype = "GOID",
                                   columns = c("GOID", "TERM"))
go_label <- setNames(go_terms$TERM, go_terms$GOID)

# ----- 2. MSigDB Hallmark gene sets -----
hallmark_gmt <- file.path(GMT_CACHE, "h.all.v2024.1.Hs.symbols.gmt")
if (!file.exists(hallmark_gmt)) {
    url <- "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt"
    cat("Downloading MSigDB Hallmark from:", url, "\n")
    options(timeout = 300)
    res <- tryCatch(download.file(url, hallmark_gmt, quiet = TRUE),
                    error = function(e) {cat("Download failed:", conditionMessage(e), "\n"); NULL})
    if (is.null(res) || !file.exists(hallmark_gmt)) {
        cat("WARNING: Could not download Hallmark; skipping.\n")
        hallmark <- NULL
    }
}
if (file.exists(hallmark_gmt)) {
    hallmark_lines <- readLines(hallmark_gmt)
    hallmark <- lapply(hallmark_lines, function(L) {
        f <- strsplit(L, "\t", fixed = TRUE)[[1]]
        toupper(f[-(1:2)])  # gene symbols, drop name + URL
    })
    names(hallmark) <- sapply(hallmark_lines, function(L)
        strsplit(L, "\t", fixed = TRUE)[[1]][1])
    cat("  Hallmark gene sets:", length(hallmark), "\n")
}

# ----- 3. Helper: build ranked vector with case-folded symbol -----
# Strip Xla suffix (.L / .S) and uppercase to match human symbol space.
strip_xla <- function(s) toupper(sub("\\.[LS]$", "", s))

build_rank <- function(res_df) {
    rk <- sign(res_df$logFC) * (-log10(res_df$P.Value))
    names(rk) <- strip_xla(res_df$gene_id)
    # Drop NAs / non-finite up front (DESeq2 sets padj+pvalue to NA for genes
    # filtered by independent filtering; those would break which.max downstream)
    keep <- !is.na(rk) & is.finite(rk) & !is.na(names(rk)) & nzchar(names(rk))
    rk <- rk[keep]
    if (length(rk) == 0) return(numeric(0))
    # Collapse duplicates (Xla L/S subgenome alleles → one human ortholog) by max abs
    rk_split <- split(rk, names(rk))
    rk_collapse <- vapply(rk_split, function(x) {
        if (length(x) == 0) return(NA_real_)
        x[which.max(abs(x))]
    }, numeric(1))
    rk_collapse <- rk_collapse[!is.na(rk_collapse) & is.finite(rk_collapse)]
    sort(rk_collapse, decreasing = TRUE)
}

# ----- 4. Run fgsea per contrast per gene set -----
all_results <- list()

run_fgsea <- function(rank_vec, pathways, db_name, contrast_label) {
    set.seed(42)
    res <- tryCatch({
        fgsea(pathways = pathways, stats = rank_vec, minSize = 15, maxSize = 500,
              eps = 0)
    }, error = function(e) {cat("    fgsea error:", conditionMessage(e), "\n"); NULL})
    if (is.null(res) || nrow(res) == 0) return(NULL)
    res <- as.data.frame(res)
    res$contrast <- contrast_label
    res$db <- db_name
    if (db_name == "GO_BP") {
        res$Description <- go_label[res$pathway]
    } else {
        res$Description <- res$pathway
    }
    res$leadingEdge <- vapply(res$leadingEdge, paste, character(1), collapse = ";")
    res <- res[order(res$padj), ]
    res
}

for (run in c("deseq2_BRAIN_primary", "deseq2_BLAST_primary",
              "voom_BRAIN_primary",   "voom_BLAST_primary")) {
    cat("\n==========================================\n")
    cat("Run:", run, "\n")
    cat("==========================================\n")
    obj <- readRDS(file.path(DEG_DIR, paste0(run, ".rds")))
    for (cn in names(obj$results)) {
        res_df <- obj$results[[cn]]
        rank_vec <- build_rank(res_df)
        cat(sprintf("\n--- %s :: %s (n=%d ranked, after Xla collapse) ---\n",
                    run, cn, length(rank_vec)))

        if (length(rank_vec) >= 100) {
            # GO BP
            r_go <- run_fgsea(rank_vec, go_to_sym, "GO_BP", paste0(run, "::", cn))
            if (!is.null(r_go)) {
                sig <- sum(r_go$padj < 0.05, na.rm = TRUE)
                cat(sprintf("    GO BP: %d sig (q<0.05) / %d total\n", sig, nrow(r_go)))
                all_results[[paste0(run, "_", cn, "_GO_BP")]] <- r_go
            }
            # Hallmark
            if (!is.null(hallmark)) {
                r_h <- run_fgsea(rank_vec, hallmark, "MSigDB_Hallmark", paste0(run, "::", cn))
                if (!is.null(r_h)) {
                    sig <- sum(r_h$padj < 0.05, na.rm = TRUE)
                    cat(sprintf("    Hallmark: %d sig (q<0.05) / %d total\n", sig, nrow(r_h)))
                    all_results[[paste0(run, "_", cn, "_Hallmark")]] <- r_h
                }
            }
        } else {
            cat("    Skipping — too few ranked genes\n")
        }
    }
}

# ----- 5. Aggregate & save -----
if (length(all_results) > 0) {
    combined <- do.call(rbind, lapply(all_results, function(d) {
        d[, c("contrast", "db", "pathway", "Description", "NES", "size",
              "pval", "padj", "leadingEdge")]
    }))
    rownames(combined) <- NULL
    write.table(combined, file.path(PATH_DIR, "all_enrichment_summary.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("\nWrote", file.path(PATH_DIR, "all_enrichment_summary.tsv"),
        "(", nrow(combined), "rows)\n")

    # Excel workbook per run
    for (g in unique(sub("_(GO_BP|Hallmark)$", "", names(all_results)))) {
        wb <- createWorkbook()
        for (k in c("GO_BP", "Hallmark")) {
            full <- paste0(g, "_", k)
            if (full %in% names(all_results)) {
                sheet <- substr(k, 1, 31)
                addWorksheet(wb, sheet)
                writeData(wb, sheet, all_results[[full]])
            }
        }
        saveWorkbook(wb, file.path(PATH_DIR, paste0(g, "_enrichment.xlsx")), overwrite = TRUE)
    }

    saveRDS(all_results, file.path(PATH_DIR, "all_enrichment_results.rds"))

    # Summary count table
    sumtab <- do.call(rbind, lapply(all_results, function(d) {
        data.frame(contrast = d$contrast[1], db = d$db[1],
                   sig_q05 = sum(d$padj < 0.05, na.rm = TRUE),
                   sig_q25 = sum(d$padj < 0.25, na.rm = TRUE),
                   total   = nrow(d))
    }))
    rownames(sumtab) <- NULL
    write.table(sumtab, file.path(PATH_DIR, "enrichment_count_summary.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat("\nEnrichment count summary:\n")
    print(sumtab)
}

cat("\n[", format(Sys.time(), "%FT%T%z"), "] Stage E (fgsea fallback) done.\n")
