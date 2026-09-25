#!/usr/bin/env Rscript
# Stage E: enrichment analysis (GSEA + ORA)
# - Per-contrast ranked GSEA via clusterProfiler (GO BP, KEGG, MSigDB Hallmark)
# - ORA on q<0.05 DEG sets where they exist (BLAST mostly)
# - Xenopus → human ortholog mapping for MSigDB/Reactome
# - Outputs: per-contrast {enrichGO, gseGO, gseKEGG, hallmark} RDS + xlsx

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Xl.eg.db)
    library(AnnotationDbi)
    library(msigdbr)
    library(openxlsx)
})

PROJ <- "/project/def-srussell/nirosha_xeno"
DEG_DIR <- file.path(PROJ, "deg")
PATH_DIR <- file.path(PROJ, "pathway")
FIG_DIR  <- file.path(PROJ, "figures", "stage_E")
dir.create(PATH_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("[", format(Sys.time(), "%FT%T%z"), "] Loading DE results...\n")

# ----- Symbol → Entrez mapping for X. laevis -----
# org.Xl.eg.db keys by Entrez Gene ID. Our gene_ids from STAR/GFF3 are NCBI symbols
# (e.g. dok1.L, gdf15.L, LOC108704873).
symbol_to_entrez <- function(symbols) {
    m <- AnnotationDbi::select(org.Xl.eg.db,
                               keys = symbols, keytype = "SYMBOL",
                               columns = c("SYMBOL", "ENTREZID"),
                               multiVals = "first")
    setNames(m$ENTREZID, m$SYMBOL)
}

# ----- MSigDB Hallmark for Xenopus via human ortholog mapping -----
# Need human ortholog Entrez per Xenopus gene. Use msigdbr's human Hallmark set + map.
hallmark_human <- msigdbr(species = "Homo sapiens", collection = "H")
# msigdbr 26.x uses ncbi_gene (formerly human_entrez_gene)
entrez_col <- intersect(c("ncbi_gene", "human_entrez_gene", "entrez_gene"),
                        colnames(hallmark_human))[1]
cat("Using msigdbr Entrez column:", entrez_col, "\n")
hallmark_h <- as.data.frame(hallmark_human[, c("gs_name", entrez_col)])
colnames(hallmark_h) <- c("gs_name", "entrez")
hallmark_h <- unique(hallmark_h)
hallmark_h <- hallmark_h[!is.na(hallmark_h$entrez), ]
hallmark_h_split <- split(hallmark_h$entrez, hallmark_h$gs_name)
cat("Hallmark gene sets:", length(hallmark_h_split), "\n")

# Map Xenopus Entrez → human Entrez via org.Xl.eg.db -> SYMBOL -> human (case-insensitive symbol match)
# This is a simplified ortholog approach; for high-confidence ortholog mapping use Ensembl Compara.
# For the manuscript we'll cite Xenbase + use HGNC-aligned symbols which match between Xla and human.

# Strategy: rank by symbol in upper case, use Hallmark gene sets keyed on human Entrez.
# Build a Xenopus_Entrez -> human_Entrez map via SYMBOL.
# For symbols with Xla-specific suffixes (.L, .S), strip them.
strip_xla_suffix <- function(s) sub("\\.[LS]$", "", s)
sym_xla <- keys(org.Xl.eg.db, keytype = "SYMBOL")
sym_to_xla_entrez <- AnnotationDbi::select(org.Xl.eg.db,
                                            keys = sym_xla, keytype = "SYMBOL",
                                            columns = c("SYMBOL", "ENTREZID"))

# Use human gene info — try org.Hs.eg.db
suppressWarnings(library(org.Hs.eg.db))
sym_hs <- keys(org.Hs.eg.db, keytype = "SYMBOL")
sym_to_hs_entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                           keys = sym_hs, keytype = "SYMBOL",
                                           columns = c("SYMBOL", "ENTREZID"))
hs_lookup <- setNames(sym_to_hs_entrez$ENTREZID,
                      toupper(sym_to_hs_entrez$SYMBOL))

xla_to_hs_entrez <- function(xla_symbols) {
    base_sym <- toupper(strip_xla_suffix(xla_symbols))
    hs_lookup[base_sym]
}

# ----- Helper: build ranked vector from a topTable result -----
build_rank <- function(res_df) {
    rk <- sign(res_df$logFC) * (-log10(res_df$P.Value))
    names(rk) <- res_df$gene_id
    rk <- rk[!is.na(rk) & is.finite(rk)]
    sort(rk, decreasing = TRUE)
}

# ----- Run GSEA against MSigDB Hallmark (human ortholog) -----
run_gsea_hallmark <- function(rank_vec, contrast_name) {
    hs_entrez <- xla_to_hs_entrez(names(rank_vec))
    keep <- !is.na(hs_entrez) & !duplicated(hs_entrez)
    rk_h <- rank_vec[keep]
    names(rk_h) <- hs_entrez[keep]
    rk_h <- sort(rk_h, decreasing = TRUE)

    if (length(rk_h) < 100) {
        cat("  Skipping Hallmark GSEA for ", contrast_name, ": ", length(rk_h),
            " mappable genes\n")
        return(NULL)
    }

    set.seed(42)
    # clusterProfiler::GSEA expects TERM2GENE as data.frame(term, gene)
    term2gene <- hallmark_h[, c("gs_name", "entrez")]
    res <- tryCatch({
        GSEA(rk_h, TERM2GENE = term2gene, pvalueCutoff = 0.25,
             minGSSize = 15, maxGSSize = 500, eps = 0)
    }, error = function(e) {cat("    GSEA error:", conditionMessage(e), "\n"); NULL})

    if (is.null(res) || nrow(res@result) == 0) {
        cat("  Hallmark GSEA: 0 enriched (", contrast_name, ")\n")
        return(NULL)
    }

    df <- res@result[, c("ID", "NES", "pvalue", "p.adjust", "qvalue", "core_enrichment")]
    df$contrast <- contrast_name
    df$db <- "MSigDB_Hallmark"
    cat("  Hallmark GSEA: ", nrow(df), " enriched (q<0.25) for ", contrast_name, "\n")
    df
}

# ----- Run GO BP GSEA via clusterProfiler -----
run_gsea_GO <- function(rank_vec, contrast_name) {
    # Convert symbol → Xla Entrez
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
        cat("  Skipping GO GSEA for ", contrast_name, ": ", length(rk_e),
            " mappable genes\n")
        return(NULL)
    }

    set.seed(42)
    res <- tryCatch({
        gseGO(rk_e, OrgDb = org.Xl.eg.db, keyType = "ENTREZID",
              ont = "BP", pvalueCutoff = 0.25,
              minGSSize = 15, maxGSSize = 500, eps = 0,
              verbose = FALSE)
    }, error = function(e) NULL)

    if (is.null(res) || nrow(res@result) == 0) {
        cat("  GO BP GSEA: 0 enriched (", contrast_name, ")\n")
        return(NULL)
    }

    df <- res@result[, c("ID", "Description", "NES", "pvalue", "p.adjust", "qvalue", "core_enrichment")]
    df$contrast <- contrast_name
    df$db <- "GO_BP"
    cat("  GO BP GSEA: ", nrow(df), " enriched (q<0.25) for ", contrast_name, "\n")
    df
}

# ----- Driver: per-tissue, per-contrast -----
all_enrich <- list()

runs_to_enrich <- c(
    "deseq2_BRAIN_primary",   # PRIMARY (per SO-02, 2026-05-16)
    "deseq2_BLAST_primary",
    "voom_BRAIN_primary",     # confirmatory
    "voom_BLAST_primary"
)
for (run in runs_to_enrich) {
    cat("\n==========================================\n")
    cat("Run: ", run, "\n")
    cat("==========================================\n")
    rds <- file.path(DEG_DIR, paste0(run, ".rds"))
    obj <- readRDS(rds)
    for (cn in names(obj$results)) {
        res_df <- obj$results[[cn]]
        rank_vec <- build_rank(res_df)
        cat("\n--- ", run, "::", cn, " (n=", length(rank_vec), " ranked) ---\n", sep = "")

        # Hallmark
        h_res <- run_gsea_hallmark(rank_vec, paste0(run, "::", cn))
        if (!is.null(h_res)) all_enrich[[paste0(run, "_", cn, "_hallmark")]] <- h_res

        # GO BP
        go_res <- run_gsea_GO(rank_vec, paste0(run, "::", cn))
        if (!is.null(go_res)) all_enrich[[paste0(run, "_", cn, "_goBP")]] <- go_res
    }
}

# ----- Combine + write -----
# Ensure Description column exists (Hallmark uses ID as description)
all_enrich <- lapply(all_enrich, function(d) {
    if (!"Description" %in% names(d)) d$Description <- d$ID
    d
})
combined_enrich <- do.call(rbind, lapply(all_enrich, function(d) {
    d[, c("contrast", "db", "ID", "Description", "NES", "pvalue", "p.adjust", "qvalue")]
}))
write.table(combined_enrich, file.path(PATH_DIR, "all_enrichment_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nWrote ", PATH_DIR, "/all_enrichment_summary.tsv (", nrow(combined_enrich), " rows)\n", sep = "")

# Save per-run enrichment as Excel workbooks for downstream figures
xl_groups <- split(seq_along(all_enrich), sub("(_hallmark|_goBP)$", "", names(all_enrich)))
for (g in names(xl_groups)) {
    wb <- createWorkbook()
    for (i in xl_groups[[g]]) {
        nm <- names(all_enrich)[i]
        sheet <- substr(sub(paste0("^", g, "_"), "", nm), 1, 31)
        addWorksheet(wb, sheet)
        writeData(wb, sheet, all_enrich[[i]])
    }
    saveWorkbook(wb, file.path(PATH_DIR, paste0(g, "_enrichment.xlsx")), overwrite = TRUE)
}

saveRDS(all_enrich, file.path(PATH_DIR, "all_enrichment_results.rds"))

cat("\n[", format(Sys.time(), "%FT%T%z"), "] Stage E done.\n")
