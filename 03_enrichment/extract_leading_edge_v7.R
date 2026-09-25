#!/usr/bin/env Rscript
# Extract leading-edge gene sets for v7 figure heatmaps.
# Source: pathway/homeolog_sum_strand_corrected_results.rds (clusterProfiler gseaResult list).
# Output: pathway/leading_edge_v7.tsv (long format) + console summary.

suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

PROJ <- "/home/srussell/projects/def-srussell/nirosha_xeno"
RDS  <- file.path(PROJ, "pathway/homeolog_sum_strand_corrected_results.rds")
OUT  <- file.path(PROJ, "pathway/leading_edge_v7.tsv")

stopifnot(file.exists(RDS))
res <- readRDS(RDS)

# Build Entrez → Xenopus parent-symbol map.
# t13 GSEA pipeline mapped Xenopus parent symbol (lowercase) → human Entrez via
#   toupper(xla_symbol) -> human SYMBOL -> ENTREZID
# So core_enrichment Entrez IDs map back via:
#   ENTREZID -> human SYMBOL (org.Hs.eg.db) -> tolower() -> Xenopus parent symbol.
# Validate against voom_*_homeolog_sum gene_id columns at heatmap time.
hs_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
  keytype = "ENTREZID",
  columns = c("ENTREZID", "SYMBOL")
)
hs_map <- hs_map[!is.na(hs_map$SYMBOL), ]
entrez_to_hs_symbol <- setNames(hs_map$SYMBOL, hs_map$ENTREZID)
to_xla_parent <- function(entrez_ids) {
  sym <- entrez_to_hs_symbol[entrez_ids]
  tolower(sym)
}

cat("Loaded RDS class:", class(res)[1], "\n")
cat("Top-level length:", length(res), "\n")
if (is.list(res) && length(res) > 0) {
  cat("First few names:\n"); print(utils::head(names(res), 6))
}

targets <- data.frame(
  panel = c("fig2_B", "fig3_B", "fig4_B", "fig4_C"),
  contrast_key = c(
    "voom_BLAST_homeolog_sum::CT_vs_ND_7DPA",
    "voom_BRAIN_homeolog_sum::CT_vs_ND_24HPA",
    "voom_BRAIN_homeolog_sum::CT_vs_ND_7DPA",
    "voom_BRAIN_homeolog_sum::CT_vs_ND_7DPA"
  ),
  db = rep("MSigDB_Hallmark", 4),
  pathway_id = c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_ALLOGRAFT_REJECTION",
    "HALLMARK_HEME_METABOLISM"
  ),
  stringsAsFactors = FALSE
)

db_suffix <- function(db) {
  switch(db,
         "MSigDB_Hallmark" = "hallmark",
         "GO_BP"           = "goBP",
         stop("Unknown db: ", db))
}

get_gsea_result <- function(contrast_key, db) {
  # RDS keys look like: voom_BRAIN_homeolog_sum_CT_vs_ND_24HPA_hallmark
  # Map TSV-style "voom_BRAIN_homeolog_sum::CT_vs_ND_24HPA" → that form.
  rds_contrast <- gsub("::", "_", contrast_key, fixed = TRUE)
  full_key <- paste0(rds_contrast, "_", db_suffix(db))
  if (full_key %in% names(res)) return(res[[full_key]])
  hits <- grep(rds_contrast, names(res), fixed = TRUE, value = TRUE)
  for (h in hits) if (grepl(db_suffix(db), h, fixed = TRUE)) return(res[[h]])
  NULL
}

extract_one <- function(panel, contrast_key, db, pathway_id) {
  gres <- get_gsea_result(contrast_key, db)
  if (is.null(gres)) {
    warning(sprintf("[%s] No gseaResult for %s :: %s", panel, contrast_key, db))
    return(NULL)
  }
  df <- as.data.frame(gres)
  row <- df[df$ID == pathway_id, , drop = FALSE]
  if (nrow(row) == 0) {
    warning(sprintf("[%s] Pathway %s missing in %s :: %s", panel, pathway_id, contrast_key, db))
    return(NULL)
  }
  ce <- row$core_enrichment[1]
  entrez_ids <- strsplit(ce, "/", fixed = TRUE)[[1]]
  hs_symbols <- entrez_to_hs_symbol[entrez_ids]
  xla_symbols <- tolower(hs_symbols)
  data.frame(
    panel = panel,
    contrast = contrast_key,
    db = db,
    pathway_id = pathway_id,
    NES = row$NES[1],
    qvalue = row$qvalue[1],
    n_leading = length(entrez_ids),
    le_rank = seq_along(entrez_ids),
    entrez = entrez_ids,
    hs_symbol = unname(hs_symbols),
    xla_symbol = unname(xla_symbols),
    stringsAsFactors = FALSE
  )
}

parts <- vector("list", nrow(targets))
for (i in seq_len(nrow(targets))) {
  parts[[i]] <- extract_one(targets$panel[i], targets$contrast_key[i],
                            targets$db[i], targets$pathway_id[i])
}
out <- do.call(rbind, parts)

if (is.null(out) || nrow(out) == 0) {
  stop("No leading-edge sets extracted; check RDS structure.")
}

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
write.table(out, OUT, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n=== Leading-edge summary ===\n")
summary_df <- aggregate(
  cbind(NES, qvalue, n_leading) ~ panel + contrast + pathway_id,
  data = out,
  FUN = function(x) x[1]
)
print(summary_df)

cat(sprintf("\nWrote %s (%d rows across %d panels)\n",
            OUT, nrow(out), length(unique(out$panel))))
