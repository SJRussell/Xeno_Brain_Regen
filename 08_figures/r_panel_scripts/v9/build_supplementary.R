#!/usr/bin/env Rscript
# Build v9 supplementary data tables:
#  (1) per-homeolog primary DE topTables for 5 contrasts -> TSV
#  (2) Hallmark + GO BP GSEA at pvalueCutoff=1 for 5 contrasts -> TSV
#  (3) Walker per-cell cluster + treatment assignment -> TSV
# WGCNA / MOFA / Liao / 50-Hallmark v9 fix TSVs are gathered as-is by the tar wrapper.
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(readr); library(tibble); library(limma)
  library(clusterProfiler); library(msigdbr); library(org.Hs.eg.db); library(org.Xl.eg.db)
  library(AnnotationDbi)
})

PROJ <- "/project/def-srussell/nirosha_xeno"
OUT  <- file.path(PROJ, "supplementary_v9_staging")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "de_topTables"),        showWarnings = FALSE)
dir.create(file.path(OUT, "gsea_full_pcutoff1"),  showWarnings = FALSE)
dir.create(file.path(OUT, "walker"),              showWarnings = FALSE)

# ---- (1) DE topTables (per-homeolog primary) ----
contrasts_spec <- list(
  BLAST_CT_vs_ND_7DPA       = list(voom = "voom_BLAST_primary.rds", slot = "CT_vs_ND_7DPA"),
  BLAST_interaction_7Dv24H  = list(voom = "voom_BLAST_primary.rds", slot = "interaction_7Dv24H"),
  BRAIN_CT_vs_ND_24HPA      = list(voom = "voom_BRAIN_primary.rds", slot = "CT_vs_ND_24HPA"),
  BRAIN_CT_vs_ND_7DPA       = list(voom = "voom_BRAIN_primary.rds", slot = "CT_vs_ND_7DPA"),
  BRAIN_interaction_7Dv24H  = list(voom = "voom_BRAIN_primary.rds", slot = "interaction_7Dv24H")
)
de_tables <- list()
for (cn in names(contrasts_spec)) {
  s <- contrasts_spec[[cn]]
  v <- readRDS(file.path(PROJ, "deg", s$voom))
  d <- v$results[[s$slot]]
  if (is.null(d)) { cat("MISSING:", cn, "\n"); next }
  out_path <- file.path(OUT, "de_topTables", sprintf("de_%s.tsv", cn))
  write_tsv(as_tibble(d), out_path)
  cat(sprintf("DE TSV: %s  rows=%d  cols=%s\n", cn, nrow(d), paste(colnames(d), collapse=",")))
  de_tables[[cn]] <- list(voom = s$voom, slot = s$slot, rds = v)
}

# ---- (2) GSEA Hallmark + GO BP at pvalueCutoff=1 ----
hm_hs <- msigdbr(species = "Homo sapiens", collection = "H")
entrez_col <- intersect(c("ncbi_gene","human_entrez_gene","entrez_gene"), colnames(hm_hs))[1]
hm_t2g <- unique(as.data.frame(hm_hs[, c("gs_name", entrez_col)]))
colnames(hm_t2g) <- c("gs_name", "entrez")
hm_t2g <- hm_t2g[!is.na(hm_t2g$entrez), ]
cat("Hallmark gene sets:", length(unique(hm_t2g$gs_name)), "\n")

sym_hs <- keys(org.Hs.eg.db, keytype = "SYMBOL")
sym_to_hs <- AnnotationDbi::select(org.Hs.eg.db, keys = sym_hs, keytype = "SYMBOL",
                                    columns = c("SYMBOL","ENTREZID"))
hs_lookup <- setNames(sym_to_hs$ENTREZID, toupper(sym_to_hs$SYMBOL))
strip_xla <- function(s) sub("\\.[LS]$", "", s)
xla_to_hs <- function(xla_symbols) hs_lookup[toupper(strip_xla(xla_symbols))]

run_full_hallmark <- function(d, cn) {
  rk <- sign(d$logFC) * (-log10(d$P.Value))
  names(rk) <- d$gene_id
  rk <- rk[is.finite(rk)]
  rk <- sort(rk, decreasing = TRUE)
  hs_e <- xla_to_hs(names(rk))
  keep <- !is.na(hs_e) & !duplicated(hs_e)
  rk_h <- rk[keep]; names(rk_h) <- hs_e[keep]
  rk_h <- sort(rk_h, decreasing = TRUE)
  set.seed(42)
  res <- GSEA(rk_h, TERM2GENE = hm_t2g, pvalueCutoff = 1, minGSSize = 5, maxGSSize = 5000,
              eps = 0, verbose = FALSE)
  out <- as.data.frame(res@result) %>%
    transmute(term = ID, NES, pvalue, padj = p.adjust, qvalue,
              core_enrichment = core_enrichment)
  out_path <- file.path(OUT, "gsea_full_pcutoff1", sprintf("hallmark_%s_full.tsv", cn))
  write_tsv(out, out_path)
  cat(sprintf("Hallmark TSV: %s  rows=%d\n", cn, nrow(out)))
}
run_full_gobp <- function(d, cn) {
  xla_e <- AnnotationDbi::select(org.Xl.eg.db, keys = d$gene_id, keytype = "SYMBOL",
                                  columns = c("SYMBOL","ENTREZID"))
  xla_e <- xla_e[!is.na(xla_e$ENTREZID) & !duplicated(xla_e$SYMBOL), ]
  rk <- sign(d$logFC) * (-log10(d$P.Value))
  names(rk) <- d$gene_id
  rk_e <- rk[match(xla_e$SYMBOL, names(rk))]
  names(rk_e) <- xla_e$ENTREZID
  rk_e <- rk_e[is.finite(rk_e)]
  rk_e <- sort(rk_e, decreasing = TRUE)
  set.seed(42)
  res <- tryCatch(
    gseGO(rk_e, OrgDb = org.Xl.eg.db, keyType = "ENTREZID", ont = "BP",
          pvalueCutoff = 1, minGSSize = 5, maxGSSize = 5000, eps = 0, verbose = FALSE),
    error = function(e) { cat("gseGO error:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) return(invisible(NULL))
  out <- as.data.frame(res@result) %>%
    transmute(term = ID, Description, NES, pvalue, padj = p.adjust, qvalue,
              core_enrichment = core_enrichment)
  out_path <- file.path(OUT, "gsea_full_pcutoff1", sprintf("gobp_%s_full.tsv", cn))
  write_tsv(out, out_path)
  cat(sprintf("GO BP TSV: %s  rows=%d\n", cn, nrow(out)))
}

for (cn in names(contrasts_spec)) {
  rds <- de_tables[[cn]]$rds
  d   <- rds$results[[ contrasts_spec[[cn]]$slot ]]
  if (is.null(d) || nrow(d) == 0) { cat("skip:", cn, "\n"); next }
  cat("=== GSEA full pvalueCutoff=1 for", cn, "===\n")
  run_full_hallmark(d, cn)
  run_full_gobp(d, cn)
}

# ---- (3) Walker per-cell cluster + treatment assignment ----
cat("Loading Walker Seurat object ...\n")
seur <- readRDS(file.path(PROJ, "comparators/walker/walker_integrated_v2.rds"))
wcells <- tibble(
  cell      = colnames(seur),
  cluster   = as.integer(as.character(seur$seurat_clusters)),
  treatment = as.character(seur$treatment)
)
write_tsv(wcells, file.path(OUT, "walker", "walker_per_cell_cluster_treatment.tsv"))
cat(sprintf("Walker per-cell TSV: rows=%d\n", nrow(wcells)))

cat("DONE.\n")
