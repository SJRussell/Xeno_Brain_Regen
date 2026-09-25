#!/usr/bin/env Rscript
# RED-2 (iter-19): re-run UCell scoring on walker_integrated.rds with corrected
# axolotl marker IDs.
#
# Diagnosis:
#   - tx2gene.tsv extracted gene_name from GTF only when the [hs] homolog tag was
#     present. Walker's signature genes (etv1, slc17a7, gad1, gfap, etc.) are
#     direct axolotl annotations in the GTF (gene_name "ETV1|AMEX60DD..."), with
#     no [hs] tag — so tx2gene wrote the gene_id as the third-col value.
#   - In the Seurat object, those genes appear under rowname AMEX60DD<gene_id>,
#     not under "ETV1" / "GAD1" / etc. Passing the toupper symbol to UCell
#     dropped them silently → glut_etv1plus_regen UCell score ~0.0044 across
#     all clusters.
#
# Fix:
#   - Build marker lists in the AMEX60DD gene-id namespace using direct GTF lookup.
#   - Where Walker's marker is absent from AmexT_v47 (e.g. SLC17A7=vGlut1), fall
#     back to closest paralog (SLC17A6=vGlut2, SLC17A8=vGlut3) and document.
#   - Re-score, re-derive cluster_annotations.tsv, and re-run treatment-population
#     enrichment with the new etv1+ cluster.

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(data.table)
    library(UCell)
    library(future)
})
set.seed(20260502)
options(future.globals.maxSize = 8 * 1024^3)
plan("sequential")

PROJ <- "/project/def-srussell/nirosha_xeno"
WORK <- "/scratch/srussell/srussell/nirosha_xeno_comparators/walker/work"
PROJ_OUT <- file.path(PROJ, "comparators", "walker")

cat("[", format(Sys.time()), "] loading walker_integrated.rds (1.2 GB) ...\n")
integrated <- readRDS(file.path(PROJ_OUT, "walker_integrated.rds"))
cat("  loaded:", ncol(integrated), "cells x", nrow(integrated), "features\n")
cat("  treatments:", paste(unique(integrated$treatment), collapse = ", "), "\n")
cat("  clusters:", paste(sort(unique(as.character(integrated$seurat_clusters))), collapse = ", "), "\n")

DefaultAssay(integrated) <- "SCT"
rn <- rownames(integrated)
cat("  rownames sample:", head(rn, 5), "...\n")

# ----------------------------------------------------------------------
# Diagnostic: per-marker non-zero-cell percentage with the OLD names
# ----------------------------------------------------------------------
walker_old <- c("ETV1","SLC17A7","CDH13","BCL11B","SOX5","NETO1","STK39","NTNG1",
                "GFAP","SOX2","GAD1","GAD2")
cat("\nDIAGNOSTIC: which Walker markers (passed as upper-case symbol) are present?\n")
for (m in walker_old) cat(sprintf("  %-12s present=%s\n", m, m %in% rn))

# Per-marker non-zero-cell %
present_old <- intersect(walker_old, rn)
if (length(present_old) > 0) {
    fd <- FetchData(integrated, vars = present_old, layer = "counts")
    pct_nz_old <- sapply(fd, function(x) round(100 * mean(x > 0), 2))
    cat("\nDIAGNOSTIC: per-marker non-zero-cell % (using OLD symbol-based names):\n")
    print(pct_nz_old)
}

# ----------------------------------------------------------------------
# Corrected marker mapping (from GTF inspection 2026-05-02)
# ----------------------------------------------------------------------
# Each marker maps to one or more rownames in the Seurat object. For genes
# with both direct ("GENE|AMEX...") and [hs] ("X [nr]|GENE [hs]|AMEX...")
# annotations we include both, since they may correspond to different
# paralogs/homeologs in the axolotl genome.
markers_amex <- list(
    glut_etv1plus_regen = c(
        "AMEX60DD022427",                  # ETV1 (direct only)
        "AMEX60DD004989", "AMEX60DD016391", # SLC17A6 = vGlut2 (direct + [hs]) — vGlut1=SLC17A7 absent
        "AMEX60DD007749",                  # SLC17A8 = vGlut3 (direct only)
        "AMEX60DD017526", "CDH13",          # CDH13 (direct + [hs] symbol)
        "AMEX60DD011386", "BCL11B",         # BCL11B
        "AMEX60DD006750",                  # SOX5 (direct only)
        "AMEX60DD039135", "NETO1",          # NETO1
        "STK39",                            # STK39 (same gene_id direct + [hs])
        "AMEX60DD019201"                   # NTNG1 (direct only)
    ),
    glutamatergic_broad = c(
        # Pan-glutamatergic — what slc17a7 was supposed to score; use any vGlut.
        "AMEX60DD004989", "AMEX60DD016391", "AMEX60DD007749"  # SLC17A6 + SLC17A8
    ),
    gabaergic = c(
        "AMEX60DD055654", "AMEX60DD038559", "GAD1",           # GAD1 direct + [hs]+ symbol
        "AMEX60DD021878", "AMEX60DD021875", "GAD2",           # GAD2 direct + [hs] + symbol
        "AMEX60DD028674",                                    # SLC32A1 (direct only)
        "AMEX60DD023994"                                     # SLC6A1 ([hs] only)
    ),
    ependymoglia = c(
        "AMEX60DD009861",                                    # GFAP (direct only)
        "AMEX60DD001640", "AMEX60DD012772", "SOX2",           # SOX2 direct + [hs] + symbol
        "AMEX60DD034100", "AMEX60DD005433", "FABP7"           # FABP7 direct + [hs] + symbol
    ),
    radial_glia = c(
        "AMEX60DD016851",                                    # NES (direct)
        "AMEX60DD022108", "AMEX60DD033837", "VIM",            # VIM direct + [hs] + symbol
        "AMEX60DD009861",                                    # GFAP
        "AMEX60DD001640", "AMEX60DD012772", "SOX2"            # SOX2
    ),
    oligodendrocyte = c(
        "AMEX60DD047308",                                    # OLIG1 (direct only)
        "AMEX60DD047306", "OLIG2",                            # OLIG2 ([hs] only) + symbol
        "AMEX60DD029302",                                    # SOX10 (direct only)
        "AMEX60DD039196", "MBP"                              # MBP direct + symbol
    ),
    microglia_macro = c(
        "AMEX60DD012740",                                    # CD68 (direct)
        "AMEX60DD010609",                                    # AIF1 (direct)
        "AMEX60DD030145",                                    # CSF1R (direct)
        "AMEX60DD001912",                                    # P2RY12 (direct)
        "AMEX60DD028004", "ITGAM"                             # ITGAM ([hs])
    ),
    neuroendocrine = c(
        "AMEX60DD046200", "AVP",                              # AVP ([hs])
        "AMEX60DD046199", "OXT",                              # OXT ([hs])
        "AMEX60DD039827",                                    # CRH (direct)
        "AMEX60DD032740"                                     # POMC (direct)
    )
)

# Validate: which entries actually exist as rownames?
cat("\nMarker presence in Seurat rownames (after correction):\n")
for (set_name in names(markers_amex)) {
    ms <- markers_amex[[set_name]]
    present <- ms %in% rn
    cat(sprintf("  %-22s n_total=%d  present=%d  missing=%d\n",
                set_name, length(ms), sum(present), sum(!present)))
    if (any(!present)) {
        cat("    missing entries (unexpected — please review):", paste(ms[!present], collapse = ", "), "\n")
    }
}

# Drop any missing entries before scoring (UCell handles missing but we want clean accounting)
markers_amex <- lapply(markers_amex, function(ms) ms[ms %in% rn])

# Per-marker non-zero-cell % for the corrected set
all_corrected <- unique(unlist(markers_amex))
cat("\nDIAGNOSTIC: per-marker non-zero-cell % (corrected gene-id namespace):\n")
fd_new <- FetchData(integrated, vars = all_corrected, layer = "counts")
pct_nz_new <- round(100 * colMeans(fd_new > 0), 2)
print(pct_nz_new)

# ----------------------------------------------------------------------
# Re-score with UCell using corrected markers
# ----------------------------------------------------------------------
cat("\nRunning AddModuleScore_UCell with corrected markers ...\n")
# Drop any prior _UCell columns to avoid collision
old_cols <- grep("_UCell$", colnames(integrated@meta.data), value = TRUE)
if (length(old_cols) > 0) {
    cat("  removing", length(old_cols), "stale _UCell columns from prior scoring\n")
    integrated@meta.data[, old_cols] <- NULL
}
integrated <- AddModuleScore_UCell(integrated, features = markers_amex, ncores = 4,
                                    name = "_UCellv2")

# ----------------------------------------------------------------------
# Cluster-level mean scores → re-derive cluster_annotations.tsv
# ----------------------------------------------------------------------
score_cols <- paste0(names(markers_amex), "_UCellv2")
cluster_scores <- sapply(names(markers_amex), function(set_name) {
    score_col <- paste0(set_name, "_UCellv2")
    tapply(integrated[[score_col, drop = TRUE]],
           integrated$seurat_clusters, mean, na.rm = TRUE)
})
cluster_anno <- data.frame(
    cluster = rownames(cluster_scores),
    cluster_scores,
    candidate_label = colnames(cluster_scores)[apply(cluster_scores, 1, which.max)],
    candidate_score = apply(cluster_scores, 1, max)
)
cat("\nCorrected cluster annotations:\n")
print(cluster_anno)
fwrite(cluster_anno, file.path(WORK, "cluster_annotations_v2.tsv"), sep = "\t")
file.copy(file.path(WORK, "cluster_annotations_v2.tsv"),
          file.path(PROJ_OUT, "cluster_annotations_v2.tsv"), overwrite = TRUE)

# Per-cell scores by treatment
per_cell <- as.data.table(integrated@meta.data[, c("sample","treatment","seurat_clusters", score_cols)],
                          keep.rownames = "cell")
fwrite(per_cell, file.path(WORK, "ucell_signatures_per_cell_v2.tsv"), sep = "\t")
file.copy(file.path(WORK, "ucell_signatures_per_cell_v2.tsv"),
          file.path(PROJ_OUT, "ucell_signatures_per_cell_v2.tsv"), overwrite = TRUE)

per_cluster_treat <- per_cell[, lapply(.SD, mean, na.rm = TRUE),
                              by = .(seurat_clusters, treatment),
                              .SDcols = score_cols]
fwrite(per_cluster_treat, file.path(WORK, "ucell_signatures_per_cluster_v2.tsv"), sep = "\t")
file.copy(file.path(WORK, "ucell_signatures_per_cluster_v2.tsv"),
          file.path(PROJ_OUT, "ucell_signatures_per_cluster_v2.tsv"), overwrite = TRUE)

# ----------------------------------------------------------------------
# Re-run treatment-population enrichment for new etv1+ cluster
# ----------------------------------------------------------------------
etv1_cluster_v2 <- cluster_anno$cluster[which.max(cluster_anno$glut_etv1plus_regen)]
cat("\netv1+ regen-required cluster (corrected best-match):", etv1_cluster_v2, "\n")
cat("  glut_etv1plus_regen score for this cluster:",
    cluster_anno$glut_etv1plus_regen[which.max(cluster_anno$glut_etv1plus_regen)], "\n")

per_cell[, is_etv1_pop := seurat_clusters == etv1_cluster_v2]
prop_dt <- per_cell[, .(n_etv1 = sum(is_etv1_pop), n_total = .N,
                         prop = sum(is_etv1_pop) / .N), by = treatment]
fwrite(prop_dt, file.path(WORK, "etv1_population_proportions_v2.tsv"), sep = "\t")
file.copy(file.path(WORK, "etv1_population_proportions_v2.tsv"),
          file.path(PROJ_OUT, "etv1_population_proportions_v2.tsv"), overwrite = TRUE)

fisher_test <- function(t1, t2) {
    a <- per_cell[treatment == t1, sum(is_etv1_pop)]
    b <- per_cell[treatment == t1, sum(!is_etv1_pop)]
    c <- per_cell[treatment == t2, sum(is_etv1_pop)]
    d <- per_cell[treatment == t2, sum(!is_etv1_pop)]
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
    data.table(comparison = paste(t1, "vs", t2),
               OR = ft$estimate, p = ft$p.value,
               lo95 = ft$conf.int[1], hi95 = ft$conf.int[2])
}
treat_enrich <- rbind(
    fisher_test("tail", "uninjured"),
    fisher_test("limb", "uninjured"),
    fisher_test("tail", "limb")
)
fwrite(treat_enrich, file.path(WORK, "treatment_pop_enrichment_v2.tsv"), sep = "\t")
file.copy(file.path(WORK, "treatment_pop_enrichment_v2.tsv"),
          file.path(PROJ_OUT, "treatment_pop_enrichment_v2.tsv"), overwrite = TRUE)
cat("\nCorrected treatment-population enrichment:\n"); print(treat_enrich)

# ----------------------------------------------------------------------
# Save the re-scored Seurat object (preserves prior _UCell columns dropped above)
# ----------------------------------------------------------------------
saveRDS(integrated, file.path(PROJ_OUT, "walker_integrated_v2.rds"))

cat("\n[", format(Sys.time(), "%FT%T%z"), "] walker_marker_remap.R complete.\n", sep = "")
