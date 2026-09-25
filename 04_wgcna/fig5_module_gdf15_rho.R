#!/usr/bin/env Rscript
# P1: re-derive Fig 5C and Fig 5D module-GDF15 rhos cleanly.
#
# Inputs:
#   /project/def-srussell/nirosha_xeno/wgcna/wgcna_BRAIN.rds
#   /project/def-srussell/nirosha_xeno/wgcna/module_assignments_BRAIN.tsv
#   /project/def-srussell/nirosha_xeno/deg/vsd_blind.rds
#   /project/def-srussell/nirosha_xeno/deg/sample_metadata.tsv
#
# Outputs:
#   /project/def-srussell/nirosha_xeno/wgcna/brain_module_hallmark_enrichment.tsv
#   /project/def-srussell/nirosha_xeno/wgcna/brain_module_eigengene_trait_lm.tsv
#   /project/def-srussell/nirosha_xeno/wgcna/brain_module_gdf15L_correlations.tsv
#   /project/def-srussell/nirosha_xeno/figures/fig5/C_panel.pdf  (+ .png + metadata)
#   /project/def-srussell/nirosha_xeno/figures/fig5/D_panel.pdf  (+ .png + metadata)
#   /project/def-srussell/nirosha_xeno/figures/fig5/p1_log.txt

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(clusterProfiler)
  library(msigdbr)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(patchwork)
  library(yaml)
})

PROJ <- "/project/def-srussell/nirosha_xeno"
WGCNA_RDS <- file.path(PROJ, "wgcna/wgcna_BRAIN.rds")
MOD_TSV   <- file.path(PROJ, "wgcna/module_assignments_BRAIN.tsv")
VSD_RDS   <- file.path(PROJ, "deg/vsd_blind.rds")
META_TSV  <- file.path(PROJ, "deg/sample_metadata.tsv")
FIG5_DIR  <- file.path(PROJ, "figures/fig5")
WG_OUT    <- file.path(PROJ, "wgcna")
dir.create(FIG5_DIR, showWarnings = FALSE, recursive = TRUE)

TSTAMP <- format(Sys.time(), "%Y%m%dT%H%M%S")
LOG_PATH <- file.path(FIG5_DIR, "p1_log.txt")
log_con <- file(LOG_PATH, open = "w")
log <- function(...) { msg <- paste0(..., collapse = ""); cat(msg, "\n"); cat(msg, "\n", file = log_con, append = TRUE) }
log("== P1 Fig 5C/5D rho re-derivation — ", TSTAMP, " ==")

# ---- Preserve prior C_panel.pdf and D_panel.pdf per SO-12 ----
for (p in c("C_panel.pdf", "D_panel.pdf", "C_panel.png", "D_panel.png", "panel_metadata_C.yaml", "panel_metadata_D.yaml")) {
  src <- file.path(FIG5_DIR, p)
  if (file.exists(src)) {
    ext <- tools::file_ext(p); base <- tools::file_path_sans_ext(p)
    dst <- file.path(FIG5_DIR, paste0(base, "_", TSTAMP, ".", ext))
    file.rename(src, dst); log("SO-12 preserved: ", p, " -> ", basename(dst))
  }
}

# ---- Load WGCNA + module assignments ----
wg <- readRDS(WGCNA_RDS)
log("wgcna_BRAIN.rds top-level slots: ", paste(names(wg), collapse=","))

MEs <- NULL
if (!is.null(wg$MEs))      MEs <- wg$MEs
if (is.null(MEs) && !is.null(wg$net) && !is.null(wg$net$MEs)) MEs <- wg$net$MEs
if (is.null(MEs)) stop("Could not locate moduleEigengenes table in wgcna_BRAIN.rds")
log("MEs dim: ", paste(dim(MEs), collapse=" x "))
log("ME columns: ", paste(colnames(MEs), collapse=","))

mod_assign <- read_tsv(MOD_TSV, show_col_types = FALSE)
log("module_assignments rows=", nrow(mod_assign), " unique colors=", length(unique(mod_assign$module_color)))

# Label <-> color crosswalk: MEs are label-keyed (ME0, ME17...), pathway picks are color-keyed.
label_color <- mod_assign %>% distinct(module_label, module_color)
color_to_label <- setNames(as.character(label_color$module_label), label_color$module_color)
label_to_color <- setNames(label_color$module_color, as.character(label_color$module_label))
log("label_color crosswalk built: ", nrow(label_color), " modules")

# ---- Load VST + metadata, subset to 12 brain samples ----
vsd <- readRDS(VSD_RDS)
meta <- read_tsv(META_TSV, show_col_types = FALSE)
brain_meta <- meta %>% filter(tissue == "BRAIN")
log("brain samples in metadata: n=", nrow(brain_meta))
stopifnot(nrow(brain_meta) == 12)

vst_mat <- assay(vsd)
brain_ids <- brain_meta$sample
brain_in_vsd <- intersect(brain_ids, colnames(vst_mat))
log("brain samples found in vsd_blind: ", length(brain_in_vsd), " / 12")
stopifnot(length(brain_in_vsd) == 12)
vst_brain <- vst_mat[, brain_ids]

# ---- Match MEs row order to brain_ids ----
me_rn <- rownames(MEs)
if (length(intersect(me_rn, brain_ids)) < 12) {
  alt <- intersect(sub("_S\\d+$","", me_rn), sub("_S\\d+$","", brain_ids))
  log("Direct sample-ID match insufficient (intersect=", length(intersect(me_rn, brain_ids)),
      "); attempting stripped-suffix match")
  stopifnot(length(alt) >= 12)
  rownames(MEs) <- sub("_S\\d+$","", me_rn)
  brain_ids_stripped <- sub("_S\\d+$","", brain_ids)
  MEs <- MEs[brain_ids_stripped, , drop = FALSE]
  rownames(MEs) <- brain_ids
} else {
  MEs <- MEs[brain_ids, , drop = FALSE]
}
log("MEs aligned to 12 brain samples in order: ", paste(rownames(MEs), collapse=", "))

# ---- gdf15.L VST per sample ----
g_rn <- rownames(vst_brain)
gdf15_candidates <- g_rn[grepl("^gdf15\\.L$", g_rn, ignore.case = TRUE)]
if (length(gdf15_candidates) == 0) {
  gdf15_candidates <- g_rn[grepl("gdf15", g_rn, ignore.case = TRUE)]
  log("Strict gdf15.L not found; loose matches: ", paste(gdf15_candidates, collapse=", "))
}
stopifnot(length(gdf15_candidates) >= 1)
gdf15L_id <- gdf15_candidates[grepl("\\.L$", gdf15_candidates, ignore.case = TRUE)][1]
if (is.na(gdf15L_id)) gdf15L_id <- gdf15_candidates[1]
log("Using gdf15 ID: ", gdf15L_id)
gdf15_vst <- as.numeric(vst_brain[gdf15L_id, ]); names(gdf15_vst) <- colnames(vst_brain)
log("gdf15.L VST per sample:")
for (s in names(gdf15_vst)) log("  ", s, "  ", sprintf("%.4f", gdf15_vst[s]))

# ---- ME-trait linear models (treatment × time interaction) ----
brain_meta <- brain_meta %>% mutate(
  treatment = factor(treatment, levels = c("ND", "CT")),
  time      = factor(time,      levels = c("24HPA", "7DPA"))
) %>% arrange(match(sample, rownames(MEs)))
stopifnot(all(brain_meta$sample == rownames(MEs)))

me_traits <- lapply(colnames(MEs), function(mc) {
  fit <- lm(MEs[[mc]] ~ treatment * time, data = brain_meta)
  s <- summary(fit)$coefficients
  trt  <- if ("treatmentCT" %in% rownames(s))        s["treatmentCT",         c("Estimate","Pr(>|t|)")] else c(NA,NA)
  tm   <- if ("time7DPA" %in% rownames(s))           s["time7DPA",            c("Estimate","Pr(>|t|)")] else c(NA,NA)
  intc <- if ("treatmentCT:time7DPA" %in% rownames(s)) s["treatmentCT:time7DPA", c("Estimate","Pr(>|t|)")] else c(NA,NA)
  lab <- sub("^ME","", mc)
  tibble(module = mc,
         module_label = lab,
         module_color = unname(label_to_color[lab]),
         beta_treatment = trt[1], p_treatment = trt[2],
         beta_time = tm[1], p_time = tm[2],
         beta_interaction = intc[1], p_interaction = intc[2])
}) %>% bind_rows()
write_tsv(me_traits, file.path(WG_OUT, "brain_module_eigengene_trait_lm.tsv"))
log("Wrote brain_module_eigengene_trait_lm.tsv (n_modules=", nrow(me_traits), ")")

# ---- Per-module Hallmark over-representation ----
hm <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
hm_t2g <- hm %>% distinct(gs_name, gene_symbol) %>% dplyr::select(gs_name, gene_symbol)

# Build module gene lists keyed by module_color (matches ME naming "ME<color>")
mod_groups <- split(mod_assign$gene_id, mod_assign$module_color)
universe <- unique(mod_assign$gene_id)
# Strip homeolog suffix and uppercase for Hallmark symbol matching
norm_sym <- function(x) toupper(sub("\\.[LSls]$", "", x))
universe_sym <- unique(norm_sym(universe))

enrich_rows <- list()
for (col in names(mod_groups)) {
  module_genes <- unique(norm_sym(mod_groups[[col]]))
  if (length(module_genes) < 5) next
  res <- tryCatch(
    clusterProfiler::enricher(gene = module_genes, TERM2GENE = hm_t2g, universe = universe_sym,
                              pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 5, maxGSSize = 5000),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0) next
  df <- as.data.frame(res)
  df$module_color <- col
  df$ME_id <- paste0("ME", col)
  enrich_rows[[col]] <- df
}
enr_df <- bind_rows(enrich_rows)
write_tsv(enr_df, file.path(WG_OUT, "brain_module_hallmark_enrichment.tsv"))
log("Wrote brain_module_hallmark_enrichment.tsv (rows=", nrow(enr_df), ")")

pick_min_padj <- function(df, term) {
  sub <- df %>% filter(ID == term) %>% arrange(p.adjust)
  if (nrow(sub) == 0) return(NULL); sub
}

# ---- Pick anchor modules ----
oxphos_rank <- pick_min_padj(enr_df, "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
fa_rank     <- pick_min_padj(enr_df, "HALLMARK_FATTY_ACID_METABOLISM")
mtorc1_rank <- pick_min_padj(enr_df, "HALLMARK_MTORC1_SIGNALING")
emt_rank    <- pick_min_padj(enr_df, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
tgfb_rank   <- pick_min_padj(enr_df, "HALLMARK_TGF_BETA_SIGNALING")

log("\n== Mitochondrial-energetic candidates (top 5 by OXPHOS padj) ==")
if (!is.null(oxphos_rank)) {
  for (i in seq_len(min(5, nrow(oxphos_rank)))) log("  ", oxphos_rank$module_color[i], " padj=", signif(oxphos_rank$p.adjust[i],3), " geneRatio=", oxphos_rank$GeneRatio[i])
}
log("\n== Anti-fibrotic-ECM candidates (top 5 by EMT padj) ==")
if (!is.null(emt_rank)) {
  for (i in seq_len(min(5, nrow(emt_rank)))) log("  ", emt_rank$module_color[i], " padj=", signif(emt_rank$p.adjust[i],3), " geneRatio=", emt_rank$GeneRatio[i])
}
log("\n== TGF_BETA top 5 (cross-check for anti-fibrotic) ==")
if (!is.null(tgfb_rank)) {
  for (i in seq_len(min(5, nrow(tgfb_rank)))) log("  ", tgfb_rank$module_color[i], " padj=", signif(tgfb_rank$p.adjust[i],3), " geneRatio=", tgfb_rank$GeneRatio[i])
}

# Combined criterion: OXPHOS primary + interaction-coefficient absolute value among top 3 OXPHOS candidates
top_n <- 3
oxphos_top <- if (!is.null(oxphos_rank)) head(oxphos_rank, top_n) else NULL
if (is.null(oxphos_top) || nrow(oxphos_top) == 0) stop("No module had any OXPHOS enrichment")
mito_pick <- oxphos_top %>%
  left_join(me_traits, by = "module_color") %>%
  arrange(desc(abs(beta_interaction))) %>%
  slice(1)
log("\nMITO-ENERGETIC anchor module: ", mito_pick$module_color,
    "  OXPHOS_padj=", signif(mito_pick$p.adjust,3),
    "  interaction_beta=", signif(mito_pick$beta_interaction,3),
    "  interaction_p=", signif(mito_pick$p_interaction,3))

# Anti-fibrotic ECM: rank by combined rank of EMT_padj + TGFB_padj across all enriched modules
combined_ecm <- enr_df %>%
  filter(ID %in% c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TGF_BETA_SIGNALING")) %>%
  group_by(module_color) %>%
  summarise(min_padj_pair = min(p.adjust),
            n_hits = n(),
            terms = paste(ID, collapse="|"),
            padjs = paste(signif(p.adjust,3), collapse="|"),
            .groups = "drop") %>%
  arrange(min_padj_pair, desc(n_hits))
log("\nAnti-fibrotic candidate ranking (combined EMT+TGFB):")
for (i in seq_len(min(5, nrow(combined_ecm)))) log("  ", combined_ecm$module_color[i], " min_padj=", signif(combined_ecm$min_padj_pair[i],3), " hits=", combined_ecm$n_hits[i])
ecm_pick <- combined_ecm %>% slice_head(n = top_n) %>%
  left_join(me_traits, by = "module_color") %>%
  arrange(desc(abs(beta_interaction))) %>%
  slice(1)
log("\nANTI-FIBROTIC-ECM anchor module: ", ecm_pick$module_color,
    "  combined_min_padj=", signif(ecm_pick$min_padj_pair,3),
    "  interaction_beta=", signif(ecm_pick$beta_interaction,3),
    "  interaction_p=", signif(ecm_pick$p_interaction,3))

# ---- Spearman rho gdf15.L vs ME ----
spearman_all <- bind_rows(lapply(colnames(MEs), function(mc) {
  v <- MEs[[mc]]
  ct <- suppressWarnings(cor.test(gdf15_vst, v, method = "spearman", exact = FALSE))
  lab <- sub("^ME","", mc)
  tibble(module = mc, module_label = lab, module_color = unname(label_to_color[lab]),
         rho_spearman = ct$estimate, p_spearman = ct$p.value)
})) %>% arrange(p_spearman)
write_tsv(spearman_all, file.path(WG_OUT, "brain_module_gdf15L_correlations.tsv"))
log("\nWrote brain_module_gdf15L_correlations.tsv")

mito_ME <- paste0("ME", color_to_label[mito_pick$module_color])
ecm_ME  <- paste0("ME", color_to_label[ecm_pick$module_color])
rho_mito <- spearman_all %>% filter(module == mito_ME)
rho_ecm  <- spearman_all %>% filter(module == ecm_ME)
log("\nFig 5.C anchor: ", mito_ME, " (color=", mito_pick$module_color, ")  rho=", signif(rho_mito$rho_spearman,3), "  p=", signif(rho_mito$p_spearman,3))
log("Fig 5.D anchor: ", ecm_ME,  " (color=", ecm_pick$module_color,  ")  rho=", signif(rho_ecm$rho_spearman,3),  "  p=", signif(rho_ecm$p_spearman,3))

# ---- Scatter panels ----
plot_df <- function(mod_color, label) {
  me_col <- paste0("ME", color_to_label[mod_color])
  if (!me_col %in% colnames(MEs)) stop("ME column not found: ", me_col, " for color ", mod_color)
  tibble(
    sample = brain_meta$sample,
    treatment = brain_meta$treatment,
    time = brain_meta$time,
    gdf15L_vst = gdf15_vst[brain_meta$sample],
    eigengene = MEs[[me_col]],
    module = label
  )
}
df_C <- plot_df(mito_pick$module_color, "mito-energetic")
df_D <- plot_df(ecm_pick$module_color,  "anti-fibrotic ECM")

panel_theme <- theme_classic(base_size = 10) +
  theme(axis.text  = element_text(size = 9, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        plot.title = element_text(size = 10, face = "plain", hjust = 0),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9))

mk_panel <- function(df, mod_color, rho, p, title) {
  ggplot(df, aes(x = gdf15L_vst, y = eigengene)) +
    geom_point(aes(color = treatment, shape = time), size = 3, stroke = 0.6) +
    scale_color_manual(values = c(ND = "#3B7DBF", CT = "#D9614C")) +
    scale_shape_manual(values = c(`24HPA` = 16, `7DPA` = 17)) +
    labs(title = title,
         subtitle = sprintf("ME%s — Spearman rho = %+.2f, p = %.3g (n=12)", mod_color, rho, p),
         x = expression(italic("gdf15.L")~"VST (brain)"),
         y = sprintf("ME%s eigengene", mod_color),
         color = "Treatment", shape = "Timepoint") +
    panel_theme
}

panel_C <- mk_panel(df_C, mito_pick$module_color, rho_mito$rho_spearman, rho_mito$p_spearman,
                    "Mitochondrial-energetic module vs brain gdf15.L")
panel_D <- mk_panel(df_D, ecm_pick$module_color,  rho_ecm$rho_spearman,  rho_ecm$p_spearman,
                    "Anti-fibrotic ECM module vs brain gdf15.L")

ggsave(file.path(FIG5_DIR, "C_panel.pdf"), panel_C, width = 90, height = 80, units = "mm")
ggsave(file.path(FIG5_DIR, "C_panel.png"), panel_C, width = 90, height = 80, units = "mm", dpi = 300)
ggsave(file.path(FIG5_DIR, "D_panel.pdf"), panel_D, width = 90, height = 80, units = "mm")
ggsave(file.path(FIG5_DIR, "D_panel.png"), panel_D, width = 90, height = 80, units = "mm", dpi = 300)

# panel_metadata.yaml per SO-12
mk_meta <- function(letter, mod_color, rho, p, role, top_hallmarks) {
  list(
    panel = paste0("Fig5.", letter),
    role = role,
    generated_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    script = "scripts/fig5_module_gdf15_rho.R",
    inputs = list(
      wgcna_rds = WGCNA_RDS,
      module_assignments_tsv = MOD_TSV,
      vsd_rds = VSD_RDS,
      sample_metadata_tsv = META_TSV
    ),
    anchor_module = paste0("ME", mod_color),
    module_color = mod_color,
    spearman_rho = signif(rho, 4),
    spearman_p = signif(p, 4),
    n_samples = 12,
    top_hallmark_hits = top_hallmarks,
    styling = list(
      width_mm = 90, height_mm = 80,
      treatment_colors = list(ND = "#3B7DBF", CT = "#D9614C"),
      time_shapes = list(`24HPA` = "circle", `7DPA` = "triangle"),
      regression_line = FALSE
    )
  )
}
mito_top <- enr_df %>% filter(module_color == mito_pick$module_color) %>% arrange(p.adjust) %>% head(5) %>% transmute(term = ID, padj = signif(p.adjust,3))
ecm_top  <- enr_df %>% filter(module_color == ecm_pick$module_color)  %>% arrange(p.adjust) %>% head(5) %>% transmute(term = ID, padj = signif(p.adjust,3))
write_yaml(mk_meta("C", mito_pick$module_color, rho_mito$rho_spearman, rho_mito$p_spearman, "mitochondrial-energetic module x brain gdf15.L", as.list(mito_top$term)),
           file.path(FIG5_DIR, "panel_metadata_C.yaml"))
write_yaml(mk_meta("D", ecm_pick$module_color,  rho_ecm$rho_spearman,  rho_ecm$p_spearman,  "anti-fibrotic ECM module x brain gdf15.L",   as.list(ecm_top$term)),
           file.path(FIG5_DIR, "panel_metadata_D.yaml"))

log("\nP1 COMPLETE. Outputs in ", FIG5_DIR, " and ", WG_OUT)
close(log_con)
