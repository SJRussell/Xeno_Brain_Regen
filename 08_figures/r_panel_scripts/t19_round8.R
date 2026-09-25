#!/usr/bin/env Rscript
# t19_round8.R — Round 8 brain plasticity arc deliverables.
# Eight jobs (analyses + panels + data files). Strict no-subtitle policy.
# Per-homeolog throughout. Overwrites prior versions in figures_current/v10_pulls/brain/.

suppressPackageStartupMessages({
    library(DESeq2); library(dplyr); library(readr); library(tibble);
    library(ggplot2); library(ggrepel); library(yaml); library(tidyr);
    library(clusterProfiler); library(org.Hs.eg.db); library(org.Xl.eg.db);
    library(msigdbr); library(AnnotationDbi); library(edgeR); library(limma)
})

PROJ    <- "/project/def-srussell/nirosha_xeno"
OUT     <- file.path(PROJ, "figures_current/v10_pulls/brain")
DEG     <- file.path(PROJ, "deg")
PATH    <- file.path(PROJ, "pathway")
WGCNA   <- file.path(PROJ, "wgcna")
MOFA    <- file.path(PROJ, "mofa")
COUNTS  <- file.path(PROJ, "counts", "star_gene_counts.tsv")
META_RDS <- file.path(DEG, "sample_metadata.rds")
META_TSV <- file.path(DEG, "sample_metadata.tsv")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("[", format(Sys.time(), "%FT%T%z"), "] t19 Round 8 start ->", OUT, "\n")

# ============================================================
# Universal styling (no subtitles ANYWHERE)
# ============================================================
panel_theme <- theme_classic(base_size = 11) +
    theme(axis.text.x  = element_text(size = 9.5, color = "black"),
          axis.text.y  = element_text(size = 9,   color = "black"),
          axis.title   = element_text(size = 10,  color = "black"),
          plot.title   = element_text(size = 11,  face = "plain", family = "sans", hjust = 0),
          plot.subtitle = element_blank(),
          legend.text  = element_text(size = 8.5),
          legend.title = element_text(size = 9),
          legend.position = "top",
          plot.margin  = margin(8, 12, 6, 8))

W_MM <- 170; H_MM <- 115

save_panel <- function(p, stem, w = W_MM, h = H_MM) {
    ggsave(file.path(OUT, paste0(stem, ".png")), p, width = w, height = h,
           units = "mm", dpi = 300, device = grDevices::png, type = "cairo")
    ggsave(file.path(OUT, paste0(stem, ".pdf")), p, width = w, height = h,
           units = "mm", device = cairo_pdf)
}
write_meta <- function(stem, info) {
    info$generated_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    info$script <- "scripts/t19_round8.R"
    write_yaml(info, file.path(OUT, paste0(stem, "_meta.yaml")))
}

sig_star <- function(p) ifelse(is.na(p), "ns",
                        ifelse(p < 0.001, "***",
                        ifelse(p < 0.01,  "**",
                        ifelse(p < 0.05,  "*", "ns"))))

# ============================================================
# Load shared inputs
# ============================================================
meta_full <- readRDS(META_RDS)
brain_pri <- readRDS(file.path(DEG, "deseq2_BRAIN_primary.rds"))
blast_pri <- readRDS(file.path(DEG, "deseq2_BLAST_primary.rds"))

# ============================================================
# JOB 1 — LRT per-homeolog (overwrite brain_lrt_trtxtime_deseq2)
# ============================================================
cat("\n[Job 1] LRT per-homeolog brain trt x time\n")
mat_full <- as.matrix(read.table(COUNTS, header = TRUE, sep = "\t",
                                  row.names = 1, check.names = FALSE,
                                  stringsAsFactors = FALSE))
storage.mode(mat_full) <- "integer"
brain_late_idx <- meta_full$tissue == "BRAIN" &
                  meta_full$time %in% c("24HPA", "7DPA")
mat_b <- mat_full[, brain_late_idx]
meta_b <- droplevels(meta_full[brain_late_idx, ])
# CPM>1 in >=3 samples filter (matches t05c)
keep_b <- rowSums(cpm(mat_b) > 1) >= 3
mat_b <- mat_b[keep_b, ]
cat("  Per-homeolog features after filter:", nrow(mat_b), ", samples:", ncol(mat_b), "\n")

dds_lrt <- DESeqDataSetFromMatrix(countData = mat_b, colData = meta_b,
                                   design = ~ time + treatment + treatment:time)
dds_lrt <- DESeq(dds_lrt, test = "LRT",
                 reduced = ~ time + treatment)
res_lrt <- results(dds_lrt)
df_lrt <- as.data.frame(res_lrt)
df_lrt$gene_id <- rownames(df_lrt)
n_sig_lrt <- sum(df_lrt$padj < 0.05, na.rm = TRUE)
cat("  LRT per-homeolog sig at q<0.05:", n_sig_lrt, "\n")

# Save LRT object replacing the homeolog-sum version's data slot for downstream
saveRDS(list(dds = dds_lrt, res = res_lrt,
             label = "BRAIN_late_7v24_interaction_LRT_per_homeolog",
             tissue = "BRAIN", timepoints = c("24HPA", "7DPA")),
        file.path(DEG, "deseq2_lrt_BRAIN_late_7v24_per_homeolog.rds"))

# ============================================================
# Volcano helper (Round 7 + 8 conventions, NO subtitle)
# ============================================================
make_volcano <- function(df, stem, title,
                         up_label = "Up in CT", down_label = "Down in CT",
                         top_n_labels = 10, clip_x = c(-5, 5),
                         color_up = "#B2182B", color_down = "#2166AC",
                         x_axis_label = "log\u2082 FC (CT vs ND)",
                         legend_ns_label = "n.s.") {
    df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
    df$neg_log_padj <- -log10(df$padj)
    df$direction <- factor(
        ifelse(df$padj < 0.05 & df$log2FoldChange > 0, "up",
        ifelse(df$padj < 0.05 & df$log2FoldChange < 0, "down", "ns")),
        levels = c("up", "down", "ns"))

    top <- df[df$direction != "ns", ]
    if ("stat" %in% colnames(top) && any(!is.na(top$stat))) {
        top <- top[order(-abs(top$stat))[seq_len(min(top_n_labels, nrow(top)))], ]
    } else {
        top <- top[order(top$padj)[seq_len(min(top_n_labels, nrow(top)))], ]
    }
    if (grepl("7dpt", stem) && !grepl("7v24", stem) &&
        "gdf15.L" %in% df$gene_id && !"gdf15.L" %in% top$gene_id) {
        top <- rbind(top, df[df$gene_id == "gdf15.L", ])
    }

    df$display_lfc  <- pmax(pmin(df$log2FoldChange, clip_x[2]), clip_x[1])
    df$off_scale    <- df$log2FoldChange < clip_x[1] | df$log2FoldChange > clip_x[2]
    top$display_lfc <- pmax(pmin(top$log2FoldChange, clip_x[2]), clip_x[1])
    top$off_scale   <- top$log2FoldChange < clip_x[1] | top$log2FoldChange > clip_x[2]

    dir_palette <- c("up" = color_up, "down" = color_down, "ns" = "grey75")
    dir_labels  <- c("up" = up_label, "down" = down_label, "ns" = legend_ns_label)

    df_ordered <- df[order(df$direction != "ns"), ]
    p <- ggplot(df_ordered, aes(x = display_lfc, y = neg_log_padj)) +
        geom_point(aes(color = direction, shape = off_scale),
                   alpha = 0.75, size = 1.5, stroke = 0.25) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                   color = "grey50", linewidth = 0.35) +
        geom_text_repel(data = top, aes(x = display_lfc, y = neg_log_padj,
                                         label = gene_id),
                        size = 2.7, color = "black", fontface = "italic",
                        box.padding = 0.35, point.padding = 0.2,
                        segment.size = 0.25, segment.color = "grey55",
                        max.overlaps = Inf, min.segment.length = 0) +
        scale_color_manual(values = dir_palette, labels = dir_labels,
                            name = NULL, drop = FALSE,
                            guide = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
        scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 18), guide = "none") +
        coord_cartesian(xlim = clip_x) +
        labs(title = title, x = x_axis_label,
             y = expression(-log[10]~"(padj)")) +
        panel_theme
    save_panel(p, stem)

    list(n_sig = sum(df$direction != "ns"),
         n_up  = sum(df$direction == "up"),
         n_down = sum(df$direction == "down"),
         n_off = sum(df$off_scale & df$direction != "ns"),
         top_hits = lapply(seq_len(nrow(top)), function(i) list(
             gene_id = top$gene_id[i],
             log2FC = signif(top$log2FoldChange[i], 3),
             padj   = signif(top$padj[i], 3))))
}

cat("  Rendering brain_lrt_trtxtime_deseq2 (per-homeolog)...\n")
info_lrt <- make_volcano(df_lrt,
    "brain_lrt_trtxtime_deseq2",
    "BRAIN treatment \u00d7 time interaction",
    up_label = "CT trajectory \u2191 vs ND",
    down_label = "CT trajectory \u2193 vs ND",
    color_up = "#1B7837", color_down = "#762A83",
    x_axis_label = "log\u2082 FC interaction (CT_7DPA \u2212 CT_24HPA) \u2212 (ND_7DPA \u2212 ND_24HPA)")
write_meta("brain_lrt_trtxtime_deseq2", list(
    panel = "Fig 2C", title = "BRAIN treatment \u00d7 time interaction",
    analysis_engine = "DESeq2 LRT, per-homeolog (full ~ time + treatment + treatment:time vs reduced ~ time + treatment)",
    framing = "per-homeolog (Round 8 override of homeolog-sum)",
    fdr_threshold = 0.05, n_sig = info_lrt$n_sig,
    n_up = info_lrt$n_up, n_down = info_lrt$n_down, n_off_scale = info_lrt$n_off,
    input_table = "deg/deseq2_lrt_BRAIN_late_7v24_per_homeolog.rds",
    top_hits = info_lrt$top_hits))
cat("  -> brain_lrt_trtxtime_deseq2: ", info_lrt$n_sig, " sig (",
    info_lrt$n_up, " up, ", info_lrt$n_down, " down)\n", sep = "")

# ============================================================
# JOB 2 — Brain ECM GO BP at 24 HPT
# ============================================================
cat("\n[Job 2] Brain ECM GO BP at 24 HPT\n")
enrich <- readRDS(file.path(PATH, "all_enrichment_results.rds"))
x_ecm <- enrich[["deseq2_BRAIN_primary_CT_vs_ND_24HPA_goBP"]]
df_ecm <- if (isS4(x_ecm) && "result" %in% slotNames(x_ecm)) as.data.frame(x_ecm@result) else as.data.frame(x_ecm)
ecm_pat <- "extracellular matrix|collagen|basement membrane|cell-substrate adhesion|extracellular structure|connective tissue|matrix remodel|fibronectin|laminin|proteoglycan"
ecm <- df_ecm[grepl(ecm_pat, df_ecm$Description, ignore.case = TRUE), ]
ecm <- ecm[order(ecm$p.adjust), ][seq_len(min(10, nrow(ecm))), ]
ecm$term_lbl <- factor(ecm$Description, levels = rev(ecm$Description))
cat("  ECM-keyword GO BP terms shown:", nrow(ecm), "\n")
x_range <- max(abs(ecm$NES))
x_lim <- c(-x_range * 1.20, x_range * 1.20)
p_ecm <- ggplot(ecm, aes(x = NES, y = term_lbl, fill = NES)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_col(color = "grey25", linewidth = 0.2, width = 0.7) +
    geom_text(aes(x = 0,
                  label = sprintf("q = %.2g", p.adjust),
                  hjust = ifelse(NES > 0, 1.05, -0.05)),
              size = 2.7, color = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "NES") +
    coord_cartesian(xlim = x_lim) +
    labs(title = "BRAIN ECM GO BP NES at 24 HPT",
         x = "NES (CT vs ND, 24 HPT)", y = NULL) +
    panel_theme + theme(legend.position = "right")
save_panel(p_ecm, "brain_ecm_go_24hpt_deseq2", h = max(85, 10 * nrow(ecm) + 40))
write_meta("brain_ecm_go_24hpt_deseq2", list(
    panel = "Fig 3B (paired with synaptic 7DPT)",
    title = "BRAIN ECM GO BP NES at 24 HPT",
    analysis_engine = "clusterProfiler GSEA on GO BP; DESeq2 Wald ranks; ECM-keyword filter",
    keyword_filter = ecm_pat, fdr_threshold = 0.05, n_sig = nrow(ecm),
    input_table = "pathway/all_enrichment_results.rds [deseq2_BRAIN_primary_CT_vs_ND_24HPA_goBP]",
    terms = lapply(seq_len(nrow(ecm)), function(i) list(
        ID = ecm$ID[i], Description = ecm$Description[i],
        NES = signif(ecm$NES[i], 3),
        padj = signif(ecm$p.adjust[i], 3)))))
cat("  -> brain_ecm_go_24hpt_deseq2: ", nrow(ecm), " ECM terms\n", sep = "")

# ============================================================
# JOB 3 — Salmon4 module Hallmark over-rep
# ============================================================
cat("\n[Job 3] Salmon4 Hallmark over-rep\n")
hm_all <- read_tsv(file.path(WGCNA, "brain_module_hallmark_enrichment.tsv"),
                    show_col_types = FALSE)
salmon4 <- hm_all[hm_all$module_color == "salmon4", ]
salmon4 <- salmon4[order(salmon4$p.adjust), ][seq_len(min(5, nrow(salmon4))), ]
salmon4$term_short <- sub("^HALLMARK_", "", salmon4$ID)
salmon4$term_lbl   <- factor(salmon4$term_short, levels = rev(salmon4$term_short))
salmon4$neg_log_q  <- -log10(salmon4$p.adjust)
salmon4_n_sig <- sum(salmon4$p.adjust < 0.05)
cat("  salmon4 Hallmarks at q<0.05:", salmon4_n_sig, "\n")
print(salmon4[, c("term_short", "Count", "p.adjust")])

p_salmon <- ggplot(salmon4, aes(x = Count, y = term_lbl, fill = neg_log_q)) +
    geom_col(color = "grey25", linewidth = 0.2, width = 0.65) +
    geom_text(aes(label = sprintf("q=%.1e", p.adjust)),
              hjust = -0.1, size = 2.7, color = "black") +
    scale_fill_gradient(low = "#FEE0B6", high = "#7F3B08",
                        name = expression(-log[10]~q)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.32))) +
    labs(title = "Brain anti-fibrotic candidate module (salmon4, n = 47): Hallmark over-rep",
         x = "Overlap count", y = NULL) +
    panel_theme + theme(legend.position = "right")
save_panel(p_salmon, "brain_salmon4_hallmark_overrep", h = 80)
write_meta("brain_salmon4_hallmark_overrep", list(
    panel = "Fig 4D paired",
    title = "Brain anti-fibrotic candidate module (salmon4, n = 47): Hallmark over-rep",
    analysis_engine = "Hypergeometric over-rep (precomputed in wgcna/brain_module_hallmark_enrichment.tsv)",
    fdr_threshold = 0.05, n_sig_at_q05 = salmon4_n_sig,
    module_n_genes = 47, module = "salmon4",
    input_table = "wgcna/brain_module_hallmark_enrichment.tsv [filter module_color==salmon4]",
    terms = lapply(seq_len(nrow(salmon4)), function(i) list(
        ID = salmon4$ID[i], Count = salmon4$Count[i],
        padj = signif(salmon4$p.adjust[i], 4),
        fold_enrich = signif(salmon4$FoldEnrichment[i], 3)))))
cat("  -> brain_salmon4_hallmark_overrep: ", nrow(salmon4),
    " Hallmarks shown, ", salmon4_n_sig, " sig at q<0.05\n", sep = "")

# Salmon4 gene list TSV (for speaker notes)
salmon4_genes <- read_tsv(file.path(WGCNA, "module_assignments_BRAIN.tsv"),
                          show_col_types = FALSE)
salmon4_genes <- salmon4_genes[salmon4_genes$module_color == "salmon4", ]
write_tsv(salmon4_genes, file.path(OUT, "brain_salmon4_gene_list.tsv"))
cat("  -> brain_salmon4_gene_list.tsv: ", nrow(salmon4_genes), " genes\n", sep = "")

# ============================================================
# JOBS 4-5 — Hallmark GSEA on CT_7v24 and ND_7v24
# ============================================================
cat("\n[Jobs 4-5] Hallmark GSEA on CT_7v24 and ND_7v24\n")

# Build Xla -> Hs ortholog map via case-folded symbols (same as t06)
strip_xla <- function(s) sub("\\.[LS]$", "", s)
hs_sym <- keys(org.Hs.eg.db, keytype = "SYMBOL")
hs_lookup <- AnnotationDbi::select(org.Hs.eg.db, keys = hs_sym, keytype = "SYMBOL",
                                    columns = c("SYMBOL", "ENTREZID"))
hs_lookup <- setNames(hs_lookup$ENTREZID, toupper(hs_lookup$SYMBOL))

hallmark_human <- msigdbr(species = "Homo sapiens", collection = "H")
entrez_col <- intersect(c("ncbi_gene", "human_entrez_gene", "entrez_gene"),
                        colnames(hallmark_human))[1]
hallmark_h <- unique(as.data.frame(hallmark_human[, c("gs_name", entrez_col)]))
colnames(hallmark_h) <- c("gs_name", "entrez")
hallmark_h <- hallmark_h[!is.na(hallmark_h$entrez), ]

build_rank_hs <- function(res_df) {
    rk <- sign(res_df$log2FoldChange) * (-log10(res_df$pvalue))
    names(rk) <- toupper(strip_xla(res_df$gene_id))
    rk <- rk[!is.na(rk) & is.finite(rk)]
    hs_entrez <- hs_lookup[names(rk)]
    keep <- !is.na(hs_entrez) & !duplicated(hs_entrez)
    rk_h <- rk[keep]; names(rk_h) <- hs_entrez[keep]
    sort(rk_h, decreasing = TRUE)
}

run_hallmark_gsea <- function(res_df, label) {
    rk <- build_rank_hs(res_df)
    set.seed(42)
    res <- tryCatch({
        GSEA(rk, TERM2GENE = hallmark_h[, c("gs_name", "entrez")],
             pvalueCutoff = 1, minGSSize = 15, maxGSSize = 500, eps = 0)
    }, error = function(e) {cat("    GSEA error (", label, "):", conditionMessage(e), "\n"); NULL})
    if (is.null(res) || nrow(res@result) == 0) return(NULL)
    res@result
}

render_hallmark_dot <- function(df, stem, title, time_label, n_show_cap = NULL) {
    sig <- df[df$p.adjust < 0.05, ]
    n_sig_total <- nrow(sig)
    if (n_sig_total == 0) {
        cat("  [WARN] no sig at q<0.05; showing top 10 by p for", stem, "\n")
        sig <- df[order(df$p.adjust), ][seq_len(min(10, nrow(df))), ]
    } else {
        sig <- sig[order(-sig$NES), ]
        if (!is.null(n_show_cap)) sig <- sig[seq_len(min(n_show_cap, nrow(sig))), ]
    }
    sig$term_short <- sub("^HALLMARK_", "", sig$ID)
    sig$term_lbl   <- factor(sig$term_short, levels = rev(sig$term_short))
    sig$neg_log_q  <- -log10(sig$p.adjust)
    p <- ggplot(sig, aes(x = NES, y = term_lbl, color = NES, size = neg_log_q)) +
        geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70") +
        geom_point(stroke = 0.4, alpha = 0.9) +
        scale_color_gradient2(low = "#2166AC", mid = "#F0F0F0", high = "#B2182B",
                              midpoint = 0, name = "NES") +
        scale_size_continuous(range = c(2.5, 7), name = expression(-log[10]~q)) +
        labs(title = title, x = sprintf("NES (%s)", time_label), y = NULL) +
        panel_theme + theme(legend.position = "right")
    save_panel(p, stem, h = max(80, 8 * nrow(sig) + 40))
    list(n_sig_total = n_sig_total,
         terms = lapply(seq_len(nrow(sig)), function(i) list(
             ID = sig$ID[i], NES = signif(sig$NES[i], 3),
             padj = signif(sig$p.adjust[i], 3))))
}

cat("  CT_7v24 Hallmark GSEA...\n")
gsea_ct <- run_hallmark_gsea(brain_pri$results[["CT_7DPA_vs_24HPA"]],
                              "BRAIN CT_7v24")
if (!is.null(gsea_ct)) {
    info_ct <- render_hallmark_dot(gsea_ct, "brain_hallmark_CT_7v24_deseq2",
        "BRAIN Hallmark NES: CT 7 DPT vs 24 HPT",
        "7 DPT vs 24 HPT, CT animals")
    write_meta("brain_hallmark_CT_7v24_deseq2", list(
        panel = "Fig 2H paired",
        title = "BRAIN Hallmark NES: CT 7 DPT vs 24 HPT",
        analysis_engine = "clusterProfiler GSEA Hallmark on signed -log10(p)*sign(LFC); DESeq2 Wald ranks",
        contrast = "CT_7DPA vs CT_24HPA, per-homeolog brain",
        fdr_threshold = 0.05, n_sig = info_ct$n_sig_total,
        terms = info_ct$terms))
    cat("  -> brain_hallmark_CT_7v24_deseq2: ", info_ct$n_sig_total, " sig\n")
}

cat("  ND_7v24 Hallmark GSEA...\n")
gsea_nd <- run_hallmark_gsea(brain_pri$results[["ND_7DPA_vs_24HPA"]],
                              "BRAIN ND_7v24")
if (!is.null(gsea_nd)) {
    info_nd <- render_hallmark_dot(gsea_nd, "brain_hallmark_ND_7v24_deseq2",
        "BRAIN Hallmark NES: ND 7 DPT vs 24 HPT",
        "7 DPT vs 24 HPT, ND animals")
    write_meta("brain_hallmark_ND_7v24_deseq2", list(
        panel = "Fig 2G paired",
        title = "BRAIN Hallmark NES: ND 7 DPT vs 24 HPT",
        analysis_engine = "clusterProfiler GSEA Hallmark on signed -log10(p)*sign(LFC); DESeq2 Wald ranks",
        contrast = "ND_7DPA vs ND_24HPA, per-homeolog brain",
        fdr_threshold = 0.05, n_sig = info_nd$n_sig_total,
        terms = info_nd$terms))
    cat("  -> brain_hallmark_ND_7v24_deseq2: ", info_nd$n_sig_total, " sig\n")
}

# ============================================================
# JOB 6 — MOFA Factor 3 BRAIN-view loadings
# ============================================================
cat("\n[Job 6] MOFA Factor 3 BRAIN-view loadings\n")
ml <- read_tsv(file.path(MOFA, "factor_loadings_top.tsv"), show_col_types = FALSE)
f3 <- ml[ml$factor == "Factor3" & ml$view == "BRAIN", ]
f3$abs_w <- abs(f3$weight)
f3 <- f3[order(-f3$abs_w), ][seq_len(min(10, nrow(f3))), ]
if (nrow(f3) > 0) {
    f3$gene_clean <- sub("_BRAIN$", "", f3$gene)
    f3$gene_lbl   <- factor(f3$gene_clean, levels = rev(f3$gene_clean))
    p_f3 <- ggplot(f3, aes(x = weight, y = gene_lbl, fill = weight)) +
        geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
        geom_col(color = "grey25", linewidth = 0.2, width = 0.7) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                             midpoint = 0, name = "weight") +
        labs(title = "MOFA+ Factor 3: top 10 BRAIN-view loadings",
             x = "Factor 3 weight (BRAIN view)", y = NULL) +
        panel_theme +
        theme(axis.text.y = element_text(face = "italic"),
              legend.position = "right")
    save_panel(p_f3, "mofa_factor3_brain_loadings", h = 100)
    write_meta("mofa_factor3_brain_loadings", list(
        panel = "Fig 6B (paired with mofa_factor3_blast_loadings)",
        title = "MOFA+ Factor 3 top 10 BRAIN-view loadings",
        analysis_engine = "MOFA+ factor model (precomputed)",
        input_table = file.path(MOFA, "factor_loadings_top.tsv"),
        top_hits = lapply(seq_len(nrow(f3)), function(i) list(
            gene_id = f3$gene_clean[i],
            weight = signif(f3$weight[i], 4)))))
    cat("  -> mofa_factor3_brain_loadings: ", nrow(f3), " loadings\n", sep = "")
}

# ============================================================
# JOB 7 — GDF15 significance grid
# ============================================================
cat("\n[Job 7] gdf15 significance grid\n")
extract_sig <- function(res_obj, contrast_keys, tissue, contrast_label_map) {
    rows <- list()
    for (gene in c("gdf15.L", "gdf15.S")) {
        for (cn in contrast_keys) {
            df <- res_obj$results[[cn]]
            if (is.null(df) || !(gene %in% df$gene_id)) {
                rows[[length(rows)+1]] <- tibble(
                    tissue = tissue, homeolog = gene,
                    contrast = contrast_label_map[[cn]],
                    logFC = NA_real_, padj = NA_real_, sig = "filtered_or_absent")
                next
            }
            i <- which(df$gene_id == gene)[1]
            p <- df$padj[i]; lfc <- df$log2FoldChange[i]
            rows[[length(rows)+1]] <- tibble(
                tissue = tissue, homeolog = gene,
                contrast = contrast_label_map[[cn]],
                logFC = signif(lfc, 4),
                padj  = signif(p, 4),
                sig   = sig_star(p))
        }
    }
    bind_rows(rows)
}

brain_map <- list(CT_vs_ND_24HPA="CT_vs_ND_24HPT", CT_vs_ND_7DPA="CT_vs_ND_7DPT",
                  CT_7DPA_vs_24HPA="CT_7v24", ND_7DPA_vs_24HPA="ND_7v24")
blast_map <- list(CT_vs_ND_11HPA="CT_vs_ND_11HPT", CT_vs_ND_24HPA="CT_vs_ND_24HPT",
                  CT_vs_ND_7DPA="CT_vs_ND_168HPT")

gdf_grid <- bind_rows(
    extract_sig(brain_pri, names(brain_map), "brain", brain_map),
    extract_sig(blast_pri, names(blast_map), "blastema", blast_map))
write_tsv(gdf_grid, file.path(OUT, "gdf15_significance_grid.tsv"))
cat("  -> gdf15_significance_grid.tsv (", nrow(gdf_grid), " rows)\n", sep = "")

# ============================================================
# JOB 8 — Trajectory gene significance grid
# ============================================================
cat("\n[Job 8] trajectory gene significance grid\n")
panel3b_genes <- c("slc17a7.L","gria1.L","syn2.S","grin2b.L","syp.S","synm.L","synm.S")
panel5d_genes <- c("tgfb2.L","tgfb2.S","sox9.L","mmp13.S","col17a1.L","col1a1.L","col3a1.S","fn1.S")
panel5c_genes <- c("gdf15.L","gdf15.S")
panel4a_genes <- c("gdf15.L")

extract_genes <- function(res_obj, genes, contrast_keys, tissue,
                          contrast_label_map, panel_tag) {
    rows <- list()
    for (gene in genes) {
        for (cn in contrast_keys) {
            df <- res_obj$results[[cn]]
            if (is.null(df) || !(gene %in% df$gene_id)) {
                rows[[length(rows)+1]] <- tibble(
                    panel = panel_tag, tissue = tissue, gene = gene,
                    contrast = contrast_label_map[[cn]],
                    logFC = NA_real_, padj = NA_real_, sig = "filtered_or_absent")
                next
            }
            i <- which(df$gene_id == gene)[1]
            p <- df$padj[i]; lfc <- df$log2FoldChange[i]
            rows[[length(rows)+1]] <- tibble(
                panel = panel_tag, tissue = tissue, gene = gene,
                contrast = contrast_label_map[[cn]],
                logFC = signif(lfc, 4),
                padj  = signif(p, 4),
                sig   = sig_star(p))
        }
    }
    bind_rows(rows)
}

brain_ctnd <- c("CT_vs_ND_24HPA","CT_vs_ND_7DPA")
brain_ctnd_lbl <- list(CT_vs_ND_24HPA="CT_vs_ND_24HPT", CT_vs_ND_7DPA="CT_vs_ND_7DPT")
blast_ctnd <- c("CT_vs_ND_11HPA","CT_vs_ND_24HPA","CT_vs_ND_7DPA")
blast_ctnd_lbl <- list(CT_vs_ND_11HPA="CT_vs_ND_11HPT",
                       CT_vs_ND_24HPA="CT_vs_ND_24HPT",
                       CT_vs_ND_7DPA="CT_vs_ND_168HPT")

trj_grid <- bind_rows(
    extract_genes(brain_pri, panel3b_genes, brain_ctnd, "brain", brain_ctnd_lbl, "3B"),
    extract_genes(brain_pri, panel4a_genes, brain_ctnd, "brain", brain_ctnd_lbl, "4A"),
    extract_genes(blast_pri, panel5c_genes, blast_ctnd, "blastema", blast_ctnd_lbl, "5C"),
    extract_genes(blast_pri, panel5d_genes, blast_ctnd, "blastema", blast_ctnd_lbl, "5D"))
write_tsv(trj_grid, file.path(OUT, "trajectory_gene_significance_grid.tsv"))
cat("  -> trajectory_gene_significance_grid.tsv (", nrow(trj_grid), " rows)\n", sep = "")

# ============================================================
# JOB 9 — Top-10 hits per brain DE contrast (text dump)
# ============================================================
cat("\n[Job 9] brain_de_top10_hits.txt\n")
emit_top10 <- function(res_df, label) {
    df <- res_df[!is.na(res_df$padj), ]
    df <- df[order(df$padj), ][seq_len(min(10, nrow(df))), ]
    lines <- sprintf("%s (%+.2f, %.2g)", df$gene_id, df$log2FoldChange, df$padj)
    c(label, lines, "")
}

top10_lines <- c(
    emit_top10(brain_pri$results[["CT_vs_ND_24HPA"]], "BRAIN_CT_vs_ND_24HPT"),
    emit_top10(brain_pri$results[["CT_vs_ND_7DPA"]],  "BRAIN_CT_vs_ND_7DPT"),
    emit_top10(brain_pri$results[["ND_7DPA_vs_24HPA"]], "BRAIN_ND_7v24"),
    emit_top10(brain_pri$results[["CT_7DPA_vs_24HPA"]], "BRAIN_CT_7v24"),
    emit_top10({tmp <- df_lrt; tmp[order(tmp$padj), ]}, "BRAIN_LRT_trt_x_time"))
writeLines(top10_lines, file.path(OUT, "brain_de_top10_hits.txt"))
cat("  -> brain_de_top10_hits.txt (", length(top10_lines), " lines)\n", sep = "")

# ============================================================
# JOB 10 — Trajectory panels with sig stars
# ============================================================
cat("\n[Job 10] Trajectory panels with sig stars\n")
voom_BRAIN <- readRDS(file.path(DEG, "voom_BRAIN_primary.rds"))
voom_BLAST <- readRDS(file.path(DEG, "voom_BLAST_primary.rds"))
expr_brain <- voom_BRAIN$v$E
expr_blast <- voom_BLAST$v$E
meta_t <- read_tsv(META_TSV, show_col_types = FALSE)
meta_brain <- meta_t %>% filter(tissue == "BRAIN") %>%
    mutate(treatment = factor(treatment, levels = c("ND", "CT")),
           time = factor(time, levels = c("24HPA", "7DPA")),
           time_lbl = factor(ifelse(time == "24HPA", "24 HPT", "7 DPT"),
                             levels = c("24 HPT", "7 DPT")))
meta_blast <- meta_t %>% filter(tissue == "BLAST") %>%
    mutate(treatment = factor(treatment, levels = c("ND", "CT")),
           time = factor(time, levels = c("11HPA","24HPA","7DPA")),
           time_lbl = factor(case_when(time == "11HPA" ~ "11 HPT",
                                        time == "24HPA" ~ "24 HPT",
                                        time == "7DPA"  ~ "7 DPT"),
                              levels = c("11 HPT","24 HPT","7 DPT")))

TRT_COL <- c(ND = "#3B7DBF", CT = "#D9614C")

make_trajectory <- function(genes, expr_mat, meta_sub, panel_tag, stem, title) {
    genes_present <- genes[genes %in% rownames(expr_mat)]
    missing <- setdiff(genes, genes_present)
    if (length(missing) > 0) cat("  [", panel_tag, "] missing: ",
                                  paste(missing, collapse = ", "), "\n", sep = "")
    if (length(genes_present) == 0) {
        cat("  [", panel_tag, "] no genes available; skipping\n"); return(NULL)
    }
    samp_ids <- meta_sub$sample
    expr_g <- expr_mat[genes_present, samp_ids, drop = FALSE]
    dat <- do.call(rbind, lapply(genes_present, function(g)
        tibble(sample = samp_ids, gene = g, expression = as.numeric(expr_g[g, ]))))
    dat <- dat %>% left_join(
        meta_sub %>% dplyr::select(sample, treatment, time_lbl), by = "sample") %>%
        mutate(gene = factor(gene, levels = genes_present))
    sumtbl <- dat %>% group_by(gene, treatment, time_lbl) %>%
        summarise(mean = mean(expression), sem = sd(expression)/sqrt(n()),
                  .groups = "drop")

    # Sig-star annotations per gene per timepoint (from CT-vs-ND-at-time contrast)
    sig_anno <- trj_grid %>% filter(panel == panel_tag,
                                    gene %in% genes_present) %>%
        filter(sig %in% c("*","**","***"))
    cat("  [", panel_tag, "] sig stars: ", nrow(sig_anno), "\n", sep = "")

    sig_anno <- sig_anno %>%
        mutate(time_lbl = case_when(
            grepl("11HPT", contrast) ~ "11 HPT",
            grepl("24HPT", contrast) ~ "24 HPT",
            grepl("7DPT|168HPT", contrast) ~ "7 DPT",
            TRUE ~ NA_character_)) %>%
        filter(!is.na(time_lbl)) %>%
        left_join(dat %>% group_by(gene, time_lbl) %>%
                      summarise(y_top = max(expression, na.rm = TRUE),
                                .groups = "drop"),
                  by = c("gene","time_lbl")) %>%
        mutate(gene = factor(gene, levels = genes_present))

    p <- ggplot(dat, aes(x = time_lbl, y = expression,
                          color = treatment, group = treatment)) +
        geom_line(data = sumtbl, aes(y = mean, group = treatment),
                  linewidth = 0.8, alpha = 0.85) +
        geom_errorbar(data = sumtbl,
                      aes(x = time_lbl, ymin = mean - sem, ymax = mean + sem,
                          color = treatment, group = treatment),
                      width = 0.12, linewidth = 0.5, alpha = 0.85,
                      inherit.aes = FALSE) +
        geom_jitter(width = 0.08, height = 0, size = 2.0, alpha = 0.85, stroke = 0) +
        scale_color_manual(values = TRT_COL, name = "Treatment") +
        facet_wrap(~ gene, ncol = ifelse(length(genes_present) <= 4, length(genes_present),
                                          ceiling(sqrt(length(genes_present)))),
                   scales = "free_y") +
        labs(title = title, x = "Timepoint", y = "Normalized expression") +
        panel_theme +
        theme(strip.background = element_blank(),
              strip.text = element_text(size = 10, face = "italic"),
              panel.spacing = unit(10, "pt"),
              legend.position = "right")

    if (nrow(sig_anno) > 0) {
        p <- p + geom_text(data = sig_anno,
                            aes(x = time_lbl, y = y_top, label = sig),
                            inherit.aes = FALSE,
                            vjust = -0.3, size = 4.0, color = "black",
                            fontface = "bold")
    }
    h <- max(95, 55 * length(genes_present) %/% 4 + 60)
    save_panel(p, stem, h = h)

    list(genes = genes_present, missing = missing,
         n_stars = nrow(sig_anno))
}

cat("  Panel 3B brain synaptic gene trajectories...\n")
info_3b <- make_trajectory(panel3b_genes, expr_brain, meta_brain, "3B",
    "brain_synaptic_gene_trajectories_v10",
    "BRAIN synaptic genes by treatment \u00d7 time")
write_meta("brain_synaptic_gene_trajectories_v10", list(
    panel = "Fig 3B (v10 — sig stars added)",
    title = "BRAIN synaptic genes by treatment x time",
    analysis_engine = "voom-normalized expression; CT-vs-ND DESeq2 Wald sig stars (ns/*/**/***)",
    genes_requested = panel3b_genes, genes_rendered = info_3b$genes,
    genes_missing = info_3b$missing, n_sig_stars = info_3b$n_stars,
    input_table = "deg/voom_BRAIN_primary.rds $v$E (expression); deg/deseq2_BRAIN_primary.rds (sig)"))
cat("  -> brain_synaptic_gene_trajectories_v10: ", length(info_3b$genes),
    " genes, ", info_3b$n_stars, " stars\n", sep = "")

cat("  Panel 4A brain gdf15.L trajectory...\n")
info_4a <- make_trajectory("gdf15.L", expr_brain, meta_brain, "4A",
    "brain_gdf15L_trajectory_v10",
    "BRAIN gdf15.L by treatment \u00d7 time")
write_meta("brain_gdf15L_trajectory_v10", list(
    panel = "Fig 4A (v10 — sig stars added)",
    title = "BRAIN gdf15.L by treatment x time",
    analysis_engine = "voom-normalized expression; CT-vs-ND DESeq2 Wald sig stars (ns/*/**/***)",
    genes_rendered = info_4a$genes,
    n_sig_stars = info_4a$n_stars,
    input_table = "deg/voom_BRAIN_primary.rds $v$E (expression); deg/deseq2_BRAIN_primary.rds (sig)"))
cat("  -> brain_gdf15L_trajectory_v10: ", info_4a$n_stars, " stars\n", sep = "")

cat("  Panel 5C blastema gdf15 trajectories...\n")
info_5c <- make_trajectory(panel5c_genes, expr_blast, meta_blast, "5C",
    "blastema_gdf15_trajectories_v10",
    "BLASTEMA gdf15 (.L/.S) by treatment \u00d7 time")
write_meta("blastema_gdf15_trajectories_v10", list(
    panel = "Fig 5C (v10 — sig stars added)",
    title = "BLASTEMA gdf15 (.L/.S) by treatment x time",
    analysis_engine = "voom-normalized expression; CT-vs-ND DESeq2 Wald sig stars (ns/*/**/***)",
    genes_rendered = info_5c$genes,
    genes_missing = info_5c$missing, n_sig_stars = info_5c$n_stars,
    input_table = "deg/voom_BLAST_primary.rds $v$E (expression); deg/deseq2_BLAST_primary.rds (sig)"))
cat("  -> blastema_gdf15_trajectories_v10: ", length(info_5c$genes),
    " genes, ", info_5c$n_stars, " stars\n", sep = "")

cat("  Panel 5D blastema anti-fibrotic trajectories...\n")
info_5d <- make_trajectory(panel5d_genes, expr_blast, meta_blast, "5D",
    "blastema_anti_fibrotic_trajectories_v10",
    "BLASTEMA anti-fibrotic candidates by treatment \u00d7 time")
write_meta("blastema_anti_fibrotic_trajectories_v10", list(
    panel = "Fig 5D (v10 — sig stars added)",
    title = "BLASTEMA anti-fibrotic candidates by treatment x time",
    analysis_engine = "voom-normalized expression; CT-vs-ND DESeq2 Wald sig stars (ns/*/**/***)",
    genes_requested = panel5d_genes, genes_rendered = info_5d$genes,
    genes_missing = info_5d$missing, n_sig_stars = info_5d$n_stars,
    input_table = "deg/voom_BLAST_primary.rds $v$E (expression); deg/deseq2_BLAST_primary.rds (sig)"))
cat("  -> blastema_anti_fibrotic_trajectories_v10: ", length(info_5d$genes),
    " genes, ", info_5d$n_stars, " stars\n", sep = "")

cat("\n[", format(Sys.time(), "%FT%T%z"), "] t19 Round 8 done\n")
