#!/usr/bin/env Rscript
# t18_v10_brain_panels.R (Round 7 restyle)
# Round 6/iter-3 outputs are overwritten in place per Round 7 brief:
#   - simple titles (no engine parenthetical)
#   - no subtitle (or tiny n_sig line)
#   - up/down direction colours + discrete legend on all volcanos
#   - bigger dots
#   - x-axis tightened to ~±5 (with off-scale diamond markers)
#   - Hallmark dot plots show ALL sig at q<0.05 (no top 10 cap)
# Plus 2 new temporal contrasts: ND 7v24 and CT 7v24 brain volcanos.

suppressPackageStartupMessages({
    library(dplyr); library(readr); library(tibble); library(ggplot2);
    library(ggrepel); library(yaml); library(grid)
})

PROJ    <- "/project/def-srussell/nirosha_xeno"
OUT     <- file.path(PROJ, "figures_current/v10_pulls/brain")
DEG     <- file.path(PROJ, "deg")
PATHRDS <- file.path(PROJ, "pathway/all_enrichment_results.rds")
WGCNA   <- file.path(PROJ, "wgcna")
MOFA    <- file.path(PROJ, "mofa")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

cat("[", format(Sys.time(), "%FT%T%z"), "] v10 brain panel render (R7) start ->", OUT, "\n")

panel_theme <- theme_classic(base_size = 11) +
    theme(axis.text.x  = element_text(size = 9.5, color = "black"),
          axis.text.y  = element_text(size = 9,   color = "black"),
          axis.title   = element_text(size = 10,  color = "black"),
          plot.title   = element_text(size = 11,  face = "plain", family = "sans", hjust = 0),
          plot.subtitle= element_text(size = 8,   color = "grey35"),
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
    info$script <- "scripts/t18_v10_brain_panels.R"
    write_yaml(info, file.path(OUT, paste0(stem, "_meta.yaml")))
}

# ============================================================
# Volcano helper (with up/down colouring + legend, narrower clip)
# ============================================================
make_volcano <- function(df, stem, title,
                         up_label = "Up in CT", down_label = "Down in CT",
                         top_n_labels = 10, clip_x = c(-5, 5),
                         color_up = "#B2182B", color_down = "#2166AC",
                         x_axis_label = "log\u2082 FC (CT vs ND)") {
    df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
    df$neg_log_padj <- -log10(df$padj)
    df$direction <- factor(
        ifelse(df$padj < 0.05 & df$log2FoldChange > 0, "up",
        ifelse(df$padj < 0.05 & df$log2FoldChange < 0, "down", "ns")),
        levels = c("up", "down", "ns")
    )
    n_sig <- sum(df$direction != "ns")
    n_up  <- sum(df$direction == "up")
    n_dn  <- sum(df$direction == "down")

    # Top labels by |stat| within sig
    top <- df[df$direction != "ns", ]
    top <- top[order(-abs(top$stat))[seq_len(min(top_n_labels, nrow(top)))], ]
    # Ensure *gdf15.L* labelled on 7DPT CT-vs-ND
    if (grepl("7dpt", stem) && !grepl("7v24", stem) &&
        "gdf15.L" %in% df$gene_id && !"gdf15.L" %in% top$gene_id) {
        top <- rbind(top, df[df$gene_id == "gdf15.L", ])
    }

    df$display_lfc  <- pmax(pmin(df$log2FoldChange, clip_x[2]), clip_x[1])
    df$off_scale    <- df$log2FoldChange < clip_x[1] | df$log2FoldChange > clip_x[2]
    top$display_lfc <- pmax(pmin(top$log2FoldChange, clip_x[2]), clip_x[1])
    top$off_scale   <- top$log2FoldChange < clip_x[1] | top$log2FoldChange > clip_x[2]
    n_off <- sum(df$off_scale & df$direction != "ns")

    dir_palette <- c("up" = color_up, "down" = color_down, "ns" = "grey75")
    dir_labels  <- c("up" = up_label, "down" = down_label, "ns" = "n.s.")

    # Plot order: ns first (background), then sig on top
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
        labs(title = title,
             subtitle = sprintf("n = %d sig at FDR < 0.05%s",
                                n_sig,
                                if (n_off > 0) sprintf(" \u2022 %d off-scale (\u2666 clipped at \u00b1%g)",
                                                       n_off, clip_x[2]) else ""),
             x = x_axis_label,
             y = expression(-log[10]~"(padj)")) +
        panel_theme

    save_panel(p, stem)

    top_hits <- lapply(seq_len(nrow(top)), function(i) list(
        gene_id = top$gene_id[i],
        log2FC  = signif(top$log2FoldChange[i], 3),
        stat    = signif(top$stat[i], 3),
        padj    = signif(top$padj[i], 3)
    ))
    list(n_sig = n_sig, n_up = n_up, n_down = n_dn, top_hits = top_hits,
         n_total = nrow(df), n_off = n_off)
}

# ============================================================
# Load DESeq2 brain primary
# ============================================================
brain_pri <- readRDS(file.path(DEG, "deseq2_BRAIN_primary.rds"))

# Panel 1 — Brain volcano 7DPT (CT vs ND)
cat("\n[P1] brain_volcano_7dpt_deseq2\n")
info1 <- make_volcano(brain_pri$results[["CT_vs_ND_7DPA"]],
                      "brain_volcano_7dpt_deseq2",
                      "BRAIN CT vs ND at 7 DPT")
write_meta("brain_volcano_7dpt_deseq2", list(
    panel = "Fig 2A", title = "BRAIN CT vs ND at 7 DPT",
    analysis_engine = "DESeq2 Wald, per-homeolog",
    fdr_threshold = 0.05, n_sig = info1$n_sig,
    n_up = info1$n_up, n_down = info1$n_down,
    n_total = info1$n_total,
    input_table = "deg/deseq2_BRAIN_primary.rds (results$CT_vs_ND_7DPA)",
    top_hits = info1$top_hits))
cat(sprintf("  -> %d sig (%d up, %d down)%s\n",
            info1$n_sig, info1$n_up, info1$n_down,
            if (info1$n_off > 0) sprintf("; %d off-scale", info1$n_off) else ""))

# Panel 2 — Brain volcano 24HPT (CT vs ND)
cat("[P2] brain_volcano_24hpt_deseq2\n")
info2 <- make_volcano(brain_pri$results[["CT_vs_ND_24HPA"]],
                      "brain_volcano_24hpt_deseq2",
                      "BRAIN CT vs ND at 24 HPT")
write_meta("brain_volcano_24hpt_deseq2", list(
    panel = "Fig 2B", title = "BRAIN CT vs ND at 24 HPT",
    analysis_engine = "DESeq2 Wald, per-homeolog",
    fdr_threshold = 0.05, n_sig = info2$n_sig,
    n_up = info2$n_up, n_down = info2$n_down,
    n_total = info2$n_total,
    input_table = "deg/deseq2_BRAIN_primary.rds (results$CT_vs_ND_24HPA)",
    top_hits = info2$top_hits))
cat(sprintf("  -> %d sig (%d up, %d down)\n", info2$n_sig, info2$n_up, info2$n_down))

# Panel 3 — Brain LRT volcano (treatment × time, homeolog-sum)
cat("[P3] brain_lrt_trtxtime_deseq2\n")
lrt_obj <- readRDS(file.path(DEG, "deseq2_lrt_BRAIN_late_7v24_homeolog_sum.rds"))
lrt_df <- as.data.frame(lrt_obj$res)
lrt_df$gene_id <- rownames(lrt_df)
info3 <- make_volcano(lrt_df,
                      "brain_lrt_trtxtime_deseq2",
                      "BRAIN treatment \u00d7 time interaction",
                      up_label = "CT trajectory \u2191 vs ND",
                      down_label = "CT trajectory \u2193 vs ND",
                      color_up = "#1B7837", color_down = "#762A83",
                      x_axis_label = "log\u2082 FC interaction (CT_7DPA \u2212 CT_24HPA) \u2212 (ND_7DPA \u2212 ND_24HPA)")
write_meta("brain_lrt_trtxtime_deseq2", list(
    panel = "Fig 2C", title = "BRAIN treatment \u00d7 time interaction",
    analysis_engine = "DESeq2 LRT, homeolog-sum",
    fdr_threshold = 0.05, n_sig = info3$n_sig,
    n_up = info3$n_up, n_down = info3$n_down,
    n_total = info3$n_total,
    input_table = "deg/deseq2_lrt_BRAIN_late_7v24_homeolog_sum.rds",
    top_hits = info3$top_hits))
cat(sprintf("  -> %d sig (%d up, %d down)%s\n",
            info3$n_sig, info3$n_up, info3$n_down,
            if (info3$n_off > 0) sprintf("; %d off-scale", info3$n_off) else ""))

# Panel 9 (NEW) — Brain volcano ND 7v24 (temporal, ND only)
cat("[P9 NEW] brain_volcano_ND_7v24_deseq2\n")
info9 <- make_volcano(brain_pri$results[["ND_7DPA_vs_24HPA"]],
                      "brain_volcano_ND_7v24_deseq2",
                      "BRAIN ND: 7 DPT vs 24 HPT (temporal, no cocktail)",
                      up_label = "Up at 7 DPT", down_label = "Down at 7 DPT",
                      x_axis_label = "log\u2082 FC (7 DPT vs 24 HPT, ND animals)")
write_meta("brain_volcano_ND_7v24_deseq2", list(
    panel = "Round 7 new",
    title = "BRAIN ND: 7 DPT vs 24 HPT (temporal, no cocktail)",
    analysis_engine = "DESeq2 Wald, per-homeolog, ND-only contrast",
    fdr_threshold = 0.05, n_sig = info9$n_sig,
    n_up = info9$n_up, n_down = info9$n_down, n_total = info9$n_total,
    input_table = "deg/deseq2_BRAIN_primary.rds (results$ND_7DPA_vs_24HPA)",
    top_hits = info9$top_hits))
cat(sprintf("  -> %d sig (%d up, %d down)\n",
            info9$n_sig, info9$n_up, info9$n_down))

# Panel 10 (NEW) — Brain volcano CT 7v24 (temporal, CT only)
cat("[P10 NEW] brain_volcano_CT_7v24_deseq2\n")
info10 <- make_volcano(brain_pri$results[["CT_7DPA_vs_24HPA"]],
                       "brain_volcano_CT_7v24_deseq2",
                       "BRAIN CT: 7 DPT vs 24 HPT (temporal, cocktail-treated)",
                       up_label = "Up at 7 DPT", down_label = "Down at 7 DPT",
                       x_axis_label = "log\u2082 FC (7 DPT vs 24 HPT, CT animals)")
write_meta("brain_volcano_CT_7v24_deseq2", list(
    panel = "Round 7 new",
    title = "BRAIN CT: 7 DPT vs 24 HPT (temporal, cocktail-treated)",
    analysis_engine = "DESeq2 Wald, per-homeolog, CT-only contrast",
    fdr_threshold = 0.05, n_sig = info10$n_sig,
    n_up = info10$n_up, n_down = info10$n_down, n_total = info10$n_total,
    input_table = "deg/deseq2_BRAIN_primary.rds (results$CT_7DPA_vs_24HPA)",
    top_hits = info10$top_hits))
cat(sprintf("  -> %d sig (%d up, %d down)\n",
            info10$n_sig, info10$n_up, info10$n_down))

# ============================================================
# Hallmark dot plots — show ALL sig (no top 10 cap)
# ============================================================
enrich <- readRDS(PATHRDS)

make_hallmark_dot <- function(key, stem, title, time_label) {
    x <- enrich[[key]]
    if (is.null(x)) { cat("  [WARN] key not found:", key, "\n"); return(NULL) }
    df <- if (isS4(x) && "result" %in% slotNames(x)) as.data.frame(x@result) else as.data.frame(x)
    sig <- df[df$p.adjust < 0.05, ]
    n_sig <- nrow(sig)
    if (n_sig == 0) {
        cat("  [WARN] no sig Hallmarks for", key, "; showing top 10 by p\n")
        sig <- df[order(df$p.adjust), ][seq_len(min(10, nrow(df))), ]
    } else {
        # Round 7: show ALL sig at q<0.05, ordered by NES descending
        sig <- sig[order(-sig$NES), ]
    }
    sig$term_short <- sub("^HALLMARK_", "", sig$ID)
    sig$term_lbl   <- factor(sig$term_short, levels = rev(sig$term_short))
    sig$neg_log_q  <- -log10(sig$p.adjust)
    sig$set_size   <- vapply(strsplit(sig$core_enrichment, "/"), length, integer(1))

    p <- ggplot(sig, aes(x = NES, y = term_lbl, color = NES, size = neg_log_q)) +
        geom_vline(xintercept = 0, linewidth = 0.4, color = "grey70") +
        geom_point(stroke = 0.4, alpha = 0.9) +
        scale_color_gradient2(low = "#2166AC", mid = "#F0F0F0", high = "#B2182B",
                              midpoint = 0, name = "NES") +
        scale_size_continuous(range = c(2.5, 7), name = expression(-log[10]~q)) +
        labs(title = title,
             subtitle = sprintf("n = %d sig Hallmarks at q < 0.05", n_sig),
             x = sprintf("NES (CT vs ND, %s)", time_label), y = NULL) +
        panel_theme + theme(legend.position = "right")

    save_panel(p, stem, h = max(80, 8 * nrow(sig) + 40))

    list(panel_title = title, n_sig_total = n_sig,
         terms = lapply(seq_len(nrow(sig)), function(i) list(
             ID = sig$ID[i], NES = signif(sig$NES[i], 3),
             padj = signif(sig$p.adjust[i], 3), n_le = sig$set_size[i])),
         key = key)
}

cat("[P4] brain_hallmark_24hpt_deseq2 (ALL sig)\n")
info4 <- make_hallmark_dot("deseq2_BRAIN_primary_CT_vs_ND_24HPA_hallmark",
                            "brain_hallmark_24hpt_deseq2",
                            "BRAIN Hallmark NES at 24 HPT", "24 HPT")
if (!is.null(info4)) {
    write_meta("brain_hallmark_24hpt_deseq2", list(
        panel = "Fig 2D", title = info4$panel_title,
        analysis_engine = "clusterProfiler GSEA on signed -log10(p) * sign(logFC); DESeq2 Wald ranks",
        fdr_threshold = 0.05, n_sig = info4$n_sig_total,
        input_table = paste0("pathway/all_enrichment_results.rds [", info4$key, "]"),
        terms = info4$terms))
    cat(sprintf("  -> %d Hallmarks shown\n", info4$n_sig_total))
}

cat("[P5] brain_hallmark_7dpt_deseq2 (ALL sig)\n")
info5 <- make_hallmark_dot("deseq2_BRAIN_primary_CT_vs_ND_7DPA_hallmark",
                            "brain_hallmark_7dpt_deseq2",
                            "BRAIN Hallmark NES at 7 DPT", "7 DPT")
if (!is.null(info5)) {
    write_meta("brain_hallmark_7dpt_deseq2", list(
        panel = "Fig 2E", title = info5$panel_title,
        analysis_engine = "clusterProfiler GSEA on signed -log10(p) * sign(logFC); DESeq2 Wald ranks",
        fdr_threshold = 0.05, n_sig = info5$n_sig_total,
        input_table = paste0("pathway/all_enrichment_results.rds [", info5$key, "]"),
        terms = info5$terms))
    cat(sprintf("  -> %d Hallmarks shown\n", info5$n_sig_total))
}

# ============================================================
# Panel 6 — Synaptic GO BP bar at 7 DPT
# ============================================================
cat("[P6] brain_synaptic_go_7dpt_deseq2\n")
x6 <- enrich[["deseq2_BRAIN_primary_CT_vs_ND_7DPA_goBP"]]
df6 <- if (isS4(x6) && "result" %in% slotNames(x6)) as.data.frame(x6@result) else as.data.frame(x6)
syn_pat <- "synap|presynaptic|postsynaptic|neurotransmitter|axon[^o]|dendrit"
syn <- df6[grepl(syn_pat, df6$Description, ignore.case = TRUE), ]
syn <- syn[order(syn$p.adjust), ][seq_len(min(10, nrow(syn))), ]
syn$term_lbl <- factor(syn$Description, levels = rev(syn$Description))

x_range <- max(abs(syn$NES))
x_lim   <- c(-x_range * 1.20, x_range * 1.20)

p6 <- ggplot(syn, aes(x = NES, y = term_lbl, fill = NES)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
    geom_col(color = "grey25", linewidth = 0.2, width = 0.7) +
    geom_text(aes(x = 0,
                  label = sprintf("q = %.2g", p.adjust),
                  hjust = ifelse(NES > 0, 1.05, -0.05)),
              size = 2.7, color = "black") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "NES") +
    coord_cartesian(xlim = x_lim) +
    labs(title = "BRAIN synaptic GO BP NES at 7 DPT",
         subtitle = sprintf("n = %d synaptic-keyword terms", nrow(syn)),
         x = "NES (CT vs ND, 7 DPT)", y = NULL) +
    panel_theme + theme(legend.position = "right")
save_panel(p6, "brain_synaptic_go_7dpt_deseq2", h = max(85, 10 * nrow(syn) + 40))
write_meta("brain_synaptic_go_7dpt_deseq2", list(
    panel = "Fig 4A", title = "BRAIN synaptic GO BP NES at 7 DPT",
    analysis_engine = "clusterProfiler GSEA on GO BP; DESeq2 Wald ranks; synaptic-keyword filter",
    keyword_filter = syn_pat, fdr_threshold = 0.05, n_sig = nrow(syn),
    input_table = "pathway/all_enrichment_results.rds [deseq2_BRAIN_primary_CT_vs_ND_7DPA_goBP]",
    terms = lapply(seq_len(nrow(syn)), function(i) list(
        ID = syn$ID[i], Description = syn$Description[i],
        NES = signif(syn$NES[i], 3),
        padj = signif(syn$p.adjust[i], 3)))))
cat(sprintf("  -> %d synaptic terms shown\n", nrow(syn)))

# ============================================================
# Panel 7 — ME1 turquoise WGCNA x Hallmark over-rep
# ============================================================
cat("[P7] brain_me1_hallmark_overrep\n")
hm <- read_tsv(file.path(WGCNA, "brain_module_hallmark_enrichment.tsv"),
               show_col_types = FALSE)
turq <- hm[hm$module_color == "turquoise", ]
turq <- turq[order(turq$p.adjust), ][seq_len(min(5, nrow(turq))), ]
turq$term_short <- sub("^HALLMARK_", "", turq$ID)
turq$term_lbl   <- factor(turq$term_short, levels = rev(turq$term_short))
turq$neg_log_q  <- -log10(turq$p.adjust)

p7 <- ggplot(turq, aes(x = Count, y = term_lbl, fill = neg_log_q)) +
    geom_col(color = "grey25", linewidth = 0.2, width = 0.65) +
    geom_text(aes(label = sprintf("q=%.1e", p.adjust)),
              hjust = -0.1, size = 2.7, color = "black") +
    scale_fill_gradient(low = "#FEE0B6", high = "#7F3B08",
                        name = expression(-log[10]~q)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.32))) +
    labs(title = "ME1 turquoise (n = 1,372 genes): top 5 Hallmark over-rep",
         x = "Overlap count", y = NULL) +
    panel_theme + theme(legend.position = "right")
save_panel(p7, "brain_me1_hallmark_overrep", h = 80)
write_meta("brain_me1_hallmark_overrep", list(
    panel = "Fig 5*", title = "ME1 turquoise module Hallmark over-representation",
    analysis_engine = "Hypergeometric over-rep (precomputed)",
    fdr_threshold = 0.05, n_sig = sum(turq$p.adjust < 0.05),
    input_table = "wgcna/brain_module_hallmark_enrichment.tsv [filter module_color==turquoise]",
    module = "ME1 turquoise", module_n_genes = 1372,
    terms = lapply(seq_len(nrow(turq)), function(i) list(
        ID = turq$ID[i], Count = turq$Count[i],
        padj = signif(turq$p.adjust[i], 3),
        fold_enrich = signif(turq$FoldEnrichment[i], 3)))))
cat(sprintf("  -> %d Hallmarks shown\n", nrow(turq)))

# ============================================================
# Panel 8 — MOFA F3 BLAST loadings
# ============================================================
cat("[P8] mofa_factor3_blast_loadings\n")
ml <- read_tsv(file.path(MOFA, "factor_loadings_top.tsv"),
               show_col_types = FALSE)
f3b <- ml[ml$factor == "Factor3" & ml$view == "BLAST", ]
f3b$abs_w <- abs(f3b$weight)
f3b <- f3b[order(-f3b$abs_w), ][seq_len(min(10, nrow(f3b))), ]
if (nrow(f3b) > 0) {
    f3b$gene_clean <- sub("_BLAST$", "", f3b$gene)
    f3b$gene_lbl   <- factor(f3b$gene_clean, levels = rev(f3b$gene_clean))

    p8 <- ggplot(f3b, aes(x = weight, y = gene_lbl, fill = weight)) +
        geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
        geom_col(color = "grey25", linewidth = 0.2, width = 0.7) +
        scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                             midpoint = 0, name = "weight") +
        labs(title = "MOFA+ Factor 3: top 10 BLAST-view loadings",
             x = "Factor 3 weight (BLAST view)", y = NULL) +
        panel_theme + theme(legend.position = "right")
    save_panel(p8, "mofa_factor3_blast_loadings", h = 100)
    write_meta("mofa_factor3_blast_loadings", list(
        panel = "Fig 5* (optional)",
        title = "MOFA+ Factor 3 top 10 BLAST-view loadings",
        analysis_engine = "MOFA+ factor model (precomputed)",
        input_table = file.path(MOFA, "factor_loadings_top.tsv"),
        top_hits = lapply(seq_len(nrow(f3b)), function(i) list(
            gene_id = f3b$gene_clean[i],
            weight = signif(f3b$weight[i], 4)))))
    cat(sprintf("  -> %d loadings shown\n", nrow(f3b)))
}

# Summary
cat("\n========================================\n")
cat("Output files in", OUT, ":\n")
cat("========================================\n")
fs <- sort(list.files(OUT, full.names = FALSE))
for (f in fs) cat("  ", f, "\n")
cat("\n[", format(Sys.time(), "%FT%T%z"), "] t18 (R7) done\n")
