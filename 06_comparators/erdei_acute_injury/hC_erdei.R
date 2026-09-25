#!/usr/bin/env Rscript
# H-C comparator: Erdei 2025 IJMS mouse S1 cortex bulk RNA-seq (PRJNA1137110, cortex subset).
# DEC-016 — re-engaged after erroneous drop on 2026-04-30.
# 3 samples × 3 conditions × 1 rep:
#   SRX25371062 CTX_nt (non-treated)
#   SRX25371064 CTX_bi (burn-injury, 1 h post-burn)
#   SRX25371066 CTX_fa (formalin, 1 h post-formalin)
# Primary contrast: bi vs nt (the classical pain-vs-baseline analogue).
# Secondary: bi vs fa (specificity of burn vs alternative noxious stimulus).
#
# Cross-species concordance: Xenopus 24HPA CT brain (NOT 7DPA — Erdei is acute
# 1h, matches our acute timepoint).
#
# n=1 per condition is power-limited; we use limma-voom with no replicate
# (formula ~condition, single-pair contrast). Results are descriptive at the
# gene level and directional at the GSEA/Hallmark level. We frame §2.10
# accordingly: pathway-level cross-species alignment with Xenopus 24HPA CT.

suppressPackageStartupMessages({
    library(tximport)
    library(edgeR)
    library(limma)
    library(fgsea)
    library(msigdbr)
    library(data.table)
    library(ggplot2)
})

PROJ <- "/project/def-srussell/nirosha_xeno"
DEG_DIR     <- file.path(PROJ, "deg")
PATHWAY_DIR <- file.path(PROJ, "pathway")
WORK <- "/scratch/srussell/srussell/nirosha_xeno_comparators/erdei/work"
QUANT <- file.path(WORK, "salmon")
FIG  <- file.path(WORK, "figures")
REF  <- "/scratch/srussell/srussell/refs/mm39_gencode_m33"
TX2GENE <- file.path(REF, "tx2gene.gencode.vM33.tsv")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. tximport: load salmon quant for the 3 cortex samples
# ------------------------------------------------------------
samples <- data.frame(
    sample    = c("CTX_nt", "CTX_bi", "CTX_fa"),
    srx       = c("SRX25371062", "SRX25371064", "SRX25371066"),
    condition = c("nt", "bi", "fa"),
    stringsAsFactors = FALSE
)
samples$quant_path <- file.path(QUANT, samples$sample, "quant.sf")
stopifnot(all(file.exists(samples$quant_path)))

# Load tx2gene (GENCODE: tx_id, gene_id, gene_symbol)
t2g <- fread(TX2GENE, header = FALSE, col.names = c("tx_id", "gene_id", "gene_symbol"))
# Use gene_symbol as identifier for cross-species comparison. Strip ENSMUSG version.
t2g[, tx_id := gsub("\\.\\d+$", "", tx_id)]
t2g[, gene_id := gsub("\\.\\d+$", "", gene_id)]
# tx2gene map for tximport: tx_id -> gene_symbol (preferred) with fallback to gene_id
t2g[, gene_label := ifelse(gene_symbol == "" | is.na(gene_symbol), gene_id, gene_symbol)]
tx2gene_for_tximport <- as.data.frame(t2g[, .(tx_id, gene_label)])

cat("[", format(Sys.time()), "] Running tximport on", nrow(samples), "samples...\n")
txi <- tximport(samples$quant_path, type = "salmon", tx2gene = tx2gene_for_tximport,
                ignoreAfterBar = TRUE, ignoreTxVersion = TRUE,
                countsFromAbundance = "lengthScaledTPM")
counts <- txi$counts
colnames(counts) <- samples$sample
cat("Counts matrix:", nrow(counts), "genes ×", ncol(counts), "samples\n")
cat("Median counts/sample:", round(median(colSums(counts))/1e6, 2), "M\n")

# ------------------------------------------------------------
# 2. n=1 descriptive DE (no replicate -> 0 residual df, can't run eBayes)
# ------------------------------------------------------------
# Erdei 2025 cortex: 1 sample per condition. Standard limma-voom + eBayes
# fails ("No residual degrees of freedom in linear model fits"). Use
# log2-CPM ratios as descriptive logFC; rank by logFC for fgsea. This is
# the appropriate methodology for unreplicated designs (cf. limma user
# guide §17 "DE without replicates"). Cross-species concordance and
# pathway-level GSEA remain valid; gene-level inference is suspended.
dge <- DGEList(counts = round(counts))
# Filter on expression in at least one sample; with n=1/condition there is
# no group-aware filter. Use a simple cpm threshold.
keep <- rowSums(cpm(dge) > 1) >= 1
cat("Genes kept (cpm>1 in >=1 sample):", sum(keep), "/", nrow(counts), "\n")
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

logcpm <- cpm(dge, log = TRUE, prior.count = 2)  # log2(CPM + prior)
# Sample order in logcpm columns matches samples$sample after tximport
stopifnot(all(colnames(logcpm) == samples$sample))

contrast_pairs <- list(
    bi_vs_nt = c("CTX_bi", "CTX_nt"),
    fa_vs_nt = c("CTX_fa", "CTX_nt"),
    bi_vs_fa = c("CTX_bi", "CTX_fa")
)

de_results <- list()
for (ctr in names(contrast_pairs)) {
    pair <- contrast_pairs[[ctr]]
    logFC <- logcpm[, pair[1]] - logcpm[, pair[2]]
    AveExpr <- (logcpm[, pair[1]] + logcpm[, pair[2]]) / 2
    de <- data.table(gene = names(logFC),
                     logFC = as.numeric(logFC),
                     AveExpr = as.numeric(AveExpr),
                     contrast = ctr)
    setorder(de, -logFC)  # ranked for fgsea
    fwrite(de, file.path(WORK, paste0("de_", ctr, ".tsv")), sep = "\t")
    de_results[[ctr]] <- de
    cat(ctr, ": top up:", head(de$gene, 5),
        "/ top down:", tail(de$gene, 5), "\n")
}

# ------------------------------------------------------------
# 3. fgsea on Hallmark + GO:BP (mouse)
# ------------------------------------------------------------
H_msig <- msigdbr(species = "Mus musculus", category = "H")
H_pw <- split(H_msig$gene_symbol, H_msig$gs_name)

# GO:BP — use C5 GO:BP subcategory
GOBP_msig <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
GOBP_pw <- split(GOBP_msig$gene_symbol, GOBP_msig$gs_name)

gsea_results <- list()
for (ctr in names(de_results)) {
    de <- de_results[[ctr]]
    # Rank by raw logFC (no p-values available at n=1)
    rank_metric <- de$logFC
    names(rank_metric) <- de$gene
    rank_metric <- sort(rank_metric, decreasing = TRUE)
    rank_metric <- rank_metric[!is.na(rank_metric) & is.finite(rank_metric)]

    for (db_name in c("hallmark", "gobp")) {
        pw <- if (db_name == "hallmark") H_pw else GOBP_pw
        cat(ctr, db_name, ": running fgsea on", length(rank_metric), "ranked genes vs", length(pw), "pathways\n")
        fgsea_res <- tryCatch({
            fgsea(pathways = pw, stats = rank_metric, eps = 0,
                  minSize = 15, maxSize = 500)
        }, error = function(e) {
            cat("  fgsea error:", conditionMessage(e), "\n"); NULL
        })
        if (!is.null(fgsea_res) && nrow(fgsea_res) > 0) {
            fgsea_res[, contrast := ctr]
            fgsea_res[, db := db_name]
            # leadingEdge is a list-col; collapse for tsv
            fgsea_res[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
            fwrite(fgsea_res, file.path(WORK, paste0("gsea_", db_name, "_", ctr, ".tsv")), sep = "\t")
            gsea_results[[paste0(ctr, "_", db_name)]] <- fgsea_res
            cat("  q<0.05:", sum(fgsea_res$padj < 0.05), "; q<0.25:", sum(fgsea_res$padj < 0.25), "\n")
        }
    }
}

# ------------------------------------------------------------
# 4. Cross-species concordance: Xenopus 24HPA CT brain vs mouse 1h cortex burn
# ------------------------------------------------------------
xen_de_path <- file.path(DEG_DIR, "voom_BRAIN_primary.rds")
xen_table <- NULL
if (file.exists(xen_de_path)) {
    xen_de <- readRDS(xen_de_path)
    if (is.list(xen_de) && "results" %in% names(xen_de)) {
        nm <- names(xen_de$results)
        target <- grep("CT_vs_ND.*24HPA|24HPA.*CT_vs_ND", nm, value = TRUE, ignore.case = TRUE)
        if (length(target) >= 1) {
            xt <- xen_de$results[[target[1]]]
            xen_table <- data.table(gene = rownames(xt), logFC = xt$logFC,
                                    P = xt$P.Value, adj.P = xt$adj.P.Val)
            xen_table[, gene_clean := toupper(gsub("\\.[LS]$", "", gene))]
            xen_table[, gene_clean := gsub("\\.\\d+$", "", gene_clean)]  # strip homeolog numeric suffix
            cat("Xenopus 24HPA CT BRAIN signature loaded:", nrow(xen_table), "genes\n")
        } else {
            cat("WARN: no CT_vs_ND_24HPA in voom_BRAIN_primary$results (names: ", paste(nm, collapse = ", "), ")\n")
        }
    }
}

if (!is.null(xen_table) && length(de_results) > 0) {
    primary_de <- de_results[["bi_vs_nt"]]
    primary_de[, gene_upper := toupper(gene)]
    m <- merge(xen_table[, .(gene_clean, xen_logFC = logFC, xen_p = P, xen_q = adj.P)],
               primary_de[, .(gene_upper, mouse_logFC = logFC)],
               by.x = "gene_clean", by.y = "gene_upper")
    if (nrow(m) > 100) {
        cat("\n=== Cross-species concordance (Xenopus 24HPA CT BRAIN vs Erdei mouse cortex bi_vs_nt) ===\n")
        cat("Overlap n =", nrow(m), "\n")
        ct <- cor.test(m$xen_logFC, m$mouse_logFC, method = "spearman")
        cat("Spearman ρ =", round(ct$estimate, 3), "p =", format(ct$p.value, digits = 3), "\n")
        sig_xen <- m[xen_p < 0.05]
        if (nrow(sig_xen) > 0) {
            agree <- sum(sign(sig_xen$xen_logFC) == sign(sig_xen$mouse_logFC))
            cat("Direction agreement among Xen-sig genes:", agree, "/", nrow(sig_xen),
                " (", round(100 * agree / nrow(sig_xen), 1), "%)\n")
            conc_dt <- data.table(
                comparison = "xenopus_24HPA_CT_BRAIN_vs_erdei_bi_vs_nt",
                n_overlap = nrow(m),
                spearman_rho = ct$estimate, spearman_p = ct$p.value,
                xen_sig_n = nrow(sig_xen), agree = agree,
                pct_agree = round(100 * agree / nrow(sig_xen), 1)
            )
            fwrite(conc_dt, file.path(WORK, "xenopus_concordance.tsv"), sep = "\t")
        }
    }
}

# ------------------------------------------------------------
# 5. Hallmark NES correlation: Erdei vs Xenopus brain
# ------------------------------------------------------------
xen_gsea_path <- file.path(PATHWAY_DIR, "all_enrichment_summary.tsv")
if (file.exists(xen_gsea_path) &&
    "bi_vs_nt_hallmark" %in% names(gsea_results)) {
    xen_gsea <- fread(xen_gsea_path)
    xen_h <- xen_gsea[grepl("BRAIN", contrast, ignore.case = TRUE) &
                      grepl("CT_vs_ND.*24HPA|24HPA.*CT_vs_ND", contrast, ignore.case = TRUE) &
                      grepl("hallmark", db, ignore.case = TRUE)]
    if (nrow(xen_h) > 0) {
        mouse_h <- gsea_results[["bi_vs_nt_hallmark"]]
        # GENCODE Hallmark IDs match between species (msigdbr Hallmark gs_name is identical)
        joined <- merge(
            xen_h[, .(pathway = ID, xen_NES = NES)],
            mouse_h[, .(pathway, mouse_NES = NES)],
            by = "pathway"
        )
        if (nrow(joined) >= 5) {
            ct <- cor.test(joined$xen_NES, joined$mouse_NES, method = "spearman")
            cat("\nHallmark NES correlation: Xen-BRAIN-24HPA-CT vs Erdei-bi_vs_nt: ρ=",
                round(ct$estimate, 3), " p=", format(ct$p.value, digits = 3),
                " (n=", nrow(joined), ")\n", sep = "")
            fwrite(joined, file.path(WORK, "hallmark_nes_xen_vs_erdei.tsv"), sep = "\t")
        }
    }
}

# ------------------------------------------------------------
# 6. Volcano figures (bi_vs_nt + bi_vs_fa)
# ------------------------------------------------------------
for (ctr in c("bi_vs_nt", "bi_vs_fa")) {
    if (!is.null(de_results[[ctr]])) {
        de <- de_results[[ctr]]
        de[, abs_lfc := abs(logFC)]
        de[, top := abs_lfc > quantile(abs_lfc, 0.99, na.rm = TRUE)]
        pdf(file.path(FIG, paste0("rankplot_", ctr, ".pdf")), width = 6, height = 5)
        # MA-style: AveExpr vs logFC. n=1 so no p-value coloring.
        print(ggplot(de, aes(AveExpr, logFC)) +
                  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
                  geom_point(alpha = 0.3, size = 0.5, color = "grey40") +
                  geom_point(data = de[top == TRUE], color = "red", size = 0.8, alpha = 0.6) +
                  theme_minimal() +
                  labs(title = paste0("Erdei mouse cortex ", ctr, " (n=1/condition; descriptive log2-CPM ratio)"),
                       x = "Mean log2-CPM", y = "log2 fold-change"))
        dev.off()
    }
}

# ------------------------------------------------------------
# 7. Mirror outputs to /project (small TSVs and PDFs only — DEC-014 storage policy)
# ------------------------------------------------------------
PROJ_OUT <- file.path(PROJ, "comparators", "erdei")
dir.create(PROJ_OUT, recursive = TRUE, showWarnings = FALSE)
for (f in list.files(WORK, pattern = "\\.tsv$", full.names = TRUE)) {
    file.copy(f, PROJ_OUT, overwrite = TRUE)
}
PROJ_FIG <- file.path(PROJ, "figures", "stage_H", "erdei")
dir.create(PROJ_FIG, recursive = TRUE, showWarnings = FALSE)
for (f in list.files(FIG, pattern = "\\.pdf$", full.names = TRUE)) {
    file.copy(f, PROJ_FIG, overwrite = TRUE)
}

cat("[", format(Sys.time()), "] hC_erdei.R complete.\n", sep = "")
