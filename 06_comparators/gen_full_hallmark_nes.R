#!/usr/bin/env Rscript
# v9 Group C fix — generate full 50-Hallmark NES TSVs for X. laevis, Descalzi, Erdei.
# Outputs:
#   pathway/hallmark_nes_BRAIN_24HPT_full.tsv  (X. laevis 24 HPT, voom_BRAIN_primary, all Hallmarks)
#   pathway/hallmark_nes_BRAIN_7DPT_full.tsv   (X. laevis 7 DPT,  voom_BRAIN_primary, all Hallmarks)
#   comparators/descalzi/gsea_<region>_full.tsv  (mPFC, NAc, PAG; fgsea, no filter)
#   comparators/erdei/gsea_hallmark_<contrast>_full.tsv  (bi_vs_nt, fa_vs_nt, bi_vs_fa; fgsea, no filter)
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(fgsea); library(msigdbr)
})

PROJ <- "/project/def-srussell/nirosha_xeno"

# ---- 1. X. laevis full Hallmark NES from gseaResult @result ----
RDS <- file.path(PROJ, "pathway/all_enrichment_results.rds")
r <- readRDS(RDS)
cat("all_enrichment_results.rds: ", length(r), " entries\n")

extract_xen <- function(slot_name, out_tag, out_path) {
  stopifnot(slot_name %in% names(r))
  x <- r[[slot_name]]
  df <- if (isS4(x) && "result" %in% slotNames(x)) as.data.frame(x@result)
        else as.data.frame(x)
  out <- df %>%
    transmute(term = ID, NES = NES, pvalue = pvalue, padj = p.adjust)
  write_tsv(out, out_path)
  cat(sprintf("Wrote %s : %d Hallmarks (%s)\n",
              out_path, nrow(out), out_tag))
  invisible(out)
}
extract_xen("voom_BRAIN_primary_CT_vs_ND_24HPA_hallmark",
            "X. laevis 24 HPT BRAIN per-homeolog",
            file.path(PROJ, "pathway/hallmark_nes_BRAIN_24HPT_full.tsv"))
extract_xen("voom_BRAIN_primary_CT_vs_ND_7DPA_hallmark",
            "X. laevis 7 DPT BRAIN per-homeolog",
            file.path(PROJ, "pathway/hallmark_nes_BRAIN_7DPT_full.tsv"))

# ---- 2. Mouse Hallmark gene sets ----
mm_h <- msigdbr(species = "Mus musculus", category = "H")
gs_list <- split(mm_h$gene_symbol, mm_h$gs_name)
gs_list <- lapply(gs_list, unique)
cat(sprintf("Mouse Hallmark gene-set list: %d sets, total %d (unique) symbols\n",
            length(gs_list), length(unique(unlist(gs_list)))))

set.seed(42)

run_fgsea <- function(rank_df, out_path, label) {
  rank_vec <- setNames(rank_df$rank, rank_df$gene)
  rank_vec <- rank_vec[!is.na(rank_vec)]
  rank_vec <- rank_vec[!duplicated(names(rank_vec))]
  rank_vec <- sort(rank_vec, decreasing = TRUE)
  cat(sprintf("[%s] ranked genes: %d (range %+.3f .. %+.3f)\n",
              label, length(rank_vec), min(rank_vec), max(rank_vec)))
  res <- fgseaMultilevel(pathways = gs_list, stats = rank_vec,
                         minSize = 5, maxSize = 5000,
                         eps = 0, nproc = 1)
  out <- as_tibble(res) %>%
    transmute(term = pathway, NES, pvalue = pval, padj) %>%
    arrange(term)
  write_tsv(out, out_path)
  cat(sprintf("Wrote %s : %d Hallmarks\n", out_path, nrow(out)))
  invisible(out)
}

# ---- 3. Descalzi re-GSEA (t-stat ranking) ----
for (region in c("mPFC","NAc","PAG")) {
  de  <- read_tsv(file.path(PROJ, sprintf("comparators/descalzi/de_%s.tsv", region)),
                  show_col_types = FALSE)
  rdf <- de %>% transmute(gene = gene, rank = t)
  run_fgsea(rdf,
            file.path(PROJ, sprintf("comparators/descalzi/gsea_%s_full.tsv", region)),
            sprintf("Descalzi %s (t-stat)", region))
}

# ---- 4. Erdei re-GSEA (logFC ranking; n=1 per condition so no t-stat) ----
for (con in c("bi_vs_nt","fa_vs_nt","bi_vs_fa")) {
  de  <- read_tsv(file.path(PROJ, sprintf("comparators/erdei/de_%s.tsv", con)),
                  show_col_types = FALSE)
  rdf <- de %>% transmute(gene = gene, rank = logFC)
  run_fgsea(rdf,
            file.path(PROJ, sprintf("comparators/erdei/gsea_hallmark_%s_full.tsv", con)),
            sprintf("Erdei %s (logFC)", con))
}

cat("DONE.\n")
