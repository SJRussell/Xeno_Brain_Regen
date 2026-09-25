#!/usr/bin/env Rscript
# Fig 6.e (v9) — X. laevis BRAIN bulk UCell against Walker 8-marker etv1+/glut signature.
# PER-HOMEOLOG PRIMARY (manuscript of record §4.6): use voom_BRAIN_primary.rds, not vsd_blind.
# Time labels: HPT/DPT (not HPA/DPA). 1:1. Plain 11pt sans title, ≤ 80 chars, no subtitle.
suppressPackageStartupMessages({
  library(UCell); library(dplyr); library(readr); library(tibble); library(ggplot2); library(yaml)
})

PROJ  <- "/project/def-srussell/nirosha_xeno"
VOOM  <- file.path(PROJ, "deg/voom_BRAIN_primary.rds")        # per-homeolog primary
META  <- file.path(PROJ, "deg/sample_metadata.tsv")
FIG_D <- file.path(PROJ, "figures_v9/fig6")
dir.create(FIG_D, showWarnings = FALSE, recursive = TRUE)

TS  <- format(Sys.time(), "%Y%m%dT%H%M%S")
LOG <- file(file.path(FIG_D, "fig6e_log.txt"), "w")
log <- function(...) { m <- paste0(..., collapse=""); cat(m, "\n"); cat(m, "\n", file = LOG, append = TRUE) }
log("== Fig 6.e (v9) X. laevis BRAIN UCell vs Walker etv1+ (per-homeolog primary) — ", TS, " ==")

WALKER <- c("etv1","slc17a7","cdh13","bcl11b","sox5","neto1","stk39","ntng1")

log("Loading voom_BRAIN_primary.rds (nested: list with $v EList, $fit MArrayLM, etc.) ...")
v_obj <- readRDS(VOOM)
# Structure is list(v=EList, fit=MArrayLM, results, summary, contrasts, design)
elist <- if (!is.null(v_obj$v)) v_obj$v else if (inherits(v_obj, "EList")) v_obj else stop("cannot locate EList")
expr  <- if (!is.null(elist$E)) elist$E else stop("EList missing $E")
rn    <- rownames(expr)
log("Expression matrix dim (genes x samples): ", paste(dim(expr), collapse=" x "))

# Resolve Walker markers to X. laevis homeologs (.L / .S)
resolve_xlaevis <- function(m) {
  patterns <- c(paste0("^", m, "\\.L$"), paste0("^", m, "\\.S$"),
                paste0("^", m, "\\.\\d+\\.[LS]$"), paste0("^", m, "$"),
                paste0("^", m, "[\\.\\-_]"))
  found <- character(0)
  for (pat in patterns) {
    hits <- rn[grepl(pat, rn, ignore.case = TRUE, perl = TRUE)]
    found <- unique(c(found, hits))
  }
  found
}
hits_per_marker <- setNames(lapply(WALKER, resolve_xlaevis), WALKER)
log("X. laevis homolog resolution per marker:")
for (m in WALKER) log("  ", m, " -> [", length(hits_per_marker[[m]]), "] ",
                      paste(hits_per_marker[[m]], collapse = ", "))

walker_signature <- unique(unlist(hits_per_marker, use.names = FALSE))
miss <- WALKER[sapply(hits_per_marker, length) == 0]
log("Total signature genes: ", length(walker_signature),
    "  (missing markers: ", paste(miss, collapse=", "), ")")

# BRAIN samples
meta <- read_tsv(META, show_col_types = FALSE) %>%
  filter(tissue == "BRAIN") %>%
  mutate(treatment = factor(treatment, levels = c("ND","CT")),
         time      = factor(time,      levels = c("24HPA","7DPA")),
         time_lbl  = factor(ifelse(time == "24HPA", "24 HPT", "7 DPT"),
                            levels = c("24 HPT","7 DPT")),
         rep       = factor(rep))
brain_ids <- meta$sample
expr_brain <- expr[, intersect(brain_ids, colnames(expr))]
log("BRAIN per-homeolog expression dim: ", paste(dim(expr_brain), collapse=" x "))

sig_list  <- list(walker_etv1_glut_regen = walker_signature)
ucell_mat <- ScoreSignatures_UCell(matrix = expr_brain, features = sig_list, ncores = 1)
log("UCell output dim: ", paste(dim(ucell_mat), collapse=" x "))

dat <- meta %>% mutate(ucell = ucell_mat[match(sample, rownames(ucell_mat)), 1])
log("Per-sample UCell scores:")
log(capture.output(print(dat %>% select(sample, treatment, time_lbl, rep, ucell) %>%
                          arrange(treatment, time_lbl, rep))))

# Linear model: trt x time interaction
fit  <- lm(ucell ~ treatment * time + rep, data = dat)
s    <- summary(fit)$coefficients
log("lm coefficients:"); log(capture.output(print(s)))

int_row <- grep("^treatmentCT:time7DPA$", rownames(s))
stopifnot(length(int_row) == 1)
beta <- s[int_row, "Estimate"]; se <- s[int_row, "Std. Error"]
tval <- s[int_row, "t value"];  pval <- s[int_row, "Pr(>|t|)"]
fdr  <- p.adjust(pval, method = "BH")

log(sprintf("trt x time interaction: beta = %+.4f, p = %.4g", beta, pval))

summary_tbl <- dat %>% group_by(treatment, time_lbl) %>%
  summarise(mean = mean(ucell), sem = sd(ucell)/sqrt(n()), n = n(), .groups = "drop")

write_tsv(summary_tbl, file.path(FIG_D, "fig6e_summary.tsv"))
write_tsv(dat %>% select(sample, treatment, time_lbl, rep, ucell),
          file.path(FIG_D, "fig6e_per_sample.tsv"))

TRT_COL <- c(ND = "#3B7DBF", CT = "#D9614C")

panel_theme <- theme_classic(base_size = 11) +
  theme(axis.text       = element_text(size = 9.5, color = "black"),
        axis.title      = element_text(size = 10,  color = "black"),
        plot.title      = element_text(size = 11,  face = "plain", family = "sans", hjust = 0),
        legend.text     = element_text(size = 9),
        legend.title    = element_text(size = 9.5),
        plot.margin     = margin(8, 14, 6, 8))

stat_lbl <- sprintf("trt x time beta = %+.3f, p = %.3f, FDR = %.3f",
                    beta, pval, fdr)
y_range  <- range(dat$ucell); y_pad <- diff(y_range) * 0.12

panel_e <- ggplot(dat, aes(x = time_lbl, y = ucell, fill = treatment, color = treatment)) +
  geom_boxplot(aes(group = interaction(time_lbl, treatment)),
               position = position_dodge(0.7), width = 0.5,
               linewidth = 0.45, outlier.shape = NA, alpha = 0.35) +
  geom_point(aes(group = interaction(time_lbl, treatment)),
             position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
             size = 2.0, alpha = 0.9) +
  annotate("text", x = 2.5, y = y_range[2] + y_pad * 0.6,
           label = stat_lbl, hjust = 1, vjust = 1, size = 3.0, color = "black") +
  coord_cartesian(xlim = c(0.6, 2.5), clip = "off") +
  scale_color_manual(values = TRT_COL, name = "Treatment") +
  scale_fill_manual(values  = TRT_COL, name = "Treatment") +
  labs(title = "X. laevis BRAIN UCell against Walker etv1+ signature",
       x = "Timepoint",
       y = "UCell score") +
  panel_theme

ggsave(file.path(FIG_D, "e_panel.pdf"), panel_e, width = 140, height = 130, units = "mm")
ggsave(file.path(FIG_D, "e_panel.png"), panel_e, width = 140, height = 130, units = "mm", dpi = 300)

key_claim <- sprintf("UCell trt x time interaction beta = %+.3f (p = %.3f, FDR = %.3f) on n = %d BRAIN samples; signature = %d X. laevis homeologs from %d Walker markers.",
                     beta, pval, fdr, nrow(dat), length(walker_signature), length(WALKER))

write_yaml(list(
  panel              = "Fig 6e",
  title              = "X. laevis BRAIN UCell against Walker etv1+ signature",
  axis_x_label       = "Timepoint",
  axis_y_label       = "UCell score",
  data_source_path   = VOOM,
  analysis_layer     = "per-homeolog (primary, manuscript of record per Methods §4.6)",
  n_obs              = nrow(dat),
  walker_markers_axolotl = as.list(WALKER),
  xlaevis_signature_genes = as.list(walker_signature),
  xlaevis_signature_n = length(walker_signature),
  walker_markers_with_no_homolog = if (length(miss) == 0) NULL else as.list(miss),
  model              = "ucell ~ treatment * time + rep",
  interaction        = list(beta = signif(beta, 4), se = signif(se, 4),
                            t = signif(tval, 4), p = signif(pval, 4), fdr = signif(fdr, 4)),
  key_numerical_claim = key_claim,
  group_summary      = lapply(seq_len(nrow(summary_tbl)), function(i) as.list(summary_tbl[i, ])),
  generated_utc      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  script             = "scripts/v9/fig6e_xlaevis_walker_etv1_ucell.R",
  styling            = list(width_mm = 140, height_mm = 130,
                            treatment_colors = list(ND = "#3B7DBF", CT = "#D9614C"),
                            time_axis_labels = list("24 HPT","7 DPT"))
), file.path(FIG_D, "panel_meta_e.yaml"))

log("Fig 6.e render COMPLETE.")
close(LOG)
