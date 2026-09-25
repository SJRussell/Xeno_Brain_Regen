#!/usr/bin/env Rscript
# Stage G — Sub-step (a): MOFA+ on n=12 paired blastema+brain.
#
# Per DEC-012:
#   - Two views (BRAIN, BLAST), VST-transformed expression, samples = paired
#     animal × treatment × timepoint at 24HPA + 7DPA only (11HPA BLAST excluded —
#     no BRAIN match by design).
#   - Top 8000 variable genes per view (matches WGCNA convention from Stage F).
#   - Latent factors K = 10 (drop degenerate at convergence).
#   - Outputs: factor variance per view, factor scores by sample-pair, factor
#     loadings (top 200 by |weight| per view per factor), trained model HDF5.
#
# Input:
#   /project/.../deg/dgelist_filtered_tmm.rds
#   /project/.../deg/vsd_blind.rds
#   /project/.../deg/sample_metadata.rds
# Output dir:  /scratch/.../nirosha_xeno_stage_G/mofa/
#   mofa_model.hdf5
#   factor_variance_explained.tsv
#   factor_scores.tsv
#   factor_loadings_top.tsv
#   factor_treatment_associations.tsv
#   figures/{factor_var_explained,factor_score_heatmap,factor_loadings_volcano}.pdf

suppressPackageStartupMessages({
    library(MOFA2)
    library(DESeq2)
    library(SummarizedExperiment)
    library(matrixStats)
    library(data.table)
    library(ggplot2)
})

PROJ <- "/project/def-srussell/nirosha_xeno"
DEG_DIR <- file.path(PROJ, "deg")
WORK <- "/scratch/srussell/srussell/nirosha_xeno_stage_G/mofa"
FIG  <- file.path(WORK, "figures")
dir.create(WORK, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

cat("[", format(Sys.time()), "] Loading inputs...\n", sep = "")
y      <- readRDS(file.path(DEG_DIR, "dgelist_filtered_tmm.rds"))
vsd    <- readRDS(file.path(DEG_DIR, "vsd_blind.rds"))
meta   <- readRDS(file.path(DEG_DIR, "sample_metadata.rds"))

stopifnot(all(colnames(vsd) == meta$sample) ||
          all(colnames(vsd) == rownames(meta)))
# vsd colnames have suffix like "_S3" — use meta$sample column (matches vsd colnames exactly)

vst_mat <- assay(vsd)  # genes × all 30 samples

# ------------------------------------------------------------
# Build paired sample table at matched (animal × treatment × timepoint)
# ------------------------------------------------------------
# Restrict to 24HPA + 7DPA timepoints — 11HPA BLAST has no BRAIN match.
meta$pair_id <- paste(meta$rep, meta$time, meta$treatment, sep = "_")
keep_time <- meta$time %in% c("24HPA", "7DPA")
meta_paired <- meta[keep_time, ]
# Identify pair_ids that have BOTH BRAIN and BLAST samples (the paired set).
pair_tab <- table(meta_paired$pair_id, meta_paired$tissue)
paired_ids <- rownames(pair_tab)[rowSums(pair_tab > 0) == 2]
cat("Paired pair_ids (have both BRAIN+BLAST at 24HPA/7DPA):", length(paired_ids), "\n")
print(paired_ids)

# n=12 expected: 3 reps × 2 treatments × 2 timepoints
if (length(paired_ids) != 12) {
    cat("WARN: expected 12 paired pair_ids, got", length(paired_ids), "\n")
    cat("     Continuing with whatever is paired; document n in DEC-012 followup.\n")
}

meta_paired <- meta_paired[meta_paired$pair_id %in% paired_ids, ]
meta_paired <- meta_paired[order(meta_paired$pair_id, meta_paired$tissue), ]
cat("Paired samples (rows in meta_paired):", nrow(meta_paired), "\n")

# ------------------------------------------------------------
# Build per-view matrices: features × samples (samples = pair_id, indexed by tissue)
# ------------------------------------------------------------
build_view <- function(tissue_label, top_n = 8000) {
    sel <- meta_paired$tissue == tissue_label
    sample_ids <- meta_paired$sample[sel]
    pair_ids_v <- meta_paired$pair_id[sel]
    stopifnot(length(sample_ids) == length(paired_ids))

    # Order columns by pair_id for consistent sample ordering across views
    ord <- match(paired_ids, pair_ids_v)
    sample_ids <- sample_ids[ord]
    pair_ids_v <- pair_ids_v[ord]

    M <- vst_mat[, sample_ids, drop = FALSE]
    colnames(M) <- pair_ids_v

    # Top-N variable genes within this view
    vars <- rowVars(M)
    top  <- order(vars, decreasing = TRUE)[seq_len(min(top_n, nrow(M)))]
    M[top, , drop = FALSE]
}

view_BRAIN <- build_view("BRAIN", top_n = 8000)
view_BLAST <- build_view("BLAST", top_n = 8000)
cat("BRAIN view:", dim(view_BRAIN)[1], "genes ×", dim(view_BRAIN)[2], "samples\n")
cat("BLAST view:", dim(view_BLAST)[1], "genes ×", dim(view_BLAST)[2], "samples\n")
stopifnot(identical(colnames(view_BRAIN), colnames(view_BLAST)))

# ------------------------------------------------------------
# MOFA+ object construction
# ------------------------------------------------------------
data_list <- list(BRAIN = view_BRAIN, BLAST = view_BLAST)
MOFAobject <- create_mofa(data_list)

# Sample-level metadata for downstream interpretation
samp_meta <- meta_paired[!duplicated(meta_paired$pair_id), c("pair_id", "rep", "time", "treatment")]
samp_meta <- as.data.frame(samp_meta)
samp_meta <- samp_meta[match(colnames(view_BRAIN), samp_meta$pair_id), ]
samp_meta$sample <- samp_meta$pair_id
samples_metadata(MOFAobject) <- samp_meta

# ------------------------------------------------------------
# Model options
# ------------------------------------------------------------
data_opts  <- get_default_data_options(MOFAobject)
data_opts$scale_views <- TRUE  # equalize variance contribution between BRAIN & BLAST

model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 10
model_opts$likelihoods <- c(BRAIN = "gaussian", BLAST = "gaussian")
model_opts$ard_factors <- TRUE   # automatic relevance determination — drops empty factors
model_opts$ard_weights <- TRUE
model_opts$spikeslab_weights <- TRUE  # sparsity in loadings

train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "medium"
train_opts$seed <- 20260430
train_opts$drop_factor_threshold <- 0.01  # drop factor if < 1% variance in any view
train_opts$verbose <- TRUE
train_opts$maxiter <- 2000

MOFAobject <- prepare_mofa(MOFAobject,
                           data_options = data_opts,
                           model_options = model_opts,
                           training_options = train_opts)

# ------------------------------------------------------------
# Train (writes HDF5)
# ------------------------------------------------------------
hdf5_out <- file.path(WORK, "mofa_model.hdf5")
cat("[", format(Sys.time()), "] Training MOFA+...\n", sep = "")
MOFAobject <- run_mofa(MOFAobject, outfile = hdf5_out, use_basilisk = TRUE)
cat("[", format(Sys.time()), "] Training complete. Model written to", hdf5_out, "\n", sep = "")

# ------------------------------------------------------------
# Extract & write outputs
# ------------------------------------------------------------
# 1. Variance explained per factor per view
ve <- get_variance_explained(MOFAobject)
ve_per_factor <- ve$r2_per_factor[[1]]  # group 1 (single group)
ve_dt <- as.data.table(ve_per_factor, keep.rownames = "factor")
fwrite(ve_dt, file.path(WORK, "factor_variance_explained.tsv"), sep = "\t")
cat("Factor variance explained per view:\n"); print(ve_dt)

# 2. Factor scores per sample
Z <- get_factors(MOFAobject)[[1]]  # group 1
Z_dt <- as.data.table(Z, keep.rownames = "pair_id")
Z_dt <- merge(Z_dt, samp_meta[, c("pair_id", "rep", "time", "treatment")], by = "pair_id")
fwrite(Z_dt, file.path(WORK, "factor_scores.tsv"), sep = "\t")

# 3. Top factor loadings per factor per view (top 200 by |weight|)
W_list <- get_weights(MOFAobject)
top_load <- list()
for (v in names(W_list)) {
    W <- W_list[[v]]
    for (f in colnames(W)) {
        ord <- order(abs(W[, f]), decreasing = TRUE)[1:200]
        top_load[[length(top_load) + 1]] <- data.table(
            view   = v,
            factor = f,
            gene   = rownames(W)[ord],
            weight = W[ord, f]
        )
    }
}
top_load_dt <- rbindlist(top_load)
fwrite(top_load_dt, file.path(WORK, "factor_loadings_top.tsv"), sep = "\t")

# 4. Treatment & time × factor associations (Wald test on factor scores)
assoc <- list()
for (f in colnames(Z)) {
    df <- data.frame(
        z = Z[, f],
        treatment = factor(samp_meta$treatment, levels = c("ND", "CT")),
        time      = factor(samp_meta$time, levels = c("24HPA", "7DPA")),
        rep       = samp_meta$rep
    )
    fit <- tryCatch(lm(z ~ treatment * time + rep, data = df),
                    error = function(e) NULL)
    if (is.null(fit)) next
    cf  <- summary(fit)$coefficients
    for (term in rownames(cf)) {
        assoc[[length(assoc) + 1]] <- data.table(
            factor = f,
            term   = term,
            estimate = cf[term, 1],
            se       = cf[term, 2],
            t        = cf[term, 3],
            p        = cf[term, 4]
        )
    }
}
assoc_dt <- rbindlist(assoc)
assoc_dt[, fdr := p.adjust(p, method = "BH")]
fwrite(assoc_dt, file.path(WORK, "factor_treatment_associations.tsv"), sep = "\t")

# 5. Quick figures (each wrapped in tryCatch — MOFA2 plot helpers can hit S7/ggplot2
#    compatibility issues on some R versions; failure here must not block the TSV outputs
#    or the file.copy step at the end of the script.)
tryCatch({
    pdf(file.path(FIG, "factor_var_explained.pdf"), width = 7, height = 5)
    print(plot_variance_explained(MOFAobject, plot_total = TRUE)[[1]])
    print(plot_variance_explained(MOFAobject)[[1]])
    dev.off()
}, error = function(e) {
    cat("WARN factor_var_explained.pdf:", conditionMessage(e), "\n")
    try(dev.off(), silent = TRUE)
})

tryCatch({
    pdf(file.path(FIG, "factor_score_heatmap.pdf"), width = 8, height = 6)
    print(plot_factor(MOFAobject, factors = "all", color_by = "treatment", shape_by = "time"))
    dev.off()
}, error = function(e) {
    cat("WARN factor_score_heatmap.pdf:", conditionMessage(e), "\n")
    try(dev.off(), silent = TRUE)
})

# 6. Mark the treatment-coupled factors (FDR < 0.10 on treatment, treatment:time, or time)
treat_assoc <- assoc_dt[grep("treatment|time", term)]
sig_factors <- unique(treat_assoc[fdr < 0.10, factor])
cat("Treatment/time-associated factors (FDR<0.10):\n")
print(treat_assoc[fdr < 0.10])
fwrite(data.table(sig_factors = sig_factors),
       file.path(WORK, "sig_factors.tsv"), sep = "\t")

# Copy small final outputs back to /project for the manuscript record.
PROJ_OUT <- file.path(PROJ, "mofa")
dir.create(PROJ_OUT, recursive = TRUE, showWarnings = FALSE)
file.copy(c(file.path(WORK, "factor_variance_explained.tsv"),
            file.path(WORK, "factor_scores.tsv"),
            file.path(WORK, "factor_loadings_top.tsv"),
            file.path(WORK, "factor_treatment_associations.tsv"),
            file.path(WORK, "sig_factors.tsv")),
          PROJ_OUT, overwrite = TRUE)

PROJ_FIG <- file.path(PROJ, "figures", "stage_G")
dir.create(PROJ_FIG, recursive = TRUE, showWarnings = FALSE)
file.copy(list.files(FIG, full.names = TRUE), PROJ_FIG, overwrite = TRUE)

cat("[", format(Sys.time()), "] t08a_stage_G_mofa.R complete.\n", sep = "")
