#!/usr/bin/env Rscript
# Stage F: WGCNA per tissue
# - Soft-threshold β by scale-free fit + connectivity
# - Block-wise modules, min size 30
# - Module-trait correlations (treatment, time, treatment×time interaction)
# - Hub-gene identification via kME
# - Outputs: module assignment, eigengene matrix, topology plots, kME tables

suppressPackageStartupMessages({
    library(WGCNA)
    library(DESeq2)
    library(edgeR)
    library(ggplot2)
    library(matrixStats)
})

# Fix the well-known WGCNA + newer R cor() incompatibility:
# WGCNA's blockwiseModules calls cor(...) with weights.x/weights.y/cosine
# args, which stats::cor() doesn't accept. Aliasing to WGCNA::cor fixes it.
cor <- WGCNA::cor

allowWGCNAThreads(nThreads = 4)
options(stringsAsFactors = FALSE)

PROJ <- "/project/def-srussell/nirosha_xeno"
DEG_DIR <- file.path(PROJ, "deg")
WGCNA_DIR <- file.path(PROJ, "wgcna")
FIG_DIR  <- file.path(PROJ, "figures", "stage_F")
dir.create(WGCNA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("[", format(Sys.time(), "%FT%T%z"), "] Loading filtered+TMM DGEList from Stage C...\n")
y <- readRDS(file.path(DEG_DIR, "dgelist_filtered_tmm.rds"))
meta <- readRDS(file.path(DEG_DIR, "sample_metadata.rds"))

# Use VST-transformed counts as input to WGCNA per recommendation
vsd <- readRDS(file.path(DEG_DIR, "vsd_blind.rds"))
vst_mat <- assay(vsd)

# Per-tissue runs
for (tiss in c("BRAIN", "BLAST")) {
    cat("\n==========================================\n")
    cat("Tissue:", tiss, "\n")
    cat("==========================================\n")

    sub_idx <- meta$tissue == tiss
    expr <- t(vst_mat[, sub_idx])  # WGCNA expects samples in rows
    meta_sub <- meta[sub_idx, ]

    # Filter to top 8000 most variable genes per tissue (standard WGCNA practice
    # for n<20 samples to avoid noise dominating modules)
    vars <- colVars(as.matrix(expr))
    top <- order(vars, decreasing = TRUE)[seq_len(min(8000, ncol(expr)))]
    expr <- expr[, top]
    cat("  WGCNA input:", nrow(expr), "samples × ", ncol(expr), "genes\n")

    # Goodness-of-samples-and-genes check
    gsg <- goodSamplesGenes(expr, verbose = 1)
    if (!gsg$allOK) {
        cat("  Filtering bad genes/samples per goodSamplesGenes...\n")
        expr <- expr[gsg$goodSamples, gsg$goodGenes]
    }
    cat("  Post-QC:", nrow(expr), "samples × ", ncol(expr), "genes\n")

    # Soft-threshold
    cat("\n  Picking soft threshold...\n")
    powers <- c(1:10, seq(12, 30, 2))
    sft <- pickSoftThreshold(expr, powerVector = powers, verbose = 0,
                              networkType = "signed",
                              moreNetworkConcepts = FALSE)
    cat("  Power table:\n")
    print(sft$fitIndices[, 1:6])
    chosen_power <- sft$powerEstimate
    if (is.na(chosen_power) || is.null(chosen_power)) chosen_power <- 14
    cat("  Chosen power:", chosen_power, "\n")

    # Save sft plot
    pdf(file.path(FIG_DIR, paste0("sft_scale_free_", tiss, ".pdf")), width = 10, height = 5)
    par(mfrow = c(1, 2))
    plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
         xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
         type = "n", main = paste(tiss, "scale-free fit"))
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
         labels = powers, cex = 0.8, col = "red")
    abline(h = 0.85, col = "red", lty = 2)
    plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
         xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
         type = "n", main = paste(tiss, "mean connectivity"))
    text(sft$fitIndices[, 1], sft$fitIndices[, 5],
         labels = powers, cex = 0.8, col = "red")
    dev.off()

    # Block-wise modules
    cat("\n  Computing block-wise modules (power=", chosen_power, ", min size 30)...\n")
    bwnet <- blockwiseModules(expr,
        power = chosen_power,
        networkType = "signed",
        TOMType = "signed",
        minModuleSize = 30,
        reassignThreshold = 0,
        mergeCutHeight = 0.25,
        numericLabels = TRUE,
        pamRespectsDendro = FALSE,
        saveTOMs = FALSE,
        verbose = 0,
        maxBlockSize = 8000)

    cat("  Module sizes (top 12):\n")
    print(sort(table(bwnet$colors), decreasing = TRUE)[1:12])

    # Module eigengenes
    MEs <- bwnet$MEs
    moduleColors <- labels2colors(bwnet$colors)

    # Module-trait correlations
    traits <- data.frame(
        treatment_CT = ifelse(meta_sub$treatment == "CT", 1, 0),
        time_late    = ifelse(meta_sub$time == "7DPA", 1, 0),
        time_24h     = ifelse(meta_sub$time == "24HPA", 1, 0)
    )
    if ("11HPA" %in% as.character(meta_sub$time)) {
        traits$time_11h <- ifelse(meta_sub$time == "11HPA", 1, 0)
    }

    moduleTraitCor <- cor(MEs, traits, use = "pairwise.complete.obs")
    moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(expr))

    cat("  Module-trait correlation matrix (head):\n")
    print(head(round(moduleTraitCor, 3)))

    # Hub genes (top 5 kME per module)
    kME <- signedKME(expr, MEs)
    colnames(kME) <- sub("^kME", "", colnames(kME))
    hub_table <- list()
    for (m in unique(bwnet$colors)) {
        mod_genes <- colnames(expr)[bwnet$colors == m]
        if (length(mod_genes) < 3) next
        kme_for_module <- kME[mod_genes, paste0(m), drop = TRUE]
        names(kme_for_module) <- mod_genes
        top5 <- sort(kme_for_module, decreasing = TRUE)[1:min(10, length(kme_for_module))]
        hub_table[[as.character(m)]] <- data.frame(
            module = m,
            color = labels2colors(m),
            gene = names(top5),
            kME = round(top5, 4)
        )
    }
    hub_df <- do.call(rbind, hub_table)

    # Save outputs
    saveRDS(list(bwnet = bwnet, MEs = MEs, moduleColors = moduleColors,
                 moduleTraitCor = moduleTraitCor, moduleTraitP = moduleTraitP,
                 traits = traits, kME = kME, sft = sft, chosen_power = chosen_power),
            file.path(WGCNA_DIR, paste0("wgcna_", tiss, ".rds")))

    # Module assignment table: gene → module
    gene_assign <- data.frame(
        gene_id = colnames(expr),
        module_label = bwnet$colors,
        module_color = moduleColors,
        stringsAsFactors = FALSE
    )
    write.table(gene_assign,
                file.path(WGCNA_DIR, paste0("module_assignments_", tiss, ".tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(hub_df,
                file.path(WGCNA_DIR, paste0("hub_genes_", tiss, ".tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)

    # Module-trait heatmap
    pdf(file.path(FIG_DIR, paste0("module_trait_", tiss, ".pdf")), width = 10, height = 8)
    textMatrix <- paste(round(moduleTraitCor, 2), "\n(p=",
                        signif(moduleTraitP, 2), ")", sep = "")
    dim(textMatrix) <- dim(moduleTraitCor)
    par(mar = c(8, 9, 3, 2))
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = colnames(traits),
                   yLabels = rownames(moduleTraitCor),
                   ySymbols = rownames(moduleTraitCor),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE, cex.text = 0.6,
                   cex.lab = 0.7,
                   zlim = c(-1, 1),
                   main = paste("Module-trait correlations:", tiss))
    dev.off()

    # Tissue-level summary
    n_modules <- length(unique(bwnet$colors))
    cat(sprintf("\n  %s WGCNA: %d modules from %d genes; chosen power=%d\n",
                tiss, n_modules, ncol(expr), chosen_power))
}

cat("\n[", format(Sys.time(), "%FT%T%z"), "] Stage F done.\n")
