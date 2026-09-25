#!/usr/bin/env Rscript
# One-shot export for the local figure rebuild: per-gene DE for every brain+blast
# contrast, the voom E matrices (homeolog-sum strand-corrected), sample metadata.
# Writes to scratch and tars for a single SFTP pull.
# Unlocks Fig 2C/2D/2F, 3B/3C/3D, 4D, 6B, SuppF1 EMT panels.
suppressPackageStartupMessages({ library(limma) })
DEG <- "/project/def-srussell/nirosha_xeno/deg"
OUT <- "/scratch/srussell/srussell/fig_export"
unlink(OUT, recursive = TRUE); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

KEEP <- c("gene", "gene_id", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "baseMean")

dump_de <- function(rds, tissue) {
  x <- readRDS(file.path(DEG, rds))
  stopifnot(is.list(x), "results" %in% names(x))
  for (nm in names(x$results)) {
    d <- as.data.frame(x$results[[nm]]); d$gene <- rownames(d)
    d <- d[, intersect(KEEP, colnames(d))]
    write.table(d, file.path(OUT, sprintf("de_%s_%s.tsv", tissue, nm)),
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat("de", tissue, nm, nrow(d), "genes\n")
  }
}

dump_voom <- function(rds, tissue) {
  v <- readRDS(file.path(DEG, rds))
  E <- if (!is.null(v$v)) v$v$E else if (!is.null(v$E)) v$E else stop("no E in ", rds)
  E <- as.data.frame(as.matrix(E), check.names = FALSE)
  write.table(cbind(gene = rownames(E), E), file.path(OUT, sprintf("voom_E_%s.tsv", tissue)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("voom_E", tissue, nrow(E), "x", ncol(E), "\n")
}

dump_de("deseq2_BRAIN_primary.rds", "BRAIN")
dump_de("deseq2_BLAST_primary.rds", "BLAST")
dump_voom("voom_BRAIN_homeolog_sum_strand_corrected.rds", "BRAIN")
dump_voom("voom_BLAST_homeolog_sum_strand_corrected.rds", "BLAST")

md <- readRDS(file.path(DEG, "sample_metadata.rds"))
write.table(as.data.frame(md), file.path(OUT, "sample_metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("metadata", nrow(as.data.frame(md)), "rows\n")

system(sprintf("tar czf %s.tar.gz -C %s %s", OUT, dirname(OUT), basename(OUT)))
cat("DONE_EXPORT tarball:", paste0(OUT, ".tar.gz"), "\n")
print(list.files(OUT))
