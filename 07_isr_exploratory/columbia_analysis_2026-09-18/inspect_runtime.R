cat(R.version.string, '\n')
pkgs <- c('limma','edgeR','DESeq2','fgsea','jsonlite','data.table','readxl','msigdbr','AnnotationDbi','org.Xl.eg.db','BiocManager')
for (p in pkgs) cat(p, if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else 'MISSING', '\n')
cat('Library paths:', .libPaths(), '\n')
