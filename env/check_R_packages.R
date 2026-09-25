cat("Searching for available R packages...\n")
ip <- installed.packages()
cat("Total packages:", nrow(ip), "\n\n")

wanted <- c("clusterProfiler", "msigdbr", "org.Xl.eg.db", "org.Hs.eg.db",
            "ReactomePA", "reactome.db", "fgsea", "enrichplot",
            "GO.db", "GSEABase", "DOSE")

cat("Wanted packages:\n")
for (p in wanted) {
    if (p %in% rownames(ip)) {
        cat(sprintf("  [OK]  %s %s\n", p, ip[p, "Version"]))
    } else {
        cat(sprintf("  [NA]  %s -- NOT INSTALLED\n", p))
    }
}

cat("\nAll Bioc packages with 'enrich' or 'GSEA' in name:\n")
matches <- grep("enrich|gsea|GSEA|cluster|profil|msig|GO|reactome|DOSE",
                rownames(ip), value = TRUE)
for (p in sort(matches)) {
    cat(sprintf("  %s %s\n", p, ip[p, "Version"]))
}
