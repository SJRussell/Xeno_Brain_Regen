#!/bin/bash
#SBATCH --job-name=v9_build_suppl
#SBATCH --account=def-srussell
#SBATCH --time=00:30:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2
#SBATCH --output=/scratch/srussell/srussell/logs/v9_build_suppl_%j.out
#SBATCH --error=/scratch/srussell/srussell/logs/v9_build_suppl_%j.err

set -euo pipefail
MANIFEST="/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/${SLURM_JOB_NAME}_${SLURM_JOB_ID}.yaml"
cat > "$MANIFEST" <<EOF
job_id: ${SLURM_JOB_ID}
job_name: ${SLURM_JOB_NAME}
task: v9 supplementary build - DE topTables + full GSEA + Walker per-cell + tar package
script: /scratch/srussell/srussell/scripts/v9/build_supplementary.R
submitted: $(date -Is)
EOF

module load StdEnv/2023 r/4.5.0 r-bundle-bioconductor/3.21
export R_LIBS_USER="${HOME}/Rlibs/4.5:${R_LIBS_USER:-}"
cd /project/def-srussell/nirosha_xeno
Rscript /scratch/srussell/srussell/scripts/v9/build_supplementary.R

# --- Stage all supplementary content into a working dir ---
STAGE=/project/def-srussell/nirosha_xeno/supplementary_v9_staging
mkdir -p "$STAGE/wgcna" "$STAGE/mofa" "$STAGE/liao" "$STAGE/walker_ucell" "$STAGE/hallmark_50_v9redfix"

# WGCNA
cp /project/def-srussell/nirosha_xeno/wgcna/module_assignments_BRAIN.tsv               "$STAGE/wgcna/"
cp /project/def-srussell/nirosha_xeno/wgcna/brain_module_eigengene_trait_lm.tsv        "$STAGE/wgcna/"
cp /project/def-srussell/nirosha_xeno/wgcna/brain_module_gdf15L_correlations.tsv       "$STAGE/wgcna/"
cp /project/def-srussell/nirosha_xeno/wgcna/brain_module_hallmark_enrichment.tsv       "$STAGE/wgcna/"
cp /project/def-srussell/nirosha_xeno/wgcna/hub_genes_BRAIN.tsv                        "$STAGE/wgcna/"
cp /project/def-srussell/nirosha_xeno/wgcna/module_assignments_BLAST.tsv               "$STAGE/wgcna/"
cp /project/def-srussell/nirosha_xeno/wgcna/hub_genes_BLAST.tsv                        "$STAGE/wgcna/"
# eigengenes per sample: extracted from wgcna_BRAIN.rds at session close
Rscript -e '
w <- readRDS("/project/def-srussell/nirosha_xeno/wgcna/wgcna_BRAIN.rds")
MEs <- if (!is.null(w$MEs)) w$MEs else w$net$MEs
write.table(data.frame(sample = rownames(MEs), MEs, check.names = FALSE),
            "'"$STAGE"'/wgcna/brain_module_eigengenes_per_sample.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
'

# MOFA
cp /project/def-srussell/nirosha_xeno/mofa/factor_loadings_top.tsv             "$STAGE/mofa/"
cp /project/def-srussell/nirosha_xeno/mofa/factor_scores.tsv                   "$STAGE/mofa/"
cp /project/def-srussell/nirosha_xeno/mofa/factor_treatment_associations.tsv   "$STAGE/mofa/"
cp /project/def-srussell/nirosha_xeno/mofa/factor_variance_explained.tsv       "$STAGE/mofa/"
cp /project/def-srussell/nirosha_xeno/mofa/factor_vs_brain_gdf15L_spearman.tsv "$STAGE/mofa/"

# Liao
cp /project/def-srussell/nirosha_xeno/comparators/liao/per_celltype_interaction.tsv "$STAGE/liao/"
cp /project/def-srussell/nirosha_xeno/comparators/liao/per_celltype_scores.tsv      "$STAGE/liao/"
cp /project/def-srussell/nirosha_xeno/comparators/liao/per_celltype_markers.tsv     "$STAGE/liao/"

# Walker UCell + Fisher
cp /project/def-srussell/nirosha_xeno/comparators/walker/ucell_signatures_per_cluster_v2.tsv  "$STAGE/walker_ucell/"
cp /project/def-srussell/nirosha_xeno/figures_v9/fig6/fig6b_cluster4_fractions.tsv            "$STAGE/walker_ucell/"
cp /project/def-srussell/nirosha_xeno/figures_v9/fig6/fig6c_cluster_ranking.tsv               "$STAGE/walker_ucell/"
cp /project/def-srussell/nirosha_xeno/figures_v9/fig6/fig6d_cluster4_proportion_fisher.tsv    "$STAGE/walker_ucell/"

# Full 50-Hallmark v9 RED-fix
cp /project/def-srussell/nirosha_xeno/pathway/hallmark_nes_BRAIN_24HPT_full.tsv                       "$STAGE/hallmark_50_v9redfix/"
cp /project/def-srussell/nirosha_xeno/pathway/hallmark_nes_BRAIN_7DPT_full.tsv                        "$STAGE/hallmark_50_v9redfix/"
cp /project/def-srussell/nirosha_xeno/comparators/descalzi/gsea_mPFC_full.tsv                        "$STAGE/hallmark_50_v9redfix/"
cp /project/def-srussell/nirosha_xeno/comparators/descalzi/gsea_NAc_full.tsv                         "$STAGE/hallmark_50_v9redfix/"
cp /project/def-srussell/nirosha_xeno/comparators/descalzi/gsea_PAG_full.tsv                         "$STAGE/hallmark_50_v9redfix/"
cp /project/def-srussell/nirosha_xeno/comparators/erdei/gsea_hallmark_bi_vs_nt_full.tsv              "$STAGE/hallmark_50_v9redfix/"
cp /project/def-srussell/nirosha_xeno/comparators/erdei/gsea_hallmark_fa_vs_nt_full.tsv              "$STAGE/hallmark_50_v9redfix/"
cp /project/def-srussell/nirosha_xeno/comparators/erdei/gsea_hallmark_bi_vs_fa_full.tsv              "$STAGE/hallmark_50_v9redfix/"

# README inventory
cat > "$STAGE/README.md" <<'EOF_README'
# Supplementary data tables for v9 figure set

Generated 2026-05-13 alongside the v9 scorial-rendered main-figure panels at
`figures_v9/fig{1..7}/{a..e}_panel.{pdf,png}`. All X. laevis analyses use the
per-homeolog primary layer (Methods §4.6 manuscript-of-record).

## Contents

### de_topTables/
limma topTables exported per contrast from `voom_BRAIN_primary.rds` and
`voom_BLAST_primary.rds`. Columns: gene_id, logFC, AveExpr, t, P.Value,
adj.P.Val, B, contrast. Per-homeolog (.L / .S resolved separately).

- `de_BLAST_CT_vs_ND_7DPA.tsv`        — Fig 2.a anchor; tgfb2.L logFC=-2.04 q=0.001
- `de_BLAST_interaction_7Dv24H.tsv`   — interaction (CT_late - CT_early) - (ND_late - ND_early)
- `de_BRAIN_CT_vs_ND_24HPA.tsv`       — Fig 3.a anchor; TNFA leading-edge source
- `de_BRAIN_CT_vs_ND_7DPA.tsv`        — Fig 3.c + Fig 4.a anchor
- `de_BRAIN_interaction_7Dv24H.tsv`   — interaction

### gsea_full_pcutoff1/
clusterProfiler GSEA results with `pvalueCutoff = 1` (full untrimmed Hallmark
and GO BP NES vectors). Columns: term, NES, pvalue, padj, qvalue, core_enrichment.
For each of 5 contrasts:

- `hallmark_<contrast>_full.tsv`  — MSigDB Hallmark (50 sets, all reported)
- `gobp_<contrast>_full.tsv`      — GO Biological Process (all sets with size >=5)

### wgcna/
- `module_assignments_BRAIN.tsv`              — gene_id -> module_label / module_color
- `module_assignments_BLAST.tsv`
- `brain_module_eigengenes_per_sample.tsv`    — sample × ME matrix
- `brain_module_eigengene_trait_lm.tsv`       — lm(ME ~ treatment * time) per module
- `brain_module_gdf15L_correlations.tsv`      — Spearman ME vs gdf15.L per module
- `brain_module_hallmark_enrichment.tsv`      — Hallmark over-representation per module
- `hub_genes_BRAIN.tsv`, `hub_genes_BLAST.tsv`

### mofa/
- `factor_loadings_top.tsv`            — top loadings per factor (NOT all genes;
                                          MOFA HDF5 reload required for full loadings)
- `factor_scores.tsv`                  — per-sample factor scores
- `factor_variance_explained.tsv`      — variance per factor per view
- `factor_treatment_associations.tsv`  — lm(factor ~ treatment * time)
- `factor_vs_brain_gdf15L_spearman.tsv` — Spearman vs gdf15.L

### liao/
- `per_celltype_interaction.tsv`  — lm interaction by Liao cell type; multi-celltype
                                     FDR for Fig 4.c (zbtb18+ excitatory FDR_int = 0.081)
- `per_celltype_scores.tsv`       — UCell scores per (sample, cell type)
- `per_celltype_markers.tsv`      — Liao marker panels per cell type

### walker/
- `walker_per_cell_cluster_treatment.tsv`  — 19,827 cells × (cluster, treatment)

### walker_ucell/
- `ucell_signatures_per_cluster_v2.tsv`         — per (cluster, treatment) UCell, multiple panels
- `fig6b_cluster4_fractions.tsv`                 — cluster 4 cell counts/pct per treatment
- `fig6c_cluster_ranking.tsv`                    — cluster ranking by mean UCell
- `fig6d_cluster4_proportion_fisher.tsv`         — Fisher exact OR vs uninjured

### hallmark_50_v9redfix/
Group C RED-fix outputs (2026-05-12): full 50-Hallmark NES vectors with NO
sig-filter. Use these for any cross-species correlation. Earlier biased subsets
(stored in pathway/all_enrichment_summary.tsv with clusterProfiler pvalueCutoff=0.25
and Descalzi/Erdei sig-filtered TSVs) are NOT suitable for cross-comparison.

- `hallmark_nes_BRAIN_24HPT_full.tsv`     — X. laevis voom_BRAIN_primary 24HPA, 50 Hallmarks
- `hallmark_nes_BRAIN_7DPT_full.tsv`      — X. laevis 7DPA
- `gsea_mPFC_full.tsv`                    — Descalzi mPFC fgsea on t-stat, no filter
- `gsea_NAc_full.tsv`                     — Descalzi NAc
- `gsea_PAG_full.tsv`                     — Descalzi PAG
- `gsea_hallmark_bi_vs_nt_full.tsv`       — Erdei burn vs nt fgsea on logFC, no filter
- `gsea_hallmark_fa_vs_nt_full.tsv`       — Erdei formalin vs nt
- `gsea_hallmark_bi_vs_fa_full.tsv`       — Erdei burn vs formalin

## Provenance

- voom EList objects: `/project/def-srussell/nirosha_xeno/deg/voom_{BRAIN,BLAST}_primary.rds`
- GSEA gene-set library: msigdbr Homo sapiens Hallmark + GO BP via org.Xl.eg.db
- WGCNA: signed network, blockwiseModules; module_eigengenes from wgcna_BRAIN.rds
- MOFA: HDF5 model at `mofa/` (not included in tarball; reach for full loadings if needed)
- Walker: walker_integrated_v2.rds (1.6 GB, 4-sample re-integration; cluster 4 = etv1+ rank 1)
- Container: R 4.5.0 + r-bundle-bioconductor 3.21 via Nibi StdEnv/2023
- Seed: 42 for all GSEA / fgsea runs

Builder scripts:
- `scripts/v9/build_supplementary.R` (DE topTables + full pcutoff=1 GSEA + Walker per-cell)
- `scripts/v9/build_supplementary.sh` (wraps Rscript + stages all files + tars)
- `scripts/v9/gen_full_hallmark_nes.R` (Descalzi + Erdei full 50-Hallmark fgsea, 2026-05-12)
- `scripts/v9/gen_xen_full_hallmark.R` (X. laevis full 50-Hallmark clusterProfiler pcutoff=1)
EOF_README

# --- Tar package ---
TARBALL=/project/def-srussell/nirosha_xeno/supplementary_v9.tar.gz
cd "$(dirname "$STAGE")"
tar -czf "$TARBALL" -C "$(dirname "$STAGE")" "$(basename "$STAGE")"
ls -la "$TARBALL"

EXIT=$?
echo "exit_code: ${EXIT}"    >> "$MANIFEST"
echo "finished: $(date -Is)"  >> "$MANIFEST"
exit ${EXIT}
