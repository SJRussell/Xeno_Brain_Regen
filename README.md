# Code for "Limb-localized induction of regeneration reshapes brain and blastema transcriptional states in adult Xenopus laevis"

Stewart J. Russell, Kelsie A. Miller, David L. Kaplan, Michael Levin, Nirosha J. Murugan

Submitted to npj Regenerative Medicine. TODO: add the article DOI on publication.

## What this repository contains

This repository holds the analysis and figure-generation code for a reanalysis of bulk RNA-seq from adult *Xenopus laevis* after hindlimb amputation, comparing multidrug treatment (MDT) with no device (ND), n=3 per group: whole brain at 24 hours and 7 days post-treatment (HPT, DPT) and limb blastema at 11 HPT, 24 HPT and 7 DPT. The cluster scripts (bash SLURM wrappers with R and Python) ran on the Nibi cluster of the Digital Research Alliance of Canada and cover reference building, read trimming, alignment, quantification, filtering, differential expression, gene-set enrichment, WGCNA, MOFA+ and the comparator datasets. The local Python scripts draw the figure panels from tables exported from the cluster. A separate folder holds the exploratory integrated stress response (ISR) scoring. The repository holds code only: no reads, count matrices, result tables or rendered figures. Scripts keep their original file names so they match the job records on the cluster.

## Label conventions

- The code and tables call the treatment `CT`. The manuscript calls the same group MDT. `ND` (no device) is the same in both.
- The code uses the original study's labels in hours or days post-amputation: `11HPA`, `24HPA` and `7DPA` are the collection points the manuscript calls 11 HPT, 24 HPT and 7 DPT.
- Tissues are `BRAIN` and `BLAST` (blastema, the sampled amputation-site tissue).
- Sample IDs follow treatment, replicate, timepoint, tissue and sequencing sample number, for example `CT_T1_24HPA_BRAIN_S21`.
- Gene symbols ending in `.L` or `.S` are the two homeologs of the allotetraploid *X. laevis* genome.

## Folder guide

| Folder | Contents |
|---|---|
| `01_preprocessing_alignment/` | Reference download and indices (`01_reference/`), FastQC and fastp (`02_qc_trim/`), STAR, featureCounts and Salmon (`03_alignment/`), filtering, TMM, VST, QC PCA and outlier check (`04_qc_normalization/`) |
| `02_differential_expression/` | limma-voom (per homeolog and homeolog-summed), DESeq2 Wald tests (per homeolog and homeolog-summed), temporal DEG overlap lists |
| `03_enrichment/` | clusterProfiler GSEA and gseGO, direct fgsea, homeolog-summed GSEA and leading-edge extraction |
| `04_wgcna/` | WGCNA per tissue; module eigengene models, module Hallmark over-representation and eigengene correlation with gdf15.L |
| `05_mofa/` | MOFA+ on the 12 paired brain and blastema samples |
| `06_comparators/` | Axolotl brain single-cell data (`walker_axolotl/`), *X. laevis* brain atlas (`liao_xenopus_atlas/`), mouse spared-nerve injury (`descalzi_sni/`), mouse acute cortical injury (`erdei_acute_injury/`), and `gen_full_hallmark_nes.R`, which writes all 50 Hallmark NES values for the mouse contrasts. No scripts for the GSE73419 uninjured baseline are included (see the TODO list below). |
| `07_isr_exploratory/` | Exploratory ISR scoring (`columbia_analysis_2026-09-18/`) and the script that draws Supplementary Figure 3 (`revision_2026-09-18_ISR/insert_isr.py`) |
| `08_figures/figure_rebuild_2026-07/` | Python panel scripts for Figures 1 to 9 and Supplementary Figures 1 and 2, shared helpers (`_style.py`, `_traj.py`, `_panels.py`, `fig3f_layout.py`), the font-scale driver `rebuild_bigfonts_2026-09-06.py`, and `nibi_export/export_fig_data.R`, which exported DE tables and expression matrices from the cluster |
| `08_figures/figures_src/` | Heatmap and running-enrichment panel scripts with their own `_style.py`, `_heatmap.py` and `_mountain.py` |
| `08_figures/revision_2026-09-17/` | Scripts that rebuild and compose Figures 2, 6, 7 and 9 as submitted |
| `08_figures/r_panel_scripts/` | Cluster R scripts: PCA for Figure 1e and 1f, the axolotl panels and tests in their R form (`v9/`), the script that staged the supplementary and comparator tables (`v9/build_supplementary.*`), and R versions of earlier brain and blastema panels (`t18_*`, `t19_*`) that are not used to draw the submitted figures |
| `env/` | R package checks and the MOFA2, UCell and BayesPrism install script for the cluster |

## Figure panels and the scripts that draw them

Script paths in this table are relative to `08_figures/` unless they start with a numbered stage folder. "Typed in" means the script holds fixed values (term lists, gene lists or test statistics) copied from records written by cluster scripts that are not in this repository.

| Panel | What it shows | Script | Inputs and notes |
|---|---|---|---|
| 1A to 1D | Design, sampling, workflow and comparator schematic | none | Hand-drawn schematic. |
| 1E | Brain PCA | `figure_rebuild_2026-07/render_fig1_pca.py` | PC scores from `r_panel_scripts/fig1D_brain_pca.R` (homeolog-summed voom expression, 5,000 most variable genes). The axis-label percentages are typed into `render_fig1_pca.py`. |
| 1F | Blastema PCA | `figure_rebuild_2026-07/render_fig1_pca.py` | PC scores from `r_panel_scripts/fig2C_blast_pca.R`; as 1E. |
| 2A | Volcano plots, MDT vs ND, 24 HPT and 7 DPT | `figure_rebuild_2026-07/render_fig2.py` (`volcano`) | `de_BRAIN_CT_vs_ND_{24HPA,7DPA}.tsv`, DESeq2 results from `02_differential_expression/t05c_deseq2_wald.R`, exported by `figure_rebuild_2026-07/nibi_export/export_fig_data.R` |
| 2B | Within-condition volcano plots, 7 DPT vs 24 HPT | `figure_rebuild_2026-07/render_fig2.py` (`volcano`) | `de_BRAIN_{CT,ND}_7DPA_vs_24HPA.tsv` |
| 2C | Heatmap of MDT temporal DEGs | `revision_2026-09-17/rebuild_figure2.py`, which calls `figures_src/_heatmap.py` | Row order from `brain_heatmap_CTtemporal_DEGs_gene_order.tsv` (written on the cluster, TODO); column order fixed in `_heatmap.CANON_2F_COLS`. `figures_src/fig2f_cttemporal_heatmap.py` is the stand-alone version of this panel. |
| 2D | Overlap of the two temporal DEG sets | `figure_rebuild_2026-07/render_fig2.py` (`temporal_venn`) | Recomputed from the two temporal DE tables |
| Figure 2 layout | Composition | `revision_2026-09-17/rebuild_figure2.py` (imports `rebuild_figures.py`) | |
| 3A | Hallmark GSEA, 24 HPT | `figure_rebuild_2026-07/render_fig3.py` | `hallmark_nes_BRAIN_24HPT_v10method.tsv` (TODO) |
| 3B | vwf.L, col4a1.L and col4a2.L expression | `figure_rebuild_2026-07/render_fig3b.py` | Per-homeolog voom expression (`voom_E_BRAIN_primary.tsv`) |
| 3C | mTORC1 leading-edge heatmap | `figures_src/fig3c_mtorc1_le_heatmap.py` with `_heatmap.py` | Gene list typed in (TODO) |
| 3D | Neurovascular signature, log2 fold change per timepoint | `figure_rebuild_2026-07/render_fig3d.py` | `nibi_export/neurovascular_signature_table.tsv` (TODO) |
| 3E | Reactome over-representation | `figure_rebuild_2026-07/render_fig3.py` with `_panels.py` | `reactome_overrep_CTvsND_24HPT.tsv` (TODO) |
| 3F | Gene-pathway network | `figure_rebuild_2026-07/render_fig3.py`, which calls `fig3f_layout.py`; `drag_fig3f.py` is an interactive layout editor | `reactome_cnet_edges_CTvsND_24HPT.tsv` (TODO); node positions in `panels_png/fig3F_positions.json` when present |
| 4A | Hallmark GSEA, 7 DPT | `figure_rebuild_2026-07/render_fig4.py` (`panel_4A`) | Terms, NES and adjusted P typed in (TODO) |
| 4B | Heme metabolism running-enrichment plot | `figures_src/fig4b_heme_mountain.py` with `_mountain.py` | Ranks from `de_BRAIN_CT_vs_ND_7DPA.tsv`; gene sets and Entrez map from `tables/brain/gsea_mountain_inputs/` (TODO); NES and FDR from the v10method Hallmark table |
| 4C | IFN-α running-enrichment plots, 24 HPT and 7 DPT | `figures_src/fig4c_ifn_inversion.py` with `_mountain.py` | As 4B |
| 4D | MYC targets leading-edge heatmap | `figures_src/fig4d_myc_le_heatmap.py` with `_heatmap.py` | Gene list typed in (TODO) |
| 4E | Hallmark over-representation in the turquoise module | `figure_rebuild_2026-07/render_fig4.py` (`panel_4E`) | Values typed in (TODO); the test itself is in `04_wgcna/fig5_module_gdf15_rho.R` |
| 4F | Per-Hallmark NES, 24 HPT vs 7 DPT | `figure_rebuild_2026-07/render_fig4.py` (`panel_4F`) | v10method Hallmark tables (TODO) |
| 5A | Brain gdf15.L expression | `figure_rebuild_2026-07/render_fig5.py` (`panel_5A`) | voom expression; DESeq2 adjusted P for the stars; interaction estimate typed in (TODO) |
| 5B | Blastema gdf15.L and gdf15.S expression | `figure_rebuild_2026-07/render_fig5.py` (`panel_5B`) | voom expression; DESeq2 adjusted P for the stars |
| 5C | Turquoise eigengene vs gdf15.L | `figure_rebuild_2026-07/render_fig5.py` (`panel_5C`) | Eigengenes from `04_wgcna/t07_stage_F_wgcna.R`, exported by `r_panel_scripts/v9/build_supplementary.sh`; the fit is recomputed and checked against typed-in values |
| 5D | Paired blastema vs brain gdf15.L (VST) | `figure_rebuild_2026-07/render_fig5.py` (`panel_5D`) | `persample_gdf15_{BRAIN,BLAST}.tsv` (TODO) |
| 6A | Synaptic and axon GO:BP enrichment, 24 HPT and 7 DPT | `figure_rebuild_2026-07/render_fig6.py` (`panel_6A_expanded`), run by `revision_2026-09-17/rebuild_figures.py` | `go_bp_rds/goBP_FULL_BRAIN_CT_vs_ND_{24HPA,7DPA}.tsv` (TODO) |
| 6B | Synaptic and structural gene expression | `figure_rebuild_2026-07/render_fig6.py` (`panel_6B`) | voom expression; interaction estimates typed in (TODO) |
| 6C | zbtb18+ excitatory-neuron UCell score | `figure_rebuild_2026-07/render_fig6.py` (`panel_6C`) | `per_celltype_{scores,interaction}.tsv` from `06_comparators/liao_xenopus_atlas/hD_liao_analyze.R` |
| Figure 6 layout | Composition | `revision_2026-09-17/rebuild_figures.py` | |
| 7A | Process-state tiles, 24 HPT and 7 DPT | `figure_rebuild_2026-07/render_fig7.py` (`panel_7A`), run by `rebuild_figures.py` | `go_process_clusters_BRAIN.tsv`, `gobp_BRAIN_CTvND_{24HPT,7DPT}.tsv` (TODO) |
| 7B | Hallmark shift from the uninjured baseline | `revision_2026-09-17/rebuild_figures.py` (`grouped7b`, using `render_fig7._hallmark_vs_baseline`) | `gsea_rank_shift_{CT,ND}.tsv` (TODO) |
| 7C | Concordance of the two conditions' shifts | `figure_rebuild_2026-07/render_fig7.py` (`panel_7C`) | `rank_shift_genome_wide.tsv` (TODO) |
| 7D | gdf15.L and gdf15.S percentile vs baseline | `figure_rebuild_2026-07/render_fig7.py` (`panel_7D`) | `rank_shift_genome_wide.tsv`, `tpm_brain_combined.tsv` (TODO) |
| 7E | Axolotl etv1+ signature score vs baseline | `figure_rebuild_2026-07/render_fig7.py` (`panel_7E`) | `tpm_brain_combined.tsv` (TODO) |
| Figure 7 layout | Composition | `revision_2026-09-17/rebuild_figures.py` | |
| 8A | MOFA variance explained and factor-covariate associations | `figure_rebuild_2026-07/render_fig8.py` (`panel_8A`) | `factor_variance_explained.tsv`, `factor_treatment_associations.tsv` from `05_mofa/t08a_stage_G_mofa.R` |
| 8B | Factor 3 scores | `figure_rebuild_2026-07/render_fig8.py` (`panel_8B`) | `factor_scores.tsv` |
| 8C, 8D | Top Factor 3 loadings, brain and blastema | `figure_rebuild_2026-07/render_fig8.py` (`_loadings`) | `factor_loadings_top.tsv` |
| 9A | Axolotl brain UMAP | `figure_rebuild_2026-07/render_fig9.py` (`panel_9A`), run by `rebuild_figures.py` | `walker_umap_coords_v2.tsv` (TODO); `r_panel_scripts/v9/fig6a_walker_umap_clusters.R` draws the same UMAP in R from the Seurat object |
| 9B | etv1+ signature score by cluster | `figure_rebuild_2026-07/render_fig9.py` (`panel_9B`) | `ucell_signatures_per_cluster_v2.tsv` from `06_comparators/walker_axolotl/walker_marker_remap.R`; R version in `r_panel_scripts/v9/fig6c_walker_ucell_ranked_bar.R` |
| 9C | etv1+ signature in *X. laevis* brain | `figure_rebuild_2026-07/render_fig9.py` (`panel_9C`) | UCell recomputed in Python from voom expression; R version in `r_panel_scripts/v9/fig6e_xlaevis_walker_etv1_ucell.R` |
| 9D | Hallmark NES across comparator contrasts | `revision_2026-09-17/rebuild_figures.py` (`grouped9d`, using `render_fig9.COMP`) | *X. laevis* columns from the v10method Hallmark tables (TODO); mouse columns from `supplementary_v9/supplementary_v9_staging/hallmark_50_v9redfix/*_full.tsv`, written by `06_comparators/gen_full_hallmark_nes.R` |
| Figure 9 layout | Composition | `revision_2026-09-17/rebuild_figures.py` | |
| S1A | Blastema Hallmark GSEA, 7 DPT | `figure_rebuild_2026-07/render_supp1.py` (`panel_A`) | Terms typed in (TODO) |
| S1B | Blastema volcano plot with EMT and anti-fibrotic genes | `figure_rebuild_2026-07/render_supp1.py` (`panel_B`) | `de_BLAST_CT_vs_ND_7DPA.tsv`; `blast_LE_genes_per_hallmark_CT_vs_ND_7DPT.tsv` (TODO) |
| S1C | Expression of EMT-intersected DEGs | `figure_rebuild_2026-07/render_supp1.py` (`panel_C`) | Gene list typed in (TODO); voom expression |
| S2A | EMT leading-edge gene expression | `figure_rebuild_2026-07/render_supp1.py` (`panel_traj`, saved as `figS1D_*`) | voom expression; DESeq2 adjusted P for the stars |
| S2B | Fibrosis and EMT-associated gene expression | `figure_rebuild_2026-07/render_supp1.py` (`panel_traj`, saved as `figS1F_*`) | As S2A |
| S2C | EMT running-enrichment plot | `figure_rebuild_2026-07/render_supp1.py` (`panel_E`, saved as `figS1E_*`), using `figures_src/_mountain.py` | Blastema ranks from `de_BLAST_CT_vs_ND_7DPA.tsv` |
| S3 | Limb ISR score by condition and timepoint | `07_isr_exploratory/revision_2026-09-18_ISR/insert_isr.py`, lines 35 to 99 | `results/sample_scores.tsv` and `results/signature_tests.tsv` from `run_analyses.R` |
| Layout of Figures 1, 3, 4, 5, 8 and Supplementary Figures 1 and 2 | Composition | none | Panels were arranged by hand in PowerPoint. |

## Methods analyses and their scripts

| Methods section | Scripts | Not in this repository |
|---|---|---|
| Animals and multidrug treatment | none (original study) | |
| Tissue collection and RNA sequencing | `01_preprocessing_alignment/02_qc_trim/t02a_unzip_inventory.sh` unzips and inventories the delivered FASTQ files | |
| Read processing and quantification | `01_reference/`: `t01a_download_xenla_ref.sh`; `t01b2_rebuild_star_index.sh` (rebuild of the first STAR index from `t01b_build_star_index.sh`); `t01c_fix_and_build_salmon.sh` (fixed version of `t01c_build_salmon_index.sh`). `02_qc_trim/`: `t02b_fastqc_raw.sh`, `t02c_fastp.sh`, `t02d_fastqc_trim.sh` or `t02d2_fastqc_trim.sh`. `03_alignment/`: `t03a_star_align.sh` (two-pass, GeneCounts), `aggregate_star_counts.py` (reverse-strand column into a gene-by-sample matrix), `t03b_featurecounts.sh` (featureCounts count table), `t03c_salmon_quant.sh` | Joint Salmon quantification with GSE73419 (TODO) |
| Filtering and normalization | `04_qc_normalization/t04_stage_C_qc.R` (CPM filter, TMM, VST, QC PCA on the 2,000 most variable genes, outlier check); `02_differential_expression/t05_stage_D_de.R` (limma-voom with quality weights, per homeolog); `t10_homeolog_sum_de.R` (homeolog-summed limma-voom); `08_figures/r_panel_scripts/fig1D_brain_pca.R` and `fig2C_blast_pca.R` (PCA plotted in Figure 1e and 1f: homeolog-summed voom expression, 5,000 most variable genes) | The script that wrote `voom_*_homeolog_sum_strand_corrected.rds`, the input to both PCA scripts (TODO; `t10_homeolog_sum_de.R` writes `voom_*_homeolog_sum.rds`) |
| Differential expression | `02_differential_expression/t05c_deseq2_wald.R` (per homeolog), `t05d_deseq2_wald_homeolog.R` (homeolog-summed), `build_temporal_venn.py` (temporal DEG partitions); `08_figures/figure_rebuild_2026-07/nibi_export/export_fig_data.R` (exports DE tables); `render_fig2.py` (overlap); `render_fig5.py` (paired Spearman correlations: pooled, within MDT and 7 DPT only); `07_isr_exploratory/columbia_analysis_2026-09-18/correlations_and_audit.py` (further cross-tissue checks) | Within-signature limma test for the neurovascular panel (TODO); homeolog-summed coupling check (TODO) |
| Ortholog mapping, gene-set enrichment and over-representation | `03_enrichment/t06_stage_E_enrichment.R` (clusterProfiler GSEA on Hallmark, gseGO on GO:BP); `t06b_stage_E_fgsea.R` (direct fgsea); `t13_homeolog_sum_gsea.R` and `extract_leading_edge_v7.R` (homeolog-summed GSEA and leading-edge lists for an earlier figure version); `08_figures/figures_src/_mountain.py` (running-enrichment curves) | Re-runs of the Hallmark GSEA and gseGO keeping all tested sets (TODO); Reactome over-representation and gene-concept network (TODO); leading-edge tables used by the panels (TODO) |
| Co-expression networks (WGCNA) | `04_wgcna/t07_stage_F_wgcna.R`; `04_wgcna/fig5_module_gdf15_rho.R` (eigengene linear models, Hallmark over-representation of module genes, eigengene correlation with gdf15.L) | |
| Multi-omic factor analysis (MOFA+) | `05_mofa/t08a_stage_G_mofa.R`; `env/setup_install_stage_G_libs.sh` installs MOFA2 | |
| Brain process-state timeline | `08_figures/figure_rebuild_2026-07/render_fig7.py` (`panel_7A`: mean member NES per cluster and the 25% flag) | GO:BP fgsea across 3,473 terms and GOSemSim clustering (TODO) |
| Cell-type signature projection | `06_comparators/liao_xenopus_atlas/hD_liao_download.sh`, `hD_liao_analyze.R` | |
| Axolotl single-cell comparator | `06_comparators/walker_axolotl/`: `hA_walker_download.sh`, `hA_walker_refbuild.sh`, `hA_walker_alevin.sh`; muscimol sample SRR31800874: `hA_walker_dl_muscimol.sh`, `walker_muscimol_recover.sh`, `walker_muscimol_finish.sh`, `walker_alevin_muscimol_only.sh`; `hA_walker_seurat.R` with `walker_seurat_4sample_recluster.sh` (four-sample integration); `walker_marker_remap.R` (marker IDs, UCell scores, Fisher tests). `08_figures/r_panel_scripts/v9/`: `fig6d_walker_cluster4_proportion.R` (cluster proportions and Fisher tests), `fig6e_xlaevis_walker_etv1_ucell.R` (signature on the *X. laevis* samples) | UMAP coordinate export for Figure 9a (TODO) |
| Mammalian injury comparators | `06_comparators/descalzi_sni/hB_descalzi.R`; `06_comparators/erdei_acute_injury/hC_erdei.R`; `06_comparators/gen_full_hallmark_nes.R` (all 50 Hallmark NES per mouse contrast); `08_figures/r_panel_scripts/v9/build_supplementary.sh` (copies those tables into `supplementary_v9_staging/`); `08_figures/revision_2026-09-17/rebuild_figures.py` (Ward clustering and heatmap) | |
| Uninjured baseline anchor | `08_figures/figure_rebuild_2026-07/render_fig7.py` (`panel_7C`, `panel_7D`, `panel_7E`); `08_figures/revision_2026-09-17/rebuild_figures.py` (`grouped7b`) | All GSE73419 processing (TODO) |
| Statistics and reproducibility | Random seeds set in the code: 42 in `t06_stage_E_enrichment.R`, `t06b_stage_E_fgsea.R`, `t13_homeolog_sum_gsea.R`, `gen_full_hallmark_nes.R`, `build_supplementary.R` and `t19_round8.R`; 20260430 for MOFA+ training (`t08a_stage_G_mofa.R`); 20260501 in `hA_walker_seurat.R`; 20260502 in `walker_marker_remap.R`; 20260918 in the ISR `run_analyses.R` | sessionInfo() record for the main R pipeline (TODO) |
| Exploratory integrated stress response scoring | `07_isr_exploratory/columbia_analysis_2026-09-18/`: `get_public_sources.py` and `inspect_sources.py` (download the published ISR sources), `prepare_inputs.py` (Xenbase ortholog mapping, reconstruction of the published scores, adapted coefficients), `run_analyses.R` (scores, ordinary linear models, BH adjustment, sensitivity variants), `correlations_and_audit.py`, `verify_and_plot.py` (independent recomputation of every estimate); `revision_2026-09-18_ISR/insert_isr.py` lines 35 to 99 (Supplementary Figure 3) | Preparation of `isr_training_120.csv`, `published_scores.csv` and the Xenbase ortholog table in `sources/` (TODO) |

## Analyses whose scripts are not in this repository

These scripts ran on the cluster under `/project/def-srussell/nirosha_xeno/` and have not been added here. TODO: add each script, or state where it is available.

1. The GSE73419 uninjured-baseline analysis: trimming and Salmon quantification of GSE73419, joint Salmon quantification of this study's brain libraries, TPM filtering, percentile rank shifts, fgsea on the shift vectors and the eye-sample check (Figure 7b to 7e; tables `gsea_rank_shift_{CT,ND}.tsv`, `rank_shift_genome_wide.tsv`, `tpm_brain_combined.tsv`).
2. The neurovascular signature table with its within-signature limma test on homeolog-collapsed genes (Figure 3d).
3. Reactome over-representation and gene-concept network tables (Figure 3e and 3f).
4. The re-runs of the `t06_stage_E_enrichment.R` Hallmark GSEA and gseGO that keep all tested gene sets: `hallmark_nes_BRAIN_{24HPT,7DPT}_v10method.tsv` (Figures 3a, 4b, 4c, 4f and 9d) and `goBP_FULL_BRAIN_CT_vs_ND_{24HPA,7DPA}.tsv` (Figure 6a).
5. The GO:BP fgsea across 3,473 terms and the GOSemSim semantic-similarity clustering behind `gobp_BRAIN_CTvND_*.tsv` and `go_process_clusters_BRAIN.tsv` (Figure 7a).
6. The cluster panel scripts `t20_round9.R`, `t21_round10.R`, `t22_round11.R`, `t22_round12.R`, `t22_round12_fix.R` and `t23_brain_heatmap_titlefix.R`, which recorded the values typed into the local scripts (Figures 2c row order, 3c, 4a, 4d, 4e, 5a, 5c, 6b, Supplementary Figures 1a and 1c, 2a and 2c) and the leading-edge gene tables.
7. Exports of the per-homeolog voom expression matrices (`voom_E_*_primary.tsv`), the per-sample gdf15 tables (`persample_gdf15_*.tsv`) and the running-enrichment inputs (`gsea_mountain_inputs/`). `export_fig_data.R` as saved exports only the DE tables and the homeolog-summed matrices.
8. The script that wrote `voom_*_homeolog_sum_strand_corrected.rds`, the input to the PCA scripts and to `export_fig_data.R`.
9. The export of UMAP coordinates from the integrated axolotl Seurat object (`walker_umap_coords_v2.tsv`, Figure 9a).
10. The homeolog-summed gdf15 coupling check named in the Methods.
11. Preparation of the ISR reference files in `07_isr_exploratory/columbia_analysis_2026-09-18/sources/`.

## Software

- Stated in the Methods: R 4.5.0 with Bioconductor 3.21; msigdbr 26.1.0; fastp 0.24.0; STAR 2.7.11b; Salmon 1.10.2; DESeq2, limma-voom, edgeR, clusterProfiler, fgsea, WGCNA, MOFA+ (R package MOFA2), UCell, Seurat and GOSemSim (TODO: versions). The ISR analysis used R 4.6.1 with limma 3.68.5; its full package list is in `07_isr_exploratory/columbia_analysis_2026-09-18/R_sessionInfo.txt`. Figure panels use Python 3 with matplotlib (TODO: Python and package versions).
- Also loaded by the SLURM wrappers: the Nibi module environment `StdEnv/2023` with `r/4.5.0` and `r-bundle-bioconductor/3.21`, FastQC 0.12.1, Subread 2.0.6 (featureCounts), SAMtools 1.20 and SRA Toolkit 3.0.9.
- R packages loaded by the scripts: AnnotationDbi, biomaRt, clusterProfiler, data.table, DESeq2, dplyr, edgeR, fgsea, future, GEOquery, ggplot2, ggrepel, ggsci, GO.db, GSEABase, jsonlite, limma, Matrix, matrixStats, MOFA2, msigdbr, openxlsx, org.Hs.eg.db, org.Xl.eg.db, patchwork, pheatmap, RColorBrewer, readr, scales, scran, Seurat, SingleCellExperiment, SummarizedExperiment, tibble, tidyr, tximport, UCell, WGCNA, yaml. No script in this repository loads GOSemSim (see item 5 above).
- Python packages imported by the scripts: numpy, pandas, scipy, matplotlib, matplotlib-venn, adjustText, networkx, Pillow, python-pptx, lxml, openpyxl.
- Environment files: `env/check_R.sh` and `env/check_R_packages.R` list the installed enrichment packages on Nibi; `env/setup_install_stage_G_libs.sh` installs MOFA2, UCell, BayesPrism and the mofapy2 backend into a user library; `07_isr_exploratory/columbia_analysis_2026-09-18/install_analysis_packages.R` and `inspect_runtime.R` set up and check the ISR library.
- TODO: a sessionInfo() and package-version record for the main R pipeline.

## Data inputs

- This study: raw reads and processed count matrices will be deposited at GEO/SRA. TODO: add the accession.
- Reference: *X. laevis* v10.1 assembly (XENLA_10.1) with the NCBI RefSeq annotation from Xenbase, downloaded by `t01a_download_xenla_ref.sh`.
- Comparator datasets: Walker et al. 2025 axolotl brain single-cell RNA-seq (SRR31800874 to SRR31800877); Liao et al. 2022 *X. laevis* brain atlas (GSE195790); Descalzi et al. 2017 mouse spared-nerve injury (GSE91396); Erdei et al. 2025 mouse acute cortical injury (PRJNA1137110); uninjured adult *X. laevis* brain baseline (GSE73419).
- Gene sets: MSigDB collections through msigdbr in the cluster scripts. The ISR scripts download the MSigDB 2024.1 Hallmark GMT, the published ISR code and the supplementary tables of Smith et al.
- The intermediate tables that the figure and ISR scripts read (DE tables, expression matrices, enrichment and comparator tables) are not in this repository. TODO: state where they will be available.

## Running the scripts

### Cluster scripts

- Each R script has a SLURM wrapper (`.sh`) written for Nibi. The wrappers use the account `def-srussell` and run the R scripts from a copy under `/scratch/srussell/srussell/scripts/`.
- The scripts use absolute paths: `/project/def-srussell/nirosha_xeno/` for project data and outputs, `/scratch/srussell/srussell/` for working files and logs, and `/home/srussell/Rlibs/4.5` for the user R library. Re-path these before running anywhere else.
- Many wrappers also write a short job record to `/project/def-srussell/nirosha_xeno/scorial/01_state/slurm_manifests/`. These lines only record the job and can be removed.
- Run order: `01_preprocessing_alignment`, then `02_differential_expression`, `03_enrichment`, `04_wgcna` and `05_mofa`. The comparator scripts in `06_comparators` process their datasets independently; `gen_full_hallmark_nes.R` also reads the *X. laevis* enrichment results from `03_enrichment`. The scripts in `08_figures/r_panel_scripts` read the outputs of the earlier stages.

### Local figure scripts

- Run each script from its own folder, for example `python render_fig5.py` inside `08_figures/figure_rebuild_2026-07/`.
- The scripts find their inputs relative to `08_figures/`. Before running, place these data folders there: `tables/` (with `brain/`, `blastema/`, `wgcna/`, `mofa/`, `comparators/{liao,walker,gse73419}/` and `sample_metadata.tsv`), `supplementary_v9/supplementary_v9_staging/hallmark_50_v9redfix/`, `manuscript_lock_2026-06-08/tables/{brain,blastema}/` (PCA scores), and `figure_rebuild_2026-07/nibi_export/`, which holds `fig_export/` and `neurovascular_signature_table.tsv`.
- `figures_src/_style.py` and `figure_rebuild_2026-07/_style.py` are different modules with the same name. Keep each in its folder and run scripts from the two folders as separate processes.
- `figure_rebuild_2026-07/rebuild_bigfonts_2026-09-06.py render` runs each panel script in its own process with the font scale used for the submitted panels (environment variable `XENO_FS_SCALE`). Running a panel script directly uses scale 1.
- Figures 2, 6, 7 and 9 as submitted come from `revision_2026-09-17/rebuild_figure2.py` (Figure 2) and `revision_2026-09-17/rebuild_figures.py` (Figures 6, 7 and 9). These load the `render_fig*.py` modules from `figure_rebuild_2026-07/`, change some label text and sizes, and write `panels/` and `figures/` inside `revision_2026-09-17/`. They load Arial from `C:/Windows/Fonts`.
- `figures_src/render_all.py` runs every `fig*.py` in that folder.
- Panel images are written to `figure_rebuild_2026-07/panels_png/` and `08_figures/figures_png/`.
- `02_differential_expression/build_temporal_venn.py` has a hard-coded Windows project path and reads temporal DEG lists extracted on the cluster; set `ROOT` before running.

### ISR scripts

- Keep the folder names `columbia_analysis_2026-09-18` and `revision_2026-09-18_ISR`: `run_analyses.R`, `install_analysis_packages.R` and `insert_isr.py` refer to them by name.
- Run from `07_isr_exploratory/`, which the scripts treat as the project root. They also read `tables/sample_metadata.tsv`, `tables/{brain,blastema}/voom_E_*_primary.tsv`, `tables/brain/brain_LE_genes_per_hallmark_*.tsv`, `supplementary_v9/supplementary_v9_staging/hallmark_50_v9redfix/` and `figure_rebuild_2026-07/nibi_export/fig_export/de_*.tsv` relative to that folder.
- Order: `get_public_sources.py`, `inspect_sources.py`, preparation of the remaining `sources/` files (TODO), `prepare_inputs.py`, `run_analyses.R`, `correlations_and_audit.py`, `verify_and_plot.py`, then `insert_isr.py`.
- `inputs/panels.json` holds the marker-gene panels that `prepare_inputs.py` writes and `run_analyses.R` reads.
- `verify_and_plot.py` checks a checksum of a manuscript draft that is not in this repository and stops at that step without it.
- `insert_isr.py`: lines 35 to 99 check the primary contrasts and draw Supplementary Figure 3 into `revision_2026-09-18_ISR/ISR_assets/`. The code from line 100 onward inserts the result text into a Word draft that is not in this repository, so the script stops with a file-not-found error after the figure is written.
