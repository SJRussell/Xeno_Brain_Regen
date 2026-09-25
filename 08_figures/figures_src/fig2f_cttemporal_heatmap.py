"""Fig 2C (was 2F) - Brain MDT-temporal DEG heatmap: 932 DEGs detected, 924 displayed.

LOCKED to the canonical render (scripts/t23_brain_heatmap_titlefix.R ->
brain_heatmap_CTtemporal_DEGs_sampleclust.png):
  rows    = tables/brain/brain_heatmap_CTtemporal_DEGs_gene_order.tsv  (R ward.D2 leaf order)
  columns = CANON_2F_COLS (7 DPT block | 24 HPT block, timepoint-first)
  z per row across 12 brain samples; RdBu_r clipped +-4; origin='upper' (rank1 top).
Only the presentation changes vs the canonical: MDT relabel, concise sample labels, Arial.
See tables/_meta/brain_heatmap_CTtemporal_DEGs_sampleclust_meta.yaml.
"""
import shutil
from pathlib import Path
import pandas as pd
import _style as S            # noqa: F401  (applies rcParams on import)
from _heatmap import ROOT, load_brain, render_heatmap, CANON_2F_COLS

OUT = ROOT / "figures_png"; OUT.mkdir(exist_ok=True)
MIRROR = ROOT / "figure_rebuild_2026-07" / "panels_png" / "fig2F_temporal_heatmap.png"

E, md = load_brain()
go = pd.read_csv(ROOT / "tables" / "brain" / "brain_heatmap_CTtemporal_DEGs_gene_order.tsv", sep="\t")

fp = OUT / "fig2f_cttemporal_heatmap.png"
# Title reports what is DISPLAYED (924), not the DEG count (932). Verified
# 2026-07-16: de_BRAIN_CT_7DPA_vs_24HPA padj<0.05 = 932 DEGs (573 up / 359 down,
# == Fig 2B and the 2D Venn); 924 of them carry voom expression values and are all
# present here; the other 8 (znf577.S, acer1.L, arl11.L, krt18.3.S + 4 LOC, all
# padj 0.010-0.043) are absent from voom_E_BRAIN_primary entirely -- they cleared
# DESeq2's independent filter but fell below the voom expression filter, so they
# cannot be z-scored. Every displayed row is a DEG (0 non-DEGs in the set).
# Legend must state: 932 detected, 924 displayed.
render_heatmap(E, md, go["gene_id"].tolist(),
               "BRAIN MDT-Temporal DEGs (n = 924 displayed)", fp,
               lock_cols=CANON_2F_COLS, cluster_rows=False, gene_labels=False,
               figsize=(5.2, 6.0))
shutil.copyfile(fp, MIRROR)
print("mirrored ->", MIRROR.name)
