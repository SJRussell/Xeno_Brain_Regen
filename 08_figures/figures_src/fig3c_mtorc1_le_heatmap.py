"""Fig 3C - BRAIN mTORC1 top leading-edge heatmap.

Canonical render: scripts/t22_round11.R -> brain_MTORC1_topLE_heatmap.png
(tables/_meta/brain_MTORC1_topLE_heatmap_meta.yaml). Pathway HALLMARK_MTORC1_SIGNALING,
NES 1.41 padj 0.048 (CT vs ND, 24 HPT). n_top_shown = 15 leading-edge homeologs ranked by
AveExpr; recorded top_genes hardcoded below (authoritative from the meta). voom z per row,
both axes clustered (ward/euclidean), RdBu_r +-4. Presentation: MDT relabel, concise labels, Arial.
"""
import shutil
import _style as S            # noqa: F401
from _heatmap import ROOT, load_brain, render_heatmap

OUT = ROOT / "figures_png"; OUT.mkdir(exist_ok=True)
MIRROR = ROOT / "figure_rebuild_2026-07" / "panels_png" / "fig3C_mtorc1_heatmap.png"

# 15 recorded top leading-edge genes (brain_MTORC1_topLE_heatmap_meta.yaml : top_genes.genes)
GENES = ["ppia.S", "ppia.L", "canx.L", "ccng1.L", "ldha.S", "tpi1.L", "slc2a1.L",
         "slc2a1.S", "aldoa.L", "acly.L", "slc7a5.L", "acly.S", "eprs1.L",
         "hsp90b1.L", "canx.S"]

E, md = load_brain()
fp = OUT / "fig3c_mtorc1_le_heatmap.png"
# 2026-09-06 (Levin round 2, Stewart: "3C got really small"): rendered at the deck slot size
# (4.7 x 3.9 in, placed 1:1 on slide 3) with a wider heatmap (0.50 of the width) and explicit
# sizes: sample labels 7x, gene labels 6x (15 rows must fit the row pitch), legend 7x,
# title 8.5x the XENO_FS_SCALE (1.35 for Fig 3 -> 9.5 / 8.1 / 9.5 / 11.5 pt).
render_heatmap(E, md, GENES, "BRAIN mTORC1 top leading-edge", fp,
               lock_cols=None, cluster_rows=True, gene_labels=True,
               figsize=(4.7, 3.9), hm_w=0.50,
               title_size=8.5 * S.FS_SCALE, xtick_size=7 * S.FS_SCALE,
               ytick_size=6 * S.FS_SCALE, legend_size=7 * S.FS_SCALE)
shutil.copyfile(fp, MIRROR)
print("mirrored ->", MIRROR.name)
