"""Fig 4D - BRAIN MYC targets top leading-edge heatmap.

Canonical render: scripts/t22_round11.R -> brain_MYC_TARGETS_topLE_heatmap.png
(tables/_meta/brain_MYC_TARGETS_topLE_heatmap_meta.yaml). Pathway HALLMARK_MYC_TARGETS_V1,
NES 1.79 padj 3.3e-06 (CT vs ND, 7 DPT). n_top_shown = 15 leading-edge homeologs ranked by
AveExpr; recorded top_genes hardcoded below (authoritative from the meta). voom z per row,
both axes clustered (ward/euclidean), RdBu_r +-4. Presentation: MDT relabel, concise labels, Arial.
"""
import shutil
import _style as S            # noqa: F401
from _heatmap import ROOT, load_brain, render_heatmap

OUT = ROOT / "figures_png"; OUT.mkdir(exist_ok=True)
MIRROR = ROOT / "figure_rebuild_2026-07" / "panels_png" / "fig4D_myc_heatmap.png"

# 15 recorded top leading-edge genes (brain_MYC_TARGETS_topLE_heatmap_meta.yaml : top_genes.genes)
GENES = ["ppia.S", "ppia.L", "pabpc1.S", "serbp1.S", "canx.L", "rplp0.L", "rpl6.L",
         "rack1.L", "serbp1.L", "rpl6.S", "ddx21.S", "ddx21.L", "rps6.L",
         "rps10.L", "eif3d.S"]

E, md = load_brain()
fp = OUT / "fig4d_myc_le_heatmap.png"
# 2026-09-06 (Levin round 2): Fig 4 is a portrait 2-column grid with every panel placed 1:1
# at 4.4 in wide (XENO_FS_SCALE = 1.15 for the whole figure). Sizes follow the rebuild house
# style: title DECK_FS_TITLE (regular), sample labels DECK_FS_AXIS, legend DECK_FS_LEGEND;
# gene labels 7.4x so 15 rows fit the row pitch at this height.
render_heatmap(E, md, GENES, "BRAIN MYC targets top leading-edge", fp,
               lock_cols=None, cluster_rows=True, gene_labels=True,
               figsize=(4.4, 4.3), hm_w=0.50,
               title_size=S.DECK_FS_TITLE, xtick_size=S.DECK_FS_AXIS,
               ytick_size=7.4 * S.FS_SCALE, legend_size=S.DECK_FS_LEGEND)
shutil.copyfile(fp, MIRROR)
print("mirrored ->", MIRROR.name)
