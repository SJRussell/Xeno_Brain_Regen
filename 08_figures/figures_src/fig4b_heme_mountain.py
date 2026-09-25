"""Fig 4B - HALLMARK_HEME_METABOLISM fgsea running-enrichment mountain (BRAIN MDT vs ND, 7 DPT).

Local fgsea reproduction via _mountain (validated identical to the Nibi rk). Fixed plot-area
size (3.52 x 1.51 in) shared with the two 4C IFN mountains so all three are the same size.

2026-09-06 (Levin round 2): Fig 4 became a portrait 2-column grid with every panel placed 1:1
at 4.4 in wide (XENO_FS_SCALE = 1.15). Sizes now follow the rebuild house style explicitly:
title DECK_FS_TITLE regular weight, ticks/axis DECK_FS_AXIS, annotation DECK_FS_LEGEND.
"""
import shutil
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from _mountain import ROOT, render_mountain
import _style as S

FIGSIZE = (4.4, 2.6)

OUT = ROOT / "figures_png"; OUT.mkdir(exist_ok=True)
MIR = ROOT / "figure_rebuild_2026-07" / "panels_png" / "fig4B_heme_mountain.png"

fig = plt.figure(figsize=FIGSIZE)
ax = fig.add_axes([0.15, 0.25, 0.80, 0.58])                # plot area 3.52 x 1.51 in
render_mountain(ax, "CT_vs_ND_7DPA", "HALLMARK_HEME_METABOLISM",
                "Heme metabolism: brain MDT vs ND, 7 DPT",
                title_size=S.DECK_FS_TITLE, tick_size=S.DECK_FS_AXIS,
                label_size=S.DECK_FS_AXIS, annot_size=S.DECK_FS_LEGEND,
                title_weight="normal")
fp = OUT / "fig4b_heme_mountain.png"
fig.savefig(fp, dpi=S.DEFAULT_DPI, bbox_inches=Bbox.from_bounds(0, 0, *FIGSIZE))
plt.close(fig)
shutil.copyfile(fp, MIR)
print("wrote", fp.name, "->", MIR.name)
