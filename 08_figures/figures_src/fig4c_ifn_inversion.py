"""Fig 4C - HALLMARK_INTERFERON_ALPHA_RESPONSE sign-flip: suppressed at 24 HPT, induced at 7 DPT.

Two fgsea mountains (BRAIN MDT vs ND), each with the SAME plot-area size as 4B (3.52 x 1.5 in).
Line colour encodes direction (down/blue at 24 HPT -> up/coral at 7 DPT), so the sign-flip
is visible in the palette. Local fgsea reproduction via _mountain.

2026-09-06 (Levin round 2): the two mountains are STACKED (24 HPT over 7 DPT) so the panel
fits one column of the portrait Fig 4 grid at 4.4 in wide, placed 1:1 (XENO_FS_SCALE = 1.15).
Sizes follow the rebuild house style: title DECK_FS_TITLE regular, ticks DECK_FS_AXIS.
"""
import shutil
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from _mountain import ROOT, render_mountain
import _style as S

FIGSIZE = (4.4, 4.6)

OUT = ROOT / "figures_png"; OUT.mkdir(exist_ok=True)
MIR = ROOT / "figure_rebuild_2026-07" / "panels_png" / "fig4C_ifn_inversion.png"

fig = plt.figure(figsize=FIGSIZE)
ax1 = fig.add_axes([0.15, 0.565, 0.80, 0.31])              # each plot area 3.52 x 1.43 in
ax2 = fig.add_axes([0.15, 0.105, 0.80, 0.31])
kw = dict(title_size=S.DECK_FS_AXIS + 0.5, tick_size=S.DECK_FS_AXIS, label_size=S.DECK_FS_AXIS,
          annot_size=S.DECK_FS_LEGEND, title_weight="normal")
render_mountain(ax1, "CT_vs_ND_24HPA", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "24 HPT", **kw)
render_mountain(ax2, "CT_vs_ND_7DPA", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "7 DPT", **kw)
fig.suptitle("Interferon alpha response: brain MDT vs ND", fontsize=S.DECK_FS_TITLE,
             fontweight="normal", y=0.985)
fp = OUT / "fig4c_ifn_inversion.png"
fig.savefig(fp, dpi=S.DEFAULT_DPI, bbox_inches=Bbox.from_bounds(0, 0, *FIGSIZE))
plt.close(fig)
shutil.copyfile(fp, MIR)
print("wrote", fp.name, "->", MIR.name)
