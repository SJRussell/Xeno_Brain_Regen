"""Fig 3D — neurovascular gene signature (brain, MDT vs ND), compartment-grouped
log2FC tile map with significance stars. Rebuilt EXACTLY from the curated Nibi
table (nibi_export/neurovascular_signature_table.tsv): logFC_clip = colour, sig =
star. Significance is the within-signature limma FDR from that table, NOT recomputed.

Stewart's edits: drop the Pericyte compartment, CT -> MDT in title/legend, no subtitle.
Gene order within each compartment matches the approved v13 panel (read off the PNG).

PERICYTE, RESOLVED (Stewart 2026-07-17): stays OUT of the panel, and comes OUT of the
Methods and the legend, which currently promise "endothelial + pericyte + tight-junction
+ basement-membrane". Supported by the source table: all 10 pericyte rows (rgs5, cspg4,
mcam, kcnj8, pdgfrb x 2 timepoints) are ns at adj.P.Val 0.688-1.000. Nothing to show.
"""
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Rectangle
import _style as S

HERE = Path(__file__).resolve().parent
TBL  = HERE / "nibi_export" / "neurovascular_signature_table.tsv"
OUT  = HERE / "panels_png"; OUT.mkdir(exist_ok=True)

# compartment + gene order from the approved panel (Pericyte dropped)
ORDER = {
    "Endothelial":       ["vwf", "tek", "pecam1", "kdr", "flt1", "esam", "cdh5"],
    "Tight junction":    ["tjp1", "ocln", "f11r", "cldn5"],
    "Basement membrane": ["nid1", "lamc1", "lamb1", "lama4", "hspg2", "col4a2", "col4a1"],
}
COMPARTS = ["Endothelial", "Tight junction", "Basement membrane"]
COLS = ["24 HPT", "7 DPT"]
ZMAX = 2.0

df = pd.read_csv(TBL, sep="\t")
df["gene_root"] = df["gene_root"].astype(str)

def cell(gene, contrast, field):
    r = df[(df["gene_root"] == gene) & (df["contrast"] == contrast)]
    return r[field].iloc[0] if len(r) else np.nan

cmap = mpl.colormaps[S.CMAP_DIVERGING]
norm = mpl.colors.Normalize(vmin=-ZMAX, vmax=ZMAX)

fig = plt.figure(figsize=(4.0, 5.4))
TILES_L, TILES_W = 0.42, 0.34          # tile block x-span (2 columns)
TV, BV = 0.90, 0.10                    # vertical plot span (room for title / col labels)
n_total = sum(len(ORDER[c]) for c in COMPARTS)
GAP_ROWS = 0.75
u = (TV - BV) / (n_total + GAP_ROWS * (len(COMPARTS) - 1))   # figure-fraction per gene row

y = TV
axes_bottom = None
for ci, comp in enumerate(COMPARTS):
    genes = ORDER[comp]
    n = len(genes)
    h = n * u
    yb = y - h
    ax = fig.add_axes([TILES_L, yb, TILES_W, h])
    data = np.array([[cell(g, c, "logFC_clip") for c in COLS] for g in genes], dtype=float)
    ax.imshow(data, aspect="auto", cmap=cmap, norm=norm, origin="upper",
              extent=[-0.5, 1.5, n - 0.5, -0.5])
    # significance / ns text per tile
    for gi, g in enumerate(genes):
        for xi, c in enumerate(COLS):
            s = cell(g, c, "sig")
            ax.text(xi, gi, "" if pd.isna(s) else str(s), ha="center", va="center",
                    fontsize=S.FS_LEGEND, color="#222222")
    ax.set_xticks([]); ax.set_yticks(range(n))
    ax.set_yticklabels(genes, fontsize=S.FS_TICK - 1, fontstyle="italic")
    ax.tick_params(axis="y", length=0)
    for sp in ax.spines.values(): sp.set_edgecolor("#888888"); sp.set_linewidth(0.6)
    # grey compartment strip + label in the left gutter
    fig.patches.append(Rectangle((0.03, yb), 0.16, h, transform=fig.transFigure,
                                 facecolor="#ECECEC", edgecolor="none", zorder=0))
    clab = comp.replace(" ", "\n") if (" " in comp and len(comp) > 12) else comp
    fig.text(0.11, yb + h / 2, clab, ha="center", va="center", fontsize=S.FS_LEGEND)
    if ci == len(COMPARTS) - 1:
        axes_bottom = ax
    y = yb - GAP_ROWS * u

# column headers under the bottom compartment
axes_bottom.set_xticks([0, 1]); axes_bottom.set_xticklabels(COLS, fontsize=S.FS_TICK)
axes_bottom.tick_params(axis="x", length=0)

# colorbar
cax = fig.add_axes([TILES_L + TILES_W + 0.06, BV + 0.10, 0.03, 0.34])
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
cb = fig.colorbar(sm, cax=cax); cb.set_ticks([-2, -1, 0, 1, 2])
cb.ax.tick_params(labelsize=S.FS_LEGEND)
cb.ax.set_title("log2 FC\n(MDT vs ND)\nclip ±2", fontsize=S.FS_LEGEND, pad=4, loc="left")

fig.suptitle("Neurovascular signature, MDT vs ND",
             fontsize=S.FS_TITLE, color=S.COL_TEXT, y=0.965)
# Star key. TWO fixes, 2026-07-17:
#  1. Said "p", but `sig` in the source table is the within-signature limma FDR
#     (adj.P.Val), not a raw p. That mislabel is what starred pecam1 at adj.P=0.067.
#     Verified: '**' = adj.P <= 0.0217, '*' = 0.0548-0.0668, 'ns' >= 0.6875 (nothing
#     lands between 0.07 and 0.69, so the thresholds are not marginal).
#  2. Was at fig coords (0.03, 0.03), i.e. the figure's bottom-left corner. The deck's
#     3D is three HAND-CROPS of the tile blocks plus a separately-placed colorbar, so
#     the corner footnote was cropped away and the deck showed */** with no definition
#     anywhere on the slide. Anchored under the tile block now so it travels with the
#     panel. PLACE 3D AS ONE IMAGE -- re-cropping into pieces will drop this again.
fig.text(TILES_L, 0.035, "*  FDR < 0.10      **  FDR < 0.05",
         fontsize=S.FS_LEGEND - 1, color="#666666", ha="left")
fig.savefig(OUT / "fig3D_neurovascular_signature.png", dpi=S.DEFAULT_DPI)
plt.close(fig)
print("wrote fig3D_neurovascular_signature.png |", n_total, "genes (pericyte dropped)")
