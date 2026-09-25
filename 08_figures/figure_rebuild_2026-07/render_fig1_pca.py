"""Fig 1E/1F — Brain and Blastema PCA, matplotlib reproduction of the approved
deck codex panels, with Stewart's 2026-07-07 edits:
  - title centered, Arial 10 pt, no subtitle
  - x-stretched ~1.25x (landscape) to tuck under the schematic
  - MDT relabel (was CT)
  - convex hulls per (treatment, timepoint); ND blue / MDT coral; time by shape
Data: manuscript_lock_2026-06-08/tables/{brain,blastema}/pca_*_strand_corrected.tsv
Variance % is computed from the same PCA (2026-09-25): prcomp-equivalent on voom E of the
homeolog-summed counts (nibi_export/fig_export/voom_E_{BRAIN,BLAST}.tsv), 5,000 most variable
genes, centered, unscaled, as in fig1D_brain_pca.R / fig2C_blast_pca.R. The recomputed PC1/PC2
scores must match the plotted score tables. The labels used before this date (brain 24.5/12.6,
blast 46.9/15.9) were typed in from an older deck and did not match the plotted PCA
(true values: brain 23.7/12.5, blast 46.5/14.8).
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Polygon
from matplotlib.lines import Line2D
from scipy.spatial import ConvexHull
import _style as S

BASE = Path(__file__).resolve().parent.parent
TBL  = BASE / "manuscript_lock_2026-06-08" / "tables"
VOOM = Path(__file__).resolve().parent / "nibi_export" / "fig_export"
OUT  = Path(__file__).resolve().parent / "panels_png"; OUT.mkdir(exist_ok=True)


def variance_explained(scores, voom_tsv, n_top=5000):
    """Recompute the plotted PCA (R: prcomp(t(E[top, ]), center=TRUE, scale.=FALSE) on the n_top genes
    with the largest sample variance) and return % variance of PC1 and PC2 after checking its scores."""
    E = pd.read_csv(voom_tsv, sep="\t", index_col=0)[scores["sample"]]
    top = E.var(axis=1, ddof=1).sort_values(ascending=False, kind="mergesort").index[:n_top]
    X = E.loc[top].T.values
    X = X - X.mean(axis=0)
    U, s, _ = np.linalg.svd(X, full_matrices=False)
    pcs = U * s
    for k in (0, 1):
        a, b = pcs[:, k], scores[f"PC{k + 1}"].values
        assert np.allclose(np.sign(a @ b) * a, b, atol=1e-6), f"PC{k + 1} scores do not match {voom_tsv.name}"
    ve = s ** 2 / np.sum(s ** 2) * 100
    return ve[0], ve[1]

# PCA-only palette: Stewart wants a clean red/blue here (the global MDT coral read
# as "orangeish" on these panels). Blue stays the house ND blue.
PCA_COL = {"ND": S.COL_ND, "CT": "#D62728"}


def pca_panel(df, title, varexp, time_shapes, out_fp, figsize):
    treats = ["ND", "CT"]                       # plot order (ND under MDT)
    times  = list(time_shapes.keys())
    fig, ax = plt.subplots(figsize=figsize)
    # reserve a bottom band: axis label gets its own row above the single-row legend
    fig.subplots_adjust(left=0.15, right=0.97, top=0.91, bottom=0.27)

    for trt in treats:
        col = PCA_COL[trt]
        for tp in times:
            sub = df[(df["treatment"] == trt) & (df["time"] == tp)]
            if not len(sub):
                continue
            pts = sub[["PC1", "PC2"]].values
            # convex hull (needs >=3 non-collinear points)
            if len(pts) >= 3:
                try:
                    h = ConvexHull(pts)
                    poly = Polygon(pts[h.vertices], closed=True, facecolor=col,
                                   edgecolor=S.COL_HULL_EDGE, alpha=0.18, lw=0.6,
                                   zorder=1)
                    ax.add_patch(poly)
                except Exception:
                    pass
            ax.scatter(sub["PC1"], sub["PC2"], color=col,
                       marker=time_shapes[tp], s=55, edgecolors="black",
                       linewidths=0.6, zorder=3)

    ax.set_xlabel(f"PC1 ({varexp[0]:.1f}%)", labelpad=2)
    ax.set_ylabel(f"PC2 ({varexp[1]:.1f}%)", labelpad=2)
    ax.set_title(title, loc="center", fontsize=S.FS_TITLE, color=S.COL_TEXT, pad=4)
    ax.axhline(0, color="grey", lw=0.4, alpha=0.4, zorder=0)
    ax.axvline(0, color="grey", lw=0.4, alpha=0.4, zorder=0)

    # ---- single-row legend: treatment by colour, timepoint by grey shape ----
    treat_handles = [Patch(facecolor=PCA_COL[t], edgecolor="black", lw=0.5,
                           label=S.TREAT_LABEL[t]) for t in treats]
    time_handles = [Line2D([0], [0], marker=time_shapes[tp], color="none",
                           markerfacecolor="#555555", markeredgecolor="black",
                           markersize=7, label=S.TIME_LABEL[tp]) for tp in times]
    handles = treat_handles + time_handles
    fig.legend(handles=handles, ncol=len(handles), loc="lower center",
               bbox_to_anchor=(0.5, 0.015), frameon=False, fontsize=S.FS_LEGEND,
               handletextpad=0.3, columnspacing=1.0)

    fig.savefig(out_fp, dpi=S.DEFAULT_DPI)   # exact figsize, no bbox='tight'
    plt.close(fig)
    print("wrote", out_fp.name, f"{figsize[0]}x{figsize[1]}in AR={figsize[0]/figsize[1]:.2f}")


# final panel size for a 180 mm (7.09 in) double-column Fig 1: two PCAs share the
# bottom row -> 3.47 in wide each (7.09 total, matches schematic width). Height set
# for the bottom 2-row legend to breathe; still wider than the codex portrait panel.
FIGSIZE = (3.47, 2.88)

brain = pd.read_csv(TBL / "brain/pca_BRAIN_homeolog_sum_strand_corrected.tsv", sep="\t")
ve_brain = variance_explained(brain, VOOM / "voom_E_BRAIN.tsv")
print(f"brain PC1 {ve_brain[0]:.2f}%, PC2 {ve_brain[1]:.2f}%")
pca_panel(brain, "Brain", ve_brain, S.TIME_SHAPE_BRAIN,
          OUT / "fig1E_brain_pca.png", FIGSIZE)

blast = pd.read_csv(TBL / "blastema/pca_BLAST_homeolog_sum_strand_corrected.tsv", sep="\t")
print("blast times:", sorted(blast["time"].unique()), "treatments:", sorted(blast["treatment"].unique()))
ve_blast = variance_explained(blast, VOOM / "voom_E_BLAST.tsv")
print(f"blastema PC1 {ve_blast[0]:.2f}%, PC2 {ve_blast[1]:.2f}%")
pca_panel(blast, "Blastema", ve_blast, S.TIME_SHAPE,
          OUT / "fig1F_blastema_pca.png", FIGSIZE)
