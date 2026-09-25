"""Shared clustered/locked heatmap renderer for the brain voom z-score panels
(2F CT/MDT-temporal DEGs, 3C mTORC1 LE, 4D MYC LE).

Design decision (2026-07 consolidation, ends the iter-14 drift):
- 2F is LOCKED to recorded artifacts: rows = tables/brain/brain_heatmap_CTtemporal_DEGs_gene_order.tsv
  (R ward.D2 leaf order), columns = the canonical 12-vector read off the canonical
  sampleclust reference. scipy re-clustering does NOT reproduce R's leaf order, so we do
  not re-cluster; we lock the order and draw order-honoring guide trees (ordered_linkage).
- 3C/4D LE heatmaps show the 15 recorded top_genes (from the panel meta.yaml) and cluster
  live (ward/euclidean); no fixed order was ever pinned for them, so a deterministic scipy
  clustering is faithful to the recorded method and stable across iterations.

Expression matrix = voom_E_BRAIN_primary (per-homeolog .L/.S), z-scored per row across the
12 brain samples. Palette RdBu_r clipped +-ZMAX_HEATMAP. Everything reads from tables/ with a
fallback to the figure_rebuild pull so it runs regardless of copy timing.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.colors import ListedColormap
from matplotlib.transforms import Bbox
from scipy.cluster.hierarchy import linkage, dendrogram
import _style as S

ROOT = Path(__file__).resolve().parent.parent
_PULL = ROOT / "figure_rebuild_2026-07" / "nibi_export" / "fig_export"

# canonical 2F column order, read off brain_heatmap_CTtemporal_DEGs_sampleclust.png
# (7 DPT block then 24 HPT block; timepoint-first, from R ward.D2/euclidean)
CANON_2F_COLS = ["S19", "S27", "S23", "S20", "S24", "S28",
                 "S30", "S26", "S22", "S29", "S25", "S21"]


def _find(name, subdir="brain"):
    for p in (ROOT / "tables" / subdir / name, ROOT / "tables" / name, _PULL / name):
        if p.exists():
            return p
    raise FileNotFoundError(f"{name} not found in tables/{subdir}, tables/, or the pull")


def load_brain():
    """voom E (columns renamed to S##) and metadata indexed by S## (trt = MDT/ND)."""
    E = pd.read_csv(_find("voom_E_BRAIN_primary.tsv"), sep="\t", index_col=0)
    E.columns = [c.split("_")[-1] for c in E.columns]              # -> S##
    md = pd.read_csv(_find("sample_metadata.tsv", subdir="."), sep="\t")
    md = md[md.tissue == "BRAIN"].copy()
    md["trt"] = md.treatment.map({"CT": "MDT", "ND": "ND"})
    md["sid"] = md["sample"].str.split("_").str[-1]
    md = md.set_index("sid")
    return E, md


def concise(sid, md):
    """S19 -> MDT_T1_7DPT (drop tissue + S#; CT->MDT; HPA/DPA -> HPT/DPT)."""
    r = md.loc[sid]
    tp = str(r["time"]).replace("HPA", "HPT").replace("DPA", "DPT")
    return f"{r['trt']}_{r['rep']}_{tp}"


def ordered_linkage(X):
    """Adjacent-merge agglomeration that PRESERVES input order. Returns a scipy linkage.
    Merge the adjacent cluster pair with min centroid-euclidean distance each step; heights
    are real distances (monotone-enforced). Leaves come out in input order 0..n-1, so the
    tree aligns with a locked row/column order. A legit guide-tree for a pre-ordered heatmap."""
    n = X.shape[0]
    segs = [[i] for i in range(n)]
    cid = list(range(n))
    cent = [X[i].astype(float).copy() for i in range(n)]
    Z, nxt, hprev = [], n, 0.0
    while len(segs) > 1:
        d = [float(np.linalg.norm(cent[i] - cent[i + 1])) for i in range(len(segs) - 1)]
        k = int(np.argmin(d))
        hprev = max(hprev, d[k]) + 1e-9
        Z.append([cid[k], cid[k + 1], hprev, len(segs[k]) + len(segs[k + 1])])
        merged = segs[k] + segs[k + 1]
        cnew = (cent[k] * len(segs[k]) + cent[k + 1] * len(segs[k + 1])) / len(merged)
        segs[k:k + 2] = [merged]; cid[k:k + 2] = [nxt]; cent[k:k + 2] = [cnew]
        nxt += 1
    return np.array(Z, dtype=float)


def render_heatmap(E, md, genes, title, out_fp, *, lock_cols=None, cluster_rows=True,
                   gene_labels=False, figsize=(5.2, 6.0), zmax=None,
                   title_size=None, xtick_size=None, ytick_size=None, legend_size=None,
                   hm_w=None):
    """Render one z-score heatmap.
    lock_cols   : list of S## for a LOCKED column order (2F). None -> cluster columns live.
    cluster_rows: True -> scipy ward/euclidean row clustering (3C/4D). False -> keep `genes`
                  order and draw an order-honoring guide tree (2F, locked to gene_order.tsv).
    gene_labels : show italic gene names on the right (LE heatmaps).
    title_size / xtick_size (sample labels) / ytick_size (gene labels) / legend_size:
                  explicit point sizes (2026-09-06, Levin round 2). Defaults are the DECK_FS_*
                  sizes, i.e. the previous behaviour. The gutters follow the sizes actually
                  used, so 3C/4D can carry larger sample labels at the deck slot size.
    hm_w        : heatmap width as a figure fraction (default 0.42 with gene labels, 0.50 for 2F).
    """
    zmax = S.ZMAX_HEATMAP if zmax is None else zmax
    xtick_size = S.DECK_FS_TICK if xtick_size is None else xtick_size
    ytick_size = xtick_size if ytick_size is None else ytick_size
    legend_size = S.DECK_FS_LEGEND if legend_size is None else legend_size
    title_size = S.DECK_FS_TITLE if title_size is None else title_size
    genes = [g for g in genes if g in E.index]
    cols0 = lock_cols if lock_cols is not None else list(E.columns)
    sub = E.loc[genes, cols0]
    z = sub.sub(sub.mean(1), axis=0).div(sub.std(1).replace(0, np.nan), axis=0).dropna()
    genes = list(z.index)

    # ---- linkages + leaf orders ----
    if cluster_rows:
        Zr = linkage(z.values, method="ward", metric="euclidean")
    else:
        Zr = ordered_linkage(z.values)                       # rows locked to `genes`
    if lock_cols is not None:
        Zc = ordered_linkage(z.values.T)                     # cols locked to CANON
    else:
        Zc = linkage(z.values.T, method="ward", metric="euclidean")
    rl = dendrogram(Zr, no_plot=True)["leaves"]
    cl = dendrogram(Zc, no_plot=True)["leaves"]

    z = z.iloc[rl, cl]                                        # reorder data to leaves
    genes = [genes[i] for i in rl]
    cols = [z.columns[j] for j in range(len(cl))]            # already reordered
    trt_of = [md.loc[s, "trt"] for s in cols]
    tp_of  = [md.loc[s, "time"] for s in cols]

    # ---- layout ----
    fig = plt.figure(figsize=figsize)
    # 2026-09-06: gutters scale with FS_SCALE so the larger review-round labels (sample
    # names along the bottom, gene names on the right) neither clip nor run into the legends.
    sc = S.FS_SCALE
    tsc = xtick_size / 7.0     # sample-label size relative to the 7 pt base (== sc by default)
    lsc = legend_size / 8.0    # legend size relative to the 8 pt base
    if gene_labels:            # LE heatmaps: room on the right for gene names
        HM_L = 0.14 + 0.02 * (tsc - 1)
        HM_W = 0.42 if hm_w is None else hm_w
        HM_B = 0.25 + 0.10 * (tsc - 1)           # rotated sample labels run ~7.6 em long
        HM_H = 0.52 - max(0.0, HM_B - 0.25)      # the top edge stays put as the gutter grows
        RD_L, RD_W = 0.035, 0.09
        # right gutter = longest gene label (~4.8 em at ytick_size) + a gap, figure fraction
        LEG_X = HM_L + HM_W + 0.03 + (4.8 * ytick_size / 72.0) / figsize[0]
    else:                      # 2F: no gene labels, x labels at bottom
        HM_L, HM_B, HM_H = 0.16, 0.135 + 0.1 * (sc - 1), 0.62 - 0.1 * (sc - 1)
        HM_W = 0.50 if hm_w is None else hm_w
        RD_L, RD_W = 0.045, 0.105
        LEG_X = HM_L + HM_W + 0.025   # snug gap from the heatmap edge
    ANN_H = 0.017 + 0.02 * (tsc - 1)         # annotation bars grow with the bar-label font
    BAR_FS = legend_size - (1.5 if gene_labels else 1.0)   # Treatment/Timepoint bar labels
    TP_B  = HM_B + HM_H + 0.004
    TRT_B = TP_B + ANN_H + 0.003
    CD_B  = TRT_B + ANN_H + 0.006

    ax    = fig.add_axes([HM_L, HM_B, HM_W, HM_H])
    ax_rd = fig.add_axes([RD_L, HM_B, RD_W, HM_H]); ax_rd.axis("off")
    ax_cd = fig.add_axes([HM_L, CD_B, HM_W, 0.075]); ax_cd.axis("off")
    ax_trt = fig.add_axes([HM_L, TRT_B, HM_W, ANN_H])
    ax_tp  = fig.add_axes([HM_L, TP_B,  HM_W, ANN_H])

    dendrogram(Zr, orientation="left", no_labels=True, ax=ax_rd,
               color_threshold=0, above_threshold_color="#555555")
    dendrogram(Zc, orientation="top", no_labels=True, ax=ax_cd,
               color_threshold=0, above_threshold_color="#555555")
    ax_rd.invert_yaxis()                                      # leaf0 -> top, matches origin='upper'

    im = ax.imshow(z.values, aspect="auto", cmap=S.CMAP_DIVERGING, vmin=-zmax, vmax=zmax,
                   interpolation="nearest", origin="upper")
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels([concise(s, md) for s in cols], rotation=90,
                       fontsize=xtick_size, va="top")
    ax.tick_params(axis="x", length=0, pad=2)
    if gene_labels:
        ax.yaxis.tick_right()
        ax.set_yticks(np.arange(len(genes)))
        ax.set_yticklabels(genes, fontsize=ytick_size, fontstyle="italic")
        ax.tick_params(axis="y", length=0, pad=2)
    else:
        ax.set_yticks([])

    # ---- annotation bars ----
    ax_trt.imshow([[0 if t == "MDT" else 1 for t in trt_of]], aspect="auto",
                  cmap=ListedColormap([S.COL_MDT, S.COL_ND]))
    ax_tp.imshow([[0 if t == "24HPA" else 1 for t in tp_of]], aspect="auto",
                 cmap=ListedColormap([S.TP_COL["24HPA"], S.TP_COL["7DPA"]]))
    for a in (ax_trt, ax_tp): a.set_xticks([]); a.set_yticks([])
    fig.text(HM_L - 0.006, TRT_B + ANN_H / 2, "Treatment", ha="right", va="center",
             fontsize=BAR_FS)
    fig.text(HM_L - 0.006, TP_B + ANN_H / 2, "Timepoint", ha="right", va="center",
             fontsize=BAR_FS)

    # ---- colorbar + legends (right), stacked top->bottom: Treatment, Timepoint, z-score.
    # Drawn MANUALLY (not via Legend) so the swatch left edges, the titles, and the colorbar
    # bar all sit at exactly LEG_X. A Legend insets its swatches by borderpad (~0.4*fontsize),
    # which put them right of the colorbar bar -> the "top two too far out" misalignment.
    if gene_labels:
        cb_h, trt_yf, tp_yf, zs_yf = 0.14, 1.00, 0.66, 0.34
    else:
        cb_h, trt_yf, tp_yf, zs_yf = 0.20, 0.90, 0.60, 0.30
    # legend row pitch follows the font scale so bigger review-round text does not overlap
    SW, SH = 0.018, 0.022 * lsc
    ROW_DY, TITLE_GAP = 0.033 * lsc, 0.036 * lsc   # SW shared with the colorbar width

    def swatch_legend(title, entries, y):
        fig.text(LEG_X, y, title, ha="left", va="center", fontsize=legend_size,
                 fontweight="bold", color="#222222")
        for i, (color, label) in enumerate(entries):
            ry = y - TITLE_GAP - i * ROW_DY                  # swatch bottom (figure fraction)
            fig.add_artist(Rectangle((LEG_X, ry), SW, SH, transform=fig.transFigure,
                                     facecolor=color, edgecolor="none"))
            fig.text(LEG_X + SW + 0.008, ry + SH / 2, label, ha="left", va="center",
                     fontsize=legend_size, color="#222222")

    swatch_legend("Treatment", [(S.COL_MDT, "MDT"), (S.COL_ND, "ND")], HM_B + HM_H * trt_yf)
    swatch_legend("Timepoint", [(S.TP_COL["24HPA"], "24 HPT"), (S.TP_COL["7DPA"], "7 DPT")],
                  HM_B + HM_H * tp_yf)

    zs_y = HM_B + HM_H * zs_yf
    fig.text(LEG_X, zs_y, "z-score", ha="left", va="center", fontsize=legend_size,
             fontweight="bold", color="#222222")
    cax = fig.add_axes([LEG_X, zs_y - TITLE_GAP - cb_h, SW, cb_h])
    cb = fig.colorbar(im, cax=cax)
    cb.set_ticks([-zmax, -zmax / 2, 0, zmax / 2, zmax])
    cb.ax.tick_params(labelsize=legend_size, length=2, pad=1.5)
    cb.outline.set_linewidth(0.5)

    fig.suptitle(title, fontsize=title_size, color="#222222", x=HM_L + HM_W / 2, y=0.965)
    # explicit full-figure Bbox: bbox_inches=None would fall back to rcParams "tight",
    # which re-maps the transFigure-anchored legends off the colorbar. Force fixed figsize.
    fig.savefig(out_fp, dpi=S.DEFAULT_DPI, bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    plt.close(fig)
    print(f"wrote {Path(out_fp).name} | {z.shape[0]} genes x {z.shape[1]} samples "
          f"| {'locked' if not cluster_rows else 'clustered'} rows, "
          f"{'locked' if lock_cols else 'clustered'} cols")
    return list(z.index), cols
