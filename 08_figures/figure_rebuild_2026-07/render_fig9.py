"""Fig 9 - Cross-dataset comparators: the treated frog brain is not a generic
distal-injury / chronic-pain response.

  9A  Walker axolotl brain scRNA UMAP, 23 Louvain clusters (cluster 4 = etv1+ winner)
  9B  Walker etv1+/glutamatergic UCell score by cluster (ranked; cluster 4 leads)
  9C  X. laevis brain bulk UCell against the Walker etv1+ signature (MDT vs ND x time)
  9D  Full 50 MSigDB Hallmark NES across 8 comparator contrasts (Ward.D2 dendrogram)

House-style rebuild (Arial, FS_*, MDT palette, RdBu_r) of the v9 canonical R panels:
  fig6a_walker_umap_clusters.R / fig6c_walker_ucell_ranked_bar.R /
  fig6e_xlaevis_walker_etv1_ucell.R / round3_fig7d.R.
Values are reproduced from LOCAL tables and VERIFIED against the recorded metas:
  9C  trt x time beta = -0.028, p = 0.024 (panel_meta_e.yaml). NO FDR is shown: this is a
      single prespecified test, so BH returns p unchanged and printing "FDR" would imply a
      multiple-testing correction that was never performed. Contrast 6C, where 26 cell-type
      signatures ARE corrected and the FDR is real and reported (2026-08-07).
  9D  Ward.D2 column order + xenopus-columns-sister (panel_meta_d.yaml, round3)
  9B  cluster-4 mean 0.0777 rank 1, cluster-7 rank 2 (fig6c_cluster_ranking.tsv)
9A UMAP coordinates were pulled once from the Walker Seurat .rds on Nibi to
tables/comparators/walker/walker_umap_coords_v2.tsv (19827 cells, cluster per cell).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.transforms import Bbox
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from scipy.stats import rankdata, t as tdist
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import _style as S

BASE = Path(__file__).resolve().parent.parent
WALK = BASE / "tables" / "comparators" / "walker"
HALL = BASE / "supplementary_v9" / "supplementary_v9_staging" / "hallmark_50_v9redfix"
VOOM = BASE / "tables" / "brain" / "voom_E_BRAIN_primary.tsv"
META = BASE / "tables" / "sample_metadata.tsv"
OUT = Path(__file__).resolve().parent / "panels_png"; OUT.mkdir(exist_ok=True)
FPNG = BASE / "figures_png"; FPNG.mkdir(exist_ok=True)

COL_ND, COL_MDT = S.COL_ND, S.COL_MDT
COL_RANK2 = "#E39385"     # light coral = runner-up cluster (on-palette, clearly secondary)
COL_OTHER = "#CBCBCB"
TIME = {"24HPA": "24 HPT", "7DPA": "7 DPT"}
WALK_UMAP = WALK / "walker_umap_coords_v2.tsv"
_verify = {}


def _save(fig, name, figsize):
    for d in (OUT, FPNG):
        fig.savefig(d / f"{name}.png", dpi=S.DEFAULT_DPI, bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    plt.close(fig)
    print("wrote", name)


# ----------------------------------------------------------------- 9A Walker UMAP
def panel_9A(figsize=(5.9, 2.7)):
    if not WALK_UMAP.exists():                          # placeholder until coords pulled
        fig = plt.figure(figsize=figsize)
        ax = fig.add_axes([0.02, 0.04, 0.96, 0.88]); ax.axis("off")
        ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor="#F4F4F4",
                               edgecolor="#BBBBBB", linewidth=1.0, linestyle=(0, (4, 3))))
        ax.text(0.5, 0.58, "Walker axolotl brain scRNA UMAP\n23 Louvain clusters (cluster 4 enlarged)",
                ha="center", va="center", fontsize=S.FS_AXIS, color="#777777", linespacing=1.5)
        ax.text(0.5, 0.24, "[embedding coordinates pending one-time Nibi pull]",
                ha="center", va="center", fontsize=S.FS_ANNOT, color="#999999", style="italic")
        fig.suptitle("Walker axolotl brain scRNA, 23 clusters", fontsize=S.FS_TITLE,
                     color=S.COL_TEXT, y=0.985)
        _save(fig, "fig9A_walker_umap_placeholder", figsize)
        return

    df = pd.read_csv(WALK_UMAP, sep="\t")
    clusters = sorted(df["cluster"].unique())          # 0..22
    _verify["9A_ncells"] = len(df); _verify["9A_nclusters"] = len(clusters)
    base = mpl.colormaps["tab20"](np.linspace(0, 1, 20))   # == ggsci d3 category20
    ramp = mpl.colors.LinearSegmentedColormap.from_list("d3ext", base)  # colorRampPalette
    pal = {c: ramp(i / (len(clusters) - 1)) for i, c in enumerate(clusters)}
    TARGET = 4
    bg = df[df["cluster"] != TARGET]; hl = df[df["cluster"] == TARGET]

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.105, 0.15, 0.50, 0.74])         # UMAP plot; legend needs ~0.38 at review fonts
    ax.scatter(bg["UMAP_1"], bg["UMAP_2"], c=[pal[c] for c in bg["cluster"]],
               s=1.3, alpha=0.7, linewidths=0)
    ax.scatter(hl["UMAP_1"], hl["UMAP_2"], c=[pal[TARGET]] * len(hl),
               s=4.5, alpha=0.9, linewidths=0)
    ax.text(hl["UMAP_1"].median(), hl["UMAP_2"].median(), "4", ha="center", va="center",
            fontsize=S.FS_TITLE, fontweight="bold", color="black")
    # UMAP axes are unitless: let the embedding fill the rect width (aspect not locked)
    # so the x-axis matches 9B/9C instead of boxing to a narrow square.
    ax.set_xlabel("UMAP 1", fontsize=S.FS_AXIS, labelpad=1)
    ax.set_ylabel("UMAP 2", fontsize=S.FS_AXIS, labelpad=1)
    ax.tick_params(labelsize=S.FS_TICK - 1)
    handles = [Line2D([], [], marker="o", linestyle="", markersize=3.2,
                      markeredgewidth=0, color=pal[c], label=str(c)) for c in clusters]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.015, 0.50), ncol=2,
              frameon=False, fontsize=S.FS_LEGEND - 1.5, title="Cluster",
              title_fontsize=S.FS_LEGEND, handletextpad=0.3, columnspacing=0.7,
              labelspacing=0.28, borderaxespad=0.0)
    fig.suptitle("Walker axolotl brain scRNA, 23 clusters", fontsize=S.FS_TITLE,
                 color=S.COL_TEXT, x=(0.105 + 0.605) / 2 + 0.03, y=0.975)
    _save(fig, "fig9A_walker_umap", figsize)


# ----------------------------------------------------------------- 9B ranked UCell bar
def panel_9B(figsize=(4.0, 2.7)):
    uc = pd.read_csv(WALK / "ucell_signatures_per_cluster_v2.tsv", sep="\t")
    g = uc.groupby("seurat_clusters")["glut_etv1plus_regen_UCellv2"]
    rk = pd.DataFrame({"mean": g.mean(), "sem": g.std(ddof=1) / np.sqrt(g.count())})
    rk = rk.sort_values("mean", ascending=False).reset_index()
    rk["rank"] = np.arange(1, len(rk) + 1)
    _verify["9B_c4_mean"] = float(rk.loc[rk["seurat_clusters"] == 4, "mean"].iloc[0])
    _verify["9B_c4_sem"] = float(rk.loc[rk["seurat_clusters"] == 4, "sem"].iloc[0])
    _verify["9B_rank1"] = int(rk.iloc[0]["seurat_clusters"])
    _verify["9B_rank2"] = int(rk.iloc[1]["seurat_clusters"])

    def col(c):
        return COL_MDT if c == 4 else COL_RANK2 if c == 7 else COL_OTHER
    colors = [col(c) for c in rk["seurat_clusters"]]

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.125, 0.235, 0.85, 0.60])
    x = np.arange(len(rk))
    ax.bar(x, rk["mean"], width=0.78, color=colors, edgecolor="#5A5A5A", linewidth=0.3, zorder=2)
    ax.errorbar(x, rk["mean"], yerr=rk["sem"], fmt="none", ecolor="#4A4A4A",
                elinewidth=0.5, capsize=1.6, capthick=0.5, zorder=3)
    ax.set_xticks(x); ax.set_xticklabels(rk["seurat_clusters"].astype(int), fontsize=S.FS_TICK - 2)
    ax.set_xlim(-0.7, len(rk) - 0.3)
    ax.set_ylim(0, (rk["mean"] + rk["sem"]).max() * 1.10)
    ax.set_xlabel("Cluster (ranked by mean UCell, descending)", fontsize=S.FS_AXIS, labelpad=2)
    ax.set_ylabel("Mean UCell score", fontsize=S.FS_AXIS, labelpad=3)
    ax.tick_params(axis="y", labelsize=S.FS_TICK)
    fig.suptitle("Walker etv1+ UCell score by cluster", fontsize=S.FS_TITLE - 0.5,
                 color=S.COL_TEXT, x=(0.125 + 0.975) / 2, y=0.965)
    _save(fig, "fig9B_walker_ucell_ranked_bar", figsize)


# ----------------------------------------------------------------- 9C X. laevis UCell boxplot
WALKER_SIG = ["etv1.S", "slc17a7.L", "slc17a7.S", "cdh13.L", "cdh13.S", "bcl11b.L",
              "bcl11b.S", "sox5.L", "sox5.S", "neto1.L", "neto1.S", "stk39.L", "stk39.S"]


def _ucell(expr, sig_rows, max_rank=1500):
    """UCell (Andreatta & Carmona 2021) per column: rank genes descending (1 = highest),
    cap at maxRank+1, AUC = 1 - (sum_rank - n(n+1)/2)/(n*maxRank). Deterministic."""
    scores = []
    n = len(sig_rows)
    for j in range(expr.shape[1]):
        r = rankdata(-expr[:, j], method="average")[sig_rows].astype(float)
        r[r > max_rank] = max_rank + 1
        u = r.sum() - n * (n + 1) / 2.0
        scores.append(1.0 - u / (n * max_rank))
    return np.array(scores)


def panel_9C(figsize=(5.6, 2.7)):
    ev = pd.read_csv(VOOM, sep="\t").set_index("gene")
    assert all(g in ev.index for g in WALKER_SIG), "signature gene missing from voom E"
    sig_rows = np.array([ev.index.get_loc(g) for g in WALKER_SIG])
    expr = ev.to_numpy()
    ucell = _ucell(expr, sig_rows)
    dat = pd.DataFrame({"sample": ev.columns, "ucell": ucell}).merge(
        pd.read_csv(META, sep="\t"), on="sample", how="left")
    dat = dat[dat["tissue"] == "BRAIN"].copy()

    # lm: ucell ~ treatment*time + rep  (ND / 24HPA / T1 reference), numpy OLS == R lm
    X = pd.DataFrame({
        "Intercept": 1.0,
        "treatmentCT": (dat["treatment"] == "CT").astype(float).values,
        "time7DPA": (dat["time"] == "7DPA").astype(float).values,
        "treatmentCT:time7DPA": ((dat["treatment"] == "CT") & (dat["time"] == "7DPA")).astype(float).values,
        "repT2": (dat["rep"] == "T2").astype(float).values,
        "repT3": (dat["rep"] == "T3").astype(float).values,
    })
    Xm = X.to_numpy(); y = dat["ucell"].to_numpy()
    beta, *_ = np.linalg.lstsq(Xm, y, rcond=None)
    resid = y - Xm @ beta
    dof = len(y) - Xm.shape[1]
    sigma2 = (resid @ resid) / dof
    covb = sigma2 * np.linalg.inv(Xm.T @ Xm)
    ix = list(X.columns).index("treatmentCT:time7DPA")
    b_int = beta[ix]; se_int = np.sqrt(covb[ix, ix])
    t_int = b_int / se_int
    p_int = 2 * tdist.sf(abs(t_int), dof)
    _verify["9C_beta"] = float(b_int); _verify["9C_p"] = float(p_int)

    order = [("ND", "24HPA"), ("ND", "7DPA"), ("CT", "24HPA"), ("CT", "7DPA")]
    grp = {(t, tm): dat[(dat["treatment"] == t) & (dat["time"] == tm)]["ucell"].values
           for t, tm in order}
    _verify["9C_means"] = {f"{t}_{tm}": float(v.mean()) for (t, tm), v in grp.items()}

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.135, 0.175, 0.54, 0.66])
    xpos = {"24HPA": 0, "7DPA": 1}; dx = {"ND": -0.19, "CT": 0.19}
    for trt, col in (("ND", COL_ND), ("CT", COL_MDT)):
        for tm in ("24HPA", "7DPA"):
            xc = xpos[tm] + dx[trt]; v = grp[(trt, tm)]
            bp = ax.boxplot(v, positions=[xc], widths=0.30, patch_artist=True,
                            showfliers=False, zorder=2,
                            medianprops=dict(color=col, linewidth=1.1),
                            whiskerprops=dict(color=col, linewidth=0.8),
                            capprops=dict(color=col, linewidth=0.8),
                            boxprops=dict(facecolor=col, alpha=0.28, edgecolor=col, linewidth=0.9))
            ax.scatter(np.full(len(v), xc) + np.linspace(-0.05, 0.05, len(v)), v,
                       s=15, color=col, edgecolors="white", linewidths=0.3, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels([TIME[t] for t in ("24HPA", "7DPA")])
    ax.set_xlim(-0.5, 1.5)
    ax.set_xlabel("Timepoint", fontsize=S.FS_AXIS, labelpad=2)
    ax.set_ylabel("UCell score", fontsize=S.FS_AXIS, labelpad=3)
    ax.tick_params(labelsize=S.FS_TICK)
    handles = [Patch(facecolor=COL_ND, alpha=0.5, edgecolor=COL_ND, label="ND"),
               Patch(facecolor=COL_MDT, alpha=0.5, edgecolor=COL_MDT, label="MDT")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              frameon=False, fontsize=S.FS_LEGEND, title="Treatment",
              title_fontsize=S.FS_LEGEND, handlelength=1.1)
    # No FDR is printed. This is a single prespecified test, so BH returns p unchanged and
    # an "FDR" line would imply a correction that was never performed. 6C is the contrast:
    # there 26 cell-type signatures are corrected, so its FDR is real and is shown.
    ax.text(1.02, 0.34, f"trt × time\nβ = {b_int:+.3f}\np = {p_int:.3f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=S.FS_ANNOT,
            color="#333333", linespacing=1.4)
    fig.suptitle("Brain UCell vs Walker etv1+ signature", fontsize=S.FS_TITLE,
                 color=S.COL_TEXT, x=(0.135 + 0.675) / 2 + 0.06, y=0.965)
    _save(fig, "fig9C_xlaevis_walker_ucell", figsize)


# ----------------------------------------------------------------- 9D full-50 hallmark heatmap
XEN24 = "X. laevis brain 24 HPT (this study)"
XEN7 = "X. laevis brain 7 DPT (this study)"
BRAIN_T = BASE / "tables" / "brain"

# X. laevis columns read the v10method tables, NOT the hallmark_50_v9redfix
# `hallmark_nes_BRAIN_*_full.tsv` pair, which is quarantined as CORRUPT
# (tables/brain/_quarantine_corrupt/WHY_CORRUPT.txt, verified on Nibi 2026-06-07):
# that pair flips 5 hallmarks against every legitimate run, including
# HEME_METABOLISM 24 HPT (+0.861 n.s. -> really -1.381 sig) and MYOGENESIS
# 7 DPT (+0.915 n.s. -> really -1.437 sig), i.e. it contradicted Fig 3A / 4A
# inside the same deck. v10method is canonical per D1 (2026-07-16).
#
# The mouse/Erdei columns stay on the v9redfix run: the corruption is specific
# to the two BRAIN tables. Checked against the independent clean descalzi
# tables on every term they share (mPFC 37/50, PAG 27/50, NAc 9/50): ZERO sign
# flips, |dNES| <= 0.22 (v9redfix runs ~10% hotter). Direction is what 9D reads,
# and direction is intact. Erdei has no local per-gene DE, so its columns can
# only come from this run in any case.
COMP = [  # (path, column label)  -- xen columns flagged "(this study)", brain lowercase
    (BRAIN_T / "hallmark_nes_BRAIN_24HPT_v10method.tsv", XEN24),
    (BRAIN_T / "hallmark_nes_BRAIN_7DPT_v10method.tsv", XEN7),
    (HALL / "gsea_mPFC_full.tsv", "Mouse mPFC SNI"),
    (HALL / "gsea_NAc_full.tsv", "Mouse NAc SNI"),
    (HALL / "gsea_PAG_full.tsv", "Mouse PAG SNI"),
    (HALL / "gsea_hallmark_bi_vs_nt_full.tsv", "Mouse cortex burn vs nt"),
    (HALL / "gsea_hallmark_fa_vs_nt_full.tsv", "Mouse cortex formalin vs nt"),
    (HALL / "gsea_hallmark_bi_vs_fa_full.tsv", "Mouse cortex burn vs formalin"),
]
# NOTE: the recorded `ward_d2_order_full50` lock (panel_meta_d.yaml) is VOID -- it
# was fitted to the corrupt X. laevis columns. Column order is now re-clustered
# from the corrected matrix on every render. Re-read the dendrogram before
# writing the caption; it is the panel's argument.


def _short(h):
    # title case (2026-09-06): narrower than caps, so the 50 labels fit the gutter at the
    # larger font; matches the 3A/4A dot-plot label style
    return h.replace("HALLMARK_", "").replace("_", " ").title()


def panel_9D(figsize=(6.0, 10.0)):
    # 2026-09-06: taller (8.6 -> 10.0 in) so 50 rows carry ~10 pt labels at the review-round
    # font sizes; slide 9's left column (A/B/C) grows to 3.15 in each to match.
    frames = []
    for path, lab in COMP:
        d = pd.read_csv(path, sep="\t")[["term", "NES", "padj"]].copy()
        d["column"] = lab
        frames.append(d)
    long = pd.concat(frames, ignore_index=True)
    hallmarks = sorted(long["term"].unique())
    assert len(hallmarks) == 50, f"expected 50 hallmarks, got {len(hallmarks)}"

    nes = long.pivot_table(index="term", columns="column", values="NES", aggfunc="mean").reindex(hallmarks)
    padj = long.pivot_table(index="term", columns="column", values="padj", aggfunc="mean").reindex(hallmarks)

    # Ward.D2 (== scipy 'ward') on column NES profiles, NA -> 0.
    # Re-clustered live: the old locked order was fitted to the corrupt inputs.
    cols_all = [lab for _, lab in COMP]
    M = nes[cols_all].to_numpy().T
    M0 = np.nan_to_num(M, nan=0.0)
    Z = linkage(pdist(M0, metric="euclidean"), method="ward")
    col_order = [cols_all[i] for i in leaves_list(Z)]
    _verify["9D_col_order"] = col_order
    _verify["9D_xen_are_sisters"] = abs(col_order.index(XEN24) - col_order.index(XEN7)) == 1
    _verify["9D_nsig"] = int(((padj < 0.05) & nes.notna()).to_numpy().sum())
    for _t, _tp in [("HALLMARK_HEME_METABOLISM", XEN24), ("HALLMARK_MYOGENESIS", XEN7)]:
        _verify[f"9D_{_t.replace('HALLMARK_', '')}_{_tp[15:21].strip()}"] = (
            round(float(nes.loc[_t, _tp]), 3), round(float(padj.loc[_t, _tp]), 4))

    # dendrogram segments a la round3_fig7d.R, leaf x-positions = locked order
    pos = {lab: i + 1 for i, lab in enumerate(col_order)}       # 1..8 left->right
    leaf_x = [pos[cols_all[i]] for i in range(len(cols_all))]   # x of scipy leaf i
    n = len(cols_all)
    node_x = {}
    seg = []
    for step in range(Z.shape[0]):
        a, b = int(Z[step, 0]), int(Z[step, 1]); h = Z[step, 2]
        xa = leaf_x[a] if a < n else node_x[a]
        xb = leaf_x[b] if b < n else node_x[b]
        ya = 0.0 if a < n else Z[a - n, 2]
        yb = 0.0 if b < n else Z[b - n, 2]
        seg += [((xa, ya), (xa, h)), ((xb, yb), (xb, h)), ((xa, h), (xb, h))]
        node_x[n + step] = (xa + xb) / 2.0

    rows = list(hallmarks)                     # alphabetical; row 0 at top
    nrow, ncol = len(rows), len(col_order)
    Mp = nes.reindex(index=rows, columns=col_order).to_numpy()
    Pp = padj.reindex(index=rows, columns=col_order).to_numpy()

    cmap = mpl.colormaps[S.CMAP_DIVERGING].copy(); cmap.set_bad("#EBEBEB")
    norm = mpl.colors.Normalize(vmin=-3, vmax=3)

    fig = plt.figure(figsize=figsize)
    L, Wd = 0.35, 0.525                        # heatmap left / width (label gutter + cbar)
    axd = fig.add_axes([L, 0.888, Wd, 0.050])  # dendrogram band
    # bottom margin 0.155 -> 0.21 (2026-09-06, round 2): the longest 45-degree column label
    # ("X. laevis brain 24 HPT (this study)", ~2.6 in at this size) needs ~1.9 in of drop.
    axh = fig.add_axes([L, 0.21, Wd, 0.672])   # heatmap (deep bottom margin for long angled labels)

    for (x0, y0), (x1, y1) in seg:
        axd.plot([x0, x1], [y0, y1], color="#4A4A4A", lw=0.6)
    axd.set_xlim(0.5, ncol + 0.5); axd.set_ylim(0, Z[:, 2].max() * 1.05); axd.axis("off")

    axh.imshow(np.ma.masked_invalid(Mp), cmap=cmap, norm=norm, aspect="auto",
               extent=[0.5, ncol + 0.5, nrow + 0.5, 0.5])
    for i in range(nrow):
        for j in range(ncol):
            axh.add_patch(Rectangle((j + 0.5, i + 0.5), 1, 1, fill=False,
                                    edgecolor="#DADADA", linewidth=0.25))
            v, p = Mp[i, j], Pp[i, j]
            if np.isfinite(v) and np.isfinite(p) and p < 0.05:
                axh.text(j + 1, i + 1, f"{v:+.1f}", ha="center", va="center",
                         fontsize=S.FS_LEGEND - 1.5, color="black")
    axh.set_xlim(0.5, ncol + 0.5); axh.set_ylim(nrow + 0.5, 0.5)
    axh.set_xticks(range(1, ncol + 1))
    axh.set_xticklabels(col_order, rotation=45, ha="right", rotation_mode="anchor",
                        fontsize=S.FS_LEGEND)
    axh.set_yticks(range(1, nrow + 1))
    axh.set_yticklabels([_short(r) for r in rows], fontsize=S.FS_LEGEND - 0.5)
    axh.tick_params(length=0)
    for s in axh.spines.values():
        s.set_visible(False)

    cax = fig.add_axes([L + Wd + 0.022, 0.45, 0.017, 0.19])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_ticks([-3, 0, 3]); cb.set_label("NES (clamp ±3)", fontsize=S.FS_LEGEND)
    cb.ax.tick_params(labelsize=S.FS_LEGEND)
    fig.text(L + Wd / 2, 0.982, "Full 50 Hallmark NES across comparator contrasts",
             ha="center", fontsize=S.FS_TITLE, color=S.COL_TEXT)
    _save(fig, "fig9D_hallmark_comparator_heatmap", figsize)


if __name__ == "__main__":
    panel_9A(figsize=(5.9, 3.15))
    panel_9B(figsize=(4.0, 3.15))
    panel_9C(figsize=(5.6, 3.15))
    panel_9D()
    for k, v in _verify.items():
        print(f"  VERIFY {k}: {v}")
    print("\n=== VERIFY ===")
    print(f"9B cluster-4 mean {_verify['9B_c4_mean']:.5f} sem {_verify['9B_c4_sem']:.5f} "
          f"(meta 0.07768 / 0.001123); rank1={_verify['9B_rank1']} rank2={_verify['9B_rank2']} (meta 4 / 7)")
    print(f"9C beta {_verify['9C_beta']:+.5f} p {_verify['9C_p']:.5f} "
          f"(meta beta -0.02797 p 0.02413)")
    for k, v in _verify["9C_means"].items():
        print(f"     mean {k}: {v:.4f}")
    print("     meta means ND_24=0.1710 ND_7=0.1691 CT_24=0.1765 CT_7=0.1466")
    print(f"9D n_sig cells {_verify['9D_nsig']} (meta 74); col order {_verify['9D_col_order']}; "
          f"xen sisters: {_verify['9D_xen_are_sisters']}")
