"""Supplementary Figure 1 - Paired blastema reproduces the peripheral treatment
phenotype (EMT/fibrotic suppression + cell-cycle/biosynthetic/metabolic activation),
confirming the brain changes occur in a regenerative-permissive, anti-fibrotic local state.

  S1A  Blastema Hallmark GSEA dot plot at 7 DPT (17 FDR<0.05 terms)
  S1B  Blastema MDT-vs-ND volcano at 7 DPT (EMT / anti-fibrotic genes overlaid)
  S1C  EMT DEG-intersected z-score heatmap (18 genes, HALLMARK_EMT x q<0.05 7 DPT)
  S1D  EMT leading-edge trajectories (6 genes, 2x3)
  S1E  HALLMARK_EMT running-enrichment (fgsea mountain, blastema 7 DPT)
  S1F  Anti-fibrotic gene trajectories (6 genes, 2x3)

House-style rebuild (Arial, FS_*, MDT palette, RdBu_r) of the manuscript-lock blastema
panels; values reproduced from LOCAL tables and VERIFIED against the recorded metas:
  A  17 sig terms, NES/padj from blast_hallmark_7dpt_deseq2_meta.yaml
  C  18 genes from blastema_emt_heatmap_v10_meta.yaml (fixed ND|MDT x 11/24/7 columns)
  D  col1a2.S/edil3.S/lum.L/thbs2.L/fbln2.L/loxl2.S; per-timepoint MDT-vs-ND stars
  E  NES = -2.11, FDR = 6.34e-09 (blastema_mountain_EMT_7dpt_meta.yaml)
  F  6 of 8 anti-fibrotic genes (dropped tgfb2.S + col3a1.S per Stewart 2026-07-15);
     per-timepoint MDT-vs-ND stars, NOT a single trt x time beta (3 timepoints)

Blastema has 3 timepoints, so D/F use per-timepoint DESeq2 Wald stars (like Fig 5B),
not the single trt x time beta of the 2-timepoint brain panels (Figs 5A/6B).
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.transforms import Bbox
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import _style as S
import _traj as T

BASE = Path(__file__).resolve().parent.parent
sys.path.append(str(BASE / "figures_src"))       # for _mountain (running_es / load_pathway)
import _mountain as M                              # noqa: E402
S.apply()                                          # figure_rebuild style wins (Arial/FS_*/MDT)

BL = BASE / "tables" / "blastema"
META = BASE / "tables" / "sample_metadata.tsv"
VOOM = BL / "voom_E_BLAST_primary.tsv"
OUT = Path(__file__).resolve().parent / "panels_png"; OUT.mkdir(exist_ok=True)
FPNG = BASE / "figures_png"; FPNG.mkdir(exist_ok=True)

TIMES = ["11HPA", "24HPA", "7DPA"]
TP_COL = {"11HPA": "#8C6BB1", "24HPA": "#41AE76", "7DPA": "#E6842A"}   # timepoint annotation
COL_ND, COL_MDT = S.COL_ND, S.COL_MDT
_verify = {}
DE_BL = {tp: pd.read_csv(BL / f"de_BLAST_CT_vs_ND_{tp}.tsv", sep="\t").set_index("gene")
         for tp in TIMES}


def _save(fig, name, figsize):
    for d in (OUT, FPNG):
        fig.savefig(d / f"{name}.png", dpi=S.DEFAULT_DPI, bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    plt.close(fig)
    print("wrote", name)


def de_padj_bl(g, tp):
    d = DE_BL[tp]
    return float(d.loc[g, "padj"]) if g in d.index else np.nan


# =============================================================== S1A hallmark dot plot
# 17 FDR<0.05 terms locked to blast_hallmark_7dpt_deseq2_meta.yaml (clusterProfiler GSEA)
HALLMARK_A = [
    ("HALLMARK_CHOLESTEROL_HOMEOSTASIS", 2.19, 2.62e-10),
    ("HALLMARK_G2M_CHECKPOINT", 2.09, 1.75e-14),
    ("HALLMARK_E2F_TARGETS", 2.02, 2.07e-12),
    ("HALLMARK_MYC_TARGETS_V2", 1.99, 3.24e-06),
    ("HALLMARK_MTORC1_SIGNALING", 1.98, 1.25e-11),
    ("HALLMARK_MITOTIC_SPINDLE", 1.71, 6.5e-06),
    ("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", 1.66, 0.00134),
    ("HALLMARK_MYC_TARGETS_V1", 1.60, 8.9e-05),
    ("HALLMARK_TNFA_SIGNALING_VIA_NFKB", 1.48, 0.00518),
    ("HALLMARK_ANDROGEN_RESPONSE", 1.45, 0.0372),
    ("HALLMARK_DNA_REPAIR", 1.45, 0.0326),
    ("HALLMARK_ESTROGEN_RESPONSE_LATE", 1.44, 0.0241),
    ("HALLMARK_PI3K_AKT_MTOR_SIGNALING", 1.43, 0.0407),
    ("HALLMARK_IL2_STAT5_SIGNALING", 1.38, 0.0413),
    ("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", -1.76, 0.000104),
    ("HALLMARK_INTERFERON_GAMMA_RESPONSE", -2.04, 6.5e-06),
    ("HALLMARK_INTERFERON_ALPHA_RESPONSE", -2.61, 3.17e-09),
]


def panel_A(figsize=(5.0, 4.7)):
    # 2026-09-06: 3.95 -> 5.0 in wide with a 0.52 label gutter so the longest Hallmark name
    # fits at the review-round font size (Supplementary Fig. 1 now = A/B over C).
    d = pd.DataFrame(HALLMARK_A, columns=["ID", "NES", "padj"])
    d["label"] = d["ID"].str.replace("HALLMARK_", "", regex=False).str.replace("_", " ").str.title()
    d = d.sort_values("NES").reset_index(drop=True)             # ascending -> highest NES at top
    d["mlogq"] = -np.log10(d["padj"].clip(lower=1e-300))
    d["size"] = (16 + np.sqrt(d["mlogq"]) * 15).clip(upper=145)
    _verify["A_nterms"] = len(d)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.52, 0.10, 0.34, 0.82])
    sc = ax.scatter(d["NES"], np.arange(len(d)), s=d["size"], c=d["NES"],
                    cmap=S.CMAP_DIVERGING, vmin=-2.6, vmax=2.6,
                    edgecolors="black", linewidths=0.5, zorder=3)
    ax.axvline(0, color="grey", ls="--", lw=0.6)
    ax.set_yticks(range(len(d))); ax.set_yticklabels(d["label"], fontsize=S.FS_TICK - 1)
    ax.set_ylim(-0.6, len(d) - 0.4)
    ax.set_xlabel("NES (MDT vs ND)", fontsize=S.FS_AXIS, labelpad=2)
    ax.tick_params(axis="x", labelsize=S.FS_TICK)
    ax.grid(True, axis="x", alpha=0.3)
    cax = fig.add_axes([0.89, 0.30, 0.02, 0.42])
    cb = fig.colorbar(sc, cax=cax); cb.set_label("NES", fontsize=S.FS_LEGEND)
    cb.ax.tick_params(labelsize=S.FS_LEGEND)
    fig.suptitle("Blastema Hallmark GSEA at 7 DPT", fontsize=S.FS_TITLE,
                 color=S.COL_TEXT, x=(0.52 + 0.86) / 2, y=0.965)
    _save(fig, "figS1A_blast_hallmark", figsize)


# =============================================================== S1B volcano
def panel_B(figsize=(4.4, 4.7)):
    de = pd.read_csv(BL / "de_BLAST_CT_vs_ND_7DPA.tsv", sep="\t").dropna(
        subset=["padj", "pvalue", "log2FoldChange"]).copy()
    de["mlogp"] = -np.log10(de["pvalue"].clip(lower=1e-300))
    LFC, Q, CAP = 1.0, 0.05, 7.0
    de["x"] = de["log2FoldChange"].clip(-CAP, CAP)       # clamp outliers onto the +-CAP axes
    up = (de["log2FoldChange"] > LFC) & (de["padj"] < Q)
    dn = (de["log2FoldChange"] < -LFC) & (de["padj"] < Q)
    ns = ~(up | dn)

    le = pd.read_csv(BL / "blast_LE_genes_per_hallmark_CT_vs_ND_7DPT.tsv", sep="\t")
    emt = set(le.loc[le["pathway"] == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "gene_xla"])
    af = {"tgfb2.L", "tgfb2.S", "sox9.L", "mmp13.S", "col17a1.L", "col1a1.L", "col3a1.S", "fn1.S"}
    overlay = emt | af
    ov = de[de["gene"].isin(overlay)].copy()
    ov_sig = ov[ov["padj"] < Q]; ov_ns = ov[ov["padj"] >= Q]
    _verify["B_up"] = int(up.sum()); _verify["B_dn"] = int(dn.sum())
    _verify["B_overlay_sig"] = int(len(ov_sig))

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.14, 0.15, 0.83, 0.72])
    ax.scatter(de.loc[ns, "x"], de.loc[ns, "mlogp"], s=5, color="#CFCFCF",
               alpha=0.45, linewidths=0, zorder=1)
    ax.scatter(de.loc[dn, "x"], de.loc[dn, "mlogp"], s=7, color=COL_ND,
               alpha=0.8, linewidths=0, label="Down in MDT vs ND", zorder=2)
    ax.scatter(de.loc[up, "x"], de.loc[up, "mlogp"], s=7, color=COL_MDT,
               alpha=0.8, linewidths=0, label="Up in MDT vs ND", zorder=2)
    ax.scatter(ov_ns["x"], ov_ns["mlogp"], s=24, marker="^", facecolor="none",
               edgecolor="#5E3C99", linewidths=0.8, label="Anti-fibrotic / EMT", zorder=4)
    ax.scatter(ov_sig["x"], ov_sig["mlogp"], s=28, marker="^", color="#F5B301",
               edgecolor="black", linewidths=0.5, label="Anti-fibrotic / EMT (sig)", zorder=5)
    ax.axvline(0, color="#999999", lw=0.6, zorder=0)
    for xv in (-LFC, LFC):
        ax.axvline(xv, color="#CCCCCC", ls="--", lw=0.6, zorder=0)

    # label the strongest overlay genes (at clamped x)
    lab = ov_sig.reindex(ov_sig["log2FoldChange"].abs().sort_values(ascending=False).index).head(8)
    try:
        from adjustText import adjust_text
        txts = [ax.text(r["x"], r["mlogp"], r["gene"], fontsize=S.FS_ANNOT - 1.5,
                        fontstyle="italic", color="#222222") for _, r in lab.iterrows()]
        adjust_text(txts, ax=ax, only_move={"points": "y", "text": "xy"},
                    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4), ensure_inside_axes=True)
    except Exception:
        for _, r in lab.iterrows():
            ax.text(r["x"], r["mlogp"] + 0.1, r["gene"], fontsize=S.FS_ANNOT - 1.5,
                    fontstyle="italic", ha="center", color="#222222")

    ax.set_xlim(-CAP - 0.4, CAP + 0.4)
    ax.set_xlabel(r"$\log_2$ fold change (MDT vs ND)", fontsize=S.FS_AXIS, labelpad=2)
    ax.set_ylabel(r"$-\log_{10}\ p$", fontsize=S.FS_AXIS, labelpad=3)
    ax.tick_params(labelsize=S.FS_TICK)
    ax.legend(loc="upper left", bbox_to_anchor=(-0.02, 1.0), frameon=False,
              fontsize=S.FS_LEGEND - 1.5, handletextpad=0.3, labelspacing=0.25, borderaxespad=0.1)
    fig.suptitle("Blastema MDT vs ND at 7 DPT", fontsize=S.FS_TITLE,
                 color=S.COL_TEXT, x=(0.14 + 0.97) / 2, y=0.965)
    _save(fig, "figS1B_blast_volcano", figsize)


# =============================================================== S1C EMT DEG heatmap
HEAT18 = ["mmp1.S", "tnfrsf11b.L", "cxcl8.L", "loxl2.S", "lum.L", "serpine2.L", "thbs2.L",
          "fbln2.L", "cdh11.S", "efemp2.L", "igfbp2.S", "cdh11.L", "plod1.S", "dpysl3.L",
          "pmp22.L", "p3h1.L", "gas1.S", "calu.L"]


def panel_C(figsize=(5.5, 3.7)):
    ev = pd.read_csv(VOOM, sep="\t").set_index("gene")
    mb = pd.read_csv(META, sep="\t"); mb = mb[mb["tissue"] == "BLAST"]
    order = []                                            # fixed ND|MDT x 11/24/7 x reps
    for trt in ["ND", "CT"]:
        for tp in TIMES:
            order += list(mb[(mb["treatment"] == trt) & (mb["time"] == tp)].sort_values("rep")["sample"])
    X = ev.loc[HEAT18, order].to_numpy(dtype=float)
    Z = (X - X.mean(axis=1, keepdims=True)) / X.std(axis=1, ddof=0, keepdims=True)
    Z = np.clip(Z, -3, 3)
    trt_of = {s: mb.set_index("sample").loc[s, "treatment"] for s in order}
    tp_of = {s: mb.set_index("sample").loc[s, "time"] for s in order}
    _verify["C_ngenes"] = len(HEAT18); _verify["C_nsamples"] = len(order)

    fig = plt.figure(figsize=figsize)
    axh = fig.add_axes([0.145, 0.045, 0.60, 0.80])
    axt = fig.add_axes([0.145, 0.895, 0.60, 0.028])       # treatment bar
    axp = fig.add_axes([0.145, 0.855, 0.60, 0.028])       # timepoint bar
    axh.imshow(Z, cmap=S.CMAP_DIVERGING, vmin=-3, vmax=3, aspect="auto")
    axh.set_xticks([]); axh.set_yticks(range(len(HEAT18)))
    axh.set_yticklabels(HEAT18, fontsize=S.FS_TICK - 2, fontstyle="italic")
    for s in axh.spines.values():
        s.set_visible(False)
    axt.imshow([[0 if trt_of[s] == "ND" else 1 for s in order]], aspect="auto",
               cmap=mpl.colors.ListedColormap([COL_ND, COL_MDT]))
    axp.imshow([[TIMES.index(tp_of[s]) for s in order]], aspect="auto",
               cmap=mpl.colors.ListedColormap([TP_COL[t] for t in TIMES]))
    for a, lab in ((axt, "Treatment"), (axp, "Timepoint")):
        a.set_xticks([]); a.set_yticks([])
        a.text(-0.012, 0.5, lab, transform=a.transAxes, ha="right", va="center", fontsize=S.FS_LEGEND)
        for s in a.spines.values():
            s.set_visible(False)

    cax = fig.add_axes([0.765, 0.30, 0.02, 0.42])
    cb = fig.colorbar(ScalarMappable(norm=mpl.colors.Normalize(-3, 3), cmap=S.CMAP_DIVERGING), cax=cax)
    cb.set_label("Row z-score (clamp ±3)", fontsize=S.FS_LEGEND); cb.ax.tick_params(labelsize=S.FS_LEGEND)
    th = [Patch(facecolor=COL_ND, label="ND"), Patch(facecolor=COL_MDT, label="MDT")]
    ph = [Patch(facecolor=TP_COL[t], label=S.TIME_LABEL[t]) for t in TIMES]
    leg1 = fig.legend(handles=th, title="Treatment", loc="upper left", bbox_to_anchor=(0.845, 0.86),
                      frameon=False, fontsize=S.FS_LEGEND, title_fontsize=S.FS_LEGEND, handlelength=1.0)
    fig.add_artist(leg1)
    fig.legend(handles=ph, title="Timepoint", loc="upper left", bbox_to_anchor=(0.845, 0.65),
               frameon=False, fontsize=S.FS_LEGEND, title_fontsize=S.FS_LEGEND, handlelength=1.0)
    fig.suptitle("Blastema EMT-intersected DEG expression", fontsize=S.FS_TITLE,
                 color=S.COL_TEXT, x=(0.145 + 0.745) / 2, y=0.975)
    _save(fig, "figS1C_blast_emt_heatmap", figsize)


# =============================================================== S1D / S1F trajectories (2x3)
def panel_traj(genes, title, name, figsize):
    ev = pd.read_csv(VOOM, sep="\t", index_col=0)
    mb = pd.read_csv(META, sep="\t"); mb = mb[mb["tissue"] == "BLAST"]
    lm, rm, tm, bm, gx, gy = 0.135, 0.02, 0.165, 0.09, 0.10, 0.115
    ncol = 2
    nrow = int(np.ceil(len(genes) / ncol))               # 3 rows x 2 cols (portrait)
    aw = (1 - lm - rm - (ncol - 1) * gx) / ncol
    ah = (1 - tm - bm - (nrow - 1) * gy) / nrow
    fig = plt.figure(figsize=figsize)
    handles = None
    for k, g in enumerate(genes):
        r, c = divmod(k, ncol)
        x0 = lm + c * (aw + gx)
        y0 = 1 - tm - (r + 1) * ah - r * gy
        ax = fig.add_axes([x0, y0, aw, ah])
        star_by = {tp: T.stars(de_padj_bl(g, tp)) for tp in TIMES}
        handles = T.draw_trajectory(ax, ev.loc[g], mb, TIMES, star_by, ylabel_text=None)
        ax.set_xlabel("")
        ax.set_title(g, fontsize=S.FS_TITLE - 1, fontstyle="italic", color=S.COL_TEXT, pad=3)
        ax.tick_params(labelsize=S.FS_TICK - 1)
    fig.text((lm + (1 - rm)) / 2, 0.022, "Timepoint", ha="center", va="bottom", fontsize=S.FS_AXIS)
    fig.text(0.03, (bm + (1 - tm)) / 2, "Normalized expression", rotation=90, ha="left",
             va="center", fontsize=S.FS_AXIS)
    # legend centered ABOVE the grid (between suptitle and facet titles) so it never overlaps a facet
    fig.legend([handles["ND"], handles["MDT"]], ["ND", "MDT"], title="Treatment",
               loc="upper center", bbox_to_anchor=(0.5, 0.945), frameon=False,
               fontsize=S.FS_LEGEND, title_fontsize=S.FS_LEGEND, handlelength=1.6,
               ncol=2, columnspacing=1.1)
    fig.suptitle(title, fontsize=S.FS_TITLE, color=S.COL_TEXT, x=(lm + (1 - rm)) / 2, y=0.99)
    _save(fig, name, figsize)


# =============================================================== S1E EMT mountain
def _blast_rank():
    de = pd.read_csv(BL / "de_BLAST_CT_vs_ND_7DPA.tsv", sep="\t")
    rk = np.sign(de["log2FoldChange"]) * -np.log10(de["pvalue"])
    sym = de["gene_id"].str.replace(r"\.[LS]$", "", regex=True).str.upper()
    d = pd.DataFrame({"sym": sym.values, "rk": rk.values}).replace([np.inf, -np.inf], np.nan).dropna(subset=["rk"])
    d["entrez"] = d["sym"].map(M._sym2ez())
    d = d.dropna(subset=["entrez"]).drop_duplicates("entrez", keep="first")
    return pd.Series(d["rk"].values, index=d["entrez"].values).sort_values(ascending=False)


def panel_E(figsize=(4.1, 3.5)):
    # NES/FDR read from the SAME locked list S1A plots, so the two panels cannot
    # disagree about one result. Was hardcoded -2.11 / 6.34e-09 from
    # blastema_mountain_EMT_7dpt_meta.yaml (the direct-fgsea run) while S1A showed
    # -1.76 / 1.04e-04 for the same pathway+contrast. Canonical = v10method (D1).
    # The ES curve below is unaffected: same ranks, same running-ES math.
    NES, FDR = next((n, p) for i, n, p in HALLMARK_A
                    if i == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
    rk = _blast_rank()
    xs, ys, es, S_pos, N = M.running_es(rk, M.load_pathway("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"))
    _verify["E_es"] = float(es); _verify["E_le_positions"] = int((ys.min() == ys).sum())
    color = COL_MDT if NES >= 0 else COL_ND
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.155, 0.15, 0.815, 0.72])
    ax.axhline(0, color="#BBBBBB", lw=0.6, zorder=1)
    ax.plot(xs, ys, color=color, lw=1.4, zorder=3)
    y0, y1 = min(ys.min(), 0), max(ys.max(), 0)
    span = (y1 - y0) or 1.0
    ax.vlines(S_pos, y0 - 0.13 * span, y0 - 0.04 * span, color="#444444", lw=0.35, zorder=2)
    ax.set_xlim(0, N); ax.set_ylim(y0 - 0.17 * span, y1 + 0.10 * span)
    ax.set_xlabel("Rank in ordered gene list", fontsize=S.FS_AXIS, labelpad=2)
    ax.set_ylabel("Running enrichment score", fontsize=S.FS_AXIS, labelpad=3)
    ax.tick_params(labelsize=S.FS_TICK)
    ax.annotate(f"NES = {NES:+.2f}\nFDR = {FDR:.2g}", xy=(0.96, 0.96), xycoords="axes fraction",
                ha="right", va="top", fontsize=S.FS_ANNOT, color="#333333")
    fig.suptitle("EMT Hallmark: blastema MDT vs ND, 7 DPT", fontsize=S.FS_TITLE - 1,
                 color=S.COL_TEXT, x=(0.155 + 0.97) / 2, y=0.965)
    _save(fig, "figS1E_blast_emt_mountain", figsize)


if __name__ == "__main__":
    panel_A()
    panel_B()
    panel_C()
    panel_traj(["col1a2.S", "edil3.S", "lum.L", "thbs2.L", "fbln2.L", "loxl2.S"],
               "Blastema EMT leading-edge trajectories", "figS1D_blast_emt_trajectories", (4.7, 6.6))
    panel_E()
    panel_traj(["tgfb2.L", "sox9.L", "mmp13.S", "col17a1.L", "col1a1.L", "fn1.S"],
               # Only tgfb2.L (*7), mmp13.S (**7) and col17a1.L (*11/***24) are ever
               # significant; sox9.L / col1a1.L / fn1.S never are. "Anti-fibrotic"
               # outruns that, so the title names the gene class, not a conclusion.
               "Selected fibrosis/EMT-associated gene trajectories",
               "figS1F_blast_antifibrotic", (4.7, 6.6))
    print("\n=== VERIFY ===")
    for k, v in _verify.items():
        print(f"  {k}: {v}")
