"""Fig 4 local panels (late brain, 7 DPT) -- 2026-09-06, Levin review round 2.

Stewart: Fig 4's text read smaller than every other figure (its 2 x 3 grid was landscape, so
the page placed it at 0.52x) and the six panel titles came out in different sizes, weights
and fonts, because half the panels were rendered by figures_src (bold 14 pt titles, tight
bbox) and half by this rebuild. Fix: a PORTRAIT 2-column grid (9.4 in wide, page scale ~0.67)
with every panel rendered at its deck slot size and placed 1:1, all under XENO_FS_SCALE=1.15,
so titles are FS_TITLE regular, ticks FS_TICK, legends FS_LEGEND everywhere in the figure
(~7.7 / 6.9 / 6.2 pt on the page). 4B/4C/4D stay in figures_src (fig4b/fig4c/fig4d) but now
take the same sizes explicitly.

  4A  Hallmark NES dot plot, MDT vs ND, 7 DPT   LOCKED to the 11 recorded terms
      (figures_current/v10_pulls/brain/brain_hallmark_7dpt_deseq2_meta.yaml), as fig4_hallmark_7dpt.py.
  4E  ME1 (turquoise) module Hallmark over-representation   LOCKED to the 5 recorded terms
      (brain_me1_hallmark_overrep_meta.yaml), as fig4_me1_hallmark_overrep.py (local wgcna table drifted).
  4F  Hallmark NES 24 HPT vs 7 DPT state-transition scatter, all 50 Hallmarks from the
      v10method tables (tables/brain/hallmark_nes_BRAIN_{24HPT,7DPT}_v10method.tsv),
      as fig4f_brain_state_transition.py. Labels: the four sign-flip pathways sit in the empty
      left half (no Hallmark has NES_24 < -1.9), the sustained-early set in the right margin.
The figures_src originals are kept for provenance and are no longer in the driver.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.transforms import Bbox
import _style as S
from _panels import nes_dotplot

BASE = Path(__file__).resolve().parent.parent
TBL = BASE / "tables" / "brain"
OUT = Path(__file__).resolve().parent / "panels_png"; OUT.mkdir(exist_ok=True)
FPNG = BASE / "figures_png"; FPNG.mkdir(exist_ok=True)


def _save(fig, name, figsize):
    for d in (OUT, FPNG):
        fig.savefig(d / f"{name}.png", dpi=S.DEFAULT_DPI, bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- 4A Hallmark dot plot (locked)
# (ID, NES, padj) from brain_hallmark_7dpt_deseq2_meta.yaml : terms
TERMS_4A = [
    ("HALLMARK_HEME_METABOLISM",           2.09, 2.97e-11),
    ("HALLMARK_ALLOGRAFT_REJECTION",       2.00, 6.01e-07),
    ("HALLMARK_INTERFERON_GAMMA_RESPONSE", 1.87, 1.32e-05),
    ("HALLMARK_MYC_TARGETS_V1",            1.79, 3.34e-06),
    ("HALLMARK_INTERFERON_ALPHA_RESPONSE", 1.78, 1.63e-03),
    ("HALLMARK_IL2_STAT5_SIGNALING",       1.62, 8.86e-04),
    ("HALLMARK_INFLAMMATORY_RESPONSE",     1.55, 1.35e-02),
    ("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", 1.52, 1.87e-02),
    ("HALLMARK_IL6_JAK_STAT3_SIGNALING",   1.52, 4.85e-02),
    ("HALLMARK_P53_PATHWAY",               1.44, 1.87e-02),
    ("HALLMARK_MYOGENESIS",               -1.44, 1.35e-02),
]


def panel_4A(figsize=(4.4, 3.3)):
    df = pd.DataFrame(TERMS_4A, columns=["term", "NES", "padj"])
    nes_dotplot(df, "Brain Hallmark GSEA, 7 DPT", OUT / "fig4A_hallmark_7dpt.png",
                top_n=20, figsize=figsize, left=0.56, label_size=S.FS_TICK - 1)
    (FPNG / "fig4A_hallmark_7dpt.png").write_bytes((OUT / "fig4A_hallmark_7dpt.png").read_bytes())


# ---------------------------------------------------------------- 4E ME1 module over-rep (locked)
# (ID, overlap Count, padj, fold_enrich) from brain_me1_hallmark_overrep_meta.yaml : terms,
# ordered least -> most significant so barh puts the most significant (G2M) on top.
TERMS_4E = [
    ("MYC_TARGETS_V1",  11, 3.25e-02, 2.15),
    ("MITOTIC_SPINDLE", 24, 8.57e-07, 2.62),
    ("HEME_METABOLISM", 36, 5.07e-09, 2.48),
    ("E2F_TARGETS",     46, 1.61e-11, 2.48),
    ("G2M_CHECKPOINT",  44, 8.35e-12, 2.59),
]
ME1_N = 1372


def panel_4E(figsize=(3.8, 3.0)):
    d = pd.DataFrame(TERMS_4E, columns=["label", "Count", "padj", "fe"])
    d["log2fe"] = np.log2(d["fe"])
    d["pretty"] = d["label"].str.replace("_", " ").str.title()
    norm = mpl.colors.Normalize(vmin=d["log2fe"].min(), vmax=d["log2fe"].max())
    cmap = mpl.colormaps[S.CMAP_ENRICHMENT]

    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=0.38, right=0.86, top=0.87, bottom=0.17)
    ax.barh(d["pretty"], d["Count"], color=[cmap(norm(v)) for v in d["log2fe"]],
            edgecolor="black", linewidth=0.5, height=0.7)
    for y, (cnt, q) in enumerate(zip(d["Count"], d["padj"])):
        ax.text(cnt + d["Count"].max() * 0.02, y, f"q = {q:.1e}", va="center", ha="left",
                fontsize=S.FS_LEGEND - 1, color="#333333")
    ax.set_xlim(0, d["Count"].max() * 1.7)
    ax.set_xlabel("Overlap count", labelpad=2)
    ax.tick_params(axis="y", labelsize=S.FS_TICK - 1)
    ax.grid(True, axis="x", alpha=0.3)
    fig.suptitle(f"ME1 module Hallmark enrichment (n = {ME1_N:,})", fontsize=S.FS_TITLE,
                 color=S.COL_TEXT, y=0.965)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.05, pad=0.03)
    cb.set_label(r"$\log_2$ fold enrichment", fontsize=S.FS_LEGEND)
    cb.ax.tick_params(labelsize=S.FS_LEGEND)
    _save(fig, "fig4E_me1_hallmark_overrep", figsize)


# ---------------------------------------------------------------- 4F state-transition scatter
PRETTY = {
    "OXIDATIVE_PHOSPHORYLATION": "OXPHOS",
    "IL6_JAK_STAT3_SIGNALING":   "IL6/JAK/STAT3",
    "TNFA_SIGNALING_VIA_NFKB":   "TNFA via NF-kB",
    "EPITHELIAL_MESENCHYMAL_TRANSITION": "EMT",
    "MYC_TARGETS_V1":            "MYC targets V1",
    "INTERFERON_ALPHA_RESPONSE": "Interferon\nalpha",
    "INTERFERON_GAMMA_RESPONSE": "Interferon\ngamma",
    "ALLOGRAFT_REJECTION":       "Allograft\nrejection",
    "HEME_METABOLISM":           "Heme\nmetabolism",
    "HYPOXIA":                   "Hypoxia",
    "COAGULATION":               "Coagulation",
    "MYOGENESIS":                "Myogenesis",
}
# label anchor (data coords) + alignment; leader lines run from the point to the anchor.
# Left column: two-line labels right-aligned at x = -2.15, ordered by the point's y so the
# leaders never cross (the empty half-plane: no Hallmark has NES_24 below -1.9).
# Right column: single-line labels left-aligned at x = 2.3, same rule.
X_L, X_R = -2.15, 2.3
LABEL_POS = {
    "HEME_METABOLISM":            (X_L, 2.6, "right"),
    "ALLOGRAFT_REJECTION":        (X_L, 1.9, "right"),
    "INTERFERON_GAMMA_RESPONSE":  (X_L, 1.2, "right"),
    "INTERFERON_ALPHA_RESPONSE":  (X_L, 0.5, "right"),
    "MYC_TARGETS_V1":             (X_R, 2.4, "left"),
    "IL6_JAK_STAT3_SIGNALING":    (X_R, 2.0, "left"),
    "COAGULATION":                (X_R, 1.6, "left"),
    "HYPOXIA":                    (X_R, 1.2, "left"),
    "TNFA_SIGNALING_VIA_NFKB":    (X_R, 0.8, "left"),
    "EPITHELIAL_MESENCHYMAL_TRANSITION": (X_R, 0.4, "left"),
    "OXIDATIVE_PHOSPHORYLATION":  (X_R, 0.0, "left"),
    "MYOGENESIS":                 (X_R, -1.45, "left"),
}


def panel_4F(figsize=(5.0, 4.0)):
    a = pd.read_csv(TBL / "hallmark_nes_BRAIN_24HPT_v10method.tsv", sep="\t").rename(
        columns={"NES": "NES_24", "padj": "padj_24"})
    b = pd.read_csv(TBL / "hallmark_nes_BRAIN_7DPT_v10method.tsv", sep="\t").rename(
        columns={"NES": "NES_7", "padj": "padj_7"})
    m = a.merge(b, on="term")[["term", "NES_24", "padj_24", "NES_7", "padj_7"]].copy()
    m["short"] = m["term"].str.replace("HALLMARK_", "", regex=False)
    m["delta"] = m["NES_7"] - m["NES_24"]
    m["sig_any"] = m[["padj_24", "padj_7"]].min(axis=1) < 0.05
    assert m["NES_24"].min() > -1.95, "left label column assumes no point left of -1.95"

    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=0.125, right=0.80, top=0.90, bottom=0.14)
    xlo, xhi, ylo, yhi = -4.8, 5.9, -2.3, 3.0
    ax.axhline(0, color="#999999", lw=0.7)
    ax.axvline(0, color="#999999", lw=0.7)
    ax.plot([-2.2, 2.9], [-2.2, 2.9], ls="--", color="#CCCCCC", lw=0.9, zorder=0)

    cmap = mpl.colormaps[S.CMAP_DIVERGING]            # red = positive delta = gained late
    norm = mpl.colors.Normalize(vmin=-2.5, vmax=2.5)
    ns, sig = m[~m["sig_any"]], m[m["sig_any"]]
    ax.scatter(ns["NES_24"], ns["NES_7"], s=14, c="#DDDDDD", edgecolors="#BBBBBB",
               linewidths=0.4, zorder=2)
    sc = ax.scatter(sig["NES_24"], sig["NES_7"], s=40, c=sig["delta"], cmap=cmap, norm=norm,
                    edgecolors="black", linewidths=0.6, zorder=3)

    idx = m.set_index("short")
    for short, (xt, yt, ha) in LABEL_POS.items():
        if short not in idx.index or not idx.loc[short, "sig_any"]:
            continue
        r = idx.loc[short]
        ax.plot([r["NES_24"], xt], [r["NES_7"], yt], color="#AAAAAA", lw=0.6, zorder=2)
        ax.annotate(PRETTY.get(short, short.replace("_", " ").capitalize()), (xt, yt),
                    fontsize=S.FS_LEGEND, ha=ha, va="center", color="#222222", zorder=4,
                    linespacing=0.95, xytext=(-2 if ha == "right" else 2, 0),
                    textcoords="offset points")

    ax.set_xlabel("Hallmark NES (MDT vs ND), 24 HPT (early)", labelpad=2)
    ax.set_ylabel("Hallmark NES (MDT vs ND), 7 DPT (late)", labelpad=2)
    ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi)
    ax.grid(True, alpha=0.2)
    # contrast named on the axes, so the title stays clear of the colorbar's "gained late"
    fig.suptitle("Brain state transition, early to late", fontsize=S.FS_TITLE,
                 color=S.COL_TEXT, x=0.46, y=0.965)

    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.04)
    cb.set_label("NES shift (7 DPT minus 24 HPT)", fontsize=S.FS_LEGEND)
    cb.ax.tick_params(labelsize=S.FS_LEGEND)
    cb.ax.text(0.5, 1.03, "gained late", transform=cb.ax.transAxes, ha="center", va="bottom",
               fontsize=S.FS_LEGEND - 1, color=S.COL_MDT)
    cb.ax.text(0.5, -0.03, "lost late", transform=cb.ax.transAxes, ha="center", va="top",
               fontsize=S.FS_LEGEND - 1, color=S.COL_ND)
    _save(fig, "fig4F_brain_state_transition", figsize)
    print(f"   {int(m['sig_any'].sum())} sig pathways highlighted, "
          f"{sum(1 for s in LABEL_POS if s in idx.index and idx.loc[s, 'sig_any'])} labelled")


if __name__ == "__main__":
    panel_4A()
    panel_4E()
    panel_4F()
