"""Shared panel renderers for the 2026-07 rebuild (MDT-relabelled, Arial, fixed size).
Reused across figures: nes_dotplot (Hallmark/Reactome/GO NES), reactome_barh.
"""
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import _style as S


def nes_dotplot(df, title, out_fp, top_n=20, figsize=(3.9, 4.6),
                padj_col="padj", nes_col="NES", term_col="term",
                strip_prefix="HALLMARK_", xlabel="NES (MDT vs ND)", left=0.60,
                label_size=None):
    """Hallmark/Reactome/GO NES dot plot: y=pathway, x=NES, colour=NES, size=-log10(q).
    Fixed figsize (no bbox tight) for clean rectangle tiling.
    2026-09-06: default widened 3.4 -> 3.9 in and the label gutter to 0.60 so the longest
    Hallmark name still fits at the larger review-round font sizes."""
    d = df.copy().dropna(subset=[nes_col, padj_col])
    d = d.sort_values(padj_col).head(top_n)
    d["label"] = (d[term_col].str.replace(strip_prefix, "", regex=False)
                  .str.replace("_", " ").str.title())
    d["mlogq"] = -d[padj_col].apply(lambda x: math.log10(max(x, 1e-300)))
    d["size"] = (14 + d["mlogq"].pow(0.5) * 16).clip(upper=130)
    d = d.iloc[::-1].reset_index(drop=True)

    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=left, right=0.86, top=0.88, bottom=0.13)
    sc = ax.scatter(d[nes_col], d["label"], s=d["size"], c=d[nes_col],
                    cmap=S.CMAP_DIVERGING, vmin=-2, vmax=2,
                    edgecolors="black", linewidths=0.5)
    ax.axvline(0, color="grey", linestyle="--", linewidth=0.6)
    ax.set_xlabel(xlabel, labelpad=2)
    ax.set_ylabel("")
    fig.suptitle(title, fontsize=S.FS_TITLE, color=S.COL_TEXT, y=0.965)
    ax.tick_params(axis="y", labelsize=(S.FS_TICK - 1) if label_size is None else label_size)
    ax.grid(True, axis="x", alpha=0.3)
    cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("NES", fontsize=S.FS_LEGEND)
    cb.ax.tick_params(labelsize=S.FS_LEGEND)
    fig.savefig(out_fp, dpi=S.DEFAULT_DPI)
    plt.close(fig)
    print("wrote", out_fp.name)


def reactome_barh(df, title, out_fp, must_include, figsize=(5.0, 3.4), ns_q=0.05, left=0.56):
    """Horizontal bar of -log10(q) for a curated ECM/collagen Reactome set.

    `ns_q`: pathways with p.adjust >= ns_q are drawn hollow with a grey dashed edge and
    an italic grey label, so a curated-but-non-significant pathway cannot be read as a
    finding. The set is curated, so these are shown for context, not because they passed.
    """
    d = df[df["direction"] == "ANY"].copy()
    d = d[d["ID"].isin(must_include)].copy()
    d["mlogq"] = -d["p.adjust"].apply(lambda x: math.log10(max(x, 1e-300)))
    d["sig"] = d["p.adjust"] < ns_q
    d["label"] = (d["ID"].str.replace("REACTOME_", "", regex=False)
                  .str.replace("_", " ").str.capitalize())
    d["label"] = d["label"].apply(lambda s: s if len(s) <= 36 else s[:35].rstrip() + "...")
    d = d.sort_values("mlogq").reset_index(drop=True)

    norm = mpl.colors.Normalize(vmin=0.5, vmax=max(d["mlogq"].max(), 2.0))
    cmap = mpl.colormaps[S.CMAP_ENRICHMENT]
    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=left, right=0.88, top=0.88, bottom=0.16)
    bars = ax.barh(d["label"], d["mlogq"],
                   color=[cmap(norm(x)) if s else "none" for x, s in zip(d["mlogq"], d["sig"])],
                   edgecolor=["black" if s else "#9A9A9A" for s in d["sig"]],
                   linewidth=0.5)
    for bar, s in zip(bars, d["sig"]):        # barh takes one linestyle, so set per-bar
        if not s:
            bar.set_linestyle((0, (2, 1.5)))
    for tick, s in zip(ax.get_yticklabels(), d["sig"]):
        if not s:
            tick.set_color("#8A8A8A"); tick.set_fontstyle("italic")
    ax.axvline(-math.log10(ns_q), color="grey", linestyle="--", linewidth=1)
    ax.set_xlabel(r"$-\log_{10}$(adjusted p)", labelpad=2)
    fig.suptitle(title, fontsize=S.FS_TITLE, color=S.COL_TEXT, y=0.965)
    ax.set_xlim(0, max(d["mlogq"].max() * 1.12, 2.2))
    ax.tick_params(axis="y", labelsize=S.FS_TICK - 1)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label(r"$-\log_{10}(q)$", fontsize=S.FS_LEGEND)
    cb.ax.tick_params(labelsize=S.FS_LEGEND)
    ax.grid(True, axis="x", alpha=0.3)
    fig.savefig(out_fp, dpi=S.DEFAULT_DPI)
    plt.close(fig)
    print("wrote", out_fp.name)
