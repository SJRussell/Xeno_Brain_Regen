"""Fig 2 local panels — brain treatment-transition.
  2A brain volcano MDT vs ND @ 24 HPT   (local: de_BRAIN_CTvND_24HPT.tsv)
  2B brain volcano MDT vs ND @ 7 DPT    (local: de_BRAIN_CTvND_7DPT.tsv)
  2E temporal DEG Venn (MDT 7v24 vs ND 7v24), counts moved OFF the circles
     (local: pipeline_pulls/brain_{CT,ND}_7v24_sig.tsv)
2C/2D (temporal volcanoes) + 2F (932-gene heatmap) are Nibi-origin -> reused PNGs.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib_venn import venn2, venn2_circles
import _style as S

BASE = Path(__file__).resolve().parent.parent
TBL  = BASE / "tables"          # was manuscript_lock_2026-06-08/tables (absent here)
OUT  = Path(__file__).resolve().parent / "panels_png"; OUT.mkdir(exist_ok=True)


def _cap(df, padj_thresh=0.05):
    """Symmetric x-cap = max(4, 95th pct of |log2FC| among significant genes)."""
    d = df[df["padj"].notna() & df["log2FoldChange"].notna()]
    a = d.loc[d["padj"] < padj_thresh, "log2FoldChange"].abs()
    return max(4.0, float(np.percentile(a, 95))) if len(a) else 5.0


def _ymax(df):
    """Max -log10(padj) in a panel, for a shared y-frame across panels."""
    d = df[df["padj"].notna()]
    return float((-np.log10(d["padj"].clip(lower=1e-300))).max()) if len(d) else 1.0


def volcano(df, title, out_fp, figsize=(3.4, 3.25), padj_thresh=0.05,
            auto_label_n=8, label_genes=None, cap=None, ylim=None,
            up_label="up in MDT", down_label="down in MDT",
            xlabel="log2 fold change (MDT vs ND)", deck=False):
    """deck=True renders the DECK variant: no title, no y-axis label (the slide carries
    coloured title bars and a bracketed y-label per column), x-label only when `xlabel`
    is given. Stewart used to hand-crop these off the full panels (2026-09-06)."""
    d = df.copy()
    d = d[d["padj"].notna() & d["log2FoldChange"].notna()].copy()
    d["nlp"] = -np.log10(d["padj"].clip(lower=1e-300))
    sig = d["padj"] < padj_thresh
    up = sig & (d["log2FoldChange"] > 0)
    dn = sig & (d["log2FoldChange"] < 0)
    ns = ~(up | dn)
    if cap is None:                       # per-panel cap unless a shared one is passed
        cap = _cap(df, padj_thresh)
    d["xplot"] = d["log2FoldChange"].clip(-cap, cap)

    fig, ax = plt.subplots(figsize=figsize)
    # right margin leaves runway for right-anchored gene labels (up genes pile up at the
    # fold-change cap); with a fixed-size crop they'd otherwise clip at the figure edge.
    if deck:
        fig.subplots_adjust(left=0.13, right=0.90, top=0.97, bottom=(0.17 if xlabel else 0.09))
    else:
        fig.subplots_adjust(left=0.18, right=0.86, top=0.90, bottom=0.17)
    ax.scatter(d.loc[ns, "xplot"], d.loc[ns, "nlp"], s=5, c="#D0D0D0",
               edgecolors="none", rasterized=True)
    # deck variant: short legend entries ("up (n=573)"); the slide's coloured title bar and
    # the x-label already name the contrast, and the long entries were what collided with
    # the top genes of the temporal volcanoes (2026-09-06, round 2)
    dn_lab = f"down (n={int(dn.sum())})" if deck else f"{down_label} (n={int(dn.sum())})"
    up_lab = f"up (n={int(up.sum())})" if deck else f"{up_label} (n={int(up.sum())})"
    ax.scatter(d.loc[dn, "xplot"], d.loc[dn, "nlp"], s=11, c=S.COL_ND,
               edgecolors="none", label=dn_lab)
    ax.scatter(d.loc[up, "xplot"], d.loc[up, "nlp"], s=11, c=S.COL_MDT,
               edgecolors="none", label=up_lab)
    ax.axhline(-np.log10(padj_thresh), color="grey", ls="--", lw=0.6)
    ax.axvline(0, color="grey", lw=0.5)
    ax.set_xlim(-cap * 1.1, cap * 1.1)
    if ylim is not None:                   # shared y-frame across panels
        ax.set_ylim(*ylim)
    elif deck:
        # 12% headroom (2026-09-06, round 2): the two most significant temporal DEGs sit at
        # the fold-change cap AND the top of 2C, so their labels had nowhere to go and overlapped
        ax.set_ylim(-0.5, float(d["nlp"].max()) * 1.12)

    to_label = {}
    for _, r in d[sig].sort_values("padj").head(auto_label_n).iterrows():
        to_label[r["gene"]] = (r["xplot"], r["nlp"])
    if label_genes:
        for g in label_genes:
            row = d[d["gene"].str.upper().str.replace(r"\.[LS]$", "", regex=True) == g.upper()]
            if len(row):
                r = row.sort_values("padj").iloc[0]
                to_label[r["gene"]] = (r["xplot"], r["nlp"])
    texts = [ax.text(gx, gy, gn, fontsize=S.FS_TICK - 2, fontstyle="italic")
             for gn, (gx, gy) in to_label.items()]
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="#999999", lw=0.4),
                    expand=(1.2, 1.5), ensure_inside_axes=True)
    except Exception as e:
        print("  (adjustText:", e, ")")

    if xlabel:
        ax.set_xlabel(xlabel, labelpad=2)
    if not deck:
        ax.set_ylabel(r"$-\log_{10}$(adjusted p)", labelpad=2)
        ax.set_title(title, fontsize=S.FS_TITLE, color=S.COL_TEXT, pad=4)
    # Legend corner chosen per panel (2026-09-06, Levin round 2): a fixed upper-left legend
    # sat on the most significant genes of the temporal volcanoes (2C/2D). Score the four
    # corners by what the legend box would cover (gene labels 10 each, with a small margin;
    # significant points 5; grey points 0.02) and take the emptiest; ties keep the upper left.
    lkw = dict(fontsize=S.FS_LEGEND - 1, frameon=False, handletextpad=0.3, borderaxespad=0.3)
    leg = ax.legend(loc="upper left", **lkw)
    fig.canvas.draw()
    bb = leg.get_window_extent().transformed(ax.transAxes.inverted())
    lw, lh = bb.width, bb.height
    pts = ax.transAxes.inverted().transform(ax.transData.transform(np.c_[d["xplot"], d["nlp"]]))
    sigm = sig.to_numpy()
    lab = [t.get_window_extent().transformed(ax.transAxes.inverted()) for t in texts]

    def cover(x0, y0):
        x1, y1 = x0 + lw, y0 + lh
        inside = (pts[:, 0] >= x0) & (pts[:, 0] <= x1) & (pts[:, 1] >= y0) & (pts[:, 1] <= y1)
        c = 5.0 * inside[sigm].sum() + 0.02 * inside[~sigm].sum()
        e = 0.015                                             # label margin (axes fraction)
        c += 10.0 * sum(1 for b in lab
                        if b.x1 + e > x0 and b.x0 - e < x1 and b.y1 + e > y0 and b.y0 - e < y1)
        return c

    mg = 0.015
    corners = {"upper left": (mg, 1 - mg - lh), "upper right": (1 - mg - lw, 1 - mg - lh),
               "lower left": (mg, mg), "lower right": (1 - mg - lw, mg)}
    best = min(corners, key=lambda k: cover(*corners[k]))     # dict order breaks ties
    leg.remove()
    ax.legend(loc=best, **lkw)
    print(f"  legend -> {best}")
    ax.grid(True, alpha=0.18)
    # fixed full-figure bbox so all four volcanoes come out at identical dimensions.
    # bbox_inches="tight" was cropping to content, and 2C's auto-labels pushed its width
    # out (1049 vs 984 px). ensure_inside_axes keeps labels off the crop edge.
    fig.savefig(out_fp, dpi=S.DEFAULT_DPI, bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    plt.close(fig)
    print("wrote", out_fp.name)


def temporal_venn(out_fp, figsize=(3.7, 3.0), scale=1.25):
    # Derived from the same temporal DE tables the 2B volcanoes use, rather than the
    # pipeline_pulls/brain_{CT,ND}_7v24_sig.tsv snapshots (not on this machine).
    # Verified 2026-07-17 to reproduce the recorded Venn exactly: 925 / 7 / 2.
    # 2026-09-06 (Levin review): the 7 shared and 2 ND-only gene names are listed beside
    # their callouts, in the panel's unused right-hand space.
    def _sig(fn):
        d = pd.read_csv(TBL / "brain" / fn, sep="\t")
        return set(d.loc[d["padj"] < 0.05, "gene_id"])
    ct = _sig("de_BRAIN_CT_7DPA_vs_24HPA.tsv")
    nd = _sig("de_BRAIN_ND_7DPA_vs_24HPA.tsv")
    a_only, b_only, shared = len(ct - nd), len(nd - ct), len(ct & nd)
    shared_genes, nd_genes = sorted(ct & nd, key=str.lower), sorted(nd - ct, key=str.lower)
    print(f"venn: MDT-only={a_only} ND-only={b_only} shared={shared}")
    print("  shared:", shared_genes, "| ND-only:", nd_genes)

    def _wrap(genes, per_line):
        return "\n".join(", ".join(genes[i:i + per_line]) for i in range(0, len(genes), per_line))

    # scale the GRAPHIC only: bigger canvas -> circles grow `scale`x (they fill the axes),
    # but fonts stay at deck sizes so labels/title match the other panels (scale=1.0 = original).
    figsize = (figsize[0] * scale, figsize[1] * scale)
    fs_axis, fs_legend, fs_title = S.FS_AXIS, S.FS_LEGEND, S.FS_TITLE
    fig, ax = plt.subplots(figsize=figsize)
    fig.subplots_adjust(left=0.02, right=0.98, top=0.86, bottom=0.06)
    v = venn2(subsets=(a_only, b_only, shared), set_labels=("MDT", "ND"), ax=ax)
    # circle fills = treatment colours
    if v.get_patch_by_id("10"): v.get_patch_by_id("10").set_color(S.COL_MDT); v.get_patch_by_id("10").set_alpha(0.55)
    if v.get_patch_by_id("01"): v.get_patch_by_id("01").set_color(S.COL_ND);  v.get_patch_by_id("01").set_alpha(0.55)
    if v.get_patch_by_id("11"): v.get_patch_by_id("11").set_color("#7A5C8A"); v.get_patch_by_id("11").set_alpha(0.6)
    # capture the LIBRARY-COMPUTED region centroids before blanking the labels, so
    # leader lines point at the true region (with 925 vs 2 the circles are wildly
    # asymmetric; hand-guessed target coords land inside the big MDT circle).
    cen = {}
    for cid in ("10", "01", "11"):
        lbl = v.get_label_by_id(cid)
        cen[cid] = np.array(lbl.get_position()) if lbl else None
        if lbl: lbl.set_text("")            # hide default in-circle counts
    # set labels (MDT / ND) styling
    for sid, col in (("A", S.COL_MDT), ("B", S.COL_ND)):
        t = v.get_label_by_id(sid)
        if t: t.set_fontsize(fs_axis); t.set_fontweight("bold"); t.set_color(col)
    venn2_circles(subsets=(a_only, b_only, shared), ax=ax, lw=0.8, color="#404040")

    # MDT-only count sits inside the big circle (no leader needed)
    if cen["10"] is not None:
        ax.text(cen["10"][0], cen["10"][1], f"{a_only}\nMDT-only", ha="center",
                va="center", fontsize=fs_axis, fontweight="bold", color="#222222")

    # shared + ND-only: callouts on the white background, leaders to the true centroids;
    # the gene names sit under each callout (italic), left-aligned in the free right margin.
    def callout(text, target, xytext, genes):
        ax.annotate(text, xy=tuple(target), xytext=xytext, fontsize=fs_legend,
                    ha="left", va="bottom", fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color="#666666", lw=0.6,
                                    relpos=(0.0, 0.0)),
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none"))
        ax.text(xytext[0], xytext[1] - 0.03, _wrap(genes, 2), fontsize=fs_legend - 1,
                fontstyle="italic", ha="left", va="top", color="#222222", linespacing=1.15)
    if cen["11"] is not None:
        callout(f"{shared} shared", cen["11"], (cen["11"][0] + 0.30, cen["11"][1] + 0.62),
                shared_genes)
    if cen["01"] is not None:
        callout(f"{b_only} ND-only", cen["01"], (cen["01"][0] + 0.30, cen["01"][1] - 0.30),
                nd_genes)

    ax.set_title("Temporal DEGs, 7 DPT vs 24 HPT", fontsize=fs_title,
                 color=S.COL_TEXT, pad=6)
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    ax.set_xlim(x0 - 0.12, x1 + 1.05); ax.set_ylim(y0 - 0.55, y1 + 0.40)
    fig.savefig(out_fp, dpi=S.DEFAULT_DPI)
    plt.close(fig)
    print("wrote", out_fp.name)


# Repointed 2026-07-17 from manuscript_lock_2026-06-08 + nibi_export (neither present on
# this machine) to the restored canonical tables/. The lock used `de_BRAIN_CTvND_24HPT`,
# the restore uses `de_BRAIN_CT_vs_ND_24HPA` -- same data under a different filename.
# VERIFIED equivalent against the recorded panel counts before repointing (padj<0.05,
# which is what volcano() uses; no LFC filter):
#   24 HPT  29 DEGs = 16 up / 13 down     7 DPT  32 DEGs = 20 up / 12 down
#   MDT 7v24  932 = 573 up / 359 down     ND 7v24  9 = 2 up / 7 down
# and the 2D Venn derives to 925 MDT-only / 7 shared / 2 ND-only. All match the deck.
v24 = pd.read_csv(TBL / "brain/de_BRAIN_CT_vs_ND_24HPA.tsv", sep="\t")
v7  = pd.read_csv(TBL / "brain/de_BRAIN_CT_vs_ND_7DPA.tsv", sep="\t")
# 2C/2D temporal volcanoes (7 DPT vs 24 HPT within treatment)
mdt_t = pd.read_csv(TBL / "brain/de_BRAIN_CT_7DPA_vs_24HPA.tsv", sep="\t")
nd_t  = pd.read_csv(TBL / "brain/de_BRAIN_ND_7DPA_vs_24HPA.tsv", sep="\t")

# Each panel auto-scales its own x/y range (Stewart: data ranges do NOT need to match).
# What must match is the plot-box geometry, which is guaranteed by the fixed figsize +
# subplots_adjust + fixed-Bbox save (all four -> box [162,878]x[97,830]); so stacked in a
# column they align without per-image scaling. `cap`/`ylim` left at defaults = per-panel.
volcano(v24, "MDT vs ND, 24 HPT", OUT / "fig2A_brain_volcano_24hpt.png")
volcano(v7,  "MDT vs ND, 7 DPT",  OUT / "fig2B_brain_volcano_7dpt.png",
        label_genes=["gdf15"])
tkw = dict(up_label="up at 7 DPT", down_label="down at 7 DPT",
           xlabel="log2 fold change (7 DPT vs 24 HPT)")
volcano(mdt_t, "MDT: 7 DPT vs 24 HPT", OUT / "fig2C_brain_volcano_MDT_7v24.png", **tkw)
volcano(nd_t,  "ND: 7 DPT vs 24 HPT",  OUT / "fig2D_brain_volcano_ND_7v24.png", **tkw)

# deck variants (slide 2 column: title bars + bracket y-labels live on the slide; the
# x-label appears once per stacked pair, under the lower panel)
volcano(v24, "", OUT / "fig2A_brain_volcano_24hpt_deck.png", deck=True, xlabel=None,
        figsize=(3.4, 2.95))
volcano(v7,  "", OUT / "fig2B_brain_volcano_7dpt_deck.png", deck=True, label_genes=["gdf15"])
volcano(mdt_t, "", OUT / "fig2C_brain_volcano_MDT_7v24_deck.png", deck=True,
        up_label="up at 7 DPT", down_label="down at 7 DPT", xlabel=None, figsize=(3.4, 2.95),
        auto_label_n=6)   # 8 labels pile up at the fold-change cap at the larger font
volcano(nd_t,  "", OUT / "fig2D_brain_volcano_ND_7v24_deck.png", deck=True,
        up_label="up at 7 DPT", down_label="down at 7 DPT",
        xlabel="log2 fold change (7 DPT vs 24 HPT)")

temporal_venn(OUT / "fig2E_temporal_venn.png")
