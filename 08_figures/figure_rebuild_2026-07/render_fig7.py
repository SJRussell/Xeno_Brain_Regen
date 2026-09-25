"""Fig 7 - Brain process timeline (A) + uninjured-baseline anchors (B/C/D).

7A  Brain process state, early to late.
    Reproduces manuscript_lock figProc_go_timeline.py panel B EXACTLY (verified 2026-07-14):
      - fgsea over all 3473 GO:BP terms at 24 HPT and 7 DPT (MDT vs ND), locked ranks
        -> tables/brain/gobp_BRAIN_CTvND_{24HPT,7DPT}.tsv (term, NES, padj)
      - 193 terms sig at either timepoint -> Wang semantic similarity + average linkage, k=12
        -> tables/brain/go_process_clusters_BRAIN.tsv (cluster members, GOBP_-prefixed)
      - process score = MEAN NES of member terms; FLAGGED when >=25% of member terms
        reach FDR<0.05. Synaptic transmission = clusters 11+5 merged.
    Restyle only (values unchanged): figure_rebuild _style (Arial, FS_*, RdBu_r, MDT not CT).

7B/7C/7D/7E  GSE73419 uninjured-baseline anchors (D3, Stewart 2026-07-16: build all three,
    he keeps what works). All data LOCAL in tables/comparators/gse73419/, no Nibi.

    READ THIS BEFORE WRITING ANY CAPTION. These are QUALITATIVE anchors, hard
    pre-registered as "no DE statistics between GSE73419 and Murugan". The quantity
    plotted is a percentile RANK SHIFT: each gene's expression percentile within our
    dataset minus its percentile within GSE73419. It is an effect-size vector, not a
    tested contrast. Caveats that must survive into the legend: n=2 baseline animals,
    different clutch/strain, library-prep mismatch (poly-A vs Ribo-Zero), read-layout
    mismatch (PE 2x100 vs SE 51).

    7B  Hallmark NES on the rank-shift vector, ND and MDT arms (dumbbell).
        Verified 2026-07-17: 50 Hallmark terms; 19 sig in EITHER arm, 14 in BOTH, and
        only 2 of 50 differ in sign between arms. The arms sit on top of each other.
    7C  Genome-wide concordance of the two arms' rank shifts.
        r = 0.987 over the 26,481 genes scored in BOTH arms (997 of the 27,478 table
        rows lack a shift in one arm and are dropped). Movers (|shift| > 0.20):
        6,119 shared of 7,109 union.
        NOTE for the legend: the shipped `mover_flag` column says 7,119, not 7,109. The
        10-gene gap is genes flagged on ONE arm whose other arm is NaN (hoxb8.L,
        hoxb6.S, rpl35a.L, gkn1.L, ...), so they cannot be plotted in 2-D or counted as
        shared. Quote 7,109 for this panel; 7,119 is the table's own flag count.
        gdf15.L is rank 39 of 26,481 by off-diagonal distance |MDT shift - ND shift|
        (NOT by |MDT shift|, where it is only rank 3,377).
    7D  gdf15 three-point anchor, expression percentile in baseline -> ND -> MDT.
        gdf15.L 0.073 -> 0.136 -> 0.370; gdf15.S 0.053 -> 0.102 -> 0.106 (flat control).
    7E  Walker etv1+ signature (UCell) by arm and timepoint vs the same baseline.
        Added 2026-08-07. baseline 0.049; ND 0.114 -> 0.110; MDT 0.115 -> 0.097, i.e. BOTH
        arms ~2.3x baseline at 24 HPT. Sibling of 7D. Scored on the combined-TPM table,
        so values are NOT comparable with Fig 9C (voom, this study only).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from matplotlib.transforms import Bbox
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
import _style as S

BASE = Path(__file__).resolve().parent.parent
TBL = BASE / "tables" / "brain"
GSE = BASE / "tables" / "comparators" / "gse73419"
OUT = Path(__file__).resolve().parent / "panels_png"; OUT.mkdir(exist_ok=True)
FPNG = BASE / "figures_png"; FPNG.mkdir(exist_ok=True)

CL = pd.read_csv(TBL / "go_process_clusters_BRAIN.tsv", sep="\t")
G24 = pd.read_csv(TBL / "gobp_BRAIN_CTvND_24HPT.tsv", sep="\t").set_index("term")
G7 = pd.read_csv(TBL / "gobp_BRAIN_CTvND_7DPT.tsv", sep="\t").set_index("term")
FRAC, FDR = 0.25, 0.05
# Label fixes 2026-07-16/17:
#  - cluster 2 was "Heme / myeloid / wound response": "myeloid" asserts a cell-type call
#    that bulk brain cannot make. -> "Heme / iron and wound response".
#  - cluster 6 was "Muscle / contraction (early)": labelled early but n.s. early and
#    flagged LATE (see the VERIFY block). "(early)" dropped.
NAME = {3: "Protein biosynthesis\n(translation, ribosome)", 1: "Immune / defense response",
        2: "Heme / iron and\nwound response", 7: "Endocrine /\noxidative-stress response",
        6: "Muscle / contraction", 8: "Cholesterol / sterol\nbiosynthesis"}
# (two-line names added 2026-09-06 so the labels fit at the larger review-round font size)
MERGE = {"Synaptic transmission": [11, 5]}
MOVER = 0.20        # |percentile rank shift| defining a "mover" (== the mover_flag column)


def agg(sub):
    terms = [f"GOBP_{t}" for m in sub["members"] for t in str(m).split("; ")]
    n24, p24 = G24["NES"].reindex(terms).dropna(), G24["padj"].reindex(terms).dropna()
    n7, p7 = G7["NES"].reindex(terms).dropna(), G7["padj"].reindex(terms).dropna()
    nt = len(terms)
    return dict(NES24=n24.mean(), NES7=n7.mean(),
                frac24=(p24 < FDR).sum() / nt, frac7=(p7 < FDR).sum() / nt)


rows = [{"name": nm, **agg(CL[CL["cluster"].isin(ids)])} for nm, ids in MERGE.items()]
rows += [{"name": nm, **agg(CL[CL["cluster"] == cid])} for cid, nm in NAME.items()]
d = pd.DataFrame(rows).sort_values("NES7", ascending=False).reset_index(drop=True)
n = len(d)

cmap = mpl.colormaps[S.CMAP_DIVERGING]
norm = mpl.colors.Normalize(vmin=-2.0, vmax=2.0)
COLS = [("NES24", "frac24", "24 HPT", "(early)"), ("NES7", "frac7", "7 DPT", "(late)")]


def panel_7A():
    FS = (7.4, 5.2)
    fig = plt.figure(figsize=FS)
    # 2026-09-06 (round 2, Stewart): the "number = mean member NES ..." footnote is dropped
    # (that detail belongs in the legend / Methods); the grid takes the freed height.
    ax = fig.add_axes([0.40, 0.19, 0.43, 0.59])    # cell grid; labels spill left, cbar right
    ax.set_xlim(0, 2); ax.set_ylim(0, n); ax.axis("off")

    for i, (_, r) in enumerate(d.iterrows()):
        y = n - 1 - i
        ax.text(-0.06, y + 0.5, r["name"], ha="right", va="center",
                fontsize=S.FS_TICK, color=S.COL_TEXT, clip_on=False, linespacing=1.05)
        for j, (col, fcol, _, _) in enumerate(COLS):
            val = r[col]; sig = r[fcol] >= FRAC
            base = cmap(norm(val))
            if sig:
                ax.add_patch(mpatches.Rectangle((j + 0.03, y + 0.06), 0.94, 0.88,
                             facecolor=base, edgecolor="#111111", linewidth=1.8, zorder=2))
                txt = f"{val:+.1f} *"
                tc = "white" if abs(val) > 1.15 else S.COL_TEXT
                ax.text(j + 0.5, y + 0.5, txt, ha="center", va="center",
                        fontsize=S.FS_TICK, fontweight="bold", color=tc, zorder=3)
            else:
                ax.add_patch(mpatches.Rectangle((j + 0.03, y + 0.06), 0.94, 0.88,
                             facecolor=base, alpha=0.28, edgecolor="#9E9E9E",
                             linewidth=0.8, hatch="////", zorder=2))
                ax.text(j + 0.5, y + 0.5, f"{val:+.1f}", ha="center", va="center",
                        fontsize=S.FS_TICK, color="#777777", zorder=3)

    # column headers (two-line, bold; wide vertical gap because 1 data unit is a thin slice)
    for j, (_, _, hdr, sub) in enumerate(COLS):
        ax.text(j + 0.5, n + 0.62, hdr, ha="center", va="bottom",
                fontsize=S.FS_AXIS, fontweight="bold", color=S.COL_TEXT, clip_on=False)
        ax.text(j + 0.5, n + 0.20, sub, ha="center", va="bottom",
                fontsize=S.FS_ANNOT, color="#666666", clip_on=False)

    # title (centered over the panel, house style)
    fig.text(0.5, 0.945, "Brain process state, early to late",
             ha="center", fontsize=S.FS_TITLE, color=S.COL_TEXT)

    # colorbar
    cax = fig.add_axes([0.855, 0.30, 0.022, 0.39])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_ticks([-2, -1, 0, 1, 2])
    cb.set_label("NES (red = up in MDT, blue = down)", fontsize=S.FS_LEGEND)
    cb.ax.tick_params(labelsize=S.FS_LEGEND)

    # Encoding legend. "flagged", not "significant": no process-level test is performed,
    # the star is a fraction-of-members rule. Spelling out that the number and the star
    # measure DIFFERENT things is what stops "+0.2 *" beside "+0.9" reading as a
    # contradiction (it is a small mean NES whose members are mostly individually sig,
    # sitting next to a larger mean NES whose members mostly are not).
    sig_h = mpatches.Patch(facecolor="#C24A38", edgecolor="#111111", linewidth=1.6,
                           label="flagged  (≥25% of member GO terms at FDR < 0.05)")
    ns_h = mpatches.Patch(facecolor="#C24A38", alpha=0.28, edgecolor="#9E9E9E",
                          linewidth=0.8, hatch="////", label="not flagged")
    fig.legend(handles=[sig_h, ns_h], loc="lower center", ncol=2, frameon=False,
               fontsize=S.FS_LEGEND, bbox_to_anchor=(0.50, 0.06), handlelength=1.6,
               columnspacing=1.4)

    for d_ in (OUT, FPNG):
        fig.savefig(d_ / "fig7a_process_state.png", dpi=S.DEFAULT_DPI,
                    bbox_inches=Bbox.from_bounds(0, 0, *FS))
    plt.close(fig)
    print("wrote fig7a_process_state.png")
    print("\n=== VERIFY 7A (NES7-sorted) ===")
    for _, r in d.iterrows():
        s24 = "*" if r["frac24"] >= FRAC else " "
        s7 = "*" if r["frac7"] >= FRAC else " "
        print(f"{r['name']:<44} 24HPT {r.NES24:+.1f}{s24}  7DPT {r.NES7:+.1f}{s7}")


# ===================================================================== baseline anchors
def _hallmark_vs_baseline():
    """Hallmark rank-shift GSEA for both arms vs the GSE73419 uninjured baseline."""
    def h(fn):
        d = pd.read_csv(GSE / fn, sep="\t")
        return d[d["collection"] == "Hallmark"].set_index("pathway")
    c, nd = h("gsea_rank_shift_CT.tsv"), h("gsea_rank_shift_ND.tsv")
    m = pd.DataFrame({"NES_ND": nd["NES"], "padj_ND": nd["padj"],
                      "NES_MDT": c["NES"], "padj_MDT": c["padj"]}).dropna()
    m["sig_either"] = (m["padj_ND"] < FDR) | (m["padj_MDT"] < FDR)
    return m


SHORT_7B = {"Epithelial mesenchymal transition": "Epithelial-mesenchymal (EMT)"}


def panel_7B(figsize=(4.6, 5.2)):
    """Dumbbell: each Hallmark's NES vs uninjured baseline in the ND and MDT conditions.

    The dumbbells come out near-zero-length, which IS the result: amputation moves the
    brain away from uninjured, and treatment does not redirect that at the process
    level, it mostly changes magnitude. MTORC1 is the one that separates.
    2026-09-06: figsize 5.6 -> 4.6 (the deck placed it at 4.36 in with a 13% left crop
    that clipped the long labels); label gutter widened to 0.57 so nothing clips.
    """
    m = _hallmark_vs_baseline()
    d = m[m["sig_either"]].sort_values("NES_MDT").reset_index()
    y = np.arange(len(d))

    fig = plt.figure(figsize=figsize)
    # 2026-09-06 (round 2, Stewart): "open symbol = n.s." footnote dropped (legend text);
    # the axes take the freed height.
    ax = fig.add_axes([0.57, 0.16, 0.40, 0.755])
    ax.axvline(0, color="#555555", lw=0.6, zorder=1)
    for yi, r in zip(y, d.itertuples()):
        ax.plot([r.NES_ND, r.NES_MDT], [yi, yi], color="#BBBBBB", lw=1.2, zorder=2)
    # ND is drawn LARGER and BEHIND, MDT smaller and in front. Several terms have the two
    # arms at effectively the same NES (TNFA -2.05/-2.05), and at equal size the blue
    # simply vanished under the coral, reading as a missing value rather than as exact
    # agreement. Concentric = the arms agree; separated = they differ.
    for col, pcol, lab, sz, zo in (("NES_ND", "padj_ND", "ND", 78, 3),
                                   ("NES_MDT", "padj_MDT", "MDT", 26, 4)):
        colr = S.COL_ND if lab == "ND" else S.COL_MDT
        sig = d[pcol] < FDR
        ax.scatter(d.loc[sig, col], y[sig.values], s=sz, facecolor=colr,
                   edgecolor="black", linewidths=0.6, zorder=zo)
        ax.scatter(d.loc[~sig, col], y[(~sig).values], s=sz, facecolor="white",
                   edgecolor=colr, linewidths=1.1, zorder=zo)
    ax.set_yticks(y)
    labs = [p.replace("HALLMARK_", "").replace("_", " ").capitalize() for p in d["pathway"]]
    ax.set_yticklabels([SHORT_7B.get(l, l) for l in labs], fontsize=S.FS_TICK - 1)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.7, len(d) - 0.3)
    ax.set_xlabel("NES vs baseline", labelpad=2)
    ax.grid(True, axis="x", alpha=0.15)
    fig.suptitle("Hallmark shift from uninjured baseline", fontsize=S.FS_TITLE,
                 color=S.COL_TEXT, x=0.5, y=0.965)
    # Legend + footnote sit BELOW the axes (2026-09-06): at the review-round font size an
    # in-axes legend covered the bottom rows on this narrower panel.
    handles = [Line2D([0], [0], marker="o", ls="none", mfc=S.COL_ND, mec="black",
                      mew=0.6, ms=8.5, label="ND vs baseline"),
               Line2D([0], [0], marker="o", ls="none", mfc=S.COL_MDT, mec="black",
                      mew=0.6, ms=5, label="MDT vs baseline")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=S.FS_LEGEND, handletextpad=0.3, columnspacing=1.2,
               bbox_to_anchor=(0.55, 0.012))
    for d_ in (OUT, FPNG):
        fig.savefig(d_ / "fig7B_hallmark_vs_baseline.png", dpi=S.DEFAULT_DPI,
                    bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    plt.close(fig)
    print("wrote fig7B_hallmark_vs_baseline.png")
    print(f"   {len(m)} Hallmarks | {int(m.sig_either.sum())} sig in either arm | "
          f"{int(((m.padj_ND < FDR) & (m.padj_MDT < FDR)).sum())} in both | "
          f"{int((m.NES_ND * m.NES_MDT < 0).sum())} sign-discordant")


def panel_7C(figsize=(3.9, 3.9)):
    """Genome-wide concordance of the two arms' shifts from baseline."""
    g = pd.read_csv(GSE / "rank_shift_genome_wide.tsv", sep="\t")
    n_rows = len(g)
    s = g.dropna(subset=["rank_shift_ND", "rank_shift_CT"]).copy()
    r = float(np.corrcoef(s["rank_shift_ND"], s["rank_shift_CT"])[0, 1])
    nd_mv, ct_mv = s["rank_shift_ND"].abs() > MOVER, s["rank_shift_CT"].abs() > MOVER
    shared, union = int((nd_mv & ct_mv).sum()), int((nd_mv | ct_mv).sum())
    s["mv"] = nd_mv | ct_mv

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.215, 0.145, 0.745, 0.735])
    lim = 0.85
    ax.plot([-lim, lim], [-lim, lim], color="#999999", lw=0.7, ls="--", zorder=2)
    ax.scatter(s.loc[~s.mv, "rank_shift_ND"], s.loc[~s.mv, "rank_shift_CT"], s=1.4,
               color="#D2D2D2", linewidths=0, rasterized=True, zorder=1)
    ax.scatter(s.loc[s.mv, "rank_shift_ND"], s.loc[s.mv, "rank_shift_CT"], s=1.6,
               color="#8FB4D9", linewidths=0, rasterized=True, zorder=1)
    gd = s[s["gene_symbol"] == "gdf15.L"]
    gx, gy = float(gd["rank_shift_ND"].iloc[0]), float(gd["rank_shift_CT"].iloc[0])
    ax.scatter([gx], [gy], s=48, facecolor=S.COL_MDT, edgecolor="black",
               linewidths=0.8, zorder=4)
    ax.annotate("gdf15.L", xy=(gx, gy), xytext=(-0.38, 0.55), fontsize=S.FS_ANNOT,
                fontstyle="italic",
                arrowprops=dict(arrowstyle="-", color="#666666", lw=0.6), zorder=5)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("ND percentile shift vs baseline", labelpad=2)
    ax.set_ylabel("MDT percentile shift vs baseline", labelpad=2)
    ax.grid(True, alpha=0.15)
    # stats sit lower-right (2026-09-06, round 2): the diagonal cloud leaves that corner empty,
    # whereas top-left collided with the gdf15.L callout at the review-round font size
    ax.annotate(f"r = {r:.3f}  (n = {len(s):,})\n{shared:,} of {union:,} movers shared",
                xy=(0.96, 0.04), xycoords="axes fraction", ha="right", va="bottom",
                fontsize=S.FS_LEGEND, color="#333333")
    # short title: the axis labels already name both arms and the baseline
    fig.suptitle("Genome-wide shift vs baseline",
                 fontsize=S.FS_TITLE, color=S.COL_TEXT, x=0.59, y=0.965)
    for d_ in (OUT, FPNG):
        fig.savefig(d_ / "fig7C_baseline_concordance.png", dpi=S.DEFAULT_DPI,
                    bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    plt.close(fig)
    print("wrote fig7C_baseline_concordance.png")
    print(f"   r={r:.4f} over {len(s):,} scored genes ({n_rows - len(s):,} of {n_rows:,} "
          f"dropped for a missing shift) | movers |shift|>{MOVER}: shared {shared:,} / union {union:,}")
    gg = s.assign(k=(s["rank_shift_CT"] - s["rank_shift_ND"]).abs()) \
          .sort_values("k", ascending=False).reset_index(drop=True)
    print(f"   gdf15.L off-diagonal rank {int(gg.index[gg.gene_symbol == 'gdf15.L'][0]) + 1} "
          f"of {len(gg):,} by |MDT-ND|")


def _pctile_recipe(sub, cols, sym):
    """Expression percentile of `sym` within one dataset, per the REGISTERED pipeline.

    Recipe reverse-engineered and validated 2026-07-17: restrict to the 27,478 gene_ids
    in rank_shift_genome_wide.tsv, mean the dataset's columns, drop genes undetected in
    THAT dataset (mean == 0), then rank(pct=True). Reproduces the shipped columns for
    gdf15.L EXACTLY -- gse_pctile 0.0734, mur_nd_pctile 0.1359, mur_ct_pctile 0.3701 --
    so the per-timepoint points below are on the same footing as 7C, not a new method.
    """
    d = sub.assign(m=sub[cols].mean(axis=1))
    d = d[d["m"] > 0]
    p = d["m"].rank(pct=True)
    i = d.index[d["gene_symbol"] == sym]
    return float(p.loc[i[0]]) if len(i) else np.nan


def panel_7D(figsize=(5.0, 3.9)):
    """gdf15 percentile by arm and timepoint, against an uninjured-baseline reference.

    REBUILT 2026-07-17 (Stewart's two objections, both correct):

    1. The old panel's ND and MDT points were POOLED over 24 HPT + 7 DPT
       (mur_{nd,ct}_mean average all six samples per arm), so they were at no timepoint
       at all. That matters more here than anywhere: gdf15.L is the most time-dependent
       gene in the brain (top 7 DPT DEG, trt x time beta +1.41, Fig 5A). Pooling blends
       MDT 24 HPT 0.159 with MDT 7 DPT 0.516 into a meaningless 0.370 and averages away
       the entire finding. Split by timepoint, the story matches Fig 5A exactly.
    2. The old panel drew a LINE baseline -> ND -> MDT. Those three are not a series:
       an external dataset and two arms pooled over time. The line manufactured a
       trajectory. Lines now connect ONLY 24 HPT -> 7 DPT within an arm, which is a real
       time course; the baseline is a horizontal REFERENCE LEVEL (it has no timepoint),
       not a waypoint.

    gdf15.S is the internal control: same locus, other homeolog, flat everywhere.
    """
    t = pd.read_csv(GSE / "tpm_brain_combined.tsv", sep="\t")
    keep = set(pd.read_csv(GSE / "rank_shift_genome_wide.tsv", sep="\t")["gene_id"])
    sub = t[t["gene_id"].isin(keep)].copy()
    BASE_C = ["GSE_Taira_brain", "GSE_Ueno_brain"]
    ARM = {"ND": "ND", "MDT": "CT"}          # display label -> column token
    TP = [("24HPA", "24 HPT"), ("7DPA", "7 DPT")]
    x = np.arange(2)

    fig = plt.figure(figsize=figsize)
    axes = [fig.add_axes([0.135 + i * 0.455, 0.150, 0.375, 0.655]) for i in range(2)]
    vals = {}
    for ax, gene in zip(axes, ("gdf15.L", "gdf15.S")):
        base = _pctile_recipe(sub, BASE_C, gene)
        ax.axhline(base, color="#8A8A8A", lw=1.1, ls="--", zorder=2)
        # Caption goes below the baseline line, on whichever side is clear: gdf15.L keeps
        # rising so its right side is busy (anchor left); gdf15.S the arms converge on the
        # right (anchor left there too, where 24 HPT still has headroom under the line).
        # centred between the two timepoints, under the line: that strip is empty for both
        # homeologs, whereas either end collides with a point label at the larger font size
        # 2026-09-06: the baseline is named in each panel's legend instead (an in-plot caption
        # collided with a point label wherever it sat at the larger font size)
        ax.lines[-1].set_label(f"uninjured baseline ({base:.2f})")
        series = {}
        for lab, tok in ARM.items():
            v = [_pctile_recipe(sub, [f"MUR_{tok}_T{i}_{tp}_BRAIN" for i in (1, 2, 3)], gene)
                 for tp, _ in TP]
            series[lab] = v
            vals[(gene, lab)] = v
        for lab, v in series.items():
            colr = S.COL_ND if lab == "ND" else S.COL_MDT
            ax.plot(x, v, color=colr, lw=1.7, marker="o", ms=7, mec="black", mew=0.6,
                    label=lab, zorder=3)
            other = series["MDT" if lab == "ND" else "ND"]
            for xi, vi, oi in zip(x, v, other):
                # label above when this series is the higher one at this x, else below,
                # so the two arms' values never collide (they nearly touch on gdf15.S)
                dy = 9 if vi >= oi else -15
                ax.annotate(f"{vi:.2f}", xy=(xi, vi), xytext=(0, dy),
                            textcoords="offset points", ha="center",
                            fontsize=S.FS_ANNOT - 0.5, color=colr)
        ax.set_xticks(x); ax.set_xticklabels([n for _, n in TP], fontsize=S.FS_TICK)
        ax.set_xlim(-0.4, 1.35); ax.set_ylim(-0.03, 0.66)
        ax.set_title(gene, fontsize=S.FS_AXIS, fontstyle="italic", color="#333333", pad=3)
        ax.grid(True, axis="y", alpha=0.15)
        ax.legend(loc="upper left", frameon=False, fontsize=S.FS_LEGEND - 1.5,
                  handlelength=1.6, borderaxespad=0.3, labelspacing=0.3)
    axes[0].set_ylabel("Expression percentile", labelpad=2)
    axes[1].set_yticklabels([])
    fig.suptitle("gdf15 percentile vs uninjured baseline",
                 fontsize=S.FS_TITLE, color=S.COL_TEXT, y=0.955)
    for d_ in (OUT, FPNG):
        fig.savefig(d_ / "fig7D_gdf15_baseline_anchor.png", dpi=S.DEFAULT_DPI,
                    bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    plt.close(fig)
    print("wrote fig7D_gdf15_baseline_anchor.png")
    for (gene, lab), v in vals.items():
        print(f"   {gene:8s} {lab:3s}  24 HPT {v[0]:.3f} -> 7 DPT {v[1]:.3f}")


def panel_7E(figsize=(3.5, 3.9)):
    """Walker etv1+ signature by arm and timepoint, against the uninjured baseline.

    ADDED 2026-08-07. This is the sibling of 7D: one feature scored per arm per
    timepoint against the GSE73419 anchor. 7D shows a gene where the arms separate
    (gdf15.L); this shows a signature where they barely do, but where BOTH sit far
    above uninjured.

    Why it lives in Fig 7 and not Fig 9: the evidence type is cross-dataset baseline
    anchoring, which is what all of Fig 7 is, and it carries the Fig 7 caveats (n = 2
    baseline, different clutch and library chemistry, rank-based, NO between-dataset
    statistics). Fig 9c/9d are within-experiment statistics. Keeping the two apart
    keeps one evidence standard per figure.

    FOOTING WARNING: scored on the combined-TPM table, like 7B/7C/7D. Fig 9C scores the
    same signature on voom E for the Murugan samples only, so absolute UCell values here
    are NOT comparable to 9C's. Only the arm/timepoint pattern and the distance from
    baseline are being read.

    Signature is the 13 Walker etv1+/glutamatergic genes used in 9C.
    """
    from scipy.stats import rankdata
    SIG = ["etv1.S", "slc17a7.L", "slc17a7.S", "cdh13.L", "cdh13.S", "bcl11b.L",
           "bcl11b.S", "sox5.L", "sox5.S", "neto1.L", "neto1.S", "stk39.L", "stk39.S"]
    BASE_C = ["GSE_Taira_brain", "GSE_Ueno_brain"]
    ARM = {"ND": "ND", "MDT": "CT"}
    TP = [("24HPA", "24 HPT"), ("7DPA", "7 DPT")]

    t = pd.read_csv(GSE / "tpm_brain_combined.tsv", sep="\t")
    rows = np.array([i for i, s in enumerate(t["gene_symbol"]) if s in SIG])
    assert len(rows) == len(SIG), f"signature genes missing: {len(rows)} of {len(SIG)}"

    def ucell(cols, max_rank=1500):
        """UCell (Andreatta & Carmona 2021), identical recipe to 9C: rank descending
        within a sample, cap at maxRank+1, AUC = 1 - (sum_rank - n(n+1)/2)/(n*maxRank)."""
        v = t[cols].mean(axis=1).to_numpy()
        n = len(rows)
        r = rankdata(-v, method="average")[rows].astype(float)
        r[r > max_rank] = max_rank + 1
        return 1.0 - (r.sum() - n * (n + 1) / 2.0) / (n * max_rank)

    base = ucell(BASE_C)
    x = np.arange(len(TP))
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.225, 0.135, 0.745, 0.75])
    ax.axhline(base, color="#8A8A8A", lw=1.1, ls="--", zorder=2,
               label=f"uninjured baseline ({base:.3f})")   # named in the legend (2026-09-06)

    series = {}
    for lab, tok in ARM.items():
        series[lab] = [ucell([f"MUR_{tok}_T{i}_{tp}_BRAIN" for i in (1, 2, 3)])
                       for tp, _ in TP]
    for lab, v in series.items():
        colr = S.COL_ND if lab == "ND" else S.COL_MDT
        ax.plot(x, v, color=colr, lw=1.7, marker="o", ms=7, mec="black", mew=0.6,
                label=lab, zorder=3)
        other = series["MDT" if lab == "ND" else "ND"]
        for xi, vi, oi in zip(x, v, other):
            ax.annotate(f"{vi:.3f}", xy=(xi, vi), xytext=(0, 9 if vi >= oi else -15),
                        textcoords="offset points", ha="center",
                        fontsize=S.FS_ANNOT - 0.5, color=colr)

    ax.set_xticks(x); ax.set_xticklabels([n for _, n in TP], fontsize=S.FS_TICK)
    ax.set_xlim(-0.38, 1.28)
    lo = base - 0.022
    ax.set_ylim(lo, max(max(v) for v in series.values()) + 0.016)
    ax.set_ylabel("UCell score", labelpad=3)
    ax.grid(True, axis="y", alpha=0.15)
    ax.legend(loc="lower left", frameon=False, fontsize=S.FS_LEGEND - 1.5, borderaxespad=0.4,
              handlelength=1.6, labelspacing=0.3)
    # title kept short: the panel is slim, and the dashed line is captioned in-plot
    fig.suptitle("etv1+ signature vs baseline",
                 fontsize=S.FS_TITLE - 1, color=S.COL_TEXT, x=0.6, y=0.965)
    for d_ in (OUT, FPNG):
        fig.savefig(d_ / "fig7E_walker_etv1_baseline_anchor.png", dpi=S.DEFAULT_DPI,
                    bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    plt.close(fig)
    print("wrote fig7E_walker_etv1_baseline_anchor.png")
    print(f"   uninjured baseline {base:.4f} ({len(SIG)} genes, n = {len(BASE_C)})")
    for lab, v in series.items():
        print(f"   {lab:3s}  24 HPT {v[0]:.4f} -> 7 DPT {v[1]:.4f}   "
              f"({v[0]/base:.2f}x -> {v[1]/base:.2f}x baseline)")


if __name__ == "__main__":
    panel_7A()
    print()
    panel_7B()
    panel_7C()
    panel_7D()
    panel_7E()
