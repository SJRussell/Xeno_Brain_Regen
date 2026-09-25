"""Fig 8 - Paired brain-blastema MOFA (Factor 3 = treatment x time cross-tissue axis).

  8A  MOFA factor variance explained per view + factor x design-term associations (two sub-panels)
  8B  Factor 3 sample-score trajectory, ND vs CT crossing (treatment x time FDR = 0.029)
  8C  Factor 3 top brain-view loadings
  8D  Factor 3 top blastema-view loadings
  8E  [placeholder] Nirosha working-model schematic

House style (Arial, FS_*), BUT Fig 8 KEEPS "CT" not "MDT" (documented exception:
the MOFA design/source uses CT). All values reproduce from tables/mofa/* and the
recorded mofa_* metas (verified 2026-07-14).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.transforms import Bbox
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
import _style as S

BASE = Path(__file__).resolve().parent.parent
MOFA = BASE / "tables" / "mofa"
OUT = Path(__file__).resolve().parent / "panels_png"; OUT.mkdir(exist_ok=True)
FPNG = BASE / "figures_png"; FPNG.mkdir(exist_ok=True)
TIME = {"24HPA": "24 HPT", "7DPA": "7 DPT"}
COL_CT, COL_ND = S.COL_MDT, S.COL_ND          # CT drawn coral, ND blue (Fig 8 keeps "CT")


def _save(fig, name, figsize):
    for d in (OUT, FPNG):
        fig.savefig(d / f"{name}.png", dpi=S.DEFAULT_DPI, bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    plt.close(fig)
    print("wrote", name)


def _heat(ax, M, rows, cols, cmap, norm, fmt="{:.1f}", stars=None, dark=1.4):
    ax.imshow(M, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, fontsize=S.FS_TICK - 1)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=S.FS_TICK - 1)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            star = "" if stars is None else ("*" if stars[i, j] else "")
            tc = "white" if norm(v) > 0.62 else S.COL_TEXT
            ax.text(j, i, fmt.format(v) + star, ha="center", va="center",
                    fontsize=S.FS_ANNOT, color=tc)


# ---------------------------------------------------------------- 8A variance + associations
def panel_8A():
    var = pd.read_csv(MOFA / "factor_variance_explained.tsv", sep="\t").set_index("factor")
    order = var["BRAIN"].sort_values(ascending=False).index.tolist()   # F1,F2,F4,F3,F5,F6
    var = var.reindex(order)
    Mv = var[["BRAIN", "BLAST"]].to_numpy()

    A = pd.read_csv(MOFA / "factor_treatment_associations.tsv", sep="\t")
    TERMS = [("time7DPA", "7 DPT"), ("treatmentCT", "MDT"),
             ("treatmentCT:time7DPA", "MDT:7DPT"), ("repT2", "repT2"), ("repT3", "repT3")]
    p = A.pivot(index="factor", columns="term", values="p").reindex(order)
    fdr = A.pivot(index="factor", columns="term", values="fdr").reindex(order)
    Ma = np.column_stack([-np.log10(p[t].to_numpy()) for t, _ in TERMS])
    Sa = np.column_stack([(fdr[t].to_numpy() < 0.05) for t, _ in TERMS])

    warm = mpl.colormaps["YlOrRd"]
    nV = mpl.colors.Normalize(0, Mv.max())
    nA = mpl.colors.Normalize(0, np.nanmax(Ma))
    # 2026-09-06 (round 2, Stewart: "the two heatmaps are not aligned, too much whitespace"):
    # one tight image placed UNCROPPED in the deck (the two hand-crops of the old wide image
    # were scaled differently and drifted apart). Rows are the same six factors in the same
    # order, so the associations panel drops its duplicate row labels and sits right after
    # the variance colorbar; the variance panel's "MOFA factor" label now has its margin.
    FS = (6.5, 3.4)
    fig = plt.figure(figsize=FS)
    axV = fig.add_axes([0.155, 0.27, 0.245, 0.55])
    axA = fig.add_axes([0.535, 0.27, 0.34, 0.55])

    # -- variance (numbers) --
    axV.imshow(Mv, cmap=warm, norm=nV, aspect="auto")
    axV.set_xticks([0, 1]); axV.set_xticklabels(["Brain", "Blast"], fontsize=S.FS_TICK - 1)
    axV.set_yticks(range(len(order))); axV.set_yticklabels(order, fontsize=S.FS_TICK - 1)
    axV.set_xlabel("Tissue", fontsize=S.FS_AXIS); axV.set_ylabel("MOFA factor", fontsize=S.FS_AXIS)
    axV.tick_params(length=0)
    for s in axV.spines.values():
        s.set_visible(False)
    for i in range(Mv.shape[0]):
        for j in range(2):
            v = Mv[i, j]
            axV.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=S.FS_ANNOT,
                     color="white" if nV(v) > 0.6 else S.COL_TEXT)
    axV.set_title("Variance explained (%)", fontsize=S.FS_AXIS + 0.5,
                  fontweight="bold", color=S.COL_TEXT, pad=6)
    cbV = fig.colorbar(ScalarMappable(norm=nV, cmap=warm), cax=fig.add_axes([0.412, 0.30, 0.015, 0.48]))
    cbV.set_label("% variance", fontsize=S.FS_LEGEND); cbV.ax.tick_params(labelsize=S.FS_LEGEND)

    # -- associations (stars only, no numbers; rows shared with the variance panel) --
    axA.imshow(Ma, cmap=warm, norm=nA, aspect="auto")
    axA.set_xticks(range(len(TERMS))); axA.set_xticklabels([lab for _, lab in TERMS],
                                                           fontsize=S.FS_TICK - 2,
                                                           rotation=35, ha="right",
                                                           rotation_mode="anchor")
    axA.set_yticks(range(len(order))); axA.set_yticklabels([])
    axA.set_xlabel("Covariate", fontsize=S.FS_AXIS)
    axA.tick_params(length=0)
    for s in axA.spines.values():
        s.set_visible(False)
    for i in range(Ma.shape[0]):
        for j in range(Ma.shape[1]):
            if Sa[i, j]:
                axA.text(j, i, "*", ha="center", va="center", fontsize=S.FS_TITLE + 3,
                         color="white", fontweight="bold")
    axA.set_title("Factor-covariate associations", fontsize=S.FS_AXIS + 0.5,
                  fontweight="bold", color=S.COL_TEXT, pad=6)
    cbA = fig.colorbar(ScalarMappable(norm=nA, cmap=warm), cax=fig.add_axes([0.887, 0.30, 0.015, 0.48]))
    cbA.set_label(r"$-\log_{10}(p)$", fontsize=S.FS_LEGEND); cbA.ax.tick_params(labelsize=S.FS_LEGEND)
    _save(fig, "fig8A_mofa_variance_assoc", FS)


# ---------------------------------------------------------------- 8B Factor 3 trajectory
def panel_8B():
    sc = pd.read_csv(MOFA / "factor_scores.tsv", sep="\t")
    times = ["24HPA", "7DPA"]
    dx = {"ND": -0.06, "CT": 0.06}
    FS = (3.8, 3.4)
    fig = plt.figure(figsize=FS)
    ax = fig.add_axes([0.23, 0.14, 0.745, 0.72])
    handles = {}
    for trt, col in (("ND", COL_ND), ("CT", COL_CT)):
        sub = sc[sc["treatment"] == trt]
        means, sems, xs = [], [], []
        for i, tp in enumerate(times):
            g = sub[sub["time"] == tp]["Factor3"].values
            xs.append(i + dx[trt])
            means.append(g.mean()); sems.append(g.std(ddof=1) / np.sqrt(len(g)))
            ax.scatter(np.full(len(g), i + dx[trt]) + np.linspace(-0.03, 0.03, len(g)), g,
                       s=16, color=col, edgecolors="white", linewidths=0.3, alpha=0.55, zorder=3)
        ln, = ax.plot(xs, means, color=col, lw=1.6, zorder=2, marker="o", ms=5)
        ax.errorbar(xs, means, yerr=sems, fmt="none", ecolor=col, elinewidth=1.0,
                    capsize=2.5, capthick=1.0, zorder=2)
        handles[trt] = ln
    ax.axhline(0, color="#999999", lw=0.6, ls="--")
    ax.set_xticks(range(len(times))); ax.set_xticklabels([TIME[t] for t in times])
    ax.set_xlim(-0.5, len(times) - 0.5)
    ax.set_ylabel("Factor 3 score", labelpad=2)
    ax.set_title("MOFA Factor 3 score", fontsize=S.FS_TITLE - 1,
                 color=S.COL_TEXT, pad=6)
    ax.text(0.03, 0.03, "treatment × time\nFDR = 0.029", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=S.FS_ANNOT, color="#333333", linespacing=1.3)
    ax.legend([handles["ND"], handles["CT"]], ["ND (amputation only)", "MDT (treated)"],
              loc="upper right", frameon=False, fontsize=S.FS_LEGEND, handlelength=1.6)
    _save(fig, "fig8B_mofa_factor3_trajectory", FS)


# ---------------------------------------------------------------- 8C/8D loadings lollipops
def _loadings(view, title, name):
    L = pd.read_csv(MOFA / "factor_loadings_top.tsv", sep="\t")
    sub = L[(L["view"] == view) & (L["factor"] == "Factor3")].copy()
    sub["abs"] = sub["weight"].abs()
    sub = sub.sort_values("abs", ascending=False).head(12).sort_values("weight")
    genes = [g.replace("_BRAIN", "").replace("_BLAST", "") for g in sub["gene"]]
    w = sub["weight"].to_numpy()
    y = np.arange(len(w))
    FS = (3.8, 3.4)
    fig = plt.figure(figsize=FS)
    ax = fig.add_axes([0.41, 0.15, 0.55, 0.705])
    for yi, wi in zip(y, w):
        c = COL_CT if wi > 0 else COL_ND
        ax.plot([0, wi], [yi, yi], color=c, lw=1.4, zorder=1)
        ax.scatter([wi], [yi], s=34, color=c, edgecolors="white", linewidths=0.4, zorder=2)
    ax.axvline(0, color="#999999", lw=0.6)
    ax.set_yticks(y); ax.set_yticklabels(genes, fontsize=S.FS_TICK - 1, fontstyle="italic")
    ax.set_ylim(-0.6, len(w) - 0.4)
    ax.set_xlabel("Factor 3 loading weight", labelpad=2)
    ax.set_title(title, fontsize=S.FS_TITLE - 1, color=S.COL_TEXT, pad=6)
    ax.grid(True, axis="x", alpha=0.18)
    _save(fig, name, FS)


# ---------------------------------------------------------------- 8E placeholder
def panel_8E():
    FS = (3.8, 3.4)
    fig = plt.figure(figsize=FS)
    ax = fig.add_axes([0.06, 0.08, 0.88, 0.84]); ax.axis("off")
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, transform=ax.transAxes, facecolor="#F4F4F4",
                               edgecolor="#BBBBBB", linewidth=1.0, linestyle=(0, (4, 3))))
    ax.text(0.5, 0.5, "Working model of coordinated\nbrain–blastema regenerative trajectory\n\n"
            "[schematic to come]", ha="center", va="center", fontsize=S.FS_AXIS,
            color="#888888", linespacing=1.5)
    _save(fig, "fig8E_working_model_placeholder", FS)


if __name__ == "__main__":
    panel_8A()
    panel_8B()
    _loadings("BRAIN", "Factor 3 loadings: brain", "fig8C_factor3_brain_loadings")
    _loadings("BLAST", "Factor 3 loadings: blastema", "fig8D_factor3_blast_loadings")
    panel_8E()
    print("\n=== VERIFY ===")
    var = pd.read_csv(MOFA / "factor_variance_explained.tsv", sep="\t").set_index("factor")
    print("F3 variance BRAIN/BLAST:", round(var.loc["Factor3", "BRAIN"], 1),
          "/", round(var.loc["Factor3", "BLAST"], 1))
    A = pd.read_csv(MOFA / "factor_treatment_associations.tsv", sep="\t")
    r = A[(A["factor"] == "Factor3") & (A["term"] == "treatmentCT:time7DPA")].iloc[0]
    print(f"F3 treatment×time: beta={r['estimate']:+.2f} p={r['p']:.4f} fdr={r['fdr']:.4f}")
