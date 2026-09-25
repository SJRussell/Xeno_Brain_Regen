"""Shared treatment x time trajectory helpers for the 2026-07 rebuild (Figs 5, 6).

One place for the mean-line + SEM + jittered-point interaction plot, its treatment legend,
the padj->star map, and the voom/DE readers. Both render_fig5.py and render_fig6.py import
from here so the panels stay identical (anti-drift).
"""
import numpy as np
import pandas as pd
import _style as S

TRT_ORDER = ["ND", "CT"]                 # draw ND first so MDT sits on top
TRT_DX    = {"ND": -0.09, "CT": 0.09}    # horizontal split so clouds don't overlap


def stars(padj):
    if padj is None or (isinstance(padj, float) and np.isnan(padj)): return "ns"
    if padj < 1e-3: return "***"
    if padj < 1e-2: return "**"
    if padj < 5e-2: return "*"
    return "ns"


def voom_row(fp, gene):
    """Per-sample voom E for one gene, as a Series indexed by sample."""
    d = pd.read_csv(fp, sep="\t", index_col=0)
    return d.loc[gene]


def de_padj(fp, gene):
    d = pd.read_csv(fp, sep="\t")
    gc = "gene" if "gene" in d.columns else d.columns[0]
    r = d[d[gc] == gene]
    return float(r["padj"].iloc[0]) if len(r) else np.nan


def draw_trajectory(ax, expr, meta_sub, times, star_by_time, ylabel_text=None):
    """One trajectory sub-plot: mean line + SEM + jittered points per treatment.
    `expr` maps sample -> value; `meta_sub` has sample/treatment/time for this tissue;
    `star_by_time` maps a data-time (e.g. '24HPA') -> star string above that timepoint;
    `ylabel_text` is the y-axis label (None = no label, e.g. right-column facets)."""
    m = meta_sub.copy()
    m["val"] = m["sample"].map(expr)
    m = m.dropna(subset=["val"])
    handles = {}
    ymax_at = {}
    for trt in TRT_ORDER:
        col = S.TREAT_COL[trt]
        sub = m[m["treatment"] == trt]
        means, sems, xs = [], [], []
        for i, tp in enumerate(times):
            g = sub[sub["time"] == tp]["val"].values
            xs.append(i + TRT_DX[trt])
            means.append(np.mean(g)); sems.append(np.std(g, ddof=1) / np.sqrt(len(g)))
            jit = np.linspace(-0.035, 0.035, len(g))
            ax.scatter(i + TRT_DX[trt] + jit, g, s=17, color=col, edgecolors="white",
                       linewidths=0.3, alpha=0.9, zorder=3)
            ymax_at[i] = max(ymax_at.get(i, -np.inf), (np.mean(g) + sems[-1]), g.max())
        ln, = ax.plot(xs, means, color=col, lw=1.5, zorder=2)
        ax.errorbar(xs, means, yerr=sems, fmt="none", ecolor=col, elinewidth=1.0,
                    capsize=2.5, capthick=1.0, zorder=2)
        handles[S.TREAT_LABEL[trt]] = ln

    ax.set_xticks(range(len(times)))
    ax.set_xticklabels([S.TIME_LABEL[t] for t in times])
    ax.set_xlim(-0.5, len(times) - 0.5)
    ax.set_xlabel("Timepoint", labelpad=2)
    if ylabel_text:
        ax.set_ylabel(ylabel_text, labelpad=2)
    ax.grid(True, axis="y", alpha=0.18)

    y0, y1 = ax.get_ylim(); pad = 0.06 * (y1 - y0)
    for i, tp in enumerate(times):
        s = star_by_time.get(tp, "ns")
        if s != "ns":
            ax.text(i, ymax_at[i] + pad, s, ha="center", va="bottom",
                    fontsize=S.FS_TITLE, color="#222222")
    ax.set_ylim(y0, y1 + 2.2 * pad)
    return handles


def treatment_legend(ax, handles, loc="upper left"):
    ax.legend([handles[k] for k in ("ND", "MDT")], ["ND", "MDT"], title="Treatment",
              loc=loc, frameon=False, fontsize=S.FS_LEGEND,
              title_fontsize=S.FS_LEGEND, handlelength=1.6)
