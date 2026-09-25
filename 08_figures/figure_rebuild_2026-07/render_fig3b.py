"""Fig 3B - vwf.L / col4a1.L / col4a2.L brain expression trajectories (mean +/- s.e.m.).

Split out of render_voom.py on 2026-09-06 so 3B can be re-rendered on its own: running
render_voom.py would also overwrite the LOCKED 2F/3C/4D heatmap mirrors with its older,
drifting clustering. Same data and look as the deck panel (means + s.e.m. per condition,
no per-animal points, MDT coral / ND blue); reads the canonical tables/ copies, and
saves at a FIXED figsize (3.78 x 2.12 in, the deck slot) instead of bbox='tight'.
"""
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox
import _style as S

BASE = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent / "panels_png"; OUT.mkdir(exist_ok=True)
FPNG = BASE / "figures_png"; FPNG.mkdir(exist_ok=True)

E = pd.read_csv(BASE / "tables" / "brain" / "voom_E_BRAIN_primary.tsv", sep="\t", index_col=0)
md = pd.read_csv(BASE / "tables" / "sample_metadata.tsv", sep="\t")
md = md[md.tissue == "BRAIN"].copy()
md["trt"] = md.treatment.map({"CT": "MDT", "ND": "ND"})
md = md.set_index("sample")


def trajectories(genes, title, name, ncol=3, figsize=(3.78, 2.12)):
    genes = [g for g in genes if g in E.index]
    nrow = int(np.ceil(len(genes) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize, squeeze=False)
    tps = ["24HPA", "7DPA"]; xlab = ["24 HPT", "7 DPT"]
    for i, g in enumerate(genes):
        ax = axes[i // ncol][i % ncol]
        for trt, col in [("MDT", S.COL_MDT), ("ND", S.COL_ND)]:
            means, ses = [], []
            for tp in tps:
                cols = md[(md.trt == trt) & (md.time == tp)].index
                vals = E.loc[g, E.columns.intersection(cols)].astype(float)
                means.append(vals.mean()); ses.append(vals.std() / np.sqrt(len(vals)))
            ax.errorbar([0, 1], means, yerr=ses, color=col, marker="o", ms=3,
                        lw=1.2, capsize=2)
        ax.set_title(g, fontsize=S.FS_TICK, style="italic", pad=2)
        ax.set_xticks([0, 1]); ax.set_xticklabels(xlab, fontsize=S.FS_LEGEND - 1.5)
        ax.tick_params(axis="y", labelsize=S.FS_LEGEND)
        ax.set_xlim(-0.3, 1.3)
    for j in range(len(genes), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    axes[0][0].set_ylabel("norm. expr.", fontsize=S.FS_LEGEND)
    handles = [Line2D([0], [0], color=S.COL_MDT, marker="o", ms=3, label="MDT"),
               Line2D([0], [0], color=S.COL_ND, marker="o", ms=3, label="ND")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=S.FS_LEGEND, bbox_to_anchor=(0.5, 0.0), handlelength=1.6)
    fig.suptitle(title, fontsize=S.FS_TITLE, color=S.COL_TEXT, y=0.985)
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.34, top=0.76, wspace=0.55)
    for d in (OUT, FPNG):
        fig.savefig(d / f"{name}.png", dpi=S.DEFAULT_DPI, bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    plt.close(fig)
    print("wrote", name)


if __name__ == "__main__":
    trajectories(["vwf.L", "col4a1.L", "col4a2.L"], "ECM gene trajectories", "fig3B_ecm_trajectories")
