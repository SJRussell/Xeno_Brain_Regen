"""Fig 5 local panels - GDF15 brain-body coupling (NPJ Regen Med rebuild).

  5A  brain gdf15.L trajectory (voom E; MDT vs ND x time; *** at 7 DPT; trt:time beta)
  5B  blastema gdf15.L/.S trajectories (voom E; * on gdf15.L at 24 HPT)
  5C  brain ME1(turquoise) eigengene vs gdf15.L (lm + 95% CI band; Spearman)
  5D  cross-tissue gdf15.L coupling (blastema VST vs brain VST, paired per animal)

House-style consistency with Fig 2-4 (figure_rebuild _style): Arial, FS_TITLE=10 (centered,
matplotlib default), MDT/ND palette, timepoints shown as HPT/DPT, panel letters added later
as movable pptx textboxes. "CT" in the data is displayed as "MDT" everywhere.

All inputs are LOCAL (no Nibi dependency at render time):
  tables/brain/voom_E_BRAIN_primary.tsv                    -> 5A y, 5C x   (gdf15.L)
  tables/blastema/voom_E_BLAST_primary.tsv                 -> 5B          (gdf15.L/.S)
  tables/wgcna/brain_module_eigengenes_per_sample.tsv      -> 5C y        (ME1 = turquoise)
  tables/{brain,blastema}/persample_gdf15_{BRAIN,BLAST}.tsv-> 5D          (vst_gdf15.L)
  tables/brain/de_BRAIN_CT_vs_ND_{24HPA,7DPA}.tsv          -> 5A stars    (padj)
  tables/blastema/de_BLAST_CT_vs_ND_{11HPA,24HPA,7DPA}.tsv -> 5B stars    (padj)
  tables/sample_metadata.tsv                               -> treatment/time/rep

Annotations are LOCKED to the recorded v10 meta (tables/_meta/*_meta.yaml) and cross-checked
against a local recompute in the VERIFY block below:
  5A trt:time beta = +1.41, p = 0.0044   (brain_gdf15L_trajectory_v10_meta)
  5C lm: ME1 = +0.149*gdf15.L - 0.245 ; R2 = 0.470 ; Spearman rho = +0.70, p = 0.0114
  Stars from CT-vs-ND padj: *** < 1e-3, ** < 1e-2, * < 0.05.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.lines import Line2D
from scipy.stats import spearmanr, t as tdist
import _style as S

BASE = Path(__file__).resolve().parent.parent
TBL  = BASE / "tables"
META = pd.read_csv(TBL / "sample_metadata.tsv", sep="\t")
OUT  = Path(__file__).resolve().parent / "panels_png"; OUT.mkdir(exist_ok=True)
FPNG = BASE / "figures_png"; FPNG.mkdir(exist_ok=True)

# locked recorded-meta values (v10)
LOCK_5A_BETA, LOCK_5A_P = 1.406, 0.004373
LOCK_5C = dict(slope=0.1487, intercept=-0.245, r2=0.4701, rho=0.6993, p=0.01137)

# trajectory helpers are shared with render_fig6 (anti-drift); see _traj.py
from _traj import (TRT_ORDER, TRT_DX, stars, voom_row, de_padj,   # noqa: F401
                   draw_trajectory, treatment_legend)


# ------------------------------------------------------------------------ Panel 5A
def panel_5A():
    gene = "gdf15.L"
    expr = voom_row(TBL / "brain/voom_E_BRAIN_primary.tsv", gene).to_dict()
    msub = META[META["tissue"] == "BRAIN"][["sample", "treatment", "time", "rep"]]
    times = ["24HPA", "7DPA"]
    star_by = {t: stars(de_padj(TBL / f"brain/de_BRAIN_CT_vs_ND_{t}.tsv", gene)) for t in times}

    fig = plt.figure(figsize=(3.4, 3.3))
    ax = fig.add_axes([0.165, 0.135, 0.80, 0.75])
    h = draw_trajectory(ax, expr, msub, times, star_by, ylabel_text="Normalized expression")
    ax.set_title("Brain gdf15.L by treatment × time", fontsize=S.FS_TITLE,
                 color=S.COL_TEXT, pad=6)
    treatment_legend(ax, h, loc="upper left")
    ax.text(0.97, 0.03, f"β = +{LOCK_5A_BETA:.2f}, p = {LOCK_5A_P:.4f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=S.FS_ANNOT,
            color="#555555")
    _save(fig, "fig5A_brain_gdf15_trajectory", (3.4, 3.3))
    return star_by


# ------------------------------------------------------------------------ Panel 5B
def panel_5B():
    msub = META[META["tissue"] == "BLAST"][["sample", "treatment", "time", "rep"]]
    times = ["11HPA", "24HPA", "7DPA"]
    genes = ["gdf15.L", "gdf15.S"]
    vdf = pd.read_csv(TBL / "blastema/voom_E_BLAST_primary.tsv", sep="\t", index_col=0)
    star_maps = {}
    for g in genes:
        star_maps[g] = {t: stars(de_padj(TBL / f"blastema/de_BLAST_CT_vs_ND_{t}.tsv", g))
                        for t in times}

    fig = plt.figure(figsize=(4.9, 3.3))
    axes = [fig.add_axes([0.105, 0.15, 0.355, 0.67]),
            fig.add_axes([0.590, 0.15, 0.355, 0.67])]
    h = None
    for ax, g in zip(axes, genes):
        expr = vdf.loc[g].to_dict()
        h = draw_trajectory(ax, expr, msub, times, star_maps[g],
                            ylabel_text=("Normalized expression" if g == "gdf15.L" else None))
        ax.set_title(g, fontsize=S.FS_AXIS, fontstyle="italic", color="#333333", pad=3)
    treatment_legend(axes[1], h, loc="upper right")
    fig.suptitle("Blastema gdf15 (.L/.S) by treatment × time", fontsize=S.FS_TITLE,
                 color=S.COL_TEXT, y=0.98)
    _save(fig, "fig5B_blastema_gdf15_trajectories", (4.9, 3.3))
    return star_maps


# ------------------------------------------------------------------------ Panel 5C
def panel_5C():
    x = voom_row(TBL / "brain/voom_E_BRAIN_primary.tsv", "gdf15.L")
    me = pd.read_csv(TBL / "wgcna/brain_module_eigengenes_per_sample.tsv", sep="\t",
                     index_col=0)["ME1"]
    msub = META[META["tissue"] == "BRAIN"].set_index("sample")
    samples = [s for s in x.index if s in me.index]
    df = pd.DataFrame({"x": x[samples].values, "y": me[samples].values,
                       "trt": msub.loc[samples, "treatment"].values,
                       "time": msub.loc[samples, "time"].values})
    xv, yv = df["x"].values, df["y"].values
    a, b = np.polyfit(xv, yv, 1)
    rho, pr = spearmanr(xv, yv)
    n = len(xv); xbar = xv.mean(); Sxx = np.sum((xv - xbar) ** 2)
    s2 = np.sum((yv - (a * xv + b)) ** 2) / (n - 2)
    xs = np.linspace(xv.min(), xv.max(), 100)
    band = tdist.ppf(0.975, n - 2) * np.sqrt(s2 * (1 / n + (xs - xbar) ** 2 / Sxx))

    fig = plt.figure(figsize=(3.5, 3.3))
    ax = fig.add_axes([0.165, 0.135, 0.80, 0.75])
    ax.fill_between(xs, a * xs + b - band, a * xs + b + band, color="#CCCCCC", alpha=0.45, lw=0)
    ax.plot(xs, a * xs + b, color="#444444", lw=1.2, zorder=2)
    for trt in TRT_ORDER:
        for tp in ["24HPA", "7DPA"]:
            d = df[(df["trt"] == trt) & (df["time"] == tp)]
            ax.scatter(d["x"], d["y"], s=42, color=S.TREAT_COL[trt],
                       marker=S.TIME_SHAPE_BRAIN[tp], edgecolors="black", linewidths=0.5,
                       alpha=0.92, zorder=3)
    ax.set_xlabel("gdf15.L expression (voom log2)", labelpad=2)
    ax.set_ylabel("ME1 eigengene (turquoise module)", labelpad=2)
    ax.set_title("Brain ME1 eigengene vs gdf15.L", fontsize=S.FS_TITLE, color=S.COL_TEXT, pad=6)
    ax.grid(True, alpha=0.18)
    ann = (f"ME1 = +{LOCK_5C['slope']:.3f} × gdf15.L − {abs(LOCK_5C['intercept']):.3f}\n"
           f"R² = {LOCK_5C['r2']:.3f}\n"
           f"Spearman ρ = +{LOCK_5C['rho']:.2f}, p = {LOCK_5C['p']:.4f}")
    ax.text(0.04, 0.96, ann, transform=ax.transAxes, ha="left", va="top",
            fontsize=S.FS_ANNOT, color="#333333")
    leg1 = [Line2D([0], [0], marker="s", ls="none", mfc=S.COL_ND, mec="black", ms=6, label="ND"),
            Line2D([0], [0], marker="s", ls="none", mfc=S.COL_MDT, mec="black", ms=6, label="MDT")]
    leg2 = [Line2D([0], [0], marker="o", ls="none", mfc="#BBBBBB", mec="black", ms=6, label="24 HPT"),
            Line2D([0], [0], marker="^", ls="none", mfc="#BBBBBB", mec="black", ms=6, label="7 DPT")]
    l1 = ax.legend(handles=leg1, title="Treatment", loc="lower right", frameon=False,
                   fontsize=S.FS_LEGEND, title_fontsize=S.FS_LEGEND,
                   handletextpad=0.4, bbox_to_anchor=(1.0, 0.22))
    ax.add_artist(l1)
    ax.legend(handles=leg2, title="Timepoint", loc="lower right", frameon=False,
              fontsize=S.FS_LEGEND, title_fontsize=S.FS_LEGEND,
              handletextpad=0.4, bbox_to_anchor=(1.0, 0.0))
    _save(fig, "fig5C_me1_vs_gdf15", (3.5, 3.3))
    return dict(slope=a, intercept=b, r2=np.corrcoef(xv, yv)[0, 1] ** 2, rho=rho, p=pr, n=n)


# ------------------------------------------------------------------------ Panel 5D
def _load_persample(fp, col):
    d = pd.read_csv(fp, sep="\t")
    d[["trt", "animal", "tp"]] = d["sample"].str.split("_", expand=True).iloc[:, :3]
    return d.rename(columns={col: "g"})[["trt", "animal", "tp", "g"]]


def panel_5D():
    SHARED = ["24HPA", "7DPA"]
    b = _load_persample(TBL / "brain/persample_gdf15_BRAIN.tsv", "vst_gdf15.L").rename(columns={"g": "bg"})
    p = _load_persample(TBL / "blastema/persample_gdf15_BLAST.tsv", "vst_gdf15.L").rename(columns={"g": "pg"})
    m = b.merge(p, on=["trt", "animal", "tp"]).query("tp in @SHARED").reset_index(drop=True)
    x, y = m["pg"].values, m["bg"].values
    a, b0 = np.polyfit(x, y, 1)
    rho, pv = spearmanr(x, y)
    rc, _ = spearmanr(*[m[m["trt"] == "CT"][c] for c in ("pg", "bg")])
    r7, _ = spearmanr(*[m[m["tp"] == "7DPA"][c] for c in ("pg", "bg")])

    fig = plt.figure(figsize=(3.8, 3.3))
    ax = fig.add_axes([0.155, 0.135, 0.81, 0.77])
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, a * xs + b0, color="#888888", lw=1.2, ls="--", zorder=1)
    mark = {"24HPA": "o", "7DPA": "^"}
    for _, r in m.iterrows():
        ax.scatter(r["pg"], r["bg"], s=95, color=S.TREAT_COL[r["trt"]], marker=mark[r["tp"]],
                   edgecolors="black", linewidths=0.7, zorder=3)
    ax.set_xlabel("Blastema gdf15 (VST)", labelpad=2)
    ax.set_ylabel("Brain gdf15 (VST)", labelpad=2)
    # "Coupling" implies a mechanism a pooled correlation cannot show: rho=0.83 is
    # structured by treatment and time (within-ND 0.64, 24 HPT only 0.38; n=3 per
    # subgroup cell). "Association" is what the panel supports.
    ax.set_title("Paired cross-tissue gdf15.L association", fontsize=S.FS_TITLE,
                 color=S.COL_TEXT, pad=6)
    ax.grid(True, alpha=0.2)
    stat = (f"Spearman ρ = +{rho:.2f}, p = {pv:.3f}  (n = {len(m)})\n"
            f"within MDT: ρ = +{rc:.2f}   ·   7 DPT: ρ = +{r7:.2f}")
    ax.text(1.0, 0.05, stat, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=S.FS_ANNOT, color="#222222",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#CCCCCC", lw=0.8))
    handles = [
        Line2D([0], [0], marker="o", ls="none", mfc=S.COL_MDT, mec="black", ms=8, label="MDT (treated)"),
        Line2D([0], [0], marker="o", ls="none", mfc=S.COL_ND, mec="black", ms=8, label="ND (amputation only)"),
        Line2D([0], [0], marker="o", ls="none", mfc="#BBBBBB", mec="black", ms=8, label="24 HPT"),
        Line2D([0], [0], marker="^", ls="none", mfc="#BBBBBB", mec="black", ms=8, label="7 DPT"),
    ]
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=S.FS_LEGEND)
    _save(fig, "fig5D_crosstissue_coupling", (3.8, 3.3))
    return dict(rho=rho, within_mdt=rc, at_7dpt=r7, n=len(m))


def _save(fig, name, figsize):
    fp = OUT / f"{name}.png"
    fig.savefig(fp, dpi=S.DEFAULT_DPI, bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    fig.savefig(FPNG / f"{name}.png", dpi=S.DEFAULT_DPI, bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    plt.close(fig)
    print("wrote", fp.name)


if __name__ == "__main__":
    s5a = panel_5A()
    s5b = panel_5B()
    v5c = panel_5C()
    v5d = panel_5D()

    print("\n=== VERIFY (local recompute vs recorded v10 meta) ===")
    print(f"5A stars   24HPA={s5a['24HPA']}  7DPA={s5a['7DPA']}   (expect ns / ***)")
    print(f"5B.L stars 11={s5b['gdf15.L']['11HPA']} 24={s5b['gdf15.L']['24HPA']} 7={s5b['gdf15.L']['7DPA']}"
          f"   (expect ns / * / ns)")
    print(f"5B.S stars 11={s5b['gdf15.S']['11HPA']} 24={s5b['gdf15.S']['24HPA']} 7={s5b['gdf15.S']['7DPA']}"
          f"   (expect ns / ns / ns)")
    print(f"5C lm slope={v5c['slope']:+.4f} (meta {LOCK_5C['slope']:+.4f}) "
          f"int={v5c['intercept']:+.4f} (meta {LOCK_5C['intercept']:+.4f}) "
          f"R2={v5c['r2']:.4f} (meta {LOCK_5C['r2']:.4f})")
    print(f"5C spearman rho={v5c['rho']:+.4f} (meta {LOCK_5C['rho']:+.4f}) "
          f"p={v5c['p']:.5f} (meta {LOCK_5C['p']:.5f})")
    print(f"5D coupling rho={v5d['rho']:+.3f} within_MDT={v5d['within_mdt']:+.3f} "
          f"7DPT={v5d['at_7dpt']:+.3f} n={v5d['n']}   (expect +0.83 / +0.83 / +0.94)")
