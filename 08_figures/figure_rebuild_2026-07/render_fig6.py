"""Fig 6 local panels - brain suppresses synaptic-transmission programs by 7 DPT.

  6A  synaptic GO:BP NES bars at 7 DPT   (LOCKED to 10 recorded meta terms)
  6B  synaptic + structural gene trajectories, 3x3 (voom E; trt:time beta LOCKED to meta)
  6C  Liao zbtb18+ excitatory-neuron UCell trajectory (per_celltype_scores.tsv)
  6D  synaptic GO:BP 24 HPT vs 7 DPT side-by-side (full GO:BP GSEA, genome-wide clusterProfiler BH)

House style (figure_rebuild _style): Arial, FS_TITLE=10 centered, MDT/ND palette, HPT/DPT.
"CT" in the data displays as "MDT".

Inputs (all local):
  6A  tables/_meta/brain_synaptic_go_7dpt_deseq2_meta.yaml  (locked; from clusterProfiler GO:BP)
  6B  tables/brain/voom_E_BRAIN_primary.tsv (9 genes) + beta/p locked to
      tables/_meta/brain_synaptic_gene_trajectories_v10_meta.yaml
  6C  tables/comparators/liao/per_celltype_scores.tsv ["Excitatory neuron_zbtb18 high"]
      + per_celltype_interaction.tsv (est_int/p_int)
  6D  tables/brain/go_bp_rds/goBP_FULL_BRAIN_CT_vs_ND_{24HPA,7DPA}.tsv
      (re-run of the canonical t06 gseGO with pvalueCutoff=1 so all terms carry NES+padj;
       Nibi pull 2026-07-11)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
import _style as S
from _traj import draw_trajectory, treatment_legend

BASE = Path(__file__).resolve().parent.parent
TBL  = BASE / "tables"
META = pd.read_csv(TBL / "sample_metadata.tsv", sep="\t")
OUT  = Path(__file__).resolve().parent / "panels_png"; OUT.mkdir(exist_ok=True)
FPNG = BASE / "figures_png"; FPNG.mkdir(exist_ok=True)


def _save(fig, name, figsize):
    for d in (OUT, FPNG):
        fig.savefig(d / f"{name}.png", dpi=S.DEFAULT_DPI, bbox_inches=Bbox.from_bounds(0, 0, *figsize))
    plt.close(fig)
    print("wrote", name)


# ---------------------------------------------------------------- 6A synaptic GO bars
# LOCKED to brain_synaptic_go_7dpt_deseq2_meta.yaml (clusterProfiler GSEA GO:BP, DESeq2 ranks).
GO_7DPT = [  # (label, NES, q-label)  most-significant first
    ("synaptic signaling", -1.91, "1.6e-05"),
    ("axon development", -1.88, "2.5e-05"),
    ("chemical synaptic transmission", -1.92, "4.2e-05"),
    ("anterograde trans-synaptic signaling", -1.92, "4.2e-05"),
    ("trans-synaptic signaling", -1.91, "7.8e-05"),
    ("modulation of chemical synaptic transmission", -2.10, "0.00015"),
    ("regulation of trans-synaptic signaling", -2.10, "0.00015"),
    ("regulation of synapse structure or activity", -1.94, "0.0032"),
    ("synapse organization", -1.83, "0.0053"),
    ("axon guidance", -1.70, "0.0075"),
]


def panel_6A(figsize=(5.2, 3.3), axleft=0.505):
    labels = [t[0] for t in GO_7DPT]
    nes = np.array([t[1] for t in GO_7DPT])
    qlab = [t[2] for t in GO_7DPT]
    y = np.arange(len(labels))[::-1]                       # first item at top
    norm = Normalize(vmin=-2.15, vmax=-1.65)
    cmap = plt.cm.Blues_r
    colors = cmap(norm(nes))

    FS = figsize
    fig = plt.figure(figsize=FS)
    ax = fig.add_axes([axleft, 0.14, 0.845 - axleft, 0.78])
    ax.barh(y, nes, color=colors, edgecolor="#2C4A6E", linewidth=0.4, height=0.72)
    for yi, q in zip(y, qlab):
        ax.text(0.06, yi, f"q = {q}", va="center", ha="left", fontsize=S.FS_ANNOT - 1,
                color="#333333")
    ax.axvline(0, color="grey", lw=0.6)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=S.FS_TICK - 1.5)
    ax.set_xlim(-2.6, 1.1)
    ax.set_xlabel("NES (MDT vs ND, 7 DPT)", labelpad=2)
    # title as a figure-level text at a fixed height so it lines up with 6B's title
    fig.suptitle("Brain synaptic GO BP NES at 7 DPT", fontsize=S.FS_TITLE, color=S.COL_TEXT,
                 x=(axleft + 0.845) / 2, y=0.955)
    ax.grid(True, axis="x", alpha=0.18)
    cax = fig.add_axes([0.885, 0.30, 0.022, 0.44])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.set_label("NES", fontsize=S.FS_LEGEND); cb.set_ticks([-1.8, -1.9, -2.0])
    cb.ax.tick_params(labelsize=S.FS_TICK - 1)
    _save(fig, "fig6A_synaptic_go_7dpt", FS)


# ---------------------------------------------------------------- 6B gene trajectories
# 9 genes; trt:time beta/p LOCKED to brain_synaptic_gene_trajectories_v10_meta.yaml.
GENES_6B = [  # (gene, beta, p)  row-major 3x3 as in the figure
    ("slc17a7.L", -0.4219, 0.001153), ("gria1.L", -0.2163, 0.01418), ("syn2.S", -0.3499, 0.00236),
    ("grin2b.L", -0.117, 0.2164),     ("syp.S", -0.2756, 0.004496),  ("synm.L", 0.1008, 0.6227),
    ("synm.S", 3.165, 0.01585),       ("map2.L", -0.05555, 0.4345),  ("mapt.L", -0.02507, 0.8979),
]


def panel_6B(figsize=(5.84, 5.0)):
    # 2026-09-06: taller (4.14 -> 5.0 in), wider gaps, figure-level legend, so the 3x3 grid
    # stays legible at the review-round font sizes (6C moves down on the slide).
    vdf = pd.read_csv(TBL / "brain/voom_E_BRAIN_primary.tsv", sep="\t", index_col=0)
    msub = META[META["tissue"] == "BRAIN"][["sample", "treatment", "time"]]
    times = ["24HPA", "7DPA"]
    fig = plt.figure(figsize=figsize)
    lm, bm, gx, gy = 0.115, 0.125, 0.075, 0.125          # margins / gaps (bm clears the bottom legend)
    rm, tm = 0.02, 0.135                                 # right / top margins (tm clears suptitle + legend)
    pw = (1 - lm - rm - 2 * gx) / 3
    ph = (1 - bm - tm - 2 * gy) / 3
    h = None
    for idx, (gene, beta, p) in enumerate(GENES_6B):
        r, c = divmod(idx, 3)
        ax = fig.add_axes([lm + c * (pw + gx), bm + (2 - r) * (ph + gy), pw, ph])
        expr = vdf.loc[gene].to_dict()
        h = draw_trajectory(ax, expr, msub, times, {}, ylabel_text=None)
        ax.set_xlabel("")                                # single shared labels added below
        ax.set_title(gene, fontsize=S.FS_AXIS - 2, fontstyle="italic", color="#333333", pad=2)
        ax.tick_params(labelsize=S.FS_TICK - 2)
        ax.text(0.03, 0.04, f"β = {beta:+.2f}, p = {p:.2g}", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=S.FS_ANNOT - 2.5, color="#555555")
    # 2026-09-06 (round 2, Stewart: "6B title slightly clipped"): the title was centred at
    # x=0.33 to clear the top-right legend and ran off the left edge of the PNG. Title now
    # centred; the treatment legend sits bottom-right beside the shared "Timepoint" label.
    fig.legend([h["ND"], h["MDT"]], ["ND", "MDT"], loc="lower right",
               bbox_to_anchor=(0.985, 0.0), frameon=False, ncol=2, fontsize=S.FS_LEGEND,
               handlelength=1.4, columnspacing=0.9, borderaxespad=0.2)
    cx, cy = (lm + (1 - rm)) / 2, (bm + (1 - tm)) / 2    # centre of the 3x3 block
    fig.text(cx, 0.02, "Timepoint", ha="center", va="bottom", fontsize=S.FS_AXIS)
    fig.text(0.022, cy, "Normalized expression", rotation=90, ha="left", va="center",
             fontsize=S.FS_AXIS)
    fig.suptitle("Brain synaptic and structural gene trajectories", fontsize=S.FS_TITLE,
                 color=S.COL_TEXT, x=0.5, y=0.965)
    _save(fig, "fig6B_synaptic_gene_trajectories", figsize)


# ---------------------------------------------------------------- 6C Liao UCell
def panel_6C(figsize=(5.84, 3.61)):
    sc = pd.read_csv(TBL / "comparators/liao/per_celltype_scores.tsv", sep="\t")
    col = "Excitatory neuron_zbtb18 high"
    expr = dict(zip(sc["sample"], sc[col]))
    inter = pd.read_csv(TBL / "comparators/liao/per_celltype_interaction.tsv", sep="\t")
    row = inter[inter["celltype"] == col].iloc[0]
    beta, p, fdr = float(row["est_int"]), float(row["p_int"]), float(row["fdr_int"])
    n_tested = int(inter["p_int"].notna().sum())
    msub = META[META["tissue"] == "BRAIN"][["sample", "treatment", "time"]]

    FS = figsize
    fig = plt.figure(figsize=FS)
    ax = fig.add_axes([0.165, 0.135, 0.80, 0.78])
    h = draw_trajectory(ax, expr, msub, ["24HPA", "7DPA"], {}, ylabel_text="UCell score")
    ax.set_title("Liao zbtb18+ excitatory-neuron signature in brain samples",
                 fontsize=S.FS_TITLE - 1.5, color=S.COL_TEXT, pad=6)
    # correct p + BH FDR (across the cell-type signatures tested) as a bottom-left inset
    ax.text(0.03, 0.03, f"β = {beta:+.3f}\np = {p:.3f}\nFDR = {fdr:.3f}", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=S.FS_ANNOT, color="#333333", linespacing=1.35)
    treatment_legend(ax, h, loc="upper right")
    _save(fig, "fig6C_liao_zbtb18_ucell", FS)
    return beta, p, fdr, n_tested


# ---------------------------------------------------------------- 6D synaptic GO 24v7
SYN_KEY = r"synap|neurotransmitter|nerve impulse"   # synapse-specific (matches the panel terms)
# The manuscript defines gene-set significance at FDR < 0.05; 6D used to render at
# 0.10, i.e. the panel was more permissive than the paper's own stated threshold.
# Verified 2026-07-16 from go_bp_rds/: of 35 synaptic-filtered terms, 7 DPT has 17 at
# padj<0.10 and 15 at padj<0.05, all 15 with negative NES; 24 HPT has 0 at BOTH
# thresholds. The late-only suppression conclusion is unchanged, so this only
# tightens the claim. (Stewart's older "45 of 98" used a within-synaptic-set BH that
# was rejected -- it must not reach the manuscript.)
SYN_FDR = 0.05


def panel_6D(figsize=(7.1, 5.4)):
    D = TBL / "brain/go_bp_rds"
    d24 = pd.read_csv(D / "goBP_FULL_BRAIN_CT_vs_ND_24HPA.tsv", sep="\t")
    d7 = pd.read_csv(D / "goBP_FULL_BRAIN_CT_vs_ND_7DPA.tsv", sep="\t")
    for d in (d24, d7):
        d["dl"] = d["Description"].str.lower()
    ids = (set(d24[d24["dl"].str.contains(SYN_KEY, regex=True)]["Description"]) &
           set(d7[d7["dl"].str.contains(SYN_KEY, regex=True)]["Description"]))

    # Significance = clusterProfiler genome-wide BH p.adjust (same run as 6A). This gives
    # 0 synaptic terms sig at 24 HPT (the "late-only" story) and is defensible in review.
    def prep(d):
        return d[d["Description"].isin(ids)].set_index("Description")
    s24, s7 = prep(d24), prep(d7)
    n_sig7 = int((s7["p.adjust"] < SYN_FDR).sum()); n_sig24 = int((s24["p.adjust"] < SYN_FDR).sum())
    m = pd.DataFrame({"n24": s24["NES"], "p24": s24["p.adjust"],
                      "n7": s7["NES"], "p7": s7["p.adjust"]})
    # honest framing (Stewart): show only terms significant at 7 DPT; drop the always-n.s. ones.
    top = m[m["p7"] < SYN_FDR].sort_values("p7").iloc[::-1]   # most-significant at top
    assert n_sig24 == 0, f"expected 0 sig synaptic terms at 24 HPT, got {n_sig24}"
    y = np.arange(len(top))

    def barcol(nes):
        return S.COL_MDT if nes > 0 else S.COL_ND         # MDT>ND red, MDT<ND blue

    FS = figsize
    fig = plt.figure(figsize=FS)
    axL = fig.add_axes([0.375, 0.14, 0.29, 0.74])
    axR = fig.add_axes([0.695, 0.14, 0.29, 0.74])
    for ax, ncol, pcol, ttl in ((axL, "n24", "p24", "24 HPT"), (axR, "n7", "p7", "7 DPT")):
        for yi, (_, r) in zip(y, top.iterrows()):
            sig = r[pcol] < SYN_FDR
            ax.barh(yi, r[ncol], color=barcol(r[ncol]), edgecolor="#333333", linewidth=0.3,
                    height=0.74, alpha=1.0 if sig else 0.30)
        ax.axvline(0, color="#555555", lw=0.6)
        ax.set_xlim(-2.5, 2.5); ax.set_ylim(-0.6, len(top) - 0.4)
        ax.set_title(ttl, fontsize=S.FS_AXIS, color="#333333", pad=3)
        ax.set_xlabel("NES (MDT vs ND)", labelpad=2)
        ax.grid(True, axis="x", alpha=0.15)
    labels = [d[:1].upper() + d[1:] for d in top.index]
    axL.set_yticks(y); axL.set_yticklabels(labels, fontsize=S.FS_TICK - 2)
    axR.set_yticks(y); axR.set_yticklabels([])
    fig.text(0.5, 0.965, "Synaptic GO BP: 24 HPT vs 7 DPT (brain, MDT vs ND)",
             ha="center", fontsize=S.FS_TITLE, color=S.COL_TEXT)
    fig.text(0.5, 0.935, f"{n_sig7} synaptic GO:BP terms significant at 7 DPT (FDR < {SYN_FDR:g}); "
             "all n.s. at 24 HPT.", ha="center", fontsize=S.FS_ANNOT - 1, color="#666666")
    handles = [Line2D([0], [0], marker="s", ls="none", mfc=S.COL_ND, mec="none", ms=8, label="MDT < ND"),
               Line2D([0], [0], marker="s", ls="none", mfc=S.COL_MDT, mec="none", ms=8, label="MDT > ND")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=S.FS_LEGEND, bbox_to_anchor=(0.66, 0.005))
    _save(fig, "fig6D_synaptic_go_24v7", FS)
    return len(ids), n_sig7, n_sig24


# ------------------------------------------------- new 6A: synaptic + axon/dendrite, 24v7
# 2026-09-01 (Stewart): 6D promoted to panel A with the BROADER keyword filter that the old
# locked 6A meta used (adds axon development, axon guidance, dendrite morphogenesis -> 18
# terms sig at 7 DPT, still 0 at 24 HPT). Old panel_6A / panel_6D kept for provenance.
BROAD_KEY = r"synap|presynaptic|postsynaptic|neurotransmitter|axon(?!o)|dendrit|nerve impulse"


def panel_6A_expanded(figsize=(6.4, 7.9)):
    import textwrap
    D = TBL / "brain/go_bp_rds"
    d24 = pd.read_csv(D / "goBP_FULL_BRAIN_CT_vs_ND_24HPA.tsv", sep="\t")
    d7 = pd.read_csv(D / "goBP_FULL_BRAIN_CT_vs_ND_7DPA.tsv", sep="\t")
    for d in (d24, d7):
        d["dl"] = d["Description"].str.lower()
    ids = (set(d24[d24["dl"].str.contains(BROAD_KEY, regex=True)]["Description"]) &
           set(d7[d7["dl"].str.contains(BROAD_KEY, regex=True)]["Description"]))

    def prep(d):
        return d[d["Description"].isin(ids)].set_index("Description")
    s24, s7 = prep(d24), prep(d7)
    n_sig7 = int((s7["p.adjust"] < SYN_FDR).sum()); n_sig24 = int((s24["p.adjust"] < SYN_FDR).sum())
    m = pd.DataFrame({"n24": s24["NES"], "p24": s24["p.adjust"],
                      "n7": s7["NES"], "p7": s7["p.adjust"]})
    top = m[m["p7"] < SYN_FDR].sort_values("p7").iloc[::-1]   # most-significant at top
    assert n_sig24 == 0, f"expected 0 sig terms at 24 HPT, got {n_sig24}"
    y = np.arange(len(top))

    def barcol(nes):
        return S.COL_MDT if nes > 0 else S.COL_ND

    FS = figsize
    fig = plt.figure(figsize=FS)
    # 2026-09-06 (round 2, Stewart): subtitle dropped (house rule, no subtitles); the term count
    # lives in the legend. The bars take the freed height.
    axL = fig.add_axes([0.40, 0.125, 0.27, 0.78])
    axR = fig.add_axes([0.70, 0.125, 0.27, 0.78])
    for ax, ncol, pcol, ttl in ((axL, "n24", "p24", "24 HPT"), (axR, "n7", "p7", "7 DPT")):
        for yi, (_, r) in zip(y, top.iterrows()):
            sig = r[pcol] < SYN_FDR
            ax.barh(yi, r[ncol], color=barcol(r[ncol]), edgecolor="#333333", linewidth=0.3,
                    height=0.74, alpha=1.0 if sig else 0.30)
        ax.axvline(0, color="#555555", lw=0.6)
        ax.set_xlim(-2.5, 2.5); ax.set_ylim(-0.6, len(top) - 0.4)
        ax.set_title(ttl, fontsize=S.FS_AXIS, color="#333333", pad=3)
        ax.set_xlabel("NES (MDT vs ND)", labelpad=2)
        ax.grid(True, axis="x", alpha=0.15)
    labels = [textwrap.fill(d[:1].upper() + d[1:], 26) for d in top.index]   # 2 lines max
    axL.set_yticks(y); axL.set_yticklabels(labels, fontsize=S.FS_TICK - 1.5, linespacing=1.0)
    axR.set_yticks(y); axR.set_yticklabels([])
    fig.text(0.55, 0.965, "Synaptic and axon GO BP: 24 HPT vs 7 DPT (brain)",
             ha="center", fontsize=S.FS_TITLE, color=S.COL_TEXT)
    handles = [Line2D([0], [0], marker="s", ls="none", mfc=S.COL_ND, mec="none", ms=8, label="MDT < ND"),
               Line2D([0], [0], marker="s", ls="none", mfc=S.COL_MDT, mec="none", ms=8, label="MDT > ND")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=S.FS_LEGEND, bbox_to_anchor=(0.685, 0.004))
    _save(fig, "fig6A_synaptic_axon_24v7", FS)
    return len(ids), n_sig7, n_sig24, list(top.index)


if __name__ == "__main__":
    # current Fig 6 (2026-09-01 layout): A = expanded synaptic+axon 24v7, B = trajectories,
    # C = Liao UCell. Old panel_6A / panel_6D kept for provenance, not rendered by default.
    n_ids, sig7, sig24, terms = panel_6A_expanded()
    panel_6B()
    b, p, fdr, n_tested = panel_6C()
    print("\n=== VERIFY ===")
    print(f"6C Liao zbtb18: beta={b:+.4f} p={p:.4f} FDR={fdr:.4f} "
          f"(BH across {n_tested} tested cell-type signatures)")
    print(f"6A synaptic+axon set n={n_ids}: showing {sig7} terms sig@7DPT (padj<0.05); sig@24HPT={sig24}")
