"""Proportional 2-circle Venns of temporal DEGs (CT vs ND, 7 DPA vs 24 HPA),
one panel per tissue. Also writes a multi-tab xlsx of CT-only / ND-only /
shared gene partitions.

Source: pipeline_pulls/{brain,blast}_{CT,ND}_7v24_sig.tsv
  (extracted from deg/deseq2_{BRAIN,BLAST}_primary.rds on Nibi via inline
   R, already filtered to padj <= 0.05)

Outputs (in figures_current/v10_pulls/brain/):
  temporal_venn_CT_vs_ND_brain.{pdf,png}
  temporal_venn_CT_vs_ND_blastema.{pdf,png}
  temporal_venn_CT_vs_ND_both_tissues.{pdf,png}
  temporal_venn_CT_vs_ND_sets.xlsx
"""
import csv
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

ROOT = Path(
    r"C:\Users\srussell\OneDrive - Wilfrid Laurier University"
    r"\Research-related\Nirosha_Xeno_Apr2026"
)
PULLS = ROOT / "pipeline_pulls"
OUT_DIR = ROOT / "figures_current" / "v10_pulls" / "brain"

FDR = 0.05

TISSUES = [
    ("brain",    "Brain",    PULLS / "brain_CT_7v24_sig.tsv",  PULLS / "brain_ND_7v24_sig.tsv"),
    ("blastema", "Blastema", PULLS / "blast_CT_7v24_sig.tsv",  PULLS / "blast_ND_7v24_sig.tsv"),
]


def load_tsv(path):
    out = {}
    with open(path, encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            g = row["gene_id"]
            try:
                out[g] = (float(row["log2FoldChange"]), float(row["padj"]))
            except (KeyError, TypeError, ValueError):
                pass
    return out


def render_single_panel(ax, ax_txt, tissue_label, ct, nd, with_shared_names=False):
    ct_set = set(ct); nd_set = set(nd)
    inter = ct_set & nd_set
    ct_only = ct_set - nd_set
    nd_only = nd_set - ct_set
    pct_nd = 100.0 * len(inter) / max(1, len(nd_set))
    pct_ct = 100.0 * len(inter) / max(1, len(ct_set))

    v = venn2(
        subsets=(len(ct_only), len(nd_only), len(inter)),
        set_labels=("", ""),
        set_colors=("#C0392B", "#1F4E78"),
        alpha=0.55, ax=ax,
    )
    venn2_circles(
        subsets=(len(ct_only), len(nd_only), len(inter)),
        ax=ax, linewidth=1.2, linestyle="solid", color="black",
    )
    for lbl in v.subset_labels:
        if lbl is not None:
            lbl.set_fontsize(13)
            lbl.set_fontweight("bold")
    ax.set_title(
        f"{tissue_label}: temporal DEGs (7 DPA vs 24 HPA, padj < {FDR})",
        fontsize=13, pad=10,
    )

    # Right text panel
    if ax_txt is not None:
        ax_txt.axis("off")
        y = 0.94
        ax_txt.text(0.0, y, f"{tissue_label}: counts",
                    fontsize=12, fontweight="bold", color="#222")
        y -= 0.07
        rows = [
            ("CT only",            f"{len(ct_only):>5,}"),
            ("ND only",            f"{len(nd_only):>5,}"),
            ("Shared",             f"{len(inter):>5,}"),
            ("",                   ""),
            ("CT total",           f"{len(ct_set):>5,}"),
            ("ND total",           f"{len(nd_set):>5,}"),
        ]
        for label, val in rows:
            if not label:
                y -= 0.03; continue
            ax_txt.text(0.0, y, label, fontsize=11, color="#333")
            ax_txt.text(0.95, y, val, fontsize=11, color="#333",
                        ha="right", family="monospace")
            y -= 0.055
        y -= 0.04
        if len(nd_set) > 0:
            ax_txt.text(0.0, y,
                f"{pct_nd:.0f}% of ND-temporal genes\nare also CT-temporal.",
                fontsize=10, color="#444", style="italic")
            y -= 0.10
        if with_shared_names and 0 < len(inter) <= 20:
            ax_txt.text(0.0, y, "Shared genes:",
                        fontsize=11, fontweight="bold", color="#6C3483")
            y -= 0.05
            for g in sorted(inter):
                ax_txt.text(0.04, y, g, fontsize=10, color="#6C3483",
                            family="monospace")
                y -= 0.045

    return dict(
        ct_set=ct_set, nd_set=nd_set, inter=inter,
        ct_only=ct_only, nd_only=nd_only,
        pct_nd=pct_nd, pct_ct=pct_ct,
    )


def render_tissue_panel(key, label, ct_path, nd_path):
    ct = load_tsv(ct_path); nd = load_tsv(nd_path)
    print(f"\n{label}: CT n={len(ct):,}, ND n={len(nd):,}")
    fig = plt.figure(figsize=(11.5, 6.0), dpi=220)
    gs = fig.add_gridspec(1, 2, width_ratios=[2.4, 1.0])
    ax = fig.add_subplot(gs[0, 0])
    ax_txt = fig.add_subplot(gs[0, 1])
    res = render_single_panel(
        ax, ax_txt, label, ct, nd,
        with_shared_names=(len(ct & nd if False else (set(ct) & set(nd))) <= 20),
    )
    print(f"  CT only={len(res['ct_only']):,}  ND only={len(res['nd_only']):,}  shared={len(res['inter']):,}")
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"temporal_venn_CT_vs_ND_{key}.{ext}",
                    bbox_inches="tight")
    plt.close(fig)
    return ct, nd, res


def render_both_tissues(brain_pkg, blast_pkg):
    """Single 1x2 panel with brain on left, blastema on right (no side text)."""
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 6.0), dpi=220)
    for ax, (label, (ct, nd, res)) in zip(axes, [("Brain", brain_pkg), ("Blastema", blast_pkg)]):
        ct_set, nd_set, inter = res["ct_set"], res["nd_set"], res["inter"]
        ct_only, nd_only = res["ct_only"], res["nd_only"]
        v = venn2(
            subsets=(len(ct_only), len(nd_only), len(inter)),
            set_labels=("CT (regenerating)", "ND (control)"),
            set_colors=("#C0392B", "#1F4E78"),
            alpha=0.55, ax=ax,
        )
        venn2_circles(
            subsets=(len(ct_only), len(nd_only), len(inter)),
            ax=ax, linewidth=1.2, color="black",
        )
        for lbl in v.subset_labels:
            if lbl is not None:
                lbl.set_fontsize(12)
                lbl.set_fontweight("bold")
        for lbl in v.set_labels:
            if lbl is not None:
                lbl.set_fontsize(11)
        ax.set_title(
            f"{label}: temporal DEGs (7 DPA vs 24 HPA, padj < {FDR})\n"
            f"CT={len(ct_set):,}   ND={len(nd_set):,}   shared={len(inter):,} "
            f"({100*len(inter)/max(1,len(nd_set)):.0f}% of ND)",
            fontsize=11, pad=10,
        )
    plt.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT_DIR / f"temporal_venn_CT_vs_ND_both_tissues.{ext}",
                    bbox_inches="tight")
    plt.close(fig)


def build_xlsx(brain_ct, brain_nd, brain_res, blast_ct, blast_nd, blast_res, out_path):
    wb = Workbook()
    ws = wb.active
    ws.title = "README"
    for line in [
        "Temporal DEGs (CT vs ND, 7 DPA vs 24 HPA) per tissue, with set partitions for the Venns.",
        "",
        "Engine     : DESeq2 Wald, per-gene (homeolog-summed primary fit)",
        f"FDR cutoff : padj <= {FDR}",
        "Source     : deg/deseq2_{BRAIN,BLAST}_primary.rds on Nibi",
        "             contrasts CT_7DPA_vs_24HPA and ND_7DPA_vs_24HPA",
        "",
        "Brain",
        f"  CT only          : {len(brain_res['ct_only']):>5,}",
        f"  ND only          : {len(brain_res['nd_only']):>5,}",
        f"  Shared           : {len(brain_res['inter']):>5,}",
        f"  CT total         : {len(brain_res['ct_set']):>5,}",
        f"  ND total         : {len(brain_res['nd_set']):>5,}",
        "",
        "Blastema",
        f"  CT only          : {len(blast_res['ct_only']):>5,}",
        f"  ND only          : {len(blast_res['nd_only']):>5,}",
        f"  Shared           : {len(blast_res['inter']):>5,}",
        f"  CT total         : {len(blast_res['ct_set']):>5,}",
        f"  ND total         : {len(blast_res['nd_set']):>5,}",
        "",
        "Each set tab below lists gene_id, log2FoldChange and padj for both arms",
        "(value blank where the gene is not in that arm's sig list).",
    ]:
        ws.append([line])
    ws.column_dimensions["A"].width = 100
    ws["A1"].font = Font(bold=True, size=13)

    header_font = Font(bold=True, color="FFFFFF")
    center = Alignment(horizontal="center", vertical="center")

    def add_tab(name, gene_iter, color, ct, nd):
        sh = wb.create_sheet(name)
        sh.append(["gene_id", "log2FC_CT_7v24", "padj_CT_7v24",
                   "log2FC_ND_7v24", "padj_ND_7v24"])
        for c in range(1, 6):
            cell = sh.cell(row=1, column=c)
            cell.font = header_font
            cell.fill = PatternFill("solid", start_color=color)
            cell.alignment = center

        def sort_key(g):
            return min(ct.get(g, (None, 1.0))[1],
                       nd.get(g, (None, 1.0))[1])

        for g in sorted(gene_iter, key=sort_key):
            ct_lfc, ct_p = ct.get(g, (None, None))
            nd_lfc, nd_p = nd.get(g, (None, None))
            sh.append([g, ct_lfc, ct_p, nd_lfc, nd_p])
        sh.freeze_panes = "A2"
        for i, w in {1: 22, 2: 16, 3: 14, 4: 16, 5: 14}.items():
            sh.column_dimensions[get_column_letter(i)].width = w
        for r in range(2, sh.max_row + 1):
            sh.cell(row=r, column=2).number_format = "0.000"
            sh.cell(row=r, column=3).number_format = "0.00E+00"
            sh.cell(row=r, column=4).number_format = "0.000"
            sh.cell(row=r, column=5).number_format = "0.00E+00"

    # Brain tabs
    add_tab("Brain_CT_only",         brain_res["ct_only"], "C0392B", brain_ct, brain_nd)
    add_tab("Brain_ND_only",         brain_res["nd_only"], "1F4E78", brain_ct, brain_nd)
    add_tab("Brain_Shared",          brain_res["inter"],   "6C3483", brain_ct, brain_nd)
    # Blastema tabs
    add_tab("Blastema_CT_only",      blast_res["ct_only"], "C0392B", blast_ct, blast_nd)
    add_tab("Blastema_ND_only",      blast_res["nd_only"], "1F4E78", blast_ct, blast_nd)
    add_tab("Blastema_Shared",       blast_res["inter"],   "6C3483", blast_ct, blast_nd)

    wb.save(out_path)


def main():
    pkgs = {}
    for key, label, ct_p, nd_p in TISSUES:
        pkgs[key] = render_tissue_panel(key, label, ct_p, nd_p)
    render_both_tissues(pkgs["brain"], pkgs["blastema"])
    build_xlsx(
        pkgs["brain"][0], pkgs["brain"][1], pkgs["brain"][2],
        pkgs["blastema"][0], pkgs["blastema"][1], pkgs["blastema"][2],
        OUT_DIR / "temporal_venn_CT_vs_ND_sets.xlsx",
    )
    print(f"\nWrote PNG/PDF and xlsx to {OUT_DIR}")


if __name__ == "__main__":
    main()
