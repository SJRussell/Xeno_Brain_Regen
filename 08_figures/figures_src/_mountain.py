"""Shared local fgsea-style Hallmark mountain (running-enrichment) renderer.

Fully local + provenance-tracked. Reproduces t20_round9.R's build_rank_hs and
fgsea::plotEnrichment exactly, from data we hold locally:
  - ranked list  : tables/brain/de_BRAIN_CT_vs_ND_{7DPA,24HPA}.tsv  (DESeq2 topTables)
  - symbol->entrez: tables/brain/gsea_mountain_inputs/hs_symbol_entrez.tsv  (org.Hs.eg.db, pulled once)
  - hallmark sets : tables/brain/gsea_mountain_inputs/hallmark_entrez.tsv   (msigdbr H, pulled once)
  - NES / FDR     : tables/brain/hallmark_nes_BRAIN_{24HPT,7DPT}_v10method.tsv (canonical, D1)
                    NOT mountain_fgsea_stats.tsv -- that is the direct-fgsea run, whose
                    normalization disagrees with the 3A/4A dot plots for the same result.

Ranking:  rk = sign(log2FC) * -log10(pvalue), keyed by upper(strip .L/.S of gene_id),
mapped symbol->entrez (first per symbol), dedup entrez (first in de-table order), sort desc.
Running ES = fgsea calcGseaStat: tops[i]=cumsum(|adj[S]|)/NR - (S[i]-i)/(N-m);
bottoms[i]=tops[i]-|adj[S[i]]|/NR, with adj = stat/max|stat| (gseaParam=1).

Line colour encodes direction with the deck palette: MDT-up (NES>0) coral, MDT-down blue,
so the IFN-alpha sign-flip reads as blue(24HPT) -> coral(7DPT).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import _style as S

S.apply()
ROOT = Path(__file__).resolve().parent.parent
TBL = ROOT / "tables" / "brain"
INP = TBL / "gsea_mountain_inputs"

_SYM2EZ = None
_STATS = None


def _sym2ez():
    global _SYM2EZ
    if _SYM2EZ is None:
        m = pd.read_csv(INP / "hs_symbol_entrez.tsv", sep="\t", dtype=str)
        _SYM2EZ = m.drop_duplicates("SYMBOL").set_index("SYMBOL")["ENTREZID"]  # first per symbol
    return _SYM2EZ


def build_rank(contrast):
    """Reproduce t20_round9.R build_rank_hs -> Series(index=entrez, value=stat) sorted desc."""
    de = pd.read_csv(TBL / f"de_BRAIN_{contrast}.tsv", sep="\t")
    rk = np.sign(de["log2FoldChange"]) * -np.log10(de["pvalue"])
    sym = de["gene_id"].str.replace(r"\.[LS]$", "", regex=True).str.upper()
    d = pd.DataFrame({"sym": sym.values, "rk": rk.values})
    d = d.replace([np.inf, -np.inf], np.nan).dropna(subset=["rk"])         # finite only
    d["entrez"] = d["sym"].map(_sym2ez())
    d = d.dropna(subset=["entrez"])
    d = d.drop_duplicates("entrez", keep="first")                          # first in de-table order
    s = pd.Series(d["rk"].values, index=d["entrez"].values)
    return s.sort_values(ascending=False)


def load_pathway(gs_name):
    hm = pd.read_csv(INP / "hallmark_entrez.tsv", sep="\t", dtype=str)
    return set(hm.loc[hm["gs_name"] == gs_name, "entrez"])


# contrast -> canonical v10method Hallmark table (clusterProfiler; D1, 2026-07-16)
_V10 = {"CT_vs_ND_24HPA": TBL / "hallmark_nes_BRAIN_24HPT_v10method.tsv",
        "CT_vs_ND_7DPA": TBL / "hallmark_nes_BRAIN_7DPT_v10method.tsv"}


def fgsea_stat(contrast, gs_name):
    """Canonical NES/FDR for the mountain annotation.

    The ES *curve* is engine-agnostic (clusterProfiler calls fgsea internally via
    DOSE on the same rank vector), so only the NES/FDR annotation depends on the
    normalization. Per D1 the annotation reads the v10method tables, which is what
    the dot plots (3A/4A) and Nirosha's v10 numbers use. Reading the old
    `mountain_fgsea_stats.tsv` here is what made 4B report HEME +2.39 while 4A
    reported +2.09 for the same pathway, contrast and timepoint.
    """
    global _STATS
    if contrast not in _V10:
        raise KeyError(f"no canonical v10method table mapped for contrast {contrast!r}")
    if _STATS is None:
        _STATS = {}
    if contrast not in _STATS:
        _STATS[contrast] = pd.read_csv(_V10[contrast], sep="\t")
    d = _STATS[contrast]
    r = d[d["term"] == gs_name]
    return (float(r["NES"].iloc[0]), float(r["padj"].iloc[0])) if len(r) else (np.nan, np.nan)


def running_es(rk, pathway_entrez):
    """fgsea plotEnrichment curve. rk = Series sorted desc. Returns xs, ys, es."""
    stat = rk.values.astype(float)
    N = len(stat)
    adj = stat / np.max(np.abs(stat))
    idx = {ez: i + 1 for i, ez in enumerate(rk.index)}                     # 1-indexed position
    S_pos = np.array(sorted(idx[e] for e in pathway_entrez if e in idx))
    m = len(S_pos)
    rAdj = np.abs(adj[S_pos - 1])                                          # gseaParam=1
    NR = rAdj.sum()
    rcum = np.cumsum(rAdj) / NR
    i = np.arange(1, m + 1)
    tops = rcum - (S_pos - i) / (N - m)
    bottoms = tops - rAdj / NR
    xs = np.empty(2 * m); xs[0::2] = S_pos - 1; xs[1::2] = S_pos
    ys = np.empty(2 * m); ys[0::2] = bottoms; ys[1::2] = tops
    xs = np.concatenate([[0], xs, [N + 1]])
    ys = np.concatenate([[0.0], ys, [0.0]])
    allext = np.concatenate([tops, bottoms])
    es = allext[np.argmax(np.abs(allext))]
    return xs, ys, es, S_pos, N


def render_mountain(ax, contrast, gs_name, title, *, title_size=None, tick_size=None,
                    label_size=None, annot_size=None, title_weight=None):
    """2026-09-06 (Levin round 2): explicit size/weight overrides so the Fig 4 mountains
    can match the figure_rebuild house sizes (regular-weight titles, 9-pt-class ticks).
    Defaults reproduce the previous behaviour (DECK_FS_* sizes, rcParams title weight)."""
    tick_size = S.DECK_FS_TICK if tick_size is None else tick_size
    label_size = (tick_size + 1) if label_size is None else label_size
    annot_size = S.DECK_FS_LEGEND if annot_size is None else annot_size
    rk = build_rank(contrast)
    xs, ys, es, S_pos, N = running_es(rk, load_pathway(gs_name))
    nes, padj = fgsea_stat(contrast, gs_name)
    color = S.COL_MDT if nes >= 0 else S.COL_ND                            # up in MDT coral, down blue
    ax.axhline(0, color="#BBBBBB", lw=0.6, zorder=1)
    ax.plot(xs, ys, color=color, lw=1.3, zorder=3)
    # leading-edge rug (barcode) at the bottom band
    y0, y1 = min(ys.min(), 0), max(ys.max(), 0)
    span = (y1 - y0) or 1.0
    rug_hi = y0 - 0.04 * span
    rug_lo = y0 - 0.13 * span
    ax.vlines(S_pos, rug_lo, rug_hi, color="#444444", lw=0.4, zorder=2)
    ax.set_xlim(0, N)
    ax.set_ylim(rug_lo - 0.04 * span, y1 + 0.10 * span)
    ax.set_xlabel("Rank in ordered gene list", fontsize=label_size)
    ax.set_ylabel("Running ES", fontsize=label_size)
    ax.tick_params(labelsize=tick_size)
    tkw = {} if title_weight is None else {"fontweight": title_weight}
    ax.set_title(title, fontsize=title_size or S.DECK_FS_TITLE, loc="center", **tkw)
    ax.annotate(f"NES = {nes:+.2f}\nFDR = {padj:.2g}", xy=(0.96, 0.92),
                xycoords="axes fraction", ha="right", va="top",
                fontsize=annot_size, color="#333333")
    return es, nes


def _validate():
    """Confirm the local rank reproduction matches the Nibi-exported rk exactly."""
    for ct in ("CT_vs_ND_7DPA", "CT_vs_ND_24HPA"):
        mine = build_rank(ct)
        ref = pd.read_csv(INP / f"rk_{ct}.tsv", sep="\t", dtype={"entrez": str})
        ref_s = pd.Series(ref["stat"].values, index=ref["entrez"].values)
        same_len = len(mine) == len(ref_s)
        # align on entrez, compare stat
        j = mine.reindex(ref_s.index)
        max_abs = float(np.nanmax(np.abs(j.values - ref_s.values))) if same_len else np.nan
        order_ok = list(mine.index) == list(ref_s.index)
        print(f"{ct}: local={len(mine)} ref={len(ref_s)} same_len={same_len} "
              f"order_match={order_ok} max_abs_dstat={max_abs:.3e}")


if __name__ == "__main__":
    _validate()
