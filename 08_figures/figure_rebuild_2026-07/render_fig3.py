"""Fig 3 local panels (early brain, 24 HPT):
  3A Hallmark NES dotplot (MDT vs ND)      local: hallmark_BRAIN_CTvND_24HPT.tsv
  3E Reactome ECM/collagen over-rep barh   local: reactome_overrep_BRAIN_24HPT.tsv
  3F Reactome ECM/collagen cnetplot (MDT)  local: reactome_cnet_edges_BRAIN_24HPT.tsv
3B (ecm trajectories), 3C (mTORC1 LE heatmap), 3D (neurovascular, drop pericyte)
are VOOM -> reused deck placeholders until the Nibi pull.
"""
from pathlib import Path
import pandas as pd, numpy as np
import matplotlib.pyplot as plt, matplotlib as mpl, networkx as nx
import matplotlib.patheffects as pe
import _style as S
from _panels import nes_dotplot, reactome_barh

BASE = Path(__file__).resolve().parent.parent
# Repointed 2026-07-17 from manuscript_lock_2026-06-08 (absent from this machine) to the
# restored canonical tables/. Same data, different filenames: the lock used
# `*_BRAIN_24HPT.tsv`, the restore uses `*_CTvsND_24HPT.tsv`.
TBL  = BASE / "tables" / "brain"
OUT  = Path(__file__).resolve().parent / "panels_png"; OUT.mkdir(exist_ok=True)

# 6 of the 14 curated ECM/collagen pathways do NOT pass adj P < 0.05 and must be
# drawn as non-significant in BOTH 3E and 3F. Verified against
# reactome_overrep_CTvsND_24HPT.tsv (direction=="ANY") 2026-07-17:
#   collagen chain trimerization 0.0600 | collagen degradation 0.1069
#   assembly of collagen fibrils 0.1236 | collagen biosynth/modifying 0.1442
#   collagen formation 0.1936          | degradation of ECM 0.2006
# The 8 that do pass are q 0.0145-0.0417 on Counts of only 2-5 genes.
NS_Q = 0.05

MUST = [
    "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION", "REACTOME_COLLAGEN_FORMATION",
    "REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES",
    "REACTOME_ASSEMBLY_OF_COLLAGEN_FIBRILS_AND_OTHER_MULTIMERIC_STRUCTURES",
    "REACTOME_COLLAGEN_CHAIN_TRIMERIZATION", "REACTOME_COLLAGEN_DEGRADATION",
    "REACTOME_CROSSLINKING_OF_COLLAGEN_FIBRILS",
    "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX",
    "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
    "REACTOME_FIBRONECTIN_MATRIX_FORMATION",
    "REACTOME_NON_INTEGRIN_MEMBRANE_ECM_INTERACTIONS", "REACTOME_LAMININ_INTERACTIONS",
    "REACTOME_PLATELET_ADHESION_TO_EXPOSED_COLLAGEN", "REACTOME_ANCHORING_FIBRIL_FORMATION",
]

# 3A Hallmark dotplot.
# Reads the v10method table (clusterProfiler), canonical per D1 2026-07-16. Previously
# read the fgsea variant `hallmark_BRAIN_CTvND_24HPT.tsv` out of manuscript_lock, which
# is a THIRD engine relative to the 4A dot plot and the mountains -- the drift PANEL_INDEX
# exists to prevent. PANEL_INDEX's 3A row already claimed v10method; this makes that true.
hm = pd.read_csv(TBL / "hallmark_nes_BRAIN_24HPT_v10method.tsv", sep="\t")
hm_sig = hm[hm["padj"] < 0.05].copy()
if len(hm_sig) < 5:
    hm_sig = hm.sort_values("padj").head(15)
nes_dotplot(hm_sig, "Brain Hallmark GSEA, 24 HPT",
            OUT / "fig3A_hallmark_24hpt.png", top_n=20, label_size=S.FS_TICK - 2.2)

# 3E Reactome ECM/collagen barh, non-significant pathways drawn hollow/grey.
rx = pd.read_csv(TBL / "reactome_overrep_CTvsND_24HPT.tsv", sep="\t")
reactome_barh(rx, "ECM/collagen Reactome pathways",
              OUT / "fig3E_reactome_overrep.png", MUST, ns_q=NS_Q)

# 3F Reactome cnetplot (MDT relabel) — WIDE short strip. Layout + draw code live in
# fig3f_layout.py so the interactive dragger (drag_fig3f.py) and this static build
# share one pipeline. If panels_png/fig3F_positions.json exists (Stewart's hand
# arrangement) it wins; otherwise the computed wide layout fills the row.
from fig3f_layout import render_panel as _render_3f
_render_3f(out_png=OUT / "fig3F_reactome_cnetplot.png")
plt.close("all")
print("wrote fig3F_reactome_cnetplot.png")
