"""Shared layout + draw code for Figure 3F (ECM/collagen Reactome cnetplot).

Imported by BOTH the static renderer (render_fig3.py) and the interactive
dragger (drag_fig3f.py) so the on-screen preview and the final publication PNG
are drawn by identical code.

Wide-strip default layout: the two collagen hubs (COL4A1/COL4A2) are separated
horizontally and the far endpoints (COL4A1 left, GP1BB right) are pinned, so the
14 pathway nodes string out across the full width instead of stacking in a fan.
That lets the panel be short and wide.

Positions can be overridden by a saved hand-arranged layout: if
panels_png/fig3F_positions.json exists it wins over the computed spring layout,
so Stewart's drag edits are reproducible on every re-render.
"""
from pathlib import Path
import json
import pandas as pd, numpy as np
import matplotlib as mpl, matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import networkx as nx
import _style as S

BASE = Path(__file__).resolve().parent.parent
# Repointed 2026-07-17 from manuscript_lock_2026-06-08 (absent from this machine) to the
# restored canonical tables/ (same data, `_CTvsND_` rather than `_BRAIN_` filenames).
TBL  = BASE / "tables" / "brain"
EDGES = TBL / "reactome_cnet_edges_CTvsND_24HPT.tsv"
OVERREP = TBL / "reactome_overrep_CTvsND_24HPT.tsv"
POS_JSON = Path(__file__).resolve().parent / "panels_png" / "fig3F_positions.json"

# Pinned nodes: the two collagen hubs held apart (shared pathways spread between
# them) plus GP1BB as a moderate right anchor. GP1BB's only edge is to "Platelet
# adhesion", so pinning GP1BB keeps that leaf (and its label) on-frame instead of
# being flung to the edge.
# Hub separation is the real lever on horizontal fill: the collagen fan hangs off
# COL4A1/COL4A2, so holding them 1.2 apart (was 0.6) spreads the dense left cluster
# across the width instead of scaling the whole blob. GP1BB anchors the thin
# platelet tail at the right.
PIN = {
    "COL4A1": (-1.20,  0.00),
    "COL4A2": ( 0.00,  0.00),
    "GP1BB":  ( 1.20,  0.00),
}
# soft starting spots for the free genes (not fixed), to steer the tail rightward
SEED = {
    "PECAM1": ( 0.25,  0.35),
    "VWF":    ( 0.60, -0.05),
}
# a couple of pathway names are too long even truncated; alias them shorter so the
# rightmost label doesn't collide with the colorbar
ALIAS = {
    "REACTOME_PLATELET_ADHESION_TO_EXPOSED_COLLAGEN": "Platelet adhesion",
}

DEFAULT_FIGSIZE = (7.4, 3.3)   # wide short strip (~20% narrower than the 9.2 draft);
                               # height 2.95 -> 3.3 on 2026-09-06 for label room at the
                               # larger review-round font sizes
# Node spread inside FRAME. Normalized coords are +-0.5 before scaling, so
# X_SPREAD=2.30 -> x in +-1.15 (fills the +-1.20 frame, leaving a little label runway)
# and Y_SQUASH=0.95 -> y in +-0.475 (was 1.15 -> +-0.575). Stewart 2026-07-17:
# "more spread out on the x-axis, ideally shorter on the y axis".
X_SPREAD = 2.05
Y_SQUASH = 0.95
AXES_RECT = [0.035, 0.06, 0.83, 0.84]
# fixed drawing frame shared by the static render and the dragger, so a dragged
# arrangement maps 1:1 to the saved PNG. The wide_layout normalizer keeps every
# node inside x in [-1,1], y in ~[-0.6,0.6]; this frame pads that.
FRAME_X = (-1.20, 1.20)
# Y frame tightened 0.80 -> 0.65: nodes span +-0.475 (Y_SQUASH), so the old frame left
# ~40% of the height as dead band. Keep some slack for label runway under ensure_inside_axes.
FRAME_Y = (-0.65, 0.65)


def short_label(pathway):
    if pathway in ALIAS:
        return ALIAS[pathway]
    s = pathway.replace("REACTOME_", "").replace("_", " ").capitalize()
    return s if len(s) <= 28 else s[:27].rstrip() + "..."


def path_significance(q=0.05):
    """short_label -> True if that curated Reactome pathway passes adj P < q.

    The cnet shows all 14 curated pathways, 6 of which are not significant (q 0.060
    to 0.201). Without this the network implies all 14 were identified as enriched.
    """
    o = pd.read_csv(OVERREP, sep="\t")
    o = o[o["direction"] == "ANY"].drop_duplicates("ID")
    return {short_label(r["ID"]): bool(r["p.adjust"] < q) for _, r in o.iterrows()}


def build_graph():
    """Return (G, df, gene_fc, path_deg, gene_nodes, path_nodes, path_sig)."""
    df = pd.read_csv(EDGES, sep="\t")
    df["path_label"] = df["pathway"].apply(short_label)
    G = nx.Graph()
    path_nodes = df["path_label"].unique().tolist()
    gene_nodes = df["gene"].unique().tolist()
    for p in path_nodes:
        G.add_node(p, kind="pathway")
    for g in gene_nodes:
        G.add_node(g, kind="gene")
    for _, r in df.iterrows():
        G.add_edge(r["path_label"], r["gene"])
    path_deg = {p: G.degree(p) for p in path_nodes}
    gene_fc = df.drop_duplicates("gene").set_index("gene")["log2FC"].to_dict()
    sig = path_significance()
    path_sig = {p: sig.get(p, True) for p in path_nodes}
    return G, df, gene_fc, path_deg, gene_nodes, path_nodes, path_sig


def wide_layout(G, seed=7):
    """Spring layout with the two hubs pinned apart, so the graph fills a wide short
    box. Returns {node: (x, y)} normalized to ~[-1,1] x."""
    init = {**SEED, **PIN}
    init = {k: np.array(v, float) for k, v in init.items() if k in G}
    fixed = [n for n in PIN if n in G]
    pos = nx.spring_layout(G, pos=init, fixed=fixed, k=0.62,
                           iterations=600, seed=seed)
    pos = {k: np.asarray(v, float) for k, v in pos.items()}
    # normalize to a stable frame: span x by midpoint (not mean) so a right-hand
    # outlier can't push everything left, and gently compress vertically
    xs = np.array([p[0] for p in pos.values()])
    ys = np.array([p[1] for p in pos.values()])
    cx = (xs.max() + xs.min()) / 2
    cy = (ys.max() + ys.min()) / 2
    sx = (xs.max() - xs.min()) or 1.0
    sy = (ys.max() - ys.min()) or 1.0
    for k in pos:
        pos[k] = np.array([(pos[k][0] - cx) / sx * X_SPREAD,
                           (pos[k][1] - cy) / sy * Y_SQUASH])
    return pos


def load_positions(G):
    """Saved hand layout if present and complete, else the computed wide layout."""
    if POS_JSON.exists():
        raw = json.loads(POS_JSON.read_text())
        if set(raw) >= set(G.nodes):
            return {n: np.asarray(raw[n], float) for n in G.nodes}
        print("  (saved positions incomplete -> recomputing wide layout)")
    return wide_layout(G)


def save_positions(pos):
    POS_JSON.parent.mkdir(parents=True, exist_ok=True)
    POS_JSON.write_text(json.dumps(
        {k: [float(v[0]), float(v[1])] for k, v in pos.items()}, indent=2))


def _psize(deg):
    return 120 + 110 * deg


def _gsize(fc):
    return 170 + 150 * abs(fc)


def draw_network(ax, G, pos, gene_fc, path_deg, gene_nodes, path_nodes,
                 norm, cmap, label_pathways=True, adjust=True, path_sig=None):
    """Draw edges + nodes + labels onto ax. Returns dict of artists for reuse.

    `path_sig`: label -> bool. Non-significant curated pathways (adj P >= 0.05) are
    drawn hollow with a dashed grey rim so they read as context, not as findings.
    """
    ax.clear()
    ax.axis("off")
    if path_sig is None:
        path_sig = {p: True for p in path_nodes}
    edge_art = nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#BBBBBB",
                                      width=0.7, alpha=0.7)
    # Significant = tan + rim (unchanged from the original look). Non-significant =
    # same tan fill, rim simply dropped. Subtle on purpose (Stewart 2026-07-17: the
    # earlier dashed-hollow rings were too loud).
    path_art = nx.draw_networkx_nodes(
        G, pos, nodelist=path_nodes, ax=ax, node_color="#E8D9C0",
        edgecolors=["#9A8C6A" if path_sig[p] else "none" for p in path_nodes],
        linewidths=[0.8 if path_sig[p] else 0.0 for p in path_nodes],
        node_size=[_psize(path_deg[p]) for p in path_nodes])
    gene_art = nx.draw_networkx_nodes(
        G, pos, nodelist=gene_nodes, ax=ax,
        node_color=[cmap(norm(gene_fc[g])) for g in gene_nodes],
        edgecolors="black", linewidths=0.9,
        node_size=[_gsize(gene_fc[g]) for g in gene_nodes])
    gtxt = nx.draw_networkx_labels(
        G, pos, labels={g: g for g in gene_nodes}, ax=ax,
        font_size=S.FS_TICK - 2, font_weight="bold")
    for t in gtxt.values():
        t.set_path_effects([pe.withStroke(linewidth=2.2, foreground="white")])

    ptexts = []
    if label_pathways:
        xs = [pos[n][0] for n in G.nodes]
        ys = [pos[n][1] for n in G.nodes]
        ptexts = [ax.text(pos[p][0], pos[p][1], p, fontsize=S.FS_LEGEND - 1.6,
                          color="#333333" if path_sig[p] else "#8E8E8E",
                          ha="center", va="center") for p in path_nodes]
        for t in ptexts:
            t.set_path_effects([pe.withStroke(linewidth=1.8, foreground="white")])
        if adjust:
            try:
                from adjustText import adjust_text
                # ensure_inside_axes: without it the fixed Bbox CLIPS left-anchored
                # labels ("Anchoring fibril formation" sheared off the left edge once
                # X_SPREAD grew). Same lesson as the Fig 2 volcanoes.
                # gene labels are passed as `objects` so pathway labels also avoid THEM
                # (at the 2026-09 font sizes they collided with COL4A1 / VWF / GP1BB)
                adjust_text(ptexts, x=xs, y=ys, ax=ax, objects=list(gtxt.values()),
                            expand=(1.6, 2.1), force_text=(0.8, 1.2),
                            force_static=(0.9, 1.3), force_explode=(0.3, 0.5),
                            ensure_inside_axes=True, time_lim=30,
                            arrowprops=dict(arrowstyle="-", color="#999999", lw=0.4))
            except Exception as e:
                print("  (adjustText:", e, ")")

    # fixed frame (shared with the dragger) so drags map 1:1 to the saved PNG
    ax.set_xlim(*FRAME_X)
    ax.set_ylim(*FRAME_Y)
    return dict(edges=edge_art, path=path_art, gene=gene_art,
                gene_txt=gtxt, path_txt=ptexts)


def make_norm_cmap(gene_fc):
    cmap = mpl.colormaps[S.CMAP_DIVERGING]
    vmax = max(abs(v) for v in gene_fc.values())
    norm = mpl.colors.Normalize(vmin=-vmax, vmax=vmax)
    return norm, cmap


def render_panel(pos=None, figsize=DEFAULT_FIGSIZE, out_png=None, adjust=True):
    """Full publication panel: title + colorbar + MDT/ND cues. Saves PNG."""
    S.apply()
    G, df, gene_fc, path_deg, gene_nodes, path_nodes, path_sig = build_graph()
    if pos is None:
        pos = load_positions(G)
    norm, cmap = make_norm_cmap(gene_fc)

    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(AXES_RECT)
    draw_network(ax, G, pos, gene_fc, path_deg, gene_nodes, path_nodes,
                 norm, cmap, adjust=adjust, path_sig=path_sig)
    ax.set_title("ECM/collagen Reactome network, 24 HPT",
                 fontsize=S.FS_TITLE, color=S.COL_TEXT, pad=6)

    cax = fig.add_axes([0.895, 0.22, 0.017, 0.50])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("log2FC (MDT vs ND)", fontsize=S.FS_LEGEND)
    cb.ax.tick_params(labelsize=S.FS_LEGEND)
    cb.ax.text(0.5, 1.06, "up in MDT", transform=cb.ax.transAxes, ha="center",
               va="bottom", fontsize=S.FS_LEGEND, color=S.COL_MDT)
    cb.ax.text(0.5, -0.06, "up in ND", transform=cb.ax.transAxes, ha="center",
               va="top", fontsize=S.FS_LEGEND, color=S.COL_ND)

    if out_png:
        fig.savefig(out_png, dpi=S.DEFAULT_DPI)
        print("wrote", Path(out_png).name, f"({figsize[0]}x{figsize[1]} in)")
    return fig


if __name__ == "__main__":
    out = Path(__file__).resolve().parent / "panels_png" / "fig3F_reactome_cnetplot.png"
    render_panel(out_png=out)
    plt.close("all")
