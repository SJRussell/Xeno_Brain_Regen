"""Interactive node dragger for Figure 3F (ECM/collagen Reactome network).

Run it locally:   python drag_fig3f.py

A window opens with the network seeded from the wide layout (or your last saved
arrangement). Drag any node with the mouse; edges and labels follow live. When it
looks right, press 's' to save the arrangement AND re-render the publication PNG
(panels_png/fig3F_reactome_cnetplot.png) through the exact figure pipeline, so the
deck panel matches what you dragged. Positions persist in
panels_png/fig3F_positions.json and are reused on every future re-render.

Keys
  drag      move the node under the cursor
  s         save positions + re-render the final publication PNG
  r         reset to the computed wide layout
  l         toggle pathway labels (declutter while arranging)
  q / close save positions and quit
"""
import matplotlib
matplotlib.use("TkAgg")
import numpy as np
import matplotlib.pyplot as plt
import fig3f_layout as L

PANEL_PNG = L.POS_JSON.parent / "fig3F_reactome_cnetplot.png"

G, df, gene_fc, path_deg, gene_nodes, path_nodes, path_sig = L.build_graph()
norm, cmap = L.make_norm_cmap(gene_fc)
pos = L.load_positions(G)
all_nodes = list(G.nodes)

# same figure aspect + axes rect + frame as the output panel, scaled up ~1.9x for a
# comfortable editing window, so what you drag is what the saved PNG shows (WYSIWYG).
WIN = 1.9
fig = plt.figure(figsize=(L.DEFAULT_FIGSIZE[0] * WIN, L.DEFAULT_FIGSIZE[1] * WIN))
ax = fig.add_axes(L.AXES_RECT)
fig.canvas.manager.set_window_title("Fig 3F  -  drag nodes, 's' to save+render")
ax.set_title("Drag nodes to arrange. Keys: s=save+render  r=reset  l=labels  q=quit",
             fontsize=10)

# --- edges (one Line2D each so we can update endpoints on drag) ---
edge_lines = []
for u, v in G.edges:
    (ln,) = ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                    color="#BBBBBB", lw=0.7, alpha=0.7, zorder=1)
    edge_lines.append((u, v, ln))

# --- nodes as two scatter collections (pathway beige, genes fc-colored) ---
path_xy = np.array([pos[n] for n in path_nodes])
gene_xy = np.array([pos[n] for n in gene_nodes])
path_sizes = [L._psize(path_deg[p]) for p in path_nodes]
gene_sizes = [L._gsize(gene_fc[g]) for g in gene_nodes]
# no rim = curated but NOT significant (adj P >= 0.05), matching render_panel
path_coll = ax.scatter(path_xy[:, 0], path_xy[:, 1], s=path_sizes,
                       facecolor="#E8D9C0",
                       edgecolor=["#9A8C6A" if path_sig[p] else "none" for p in path_nodes],
                       linewidths=[0.8 if path_sig[p] else 0.0 for p in path_nodes], zorder=2)
gene_coll = ax.scatter(gene_xy[:, 0], gene_xy[:, 1], s=gene_sizes,
                       facecolor=[cmap(norm(gene_fc[g])) for g in gene_nodes],
                       edgecolor="black", linewidths=0.9, zorder=3)
path_idx = {n: i for i, n in enumerate(path_nodes)}
gene_idx = {n: i for i, n in enumerate(gene_nodes)}

# --- labels (follow their node; pathway labels sit just above the node) ---
show_labels = [True]
labels = {}
for g in gene_nodes:
    labels[g] = ax.text(pos[g][0], pos[g][1], g, fontsize=7, fontweight="bold",
                        ha="center", va="center", zorder=4)
for p in path_nodes:
    labels[p] = ax.text(pos[p][0], pos[p][1] + 0.06, p, fontsize=6.3, color="#333333",
                        ha="center", va="bottom", zorder=4)

ax.set_xlim(*L.FRAME_X)
ax.set_ylim(*L.FRAME_Y)
ax.axis("off")

drag = {"node": None}


def _refresh_node(n):
    if n in path_idx:
        off = path_coll.get_offsets().copy(); off[path_idx[n]] = pos[n]
        path_coll.set_offsets(off)
    else:
        off = gene_coll.get_offsets().copy(); off[gene_idx[n]] = pos[n]
        gene_coll.set_offsets(off)
    dy = 0.06 if n in path_idx else 0.0
    labels[n].set_position((pos[n][0], pos[n][1] + dy))
    for u, v, ln in edge_lines:
        if u == n or v == n:
            ln.set_data([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]])


def _nearest(event):
    if event.xdata is None:
        return None
    # pick in pixel space so hit radius is uniform regardless of zoom
    p = np.array([event.x, event.y])
    best, bd = None, 1e9
    for n in all_nodes:
        px, py = ax.transData.transform(pos[n])
        d = np.hypot(px - p[0], py - p[1])
        if d < bd:
            best, bd = n, d
    return best if bd < 22 else None


def on_press(event):
    if event.inaxes != ax:
        return
    drag["node"] = _nearest(event)


def on_motion(event):
    n = drag["node"]
    if n is None or event.inaxes != ax or event.xdata is None:
        return
    pos[n] = np.array([event.xdata, event.ydata])
    _refresh_node(n)
    fig.canvas.draw_idle()


def on_release(event):
    drag["node"] = None


def on_key(event):
    if event.key == "s":
        L.save_positions(pos)
        L.render_panel(pos={k: np.asarray(v, float) for k, v in pos.items()},
                       out_png=PANEL_PNG)
        print(f"  saved positions + re-rendered {PANEL_PNG.name}")
    elif event.key == "r":
        newp = L.wide_layout(G)
        for n in all_nodes:
            pos[n] = newp[n]; _refresh_node(n)
        fig.canvas.draw_idle(); print("  reset to wide layout")
    elif event.key == "l":
        show_labels[0] = not show_labels[0]
        for p in path_nodes:
            labels[p].set_visible(show_labels[0])
        fig.canvas.draw_idle()


fig.canvas.mpl_connect("button_press_event", on_press)
fig.canvas.mpl_connect("motion_notify_event", on_motion)
fig.canvas.mpl_connect("button_release_event", on_release)
fig.canvas.mpl_connect("key_press_event", on_key)


def _on_close(event):
    L.save_positions(pos)
    print(f"  positions saved to {L.POS_JSON.name}")


fig.canvas.mpl_connect("close_event", _on_close)

print(__doc__)
plt.show()
