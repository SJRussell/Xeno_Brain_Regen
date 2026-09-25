"""Shared style for the 2026-07 figure rebuild (NPJ Regen Med).

Global rules from Stewart 2026-07-07:
  - Arial everywhere. Figure/panel titles ~10 pt. No subtitles.
  - "MDT" (multi-drug treatment) everywhere, never "CT" (cocktail).
  - Panel letters are added as MOVABLE pptx textboxes (Arial 12 pt caps),
    NOT baked into the PNGs, so Stewart can nudge them.
  - Panels render at an EXACT figsize (no bbox='tight') so each fills a clean
    rectangle in the deck grid. Font points are therefore consistent panel-to-panel
    when panels are placed at 1:1.

One change here propagates to every figure on re-render.
"""
import os
import matplotlib as mpl

# ---- toggles ----
SHOW_SUBTITLE = False
DEFAULT_DPI = 300

# ---- font ----
# XENO_FS_SCALE (env, default 1) multiplies every point size below. Added 2026-09-06 for the
# Levin review round: the composed figures are placed at 0.5-0.8x in the manuscript page, so
# the 9/10 pt deck sizes landed at 4-6 pt in Word. The driver (figure_rebuild_2026-07/
# rebuild_bigfonts_2026-09-06.py) sets a per-panel scale so tick text lands near 7 pt at page
# width. Everything stays relative (FS_TICK - 1 etc.), so layouts scale together.
FS_SCALE = float(os.environ.get("XENO_FS_SCALE", "1"))
FONT_FAMILY = "Arial"
FS_TITLE  = 10 * FS_SCALE
FS_AXIS   = 9 * FS_SCALE
FS_TICK   = 9 * FS_SCALE
FS_LEGEND = 8 * FS_SCALE
FS_ANNOT  = 8 * FS_SCALE

# ---- treatment palette (matches the approved deck / _nibi_refs R panels) ----
COL_ND  = "#3B7DBF"   # blue  (no-drug control)
COL_MDT = "#D9614C"   # coral (multi-drug treatment; was "CT")
COL_UNINJURED = "#7F7F7F"

# ---- direction (volcanoes / heatmaps) ----
COL_UP   = "#D9614C"
COL_DOWN = "#3B7DBF"
COL_NS   = "#BBBBBB"

CMAP_DIVERGING  = "RdBu_r"
CMAP_SEQUENTIAL = "YlOrRd"
CMAP_ENRICHMENT = "Reds"

COL_HULL_EDGE = "#404040"
COL_TEXT      = "#222222"

# ---- treatment / timepoint display labels ----
# data uses CT/ND and 24HPA/7DPA/11HPA; the figures show MDT/ND and HPT/DPT.
TREAT_LABEL = {"CT": "MDT", "ND": "ND", "U": "Uninjured"}
TREAT_COL   = {"CT": COL_MDT, "ND": COL_ND, "U": COL_UNINJURED}
TIME_LABEL  = {"11HPA": "11 HPT", "24HPA": "24 HPT", "7DPA": "7 DPT"}
TIME_SHAPE  = {"11HPA": "o", "24HPA": "^", "7DPA": "s"}   # blastema (3 timepoints)
TIME_SHAPE_BRAIN = {"24HPA": "o", "7DPA": "^"}            # brain (2 timepoints)


def apply():
    mpl.rcParams["font.family"]       = FONT_FAMILY
    mpl.rcParams["font.size"]         = FS_AXIS
    mpl.rcParams["axes.titlesize"]    = FS_TITLE
    mpl.rcParams["axes.titleweight"]  = "normal"
    mpl.rcParams["axes.titlecolor"]   = COL_TEXT
    mpl.rcParams["axes.labelsize"]    = FS_AXIS
    mpl.rcParams["xtick.labelsize"]   = FS_TICK
    mpl.rcParams["ytick.labelsize"]   = FS_TICK
    mpl.rcParams["legend.fontsize"]   = FS_LEGEND
    mpl.rcParams["figure.dpi"]        = DEFAULT_DPI
    mpl.rcParams["savefig.dpi"]       = DEFAULT_DPI
    mpl.rcParams["axes.spines.top"]   = False
    mpl.rcParams["axes.spines.right"] = False
    mpl.rcParams["axes.edgecolor"]    = "black"
    mpl.rcParams["pdf.fonttype"]      = 42   # editable text in Illustrator
    mpl.rcParams["ps.fonttype"]       = 42

apply()
