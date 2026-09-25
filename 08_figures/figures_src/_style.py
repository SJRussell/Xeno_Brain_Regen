"""Shared style for all Murugan/Russell brain RNA-seq figures.

Every figure script imports from here. To globally restyle:
- Toggle SHOW_SUBTITLE / SHOW_CAPTIONS
- Change FONT_FAMILY, FONT_SIZE_*
- Change PALETTE_* values

One change here propagates to every figure on re-render.

2026-07 DECK CONVENTIONS (rebuild consolidation):
- Font is Arial (deck house style). Falls back to DejaVu Sans if Arial is absent.
- Treatment palette is MDT (coral) / ND (blue). "CT" is the deprecated old name for
  the multi-drug-treatment group; COL_CT is kept ONLY as a back-compat alias == COL_MDT
  so legacy figXX scripts don't crash. New scripts use COL_MDT and label it "MDT".
- Timepoint annotation palette TP_COL (24 HPT orange / 7 DPT green) and the heatmap
  z-range ZMAX_HEATMAP (+-4) live here so every heatmap matches the canonical renders.

Render loop:
  for f in figures_src/fig*.py: python "$f"
"""
import os
import matplotlib as mpl

# ---- global toggles ----
SHOW_SUBTITLE = False        # Stewart 2026-06-05: drop subtitles globally
SHOW_FIGURE_CAPTIONS = False
TITLE_BOLD = True
DEFAULT_DPI = 300            # 2026-07: print-quality panels for the rebuild deck

# ---- font ----
# Arial is the deck house style; keep DejaVu Sans as the fallback so renders never fail.
# XENO_FS_SCALE (env, default 1) multiplies every point size (2026-09-06, Levin review round:
# composed figures land at 0.5-0.8x on the manuscript page, so deck sizes read at 4-6 pt).
# DECK_FS_TICK baseline raised 6 -> 7 at the same time (6 pt sample/gene labels were the
# smallest text in the deck).
FS_SCALE = float(os.environ.get("XENO_FS_SCALE", "1"))
FONT_FAMILY = "Arial"
FONT_FALLBACK = ["Arial", "Helvetica", "DejaVu Sans"]
FONT_SIZE_TITLE   = 14 * FS_SCALE   # legacy panels; rebuild heatmaps set ~10 pt locally
FONT_SIZE_AXIS    = 11 * FS_SCALE
FONT_SIZE_TICK    = 10 * FS_SCALE
FONT_SIZE_LEGEND  = 10 * FS_SCALE
FONT_SIZE_LABEL   = 11 * FS_SCALE

# rebuild deck sizes (used by the 2026-07 per-figure heatmap generators)
DECK_FS_TITLE   = 10 * FS_SCALE
DECK_FS_LEGEND  = 8 * FS_SCALE
DECK_FS_TICK    = 7 * FS_SCALE
DECK_FS_AXIS    = 9 * FS_SCALE   # standard (non-dense) tick/axis size; == figure_rebuild FS_TICK/FS_AXIS

# ---- colour palettes (kept consistent across figures) ----
# Treatments (2026-07 approved rebuild palette; matches the _nibi_refs R panels)
COL_ND  = "#3B7DBF"   # no-drug control, blue
COL_MDT = "#D9614C"   # multi-drug treatment, coral
COL_CT  = COL_MDT     # DEPRECATED alias: "CT" == MDT. Do not use in new scripts.
COL_UNINJURED = "#7F7F7F"

# Timepoint annotation bars (heatmaps)
TP_COL = {"24HPA": "#F0A860", "7DPA": "#4C9E7C"}   # 24 HPT orange, 7 DPT green

# Direction (volcanoes)
COL_UP   = COL_MDT     # up = coral
COL_DOWN = COL_ND      # down = blue
COL_NS   = "#BBBBBB"

# Diverging colormap (NES, log2FC, z-score heatmaps): RdBu_r
CMAP_DIVERGING = "RdBu_r"
CMAP_SEQUENTIAL = "YlOrRd"
CMAP_ENRICHMENT = "Reds"

# Heatmap z-score display range (canonical: RdBu_r clipped +-4, ticks at +-2/+-4)
ZMAX_HEATMAP = 4.0

# Title text colour
COL_TEXT_TITLE = "#333333"

# ---- rcParams ----
def apply():
    mpl.rcParams["font.family"]          = FONT_FAMILY
    mpl.rcParams["font.sans-serif"]      = FONT_FALLBACK
    mpl.rcParams["axes.unicode_minus"]   = False
    mpl.rcParams["font.size"]            = FONT_SIZE_AXIS
    mpl.rcParams["axes.titlesize"]       = FONT_SIZE_TITLE
    mpl.rcParams["axes.titleweight"]     = "bold" if TITLE_BOLD else "normal"
    mpl.rcParams["axes.titlecolor"]      = COL_TEXT_TITLE
    mpl.rcParams["axes.labelsize"]       = FONT_SIZE_LABEL
    mpl.rcParams["xtick.labelsize"]      = FONT_SIZE_TICK
    mpl.rcParams["ytick.labelsize"]      = FONT_SIZE_TICK
    mpl.rcParams["legend.fontsize"]      = FONT_SIZE_LEGEND
    mpl.rcParams["figure.dpi"]           = DEFAULT_DPI
    mpl.rcParams["savefig.dpi"]          = DEFAULT_DPI
    mpl.rcParams["savefig.bbox"]         = "tight"
    mpl.rcParams["axes.spines.top"]      = False
    mpl.rcParams["axes.spines.right"]    = False

apply()
