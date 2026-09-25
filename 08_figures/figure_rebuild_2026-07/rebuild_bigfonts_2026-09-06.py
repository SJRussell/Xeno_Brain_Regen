"""Levin review round (2026-09-06): re-render every pipeline panel with larger fonts.

Why: the composed figures are placed at 0.5-0.8x on the manuscript page, so the 9/10 pt
deck fonts read at 4-6 pt in Word (Mike: "text labels super tiny", Figs 2/3/7). Each
generator is run with XENO_FS_SCALE chosen per panel so tick-size text lands near 7 pt at
page width (scale = 7 / (base tick pt x deck placement scale x Word placement scale)).

Usage:  python rebuild_bigfonts_2026-09-06.py render      # stage 1: panels -> panels_png/
        (deck + docx stages live in the scratchpad driver for this round)
"""
import os, subprocess, sys, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent / "figures_src"
PY = sys.executable

# (script, scale)  -- see the per-figure derivation in the 2026-09-06 session notes
JOBS = [
    (HERE / "render_fig2.py", 1.35),
    (SRC / "fig2f_cttemporal_heatmap.py", 1.15),
    (HERE / "render_fig3.py", 1.35),
    (HERE / "render_fig3b.py", 1.35),
    (SRC / "fig3c_mtorc1_le_heatmap.py", 1.35),
    (HERE / "render_fig3d.py", 1.35),
    # Fig 4 (round 2, 2026-09-06): portrait grid, every panel at its slot size placed 1:1,
    # ONE scale for the whole figure. 4A/4E/4F moved to render_fig4.py (house style).
    (HERE / "render_fig4.py", 1.15),
    (SRC / "fig4b_heme_mountain.py", 1.15),
    (SRC / "fig4c_ifn_inversion.py", 1.15),
    (SRC / "fig4d_myc_le_heatmap.py", 1.15),
    (HERE / "render_fig6.py", 1.5),
    (HERE / "render_fig7.py", 1.5),
    (HERE / "render_fig8.py", 1.5),
    (HERE / "render_fig9.py", 1.35),
    (HERE / "render_supp1.py", 1.3),
]


def render(only=None):
    fails = []
    for script, scale in JOBS:
        if only and script.name not in only:
            continue
        env = dict(os.environ, XENO_FS_SCALE=str(scale), PYTHONIOENCODING="utf-8", MPLBACKEND="Agg")
        t0 = time.time()
        r = subprocess.run([PY, str(script)], cwd=str(script.parent), env=env,
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        tail = (r.stdout.strip().splitlines() or [""])[-3:]
        status = "OK " if r.returncode == 0 else "FAIL"
        print(f"[{status}] {script.name:36s} x{scale:<4} {time.time()-t0:5.1f}s | " + " || ".join(tail))
        if r.returncode != 0:
            fails.append(script.name)
            print("   STDERR:", r.stderr.strip()[-1500:])
    print("failures:", fails)
    return fails


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "render"
    if stage == "render":
        render(only=set(sys.argv[2:]) or None)
