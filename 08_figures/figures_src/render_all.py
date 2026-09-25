"""Render every fig*.py in figures_src/ in alphabetical order.

Run from repo root or from figures_src/:
    py figures_src/render_all.py
"""
import sys, subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

scripts = sorted(HERE.glob("fig*.py"))
print(f"Found {len(scripts)} figure scripts.\n")

failed = []
for s in scripts:
    rel = s.relative_to(ROOT)
    print(f"  [{s.name:50}] ", end="", flush=True)
    r = subprocess.run([sys.executable, str(s)], capture_output=True, text=True)
    if r.returncode == 0:
        out_line = (r.stdout or "").strip().splitlines()[-1] if r.stdout else "ok"
        print(out_line)
    else:
        print(f"FAILED ({r.returncode})")
        if r.stderr: print("    stderr:", r.stderr.strip()[:400])
        failed.append(s.name)

print(f"\n{len(scripts) - len(failed)}/{len(scripts)} rendered.")
if failed:
    print("Failed:")
    for f in failed: print(f"  - {f}")
    sys.exit(1)
