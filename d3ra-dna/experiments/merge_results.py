"""Merge parallel baseline parts into results/baselines.csv and render figures + tables."""
import glob
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from d3ra.visualization import plots, reports  # noqa: E402

res = ROOT / "results"
parts = sorted(glob.glob(str(res / "baselines_part*.csv")))
if parts:
    pd.concat([pd.read_csv(p) for p in parts]).to_csv(res / "baselines.csv", index=False)
figs = plots.make_all(res, ROOT / "reports" / "figures")
reports.write_results_tables(res, ROOT / "reports" / "results_tables.md")
print("\n".join(map(str, figs)))
