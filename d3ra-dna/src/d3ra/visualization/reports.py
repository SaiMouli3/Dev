"""Markdown rendering of the result CSVs (tables + figure index)."""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

try:                                   # pandas.to_markdown needs tabulate
    import tabulate  # noqa: F401
    _HAS_TABULATE = True
except ImportError:                    # pragma: no cover
    _HAS_TABULATE = False

SECTIONS = [
    ("baselines.csv", "Baselines A-F per scenario",
     ["scenario", "design", "variant", "B_ref", "cap@0.99", "mean_spend@0.99", "cap@0.999", "mean_spend@0.999",
      "p_exact@ref", "mean_cost_success", "mean_rounds", "mean_containers_accessed", "redundancy_overhead", "runtime_s"],
     lambda d: d[~d["design"].astype(str).str.startswith("cmp:")] if "design" in d else d),
    ("baselines.csv", "Design comparisons (bootstrap mean cost difference, Cliff's delta, Holm-adjusted McNemar p)",
     ["scenario", "design", "mean_cost_diff", "mean_cost_diff_lo", "mean_cost_diff_hi", "cliffs_delta", "p_holm"],
     lambda d: d[d["design"].astype(str).str.startswith("cmp:")] if "design" in d else d.iloc[:0]),
    ("h1_model_accuracy.csv", "H1 - predictive accuracy of channel models",
     ["model", "mae_p", "max_abs_err", "B99_pred", "B99_mc", "B99_rel_err"],
     lambda d: d[d["experiment"] == "H1-summary"] if "experiment" in d else d),
    ("h5_robustness.csv", "H5 - robustness to mis-estimated channel parameters",
     ["scenario", "controller", "cap@0.99", "mean_spend@0.99", "cap@0.999", "mean_spend@0.999", "mean_cost_success",
      "p_success_uncapped"], None),
    ("ablations.csv", "Ablations",
     ["scenario", "ablation", "variant", "cap@0.99", "mean_spend@0.99", "cap@0.999", "mean_spend@0.999",
      "mean_cost_success"], None),
    ("redundancy.csv", "Redundancy allocation (uniform vs risk-adaptive parity)",
     ["experiment", "reads_per_fragment", "method", "parity_total", "overhead", "parity_min", "parity_max",
      "p_file_pred", "p_file_mc"], None),
    ("replication_vs_parity.csv", "Replication vs RS parity under shared bias",
     ["shared_variance_fraction", "scheme", "stored_molecules", "cap@0.99", "mean_cost_success",
      "p_success_uncapped"], None),
    ("placement_opt.csv", "Placement optimisers",
     ["scenario", "stripes", "method", "objective", "gap_vs_ilp", "time"], None),
    ("scalability.csv", "Scalability",
     ["containers", "stripes", "fragments", "placement_s", "predictor_bisection_s", "mc_static_2000_s", "adaptive_s",
      "peak_mem_mb", "B99_pred", "B99_mc"], None),
]


def _fmt(v) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        if np.isinf(v):
            return "inf" if v > 0 else "-inf"
        a = abs(v)
        if a == 0:
            return "0"
        if a >= 1000:
            return f"{v:,.0f}"
        if a >= 10:
            return f"{v:.1f}"
        if a >= 0.01:
            return f"{v:.3f}"
        return f"{v:.2e}"
    return str(v).replace("|", "\\|")


def _add_ci(d: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Render `x` with its `x_ci_lo/x_ci_hi` columns as 'x [lo, hi]'."""
    out = pd.DataFrame(index=d.index)
    for c in cols:
        if c not in d.columns:
            continue
        lo, hi = f"{c}_ci_lo", f"{c}_ci_hi"
        if lo in d.columns and hi in d.columns:
            out[c] = [f"{_fmt(v)} [{_fmt(l)}, {_fmt(h)}]" if _fmt(v) else ""
                      for v, l, h in zip(d[c], d[lo], d[hi])]
        else:
            out[c] = [_fmt(v) for v in d[c]]
    return out


def df_to_markdown(df: pd.DataFrame) -> str:
    if df.empty:
        return "_(no rows)_\n"
    if _HAS_TABULATE:
        return df.to_markdown(index=False) + "\n"
    cols = [str(c) for c in df.columns]
    rows = [[str(x) for x in r] for r in df.itertuples(index=False)]
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(lines) + "\n"


def write_results_tables(results_dir="results", out_md="results/RESULTS.md") -> str:
    """Render the key CSVs in ``results_dir`` to one markdown file; returns its path."""
    parts = ["# D3RA results tables\n",
             "Costs are in cost units (reads x cost per read); `cap@t` is the budget cap reaching "
             "reliability t, `mean_spend@t` the mean spend under that cap. Brackets: 95% CI.\n"]
    for csv, title, cols, filt in SECTIONS:
        path = os.path.join(results_dir, csv)
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            continue
        try:
            d = pd.read_csv(path)
        except Exception:
            continue
        if filt is not None:
            try:
                d = filt(d)
            except Exception:
                pass
        if d.empty:
            continue
        t = _add_ci(d, [c for c in cols if c in d.columns])
        if t.shape[1] == 0:
            continue
        parts.append(f"\n## {title}\n\n_source: `{csv}`_\n\n")
        parts.append(df_to_markdown(t.dropna(axis=1, how="all")))
    out = Path(out_md)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(parts))
    return str(out)


def figure_index(figures, relative_to=None) -> str:
    """Markdown list (with embedded images) of figure paths."""
    lines = ["## Figures\n"]
    for f in figures:
        p = os.path.relpath(f, relative_to) if relative_to else str(f)
        name = Path(f).stem.replace("_", " ")
        lines.append(f"- **{name}**: `{p}`\n\n  ![{name}]({p})\n")
    return "\n".join(lines) + "\n"
