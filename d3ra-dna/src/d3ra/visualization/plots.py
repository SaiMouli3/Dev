"""Publication figures for the D3RA experiments.

Every ``plot_*`` function takes a pandas DataFrame (as written to ``results/``
by ``experiments/run_all.py``) and an output path, and returns the written path
or ``None`` when the frame lacks the data it needs (it never raises on missing
columns). ``make_all`` renders every figure for which a CSV exists.

Cost units: total sequencing spend = sum over containers of reads x cost per read.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# Okabe-Ito colour-blind-safe palette
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#000000"]
DESIGN_ORDER = ["A_uniform_fixed", "B_reliability_fixed", "C_uniform_adaptive", "D_reliability_adaptive",
                "E_joint_adaptive", "F_oracle"]
DESIGN_LABEL = {"A_uniform_fixed": "A uniform/fixed", "B_reliability_fixed": "B reliab./fixed",
                "C_uniform_adaptive": "C uniform/adaptive", "D_reliability_adaptive": "D reliab./adaptive",
                "E_joint_adaptive": "E joint/adaptive", "F_oracle": "F oracle"}
COST_UNITS = "cost units (reads x cost/read)"
HATCHES = ["", "//", "", "\\\\", "", ".."]

STYLE = {
    "figure.dpi": 110, "figure.constrained_layout.use": True, "savefig.dpi": 160, "savefig.bbox": "tight",
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "legend.fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.alpha": 0.3, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "axes.prop_cycle": matplotlib.cycler(color=PALETTE),
}


def _colour(i: int) -> str:
    return PALETTE[i % len(PALETTE)]


def _has(df: pd.DataFrame | None, *cols: str) -> bool:
    return df is not None and not df.empty and all(c in df.columns for c in cols)


def _num(s) -> np.ndarray:
    """Numeric array with +-inf mapped to nan (unplottable)."""
    a = pd.to_numeric(pd.Series(s), errors="coerce").to_numpy(dtype=float, copy=True)
    a[~np.isfinite(a)] = np.nan
    return a


def _yerr(df: pd.DataFrame, col: str):
    lo, hi = f"{col}_ci_lo", f"{col}_ci_hi"
    if lo not in df.columns or hi not in df.columns:
        return None
    v, l, h = _num(df[col]), _num(df[lo]), _num(df[hi])
    e = np.vstack([np.clip(v - l, 0, None), np.clip(h - v, 0, None)])
    return np.nan_to_num(e, nan=0.0)


def _save(fig, path) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    return str(path)


def _grouped_bars(ax, df: pd.DataFrame, group_col: str, series_col: str, value_col: str,
                  series_order=None, labels=None, err=True, hatch=False):
    groups = list(dict.fromkeys(df[group_col]))
    series = [s for s in (series_order or []) if s in set(df[series_col])]
    series += [s for s in dict.fromkeys(df[series_col]) if s not in series]
    n = max(len(series), 1)
    width = 0.8 / n
    x = np.arange(len(groups))
    for i, s in enumerate(series):
        sub = df[df[series_col] == s].set_index(group_col).reindex(groups)
        vals = _num(sub[value_col])
        ye = _yerr(sub.reset_index(), value_col) if err else None
        ax.bar(x + (i - (n - 1) / 2) * width, np.nan_to_num(vals), width, yerr=ye, capsize=1.5,
               color=_colour(i), edgecolor="white" if not hatch else "black", linewidth=0.4,
               hatch=HATCHES[i % len(HATCHES)] if hatch else None,
               label=(labels or {}).get(s, s), error_kw={"elinewidth": 0.7, "ecolor": "#333333"})
        for xi, v in zip(x + (i - (n - 1) / 2) * width, vals):
            if np.isnan(v):          # unreachable target (inf) -> marker
                ax.plot(xi, 0, marker="x", color=_colour(i), ms=4)
    ax.set_xticks(x)
    ax.set_xticklabels([str(g).replace("_", "\n", 1) for g in groups], rotation=0)
    return series


# --------------------------------------------------------------------------- baselines

def plot_baselines(df: pd.DataFrame, out) -> str | None:
    if not _has(df, "scenario", "design"):
        return None
    d = df[df["design"].isin(DESIGN_ORDER)].copy()
    panels = [(c, t) for c, t in [("mean_spend@0.99", "Mean spend at 0.99 reliability cap"),
                                  ("cap@0.99", "Budget cap for 0.99 reliability"),
                                  ("p_exact@ref", "P(exact recovery) at reference budget")] if c in d.columns]
    if d.empty or not panels:
        return None
    with plt.rc_context(STYLE):
        nsc = d["scenario"].nunique()
        fig, axes = plt.subplots(len(panels), 1, figsize=(max(7, 1.1 * nsc + 2), 2.6 * len(panels)), squeeze=False)
        for ax, (col, title) in zip(axes[:, 0], panels):
            _grouped_bars(ax, d, "scenario", "design", col, DESIGN_ORDER, DESIGN_LABEL, hatch=True)
            ax.set_title(title)
            if col.startswith("p_exact"):
                ax.set_ylabel("P(exact recovery)")
                ax.set_ylim(0, 1.05)
            else:
                ax.set_ylabel(COST_UNITS)
            ax.tick_params(axis="x", labelsize=7)
        h, l = axes[0, 0].get_legend_handles_labels()
        fig.legend(h, l, ncol=6, loc="outside upper center", frameon=False)
        fig.supxlabel("x marker = target unreachable (infinite budget)", fontsize=7, color="#555555")
        return _save(fig, out)


def plot_baseline_comparisons(df: pd.DataFrame, out) -> str | None:
    """Forest plot of paired/unpaired mean-cost differences (cmp: rows)."""
    if not _has(df, "design", "mean_cost_diff"):
        return None
    d = df[df["design"].astype(str).str.startswith("cmp:")].copy()
    if d.empty:
        return None
    d["label"] = d["scenario"].astype(str) + "  " + d["design"].str[4:].str.replace("_", " ")
    v = _num(d["mean_cost_diff"])
    lo = _num(d.get("mean_cost_diff_lo", d["mean_cost_diff"]))
    hi = _num(d.get("mean_cost_diff_hi", d["mean_cost_diff"]))
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(7.5, 0.22 * len(d) + 1.2))
        y = np.arange(len(d))[::-1]
        comps = list(dict.fromkeys(d["design"]))
        cols = [_colour(comps.index(c)) for c in d["design"]]
        ax.errorbar(v, y, xerr=np.nan_to_num(np.vstack([v - lo, hi - v]).clip(0)), fmt="none", ecolor="#555", lw=0.8)
        ax.scatter(v, y, c=cols, s=18, zorder=3)
        if "p_holm" in d.columns:
            for yi, vi, p in zip(y, v, _num(d["p_holm"])):
                if np.isfinite(p) and p < 0.05:
                    ax.annotate("*", (vi, yi), xytext=(4, -3), textcoords="offset points", fontsize=8)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(d["label"], fontsize=6.5)
        ax.set_xlabel(f"Mean capped-cost difference (first - second), {COST_UNITS}")
        ax.set_title("Design comparisons (95% bootstrap CI; * Holm p<0.05 McNemar)")
        return _save(fig, out)


# --------------------------------------------------------------------------- H1

def plot_h1(df: pd.DataFrame, out) -> str | None:
    if not _has(df, "model", "p_mc", "p_pred"):
        return None
    d = df[df["experiment"] == "H1"] if "experiment" in df.columns else df
    d = d.dropna(subset=["p_mc", "p_pred"])
    if d.empty:
        return None
    summ = df[df["experiment"] == "H1-summary"].set_index("model") if "experiment" in df.columns else pd.DataFrame()
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(10, 4))
        ax = axes[0]
        ax.plot([0, 1], [0, 1], color="#888", ls="--", lw=0.9, label="y = x")
        markers = "osD^v<>"
        for i, (m, g) in enumerate(d.groupby("model", sort=False)):
            g = g.sort_values("p_mc")
            lab = m
            if "mae_p" in summ.columns and m in summ.index:
                lab += f" (MAE {float(summ.loc[m, 'mae_p']):.3f})"
            ax.plot(g["p_mc"], g["p_pred"], marker=markers[i % len(markers)], ms=4, color=_colour(i), lw=1.2, label=lab)
        ax.set_xlabel("Monte-Carlo P(file recovered)")
        ax.set_ylabel("Predicted P(file recovered)")
        ax.set_xlim(-0.02, 1.02)
        ax.set_ylim(-0.02, 1.02)
        ax.set_aspect("equal")
        ax.set_title("H1: predictor calibration")
        ax.legend(frameon=False, fontsize=7)
        ax = axes[1]
        if "budget" in d.columns:
            mc = d.drop_duplicates("budget").sort_values("budget")
            ax.plot(mc["budget"], mc["p_mc"], color="black", lw=2, label="Monte-Carlo (truth)")
            for i, (m, g) in enumerate(d.groupby("model", sort=False)):
                g = g.sort_values("budget")
                ax.plot(g["budget"], g["p_pred"], marker=markers[i % len(markers)], ms=3.5, color=_colour(i),
                        lw=1, ls="--", label=m)
            ax.set_xlabel("Total budget (reads)")
            ax.set_ylabel("P(file recovered)")
            ax.set_title("Recovery curve: predicted vs MC")
            ax.legend(frameon=False, fontsize=7)
        else:
            ax.axis("off")
        return _save(fig, out)


# --------------------------------------------------------------------------- H5

def plot_h5(df: pd.DataFrame, out) -> str | None:
    if not _has(df, "scenario", "controller"):
        return None
    panels = [c for c in ("mean_spend@0.99", "cap@0.99") if c in df.columns]
    if not panels:
        return None
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, len(panels), figsize=(5.5 * len(panels), 3.8), squeeze=False)
        for ax, col in zip(axes[0], panels):
            _grouped_bars(ax, df, "scenario", "controller", col)
            ax.set_ylabel(COST_UNITS)
            ax.set_title({"mean_spend@0.99": "H5: mean spend (0.99 cap)", "cap@0.99": "H5: budget cap for 0.99"}[col])
        axes[0, 0].legend(frameon=False, fontsize=7, loc="upper left")
        return _save(fig, out)


# --------------------------------------------------------------------------- ablations

def plot_ablations(df: pd.DataFrame, out, metric: str = "mean_spend@0.99") -> str | None:
    if not _has(df, "ablation") or metric not in df.columns:
        return None
    ref = df[df["ablation"] == "full (E)"]
    if ref.empty:
        return None
    figs = []
    with plt.rc_context(STYLE):
        scen = list(dict.fromkeys(df["scenario"])) if "scenario" in df.columns else [None]
        fig, axes = plt.subplots(len(scen), 1, figsize=(8, 0.35 * len(df) / max(len(scen), 1) + 1.5) if len(scen) == 1
                                 else (8, 3.2 * len(scen)), squeeze=False)
        for ax, s in zip(axes[:, 0], scen):
            d = df if s is None else df[df["scenario"] == s]
            r = d[d["ablation"] == "full (E)"]
            if r.empty:
                ax.axis("off")
                continue
            base = _num(r[metric])[0]
            d = d[d["ablation"] != "full (E)"]
            rel = (_num(d[metric]) / base - 1) * 100
            y = np.arange(len(d))[::-1]
            cols = [_colour(3) if v > 0 else _colour(2) for v in np.nan_to_num(rel)]
            ax.barh(y, np.nan_to_num(rel), color=cols, edgecolor="white")
            for yi, v in zip(y, rel):
                ax.annotate("unreachable (inf)" if np.isnan(v) else f"{v:+.1f}%", (0 if np.isnan(v) else v, yi),
                            xytext=(3 if (np.isnan(v) or v >= 0) else -3, 0), textcoords="offset points",
                            va="center", ha="left" if (np.isnan(v) or v >= 0) else "right", fontsize=7)
            ax.axvline(0, color="black", lw=0.8)
            ax.margins(x=0.18)
            ax.set_yticks(y)
            ax.set_yticklabels(d["ablation"], fontsize=7)
            ax.set_xlabel(f"Change in {metric} vs full (E) = {base:,.0f} {COST_UNITS} [%]")
            ax.set_title(f"Ablations{'' if s is None else ' - ' + str(s)} (positive = worse)")
            figs.append(ax)
        return _save(fig, out)


# --------------------------------------------------------------------------- redundancy

def plot_redundancy(df: pd.DataFrame, out) -> str | None:
    if not _has(df, "method", "reads_per_fragment"):
        return None
    pred = df[df["experiment"] == "redundancy"] if "experiment" in df.columns else df
    if "overhead" not in pred.columns or pred.dropna(subset=["overhead"]).empty:
        return None
    pred = pred.dropna(subset=["overhead"])
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
        _grouped_bars(axes[0], pred.assign(rpf=pred["reads_per_fragment"].astype(str) + " reads/frag"),
                      "rpf", "method", "overhead", err=False)
        axes[0].set_ylabel("Parity overhead (parity / data fragments)")
        axes[0].set_title("Storage redundancy at equal reliability target")
        axes[0].legend(frameon=False)
        ax = axes[1]
        cols = [c for c in ("p_file_pred", "p_file_mc") if c in df.columns]
        if cols:
            rows = []
            for (rpf, meth), g in df.groupby(["reads_per_fragment", "method"]):
                for c in cols:
                    v = _num(g[c])
                    v = v[np.isfinite(v)]
                    if len(v):
                        rows.append({"grp": f"{rpf:g} r/f\n{meth}", "kind": "predicted" if c == "p_file_pred" else "Monte-Carlo",
                                     "p": v[0]})
            if rows:
                _grouped_bars(ax, pd.DataFrame(rows), "grp", "kind", "p", err=False)
                lo = min(r["p"] for r in rows)
                ax.set_ylim(max(0, lo - 0.01), 1.0005)
                ax.legend(frameon=False, loc="lower right", ncol=2)
            ax.set_ylabel("P(file recovered)")
            ax.set_title("Predicted vs MC reliability")
            ax.tick_params(axis="x", labelsize=7)
        else:
            ax.axis("off")
        return _save(fig, out)


def plot_replication_vs_parity(df: pd.DataFrame, out) -> str | None:
    if not _has(df, "shared_variance_fraction", "scheme"):
        return None
    cols = [c for c in ("cap@0.99", "mean_cost_success") if c in df.columns]
    if not cols:
        return None
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, len(cols), figsize=(4.6 * len(cols), 3.4), squeeze=False)
        for ax, col in zip(axes[0], cols):
            for i, (s, g) in enumerate(df.groupby("scheme", sort=False)):
                g = g.sort_values("shared_variance_fraction")
                x, v = _num(g["shared_variance_fraction"]), _num(g[col])
                ax.plot(x, v, marker="os"[i % 2], color=_colour(i), label=s)
                e = _yerr(g, col)
                if e is not None:
                    ax.fill_between(x, v - e[0], v + e[1], color=_colour(i), alpha=0.18, lw=0)
            ax.set_xlabel("Shared (sequence-intrinsic) share of abundance log-variance")
            ax.set_ylabel("Reads (unit cost/read)")
            ax.set_title({"cap@0.99": "Budget for 0.99 reliability", "mean_cost_success": "Mean reads to recover"}[col])
        axes[0, 0].legend(frameon=False)
        fig.suptitle("Replication vs RS parity at equal (2x) stored molecules", fontsize=10)
        return _save(fig, out)


# --------------------------------------------------------------------------- placement / scalability

def plot_placement_gap(df: pd.DataFrame, out) -> str | None:
    if not _has(df, "scenario", "method", "gap_vs_ilp"):
        return None
    d = df.copy()
    d["gap_pct"] = _num(d["gap_vs_ilp"]) * 100
    with plt.rc_context(STYLE):
        has_t = "time" in d.columns
        fig, axes = plt.subplots(1, 2 if has_t else 1, figsize=(10 if has_t else 5.5, 3.6), squeeze=False)
        _grouped_bars(axes[0, 0], d, "scenario", "method", "gap_pct", err=False)
        axes[0, 0].set_ylabel("Objective gap vs ILP [%]")
        axes[0, 0].set_title("Placement optimiser optimality gap")
        axes[0, 0].axhline(0, color="black", lw=0.8)
        axes[0, 0].legend(frameon=False, fontsize=7)
        if has_t:
            ax = axes[0, 1]
            _grouped_bars(ax, d.assign(time=np.clip(_num(d["time"]), 1e-4, None)), "scenario", "method", "time", err=False)
            ax.set_yscale("log")
            ax.set_ylabel("Wall-clock time [s]")
            ax.set_title("Optimiser runtime")
        return _save(fig, out)


TIME_LABEL = {"placement_s": "placement", "predictor_bisection_s": "predictor + bisection (B99)",
              "mc_static_2000_s": "static MC (2000 trials)", "adaptive_s": "adaptive controller"}


def plot_scalability(df: pd.DataFrame, out) -> str | None:
    if not _has(df, "fragments"):
        return None
    tcols = [c for c in ("placement_s", "predictor_bisection_s", "mc_static_2000_s", "adaptive_s") if c in df.columns]
    if not tcols:
        return None
    d = df.sort_values("fragments")
    x = _num(d["fragments"])
    with plt.rc_context(STYLE):
        extra = [c for c in ("B99_pred", "B99_mc") if c in d.columns]
        fig, axes = plt.subplots(1, 2 if extra else 1, figsize=(10 if extra else 5.5, 3.8), squeeze=False)
        ax = axes[0, 0]
        for i, c in enumerate(tcols):
            ax.plot(x, np.clip(_num(d[c]), 1e-4, None), marker="osD^"[i], color=_colour(i), label=TIME_LABEL.get(c, c))
        if len(x) >= 2:   # linear-scaling guide
            ref = np.nanmax(np.clip(_num(d[tcols[0]]), 1e-4, None))
            ax.plot(x, ref * x / x.max(), color="#999", ls=":", lw=1, label="O(n) guide")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Encoded fragments")
        ax.set_ylabel("Runtime [s]")
        ax.set_title("Scalability (log-log)")
        ax.legend(frameon=False, fontsize=7)
        if "containers" in d.columns:
            for xi, m, s in zip(x, d["containers"], d.get("stripes", d["containers"])):
                ax.annotate(f"M={m}", (xi, ax.get_ylim()[0]), xytext=(0, 2), textcoords="offset points",
                            fontsize=6, ha="center", color="#555")
        if extra:
            ax = axes[0, 1]
            for i, c in enumerate(extra):
                ax.plot(x, _num(d[c]), marker="os"[i], color=_colour(i + 4), label={"B99_pred": "predicted", "B99_mc": "Monte-Carlo"}[c])
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_xlabel("Encoded fragments")
            ax.set_ylabel("Budget for 0.99 reliability (reads)")
            ax.set_title("Predictor vs MC at scale")
            ax.legend(frameon=False)
        return _save(fig, out)


# --------------------------------------------------------------------------- recovery curve helper

def plot_recovery_curve(layout, ch, out, weights=None, budgets=None, mc_samples=None, label="predicted",
                        ax=None, n_points: int = 40) -> str | None:
    """Predicted P(file recovered) vs total budget in cost units for a static
    allocation ``weights`` (read fractions; default uniform over used
    containers). ``mc_samples`` (per-trial costs) adds an empirical ECDF."""
    from ..channel.surrogate import Predictor
    from ..retrieval.fixed_policy import uniform_weights

    w = np.asarray(uniform_weights(layout) if weights is None else weights, float)
    w = w / w.sum()
    cpr = float(ch.cost @ w)                        # cost per read under w
    pred = Predictor(layout, ch)
    if budgets is None:
        hi = 30.0 * layout.n * cpr
        if mc_samples is not None and np.isfinite(mc_samples).any():
            fin = np.asarray(mc_samples, float)
            hi = 1.3 * np.quantile(fin[np.isfinite(fin)], 0.999)
        budgets = np.linspace(0, hi, n_points)[1:]
    budgets = np.asarray(budgets, float)
    p = np.array([pred.p_file(B * w / cpr) for B in budgets])
    own = ax is None
    with plt.rc_context(STYLE):
        if own:
            fig, ax = plt.subplots(figsize=(5.5, 3.6))
        ax.plot(budgets, p, color=_colour(0), lw=1.6, label=label)
        if mc_samples is not None:
            s = np.sort(np.asarray(mc_samples, float))
            ax.step(s, np.arange(1, len(s) + 1) / len(s), where="post", color=_colour(1), lw=1.2, label="Monte-Carlo ECDF")
        for t in (0.99, 0.999):
            ax.axhline(t, color="#999", ls=":", lw=0.8)
        ax.set_xlabel(f"Total budget, {COST_UNITS}")
        ax.set_ylabel("P(file recovered)")
        ax.set_ylim(0, 1.01)
        ax.legend(frameon=False)
        if own:
            ax.set_title("Recovery probability vs budget")
            return _save(fig, out)
    return None


# --------------------------------------------------------------------------- driver

FIGURES = [
    ("baselines.csv", plot_baselines, "baselines.png"),
    ("baselines.csv", plot_baseline_comparisons, "baseline_comparisons.png"),
    ("h1_model_accuracy.csv", plot_h1, "h1_model_accuracy.png"),
    ("h5_robustness.csv", plot_h5, "h5_robustness.png"),
    ("ablations.csv", plot_ablations, "ablations.png"),
    ("redundancy.csv", plot_redundancy, "redundancy.png"),
    ("replication_vs_parity.csv", plot_replication_vs_parity, "replication_vs_parity.png"),
    ("placement_opt.csv", plot_placement_gap, "placement_gap.png"),
    ("scalability.csv", plot_scalability, "scalability.png"),
]


def load_csv(path) -> pd.DataFrame | None:
    try:
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            return None
        return pd.read_csv(path)
    except Exception:
        return None


def make_all(results_dir="results", figures_dir="figures") -> list[str]:
    """Render every figure whose CSV exists; returns the list of written paths."""
    written = []
    for csv, fn, png in FIGURES:
        df = load_csv(os.path.join(results_dir, csv))
        if df is None:
            continue
        try:
            p = fn(df, os.path.join(figures_dir, png))
        except Exception as e:          # one bad CSV must not stop the rest
            print(f"[plots] {fn.__name__} failed on {csv}: {e!r}")
            plt.close("all")
            p = None
        if p:
            written.append(p)
    return written


if __name__ == "__main__":
    import sys
    rd = sys.argv[1] if len(sys.argv) > 1 else "results"
    fd = sys.argv[2] if len(sys.argv) > 2 else "figures"
    for p in make_all(rd, fd):
        print(p)
