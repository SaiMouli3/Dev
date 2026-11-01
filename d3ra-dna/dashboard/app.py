"""D3RA benchmarking dashboard.

Run from the repository root:

    streamlit run dashboard/app.py
"""
from __future__ import annotations

import os
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import streamlit as st  # noqa: E402

from d3ra.visualization import plots  # noqa: E402

st.set_page_config(page_title="D3RA benchmark dashboard", layout="wide")

COST = plots.COST_UNITS
HEADLINE_COLS = ["design", "variant", "B_ref", "cap@0.99", "mean_spend@0.99", "cap@0.999", "mean_spend@0.999",
                 "p_exact@ref", "mean_rounds", "mean_containers_accessed", "redundancy_overhead"]


# --------------------------------------------------------------------------- helpers

@st.cache_data(show_spinner=False)
def load(results_dir: str, name: str, mtime: float) -> pd.DataFrame | None:   # mtime busts the cache
    return plots.load_csv(os.path.join(results_dir, name))


def get(results_dir: str, name: str) -> pd.DataFrame | None:
    p = os.path.join(results_dir, name)
    return load(results_dir, name, os.path.getmtime(p)) if os.path.exists(p) else None


def show_figure(fn, df, **kw):
    """Render a plots.* function into a temp PNG and display it."""
    with tempfile.TemporaryDirectory() as d:
        try:
            out = fn(df, os.path.join(d, "fig.png"), **kw)
        except Exception as e:                        # keep the dashboard alive
            st.warning(f"{fn.__name__} failed: {e!r}")
            return
        if out:
            st.image(out, width="stretch")
        else:
            st.info("Not enough columns in this CSV to draw the figure.")


def missing(name: str):
    st.info(f"`{name}` not found in the selected results directory yet. Run `experiments/run_all.py` "
            "(or wait for the running experiments) and press *Reload*.")


def fmt_table(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    keep = [c for c in cols if c in df.columns]
    out = df[keep].copy()
    for c in keep:
        lo, hi = f"{c}_ci_lo", f"{c}_ci_hi"
        if lo in df.columns and hi in df.columns:
            out[c] = [_ci(v, l, h) for v, l, h in zip(df[c], df[lo], df[hi])]
    return out


def _f(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return str(v)
    if np.isnan(v):
        return ""
    if np.isinf(v):
        return "inf"
    return f"{v:,.0f}" if abs(v) >= 100 else f"{v:.3f}"


def _ci(v, lo, hi):
    return f"{_f(v)} [{_f(lo)}, {_f(hi)}]" if _f(v) else ""


# --------------------------------------------------------------------------- sidebar

st.sidebar.title("D3RA dashboard")
default_dir = str(ROOT / "results")
results_dir = st.sidebar.text_input("Results directory", value=default_dir)
if st.sidebar.button("Reload"):
    load.clear()
if not os.path.isdir(results_dir):
    st.sidebar.error("Directory does not exist.")
else:
    files = sorted(f for f in os.listdir(results_dir) if f.endswith(".csv"))
    st.sidebar.caption(f"{len(files)} CSV file(s): " + ", ".join(files) if files else "No CSV files yet.")
st.sidebar.markdown(f"Costs are in **{COST}**.")

tabs = st.tabs(["Overview", "Baselines", "Hypotheses (H1, H5)", "Ablations", "Redundancy & Replication",
                "Placement optimisers", "Scalability", "Interactive simulator"])

# --------------------------------------------------------------------------- Overview
with tabs[0]:
    st.header("Overview")
    st.markdown("Designs A-F share K, parity, stripe count and file. `cap@t`: budget cap reaching reliability t; "
                "`mean_spend@t`: mean spend under that cap; `p_exact@ref`: P(exact recovery) at the reference "
                "budget (A's 0.99 budget x scenario scale). Brackets are 95% CIs.")
    b = get(results_dir, "baselines.csv")
    if b is None or "design" not in b.columns:
        missing("baselines.csv")
    else:
        d = b[~b["design"].astype(str).str.startswith("cmp:")]
        if "mean_spend@0.99" in d.columns:
            best = d[d["design"] != "F_oracle"].copy()
            best["_v"] = pd.to_numeric(best["mean_spend@0.99"], errors="coerce")
            win = best.loc[best.groupby("scenario")["_v"].idxmin().dropna()]
            st.subheader("Cheapest non-oracle design per scenario (mean spend @ 0.99)")
            st.dataframe(fmt_table(win, ["scenario", "design", "mean_spend@0.99", "cap@0.99"]), hide_index=True)
        for sc, g in d.groupby("scenario", sort=False):
            st.subheader(sc)
            st.dataframe(fmt_table(g, HEADLINE_COLS), hide_index=True, width="stretch")

# --------------------------------------------------------------------------- Baselines
with tabs[1]:
    st.header("Baselines A-F")
    b = get(results_dir, "baselines.csv")
    if b is None:
        missing("baselines.csv")
    else:
        scen = sorted(b["scenario"].dropna().unique()) if "scenario" in b.columns else []
        pick = st.multiselect("Scenarios", scen, default=scen)
        bb = b[b["scenario"].isin(pick)] if pick else b
        show_figure(plots.plot_baselines, bb)
        st.subheader("Pairwise comparisons")
        show_figure(plots.plot_baseline_comparisons, bb)
        cmp_ = bb[bb["design"].astype(str).str.startswith("cmp:")] if "design" in bb.columns else bb.iloc[:0]
        if not cmp_.empty:
            st.dataframe(cmp_.dropna(axis=1, how="all"), hide_index=True, width="stretch")
        with st.expander("Raw table"):
            st.dataframe(bb, width="stretch")

# --------------------------------------------------------------------------- Hypotheses
with tabs[2]:
    st.header("H1 - predictive accuracy of channel models")
    h1 = get(results_dir, "h1_model_accuracy.csv")
    if h1 is None:
        missing("h1_model_accuracy.csv")
    else:
        show_figure(plots.plot_h1, h1)
        if "experiment" in h1.columns:
            st.dataframe(h1[h1["experiment"] == "H1-summary"].dropna(axis=1, how="all"), hide_index=True)
    st.header("H5 - robustness to mis-estimated parameters")
    h5 = get(results_dir, "h5_robustness.csv")
    if h5 is None:
        missing("h5_robustness.csv")
    else:
        show_figure(plots.plot_h5, h5)
        st.dataframe(fmt_table(h5, ["scenario", "controller", "cap@0.99", "mean_spend@0.99", "cap@0.999",
                                    "mean_spend@0.999", "mean_cost_success"]), hide_index=True, width="stretch")

# --------------------------------------------------------------------------- Ablations
with tabs[3]:
    st.header("Ablations (relative to full design E)")
    ab = get(results_dir, "ablations.csv")
    if ab is None:
        missing("ablations.csv")
    else:
        metric = st.radio("Metric", [c for c in ("mean_spend@0.99", "cap@0.99", "mean_spend@0.999", "cap@0.999")
                                     if c in ab.columns] or ["mean_spend@0.99"], horizontal=True)
        show_figure(plots.plot_ablations, ab, metric=metric)
        st.dataframe(fmt_table(ab, ["scenario", "ablation", "variant", "cap@0.99", "mean_spend@0.99", "cap@0.999",
                                    "mean_spend@0.999"]), hide_index=True, width="stretch")

# --------------------------------------------------------------------------- Redundancy
with tabs[4]:
    st.header("Redundancy allocation")
    rd = get(results_dir, "redundancy.csv")
    if rd is None:
        missing("redundancy.csv")
    else:
        show_figure(plots.plot_redundancy, rd)
        st.dataframe(rd.dropna(axis=1, how="all"), hide_index=True, width="stretch")
    st.header("Replication vs RS parity under shared sequence bias")
    rp = get(results_dir, "replication_vs_parity.csv")
    if rp is None:
        missing("replication_vs_parity.csv")
    else:
        show_figure(plots.plot_replication_vs_parity, rp)
        st.dataframe(fmt_table(rp, ["shared_variance_fraction", "scheme", "stored_molecules", "cap@0.99",
                                    "mean_cost_success"]), hide_index=True, width="stretch")

# --------------------------------------------------------------------------- Placement
with tabs[5]:
    st.header("Placement optimisers vs ILP")
    po = get(results_dir, "placement_opt.csv")
    if po is None:
        missing("placement_opt.csv")
    else:
        show_figure(plots.plot_placement_gap, po)
        st.dataframe(po.dropna(axis=1, how="all"), hide_index=True, width="stretch")

# --------------------------------------------------------------------------- Scalability
with tabs[6]:
    st.header("Scalability")
    sc_df = get(results_dir, "scalability.csv")
    if sc_df is None:
        missing("scalability.csv")
    else:
        show_figure(plots.plot_scalability, sc_df)
        st.dataframe(sc_df.dropna(axis=1, how="all"), hide_index=True, width="stretch")

# --------------------------------------------------------------------------- Simulator
PLACEMENTS = {"round robin": "round_robin", "reliability-weighted": "reliability",
              "cost-tiered (slack 1)": "tiered1", "domain-aware": "domain"}
CONTROLLERS = ["pivotal", "thompson", "coverage", "uniform", "oracle"]


@st.cache_resource(show_spinner=False)
def scenario(name: str, seed: int):
    from d3ra.evaluation import scenarios as S
    return S.make(name, seed, calibrated=False)


def build_layout(sc, placement: str):
    from d3ra.evaluation.baselines import reliability_layout, uniform_layout
    from d3ra.optimization.joint_optimizer import _layout
    if placement == "round_robin":
        return uniform_layout(sc)
    if placement == "reliability":
        return reliability_layout(sc)
    return _layout(placement, sc.n_stripes, sc.k, sc.parity, sc.belief, sc.domain)


@st.cache_data(show_spinner=False)
def simulate(name: str, placement: str, controller: str, trials: int, seed: int, cap_mult: float, prior_sd: float):
    from d3ra.channel.surrogate import sample_T
    from d3ra.evaluation.metrics import adaptive_summary, static_summary
    from d3ra.retrieval.adaptive_policy import AdaptiveConfig, AdaptiveController
    from d3ra.retrieval.fixed_policy import StaticPolicy, budget_for_target, uniform_weights

    sc = scenario(name, seed)
    L = build_layout(sc, placement)
    t0 = time.time()
    # static uniform baseline on the same layout (more trials: it is cheap)
    static = StaticPolicy(uniform_weights(L)).budget_samples(L, sc.truth, max(4 * trials, 2000),
                                                             np.random.default_rng(seed))
    t_static = time.time() - t0
    fin = static[np.isfinite(static)]
    # adaptive cost cap: cap_mult x static 0.99 budget (or its finite 0.97 quantile when 0.99 is unreachable,
    # e.g. under container failures) -- trials that cannot complete stop at the cap
    ref = budget_for_target(static, 0.99)
    if not np.isfinite(ref):
        ref = float(np.quantile(fin, 0.97)) if len(fin) else 1e5
    cap = cap_mult * ref
    rng = np.random.default_rng(seed + 1)
    T, failed = sample_T(L, sc.truth, trials, rng, return_state=True)
    ctl = AdaptiveController(L, sc.belief, AdaptiveConfig(controller, cap=cap, prior_sd=prior_sd), truth=sc.truth)
    t1 = time.time()
    r = ctl.run(trials, np.random.default_rng(seed + 2), T=T, failed=failed)
    t_ad = time.time() - t1
    adaptive = np.where(r.success, r.cost, np.inf)
    s_ad = adaptive_summary(adaptive, r.success, (0.99, 0.999), ref, r.rounds, r.containers_used, r.extra_after_probe)
    s_st = static_summary(static, (0.99, 0.999), ref)
    info = {"description": sc.description, "fragments": int(L.n), "containers": int(L.m),
            "stored_molecules": float(L.placement.sum()), "redundancy_overhead": float(L.placement.sum() / L.k.sum() - 1),
            "cap": cap, "ref_budget": ref, "t_static": t_static, "t_adaptive": t_ad,
            "reads_per_container": r.reads_per_container.mean(0), "stored_per_container": L.stored_per_container(),
            "cost_per_read": sc.truth.cost, "yield": sc.truth.yield_}
    return adaptive, static, s_ad, s_st, info


def _metric_rows(s_ad: dict, s_st: dict, ctl: str) -> pd.DataFrame:
    keys = ["p_success_uncapped", "cap@0.99", "mean_spend@0.99", "cap@0.999", "mean_spend@0.999",
            "mean_cost_success", "p_exact@ref", "mean_rounds", "mean_containers_accessed"]
    rows = []
    for lab, s in [(f"adaptive ({ctl})", s_ad), ("static uniform", s_st)]:
        row = {"policy": lab}
        for k in keys:
            if k in s:
                ci = s.get(f"{k}_ci")
                row[k] = _ci(s[k], *ci) if isinstance(ci, tuple) else _f(s[k])
        rows.append(row)
    return pd.DataFrame(rows)


with tabs[7]:
    st.header("Interactive simulator")
    st.caption("Uncalibrated (analytic) consensus curves for speed. The adaptive controller and a static uniform "
               "allocation are evaluated on the same layout and true channel; the controller only sees the "
               "scenario's *belief* parameters (except `oracle`).")
    from d3ra.evaluation.scenarios import SCENARIOS

    c1, c2, c3 = st.columns(3)
    name = c1.selectbox("Scenario", SCENARIOS, index=SCENARIOS.index("s2_heterogeneous"))
    placement_lab = c2.selectbox("Placement", list(PLACEMENTS), index=2)
    controller = c3.selectbox("Controller", CONTROLLERS)
    c4, c5, c6, c7 = st.columns(4)
    trials = c4.slider("Trials", 20, 500, 150, step=10)
    seed = c5.number_input("Seed", 0, 10_000, 0)
    cap_mult = c6.slider("Cost cap (x static reference budget)", 1.0, 8.0, 4.0, step=0.5)
    prior_sd = c7.slider("Prior sd (estimation uncertainty)", 0.02, 1.5, 0.5, step=0.02)
    st.caption("Runtime grows with trials and with the share of trials that hit the cap (e.g. container failures).")
    if st.button("Run simulation", type="primary"):
        st.session_state["sim_args"] = (name, PLACEMENTS[placement_lab], controller, int(trials), int(seed),
                                        float(cap_mult), float(prior_sd))
    args = st.session_state.get("sim_args")
    if args:
        with st.spinner("Simulating..."):
            try:
                adaptive, static, s_ad, s_st, info = simulate(*args)
            except Exception as e:
                st.error(f"Simulation failed: {e!r}")
                adaptive = None
        if adaptive is not None:
            st.markdown(f"**{args[0]}** - {info['description']}  \n"
                        f"{info['fragments']} fragments on {info['containers']} containers, redundancy overhead "
                        f"{info['redundancy_overhead']:.1%}; reference budget {info['ref_budget']:,.0f}, cap "
                        f"{info['cap']:,.0f} {COST}. Runtime: static {info['t_static']:.1f}s, adaptive "
                        f"{info['t_adaptive']:.1f}s.")
            m1, m2, m3, m4 = st.columns(4)
            fa, fs = adaptive[np.isfinite(adaptive)], static[np.isfinite(static)]
            m1.metric("Adaptive mean cost (successes)", _f(fa.mean()) if len(fa) else "-")
            m2.metric("Static mean cost (successes)", _f(fs.mean()) if len(fs) else "-",
                      delta=(f"{(fa.mean() / fs.mean() - 1):+.1%} adaptive vs static" if len(fa) and len(fs) else None),
                      delta_color="off")
            m3.metric("Adaptive success rate (within cap)", f"{np.isfinite(adaptive).mean():.1%}")
            m4.metric("Static success rate (unbounded)", f"{np.isfinite(static).mean():.1%}")
            st.dataframe(_metric_rows(s_ad, s_st, args[2]), hide_index=True, width="stretch")

            with plt.rc_context(plots.STYLE):
                fig, axes = plt.subplots(1, 3, figsize=(15, 3.8))
                hi = np.quantile(np.concatenate([fa, fs]), 0.995) if len(fa) + len(fs) else 1
                bins = np.linspace(0, hi, 50)
                axes[0].hist(fs, bins=bins, color=plots.PALETTE[1], alpha=0.65, density=True, label="static uniform")
                axes[0].hist(fa, bins=bins, color=plots.PALETTE[0], alpha=0.65, density=True, label=f"adaptive ({args[2]})")
                axes[0].set_xlabel(f"Cost to recover, {COST}")
                axes[0].set_ylabel("Density")
                axes[0].set_title("Cost distribution (successful trials)")
                axes[0].legend(frameon=False)
                for x, lab, col in [(static, "static uniform", plots.PALETTE[1]),
                                    (adaptive, f"adaptive ({args[2]})", plots.PALETTE[0])]:
                    s = np.sort(x)
                    axes[1].step(s, np.arange(1, len(s) + 1) / len(s), where="post", color=col, label=lab)
                for t in (0.99, 0.999):
                    axes[1].axhline(t, color="#999", ls=":", lw=0.8)
                axes[1].axvline(info["cap"], color=plots.PALETTE[3], ls="--", lw=0.9, label="adaptive cap")
                axes[1].set_xlim(0, max(hi, info["cap"]) * 1.05)
                axes[1].set_ylim(0, 1.01)
                axes[1].set_xlabel(f"Budget, {COST}")
                axes[1].set_ylabel("P(recovered within budget)")
                axes[1].set_title("ECDF of cost (failures never reach 1)")
                axes[1].legend(frameon=False, loc="lower right")
                x = np.arange(info["containers"])
                axes[2].bar(x - 0.2, info["reads_per_container"], 0.4, color=plots.PALETTE[0], label="mean reads (adaptive)")
                ax2 = axes[2].twinx()
                ax2.bar(x + 0.2, info["stored_per_container"], 0.4, color=plots.PALETTE[2], label="stored fragments")
                ax2.grid(False)
                ax2.set_ylabel("Stored fragment copies")
                axes[2].set_xticks(x)
                axes[2].set_xticklabels([f"C{j}\n{c:.1f}" for j, c in enumerate(info["cost_per_read"])], fontsize=7)
                axes[2].set_xlabel("Container (cost/read)")
                axes[2].set_ylabel("Reads per trial")
                axes[2].set_title("Where the controller reads")
                h1_, l1_ = axes[2].get_legend_handles_labels()
                h2_, l2_ = ax2.get_legend_handles_labels()
                axes[2].legend(h1_ + h2_, l1_ + l2_, frameon=False, fontsize=7)
                st.pyplot(fig)
                plt.close(fig)
    else:
        st.info("Choose settings and press **Run simulation**.")
