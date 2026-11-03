"""Figures of the digital-twin (DT4DDS) end-to-end study and of the calibration study.

All numbers come from d3ra-dna/results/ (dt_summary.csv is written by
statistical_analysis/dt_analysis.py; run that first).

    python research_paper/experiments/scripts/make_dt_figures.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PROJ = Path(__file__).resolve().parents[3]
RES = PROJ / "results"
OUT = PROJ / "research_paper" / "manuscript" / "figures"
COL = 3.45
OK = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#000000", "#999999"]
plt.rcParams.update({"font.size": 7.5, "axes.titlesize": 7.5, "axes.labelsize": 7.5, "legend.fontsize": 6.3,
                     "xtick.labelsize": 6.8, "ytick.labelsize": 6.8, "font.family": "serif",
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
                     "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.4})
NAMES = {"round_robin_stop": "round-robin + verified stop", "coverage": "coverage-greedy", "pivotal": "pivotal",
         "thompson": "Thompson sampling", "ucb": "UCB ($q{=}0.9$)", "learn_then_commit": "learn-then-commit",
         "info_greedy": "info-greedy", "tail_aware": "tail-aware (LCB)"}
ORDER = list(NAMES)


def save(fig, name):
    fig.tight_layout(pad=0.3)
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)
    print("wrote", name)


def fig_dt_policies(S):
    """(a) baseline: mean reads (CI) per adaptive policy, q99 as open markers, static plans as lines;
    (b) relative cost vs round-robin on jointly successful retrievals in every scenario."""
    fig, ax = plt.subplots(1, 2, figsize=(2 * COL, 2.0), gridspec_kw={"width_ratios": [1.1, 1]})
    b = S[S.scenario == "baseline"].set_index("policy")
    y = np.arange(len(ORDER))[::-1]
    for k, p in enumerate(ORDER):
        if p not in b.index:
            continue
        r = b.loc[p]
        ax[0].errorbar(r.mean_reads / 1e3, y[k], xerr=[[(r.mean_reads - r.mean_lo) / 1e3], [(r.mean_hi - r.mean_reads) / 1e3]],
                       fmt="o", color=OK[0], ms=3.5, lw=0.9, capsize=1.5, label="mean (pay per read)" if k == 0 else None)
        ax[0].plot(r.q99_reads / 1e3, y[k], "s", mfc="none", color=OK[3], ms=3.5, label="q99 (reserved budget)" if k == 0 else None)
    st = b[b.index.str.startswith("static_")]
    if len(st):
        best = st.mean_reads.idxmin()
        ax[0].axvline(b.loc["static_a0.0"].mean_reads / 1e3, color=OK[7], ls=":", lw=0.9, label="fixed depth (uniform, hindsight)")
        ax[0].axvline(st.loc[best].mean_reads / 1e3, color=OK[2], ls="--", lw=0.9, label="best static (hindsight)")
    ax[0].set_yticks(y)
    ax[0].set_yticklabels([NAMES[p] for p in ORDER])
    ax[0].set_xscale("log")
    ax[0].set(xlabel="reads per retrieval (thousands, log)", title="(a) baseline scenario, 300 retrievals per policy")
    ax[0].legend(loc="lower right", frameon=False, fontsize=5.6)
    scen = [s for s in ["baseline", "fail1", "photolitho", "homog"] if s in set(S.scenario)]
    pols = [p for p in ORDER if p != "round_robin_stop"]
    w = 0.8 / max(len(scen), 1)
    for i, s in enumerate(scen):
        d = S[(S.scenario == s)].set_index("policy")
        vals = [d.loc[p].rel_vs_rr * 100 if p in d.index and pd.notna(d.loc[p].get("rel_vs_rr", np.nan)) else np.nan for p in pols]
        lo = [d.loc[p].rel_vs_rr_lo * 100 if p in d.index else np.nan for p in pols]
        hi = [d.loc[p].rel_vs_rr_hi * 100 if p in d.index else np.nan for p in pols]
        yy = np.arange(len(pols))[::-1] + (i - (len(scen) - 1) / 2) * w
        ax[1].barh(yy, vals, height=w, color=OK[i], label=s)
        ax[1].errorbar(vals, yy, xerr=[np.subtract(vals, lo), np.subtract(hi, vals)], fmt="none", ecolor="k", lw=0.5, capsize=1)
    ax[1].axvline(0, color="k", lw=0.6)
    ax[1].set_yticks(np.arange(len(pols))[::-1])
    ax[1].set_yticklabels([NAMES[p] for p in pols])
    ax[1].set_xlim(-50, 150)
    ax[1].set(xlabel="mean reads vs round-robin + verified stop (%)", title="(b) paired, jointly successful retrievals")
    ax[1].legend(frameon=False, loc="lower right", fontsize=5.6)
    save(fig, "fig_dt_policies")


def fig_dt_sweep(S):
    sw = S[S.scenario.str.startswith("sweep_cov")].copy()
    if sw.empty:
        return
    sw["cov"] = sw.scenario.str.replace("sweep_cov", "").astype(float)
    fig, ax = plt.subplots(1, 2, figsize=(2 * COL, 1.7))
    for k, p in enumerate(["round_robin_stop", "coverage", "thompson", "tail_aware"]):
        d = sw[sw.policy == p].sort_values("cov")
        if d.empty:
            continue
        ax[0].errorbar(d["cov"], d.mean_reads / 1e3, yerr=[(d.mean_reads - d.mean_lo) / 1e3, (d.mean_hi - d.mean_reads) / 1e3],
                       fmt="o-", ms=3, lw=0.9, capsize=1.5, color=OK[k], label=NAMES[p])
        ax[1].plot(d["cov"], d.success, "o-", ms=3, lw=0.9, color=OK[k], label=NAMES[p])
    for a in ax:
        a.set_xscale("log")
        a.set_xlabel("copies per strand in the degraded container")
    ax[0].set(ylabel="mean reads (thousands)", title="(a) cost (failures counted at the 72k-read stream end)")
    ax[1].set(ylabel="P(file byte-identical)", ylim=(-0.03, 1.03), title="(b) success")
    ax[0].legend(frameon=False, fontsize=5.6)
    save(fig, "fig_dt_sweep")


def fig_calibration():
    f = RES / "dt_calibration_select_baseline.csv"
    if not f.exists():
        return
    d = pd.read_csv(f)
    names = ["fresh", "aged 1 hl", "aged 2 hl, 10 copies", "PCR 30 cyc.", "PCR 45 cyc., 20 copies", "near-dead (3 copies)"]
    fig, ax = plt.subplots(1, 2, figsize=(2 * COL, 1.7))
    rhos = sorted(d.rho.unique())
    xs = np.arange(len(rhos))
    for j in sorted(d.container.unique()):
        c = d[d.container == j].groupby("rho").covered.mean().reindex(rhos)
        ax[0].plot(xs, c.values, "o-", ms=2.5, lw=0.8, color=OK[j], label=names[j])
    ax[0].plot(xs, d.groupby("rho").covered.mean().reindex(rhos).values, "k-", lw=1.8, label="mean")
    ax[0].axhline(0.9, color="k", ls=":", lw=0.8)
    ax[0].set_xticks(xs)
    ax[0].set_xticklabels([str(r) for r in rhos])
    ax[0].set(xlabel=r"beta-binomial ICC $\rho_U$", ylabel="90% CI coverage of $\\theta^\\ast$",
              title="(a) selection replicates 1000-1009", ylim=(-0.03, 1.03))
    ax[0].legend(frameon=False, fontsize=5.2, ncol=2, loc="lower right")
    fe = RES / "dt_calibration_eval_baseline.csv"
    if fe.exists():
        e = pd.read_csv(fe)
        c = e.groupby("container").covered.mean()
        ax[1].bar(np.arange(len(c)), c.values, color=[OK[j] for j in c.index])
        ax[1].axhline(0.9, color="k", ls=":", lw=0.8)
        ax[1].set_xticks(np.arange(len(c)))
        ax[1].set_xticklabels([n.split(" (")[0].replace(", ", "\n") for n in names], fontsize=5.4)
        rho = e.rho.iloc[0]
        ax[1].set(ylabel="coverage", ylim=(0, 1.05), title=f"(b) held-out replicates 0-29, $\\rho_U={rho}$")
    save(fig, "fig_dt_calibration")


def _rel(df, bins=10):
    edges = np.linspace(0, 1, bins + 1)
    ix = np.clip(np.digitize(df.p_complete, edges) - 1, 0, bins - 1)
    out = []
    for b in range(bins):
        m = ix == b
        if m.sum() >= 5:
            out.append((df.p_complete[m].mean(), df.completed[m].mean(), int(m.sum())))
    return np.array(out)


def fig_forecast():
    f = RES / "dt_forecast_baseline.csv"
    if not f.exists():
        return
    d = pd.read_csv(f)
    fig, ax = plt.subplots(1, 1, figsize=(COL, 2.0))
    for k, fc in enumerate(["normal", "mc"]):
        r = _rel(d[d.forecast == fc].reset_index(drop=True))
        if len(r):
            ax.plot(r[:, 0], r[:, 1], "o-", ms=3, lw=0.9, color=OK[k],
                    label={"normal": "analytic (normal)", "mc": "posterior simulation"}[fc])
    ax.plot([0, 1], [0, 1], "k:", lw=0.8)
    ax.set(xlabel=r"forecast $\hat P_{\mathrm{complete}}$", ylabel="observed completion frequency", xlim=(0, 1), ylim=(0, 1),
           title="DT baseline, checkpoints 30/50/70/85% pooled")
    ax.legend(frameon=False, loc="upper left")
    save(fig, "fig_forecast")


if __name__ == "__main__":
    S = pd.read_csv(RES / "dt_summary.csv")
    fig_dt_policies(S)
    fig_dt_sweep(S)
    fig_calibration()
    fig_forecast()
