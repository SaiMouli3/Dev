"""Figures of the digital-twin (DT4DDS) end-to-end study and of the calibration study.

All numbers come from d3ra-dna/results/ (dt_summary.csv is written by
statistical_analysis/dt_analysis.py; run that first).

    python research_paper/experiments/scripts/make_dt_figures.py
"""
from __future__ import annotations

import json
import os
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
    # bbox_inches="tight" so legends or labels placed outside the axes are never clipped
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("wrote", name)


def fig_dt_policies(S):
    """(a) baseline: reads consumed per policy, mean with CI and q99, against the static reservations;
    (b) paired relative cost against round-robin on jointly successful retrievals, per scenario."""
    from matplotlib.ticker import FixedLocator, ScalarFormatter
    fig, ax = plt.subplots(1, 2, figsize=(2 * COL, 2.5), gridspec_kw={"width_ratios": [1.12, 1]})
    b = S[S.scenario == "baseline"].set_index("policy")
    y = np.arange(len(ORDER))[::-1]
    for k, p in enumerate(ORDER):
        if p not in b.index:
            continue
        r = b.loc[p]
        ax[0].errorbar(r.mean_reads / 1e3, y[k],
                       xerr=[[(r.mean_reads - r.mean_lo) / 1e3], [(r.mean_hi - r.mean_reads) / 1e3]],
                       fmt="o", color=OK[0], ms=3.5, lw=0.9, capsize=1.5,
                       label="mean reads consumed" if k == 0 else None)
        ax[0].plot(r.q99_reads / 1e3, y[k], "s", mfc="none", color=OK[3], ms=3.5,
                   label="q99 of reads consumed" if k == 0 else None)
    st = b[b.index.str.startswith("static_")]
    if len(st):
        best = st.q99_reads.idxmin()
        ax[0].axvline(b.loc["static_a0.0"].q99_reads / 1e3, color=OK[7], ls=":", lw=0.9,
                      label="reserve: fixed depth")
        ax[0].axvline(st.loc[best].q99_reads / 1e3, color=OK[2], ls="--", lw=0.9,
                      label="reserve: best static plan")
    ax[0].set_yticks(y)
    ax[0].set_yticklabels([NAMES[p] for p in ORDER])
    ax[0].set_xscale("log")
    ax[0].xaxis.set_major_locator(FixedLocator([3, 4, 5, 7, 10, 15, 20, 30, 40]))
    ax[0].xaxis.set_minor_locator(FixedLocator([]))
    ax[0].xaxis.set_major_formatter(ScalarFormatter())
    ax[0].set_xlim(2.8, 42)
    ax[0].set_ylim(-0.6, len(ORDER) - 0.4)
    ax[0].set(xlabel="reads per retrieval (thousands, log scale)",
              title="(a) baseline: consumed vs. reserved")
    ax[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), frameon=False, fontsize=6.0, ncol=2,
                 handletextpad=0.4, columnspacing=1.0, borderaxespad=0.0)

    # (b) only scenarios in which policies differ measurably (a scenario with no successful
    # retrieval, such as photolitho, has no paired statistic and is omitted)
    scen = [s_ for s_ in ["baseline", "fail1", "homog"]
            if s_ in set(S.scenario) and S[S.scenario == s_].rel_vs_rr.notna().any()]
    pols = [p for p in ORDER if p != "round_robin_stop"]
    w = 0.8 / max(len(scen), 1)
    for i_, s_ in enumerate(scen):
        d = S[S.scenario == s_].set_index("policy")
        vals = [d.loc[p].rel_vs_rr * 100 if p in d.index and pd.notna(d.loc[p].get("rel_vs_rr", np.nan)) else np.nan
                for p in pols]
        lo = [d.loc[p].rel_vs_rr_lo * 100 if p in d.index else np.nan for p in pols]
        hi = [d.loc[p].rel_vs_rr_hi * 100 if p in d.index else np.nan for p in pols]
        yy = np.arange(len(pols))[::-1] + (i_ - (len(scen) - 1) / 2) * w
        ax[1].barh(yy, vals, height=w, color=OK[i_], label={"fail1": "1 failed", "homog": "homogeneous"}.get(s_, s_))
        ax[1].errorbar(vals, yy, xerr=[np.subtract(vals, lo), np.subtract(hi, vals)], fmt="none",
                       ecolor="k", lw=0.5, capsize=1)
    ax[1].axvline(0, color="k", lw=0.6)
    ax[1].set_yticks(np.arange(len(pols))[::-1])
    ax[1].set_yticklabels([NAMES[p] for p in pols])
    ax[1].set_xlim(-55, 120)
    ax[1].set_ylim(-0.6, len(pols) - 0.4)
    ax[1].set(xlabel="paired cost vs. round-robin (%), negative is cheaper",
              title="(b) allocation effect, jointly successful")
    ax[1].annotate("info-greedy runs off scale:\n+441% baseline, +81% failed", xy=(0.97, 0.52),
                   xycoords="axes fraction", ha="right", va="center", fontsize=5.6, style="italic",
                   bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))
    ax[1].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.30), fontsize=6.0, ncol=3,
                 handletextpad=0.4, columnspacing=1.0, borderaxespad=0.0)
    save(fig, "fig_dt_policies")


def fig_dt_sweep(S):
    """Severity sweep, single column: success is 1.00 at every level, so only cost is plotted and the
    success result is stated in the caption instead of a flat panel."""
    from matplotlib.ticker import FixedLocator, ScalarFormatter
    sw = S[S.scenario.str.startswith("sweep_cov")].copy()
    if sw.empty:
        return
    sw["cov"] = sw.scenario.str.replace("sweep_cov", "").astype(float)
    fig, ax = plt.subplots(1, 1, figsize=(COL, 2.1))
    for k, p in enumerate(["round_robin_stop", "coverage", "thompson", "tail_aware"]):
        d = sw[sw.policy == p].sort_values("cov")
        if d.empty:
            continue
        ax.errorbar(d["cov"], d.mean_reads / 1e3,
                    yerr=[(d.mean_reads - d.mean_lo) / 1e3, (d.mean_hi - d.mean_reads) / 1e3],
                    fmt="o-", ms=3, lw=0.9, capsize=1.5, color=OK[k], label=NAMES[p])
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator([1, 3, 10, 30, 200]))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.invert_xaxis()
    ax.set(xlabel="copies per strand in the degraded container (fewer = more degraded)",
           ylabel="mean reads consumed (thousands)")
    ax.legend(frameon=False, fontsize=6.0, loc="upper left", handletextpad=0.4)
    save(fig, "fig_dt_sweep")


def fig_calibration():
    f = RES / "dt_calibration_select_baseline.csv"
    if not f.exists():
        return
    d = pd.read_csv(f)
    names = ["fresh", "aged 1 hl", "aged 2 hl, 10 copies", "PCR 30 cyc.", "PCR 45 cyc., 20 copies", "near-dead (3 copies)"]
    fig, ax = plt.subplots(1, 2, figsize=(2 * COL, 2.1))
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
        # the held-out panel must show the SELECTED rho only; the file holds every rho on the grid and
        # pooling them silently reports a different quantity from the one the text reports
        rho = float(os.environ.get("D3RA_FIG_RHO", 0.1))
        e = e[np.isclose(e.rho, rho)]
        c = e.groupby("container").covered.mean()
        ax[1].bar(np.arange(len(c)), c.values, color=[OK[j] for j in c.index])
        ax[1].axhline(0.9, color="k", ls=":", lw=0.8)
        ax[1].set_xticks(np.arange(len(c)))
        ax[1].set_xticklabels([n.split(" (")[0].replace(", ", "\n") for n in names], fontsize=5.6,
                              rotation=20, ha="right")
        ax[1].set(ylabel="90% CI coverage of $\\theta^\\ast$", ylim=(0, 1.08),
                  title=f"(b) held-out replicates 0–29 at the selected $\\rho_U={rho:g}$")
        for k_, v_ in enumerate(c.values):
            ax[1].text(k_, v_ + 0.03, f"{v_:.2f}", ha="center", fontsize=5.6)
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
    fig, ax = plt.subplots(1, 1, figsize=(COL, 2.1))
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
