"""Generate all data figures of the IEEE manuscript from the committed result files.

Every figure here is computed from files in d3ra-dna/results/ (no hand-entered
numbers). Conceptual diagrams are drawn in TikZ inside main.tex.

    python research_paper/experiments/scripts/make_figures.py
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
OUT.mkdir(parents=True, exist_ok=True)

COL = 3.45      # IEEE single-column width (in)
OK = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#000000"]
plt.rcParams.update({"font.size": 7.5, "axes.titlesize": 7.5, "axes.labelsize": 7.5, "legend.fontsize": 6.3,
                     "xtick.labelsize": 6.8, "ytick.labelsize": 6.8, "font.family": "serif",
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42,
                     "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.4})


def save(fig, name):
    fig.tight_layout(pad=0.3)
    fig.savefig(OUT / f"{name}.pdf")
    plt.close(fig)
    print("wrote", name)


def fig_consensus():
    ch = pd.read_csv(RES / "realdata_channel.csv")
    grid = [1, 2, 3, 4, 5, 6, 8, 10, 13, 16]
    val = pd.read_csv(RES / "q3_model_validation.csv")
    fig, ax = plt.subplots(1, 2, figsize=(2 * COL, 1.75))
    for i, r in ch.iterrows():
        ax[0].plot(grid, [r[f"phi_exact_n{n}"] for n in grid], "o-", ms=2.5, lw=0.9, color=OK[i % 7], label=r["condition"].replace("_", " "))
    ax[0].set(xlabel="reads per strand $n$", ylabel=r"$\hat\varphi_j(n)$ (exact consensus)", xlim=(0.5, 16), ylim=(0, 1.02),
              title="(a) Illumina, 6 handling conditions (PRJEB65931)")
    ax[0].legend(ncol=2, loc="lower right", frameon=False)
    ax[1].plot(val.n_reads, val.real_FC2, "o-", color=OK[3], ms=3, label="real reads (flowcell 2)")
    ax[1].plot(val.n_reads, val.iid, "s--", color=OK[0], ms=3, label="i.i.d. errors (rates fitted on FC1)")
    ax[1].plot(val.n_reads, val.context_only, "^--", color=OK[5], ms=3, label="context-only model")
    ax[1].plot(val.n_reads, val.hierarchical, "D-", color=OK[2], ms=3, label="hierarchical beta-binomial")
    ax[1].set(xlabel="reads per strand $n$", ylabel="P(exact consensus)", ylim=(0, 1.02),
              title="(b) Nanopore, held-out strands (DNAformer)")
    ax[1].legend(loc="lower right", frameon=False)
    save(fig, "fig_consensus")


# v2 = policy tables re-run on the current real-data caches (the v1 files are superseded, kept for transparency)
TW_POL = RES / "realdata_policies_v2.csv"
GS_POL = [RES / "realdata_genscript_policies_rr_v2.csv", RES / "realdata_genscript_policies_sb_v2.csv"]


def fig_cost_real():
    tw = pd.read_csv(TW_POL)
    gs = pd.concat([pd.read_csv(f) for f in GS_POL])
    pols = [("static_uniform", "fixed depth"), ("static_oracle_weights", "oracle-tuned static"),
            ("adaptive_uniform", "adaptive stop, round-robin"), ("adaptive_pivotal", "adaptive (proposed)")]
    fig, ax = plt.subplots(1, 2, figsize=(2 * COL, 2.05))
    for a, (d, title, places) in zip(ax, [(tw, "(a) Twist pool", ["round_robin", "reliability", "tiered2", "single_best"]),
                                          (gs, "(b) Genscript pool (independent)", ["round_robin", "single_best"])]):
        x = np.arange(len(places))
        w = 0.2
        for k, (p, lab) in enumerate(pols):
            sub = d[d.policy == p].set_index("placement").reindex(places)
            y = sub["mean_spend@0.99"].values / 1e3
            lo = y - sub["mean_spend@0.99_ci_lo"].values / 1e3
            hi = sub["mean_spend@0.99_ci_hi"].values / 1e3 - y
            a.bar(x + (k - 1.5) * w, y, w, yerr=[np.clip(lo, 0, None), np.clip(hi, 0, None)], color=OK[k], label=lab,
                  capsize=1.5, error_kw={"lw": 0.6})
        a.set_xticks(x, [p.replace("_", " ").replace("tiered2", "tiered") for p in places])
        a.set(ylabel="reads to reach $P_{exact}\\geq0.99$ ($\\times10^3$)", title=title)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.0))
    fig.tight_layout(pad=0.3, rect=(0, 0, 1, 0.9))
    fig.savefig(OUT / "fig_cost_real.pdf")
    plt.close(fig)
    print("wrote fig_cost_real")


def fig_failure():
    fi = pd.read_csv(RES / "realdata_failure_injection_v2.csv")
    fig, ax = plt.subplots(figsize=(COL, 1.9))
    names = list(dict.fromkeys(fi.lost_container))
    x = np.arange(len(names))
    for k, (pl, g) in enumerate(fi.groupby("placement", sort=False)):
        g = g.set_index("lost_container").reindex(names)
        cost = g.adaptive_mean_cost_success.values / 1e3
        ok = g.adaptive_p_success_within_cap.values
        ax.plot(x, np.where(np.isfinite(cost), cost, np.nan), "o-", ms=3, lw=0.9, color=OK[k], label=pl.replace("_", " "))
        for xi, c, p in zip(x, cost, ok):
            if not np.isfinite(c) or p < 0.999:
                ax.annotate("lost" if not np.isfinite(c) else f"{p:.2f}", (xi, (c if np.isfinite(c) else 40)), fontsize=5.5,
                            color=OK[k], xytext=(2, 3), textcoords="offset points")
    ax.set_xticks(x, [n.replace("_", " ") for n in names], rotation=30, ha="right")
    ax.set(ylabel="mean reads to recover ($\\times10^3$)", xlabel="container removed")
    ax.legend(frameon=False, fontsize=5.8)
    save(fig, "fig_failure")


def fig_confidence():
    f = RES / "paper_ec_confidence_bins.csv"
    if not f.exists():
        return
    b = pd.read_csv(f)
    allb = b[(b.spend_fraction.astype(str) == "all") & (b.n >= 5)]  # bins with <5 trials omitted
    fig, ax = plt.subplots(figsize=(COL, 2.0))
    ax.plot([0, 1], [0, 1], "k:", lw=0.8, label="perfect calibration")
    for k, (name, g) in enumerate(allb.groupby("belief")):
        ax.plot(g.mean_p, g.empirical, "o-", color=OK[k], ms=3, label=name)
        for _, r in g.iterrows():
            ax.annotate(int(r.n), (r.mean_p, r.empirical), fontsize=5, xytext=(2, -6), textcoords="offset points", color=OK[k])
    ax.set(xlabel="predicted P(complete within budget)", ylabel="observed completion frequency", xlim=(0, 1), ylim=(0, 1))
    ax.legend(frameon=False, loc="lower right")
    save(fig, "fig_confidence")


def fig_limited():
    f = RES / "paper_e4_limited_calibration.csv"
    if not f.exists():
        return
    d = pd.read_csv(f)
    order = list(dict.fromkeys(d.belief))
    labels = []
    for s in order:
        if s.startswith("generic"):
            labels.append("0")
        elif s.startswith("counts"):
            labels.append(s.split()[2].replace(",", "")[:-3] + "k")
        else:
            labels.append("300k+$\\hat\\varphi$")
    tw = pd.read_csv(TW_POL)
    ref_or = tw[(tw.placement == "round_robin") & (tw.policy == "static_oracle_weights")]["mean_spend@0.99"].iloc[0] / 1e3
    ref_fx = tw[(tw.placement == "round_robin") & (tw.policy == "static_uniform")]["mean_spend@0.99"].iloc[0] / 1e3
    fig, ax = plt.subplots(figsize=(COL, 1.9))
    x = np.arange(len(order))
    for k, psd in enumerate(sorted(d.prior_sd.unique())):
        g = d[d.prior_sd == psd].set_index("belief").reindex(order)
        y = g["mean_spend@0.99"].values / 1e3
        ax.errorbar(x + (k - 0.5) * 0.12, y, yerr=[y - g["mean_spend@0.99_ci_lo"].values / 1e3, g["mean_spend@0.99_ci_hi"].values / 1e3 - y],
                    fmt="o-", ms=3, lw=0.9, capsize=1.5, color=OK[k], label=f"adaptive, prior sd {psd}")
    ax.axhline(ref_or, color=OK[3], ls="--", lw=0.8, label="oracle-tuned static")
    ax.axhline(ref_fx, color="grey", ls=":", lw=0.8, label="fixed depth")
    ax.set_xticks(x, labels)
    ax.set(xlabel="calibration reads per container used to build the prior", ylabel="reads @0.99 ($\\times10^3$)")
    ax.set_ylim(53, 70)
    ax.legend(frameon=False, fontsize=5.8, loc="center right", ncol=2)
    save(fig, "fig_limited_calibration")


def fig_ablation():
    d = pd.read_csv(RES / "ablations.csv")
    d = d[d.scenario == "s4_correlated_failures"]
    full = d[d.ablation == "full (E)"]["mean_cost_success"].iloc[0]
    fig, ax = plt.subplots(figsize=(COL, 2.1))
    lab = {"full (E)": "full design", "A: -reliability-aware placement": "uniform placement",
           "B: -adaptive sequencing (static optimised weights)": "static reads (no adaptivity)",
           "D: -correlated failure modelling": "ignore correlated failures",
           "E: fixed stopping (controller's mean allocation, committed)": "fixed stopping",
           "F: perfect channel information (oracle)": "oracle channel knowledge",
           "G: independent instead of joint optimisation": "independent optimisation",
           "adaptive stopping only (uniform round-robin allocation)": "adaptive stop, round-robin alloc.",
           "coverage-greedy instead of pivotal": "coverage-greedy score"}
    d = d[d.ablation.isin(lab)]
    y = np.arange(len(d))
    rel = d.mean_cost_success.values / full - 1
    succ = d.p_success_uncapped.values
    colors = [OK[2] if s >= 0.99 else OK[3] for s in succ]
    ax.barh(y, rel * 100, color=colors)
    for yi, r, s in zip(y, rel, succ):
        ax.text(max(r * 100, 0) + 2, yi, f"{r*100:+.1f}%, $P$={s:.3f}", va="center", ha="left", fontsize=5.5)
    ax.set_xlim(-8, 105)
    ax.set_yticks(y, [lab[a] for a in d.ablation])
    ax.axvline(0, color="k", lw=0.6)
    ax.set(xlabel="change in mean reads of successful retrievals vs full design (%)")
    ax.invert_yaxis()
    save(fig, "fig_ablation")


if __name__ == "__main__":
    for f in (fig_consensus, fig_cost_real, fig_failure, fig_confidence, fig_limited, fig_ablation):
        f()
