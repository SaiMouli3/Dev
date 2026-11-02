"""Experiments added for the IEEE manuscript (research questions RQ1-RQ3 and
confidence calibration). All use real read streams (empirical recovery-time
matrices) unless stated.

  E4  limited-observation self-calibration (Twist pool)        -> results/paper_e4_limited_calibration.csv
  EC  posterior-predictive confidence calibration (Twist pool)  -> results/paper_ec_confidence_calibration.csv / _bins.csv
  ED  change detection: simulated (s7) and real (Twist)         -> results/paper_ed_change_detection.csv
  E10 unfamiliar pool: Genscript with in-pool / Twist / generic beliefs -> results/paper_e10_cross_pool.csv

Usage:
  D3RA_POOL=twist     python experiments/paper_experiments.py e4 ec ed
  D3RA_POOL=genscript python experiments/paper_experiments.py e10
"""
from __future__ import annotations

import json
import os
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

import real_data as RD  # noqa: E402  (pool selected by D3RA_POOL at import)
from d3ra.channel.failures import FailureModel  # noqa: E402
from d3ra.channel.surrogate import ChannelParams, Layout  # noqa: E402
from d3ra.evaluation.metrics import adaptive_summary, static_summary  # noqa: E402
from d3ra.realdata import empirical as EM  # noqa: E402
from d3ra.retrieval.adaptive_policy import THETA, AdaptiveConfig, AdaptiveController  # noqa: E402
from d3ra.storage.placement import round_robin  # noqa: E402

RES = ROOT / "results"
SEED = 0


def generic_prior(m: int) -> ChannelParams:
    """Belief for a completely unfamiliar pool: literature-level Illumina defaults
    (yield 0.9, NB shape 10, dropout 0.5%, phi from the simulator at the
    published Illumina error profile), identical for every container."""
    from d3ra.evaluation.scenarios import phi_for_scale
    phi = np.tile(phi_for_scale(1.0), (m, 1))
    return ChannelParams(np.full(m, 0.9), np.full(m, 10.0), np.full(m, 0.005), np.ones(m), phi,
                         FailureModel(np.full(m, 1e-3)), sigma_g=0.3, d_syn=0.0)


def counts_belief(C, n_cal: int, rng) -> ChannelParams:
    """Belief estimated from only n_cal calibration reads per container (binomial
    thinning of the 300k-read calibration counts); phi stays generic because
    n_cal << strands makes per-strand consensus statistics unavailable."""
    runs = [C[r]["run"] for r in RD.RUNS]
    m = len(runs)
    g = generic_prior(m)
    shapes, drops, yields = [], [], []
    for rd in runs:
        p = min(1.0, n_cal / max(rd.calib_counts.sum(), 1))
        c = rng.binomial(rd.calib_counts, p).astype(float)
        mu, var = c.mean(), c.var()
        # method of moments; with very few reads the estimate is noisy, so shrink to the prior
        k_hat = mu ** 2 / max(var - mu, 1e-9) if var > mu else 200.0
        w = min(1.0, mu / 5.0)                      # trust the estimate once ~5 reads/strand
        shapes.append(float(np.exp(w * np.log(np.clip(k_hat, 0.5, 200)) + (1 - w) * np.log(10.0))))
        exp_zero = np.exp(-mu) if mu > 0 else 1.0
        drops.append(float(np.clip(w * ((c == 0).mean() - exp_zero) + (1 - w) * 0.005, 1e-4, 0.9)))
        yields.append(float(rd.calib_yield))
    return replace(g, yield_=np.array(yields), shape=np.array(shapes), dropout=np.array(drops))


def _setup_layout(n_stripes=250, k=32, parity=16):
    belief_full, runs_v, phi_v, stripe_of, kk = RD._real_setup(n_stripes, k, parity, SEED)
    L = Layout(stripe_of, kk, round_robin(stripe_of, belief_full.m))
    return belief_full, runs_v, phi_v, L


def _run(L, belief, T, trials, controller="pivotal", cap=np.inf, prior_sd=0.5, checkpoints=None, batch=200):
    ctl = AdaptiveController(L, belief, AdaptiveConfig(controller, cap=cap, batch=batch, prior_sd=prior_sd))
    return ctl.run(trials, np.random.default_rng(SEED + 11), T=T, failed=np.zeros((trials, L.m), bool),
                   checkpoints=checkpoints)


def _flat(d):
    out = {}
    for k, v in d.items():
        if isinstance(v, tuple):
            out[k + "_lo"], out[k + "_hi"] = v
        else:
            out[k] = v
    return out


# ----------------------------------------------------------------------------- E4

def e4(trials=100):
    C = RD.load_cached()
    belief_full, runs_v, phi_v, L = _setup_layout()
    T = EM.empirical_T(L, runs_v, phi_v, trials, np.random.default_rng(SEED + 7))
    rng = np.random.default_rng(SEED + 3)
    levels = [("generic prior (0 reads)", generic_prior(belief_full.m))]
    for n in [1_000, 3_000, 10_000, 30_000, 300_000]:
        levels.append((f"counts from {n:,} reads/container, generic phi", counts_belief(C, n, rng)))
    levels.append(("full calibration (300k reads + empirical phi)", belief_full))
    rows = []
    for name, b in levels:
        for psd in (0.5, 1.0):
            r = _run(L, b, T, trials, prior_sd=psd)
            c = np.where(r.success, r.cost, np.inf)
            s = adaptive_summary(c, r.success, (0.95, 0.99))
            rows.append({"belief": name, "prior_sd": psd, "trials": trials, **_flat(s)})
            print(rows[-1]["belief"], psd, round(s["mean_spend@0.99"]), s["p_success_uncapped"], flush=True)
    pd.DataFrame(rows).to_csv(RES / "paper_e4_limited_calibration.csv", index=False)
    (RES / "paper_twist_belief_mean.json").write_text(json.dumps({
        "yield": float(belief_full.yield_.mean()), "shape": float(np.exp(np.log(belief_full.shape).mean())),
        "dropout": float(belief_full.dropout.mean()), "phi": belief_full.phi.mean(0).tolist(),
        "sigma_g": belief_full.sigma_g}))


# ----------------------------------------------------------------------------- EC

def ec(trials=200):
    belief_full, runs_v, phi_v, L = _setup_layout()
    T = EM.empirical_T(L, runs_v, phi_v, trials, np.random.default_rng(SEED + 7))
    # budget set near the median adaptive cost so that outcomes are mixed (informative calibration)
    ref = _run(L, belief_full, T, trials)
    cap = float(np.quantile(ref.cost[ref.success], 0.5))
    cps = [f * cap for f in (0.3, 0.5, 0.7, 0.85)]
    rows, bins_rows = [], []
    for name, b in [("full calibration", belief_full), ("generic prior", generic_prior(belief_full.m))]:
        r = _run(L, b, T, trials, cap=cap, checkpoints=cps)
        outcome = r.success
        for tr, lvl, p in r.meta["checkpoints"]:
            rows.append({"belief": name, "trial": tr, "spend_fraction": round(lvl / cap, 2), "p_complete": p,
                         "completed": bool(outcome[tr])})
    df = pd.DataFrame(rows)
    df.to_csv(RES / "paper_ec_confidence_calibration.csv", index=False)
    edges = np.linspace(0, 1, 11)
    for (name, frac), g in df.groupby(["belief", "spend_fraction"]):
        p, y = g.p_complete.values, g.completed.values.astype(float)
        idx = np.clip(np.digitize(p, edges) - 1, 0, 9)
        ece = sum(abs(p[idx == b].mean() - y[idx == b].mean()) * (idx == b).sum() for b in range(10) if (idx == b).any()) / len(p)
        bins_rows.append({"belief": name, "spend_fraction": frac, "n": len(p), "mean_p": p.mean(), "empirical": y.mean(),
                          "brier": float(np.mean((p - y) ** 2)), "ece": float(ece),
                          "brier_climatology": float(np.mean((y.mean() - y) ** 2))})
    for name, g in df.groupby("belief"):
        p, y = g.p_complete.values, g.completed.values.astype(float)
        idx = np.clip(np.digitize(p, edges) - 1, 0, 9)
        for b in range(10):
            if (idx == b).any():
                bins_rows.append({"belief": name, "spend_fraction": "all", "bin": f"{edges[b]:.1f}-{edges[b+1]:.1f}",
                                  "n": int((idx == b).sum()), "mean_p": float(p[idx == b].mean()),
                                  "empirical": float(y[idx == b].mean())})
    pd.DataFrame(bins_rows).to_csv(RES / "paper_ec_confidence_bins.csv", index=False)
    print(pd.DataFrame(bins_rows).round(3).to_string(), flush=True)
    (RES / "paper_ec_meta.json").write_text(json.dumps({"cap": cap, "checkpoints": cps, "trials": trials}))


# ----------------------------------------------------------------------------- ED

def _rel_eff_and_rstar(L, T_truth, T_belief, j, frac=0.4):
    """True relative efficiency of container j = (median recovery time under the
    belief model) / (median under the truth); R* = reads at which the BELIEF model
    predicts `frac` of the container's pairs to be recovered."""
    sel = L.pair_cont == j
    tb, tt = T_belief[:, sel], T_truth[:, sel]
    mb = np.median(tb[np.isfinite(tb)])
    mt = np.median(tt[np.isfinite(tt)])
    return float(mb / mt), float(np.quantile(tb[np.isfinite(tb)], frac))


def _posterior_at(ctl, L, T, j, R_j):
    R = np.zeros((T.shape[0], L.m))
    R[:, j] = R_j
    U = L.pair_to_cont_sum((T <= R[:, L.pair_cont]).astype(np.float32))
    return ctl.posterior(R, U)[:, j, :]            # (trials, G)


def _detect_rows(setting, names, ctl, L, T_truth, T_belief, thr=0.7):
    rows = []
    cdf_idx = np.cumsum
    for j in range(L.m):
        eff, rstar = _rel_eff_and_rstar(L, T_truth, T_belief, j)
        for mult in (1, 2):
            P = _posterior_at(ctl, L, T_truth, j, mult * rstar)
            mean = P @ THETA
            cdf = np.cumsum(P, 1)
            lo = THETA[(cdf < 0.05).sum(1).clip(max=len(THETA) - 1)]
            hi = THETA[(cdf < 0.95).sum(1).clip(max=len(THETA) - 1)]
            flag = P[:, THETA < thr].sum(1) > 0.9
            rows.append({"setting": setting, "container": names[j], "reads_per_stored_strand": mult * rstar / (L.pair_cont == j).sum(),
                         "operating_point": f"{mult}x R*", "true_relative_efficiency": eff,
                         "truly_degraded": eff < thr, "posterior_mean_median": float(np.median(mean)),
                         "ci90_coverage": float(((lo <= eff) & (eff <= hi)).mean()),
                         "flag_rate": float(flag.mean())})
            print(rows[-1], flush=True)
    return rows


def ed(trials=200, od=0.0, out="paper_ed_change_detection.csv", sim=True):
    from d3ra.channel.surrogate import sample_T
    from d3ra.evaluation import scenarios as S
    from d3ra.evaluation.baselines import uniform_layout
    rows = []
    # (i) simulation with known degradation (s7: containers 0 and 2 degrade after placement)
    sc = S.make("s7_channel_change") if sim else None
    if sim:
        Ls = uniform_layout(sc)
        T_truth = sample_T(Ls, replace(sc.truth, failure=FailureModel(np.zeros(Ls.m))), trials, np.random.default_rng(SEED + 5))
        T_bel = sample_T(Ls, replace(sc.belief, failure=FailureModel(np.zeros(Ls.m))), trials, np.random.default_rng(SEED + 6))
        ctl = AdaptiveController(Ls, sc.belief, AdaptiveConfig("pivotal", prior_sd=1.0, overdispersion=od))
        rows += _detect_rows("simulated s7 (containers 0,2 degraded)", [f"C{j}" for j in range(Ls.m)], ctl, Ls, T_truth, T_bel)
    # (ii) real Twist streams; belief = every container behaves like the fresh pool
    belief_full, runs_v, phi_v, L = _setup_layout()
    T = EM.empirical_T(L, runs_v, phi_v, trials, np.random.default_rng(SEED + 7))
    f = 0
    b_fresh = replace(belief_full, yield_=np.full(L.m, belief_full.yield_[f]), shape=np.full(L.m, belief_full.shape[f]),
                      dropout=np.full(L.m, belief_full.dropout[f]), phi=np.tile(belief_full.phi[f], (L.m, 1)))
    # belief-model times: the fresh container's real stream used for every container (same layout share)
    T_bel = np.empty_like(T)
    for j in range(L.m):
        T_bel[:, L.pair_cont == j] = T[:, L.pair_cont == f][:, : (L.pair_cont == j).sum()]
    ctl = AdaptiveController(L, b_fresh, AdaptiveConfig("pivotal", prior_sd=1.0, overdispersion=od))
    rows += _detect_rows(f"real {RD.POOL} (belief = fresh pool)", list(RD.RUNS.values()), ctl, L, T, T_bel)
    for r in rows:
        r["overdispersion"] = od
    pd.DataFrame(rows).to_csv(RES / out, index=False)
    return rows


def ed_od(trials=200):
    """Exploratory: same as ED with a beta-binomial (overdispersed) likelihood, ICC in {0 (reference), 0.01, 0.05}; run with PYTHONHASHSEED=0."""
    rows = []
    for od in (0.0, 0.01, 0.05):
        rows += ed(trials, od=od, out=f"_tmp_ed_od{od}.csv")
        (RES / f"_tmp_ed_od{od}.csv").unlink()
    pd.DataFrame(rows).to_csv(RES / "paper_ed_overdispersion.csv", index=False)


def ed_holdout(trials=200):
    """Held-out check of the ICC chosen on Twist (0.01) on the independent Genscript pool (real part only)."""
    assert RD.POOL == "genscript", "run with D3RA_POOL=genscript"
    rows = []
    for od in (0.0, 0.01):
        rows += ed(trials, od=od, out=f"_tmp_edgs_od{od}.csv", sim=False)
        (RES / f"_tmp_edgs_od{od}.csv").unlink()
    pd.DataFrame(rows).to_csv(RES / "paper_ed_holdout_genscript.csv", index=False)


# ----------------------------------------------------------------------------- E10

def e10(trials=60):
    assert RD.POOL == "genscript", "run with D3RA_POOL=genscript"
    belief_full, runs_v, phi_v, L = _setup_layout()
    T = EM.empirical_T(L, runs_v, phi_v, trials, np.random.default_rng(SEED + 7))
    m = belief_full.m
    tw = json.loads((RES / "paper_twist_belief_mean.json").read_text())
    b_twist = ChannelParams(np.full(m, tw["yield"]), np.full(m, tw["shape"]), np.full(m, tw["dropout"]), np.ones(m),
                            np.tile(np.array(tw["phi"]), (m, 1)), FailureModel(np.full(m, 1e-3)), tw["sigma_g"], 0.0)
    rows = []
    for name, b in [("in-pool calibration", belief_full), ("transferred from Twist pool", b_twist),
                    ("generic prior", generic_prior(m))]:
        for psd in (0.5, 1.0):
            r = _run(L, b, T, trials, prior_sd=psd)
            c = np.where(r.success, r.cost, np.inf)
            s = adaptive_summary(c, r.success, (0.95, 0.99))
            share = r.reads_per_container.mean(0) / r.reads_per_container.mean(0).sum()
            rows.append({"belief": name, "prior_sd": psd, "trials": trials, **_flat(s),
                         "read_share_broken_container": float(share[list(RD.RUNS.values()).index("PCR_30c")])})
            print(name, psd, round(s["mean_spend@0.99"]), rows[-1]["read_share_broken_container"], flush=True)
            pd.DataFrame(rows).to_csv(RES / "paper_e10_cross_pool.csv", index=False)


def e10_static(trials=60):
    """Static references on the SAME 60 Genscript recovery-time matrices as E10."""
    assert RD.POOL == "genscript", "run with D3RA_POOL=genscript"
    from d3ra.channel.surrogate import file_time, fragment_times
    from d3ra.evaluation.metrics import static_summary
    from d3ra.retrieval.fixed_policy import uniform_weights
    belief_full, runs_v, phi_v, L = _setup_layout()
    T = EM.empirical_T(L, runs_v, phi_v, trials, np.random.default_rng(SEED + 7))
    rows = []
    for name, w in [("fixed depth (uniform)", uniform_weights(L)), ("oracle-tuned static", RD._oracle_static_weights(L, T))]:
        s = static_summary(file_time(L, fragment_times(L, T, w)), (0.95, 0.99))
        rows.append({"policy": name, "trials": trials, **_flat(s)})
        print(name, round(s["mean_spend@0.99"]), flush=True)
    pd.DataFrame(rows).to_csv(RES / "paper_e10_static_reference.csv", index=False)


REPLAY_POLICIES = ["uniform", "coverage", "pivotal", "ts", "ucb", "ltc", "info", "tail_lcb"]


def replay_v3(trials=200):
    """Leak-free replay (v3 caches): every adaptive controller starts from the GENERIC prior.
    Static references: fixed depth (uniform) and hindsight-tuned static weights. The leak-free
    calibrated belief (calibration reads only) is run as a secondary condition. Per-trial
    costs are saved so that every comparison can be re-computed (paired over read orders)."""
    from d3ra.channel.surrogate import fragment_times, file_time
    from d3ra.retrieval.fixed_policy import uniform_weights, budget_for_target
    belief_cal, runs_v, phi_v, L = _setup_layout()
    T = EM.empirical_T(L, runs_v, phi_v, trials, np.random.default_rng(SEED + 7))
    gen = generic_prior(L.m)
    per_trial, rows = {}, []
    for name, w in [("static_uniform", uniform_weights(L)), ("static_hindsight", RD._oracle_static_weights(L, T))]:
        c = file_time(L, fragment_times(L, T, w))
        per_trial[name] = c
        s = static_summary(c, (0.95, 0.99))
        rows.append({"policy": name, "belief": "--", **_flat(s), "mean_reads": float(budget_for_target(c, 0.99)),
                     "q99_reads": float(budget_for_target(c, 0.99)), "weights": json.dumps(np.round(w, 4).tolist())})
        print(RD.POOL, name, round(rows[-1]["q99_reads"]), flush=True)
    for belief_name, b in [("generic", gen), ("calibrated_leakfree", belief_cal)]:
        for pol in (REPLAY_POLICIES if belief_name == "generic" else ["coverage"]):
            ctl = AdaptiveController(L, b, AdaptiveConfig(pol, batch=200))
            r = ctl.run(trials, np.random.default_rng(SEED + 11), T=T, failed=np.zeros((trials, L.m), bool))
            c = np.where(r.success, r.cost, np.inf)
            per_trial[f"{pol}|{belief_name}"] = c
            s = adaptive_summary(c, r.success, (0.95, 0.99))
            d = RD._paired_spend_diff(c, per_trial["static_hindsight"])
            rows.append({"policy": pol, "belief": belief_name, **_flat(s), "mean_reads": float(np.mean(c)),
                         "q99_reads": float(budget_for_target(c, 0.99)),
                         "diff_vs_hindsight_static": d[0], "diff_lo": d[1], "diff_hi": d[2],
                         "read_share": json.dumps((r.reads_per_container.mean(0) / r.reads_per_container.mean(0).sum()).round(4).tolist())})
            print(RD.POOL, pol, belief_name, round(rows[-1]["mean_reads"]), round(rows[-1]["q99_reads"]), flush=True)
    suf = "" if RD.POOL == "twist" else "_genscript"
    pd.DataFrame(rows).to_csv(RES / f"replay_v3{suf}.csv", index=False)
    pd.DataFrame({k: v for k, v in per_trial.items()}).to_csv(RES / f"replay_v3{suf}_per_trial.csv", index=False)


if __name__ == "__main__":
    for stage in sys.argv[1:]:
        {"e4": e4, "ec": ec, "ed": ed, "ed_od": ed_od, "ed_holdout": ed_holdout, "e10": e10, "e10_static": e10_static, "replay_v3": replay_v3}[stage]()
