"""Experiment runners: baselines across scenarios, hypotheses H1-H5,
ablations A-G, placement-optimiser comparison, redundancy trade-off,
replication-vs-parity under shared bias, and scalability.

All randomness is seeded; results are returned as lists of flat dicts and
written to results/ by experiments/run_all.py.
"""
from __future__ import annotations

import time
import tracemalloc

import numpy as np

from ..channel.failures import FailureModel
from ..channel.surrogate import ChannelParams, Layout, Predictor, sample_T, static_budget_samples
from ..optimization.budget_optimizer import min_budget, optimise_weights
from ..optimization.joint_optimizer import _layout, independent_optimise, joint_optimise, joint_optimise_adaptive
from ..optimization.placement_optimizer import PlacementProblem, compare_methods
from ..optimization.redundancy_optimizer import build_layout, risk_adaptive, uniform_parity
from ..retrieval.adaptive_policy import AdaptiveConfig, AdaptiveController
from ..retrieval.fixed_policy import StaticPolicy, budget_for_target, coverage_equalising_weights, uniform_weights
from ..storage.placement import cost_tiered, domain_aware, round_robin, unit_retrieval_cost
from . import scenarios as S
from .baselines import Design, designs, uniform_layout
from .metrics import adaptive_summary, static_summary
from .statistical_tests import (bootstrap_ci, cliffs_delta, cohens_dz, fisher_two_proportions, holm, mcnemar, paired_bootstrap_diff,
                                unpaired_bootstrap_diff, wilcoxon_paired)

TARGETS = (0.99, 0.999)


def _flat(d: dict, prefix: str = "") -> dict:
    out = {}
    for k, v in d.items():
        if isinstance(v, tuple):
            out[f"{prefix}{k}_lo"], out[f"{prefix}{k}_hi"] = v
        else:
            out[prefix + k] = v
    return out


def run_design(d: Design, sc: S.Scenario, trials: int, seed: int, cap: float, static_trials: int | None = None):
    """Return (summary dict, per-trial cost array with inf for failures, success array)."""
    rng = np.random.default_rng(seed)
    if d.policy == "static":
        pol = StaticPolicy(d.weights)
        s = pol.budget_samples(d.layout, sc.truth, static_trials or trials, rng)
        return s, np.isfinite(s), None
    ctl = AdaptiveController(d.layout, sc.belief, AdaptiveConfig(controller=d.policy, cap=cap), truth=sc.truth)
    T, failed = sample_T(d.layout, sc.truth, trials, rng, return_state=True)
    r = ctl.run(trials, rng, T=T, failed=failed)
    cost = np.where(r.success, r.cost, np.inf)
    return cost, r.success, r


def baseline_experiment(name: str, trials: int = 2000, static_trials: int = 20000, seed: int = 0,
                        calibrated: bool = True) -> list[dict]:
    sc = S.make(name, seed, calibrated)
    t0 = time.time()
    ds = designs(sc)
    t_design = time.time() - t0
    # reference budget: baseline A's 0.99 budget (x budget_scale for the tight-budget scenario)
    A = ds[0]
    sA, _, _ = run_design(A, sc, trials, seed, np.inf, static_trials)
    b99 = budget_for_target(sA, 0.99)
    if not np.isfinite(b99):   # A cannot reach 0.99: reference = its 0.99 quantile among recoverable trials
        b99 = float(np.quantile(sA[np.isfinite(sA)], 0.99))
    B_ref = b99 * sc.budget_scale
    cap = 4 * b99
    rows, per_trial = [], {}
    for d in ds:
        t = time.time()
        cost, succ, r = (sA, np.isfinite(sA), None) if d is A else run_design(d, sc, trials, seed, cap, static_trials)
        el = time.time() - t
        if d.policy == "static":
            summ = static_summary(cost, TARGETS, B_ref)
        else:
            summ = adaptive_summary(cost, succ, TARGETS, B_ref, r.rounds, r.containers_used, r.extra_after_probe)
            summ["mean_latency"] = float(r.latency.mean())
        per_trial[d.name] = cost
        n_mol = d.layout.placement.sum()
        rows.append({"scenario": name, "design": d.name, "variant": d.variant, "policy": d.policy,
                     "parity": sc.parity, "B_ref": B_ref, "runtime_s": el, "design_time_s": t_design,
                     "redundancy_overhead": n_mol / d.layout.k.sum() - 1, **_flat(summ)})
    # paired comparisons on common layouts (C vs A share layout but A is static -> unpaired)
    comps = [("C_uniform_adaptive", "A_uniform_fixed"), ("D_reliability_adaptive", "B_reliability_fixed"),
             ("E_joint_adaptive", "C_uniform_adaptive"), ("E_joint_adaptive", "D_reliability_adaptive"),
             ("F_oracle", "E_joint_adaptive")]
    pv = []
    for a, b in comps:
        ca, cb = per_trial[a], per_trial[b]
        pa = np.isfinite(ca) & (ca <= B_ref)
        pb = np.isfinite(cb) & (cb <= B_ref)
        paired = a == "F_oracle"          # F replays E's exact T matrices (same layout and seed)
        p = mcnemar(pa, pb) if paired else fisher_two_proportions(pa, pb)
        pv.append(p)
    adj = holm([p if np.isfinite(p) else 1.0 for p in pv])
    for (a, b), p, pa in zip(comps, pv, adj):
        ca, cb = per_trial[a], per_trial[b]
        capa = np.where(np.isfinite(ca), ca, cap)
        capb = np.where(np.isfinite(cb), cb, cap)
        m = min(len(capa), len(capb))
        diff = unpaired_bootstrap_diff(capa[:m], capb[:m]) if a[0] != "F" else paired_bootstrap_diff(capa[:m], capb[:m])
        rows.append({"scenario": name, "design": f"cmp:{a}-vs-{b}", "mean_cost_diff": diff[0], "mean_cost_diff_lo": diff[1],
                     "mean_cost_diff_hi": diff[2], "cliffs_delta": cliffs_delta(capa, capb),
                     "p_exact_ref_test_p": p, "test": "mcnemar" if a == "F_oracle" else "fisher", "p_holm": pa})
    return rows


# --------------------------------------------------------------------------- H1

def h1_model_accuracy(name: str = "s3_nonuniform_abundance", trials: int = 20000, seed: int = 0) -> list[dict]:
    """Predictive accuracy of channel models vs Monte-Carlo ground truth."""
    sc = S.make(name, seed)
    L = uniform_layout(sc)
    w = uniform_weights(L)
    rng = np.random.default_rng(seed)
    truth_samples = static_budget_samples(L, sc.truth, w, trials, rng)
    ch = sc.truth
    m = ch.m
    homog = ChannelParams(np.full(m, ch.yield_.mean()), np.full(m, ch.shape.mean()), np.full(m, ch.dropout.mean()),
                          ch.cost, np.tile(ch.phi.mean(0), (m, 1)), FailureModel(np.full(m, ch.failure.f.mean())),
                          ch.sigma_g, ch.d_syn)
    poisson = ChannelParams(ch.yield_, np.full(m, 1e4), ch.dropout, ch.cost, ch.phi, ch.failure, 0.0, ch.d_syn)
    nofail = ChannelParams(ch.yield_, ch.shape, ch.dropout, ch.cost, ch.phi, FailureModel(np.zeros(m)), ch.sigma_g, ch.d_syn)
    models = {"heterogeneous(full)": ch, "homogeneous-avg": homog, "poisson-no-overdispersion": poisson,
              "no-failure-model": nofail}
    budgets = np.quantile(truth_samples[np.isfinite(truth_samples)], [0.1, 0.3, 0.5, 0.7, 0.9, 0.97, 0.99])
    rows = []
    for mname, mch in models.items():
        pred = Predictor(L, mch)
        errs = []
        for B in budgets:
            p_mc = float((truth_samples <= B).mean())
            p_hat = pred.p_file(B * w)
            errs.append(abs(p_hat - p_mc))
            rows.append({"experiment": "H1", "model": mname, "budget": float(B), "p_mc": p_mc, "p_pred": float(p_hat)})
        B99_pred = min_budget(pred, w, 0.99)
        rows.append({"experiment": "H1-summary", "model": mname, "mae_p": float(np.mean(errs)), "max_abs_err": float(np.max(errs)),
                     "B99_pred": B99_pred, "B99_mc": budget_for_target(truth_samples, 0.99),
                     "B99_rel_err": B99_pred / budget_for_target(truth_samples, 0.99) - 1})
    return rows


# --------------------------------------------------------------------------- H5

def h5_robustness(trials: int = 2000, seed: int = 0) -> list[dict]:
    """Uncertainty-aware controllers vs a controller that trusts a point estimate."""
    rows = []
    for name in ["s7_channel_change", "s9_uncertain_estimates", "s10_mixed"]:
        sc = S.make(name, seed)
        L = _layout("tiered1", sc.n_stripes, sc.k, sc.parity, sc.belief, sc.domain)
        rng = np.random.default_rng(seed + 1)
        T, failed = sample_T(L, sc.truth, trials, rng, return_state=True)
        cap = 1e6
        for label, cfg, belief in [
            ("point-estimate (prior_sd=0.02)", AdaptiveConfig("pivotal", prior_sd=0.02, cap=cap), sc.belief),
            ("bayes (prior_sd=0.5)", AdaptiveConfig("pivotal", prior_sd=0.5, cap=cap), sc.belief),
            ("thompson (prior_sd=0.5)", AdaptiveConfig("thompson", prior_sd=0.5, cap=cap), sc.belief),
            ("bayes (prior_sd=1.0)", AdaptiveConfig("pivotal", prior_sd=1.0, cap=cap), sc.belief),
            ("oracle-params", AdaptiveConfig("pivotal", prior_sd=0.02, cap=cap), sc.truth),
        ]:
            ctl = AdaptiveController(L, belief, cfg, truth=sc.truth)
            r = ctl.run(trials, np.random.default_rng(seed + 2), T=T, failed=failed)
            c = np.where(r.success, r.cost, np.inf)
            rows.append({"experiment": "H5", "scenario": name, "controller": label,
                         **_flat(adaptive_summary(c, r.success, TARGETS))})
        # static allocation optimised on the (wrong) belief, evaluated on truth
        pred = Predictor(L, sc.belief)
        w = uniform_weights(L)
        B = min_budget(pred, w, 0.99)
        w = optimise_weights(pred, B, w0=w, iters=15)
        s = static_budget_samples(L, sc.truth, w, trials * 5, np.random.default_rng(seed + 3)) * float(sc.truth.cost @ w)
        rows.append({"experiment": "H5", "scenario": name, "controller": "static-optimised-on-belief",
                     "predicted_cap@0.99": float(min_budget(pred, w, 0.99)), **_flat(static_summary(s, TARGETS))})
    return rows


# --------------------------------------------------------------------------- ablations

def ablations(name: str = "s10_mixed", trials: int = 2000, seed: int = 0) -> list[dict]:
    sc = S.make(name, seed)
    rows = []
    res = joint_optimise_adaptive(sc.n_stripes, sc.k, sc.belief, 0.99, sc.domain, parities=[sc.parity])
    full_layout = res.layout
    variant = res.variant
    res_static = joint_optimise(sc.n_stripes, sc.k, sc.belief, 0.99, sc.domain, parities=[sc.parity], iters=3, p0=sc.parity)
    rng_seed = seed + 11

    def adaptive(layout, controller="pivotal", belief=None):
        ctl = AdaptiveController(layout, belief or sc.belief, AdaptiveConfig(controller, cap=1e6), truth=sc.truth)
        r = ctl.run(trials, np.random.default_rng(rng_seed), T=None)
        return np.where(r.success, r.cost, np.inf), r

    def static(layout, w):
        return static_budget_samples(layout, sc.truth, w, trials * 5, np.random.default_rng(rng_seed)) * float(sc.truth.cost @ w)

    c, r = adaptive(full_layout)
    rows.append({"ablation": "full (E)", "variant": variant, **_flat(adaptive_summary(c, r.success, TARGETS))})
    c2, r2 = adaptive(uniform_layout(sc))
    rows.append({"ablation": "A: -reliability-aware placement", **_flat(adaptive_summary(c2, r2.success, TARGETS))})
    from ..optimization.budget_optimizer import optimise_weights as _ow
    _pred = Predictor(full_layout, sc.belief)
    _w = uniform_weights(full_layout)
    _B = min_budget(_pred, _w, 0.99)
    w_full = _ow(_pred, _B, w0=_w, iters=15) if np.isfinite(_B) else _w
    s = static(full_layout, w_full)
    rows.append({"ablation": "B: -adaptive sequencing (static optimised weights)", **_flat(static_summary(s, TARGETS))})
    # D: placement optimised ignoring correlated failures (independent marginals only)
    if sc.truth.failure.domains:
        nodom = ChannelParams(sc.belief.yield_, sc.belief.shape, sc.belief.dropout, sc.belief.cost, sc.belief.phi,
                              FailureModel(sc.belief.failure.marginal()), sc.belief.sigma_g, sc.belief.d_syn)
        res_nd = joint_optimise_adaptive(sc.n_stripes, sc.k, nodom, 0.99, None, parities=[sc.parity])
        c3, r3 = adaptive(res_nd.layout)
        rows.append({"ablation": "D: -correlated failure modelling", "variant": res_nd.variant,
                     **_flat(adaptive_summary(c3, r3.success, TARGETS))})
    # E: fixed stopping -- commit the controller's mean allocation as a fixed budget
    w_ad = r.reads_per_container.mean(0)
    w_ad = w_ad / w_ad.sum()
    s2 = static(full_layout, w_ad)
    rows.append({"ablation": "E: fixed stopping (controller's mean allocation, committed)", **_flat(static_summary(s2, TARGETS))})
    c4, r4 = adaptive(full_layout, "oracle", sc.truth)
    rows.append({"ablation": "F: perfect channel information (oracle)", **_flat(adaptive_summary(c4, r4.success, TARGETS))})
    ind = independent_optimise(sc.n_stripes, sc.k, sc.belief, 0.99, sc.domain, parities=[sc.parity])
    c5, r5 = adaptive(ind.layout)
    rows.append({"ablation": "G: independent instead of joint optimisation", **_flat(adaptive_summary(c5, r5.success, TARGETS))})
    cs_, rs_ = adaptive(res_static.layout)
    rows.append({"ablation": "joint with static-budget proxy (v1) instead of adaptive-in-the-loop",
                 "variant": res_static.history[-1]["variant"], **_flat(adaptive_summary(cs_, rs_.success, TARGETS))})
    c6, r6 = adaptive(full_layout, "uniform")
    rows.append({"ablation": "adaptive stopping only (uniform round-robin allocation)", **_flat(adaptive_summary(c6, r6.success, TARGETS))})
    c7, r7 = adaptive(full_layout, "coverage")
    rows.append({"ablation": "coverage-greedy instead of pivotal", **_flat(adaptive_summary(c7, r7.success, TARGETS))})
    for row in rows:
        row["scenario"] = name
    return rows


# --------------------------------------------------------------------------- redundancy (Algorithm 3, H4 storage)

def redundancy_experiment(seed: int = 0, target: float = 0.999) -> list[dict]:
    """Uniform vs risk-adaptive parity with stripes spread over *subsets* of a
    heterogeneous population (stripe risk differs)."""
    sc = S.scale_scenario(24, 40, k=32, parity=0, seed=seed)
    ch = sc.truth
    m = ch.m
    groups = np.array_split(np.random.default_rng(seed).permutation(m), 6)   # each stripe lives in one 4-container group

    def place(stripe_of, k):
        X = np.zeros((len(stripe_of), m))
        for s in np.unique(stripe_of):
            idx = np.flatnonzero(stripe_of == s)
            g = groups[s % len(groups)]
            X[idx, g[np.arange(len(idx)) % len(g)]] = 1
        return X

    k = np.full(40, 32)
    rows = []
    for rpf in [6.0, 10.0]:
        p, Lu, pu = uniform_parity(k, ch, place, rpf, target)
        n_ad, La, pa = risk_adaptive(k, ch, place, rpf, target, start_parity=0)
        rows.append({"experiment": "redundancy", "reads_per_fragment": rpf, "method": "uniform", "parity_total": int(p * len(k)),
                     "overhead": p * len(k) / k.sum(), "p_file_pred": pu})
        rows.append({"experiment": "redundancy", "reads_per_fragment": rpf, "method": "risk-adaptive",
                     "parity_total": int((n_ad - k).sum()), "overhead": (n_ad - k).sum() / k.sum(), "p_file_pred": pa,
                     "parity_min": int((n_ad - k).min()), "parity_max": int((n_ad - k).max())})
        # MC validation of both
        for lab, L in [("uniform", Lu), ("risk-adaptive", La)]:
            R = rpf * L.placement.sum(0)
            s = static_budget_samples(L, ch, R / R.sum(), 20000, np.random.default_rng(seed))
            rows.append({"experiment": "redundancy-mc", "reads_per_fragment": rpf, "method": lab,
                         "p_file_mc": float((s <= R.sum()).mean())})
    return rows


# --------------------------------------------------------------------------- replication vs parity under shared bias

def replication_vs_parity(seed: int = 0, trials: int = 20000) -> list[dict]:
    """At equal stored molecules (2x), compare 2-way cross-container replication
    vs RS parity, as the share of abundance variance that is sequence-intrinsic
    (sigma_g) grows at fixed total log-variance."""
    rows = []
    base = S.make("s2_heterogeneous", seed)
    m = base.truth.m
    total_var = 0.6 ** 2
    for frac in [0.0, 0.25, 0.5, 0.75, 1.0]:
        sg = np.sqrt(frac * total_var)
        sh_var = (1 - frac) * total_var
        shape = 1 / max(np.exp(sh_var) - 1, 1e-4)        # gamma shape with same log-variance (approx)
        ch = ChannelParams(base.truth.yield_, np.full(m, shape), base.truth.dropout, np.ones(m), base.truth.phi,
                           base.truth.failure, sg, 0.0)
        K, S_ = 32, 20
        # replication: N=K data fragments, each copied into 2 containers
        so = np.repeat(np.arange(S_), K)
        Xr = np.zeros((len(so), m))
        for i in range(len(so)):
            Xr[i, i % m] = 1
            Xr[i, (i + 3) % m] = 1
        Lr = Layout(so, np.full(S_, K), Xr)
        so2 = np.repeat(np.arange(S_), 2 * K)
        Lp = Layout(so2, np.full(S_, K), round_robin(so2, m))
        for lab, L in [("replication-2x", Lr), ("rs-parity-2x", Lp)]:
            w = uniform_weights(L)
            s = static_budget_samples(L, ch, w, trials, np.random.default_rng(seed))
            rows.append({"experiment": "replication_vs_parity", "shared_variance_fraction": frac, "scheme": lab,
                         "stored_molecules": float(L.placement.sum()), **_flat(static_summary(s, (0.99,)))})
    return rows


# --------------------------------------------------------------------------- placement optimisers (gap)

def placement_optimizer_comparison(seed: int = 0) -> list[dict]:
    rows = []
    for name, ns in [("s2_heterogeneous", 4), ("s4_correlated_failures", 4), ("s10_mixed", 6)]:
        sc = S.make(name, seed)
        ch = sc.belief
        e = unit_retrieval_cost(ch)
        so = np.repeat(np.arange(ns), sc.k + sc.parity)
        states = ch.failure.enumerate_states(max_states=2 ** 12)
        if states is None:
            st = ch.failure.sample(512, np.random.default_rng(seed))
            states, weights = st.astype(float), np.ones(len(st))
        else:
            weights = np.array([p for _, p in states])
            states = np.stack([s for s, _ in states]).astype(float)
        prob = PlacementProblem(so, np.full(ns, sc.k), e, states, weights, delta=1.0, eps=1e-2,
                                domain=sc.domain)
        for r in compare_methods(prob, np.random.default_rng(seed), ilp_limit=60):
            rows.append({"experiment": "placement_opt", "scenario": name, "stripes": ns, **r})
    return rows


# --------------------------------------------------------------------------- scalability (scenario 8)

def scalability(seed: int = 0, adaptive_trials: int = 200) -> list[dict]:
    rows = []
    for m, ns in [(10, 20), (50, 100), (100, 200), (100, 600)]:
        sc = S.scale_scenario(m, ns, seed=seed)
        row = {"experiment": "scalability", "containers": m, "stripes": ns, "fragments": ns * (sc.k + sc.parity)}
        tracemalloc.start()
        t = time.time()
        L = _layout("tiered1", ns, sc.k, sc.parity, sc.truth, None)
        row["placement_s"] = time.time() - t
        t = time.time()
        pred = Predictor(L, sc.truth, n_states=32)
        w = uniform_weights(L)
        row["B99_pred"] = min_budget(pred, w, 0.99)
        row["predictor_bisection_s"] = time.time() - t
        t = time.time()
        s = static_budget_samples(L, sc.truth, w, 2000, np.random.default_rng(seed))
        row["mc_static_2000_s"] = time.time() - t
        row["B99_mc"] = budget_for_target(s, 0.99)
        t = time.time()
        ctl = AdaptiveController(L, sc.truth, AdaptiveConfig("pivotal", cap=1e7))
        r = ctl.run(adaptive_trials, np.random.default_rng(seed))
        row["adaptive_s"] = time.time() - t
        row["adaptive_trials"] = adaptive_trials
        row["adaptive_mean_cost"] = float(r.cost[r.success].mean()) if r.success.any() else np.inf
        row["peak_mem_mb"] = tracemalloc.get_traced_memory()[1] / 1e6
        tracemalloc.stop()
        rows.append(row)
    return rows


# --------------------------------------------------------------------------- spreading threshold

def spreading_threshold(name: str = "s4_correlated_failures", parities=(4, 8, 12, 16, 20, 24, 28, 32),
                        seed: int = 0, trials: int = 20000) -> list[dict]:
    """For each parity, exact stripe-survival probability (failures + dropout)
    of each placement family, plus a Monte-Carlo check of file-level infeasibility."""
    from ..channel.surrogate import file_time, fragment_times
    from ..storage.placement import survival_prob
    sc = S.make(name, seed)
    ch = sc.truth
    ex = ch.failure.enumerate_states(max_states=2 ** 12)
    states = np.stack([s_ for s_, _ in ex])
    weights = np.array([p_ for _, p_ in ex])
    presence = (1 - ch.dropout) * (1 - ch.d_syn)
    D = len(np.unique(sc.domain)) if sc.domain is not None else ch.m
    rows = []
    for p in parities:
        for v in ["round_robin", "domain", "tiered1", "failure0"]:
            if v == "round_robin":
                so = np.repeat(np.arange(sc.n_stripes), sc.k + p)
                L = Layout(so, np.full(sc.n_stripes, sc.k), round_robin(so, ch.m))
            else:
                L = _layout(v, sc.n_stripes, sc.k, p, ch, sc.domain)
            q = L.placement[L.stripe_of == 0].sum(0).astype(int)
            ps = survival_prob(q, sc.k, states, weights, presence)
            T = sample_T(L, ch, trials, np.random.default_rng(seed))
            infeas = float(np.isinf(file_time(L, fragment_times(L, T, np.ones(ch.m)))).mean())
            rows.append({"experiment": "spreading_threshold", "scenario": name, "parity": p,
                         "parity_fraction": p / (sc.k + p), "one_over_D": 1 / D, "placement": v,
                         "counts_stripe0": " ".join(map(str, q)), "stripe_survival_exact": ps,
                         "file_infeasible_mc": infeas})
    return rows
