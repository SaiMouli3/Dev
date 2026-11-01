"""Joint placement / redundancy / retrieval optimisation (alternating).

Step 1  placement given retrieval policy  (capped water-filling seed + SA on the
        decomposable objective with unit costs re-weighted by the current
        allocation's realised read shares)
Step 2  retrieval policy given placement  (static weights by exponentiated
        gradient; adaptive controller configured with those weights' horizon)
Step 3  redundancy given placement & policy (smallest uniform parity, or
        risk-adaptive per-stripe parity, meeting the target)
Repeat until the predicted objective J = alpha*C_seq + beta*C_storage stops
improving (reliability enters as a hard constraint P_file >= target).

Independent optimisation (the comparison baseline) performs each step once,
each with the *default* assumptions of the others (uniform placement for
budget optimisation, uniform reads for redundancy, etc.).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..channel.surrogate import ChannelParams, Layout, Predictor
from ..storage.placement import cost_tiered, domain_aware, failure_aware, unit_retrieval_cost
from .budget_optimizer import min_budget, optimise_weights


@dataclass
class JointResult:
    layout: Layout
    weights: np.ndarray
    parity: int
    budget: float
    p_file: float
    storage_overhead: float
    objective: float
    history: list = field(default_factory=list)


def _objective(budget, n_total, k_total, alpha, beta, ref_budget):
    return alpha * budget / ref_budget + beta * (n_total / k_total - 1)


def evaluate(layout: Layout, ch: ChannelParams, w: np.ndarray, target: float):
    pred = Predictor(layout, ch)
    B = min_budget(pred, w, target)
    return B, pred


PLACEMENT_VARIANTS = ("tiered0", "tiered1", "tiered2", "domain", "failure0", "failure1")


def _layout(variant: str, n_stripes: int, k: int, p: int, ch: ChannelParams, domain) -> Layout:
    m = ch.m
    stripe_of = np.repeat(np.arange(n_stripes), k + p)
    dom = domain if domain is not None else np.full(m, -1)
    if variant.startswith("tiered"):
        X = cost_tiered(stripe_of, np.full(n_stripes, k), unit_retrieval_cost(ch), dom,
                        slack=min(int(variant[-1]), max(p - 1, 0)))
    elif variant.startswith("failure"):
        presence = (1 - ch.dropout) * (1 - ch.d_syn)
        X = failure_aware(stripe_of, np.full(n_stripes, k), ch.failure, presence, unit_retrieval_cost(ch),
                          lam=0.0 if variant == "failure0" else 0.05)
    else:
        X = domain_aware(stripe_of, k, dom, m)
    return Layout(stripe_of, np.full(n_stripes, k), X)


def joint_optimise(n_stripes: int, k: int, ch: ChannelParams, target: float, domain=None,
                   parities=range(2, 25, 2), alpha: float = 1.0, beta: float = 0.5,
                   iters: int = 4, p0: int | None = None) -> JointResult | None:
    """Alternating (block-coordinate) optimisation over placement variant v,
    read weights w and parity p, minimising J subject to P_file >= target."""
    parities = list(parities)
    p = p0 if p0 is not None else parities[len(parities) // 2]
    v = "domain"
    w = None
    hist = []
    ref = [None]

    def score(v_, p_, w_):
        L = _layout(v_, n_stripes, k, p_, ch, domain)
        pred = Predictor(L, ch)
        ww = w_ if w_ is not None else (L.stored_per_container() > 0) / max((L.stored_per_container() > 0).sum(), 1)
        B = min_budget(pred, ww, target)
        if ref[0] is None and np.isfinite(B):
            ref[0] = B
        J = _objective(B, L.n, n_stripes * k, alpha, beta, ref[0] or 1.0) if np.isfinite(B) else np.inf
        return J, B, L, pred, ww

    J, B, L, pred, w = score(v, p, w)
    for it in range(iters):
        changed = False
        # step 1: placement | (p, w)
        cand = [(score(vv, p, w)[0], vv) for vv in PLACEMENT_VARIANTS]
        Jv, vv = min(cand)
        if Jv < J - 1e-9:
            v, changed = vv, True
        # step 2: retrieval weights | (v, p)
        J, B, L, pred, _ = score(v, p, w)
        if np.isfinite(B):
            w_new = optimise_weights(pred, B, w0=w, iters=15)
            Jw, Bw, *_ = score(v, p, w_new)
            if Jw < J - 1e-9:
                w, changed = w_new, True
        # step 3: redundancy | (v, w)
        cand = [(score(v, pp, w)[0], pp) for pp in parities]
        Jp, pp = min(cand)
        if Jp < score(v, p, w)[0] - 1e-9:
            p, changed = pp, True
        J, B, L, pred, w = score(v, p, w)
        hist.append({"iter": it, "variant": v, "parity": p, "budget": B, "J": J})
        if not changed:
            break
    if not np.isfinite(B):
        return None
    return JointResult(L, w, p, B, pred.p_file(B * w / (ch.cost @ w)), L.n / (n_stripes * k), J, hist)


def independent_optimise(n_stripes: int, k: int, ch: ChannelParams, target: float, domain=None,
                         parities=range(2, 25, 2), alpha: float = 1.0, beta: float = 0.5) -> JointResult:
    """Each component optimised under default assumptions of the others:
    redundancy chosen for uniform placement + uniform reads, placement chosen by
    reliability only, read weights optimised last."""
    m = ch.m
    dom = domain if domain is not None else np.full(m, -1)
    chosen = None
    for p in parities:
        stripe_of = np.repeat(np.arange(n_stripes), k + p)
        L = Layout(stripe_of, np.full(n_stripes, k), domain_aware(stripe_of, k, dom, m))
        pred = Predictor(L, ch)
        w = np.ones(m) / m
        B = min_budget(pred, w, target)
        if np.isfinite(B):
            J = _objective(B, L.n, n_stripes * k, alpha, beta, B)
            if chosen is None or J < chosen[1]:
                chosen = (p, J, B)
    if chosen is None:
        return None
    p = chosen[0]
    stripe_of = np.repeat(np.arange(n_stripes), k + p)
    L = Layout(stripe_of, np.full(n_stripes, k), domain_aware(stripe_of, k, dom, m))
    pred = Predictor(L, ch)
    w = np.ones(m) / m
    B = min_budget(pred, w, target)
    w = optimise_weights(pred, B, w0=w, iters=15)
    B = min_budget(pred, w, target)
    return JointResult(L, w, p, B, pred.p_file(B * w / (ch.cost @ w)), L.n / (n_stripes * k),
                       _objective(B, L.n, n_stripes * k, alpha, beta, chosen[2]))


# --------------------------------------------------------------------------- simulation-in-the-loop joint optimisation

ADAPTIVE_VARIANTS = ("reliability", "tiered0", "tiered1", "tiered2", "domain", "failure0", "failure1")


def _reliability_scores(ch: ChannelParams) -> np.ndarray:
    return (1 - ch.failure.marginal()) * (1 - ch.dropout) * ch.yield_ * ch.phi[:, 3]


def layout_for(variant: str, n_stripes: int, k: int, p: int, ch: ChannelParams, domain) -> Layout:
    if variant == "reliability":
        from ..storage.placement import reliability_weighted
        stripe_of = np.repeat(np.arange(n_stripes), k + p)
        return Layout(stripe_of, np.full(n_stripes, k), reliability_weighted(stripe_of, _reliability_scores(ch)))
    return _layout(variant, n_stripes, k, p, ch, domain)


def joint_optimise_adaptive(n_stripes: int, k: int, belief: ChannelParams, target: float, domain=None,
                            parities=(16,), controller: str = "pivotal", trials: int = 300,
                            feas_trials: int = 4000, beta: float = 0.0, seed: int = 12345,
                            variants=ADAPTIVE_VARIANTS) -> JointResult | None:
    """Joint placement x redundancy optimisation with the *adaptive retrieval
    policy in the loop*: each candidate (variant, parity) is scored by simulating
    the controller on the **believed** channel (never the truth):

        feasible  iff  P_infeasible(belief) <= 1 - target   (MC on the T matrix)
        J = E[adaptive spend | belief] / ref  +  beta * storage overhead

    Infeasible candidates are only chosen if no candidate is feasible (then the
    one with the smallest infeasibility wins). Common random numbers are used
    across candidates."""
    from ..channel.surrogate import file_time, fragment_times, sample_T
    from ..retrieval.adaptive_policy import AdaptiveConfig, AdaptiveController
    hist = []
    best = None
    ref = None
    for p in parities:
        for v in variants:
            L = layout_for(v, n_stripes, k, p, belief, domain)
            T = sample_T(L, belief, feas_trials, np.random.default_rng(seed))
            p_inf = float(np.isinf(file_time(L, fragment_times(L, T, np.ones(L.m)))).mean())
            ctl = AdaptiveController(L, belief, AdaptiveConfig(controller, cap=np.inf))
            r = ctl.run(trials, np.random.default_rng(seed + 1))
            spend = float(r.cost[r.success].mean()) if r.success.any() else np.inf
            ref = ref or spend
            J = spend / ref + beta * (L.placement.sum() / (n_stripes * k) - 1)
            feasible = p_inf <= 1 - target
            hist.append({"variant": v, "parity": p, "p_infeasible": p_inf, "spend": spend, "J": J, "feasible": feasible})
            key = (not feasible, p_inf if not feasible else 0.0, J)
            if best is None or key < best[0]:
                best = (key, v, p, L, spend, p_inf)
    if best is None:
        return None
    _, v, p, L, spend, p_inf = best
    res = JointResult(L, (L.stored_per_container() > 0) / max((L.stored_per_container() > 0).sum(), 1), p, spend,
                      1 - p_inf, L.placement.sum() / (n_stripes * k), best[0][2], hist)
    res.variant = v
    return res
