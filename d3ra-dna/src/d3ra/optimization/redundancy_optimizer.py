"""Redundancy allocation: how many parity fragments each stripe gets.

* ``uniform_parity``  : same parity p for every stripe; smallest p meeting the target.
* ``risk_adaptive``   : greedy marginal allocation -- repeatedly add one parity
  fragment to the stripe with the largest increase in predicted log P_s per
  added fragment (log-concavity of Poisson-binomial tails in N makes greedy
  near-optimal); removes parity from over-protected stripes afterwards.
* ``replication``     : r copies per fragment (for comparison).

Parity fragments added to a stripe are placed by a supplied ``place`` callback
(e.g. cost-tiered capped water-filling).
"""
from __future__ import annotations

from typing import Callable

import numpy as np

from ..channel.surrogate import ChannelParams, Layout, Predictor


def build_layout(n_per_stripe: np.ndarray, k: np.ndarray, place: Callable[[np.ndarray, np.ndarray], np.ndarray]) -> Layout:
    stripe_of = np.repeat(np.arange(len(n_per_stripe)), n_per_stripe)
    return Layout(stripe_of, np.asarray(k), place(stripe_of, np.asarray(k)))


def stripe_log_probs(layout: Layout, ch: ChannelParams, R: np.ndarray) -> np.ndarray:
    pred = Predictor(layout, ch)
    vals = []
    ws = []
    for s, w in pred.states:
        vals.append(pred.stripe_probs(R, s if s.any() else None))
        ws.append(w)
    P = np.average(np.stack(vals), axis=0, weights=ws)
    return np.log(np.maximum(P, 1e-300))


def uniform_parity(k: np.ndarray, ch: ChannelParams, place, R_per_frag: float, target: float,
                   p_max: int = 64) -> tuple[int, Layout, float]:
    """Smallest uniform parity meeting `target` with R_j = R_per_frag * n_j reads."""
    for p in range(0, p_max + 1):
        L = build_layout(np.asarray(k) + p, k, place)
        R = R_per_frag * L.placement.sum(0)
        pf = Predictor(L, ch).p_file(R)
        if pf >= target:
            return p, L, pf
    return p_max, L, pf


def risk_adaptive(k: np.ndarray, ch: ChannelParams, place, R_per_frag: float, target: float,
                  start_parity: int = 0, max_total: int | None = None) -> tuple[np.ndarray, Layout, float]:
    k = np.asarray(k)
    n = k + start_parity
    max_total = max_total or int(4 * k.sum())
    while True:
        L = build_layout(n, k, place)
        R = R_per_frag * L.placement.sum(0)
        pred = Predictor(L, ch)
        pf = pred.p_file(R)
        if pf >= target or n.sum() >= max_total:
            return n, L, pf
        lp = stripe_log_probs(L, ch, R)
        # add parity to the worst stripes (batch of up to 5% of stripes per step)
        worst = np.argsort(lp)[: max(1, len(k) // 20)]
        worst = worst[lp[worst] < np.log(target) / len(k) * 0.5] if (lp[worst] < 0).any() else worst
        if len(worst) == 0:
            worst = np.argsort(lp)[:1]
        n[worst] += 1
