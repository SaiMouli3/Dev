"""Static (non-adaptive) sequencing allocations.

A static policy commits read fractions w_j (summing to 1 over reads) and a
total budget B before observing anything. Its exact per-trial budget-to-
recover is computed on the recovery-time matrix (see channel.surrogate)."""
from __future__ import annotations

import numpy as np

from ..channel.surrogate import ChannelParams, Layout, static_budget_samples


def uniform_weights(layout: Layout) -> np.ndarray:
    m = layout.m
    used = layout.stored_per_container() > 0
    w = used / used.sum()
    return w.astype(float)


def coverage_equalising_weights(layout: Layout, ch: ChannelParams) -> np.ndarray:
    """w_j ∝ n_j * E[t_j] / y_j : equal expected *decodable* coverage per
    fragment in every container (reliability-weighted static allocation)."""
    n_j = layout.placement.sum(axis=0)
    pmf = np.diff(ch.phi, axis=1, prepend=0.0)
    et = (pmf * np.arange(ch.phi.shape[1])).sum(axis=1) + (1 - ch.phi[:, -1]) * ch.phi.shape[1]
    w = n_j * et / np.maximum(ch.yield_, 1e-6)
    w = np.where(n_j > 0, w, 0)
    return w / w.sum()


class StaticPolicy:
    def __init__(self, w: np.ndarray, name: str = "static"):
        self.w = np.asarray(w, float) / np.sum(w)
        self.name = name

    def budget_samples(self, layout: Layout, ch: ChannelParams, trials: int, rng, cost: bool = True):
        """Per-trial minimal spend (cost units if cost=True else reads)."""
        reads = static_budget_samples(layout, ch, self.w, trials, rng)
        return reads * (float(ch.cost @ self.w) if cost else 1.0)


def reliability_at(samples: np.ndarray, budget: float) -> float:
    return float((samples <= budget).mean())


def budget_for_target(samples: np.ndarray, target: float) -> float:
    """Smallest budget whose empirical reliability >= target (inf if unreachable)."""
    s = np.sort(samples)
    k = int(np.ceil(target * len(s))) - 1
    return float(s[min(max(k, 0), len(s) - 1)]) if np.isfinite(s[min(max(k, 0), len(s) - 1)]) else np.inf
