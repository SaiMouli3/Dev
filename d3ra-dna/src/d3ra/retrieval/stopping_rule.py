"""Stopping rules.

* ``FixedBudgetStop``  : stop after a predetermined budget (classical fixed depth).
* ``DecodedStop``      : stop as soon as every stripe holds K inner-verified
  fragments. Confidence model: each fragment passed inner RS + CRC-c, so the
  probability that a stripe decodes to wrong data is bounded by the union of
  per-fragment false-accept probabilities (~K 2^-(8c)); the file-level
  SHA-256 then catches any residual error, in which case retrieval continues.
* ``FutilityStop``     : additionally abandon when the predicted probability of
  completing within the remaining budget falls below ``eps``.
"""
from __future__ import annotations

import numpy as np
from scipy.special import log_ndtr


def false_accept_bound(k_total: int, crc_bytes: int) -> float:
    return min(1.0, k_total * 2.0 ** (-8 * crc_bytes))


def completion_probability(rho: np.ndarray, need: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Normal-approx P(every stripe gets `need` more of its unrecovered
    fragments) where rho (..., S, N) are per-fragment probabilities of
    recovery within the remaining budget."""
    r = np.where(mask, rho, 0.0)
    mu = r.sum(-1)
    var = (r * (1 - r)).sum(-1) + 1e-9
    z = (mu - (need - 0.5)) / np.sqrt(var)
    logp = np.where(need > 0, log_ndtr(z), 0.0)
    return np.exp(logp.sum(-1))


class FixedBudgetStop:
    def __init__(self, budget: float):
        self.budget = budget

    def __call__(self, cost: np.ndarray, done: np.ndarray, **_) -> np.ndarray:
        return cost >= self.budget


class DecodedStop:
    def __init__(self, cap: float = np.inf):
        self.cap = cap

    def __call__(self, cost: np.ndarray, done: np.ndarray, **_) -> np.ndarray:
        return done | (cost >= self.cap)


class FutilityStop(DecodedStop):
    def __init__(self, cap: float, eps: float = 1e-3):
        super().__init__(cap)
        self.eps = eps

    def __call__(self, cost, done, p_complete=None, **_):
        stop = super().__call__(cost, done)
        if p_complete is not None:
            stop |= (~done) & (p_complete < self.eps)
        return stop
