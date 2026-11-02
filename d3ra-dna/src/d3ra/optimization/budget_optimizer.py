"""Static sequencing-budget optimisation on the analytical predictor.

maximise  log P_file(R)   s.t.  sum_j c_j R_j = B,  R >= 0
solved by exponentiated-gradient (mirror ascent) on the cost simplex; and the
constrained dual problem  min B  s.t. P_file(B w*) >= target  by bisection.
"""
from __future__ import annotations

import numpy as np

from ..channel.surrogate import ChannelParams, Layout, Predictor


def optimise_weights(pred: Predictor, budget: float, iters: int = 40, lr: float = 0.5,
                     w0: np.ndarray | None = None, nominal_only: bool = True) -> np.ndarray:
    """Return read fractions w (sum 1) maximising predicted P_file at cost `budget`."""
    ch, L = pred.ch, pred.layout
    used = L.stored_per_container() > 0
    c = ch.cost
    # work with cost shares v (sum 1): R_j = budget * v_j / c_j
    v = (np.ones(L.m) * used) if w0 is None else np.asarray(w0) * c
    v = v / v.sum()
    states = pred.states
    if nominal_only:
        pred.states = [(np.zeros(L.m, bool), 1.0)]
    try:
        def f(vv):
            return np.log(max(pred.p_file(budget * vv / c), 1e-300))
        cur = f(v)
        for _ in range(iters):
            g = np.zeros(L.m)
            for j in np.flatnonzero(used):
                d = 0.02 * max(v[j], 0.01)
                vj = v.copy()
                vj[j] += d
                vj /= vj.sum()
                g[j] = (f(vj) - cur) / d
            g = np.where(used, g, -np.inf)
            nv = v * np.exp(lr * (g - g[used].max()) / (np.abs(g[used]).max() + 1e-12))
            nv = np.where(used, nv, 0)
            nv /= nv.sum()
            nf = f(nv)
            if nf >= cur:
                v, cur = nv, nf
            else:
                lr *= 0.5
                if lr < 1e-3:
                    break
    finally:
        pred.states = states
    R = budget * v / c
    return R / R.sum()


def min_budget(pred: Predictor, w: np.ndarray, target: float, lo: float = 1.0, hi: float | None = None,
               tol: float = 0.01) -> float:
    """Smallest cost B with predicted P_file(B * w / (c.w)) >= target (inf if unreachable)."""
    c = pred.ch.cost
    cw = float(c @ w)
    f = lambda B: pred.p_file(B * w / cw)
    hi = hi or 1000.0 * pred.layout.n
    if f(hi) < target:
        return np.inf
    while hi / lo > 1 + tol:
        mid = np.sqrt(lo * hi)
        if f(mid) >= target:
            hi = mid
        else:
            lo = mid
    return hi


def optimise_for_target(pred: Predictor, target: float, rounds: int = 2) -> tuple[np.ndarray, float]:
    """Alternate weight optimisation and budget bisection."""
    L = pred.layout
    w = (L.stored_per_container() > 0).astype(float)
    w /= w.sum()
    B = min_budget(pred, w, target)
    for _ in range(rounds):
        if not np.isfinite(B):
            break
        w = optimise_weights(pred, B, w0=w)
        B = min_budget(pred, w, target)
    return w, B
