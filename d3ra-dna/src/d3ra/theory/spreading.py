"""Spreading vs concentration of an erasure-coded stripe across failure domains.

Setting: one (N, K) stripe; D failure domains; domain d fails (loses all its
fragments) independently with probability q_d; fragments may additionally be
erased independently with probability delta (dropout). Placement
n = (n_1..n_D), sum n_d = N. MDS: the stripe survives iff the surviving,
non-dropped fragments number >= K.

Results implemented here (proofs in docs/theory_spreading.md):
  * loss_probability(n, K, q, delta)        exact (enumeration over failure sets)
  * homogeneous_crossover(D, t)            unique q* (Theorem 1b)
  * tolerance(n, K)                        max #domain failures always survived
  * no_tolerance_dominated(...)            Proposition 2 check
"""
from __future__ import annotations

import itertools
import math

import numpy as np
from scipy import optimize, stats


def loss_probability(n, K: int, q, delta: float = 0.0, k_eff: float | None = None) -> float:
    """Exact P(stripe lost). q scalar or per-domain array. k_eff overrides K
    (e.g. K(1+eps) for non-MDS fountain codes)."""
    n = np.asarray(n, int)
    D = len(n)
    q = np.broadcast_to(np.asarray(q, float), (D,))
    need = K if k_eff is None else k_eff
    N = int(n.sum())
    total = 0.0
    for fail in itertools.product([0, 1], repeat=D):
        f = np.array(fail, bool)
        p = float(np.prod(np.where(f, q, 1 - q)))
        if p == 0:
            continue
        alive = int(n[~f].sum())
        if delta > 0:
            surv = 1 - stats.binom.cdf(math.ceil(need) - 1, alive, 1 - delta)
        else:
            surv = 1.0 if alive >= need else 0.0
        total += p * (1 - surv)
    return float(total)


def even_split(N: int, D: int) -> np.ndarray:
    base = np.full(D, N // D)
    base[: N % D] += 1
    return base


def tolerance(n, K: int) -> int:
    """Largest t such that ANY t domain failures leave >= K fragments."""
    s = np.sort(np.asarray(n))[::-1]
    N = s.sum()
    t = 0
    while t < len(s) and N - s[: t + 1].sum() >= K:
        t += 1
    return t


def homogeneous_spread_loss(D: int, t: int, q):
    """P(Bin(D,q) >= t+1): loss of an even split tolerating exactly t failures."""
    return stats.binom.sf(t, D, q)


def homogeneous_crossover(D: int, t: int) -> float:
    """Unique q* in (0,1] with P(Bin(D,q) >= t+1) = q (Theorem 1b); spreading is
    strictly better than concentration iff q < q*. t = 0 -> q* = 0 (never better)."""
    if t <= 0:
        return 0.0
    if t >= D - 1:
        return 1.0
    g = lambda q: homogeneous_spread_loss(D, t, q) - q
    # g < 0 near 0 (slope -1), g > 0 just below 1 (h > 1 there); bracket the sign change
    grid = np.linspace(1e-9, 1 - 1e-9, 20001)
    vals = g(grid)
    idx = np.flatnonzero((vals[:-1] < 0) & (vals[1:] >= 0))
    return float(optimize.brentq(g, grid[idx[0]], grid[idx[0] + 1]))


def threshold_parity_fraction(D: int) -> float:
    """Minimum parity fraction rho = (N-K)/N for an even split to tolerate one failure."""
    return 1.0 / D


def dropout_adjusted_parity(N: int, K: int, D: int, delta: float, z: float = 3.0) -> float:
    """Effective parity fraction after reserving a dropout margin: the even split
    tolerates one domain loss with high probability iff
    N(1 - 1/D)(1 - delta) - z*sqrt(N(1-1/D) delta(1-delta)) >= K."""
    alive = N * (1 - 1 / D)
    return (alive * (1 - delta) - z * math.sqrt(alive * delta * (1 - delta)) - K)


def no_tolerance_dominated(n, K: int, q) -> tuple[bool, float, float]:
    """Proposition 2: if some domain holds more than N-K fragments, the placement's
    loss probability is >= q of that domain >= min_d q_d (= concentration in the
    most reliable domain). Returns (holds, P_loss(n), min q)."""
    n = np.asarray(n)
    q = np.broadcast_to(np.asarray(q, float), (len(n),))
    N = n.sum()
    heavy = np.flatnonzero(n > N - K)
    pl = loss_probability(n, K, q)
    if heavy.size == 0:
        return True, pl, float(q.min())
    return bool(pl >= q[heavy].max() - 1e-12 and pl >= q.min() - 1e-12), pl, float(q.min())


def best_tolerant_vs_concentration(N: int, K: int, q, capacity=None, max_support: int | None = None):
    """Heterogeneous generalisation: among placements that tolerate every single
    domain failure (max n_d <= N-K, capacities respected), find the best by
    greedy-then-local search on exact loss; compare with concentration in the
    most reliable domain. Returns dict."""
    q = np.asarray(q, float)
    D = len(q)
    cap = np.full(D, N) if capacity is None else np.asarray(capacity)
    c = N - K
    order = np.argsort(q)
    best = None
    for r in range(1, D + 1):
        S = order[:r]
        lim = np.minimum(cap[S], c)
        if lim.sum() < N:
            continue
        n = np.zeros(D, int)
        # water-fill evenly inside S respecting limits
        rem = N
        while rem > 0:
            for j in S:
                if rem and n[j] < min(cap[j], c):
                    n[j] += 1
                    rem -= 1
        pl = loss_probability(n, K, q)
        if best is None or pl < best[1]:
            best = (n.copy(), pl, r)
    # most-concentrated feasible placement: fill the most reliable domains to capacity
    nm = np.zeros(D, int)
    rem = N
    for j in order:
        take = min(cap[j], rem)
        nm[j] = take
        rem -= take
    conc = loss_probability(nm, K, q) if rem == 0 else float("nan")
    return {"best_tolerant": None if best is None else best[0].tolist(),
            "p_loss_tolerant": None if best is None else best[1],
            "support": None if best is None else best[2],
            "p_loss_concentration": conc, "concentrated_placement": nm.tolist(),
            "min_support_needed": math.ceil(N / c) if c > 0 else None}
