"""Fragment placement strategies across containers.

All strategies return an (n_fragments, M) copy-count matrix X with
sum_j X_ij >= 1, respecting container capacities.

* ``random``              : uniform random container per fragment.
* ``round_robin``         : fragment i -> container (i mod M), stripe-interleaved.
* ``reliability_weighted``: per stripe, fragments allocated in proportion to a
                            scalar reliability score (no failure-domain cap) --
                            the "intuitive" baseline that concentrates data in
                            the best containers.
* ``domain_aware``        : spreads each stripe as evenly as possible across failure
                            domains, then containers; any single domain loses <= N-K
                            fragments whenever feasible.
* ``cost_tiered``         : (proposed) capped water-filling -- per stripe, fill the
                            cheapest-to-read containers first (unit retrieval cost
                            e_j), but never put more than ``cap`` fragments of a
                            stripe into one container or failure domain (cap = N-K
                            by default, so any single container/domain failure is
                            survivable). Cheap containers become a primary read set;
                            expensive ones hold the fallback remainder.
"""
from __future__ import annotations

import numpy as np


def unit_retrieval_cost(ch) -> np.ndarray:
    """e_j = c_j * E[t_j] * (1 + 1/k_j) / (y_j (1-d_j)): expected cost to make one
    stored fragment decodable from container j (dispersion-penalised)."""
    pmf = np.diff(ch.phi, axis=1, prepend=0.0)
    et = (pmf * np.arange(ch.phi.shape[1])).sum(1) + (1 - ch.phi[:, -1]) * ch.phi.shape[1]
    return ch.cost * et * (1 + 1 / ch.shape) / (np.maximum(ch.yield_, 1e-6) * (1 - ch.dropout))


def _check_capacity(X, capacity):
    if capacity is not None and ((X > 0).sum(0) > np.asarray(capacity)).any():
        raise OverflowError("placement exceeds container capacity")
    return X


def random_placement(stripe_of, m, rng, capacity=None):
    n = len(stripe_of)
    X = np.zeros((n, m))
    cap = np.full(m, np.inf) if capacity is None else np.asarray(capacity, float).copy()
    for i in range(n):
        ok = np.flatnonzero(cap > 0)
        j = rng.choice(ok)
        X[i, j] = 1
        cap[j] -= 1
    return X


def round_robin(stripe_of, m, capacity=None, offset_per_stripe: bool = True):
    n = len(stripe_of)
    X = np.zeros((n, m))
    pos = np.zeros(stripe_of.max() + 1, dtype=int)
    for i, s in enumerate(stripe_of):
        j = (pos[s] + (s if offset_per_stripe else 0)) % m
        X[i, j] = 1
        pos[s] += 1
    return _check_capacity(X, capacity)


def _largest_remainder(total: int, weights: np.ndarray, cap: np.ndarray | None = None) -> np.ndarray:
    w = np.asarray(weights, float)
    w = w / w.sum()
    q = np.floor(total * w).astype(int)
    if cap is not None:
        q = np.minimum(q, cap)
    rem = total - q.sum()
    order = np.argsort(-(total * w - q))
    k = 0
    while rem > 0:
        j = order[k % len(order)]
        if cap is None or q[j] < cap[j]:
            q[j] += 1
            rem -= 1
        k += 1
        if k > 100 * len(order):
            raise ValueError("infeasible caps")
    return q


def _assign_counts(stripe_of, counts_per_stripe):
    n = len(stripe_of)
    m = counts_per_stripe.shape[1]
    X = np.zeros((n, m))
    for s in range(counts_per_stripe.shape[0]):
        idx = np.flatnonzero(stripe_of == s)
        cols = np.repeat(np.arange(m), counts_per_stripe[s])
        X[idx, cols[: len(idx)]] = 1
    return X


def reliability_weighted(stripe_of, scores, capacity=None):
    S = stripe_of.max() + 1
    n_s = np.bincount(stripe_of)
    counts = np.stack([_largest_remainder(n_s[s], scores) for s in range(S)])
    return _check_capacity(_assign_counts(stripe_of, counts), capacity)


def _capped_fill(n: int, order: np.ndarray, cap_c: np.ndarray, domain: np.ndarray, cap_d: int) -> np.ndarray:
    """Place n fragments one at a time into the cheapest container whose
    container cap and failure-domain cap are not yet reached; relax the caps
    only if the instance is otherwise infeasible."""
    m = len(cap_c)
    cap_c = np.asarray(cap_c, dtype=int).copy()
    q = np.zeros(m, dtype=int)
    dom_used: dict[int, int] = {}
    for _ in range(n):
        for _relax in range(4 * n + 4):
            pick = -1
            for j in order:
                d = domain[j]
                if q[j] < cap_c[j] and (d < 0 or dom_used.get(d, 0) < cap_d):
                    pick = j
                    break
            if pick >= 0:
                break
            cap_d += 1
            if cap_d > n:
                cap_c += 1
        q[pick] += 1
        if domain[pick] >= 0:
            dom_used[domain[pick]] = dom_used.get(domain[pick], 0) + 1
    return q


def domain_aware(stripe_of, k, domain, m, capacity=None):
    """Even spread over domains then containers (single-domain-failure tolerant when feasible)."""
    domain = np.asarray(domain)
    dom = np.where(domain < 0, np.arange(m) + 10**6, domain)
    S = stripe_of.max() + 1
    n_s = np.bincount(stripe_of)
    counts = []
    uniq = np.unique(dom)
    for s in range(S):
        per_dom = _largest_remainder(n_s[s], np.ones(len(uniq)))
        q = np.zeros(m, dtype=int)
        rot = s % len(uniq)
        per_dom = np.roll(per_dom, rot)
        for d, c in zip(uniq, per_dom):
            members = np.flatnonzero(dom == d)
            members = np.roll(members, s % len(members))
            q[members] = _largest_remainder(c, np.ones(len(members))) if c else 0
        counts.append(q)
    return _check_capacity(_assign_counts(stripe_of, np.stack(counts)), capacity)


def cost_tiered(stripe_of, k, unit_cost, domain=None, cap: np.ndarray | int | None = None,
                capacity=None, slack: int = 0):
    """Capped water-filling by unit retrieval cost (proposed placement).

    cap: per-stripe per-container cap (default N_s-K_s, i.e. single-failure tolerance).
    slack: reserve this many extra fragments beyond the cap (cap = N-K-slack), increasing
           robustness to a failure combined with coverage shortfalls.
    """
    m = len(unit_cost)
    domain = np.full(m, -1) if domain is None else np.asarray(domain)
    S = stripe_of.max() + 1
    n_s = np.bincount(stripe_of)
    k = np.broadcast_to(np.asarray(k), (S,))
    order = np.argsort(unit_cost, kind="stable")
    load = np.zeros(m, dtype=int)
    capv = np.full(m, np.inf) if capacity is None else np.asarray(capacity, float)
    counts = []
    for s in range(S):
        c = max(1, int(n_s[s] - k[s] - slack)) if cap is None else int(np.broadcast_to(cap, (S,))[s])
        cap_c = np.minimum(np.full(m, c), (capv - load)).astype(int)
        q = _capped_fill(n_s[s], order, np.maximum(cap_c, 0), domain, c)
        load += q
        counts.append(q)
    return _check_capacity(_assign_counts(stripe_of, np.stack(counts)), capacity)


def add_replicas(X: np.ndarray, frag_idx: np.ndarray, unit_cost: np.ndarray, copies: int = 1) -> np.ndarray:
    """Replicate selected fragments into the cheapest containers that don't hold them."""
    X = X.copy()
    order = np.argsort(unit_cost)
    for i in frag_idx:
        added = 0
        for j in order:
            if X[i, j] == 0:
                X[i, j] = 1
                added += 1
                if added == copies:
                    break
    return X


def stripe_failure_exposure(X, stripe_of, k, failure_states, weights) -> float:
    """P(some stripe loses > N-K fragments) under (sampled or enumerated) failure states."""
    S = stripe_of.max() + 1
    present = (X > 0).astype(float)
    # alive fragments per stripe per state: fragment alive if any copy in an alive container
    F = np.asarray(failure_states, float)              # (F, M)
    alive_frag = (present @ (1 - F).T) > 0             # (n, F)
    alive_s = np.zeros((S, F.shape[0]))
    np.add.at(alive_s, stripe_of, alive_frag)
    ok = (alive_s >= np.asarray(k)[:, None]).all(0)
    w = np.asarray(weights, float)
    return float(1 - (ok * w).sum() / w.sum())


def survival_prob(q: np.ndarray, k: int, states: np.ndarray, weights: np.ndarray, presence: np.ndarray) -> float:
    """P(stripe keeps >= k present fragments) for per-container counts q,
    averaged over failure states (rows of ``states``, True = failed); a fragment
    in an alive container j is present w.p. presence_j. Exact: DP over
    containers with binomial blocks, vectorised across states."""
    from scipy.stats import binom
    q = np.asarray(q, dtype=int)
    N = int(q.sum())
    F = states.shape[0]
    dp = np.zeros((F, N + 1))
    dp[:, 0] = 1.0
    for j in np.flatnonzero(q > 0):
        b = binom.pmf(np.arange(q[j] + 1), q[j], presence[j])
        new = np.zeros_like(dp)
        for t, bt in enumerate(b):
            new[:, t:] += bt * dp[:, : N + 1 - t]
        alive = ~states[:, j].astype(bool)
        dp = np.where(alive[:, None], new, dp)
    ok = dp[:, k:].sum(1)
    return float((ok * weights).sum() / weights.sum())


def failure_aware(stripe_of, k, failure, presence: np.ndarray, unit_cost: np.ndarray | None = None,
                  lam: float = 0.0, capacity=None, max_states: int = 4096, n_samples: int = 1024, seed: int = 0,
                  max_evals: int = 600):
    """Failure-exposure-minimising placement (proposed).

    Chooses per-stripe container counts q maximising
        log P(stripe survives failures and dropout)  -  lam * mean unit retrieval cost
    (survival exact over enumerated container + domain failure states, or a fixed
    sample when enumeration is too large). Search: start from 'spread evenly over
    the r most reliable containers' for every r, keep the best, then
    best-improvement moves of single fragments between containers."""
    m = len(presence)
    ex = failure.enumerate_states(max_states=max_states)
    if ex is None:
        st = failure.sample(n_samples, np.random.default_rng(seed))
        states, weights = st, np.ones(len(st))
    else:
        states = np.stack([s for s, _ in ex])
        weights = np.array([p for _, p in ex])
    uc = np.ones(m) if unit_cost is None else np.asarray(unit_cost, float)
    uc = uc / uc.min()
    marg = failure.marginal()
    rel_order = np.lexsort((uc, -(1 - marg) * presence))
    S = stripe_of.max() + 1
    n_s = np.bincount(stripe_of)
    k = np.broadcast_to(np.asarray(k), (S,))

    def score(q, N, K):
        return np.log(max(survival_prob(q, K, states, weights, presence), 1e-300)) - lam * float(q @ uc) / N

    cache = {}
    counts = []
    for s in range(S):
        key = (int(n_s[s]), int(k[s]))
        if key not in cache:
            N, K = key
            best_q, best_v = None, -np.inf
            for r in range(1, m + 1):
                q = np.zeros(m, dtype=int)
                q[rel_order[:r]] = _largest_remainder(N, np.ones(r))
                v = score(q, N, K)
                if v > best_v:
                    best_q, best_v = q, v
            q, v, evals = best_q.copy(), best_v, 0
            improved = True
            while improved and evals < max_evals:
                improved = False
                cand = None
                for a_ in np.flatnonzero(q > 0):
                    for b_ in range(m):
                        if a_ == b_:
                            continue
                        q[a_] -= 1
                        q[b_] += 1
                        vv = score(q, N, K)
                        evals += 1
                        if vv > v + 1e-12 and (cand is None or vv > cand[0]):
                            cand = (vv, a_, b_)
                        q[a_] += 1
                        q[b_] -= 1
                if cand is not None:
                    v, a_, b_ = cand
                    q[a_] -= 1
                    q[b_] += 1
                    improved = True
            cache[key] = q
        counts.append(cache[key])
    return _check_capacity(_assign_counts(stripe_of, np.stack(counts)), capacity)
