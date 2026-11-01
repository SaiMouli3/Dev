"""Placement optimisation: exact ILP (small instances) and heuristics
(greedy, local search, simulated annealing, genetic algorithm).

Objective (fast, decomposable surrogate of the full predictor):

    J(X) = sum_s  [ cost of the K_s cheapest fragments of stripe s ]  / C_ref
         + delta * P_loss(X) / eps

* first term: expected retrieval cost when the reader collects the cheapest
  K_s fragments of every stripe (unit costs e_j, see placement.unit_retrieval_cost);
* P_loss(X): probability that a failure state leaves some stripe with < K_s
  surviving fragments (exact over enumerated states, or sample-average);
* eps: tolerated failure probability (1 - target).

The full analytical/MC predictor is used afterwards to *validate* placements.
"""
from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix

from ..storage.placement import cost_tiered, stripe_failure_exposure


@dataclass
class PlacementProblem:
    stripe_of: np.ndarray
    k: np.ndarray
    unit_cost: np.ndarray
    states: np.ndarray            # (F, M) failure states
    weights: np.ndarray           # (F,)
    capacity: np.ndarray | None = None
    delta: float = 1.0
    eps: float = 1e-3
    domain: np.ndarray | None = None

    @property
    def m(self):
        return len(self.unit_cost)

    def c_ref(self) -> float:
        return float(self.k.sum() * self.unit_cost.min())

    def counts(self, assign: np.ndarray) -> np.ndarray:
        S = self.stripe_of.max() + 1
        C = np.zeros((S, self.m), dtype=int)
        np.add.at(C, (self.stripe_of, assign), 1)
        return C

    def objective(self, assign: np.ndarray) -> float:
        C = self.counts(assign)
        return self.objective_counts(C)

    def objective_counts(self, C: np.ndarray) -> float:
        # retrieval-cost term
        order = np.argsort(self.unit_cost)
        cost = 0.0
        for s in range(C.shape[0]):
            need = self.k[s]
            for j in order:
                take = min(need, C[s, j])
                cost += take * self.unit_cost[j]
                need -= take
                if need == 0:
                    break
        # failure term: alive counts per stripe per state
        alive = C @ (1 - self.states.T)               # (S, F)
        ok = (alive >= self.k[:, None]).all(0)
        ploss = 1 - (ok * self.weights).sum() / self.weights.sum()
        pen = 0.0
        if self.capacity is not None:
            over = np.maximum(C.sum(0) - self.capacity, 0).sum()
            pen = 1e3 * over
        return cost / self.c_ref() + self.delta * ploss / self.eps + pen


def to_matrix(assign: np.ndarray, m: int) -> np.ndarray:
    X = np.zeros((len(assign), m))
    X[np.arange(len(assign)), assign] = 1
    return X


def from_matrix(X: np.ndarray) -> np.ndarray:
    return X.argmax(1)


# ----------------------------------------------------------------------------- heuristics

def greedy(prob: PlacementProblem, cap: int | None = None) -> np.ndarray:
    return from_matrix(cost_tiered(prob.stripe_of, prob.k, prob.unit_cost, prob.domain, cap, prob.capacity))


def local_search(prob: PlacementProblem, assign: np.ndarray, iters: int = 5000, rng=None) -> np.ndarray:
    """First-improvement moves on per-stripe container counts (works on counts)."""
    rng = rng or np.random.default_rng(0)
    C = prob.counts(assign)
    best = prob.objective_counts(C)
    S, M = C.shape
    for _ in range(iters):
        s = rng.integers(S)
        a, b = rng.choice(M, 2, replace=False)
        if C[s, a] == 0:
            continue
        C[s, a] -= 1
        C[s, b] += 1
        v = prob.objective_counts(C)
        if v <= best:
            best = v
        else:
            C[s, a] += 1
            C[s, b] -= 1
    return counts_to_assign(prob.stripe_of, C)


def simulated_annealing(prob: PlacementProblem, assign: np.ndarray, iters: int = 20000, t0: float = 0.05,
                        rng=None) -> np.ndarray:
    rng = rng or np.random.default_rng(0)
    C = prob.counts(assign)
    cur = prob.objective_counts(C)
    best, bestC = cur, C.copy()
    S, M = C.shape
    for it in range(iters):
        T = t0 * (1 - it / iters) + 1e-6
        s = rng.integers(S)
        a, b = rng.choice(M, 2, replace=False)
        if C[s, a] == 0:
            continue
        C[s, a] -= 1
        C[s, b] += 1
        v = prob.objective_counts(C)
        if v <= cur or rng.random() < np.exp(-(v - cur) / T):
            cur = v
            if v < best:
                best, bestC = v, C.copy()
        else:
            C[s, a] += 1
            C[s, b] -= 1
    return counts_to_assign(prob.stripe_of, bestC)


def genetic(prob: PlacementProblem, pop: int = 30, gens: int = 60, rng=None,
            seed_assign: np.ndarray | None = None) -> np.ndarray:
    """GA over per-stripe count rows (crossover = stripe-wise mixing)."""
    rng = rng or np.random.default_rng(0)
    S = prob.stripe_of.max() + 1
    M = prob.m
    n_s = np.bincount(prob.stripe_of)

    def rand_counts():
        return np.stack([rng.multinomial(n_s[s], np.ones(M) / M) for s in range(S)])

    P = [rand_counts() for _ in range(pop)]
    if seed_assign is not None:
        P[0] = prob.counts(seed_assign)
    fit = np.array([prob.objective_counts(c) for c in P])
    for _ in range(gens):
        children = []
        for _c in range(pop):
            i, j = rng.choice(pop, 2, replace=False)
            a = P[i] if fit[i] < fit[j] else P[j]
            i, j = rng.choice(pop, 2, replace=False)
            b = P[i] if fit[i] < fit[j] else P[j]
            mask = rng.random(S) < 0.5
            ch = np.where(mask[:, None], a, b).copy()
            for _m in range(2):  # mutation: move one fragment
                s = rng.integers(S)
                x, y = rng.choice(M, 2, replace=False)
                if ch[s, x] > 0:
                    ch[s, x] -= 1
                    ch[s, y] += 1
            children.append(ch)
        cf = np.array([prob.objective_counts(c) for c in children])
        allP = P + children
        allf = np.concatenate([fit, cf])
        keep = np.argsort(allf)[:pop]
        P = [allP[k] for k in keep]
        fit = allf[keep]
    return counts_to_assign(prob.stripe_of, P[0])


def counts_to_assign(stripe_of: np.ndarray, C: np.ndarray) -> np.ndarray:
    assign = np.zeros(len(stripe_of), dtype=int)
    for s in range(C.shape[0]):
        idx = np.flatnonzero(stripe_of == s)
        assign[idx] = np.repeat(np.arange(C.shape[1]), C[s])[: len(idx)]
    return assign


# ----------------------------------------------------------------------------- exact

def ilp(prob: PlacementProblem, time_limit: float = 60.0) -> tuple[np.ndarray, dict]:
    """Exact MILP on per-stripe counts for small instances.

    Variables: y_sj (integer count of stripe s in container j), u_sj (fragments of
    stripe s *read* from j, u <= y, sum_j u_sj = K_s), z_f (state f is a loss state).
    Loss linearisation: for each stripe s and state f, alive_sf = sum_j y_sj (1-F_fj);
    alive_sf >= K_s - K_s z_f. Objective: sum u_sj e_j / C_ref + delta/eps * sum_f w_f z_f.
    """
    S = prob.stripe_of.max() + 1
    M = prob.m
    F = prob.states.shape[0]
    n_s = np.bincount(prob.stripe_of)
    w = prob.weights / prob.weights.sum()
    ny, nu, nz = S * M, S * M, F
    nv = ny + nu + nz
    c = np.zeros(nv)
    c[ny:ny + nu] = np.tile(prob.unit_cost, S) / prob.c_ref()
    c[ny + nu:] = prob.delta / prob.eps * w
    rows, lo, hi = [], [], []
    A = lil_matrix((S + S + S * M + S * F + (M if prob.capacity is not None else 0), nv))
    r = 0
    for s in range(S):                      # sum_j y_sj = N_s
        A[r, s * M:(s + 1) * M] = 1
        lo.append(n_s[s]); hi.append(n_s[s]); r += 1
    for s in range(S):                      # sum_j u_sj = K_s
        A[r, ny + s * M: ny + (s + 1) * M] = 1
        lo.append(prob.k[s]); hi.append(prob.k[s]); r += 1
    for s in range(S):                      # u <= y
        for j in range(M):
            A[r, ny + s * M + j] = 1
            A[r, s * M + j] = -1
            lo.append(-np.inf); hi.append(0); r += 1
    for s in range(S):                      # alive + K z >= K
        for f in range(F):
            for j in range(M):
                if prob.states[f, j] == 0:
                    A[r, s * M + j] = 1
            A[r, ny + nu + f] = prob.k[s]
            lo.append(prob.k[s]); hi.append(np.inf); r += 1
    if prob.capacity is not None:
        for j in range(M):
            for s in range(S):
                A[r, s * M + j] = 1
            lo.append(-np.inf); hi.append(prob.capacity[j]); r += 1
    ub = np.concatenate([np.repeat(n_s, M), np.repeat(prob.k, M), np.ones(nz)])
    t = time.time()
    res = milp(c, constraints=LinearConstraint(A.tocsr(), lo, hi), integrality=np.ones(nv),
               bounds=Bounds(np.zeros(nv), ub), options={"time_limit": time_limit, "disp": False})
    info = {"status": res.status, "message": res.message, "time": time.time() - t,
            "objective": float(res.fun) if res.x is not None else np.nan,
            "mip_gap": getattr(res, "mip_gap", np.nan)}
    if res.x is None:
        return greedy(prob), info
    Y = np.rint(res.x[:ny]).astype(int).reshape(S, M)
    return counts_to_assign(prob.stripe_of, Y), info


def compare_methods(prob: PlacementProblem, rng=None, ilp_limit: float = 30.0) -> list[dict]:
    rng = rng or np.random.default_rng(0)
    out = []
    M = prob.m
    rr = np.array([(i + s) % M for i, s in zip(np.arange(len(prob.stripe_of)), prob.stripe_of)])
    candidates = {
        "random": lambda: rng.integers(0, M, len(prob.stripe_of)),
        "round_robin": lambda: rr,
        "greedy_tiered": lambda: greedy(prob),
        "local_search": lambda: local_search(prob, greedy(prob), rng=rng),
        "simulated_annealing": lambda: simulated_annealing(prob, greedy(prob), rng=rng),
        "genetic": lambda: genetic(prob, rng=rng, seed_assign=greedy(prob)),
    }
    for name, fn in candidates.items():
        t = time.time()
        a = fn()
        out.append({"method": name, "objective": prob.objective(a), "time": time.time() - t})
    t = time.time()
    a, info = ilp(prob, ilp_limit)
    out.append({"method": "ilp", "objective": prob.objective(a), "time": time.time() - t,
                "ilp_status": info["status"]})
    best = min(r["objective"] for r in out if r["method"] == "ilp")
    for r in out:
        r["gap_vs_ilp"] = (r["objective"] - best) / max(abs(best), 1e-12)
    return out
