"""Correlation-aware (hierarchical) read-error model.

Each (strand, position) has a persistent error propensity
    p_sp ~ Beta(mean = mu_ctx(s,p), intra-class correlation = icc)
where mu_ctx is the context-dependent mean (5-mer + homopolymer run) and icc
captures how strongly errors repeat across reads of the same strand beyond
context. Every read of strand s errs at p independently w.p. p_sp; the error
type is drawn from the context's type mix. icc = 0 recovers the context-only
model; a constant mu recovers the i.i.d. model.
"""
from __future__ import annotations

import numpy as np

from ..realdata.mapping import parse_cigar

BASES = "ACGT"


def context_key(ref: str, p: int, k: int = 2) -> tuple:
    left = ref[max(0, p - k):p].rjust(k, "N")
    right = ref[p + 1:p + 1 + k].ljust(k, "N")
    # homopolymer run length containing p (capped)
    i, j = p, p
    while i > 0 and ref[i - 1] == ref[p]:
        i -= 1
    while j < len(ref) - 1 and ref[j + 1] == ref[p]:
        j += 1
    return (left + ref[p] + right, min(j - i + 1, 4))


def read_error_vector(ref: str, read: str):
    """Per reference position error type for one read: 0 ok, 1 sub, 2 del,
    3 insertion immediately before this position (first event wins)."""
    import edlib
    a = edlib.align(ref, read, mode="NW", task="path")
    v = np.zeros(len(ref), np.int8)
    if a["editDistance"] < 0:
        return None
    pos = 0
    for n, op in parse_cigar(a["cigar"]):
        if op == "=":
            pos += n
        elif op == "X":
            v[pos:pos + n] = np.where(v[pos:pos + n] == 0, 1, v[pos:pos + n])
            pos += n
        elif op == "I":
            v[pos:pos + n] = np.where(v[pos:pos + n] == 0, 2, v[pos:pos + n])
            pos += n
        elif op == "D":
            if pos < len(ref) and v[pos] == 0:
                v[pos] = 3
    return v


class HierarchicalErrorModel:
    def __init__(self, ctx_rate: dict, ctx_types: dict, global_rate: float, global_types: np.ndarray,
                 icc: float, min_count: int = 50):
        self.ctx_rate, self.ctx_types = ctx_rate, ctx_types
        self.global_rate, self.global_types = global_rate, global_types
        self.icc = icc
        self.min_count = min_count

    @classmethod
    def fit(cls, vectors_by_strand: dict, refs: dict, icc: float | None = None):
        """vectors_by_strand: {strand_id: [error vectors]}; refs: {strand_id: seq}."""
        num, den, types = {}, {}, {}
        g_err = g_n = 0
        g_types = np.zeros(4)
        for s, vecs in vectors_by_strand.items():
            ref = refs[s]
            V = np.stack(vecs)
            err = (V > 0)
            for p in range(len(ref)):
                key = context_key(ref, p)
                e = int(err[:, p].sum())
                num[key] = num.get(key, 0) + e
                den[key] = den.get(key, 0) + V.shape[0]
                t = types.setdefault(key, np.zeros(4))
                for k in (1, 2, 3):
                    t[k] += int((V[:, p] == k).sum())
            g_err += int(err.sum())
            g_n += err.size
            for k in (1, 2, 3):
                g_types[k] += int((V == k).sum())
        ctx_rate = {k: num[k] / den[k] for k in num if den[k] >= 50}
        ctx_types = {k: types[k] / max(types[k].sum(), 1) for k in types if den.get(k, 0) >= 50}
        m = cls(ctx_rate, ctx_types, g_err / g_n, g_types / g_types.sum(), 0.0)
        m.icc = icc if icc is not None else m.estimate_icc(vectors_by_strand, refs)
        return m

    def mean_rate(self, ref: str, p: int) -> float:
        return self.ctx_rate.get(context_key(ref, p), self.global_rate)

    def estimate_icc(self, vectors_by_strand: dict, refs: dict) -> float:
        """Method-of-moments intra-class correlation of error indicators across
        reads of the same (strand, position), after context normalisation:
        Var_between(rate) beyond binomial / (mu (1-mu))."""
        excess, base = 0.0, 0.0
        for s, vecs in vectors_by_strand.items():
            V = (np.stack(vecs) > 0).astype(float)
            n = V.shape[0]
            if n < 4:
                continue
            ref = refs[s]
            mu = np.array([self.mean_rate(ref, p) for p in range(len(ref))])
            obs = V.mean(0)
            # E[(obs-mu)^2] = mu(1-mu)/n + icc*mu(1-mu)(1-1/n)
            excess += ((obs - mu) ** 2 - mu * (1 - mu) / n).sum()
            base += (mu * (1 - mu) * (1 - 1 / n)).sum()
        return float(max(excess / base, 0.0))

    def strand_propensities(self, ref: str, rng: np.random.Generator, icc: float | None = None) -> np.ndarray:
        icc = self.icc if icc is None else icc
        mu = np.clip(np.array([self.mean_rate(ref, p) for p in range(len(ref))]), 1e-5, 0.9)
        if icc <= 1e-9:
            return mu
        conc = (1 - icc) / icc            # alpha + beta
        return rng.beta(mu * conc, (1 - mu) * conc)

    def reads(self, ref: str, n: int, rng: np.random.Generator, p_sp: np.ndarray | None = None,
              icc: float | None = None) -> list[str]:
        p_sp = self.strand_propensities(ref, rng, icc) if p_sp is None else p_sp
        out = []
        types = [self.ctx_types.get(context_key(ref, p), self.global_types) for p in range(len(ref))]
        for _ in range(n):
            hit = rng.random(len(ref)) < p_sp
            s = []
            for p, ch in enumerate(ref):
                if not hit[p]:
                    s.append(ch)
                    continue
                t = types[p]
                k = rng.choice(4, p=t / t.sum()) if t.sum() > 0 else 1
                if k == 1 or k == 0:
                    s.append(BASES[(BASES.index(ch) + int(rng.integers(1, 4))) % 4])
                elif k == 2:
                    pass
                else:
                    s.append(BASES[int(rng.integers(4))])
                    s.append(ch)
            out.append("".join(s))
        return out
