"""Container-level failures with shared failure domains (common-mode causes).

A failure state is a boolean vector over containers (True = failed). Each
container fails independently w.p. f_j, and each failure domain d (a set of
containers sharing e.g. a freezer, rack, buffer batch or site) fails as a
whole w.p. q_d.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np


@dataclass
class FailureModel:
    f: np.ndarray                                   # (M,) independent failure probs
    domains: list[list[int]] = field(default_factory=list)
    q: list[float] = field(default_factory=list)   # domain failure probs

    @property
    def m(self) -> int:
        return len(self.f)

    def sample(self, n: int, rng: np.random.Generator) -> np.ndarray:
        failed = rng.random((n, self.m)) < np.asarray(self.f)[None, :]
        for dom, q in zip(self.domains, self.q):
            hit = rng.random(n) < q
            if hit.any():
                failed[np.ix_(hit, dom)] = True
        return failed

    def enumerate_states(self, max_states: int = 4096):
        """Exact enumeration of (state, prob) when feasible; otherwise None.

        Enumerates domain events and independent failures jointly; prunes
        states with negligible probability (< 1e-12)."""
        m, nd = self.m, len(self.domains)
        if 2 ** (m + nd) > max_states:
            return None
        out: dict[tuple, float] = {}
        f = np.asarray(self.f)
        for dstate in itertools.product([0, 1], repeat=nd):
            pd = np.prod([q if s else 1 - q for s, q in zip(dstate, self.q)]) if nd else 1.0
            forced = np.zeros(m, bool)
            for s, dom in zip(dstate, self.domains):
                if s:
                    forced[dom] = True
            for istate in itertools.product([0, 1], repeat=m):
                ist = np.array(istate, bool)
                pi = np.prod(np.where(ist, f, 1 - f))
                st = tuple(forced | ist)
                out[st] = out.get(st, 0.0) + pd * pi
        return [(np.array(s), p) for s, p in out.items() if p > 1e-12]

    def marginal(self) -> np.ndarray:
        """P(container j failed)."""
        surv = 1 - np.asarray(self.f, float)
        for dom, q in zip(self.domains, self.q):
            surv[dom] *= 1 - q
        return 1 - surv

    def domain_of(self) -> np.ndarray:
        """Primary domain label per container (-1 if none)."""
        lab = np.full(self.m, -1)
        for d, dom in enumerate(self.domains):
            lab[dom] = np.where(lab[dom] < 0, d, lab[dom])
        return lab
