"""Track recovered fragments per stripe and per container."""
from __future__ import annotations

from collections import defaultdict

import numpy as np


class FragmentTracker:
    def __init__(self, k_per_stripe: np.ndarray, n_per_stripe: np.ndarray):
        self.k = np.asarray(k_per_stripe)
        self.n = np.asarray(n_per_stripe)
        self.rows: dict[tuple[int, int], bytes] = {}
        self.by_stripe = np.zeros(len(self.k), dtype=int)
        self.source: dict[tuple[int, int], int] = {}
        self.attempts: dict = defaultdict(int)

    def add(self, stripe: int, index: int, payload: bytes, container: int = -1) -> bool:
        key = (stripe, index)
        if key in self.rows:
            return False
        self.rows[key] = payload
        self.by_stripe[stripe] += 1
        self.source[key] = container
        return True

    def need(self) -> np.ndarray:
        return np.maximum(self.k - self.by_stripe, 0)

    def complete(self) -> bool:
        return bool((self.by_stripe >= self.k).all())

    def missing(self) -> list[tuple[int, int]]:
        return [(s, i) for s in range(len(self.k)) for i in range(self.n[s]) if (s, i) not in self.rows]

    def summary(self) -> dict:
        return {"recovered": len(self.rows), "stripes_complete": int((self.by_stripe >= self.k).sum()),
                "stripes": len(self.k), "min_margin": int((self.by_stripe - self.k).min())}
