"""Redundancy plans: how many parity fragments each stripe receives and how
fragments are replicated across containers.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class RedundancyPlan:
    """Per-stripe (K_s, N_s) plus per-fragment replication counts."""

    k: list[int]
    n: list[int]
    replication: list[np.ndarray] = field(default_factory=list)  # per stripe, (N_s,) copies
    code: str = "rs-mds"

    def __post_init__(self):
        if not self.replication:
            self.replication = [np.ones(n, dtype=int) for n in self.n]

    @classmethod
    def uniform(cls, n_stripes: int, k: int, parity: int, copies: int = 1, code: str = "rs-mds"):
        return cls([k] * n_stripes, [k + parity] * n_stripes,
                   [np.full(k + parity, copies, dtype=int) for _ in range(n_stripes)], code)

    @property
    def n_stripes(self) -> int:
        return len(self.k)

    def total_fragments(self) -> int:
        return int(sum(r.sum() for r in self.replication))

    def unique_fragments(self) -> int:
        return int(sum(self.n))

    def storage_overhead(self) -> float:
        """Stored molecules (distinct sequences x copies) per information fragment."""
        return self.total_fragments() / sum(self.k)

    def code_rate(self) -> float:
        return sum(self.k) / sum(self.n)
