"""Substitution channel with optional context/positional dependence."""
from __future__ import annotations

import numpy as np

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


def position_profile(L: int, slope: float = 0.0) -> np.ndarray:
    """Multiplicative per-position error profile (1 at 5' end, 1+slope at 3')."""
    return 1.0 + slope * np.linspace(0, 1, L)


def substitute(seq: np.ndarray, p: float | np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Replace each base with a uniformly random *different* base w.p. p."""
    seq = seq.copy()
    hit = rng.random(seq.size) < p
    if hit.any():
        idx = np.flatnonzero(hit)
        cur = np.searchsorted(BASES, seq[idx])
        seq[idx] = BASES[(cur + rng.integers(1, 4, idx.size)) % 4]
    return seq
