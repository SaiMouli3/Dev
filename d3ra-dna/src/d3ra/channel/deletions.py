"""Deletion channel."""
from __future__ import annotations

import numpy as np


def delete(seq: np.ndarray, p: float | np.ndarray, rng: np.random.Generator) -> np.ndarray:
    keep = rng.random(seq.size) >= p
    return seq[keep]


def burst_delete(seq: np.ndarray, p_burst: float, mean_len: float, rng: np.random.Generator) -> np.ndarray:
    """Correlated (burst) deletions: geometric-length bursts starting w.p. p_burst."""
    if p_burst <= 0:
        return seq
    keep = np.ones(seq.size, dtype=bool)
    for s in np.flatnonzero(rng.random(seq.size) < p_burst):
        keep[s:s + rng.geometric(1 / mean_len)] = False
    return seq[keep]
