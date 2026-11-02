"""Insertion channel (homopolymer-elevated)."""
from __future__ import annotations

import numpy as np

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


def homopolymer_factor(seq: np.ndarray, boost: float = 1.0) -> np.ndarray:
    """Per-position multiplier boost**(run_length-1) for indels inside runs."""
    if boost == 1.0 or seq.size == 0:
        return np.ones(seq.size)
    run = np.ones(seq.size)
    same = seq[1:] == seq[:-1]
    for i in np.flatnonzero(same) + 1:
        run[i] = run[i - 1] + 1
    return boost ** (run - 1)


def insert(seq: np.ndarray, p: float | np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """After each position insert a random base w.p. p (half are duplications)."""
    hit = rng.random(seq.size) < p
    if not hit.any():
        return seq
    idx = np.flatnonzero(hit)
    ins = np.where(rng.random(idx.size) < 0.5, seq[idx], BASES[rng.integers(0, 4, idx.size)])
    return np.insert(seq, idx + 1, ins)
