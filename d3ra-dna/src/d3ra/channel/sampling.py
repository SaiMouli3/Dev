"""Read sampling models: which fragment each read comes from."""
from __future__ import annotations

import numpy as np


def sample_counts(abundance: np.ndarray, reads: int, model: str, rng: np.random.Generator,
                  dispersion: float = 10.0) -> np.ndarray:
    """Read counts per fragment for one container.

    poisson     : n_i ~ Poisson(R a_i) independently
    multinomial : (n_i) ~ Multinomial(R, a)      (exact fixed-total)
    negbin      : n_i ~ Poisson(R a_i G_i), G_i ~ Gamma(k, 1/k)  (overdispersed)
    """
    a = np.asarray(abundance, dtype=float)
    if reads <= 0 or a.sum() == 0:
        return np.zeros_like(a, dtype=np.int64)
    a = a / a.sum()
    if model == "poisson":
        return rng.poisson(reads * a)
    if model == "multinomial":
        return rng.multinomial(reads, a)
    if model in ("negbin", "negative_binomial"):
        g = rng.gamma(dispersion, 1.0 / dispersion, a.shape)
        return rng.poisson(reads * a * g)
    raise ValueError(model)


def apply_duplicates(counts: np.ndarray, dup_rate: float, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """PCR/optical duplicates: a fraction of reads are copies of another read
    of the same molecule and add no independent information.
    Returns (total_reads, independent_reads)."""
    if dup_rate <= 0:
        return counts, counts
    indep = rng.binomial(counts, 1 - dup_rate)
    indep = np.where((counts > 0) & (indep == 0), 1, indep)
    return counts, indep


def apply_missing_reads(counts: np.ndarray, loss: float, rng: np.random.Generator) -> np.ndarray:
    """Reads lost to filtering / failed demultiplexing (thinning)."""
    return rng.binomial(counts, 1 - loss) if loss > 0 else counts
