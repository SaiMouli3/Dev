"""Fragment abundance models.

Abundance of fragment i in container j is modelled as

    a_ij  ∝  c_ij * g_i * h_ij            (normalised within container j)

* g_i  : *sequence-intrinsic* bias shared by every copy of the same sequence
         (synthesis yield, GC-dependent PCR efficiency). Because replicas of a
         fragment are the same molecule species, g_i is shared across
         containers -- replicas do NOT see independent coverage.
* h_ij : container-specific variation (aliquoting, storage decay, local PCR).
* c_ij : number of copies (placement replication) of fragment i in j.

Both factors are log-normal or gamma (overdispersion); a gamma-Poisson
mixture gives negative-binomial read counts.
"""
from __future__ import annotations

import numpy as np

from ..encoding.constraints import gc_content


def pcr_efficiency(gc: np.ndarray, eff0: float = 0.95, kappa: float = 0.6, sd: float = 0.02,
                   rng: np.random.Generator | None = None) -> np.ndarray:
    """Per-sequence PCR efficiency decreasing with distance from 50% GC."""
    rng = rng or np.random.default_rng(0)
    eff = eff0 - kappa * np.abs(np.asarray(gc) - 0.5) + rng.normal(0, sd, np.shape(gc))
    return np.clip(eff, 0.3, 1.0)


def pcr_bias(gc: np.ndarray, cycles: int = 20, **kw) -> np.ndarray:
    """Relative abundance after ``cycles`` of PCR: (1+eff_i)^cycles, mean-normalised."""
    eff = pcr_efficiency(gc, **kw)
    a = np.exp(cycles * np.log1p(eff))
    return a / a.mean()


def intrinsic_bias(strands: list[str] | None, n: int, sigma: float, cycles: int,
                   rng: np.random.Generator) -> np.ndarray:
    """Shared, sequence-intrinsic multiplicative bias g_i (mean 1)."""
    g = rng.lognormal(-sigma ** 2 / 2, sigma, n) if sigma > 0 else np.ones(n)
    if strands is not None and cycles > 0:
        gc = np.array([gc_content(s) for s in strands])
        g = g * pcr_bias(gc, cycles, rng=rng)
    return g / g.mean()


def container_variation(n: int, m: int, shape: np.ndarray | float, rng: np.random.Generator) -> np.ndarray:
    """h_ij ~ Gamma(shape_j, 1/shape_j) (mean 1, CV = 1/sqrt(shape_j))."""
    shape = np.broadcast_to(np.asarray(shape, dtype=float), (m,))
    return rng.gamma(shape[None, :], 1.0 / shape[None, :], size=(n, m))


def abundance_matrix(g: np.ndarray, h: np.ndarray, placement: np.ndarray) -> np.ndarray:
    """Combine shared/container effects; zero where not placed; columns sum to 1
    over stored copies. ``placement``: (n, m) copy counts."""
    a = g[:, None] * h * placement
    s = a.sum(axis=0, keepdims=True)
    s[s == 0] = 1
    return a / s
