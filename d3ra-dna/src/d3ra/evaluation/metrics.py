"""Metrics for recovery reliability, sequencing / storage efficiency,
retrieval performance and computation."""
from __future__ import annotations

import numpy as np

from ..retrieval.fixed_policy import budget_for_target
from .statistical_tests import bootstrap_ci, wilson_ci


def quantile_ci(samples, target, n_boot=500, rng=None):
    rng = rng or np.random.default_rng(0)
    s = np.asarray(samples)
    vals = [budget_for_target(s[rng.integers(0, len(s), len(s))], target) for _ in range(n_boot)]
    vals = np.asarray(vals)
    fin = vals[np.isfinite(vals)]
    lo = float(np.quantile(fin, 0.025)) if len(fin) > 0.975 * len(vals) else np.inf
    hi = float(np.quantile(vals, 0.975)) if np.isfinite(vals).mean() >= 0.975 else np.inf
    return budget_for_target(s, target), lo, hi


def adaptive_summary(cost: np.ndarray, success: np.ndarray, targets=(0.99, 0.999), budget_ref: float | None = None,
                     rounds=None, containers_used=None, extra=None) -> dict:
    """cost: uncapped stopping cost per trial (inf if never succeeded within the run cap)."""
    c = np.where(success, cost, np.inf)
    out = {"trials": len(c), "p_success_uncapped": float(np.isfinite(c).mean())}
    for t in targets:
        cap, lo, hi = quantile_ci(c, t)
        out[f"cap@{t}"] = cap
        out[f"cap@{t}_ci"] = (lo, hi)
        spend = np.minimum(np.where(np.isfinite(c), c, cap), cap) if np.isfinite(cap) else np.full(len(c), np.inf)
        m, l, h = bootstrap_ci(spend) if np.isfinite(cap) else (np.inf, np.inf, np.inf)
        out[f"mean_spend@{t}"] = m
        out[f"mean_spend@{t}_ci"] = (l, h)
    fin = c[np.isfinite(c)]
    out["mean_cost_success"] = float(fin.mean()) if len(fin) else np.inf
    if budget_ref is not None:
        k = int((c <= budget_ref).sum())
        p, lo, hi = wilson_ci(k, len(c))
        out["p_exact@ref"] = p
        out["p_exact@ref_ci"] = (lo, hi)
    if rounds is not None:
        out["mean_rounds"] = float(np.mean(rounds))
    if containers_used is not None:
        out["mean_containers_accessed"] = float(np.mean(containers_used))
    if extra is not None:
        out["mean_extra_after_probe"] = float(np.mean(extra))
    return out


def static_summary(samples: np.ndarray, targets=(0.99, 0.999), budget_ref: float | None = None) -> dict:
    out = {"trials": len(samples), "p_success_uncapped": float(np.isfinite(samples).mean())}
    for t in targets:
        cap, lo, hi = quantile_ci(samples, t)
        out[f"cap@{t}"] = cap
        out[f"cap@{t}_ci"] = (lo, hi)
        out[f"mean_spend@{t}"] = cap          # a static policy commits its whole budget
        out[f"mean_spend@{t}_ci"] = (lo, hi)
    fin = samples[np.isfinite(samples)]
    out["mean_cost_success"] = float(fin.mean()) if len(fin) else np.inf
    if budget_ref is not None:
        k = int((samples <= budget_ref).sum())
        p, lo, hi = wilson_ci(k, len(samples))
        out["p_exact@ref"] = p
        out["p_exact@ref_ci"] = (lo, hi)
    return out


def storage_metrics(layout, info_bytes_per_fragment: int, strand_nt: int) -> dict:
    n_mol = float(layout.placement.sum())
    k = float(layout.k.sum())
    return {
        "encoded_fragments": int(layout.n),
        "stored_molecule_species": n_mol,
        "redundancy_overhead": n_mol / k - 1,
        "replication_factor": n_mol / layout.n,
        "encoded_bases": n_mol * strand_nt,
        "info_density_bits_per_nt": 8 * info_bytes_per_fragment * k / (n_mol * strand_nt),
    }
