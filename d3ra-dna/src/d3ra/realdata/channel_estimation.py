"""Estimate channel parameters from real reads with known references:
per-base substitution / insertion / deletion rates (overall, by position, in
homopolymers), per-read error-count dispersion, coverage distribution
(mean, CV, negative-binomial shape, log-normal sigma), dropout, GC effect,
and cross-container decomposition of coverage variance into a shared
(sequence-intrinsic) and a container-specific part.
"""
from __future__ import annotations

import numpy as np

from ..encoding.constraints import gc_content
from .mapping import error_events


def homopolymer_mask(seq: str) -> np.ndarray:
    m = np.zeros(len(seq), dtype=bool)
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            m[i] = m[i - 1] = True
    return m


def error_profile(pairs: list[tuple[str, str]], L: int | None = None) -> dict:
    """pairs: (reference, read) for mapped reads."""
    L = L or max(len(r) for r, _ in pairs)
    pos = {k: np.zeros(L) for k in ("sub", "ins", "del")}
    hp = {k: 0 for k in ("sub", "ins", "del")}
    hp_bases = nonhp_bases = 0
    per_read = []
    bases = 0
    for ref, read in pairs:
        ev = error_events(ref, read)
        if ev is None:
            continue
        m = homopolymer_mask(ref)
        hp_bases += m.sum()
        nonhp_bases += (~m).sum()
        bases += len(ref)
        per_read.append(len(ev))
        for p, kind in ev:
            if p < L:
                pos[kind][p] += 1
            if p < len(ref) and m[p]:
                hp[kind] += 1
    n = len(per_read)
    tot = {k: float(v.sum()) for k, v in pos.items()}
    out = {f"rate_{k}": tot[k] / max(bases, 1) for k in tot}
    out["rate_total"] = sum(tot.values()) / max(bases, 1)
    out["reads_aligned"] = n
    for k in hp:
        hp_rate = hp[k] / max(hp_bases, 1)
        nonhp_rate = (tot[k] - hp[k]) / max(nonhp_bases, 1)
        out[f"hp_ratio_{k}"] = hp_rate / max(nonhp_rate, 1e-12)
    pr = np.array(per_read)
    out["errors_per_read_mean"] = float(pr.mean()) if n else np.nan
    out["errors_per_read_var"] = float(pr.var()) if n else np.nan
    out["frac_error_free_reads"] = float((pr == 0).mean()) if n else np.nan
    out["position_profile"] = {k: (v / max(n, 1)).tolist() for k, v in pos.items()}
    return out


def coverage_stats(counts: np.ndarray, refs: list[str] | None = None) -> dict:
    c = np.asarray(counts, float)
    mu, var = c.mean(), c.var()
    nz = c[c > 0]
    out = {"strands": len(c), "mean_coverage": mu, "cv": np.sqrt(var) / max(mu, 1e-12),
           "nb_shape": mu ** 2 / max(var - mu, 1e-12) if var > mu else np.inf,
           "lognormal_sigma": float(np.log(nz).std()) if len(nz) > 2 else np.nan,
           "dropout_frac": float((c == 0).mean()),
           "p10_over_mean": float(np.quantile(c, 0.1) / max(mu, 1e-12)),
           "p90_over_mean": float(np.quantile(c, 0.9) / max(mu, 1e-12))}
    if refs is not None:
        gc = np.array([gc_content(r) for r in refs])
        ok = c > 0
        out["gc_coverage_spearman"] = float(_spearman(gc[ok], c[ok]))
    return out


def _spearman(a, b):
    ra = np.argsort(np.argsort(a))
    rb = np.argsort(np.argsort(b))
    return np.corrcoef(ra, rb)[0, 1]


def shared_variance_decomposition(count_matrix: np.ndarray) -> dict:
    """count_matrix: (strands, containers) read counts at comparable depth.
    Decompose log-coverage variance into a strand effect shared by all
    containers and a residual (container-specific) effect, after removing
    Poisson sampling noise (approx. var of log count ~ 1/mean)."""
    C = np.asarray(count_matrix, float)
    keep = (C > 0).all(1)
    Lc = np.log(C[keep] / C[keep].mean(0, keepdims=True))
    m = C[keep].mean(0)
    poisson_var = (1.0 / m).mean()            # delta-method sampling variance of log count
    strand_mean = Lc.mean(1, keepdims=True)
    resid = Lc - strand_mean
    J = C.shape[1]
    var_resid = resid.var(ddof=0) * J / (J - 1)
    var_shared = max(strand_mean.var() - var_resid / J, 0.0)
    var_specific = max(var_resid - poisson_var, 0.0)
    tot = var_shared + var_specific
    pair_corr = np.corrcoef(Lc.T)
    iu = np.triu_indices(J, 1)
    return {"var_shared": var_shared, "var_container_specific": var_specific, "poisson_var": poisson_var,
            "shared_fraction": var_shared / tot if tot > 0 else np.nan,
            "mean_pairwise_log_corr": float(pair_corr[iu].mean()), "strands_used": int(keep.sum())}
