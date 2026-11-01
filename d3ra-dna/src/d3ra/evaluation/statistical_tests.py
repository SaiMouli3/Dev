"""Statistics: bootstrap CIs (paired and unpaired), paired tests, effect
sizes and Holm-Bonferroni correction."""
from __future__ import annotations

import numpy as np
from scipy import stats


def bootstrap_ci(x, stat=np.mean, n_boot: int = 2000, alpha: float = 0.05, rng=None):
    rng = rng or np.random.default_rng(0)
    x = np.asarray(x)
    idx = rng.integers(0, len(x), (n_boot, len(x)))
    vals = np.array([stat(x[i]) for i in idx])
    return float(stat(x)), float(np.nanquantile(vals, alpha / 2)), float(np.nanquantile(vals, 1 - alpha / 2))


def paired_bootstrap_diff(a, b, stat=np.mean, n_boot: int = 2000, alpha: float = 0.05, rng=None):
    """CI for stat(a) - stat(b) resampling *pairs* (common random numbers)."""
    rng = rng or np.random.default_rng(0)
    a, b = np.asarray(a), np.asarray(b)
    idx = rng.integers(0, len(a), (n_boot, len(a)))
    vals = np.array([stat(a[i]) - stat(b[i]) for i in idx])
    return float(stat(a) - stat(b)), float(np.nanquantile(vals, alpha / 2)), float(np.nanquantile(vals, 1 - alpha / 2))


def unpaired_bootstrap_diff(a, b, stat=np.mean, n_boot: int = 2000, alpha: float = 0.05, rng=None):
    rng = rng or np.random.default_rng(0)
    a, b = np.asarray(a), np.asarray(b)
    va = np.array([stat(a[rng.integers(0, len(a), len(a))]) for _ in range(n_boot)])
    vb = np.array([stat(b[rng.integers(0, len(b), len(b))]) for _ in range(n_boot)])
    d = va - vb
    return float(stat(a) - stat(b)), float(np.nanquantile(d, alpha / 2)), float(np.nanquantile(d, 1 - alpha / 2))


def wilcoxon_paired(a, b) -> float:
    d = np.asarray(a) - np.asarray(b)
    d = d[np.isfinite(d) & (d != 0)]
    if len(d) < 10:
        return float("nan")
    return float(stats.wilcoxon(d).pvalue)


def mcnemar(success_a, success_b) -> float:
    """Exact McNemar test on paired binary outcomes."""
    a, b = np.asarray(success_a, bool), np.asarray(success_b, bool)
    n01, n10 = int((~a & b).sum()), int((a & ~b).sum())
    if n01 + n10 == 0:
        return 1.0
    return float(stats.binomtest(n01, n01 + n10, 0.5).pvalue)


def cohens_dz(a, b) -> float:
    d = np.asarray(a) - np.asarray(b)
    d = d[np.isfinite(d)]
    return float(d.mean() / (d.std(ddof=1) + 1e-12))


def cliffs_delta(a, b, max_n: int = 3000, rng=None) -> float:
    rng = rng or np.random.default_rng(0)
    a, b = np.asarray(a), np.asarray(b)
    if len(a) > max_n:
        a = rng.choice(a, max_n, replace=False)
    if len(b) > max_n:
        b = rng.choice(b, max_n, replace=False)
    gt = (a[:, None] > b[None, :]).mean()
    lt = (a[:, None] < b[None, :]).mean()
    return float(gt - lt)


def holm(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    run = 0.0
    for r, i in enumerate(order):
        run = max(run, (m - r) * p[i])
        adj[i] = min(1.0, run)
    return adj.tolist()


def wilson_ci(successes: int, n: int, alpha: float = 0.05):
    if n == 0:
        return (np.nan, np.nan, np.nan)
    z = stats.norm.ppf(1 - alpha / 2)
    p = successes / n
    den = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return p, max(0.0, c - h), min(1.0, c + h)
