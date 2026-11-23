"""Shared cluster bootstrap used by every analysis script and by the manuscript verifier.

The resampling indices depend only on the number of clusters (synthesis replicates) and the number of
resamples, never on the order in which a script happens to compute its statistics. The same statistic on
the same retrievals therefore gets the same 95% interval in every table, in the text and in the verifier.
(Earlier revisions drew from one sequential random stream per script, so a quantity computed in two scripts
could be printed with two slightly different Monte-Carlo intervals.)
"""
from __future__ import annotations

import zlib

import numpy as np
import pandas as pd

B = 2000
_IDX: dict = {}


def index_sets(n_clusters: int, n_resamples: int = B) -> np.ndarray:
    """(n_resamples, n_clusters) array of cluster indices, fixed for a given (n_clusters, n_resamples)."""
    key = (int(n_clusters), int(n_resamples))
    if key not in _IDX:
        rng = np.random.default_rng(zlib.crc32(f"cluster-bootstrap|{key[0]}|{key[1]}".encode()))
        _IDX[key] = rng.integers(0, key[0], size=(key[1], key[0]))
    return _IDX[key]


def cluster_boot(df: pd.DataFrame, stat, by: str = "rep", n: int = B):
    """Percentile interval over clusters (sorted by `by`); returns (estimate, lo, hi)."""
    groups = [g for _, g in df.groupby(by, sort=True)]
    est = stat(df)
    vals = [stat(pd.concat([groups[i] for i in ix])) for ix in index_sets(len(groups), n)]
    return est, float(np.nanquantile(vals, 0.025)), float(np.nanquantile(vals, 0.975))


def matrix_boot(stat, mats, n: int = B):
    """Same resampling for (cluster x read-order) matrices: rows are clusters in sorted order."""
    vals = [stat(*[m[ix].ravel() for m in mats]) for ix in index_sets(mats[0].shape[0], n)]
    return np.nanpercentile(vals, [2.5, 97.5])
