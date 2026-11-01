"""Fragment dropout: molecules absent from a container.

* synthesis dropout  : sequence-intrinsic, shared by every replica (a failed
                       synthesis of sequence i is missing everywhere).
* container dropout  : independent per (fragment, container) with rate d_j.
* correlated dropout : probability rises with predicted bias (e.g. GC extremes
                       or low intrinsic abundance), inducing correlation with
                       coverage.
"""
from __future__ import annotations

import numpy as np


def sample_dropout(n: int, m: int, d_container: np.ndarray, d_synthesis: float,
                   rng: np.random.Generator, g: np.ndarray | None = None,
                   correlation: float = 0.0) -> np.ndarray:
    """Boolean (n, m) mask, True = molecule present."""
    d_container = np.broadcast_to(np.asarray(d_container, float), (m,))
    p_syn = np.full(n, d_synthesis)
    if g is not None and correlation > 0:
        # low-abundance sequences are more likely to drop out entirely
        w = np.exp(-correlation * np.log(np.maximum(g, 1e-9)))
        p_syn = np.clip(d_synthesis * w / w.mean(), 0, 0.5)
    syn_ok = rng.random(n) >= p_syn
    cont_ok = rng.random((n, m)) >= d_container[None, :]
    return syn_ok[:, None] & cont_ok
