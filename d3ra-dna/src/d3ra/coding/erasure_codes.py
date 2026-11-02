"""Erasure codes for the outer (cross-fragment) layer.

* ``MDSStripeCode`` wraps ``CauchyRS`` (any K of N).
* ``LTFountainCode``: Luby-transform code with robust-soliton degrees and a
  GF(2) Gaussian-elimination decoder (the family used by DNA Fountain). It is
  not MDS: decoding needs K(1+eps) symbols with a random eps.
* ``ReplicationCode``: r copies of each data fragment (N = rK).
"""
from __future__ import annotations

import math

import numpy as np

from .reed_solomon import CauchyRS


class MDSStripeCode:
    name = "rs-mds"

    def __init__(self, k: int, n: int):
        self.k, self.n = k, n
        self.rs = CauchyRS(k, n)

    def encode(self, data: np.ndarray) -> np.ndarray:
        return self.rs.encode(data)

    def can_decode(self, ids) -> bool:
        return len(set(ids)) >= self.k

    def decode(self, rows: dict[int, np.ndarray]) -> np.ndarray:
        return self.rs.decode(rows)


class ReplicationCode:
    name = "replication"

    def __init__(self, k: int, r: int):
        self.k, self.r, self.n = k, r, k * r

    def encode(self, data: np.ndarray) -> np.ndarray:
        return np.concatenate([data] * self.r, axis=0)

    def can_decode(self, ids) -> bool:
        return len({i % self.k for i in ids}) == self.k

    def decode(self, rows: dict[int, np.ndarray]) -> np.ndarray:
        out = {}
        for i, row in rows.items():
            out.setdefault(i % self.k, row)
        if len(out) < self.k:
            raise ValueError("missing data fragments")
        return np.stack([out[i] for i in range(self.k)])


def robust_soliton(k: int, c: float = 0.1, delta: float = 0.05) -> np.ndarray:
    R = c * math.log(k / delta) * math.sqrt(k)
    rho = np.array([1 / k] + [1 / (d * (d - 1)) for d in range(2, k + 1)])
    tau = np.zeros(k)
    pivot = max(1, min(k, int(round(k / R))))
    for d in range(1, pivot):
        tau[d - 1] = R / (d * k)
    tau[pivot - 1] = R * math.log(R / delta) / k if R / delta > 1 else 0
    mu = rho + tau
    return mu / mu.sum()


class LTFountainCode:
    name = "lt-fountain"

    def __init__(self, k: int, n: int, seed: int = 7):
        self.k, self.n, self.seed = k, n, seed
        dist = robust_soliton(k)
        rng = np.random.default_rng(seed)
        self.masks = np.zeros((n, k), dtype=bool)
        for i in range(n):
            d = int(rng.choice(np.arange(1, k + 1), p=dist))
            self.masks[i, rng.choice(k, size=d, replace=False)] = True

    def encode(self, data: np.ndarray) -> np.ndarray:
        out = np.zeros((self.n, data.shape[1]), dtype=np.uint8)
        for i in range(self.n):
            out[i] = np.bitwise_xor.reduce(data[self.masks[i]], axis=0)
        return out

    def _eliminate(self, ids, rows=None):
        A = self.masks[list(ids)].copy()
        Y = None if rows is None else np.stack([rows[i] for i in ids]).copy()
        r = 0
        piv_cols = []
        for c in range(self.k):
            nz = np.flatnonzero(A[r:, c])
            if nz.size == 0:
                continue
            p = r + nz[0]
            if p != r:
                A[[r, p]] = A[[p, r]]
                if Y is not None:
                    Y[[r, p]] = Y[[p, r]]
            others = np.flatnonzero(A[:, c])
            others = others[others != r]
            A[others] ^= A[r]
            if Y is not None:
                Y[others] ^= Y[r]
            piv_cols.append(c)
            r += 1
            if r == A.shape[0]:
                break
        return r, A, Y

    def can_decode(self, ids) -> bool:
        ids = sorted(set(ids))
        if len(ids) < self.k:
            return False
        return self._eliminate(ids)[0] == self.k

    def decode(self, rows: dict[int, np.ndarray]) -> np.ndarray:
        ids = sorted(rows)
        rank, A, Y = self._eliminate(ids, rows)
        if rank < self.k:
            raise ValueError("LT system not full rank")
        return Y[: self.k]


def make_code(kind: str, k: int, n: int):
    if kind in ("rs", "rs-mds", "mds"):
        return MDSStripeCode(k, n)
    if kind == "replication":
        if n % k:
            raise ValueError("replication needs n multiple of k")
        return ReplicationCode(k, n // k)
    if kind in ("lt", "lt-fountain", "fountain"):
        return LTFountainCode(k, n)
    raise ValueError(kind)
