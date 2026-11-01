"""Reed-Solomon codes.

* ``GF256``: vectorised GF(2^8) arithmetic (primitive polynomial 0x11d).
* ``CauchyRS``: systematic MDS erasure code over GF(256) used as the *outer*
  code across the N fragments of a stripe (any K of N recover the stripe).
* ``InnerRS``: per-strand error-correcting code (``reedsolo``) with a CRC
  guard; strands that fail decoding or the CRC are declared erasures.
"""
from __future__ import annotations

import numpy as np
from reedsolo import ReedSolomonError, RSCodec

from ..ingestion.checksum import crc32


class GF256:
    PRIM = 0x11D
    exp = np.zeros(512, dtype=np.int32)
    log = np.zeros(256, dtype=np.int32)
    _x = 1
    for _i in range(255):
        exp[_i] = _x
        log[_x] = _i
        _x <<= 1
        if _x & 0x100:
            _x ^= PRIM
    exp[255:510] = exp[:255]
    MUL = np.zeros((256, 256), dtype=np.uint8)
    _a = np.arange(1, 256)
    MUL[1:, 1:] = exp[(log[_a][:, None] + log[_a][None, :]) % 255]
    INV = np.zeros(256, dtype=np.uint8)
    INV[1:] = exp[(255 - log[_a]) % 255]

    @classmethod
    def matmul(cls, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        """(r x k) @ (k x n) over GF(256)."""
        A = np.asarray(A, dtype=np.uint8)
        B = np.asarray(B, dtype=np.uint8)
        out = np.zeros((A.shape[0], B.shape[1]), dtype=np.uint8)
        for k in range(A.shape[1]):
            out ^= cls.MUL[A[:, k]][:, B[k]]
        return out

    @classmethod
    def inv_matrix(cls, M: np.ndarray) -> np.ndarray:
        M = np.array(M, dtype=np.uint8)
        n = M.shape[0]
        A = np.concatenate([M, np.eye(n, dtype=np.uint8)], axis=1)
        for c in range(n):
            piv = c + int(np.flatnonzero(A[c:, c])[0]) if A[c:, c].any() else -1
            if piv < 0:
                raise np.linalg.LinAlgError("singular matrix over GF(256)")
            if piv != c:
                A[[c, piv]] = A[[piv, c]]
            A[c] = cls.MUL[cls.INV[A[c, c]]][A[c]]
            nz = np.flatnonzero(A[:, c])
            nz = nz[nz != c]
            if nz.size:
                A[nz] ^= cls.MUL[A[nz, c][:, None], A[c][None, :]]
        return A[:, n:]


class CauchyRS:
    """Systematic (N, K) MDS erasure code. Requires N <= 256."""

    def __init__(self, k: int, n: int):
        if not (0 < k <= n <= 256):
            raise ValueError("need 0 < k <= n <= 256")
        self.k, self.n = k, n
        p = n - k
        x = np.arange(k, k + p)
        y = np.arange(k)
        C = GF256.INV[(x[:, None] ^ y[None, :]).astype(np.uint8)]
        self.G = np.concatenate([np.eye(k, dtype=np.uint8), C.astype(np.uint8)], axis=0)

    @property
    def rate(self) -> float:
        return self.k / self.n

    def encode(self, data: np.ndarray) -> np.ndarray:
        """data: (K, B) uint8 -> (N, B) codeword rows (first K are systematic)."""
        data = np.asarray(data, dtype=np.uint8)
        assert data.shape[0] == self.k
        if self.n == self.k:
            return data.copy()
        return np.concatenate([data, GF256.matmul(self.G[self.k:], data)], axis=0)

    def decode(self, rows: dict[int, np.ndarray]) -> np.ndarray:
        """rows: {codeword index -> (B,) row}. Needs >= K distinct indices."""
        idx = sorted(rows)[: self.k]
        if len(idx) < self.k:
            raise ValueError(f"need {self.k} fragments, have {len(idx)}")
        if idx == list(range(self.k)):
            return np.stack([rows[i] for i in idx])
        sub = self.G[idx]
        Y = np.stack([np.asarray(rows[i], dtype=np.uint8) for i in idx])
        return GF256.matmul(GF256.inv_matrix(sub), Y)


class InnerRS:
    """Per-strand RS(nsym) over bytes + 2-byte CRC guard on the message."""

    def __init__(self, nsym: int = 4, crc_bytes: int = 2):
        self.nsym = nsym
        self.crc_bytes = crc_bytes
        self.codec = RSCodec(nsym) if nsym > 0 else None

    def overhead(self) -> int:
        return self.nsym + self.crc_bytes

    def _crc(self, msg: bytes) -> bytes:
        return crc32(msg).to_bytes(4, "big")[: self.crc_bytes]

    def encode(self, msg: bytes) -> bytes:
        body = msg + self._crc(msg)
        return bytes(self.codec.encode(body)) if self.codec else body

    def decode(self, word: bytes) -> bytes | None:
        """Return message bytes, or None if uncorrectable / CRC mismatch."""
        try:
            body = bytes(self.codec.decode(word)[0]) if self.codec else word
        except ReedSolomonError:
            return None
        msg, crc = body[: -self.crc_bytes or None], body[-self.crc_bytes:] if self.crc_bytes else b""
        if self.crc_bytes and self._crc(msg) != crc:
            return None
        return msg
