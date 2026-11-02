"""Binary <-> DNA codecs.

* ``TwoBitCodec``       : 2 bits/nt, no constraint guarantees.
* ``ScrambledTwoBitCodec``: 2 bits/nt after keystream whitening (probabilistic constraints).
* ``RotatingTernaryCodec``: Goldman-style rotation code; base-3 digits mapped
  relative to the previous base, so homopolymers never occur
  (log2(3)~1.58 bits/nt raw; 40 bits -> 26 trits = 1.538 bits/nt here).

All codecs are deterministic and length-preserving for a fixed byte length,
which is what the fragment addressing layer requires.
"""
from __future__ import annotations

import numpy as np

from .binary_encoder import unwhiten, whiten

BASES = "ACGT"
_B2I = {b: i for i, b in enumerate(BASES)}
_LUT = np.frombuffer(b"ACGT", dtype=np.uint8)
_REV = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _REV[_b] = _i


class Codec:
    name = "abstract"
    bits_per_nt = 0.0

    def encode(self, data: bytes) -> str:  # pragma: no cover - interface
        raise NotImplementedError

    def decode(self, seq: str, n_bytes: int) -> bytes:  # pragma: no cover
        raise NotImplementedError

    def length_for(self, n_bytes: int) -> int:
        return len(self.encode(bytes(n_bytes)))


class TwoBitCodec(Codec):
    name = "2bit"
    bits_per_nt = 2.0

    def encode(self, data: bytes) -> str:
        a = np.frombuffer(data, dtype=np.uint8)
        q = np.stack([(a >> 6) & 3, (a >> 4) & 3, (a >> 2) & 3, a & 3], axis=1).ravel()
        return _LUT[q].tobytes().decode()

    def decode(self, seq: str, n_bytes: int) -> bytes:
        s = seq[: 4 * n_bytes].ljust(4 * n_bytes, "A")
        q = _REV[np.frombuffer(s.encode(), dtype=np.uint8)]
        q = np.where(q == 255, 0, q).astype(np.uint8).reshape(-1, 4)
        a = (q[:, 0] << 6) | (q[:, 1] << 4) | (q[:, 2] << 2) | q[:, 3]
        return a.astype(np.uint8).tobytes()


class ScrambledTwoBitCodec(TwoBitCodec):
    name = "scrambled-2bit"

    def __init__(self, seed: int = 0xD3A):
        self.seed = seed

    def encode(self, data: bytes) -> str:
        return super().encode(whiten(data, self.seed))

    def decode(self, seq: str, n_bytes: int) -> bytes:
        return unwhiten(super().decode(seq, n_bytes), self.seed)


class RotatingTernaryCodec(Codec):
    """Goldman et al. (2013)-style rotating code (no homopolymers by construction)."""

    name = "rotating-ternary"
    BLOCK_BYTES = 5
    BLOCK_TRITS = 26  # 3**26 > 2**40
    bits_per_nt = 40 / 26

    def _trits(self, data: bytes) -> list[int]:
        pad = (-len(data)) % self.BLOCK_BYTES
        data = data + bytes(pad)
        out: list[int] = []
        for k in range(0, len(data), self.BLOCK_BYTES):
            v = int.from_bytes(data[k:k + self.BLOCK_BYTES], "big")
            digs = []
            for _ in range(self.BLOCK_TRITS):
                v, r = divmod(v, 3)
                digs.append(r)
            out.extend(reversed(digs))
        return out

    def encode(self, data: bytes) -> str:
        prev = "A"  # virtual previous base (fixed so decoding is stateless)
        out = []
        for t in self._trits(data):
            choices = [b for b in BASES if b != prev]
            prev = choices[t]
            out.append(prev)
        return "".join(out)

    def decode(self, seq: str, n_bytes: int) -> bytes:
        nblocks = -(-n_bytes // self.BLOCK_BYTES)
        need = nblocks * self.BLOCK_TRITS
        seq = seq[:need].ljust(need, "A")
        prev = "A"
        trits = []
        for b in seq:
            choices = [c for c in BASES if c != prev]
            trits.append(choices.index(b) if b in choices else 0)  # homopolymer => error
            prev = b
        out = bytearray()
        for k in range(nblocks):
            v = 0
            for t in trits[k * self.BLOCK_TRITS:(k + 1) * self.BLOCK_TRITS]:
                v = v * 3 + t
            v &= (1 << (8 * self.BLOCK_BYTES)) - 1  # out-of-range (corrupted) values wrap
            out += v.to_bytes(self.BLOCK_BYTES, "big")
        return bytes(out[:n_bytes])


CODECS = {c.name: c for c in (TwoBitCodec, ScrambledTwoBitCodec, RotatingTernaryCodec)}


def get_codec(name: str) -> Codec:
    return CODECS[name]()
