"""File fragmentation into stripes of equal-sized chunks + logical metadata."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Fragment:
    """Logical fragment. ``payload`` is the outer-codeword row; ``strand`` is
    the physical DNA sequence (address + payload + inner ECC, encoded)."""

    obj: int
    stripe: int
    index: int           # position within stripe codeword (0..N_s-1)
    n_in_stripe: int     # total fragment count of the stripe (logical)
    k_in_stripe: int
    payload: bytes = field(repr=False)
    strand: str = field(default="", repr=False)
    encoding_version: int = 1
    code: str = "rs-mds"
    crc: int = 0

    @property
    def is_parity(self) -> bool:
        return self.index >= self.k_in_stripe

    @property
    def uid(self) -> tuple[int, int, int]:
        return (self.obj, self.stripe, self.index)


def split_into_stripes(blob: bytes, payload_bytes: int, k: int) -> list[np.ndarray]:
    """Return list of (k, payload_bytes) uint8 arrays; last stripe zero-padded."""
    chunk = payload_bytes * k
    n_stripes = max(1, -(-len(blob) // chunk))
    padded = blob + bytes(n_stripes * chunk - len(blob))
    arr = np.frombuffer(padded, dtype=np.uint8).reshape(n_stripes, k, payload_bytes)
    return [arr[s].copy() for s in range(n_stripes)]


def join_stripes(stripes: list[np.ndarray], length: int) -> bytes:
    return b"".join(s.tobytes() for s in stripes)[:length]
