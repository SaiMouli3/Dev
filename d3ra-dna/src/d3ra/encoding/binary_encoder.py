"""Bit-level helpers and a fixed keystream whitener (scrambler).

Whitening with a fixed pseudo-random keystream (as in randomised DNA storage
encoders) makes long homopolymers and GC excursions unlikely without
reducing the 2 bits/nt rate; it is a *probabilistic* constraint, unlike the
rotating ternary code which guarantees no homopolymers.
"""
from __future__ import annotations

import numpy as np

_KEYSTREAM_CACHE: dict[tuple[int, int], np.ndarray] = {}


def bytes_to_bits(data: bytes) -> np.ndarray:
    return np.unpackbits(np.frombuffer(data, dtype=np.uint8))


def bits_to_bytes(bits: np.ndarray) -> bytes:
    return np.packbits(bits.astype(np.uint8)).tobytes()


def keystream(n: int, seed: int = 0xD3A) -> np.ndarray:
    key = (n, seed)
    if key not in _KEYSTREAM_CACHE:
        _KEYSTREAM_CACHE[key] = np.random.default_rng(seed).integers(0, 256, n, dtype=np.uint8)
    return _KEYSTREAM_CACHE[key]


def whiten(data: bytes, seed: int = 0xD3A) -> bytes:
    a = np.frombuffer(data, dtype=np.uint8)
    return (a ^ keystream(len(a), seed)).tobytes()


unwhiten = whiten  # XOR is an involution
