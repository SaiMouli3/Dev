"""Optional lossless compression. The chosen method is recorded in metadata."""
from __future__ import annotations

import bz2
import lzma
import zlib

_METHODS = {
    "none": (lambda b: b, lambda b: b),
    "zlib": (lambda b: zlib.compress(b, 9), zlib.decompress),
    "bz2": (lambda b: bz2.compress(b, 9), bz2.decompress),
    "lzma": (lambda b: lzma.compress(b, preset=9), lzma.decompress),
}
METHOD_IDS = {name: i for i, name in enumerate(_METHODS)}
METHOD_NAMES = {i: n for n, i in METHOD_IDS.items()}


def compress(data: bytes, method: str = "auto") -> tuple[bytes, str]:
    """Compress ``data``. ``auto`` picks the smallest output (incl. 'none')."""
    if method == "auto":
        best = min(((_METHODS[m][0](data), m) for m in _METHODS), key=lambda t: len(t[0]))
        return best
    if method not in _METHODS:
        raise ValueError(f"unknown compression method {method!r}")
    return _METHODS[method][0](data), method


def decompress(data: bytes, method: str) -> bytes:
    return _METHODS[method][1](data)
