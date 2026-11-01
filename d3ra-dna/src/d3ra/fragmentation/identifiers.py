"""Stable identifiers for files/objects and fragments."""
from __future__ import annotations

import hashlib


def object_id(sha256_hex: str, bits: int = 8) -> int:
    """Short object id derived from content hash (collisions resolved by manifest)."""
    return int(sha256_hex[:8], 16) % (1 << bits)


def fragment_uid(obj: int, stripe: int, index: int) -> str:
    return f"{obj:x}:{stripe}:{index}"


def manifest_id(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()[:16]
