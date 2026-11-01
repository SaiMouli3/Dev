"""Integrity primitives: SHA-256 for files, CRC-32 for fragments."""
from __future__ import annotations

import hashlib
import zlib


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_bytes(data: bytes) -> bytes:
    return hashlib.sha256(data).digest()


def crc32(data: bytes) -> int:
    return zlib.crc32(data) & 0xFFFFFFFF


def crc32_bytes(data: bytes) -> bytes:
    return crc32(data).to_bytes(4, "big")


def verify_sha256(data: bytes, expected_hex: str) -> bool:
    return sha256(data) == expected_hex
