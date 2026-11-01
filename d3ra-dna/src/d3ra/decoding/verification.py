"""Recovery outcome classification. A file counts as recovered only if its
SHA-256 matches -- a plausible-looking file is never treated as success."""
from __future__ import annotations

from enum import Enum

from ..ingestion.checksum import sha256


class RecoveryStatus(str, Enum):
    EXACT = "exact"
    PARTIAL = "partial"                # some stripes / files recovered, not all
    FAILED = "failed"                  # nothing decodable
    INTEGRITY_MISMATCH = "integrity_mismatch"  # decoded but hash mismatch


def classify(data: bytes | None, expected_sha: str, stripes_ok: int, stripes_total: int) -> RecoveryStatus:
    if data is not None:
        return RecoveryStatus.EXACT if sha256(data) == expected_sha else RecoveryStatus.INTEGRITY_MISMATCH
    if stripes_ok == 0:
        return RecoveryStatus.FAILED
    return RecoveryStatus.PARTIAL
