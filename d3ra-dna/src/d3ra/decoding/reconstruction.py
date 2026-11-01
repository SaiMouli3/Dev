"""Outer decoding, file reassembly, decompression and verification."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..coding.erasure_codes import make_code
from ..fragmentation.splitter import join_stripes
from ..ingestion.checksum import sha256
from ..ingestion.compressor import decompress
from ..ingestion.file_loader import FileRecord, unbundle
from ..storage.metadata import Manifest
from .verification import RecoveryStatus, classify


@dataclass
class RecoveryResult:
    status: RecoveryStatus
    files: list[FileRecord] = field(default_factory=list)
    file_status: dict = field(default_factory=dict)
    stripes_ok: int = 0
    stripes_total: int = 0


def reconstruct(manifest: Manifest, rows: dict[tuple[int, int], bytes]) -> RecoveryResult:
    """rows: {(stripe, index): payload bytes} of verified (inner-decoded) fragments."""
    S = len(manifest.stripes)
    stripes: list[np.ndarray | None] = []
    ok = 0
    for s, meta in enumerate(manifest.stripes):
        code = make_code(manifest.code, meta["k"], meta["n"])
        got = {i: np.frombuffer(p, dtype=np.uint8) for (ss, i), p in rows.items() if ss == s}
        if code.can_decode(got.keys()):
            try:
                stripes.append(code.decode(got))
                ok += 1
                continue
            except (ValueError, np.linalg.LinAlgError):
                pass
        stripes.append(None)
    if ok < S:
        return RecoveryResult(classify(None, manifest.blob_sha256, ok, S), stripes_ok=ok, stripes_total=S)
    blob = join_stripes(stripes, manifest.blob_length)
    if sha256(blob) != manifest.blob_sha256:
        return RecoveryResult(RecoveryStatus.INTEGRITY_MISMATCH, stripes_ok=ok, stripes_total=S)
    files = unbundle(decompress(blob, manifest.compression))
    fstat = {f.name: classify(f.data, m["sha256"], 1, 1) for f, m in zip(files, manifest.files)}
    status = RecoveryStatus.EXACT if all(v == RecoveryStatus.EXACT for v in fstat.values()) \
        else RecoveryStatus.INTEGRITY_MISMATCH
    return RecoveryResult(status, files, fstat, ok, S)
