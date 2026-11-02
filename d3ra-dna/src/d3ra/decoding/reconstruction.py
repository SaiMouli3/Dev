"""Outer decoding, file reassembly, decompression and verification."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..coding.erasure_codes import make_code
from ..fragmentation.splitter import join_stripes
from itertools import combinations

from ..ingestion.checksum import crc32, sha256
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
    repaired: list = field(default_factory=list)       # (stripe, excluded fragment indices)
    unverified_stripes: list = field(default_factory=list)


def decode_stripe_verified(code, got: dict, crc: int | None, max_exclude: int = 2):
    """Decode a stripe; if a per-stripe CRC is known and fails, search for the
    corrupted fragment(s) by excluding up to `max_exclude` of the available
    fragments (needs |got| >= K + #excluded). Returns (data or None, excluded)."""
    ids = sorted(got)
    def attempt(use):
        try:
            return code.decode({i: got[i] for i in use})
        except (ValueError, np.linalg.LinAlgError):
            return None
    d = attempt(ids)
    if d is not None and (crc is None or crc32(d.tobytes()) == crc):
        return d, ()
    if crc is None:
        return None, ()
    for r in range(1, max_exclude + 1):
        if len(ids) - r < code.k:
            break
        for ex in combinations(ids, r):
            use = [i for i in ids if i not in ex]
            d = attempt(use)
            if d is not None and crc32(d.tobytes()) == crc:
                return d, ex
    return None, ()


def reconstruct(manifest: Manifest, rows: dict[tuple[int, int], bytes]) -> RecoveryResult:
    """rows: {(stripe, index): payload bytes} of verified (inner-decoded) fragments."""
    S = len(manifest.stripes)
    stripes: list[np.ndarray | None] = []
    ok = 0
    repaired, bad = [], []
    for s, meta in enumerate(manifest.stripes):
        code = make_code(manifest.code, meta["k"], meta["n"])
        got = {i: np.frombuffer(p, dtype=np.uint8) for (ss, i), p in rows.items() if ss == s}
        if code.can_decode(got.keys()):
            d, ex = decode_stripe_verified(code, got, meta.get("crc32"))
            if d is not None:
                stripes.append(d)
                ok += 1
                if ex:
                    repaired.append((s, ex))
                continue
            bad.append(s)
        stripes.append(None)
    if ok < S:
        r = RecoveryResult(classify(None, manifest.blob_sha256, ok, S), stripes_ok=ok, stripes_total=S)
        r.unverified_stripes, r.repaired = bad, repaired
        return r
    blob = join_stripes(stripes, manifest.blob_length)
    if sha256(blob) != manifest.blob_sha256:
        return RecoveryResult(RecoveryStatus.INTEGRITY_MISMATCH, stripes_ok=ok, stripes_total=S)
    files = unbundle(decompress(blob, manifest.compression))
    fstat = {f.name: classify(f.data, m["sha256"], 1, 1) for f, m in zip(files, manifest.files)}
    status = RecoveryStatus.EXACT if all(v == RecoveryStatus.EXACT for v in fstat.values()) \
        else RecoveryStatus.INTEGRITY_MISMATCH
    return RecoveryResult(status, files, fstat, ok, S, repaired)
