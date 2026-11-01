"""Calibrate the surrogate's consensus-success curve phi(n) from the
sequence-level simulator: for n noisy reads of a random strand, run
two-sided BMA consensus and inner RS+CRC decoding, and record success.

Results are cached as JSON (data/processed/phi_calibration.json) keyed by the
error model and strand format, so experiments are reproducible and cheap.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from ..decoding.dna_decoder import StrandFormat
from ..fragmentation.addressing import AddressLayout
from ..retrieval.consensus import consensus
from .sequencer import ReadErrorModel, noisy_reads
from .surrogate import NMAX

CACHE = Path(__file__).resolve().parents[3] / "data" / "processed" / "phi_calibration.json"


def _key(err: ReadErrorModel, fmt_desc: dict, n_list, trials) -> str:
    d = {**err.__dict__, **fmt_desc, "n": list(n_list), "trials": trials}
    return hashlib.sha1(json.dumps(d, sort_keys=True).encode()).hexdigest()[:16]


def measure_phi(err: ReadErrorModel, fmt: StrandFormat | None = None, n_list=(1, 2, 3, 4, 5, 6, 8, 10, 13, 16, 20),
                trials: int = 60, seed: int = 0, use_cache: bool = True) -> dict:
    fmt = fmt or StrandFormat.build(AddressLayout(1, 2, 1), 32, 4, "scrambled-2bit")
    desc = {"payload": fmt.payload_bytes, "nsym": fmt.inner.nsym, "codec": fmt.codec.name}
    key = _key(err, desc, n_list, trials)
    cache = {}
    if use_cache and CACHE.exists():
        cache = json.loads(CACHE.read_text())
        if key in cache:
            return cache[key]
    rng = np.random.default_rng(seed)
    L = fmt.strand_nt
    succ = []
    for n in n_list:
        ok = 0
        for t in range(trials):
            payload = rng.integers(0, 256, fmt.payload_bytes, dtype=np.uint8).tobytes()
            strand = fmt.encode(0, int(t), int(n) % 256, payload)
            reads = noisy_reads(strand, n, err, rng)
            est = consensus(reads, L)
            dec = fmt.decode(est)
            ok += dec is not None and dec[1] == payload
        succ.append(ok / trials)
    res = {"n": list(n_list), "success": succ, "strand_nt": L, "error": err.__dict__, "format": desc}
    if use_cache:
        CACHE.parent.mkdir(parents=True, exist_ok=True)
        cache[key] = res
        CACHE.write_text(json.dumps(cache, indent=1))
    return res


def phi_from_measurement(res: dict, n_max: int = NMAX) -> np.ndarray:
    """Monotone interpolation of measured success to n = 0..n_max."""
    n = np.asarray(res["n"], float)
    s = np.maximum.accumulate(np.asarray(res["success"], float))
    grid = np.arange(n_max + 1, dtype=float)
    phi = np.interp(grid, np.concatenate([[0], n]), np.concatenate([[0], s]))
    # extrapolate beyond the last measured n with a geometric approach to 1
    last_n, last_s = n[-1], s[-1]
    if last_s < 1:
        rate = -np.log(max(1 - last_s, 1e-6)) / last_n
        tail = 1 - np.exp(-rate * grid)
        phi = np.where(grid > last_n, np.maximum(phi, tail), phi)
    phi[0] = 0
    return np.clip(np.maximum.accumulate(phi), 0, 1)


def calibrate_containers(containers, trials: int = 40, fmt: StrandFormat | None = None) -> np.ndarray:
    return np.stack([phi_from_measurement(measure_phi(ReadErrorModel.from_container(c), fmt, trials=trials))
                     for c in containers])
