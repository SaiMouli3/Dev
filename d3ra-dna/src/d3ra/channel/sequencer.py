"""Sequence-level sequencing simulator: molecules -> noisy reads.

Combines sampling, dropout, PCR/intrinsic bias, substitution / insertion /
deletion (homopolymer- and position-dependent), per-read quality variation,
burst errors, duplicates, missing reads and container-level failure.
All randomness flows from an explicit ``numpy.random.Generator``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .deletions import burst_delete, delete
from .insertions import homopolymer_factor, insert
from .sampling import apply_duplicates, apply_missing_reads, sample_counts
from .substitutions import position_profile, substitute

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


@dataclass
class ReadErrorModel:
    e_sub: float = 0.005
    e_ins: float = 0.001
    e_del: float = 0.003
    homopolymer_boost: float = 1.5
    quality_sigma: float = 0.3
    position_slope: float = 0.5
    burst_rate: float = 0.0
    burst_len: float = 3.0

    @classmethod
    def from_container(cls, c) -> "ReadErrorModel":
        return cls(c.e_sub, c.e_ins, c.e_del, c.homopolymer_boost, c.quality_sigma,
                   burst_rate=c.burst_rate)

    def corrupt(self, seq: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        q = rng.lognormal(-self.quality_sigma ** 2 / 2, self.quality_sigma) if self.quality_sigma else 1.0
        prof = position_profile(seq.size, self.position_slope) * q
        hp = homopolymer_factor(seq, self.homopolymer_boost)
        s = substitute(seq, np.minimum(self.e_sub * prof, 0.5), rng)
        hp = homopolymer_factor(s, self.homopolymer_boost)
        s2 = delete(s, np.minimum(self.e_del * prof * hp, 0.5), rng)
        s3 = insert(s2, np.minimum(self.e_ins * position_profile(s2.size, self.position_slope) * q, 0.5), rng)
        if self.burst_rate:
            s3 = burst_delete(s3, self.burst_rate, self.burst_len, rng)
        return s3


def to_array(seq: str) -> np.ndarray:
    return np.frombuffer(seq.encode(), dtype=np.uint8).copy()


def to_str(a: np.ndarray) -> str:
    return a.tobytes().decode()


def noisy_reads(strand: str, n: int, err: ReadErrorModel, rng: np.random.Generator) -> list[str]:
    a = to_array(strand)
    return [to_str(err.corrupt(a, rng)) for _ in range(n)]


def junk_read(length: int, rng: np.random.Generator) -> str:
    return to_str(BASES[rng.integers(0, 4, length)])


def sequence_container(strands: list[str], abundance: np.ndarray, present: np.ndarray, reads: int,
                       err: ReadErrorModel, rng: np.random.Generator, model: str = "multinomial",
                       yield_: float = 1.0, dup_rate: float = 0.0, loss: float = 0.0,
                       failed: bool = False, dispersion: float = 10.0) -> list[tuple[int, str]]:
    """Return list of (true_fragment_index or -1 for junk, read string)."""
    if failed or reads <= 0:
        L = len(strands[0]) if strands else 100
        return [(-1, junk_read(L, rng)) for _ in range(reads if not failed else reads)]
    usable = rng.binomial(reads, yield_)
    a = np.asarray(abundance, float) * present
    counts = sample_counts(a, usable, model, rng, dispersion) if a.sum() > 0 else np.zeros(len(strands), int)
    counts = apply_missing_reads(counts, loss, rng)
    _, indep = apply_duplicates(counts, dup_rate, rng)
    out: list[tuple[int, str]] = []
    for i in np.flatnonzero(indep):
        rs = noisy_reads(strands[i], int(indep[i]), err, rng)
        # duplicates: exact copies of an independent read
        extra = int(counts[i] - indep[i])
        rs += [rs[rng.integers(len(rs))] for _ in range(extra)]
        out += [(int(i), r) for r in rs]
    L = len(strands[0])
    out += [(-1, junk_read(L, rng)) for _ in range(reads - usable)]
    order = rng.permutation(len(out))
    return [out[k] for k in order]
