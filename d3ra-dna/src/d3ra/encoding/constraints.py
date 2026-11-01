"""Biochemical sequence constraints: GC content, homopolymers, validation."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

BASES = "ACGT"


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def max_homopolymer(seq: str) -> int:
    best = run = 0
    prev = ""
    for ch in seq:
        run = run + 1 if ch == prev else 1
        prev = ch
        best = max(best, run)
    return best


def homopolymer_runs(seq: str, min_len: int = 4) -> list[tuple[int, int, str]]:
    """Return (start, length, base) for runs of length >= min_len."""
    runs, i = [], 0
    while i < len(seq):
        j = i
        while j < len(seq) and seq[j] == seq[i]:
            j += 1
        if j - i >= min_len:
            runs.append((i, j - i, seq[i]))
        i = j
    return runs


def windowed_gc(seq: str, window: int = 20) -> np.ndarray:
    gc = np.array([c in "GC" for c in seq], dtype=float)
    if len(gc) < window:
        return np.array([gc.mean()]) if len(gc) else np.array([0.0])
    k = np.ones(window) / window
    return np.convolve(gc, k, mode="valid")


@dataclass
class ConstraintSpec:
    gc_min: float = 0.40
    gc_max: float = 0.60
    max_run: int = 3
    window: int | None = None  # optional local-GC window
    window_gc_min: float = 0.25
    window_gc_max: float = 0.75


@dataclass
class ConstraintReport:
    gc: float
    max_run: int
    ok_gc: bool
    ok_run: bool
    ok_window: bool

    @property
    def ok(self) -> bool:
        return self.ok_gc and self.ok_run and self.ok_window


def validate(seq: str, spec: ConstraintSpec | None = None) -> ConstraintReport:
    spec = spec or ConstraintSpec()
    if set(seq) - set(BASES):
        raise ValueError("sequence contains non-ACGT symbols")
    gc = gc_content(seq)
    mr = max_homopolymer(seq)
    ok_w = True
    if spec.window:
        w = windowed_gc(seq, spec.window)
        ok_w = bool((w.min() >= spec.window_gc_min) and (w.max() <= spec.window_gc_max))
    return ConstraintReport(gc, mr, spec.gc_min <= gc <= spec.gc_max, mr <= spec.max_run, ok_w)
