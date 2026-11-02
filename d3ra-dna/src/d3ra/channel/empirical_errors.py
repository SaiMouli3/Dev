"""Empirical error transplant.

Edit scripts (substitution / insertion / deletion events with their relative
positions and inserted/substituted bases) are extracted from real reads
aligned to their designed references. Simulated reads of *our* strands are
corrupted by replaying a randomly chosen real edit script, rescaled to the
strand length. This preserves real per-type rates, positional profiles,
within-read error correlation and the heavy tail of per-read error counts.
"""
from __future__ import annotations

import numpy as np

from ..realdata.mapping import parse_cigar

import edlib


def edit_script(ref: str, read: str):
    a = edlib.align(ref, read, mode="SHW", task="path")
    if a["editDistance"] < 0:
        return None
    ops, i, j = [], 0, 0
    L = len(ref)
    for n, op in parse_cigar(a["cigar"]):
        if op == "=":
            i += n
            j += n
        elif op == "X":
            for t in range(n):
                ops.append((i / L, "S", read[j]))
                i += 1
                j += 1
        elif op == "I":
            for t in range(n):
                ops.append((i / L, "D", ""))
                i += 1
        elif op == "D":
            for t in range(n):
                ops.append((i / L, "I", read[j]))
                j += 1
    return ops


class ErrorTransplant:
    name = "transplant"

    def __init__(self, pairs: list[tuple[str, str]], max_scripts: int = 50000):
        self.scripts = []
        for ref, read in pairs[:max_scripts]:
            s = edit_script(ref, read)
            if s is not None:
                self.scripts.append(s)
        if not self.scripts:
            raise ValueError("no alignable pairs")
        n_ev = np.array([len(s) for s in self.scripts])
        self.mean_errors = float(n_ev.mean())
        self.frac_clean = float((n_ev == 0).mean())

    def corrupt(self, seq: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        s = list(seq.tobytes().decode())
        L = len(s)
        script = self.scripts[rng.integers(len(self.scripts))]
        # apply from the 3' end so earlier positions stay valid
        for rel, op, base in sorted(script, key=lambda x: -x[0]):
            p = min(int(rel * L), L - 1 if op != "I" else L)
            if op == "S" and p < len(s):
                s[p] = base if base and base != s[p] else "ACGT"[("ACGT".index(s[p]) + 1) % 4]
            elif op == "D" and p < len(s):
                del s[p]
            elif op == "I":
                s.insert(p, base or "A")
        return np.frombuffer("".join(s).encode(), dtype=np.uint8).copy()
