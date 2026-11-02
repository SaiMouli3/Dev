"""Map raw reads to designed reference strands.

* adapter trimming (Illumina TruSeq read-1 adapter read-through)
* k-mer voting index over the references (unique k-mers only) to propose a
  reference, verified by edit distance (edlib, prefix-global alignment) with a
  relative threshold. Unassigned reads are 'junk' (adapter dimers, PhiX,
  chimeras, heavily corrupted reads).
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import edlib
import numpy as np

ADAPTER_R1 = "AGATCGGAAGAGC"


def trim_adapter(seq: str, adapter: str = ADAPTER_R1, min_overlap: int = 6) -> str:
    i = seq.find(adapter)
    if i >= 0:
        return seq[:i]
    # partial adapter at the 3' end
    for k in range(len(adapter) - 1, min_overlap - 1, -1):
        if seq.endswith(adapter[:k]):
            return seq[:-k]
    return seq


@dataclass
class MapResult:
    ref: int          # -1 if unmapped
    dist: int         # edit distance to the reference (or -1)


class KmerMapper:
    def __init__(self, refs: list[str], k: int = 14, stride: int = 7, max_rel_dist: float = 0.25):
        self.refs, self.k, self.stride, self.max_rel = refs, k, stride, max_rel_dist
        counts = Counter()
        for r in refs:
            for p in range(0, len(r) - k + 1):
                counts[r[p:p + k]] += 1
        self.index = {}
        for i, r in enumerate(refs):
            for p in range(0, len(r) - k + 1):
                km = r[p:p + k]
                if counts[km] == 1:
                    self.index[km] = i

    def propose(self, read: str) -> int:
        votes = Counter()
        k = self.k
        for p in range(0, max(1, len(read) - k + 1), self.stride):
            i = self.index.get(read[p:p + k])
            if i is not None:
                votes[i] += 1
        if not votes:
            # dense fallback (indels can break all strided k-mers)
            for p in range(0, max(1, len(read) - k + 1)):
                i = self.index.get(read[p:p + k])
                if i is not None:
                    votes[i] += 1
        return votes.most_common(1)[0][0] if votes else -1

    def map(self, read: str) -> MapResult:
        i = self.propose(read)
        if i < 0:
            return MapResult(-1, -1)
        ref = self.refs[i]
        d = edlib.align(ref, read, mode="SHW", task="distance")["editDistance"]
        if d < 0 or d > self.max_rel * len(ref):
            return MapResult(-1, d)
        return MapResult(i, d)


def parse_cigar(cigar: str):
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
        else:
            yield int(num), ch
            num = ""


def error_events(ref: str, read: str):
    """Edit operations of the read relative to the reference, as a list of
    (ref_position, kind) with kind in {'sub','ins','del'} (SHW alignment:
    reference globally, read prefix)."""
    a = edlib.align(ref, read, mode="SHW", task="path")
    if a["editDistance"] < 0:
        return None
    ev = []
    pos = 0
    for n, op in parse_cigar(a["cigar"]):
        if op == "=":
            pos += n
        elif op == "X":
            ev += [(pos + t, "sub") for t in range(n)]
            pos += n
        elif op == "I":       # base in reference (query) missing from read -> deletion
            ev += [(pos + t, "del") for t in range(n)]
            pos += n
        elif op == "D":       # extra base in read -> insertion before pos
            ev += [(pos, "ins")] * n
    return ev
