"""Read clustering and consensus (trace reconstruction).

Clustering is *inventory-aware*: the container inventory provides a codebook
of expected address prefixes; each read is assigned to the nearest codebook
entry (Hamming distance on the prefix, also trying +-1 shifts to absorb an
early indel) if within a tolerance, otherwise discarded.

Consensus uses two-sided Bitwise Majority Alignment (Batu et al., 2004) with
a lookahead that classifies disagreements as substitution / insertion /
deletion; the forward pass supplies the first half and a reversed pass the
second half of the estimate.
"""
from __future__ import annotations

from collections import Counter, defaultdict

import numpy as np


def _prefix_array(seqs: list[str], length: int) -> np.ndarray:
    buf = np.full((len(seqs), length), ord("N"), dtype=np.uint8)
    for r, s in enumerate(seqs):
        b = s[:length].encode()
        buf[r, : len(b)] = np.frombuffer(b, dtype=np.uint8)
    return buf


class InventoryClusterer:
    def __init__(self, keys: list, prefixes: list[str], prefix_len: int, max_dist: int | None = None):
        self.keys = list(keys)
        self.prefix_len = prefix_len
        self.book = _prefix_array(prefixes, prefix_len)
        self.max_dist = max_dist if max_dist is not None else max(2, prefix_len // 5)

    def assign(self, reads: list[str]) -> list:
        if not reads or not self.keys:
            return [None] * len(reads)
        out = []
        P = self.prefix_len
        best_d = np.full(len(reads), 10**6)
        best_k = np.zeros(len(reads), dtype=int)
        for shift in (0, 1, -1):
            if shift >= 0:
                arr = _prefix_array([r[shift:] for r in reads], P)
            else:
                arr = _prefix_array(["N" + r for r in reads], P)
            # (reads x book) Hamming distances, chunked for memory
            for s in range(0, len(reads), 512):
                d = (arr[s:s + 512, None, :] != self.book[None, :, :]).sum(axis=2) + abs(shift)
                k = d.argmin(axis=1)
                dv = d[np.arange(d.shape[0]), k]
                upd = dv < best_d[s:s + 512]
                best_d[s:s + 512][upd] = dv[upd]
                best_k[s:s + 512][upd] = k[upd]
        for d, k in zip(best_d, best_k):
            out.append(self.keys[k] if d <= self.max_dist else None)
        return out


def cluster_reads(reads: list[str], clusterer: InventoryClusterer) -> dict:
    groups: dict = defaultdict(list)
    for key, r in zip(clusterer.assign(reads), reads):
        if key is not None:
            groups[key].append(r)
    return groups


def _bma(reads: list[str], L: int, w: int = 3) -> str:
    ptr = [0] * len(reads)
    out = []
    for _ in range(L):
        chars = [r[p] if p < len(r) else None for r, p in zip(reads, ptr)]
        cnt = Counter(c for c in chars if c is not None)
        if not cnt:
            out.append("A")
            continue
        cand = cnt.most_common(1)[0][0]
        out.append(cand)
        agree = [k for k, c in enumerate(chars) if c == cand]
        fut_cnt = [Counter() for _ in range(w)]
        for k in agree:
            r, p = reads[k], ptr[k]
            for t in range(w):
                if p + 1 + t < len(r):
                    fut_cnt[t][r[p + 1 + t]] += 1
        fut = "".join(fc.most_common(1)[0][0] if fc else "N" for fc in fut_cnt)
        for k, c in enumerate(chars):
            if c is None:
                continue
            if c == cand:
                ptr[k] += 1
                continue
            r, p = reads[k], ptr[k]

            def score(seg: str) -> int:
                return sum(a == b for a, b in zip(seg, fut))
            s_sub = score(r[p + 1:p + 1 + w])
            s_del = score(r[p:p + w])            # read skipped cand: current char is the next symbol
            s_ins = score(r[p + 2:p + 2 + w]) + (1 if p + 1 < len(r) and r[p + 1] == cand else -1)
            best = max((s_sub, 0), (s_ins, 1), (s_del, 2))[1]
            ptr[k] += (1, 2, 0)[best]
    return "".join(out)


def consensus(reads: list[str], L: int, w: int = 3) -> str:
    """Two-sided BMA consensus of length L."""
    if not reads:
        return ""
    if len(reads) == 1:
        return reads[0][:L].ljust(L, "A")
    fwd = _bma(reads, L, w)
    bwd = _bma([r[::-1] for r in reads], L, w)[::-1]
    h = L // 2
    return fwd[:h] + bwd[h:]


def majority_vote(reads: list[str], L: int) -> str:
    """Naive position-wise majority (baseline that ignores indels)."""
    out = []
    for t in range(L):
        c = Counter(r[t] for r in reads if t < len(r))
        out.append(c.most_common(1)[0][0] if c else "A")
    return "".join(out)
