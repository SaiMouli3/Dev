"""FASTQ / FASTA I/O (plain or gzip), streaming.

Quality strings for simulated reads are derived from the per-read error
probability (Phred = -10 log10 p), so simulated FASTQ files are consumable by
standard tools.
"""
from __future__ import annotations

import gzip
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


def _open(path, mode="rt"):
    path = str(path)
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


@dataclass
class FastqRecord:
    name: str
    seq: str
    qual: str

    def mean_phred(self) -> float:
        return sum(ord(c) - 33 for c in self.qual) / max(len(self.qual), 1)


def read_fastq(path, limit: int | None = None) -> Iterator[FastqRecord]:
    with _open(path) as fh:
        n = 0
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().strip()
            fh.readline()
            q = fh.readline().strip()
            yield FastqRecord(h[1:].strip(), s, q)
            n += 1
            if limit is not None and n >= limit:
                return


def write_fastq(path, records: Iterable[FastqRecord]) -> int:
    n = 0
    with _open(path, "wt") as fh:
        for r in records:
            fh.write(f"@{r.name}\n{r.seq}\n+\n{r.qual}\n")
            n += 1
    return n


def phred_string(length: int, p_error: float) -> str:
    q = int(min(41, max(2, round(-10 * math.log10(max(p_error, 1e-5))))))
    return chr(q + 33) * length


def read_fasta(path) -> tuple[list[str], list[str]]:
    names, seqs, cur = [], [], []
    with _open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
                names.append(line[1:])
            else:
                cur.append(line.upper())
    if cur:
        seqs.append("".join(cur))
    return names, seqs


def write_fasta(path, names, seqs) -> None:
    with _open(path, "wt") as fh:
        for n, s in zip(names, seqs):
            fh.write(f">{n}\n{s}\n")
