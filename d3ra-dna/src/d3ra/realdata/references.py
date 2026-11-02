"""Loaders for public DNA-storage datasets with ground-truth references.

* binned format (Zenodo 14296588, e.g. Srinivasavaradhan et al. clustered
  nanopore reads): reference line, a line of '*', the reads of that cluster,
  blank line(s) between clusters.
"""
from __future__ import annotations

from pathlib import Path


def read_binned(path, max_clusters: int | None = None) -> list[tuple[str, list[str]]]:
    out = []
    ref, reads, state = None, [], 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                if ref is not None and state == 2:
                    out.append((ref, reads))
                    ref, reads, state = None, [], 0
                    if max_clusters and len(out) >= max_clusters:
                        return out
                continue
            if state == 0:
                ref, state = line, 1
            elif state == 1 and set(line) == {"*"}:
                state = 2
            else:
                state = 2
                reads.append(line)
    if ref is not None:
        out.append((ref, reads))
    return out
