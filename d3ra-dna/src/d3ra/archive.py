"""On-disk archive: everything needed to retrieve without the original files.

archive/
  manifest.json        logical metadata (codec, code parameters, sizes, hashes)
  placement.npy        (n_fragments, M) copy counts
  containers.json      container descriptors (belief parameters)
  synthesis/container_XX.fasta   strands to synthesise per container
Fragment order is the library order (stripe-major), so the address of
fragment f is recoverable from the manifest alone.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .decoding.dna_decoder import StrandFormat
from .fragmentation.addressing import AddressLayout
from .io.fastq import write_fasta
from .storage.metadata import Manifest


def save_archive(lib, placement: np.ndarray, containers, path) -> Path:
    p = Path(path)
    (p / "synthesis").mkdir(parents=True, exist_ok=True)
    (p / "manifest.json").write_text(lib.manifest.to_json())
    np.save(p / "placement.npy", placement)
    (p / "containers.json").write_text(json.dumps([c.as_dict() for c in containers], indent=1))
    for j in range(placement.shape[1]):
        idx = np.flatnonzero(placement[:, j] > 0)
        write_fasta(p / "synthesis" / f"container_{j:02d}.fasta",
                    [f"frag{i}_s{lib.fragments[i].stripe}_p{lib.fragments[i].index}" for i in idx],
                    [lib.fragments[i].strand for i in idx])
    return p


def load_archive(path):
    from .storage.container import Container
    p = Path(path)
    man = Manifest.from_json((p / "manifest.json").read_text())
    placement = np.load(p / "placement.npy")
    containers = [Container(**c) for c in json.loads((p / "containers.json").read_text())]
    fmt = StrandFormat.build(AddressLayout(**man.layout), man.payload_bytes, man.inner_nsym, man.codec,
                             man.crc_bytes, man.addr_check)
    addresses = [(man.obj, s, i) for s, st in enumerate(man.stripes) for i in range(st["n"])]
    return man, fmt, placement, containers, addresses


def address_prefix(fmt: StrandFormat, addr) -> str:
    """Expected DNA prefix of a strand, computed from its address alone."""
    header = fmt._header(*addr)
    return fmt.codec.encode(header + bytes(fmt.message_bytes + fmt.inner.overhead() - len(header)))[: fmt.prefix_nt]
