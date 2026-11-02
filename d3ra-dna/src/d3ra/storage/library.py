"""Build a DNA library from files: bundle -> compress -> stripes -> outer
code -> strands. The library is the bridge between the digital pipeline and
the distributed storage simulator."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .. import ENCODING_VERSION
from ..coding.erasure_codes import make_code
from ..coding.redundancy import RedundancyPlan
from ..decoding.dna_decoder import StrandFormat
from ..fragmentation.addressing import AddressLayout
from ..fragmentation.splitter import Fragment, split_into_stripes
from ..ingestion.checksum import crc32, sha256
from ..ingestion.compressor import compress
from ..ingestion.file_loader import FileRecord, bundle
from .metadata import Manifest


@dataclass
class DNALibrary:
    manifest: Manifest
    fmt: StrandFormat
    fragments: list[Fragment] = field(default_factory=list)

    @property
    def strands(self) -> list[str]:
        return [f.strand for f in self.fragments]

    @property
    def stripe_of(self) -> np.ndarray:
        return np.array([f.stripe for f in self.fragments])

    def k_per_stripe(self) -> np.ndarray:
        return np.array([s["k"] for s in self.manifest.stripes])

    def index_of(self) -> dict:
        return {f.uid: i for i, f in enumerate(self.fragments)}

    def total_nt(self) -> int:
        return sum(len(s) for s in self.strands)

    def stats(self) -> dict:
        info_bytes = self.manifest.blob_length
        nt = self.total_nt()
        return {
            "fragments": len(self.fragments),
            "stripes": len(self.manifest.stripes),
            "strand_nt": self.fmt.strand_nt,
            "total_nt": nt,
            "info_bytes": info_bytes,
            "bits_per_nt": 8 * info_bytes / max(nt, 1),
            "code_rate": sum(s["k"] for s in self.manifest.stripes) / sum(s["n"] for s in self.manifest.stripes),
        }


def build_library(files: list[FileRecord], payload_bytes: int = 32, k: int = 32, parity: int = 8,
                  inner_nsym: int = 4, codec: str = "scrambled-2bit", compression: str = "auto",
                  code: str = "rs-mds", plan: RedundancyPlan | None = None, obj: int = 0,
                  crc_bytes: int = 2) -> DNALibrary:
    raw = bundle(files)
    blob, method = compress(raw, compression)
    data_stripes = split_into_stripes(blob, payload_bytes, k)
    S = len(data_stripes)
    if plan is None:
        plan = RedundancyPlan.uniform(S, k, parity, code=code)
    assert plan.n_stripes == S, "redundancy plan must match stripe count"
    max_n = max(plan.n)
    layout = AddressLayout.compact(1, S, max_n)
    fmt = StrandFormat.build(layout, payload_bytes, inner_nsym, codec, crc_bytes=crc_bytes)
    frags: list[Fragment] = []
    for s, data in enumerate(data_stripes):
        n_s = plan.n[s]
        cw = make_code(plan.code, k, n_s).encode(data)
        for i in range(n_s):
            payload = cw[i].tobytes()
            frags.append(Fragment(obj, s, i, n_s, k, payload, fmt.encode(obj, s, i, payload),
                                  ENCODING_VERSION, plan.code, crc32(payload)))
    man = Manifest(ENCODING_VERSION, codec, method, len(blob), sha256(blob), payload_bytes, inner_nsym,
                   fmt.inner.crc_bytes, fmt.addr_check, layout.__dict__, plan.code,
                   [{"k": k, "n": n, "crc32": crc32(data_stripes[i].tobytes())} for i, n in enumerate(plan.n)],
                   [f.metadata() for f in files], obj)
    return DNALibrary(man, fmt, frags)
