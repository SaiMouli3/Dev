"""Command-line interface: encode/store/retrieve files through the simulator."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .ingestion.file_loader import FileRecord
from .pipeline import store_and_retrieve
from .storage.container import make_heterogeneous


def main(argv=None):
    ap = argparse.ArgumentParser(prog="d3ra", description="D3RA-DNA simulated store-and-retrieve")
    ap.add_argument("files", nargs="+")
    ap.add_argument("--containers", type=int, default=6)
    ap.add_argument("--placement", default="tiered", choices=["tiered", "domain", "round_robin"])
    ap.add_argument("--controller", default="pivotal",
                    choices=["uniform", "coverage", "pivotal", "thompson"])
    ap.add_argument("--parity", type=int, default=8)
    ap.add_argument("--k", type=int, default=24)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="recovered")
    a = ap.parse_args(argv)
    files = [FileRecord.from_path(p) for p in a.files]
    cs = make_heterogeneous(a.containers, np.random.default_rng(a.seed), 1.0)
    out = store_and_retrieve(files, cs, a.placement, a.controller, k=a.k, parity=a.parity, seed=a.seed)
    print(f"status={out['status']} reads={out['reads']:.0f} cost={out['cost']:.0f} actions={out['actions']}")
    print("library:", out["library"])
    od = Path(a.out)
    od.mkdir(exist_ok=True)
    for f in out["result"].files:
        (od / f.name).write_bytes(f.data)
        print(f"  wrote {od / f.name} ({len(f.data)} bytes) status={out['result'].file_status[f.name].value}")


if __name__ == "__main__":
    main()
