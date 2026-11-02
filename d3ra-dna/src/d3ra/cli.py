"""D3RA-DNA command line.

  d3ra store FILES... --archive DIR [--containers 6] [--placement tiered|domain|failure|round_robin]
  d3ra simulate --archive DIR --out FASTQ_DIR [--reads 20000] [--transplant illumina|nanopore] [--fail 0,3]
  d3ra retrieve --archive DIR --fastq FASTQ_DIR [--controller pivotal|uniform|coverage] --out OUT_DIR
  d3ra demo FILES...            (in-memory store + simulated retrieve, as before)
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .ingestion.file_loader import FileRecord


def _containers(n, seed):
    from .storage.container import make_heterogeneous
    return make_heterogeneous(n, np.random.default_rng(seed), 1.0)


def cmd_store(a):
    from .archive import save_archive
    from .channel.failures import FailureModel
    from .channel.surrogate import ChannelParams
    from .storage.library import build_library
    from .storage.placement import cost_tiered, domain_aware, failure_aware, round_robin, unit_retrieval_cost
    files = [FileRecord.from_path(p) for p in a.files]
    lib = build_library(files, a.payload, a.k, a.parity, a.nsym, a.codec)
    cs = _containers(a.containers, a.seed)
    ch = ChannelParams.from_containers(cs, FailureModel(np.array([c.fail_prob for c in cs])), L=lib.fmt.strand_nt)
    so, k = lib.stripe_of, lib.k_per_stripe()
    dom = np.array([c.domain for c in cs])
    X = {"tiered": lambda: cost_tiered(so, k, unit_retrieval_cost(ch), dom, slack=2),
         "domain": lambda: domain_aware(so, k, dom, len(cs)),
         "failure": lambda: failure_aware(so, k, ch.failure, (1 - ch.dropout) * (1 - ch.d_syn), unit_retrieval_cost(ch), 0.05),
         "round_robin": lambda: round_robin(so, len(cs))}[a.placement]()
    save_archive(lib, X, cs, a.archive)
    print(json.dumps({"archive": a.archive, **lib.stats(), "per_container": X.sum(0).astype(int).tolist()}, indent=1))


def _transplant(kind):
    import pickle
    root = Path(__file__).resolve().parents[2] / "data" / "processed" / "realdata"
    f = root / f"transplant_{kind}.pkl"
    if not f.exists():
        raise SystemExit(f"{f} missing: run experiments/real_data.py transplant first")
    with open(f, "rb") as fh:
        return pickle.load(fh)


def cmd_simulate(a):
    from .retrieval.fastq_retrieval import simulate_fastq
    em = None
    if a.transplant:
        t = _transplant(a.transplant)
        em = {j: t for j in range(1000)}
    failed = None
    if a.fail:
        from .archive import load_archive
        m = load_archive(a.archive)[2].shape[1]
        failed = np.zeros(m, bool)
        failed[[int(x) for x in a.fail.split(",")]] = True
    paths = simulate_fastq(a.archive, a.out, a.reads, em, a.seed, failed)
    print("\n".join(str(p) for p in paths))


def cmd_retrieve(a):
    from .retrieval.fastq_retrieval import FastqRetriever
    from .archive import load_archive
    m = load_archive(a.archive)[2].shape[1]
    fq = [Path(a.fastq) / f"container_{j:02d}.fq.gz" for j in range(m)]
    r = FastqRetriever(a.archive, [f if f.exists() else None for f in fq], a.controller).run()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for f in r["result"].files:
        (out / f.name).write_bytes(f.data)
    print(json.dumps({k: v for k, v in r.items() if k != "result"}, indent=1))
    for f in r["result"].files:
        print(f"  {f.name}: {r['result'].file_status[f.name].value} ({len(f.data)} bytes) -> {out / f.name}")


def cmd_demo(a):
    from .pipeline import store_and_retrieve
    files = [FileRecord.from_path(p) for p in a.files]
    out = store_and_retrieve(files, _containers(a.containers, a.seed), "tiered", a.controller, seed=a.seed)
    print(out["status"], out["reads"], out["library"])


def main(argv=None):
    ap = argparse.ArgumentParser(prog="d3ra")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("store")
    s.add_argument("files", nargs="+")
    s.add_argument("--archive", required=True)
    s.add_argument("--containers", type=int, default=6)
    s.add_argument("--placement", default="tiered", choices=["tiered", "domain", "failure", "round_robin"])
    s.add_argument("--k", type=int, default=24)
    s.add_argument("--parity", type=int, default=12)
    s.add_argument("--payload", type=int, default=32)
    s.add_argument("--nsym", type=int, default=4)
    s.add_argument("--codec", default="scrambled-2bit")
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(fn=cmd_store)
    s = sub.add_parser("simulate")
    s.add_argument("--archive", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--reads", type=int, default=20000)
    s.add_argument("--transplant", choices=["illumina", "nanopore"])
    s.add_argument("--fail", default="")
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(fn=cmd_simulate)
    s = sub.add_parser("retrieve")
    s.add_argument("--archive", required=True)
    s.add_argument("--fastq", required=True)
    s.add_argument("--controller", default="pivotal", choices=["pivotal", "uniform", "coverage"])
    s.add_argument("--out", default="recovered")
    s.set_defaults(fn=cmd_retrieve)
    s = sub.add_parser("demo")
    s.add_argument("files", nargs="+")
    s.add_argument("--containers", type=int, default=6)
    s.add_argument("--controller", default="pivotal")
    s.add_argument("--seed", type=int, default=0)
    s.set_defaults(fn=cmd_demo)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
