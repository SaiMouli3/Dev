"""Run the full D3RA-DNA experimental programme and write results/*.csv.

Usage:  python experiments/run_all.py [--quick] [--only baselines,h1,...]
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from d3ra.evaluation import experiments as E  # noqa: E402
from d3ra.evaluation.scenarios import SCENARIOS  # noqa: E402

RES = ROOT / "results"


def save(rows, name):
    RES.mkdir(exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(RES / f"{name}.csv", index=False)
    try:
        df.to_parquet(RES / f"{name}.parquet", index=False)
    except Exception:
        pass
    print(f"[saved] {name}: {len(df)} rows", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="small trial counts (smoke run)")
    ap.add_argument("--only", default="", help="comma list: baselines,h1,h5,ablations,redundancy,replication,placement,scalability")
    ap.add_argument("--scenarios", default=",".join(SCENARIOS))
    ap.add_argument("--tag", default="", help="suffix for output file names (parallel runs)")
    a = ap.parse_args(argv)
    q = a.quick
    tr, st = (300, 3000) if q else (2000, 20000)
    only = set(filter(None, a.only.split(",")))
    run = lambda k: not only or k in only
    meta = {"started": time.ctime(), "quick": q, "adaptive_trials": tr, "static_trials": st}
    if run("baselines"):
        rows = []
        for sc in a.scenarios.split(","):
            t = time.time()
            rows += E.baseline_experiment(sc, trials=tr, static_trials=st)
            print(f"  {sc}: {time.time() - t:.0f}s", flush=True)
            save(rows, "baselines" + a.tag)
    if run("h1"):
        save(E.h1_model_accuracy(trials=st), "h1_model_accuracy")
    if run("h5"):
        save(E.h5_robustness(trials=tr), "h5_robustness")
    if run("ablations"):
        rows = []
        for sc in ["s10_mixed", "s4_correlated_failures"]:
            rows += E.ablations(sc, trials=tr)
        save(rows, "ablations")
    if run("redundancy"):
        save(E.redundancy_experiment(), "redundancy")
    if run("replication"):
        save(E.replication_vs_parity(trials=st), "replication_vs_parity")
    if run("placement"):
        save(E.placement_optimizer_comparison(), "placement_opt")
    if run("threshold"):
        save(E.spreading_threshold(trials=2000 if q else 5000), "spreading_threshold")
    if run("scalability"):
        save(E.scalability(adaptive_trials=100 if q else 200), "scalability")
    meta["finished"] = time.ctime()
    (RES / f"run_meta{a.tag}.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
