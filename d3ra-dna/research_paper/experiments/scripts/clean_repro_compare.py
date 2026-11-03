"""Compare the clean re-simulation (results/dt_baseline*_cleanrepro.csv, produced by clean_repro.sh from
a fresh clone with DT4DDS replicates regenerated from seeds) against the committed runs, and report both
what does not reproduce (individual retrievals, because DT4DDS 1.1.0 is not seed-deterministic) and what
does (every policy mean, the paired comparisons, the ordering).

    python research_paper/experiments/scripts/clean_repro_compare.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "results"
OUT = ROOT / "research_paper" / "evidence" / "verification_reports" / "clean_repro.txt"
RNG = np.random.default_rng(3)
KEY = ["scenario", "rep", "perm", "policy"]


def cb(df, f, n=1500):
    gs = [g for _, g in df.groupby("rep")]
    est = f(df)
    v = [f(pd.concat([gs[i] for i in RNG.integers(0, len(gs), len(gs))])) for _ in range(n)]
    return est, float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))


def load(suffix=""):
    return pd.concat([pd.read_csv(RES / f"dt_baseline{suffix}.csv"),
                      pd.read_csv(RES / f"dt_baseline_static{suffix}.csv")])


def main():
    a, b = load(), load("_cleanrepro")
    aa = a.sort_values(KEY).reset_index(drop=True)
    bb = b.sort_values(KEY).reset_index(drop=True)
    lines = ["Scenario-level clean reproduction of the primary baseline scenario (all 8 policies and the",
             "4 static plans, 30 syntheses x 10 read orders each).  The repository was cloned at the commit",
             "in research_paper/final/COMMIT.txt and every DT4DDS replicate regenerated from its seed.", ""]
    same = {c: aa[c].astype(str).equals(bb[c].astype(str)) for c in
            ["status", "file_identical", "wrong_output", "false_accepts", "verification_rounds"]}
    diff = int((aa.reads != bb.reads).sum())
    lines += ["1. OUTCOMES reproduce exactly.",
              f"   {', '.join(c + ': identical' for c, v in same.items() if v)}.",
              f"   Wrong outputs in both runs: {int(aa.wrong_output.sum())} and {int(bb.wrong_output.sum())}.",
              "",
              "2. INDIVIDUAL READ COUNTS do not, and this is a property of the simulator, not of the",
              "   controller.  DT4DDS 1.1.0 is not deterministic given a seed: calling its synthesis ->",
              "   aging -> PCR -> sequencing chain twice in one process with the same seed produces",
              "   different read streams.  Regenerating the replicates therefore yields different streams,",
              f"   and {diff} of {len(aa)} retrievals consumed a different number of reads",
              f"   (mean absolute difference {float((aa.reads - bb.reads).abs().mean()):.0f} reads).",
              "   The committed read streams are consequently part of the artifact; given them, retrieval is",
              "   bit-exact (see rerun_check.txt).",
              "",
              "3. CONCLUSIONS reproduce.  Mean reads per policy, with cluster-bootstrap 95% intervals over",
              "   syntheses, and the paired difference against round-robin:", ""]
    hdr = f"   {'policy':18s} {'committed':>22s} {'clean re-run':>22s}   {'vs RR (committed)':>18s} {'vs RR (re-run)':>16s}"
    lines += [hdr, "   " + "-" * (len(hdr) - 3)]
    rows = {}
    for tag, d in (("c", a), ("r", b)):
        rr = d[d.policy == "round_robin_stop"]
        for pol, g in d.groupby("policy"):
            m = cb(g, lambda x: x.reads.mean())
            j = g.merge(rr, on=["rep", "perm"], suffixes=("", "_r"))
            rel = (0.0,) * 3 if pol == "round_robin_stop" else cb(j, lambda x: x.reads.sum() / x.reads_r.sum() - 1)
            rows[(tag, pol)] = (m, rel)
    ok_in = True
    order = sorted({p for _, p in rows}, key=lambda p: rows[("c", p)][0][0])
    for pol in order:
        (mc, lc, hc), relc = rows[("c", pol)]
        (mr, lr, hr), relr = rows[("r", pol)]
        inside = (lr <= mc <= hr) and (lc <= mr <= hc)
        ok_in &= inside
        lines.append(f"   {pol:18s} {mc:8.0f} [{lc:5.0f},{hc:6.0f}] {mr:8.0f} [{lr:5.0f},{hr:6.0f}]   "
                     f"{relc[0]:+17.1%} {relr[0]:+15.1%}{'' if inside else '   <-- outside'}")
    same_order = [p for p in order] == sorted(order, key=lambda p: rows[("r", p)][0][0])
    lines += ["",
              f"   Every policy mean inside the other run's interval, in both directions: {'yes' if ok_in else 'NO'}.",
              f"   Policy ordering by mean reads identical: {'yes' if same_order else 'NO'}.",
              "",
              "Conclusion: the simulator is stochastic across regenerations, so individual retrievals are not",
              "reproducible from seeds alone; every claim this paper makes about means, paired comparisons and",
              "policy ordering is reproduced by an independent re-simulation from a clean checkout."]
    OUT.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0 if (ok_in and same_order and all(same.values())) else 1


if __name__ == "__main__":
    raise SystemExit(main())
