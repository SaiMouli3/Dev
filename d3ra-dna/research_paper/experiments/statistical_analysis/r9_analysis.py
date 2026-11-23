"""Final-revision analyses: stopping benchmarks and model-free allocation baselines.

1. Reservation baselines on the committed primary retrievals (results/dt_baseline*.csv, 30 x 10): a plan that
   does not stop and commits a fixed budget, set at the q99, at 1.2 x the median or at 1.2 x the mean of the
   fixed-depth read requirement.  Reported: the budget, the share of retrievals it completes, and the expected
   reads per completed retrieval if a failed attempt is repeated with a fresh budget (budget / success).
2. Exploratory baselines on regenerated read streams (results/dt_{baseline,homog,fail1}_r9regen.csv):
   round-robin with verified stopping (the reference, re-run on the same streams), the same schedule with the
   stopping test applied after every read (per-read stop), successive elimination, and a yield-threshold rule.
   DT4DDS is not seed-deterministic, so these streams differ from the committed ones; every comparison is paired
   within the regenerated set.

    python research_paper/experiments/statistical_analysis/r9_analysis.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bootstrap import cluster_boot  # noqa: E402

RES = Path(__file__).resolve().parents[3] / "results"


def reservations():
    st = pd.read_csv(RES / "dt_baseline_static.csv")
    rows = []
    for pol in ("static_a0.0", "static_a1.0"):
        R = st[st.policy == pol]
        for rule, f in (("q99", lambda x: np.quantile(x, .99)), ("1.2 x median", lambda x: 1.2 * np.median(x)),
                        ("1.2 x mean", lambda x: 1.2 * np.mean(x))):
            bud = f(R.reads.values)
            succ = cluster_boot(R, lambda x, b=bud: float((x.reads <= b).mean()))
            rows.append({"plan": pol, "rule": rule, "budget": bud, "success": succ[0], "success_lo": succ[1],
                         "success_hi": succ[2], "reads_per_completed": bud / succ[0]})
    return pd.DataFrame(rows)


def extra_baselines():
    rows = []
    for scen in ("baseline", "homog", "fail1"):
        f = RES / f"dt_{scen}_r9regen.csv"
        if not f.exists():
            continue
        d = pd.read_csv(f)
        rr = d[d.policy == "round_robin_stop"]
        for pol, g in d.groupby("policy"):
            ok = g.status == "exact"
            rpc = np.stack(g.reads_per_container.map(json.loads).values)
            m = cluster_boot(g[ok], lambda x: x.reads.mean()) if ok.any() else (np.nan,) * 3
            row = {"scenario": scen, "policy": pol, "n": len(g), "success": ok.mean(),
                   "wrong_outputs": int(g.wrong_output.sum()), "false_accepts": int(g.false_accepts.sum()),
                   "mean_success": m[0], "mean_lo": m[1], "mean_hi": m[2],
                   "q99_all": np.quantile(g.reads, .99), "switches": g.switches.mean(),
                   "share_last_container": float((rpc[:, -1] / rpc.sum(1)).mean())}
            if pol != "round_robin_stop":
                j = g.merge(rr, on=["rep", "perm"], suffixes=("", "_r"))   # every synthesis is a cluster
                both = (j.status == "exact") & (j.status_r == "exact")
                e = cluster_boot(j, lambda x: x.reads[(x.status == "exact") & (x.status_r == "exact")].sum()
                                 / x.reads_r[(x.status == "exact") & (x.status_r == "exact")].sum() - 1)
                row.update({"rel_vs_rr": e[0], "lo": e[1], "hi": e[2], "n_joint": int(both.sum())})
            rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    a = reservations()
    a.to_csv(RES / "r9_reservations.csv", index=False)
    b = extra_baselines()
    b.to_csv(RES / "r9_extra_baselines.csv", index=False)
    pd.set_option("display.width", 220)
    print(a.round(3).to_string())
    print(b.round(4).to_string())
