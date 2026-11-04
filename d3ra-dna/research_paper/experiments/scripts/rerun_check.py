"""Bit-for-bit re-run check: re-execute a sample of committed digital-twin retrievals in a fresh
process (single worker, different job order) and compare every recorded field.

    python research_paper/experiments/scripts/rerun_check.py [scenario] [n]

Writes research_paper/evidence/verification_reports/rerun_check.txt.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJ = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PROJ / "experiments"))
import dt_experiments as DT  # noqa: E402

FIELDS = ["status", "reads", "reads_per_container", "rounds", "decode_attempts", "false_accepts",
          "verification_rounds", "file_identical", "wrong_output"]


def main(scen="baseline", n=24):
    df = pd.read_csv(PROJ / "results" / f"dt_{scen}.csv")
    sample = df.sample(n=n, random_state=0).iloc[::-1]          # reversed order on purpose
    lines, mism = [], 0
    for _, row in sample.iterrows():
        new = DT.one_retrieval((scen, int(row.rep), int(row.perm), row.policy, {}, np.inf))
        new["reads_per_container"] = json.dumps(new["reads_per_container"])
        diff = [f for f in FIELDS if str(new[f]) != str(row[f])]
        mism += bool(diff)
        lines.append(f"{scen} rep={row.rep} perm={row.perm} {row.policy:18s} reads={row.reads:.0f} "
                     f"-> {new['reads']:.0f} {'IDENTICAL' if not diff else 'DIFF ' + ','.join(diff)}")
    head = (f"Re-run check of {n} randomly sampled committed retrievals of results/dt_{scen}.csv in a fresh, "
            f"single-process run (reverse order). Fields compared: {', '.join(FIELDS)}.\n"
            f"Mismatching retrievals: {mism} of {n}.\n")
    out = PROJ / "research_paper" / "evidence" / "verification_reports" / "rerun_check.txt"
    out.write_text(head + "\n".join(lines) + "\n")
    print(head)
    return mism


if __name__ == "__main__":
    sys.exit(main(*(sys.argv[1:2] or ["baseline"]), *(int(a) for a in sys.argv[2:3])) > 0)
