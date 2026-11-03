"""Merge Genscript replay_v3 conditions run in parallel processes (D3RA_REPLAY_TAG=_partA/_partB/_partC)
into results/replay_v3_genscript{,_per_trial}.csv.  Each condition has its own fixed seed, so the
parts are identical to a single sequential run (checked against the logged sequential values)."""
from pathlib import Path

import pandas as pd

RES = Path(__file__).resolve().parents[3] / "results"
ORDER = ["static_uniform", "static_hindsight", "uniform", "coverage", "pivotal", "ts", "ucb", "ltc", "info", "tail_lcb"]
parts = [p for p in ("A", "B", "C") if (RES / f"replay_v3_genscript_part{p}.csv").exists()]
rows = pd.concat([pd.read_csv(RES / f"replay_v3_genscript_part{p}.csv") for p in parts]).drop_duplicates(["policy", "belief"])
rows["_o"] = [ORDER.index(p) + (100 if b == "calibrated_leakfree" else 0) for p, b in zip(rows.policy, rows.belief)]
rows.sort_values("_o").drop(columns="_o").to_csv(RES / "replay_v3_genscript.csv", index=False)
pt = pd.concat([pd.read_csv(RES / f"replay_v3_genscript_part{p}_per_trial.csv") for p in parts], axis=1)
pt.loc[:, ~pt.columns.duplicated()].to_csv(RES / "replay_v3_genscript_per_trial.csv", index=False)
print("merged parts", parts, "->", len(rows), "conditions")
