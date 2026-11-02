"""Q5: sensitivity of the key simulation comparisons to the scenario seed
(container randomness, trial randomness)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from d3ra.evaluation.experiments import baseline_experiment  # noqa: E402

rows = []
for sc in ["s2_heterogeneous", "s4_correlated_failures"]:
    for seed in range(5):
        r = baseline_experiment(sc, trials=1000, static_trials=10000, seed=seed)
        d = {x["design"]: x for x in r}
        rows.append({"scenario": sc, "seed": seed,
                     **{f"{k}_spend99": d[k]["mean_spend@0.99"] for k in ["A_uniform_fixed", "C_uniform_adaptive", "D_reliability_adaptive", "E_joint_adaptive"]},
                     **{f"{k}_pexact": d[k]["p_exact@ref"] for k in ["A_uniform_fixed", "C_uniform_adaptive", "D_reliability_adaptive", "E_joint_adaptive"]},
                     "E_variant": d["E_joint_adaptive"]["variant"]})
        print(rows[-1], flush=True)
        pd.DataFrame(rows).to_csv(ROOT / "results/q5_seed_sensitivity.csv", index=False)
