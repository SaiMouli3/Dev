"""AUC (Mann-Whitney) of the posterior-predictive P(complete) for predicting completion within the cap,
with bootstrap 95% CIs (1,000 resamples, seed 0). Input: results/paper_ec_confidence_calibration.csv."""
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata

RES = Path(__file__).resolve().parents[3] / "results"

def auc(y, p):
    r = rankdata(p); n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

d = pd.read_csv(RES / "paper_ec_confidence_calibration.csv")
rows = []
for (b, s), g in d.groupby(["belief", "spend_fraction"]):
    y = g.completed.values.astype(int); p = g.p_complete.values; rng = np.random.default_rng(0)
    bs = []
    for _ in range(1000):
        i = rng.integers(0, len(y), len(y))
        if 0 < y[i].sum() < len(i): bs.append(auc(y[i], p[i]))
    lo, hi = np.percentile(bs, [2.5, 97.5])
    rows.append(dict(belief=b, spend_fraction=s, n=len(g), completion_rate=y.mean(), mean_p=p.mean(), auc=auc(y, p), auc_lo=lo, auc_hi=hi))
out = pd.DataFrame(rows); out.to_csv(RES / "paper_ec_discrimination.csv", index=False); print(out.round(3).to_string())
