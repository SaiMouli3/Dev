"""Numerical verification of the spreading-threshold results (docs/theory_spreading.md)."""
import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from d3ra.theory.spreading import (best_tolerant_vs_concentration, even_split, homogeneous_crossover,  # noqa: E402
                                   homogeneous_spread_loss, loss_probability, no_tolerance_dominated, tolerance)

RES = ROOT / "results"
rng = np.random.default_rng(0)
rows = []

# T1: homogeneous theorem -- closed form vs exact enumeration vs Monte Carlo; sign prediction
agree = total = 0
for D in range(2, 9):
    for N in [24, 36, 48, 60]:
        for K in range(N // 3, N, max(1, N // 12)):
            n = even_split(N, D)
            t = tolerance(n, K)
            qstar = homogeneous_crossover(D, t)
            for q in [0.001, 0.01, 0.03, 0.1, 0.2, 0.4]:
                p_sp = loss_probability(n, K, q)
                p_cf = homogeneous_spread_loss(D, t, q) if N % D == 0 else np.nan
                sims = rng.random((20000, D)) < q
                mc = float(((~sims) @ n < K).mean())
                pred_better = q < qstar
                actual_better = p_sp < q
                if not np.isclose(p_sp, q):
                    total += 1
                    agree += pred_better == actual_better
                rows.append({"test": "T1", "D": D, "N": N, "K": K, "rho": (N - K) / N, "t": t, "q": q,
                             "q_star": qstar, "p_spread_exact": p_sp, "p_spread_closed": p_cf, "p_spread_mc": mc,
                             "p_concentrate": q, "spread_better": actual_better, "predicted_better": pred_better})
df = pd.DataFrame(rows)
summary = {"T1_sign_agreement": agree / total, "T1_cases": total,
           "T1_closed_form_max_abs_err": float(np.nanmax(np.abs(df.p_spread_exact - df.p_spread_closed))),
           "T1_mc_max_abs_err": float(np.max(np.abs(df.p_spread_exact - df.p_spread_mc))),
           "T1a_rho_below_1_over_D_spread_never_better": bool(not df[(df.t == 0)].spread_better.any())}

# T2: q* table
qs = [{"D": D, "t": t, "q_star": homogeneous_crossover(D, t)} for D in range(2, 9) for t in range(0, D)]
pd.DataFrame(qs).to_csv(RES / "theory_qstar.csv", index=False)

# T3: Proposition 2 (heterogeneous, random instances): no single-failure tolerance => no better than concentration
viol = 0
cnt = 0
for _ in range(3000):
    D = int(rng.integers(2, 7))
    N = int(rng.integers(8, 40))
    K = int(rng.integers(1, N))
    q = np.exp(rng.uniform(np.log(1e-3), np.log(0.4), D))
    n = rng.multinomial(N, rng.dirichlet(np.ones(D)))
    if (n > N - K).any():
        ok, pl, qmin = no_tolerance_dominated(n, K, q)
        cnt += 1
        viol += not ok
summary["P2_instances"] = cnt
summary["P2_violations"] = viol

# T4: dropout margin -- crossover parity for the s4-like scenario (D=4, K=32, delta=2%) by exact computation
trows = []
for p in range(4, 40, 2):
    N = 32 + p
    n = even_split(N, 4)
    for q in [0.01, 0.03, 0.05]:
        trows.append({"parity": p, "rho": p / N, "q": q, "spread": loss_probability(n, 32, q, delta=0.02),
                      "concentrate": loss_probability([N, 0, 0, 0], 32, q, delta=0.02)})
t4 = pd.DataFrame(trows)
t4["spread_better"] = t4.spread < t4.concentrate
t4.to_csv(RES / "theory_dropout_margin.csv", index=False)
cross = t4[t4.spread_better].groupby("q").rho.min().to_dict()
summary["T4_min_rho_where_spreading_wins_delta2pct"] = {str(k): v for k, v in cross.items()}

# T5: heterogeneous generalisation with capacities
hrows = []
for _ in range(200):
    D = 6
    q = np.sort(np.exp(rng.uniform(np.log(0.002), np.log(0.1), D)))
    N, K = 48, int(rng.integers(24, 44))
    cap = rng.integers(8, 49, D)
    r = best_tolerant_vs_concentration(N, K, q, cap)
    hrows.append({"K": K, "rho": (N - K) / N, **r})
h = pd.DataFrame(hrows)
h["tolerant_better"] = h.p_loss_tolerant < h.p_loss_concentration
h["comparable"] = h.p_loss_tolerant.notna() & h.p_loss_concentration.notna()
h.to_csv(RES / "theory_heterogeneous.csv", index=False)
summary["T5_share_tolerant_better_when_feasible"] = float(h[h.comparable].tolerant_better.mean())
summary["T5_rho_tolerant_better_median"] = float(h[h.comparable & h.tolerant_better].rho.median()) if (h.comparable & h.tolerant_better).any() else None
summary["T5_rho_concentrated_better_median"] = float(h[h.comparable & ~h.tolerant_better].rho.median()) if (h.comparable & ~h.tolerant_better).any() else None
summary["T5_infeasible_tolerance_share"] = float(h.p_loss_tolerant.isna().mean())

# T6: non-MDS codes: fountain needing K(1+eps) shifts the threshold
frows = []
for eps in [0.0, 0.05, 0.1, 0.2]:
    for p in range(4, 40, 4):
        N = 32 + p
        n = even_split(N, 4)
        frows.append({"eps": eps, "parity": p, "rho": p / N, "t": tolerance(n, int(np.ceil(32 * (1 + eps)))),
                      "spread": loss_probability(n, 32, 0.02, k_eff=32 * (1 + eps)),
                      "concentrate": loss_probability([N, 0, 0, 0], 32, 0.02, k_eff=32 * (1 + eps))})
f = pd.DataFrame(frows)
f["spread_better"] = f.spread < f.concentrate
summary["T6_min_rho_spread_wins_by_eps"] = {str(k): float(v) for k, v in f[f.spread_better].groupby("eps").rho.min().items()}
f.to_csv(RES / "theory_nonmds.csv", index=False)

df.to_csv(RES / "theory_homogeneous.csv", index=False)
(RES / "theory_summary.json").write_text(json.dumps(summary, indent=1))
print(json.dumps(summary, indent=1))
