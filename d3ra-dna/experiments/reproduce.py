"""Independent reproduction of headline claims with explicit tolerances.

  python experiments/reproduce.py            # simulation + theory claims (~2-5 min, CPU)
  python experiments/reproduce.py --real     # + real-data claims (needs `real_data.py download prep channel`)

Each check prints PASS/FAIL with the measured value and the claimed value."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

results = []


def check(name, value, ok, claim):
    results.append({"check": name, "value": value, "pass": bool(ok), "claim": claim})
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {value}  (claim: {claim})", flush=True)


def theory():
    from d3ra.theory.spreading import even_split, homogeneous_crossover, loss_probability, tolerance
    agree = tot = 0
    for D in range(2, 9):
        for N in range(D * 4, 61, D):
            for K in range(N // 3, N):
                n = even_split(N, D)
                t = tolerance(n, K)
                qs = homogeneous_crossover(D, t)
                for q in (0.005, 0.02, 0.08, 0.3):
                    p = loss_probability(n, K, q)
                    if abs(p - q) < 1e-12:
                        continue
                    tot += 1
                    agree += (p < q) == (q < qs)
    check("Theorem 1 sign prediction (D | N)", f"{agree}/{tot}", agree == tot, "100% agreement")


def predictor_vs_mc():
    from d3ra.channel.surrogate import ChannelParams, Layout, Predictor, static_budget_samples
    from d3ra.storage.container import make_heterogeneous
    rng = np.random.default_rng(0)
    ch = ChannelParams.from_containers(make_heterogeneous(6, rng, 1.0))
    S, K, N = 20, 30, 40
    pl = np.zeros((S * N, 6)); pl[np.arange(S * N), np.arange(S * N) % 6] = 1
    L = Layout(np.repeat(np.arange(S), N), np.full(S, K), pl)
    w = np.ones(6) / 6
    bs = static_budget_samples(L, ch, w, 10000, rng)
    pr = Predictor(L, ch)
    err = max(abs(pr.p_file(B * w) - (bs <= B).mean()) for B in (5000, 8000, 15000))
    check("analytic predictor vs Monte Carlo", round(err, 4), err < 0.02, "max |dP| < 0.02")


def adaptive_vs_fixed():
    from d3ra.evaluation.experiments import baseline_experiment
    rows = baseline_experiment("s1_homogeneous", trials=600, static_trials=6000)
    d = {r["design"]: r for r in rows if not r["design"].startswith("cmp")}
    a, c = d["A_uniform_fixed"]["mean_spend@0.99"], d["C_uniform_adaptive"]["mean_spend@0.99"]
    check("s1 adaptive stopping saving vs fixed depth", f"{1 - c / a:.3f}", 0.3 < 1 - c / a < 0.65, "0.48 at full trial counts (0.44 at these reduced counts)")


def threshold():
    from d3ra.evaluation.experiments import spreading_threshold
    import pandas as pd
    df = pd.DataFrame(spreading_threshold(parities=(12, 24), trials=3000))
    p = df.pivot(index="parity", columns="placement", values="file_infeasible_mc")
    ok = p.loc[12, "failure0"] < p.loc[12, "domain"] and p.loc[24, "domain"] < 0.03
    check("spreading threshold direction (s4)", p.round(4).to_dict(), ok,
          "parity 12: concentration safer than spreading; parity 24: spreading <3% loss")


def real():
    import pandas as pd
    ch = pd.read_csv(ROOT / "results/realdata_channel.csv")
    v = float(ch.loc[ch.condition == "Aging_0d", "phi_exact_n3"].iloc[0])
    check("real Illumina phi(3), fresh pool", round(v, 3), abs(v - 0.977) < 0.02, "0.977")
    sv = json.loads((ROOT / "results/realdata_shared_variance.json").read_text())["all_runs"]["shared_fraction"]
    check("shared coverage-variance fraction", round(sv, 3), abs(sv - 0.69) < 0.05, "0.69")
    pol = ROOT / "results/realdata_policies.csv"
    if pol.exists():
        d = pd.read_csv(pol)
        a = d[(d.placement == "round_robin") & (d.policy == "static_uniform")]["mean_spend@0.99"].iloc[0]
        e = d[(d.placement == "single_best") & (d.policy == "adaptive_pivotal")]["mean_spend@0.99"].iloc[0]
        check("real-data placement+adaptive saving", f"{1 - e / a:.3f}", 0.25 < 1 - e / a < 0.42, "0.34")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", action="store_true")
    a = ap.parse_args()
    theory()
    predictor_vs_mc()
    threshold()
    adaptive_vs_fixed()
    if a.real:
        real()
    (ROOT / "results/reproduce_report.json").write_text(json.dumps(results, indent=1, default=str))
    print(f"{sum(r['pass'] for r in results)}/{len(results)} checks passed")
    sys.exit(0 if all(r["pass"] for r in results) else 1)
