"""Analysis of the digital-twin end-to-end experiments (frozen protocol, dt_protocol.yaml).

For every scenario and policy:
  * success = output file byte-identical to the stored file; wrong outputs (file returned but not
    identical) with an exact Clopper-Pearson 95% upper bound;
  * false accepts per decode attempt (exact 95% upper bound);
  * mean reads per retrieval (pay-per-read cost) and q99 of reads (reserved-capacity budget);
    static plans pay their q99 budget on every retrieval;
  * CIs: cluster bootstrap over synthesis replicates (2,000 resamples);
  * paired differences vs uniform fixed depth and vs the best (hindsight) static plan, resampling
    replicates jointly for both policies (same syntheses and read orders).

    python research_paper/experiments/statistical_analysis/dt_analysis.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

RES = Path(__file__).resolve().parents[3] / "results"
B = 2000


def cp_upper(k, n, a=0.05):
    return float(beta.ppf(1 - a / 2, k + 1, n - k)) if n > 0 and k < n else 1.0


def q99(x):
    x = np.sort(np.asarray(x, float))
    return float(x[min(len(x) - 1, int(np.ceil(0.99 * len(x))) - 1)])


def _matrix(df, pol):
    """reads as (rep, perm) array for one policy (inf where the file was not recovered)."""
    d = df[df.policy == pol]
    v = np.where(d.file_identical, d.reads, np.inf)
    reps = sorted(d.rep.unique())
    perms = sorted(d.perm.unique())
    M = np.full((len(reps), len(perms)), np.inf)
    for (r, p), x in zip(zip(d.rep, d.perm), v):
        M[reps.index(r), perms.index(p)] = x
    return M


def boot(stat, mats, rng):
    n = mats[0].shape[0]
    out = []
    for _ in range(B):
        ix = rng.integers(0, n, n)
        out.append(stat(*[m[ix].ravel() for m in mats]))
    return np.percentile(out, [2.5, 97.5])


def analyse(scen):
    f = RES / f"dt_{scen}.csv"
    if not f.exists():
        return []
    df = pd.read_csv(f)
    fs = RES / f"dt_{scen}_static.csv"
    if fs.exists():
        st = pd.read_csv(fs)
        df = pd.concat([df, st], ignore_index=True)
    rng = np.random.default_rng(0)
    mats = {p: _matrix(df, p) for p in df.policy.unique()}
    statics = [p for p in mats if p.startswith("static_")]
    best_static = min(statics, key=lambda p: q99(mats[p].ravel())) if statics else None
    rows = []
    for pol, Mx in mats.items():
        d = df[df.policy == pol]
        x = Mx.ravel()
        n = len(d)
        succ = int(d.file_identical.sum())
        wrong = int(d.wrong_output.sum())
        fa, att = int(d.false_accepts.sum()), int(d.decode_attempts.sum())
        is_static = pol.startswith("static_")
        mean = q99(x) if is_static else float(np.mean(x))
        mean_ci = boot(q99 if is_static else np.mean, [Mx], rng)
        row = {"scenario": scen, "policy": pol, "retrievals": n, "success": succ / n, "wrong_outputs": wrong,
               "wrong_output_upper95": cp_upper(wrong, n), "false_accepts": fa, "decode_attempts": att,
               "false_accept_upper95": cp_upper(fa, att), "mean_reads": mean, "mean_lo": mean_ci[0],
               "mean_hi": mean_ci[1], "mean_success": float(np.mean(x[np.isfinite(x)])) if np.isfinite(x).any() else np.nan,
               "q99_reads": q99(x), "q99_lo": boot(q99, [Mx], rng)[0],
               "q99_hi": boot(q99, [Mx], rng)[1], "best_static": best_static}
        for ref in ("static_a0.0", best_static):
            if ref is None or ref not in mats or ref == pol:
                continue
            R = mats[ref]
            if is_static:
                st_ = lambda a, r: q99(a) - q99(r)
            else:
                st_ = lambda a, r: float(np.mean(a)) - q99(r)
            key = "fixed" if ref == "static_a0.0" else "beststatic"
            row[f"diff_vs_{key}"] = st_(x, R.ravel())
            lo, hi = boot(st_, [Mx, R], rng)
            row[f"diff_vs_{key}_lo"], row[f"diff_vs_{key}_hi"] = lo, hi
            row[f"q99diff_vs_{key}"] = q99(x) - q99(R.ravel())
            lo, hi = boot(lambda a, r: q99(a) - q99(r), [Mx, R], rng)
            row[f"q99diff_vs_{key}_lo"], row[f"q99diff_vs_{key}_hi"] = lo, hi
        # paired vs round-robin + verified stop, on cells where both recovered the file
        # (success is identical across policies in every scenario we ran; cost only compared where both paid to success)
        if not is_static and pol != "round_robin_stop" and "round_robin_stop" in mats:
            R = mats["round_robin_stop"]
            both = np.isfinite(Mx) & np.isfinite(R)
            row["n_joint_success"] = int(both.sum())
            if both.any():
                def st_rr(a, r):
                    m = np.isfinite(a) & np.isfinite(r)
                    return float(np.mean(a[m]) / np.mean(r[m]) - 1) if m.any() else np.nan
                row["rel_vs_rr"] = st_rr(Mx.ravel(), R.ravel())
                lo, hi = np.nanpercentile([st_rr(*(m[ix].ravel() for m in (Mx, R)))
                                           for ix in (rng.integers(0, Mx.shape[0], Mx.shape[0]) for _ in range(B))], [2.5, 97.5])
                row["rel_vs_rr_lo"], row["rel_vs_rr_hi"] = lo, hi
        rows.append(row)
    return rows


if __name__ == "__main__":
    scens = ["baseline", "fail1", "fail2", "homog", "photolitho",
             "sweep_cov200", "sweep_cov30", "sweep_cov10", "sweep_cov3", "sweep_cov1"]
    rows = []
    for s in scens:
        rows += analyse(s)
    out = pd.DataFrame(rows)
    out.to_csv(RES / "dt_summary.csv", index=False)
    pd.set_option("display.width", 250)
    cols = ["scenario", "policy", "success", "wrong_outputs", "mean_reads", "q99_reads",
            "diff_vs_beststatic", "rel_vs_rr", "rel_vs_rr_lo", "rel_vs_rr_hi"]
    print(out[[c for c in cols if c in out]].round(3).to_string())
