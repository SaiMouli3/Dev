"""Revision-3 claim audit: recompute every headline number of main.tex directly from the committed
result files (no cached summaries other than results/dt_summary.csv, which dt_analysis.py
regenerates from the per-retrieval CSVs) and write audit/claim_recompute_rev3.txt.

    python research_paper/audit/recompute_claims.py
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, rankdata

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results"
out = []


def p(*a):
    out.append(" ".join(str(x) for x in a))


def cp_upper(k, n, a=0.05):
    return 1.0 if k >= n else float(beta.ppf(1 - a / 2, k + 1, n - k))  # upper limit of two-sided 95% CP interval


S = pd.read_csv(R / "dt_summary.csv")
g = lambda sc, po, col: float(S[(S.scenario == sc) & (S.policy == po)][col].iloc[0])

p("== Safety (all dt_*.csv except forecast/calibration runs) ==")
n = w = fa = att = 0
for f in sorted(glob.glob(str(R / "dt_*.csv"))):
    if any(k in f for k in ("summary", "forecast", "calibration")):
        continue
    d = pd.read_csv(f)
    n += len(d); w += int(d.wrong_output.sum()); fa += int(d.false_accepts.sum()); att += int(d.decode_attempts.sum())
p(f"retrievals={n} wrong={w} UB95={cp_upper(w, n):.2e} decode_attempts={att} false_accepts={fa} UB95={cp_upper(fa, att):.2e}")

p("\n== Baseline DT ==")
rr = g("baseline", "round_robin_stop", "mean_reads")
st = {a: g("baseline", f"static_a{a}", "q99_reads") for a in ("0.0", "0.5", "1.0", "2.0")}
p(f"RR mean={rr:.0f} CI=[{g('baseline','round_robin_stop','mean_lo'):.0f},{g('baseline','round_robin_stop','mean_hi'):.0f}] q99={g('baseline','round_robin_stop','q99_reads'):.0f}")
p("static q99 (= capacity a non-stopping plan must reserve):", st, f"best={min(st.values()):.0f}")
p(f"  consumed-vs-reserved saving of stopping: RR consumed {rr:.0f} vs best reservation {min(st.values()):.0f} "
  f"= {1 - rr / min(st.values()):.1%}; vs uniform reservation {st['0.0']:.0f} = {1 - rr / st['0.0']:.1%}")
p("  NOTE: this is a consumed-vs-reserved comparison, not an allocation comparison; for allocation see "
  "r4_actual_vs_reserved.csv, where the static plans also stop.")
for sc in ("baseline", "fail1", "homog"):
    for po in ("coverage", "pivotal", "thompson", "ucb", "learn_then_commit", "info_greedy", "tail_aware"):
        r = S[(S.scenario == sc) & (S.policy == po)].iloc[0]
        p(f"{sc:9s} {po:18s} success={r.success:.2f} mean={r.mean_reads:.0f} q99={r.q99_reads:.0f} "
          f"rel_vs_RR={r.rel_vs_rr:+.1%} [{r.rel_vs_rr_lo:+.1%},{r.rel_vs_rr_hi:+.1%}] n_joint={r.n_joint_success:.0f}")
    p(f"{sc:9s} RR mean={g(sc,'round_robin_stop','mean_reads'):.0f} q99={g(sc,'round_robin_stop','q99_reads'):.0f} "
      f"success={g(sc,'round_robin_stop','success'):.2f}; static q99 best="
      f"{min(g(sc, f'static_a{a}', 'q99_reads') for a in ('0.0','0.5','1.0','2.0')):.0f}")

d = pd.read_csv(R / "dt_baseline.csv")
p("\nread share by policy (baseline; mean of per-retrieval shares; containers 0..5, 5 = near-dead):")
for po, grp in d.groupby("policy"):
    a = np.stack(grp.reads_per_container.apply(lambda x: np.array(json.loads(x))))
    p(f"  {po:18s}", np.round((a / a.sum(1, keepdims=True)).mean(0), 3).tolist())

p("\n== fail2 / photolitho success ==")
for sc in ("fail2", "photolitho"):
    sub = S[S.scenario == sc]; p(sc, "max success over all policies and static plans =", sub.success.max(),
                                 "wrong outputs =", int(sub.wrong_outputs.sum()))

p("\n== Sweep (mean, q99) ==")
for c in (200, 30, 10, 3, 1):
    sub = S[S.scenario == f"sweep_cov{c}"]
    p(f"cov{c}:", {r.policy: (round(r.mean_reads), round(r.q99_reads)) for r in sub.itertuples()})

p("\n== Post-hoc rho=0 ==")
for sc in ("baseline", "fail1"):
    d0 = pd.read_csv(R / f"dt_{sc}_posthoc_rho0.csv"); dr = pd.read_csv(R / f"dt_{sc}.csv")
    rrm = dr[(dr.policy == "round_robin_stop") & (dr.status == "exact")]
    for po, grp in d0.groupby("policy"):
        ok = grp[grp.status == "exact"]
        j = ok.merge(rrm, on=["rep", "perm"], suffixes=("", "_rr"))
        p(f"{sc} {po:12s} mean(success)={ok.reads.mean():.0f} paired vs RR={j.reads.sum() / j.reads_rr.sum() - 1:+.1%} n={len(j)} q99={np.quantile(grp.reads, .99):.0f}")

p("\n== Calibration ==")
ce = pd.read_csv(R / "dt_calibration_eval_baseline.csv")
for rho, grp in ce.groupby("rho"):
    p(f"rho={rho}: coverage={grp.covered.mean():.3f}", grp.groupby("op").covered.mean().round(3).to_dict(),
      "by container", grp.groupby("container").covered.mean().round(2).to_dict())
cs = pd.read_csv(R / "dt_calibration_select_baseline.csv")
p("selection:", cs.groupby("rho").covered.mean().round(3).to_dict() if "rho" in cs else cs.head())

p("\n== Forecast ==")
fc = pd.read_csv(R / "dt_forecast_baseline.csv")
for f_, grp in fc.groupby("forecast"):
    y = grp.completed.astype(float).values; q = grp.p_complete.values
    bins = np.clip((q * 10).astype(int), 0, 9)
    ece = sum(abs(q[bins == b].mean() - y[bins == b].mean()) * (bins == b).mean() for b in range(10) if (bins == b).any())
    r = rankdata(q); npos = y.sum(); auc = (r[y == 1].sum() - npos * (npos + 1) / 2) / (npos * (len(y) - npos))
    p(f"{f_}: n={len(y)} base={y.mean():.2f} ECE={ece:.2f} Brier={np.mean((q - y) ** 2):.3f} AUC={auc:.2f} "
      f"climatology Brier={np.mean((y.mean() - y) ** 2):.3f}")
    for fr, gg in grp.groupby("frac"):
        p(f"   frac={fr}: mean forecast={gg.p_complete.mean():.2f} observed={gg.completed.mean():.2f}")

p("\n== Replay ==")
for name in ("replay_v3.csv", "replay_v3_genscript.csv"):
    fp = R / name
    if not fp.exists():
        p(name, "MISSING"); continue
    rp = pd.read_csv(fp)
    su = float(rp[rp.policy == "static_uniform"].mean_reads.iloc[0]); sh = float(rp[rp.policy == "static_hindsight"].mean_reads.iloc[0])
    rrr = rp[(rp.policy == "round_robin_stop") | (rp.policy == "round_robin")]
    p(f"{name}: static_uniform={su:.0f} hindsight={sh:.0f}")
    for r in rp.itertuples():
        if r.belief == "--":
            continue
        p(f"   {r.policy:14s} {r.belief:20s} mean={r.mean_reads:.0f} q99={r.q99_reads:.0f} vs uniform={r.mean_reads / su - 1:+.1%} "
          f"vs hindsight={r.mean_reads / sh - 1:+.1%} diff_vs_hindsight={r.diff_vs_hindsight_static:.0f} [{r.diff_lo:.0f},{r.diff_hi:.0f}]")

rc = (ROOT / "research_paper/evidence/verification_reports/rerun_check.txt").read_text().splitlines()[1]
p("\n== Re-run check ==\n" + rc)

txt = "\n".join(out) + "\n"
(ROOT / "research_paper/audit/claim_recompute_rev3.txt").write_text(txt)
print(txt)
