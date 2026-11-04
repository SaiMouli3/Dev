"""Final numeric verification: recompute every quantitative claim in the manuscript directly from the
committed result files and compare with the value printed in the PDF.

    python research_paper/audit/verify_manuscript.py

Each check names the section it guards, the recomputed value, the value asserted in the manuscript and
the tolerance.  Exit status is non-zero if any check fails.  Output: audit/verify_manuscript.txt
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results"
PDF = ROOT / "research_paper" / "manuscript" / "main.pdf"
OUT = Path(__file__).resolve().parent / "verify_manuscript.txt"
RNG = np.random.default_rng(20261104)
rows: list[tuple] = []


def load(n):
    p = RES / n
    return pd.read_csv(p) if p.exists() else None


def cb(df, f, n=2000):
    gs = [g for _, g in df.groupby("rep")]
    v = [f(pd.concat([gs[i] for i in RNG.integers(0, len(gs), len(gs))])) for _ in range(n)]
    return f(df), float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))


def paired(g, rr, stat=lambda x: x.reads.sum() / x.reads_r.sum() - 1, ok_only=True):
    a = g[g.status == "exact"] if ok_only else g
    b = rr[rr.status == "exact"] if ok_only else rr
    j = a.merge(b, on=["rep", "perm"], suffixes=("", "_r"))
    return cb(j, stat) if len(j) else (np.nan, np.nan, np.nan)


def check(section, what, recomputed, claimed, tol, unit=""):
    ok = (recomputed is not None and claimed is not None
          and abs(recomputed - claimed) <= tol)
    rows.append((ok, section, what, recomputed, claimed, tol, unit))


def cp_upper(k, n, a=0.05):
    return 1.0 if k >= n else float(beta.ppf(1 - a / 2, k + 1, n - k))


# ---------------------------------------------------------------- safety (abstract, Sec. VI-A, Table IV)
files = [f for f in sorted(RES.glob("dt_*.csv"))
         if not any(k in f.name for k in ("summary", "forecast", "calibration", "cleanrepro"))]
n = w = fa = att = 0
units = set()
for f in files:
    d = pd.read_csv(f)
    if "wrong_output" not in d:
        continue
    geom = d["geom"].iloc[0] if "geom" in d and d["geom"].notna().any() else "main"
    n += len(d); w += int(d.wrong_output.sum()); fa += int(d.false_accepts.sum()); att += int(d.decode_attempts.sum())
    units |= {(geom, sc, int(r)) for sc, r in zip(d.scenario, d.rep)}
check("Abstract / VI-A", "retrievals", n, 55110, 0)
check("Abstract / VI-A", "wrong outputs", w, 0, 0)
check("Abstract / VI-A", "false accepts", fa, 0, 0)
check("Abstract / VI-A", "decode attempts (millions)", att / 1e6, 69.5, 0.05)
check("Abstract / VI-A", "synthesis replicates", len(units), 580, 0)
check("Abstract / VI-A", "wrong-output UL per retrieval (1e-5)", cp_upper(w, n) * 1e5, 6.7, 0.05)
check("Abstract / VI-A", "wrong-output UL per synthesis (1e-3)", cp_upper(w, len(units)) * 1e3, 6.3, 0.05)
check("Abstract / VI-A", "false-accept UL per attempt (1e-8)", cp_upper(fa, att) * 1e8, 5.3, 0.05)

# ---------------------------------------------------------------- baseline consumed vs reserved (Sec. VI-B, Table V)
b = load("dt_baseline.csv")
stat = pd.concat([x for x in (load("dt_baseline_static.csv"), load("dt_baseline_static_r4fine.csv"),
                              load("dt_baseline_static_r4cal.csv")) if x is not None])
allb = pd.concat([b, stat])
rr = allb[allb.policy == "round_robin_stop"]
m, lo, hi = cb(rr, lambda x: x.reads.mean())
check("VI-B", "round-robin mean reads (k)", m / 1e3, 3.49, 0.005)
check("VI-B", "round-robin mean CI lo (k)", lo / 1e3, 3.37, 0.02)
check("VI-B", "round-robin mean CI hi (k)", hi / 1e3, 3.60, 0.02)
check("VI-B", "uniform reservation q99 (k)", np.quantile(allb[allb.policy == "static_a0.0"].reads, .99) / 1e3, 5.15, 0.005)
check("VI-B", "best pre-specified reservation q99 (k)",
      min(np.quantile(allb[allb.policy == f"static_a{a}"].reads, .99) for a in ("0.0", "0.5", "1.0", "2.0")) / 1e3,
      4.80, 0.005)
check("VI-B", "stopping saving vs best reservation (%)",
      (1 - m / min(np.quantile(allb[allb.policy == f"static_a{a}"].reads, .99)
                   for a in ("0.0", "0.5", "1.0", "2.0"))) * 100, 27.3, 0.3)
check("VI-B", "stopping saving vs uniform reservation (%)",
      (1 - m / np.quantile(allb[allb.policy == "static_a0.0"].reads, .99)) * 100, 32.2, 0.3)
for pol, claim, tolv in (("static_a1.0", -4.21, 0.05), ("static_a1.5", -4.57, 0.05), ("static_cal_a1.0", -4.35, 0.05)):
    e, l_, h_ = paired(allb[allb.policy == pol], rr)
    check("VI-B", f"{pol} vs round-robin (%)", e * 100, claim, tolv)
check("VI-B", "static alpha=1 mean (k)", allb[allb.policy == "static_a1.0"].reads.mean() / 1e3, 3.343, 0.002)
check("VI-B", "calibration-estimated mean (k)", allb[allb.policy == "static_cal_a1.0"].reads.mean() / 1e3, 3.338, 0.002)

# ---------------------------------------------------------------- policy comparison (Sec. VI-C, Table VI)
for pol, claim in (("thompson", 7), ("tail_aware", 13), ("coverage", 15), ("ucb", 15),
                   ("learn_then_commit", 45), ("pivotal", 49), ("info_greedy", 441)):
    e, _, _ = paired(b[b.policy == pol], b[b.policy == "round_robin_stop"])
    check("VI-C", f"{pol} vs round-robin (%)", round(e * 100), claim, 0.5)
rpc = lambda g: np.stack(g.reads_per_container.map(json.loads).values)
sh = lambda g: (rpc(g) / rpc(g).sum(1, keepdims=True)).mean(0)[-1]
check("VI-C", "round-robin near-dead read share", sh(b[b.policy == "round_robin_stop"]), 0.16, 0.005)
check("VI-C", "info-greedy near-dead read share", sh(b[b.policy == "info_greedy"]), 0.66, 0.005)

# ---------------------------------------------------------------- failed container (Sec. VI-D, Table VII)
f1 = load("dt_fail1.csv")
f1s = load("dt_fail1_static.csv")
a1 = pd.concat([f1, f1s])
check("VI-D", "fail1 success (all policies)", a1.groupby("policy").apply(
    lambda g: (g.status == "exact").mean(), include_groups=False).unique()[0], 0.34, 0.005)
check("VI-D", "fail1 distinct success values", a1.groupby("policy").apply(
    lambda g: round((g.status == "exact").mean(), 4), include_groups=False).nunique(), 1, 0)
rr1 = f1[f1.policy == "round_robin_stop"]
check("VI-D", "fail1 round-robin mean on successes (k)", rr1[rr1.status == "exact"].reads.mean() / 1e3, 20.31, 0.01)
u1 = f1[f1.policy == "ucb"]
check("VI-D", "fail1 UCB mean on successes (k)", u1[u1.status == "exact"].reads.mean() / 1e3, 15.49, 0.01)
check("VI-D", "fail1 UCB vs round-robin, conditional (%)", paired(u1, rr1)[0] * 100, -24, 0.6)
check("VI-D", "fail1 UCB mean over all attempts (k)", u1.reads.mean() / 1e3, 52.79, 0.01)
check("VI-D", "fail1 round-robin mean over all attempts (k)", rr1.reads.mean() / 1e3, 54.43, 0.01)
check("VI-D", "fail1 UCB vs round-robin, unconditional (%)", (u1.reads.mean() / rr1.reads.mean() - 1) * 100, -3, 0.2)

# ---------------------------------------------------------------- ablation (Sec. VI-E, Table VIII)
abl = load("r4_ablation.csv")
g = abl[(abl.scenario == "baseline") & (abl.policy == "coverage")]
for rho, sd, claim in (("0.1", "0.5", 15.1), ("0.0", "0.5", -4.4), ("0.1", "1.0", -0.6)):
    v = g[(g.model == "theta") & (g.rho_U.astype(str) == rho) & (g.prior_sd.astype(str) == sd)]
    check("VI-E", f"coverage theta rho={rho} s={sd} (%)", v.rel_vs_rr.iloc[0] * 100, claim, 0.1)
v = g[(g.model == "theta+d") & (g.rho_U.astype(str) == "0.1") & (g.prior_sd.astype(str) == "0.5")]
check("VI-E", "coverage theta+d rho=0.1 s=0.5 (%)", v.rel_vs_rr.iloc[0] * 100, 1.6, 0.1)

# ---------------------------------------------------------------- dropout model (Sec. VI-F)
drop = load("dt_baseline_r4drop.csv")
rrb = b[b.policy == "round_robin_stop"]
e, l_, h_ = paired(drop[drop.policy == "ucb"], rrb)
check("VI-F", "UCB+dropout vs round-robin (%)", e * 100, -3.5, 0.1)
check("VI-F", "UCB+dropout CI lo (%)", l_ * 100, -4.9, 0.3)
check("VI-F", "UCB+dropout CI hi (%)", h_ * 100, -2.2, 0.3)
check("VI-F", "coverage+dropout near-dead share", sh(drop[drop.policy == "coverage"]), 0.20, 0.006)
check("VI-F", "UCB+dropout near-dead share", sh(drop[drop.policy == "ucb"]), 0.17, 0.006)
f1d = load("dt_fail1_r4drop.csv")
if f1d is not None:
    check("VI-F", "fail1 coverage+dropout vs round-robin (%)", paired(f1d[f1d.policy == "coverage"], rr1)[0] * 100,
          -35.0, 0.6)

# ---------------------------------------------------------------- tails (Sec. VI-H, Table IX)
ext = load("dt_baseline_r4ext.csv")
extd = load("dt_baseline_r4ext_drop.csv")
ste = load("dt_baseline_static_r4ext.csv")
big = pd.concat([b[b.perm < 3], ext, extd, ste])
rrx = big[big.policy == "round_robin_stop"]
check("VI-H", "100-synthesis replicate count", rrx.rep.nunique(), 100, 0)
check("VI-H", "100-synthesis round-robin q99 (k)", np.quantile(rrx.reads, .99) / 1e3, 5.15, 0.005)
stb = load("dt_baseline_static.csv")
sa1 = pd.concat([ste[ste.policy == "static_a1.0"],
                 stb[(stb.perm < 3) & (stb.policy == "static_a1.0")]])   # all 100 replicates
check("VI-H", "100-synthesis static alpha=1 q99 (k)", np.quantile(sa1.reads, .99) / 1e3, 4.70, 0.005)
e, l_, h_ = paired(sa1, rrx, stat=lambda x: np.quantile(x.reads, .99) - np.quantile(x.reads_r, .99), ok_only=False)
check("VI-H", "static alpha=1 paired q99 difference (reads)", e, -450, 30)
check("VI-H", "static q99 difference CI excludes zero (hi<0)", 1.0 if h_ < 0 else 0.0, 1.0, 0)
ud = drop[drop.policy == "ucb"]
e, l_, h_ = paired(ud, rrb, stat=lambda x: np.quantile(x.reads, .99) - np.quantile(x.reads_r, .99), ok_only=False)
check("VI-H", "UCB+dropout paired q99 difference (reads)", e, -450, 40)
check("VI-H", "UCB+dropout q99 CI contains zero (hi>0)", 1.0 if h_ > 0 else 0.0, 1.0, 0)
check("VI-H", "UCB+dropout q99 numerical change (%)",
      (np.quantile(ud.reads, .99) / np.quantile(rrb.reads, .99) - 1) * 100, -8.7, 0.15)

# ---------------------------------------------------------------- calibration (Sec. VI-I, Table X)
ce = load("dt_calibration_eval_baseline.csv")
sel = ce[np.isclose(ce.rho, 0.1)]
check("VI-I", "held-out coverage at rho=0.1", sel.covered.mean(), 0.87, 0.01)
check("VI-I", "near-dead coverage at rho=0.1", sel[sel.container == 5].covered.mean(), 0.20, 0.01)
r4c = load("r4_recalibration.csv")
v = r4c[r4c.split.str.startswith("validate")]
check("VI-I", "best re-selected coverage (s=0.5)", v.coverage.max(), 0.831, 0.005)
check("VI-I", "re-selected best rho", float(v.loc[v.coverage.idxmax(), "rho_U"]), 0.03, 0.001)
check("VI-I", "near-dead coverage max over extended grid", v.near_dead.max(), 0.033, 0.005)

# ---------------------------------------------------------------- forecast (Sec. VI-I)
fc = load("dt_forecast_baseline.csv")
check("VI-I", "forecast checkpoints", len(fc) // 2, 1127, 0)
check("VI-I", "forecast base rate", fc[fc.forecast == "mc"].completed.mean(), 0.49, 0.02)

# ---------------------------------------------------------------- replay (Sec. VI-J, Table XI)
tw, gs = load("replay_v3.csv"), load("replay_v3_genscript.csv")
for nm, d_, u_claim in (("Twist", tw, 13.9), ("Genscript", gs, 25.3)):
    uni = float(d_[d_.policy == "static_uniform"].mean_reads.iloc[0])
    rrr = float(d_[d_.policy == "uniform"].mean_reads.iloc[0])
    check("VI-J", f"{nm} stopping vs fixed depth (%)", (1 - rrr / uni) * 100, u_claim, 0.1)
check("VI-J", "Twist round-robin mean (k)", float(tw[tw.policy == "uniform"].mean_reads.iloc[0]) / 1e3, 58.6, 0.05)
check("VI-J", "Genscript round-robin mean (k)", float(gs[gs.policy == "uniform"].mean_reads.iloc[0]) / 1e3, 637, 0.5)
cal = load("r4_replay_calcost.csv")
full = cal[cal.belief == "full"].iloc[0]
check("VI-J", "calibrated replay mean (k)", full.mean_reads / 1e3, 533, 0.5)
check("VI-J", "calibration reads total (M)", full.calib_reads_total / 1e6, 1.8, 0.01)
check("VI-J", "break-even retrievals (full)", full.break_even_retrievals, 17.3, 0.2)
c100 = cal[cal.calib_reads_per_container == 100000].iloc[0]
check("VI-J", "break-even retrievals (100k/container)", c100.break_even_retrievals, 10.0, 0.2)
c30 = cal[cal.calib_reads_per_container == 30000].iloc[0]
check("VI-J", "30k/container is worse than generic (reads k)", c30.mean_reads / 1e3, 658, 1.0)

# ---------------------------------------------------------------- photolithographic (Sec. VI-C/VI-K)
ph = pd.concat([x for x in (load("dt_photolitho.csv"), load("dt_photolitho_static.csv"),
                            load("dt_photolitho_r4full.csv")) if x is not None])
check("VI-C", "photolithographic retrievals", len(ph), 660, 0)
check("VI-C", "photolithographic successes", int((ph.status == "exact").sum()), 0, 0)

# ---------------------------------------------------------------- consensus cap (Sec. VI-K, Table XIII)
cap = load("r4_consensus_cap.csv")
f1cap = cap[cap.scenario == "fail1"]
for c, claim in ((12, 0.313), (24, 0.340), (48, 0.357)):
    r_ = f1cap[f1cap.cap == c]
    if len(r_):
        check("VI-K", f"fail1 round-robin success at cap {c}", r_.success.iloc[0], claim, 0.002)

# ---------------------------------------------------------------- failure modes & geometries (Tables XIV, XV)
fm = load("r4_failure_modes.csv")
for setting, pol, claim in (("partial_c0_50", "coverage", 76.6), ("partial_c0_50", "coverage+dropout", -6.0),
                            ("common_cause", "coverage", 46.9), ("common_cause", "coverage+dropout", -33.6)):
    r_ = fm[(fm.setting == setting) & (fm.policy == pol)]
    check("VI-K", f"{setting} {pol} vs round-robin (%)", r_.rel_vs_rr_success.iloc[0] * 100, claim, 0.2)
ge = load("r4_geometries.csv")
cov = ge[ge.policy == "coverage"]
check("VI-K", "geometry min coverage penalty (%)", cov.rel_vs_rr_success.min() * 100, 1.1, 0.2)
check("VI-K", "geometry max coverage penalty (%)", cov.rel_vs_rr_success.max() * 100, 93.5, 0.3)

# ---------------------------------------------------------------- cost model (Table XVI)
cm = load("r4_cost_model.csv")
cmb = cm[(cm.scenario == "baseline") & cm.switches.notna()].drop_duplicates("policy")
ad = cmb[~cmb.policy.str.startswith("Static") & ~cmb.policy.str.contains("static")
          & (cmb.policy != "round_robin_stop")]
check("VI-L", "min adaptive container switches", ad.switches.min(), 38, 1)
check("VI-L", "max adaptive container switches", ad.switches.max(), 61, 1)
check("VI-L", "round-robin container switches", cmb[cmb.policy == "round_robin_stop"].switches.iloc[0], 63, 1)
check("VI-L", "max static container switches", cmb[cmb.policy.str.contains("static")].switches.max(), 69, 1)

# ---------------------------------------------------------------- what the PDF actually prints
txt = subprocess.run(["pdftotext", str(PDF), "-"], capture_output=True, text=True).stdout
present = {
    "abstract 3.49k": "3.49k" in txt,
    "abstract 27--32%": "27–32%" in txt,
    "no bit-for-bit overclaim": "bit-for-bit reproducible given the committed" in txt.replace("\n", " "),
    "no wet-lab claim": "No DNA was" in txt,
    "oracle decodability stated": "oracle" in txt,
    "near-dead 0.20 preserved": "0.20" in txt,
    "17 retrievals break-even": "17 retrievals" in txt,
    "per-synthesis bound present": "6.3" in txt,
    "no TBD markers": "TBD" not in txt,
    "no unresolved refs": "??" not in txt,
}
for k, v in present.items():
    rows.append((v, "PDF text", k, 1.0 if v else 0.0, 1.0, 0, ""))

bad = [r for r in rows if not r[0]]
lines = [f"Numeric verification of the manuscript against committed result files.",
         f"Checks: {len(rows)}, failures: {len(bad)}.", ""]
lines.append(f"{'':3s} {'section':14s} {'quantity':52s} {'recomputed':>12s} {'in paper':>10s}")
for ok, sec, what, got, claimed, tol, unit in rows:
    g = "" if got is None else (f"{got:12.4g}" if isinstance(got, (int, float, np.floating)) else str(got))
    c = "" if claimed is None else f"{claimed:10.4g}"
    lines.append(f"{'ok ' if ok else 'FAIL'} {sec:14s} {what:52s} {g} {c}")
lines += ["", "ALL CHECKS PASS" if not bad else f"{len(bad)} CHECK(S) FAILED"]
OUT.write_text("\n".join(lines) + "\n")
print("\n".join(lines))
sys.exit(1 if bad else 0)
