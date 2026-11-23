"""Tables added in the final revision, built only from committed result files.

    python research_paper/experiments/scripts/make_r9_tables.py

  tables/r9_baselines.tex     stopping benchmarks and model-free allocation baselines (main text)
  tables/supp_infogreedy.tex  every information-greedy result, moved out of the main tables (supplement)
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "results"
TAB = ROOT / "research_paper" / "manuscript" / "tables"


def k(x, d=2):
    return "--" if not np.isfinite(x) else f"{x / 1e3:.{d}f}"


def pct(x, d=1):
    if not np.isfinite(x):
        return "--"
    v = round(x * 100, d)
    return "$0$" if v == 0 else f"${v:+.{d}f}$"


def baselines():
    res = pd.read_csv(RES / "r9_reservations.csv")
    ext = pd.read_csv(RES / "r9_extra_baselines.csv")
    rr_all = pd.read_csv(RES / "dt_baseline.csv")
    rr = rr_all[rr_all.policy == "round_robin_stop"]
    plan = {"static_a0.0": "fixed depth", "static_a1.0": r"$w\propto e$"}
    rule = {"q99": "q99 budget", "1.2 x median": r"$1.2\times$ median", "1.2 x mean": r"$1.2\times$ mean"}
    out = [r"\begin{table}[!t]", r"\centering",
           r"\caption{Stopping benchmarks and model-free allocation baselines (digital twin).  "
           r"\textbf{(a)} On the same 300 retrievals as Table~\ref{tab:actual}: verified stopping against plans "
           r"that do not stop and commit a fixed budget.  Success is the share of retrievals the budget completes "
           r"(cluster-bootstrap 95\% CI); reads per completed retrieval assume that a failed attempt is repeated "
           r"with a fresh budget (budget/success).  \textbf{(b)} Exploratory baselines on \emph{regenerated} read "
           r"streams (DT4DDS is not seed-deterministic), each paired against round-robin re-run on the same streams; "
           r"mean reads over successful retrievals, relative difference on jointly successful retrievals.  "
           r"Per-read stop applies the stopping test after every read instead of every 50.  No retrieval in (b) "
           r"returned a wrong file.  Thousands of reads.}",
           r"\label{tab:r9base}", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
           r"\begin{tabular}{@{}llrrr@{}}", r"\toprule",
           r"\multicolumn{5}{@{}l}{\textbf{(a) Committed retrievals, baseline (30 syntheses $\times$ 10 orders)}} \\",
           r"Plan & Budget rule & Reads (budget) & Success [95\% CI] & Reads per completed \\", r"\midrule",
           f"Round-robin + verified stop & stops when verified & {k(rr.reads.mean())} (mean) & 1.00 & "
           f"{k(rr.reads.mean())} \\\\"]
    for p in ("static_a0.0", "static_a1.0"):
        for r in res[res.plan == p].itertuples():
            out.append(f"Static, {plan[p]}, no stop & {rule[r.rule]} & {k(r.budget)} & "
                       f"{r.success:.2f} [{r.success_lo:.2f}, {r.success_hi:.2f}] & {k(r.reads_per_completed)} \\\\")
    out += [r"\midrule",
            r"\multicolumn{5}{@{}l}{\textbf{(b) Regenerated streams (30 syntheses $\times$ 10 orders per scenario)}} \\",
            r"Scenario & Policy & Success & Mean reads [95\% CI] & vs.\ round-robin, \% [95\% CI] \\", r"\midrule"]
    pn = {"round_robin_stop": "Round-robin + verified stop", "per_read_stop": "Round-robin, per-read stop",
          "succ_elim": "Successive elimination", "yield_threshold": "Yield threshold (25\\%)"}
    sn = {"baseline": "Baseline", "homog": "Homogeneous", "fail1": "One failed"}
    for s in ("baseline", "homog", "fail1"):
        d = ext[ext.scenario == s]
        first = True
        for p in pn:
            r = d[d.policy == p]
            if r.empty:
                continue
            r = r.iloc[0]
            rel = "" if p == "round_robin_stop" else f"{pct(r.rel_vs_rr)} [{pct(r.lo)}, {pct(r.hi)}]"
            out.append(f"{sn[s] if first else ''} & {pn[p]} & {r.success:.2f} & "
                       f"{k(r.mean_success)} [{k(r.mean_lo)}, {k(r.mean_hi)}] & {rel} \\\\")
            first = False
        out.append(r"\midrule")
    out[-1] = r"\bottomrule"
    out += [r"\end{tabular}", r"\end{table}"]
    (TAB / "r9_baselines.tex").write_text("\n".join(out) + "\n")
    print("wrote r9_baselines.tex")


def infogreedy():
    S = pd.read_csv(RES / "dt_summary.csv")
    lik = pd.read_csv(RES / "r6_likelihood_sweep.csv")
    abl = pd.read_csv(RES / "r4_ablation.csv")
    tails = pd.read_csv(RES / "r4_tails.csv")
    tw = pd.read_csv(RES / "replay_v3.csv")
    rows = []
    for s, sname in (("baseline", "Baseline"), ("homog", "Homogeneous"), ("fail1", "One failed")):
        r = S[(S.scenario == s) & (S.policy == "info_greedy")].iloc[0]
        mean = r.mean_reads if np.isfinite(r.mean_reads) else r.mean_success
        rows.append((sname, r"primary ($\theta$, $\rho_U{=}0.1$, $s{=}0.5$)", f"{r.success:.2f}", k(mean) +
                     ("" if np.isfinite(r.mean_reads) else "$^{s}$"),
                     f"{pct(r.rel_vs_rr, 0)} [{pct(r.rel_vs_rr_lo, 0)}, {pct(r.rel_vs_rr_hi, 0)}]"))
        l = lik[(lik.scenario == s) & (lik.policy == "info_greedy") & lik.rho_U.str.startswith("0.03")]
        if len(l):
            l = l.iloc[0]
            rows.append(("", r"$\rho_U{=}0.03$", f"{l.success:.2f}", k(l.mean_success),
                         f"{pct(l.rel_vs_rr, 0)} [{pct(l.lo, 0)}, {pct(l.hi, 0)}]"))
        a = abl[(abl.scenario == s) & (abl.policy == "info_greedy") & (abl.model == "theta+d")]
        if len(a):
            a = a.iloc[0]
            rows.append(("", r"$\theta{+}d$ (dropout model)", f"{a.success:.2f}", k(a.mean_success),
                         f"{pct(a.rel_vs_rr, 0)} [{pct(a.lo, 0)}, {pct(a.hi, 0)}]"))
    t = tails[(tails.design.str.startswith("30")) & (tails.policy == "info_greedy")].iloc[0]
    ti = tw[tw.policy == "info"].iloc[0]
    out = [r"\begin{table}[!ht]", r"\centering",
           r"\caption{Information-greedy (Bayesian experimental design) results, removed from the main tables.  "
           r"Paired relative difference in mean reads against round-robin with verified stopping on jointly "
           r"successful retrievals (cluster-bootstrap 95\% CI).  $^{s}$: mean over successful retrievals.  "
           r"Thousands of reads.}",
           r"\label{tab:suppinfo}", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
           r"\begin{tabular}{@{}llrrr@{}}", r"\toprule",
           r"Scenario & Configuration & Success & Mean reads & vs.\ round-robin, \% [95\% CI] \\", r"\midrule"]
    out += [" & ".join(r) + r" \\" for r in rows]
    out += [r"\bottomrule", r"\end{tabular}", "", r"\vspace{2pt}",
            r"\parbox{\columnwidth}{\scriptsize Baseline tail (frozen configuration): q99 "
            f"{k(t.q99)}k [{k(t.q99_lo)}, {k(t.q99_hi)}].  "
            r"Twist replay: after the first probes information-greedy reduced to coverage-greedy "
            f"(mean {ti.mean_reads / 1e3:.1f}k, q99 {ti.q99_reads / 1e3:.1f}k).  "
            r"Genscript replay: not completed; more than 1\% of its 200 trials never recovered the file, so the "
            r"condition is censored (Section~\ref{app:devfull}).}",
            r"\end{table}"]
    (TAB / "supp_infogreedy.tex").write_text("\n".join(out) + "\n")
    print("wrote supp_infogreedy.tex")


if __name__ == "__main__":
    baselines()
    infogreedy()
