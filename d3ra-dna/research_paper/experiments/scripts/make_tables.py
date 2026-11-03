"""Generate the LaTeX tables of the revised manuscript from committed result files
(no hand-entered numbers).  Run dt_analysis.py first.

    python research_paper/experiments/scripts/make_tables.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

PROJ = Path(__file__).resolve().parents[3]
RES = PROJ / "results"
TAB = PROJ / "research_paper" / "manuscript" / "tables"
NAMES = {"round_robin_stop": "Round-robin + verified stop", "coverage": "Coverage-greedy", "pivotal": "Pivotal",
         "thompson": "Thompson sampling", "ucb": "UCB ($q{=}0.9$)", "learn_then_commit": "Learn-then-commit",
         "info_greedy": "Info-greedy (exp.\\ design)", "tail_aware": "Tail-aware (LCB)",
         "static_a0.0": "Static, uniform (fixed depth)", "static_a0.5": "Static, $w\\propto e^{0.5}$",
         "static_a1.0": "Static, $w\\propto e$ (equal coverage)", "static_a2.0": "Static, $w\\propto e^{2}$"}
SCEN = {"baseline": "Baseline", "fail1": "One container failed", "fail2": "Two containers failed",
        "homog": "Homogeneous (control)", "photolitho": "Photolithographic synthesis$^{\\ddagger}$"}


def cp_upper(k, n, a=0.05):
    return float(beta.ppf(1 - a / 2, k + 1, n - k)) if n > 0 and k < n else 1.0


def k(x, d=2):
    return "--" if not np.isfinite(x) else f"{x / 1e3:.{d}f}"


def pct(x):
    return "--" if not np.isfinite(x) else f"{100 * x:+.0f}"


def dt_main(S):
    lines = [r"\begin{table*}[!t]", r"\centering",
             r"\caption{Digital-twin end-to-end retrieval (DT4DDS; primary evidence).  Per scenario and policy: share of retrievals whose output file is byte-identical to the stored file (no retrieval returned a wrong file), mean reads per retrieval (thousands; cluster-bootstrap 95\% CI over syntheses), the 0.99-quantile q99 (the budget to reserve), and the paired relative difference in mean reads against round-robin with verified stopping on retrievals where both succeeded (95\% CI).  Static plans read fixed proportions until decodable and must reserve their q99 on every retrieval, so their mean cost equals q99; their proportions use per-container efficiencies measured in hindsight on the same scenario.}",
             r"\label{tab:dt}", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
             r"\begin{tabular}{@{}llrrrr@{}}", r"\toprule",
             r"Scenario & Policy & Success & Mean reads (k) [95\% CI] & q99 (k) & vs.\ round-robin, \% [95\% CI] \\", r"\midrule"]
    for s, sname in SCEN.items():
        d = S[S.scenario == s]
        if d.empty:
            continue
        if (d.success == 0).all():
            lines.append(f"{sname} & all {len(d)} policies and static plans & 0.00 & -- (stream exhausted, 72k reads) & -- & \\\\")
            lines.append(r"\midrule")
            continue
        first = True
        order = [p for p in NAMES if p in set(d.policy)]
        for p in order:
            r = d[d.policy == p].iloc[0]
            mean = (f"{k(r.mean_reads)} [{k(r.mean_lo)}, {k(r.mean_hi)}]" if np.isfinite(r.mean_reads)
                    else f"{k(r.mean_success)}$^{{s}}$")
            if p.startswith("static_"):
                mean = f"{k(r.mean_reads)} (= q99)" if np.isfinite(r.mean_reads) else "--"
            rel = ""
            if "rel_vs_rr" in r and pd.notna(r.get("rel_vs_rr", np.nan)):
                rel = f"{pct(r.rel_vs_rr)} [{pct(r.rel_vs_rr_lo)}, {pct(r.rel_vs_rr_hi)}]"
            lines.append(f"{sname if first else ''} & {NAMES[p]} & {r.success:.2f} & {mean} & {k(r.q99_reads)} & {rel} \\\\")
            first = False
        lines.append(r"\midrule")
    lines[-1] = r"\bottomrule"
    lines += [r"\end{tabular}", "", r"\vspace{2pt}",
              r"\parbox{\textwidth}{\scriptsize 30 independent syntheses of all six containers $\times$ 10 read orders = 300 retrievals per policy and scenario, except $^{\ddagger}$10 $\times$ 3 (protocol deviation, compute).  A failed retrieval read the whole stream (72{,}000 reads) without verifying; $^{s}$: mean over successful retrievals only (a mean or quantile including failures is undefined, ``--'').  When success $<1$, success was identical for every policy on every synthesis and read order, so the paired comparison uses the same retrievals for all policies.}",
              r"\end{table*}"]
    (TAB / "dt_main.tex").write_text("\n".join(lines) + "\n")
    print("wrote dt_main.tex")


def safety():
    rows = []
    files = sorted(RES.glob("dt_*.csv"))
    tot_n = tot_w = tot_fa = tot_att = tot_ok = 0
    for f in files:
        if any(t in f.name for t in ("calibration", "forecast", "summary")):
            continue
        d = pd.read_csv(f)
        if "wrong_output" not in d:
            continue
        tot_n += len(d)
        tot_ok += int(d.file_identical.sum())
        tot_w += int(d.wrong_output.sum())
        tot_fa += int(d.false_accepts.sum())
        tot_att += int(d.decode_attempts.sum())
        rows.append((f.name, len(d)))
    body = [r"\begin{table}[!t]", r"\centering",
            r"\caption{Safety of verification-gated stopping over every digital-twin retrieval run for this paper (all scenarios, policies, static plans and ablations; the forecast-evaluation runs are excluded).  Exact 95\% upper bounds (Clopper--Pearson).}",
            r"\label{tab:safety}", r"\footnotesize", r"\begin{tabular}{@{}lr@{}}", r"\toprule",
            f"End-to-end retrievals & {tot_n:,} \\\\".replace(",", "{,}"),
            f"\\quad output byte-identical to stored file & {tot_ok:,} \\\\".replace(",", "{,}"),
            f"\\quad reported as not recovered (no output) & {tot_n - tot_ok - tot_w:,} \\\\".replace(",", "{,}"),
            f"\\quad wrong output (file returned, not identical) & {tot_w} \\\\",
            f"Wrong-output rate, 95\\% upper bound & {cp_upper(tot_w, tot_n):.1e} \\\\",
            f"Fragment decode attempts & {tot_att:,} \\\\".replace(",", "{,}"),
            f"\\quad false accepts (wrong payload passed inner checks) & {tot_fa} \\\\",
            f"False-accept rate per attempt, 95\\% upper bound & {cp_upper(tot_fa, tot_att):.1e} \\\\",
            r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    (TAB / "safety.tex").write_text("\n".join(body) + "\n")
    print("wrote safety.tex", tot_n, tot_w, tot_fa, tot_att, [r for r in rows])
    return tot_n, tot_w, tot_fa, tot_att


def replay():
    out = [r"\begin{table}[!t]", r"\centering",
           r"\caption{Secondary evidence: semi-empirical replay of real Illumina read streams (ENA PRJEB65931), generic prior, no evaluation-split statistics.  Mean reads (thousands) and q99; $\Delta$: mean minus the hindsight-tuned static budget (paired bootstrap 95\% CI).  Decodability is an oracle event here (consensus equals the designed strand), so this is \emph{not} exact retrieval.}",
           r"\label{tab:replay}", r"\scriptsize", r"\setlength{\tabcolsep}{2.5pt}", r"\begin{tabular}{@{}llrrr@{}}", r"\toprule",
           r"Pool & Policy & Mean & q99 & $\Delta$ vs.\ static [95\% CI] \\", r"\midrule"]
    pn = {"static_uniform": "Static, uniform (hindsight budget)", "static_hindsight": "Static, hindsight-tuned",
          "uniform": "Round-robin + verified stop", "coverage": "Coverage-greedy", "pivotal": "Pivotal",
          "ts": "Thompson sampling", "ucb": "UCB", "ltc": "Learn-then-commit", "info": "Info-greedy",
          "tail_lcb": "Tail-aware (LCB)"}
    for pool, f in (("Twist", "replay_v3.csv"), ("Genscript", "replay_v3_genscript.csv")):
        if not (RES / f).exists():
            continue
        d = pd.read_csv(RES / f)
        first = True
        for _, r in d.iterrows():
            name = pn.get(r.policy, r.policy)
            if r.belief == "calibrated_leakfree":
                name = "Coverage-greedy, leak-free calibrated"
            if r.policy.startswith("static"):
                m, q, dl = r["cap@0.99"], r["cap@0.99"], ""
            else:
                m, q = r.mean_reads, r.q99_reads
                dl = f"{r.diff_vs_hindsight_static / 1e3:+.1f} [{r.diff_lo / 1e3:+.1f}, {r.diff_hi / 1e3:+.1f}]"
            dig = 1 if pool == "Twist" else 0
            out.append(f"{pool if first else ''} & {name} & {m / 1e3:.{dig}f} & {q / 1e3:.{dig}f} & {dl} \\\\")
            first = False
        out.append(r"\midrule")
    out[-1] = r"\bottomrule"
    out += [r"\end{tabular}", r"\end{table}"]
    (TAB / "replay.tex").write_text("\n".join(out) + "\n")
    print("wrote replay.tex")


if __name__ == "__main__":
    S = pd.read_csv(RES / "dt_summary.csv")
    dt_main(S)
    safety()
    replay()
