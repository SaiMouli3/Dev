"""Build the revision-4 LaTeX tables and macros from results/r4_*.csv (written by
research_paper/experiments/statistical_analysis/r4_analysis.py).  Each table is skipped if its
input file is missing, so this can be re-run while runs are still finishing.

    python research_paper/experiments/scripts/make_r4_tables.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "results"
TAB = ROOT / "research_paper" / "manuscript" / "tables"

PN = {"round_robin_stop": "Round-robin + verified stop", "coverage": "Coverage-greedy", "pivotal": "Pivotal",
      "thompson": "Thompson sampling", "ucb": "UCB", "learn_then_commit": "Learn-then-commit",
      "info_greedy": "Info-greedy", "tail_aware": "Tail-aware", "static_a0.0": "Static $\\alpha{=}0$ (fixed depth)",
      "static_a0.5": "Static $\\alpha{=}0.5$", "static_a1.0": "Static $\\alpha{=}1$", "static_a2.0": "Static $\\alpha{=}2$",
      "static_a0.25": "Static $\\alpha{=}0.25$", "static_a0.75": "Static $\\alpha{=}0.75$",
      "static_a1.25": "Static $\\alpha{=}1.25$", "static_a1.5": "Static $\\alpha{=}1.5$",
      "static_a3.0": "Static $\\alpha{=}3$", "static_cal_a0.5": "Static $\\alpha{=}0.5$, calibration-estimated",
      "static_cal_a1.0": "Static $\\alpha{=}1$, calibration-estimated",
      "static_cal_a2.0": "Static $\\alpha{=}2$, calibration-estimated"}


def name(p):
    base = p.replace("+dropout", "")
    return PN.get(base, base.replace("_", " ")) + ("\\,+\\,dropout model" if p.endswith("+dropout") else "")


def load(n):
    p = RES / n
    if not p.exists() or p.stat().st_size == 0:
        return None
    try:
        d = pd.read_csv(p)
    except pd.errors.EmptyDataError:
        return None
    return None if d.empty else d


def k(x, d=2):
    return "--" if not np.isfinite(x) else f"{x / 1e3:.{d}f}"


def pct(x, d=0):
    return "--" if not np.isfinite(x) else f"${x * 100:+.{d}f}$"


def write(fn, lines, msg):
    (TAB / fn).write_text("\n".join(lines) + "\n")
    print("wrote", fn, msg)


# ----------------------------------------------------------------------------- R4.1 actual vs reserved
def actual_vs_reserved():
    d = load("r4_actual_vs_reserved.csv")
    if d is None:
        return
    order = ["static_a0.0", "static_a0.25", "static_a0.5", "static_a0.75", "static_a1.0", "static_a1.25",
             "static_a1.5", "static_a2.0", "static_a3.0", "static_cal_a0.5", "static_cal_a1.0", "static_cal_a2.0",
             "round_robin_stop", "coverage", "pivotal", "thompson", "ucb", "learn_then_commit", "info_greedy",
             "tail_aware"]
    out = [r"\begin{table*}[!t]", r"\centering",
           r"\caption{Reads \emph{consumed} with verified stopping versus capacity that must be \emph{reserved} "
           r"without it (digital twin, 300 retrievals per row).  All rows stop on verification, so the comparison "
           r"between rows is allocation only; the last column is the budget the same read proportions would need "
           r"without stopping, which coincides with q99 wherever every retrieval succeeds.  "
           r"Thousands of reads; mean over all attempts with a cluster-bootstrap 95\% interval.}",
           r"\label{tab:actual}", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
           r"\begin{tabular}{@{}llrrrrr@{}}", r"\toprule",
           r"Scen. & Allocation & $p_{\text{ok}}$ & $\bar r_{\text{all}}$ [95\% CI] & q95 & q99 & $B_{0.99}$ \\",
           r"\midrule"]
    for sc, lab in (("baseline", "Base"), ("homog", "Homog.")):
        g = d[d.scenario == sc]
        if g.empty:
            continue
        first = True
        for pol in order:
            r = g[g.policy == pol]
            if r.empty:
                continue
            r = r.iloc[0]
            out.append(f"{lab if first else ''} & {name(pol)} & {r.success:.2f} & {k(r.mean_all)} "
                       f"[{k(r.mean_all_lo)}, {k(r.mean_all_hi)}] & {k(r.q95)} & {k(r.q99)} & {k(r.reservation_q99)} \\\\")
            first = False
        out.append(r"\midrule")
    out[-1] = r"\bottomrule"
    out += [r"\end{tabular}", r"\end{table*}"]
    write("r4_actual.tex", out, "")


# ----------------------------------------------------------------------------- R4.2 unconditional fail1
def unconditional():
    d = load("r4_actual_vs_reserved.csv")
    if d is None or "fail1" not in set(d.scenario):
        return
    g = d[d.scenario == "fail1"].copy()
    g["o"] = [0 if p == "round_robin_stop" else 1 for p in g.policy]
    g = g.sort_values(["o", "mean_success"])
    out = [r"\begin{table}[!t]", r"\centering",
           r"\caption{One failed container: cost conditional on success and over all attempts (300 retrievals per "
           r"row, 72k-read stream).  Every policy has the same success probability, so the conditional saving does "
           r"not buy a better chance of recovery, and failures dominate the unconditional mean.  Thousands of reads.}",
           r"\label{tab:uncond}", r"\scriptsize", r"\begin{tabular}{@{}lrrrr@{}}", r"\toprule",
           r"Allocation & $p_{\text{ok}}$ & $\bar r_{\text{ok}}$ & vs.\ RR & $\bar r_{\text{all}}$ \\", r"\midrule"]
    for r in g.itertuples():
        out.append(f"{name(r.policy)} & {r.success:.2f} & {k(r.mean_success)} & "
                   f"{pct(r.mean_success_vs_rr) if r.policy != 'round_robin_stop' else '--'} & {k(r.mean_all)} \\\\")
    out += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    write("r4_uncond.tex", out, "")


# ----------------------------------------------------------------------------- R4.3 safety
def safety():
    d = load("r4_safety.csv")
    if d is None:
        return
    r = d.iloc[0]

    def sci(x):
        e = int(np.floor(np.log10(x))); return f"${x / 10 ** e:.1f}\\times10^{{{e}}}$"
    out = [r"\begin{table}[!t]", r"\centering",
           r"\caption{Safety of verification-gated stopping over every digital-twin retrieval of this paper, with "
           r"the dependence between retrievals made explicit.  Retrievals that share a synthesis replicate are not "
           r"independent, so we also report the bound that treats one synthesis as one unit.  Upper limits of exact "
           r"two-sided 95\% Clopper--Pearson intervals.}",
           r"\label{tab:safety}", r"\scriptsize", r"\begin{tabular}{@{}lr@{}}", r"\toprule",
           r"Quantity & Value \\", r"\midrule",
           f"Retrievals (all scenarios, policies, geometries) & {int(r.retrievals):,} \\\\".replace(",", "{,}"),
           f"\\quad wrong output files & {int(r.wrong_outputs)} \\\\",
           f"Independent synthesis replicates & {int(r.synthesis_units)} \\\\",
           f"\\quad retrievals per replicate (mean) & {r.retrievals_per_unit:.0f} \\\\",
           f"\\quad synthesis replicates with a wrong output & {int(r.units_with_wrong_output)} \\\\",
           f"Wrong-output rate per retrieval, 95\\% upper limit & {sci(r.wrong_ub_retrieval)} \\\\",
           f"\\quad per synthesis replicate (conservative) & {sci(r.wrong_ub_synthesis)} \\\\",
           r"\midrule",
           f"Fragment decode attempts & {int(r.decode_attempts):,} \\\\".replace(",", "{,}"),
           f"\\quad false accepts (wrong payload passed inner checks) & {int(r.false_accepts)} \\\\",
           f"\\quad attempts per synthesis replicate (mean) & {r.attempts_per_unit:,.0f} \\\\".replace(",", "{,}"),
           f"False-accept rate per attempt, 95\\% upper limit & {sci(r.fa_ub_attempt)} \\\\",
           f"\\quad per synthesis replicate (conservative) & {sci(r.fa_ub_synthesis)} \\\\",
           r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    write("safety.tex", out, "(cluster-aware)")
    m = [r"% generated by make_r4_tables.py from results/r4_safety.csv",
         f"\\newcommand{{\\SafeN}}{{{int(r.retrievals):,}}}".replace(",", "{,}"),
         r"\newcommand{\SafeWrong}{" + str(int(r.wrong_outputs)) + "}",
         r"\newcommand{\SafeWUB}{" + sci(r.wrong_ub_retrieval) + "}",
         r"\newcommand{\SafeWUBc}{" + sci(r.wrong_ub_synthesis) + "}",
         r"\newcommand{\SafeUnits}{" + f"{int(r.synthesis_units)}" + "}",
         r"\newcommand{\SafePerUnit}{" + f"{r.retrievals_per_unit:.0f}" + "}",
         r"\newcommand{\SafeAtt}{" + f"{r.decode_attempts / 1e6:.1f} million" + "}",
         r"\newcommand{\SafeFA}{" + str(int(r.false_accepts)) + "}",
         r"\newcommand{\SafeFUB}{" + sci(r.fa_ub_attempt) + "}",
         r"\newcommand{\SafeFUBc}{" + sci(r.fa_ub_synthesis) + "}"]
    write("macros.tex", m, "")


# ----------------------------------------------------------------------------- R4.4 tails
def tails():
    d = load("r4_tails.csv")
    if d is None:
        return
    out = [r"\begin{table*}[!t]", r"\centering",
           r"\caption{Tail metrics with cluster-bootstrap 95\% intervals over synthesis replicates (baseline).  "
           r"Extreme quantiles estimated from a few hundred retrievals are imprecise, and the intervals of most "
           r"policies overlap; the ordering of the means does not carry over to q99.  Thousands of reads.}",
           r"\label{tab:r4tails}", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
           r"\begin{tabular}{@{}llrrr@{}}", r"\toprule",
           r"Design & Allocation & $\bar r_{\text{all}}$ & q95 [95\% CI] & q99 [95\% CI] \\", r"\midrule"]
    for des, g in d.groupby("design", sort=False):
        first = True
        for r in g.sort_values("mean").itertuples():
            out.append(f"{des.replace('x', '$\\times$') if first else ''} & {name(r.policy)} & {k(r.mean)} & "
                       f"{k(r.q95)} [{k(r.q95_lo)}, {k(r.q95_hi)}] & {k(r.q99)} [{k(r.q99_lo)}, {k(r.q99_hi)}] \\\\")
            first = False
        out.append(r"\midrule")
    out[-1] = r"\bottomrule"
    out += [r"\end{tabular}", r"\end{table*}"]
    write("r4_tails.tex", out, "")


# ----------------------------------------------------------------------------- R4.6 ablation
def ablation():
    """Pivoted: one row per scenario and policy, one column per ablation cell.  The long form with
    cluster-bootstrap intervals is results/r4_ablation.csv; intervals for the cells the text relies on
    are quoted in the text."""
    d = load("r4_ablation.csv")
    if d is None:
        return
    cells = [("theta", "0.1", "0.5"), ("theta", "0.0", "0.5"), ("theta", "0.1", "1.0"), ("theta", "0.0", "1.0"),
             ("theta+d", "0.1", "0.5"), ("theta+d", "0.0", "0.5"), ("theta+d", "0.1", "1.0"),
             ("theta+d", "0.0", "1.0")]
    d["key"] = list(zip(d.model, d.rho_U.astype(str), d.prior_sd.astype(str)))
    pols = ["coverage", "thompson", "tail_aware", "pivotal", "ucb", "learn_then_commit", "info_greedy"]
    out = [r"\begin{table*}[!t]", r"\centering",
           r"\caption{Ablation of the three modelling choices (digital twin; paired mean reads against "
           r"round-robin on jointly successful retrievals, per cent, negative is cheaper than round-robin).  "
           r"$\theta$: one efficiency multiplier per container (the frozen model); $\theta{+}d$: efficiency and "
           r"a dropout fraction.  $\rho_U$ is the beta-binomial correlation and $s$ the prior width; the frozen "
           r"configuration is the first column.  The penalty is large only in that column: relaxing any one of "
           r"the three choices moves every rule into a band a few per cent below round-robin.  Blank: not run.  "
           r"Cluster-bootstrap intervals for every cell are in \texttt{results/r4\_ablation.csv}; those the text "
           r"relies on are quoted there.}",
           r"\label{tab:r4abl}", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
           r"\begin{tabular}{@{}ll" + "r" * len(cells) + r"@{}}", r"\toprule",
           r" & & \multicolumn{4}{c}{$\theta$ only} & \multicolumn{4}{c}{$\theta$ and dropout $d$} \\",
           r"\cmidrule(lr){3-6}\cmidrule(lr){7-10}",
           r"Scenario & Policy ($\rho_U$, $s$) & " + " & ".join(
               f"${c[1]}$, ${c[2]}$" for c in cells) + r" \\", r"\midrule"]
    for sc in ["baseline", "fail1", "homog", "sweep_cov3", "sweep_cov1"]:
        g = d[d.scenario == sc]
        if g.empty:
            continue
        first = True
        for pol in pols:
            gg = g[g.policy == pol]
            if gg.empty:
                continue
            vals = []
            for c in cells:
                r_ = gg[gg.key == c]
                vals.append(f"{r_.rel_vs_rr.iloc[0] * 100:+.0f}" if len(r_) else "")
            out.append(f"{sc.replace('_', '-') if first else ''} & {name(pol)} & " + " & ".join(vals) + r" \\")
            first = False
        out.append(r"\midrule")
    out[-1] = r"\bottomrule"
    out += [r"\end{tabular}", r"\end{table*}"]
    write("r4_ablation.tex", out, "")


# ----------------------------------------------------------------------------- R4.9/4.10 settings table
def settings_table(src, fn, label, caption):
    d = load(src)
    if d is None or "setting" not in d:
        return
    out = [r"\begin{table*}[!t]", r"\centering", r"\caption{" + caption + "}", r"\label{" + label + "}",
           r"\scriptsize", r"\setlength{\tabcolsep}{3pt}", r"\begin{tabular}{@{}llrrrrr@{}}", r"\toprule",
           r"Setting & Allocation & $n$ & $p_{\text{ok}}$ & $\bar r_{\text{ok}}$ & vs.\ RR & q99 \\", r"\midrule"]
    for st, g in d.groupby("setting", sort=False):
        first = True
        g = g.copy()
        g["o"] = [0 if p == "round_robin_stop" else (2 if p.startswith("static") else 1) for p in g.policy]
        for r in g.sort_values(["o", "policy"]).itertuples():
            out.append(f"{st.replace('_', '-') if first else ''} & {name(r.policy)} & {r.n} & {r.success:.2f} & "
                       f"{k(r.mean_success)} & {pct(r.rel_vs_rr_success) if r.policy != 'round_robin_stop' else '--'} & "
                       f"{k(r.q99)} \\\\")
            first = False
        out.append(r"\midrule")
    out[-1] = r"\bottomrule"
    out += [r"\end{tabular}", r"\end{table*}"]
    write(fn, out, "")


# ----------------------------------------------------------------------------- R4.11 consensus cap
def cap_table():
    d = load("r4_consensus_cap.csv")
    if d is None:
        return
    out = [r"\begin{table}[!t]", r"\centering",
           r"\caption{Sensitivity to the consensus cap, the number of reads per cluster used to form a "
           r"consensus (frozen value 24).  On the baseline the cap is not binding: clusters rarely reach a "
           r"dozen reads before the fragment decodes.  Thousands of reads.}",
           r"\label{tab:r4cap}", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
           r"\begin{tabular}{@{}llrrrr@{}}", r"\toprule",
           r"Scen. & Allocation & Cap & $p_{\text{ok}}$ & $\bar r_{\text{ok}}$ & q99 \\", r"\midrule"]
    for (sc, pol), g in d.groupby(["scenario", "policy"], sort=False):
        first = True
        for r in g.sort_values("cap").itertuples():
            out.append(f"{sc if first else ''} & {name(pol) if first else ''} & {r.cap} & {r.success:.2f} & "
                       f"{k(r.mean_success)} & {k(r.q99)} \\\\")
            first = False
        out.append(r"\midrule")
    out[-1] = r"\bottomrule"
    out += [r"\end{tabular}", r"\end{table}"]
    write("r4_cap.tex", out, "")


# ----------------------------------------------------------------------------- R6 likelihood sweep
def likelihood_table():
    d = load("r6_likelihood_sweep.csv")
    if d is None:
        return
    order = ["0.1 (frozen)", "0.03 (re-selected)", "0 (binomial)"]
    pols = ["coverage", "thompson", "ucb", "tail_aware", "pivotal", "learn_then_commit", "info_greedy"]
    out = [r"\begin{table*}[!t]", r"\centering",
           r"\caption{The complete primary comparison repeated under each likelihood setting: the frozen "
           r"$\rho_U{=}0.1$, the value the corrected selection picks ($\rho_U{=}0.03$, Section~\ref{sec:res-calib}) "
           r"and the binomial ($\rho_U{=}0$).  All three use the full design (30 syntheses $\times$ 10 read "
           r"orders per policy) and the frozen prior width $s{=}0.5$.  Entries are the paired difference in mean "
           r"reads against round-robin on jointly successful retrievals, per cent, with cluster-bootstrap 95\% "
           r"intervals; negative is cheaper than round-robin.  Round-robin never evaluates the posterior, so it is "
           r"bit-identical across the three settings and is the common reference.  Success probability is "
           r"unchanged by $\rho_U$ in every scenario.}",
           r"\label{tab:r6lik}", r"\scriptsize", r"\setlength{\tabcolsep}{4pt}",
           r"\begin{tabular}{@{}ll" + "r" * len(order) + r"@{}}", r"\toprule",
           r"Scenario & Policy & " + " & ".join(
               {"0.1 (frozen)": r"$\rho_U{=}0.1$ (frozen)",
                "0.03 (re-selected)": r"$\rho_U{=}0.03$ (re-selected)",
                "0 (binomial)": r"$\rho_U{=}0$ (binomial)"}[o] for o in order) + r" \\", r"\midrule"]
    for sc in ["baseline", "homog", "fail1"]:
        g = d[d.scenario == sc]
        if g.empty:
            continue
        first = True
        for pol in pols:
            cells = []
            for o in order:
                r_ = g[(g.policy == pol) & (g.rho_U == o)]
                cells.append(f"{r_.rel_vs_rr.iloc[0] * 100:+.0f} [{r_.lo.iloc[0] * 100:+.0f}, "
                             f"{r_.hi.iloc[0] * 100:+.0f}]" if len(r_) else "")
            if not any(cells):
                continue
            out.append(f"{sc if first else ''} & {name(pol)} & " + " & ".join(cells) + r" \\")
            first = False
        out.append(r"\midrule")
    out[-1] = r"\bottomrule"
    out += [r"\end{tabular}", r"\end{table*}"]
    write("r6_likelihood.tex", out, "")


# ----------------------------------------------------------------------------- R4.7 recalibration
def recalibration():
    d = load("r4_recalibration.csv")
    if d is None:
        return
    cols = ["fresh", "aged_1hl", "aged_2hl_lowcopy", "pcr_30c", "pcr_45c_lowcopy", "near_dead"]
    short = {"fresh": "fresh", "aged_1hl": "aged 1\\,hl", "aged_2hl_lowcopy": "aged 2\\,hl, 10 cop.",
             "pcr_30c": "30 cyc.", "pcr_45c_lowcopy": "45 cyc., 20 cop.", "near_dead": "near-dead"}
    hdr = " & ".join(short[c] for c in cols)
    out = [r"\begin{table*}[!t]", r"\centering",
           r"\caption{Re-run of the $\rho_U$ selection with the prior actually used in the policy runs "
           r"($s{=}0.5$) and an extended grid: coverage of the 90\% credible interval of $\theta^\ast_j$, overall "
           r"and per container.  No value in the grid reaches 0.90 overall, and no value covers the near-dead "
           r"container.}",
           r"\label{tab:r4calib}", r"\scriptsize", r"\setlength{\tabcolsep}{2.5pt}",
           r"\begin{tabular}{@{}llr" + "r" * len(cols) + r"@{}}", r"\toprule",
           r"Split & $\rho_U$ & All & " + hdr + r" \\", r"\midrule"]
    for sp, g in d.groupby("split", sort=False):
        first = True
        for r in g.itertuples():
            vals = " & ".join(f"{getattr(r, c, float('nan')):.2f}" for c in cols)
            out.append(f"{sp if first else ''} & {r.rho_U} & {r.coverage:.2f} & {vals} \\\\")
            first = False
        out.append(r"\midrule")
    out[-1] = r"\bottomrule"
    out += [r"\end{tabular}", r"\end{table*}"]
    write("r4_calib.tex", out, "")


# ----------------------------------------------------------------------------- R4.14 calibration cost
def calcost():
    d = load("r4_replay_calcost.csv")
    if d is None:
        return
    out = [r"\begin{table}[!t]", r"\centering",
           r"\caption{What the calibrated prior costs (Genscript replay, coverage-greedy).  Calibration reads are "
           r"charged explicitly: they are the first reads of the same six sequencing runs, disjoint from the "
           r"evaluation reads.  Break-even is the number of retrievals over which the calibration reads are repaid "
           r"by the per-retrieval saving against round-robin.  Thousands of reads.}",
           r"\label{tab:r4calcost}", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
           r"\begin{tabular}{@{}lrrrrr@{}}", r"\toprule",
           r"Belief & Calib./cont. & Calib.\ total & Retrieval & Saving vs.\ RR & Break-even \\", r"\midrule"]
    for r in d.itertuples():
        be = "--" if not np.isfinite(r.break_even_retrievals) else f"{r.break_even_retrievals:.1f}"
        out.append(f"{r.belief} & {k(r.calib_reads_per_container)} & {k(r.calib_reads_total)} & {k(r.mean_reads)} & "
                   f"{k(r.retrieval_saving_vs_rr)} & {be} \\\\")
    out += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    write("r4_calcost.tex", out, "")


# ----------------------------------------------------------------------------- R4.15 cost model
def cost_model():
    d = load("r4_cost_model.csv")
    if d is None:
        return
    d = d[d.scenario == "baseline"].copy()
    # prefer the runs that record container switches (added in revision 4)
    if d.switches.notna().any():
        d = d[d.switches.notna()]
    d = d.drop_duplicates("policy")
    cols = [c for c in d.columns if c.startswith("cost[") and "reads only" not in c]
    out = [r"\begin{table}[!t]", r"\centering",
           r"\caption{Reads are not the only cost.  Two further charging models applied to the same baseline "
           r"retrievals; the reads-only model is the Reads column itself.  A multiplexed run charges 50 "
           r"read-equivalents per batch of 50 reads and 500 per container library preparation; "
           r"and a separate sequencing run per batch (2{,}000 per batch).  Costs in read-equivalents (thousands).  "
           r"Switching containers every batch is free in the first model and dominant in the third.}",
           r"\label{tab:r4cost}", r"\scriptsize", r"\setlength{\tabcolsep}{3pt}",
           r"\begin{tabular}{@{}p{26mm}rrr" + "r" * len(cols) + r"@{}}", r"\toprule",
           r"Allocation & Reads & Batch. & Sw. & " +
           " & ".join({"reads only": "reads", "multiplexed run (shared flow cell)": "multiplexed",
                       "separate run per batch": "run/batch"}.get(c[5:-1], c[5:-1]) for c in cols) + r" \\",
           r"\midrule"]
    for _, r in d.sort_values("reads").iterrows():
        sw = "--" if not np.isfinite(r.switches) else f"{r.switches:.0f}"
        vals = " & ".join(k(r[c]) for c in cols)
        out.append(f"{name(r.policy)} & {k(r.reads)} & {r.batches:.0f} & {sw} & {vals} \\\\")
    out += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    write("r4_cost.tex", out, "")


if __name__ == "__main__":
    actual_vs_reserved()
    unconditional()
    safety()
    tails()
    ablation()
    settings_table("r4_failure_modes.csv", "r4_failures.tex", "tab:r4fail",
                   r"Failure modes beyond total loss (20 syntheses $\times$ 3 read orders).  "
                   r"\texttt{partial}: a fraction of a container's designs is physically lost; \texttt{lowyield}: "
                   r"90\% of a container's reads do not derive from stored strands; \texttt{contam}: 50\% do not; "
                   r"\texttt{common-cause}: every container aged two further half-lives and amplified harder.  "
                   r"Paired against round-robin on jointly successful retrievals.  Thousands of reads.")
    settings_table("r4_geometries.csv", "r4_geoms.tex", "tab:r4geom",
                   r"Other archive geometries, baseline container conditions (20 syntheses $\times$ 3 read orders; "
                   r"\texttt{big}: 10 $\times$ 3).  \texttt{m3}/\texttt{m12}: 3 and 12 containers; \texttt{n40}/"
                   r"\texttt{n64}: $(32,40)$ and $(32,64)$; \texttt{k16}: $(16,24)$; \texttt{big}: a four times "
                   r"larger file.  Thousands of reads.")
    cap_table()
    recalibration()
    likelihood_table()
    calcost()
    cost_model()
