"""Fill the @@TOKEN@@ placeholders of the manuscript text from committed result files.

    python research_paper/experiments/scripts/fill_text_numbers.py [--check]

Every number that the final revision added to the prose, and every interval the prose quotes from a table, is
written here from the same derived CSVs that build the tables, so text and tables cannot disagree.  Writes
manuscript/text_numbers.tex (a record of every value) and substitutes the tokens in main.tex and the
tables/*.tex files that carry them.  --check only reports tokens without a value.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "statistical_analysis"))
from bootstrap import cluster_boot  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "results"
MAN = ROOT / "research_paper" / "manuscript"


def sgn(x, d=1):
    """signed per cent with LaTeX minus, e.g. $-$4.9"""
    v = round(x * 100, d)
    if v == 0:
        v = 0.0
    s = f"{abs(v):.{d}f}"
    return ("$-$" if v < 0 else "$+$") + s


def ci(lo, hi, d=1):
    return f"{sgn(lo, d)} to {sgn(hi, d)}"


def tok():
    T = {}
    # --- Table V / round-robin
    act = pd.read_csv(RES / "r4_actual_vs_reserved.csv")
    r = act[(act.scenario == "baseline") & (act.policy == "round_robin_stop")].iloc[0]
    T["RR_CI"] = f"{r.mean_all_lo / 1e3:.2f}--{r.mean_all_hi / 1e3:.2f}"
    b0 = pd.read_csv(RES / "dt_baseline.csv")
    rr0 = b0[b0.policy == "round_robin_stop"]
    st = pd.concat([pd.read_csv(RES / f) for f in ("dt_baseline_static.csv", "dt_baseline_static_r4fine.csv",
                                                   "dt_baseline_static_r4cal.csv")])
    for key, pol in (("SA1_CI", "static_a1.0"), ("SA15_CI", "static_a1.5"), ("SCAL_CI", "static_cal_a1.0")):
        j = st[(st.policy == pol) & (st.status == "exact")].merge(rr0[rr0.status == "exact"], on=["rep", "perm"],
                                                                  suffixes=("", "_r"))
        e = cluster_boot(j, lambda x: x.reads.sum() / x.reads_r.sum() - 1)
        T[key] = ci(e[1], e[2])
    # --- reservations
    res = pd.read_csv(RES / "r9_reservations.csv")
    med = res[(res.plan == "static_a0.0") & (res.rule == "1.2 x median")].iloc[0]
    mea = res[(res.plan == "static_a0.0") & (res.rule == "1.2 x mean")].iloc[0]
    T["RES_MED_BUD"] = f"{med.budget / 1e3:.2f}"
    T["RES_MED_SUCC"] = f"{med.success * 100:.0f}"
    T["RES_MED_SUCC_CI"] = f"{med.success_lo * 100:.0f}--{med.success_hi * 100:.0f}"
    T["RES_MED_RPC"] = f"{med.reads_per_completed / 1e3:.2f}"
    T["RES_MEAN_BUD"] = f"{mea.budget / 1e3:.2f}"
    T["RES_MEAN_SUCC"] = f"{mea.success * 100:.0f}"
    T["RES_MEAN_RPC"] = f"{mea.reads_per_completed / 1e3:.2f}"
    # --- regenerated-stream baselines
    ex = pd.read_csv(RES / "r9_extra_baselines.csv")
    p = ex[(ex.scenario == "baseline") & (ex.policy == "per_read_stop")].iloc[0]
    T["PRS_PCT"] = f"{-p.rel_vs_rr * 100:.1f}"
    T["PRS_CI"] = ci(p.lo, p.hi)
    T["EXTRA_TEXT"] = extra_text(ex)
    # --- likelihood sweep
    lik = pd.read_csv(RES / "r6_likelihood_sweep.csv")
    b03 = lik[(lik.scenario == "baseline") & lik.rho_U.str.startswith("0.03")].set_index("policy")
    for key, pol in (("LIK_COV", "coverage"), ("LIK_TS", "thompson"), ("LIK_UCB", "ucb"), ("LIK_TAIL", "tail_aware")):
        r = b03.loc[pol]
        T[key] = f"{sgn(r.rel_vs_rr)}\\% ({ci(r.lo, r.hi)}\\%)"
    T["LIK_PIV_PT"] = sgn(b03.loc["pivotal"].rel_vs_rr)
    T["LIK_LTC_PT"] = sgn(b03.loc["learn_then_commit"].rel_vs_rr)
    hm = lik[(lik.scenario == "homog") & (lik.policy != "info_greedy")]
    T["HOMOG_RANGE"] = f"between {sgn(hm.rel_vs_rr.min())}\\% and {sgn(hm.rel_vs_rr.max())}\\%"
    # --- Table VI (dt_summary): baseline and one failed container
    S = pd.read_csv(RES / "dt_summary.csv").set_index(["scenario", "policy"])

    def pos(r):   # "7\% (4--11\%)" for a penalty
        return f"{r.rel_vs_rr * 100:.0f}\\% ({r.rel_vs_rr_lo * 100:.0f}--{r.rel_vs_rr_hi * 100:.0f}\\%)"

    def neg(r):   # "24\% (17--28\%)" for a saving
        return f"{-r.rel_vs_rr * 100:.0f}\\% ({-r.rel_vs_rr_hi * 100:.0f}--{-r.rel_vs_rr_lo * 100:.0f}\\%)"

    for key, pol in (("DT_TS", "thompson"), ("DT_TAIL", "tail_aware"), ("DT_COV", "coverage"), ("DT_UCB", "ucb"),
                     ("DT_LTC", "learn_then_commit"), ("DT_PIV", "pivotal")):
        T[key] = pos(S.loc[("baseline", pol)])
    for key, pol in (("F1_UCB", "ucb"), ("F1_COV", "coverage"), ("F1_PIV", "pivotal"), ("F1_TS", "thompson")):
        T[key] = neg(S.loc[("fail1", pol)])
    r = S.loc[("fail1", "tail_aware")]
    T["F1_TAIL"] = f"{sgn(r.rel_vs_rr, 0)}\\%, {ci(r.rel_vs_rr_lo, r.rel_vs_rr_hi, 0)}\\%"
    T["F1_LTC"] = pos(S.loc[("fail1", "learn_then_commit")])
    # --- ablation
    ab = pd.read_csv(RES / "r4_ablation.csv")
    ab = ab[ab.scenario == "baseline"]

    def cell(model, rho, sd, pol):
        r = ab[(ab.model == model) & np.isclose(ab.rho_U, rho) & np.isclose(ab.prior_sd, sd) & (ab.policy == pol)].iloc[0]
        return r

    for key, args in (("ABL_COV_PRIM", ("theta", .1, .5)), ("ABL_COV_RHO0", ("theta", 0, .5)),
                      ("ABL_COV_SD1", ("theta", .1, 1.0)), ("ABL_COV_DROP", ("theta+d", .1, .5))):
        r = cell(*args, "coverage")
        T[key] = f"{sgn(r.rel_vs_rr)}\\% (95\\% CI {ci(r.lo, r.hi)}\\%)"
    r = cell("theta+d", .1, .5, "ucb")
    T["UCBD_CI"] = ci(r.lo, r.hi)
    # --- tails
    tl = pd.read_csv(RES / "r4_tails.csv")

    def dq(design, pol, style="paren"):
        r = tl[tl.design.str.startswith(design) & (tl.policy == pol)].iloc[0]
        f = lambda x: f"${x:+,.0f}$".replace(",", "{,}")
        body = f"{f(r.dq99_vs_rr)} reads"
        c = f"95\\% CI {f(r.dq99_lo)} to {f(r.dq99_hi)}"
        return (f"{body} ({c})" if style == "paren" else f"{body}, {c}"), r

    T["TQ_SA1_30"], r30 = dq("30", "static_a1.0")
    T["TQ_SA1_30_PCT"] = sgn(r30.dq99_vs_rr / tl[tl.design.str.startswith("30") & (tl.policy == "round_robin_stop")].q99.iloc[0])
    T["TQ_SA1_100"], r100 = dq("100", "static_a1.0")
    T["TQ_SA1_100_PCT"] = sgn(r100.dq99_vs_rr / tl[tl.design.str.startswith("100") & (tl.policy == "round_robin_stop")].q99.iloc[0])
    T["TQ_UCBD_30"] = dq("30", "ucb+dropout", "comma")[0]
    T["TQ_COVD_30"] = dq("30", "coverage+dropout", "comma")[0]
    T["TQ_COVD_100"] = dq("100", "coverage+dropout", "comma")[0]
    T["TQ_TSD_100"] = dq("100", "thompson+dropout", "comma")[0]
    # --- container switches on the 300 retrievals
    cm = pd.read_csv(RES / "r4_cost_model.csv")
    cm = cm[cm.scenario == "baseline"]
    T["SW_RR"] = f"{cm[(cm.source == 'dt_baseline_rho003.csv') & (cm.policy == 'round_robin_stop@rho0.03')].switches.iloc[0]:.0f}"
    d = cm[(cm.source == "dt_baseline_r4drop.csv") & ~cm.policy.str.startswith("info_greedy")].switches
    T["SW_DROP"] = f"{d.min():.0f}--{d.max():.0f}"
    # --- Section VII-A/B
    b = pd.read_csv(RES / "dt_baseline.csv")
    rr = b[b.policy == "round_robin_stop"]
    T["CV_MAIN"] = f"{rr.reads.std() / rr.reads.mean():.2f}"
    T["GAP_MAIN"] = f"{(np.quantile(rr.reads, .99) / rr.reads.mean() - 1) * 100:.0f}"
    g = rr.groupby("rep").reads
    bv, wv = g.mean().var(ddof=1), g.var(ddof=1).mean()
    T["BETWEEN_SHARE"] = f"{bv / (bv + wv) * 100:.0f}"
    T["WITHIN_SHARE"] = f"{wv / (bv + wv) * 100:.0f}"
    e = np.array(json.load(open(RES / "dt_baseline_hindsight_eff.json")))
    T["EFF_RANGE"] = f"{e.min():.3f}--{e.max():.3f}"
    T["FLUID_BOUND"] = f"{(1 - e.mean() / e.max()) * 100:.0f}"
    for key, geom in (("BIG", "big"), ("N64", "n64"), ("M12", "m12")):
        dd = pd.read_csv(RES / f"dt_baseline_geom_{geom}_r4.csv")
        ss = pd.read_csv(RES / f"dt_baseline_geom_{geom}_static_r4.csv")
        r_ = dd[dd.policy == "round_robin_stop"].reads
        q = np.quantile(ss[ss.policy == "static_a0.0"].reads, .99)
        T[f"{key}_SAVING"] = f"{(1 - r_.mean() / q) * 100:.0f}"
        if key == "BIG":
            T["BIG_CV"] = f"{r_.std() / r_.mean():.2f}"
    return T


def extra_text(ex):
    def g(s, p):
        r = ex[(ex.scenario == s) & (ex.policy == p)]
        return None if r.empty else r.iloc[0]

    def eff(r):
        if abs(r.rel_vs_rr) < 1e-12 and abs(r.lo) < 1e-12 and abs(r.hi) < 1e-12:
            return "reproduced round-robin exactly"
        return f"gave {sgn(r.rel_vs_rr)}\\% ({ci(r.lo, r.hi)}\\%)"

    rrb = g("baseline", "round_robin_stop")
    parts = [f"On the regenerated baseline streams round-robin used {rrb.mean_success / 1e3:.2f}k reads, against "
             f"3.49k on the committed streams.  On the baseline, successive elimination never eliminated a container "
             f"and {eff(g('baseline', 'succ_elim'))}, and the yield threshold {eff(g('baseline', 'yield_threshold'))}: "
             f"the near-dead container's early yield is too close to the others' for either rule to set it aside."]
    se, yt = g("homog", "succ_elim"), g("homog", "yield_threshold")
    parts.append(f"In the homogeneous control successive elimination {eff(se)} and the yield threshold {eff(yt)}.")
    se, yt, rr1 = g("fail1", "succ_elim"), g("fail1", "yield_threshold"), g("fail1", "round_robin_stop")
    parts.append(f"With one failed container, successive elimination {eff(se)} and the yield threshold {eff(yt)} on "
                 f"the retrievals that succeeded.  On these regenerated streams {rr1.success * 100:.0f}\\% of "
                 f"one-failed retrievals succeeded, against 34\\% on the committed streams, a reminder that "
                 f"recoverability after a failure depends on which strands each synthesis loses.")
    w = int(ex.wrong_outputs.sum())
    parts.append((f"None of the {int(ex.n.sum()):,} retrievals on regenerated streams returned a wrong file."
                  .replace(",", "{,}")) if w == 0 else f"{w} retrievals returned a wrong file.")
    return "  ".join(parts)


if __name__ == "__main__":
    T = tok()
    files = [MAN / "main.tex", MAN / "supplement.tex"] + sorted((MAN / "tables").glob("*.tex"))
    missing = set()
    for f in files:
        s = f.read_text()
        for k in set(re.findall(r"@@([A-Z0-9_]+)@@", s)):
            if k not in T:
                missing.add(k)
        if "--check" not in sys.argv:
            s2 = re.sub(r"@@([A-Z0-9_]+)@@", lambda m: T.get(m.group(1), m.group(0)), s)
            if s2 != s:
                f.write_text(s2)
                print("filled", f.name)
    (MAN / "text_numbers.tex").write_text("".join(f"% {k} = {v}\n" for k, v in sorted(T.items())))
    print("missing:", sorted(missing) or "none")
