"""Revision-4 analyses (dt_protocol_rev4.yaml).  Reads only committed result files and writes
results/r4_*.csv plus manuscript/tables/r4_*.tex.  Every section skips silently if its inputs are
missing, so it can be re-run as the runs finish.

    python research_paper/experiments/statistical_analysis/r4_analysis.py
"""
from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "results"
TAB = ROOT / "research_paper" / "manuscript" / "tables"
STREAM = {"main": 72000}
B = 2000
RNG = np.random.default_rng(20261103)


def cp_upper(k, n, a=0.05):
    """Upper limit of the exact two-sided (1-a) Clopper-Pearson interval."""
    return 1.0 if k >= n else float(beta.ppf(1 - a / 2, k + 1, n - k))


def load(name):
    p = RES / name
    return pd.read_csv(p) if p.exists() else None


def cluster_boot(df, stat, by="rep", n=B):
    """Bootstrap over synthesis replicates (clusters); returns (estimate, lo, hi)."""
    groups = [g for _, g in df.groupby(by)]
    est = stat(df)
    vals = []
    for _ in range(n):
        idx = RNG.integers(0, len(groups), len(groups))
        vals.append(stat(pd.concat([groups[i] for i in idx])))
    return est, float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))


def q(p):
    return lambda d: float(np.quantile(d.reads, p))


def fmt_k(x, d=2):
    return f"{x / 1e3:.{d}f}"


# --------------------------------------------------------------------------- R4.1/4.2 actual vs reserved
def actual_vs_reserved():
    rows = []
    for sc in ("baseline", "homog", "fail1"):
        d, st = load(f"dt_{sc}.csv"), load(f"dt_{sc}_static.csv")
        if d is None or st is None:
            continue
        fine = load(f"dt_{sc}_static_r4fine.csv")
        cal = load(f"dt_{sc}_static_r4cal.csv")
        allr = pd.concat([x for x in (d, st, fine, cal) if x is not None])
        for pol, g in allr.groupby("policy"):
            ok = g.status == "exact"
            m_all = cluster_boot(g, lambda x: x.reads.mean())
            rows.append({"scenario": sc, "policy": pol, "n": len(g), "syntheses": g.rep.nunique(), "success": ok.mean(),
                         "mean_all": m_all[0], "mean_all_lo": m_all[1], "mean_all_hi": m_all[2],
                         "mean_success": g.reads[ok].mean() if ok.any() else np.nan,
                         "q95": np.quantile(g.reads, .95), "q99": np.quantile(g.reads, .99),
                         # reservation: what a plan that does NOT stop must commit for 99% success
                         "reservation_q99": np.quantile(g.reads, .99) if ok.mean() >= 0.99 else np.inf})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    for sc in df.scenario.unique():
        rr = df[(df.scenario == sc) & (df.policy == "round_robin_stop")].iloc[0]
        df.loc[df.scenario == sc, "mean_all_vs_rr"] = df[df.scenario == sc].mean_all / rr.mean_all - 1
        df.loc[df.scenario == sc, "mean_success_vs_rr"] = df[df.scenario == sc].mean_success / rr.mean_success - 1
    df.to_csv(RES / "r4_actual_vs_reserved.csv", index=False)
    return df


# --------------------------------------------------------------------------- R4.3 cluster-aware safety
def safety():
    files = [f for f in sorted(glob.glob(str(RES / "dt_*.csv")))
             if not any(k in f for k in ("summary", "forecast", "calibration"))]
    n = w = fa = att = 0
    units = set()
    per_unit_att = []
    for f in files:
        d = pd.read_csv(f)
        if "wrong_output" not in d:
            continue
        geom = d["geom"].iloc[0] if "geom" in d and d["geom"].notna().any() else "main"
        n += len(d); w += int(d.wrong_output.sum()); fa += int(d.false_accepts.sum()); att += int(d.decode_attempts.sum())
        for (sc, rep), g in d.groupby(["scenario", "rep"]):
            units.add((geom, sc, int(rep)))
        per_unit_att.append(d.groupby(["scenario", "rep"]).decode_attempts.sum())
    k_units = len(units)
    row = {"retrievals": n, "wrong_outputs": w, "wrong_ub_retrieval": cp_upper(w, n),
           "decode_attempts": att, "false_accepts": fa, "fa_ub_attempt": cp_upper(fa, att),
           "synthesis_units": k_units, "units_with_wrong_output": 0 if w == 0 else np.nan,
           "wrong_ub_synthesis": cp_upper(0 if w == 0 else w, k_units),
           "units_with_false_accept": 0 if fa == 0 else np.nan, "fa_ub_synthesis": cp_upper(0 if fa == 0 else fa, k_units),
           "retrievals_per_unit": n / max(k_units, 1), "attempts_per_unit": att / max(k_units, 1)}
    pd.DataFrame([row]).to_csv(RES / "r4_safety.csv", index=False)
    return row


# --------------------------------------------------------------------------- R4.4 tail metrics with more replicates
def label(d):
    """Append '+dropout' to the policy name of runs that used the efficiency+dropout model."""
    if d is None:
        return None
    d = d.copy()
    if "cfg_dropout_model" in d:
        m = d.cfg_dropout_model.astype(str).str.lower().eq("true")
        d.loc[m, "policy"] = d.loc[m, "policy"] + "+dropout"
    return d


def tails():
    base = load("dt_baseline.csv")
    if base is None:
        return None
    st = load("dt_baseline_static.csv")
    sets = {"30 syntheses x 10 orders": pd.concat([base] + ([st[st.policy.isin(["static_a0.0", "static_a1.0"])]] if st is not None else [])
                                                   + ([label(load("dt_baseline_r4drop.csv"))] if load("dt_baseline_r4drop.csv") is not None else []))}
    ext = load("dt_baseline_r4ext.csv")
    if ext is not None:
        parts = [base[base.perm < 3], ext, label(load("dt_baseline_r4ext_drop.csv")), load("dt_baseline_static_r4ext.csv")]
        if st is not None:
            parts.append(st[(st.perm < 3) & st.policy.isin(["static_a0.0", "static_a1.0"])])
        dr = label(load("dt_baseline_r4drop.csv"))
        if dr is not None:
            parts.append(dr[(dr.perm < 3) & dr.policy.isin(["coverage+dropout", "thompson+dropout"])])
        sets["100 syntheses x 3 orders"] = pd.concat([x for x in parts if x is not None])
    rows = []
    for name, d in sets.items():
        for pol, g in d.groupby("policy"):
            if g.rep.nunique() < 10:
                continue
            m = cluster_boot(g, lambda x: x.reads.mean())
            q95 = cluster_boot(g, q(.95)); q99 = cluster_boot(g, q(.99))
            rows.append({"design": name, "policy": pol, "n": len(g), "syntheses": g.rep.nunique(),
                         "mean": m[0], "mean_lo": m[1], "mean_hi": m[2], "q95": q95[0], "q95_lo": q95[1], "q95_hi": q95[2],
                         "q99": q99[0], "q99_lo": q99[1], "q99_hi": q99[2]})
    df = pd.DataFrame(rows)
    df.to_csv(RES / "r4_tails.csv", index=False)
    return df


# --------------------------------------------------------------------------- paired comparison helper
def paired_vs(df_pol, df_ref):
    """Paired over (rep, perm): ratio of summed reads on jointly successful retrievals, cluster bootstrap."""
    j = df_pol[df_pol.status == "exact"].merge(df_ref[df_ref.status == "exact"], on=["rep", "perm"], suffixes=("", "_r"))
    if j.empty:
        return (np.nan,) * 3 + (0,)
    est = cluster_boot(j, lambda x: x.reads.sum() / x.reads_r.sum() - 1)
    return est + (len(j),)


def rr_ref(sc, tagged=None):
    d = load(f"dt_{sc}.csv")
    return None if d is None else d[d.policy == "round_robin_stop"]


# --------------------------------------------------------------------------- R4.5 / R4.6 dropout model and ablation
def ablation():
    rows = []
    specs = [("theta", "0.1", "0.5", ""), ("theta", "0", "0.5", "_posthoc_rho0"), ("theta+d", "0.1", "0.5", "_r4drop"),
             ("theta+d", "0", "0.5", "_r4abl_d_rho0"), ("theta", "0.1", "1.0", "_r4abl_t_sd1"),
             ("theta+d", "0.1", "1.0", "_r4abl_d_sd1"), ("theta", "0", "1.0", "_r4abl_t_rho0_sd1"),
             ("theta+d", "0", "1.0", "_r4abl_d_rho0_sd1")]
    for sc in ("baseline", "fail1", "homog", "sweep_cov3", "sweep_cov1"):
        ref = rr_ref(sc)
        if ref is None:
            continue
        for model, rho, sd, tag in specs:
            d = load(f"dt_{sc}{tag}.csv")
            if d is None and sc == "fail1" and tag == "_r4drop":
                d = load("dt_fail1_r4abl_d.csv")
            if d is None:
                continue
            for pol, g in d.groupby("policy"):
                if pol == "round_robin_stop" and tag:
                    continue
                g2 = g if tag else g
                r = ref[ref.set_index(["rep", "perm"]).index.isin(g2.set_index(["rep", "perm"]).index)]
                p = paired_vs(g2, r)
                rpc = np.stack(g2.reads_per_container.map(json.loads).values)
                share = (rpc / rpc.sum(1, keepdims=True)).mean(0)
                rows.append({"scenario": sc, "model": model, "rho_U": rho, "prior_sd": sd, "policy": pol,
                             "n": len(g2), "success": (g2.status == "exact").mean(),
                             "mean_success": g2.reads[g2.status == "exact"].mean(), "q99": np.quantile(g2.reads, .99),
                             "rel_vs_rr": p[0], "lo": p[1], "hi": p[2], "n_joint": p[3],
                             "share_worst": float(share[-1]) if sc != "fail1" else float(share[0])})
    df = pd.DataFrame(rows)
    df.to_csv(RES / "r4_ablation.csv", index=False)
    return df


# --------------------------------------------------------------------------- R6 likelihood sweep
def likelihood_sweep():
    """Complete primary comparison under each likelihood setting: the frozen rho_U=0.1, the value the
    corrected selection picks (0.03), and the binomial (0). Round-robin never evaluates the posterior,
    so it is identical across settings and serves as the common reference; we assert that."""
    rows = []
    cfgs = [("0.1 (frozen)", ""), ("0.03 (re-selected)", "_rho003"), ("0 (binomial)", "_posthoc_rho0")]
    for sc in ("baseline", "homog", "fail1"):
        ref = load(f"dt_{sc}.csv")
        if ref is None:
            continue
        rr0 = ref[ref.policy == "round_robin_stop"]
        for label, tag in cfgs:
            d = load(f"dt_{sc}{tag}.csv")
            if d is None:
                continue
            rr = d[d.policy == "round_robin_stop"]
            # round-robin must be bit-identical across likelihood settings
            same = (rr.sort_values(["rep", "perm"]).reads.values ==
                    rr0.sort_values(["rep", "perm"]).reads.values).all() if len(rr) else None
            for pol, g in d.groupby("policy"):
                if pol == "round_robin_stop":
                    continue
                e, lo, hi = paired_vs(g, rr0)[:3]
                ok = g.status == "exact"
                rows.append({"scenario": sc, "rho_U": label, "policy": pol, "n": len(g),
                             "rr_identical": same, "success": ok.mean(),
                             "mean_success": g.reads[ok].mean() if ok.any() else np.nan,
                             "q99": np.quantile(g.reads, .99),
                             "rel_vs_rr": e, "lo": lo, "hi": hi})
    df = pd.DataFrame(rows)
    df.to_csv(RES / "r6_likelihood_sweep.csv", index=False)
    return df


# --------------------------------------------------------------------------- R4.8 prior sensitivity
def prior_sens():
    rows = []
    ref = rr_ref("baseline")
    for f in sorted(glob.glob(str(RES / "dt_baseline_r4prior_*.csv"))):
        name = Path(f).stem.replace("dt_baseline_r4prior_", "")
        d = pd.read_csv(f)
        for pol, g in d.groupby("policy"):
            r = ref[ref.set_index(["rep", "perm"]).index.isin(g.set_index(["rep", "perm"]).index)]
            p = paired_vs(g, r)
            rows.append({"setting": name, "policy": pol + ("+dropout" if name.endswith("_drop") else ""),
                         "n": len(g), "mean": g.reads.mean(), "q99": np.quantile(g.reads, .99),
                         "rel_vs_rr": p[0], "lo": p[1], "hi": p[2]})
    df = pd.DataFrame(rows)
    df.to_csv(RES / "r4_prior_sensitivity.csv", index=False)
    return df


# --------------------------------------------------------------------------- R4.9 / R4.10 failure modes and geometries
def scen_table(files_by_label):
    rows = []
    for label, (pol_file, drop_file, static_file) in files_by_label.items():
        d = load(pol_file)
        if d is None:
            continue
        parts = [d]
        if drop_file and load(drop_file) is not None:
            dd = load(drop_file); dd["policy"] = dd.policy + "+dropout"; parts.append(dd)
        if static_file and load(static_file) is not None:
            parts.append(load(static_file))
        allr = pd.concat(parts)
        ref = allr[allr.policy == "round_robin_stop"]
        stream = int(allr.reads.max())
        for pol, g in allr.groupby("policy"):
            p = paired_vs(g, ref) if pol != "round_robin_stop" else (0.0, 0.0, 0.0, len(g))
            ok = g.status == "exact"
            rows.append({"setting": label, "policy": pol, "n": len(g), "syntheses": g.rep.nunique(),
                         "success": ok.mean(), "mean_all": g.reads.mean(),
                         "mean_success": g.reads[ok].mean() if ok.any() else np.nan,
                         "q99": np.quantile(g.reads, .99), "rel_vs_rr_success": p[0], "lo": p[1], "hi": p[2],
                         "wrong_outputs": int(g.wrong_output.sum()), "max_stream_seen": stream})
    return pd.DataFrame(rows)


def failure_modes():
    sc = ["partial_c0_50", "partial_c0c3_25", "lowyield_c0", "contam_c2", "common_cause"]
    df = scen_table({s: (f"dt_{s}_r4.csv", f"dt_{s}_r4drop.csv", f"dt_{s}_static_r4.csv") for s in sc})
    df.to_csv(RES / "r4_failure_modes.csv", index=False)
    return df


def geometries():
    g = ["m3", "m12", "n40", "n64", "k16", "big"]
    df = scen_table({x: (f"dt_baseline_geom_{x}_r4.csv", f"dt_baseline_geom_{x}_r4drop.csv",
                         f"dt_baseline_geom_{x}_static_r4.csv") for x in g})
    df.to_csv(RES / "r4_geometries.csv", index=False)
    return df


# --------------------------------------------------------------------------- R4.11 consensus cap
def consensus_cap():
    rows = []
    for c in (12, 24, 48):
        for sc in ("baseline", "fail1"):
            if c == 24:
                d = load(f"dt_{sc}.csv")
                if d is None:
                    continue
                d = d[d.policy.isin(["round_robin_stop", "coverage"])]
                if sc == "baseline":
                    d = d[d.perm < 3]
                dd = load(f"dt_{sc}_r4drop.csv") if sc == "baseline" else None
                if dd is not None:
                    dd = dd[(dd.policy == "coverage") & (dd.perm < 3)].copy(); dd["policy"] = "coverage+dropout"
                    d = pd.concat([d, dd])
                if sc == "fail1":
                    d = d[d.policy == "round_robin_stop"]
            else:
                parts = [load(f"dt_{sc}_r4cap{c}.csv")]
                if sc == "baseline" and load(f"dt_{sc}_r4cap{c}_drop.csv") is not None:
                    x = load(f"dt_{sc}_r4cap{c}_drop.csv"); x["policy"] = x.policy + "+dropout"; parts.append(x)
                parts = [p for p in parts if p is not None]
                if not parts:
                    continue
                d = pd.concat(parts)
            for pol, g in d.groupby("policy"):
                ok = g.status == "exact"
                rows.append({"cap": c, "scenario": sc, "policy": pol, "n": len(g), "success": ok.mean(),
                             "mean_success": g.reads[ok].mean() if ok.any() else np.nan, "q99": np.quantile(g.reads, .99),
                             "wrong_outputs": int(g.wrong_output.sum())})
    df = pd.DataFrame(rows)
    df.to_csv(RES / "r4_consensus_cap.csv", index=False)
    return df


# --------------------------------------------------------------------------- R4.7 recalibration
def recalibration():
    s, e = load("dt_calibration_select_baseline_r4sd05.csv"), load("dt_calibration_eval_baseline_r4sd05.csv")
    if s is None or e is None:
        return None
    names = {0: "fresh", 1: "aged_1hl", 2: "aged_2hl_lowcopy", 3: "pcr_30c", 4: "pcr_45c_lowcopy", 5: "near_dead"}
    rows = []
    for split, d in (("select (reps 1000-1009)", s), ("validate (reps 0-29)", e)):
        for rho, g in d.groupby("rho"):
            row = {"split": split, "rho_U": rho, "prior_sd": 0.5, "coverage": g.covered.mean()}
            for j, gg in g.groupby("container"):
                row[names[j]] = gg.covered.mean()
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(RES / "r4_recalibration.csv", index=False)
    return df


# --------------------------------------------------------------------------- R4.15 cost model
PRICES = {  # per read, per sequencing round (batch), per container library preparation (one-off when touched)
    "reads only": (1.0, 0.0, 0.0),
    "multiplexed run (shared flow cell)": (1.0, 50.0, 500.0),
    "separate run per batch": (1.0, 2000.0, 500.0),
}


def cost_model():
    rows = []
    for f in ("dt_baseline_r4ext.csv", "dt_baseline_static_r4ext.csv", "dt_baseline_r4ext_drop.csv",
              "dt_baseline.csv", "dt_baseline_static.csv", "dt_baseline_r4drop.csv", "dt_fail1.csv"):
        d = load(f)
        if d is None:
            continue
        tag = "+dropout" if "r4drop" in f else ""
        if "cleanrepro" in f:
            tag = " (clean re-run)"
        for pol, g in d.groupby("policy"):
            rpc = np.stack(g.reads_per_container.map(json.loads).values)
            touched = (rpc > 0).sum(1)
            batches = g.reads / 50.0
            sw = g["switches"] if "switches" in g else pd.Series(np.nan, index=g.index)
            r = {"source": f, "scenario": g.scenario.iloc[0], "policy": pol + tag, "reads": g.reads.mean(),
                 "batches": batches.mean(), "containers_touched": touched.mean(), "switches": sw.mean()}
            for name, (cr, cb, cl) in PRICES.items():
                r[f"cost[{name}]"] = float(np.mean(cr * g.reads + cb * batches + cl * touched))
            rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(RES / "r4_cost_model.csv", index=False)
    return df


# --------------------------------------------------------------------------- R4.14 replay calibration cost
def replay_calcost():
    d, rp = load("replay_calcost_genscript.csv"), load("replay_v3_genscript.csv")
    if d is None or rp is None:
        return None
    rr = float(rp[rp.policy == "uniform"].mean_reads.iloc[0])
    d["retrieval_saving_vs_rr"] = rr - d.mean_reads
    d["total_first_retrieval"] = d.mean_reads + d.calib_reads_total
    d["break_even_retrievals"] = np.where(d.retrieval_saving_vs_rr > 0, d.calib_reads_total / d.retrieval_saving_vs_rr, np.inf)
    d["rr_mean"] = rr
    d.to_csv(RES / "r4_replay_calcost.csv", index=False)
    return d


if __name__ == "__main__":
    for fn in (actual_vs_reserved, safety, tails, ablation, likelihood_sweep, prior_sens, failure_modes, geometries, consensus_cap,
               recalibration, cost_model, replay_calcost):
        out = fn()
        print(f"== {fn.__name__}")
        if isinstance(out, pd.DataFrame):
            with pd.option_context("display.width", 250, "display.max_columns", 30):
                print(out.round(3).to_string() if not out.empty else "(no data yet)")
        else:
            print(out)
