"""Real-data validation of D3RA-DNA on public DNA-storage sequencing data.

Data (downloaded by --download; not committed):
  * ENA PRJEB65931 (Gimpel et al. 2023, DT4DDS): Illumina iSeq reads of one
    12,000-strand Twist pool after different handling (aging 0/2/7 days,
    PCR 15/60/90 cycles); designed references from the DT4DDS release
    (Zenodo 12666365, scripts/internal_validation/Twist_GCall/design_files.fasta).
  * Zenodo 14296588: binned (clustered) real nanopore reads with ground truth
    (Srinivasavaradhan et al.).

Stages:
  prep      map reads, estimate per-strand empirical consensus curves (cached)
  channel   error profiles, coverage stats, shared-variance decomposition
  nanopore  real nanopore phi(n) vs simulated
  policies  real recovery-time matrices: placement x policy comparison, H1 on real data
"""
from __future__ import annotations

import argparse
import json
import pickle
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from d3ra.io.fastq import read_fasta  # noqa: E402
from d3ra.realdata import empirical as EM  # noqa: E402
from d3ra.realdata.channel_estimation import coverage_stats, error_profile, shared_variance_decomposition  # noqa: E402
from d3ra.realdata.mapping import KmerMapper  # noqa: E402

RAW = ROOT / "data" / "raw"
CACHE = ROOT / "data" / "processed" / "realdata"
RES = ROOT / "results"
REF_FASTA = RAW / "zenodo/fml-ethz-dt4dds-21968d9/scripts/internal_validation/Twist_GCall/design_files.fasta"
RUNS = {  # ENA run -> condition
    "ERR12033805": "Aging_0d", "ERR12033813": "Aging_2d", "ERR12033821": "Aging_7d",
    "ERR12033826": "PCR_15c", "ERR12033837": "PCR_60c", "ERR12033841": "PCR_90c",
}
N_CALIB, N_EVAL = 300_000, 900_000


def download():
    (RAW / "ena").mkdir(parents=True, exist_ok=True)
    (RAW / "zenodo").mkdir(parents=True, exist_ok=True)
    for run in RUNS:
        out = RAW / "ena" / f"{run}_R1.fq.gz"
        if not out.exists():
            url = f"https://ftp.sra.ebi.ac.uk/vol1/run/{run[:6]}/{run}/R1.fq.gz"
            subprocess.run(["curl", "-sSL", "-o", str(out), url], check=True)
    z = RAW / "zenodo"
    for rec, f in [("12666365", "fml-ethz/dt4dds-v1.1.0.zip"), ("14296588", "Srinivasavaradhan.txt")]:
        out = z / Path(f).name
        if not out.exists():
            subprocess.run(["curl", "-sSL", "-o", str(out), f"https://zenodo.org/records/{rec}/files/{f}?download=1"], check=True)
    if not REF_FASTA.exists():
        subprocess.run(["unzip", "-qo", str(z / "dt4dds-v1.1.0.zip"),
                        "fml-ethz-dt4dds-21968d9/scripts/internal_validation/*", "-d", str(z)], check=True)


def _prep_one(run):
    out = CACHE / f"{run}.pkl"
    if out.exists():
        return run, "cached"
    _, refs = read_fasta(REF_FASTA)
    mapper = KmerMapper(refs)
    t = time.time()
    rd = EM.load_run(RUNS[run], RAW / "ena" / f"{run}_R1.fq.gz", mapper, N_CALIB, N_EVAL, len(refs[0]))
    t_map = time.time() - t
    EM.estimate_phi(rd, refs, seed=hash(run) % 2**32)
    prof = error_profile(rd.pairs_sample)
    rd.reads_by_ref = {}
    rd.pairs_sample = []
    with open(out, "wb") as fh:
        pickle.dump({"run": rd, "error_profile": prof, "t_map": t_map, "t_total": time.time() - t}, fh)
    return run, f"{time.time() - t:.0f}s"


def prep(workers=4):
    CACHE.mkdir(parents=True, exist_ok=True)
    with Pool(workers) as p:
        for run, msg in p.imap_unordered(_prep_one, list(RUNS)):
            print("prep", run, RUNS[run], msg, flush=True)


def load_cached():
    out = {}
    for run in RUNS:
        with open(CACHE / f"{run}.pkl", "rb") as fh:
            out[run] = pickle.load(fh)
    return out


def nanopore_diagnostics(n=20, max_clusters=500, seed=1):
    """Why does real consensus fail where i.i.d. simulation succeeds?
    (a) transplant real edit scripts from OTHER strands' reads (realistic
        within-read structure, independent across reads);
    (b) disjoint-subset test: two consensuses from non-overlapping read sets of
        the same strand; identical wrong answers => strand-systematic error."""
    import edlib
    from d3ra.channel.empirical_errors import ErrorTransplant
    from d3ra.channel.sequencer import noisy_reads
    from d3ra.realdata.references import read_binned
    from d3ra.retrieval.consensus import consensus
    rng = np.random.default_rng(seed)
    clusters = read_binned(RAW / "zenodo" / "Srinivasavaradhan.txt")
    big = [c for c in clusters if len(c[1]) >= 2 * n]
    rng.shuffle(big)
    donors = [c for c in clusters if 10 <= len(c[1]) < 2 * n][:800]
    tp = ErrorTransplant([(ref, rd) for ref, reads in donors for rd in reads[:5]])
    big = big[:max_clusters]
    real_ok = tp_ok = both_wrong = same_wrong = 0
    wrong_dists = []
    for ref, reads in big:
        L = len(ref)
        idx = rng.permutation(len(reads))
        a = consensus([reads[i] for i in idx[:n]], L)
        b = consensus([reads[i] for i in idx[n:2 * n]], L)
        real_ok += a == ref
        if a != ref and b != ref:
            both_wrong += 1
            same_wrong += a == b
            wrong_dists.append(edlib.align(a, ref)["editDistance"])
        simr = noisy_reads(ref, n, tp, rng)
        tp_ok += consensus(simr, L) == ref
    m = len(big)
    out = {"n_reads": n, "clusters": m, "real_exact": real_ok / m, "transplant_independent_exact": tp_ok / m,
           "transplant_mean_errors_per_read": tp.mean_errors, "frac_both_disjoint_wrong": both_wrong / m,
           "frac_identical_wrong_given_both_wrong": same_wrong / max(both_wrong, 1),
           "median_edit_distance_of_wrong_consensus": float(np.median(wrong_dists)) if wrong_dists else 0.0}
    (RES / "realdata_nanopore_diagnostics.json").write_text(json.dumps(out, indent=1))
    print(out, flush=True)
    import pickle
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(CACHE / "transplant_nanopore.pkl", "wb") as fh:
        pickle.dump(tp, fh)
    return out


def transplant(n_reads=40000):
    """Build the Illumina error-transplant model from real aligned reads (fresh pool)."""
    import pickle
    from d3ra.channel.empirical_errors import ErrorTransplant
    from d3ra.io.fastq import read_fastq
    from d3ra.realdata.mapping import trim_adapter
    _, refs = read_fasta(REF_FASTA)
    mp = KmerMapper(refs)
    pairs = []
    for rec in read_fastq(RAW / "ena" / "ERR12033805_R1.fq.gz", limit=n_reads):
        s = trim_adapter(rec.seq)
        m = mp.map(s)
        if m.ref >= 0:
            pairs.append((refs[m.ref], s[: len(refs[m.ref]) + 8]))
    tp = ErrorTransplant(pairs)
    CACHE.mkdir(parents=True, exist_ok=True)
    with open(CACHE / "transplant_illumina.pkl", "wb") as fh:
        pickle.dump(tp, fh)
    print({"illumina_transplant_scripts": len(tp.scripts), "mean_errors_per_read": tp.mean_errors,
           "frac_error_free": tp.frac_clean}, flush=True)
    return tp


def channel():
    """Real channel parameters per run + cross-run shared-variance decomposition."""
    _, refs = read_fasta(REF_FASTA)
    C = load_cached()
    rows, prof_out = [], {}
    for run, d in C.items():
        rd, prof = d["run"], d["error_profile"]
        cov = coverage_stats(rd.calib_counts, refs)
        eval_counts = np.bincount(rd.ref_of_read[rd.ref_of_read >= 0], minlength=len(refs))
        phi = rd.phi_strand.mean(0)
        rows.append({"run": run, "condition": RUNS[run], "yield_mapped": rd.calib_yield,
                     **{k: v for k, v in prof.items() if k != "position_profile"}, **cov,
                     "dropout_at_full_depth": float(((rd.calib_counts + eval_counts) == 0).mean()),
                     "mean_coverage_full": float((rd.calib_counts + eval_counts).mean()),
                     **{f"phi_exact_n{n}": float(v) for n, v in zip(EM.N_GRID, phi)},
                     **{f"phi_h2_n{n}": float(v) for n, v in zip(EM.N_GRID, rd.phi_strand_h2.mean(0))},
                     "strands_unrecoverable_at_16": float((rd.phi_strand[:, -1] == 0).mean()),
                     "map_seconds": d["t_map"], "prep_seconds": d["t_total"]})
        prof_out[RUNS[run]] = prof["position_profile"]
        print(RUNS[run], {k: round(v, 4) if isinstance(v, float) else v for k, v in rows[-1].items()
                          if k in ("yield_mapped", "rate_sub", "rate_ins", "rate_del", "cv", "dropout_at_full_depth",
                                   "phi_exact_n1", "phi_exact_n3", "phi_exact_n5")}, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(RES / "realdata_channel.csv", index=False)
    M = np.stack([C[r]["run"].calib_counts for r in RUNS], 1)
    dec = shared_variance_decomposition(M)
    dec_fresh = shared_variance_decomposition(M[:, [0, 1, 3]])   # aging 0/2 d + PCR 15c (mild)
    out = {"all_runs": dec, "mild_runs": dec_fresh, "position_profiles": prof_out}
    (RES / "realdata_shared_variance.json").write_text(json.dumps(out, indent=1, default=float))
    print("shared-variance:", {k: round(v, 4) for k, v in dec.items()}, flush=True)
    return df, dec


def belief_from_calibration(C, refs):
    """ChannelParams estimated ONLY from calibration-split reads (+ pooled phi from
    even-indexed strands of the evaluation split; T uses all strands)."""
    from d3ra.channel.failures import FailureModel
    from d3ra.channel.surrogate import ChannelParams
    runs = [C[r]["run"] for r in RUNS]
    M = np.stack([rd.calib_counts for rd in runs], 1).astype(float)
    dec = shared_variance_decomposition(M)
    yields, shapes, drops, phis = [], [], [], []
    for rd in runs:
        c = rd.calib_counts.astype(float)
        mu, var = c.mean(), c.var()
        shapes.append(float(np.clip(mu ** 2 / max(var - mu, 1e-9), 0.5, 200)))
        yields.append(rd.calib_yield)
        drops.append(float(max((c == 0).mean() - np.mean(np.exp(-c.mean())), 1e-4)))
        phis.append(EM.grid_to_full(rd.phi_strand[::2].mean(0)))
    m = len(runs)
    return ChannelParams(np.array(yields), np.array(shapes), np.array(drops), np.ones(m), np.stack(phis),
                         FailureModel(np.full(m, 1e-3)), sigma_g=float(np.sqrt(max(dec["var_shared"], 1e-6))), d_syn=0.0), dec


def policies(trials=200, n_stripes=250, k=32, parity=16, seed=0):
    from d3ra.channel.surrogate import Layout, file_time, fragment_times, sample_T
    from d3ra.optimization.joint_optimizer import joint_optimise_adaptive, layout_for
    from d3ra.retrieval.adaptive_policy import AdaptiveConfig, AdaptiveController
    from d3ra.retrieval.fixed_policy import budget_for_target, coverage_equalising_weights, uniform_weights
    from d3ra.evaluation.metrics import adaptive_summary, static_summary
    _, refs = read_fasta(REF_FASTA)
    C = load_cached()
    runs = [C[r]["run"] for r in RUNS]
    belief, dec = belief_from_calibration(C, refs)
    rng = np.random.default_rng(seed)
    N = k + parity
    perm = rng.permutation(len(refs))[: n_stripes * N]        # strand ids used, grouped into stripes
    phi_full = [EM.grid_to_full(rd.phi_strand) for rd in runs]
    # strands are re-indexed: fragment f <-> design strand perm[f]
    runs_v = []
    for rd, ph in zip(runs, phi_full):
        inv = np.full(len(refs), -1)
        inv[perm] = np.arange(perm.size)
        r = rd.ref_of_read.copy()
        r = np.where(r >= 0, inv[np.maximum(r, 0)], -1)
        # reads of design strands not used are dropped (they are not in any container)
        unused = (rd.ref_of_read >= 0) & (r < 0)
        r = r[~unused]
        runs_v.append(EM.RunData(rd.name, r.astype(np.int32), rd.calib_counts[perm], rd.calib_yield))
    phi_v = [ph[perm] for ph in phi_full]
    stripe_of = np.repeat(np.arange(n_stripes), N)
    kk = np.full(n_stripes, k)
    variants = ["round_robin", "reliability", "failure1", "tiered2", "single_best"]
    best = int(np.argmax(belief.yield_ * belief.phi[:, 3] * (1 - belief.dropout)))
    js = joint_optimise_adaptive(n_stripes, k, belief, 0.99, None, parities=[parity], trials=200)
    print("joint (belief) chose", js.variant, js.history, flush=True)
    rows = []
    for v in variants + ["joint"]:
        if v == "round_robin":
            from d3ra.storage.placement import round_robin
            L = Layout(stripe_of, kk, round_robin(stripe_of, belief.m))
        elif v == "single_best":
            X = np.zeros((stripe_of.size, belief.m))
            X[:, best] = 1
            L = Layout(stripe_of, kk, X)
        elif v == "joint":
            L = js.layout
        else:
            L = layout_for(v, n_stripes, k, parity, belief, None)
        counts = L.placement.sum(0).astype(int).tolist()
        t0 = time.time()
        T = EM.empirical_T(L, runs_v, phi_v, trials, np.random.default_rng(seed + 7))
        t_T = time.time() - t0
        for pol in ["static_uniform", "static_coverage_eq", "adaptive_uniform", "adaptive_coverage", "adaptive_pivotal"]:
            if pol.startswith("static"):
                w = uniform_weights(L) if pol == "static_uniform" else coverage_equalising_weights(L, belief)
                cost = file_time(L, fragment_times(L, T, w))
                summ = static_summary(cost, (0.95, 0.99))
                # model prediction of the same quantity on the believed channel
                Ts = sample_T(L, belief, 2000, np.random.default_rng(seed + 9))
                pred = file_time(L, fragment_times(L, Ts, w))
                summ.update({"pred_cap@0.95": budget_for_target(pred, 0.95), "pred_cap@0.99": budget_for_target(pred, 0.99),
                             "pred_mean": float(pred[np.isfinite(pred)].mean()) if np.isfinite(pred).any() else np.inf})
            else:
                ctl = AdaptiveController(L, belief, AdaptiveConfig(pol.split("_")[1], cap=np.inf, batch=200))
                r = ctl.run(trials, np.random.default_rng(seed + 11), T=T, failed=np.zeros((trials, belief.m), bool))
                c = np.where(r.success, r.cost, np.inf)
                summ = adaptive_summary(c, r.success, (0.95, 0.99), rounds=r.rounds, containers_used=r.containers_used)
                Ts, fs = sample_T(L, belief, trials, np.random.default_rng(seed + 9), return_state=True)
                rp = ctl.run(trials, np.random.default_rng(seed + 12), T=Ts, failed=fs)
                summ["pred_mean_cost"] = float(rp.cost[rp.success].mean()) if rp.success.any() else np.inf
            flat = {}
            for kk2, vv in summ.items():
                if isinstance(vv, tuple):
                    flat[kk2 + "_lo"], flat[kk2 + "_hi"] = vv
                else:
                    flat[kk2] = vv
            rows.append({"placement": v, "counts": " ".join(map(str, counts)), "policy": pol, "trials": trials,
                         "T_build_s": t_T, **flat})
            print(v, pol, {x: round(flat.get(x, np.nan), 1) for x in ("mean_cost_success", "cap@0.99", "mean_spend@0.99",
                                                                        "pred_mean", "pred_mean_cost", "p_success_uncapped")}, flush=True)
        pd.DataFrame(rows).to_csv(RES / "realdata_policies.csv", index=False)
    (RES / "realdata_belief.json").write_text(json.dumps({"yield": belief.yield_.tolist(), "shape": belief.shape.tolist(),
        "dropout": belief.dropout.tolist(), "sigma_g": belief.sigma_g, "runs": [RUNS[r] for r in RUNS],
        "joint_choice": js.variant, "joint_history": js.history, "shared_variance": dec}, indent=1, default=float))
    return rows


def _real_setup(n_stripes=250, k=32, parity=16, seed=0):
    _, refs = read_fasta(REF_FASTA)
    C = load_cached()
    runs = [C[r]["run"] for r in RUNS]
    belief, dec = belief_from_calibration(C, refs)
    rng = np.random.default_rng(seed)
    N = k + parity
    perm = rng.permutation(len(refs))[: n_stripes * N]
    phi_full = [EM.grid_to_full(rd.phi_strand) for rd in runs]
    inv = np.full(len(refs), -1)
    inv[perm] = np.arange(perm.size)
    runs_v = []
    for rd in runs:
        r = np.where(rd.ref_of_read >= 0, inv[np.maximum(rd.ref_of_read, 0)], -1)
        unused = (rd.ref_of_read >= 0) & (r < 0)
        runs_v.append(EM.RunData(rd.name, r[~unused].astype(np.int32), rd.calib_counts[perm], rd.calib_yield))
    phi_v = [ph[perm] for ph in phi_full]
    return belief, runs_v, phi_v, np.repeat(np.arange(n_stripes), N), np.full(n_stripes, k)


def failure_injection(trials=60, n_stripes=250, k=32, parity=16, seed=0):
    """Knock out each container's real read stream in turn (container loss) and
    measure static feasibility and adaptive cost/success per placement."""
    from d3ra.channel.surrogate import Layout, file_time, fragment_times
    from d3ra.optimization.joint_optimizer import layout_for
    from d3ra.retrieval.adaptive_policy import AdaptiveConfig, AdaptiveController
    from d3ra.retrieval.fixed_policy import budget_for_target, uniform_weights
    from d3ra.storage.placement import round_robin
    belief, runs_v, phi_v, stripe_of, kk = _real_setup(n_stripes, k, parity, seed)
    m = belief.m
    layouts = {"round_robin": Layout(stripe_of, kk, round_robin(stripe_of, m)),
               "reliability": layout_for("reliability", n_stripes, k, parity, belief, None),
               "tiered2": layout_for("tiered2", n_stripes, k, parity, belief, None),
               "failure1 (joint choice)": layout_for("failure1", n_stripes, k, parity, belief, None)}
    rows = []
    for name, L in layouts.items():
        T = EM.empirical_T(L, runs_v, phi_v, trials, np.random.default_rng(seed + 7))
        w = uniform_weights(L)
        base = file_time(L, fragment_times(L, T, w))
        cap = 3 * budget_for_target(base, 0.95)
        for lost in [None] + list(range(m)):
            Tf = T.copy()
            if lost is not None:
                Tf[:, L.pair_cont == lost] = np.inf
            st = file_time(L, fragment_times(L, Tf, w))
            ctl = AdaptiveController(L, belief, AdaptiveConfig("coverage", cap=cap, batch=200))
            failed = np.zeros((trials, m), bool)
            if lost is not None:
                failed[:, lost] = True
            r = ctl.run(trials, np.random.default_rng(seed + 11), T=Tf, failed=failed)
            rows.append({"placement": name, "counts": " ".join(map(str, L.placement.sum(0).astype(int))),
                         "lost_container": "none" if lost is None else RUNS[list(RUNS)[lost]],
                         "static_p_recoverable": float(np.isfinite(st).mean()),
                         "adaptive_p_success_within_cap": float(r.success.mean()),
                         "adaptive_mean_cost_success": float(r.cost[r.success].mean()) if r.success.any() else np.inf,
                         "cap": cap})
            print(rows[-1], flush=True)
        pd.DataFrame(rows).to_csv(RES / "realdata_failure_injection.csv", index=False)
    return rows


def plots():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from d3ra.evaluation.scenarios import ERR_SCALES, phi_grid
    out = ROOT / "reports" / "figures"
    out.mkdir(parents=True, exist_ok=True)
    okabe = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#000000"]
    # 1) Illumina real phi per container vs simulated grid
    ch = pd.read_csv(RES / "realdata_channel.csv")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for i, row in ch.iterrows():
        ax[0].plot(EM.N_GRID, [row[f"phi_exact_n{n}"] for n in EM.N_GRID], "o-", color=okabe[i % 7], label=row["condition"])
    try:
        G = phi_grid()
        for gi, sc in enumerate(ERR_SCALES[:3]):
            ax[0].plot(range(0, 17), G[gi, :17], "--", color="grey", alpha=0.4 + 0.2 * gi, lw=1)
        ax[0].plot([], [], "--", color="grey", label="simulated (error x0.5/1/2)")
    except Exception:
        pass
    ax[0].set(xlabel="reads per strand (n)", ylabel="P(exact consensus)", title="Illumina (ENA PRJEB65931): real consensus curves", xlim=(0.5, 16))
    ax[0].legend(fontsize=7)
    nano = pd.read_csv(RES / "realdata_nanopore_phi.csv")
    ax[1].plot(nano.n_reads, nano.real_exact, "o-", color=okabe[3], label="real nanopore reads")
    ax[1].plot(nano.n_reads, nano.simulated_exact_at_measured_rates, "s--", color=okabe[0], label="simulated, i.i.d. errors at measured rates")
    ax[1].plot(nano.n_reads, nano.real_majority_vote_exact, "x:", color=okabe[6], label="real, position-wise majority vote")
    ax[1].set(xlabel="reads per strand (n)", ylabel="P(exact consensus)", title="Nanopore (Srinivasavaradhan et al.): real vs model")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "realdata_consensus_curves.png", dpi=150)
    plt.close(fig)
    # 2) policies on real recovery-time matrices
    d = pd.read_csv(RES / "realdata_policies.csv")
    d = d[d.placement.isin(["round_robin", "reliability", "tiered2", "single_best"])]
    pols = ["static_uniform", "adaptive_uniform", "adaptive_coverage", "adaptive_pivotal"]
    fig, ax = plt.subplots(figsize=(9, 4.2))
    pl = list(dict.fromkeys(d.placement))
    wdt = 0.2
    for k, pol in enumerate(pols):
        sub = d[d.policy == pol].set_index("placement").reindex(pl)
        y = sub["mean_spend@0.99"].values
        lo = y - sub["mean_spend@0.99_ci_lo"].values
        hi = sub["mean_spend@0.99_ci_hi"].values - y
        ax.bar(np.arange(len(pl)) + (k - 1.5) * wdt, y, wdt, yerr=[lo, hi], color=okabe[k], label=pol, capsize=2)
        pred = (sub["pred_mean"].fillna(sub["pred_mean_cost"])).values
        ax.plot(np.arange(len(pl)) + (k - 1.5) * wdt, pred, "k_", ms=12)
    ax.plot([], [], "k_", ms=12, label="surrogate prediction (calibrated belief)")
    ax.set_xticks(range(len(pl)), pl)
    ax.set(ylabel="reads to reach 0.99 (expected spend)", title="Real recovery-time matrices: 12,000 strands, 6 real containers")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(out / "realdata_policies.png", dpi=150)
    plt.close(fig)
    # 3) failure injection
    f = RES / "realdata_failure_injection.csv"
    if f.exists():
        fi = pd.read_csv(f)
        fig, ax = plt.subplots(figsize=(9, 4))
        for k, (name, g) in enumerate(fi.groupby("placement", sort=False)):
            ax.plot(g.lost_container, g.adaptive_p_success_within_cap, "o-", color=okabe[k], label=name)
        ax.set(ylabel="P(exact recovery within cap)", xlabel="container lost", title="Real read streams with one container removed")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out / "realdata_failure_injection.png", dpi=150)
        plt.close(fig)
    print("plots written")


def nanopore(max_clusters=600, min_reads=20, draws=4, seed=0):
    """Real nanopore phi(n) (clustered reads with ground truth) vs our simulator
    at the measured per-type error rates on the same reference strands."""
    from d3ra.channel.sequencer import ReadErrorModel, noisy_reads
    from d3ra.realdata.references import read_binned
    from d3ra.retrieval.consensus import consensus, majority_vote
    rng = np.random.default_rng(seed)
    clusters = read_binned(RAW / "zenodo" / "Srinivasavaradhan.txt")
    sizes = np.array([len(r) for _, r in clusters])
    print(f"nanopore clusters={len(clusters)} reads={sizes.sum()} median_cluster={np.median(sizes)}", flush=True)
    big = [c for c in clusters if len(c[1]) >= min_reads]
    rng.shuffle(big)
    big = big[:max_clusters]
    pairs = [(ref, rd) for ref, reads in big[:400] for rd in reads[:10]]
    prof = error_profile(pairs)
    err = ReadErrorModel(prof["rate_sub"], prof["rate_ins"], prof["rate_del"], homopolymer_boost=1.0,
                         quality_sigma=0.3, position_slope=0.0)
    n_grid = [1, 2, 3, 5, 7, 10, 15, 20]
    rows = []
    for n in n_grid:
        real_ex = real_h = sim_ex = mv_ex = 0
        tot = 0
        for ref, reads in big:
            L = len(ref)
            for _ in range(draws):
                sub = [reads[i] for i in rng.choice(len(reads), n, replace=False)]
                est = consensus(sub, L) if n > 1 else sub[0][:L]
                real_ex += est == ref
                real_h += EM._hamming_ok(est, ref, 2)
                if n > 2:
                    mv_ex += majority_vote(sub, L) == ref
                simr = noisy_reads(ref, n, err, rng)
                se = consensus(simr, L) if n > 1 else simr[0][:L]
                sim_ex += se == ref
                tot += 1
        rows.append({"n_reads": n, "real_exact": real_ex / tot, "real_hamming_le2": real_h / tot,
                     "real_majority_vote_exact": mv_ex / tot if n > 2 else np.nan,
                     "simulated_exact_at_measured_rates": sim_ex / tot, "trials": tot})
        print(rows[-1], flush=True)
    prof_small = {k: v for k, v in prof.items() if k != "position_profile"}
    pd.DataFrame(rows).to_csv(RES / "realdata_nanopore_phi.csv", index=False)
    (RES / "realdata_nanopore_channel.json").write_text(json.dumps(
        {"clusters": len(clusters), "reads": int(sizes.sum()), "median_cluster_size": float(np.median(sizes)),
         "cluster_size_cv": float(sizes.std() / sizes.mean()), "empty_or_singleton_clusters": int((sizes <= 1).sum()),
         "strand_length": int(np.median([len(r) for r, _ in clusters[:1000]])), **prof_small}, indent=1))
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stages", nargs="*", default=["prep"])
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if "download" in a.stages:
        download()
    if "prep" in a.stages:
        prep(a.workers)
    if "nanopore" in a.stages:
        nanopore()
    if "transplant" in a.stages:
        transplant()
    if "nanodiag" in a.stages:
        nanopore_diagnostics()
    if "channel" in a.stages:
        channel()
    if "policies" in a.stages:
        policies()
    if "failures" in a.stages:
        failure_injection()
    if "plots" in a.stages:
        plots()
