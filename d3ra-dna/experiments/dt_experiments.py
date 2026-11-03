"""Digital-twin end-to-end experiments (primary evidence of the revised paper).

Every container is an independent DT4DDS synthesis + aging + PCR + iSeq100 run
of the strands placed in it. Retrievals run the full decode chain (clustering,
consensus, inner RS + CRC-16 + address check, outer RS, per-stripe CRC-32,
SHA-256) with verification-gated stopping; the reported outcome of a retrieval
is whether the output file is byte-identical to the stored file. Controllers
use a generic literature prior (no calibration data of any kind).

    python experiments/dt_experiments.py gen          # synthesise replicates (cached)
    python experiments/dt_experiments.py run SCEN     # SCEN in SCENARIOS
    python experiments/dt_experiments.py calib        # rho_U selection on calibration replicates

The protocol (scenarios, replicates, permutations, policies, metrics) is frozen in
research_paper/experiments/configurations/dt_protocol.yaml before the final runs.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import pickle
import sys
import time
import zlib
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from d3ra.channel.digital_twin import ContainerCondition, simulate_container  # noqa: E402
from d3ra.ingestion.file_loader import FileRecord  # noqa: E402
from d3ra.retrieval.adaptive_policy import AdaptiveConfig  # noqa: E402
from d3ra.retrieval.stream_retrieval import StreamRetriever  # noqa: E402
from d3ra.storage.library import build_library  # noqa: E402
from d3ra.storage.placement import round_robin  # noqa: E402

RES = ROOT / "results"
CACHE = ROOT / "data" / "processed" / "dt"
if os.environ.get("D3RA_GEOM", "main") != "main":
    CACHE = CACHE / f"geom_{os.environ['D3RA_GEOM']}"
SEED = 20261102

# ---------------------------------------------------------------- frozen design
PAYLOAD, K, PARITY, NSYM, M = 24, 32, 16, 4, 6
# 18,432 payload bytes; the library bundles the file with a 162-byte header (name, size, hash) into an
# 18,594-byte blob, which fills ceil(18594 / (32*24)) = 25 stripes (the last one zero-padded):
# 25 x 48 = 1,200 strands, 200 per container, and a 60 x 1,200 = 72,000-read maximum stream.
FILE_BYTES = 24 * K * 24
READS_PER_STRAND = 60             # reads simulated per stored strand (stream length; >> reads needed)

# Revision-4 archive geometries (D3RA_GEOM); "main" is the frozen primary design above.
GEOMS = {
    "main": dict(K=32, PARITY=16, M=6, FILE_BYTES=FILE_BYTES),
    "m3": dict(K=32, PARITY=16, M=3, FILE_BYTES=FILE_BYTES),             # 3 containers, 16 fragments/stripe each
    "m12": dict(K=32, PARITY=16, M=12, FILE_BYTES=FILE_BYTES),           # 12 containers, 4 fragments/stripe each
    "n40": dict(K=32, PARITY=8, M=6, FILE_BYTES=FILE_BYTES),             # lower redundancy (32,40)
    "n64": dict(K=32, PARITY=32, M=6, FILE_BYTES=FILE_BYTES),            # higher redundancy (32,64)
    "k16": dict(K=16, PARITY=8, M=6, FILE_BYTES=FILE_BYTES),             # short stripes (16,24)
    "big": dict(K=32, PARITY=16, M=6, FILE_BYTES=4 * FILE_BYTES),        # 4x file (~100 stripes, 4,800 strands)
}
GEOM = os.environ.get("D3RA_GEOM", "main")
K, PARITY, M, FILE_BYTES = (GEOMS[GEOM][x] for x in ("K", "PARITY", "M", "FILE_BYTES"))
GTAG = "" if GEOM == "main" else f"_geom_{GEOM}"

# Handling histories mirror the six conditions of the PRJEB65931 study (aging / PCR),
# but every container is an independent synthesis.
BASE = [
    ContainerCondition("fresh", halflives=0.0, pcr_cycles=15),
    ContainerCondition("aged_1hl", halflives=1.0, pcr_cycles=15),
    ContainerCondition("aged_2hl_lowcopy", halflives=2.0, pcr_cycles=30, coverage=10),
    ContainerCondition("pcr_30c", halflives=0.0, pcr_cycles=30, pcr_eff=0.93),
    ContainerCondition("pcr_45c_lowcopy", halflives=0.0, pcr_cycles=45, pcr_eff=0.90, eff_sd=0.02, coverage=20),
    ContainerCondition("near_dead", halflives=2.0, pcr_cycles=40, pcr_eff=0.90, sigma=0.5, coverage=3),
]


def _with(conds, **changes_by_index):
    out = list(conds)
    for k, kw in changes_by_index.items():
        i = int(k[1:])
        out[i] = ContainerCondition(**{**out[i].__dict__, **kw})
    return out


SCENARIOS = {
    "baseline": BASE,                                                  # heterogeneous, no failure
    "fail1": _with(BASE, c0={"failed": True, "name": "fresh(failed)"}),  # best container lost
    "fail2": _with(BASE, c0={"failed": True, "name": "fresh(failed)"}, c3={"failed": True, "name": "pcr_30c(failed)"}),
    "homog": [ContainerCondition(f"fresh{j}") for j in range(M)],     # homogeneous control
    # deviation from the frozen protocol: DT4DDS's CustomArray model has a 117-nt position-bias table and
    # cannot synthesise 136-nt strands, so the alternative synthesis model is photolithographic
    "photolitho": [ContainerCondition(c.name + "_pl", synthesis="photolithographic", halflives=c.halflives,
                                       pcr_cycles=c.pcr_cycles, pcr_eff=c.pcr_eff, eff_sd=c.eff_sd,
                                       sigma=c.sigma, coverage=c.coverage) for c in BASE],     # different synthesis vendor model
}
# revision-4 failure modes (pre-specified in dt_protocol_rev4.yaml)
SCENARIOS.update({
    "partial_c0_50": _with(BASE, c0={"loss_frac": 0.5, "name": "fresh(50% strands lost)"}),
    "partial_c0c3_25": _with(BASE, c0={"loss_frac": 0.25, "name": "fresh(25% lost)"},
                             c3={"loss_frac": 0.25, "name": "pcr_30c(25% lost)"}),
    "lowyield_c0": _with(BASE, c0={"junk_frac": 0.9, "name": "fresh(90% junk reads)"}),
    "contam_c2": _with(BASE, c2={"junk_frac": 0.5, "name": "aged_2hl_lowcopy(50% contamination)"}),
    "common_cause": [ContainerCondition(c.name + "_cc", halflives=c.halflives + 2.0, pcr_cycles=max(c.pcr_cycles, 30),
                                        pcr_eff=min(c.pcr_eff, 0.93), eff_sd=c.eff_sd, sigma=c.sigma, coverage=c.coverage)
                     for c in BASE],                                   # correlated degradation of every container
})
# degradation sweep: severity of the last container (copies per design after handling)
for cov in (200, 30, 10, 3, 1):
    SCENARIOS[f"sweep_cov{cov}"] = _with(BASE, c5={"coverage": float(cov), "name": f"c5_cov{cov}"})

POLICIES = {   # name -> (controller, static weights or None)
    "round_robin_stop": ("uniform", None),
    "coverage": ("coverage", None),
    "pivotal": ("pivotal", None),
    "thompson": ("ts", None),
    "ucb": ("ucb", None),
    "learn_then_commit": ("ltc", None),
    "info_greedy": ("info", None),
    "tail_aware": ("tail_lcb", None),
}
N_REP, N_PERM = 30, 10             # 30 independent syntheses x 10 read orders = 300 retrievals per policy
CALIB_REPS = range(1000, 1010)     # disjoint replicates used only to select rho_U


def generic_belief(m, over=None):
    """Generic literature prior; `over` (prior-sensitivity runs) may override failure, dropout,
    shape, yield, sigma_g or the consensus-success curve scale (phi_scale)."""
    from paper_experiments import generic_prior
    b = generic_prior(m)
    for k_, v in (over or {}).items():
        if k_ == "failure":
            from d3ra.channel.surrogate import FailureModel
            b.failure = FailureModel(np.full(m, float(v)))
        elif k_ == "phi_scale":
            from d3ra.evaluation.scenarios import phi_for_scale
            b.phi = np.tile(phi_for_scale(float(v)), (m, 1))
        elif k_ == "sigma_g":
            b.sigma_g = float(v)
        else:
            setattr(b, {"yield": "yield_"}.get(k_, k_), np.full(m, float(v)))
    return b


def conditions(scen):
    """Container conditions of `scen` for the current geometry (M containers)."""
    base = SCENARIOS[scen]
    if M == len(base):
        return base
    if M == 3:
        return [base[0], base[2], base[5]]       # fresh, aged low-copy, near-dead
    return [ContainerCondition(**{**c.__dict__, "name": f"{c.name}_{r}"}) for r in range(M // len(base)) for c in base]


def archive():
    rng = np.random.default_rng(SEED)
    data = rng.integers(0, 256, FILE_BYTES, dtype=np.uint8).tobytes()
    f = FileRecord.from_bytes("dt_payload.bin", data)
    f.created = 0.0                      # deterministic archive: identical strands in every process
    lib = build_library([f], PAYLOAD, K, PARITY, NSYM, compression="none")
    X = round_robin(lib.stripe_of, M)
    from d3ra.archive import load_archive  # noqa: F401  (format objects come from the library)
    addresses = [(fr.obj, fr.stripe, fr.index) for fr in lib.fragments]
    truth = {(fr.obj, fr.stripe, fr.index): fr.payload for fr in lib.fragments}
    return lib, X, addresses, truth, data


def _rep_seed(scen: str, rep: int, j: int) -> int:
    return zlib.crc32(f"{scen}|{rep}|{j}".encode()) ^ SEED


def gen_rep(args):
    scen, rep = args
    out = CACHE / scen / f"rep{rep:04d}.pkl.gz"
    if out.exists():
        return str(out), "cached"
    lib, X, *_ = archive()
    conds = conditions(scen)
    streams = []
    for j in range(M):
        idx = np.flatnonzero(X[:, j] > 0)
        seqs = [lib.fragments[i].strand for i in idx]
        streams.append(simulate_container(seqs, conds[j], READS_PER_STRAND * len(seqs), _rep_seed(scen, rep, j)))
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp")
    with gzip.open(tmp, "wb") as fh:   # write-then-rename: an interrupted run never leaves a truncated cache
        pickle.dump(streams, fh)
    tmp.replace(out)
    return str(out), "new"


def load_rep(scen, rep):
    with gzip.open(CACHE / scen / f"rep{rep:04d}.pkl.gz", "rb") as fh:
        return pickle.load(fh)


def one_retrieval(args):
    scen, rep, perm, pol, cfg_kw, cap = args
    lib, X, addresses, truth, data = archive()
    streams = load_rep(scen, rep)
    prng = np.random.default_rng(zlib.crc32(f"perm|{scen}|{rep}|{perm}".encode()))
    streams = [[s[i] for i in prng.permutation(len(s))] for s in streams]
    cfg_kw = dict(cfg_kw)
    meta = {k_: v for k_, v in cfg_kw.items() if k_.startswith("_") and k_ != "_checkpoints"}
    cps = cfg_kw.pop("_checkpoints", None)
    belief_over = cfg_kw.pop("_belief", None)
    mcr = cfg_kw.pop("_max_cluster_reads", 24)
    w = cfg_kw.pop("static_weights", None)
    ctl = POLICIES[pol][0] if pol in POLICIES else "uniform"
    cfg = AdaptiveConfig(controller=ctl, batch=50, cap=cap, **cfg_kw)
    r = StreamRetriever(lib.manifest, lib.fmt, X, addresses, streams, generic_belief(M, belief_over), cfg, truth=truth,
                        static_weights=w, cap_reads=cap, checkpoints=cps, max_cluster_reads=mcr)
    out = r.run(np.random.default_rng(zlib.crc32(f"ctl|{scen}|{rep}|{perm}|{pol}".encode())))
    files = out.pop("files")
    exact_bytes = bool(files) and files[0].data == data
    return {"scenario": scen, "rep": rep, "perm": perm, "policy": pol, **{k: v for k, v in out.items() if k != "checkpoints"},
            "file_identical": exact_bytes, "wrong_output": bool(files) and not exact_bytes,
            "checkpoints": json.dumps(out["checkpoints"]), "geom": GEOM,
            **{f"cfg_{k}": v for k, v in cfg_kw.items()}, **{f"cfg{k}": json.dumps(v) for k, v in meta.items()}}


# ---------------------------------------------------------------- static baselines
def static_family(scen, workers=4):
    """Static allocations with stopping disabled are characterised by the reads each one needs
    to decode (reading in fixed proportions w until decodable); a static plan must commit the
    0.99-quantile of that requirement as its budget. Family: w_j ~ e_j^alpha with e_j the
    hindsight per-read recovery rate of container j measured on the same replicates
    (alpha = 0 is uniform fixed depth; the best alpha is a hindsight-tuned static plan)."""
    eff = np.array([0.0] * M)
    # hindsight efficiency: recovered fragments per read in a single-container pass, measured on the
    # evaluation replicates; with D3RA_STATIC_EFF=calib it is measured on the disjoint calibration
    # replicates instead (a practical static plan that sees no evaluation data)
    src = os.environ.get("D3RA_STATIC_EFF", "hindsight")
    eff_reps = list(CALIB_REPS) if src == "calib" else list(range(min(N_REP, 10)))
    for r_ in eff_reps:
        gen_rep((scen, r_))
    lib, X, addresses, truth, data = archive()
    for j in range(M):
        tot = rec = 0
        for rep in eff_reps:
            s = load_rep(scen, rep)[j]
            cur = container_curve(lib, X, addresses, s, j, upto=8 * int(X[:, j].sum()))
            tot += cur[-1][0]
            rec += cur[-1][1]
        eff[j] = rec / max(tot, 1)
    rows = []
    alphas = [float(a) for a in os.environ.get("D3RA_STATIC_ALPHAS", "0,0.5,1,2").split(",")]
    tag = os.environ.get("D3RA_STATIC_TAG", "")
    pfx = "static_cal_a" if src == "calib" else "static_a"
    rep0 = int(os.environ.get("D3RA_REP0", 0))
    for alpha in alphas:
        w = np.where(eff > 0, np.maximum(eff, 1e-6) ** alpha, 0.0)
        w = w / w.sum()
        n_rep = int(os.environ.get("D3RA_NREP", N_REP))
        n_perm = int(os.environ.get("D3RA_NPERM", N_PERM))
        jobs = [(scen, r, p, f"{pfx}{alpha}", {"static_weights": w}, np.inf)
                for r in range(rep0, rep0 + n_rep) for p in range(n_perm)]
        with Pool(workers) as pool:
            for row in pool.imap_unordered(one_retrieval, jobs, chunksize=2):
                row["reads_per_container"] = json.dumps(row["reads_per_container"])
                row["alpha"] = alpha
                row["eff_source"] = src
                row["weights"] = json.dumps(w.round(4).tolist())
                rows.append(row)
        print("static", scen, alpha, flush=True)
    pd.DataFrame(rows).to_csv(RES / f"dt_{scen}{GTAG}_static{tag}.csv", index=False)
    (RES / f"dt_{scen}{GTAG}_{'calib' if src == 'calib' else 'hindsight'}_eff{tag}.json").write_text(json.dumps(eff.tolist()))


# ---------------------------------------------------------------- calibration (rho_U) on disjoint replicates
def container_curve(lib, X, addresses, stream, j, upto, batch=50):
    """Read container j alone in batches; return [(reads, own recoveries)]."""
    from d3ra.channel.surrogate import Layout  # noqa: F401
    one = [[] for _ in range(M)]
    one[j] = stream
    r = StreamRetriever(lib.manifest, lib.fmt, X, addresses, one, generic_belief(M),
                        AdaptiveConfig(controller="uniform", batch=batch))
    out = []
    while r.reads_used[j] < min(upto, len(stream)):
        r.read_container(j, batch)
        out.append((r.reads_used[j], sum(1 for (i, jj), v in r.rec_pair.items() if jj == j and v)))
    return out


def _theta_star(ctl, j, curve):
    """Efficiency multiplier that best reproduces the container's full recovery curve
    under the controller's own (generic) model: the reference value for calibration."""
    from d3ra.retrieval.adaptive_policy import THETA  # noqa: F401
    grid = np.exp(np.linspace(np.log(0.01), np.log(5), 400))
    R = np.array([c[0] for c in curve])
    u = np.array([c[1] for c in curve]) / ctl.n_j[j]
    best, val = None, np.inf
    for g in grid:
        q = ctl.q.per_container(np.outer(R * ctl.pi_cont[j] * g, np.ones(M)))[:, j]
        e = ((q - u) ** 2).sum()
        if e < val:
            best, val = g, e
    return best


def calib_job(args):
    scen, rep, rhos = args
    from d3ra.retrieval.adaptive_policy import THETA, AdaptiveController
    from d3ra.channel.surrogate import Layout
    lib, X, addresses, truth, data = archive()
    streams = load_rep(scen, rep)
    k = np.array([s["k"] for s in lib.manifest.stripes]); n = np.array([s["n"] for s in lib.manifest.stripes])
    L = Layout(np.repeat(np.arange(len(k)), n), k, X)
    ctl0 = AdaptiveController(L, generic_belief(M), AdaptiveConfig())
    rows = []
    for j in range(M):
        curve = container_curve(lib, X, addresses, streams[j], j, upto=10 * int(X[:, j].sum()))
        th = _theta_star(ctl0, j, curve)
        # R*: reads at which the generic belief predicts 40% recovery (theta = 1)
        Rgrid = np.arange(50, 20000, 50)
        q = ctl0.q.per_container(np.outer(Rgrid * ctl0.pi_cont[j], np.ones(M)))[:, j]
        rstar = Rgrid[np.searchsorted(q, 0.4)] if (q >= 0.4).any() else Rgrid[-1]
        for mult in (1, 2):
            Rj = mult * rstar
            U = [c[1] for c in curve if c[0] <= Rj]
            U = U[-1] if U else 0
            for rho in rhos:
                ctl = AdaptiveController(L, generic_belief(M), AdaptiveConfig(prior_sd=CALIB_SD, overdispersion=rho))
                Rm = np.zeros((1, M)); Rm[0, j] = Rj
                Um = np.zeros((1, M)); Um[0, j] = U
                P = ctl.posterior(Rm, Um)[0, j]
                cdf = np.cumsum(P)
                lo = THETA[min((cdf < 0.05).sum(), len(THETA) - 1)]
                hi = THETA[min((cdf < 0.95).sum(), len(THETA) - 1)]
                rows.append({"scenario": scen, "rep": rep, "container": j, "op": f"{mult}xR*", "reads": Rj, "U": U,
                             "theta_star": th, "rho": rho, "post_mean": float(P @ THETA), "lo90": lo, "hi90": hi,
                             "covered": bool(lo <= th <= hi), "p_below_0.7": float(P[THETA < 0.7].sum())})
    return rows


RHOS = tuple(float(x) for x in os.environ.get("D3RA_RHOS", "0,0.003,0.01,0.03,0.1").split(","))
CALIB_SD = float(os.environ.get("D3RA_CALIB_SD", 1.0))   # frozen selection used 1.0; revision 4 repeats it at 0.5


def calib(scen="baseline", reps=CALIB_REPS, out="dt_calibration_select", workers=4):
    jobs = [(scen, r, RHOS) for r in reps]
    for r in reps:
        gen_rep((scen, r))
    rows = []
    with Pool(workers) as p:
        for rr_ in p.imap_unordered(calib_job, jobs):
            rows += rr_
    df = pd.DataFrame(rows)
    df["prior_sd"] = CALIB_SD
    df.to_csv(RES / f"{out}_{scen}{os.environ.get('D3RA_CALIB_TAG', '')}.csv", index=False)
    s = df.groupby("rho").covered.mean()
    print(s, flush=True)
    return df


# ---------------------------------------------------------------- completion-forecast calibration
FRACS = (0.3, 0.5, 0.7, 0.85)


def forecast_eval(scen="baseline", workers=4):
    """Budget = median cost of the coverage controller in `scen`; record P_complete at fixed
    fractions of it under the analytic (normal) and the posterior-simulation (mc) forecast."""
    base = pd.read_csv(RES / f"dt_{scen}.csv")
    cap = float(base[base.policy == "coverage"].reads.median())
    cps = [f * cap for f in FRACS]
    rows = []
    for fc in ("normal", "mc"):
        jobs = [(scen, r, p, "coverage", {"forecast": fc, "_checkpoints": list(cps)}, cap)
                for r in range(N_REP) for p in range(N_PERM)]
        with Pool(workers) as pool:
            for row in pool.imap_unordered(one_retrieval, jobs, chunksize=2):
                done = row["file_identical"]
                for level, pc in json.loads(row["checkpoints"]):
                    rows.append({"forecast": fc, "rep": row["rep"], "perm": row["perm"], "frac": round(level / cap, 2),
                                 "p_complete": pc, "completed": done, "cap": cap})
    df = pd.DataFrame(rows)
    df.to_csv(RES / f"dt_forecast_{scen}.csv", index=False)
    return df


def gen(scens=None, workers=4):
    rep0 = int(os.environ.get("D3RA_REP0", 0))
    n_rep = int(os.environ.get("D3RA_NREP", N_REP))
    jobs = [(s, r) for s in (scens or SCENARIOS) for r in range(rep0, rep0 + n_rep)]
    with Pool(workers) as p:
        for path, msg in p.imap_unordered(gen_rep, jobs):
            print("gen", path, msg, flush=True)


def run(scen, policies=None, workers=4, tag="", cfg_kw=None):
    pols = policies or list(POLICIES)
    n_rep = int(os.environ.get("D3RA_NREP", N_REP))     # reduced designs are logged as protocol deviations
    n_perm = int(os.environ.get("D3RA_NPERM", N_PERM))
    rep0 = int(os.environ.get("D3RA_REP0", 0))
    jobs = [(scen, r, p, pol, dict(cfg_kw or {}), np.inf) for pol in pols for r in range(rep0, rep0 + n_rep)
            for p in range(n_perm)]
    rows = []
    t0 = time.time()
    with Pool(workers) as pool:
        for k, row in enumerate(pool.imap_unordered(one_retrieval, jobs, chunksize=2)):
            row["reads_per_container"] = json.dumps(row["reads_per_container"])
            rows.append(row)
            if k % 50 == 0:
                print(scen, k, len(jobs), f"{time.time() - t0:.0f}s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(RES / f"dt_{scen}{GTAG}{tag}.csv", index=False)
    return df


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "gen":
        gen(sys.argv[2:] or None, int(os.environ.get("D3RA_WORKERS", 4)))
    elif cmd == "run":
        run(sys.argv[2], sys.argv[3:] or None, int(os.environ.get("D3RA_WORKERS", 4)))
    elif cmd == "posthoc_rho0":     # POST-HOC ablation (not in the frozen protocol): binomial likelihood
        run(sys.argv[2], ["coverage", "thompson", "tail_aware", "pivotal"], int(os.environ.get("D3RA_WORKERS", 4)),
            tag="_posthoc_rho0", cfg_kw={"overdispersion": 0.0})
    elif cmd == "runx":             # revision-4 runs: runx SCEN TAG JSON_CFG POLICY...
        run(sys.argv[2], sys.argv[5:] or None, int(os.environ.get("D3RA_WORKERS", 4)), tag=sys.argv[3],
            cfg_kw=json.loads(sys.argv[4]))
    elif cmd == "static":
        static_family(sys.argv[2], int(os.environ.get("D3RA_WORKERS", 4)))
    elif cmd == "forecast":
        forecast_eval(sys.argv[2] if len(sys.argv) > 2 else "baseline", int(os.environ.get("D3RA_WORKERS", 4)))
    elif cmd == "calib":
        calib(sys.argv[2] if len(sys.argv) > 2 else "baseline")
    elif cmd == "calib_eval":       # held-out replicates (the experiment replicates) after selection
        calib(sys.argv[2], range(N_REP), out="dt_calibration_eval")
