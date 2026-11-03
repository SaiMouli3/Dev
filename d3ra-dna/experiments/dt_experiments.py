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
SEED = 20261102

# ---------------------------------------------------------------- frozen design
PAYLOAD, K, PARITY, NSYM, M = 24, 32, 16, 4, 6
FILE_BYTES = 24 * K * 24          # 24 data stripes -> 24 stripes x 48 = 1,152 strands, 192 per container
READS_PER_STRAND = 60             # reads simulated per stored strand (stream length; >> reads needed)

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


def generic_belief(m):
    from paper_experiments import generic_prior
    return generic_prior(m)


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
    conds = SCENARIOS[scen]
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
    cps = cfg_kw.pop("_checkpoints", None)
    w = cfg_kw.pop("static_weights", None)
    ctl = POLICIES[pol][0] if pol in POLICIES else "uniform"
    cfg = AdaptiveConfig(controller=ctl, batch=50, cap=cap, **cfg_kw)
    r = StreamRetriever(lib.manifest, lib.fmt, X, addresses, streams, generic_belief(M), cfg, truth=truth,
                        static_weights=w, cap_reads=cap, checkpoints=cps)
    out = r.run(np.random.default_rng(zlib.crc32(f"ctl|{scen}|{rep}|{perm}|{pol}".encode())))
    files = out.pop("files")
    exact_bytes = bool(files) and files[0].data == data
    return {"scenario": scen, "rep": rep, "perm": perm, "policy": pol, **{k: v for k, v in out.items() if k != "checkpoints"},
            "file_identical": exact_bytes, "wrong_output": bool(files) and not exact_bytes,
            "checkpoints": json.dumps(out["checkpoints"]), **{f"cfg_{k}": v for k, v in cfg_kw.items()}}


# ---------------------------------------------------------------- static baselines
def static_family(scen, workers=4):
    """Static allocations with stopping disabled are characterised by the reads each one needs
    to decode (reading in fixed proportions w until decodable); a static plan must commit the
    0.99-quantile of that requirement as its budget. Family: w_j ~ e_j^alpha with e_j the
    hindsight per-read recovery rate of container j measured on the same replicates
    (alpha = 0 is uniform fixed depth; the best alpha is a hindsight-tuned static plan)."""
    rr = pd.read_csv(RES / f"dt_{scen}.csv")
    rr = rr[rr.policy == "round_robin_stop"]
    rpc = np.stack(rr.reads_per_container.map(json.loads).values)
    eff = np.array([0.0] * M)
    # hindsight efficiency: recovered fragments per read in a single-container pass
    lib, X, addresses, truth, data = archive()
    for j in range(M):
        tot = rec = 0
        for rep in range(min(N_REP, 10)):
            s = load_rep(scen, rep)[j]
            cur = container_curve(lib, X, addresses, s, j, upto=8 * int(X[:, j].sum()))
            tot += cur[-1][0]
            rec += cur[-1][1]
        eff[j] = rec / max(tot, 1)
    rows = []
    for alpha in (0.0, 0.5, 1.0, 2.0):
        w = np.where(eff > 0, np.maximum(eff, 1e-6) ** alpha, 0.0)
        w = w / w.sum()
        jobs = [(scen, r, p, f"static_a{alpha}", {"static_weights": w}, np.inf) for r in range(N_REP) for p in range(N_PERM)]
        with Pool(workers) as pool:
            for row in pool.imap_unordered(one_retrieval, jobs, chunksize=2):
                row["reads_per_container"] = json.dumps(row["reads_per_container"])
                row["alpha"] = alpha
                row["weights"] = json.dumps(w.round(4).tolist())
                rows.append(row)
        print("static", scen, alpha, flush=True)
    pd.DataFrame(rows).to_csv(RES / f"dt_{scen}_static.csv", index=False)
    (RES / f"dt_{scen}_hindsight_eff.json").write_text(json.dumps(eff.tolist()))


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
                ctl = AdaptiveController(L, generic_belief(M), AdaptiveConfig(prior_sd=1.0, overdispersion=rho))
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


RHOS = (0.0, 0.003, 0.01, 0.03, 0.1)


def calib(scen="baseline", reps=CALIB_REPS, out="dt_calibration_select", workers=4):
    jobs = [(scen, r, RHOS) for r in reps]
    for r in reps:
        gen_rep((scen, r))
    rows = []
    with Pool(workers) as p:
        for rr_ in p.imap_unordered(calib_job, jobs):
            rows += rr_
    df = pd.DataFrame(rows)
    df.to_csv(RES / f"{out}_{scen}.csv", index=False)
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
    jobs = [(s, r) for s in (scens or SCENARIOS) for r in range(N_REP)]
    with Pool(workers) as p:
        for path, msg in p.imap_unordered(gen_rep, jobs):
            print("gen", path, msg, flush=True)


def run(scen, policies=None, workers=4, tag="", cfg_kw=None):
    pols = policies or list(POLICIES)
    jobs = [(scen, r, p, pol, dict(cfg_kw or {}), np.inf) for pol in pols for r in range(N_REP) for p in range(N_PERM)]
    rows = []
    t0 = time.time()
    with Pool(workers) as pool:
        for k, row in enumerate(pool.imap_unordered(one_retrieval, jobs, chunksize=2)):
            row["reads_per_container"] = json.dumps(row["reads_per_container"])
            rows.append(row)
            if k % 50 == 0:
                print(scen, k, len(jobs), f"{time.time() - t0:.0f}s", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(RES / f"dt_{scen}{tag}.csv", index=False)
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
    elif cmd == "static":
        static_family(sys.argv[2], int(os.environ.get("D3RA_WORKERS", 4)))
    elif cmd == "forecast":
        forecast_eval(sys.argv[2] if len(sys.argv) > 2 else "baseline", int(os.environ.get("D3RA_WORKERS", 4)))
    elif cmd == "calib":
        calib(sys.argv[2] if len(sys.argv) > 2 else "baseline")
    elif cmd == "calib_eval":       # held-out replicates (the experiment replicates) after selection
        calib(sys.argv[2], range(N_REP), out="dt_calibration_eval")
