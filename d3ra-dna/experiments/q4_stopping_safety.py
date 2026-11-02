"""Q4: safety of verification-gated adaptive stopping.

A. Analytic per-attempt false-accept bound for the inner checks.
B. Empirical false-accept rates of the inner code on realistic wrong consensus
   words (nanopore-error clusters with few reads), for weak and default codes.
C. End-to-end retrievals with a deliberately weak inner code: false accepts
   occur, are caught by the per-stripe CRC / file SHA-256, and repaired;
   wrong files are never output.
D. Futility (give-up) rule: how often it abandons a retrieval that would have
   succeeded within the cap.
"""
import json
import math
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from d3ra.channel.sequencer import ReadErrorModel, noisy_reads  # noqa: E402
from d3ra.decoding.dna_decoder import StrandFormat  # noqa: E402
from d3ra.fragmentation.addressing import AddressLayout  # noqa: E402
from d3ra.retrieval.consensus import consensus  # noqa: E402

RES = ROOT / "results"
rng = np.random.default_rng(0)


def ball(n, t, q=256):
    return sum(math.comb(n, i) * (q - 1) ** i for i in range(t + 1))


def analytic_bound(n_bytes, nsym, crc_bytes, addr_check=1):
    """Upper-bound-style estimate for a *random* wrong word: P(RS decoder lands in a
    codeword's radius-t ball) * P(CRC passes) * P(address check passes | address damaged).
    The address factor is omitted (conservative) because errors may spare the header."""
    t = nsym // 2
    p_rs = min(1.0, ball(n_bytes, t) / 256 ** nsym) if nsym > 0 else 1.0
    return p_rs * 2.0 ** (-8 * crc_bytes)


def empirical(fmt, n_attempts, n_reads, err, seed=0):
    r = np.random.default_rng(seed)
    wrong = accepted_wrong = accepted_right = 0
    L = fmt.strand_nt
    for t in range(n_attempts):
        payload = r.integers(0, 256, fmt.payload_bytes, dtype=np.uint8).tobytes()
        addr = (0, int(r.integers(0, 2 ** 15)), int(r.integers(0, 255)))
        strand = fmt.encode(*addr, payload)
        est = consensus(noisy_reads(strand, n_reads, err, r), L) if n_reads > 1 else noisy_reads(strand, 1, err, r)[0][:L]
        dec = fmt.decode(est)
        if est != strand:
            wrong += 1
        if dec is not None:
            if dec[0] == addr and dec[1] == payload:
                accepted_right += 1
            else:
                accepted_wrong += 1
    return wrong, accepted_wrong, accepted_right


def main():
    out = {}
    # ---- A + B
    tp = None
    tpf = ROOT / "data/processed/realdata/transplant_nanopore.pkl"
    if tpf.exists():
        with open(tpf, "rb") as fh:
            tp = pickle.load(fh)
    err = tp if tp is not None else ReadErrorModel(0.015, 0.015, 0.02, quality_sigma=0.3)
    rows = []
    for nsym, crc in [(0, 0), (2, 0), (2, 1), (4, 0), (4, 2), (6, 2)]:
        fmt = StrandFormat.build(AddressLayout(1, 2, 1), 32, nsym, "scrambled-2bit", crc_bytes=crc)
        n_att = 4000
        w, aw, ar = empirical(fmt, n_att, 3, err, seed=nsym * 10 + crc)
        ub = 3.0 / w if aw == 0 and w else np.nan            # rule of three, 95%
        rows.append({"inner_nsym": nsym, "crc_bytes": crc, "strand_bytes": fmt.strand_bytes, "attempts": n_att,
                     "wrong_consensus": w, "false_accepts": aw, "correct_accepts": ar,
                     "false_accept_per_wrong_word": aw / max(w, 1), "upper95_if_zero": ub,
                     "analytic_estimate": analytic_bound(fmt.strand_bytes, nsym, crc)})
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(RES / "q4_inner_false_accepts.csv", index=False)
    # ---- C: end-to-end with a weak inner code (nsym=2, crc=0) vs default
    import os
    from d3ra.channel.failures import FailureModel
    from d3ra.channel.surrogate import ChannelParams
    from d3ra.ingestion.file_loader import FileRecord
    from d3ra.retrieval.adaptive_policy import AdaptiveConfig
    from d3ra.retrieval.sequencing_env import SequencingEnv, run_policy
    from d3ra.storage.container import make_heterogeneous
    from d3ra.storage.library import build_library
    from d3ra.storage.placement import round_robin
    e2e = []
    data = os.urandom(6000)
    for nsym, crc, scale in [(2, 0, 3.0), (4, 2, 3.0)]:
        for seed in range(12):
            lib = build_library([FileRecord.from_bytes("f.bin", data)], 32, 16, 8, nsym, crc_bytes=crc, compression="none")
            cs = make_heterogeneous(4, np.random.default_rng(seed), 0.5)
            for c in cs:
                c.e_sub, c.e_ins, c.e_del = 0.004 * scale, 0.001 * scale, 0.003 * scale
                c.fail_prob = 0.0
            X = round_robin(lib.stripe_of, 4)
            ch = ChannelParams.from_containers(cs, FailureModel(np.zeros(4)), L=lib.fmt.strand_nt, nsym=nsym)
            env = SequencingEnv(lib, cs, X, np.random.default_rng(100 + seed), failed=np.zeros(4, bool))
            r = run_policy(env, ch, AdaptiveConfig("pivotal", batch=60), max_actions=3000)
            out_ok = r["status"] == "exact" and r["result"].files[0].data == data
            e2e.append({"inner_nsym": nsym, "crc_bytes": crc, "seed": seed, "status": r["status"],
                        "output_correct": out_ok, "wrong_file_output": (r["status"] == "exact") and not out_ok,
                        "decode_attempts": r["decode_attempts"], "false_accepts": r["false_accepts"],
                        "verification_rounds": r["verification_rounds"], "stripes_repaired": len(r["repaired"]),
                        "reads": r["reads"]})
            print(e2e[-1], flush=True)
    E = pd.DataFrame(e2e)
    E.to_csv(RES / "q4_end_to_end_weak_code.csv", index=False)
    # ---- D: futility rule on the surrogate (channel-change scenario)
    from d3ra.channel.surrogate import sample_T
    from d3ra.evaluation import scenarios as S
    from d3ra.evaluation.baselines import uniform_layout
    from d3ra.retrieval.adaptive_policy import AdaptiveController
    sc = S.make("s7_channel_change")
    L = uniform_layout(sc)
    T, f = sample_T(L, sc.truth, 600, np.random.default_rng(5), return_state=True)
    cap = 150000
    base = AdaptiveController(L, sc.belief, AdaptiveConfig("pivotal", cap=cap, prior_sd=1.0), truth=sc.truth).run(600, np.random.default_rng(1), T=T, failed=f)
    fut = []
    for eps in [1e-3, 1e-2, 5e-2]:
        r = AdaptiveController(L, sc.belief, AdaptiveConfig("pivotal", cap=cap, prior_sd=1.0, futility_eps=eps),
                               truth=sc.truth).run(600, np.random.default_rng(1), T=T, failed=f)
        gave = r.meta["gave_up"]
        fut.append({"futility_eps": eps, "gave_up": float(gave.mean()),
                    "wrongful_give_up": float((gave & base.success).mean()),
                    "success": float(r.success.mean()), "success_without_futility": float(base.success.mean()),
                    "mean_cost": float(r.cost.mean()), "mean_cost_without_futility": float(base.cost.mean())})
        print(fut[-1], flush=True)
    pd.DataFrame(fut).to_csv(RES / "q4_futility.csv", index=False)
    out = {"inner": rows, "end_to_end": {
        f"nsym{n}_crc{c}": {"runs": int(len(g)), "exact": int((g.status == "exact").sum()),
                            "wrong_file_output": int(g.wrong_file_output.sum()),
                            "false_accepts_total": int(g.false_accepts.sum()),
                            "decode_attempts_total": int(g.decode_attempts.sum()),
                            "runs_with_repair": int((g.stripes_repaired > 0).sum())}
        for (n, c), g in E.groupby(["inner_nsym", "crc_bytes"])}, "futility": fut}
    (RES / "q4_summary.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out["end_to_end"], indent=1))


if __name__ == "__main__":
    main()
