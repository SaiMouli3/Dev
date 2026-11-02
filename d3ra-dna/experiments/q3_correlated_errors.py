"""Q3: correlated nanopore errors -- persistence across two flowcells,
context vs strand decomposition, and a validated hierarchical error model.
Data: Zenodo 17473983 (DNAformer), binned nanopore reads of 55,000 strands
sequenced on two separate flowcells."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from d3ra.channel.correlated_errors import HierarchicalErrorModel, context_key, read_error_vector  # noqa: E402
from d3ra.realdata.references import read_binned  # noqa: E402
from d3ra.retrieval.consensus import consensus  # noqa: E402

Z = ROOT / "data/raw/zenodo"
RES = ROOT / "results"
rng = np.random.default_rng(0)


def load():
    a = dict(read_binned(Z / "BinnedNanoporeFirstFlowcell_Random.txt"))
    b = dict(read_binned(Z / "BinnedNanoporeSecondFlowcell_Random.txt"))
    return a, b


def main(n_fit=2500, n_test=1200, min_reads=12):
    A, B = load()
    both = [s for s in A if s in B and len(A[s]) >= min_reads and len(B[s]) >= min_reads]
    rng.shuffle(both)
    fit_ids, test_ids = both[:n_fit], both[n_fit:n_fit + n_test]
    print("strands with >=%d reads in both flowcells: %d" % (min_reads, len(both)), flush=True)
    refs = {s: s for s in both}

    def vectors(D, ids):
        out = {}
        for s in ids:
            vs = [v for v in (read_error_vector(s, r) for r in D[s]) if v is not None]
            if len(vs) >= 4:
                out[s] = vs
        return out

    VA_fit, VB_fit = vectors(A, fit_ids), vectors(B, fit_ids)
    # ---- 1. context model + ICC on FC1 (fit strands)
    model = HierarchicalErrorModel.fit(VA_fit, refs)
    icc_A = model.icc
    icc_B = model.estimate_icc(VB_fit, refs)
    print("global error rate FC1 %.4f  ICC FC1 %.4f  ICC FC2 %.4f" % (model.global_rate, icc_A, icc_B), flush=True)
    # ---- 2. persistence: per (strand,position) rates FC1 vs FC2 on the same strands
    xa, xb, mu = [], [], []
    for s in VA_fit:
        if s not in VB_fit:
            continue
        ra = (np.stack(VA_fit[s]) > 0).mean(0)
        rb = (np.stack(VB_fit[s]) > 0).mean(0)
        m = np.array([model.mean_rate(s, p) for p in range(len(s))])
        xa.append(ra); xb.append(rb); mu.append(m)
    xa, xb, mu = map(np.concatenate, (xa, xb, mu))
    raw_corr = float(np.corrcoef(xa, xb)[0, 1])
    resid_corr = float(np.corrcoef(xa - mu, xb - mu)[0, 1])
    # null: residuals of FC2 at the same context but a different strand
    keys = []
    for s in VA_fit:
        if s in VB_fit:
            keys += [context_key(s, p) for p in range(len(s))]
    keys = np.array([hash(k) for k in keys])
    perm = np.arange(len(keys))
    for k in np.unique(keys):
        idx = np.flatnonzero(keys == k)
        perm[idx] = rng.permutation(idx)
    null_corr = float(np.corrcoef(xa - mu, (xb - mu)[perm])[0, 1])
    ctx_corr = float(np.corrcoef(mu, xb)[0, 1])
    print("rate corr FC1-FC2 raw %.3f | context-only vs FC2 %.3f | residual (strand-specific) %.3f | null %.3f" %
          (raw_corr, ctx_corr, resid_corr, null_corr), flush=True)
    # ---- 3. consensus failure persistence across flowcells (n = 10)
    n = 10
    fa, fb, same_wrong = [], [], 0
    for s in fit_ids[:1500]:
        if len(A[s]) < n or len(B[s]) < n:
            continue
        ca = consensus([A[s][i] for i in rng.choice(len(A[s]), n, replace=False)], len(s))
        cb = consensus([B[s][i] for i in rng.choice(len(B[s]), n, replace=False)], len(s))
        fa.append(ca != s); fb.append(cb != s)
        same_wrong += (ca != s) and (cb != s) and ca == cb
    fa, fb = np.array(fa), np.array(fb)
    p_b = fb.mean(); p_b_given_a = fb[fa].mean() if fa.any() else np.nan
    print("consensus fail FC1 %.3f FC2 %.3f  P(FC2 fail | FC1 fail) %.3f  lift %.2f  identical wrong %d/%d" %
          (fa.mean(), p_b, p_b_given_a, p_b_given_a / p_b, same_wrong, int((fa & fb).sum())), flush=True)
    # ---- 4. model validation on held-out strands (FC2 real vs models fitted on FC1 fit strands)
    grid = [3, 5, 10]
    rows = []
    test = [s for s in test_ids if len(B[s]) >= max(grid)]
    sub_rate, ins_rate, del_rate = model.global_types[1:] * model.global_rate
    from d3ra.channel.sequencer import ReadErrorModel, noisy_reads
    iid = ReadErrorModel(sub_rate, ins_rate, del_rate, homopolymer_boost=1.0, quality_sigma=0.0, position_slope=0.0)
    VA_test = vectors(A, test)
    for nn in grid:
        acc = {"real_FC2": 0, "iid": 0, "context_only": 0, "hierarchical": 0, "strand_plugin_FC1": 0}
        cnt = 0
        for s in test:
            L = len(s)
            real = consensus([B[s][i] for i in rng.choice(len(B[s]), nn, replace=False)], L)
            acc["real_FC2"] += real == s
            acc["iid"] += consensus(noisy_reads(s, nn, iid, rng), L) == s
            acc["context_only"] += consensus(model.reads(s, nn, rng, icc=0.0), L) == s
            acc["hierarchical"] += consensus(model.reads(s, nn, rng), L) == s
            if s in VA_test:   # strand-specific plug-in: propensities estimated from this strand's FC1 reads
                V = (np.stack(VA_test[s]) > 0).mean(0)
                m = np.array([model.mean_rate(s, p) for p in range(L)])
                w = len(VA_test[s]) / (len(VA_test[s]) + (1 - icc_A) / max(icc_A, 1e-6))
                p_sp = np.clip(w * V + (1 - w) * m, 1e-5, 0.9)
                acc["strand_plugin_FC1"] += consensus(model.reads(s, nn, rng, p_sp=p_sp), L) == s
            cnt += 1
        rows.append({"n_reads": nn, "strands": cnt, **{k: v / cnt for k, v in acc.items()}})
        print(rows[-1], flush=True)
    val = pd.DataFrame(rows)
    val.to_csv(RES / "q3_model_validation.csv", index=False)
    mae = {k: float((val[k] - val["real_FC2"]).abs().mean()) for k in ["iid", "context_only", "hierarchical", "strand_plugin_FC1"]}
    out = {"strands_both_flowcells": len(both), "global_error_rate_FC1": model.global_rate,
           "type_mix_sub_del_ins": model.global_types[1:].tolist(),
           "icc_FC1": icc_A, "icc_FC2": icc_B, "rate_corr_raw": raw_corr, "rate_corr_context_vs_FC2": ctx_corr,
           "rate_corr_residual_strand_specific": resid_corr, "rate_corr_residual_null": null_corr,
           "consensus_fail_FC1": float(fa.mean()), "consensus_fail_FC2": float(p_b),
           "p_fail_FC2_given_FC1": float(p_b_given_a), "identical_wrong_across_flowcells": same_wrong,
           "both_fail": int((fa & fb).sum()), "model_mae_vs_real": mae, "contexts": len(model.ctx_rate)}
    (RES / "q3_summary.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
