"""Q3b: is the 'correlated nanopore error' effect an artefact of wrong references?

A wrong reference (or a majority synthesis variant) makes every read 'err' at the
same position -- indistinguishable from correlated errors if the design is used
as truth. Tests (DNAformer two-flowcell data, and the Srinivasavaradhan set):
  1. reference audit: consensus of ALL flowcell-2 reads vs the design;
  2. per-position error-rate distribution: share of positions with rate > 0.5
     (reference-mismatch signature) vs moderately elevated rates;
  3. re-run real-vs-i.i.d. consensus curves and the edit-script transplant test
     using an independent truth (FC2 all-read consensus), restricted to strands
     whose design agrees with that truth;
  4. strand-specific residual correlation FC1 vs FC2 excluding rate>0.5 positions.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from d3ra.channel.correlated_errors import HierarchicalErrorModel, read_error_vector  # noqa: E402
from d3ra.channel.empirical_errors import ErrorTransplant  # noqa: E402
from d3ra.channel.sequencer import ReadErrorModel, noisy_reads  # noqa: E402
from d3ra.realdata.channel_estimation import error_profile  # noqa: E402
from d3ra.realdata.references import read_binned  # noqa: E402
from d3ra.retrieval.consensus import consensus  # noqa: E402

Z = ROOT / "data/raw/zenodo"
rng = np.random.default_rng(3)


def audit(A, B, ids, L, min_truth_reads=15):
    agree, total = [], 0
    truth = {}
    for s in ids:
        if len(B[s]) < min_truth_reads:
            continue
        t = consensus(B[s][:40], L)
        truth[s] = t
        agree.append(t == s)
    return truth, float(np.mean(agree)), len(agree)


def curves(A, truth, ids, n_grid, err_iid, tp, L):
    out = []
    for n in n_grid:
        real = iid = tpo = cnt = 0
        for s in ids:
            if len(A[s]) < n:
                continue
            t = truth[s]
            real += consensus([A[s][i] for i in rng.choice(len(A[s]), n, replace=False)], L) == t
            iid += consensus(noisy_reads(t, n, err_iid, rng), L) == t
            tpo += consensus(noisy_reads(t, n, tp, rng), L) == t
            cnt += 1
        out.append({"n": n, "strands": cnt, "real": real / cnt, "iid": iid / cnt, "transplant_independent": tpo / cnt})
        print(out[-1], flush=True)
    return out


def main():
    A = dict(read_binned(Z / "BinnedNanoporeFirstFlowcell_Random.txt"))
    B = dict(read_binned(Z / "BinnedNanoporeSecondFlowcell_Random.txt"))
    ids = [s for s in A if s in B and len(A[s]) >= 12 and len(B[s]) >= 15]
    rng.shuffle(ids)
    ids = ids[:3000]
    L = len(ids[0])
    truth, agree_rate, n_aud = audit(A, B, ids, L)
    print(f"design == FC2 all-read consensus for {agree_rate:.3f} of {n_aud} strands", flush=True)
    good = [s for s in ids if s in truth and truth[s] == s]
    # per-position error rates FC1 vs FC2 on good strands; reference-mismatch signature
    xa, xb = [], []
    for s in good[:600]:
        va = [v for v in (read_error_vector(s, r) for r in A[s]) if v is not None]
        vb = [v for v in (read_error_vector(s, r) for r in B[s]) if v is not None]
        xa.append((np.stack(va) > 0).mean(0)); xb.append((np.stack(vb) > 0).mean(0))
    xa, xb = np.concatenate(xa), np.concatenate(xb)
    high = (xa > 0.5) | (xb > 0.5)
    corr_all = float(np.corrcoef(xa, xb)[0, 1])
    corr_excl = float(np.corrcoef(xa[~high], xb[~high])[0, 1])
    # i.i.d. and transplant models from FC1 reads of good strands (truth = design, verified)
    pairs = [(s, r) for s in good[:600] for r in A[s][:5]]
    prof = error_profile(pairs)
    iid = ReadErrorModel(prof["rate_sub"], prof["rate_ins"], prof["rate_del"], homopolymer_boost=1.0,
                         quality_sigma=0.3, position_slope=0.0)
    tp = ErrorTransplant(pairs)
    cv = curves(A, truth, good, [5, 10], iid, tp, L)
    out = {"dnaformer": {"strands_audited": n_aud, "design_matches_independent_consensus": agree_rate,
                         "positions_rate_gt_0.5_share": float(high.mean()),
                         "fc1_fc2_rate_corr_all": corr_all, "fc1_fc2_rate_corr_excluding_gt_0.5": corr_excl,
                         "error_rate": prof["rate_total"], "curves_on_verified_strands": cv}}
    # Srinivasavaradhan set: audit with all-read consensus of the cluster's other half
    clusters = [c for c in read_binned(Z / "Srinivasavaradhan.txt") if len(c[1]) >= 40]
    rng.shuffle(clusters)
    ag = []
    for ref, reads in clusters[:400]:
        ag.append(consensus(reads[20:60], len(ref)) == ref)
    out["srinivasavaradhan"] = {"clusters_audited": len(ag), "design_matches_consensus_of_20_40_reads": float(np.mean(ag))}
    print(json.dumps(out, indent=1))
    (ROOT / "results/q3b_reference_audit.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
