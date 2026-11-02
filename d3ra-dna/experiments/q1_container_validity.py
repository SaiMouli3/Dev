"""Q1: are the six ENA runs valid stand-ins for independent containers?
Decompose between-sample coverage and error differences into technical
(replicate library on the same flowcell), handling (aging/PCR, same flowcell)
and sequencing-run (different flowcell) components."""
import json
import pickle
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from d3ra.io.fastq import read_fasta, read_fastq  # noqa: E402
from d3ra.realdata.channel_estimation import error_profile  # noqa: E402
from d3ra.realdata.mapping import KmerMapper, trim_adapter  # noqa: E402

RAW = ROOT / "data/raw"
REF = RAW / "zenodo/fml-ethz-dt4dds-21968d9/scripts/internal_validation/Twist_GCall/design_files.fasta"
RUNS = {"ERR12033805": ("Aging_0d_repA", 123), "ERR12033809": ("Aging_0d_repB", 123),
        "ERR12033813": ("Aging_2d", 123), "ERR12033821": ("Aging_7d", 123),
        "ERR12033826": ("PCR_15c", 121), "ERR12033837": ("PCR_60c", 121), "ERR12033841": ("PCR_90c", 121)}
N_READS = 600_000


def profile(run, mapper, refs):
    counts = np.zeros(len(refs), np.int64)
    pairs = []
    rng = np.random.default_rng(0)
    for rec in read_fastq(RAW / "ena" / f"{run}_R1.fq.gz", limit=N_READS):
        s = trim_adapter(rec.seq)
        m = mapper.map(s)
        if m.ref >= 0:
            counts[m.ref] += 1
            if len(pairs) < 15000 and rng.random() < 0.05:
                pairs.append((refs[m.ref], s[:116]))
    prof = error_profile(pairs)
    return counts, {k: v for k, v in prof.items() if k != "position_profile"}


def main():
    _, refs = read_fasta(REF)
    mapper = KmerMapper(refs)
    C, E = {}, {}
    for run in RUNS:
        C[run], E[run] = profile(run, mapper, refs)
        print(run, RUNS[run][0], "rate_sub", round(E[run]["rate_sub"], 5), "mapped", C[run].sum(), flush=True)
    rows = []
    for a, b in combinations(RUNS, 2):
        ca, cb = C[a].astype(float), C[b].astype(float)
        ok = (ca > 0) & (cb > 0)
        la, lb = np.log(ca[ok] / ca[ok].mean()), np.log(cb[ok] / cb[ok].mean())
        poisson = (1 / ca[ok]).mean() + (1 / cb[ok]).mean()
        excess = float(np.var(la - lb) - poisson)
        same_fc = RUNS[a][1] == RUNS[b][1]
        na, nb = RUNS[a][0], RUNS[b][0]
        kind = ("technical replicate" if {na, nb} == {"Aging_0d_repA", "Aging_0d_repB"}
                else "same flowcell, different handling" if same_fc else "different flowcell")
        rows.append({"a": na, "b": nb, "comparison": kind, "excess_logratio_var": max(excess, 0.0),
                     "log_corr": float(np.corrcoef(la, lb)[0, 1]),
                     "sub_rate_diff": abs(E[a]["rate_sub"] - E[b]["rate_sub"]),
                     "del_rate_diff": abs(E[a]["rate_del"] - E[b]["rate_del"])})
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "results/q1_pairwise.csv", index=False)
    summ = df.groupby("comparison")[["excess_logratio_var", "log_corr", "sub_rate_diff", "del_rate_diff"]].agg(["mean", "min", "max", "count"])
    print(summ.round(5).to_string())
    errs = pd.DataFrame({RUNS[r][0]: E[r] for r in RUNS}).T
    errs.to_csv(ROOT / "results/q1_error_rates.csv")
    out = {"pairwise_summary": json.loads(summ.round(6).to_json()),
           "headers": {RUNS[r][0]: f"instrument FS10000169 run {RUNS[r][1]}" for r in RUNS}}
    (ROOT / "results/q1_summary.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
