#!/bin/bash
# Scenario-level clean reproduction (revision 4, reviewer item 33).
#
# Clones the repository at the current HEAD into a scratch directory, regenerates the primary
# baseline scenario (all 8 policies x 30 syntheses x 10 read orders) and its static family from
# seeds alone, and compares every recorded field against the committed result files.
#
#   bash research_paper/experiments/scripts/clean_repro.sh [workdir]
#
# Writes research_paper/evidence/verification_reports/clean_repro.txt.
set -u
SRC=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
PROJ=$SRC/d3ra-dna
[ -d "$PROJ" ] || PROJ=$SRC
WORK=${1:-/tmp/clean_repro}
OUT=$PROJ/research_paper/evidence/verification_reports/clean_repro.txt
COMMIT=$(git -C "$SRC" rev-parse HEAD)

rm -rf "$WORK"; mkdir -p "$WORK"
git -C "$SRC" clone -q --no-hardlinks "$SRC" "$WORK/repo"
git -C "$WORK/repo" checkout -q "$COMMIT"
CLONE=$WORK/repo/d3ra-dna
[ -d "$CLONE" ] || CLONE=$WORK/repo

# The DT4DDS synthesis cache is regenerated from seeds in the clone (nothing is copied).
cd "$CLONE" || exit 1
export D3RA_WORKERS=${D3RA_WORKERS:-3} PYTHONHASHSEED=0
python3 experiments/dt_experiments.py gen baseline  || exit 1
python3 experiments/dt_experiments.py run baseline  || exit 1
python3 experiments/dt_experiments.py static baseline || exit 1

python3 - "$PROJ" "$CLONE" "$COMMIT" "$OUT" <<'PY'
import sys
import pandas as pd

proj, clone, commit, out = sys.argv[1:5]
KEY = ["scenario", "rep", "perm", "policy"]
FIELDS = ["status", "reads", "reads_per_container", "rounds", "decode_attempts", "false_accepts",
          "verification_rounds", "file_identical", "wrong_output"]
lines = [f"Scenario-level clean reproduction of the primary baseline scenario.",
         f"Repository cloned at commit {commit} into {clone}; DT4DDS replicates regenerated from seeds.",
         f"Compared fields: {', '.join(FIELDS)}.", ""]
bad = 0
for f in ("dt_baseline.csv", "dt_baseline_static.csv"):
    a = pd.read_csv(f"{proj}/results/{f}").sort_values(KEY).reset_index(drop=True)
    b = pd.read_csv(f"{clone}/results/{f}").sort_values(KEY).reset_index(drop=True)
    if len(a) != len(b):
        lines.append(f"{f}: ROW COUNT {len(a)} -> {len(b)}"); bad += 1; continue
    diff = [c for c in FIELDS if not a[c].astype(str).equals(b[c].astype(str))]
    bad += bool(diff)
    lines.append(f"{f}: {len(a)} retrievals, {'IDENTICAL' if not diff else 'DIFFER in ' + ','.join(diff)}")
lines.append("")
lines.append(f"Result: {'all compared files reproduce exactly' if not bad else str(bad) + ' file(s) differ'}.")
open(out, "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
PY
echo CLEAN_REPRO_DONE
