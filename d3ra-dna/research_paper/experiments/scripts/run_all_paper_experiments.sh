#!/usr/bin/env bash
# Re-runs every experiment behind the manuscript (several hours on 4 cores). Run from d3ra-dna/.
set -euo pipefail
python experiments/real_data.py download prep channel
D3RA_TAG=_v2 python experiments/real_data.py policies failures
D3RA_POOL=genscript python experiments/real_data.py prep channel
D3RA_POOL=genscript D3RA_SKIP_JOINT=1 D3RA_VARIANTS=round_robin D3RA_TAG=_rr_v2 python experiments/real_data.py policies
D3RA_POOL=genscript D3RA_SKIP_JOINT=1 D3RA_VARIANTS=single_best D3RA_TRIALS=100 D3RA_TAG=_sb_v2 python experiments/real_data.py policies
python experiments/run_all.py
python experiments/theory_checks.py
python experiments/q3_correlated_errors.py
python experiments/q3b_reference_audit.py
python experiments/q4_stopping_safety.py
python experiments/q5_seed_sensitivity.py
python experiments/paper_experiments.py e4 ec ed
PYTHONHASHSEED=0 python experiments/paper_experiments.py ed_od
D3RA_POOL=genscript PYTHONHASHSEED=0 python experiments/paper_experiments.py ed_holdout
D3RA_POOL=genscript python experiments/paper_experiments.py e10 e10_static
python research_paper/experiments/statistical_analysis/ec_discrimination.py
python research_paper/experiments/statistical_analysis/ed_overdispersion_summary.py
python research_paper/experiments/scripts/make_figures.py
(cd research_paper/manuscript && latexmk -pdf main.tex)
