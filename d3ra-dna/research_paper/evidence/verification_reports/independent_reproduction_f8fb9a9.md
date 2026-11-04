> Note (added at packaging): this independent clean-environment reproduction was run at commit f8fb9a9, **before** the paper experiments E4/EC/ED/E10, the overdispersion option and the v2 policy re-runs. It does not cover them. The missing `edlib` requirement it reports has since been added to `requirements.txt`.

# D3RA-DNA independent reproduction report

Date: 2026-11-02. The original working copy at /home/user/Dev was not modified.

## Environment
- Clone: `git clone --branch claude/fervent-faraday-92s7l6 /home/user/Dev repro_clone`. HEAD is `f8fb9a9` ("Q1-Q4 work in progress: ...").
- Fresh venv: `python3 -m venv repro_venv`. Python 3.13.18, Linux, 4 cores. Other heavy jobs were running at the same time (`q4_stopping_safety.py` and `real_data.py prep channel policies --workers 3` from the original copy, load average about 5.5).
- `pip install -r d3ra-dna/requirements.txt` succeeded in 41 s with no failures.
- Full freeze: `scratchpad/pip_freeze.txt` (before edlib was added). Key versions: numpy 2.5.4, scipy 1.18.1, pandas 3.1.0, matplotlib 3.11.2, networkx 3.7, pyarrow 26.0.0, reedsolo 1.7.0, PyYAML 6.0.3, pytest 9.2.0, hypothesis 6.170.14, streamlit 1.66.1, tabulate 0.10.0. edlib 1.3.9.post1 was installed afterwards; see below.

## 1. Test suite (`python -m pytest -q`)
- **With requirements.txt only: FAILS at collection.** `tests/test_realdata.py` imports `d3ra.channel.empirical_errors`, which imports `d3ra.realdata.mapping`, which runs `import edlib`, and that raises `ModuleNotFoundError: No module named 'edlib'`. pytest stops with "Interrupted: 1 error during collection", so 0 tests run.
- With `--ignore=tests/test_realdata.py`: **189 passed** in 7.3 s.
- After `pip install edlib`: **201 passed** in 11.1 s, with 0 failures and 0 skips.
- `edlib` is imported in `src/d3ra/realdata/mapping.py`, `src/d3ra/channel/empirical_errors.py`, `src/d3ra/channel/correlated_errors.py` and `experiments/real_data.py`. It is declared in neither `requirements.txt` nor `pyproject.toml`. It is the only undeclared third-party import found; I scanned every import in src/, experiments/, tests/ and dashboard/.
- Cosmetic issue: `addopts = "-q"` in pyproject plus the README's `-q` gives `-qq`, which suppresses the "N passed" summary line. A new user sees only dots.

## 2. `python experiments/reproduce.py`
```
[PASS] Theorem 1 sign prediction (D | N): 8148/8148  (claim: 100% agreement)
[PASS] analytic predictor vs Monte Carlo: 0.0044  (claim: max |dP| < 0.02)
[PASS] spreading threshold direction (s4): {'domain': {12: 0.1453, 24: 0.0113}, 'failure0': {12: 0.0347, 24: 0.0113}, 'round_robin': {12: 0.1427, 24: 0.0103}, 'tiered1': {12: 0.174, 24: 0.061}}  (claim: parity 12: concentration safer than spreading; parity 24: spreading <3% loss)
[PASS] s1 adaptive stopping saving vs fixed depth: 0.439  (claim: ~0.48)
4/4 checks passed
real 2m6s, exit 0
```
- **4/4 PASS.** A second run gave identical output, and `results/reproduce_report.json` was byte-identical, so the script is deterministic.
- Runtime was about 2 min, not the "~15-25 min" in its docstring. That is good, but the docstring is stale.
- Note: the s1 saving measured 0.439, while the claim is ~0.48 and the README headline is "48–79%". It passes only because the tolerance band is wide (0.3–0.65), and it uses reduced trials (600 adaptive / 6000 static). A reader should not treat this check as confirming 0.48 exactly.
- `--real` checks were not run. They need the ~0.5 GB ENA/Zenodo download plus `real_data.py prep channel`.
- The script writes an untracked `results/reproduce_report.json` into the tree.

## 3. `python experiments/theory_checks.py` vs committed `results/theory_summary.json`
- Ran in 4 s. The regenerated `theory_summary.json` is **byte-identical** to the committed one: `git status` showed no change, and a sorted-key JSON diff was empty.
- All regenerated `results/theory_*.csv` (qstar, homogeneous, heterogeneous, nonmds, dropout_margin) are byte-identical to the committed files.
- Values: T1_sign_agreement 0.99554 (1344 cases), closed-form max err 4.0e-15, MC max err 0.00928, P2 violations 0/1680, T5 share 0.911. This matches `docs/theory_spreading.md`, which documents the 6 mismatches as uneven splits (D does not divide N).
- **No discrepancies.**

## 4. README quick-start check
| Command | Result |
|---|---|
| `pip install -r requirements.txt` | Works, but is **missing edlib**. Also `streamlit>=1.30` is listed twice, which is harmless. |
| `python -m pytest -q` | **Fails** without edlib (see 1). |
| `PYTHONPATH=src python -m d3ra.cli README.md --containers 6` | **Wrong command.** It errors with `invalid choice: 'README.md' (choose from 'store', 'simulate', 'retrieve', 'demo')`. The `demo` subcommand is missing. |
| `PYTHONPATH=src python -m d3ra.cli demo README.md --containers 4` | Works in 8.7 s: `exact 900.0 {... 'code_rate': 0.75}`. But (a) README step 3 says "Recovered files are written to `./recovered`", and `demo` writes no files; it only prints status. (b) **Side effect: it appends a new entry to the tracked file `data/processed/phi_calibration.json`** (a calibration cache keyed by parameter hash), which leaves the working tree dirty. I reverted this in the clone. |
| `store my.txt --archive arch` / `simulate --archive arch --out fq2` / `retrieve --archive arch --fastq fq2 --out recovered` | Works; a 3000-byte file was recovered **byte-exact**. |
| `simulate ... --transplant illumina` | **Fails on a fresh clone**: `data/processed/realdata/transplant_illumina.pkl missing: run experiments/real_data.py transplant first`. The pickle is gitignored. The README lists `real_data.py` as the *last* command, after this one, so the order is wrong. The README `retrieve --fastq fq` step then has no input. |
| `python -m d3ra.visualization.plots results <out>` | Works in 5 s and produces 9 PNGs. Committed `reports/figures/` has 12; the other 3 come from `real_data.py plots`. |
| `run_all.py --quick` / full, `streamlit run dashboard/app.py`, `real_data.py download ...` | Not run: they take too long or need network/data. Entry points exist and `--help` works. |

- The README does not mention `experiments/reproduce.py` or `experiments/theory_checks.py`, even though these are the fastest way to check the headline claims.
- PYTHONPATH=src is needed because the package is not installed. `pip install -e .` would also work through pyproject, but the README does not say so.

## Reproducibility problems found
1. **Undeclared dependency `edlib`.** The test suite fails at collection on a clean install.
2. **README demo command is missing the `demo` subcommand**, and its description (files written to ./recovered) does not match what `demo` does.
3. **Order problem in the README:** `simulate --transplant illumina` needs gitignored data from `real_data.py transplant`, which is listed afterwards.
4. **The demo/pipeline writes to the tracked cache `data/processed/phi_calibration.json`.** Runs mutate the committed state, and results may depend on which cache entries already exist.
5. Minor: the reproduce.py docstring runtime is stale (2 min, not 15–25); the s1 check tolerance is loose relative to the claim (0.439 measured vs ~0.48); `-q` twice hides the pytest summary; streamlit is duplicated in requirements.txt.
6. No nondeterminism was observed. reproduce.py (two runs) and theory_checks.py reproduce exactly.

Logs are in the scratchpad: pip_install.log, pip_freeze.txt, pytest.log, pytest_noedlib.log, pytest_edlib.log, reproduce.log, reproduce_run2.log, theory_checks.log, theory_summary_committed.json, demo.log, plots.log, phi_after_demo.json, cli_work/.
