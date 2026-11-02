# D3RA-DNA: Distributed DNA Reliability and Retrieval Architecture

**Reliability-aware coding, fragment placement and adaptive sequencing-budget allocation for distributed DNA data storage under heterogeneous, non-uniform and noisy channels.**

This is a software-only research prototype. All synthesis, storage and sequencing behaviour is simulated. Each parameter is tagged with its provenance in `docs/channel_calibration.md`: measured, estimated from published evidence, or hypothetical. See *Limitations* below.

## What is here

| Layer | Modules (`src/d3ra/`) | Highlights |
|---|---|---|
| Ingestion | `ingestion/` | Arbitrary binaries; PDF/PNG/JPEG/ZIP/text sniffing; multi-file versioned bundles; optional zlib/bz2/lzma (`auto` picks the smallest); SHA-256 |
| Encoding | `encoding/` | 2-bit, scrambled 2-bit and Goldman-style rotating ternary codecs; GC, homopolymer and windowed-GC validation |
| Fragmentation | `fragmentation/` | Compact physical addressing kept separate from logical manifest metadata; overhead trade-off analysis |
| Coding | `coding/` | Vectorised GF(256) Cauchy Reed–Solomon (outer, MDS); `reedsolo` inner RS + CRC; LT fountain; replication; per-stripe redundancy plans |
| Storage | `storage/` | Container objects with capacity, inventory, error/dropout/failure/cost and history; six placement strategies; library builder; manifest; SQLite results store |
| Channel | `channel/` | Abundance (shared sequence-intrinsic + container-specific, gamma/lognormal, PCR); Poisson/multinomial/NB sampling; context- and position-dependent substitution/insertion/deletion; bursts; dropout; correlated failure domains; **recovery-time-matrix surrogate** and analytical predictor; calibration of φ(n) from the sequence-level simulator |
| Retrieval | `retrieval/` | Inventory-aware clustering; two-sided BMA consensus; fixed and adaptive (uniform / coverage / pivotal / Thompson / oracle / clairvoyant) controllers; stopping rules; sequence-level environment |
| Optimisation | `optimization/` | Placement: MILP (HiGHS), greedy, local search, SA, GA. Budget: exponentiated gradient + bisection. Redundancy: uniform vs risk-adaptive. Joint: alternating optimisation |
| Decoding | `decoding/` | Strand decoding; outer erasure decoding; reassembly; decompression; verification into exact / partial / failed / integrity-mismatch |
| Evaluation | `evaluation/` | Scenarios 1–10; baselines A–F; H1–H5; ablations A–G; bootstrap CIs, McNemar, Wilcoxon, Cliff's δ, Holm |
| Visualisation | `visualization/`, `dashboard/` | Matplotlib figures, markdown tables, Streamlit dashboard |

Documents:
- `docs/system_model.md`: the mathematical specification.
- `docs/literature_review.md`: related work and where this project's novelty claims stand.
- `docs/channel_calibration.md`: parameter provenance and public datasets.
- `reports/technical_report.md`: the technical report.
- `reports/paper_draft.md`: a paper draft.

## Quick start

```bash
pip install -r requirements.txt
```

```bash
python -m pytest -q
```

```bash
PYTHONPATH=src python -m d3ra.cli demo README.md --containers 4
```

```bash
python experiments/run_all.py --quick
```

```bash
python experiments/run_all.py
```

```bash
PYTHONPATH=src python -m d3ra.visualization.plots results reports/figures
```

```bash
streamlit run dashboard/app.py
```

```bash
PYTHONPATH=src python -m d3ra.cli store my.pdf --archive arch
```

```bash
PYTHONPATH=src python -m d3ra.cli simulate --archive arch --out fq
```

Real-error transplant (`--transplant illumina` or `--transplant nanopore`) requires first running `python experiments/real_data.py download transplant` (nanopore: `nanodiag`).

```bash
PYTHONPATH=src python -m d3ra.cli retrieve --archive arch --fastq fq --out recovered
```

```bash
python experiments/real_data.py download prep transplant nanopore nanodiag channel policies failures plots
```

What each command does, in order:
1. Installs dependencies.
2. Runs the test suite.
3. Stores and retrieves a real file end-to-end through simulated containers, in memory, and prints the status. To write recovered files to disk, use the `store` / `simulate` / `retrieve` workflow below.
4. Smoke-runs every experiment, about 10–20 minutes.
5. Runs the full programme, which writes `results/*.csv` (about 1–2 h on a laptop CPU).
6. Renders the figures.
7. Starts the dashboard.

The last four commands are the on-disk workflow: store the archive, write FASTQ files with real error profiles, retrieve by streaming FASTQ adaptively, then run the real-data validation (downloads about 0.5 GB from ENA and Zenodo).

Everything runs on a CPU. No GPU is needed.

## Verifying the results

```bash
python experiments/reproduce.py
```

```bash
python experiments/theory_checks.py
```

`reproduce.py` re-derives the headline claims with tolerances in about 2–5 minutes; add `--real` once the real data has been prepared. `theory_checks.py` regenerates `results/theory_*`, which should match the committed files exactly.

## Core idea in one paragraph

Each sequencing read from container j is treated as an event on that container's own read clock. A (fragment, container) cluster becomes decodable after t_ij reads. t_ij is drawn from a consensus-success curve, and that curve is calibrated by running real BMA consensus and the inner RS decoder on simulated noisy reads. So the read count at which fragment i becomes recoverable from container j is a random variable T_ij, and one sampled T-matrix determines how *every* policy performs on that trial:
- For static allocations, the minimal budget for exact recovery is closed-form per trial.
- Adaptive controllers are replayed on the same T, which gives common random numbers across policies.

This makes 10⁴–10⁵-trial reliability estimates (targets 0.99 and 0.999) cheap. An independent analytical predictor (negative-binomial counts, Gauss–Hermite quadrature over the shared bias, exact Poisson-binomial tails, failure-state averaging) matches the Monte Carlo to within about 0.002–0.02 in success probability.

## Results

All results are in simulation, under the documented assumptions. Full tables, CIs, negative results and failure cases are in `reports/technical_report.md`.

- **Adaptive stopping vs fixed depth.** Stopping once every stripe holds K verified fragments cuts expected sequencing spend at 0.99 reliability by 48–79% compared with the fixed depth that reaches 0.99.
- **Placement must be optimised with the adaptive reader in the loop. It is often *not* distributed.** When one container meets the failure budget, putting each stripe there is 22–37% cheaper than reliability-weighted spreading. With correlated failure domains, the search spreads over five containers and is the only design to reach 0.99. Without adaptive allocation, that design costs +79%.
- **Spreading threshold.** Spreading stripes over D failure domains helps only once the parity fraction exceeds about 1/D. Below it, concentration was 3.9× safer.
- **Channel models (H1).** Container-specific models predict recovery with mean absolute error 0.007, versus 0.073 for homogeneous models.
- **Robustness (H5).** Wide Bayesian priors match oracle completion under post-placement channel change (97.6% vs 74% for an over-confident controller).
- **Real data** (ENA PRJEB65931 Illumina runs; Microsoft nanopore clusters; `reports/real_data_validation.md`):
  - The real coverage bias is 69% shared across containers.
  - Independent-error models overstate nanopore consensus by 20–40 points, because of cross-read correlated errors.
  - On real read streams, placement plus adaptive stopping saves 34% (stopping alone 15–17%).
  - After a container loss, cost is governed by the post-failure margin.
  - Real files are recovered byte-exact from FASTQ with real error profiles.
- **Negative results.** Thompson sampling; stripe-importance weighting vs plain hazard-greedy allocation; risk-adaptive parity (predictor tail too optimistic); shared-bias vs replication; no design meets 0.99 under the compound-stress scenario.

## Limitations

- No wet-lab validation.
- Container failure probabilities and domain structures are hypothetical.
- Per-container consensus: reads of the same fragment are not pooled across containers.
- Costs are in abstract units.

None of the numbers here should be read as real-world guarantees.
