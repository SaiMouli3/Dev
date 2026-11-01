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
PYTHONPATH=src python -m d3ra.cli README.md --containers 6
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

What each command does, in order:
1. Installs dependencies.
2. Runs the test suite.
3. Stores and retrieves a real file end-to-end through simulated containers. Recovered files are written to `./recovered`.
4. Smoke-runs every experiment, about 10–20 minutes.
5. Runs the full programme, which writes `results/*.csv` (about 1–2 h on a laptop CPU).
6. Renders the figures.
7. Starts the dashboard.

Everything runs on a CPU. No GPU is needed.

## Core idea in one paragraph

Each sequencing read from container j is treated as an event on that container's own read clock. A (fragment, container) cluster becomes decodable after t_ij reads. t_ij is drawn from a consensus-success curve, and that curve is calibrated by running real BMA consensus and the inner RS decoder on simulated noisy reads. So the read count at which fragment i becomes recoverable from container j is a random variable T_ij, and one sampled T-matrix determines how *every* policy performs on that trial:
- For static allocations, the minimal budget for exact recovery is closed-form per trial.
- Adaptive controllers are replayed on the same T, which gives common random numbers across policies.

This makes 10⁴–10⁵-trial reliability estimates (targets 0.99 and 0.999) cheap. An independent analytical predictor (negative-binomial counts, Gauss–Hermite quadrature over the shared bias, exact Poisson-binomial tails, failure-state averaging) matches the Monte Carlo to within about 0.002–0.02 in success probability.

## Results

See `reports/technical_report.md` for all numbers, confidence intervals, negative results and failure cases.

## Limitations

- No wet-lab validation.
- Container failure probabilities and domain structures are hypothetical.
- Per-container consensus: reads of the same fragment are not pooled across containers.
- Costs are in abstract units.

None of the numbers here should be read as real-world guarantees.
