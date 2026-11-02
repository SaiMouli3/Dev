# D3RA-DNA — Validation on Real DNA-Storage Sequencing Data

This report covers what changed when the D3RA-DNA retrieval stack and models met real reads, including the bugs that real data exposed.

Reproduce everything with:

```bash
python experiments/real_data.py download prep transplant nanopore nanodiag channel policies failures plots
```

Raw data (~0.5 GB) is not committed. Derived results are in `results/realdata_*`.

## 1. Data

| Dataset | What it is | Ground truth | Use here |
|---|---|---|---|
| ENA **PRJEB65931** (Gimpel et al., *Nat. Commun.* 2023, DT4DDS), 6 runs | Illumina iSeq 100 reads of **one 12,000-strand Twist pool** (108 nt) after different handling: aging 0 / 2 / 7 days, PCR 15 / 60 / 90 cycles. 1.5–2.5 M reads per run | Designed sequences from the DT4DDS release (Zenodo 12666365, `Twist_GCall/design_files.fasta`) | Error / coverage / dropout estimation; **six physically distinct real "containers"** holding the same strands |
| Zenodo **14296588** (Srinivasavaradhan et al., clustered nanopore reads) | 9,954 clusters, 166,296 nanopore reads, 110-nt strands | Reference strand per cluster | Real consensus curves under nanopore noise |

**Why one pool with different handling is a valid multi-container testbed.** Reads are independent draws from each sample. A container that holds only a subset of strands therefore produces exactly the run's read stream restricted to that subset, plus its share of junk reads. We use this restriction to emulate any placement on real reads (`realdata/empirical.py`). One caveat: no real container *fails* in this data, so failures are injected (§6).

## 2. Pipeline on real reads

1. **Stream FASTQ.** `io/fastq.py`.
2. **Trim the adapter.** About 99.5% of reads map to a designed strand: k-mer voting, then edlib prefix-global verification with relative edit distance ≤ 0.25.
3. **Cluster per strand, then take a consensus.** Exact success means the consensus equals the design.
4. **Run everything with a held-out split.** The first 300k reads of each run are **calibration** and estimate the controller's beliefs. The next 900k are **evaluation** reads and are never seen by the beliefs. The one exception: the pooled consensus curve in the belief uses even-indexed strands of the evaluation split, while evaluation uses all strands. This is a small leak, noted here.

### A bug found only on real data

Two-sided bitwise-majority-alignment (BMA) consensus combines a forward pass with a reversed pass. It assumes every read ends exactly at the strand end. Raw iSeq reads continue past the 108-nt insert, mostly as read-through noise and poly-T. Only about 3% end in a clean adapter. The reversed pass therefore anchors on garbage.

| Consensus | 1 read | 2 reads | 3 reads | 5 reads | 8 reads |
|---|---|---|---|---|---|
| two-sided (original) | 0.52 | 0.05 | 0.01 | 0.02 | 0.00 |
| **forward-only (fix)** | 0.52 | 0.70 | 0.97 | 0.995 | 1.00 |

The fix is `consensus(..., two_sided=False)` for raw reads, plus a `--consensus forward` option in the CLI retriever. Strand-trimmed reads, such as our own strands or the nanopore set, keep the two-sided mode, which is better there (§4). The invalid first channel table is kept as `results/realdata_channel_INVALID_two_sided_consensus.csv`.

## 3. Real Illumina channel vs. model assumptions

From `results/realdata_channel.csv` and `results/realdata_shared_variance.json`, with about 25× calibration coverage per strand.

**Error rates and coverage:**

| condition | sub /nt | del /nt | ins /nt | error-free reads | coverage CV | NB shape | dropout (≈100×) |
|---|---|---|---|---|---|---|---|
| PCR 15 c | 0.35% | 0.06% | 0.01% | 67.7% | 0.35 | 12.4 | 0.01% |
| Aging 0 d (fresh) | 0.60% | 0.07% | 0.01% | 52.4% | 0.35 | 12.0 | 0.08% |
| Aging 2 d | 0.71% | 0.08% | 0.01% | 46.1% | 0.37 | 10.7 | 0.13% |
| Aging 7 d | 0.74% | 0.09% | 0.01% | 43.3% | 0.37 | 10.1 | 0.09% |
| PCR 60 c | 0.87% | 0.09% | 0.01% | 37.7% | 0.39 | 8.9 | 0.23% |
| PCR 90 c | 1.15% | 0.11% | 0.01% | 27.4% | 0.42 | 7.1 | 0.43% |

**Consensus success and unrecoverable strands:**

| condition | P(exact) 1 / 3 / 5 reads | never exact (16 reads) |
|---|---|---|
| PCR 15 c | 0.68 / 0.987 / 0.999 | 0.07% |
| Aging 0 d (fresh) | 0.52 / 0.977 / 0.997 | 0.30% |
| Aging 2 d | 0.46 / 0.965 / 0.996 | 0.42% |
| Aging 7 d | 0.43 / 0.962 / 0.996 | 0.39% |
| PCR 60 c | 0.36 / 0.928 / 0.993 | 0.59% |
| PCR 90 c | 0.27 / 0.873 / 0.983 | 1.28% |

Findings:

- **The containers really are heterogeneous.** Error rates differ 3.3× between the best and worst samples. The over-amplified sample needs about 2× more reads per strand, and has 4.3× the dropout and 18× the never-recoverable strands.
- **Errors are dominated by substitutions.** My literature-based default (sub 0.4%, del 0.2%, ins 0.04%; Organick 2018 profile) overstated deletions 2–3× and insertions about 4× for this platform.
  - Deletions are 1.1–1.5× more likely inside homopolymers, consistent with the model's homopolymer boost.
  - Per-read error counts are over-dispersed: variance 1.2–2.1 against a mean of 0.46–1.37. This supports per-read quality variation (`quality_sigma`) rather than a pure Poisson model.
- **A shared, sequence-intrinsic abundance component exists and is large.** After removing sampling noise, **69%** of log-coverage variance is shared by a strand across all six containers, and 31% is container-specific (mean pairwise log-coverage correlation 0.54).
  - The shared standard deviation is σ_g = 0.31, matching the model default of 0.30 taken from the literature.
  - This directly supports C1 (correlated coverage across replicas or containers from one synthesis batch). It also confirms the warning in the literature review: the correlation is a property of the synthesis batch.
- **Coverage barely depends on GC.** Spearman ρ between GC and coverage is −0.19 to +0.03, consistent with Chen et al. 2020: PCR bias is stochastic rather than GC-driven.
- **The φ(n) curves match the simulator.** At error scale ×1 the simulated curve is 0.55 / 0.98 / 1.00 for 1 / 3 / 4 reads. The real fresh pool gives 0.52 / 0.977 / 0.994.

## 4. Real nanopore consensus vs. the independent-error model

`results/realdata_nanopore_phi.csv`, `results/realdata_nanopore_diagnostics.json`.

- **Measured channel.** Substitutions 1.5%, insertions 1.5%, deletions 2.0% (5.0%/nt in total). 3.2% of reads are error-free. Deletions are 1.6× more likely inside homopolymers.
- **Consensus success on real reads:**

  | reads | real (two-sided BMA) | simulated, independent errors at measured rates | real, position-wise majority vote |
  |---|---|---|---|
  | 3 | 0.17 | 0.29 | 0.01 |
  | 5 | 0.41 | 0.81 | 0.01 |
  | 7 | 0.60 | 0.96 | 0.01 |
  | 10 | 0.75 | 0.99 | 0.01 |
  | 20 | 0.87 | 1.00 | 0.01 |

- **Diagnosis** (500 strands with at least 40 reads, 20 reads per consensus):
  - Real reads give 89.2% exact consensus.
  - Real *edit scripts* from other strands' reads, replayed independently across reads, give **100%**. These keep the real within-read error structure and the real rate.
  - **The whole gap therefore comes from errors correlated across reads of the same strand**, not from the error rate or within-read structure.
  - 3–5% of strands fail with two disjoint read sets of 20 reads. Whether both sets produce the *same* wrong sequence depends on the consensus algorithm: 67% with two-sided, 8% with forward-only. We therefore claim only correlation, not that the molecules truly deviate from the design.
- **Implication.**
  - Independent-error simulators overstate nanopore consensus success by 20–40 percentage points at 5–10 reads, and miss a plateau of about 10% strand failures at 20 reads.
  - For nanopore retrieval, φ must be *calibrated from real reads*, as D3RA-DNA does. The outer code must also budget parity for about 5–10% **strand-systematic erasures** that no extra coverage removes.
  - Position-wise majority voting is useless under nanopore indels (≈1%). Alignment-aware consensus is required.

## 5. Real-read evaluation of placement × sequencing policy

`results/realdata_policies.csv`, `reports/figures/realdata_policies.png`.

**Setup.**
- 12,000 real strands in 250 stripes of (K = 32, N = 48), so any 32 exactly-recovered strands of 48 rebuild a stripe.
- The 6 real runs act as containers. Each Monte-Carlo trial is a random permutation of each container's real evaluation reads; there are 200 trials.
- The recovery time of strand i in container j is the stream position of the read that first makes i's real reads give an exact consensus. The consensus threshold is drawn from that strand's own empirical φ_ij.
- The controller's beliefs come from the calibration split only.

**Reads needed to reach 0.99 reliability.** For fixed depth this is the committed budget; for adaptive policies it is the expected spend. 95% bootstrap CIs in brackets.

| placement (strands per container) | fixed depth (uniform) | adaptive, uniform allocation | adaptive, coverage-greedy | adaptive, pivotal | surrogate prediction (pivotal) |
|---|---|---|---|---|---|
| round-robin (2000 ×6) | 67,938 [64,074; 68,424] | 57,787 | 55,744 | 55,810 | 48,930 |
| reliability-weighted | 65,400 | 56,019 | 54,172 | 54,059 | 48,279 |
| tiered, slack 2 (3500/3500/1500/3500/0/0) | 61,908 (54,518 coverage-equalised) | 51,347 | 46,373 | 46,186 | 43,786 |
| **single best container** (joint optimiser's choice) | 53,867 | **44,753 [44,315; 45,227]** | 44,753 | 44,753 | 36,779 |

Findings on real reads:

1. **Joint placement plus adaptive stopping** needs 44,753 reads against 67,938 for uniform placement with fixed depth: **−34%**.
   - About half of that comes from placement (concentrating in the best container).
   - About half comes from adaptive stopping (−17% at equal placement).
2. **Adaptive stopping is worth 15–17% on real reads, much less than the 48–79% seen in simulation.** The simulated gains were dominated by container failures and heavy tails. These real samples never fail and have mild tails (coverage CV 0.35–0.42). The real benefit is real but scenario-dependent.
3. **Adaptive allocation** helps when data is spread (tiered: 51,347 → 46,186, −10%; round-robin −3.5%) and is irrelevant for one container. **Pivotal and coverage-greedy scoring are indistinguishable**, as in simulation.
4. **The surrogate is 5–18% optimistic on real reads.** It predicts fewer reads than needed: 36.8k vs 44.8k for the single container, 43.8k vs 46.2k for tiered. Model-based budgets need a safety margin. Rankings were preserved.

## 6. Container loss on real read streams

`results/realdata_failure_injection.csv`, `reports/figures/realdata_failure_injection.png`. One container's real read stream is removed (60 permutations), and the adaptive controller learns this from an empty probe.

| placement | no loss | one container lost |
|---|---|---|
| single best container | 44,473 | **total loss** if that container is lost (P = 0) |
| tiered, 4 containers (loses ≤ 14 of 48 per stripe) | 46,530 | always recoverable in principle; after losing a main container only 2 spare fragments per stripe remain, so cost rises to 102–137k and only **45–48%** of retrievals finish within 3× the budget |
| round-robin, 6 containers (loses 8 of 48 per stripe) | 55,610 | +5% to +32% (58–73k); always finishes |
| reliability-weighted | 53,993 | +3% to +41%; always finishes |

**Real-data form of the spreading rule: the post-failure cost is governed by the post-failure margin.** The margin is the number of fragments left per stripe minus K. With 6 spare fragments (round-robin), losing a container costs about +30%. With 2 spare (tiered), it costs about +200% and the retrieval often stalls. With none (single container), the data is gone.

So the cheapest layout is only the right one if container loss is acceptably unlikely, which is exactly what the joint objective encodes through ε. Without measured failure rates for real vials, the choice of ε, not the optimiser, decides between these layouts.

## 7. End-to-end byte-level retrieval with real error profiles

Real files: a 168 KB PNG figure and a 6.6 KB Markdown file, stored as 7,200 strands in 150 stripes (K = 32, N = 48), 168 nt each, across 6 containers.

```bash
d3ra store fig.png readme.md --archive arch --placement failure
```

```bash
d3ra simulate --archive arch --out fq --reads 60000 --transplant illumina|nanopore
```

```bash
d3ra retrieve --archive arch --fastq fq --out rec
```

The `--transplant` option replays real edit scripts:
- *illumina*: 39,867 scripts from ENA ERR12033805, 0.75 errors per 108-nt read;
- *nanopore*: Srinivasavaradhan reads, 5.5 errors per read.

| error source | reads consumed (adaptive stop) | reads per strand | result |
|---|---|---|---|
| real Illumina scripts | 14,832 of 60,000 | 2.1 | **exact** (SHA-256 verified, `cmp` identical) |
| real nanopore scripts | 48,672 of 60,000 | 6.8 | **exact** |

Caveat: transplanted scripts are independent across reads, so the nanopore run does not include the cross-read correlation found in §4. Under real nanopore correlation, expect about 5–10% strand erasures, which K = 32 of N = 48 (33% parity) would absorb.

## 8. What this validation does and does not establish

**Established on real reads:**
- the retrieval stack (mapping/clustering, consensus, the exact-recovery criterion);
- container heterogeneity;
- a shared sequence-intrinsic coverage component (C1);
- the placement × adaptive-stopping gains, with smaller magnitudes than simulated;
- the post-failure-margin rule;
- the inadequacy of independent-error models for nanopore consensus.

**Not established:**
- *Real container-failure statistics.* Failures here are injected.
- *Byte-level decoding of a real synthesized D3RA-format pool.* That needs a wet lab: our encoding has not been synthesized. The real-read experiments use strand-level exact recovery as a proxy for inner-code decoding, which is conservative because an inner RS could fix residual errors.
- *Correlated nanopore errors in the byte-level transplant run.*
