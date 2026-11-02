# D3RA-DNA: Answers to the Six Critical Research Questions

Every statement below is backed by a script and a result file in this repository. Negative results and corrections are kept. Where an answer weakens an earlier claim, the earlier claim has been corrected in the other documents.

| Q | Short answer |
|---|---|
| 1 Container validity | The six runs are **aliquots of one synthesis batch on two flowcells**, not independent containers. Handling explains about 5× the variation of a technical replicate. Conclusions transfer to *same-batch aliquot* containers, not to independently synthesised ones. Real failure rates are unmeasured. |
| 2 Spreading threshold | Formally proved (Theorem 1, Proposition 2) and checked to 100% where the theorem applies. It is **a special case of known allocation theory** (Leong, Dimakis and Ho 2012), not a new theorem. |
| 3 Correlated nanopore errors | **Real, but smaller than first reported.** Persistent across two flowcells, mostly strand-specific beyond 5-mer context, 6–15 points on verified references. A hierarchical model cuts prediction error 2.2×. |
| 4 Stopping safety | Verification-gated stopping **cannot output a wrong file** except through a SHA-256 collision. Zero false accepts in about 10,500 decode attempts with the default code. The futility (give-up) rule is **harmful** and should be off. |
| 5 Reproducibility | A clean-environment reproduction passed after fixing four reproducibility defects. Generalisation to an independent pool and seed sensitivity: §5. |
| 6 Novelty | Most components are known. The strongest defensible contribution is a **real-data methods study of placement plus adaptive stopping across heterogeneous containers**, with honest effect sizes. |

---

## Q1. Are the six Illumina runs valid stand-ins for containers?

Scripts: `experiments/q1_container_validity.py`. Results: `results/q1_pairwise.csv`, `q1_error_rates.csv`, `q1_summary.json`.

**What the runs physically are.** The FASTQ headers show:
- **Aging series** (0 d replicate A, 0 d replicate B, 2 d, 7 d): instrument FS10000169, run 123, **one shared flowcell**, separated only by index.
- **PCR series** (15, 60, 90 cycles): run 121, a second flowcell.

All samples come from **one Twist synthesis batch**. The PCR samples are serial re-amplifications, and the aged samples are aliquots of the same pool.

**Variance decomposition.** Excess log-coverage variance beyond Poisson, between pairs of samples:

| comparison | excess variance | log-coverage r | Δ substitution rate |
|---|---|---|---|
| technical replicate (same handling, same flowcell) | **0.011** | 0.80 | 0.07%/nt |
| same flowcell, different handling | 0.052 (0.009–0.190) | 0.72 | 0.25%/nt |
| different flowcell (and different handling) | 0.061 | 0.68 | 0.33%/nt |

**Answers:**
- **Independent physical containers? No.** They are physically separate aliquots with different histories, so they are real *containers* in the sense of separately stored samples. But they share synthesis, so their coverage is correlated: 69% of variance is shared, which is an upper bound for same-batch aliquots.
- **Run vs handling.** Handling adds about 5× the variation of a technical replicate (0.052 vs 0.011). The flowcell adds a little more (0.061), but **the PCR series is entirely on one flowcell**, so PCR treatment and sequencing run are partly confounded. Error-rate differences follow the same pattern.
- **Generalisation.** The conclusions cover containers made by *aliquoting one synthesis batch*, which is a realistic replication strategy. They do **not** establish behaviour for containers synthesised independently (other chips, vendors or batches). There the shared component should be smaller, which helps replication-style redundancy and changes the placement trade-offs.
  - The Genscript pool (§5) is an independent synthesis but again *one batch per pool*.
- **Validating real container failures.** No public data exists. A credible experiment:
  1. Synthesise ≥ 2 independent batches, with ≥ 2 vendors. Aliquot into vials grouped into storage domains (different freezers, buffers, sites).
  2. Apply controlled stresses to subsets (accelerated aging at 60–70 °C as in Grass 2015, freeze–thaw cycles, humidity) to raise the failure probability into a measurable range.
  3. Sequence every vial at staged depths, defining "failure" as loss of decodability at a fixed depth.
  4. **Power:** estimating a per-vial failure probability f with relative precision r at 95% needs n ≈ 1.96²(1−f)/(f r²) vials. For f = 1% and r = 0.5 that is ≈ 1,520 vials, which is why accelerated stress is required. At f = 10% under stress, ≈ 140 vials.
  5. Domain correlation would be estimated by pairing vials within and across domains.

## Q2. Mathematical derivation of the spreading threshold

Proof: `docs/theory_spreading.md`. Code: `src/d3ra/theory/spreading.py`. Checks: `experiments/theory_checks.py`, giving `results/theory_*.csv` and `results/theory_summary.json`.

- **Theorem 1 (homogeneous domains, D | N, MDS code, independent domain failures with probability q).** Let ρ = (N−K)/N and t = ⌊Dρ⌋.
  - (a) If ρ < 1/D, even spreading is strictly worse than concentration for **every** q ∈ (0,1).
  - (b) If t ≥ 1, there is a unique q* with spreading better ⟺ q < q*. The proof uses the running mean of a unimodal Beta density.
  - Example: q* = 0.232 for D = 4, t = 1, and q* = 0.042 for D = 8, t = 1.
- **Checks:**
  - The closed form matches exact enumeration to 4×10⁻¹⁵.
  - The sign prediction is correct in **100%** of covered cases (1,338 here; 8,148 in the independent reproduction).
  - The only mismatches (6) are uneven splits, outside the theorem's assumptions.
- **Dropout margin.** Fragment dropout δ raises the threshold to about 1/D + δ + a dispersion term. With D = 4 and δ = 2%, the exact crossover is ρ = 0.273. Overdispersed simulated coverage moved it to about 0.33.
- **Other codes.** Fountain codes needing K(1+ε) shift the threshold to 0.33 at ε = 0.05–0.1 and to 0.38 at ε = 0.2.
- **Heterogeneous domains and capacities.**
  - **Proposition 2:** a placement without single-domain tolerance can never beat concentration in the most reliable domain. 0 violations in 1,680 random instances.
  - **Corollary 3:** single-domain tolerance requires Σ_d min(C_d, N−K) ≥ N.
  - In 200 random heterogeneous instances, the best tolerant placement won in 91% where both options were feasible. The losses were at low ρ (median 0.19).
- **Assumptions.** Independent domain failures, MDS (or fountain with overhead), and a single stripe. The multi-stripe file level has the same sign but different magnitudes.
- **Novelty.** Not new. It is a special case of the symmetric-allocation analysis of Leong, Dimakis and Ho (IEEE Trans. IT 2012), comparing m = D with m = 1. It is presented as a design rule.

## Q3. Correlated nanopore sequencing errors

Scripts: `experiments/q3_correlated_errors.py`, `experiments/q3b_reference_audit.py`. Model: `src/d3ra/channel/correlated_errors.py`. Results: `results/q3_summary.json`, `q3_model_validation.csv`, `q3b_reference_audit.json`.

The data is Zenodo 17473983 (DNAformer): 55,000 strands, **two separate flowcells** (634k and 453k reads).

- **Is it real, or an artefact of wrong references? Real, but smaller than first reported.**
  - For 86% of 2,393 strands, the published design equals the consensus of the independent flowcell's reads.
  - Wrong-reference signature positions (most reads deviate) are 0.1% of positions, and excluding them leaves the cross-flowcell correlation unchanged (0.47 → 0.45).
  - On the 2,057 verified strands, a conservative selection, real 5-read and 10-read exact-consensus rates are 0.729 and 0.937. Independent errors give 0.880 and 0.997; independently replayed real edit scripts give 0.856 and 0.996.
  - The earlier "20–40 point" gap (Microsoft dataset) is partly explained by that dataset's documented reference problems. **The defensible magnitude is 6–15 points.**
- **Why do errors repeat? Mostly strand-specific, persistent across runs.**
  - Per-(strand, position) error rates correlate between flowcells at r = 0.60. The 5-mer + homopolymer context mean explains only r = 0.29.
  - The *residual* (strand-specific) component correlates at r = 0.57 across flowcells, against r = 0.002 for a random strand with the same context.
  - A failure on flowcell 1 raises the failure probability on flowcell 2 from 17% to 45%.
  - Intra-class correlation of errors across reads: 0.07 on flowcell 1 and 0.12 on flowcell 2.
  - Short-context basecalling effects therefore explain only part. The rest is tied to the specific molecule or strand: synthesis variants in part of the molecules, longer-range basecaller context, or residual clustering errors. **This data cannot separate those**; resolving it needs a re-basecalled or orthogonal-platform replicate.
- **Effect on retrieval.** Correlated errors put a floor under consensus success: about 6% of verified strands still fail at 10 reads, against ≤ 0.4% for independent errors. These failures behave as erasures, and the outer code must budget parity for them. The default K = 32 of N = 48 (33% parity) absorbs them in the end-to-end runs.
- **Correlation-aware model.** Hierarchical beta-binomial: context mean × strand-position propensity. Fitted on flowcell 1 and evaluated on held-out strands of flowcell 2:

  | model | mean absolute error |
  |---|---|
  | independent errors | 0.181 |
  | context-only | 0.153 |
  | **hierarchical** | **0.083** (2.2× better than independent) |

  Remaining error: flowcell 2 is harder (17% vs 9% consensus failure; ICC 0.12 vs 0.07), so a model fitted on one run transfers only partly, and run-level calibration is needed.

## Q4. Mathematical safety of adaptive stopping

Script: `experiments/q4_stopping_safety.py`. Code: `decoding/reconstruction.py`, `retrieval/sequencing_env.py`. Results: `results/q4_*.csv`, `q4_summary.json`.

**Definition: provable decodability.** Every stripe holds ≥ K fragments, each of which passed (i) inner Reed–Solomon decoding, (ii) the CRC-c check and (iii) the address check byte. In addition, the decoded stripe matches the per-stripe CRC-32 (logical metadata in the manifest), **and** the assembled file matches its SHA-256.

**Proposition (verification-gated stopping).** The controller stops with "exact" only after the SHA-256 check passes. So P(output is wrong) ≤ 2⁻²⁵⁶ per decode attempt, a hash collision, independent of sequencing quality.
- If a wrong fragment passes the inner checks (a *false accept*), it is caught by the per-stripe CRC or the file hash. Retrieval then continues: the controller demands K+1 fragments for the flagged stripe, and the corrupted fragment is found by leave-one-out decoding.
- **Premature stopping, meaning stopping with a wrong file, is therefore impossible apart from that collision.** The only way to stop without success is to exhaust the budget, or the futility rule (below).

**Per-attempt false-accept probability.** For a random wrong word: ≈ V_t(n)/256^nsym × 2^{−8c}, where V_t is the radius-t Hamming ball, t = nsym/2 and c is the number of CRC bytes. Measured on realistic wrong consensus words (3-read nanopore-error clusters):

| inner code | false accepts per wrong word | analytic estimate |
|---|---|---|
| none | 0.98 | 1.0 |
| RS nsym 2, no CRC | 0.122 | 0.152 |
| RS nsym 4, no CRC | 0.0074 | 0.012 |
| RS nsym 2, CRC-8 | 0 / 2,746 (≤ 1.1×10⁻³) | 6×10⁻⁴ |
| **RS nsym 4, CRC-16 (default)** | **0 / 2,682 (≤ 1.1×10⁻³)** | **2×10⁻⁷** |

**End-to-end retrievals** (12 each, heavy errors):
- **Default code:** 12/12 exact, **0 false accepts in 10,567 decode attempts**, 0 wrong outputs.
- **Weak code (nsym 2, no CRC):** 54–60 false accepts per retrieval. All were caught: **0 wrong files output**. Retrievals ended *partial* because bad fragments outnumbered what leave-two-out repair can isolate. The weak code is safe but not live, which justifies the CRC.

**Futility rule, a negative result.** Giving up when the estimated completion probability drops below ε (10⁻³–5×10⁻²) abandoned 3.7–4.5% of retrievals in the channel-change scenario. **All of them would have succeeded** within the budget. The completion-probability model is miscalibrated under channel change, so the rule is off by default and is not recommended without calibration.

**Remaining assumption.** The manifest (per-stripe CRCs, SHA-256, code parameters) must itself survive. It is small and should be stored with heavy replication. This is standard but outside this analysis.

## Q5. Reproducibility and generalisation

See §5 below, which is filled from `results/q5_*` and `results/realdata_genscript_*`.

## Q6. Novelty and scientific contribution

Full verified assessment: `docs/novelty_assessment.md`.

| component | status | closest prior work |
|---|---|---|
| spreading threshold | **known** (special case) | Leong, Dimakis and Ho 2012; Nath et al. 2006 |
| post-failure margin rule | **known** (corollary) | Bar-Lev et al. (E = n(H_n − H_{n−k})); Levy et al. 2026 |
| multi-container DNA storage with failures | **known setting** | Levy, Con, Yaakobi and Kiah 2026; Kong et al. 2026 (uniform, no failures) |
| stop when decodable | **known practice** | Lopez et al. 2019; Zhang et al. 2026; readfish |
| verification-gated stopping | **standard** | inner RS/CRC rejection in DNA Fountain and others |
| choosing *which container* to read to complete K-of-N stripes | **apparently new for DNA storage** | Good-UCB (Bubeck et al. 2013); BOSS-RUNS (2023) |
| systematic nanopore errors | **known qualitatively** | Delahaye and Nicolas 2021; Wick et al. 2019; Hamoum et al. |
| cross-flowcell persistence of strand-specific errors, and a validated hierarchical model for consensus prediction | **partially new** (quantification) | — |
| real-read evaluation of placement × adaptive stopping across heterogeneous containers | **apparently new** | — |

**Strongest defensible contribution:** a **methods paper with real-data evaluation**. It would contain:
- placement plus adaptive retrieval evaluated on real recovery-time matrices from two independent pools;
- honest effect sizes: −34% reads on real Twist data, adaptive stopping alone 15–17%, against 48–79% in simulation;
- the model–reality gap (surrogate 5–18% optimistic);
- container-loss behaviour on real streams;
- the quantified, validated correlated-error model as the supporting channel result.

The combination alone is a *systems/methods* contribution, not a theoretical one. The only algorithmic element without clear precedent in DNA storage, container selection, gave the smallest measured gain (≤ 10%, and pivotal scoring tied plain coverage-greedy). A theory-venue paper would need a new result, for example optimality or regret bounds for adaptive container selection with heterogeneous efficiencies. That is not established here.

**Claims withdrawn or reworded:** the spreading threshold and margin rule as new results; "sequence-intrinsic" shared bias, now "batch-inherited, upper bound"; the 20–40-point nanopore gap, now 6–15 on verified references; Kong et al. as the closest prior work, Levy et al. added; "first multi-container" (never claimed, explicitly excluded).
