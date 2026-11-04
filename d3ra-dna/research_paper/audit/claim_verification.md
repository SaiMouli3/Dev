# Claim Verification Audit (manuscript `main.tex` + `tables/*.tex`)

- **Audited version:** `research_paper/manuscript/main.tex` at commit `89e439f` (756 lines; line numbers below refer to it), together with every table it `\input`s.
- **Evidence used:**
  - `results/*.csv|json`. These are byte-identical to `research_paper/evidence/experiment_results/`, and all SHA-256 sums in `evidence/verification_reports/results_sha256.txt` check OK.
  - `research_paper/evidence/evidence_traceability.csv`.
  - `reports/*.md`.
  - The code in `src/d3ra/` and `experiments/`.
- **What was recomputed:** every percentage, ratio and range, using pandas directly on the result files.
- **Long experiments:** none were re-run. The only script executed was the Theorem-1 check in `experiments/reproduce.py` (under 2 min).

**Verdicts**

| Verdict | Meaning |
|---|---|
| VERIFIED | Matches the evidence within rounding. |
| MISMATCH | The number differs from the evidence. The correct value is given. |
| UNSUPPORTED | No artefact in the repository supports the claim. |
| OVERCLAIM | The number is correct, but the wording overstates what the evidence shows (e.g. simulation vs real data, absolute words, causal language, selective reporting). |

---

## 1. Claim table

### 1.1 Abstract (l. 38–71)

| # | Line | Claim (abridged) | Evidence (file → value) | Verdict | Note |
|---|---|---|---|---|---|
| A1 | 53–56 | Twist: controller needs 11–17% fewer reads than a static allocation tuned in hindsight | `realdata_policies.csv`, `adaptive_pivotal` vs `static_oracle_weights`, `mean_spend@0.99`: RR 55,810/63,142 = −11.6%; reliability 54,059/61,947 = −12.7%; tiered 46,186/53,088 = −13.0%; single-best 44,753/53,867 = −16.9% | VERIFIED (numbers) / OVERCLAIM (attribution) | (a) The 17% end is the single-container placement. There, allocation is trivial and "static tuned in hindsight" is just fixed depth; the gain is from stopping only (table caption admits this, the abstract does not). (b) Adaptive stopping with round-robin allocation already gets 57.8k vs 63.1k on Twist RR, so allocation itself contributes about 3.5 percentage points. (c) The adaptive policy's own 0.99-quantile cost (`cap@0.99`) is 63.2k on Twist RR, equal to the oracle static budget. The saving is in *expected* spend only. |
| A2 | 55–56 | Paired bootstrap 95% CIs exclude zero | `diff_vs_static_oracle_lo/hi`: all upper bounds negative (−7.0k, −7.6k, −6.2k, −6.7k) | VERIFIED | Bootstrap of an extreme-quantile statistic: see the statistical review. |
| A3 | 56–57 | Combined with placement, 34% fewer than uniform fixed-depth | 44,753 (single-best, adaptive) / 67,938 (round-robin, fixed depth) = −34.1% | VERIFIED (number) / OVERCLAIM | "Placement" means putting all data in one container, chosen with a belief whose φ comes partly from evaluation reads. That layout loses the file if the container fails (Fig. 6, l. 634–635), and the abstract does not mention this trade-off. The figure also mixes a placement change with a stopping change. |
| A4 | 57–58 | Genscript: 11% and 63% | `realdata_genscript_policies_roundrobin.csv`: 577,653/650,709 = −11.2%. `…_part2.csv` single-best adaptive 306,210 vs RR fixed depth 831,222 = −63.2% | VERIFIED (numbers) / OVERCLAIM | The 63% compares different placements and different trial sets (100 vs 200 permutations). Same-placement single-best: 306,210 vs 409,341 = −25.2%. On Genscript RR the adaptive 0.99-quantile (`cap@0.99` = 974,400) is **50% higher** than the oracle static budget (650,710) and 17% higher than fixed depth. This tail cost is not reported anywhere. |
| A5 | 59–60 | Generic prior (no calibration reads) matches full calibration on pool 1 | `paper_e4_limited_calibration.csv`: generic 55,584 [55,134–56,026] vs full 55,284 [54,800–55,778] (+0.5%) | VERIFIED / minor OVERCLAIM | "Matches" is inferred from overlapping *unpaired* CIs. No equivalence test or paired difference is given, although all configurations share the same T matrices. |
| A6 | 60–61 | Pool 2: +4–5% vs in-pool calibration, still 10% below hindsight static | `paper_e10_cross_pool.csv`: generic 551.7/555.0 vs in-pool 528.9 → +4.3/+4.9%; vs oracle static 621.0 (`paper_e10_static_reference.csv`) → −11.2/−10.6% | VERIFIED | Only 60 permutations, so the 0.99-quantile is the sample maximum. In-pool E10 (528.9k) disagrees with Table II (577.7k, same pool, placement and belief type; CIs do not overlap). See F-1. |
| A7 | 61–63 | Flags strongly degraded containers after about 1–2 reads per stored strand | `paper_ed_change_detection.csv`: flag 1.00 for sim C0 (1.17 reads/strand), C2 (1.66) and real PCR_90c (1.08) | VERIFIED | Genscript PCR_30c needs 2.48 (`paper_ed_holdout_genscript.csv`). |
| A8 | 62–63 | "with no false alarms on near-nominal containers" | Containers within 5% of nominal: flag 0.00 | VERIFIED (narrowly) / OVERCLAIM | By the script's own ground truth (`truly_degraded = eff < 0.7`), PCR_60c (eff 0.75, *not* degraded) is flagged in **73%** of trials and Aging_2d (0.76) in 10% (82% for PCR_60c in the re-run). The abstract leaves out this false-alarm rate. |
| A9 | 63–64 | Stripe-importance rule is no better than hazard-greedy | Pivotal vs coverage: 55.8/55.7, 54.1/54.2, 46.2/46.4, 577.7/577.7k | VERIFIED (descriptive) | No paired test of pivotal vs coverage is in any result file (see R6). |
| A10 | 64–66 | Credible intervals over-confident; overdispersed likelihood largely fixes this on a held-out pool | `paper_ed_overdispersion_summary.csv`: Genscript coverage 0.0065 → 0.849 (min 0.22) | VERIFIED / minor OVERCLAIM | At ρ_U = 0.01, Twist coverage becomes 1.00 for a nominal 0.90, i.e. too wide. That is the same criterion the paper uses to reject ρ_U = 0.05 (l. 583–584). ρ_U was chosen post hoc. |
| A11 | 66–67 | Completion probability under-confident; early-abandonment rule harmful | `paper_ec_*` (real Twist); `q4_futility.csv` (simulated s7) | VERIFIED / OVERCLAIM (setting) | The harm of the futility rule was measured only in the **simulated** channel-change scenario s7, not on real reads. The abstract does not say so. |
| A12 | 67–68 | Surrogate 5–18% optimistic on one pool, 3.7× pessimistic on the other | `realdata_policies.csv` `pred_mean_cost`/`mean_spend@0.99` (pivotal): −5.2% to −17.8%. Genscript single-best 1,128,032/306,296 = 3.68× | VERIFIED / minor OVERCLAIM | On Genscript round-robin the surrogate predicts **infinite** cost (`pred_mean_cost = inf`), i.e. it calls the layout infeasible. "3.7×" understates the failure. |
| A13 | 68–70 | All results computational; no new synthesis or sequencing | — | VERIFIED | |
| A14 | 50–52 | "To evaluate policies on *real* reads rather than simulated ones" | `src/d3ra/realdata/empirical.py::empirical_T` | OVERCLAIM | Read *arrival order* is a permutation of real reads. Whether a cluster's consensus is correct is a random draw t_ij ~ φ̂_ij from a per-strand curve estimated from ≤ 24 reads × 2 draws (`strand_phi`). Containers are emulated by subsetting runs of the full pool. "Semi-empirical replay of real read streams" is more accurate. |

### 1.2 Introduction and contributions (l. 78–129)

| # | Line | Claim | Evidence | Verdict | Note |
|---|---|---|---|---|---|
| I1 | 80–81 | DNA demonstrated "at megabyte to hundreds-of-megabyte scale" [church, goldman, grass, erlich, organick] | literature_review.md | OVERCLAIM (minor) | Grass 2015 stored 83 kB; Church about 0.65 MB; Goldman about 0.74 MB. Use "kilobyte to hundreds-of-megabyte". |
| I2 | 83–84 | "every read carries substitution, insertion and deletion errors" | `realdata_channel.csv` `frac_error_free_reads` = 0.27–0.68 | OVERCLAIM | 27–68% of real Twist reads are error-free. Use "reads carry …". |
| I3 | 112–115 | Container selection "has no direct precedent we could find" | literature_review §14 | VERIFIED (hedged) | Acceptable hedged novelty statement. |
| I4 | 118 | "a strongest-possible static baseline tuned in hindsight" | `experiments/real_data.py::_oracle_static_weights` | OVERCLAIM | Nelder–Mead, 2 starts, maxiter 300, on a non-smooth order-statistic objective: a local optimum, not a proven optimum. It is also strongest only within the commit-whole-budget class. Static weights plus verified stopping (a deployable and stronger static allocation) is not evaluated. |
| I5 | 120–123 | Two independent real pools | critical_questions.md Q5 | OVERCLAIM (minor) | Both pools come from one study and one lab (PRJEB65931), and each is a single synthesis batch. l. 711 admits this. |
| I6 | 124–126 | "on two independent nanopore flowcells, a strand-specific error component persists across runs", 2.2× more accurate | `q3_summary.json`: MAE 0.1808/0.0833 = 2.17×; residual r = 0.569 | VERIFIED / minor OVERCLAIM | Both flowcells read the same synthesised molecules. A persistent strand-specific component is equally consistent with synthesis variants or reference errors (the q3b audit finds 14% design mismatch; real_data_validation.md §4 says so). The cause is not identified. The MAE comes from only 3 points (n = 3, 5, 10). |
| I7 | 128–129 | Everything computational | — | VERIFIED | |

### 1.3 Model and method (l. 176–381)

| # | Line | Claim | Evidence | Verdict | Note |
|---|---|---|---|---|---|
| M1 | 235 | 168–172 nt strands | not checked against a result file | UNSUPPORTED (not re-derived) | Plausible from the format, but no artefact states it. |
| M2 | 238 | Default K = 32, N = 48 | `real_data.py policies(k=32, parity=16)` | VERIFIED | |
| M3 | 305–306 | Grid {0} ∪ 25 log-spaced points in [0.01, 5] | `adaptive_policy.py` l. 41 | VERIFIED | |
| M4 | 315–317 | A floor of 0.2 caused abandonment in an earlier version | `results/archive_v2_narrow_theta/`, commit d991e8c | VERIFIED | |
| M5 | 329 | H = 6 batches per container | `AdaptiveConfig.horizon = 6.0` | VERIFIED | |
| M6 | 377–378 | q* = 0.232 for D = 4, t = 1 | `theory_qstar.csv` 0.232408; recomputed root 0.232408 | VERIFIED | |
| M7 | 380 | "sign prediction correct in all 8,148 covered cases of an independent numerical reproduction" | `experiments/reproduce.py::theory()` re-run by auditor → 8148/8148 PASS. `theory_summary.json`: T1 agreement 0.9955 over 1,344 cases (6 uneven-split mismatches) | VERIFIED / minor OVERCLAIM | The 8,148 count is not saved in any result file (it is computed at run time), and the script is by the same authors ("independent" = separate script). State the D \| N restriction. |
| M8 | 386 | 202 unit and integration tests | `evidence/test_reports/pytest.txt`: 202 passed | VERIFIED | Updated during the audit (was 201). Only 135 `def test_` functions: the rest are parametrisations. |
| M9 | 346 | "A wrong file is therefore output only on a SHA-256 collision" | `decoding/reconstruction.py` | VERIFIED (logic) | Requires an intact out-of-band manifest (critical_questions.md "remaining assumption"). Not stated in the paper. |

### 1.4 Evaluation methodology (l. 399–453)

| # | Line | Claim | Evidence | Verdict | Note |
|---|---|---|---|---|---|
| E1 | 402–404 | A container holding a subset "produces **exactly** the sample's read stream restricted to that subset" | modelling assumption (`empirical.py` docstring) | OVERCLAIM | This holds only if physically subsetting a pool leaves relative amplification and sampling unchanged, and junk reads scale with the subset share. It is an assumption, not an identity. |
| E2 | 405 | 99.5% of reads mapped | `realdata_channel.csv` `yield_mapped`: Twist 99.48–99.75% (mean 99.62%); `realdata_genscript_channel.csv`: 98.15–99.37% (mean 98.92%) | MISMATCH | Say "99.5–99.7% (Twist), 98.2–99.4% (Genscript)". |
| E3 | 406–409 | Strand-specific φ̂_ij from repeated consensus of n real reads | `empirical.py::strand_phi` (draws = 2, ≤ 24 reads/strand, early stop after two successes) | VERIFIED / note | Each grid value comes from only 2 subsamples. Disclose this. |
| E4 | 414–415 | Beliefs estimated from a disjoint calibration split (first 300k; evaluation next 900k), "except where stated" | `real_data.py::belief_from_calibration`: `phis.append(EM.grid_to_full(rd.phi_strand[::2].mean(0)))`, where `phi_strand` is estimated from **evaluation-split** reads (`load_run` stores eval reads in `reads_by_ref`) | MISMATCH (undisclosed leakage) | The "full calibration" belief's consensus curve uses evaluation reads of half of the evaluated strands. It is never "stated" in the manuscript. real_data_validation.md l. 27 calls it "a small leak". supplementary/limitations.md l. 19 says the opposite. The leaked belief drives Table II (adaptive rows), E4 "full calibration", EC, ED (fresh-pool φ), E10 in-pool and the single-best container choice. |
| E5 | 441 | 12,000 strands in 250 stripes of (32,48) | `real_data.py` (n_stripes = 250) | VERIFIED | For Genscript, a random 12,000 of 12,472 strands are used. Say so. |
| E6 | 444–448 | Oracle-tuned static: "no deployable static plan can do better" | Nelder–Mead local search | OVERCLAIM | See I4. |
| E7 | 450–451 | Paired bootstrap, 2,000 resamples [efron]; effect sizes as Cliff's δ [cliff] | `_paired_spend_diff(n_boot=2000)`; no Cliff's δ in any paper result file (`grep -i cliff results/paper_* results/realdata_*` → none) | MISMATCH / UNSUPPORTED | Cliff's δ is announced but never reported. Remove the sentence and citation, or report it. |
| E8 | Table `config` | Prior sd 0.5 (default), 1.0 in sensitivity runs | `paper_experiments.py::ed` uses `prior_sd=1.0`; `q4` futility uses 1.0 | MISMATCH (minor) | The change-detection results (Table `detect`, overdispersion, Genscript hold-out) use s = 1.0. Say so in the table and text. |
| E9 | Table `config` | Simulation scenarios s1–s7, 2,000 trials each | `ablations.csv`: B and E variants use 10,000 trials | MISMATCH (minor) | |
| E10 | Table `config` | Trials: 60 for failure injection and E10 | `failure_injection(trials=60)`, `e10(trials=60)` | VERIFIED | |

### 1.5 Results: real channels (l. 460–475, Fig. `consensus`)

| # | Line | Claim | Evidence | Verdict | Note |
|---|---|---|---|---|---|
| R1 | 461–462 | Substitution-dominated, 0.35–1.15%/nt | `realdata_channel.csv` `rate_sub` 0.352–1.150% | VERIFIED | Total error 0.42–1.27%/nt. |
| R2 | 462–463 | 3 reads give exact consensus for 87–99% of strands | `phi_exact_n3` 0.873–0.987 | VERIFIED | |
| R3 | 463–464 | Over-amplified samples need "about twice as many reads" | `phi_exact_n*`: to reach 0.95, Aging_0d needs n = 3 and PCR_90c needs n = 4; at n = 1, 0.52 vs 0.27 | UNSUPPORTED (vague) | No metric gives a factor of 2. Define the criterion or drop the claim. |
| R4 | 464–465 | 69% of excess coverage variance shared (52% Genscript) | `realdata_shared_variance.json` 0.693; `realdata_genscript_shared_variance.json` 0.518 | VERIFIED / note | The Genscript 52% uses only 4,940 of 12,472 strands (zero-coverage strands in the near-dead run are excluded). Without PCR_30c ("mild runs") it is 91%. "An upper bound for aliquots of one synthesis batch" is muddled: it is an upper bound for what *independently synthesised* containers would share. |
| R5 | 466–467 | Genscript deletion-dominated, CV ≈ 1.0 (Twist 0.35), one container with 54% dropout | `rate_del` > `rate_sub`; CV 0.98–1.06 (PCR_30c 2.0); Twist 0.35–0.42; PCR_30c `dropout_at_full_depth` 0.538 | VERIFIED | |
| R6 | 470–471 | i.i.d. MAE 0.181 → hierarchical 0.083 | `q3_summary.json` 0.1808 / 0.0833 | VERIFIED | The i.i.d. model under-predicts at n = 3 (0.285 vs 0.298). "Over-predicts" holds for n = 5, 10 only. |
| R7 | 472–473 | r = 0.57 between flowcells (null 0.002) | `rate_corr_residual_strand_specific` 0.569; null 0.0016 | VERIFIED | The null depends on Python `hash()` ordering (`q3_correlated_errors.py` l. 68), so it is not bit-reproducible. |
| R8 | 473–475 | After excluding 14% mismatching designs, 10-read exact consensus 0.937 vs 0.997 i.i.d. | `q3b_reference_audit.json`: 1 − 0.8596 = 14.0%; 0.9373 vs 0.9971 | VERIFIED | |

### 1.6 Results: cost (l. 486–509, Table `baselines`, Fig. `cost`)

| # | Line | Claim | Evidence | Verdict | Note |
|---|---|---|---|---|---|
| C1 | 488–489 | 6,902–9,114 fewer reads (11–17%) | `diff_vs_static_oracle` −6,901.9 … −9,114 | VERIFIED | See A1. |
| C2 | 489 | Paired 95% CIs exclude zero for every placement | all `_hi` < 0 | VERIFIED | |
| C3 | 491–492 | 44,753 vs 67,938 (−34%) | as A3 | VERIFIED / OVERCLAIM | as A3 |
| C4 | 492–494 | Genscript 577,653 vs 650,709 (−11%); 306,210 vs 831,222 (−63%) | as A4 | VERIFIED / OVERCLAIM | as A4 |
| C5 | 494–497 | Pivotal and coverage-greedy "statistically indistinguishable in all cases" | No test result for pivotal − coverage in any file | UNSUPPORTED (statistic) | Point estimates differ by ≤ 0.4%. Add a paired bootstrap of the difference. |
| C6 | 496–497 | "the gain comes from hazard-based allocation and stopping" | `adaptive_uniform` rows | VERIFIED (qualitative) | Twist RR: stopping −5.4k, allocation −2.0k. Genscript RR: stopping −9.0k (CI −58.5k to +1.1k, includes 0), allocation −64k. Report this decomposition explicitly. |
| C7 | 497–499 | Coverage-equalising static allocator 2.3× worse than uniform on Genscript | 1,938,018/831,222 = 2.33× | VERIFIED | |
| C8 | 499 | "the adaptive controller learns to avoid it" | E10 read share 3.4% vs 16.7% uniform | VERIFIED | |
| C9 | Table `baselines` | All 30 cells (fixed / oracle / stop-RR / cov / pivotal / Δ and CIs) | `realdata_policies.csv`, `realdata_genscript_policies_*.csv` | VERIFIED | All values match to 0.1k. |
| C10 | Fig. `cost` caption | "Reads needed for exact recovery with probability ≥ 0.99" | `make_figures.py::fig_cost_real` plots `mean_spend@0.99` | OVERCLAIM (wording) | For adaptive bars this is the *expected* spend, not the budget needed for 0.99. Say "expected reads (spend capped at the 0.99-quantile)". |

### 1.7 Results: limited calibration (l. 511–535, Fig. `limited`)

| # | Line | Claim | Evidence | Verdict | Note |
|---|---|---|---|---|---|
| L1 | 513–514 | Generic prior: yield 0.9, NB shape 10, dropout 0.5%, generic φ | `paper_experiments.py::generic_prior` | VERIFIED | Also σ_g = 0.3, close to the Twist estimate 0.31 (`realdata_belief.json`). This is not stated, and it makes the prior less "generic" than claimed. |
| L2 | 517–518 | All 14 configurations reach the target in every trial; 55,254–55,726 reads | `paper_e4_limited_calibration.csv`: `p_success_uncapped` = 1.0; min 55,254, max 55,726 | VERIFIED | |
| L3 | 518–519 | CIs overlap; largest difference 0.9% | 55,726/55,254 − 1 = 0.85% | VERIFIED | Unpaired overlap on the same T matrices: report paired differences. |
| L4 | 519–521 | 55,584 [55,134–56,026] vs 55,284 [54,800–55,778]; oracle static 63,142 | E4 csv; `realdata_policies.csv` | VERIFIED | The oracle reference comes from a different (200-permutation) trial set. The Fig. `limited` caption says so; the text does not. |
| L5 | 524 | Twist relative efficiencies 0.34–1.03 | `paper_ed_change_detection.csv` 0.337–1.032 | VERIFIED | |

### 1.8 Results: confidence calibration (l. 537–558, Fig. `confidence`)

| # | Line | Claim | Evidence | Verdict | Note |
|---|---|---|---|---|---|
| P1 | 538 | Budget = median adaptive cost 55,600 | `paper_ec_meta.json` cap 55,600 | VERIFIED | |
| P2 | 542 | Mean P̂_complete is 0.04–0.22 | `paper_ec_confidence_bins.csv` `mean_p`: full 0.043–0.216; generic **0.0005**–0.199 | MISMATCH | Correct range: 0.001–0.22 (0.04–0.22 for full calibration only). |
| P3 | 542–543 | 52% / 60% complete | `empirical` 0.515 / 0.600 | VERIFIED | |
| P4 | 543–544 | ECE 0.30–0.60; Brier 0.31–0.60 vs constant 0.24–0.25 | ECE 0.299–0.599; Brier 0.312–0.599; climatology 0.240/0.250 (recomputed) | VERIFIED | |
| P5 | 545–547 | AUC 0.49–0.56 at 30–50% (CIs include 0.5), 0.70–0.71 at 85% (CI 0.63–0.77) | `paper_ec_discrimination.csv`: 0.491–0.556; 0.706/0.702; CI 0.629–0.776 | VERIFIED | Upper bound rounds to 0.78 for the generic prior. |
| P6 | 547–548 | "this explains why the futility rule built on it abandoned recoverable retrievals" | futility measured in simulated s7 (`q4_stopping_safety.py` part D); EC measured on real Twist | OVERCLAIM (causal) | Different settings, and no analysis links the two. Use "is consistent with". |

### 1.9 Results: change detection (l. 560–592, Table `detect`)

| # | Line | Claim | Evidence | Verdict | Note |
|---|---|---|---|---|---|
| D1 | 561–562 | R* is about 1–2 reads per stored strand for most containers | Twist 1.08 (all 6); sim C0/C1/C2 1.17–1.66; sim C3–C7 2.12–8.97 | VERIFIED (borderline) | 9 of 14 containers. Five simulated containers need 2.1–9.0. |
| D2 | 565–566 | Strongly degraded flagged in 100% of trials | flag 1.00 (C0, C2, PCR_90c at 1× and 2× R*) | VERIFIED | |
| D3 | 566 | Containers within 5% of nominal never flagged | flag 0.00 | VERIFIED | See A8 for the false alarm on PCR_60c. |
| D4 | 566–568 | Posterior means 0.28–0.29 vs 0.34; 1.06 vs 1.00–1.03; sim 0.05–0.06 vs 0.14 | csv | VERIFIED | Selective: Aging_7d 0.82 vs 0.90 and PCR_60c 0.63 vs 0.75 are not mentioned. |
| D5 | 568–570 | Moderately degraded "flagged occasionally": 2 d 10%, 60 c 73% at R*, 0% at 2R* | csv | VERIFIED (numbers) / OVERCLAIM | 73% is not "occasionally". The container is not degraded by the script's own definition (eff 0.75 > 0.7), so this is a false alarm. |
| D6 | 570–571 | 90% CIs cover truth in only 0–83% | csv `ci90_coverage` 0.00–0.83 | VERIFIED | |
| D7 | 579–580 | Coverage 0.08 → 1.00 (Twist), 0.51 → 0.79 (sim) | `paper_ed_overdispersion_summary.csv` 0.081 → 1.000; 0.507 → 0.790 | VERIFIED | |
| D8 | 581–582 | 60-cycle 0.82 → 0; binomial re-run 82% vs 73% | 0.815 → 0.000; 0.73 | VERIFIED | |
| D9 | 582–583 | MAE 0.060 → 0.023 real; 0.038 → 0.027 sim | 0.0600 → 0.0229; 0.0381 → 0.0269 | VERIFIED | |
| D10 | 583–584 | ρ_U = 0.05 gives 100% coverage everywhere | 1.00 (selection set) | VERIFIED | ρ_U = 0.01 also gives 1.00 on Twist, i.e. also too wide there (A10). |
| D11 | 585–586 | Genscript: 0.01 → 0.85 (range 0.22–1.00), no false flags | `paper_ed_overdispersion_summary.csv` 0.0065 → 0.849 (0.22–1.00); `flag_healthy_max` 0 | VERIFIED | |
| D12 | 586–588 | PCR_30c flagged in all trials under both likelihoods | flag 1.0 | VERIFIED | |
| D13 | 589–590 | Genscript posterior MAE 0.34 | 0.340 / 0.335 | VERIFIED | |
| D14 | Table `detect` | All 10 rows | `paper_ed_change_detection.csv` (1× R* rows) | VERIFIED | Shows only C0, C1, C2, C5 of the 8 simulated containers. Mention that C3, C4, C6, C7 are omitted. |

### 1.10 Results: generalisation (l. 594–630, Table `transfer`)

| # | Line | Claim | Evidence | Verdict | Note |
|---|---|---|---|---|---|
| G1 | 600–602 | 551.7–557.1k; 4–5% above in-pool 528.9k; unpaired CIs overlap | `paper_e10_cross_pool.csv`: +4.3% to +5.3% | VERIFIED | Paired comparisons are possible (same 60 matrices) but not done. |
| G2 | 602–603 | 10–11% below oracle static 621.0k; 30–31% below fixed 798.2k | −10.3% to −11.2%; −30.2% to −30.9% | VERIFIED | |
| G3 | 604–605 | 5.3–5.9% of reads to the near-dead container vs 3.4% | `read_share_broken_container` 0.053–0.059; 0.034 | VERIFIED | |
| G4 | 605–606 | Prior width changes cost by < 1% | max 0.6% | VERIFIED | |
| G5 | Table `transfer` | All 8 rows incl. CIs and shares | `paper_e10_*` | VERIFIED | Fixed-depth share 16.7% = 1/6 (uniform). |
| G6 | 600–602 vs 492–493 | (implicit) consistency of in-pool Genscript cost | E10 in-pool 528,873 [513,860–544,384] (60 perms) vs Table II pivotal 577,653 [560,581–595,229] (200 perms). Static refs also differ: 621.0 vs 650.7k, 798.2 vs 831.2k | MISMATCH (internal inconsistency, not discussed) | The CIs do not overlap: a 9% gap, larger than the 4–5% transfer effect being claimed. Likely cause: the Genscript φ̂ cache was regenerated (pickles dated 2026-11-02 15:10–15:15) with process-randomised `hash(run)` seeds after the Table II results were produced (committed earlier, 0c322f9 had deleted the caches). Report and explain, or re-run both on one cache. |

### 1.11 Results: failures and stopping safety (l. 632–655, Fig. `failure`)

| # | Line | Claim | Evidence | Verdict | Note |
|---|---|---|---|---|---|
| F1 | 634–635 | Single-container layout 44,473 reads; lost when that container fails | `realdata_failure_injection.csv` 44,473; PCR_15c lost → p = 0 | VERIFIED | |
| F2 | 636–638 | Tiered 46,530; 127–137k after losing a main container; 45–48% within 3× budget; third container 102k, always completes | 126,586 / 136,681; 0.483 / 0.450; Aging_2d 101,797, p = 1 | VERIFIED | |
| F3 | 638–639 | Round-robin 55,610; +5 to +32%; always completes | 58,347/55,610 = +4.9%; 73,363/55,610 = +31.9% | VERIFIED | The failure-injection controller is **coverage-greedy**, not pivotal (`failure_injection` uses `AdaptiveConfig("coverage")`). Not stated. |
| F4 | 640 | "consistent with coupon-collector theory [barlev2025cover, levy2026expected]" | no comparison computed | UNSUPPORTED | No theoretical prediction is evaluated against the measured +5…+32%. |
| F5 | 641–642 | Simulation, 12 retrievals per configuration, heavy errors | `q4_end_to_end_weak_code.csv`, 12 seeds × 2 | VERIFIED | The file is 6,000 random bytes over 4 containers. That is tiny and should be stated. |
| F6 | 643 | 0 false accepts in **10,567** decode attempts | `q4_summary.json` `decode_attempts_total` = **10,557** (sum of the CSV column) | MISMATCH | Correct value: 10,557. Same error in `tables/evidence.tex` l. 20 and `evidence_traceability.csv` C14. |
| F7 | 643–644 | 95% upper bound 2.8 × 10⁻⁴ | 3/10,557 = 2.84 × 10⁻⁴ | VERIFIED / note | This is a bound per decode attempt. Per *wrong* consensus word the relevant bound is 3/2,682 = 1.1 × 10⁻³ (`q4_inner_false_accepts.csv`). State which. |
| F8 | 644–645 | Weak code: "54–60 wrong fragments per retrieval passed the inner check" | per-run `false_accepts` = 45, 45, 49, 50, 51, 51, 54, 57, 59, 60, 63, 65 (mean 54.1) | MISMATCH | Correct: 45–65 per retrieval (mean 54). |
| F9 | 645–646 | All were caught; retrievals ended *partial* rather than wrong | 11 partial, **1 failed** (seed 3); `wrong_file_output` all False | MISMATCH (minor) | Write "11 ended partial and 1 failed; none output a wrong file". |

### 1.12 Results: ablation (l. 657–673, Fig. `ablation`)

| # | Line | Claim | Evidence | Verdict | Note |
|---|---|---|---|---|---|
| B1 | 658–659 | s4: 8 containers in 4 domains; 2,000 trials | `scenarios.py` s4 (domains [0,1],[2,3],[4,5],[6,7]) | VERIFIED | Variants B and E use 10,000 trials. |
| B2 | 659–660 | Removing adaptive allocation costs **+79%** | `ablations.csv`: `mean_spend@0.99` 9,322.9/5,201.8 = +79.2%, but the figure (`mean_cost_success`, per its caption) shows 8,559.2/5,045.9 = **+69.6%** | MISMATCH (text vs figure metric) | The text and Fig. `ablation` use different metrics. Use +70% (the caption's metric) or change the caption. |
| B3 | 660 | Static reads cost +68% | `mean_cost_success` 8,426.4/5,045.9 = **+67.0%** (spend@0.99 infinite: p = 0.987) | MISMATCH (minor) | +67%. It is also a conditional-on-success metric with unequal success rates. |
| B4 | 660–661 | Uniform or independent placement +37–38% and lower reliability | +37.8% (p 0.984), +37.4% (p 0.9825) | VERIFIED | |
| B5 | 661–662 | Ignoring correlated failures drops success to 0.953 | p 0.9525 | VERIFIED (rounding; 0.952 also defensible) | |
| B6 | 662 | Perfect channel knowledge changes cost by < 1% | 5,032.9 vs 5,045.9 (−0.26%) | VERIFIED / OVERCLAIM | In s4 the belief **equals** the truth (`Scenario(..., ch, ch)`). The "oracle" variant only adds the true failure state, so the result is uninformative about channel-knowledge value. |
| B7 | 663 | Coverage-greedy 1.5% cheaper | spend −1.49%; figure metric −1.57% | VERIFIED | Same metric inconsistency as B2. |
| B8 | 664–665 | Seed sweep: full design most reliable in **3/5** seeds; does not reliably reach 0.99 | `q5_seed_sensitivity.csv` s4 `E_joint_adaptive_pexact` 0.978, 0.989, 0.977, 0.985, 0.980 vs `C_uniform_adaptive_pexact` 0.973, 0.983, **0.987, 0.989, 0.982** | MISMATCH | E is most reliable in **2/5** seeds (seeds 0, 1). C wins seeds 2–4. E reaches 0.99 in **0/5** seeds (max 0.989), which is stronger than "does not reliably reach". The ablation figure is from the one run where E hit exactly 0.990. |
| B9 | Fig. `ablation` caption | "green: P ≥ 0.99" | full design P = 0.9900 | VERIFIED | Marked green only because the value sits exactly at the threshold. |

### 1.13 Discussion, limitations, reproducibility, conclusion (l. 676–747)

| # | Line | Claim | Evidence | Verdict | Note |
|---|---|---|---|---|---|
| X1 | 678–683 | Two mechanisms (stopping; directing reads away from slow containers) explain most of the gain; stripe weighting adds nothing measurable | see C6, A9 | VERIFIED (qualitative) | |
| X2 | 679 | "stopping as soon as the data are **provably** decodable" | in the real-data replay, decodability is an oracle check (`T <= R`, consensus == design), not a CRC/SHA verification | OVERCLAIM (minor) | |
| X3 | 691–692 | Surrogate −5 to −18% (Twist), 3.7× (Genscript) | as A12 | VERIFIED | |
| X4 | 693–694 | Controller "is robust to this because it uses the model only to rank containers" | no experiment isolates this | OVERCLAIM (causal) | The controller uses the same Q-tables for the posterior and hazards, not only for ranking. E10 shows robustness empirically but not the stated reason. |
| X5 | 700–702 | Same-batch aliquots on two flowcells; 69% shared variance | `q1_summary.json` (different-flowcell pairs); `realdata_shared_variance.json` | VERIFIED | |
| X6 | 707–709 | Futility rule abandoned 3.7–4.5% of retrievals that would have succeeded; disabled | `q4_futility.csv` 0.0367–0.045, `wrongful_give_up` = `gave_up` | VERIFIED / OVERCLAIM (setting) | Simulated s7 only. Say "in simulation". |
| X7 | 711 | External validity: two Illumina pools from one study, two nanopore datasets | — | VERIFIED | Only one nanopore dataset (DNAformer) contributes results to the paper. The Srinivasavaradhan data appear only in the q3b audit (0.91 match), which is not reported. |
| X8 | 719–721 | Public data fetched by `experiments/real_data.py download` (ENA PRJEB65931; Zenodo 12666365, 17473983, 14296588) | `real_data.py::download` fetches ENA runs, 12666365 and 14296588 **only** | MISMATCH | Zenodo 17473983 (DNAformer) is not downloaded by that command. |
| X9 | 721–722 | Independent clean-environment reproduction passed all checks, byte-identical theory outputs | critical_questions.md Q5 cites `scratchpad/repro_report.md`, which is not in the repository; it was done on commit f8fb9a9 (before E4/EC/ED/E10 existed) and after 4 fixes | UNSUPPORTED (artefact) / OVERCLAIM | Commit the report. State that it covers the theory and simulation checks plus the earlier real-data headline, not the paper's RQ1–RQ3 experiments. |
| X10 | 723–726 | Not bit-for-bit reproducible because of process-randomised hash seeding; 73% vs 82% | `real_data.py` l. 88; `q3_correlated_errors.py` l. 68 (second hash use, not mentioned) | VERIFIED / incomplete | This also probably explains G6 (a 9% shift in Genscript costs), which is far more than "small run-to-run differences". |
| X11 | 734–737 | Conclusion: 11–17% vs hindsight static; 34–63% with placement | as A1, A3, A4 | VERIFIED / OVERCLAIM | as A1, A3, A4 |
| X12 | 737–738 | "in simulation, verification-gated stopping **never** produced a wrong file" | 0/24 end-to-end runs (12 default and 12 weak code), 6 kB file | OVERCLAIM | With 0/24 the 95% upper bound on the wrong-file rate is about 12%. Write "produced no wrong file in 24 simulated retrievals". |
| X13 | 738–739 | Works from a generic prior (cost unchanged on one pool, +4–5% on the other) | A5, A6 | VERIFIED | |
| X14 | 739–740 | Detects strongly degraded containers within 1–2 reads per stored strand | A7 | VERIFIED | |

### 1.14 Tables `literature_comparison`, `limitations`, `evidence`, `notation`

| # | Location | Claim | Evidence | Verdict | Note |
|---|---|---|---|---|---|
| T1 | `literature_comparison.tex` l. 29 | "This work": Correlated errors = **Y** | l. 368–370: the beta-binomial model "is not yet coupled into the controller" | OVERCLAIM | Should be P. Correlated errors enter only through the calibrated φ_j. |
| T2 | `limitations.tex` l. 15–17 | ECE 0.30–0.60; futility 3.7–4.5%; surrogate −5…−18% / 3.7× | as above | VERIFIED | The futility row should say "simulated". |
| T3 | `evidence.tex` l. 20 | "0/10,567 false accepts" | 10,557 | MISMATCH | |
| T4 | `evidence.tex` l. 12–23 | Status "verified, real reanalysis" etc. | traceability CSV | VERIFIED | Row 15 ("Zero-calibration prior matches full calibration") inherits the E4 leakage caveat. |
| T5 | `notation.tex` l. 28 | ε, B_ε implemented in `realdata.empirical` | `budget_for_target` is in `retrieval/fixed_policy.py`; spend in `evaluation/metrics.py` | MISMATCH (minor) | |
| T6 | `notation.tex` | h_ij is both container-specific abundance (Eq. 1) and recovery hazard (Eq. 6) | — | clarity | Symbol clash. ρ is likewise overloaded: ICC, ρ_i, ρ_par, ρ_U. |

---

## 2. Summary counts

| Verdict | Count (primary verdict per row) |
|---|---|
| VERIFIED (some with notes) | about 95 |
| MISMATCH | 16: E2, E4, E7, E8, E9, P2, G6, F6, F8, F9, B2, B3, B8, X8, T3, T5 |
| UNSUPPORTED | 6: M1, R3, C5, F4, X9, and E7 (Cliff's δ never reported) |
| OVERCLAIM (alone or alongside a correct number) | about 25: A1, A3, A4, A5, A8, A10, A11, A12, A14, I1, I2, I4, I5, I6, E1, E6, C10, P6, D5, B6, X2, X4, X6, X11, X12, T1 |

---

## 3. Required corrections

### 3.1 Hard number fixes (MISMATCH)

1. **l. 643, `tables/evidence.tex` l. 20, `evidence_traceability.csv` C14:** "10,567" → **10,557** decode attempts. The bound stays at 2.8 × 10⁻⁴. Also give the per-wrong-word bound, 1.1 × 10⁻³ (0/2,682).
2. **l. 644–645:** "54–60 wrong fragments per retrieval" → **45–65 (mean 54)**.
3. **l. 645–646:** "the retrievals ended partial rather than wrong" → **"11 ended partial and 1 failed; none output a wrong file"**.
4. **l. 659–660 (ablation):** "+79%" contradicts the figure's metric (mean reads of successful retrievals = **+69.6%**). "+68%" → **+67%**. Either quote the figure's metric throughout, or state that the +79% / −1.5% figures are expected spend at 0.99.
5. **l. 664–665:** "most reliable in 3/5 seeds" → **2/5 seeds** (uniform placement with adaptive reading is most reliable in seeds 2–4). Add that the full design reaches 0.99 in **none** of the five seeds (0.977–0.989). Also fix README.md l. 119, reports/technical_report.md l. 41 and l. 79, reports/paper_draft.md l. 15 and l. 148, and reports/critical_questions.md l. 149.
6. **l. 542:** "its mean is 0.04–0.22" → **0.001–0.22** (generic prior 0.0005–0.20; full calibration 0.04–0.22).
7. **l. 405:** "99.5% mapped" → **99.5–99.7% (Twist), 98.2–99.4% (Genscript)**.
8. **l. 414–415 (leakage):** disclose that the belief's consensus curve φ_j (in "full calibration", the fresh-pool belief of ED, the in-pool belief of E10, and the single-best container choice) is averaged over even-indexed strands of the **evaluation** split (`belief_from_calibration`). Better: re-estimate φ_j from calibration reads and re-run. Also correct `supplementary/limitations.md` l. 19, which states the opposite.
9. **Genscript inconsistency (G6):** Table II (577.7k / 650.7k / 831.2k, 200 perms) and Table `transfer` (528.9k / 621.0k / 798.2k, 60 perms) disagree beyond their CIs. Re-run both from one φ̂ cache with a stable seed, or explain the gap. The likely cause is regeneration of the Genscript cache with `hash(run)` seeding.
10. **l. 450–451:** Cliff's δ is announced but never reported. Report it, or delete the clause and `\cite{cliff1993dominance}`.
11. **Table `config`:** prior sd 1.0 is used in all change-detection experiments, not only in "sensitivity runs". Variants B and E of the ablation use 10,000 trials.
12. **l. 719–721:** `real_data.py download` does not fetch Zenodo 17473983. Fix the script or the sentence, and add the download step for `q3_correlated_errors.py`.
13. **`tables/notation.tex` l. 28:** ε, B_ε live in `retrieval/fixed_policy.budget_for_target` and `evaluation/metrics`, not `realdata.empirical`.

### 3.2 Supply missing evidence (UNSUPPORTED)

1. **l. 494–496:** add a paired bootstrap CI for pivotal − coverage-greedy on every placement before writing "statistically indistinguishable".
2. **l. 640:** either compute the coupon-collector prediction (e.g. Bar-Lev E[reads] = n(H_n − H_{n−k}) on surviving fragments, or Levy's 1/(surplus+1) scaling) and compare it with the measured +5…+32%, or delete "consistent with coupon-collector theory".
3. **l. 463–464:** define "about twice as many reads" (e.g. reads to reach φ = 0.95) or remove it.
4. **l. 721–722:** commit the clean-environment reproduction report (currently `scratchpad/repro_report.md`, not in the repository). Say that it predates the RQ1–RQ3 experiments (commit f8fb9a9).
5. **l. 380:** save the 8,148-case output of `reproduce.py` as a result file. State the D | N restriction, and that 6 of 1,344 cases with uneven splits disagree (`theory_summary.json`).
6. **l. 235:** add the 168–172 nt strand length to an artefact, or cite the format code.

### 3.3 Wording fixes (OVERCLAIM)

1. **Abstract l. 53–58 and Conclusion l. 734–737.**
   - Say that most of the 11–17% comes from verified stopping, and that the 17% case is a single container where allocation plays no role.
   - Report the stopping-only vs allocation decomposition: Twist RR −5.4k / −2.0k; Genscript RR −9.0k / −64k.
   - Report that the adaptive 0.99-quantile cost is not lower than the oracle static budget (Twist RR 63.2k vs 63.1k; Genscript RR 974k vs 651k).
   - Qualify the 34% / 63% "combined with placement" figures. They put all data in one container, which loses the file if that container fails, and they compare across placements.
2. **Abstract l. 62–63:** "no false alarms on near-nominal containers" → add "but a container at 0.75 relative efficiency (above the 0.7 threshold) is flagged in 73% of trials at R*". Also fix l. 569 "occasionally".
3. **Abstract l. 66–67, l. 707–708, `limitations.tex` l. 16:** mark the futility result as **simulated** (scenario s7).
4. **l. 547–548:** "this explains why" → "this is consistent with".
5. **Abstract l. 50–52, l. 402–404, l. 458:** describe the replay as semi-empirical. Real read order and abundance are used; consensus success is sampled from per-strand empirical curves; containers are emulated by subsetting runs of one pool. "Exactly" (l. 403) → "under the assumption that subsetting does not change relative abundance".
6. **l. 118, l. 447–448:** "strongest-possible", "no deployable static plan can do better" → "a hindsight-tuned (Nelder–Mead) static plan, an optimistic reference for fixed-budget static plans". Add the missing baseline: oracle-tuned static weights **with** verified stopping.
7. **l. 737–738:** "never produced a wrong file" → "produced no wrong file in 24 simulated retrievals of a 6 kB file".
8. **l. 662:** "Perfect channel knowledge changes the cost by < 1%" → explain that in s4 the belief equals the truth, so this variant only tests knowledge of the failure state.
9. **`literature_comparison.tex` "This work / Correlated errors":** Y → P (footnote: enters only through calibrated φ_j; not coupled to the controller).
10. **Abstract l. 67–68:** add that on Genscript round-robin the surrogate predicted the layout infeasible (infinite cost).
11. **l. 80–81:** "megabyte to hundreds-of-megabyte" → "kilobyte to hundreds-of-megabyte".
12. **l. 83:** "every read carries" → "reads carry".
13. **l. 124–126, l. 472–473:** state that the cross-flowcell strand-specific component may reflect synthesis variants or reference errors in the shared molecules, not only correlated sequencing errors.
14. **l. 693–694:** drop or soften the causal "robust … because it uses the model only to rank".
15. **l. 679:** "provably decodable" → "verified decodable" (in the real-data replay, decodability is an oracle check against the design).
16. **l. 583–584 / abstract l. 65:** acknowledge that ρ_U = 0.01 also gives 100% coverage (too wide) on Twist, so "largely fixes" holds on Genscript only.

---

## Author response and corrections applied (after the audit)

**MISMATCH**
- **1–3 (decode attempts, weak-code counts, partial/failed).** Corrected to 10,557 attempts, 45–65 wrong fragments per retrieval (mean 54), and 11 partial / 1 failed. Fixed in `main.tex`, `tables/evidence.tex`, the traceability CSV, README and reports.
- **4 (ablation).** Now reported as +70% in mean cost of successful retrievals (+79% in expected spend) and +67% for static reads. Coverage-greedy is −1.6%.
- **5 (seed sweep).** Corrected to "most reliable in 2 of 5 seeds; reaches 0.99 in none". Fixed in the manuscript, README, `paper_draft.md`, `technical_report.md` and `critical_questions.md`.
- **6 (forecast mean).** Corrected to 0.001–0.22.
- **7 (mapping rate).** Corrected to 99.5–99.7% (Twist) and 98.2–99.4% (Genscript).
- **8 (evaluation leak).** Confirmed in `belief_from_calibration`.
  - Now disclosed in Section VI, the limitations, the literature-table footnote, S4/S5 and the traceability CSV (C32).
  - Its effect is bounded by the leak-free E4 generic-φ rows (≤0.9%) and the E10 generic/transfer priors (4–5%).
  - The code was not changed, so all reported numbers are reproducible as they stand.
- **9 (Genscript inconsistency).** Root cause confirmed: the Genscript φ̂ cache was regenerated (15:10–15:15) after the first policy run.
  - All real-data policy tables and the failure injection were re-run on the current caches (`*_v2.csv`).
  - Twist v2 is identical to v1. Genscript changed: round-robin oracle-static 617.9k, adaptive 531.9k; single best 279.4k → 219.0k.
  - E10 now agrees with Table II. The superseded files are kept and labelled.
- **10 (Cliff's δ).** The mention was removed.
- **11 (config table).** Corrected: prior sd 1.0 in change detection; 2,000 adaptive / 20,000 static trials, and 10,000 for the static ablation variants.
- **12 (DNAformer download).** `real_data.py download` now fetches Zenodo 17473983 (both flowcell files).
- **13 (notation).** ε and B_ε now map to `fixed_policy` / `metrics`, in both the notation table and S3.

**UNSUPPORTED**
- "statistically indistinguishable" → "differ by at most 0.4% with overlapping CIs (no paired test was run)".
- The coupon-collector remark was removed.
- "twice as many reads" → concrete values: 52% vs 27% at 1 read, 98% vs 87% at 3 reads.
- The clean-environment reproduction is now stated as covering an earlier commit, before the RQ1–RQ3 and confidence experiments. Its report is copied to `evidence/verification_reports/independent_reproduction_f8fb9a9.md`.
- 8,148: the D|N restriction is now stated, with the unrestricted figure (99.6% of 1,344 cases).
- Strand length: replaced by a measured value (168 nt for a 3 kB test archive). The 168–172 nt range was removed.

**OVERCLAIM**
- **Abstract and conclusion.** Rewritten:
  - The gain is now reported as a gain in mean cost on multi-container placements (11.6–13.9%), with single-container cases attributed to stopping only.
  - The stopping-versus-allocation decomposition is reported per pool.
  - The q99 budget is reported: not lower on Twist, +12% on Genscript round-robin. A q99 column was added to Table II.
- **Placement savings (34% / 63%).** Removed as headline numbers. The single-container layout is described as offering no failure protection and as chosen with the leaky belief.
- **False alarms.** PCR_60c (0.75) is now explicitly called a false alarm under our own 0.7 threshold.
- **Futility rule.** Now marked "simulated s7; not tested on real reads", and the causal "this explains why" was removed.
- **"Real reads".** Now "replayed / semi-empirical", with a paragraph on what is real and what is modelled. Decodability in the replay is stated to be an oracle event, and the CRC/SHA path to run only in the FASTQ simulation.
- **Oracle-static baseline.** Now "a local optimum; not deployable; static weights + verified stopping not tested".
- **"Never produced a wrong file".** Now "24 simulated retrievals of a 6 kB file; guarantee rests on SHA-256".
- **"Oracle channel <1%".** Now described as uninformative, because belief equals truth in s4.
- **Literature table.** The correlated-errors cell for this work is now P, with a footnote.
- **Surrogate.** Recomputed from `pred_cap@0.99` vs `cap@0.99`: −16 to +6% on Twist. On Genscript, round-robin is declared infeasible and the single-container budget is 5.0× too high.
- **ρ_U = 0.01.** Now states that the intervals become too wide on Twist (coverage 1.00).
- **Other wording.**
  - "kilobyte to hundreds-of-megabyte".
  - "reads carry errors".
  - "two pools from different vendors (one study)".
  - "decodable" instead of "provably decodable".
  - "plausibly less sensitive … consistent with but not isolated" instead of "robust".

**Math vs code**
- **Eq. 4.** The vacuous constraint and the sample-maximum caveat are stated, and q99 is reported.
- **Eq. 7.** The (1−ρ+h) factor is described as a heuristic, with its first-order form given. The ρ_i horizon assumption is stated.
- **Eq. 8.** It is now stated that the remaining budget is spread over all M containers.
- **Stopping path.** Verified stopping runs in the FASTQ pipeline; the replay uses the oracle event T ≤ R.
- **Seeding.** The hash-based seeding of the q3 null permutation is now mentioned in the reproducibility section.

**Citation usage**
- Bornholt/Organick: the containers motivation is now stated as our assumption.
- Checksum gating is now cited to Erlich.
- Nath and Leong wording softened.
- The "Birnbaum-type" name was removed.
- Cliff's δ was removed.
- The two-sided BMA variant is attributed to us.
- The Cauchy-RS citation to Reed–Solomon (1960) is kept as the RS origin. No Cauchy-specific reference was verified, so it remains a known gap.

**Not changed (and why)**
- The leak and hash seeding were disclosed rather than fixed. Fixing them requires regenerating all caches and re-running every real-data experiment, which was not done in this revision. This is listed in `missing_evidence.md`.
- Eq. 3's per-pair synthesis dropout and the π̄ approximation are documented in S3 rather than in the paper.
