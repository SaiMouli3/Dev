# Scientific Quality Report: Reviewer-Style Assessment

**Manuscript:** "Self-Calibrating Adaptive Sequencing for Exact Retrieval from Heterogeneous Multi-Container DNA Storage: Method and Evaluation on Real Read Streams"

- **Audited version:** `research_paper/manuscript/main.tex` at commit `89e439f`. Line numbers below refer to that version.
- **Companion files:**
  - `audit/claim_verification.md`: the claim-by-claim table.
  - `audit/citation_verification.md`, section "Usage audit (manuscript)".

---

## 1. Summary of the submission

The paper studies retrieval of an erasure-coded file (Cauchy RS, K = 32, N = 48) stored across M containers whose sequencing efficiency is unknown and heterogeneous.

**Method.** A controller does the following:
- keeps a discrete posterior over a per-container efficiency multiplier θ_j, including a failure atom;
- maps it to per-pair recovery hazards through a Gamma–Poisson / log-normal coverage model with an empirical consensus curve φ_j;
- reads the container with the best (stripe-weighted) hazard per cost;
- stops when every stripe has K verified fragments and the SHA-256 matches.

**Evaluation.** Public Illumina runs (ENA PRJEB65931; Twist and Genscript pools, six handling conditions each) are turned into "empirical recovery-time matrices". Each run emulates one container, and policies are replayed over 60–200 read permutations.

**Headline results:**
- 11–17% fewer reads than a hindsight-tuned static plan (Twist) and 11% (Genscript);
- 34% / 63% fewer than fixed depth when combined with placement;
- the generic-prior controller is about as good as a calibrated one.

**Negative results:**
- pivotal weighting gives no measurable gain;
- credible intervals are over-confident;
- the completion forecast is under-confident;
- the futility rule is harmful;
- the surrogate model is biased.

**Supporting results:** a nanopore correlated-error model, and a placement threshold that is a special case of Leong et al. (2012).

The paper is unusually candid about novelty and negative results, and every number in it can be traced to a committed file. The core empirical claim is numerically correct, but it is framed more strongly than the evidence supports (see W1–W4). Several methodological issues also affect its interpretation: evaluation-split leakage, an internal inconsistency on Genscript, and an extreme-quantile metric with small trial counts.

---

## 2. Strengths

1. **Honest positioning.**
   - Stopping-when-decodable is credited to prior work (Lopez 2019).
   - The placement rule is credited to Leong et al.
   - Multi-container coverage and container failure are credited to Kong et al. and Levy et al.
   - No "first" claim is made. The hedged novelty statement (l. 113–115) is appropriate.
2. **Negative results are reported prominently:** pivotal ≈ coverage-greedy, over-confident posteriors, an under-confident P̂_complete, a harmful futility rule, and surrogate bias. This is rare and valuable.
3. **Strong reference baselines.**
   - The static baseline is tuned in hindsight on the evaluation matrices.
   - The fixed-depth baseline is given its hindsight budget.
   - All policies share common random numbers (identical read permutations).
4. **A reusable evaluation idea.** Empirical recovery-time matrices from public runs could be reused by others, and the closed-form B*(w) for static plans is neat.
5. **A second pool.** The Genscript replication comes from a different vendor and design, with a harsher channel, and includes a near-dead container. Cross-pool transfer (RQ1) is tested with beliefs from another pool and from a generic prior.
6. **Traceability.**
   - The claim-to-evidence CSV, the SHA-256 manifest of result files (all verified OK) and the figure scripts reading only committed results are good practice.
   - Results are byte-identical between `results/` and `evidence/experiment_results/`.
7. **Clean verification-gated decoding pipeline.** Inner RS + CRC-16 + address check, then per-stripe CRC-32 with leave-one-out repair, then SHA-256. This is implemented and tested (202 passing tests).

---

## 3. Weaknesses

### 3.1 Major

**W1. The headline saving is mostly verified stopping, and it is reported only on expected spend, not on the reliability budget.**
- On Twist round-robin, adaptive *stopping* with round-robin allocation already saves 5.4k of the 7.3k reads (63.1k → 57.8k). Allocation adds 2.0k (about 3.5 percentage points).
- The upper end of "11–17%" is the single-container placement. There, allocation is trivial and "static tuned in hindsight" equals fixed depth.
- On Genscript the picture reverses: stopping alone saves 9.0k (CI includes 0) and allocation saves 64k. This decomposition is the paper's most informative result but is not stated in the abstract or conclusion.
- The adaptive policy's 0.99-quantile cost (`cap@0.99`) is not better than the static budget:

| Pool / placement | Adaptive cap@0.99 | Oracle static budget |
|---|---|---|
| Twist RR | 63.2k | 63.1k |
| Twist tiered | 53.0k | 53.1k |
| Twist single | 54.0k | 53.9k |
| Genscript RR | **974k** | **651k** (+50%) |

- Anyone who must provision a budget for 0.99 reliability (e.g. buy a flowcell) gains nothing, and on Genscript loses. Report both expected spend and the 0.99 budget, and explain the trade-off.

**W2. A missing baseline: hindsight-tuned static weights *with* verified stopping.**
- The paper compares an adaptive-allocation + adaptive-stopping policy against a static-allocation + fixed-budget plan, and against uniform allocation + adaptive stopping.
- The natural deployable competitor is oracle-tuned (or calibration-tuned) static weights + verified stopping. It is missing, so the value of *online* allocation over a good fixed allocation is not isolated.
- Other expected references:
  - an explore-then-commit heuristic (probe, then drop slow containers);
  - Thompson sampling (implemented in `adaptive_policy.py` but never evaluated);
  - the clairvoyant lower bound (implemented but not reported), which would show the optimality gap;
  - the Kong et al. best-partition static reference, promised in literature_review.md §4 but not used.

**W3. Evaluation-split leakage, not disclosed in the manuscript.**
- `belief_from_calibration` (`experiments/real_data.py`) builds each container's believed consensus curve from `rd.phi_strand[::2]`.
- `phi_strand` is estimated from **evaluation-split** reads (`load_run` stores evaluation reads in `reads_by_ref`; `estimate_phi` uses them).
- The manuscript says beliefs come from a disjoint calibration split "except where stated" (l. 414–415), but no exception is stated. `supplementary/limitations.md` l. 19 asserts the opposite; `reports/real_data_validation.md` l. 27 admits "a small leak".
- **Affected results:** Table II (adaptive rows), E4 "full calibration", EC, ED (fresh-pool φ), E10 in-pool, and the choice of the "single best" container.
- The generic-prior results (E4, E10) are leak-free and give similar costs, so the cost impact is probably small. The change-detection and confidence results are not covered by that argument.

**W4. Internal inconsistency on Genscript, larger than the claimed effects.**
- In-pool calibrated adaptive retrieval on Genscript round-robin costs 577.7k [560.6–595.2k] in Table II (200 permutations) but 528.9k [513.9–544.4k] in Table `transfer` (60 permutations).
- The CIs do not overlap, and the static references also shift (650.7 → 621.0k; 831.2 → 798.2k).
- The 9% gap exceeds the 4–5% "price of an unfamiliar pool" (RQ1) that the paper reports.
- Likely cause, from file timestamps and git history: the Genscript per-strand φ̂ cache was deleted (commit 0c322f9) and regenerated with seeds from Python's process-randomised `hash(run)` before E10 and the Genscript hold-out ran.
- This also means "small run-to-run differences" (l. 725) understates the reproducibility problem. Re-run Table II and E10 from one cache with a stable seed before relying on cross-experiment comparisons.

**W5. The semi-empirical replay is presented as "real reads".**
Real are: read arrival order, i.e. a permutation of 900k mapped reads, and per-strand abundance in each run. Not real are:
- (a) whether a cluster's consensus succeeds: a random draw t_ij ~ φ̂_ij from a curve estimated from ≤ 24 reads with only 2 subsamples per n (`strand_phi`);
- (b) the containers: emulated by restricting a full-pool run to a subset of strands, plus a proportional share of junk reads. The paper calls this "exactly" the container's stream (l. 403), but it is an assumption that subsetting leaves relative amplification and sampling unchanged;
- (c) trial-to-trial variance: only read order varies. Library prep, cluster generation and per-strand abundance are fixed in every permutation, so the bootstrap CIs reflect order randomness only, not run-to-run variability;
- (d) the "verified" stopping in the replay: an oracle (consensus == design). The CRC/SHA path is exercised only in the 24-run simulation.

Please describe the methodology as semi-empirical, and discuss how (c) narrows the CIs.

**W6. Container heterogeneity is limited and partly artificial.**
- The six "containers" per pool are one synthesis batch under different handling, sequenced on two flowcells. Relative efficiencies on Twist are 0.34–1.03.
- The method's value depends on heterogeneity and failures that are either mild (Twist), extreme in one container (Genscript PCR_30c) or injected.
- Container-failure rates are not measured. The paper acknowledges this (l. 699–703); still, the external validity of the cost claims is narrow, and the abstract should say "handling conditions of one synthesis batch" rather than imply independent containers.

**W7. The objective (Eq. 4) is ill-posed as written.**
- B_ε is defined as the (1−ε)-quantile of C_π. The constraint P(exact recovery within B_ε) ≥ 1 − ε then holds by construction for any policy that eventually succeeds, so the problem reduces to min E[min(C_π, B_ε(π))].
- This is not the stated chance-constrained problem. It also rewards heavy-tailed policies, since the worst 1% is truncated at the policy's own quantile (cf. Genscript, cap 974k).
- **Suggested formulation:** min_π E[C_π] subject to P(C_π ≤ B) ≥ 1 − ε for a *given* B. Or report the Pareto front (B, E[C]).

**W8. Statistical treatment of an extreme quantile with 60–200 trials.**
- The 0.99-quantile from 200 trials is the 2nd-largest value; from 60 or 100 trials it is the maximum (`budget_for_target`).
- The oracle static weights are optimised to minimise this order statistic on the same trials. They are not re-optimised inside the bootstrap, so the bootstrap ignores the selection step.
- The percentile bootstrap is known to be inconsistent for extreme order statistics. Symptom: Δ = −7.3k has CI (−12.1, −7.0), with the point estimate at the edge of its own interval.
- The paper should:
  - use more permutations (at least 1,000);
  - or report the 0.95 target alongside 0.99;
  - and/or use a m-out-of-n or subsampling bootstrap;
  - and nest the oracle optimisation inside the resampling.
- E4 and E10 compare configurations through overlapping *unpaired* CIs although paired comparisons on the same matrices are available.
- No multiplicity control is applied across the many placements, configurations and checkpoints. Holm was excluded for citation reasons; Bonferroni would do.
- Cliff's δ is announced (l. 451) but never reported.

### 3.2 Minor

- **m1. Number errors** (see claim_verification §3.1):
  - 10,567 → 10,557 decode attempts;
  - 54–60 → 45–65 false accepts per weak-code retrieval;
  - 1 of 12 weak-code retrievals "failed", not "partial";
  - ablation +79% / +68% use a different metric from the figure (+69.6% / +67.0%);
  - "most reliable in 3/5 seeds" → 2/5;
  - P̂_complete mean 0.04–0.22 → 0.001–0.22;
  - 99.5% mapped → 98.2–99.7% across pools;
  - `real_data.py download` does not fetch Zenodo 17473983.
- **m2. Ablation (simulated s4):**
  - The comparisons use mean reads of *successful* retrievals across variants with different success probabilities (0.9525–0.999), which favours less reliable variants.
  - "Perfect channel knowledge" is uninformative because the belief equals the truth in s4.
  - The full design's 0.990 in the figure is the single lucky seed. Across five seeds it never reaches 0.99.
- **m3. Futility rule.** It is evaluated only in simulated s7, but is presented as a general negative result (abstract l. 66–67, l. 707–708). The causal link to the under-confident P̂_complete on real Twist streams (l. 547–548) is asserted, not shown.
- **m4. The overdispersion fix is chosen post hoc.**
  - ρ_U = 0.01 is chosen on Twist + sim.
  - On Twist it yields 100% coverage for nominal 90% intervals, i.e. too wide. That is the same criterion used to reject ρ_U = 0.05.
  - The Genscript hold-out (n = 5 containers) is the only out-of-sample evidence.
- **m5. Change detection.**
  - The flag threshold (θ < 0.7, P > 0.9) is arbitrary.
  - A container at 0.75 efficiency is flagged 73–82% of the time. By the script's own ground truth (`truly_degraded = eff < 0.7`) this is a false-alarm rate, yet the abstract claims "no false alarms on near-nominal containers".
  - Detection is evaluated statically (posterior at a fixed R* with only container j read), not inside the closed loop.
  - "True relative efficiency" is a ratio of median recovery times; for Genscript PCR_30c it is 4.87 (meaningless).
- **m6. Practicality.** A batch of 200 reads with zero setup cost (≈ 280 decision rounds per retrieval) is unrealistic for Illumina. The paper acknowledges this (iv), but no sensitivity to batch size, setup cost or horizon H is reported on real data, so the savings could vanish under realistic flowcell granularity.
- **m7. The "generic" prior is not fully generic.** It uses σ_g = 0.3, close to the Twist estimate (0.31). This is not stated.
- **m8. Nanopore correlated-error claim.**
  - The cross-flowcell strand-specific correlation (r = 0.57) is measured on the same synthesised molecules. It may reflect synthesis variants or reference errors (14% design mismatch), not correlated sequencing errors.
  - The 2.2× MAE improvement comes from 3 coverage values only (n = 3, 5, 10).
  - The model is not coupled to the controller, so the comparison-table entry "Correlated errors: Y" overstates this.
- **m9. Clarity.**
  - Symbols are overloaded: h_ij is both abundance and hazard; ρ is the ICC, ρ_i, ρ_par and ρ_U; π is both the policy and the read probability.
  - Eq. (3) uses π̄ while the model uses π with g normalisation; the approximation is not discussed.
  - The figure "Reads needed for exact recovery with probability ≥ 0.99" (Fig. `cost`) actually plots expected spend for the adaptive bars.
  - Statements about decodability in the real-data replay (l. 679 "provably") overstate an oracle check.
- **m10. Length.** The manuscript is about 11 pages in IEEEtran conference format; most IEEE conferences allow 5–8. Move the overdispersion study, the placement rule and parts of the channel section to supplementary material once a venue is chosen.
- **m11. Unreported sources.** Five of eight simulated containers are omitted from Table `detect`, and the Srinivasavaradhan dataset, although listed in Data (l. 440), contributes no reported result.

---

## 4. Statistical review

| Aspect | Assessment |
|---|---|
| **Unit of replication** | One read permutation of a fixed run. Trials are not independent sequencing experiments. CIs capture ordering randomness only (W5c). |
| **Sample sizes** | 200 (Twist, Genscript RR), 100 (Genscript single-best, E4), 60 (E10, failure injection). The 0.99-quantile equals the maximum for n ≤ 100. This is too few for the target statistic (W8). |
| **Paired design** | Correctly used for adaptive vs oracle static (`_paired_spend_diff`). Not used for pivotal vs coverage, E4 configurations or E10 beliefs, although it was available. |
| **Bootstrap** | Percentile, 2,000 resamples. For extreme quantiles it is inconsistent, and the oracle tuning is not nested. The adaptive `mean_spend` CI (`metrics.adaptive_summary`) holds the cap fixed, so cap uncertainty is ignored. |
| **Effect sizes** | Cliff's δ is promised but absent. Relative differences are reported without CIs, except via Δ. |
| **Multiplicity** | None. At least 6 placements × 3 policies, 14 E4 configurations, 4 checkpoints × 2 beliefs, and 14 containers × 2 operating points in ED. |
| **Calibration metrics** | ECE (10 bins, 200 per checkpoint), Brier against climatology, and AUC with CIs are appropriate and correctly recomputed. Bins with n < 5 are dropped in the figure (stated). |
| **Rare-event bounds** | Rule of three, 3/10,557 = 2.8 × 10⁻⁴ per decode attempt. The more relevant denominator is wrong consensus words (2,682 → 1.1 × 10⁻³). 0 wrong files in 24 runs gives a 95% upper bound of about 12% on the per-retrieval wrong-file rate. "Never" is not justified. |
| **Post-hoc choices** | ρ_U was chosen on the selection set (acknowledged); the flag threshold of 0.7 and R* at 40% are not justified; the EC budget is the median adaptive cost from the same run (circular, but declared). |
| **Seed robustness** | One seed for every real-data experiment. The simulated seed sweep (q5) shows that conclusions can flip (the "3/5" claim is actually 2/5). |

---

## 5. Reproducibility review

**Positives**
- Fixed seed offsets.
- A committed result-file SHA-256 manifest (verified OK).
- Figure scripts read only committed files.
- 202 passing tests (`evidence/test_reports/pytest.txt`).
- `experiments/reproduce.py`: its Theorem-1 check was re-run by this audit and gives 8,148/8,148 PASS.
- Claim-to-evidence CSV covering 30 claims.

**Problems**
1. **Process-randomised seeding.**
   - `real_data.py` l. 88 seeds φ̂ estimation with `hash(run) % 2**32`.
   - `q3_correlated_errors.py` l. 68 uses `hash()` to order the null permutation. This second use is not mentioned in the paper.
   - With caches deleted and regenerated, headline Genscript costs shift by about 9% (W4).
   - Fix: use a stable hash (`zlib.crc32`), set PYTHONHASHSEED in all runners, re-run everything, and commit the φ̂ caches or their checksums.
2. **Unpinned dependencies.** `requirements.txt` uses only lower bounds (`numpy>=1.26`, …). Provide a lock file or `pip freeze`, and record the Python version used (≥ 3.11 stated; 3.13 used in the reproduction).
3. **Data acquisition is incomplete.** `real_data.py download` does not fetch Zenodo 17473983, which Q3 needs, although l. 719–721 says it does.
4. **The clean-environment reproduction report is not in the repository.** critical_questions.md points to `scratchpad/repro_report.md`. It was run on commit f8fb9a9, before E4, EC, ED and E10 existed, so the RQ1–RQ3 results have never been independently reproduced.
5. **Compute cost.** The reported runtimes total about 7–9 h on 4 cores and the joint optimiser was skipped on Genscript for compute reasons. That is reasonable, but a "smoke" configuration with fewer trials and published tolerances would help reviewers.
6. **Evidence documents disagree.**
   - `supplementary/limitations.md` l. 19 (φ̂ from calibration reads) vs `real_data_validation.md` l. 27 (leak).
   - README and reports still say "3 of 5 seeds".
   - `evidence_traceability.csv` C14 says 10,567.

---

## 6. Math vs code (Eqs. 1–9)

Equation numbers follow the order in `main.tex`:

| Eq. | Label | Content |
|---|---|---|
| (1) | `eq:abundance` | read probability |
| (2) | `eq:T` | recovery time |
| (3) | `eq:Q` | pair decodability |
| (4) | `eq:objective` | objective |
| (5) | `eq:posterior` | posterior |
| (6) | `eq:hazard` | hazard |
| (7) | `eq:score` | allocation score |
| (8) | `eq:pcomplete` | completion probability |
| (9) | `eq:betabin` | beta-binomial errors |

**Code checked:**
- `src/d3ra/channel/surrogate.py` (`sample_T`, `SuccessTable`, `Layout.pi_bar`)
- `src/d3ra/retrieval/adaptive_policy.py` (`_QTable`, `posterior`, `_hazard`, `scores`, `p_complete`, `run`)
- `src/d3ra/retrieval/stopping_rule.py::completion_probability`
- `src/d3ra/channel/correlated_errors.py`
- `src/d3ra/realdata/empirical.py`
- `src/d3ra/retrieval/fixed_policy.py::budget_for_target`
- `src/d3ra/evaluation/metrics.py`

| Eq. | Paper | Code | Finding |
|---|---|---|---|
| (1) | π_ij = y_j x_ij g_i h_ij / Σ x g h; g ~ LogN(−σ²/2, σ²), h ~ Γ(k, 1/k) | `sample_T`: `lognormal(-sg**2/2, sg)`, `gamma(k, 1/k)`, `pi = yield*w/tot` | **Match.** x_ij may be a copy count (`pair_copies`), not only {0, 1} as l. 237 states. |
| (2) | T_ij = Γ(t_ij, 1)/π_ij; ∞ if absent or failed | `T = gamma(t,1)/pi`; dropout per pair (d_j) and per fragment (d_syn); failure per container | **Match.** For real data `empirical_T` replaces it with stream positions, and t_ij is drawn independently of *which* reads arrive (consensus success is not tied to the actual reads). |
| (3) | Q_j(μ) = (1−d_j) E_g[Σ_n NB(n; k_j, gμ) φ_j(n)], 9-node Gauss–Hermite | `SuccessTable` + `_QTable`: NB truncated at n = 48 with the tail assigned φ(48), monotonised by `maximum.accumulate`, multiplied by (1−d_j)(1−d_syn) | **Match** up to the stated d_syn factor (d_syn = 0 in all real-data beliefs). (a) Applying the *shared* synthesis dropout per pair treats copies in different containers as independent, contrary to the text's own description; harmless with one copy. (b) Eq. (3) evaluates at π̄_ij (layout-average share), ignoring the g normalisation in Eq. (1); undocumented approximation. |
| (4) | min_π E[min(C_π, B_ε)] s.t. P(exact within B_ε) ≥ 1−ε; B_ε = (1−ε)-quantile of C_π | `metrics.adaptive_summary`: `cap = quantile(c, t)`, `spend = mean(min(c, cap))`; static: spend = budget | Code matches the formula, but **the constraint is vacuous** (it holds by construction once B_ε is the policy's own quantile). See W7. `budget_for_target` uses the ⌈t·n⌉-th order statistic, i.e. the sample max for n ≤ 100 at t = 0.99. |
| prior | p_0(0) = f_j; p_0(θ^(ℓ)) ∝ (1−f_j) exp(−½(ln θ/s)²) on 25 log-spaced points in [0.01, 5] | `THETA` (l. 41), `log_prior` (l. 111–115) | **Match.** f_j = 10⁻³ in all real-data beliefs; not stated in the paper. |
| (5) | p(θ \| R, U) ∝ p_0(θ) Bin(U_j; n_j, Q_j(θ R_j π̄_j)); p(0) = 0 if U_j > 0 | `posterior`: `mu = R*pi_cont*THETA` with `pi_cont = y_j/Σ x`, binomial log-likelihood, `ll[...,0] = -inf if U>0`; optional beta-binomial with a = q(1−r)/r, b = (1−q)(1−r)/r | **Match.** n_j counts distinct fragments while `pi_cont` divides by total copies; identical with one copy. The beta-binomial option ("ρ_U") is described in text only; consider writing its likelihood. |
| (6) | h_ij(Δ) = [Q(θ̂(R+Δ)π̄) − Q(θ̂Rπ̄)] / [1 − Q(θ̂Rπ̄)], θ̂ = posterior mean | `_hazard` with `theta = post @ THETA` | **Match** (clipped to [0, 1]). |
| W_s | μ_s = Σρ_i, σ_s² = Σρ(1−ρ) + ¼, z_s = (μ_s − ν_s + ½)/σ_s, W_s = φ(z)/(σΦ(z)) | `scores` l. 166–177 | **Match.** ρ_i is computed with every container receiving H = 6 batches simultaneously, not one container at a time; the text's "within a look-ahead horizon H" should say "if every container received H batches". |
| (7) | S_j = Σ_open W_s h_ij(b)(1 − ρ_i + h_ij(b)) / (c_j b + c_j^setup) | `contrib = W*h1*(1-rho+h1)`; cost = c·b + setup | **Match**, but the factor (1 − ρ_i + h) is unexplained. The code comment says the marginal effect is h(1 − ρ_i^{−j}) ≈ h(1 − ρ_i), so the "+h" term is an ad-hoc correction. Justify or drop it. The exact Poisson-binomial `pivotal_weights` in `surrogate.py` is not used by the controller. |
| cov.-greedy | S_j = Σ h_ij(b)/(c_j b + c^setup) over fragments of unfinished stripes | `scores` l. 161–162 | **Match.** |
| (8) | P̂_complete = ∏_s Φ((Σρ^rem − ν_s + ½)/√Σρ^rem(1−ρ^rem)); remaining budget spread over containers | `p_complete` + `completion_probability` (+10⁻⁹ in the variance) | **Discrepancy:** the remaining budget is spread over **all M containers** (`np.repeat(rem/cost.sum(), L.m)`), including containers that store no fragments (e.g. tiered, single-best). Those layouts systematically under-estimate completion, a candidate cause the paper does not list. EC uses round-robin, where all containers are used, so the published EC numbers are unaffected. Also, P̂ uses the plug-in θ̂ and ignores posterior uncertainty; together with stripe independence this should cause over-confidence, so the observed under-confidence is probably the even-spread assumption. |
| (9) | p_sp ~ Beta(μ(1−ρ)/ρ, (1−μ)(1−ρ)/ρ); context = 5-mer + homopolymer run; MoM E[(x̄−μ)²] = μ(1−μ)/n + ρμ(1−μ)(1−1/n) | `HierarchicalErrorModel.strand_propensities`, `context_key(k=2)`, `estimate_icc` | **Match.** Undocumented details: homopolymer run capped at 4; contexts with < 50 observations fall back to the global rate; ICC is a pooled ratio estimator truncated at 0; strands with < 4 reads are skipped. |
| placement rule | t = ⌊Dρ_par⌋; q* root of P(Bin(D, q) ≥ t+1) = q | `theory/spreading.py`; recomputed q* = 0.232408 | **Match.** The theorem requires D \| N (even split); this restriction is stated in critical_questions.md but not in the paper. |
| stopping (§IV-C) | Stop only after inner RS + CRC-16 + address check, stripe CRC-32 and SHA-256; leave-one-out repair | `decoding/reconstruction.py`, `sequencing_env.py` | **Implemented**, but only in the FASTQ/simulation path. In the real-data replay (`AdaptiveController.run` with an empirical T), "verified" recovery is the oracle event T_ij ≤ R_j (consensus equals design), and no CRC/SHA logic runs. Say so in §VI. |
| run loop | — | `run` l. 230–237: trials whose *realised* T makes recovery impossible are cut short and charged the cap | Uses oracle knowledge of T for accounting only. It does not change success rates, but it does affect cost when the cap is finite (failure injection). Document it. |
| replay φ̂ | φ̂_ij from repeated consensus of n reads | `strand_phi`: forward-only BMA, `draws=2`, early stop after two consecutive exact successes, forward-fill | Each grid value is an average of only 2 subsamples. The early stop sets φ = 1 for all larger n after 2 × 2 successes. The paper describes forward-only consensus for raw reads (l. 392–393) but not these estimator details. |

---

## 7. Threats to validity not (fully) discussed

1. **Leakage** of evaluation reads into the believed φ_j (W3).
2. **Permutation-only variance.** CIs ignore library-prep and run-to-run variability (W5c).
3. **Emulated containers.** Restricting a full-pool run to a subset assumes subsetting is neutral (W5b).
4. **Seed and cache dependence** of headline Genscript numbers (W4).
5. **Granularity of sequencing.** 200-read batches and zero setup cost (m6). No sensitivity analysis.
6. **Choice of operating point.** The 0.99 target with ≤ 200 trials makes conclusions hinge on 1–2 extreme permutations (W8).
7. **Selection of the single-best container** uses the (leaky) belief, and its failure risk is ignored in the "combined with placement" headline numbers.
8. **Manifest dependence.** Exactness requires the out-of-DNA manifest (CRCs, SHA-256, parameters) to survive. Not stated.
9. **Generic prior's σ_g** close to the Twist estimate (m7).

---

## 8. Recommendation

**Major revision.**

The work is careful, transparent and potentially useful, especially the empirical recovery-time-matrix methodology and the explicit negative results. Before acceptance the authors must:
- (i) disclose or remove the evaluation-split leakage;
- (ii) resolve the Genscript inconsistency by re-running with stable seeding;
- (iii) report the 0.99-budget (quantile) comparison alongside expected spend, and decompose the savings into stopping vs allocation in the abstract;
- (iv) add the "tuned static weights + verified stopping" baseline;
- (v) fix the objective formulation and the extreme-quantile statistics (more permutations, nested bootstrap, paired tests);
- (vi) correct the number mismatches listed in `claim_verification.md` §3.1 and soften the overclaims in §3.3.

The supporting nanopore and placement results are adequately positioned but should be shortened to fit a conference page limit.

---

## 9. Process notes

- **Template.** The official IEEE template download from ieee.org returned an empty response to automated download. The manuscript therefore uses the CTAN copy of the official IEEEtran class v1.8b (`manuscript/IEEEtran.cls`, `IEEEtran.bst`), with `bare_conf.tex` as the structural reference (`manuscript/ieee_template/`).
- **Venue.** No venue is targeted, so page limits and venue-specific formatting were not checked; the current length of about 11 pages is noted in m10.
- **Author block.** The affiliation is a placeholder ("[Affiliation to be completed by the authors]", city/country and e-mail; l. 30–33) and must be completed by the author.
- **Scope of this audit.**
  - Read-only. No manuscript, table, code or result file was modified.
  - The only computation run was the Theorem-1 part of `experiments/reproduce.py` (8,148/8,148 PASS, under 2 min), plus pandas recomputations of statistics from committed result files.
  - Long experiments were not re-run.
  - `main.tex` changed during the audit (test count 201 → 202, commit 89e439f); line numbers refer to that version.
