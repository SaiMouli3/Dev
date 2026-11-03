# Supplementary S4: Methodology and Experimental Protocol (revision 3)

All work is computational. The frozen protocol is `experiments/configurations/dt_protocol.yaml`; changes made after freezing are listed there under `deviations`, and analyses that were not pre-specified are listed under `post_hoc_analyses`.

## Primary evidence: DT4DDS digital twin, end to end
- **Archive:** an 18,432-byte random file (seed 20261102), 24-byte payloads, Cauchy RS (K, N) = (32, 48), 25 stripes, 1,200 strands of 136 nt, six containers, round-robin placement.
- **Per container:** an independent DT4DDS 1.1.0 synthesis (Twist model, or the photolithographic model in one scenario), then aging, HiFi PCR and iSeq100 single-end 150 nt reads, with 60 reads per stored strand.
- **Decoding chain:** adapter trim, inventory clustering, consensus, inner RS(4) + CRC-16 + address check, outer RS with CRC-32 leave-one-out repair, SHA-256, then byte comparison.
- **No oracle:** true payloads are used only afterwards, to count wrong outputs and false accepts.
- **Design:** 30 independent syntheses × 10 read orders = 300 retrievals per policy and scenario. The photolithographic scenario uses 10 × 3 (deviation).
- **Scenarios:** baseline; fail1; fail2; homog (control); photolitho; and a sweep of the worst container at 200, 30, 10, 3 or 1 copies per strand.

## Policies and baselines
- **Verification-gated policies:**
  - round-robin;
  - coverage-greedy;
  - pivotal;
  - Thompson sampling;
  - UCB (posterior 0.9-quantile);
  - learn-then-commit (5 exploration batches per container);
  - information-greedy (largest posterior variance of ln θ until below 1e-3);
  - tail-aware (posterior 0.1-quantile and the 10% of stripes with the largest remaining need).
- **Static family:** read proportions w_j ∝ e_j^α, α ∈ {0, 0.5, 1, 2}.
  - e_j is measured in hindsight on the same scenario, so the best α is a hindsight-tuned ("oracle") static plan.
  - α = 0 is uniform fixed depth.
  - A static plan must reserve its q99, so its mean cost equals its q99.
- **Prior:** the generic prior for every container. No calibration data from the evaluated system is used.

## Metrics and statistics
- **Mean reads per retrieval:** the pay-per-read cost.
- **q99:** the reserved-capacity budget.
- **Success:** the output file is byte-identical to the stored one.
- **Safety:** wrong outputs, and false accepts per decode attempt, with exact Clopper–Pearson 95% upper bounds.
- **Confidence intervals:** cluster bootstrap over syntheses, 2,000 resamples.
- **Paired comparisons:** against round-robin + verified stop, on retrievals where both policies succeeded.
- **Reporting:** no significance stars; every pre-specified comparison (C1–C4) is reported.

## Calibration and forecasts
- **Choosing ρ_U:**
  - selected from {0, 0.003, 0.01, 0.03, 0.1} on disjoint replicates 1000–1009;
  - the chosen value is the one whose 90% credible-interval coverage of θ* is closest to 0.90;
  - it is then evaluated on held-out replicates 0–29.
- **Forecasts:** the analytic forecast is compared with posterior simulation (200 draws).
  - The budget is the median coverage-greedy cost.
  - Checkpoints are at 30/50/70/85% of the budget.
  - Metrics: ECE (10 bins), Brier score against climatology, and AUC with cluster-bootstrap CIs.
- **Futility:** disabled throughout.

## Secondary evidence: real-read replay (leak-free)
- **Data:** ENA PRJEB65931, Twist and Genscript, six runs each; each run emulates one container.
- **Decodability:** an oracle event (exact consensus).
- **Prior:** the generic prior is primary.
  - The optional calibrated belief uses only the calibration split (the first 300k reads).
  - Evaluation uses the next 900k reads.
- **Design:** 200 read permutations.
- **Earlier leaky numbers:** kept only as an upper bound (manuscript Appendix A).

## What is not done
- No wet lab.
- No observed failure rates.
- No currency cost model.
- No regret or optimality theory.

## Revision 2 protocol (superseded; real-read replay with the leaky calibrated belief)

### Research questions and how each is tested
| RQ | Question | Experiment | Data | Primary metric | Statistic |
|---|---|---|---|---|---|
| RQ1 | Does retrieval generalise to an unfamiliar pool? | E10: run the controller on the Genscript pool with (a) in-pool calibration, (b) a prior built from the Twist pool, (c) a generic literature prior; prior sd 0.5 and 1.0 | real Genscript reads (60 permutations) | expected spend at 0.99; read share sent to the near-dead PCR_30c container | bootstrap 95% CI |
| RQ2 | Can it self-calibrate from limited observations? | E4: build the controller's belief from 0, 1k, 3k, 10k, 30k or 300k calibration reads per container (generic φ), or from full calibration (300k reads + per-run φ̂) | real Twist reads (100 permutations) | expected spend at 0.99 | bootstrap 95% CI; compared with oracle-tuned static |
| RQ3 | Does it detect changes in efficiency? | ED: posterior at the read count R* where the belief predicts 40% recovery, and at 2R* | simulated s7 (2 containers degraded); real Twist with belief = fresh pool | flag rate P(θ<0.7) > 0.9; posterior mean; 90% CI coverage | proportions over 200 trials |
| RQ4 | Does adaptive allocation reduce cost? | policy comparison | real Twist (200) and Genscript (200/100) | expected spend at 0.99 | paired bootstrap difference vs oracle-tuned static (2,000 resamples) |
| RQ5 | Is recovery exact under failures and heterogeneity? | failure injection (remove one container's stream); Q4 decode-safety runs | real Twist; simulated heavy-error channel | completion within 3× cap; false accepts; wrong outputs | counts; 95% upper bound for zero events (rule of three, 3/n) |
| — | Calibration of the posterior predictive | EC: P_complete at 30/50/70/85% of a budget set to the median adaptive cost | real Twist (200) | ECE (10 bins), Brier vs climatology, AUC | bootstrap CI for AUC |

### Baselines
- **Fixed depth.** Uniform read fractions. The budget is the 0.99-quantile of its own per-trial requirement on the evaluation data. This is favourable to the baseline, because the budget is known in hindsight.
- **Coverage-equalising static.** Fractions proportional to the expected per-fragment deficit under the belief. This is a heuristic.
- **Oracle-tuned static.** Fractions optimised by Nelder–Mead on the evaluation matrices themselves. This is an upper bound on any static plan.
- **Adaptive stop, round-robin.** Our verified stopping rule with uniform allocation. This isolates the stopping effect from the allocation effect.
- **Adaptive, coverage-greedy / pivotal.** The proposed controller with its two scores.

### Metric definition
- **Expected spend at 1−ε.** For a static plan, this is its budget B_ε, the (1−ε)-quantile of the per-trial minimal budget.
- **For an adaptive policy:** E[min(C, B_ε)] with B_ε the (1−ε)-quantile of its costs C.
- **Success probability:** reported separately (uncapped).
- **Bootstrap CIs:** resample trials. For paired differences, the same resample indices are used for both policies, because they share read permutations (common random numbers).

### Train/test separation
- **Twist and Genscript runs:**
  - The first 300k reads of each run form the calibration split. Beliefs take yield, NB shape and dropout from this split.
  - **Exception (leak):** the belief's consensus curve φ_j is pooled from per-strand curves estimated on evaluation-split reads. See S5 for the leak-free conditions.
  - The next 900k reads form the evaluation split. Recovery-time matrices come only from this split.
- **E4:** the 1k–30k beliefs subsample the calibration split by binomial thinning.
- **E10, "Twist prior" condition:** uses only Twist calibration data.
- **Q3 nanopore models:** fitted on flowcell 1 and evaluated on held-out strands of flowcell 2.

### Multiple comparisons
- The headline comparison is pre-specified: adaptive pivotal vs oracle-tuned static, one per placement.
- Secondary comparisons are descriptive and are not corrected. We do not report significance stars for them.

### What is *not* done
- No new wet-lab synthesis or sequencing.
- No observed container-failure data. Failures are injected.
- No currency cost model.
- No regret analysis.
