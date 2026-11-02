# Supplementary S4: Methodology and Experimental Protocol

## Research questions and how each is tested
| RQ | Question | Experiment | Data | Primary metric | Statistic |
|---|---|---|---|---|---|
| RQ1 | Does retrieval generalise to an unfamiliar pool? | E10: run the controller on the Genscript pool with (a) in-pool calibration, (b) a prior built from the Twist pool, (c) a generic literature prior; prior sd 0.5 and 1.0 | real Genscript reads (60 permutations) | expected spend at 0.99; read share sent to the near-dead PCR_30c container | bootstrap 95% CI |
| RQ2 | Can it self-calibrate from limited observations? | E4: build the controller's belief from 0, 1k, 3k, 10k, 30k or 300k calibration reads per container (generic φ), or from full calibration (300k reads + per-run φ̂) | real Twist reads (100 permutations) | expected spend at 0.99 | bootstrap 95% CI; compared with oracle-tuned static |
| RQ3 | Does it detect changes in efficiency? | ED: posterior at the read count R* where the belief predicts 40% recovery, and at 2R* | simulated s7 (2 containers degraded); real Twist with belief = fresh pool | flag rate P(θ<0.7) > 0.9; posterior mean; 90% CI coverage | proportions over 200 trials |
| RQ4 | Does adaptive allocation reduce cost? | policy comparison | real Twist (200) and Genscript (200/100) | expected spend at 0.99 | paired bootstrap difference vs oracle-tuned static (2,000 resamples) |
| RQ5 | Is recovery exact under failures and heterogeneity? | failure injection (remove one container's stream); Q4 decode-safety runs | real Twist; simulated heavy-error channel | completion within 3× cap; false accepts; wrong outputs | counts; 95% upper bound for zero events (rule of three, 3/n) |
| — | Calibration of the posterior predictive | EC: P_complete at 30/50/70/85% of a budget set to the median adaptive cost | real Twist (200) | ECE (10 bins), Brier vs climatology, AUC | bootstrap CI for AUC |

## Baselines
- **Fixed depth.** Uniform read fractions. The budget is the 0.99-quantile of its own per-trial requirement on the evaluation data. This is favourable to the baseline, because the budget is known in hindsight.
- **Coverage-equalising static.** Fractions proportional to the expected per-fragment deficit under the belief. This is a heuristic.
- **Oracle-tuned static.** Fractions optimised by Nelder–Mead on the evaluation matrices themselves. This is an upper bound on any static plan.
- **Adaptive stop, round-robin.** Our verified stopping rule with uniform allocation. This isolates the stopping effect from the allocation effect.
- **Adaptive, coverage-greedy / pivotal.** The proposed controller with its two scores.

## Metric definition
- **Expected spend at 1−ε.** For a static plan, this is its budget B_ε, the (1−ε)-quantile of the per-trial minimal budget.
- **For an adaptive policy:** E[min(C, B_ε)] with B_ε the (1−ε)-quantile of its costs C.
- **Success probability:** reported separately (uncapped).
- **Bootstrap CIs:** resample trials. For paired differences, the same resample indices are used for both policies, because they share read permutations (common random numbers).

## Train/test separation
- **Twist and Genscript runs:**
  - The first 300k reads of each run form the calibration split. Beliefs (yield, NB shape, dropout, φ) come only from this split.
  - The next 900k reads form the evaluation split. Recovery-time matrices come only from this split.
- **E4:** the 1k–30k beliefs subsample the calibration split by binomial thinning.
- **E10, "Twist prior" condition:** uses only Twist calibration data.
- **Q3 nanopore models:** fitted on flowcell 1 and evaluated on held-out strands of flowcell 2.

## Multiple comparisons
- The headline comparison is pre-specified: adaptive pivotal vs oracle-tuned static, one per placement.
- Secondary comparisons are descriptive and are not corrected. We do not report significance stars for them.

## What is *not* done
- No new wet-lab synthesis or sequencing.
- No observed container-failure data. Failures are injected.
- No currency cost model.
- No regret analysis.
