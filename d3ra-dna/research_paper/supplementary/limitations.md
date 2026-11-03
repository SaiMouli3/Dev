# Supplementary S5: Limitations and Threats to Validity (revision 3)

## External validity
- **Computational only.**
  - The primary evidence is the DT4DDS digital twin with its default parameters.
  - Conclusions are bounded by DT4DDS's fidelity and by one archive geometry: 1,200 strands of 136 nt, (32, 48), six containers.
- **Failures are simulated.** A failed container returns non-informative reads. Real failure modes and rates are unknown.
- **Replay data.** The replay uses real reads but oracle decodability. Its containers are runs of one synthesis batch.

## Internal validity
- **ρ_U selection:**
  - ρ_U was selected with prior sd s = 1 but used with s = 0.5, a mismatch noticed after freezing.
  - The selected ρ_U = 0.1 is the largest grid value, and its selection coverage of 0.87 is below 0.90.
- **Static plans** use efficiencies measured in hindsight on the evaluated syntheses, which favours them.
- **Deviations:**
  - the photolithographic scenario uses a reduced 10 × 3 design;
  - the static allocator was fixed after freezing (it affects failed retrievals only);
  - one truncated replicate was regenerated;
  - the first Twist replay ran with ρ_U = 0 (see `dt_protocol.yaml`).

## Construct validity
- Reads are an abstract cost: no set-up cost and no flowcell granularity.
- The mean cost serves pay-per-read users. Users who must reserve capacity pay the q99, which no rule lowers.

## Statistical conclusion validity
- Safety bounds come from zero events and scale as about 3/n.
- Comparisons are paired, but all within one archive design.
- In fail1, comparisons use the 102 jointly successful retrievals.

## Measured weaknesses of the method
| Weakness | Evidence |
|---|---|
| A single efficiency multiplier cannot represent dropout; every model-based rule is costlier than round-robin with a near-dead container | `results/dt_baseline.csv`, `results/dt_sweep_cov*.csv` |
| Better-calibrated likelihood (ρ_U = 0.1) gives costlier allocation than the binomial (post hoc) | `results/dt_*_posthoc_rho0.csv` |
| Held-out 90% CI coverage 0.87 overall, 0.20 for the near-dead container | `results/dt_calibration_eval_baseline.csv` |
| Completion forecasts over-confident; Brier worse than climatology | `results/dt_forecast_baseline.csv` |
| No rule, including tail-aware, lowers q99 | `results/dt_summary.csv` |
| Information-greedy spends 66% of reads on the near-dead container (+441%) | `results/dt_baseline.csv` |

## Revision 2 limitations (real-read replay; kept for reference)

### Construct validity
1. **Containers.** Each "container" is one sequencing run of an aliquot from a single synthesis batch.
   - Coverage is strongly correlated across containers: 69% of excess variance is shared on Twist, 52% on Genscript.
   - Independently synthesised containers would be less correlated. This would change how much diversity spreading provides, and so how useful adaptive allocation is. We have not tested it.
2. **Recovery criterion.** A strand counts as recovered when its consensus matches the design exactly.
   - Our own inner code (RS + CRC) was never synthesised, so its correction gain on real reads is not measured.
   - This makes the criterion conservative for recovery, but it cannot reveal inner-code mis-corrections on real reads.
3. **Failures are injected** by deleting a container's read stream. Real vial-level failure rates are unknown, so ε and f_j are not calibrated.
4. **Costs.**
   - Read costs are abstract units and setup costs are zero.
   - Illumina flowcells are bought in fixed sizes, so batch-level adaptivity matches streaming nanopore sequencing better.
   - Savings in reads do not translate one-to-one into money.

### Internal validity
- All policies share the same read permutations.
- Baselines receive budgets tuned in hindsight, which is favourable to them.
- The φ̂ curves and recovery matrices use the evaluation split for the truth.
- **Evaluation-data leak (found by the independent audit).** "Calibrated" beliefs take yield, NB shape and dropout from the calibration split. Their per-container consensus curve φ_j, however, is the mean of the per-strand curves of even-indexed strands, and those curves are estimated from **evaluation-split** reads (`experiments/real_data.py`, `belief_from_calibration`).
  - **Affected:** the adaptive rows of the policy tables, E4 "full calibration", EC "full calibration", the ED fresh-pool belief, the E10 in-pool belief, and the choice of the single-best container.
  - **Leak-free conditions:** E4 "generic φ" and generic-prior rows (Twist: within 0.9% of full calibration), and E10 generic and Twist-derived priors (Genscript: +4–5% vs in-pool). These bound the effect.
- **Semi-empirical replay.**
  - Read orders and per-strand read counts are real.
  - Consensus success is drawn from per-strand curves estimated from real reads (two subsamples per depth).
  - Containers are emulated by subsetting runs of one pool, and only the read order varies between trials.
  - Decodability is an oracle event (exact consensus). The CRC/SHA verification path runs only in the simulated FASTQ pipeline.
- **Cache regeneration.** The real-data caches (per-strand φ̂) were regenerated between experiment runs, which shifted Genscript costs by ~9%. The policy tables were re-run on the current caches (`*_v2` result files). The superseded files are kept for transparency.
- Per-run φ̂ seeding used Python's randomised `hash()`, so re-runs are not bit-identical (see S2).

### Statistical conclusion validity
- Bootstrap CIs of a 0.99-quantile are asymmetric and can be anti-conservative when trials are few (60–100 in E4/E10/failure injection).
- The simulation ablation uses 2,000 trials, but its seed sweep showed the 0.99 target is not reliably met in s4.

### External validity
- Two Illumina pools from one study (PRJEB65931) and two nanopore datasets.
- Strand lengths are 102–108 nt.
- No long-strand, enzymatic-synthesis or cold-storage-aged pools.

### Known weaknesses of the method (measured)
| Weakness | Evidence |
|---|---|
| Posterior credible intervals are over-confident (coverage 0–83% for nominal 90%) | `paper_ed_change_detection.csv` |
| Posterior predictive P_complete is under-confident (ECE 0.30–0.60) and discriminates only late (AUC ≈0.70 at 85% of budget) | `paper_ec_*.csv` |
| Futility rule abandoned 3.7–4.5% of recoverable retrievals in simulated scenario s7 (not tested on real reads), so it is disabled | `q4_futility.csv` |
| Pivotal (Birnbaum-weighted) score gives no gain over coverage-greedy | `realdata_policies.csv`, `ablations.csv` |
| Analytical surrogate budget error −16 to +6% (Twist); Genscript: round-robin declared infeasible, single-container budget 5.0× too high | `realdata_*policies*_v2.csv` (`pred_cap@0.99`) |
| Adaptive gain is in the mean only: the q99 budget is not lower (Twist) and is 12% higher (Genscript round-robin) | `realdata_*policies*_v2.csv` (`cap@0.99`) |
| Correlated-error model not used in the controller's likelihood | code inspection |
