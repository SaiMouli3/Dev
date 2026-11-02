# D3RA-DNA Technical Report

**Reliability-aware coding, fragment placement and adaptive sequencing-budget allocation for distributed DNA data storage under heterogeneous, non-uniform and noisy channels**

Status: computational study, version 0.1 (2026-11-01). Everything here is simulated, and no wet-lab claims are made. Parameter provenance is in `docs/channel_calibration.md`, the model in `docs/system_model.md`, and the positioning against prior work in `docs/literature_review.md`.

---

> **Update.** Several claims below were tested further and corrected in `reports/critical_questions.md`:
> - the spreading threshold is a known allocation-theory result;
> - the s4 "only design to reach 0.99" claim was not robust across seeds;
> - the nanopore gap is 6–15 points on verified references;
> - the shared coverage variance is batch-inherited, an upper bound.
>
> That document takes precedence.

## 0. Executive summary

We built a complete simulation pipeline. It takes files through ingestion and DNA encoding, places the strands in multiple containers, sequences them through a noisy channel, retrieves them adaptively and checks the result byte-for-byte. We then evaluated placement, redundancy and sequencing policies across 9 scenarios, with 2,000 adaptive trials and 20,000 static trials per design.

The main findings below are all *in simulation, under stated assumptions*.

1. **Adaptive stopping is the largest single lever (H3, supported).**
   - With identical placement, stopping as soon as every stripe holds K verified fragments cuts the expected spend at 0.99 reliability by 48% against the fixed depth that reaches 0.99, even in the homogeneous scenario (5,896 vs 11,281 cost units).
   - At 0.999 the cut is 69% (5,913 vs 19,010).
   - This works because the reliability tail is dominated by rare slow trials. A fixed policy must pay for them on every retrieval; an adaptive one pays only when they happen.

2. **The cost-optimal "distributed" layout is often not distributed (H4, with an important caveat).**
   - The joint optimiser scores candidate placements by simulating the adaptive controller on the *believed* channel. In s2, s3, s5 and s7 it placed **all 48 fragments of every stripe in the single best container**. That container's failure probability (0.5%) is below the 1% the 0.99 target allows, and its yield and error rate are the best available.
   - Against reliability-weighted spreading with the same adaptive controller (D), this cut expected spend at 0.99:

     | Scenario | E vs D (expected spend) | 95% CI |
     |---|---|---|
     | s2 | −21.6% | [−14%, −24%] |
     | s3 | −28.5% | — |
     | s5 | −37.1% | — |

   - **The saving comes from *not* distributing.** Spreading dilutes reads across worse containers. The price is reliability margin: P(exact) at the reference budget is 0.996 against 0.998.
   - **Design rule:** at a target of 1 − ε, if some container already meets ε by itself, distribution costs reads. Distribution earns its place only when no single container meets the target.
   - **Where distribution was chosen.**
     - *s4, correlated failure domains:* counts 9/8/15/2/0/18/0/0 across five containers. With seed 0 and 2,000 trials E reached 0.990 at 5,217. A 5-seed re-run shows it is the most reliable design in s4 in 3 of 5 seeds (0.977–0.989), but it does **not** reliably reach 0.99 (`results/q5_seed_sensitivity.csv`); the single-seed 0.990 reported earlier was within Monte-Carlo noise.
     - *s9, wrong estimates:* three containers; −11% vs D.
   - **Static proxy (negative result).** Scoring placements by their static budget-to-target instead gave designs 12–43% *worse* than this simulation-in-the-loop search. In s4 it found no feasible design at all.

3. **Failure-aware placement is non-monotone in parity.** This is an application of known allocation theory (Leong, Dimakis and Ho 2012): see `docs/theory_spreading.md` and `docs/novelty_assessment.md`.
   - With D failure domains, spreading a stripe evenly over domains beats concentrating it **only when the parity fraction (N−K)/N clears ≈ 1/D plus a dropout margin.**
   - Below the threshold, spreading makes every domain a single point of failure. At parity 12/44 with 4 domains, the file is unrecoverable 14.0% of the time with domain-spreading versus 3.6% with exposure-minimising concentration.
   - Above it (parity ≥ 16/48), spreading wins.
   - This and point 2 are the same principle at two scales: distribute only as far as the parity budget can absorb the resulting exposure.

4. **Container-specific channel models predict recovery far better (H1, supported).**
   - Mean absolute error in P(exact recovery): 0.0067 for the heterogeneous model, 0.073 for a homogeneous-average model, 0.36 for Poisson without overdispersion.
   - Caveat: the predicted 0.99 budget still undershoots by 23%. Near the tail the P(budget) curve is very flat, so small errors in probability become large errors in budget.

5. **Robustness needs wide, degradation-aware beliefs (H5, supported for post-placement channel change).**
   - In s7, two containers degrade after placement. A controller that trusts its point estimate completes only 74% of retrievals.
   - A Bayesian controller with prior sd 1.0 completes 97.6%, the same as an oracle that knows the new channel. It costs 12% more than the oracle and 3.1× less than the point-estimate controller.
   - Under moderate estimate noise (s9, s10), all controllers are equivalent.
   - **Thompson sampling is a negative result:** 78% completion in s7.
   - **Failure case:** adaptive controllers combined with naive (uniform or reliability) placement still miss 0.99 under channel change.
   - Under the compound-stress scenario s10, *no* design reaches 0.99 at parity 20.

6. **Negative or qualified results.**
   - Risk-adaptive per-stripe parity saves 16–32% of storage at a predicted reliability of 0.999. Monte Carlo shows the predictor's sampled failure tail is too optimistic for M=24, so neither design truly reaches 0.999.
   - Making abundance bias shared across replicas, rather than independent, did not make replication markedly worse. Replication loses to parity for the classic reason: every fragment must be recovered.
   - The *oracle* controller with perfect channel parameters was never meaningfully better than the Bayesian controller. In s9 it was slightly worse: +17 cost units, CI [5, 28].

---

## 1. Problem and hypotheses

A file is split into S stripes. Each stripe carries K information fragments plus N−K parity fragments of a Cauchy-RS MDS code. The fragments are distributed across M containers, which differ in yield, abundance dispersion, error rate, dropout, failure probability and cost per read. The research question is whether jointly choosing *placement*, *redundancy* and the *sequencing policy* recovers the file exactly, verified by SHA-256, at lower sequencing cost or higher reliability than conventional choices.

| Hypothesis | Statement | Outcome |
|---|---|---|
| H1 | Container-specific channel models predict retrieval better than homogeneous ones | Supported (§6.1) |
| H2 | Reliability-aware placement improves recovery under uneven coverage and container failures | Supported, with the spreading threshold (§6.2) |
| H3 | Marginal-benefit adaptive allocation reduces reads at equal reliability | Supported. Most of the gain comes from adaptive stopping; allocation adds 0–35% depending on placement (§6.3) |
| H4 | Joint optimisation beats independent optimisation | Supported with the adaptive policy in the loop: −65% vs independent (s2), and the most reliable design in s4 in 3 of 5 seeds, though below 0.99 in all of them (seed sweep). In 4/8 scenarios the optimum is a single container. Negative with a static proxy (§6.4, §7) |
| H5 | Uncertainty-aware control generalises better to unseen channels | Supported for post-placement channel change (prior sd 1.0 matches oracle completion; point estimate loses 23% of retrievals); no difference under moderate noise; Thompson sampling worse (§6.5) |

## 2. System

The full mathematical specification is `docs/system_model.md`. In short:

- **Digital layer.**
  - Ingestion: multi-file bundles, optional compression.
  - Strands: compact address, address check byte, 32-byte payload, inner RS (4 parity bytes), CRC-16, encoded with scrambled 2-bit mapping. That gives 172 nt; the rotating-ternary codec is available for homopolymer-free output.
  - Outer code: Cauchy-RS across fragments. Logical metadata lives in the manifest and is not counted as DNA overhead.
- **Sequence-level channel.**
  - Abundance: shared log-normal sequence bias times container-level gamma variation.
  - Sampling: multinomial, Poisson or negative binomial.
  - Errors: substitution, insertion and deletion, with position, homopolymer and read-quality effects, plus bursts.
  - Losses: dropout, duplicates and container failures, including common-mode failure domains.
- **Retrieval.**
  - Reads are clustered against each container's inventory, which acts as a codebook of expected address prefixes.
  - Consensus: two-sided bitwise majority alignment.
  - Decoding: inner RS + CRC, then outer erasure decoding, then SHA-256.
  - **End-to-end check.** Two real files (text and random binary) are recovered byte-exact through 5 heterogeneous containers (`tests/test_reconstruction.py`, `pipeline.store_and_retrieve`).
- **Surrogate (recovery-time matrix).**
  - T_ij = Γ(t_ij, 1)/π_ij is the read count of container j at which fragment i becomes decodable. t_ij ~ φ_j, calibrated from the sequence-level simulator.
  - From T, the exact per-trial minimal budget is closed-form for static policies, and adaptive policies are replayed with common random numbers.
  - One parity configuration with 2,000 adaptive trials takes 1–25 minutes on one CPU core.
- **Calibrated φ(n).** The probability that a cluster of n reads decodes correctly. Measured values at error scale ×1 (Illumina-like, 0.64%/nt): φ(1) = 0.55, φ(3) = 0.98, φ(4) = 1.0. At ×12 (≈7.7%/nt, nanopore-like): φ(6) = 0.50, φ(10) = 0.75, φ(20) = 0.91, φ(40) = 0.96. The ×1 values are consistent in order of magnitude with Organick et al.'s median decoding coverage of about 5 reads per sequence.

## 3. Algorithms

| Component | Methods implemented |
|---|---|
| Placement | random, round-robin, reliability-weighted, domain-aware, cost-tiered capped water-filling (slack 0/1/2), failure-exposure-minimising (λ = 0 or 0.05), plus MILP (HiGHS), greedy, local search, SA and GA on a decomposable cost+loss objective |
| Sequencing | fixed uniform, coverage-equalising static, optimised static (exponentiated gradient + bisection), adaptive uniform, coverage-greedy, pivotal (Birnbaum-importance) greedy, Thompson, oracle, clairvoyant |
| Redundancy | uniform parity, risk-adaptive greedy parity, replication, LT fountain |
| Joint | (v1) block-coordinate descent with a static-budget proxy; (v2) **adaptive-in-the-loop** search over (placement variant × parity) scored on the believed channel |

**Pivotal weights.** The adaptive controller scores container j by Σ_{open pairs in j} W_s · h_ij(b) · (1−ρ_i) / cost_j. Here:
- W_s = ∂ log P(S_s ≥ ν_s)/∂ρ, the Birnbaum importance (Birnbaum 1969) normalised by the stripe's success probability;
- h_ij(b) is the posterior hazard of recovering pair (i,j) in the next b reads;
- ρ_i is the probability that fragment i is recovered within the look-ahead horizon.

**Belief.** The posterior is over a per-container efficiency multiplier θ ∈ {0} ∪ [0.01, 5], where θ = 0 means failed. It is updated from the binomial count of recovered fragments.

## 4. Experimental protocol

- **Geometry:** 25 stripes × (K=32, N=48), i.e. 1,200 strands of 172 nt, across 8 containers. N=52 in s4 and s10.
- **Seeds:** seed 0 is the test seed. No hyper-parameter was tuned on any reported scenario; the controller's batch size, horizon, prior sd and the placement variants were fixed a priori. Seeds 100–104 are reserved for validation.
- **Static designs (A, B):** per-trial minimal budget from 20,000 trials. A static policy's spend equals its committed budget B, the 0.99 quantile.
  - **Note:** the static baselines are handed the *true* quantile of their own budget distribution. That is an optimistic, oracle-tuned budget no operator would know, so the comparison favours them.
- **Adaptive designs (C–F):** 2,000 trials with a cap of 4× A's 0.99 budget.
  - Expected spend at a target = E[min(C_stop, cap_target)], where cap_target is the target-quantile of the stopping cost.
  - P(exact) is also reported at the reference budget B_ref = A's 0.99 budget (×0.6 in s6).
- **Statistics:**
  - Bootstrap 95% CIs for quantiles and mean spend.
  - Unpaired bootstrap differences between designs with different layouts.
  - Fisher exact test (unpaired) or McNemar test (F vs E, paired on identical T) for P(exact).
  - Holm correction within each scenario; Cliff's δ as the effect size.
- **Reliability targets:** 0.99, and 0.999 where reachable. 0.9999 would need ≥ 10⁵ adaptive trials per design and was evaluated only for static designs, where the closed form makes it cheap (see `results/`).

## 5. Baseline results

Table entries: expected spend to reach 0.99 reliability (cost units = reads × cost per read), with P(exact) at the reference budget in brackets. **Bold** marks the best design that meets 0.99. "—" means 0.99 was not reached.

| Scenario | A uniform/fixed | B reliab./fixed | C uniform/adaptive | D reliab./adaptive | E joint/adaptive | F oracle |
|---|---|---|---|---|---|---|
| s1 homogeneous | 11,281 [.990] | 11,315 [.990] | **5,896** [.998] | 5,909 [.999] | 5,909 [.999] | 5,906 [1.00] |
| s2 heterogeneous | 20,755 [.990] | 9,893 [.999] | 7,195 [.998] | 5,526 [.998] | **4,330** [.996] | 4,330 [.996] |
| s3 non-uniform abundance | 73,974 [.990] | 29,659 [.998] | 15,559 [.997] | 10,857 [.998] | **7,760** [.996] | 7,760 [.996] |
| s4 correlated failures | — [.972] | — [.951] | — [.978] | — [.954] | **5,217** [.990] | 5,214 [.990] |
| s5 high errors | 55,573 [.990] | 13,904 [.994] | 15,019 [.997] | 8,894 [.996] | **5,595** [.994] | 5,595 [.994] |
| s6 tight budget (0.6×) | 20,755 [.965] | 9,893 [.997] | 7,195 [.986] | 5,526 [.998] | **4,330** [.996] | 4,330 [.996] |
| s7 channel change | 36,433 [.990] | 36,443 [.990] | — [.898] | — [.756] | **25,102** [.996] | 25,102 [.996] |
| s9 uncertain estimates | 20,755 [.990] | 12,454 [.998] | 9,023 [.986] | 5,363 [.998] | **4,753** [.996] | 4,771 [.996] |
| s10 mixed (parity 20) | — [.964] | — [.937] | — [.972] | — [.940] | — [.968] | — [.968] |

Figure: `reports/figures/baselines.png`. Pairwise comparisons with CIs: `reports/figures/baseline_comparisons.png` and the `cmp:` rows of `results/baselines.csv`.

What the table shows:

- **Homogeneous (s1).** Placement makes no difference: E chose the reliability-weighted variant, identical to D, and C ≈ D ≈ E within CI. The architecture adds **no cost** when it isn't needed, which was the stated objective of scenario 1.
- **Adaptive vs fixed with the same placement (C vs A, D vs B).**
  - C is cheaper than A wherever both reach 0.99: −48% in s1, −65% in s2, −79% in s3, −73% in s5.
  - D vs B: −44% in s2, −63% in s3, −36% in s5, −57% in s9.
  - Exception: in s9, C needs a much higher *cap* than A (28,120 vs 20,755) even though its mean spend is lower. Uniform placement plus wrong channel estimates makes the adaptive tail heavy.
- **Joint vs reliability-weighted placement with the same adaptive controller (E vs D).** E is cheaper in every heterogeneous scenario; the reductions are listed in the executive summary (§0, point 2). All differences hold under bootstrap and Cliff's |δ| ≥ 0.6.
  - **Important:** in s2, s3, s5 and s7, E's layout puts *every* fragment in container 0, the most reliable and highest-yield container, whose failure probability of 0.005 is below ε = 0.01. The comparison there is effectively "no distribution" against "reliability-weighted distribution". With one container, all controllers behave identically, which is why E = F in those rows.
  - E used several containers only in s4 (five) and s9 (three).
  - **Trade-off:** P(exact) at the reference budget is 0.996 for E against 0.998 for D. The difference is not significant after Holm correction (Fisher p ≥ 0.35), but it is consistently in D's favour. E sits closer to the 0.99 boundary by design.
- **Correlated failures (s4).** E reached 0.990 for seed 0 but is below 0.99 in all 5 seeds of the sensitivity sweep (0.977–0.989); no design reliably meets 0.99 at parity 20. The joint search picked a failure-aware placement that partly concentrates fragments in the containers and domains least likely to fail.
- **Channel change (s7) — failure case.**
  - E reaches 0.99 with 31% less spend than the oracle-tuned static budget.
  - C and D *fail* the target: only 0.969 / 0.980 of trials finish within 4× the reference budget.
  - The controller under-reacts when two containers degrade after placement. It keeps directing reads to stripes whose fragments sit in the degraded containers.
  - E's placement depended less on those containers.
  - Robust adaptive control under naive placement remains an open problem.
- **Compound stress (s10) — failure case for every design.** Outages, domains, overdispersion, high error, loss and ±30% estimate noise act together. No design reaches 0.99 at parity 20; the best fraction recoverable at all is 0.974.
  - Reliability-weighted designs (B, D) are cheapest when they succeed (about 10.0–10.4k) but the least reliable (0.94).
  - E's search, misled by the noisy beliefs, picked plain domain-spreading (0.968).
  - Meeting 0.99 under compound stress needs more parity, which the redundancy dimension of the joint search would add.
- **Oracle (F).** Within ±0.4% of E in s1–s9. In s10 it cost 2.3% more, CI [+1.8%, +2.7%]. Knowing the true parameters bought nothing measurable, which suggests the Bayesian posterior learns containers quickly enough at this scale. The clairvoyant reference (knowing realised T), measured in pilots, was about 5% cheaper than the pivotal controller.

## 6. Hypothesis results

### 6.1 H1 — channel-model accuracy

Scenario s3. Prediction of P(exact) for a static uniform allocation at 7 budgets, compared with 20,000-trial Monte Carlo on the true channel (`results/h1_model_accuracy.csv`, `reports/figures/h1_model_accuracy.png`):

| Model | MAE in P | max abs err | 0.99-budget relative error |
|---|---|---|---|
| heterogeneous (container-specific) | **0.0067** | **0.015** | −23% |
| homogeneous average | 0.073 | 0.168 | −17% |
| Poisson (no overdispersion) | 0.360 | 0.892 | −88% |
| ignores container failures | 0.067 | 0.148 | −75% |

H1 is supported for probability prediction.

The tail budget is still hard. The P(B) curve flattens near 0.99 because rare failure-dominated trials set the 0.99 quantile, so a 0.002 error in P moves B by about 20%. The homogeneous model's *smaller* budget error is a cancellation of errors: its MAE in P is 11× worse.

### 6.2 H2 — placement under failures; spreading threshold

`results/spreading_threshold.csv` (s4: 4 domains of 2 containers, domain failure 3–5%, container failure 0.25–4%). File-level unrecoverable fraction from 5,000 trials:

| parity (N−K)/N | round-robin | domain-aware | tiered (cost) | failure-aware |
|---|---|---|---|---|
| 4/36 (0.11) | 0.288 | 0.287 | 0.284 | **0.039** |
| 8/40 (0.20) | 0.164 | 0.162 | 0.181 | **0.035** |
| 12/44 (0.27) | 0.147 | 0.140 | 0.173 | **0.036** |
| 16/48 (0.33) | 0.033 | **0.031** | 0.099 | 0.039 |
| 20/52 (0.38) | 0.021 | 0.021 | 0.063 | **0.015** |
| 24/56 (0.43) | 0.010 | 0.013 | 0.066 | **0.009** |
| 32/64 (0.50) | 0.008 | 0.009 | 0.061 | **0.002** |

**Spreading threshold.** Below a parity fraction of about 1/D + margin (here about 0.3), even spreading turns *each* of the D domains into a single point of failure: losing any domain removes about N/D > N−K fragments. Concentrating in the most reliable domain is then 3.9× safer. Above the threshold, spreading wins, and the exposure-minimising search mixes the two strategies (e.g. counts 9 8 9 8 7 10 1 0 at parity 20).

The search optimises *stripe*-level survival, so at parity 16 it is slightly worse at file level than domain spreading (0.039 vs 0.031). Correlated dropout across 25 stripes is not part of its objective; this is a known limitation.

The cost-tiered placement, which concentrates on cheap containers up to the N−K cap, is **never** best for reliability here. Cheap containers were also the failure-prone ones in this scenario. It is useful only when combined with failure awareness (the λ > 0 variant, "failure1"), which the joint search selected in 8 of 9 scenarios.

Both parts of the threshold rule follow from erasure-coding first principles. The rule is a direct special case of the symmetric-allocation analysis of Leong, Dimakis and Ho (IEEE Trans. IT 2012), comparing spreading over m = D nodes with m = 1. Nath et al. (NSDI 2006) observed the same phenomenon empirically. We present it as an explicit design rule for DNA containers, with a dropout margin and a heterogeneous extension (Theorem 1 and Proposition 2 in `docs/theory_spreading.md`), not as a new result.

### 6.3 H3 — adaptive sequencing

Decomposition on the full design (`results/ablations.csv`, s10 and s4): the full controller vs "adaptive stopping with uniform round-robin allocation" vs "coverage-greedy" vs "pivotal". See §7. In the baseline table, the C-vs-A and D-vs-B gaps measure adaptive stopping plus allocation together.

An early pilot was run before the final configuration: 6 random heterogeneous containers, 20×(30,40) stripes, 2,000 trials, with the earlier θ ≥ 0.2 posterior. The pilot script is not part of `results/`, so treat these numbers as indicative:
- uniform-allocation adaptive stopping: mean 4,524, q99 11,543;
- coverage-greedy: mean 4,161, q99 8,077;
- pivotal: mean 4,154, q99 7,792;
- Thompson: mean 4,150, q99 7,639.

Compared with the static uniform allocation (q99 11,498 on the same T), stopping gave most of the *mean* saving and allocation mostly trimmed the *tail*: roughly −30% at q99, because allocation steers away from failed containers. Pivotal and coverage-greedy differed by only about 3%. The reproducible version of this decomposition is in the ablations (§7).

### 6.4 H4 — joint vs independent optimisation

See §5 (E vs D) and §7 (ablation G and the static-proxy joint). The v1 negative result is kept in `results/archive_v1_static_proxy_joint/`: in s2, the static-proxy joint design cost 6,202 vs D's 5,526 (+12%). The v2 adaptive-in-the-loop joint design cost 4,330 (−22%).

**Lesson:** joint optimisation is only as good as the retrieval model inside it. A placement tuned for fixed-depth retrieval is the wrong placement for adaptive retrieval.

### 6.5 H5 — robustness to channel misspecification

`results/h5_robustness.csv`. Every controller uses the same tiered placement (slack 1), so the infeasible fraction caps reliability for all of them; we therefore compare *completion within the cap* (4× static 0.99 budget) and the conditional cost.

| Scenario | Controller | P(complete) | Mean cost when complete |
|---|---|---|---|
| s7 channel change (ceiling 0.976) | point estimate (prior sd 0.02) | **0.742** | 39,002 |
| | Bayes, prior sd 0.5 | 0.974 | 22,483 |
| | Thompson, prior sd 0.5 | 0.780 | 16,625 |
| | **Bayes, prior sd 1.0** | **0.976** | **12,428** |
| | oracle parameters | 0.976 | 11,072 |
| s9 estimate noise (ceiling 0.964) | all five controllers | 0.962–0.964 | 5,033–5,150 |
| s10 mixed (ceiling 0.920) | all five controllers | 0.920 | 9,885–10,790 |

What the table shows:

- **H5 is supported when the channel shifts after placement (s7).**
  - The controller that trusts its point estimate loses 23% of retrievals and costs 3.1× more when it does succeed.
  - Widening the prior to sd 1.0 recovers the oracle's completion rate (0.976) at only 12% higher cost than the oracle.
- **Moderate estimate noise (s9, s10).** All controllers are equivalent. The posterior corrects ±30–50% errors within the first probes.
- **Negative result: Thompson sampling.** It is markedly *worse* than the posterior mean in s7 (0.780). Sampling θ from a posterior with a failure atom occasionally draws "failed" for a degraded but live container, and the resulting avoidance compounds.
- **Static allocation optimised on the wrong belief: not applicable.** The believed channel already declared 0.99 unreachable for this placement, so no budget could be optimised. A static plan cannot recover from that.

**Development finding.** With posterior support θ ≥ 0.2, the controller abandoned heavily degraded containers after one empty probe in s7: 4% of trials did not finish (archived in `results/archive_v2_narrow_theta/`). Support down to 0.01 cut this to 0.7% in a 300-trial check. In the 2,000-trial baseline run, C and D in s7 still finish within the cap in only 96.9–98.0% of trials, so the robustness problem is reduced, not solved.

## 7. Ablations

Source: `results/ablations.csv`. Each ablation removes one component from the full design E (adaptive-in-the-loop joint placement + pivotal Bayesian controller), with 2,000 adaptive or 10,000 static trials and the same 4× cap. Values are expected spend at 0.99, or "—" with the achieved success rate where 0.99 is not reached.

**s2 (heterogeneous). E = all fragments in container 0.**

| Ablation | Expected spend @0.99 | vs full | P(complete) |
|---|---|---|---|
| full (E) | 4,347 [4,328; 4,364] | — | 0.990 |
| A: uniform (round-robin) placement | 7,232 | +66% | 0.996 |
| B: static optimised weights (fixed depth) | 5,739 (committed) | +32% | 0.994 |
| E: fixed stopping at the controller's mean allocation | 5,739 (committed) | +32% | 0.994 |
| F: oracle channel knowledge | 4,347 | 0% | 0.990 |
| G: independent instead of joint optimisation | 7,176 | +65% | 0.998 |
| static-budget-proxy joint (v1) | 6,195 | +43% | 0.999 |
| uniform allocation / coverage-greedy controller | 4,347 | 0% (one container, nothing to allocate) | 0.990 |

**s4 (correlated failure domains). E = 5 containers.**

| Ablation | Expected spend @0.99 | vs full | P(complete) |
|---|---|---|---|
| full (E) | 5,202 [5,106; 5,309] | — | 0.990 |
| A: uniform placement | — | target missed | 0.984 |
| B: static optimised weights (fixed depth) | — | target missed | 0.987 |
| D: placement optimised *without* correlated-failure modelling | — | target missed | 0.952 |
| E: fixed stopping at the controller's mean allocation | — | target missed | 0.987 |
| F: oracle channel knowledge | 5,184 | −0.3% (CI overlaps) | 0.990 |
| G: independent instead of joint optimisation | — | target missed | 0.982 |
| static-budget-proxy joint (v1) | no feasible design found | — | — |
| adaptive stopping only (uniform round-robin allocation) | 9,323 | **+79%** | 0.990 |
| coverage-greedy instead of pivotal | 5,124 | **−1.5%** | 0.990 |

What the ablations show:

1. **Placement matters most.**
   - In s2, removing it costs +66%.
   - In s4, every ablation that changes placement misses 0.99: uniform, independent, and no correlated-failure model. Ignoring the domains is the worst of these (0.952).
2. **Adaptive stopping vs fixed depth.**
   - Committing a fixed budget instead of stopping adaptively costs +32% in s2.
   - In s4 a fixed budget cannot reach 0.99 at any affordable depth (0.987).
3. **Adaptive allocation matters only when there is something to allocate.**
   - In s4, where data spans 5 containers, the uniform round-robin controller costs +79% over the pivotal controller. Allocation steers reads away from containers whose domain has failed and towards the cheapest live holders of the needed fragments.
   - **Negative result for the pivotal weighting.** Coverage-greedy (no Birnbaum weights) was 1.5% *cheaper* than pivotal, with overlapping CIs. The stripe-importance weighting brought no measurable benefit here. The gain comes from hazard-based allocation itself.
4. **Perfect channel knowledge** (oracle) changed costs by ≤ 0.3%. The posterior learns containers quickly at this scale.

## 8. Other experiments

**Placement optimiser quality** (`results/placement_opt.csv`; 4–6 stripes, so that the exact MILP solves):

| Scenario | random | round-robin | greedy / LS / SA / GA | MILP | MILP time |
|---|---|---|---|---|---|
| s2 | +96% | +24% | +1.4% | 0 | 0.4 s |
| s4 | +100% | +56% | 0% | 0 | 3.3 s |
| s10 | +307% | +31% | +8.0% | 0 | 11.3 s |

The local-search, SA and GA heuristics never improved on the capped greedy start. It is a strong local optimum for this objective, but not always globally optimal (8% gap in s10). The MILP is practical up to about 6 stripes × 8 containers × 4,096 failure states. Larger instances should use greedy plus MILP on stripe templates, since all stripes share one count vector.

**Redundancy** (`results/redundancy.csv`). 24 containers; each of 40 stripes lives in one group of 4 containers.
- Risk-adaptive parity used 16% less parity at 6 reads per fragment and 32% less at 10, at an equal *predicted* reliability of 0.999.
- Monte Carlo reliability: 0.987 vs 0.984 at 6 reads per fragment, and 0.991 vs 0.983 at 10. Both miss 0.999, and the adaptive allocation loses more.
- **Cause:** with M = 24 the predictor samples 64 failure states, which under-represents multi-failure states in small groups.
- **Verdict on C4: not established.** The direction is promising, but it needs an exact or importance-sampled failure predictor.

**Replication vs parity under shared bias** (`results/replication_vs_parity.csv`). At equal stored molecules (2×) with total log-variance fixed at 0.36:

| shared fraction | 2× replication P(recoverable) | RS parity P(recoverable) | RS 0.99 budget |
|---|---|---|---|
| 0.00 | 0.537 | 1.000 | 10,742 |
| 0.50 | 0.541 | 1.000 | 10,200 |
| 1.00 | 0.539 | 1.000 | 10,152 |

The hypothesis that a shared sequence-intrinsic bias would particularly hurt replication was **not** supported. Replication fails because it needs *every* fragment: the probability that both copies of some fragment are missing is about 46% over 640 fragments. Parity has no such requirement, as the coverage-depth literature also shows. The shared fraction changed RS budgets by only about 5%.

**Scalability** (`results/scalability.csv`, one CPU core):

| containers | fragments | placement | predictor + bisection | static MC (2,000 trials) | adaptive (200 trials) | peak memory |
|---|---|---|---|---|---|---|
| 10 | 880 | 0.02 s | 2.1 s | 0.5 s | 2.5 s | 0.17 GB |
| 50 | 4,400 | 0.15 s | 0.2 s | 4.7 s | 46 s | 0.86 GB |
| 100 | 8,800 | 0.24 s | 2.5 s | 8.2 s | 110 s | 1.7 GB |
| 100 | 26,400 | 1.3 s | 4.9 s | 28.7 s | 473 s | 5.2 GB |

Cost scales roughly linearly with fragments × trials. Memory is dominated by the dense T matrix (trials × pairs) and is the binding constraint. Chunking over trials is implemented for static evaluation and would be needed for adaptive evaluation beyond about 30k fragments.

At 50–100 containers the predictor's sampled failure tail disagrees with Monte Carlo on whether 0.99 is reachable at all; the predicted 0.99 budget was finite while Monte Carlo found it unreachable. This is the same limitation seen in the redundancy experiment.

## 9. Threats to validity and limitations

- **Simulation only.** No synthesis, storage or sequencing was performed. Container failure rates and domain structures are hypothetical; no public dataset measures vial-level failures. Error profiles are anchored to published Illumina values, and the ×12 grid point is a nanopore-like assumption.
- **Static baselines get the true quantile**, which is optimistic for them. The adaptive designs work from beliefs, which are wrong in s7, s9 and s10.
- **Pooling.** Consensus is per (fragment, container); reads from different containers are not pooled. This understates replication and overstates the value of placement.
- **Failure predictor.** Sampled failure states are inaccurate for large M. Exact enumeration is used only up to 2^12 states.
- **Reliability targets.** 0.999 was not reachable by any design in most scenarios, because the hypothetical failure rates (0.5–8% per container) leave an infeasibility floor at parity 16. Targets of 0.999 and above need more parity or lower failure rates. Their sequencing comparison remains future work.
- **Abstract costs.** Real sequencing has large per-run fixed costs (`setup_cost`, 0 by default). Positive setup costs would penalise the many small adaptive batches, which average 100–750 rounds. On a real nanopore with Read Until, batches are effectively free; on Illumina they are not. The adaptive gains therefore apply most directly to nanopore-style streaming sequencing.
- **Novelty.** See `docs/literature_review.md` §9. The closest prior work is Levy, Con, Yaakobi and Kiah, *Expected recovery time in DNA-based distributed storage systems* (arXiv:2602.07601, Feb 2026), which covers containers with failures; and Kong, Elishco, Wang and Duman, "The coverage depth problem in distributed DNA data storage" (arXiv:2610.02931, Oct 2026), which covers uniform, failure-free, non-adaptive sampling. The pivotal weight is Birnbaum importance (1969). Failure-domain caps are standard in erasure-coded storage (CRUSH, Copysets). Our claims are therefore about the **combination**: adaptive closed-loop sequencing, failure-aware placement selected with the adaptive policy in the loop, and the spreading-threshold rule, studied in multi-container DNA archives.


## 11. Real-data validation (summary)

Full details are in `reports/real_data_validation.md`. Data: ENA PRJEB65931 (6 real Illumina runs of one 12,000-strand pool) and Zenodo 14296588 (real nanopore clusters).

- **Retrieval bug found by real data, now fixed.** Two-sided consensus fails on raw reads with read-through tails. Forward-only consensus raises 3-read exact recovery from 1% to 97%.
- **Real channel.**
  - Substitution-dominated errors, 0.35–1.15%/nt. Deletions 2–3× and insertions about 4× lower than my default.
  - Real container heterogeneity: up to 2× the reads per strand and 4× the dropout.
  - **69% of coverage variance is shared across containers** (σ_g = 0.31, matching the assumed 0.30).
  - The real Illumina φ(n) matches the simulator.
- **Nanopore.** On references verified against an independent flowcell, independent-error simulation overstates consensus success by 6–15 points (20–40 on the Microsoft dataset, which has a reference caveat). Replaying real edit scripts independently gives ≥ 99.6% at 10 reads, so reads of a strand are dependent. The strand-specific error component persists across flowcells (r = 0.57). A hierarchical model cuts prediction error 2.2×.
- **Placement × policy on real recovery-time matrices.**
  - Single best container with adaptive stopping is −34% vs uniform placement with fixed depth.
  - Adaptive stopping alone is worth 15–17%, much less than in simulation, because these real containers never fail.
  - Adaptive allocation is worth −10% when data is spread over 4 containers; pivotal and coverage-greedy scoring tie.
  - The surrogate is 5–18% optimistic.
- **Container loss on real streams.** The single-container layout loses everything if its container fails. Tiered (4 containers, 2 spare fragments per stripe after a loss) costs +120–190% and often stalls. Round-robin (6 spare) costs +5–32% and always finishes. **The post-failure margin governs degradation.**
- **End-to-end.** Real files (175 KB) were recovered byte-exact through FASTQ with real Illumina and nanopore error scripts transplanted onto our strands.

## 10. Reproducibility

```bash
pip install -r requirements.txt
```

```bash
python -m pytest -q
```

```bash
python experiments/run_all.py
```

```bash
python experiments/merge_results.py
```

```bash
streamlit run dashboard/app.py
```

- `python -m pytest -q` runs 199 tests in about 20 s.
- Real data: `python experiments/real_data.py download prep transplant nanopore nanodiag channel policies failures plots` (about 1 h).
- `python experiments/run_all.py` takes about 1–2 h on 4 cores when run with `--tag` parts in parallel.
- `python experiments/merge_results.py` writes `results/baselines.csv`, `reports/figures/*.png` and `reports/results_tables.md`.
- All random number generators are seeded. Archived intermediate runs are kept in `results/archive_*` with the reason they were superseded.
