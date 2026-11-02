# Retrieval-in-the-Loop Placement and Adaptive Sequencing for Multi-Container DNA Archives

*Draft. Simulation study; not peer reviewed. Author list TBD.*

## Abstract

Large DNA archives will be physically split across many containers: vials, capsules or pools. These containers differ in yield, abundance skew, sequencing error, dropout, failure risk and cost. We study how to place erasure-coded strands across such containers and how to spend a sequencing budget to recover a file exactly.

The basis of the study is a *recovery-time-matrix* surrogate. Each (fragment, container) pair has a read count at which it becomes decodable. Its distribution is calibrated from a sequence-level simulator: noisy reads, bitwise-majority-alignment consensus, inner Reed–Solomon decoding with a CRC. The surrogate lets any placement and any sequencing policy be evaluated exactly over 10⁴–10⁵ Monte-Carlo trials on a laptop CPU. An independent analytical predictor matches it to a mean absolute error of 0.007 in recovery probability.

Across nine scenarios we find four things.
1. Closed-loop sequencing that stops when every stripe holds K verified fragments needs 48–79% less expected spend than the fixed depth that reaches the same 0.99 reliability.
2. Placement must be optimised *with the adaptive retrieval policy in the loop*, and the optimum is often *not* distributed.
   - When a single container already meets the failure budget, the search concentrates each stripe there, cutting expected spend by 22–37% against reliability-weighted spreading.
   - When no container does (correlated failure domains), the search spreads data over five containers and is the only design to reach 0.99.
   - A placement chosen for fixed-depth retrieval cost 12–43% more.
3. Spreading stripes over failure domains is beneficial only above a parity-fraction threshold of about 1/D. Below it, exposure-minimising concentration is 3.9× safer.
4. Bayesian container beliefs must be able to represent severe degradation. Otherwise controllers abandon degraded containers and fail under channel change, which remains the main open failure case.

We release a reproducible framework: code, tests, scenarios and statistics.

## 1. Introduction

Coding for DNA storage has matured. Inner codes correct substitutions and indels; outer erasure codes absorb strand loss. The cost of *reading* the data — the coverage-depth problem — is now an active topic [Bar-Lev et al. 2023/2025; Cao & Chen 2025; Kong et al. 2026].

Real archives will not be one homogeneous pool. Capacity limits, random access and risk management push towards many physical containers, each with its own history. Kong et al. [2026] take the first step: they analyse coverage depth when an MDS codeword is partitioned over containers that are sampled uniformly.

We study the operational setting that follows from it, in which containers are **heterogeneous**, can **fail** (alone or through shared causes), and are read by a controller that **adapts** online. Three coupled decisions arise:

1. **Placement.** Which container holds which fragments, including how many fragments of each stripe per container.
2. **Redundancy.** How much parity each stripe gets.
3. **Sequencing.** Where to send the next batch of reads, and when to stop.

**Contributions.** The ingredients are known; our claims concern their combination in this setting.

- **C-sim. Recovery-time-matrix surrogate.** It is calibrated from sequence-level consensus decoding. It turns policy evaluation into order statistics on one sampled matrix and gives common random numbers across policies (§3).
- **C-ctrl. Closed-loop sequencing controller.** It uses a discrete Bayesian posterior over container efficiency, with an explicit failure atom, and scores containers by posterior hazards weighted with stripe-level Birnbaum importance [Birnbaum 1969]. It stops on verified decodability (§4).
- **C-joint. Retrieval-in-the-loop placement.** It is a joint search over placement families (failure-exposure-minimising, cost-tiered capped water-filling, domain-spreading, reliability-weighted) and parity, scored by simulating the controller on the *believed* channel (§5). We show that the static-proxy alternative is a trap, and that the optimum is often a single container: distribute only when no container meets the failure budget alone.
- **C-rule. Spreading threshold.** This is *not* a new result: it is a direct application of Leong, Dimakis and Ho's (2012) allocation analysis. We state it as a design rule for DNA containers, with a dropout margin and a heterogeneous extension (§6.2).
- **Honest evaluation.** Nine scenarios, six baselines, Holm-corrected tests, and reported negative results (§6–7).

## 2. Related work

See `docs/literature_review.md`. The main lines are:
- DNA storage systems [Church 2012; Goldman 2013; Grass 2015; Erlich & Zielinski 2017; Organick 2018];
- channel characterisation [Heckel 2019; Chen 2020; Gimpel 2023];
- the coverage-depth line [Bar-Lev 2023/2025; Gruica 2024; Cao & Chen 2025; Bertuzzo 2025; Kong 2026];
- distributed-storage placement [CRUSH; Copysets; LRC];
- adaptive sampling [Loose 2016].

**Differences from Levy et al. [2026] and Kong et al. [2026].** Levy, Con, Yaakobi and Kiah study expected recovery time in DNA distributed storage with container failures; Kong et al. study coverage depth with uniform per-container sampling. We add container heterogeneity and failure, sequence-level calibration, adaptive allocation and stopping, and placement chosen under failure risk and retrieval cost.

**Differences from the coverage-depth line.** That work assumes a fixed sampling law. Bar-Lev et al. show that the uniform law is optimal when the sampler *cannot* be controlled; we study a controllable sampler.

## 3. Model

**Abundance.** Container j holds fragments with relative abundance a_ij ∝ x_ij g_i h_ij. Here g_i is a log-normal sequence-intrinsic bias shared by replicas, and h_ij ~ Γ(k_j, 1/k_j) is container-specific variation.

**Reads.** A fraction y_j of reads is usable. A cluster of n reads decodes with probability φ_j(n). φ_j is measured by simulating noisy reads with position-, homopolymer- and quality-dependent substitutions, insertions and deletions, followed by two-sided BMA consensus and the inner RS + CRC decoder.

**Losses.** Containers fail independently with probability f_j, and failure domains fail jointly with probability q_d. Molecules drop out with probabilities d_j and d_syn.

**Recovery-time matrix.** Treat reads as a Poisson process on each container's read clock. Then the read count at which pair (i,j) becomes decodable is

  T_ij = Γ(t_ij, 1)/π_ij,  P(t_ij ≤ n) = φ_j(n),

and T_ij = ∞ under dropout or failure. For read fractions w, the minimal total budget for exact recovery is

  B* = max_s ord_{K_s}{ min_j T_ij / w_j : i ∈ s }.

It is closed-form per trial. Adaptive policies are replayed on the same T.

**Analytical predictor.** Negative-binomial counts, Gauss–Hermite quadrature over g (which preserves correlation between replicas), exact Poisson-binomial stripe tails, and failure-state averaging.

## 4. Controller

**Belief.** For each container, a posterior over θ_j ∈ {0} ∪ [0.01, 5]. It is updated from the binomial count of fragments recovered so far.

**Score.** At each epoch the controller picks the container maximising

  Σ_{i ∈ j, open} W_{σ(i)} · h_ij(b) · (1 − ρ_i) / cost_j,

where:
- h_ij(b) is the posterior hazard that pair (i,j) is recovered within the next b reads;
- ρ_i is the probability that fragment i is recovered within a look-ahead horizon;
- W_s = ∂ log P(S_s ≥ ν_s)/∂ρ, the Birnbaum importance of stripe s normalised by its success probability, computed with a normal approximation.

**Stopping.** The controller stops when every stripe has K inner-verified fragments. The residual-error probability is at most K·2^{−16} before the outer code, and the file-level SHA-256 catches the rest. Variants: uniform round-robin (adaptive stopping only), coverage-greedy, Thompson sampling, an oracle with the true parameters, and a clairvoyant reference.

## 5. Placement and joint optimisation

Placement families:
- **Domain-spreading.**
- **Reliability-weighted.**
- **Cost-tiered capped water-filling.** Fill the cheapest containers first, with at most N−K−slack fragments of a stripe per container or domain.
- **Failure-exposure-minimising.** Choose per-stripe counts q to maximise log P_survive(q) − λ·cost, with P_survive exact over enumerated failure states including dropout. The search uses r-reliable-spread starts followed by best-improvement moves.

An exact MILP on a decomposable cost+loss objective gives optimality gaps. The heuristics come within 0–8% of the MILP optimum.

**Retrieval-in-the-loop joint search.** For each (family, parity), simulate the controller on the believed channel. Among the candidates whose believed infeasibility is at most 1 − target, pick the one with the smallest expected spend.

## 6. Results

Setup: 25 stripes of (K = 32, N = 48), 1,200 strands of 172 nt, 8 containers. Each adaptive design runs 2,000 trials; each static design runs 20,000. Static baselines receive the *true* quantile of their own budget, which favours them. Full tables are in `reports/technical_report.md` and `results/`.

### 6.1 Prediction (H1)

| Model | MAE in P(exact) |
|---|---|
| Container-specific | 0.0067 |
| Homogeneous average | 0.073 |
| Poisson | 0.36 |
| No failure model | 0.067 |

### 6.2 The spreading threshold (H2)

With 4 failure domains, parity fraction 0.27 gives an unrecoverable fraction of 14.0% for domain-spreading versus 3.6% for exposure-minimising placement. At 0.33, spreading is slightly better (3.1% vs 3.9%). At 0.38 and above, mixed placements found by the search are best.

The explanation: spreading over D domains loses about N/D fragments per domain failure. It therefore helps only if N/D plus the dropout margin is at most N − K, i.e. a parity fraction ≳ 1/D.

### 6.3 Adaptive vs fixed sequencing (H3)

The table gives the expected spend at 0.99 reliability (cost units) for the same placement under fixed and adaptive sequencing.

| Scenario | A (fixed) → C (adaptive), uniform placement | B (fixed) → D (adaptive), reliability-weighted placement |
|---|---|---|
| homogeneous | 11,281 → 5,896 (−48%) | — |
| heterogeneous | 20,755 → 7,195 | 9,893 → 5,526 (−44%) |
| non-uniform abundance | 73,974 → 15,559 | 29,659 → 10,857 |
| high error | 55,573 → 15,019 | 13,904 → 8,894 |

### 6.4 Retrieval-in-the-loop placement (H4)

Expected spend at 0.99, E (joint) vs D (reliability-weighted + adaptive):

| Scenario | E | D | Change |
|---|---|---|---|
| s2 heterogeneous | 4,330 | 5,526 | −22% |
| s3 non-uniform abundance | 7,760 | 10,857 | −29% |
| s5 high error | 5,595 | 8,894 | −37% |
| s9 uncertain estimates | 4,753 | 5,363 | −11% |

Bootstrap 95% CIs exclude 0 in all four. P(exact) at the reference budget is 0.996 for E against 0.998 for D: a consistent but non-significant reliability margin traded for cost.

**Caveat that shapes the interpretation.** In s2, s3 and s5, E puts *all* fragments in the single best container, whose failure probability of 0.5% is below ε = 1%. Its advantage is the advantage of not distributing; with only one container, the controllers tie.

Distribution was chosen where no container meets ε alone:
- **s4 (correlated failure domains):** five containers. E is the only design reaching 0.99. Removing its adaptive *allocation* (uniform round-robin) costs +79%. Placement without the correlated-failure model misses the target (0.952).
- **s9 (wrong estimates):** three containers.

The stripe-importance (Birnbaum) weighting gave no measurable benefit over plain hazard-greedy allocation (−1.5%, CIs overlap). The static-proxy joint optimiser (v1) was 12–43% *worse*, and in s4 it found no feasible design.

### 6.5 Robustness (H5) and failure cases

- **Channel change (s7).** Two containers degrade after placement. E still meets 0.99, at 31% below the oracle-tuned static budget. C and D finish within 4× the reference budget in only 96.9% and 98.0% of trials.
- **Narrow posterior.** A posterior restricted to θ ≥ 0.2 abandoned degraded containers after empty probes.
- **Compound stress (s10).** Under compound outages, overdispersion, high error and estimate noise, no design reached 0.99 at parity 20. The best recoverable fraction was 0.974.
- **Oracle.** Knowing the true parameters helped by at most 0.4% in s1–s9.
- **Prior width (H5).** With all controllers sharing the same placement in s7, an over-confident point-estimate controller completed 74% of retrievals. A Bayesian controller with prior sd 1.0 completed 97.6% (the oracle's ceiling) at 1.12× the oracle's cost. Thompson sampling completed 78%, because draws of a container's failure state caused avoidance.

H5 details and ablations: `results/h5_robustness.csv` and `results/ablations.csv`.

## 6.6 Validation on real sequencing data

The data are six Illumina runs of one 12,000-strand pool after different handling (ENA PRJEB65931), used as six real containers, plus real nanopore clusters with ground truth.

- **Shared coverage bias.** 69% of coverage variance is shared across containers (σ_g = 0.31).
- **Nanopore consensus.** On references verified against an independent flowcell, independent-error models overstate it by 6–15 points; real edit scripts replayed independently do not reproduce the gap. The strand-specific error component persists across flowcells (r = 0.57). A hierarchical beta-binomial model reduces prediction error 2.2×.
- **Real recovery-time matrices** (12,000 strands, 200 permutations of real reads):
  - Placement plus adaptive stopping is −34% versus uniform placement with fixed depth.
  - Adaptive stopping alone is −17%.
  - The surrogate is 5–18% optimistic.
- **Container loss.** Removing one container's real stream leaves round-robin placement +5–32% more costly. The 4-container tiered placement becomes +120–190% more costly and frequently stalls, and the single-container placement loses everything. The post-failure margin governs the degradation.

## 7. Negative results

- **Risk-adaptive parity.** It saved 16–32% of parity at an equal *predicted* reliability. Monte Carlo revealed that the predictor's sampled failure tail is optimistic at M = 24, so the claim is not established.
- **Shared abundance bias.** It did not specifically penalise cross-container replication. Replication loses to MDS parity for the classical need-every-fragment reason.

## 8. Limitations

- **Simulation only.** Failure rates are hypothetical, and error profiles are anchored to published Illumina measurements.
- **No cross-container pooling of reads.**
- **Abstract costs with no per-run setup cost.** Adaptive control is therefore most directly applicable to streaming (nanopore Read-Until-style) sequencing.
- **Predictor tail.** It needs exact or importance-sampled failure enumeration for large M.

## 9. Conclusion

The cheapest reliable way to read a multi-container DNA archive is to stop as soon as the data is provably decodable. For that to work, the archive must have been laid out with that adaptive reader in mind.

**Counter-intuitive consequence:** distribution is not free. At a given failure budget ε, placing data in many containers dilutes sequencing effort. The search distributes only when no single container — or, for correlated failures, no single domain — meets ε, and only as far as the parity fraction can absorb the resulting exposure (the spreading threshold).

When data is genuinely distributed, adaptive allocation matters a great deal: +79% cost without it in our correlated-failure scenario. The specific choice of allocation score matters little.

**Wet-lab validation needed.** Real container failure statistics and real per-run sequencing costs. In particular, a split-pool experiment with deliberately degraded aliquots would test both the adaptive controller and the spreading threshold.
