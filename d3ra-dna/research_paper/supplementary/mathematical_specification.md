# Supplementary S3: Mathematical Specification

This file combines the project's mathematical specification (`docs/system_model.md`) and the spreading-threshold theory (`docs/theory_spreading.md`) with the equations of the manuscript.

## Manuscript equation → code map

| Eq. | Quantity | Code |
|---|---|---|
| (1) | Read-origin probability π_ij with shared g_i (LogNormal) and container-specific h_ij (Gamma) | `channel/surrogate.py` (`sample_T`), `channel/abundance.py` |
| (2) | Recovery time T_ij = Γ(t_ij, 1)/π_ij | `channel/surrogate.py` (`sample_T`); empirical version in `realdata/empirical.py` (`empirical_T`) |
| (3) | Q_j(μ): NB-mixed, Gauss–Hermite-averaged decodability | `retrieval/adaptive_policy.py` (`_QTable`) |
| (4) | Objective: expected spend at 1−ε | `realdata/empirical.py`, `experiments/real_data.py` |
| (5) | Discrete posterior over θ_j with a failure atom | `AdaptiveController.posterior` |
| (6) | Hazard h_ij(Δ) | `AdaptiveController._hazard` |
| (7) | Pivotal score S_j; coverage-greedy score | `AdaptiveController.scores` |
| (8) | Posterior-predictive P_complete | `AdaptiveController.p_complete` |
| (9) | Hierarchical beta-binomial error propensity | `channel/correlated_errors.py` |
| Thm. | Spreading threshold q*(D, t) | `theory/spreading.py` |

## Notes on implementation details not shown in the paper
- **Q_j(μ):** also multiplies by (1 − d_syn), the synthesis dropout shared by all copies of a fragment.
- **Posterior grid:** Θ = {0} ∪ exp(linspace(log 0.01, log 5, 25)). After the first verified recovery in container j, the failure atom is set to zero.
- **Look-ahead horizon:** H = 6 batches per container.
- **Batch size:** b = 200 reads for real-data runs. If not set, it defaults to max(50, 0.02·n).
- **Overdispersion is not modelled.** The likelihood is binomial in U_j given the plug-in Q_j and ignores overdispersion. This is the over-confidence source discussed in the paper.

---

# D3RA-DNA — Mathematical System Specification

This document specifies the models that are actually implemented. File and function names refer to `src/d3ra/`.

## 1. Objects

| Symbol | Meaning | Code |
|---|---|---|
| M | number of containers (vials, compartments) | `ChannelParams.m` |
| S | number of stripes; stripe s has K_s information and N_s encoded fragments | `Layout.k`, `Layout.n_per_stripe` |
| n = Σ_s N_s | encoded fragments; fragment i belongs to stripe σ(i) | `Layout.stripe_of` |
| L | strand length in nt (address + check byte + payload + inner RS + CRC, after constrained mapping) | `StrandFormat.strand_nt` |
| X ∈ ℕ^{n×M} | placement: x_ij = copies of fragment i stored in container j | `Layout.placement` |
| y_j | usable read yield of container j | `yield_` |
| k_j | gamma shape of container-specific abundance variation (CV = k_j^{-1/2}) | `shape` |
| σ_g | log-sd of the sequence-intrinsic bias, which is shared by all replicas of a sequence | `sigma_g` |
| d_j, d_syn | per-molecule container dropout; shared synthesis dropout | `dropout`, `d_syn` |
| f_j, (D_d, q_d) | independent container failure; failure domain D_d failing as a whole with probability q_d | `FailureModel` |
| c_j, c^setup_j | cost per read; fixed cost per sequencing run | `cost`, `setup_cost` |
| φ_j(n) | probability that a cluster of n reads from container j decodes correctly | `phi` |

Each fragment s_i ∈ {A,C,G,T}^L is `codec(innerRS(addr ‖ chk(addr) ‖ payload_i ‖ CRC))`.

Codecs:
- `2bit`: 2 bits/nt.
- `scrambled-2bit`: 2 bits/nt, with a fixed keystream XOR applied first.
- `rotating-ternary`: 40/26 ≈ 1.54 bits/nt, with no homopolymers by construction.

The outer code is a systematic Cauchy Reed–Solomon (N_s, K_s) MDS code over GF(2⁸), applied column-wise across the payloads of a stripe. LT fountain codes and replication are implemented for comparison.

Logical metadata is kept in a manifest, not in the DNA: total counts, code parameters, encoding version and SHA-256 hashes. Physical overhead is therefore only the address bytes, the check byte, inner ECC and CRC (`fragmentation/addressing.py`).

## 2. Sequencing channel

### 2.1 Sampling (sequence level, `channel/sequencer.py`)

The relative abundance of fragment i in container j is

  a_ij ∝ x_ij · g_i · h_ij,  g_i ~ LogNormal(−σ_g²/2, σ_g),  h_ij ~ Gamma(k_j, 1/k_j),

normalised over the molecules stored in j. A run of R reads on container j works as follows:
1. Bin(R, y_j) reads are usable; the rest are junk.
2. Usable reads are distributed by a Poisson, multinomial or negative-binomial model over present molecules. A molecule is present unless container dropout, synthesis dropout or container failure removed it.
3. Each read is corrupted independently:
   - substitution p_sub·π(t)·Q, where π(t) is a linear 5'→3' ramp and Q ~ LogNormal is per-read quality;
   - deletion p_del·π(t)·Q·β^{run(t)−1}, which is higher inside homopolymer runs;
   - insertion p_ins·π(t)·Q, half of which are duplications of the preceding base;
   - optional burst deletions.
4. Optional extras: duplicates, which add no information, and read loss.

### 2.2 Retrieval

1. **Inventory-aware clustering.** Each read's address prefix is matched by Hamming distance (shifts 0 and ±1) against the container's inventory codebook.
2. **Two-sided bitwise-majority-alignment consensus.** A lookahead classifies each disagreement as a substitution, insertion or deletion.
3. **Inner decoding.** Inner RS decoding, then CRC and address checks.
4. **Outer decoding.** Outer MDS erasure decoding once K_s distinct fragments are verified.
5. **File check.** SHA-256 of the blob and of each file.

The outcome is classified as one of exact, partial, failed or integrity_mismatch (`decoding/verification.py`).

### 2.3 Calibrated consensus curve

φ_j(n) is measured from the sequence-level simulator: random strands, n noisy reads, BMA consensus and the inner decoder (`channel/calibration.py`). It is measured on an error-scale grid relative to the Illumina profile (sub 4·10⁻³, del 2·10⁻³, ins 4·10⁻⁴; Organick 2018) and interpolated. This is the bridge between sequence-level fidelity and fragment-level scale.

## 3. Recovery-time-matrix surrogate (`channel/surrogate.py`)

Treat the reads of container j as a Poisson process on that container's read clock. A read hits fragment i with probability

  π_ij = y_j x_ij g_i h_ij / Σ_{i'} x_{i'j} g_{i'} h_{i'j}.

Cluster (i,j) becomes decodable once it holds t_ij reads, where t_ij is drawn so that P(t_ij ≤ n) = φ_j(n). The read count at which pair (i,j) becomes recoverable is therefore

  **T_ij = Γ(t_ij, 1) / π_ij**,  and T_ij = ∞ if the molecule is absent or container j failed.

Given cumulative container reads R = (R_j), fragment i is recovered iff min_j (T_ij − R_j) ≤ 0. The file is recovered iff every stripe has at least K_s recovered fragments.

**Static allocations are closed-form per trial.** For read fractions w (Σ w_j = 1), fragment i is recovered once the total read count reaches

  τ_i = min_j T_ij / w_j.

So the minimal total reads for exact recovery is

  B*(w) = max_s K_s-th order statistic of {τ_i : σ(i)=s}.

This gives the whole distribution of budget-to-recovery from one sample of T. Reliability targets such as 0.999 are its quantiles, and 10⁴–10⁵ trials run in seconds (`static_budget_samples`).

**Adaptive policies** are simulated by advancing R batch by batch on the same T. This gives common random numbers across policies.

## 4. Analytical predictor (H1)

Given static reads R, the pair success probability conditional on the shared bias g is

  q_ij(g) = (1−d_j) · Σ_n NB(n; k_j, μ = g R_j π̄_ij) φ_j(n),  with π̄_ij = y_j x_ij / Σ_{i'} x_{i'j}.

The fragment success probability is

  p_i = (1 − d_syn) · E_g[1 − Π_{j alive} (1 − q_ij(g))].

The expectation over g uses 9-node Gauss–Hermite quadrature, which keeps replicas correlated through g. Then:
- **Stripe:** P_s = P(PoissonBinomial(p_{i∈s}) ≥ K_s), computed exactly by dynamic programming.
- **File:** P_file(R) = Σ_states w(state) · Π_s P_s(state). Failure states are enumerated exactly when 2^{M+D} ≤ 4096. Otherwise the nominal state gets its exact weight and failed states are Monte-Carlo sampled.

**Pivotal weight.** ∂P_s/∂p_i = P(S_{−i} = K_s − 1). This is the Birnbaum importance of component i in a K-out-of-N system. It is computed exactly with prefix/suffix convolutions (`pivotal_weights`). Inside the controller it is approximated by the normal ratio φ(z)/Φ(z).

## 5. Optimisation problems

### 5.1 Refined objective

The proposal's objective J = αC_seq + βC_storage + γC_latency + δP_failure double-counts. A failed retrieval also incurs the full sequencing spend, and latency is proportional to the number of batches. We therefore use the constrained form

  minimise_{X, N, π}  α E[C_seq(π; X)] + β (Σ_i Σ_j x_ij / Σ_s K_s − 1)
  subject to  P(exact recovery within cap B) ≥ 1 − ε,

where π is the retrieval policy and B its budget cap. For a static policy the expected spend equals B. For an adaptive policy it is E[min(C_stop, B)], with B the (1−ε)-quantile of the stopping cost. Latency (rounds) is reported separately rather than mixed into J.

### 5.2 Placement (Algorithm 1, `storage/placement.py`, `optimization/placement_optimizer.py`)

Each strategy satisfies Σ_j x_ij ≥ 1 and the capacity limits.

| Strategy | Code | How it places fragments |
|---|---|---|
| uniform random | `random_placement` | random |
| round-robin | `round_robin` | cyclic |
| reliability-weighted | `reliability_weighted` | proportional to a scalar score; no caps |
| failure-domain-aware | `domain_aware` | even spread over domains, then over containers |
| cost-tiered capped water-filling | `cost_tiered` | per stripe, cheapest unit retrieval cost first, with at most N_s−K_s−slack fragments per container and per domain |
| failure-exposure minimising | `failure_aware` | per-stripe counts q chosen to maximise log P_survive(q) − λ·cost (see below) |

For `cost_tiered`, the unit retrieval cost is e_j = c_j E[t_j](1+1/k_j) / (y_j(1−d_j)).

For `failure_aware`, P_survive(q) is exact. It averages over enumerated container and domain failure states the probability that the binomially-thinned (dropout) surviving count is at least K. The search starts from even spreads over the r most reliable containers for every r, then applies best-improvement single moves.

Exact MILP (`placement_optimizer.ilp`, HiGHS through `scipy.optimize.milp`):
- Variables: y_sj (fragments of stripe s in container j), u_sj (fragments read from j), z_f (state f is a loss state).
- Constraints: Σ_j y_sj = N_s; Σ_j u_sj = K_s; u ≤ y; Σ_j y_sj(1−F_fj) + K_s z_f ≥ K_s.
- Objective: Σ u_sj e_j / C_ref + (δ/ε) Σ_f w_f z_f.

The MILP gives optimality gaps for the greedy, local search, simulated annealing and genetic heuristics.

### 5.3 Static budget (Algorithm 2 baselines, `optimization/budget_optimizer.py`)

- Maximise log P_file(B v/c) over the cost simplex v by exponentiated gradient.
- Bisect on B to find the smallest budget that meets the target.

### 5.4 Adaptive sequencing (Algorithm 2, `retrieval/adaptive_policy.py`)

At each epoch, every unfinished trial gives b reads to one container.

**Belief.** A discrete posterior over a container efficiency multiplier θ_j ∈ {0} ∪ [0.01, 5]:
- The prior is a failure atom f_j plus a log-normal with sd `prior_sd`, which encodes estimation uncertainty.
- The likelihood is U_j ~ Bin(n_j, Q_j(θ R_j π̄_j)), where U_j is the number of recovered fragments and Q_j is the g-averaged decodability curve.

θ absorbs mis-estimated yield, abundance and error, and container failure, in one quantity.

**Hazard.** h_ij(Δ) = [Q(θ(R_j+Δ)π̄) − Q(θR_jπ̄)] / [1 − Q(θR_jπ̄)]. This is the probability that a still-unrecovered pair is recovered within Δ more reads. It is exact for the monotone threshold model and folds in the posterior probability that the molecule is absent.

**Controllers:**
- *uniform*: round-robin; adaptive stopping only.
- *coverage*: maximise Σ_{open pairs in j} h_ij(b) / cost_j.
- *pivotal*: maximise Σ_{pairs in j} W_σ(i) · h_ij(b) · (1−ρ_i) / cost_j. Here ρ_i is the horizon recovery probability (H = 6 batches) and W_s = ∂ log P(S_s ≥ ν_s)/∂ρ, where ν_s is the stripe's remaining need.
- *thompson*: pivotal with θ sampled from the posterior.
- *oracle*: true parameters and true failure state.
- *clairvoyant*: knows the realised T; reference only, not a bound.

**Stopping.** Stop as soon as every stripe holds K_s verified fragments. The confidence model, stated in `stopping_rule.py`:
- the undetected-error probability is at most K·2^{−8c} with a c-byte CRC per fragment, before the outer code;
- the file SHA-256 catches residual errors, and in that case sequencing continues.

Optionally, stop for futility when the normal-approximation completion probability within the remaining budget falls below ε.

### 5.5 Redundancy (Algorithm 3, `optimization/redundancy_optimizer.py`)

- Uniform parity: the smallest p meeting the target.
- Risk-adaptive: greedily add parity to the stripes with the lowest predicted log P_s until P_file ≥ target.

### 5.6 Joint optimisation (Algorithm 4, `optimization/joint_optimizer.py`)

Block-coordinate descent over (placement variant, read weights w, parity p), minimising J subject to the reliability constraint:
1. Choose the variant given (p, w).
2. Run exponentiated gradient on w given (variant, p).
3. Grid-search p given (variant, w).
4. Repeat until no block changes.

The variants are tiered with slack 0/1/2, domain-aware, and failure-aware with λ = 0 or 0.05. The comparison baseline, *independent* optimisation, chooses parity assuming uniform placement and uniform reads, places fragments domain-aware, and optimises weights last.

## 6. Modelling limitations (honest list)

- Consensus is per (fragment, container). Reads of the same fragment from different containers are not pooled. Pooling would help replication only.
- The sequence-intrinsic bias is fully shared across replicas, and the container-specific part is independent. The literature (Chen 2020; Gimpel 2023) suggests the shared fraction depends on whether replicas come from the same synthesis batch. `replication_vs_parity` sweeps this fraction.
- Failure probabilities and domain structures are hypothetical. No public data measures vial-level DNA archive failures.
- The surrogate assumes a fixed threshold t_ij per cluster. In reality consensus success is not exactly monotone in n; measured curves show a dip at n=2 from tie-breaking, which we remove with a cumulative maximum.
- Costs are in abstract units. Real sequencing cost has strong per-run fixed components, modelled only through `setup_cost`, which is 0 by default.

---

# Spreading vs concentration across failure domains: formal results

Code: `src/d3ra/theory/spreading.py`. Verification: `experiments/theory_checks.py`, which writes `results/theory_*.csv` and `results/theory_summary.json`. Prior-art status is covered in §6 and `docs/novelty_assessment.md`.

## 1. Setting

- **Stripe.** One stripe of an (N, K) erasure code. For an MDS code, any K of the N fragments recover the stripe. Parity fraction ρ = (N−K)/N.
- **Domains.** D failure domains, for example containers or groups of containers sharing a common cause. Domain d fails, losing all its fragments, independently with probability q_d.
- **Placement.** n = (n_1, …, n_D) with Σ n_d = N. Optionally, capacities n_d ≤ C_d.
- **Survival.** Without dropout, the stripe survives a failure set F iff Σ_{d∈F} n_d ≤ N − K.
- **Concentration.** All N fragments in the single most reliable domain d* = argmin q_d. Its loss probability is q_{d*}.
- **Tolerance.** τ(n) is the largest t such that *every* set of t failed domains is survived. For an even split with D | N, τ = ⌊Dρ⌋.

## 2. Homogeneous domains

**Theorem 1.** Let q_d = q ∈ (0,1) for all d, D ≥ 2, D | N, an even split n_d = N/D, and an MDS code. Let t = ⌊Dρ⌋. Then P_spread(q) = P(Bin(D,q) ≥ t+1).

**(a) Below the threshold.** If ρ < 1/D (so t = 0), then P_spread = 1 − (1−q)^D > q = P_conc for every q ∈ (0,1). Even spreading is strictly worse than concentration.

*Proof.* With t = 0 any single failure loses the stripe. 1 − (1−q)^D > q for D ≥ 2 and q ∈ (0,1), because (1−q)^D < 1−q. ∎

**(b) Above the threshold.** If 1 ≤ t ≤ D−2, there is a unique q* ∈ (0,1) such that P_spread(q) < q ⟺ q < q*. If t = D−1, then P_spread = q^D < q for all q ∈ (0,1).

*Proof.*
1. Write m = t+1 ≥ 2. Then P(Bin(D,q) ≥ m) = I_q(m, D−m+1) = ∫_0^q f(x) dx, where f is the Beta(m, D−m+1) density.
2. Let h(q) = I_q/q, the running mean of f on [0,q]. Then h′ = (f − h)/q.
3. Because m ≥ 2, f(0) = 0. f is unimodal with mode x_0 = (m−1)/(D−1). Because m ≤ D−1, f(1) = 0.
4. On (0, x_0], f is increasing, so f > h and h increases.
5. Let q_1 be the first point with f(q_1) = h(q_1). It lies after the mode. There h′ = 0 and (f−h)′ = f′ < 0, so f − h crosses zero downward. It cannot cross upward later, since that would need f′ ≥ 0 on the decreasing branch.
6. Hence h increases on (0, q_1) and decreases on (q_1, 1), with h(1) = ∫f = 1. So h(q_1) > 1 and h(0⁺) = 0.
7. Therefore h − 1 has exactly one root q* ∈ (0, q_1), with h < 1 on (0, q*) and h > 1 on (q*, 1).
8. For t = D−1, P_spread = q^D < q. ∎

**Interpretation.** ρ = 1/D is a *necessary* condition for even spreading to help, whatever the failure rate. Above it, spreading helps exactly when domains are reliable enough: q < q*(D, t).

Values of q* (`results/theory_qstar.csv`):

| D | t = 1 | t = 2 | t = 3 |
|---|---|---|---|
| 3 | 0.5000 | 1 | — |
| 4 | 0.2324 | 0.7676 | 1 |
| 6 | 0.0836 | 0.3471 | 0.6529 |
| 8 | 0.0423 | 0.1974 | 0.3955 |

For realistic container failure rates (q ≤ 0.04), any t ≥ 1 with D ≤ 8 makes spreading better. At D = 8 and t = 1 the crossover is q* = 0.042.

**Numerical check (T1).** 1,344 cases (D = 2…8, N ∈ {24, 36, 48, 60}, K from N/3 to N−1, q ∈ {0.001, …, 0.4}).
- The closed form matches exact enumeration to 4×10⁻¹⁵.
- Monte Carlo (20,000 draws) matches within 0.009.
- The sign prediction "spread better ⟺ q < q*" is correct in **all** cases with D | N.
- The 6 mismatches (0.45%) are all uneven splits (D ∤ N). The theorem does not cover them, because there the tolerance is set by the heaviest domains.

## 3. Dropout margin

With independent fragment dropout δ, losing one domain leaves N(1−1/D) fragments in expectation, of which Bin(N(1−1/D), 1−δ) survive. So the even split tolerates one loss with probability at least 1 − α only if

  N(1−1/D)(1−δ) − z_α √(N(1−1/D)δ(1−δ)) ≥ K,

that is, ρ ≳ 1/D + δ + (dispersion term).

**Check (T4)** with D = 4, K = 32, δ = 2%, q ∈ {1%, 3%, 5%}. Spreading first beats concentration at ρ = 0.273 (parity 12/44), against a plain threshold of 1/D = 0.25.

In the earlier scenario-s4 simulation, extra overdispersion and unequal domains moved the empirical crossover to parity 16/48 (ρ = 0.33). This is consistent with the margin term growing with heterogeneity. The formula is a design guide, not an exact boundary under overdispersed coverage.

## 4. Heterogeneous domains and capacities

**Proposition 2 (no tolerance ⇒ no better than concentration).** If some domain d has n_d > N − K, then P_loss(n) ≥ q_d ≥ q_{d*} = P_conc. Equality requires d = d* and no other failure set to matter.

*Proof.* The failure of d alone loses the stripe, so P_loss ≥ P(d fails) = q_d ≥ min_j q_j. ∎

- Check (T3): 1,680 random instances, 0 violations.

**Corollary 3 (general threshold).** Any placement that strictly beats concentration (generically) must tolerate every single domain failure. Single-failure tolerance is feasible iff the domains can hold N fragments with at most N−K each:

  Σ_d min(C_d, N−K) ≥ N.

Without capacity limits this means at least ⌈N/(N−K)⌉ = ⌈1/ρ⌉ domains, which is the heterogeneous form of "ρ ≥ 1/D_eff". Among tolerant placements, P_loss ≈ Σ_{pairs (a,b): n_a+n_b > N−K} q_a q_b, to first order in q. Spreading beats concentration roughly iff this is less than q_{d*}.

**Check (T5).** 200 random heterogeneous instances: D = 6, q_d ∈ [0.002, 0.1], random capacities 8–48, N = 48, random K.
- Single-failure tolerance was infeasible in 15.5% of instances, because of capacities.
- Where both a tolerant placement and the most-concentrated feasible placement exist, the best tolerant one won in **91%** of instances.
- The median ρ was 0.33 where tolerant placement won, against 0.19 where concentration won.

## 5. Dependence on the erasure code

- **Non-MDS codes** (LT/fountain needing K(1+ε) symbols) replace K by K(1+ε), so the threshold rises to roughly ρ ≥ 1/D + ε(1−ρ).
  - Check (T6), D = 4, q = 0.02: the minimum ρ where spreading wins rises from 0.27 (ε = 0) to 0.33 (ε = 0.05–0.1) and 0.38 (ε = 0.2).
- **Replication** with r copies in r distinct domains is lost when all r domains holding some fragment fail. With copies covering all domain pairs, P_loss ≈ C(D,2)q² for r = 2. This is the analogue of t = 1.

## 6. Limitations and relation to prior work

- **Independence assumption.** Domains are assumed to fail independently; correlation *within* a domain is what defines it. Correlation *between* domains (shared site, shared synthesis batch) is not covered and would raise every threshold.
- **Single stripe.** The results concern one stripe. Several stripes sharing domains fail together (the file is lost if any stripe is lost). Because all stripes use the same count vector, the file-level comparison has the same sign but different magnitudes.
- **Prior work.** These are elementary results for k-out-of-n reliability. The "spread maximally vs minimally" phenomenon is studied in distributed storage allocation theory (Leong, Dimakis and Ho, *Distributed storage allocations*, IEEE Trans. IT 2012), where the optimal allocation shifts from minimal to maximal spreading as the per-node access probability and budget grow. Our Theorem 1 and Proposition 2 should be presented as a **simple, explicit design rule for DNA containers** (with dropout margins and real-data evidence), **not** as a new theoretical result. See `docs/novelty_assessment.md` for the verified comparison.
