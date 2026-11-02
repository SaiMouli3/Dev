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
