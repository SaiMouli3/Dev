# Revision-3 claim verification

**Method.** `audit/recompute_claims.py` recomputes every headline number of `manuscript/main.tex` directly from the committed result files. Its output is in `audit/claim_recompute_rev3.txt`. Each number was then compared with the text.

| Claim in main.tex | Recomputed | Verdict |
|---|---|---|
| Safety: 23,160 retrievals, 0 wrong outputs (upper 95% limit 1.6e-4) | 23,160 / 0 / 1.59e-4 | VERIFIED |
| Safety: 30.2M decode attempts, 0 false accepts (upper 95% limit 1.2e-7) | 30,164,044 / 0 / 1.22e-7 | VERIFIED |
| Baseline round-robin (RR) mean 3.49k [3.37, 3.60], q99 5.15k | 3490 [3367, 3596], 5150 | VERIFIED |
| Verified stopping 27% below the best static plan (4.80k) | 27.3% | VERIFIED |
| Baseline policy deltas vs RR (+7, +13, +15, +15, +45, +49, +441%) with CIs | match | VERIFIED |
| Near-dead container read share 16 / 23 / 25–27 / 40–41 / 66% | mean of per-retrieval shares: .16 / .23 / .25–.27 / .40–.41 / .656 | VERIFIED (definition: mean of per-retrieval shares) |
| fail1: success 0.34; UCB −24%, coverage −17%, pivotal −15%, TS −9% (n = 102 paired) | match | VERIFIED |
| fail2 and photolitho: 0 success for every policy and static plan; no wrong output | max success 0.0; 0 wrong outputs | VERIFIED |
| Homogeneous scenario: model-based rules 2–4% cheaper than RR; info-greedy ≈ RR | −2.4% to −3.6%; info −0.2% | VERIFIED |
| Sweep table values | match | VERIFIED |
| Post-hoc ρ=0: coverage −4.4%, TS/tail −3.5%, pivotal +4.4%; fail1: tail −25%, TS −20%, coverage −15%, pivotal −8% | match | VERIFIED |
| Calibration: held-out coverage 0.867 (0.84 / 0.89); 1.0 for five containers, 0.20 for the near-dead one; ρ=0 gives 0.54 | match | VERIFIED |
| Forecast: ECE 0.29 / 0.34, Brier 0.30 / 0.34, AUC 0.73 / 0.69; climatology Brier 0.25; mean forecast 0.86 vs 0.52 observed at 30% of budget | match | VERIFIED |
| Twist replay: RR 58.6k (−14% vs fixed depth, −11% vs hindsight static); rules 56.6–57.2k; Δ coverage −9.0k [−19.3, −8.5] | match | VERIFIED |
| Genscript replay: RR 637k (−25% vs fixed depth 852k; +1% vs hindsight static 631k); generic-prior rules 0–3% above RR; calibrated 533k (−16% vs RR) | UCB −0.03%, coverage +0.4%, tail +1.1%, LTC +2.2%, TS +2.7%; calibrated −16.3% | VERIFIED ("0–3%": UCB is equal to RR) |
| Re-run check: 24/24 identical | 0 mismatches | VERIFIED |

**Corrections made during this audit**
- The safety bounds are upper limits of two-sided 95% Clopper–Pearson intervals. The wording now says so.
- The abstract and findings now state the calibrated-prior exception on Genscript. Before this revision, the claim was that allocation changed reads by at most a few percent.

**Not completed:** info-greedy on the Genscript replay. The sequential run crashed after finishing (an overflow on infinite costs). The fix has been made and a rerun is in progress. The gap is disclosed as Appendix B deviation (vi).
