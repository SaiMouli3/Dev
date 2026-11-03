# Claim Verification (revision 4)

**Method.** Every number in `manuscript/main.tex` was recomputed from the committed result files by
`audit/recompute_claims.py` (revision-3 claims; output `claim_recompute_rev3.txt`) and
`experiments/statistical_analysis/r4_analysis.py` (revision-4 claims; outputs `results/r4_*.csv`).
Tables and figures are generated from those files by `make_tables.py`, `make_r4_tables.py` and
`make_dt_figures.py`, so no number in the manuscript is typed in by hand. Verdicts: **VERIFIED** matches
the evidence within rounding; **CORRECTED** the earlier revision said something different and the
manuscript now states the recomputed value.

## 1. Headline claims

| # | Claim in the manuscript | Recomputed | Verdict |
|---|---|---|---|
| 1 | Verified stopping: 3.49k reads consumed [3.37, 3.60] against 4.80–5.15k reserved, i.e. −27% to −32% | 3490 [3370, 3601]; reservations 4803 (best static) and 5153 (uniform); −27.3% and −32.2% | VERIFIED |
| 2 | With stopping on both sides, fixed proportions beat round-robin by 4.2% (pre-specified α=1) and the calibration-estimated plan does as well | α=1 (pre-specified, primary): 3343.2, −4.21% [−4.88, −3.51]; extended-grid optimum α=1.5: 3330.8, −4.57% [−5.51, −3.51]; calibration-estimated α=1: 3338.3, −4.35% [−5.06, −3.62]. All three intervals overlap, so they are one result reported three ways, not three estimates in conflict; the manuscript labels which is primary | VERIFIED (**CORRECTED**: the abstract previously quoted the extended-grid optimum) |
| 3 | Frozen configuration: every online rule costs more than round-robin, +7% to +441% | TS +7.4%, tail +12.9%, coverage +15.1%, UCB +15.3%, LTC +44.6%, pivotal +49.3%, info +440.6% | VERIFIED |
| 4 | Ablation: the penalty needs the frozen *combination* of container model, likelihood and prior width; relaxing any one substantially reduces it; floor ≈ −4.5%. Adding dropout is the physically meaningful change and helps most where loss dominates | coverage-greedy: frozen +15.1% [+10.3, +20.6]; ρ_U=0 −4.43% [−5.92, −2.91]; s=1.0 −0.56% [−2.47, +1.63]; dropout +1.61% [−0.12, +3.42]; best combination −4.70% | VERIFIED (**CORRECTED** from "the container model is the cause"; the claim-to-evidence table carried the old wording and now matches) |
| 5 | Dropout model: UCB −3.5% [−4.9, −2.2]; near-dead read share 0.27 → 0.17 | −3.5% [−4.9, −2.2]; shares 0.274 → 0.166 (round-robin 0.160) | VERIFIED |
| 6 | After an outright failure: −24% conditional, −3% over all attempts, success identical across policies | UCB 15,487 vs RR 20,309 on 102 joint successes (−23.7%); all attempts 52,786 vs 54,425 (−3.0%); success 0.34 for all 12 policies and static plans | VERIFIED |
| 7 | No adaptive rule lowers q99; hindsight proportions lower it ≈9% | 100 syntheses: RR 5150 [4853, 5202], static α=1 4700 [4550, 4753] (disjoint, −8.7%); every adaptive rule worse or overlapping | VERIFIED (**CORRECTED** from a blanket "no rule lowers q99") |
| 8 | Safety: 55,110 retrievals, 0 wrong outputs; 69.5M decode attempts, 0 false accepts; limits 6.7e-5 per retrieval and 6.3e-3 per synthesis replicate | 55,110 / 0 / 6.69e-5; 69,520,219 / 0 / 5.31e-8; 580 synthesis replicates, 95 retrievals each, per-replicate limit 6.34e-3 | VERIFIED |
| 9 | Replay: stopping −14% (Twist) and −25% (Genscript) against fixed depth | Twist 58,627 vs 68,130 (−13.9%); Genscript 636,940 vs 852,276 (−25.3%) | VERIFIED |
| 10 | Calibrated replay prior: −16% per retrieval, 1.8M calibration reads, break-even ≈17 retrievals (≈10 at 100k per container; lower levels save nothing and 30k is worse than none) | 533,179 vs 636,940 (−16.3%); 6 × 300k = 1.8M; 1.8M / 103,761 = 17.3; at 100k/container 600k / 60,054 = 10.0; at 1k/3k/10k/30k the means are 640,335 / 637,736 / 635,857 / 658,183 against round-robin's 636,940 | VERIFIED (**CORRECTED**: the claim-to-evidence table previously said "≈1 retrieval") |
| 11 | ρ_U re-selection at s=0.5: best 0.03, no value reaches 0.90, near-dead covered 0.00–0.03 | select 0.433/0.642/0.742/0.808/0.808/0.800/0.800/0.600; validate best 0.831 at ρ=0.03; near-dead 0.00–0.03 at every ρ | VERIFIED |
| 12 | Prior sensitivity: one-parameter penalty +10.9% to +20.6% (+56.4% at s=0.25); two-parameter −0.3% to +3.4% | coverage-greedy across 8 settings: +10.9% to +20.6%, s=0.25 +56.4%; with dropout −0.3% to +3.4% | VERIFIED |
| 13 | Failure modes: partial loss +76.6% → −6.0%; common cause +46.9% → −33.6% | as stated; success 1.00/0.95/1.00/1.00/0.87, identical across policies; 0 wrong outputs in 2,700 retrievals | VERIFIED |
| 14 | Six geometries: one-parameter never cheaper than round-robin (+1.1% to +93.5%); dropout model never materially worse (−27.9% to +2.4%) | as stated; static α=1 −4.0% to −10.7% everywhere; 0 wrong outputs | VERIFIED |
| 15 | Consensus cap: not binding on the baseline; fail1 success 0.313 / 0.340 / 0.357 at caps 12 / 24 / 48 | baseline means within 10 reads across caps; fail1 as stated | VERIFIED |
| 16 | Photolithographic: 0 of 660 retrievals recovered the file, 0 wrong outputs | 300 (full-design round-robin) + 240 (8 policies, reduced) + 120 (4 static) = 660, 0 successes | VERIFIED |
| 17 | Geometry arithmetic: 18,594-byte blob → 25 stripes → 1,200 strands → 72,000-read stream | `archive()` returns 1,200 fragments in 25 stripes, 200 per container, blob_length 18,594 | VERIFIED |
| 18 | Genscript: 12,000 of 12,472 designs used; 4.23% of mapped evaluation reads discarded | 12,472 designs; 12,000 selected; 4,658,997 mapped, 4,461,804 used (4.23% discarded) | VERIFIED |

## 2. Reproducibility checks

| Check | Scope | Result |
|---|---|---|
| `rerun_check.py` | 24 sampled baseline retrievals, fresh single process, reverse order | 0 of 24 differ in any recorded field |
| Replay parallelisation | Genscript conditions re-run in separate processes | identical to the sequential run |
| `clean_repro.sh` + `clean_repro_compare.py` | baseline scenario regenerated from a clean clone at the final commit | outcomes identical; 2,350 of 2,400 read counts differ because **DT4DDS 1.1.0 is not seed-deterministic**; every policy mean inside the other run's cluster-bootstrap interval in both directions; policy ordering identical |

The third check corrects a revision-3 claim that every run was bit-for-bit reproducible. It holds for the
retrieval layer given the committed read streams; it does not hold for the simulator, so the streams are
part of the artifact and the manuscript says so.

## 3. Claims deliberately **not** made

- No wet-lab evidence of any kind; the end-to-end exactness and safety results are simulated.
- No claim that the real-read replay is exact retrieval (decodability is an oracle event there).
- No claim that any policy improves the probability of successful recovery.
- No claim that adaptive allocation reduces reserved capacity.
- No regret or optimality result.
- No claim about information-greedy on the Genscript replay (censored condition).
- No claim that photolithographic DNA archives are unrecoverable in general.

## 4. Figure corrections made in the final pass

| Figure | Defect | Fix |
|---|---|---|
| 2 | x-axis label of panel (b) clipped at the page edge; overlapping log-scale ticks and a legend drawn over the data in panel (a); an empty legend entry for a scenario with no paired statistic | explicit log ticks, legends moved below the axes, `bbox_inches="tight"` so nothing is clipped, scenarios without a paired statistic dropped, caption rewritten to define every plotted quantity |
| 3 | generated at two-column width but included at one column, halving every font; a second panel showing a flat line at 1.00 | regenerated at column width as a single panel; the constant success probability is stated in the caption instead |
| 4 | **panel (b) pooled every ρ on the grid instead of using the selected ρ_U, and its title read ρ_U = 0.0** — the bars did not show the quantity the text reports | filtered to the selected ρ_U = 0.1, values annotated, title corrected; the panel now shows 1.00 for five containers and 0.20 for the near-dead one, matching Section VI-I |
| 5 | was crammed into the same float as Figure 4 at one third of the width | split into its own single-column figure with a caption that states what calibration would look like |
| 1 | the routing arrow was drawn across the text of the controller box | arrow rerouted outside the box, box widened |
