# Missing Evidence and Pending Work (revision 4)

What the manuscript does **not** establish, and what would be needed to establish it. Ordered by
importance to the paper's claims. The revision-3 version of this file is kept below for comparison.

## Still open

| # | Gap | Why it matters | What would close it | Possible without a laboratory? |
|---|---|---|---|---|
| 1 | **No physically independent containers.** Every real-read "container" is a run of one synthesis batch; every independently synthesised container is simulated by DT4DDS. | The size of every allocation effect depends on how different containers really are, and on whether real containers show the dropout-dominated behaviour that drives our results. | Synthesise one design in ≥4 independent batches, age and amplify them differently, sequence each separately, and measure each container's fraction of permanently absent strands. | No |
| 2 | **No end-to-end wet-lab retrieval.** Our codes were never synthesised; on real reads, decodability is an oracle event (consensus equals the designed strand). | The inner-code correction gain and the real false-accept rate are unmeasured, so the safety result is simulated. | Encode a known file with these codes, synthesise, split into containers, sequence, feed the FASTQ into this decoder, compare bytes. | No |
| 3 | **No observed container-failure statistics.** All failure modes (total loss, partial loss, low yield, contamination, common cause) are injected at rates we chose. | The failure prior and the redundancy a real archive needs cannot be calibrated from our work. | Long-term storage studies reporting vial or capsule loss rates and partial-degradation profiles. | No |
| 4 | **No regret or optimality result.** Every claim is empirical, about the configurations tested. | A theory-oriented venue will want a bound; our negative result is also not proved to be unavoidable. | A bandit formulation with K-of-N stripe completion, or a lower bound for the dropout channel. | Yes (theory) |
| 5 | **Nanopore correlation evidence is two flowcells.** Supporting only. | It motivates over-dispersed models but does not validate our container model. | More independent flowcells or datasets. | Partly (depends on public data) |
| 6 | **Information-greedy on the Genscript replay is censored.** The condition did not complete; >1% of its trials never recovered the file within the available reads. | One cell of one secondary table is absent. | Re-run it with a longer budget; the fix for the overflow that lost the first run is in the artifact. | Yes (compute) |
| 7 | **Cost model uses assumed prices.** Three charging models, not quoted vendor prices; instrument time, minimum run sizes and reagent costs are not modelled. | Monetary conclusions are illustrative. | Real quotations for a specific instrument and service. | Partly |
| 8 | **Target venue not chosen.** The paper is 25 pages with appendices. | Page limits and formatting are venue-specific. | Choose a venue and cut to its limit. | Yes |

## Closed since revision 3

| # | Former gap | How it was closed in revision 4 |
|---|---|---|
| R3-4 | Posterior over-confidence; ICC chosen post hoc | ρ_U re-selected with the prior actually used (s = 0.5) on an extended grid {0 … 0.5}; reported per container. No value reaches 0.90 coverage and the near-dead container is covered 0.00–0.03 of the time at every value, so the failure is now characterised rather than hidden (`results/dt_calibration_*_r4sd05.csv`). |
| R3-6 | Currency cost / run economics | Cost model with per-batch and per-container charges under three charging models; container switches recorded (`results/r4_cost_model.csv`). Vendor prices still assumed. |
| R3-8 | Ablation of the allocation components | 2×2×2 ablation over {container model} × {ρ_U} × {prior width} on baseline and fail1, plus an 8-setting prior sweep (`results/r4_ablation.csv`, `results/r4_prior_sensitivity.csv`). |
| R3-9 | Bit-identical reproducibility | All seeds derive from `zlib.crc32`; 24 sampled retrievals reproduce every field, the replay conditions reproduce across process layouts, and a scenario-level clean reproduction from a fresh clone is in `evidence/verification_reports/`. |
| R3-13 | Evaluation-split leak in calibrated beliefs | Leak-free split documented and audited: calibration is the first 300k reads of each run, evaluation the next 900k; the generic prior sees neither. The leaky analysis is retained only as an upper bound in Appendix A. |
| R3-14 | No tail reduction | Measured properly with cluster-bootstrap intervals and 100 synthesis replicates: no adaptive rule lowers q99, the hindsight static plan lowers it by about 9%. Reported as a negative result rather than a limitation to be fixed later. |
| R3-15 | Static weights with verified stopping never evaluated | Now the central baseline: the whole static family runs with verified stopping, with a finer α grid and a calibration-estimated variant (`results/r4_actual_vs_reserved.csv`). |
| R3-1/2/3 | Single container model; failures only as total loss; one geometry | A two-parameter (efficiency + dropout) container model; four further failure modes; six further archive geometries. The underlying wet-lab gaps (items 1–3 above) remain open. |

---

