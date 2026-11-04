# Change log — final revisions (revisions 6 and 7)

Scope: only the corrections requested, plus defects found by the verification they required. No
experiment was re-run, no result file was altered, and no finding was added or removed. Section
letters refer to the rebuilt PDF (25 pages); page numbers are from that build.

## 1. Corrections made

| # | Item | Where (page / section) | What changed | Why | Numbers changed? |
|---|---|---|---|---|---|
| 1 | q99 claim | p1 Abstract; p2 Findings (iii); p11 §VI-H; p18 Discussion; p21 Table XVIII | "No adaptive rule lowered the 0.99-quantile" → no adaptive rule lowered it **reliably**. §VI-H now separates numerical from statistically supported reductions and lists all four numerical reductions with their intervals, all of which contain zero, against the one reduction that excludes zero (fixed proportions) | The dropout-aware runs do show numerical q99 reductions (UCB −8.7%, coverage-greedy −1.0% in the frozen design; coverage-greedy −9.0%, Thompson −7.8% at 100 syntheses). The old wording denied them | No. Values were already in `r4_tails.csv`; they are now stated rather than summarised away |
| 2 | 4.2% vs 3.5% | p1 Abstract; p2 Findings (ii); p8 §VI-B; p18 Discussion; p19 recommendation box; p20 Conclusion | Both figures are now labelled by experiment: 4.2% is the **pre-specified fixed allocation** (no controller, baseline, both sides stopping); 3.5% is **UCB under the dropout-aware model** (an online rule, secondary non-frozen configuration). Text states they are alternatives, not additive | They came from different configurations and could be read as one result or as cumulative | No |
| 3 | Calibration protocol | p13 §VI-I | Added an explicit statement that **every policy run, primary and revision-4, uses the frozen ρ_U = 0.1**, that the s = 0.5 re-selection is a sensitivity experiment which changes no policy run, and that **no primary conclusion changes**, with the reason | The two calibration exercises could be read as one superseding the other | No |
| 4 | Real-read validation | p15 §VI-J (new paragraph "What is real here and what is not") | Itemises what comes from the sequencing data (reads, counts, strand identity, yield, unmapped fraction, arrival order) and what we supply (container partition, the coded layout, and the **oracle** recovery decision). States that the inner RS, CRC-16, address check, outer code, CRC-32 and SHA-256 **are not executed** on this path and no bytes are compared, so the replay cannot support any statement about wrong outputs or byte-exact retrieval | Item 4 | No |
| 5 | Safety claims | p8 §VI-A | Added that zero observed errors is **not** evidence of zero error probability, only of a rate below the stated bounds, and that the bounds apply to this archive configuration only | Item 5 | No |
| 6 | Reads vs money | p19 §VI-L | Added an explicit three-way separation of reads consumed, capacity reserved and laboratory expenditure, and a statement that 27–32% is **not** a demonstrated reduction in sequencing cost | Item 6 | No |
| 7 | Title | p1 | "Verified Stopping, **Not Adaptive Allocation**, Saves Reads …" → "Verified Stopping **Dominates** Adaptive Allocation …", subtitle unchanged | The old title implies allocation never helps, contradicted by our own results: −6% under partial loss, −34% under common-cause degradation, −35% after a container failure. "Dominates" states the size relation (27–32% vs ≤4.5% in the baseline) without the false absolute | No |
| 8a | **Container switches (error found by verification)** | p19 §VI-L | "adaptive rules switch container 54–61 times per retrieval and round-robin 63" → **38–61** for adaptive rules (tail-aware fewest at 38, coverage-greedy most at 61), round-robin 63, and the fixed-proportion plans **most** of all (65 and 69) | The stated range omitted tail-aware (37.7) and the claim "no policy is switch-frugal" ignored that static plans switch most | **Yes — a wrong range in the paper is corrected.** `r4_cost_model.csv` unchanged |
| 8b | Verification scaffolding | new file `audit/verify_manuscript.py`, output `audit/verify_manuscript.txt` | 86 numeric checks plus 10 PDF-text checks recomputing every headline figure from the committed result files | Item 8 | No |

## 2. Verification results

- **Numeric (item 8).** `audit/verify_manuscript.py` → **all 96 checks pass**. Coverage: abstract, §VI-A
  safety, §VI-B consumed vs reserved, §VI-C policy comparison and read shares, §VI-D failed container,
  §VI-E ablation, §VI-F dropout model, §VI-H tails, §VI-I calibration and forecast, §VI-J replay and
  calibration cost, §VI-K consensus cap / failure modes / geometries, §VI-L cost model, and the
  photolithographic totals. Each check names the section, the recomputed value, the value in the paper
  and the tolerance.
- **References (item 9).** All 54 entries checked. 43 carry DOIs that resolve with matching title and
  year (Crossref; the two Zenodo datasets resolve on DataCite). The 11 without DOIs are preprints
  (Kong, Bertuzzo), venues that issue none (JMLR, OSDI, USENIX ATC, NeurIPS, COLT) and the ENA study
  accession. Kong et al. 2026, Levy et al. 2026 and Cao et al. 2026 were read in full and their
  described contributions match the text.
- **Formatting (item 10).** Official IEEEtran v1.8b, `conference` option, two-column. All fonts
  embedded and subsetted. No raster images: every figure is vector, so resolution is not at issue. One
  overfull hbox remains, 2.77 pt in the notation table, below the visible threshold. No float
  overflow, no unresolved references or citations, no placeholder text. Page numbers are absent, which
  is the IEEEtran conference default.
- **Reproducibility (item 11).** Repository public and reachable (HTTP 200). `final/COMMIT.txt` names
  the commit the PDF was built from and explains that the bundle commit follows it. Every script,
  protocol file, verification report and result file named in the manuscript exists, and every
  documented sub-command is present in the code. DT4DDS is pinned at 1.1.0. The manuscript states that
  DT4DDS is **not** seed-deterministic and that only aggregate conclusions, not individual retrievals,
  reproduce on re-simulation.

## 3. What did not change

The research question, the experimental design, every result file, every table of measurements, the
figures' underlying data, the frozen protocol and its deviation log, the reference list (one addition
in the previous revision, none here), and all twelve limitation disclosures of item 12 — no wet-lab
work, computational digital twin, oracle decodability in the replay, simulated container failures, no
demonstrated physical instrument switching, no regret or optimality result, uncalibrated completion
forecasts, and results that apply to the tested configurations only. The one substantive numerical
correction is the container-switch range (8a), which was wrong in the manuscript and is now right.


---

# Revision 7 — the posterior-calibration concern

**Concern raised.** The beta-binomial correlation $\rho_U$ was selected with prior width $s=1$ while the
policy runs use $s=0.5$; the corrected selection picks $\rho_U=0.03$, and no value on the grid reaches
90% coverage for the near-dead container. A reviewer could reasonably ask why the primary
adaptive-allocation comparison relies on the original configuration.

**What was done.** Rather than argue the difference is small, the **complete primary comparison was
re-run at the re-selected $\rho_U=0.03$**: all eight policies, all three primary scenarios
(baseline, homogeneous, one failed container), the full 30 syntheses x 10 read orders, with seeds,
prior width, container model, placement and decoder held at their frozen settings. New result files:
`results/dt_{baseline,homog,fail1}_rho003.csv` (7,200 retrievals). New analysis
`r4_analysis.likelihood_sweep()` -> `results/r6_likelihood_sweep.csv`; new Table XI.

**Internal check.** Round-robin never evaluates the posterior, so its retrievals must be bit-identical
across likelihood settings. Verified programmatically and asserted in the analysis; it is the common
reference for all three columns.

**Result — this changed a headline sub-claim.**

| Baseline | $\rho_U=0.1$ (frozen) | $\rho_U=0.03$ (re-selected) |
|---|---|---|
| Coverage-greedy | +15.1% [10.4, 21.0] | **−0.4% [−2.5, +1.7]** |
| Thompson sampling | +7.4% [4.3, 10.9] | **+0.3% [−1.7, +2.4]** |
| UCB | +15.3% [10.8, 20.4] | **+1.5% [−0.5, +3.8]** |
| Tail-aware | +12.9% [8.6, 17.6] | **−1.2% [−3.2, +0.6]** |
| Pivotal | +49.3% | +12.2% [8.8, 15.8] |
| Learn-then-commit | +44.6% | +34.9% [29.7, 39.6] |
| Info-greedy | +441% | +576% |

Four of seven rules move from clearly costlier than round-robin to **statistically indistinguishable
from it**; three remain clearly worse. The homogeneous control is unaffected (−2.8% to −3.7% under
both settings), as it should be with no near-dead container to misread. After an outright failure the
effect is mixed: coverage-greedy keeps less of its gain (−8.7% vs −17.0%) and UCB loses most of it
(−2.3%, n.s., vs −23.7%), while tail-aware (+1.6% -> −16.7%) and learn-then-commit (+7.9% -> −17.8%)
gain. No single $\rho_U$ is best everywhere.

**Claims rewritten** (Abstract p1; Findings (ii) p2; new Section VI-J p16 with Table XI; Discussion
"Calibration versus cost" p20; Conclusion p22; Table XVII limitations; Table XX claim-to-evidence):

- *Strong form, now scoped*: "every model-based rule read 7–441% more than round-robin" holds for the
  **pre-registered** configuration, not for the corrected one.
- *Weak form, now the load-bearing claim*: **no adaptive rule beat round-robin by more than a few per
  cent in any (scenario, likelihood) combination tested**, and three of seven were clearly worse in most.
- *Main thesis unaffected and arguably strengthened*: at the corrected likelihood, adaptive allocation
  on the baseline buys essentially nothing against stopping's 27–32%.
- *Why $\rho_U=0.1$ remains primary*: it was pre-registered, and re-tuning after seeing outcomes would
  unfreeze the protocol. Every statement about the **size** of the penalty now names its likelihood.
- *Calibration and cost separated*: no $\rho_U$ calibrates the near-dead container (a container-model
  failure, addressed only by the dropout parameter), yet $\rho_U$ strongly affects cost. A configuration
  can be badly calibrated and still allocate acceptably; our strongest negative result sat at a setting
  poor on both counts.

**Also corrected.** The clean-reproduction re-runs were being counted in the safety totals, inflating
them by 3,600 retrievals while contributing no new synthesis replicates. They are now excluded as a
reproducibility check rather than independent evidence. Safety totals are therefore **62,310
retrievals and 77.3M decode attempts** (still zero wrong outputs, zero false accepts; upper limits
5.9e-5 per retrieval, 4.8e-8 per attempt, 6.3e-3 per synthesis replicate).

**Verification.** `audit/verify_manuscript.py` extended with 20 checks over the new experiment,
including the round-robin-identity assertion. **All 125 checks pass.**

**Numbers changed?** Yes, and deliberately: a new experiment was added and one sub-claim was weakened to
match it. No existing result file was altered; no prior measurement changed value. The safety totals
changed because a double-count was removed.
