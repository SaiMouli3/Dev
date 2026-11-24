# Change log: revision 9 (reviewer document "Fix errors and inconsistencies first")

This revision applies the reviewer document, section by section. Unlike revision 8 (language only), it
changes results. Some inconsistencies were corrected, new baselines were run, and material was moved to a
supplement. Every number in the paper is regenerated from committed result files, and
`audit/verify_manuscript.py` recomputes the claims (see section 6).

Outputs: `manuscript/main.tex` + `manuscript/supplement.tex` (new). The final bundle in `final/` holds both PDFs.

## 1. Errors and inconsistencies (reviewer section A)

| # | Reviewer item | Finding when checked against the result files | Fix |
|---|---|---|---|
| A1 | Sec. VI-M says the ranking is unchanged under the cost models; Table XVIII says per-batch charges "change the ranking" | Confirmed. The cost is `reads + 50*batches + 500*libraries`, with batches proportional to reads and every policy preparing all six libraries, so the ranking cannot change | Table XVIII row now reads "leave the ranking unchanged but shrink the margin allocation competes for" |
| A2 | Cost table (XVII) disagrees with Table V | Confirmed, and worse than reported. The table took the first row per policy from `r4_cost_model.csv`, which listed the 70-synthesis extension first (`dt_baseline_r4ext*.csv`, 210 retrievals). Frozen rules and static plans came from that extension, while the dropout-model rows came from the original 300 retrievals, so one table mixed two retrieval sets | `r4_analysis.cost_model` and `make_r4_tables.cost_model` now cost the same 300 retrievals as Tables V and VI. The switch column is filled only where switches were recorded on those retrievals: dropout-model rules (38–59), and round-robin (63) from its bit-identical ρ_U=0.03 re-run. Switch counts for the frozen rules (38–61) and static plans (65, 69) are quoted in the text as coming from the extension |
| A3 | UCB+dropout q99 difference has two intervals (−1,147 to +150 vs −1,001 to +242) | Confirmed. Root cause: each analysis script drew bootstrap resamples from one sequential random stream, so the same statistic got a different Monte-Carlo interval depending on script and call order. Re-running the committed pipeline unchanged did not reproduce its own `r6_likelihood_sweep.csv` intervals | New `statistical_analysis/bootstrap.py`: resampling indices depend only on the number of syntheses. `dt_analysis.py`, `r4_analysis.py`, `r9_analysis.py`, `fill_text_numbers.py` and the verifier all use it. The tails table now has a paired Δq99 column, and the text quotes only that column: UCB+dropout −450 reads (95% CI −1,001 to +244) |
| A4 | Fixed-proportion q99 "4.70k / −450 / −8.7%" mixes designs | Confirmed. 4.70k is the 100×3 design; in the primary 30×10 design static α=1 has q99 4.80k | Text, Fig. 6 and Table XVIII now give both: −350 reads (95% CI −450 to −3, −6.8%) in 30×10, and −450 (95% CI −550 to −202, −8.7%) in 100×3. The 30×10 interval excludes zero only narrowly, and the text says so |
| A4+ | (found while fixing A4) | The text claimed coverage+dropout q99 "4.69k (−9.0%)" and Thompson+dropout "4.75k (−7.8%)" in the 100×3 design. The paper's own tails table printed 5.05k and 5.10k. The text also said "every other adaptive rule is significantly worse at q99", which is false for Thompson+dropout and pivotal+dropout | Text now uses the table values (−100 and −49 reads, intervals containing zero). The claim is restricted to the primary-model rules, which are all worse with intervals excluding zero |
| A5 | "No adaptive rule beat round-robin by more than a few per cent in any configuration" contradicts −34%, −28%, −45%, −16% | Confirmed | Restricted to the baseline everywhere (abstract, Sec. VI-C, Conclusion). The text also says larger gains appear only in dropout-dominated scenarios or with a calibrated prior |
| A6 | Failed-container intervals differ between Table VI, VI-G and Table XI | Confirmed, plus more: the VI-G text had stale intervals for coverage-greedy, Thompson, learn-then-commit and information-greedy, and Tables V and VI printed different intervals for the same means | Shared bootstrap (A3). Paired comparisons in `r4_analysis` and `r9_analysis` now also resample every synthesis (as `dt_analysis` does), not only those with a joint success. All failed-container intervals in the text are filled from `dt_summary.csv`, the source of Table VI |
| A7 | "Any one of three removes it" vs dropout-only +1.6% | Confirmed and broader. Single relaxations leave pivotal at +4 to +12%, tail-aware at +10% (dropout alone) and learn-then-commit at +13%. The abstract's "every rule 0 to 4.5% below round-robin" was false | "removes most of the penalty, leaving the rules between 4% below and 13% above round-robin". Fixed in the abstract, Sec. VI-C, the ablation caption and the evidence table |

Estimates did not change; only the bootstrap intervals moved. Interval limits moved by at most 0.8
percentage points on the baseline (excluding information-greedy), 1.3 in the failure scenarios, and 3.6 for
information-greedy (e.g. after a failure, [+50, +130] → [+51, +127]).

## 2. Headline framing (reviewer section B)

- B1: abstract, Sec. I and Conclusion now say that verified stopping *consumes 27–32% fewer reads than the
  capacity a non-stopping plan must reserve* (no more "reduces the bill by").
- B2 (new run): a **per-read stop** applies the same verified test after every read, with the same
  container schedule. On regenerated baseline streams it used 0.7% fewer reads than batched round-robin
  (95% CI −0.7 to −0.6%), so batches of 50 lose under 1%. Implemented as `StreamRetriever(per_read_stop=True)`.
- B3 (from committed data): **margin reservations**. A fixed-depth budget of 1.2× median (4.08k) completes
  85% of retrievals (95% CI 80–90%); 1.2× mean (4.19k) completes 89%. With retries these cost 4.78k and 4.69k
  reads per completed retrieval, against 3.49k with verified stopping. New Table IV(a), `r9_analysis.py`.

## 3. Allocation baselines (reviewer section C)

- C1 (new runs): **successive elimination** and a **yield-threshold** rule (skip a container while its
  recent yield is below 25% of the best's). Both were specified before running and not tuned; see
  `stream_retrieval.py`.
  - On the regenerated baseline: successive elimination never eliminated a container (identical to
    round-robin), and the yield threshold gave −1.0% (−1.6 to −0.6%).
  - Homogeneous control: both identical to round-robin.
  - One failed container: −0.5% and −3.8%.
  - No wrong output in 3,600 retrievals.
  - These rules run on regenerated streams, because DT4DDS is not seed-deterministic and the committed
    streams are not in the repository. They are therefore paired with round-robin re-run on the same
    streams, in a separate table block (Table IV(b)), not added to Tables V/XI.
  - Regenerated one-failed success is 50%, against 34% on the committed streams. The text reports this.
- C2: information-greedy removed from every main table. It keeps one sentence in the text (441%) and a full
  supplementary table (S3) with all its results.
- C3: the allocation section now opens with the ρ_U=0.03 comparison and the ablation, then gives the
  primary comparison.

## 4. Cuts and structure (reviewer section D)

- Abstract cut from 516 to about 200 words.
- New `supplement.tex` (10 pages):
  - Literature table, operational specification table, information-greedy, conditional/unconditional cost
    after a failure.
  - Posterior calibration and ρ_U re-selection (old Table X, Fig. 4), completion forecast (Fig. 5),
    severity sweep (Fig. 3).
  - Prior sensitivity, consensus cap (XIV), failure modes (XV), geometries (XVI), nanopore result (VI-N).
  - Leaky upper bound (Appendix A, Table XX), full deviations (most of Appendix B), claim-to-evidence map.
- Cross-references between the two documents use `xr`.
- VI-E + VI-F merged into "Why adaptive allocation loses, and the dropout repair". VI-D, VI-G and VI-L
  merged into "Robustness".
- Process language removed: "frozen" → "primary", "Revision 4" → "secondary analyses",
  "superseded" → "replaced", "our earlier rule" removed. "Leak-Free" was removed from the title.
- Fig. 6 cut to four recommendations.
- Main paper: 25 pages (single column), down from 34.

## 5. Additions (reviewer section E)

- E1 (Sec. VII-A): why stopping gains more than allocation.
  - The stopping saving is 1 − E[T]/q99(T). On the baseline the coefficient of variation is 0.16 and the q99
    is 48% above the mean. 69% of the read-count variance is within a synthesis (read order) and 31% between
    syntheses.
  - A fluid bound on allocation is 1 − ē/e_max = 11% with the measured container rates, against the observed
    4.2–4.6%.
  - This is stated as an approximation, not a theorem. The reviewer's suggested wording ("the max over stripes
    has a heavier tail than its mean") was not used, because the maximum of many stripes has a light tail.
- E2 (Sec. VII-B): scaling. The expected trend is argued, and checked on the existing geometries: stopping
  saving 24% with a 4× file, 18% at (32,64) and 41% with 12 containers. These geometries have 30–60
  retrievals, and the text flags that.
- E3 (Sec. V-C): three-tier statement (confirmatory / secondary / exploratory) taken from the protocol files.
  It states that intervals are per-comparison and not adjusted for multiplicity.

## 6. Verification

- `audit/verify_manuscript.py`: 167 checks, all pass. Added checks cover the reservation baselines,
  the regenerated-stream baselines (including that successive elimination is identical to round-robin on the
  baseline), the paired q99 intervals of A3/A4, the cost table using Table V's retrievals, switch counts by
  source, failed-container interval identity across tables, and the Sec. VII numbers.
- Safety totals are unchanged: 62,310 retrievals, 580 synthesis units. The regenerated-stream runs are
  excluded like the clean-reproduction runs, because they reuse replicate labels.
- Text numbers that the revision added are written by `experiments/scripts/fill_text_numbers.py` from the
  same CSVs as the tables. The values are recorded in `manuscript/text_numbers.tex`.
- Both documents compile with pdfLaTeX + BibTeX, with no undefined references. The Source zip contains both.

## 7. Not done, or done differently from the reviewer's wording

- The new baselines were not added as rows of Tables V and XI. Those tables are paired on the committed
  streams, which cannot be regenerated bit for bit, so mixing would recreate the A2 problem. They are in
  Table IV(b) instead.
- The reviewer's "keep in main text" list was followed, with the following additional tables kept:
  - notation;
  - tails (needed for the q99 claims of A3/A4);
  - cost, with the A1/A2 fixes;
  - limitations;
  - the new benchmark table.
- No real-data or wet-lab experiment was added; the paper remains computational.

## 8. Addendum: ACM journal format, one document

At the authors' request the paper was converted from IEEEtran to ACM's `acmart` class.
- **Format:** `acmsmall`, ACM's single-column journal format (e.g. ACM Transactions on Storage), with a page
  size of 6.75 x 10 in (486 x 720 pt). ACM's two-column conference format (`sigconf`, US Letter) can be
  selected by changing one class option.
- **One document:** the former supplement is now a set of appendices in the same PDF. Cross-references that
  used to point into the supplement now point to the appendices.
- **ACM front matter:** CCS concepts (to be checked with the ACM CCS tool), keywords, the ACM reference
  bibliography style, and numbered citations.
- **Placeholders, not inventions:**
  - Author affiliations and countries are visible placeholders, because ACM requires them and none was
    provided.
  - Journal, volume, article number and rights are not set. The placeholder running footer and the ACM
    reference block are switched off until ACM assigns them.
- **Tables** wider than the ACM text block are scaled down to fit. No table content changed.
- **Build and checks:** 38 pages, no LaTeX errors, no unresolved references. All 167 verifier checks pass on
  the ACM PDF. The BibTeX warnings are fields ACM's style expects but some entries lack (publisher address,
  page numbers); these were not invented.
- **IEEE version:** the IEEE sources are kept in `manuscript/ieee_version/` for reference.

## 9. Addendum: submission fixes

- **References:** publisher locations were added for standard publishers (ACM, IEEE, Cambridge University
  Press, Curran Associates for NeurIPS, SIAM for SODA). BibTeX warnings fell from 20 to 7. The remaining
  seven cannot be checked from this environment, because the network policy blocks arxiv.org, doi.org,
  papers.nips.cc, usenix.org and api.crossref.org:
  - page numbers for cao2025optimizing, cao2026modeling, ford2010availability, garivier2016explore,
    nath2006subtleties and rashtchian2017clustering;
  - the publisher city for efron1993introduction.
  They were left empty rather than filled from memory.
- **Preprint status:** kong2026coverage (arXiv:2610.02931) is still cited as a preprint, and its publication
  status could not be checked. levy2026expected already cites its ISIT 2026 version.
- **CCS concepts:** the three concepts could not be checked against the ACM CCS site (blocked). They should be
  confirmed at https://dl.acm.org/ccs.
- **Anonymous version:** `final/ACM_Research_Paper_Anonymous.pdf` is built by finalize.sh with
  `review,anonymous`. It hides authors, affiliations, acknowledgements and the repository link, and adds line
  numbers. finalize.sh fails if the anonymous PDF still contains any author name, institution or repository
  name.
