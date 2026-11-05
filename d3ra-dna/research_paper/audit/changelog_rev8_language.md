# Change log: revision 8 (language and presentation only)

This revision changes only the wording and presentation of the paper. No experiment was re-run, and no result file,
table entry, figure, equation, algorithm, parameter, dataset identifier, citation or reference was changed.
Every number in the original manuscript is still in the revised one, and no new number was added (see
section 3).

Files edited: `manuscript/main.tex` (prose only) and `manuscript/tables/limitations.tex` (caption and one
cross-reference phrase; no table data). The final bundle in `final/` was rebuilt from these files.

## 1. Edits by instruction

| # | Instruction | What changed |
|---|---|---|
| 1 | Repetitive rhetorical framing | Removed or rewrote as direct statements: "The honest summary is therefore…" (VI-A), "Three consequences, stated plainly." (VI-J), "The honest statement about every zero- or partial-success result…" (VI-L), "The right reading of our negative result is therefore not…" (VI-F), "the defensible claim is…" (abstract), "The honest claim is therefore…" (VI-H), "the honest estimate" (VI-E), "the honest statement for the coverage family" (VI-J), "the conservative reading" (VI-A). Contribution 3 "…with honest uncertainty" became "…with explicit uncertainty". |
| 2 | Reviewer-facing and meta language | "The reviewer-facing question is…" became "We next asked whether…". "What is and is not new here." became **Relation to prior work**, with the same three points stated without "We do not claim". "What we would test with wet-lab containers." became **Future experimental validation**, keeping label `sec:wetlab` and all three proposed tests plus the decisive experiment. Removed from the abstract: "Every number is regenerated from committed result files under frozen protocols with logged deviations." Reproducibility information stays in full in Section IX and Appendix B. "We state the relation precisely", "so that nothing has to be read out of the code", "One finding here is negative and we report it rather than the claim we expected to make" were removed. |
| 3 | Self-auditing language | "leak-free" is kept only where it names a condition: the table label "leak-free calibrated" and the leaky vs leak-free comparison in Appendix A. Elsewhere it reads "calibrated belief" or "calibration and evaluation reads kept disjoint". "Leakage audit of the replay" became **Separation of calibration and evaluation reads**, with the content unchanged. "We make no claim that these numbers are usable probabilities" became "the forecasts should not be used as probabilities". In the main text, "revision-4" labels became "secondary"/"additional". Appendix B keeps "Revision 4" because there it identifies `dt_protocol_rev4.yaml`. All of these are kept: pre-registration, frozen protocol, protocol deviations, cluster-bootstrap CIs, Clopper–Pearson bounds, per-synthesis bounds and the leakage methodology. |
| 4 | LLM-like scaffolding | Rewrote: "This is the clearest case in the paper where…" (VI-K), "That sentence needs care in both halves…" (VI-H), "Two things this table does not establish." (VI-A), "Three qualifications keep this from being a positive headline." (now "Three qualifications apply."), "Three quantities must stay separate and we never equate them" (now "Three quantities are distinct"), "Several distinct quantities are easy to conflate…" (III-C), "Three static numbers appear in this paper and they are easy to confuse…" (VI-B), "Before the numbers, which experiment uses which" (VI-I), "If one sentence has to carry the paper's practical content…" (VII). Shortened two subsection titles: "What a 'container switch' means here, and what it would mean in a laboratory" became "Physical interpretation of container switching", and "Robustness: what does and does not change the conclusions" became "Robustness". "What the primary comparison looks like under the corrected likelihood" became "Primary comparison under the re-selected likelihood". Labels and numbering are unchanged. |
| 5 | Practical recommendation box (Fig. 6) | All seven recommendations kept in the original order and with all their numbers (27–32%, ~9%, 24%, 0.313→0.357, 7–441%, 3.5%, 16%, 1.8M, 2.8×, ~17 and ~10 retrievals, 30k, 27%/16%, (32,48), K). Each now opens with a plain imperative heading followed by short sentences, replacing the slogan-style wording. |
| 6 | Abstract | Went from 603 to 516 words. Long sentences were split and the restatement "So adaptive allocation is worth a few percent at best, against stopping's 27–32%" was removed. Kept: problem, method, all main numbers (3.49k, 3.37–3.60k, 4.80–5.15k, 27–32%, 7–441%, ρ_U 0.1/0.03, −1.2 to +1.5%, 4.2% [3.5–4.9%], 0–4.5%, −3.5%, 1–9%, −450 / −8.7%, safety macros, 16%, 1.8M, 2.8×, 17), the digital-twin vs replay distinction, the negative findings and the limitations. |
| 7 | Transitions and repetition between sections | The one-failed-container success facts (34%, identical across policies) are now stated once in VI-D and referenced from VI-G. The two-failed-container arithmetic is stated once in VI-G (new label `sec:res-fail`) and referenced from VI-C. The duplicated second safety paragraph in VI-A was merged into the first. The duplicated oracle-decodability description in V-E was merged into one paragraph. The repeated "spread is a few per cent of a bill set by run setup" sentence in VI-M was removed. The conclusion was condensed into three paragraphs. |
| 8 | Objection pre-emption | Repeated "this is computational, not wet-lab" statements were reduced. The statement is kept in the abstract, Introduction (Scope), Limitations, the future-validation paragraph and the Conclusion's final paragraph. The Introduction's "no result should be read as a wet-lab demonstration" and the Related Work regret/optimality disclaimer were removed; the latter is still in Section VIII and Table XVIII. Every limitation remains. |
| 9 | Authors | Added J Geetha (geetha@msrit.edu) as second author next to Sai Mouli (saimouli12345678@gmail.com), using the standard IEEEtran `\and` block. No affiliation, department or title was invented. A LaTeX comment marks where each affiliation should go once the authors confirm it. No venue name was added. |
| 10 | Phrase stacking | Split the multi-claim sentences in the Findings paragraph, the stopping summary ("It needs no model, never returned a wrong file, and converts…"), and the ablation, dropout, calibration, forecast, safety and tail paragraphs. Every underlying finding is kept. |

## 2. What was deliberately not changed

- The title, including "A Leak-Free Computational Study". It is the author's choice; see the optional suggestion in section 4.
- All `tables/*.tex` data. Only the caption and one cross-reference phrase in `limitations.tex` were reworded.
- Equations (1)–(2), the TikZ pipeline figure, all `\includegraphics`, all algorithms (Table `tab:algorithms`) and `references.bib`.
- The header comment of `main.tex`.

## 3. Verification

| Check | Result |
|---|---|
| Numeric recomputation (`audit/verify_manuscript.py`, rebuilt PDF) | **125 / 125 checks pass** (`audit/verify_manuscript.txt`) |
| Numbers: original vs revised `.tex` (comments excluded) | 0 numbers removed; 0 numbers added. Some values now appear fewer times because a duplicated restatement was removed (e.g. "660 retrievals" 2→1) |
| Equations | Both `equation` environments byte-identical |
| Citations | Same 51 citation keys; none added or removed |
| Cross-references | All 31 original labels kept; 1 label added (`sec:res-fail`); no undefined references (`??` absent from PDF) |
| Tables / figures | Same 21 `\input` tables, 4 graphics files, 3 `figure` + 3 `figure*` environments, 1 TikZ figure |
| Section structure | 13 `\section`, 27 `\subsection`, appendices unchanged; numbering unchanged |
| Inline math | Nothing added. Three duplicated occurrences were removed along with redundant sentences: `$-$`, `$\rho_U$`, `$(32,48)$` |
| Compilation | pdfLaTeX + BibTeX: 26 pages (same as before), IEEEtran conference two-column layout, no errors. Font and overfull-box warnings are the same 7 as in the previous build |
| Authors | Sai Mouli / saimouli12345678@gmail.com; J Geetha / geetha@msrit.edu |

Typographic note: `No~DNA~was` in the Conclusion uses non-breaking spaces. Without them a line break splits
the phrase in the `pdftotext` output and the verifier's "no wet-lab claim" text check fails.

## 4. Items for author decision (not changed, per the final rule)

These are existing points of tension in the manuscript. Resolving any of them would change a scientific
statement or its scope, so they were left as they were:

1. **Section VI-I vs VI-J.** VI-I says that every policy run in this paper uses ρ_U = 0.1 and that no primary
   conclusion changes. VI-J reports the complete primary comparison re-run at ρ_U = 0.03 and says that the
   likelihood setting changes part of the allocation result. A possible fix: limit VI-I's statement to "every
   policy run in Sections VI-A–VI-H and the secondary scenarios" and point to VI-J.
2. **UCB + dropout q99 interval.** VI-F gives the paired q99 difference as −450 reads (95% CI −1,147 to
   +150). VI-H gives −450 (95% CI −1,001 to +242). One of the two may come from an earlier bootstrap run.
   The result files should decide which is correct.
3. **Recommendation 4 (Fig. 6)** says "every such rule we tested read 7–441% more than round-robin" without
   the qualification, given in the abstract, Findings and VI-J, that this holds at the pre-registered
   ρ_U = 0.1.
4. **Abstract length.** At 516 words it is still well above typical IEEE limits (about 150–250 words). Going
   further means dropping some secondary numbers (e.g. the 0.99-quantile or calibration-cost details), which
   the authors should choose.
5. **Title (optional).** "A Leak-Free Computational Study" could become "A Computational Study". The
   leakage control is described in the body.
6. **Affiliations.** Add once confirmed (placeholders are LaTeX comments in the `\author` block).
