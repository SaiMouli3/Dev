# Citation Verification Audit

Date: 2026-11-02. File audited: `research_paper/manuscript/references.bib` (44 entries: 41 publications and 3 datasets).

**Method.**
- **DOI items.** Metadata was fetched from `https://api.crossref.org/works/<DOI>`. The JSON fields checked were title,
  author list, container title, volume, issue, page or article number, and date. BibTeX was obtained from
  `https://api.crossref.org/works/<DOI>/transform/application/x-bibtex` (or `doi.org` content negotiation), and then
  normalised:
  - keys of the form authorYEARword;
  - Unicode converted to LaTeX/ASCII;
  - acronyms protected with braces;
  - `url`, `ISSN` and `month` fields dropped;
  - for article-number journals, the article number placed in `pages`.
- **arXiv items.** Checked through the arXiv Atom API.
- **Items without a DOI.** Checked on the publisher or proceedings landing page (USENIX, JMLR, NeurIPS) or on an
  institutional repository.
- **Datasets.** Checked through DataCite (`doi.org`) and the Zenodo or ENA APIs.
- **Test build.** All 44 entries compile with `IEEEtran.bst` (`\nocite{*}`) with no BibTeX warnings.

| bib key | DOI / arXiv | verification source (URL) | verified fields | status |
|---|---|---|---|---|
| church2012next | doi:10.1126/science.1226355 | https://api.crossref.org/works/10.1126/science.1226355 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED (single page 1628) |
| goldman2013towards | doi:10.1038/nature11875 | https://api.crossref.org/works/10.1038/nature11875 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| grass2015robust | doi:10.1002/anie.201411378 | https://api.crossref.org/works/10.1002/anie.201411378 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| bornholt2016dna | doi:10.1145/2872362.2872397 | https://api.crossref.org/works/10.1145/2872362.2872397 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| erlich2017fountain | doi:10.1126/science.aaj2038 | https://api.crossref.org/works/10.1126/science.aaj2038 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| organick2018random | doi:10.1038/nbt.4079 | https://api.crossref.org/works/10.1038/nbt.4079 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| antkowiak2020low | doi:10.1038/s41467-020-19148-3 | https://api.crossref.org/works/10.1038/s41467-020-19148-3 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| heckel2019characterization | doi:10.1038/s41598-019-45832-6 | https://api.crossref.org/works/10.1038/s41598-019-45832-6 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| chen2020quantifying | doi:10.1038/s41467-020-16958-3 | https://api.crossref.org/works/10.1038/s41467-020-16958-3 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| gimpel2023digital | doi:10.1038/s41467-023-41729-1 | https://api.crossref.org/works/10.1038/s41467-023-41729-1 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| shomorony2021dna | doi:10.1109/TIT.2021.3058966 | https://api.crossref.org/works/10.1109/TIT.2021.3058966 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| weinberger2022dna | doi:10.1109/TIT.2022.3176371 | https://api.crossref.org/works/10.1109/TIT.2022.3176371 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| lenz2020coding | doi:10.1109/TIT.2019.2961265 | https://api.crossref.org/works/10.1109/TIT.2019.2961265 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| press2020hedges | doi:10.1073/pnas.2004821117 | https://api.crossref.org/works/10.1073/pnas.2004821117 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| reed1960polynomial | doi:10.1137/0108018 | https://api.crossref.org/works/10.1137/0108018 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| barlev2025cover | doi:10.1109/TIT.2024.3496587 | https://api.crossref.org/works/10.1109/TIT.2024.3496587 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| gruica2024combinatorial | doi:10.1109/ISIT57864.2024.10619151 | https://api.crossref.org/works/10.1109/ISIT57864.2024.10619151 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| cao2025optimizing | doi:10.1109/ISIT63088.2025.11195510 | https://api.crossref.org/works/10.1109/ISIT63088.2025.11195510 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED (Crossref pages "1-6" are a placeholder and were omitted) |
| barlev2025deep | doi:10.1038/s42256-025-01003-z | https://api.crossref.org/works/10.1038/s42256-025-01003-z | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| sabary2024reconstruction | doi:10.1038/s41598-024-51730-3 | https://api.crossref.org/works/10.1038/s41598-024-51730-3 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| loose2016real | doi:10.1038/nmeth.3930 | https://api.crossref.org/works/10.1038/nmeth.3930 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| payne2021readfish | doi:10.1038/s41587-020-00746-x | https://api.crossref.org/works/10.1038/s41587-020-00746-x | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED (Crossref BibTeX year is 2020 (online); the print issue year 2021 is used) |
| weilguny2023dynamic | doi:10.1038/s41587-022-01580-z | https://api.crossref.org/works/10.1038/s41587-022-01580-z | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| lopez2019dna | doi:10.1038/s41467-019-10978-4 | https://api.crossref.org/works/10.1038/s41467-019-10978-4 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| leong2012distributed | doi:10.1109/TIT.2012.2191135 | https://api.crossref.org/works/10.1109/TIT.2012.2191135 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| weil2006crush | doi:10.1109/SC.2006.19 | https://api.crossref.org/works/10.1109/SC.2006.19 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED (Crossref "page" 31 is the article number; author middle initials not in Crossref) |
| thompson1933likelihood | doi:10.1093/biomet/25.3-4.285 | https://api.crossref.org/works/10.1093/biomet/25.3-4.285 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED (Crossref title is in all caps; converted to title case) |
| delahaye2021sequencing | doi:10.1371/journal.pone.0257521 | https://api.crossref.org/works/10.1371/journal.pone.0257521 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| wick2019performance | doi:10.1186/s13059-019-1727-y | https://api.crossref.org/works/10.1186/s13059-019-1727-y | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| cliff1993dominance | doi:10.1037/0033-2909.114.3.494 | https://api.crossref.org/works/10.1037/0033-2909.114.3.494 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED |
| srinivasavaradhan2021trellis | doi:10.1109/ISIT45174.2021.9517821 | https://api.crossref.org/works/10.1109/ISIT45174.2021.9517821 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED (DOI found by Crossref bibliographic search; the arXiv:2107.06440 record says "Extended version of paper presented at ISIT 2021") |
| kashanna2026reliability | doi:10.1109/TMBMC.2026.3690724 | https://api.crossref.org/works/10.1109/TMBMC.2026.3690724 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED (Crossref splits the surname as Hanna/"Serge Kas"; corrected to Kas Hanna (as on arXiv)) |
| efron1993introduction | doi:10.1007/978-1-4899-4541-9 | https://api.crossref.org/works/10.1007/978-1-4899-4541-9 | title, authors, venue, volume, issue, pages/article no., year (BibTeX via Crossref content negotiation) | VERIFIED (the Crossref record is the Springer US (reissued) book; the original 1993 print was Chapman & Hall (not in Crossref)) |
| levy2026expected | arXiv:2602.07601 | http://export.arxiv.org/api/query?id_list=2602.07601 | title, 4 authors (Adi Levy, Roni Con, Eitan Yaakobi, Han Mao Kiah), date 2026-02-07, v1, no journal ref | VERIFIED (preprint; note the first author is Adi Levy, not D. Levy as in docs/literature_review.md) |
| kong2026coverage | arXiv:2610.02931 | http://export.arxiv.org/api/query?id_list=2610.02931 | title, 4 authors, date 2026-10-02, comment '33 pages', abstract | VERIFIED (preprint) |
| batu2004reconstructing | none (no DOI) | https://researchonline.lse.ac.uk/id/eprint/8146 ; https://repository.upenn.edu/cis_papers/173 | title, 4 authors, SODA 2004 (New Orleans, 11-13 Jan 2004); pages 910-918 from the UPenn repository record (seen as a search snippet; the page itself returned 403) | VERIFIED (institutional repository; pages from secondary source) |
| rashtchian2017clustering | none (no DOI) | https://papers.nips.cc/paper_files/paper/2017/hash/ab7314887865c4265e896c6e209d1cd6-Abstract.html | title, 8 authors, NIPS 2017, Advances in NeurIPS 30 | VERIFIED (proceedings page; page numbers omitted because they were not shown on the page) |
| nath2006subtleties | none (no DOI) | https://www.usenix.org/conference/nsdi-06/subtleties-tolerating-correlated-failures-wide-area-storage-systems ; PDF https://www.usenix.org/legacy/event/nsdi06/tech/full_papers/nath/nath.pdf | title, venue NSDI '06, year; 4 authors read from the PDF title page | VERIFIED (the USENIX landing-page BibTeX lists only H. Yu; the author list was taken from the paper PDF) |
| ford2010availability | none (no DOI) | https://www.usenix.org/conference/osdi10/availability-globally-distributed-storage-systems ; PDF https://www.usenix.org/legacy/event/osdi10/tech/full_papers/Ford.pdf | title, 8 authors, OSDI '10, year (USENIX BibTeX and PDF title page) | VERIFIED (pages not listed by USENIX and omitted) |
| cidon2013copysets | none (no DOI) | https://www.usenix.org/conference/atc13/technical-sessions/presentation/cidon | title, 6 authors, USENIX ATC '13, pp. 37-48, year (USENIX BibTeX) | VERIFIED |
| bubeck2013optimal | arXiv:1207.5259 (no DOI) | https://www.jmlr.org/papers/v14/bubeck13a.html ; http://export.arxiv.org/api/query?id_list=1207.5259 | title, 3 authors, JMLR 14(17):601-623, 2013 | VERIFIED |
| gimpel2023prjeb65931 | ENA PRJEB65931 | https://www.ebi.ac.uk/ena/portal/api/search?result=study&query=study_accession=PRJEB65931&fields=study_accession,study_title,center_name,first_public | study title, center (FML), first public 2023-09-08 | VERIFIED (dataset) |
| technion2025dnaformerreads | doi:10.5281/zenodo.17473983 | https://doi.org/10.5281/zenodo.17473983 (DataCite BibTeX) ; https://zenodo.org/api/records/17473983 | title, creator, publication date 2025-10-29, publisher Zenodo | VERIFIED (dataset) |
| barlev2024binned | doi:10.5281/zenodo.14296588 | https://doi.org/10.5281/zenodo.14296588 (DataCite BibTeX) ; https://zenodo.org/api/records/14296588 | title, 5 creators, date 2024-12-07; the description lists Srinivasavaradhan et al. among the included datasets | VERIFIED (dataset; project docs had attributed this record to Srinivasavaradhan et al.) |
| (excluded) birnbaum1969importance | none | Crossref bibliographic search (no hit); web search (secondary citations only) | n/a | NOT VERIFIED - left out of references.bib |
| (excluded) holm1979simple | 10.2307/4615733 (asserted by secondary sources) | https://api.crossref.org/works/10.2307/4615733 returns 'Resource not found'; JSTOR page not fetchable | n/a | NOT VERIFIED - left out of references.bib |
| (dropped) efron1979bootstrap | doi:10.1214/aos/1176344552 | https://api.crossref.org/works/10.1214/aos/1176344552 | title, author, Ann. Statist. 7(1):1-26, 1979 | VERIFIED but not included; the brief allowed 'Efron 1979 or Efron and Tibshirani 1993', and the book was kept. Can be added on request |

**Corrections made relative to the earlier project notes**
- `levy2026expected`: the first author is **Adi** Levy. `docs/literature_review.md` §3 says "D. Levy".
- `barlev2024binned`: Zenodo 14296588 was published by Bar-Lev, Orr, Sabary, Etzion and Yaakobi. It packages the
  Srinivasavaradhan et al. reads, so cite both.
- `srinivasavaradhan2021trellis`: the ISIT DOI 10.1109/ISIT45174.2021.9517821 was not recorded in the project notes.
- `kashanna2026reliability`: the journal version is now published (IEEE TMBMC 12:585–598, 2026). The ISIT 2025 version
  is doi:10.1109/ISIT63088.2025.11195362, which was verified through a Crossref search but is not included.
- `payne2021readfish`: the issue year is 2021; Crossref BibTeX gives the online year 2020.

## Usage audit (manuscript)

- **Date:** 2026-11-02.
- **Scope:** every `\cite` in `manuscript/main.tex` (commit 89e439f) and `manuscript/tables/*.tex`.
- **Method:** each citation was checked against what the cited work says, as summarised in `supplementary/literature_review.md`. Where that file says "abstract only" or "snippet level", the depth of reading limits any verdict.
- **Line numbers** refer to `main.tex`.

### A. Key resolution

- All **44** distinct cite keys used in `main.tex` and the tables are present in `references.bib`.
- All 44 bib entries are cited at least once.
- No undefined or unused keys.

### B. Sentence-level usage findings

| # | Line | Sentence (abridged) | Keys | Finding | Severity / suggested fix |
|---|---|---|---|---|---|
| U1 | 80–81 | DNA "demonstrated as an archival medium at megabyte to hundreds-of-megabyte scale" | church2012next, goldman2013towards, grass2015robust, erlich2017fountain, organick2018random | Grass 2015 stored 83 kB; Church about 0.65 MB and Goldman about 0.74 MB. Only Erlich (2.14 MB) and Organick (> 200 MB) are at MB scale or above. | Minor. Write "kilobyte to hundreds-of-megabyte". |
| U2 | 83–84 | "every read carries substitution, insertion and deletion errors" | heckel2019characterization, gimpel2023digital | The cited works characterise per-base error *rates*; they do not say that every read is erroneous. The paper's own data show 27–68% error-free reads (`realdata_channel.csv`). | Minor. Write "reads carry …". |
| U3 | 91–92 | "Capacity, random access and risk management push towards many physical containers (vials, capsules, pools)" | bornholt2016dna, organick2018random | Both works concern random access *within* a pool via PCR primers. Neither argues for many separately stored physical containers or for risk management (literature_review §1: Organick "no failures across pools or containers are modelled"; Bornholt "no reasoning about … multiple containers"). Grass 2015 (silica encapsulation as a container / failure domain) supports the physical-container notion better. | Moderate. Re-attribute or rephrase as the authors' motivation. |
| U4 | 93–95 | Theory of partitioned data with uniform sampling; expected recovery time with failures | kong2026coverage, levy2026expected | Consistent with the review. Kong was read at **abstract level only**. | OK. Note that Kong was read at abstract level. |
| U5 | 112–113 | "Stopping when decodable and checksum gating are known practice" | lopez2019dna | Lopez supports *sequence-until-decodable* (full text). It is not cited for checksum or verification gating, and the review notes that automation is unclear. The precedent for verification-gated acceptance in the review is erlich2017fountain (droplets failing the RS check are rejected). | Moderate. Add `\cite{erlich2017fountain}` for checksum gating, or rephrase. |
| U6 | 135 | Systems combine constrained encodings with inner and outer codes | antkowiak2020low, goldman2013towards, grass2015robust, erlich2017fountain, organick2018random, press2020hedges | Goldman uses a rotating code with 4-fold overlap (repetition), not inner/outer codes. | Minor. Drop goldman from this sentence or write "with redundancy". |
| U7 | 139–140 | Channel modelled as sampling with replacement followed by noise | shomorony2021dna, weinberger2022dna, lenz2020coding | Shomorony and Weinberger fit. Lenz studies *coding over sets* (unordered strands), with bounds and constructions, rather than a sampling-with-replacement channel model. | Minor. Move lenz to a "coding over unordered sets" clause. |
| U8 | 142–144 | Bar-Lev: for an uncontrollable sampler, MDS and uniform sampling minimise E[reads] = n(H_n − H_{n−k}) | barlev2025cover | Matches the review. | OK. |
| U9 | 145–147 | Follow-ups on random access, real-data sampling, partitioning, failures; reliability of MDS-coded retrieval | gruica2024combinatorial, cao2025optimizing, kong2026coverage, levy2026expected, kashanna2026reliability | Consistent. Gruica, Cao and Kong were read at abstract level. "Real-data sampling distributions" for Cao means a log-normal fitted to the authors' own data (the table footnote says so). | OK. |
| U10 | 150–151 | Read Until and adaptive sampling select molecules in real time | loose2016real, payne2021readfish, weilguny2023dynamic | Consistent. Weilguny was read at abstract level. | OK. |
| U11 | 152–153 | DNA-storage readouts sequenced until sufficient coverage | lopez2019dna | Consistent (full-text quote in the review). | OK. |
| U12 | 153–155 | Container choice resembles adaptive coupon collection; analogues in bandits | thompson1933likelihood, bubeck2013optimal | Consistent as an analogy. Bubeck was read at **snippet level**. The manuscript cites Thompson, but the Thompson-sampling controller in the code is never evaluated. literature_review.md §8 also claims a "Gamma posterior", whereas the manuscript posterior is a discrete grid. | Minor. Fix the review text. Consider evaluating the Thompson variant, or cite Thompson only as an analogy. |
| U13 | 157–158 | "Erasure-coded storage limits the number of fragments of a stripe per failure domain" | weil2006crush, cidon2013copysets, ford2010availability | CRUSH and Copysets concern *replicated* data placement across failure domains; Ford covers rack-aware placement, including for stripes. Calling all three "erasure-coded" is imprecise. | Minor. Write "Distributed storage spreads replicas or erasure-coded chunks across failure domains". |
| U14 | 159 | "spreading can be worse than concentrating under correlated failures" | nath2006subtleties | The review (abstract and summary level) says rankings can *flip* under correlated failures and extra redundancy has diminishing returns. "Spreading worse than concentrating" is the authors' reading, not a stated result. | Moderate. Soften to "design rankings can reverse under correlated failures". |
| U15 | 160–161 | "the optimal allocation of coded storage across unreliable nodes was characterised by Leong et al." | leong2012distributed | Leong et al. characterise optimal (symmetric) allocations only in certain regimes (minimal vs maximal spreading). The general optimal allocation is left open (conjectured). "Was characterised" overstates. | Moderate. Write "was analysed, with optimal allocations identified in several regimes". |
| U16 | 161–162, 378 | Spreading rule is a special case of Leong et al. | leong2012distributed | Consistent with the review and an appropriately modest framing. | OK. |
| U17 | 164–166 | Trace reconstruction after clustering or learned models | batu2004reconstructing, srinivasavaradhan2021trellis, sabary2024reconstruction, rashtchian2017clustering, barlev2025deep | Consistent. | OK. |
| U18 | 166–167 | Systematic, context-dependent nanopore errors well documented | delahaye2021sequencing, wick2019performance | Consistent. Delahaye was read at abstract level. | OK. |
| U19 | 232 | "A systematic **Cauchy** Reed–Solomon code [reed1960polynomial]" | reed1960polynomial | Reed & Solomon (1960) introduce RS codes, not the Cauchy-matrix construction (Blömer et al., 1995). The citation supports "Reed–Solomon" only. | Minor. Add a Cauchy-RS reference, or cite RS for the code family only. |
| U20 | 333 | W_s is "a Birnbaum-type component importance" | (none) | Uncited attribution. Birnbaum (1969) was excluded from the bibliography because it could not be verified (§12 of the review, `missing_evidence.md` #10). | Moderate. Add the verified reference, or drop the attribution. |
| U21 | 378–381 | Placement rule "follows from the allocation analysis of [Leong]"; proof in supplementary | leong2012distributed | Consistent. | OK. |
| U22 | 392 | "two-sided bitwise-majority-alignment consensus [batu2004reconstructing]" | batu2004reconstructing | Batu et al. introduce (one-sided) BMA for deletion channels. The two-sided variant used here and its use with substitutions and insertions are not from that paper. | Minor. Write "an extension of BMA" and cite the source of the two-sided variant, if any. |
| U23 | 435–440 | Data sources | gimpel2023digital, gimpel2023prjeb65931, barlev2025deep, technion2025dnaformerreads, srinivasavaradhan2021trellis, barlev2024binned | Correct attributions (the dataset/paper split follows the corrections above). The Srinivasavaradhan / barlev2024binned data contribute no result reported in the manuscript; they appear only in the q3b audit. | Minor. Report the result, or remove these from Data. |
| U24 | 450–451 | Paired bootstrap [efron]; "effect sizes as Cliff's δ [cliff]" | efron1993introduction, cliff1993dominance | Efron is fine. **Cliff's δ is never reported** in the manuscript or the result files, so the citation supports a method that was not applied. | Moderate. Report δ, or remove the clause and the citation. |
| U25 | 639–640 | Round-robin degradation "consistent with coupon-collector theory" | barlev2025cover, levy2026expected | The works support the qualitative theory (E[reads] = n(H_n − H_{n−k}); Levy's leading-order 1/(surplus+1) scaling), but the manuscript computes no prediction to compare with the measured +5 to +32%. The citation is used as evidence of agreement without a check. | Moderate. Add the quantitative comparison, or write "qualitatively expected from". |
| U26 | `literature_comparison.tex` | Row-by-row Y/N/P | all 11 compared works | Entries match literature_review §13. Kong ("abstract only") and Bubeck ("snippet level") are characterised from limited reading; this is flagged in the table comment but not in the printed caption. The **"This work / Correlated errors = Y"** entry overstates: the correlated-error model is not coupled to the controller (l. 368–370). | Moderate. Change it to P, and add "Kong et al.: abstract only" to the caption footnote. |
| U27 | `literature_comparison.tex` caption | "Containers" in Leong are storage nodes; in Bubeck, sampling sources | leong2012distributed, bubeck2013optimal | Appropriately qualified. | OK. |

### C. Novelty and "first" claims

- **No "first", "novel" or "for the first time" claim** appears in `main.tex` or the tables.
- **l. 113–115:** "the container-selection component has no direct precedent we could find in DNA storage, but has analogues in adaptive sampling and bandits". Appropriately hedged, and consistent with literature_review §14 ("not found, not proof of absence").
- **l. 169–172 (gap statement):** a combination claim of (a)–(d). It is defensible given Levy (failures), Kong (multi-container, abstract only) and Weilguny (adaptive, other domain).
  - Kong et al. were read only at abstract level. Its full text should be checked before claiming that no prior multi-container work has heterogeneous containers.
  - The table rates Kong "P" for heterogeneity; "container sizes only" is inferred from the abstract.
- **l. 118, l. 447–448:** "strongest-possible static baseline", "no deployable static plan can do better". This is not a literature-novelty claim but an optimality overclaim: the baseline is a local Nelder–Mead optimum within the fixed-budget class. See `claim_verification.md` I4 / E6.
- **literature_review.md §14 vs manuscript.**
  - The review says "the allocation component gives a smaller, but measurable, benefit". On Twist this is about 2.0k reads (≈ 3.5 percentage points); on Genscript about 64k.
  - The manuscript abstract attributes the full 11–17% to "the controller" without this decomposition. That is a framing overclaim, not a citation problem, and is listed in `claim_verification.md` A1.

### D. Summary of required citation-usage fixes

1. **U3:** re-attribute the "many physical containers" motivation; Grass 2015 is a better fit than Bornholt / Organick.
2. **U5:** add erlich2017fountain for checksum or verification gating.
3. **U14:** soften the Nath et al. paraphrase.
4. **U15:** soften "characterised" for Leong et al.
5. **U20:** add a verified Birnbaum reference, or drop "Birnbaum-type".
6. **U24:** report Cliff's δ or remove it and its citation.
7. **U25:** add the coupon-collector comparison, or soften "consistent with".
8. **U26:** "This work / Correlated errors" Y → P; note the abstract-only reading of Kong in the caption.
9. **Minor (U1, U2, U6, U7, U12, U13, U19, U22, U23):** wording precision as listed above.
