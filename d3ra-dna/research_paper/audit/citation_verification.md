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
