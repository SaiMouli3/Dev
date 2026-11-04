# Supplementary Material S1: Annotated Literature Review

Compiled 2026-11-02 for the D3RA-DNA manuscript, which covers self-calibrating, adaptive (closed-loop, Bayesian) allocation and
stopping of sequencing reads when retrieving erasure-coded data from several heterogeneous DNA storage containers ("pools").
The setting includes container failures and correlated, strand-specific sequencing errors.

**Scope and provenance.**
- Every bibliographic record below appears in `manuscript/references.bib`. Each was checked against Crossref, the arXiv API,
  DataCite/ENA, or the publisher or proceedings page; see `audit/citation_verification.md`.
- Content summaries come from the project's earlier reading notes (`docs/literature_review.md`,
  `docs/novelty_assessment.md`, `docs/channel_calibration.md`) and from the abstracts retrieved on 2026-11-02.
- The depth of reading is stated wherever it matters. "Abstract only" means nothing beyond the abstract was read, and
  statements of the form "X does not model Y" are then limited to what the abstract shows.
- Classic works (e.g. Reed–Solomon, Thompson) are summarized from their well-known content and were not re-read.
- Two works are cited in the project notes but are **excluded** from the bibliography, because no primary record could be
  reached: Birnbaum (1969) on component importance, and Holm (1979) on step-down multiple testing. See §12.

Bib keys are given in `[brackets]`.

---

## 1. DNA data storage systems

### [church2012next] Next-generation digital information storage in DNA
- **Authors / year / venue:** G. M. Church, Y. Gao, S. Kosuri; 2012; *Science* 337(6102):1628.
- **DOI:** 10.1126/science.1226355
- **Objective:** Demonstrate that DNA can store megabit-scale digital data.
- **Method:** One bit per base, encoded in addressed oligonucleotides.
- **Setup:** Array synthesis followed by Illumina sequencing.
- **Key findings:** Recovered a book-sized file with a small number of residual errors.
- **Limitations:** No strand-level error-correcting code. Errors were corrected through redundant copies only.
- **Relation to this work:** Historical baseline. It uses a single pool and no erasure code across strands.

### [goldman2013towards] Towards practical, high-capacity, low-maintenance information storage in synthesized DNA
- **Authors / year / venue:** N. Goldman, P. Bertone, S. Chen, C. Dessimoz, E. M. LeProust, B. Sipos, E. Birney; 2013;
  *Nature* 494(7435):77–80.
- **DOI:** 10.1038/nature11875
- **Objective:** Store data robustly at scale.
- **Method:** A rotating base-3 code that avoids homopolymers, with 4-fold overlapping segments as redundancy.
- **Setup:** Strands of 117 nt and Illumina reads. The reads are in ENA PRJEB1186.
- **Key findings:** Successful recovery of several files.
- **Limitations:** Repetition-style redundancy rather than an MDS outer code. Single pool.
- **Relation to this work:** An early example of how much redundancy a single pool needs.

### [grass2015robust] Robust chemical preservation of digital information on DNA in silica with error-correcting codes
- **Authors / year / venue:** R. N. Grass, R. Heckel, M. Puddu, D. Paunescu, W. J. Stark; 2015;
  *Angew. Chem. Int. Ed.* 54(8):2552–2555.
- **DOI:** 10.1002/anie.201411378
- **Objective:** Long-term, error-tolerant storage.
- **Method:** Concatenated inner and outer Reed–Solomon (RS) codes, with silica encapsulation.
- **Setup:** 83 kB stored in 4,991 strands of 158 nt, with accelerated aging at 70 °C.
- **Key findings:** Exact recovery after the aging that was equivalent to long-term storage.
- **Limitations:** Single container. Strand loss is treated as erasures within one pool.
- **Relation to this work:** The inner/outer RS architecture we adopt. Encapsulation also gives a physical notion of a
  container as a failure domain.

### [bornholt2016dna] A DNA-based archival storage system
- **Authors / year / venue:** J. Bornholt, R. Lopez, D. M. Carmean, L. Ceze, G. Seelig, K. Strauss; 2016;
  *ASPLOS '16*, pp. 637–649.
- **DOI:** 10.1145/2872362.2872397
- **Objective:** A systems architecture for DNA storage.
- **Method:** Key-value store with random access through PCR primers, plus redundancy schemes.
- **Setup:** Small wet-lab experiments together with simulation.
- **Key findings:** Random access is feasible, and the paper quantifies the trade-off between redundancy and density.
- **Limitations:** No reasoning about read allocation across multiple containers.
- **Relation to this work:** A systems view that motivates treating a pool as a storage device.

### [erlich2017fountain] DNA Fountain enables a robust and efficient storage architecture
- **Authors / year / venue:** Y. Erlich, D. Zielinski; 2017; *Science* 355(6328):950–954.
- **DOI:** 10.1126/science.aaj2038
- **Objective:** Near-capacity coding for DNA storage.
- **Method:** An LT/fountain outer code. Each droplet carries RS bytes, and droplets are screened for GC content and
  homopolymers.
- **Setup:** 2.14 MB in 72,000 oligos of 152 nt, with Illumina reads (ENA PRJEB19305 and PRJEB19307).
- **Key findings:** High information density. Droplets that fail the inner RS check are rejected before LT decoding.
- **Limitations:** Single pool. The fountain code needs an overhead ε.
- **Relation to this work:** A precedent for verification-gated fragment acceptance, and an alternative to an MDS outer code.

### [organick2018random] Random access in large-scale DNA data storage
- **Authors / year / venue:** L. Organick et al. (22 authors); 2018; *Nat. Biotechnol.* 36(3):242–248.
- **DOI:** 10.1038/nbt.4079
- **Objective:** Random access at scale.
- **Method:** File-specific PCR primers and an RS outer code with about 15% redundancy.
- **Setup:** 35 files totalling more than 200 MB in more than 13 M oligos. Illumina and nanopore readout.
- **Key findings:**
  - Decoding at a median coverage of about 5 reads per sequence with Illumina.
  - Nanopore readout needed higher coverage.
- **Limitations:**
  - Redundancy is per file.
  - No failures across pools or containers are modelled.
  - Sequencing depth is fixed in advance.
- **Relation to this work:** Real evidence of the coverage needed for decoding, and of physically partitioned storage.

### [antkowiak2020low] Low cost DNA data storage using photolithographic synthesis and advanced information reconstruction and error correction
- **Authors / year / venue:** P. L. Antkowiak et al.; 2020; *Nat. Commun.* 11:5345.
- **DOI:** 10.1038/s41467-020-19148-3
- **Objective:** Cheap, high-error synthesis.
- **Method:** Light-directed synthesis of 60-nt strands, with stronger reconstruction and error correction.
- **Key findings:** Noisy, cheap writing can be compensated at decode time.
- **Limitations:** Single pool, with high synthesis error rates.
- **Relation to this work:** Supports making the read side adaptive to the quality of each container.

### [lopez2019dna] DNA assembly for nanopore data storage readout
- **Authors / year / venue:** R. Lopez et al.; 2019; *Nat. Commun.* 10:2933.
- **DOI:** 10.1038/s41467-019-10978-4
- **Objective:** Make nanopore readout of short-oligo storage efficient.
- **Method:** Assemble oligos into long molecules for nanopore sequencing.
- **Key findings:**
  - Real-time sequencing makes it possible to "sequence until sufficient coverage has been acquired for successful
    decoding" (full text read).
  - Fig. 1c of the paper shows a sequence-until-decoded workflow.
- **Limitations:** The full text does not show whether stopping was automated and triggered by the decoder. There is no
  allocation across containers.
- **Relation to this work:** **Prior art for decode-driven stopping.** Our stopping rule is not presented as new.

## 2. Channel characterization and molecular bias

### [heckel2019characterization] A characterization of the DNA data storage channel
- **Authors / year / venue:** R. Heckel, G. Mikutis, R. N. Grass; 2019; *Sci. Rep.* 9:9663.
- **DOI:** 10.1038/s41598-019-45832-6
- **Objective:** Characterize the errors and coverage of real DNA storage experiments.
- **Method:** Reanalysis of datasets from ETH, Goldman and Erlich.
- **Key findings:**
  - Errors within molecules come mostly from synthesis and sequencing.
  - Strand loss comes from handling and storage.
  - Coverage is approximately negative binomial.
- **Limitations:** Single-pool datasets. No analysis across containers.
- **Relation to this work:** Justifies our Gamma-Poisson (negative binomial) coverage model.

### [chen2020quantifying] Quantifying molecular bias in DNA data storage
- **Authors / year / venue:** Y.-J. Chen et al.; 2020; *Nat. Commun.* 11:3264.
- **DOI:** 10.1038/s41467-020-16958-3
- **Objective:** Separate synthesis bias from PCR bias.
- **Method:** UMIs and controlled amplification (full text read).
- **Key findings:**
  - Synthesis bias is significant and related to the oligo's spatial location on the chip.
  - PCR bias is dominated by stochasticity and is not correlated with GC content.
- **Limitations:** Does not decompose variance across replicates.
- **Relation to this work:** Explains why aliquots of the same synthesis batch share per-strand abundance. We report a
  "batch-inherited" component, not a "sequence-intrinsic" one.

### [gimpel2023digital] A digital twin for DNA data storage based on comprehensive quantification of errors and biases
- **Authors / year / venue:** A. L. Gimpel, W. J. Stark, R. Heckel, R. N. Grass; 2023; *Nat. Commun.* 14:6026.
- **DOI:** 10.1038/s41467-023-41729-1
- **Objective:** A calibrated, process-level simulator (DT4DDS).
- **Method:** Quantifies errors per process step across 40 datasets and fits a lognormal coverage distribution.
- **Setup:** Experiments covering synthesis, PCR cycling and aging. The reads are ENA PRJEB65931 [gimpel2023prjeb65931].
- **Key findings:**
  - Error rates are given for each process step.
  - The spread of PCR efficiency is narrow and does not depend on GC.
  - Aging did not significantly change the coverage distribution.
- **Limitations:**
  - Correlations with sequence properties were not investigated, as the paper itself states.
  - The main text has no decomposition of variance shared across conditions.
- **Relation to this work:** **Primary real-data source.** We treat six handling conditions as same-batch aliquot
  "containers", and the paper's caveats about lineage carry over to our results.

## 3. Coding and information theory

### [reed1960polynomial] Polynomial codes over certain finite fields
- **Authors / year / venue:** I. S. Reed, G. Solomon; 1960; *J. SIAM* 8(2):300–304.
- **DOI:** 10.1137/0108018
- **Content:** Introduces Reed–Solomon (MDS) codes.
- **Relation to this work:** Used as both the inner and the outer (stripe) code.

### [lenz2020coding] Coding over sets for DNA storage
- **Authors / year / venue:** A. Lenz, P. H. Siegel, A. Wachter-Zeh, E. Yaakobi; 2020; *IEEE Trans. Inf. Theory*
  66(4):2331–2351.
- **DOI:** 10.1109/TIT.2019.2961265
- **Content:** Bounds (Gilbert–Varshamov and sphere-packing) and near-optimal constructions for codes over unordered sets of
  strands.
- **Limitations:** Static, single-pool channel.
- **Relation to this work:** Theoretical background for coding over unordered strands.

### [press2020hedges] HEDGES error-correcting code for DNA storage corrects indels and allows sequence constraints
- **Authors / year / venue:** W. H. Press, J. A. Hawkins, S. K. Jones Jr., J. M. Schaub, I. J. Finkelstein; 2020;
  *PNAS* 117(31):18489–18496.
- **DOI:** 10.1073/pnas.2004821117
- **Content:**
  - A hash-based inner code that corrects indels, combined with an RS outer code applied diagonally across strands.
  - Error-free recovery at about 10% nucleotide error.
- **Relation to this work:** An alternative inner code. Our pipeline uses RS with a CRC instead.

### [shomorony2021dna] DNA-based storage: Models and fundamental limits
- **Authors / year / venue:** I. Shomorony, R. Heckel; 2021; *IEEE Trans. Inf. Theory* 67(6):3675–3689.
- **DOI:** 10.1109/TIT.2021.3058966
- **Content:** Defines the noisy shuffling-sampling channel and shows that index-based coding is optimal in many regimes.
- **Limitations:** Single pool, i.i.d. sampling.
- **Relation to this work:** The channel abstraction that our per-container model extends.

### [weinberger2022dna] The DNA storage channel: Capacity and error probability bounds
- **Authors / year / venue:** N. Weinberger, N. Merhav; 2022; *IEEE Trans. Inf. Theory* 68(9):5657–5700.
- **DOI:** 10.1109/TIT.2022.3176371
- **Content:**
  - Capacity and error-exponent bounds when sampling M molecules N times with replacement.
  - The dominant error event is a sampling "outage."
- **Relation to this work:** Coverage depth drives reliability. This motivates spending reads where outage risk is highest.

## 4. Coverage depth and multi-container DNA storage (closest line of work)

### [barlev2025cover] Cover your bases: How to minimize the sequencing coverage in DNA storage systems
- **Authors / year / venue:** D. Bar-Lev, O. Sabary, R. Gabrys, E. Yaakobi; 2025; *IEEE Trans. Inf. Theory*
  71(1):192–218. ISIT 2023 version; arXiv:2305.05656.
- **DOI:** 10.1109/TIT.2024.3496587
- **Objective:** Minimize the expected number of reads until decoding succeeds (the "coverage depth problem").
- **Method:** Coupon-collector analysis.
- **Key findings:**
  - Under noiseless uniform sampling, MDS codes are optimal.
  - For an [n,k] MDS code, E[reads] = n(H_n − H_{n−k}).
  - The uniform sampling distribution minimizes expected reads.
  - The paper also gives bounds on the distribution of the read count, and a random-access variant.
- **Limitations:** A single pool with a fixed, non-adaptive sampling law. No container failures.
- **Relation to this work:**
  - The uniform-allocation baseline.
  - The theory behind our post-failure margin check: with R = K+m surviving fragments, E = R(H_R − H_m).

### [gruica2024combinatorial] A combinatorial perspective on random access efficiency for DNA storage
- **Authors / year / venue:** A. Gruica, D. Bar-Lev, A. Ravagnani, E. Yaakobi; 2024; *Proc. IEEE ISIT*, pp. 675–680.
  arXiv:2401.15722.
- **DOI:** 10.1109/ISIT57864.2024.10619151
- **Content (abstract):**
  - Studies T_max(G), the worst-case expected number of reads to recover any single information strand.
  - Gives formulas based on column dependencies and on recovery sets, and introduces "recovery balanced" codes.
- **Limitations:** Uniform sampling, one pool, no failures.
- **Relation to this work:** Random-access counterpart of the full-file coverage depth problem.

### [cao2025optimizing] Optimizing sequencing coverage depth in DNA storage: Insights from DNA storage data
- **Authors / year / venue:** R. Cao, X. Chen; 2025; *Proc. IEEE ISIT*. arXiv:2501.06801.
- **DOI:** 10.1109/ISIT63088.2025.11195510
- **Content (abstract):**
  - Extends coverage depth from uniform to log-normal sampling, fitted to the authors' own PCR and sequencing data.
  - The expected minimum coverage depth increases with MDS code rate.
  - Lower bounds are extended to noisy log-normal channels.
- **Limitations:** Single pool. Allocation is not adaptive.
- **Relation to this work:** The closest work that uses non-uniform, data-fitted sampling. We add a container dimension
  and online learning.

### [kashanna2026reliability] On the reliability of information retrieval from MDS coded data in DNA storage
- **Authors / year / venue:** S. Kas Hanna; 2026; *IEEE Trans. Mol. Biol. Multi-Scale Commun.* 12:585–598.
  ISIT 2025 version; arXiv:2502.06618.
- **DOI:** 10.1109/TMBMC.2026.3690724
- **Content (abstract):**
  - The probability of successful retrieval with inner and outer MDS codes under i.i.d. substitutions.
  - It is a function of the total read count, how reads are spread across strands, and the code rates.
  - Gives the minimum number of reads and the best split of rate between inner and outer codes.
- **Limitations:** i.i.d. errors, one pool, a fixed read budget.
- **Relation to this work:** A model-based read budget. Our results show that correlated errors make i.i.d. budgets
  optimistic.

### [levy2026expected] Expected recovery time in DNA-based distributed storage systems
- **Authors / year / venue:** A. Levy, R. Con, E. Yaakobi, H. M. Kiah; 2026; arXiv:2602.07601 (preprint).
- **Objective:** Robustness against **container failures** in DNA storage spread over M containers (full text read
  except the last ~7%).
- **Method:**
  - Each read returns a uniformly random strand of a container.
  - Uses row-wise erasure codes across containers.
  - Analyses generalized coupon-collector problems.
- **Key findings:**
  - Expected time to recover a failed container's contents from the survivors.
  - The leading-order cost scales as 1/(surplus+1).
- **Limitations (from its model section):** No heterogeneous container efficiency, abundance bias, sequencing errors,
  adaptive allocation or real data.
- **Relation to this work:** **Closest prior work on container failures.** Any claim of being first to model container
  failure would be false.

### [kong2026coverage] The coverage depth problem in distributed DNA data storage
- **Authors / year / venue:** X. Kong, O. Elishco, C. Wang, T. M. Duman; 2026; arXiv:2610.02931 (preprint, 33 pp.).
- **Objective:** Coverage depth when strands are partitioned among M containers (abstract only).
- **Method:** Each round samples one strand uniformly, with replacement, from each container. Noiseless channel.
- **Key findings:**
  - Exact recovery-time distribution and expectation for arbitrary linear codes and partitions.
  - MDS codes are optimal for every fixed partition.
  - A universal lower bound on total read cost.
  - Container-size regimes in which distribution saves reads, versus regimes where it only adds parallelism.
  - A partition that attains the minimum total read cost.
- **Limitations (as far as the abstract shows):**
  - No container failures.
  - No heterogeneity beyond container size.
  - No biased sampling.
  - No adaptive allocation: the per-round sampling and the partition are fixed in advance.
- **Relation to this work:** **Closest prior work on multi-container coverage.** It is the source of the static
  best-partition baseline, and its formulas validate our simulator in the uniform case without failures.

## 5. Trace reconstruction, clustering and consensus

### [batu2004reconstructing] Reconstructing strings from random traces
- **Authors / year / venue:** T. Batu, S. Kannan, S. Khanna, A. McGregor; 2004; *SODA*, pp. 910–918.
- **Content:**
  - Introduces bitwise majority alignment (BMA) for reconstruction from deletion channels.
  - O(log n) traces suffice for random strings at deletion probability O(1/log n).
- **Relation to this work:** BMA is the consensus algorithm in our nanopore analysis.

### [rashtchian2017clustering] Clustering billions of reads for DNA data storage
- **Authors / year / venue:** C. Rashtchian et al.; 2017; *NeurIPS (NIPS) 30*.
- **Content:** Distributed approximate clustering by edit distance. The abstract reports higher accuracy and a 1000×
  speedup on three real datasets.
- **Limitations:** Relies on clusters being well separated. Clustering errors can masquerade as channel errors.
- **Relation to this work:** Clustering is part of the decoding stack, and residual clustering errors are one possible
  source of the strand-specific error correlation we observe.

### [srinivasavaradhan2021trellis] Trellis BMA: Coded trace reconstruction on IDS channels for DNA storage
- **Authors / year / venue:** S. R. Srinivasavaradhan, S. Gopi, H. D. Pfister, S. Yekhanin; 2021; *Proc. IEEE ISIT*,
  pp. 2453–2458. Extended version arXiv:2107.06440.
- **DOI:** 10.1109/ISIT45174.2021.9517821
- **Content:** Trellis-based coded trace reconstruction for insertion/deletion/substitution (IDS) channels. Releases a real
  clustered nanopore read dataset.
- **Limitations:** The arXiv note added in August 2024 states that the dataset's reference sequences have long-range
  dependencies because of a generation error, so some clusters may be malformed.
- **Relation to this work:** Data source for the nanopore consensus curves, through [barlev2024binned]. The caveat is
  why we re-checked the result on the DNAformer flowcells.

### [sabary2024reconstruction] Reconstruction algorithms for DNA-storage systems
- **Authors / year / venue:** O. Sabary, A. Yucovich, G. Shapira, E. Yaakobi; 2024; *Sci. Rep.* 14:1951.
- **DOI:** 10.1038/s41598-024-51730-3
- **Content:** Trace-reconstruction algorithms for DNA storage, evaluated at high error probabilities (up to about 0.27).
- **Relation to this work:** An alternative consensus stage, and a benchmark context.

### [barlev2025deep] Scalable and robust DNA-based storage via coding theory and deep learning (DNAformer)
- **Authors / year / venue:** D. Bar-Lev, I. Orr, O. Sabary, T. Etzion, E. Yaakobi; 2025; *Nat. Mach. Intell.*
  7(4):639–649.
- **DOI:** 10.1038/s42256-025-01003-z
- **Content (arXiv v3 full text read):**
  - A transformer reconstruction model trained on simulated data, combined with tensor-product codes.
  - Uses a confidence filter and a "safety margin."
  - Demonstrated on 3.1 MB with Illumina and nanopore readout.
- **Limitations:** Single pool. Does not quantify the gap between i.i.d.-simulated and real nanopore data.
- **Relation to this work:** The source of the two-flowcell nanopore dataset [technion2025dnaformerreads] used for our
  test of correlated errors across flowcells.

## 6. Adaptive and selective nanopore sequencing

### [loose2016real] Real-time selective sequencing using nanopore technology
- **Authors / year / venue:** M. Loose, S. Malla, M. Stout; 2016; *Nat. Methods* 13(9):751–754.
- **DOI:** 10.1038/nmeth.3930
- **Content:** "Read Until": molecules are rejected in real time on the basis of their signal.
- **Relation to this work:** The physical actuator that makes per-container allocation on a flow cell feasible.

### [payne2021readfish] Readfish enables targeted nanopore sequencing of gigabase-sized genomes
- **Authors / year / venue:** A. Payne, N. Holmes, T. Clarke, R. Munro, B. J. Debebe, M. Loose; 2021; *Nat. Biotechnol.*
  39(4):442–450.
- **DOI:** 10.1038/s41587-020-00746-x
- **Content:** A software toolkit for targeted adaptive sampling on genomic targets.
- **Relation to this work:** Practical adaptive-sampling software. Per-barcode switch-off is described in a later readfish
  preprint, which we did not cite because we did not verify it.

### [weilguny2023dynamic] Dynamic, adaptive sampling during nanopore sequencing using Bayesian experimental design (BOSS-RUNS)
- **Authors / year / venue:** L. Weilguny, N. De Maio, R. Munro, C. Manser, E. Birney, M. Loose, N. Goldman; 2023;
  *Nat. Biotechnol.* 41(7):1018–1025.
- **DOI:** 10.1038/s41587-022-01580-z
- **Content (abstract level):** Computes the expected benefit of each new read and redistributes coverage in real time.
- **Limitations:** Genomics, not erasure-coded storage. The objective is variant-calling information, not K-of-N
  decodability.
- **Relation to this work:** **Closest closed-loop Bayesian analogue** of our allocation, in another domain.

## 7. Distributed-storage placement and correlated failures

### [weil2006crush] CRUSH: Controlled, scalable, decentralized placement of replicated data
- **Authors / year / venue:** S. Weil, S. Brandt, E. Miller, C. Maltzahn (names as in Crossref); 2006; *SC '06*, art. 31.
- **DOI:** 10.1109/SC.2006.19
- **Content:** Pseudo-random placement that respects failure domains (host, rack, row).
- **Relation to this work:** The per-domain cap on chunks is standard practice, and we cite it rather than claim it.

### [nath2006subtleties] Subtleties in tolerating correlated failures in wide-area storage systems
- **Authors / year / venue:** S. Nath, H. Yu, P. B. Gibbons, S. Seshan; 2006; *NSDI '06* (USENIX).
- **Content (abstract and summaries):**
  - Under real correlated-failure traces, extra redundancy has diminishing returns.
  - Design rankings can flip between independent and correlated failures.
- **Relation to this work:** Spreading across containers can lose to concentrating in one when failures are correlated.

### [ford2010availability] Availability in globally distributed storage systems
- **Authors / year / venue:** D. Ford et al.; 2010; *OSDI '10* (USENIX).
- **Content (read in part):**
  - Availability measured over a year of production at Google.
  - Rack-aware placement increases stripe MTTF under correlated failure bursts.
- **Relation to this work:** Empirical support for placement that accounts for failure domains.

### [leong2012distributed] Distributed storage allocations
- **Authors / year / venue:** D. Leong, A. G. Dimakis, T. Ho; 2012; *IEEE Trans. Inf. Theory* 58(7):4733–4752.
- **DOI:** 10.1109/TIT.2012.2191135
- **Content (model and theorems read):**
  - Allocates a coded object over nodes that are accessed with probability p.
  - A symmetric allocation over m nodes succeeds with probability P[Bin(m,p) ≥ ⌈m/T⌉].
  - Minimal spreading versus maximal spreading regimes.
- **Limitations:** A common p, no read cost.
- **Relation to this work:** Our spreading threshold is **a special case** of this result, and is presented that way.

### [cidon2013copysets] Copysets: Reducing the frequency of data loss in cloud storage
- **Authors / year / venue:** A. Cidon, S. Rumble, R. Stutsman, S. Katti, J. Ousterhout, M. Rosenblum; 2013;
  *USENIX ATC '13*, pp. 37–48.
- **Content:** Under correlated failures, fewer distinct placement groups reduce the frequency of data loss.
- **Relation to this work:** The combinatorial trade-off between concentrating and spreading.

## 8. Sequential decisions and bandits

### [thompson1933likelihood] On the likelihood that one unknown probability exceeds another in view of the evidence of two samples
- **Authors / year / venue:** W. R. Thompson; 1933; *Biometrika* 25(3–4):285–294.
- **DOI:** 10.1093/biomet/25.3-4.285
- **Content:** Posterior-probability-matching allocation, now known as Thompson sampling.
- **Relation to this work:** One of our allocation policies samples from the Gamma posterior of each container's
  efficiency.

### [bubeck2013optimal] Optimal discovery with probabilistic expert advice: finite time analysis and macroscopic optimality
- **Authors / year / venue:** S. Bubeck, D. Ernst, A. Garivier; 2013; *JMLR* 14(17):601–623.
- **Content (snippet level):** Good-UCB chooses which of several sampling sources to query in order to maximize discovery
  of new items, with an upper confidence bound on the Good–Turing missing mass.
- **Relation to this work:** **The closest general analogue of choosing which container to sequence.** Its objective is
  unweighted discovery, whereas ours is completing every stripe to K of N.

## 9. Nanopore systematic errors

### [wick2019performance] Performance of neural network basecalling tools for Oxford Nanopore sequencing
- **Authors / year / venue:** R. R. Wick, L. M. Judd, K. E. Holt; 2019; *Genome Biol.* 20:129.
- **DOI:** 10.1186/s13059-019-1727-y
- **Content:** Consensus accuracy is about 99.5%, which shows that systematic basecalling errors remain.
- **Relation to this work:** Systematic errors are not averaged out by extra reads. This is prior evidence for our
  finding of strand-specific correlated errors.

### [delahaye2021sequencing] Sequencing DNA with nanopores: Troubles and biases
- **Authors / year / venue:** C. Delahaye, J. Nicolas; 2021; *PLoS One* 16(10):e0257521.
- **DOI:** 10.1371/journal.pone.0257521
- **Content (abstract level):**
  - Homopolymers and short repeats account for about half of the errors.
  - Errors depend on GC content.
  - Quality scores are miscalibrated.
- **Relation to this work:** Context-dependent errors. Our context-only model is the baseline that the hierarchical
  strand-level model improves on.

## 10. Statistics

### [efron1993introduction] An Introduction to the Bootstrap
- **Authors / year / venue:** B. Efron, R. J. Tibshirani; 1993; book. The Crossref record lists Springer US as publisher.
- **DOI:** 10.1007/978-1-4899-4541-9
- **Relation to this work:** Bootstrap confidence intervals for read savings.

### [cliff1993dominance] Dominance statistics: Ordinal analyses to answer ordinal questions
- **Authors / year / venue:** N. Cliff; 1993; *Psychol. Bull.* 114(3):494–509.
- **DOI:** 10.1037/0033-2909.114.3.494
- **Relation to this work:** Cliff's δ as the effect size for comparing policies.

## 11. Public datasets used

| Key | Record | Content | Used for |
|---|---|---|---|
| [gimpel2023prjeb65931] | ENA PRJEB65931, "DT4DDS – Characterization of Errors and Biases in DNA Data Storage" (FML, first public 2023-09-08) | Illumina runs of DT4DDS experiments | Six same-batch "containers", emulation of placement and stopping |
| [technion2025dnaformerreads] | Zenodo 10.5281/zenodo.17473983 (2025) | Binned DNAformer nanopore and Illumina reads | Test of correlated errors across two flowcells |
| [barlev2024binned] | Zenodo 10.5281/zenodo.14296588 (2024, Bar-Lev et al.) | Binned benchmark datasets, including the clustered nanopore reads of Srinivasavaradhan et al. | Real nanopore consensus curves |

Note: earlier project documents attribute Zenodo 14296588 to "Srinivasavaradhan et al." The record is in fact a
re-packaging by Bar-Lev, Orr, Sabary, Etzion and Yaakobi that *contains* that dataset. The paper should cite both
[barlev2024binned] and [srinivasavaradhan2021trellis].

## 12. Works excluded from the bibliography

| Work | Reason | Suggested handling |
|---|---|---|
| Z. W. Birnbaum, "On the importance of different components in a multicomponent system," in *Multivariate Analysis II* (P. R. Krishnaiah, ed.), Academic Press, 1969, pp. 581–592 | No Crossref or publisher record was reachable; the details are consistent across secondary citations only | Add after checking a library copy. The pivotal-fragment gradient must still be attributed to Birnbaum importance in the text |
| S. Holm, "A simple sequentially rejective multiple test procedure," *Scand. J. Stat.* 6(2):65–70, 1979 (JSTOR 4615733) | DOI 10.2307/4615733 is not in Crossref, and the JSTOR page could not be fetched | Add after checking JSTOR manually |
| Sokolovskii et al., bioRxiv 10.1101/2024.11.05.622081; Pic et al., arXiv:2603.05220; readfish barcode preprint; Zhang et al., *Sci. Adv.* 2026; Abraham et al. ISIT 2026; Hamoum et al. | Not re-verified in this pass (they are outside the requested scope, or the author list is incomplete) | Verify before citing. Zhang et al. 2026 is relevant prior art for real-time decode-until-success |

---

## 13. Comparison matrix (closest works)

Legend:
- **Y** = the paper does this.
- **N** = it does not, according to the paper's own model or setup.
- **P** = partially (see note).
- **n/s** = not stated, or not determinable from what was read.
- "Implicit" stopping means that the analysed quantity is the number of reads until decoding succeeds, with no
  operational stopping rule in a closed-loop system.

| Work | Multi-container | Container failures | Heterogeneous containers | Correlated errors modelled | Adaptive read allocation | Stopping rule | Real-data evaluation |
|---|---|---|---|---|---|---|---|
| Organick 2018 [organick2018random] | P (oligos from several synthesis pools, per project notes; redundancy per file) | N | n/s | N | N | n/s (fixed coverage) | Y (Illumina, nanopore) |
| Lopez 2019 [lopez2019dna] | N | N | N | n/s | N | P (sequence until decodable; automation n/s) | Y (nanopore) |
| Gimpel 2023 [gimpel2023digital] | N (one pool, many conditions) | N (strand dropout only) | P (characterizes conditions) | N | N | N | Y (40 datasets) |
| Bar-Lev 2025 [barlev2025cover] | N | N | N | N | N (uniform shown optimal for uncontrolled sampling) | Implicit | N |
| Gruica 2024 [gruica2024combinatorial] | N | N | N | N | N | Implicit | N |
| Cao 2025 [cao2025optimizing] | N | N | N (log-normal per strand) | N | N | Implicit | P (log-normal fitted to own data) |
| Kas Hanna 2026 [kashanna2026reliability] | N | N | N | N (i.i.d. substitutions) | N | N (fixed read count) | N |
| Levy 2026 [levy2026expected] | Y | Y | N | N | N (uniform per container) | Implicit (expected recovery time) | N |
| Kong 2026 [kong2026coverage] | Y | N (per abstract) | P (container sizes only) | N | N (static partition) | Implicit | N (per abstract) |
| Leong 2012 [leong2012distributed] | Y (storage nodes) | Y (node access prob. p) | N (common p) | N | N | N | N |
| Weilguny 2023 [weilguny2023dynamic] | N (genomic loci, not containers) | N | n/s | N | Y (Bayesian, real time) | n/s | Y (genomics) |
| Bubeck 2013 [bubeck2013optimal] | Y (sampling "experts") | N | Y (distributions differ) | N | Y (Good-UCB) | N (fixed horizon) | N (simulation) |
| **This work** | Y | Y (injected on real read streams; real failure rates not measured) | Y (measured across real runs) | Y (strand-specific, hierarchical beta-binomial) | Y (Gamma-Poisson posterior) | Y (verification-gated) | Y (PRJEB65931; Zenodo 17473983, 14296588) |

Caveats for the "This work" row:
- The real-data containers are aliquots of **one synthesis batch**, and they are lineage-related.
- Container failures are injected, not observed.
- The allocation component adds ≤ 10% on real reads, and pivotal scoring ties coverage-greedy scoring. Most of the
  saving comes from stopping and placement.

## 14. Research gap

**Prior work in this area.**
- *Coverage depth.* The coverage-depth literature [barlev2025cover, gruica2024combinatorial, cao2025optimizing,
  kashanna2026reliability] characterizes how many reads are needed to decode under a sampling law that is fixed in
  advance, mostly uniform. Bar-Lev et al. show that MDS codes and uniform sampling are optimal when the sampler cannot
  be controlled.
- *Multi-container storage.* Two recent preprints extend coverage depth to several containers:
  - Kong et al. [kong2026coverage] derive exact recovery-time laws for static partitions under uniform sampling.
  - Levy et al. [levy2026expected] analyse the expected time to recover a failed container.
  - **Both multi-container DNA storage and container failure are therefore already studied.**
- *Placement under failure domains.* Placement in conventional distributed storage is classical
  [leong2012distributed, weil2006crush, ford2010availability, cidon2013copysets, nath2006subtleties]. Our homogeneous
  spreading threshold is a special case of Leong, Dimakis and Ho.
- *Stopping.* Sequencing until the data decode is known in practice [lopez2019dna].
- *Closed-loop allocation in other domains.* Closed-loop Bayesian allocation of reads exists in genomics
  [weilguny2023dynamic, payne2021readfish, loose2016real] and in general adaptive discovery [bubeck2013optimal].
- *Correlated errors.* Systematic, context-dependent nanopore errors are documented
  [wick2019performance, delahaye2021sequencing].

**What the searches did not find.** This is "not found", not proof of absence. We found no work that does all of the
following at once:
1. allocates sequencing reads *online* across several DNA storage containers whose efficiency is unknown and
   heterogeneous;
2. learns the efficiency and failure state of each container from the read stream;
3. directs reads toward the erasure-coded stripes that are furthest from decodability;
4. stops on a verified K-of-N criterion;
5. evaluates the policy on real sequencing read streams with injected container loss and with errors correlated across
   reads of the same strand.

**The contribution we claim.**
- It is an integration and a real-data evaluation, not a new primitive.
- We quantify how much of the read saving comes from placement, from stopping and from allocation.
- The largest gains come from the components that are already known (stopping and placement).
- The new component (allocation) gives a smaller, but measurable, benefit on real data.
- Independent-error simulators overstate the success of nanopore consensus by 6–15 percentage points on verified
  references from one dataset.
