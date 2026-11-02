# D3RA-DNA: Literature Review and Novelty Positioning

Compiled 2026-11-01. I checked every bibliographic entry below against Crossref, the arXiv API, PubMed Central or Europe PMC, or the
publisher or proceedings page, unless it is marked **[unverified]**. Quantitative claims come from the cited full text unless
marked otherwise. "Preprint" means I found no peer-reviewed version.

---

## 1. Foundational DNA data storage systems

| # | Citation | Key points for D3RA-DNA |
|---|----------|--------------------------|
| [1] | G. M. Church, Y. Gao, S. Kosuri, "Next-generation digital information storage in DNA," *Science* 337(6102):1628, 2012. doi:10.1126/science.1226355 | First megabyte-scale demonstration. Uses one bit per base and no strand-level ECC. |
| [2] | N. Goldman, P. Bertone, S. Chen, C. Dessimoz, E. M. LeProust, B. Sipos, E. Birney, "Towards practical, high-capacity, low-maintenance information storage in synthesized DNA," *Nature* 494:77–80, 2013. doi:10.1038/nature11875 | Rotating ternary (base-3) code that avoids homopolymers. Uses 4-fold overlapping-segment redundancy. Strands are 117 nt (per [7]). Reads are in ENA PRJEB1186. |
| [3] | R. N. Grass, R. Heckel, M. Puddu, D. Paunescu, W. J. Stark, "Robust chemical preservation of digital information on DNA in silica with error-correcting codes," *Angew. Chem. Int. Ed.* 54(8):2552–2555, 2015. doi:10.1002/anie.201411378 | First inner/outer Reed–Solomon (RS) concatenated code. Stores 83 kB in 4,991 strands of 158 nt. Silica encapsulation, with accelerated aging at 70 °C for 1 week. |
| [4] | J. Bornholt, R. Lopez, D. M. Carmean, L. Ceze, G. Seelig, K. Strauss, "A DNA-based archival storage system," *ASPLOS '16*, pp. 637–649, 2016. doi:10.1145/2872362.2872397 | Systems view (Microsoft/UW): key-value store with PCR-primer random access. |
| [5] | Y. Erlich, D. Zielinski, "DNA Fountain enables a robust and efficient storage architecture," *Science* 355(6328):950–954, 2017. doi:10.1126/science.aaj2038 | LT/fountain outer code with per-droplet RS. Stores 2.14 MB in 72,000 oligos of 152 nt with about 7% redundancy. Screening: GC 45–55%, homopolymer ≤3. Reads are in ENA PRJEB19305 and PRJEB19307. |
| [6] | L. Organick, S. D. Ang, Y.-J. Chen, R. Lopez, S. Yekhanin, K. Makarychev, M. Z. Racz, G. Kamath, P. Gopalan, B. Nguyen, C. N. Takahashi, S. Newman, et al., "Random access in large-scale DNA data storage," *Nat. Biotechnol.* 36:242–248, 2018. doi:10.1038/nbt.4079 | Stores 35 files (>200 MB) in more than 13 M oligos across 9 synthesis pools, with PCR random access. 15% RS outer redundancy. Decoded at a median coverage of 5 reads per sequence (range 4–14 per file). Nanopore decoding needed 36× and 80×. |
| [7] | R. Heckel, G. Mikutis, R. N. Grass, "A characterization of the DNA data storage channel," *Sci. Rep.* 9:9663, 2019. doi:10.1038/s41598-019-45832-6 (arXiv:1803.03322) | Reanalyses datasets from ETH, Goldman and Erlich. Errors within molecules come mostly from synthesis and sequencing. Handling and storage cause strand loss. Coverage is approximately negative binomial. |
| [8] | Y.-J. Chen, C. N. Takahashi, L. Organick, C. Bee, S. D. Ang, P. Weiss, B. Peck, G. Seelig, L. Ceze, K. Strauss, "Quantifying molecular bias in DNA data storage," *Nat. Commun.* 11:3264, 2020. doi:10.1038/s41467-020-16958-3 | Uses UMIs to separate synthesis bias from PCR bias. Synthesis bias is significant and depends on the oligo's spatial location on the chip. PCR bias comes from stochasticity, not GC content. Code and data: github.com/uwmisl/storage-biasing-ncomms20 |
| [9] | W. H. Press, J. A. Hawkins, S. K. Jones Jr., J. M. Schaub, I. J. Finkelstein, "HEDGES error-correcting code for DNA storage corrects indels and allows sequence constraints," *PNAS* 117(31):18489–18496, 2020. doi:10.1073/pnas.2004821117 | Hash-based inner code that corrects indels. Its RS outer code is applied diagonally across strands. Remains error-free at about 10% nucleotide error. |
| [10] | P. L. Antkowiak, J. Lietard, M. Z. Darestani, M. M. Somoza, W. J. Stark, R. Heckel, R. N. Grass, "Low cost DNA data storage using photolithographic synthesis and advanced information reconstruction and error correction," *Nat. Commun.* 11:5345, 2020. doi:10.1038/s41467-020-19148-3 | Cheap, high-error light-directed synthesis (60 nt strands). Shows that cheap, noisy writing can be compensated at decode time. |
| [11] | A. L. Gimpel, W. J. Stark, R. Heckel, R. N. Grass, "A digital twin for DNA data storage based on comprehensive quantification of errors and biases," *Nat. Commun.* 14:6026, 2023. doi:10.1038/s41467-023-41729-1 | **Most complete calibration source.** Characterizes 40 datasets, splits error rates by process step, and fits lognormal coverage. Data: ENA PRJEB65931. Code: github.com/fml-ethz/dt4dds |
| [12] | S. M. H. Tabatabaei Yazdi, R. Gabrys, O. Milenkovic, "Portable and error-free DNA-based data storage," *Sci. Rep.* 7:5011, 2017. doi:10.1038/s41598-017-05188-1 | Early nanopore (MinION) readout of about 3 kB. |
| [13] | R. Lopez, Y.-J. Chen, S. Dumas Ang, S. Yekhanin, K. Makarychev, M. Z. Racz, G. Seelig, K. Strauss, L. Ceze, "DNA assembly for nanopore data storage readout," *Nat. Commun.* 10:2933, 2019. doi:10.1038/s41467-019-10978-4 | Concatenates oligos into long molecules for nanopore readout. |

## 2. Information-theoretic channel models

- [14] I. Shomorony, R. Heckel, "DNA-based storage: Models and fundamental limits," *IEEE Trans. Inf. Theory*
  67(6):3675–3689, 2021. doi:10.1109/TIT.2021.3058966 (arXiv:2001.06311). Defines the *noisy shuffling-sampling
  channel* and shows that index-based coding is optimal in many regimes. Earlier ISIT version: Heckel, Shomorony,
  Ramchandran, Tse, ISIT 2017 (arXiv:1705.04732).
- [15] N. Weinberger, N. Merhav, "The DNA storage channel: Capacity and error probability bounds," *IEEE Trans. Inf.
  Theory* 68(9):5657–5700, 2022. doi:10.1109/TIT.2022.3176371 (ISIT 2022, doi:10.1109/ISIT50566.2022.9834482;
  arXiv:2109.12549). The channel samples M molecules N times with replacement, with coverage depth N/M. The dominant
  error event is an "outage" caused by sampling.
- [16] N. Weinberger, "Error probability bounds for coded-index DNA storage," arXiv:2205.10077, 2022. The error
  probability decays exponentially in N, so coverage depth matters for reliability more than for capacity.
- [17] A. Lenz, P. H. Siegel, A. Wachter-Zeh, E. Yaakobi, "Coding over sets for DNA storage," *IEEE Trans. Inf.
  Theory* 66(4):2331–2351, 2020. doi:10.1109/TIT.2019.2961265. Gives Gilbert–Varshamov and sphere-packing bounds and
  near-optimal constructions for unordered sets of strands. (These authors also have papers on achievable rates for the
  DNA storage channel; I did not check their exact titles or venues, so they are **[unverified]** and not listed here.)
- [18] S. Kas Hanna, "On the reliability of information retrieval from MDS coded data in DNA storage," arXiv:2502.06618
  (the comment says "accepted in IEEE TMBMC"). Gives success probability as a function of the read count, how reads are
  spread across strands, and the inner and outer RS rates, and works out the optimal split of redundancy between inner
  and outer codes.

## 3. The coverage-depth line of work (closest to C3 and C6)

This is the most important body of prior work for D3RA-DNA.

- [19] **D. Bar-Lev, O. Sabary, R. Gabrys, E. Yaakobi, "Cover your bases: How to minimize the sequencing coverage in
  DNA storage systems,"** ISIT 2023, pp. 370–375 (doi:10.1109/ISIT54713.2023.10206882). Journal version: *IEEE Trans.
  Inf. Theory* 71(1):192–218, 2025 (doi:10.1109/TIT.2024.3496587). arXiv:2305.05656.
  Introduces the *DNA coverage depth problem*: the expected number of reads until decoding succeeds. Main results:
  - Under noiseless uniform sampling, MDS codes are optimal (a coupon-collector argument).
  - The uniform distribution minimizes the expected number of reads.
  - Bounds on the distribution of the number of reads.
  - Random-access variant: the reads needed to recover one information strand.
  The arXiv v3 cites about $120/GB as the cost of sequencing information.
- [20] A. Gruica, D. Bar-Lev, A. Ravagnani, E. Yaakobi, "A combinatorial perspective on random access efficiency for DNA
  storage," ISIT 2024, pp. 675–680 (doi:10.1109/ISIT57864.2024.10619151). The arXiv:2401.15722 version is titled
  "Reducing coverage depth in DNA storage…". Studies T_max(G), the worst-case expected reads to recover any single
  information strand under generator matrix G.
- [21] H. Abraham, R. Gabrys, E. Yaakobi, "Covering all bases: The next inning in DNA sequencing efficiency,"
  arXiv:2401.17649, 2024 (preprint). Studies the coverage depth of practical code pairings in random-access settings.
- [22] H. Abraham, (first name unverified) Feldman, E. Yaakobi, "Coverage depth analysis: Random access for DNA storage
  over noisy channels," ISIT 2026. doi:10.1109/ISIT62367.2026.11653660. Only the Crossref metadata is verified; I did
  not read the paper.
- [23] R. Cao, X. Chen, "Optimizing sequencing coverage depth in DNA storage: Insights from DNA storage data," ISIT 2025
  (doi:10.1109/ISIT63088.2025.11195510; arXiv:2501.06801). Moves from uniform to *log-normal* sampling channels fitted to
  real PCR and sequencing data. Coverage depth grows with MDS code rate. Follow-up: R. Cao, (Zhou), X. Chen, "Modeling
  DNA storage retrieval reliability via sequencing coverage depth," *Brief. Bioinform.* 27(3), 2026,
  doi:10.1093/bib/bbag303. I verified only the Crossref metadata of the follow-up.
- [24] M. Bertuzzo, A. Ravagnani, E. Yaakobi, "The coverage depth problem in DNA storage over small alphabets,"
  arXiv:2507.20639, 2025. Same authors, extended: "The DNA coverage depth problem: duality, weight distributions, and
  applications," arXiv:2603.06489, 2026. Gives closed forms for simplex, Hamming, Golay and Reed–Muller codes. For MDS
  codes, E[reads] = n(H_n − H_{n−k}).
- [25] Ş. Bodur, S. Lia, H. H. López, R. Ludhani, A. Ravagnani, L. Seccia, "The random variables of the DNA coverage depth
  problem," arXiv:2507.20645, 2025. Studies the full distribution, not just the expectation.
- [26] C. Wang, E. Yaakobi, "Random access in DNA storage: Algorithms, constructions, and bounds," arXiv:2601.07053, 2026.
- [27] R. Gross, Y. Grunbaum, M. Bertuzzo, E. Yaakobi, A. Ravagnani, "Optimal codes for the coverage depth problem and the
  performance of random codes," arXiv:2610.07370, 2026.
- [28] A. Gruica, A. Petrillo, F. Zullo, "The generalized random access problem for linear codes," arXiv:2608.20152, 2026.
  Covers recovery of several symbols at once.
- [29b] **D. Levy, R. Con, E. Yaakobi, H. M. Kiah, "Expected recovery time in DNA-based distributed storage systems," arXiv:2602.07601, Feb 2026 (preprint).** Studies multi-container DNA storage *with container failures* and expected recovery time. This is closer prior work for the failure aspect than [29]; it was added after the novelty review (`docs/novelty_assessment.md`).
- [29] **X. Kong, O. Elishco, C. Wang, T. M. Duman, "The coverage depth problem in distributed DNA data storage,"
  arXiv:2610.02931, Oct 2026 (preprint).** *Together with [29b], this is the closest prior work to the D3RA-DNA multi-container setting.*
  Model: strands are partitioned among **M containers**, and each round draws one strand per container, uniformly and
  with replacement. Results:
  - Exact recovery-time distribution and expectation for any linear code and partition.
  - MDS codes are optimal for every fixed partition.
  - A universal lower bound on total read cost.
  - Container-size regimes where distribution saves reads, versus regimes where it only adds parallelism.

  What it does **not** model, as far as the abstract shows:
  - container failure or loss, or heterogeneous container efficiency;
  - sequence-level abundance bias;
  - adaptive or non-uniform allocation of reads across containers;
  - placement driven by failure tolerance.

## 4. Consensus, trace reconstruction and clustering (inner decoding for C6's calibrator)

- [30] T. Batu, S. Kannan, S. Khanna, A. McGregor, "Reconstructing strings from random traces," *SODA 2004*,
  pp. 910–918. Introduces the bitwise majority alignment (BMA) algorithm.
- [31] S. R. Srinivasavaradhan, S. Gopi, H. D. Pfister, S. Yekhanin, "Trellis BMA: Coded trace reconstruction on IDS
  channels for DNA storage," ISIT 2021 (extended version arXiv:2107.06440). Released the Clustered Nanopore Reads
  dataset; see the calibration document.
- [32] O. Sabary, A. Yucovich, G. Shapira, E. Yaakobi, "Reconstruction algorithms for DNA-storage systems," *Sci. Rep.*
  14:1951, 2024. doi:10.1038/s41598-024-51730-3. Works at error probabilities up to about 0.27.
- [33] C. Rashtchian, K. Makarychev, M. Rácz, S. D. Ang, D. Jevdjic, S. Yekhanin, L. Ceze, K. Strauss, "Clustering
  billions of reads for DNA data storage," *NeurIPS (NIPS) 2017*. Distributed edit-distance clustering, about 1000× faster
  than the prior method.
- [34] S. Chandak et al. (Chandak, Neu, Tatwawadi, Mardia, Lau, Kubit, Hulett, Griffin, Wootters, Weissman), "Overcoming
  high nanopore basecaller error rates for DNA storage via basecaller-decoder integration and convolutional codes,"
  ICASSP 2020, pp. 8822–8826. doi:10.1109/ICASSP40776.2020.9053441. Basecalled error rate is about 10%.

## 5. Multi-pool, partitioned and "file-system" DNA storage

- [35] K. J. Tomek, K. Volkel, A. Simpson, A. G. Hass, E. W. Indermaur, J. M. Tuck, A. J. Keung, "Driving the scalability
  of DNA-based information storage systems," *ACS Synth. Biol.* 8(6):1241–1248, 2019. doi:10.1021/acssynbio.9b00100.
  Nested primer addressing and file separation in a database of about 5 TB scale.
- [36] K. N. Lin, K. Volkel, J. M. Tuck, A. J. Keung, "Dynamic and scalable DNA-based information storage," *Nat. Commun.*
  11:2981, 2020. doi:10.1038/s41467-020-16797-2 (DORIS: non-destructive, PCR-free file access).
- [37] J. L. Banal, T. R. Shepherd, J. Berleant, H. Huang, M. Reyes, C. M. Ackerman, P. C. Blainey, M. Bathe, "Random
  access DNA memory using Boolean search in an archival file storage system," *Nat. Mater.* 20(9):1272–1280, 2021.
  doi:10.1038/s41563-021-01021-3. Silica-encapsulated files with barcodes, selected by Boolean FACS. Each file
  lives in its own physical "capsule," which is a natural failure domain.
- [38] C. Bee, Y.-J. Chen, M. Queen, D. Ward, X. Liu, L. Organick, G. Seelig, K. Strauss, L. Ceze, "Molecular-level
  similarity search brings computing to DNA data storage," *Nat. Commun.* 12:4764, 2021. doi:10.1038/s41467-021-24991-z.
- [39] D. Lin, Y. Tabatabaee, Y. Pote, D. Jevdjic, "Managing reliability skew in DNA storage," *ISCA 2022*, pp. 482–494.
  doi:10.1145/3470496.3527441. Reliability varies with position inside a molecule. The paper argues that plain unequal
  error protection cannot close this gap and proposes reliability-aware data mapping instead. This is relevant to C4.
- Organick et al. [6] already used **9 physically separate synthesis pools**. Their redundancy was per file and they
  did not consider failures across pools.

## 6. Distributed-storage coding and placement theory (closest to C2)

- [40] A. G. Dimakis, P. B. Godfrey, Y. Wu, M. J. Wainwright, K. Ramchandran, "Network coding for distributed storage
  systems," *IEEE Trans. Inf. Theory* 56(9):4539–4551, 2010. doi:10.1109/TIT.2010.2054295 (regenerating codes).
- [41] P. Gopalan, C. Huang, H. Simitci, S. Yekhanin, "On the locality of codeword symbols," *IEEE Trans. Inf. Theory*
  58(11):6925–6934, 2012. doi:10.1109/TIT.2012.2208937 (locally repairable codes, LRC).
- [42] C. Huang, H. Simitci, Y. Xu, A. Ogus, B. Calder, P. Gopalan, J. Li, S. Yekhanin, "Erasure coding in Windows Azure
  Storage," *USENIX ATC 2012*. **[DOI and pages unverified.]** Production LRC. Its "degraded read" concept parallels
  C2's fallback tier.
- [43] S. A. Weil, S. A. Brandt, E. L. Miller, C. Maltzahn, "CRUSH: Controlled, scalable, decentralized placement of
  replicated data," *SC '06*. doi:10.1109/SC.2006.19. Failure-domain-aware placement rules (host, rack, row).
  Erasure-coded Ceph pools cap the chunks per failure domain, which is the standard form of C2's per-domain cap N−K.
- [44] A. Cidon, S. Rumble, R. Stutsman, S. Katti, J. Ousterhout, M. Rosenblum, "Copysets: Reducing the frequency of data
  loss in cloud storage," *USENIX ATC 2013*, pp. 37–48. Shows that the *combinatorics* of placement drives
  the probability of data loss under correlated failures.
- [45] Z. W. Birnbaum, "On the importance of different components in a multicomponent system," in P. R. Krishnaiah (ed.),
  *Multivariate Analysis II*, Academic Press, 1969, pp. 581–592. I could not access the text; the details come from
  secondary references. The Birnbaum importance ∂R/∂r_i of a k-out-of-n system equals P(exactly k−1 of the other
  components work), which is **the same identity as C3's "pivotal-fragment" gradient**.

## 7. Adaptive sequencing, adaptive sampling and bandits

- [46] M. Loose, S. Malla, M. Stout, "Real-time selective sequencing using nanopore technology," *Nat. Methods*
  13(9):751–754, 2016. doi:10.1038/nmeth.3930 (Read Until).
- [47] Sokolovskii, Ramirez-Garcia, Heinis (first names unverified), "Adaptive sampling in nanopore sequencing for PCR-free
  random access in DNA data storage," bioRxiv 2024, doi:10.1101/2024.11.05.622081 (preprint). Uses Read Until or
  adaptive sampling as PCR-free random access, which works only for sufficiently long strands.
- [48] X. Pic, N. Pinnamaneni, R. Appuswamy, "Adaptive sampling for storage of progressive images on DNA," arXiv:2603.05220,
  2026. Uses adaptive sampling to read only selected resolution layers.
- [49] W. R. Thompson, "On the likelihood that one unknown probability exceeds another in view of the evidence of two
  samples," *Biometrika* 25(3/4):285–294, 1933 **[DOI unverified]**. Tutorial: D. Russo, B. Van Roy, A. Kazerouni,
  I. Osband, Z. Wen, "A tutorial on Thompson sampling," *Found. Trends Mach. Learn.* 11(1):1–96, 2018,
  doi:10.1561/2200000070.
- **Gap.** My searches of arXiv, Crossref and the web, for terms including bandit, adaptive allocation, sequential
  stopping and multiple pools, found **no work that adaptively allocates sequencing reads across DNA storage
  containers or pools**. That includes work using a stopping rule based on per-stripe decodability. Treat this as
  "not found", not as proof that none exists.

## 8. Patents (only those I could verify)

- US 11,365,411 B2, "Efficient random access to DNA-encoded data," Microsoft Technology Licensing LLC. Google Patents
  lists the inventors as Y.-J. Chen, B. Nguyen and K. Strauss. Covers PCR random access in isolated reaction
  volumes or droplets, possibly from different pools. EP counterpart: EP 4093885 A1.
- US 10,742,233 B2, "Efficient encoding of data for storage in polymers such as DNA," Erlich Lab LLC (Y. Erlich), granted
  2020-08-11. This is the DNA Fountain family.
- US 12,562,755, "Multi-tier error correction codes for DNA data storage"; and US 12,430,202 B2, "Nested error correction
  codes for DNA data storage" (search snippets name Western Digital as assignee; **[unverified]**). Both cover
  multi-level ECC across oligos within a pool.
- US 11,755,640, "DNA-based image storage and retrieval." Unequal protection of facial-feature bits. **[assignee
  unverified]**. Relevant prior art for C4.
- I found **no patent on failure-domain-aware placement across DNA containers or on adaptive allocation of reads
  across containers**. I did not run a full patent search with claim-level review; that needs a professional
  search before any patent filing.

---

## 9. Novelty positioning

Summary: the *ingredients* of D3RA-DNA (MDS stripes, failure-domain caps, Birnbaum-style gradients, Thompson sampling,
Gamma-Poisson coverage, Poissonized simulation) are all known. The *combination*, applied to multi-container DNA
storage with failing containers and adaptive closed-loop sequencing, did not appear in any work I found. The single
closest paper is Kong et al. [29] (Oct 2026), and the paper should state the differences from it explicitly.

### C1. Heterogeneous multi-container channel: sequence-intrinsic abundance shared across replicas vs container effects

- **Closest prior work:**
  - Single-pool characterization: Heckel et al. [7], Chen et al. [8], Gimpel et al. [11].
  - Non-uniform (log-normal) sampling in coverage-depth analysis: Cao & Chen [23].
  - Multi-container sampling, but uniform and with no failures: Kong et al. [29].
- **Caution from the evidence:**
  - Chen et al. [8] found synthesis bias tied to the oligo's location on the chip and *no practically important GC
    association for PCR bias*. They concluded that PCR stochasticity dominates.
  - Gimpel et al. [11] found the PCR efficiency spread to be very narrow (σ ≈ 0.005–0.012) and to depend on
    experimental conditions, not on GC.
  - So a bias that is *sequence-intrinsic and shared across replicas* is plausible only when the containers are
    aliquots of the **same synthesis batch**. Separately synthesized replicas should be largely independent.
  - Recommendation: make the cross-container correlation a parameter (ρ ∈ [0,1]) rather than assuming it.
- **What appears new:** an explicit generative model that combines a per-sequence log-abundance shared across
  containers, a per-container efficiency, and per-container failure (total loss). It is modest modelling novelty. Its
  value is that it enables C2 and C3.

### C2. Failure-domain-aware MDS stripe placement and primary/fallback tiering

- **Closest prior work:**
  - CRUSH and Ceph [43]: the per-domain cap ≤ N−K that makes data survive the loss of one domain is standard practice in
    erasure-coded cloud storage.
  - Copysets [44]: placement combinatorics under correlated failure.
  - LRC and Azure [41, 42]: cheap normal-path reads versus degraded reads.
  - Kong et al. [29]: optimize the *partition* of an MDS codeword over containers to minimize total reads, without
    failures.
- **What appears new:**
  - Transferring failure-domain placement to physical DNA containers (vials or capsules [37]).
  - Coupling it to *read cost*. A "primary tier" that covers ≥ K fragments of every stripe lets retrieval sequence only
    the cheap containers. This resembles systematic or degraded reads, but here the cost is sequencing reads and not
    disk I/O, and the reads are sampled with replacement, which makes tier sizing a coverage-depth question.
- **Honest assessment:** the placement constraint on its own is not novel. The read-cost-aware tiering, analysed under
  the coverage-depth model with container failures, did not appear in my searches. Position it as an extension of
  [29] with failures and tiers.

### C3. Adaptive sequencing allocation (pivotal-fragment gradient, Gamma-Poisson posterior, Thompson sampling, stopping rule)

- **Closest prior work:**
  - The coverage-depth line [19–29], in particular Bar-Lev et al. [19] and Kong et al. [29]. These works compute
    expected reads under a **fixed, non-adaptive** sampling law, mostly uniform. Bar-Lev et al. show the uniform law
    minimizes expected reads *when the sampler cannot be controlled*.
  - Reliability-importance theory [45]: ∂P(S≥K)/∂p_i = P(S_{−i}=K−1) is the Birnbaum importance of a component in a
    k-out-of-n system. **Do not claim the identity as new.** Cite it and claim its *use* as an allocation signal.
  - Thompson sampling [49].
  - Read Until and adaptive sampling as the physical actuator [46–48].
- **What appears new:** closed-loop, per-container read allocation that:
  - learns container efficiency and failure online, with a Gamma-Poisson posterior;
  - directs reads to containers holding the stripes nearest their decoding threshold;
  - stops as soon as every stripe has K verified fragments.

  I found no DNA-storage paper that does this. The nearest is adaptive sampling [47, 48], which selects molecules
  by target identity, not by a decodability objective.
- **Suggested comparison baselines:**
  - uniform allocation (Bar-Lev et al. [19]);
  - the best static partition (Kong et al. [29]);
  - an oracle allocation that knows the true container efficiencies.

### C4. Risk-adaptive (unequal) parity across stripes

- **Closest prior work:**
  - Unequal error protection in DNA images (US 11,755,640; PELMI, *Brief. Bioinform.* 2024, **[details unverified]**).
  - Lin et al. [39] on reliability skew. Their argument against naive UEP should be addressed.
  - Kas Hanna [18] on redundancy split between inner and outer codes.
  - DT4DDS [11] used to size redundancy from predicted dropout.
- **What appears new:** setting the parity per *stripe* from a predicted per-stripe failure risk, which comes from
  container placement and per-sequence abundance. This is incremental, and the paper must show gains over uniform
  parity at equal total redundancy.

### C5. Alternating joint optimization (placement, parity and allocation)

- Alternating or block-coordinate optimization is generic. I found no joint optimization of placement, parity and
  sequencing allocation for DNA storage. Novelty depends on the problem formulation, not the algorithm. Report
  convergence behaviour and compare against optimizing each block separately.

### C6. Recovery-time-matrix Monte-Carlo, calibrated from a sequence-level simulator

- **Closest prior work:**
  - Kong et al. [29] derive *exact* recovery-time distributions for uniform per-container sampling. Use these to
    validate the simulator in the uniform, no-failure special case.
  - DT4DDS [11] is a calibrated sequence-level digital twin.
  - Trace-reconstruction decoders [30–32].
- **What appears new:** engineering. Poissonized sequencing gives exponential arrival times per (fragment, container).
  Any policy can then be replayed on the same matrix, which allows 10^4+ trials with common random numbers.
  Present it as a methodological contribution, and show that it agrees with (a) the Kong et al. formulas and (b) a
  full sequence-level simulation (BMA plus inner RS) on a subset of runs.
- **Caution:** with adaptive policies, arrival times depend on the policy's allocation. The matrix must be stored per
  unit of *per-container read effort*, not per unit of wall time, so that replaying it for different policies stays
  valid. Under Poissonization this holds because thinning and superposition preserve independence.

### Statements to avoid in the paper

- Do not say "first to study coverage depth" (see [19]) or "first multi-container coverage analysis" (see [29]).
- Do not say "MDS-optimal placement is new": MDS codes are optimal for every fixed partition [29].
- Do not present the pivotal-fragment gradient as a new derivation; cite Birnbaum [45].
- Do not claim sequence-intrinsic bias is established. The evidence points to synthesis and spatial effects plus PCR
  stochasticity [8, 11].
