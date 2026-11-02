# D3RA-DNA: Novelty Assessment of Candidate Contributions

Compiled 2026-11-02. This builds on `docs/literature_review.md`, `reports/technical_report.md` and
`reports/real_data_validation.md`.

**How claims were checked.**
- Every work cited here was located online (arXiv abstract pages, Europe PMC full text, Crossref metadata,
  USENIX/Ceph pages).
- Whenever I read a work's full text (or a substantial part of it), I say so.
- "Abstract only" means I read nothing beyond the abstract.
- Anything I could not confirm is marked **[unverified]**.
- No numbers are invented. Every D3RA-DNA number comes from the project's own reports.

Verdict scale:
- **KNOWN**: the claim already appears in prior work.
- **PARTIALLY KNOWN**: the core is known, and something specific is new.
- **APPARENTLY NEW**: nothing was found. The searches run are listed, and this is not proof of absence.

**Summary**

| # | Candidate contribution | Verdict |
|---|---|---|
| 1 | Spreading threshold (spread vs concentrate across D failure domains) | **KNOWN as a special case** (Leong–Dimakis–Ho 2012). The heterogeneous, multi-stripe and cost-coupled extensions are new but minor. |
| 2 | Post-failure margin rule | **KNOWN in essence** (Bar-Lev et al. 2023/2025; Levy–Con–Yaakobi–Kiah 2026). The real-data confirmation is new. |
| 3 | Adaptive closed-loop sequencing allocation and stopping | Stopping: **KNOWN in practice**. Per-container, decodability-driven allocation: **APPARENTLY NEW**, but its measured value is small. |
| 4 | Correlated nanopore errors / i.i.d.-simulation gap | **PARTIALLY KNOWN**. The quantified gap on one dataset is new, but it has a serious confound (see 4.3). |
| 5 | Shared cross-replicate coverage bias (69% of log-coverage variance) | **PARTIALLY KNOWN / implied** by Chen 2020 and Gimpel 2023. The number is new but descriptive, and the conditions are not independent. |
| 6 | Verification-gated stopping | **KNOWN / standard practice**. It is not a contribution. |

**A gap in the existing literature review.** The review names Kong et al. (arXiv:2610.02931) as the closest work. It
misses an earlier and more directly relevant paper:

- A. Levy, R. Con, E. Yaakobi, H. M. Kiah, "Expected Recovery Time in DNA-based Distributed Storage Systems,"
  arXiv:2602.07601, Feb 2026 (preprint). I read the full text through the arXiv HTML; the last ~7% was not read.

That paper explicitly defines *DNA-based distributed storage systems* in which information is encoded across multiple
DNA containers for robustness against **container failures**. Any claim that D3RA-DNA is the first to treat container
failure in multi-container DNA storage is therefore false. The paper must cite it.

---

## 1. Spreading threshold

**Our claim** (technical report §0 item 3 and §6.2):
- There are D failure domains, each failing independently with probability q.
- An (N,K) MDS stripe spread evenly over the domains tolerates t = ⌊D(N−K)/N⌋ domain failures.
- If (N−K)/N < 1/D, even spreading is strictly worse than concentrating the stripe in one domain.
- Otherwise, spreading wins iff P(Bin(D,q) ≥ t+1) < q.

### 1.1 Closest prior work

**[A] D. Leong, A. G. Dimakis, T. Ho, "Distributed storage allocations," *IEEE Trans. Inf. Theory* 58(7):4733–4752,
2012, doi:10.1109/TIT.2012.2191135** (Crossref-verified). Extended version: arXiv:1011.5287. I read the model and
theorem statements through ar5iv. Model and results:
- One unit-size object is coded across n nodes, with node i storing x_i and Σx_i ≤ T.
- Each node is accessed independently with probability p. Recovery succeeds iff the accessed nodes hold ≥ 1 in total.
- Under MDS coding, a *symmetric allocation* spreads T evenly over m nodes. It succeeds with probability
  **P_S(p,T,m) = P[Bin(m,p) ≥ ⌈m/T⌉]**, and only a small set of m values needs to be checked.
- *Minimal spreading* (m = ⌊T⌋, uncoded replication) is optimal among symmetric allocations if p ≤ 1/⌈T⌉ (their Thm 4).
- *Maximal spreading* (m = n) is asymptotically optimal when p > 1/T (Thm 1).
- The transition between the two sits roughly at p ≈ 1/T. It is not sharp, and the optimal allocation can be
  non-symmetric. Their example is (n,p,T) = (5, 2/3, 7/3).

**Mapping to our setting.** This mapping is mine; the paper treats p as an access probability, not a survival
probability.
- Node = failure domain; p = 1 − q is the domain survival probability.
- One fragment carries 1/K of the object, so T = N/K, and spreading over D domains is m = D.
- Then ⌈m/T⌉ = ⌈DK/N⌉ = D − ⌊D(N−K)/N⌋ = D − t. **Our "tolerates t failures" is exactly Leong's ⌈m/T⌉ survivor
  requirement.**
- Concentrating in one domain is m = 1, with P_S = p. For T < 2 (parity fraction < 1/2), m = 1 is Leong's minimal
  spreading ⌊T⌋.
- Our rule "spread wins iff P(Bin(D,q) ≥ t+1) < q" is precisely **P_S(p,T,D) > P_S(p,T,1)**.
- When (N−K)/N < 1/D, ⌈DK/N⌉ = D, so P_S(p,T,D) = p^D < p. The "strictly worse" part follows in one line.

**[B] Leong, Dimakis, Ho, "Distributed storage allocations for optimal delay"** (arXiv:1106.2581; abstract and search
snippets only). For symmetric allocations, minimal spreading is optimal for short deadlines and maximal spreading for
long ones. The venue was not confirmed **[venue unverified]**. Related: Leong et al., "Symmetric allocations for
distributed storage," GLOBECOM 2010 (Caltech record; not read), and D. Leong's Caltech PhD thesis, which has a section
"Optimality of Maximal Spreading" (table of contents only).

**[C] S. Nath, H. Yu, P. B. Gibbons, S. Seshan, "Subtleties in tolerating correlated failures in wide-area storage
systems," NSDI 2006, pp. 225–238** (abstract and summaries only).
- Under real correlated-failure traces, adding redundancy has strongly diminishing returns.
- Design rankings can flip. Their example: 8-of-16 erasure coding beats 1-of-4 replication under independent failures
  but achieves about 2 fewer nines under correlated failures.
- This is the qualitative "wide spreading can lose to concentration/replication" phenomenon, shown empirically.

**[D] D. Ford et al., "Availability in globally distributed storage systems," OSDI 2010** (PDF read in part).
- A rack-aware policy places no two chunks of a stripe in the same rack.
- For small and medium failure bursts and large encodings, this increases stripe MTTF "by a factor of 3 typically."
- Multi-cell schemes (e.g. RS(6,3)×3) are modelled with each cell treated as a chunk of a multi-cell stripe.
  Recovering by combining chunks across cells is explicitly not considered.

**[E] A. Cidon et al., "Copysets," USENIX ATC 2013** (search summaries only).
- Under correlated failures, *fewer* distinct placement groups (lower scatter width) reduce the frequency of data loss.
- This is a combinatorial concentrate-versus-spread trade-off, but for replication and with node-level rather than
  domain-level failure.

**[F] Ceph erasure-code documentation** (docs.ceph.com, read).
- `crush-failure-domain=rack` "ensures no two chunks are stored in the same rack."
- Most EC pools "require at least k+m CRUSH failure domains."
- M is the number of OSDs that can be lost. The general rule "per-domain chunks ≤ N−K to survive one domain loss" is
  standard practice.

### 1.2 Verdict: KNOWN, as a direct special case of [A]

| Part of our threshold | Status |
|---|---|
| t = ⌊D(N−K)/N⌋ and the binomial survival probability of even spreading | Identical to P_S(p,T,m) in [A] with m = D, T = N/K |
| Spreading is strictly worse when (N−K)/N < 1/D | Immediate from [A]'s formula (needs all D domains: p^D < p). Also folklore in failure-domain placement [F] |
| Spread beats concentrate iff P(Bin(D,q) ≥ t+1) < q | Comparison of two symmetric allocations (m = D vs m = 1) in [A]. [A] solves the more general problem of optimizing over all m |
| Spreading can lose to concentration/replication under correlated failure | Shown empirically in [C]; combinatorial analogue in [E] |

**What is not in [A] and may be new, though all of it is minor:**
1. *Heterogeneous* domain failure probabilities, where one concentrates in the *most reliable* domain. [A] assumes a
   common p, and its non-symmetric optimum is a different effect.
2. *Integer* fragment counts with hierarchical domains (containers inside domains, as in scenario s4).
3. *File-level* survival over many stripes that share the same domains. Stripe events are perfectly correlated
   through the shared domains, so the file-level comparison is not the stripe-level one raised to a power.
4. Coupling to *sequencing read cost*: the cheapest placement concentrates (technical report §0 item 2).

**Precise statement for the paper.** "The homogeneous threshold is a special case of the symmetric-allocation analysis
of Leong, Dimakis and Ho (2012). We apply it to DNA container placement and extend it numerically to heterogeneous,
hierarchical failure domains and to the joint objective with read cost."

---

## 2. Post-failure margin rule

**Our claim** (real_data_validation §6): after a container is lost, degradation is governed by the margin
m = (fragments remaining per stripe) − K. With small margins, the coupon-collector tail explodes. Measured examples:
- margin 6 (round-robin): about +5% to +32% cost;
- margin 2 (tiered): +120–190%, with 45–48% of retrievals finishing within 3× the budget;
- margin 0: total loss.

### 2.1 Closest prior work

**[G] D. Bar-Lev, O. Sabary, R. Gabrys, E. Yaakobi, "Cover your bases,"** ISIT 2023 / *IEEE TIT* 71(1), 2025
(arXiv:2305.05656; verified in the literature review).
- Under uniform sampling, recovering an [n,k] MDS-coded object needs E[reads] = n(H_n − H_{n−k}). The formula is
  restated in Bertuzzo–Ravagnani–Yaakobi, arXiv:2507.20639, per the literature review.
- The same paper also bounds the distribution of the number of reads.
- **Corollary for one stripe.** With R = K + m surviving fragments sampled uniformly, E = R(H_R − H_m).
  - For m = 0 this is R·H_R ≈ R ln R.
  - For m ≥ 1 it is ≈ R ln(R/m).
  - So the single-stripe cost is a direct function of the margin.

**[H] Levy, Con, Yaakobi, Kiah, arXiv:2602.07601 (2026)** (full text read, except the last ~7%).
- Model:
  - M containers of n strands each, with one uniform read per functioning container per step.
  - Row-wise MDS codes of length M and redundancy r.
  - The time to rebuild a *failed container* from the survivors.
- Theorem 1, as I read it: a generalized coupon collector in which each of n rows needs all but ρ of its symbols has
  E[T] = (n/(ρ+1))(ln n + ln C(m+ρ, m−1) + γ ± o(1)).
- Corollary 1 gives E[T_j] = (n/r) ln n + (n/r) ln C(M−1, r) + γn/r ± o(n).
- **The leading-order cost scales as 1/(surplus+1).** This is the asymptotic form of our margin rule, in the DNA
  multi-container setting and including container failure.
- [H] does **not** model heterogeneous container efficiency, abundance bias, sequencing errors, adaptive allocation, or
  placement against failure probabilities. These omissions come from the paper's own model section, as I read it.

**[I] X. Kong, O. Elishco, C. Wang, T. M. Duman, "The coverage depth problem in distributed DNA data storage,"
arXiv:2610.02931, 2 Oct 2026, 33 pp.** I read the abstract only and could not open the full text this session.
- **Covered** (per the abstract):
  - M containers, with one strand per container per round, under noiseless uniform sampling.
  - Exact recovery-time distribution and expectation for arbitrary linear codes and ordered partitions.
  - MDS codes are optimal for each fixed partition.
  - A universal lower bound on expected total read cost.
  - Container-size regimes where distribution saves reads, versus regimes where it only adds parallelism.
  - The q-ary simplex code uniquely minimizes single-container cost, resolving a conjecture of Bertuzzo, Ravagnani and
    Yaakobi.
  - A partition achieving minimum total read cost, plus bounds for intermediate and balanced partitions.
- **Not covered**, as far as the abstract shows:
  - container *failures* or losses;
  - *heterogeneity* in container yield, error or efficiency (only container sizes vary);
  - non-uniform or biased sampling;
  - *adaptive* allocation or stopping, since the per-round sampling is fixed and the optimized partition is a static
    design.
- Whether the full text has a failure or heterogeneity section is **[unverified]**.

**[J] Other related work.**
- Bodur et al., arXiv:2507.20645: the full distribution of coverage depth, i.e. the tail.
- Kas Hanna, arXiv:2502.06618: success probability as a function of the read count and of the inner/outer RS rates.
- Both were verified in the literature review. I did not reread them.

### 2.2 Verdict: KNOWN in essence (a corollary of [G], with the asymptotic multi-row form in [H])

- **Known:** recovery cost depends on the surplus above K through harmonic-number and coupon-collector terms
  ([G] for one stripe, [H] for many rows, including after a container failure).
- **New, but modest:**
  1. Measuring it on real read streams with heterogeneous, biased containers (PRJEB65931).
  2. The *stall* rate: the fraction of retrievals that miss a budget multiple, which the expectation results do not
     give directly.
  3. Linking it to placement choice (tiered vs round-robin) as a deployment trade-off.
- **Wording.** Call it a "consistency check with coupon-collector theory [G, H] on real reads," not a "rule" we
  discovered.

---

## 3. Adaptive (closed-loop) sequencing allocation and stopping

**Our claim** (technical report H3; real_data_validation §5):
- Per-container read allocation that learns container efficiency online (Gamma-Poisson posterior).
- Pivotal (Birnbaum-importance) scoring of stripes.
- Stopping as soon as every stripe holds K verified fragments.

Measured gains:
- Stopping saves 48–69% in simulation, but only 15–17% on real reads.
- Allocation adds −3.5% to −10% on real reads.
- Pivotal scoring is indistinguishable from coverage-greedy scoring.

### 3.1 Prior work found

**Stopping when decodable: known in practice.**
- **R. Lopez et al., *Nat. Commun.* 10:2933, 2019** (Europe PMC full text read): "real-time sequencing enables the
  ability to sequence until sufficient coverage has been acquired for successful decoding". Their Fig. 1c is a
  "sequence-until" diagram, in which the run is stopped or switched to another file after decoding. The text does not
  show whether the stopping was automated and decoder-triggered **[unverified]**.
- **Y. Zhang, R. Qin, Q. Ge, Q. Guo, W. Chen, "From spacecraft ranging to massive DNA data storage: composite ranging
  codes…," *Sci. Adv.* 2026, doi:10.1126/sciadv.aec1469** (PMC13089337, full text read in part).
  - Uses "real-time read-by-read decoding".
  - Recovered 12.87 MB in about 20 min at 3.66× coverage, with a 4.9% raw error rate, on nanopore.
  - Decoding runs simultaneously with sequencing.
- **readfish barcode-aware adaptive sampling** (Payne et al., bioRxiv 10.1101/2021.12.01.470722; follow-up Munro et al.,
  *Genome Res.* 2025) (search snippets only): sequencing of individual amplicons or barcodes can be switched off "once
  sufficient data have been generated". This is per-sample stopping on a multiplexed flow cell, which is the natural
  physical actuator for per-container allocation.
- **J. Li et al., "SUSTag-ORCtrL," *Nat. Commun.* 2025, doi:10.1038/s41467-025-64293-2** (Europe PMC abstract):
  "dynamic decision-making whilst sequencing" for PCR-free random access through nanopore signal classification. It
  selects targets; it does not optimize a read budget.
- **ACOR** (Han et al., arXiv:2610.36500, Oct 2026; abstract only): streaming consensus with decoding at "cluster
  closure" for real-time reconstruction. No file-level stopping rule is described.

**Adaptive sampling in DNA storage** (literature review [46–48], re-checked):
- Sokolovskii et al. (bioRxiv 10.1101/2024.11.05.622081; the author list was not confirmed by this search
  **[unverified]**) and Pic, Pinnamaneni, Appuswamy (arXiv:2603.05220, abstract) use Read Until to *select molecules
  by identity*.
- Neither allocates reads by decodability or stops on decoding progress, as far as the abstracts show.

**General adaptive-allocation theory.**
- **S. Bubeck, D. Ernst, A. Garivier, "Optimal discovery with probabilistic expert advice: finite time analysis and
  macroscopic optimality," *JMLR* 14:601–623, 2013** (arXiv:1207.5259; search snippets only). Good-UCB repeatedly
  chooses which of several "experts" (sampling sources) to query so as to maximize the discovery of new items. It uses
  an upper confidence bound on the Good–Turing missing mass, and comes with oracle-regret and macroscopic-optimality
  guarantees.
  - **This is the closest general analogue of choosing which container to sequence.** The difference is the
    objective: unweighted discovery there, versus completing all stripes to K-of-N here.
- **L. Weilguny et al., "BOSS-RUNS," *Nat. Biotechnol.* 41:1018–1025, 2023, doi:10.1038/s41587-022-01580-z**
  (Crossref-verified; abstract-level). Bayesian experimental design computes the expected benefit of each new read and
  redistributes coverage in real time. This is closed-loop coverage allocation in genomics.
- Coupon collection with group drawings or package choice (e.g. "Collector's problem with group drawings," *Adv.
  Appl. Probab.*; Fang & Chang 2014): **[authors and details unverified]**.
- Dobson & Tezcan, adaptive sampling with an unknown number of coupons (*Ann. Oper. Res.* 233, 2015): **[unverified;
  from a search snippet]**.

**Searches that found no per-container, decodability-driven allocation for DNA storage:**
- "DNA data storage stop sequencing once decoding succeeds real-time decoding nanopore" (extended);
- "arXiv DNA storage adaptive sequencing feedback variable-length sampling stopping time coverage depth adaptive reads
  per pool" (extended);
- "coupon collector choice of urn adaptive sampling multiple urns";
- Good-UCB / BOSS-RUNS / readfish follow-ups;
- plus the arXiv, Crossref and web searches already reported in the literature review §7.

### 3.2 Verdict

| Component | Verdict |
|---|---|
| Stop sequencing when the file decodes | **KNOWN** (Lopez 2019; Zhang 2026 Sci. Adv.; readfish per-barcode switch-off) |
| Saving from stopping vs a fixed budget at a reliability target | **PARTIALLY KNOWN.** The saving is the gap between E[T] and the (1−ε)-quantile of T, both of which follow from coverage-depth distribution results ([G], [J]). Measuring it on real multi-container reads (15–17%) is new |
| Online learning of container efficiency plus allocation across containers to complete K-of-N stripes | **APPARENTLY NEW** for DNA storage. The general idea is the same as Good-UCB and BOSS-RUNS |
| Birnbaum ("pivotal") importance as the allocation weight | Identity **KNOWN** (Birnbaum 1969, per the literature review). Using it as the allocation weight is new, but it shows **no measurable benefit** over coverage-greedy scoring in our own results |

**Honest assessment.** The novel part, allocation, is also the part with the smallest measured effect: 0–10% on real
reads, with pivotal ≈ greedy. The part with the large effect, stopping, is not novel. A paper must present this
plainly.

---

## 4. Correlated (strand- or context-specific) nanopore errors

**Our claim** (real_data_validation §4), on Zenodo 14296588 (Srinivasavaradhan et al. clustered nanopore reads):
- Real two-sided BMA consensus gives 0.41 / 0.75 / 0.87 exact at 5 / 10 / 20 reads.
- Independent errors at the measured rates give 0.81 / 0.99 / 1.00.
- Replaying real edit scripts independently across reads gives 100% at 20 reads.
- Conclusion drawn: the gap comes from errors correlated across reads of the same strand.

### 4.1 Prior work

**Nanopore errors are systematic and context-dependent: KNOWN.**
- **Delahaye & Nicolas, "Sequencing DNA with nanopores: troubles and biases," *PLoS One* 16(10):e0257521, 2021**
  (Crossref-verified; abstract-level):
  - errors in homopolymers and short repeats make up about half of all errors, mostly deletions;
  - GC content strongly affects error rates;
  - quality scores overstate the error rate.
- **Wick, Judd, Holt, "Performance of neural network basecalling tools for Oxford Nanopore sequencing," *Genome Biol.*
  20:129, 2019** (Crossref-verified; snippet-level): consensus accuracy is about 99.5%, "revealing that even the best
  basecallers still have systematic errors"; methylation-aware polishing raised it to about 99.9%.
- **Huang, Liu, Shih, "Homopolish," *Genome Biol.* 22:95, 2021** (Crossref-verified; snippet-level): many systematic
  errors remain after polishing "because the correct allele is a minority at these loci." Extra reads therefore
  cannot fix them. **This is exactly the cross-read-correlation mechanism.**
- **Hamoum, Dupraz, Conde-Canencia, Lavenier, "Channel model with memory for DNA data storage with nanopore
  sequencing,"** ISTC 2021 (HAL hal-03337117); journal extension in *IEEE Access* 2023 (abstract-level; exact title and
  DOI **[unverified]**).
  - A k-order Markov channel with memory and dependence on the input sequence, built *for DNA storage*.
  - It fits experimental data better than memoryless models (lower KL divergence).
- Chen, Milenkovic, Ribeiro, Sima (2026), on context-dependent deletions: **[unverified; search snippet only]**.

**The i.i.d.-simulation versus real-consensus gap.**
- M. Keoliya, UPenn UROP report and ISPASS 2022 poster (not peer reviewed; snippet-level): simulated per-strand
  accuracy was "consistently greater" than real accuracy, even after controlling for coverage. This is the same
  direction as our finding. It is weak, but it is prior evidence.
- DNAformer (Bar-Lev, Orr, Sabary, Etzion, Yaakobi, *Nat. Mach. Intell.* 7:639–649, 2025,
  doi:10.1038/s42256-025-01003-z) (arXiv v3 full text read):
  - It trains on simulated data whose error statistics are fitted with SOLQC.
  - On Organick et al. *Illumina* data, simulated training matched or beat real-data training (Suppl. Table 9).
  - It does not quantify an i.i.d. gap for nanopore.
  - It includes a confidence filter and a "safety margin."
- I found no paper that reports the specific "real edit scripts replayed independently versus real clusters"
  ablation. That ablation appears new.

**Hierarchical and over-dispersed error models.** Beta-binomial models with an overdispersion parameter are standard
in deep-sequencing variant calling, for example the deepSNV Bioconductor package (vignette read; the original paper
was not checked). I found no use of them for the DNA-storage consensus success curve φ(n), but this search was shallow
**[not exhaustively searched]**.

### 4.2 Verdict: PARTIALLY KNOWN

- **Known:**
  - Nanopore errors are systematic and context-dependent.
  - Consensus has a systematic floor.
  - Memory and context channel models exist for nanopore DNA storage [Hamoum et al.].
  - Simulators overstate real accuracy (weak, unreviewed evidence).
- **New:**
  - A clean ablation that separates rate and within-read structure (replayed scripts give 100%) from cross-read
    dependence (real clusters give 89%).
  - The 20–40 percentage-point gap at 5–10 reads and the ~10% plateau at 20 reads, on this dataset.

### 4.3 Serious confound that must be addressed before publication

The Trellis BMA paper (Srinivasavaradhan, Gopi, Pfister, Yekhanin, arXiv:2107.06440; full text read) contains a "Note
added on 8/12/2024":
- The 10,000 reference sequences of this very dataset "exhibit long-range dependencies instead of being uniformly
  random" because of an error in how they were generated.
- As a result, "some recovered clusters may be malformed, making the trace reconstruction problem harder."

Consequences for our claim:
- Mis-clustered reads, or a wrong reference, would produce *exactly* our signature:
  - consistent failure on two disjoint read sets;
  - a coverage-independent plateau;
  - 100% success when scripts are replayed onto the correct reference.
- So **our "cross-read correlation" may partly or largely be a dataset artefact**, not a channel property. Synthesis
  errors (wrong molecules) are a second explanation.

Before claiming anything:
- filter clusters by read-to-reference identity or by outliers in within-cluster edit distance;
- repeat the analysis on a second nanopore DNA-storage dataset with trustworthy references;
- report the result for each explanation separately.

---

## 5. Shared (cross-replicate) coverage bias

**Our claim:** in ENA PRJEB65931 (one Twist pool, 6 handling conditions), 69% of log-coverage variance is shared by a
strand across all six conditions, with σ_g = 0.31. The reports call this "sequence-intrinsic."

### 5.1 Prior work (both full texts read through Europe PMC)

**Chen et al., *Nat. Commun.* 11:3264, 2020.**
- Synthesis bias is significant and "highly related to spatial location of oligos on a synthesis chip."
- PCR bias "is not correlated with GC content", and PCR stochasticity is the main PCR effect.
- A bias tied to a molecule's position on the chip is inherited by every aliquot of that synthesis batch. This
  *implies* shared cross-replicate bias, although Chen et al. do not quantify it as a variance share.

**Gimpel et al. (DT4DDS), *Nat. Commun.* 14:6026, 2023.**
- Synthesis skew is modelled as the *initial count distribution* of the pool, and every downstream step inherits it.
- Per-sequence PCR efficiency is estimated by relating each sequence's coverage before and after PCR, which presupposes
  sequence-level persistence. The spread is σ ≈ 0.0058–0.012 and does not depend on GC.
- Aging "nor significantly affected the coverage distribution", and decay is modelled as random sampling.
- In the main text, the only occurrence of "correlat*" is a limitations sentence saying they did not investigate
  "correlations with sequence properties."
- I found **no explicit decomposition of shared versus condition-specific variance** in the main text. The
  Supplementary material was not checked **[unverified]**.
- Methods: the PCR series was produced by **sequential re-amplification with dilution** (15 → 90 cycles over six
  rounds), and the aged samples come from the same pool.

### 5.2 Verdict: PARTIALLY KNOWN / implied

The direction of the result is implied by Chen 2020 and by the structure of DT4DDS. The 69% figure is a new
descriptive statistic that DT4DDS does not report (main text). It must be reworded:

1. **The six conditions are not independent containers.**
   - The PCR 60c and 90c samples are descendants of the 15c lineage.
   - All aging samples are aliquots of one pool.
   - Part of the shared variance is inherited PCR history, so 69% is an *upper bound* for independently handled
     aliquots and says nothing about separately synthesized pools.
2. **"Sequence-intrinsic" is the wrong word.** Chen 2020 attributes synthesis bias to chip location, not to sequence
   content, and our own GC correlation is about zero. Use "batch-inherited per-strand abundance."
3. **Containers have to be independent physical aliquots for this to transfer to D3RA-DNA.** It is still a useful
   parameter estimate (ρ for same-batch aliquots), but it is not a discovery.

---

## 6. Verification-gated stopping

**Our claim:** stop only when the inner RS/CRC pass for ≥ K fragments of every stripe and the file's SHA-256
verifies.

**Prior practice:**
- DNA Fountain (Erlich & Zielinski 2017): each droplet carries 2 RS bytes, and the decoder "rejected droplets with
  errors based on the Reed-Solomon code" before LT decoding (patent US 10,742,233 and press summaries; search-level).
- Grass et al. 2015 use inner/outer RS.
- The DNAformer comparison table lists several prior schemes with an "inner code CRC."
- Lopez 2019 and Zhang 2026 sequence until decoding succeeds.
- Outside DNA, CRC-gated acknowledgement for incremental redundancy (ARQ/HARQ, rateless codes) is textbook practice.
  This is general knowledge; no specific citation was checked.

**Verdict: KNOWN / standard engineering.** A whole-file cryptographic hash as the final gate is good practice, but I
found no paper that presents it as a contribution, and it should not be presented as one. Its legitimate role in the
paper is methodological: it makes "exact recovery" an auditable endpoint. It is worth one sentence.

---

## 7. Strongest defensible contribution (ranked)

1. **Methods paper with real-data evaluation: "Multi-container DNA retrieval under heterogeneity and container loss:
   placement and adaptive stopping evaluated on real read streams."** This is the strongest option. Defensible
   content:
   - A resampling method that emulates any placement across six real sequencing runs.
   - On PRJEB65931, placement plus adaptive stopping needs 34% fewer reads than uniform placement with fixed depth.
     Half of that comes from placement and half from stopping.
   - The honest negative result: stopping gains drop from 48–69% in simulation to 15–17% on real reads.
   - Allocation is worth ≤ 10%, and pivotal scoring ties coverage-greedy scoring.
   - The surrogate is 5–18% optimistic, which is a useful warning for model-based budgeting.
   - Failure injection on real streams.
   - Byte-exact end-to-end recovery.

   Requirements before submission:
   - Baselines and theory checks against Kong et al. [I] (static optimal partition) and Levy et al. [H] (expected
     failed-container recovery time), validated in their uniform special cases.
   - Cite Lopez 2019 and Zhang 2026 for stopping.
   - A second dataset if possible.

   Suitable venues **[venue fit is my judgement]**: *IEEE TMBMC*, *Bioinformatics*, or *ACS Synth. Biol.*. The
   claimed novelty is the systematic real-data evaluation, not the individual ideas.
2. **Short empirical note on nanopore consensus under cross-read dependence.** This is viable only after resolving the
   §4.3 confound (cluster-quality filtering plus a second dataset). If the gap survives, the "replayed-scripts"
   ablation is a crisp, citable measurement.
3. **Theory note on placement under failure domains.** It is not publishable as stated, since the homogeneous result is
   in Leong et al. 2012. A note would have to contribute one of:
   - closed-form or bounded results for *heterogeneous* domain failure probabilities;
   - file-level (multi-stripe) correlated survival;
   - the joint read-cost and reliability optimum, e.g. a condition under which a single container is cost-optimal at
     target ε.

   The design rule "distribution earns its place only when no single container meets ε" (technical report §0 item 2)
   could be formalized here. I found no prior statement of it coupled to read cost, but it is close to trivial for a
   single stripe.
4. **Software and simulator artefact** (recovery-time-matrix replay with common random numbers). Only as a
   supplementary contribution or a JOSS-style paper.

---

## 8. Claims that must be dropped or reworded

| Current wording (location) | Problem | Replace with |
|---|---|---|
| "new spreading-threshold result" (tech report §0.3); "To our knowledge it has not been stated for multi-container DNA archives" (§6.2) | It is a special case of Leong–Dimakis–Ho 2012 (P[Bin(m,p) ≥ ⌈m/T⌉], comparing m = D with m = 1). Nath et al. 2006 showed the empirical phenomenon | "A direct application of the symmetric-allocation analysis of Leong et al. (2012) to DNA containers, extended numerically to heterogeneous, hierarchical domains" |
| "Real-data form of the spreading rule: the post-failure margin governs degradation" (real_data §6; tech report §11) | A corollary of coupon-collector results: Bar-Lev et al. E = n(H_n − H_{n−k}); Levy et al. 2026, E ∝ 1/(surplus+1) | "Consistent with coupon-collector theory (Bar-Lev et al.; Levy et al.), confirmed on real heterogeneous read streams" |
| Kong et al. described as the closest prior work (literature review §3, §9; tech report §10) | Levy, Con, Yaakobi, Kiah (arXiv:2602.07601, Feb 2026) already define DNA distributed storage *with container failures* | Cite both. State that ours adds heterogeneity, bias, errors, adaptivity and real data |
| "Adaptive stopping is the largest single lever" presented as a contribution | Sequence-until-decoded is known (Lopez 2019; Zhang 2026) | "We quantify the value of stopping vs fixed depth at a reliability target; it is large in simulation and 15–17% on real reads" |
| Pivotal-fragment gradient as a key algorithmic idea | The identity is Birnbaum importance, and it gave no gain over coverage-greedy scoring | Mention it as a design choice, report the tie, and do not headline it |
| "A shared, sequence-intrinsic abundance component exists and is large (69%)" (real_data §3) | The conditions are lineage-related (serial PCR, same-pool aging), and the bias is chip/batch-related, not sequence-intrinsic | "Up to 69% of log-coverage variance is shared across lineage-related samples of one synthesis batch; this is an upper bound for same-batch aliquots" |
| "The whole gap therefore comes from errors correlated across reads of the same strand" (real_data §4) | It cannot exclude malformed clusters or reference problems documented for this dataset (Trellis BMA note, 2024) or synthesis errors | "…comes from dependence among reads within a cluster, which may reflect systematic basecalling errors, synthesis variants or clustering artefacts" |
| "Independent-error simulators overstate nanopore consensus success by 20–40 pp" (general) | One dataset, one consensus algorithm, with a known dataset caveat | Restrict the claim to this dataset and BMA, and cite Keoliya (simulators overstate) and the systematic-error literature as prior evidence |
| Any implication that verification-gated stopping is novel | It is standard practice (inner RS/CRC rejection; decode-until-success) | One methods sentence |
| "first multi-container…", "first to model container failure in DNA storage" (anywhere) | False given Levy et al. 2026 (failures) and Kong et al. 2026 (multi-container) | Do not use |

---

## 9. Items still to verify before submission

- The full text of Kong et al. arXiv:2610.02931: do any sections treat failures or heterogeneity?
- The last ~7% of Levy et al. arXiv:2602.07601, and whether it cites or extends to heterogeneous containers.
- Whether the DT4DDS Supplementary Information reports cross-condition per-sequence coverage correlation.
- The exact titles, venues and DOIs for Hamoum et al. (IEEE Access 2023) and Leong et al. (optimal delay).
- A claim-level patent search on adaptive per-pool read allocation for DNA storage. None was found at search level;
  this requires a professional search.

## Sources consulted (URLs)

- Leong et al. 2012: https://arxiv.org/abs/1011.5287 , https://ar5iv.arxiv.org/html/1011.5287
- Leong et al., optimal delay: https://arxiv.org/pdf/1106.2581
- Nath et al. 2006: https://www.usenix.org/conference/nsdi-06/subtleties-tolerating-correlated-failures-wide-area-storage-systems
- Ford et al. 2010: https://www.usenix.org/legacy/event/osdi10/tech/full_papers/Ford.pdf
- Copysets: https://www.usenix.org/conference/atc13/technical-sessions/presentation/cidon
- Ceph EC documentation: https://docs.ceph.com/en/latest/rados/operations/erasure-code/
- Levy et al. 2026: https://arxiv.org/abs/2602.07601
- Kong et al. 2026: https://arxiv.org/abs/2610.02931
- Cover your bases: https://arxiv.org/abs/2305.05656
- Lopez et al. 2019: PMC6610119
- Zhang et al. 2026 (Sci. Adv.): PMC13089337
- Li et al. 2025 (SUSTag): PMC12534545
- ACOR: https://arxiv.org/abs/2610.36500
- Pic et al.: https://arxiv.org/abs/2603.05220
- Sokolovskii et al.: https://www.biorxiv.org/content/10.1101/2024.11.05.622081
- Good-UCB: https://arxiv.org/abs/1207.5259
- BOSS-RUNS: doi:10.1038/s41587-022-01580-z
- readfish barcode preprint: https://www.biorxiv.org/content/10.1101/2021.12.01.470722v2
- Delahaye & Nicolas: PMC8486125
- Wick et al.: doi:10.1186/s13059-019-1727-y
- Homopolish: doi:10.1186/s13059-021-02282-6
- Hamoum et al.: https://hal-imt-atlantique.archives-ouvertes.fr/hal-03337117
- Keoliya: https://www.engineering.upenn.edu/~mkeoliya/assets/pdf/ispass_2022_report.pdf
- DNAformer: https://arxiv.org/abs/2109.00031 ; doi:10.1038/s42256-025-01003-z
- Trellis BMA: https://arxiv.org/abs/2107.06440
- Chen et al. 2020: PMC7324401
- Gimpel et al. 2023: PMC10533828
- Sabary et al. 2024: PMC10806084
