# D3RA-DNA: Channel Parameter Calibration

Compiled 2026-11-01. Each value has one of three classifications:

- **[M] Experimentally measured**: the value is reported directly in the cited paper's text or tables.
- **[E] Estimated from published evidence**: derived by us (arithmetic on reported values, or a range across papers),
  or reported only by a secondary source.
- **[H] Hypothetically assigned for controlled simulation**: no direct measurement exists, so we choose it and sweep it.

The reference numbers ([n]) match `docs/literature_review.md`. Additional sources are cited inline.

---

## 1. Per-base error rates

### 1a. Array synthesis + PCR + Illumina sequencing-by-synthesis

| Parameter | Value / range | Source | Class |
|-----------|---------------|--------|-------|
| Substitution (whole pipeline) | 0.4% per position (4×10⁻³) | Organick 2018 [6]: Twist synthesis, NextSeq; "average error rate per position of 0.6%" overall | M |
| Deletion (whole pipeline) | 0.2% (2×10⁻³) | Organick 2018 [6] | M |
| Insertion (whole pipeline) | 0.04% (4×10⁻⁴) | Organick 2018 [6] | M |
| Mean across 40 datasets | del 6.7±6.9×10⁻³, sub 7.9±2.0×10⁻³, ins <0.3±0.2×10⁻³ per nt; max del 17.1×10⁻³, max sub 12.5×10⁻³ | Gimpel 2023 [11] (DT4DDS) | M |
| Deletion, electrochemical synthesis (CustomArray/GenScript) | 13.5±2.0×10⁻³ per nt; >5% per nt near the 5′ end | Gimpel 2023 [11] | M |
| Deletion, material-deposition synthesis (Twist) | 0.58±0.15×10⁻³ per nt | Gimpel 2023 [11] | M |
| Sequencer-only (Illumina iSeq 100, PhiX) | sub 1.8±0.8×10⁻³; del and ins each <0.1×10⁻³ per nt. Forward reads 1.1×10⁻³, reverse reads 2.5×10⁻³; strongly cycle-dependent | Gimpel 2023 [11] | M |
| Sequencer-only (Illumina, as cited by Heckel) | sub about 4×10⁻⁴ to 1.5×10⁻³; indels on the order of 10⁻⁶ | Heckel 2019 [7], citing their ref. 30 | E (secondary) |
| Total within-molecule error | about 1–2% per nt, mostly from synthesis and sequencing; synthesis errors are mainly deletions | Heckel 2019 [7] | M |
| PCR-induced substitutions | 1.09×10⁻⁴ per nt per cycle; 61% A→G/T→C | Gimpel 2023 [11] | M |
| Aging-induced substitutions | 1.64×10⁻⁴ per nt per half-life; 77% C→T/G→A | Gimpel 2023 [11] | M |
| Fraction of reads with correct length | about 57% (Goldman 56.6%, Erlich 57.7%) | Heckel 2019 [7], Fig. 3 caption | M |
| Unusable reads (length or quality filter) | 6.5% on average | Organick 2018 [6] | M |
| Photolithographic (light-directed) synthesis | much higher error rate, deliberately; per-type rates not extracted here | Antkowiak 2020 [10] | unverified |

### 1b. Nanopore sequencing

| Parameter | Value / range | Source | Class |
|-----------|---------------|--------|-------|
| Total "coordinate" error, DNA storage oligos, MinION (2017-era) | about 12% (minimum edit distance) | Organick 2018 [6] | M |
| Total basecalled error, DNA storage | about 10%, a large share of it indels | Chandak et al. ICASSP 2020 [34] | M |
| Total error on raw reads, R9.4.1, Q≥10 (genomic DNA) | about 6% mean; about 6% for low-GC and 8% for high-GC reads | Delahaye & Nicolas, *PLOS ONE* 16:e0257521, 2021, doi:10.1371/journal.pone.0257521 | M |
| Per-type split (Bonito basecaller, public raw data, 200-read sample) | del ≈2.24%, sub ≈1.74%, ins ≈0.70% (total ≈4.67%) | Quah, Sella, Heinis, arXiv:2205.05488, 2022 | E (small sample) |
| R10.4 / R10.4.1 chemistry, newer basecallers | modal read accuracy reported around Q20 (≈99%); not measured on DNA storage oligos | Sereika et al., *Nat. Methods* 19:823–826, 2022, doi:10.1038/s41592-022-01539-7 (numbers not extracted) | unverified for oligos |
| Reads of incorrect length (nanopore, concatenated oligos) | >88% of reads | Organick 2018 [6] | M |

**Suggested defaults (sequence-level simulator):**

| Platform | sub | del | ins | Basis | Class |
|----------|-----|-----|-----|-------|-------|
| Illumina | 4×10⁻³ | 2×10⁻³ | 4×10⁻⁴ | Organick [6]. For Twist + iSeq use sub 2×10⁻³, del 6×10⁻⁴ [11] | M |
| Nanopore | 2×10⁻² | 2.5×10⁻² | 1×10⁻² | Consistent with a 6–12% total [6, 34, Delahaye 2021] and a deletion-dominant split [Quah 2022] | E |

For stress tests, sweep the nanopore total from 5% to 15%.

## 2. Homopolymer and context dependence of indels

| Observation | Source | Class |
|-------------|--------|-------|
| Nanopore: homopolymer and short-repeat regions cause about **half of all errors**, mostly **deletions** (homopolymer length is mis-called) | Delahaye & Nicolas 2021 (R9.4.1) | M |
| Nanopore R10.4.1: accuracy collapses for homopolymers longer than 10 nt | Research Square preprint rs-9835208 (not peer-reviewed) | unverified |
| Illumina pipeline: error rates depend on base type and on neighboring bases; insertions and substitutions are biased toward particular bases | Organick 2018 [6], Fig. 3c,d | M (qualitative) |
| Constrained coding (homopolymer and GC limits) may be *less* efficient than randomizing the sequence and absorbing the errors with extra code redundancy, for substitution errors | Weindel, Gimpel, Grass, Heckel, arXiv:2308.05952, 2023 | E |

Suggested simulator model **[H]**:

- Nanopore: per-base deletion probability p_del·(1 + α·(h−1)) inside a homopolymer of length h, with α ∈ [0.3, 1.0].
  Default α = 0.5.
- Illumina: α ≈ 0.

No source gives α directly, so its value is hypothetical.

## 3. Strand dropout (erasures)

| Setting | Dropout | Source | Class |
|---------|---------|--------|-------|
| DNA Fountain, reads downsampled to 750k | 1.3% of oligos missing | Erlich & Zielinski [5] (bioRxiv 10.1101/074237) | M |
| DNA Fountain design redundancy | 7% (72,000 oligos for 67,088 segments) | [5] | M |
| Organick outer RS redundancy | 15% (87% of sequences carry raw data); median decode at 5 reads per sequence (range 4–14) | [6] | M |
| ETH "High PR" pool (high physical redundancy) | 1.6% not observed; **8%** after 4 half-lives of aging | Heckel 2019 [7] | M |
| ETH "Low PR" pool (about 1 copy per sequence) | **34%** not observed (theory 1−e⁻¹ ≈ 36%); **53%** after 4 half-lives | Heckel 2019 [7] | M |
| By physical copies per sequence | <2% dropout at about 10 copies per sequence; logical redundancy needed: 3% at 30 copies, 8% at 10 copies | Chen 2020 [8] | M (model fit to data) |
| Digital-twin best vs worst workflow, same physical redundancy | 98% recovered (best) vs 24% lost (worst) | Gimpel 2023 [11] | M (validated simulation) |
| Koch et al. generational experiment | designed for 80% dropout; DT4DDS shows at most 30% was needed | Gimpel 2023 [11] | E |

**Defaults** for `config/containers.yaml`:

- Baseline per-sequence dropout in a healthy container: **1–2%** [5, 7, 8] **[M]**.
- Aged or degraded container: **8%** [7] **[M]**.
- Near-single-copy container: **≈35%** [7] **[M]**.
- Whole-container failure probability (lost vial, failed extraction): **no published measurement found**. Sweep it over
  {0.01, 0.05, 0.1, 0.2} **[H]**.

## 4. Coverage (read-count) distribution shape

| Source | Fit | Reported dispersion | Class |
|--------|-----|---------------------|-------|
| Organick 2018 [6] | Gamma, moment-matched per file (so reads are negative binomial) | Overall mean 53.8, s.d. 48.7 (pooled across files, CV ≈ 0.91). Per file: mean 6.7 with s.d. 3.4 (CV ≈ 0.51), up to mean 298.6 with s.d. 139.6 (CV ≈ 0.47) | M (CV computed by us: E) |
| Heckel 2019 [7] | Erlich and Goldman data are "approximately negative binomial"; simulation uses Gamma(shape 8, scale 16) | Shape 8 means CV of the abundance = 1/√8 ≈ 0.35 | M (fit) / E (CV) |
| Chen 2020 [8] | Synthesis copies ~ Normal(10⁸, 3.2×10⁷); branching-process PCR; sampled reads | Synthesis CV = 0.32; observed coverage CV 0.41 before PCR and 0.45 after PCR (subsampled to 20×) | M |
| Gimpel 2023 [11] | **Lognormal** | σ = 0.27 (Twist, unconstrained), 0.30 (Twist, GC-constrained), 0.58 (electrochemical, GC-constrained), 1.30 (electrochemical, unconstrained) | M |
| Cao & Chen 2025 [23] | Log-normal sampling channel fitted to their PCR and sequencing data | parameters not extracted | unverified |

**Defaults:**

- Per-sequence intrinsic log-abundance: σ_seq = **0.30**, for Twist-like synthesis **[M, 11]**.
- Stress settings: σ_seq = **0.58** and **1.30** **[M, 11]**.
- Equivalent Gamma shape: k ≈ 1/(e^{σ²}−1). For σ = 0.30 this gives k ≈ 10.6 **[E]**.
- Reads given abundance: Poisson. The result is Gamma-Poisson (negative binomial), matching [6, 7].

## 5. PCR amplification bias

| Quantity | Value | Source | Class |
|----------|-------|--------|-------|
| Per-cycle success probability (model) | 0.95 | Chen 2020 [8] | M (model param fitted to data) |
| Effect of PCR on spread | s.d. of population-fraction change 0.24 before PCR vs 0.37 after (31 cycles) | Chen 2020 [8] | M |
| GC association of PCR bias | slope <0.01, "no practically important association" | Chen 2020 [8] | M |
| Per-sequence efficiency used in simulation | mean factor 1.85 per cycle, s.d. 0.07 (s.d. 0.25 in the per-cycle variant) | Heckel 2019 [7] | H (by those authors) |
| Example: 80% vs 90% efficiency over 60 cycles | ratio (1.8/1.9)^60 ≈ 0.039 | Heckel 2019 [7] | E |
| Relative efficiency (1+ε_i)/(1+ε̄), Twist pools | normal distribution, **s.d. 0.0051** (unconstrained), **0.0048** (GC-constrained) | Gimpel 2023 [11] | M |
| Same quantity re-estimated on Chen (31 cycles), Erlich (90 cycles), Koch (60 cycles) data | s.d. **0.0058–0.012** | Gimpel 2023 [11] | M |
| PCR cycles used | Goldman 22, Erlich 10 (dilution series 40), ETH High PR 65, Low PR 68 | Heckel 2019 [7], Table 1 | M |

**Defaults** for `config/sequencing.yaml`:

- Mean efficiency: ε̄ = 0.90–0.95.
- Spread of relative efficiency: s.d. = 0.0051 [11], with 0.012 as a stress value.
- Number of cycles: 15–30.
- Amplification: a branching binomial process, as in [8, 11].

Under these defaults, PCR adds only a modest spread compared with synthesis (σ_seq) unless cycle counts exceed about 60.

## 6. Sequence design constraints and oligo length

| System | Oligo length | GC constraint | Max homopolymer | Source | Class |
|--------|--------------|---------------|-----------------|--------|-------|
| Goldman 2013 | 117 nt | none explicit (rotating code) | 1 (no repeats, by construction) | [2]; length per [7] | M |
| Grass 2015 | 158 nt per [3]; [7] lists 117 nt (we assume the difference is primers, unverified) | — | — | [3], [7] | M |
| DNA Fountain 2017 | 152 nt | 45–55% | ≤3 | [5] (bioRxiv) | M |
| Organick 2018 | 150–154 nt (110 nt payload + 2×20 nt primers) | — | — | [6] | M |
| Antkowiak 2020 | 60 nt | — | — | [10] | M |
| DT4DDS pools | 143–157 nt | 50% (constrained pools) | — | [11] | M |

**Defaults:**

- Oligo length **150 nt**: 110 nt payload plus 2×20 nt primers [6] **[M]**.
- GC content **45–55%** [5] **[M]**. 40–60% is a common looser alternative; I found no single paper that defines it,
  so treat it as **[H]**.
- Maximum homopolymer run **≤3** [5] **[M]**.

## 7. Sequencing cost (approximate)

| Platform | Cost basis | Per-read order of magnitude | Source | Class |
|----------|-----------|------------------------------|--------|-------|
| Illumina NovaSeq X 25B (2×150) | about $2.0/Gb (GenoHub spec table, derived from reagent cost divided by maximum output); $3.20/Gb at an academic core (TAMU TxGen price list) | 300 bases per read pair gives about 3.3×10⁶ pairs/Gb, so **≈$0.6–1×10⁻⁶ per read pair** | genohub.com/ngs-sequencer/2/illumina-novaseq-x/ ; txgen.tamu.edu price list | E (approximate) |
| ONT PromethION (R10.4.1) | flow cell about $995 in a 4-pack (2025 lab notebook listing $3,980 per 4); 100–290 Gb per flow cell for gDNA, so about $3.4–10/Gb | for 150–200 nt oligo reads about **$10⁻⁶–10⁻⁵ per read**. Short inserts reduce yield, so this is highly uncertain | ONT product pages; third-party estimates (aseq.substack.com) | E (rough) |
| "Sequencing information" for DNA storage | about $120 per GB of stored data | — | Bar-Lev et al. [19], arXiv v3 abstract | E (secondary) |

**Defaults:**

- `cost_per_read`: Illumina = 1e-6, nanopore = 5e-6 (USD) **[E]**.
- Nanopore also offers per-read latency and adaptive sampling [46, 47]. Illumina is batch-only, so adaptive policies
  on Illumina need per-container *runs* or lanes as the decision unit **[H]**.

## 8. Multi-container parameters (D3RA-DNA specific; none measured)

| Parameter | Default | Sweep | Class | Rationale |
|-----------|---------|-------|-------|-----------|
| Cross-container correlation ρ of per-sequence log-abundance | 0.8 if the containers are aliquots of one synthesis; 0.2 if synthesized separately | 0–1 | H | Chen [8] found synthesis bias tied to chip location, so it is shared across aliquots of one batch. PCR bias is stochastic [8, 11] and so is not shared |
| Container efficiency multiplier (lognormal σ_c) | 0.3 | 0–0.8 | H | No source; plausible variation in extraction and handling |
| Container failure probability | 0.05 | 0.01–0.2 | H | No source |
| Partial degradation (extra dropout) | 8% | 2–53% | E | Aged range from Heckel [7] |
| Number of containers N_c | 6–12 | — | H | Organick [6] used 9 separate synthesis pools (M) |

## 9. Public datasets (verified through the ENA Portal API, 2026-11-01)

| Accession | Title (as deposited) | Associated paper | Platform / example runs | Status |
|-----------|----------------------|------------------|-------------------------|--------|
| **PRJEB1186** (ERP002040) | "storage of 5 computer files (739 kB) by coding into synthetic DNA oligos…" | Goldman 2013 [2] (description matches the abstract) | HiSeq 2000; ERR207813, ERR215679 | Verified |
| **PRJEB19305** (ERP021303) | "DNA storage of digital media using DNA Fountain" | Erlich & Zielinski 2017 [5]; cited as such in [11] | MiSeq; ERR1816980 | Verified |
| **PRJEB19307** (ERP021304) | "DNA fountain for DNA storage" | Erlich & Zielinski 2017 [5] | MiSeq; ERR1817032–ERR1817038 | Verified |
| **PRJEB65931** (ERP151026) | "DT4DDS – Characterization of Errors and Biases in DNA Data Storage" | Gimpel 2023 [11] | iSeq 100; ERR12033803 etc. | Verified |
| **PRJEB35217** | (title not retrieved) | Koch et al., cited in the [11] data-availability statement | — | Accession cited by [11]; not opened |
| PRJEB90546 (ERP173549) | "Codec benchmarking for DNA data storage" (ETH FML, 2025) | Six ECC benchmark; paper not identified | iSeq 100 | Accession verified; paper unverified |
| PRJEB32427 (ERP115115) | "DNA storage of digital media using composite DNA" (Technion) | Likely Anavy et al., *Nat. Biotechnol.* 37:1229, 2019, doi:10.1038/s41587-019-0240-x | HiSeq 2500 | Link to paper unverified |
| PRJEB32885 (ERP115620) | "Error profile of the DNA storage channel" (Wellcome Sanger) | not identified | MiSeq | Paper unverified |
| PRJNA978627 (SRP440729) | "DNA storage error model" (LANL) | not identified | NextSeq 500 | Paper unverified |
| GitHub `uwmisl/storage-biasing-ncomms20` | Chen 2020 data and code | [8] | — | Verified (cited in [8] and [11]) |
| Clustered Nanopore Reads (CNR) | 269,709 MinION reads clustered to centers (Twist synthesis, UW MISL) | Trellis BMA [31] | MinION | Verified on a mirror of `microsoft/clustered-nanopore-reads-dataset`. The repository notes issues with how the dataset was generated, so check it before use |
| GitHub `fml-ethz/dt4dds` (Zenodo 10.5281/zenodo.8329043) | DT4DDS simulator | [11] | — | Verified (cited in [11]) |
| Organick 2018 reads | — | [6] | NextSeq / MinION | **No accession found; unverified** |
| Heckel 2019 ETH reads | — | [7] | — | **No accession in the paper text; unverified** |

**Recommended calibration pipeline:**

1. Fit σ_seq (lognormal) and the dropout-versus-coverage curve on PRJEB19305/PRJEB19307 (Fountain) and PRJEB65931
   (DT4DDS).
2. Fit nanopore indel rates and the homopolymer slope α on the CNR dataset.
3. Cross-check the simulator against the DT4DDS tool on identical workflows.
