"""Digital-twin containers built with DT4DDS (Gimpel et al., Nat. Commun. 2023).

Every container is an *independent* synthesis of the strands placed in it,
followed by its own storage (aging) and amplification (PCR) history and an
Illumina iSeq100 sequencing run, all simulated by DT4DDS with its published
default error and bias parameters. The output is the container's read stream
(R1 sequences), which is then consumed by the same streaming retriever that
reads real FASTQ files: clustering, consensus, inner RS + CRC, outer RS,
per-stripe CRC-32 and SHA-256. Nothing in this path uses the designed strands
to decide whether a fragment is recovered.
"""
from __future__ import annotations

import gzip
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PRIMERS_0 = ("ACACGACGCTCTTCCGATCT", "AGACGTGTGCTCTTCCGATCT")
PRIMERS_2 = ("AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
             "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")


@dataclass(frozen=True)
class ContainerCondition:
    """Handling history of one physically separate container."""
    name: str
    synthesis: str = "twist"          # 'twist' | 'photolithographic' | 'customarray' (117-nt designs only)
    sigma: float = 0.30               # lognormal oligo-count spread of this synthesis batch
    halflives: float = 0.0            # aging (DT4DDS Aging, n_halflives)
    pcr_cycles: int = 15
    pcr_eff: float = 0.95             # mean PCR efficiency
    eff_sd: float = 0.0051            # per-sequence efficiency spread (DT4DDS default)
    coverage: float = 200.0           # oligo copies per design sampled into the container
    failed: bool = False              # container lost: produces no usable reads


def simulate_container(seqs: list[str], cond: ContainerCondition, n_reads: int, seed: int,
                       read_length: int = 150) -> list[str]:
    """Return the R1 read sequences of one container (DT4DDS synthesis -> aging ->
    PCR -> iSeq100 single-end sequencing). A failed container returns no reads."""
    if n_reads <= 0 or not seqs:
        return []
    if cond.failed:
        # a lost/destroyed container still returns reads when sequenced, but none of them
        # derive from the stored strands (contamination, primer dimers, adapter artefacts)
        rng = np.random.default_rng(int(seed))
        return ["".join(rng.choice(list("ACGT"), read_length)) for _ in range(int(n_reads))]
    import dt4dds
    dt4dds.config.show_progressbars = False
    dt4dds.helpers.config.set_random_seed(int(seed))
    dt4dds.properties.set_property_settings(dt4dds.settings.defaults.SequenceProperties(
        efficiency_distribution="normal", efficiency_params={"loc": 1.0, "scale": cond.eff_sd}))
    syn_settings = {"twist": dt4dds.settings.defaults.ArraySynthesis_Twist,
                    "photolithographic": dt4dds.settings.defaults.ArraySynthesis_Photolithographic,
                    "customarray": dt4dds.settings.defaults.ArraySynthesis_CustomArray}[cond.synthesis]()
    syn = dt4dds.processes.ArraySynthesis(syn_settings(oligo_distribution_type="lognormal",
                                                       oligo_distribution_params={"mean": 0, "sigma": cond.sigma}))
    syn.process(list(seqs))
    pool = syn.sample_by_counts(int(cond.coverage * len(seqs)))
    pool = dt4dds.generators.attach_primers_to_pool(pool, *PRIMERS_0)
    pool.volume = 1
    if cond.halflives > 0:
        aging = dt4dds.processes.Aging(dt4dds.settings.defaults.Aging()(), n_halflives=cond.halflives)
        pool = aging.process(pool)
        pool.volume = 1
    pcr = dt4dds.processes.PCR(dt4dds.settings.defaults.PCR_HiFi()(
        primers=list(PRIMERS_2), template_volume=1, volume=20, efficiency_mean=cond.pcr_eff, n_cycles=cond.pcr_cycles))
    pool = pcr.process(pool)
    with tempfile.TemporaryDirectory() as d:
        sbs = dt4dds.processes.SBSSequencing(dt4dds.settings.defaults.iSeq100()(
            output_directory=d, n_reads=int(n_reads), read_length=read_length, read_mode="single-end"))
        sbs.process(pool)
        with gzip.open(Path(d) / "R1.fq.gz", "rt") as fh:
            reads = [line.strip() for k, line in enumerate(fh) if k % 4 == 1]
    return reads
