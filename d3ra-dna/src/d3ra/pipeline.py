"""High-level end-to-end API: files -> DNA library -> containers -> simulated
sequencing (fixed or adaptive) -> reconstruction -> SHA-256 verification."""
from __future__ import annotations

import numpy as np

from .channel.calibration import calibrate_containers
from .channel.failures import FailureModel
from .channel.surrogate import ChannelParams
from .ingestion.file_loader import FileRecord
from .retrieval.adaptive_policy import AdaptiveConfig
from .retrieval.sequencing_env import SequencingEnv, run_policy
from .storage.library import build_library
from .storage.placement import cost_tiered, domain_aware, round_robin, unit_retrieval_cost


def store_and_retrieve(files: list[FileRecord], containers, placement: str = "tiered", controller: str = "pivotal",
                       payload_bytes: int = 32, k: int = 24, parity: int = 8, inner_nsym: int = 4,
                       codec: str = "scrambled-2bit", seed: int = 0, calibrate: bool = True,
                       failed: np.ndarray | None = None, fixed_reads_per_fragment: float | None = None,
                       cap: float = np.inf) -> dict:
    rng = np.random.default_rng(seed)
    lib = build_library(files, payload_bytes, k, parity, inner_nsym, codec)
    m = len(containers)
    phi = calibrate_containers(containers, trials=30, fmt=lib.fmt) if calibrate else None
    fm = FailureModel(np.array([c.fail_prob for c in containers]))
    ch = ChannelParams.from_containers(containers, fm, phi, L=lib.fmt.strand_nt, nsym=inner_nsym)
    stripe_of = lib.stripe_of
    kk = lib.k_per_stripe()
    dom = np.array([c.domain for c in containers])
    if placement == "tiered":
        X = cost_tiered(stripe_of, kk, unit_retrieval_cost(ch), dom)
    elif placement == "domain":
        X = domain_aware(stripe_of, kk, dom, m)
    else:
        X = round_robin(stripe_of, m)
    for j, c in enumerate(containers):
        c.inventory = [lib.fragments[i].uid for i in np.flatnonzero(X[:, j])]
    env = SequencingEnv(lib, containers, X, rng, failed=failed)
    fixed = None
    if fixed_reads_per_fragment is not None:
        fixed = fixed_reads_per_fragment * X.sum(0) / np.maximum(ch.yield_, 1e-3)
    out = run_policy(env, ch, AdaptiveConfig(controller=controller, cap=cap), fixed_reads=fixed)
    out["library"] = lib.stats()
    out["failed"] = env.failed.tolist()
    return out
