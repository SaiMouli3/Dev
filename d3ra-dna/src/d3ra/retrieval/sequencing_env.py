"""Sequence-level retrieval environment: real strands, noisy reads,
inventory-aware clustering, BMA consensus, inner decoding, verification.

The same controllers used on the surrogate drive this environment one trial
at a time (``run_policy``), which is how surrogate predictions are validated
end-to-end on actual bytes.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..channel.abundance import container_variation, intrinsic_bias
from ..channel.dropout import sample_dropout
from ..channel.sequencer import ReadErrorModel, sequence_container
from ..channel.surrogate import ChannelParams, Layout
from ..decoding.reconstruction import RecoveryResult, reconstruct
from ..storage.library import DNALibrary
from .adaptive_policy import THETA, AdaptiveConfig, AdaptiveController
from .consensus import InventoryClusterer, consensus
from .fragment_tracker import FragmentTracker


@dataclass
class SequencingEnv:
    lib: DNALibrary
    containers: list
    placement: np.ndarray                # (n, M) copies
    rng: np.random.Generator
    failed: np.ndarray | None = None
    sigma_g: float = 0.3
    d_syn: float = 0.002
    sampling: str = "multinomial"
    dup_rate: float = 0.0
    max_cluster_reads: int = 24
    pools: dict = field(default_factory=dict)
    reads_used: np.ndarray | None = None

    def __post_init__(self):
        n, m = self.placement.shape
        strands = self.lib.strands
        if self.failed is None:
            self.failed = np.array([self.rng.random() < c.fail_prob for c in self.containers])
        g = intrinsic_bias(strands, n, self.sigma_g, 0, self.rng)
        h = container_variation(n, m, np.array([c.abundance_shape for c in self.containers]), self.rng)
        self.abundance = g[:, None] * h * self.placement
        self.present = sample_dropout(n, m, np.array([c.dropout for c in self.containers]), self.d_syn, self.rng)
        self.errs = [ReadErrorModel.from_container(c) for c in self.containers]
        fmt = self.lib.fmt
        self.clusterers = []
        for j in range(m):
            idx = np.flatnonzero(self.placement[:, j] > 0)
            self.clusterers.append(InventoryClusterer(list(idx), [strands[i][: fmt.prefix_nt] for i in idx], fmt.prefix_nt))
        k = self.lib.k_per_stripe()
        n_s = np.array([s["n"] for s in self.lib.manifest.stripes])
        self.tracker = FragmentTracker(k, n_s)
        self.rec_pair = {}
        self.reads_used = np.zeros(m)
        self.cluster_reads = {}
        self.decode_attempts = 0
        self.false_accepts = 0

    def sequence(self, j: int, reads: int) -> int:
        c = self.containers[j]
        idx = np.flatnonzero(self.placement[:, j] > 0)
        strands = [self.lib.strands[i] for i in idx]
        out = sequence_container(strands, self.abundance[idx, j], self.present[idx, j], reads, self.errs[j],
                                 self.rng, self.sampling, c.yield_, self.dup_rate, failed=bool(self.failed[j]))
        self.reads_used[j] += reads
        rs = [r for _, r in out]
        keys = self.clusterers[j].assign(rs)
        touched = set()
        for key, r in zip(keys, rs):
            if key is None:
                continue
            pool = self.pools.setdefault((j, key), [])
            if len(pool) < self.max_cluster_reads:
                pool.append(r)
                touched.add(key)
        new = 0
        L = self.lib.fmt.strand_nt
        for i in touched:
            if self.rec_pair.get((i, j)):
                continue
            f = self.lib.fragments[i]
            if (f.stripe, f.index) in self.tracker.rows:
                self.rec_pair[(i, j)] = True
                continue
            est = consensus(self.pools[(j, i)], L)
            self.decode_attempts += 1
            dec = self.lib.fmt.decode(est)
            if dec is None:
                continue
            (obj, s, pos), payload = dec
            if (s, pos) != (f.stripe, f.index):
                continue   # inner code accepted a different (mis-clustered) strand: ignore
            if payload != f.payload:
                self.false_accepts += 1      # wrong fragment passed inner RS + CRC + address check
            self.rec_pair[(i, j)] = True
            new += self.tracker.add(s, pos, payload, j)
        c.record(reads, new)
        return new

    def result(self) -> RecoveryResult:
        return reconstruct(self.lib.manifest, self.tracker.rows)


def run_policy(env: SequencingEnv, ch: ChannelParams, cfg: AdaptiveConfig, fixed_reads: np.ndarray | None = None,
               max_actions: int = 2000) -> dict:
    """Drive the sequence-level env with a fixed allocation or an adaptive controller."""
    lib = env.lib
    stripe_of = lib.stripe_of
    layout = Layout(stripe_of, lib.k_per_stripe(), env.placement)
    if fixed_reads is not None:
        for j, r in enumerate(fixed_reads):
            if r > 0:
                env.sequence(j, int(r))
        res = env.result()
        return {"status": res.status.value, "reads": float(env.reads_used.sum()),
                "cost": float(env.reads_used @ ch.cost), "actions": int((fixed_reads > 0).sum()), "result": res}
    ctl = AdaptiveController(layout, ch, cfg)
    used = layout.stored_per_container() > 0
    repairs = 0
    for j in np.flatnonzero(used):
        env.sequence(j, ctl.probe)
    actions = int(used.sum())
    rr = 0
    while actions < max_actions and env.reads_used @ ch.cost < cfg.cap:
        if env.tracker.complete():
            # verification gate: stop only if every stripe verifies (per-stripe CRC) and the file hash matches
            res = env.result()
            if res.status.value == "exact" or not res.unverified_stripes:
                break
            for s_bad in res.unverified_stripes:      # need one more fragment to localise the bad one
                if env.tracker.k[s_bad] < env.tracker.n[s_bad]:
                    env.tracker.k[s_bad] += 1
            repairs += 1
            if repairs > 20:
                break
        R = env.reads_used[None, :].astype(float)
        rec_pair = np.array([env.rec_pair.get((i, j), False) for i, j in zip(layout.pair_frag, layout.pair_cont)])[None]
        rec_frag = np.zeros((1, layout.n), bool)
        uid = {(f.stripe, f.index): i for i, f in enumerate(lib.fragments)}
        for key in env.tracker.rows:
            rec_frag[0, uid[key]] = True
        by_stripe = env.tracker.by_stripe[None, :]
        U = layout.pair_to_cont_sum(rec_pair.astype(float))
        if cfg.controller == "uniform":
            cand = np.flatnonzero(used)
            j = cand[rr % len(cand)]
            rr += 1
        else:
            post = ctl.posterior(R, U)
            theta = post @ THETA
            sc = ctl.scores(R, rec_pair, rec_frag, by_stripe, theta, env.rng)
            sc[:, ~used] = -np.inf
            j = int(sc.argmax())
        env.sequence(j, ctl.batch)
        actions += 1
    res = env.result()
    return {"status": res.status.value, "reads": float(env.reads_used.sum()),
            "cost": float(env.reads_used @ ch.cost), "actions": actions, "result": res,
            "decode_attempts": env.decode_attempts, "false_accepts": env.false_accepts,
            "verification_rounds": repairs, "repaired": res.repaired}
