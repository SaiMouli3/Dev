"""Closed-loop retrieval from per-container read streams with full decoding.

Used for the digital-twin end-to-end experiments: each container is a list of
reads (in sequencing order). The controller chooses which container to read
next; a batch of reads is clustered against the expected address prefixes,
consensus is formed per cluster, and a fragment counts as recovered only if
the inner Reed-Solomon code, CRC-16 and address check accept it. Stopping is
verification-gated: when every stripe holds K fragments, the stripes are
outer-decoded and checked against their CRC-32 and the blob against its
SHA-256; a failing stripe asks for one more fragment so the bad one can be
localised. No oracle knowledge of the designed strands is used for decisions;
the true payloads are used only afterwards, to count false accepts.
"""
from __future__ import annotations

import time

import numpy as np

from ..archive import address_prefix
from ..channel.surrogate import ChannelParams, Layout
from ..decoding.reconstruction import reconstruct
from .adaptive_policy import AdaptiveConfig, AdaptiveController
from .consensus import InventoryClusterer, consensus
from .fragment_tracker import FragmentTracker

ADAPTER = "AGATCGGAAGAGC"


class StreamRetriever:
    def __init__(self, man, fmt, placement, addresses, streams: list[list[str]], belief: ChannelParams,
                 cfg: AdaptiveConfig, truth: dict | None = None, static_weights=None, cap_reads: float = np.inf,
                 max_cluster_reads: int = 24, two_sided: bool = False, checkpoints=None):
        self.man, self.fmt, self.placement, self.addresses = man, fmt, placement, addresses
        self.streams = streams
        self.pos = np.zeros(len(streams), dtype=int)
        m = placement.shape[1]
        self.clusterers = []
        for j in range(m):
            idx = np.flatnonzero(placement[:, j] > 0)
            self.clusterers.append(InventoryClusterer(list(idx), [address_prefix(fmt, addresses[i]) for i in idx],
                                                      fmt.prefix_nt))
        k = np.array([s["k"] for s in man.stripes])
        n = np.array([s["n"] for s in man.stripes])
        self.tracker = FragmentTracker(k.copy(), n)
        self.k0 = k
        stripe_of = np.repeat(np.arange(len(k)), n)
        self.layout = Layout(stripe_of, k, placement)
        self.cfg = cfg
        self.ctl = AdaptiveController(self.layout, belief, cfg)
        self.belief = belief
        self.truth = truth or {}
        self.static_w = None if static_weights is None else np.asarray(static_weights, float)
        self.cap_reads = cap_reads
        self.max_cluster_reads = max_cluster_reads
        self.two_sided = two_sided
        self.pools: dict = {}
        self.rec_pair: dict = {}
        self.reads_used = np.zeros(m)
        self.decode_attempts = 0
        self.false_accepts = 0
        self.uid = {(a[1], a[2]): i for i, a in enumerate(addresses)}
        self.checkpoints = sorted(checkpoints or [])
        self.cp_rows = []

    # ------------------------------------------------------------------ sequencing
    def _batch(self, j: int, n: int) -> list[str]:
        s = self.streams[j]
        a = self.pos[j]
        b = min(a + n, len(s))
        self.pos[j] = b
        out = []
        for seq in s[a:b]:
            k = seq.find(ADAPTER)
            out.append(seq[:k] if k >= 0 else seq)
        return out

    def exhausted(self) -> np.ndarray:
        return np.array([self.pos[j] >= len(s) for j, s in enumerate(self.streams)])

    def read_container(self, j: int, n: int) -> int:
        reads = self._batch(j, n)
        self.reads_used[j] += n            # reads are paid for even if the stream has run dry
        keys = self.clusterers[j].assign(reads) if reads else []
        touched = set()
        for key, r in zip(keys, reads):
            if key is None:
                continue
            pool = self.pools.setdefault((j, key), [])
            if len(pool) < self.max_cluster_reads:
                pool.append(r)
                touched.add(key)
        new = 0
        L = self.fmt.strand_nt
        for i in touched:
            if self.rec_pair.get((i, j)):
                continue
            addr = self.addresses[i]
            # decode even if another container already supplied this fragment: the posterior
            # over container j's efficiency must only count recoveries made from j's own reads
            self.decode_attempts += 1
            dec = self.fmt.decode(consensus(self.pools[(j, i)], L, two_sided=self.two_sided))
            if dec is None or dec[0] != addr:
                continue
            if self.truth and dec[1] != self.truth.get(addr):
                self.false_accepts += 1
            self.rec_pair[(i, j)] = True
            new += self.tracker.add(addr[1], addr[2], dec[1], j)
        return new

    # ------------------------------------------------------------------ state for the controller
    def _state(self):
        L = self.layout
        rec_pair = np.array([self.rec_pair.get((i, j), False) for i, j in zip(L.pair_frag, L.pair_cont)])[None]
        rec_frag = np.zeros((1, L.n), bool)
        for key in self.tracker.rows:
            rec_frag[0, self.uid[key]] = True
        by_stripe = np.minimum(self.tracker.by_stripe, self.k0)[None, :]
        U = L.pair_to_cont_sum(rec_pair.astype(float))
        return self.reads_used[None, :].astype(float), rec_pair, rec_frag, by_stripe, U

    def choose(self, avail: np.ndarray, rounds: int, rng) -> int:
        name = self.cfg.controller
        cand = np.flatnonzero(avail)
        if self.static_w is not None:              # static allocation: keep reads proportional to w
            w = np.where(avail, self.static_w, 0)
            ok = np.flatnonzero(w > 0)
            if ok.size == 0:                       # every positive-weight container exhausted:
                return int(cand[0])                # read out the rest (retrieval will fail)
            return int(ok[np.argmin(self.reads_used[ok] / w[ok])])
        if name == "uniform":
            return int(cand[rounds % len(cand)])
        R, rec_pair, rec_frag, by_stripe, U = self._state()
        post = self.ctl.posterior(R, U)
        if name == "ltc":
            if rounds < self._explore_end:
                self._committed = post
                return int(cand[rounds % len(cand)])
            post = self._committed
        theta = self.ctl.theta_from_post(post, rng)
        if name == "info":
            v = self.ctl.info_scores(self.ctl.posterior(R, U))
            n_open = self.layout.pair_to_cont_sum((~rec_pair).astype(np.float32))
            sc = v * (n_open > 0)
            if sc.max() <= 1e-3:
                sc = self.ctl.scores(R, rec_pair, rec_frag, by_stripe, theta, rng, base="coverage")
        else:
            sc = self.ctl.scores(R, rec_pair, rec_frag, by_stripe, theta, rng)
        sc = sc[0].copy()
        sc[~avail] = -np.inf
        return int(sc.argmax())

    def _record_checkpoint(self, rng):
        spent = self.reads_used.sum()
        while self.checkpoints and spent >= self.checkpoints[0]:
            level = self.checkpoints.pop(0)
            R, rec_pair, rec_frag, by_stripe, U = self._state()
            theta = self.ctl.theta_from_post(self.ctl.posterior(R, U), rng) if self.cfg.controller in ("coverage", "pivotal") \
                else self.ctl.posterior(R, U) @ self.ctl.grid_theta
            p = float(self.ctl.p_complete(R, rec_pair, rec_frag, by_stripe, theta, U, rng)[0])
            self.cp_rows.append((level, p))

    # ------------------------------------------------------------------ run
    def run(self, rng) -> dict:
        t0 = time.time()
        L = self.layout
        used = L.stored_per_container() > 0
        b = self.ctl.batch
        self._explore_end = int(used.sum()) * self.cfg.explore_rounds
        if self.static_w is None:
            for j in np.flatnonzero(used):
                self.read_container(j, self.ctl.probe)
        rounds, repairs = 0, 0
        switches, last = 0, -1
        status = "partial"
        while self.reads_used.sum() < self.cap_reads:
            if self.tracker.complete():
                res = reconstruct(self.man, self.tracker.rows)
                status = res.status.value
                if status == "exact" or not res.unverified_stripes or repairs > 20:
                    break
                for s_bad in res.unverified_stripes:  # one more fragment to localise the bad one
                    if self.tracker.k[s_bad] < self.tracker.n[s_bad]:
                        self.tracker.k[s_bad] += 1
                repairs += 1
            avail = used & ~self.exhausted()
            if not avail.any():
                break
            self._record_checkpoint(rng)
            j = self.choose(avail, rounds, rng)
            switches += int(last >= 0 and j != last)
            last = j
            self.read_container(j, b)
            rounds += 1
        res = reconstruct(self.man, self.tracker.rows)
        return {"status": res.status.value, "reads": float(self.reads_used.sum()),
                "reads_per_container": self.reads_used.tolist(), "rounds": rounds,
                "decode_attempts": self.decode_attempts, "false_accepts": self.false_accepts,
                "verification_rounds": repairs, "switches": switches, "files": res.files, "seconds": time.time() - t0,
                "checkpoints": self.cp_rows}
