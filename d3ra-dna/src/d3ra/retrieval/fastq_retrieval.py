"""Streaming retrieval from real or simulated FASTQ files (one per container).

The retriever reads each container's FASTQ in batches (as a sequencer would
deliver them), clusters reads against the expected address prefixes of that
container's inventory (computed from the manifest), runs BMA consensus and
inner RS+CRC decoding, and stops as soon as every stripe holds K verified
fragments; the adaptive controller chooses which container to read next.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np

from ..archive import address_prefix, load_archive
from ..channel.failures import FailureModel
from ..channel.surrogate import ChannelParams, Layout
from ..decoding.reconstruction import reconstruct
from ..io.fastq import read_fastq
from .adaptive_policy import THETA, AdaptiveConfig, AdaptiveController
from .consensus import InventoryClusterer, consensus
from .fragment_tracker import FragmentTracker


class FastqRetriever:
    def __init__(self, archive, fastqs: list, controller: str = "pivotal", batch: int | None = None,
                 max_cluster_reads: int = 24, cap_reads: float = np.inf, consensus_mode: str = "two_sided",
                 adapter: str | None = "AGATCGGAAGAGC"):
        """consensus_mode: 'two_sided' when reads end at the strand end (primer-trimmed or
        simulated reads); 'forward' for raw reads with an uncertain 3' tail."""
        self.two_sided = consensus_mode == "two_sided"
        self.adapter = adapter
        self.man, self.fmt, self.placement, self.containers, self.addresses = load_archive(archive)
        self.fastqs = [Path(f) if f else None for f in fastqs]
        m = self.placement.shape[1]
        assert len(self.fastqs) == m, "one FASTQ (or None) per container"
        self.iters = [read_fastq(f) if f is not None and f.exists() else iter(()) for f in self.fastqs]
        self.exhausted = np.array([f is None or not f.exists() for f in self.fastqs])
        self.clusterers = []
        for j in range(m):
            idx = np.flatnonzero(self.placement[:, j] > 0)
            self.clusterers.append(InventoryClusterer(list(idx), [address_prefix(self.fmt, self.addresses[i]) for i in idx],
                                                      self.fmt.prefix_nt))
        k = np.array([s["k"] for s in self.man.stripes])
        n = np.array([s["n"] for s in self.man.stripes])
        self.tracker = FragmentTracker(k, n)
        self.stripe_of = np.repeat(np.arange(len(k)), n)
        self.layout = Layout(self.stripe_of, k, self.placement)
        ch = ChannelParams.from_containers(self.containers, FailureModel(np.array([c.fail_prob for c in self.containers])),
                                           L=self.fmt.strand_nt, nsym=self.fmt.inner.nsym)
        self.ctl = AdaptiveController(self.layout, ch, AdaptiveConfig(controller=controller, batch=batch or 0))
        self.ch = ch
        self.controller = controller
        self.pools: dict = {}
        self.rec_pair: dict = {}
        self.reads_used = np.zeros(m)
        self.max_cluster_reads = max_cluster_reads
        self.cap_reads = cap_reads
        self.log = []

    def _next_batch(self, j: int, n: int) -> list[str]:
        out = []
        for _ in range(n):
            try:
                seq = next(self.iters[j]).seq
                if self.adapter:
                    k = seq.find(self.adapter)
                    seq = seq[:k] if k >= 0 else seq
                out.append(seq)
            except StopIteration:
                self.exhausted[j] = True
                break
        return out

    def read_container(self, j: int, n: int) -> int:
        reads = self._next_batch(j, n)
        self.reads_used[j] += len(reads)
        keys = self.clusterers[j].assign(reads)
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
            if (addr[1], addr[2]) in self.tracker.rows:
                self.rec_pair[(i, j)] = True
                continue
            dec = self.fmt.decode(consensus(self.pools[(j, i)], L, two_sided=self.two_sided))
            if dec is None or dec[0] != addr:
                continue
            self.rec_pair[(i, j)] = True
            new += self.tracker.add(addr[1], addr[2], dec[1], j)
        self.log.append((j, len(reads), new))
        return new

    def run(self) -> dict:
        t0 = time.time()
        L = self.layout
        used = (L.stored_per_container() > 0) & ~self.exhausted
        for j in np.flatnonzero(used):
            self.read_container(j, self.ctl.probe)
        rr = 0
        uid = {(a[1], a[2]): i for i, a in enumerate(self.addresses)}
        while not self.tracker.complete() and self.reads_used.sum() < self.cap_reads:
            avail = used & ~self.exhausted
            if not avail.any():
                break
            if self.controller == "uniform":
                cand = np.flatnonzero(avail)
                j = cand[rr % len(cand)]
                rr += 1
            else:
                R = self.reads_used[None, :]
                rec_pair = np.array([self.rec_pair.get((i, jj), False) for i, jj in zip(L.pair_frag, L.pair_cont)])[None]
                rec_frag = np.zeros((1, L.n), bool)
                for key in self.tracker.rows:
                    rec_frag[0, uid[key]] = True
                U = L.pair_to_cont_sum(rec_pair.astype(float))
                theta = self.ctl.posterior(R, U) @ THETA
                sc = self.ctl.scores(R, rec_pair, rec_frag, self.tracker.by_stripe[None, :], theta, None)
                sc[:, ~avail] = -np.inf
                j = int(sc.argmax())
            self.read_container(j, self.ctl.batch)
        res = reconstruct(self.man, self.tracker.rows)
        return {"status": res.status.value, "reads_per_container": self.reads_used.tolist(),
                "reads_total": float(self.reads_used.sum()), "actions": len(self.log), "result": res,
                "stripes_complete": int((self.tracker.by_stripe >= self.tracker.k).sum()),
                "stripes": len(self.tracker.k), "seconds": time.time() - t0}


def simulate_fastq(archive, out_dir, reads_per_container, error_models=None, seed: int = 0,
                   failed=None, sigma_g: float = 0.3, d_syn: float = 0.002, model: str = "multinomial") -> list[Path]:
    """Write one FASTQ per container from the archive's synthesis files."""
    from ..channel.abundance import container_variation, intrinsic_bias
    from ..channel.dropout import sample_dropout
    from ..channel.sequencer import ReadErrorModel, sequence_container
    from ..io.fastq import FastqRecord, phred_string, read_fasta, write_fastq
    rng = np.random.default_rng(seed)
    man, fmt, placement, containers, addresses = load_archive(archive)
    n, m = placement.shape
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    strands_all = [None] * n
    for j in range(m):
        names, seqs = read_fasta(Path(archive) / "synthesis" / f"container_{j:02d}.fasta")
        for nm, sq in zip(names, seqs):
            strands_all[int(nm.split("_")[0][4:])] = sq
    g = intrinsic_bias(None, n, sigma_g, 0, rng)
    h = container_variation(n, m, np.array([c.abundance_shape for c in containers]), rng)
    present = sample_dropout(n, m, np.array([c.dropout for c in containers]), d_syn, rng)
    failed = np.zeros(m, bool) if failed is None else np.asarray(failed, bool)
    rpc = np.broadcast_to(np.asarray(reads_per_container), (m,))
    paths = []
    for j, c in enumerate(containers):
        idx = np.flatnonzero(placement[:, j] > 0)
        if idx.size == 0:
            paths.append(None)
            continue
        err = (error_models or {}).get(j) or ReadErrorModel.from_container(c)
        p_err = getattr(err, "mean_errors", None)
        p_err = (p_err / fmt.strand_nt) if p_err is not None else c.error_rate
        reads = sequence_container([strands_all[i] for i in idx], g[idx] * h[idx, j] * placement[idx, j],
                                   present[idx, j], int(rpc[j]), err, rng, model, c.yield_, failed=bool(failed[j]))
        path = out / f"container_{j:02d}.fq.gz"
        write_fastq(path, (FastqRecord(f"c{j}_r{k}_src{src}", r, phred_string(len(r), p_err)) for k, (src, r) in enumerate(reads)))
        paths.append(path)
    return paths
