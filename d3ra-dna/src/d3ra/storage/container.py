"""DNA storage container (vial / compartment) as a simulated storage object."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Container:
    cid: int
    name: str = ""
    domain: int = -1                 # failure-domain label
    capacity: int = 10**9            # max number of distinct stored sequences
    read_share: float = 1.0          # p_j: share of a pooled run (pooled mode only)
    yield_: float = 1.0              # fraction of reads that are usable (map to a stored strand)
    abundance_shape: float = 10.0    # gamma shape of container-specific abundance variation
    e_sub: float = 0.005
    e_ins: float = 0.001
    e_del: float = 0.003
    homopolymer_boost: float = 1.5   # indel multiplier per extra base in a run
    quality_sigma: float = 0.3       # log-normal read-quality variation
    burst_rate: float = 0.0          # correlated burst-deletion start rate
    dropout: float = 0.01            # d_j: per-molecule absence probability
    fail_prob: float = 0.01          # f_j: independent container failure
    cost_per_read: float = 1.0       # c_j (arbitrary units)
    setup_cost: float = 0.0          # fixed cost per sequencing run
    latency_per_run: float = 1.0     # time units per sequencing batch
    inventory: list = field(default_factory=list)   # fragment uids stored here
    history: list = field(default_factory=list)     # (reads, new_recovered) per batch

    @property
    def error_rate(self) -> float:
        return self.e_sub + self.e_ins + self.e_del

    def reliability(self, coverage: float = 10.0) -> float:
        """Heuristic single-number reliability score used by naive strategies:
        survival x presence x per-read cleanliness at a nominal coverage."""
        clean = (1 - self.error_rate) ** 100
        return (1 - self.fail_prob) * (1 - self.dropout) * self.yield_ * (0.5 + 0.5 * clean)

    def store(self, uid) -> None:
        if len(self.inventory) >= self.capacity:
            raise OverflowError(f"container {self.cid} full")
        self.inventory.append(uid)

    def record(self, reads: int, new: int) -> None:
        self.history.append((reads, new))

    def as_dict(self) -> dict:
        d = {k: v for k, v in self.__dict__.items() if k not in ("inventory", "history")}
        return d


def containers_from_config(cfg: list[dict]) -> list[Container]:
    return [Container(cid=i, **{k: v for k, v in c.items() if k != "cid"}) for i, c in enumerate(cfg)]


def make_heterogeneous(m: int, rng: np.random.Generator, spread: float = 1.0,
                       n_domains: int = 0) -> list[Container]:
    """Random heterogeneous container population. ``spread``=0 -> homogeneous."""
    out = []
    for j in range(m):
        z = rng.normal(0, 1, 6) * spread
        out.append(Container(
            cid=j, name=f"C{j:03d}",
            domain=(j % n_domains) if n_domains else -1,
            yield_=float(np.clip(0.85 * np.exp(0.15 * z[0]), 0.2, 1.0)),
            abundance_shape=float(np.clip(8 * np.exp(0.6 * z[1]), 1.0, 100.0)),
            e_sub=float(np.clip(0.004 * np.exp(0.5 * z[2]), 1e-4, 0.05)),
            e_ins=float(np.clip(0.001 * np.exp(0.5 * z[2]), 1e-5, 0.02)),
            e_del=float(np.clip(0.003 * np.exp(0.5 * z[2]), 1e-4, 0.05)),
            dropout=float(np.clip(0.01 * np.exp(0.8 * z[3]), 1e-4, 0.3)),
            fail_prob=float(np.clip(0.01 * np.exp(0.8 * z[4]), 1e-4, 0.3)),
            cost_per_read=float(np.clip(np.exp(0.3 * z[5]), 0.2, 5.0)),
        ))
    return out
