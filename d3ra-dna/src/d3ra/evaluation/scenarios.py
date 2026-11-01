"""Experimental scenarios 1-10 (controlled simulation assumptions).

Every parameter here is *hypothetically assigned for a controlled simulation*
unless it is the base error profile, which is taken from published
Illumina-based DNA storage measurements (sub 4e-3, del 2e-3, ins 4e-4 per nt;
Organick et al. 2018, see docs/channel_calibration.md). Consensus-success
curves phi(n) are *estimated by the sequence-level simulator* on a grid of
error scales and interpolated per container.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from ..channel.calibration import measure_phi, phi_from_measurement
from ..channel.failures import FailureModel
from ..channel.sequencer import ReadErrorModel
from ..channel.surrogate import ChannelParams, default_phi
from ..storage.container import Container

BASE_ERR = (4e-3, 4e-4, 2e-3)          # (sub, ins, del) per nt, Illumina-like
ERR_SCALES = np.array([0.5, 1, 2, 3, 4, 6, 8, 12])
_PHI_GRID = None


def phi_grid(trials: int = 40) -> np.ndarray:
    global _PHI_GRID
    if _PHI_GRID is None:
        rows = []
        for s in ERR_SCALES:
            err = ReadErrorModel(BASE_ERR[0] * s, BASE_ERR[1] * s, BASE_ERR[2] * s)
            rows.append(phi_from_measurement(measure_phi(err, trials=trials)))
        _PHI_GRID = np.stack(rows)
    return _PHI_GRID


def phi_for_scale(scale: float, calibrated: bool = True) -> np.ndarray:
    if not calibrated:
        return default_phi(sum(BASE_ERR) * scale)
    G = phi_grid()
    ls = np.log(np.clip(scale, ERR_SCALES[0], ERR_SCALES[-1]))
    lg = np.log(ERR_SCALES)
    i = int(np.clip(np.searchsorted(lg, ls) - 1, 0, len(lg) - 2))
    t = (ls - lg[i]) / (lg[i + 1] - lg[i])
    return np.maximum.accumulate((1 - t) * G[i] + t * G[i + 1])


@dataclass
class Scenario:
    name: str
    description: str
    containers: list
    truth: ChannelParams
    belief: ChannelParams
    n_stripes: int = 25
    k: int = 32
    parity: int = 12
    domain: np.ndarray | None = None
    budget_scale: float = 1.0     # scenario 6 tightens the reference budget
    notes: dict = field(default_factory=dict)


def _containers(spec: dict) -> list[Container]:
    m = len(spec["yield_"])
    out = []
    for j in range(m):
        kw = {k: (v[j] if isinstance(v, (list, np.ndarray)) else v) for k, v in spec.items() if k != "err_scale"}
        s = spec["err_scale"][j]
        kw.update(e_sub=BASE_ERR[0] * s, e_ins=BASE_ERR[1] * s, e_del=BASE_ERR[2] * s)
        out.append(Container(cid=j, name=f"C{j}", **kw))
    return out


def _channel(cs, err_scale, failure=None, sigma_g=0.3, d_syn=0.002, calibrated=True):
    fm = failure or FailureModel(np.array([c.fail_prob for c in cs]))
    phi = np.stack([phi_for_scale(s, calibrated) for s in err_scale])
    return ChannelParams.from_containers(cs, fm, phi, sigma_g=sigma_g, d_syn=d_syn)


HET = dict(
    yield_=[0.95, 0.92, 0.88, 0.85, 0.75, 0.65, 0.55, 0.45],
    abundance_shape=[20, 16, 12, 10, 6, 8, 4, 3],
    err_scale=[1, 1, 1.5, 2, 3, 1.5, 4, 6],
    dropout=[0.005, 0.01, 0.01, 0.015, 0.03, 0.02, 0.05, 0.08],
    fail_prob=[0.005, 0.01, 0.01, 0.02, 0.03, 0.02, 0.05, 0.08],
    cost_per_read=[1.0, 1.0, 0.9, 0.9, 0.8, 0.6, 0.5, 0.5],
)
HOM = dict(yield_=[0.85] * 8, abundance_shape=[10] * 8, err_scale=[1.5] * 8, dropout=[0.015] * 8,
           fail_prob=[0.02] * 8, cost_per_read=[1.0] * 8)


def make(name: str, seed: int = 0, calibrated: bool = True) -> Scenario:
    rng = np.random.default_rng(seed)
    if name == "s1_homogeneous":
        cs = _containers(HOM)
        ch = _channel(cs, HOM["err_scale"], calibrated=calibrated)
        return Scenario(name, "All containers approximately equal", cs, ch, ch)
    if name == "s2_heterogeneous":
        cs = _containers(HET)
        ch = _channel(cs, HET["err_scale"], calibrated=calibrated)
        return Scenario(name, "Different abundance, dropout, error, failure and cost", cs, ch, ch)
    if name == "s3_nonuniform_abundance":
        spec = dict(HET, abundance_shape=[4, 3, 3, 2.5, 2, 2, 1.5, 1.5])
        cs = _containers(spec)
        ch = _channel(cs, spec["err_scale"], sigma_g=0.8, calibrated=calibrated)
        return Scenario(name, "Strong sequence-intrinsic (sigma_g=0.8) and container overdispersion", cs, ch, ch)
    if name == "s4_correlated_failures":
        cs = _containers(HET)
        dom = np.array([0, 0, 1, 1, 2, 2, 3, 3])
        for c, d in zip(cs, dom):
            c.domain = int(d)
        fm = FailureModel(np.array([c.fail_prob for c in cs]) * 0.5, [[0, 1], [2, 3], [4, 5], [6, 7]],
                          [0.04, 0.03, 0.03, 0.05])
        ch = _channel(cs, HET["err_scale"], failure=fm, calibrated=calibrated)
        return Scenario(name, "Pairs of containers share a common-mode failure domain", cs, ch, ch, domain=dom, parity=20)
    if name == "s5_high_errors":
        spec = dict(HET, err_scale=[2, 3, 4, 6, 8, 6, 10, 12])
        cs = _containers(spec)
        ch = _channel(cs, spec["err_scale"], calibrated=calibrated)
        return Scenario(name, "Nanopore-like error rates (up to ~7.7%/nt)", cs, ch, ch)
    if name == "s6_tight_budget":
        sc = make("s2_heterogeneous", seed, calibrated)
        return replace(sc, name=name, description="Heterogeneous, budget = 0.6x reference", budget_scale=0.6)
    if name == "s7_channel_change":
        cs = _containers(HET)
        belief = _channel(cs, HET["err_scale"], calibrated=calibrated)
        spec = dict(HET)
        spec["yield_"] = [0.30, 0.92, 0.30, 0.85, 0.75, 0.65, 0.55, 0.45]   # two cheap-good containers degrade
        spec["err_scale"] = [4, 1, 6, 2, 3, 1.5, 4, 6]
        cs2 = _containers(spec)
        truth = _channel(cs2, spec["err_scale"], calibrated=calibrated)
        return Scenario(name, "Containers 0 and 2 degrade after placement (yield, errors)", cs, truth, belief)
    if name == "s9_uncertain_estimates":
        cs = _containers(HET)
        truth = _channel(cs, HET["err_scale"], calibrated=calibrated)
        belief = truth.perturbed(rng, rel_sd=0.5)
        return Scenario(name, "Optimiser sees log-normally perturbed estimates (sd 0.5)", cs, truth, belief)
    if name == "s10_mixed":
        spec = dict(HET, err_scale=[1.5, 2, 3, 4, 6, 3, 8, 10], dropout=[0.02, 0.03, 0.03, 0.04, 0.06, 0.05, 0.08, 0.1],
                    abundance_shape=[6, 5, 4, 4, 3, 3, 2, 2])
        cs = _containers(spec)
        dom = np.array([0, 1, 2, 3, 0, 1, 2, 3])
        for c, d in zip(cs, dom):
            c.domain = int(d)
        fm = FailureModel(np.array([c.fail_prob for c in cs]), [[0, 4], [1, 5], [2, 6], [3, 7]], [0.03] * 4)
        truth = _channel(cs, spec["err_scale"], failure=fm, sigma_g=0.6, d_syn=0.005, calibrated=calibrated)
        belief = truth.perturbed(rng, rel_sd=0.3)
        return Scenario(name, "Outages + domains + overdispersion + high errors + loss + estimate noise",
                        cs, truth, belief, domain=dom, parity=20)
    raise KeyError(name)


SCENARIOS = ["s1_homogeneous", "s2_heterogeneous", "s3_nonuniform_abundance", "s4_correlated_failures",
             "s5_high_errors", "s6_tight_budget", "s7_channel_change", "s9_uncertain_estimates", "s10_mixed"]


def scale_scenario(m: int, n_stripes: int, k: int = 32, parity: int = 12, seed: int = 0) -> Scenario:
    """Scenario 8: random heterogeneous population of m containers."""
    rng = np.random.default_rng(seed)
    scale = np.exp(rng.normal(np.log(2), 0.5, m))
    spec = dict(yield_=np.clip(rng.normal(0.8, 0.12, m), 0.3, 1.0), abundance_shape=np.exp(rng.normal(np.log(8), 0.5, m)),
                err_scale=scale, dropout=np.clip(np.exp(rng.normal(np.log(0.02), 0.6, m)), 1e-3, 0.2),
                fail_prob=np.clip(np.exp(rng.normal(np.log(0.02), 0.6, m)), 1e-3, 0.2),
                cost_per_read=np.exp(rng.normal(0, 0.3, m)))
    cs = _containers(spec)
    ch = _channel(cs, scale, calibrated=False)
    return Scenario(f"s8_scale_m{m}_s{n_stripes}", "Large-scale", cs, ch, ch, n_stripes=n_stripes, k=k, parity=parity)
