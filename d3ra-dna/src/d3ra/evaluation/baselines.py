"""Baselines A-F (identical code: same K, parity, stripe count, file).

A  uniform (round-robin) placement + fixed sequencing (uniform reads)
B  reliability-weighted placement + fixed sequencing (coverage-equalising reads)
C  uniform placement + adaptive sequencing (pivotal controller)
D  reliability-weighted placement + adaptive sequencing
E  joint placement/redundancy/retrieval optimisation + adaptive sequencing
F  E's design + oracle controller (perfect channel knowledge; reference only)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..channel.surrogate import Layout
from ..optimization.joint_optimizer import _layout, joint_optimise
from ..retrieval.fixed_policy import coverage_equalising_weights, uniform_weights
from ..storage.placement import reliability_weighted, round_robin
from .scenarios import Scenario


@dataclass
class Design:
    name: str
    layout: Layout
    policy: str          # "static" | controller name
    weights: np.ndarray | None = None
    variant: str = ""
    parity: int = 0


def _stripe_of(sc: Scenario, parity: int | None = None):
    p = sc.parity if parity is None else parity
    return np.repeat(np.arange(sc.n_stripes), sc.k + p), np.full(sc.n_stripes, sc.k)


def uniform_layout(sc: Scenario) -> Layout:
    so, k = _stripe_of(sc)
    return Layout(so, k, round_robin(so, len(sc.containers)))


def reliability_layout(sc: Scenario) -> Layout:
    so, k = _stripe_of(sc)
    scores = np.array([c.reliability() for c in sc.containers]) * sc.belief.yield_ / np.maximum(sc.belief.yield_, 1e-9)
    # belief-based reliability: survival x presence x yield
    b = sc.belief
    scores = (1 - b.failure.marginal()) * (1 - b.dropout) * b.yield_ * (b.phi[:, 3])
    return Layout(so, k, reliability_weighted(so, scores))


def joint_design(sc: Scenario, target: float = 0.99, parities=None):
    res = joint_optimise(sc.n_stripes, sc.k, sc.belief, target, sc.domain,
                         parities=parities if parities is not None else [sc.parity], iters=3, p0=sc.parity)
    variant = res.history[-1]["variant"] if res and res.history else "domain"
    return res, variant


def designs(sc: Scenario, target: float = 0.99) -> list[Design]:
    U = uniform_layout(sc)
    RW = reliability_layout(sc)
    res, variant = joint_design(sc, target)
    E_layout = res.layout if res is not None else _layout("tiered1", sc.n_stripes, sc.k, sc.parity, sc.belief, sc.domain)
    return [
        Design("A_uniform_fixed", U, "static", uniform_weights(U)),
        Design("B_reliability_fixed", RW, "static", coverage_equalising_weights(RW, sc.belief)),
        Design("C_uniform_adaptive", U, "pivotal"),
        Design("D_reliability_adaptive", RW, "pivotal"),
        Design("E_joint_adaptive", E_layout, "pivotal", res.weights if res is not None else None, variant, sc.parity),
        Design("F_oracle", E_layout, "oracle", None, variant, sc.parity),
    ]
