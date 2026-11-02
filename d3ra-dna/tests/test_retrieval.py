"""Clustering, consensus, tracking, stopping rules and sequencing controllers."""
from __future__ import annotations

import numpy as np
import pytest

from d3ra.channel.failures import FailureModel
from d3ra.channel.sequencer import ReadErrorModel, noisy_reads
from d3ra.channel.surrogate import ChannelParams, Layout, sample_T
from d3ra.retrieval.adaptive_policy import THETA, AdaptiveConfig, AdaptiveController
from d3ra.retrieval.consensus import InventoryClusterer, cluster_reads, consensus, majority_vote
from d3ra.retrieval.fixed_policy import (StaticPolicy, budget_for_target, coverage_equalising_weights,
                                         reliability_at, uniform_weights)
from d3ra.retrieval.fragment_tracker import FragmentTracker
from d3ra.retrieval.stopping_rule import (DecodedStop, FixedBudgetStop, FutilityStop, completion_probability,
                                          false_accept_bound)
from d3ra.storage.container import make_heterogeneous
from d3ra.storage.placement import round_robin


def _rand_strand(n, rng):
    return "".join(rng.choice(list("ACGT"), n))


# ----------------------------------------------------------------------------- consensus

def test_consensus_recovers_strand_from_noisy_reads():
    rng = np.random.default_rng(0)
    err = ReadErrorModel()                     # ~0.9% total error per base
    ok = maj = 0
    trials = 60
    for _ in range(trials):
        s = _rand_strand(120, rng)
        reads = noisy_reads(s, 7, err, rng)
        ok += consensus(reads, 120) == s
        maj += majority_vote(reads, 120) == s
    assert ok / trials >= 0.85
    assert ok >= maj                           # indel-aware BMA beats naive voting


def test_consensus_edge_cases():
    assert consensus([], 10) == ""
    assert consensus(["ACG"], 6) == "ACGAAA"
    assert consensus(["ACGTACGTAA"], 4) == "ACGT"
    assert consensus(["ACGT", "ACGT", "ACGT"], 4) == "ACGT"
    # one read with an insertion, one with a deletion, one clean
    truth = "ACGTTGCAACGT"
    reads = [truth, truth[:5] + "A" + truth[5:], truth[:7] + truth[8:]]
    assert consensus(reads, len(truth)) == truth
    assert majority_vote(["AC", "AG", "AG"], 2) == "AG"
    assert majority_vote(["A"], 3) == "AAA"


def test_inventory_clusterer_clean_reads():
    rng = np.random.default_rng(1)
    strands = [_rand_strand(60, rng) for _ in range(40)]
    keys = [f"frag{i}" for i in range(40)]
    P = 16
    cl = InventoryClusterer(keys, [s[:P] for s in strands], P)
    order = rng.permutation(40)
    reads = [strands[i] for i in order]
    assert cl.assign(reads) == [keys[i] for i in order]
    # single substitution in the prefix, one leading insertion, one leading deletion
    r_sub = strands[3][:2] + ("A" if strands[3][2] != "A" else "C") + strands[3][3:]
    r_ins = "G" + strands[5]
    r_del = strands[7][1:]
    assert cl.assign([r_sub, r_ins, r_del]) == ["frag3", "frag5", "frag7"]
    # junk reads are rejected
    junk = [_rand_strand(60, rng) for _ in range(20)]
    assert sum(k is not None for k in cl.assign(junk)) <= 1
    assert cl.assign([]) == []
    assert InventoryClusterer([], [], P).assign(["ACGT"]) == [None]
    groups = cluster_reads(reads + reads[:5], cl)
    assert sum(len(v) for v in groups.values()) == 45
    assert all(len(groups[keys[i]]) == 2 for i in order[:5])


# ----------------------------------------------------------------------------- tracker

def test_fragment_tracker():
    t = FragmentTracker(np.array([2, 3]), np.array([3, 4]))
    assert t.need().tolist() == [2, 3] and not t.complete()
    assert t.add(0, 0, b"a", container=1)
    assert not t.add(0, 0, b"b")                       # duplicate ignored
    assert t.rows[(0, 0)] == b"a" and t.source[(0, 0)] == 1
    t.add(0, 2, b"c")
    for i in range(3):
        t.add(1, i, b"x")
    assert t.complete() and t.need().tolist() == [0, 0]
    assert t.missing() == [(0, 1), (1, 3)]
    assert t.summary() == {"recovered": 5, "stripes_complete": 2, "stripes": 2, "min_margin": 0}


# ----------------------------------------------------------------------------- stopping rules

def test_stopping_rules():
    assert false_accept_bound(10, 2) == pytest.approx(10 / 65536)
    assert false_accept_bound(10 ** 9, 1) == 1.0
    cost = np.array([5.0, 10.0, 20.0])
    done = np.array([False, True, False])
    assert FixedBudgetStop(10)(cost, done).tolist() == [False, True, True]
    assert DecodedStop()(cost, done).tolist() == [False, True, False]
    assert DecodedStop(cap=20)(cost, done).tolist() == [False, True, True]
    fs = FutilityStop(cap=100, eps=0.01)
    assert fs(cost, done).tolist() == [False, True, False]
    assert fs(cost, done, p_complete=np.array([0.001, 0.0, 0.5])).tolist() == [True, True, False]


def test_completion_probability():
    mask = np.ones((1, 2, 4), bool)
    rho = np.ones((1, 2, 4))
    assert completion_probability(rho, np.array([[2, 3]]), mask)[0] == pytest.approx(1.0)
    assert completion_probability(np.zeros((1, 2, 4)), np.array([[1, 0]]), mask)[0] < 1e-6
    assert completion_probability(np.zeros((1, 2, 4)), np.array([[0, 0]]), mask)[0] == 1.0
    lo = completion_probability(np.full((1, 2, 4), 0.3), np.array([[2, 2]]), mask)[0]
    hi = completion_probability(np.full((1, 2, 4), 0.7), np.array([[2, 2]]), mask)[0]
    assert 0 < lo < hi < 1
    # masked-out entries do not count
    m2 = mask.copy()
    m2[..., 2:] = False
    assert completion_probability(rho, np.array([[3, 0]]), m2)[0] < 1e-3


# ----------------------------------------------------------------------------- fixed policy

def test_budget_for_target_and_reliability():
    s = np.arange(1, 101, dtype=float)
    assert budget_for_target(s, 0.99) == 99
    assert budget_for_target(s, 1.0) == 100
    assert budget_for_target(s, 0.5) == 50
    assert budget_for_target(s, 0.0) == 1
    assert reliability_at(s, 99) == pytest.approx(0.99)
    assert reliability_at(s, budget_for_target(s, 0.9)) >= 0.9
    s2 = np.concatenate([s[:95], np.full(5, np.inf)])
    assert budget_for_target(s2, 0.9) == 90
    assert budget_for_target(s2, 0.99) == np.inf


# ----------------------------------------------------------------------------- controllers

CONTROLLERS = ["uniform", "coverage", "pivotal", "thompson", "oracle", "clairvoyant"]


@pytest.fixture(scope="module")
def surrogate_setup():
    cs = make_heterogeneous(4, np.random.default_rng(0))
    ch = ChannelParams.from_containers(cs, failure=FailureModel(np.zeros(4)))
    S, K, N = 4, 8, 12
    stripe_of = np.repeat(np.arange(S), N)
    L = Layout(stripe_of, np.full(S, K), round_robin(stripe_of, 4))
    return L, ch


def test_static_weights(surrogate_setup):
    L, ch = surrogate_setup
    w = uniform_weights(L)
    np.testing.assert_allclose(w, 0.25)
    wc = coverage_equalising_weights(L, ch)
    assert wc.sum() == pytest.approx(1.0) and (wc > 0).all()
    sp = StaticPolicy(np.array([2.0, 2.0, 2.0, 2.0]))
    np.testing.assert_allclose(sp.w, 0.25)
    reads = sp.budget_samples(L, ch, 300, np.random.default_rng(0), cost=False)
    cost = sp.budget_samples(L, ch, 300, np.random.default_rng(0), cost=True)
    np.testing.assert_allclose(cost, reads * float(ch.cost @ sp.w))


def test_all_controllers_succeed(surrogate_setup):
    L, ch = surrogate_setup
    static = StaticPolicy(uniform_weights(L)).budget_samples(L, ch, 3000, np.random.default_rng(3))
    b99 = budget_for_target(static, 0.99)
    assert np.isfinite(b99)
    T, failed = sample_T(L, ch, 200, np.random.default_rng(2), return_state=True)
    results = {}
    for ctl in CONTROLLERS:
        ctrl = AdaptiveController(L, ch, AdaptiveConfig(controller=ctl, batch=10, cap=1e4))
        res = ctrl.run(200, np.random.default_rng(4), T=T, failed=failed)
        assert res.success.all(), ctl
        assert np.isfinite(res.cost).all() and (res.cost > 0).all(), ctl
        assert res.reads_per_container.shape == (200, 4)
        np.testing.assert_allclose(res.reads, res.reads_per_container.sum(1))
        np.testing.assert_allclose(res.cost, res.reads_per_container @ ch.cost)   # setup_cost == 0
        results[ctl] = res.cost.mean()
        if ctl != "uniform" and ctl != "clairvoyant":
            assert results[ctl] <= b99, ctl
    assert results["uniform"] <= b99
    # adaptive stopping alone beats the 0.99 static budget; adaptive allocation beats round-robin
    assert results["pivotal"] <= results["uniform"]
    assert results["coverage"] <= results["uniform"]
    # knowing the realised T is the best reference
    assert results["clairvoyant"] <= min(v for k, v in results.items() if k != "clairvoyant")


def test_controller_probe_and_batch_defaults(surrogate_setup):
    L, ch = surrogate_setup
    ctrl = AdaptiveController(L, ch, AdaptiveConfig())
    assert ctrl.batch == 50 and ctrl.probe == 50
    res = ctrl.run(20, np.random.default_rng(0))
    assert (res.reads_per_container >= 50).all()           # every container probed once
    assert (res.containers_used == 4).all()
    assert res.meta["batch"] == 50


def test_controller_respects_cap(surrogate_setup):
    L, ch = surrogate_setup
    ctrl = AdaptiveController(L, ch, AdaptiveConfig(controller="pivotal", batch=10, cap=30.0))
    res = ctrl.run(30, np.random.default_rng(0))
    assert not res.success.any()                           # far too little budget
    # stops right after crossing the cap: probe (4 containers) + at most one extra batch
    assert (res.cost <= max(30.0, 4 * 10 * ch.cost.max()) + 10 * ch.cost.max() + 1e-9).all()


def test_posterior_is_distribution_and_detects_failure(surrogate_setup):
    L, ch = surrogate_setup
    ctrl = AdaptiveController(L, ch, AdaptiveConfig(controller="pivotal"))
    R = np.full((2, 4), 300.0)
    U = np.array([[12, 12, 12, 12], [12, 0, 12, 12]], dtype=float)
    post = ctrl.posterior(R, U)
    assert post.shape == (2, 4, len(THETA))
    np.testing.assert_allclose(post.sum(2), 1.0)
    assert post[0, :, 0].max() < 1e-6                      # recovered fragments rule out failure
    mean = post @ THETA
    assert mean[1, 1] < 0.3 * mean[0, 1]                   # silent container looks (near-)dead


def test_overdispersed_posterior_is_wider(surrogate_setup):
    L, ch = surrogate_setup
    R = np.full((1, 4), 300.0)
    U = np.array([[12, 12, 12, 12]], dtype=float)
    sd = []
    for od in (0.0, 0.05):
        post = AdaptiveController(L, ch, AdaptiveConfig(controller="pivotal", overdispersion=od)).posterior(R, U)
        np.testing.assert_allclose(post.sum(2), 1.0)
        lt = np.log(THETA[1:])
        p = post[0, 0, 1:] / post[0, 0, 1:].sum()
        sd.append(np.sqrt(p @ lt ** 2 - (p @ lt) ** 2))
    assert sd[1] > sd[0]                                   # beta-binomial spreads the posterior


@pytest.mark.parametrize("ctl_name", ["ts", "ucb", "ltc", "info", "tail", "tail_lcb"])
def test_new_controllers_complete(surrogate_setup, ctl_name):
    L, ch = surrogate_setup
    ctl = AdaptiveController(L, ch, AdaptiveConfig(controller=ctl_name, explore_rounds=2))
    r = ctl.run(20, np.random.default_rng(1))
    assert r.success.mean() > 0.8
    assert np.all(np.isfinite(r.cost[r.success]))


def test_mc_forecast_is_probability(surrogate_setup):
    L, ch = surrogate_setup
    base = AdaptiveController(L, ch, AdaptiveConfig(controller="coverage")).run(10, np.random.default_rng(2))
    cap = float(np.median(base.cost))
    ctl = AdaptiveController(L, ch, AdaptiveConfig(controller="coverage", cap=cap, forecast="mc", forecast_draws=50))
    r = ctl.run(10, np.random.default_rng(2), checkpoints=[0.6 * cap])
    vals = [v for _, _, v in r.meta["checkpoints"]]
    assert vals and all(0.0 <= v <= 1.0 for v in vals)
