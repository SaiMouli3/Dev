"""Placement strategies and the placement optimiser."""
from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from d3ra.channel.failures import FailureModel
from d3ra.channel.surrogate import ChannelParams
from d3ra.optimization.placement_optimizer import (PlacementProblem, counts_to_assign, from_matrix, genetic,
                                                   greedy, ilp, local_search, simulated_annealing, to_matrix)
from d3ra.storage.container import make_heterogeneous
from d3ra.storage.placement import (add_replicas, cost_tiered, domain_aware, failure_aware, random_placement,
                                    reliability_weighted, round_robin, stripe_failure_exposure,
                                    unit_retrieval_cost)


def _stripes(S, N):
    return np.repeat(np.arange(S), N)


def _assert_valid(X, n, m, capacity=None):
    assert X.shape == (n, m)
    assert (X >= 0).all()
    assert ((X > 0).sum(1) >= 1).all(), "every fragment must be placed"
    if capacity is not None:
        assert ((X > 0).sum(0) <= np.asarray(capacity)).all()


def _per_stripe_counts(X, stripe_of):
    S = stripe_of.max() + 1
    C = np.zeros((S, X.shape[1]))
    np.add.at(C, stripe_of, (X > 0).astype(float))
    return C


@pytest.fixture
def channel():
    cs = make_heterogeneous(6, np.random.default_rng(0))
    return ChannelParams.from_containers(cs)


# ----------------------------------------------------------------------------- strategies

def test_unit_retrieval_cost_positive_and_monotone(channel):
    e = unit_retrieval_cost(channel)
    assert e.shape == (6,) and np.isfinite(e).all() and (e > 0).all()
    ch2 = ChannelParams(channel.yield_, channel.shape, channel.dropout, channel.cost * 2, channel.phi,
                        channel.failure)
    np.testing.assert_allclose(unit_retrieval_cost(ch2), 2 * e)


def test_all_strategies_valid(channel):
    S, N, K, m = 3, 10, 6, 6
    stripe_of = _stripes(S, N)
    n = len(stripe_of)
    rng = np.random.default_rng(1)
    cap = np.full(m, 8)
    dom = np.array([0, 0, 1, 1, 2, 2])
    e = unit_retrieval_cost(channel)
    scores = np.array([c for c in np.linspace(0.5, 1.0, m)])
    mats = {
        "random": random_placement(stripe_of, m, rng, capacity=cap),
        "round_robin": round_robin(stripe_of, m, capacity=cap),
        "reliability": reliability_weighted(stripe_of, scores),
        "domain": domain_aware(stripe_of, np.full(S, K), dom, m, capacity=cap),
        "tiered": cost_tiered(stripe_of, K, e, dom, capacity=cap),
    }
    for name, X in mats.items():
        _assert_valid(X, n, m, cap if name != "reliability" else None)
        assert ((X > 0).sum(1) == 1).all(), name     # single copy each


def test_round_robin_layout_and_capacity_error():
    stripe_of = _stripes(2, 4)
    X = round_robin(stripe_of, 3)
    assert from_matrix(X).tolist() == [0, 1, 2, 0, 1, 2, 0, 1]   # stripe 1 offset by 1
    X0 = round_robin(stripe_of, 3, offset_per_stripe=False)
    assert from_matrix(X0).tolist() == [0, 1, 2, 0, 0, 1, 2, 0]
    with pytest.raises(OverflowError):
        round_robin(stripe_of, 3, capacity=[1, 1, 1])


def test_reliability_weighted_concentrates():
    stripe_of = _stripes(1, 12)
    X = reliability_weighted(stripe_of, np.array([3.0, 1.0, 0.0, 0.0]))
    np.testing.assert_array_equal((X > 0).sum(0), [9, 3, 0, 0])


@settings(max_examples=60, deadline=None)
@given(S=st.integers(1, 4), K=st.integers(1, 10), P=st.integers(1, 6), m=st.integers(2, 8),
       seed=st.integers(0, 10_000))
def test_cost_tiered_cap_n_minus_k(S, K, P, m, seed):
    N = K + P
    stripe_of = _stripes(S, N)
    e = np.random.default_rng(seed).uniform(0.5, 5, m)
    X = cost_tiered(stripe_of, K, e)
    _assert_valid(X, len(stripe_of), m)
    C = _per_stripe_counts(X, stripe_of)
    assert (C.sum(1) == N).all()
    if m * P >= N:                                    # feasible: single-container failure survivable
        assert C.max() <= P
        states = np.eye(m, dtype=bool)
        assert stripe_failure_exposure(X, stripe_of, np.full(S, K), states, np.ones(m)) == 0.0


def test_cost_tiered_fills_cheapest_first():
    stripe_of = _stripes(2, 12)
    e = np.array([5.0, 1.0, 3.0, 2.0, 4.0])
    X = cost_tiered(stripe_of, 8, e)
    C = _per_stripe_counts(X, stripe_of)
    # cap = 4: cheapest three containers (1, 3, 2) hold everything
    np.testing.assert_array_equal(C, [[0, 4, 4, 4, 0]] * 2)
    # slack lowers the cap and spreads further
    C2 = _per_stripe_counts(cost_tiered(stripe_of, 8, e, slack=1), stripe_of)
    assert C2.max() <= 3 and (C2 > 0).sum(1).min() == 4


def test_cost_tiered_respects_domains_and_capacity():
    stripe_of = _stripes(2, 12)
    e = np.array([1.0, 1.1, 1.2, 5.0, 6.0, 7.0])
    dom = np.array([0, 0, 0, 1, 2, 3])
    X = cost_tiered(stripe_of, 8, e, dom)
    C = _per_stripe_counts(X, stripe_of)
    for s in range(2):
        per_dom = np.bincount(dom, weights=C[s])
        assert per_dom.max() <= 4
    cap = np.full(6, 5)
    X = cost_tiered(stripe_of, 8, e, dom, capacity=cap)
    _assert_valid(X, 24, 6, cap)
    with pytest.raises(OverflowError):
        cost_tiered(stripe_of, 8, e, capacity=np.full(6, 3))       # 18 < 24 slots


def test_domain_aware_tolerates_single_domain_loss():
    S, N, K, m = 4, 12, 8, 8
    stripe_of = _stripes(S, N)
    dom = np.repeat(np.arange(4), 2)         # 4 domains x 2 containers
    X = domain_aware(stripe_of, np.full(S, K), dom, m)
    _assert_valid(X, S * N, m)
    C = _per_stripe_counts(X, stripe_of)
    for s in range(S):
        assert np.bincount(dom, weights=C[s]).max() <= N - K
    states = np.stack([dom == d for d in range(4)])
    assert stripe_failure_exposure(X, stripe_of, np.full(S, K), states, np.ones(4)) == 0.0
    # whereas losing two domains is fatal for some stripe
    two = np.stack([(dom == 0) | (dom == 1)])
    assert stripe_failure_exposure(X, stripe_of, np.full(S, K), two, [1.0]) == 1.0


def test_domain_aware_without_domains_spreads_evenly():
    stripe_of = _stripes(3, 10)
    X = domain_aware(stripe_of, np.full(3, 6), np.full(5, -1), 5)
    C = _per_stripe_counts(X, stripe_of)
    assert (C == 2).all()


def test_add_replicas():
    X = np.zeros((3, 4))
    X[[0, 1, 2], [0, 1, 2]] = 1
    Y = add_replicas(X, np.array([0, 2]), np.array([4.0, 3.0, 2.0, 1.0]), copies=2)
    assert Y[0].tolist() == [1, 0, 1, 1]
    assert Y[2].tolist() == [0, 1, 1, 1]
    assert Y[1].tolist() == [0, 1, 0, 0]
    assert X.sum() == 3          # input untouched


def test_stripe_failure_exposure_values():
    stripe_of = _stripes(2, 3)
    X = np.zeros((6, 3))
    X[np.arange(6), [0, 1, 2, 0, 0, 1]] = 1
    k = np.array([2, 2])
    states = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=bool)
    # stripe 1 holds 2 fragments in container 0 -> state 1 is fatal; others survive
    w = np.array([0.7, 0.1, 0.1, 0.1])
    assert stripe_failure_exposure(X, stripe_of, k, states, w) == pytest.approx(0.1)
    # replicas protect: copy stripe-1 fragment 3 into container 2
    X2 = X.copy()
    X2[3, 2] = 1
    assert stripe_failure_exposure(X2, stripe_of, k, states, w) == pytest.approx(0.0)


def test_failure_aware_avoids_risky_domain():
    S, N, K, m = 2, 8, 5, 4
    stripe_of = _stripes(S, N)
    # containers 0,1 share a failure domain that fails often; 2,3 are independent and safe
    fm = FailureModel(np.array([0.01, 0.01, 0.01, 0.01]), domains=[[0, 1]], q=[0.3])
    X = failure_aware(stripe_of, K, fm, presence=np.full(m, 0.99))
    _assert_valid(X, S * N, m)
    C = _per_stripe_counts(X, stripe_of)
    assert (C.sum(1) == N).all()
    assert (C[:, :2].sum(1) <= N - K).all()        # risky domain alone can't sink a stripe
    states = np.stack([s for s, _ in fm.enumerate_states()])
    weights = np.array([p for _, p in fm.enumerate_states()])
    rr = round_robin(stripe_of, m)
    assert stripe_failure_exposure(X, stripe_of, np.full(S, K), states, weights) <= \
        stripe_failure_exposure(rr, stripe_of, np.full(S, K), states, weights) + 1e-12


# ----------------------------------------------------------------------------- optimiser

@pytest.fixture
def tiny_problem():
    S, N, K, m = 2, 6, 4, 4
    stripe_of = _stripes(S, N)
    fm = FailureModel(np.array([0.05, 0.02, 0.1, 0.03]))
    ex = fm.enumerate_states()
    states = np.stack([s for s, _ in ex]).astype(int)
    weights = np.array([p for _, p in ex])
    unit_cost = np.array([1.0, 2.0, 1.5, 4.0])
    return PlacementProblem(stripe_of, np.full(S, K), unit_cost, states, weights, eps=1e-2)


def test_placement_problem_objective(tiny_problem):
    prob = tiny_problem
    assert prob.m == 4 and prob.c_ref() == pytest.approx(8.0)
    # all of stripe s in container 0: cheap but any failure of container 0 is fatal
    a = np.zeros(12, dtype=int)
    C = prob.counts(a)
    np.testing.assert_array_equal(C, [[6, 0, 0, 0], [6, 0, 0, 0]])
    f0 = prob.weights[prob.states[:, 0] == 1].sum() / prob.weights.sum()
    assert prob.objective(a) == pytest.approx(8.0 / 8.0 + f0 / prob.eps)
    # objective(assign) == objective_counts(counts)
    b = np.array([0, 1, 2, 3, 0, 1, 2, 0, 3, 1, 0, 2])
    assert prob.objective(b) == pytest.approx(prob.objective_counts(prob.counts(b)))
    np.testing.assert_array_equal(prob.counts(counts_to_assign(prob.stripe_of, prob.counts(b))), prob.counts(b))
    np.testing.assert_array_equal(from_matrix(to_matrix(b, 4)), b)
    # capacity penalty
    capped = PlacementProblem(prob.stripe_of, prob.k, prob.unit_cost, prob.states, prob.weights,
                              capacity=np.array([3, 12, 12, 12]), eps=prob.eps)
    assert capped.objective(a) == pytest.approx(prob.objective(a) + 1e3 * 9)


def test_ilp_not_worse_than_heuristics(tiny_problem):
    prob = tiny_problem
    a_ilp, info = ilp(prob, time_limit=20)
    assert info["status"] == 0
    j_ilp = prob.objective(a_ilp)
    assert j_ilp == pytest.approx(info["objective"], rel=1e-6, abs=1e-6)
    assert ((a_ilp >= 0) & (a_ilp < 4)).all()
    np.testing.assert_array_equal(prob.counts(a_ilp).sum(1), [6, 6])
    g = greedy(prob)
    rng = np.random.default_rng(0)
    cands = {
        "greedy": g,
        "local_search": local_search(prob, g, iters=500, rng=rng),
        "annealing": simulated_annealing(prob, g, iters=1000, rng=rng),
        "genetic": genetic(prob, pop=10, gens=10, rng=rng, seed_assign=g),
        "round_robin": np.array([(i + s) % 4 for i, s in enumerate(prob.stripe_of)]),
    }
    for name, a in cands.items():
        assert j_ilp <= prob.objective(a) + 1e-9, name
    # local search never worsens its starting point
    assert prob.objective(cands["local_search"]) <= prob.objective(g) + 1e-12


def test_greedy_matches_cost_tiered(tiny_problem):
    prob = tiny_problem
    g = greedy(prob)
    X = cost_tiered(prob.stripe_of, prob.k, prob.unit_cost)
    np.testing.assert_array_equal(g, from_matrix(X))
