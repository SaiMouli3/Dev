"""Channel models: sampling, IDS errors, failures and the fragment-level surrogate."""
from __future__ import annotations

import itertools
from dataclasses import replace

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from d3ra.channel.deletions import burst_delete, delete
from d3ra.channel.failures import FailureModel
from d3ra.channel.insertions import homopolymer_factor, insert
from d3ra.channel.sampling import apply_duplicates, apply_missing_reads, sample_counts
from d3ra.channel.sequencer import ReadErrorModel, noisy_reads, sequence_container, to_array, to_str
from d3ra.channel.substitutions import position_profile, substitute
from d3ra.channel.surrogate import (ChannelParams, Layout, Predictor, default_phi, file_time, fragment_times,
                                    pivotal_weights, poisson_binomial_tail, sample_T, sample_thresholds,
                                    static_budget_samples)
from d3ra.retrieval.fixed_policy import uniform_weights
from d3ra.storage.container import Container, make_heterogeneous
from d3ra.storage.placement import round_robin

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


def _rand_seq(n, rng):
    return BASES[rng.integers(0, 4, n)].copy()


def _binom_tol(p, n, z=4.0):
    """z standard errors of an empirical Bernoulli rate over n trials."""
    return z * np.sqrt(p * (1 - p) / n)


# ----------------------------------------------------------------------------- sampling

@pytest.mark.parametrize("model", ["poisson", "multinomial", "negbin"])
def test_sample_counts_means(model):
    rng = np.random.default_rng(0)
    a = np.array([1.0, 2.0, 3.0, 4.0])
    R = 1000
    draws = np.stack([sample_counts(a, R, model, rng) for _ in range(2000)])
    np.testing.assert_allclose(draws.mean(0), R * a / a.sum(), rtol=0.02)
    if model == "multinomial":
        assert (draws.sum(1) == R).all()
    if model == "negbin":   # overdispersed relative to Poisson: var = mu + mu^2/k
        mu = R * a / a.sum()
        np.testing.assert_allclose(draws.var(0), mu + mu ** 2 / 10, rtol=0.15)


def test_sample_counts_edge_cases():
    rng = np.random.default_rng(0)
    assert sample_counts(np.ones(3), 0, "poisson", rng).tolist() == [0, 0, 0]
    assert sample_counts(np.zeros(3), 10, "poisson", rng).tolist() == [0, 0, 0]
    with pytest.raises(ValueError):
        sample_counts(np.ones(3), 10, "zipf", rng)


def test_duplicates_and_missing_reads():
    rng = np.random.default_rng(1)
    c = rng.poisson(20, 5000)
    tot, ind = apply_duplicates(c, 0.3, rng)
    assert (tot == c).all() and (ind <= c).all()
    assert ((c > 0) <= (ind > 0)).all()                    # a sampled molecule keeps >= 1 read
    assert ind.sum() / c.sum() == pytest.approx(0.7, abs=0.02)
    same, same2 = apply_duplicates(c, 0.0, rng)
    assert same is c and same2 is c
    kept = apply_missing_reads(c, 0.25, rng)
    assert (kept <= c).all() and kept.sum() / c.sum() == pytest.approx(0.75, abs=0.02)
    assert apply_missing_reads(c, 0.0, rng) is c


# ----------------------------------------------------------------------------- IDS channels

@pytest.mark.parametrize("p", [0.01, 0.05, 0.2])
def test_substitution_rate(p):
    rng = np.random.default_rng(2)
    seq = _rand_seq(100_000, rng)
    out = substitute(seq, p, rng)
    assert out.size == seq.size
    rate = (out != seq).mean()
    assert rate == pytest.approx(p, abs=_binom_tol(p, seq.size))
    assert set(np.unique(out).tolist()) <= set(BASES.tolist())
    # replacements are uniform over the 3 other bases
    hit = out != seq
    pairs = np.unique((np.searchsorted(BASES, seq[hit]) - np.searchsorted(BASES, out[hit])) % 4,
                      return_counts=True)
    assert pairs[0].tolist() == [1, 2, 3]
    np.testing.assert_allclose(pairs[1] / hit.sum(), 1 / 3, atol=_binom_tol(1 / 3, hit.sum()))
    assert (seq == substitute(seq, 0.0, rng)).all()


def test_substitute_does_not_mutate_input():
    rng = np.random.default_rng(0)
    seq = _rand_seq(100, rng)
    copy = seq.copy()
    substitute(seq, 0.5, rng)
    assert (seq == copy).all()


@pytest.mark.parametrize("p", [0.01, 0.05, 0.2])
def test_deletion_rate(p):
    rng = np.random.default_rng(3)
    seq = _rand_seq(100_000, rng)
    out = delete(seq, p, rng)
    assert 1 - out.size / seq.size == pytest.approx(p, abs=_binom_tol(p, seq.size))


@pytest.mark.parametrize("p", [0.01, 0.05, 0.2])
def test_insertion_rate(p):
    rng = np.random.default_rng(4)
    seq = _rand_seq(100_000, rng)
    out = insert(seq, p, rng)
    assert out.size / seq.size - 1 == pytest.approx(p, abs=_binom_tol(p, seq.size))
    assert insert(seq, 0.0, rng) is seq


def test_burst_deletions():
    rng = np.random.default_rng(5)
    seq = _rand_seq(100_000, rng)
    assert burst_delete(seq, 0.0, 3.0, rng) is seq
    out = burst_delete(seq, 0.01, 3.0, rng)
    # ~1% burst starts x mean length 3, minus overlaps
    assert 0.02 < 1 - out.size / seq.size < 0.035


def test_homopolymer_factor_and_profile():
    seq = to_array("AAACGG")
    np.testing.assert_allclose(homopolymer_factor(seq, 2.0), [1, 2, 4, 1, 1, 2])
    np.testing.assert_allclose(homopolymer_factor(seq, 1.0), np.ones(6))
    prof = position_profile(5, 0.5)
    assert prof[0] == 1.0 and prof[-1] == 1.5 and np.all(np.diff(prof) > 0)


def test_read_error_model_total_rate():
    err = ReadErrorModel(e_sub=0.01, e_ins=0.0, e_del=0.0, quality_sigma=0.0, position_slope=0.0)
    rng = np.random.default_rng(6)
    strand = to_str(_rand_seq(200, rng))
    reads = noisy_reads(strand, 500, err, rng)
    assert all(len(r) == 200 for r in reads)
    rate = np.mean([sum(a != b for a, b in zip(r, strand)) for r in reads]) / 200
    assert rate == pytest.approx(0.01, rel=0.1)
    c = Container(0, e_sub=0.02, e_ins=0.003, e_del=0.004, burst_rate=0.001)
    m = ReadErrorModel.from_container(c)
    assert (m.e_sub, m.e_ins, m.e_del, m.burst_rate) == (0.02, 0.003, 0.004, 0.001)


def test_sequence_container_labels_and_yield():
    rng = np.random.default_rng(7)
    strands = [to_str(_rand_seq(60, rng)) for _ in range(5)]
    err = ReadErrorModel(0, 0, 0, quality_sigma=0)
    out = sequence_container(strands, np.ones(5), np.array([1, 1, 0, 1, 1], bool), 2000, err, rng,
                             yield_=0.8)
    assert len(out) == 2000
    labels = np.array([i for i, _ in out])
    assert 2 not in labels                                     # absent molecule never read
    assert (labels == -1).mean() == pytest.approx(0.2, abs=0.03)
    for i, r in out:
        if i >= 0:
            assert r == strands[i]                             # error-free channel
    junk = sequence_container(strands, np.ones(5), np.ones(5, bool), 50, err, rng, failed=True)
    assert len(junk) == 50 and all(i == -1 for i, _ in junk)


# ----------------------------------------------------------------------------- failures

@pytest.fixture
def fmodel():
    return FailureModel(np.array([0.1, 0.05, 0.2, 0.01, 0.15]), domains=[[0, 1], [2, 3, 4]], q=[0.05, 0.1])


def test_enumerate_states_sum_to_one(fmodel):
    ex = fmodel.enumerate_states()
    probs = np.array([p for _, p in ex])
    assert probs.sum() == pytest.approx(1.0, abs=1e-9)
    states = np.stack([s for s, _ in ex])
    assert len({tuple(s) for s in states}) == len(states) == 2 ** 5
    # marginals of the enumeration equal marginal()
    np.testing.assert_allclose(probs @ states, fmodel.marginal(), atol=1e-12)


def test_enumerate_states_no_domains_and_too_big():
    fm = FailureModel(np.array([0.3, 0.6]))
    ex = dict((tuple(s.tolist()), p) for s, p in fm.enumerate_states())
    assert ex[(False, False)] == pytest.approx(0.28)
    assert ex[(True, True)] == pytest.approx(0.18)
    assert FailureModel(np.full(13, 0.01)).enumerate_states(max_states=4096) is None


def test_marginal_matches_sampling(fmodel):
    s = fmodel.sample(200_000, np.random.default_rng(8))
    np.testing.assert_allclose(s.mean(0), fmodel.marginal(), atol=0.004)
    expected = 1 - (1 - np.array([0.1, 0.05, 0.2, 0.01, 0.15])) * (1 - np.array([0.05, 0.05, 0.1, 0.1, 0.1]))
    np.testing.assert_allclose(fmodel.marginal(), expected)
    # domain events are common-mode: containers 2,3,4 fail together more often than independently
    joint = (s[:, 2] & s[:, 3]).mean()
    assert joint > s[:, 2].mean() * s[:, 3].mean() * 1.5
    np.testing.assert_array_equal(fmodel.domain_of(), [0, 0, 1, 1, 1])


def test_sampling_matches_enumeration(fmodel):
    ex = fmodel.enumerate_states()
    s = fmodel.sample(100_000, np.random.default_rng(9))
    codes = s @ (1 << np.arange(5))
    freq = np.bincount(codes, minlength=32) / len(s)
    for st_, p in ex:
        assert freq[int(st_ @ (1 << np.arange(5)))] == pytest.approx(p, abs=0.004)


# ----------------------------------------------------------------------------- surrogate math

def _brute_tail(p, k):
    N = len(p)
    tot = 0.0
    for bits in itertools.product([0, 1], repeat=N):
        if sum(bits) >= k:
            tot += np.prod([pi if b else 1 - pi for pi, b in zip(p, bits)])
    return tot


@settings(max_examples=40, deadline=None)
@given(p=st.lists(st.floats(0, 1), min_size=1, max_size=8), k=st.integers(0, 10))
def test_poisson_binomial_tail_vs_bruteforce(p, k):
    got = poisson_binomial_tail(np.array([p]), np.array([k]))[0]
    assert got == pytest.approx(_brute_tail(p, k), abs=1e-10)


def test_poisson_binomial_tail_batched_and_masked():
    rng = np.random.default_rng(10)
    P = rng.random((3, 4, 6))
    K = rng.integers(0, 8, (3, 4))
    out = poisson_binomial_tail(P, K)
    assert out.shape == (3, 4)
    for a, b in itertools.product(range(3), range(4)):
        assert out[a, b] == pytest.approx(_brute_tail(P[a, b], K[a, b]), abs=1e-10)
    mask = np.array([True, True, False, True, False, False])
    m = poisson_binomial_tail(P[0], K[0], mask=np.broadcast_to(mask, P[0].shape))
    for b in range(4):
        assert m[b] == pytest.approx(_brute_tail(P[0, b][mask], K[0, b]), abs=1e-10)


def test_pivotal_weights_vs_finite_difference():
    rng = np.random.default_rng(11)
    S, N = 4, 7
    p = rng.uniform(0.05, 0.95, (S, N))
    need = np.array([1, 3, 5, 7])
    mask = np.ones((S, N), bool)
    mask[2, 6] = False
    W = pivotal_weights(p, need, mask)
    h = 1e-6
    for s in range(S):
        for i in range(N):
            if not mask[s, i]:
                assert W[s, i] == 0.0
                continue
            pp, pm = p.copy(), p.copy()
            pp[s, i] += h
            pm[s, i] -= h
            lp = np.log(poisson_binomial_tail(np.where(mask, pp, 0), need)[s])
            lm = np.log(poisson_binomial_tail(np.where(mask, pm, 0), need)[s])
            assert W[s, i] == pytest.approx((lp - lm) / (2 * h), rel=1e-4, abs=1e-8)


def test_sample_thresholds_inverse_cdf():
    phi = np.array([0.0, 0.5, 0.8, 0.9])
    u = np.array([0.1, 0.5, 0.6, 0.85, 0.95])
    t = sample_thresholds(phi, u)
    np.testing.assert_array_equal(t, [1, 1, 2, 3, np.inf])


def test_default_phi_is_cdf():
    for err in (0.001, 0.01, 0.05):
        phi = default_phi(err)
        assert phi[0] == 0 and (np.diff(phi) >= 0).all() and phi[-1] <= 1
    assert default_phi(0.001)[3] > default_phi(0.05)[3]


# ----------------------------------------------------------------------------- surrogate MC + predictor

def _small_setup(seed=2, S=6, K=8, N=12, m=4):
    cs = make_heterogeneous(m, np.random.default_rng(seed))
    ch = ChannelParams.from_containers(cs)
    stripe_of = np.repeat(np.arange(S), N)
    L = Layout(stripe_of, np.full(S, K), round_robin(stripe_of, m))
    return L, ch


def test_layout_bookkeeping():
    L, ch = _small_setup(S=2)
    assert L.n == 24 and L.m == 4 and L.n_stripes == 2
    assert len(L.pair_frag) == 24
    assert (L.stripe_index >= 0).all() and sorted(L.stripe_index.ravel().tolist()) == list(range(24))
    np.testing.assert_array_equal(L.stored_per_container(), [6, 6, 6, 6])
    np.testing.assert_allclose(L.pi_bar(ch), ch.yield_[L.pair_cont] / 6)


def test_pair_to_cont_min_with_empty_last_container():
    L = Layout(np.zeros(3, int), np.array([1]), np.array([[1, 0, 0], [0, 1, 0], [0, 1, 0]]))
    out = L.pair_to_cont_min(np.array([[5.0, 9.0, 1.0]]))
    np.testing.assert_array_equal(out, [[5.0, 1.0, np.inf]])


def test_pair_to_frag_min_with_unplaced_last_fragment():
    L = Layout(np.zeros(2, int), np.array([1]), np.array([[1, 1], [0, 0]]))
    out = L.pair_to_frag_min(np.array([[5.0, 1.0]]))
    np.testing.assert_array_equal(out, [[1.0, np.inf]])


def test_pair_min_aggregations_regular_layout():
    L, _ = _small_setup(S=2)
    X = np.random.default_rng(0).random((3, len(L.pair_frag)))
    fm = L.pair_to_frag_min(X)
    cm = L.pair_to_cont_min(X)
    for i in range(L.n):
        np.testing.assert_allclose(fm[:, i], X[:, L.pair_frag == i].min(1))
    for j in range(L.m):
        np.testing.assert_allclose(cm[:, j], X[:, L.pair_cont == j].min(1))


def test_sample_T_reproducible():
    L, ch = _small_setup(S=2)
    a = sample_T(L, ch, 50, np.random.default_rng(123))
    b = sample_T(L, ch, 50, np.random.default_rng(123))
    c = sample_T(L, ch, 50, np.random.default_rng(124))
    np.testing.assert_array_equal(a, b)
    assert not np.array_equal(a, c)
    assert a.shape == (50, len(L.pair_frag)) and (a > 0).all()
    T, failed = sample_T(L, ch, 50, np.random.default_rng(123), return_state=True)
    np.testing.assert_array_equal(T, a)
    assert np.isinf(T[failed[:, L.pair_cont]]).all()


def test_sample_T_forced_failure():
    L, ch = _small_setup(S=2)
    failed = np.zeros((20, 4), bool)
    failed[:, 1] = True
    T = sample_T(L, ch, 20, np.random.default_rng(0), failed=failed)
    assert np.isinf(T[:, L.pair_cont == 1]).all()


def test_file_time_matches_definition():
    L, ch = _small_setup(S=2)
    T = sample_T(L, ch, 30, np.random.default_rng(1))
    w = uniform_weights(L)
    tau = fragment_times(L, T, w)
    ft = file_time(L, tau)
    for t in range(30):
        per = [np.sort(tau[t, L.stripe_of == s])[L.k[s] - 1] for s in range(2)]
        assert ft[t] == max(per)


def test_predictor_exact_without_abundance_noise():
    """With no shared/container abundance noise the predictor's only
    approximation vanishes, so it must match Monte Carlo tightly."""
    L, ch = _small_setup(S=2)
    ch = replace(ch, sigma_g=0.0, shape=np.full(4, 1e6))
    pred = Predictor(L, ch)
    w = uniform_weights(L)
    samp = static_budget_samples(L, ch, w, 6000, np.random.default_rng(1))
    for B in np.quantile(samp, [0.2, 0.5, 0.8, 0.95]):
        assert pred.p_file(w * B) == pytest.approx((samp <= B).mean(), abs=0.015)


def test_predictor_agrees_with_monte_carlo():
    L, ch = _small_setup(S=6)
    pred = Predictor(L, ch)
    assert sum(p for _, p in pred.states) == pytest.approx(1.0, abs=1e-6)
    w = uniform_weights(L)
    samp = static_budget_samples(L, ch, w, 6000, np.random.default_rng(1))
    for B in np.quantile(samp, [0.1, 0.3, 0.5, 0.7, 0.9, 0.97]):
        assert pred.p_file(w * B) == pytest.approx((samp <= B).mean(), abs=0.025)
    # monotone in budget, and gradient is non-negative
    vals = [pred.p_file(w * B) for B in (20, 50, 100, 200, 400)]
    assert all(a <= b + 1e-12 for a, b in zip(vals, vals[1:]))
    assert (pred.gradient(w * 100) >= -1e-9).all()
    p, by = pred.p_file(w * 100, by_state=True)
    assert len(by) == len(pred.states) and 0 <= p <= 1


def test_static_budget_reproducible():
    L, ch = _small_setup(S=2)
    w = uniform_weights(L)
    a = static_budget_samples(L, ch, w, 500, np.random.default_rng(5), chunk=128)
    b = static_budget_samples(L, ch, w, 500, np.random.default_rng(5), chunk=128)
    np.testing.assert_array_equal(a, b)
    assert a.shape == (500,)
