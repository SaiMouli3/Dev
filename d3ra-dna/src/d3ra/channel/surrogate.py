"""Fragment-level surrogate channel ("recovery-time matrix" model).

Sequencing container j is modelled as a Poisson process in that container's
*read clock*: each read independently originates from stored fragment i with
probability pi_ij = y_j c_ij g_i h_ij / sum_i' c_i'j g_i' h_i'j (y_j = usable
yield). A fragment cluster in container j becomes decodable once it holds
t_ij reads, where t_ij is drawn from the consensus-success curve
phi_j(n) = P(t_ij <= n) calibrated by the sequence-level simulator
(BMA consensus + inner RS + CRC). Hence the read count at which fragment i
becomes recoverable from container j is

    T_ij = Gamma(t_ij, 1) / pi_ij      (inf if dropped out or container failed)

and *any* allocation policy can be evaluated exactly on a sampled T matrix:
fragment i is recovered after container reads R iff min_j (T_ij - R_j) <= 0.
For static allocations the minimum total budget per trial is closed-form,
which lets us estimate reliability targets such as 0.9999 with 10^4-10^5
Monte-Carlo trials on a laptop CPU.

The same parameters drive an analytical predictor (negative-binomial read
counts, Gauss-Hermite quadrature over the shared sequence-intrinsic bias,
exact Poisson-binomial stripe tails, failure-state averaging).
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy import sparse, stats
from scipy.special import gammaln

from .failures import FailureModel

NMAX = 48


def default_phi(err: float, n_max: int = NMAX, nsym: int = 4, L: int = 140) -> np.ndarray:
    """Parametric fallback consensus-success curve phi(n) (used when no
    sequence-level calibration is available). Single-read success requires
    no indel and few substitutions; with more reads consensus success
    approaches 1 at a rate governed by the per-base error rate."""
    n = np.arange(n_max + 1, dtype=float)
    indel_free = np.exp(-0.55 * err * L)
    s1 = indel_free * stats.poisson.cdf(nsym // 2 * 4, 0.45 * err * L)
    # effective "reads needed" scale grows ~ log(L)/(-log(4 err (1-err)))
    eps = min(max(err, 1e-5), 0.24)
    n0 = max(1.0, np.log(L) / -np.log(4 * eps * (1 - eps)))
    phi = 1 - (1 - s1) * np.exp(-np.maximum(n - 1, 0) / n0) ** 1.0
    phi = np.where(n >= 2, 1 - (1 - s1) * np.exp(-(np.maximum(n - 1, 0) / n0) ** 1.3), phi)
    phi[0] = 0.0
    phi[1] = s1
    return np.clip(np.maximum.accumulate(phi), 0, 1)


@dataclass
class ChannelParams:
    """Per-container parameters for the surrogate (arrays of length M)."""

    yield_: np.ndarray
    shape: np.ndarray          # gamma shape of container-specific abundance variation
    dropout: np.ndarray        # d_j
    cost: np.ndarray           # c_j per read
    phi: np.ndarray            # (M, NMAX+1) consensus success curves
    failure: FailureModel
    sigma_g: float = 0.3       # log-sd of shared sequence-intrinsic bias
    d_syn: float = 0.002       # shared synthesis dropout
    setup_cost: np.ndarray | None = None
    latency: np.ndarray | None = None

    @property
    def m(self) -> int:
        return len(self.yield_)

    def __post_init__(self):
        for name in ("yield_", "shape", "dropout", "cost"):
            setattr(self, name, np.asarray(getattr(self, name), dtype=float))
        if self.setup_cost is None:
            self.setup_cost = np.zeros(self.m)
        if self.latency is None:
            self.latency = np.ones(self.m)

    @classmethod
    def from_containers(cls, containers, failure: FailureModel | None = None, phi=None,
                        sigma_g: float = 0.3, d_syn: float = 0.002, L: int = 140, nsym: int = 4):
        if failure is None:
            failure = FailureModel(np.array([c.fail_prob for c in containers]))
        if phi is None:
            phi = np.stack([default_phi(c.error_rate, L=L, nsym=nsym) for c in containers])
        return cls(np.array([c.yield_ for c in containers]), np.array([c.abundance_shape for c in containers]),
                   np.array([c.dropout for c in containers]), np.array([c.cost_per_read for c in containers]),
                   phi, failure, sigma_g, d_syn, np.array([c.setup_cost for c in containers]),
                   np.array([c.latency_per_run for c in containers]))

    def perturbed(self, rng: np.random.Generator, rel_sd: float = 0.5) -> "ChannelParams":
        """Imperfect estimate of this channel (for robustness experiments)."""
        mul = lambda x: x * np.exp(rng.normal(0, rel_sd, np.shape(x)))
        phi = self.phi.copy()
        # mis-estimated error rate: stretch/compress the consensus curve
        for j in range(self.m):
            s = np.exp(rng.normal(0, rel_sd))
            n = np.arange(phi.shape[1])
            phi[j] = np.interp(n / s, n, self.phi[j])
        return replace(self, yield_=np.clip(mul(self.yield_), 0.01, 1.0), shape=mul(self.shape),
                       dropout=np.clip(mul(self.dropout), 1e-5, 0.9), phi=phi,
                       failure=FailureModel(np.clip(mul(self.failure.f), 1e-6, 0.9),
                                            self.failure.domains, list(self.failure.q)))


@dataclass
class Layout:
    """Placement of a striped code: (pair -> fragment, container, copies)."""

    stripe_of: np.ndarray      # (n,) stripe of each fragment
    k: np.ndarray              # (S,) K per stripe
    placement: np.ndarray      # (n, M) copy counts

    def __post_init__(self):
        self.placement = np.asarray(self.placement, dtype=float)
        self.pair_frag, self.pair_cont = np.nonzero(self.placement)
        self.pair_copies = self.placement[self.pair_frag, self.pair_cont]
        P, n, m = len(self.pair_frag), self.placement.shape[0], self.placement.shape[1]
        one = np.ones(P)
        self.A_pf = sparse.csr_matrix((one, (np.arange(P), self.pair_frag)), shape=(P, n))
        self.A_pc = sparse.csr_matrix((one, (np.arange(P), self.pair_cont)), shape=(P, m))
        self.frag_start = np.searchsorted(self.pair_frag, np.arange(n))
        self.frag_has = np.bincount(self.pair_frag, minlength=n) > 0
        self.cont_order = np.argsort(self.pair_cont, kind="stable")
        self.cont_start = np.searchsorted(self.pair_cont[self.cont_order], np.arange(m))
        self.cont_has = np.bincount(self.pair_cont, minlength=m) > 0
        self.n_stripes = len(self.k)
        self.n_per_stripe = np.bincount(self.stripe_of, minlength=self.n_stripes)
        self.A_fs = sparse.csr_matrix((np.ones(n), (np.arange(n), self.stripe_of)), shape=(n, self.n_stripes))
        self.nmax = int(self.n_per_stripe.max())
        # padded (S, Nmax) index of fragments per stripe
        order = np.argsort(self.stripe_of, kind="stable")
        self.stripe_index = np.full((self.n_stripes, self.nmax), -1)
        pos = np.zeros(self.n_stripes, dtype=int)
        for i in order:
            s = self.stripe_of[i]
            self.stripe_index[s, pos[s]] = i
            pos[s] += 1

    @property
    def n(self) -> int:
        return self.placement.shape[0]

    # fast aggregations (avoid ufunc.at)
    def pair_to_frag_sum(self, X):
        return np.asarray(X @ self.A_pf)

    def pair_to_cont_sum(self, X):
        return np.asarray(X @ self.A_pc)

    def frag_to_stripe_sum(self, X):
        return np.asarray(X @ self.A_fs)

    @staticmethod
    def _group_min(X, start, has):
        """Segment minimum over sorted pair columns; empty groups -> inf.
        An inf sentinel column is appended so empty groups never truncate a
        neighbouring group (reduceat semantics)."""
        Xp = np.concatenate([X, np.full((X.shape[0], 1), np.inf)], axis=1)
        st = np.where(has, start, X.shape[1])
        # reduceat needs non-decreasing indices for correct segments: use only non-empty groups
        out = np.full((X.shape[0], len(start)), np.inf)
        idx = np.flatnonzero(has)
        if idx.size:
            out[:, idx] = np.minimum.reduceat(Xp, np.append(st[idx], X.shape[1]), axis=1)[:, : idx.size]
        return out

    def pair_to_frag_min(self, X):
        return self._group_min(X, self.frag_start, self.frag_has)

    def pair_to_cont_min(self, X):
        return self._group_min(X[:, self.cont_order], self.cont_start, self.cont_has)

    @property
    def m(self) -> int:
        return self.placement.shape[1]

    def pi_bar(self, ch: ChannelParams) -> np.ndarray:
        """Expected per-read probability of each pair (y_j c_ij / sum_i c_ij)."""
        tot = self.placement.sum(axis=0)
        tot[tot == 0] = 1
        return ch.yield_[self.pair_cont] * self.pair_copies / tot[self.pair_cont]

    def stored_per_container(self) -> np.ndarray:
        return (self.placement > 0).sum(axis=0)


# --------------------------------------------------------------------------- sampling

def sample_thresholds(phi_rows: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Inverse-CDF sampling of reads-needed t from phi rows; inf if never."""
    t = (u[..., None] > phi_rows).sum(axis=-1).astype(float)
    t[t >= phi_rows.shape[-1]] = np.inf
    return t


def sample_T(layout: Layout, ch: ChannelParams, trials: int, rng: np.random.Generator,
             failed: np.ndarray | None = None, return_state: bool = False):
    """Recovery-time matrix: (trials, n_pairs) container-read counts at which
    each (fragment, container) pair becomes decodable."""
    P = len(layout.pair_frag)
    n, m = layout.n, layout.m
    f, c = layout.pair_frag, layout.pair_cont
    sg = ch.sigma_g
    g = rng.lognormal(-sg ** 2 / 2, sg, (trials, n)) if sg > 0 else np.ones((trials, n))
    k = ch.shape[c]
    h = rng.gamma(k[None, :], 1.0 / k[None, :], (trials, P))
    w = layout.pair_copies[None, :] * g[:, f] * h
    tot = layout.pair_to_cont_sum(w)
    pi = ch.yield_[c][None, :] * w / tot[:, c]
    phi_rows = ch.phi[c]                                    # (P, NMAX+1)
    t = sample_thresholds(phi_rows[None, :, :], rng.random((trials, P)))
    finite = np.isfinite(t)
    T = np.full((trials, P), np.inf)
    T[finite] = rng.gamma(t[finite], 1.0) / pi[finite]
    present = rng.random((trials, P)) >= ch.dropout[c][None, :]
    syn_ok = rng.random((trials, n)) >= ch.d_syn
    present &= syn_ok[:, f]
    if failed is None:
        failed = ch.failure.sample(trials, rng)
    alive = ~failed[:, c]
    T[~(present & alive)] = np.inf
    if return_state:
        return T, failed
    return T


def fragment_times(layout: Layout, T: np.ndarray, w: np.ndarray) -> np.ndarray:
    """tau_i = min_j T_ij / w_j : total-read time (under allocation fractions w)
    at which fragment i is recovered."""
    trials = T.shape[0]
    wj = np.asarray(w, float)[layout.pair_cont]
    with np.errstate(divide="ignore", invalid="ignore"):
        scaled = np.where(wj[None, :] > 0, T / wj[None, :], np.inf)
    return layout.pair_to_frag_min(scaled)


def file_time(layout: Layout, tau: np.ndarray) -> np.ndarray:
    """Minimum total budget at which every stripe has K recovered fragments."""
    idx = layout.stripe_index
    pad = np.where(idx[None] >= 0, tau[:, np.maximum(idx, 0)], np.inf)   # (trials, S, Nmax)
    kth = np.take_along_axis(np.sort(pad, axis=2), (layout.k - 1)[None, :, None], axis=2)[..., 0]
    return kth.max(axis=1)


def static_budget_samples(layout: Layout, ch: ChannelParams, w: np.ndarray, trials: int,
                          rng: np.random.Generator, chunk: int = 2000) -> np.ndarray:
    """Per-trial minimal total *reads* for exact recovery under fixed fractions w."""
    out = []
    for s in range(0, trials, chunk):
        T = sample_T(layout, ch, min(chunk, trials - s), rng)
        out.append(file_time(layout, fragment_times(layout, T, w)))
    return np.concatenate(out)


# --------------------------------------------------------------------------- analytics

GH_X, GH_W = np.polynomial.hermite_e.hermegauss(9)
GH_W = GH_W / GH_W.sum()


class SuccessTable:
    """Psucc_j(mu) = sum_n NB(n; k_j, mu) phi_j(n), tabulated on a log grid."""

    def __init__(self, ch: ChannelParams, grid=np.concatenate([[0.0], np.logspace(-3, 4, 400)])):
        self.grid = grid
        n = np.arange(ch.phi.shape[1])
        tab = np.zeros((ch.m, grid.size))
        dmu = np.zeros((ch.m, grid.size))
        for j in range(ch.m):
            k = ch.shape[j]
            mu = grid[:, None]
            # NB pmf with mean mu, shape k
            mu = np.maximum(mu, 1e-12)
            logp, log1mp = np.log(k / (k + mu)), np.log(mu / (k + mu))
            logpmf = (gammaln(n + k) - gammaln(k) - gammaln(n + 1))[None, :] + k * logp + n[None, :] * log1mp
            pmf = np.exp(logpmf)
            tail = np.clip(1 - pmf.sum(axis=1), 0, 1)   # mass beyond NMAX
            tab[j] = pmf @ ch.phi[j] + tail * ch.phi[j, -1]
        self.tab = np.maximum.accumulate(tab, axis=1)

    def __call__(self, j: np.ndarray, mu: np.ndarray) -> np.ndarray:
        out = np.empty(np.shape(mu))
        mu = np.asarray(mu)
        for jj in np.unique(j):
            sel = (j == jj)
            out[..., sel] = np.interp(mu[..., sel], self.grid, self.tab[jj])
        return out


def poisson_binomial_tail(p: np.ndarray, k: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
    """P(sum Bernoulli(p) >= k) row-wise. p: (..., N); k broadcastable (...)."""
    p = np.asarray(p, float)
    if mask is not None:
        p = np.where(mask, p, 0.0)
    N = p.shape[-1]
    dp = np.zeros(p.shape[:-1] + (N + 1,))
    dp[..., 0] = 1.0
    for i in range(N):
        pi = p[..., i:i + 1]
        dp[..., 1:] = dp[..., 1:] * (1 - pi) + dp[..., :-1] * pi
        dp[..., :1] = dp[..., :1] * (1 - pi)
    kk = np.broadcast_to(np.asarray(k), p.shape[:-1])
    cdf_rev = np.cumsum(dp[..., ::-1], axis=-1)[..., ::-1]   # P(S >= x)
    kk = np.clip(kk, 0, N + 1)
    ext = np.concatenate([cdf_rev, np.zeros(p.shape[:-1] + (1,))], axis=-1)
    return np.take_along_axis(ext, kk[..., None].astype(int), axis=-1)[..., 0]


def pivotal_weights(p: np.ndarray, need: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """d log P(S>=need) / d p_i = P(S_{-i} = need-1) / P(S >= need), exact via
    prefix/suffix Poisson-binomial convolutions. p, mask: (S, N)."""
    p = np.where(mask, p, 0.0)
    S, N = p.shape
    pre = np.zeros((S, N + 1, N + 1))
    pre[:, 0, 0] = 1
    for i in range(N):
        pre[:, i + 1, 1:] = pre[:, i, 1:] * (1 - p[:, i:i + 1]) + pre[:, i, :-1] * p[:, i:i + 1]
        pre[:, i + 1, 0] = pre[:, i, 0] * (1 - p[:, i])
    suf = np.zeros((S, N + 1, N + 1))
    suf[:, N, 0] = 1
    for i in range(N - 1, -1, -1):
        suf[:, i, 1:] = suf[:, i + 1, 1:] * (1 - p[:, i:i + 1]) + suf[:, i + 1, :-1] * p[:, i:i + 1]
        suf[:, i, 0] = suf[:, i + 1, 0] * (1 - p[:, i])
    tail = poisson_binomial_tail(p, need)
    w = np.zeros((S, N))
    for i in range(N):
        a, b = pre[:, i, :], suf[:, i + 1, :]          # pmfs excluding i
        # P(S_{-i} = need-1) = sum_t a[t] b[need-1-t]
        conv = np.zeros(S)
        for s in range(S):
            nn = int(need[s]) - 1
            if 0 <= nn <= N:
                conv[s] = np.dot(a[s, : nn + 1], b[s, nn::-1][: nn + 1])
        w[:, i] = conv
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(mask, w / np.maximum(tail[:, None], 1e-300), 0.0)


@dataclass
class Predictor:
    """Analytical recovery predictor for static read allocations."""

    layout: Layout
    ch: ChannelParams
    n_states: int = 64
    seed: int = 1234
    states: list = field(default_factory=list)

    def __post_init__(self):
        self.table = SuccessTable(self.ch)
        self.pi = self.layout.pi_bar(self.ch)
        ex = self.ch.failure.enumerate_states(max_states=2 ** 12)
        if ex is not None:
            self.states = [(s, p) for s, p in ex if p > 1e-9]
        else:
            rng = np.random.default_rng(self.seed)
            st = self.ch.failure.sample(self.n_states, rng)
            st[0] = False  # include the nominal state (variance reduction)
            p0 = np.prod(1 - self.ch.failure.marginal())
            rest = (1 - p0) / max(1, (st[1:].any(axis=1)).sum())
            self.states = [(st[0], p0)] + [(s, rest) for s in st[1:] if s.any()]
            # unbiased-ish mixture: nominal exactly weighted, failures via MC

    def pair_success(self, R: np.ndarray) -> np.ndarray:
        """(G, P) prob a pair is decodable given g quadrature node (excl. dropout)."""
        g = np.exp(self.ch.sigma_g * GH_X - self.ch.sigma_g ** 2 / 2)
        mu = g[:, None] * (np.asarray(R, float)[self.layout.pair_cont] * self.pi)[None, :]
        return self.table(self.layout.pair_cont, mu)

    def fragment_probs(self, R: np.ndarray, failed: np.ndarray | None = None) -> np.ndarray:
        L = self.layout
        q = self.pair_success(R) * (1 - self.ch.dropout[L.pair_cont])[None, :]
        if failed is not None:
            q = q * (~failed[L.pair_cont])[None, :]
        logmiss = L.pair_to_frag_sum(np.log1p(-np.clip(q, 0, 1 - 1e-15)))
        p = GH_W @ (1 - np.exp(logmiss))
        return p * (1 - self.ch.d_syn)

    def stripe_probs(self, R, failed=None) -> np.ndarray:
        p = self.fragment_probs(R, failed)
        idx = self.layout.stripe_index
        P = np.where(idx >= 0, p[np.maximum(idx, 0)], 0.0)
        return poisson_binomial_tail(P, self.layout.k)

    def p_file(self, R, by_state: bool = False):
        vals = [(np.prod(self.stripe_probs(R, s if s.any() else None)), w) for s, w in self.states]
        tot = sum(w for _, w in vals)
        p = sum(v * w for v, w in vals) / tot
        return (p, vals) if by_state else p

    def gradient(self, R, rel: float = 0.02) -> np.ndarray:
        R = np.asarray(R, float)
        base = self.p_file(R)
        g = np.zeros_like(R)
        for j in range(len(R)):
            d = max(1.0, rel * R[j])
            Rj = R.copy()
            Rj[j] += d
            g[j] = (self.p_file(Rj) - base) / d
        return g
