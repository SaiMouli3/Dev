"""Adaptive sequential sequencing allocation (vectorised over Monte-Carlo trials).

At every decision epoch the controller, for each unfinished trial, picks one
container to receive a batch of ``b`` reads, then observes which fragments
became decodable. Controllers:

* ``uniform``        : round-robin over containers (adaptive *stopping* only).
* ``coverage``       : greedy expected number of *useful* new fragments per cost
                       (fragments of stripes that are not yet complete).
* ``pivotal``        : greedy marginal increase of  sum_s log P(stripe s completes
                       within horizon H)  per cost, using pivotal-fragment weights
                       d log P(S>=nu)/d rho_i = P(S_{-i}=nu-1)/P(S>=nu).
* ``thompson``       : ``pivotal`` with container efficiencies sampled from the
                       posterior (exploration under estimation uncertainty).
* ``oracle``         : ``pivotal`` with the *true* channel parameters and the true
                       container-failure state (upper-bound reference).
* ``clairvoyant``    : knows the realised recovery-time matrix (reference only).

Container belief: a discrete posterior over an efficiency multiplier
theta_j in {0 (failed)} U {0.01..5} combined with a binomial likelihood of the
number of recovered fragments U_j after R_j reads. theta absorbs
mis-estimated yield/abundance/error and container failure in one object.
Stopping: as soon as every stripe has K verified fragments (``DecodedStop``),
or the cap is reached, or (optionally) predicted completion probability falls
below ``futility_eps``.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.special import gammaln
from scipy.special import log_ndtr

from ..channel.surrogate import GH_W, GH_X, ChannelParams, Layout, SuccessTable, file_time, sample_T
from .stopping_rule import completion_probability

# Efficiency support must reach deep degradation: with a floor of 0.2, a heavily
# degraded container that returns nothing in its probe is better explained by
# "failed" than "degraded" and is abandoned for good (observed in s7, see report).
THETA = np.concatenate([[0.0], np.exp(np.linspace(np.log(0.01), np.log(5.0), 25))])


# controller name -> base score / theta selection
_BASE_SCORE = {"thompson": "pivotal", "ts": "coverage", "ucb": "coverage", "ltc": "coverage",
               "info": "coverage", "tail": "tail", "tail_lcb": "tail"}
_THETA_MODE = {"thompson": "sample", "ts": "sample", "ucb": "upper", "tail_lcb": "lower"}


@dataclass
class AdaptiveConfig:
    controller: str = "pivotal"
    batch: int = 0                 # reads per action (0 -> auto: 2% of fragments, >=50)
    probe: int = 0                 # initial reads per container (0 -> batch)
    horizon: float = 6.0           # look-ahead in batches per container
    cap: float = np.inf            # cost cap (budget)
    prior_sd: float = 0.5          # log-sd of efficiency prior (estimation uncertainty)
    futility_eps: float = 0.0
    seed: int = 0
    overdispersion: float = 0.0    # beta-binomial ICC of recoveries in the posterior (0 -> binomial)
    explore_rounds: int = 5        # learn-then-commit: round-robin batches per container before committing
    quantile: float = 0.9          # UCB uses this upper posterior quantile, tail/LCB the lower (1-q)
    tail_frac: float = 0.1         # tail-aware: fraction of unfinished stripes treated as bottlenecks
    forecast: str = "normal"       # p_complete: 'normal' (analytic) or 'mc' (posterior simulation)
    forecast_draws: int = 200


@dataclass
class AdaptiveResult:
    cost: np.ndarray
    reads: np.ndarray
    success: np.ndarray
    rounds: np.ndarray
    containers_used: np.ndarray
    reads_per_container: np.ndarray
    extra_after_probe: np.ndarray
    latency: np.ndarray
    meta: dict = field(default_factory=dict)


class _QTable:
    """Q_j(mu) = (1-d_j)(1-d_syn) E_g Psucc_j(g mu): pair decodability vs mean reads."""

    def __init__(self, ch: ChannelParams):
        st = SuccessTable(ch)
        self.grid = st.grid
        g = np.exp(ch.sigma_g * GH_X - ch.sigma_g ** 2 / 2)
        tab = np.zeros_like(st.tab)
        for j in range(ch.m):
            vals = np.stack([np.interp(gg * self.grid, st.grid, st.tab[j]) for gg in g])
            tab[j] = GH_W @ vals
        self.tab = tab * ((1 - ch.dropout) * (1 - ch.d_syn))[:, None]

    def __call__(self, cont: np.ndarray, mu: np.ndarray) -> np.ndarray:
        out = np.empty(mu.shape)
        for j in np.unique(cont):
            sel = cont == j
            out[..., sel] = np.interp(mu[..., sel], self.grid, self.tab[j])
        return out

    def per_container(self, mu: np.ndarray) -> np.ndarray:
        """mu (..., M) -> (..., M)."""
        out = np.empty(mu.shape)
        for j in range(mu.shape[-1]):
            out[..., j] = np.interp(mu[..., j], self.grid, self.tab[j])
        return out


class AdaptiveController:
    def __init__(self, layout: Layout, belief: ChannelParams, cfg: AdaptiveConfig,
                 truth: ChannelParams | None = None):
        self.L, self.ch, self.cfg = layout, belief, cfg
        self.truth = truth if truth is not None else belief
        self.q = _QTable(belief)
        self.pi = layout.pi_bar(belief)                      # per pair expected read share
        n_j = layout.placement.sum(axis=0)
        self.n_j = (layout.placement > 0).sum(axis=0)
        self.pi_cont = np.where(n_j > 0, belief.yield_ / np.maximum(n_j, 1), 0)  # mean pi per container
        self.batch = cfg.batch or max(50, int(0.02 * layout.n))
        self.probe = cfg.probe or self.batch
        # efficiency prior: failure atom + lognormal around 1
        f = belief.failure.marginal()
        lp = -0.5 * (np.log(THETA[1:]) / cfg.prior_sd) ** 2
        lp = np.exp(lp - lp.max())
        lp /= lp.sum()
        self.log_prior = np.log(np.concatenate([f[:, None], (1 - f)[:, None] * lp[None, :]], axis=1) + 1e-300)

    # ------------------------------------------------------------------ belief
    def posterior(self, R: np.ndarray, U: np.ndarray) -> np.ndarray:
        """(trials, M, G) posterior over theta."""
        mu = R[:, :, None] * self.pi_cont[None, :, None] * THETA[None, None, :]
        q = np.empty(mu.shape)
        for gi in range(len(THETA)):
            q[:, :, gi] = self.q.per_container(mu[:, :, gi])
        q = np.clip(q, 1e-12, 1 - 1e-12)
        n = self.n_j[None, :, None]
        u = U[:, :, None]
        r = self.cfg.overdispersion
        if r > 0:
            # beta-binomial likelihood: recoveries within a container are positively correlated
            # (shared abundance / error effects), which the binomial ignores -> over-confidence
            a, b = q * (1 - r) / r, (1 - q) * (1 - r) / r
            ll = gammaln(u + a) + gammaln(n - u + b) - gammaln(n + a + b) - gammaln(a) - gammaln(b) + gammaln(a + b)
        else:
            ll = u * np.log(q) + (n - u) * np.log1p(-q)
        # theta = 0 : no fragments can be recovered
        ll[:, :, 0] = np.where(U > 0, -np.inf, 0.0)
        lpst = ll + self.log_prior[None]
        lpst -= lpst.max(axis=2, keepdims=True)
        p = np.exp(lpst)
        return p / p.sum(axis=2, keepdims=True)

    # ------------------------------------------------------------------ scores
    def _hazard(self, theta: np.ndarray, R: np.ndarray, add: np.ndarray) -> np.ndarray:
        """P(pair recovered within `add` more reads | not yet), (trials, P)."""
        c = self.L.pair_cont
        th = theta[:, c]
        mu0 = th * R[:, c] * self.pi[None, :]
        mu1 = th * (R[:, c] + add[:, c]) * self.pi[None, :]
        q0, q1 = self.q(c, mu0), self.q(c, mu1)
        return np.clip((q1 - q0) / np.maximum(1 - q0, 1e-12), 0, 1)

    def theta_from_post(self, post: np.ndarray, rng) -> np.ndarray:
        """Point value of theta per (trial, container) used for planning."""
        mode = _THETA_MODE.get(self.cfg.controller, "mean")
        if mode == "mean":
            return post @ THETA
        cdf = post.cumsum(axis=2)
        if mode == "sample":
            u = rng.random(cdf.shape[:2] + (1,))
        elif mode == "upper":
            u = np.full(cdf.shape[:2] + (1,), self.cfg.quantile)
        else:  # lower
            u = np.full(cdf.shape[:2] + (1,), 1 - self.cfg.quantile)
        return THETA[(u > cdf).sum(axis=2).clip(max=len(THETA) - 1)]

    def scores(self, R, rec_pair, rec_frag, by_stripe, theta, rng, base: str | None = None) -> np.ndarray:
        L, cfg = self.L, self.cfg
        trials, M = R.shape
        b = np.full((trials, M), float(self.batch))
        cost = self.ch.cost[None, :] * self.batch + self.ch.setup_cost[None, :]
        need = np.maximum(L.k[None, :] - by_stripe, 0)          # (trials, S)
        frag_open = (~rec_frag) & (need[:, L.stripe_of] > 0)     # useful fragments
        pair_open = frag_open[:, L.pair_frag] & (~rec_pair)
        h1 = self._hazard(theta, R, b) * pair_open
        base = base or _BASE_SCORE.get(cfg.controller, cfg.controller)
        if base == "coverage":
            return L.pair_to_cont_sum(h1) / cost
        # pivotal: horizon recovery probabilities per fragment
        H = np.full((trials, M), cfg.horizon * self.batch)
        hH = self._hazard(theta, R, H) * pair_open
        logmiss = L.pair_to_frag_sum(np.log1p(-np.clip(hH, 0, 1 - 1e-12)))
        rho = (1 - np.exp(logmiss)) * frag_open                   # (trials, n)
        idx = L.stripe_index
        mask = idx >= 0
        R_s = np.where(mask[None], rho[:, np.maximum(idx, 0)], 0.0)  # (trials, S, Nmax)
        mu = R_s.sum(-1)
        var = (R_s * (1 - R_s)).sum(-1) + 0.25
        # d log P(S >= need) / d rho_i  ~  pdf / sf  (normal approx, continuity-corrected)
        z = (mu - (need - 0.5)) / np.sqrt(var)
        logpdf = -0.5 * z ** 2 - 0.5 * np.log(2 * np.pi * var)
        wlog = np.exp(logpdf - log_ndtr(z)) * (need > 0)           # (trials, S)
        if base == "tail":
            # risk-aware: concentrate on the bottleneck stripes (lowest standardized slack z),
            # which determine the upper tail of the completion time
            zz = np.where(need > 0, z, np.inf)
            k_b = np.maximum(1, np.ceil(cfg.tail_frac * (need > 0).sum(1))).astype(int)
            thr = np.take_along_axis(np.sort(zz, axis=1), (k_b - 1)[:, None], axis=1)
            Wt = ((zz <= thr) & (need > 0)).astype(float)[:, L.stripe_of]
            contrib = Wt[:, L.pair_frag] * h1
            sc = L.pair_to_cont_sum(contrib)
            # fall back to coverage when no bottleneck pair can be advanced
            sc = np.where(sc.max(1, keepdims=True) > 0, sc, L.pair_to_cont_sum(h1))
            return sc / cost
        W = wlog[:, L.stripe_of]                                   # (trials, n)
        # marginal effect of one batch on pair (i,j): d rho_i = (1-rho_i^{-j}) * h1_ij ~ h1 * (1-rho_i)
        contrib = W[:, L.pair_frag] * h1 * (1 - rho[:, L.pair_frag] + h1)
        return L.pair_to_cont_sum(contrib) / cost

    def p_complete(self, R, rec_pair, rec_frag, by_stripe, theta, U=None, rng=None) -> np.ndarray:
        if self.cfg.forecast == "mc" and U is not None:
            return self.p_complete_mc(R, rec_pair, rec_frag, by_stripe, U, rng or np.random.default_rng(0))
        L = self.L
        rem = np.maximum(self.cfg.cap - (R * self.ch.cost[None, :]).sum(1), 0)
        H = np.repeat((rem / self.ch.cost.sum())[:, None], L.m, 1)  # spread remaining evenly
        need = np.maximum(L.k[None, :] - by_stripe, 0)
        frag_open = ~rec_frag
        hH = self._hazard(theta, R, H) * (frag_open[:, L.pair_frag] & ~rec_pair)
        rho = 1 - np.exp(L.pair_to_frag_sum(np.log1p(-np.clip(hH, 0, 1 - 1e-12))))
        idx = L.stripe_index
        return completion_probability(rho[:, np.maximum(idx, 0)], need, (idx >= 0)[None])

    def p_complete_mc(self, R, rec_pair, rec_frag, by_stripe, U, rng) -> np.ndarray:
        """Posterior-simulation forecast of completing within the remaining cap.
        Draws theta from the posterior, spreads the remaining budget in proportion to the
        reads spent so far (the controller's revealed allocation), draws fragment recoveries
        from the hazards and counts stripes that reach K."""
        L, cfg = self.L, self.cfg
        trials = R.shape[0]
        post = self.posterior(R, U)
        rem = np.maximum(cfg.cap - (R * self.ch.cost[None, :]).sum(1), 0)
        share = R / np.maximum(R.sum(1, keepdims=True), 1)
        H = share * (rem / np.maximum((share * self.ch.cost[None, :]).sum(1), 1e-12))[:, None]
        need = np.maximum(L.k[None, :] - by_stripe, 0)
        open_ = (~rec_frag)[:, L.pair_frag] & ~rec_pair
        cdf = post.cumsum(axis=2)
        out = np.zeros(trials)
        D = cfg.forecast_draws
        for _ in range(D):
            u = rng.random(cdf.shape[:2] + (1,))
            th = THETA[(u > cdf).sum(axis=2).clip(max=len(THETA) - 1)]
            hH = self._hazard(th, R, H) * open_
            miss = L.pair_to_frag_sum(np.log1p(-np.clip(hH, 0, 1 - 1e-12)))
            rho = (1 - np.exp(miss)) * (~rec_frag)
            got = rng.random(rho.shape) < rho
            gain = np.rint(L.frag_to_stripe_sum(got.astype(np.float32))).astype(int)
            out += (gain >= need).all(1)
        return out / D

    # ------------------------------------------------------------------ run
    def run(self, trials: int, rng: np.random.Generator, T: np.ndarray | None = None,
            failed: np.ndarray | None = None, checkpoints=None) -> AdaptiveResult:
        """checkpoints: optional list of spend levels; when an unfinished trial first
        reaches a level, the controller's posterior-predictive probability of completing
        within the remaining cap is recorded (meta['checkpoints']: rows of
        (trial, level, p_complete)). Used to evaluate confidence calibration."""
        L, cfg, ch = self.L, self.cfg, self.ch
        if T is None:
            T, failed = sample_T(L, self.truth, trials, rng, return_state=True)
        trials = T.shape[0]
        M = L.m
        R = np.zeros((trials, M))
        cost = np.zeros(trials)
        rounds = np.zeros(trials, dtype=int)
        active = np.ones(trials, dtype=bool)
        used = L.stored_per_container() > 0
        rr = np.zeros(trials, dtype=int)

        def observe():
            rec_pair = T <= R[:, L.pair_cont]
            rp = rec_pair.astype(np.float32)
            rec_frag = L.pair_to_frag_sum(rp) > 0
            by_stripe = np.rint(L.frag_to_stripe_sum(rec_frag.astype(np.float32))).astype(int)
            U = L.pair_to_cont_sum(rp)
            return rec_pair, rec_frag, by_stripe, U

        # initial probe of every container that stores fragments
        if cfg.controller != "clairvoyant":
            for j in np.flatnonzero(used):
                R[:, j] += self.probe
                cost += ch.cost[j] * self.probe + ch.setup_cost[j]
                rounds += 1
        probe_cost = cost.copy()
        rec_pair, rec_frag, by_stripe, U = observe()
        done = (by_stripe >= L.k[None, :]).all(1)
        active &= ~done
        # simulation shortcut: trials whose realised T makes recovery impossible would
        # sequence until the cap and fail; account for them exactly without simulating.
        tau = L.pair_to_frag_min(np.where(np.isfinite(T), 0.0, np.inf))
        feasible = np.isfinite(file_time(L, tau))
        infeasible = active & ~feasible
        if np.isfinite(cfg.cap):
            cost[infeasible] = cfg.cap
        active &= feasible
        gave_up = np.zeros(trials, dtype=bool)
        explore_end = int(used.sum()) * (1 + cfg.explore_rounds)
        committed = np.zeros((trials, M, len(THETA))) if cfg.controller == "ltc" else None
        cps = sorted(checkpoints) if checkpoints is not None else []
        next_cp = np.zeros(trials, dtype=int)
        cp_rows = []
        it = 0
        while active.any():
            it += 1
            a = np.flatnonzero(active)
            if cfg.controller == "uniform":
                cand = np.flatnonzero(used)
                choice = cand[rr[a] % len(cand)]
                rr[a] += 1
            elif cfg.controller == "clairvoyant":
                choice = self._clairvoyant(T[a], R[a], rec_pair[a], rec_frag[a], by_stripe[a])
            else:
                if cfg.controller == "oracle":
                    theta = np.where(failed[a], 0.0, 1.0)
                else:
                    post = self.posterior(R[a], U[a])
                    if cfg.controller == "ltc":
                        # learn-then-commit: keep updating during exploration, then freeze the belief
                        learning = rounds[a] <= explore_end
                        committed[a[learning]] = post[learning]
                        post = committed[a]
                    theta = self.theta_from_post(post, rng)
                if cfg.controller == "ltc":
                    explore = rounds[a] < explore_end
                    cand = np.flatnonzero(used)
                    rr_choice = cand[rr[a] % len(cand)]
                if cfg.controller == "info":
                    # Bayesian design / pure information: container whose efficiency is most uncertain
                    lt = np.log(np.maximum(THETA, 1e-3))
                    post_i = self.posterior(R[a], U[a])
                    v = (post_i * lt ** 2).sum(2) - ((post_i * lt).sum(2)) ** 2
                    n_open = L.pair_to_cont_sum((~rec_pair[a]).astype(np.float32))
                    sc = v * (n_open > 0) / (ch.cost[None, :] * self.batch + ch.setup_cost[None, :])
                    # once uncertainty is resolved everywhere, behave like coverage-greedy
                    sc_cov = self.scores(R[a], rec_pair[a], rec_frag[a], by_stripe[a], theta, rng, base="coverage")
                    sc = np.where(sc.max(1, keepdims=True) > 1e-3, sc, sc_cov)
                else:
                    sc = self.scores(R[a], rec_pair[a], rec_frag[a], by_stripe[a], theta, rng)
                sc[:, ~used] = -np.inf
                choice = sc.argmax(axis=1)
                if cfg.controller == "ltc":
                    choice = np.where(explore, rr_choice, choice)
                    rr[a[explore]] += 1
                if cfg.futility_eps > 0 and np.isfinite(cfg.cap) and it % 5 == 0:
                    pc = self.p_complete(R[a], rec_pair[a], rec_frag[a], by_stripe[a], theta, U[a], rng)
                    quit_ = pc < cfg.futility_eps
                    gave_up[a[quit_]] = True
                    active[a[quit_]] = False
            if cps and cfg.controller not in ("uniform", "clairvoyant"):
                due = np.flatnonzero((next_cp[a] < len(cps)) &
                                     (cost[a] >= np.array(cps + [np.inf])[np.minimum(next_cp[a], len(cps))]))
                if due.size:
                    th = theta if cfg.controller != "oracle" else np.where(failed[a], 0.0, 1.0)
                    sub = a[due]
                    pc = self.p_complete(R[sub], rec_pair[sub], rec_frag[sub], by_stripe[sub], th[due], U[sub], rng)
                    for tr, v in zip(sub, pc):
                        cp_rows.append((int(tr), cps[next_cp[tr]], float(v)))
                        next_cp[tr] += 1
            R[a, choice] += self.batch
            cost[a] += ch.cost[choice] * self.batch + ch.setup_cost[choice]
            rounds[a] += 1
            rec_pair, rec_frag, by_stripe, U = observe()
            done = (by_stripe >= L.k[None, :]).all(1)
            active &= ~done
            active &= cost < cfg.cap
            if it > 100000:
                break
        success = done & (cost <= cfg.cap + 1e-9)
        latency = (R / self.batch) @ ch.latency
        return AdaptiveResult(cost, R.sum(1), success, rounds, (R > 0).sum(1), R,
                              cost - probe_cost, latency, {"gave_up": gave_up, "batch": self.batch,
                                                           "checkpoints": cp_rows})

    def _clairvoyant(self, T, R, rec_pair, rec_frag, by_stripe) -> np.ndarray:
        """Greedy with knowledge of realised T: pick the container whose next
        batch yields the most useful fragments (weighted by stripe urgency)."""
        L = self.L
        need = np.maximum(L.k[None, :] - by_stripe, 0)
        frag_open = (~rec_frag) & (need[:, L.stripe_of] > 0)
        pair_open = frag_open[:, L.pair_frag] & ~rec_pair
        urg = (need / np.maximum(L.n_per_stripe[None, :] - by_stripe, 1))[:, L.stripe_of]
        hit = pair_open & (T <= R[:, L.pair_cont] + self.batch)
        gain = L.pair_to_cont_sum(hit * urg[:, L.pair_frag])
        cost = self.ch.cost[None, :] * self.batch
        sc = gain / cost
        none = sc.max(1) <= 0
        if none.any():
            # jump toward the nearest useful event
            dt = np.where(pair_open, T - R[:, L.pair_cont], np.inf)
            best = L.pair_to_cont_min(dt)
            sc[none] = -best[none] * self.ch.cost[None, :]
        return sc.argmax(1)
