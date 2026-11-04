"""Empirical recovery from real sequencing runs.

Each sequencing run of the same synthesized pool is treated as a physically
distinct *container* (aging / PCR history differ). Reads are i.i.d. draws, so
restricting a run's read stream to the strands placed in a container (plus a
proportional share of its junk reads) is an exact emulation of that container
under multinomial sampling.

For every (strand i, container j) we estimate an empirical consensus-success
curve phi_ij(n) from the strand's real reads (random subsets of size n, two-
sided BMA consensus, success = exact match to the design; or Hamming <= h as an
inner-ECC-correctable proxy). A Monte-Carlo trial is a random permutation of
the container's evaluation reads; the recovery time T_ij is the stream
position of the t_ij-th read of strand i with t_ij ~ phi_ij. The resulting T
matrix plugs into the same static/adaptive evaluation used for the surrogate.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..channel.surrogate import NMAX, Layout
from ..io.fastq import read_fastq
from ..retrieval.consensus import consensus
from .mapping import KmerMapper, trim_adapter

N_GRID = (1, 2, 3, 4, 5, 6, 8, 10, 13, 16)


@dataclass
class RunData:
    name: str
    ref_of_read: np.ndarray                  # evaluation split, -1 = junk
    calib_counts: np.ndarray                 # per-strand counts in calibration split
    calib_yield: float
    reads_by_ref: dict = field(default_factory=dict, repr=False)   # eval reads per strand (trimmed)
    pairs_sample: list = field(default_factory=list, repr=False)    # (ref, read) for error profiling
    phi_strand: np.ndarray | None = None     # (strands, len(N_GRID)) success probability
    phi_strand_h2: np.ndarray | None = None
    calib_by_ref: dict = field(default_factory=dict, repr=False)   # calibration reads per strand (trimmed)
    phi_calib: np.ndarray | None = None      # (len(N_GRID),) pooled phi from CALIBRATION reads only (leak-free)


def load_run(name: str, path, mapper: KmerMapper, n_calib: int, n_eval: int, L: int,
             keep_reads_per_ref: int = 24, n_pairs: int = 20000, seed: int = 0,
             keep_calib_per_ref: int = 16) -> RunData:
    rng = np.random.default_rng(seed)
    S = len(mapper.refs)
    calib = np.zeros(S, dtype=np.int64)
    calib_mapped = 0
    refs_eval = np.full(n_eval, -1, dtype=np.int32)
    by_ref: dict[int, list[str]] = {}
    calib_by_ref: dict[int, list[str]] = {}
    pairs = []
    for t, rec in enumerate(read_fastq(path, limit=n_calib + n_eval)):
        seq = trim_adapter(rec.seq)
        m = mapper.map(seq)
        if t < n_calib:
            if m.ref >= 0:
                calib[m.ref] += 1
                calib_mapped += 1
                cl = calib_by_ref.setdefault(m.ref, [])
                if len(cl) < keep_calib_per_ref:
                    cl.append(seq[: L + 8])
                if len(pairs) < n_pairs and rng.random() < 0.05:
                    pairs.append((mapper.refs[m.ref], seq))
            continue
        k = t - n_calib
        refs_eval[k] = m.ref
        if m.ref >= 0:
            lst = by_ref.setdefault(m.ref, [])
            if len(lst) < keep_reads_per_ref:
                lst.append(seq[: L + 8])
    return RunData(name, refs_eval, calib, calib_mapped / max(n_calib, 1), by_ref, pairs, calib_by_ref=calib_by_ref)


def _hamming_ok(a: str, b: str, h: int) -> bool:
    return len(a) == len(b) and sum(x != y for x, y in zip(a, b)) <= h


def strand_phi(reads: list[str], ref: str, rng, draws: int = 2, h: int = 2, two_sided: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Empirical P(consensus of n random reads == ref) on N_GRID (exact and Hamming<=h).
    Forward-only consensus by default: raw Illumina reads carry an uncertain 3' tail
    beyond the strand (read-through), which breaks the reversed BMA pass."""
    L = len(ref)
    ex = np.full(len(N_GRID), np.nan)
    hh = np.full(len(N_GRID), np.nan)
    if not reads:
        return np.zeros(len(N_GRID)), np.zeros(len(N_GRID))
    streak = 0
    for gi, n in enumerate(N_GRID):
        if n > len(reads):
            break
        s_ex = s_h = 0
        for _ in range(draws if n < len(reads) else 1):
            sub = [reads[i] for i in rng.choice(len(reads), n, replace=False)]
            est = consensus(sub, L, two_sided=two_sided) if n > 1 else sub[0][:L]
            s_ex += est == ref
            s_h += _hamming_ok(est, ref, h)
        d = draws if n < len(reads) else 1
        ex[gi], hh[gi] = s_ex / d, s_h / d
        streak = streak + 1 if ex[gi] == 1.0 else 0
        if streak >= 2:              # stable success: assume success for larger n
            ex[gi:], hh[gi:] = 1.0, 1.0
            break
    # forward-fill untested larger n with last value (conservative: no improvement assumed)
    for arr in (ex, hh):
        last = 0.0
        for gi in range(len(arr)):
            if np.isnan(arr[gi]):
                arr[gi] = last
            last = max(last, arr[gi])
            arr[gi] = last
    return ex, hh


def estimate_phi(run: RunData, refs: list[str], seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    S = len(refs)
    run.phi_strand = np.zeros((S, len(N_GRID)))
    run.phi_strand_h2 = np.zeros((S, len(N_GRID)))
    for i in range(S):
        r = run.reads_by_ref.get(i, [])
        run.phi_strand[i], run.phi_strand_h2[i] = strand_phi(r, refs[i], rng)
    # leak-free pooled curve for controller beliefs: calibration-split reads only, separate RNG stream
    if run.calib_by_ref:
        rng_c = np.random.default_rng([seed, 1])
        cal = [strand_phi(run.calib_by_ref.get(i, []), refs[i], rng_c)[0] for i in range(S)]
        run.phi_calib = np.mean(cal, axis=0)


def grid_to_full(phi_grid: np.ndarray, n_max: int = NMAX) -> np.ndarray:
    """Interpolate N_GRID values to n = 0..n_max (step function, monotone)."""
    full = np.zeros(phi_grid.shape[:-1] + (n_max + 1,))
    g = np.array(N_GRID)
    for n in range(1, n_max + 1):
        idx = np.searchsorted(g, n, side="right") - 1
        full[..., n] = phi_grid[..., max(idx, 0)] if idx >= 0 else 0.0
    return full


def empirical_T(layout: Layout, runs: list[RunData], phi_full: list[np.ndarray], trials: int,
                rng: np.random.Generator) -> np.ndarray:
    """(trials, pairs) recovery times in container reads from real read streams."""
    P = len(layout.pair_frag)
    T = np.full((trials, P), np.inf)
    for j, run in enumerate(runs):
        sel = np.flatnonzero(layout.pair_cont == j)
        if sel.size == 0:
            continue
        strands = layout.pair_frag[sel]
        inv = np.zeros(len(run.calib_counts), bool)
        inv[strands] = True
        r = run.ref_of_read
        mapped_mask = (r >= 0) & inv[np.maximum(r, 0)]
        junk = np.flatnonzero(r < 0)
        share = inv.mean()
        for t in range(trials):
            keep_junk = junk[rng.random(junk.size) < share]
            stream = np.concatenate([np.flatnonzero(mapped_mask), keep_junk])
            perm = rng.permutation(stream.size)
            pos = np.empty(stream.size, dtype=np.int64)
            pos[perm] = np.arange(1, stream.size + 1)
            ids = r[stream]
            order = np.lexsort((pos, ids))           # group by strand, then by position
            ids_o, pos_o = ids[order], pos[order]
            starts = np.searchsorted(ids_o, strands, side="left")
            ends = np.searchsorted(ids_o, strands, side="right")
            u = rng.random(sel.size)
            ph = phi_full[j][strands]                 # (sel, NMAX+1)
            need = (u[:, None] > ph).sum(1)            # first n with phi >= u
            ok = (need >= 1) & (need <= NMAX) & (ends - starts >= need)
            T[t, sel[ok]] = pos_o[starts[ok] + need[ok] - 1]
    return T
