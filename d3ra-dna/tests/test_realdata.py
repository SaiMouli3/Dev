"""Tests for FASTQ I/O, archives, streaming FASTQ retrieval, error transplant
and the real-data estimation utilities (no network needed)."""
import os

import numpy as np
import pytest

from d3ra.archive import address_prefix, load_archive, save_archive
from d3ra.channel.empirical_errors import ErrorTransplant, edit_script
from d3ra.channel.sequencer import noisy_reads, ReadErrorModel
from d3ra.channel.surrogate import Layout
from d3ra.ingestion.file_loader import FileRecord
from d3ra.io.fastq import FastqRecord, read_fasta, read_fastq, write_fasta, write_fastq
from d3ra.realdata import empirical as EM
from d3ra.realdata.channel_estimation import coverage_stats, error_profile, shared_variance_decomposition
from d3ra.realdata.mapping import KmerMapper, error_events, trim_adapter
from d3ra.realdata.references import read_binned
from d3ra.retrieval.fastq_retrieval import FastqRetriever, simulate_fastq
from d3ra.storage.container import make_heterogeneous
from d3ra.storage.library import build_library
from d3ra.storage.placement import round_robin


def rand_seq(rng, n):
    return "".join(rng.choice(list("ACGT"), n))


def test_fastq_roundtrip_gz(tmp_path):
    recs = [FastqRecord(f"r{i}", "ACGT" * 5, "I" * 20) for i in range(10)]
    p = tmp_path / "x.fq.gz"
    assert write_fastq(p, recs) == 10
    back = list(read_fastq(p))
    assert [r.seq for r in back] == [r.seq for r in recs] and back[3].name == "r3"
    assert list(read_fastq(p, limit=4))[-1].name == "r3"


def test_fasta_roundtrip(tmp_path):
    write_fasta(tmp_path / "a.fa", ["x", "y"], ["ACGT", "GGCC"])
    assert read_fasta(tmp_path / "a.fa") == (["x", "y"], ["ACGT", "GGCC"])


def test_binned_reader(tmp_path):
    p = tmp_path / "b.txt"
    p.write_text("AAAA\n*****\nAAAT\nAAAA\n\n\nCCCC\n*****\nCCCG\n\n")
    out = read_binned(p)
    assert out == [("AAAA", ["AAAT", "AAAA"]), ("CCCC", ["CCCG"])]


def test_trim_and_map():
    rng = np.random.default_rng(0)
    refs = [rand_seq(rng, 108) for _ in range(200)]
    mp = KmerMapper(refs)
    read = refs[17][:50] + "T" + refs[17][51:] + "AGATCGGAAGAGCACAC"
    t = trim_adapter(read)
    assert len(t) == 108
    m = mp.map(t)
    assert m.ref == 17 and m.dist == 1
    assert mp.map(rand_seq(rng, 108)).ref == -1


def test_error_events_kinds():
    ref = "ACGTACGTAC"
    assert error_events(ref, ref) == []
    assert [k for _, k in error_events(ref, "ACGAACGTAC")] == ["sub"]
    assert [k for _, k in error_events(ref, "ACGACGTAC")] == ["del"]
    assert [k for _, k in error_events(ref, "ACGTTACGTAC")] == ["ins"]


def test_error_profile_rates():
    rng = np.random.default_rng(1)
    em = ReadErrorModel(0.01, 0.002, 0.005, homopolymer_boost=1.0, quality_sigma=0.0, position_slope=0.0)
    pairs = []
    for _ in range(400):
        ref = rand_seq(rng, 100)
        pairs.append((ref, noisy_reads(ref, 1, em, rng)[0]))
    prof = error_profile(pairs)
    assert abs(prof["rate_sub"] - 0.01) < 0.004
    assert abs(prof["rate_del"] - 0.005) < 0.003
    assert prof["rate_ins"] < 0.006


def test_transplant_reproduces_rates():
    rng = np.random.default_rng(2)
    em = ReadErrorModel(0.02, 0.005, 0.01, homopolymer_boost=1.0, quality_sigma=0.0, position_slope=0.0)
    pairs = [(r, noisy_reads(r, 1, em, rng)[0]) for r in (rand_seq(rng, 110) for _ in range(500))]
    tp = ErrorTransplant(pairs)
    target = rand_seq(rng, 110)
    sim = [(target, noisy_reads(target, 1, tp, rng)[0]) for _ in range(500)]
    p0, p1 = error_profile(pairs), error_profile(sim)
    assert abs(p0["rate_total"] - p1["rate_total"]) < 0.01
    assert edit_script("ACGT", "ACGT") == []


def test_coverage_and_shared_variance():
    rng = np.random.default_rng(3)
    S, J = 4000, 4
    g = rng.normal(0, 0.5, S)                       # shared log effect, var 0.25
    h = rng.normal(0, 0.3, (S, J))                  # specific, var 0.09
    lam = 200 * np.exp(g[:, None] + h)
    C = rng.poisson(lam)
    dec = shared_variance_decomposition(C)
    assert abs(dec["var_shared"] - 0.25) < 0.04
    assert abs(dec["var_container_specific"] - 0.09) < 0.03
    cs = coverage_stats(C[:, 0])
    assert cs["dropout_frac"] < 0.01 and cs["cv"] > 0.4


def test_empirical_T_positions():
    # one container, 2 strands; read stream: strand0, junk, strand1, strand0, strand0
    run = EM.RunData("x", np.array([0, -1, 1, 0, 0], dtype=np.int32), np.array([3, 1]), 0.8)
    L = Layout(np.array([0, 0]), np.array([1]), np.array([[1.0], [1.0]]))
    phi = np.zeros((2, EM.NMAX + 1))
    phi[0, 2:] = 1.0      # strand 0 needs exactly 2 reads
    phi[1, 1:] = 1.0      # strand 1 needs 1 read
    T = EM.empirical_T(L, [run], [phi], 50, np.random.default_rng(0))
    assert np.isfinite(T).all()
    # T for strand 0 is the position of its 2nd read in a permuted 4-5 element stream
    assert (T[:, 0] >= 2).all() and (T[:, 0] <= 5).all() and (T[:, 1] >= 1).all()
    assert (T[:, 0] > 1).all()


def test_archive_prefix_and_fastq_retrieval_end_to_end(tmp_path):
    data = os.urandom(1500) + b"hello" * 200
    lib = build_library([FileRecord.from_bytes("f.bin", data)], 32, 12, 6)
    cs = make_heterogeneous(3, np.random.default_rng(0), 0.3)
    for c in cs:
        c.fail_prob = 0.0
    X = round_robin(lib.stripe_of, 3)
    save_archive(lib, X, cs, tmp_path / "arch")
    man, fmt, placement, containers, addresses = load_archive(tmp_path / "arch")
    for i in (0, 7, len(addresses) - 1):
        assert address_prefix(fmt, addresses[i]) == lib.fragments[i].strand[: fmt.prefix_nt]
    paths = simulate_fastq(tmp_path / "arch", tmp_path / "fq", 3000, seed=1)
    r = FastqRetriever(tmp_path / "arch", paths, "pivotal", batch=100).run()
    assert r["status"] == "exact"
    assert r["result"].files[0].data == data
    assert r["reads_total"] < 9000     # stopped before exhausting the files


def test_fastq_retrieval_missing_container(tmp_path):
    data = os.urandom(1200)
    lib = build_library([FileRecord.from_bytes("g", data)], 32, 12, 8)
    cs = make_heterogeneous(4, np.random.default_rng(1), 0.3)
    X = round_robin(lib.stripe_of, 4)
    save_archive(lib, X, cs, tmp_path / "a")
    paths = simulate_fastq(tmp_path / "a", tmp_path / "fq", 3000, seed=2)
    paths[2] = None          # container 2 lost: 3 of 20 fragments per stripe, parity 8
    r = FastqRetriever(tmp_path / "a", paths, "pivotal", batch=100).run()
    assert r["status"] == "exact" and r["reads_per_container"][2] == 0
