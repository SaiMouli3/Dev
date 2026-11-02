"""End-to-end: library construction, strand decoding, outer decoding, verification."""
from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from d3ra.coding.redundancy import RedundancyPlan
from d3ra.decoding.dna_decoder import StrandFormat
from d3ra.decoding.reconstruction import reconstruct
from d3ra.decoding.verification import RecoveryStatus, classify
from d3ra.encoding.constraints import max_homopolymer
from d3ra.fragmentation.addressing import AddressLayout
from d3ra.ingestion.checksum import sha256
from d3ra.ingestion.file_loader import FileRecord
from d3ra.pipeline import store_and_retrieve
from d3ra.storage.container import make_heterogeneous
from d3ra.storage.library import build_library
from d3ra.storage.metadata import Manifest, ResultStore

RNG_DATA = np.random.default_rng(5).integers(0, 256, 600, dtype=np.uint8).tobytes()
PNG = b"\x89PNG\r\n\x1a\n" + bytes(40)


def _files():
    return [FileRecord.from_bytes("a.bin", RNG_DATA), FileRecord.from_bytes("b.png", PNG),
            FileRecord.from_bytes("c.txt", b"hello DNA " * 20)]


def _rows(lib, keep=None):
    """{(stripe, index): payload} decoded from the strands themselves."""
    rows = {}
    for f in lib.fragments:
        if keep is not None and (f.stripe, f.index) not in keep:
            continue
        dec = lib.fmt.decode(f.strand)
        assert dec is not None
        (obj, s, i), payload = dec
        assert (obj, s, i) == f.uid
        rows[(s, i)] = payload
    return rows


@pytest.fixture(scope="module")
def lib():
    return build_library(_files(), payload_bytes=16, k=6, parity=3, inner_nsym=4, compression="none")


def _assert_exact(res, files):
    assert res.status == RecoveryStatus.EXACT
    assert [f.data for f in res.files] == [f.data for f in files]
    assert [f.name for f in res.files] == [f.name for f in files]
    assert all(v == RecoveryStatus.EXACT for v in res.file_status.values())
    assert res.stripes_ok == res.stripes_total


# ----------------------------------------------------------------------------- library

def test_library_structure(lib):
    man = lib.manifest
    S = len(man.stripes)
    assert S == -(-man.blob_length // (16 * 6)) and S >= 3
    assert len(lib.fragments) == 9 * S
    assert all((s["k"], s["n"]) == (6, 9) and "crc32" in s for s in man.stripes)
    assert len({len(s) for s in lib.strands}) == 1 and len(lib.strands[0]) == lib.fmt.strand_nt
    assert len(set(lib.strands)) == len(lib.strands)
    st_ = lib.stats()
    assert st_["fragments"] == 9 * S and st_["code_rate"] == pytest.approx(6 / 9)
    assert lib.total_nt() == 9 * S * lib.fmt.strand_nt
    assert lib.index_of()[(0, 1, 2)] == 9 + 2
    np.testing.assert_array_equal(lib.k_per_stripe(), np.full(S, 6))
    assert [f.is_parity for f in lib.fragments[:9]] == [False] * 6 + [True] * 3


def test_strand_format_roundtrip_and_rejects_corruption():
    fmt = StrandFormat.build(AddressLayout(1, 2, 1), 20, 4, "scrambled-2bit")
    payload = bytes(range(20))
    s = fmt.encode(3, 513, 7, payload)
    assert len(s) == fmt.strand_nt == 4 * fmt.strand_bytes
    assert fmt.decode(s) == ((3, 513, 7), payload)
    # one substituted base (<= 1 byte error) is corrected by the inner code
    s1 = s[:10] + ("A" if s[10] != "A" else "C") + s[11:]
    assert fmt.decode(s1) == ((3, 513, 7), payload)
    # garbage is rejected
    rng = np.random.default_rng(0)
    junk = "".join(rng.choice(list("ACGT"), len(s)))
    assert fmt.decode(junk) is None
    rot = StrandFormat.build(AddressLayout(1, 1, 1), 10, 4, "rotating-ternary")
    rs = rot.encode(0, 1, 2, bytes(10))
    assert max_homopolymer(rs) == 1 and rot.decode(rs) == ((0, 1, 2), bytes(10))
    assert 0 < rot.prefix_nt < rot.strand_nt and fmt.prefix_nt == 4 * fmt.header_bytes


def test_manifest_json_roundtrip(lib, tmp_path):
    m2 = Manifest.from_json(lib.manifest.to_json())
    assert m2 == lib.manifest
    _assert_exact(reconstruct(m2, _rows(lib)), _files())
    store = ResultStore(tmp_path / "r.sqlite")
    store.add("exp", "k1", {"a": 1})
    store.add("exp", "k2", {"a": 2})
    store.add("other", "k", {"a": 3})
    assert store.get("exp") == [{"a": 1}, {"a": 2}]


# ----------------------------------------------------------------------------- reconstruction

def test_reconstruct_all_fragments(lib):
    _assert_exact(reconstruct(lib.manifest, _rows(lib)), _files())


@settings(max_examples=25, deadline=None)
@given(seed=st.integers(0, 2**32 - 1))
def test_reconstruct_from_any_k_of_n(lib, seed):
    rng = np.random.default_rng(seed)
    keep = set()
    for s, meta in enumerate(lib.manifest.stripes):
        for i in rng.choice(meta["n"], meta["k"], replace=False):
            keep.add((s, int(i)))
    _assert_exact(reconstruct(lib.manifest, _rows(lib, keep)), _files())


def test_reconstruct_parity_only(lib):
    """Every stripe recovered from parity-heavy subsets (systematic rows mostly lost)."""
    keep = {(s, i) for s in range(len(lib.manifest.stripes)) for i in range(3, 9)}
    _assert_exact(reconstruct(lib.manifest, _rows(lib, keep)), _files())


def _corrupt(rows, key):
    p = bytearray(rows[key])
    p[0] ^= 0x55
    rows[key] = bytes(p)


def test_integrity_mismatch_when_payload_corrupted_legacy_manifest(lib):
    import copy
    man = copy.deepcopy(lib.manifest)
    for st in man.stripes:
        st.pop("crc32", None)          # no per-stripe CRC: only the file hash can catch it
    rows = _rows(lib, {(s, i) for s in range(len(man.stripes)) for i in range(6)})
    _corrupt(rows, (1, 2))
    res = reconstruct(man, rows)
    assert res.status == RecoveryStatus.INTEGRITY_MISMATCH
    assert res.files == [] and res.stripes_ok == res.stripes_total


def test_stripe_crc_detects_and_localises_corruption(lib):
    rows = _rows(lib, {(s, i) for s in range(len(lib.manifest.stripes)) for i in range(6)})   # exactly K
    _corrupt(rows, (1, 2))
    res = reconstruct(lib.manifest, rows)
    assert res.status == RecoveryStatus.PARTIAL and res.unverified_stripes == [1]


def test_stripe_crc_repairs_with_extra_fragments(lib):
    rows = _rows(lib, {(s, i) for s in range(len(lib.manifest.stripes)) for i in range(8)})   # K + 2
    _corrupt(rows, (1, 2))
    _corrupt(rows, (2, 7))           # corrupt but unused (decoder takes the first K): no repair needed
    _corrupt(rows, (3, 0))
    res = reconstruct(lib.manifest, rows)
    assert res.status == RecoveryStatus.EXACT
    assert dict(res.repaired) == {1: (2,), 3: (0,)}


def test_partial_and_failed(lib):
    S = len(lib.manifest.stripes)
    rows = _rows(lib)
    short = {key: v for key, v in rows.items() if not (key[0] == 1 and key[1] >= 5)}   # stripe 1: 5 < K
    res = reconstruct(lib.manifest, short)
    assert res.status == RecoveryStatus.PARTIAL
    assert (res.stripes_ok, res.stripes_total) == (S - 1, S) and res.files == []
    few = {key: v for key, v in rows.items() if key[1] < 5}
    res = reconstruct(lib.manifest, few)
    assert res.status == RecoveryStatus.FAILED and res.stripes_ok == 0
    assert reconstruct(lib.manifest, {}).status == RecoveryStatus.FAILED


def test_classify():
    assert classify(b"x", sha256(b"x"), 1, 1) == RecoveryStatus.EXACT
    assert classify(b"y", sha256(b"x"), 1, 1) == RecoveryStatus.INTEGRITY_MISMATCH
    assert classify(None, "", 0, 3) == RecoveryStatus.FAILED
    assert classify(None, "", 2, 3) == RecoveryStatus.PARTIAL


@pytest.mark.parametrize("codec", ["2bit", "scrambled-2bit", "rotating-ternary"])
@pytest.mark.parametrize("compression", ["none", "zlib", "auto"])
def test_library_codecs_and_compression(codec, compression):
    files = _files()
    lb = build_library(files, payload_bytes=20, k=5, parity=2, codec=codec, compression=compression)
    assert lb.manifest.codec == codec
    keep = {(s, i) for s in range(len(lb.manifest.stripes)) for i in (0, 2, 4, 5, 6)}
    _assert_exact(reconstruct(lb.manifest, _rows(lb, keep)), files)


def test_library_with_replication_and_lt_codes():
    files = _files()
    S = build_library(files, payload_bytes=16, k=4, parity=4, compression="none").manifest.stripes
    nS = len(S)
    rep = build_library(files, payload_bytes=16, k=4, parity=4, compression="none", code="replication")
    assert rep.manifest.code == "replication" and len(rep.fragments) == 8 * nS
    keep = {(s, i) for s in range(nS) for i in (4, 1, 6, 7)}          # one copy of each data row
    _assert_exact(reconstruct(rep.manifest, _rows(rep, keep)), files)
    bad = {(s, i) for s in range(nS) for i in (0, 4, 1, 5)}           # rows 2, 3 missing
    assert reconstruct(rep.manifest, _rows(rep, bad)).status == RecoveryStatus.FAILED
    plan = RedundancyPlan.uniform(nS, 4, 12, code="lt")
    lt = build_library(files, payload_bytes=16, k=4, compression="none", plan=plan)
    _assert_exact(reconstruct(lt.manifest, _rows(lt)), files)


def test_library_with_heterogeneous_plan():
    files = _files()
    base = build_library(files, payload_bytes=16, k=6, parity=2, compression="none")
    nS = len(base.manifest.stripes)
    plan = RedundancyPlan([6] * nS, [7 + (s % 3) for s in range(nS)])
    lb = build_library(files, payload_bytes=16, k=6, compression="none", plan=plan)
    assert [m["n"] for m in lb.manifest.stripes] == plan.n
    rng = np.random.default_rng(0)
    keep = {(s, int(i)) for s, n in enumerate(plan.n) for i in rng.choice(n, 6, replace=False)}
    _assert_exact(reconstruct(lb.manifest, _rows(lb, keep)), files)
    with pytest.raises(AssertionError):
        build_library(files, payload_bytes=16, k=6, plan=RedundancyPlan([6], [8]))


# ----------------------------------------------------------------------------- pipeline

@pytest.mark.parametrize("placement,controller", [("tiered", "pivotal"), ("domain", "uniform"),
                                                  ("round_robin", "coverage")])
def test_store_and_retrieve_end_to_end(placement, controller):
    files = [FileRecord.from_bytes("a.bin", RNG_DATA)]
    cs = make_heterogeneous(5, np.random.default_rng(1), n_domains=0)
    out = store_and_retrieve(files, cs, placement=placement, controller=controller, k=8, parity=4,
                             calibrate=False, failed=np.zeros(5, bool), seed=0)
    assert out["status"] == "exact"
    res = out["result"]
    assert res.files[0].data == RNG_DATA and res.files[0].name == "a.bin"
    assert out["reads"] > 0 and np.isfinite(out["cost"]) and out["failed"] == [False] * 5
    assert out["library"]["fragments"] == 12 * out["library"]["stripes"]
    assert sum(len(c.inventory) for c in cs) == out["library"]["fragments"]


def test_store_and_retrieve_fixed_reads_and_failed_container():
    files = [FileRecord.from_bytes("a.bin", RNG_DATA), FileRecord.from_bytes("b.png", PNG)]
    cs = make_heterogeneous(4, np.random.default_rng(2))
    failed = np.array([False, True, False, False])
    out = store_and_retrieve(files, cs, placement="round_robin", k=8, parity=4, calibrate=False,
                             failed=failed, fixed_reads_per_fragment=20, seed=1)
    # round robin puts 3 of 12 fragments per stripe in each container: one failure is survivable
    assert out["status"] == "exact"
    assert [f.data for f in out["result"].files] == [RNG_DATA, PNG]
    assert out["failed"] == failed.tolist()
