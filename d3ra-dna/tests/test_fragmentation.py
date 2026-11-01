"""Splitting, physical addressing, GF(256) arithmetic and the outer/inner codes."""
from __future__ import annotations

import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from d3ra.coding.erasure_codes import (LTFountainCode, MDSStripeCode, ReplicationCode, make_code,
                                       robust_soliton)
from d3ra.coding.redundancy import RedundancyPlan
from d3ra.coding.reed_solomon import GF256, CauchyRS, InnerRS
from d3ra.fragmentation.addressing import AddressLayout, address_overhead, overhead_tradeoff
from d3ra.fragmentation.identifiers import fragment_uid, manifest_id, object_id
from d3ra.fragmentation.splitter import Fragment, join_stripes, split_into_stripes


# ----------------------------------------------------------------------------- splitter

@settings(max_examples=60, deadline=None)
@given(blob=st.binary(max_size=500), payload=st.integers(1, 17), k=st.integers(1, 9))
def test_split_join_roundtrip(blob, payload, k):
    stripes = split_into_stripes(blob, payload, k)
    assert len(stripes) == max(1, math.ceil(len(blob) / (payload * k)))
    for s in stripes:
        assert s.shape == (k, payload) and s.dtype == np.uint8
    assert join_stripes(stripes, len(blob)) == blob
    # zero padding of the tail
    flat = b"".join(s.tobytes() for s in stripes)
    assert flat[len(blob):] == bytes(len(flat) - len(blob))


def test_split_stripes_are_writable_copies():
    stripes = split_into_stripes(b"abcdefgh", 2, 2)
    stripes[0][0, 0] = 0      # must not raise (frombuffer views are read-only)
    assert len(stripes) == 2


def test_fragment_properties():
    f = Fragment(obj=3, stripe=1, index=5, n_in_stripe=6, k_in_stripe=4, payload=b"xx")
    assert f.is_parity and f.uid == (3, 1, 5)
    assert not Fragment(0, 0, 3, 6, 4, b"").is_parity


def test_identifiers():
    assert fragment_uid(255, 3, 7) == "ff:3:7"
    assert 0 <= object_id("f" * 64, bits=8) < 256
    assert object_id("00000010" + "0" * 56, bits=8) == 16
    assert manifest_id(b"abc") == "ba7816bf8f01cfea"


# ----------------------------------------------------------------------------- addressing

@settings(max_examples=80, deadline=None)
@given(ob=st.integers(1, 3), sb=st.integers(1, 3), pb=st.integers(1, 2), data=st.data())
def test_address_pack_unpack(ob, sb, pb, data):
    lay = AddressLayout(ob, sb, pb)
    o = data.draw(st.integers(0, 256 ** ob - 1))
    s = data.draw(st.integers(0, 256 ** sb - 1))
    p = data.draw(st.integers(0, 256 ** pb - 1))
    b = lay.pack(o, s, p)
    assert len(b) == lay.size == ob + sb + pb
    assert lay.unpack(b) == (o, s, p)


def test_address_pack_overflow():
    with pytest.raises(OverflowError):
        AddressLayout(1, 1, 1).pack(0, 256, 0)


@pytest.mark.parametrize("count,nbytes", [(1, 1), (2, 1), (255, 1), (256, 1), (257, 2), (65536, 2),
                                          (65537, 3)])
def test_compact_sizing(count, nbytes):
    lay = AddressLayout.compact(1, count, 2)
    assert lay.stripe_bytes == nbytes
    # every index 0..count-1 is representable
    assert lay.unpack(lay.pack(0, count - 1, 1)) == (0, count - 1, 1)


def test_compact_layout_fields():
    lay = AddressLayout.compact(3, 1000, 40)
    assert (lay.obj_bytes, lay.stripe_bytes, lay.pos_bytes) == (1, 2, 1)
    assert lay.size == 4
    assert AddressLayout().size == 4


def test_address_overhead_and_tradeoff():
    d = address_overhead(AddressLayout(1, 2, 1), 32, 6)
    assert d["strand_bytes"] == 42
    assert d["address_fraction"] == pytest.approx(4 / 42)
    assert d["payload_fraction"] == pytest.approx(32 / 42)
    rows = overhead_tradeoff(100_000)
    assert [r["payload_bytes"] for r in rows] == [16, 24, 32, 48, 64]
    for r in rows:
        assert r["n_fragments"] == math.ceil(100_000 / r["payload_bytes"])
        assert r["address_fraction"] + r["payload_fraction"] < 1
    fr = [r["address_fraction"] for r in rows]
    assert all(a > b for a, b in zip(fr, fr[1:]))      # longer payload -> less address overhead


# ----------------------------------------------------------------------------- GF(256)

def test_gf256_tables():
    a = np.arange(1, 256)
    assert (GF256.MUL[a, GF256.INV[a]] == 1).all()
    assert (GF256.MUL[0] == 0).all() and (GF256.MUL[:, 0] == 0).all()
    assert (GF256.MUL[1, :] == np.arange(256)).all()
    assert (GF256.MUL == GF256.MUL.T).all()
    assert len(set(GF256.exp[:255].tolist())) == 255      # 2 is a generator
    # distributivity over XOR on a sample
    rng = np.random.default_rng(0)
    x, y, z = rng.integers(0, 256, (3, 200))
    assert (GF256.MUL[x, y ^ z] == (GF256.MUL[x, y] ^ GF256.MUL[x, z])).all()


@pytest.mark.parametrize("n", [1, 2, 5, 12])
def test_gf256_inverse_matrix(n):
    rng = np.random.default_rng(n)
    eye = np.eye(n, dtype=np.uint8)
    for _ in range(10):
        M = rng.integers(0, 256, (n, n), dtype=np.uint8)
        try:
            Minv = GF256.inv_matrix(M)
        except np.linalg.LinAlgError:
            continue
        np.testing.assert_array_equal(GF256.matmul(M, Minv), eye)
        np.testing.assert_array_equal(GF256.matmul(Minv, M), eye)


def test_gf256_singular_matrix():
    M = np.array([[1, 2], [1, 2]], dtype=np.uint8)
    with pytest.raises(np.linalg.LinAlgError):
        GF256.inv_matrix(M)


# ----------------------------------------------------------------------------- Cauchy RS

@settings(max_examples=60, deadline=None)
@given(k=st.integers(1, 12), parity=st.integers(0, 8), width=st.integers(1, 8), data=st.data())
def test_cauchy_rs_any_k_of_n(k, parity, width, data):
    n = k + parity
    code = CauchyRS(k, n)
    rng = np.random.default_rng(data.draw(st.integers(0, 2**32 - 1)))
    msg = rng.integers(0, 256, (k, width), dtype=np.uint8)
    cw = code.encode(msg)
    assert cw.shape == (n, width)
    np.testing.assert_array_equal(cw[:k], msg)                       # systematic
    keep = data.draw(st.lists(st.integers(0, n - 1), min_size=k, max_size=k, unique=True))
    out = code.decode({i: cw[i] for i in keep})
    np.testing.assert_array_equal(out, msg)


def test_cauchy_rs_every_erasure_pattern_small():
    import itertools
    k, n = 4, 7
    code = CauchyRS(k, n)
    msg = np.random.default_rng(1).integers(0, 256, (k, 5), dtype=np.uint8)
    cw = code.encode(msg)
    for keep in itertools.combinations(range(n), k):
        np.testing.assert_array_equal(code.decode({i: cw[i] for i in keep}), msg)


def test_cauchy_rs_limits():
    with pytest.raises(ValueError):
        CauchyRS(0, 4)
    with pytest.raises(ValueError):
        CauchyRS(5, 4)
    with pytest.raises(ValueError):
        CauchyRS(10, 257)
    code = CauchyRS(250, 256)
    msg = np.random.default_rng(0).integers(0, 256, (250, 2), dtype=np.uint8)
    cw = code.encode(msg)
    keep = list(range(6, 256))
    np.testing.assert_array_equal(code.decode({i: cw[i] for i in keep}), msg)
    with pytest.raises(ValueError):
        CauchyRS(4, 6).decode({0: cw[0, :1], 1: cw[1, :1]})
    assert CauchyRS(4, 6).rate == pytest.approx(4 / 6)


# ----------------------------------------------------------------------------- inner RS

@pytest.mark.parametrize("nsym,correctable", [(4, 2), (6, 3), (8, 4)])
def test_inner_rs_corrects_byte_errors(nsym, correctable):
    irs = InnerRS(nsym)
    rng = np.random.default_rng(nsym)
    for _ in range(20):
        msg = rng.integers(0, 256, 40, dtype=np.uint8).tobytes()
        word = bytearray(irs.encode(msg))
        assert len(word) == len(msg) + irs.overhead()
        for p in rng.choice(len(word), correctable, replace=False):
            word[p] ^= int(rng.integers(1, 256))
        assert irs.decode(bytes(word)) == msg


def test_inner_rs_rejects_heavy_corruption():
    irs = InnerRS(4)
    rng = np.random.default_rng(42)
    msg = bytes(range(32))
    for _ in range(200):
        word = bytearray(irs.encode(msg))
        for p in rng.choice(len(word), 12, replace=False):
            word[p] ^= int(rng.integers(1, 256))
        assert irs.decode(bytes(word)) is None


def test_inner_rs_without_parity_uses_crc():
    irs = InnerRS(nsym=0, crc_bytes=2)
    assert irs.overhead() == 2
    w = irs.encode(b"payload")
    assert irs.decode(w) == b"payload"
    bad = bytearray(w)
    bad[0] ^= 1
    assert irs.decode(bytes(bad)) is None


# ----------------------------------------------------------------------------- other codes

def test_mds_stripe_code_wrapper():
    code = make_code("rs", 3, 5)
    assert isinstance(code, MDSStripeCode)
    msg = np.arange(12, dtype=np.uint8).reshape(3, 4)
    cw = code.encode(msg)
    assert code.can_decode([0, 3, 4]) and not code.can_decode([0, 0, 1])
    np.testing.assert_array_equal(code.decode({0: cw[0], 3: cw[3], 4: cw[4]}), msg)


def test_robust_soliton_is_distribution():
    for k in (1, 2, 8, 50):
        mu = robust_soliton(k)
        assert mu.shape == (k,)
        assert (mu >= 0).all() and mu.sum() == pytest.approx(1.0)


def test_lt_fountain_decode():
    k, n = 8, 24
    code = make_code("lt", k, n)
    assert isinstance(code, LTFountainCode)
    rng = np.random.default_rng(3)
    msg = rng.integers(0, 256, (k, 6), dtype=np.uint8)
    cw = code.encode(msg)
    for i in range(n):        # each symbol is the XOR of its neighbours
        np.testing.assert_array_equal(cw[i], np.bitwise_xor.reduce(msg[code.masks[i]], axis=0))
    assert not code.can_decode(range(k - 1))
    decoded = 0
    for _ in range(30):
        ids = rng.permutation(n)[: rng.integers(k, n + 1)]
        if code.can_decode(ids):
            np.testing.assert_array_equal(code.decode({int(i): cw[i] for i in ids}), msg)
            decoded += 1
        else:
            with pytest.raises(ValueError):
                code.decode({int(i): cw[i] for i in ids})
    assert decoded > 0
    assert code.can_decode(range(n))
    # deterministic construction (decoder must rebuild the same graph)
    np.testing.assert_array_equal(LTFountainCode(k, n).masks, code.masks)


def test_replication_code():
    code = make_code("replication", 3, 9)
    assert isinstance(code, ReplicationCode) and code.r == 3 and code.n == 9
    msg = np.arange(6, dtype=np.uint8).reshape(3, 2)
    cw = code.encode(msg)
    assert cw.shape == (9, 2)
    assert code.can_decode([0, 4, 8]) and not code.can_decode([0, 3, 6])
    np.testing.assert_array_equal(code.decode({8: cw[8], 0: cw[0], 4: cw[4]}), msg)
    with pytest.raises(ValueError):
        code.decode({0: cw[0], 3: cw[3]})
    with pytest.raises(ValueError):
        make_code("replication", 3, 8)
    with pytest.raises(ValueError):
        make_code("turbo", 3, 8)


def test_redundancy_plan():
    p = RedundancyPlan.uniform(3, k=4, parity=2, copies=2)
    assert p.n_stripes == 3 and p.n == [6, 6, 6]
    assert p.unique_fragments() == 18
    assert p.total_fragments() == 36
    assert p.storage_overhead() == pytest.approx(36 / 12)
    assert p.code_rate() == pytest.approx(12 / 18)
    q = RedundancyPlan([4, 4], [5, 7])
    assert [r.tolist() for r in q.replication] == [[1] * 5, [1] * 7]
    assert q.total_fragments() == 12 and q.storage_overhead() == pytest.approx(12 / 8)
