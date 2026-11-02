"""Codecs, sequence constraints, whitening and the ingestion layer."""
from __future__ import annotations

import io
import zipfile

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from d3ra.encoding import binary_encoder as be
from d3ra.encoding.constraints import (ConstraintSpec, gc_content, homopolymer_runs, max_homopolymer,
                                       validate, windowed_gc)
from d3ra.encoding.dna_encoder import (CODECS, RotatingTernaryCodec, ScrambledTwoBitCodec, TwoBitCodec,
                                       get_codec)
from d3ra.ingestion import checksum
from d3ra.ingestion.compressor import METHOD_IDS, METHOD_NAMES, compress, decompress
from d3ra.ingestion.file_loader import FileRecord, bundle, sniff_type, unbundle

ALL_CODECS = [TwoBitCodec(), ScrambledTwoBitCodec(), ScrambledTwoBitCodec(seed=12345), RotatingTernaryCodec()]
CODEC_IDS = ["2bit", "scrambled", "scrambled-seed", "rotating"]


# ----------------------------------------------------------------------------- codecs

@pytest.mark.parametrize("codec", ALL_CODECS, ids=CODEC_IDS)
@settings(max_examples=60, deadline=None)
@given(data=st.binary(min_size=0, max_size=64))
def test_codec_roundtrip(codec, data):
    seq = codec.encode(data)
    assert set(seq) <= set("ACGT")
    assert codec.decode(seq, len(data)) == data


@pytest.mark.parametrize("codec", ALL_CODECS, ids=CODEC_IDS)
@settings(max_examples=30, deadline=None)
@given(n=st.integers(0, 40), seed=st.integers(0, 2**32 - 1))
def test_codec_length_preserving(codec, n, seed):
    """All inputs of the same byte length encode to the same strand length."""
    data = np.random.default_rng(seed).integers(0, 256, n, dtype=np.uint8).tobytes()
    assert len(codec.encode(data)) == codec.length_for(n)


def test_two_bit_rate_and_mapping():
    c = TwoBitCodec()
    assert c.encode(b"\x00") == "AAAA"
    assert c.encode(b"\xff") == "TTTT"
    assert c.encode(b"\x1b") == "ACGT"     # 00 01 10 11
    assert c.length_for(10) == 40
    assert c.bits_per_nt == 2.0


def test_scrambled_differs_from_plain_but_same_length():
    data = bytes(32)
    plain, scr = TwoBitCodec().encode(data), ScrambledTwoBitCodec().encode(data)
    assert len(plain) == len(scr)
    assert plain == "A" * 128 and scr != plain
    # whitening breaks the all-zero homopolymer
    assert max_homopolymer(scr) < 20


@settings(max_examples=80, deadline=None)
@given(data=st.binary(min_size=0, max_size=80))
def test_rotating_code_never_has_homopolymers(data):
    seq = RotatingTernaryCodec().encode(data)
    assert max_homopolymer(seq) <= 1
    assert all(a != b for a, b in zip(seq, seq[1:]))


def test_rotating_code_block_sizes():
    c = RotatingTernaryCodec()
    assert 3 ** c.BLOCK_TRITS > 2 ** (8 * c.BLOCK_BYTES)
    assert c.length_for(5) == 26
    assert c.length_for(6) == 52
    assert c.length_for(0) == 0
    # adversarial all-0xFF block (largest value) round trips
    assert c.decode(c.encode(b"\xff" * 15), 15) == b"\xff" * 15


def test_get_codec_registry():
    assert set(CODECS) == {"2bit", "scrambled-2bit", "rotating-ternary"}
    for name in CODECS:
        assert get_codec(name).name == name
    with pytest.raises(KeyError):
        get_codec("nope")


# ----------------------------------------------------------------------------- constraints

def test_gc_content():
    assert gc_content("") == 0.0
    assert gc_content("GGCC") == 1.0
    assert gc_content("ATAT") == 0.0
    assert gc_content("ACGT") == 0.5


def test_max_homopolymer_and_runs():
    assert max_homopolymer("") == 0
    assert max_homopolymer("A") == 1
    assert max_homopolymer("ACGT") == 1
    assert max_homopolymer("AACCCGTTTT") == 4
    assert homopolymer_runs("AACCCCGTTTTT", 4) == [(2, 4, "C"), (7, 5, "T")]
    assert homopolymer_runs("ACGT", 2) == []


@settings(max_examples=50, deadline=None)
@given(seq=st.text(alphabet="ACGT", min_size=1, max_size=60))
def test_max_homopolymer_matches_runs(seq):
    m = max_homopolymer(seq)
    runs = homopolymer_runs(seq, 1)
    assert m == max(length for _, length, _ in runs)
    assert sum(length for _, length, _ in runs) == len(seq)


def test_windowed_gc():
    w = windowed_gc("GGGGAAAA", 4)
    np.testing.assert_allclose(w, [1, 0.75, 0.5, 0.25, 0.0])
    assert windowed_gc("GA", 4).tolist() == [0.5]


def test_validate():
    rep = validate("ACGTACGTGC")
    assert rep.ok and rep.max_run == 1 and abs(rep.gc - 0.6) < 1e-12
    rep = validate("AAAAGC")
    assert not rep.ok_run and not rep.ok
    rep = validate("GGCGCC", ConstraintSpec(max_run=3))
    assert not rep.ok_gc and rep.ok_run
    spec = ConstraintSpec(gc_min=0.0, gc_max=1.0, max_run=10, window=4, window_gc_min=0.25, window_gc_max=0.75)
    assert not validate("GCGCATAT", spec).ok_window
    assert validate("GACTGACT", spec).ok_window
    with pytest.raises(ValueError):
        validate("ACGN")


# ----------------------------------------------------------------------------- whitening

@settings(max_examples=50, deadline=None)
@given(data=st.binary(max_size=200), seed=st.integers(0, 2**31 - 1))
def test_whitening_is_involution(data, seed):
    w = be.whiten(data, seed)
    assert len(w) == len(data)
    assert be.whiten(w, seed) == data
    assert be.unwhiten(w, seed) == data


def test_whitening_changes_data_and_keystream_deterministic():
    data = bytes(64)
    assert be.whiten(data) != data
    np.testing.assert_array_equal(be.keystream(16, 7), np.random.default_rng(7).integers(0, 256, 16, dtype=np.uint8))
    assert be.whiten(data, 1) != be.whiten(data, 2)


@settings(max_examples=30, deadline=None)
@given(data=st.binary(max_size=64))
def test_bits_roundtrip(data):
    bits = be.bytes_to_bits(data)
    assert bits.size == 8 * len(data)
    assert be.bits_to_bytes(bits) == data


# ----------------------------------------------------------------------------- checksum

def test_checksums_known_vectors():
    assert checksum.sha256(b"abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    assert checksum.sha256_bytes(b"abc").hex() == checksum.sha256(b"abc")
    assert checksum.crc32(b"123456789") == 0xCBF43926
    assert checksum.crc32_bytes(b"123456789") == bytes.fromhex("cbf43926")
    assert checksum.verify_sha256(b"abc", checksum.sha256(b"abc"))
    assert not checksum.verify_sha256(b"abd", checksum.sha256(b"abc"))


# ----------------------------------------------------------------------------- compression

SAMPLES = [b"", b"x", bytes(1000), b"hello world " * 50,
           np.random.default_rng(0).integers(0, 256, 500, dtype=np.uint8).tobytes()]


@pytest.mark.parametrize("method", ["none", "zlib", "bz2", "lzma", "auto"])
@pytest.mark.parametrize("data", SAMPLES, ids=["empty", "one", "zeros", "text", "random"])
def test_compress_roundtrip(method, data):
    blob, used = compress(data, method)
    assert used in METHOD_IDS
    if method != "auto":
        assert used == method
    assert decompress(blob, used) == data


@settings(max_examples=25, deadline=None)
@given(data=st.binary(max_size=300))
def test_compress_auto_roundtrip_and_minimal(data):
    blob, used = compress(data, "auto")
    assert decompress(blob, used) == data
    assert len(blob) == min(len(compress(data, m)[0]) for m in METHOD_IDS)
    assert len(blob) <= len(data)       # 'none' is always a candidate


def test_compress_auto_prefers_none_for_random_and_compresses_redundant():
    rnd = SAMPLES[-1]
    assert compress(rnd, "auto")[1] == "none"
    blob, used = compress(bytes(5000), "auto")
    assert used != "none" and len(blob) < 100


def test_compress_unknown_method():
    with pytest.raises(ValueError):
        compress(b"abc", "brotli")
    assert {METHOD_NAMES[i] for i in METHOD_IDS.values()} == set(METHOD_IDS)


# ----------------------------------------------------------------------------- file loader

PDF = b"%PDF-1.7\n%\xe2\xe3\xcf\xd3\n1 0 obj<<>>endobj\n%%EOF"
PNG = b"\x89PNG\r\n\x1a\n" + b"\x00\x00\x00\rIHDR" + bytes(13) + b"\xde\xad\xbe\xef"
JPEG = b"\xff\xd8\xff\xe0\x00\x10JFIF\x00" + bytes(20) + b"\xff\xd9"
GIF = b"GIF89a" + bytes(10)


def _zip_bytes() -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("inner.txt", "zipped content")
    return buf.getvalue()


ZIP = _zip_bytes()
EMPTY_ZIP = b"PK\x05\x06" + bytes(18)


@pytest.mark.parametrize("data,kind", [
    (PDF, "pdf"), (PNG, "png"), (JPEG, "jpeg"), (ZIP, "zip"), (EMPTY_ZIP, "zip"), (GIF, "gif"),
    (b"plain utf-8 text \xc3\xa9", "text"), (b"", "text"), (b"\xff\xfe\x00\x80\x81", "binary"),
])
def test_sniff_type(data, kind):
    assert sniff_type(data) == kind


def test_sniff_type_text_with_multibyte_char_at_4096_boundary():
    data = b"a" * 4095 + "é".encode()   # 2-byte char split by the 4096-byte window
    data.decode("utf-8")                      # it is valid UTF-8 text
    assert sniff_type(data) == "text"


def test_file_record_from_bytes_and_path(tmp_path):
    r = FileRecord.from_bytes("doc.pdf", PDF, version=3)
    assert r.kind == "pdf" and r.size == len(PDF) and r.version == 3
    assert r.sha256 == checksum.sha256(PDF)
    md = r.metadata()
    assert "data" not in md and md["name"] == "doc.pdf"
    p = tmp_path / "img.png"
    p.write_bytes(PNG)
    q = FileRecord.from_path(p)
    assert q.name == "img.png" and q.kind == "png" and q.data == PNG


def test_bundle_unbundle_multi_file():
    files = [FileRecord.from_bytes(n, d, version=i + 1) for i, (n, d) in enumerate(
        [("a.pdf", PDF), ("b.png", PNG), ("c.jpg", JPEG), ("d.zip", ZIP), ("e.txt", b"hello"),
         ("empty.bin", b""), ("rand.bin", SAMPLES[-1])])]
    blob = bundle(files)
    out = unbundle(blob)
    assert len(out) == len(files)
    for a, b in zip(files, out):
        assert (b.name, b.data, b.kind, b.size, b.sha256, b.version) == (a.name, a.data, a.kind, a.size,
                                                                         a.sha256, a.version)
        assert checksum.sha256(b.data) == b.sha256
    assert [f.kind for f in out] == ["pdf", "png", "jpeg", "zip", "text", "text", "binary"]


@settings(max_examples=25, deadline=None)
@given(payloads=st.lists(st.binary(max_size=100), min_size=0, max_size=5))
def test_bundle_roundtrip_arbitrary(payloads):
    files = [FileRecord.from_bytes(f"f{i}", d) for i, d in enumerate(payloads)]
    out = unbundle(bundle(files))
    assert [f.data for f in out] == payloads
    assert [f.name for f in out] == [f.name for f in files]


def test_bundle_survives_compression():
    files = [FileRecord.from_bytes("a.txt", b"abc" * 100), FileRecord.from_bytes("b.png", PNG)]
    blob, m = compress(bundle(files), "auto")
    out = unbundle(decompress(blob, m))
    assert [f.data for f in out] == [f.data for f in files]
