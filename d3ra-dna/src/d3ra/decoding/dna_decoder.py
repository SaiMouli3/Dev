"""Strand codec: (address, payload) <-> DNA, with inner RS + CRC and an
address check byte used by inventory-aware clustering."""
from __future__ import annotations

from dataclasses import dataclass

from ..coding.reed_solomon import InnerRS
from ..encoding.dna_encoder import Codec, get_codec
from ..fragmentation.addressing import AddressLayout
from ..ingestion.checksum import crc32


@dataclass
class StrandFormat:
    layout: AddressLayout
    payload_bytes: int
    inner: InnerRS
    codec: Codec
    addr_check: int = 1   # bytes of address check (raises codebook distance)

    @classmethod
    def build(cls, layout: AddressLayout, payload_bytes: int, inner_nsym: int = 4,
              codec: str = "scrambled-2bit", crc_bytes: int = 2, addr_check: int = 1) -> "StrandFormat":
        return cls(layout, payload_bytes, InnerRS(inner_nsym, crc_bytes), get_codec(codec), addr_check)

    @property
    def header_bytes(self) -> int:
        return self.layout.size + self.addr_check

    @property
    def message_bytes(self) -> int:
        return self.header_bytes + self.payload_bytes

    @property
    def strand_bytes(self) -> int:
        return self.message_bytes + self.inner.overhead()

    @property
    def strand_nt(self) -> int:
        return self.codec.length_for(self.strand_bytes)

    @property
    def prefix_nt(self) -> int:
        """nt that encode the address header (used for clustering)."""
        return self.codec.length_for(self.header_bytes) if self.codec.name != "rotating-ternary" \
            else int(self.header_bytes * 8 / self.codec.bits_per_nt)

    def _header(self, obj: int, stripe: int, pos: int) -> bytes:
        a = self.layout.pack(obj, stripe, pos)
        return a + crc32(a).to_bytes(4, "big")[: self.addr_check]

    def encode(self, obj: int, stripe: int, pos: int, payload: bytes) -> str:
        assert len(payload) == self.payload_bytes
        return self.codec.encode(self.inner.encode(self._header(obj, stripe, pos) + payload))

    def decode(self, seq: str) -> tuple[tuple[int, int, int], bytes] | None:
        word = self.codec.decode(seq, self.strand_bytes)
        msg = self.inner.decode(word)
        if msg is None:
            return None
        hdr, payload = msg[: self.header_bytes], msg[self.header_bytes:]
        addr = self.layout.unpack(hdr[: self.layout.size])
        if hdr != self._header(*addr):
            return None
        return addr, payload
