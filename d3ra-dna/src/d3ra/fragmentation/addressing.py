"""Physical address layout inside each strand.

The *physical* address is the minimum information a read must carry to be
placed: (object id, stripe, position-in-stripe). Everything else (total
fragment count, code parameters, encoding version, file hash) is *logical*
metadata held in the manifest, so the DNA overhead of addressing can be
measured separately from the payload.
"""
from __future__ import annotations

import math
from dataclasses import dataclass


def _bytes_for(count: int) -> int:
    return max(1, math.ceil(math.log2(max(count, 2)) / 8))


@dataclass(frozen=True)
class AddressLayout:
    obj_bytes: int = 1
    stripe_bytes: int = 2
    pos_bytes: int = 1

    @classmethod
    def compact(cls, n_objects: int, n_stripes: int, max_n: int) -> "AddressLayout":
        return cls(_bytes_for(n_objects), _bytes_for(n_stripes), _bytes_for(max_n))

    @property
    def size(self) -> int:
        return self.obj_bytes + self.stripe_bytes + self.pos_bytes

    def pack(self, obj: int, stripe: int, pos: int) -> bytes:
        return (obj.to_bytes(self.obj_bytes, "big") + stripe.to_bytes(self.stripe_bytes, "big")
                + pos.to_bytes(self.pos_bytes, "big"))

    def unpack(self, b: bytes) -> tuple[int, int, int]:
        o = int.from_bytes(b[: self.obj_bytes], "big")
        s = int.from_bytes(b[self.obj_bytes: self.obj_bytes + self.stripe_bytes], "big")
        p = int.from_bytes(b[self.obj_bytes + self.stripe_bytes: self.size], "big")
        return o, s, p


def address_overhead(layout: AddressLayout, payload_bytes: int, inner_overhead: int) -> dict:
    total = layout.size + payload_bytes + inner_overhead
    return {
        "address_bytes": layout.size,
        "payload_bytes": payload_bytes,
        "inner_ecc_bytes": inner_overhead,
        "strand_bytes": total,
        "address_fraction": layout.size / total,
        "payload_fraction": payload_bytes / total,
    }


def overhead_tradeoff(total_bytes: int, payload_options=(16, 24, 32, 48, 64), inner=6, max_n=255):
    """Address overhead vs payload size (longer strands: less overhead, but
    synthesis/sequencing error per strand grows)."""
    rows = []
    for p in payload_options:
        n_frag = math.ceil(total_bytes / p)
        layout = AddressLayout.compact(1, math.ceil(n_frag / max_n) + 1, max_n)
        rows.append({"payload_bytes": p, "n_fragments": n_frag, **address_overhead(layout, p, inner)})
    return rows
