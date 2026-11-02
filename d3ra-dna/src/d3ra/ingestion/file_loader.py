"""Digital file ingestion: arbitrary binary input, type sniffing, metadata,
multi-file bundling with versioned manifest.

Every file must be recoverable byte-for-byte; the bundle format is a simple,
self-describing length-prefixed container whose manifest carries SHA-256 of
each member so verification is independent of the storage path.
"""
from __future__ import annotations

import json
import struct
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

from .checksum import sha256

MAGIC = {
    b"%PDF": "pdf",
    b"\x89PNG\r\n\x1a\n": "png",
    b"\xff\xd8\xff": "jpeg",
    b"PK\x03\x04": "zip",
    b"PK\x05\x06": "zip",
    b"GIF8": "gif",
}


def sniff_type(data: bytes, name: str = "") -> str:
    for sig, kind in MAGIC.items():
        if data.startswith(sig):
            return kind
    import codecs
    try:
        codecs.getincrementaldecoder("utf-8")().decode(data[:4096], final=len(data) <= 4096)
        return "text"
    except UnicodeDecodeError:
        return "binary"


@dataclass
class FileRecord:
    name: str
    data: bytes = field(repr=False)
    kind: str = "binary"
    size: int = 0
    sha256: str = ""
    version: int = 1
    created: float = 0.0

    @classmethod
    def from_bytes(cls, name: str, data: bytes, version: int = 1) -> "FileRecord":
        return cls(name=name, data=data, kind=sniff_type(data, name), size=len(data),
                   sha256=sha256(data), version=version, created=time.time())

    @classmethod
    def from_path(cls, path: str | Path, version: int = 1) -> "FileRecord":
        p = Path(path)
        return cls.from_bytes(p.name, p.read_bytes(), version)

    def metadata(self) -> dict:
        d = asdict(self)
        d.pop("data")
        return d


def bundle(files: list[FileRecord]) -> bytes:
    """Serialise several files into one byte string (manifest + payloads)."""
    manifest = json.dumps([f.metadata() for f in files], separators=(",", ":")).encode()
    out = [struct.pack(">I", len(manifest)), manifest]
    out += [f.data for f in files]
    return b"".join(out)


def unbundle(blob: bytes) -> list[FileRecord]:
    (mlen,) = struct.unpack(">I", blob[:4])
    manifest = json.loads(blob[4:4 + mlen])
    pos = 4 + mlen
    files = []
    for m in manifest:
        data = blob[pos:pos + m["size"]]
        pos += m["size"]
        files.append(FileRecord(data=data, **m))
    return files
