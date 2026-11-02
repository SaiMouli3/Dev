"""Manifest (logical metadata) and SQLite persistence of experiments."""
from __future__ import annotations

import json
import sqlite3
from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass
class Manifest:
    encoding_version: int
    codec: str
    compression: str
    blob_length: int
    blob_sha256: str
    payload_bytes: int
    inner_nsym: int
    crc_bytes: int
    addr_check: int
    layout: dict
    code: str
    stripes: list[dict]           # [{"k":..,"n":..}]
    files: list[dict] = field(default_factory=list)
    obj: int = 0

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=1)

    @classmethod
    def from_json(cls, s: str) -> "Manifest":
        return cls(**json.loads(s))


class ResultStore:
    """Tiny SQLite store for experiment rows (JSON payload per row)."""

    def __init__(self, path: str | Path):
        self.conn = sqlite3.connect(str(path))
        self.conn.execute("CREATE TABLE IF NOT EXISTS results (experiment TEXT, key TEXT, payload TEXT)")

    def add(self, experiment: str, key: str, payload: dict) -> None:
        self.conn.execute("INSERT INTO results VALUES (?,?,?)", (experiment, key, json.dumps(payload)))
        self.conn.commit()

    def get(self, experiment: str) -> list[dict]:
        cur = self.conn.execute("SELECT payload FROM results WHERE experiment=?", (experiment,))
        return [json.loads(r[0]) for r in cur.fetchall()]
