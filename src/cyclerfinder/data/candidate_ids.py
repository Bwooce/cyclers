"""Stable ids for the candidate records of the two-working-body gauntlet JSONs (#1071).

The ``key`` of a record is the cycle template of the best member of a merged group, and two
groups can share a template, so ``key`` is not an identifier
(``data/998_pluto_smallmoons/pk_gauntlet.json``: 62 candidates, 56 distinct keys). New writes
carry a ``cid``; committed files do not and are addressed by their index in ``candidates``, or by
:func:`record_id` (the key plus the ideal worst ratio, which is unique in every committed file).
This module imports nothing from the literature gate.
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterable, Sequence
from typing import Any

Record = dict[str, Any]


def candidate_cid(cell: str, index: int, x_days: Sequence[float]) -> str:
    """``<cell>-<4-digit zero-based index>-<8 hex of sha1 of the dates rounded to 6 decimals>``."""
    digest = hashlib.sha1(repr(tuple(round(float(v), 6) for v in x_days)).encode()).hexdigest()
    return f"{cell}-{index:04d}-{digest[:8]}"


def record_ratio(rec: Record) -> float | None:
    """The record's ideal worst demanded/available ratio (``worst_ratio``, else the largest
    ``ratio`` of its ``flyby_table``); ``None`` when it has none."""
    if rec.get("worst_ratio") is not None:
        return float(rec["worst_ratio"])
    ratios = [t["ratio"] for t in rec.get("flyby_table", []) if t.get("ratio") is not None]
    return float(max(ratios)) if ratios else None


def record_id(rec: Record) -> str:
    """``cid`` when the record has one, else ``"<key>@<ratio to 9 decimals>"``."""
    if rec.get("cid"):
        return str(rec["cid"])
    ratio = record_ratio(rec)
    return f"{rec['key']}@{'none' if ratio is None else format(ratio, '.9f')}"


def index_records(cands: Iterable[Record]) -> dict[str, Record]:
    """Dict of records by :func:`record_id`; raises if two records share an id."""
    out: dict[str, Record] = {}
    for rec in cands:
        rid = record_id(rec)
        if rid in out:
            raise ValueError(f"duplicate record id {rid!r}")
        out[rid] = rec
    return out


def index_by_key(cands: Iterable[Record]) -> dict[str, Record]:
    """Dict of records by ``key``; raises on an ambiguous key instead of keeping the last one."""
    out: dict[str, Record] = {}
    for rec in cands:
        if rec["key"] in out:
            raise ValueError(f"ambiguous key {rec['key']!r}: use record_id / index_records")
        out[rec["key"]] = rec
    return out


def same_record(a: Record, b: Record) -> bool:
    """Same candidate: equal ``cid`` when both have one, else equal key and ratio (1e-9)."""
    if a.get("cid") and b.get("cid"):
        return bool(a["cid"] == b["cid"])
    ra, rb = record_ratio(a), record_ratio(b)
    if ra is None or rb is None:
        return bool(a["key"] == b["key"] and ra is rb)
    return bool(a["key"] == b["key"] and abs(ra - rb) < 1e-9)
