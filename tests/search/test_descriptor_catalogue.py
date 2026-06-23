"""M-ED Phase 2 gate: parse_free_return_arcs over the 12 sourced-descriptor rows.

GOLDEN DISCIPLINE (project memory feedback_golden_tests_sourced_only): this gate
asserts only the *shape/structure* parsed from the catalogue's sourced descriptor
strings (lengths consistent, revs >= 0, branches in the allowed set, seeds > 0).
No V_inf or any value our own code computed is ever asserted here.

The 12 rows are exactly the catalogue entries that carry a non-null
free_return_arcs[] (spec §16.7.7 / docs/spec.md:1018). The 3 gapped Russell rows
and the russell-ocampo-* rows carry null free_return_arcs and are therefore
absent from this set (no explicit skip needed).
"""

from __future__ import annotations

from typing import Any

import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.search.descriptor import parse_free_return_arcs
from tests._catalogue_loader import CATALOGUE_PATH

_ALLOWED_BRANCHES = {"single", "low", "high"}


def _rows_with_arcs() -> list[tuple[str, list[dict[str, Any]]]]:
    rows = yaml.safe_load(CATALOGUE_PATH.read_text())
    return [(r["id"], r["free_return_arcs"]) for r in rows if r.get("free_return_arcs")]


def test_exactly_thirteen_sourced_descriptor_rows() -> None:
    # 13 since #388 (2026-06-23): russell-ocampo-4.3.1-5 ingested its
    # McConaghy-Russell-Longuski 2005 Table 2 per-arc descriptor
    # (4 g(7-1/14,5-1/14 rev,L) f(1:1,...) h(0.5,...)) — was 12.
    assert len(_rows_with_arcs()) == 13


@pytest.mark.parametrize(
    "entry_id,arcs",
    _rows_with_arcs(),
    ids=[rid for rid, _ in _rows_with_arcs()],
)
def test_parse_free_return_arcs_well_formed(entry_id: str, arcs: list[dict[str, Any]]) -> None:
    revs, branches, seeds = parse_free_return_arcs(arcs)
    n = len(arcs)
    assert len(revs) == n
    assert len(branches) == n
    assert len(seeds) == n
    assert all(r >= 0 for r in revs)
    assert all(b in _ALLOWED_BRANCHES for b in branches)
    assert all(s > 0.0 for s in seeds)
