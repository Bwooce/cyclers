"""#1024 regression: resonant_circle and half_rev_vectors are BIT-IDENTICAL on the Keplerian cells.

The fixture was recorded from the pre-edit code (scripts/make_1024_regression_fixture.py) before the
barycentric option touched src/cyclerfinder/search/two_working_body.py (pre-registration amendment
A1, docs/notes/2026-10-11-1024-barycentric-binary-model.md sec. 3a). The ev, gc, vm and #998 pn
cells are all Keplerian, so their #973 / #998 numbers cannot move without this test failing.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

from cyclerfinder.search.two_working_body import (
    HalfRevLeg,
    ResonantLeg,
    half_rev_vectors,
    resonant_circle,
)

REPO = Path(__file__).resolve().parents[2]
FIXTURE = REPO / "tests" / "search" / "fixtures" / "1024_keplerian_regression.json"


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_keplerian_cells_bit_identical() -> None:
    maker = _load(
        "make_1024_regression_fixture", REPO / "scripts" / "make_1024_regression_fixture.py"
    )
    systems = {cell: s for cell, (s, _, _) in maker.cells().items()}
    rec = json.loads(FIXTURE.read_text())
    n = 0
    for cell, rows in rec.items():
        system = systems[cell]
        for r in rows:
            if r["kind"] == "R":
                c = resonant_circle(system, ResonantLeg(r["body"], *r["nm"]), r["t"], r["v"])
                got = None if c is None else [c[0], c[1], *[list(map(float, u)) for u in c[2:]]]
            else:
                h, k, p = r["hkp"]
                vs = half_rev_vectors(system, HalfRevLeg(r["body"], h, k, p, 0), r["t"], r["v"])
                got = [list(map(float, u)) for u in vs]
            assert got == r["out"], (cell, r)
            n += 1
    assert n == 459
