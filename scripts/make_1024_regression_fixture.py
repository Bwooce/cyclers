"""#1024: record resonant_circle and half_rev_vectors outputs of the current (Keplerian) cells.

The pinned regression test tests/search/test_1024_keplerian_regression.py asserts the edited code
returns these values bit-identically (pre-registration amendment A1,
docs/notes/2026-10-11-1024-barycentric-binary-model.md sec. 3a). Run ONCE, on the pre-edit code.

Usage: uv run python scripts/make_1024_regression_fixture.py
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

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "tests" / "search" / "fixtures" / "1024_keplerian_regression.json"
RESONANCES = [(1, 1), (2, 1), (1, 2), (3, 2), (2, 3)]
HALF_REVS = [(1, 0, True), (1, 0, False), (3, 1, True), (3, 1, False)]


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def cells() -> dict[str, tuple[Any, list[str], list[float]]]:
    enum = _load("run_942_enumerate", REPO / "scripts" / "run_942_enumerate.py")
    drv = _load("run_998_enumerate", REPO / "scripts" / "run_998_enumerate.py")
    out = {}
    for cell, vinfs in (("ev", [3.0, 6.0, 10.0]), ("gc", [1.5, 3.0, 6.0]), ("vm", [8.0, 12.0])):
        system, a, b = enum.cell_system(cell)
        out[cell] = (system, [c for c in (a, b) if c not in system.massless], vinfs)
    out["pn"] = (drv.pluto_system("Nix"), ["Charon"], [0.15, 0.25, 0.35])
    return out


def record() -> dict[str, Any]:
    rec: dict[str, Any] = {}
    for cell, (system, bodies, vinfs) in cells().items():
        rows = []
        for body in bodies:
            per = system.period_s(body)
            for t in (0.0, 0.37 * per, 1.91 * per):
                for v in vinfs:
                    for n, m in RESONANCES:
                        c = resonant_circle(system, ResonantLeg(body, n, m), t, v)
                        rows.append(
                            {"kind": "R", "body": body, "nm": [n, m], "t": t, "v": v,
                             "out": None if c is None else [c[0], c[1], *map(list, c[2:])]}
                        )  # fmt: skip
                    for h, k, p in HALF_REVS:
                        vs = half_rev_vectors(system, HalfRevLeg(body, h, k, p, 0), t, v)
                        rows.append(
                            {"kind": "H", "body": body, "hkp": [h, k, p], "t": t, "v": v,
                             "out": [list(map(float, u)) for u in vs]}
                        )  # fmt: skip
        rec[cell] = rows
    return rec


def main() -> None:
    rec = record()
    OUT.write_text(json.dumps(rec, default=float))
    print(f"wrote {OUT.relative_to(REPO)} ({sum(len(v) for v in rec.values())} rows)")


if __name__ == "__main__":
    main()
