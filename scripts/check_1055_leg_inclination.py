"""#1055 sec. 14.1 follow-up (read-only): per-leg inclination of each gc member's chain-tool
patched-conic member in the ideal model, against the catalogue's row-level inclination_deg.

Each node's outbound leg: the spacecraft state is the moon's state plus the outbound V_inf (the
two-body asymptote of the seed node at full mass); the leg inclination is the angle of r x v from
the moons' plane normal. The moon's position (not the periapsis node) is used, so an out-of-plane
periapsis offset does not tilt an in-plane leg.

    uv run python scripts/check_1055_leg_inclination.py [member ...]

Writes data/1055_leg_inclination.json when run over all members.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MEMBERS = [
    "gc4-19",
    "gc-1",
    "gc-2",
    "gc4-1",
    "gc4-2",
    "gc4-6",
    "gc4-8",
    "gc4-10",
    "gc4-25",
    "gc4-28",
    "gc4-38",
    "gc4-41",
    "gc5-3",
    "gc5-5",
    "gc5-6",
    "gc6-0",
    "gc6-1",
    "gc6-2",
]


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


C = _load("check_1055_vinf_direction", ROOT / "scripts" / "check_1055_vinf_direction.py")


def legs(name: str) -> list[dict[str, Any]]:
    p = C.B.build(name)
    z = np.asarray(p.seed)
    moons = list(p.moons)
    vv = C.vinf_vectors(moons, z, 1.0)
    out = []
    for k, moon in enumerate(moons):
        rm, vm = p.circ.state(moon, float(z[7 * k + 6]))
        h = np.cross(np.asarray(rm), np.asarray(vm) + vv[k][1])
        inc = math.degrees(math.acos(min(1.0, abs(float(h[2])) / float(np.linalg.norm(h)))))
        nxt = moons[(k + 1) % len(moons)]
        out.append({"leg": f"{moon}{k} -> {nxt}", "inclination_deg": inc})
    return out


def main() -> None:
    names = sys.argv[1:] or MEMBERS
    res = {}
    for name in names:
        lg = legs(name)
        res[name] = {"legs": lg, "max_deg": max(x["inclination_deg"] for x in lg)}
        print(
            name,
            f"max {res[name]['max_deg']:.3f}",
            [(x["leg"], round(x["inclination_deg"], 3)) for x in lg],
        )
    if not sys.argv[1:]:
        (ROOT / "data" / "1055_leg_inclination.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
