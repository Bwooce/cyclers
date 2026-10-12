"""#1055 sec. 14.2 diagnostic (NOT pre-registered; lead-approved): crank resolution of gc4-1's
down points in the ideal model. For each stored point of data/1055_vinf_limit/gc4-1.json: the
V_inf-angle effect of the Newton move along each of the three weakest column-scaled directions,
and the angles after a 10 + 10 iteration polish. Read-only.

    uv run python scripts/diag_1055_gc4_1_crank_resolution.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

from cyclerfinder.nbody.degenerate_family import column_scales, scaled_svd

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


C = _load("check_1055_vinf_direction", ROOT / "scripts" / "check_1055_vinf_direction.py")


def main() -> None:
    b = C.B
    p = b.build("gc4-1")
    moons = list(p.moons)
    chain = C.vinf_vectors(moons, np.asarray(p.seed), 1.0)

    def ang(z: Any, s: float) -> Any:
        vv = C.vinf_vectors(moons, z, s)
        return np.array([C._ang(vv[k][j], chain[k][j]) for k in range(4) for j in (0, 1)])

    rec = json.loads((ROOT / "data" / "1055_vinf_limit" / "gc4-1.json").read_text())
    for pt in rec["points"][1:]:
        s, z = pt["sigma"], np.asarray(pt["z"])
        p.set_sigma(s)
        r, j, _ = b.residual_and_jac(p, z)
        d = column_scales(j)
        u, sv, vt = scaled_svd(j, d)
        a0 = ang(z, s)
        eff = []
        for i in (1, 2, 3):
            dz = vt[-i] * (-float(u[:, -i] @ r) / float(sv[-i])) / d
            eff.append((float(sv[-i]), float(np.max(np.abs(ang(z + dz, s) - a0)))))
        zn, inf = b.solve(p, z, 10, 10, False)
        a1 = ang(zn, s)
        print(
            f"sigma {s}: |r| {np.linalg.norm(r):.2e} G2out/G3in {a0[5]:.4f}/{a0[6]:.4f} "
            f"weak-direction effects (sv, deg) {eff}; polished {a1[5]:.4f}/{a1[6]:.4f} "
            f"converged {b.converged(inf)}",
            flush=True,
        )


if __name__ == "__main__":
    main()
