"""#1007 ev positive control, take 3 (note sec. 8.2): H&M 1970 orbit 1 from its printed dates.

The 5-cycle chain of key k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s on Standish mean elements is solved
directly (Levenberg-Marquardt on the chain tool's date residual) from the Table 3 dates, then judged
by scripts/v2_chain_1007.py under the sec. 7 amended rules. The ev negative control (ev-C epoch 0
with the 2nd Lambert start + 1 d) is judged in the same run.

Usage: uv run python scripts/v2_hm_printed_1007.py
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from cyclerfinder.search.hollister_menning_1970 import load_table3

REPO = Path(__file__).resolve().parents[1]
DAY = 86400.0
KEY = "k2|RE/1:1|LE>V/0s|RV/1:1|RV/1:1|LV>E/0s"
X_IDEAL = [474.59719334719347, 1100.6183800415804]
N = 5
OUT = REPO / "data" / "1007_v2" / "hm_orbit1_printed.json"


def _load(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


V2 = _load("v2_chain_1007", REPO / "scripts" / "v2_chain_1007.py")
CHAIN = V2.CHAIN


def main() -> None:
    V2.RULES = "amended"
    real = CHAIN.real_ephemeris("ev", "mean")
    rows = load_table3()[1]
    starts = []
    for k in range(N):
        starts += [rows[5 * k + 1].date, rows[5 * k + 4].date]
    x0 = starts[0]
    y_seed = np.array([*starts[1:], 5844.0])
    circ, a, b = CHAIN.ENUM.cell_system("ev")
    _, one = CHAIN.ENUM.parse_cycle_key(KEY, circ, a, b)
    legs = one.legs * N
    r0 = CHAIN.date_chain_residual(real, legs, x0 * DAY, y_seed)
    sol = least_squares(
        lambda yy: CHAIN.date_chain_residual(real, legs, x0 * DAY, yy),
        y_seed,
        method="lm",
        xtol=1e-14,
        ftol=1e-14,
        gtol=1e-14,
        max_nfev=200 * len(y_seed),
    )
    r1 = CHAIN.date_chain_residual(real, legs, x0 * DAY, sol.x)
    shift = float(np.max(np.abs(sol.x[:-1] - y_seed[:-1])))
    print(
        f"seed residual {np.max(np.abs(r0)):.3e} km/s -> solved {np.max(np.abs(r1)):.3e} km/s; "
        f"max date shift from print {shift:.3f} d; chain period {sol.x[-1]:.3f} d (print 5844)",
        flush=True,
    )
    pos = V2.judge("ev", KEY, X_IDEAL, N, "ramp", x0, list(sol.x), real, "CONTROL+ H&M orbit 1")
    pos |= {
        "seed_residual_kms": float(np.max(np.abs(r0))),
        "solved_residual_kms": float(np.max(np.abs(r1))),
        "max_shift_from_print_d": shift,
        "chain_period_d": float(sol.x[-1]),
        "x0_jd": 2440000.0 + x0,
        "printed_starts_days": starts,
        "solved_y_days": sol.x.tolist(),
    }
    direct = json.loads((REPO / "data" / "942_direct_route_and_gm_fix.json").read_text())
    e0 = direct["evC_standish_rerun"]["epochs"][0]
    wb = json.loads((REPO / "data" / "942_943_writeback_v1.json").read_text())["ev-C"]
    y_bad = list(e0["final_y_dates"])
    y_bad[0] += 1.0
    neg = V2.judge("ev", wb["key"], wb["x_days"], 5, "ramp", e0["x0_days"], y_bad, real,
                   "CONTROL- ev-C e0 broken")  # fmt: skip
    print(
        f"CONTROL+ {'PASS' if pos['pass'] else 'FAIL'}; CONTROL- "
        f"{'FAILS as required' if not neg['pass'] else 'PASSED (BAD)'}",
        flush=True,
    )
    OUT.write_text(json.dumps({"positive": pos, "negative": neg}, indent=1, default=float))
    print("wrote", OUT.relative_to(REPO), flush=True)


if __name__ == "__main__":
    main()
