"""#1059: a trial date outside the ephemeris kernel's span is a rejected step, not a crash.

The chain tool's shoot evaluates solver trial steps that can leave DE440's span (#1057 step 1b e2
reached a date past 2650-01-25 and the run died with jplephem's OutOfRangeError). `_res` must
return the 1e3 penalty and `_safe_eval` None. The synthetic system below raises the real
jplephem exception for any date past a cutoff.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
from jplephem.exceptions import OutOfRangeError

REPO = Path(__file__).resolve().parents[2]


def _load(name: str, rel: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # dataclasses resolve string annotations through sys.modules
    spec.loader.exec_module(mod)
    return mod


class _Limited:
    """Delegates to a real circular system but raises OutOfRangeError past ``cutoff_s``."""

    def __init__(self, inner: Any, cutoff_s: float) -> None:
        self._inner = inner
        self._cutoff = cutoff_s

    def state(self, code: str, t_s: float) -> Any:
        if t_s > self._cutoff:
            raise OutOfRangeError(
                "ephemeris segment only covers dates 1549-12-31 through 2650-01-25",
                np.array([t_s]),
            )
        return self._inner.state(code, t_s)

    def __getattr__(self, name: str) -> Any:
        return getattr(self._inner, name)


def _setup() -> tuple[Any, Any, tuple[Any, ...], float, np.ndarray]:
    chain = _load("run_942_realeph_chain_1059", "scripts/run_942_realeph_chain.py")
    enum = chain.ENUM
    system, a, b = enum.cell_system("ev")
    cand = next(
        c
        for c in json.loads((REPO / "data/942_cell_ev_gauntlet.json").read_text())["candidates"]
        if c["key"] == "k2|LE>V/0s|LV>V/1h|LV>E/0s"
    )
    _, cycle = enum.parse_cycle_key(cand["key"], system, a, b)
    x = list(cand["x_days"])
    period_d = cycle.period_s / 86400.0
    y = np.array([*x[1:], period_d])
    return chain, system, cycle.legs, x[0], y


def test_in_range_evaluation_is_not_penalised() -> None:
    chain, system, legs, x0, y = _setup()
    r = chain._res(system, legs, x0, y)
    assert np.all(np.abs(r) < 1e-6)  # the stored ev-C zero closes in the ideal model


def test_out_of_range_date_returns_penalty_and_does_not_raise() -> None:
    chain, system, legs, x0, y = _setup()
    limited = _Limited(system, cutoff_s=0.5 * 86400.0)  # every encounter date is past the cutoff
    r = chain._res(limited, legs, x0, y)
    assert r.shape[0] == len(legs) - sum(1 for lg in legs if not hasattr(lg, "nrev"))
    assert np.all(r == 1e3)
    assert chain._safe_eval(limited, legs, x0, y) is None
