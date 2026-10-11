"""DOP853 re-fly of a solved two-working-body cycle (independent of the Lambert solver).

Moved out of scripts/gauntlet_942.py (#1038) so that a task that must not load
``search.literature_check`` can reuse the re-fly: this module imports nothing from it.
``scripts/gauntlet_942.py`` re-exports ``fly`` and ``cross_check`` under the same names.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from cyclerfinder.search.two_working_body import (
    Flyby,
    _blocks,
    eval_lambert_legs,
    fixed_duration_s,
    gate_cycle,
)


def _twobody(_t: float, y: np.ndarray, mu: float) -> np.ndarray:
    r = y[:3]
    return np.concatenate([y[3:], -mu * r / float(np.linalg.norm(r)) ** 3])


def fly(mu: float, r0: np.ndarray, v0: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    sol = solve_ivp(
        _twobody,
        (0.0, dt),
        np.concatenate([r0, v0]),
        args=(mu,),
        method="DOP853",
        rtol=1e-13,
        atol=1e-8,
    )
    return sol.y[:3, -1], sol.y[3:, -1]


def cross_check(system: Any, cycle: Any, x: np.ndarray, flybys: list[Flyby]) -> dict[str, Any]:
    """Re-fly every leg with DOP853 and rebuild the junctions and the gate."""
    legs = eval_lambert_legs(system, cycle, x)
    assert legs is not None
    worst_miss = 0.0
    worst_dv = 0.0
    integrated_in: dict[float, np.ndarray] = {}
    for ev, li in zip(legs, cycle.lambert_index, strict=True):
        leg = cycle.legs[li]
        r0, w0 = system.state(leg.frm, ev.t_dep)
        r1, v1 = fly(system.mu, r0, w0 + ev.vinf_dep, ev.t_arr - ev.t_dep)
        rb, wb = system.state(leg.to, ev.t_arr)
        worst_miss = max(worst_miss, float(np.linalg.norm(r1 - rb)))
        worst_dv = max(worst_dv, float(np.linalg.norm((v1 - wb) - ev.vinf_arr)))
        integrated_in[round(ev.t_arr, 3)] = v1 - wb
    # fixed legs: start from each chosen direction, compare with the next flyby's inbound
    for blk in _blocks(system, cycle, legs):
        for leg, t0 in blk.fixed:
            fb_out = next(f for f in flybys if abs(f.t_s - t0) < 1.0 and f.body == blk.body)
            dt = fixed_duration_s(system, leg)
            r0, w0 = system.state(blk.body, t0)
            r1, v1 = fly(system.mu, r0, w0 + fb_out.vinf_out, dt)
            rb, wb = system.state(blk.body, t0 + dt)
            worst_miss = max(worst_miss, float(np.linalg.norm(r1 - rb)))
            nxt = next(f for f in flybys if abs(f.t_s - (t0 + dt)) < 1.0 and f.body == blk.body)
            worst_dv = max(worst_dv, float(np.linalg.norm((v1 - wb) - nxt.vinf_in)))
            integrated_in[round(t0 + dt, 3)] = v1 - wb
    # gate on integrated inbound vectors where available
    rebuilt = []
    for f in flybys:
        vin = integrated_in.get(round(f.t_s, 3))
        if vin is None:
            # the wrap arrival is stored at t + period
            vin = integrated_in.get(round(f.t_s + cycle.period_s, 3), f.vinf_in)
        rebuilt.append(Flyby(f.body, f.t_s, np.asarray(vin), f.vinf_out))
    rep = gate_cycle(system, rebuilt)
    mism = max(
        (
            abs(float(np.linalg.norm(f.vinf_in)) - float(np.linalg.norm(f.vinf_out)))
            for f in rebuilt
        ),
        default=0.0,
    )
    return {
        "max_arrival_miss_km": worst_miss,
        "max_vinf_vector_error_kms": worst_dv,
        "max_junction_mismatch_kms": mism,
        "gate_status_integrated": rep.status,
        "worst_ratio_integrated": rep.gate.worst_ratio,
    }
