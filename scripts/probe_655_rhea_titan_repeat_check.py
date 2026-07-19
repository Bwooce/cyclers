"""#655 -- multi-cycle repeat-instrumentation control for the Rhea-Titan
inclination-extension closing basins found by
``scripts/probe_655_rhea_titan_3d_closure.py``.

Mirrors #575's own C2 discipline: a single-cycle 3D Lambert closure at the
real mutual inclination is NOT sufficient evidence of a genuine repeating
cycler (#575's own Titan-Iapetus result found 6-9 single-cycle-closing
basins that then ALL failed to repeat over 3 cycles -- see
``data/empty_regions.jsonl`` region
``saturn-titan-iapetus-symmetric-closure-inclination-empty-575``). This
script re-solves the SAME (Omega, tof_scale, n_rev, rel_offset) Lambert-leg
geometry each of #655's inclination-extension closing basins converged to,
across 3 consecutive cycles, with BOTH Rhea (kept equatorial, per
#572/#655's own "anchor stays in-plane" convention) and Titan (on the fixed
inclined circular orbit found by the node-alignment search) advanced through
their natural Keplerian motion -- exactly the #575 C2 / #574-Stage-B
"same-model re-solve over cycles" mechanism, built directly on
``probe_572_titan_iapetus_3d_closure.py``'s own ``iapetus_state_3d``/
``_leg_best`` primitives (not reimplemented) rather than the Titan-Iapetus-
specific ``v2_saturn_3d``/``titan_iapetus_corrector`` modules (those hardcode
``ANCHOR="Titan"``/``FLYBY="Iapetus"`` internals and are not a drop-in fit
for a different pair without a rewrite -- reusing the lower-level
``iapetus_state_3d``/``_leg_best`` primitives instead keeps this a reuse, not
a new corrector build).

PASS criterion (same discipline as #572/#575, NOT relaxed): every cycle's
Lambert leg must converge to the SAME requested n_rev, the V_inf-continuity
residual must stay under ``GATE_RESIDUAL_KMS`` (the #558-lineage project-wide
bar) for all 3 cycles, and the physical bend gate must keep passing at every
cycle.

Discipline: NO catalogue writeback, no new gate logic, read-only on
``data/probe_655_rhea_titan_3d_closure.jsonl``.

Run as::

    uv run python scripts/probe_655_rhea_titan_repeat_check.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import probe_572_titan_iapetus_3d_closure as p572  # noqa: E402
from scan_558_uranus_all_pairs_offset_sweep import GATE_RESIDUAL_KMS  # noqa: E402

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.discovery_campaign import (  # noqa: E402
    DAY_S,
    _mean_motion_rad_day,
    _moon_state,
)
from cyclerfinder.search.physical_sanity import (  # noqa: E402
    DEFAULT_MIN_USEFUL_BEND_DEG,
    candidate_passes_physical_gate,
)

__all__ = [
    "DEFAULT_MIN_USEFUL_BEND_DEG",
    "GATE_RESIDUAL_KMS",
    "candidate_passes_physical_gate",
    "main",
    "p572",
]

p572.PRIMARY = "Saturn"
p572.ANCHOR = "Rhea"
p572.FLYBY = "Titan"
p572.INCLINATION_DEG = 0.7

DATA_DIR = ROOT / "data"
IN_PATH = DATA_DIR / "probe_655_rhea_titan_3d_closure.jsonl"
OUT_PATH = DATA_DIR / "probe_655_rhea_titan_repeat_check.jsonl"

N_CYCLES = 3


def load_closing_basins() -> list[dict[str, Any]]:
    """The best-closing basin from each #655 candidate verdict where
    ``closure_found`` is True, loaded directly from the committed 3D-closure
    probe output (never hand-transcribed)."""
    out = []
    with IN_PATH.open(encoding="utf-8") as fh:
        for line in fh:
            d = json.loads(line)
            if d.get("kind") != "candidate_verdict":
                continue
            if not d.get("closure_found"):
                continue
            best = d["best_closing_basin"]
            out.append({"label": d["label"], **best})
    return out


def repeat_check(
    cand: dict[str, Any], rel_offset_deg: float, n_rev: tuple[int, int]
) -> dict[str, Any]:
    mu = PRIMARIES["Saturn"]
    sat_a = SATELLITES["Rhea"]
    sat_b = SATELLITES["Titan"]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    v_circ_b = math.sqrt(mu / sat_b.sma_km)

    tof_scale = cand["tof_scale"]
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b
    tof_days = tof_scale * math.sqrt(p_a * p_b)
    cycle_period_days = 2.0 * tof_days

    omega = math.radians(cand["omega_deg"])
    inc = math.radians(p572.INCLINATION_DEG)
    u0 = math.radians(rel_offset_deg)
    n0, n1 = n_rev

    per_cycle: list[dict[str, Any]] = []
    max_residual = 0.0
    max_drift_km = 0.0
    r2_cycle0 = None

    for k in range(N_CYCLES):
        t0 = k * cycle_period_days
        r0, v0 = _moon_state(0.0, n_a, t0, sat_a.sma_km, mu)
        u_titan = u0 + n_b * (t0 + tof_days)
        r1, v1 = p572.iapetus_state_3d(u_titan, v_circ_b, sat_b.sma_km, omega, inc)
        r2, v2 = _moon_state(0.0, n_a, t0 + 2.0 * tof_days, sat_a.sma_km, mu)

        tof_s = tof_days * DAY_S
        try:
            leg0 = p572._leg_best(r0, v0, r1, v1, tof_s, mu, n0)
            leg1 = p572._leg_best(r1, v1, r2, v2, tof_s, mu, n1)
        except Exception as exc:  # LambertGeometryError/ConvergenceError
            per_cycle.append({"cycle_index": k, "converged": False, "error": str(exc)})
            continue

        if leg0 is None or leg1 is None:
            per_cycle.append({"cycle_index": k, "converged": False, "error": "n_rev infeasible"})
            continue

        r_mid = abs(leg0["vinf_in"] - leg1["vinf_out"])
        r_periodic = abs(leg0["vinf_out"] - leg1["vinf_in"])
        residual = max(r_mid, r_periodic)
        max_residual = max(max_residual, residual)

        vinf_tuple = (leg0["vinf_out"], max(leg0["vinf_in"], leg1["vinf_out"]), leg1["vinf_in"])
        seq = ("Rhea", "Titan", "Rhea")
        gate_pass, gate_verdicts = candidate_passes_physical_gate(
            seq, vinf_tuple, min_useful_bend_deg=DEFAULT_MIN_USEFUL_BEND_DEG
        )
        bends = [v.max_bend_deg for v in gate_verdicts]

        if k == 0:
            r2_cycle0 = r2
        else:
            assert r2_cycle0 is not None
            import numpy as np

            drift = float(np.linalg.norm(r2 - r2_cycle0))
            max_drift_km = max(max_drift_km, drift)

        per_cycle.append(
            {
                "cycle_index": k,
                "converged": True,
                "residual_kms": residual,
                "vinf_kms": list(vinf_tuple),
                "max_bend_deg_per_encounter": bends,
                "physical_gate_pass": gate_pass,
            }
        )

    n_cycles_completed = sum(1 for c in per_cycle if c.get("converged"))
    all_residuals_ok = all(
        c.get("residual_kms", math.inf) < GATE_RESIDUAL_KMS for c in per_cycle if c.get("converged")
    )
    all_bend_ok = all(c.get("physical_gate_pass", False) for c in per_cycle if c.get("converged"))
    repeats = bool(n_cycles_completed == N_CYCLES and all_residuals_ok and all_bend_ok)

    return {
        "label": cand["label"],
        "rel_offset_deg": rel_offset_deg,
        "n_rev": list(n_rev),
        "omega_deg": cand["omega_deg"],
        "tof_scale": tof_scale,
        "cycle0_residual_kms": cand["residual_kms"],
        "n_cycles_requested": N_CYCLES,
        "n_cycles_completed": n_cycles_completed,
        "max_residual_kms": max_residual,
        "max_drift_km": max_drift_km,
        "per_cycle": per_cycle,
        "repeats_as_genuine_cycle": repeats,
    }


def main() -> int:
    basins = load_closing_basins()
    print(f"[655-repeat] {len(basins)} inclination-extension closing basins to test", flush=True)

    # rel_offset_deg / n_rev per label, parsed from the label string
    # (format "n{n}_nrev[{n0}, {n1}]_rel{rel}") -- also cross-checked against
    # the original #655 enumeration record for the same label.
    import re

    enum_by_label: dict[str, dict[str, Any]] = {}
    with (DATA_DIR / "enumerate_655_saturn_rhea_titan_symmetric_closures.jsonl").open() as fh:
        for line in fh:
            d = json.loads(line)
            if d.get("kind") == "pass" and d.get("anchor") == "Rhea":
                lbl = f"n{d['n_commensurate_int']}_nrev{d['n_rev']}_rel{d['rel_offset_deg']:.0f}"
                enum_by_label[lbl] = d

    results = []
    n_repeat = 0
    for b in basins:
        m = re.match(r"n(\d+)_nrev\[(\d+), (\d+)\]_rel(\d+)", b["label"])
        assert m is not None, b["label"]
        n_rev = (int(m.group(2)), int(m.group(3)))
        rel_offset_deg = enum_by_label[b["label"]]["rel_offset_deg"]
        res = repeat_check(b, rel_offset_deg, n_rev)
        results.append(res)
        if res["repeats_as_genuine_cycle"]:
            n_repeat += 1
        print(
            f"[655-repeat]   {res['label']}: n_cycles_completed={res['n_cycles_completed']}/"
            f"{N_CYCLES} max_residual_kms={res['max_residual_kms']:.4f} "
            f"max_drift_km={res['max_drift_km']:.1f} "
            f"repeats={res['repeats_as_genuine_cycle']}",
            flush=True,
        )

    print(f"[655-repeat] {n_repeat}/{len(basins)} basins repeat as genuine 3D cycles", flush=True)

    with OUT_PATH.open("w", encoding="utf-8") as fh:
        fh.write(
            json.dumps(
                {
                    "_meta": True,
                    "task": "#655 Rhea-Titan inclination-extension multi-cycle repeat check",
                    "n_cycles": N_CYCLES,
                    "inclination_deg": p572.INCLINATION_DEG,
                    "gate_residual_kms": GATE_RESIDUAL_KMS,
                    "n_basins_tested": len(basins),
                    "n_basins_repeat": n_repeat,
                }
            )
            + "\n"
        )
        for res in results:
            fh.write(json.dumps({"kind": "repeat_result", **res}) + "\n")
    print(f"[655-repeat] written: {OUT_PATH}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
