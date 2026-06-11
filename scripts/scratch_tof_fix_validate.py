#!/usr/bin/env python3
"""Scratch: close + V1-gate + n-body-confirm the 6 descriptor-bearing russell-ch4 rows.

For each row: close the G leg on real DE440 with the joint (epoch, ToF) free-variable
closer (the 2026-06-10 ToF-artifact fix), then
  V1: re-solve the SAME (r1, r2, tof) endpoints with lamberthub izzo + gooding and
      check agreement < 1e-3 m/s vs the in-house Lambert, plus a Kepler reprop of the
      arc back to the Mars-arrival position (independent of the Lambert that built it);
  V3: independent REBOUND/IAS15 Sun-only propagation of the found departure state to
      the Mars-arrival epoch (the #167 arbiter), recording the geometric miss + the
      emerged Mars v_inf in-band.

Emerged E/M v_inf are EVIDENCE compared to the SOURCED anchors (the golden). NO
writeback. Prints a per-row table consumed by the results note.
"""

from __future__ import annotations

import warnings
from datetime import UTC, datetime, timedelta
from pathlib import Path

import numpy as np
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import propagate as kepler_propagate
from cyclerfinder.core.lambert import lambert
from cyclerfinder.nbody.propagator import RestrictedNBody
from cyclerfinder.search.self_seeding import FamilyAnchors, joint_epoch_tof_close, on_family

J2000 = datetime(2000, 1, 1, 12, 0, 0, tzinfo=UTC)
SPD = 86400.0
IDS = [
    "russell-ch4-9.353Gg2",
    "russell-ch4-9.94Gg3",
    "russell-ch4-3.78Gg3",
    "russell-ch4-6.44Gg3",
    "russell-ch4-3.64gGg3",
    "russell-ch4-5.30ggF3",
]

# 3-Mars-SOI band, the #165/#167 constant (NEVER loosened).
_MARS_SOI_AU = PLANETS["M"].sma_au * (PLANETS["M"].mu_km3_s2 / MU_SUN_KM3_S2) ** 0.4
BAND_AU = 3.0 * _MARS_SOI_AU
V1_TOL_MPS = 1.0e-3


def v1_crosscheck(ephem: Ephemeris, t_depart: float, t_arrive: float) -> dict[str, float | bool]:
    """lamberthub izzo+gooding agreement + Kepler reprop on the closing leg endpoints."""
    from lamberthub import gooding1990, izzo2015  # type: ignore[import-untyped]

    r_e, _v_e = ephem.state("E", t_depart)
    r_m, _v_m = ephem.state("M", t_arrive)
    r1 = np.asarray(r_e, dtype=np.float64)
    r2 = np.asarray(r_m, dtype=np.float64)
    tof = t_arrive - t_depart

    sols = lambert(r1, r2, tof, mu=MU_SUN_KM3_S2, prograde=True, max_revs=0)
    mine = next(s for s in sols if s.n_revs == 0)
    mine_v1 = np.asarray(mine.v1, dtype=np.float64)
    v1_izzo, _ = izzo2015(MU_SUN_KM3_S2, r1, r2, tof, M=0, prograde=True, low_path=True)
    v1_good, _ = gooding1990(MU_SUN_KM3_S2, r1, r2, tof, M=0, prograde=True, low_path=True)
    d_izzo = float(np.linalg.norm(mine_v1 - np.asarray(v1_izzo)))
    d_good = float(np.linalg.norm(mine_v1 - np.asarray(v1_good)))
    max_diff_mps = 1000.0 * max(d_izzo, d_good)

    # Kepler reprop (not the Lambert that built it): does the arc reach Mars's position?
    r_end, _v_end = kepler_propagate(r1, mine_v1, tof)
    kepler_miss_au = float(np.linalg.norm(np.asarray(r_end) - r2) / AU_KM)

    return {
        "max_diff_mps": max_diff_mps,
        "kepler_miss_au": kepler_miss_au,
        "passed": bool(max_diff_mps < V1_TOL_MPS and kepler_miss_au < 1e-6),
    }


def nbody_confirm(
    ephem: Ephemeris, prop: RestrictedNBody, t_depart: float, t_arrive: float, vinf_vec: np.ndarray
) -> dict[str, float | bool]:
    """Independent REBOUND/IAS15 Sun-only confirm of the found departure state (#167 arbiter)."""
    r_e, v_e = ephem.state("E", t_depart)
    r0 = np.asarray(r_e, dtype=np.float64)
    v0 = np.asarray(v_e, dtype=np.float64) + np.asarray(vinf_vec, dtype=np.float64)
    out = prop.propagate(r0, v0, t_depart, t_arrive, accuracy=1e-11)
    r_m, v_m = ephem.state("M", t_arrive)
    miss_au = float(np.linalg.norm(out.r_km - np.asarray(r_m)) / AU_KM)
    vinf_m = float(np.linalg.norm(out.v_km_s - np.asarray(v_m)))
    return {
        "converged": bool(out.converged),
        "energy_drift": float(out.energy_rel_drift),
        "miss_au": miss_au,
        "vinf_m_kms": vinf_m,
        "in_band": bool(out.converged and miss_au < BAND_AU),
    }


def main() -> None:
    warnings.filterwarnings("ignore")
    catalogue = yaml.safe_load(Path("data/catalogue.yaml").read_text())
    rows = {str(r["id"]): r for r in catalogue}
    ephem = Ephemeris("astropy")
    prop = RestrictedNBody("rebound")
    t_center = (datetime(2027, 1, 1, tzinfo=UTC) - J2000).total_seconds()

    print(f"3-SOI band = {BAND_AU:.5f} AU ; V1 tol = {V1_TOL_MPS} m/s")
    for rid in IDS:
        row = rows[rid]
        vinf = {e["body"]: float(e["vinf_kms"]) for e in (row.get("vinf_kms_at_encounters") or [])}
        transit = (row.get("invariants") or {}).get("transit_times_days")
        sigs = sorted({float(t) for t in (transit or [])})
        anchors = FamilyAnchors(vinf_e=vinf["E"], vinf_m=vinf["M"], vinf_band_kms=1.5)
        print(f"\n=== {rid}  sourced E={vinf['E']} M={vinf['M']}  sig={transit} ===")
        for sig in sigs:
            r = joint_epoch_tof_close(ephem, anchors, t_center, sig, max_revs=2)
            if r is None:
                print(f"  sig={sig:.0f}d  NO LAMBERT SOLUTION (negative)")
                continue
            v = on_family(r, anchors)
            v1 = v1_crosscheck(ephem, r.t_depart_sec, r.t_arrive_sec)
            nb = nbody_confirm(ephem, prop, r.t_depart_sec, r.t_arrive_sec, r.vinf_vec)
            depart = (J2000 + timedelta(seconds=r.t_depart_sec)).date().isoformat()
            print(
                f"  sig={sig:.0f}d tof={r.tof_g_days:.1f}d depart={depart} "
                f"(t={r.t_depart_sec:.0f}s)  "
                f"E={r.vinf_e_kms:6.2f}(d{r.vinf_e_kms - vinf['E']:+.2f}) "
                f"M={r.vinf_m_kms:6.2f}(d{r.vinf_m_kms - vinf['M']:+.2f})  "
                f"on_family={v.on_family}\n"
                f"        V1: izzo/gooding={v1['max_diff_mps']:.2e}m/s "
                f"kepler_miss={v1['kepler_miss_au']:.2e}AU PASS={v1['passed']}\n"
                f"        V3(nbody): conv={nb['converged']} drift={nb['energy_drift']:.1e} "
                f"miss={nb['miss_au']:.2e}AU vinfM={nb['vinf_m_kms']:.2f} in_band={nb['in_band']}"
            )


if __name__ == "__main__":
    main()
