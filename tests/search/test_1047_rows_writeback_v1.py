"""#1047 rows A, B, D (owner decision 2026-10-11, "Rows for A, B, D only"): the V1 evidence,
recomputed from data/1047_rows/writeback_v1.json (scripts/writeback_1047_rows.py).

For each row's corrected state (planar CR3BP, the project's registry Earth-Moon mu):
- the corrector's own propagation (DOP853 with the STM, as in correct_periodic) closes to < 1e-10,
  the corrector tolerance;
- an independent Radau re-propagation keeps the Jacobi constant to < 1e-8 over one period, and
  agrees with DOP853 arc by arc (4p equal arcs, each Radau arc from the DOP853 node) to < 1e-8
  (crosscheck_periodic's tolerance). This is the multiple-shooting form of the independent check:
  these orbits have |lambda| ~ 1e2-1e3, so a whole-period closure depends on the integrator's step
  sequence at the 1e-7 to 1e-5 level (a plain DOP853 run without the STM; recorded in the evidence
  file as closure_plain_dop853), which is error amplification, not a closure failure;
- the orbit is cycler class: perigee between the Earth floor and the GEO radius, periselene between
  the Moon floor and the Hill radius (registry radii and safe altitudes; GEO 42,164 km);
- it is asymmetric (no perpendicular x-axis crossing; the mirror image is a different orbit) and its
  inertial p:q and sense are the family's;
- once the rows are in the catalogue, each row's state, C and period are the evidence values.

The thresholds are model and project conventions, not values computed by this code.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import yaml  # type: ignore[import-untyped]
from scipy.integrate import solve_ivp

from cyclerfinder.core import cr3bp
from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.satellites import SATELLITES

REPO = Path(__file__).resolve().parents[2]
EVIDENCE = json.loads((REPO / "data" / "1047_rows" / "writeback_v1.json").read_text())
SYSTEM = cr3bp.cr3bp_system("Earth", "Moon")
FLOOR_E_KM = PLANETS["E"].radius_eq_km + PLANETS["E"].safe_alt_km
FLOOR_M_KM = SATELLITES["Moon"].radius_eq_km + SATELLITES["Moon"].safe_alt_km
GEO_KM = 42164.0
HILL_KM = (SYSTEM.mu / 3.0) ** (1.0 / 3.0) * SYSTEM.l_km
EXPECTED = {"A": ((5, 2), "prograde"), "B": ((5, 3), "retrograde"), "D": ((7, 3), "prograde")}


def _base() -> Any:
    name = "run_1000_complement"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, REPO / "scripts" / f"{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("fam", sorted(EVIDENCE))
def test_closure_and_independent_integrator(fam: str) -> None:
    e = EVIDENCE[fam]
    assert e["mass_ratio"] == SYSTEM.mu
    s0 = np.array(e["state_nd"])
    period = float(e["period_nd"])
    arc = cr3bp.propagate(SYSTEM, s0, period, with_stm=True)
    assert float(np.linalg.norm(np.asarray(arc.state_f) - s0)) < 1e-10
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        s0,
        args=(SYSTEM.mu,),
        method="Radau",
        rtol=1e-11,
        atol=1e-11,
    )
    c0 = cr3bp.jacobi_constant(s0, SYSTEM.mu)
    assert abs(cr3bp.jacobi_constant(sol.y[:, -1], SYSTEM.mu) - c0) < 1e-8
    assert c0 == pytest.approx(e["jacobi_constant"], abs=1e-12)
    n = 4 * EXPECTED[fam][0][0]
    s = s0.copy()
    for _ in range(n):
        d = cr3bp.propagate(SYSTEM, s, period / n, rtol=1e-12, atol=1e-12)
        r = solve_ivp(
            cr3bp.cr3bp_eom,
            (0.0, period / n),
            s,
            args=(SYSTEM.mu,),
            method="Radau",
            rtol=1e-11,
            atol=1e-11,
        )
        assert float(np.linalg.norm(r.y[:, -1] - np.asarray(d.state_f))) < 1e-8
        s = np.asarray(d.state_f)


@pytest.mark.parametrize("fam", sorted(EVIDENCE))
def test_cycler_class_asymmetric_and_resonance(fam: str) -> None:
    e = EVIDENCE[fam]
    s0 = np.array(e["state_nd"])
    s4 = np.array([s0[0], s0[1], s0[3], s0[4]])
    period = float(e["period_nd"])
    cl = _base().classify(s4, period)
    assert FLOOR_E_KM <= cl["perigee_km"] <= GEO_KM
    assert FLOOR_M_KM <= cl["periselene_km"] <= HILL_KM
    assert not cl["symmetric"]
    n_in = cl["wind_E"] + period / (2 * math.pi)
    pq, sense = EXPECTED[fam]
    assert (abs(round(n_in)), round(period / (2 * math.pi))) == pq
    assert ("prograde" if n_in > 0 else "retrograde") == sense
    assert e["mirror_min_distance_nd"] > 1e-3


@pytest.mark.parametrize("fam", sorted(EVIDENCE))
def test_catalogue_row_matches_evidence(fam: str) -> None:
    e = EVIDENCE[fam]
    rows = {r["id"]: r for r in yaml.safe_load((REPO / "data" / "catalogue.yaml").read_text())}
    if e["id"] not in rows:
        pytest.skip("row not yet written back (catalogue window)")
    row = rows[e["id"]]
    cr = row["orbit_elements"]["cr3bp"]
    assert cr["state_nd"] == e["state_nd"]
    assert cr["jacobi_constant"] == e["jacobi_constant"]
    assert cr["period_nd"] == e["period_nd"]
    assert row["our_status"] == "known-class-member"
    assert row["validation_level"] == "V1"
