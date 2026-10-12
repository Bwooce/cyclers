"""#975: Varin's closed families i_k.

Pre-registration: docs/notes/2026-10-11-975-varin-cascade-satellite-mass-ratios.md.

Expected values come from Bruno & Varin 2009, Solar System Research 43(3):253-276, Table 1 (eight
critical orbits printed to 8 digits) and Table 2 (their mu), as transcribed in
data/sources/varin-bruno-2009-closed-families-tables.yaml, and from model identities. Tolerances are
the pre-registration's control P1 (sec. 5).

Production scale: the Jacobian and geometry checks run at the registry mass ratios and length units
of the target systems (Jupiter-Ganymede, Saturn-Titan, Neptune-Triton), so no step or tolerance can
silently assume Varin's mu or a unit length.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
import pytest
import yaml  # type: ignore[import-untyped]

from cyclerfinder.core import cr3bp
from cyclerfinder.search import varin_closed_families as vcf

TABLES = "data/sources/varin-bruno-2009-closed-families-tables.yaml"


def _table1() -> dict[str, dict[str, Any]]:
    with open(TABLES) as fh:
        d = yaml.safe_load(fh.read())
    mus = {r["k"]: (float(r["mu_p"]), float(r["mu_pp"])) for r in d["table2_bifurcation_mu"]}
    out = {}
    for r in d["table1_bifurcation_orbits"]:
        lab = r["mu"]
        mu = mus[lab[2]][1] if lab.endswith("pp") else mus[lab[2]][0]
        x1, y2 = vcf.crossing_from_varin(float(r["a"]), float(r["e"]), mu)
        out[lab] = {
            "mu": mu,
            "x1": x1,
            "y2": y2,
            "T": 2.0 * math.pi * float(r["T_over_2pi"]),
            "C": float(r["C"]),
            "Trv": float(r["Trv"]),
        }
    return out


ROWS = _table1()


@pytest.mark.parametrize("lab", sorted(ROWS))
def test_printed_orbit_closes_and_jacobi_convention(lab: str) -> None:
    r = ROWS[lab]
    _, y1 = vcf.crossing_near_half(r["x1"], r["y2"], r["T"], r["mu"])
    assert y1 <= 1e-8
    cj = cr3bp.jacobi_constant(vcf.state_from_crossing(r["x1"], r["y2"], r["mu"]), r["mu"])
    assert cj == pytest.approx(vcf.jacobi_from_printed(r["C"], r["mu"]), abs=1e-8)


@pytest.mark.parametrize("lab", ["mu1p", "mu1pp"])
def test_traces_match_printed(lab: str) -> None:
    r = ROWS[lab]
    tr, trv = vcf.traces(r["x1"], r["y2"], r["T"], r["mu"])
    assert vcf.modified_trace(tr) == pytest.approx(2.0, abs=1e-4)
    assert vcf.modified_trace(trv) == pytest.approx(r["Trv"], abs=1e-6)


@pytest.mark.parametrize(("lab", "half"), [("mu1p", False), ("mu1p", True), ("mu4pp", True)])
def test_critical_point_solver_returns_to_printed(lab: str, half: bool) -> None:
    r = ROWS[lab]
    z = np.array([r["x1"], r["y2"], r["T"], r["mu"]])
    cp = vcf.solve_critical_point(z * (1.0 + 1e-4 * np.array([1.0, -1.0, 1.0, -1.0])), half=half)
    assert cp.converged
    zc = np.array([cp.x1, cp.y2, cp.t_full, cp.mu])
    assert float(np.max(np.abs((zc - z) / z))) <= 5e-8


@pytest.mark.parametrize(
    ("primary", "secondary"), [("Jupiter", "Ganymede"), ("Saturn", "Titan"), ("Neptune", "Triton")]
)
def test_half_jacobian_at_production_mass_ratios(primary: str, secondary: str) -> None:
    """The analytic half-period Jacobian agrees with central differences at the target mu, on a
    P1-circling orbit of radius 0.9 (no periodicity needed for a derivative check)."""
    system = cr3bp.cr3bp_system(primary, secondary)
    mu = system.mu
    x1, t_full = 0.9, 2.0 * math.pi * 1.2
    y2 = math.sqrt((1.0 - mu) / x1)
    _, jac = vcf.half_residuals_jac(x1, y2, t_full, mu)
    fd = np.empty((2, 3))
    for j, h in enumerate((1e-7 * x1, 1e-7 * y2, 1e-7 * t_full)):
        zp = [x1, y2, t_full]
        zm = [x1, y2, t_full]
        zp[j] += h
        zm[j] -= h
        rp = vcf.half_residuals(zp[0], zp[1], zp[2], mu)
        rm = vcf.half_residuals(zm[0], zm[1], zm[2], mu)
        fd[:, j] = [(rp[0] - rm[0]) / (2 * h), (rp[1] - rm[1]) / (2 * h)]
    assert np.allclose(jac, fd, rtol=1e-5, atol=1e-6)


def test_member_geometry_in_km_uses_system_length() -> None:
    """A near-circular P1 orbit of radius 0.5 (far from P2) at Ganymede's mu: min and max r1 about
    0.5 l* in km."""
    system = cr3bp.cr3bp_system("Jupiter", "Ganymede")
    x1 = 0.5
    y2 = math.sqrt((1.0 - system.mu) / x1)
    g = vcf.member_geometry(x1, y2, 2.0 * math.pi, system, n=2001)
    assert g.min_r1_km == pytest.approx(0.5 * system.l_km, rel=2e-3)
    assert g.max_r1_km == pytest.approx(0.5 * system.l_km, rel=2e-3)


def test_loop_tracer_closes_on_a_small_oval() -> None:
    """Just below mu''_2 the closed family i_2 is a small oval around the elliptic orbit: the
    extremum point is corrected there and the fixed-mu trace returns to it."""
    r = ROWS["mu2pp"]
    mu = r["mu"] * (1.0 - 1e-3)
    w = vcf.correct_extremum_at_mu(
        np.array([r["x1"] * 0.995, r["y2"] * 1.005, r["T"], math.log(mu)]), mu
    )
    assert w is not None
    members, why, closure = vcf.trace_loop(
        w[0], w[1], w[2], mu, direction=np.array([1.0, 0.0, 0.0]), ds=1e-3, ds_max=2e-3
    )
    assert why == "closed"
    assert closure < 1e-6
    assert len(members) >= 20


#: Committed members 64 and 65 of data/975_varin/loops/p2_up_loop.jsonl (i_1 at the Earth-Moon mu of
#: the P2 curve end), p = (x1(0), y2(0), T/2pi); member 65 sits at the oval's hairpin (#975 note
#: sec. 11.2).
_MU_EM_P2 = 0.01215058439469525
_UP_64 = [0.33323878351423125, 2.1668537174022484, 2.035542608069052]
_UP_65 = [0.32822603452771587, 2.1856692375299507, 2.0308441981041976]


def test_a2_guard_rejects_the_p2_up_corrector_jump() -> None:
    """Tool amendment A2 (#975 note sec. 11.4): from member 65 the unguarded 0.02 step lands on
    another family more than 0.5 away (the jump that broke the lead's P2 up trace); the guarded
    step rejects it, and a smaller guarded step lands within 0.5 |step| of its predictor."""
    p = np.array(_UP_65)
    t = vcf._tangent(vcf._loop_eval(p, _MU_EM_P2)[1], p - np.array(_UP_64))
    jumped = vcf._loop_step(p, t, 0.02, _MU_EM_P2, 1e-10)
    assert jumped is not None
    assert float(np.linalg.norm(jumped - p)) > 0.5
    stats: dict[str, int] = {}
    assert vcf._accepted_step(p, t, 0.02, _MU_EM_P2, 1e-10, True, stats) is None
    assert stats == {"displacement": 1}
    acc = vcf._accepted_step(p, t, 0.0025, _MU_EM_P2, 1e-10, True, None)
    assert acc is not None
    pc, tn = acc
    assert float(np.linalg.norm(pc - (p + 0.0025 * t))) <= vcf.MAX_DISP_FRAC * 0.0025
    assert float(tn @ t) >= math.cos(math.radians(vcf.MAX_TURN_DEG))
