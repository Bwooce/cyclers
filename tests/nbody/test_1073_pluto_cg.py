"""#1073: Pluto continuous-gravity lane: field, STM and chain tests at production scale (km).

Pre-registration: docs/notes/2026-10-12-1073-pluto-continuous-gravity-lane-prereg.md sec. 8 item 1.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.nbody.pluto_cg import (
    BinaryCG,
    Chain,
    NodeSpec,
    accel_and_gradient,
    flyby_seed_offsets,
    jacobi,
    propagate,
    wrap_q6,
)
from cyclerfinder.search.two_working_body import CircularSystem
from cyclerfinder.search.two_working_body_binary import (
    BarycentricCircularSystem,
    barycentric_binary,
)

DAY = 86400.0
GM_TOT = 975.5
GM_C = 106.1
P_C = 6.38722 * DAY
P_NIX = 360.0 / 14.484 * DAY


def _kepler_a(per: float) -> float:
    return float((GM_TOT * (per / (2.0 * math.pi)) ** 2) ** (1.0 / 3.0))


def _model() -> BinaryCG:
    return BinaryCG(
        gm_total=GM_TOT,
        gm_secondary=GM_C,
        a_rel_km=_kepler_a(P_C),
        period_s=P_C,
        secondary_radius_km=606.0,
        primary_radius_km=1188.3,
        targets={"Nix": (_kepler_a(P_NIX), P_NIX, 0.0)},
    )


def _bary() -> BarycentricCircularSystem:
    off = CircularSystem(
        GM_TOT,
        {"Charon": (_kepler_a(P_C), P_C, 0.0), "Nix": (_kepler_a(P_NIX), P_NIX, 0.0)},
        frozenset({"Nix"}),
    )
    return barycentric_binary(off, "Charon", GM_C / GM_TOT)


def test_rails_equal_barycentric_binary() -> None:
    m, s = _model(), _bary()
    for t in [0.0, 1.3 * DAY, 17.9 * DAY]:
        rc, vc, _ = m.secondary_state(t)
        rp, vp, _ = m.primary_state(t)
        rn, vn, _ = m.target_state("Nix", t)
        assert np.allclose(rc, s.state("Charon", t)[0], atol=1e-9, rtol=0)
        assert np.allclose(vc, s.state("Charon", t)[1], atol=1e-13, rtol=0)
        assert np.allclose(rp, s.primary_state(t)[0], atol=1e-9, rtol=0)
        assert np.allclose(vp, s.primary_state(t)[1], atol=1e-13, rtol=0)
        assert np.allclose(rn, s.state("Nix", t)[0], atol=1e-9, rtol=0)
        assert np.allclose(vn, s.state("Nix", t)[1], atol=1e-13, rtol=0)
    # the binary is the exact two-body solution: n^2 a_rel^3 = GM_total
    n = 2.0 * math.pi / m.period_s
    assert abs(n**2 * m.a_rel_km**3 / GM_TOT - 1.0) < 1e-12


def test_sigma0_field_is_the_monopole() -> None:
    m = _model()
    for r in [np.array([306.0, 40.0, 10.0]), np.array([-30000.0, 12000.0, 3000.0])]:
        a, _ = accel_and_gradient(m, 0.0, r, 2.1 * DAY)
        mono = -GM_TOT * r / np.linalg.norm(r) ** 3
        assert np.allclose(a, mono, rtol=1e-14, atol=0)


def test_sigma1_field_is_the_cr3bp() -> None:
    m = _model()
    t = 3.3 * DAY
    r = np.array([5000.0, -9000.0, 2500.0])
    a, _ = accel_and_gradient(m, 1.0, r, t)
    rc = m.secondary_state(t)[0]
    rp = m.primary_state(t)[0]
    exp = (
        -GM_C * (r - rc) / np.linalg.norm(r - rc) ** 3
        - (GM_TOT - GM_C) * (r - rp) / np.linalg.norm(r - rp) ** 3
    )
    assert np.allclose(a, exp, rtol=1e-13, atol=0)


@pytest.mark.parametrize("sigma", [0.0, 0.5, 1.0])
def test_gradient_matches_finite_difference(sigma: float) -> None:
    m = _model()
    t = 1.7 * DAY
    r = np.array([12000.0, 9000.0, -1500.0])
    _, g = accel_and_gradient(m, sigma, r, t)
    h = 1e-2
    fd = np.zeros((3, 3))
    for k in range(3):
        e = np.zeros(3)
        e[k] = h
        fd[:, k] = (
            accel_and_gradient(m, sigma, r + e, t)[0] - accel_and_gradient(m, sigma, r - e, t)[0]
        ) / (2 * h)
    assert np.allclose(g, fd, rtol=1e-6, atol=1e-16)


def _near_charon_state(m: BinaryCG, t: float) -> np.ndarray:
    rc, vc, _ = m.secondary_state(t)
    return np.concatenate(
        [rc + np.array([2500.0, 1500.0, 800.0]), vc + np.array([0.05, -0.2, 0.06])]
    )


def test_stm_matches_finite_difference_near_charon() -> None:
    m = _model()
    t0, t1 = 0.4 * DAY, 1.4 * DAY
    x0 = _near_charon_state(m, t0)
    x1, phi = propagate(m, 1.0, x0, t0, t1)
    rng = np.random.default_rng(1073)
    for _ in range(3):
        d = rng.normal(size=6)
        d[:3] *= 1.0  # km
        d[3:] *= 1e-5  # km/s
        eps = 1e-2
        xp, _ = propagate(m, 1.0, x0 + eps * d, t0, t1, stm=False)
        xm, _ = propagate(m, 1.0, x0 - eps * d, t0, t1, stm=False)
        fd = (xp - xm) / 2.0
        lin = phi @ (eps * d)
        assert np.linalg.norm(lin - fd) / np.linalg.norm(fd) < 1e-5
    assert np.all(np.isfinite(x1))


@pytest.mark.parametrize("sigma", [0.0, 0.5, 1.0])
def test_jacobi_conserved_3d(sigma: float) -> None:
    m = _model()
    t0 = 0.2 * DAY
    x0 = _near_charon_state(m, t0)
    j0 = jacobi(m, sigma, x0, t0)
    for t1 in [0.7 * DAY, 3.0 * DAY]:
        x1, _ = propagate(m, sigma, x0, t0, t1, stm=False)
        assert abs(jacobi(m, sigma, x1, t1) - j0) / abs(j0) < 1e-10


def test_wrap_preserves_jacobi() -> None:
    m = _model()
    t, dt = 1.1 * DAY, 3.7 * DAY
    x = _near_charon_state(m, t)
    q6 = wrap_q6(m, dt)
    for s in (0.3, 1.0):
        assert abs(jacobi(m, s, q6 @ x, t + dt) - jacobi(m, s, x, t)) < 1e-12 * abs(
            jacobi(m, s, x, t)
        )


def test_flyby_seed_offsets_give_the_patched_turn() -> None:
    m = _model()
    t = 2.0 * DAY
    vin = np.array([0.2, 0.05, 0.1])
    ang = math.radians(35.0)
    c, s = math.cos(ang), math.sin(ang)
    vout = np.array([c * vin[0] - s * vin[1], s * vin[0] + c * vin[1], vin[2]])
    vout *= np.linalg.norm(vin) / np.linalg.norm(vout)
    dr, dv = flyby_seed_offsets(GM_C, vin, vout)
    # a two-body hyperbola with these periapsis offsets turns by 35 deg at |vin|
    rp = float(np.linalg.norm(dr))
    vp = float(np.linalg.norm(dv))
    assert abs(float(dr @ dv)) < 1e-9 * rp * vp
    vinf2 = vp**2 - 2.0 * GM_C / rp
    e = 1.0 + rp * vinf2 / GM_C
    assert abs(math.sqrt(vinf2) - np.linalg.norm(vin)) < 1e-12
    turn = math.degrees(
        math.acos(float(vin @ vout) / float(np.linalg.norm(vin) * np.linalg.norm(vout)))
    )
    assert turn > 20.0
    assert abs(math.degrees(2.0 * math.asin(1.0 / e)) - turn) < 1e-9
    assert np.all(np.isfinite(m.secondary_state(t)[0]))


def _small_chain(m: BinaryCG, free_bz: bool) -> tuple[Chain, np.ndarray]:
    """Two Charon flyby nodes and one Nix pin; not an orbit, a Jacobian test bed."""
    specs = [NodeSpec("flyby"), NodeSpec("pin", "Nix"), NodeSpec("flyby")]
    ch = Chain(m, specs, period_s=3.0 * P_C, sigma=0.1, free_bz=free_bz)
    z: list[float] = []
    z += [400.0, 250.0, 60.0, 0.03, -0.25, 0.02, 1.0 * DAY]
    z += [0.08, 0.02, 0.01, 6.0 * DAY]
    if free_bz:
        z += [35.0]
    z += [-300.0, 410.0, -50.0, -0.2, 0.1, 0.03, 13.0 * DAY]
    return ch, np.asarray(z, dtype=np.float64)


@pytest.mark.parametrize("free_bz", [False, True])
def test_chain_jacobian_matches_finite_difference(free_bz: bool) -> None:
    m = _model()
    ch, z = _small_chain(m, free_bz)
    r, j, _ = ch.residual_and_jac(z)
    assert j.shape == (len(r), len(z))
    assert len(r) == 6 * 3 + 2  # 3 legs, 2 gauges
    steps = ch.fd_steps(z)
    for k in range(len(z)):
        zp, zm = z.copy(), z.copy()
        zp[k] += steps[k]
        zm[k] -= steps[k]
        fd = (
            ch.residual_and_jac(zp, want_jac=False)[0] - ch.residual_and_jac(zm, want_jac=False)[0]
        ) / (2 * steps[k])
        scale = max(float(np.linalg.norm(fd)), float(np.linalg.norm(j[:, k])), 1e-12)
        assert float(np.linalg.norm(j[:, k] - fd)) / scale < 2e-4, k


def test_chain_counts() -> None:
    m = _model()
    ch, z = _small_chain(m, False)
    assert ch.n_unknowns == 7 + 4 + 7 == len(z)
    ch2, z2 = _small_chain(m, True)
    assert ch2.n_unknowns == 7 + 5 + 7 == len(z2)


def test_free_period_chain_jacobian() -> None:
    """PC1 formulation: two free nodes at fixed epochs, period free, node-0 components fixed."""
    m = _model()
    specs = [NodeSpec("free", t_fixed=0.0), NodeSpec("free", t_fixed=1.5 * DAY)]
    ch = Chain(
        m, specs, period_s=3.0 * DAY, sigma=1.0, free_period=True, extra=((0, -8000.0), (1, 0.0))
    )
    z = np.array(
        [-8100.0, 30.0, 0.0, 0.01, -0.30, 0.0, 9000.0, -200.0, 0.0, -0.02, 0.25, 0.0, 3.1 * DAY]
    )
    r, j, _ = ch.residual_and_jac(z)
    assert ch.n_unknowns == 13 and len(r) == 12 + 2 == ch.n_rows
    steps = ch.fd_steps(z)
    for k in range(len(z)):
        zp, zm = z.copy(), z.copy()
        zp[k] += steps[k]
        zm[k] -= steps[k]
        fd = (
            ch.residual_and_jac(zp, want_jac=False)[0] - ch.residual_and_jac(zm, want_jac=False)[0]
        ) / (2 * steps[k])
        scale = max(float(np.linalg.norm(fd)), float(np.linalg.norm(j[:, k])), 1e-12)
        assert float(np.linalg.norm(j[:, k] - fd)) / scale < 2e-4, k


def _ladder_point(
    sigma: float, shift: float, rp_sigma: float, t_shift_days: float
) -> dict[str, object]:
    return {
        "sigma": sigma,
        "nodes": [
            {
                "kind": "flyby",
                "vinf_kms": 0.25 + shift * sigma,
                "t_days": 10.0 + t_shift_days * sigma,
                "rp_over_sigma": rp_sigma,
                "turn_deg": 40.0,
            },
            {"kind": "pin", "vinf_kms": 0.1 + 0.5 * shift * sigma, "t_days": 20.0},
        ],
    }


REF = [
    {"kind": "flyby", "vinf_kms": 0.25, "t_days": 10.0, "rp_km": 4000.0, "turn_deg": 40.0},
    {"kind": "pin", "vinf_kms": 0.1, "t_days": 20.0},
]


def test_identity_pass_sibling_fail() -> None:
    from cyclerfinder.nbody.pluto_cg import identity_check

    lad = [_ladder_point(s, 0.2, 4000.0, 1.0) for s in (0.005, 0.0025, 0.00125)]
    assert identity_check(lad, REF)["verdict"] == "PASS"
    # a different periapsis in the limit: sibling
    lad = [_ladder_point(s, 0.2, 5000.0, 1.0) for s in (0.0025, 0.00125)]
    assert identity_check(lad, REF)["verdict"] == "SIBLING"
    # an epoch 2.46 d away (the phb 26 sibling zero, #1065 sec. 10.2): fail
    lad = [_ladder_point(s, 0.2, 4000.0, 1.0) for s in (0.0025, 0.00125)]
    far = [dict(r) for r in REF]
    far[0] = dict(far[0], t_days=12.46)
    assert identity_check(lad, far)["verdict"] == "FAIL"
    # a ladder that stops too high
    lad = [_ladder_point(s, 0.2, 4000.0, 1.0) for s in (0.01, 0.005)]
    assert identity_check(lad, REF)["verdict"] == "FAIL"
