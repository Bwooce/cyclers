"""#970: re-close Schwaniger 1963's retrograde Earth-Moon periodic free return (NASA TN D-1833
sec. III.E) with the project's CR3BP corrector, and record it as a positive control.

Source: docs/notes/2026-10-07-digest-schwaniger-1963-earth-moon-symmetrical-free-return.md.
Schwaniger defines the orbit as the counter-rotation cislunar symmetric free return whose
perigee (r = 6555 km, 100 n.mi. altitude, horizontal) lies on the Earth-Moon line. In the
one-parameter family of x-axis-symmetric periodic orbits (parameter C) that is the member whose
half-period perpendicular crossing near the Earth is at r1 = 6555 km, so the script:

1. corrects the symmetric orbit at fixed C with ``correct_symmetric_fixed_jacobi`` (plain
   DOP853 in physical time; no regularisation) from the periselenum start ``(x0, 0, 0, ydot0)``;
2. solves for C by a secant on ``r1(t_half) - 6555 km``;
3. integrates the full-period variational equations (``stm_mode="fixed_path"``) for the
   monodromy and the stability indices, and checks them against the half-period Barden form;
4. re-propagates the corrected state over one period with two other integrators: the
   Moon-centred KS propagator (``core.cr3bp_ks``, with its own transition matrix) and the
   Sundman r1*r2 time-transformed propagator (``core.cr3bp_regularized``).

Run for the digest's mu = 0.012150 (L = 384,400 km) and for the registry Earth-Moon system.
Output: data/970_schwaniger/control.json (and stdout).
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

from cyclerfinder.core.constants import PLANETS
from cyclerfinder.core.cr3bp import (
    CR3BPSystem,
    cr3bp_system,
    jacobi_constant,
    propagate,
)
from cyclerfinder.core.cr3bp_ks import MoonCentredCR3BP, propagate_ks
from cyclerfinder.core.cr3bp_regularized import (
    physical_to_regularized_span,
    propagate_regularized,
)
from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.data.method_capability import MethodCapability
from cyclerfinder.data.preflight import preflight_search
from cyclerfinder.search.cr3bp_periodic import (
    _xaxis_crossings,
    barden_stability,
    correct_symmetric_fixed_jacobi,
)

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "data" / "970_schwaniger" / "control.json"

R_PERIGEE_KM = 6555.0  # Schwaniger p.2: injection and re-entry at 100 n.mi. altitude
# Digest seed (mu = 0.012150): periselenum start, 6 digits as printed in the digest.
SEED_X0 = 0.982120
SEED_YDOT0 = -2.471279
SEED_PERIOD = 6.0015
# Schwaniger's printed values (p.6-7): "about 2150 km", "about 650 hours".
PAPER_PERISELENUM_KM = 2150.0
PAPER_PERIOD_H = 650.0


def _ts(*a: object) -> None:
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def _member(system: CR3BPSystem, x0: float, jac: float, period: float) -> dict:
    orb = correct_symmetric_fixed_jacobi(
        system, x0, jac, period, ydot0_sign=-1.0, half_crossings=2, tol=1e-12
    )
    if not orb.converged:
        raise RuntimeError(f"no convergence at C={jac}: {orb}")
    s0 = np.array([orb.x0, 0.0, 0.0, 0.0, orb.ydot0, 0.0])
    half = propagate(system, s0, orb.t_half)
    xh = half.state_f
    r1 = math.hypot(xh[0] + system.mu, xh[1]) * system.l_km
    return {"orb": orb, "s0": s0, "r1_half_km": r1, "state_half": xh}


def solve_member(system: CR3BPSystem) -> dict:
    """Secant on C so that the half-period crossing is at r1 = 6555 km."""
    c_a = jacobi_constant(np.array([SEED_X0, 0, 0, 0, SEED_YDOT0, 0.0]), system.mu)
    m_a = _member(system, SEED_X0, c_a, SEED_PERIOD)
    c_b = c_a + 2e-4
    m_b = _member(system, m_a["orb"].x0, c_b, m_a["orb"].period)
    f_a, f_b = m_a["r1_half_km"] - R_PERIGEE_KM, m_b["r1_half_km"] - R_PERIGEE_KM
    for it in range(40):
        c_n = c_b - f_b * (c_b - c_a) / (f_b - f_a)
        m_n = _member(system, m_b["orb"].x0, c_n, m_b["orb"].period)
        f_n = m_n["r1_half_km"] - R_PERIGEE_KM
        _ts(f"  secant {it}: C={c_n:.12f} r1(T/2)-6555 = {f_n:+.3e} km")
        c_a, f_a, c_b, f_b, m_b = c_b, f_b, c_n, f_n, m_n
        if abs(f_n) < 1e-6:
            break
    else:
        raise RuntimeError("secant on C did not converge")
    m_b["C"] = c_b
    return m_b


def _min_distances(system: CR3BPSystem, s0: np.ndarray, period: float) -> dict:
    """Dense sampling plus event-refined minima of r1, r2 and max r1 over one period."""
    from scipy.integrate import solve_ivp

    from cyclerfinder.core.cr3bp import cr3bp_eom

    mu = system.mu

    def dr1(t: float, y: np.ndarray, m: float) -> float:
        return float((y[0] + m) * y[3] + y[1] * y[4])

    def dr2(t: float, y: np.ndarray, m: float) -> float:
        return float((y[0] - 1 + m) * y[3] + y[1] * y[4])

    sol = solve_ivp(
        cr3bp_eom,
        (0, period),
        s0,
        args=(mu,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        events=[dr1, dr2],
    )
    r1s = [math.hypot(y[0] + mu, y[1]) for y in sol.y_events[0]]
    r2s = [math.hypot(y[0] - 1 + mu, y[1]) for y in sol.y_events[1]]
    r2s.append(math.hypot(s0[0] - 1 + mu, s0[1]))
    return {
        "r1_min_km": min(r1s) * system.l_km,
        "r1_max_km": max(r1s) * system.l_km,
        "r2_min_km": min(r2s) * system.l_km,
        "n_r1_extrema": len(r1s),
        "n_r2_extrema": len(r2s),
    }


def monodromy_report(system: CR3BPSystem, s0: np.ndarray, period: float, orb: object) -> dict:
    arc = propagate(system, s0, period, with_stm=True, stm_mode="fixed_path")
    m = arc.stm
    assert m is not None
    idx = [0, 1, 3, 4]
    m4 = m[np.ix_(idx, idx)]
    mz = m[np.ix_([2, 5], [2, 5])]
    eig4 = np.linalg.eigvals(m4)
    eigz = np.linalg.eigvals(mz)
    b_h = float(np.trace(m4) - 2.0)
    b_v = float(np.trace(mz))
    nu_barden, lam_barden = barden_stability(system, orb)  # type: ignore[arg-type]
    order = np.argsort(np.abs(eig4 - 1.0))
    closure = arc.state_f - s0
    return {
        "closure_full_period_plain": {
            "dr_km": float(np.linalg.norm(closure[:3]) * system.l_km),
            "dv_m_s": float(np.linalg.norm(closure[3:]) * system.l_km / system.t_s * 1e3),
        },
        "det_M6": float(np.linalg.det(m)),
        "det_M4": float(np.linalg.det(m4)),
        "eig_planar": [complex(e) for e in eig4],
        "eig_vertical": [complex(e) for e in eigz],
        "trivial_pair": [complex(eig4[order[0]]), complex(eig4[order[1]])],
        "b_h_trace": b_h,
        "b_v_trace": b_v,
        "lambda_max": float(np.max(np.abs(eig4))),
        "barden_nu": nu_barden,
        "barden_lambda": complex(lam_barden),
        "barden_b_h_equiv": 2.0 * nu_barden,
        "stm_planar": m4.tolist(),
    }


def ks_crosscheck(system: CR3BPSystem, s0: np.ndarray, period: float, m_plain: dict) -> dict:
    model = MoonCentredCR3BP(system.mu)
    arc = propagate_ks(model, s0, period, with_stm=True)
    assert arc.stm is not None
    idx = [0, 1, 3, 4]
    m4 = arc.stm[np.ix_(idx, idx)]
    eig = np.linalg.eigvals(m4)
    d = arc.state - s0
    return {
        "dr_km": float(np.linalg.norm(d[:3]) * system.l_km),
        "dv_m_s": float(np.linalg.norm(d[3:]) * system.l_km / system.t_s * 1e3),
        "b_h_trace": float(np.trace(m4) - 2.0),
        "lambda_max": float(np.max(np.abs(eig))),
        "hamiltonian_max": arc.hamiltonian_max,
        "bilinear_max": arc.bilinear_max,
        "energy_drift": arc.energy_drift,
        "r_min_moon_km": arc.r_min * system.l_km,
        "nfev": arc.nfev,
        "rel_b_h_vs_plain": float(
            abs(np.trace(m4) - 2.0 - m_plain["b_h_trace"]) / abs(m_plain["b_h_trace"])
        ),
    }


def sundman_crosscheck(system: CR3BPSystem, s0: np.ndarray, period: float) -> dict:
    span = physical_to_regularized_span(system, s0, (0.0, period))
    span = (span[0], span[1] * 1.5)
    arc = propagate_regularized(system, s0, span, regularization="r1r2", t_stop=period)
    sf = arc.state_at_s[:, -1]
    d = sf - s0
    return {
        "t_end": float(arc.t_at_s[-1]),
        "dr_km": float(np.linalg.norm(d[:3]) * system.l_km),
        "dv_m_s": float(np.linalg.norm(d[3:]) * system.l_km / system.t_s * 1e3),
        "nfev": arc.nfev,
    }


def run(system: CR3BPSystem, label: str) -> dict:
    mu = system.mu
    _ts(f"[{label}] mu={mu!r} L={system.l_km} km TU={system.t_s} s")
    mem = solve_member(system)
    orb = mem["orb"]
    s0 = mem["s0"]
    period = orb.period
    times, _states = _xaxis_crossings(
        system, s0, 1.01 * period, with_stm=False, rtol=1e-12, atol=1e-12
    )
    dist = _min_distances(system, s0, period)
    vu = system.l_km / system.t_s
    xh = mem["state_half"]
    v_rot = math.hypot(xh[3], xh[4]) * vu
    # Earth-relative inertial velocity at perigee: rotating velocity plus omega x (r - r_E)
    # (omega = 1; the Earth's own inertial velocity is omega x r_E).
    vx_in, vy_in = xh[3] - xh[1], xh[4] + xh[0] + mu
    v_in = math.hypot(vx_in, vy_in) * vu
    h_in = (xh[0] + mu) * vy_in - xh[1] * vx_in
    mono = monodromy_report(system, s0, period, orb)
    _ts(
        f"[{label}] monodromy: b_h={mono['b_h_trace']:.6g} "
        f"barden 2nu={mono['barden_b_h_equiv']:.6g}"
    )
    ks = ks_crosscheck(system, s0, period, mono)
    sund = sundman_crosscheck(system, s0, period)
    r_e = PLANETS["E"].radius_eq_km
    r_m = SATELLITES["Moon"].radius_eq_km
    c = mem["C"]
    res = {
        "label": label,
        "mu": mu,
        "l_km": system.l_km,
        "t_s": system.t_s,
        "x0": orb.x0,
        "ydot0": orb.ydot0,
        "state_nd": [orb.x0, 0.0, 0.0, 0.0, orb.ydot0, 0.0],
        "jacobi_project": c,
        "jacobi_with_mu1mu": c + mu * (1.0 - mu),
        "t_half": orb.t_half,
        "period_nd": period,
        "period_h": period * system.t_s / 3600.0,
        "period_d": period * system.t_s / 86400.0,
        "crossing_residual_xdot": orb.crossing_residual,
        "y_crossings_half_period": int(np.sum(times < orb.t_half + 1e-9)),
        "perigee_km": mem["r1_half_km"],
        "perigee_alt_km": mem["r1_half_km"] - r_e,
        "perigee_x_nd": float(xh[0]),
        "periselene_km": dist["r2_min_km"],
        "periselene_alt_km": dist["r2_min_km"] - r_m,
        "apogee_km": dist["r1_max_km"],
        "r1_min_km_over_period": dist["r1_min_km"],
        "perigee_speed_rot_m_s": v_rot * 1e3,
        "perigee_speed_inertial_m_s": v_in * 1e3,
        "perigee_inertial_ang_mom_sign": "retrograde" if h_in < 0 else "prograde",
        "earth_floor_alt_km": PLANETS["E"].safe_alt_km,
        "moon_floor_alt_km": SATELLITES["Moon"].safe_alt_km,
        "passes_earth_floor": bool(mem["r1_half_km"] - r_e >= PLANETS["E"].safe_alt_km),
        "passes_moon_floor": bool(dist["r2_min_km"] - r_m >= SATELLITES["Moon"].safe_alt_km),
        "monodromy": mono,
        "ks_moon_centred_crosscheck": ks,
        "sundman_r1r2_crosscheck": sund,
        "paper": {"periselenum_km_about": PAPER_PERISELENUM_KM, "period_h_about": PAPER_PERIOD_H},
    }
    res["paper_rel_diff"] = {
        "periselenum": res["periselene_km"] / PAPER_PERISELENUM_KM - 1.0,
        "period": res["period_h"] / PAPER_PERIOD_H - 1.0,
    }
    _ts(
        f"[{label}] x0={orb.x0:.15f} ydot0={orb.ydot0:.15f} C={c:.12f} "
        f"(+mu(1-mu): {res['jacobi_with_mu1mu']:.12f}) T={period:.10f} "
        f"({res['period_h']:.2f} h) perigee {res['perigee_km']:.3f} km "
        f"periselene {res['periselene_km']:.2f} km apogee {res['apogee_km']:.0f} km"
    )
    _ts(
        f"[{label}] closure plain {mono['closure_full_period_plain']}, "
        f"KS dr={ks['dr_km']:.3e} km b_h rel {ks['rel_b_h_vs_plain']:.2e}, "
        f"Sundman dr={sund['dr_km']:.3e} km"
    )
    return res


def _jsonable(o: object) -> object:
    if isinstance(o, complex):
        return [o.real, o.imag]
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, list | tuple):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.floating | np.integer | np.bool_):
        return o.item()
    return o


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.parse_args()
    preflight_search(
        task_no=970,
        region_id="schwaniger-1963-em-retrograde-free-return-control",
        method=MethodCapability(
            genome="planar x-axis-symmetric periodic orbit, CR3BP Earth-Moon",
            corrector="correct_symmetric_fixed_jacobi + secant on C (perigee 6555 km)",
            capability_tags=frozenset({"ballistic", "cr3bp", "planar"}),
            git_sha="working-tree",
        ),
        script_path=Path(__file__),
        n_points=2,
    )
    digest_sys = CR3BPSystem(
        mu=0.012150,
        primary="Earth",
        secondary="Moon",
        l_km=384400.0,
        t_s=math.sqrt(384400.0**3 / (398600.4 / (1.0 - 0.012150))),
    )
    results = [
        run(digest_sys, "digest-mu-0.012150"),
        run(cr3bp_system("Earth", "Moon"), "registry"),
    ]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(_jsonable({"task": 970, "results": results}), indent=1) + "\n")
    _ts(f"wrote {OUT}")


if __name__ == "__main__":
    main()
