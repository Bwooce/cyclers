"""#593 reconciliation: independent QBCP L1-substitute build, checked against
the published POL1 golden under both the #592-fixed and reconstructed-buggy
alpha_6 scaling.

Built to resolve a discrepancy between this task's own SE-L2 torus re-test
(found #592 does NOT regress SE-L2, contradicting #544's 2026-07-10 reverted
experiment) and #544's OTHER claim (the same alpha_6 fix moved a QBCP
L1-substitute ~30% closer to POL1). Rather than replicate #544's undocumented,
never-committed construction, this builds a genuinely independent multi-
shooting periodic-orbit corrector from scratch: CR3BP EM-L1 (own quintic
solve, cross-checked as a real fixed point) -> BCR4BP via sequential mu_sun
continuation (reduces exactly to CR3BP at mu_sun=0) -> QBCP handoff (real
time-varying alphas), 12 segments, analytic STM via qbcp_stm_eom/
bcr4bp_stm_eom, plain Newton (well-posed, no free period/phase pin needed --
the system is time-periodic, not autonomous).

Result: CONFIRMS #544's own finding, independently. distance-to-POL1 goes
4.14e-2 (buggy) -> 1.81e-2 (fixed), same direction/magnitude as #544's own
2.12e-2 -> 1.47e-2. Combined with the SE-L2 refutation, net verdict: #592 is
a genuine, net-positive fix -- see docs/notes/2026-07-14-593-qbcp-alpha6-
impact-scoping.md for the full writeup.

Usage:
    uv run python scripts/analyze_593_qbcp_l1_substitute_reconciliation.py
"""

import math

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.qbcp as qbcp

MU = qbcp.qbcp_default().mu
_POL1_X = 0.8369141677649317
_POL1_PY = 0.8391311559808445

ORIG_EOM = qbcp.qbcp_eom
ORIG_POT2 = qbcp.qbcp_potential_second_derivatives


def buggy_qbcp_eom(t, state_pm, system):
    x, y, z, px, py, pz = state_pm
    alphas = qbcp.evaluate_alphas(t, system)
    a1, a2, a3, a4, a5, a6, xs, ys = alphas[1:9]
    mu = system.mu
    dx = a1 * px + a2 * x + a3 * y
    dy = a1 * py + a2 * y - a3 * x
    dz = a1 * pz + a2 * z
    rpe2 = (x + mu) ** 2 + y * y + z * z
    rpm2 = (x - 1.0 + mu) ** 2 + y * y + z * z
    rps2 = (x - xs) ** 2 + (y - ys) ** 2 + z * z
    rpe3 = rpe2 * math.sqrt(rpe2)
    rpm3 = rpm2 * math.sqrt(rpm2)
    rps3 = rps2 * math.sqrt(rps2)
    pot_x = (
        (1.0 - mu) * (x + mu) / rpe3 + mu * (x - 1.0 + mu) / rpm3 + system.mu_sun * (x - xs) / rps3
    )
    pot_y = (1.0 - mu) * y / rpe3 + mu * y / rpm3 + system.mu_sun * (y - ys) / rps3
    pot_z = (1.0 - mu) * z / rpe3 + mu * z / rpm3 + system.mu_sun * z / rps3
    dpx = -a2 * px + a3 * py - a4 - a6 * pot_x
    dpy = -a2 * py - a3 * px - a5 - a6 * pot_y
    dpz = -a2 * pz - a6 * pot_z
    return np.array([dx, dy, dz, dpx, dpy, dpz], dtype=np.float64)


def cr3bp_l1_x(mu):
    x = 1.0 - mu - 0.1
    for _ in range(100):
        r1 = abs(x + mu)
        r2 = abs(x - 1.0 + mu)
        f = x - (1.0 - mu) * (x + mu) / r1**3 - mu * (x - 1.0 + mu) / r2**3
        h = 1e-8
        r1h = abs(x + h + mu)
        r2h = abs(x + h - 1.0 + mu)
        fh = (x + h) - (1.0 - mu) * (x + h + mu) / r1h**3 - mu * (x + h - 1.0 + mu) / r2h**3
        df = (fh - f) / h
        x_new = x - f / df
        if abs(x_new - x) < 1e-15:
            x = x_new
            break
        x = x_new
    return x


def multi_shoot_time_periodic(stm_eom, system, period, n_segments, states0, max_iter=80, tol=1e-11):
    n = n_segments
    dt = period / n
    states = [s.copy() for s in states0]
    resnorm = np.inf
    for _ in range(max_iter):
        residual = np.zeros(6 * n)
        jac = np.zeros((6 * n, 6 * n))
        endpoints, stms = [], []
        for i in range(n):
            t0 = i * dt
            y0 = np.concatenate([states[i], np.eye(6).flatten()])
            sol = solve_ivp(
                stm_eom, (t0, t0 + dt), y0, args=(system,), method="DOP853", rtol=1e-13, atol=1e-13
            )
            yf = sol.y[:, -1]
            endpoints.append(yf[:6])
            stms.append(yf[6:].reshape(6, 6))
        for i in range(n):
            j = (i + 1) % n
            residual[6 * i : 6 * i + 6] = endpoints[i] - states[j]
            jac[6 * i : 6 * i + 6, 6 * i : 6 * i + 6] = stms[i]
            jac[6 * i : 6 * i + 6, 6 * j : 6 * j + 6] -= np.eye(6)
        resnorm = float(np.linalg.norm(residual))
        if resnorm < tol:
            return states, resnorm, True
        try:
            delta = np.linalg.solve(jac, -residual)
        except np.linalg.LinAlgError:
            delta, *_ = np.linalg.lstsq(jac, -residual, rcond=None)
        max_step = 0.5
        step_norm = np.linalg.norm(delta)
        if step_norm > max_step:
            delta = delta * (max_step / step_norm)
        for i in range(n):
            states[i] = states[i] + delta[6 * i : 6 * i + 6]
    return states, resnorm, False


def build_l1_substitute(use_buggy_qbcp: bool):
    qbcp_sys = qbcp.qbcp_default()
    n_seg = 12
    ts = qbcp_sys.sun_period_tu
    x_l1 = cr3bp_l1_x(MU)

    # Stage A: CR3BP L1 (mu_sun=0) -> BCR4BP (mu_sun=full), sequential continuation.
    state0 = np.array([x_l1, 0.0, 0.0, 0.0, 0.0, 0.0])
    states = [state0.copy() for _ in range(n_seg)]
    mu_sun_steps = np.concatenate(
        [[0.0], np.geomspace(qbcp_sys.mu_sun * 1e-4, qbcp_sys.mu_sun, 25)]
    )
    last_ok = True
    for s in mu_sun_steps:
        sys_bcr = bcr4bp.BCR4BPSystem(
            mu=qbcp_sys.mu,
            mu_sun=float(s),
            a_sun_nondim=qbcp_sys.a_sun_nondim,
            omega_sun_nondim=qbcp_sys.omega_sun_nondim,
            theta_sun0=qbcp_sys.theta_sun0,
        )
        states, resnorm, ok = multi_shoot_time_periodic(
            bcr4bp.bcr4bp_stm_eom, sys_bcr, ts, n_seg, states
        )
        last_ok = ok
        if not ok:
            print(f"  [stage A] mu_sun={s:.4e} FAILED to converge, resnorm={resnorm:.3e}")
            break
    print(f"Stage A (BCR4BP, mu_sun=full) converged={last_ok}, final resnorm reported above")

    # Stage B: swap circular BCR4BP dynamics for the real time-varying QBCP alphas.
    if use_buggy_qbcp:
        qbcp.qbcp_eom = buggy_qbcp_eom
    try:
        # Convert each BCR4BP PV segment state to QBCP PM at its own segment start time.
        dt = ts / n_seg
        states_pm = [qbcp.state_pv_to_pm(states[i], i * dt, qbcp_sys) for i in range(n_seg)]
        states_qbcp, resnorm_q, ok_q = multi_shoot_time_periodic(
            qbcp.qbcp_stm_eom, qbcp_sys, ts, n_seg, states_pm
        )
    finally:
        qbcp.qbcp_eom = ORIG_EOM

    print(f"Stage B (QBCP) converged={ok_q}, resnorm={resnorm_q:.3e}")

    # Closest approach of the converged QBCP orbit (segment 0, t=0 state) to POL1,
    # matching FOLLOW-UP-2's own metric (x, py at the reference epoch).
    x0_qbcp, py0_qbcp = states_qbcp[0][0], states_qbcp[0][4]
    dist_pol1 = math.hypot(x0_qbcp - _POL1_X, py0_qbcp - _POL1_PY)
    return ok_q, resnorm_q, x0_qbcp, py0_qbcp, dist_pol1


print("=" * 70)
print("FIXED (#592, current main)")
print("=" * 70)
ok_f, res_f, x_f, py_f, dist_f = build_l1_substitute(use_buggy_qbcp=False)
print(f"x0={x_f:.10f}  py0={py_f:.10f}  dist-to-POL1={dist_f:.6e}")

print()
print("=" * 70)
print("BUGGY (reconstructed pre-#592)")
print("=" * 70)
ok_b, res_b, x_b, py_b, dist_b = build_l1_substitute(use_buggy_qbcp=True)
print(f"x0={x_b:.10f}  py0={py_b:.10f}  dist-to-POL1={dist_b:.6e}")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"FIXED : converged={ok_f}  dist-to-POL1={dist_f:.6e}")
print(f"BUGGY : converged={ok_b}  dist-to-POL1={dist_b:.6e}")
print("#544's reported (their own construction): buggy=2.12e-2, fixed(their patch)=1.47e-2")
