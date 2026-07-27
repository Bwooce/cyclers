"""#724: independent foreground re-run of #720/#723's mu_Io continuation.

Mirrors tests/search/test_variational_crnbp_torus.py::
test_continuation_reaches_physical_mu_io_at_n1_2 exactly (same seed pipeline,
same solver settings), printing every quantity the confirmation pass needs.
"""

import math

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.crnbp as crnbp
import cyclerfinder.search.variational_ccr4bp_torus as vt
import cyclerfinder.search.variational_crnbp_torus as vc
from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.genome.composed_moon_map import resonance_semimajor


def resonant_symmetric_orbit(
    mu: float,
    p_sc: int,
    q_moon: int,
    max_iter: int = 60,
    tol: float = 1e-12,
    cap: float = 0.05,
) -> tuple[np.ndarray, float, float]:
    a = resonance_semimajor(p_sc, q_moon)
    period = 2.0 * np.pi * q_moon
    th = 0.5 * period
    x0 = a - mu
    vy0 = float(np.sqrt((1.0 - mu) / a)) - x0
    res = np.inf
    for k in range(max_iter):
        s0 = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0])
        y42 = np.concatenate([s0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom,
            (0.0, th),
            y42,
            args=(mu,),
            method="DOP853",
            rtol=1e-12,
            atol=1e-12,
        )
        sf = sol.y[:, -1]
        phi = sf[6:].reshape(6, 6)
        g = np.array([sf[1], sf[3]])
        res = float(np.linalg.norm(g))
        if res < tol:
            break
        jac = np.array([[phi[1, 0], phi[1, 4]], [phi[3, 0], phi[3, 4]]])
        dz = np.linalg.solve(jac, -g)
        dz = dz * (0.3 if k < 8 else 1.0)
        norm = float(np.linalg.norm(dz))
        if norm > cap:
            dz = dz / norm * cap
        x0 += dz[0]
        vy0 += dz[1]
    return np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0]), period, res


def main() -> None:
    sys4 = ccr4bp.jupiter_europa_ganymede_default()
    target = crnbp.jupiter_europa_io_ganymede_default()
    io, gan = target.perturbers

    print("=== System parameters (current code, post-#723) ===")
    print(f"mu (Europa)      = {sys4.mu:.10e}")
    print(f"mu_io            = {io.mu:.10e}   a_io  = {io.a:.10f}")
    print(f"mu_gan           = {gan.mu:.10e}   a_gan = {gan.a:.10f}")
    print(f"omega_gan        = {gan.omega:.10f}  (registry two-body synodic rate)")
    print(f"omega_io         = {io.omega:.10f}  (= -2*omega_gan: {io.omega == -2.0 * gan.omega})")
    print(f"theta_io0        = {io.theta0:.10f}  (pi = {math.pi:.10f})")
    print(f"theta_gan0       = {gan.theta0:.10f}")
    print(
        f"Laplace phi_L    = theta_io0 + 2*theta_gan0 = "
        f"{(io.theta0 + 2 * gan.theta0) % (2 * math.pi):.10f} rad "
        f"({math.degrees((io.theta0 + 2 * gan.theta0) % (2 * math.pi)):.4f} deg)"
    )

    # Independent ephemeris cross-check of the synodic rates (IAU/JPL sidereal
    # periods, days): Io 1.769137786, Europa 3.551181041, Ganymede 7.15455296.
    p_io, p_eu, p_gan = 1.769137786, 3.551181041, 7.15455296
    omega_gan_eph = p_eu / p_gan - 1.0
    omega_io_eph = p_eu / p_io - 1.0
    print("\n=== Independent ephemeris cross-check (sidereal periods) ===")
    print(
        f"omega_gan (ephemeris) = {omega_gan_eph:.10f}  vs code {gan.omega:.10f}  "
        f"rel dev = {abs(gan.omega - omega_gan_eph) / abs(omega_gan_eph):.3e}"
    )
    print(
        f"omega_io  (ephemeris) = {omega_io_eph:.10f}  vs code {io.omega:.10f}  "
        f"rel dev = {abs(io.omega - omega_io_eph) / abs(omega_io_eph):.3e}"
    )
    print(
        f"TCP idealized omega_gan = -0.5 exactly; physical deviates by "
        f"{abs(omega_gan_eph + 0.5) / 0.5:.3e} relative"
    )
    print(f"Io/Gan synodic ratio (ephemeris) = {omega_io_eph / omega_gan_eph:.9f} (Laplace -2)")

    # Registry SMA sanity (seed orbit is interior: a_34 < 1).
    a_34 = resonance_semimajor(3, 4)
    print(f"\n3:4 resonant orbit SMA (Europa units) = {a_34:.10f}  (interior: {a_34 < 1.0})")
    print(f"a_io registry = {SATELLITES['Io'].sma_km / SATELLITES['Europa'].sma_km:.6f}")

    print("\n=== Stage 1: Kumar-class Europa 3:4 symmetric resonant CR3BP orbit ===")
    s0, period, res = resonant_symmetric_orbit(sys4.mu, 3, 4)
    print(
        f"x0 = {s0[0]:.12f}, vy0 = {s0[4]:.12f}, period = {period:.10f} "
        f"({period / (2 * np.pi):.1f} Europa periods), perp residual = {res:.3e}"
    )
    assert res < 1e-10

    print("\n=== Stage 2: #690 CCR4BP torus (n1=2, n2=20) ===")
    phys = vt.discover_ccr4bp_torus_from_resonant_orbit(
        sys4,
        s0,
        period,
        n1=2,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    print(
        f"residual_rms = {phys.residual_rms:.10e}, closure = {phys.closure_residual:.6e}, "
        f"rotation = {phys.rotation_number:.9f}, converged = {phys.converged}"
    )

    print("\n=== Stage 3: mu_io = 0 CRNBP seed (theta_io0 = pi, physical phase) ===")
    seed = vc.discover_crnbp_torus_from_ccr4bp_seed(
        phys,
        mu_io=0.0,
        a_io=io.a,
        omega_io=io.omega,
        theta_io0=io.theta0,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    print(
        f"residual_rms = {seed.residual_rms:.10e}, closure = {seed.closure_residual:.6e}, "
        f"rotation = {seed.rotation_number:.9f}"
    )

    print("\n=== Stage 4: 8-step mu_Io continuation to the physical mass ===")
    steps = vc.continue_crnbp_torus_mu_io(
        seed,
        target,
        n_steps=8,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    for i, st in enumerate(steps, 1):
        print(
            f"step {i}: mu_io = {st.system.perturbers[0].mu:.6e}  "
            f"residual_rms = {st.residual_rms:.10e}  closure = {st.closure_residual:.6e}  "
            f"rot = {st.rotation_number:.9f}"
        )
    final = steps[-1]
    print(
        f"\nFINAL: mu_io = {final.system.perturbers[0].mu:.10e} "
        f"(target {target.perturbers[0].mu:.10e})"
    )
    print(f"FINAL residual_rms = {final.residual_rms:.10e}")
    print(f"FINAL theta_io0 in system = {final.system.perturbers[0].theta0:.10f}")
    print(
        f"expected #723 value ~1.2343143649e-04; "
        f"match = {abs(final.residual_rms - 1.2343143649e-04) / 1.2343143649e-04:.3e} rel"
    )


if __name__ == "__main__":
    main()
