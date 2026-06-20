"""#412 Phase-B scoping spike: does a BCR4BP EM-libration orbit's support reach the
SE-L region (~3.9 lunar distances from Earth) as the Sun term turns on?

Feasibility gate for the full Phase-B build (a both-region-spanning synodic-resonant
cross-system cycle). For a few EM-L1 Lyapunov amplitudes: build the CR3BP orbit, cast it
as a BCR4BP-at-mu_sun=0 seed, mu_sun-continue to the full Sun value, and measure the max
Earth-relative reach (nondim = lunar distances) of the converged BCR4BP orbit. If reach
approaches ~3.9 LD (SE-L), the spanning object is reachable -> full build justified; if it
stays EM-Hill-bounded (~1-2 LD, the published outcome), that's a fast structural negative.

Reuse only (correct_symmetric_fixed_jacobi / correct_bcr4bp_periodic / continuation are
already tested). No orbit claimed beyond the reach measurement. Run:
  uv run python scripts/spike_412_bcr4bp_reach.py
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

import cyclerfinder.core.bcr4bp as bcr4bp
import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.genome.bcr4bp_continuation import continue_bcr4bp_family_in_musun
from cyclerfinder.genome.bcr4bp_genome import (
    IDX_T,
    IDX_X,
    IDX_XDOT,
    IDX_Y,
    IDX_YDOT,
    IDX_ZDOT,
    correct_bcr4bp_periodic,
)
from cyclerfinder.genome.bcr4bp_systems import SEM_ANDREU, build_bcr4bp_system
from cyclerfinder.search.cr3bp_periodic import correct_symmetric_fixed_jacobi

MU_EM = 0.012150581600000
SE_L_REACH_LD = 1.5e6 / 384400.0  # SE-L region distance from Earth, in lunar distances (~3.9)


def _max_earth_reach_ld(system: bcr4bp.BCR4BPSystem, state0: np.ndarray, period: float) -> float:
    """Max Earth-relative distance over one period, in lunar distances (nondim = LD)."""
    earth_x = -system.mu  # Earth at (-mu, 0) in EM-rotating nondim coords
    te = np.linspace(0.0, period, 1500)
    sol = solve_ivp(
        bcr4bp.bcr4bp_eom,
        (0.0, period),
        np.asarray(state0, float),
        args=(system,),
        method="DOP853",
        rtol=1e-11,
        atol=1e-11,
        t_eval=te,
    )
    if not sol.success:
        return float("nan")
    dx = sol.y[0, :] - earth_x
    dy = sol.y[1, :]
    return float(np.max(np.hypot(dx, dy)))


def main() -> None:
    cr = cr3bp.CR3BPSystem(mu=MU_EM, primary="earth", secondary="moon", l_km=384400.0, t_s=375190.0)
    sys_zero = build_bcr4bp_system(SEM_ANDREU, mu_sun_override=0.0)
    target = SEM_ANDREU.mu_sun
    free_vars = (IDX_X, IDX_YDOT, IDX_T)
    resid = (IDX_Y, IDX_XDOT, IDX_ZDOT)

    print("=== #412 BCR4BP reach spike (EM-L1 Lyapunov -> full mu_sun) ===")
    print(f"SE-L target reach = {SE_L_REACH_LD:.2f} lunar distances (LD) from Earth")
    best_reach = 0.0
    for c_j, x0 in ((3.18, 0.83), (3.1294, 0.8115), (3.05, 0.79), (2.95, 0.76)):
        try:
            orb = correct_symmetric_fixed_jacobi(cr, x0, c_j, 3.0, ydot0_sign=1.0)
            if not orb.converged:
                print(f"C={c_j}: CR3BP L1 Lyapunov did not converge (x0={x0}); skip")
                continue
            state6 = np.array([orb.x0, 0.0, 0.0, 0.0, orb.ydot0, 0.0])
            seed = correct_bcr4bp_periodic(
                sys_zero,
                state6,
                orb.period,
                sun_commensurate_n=1,
                free_vars=free_vars,
                residual_indices=resid,
                is_half_period_residual=True,
                tol=1e-12,
                independent_tol=1e-6,
            )
            if not seed.converged:
                print(f"C={c_j}: BCR4BP@mu_sun=0 seed failed to close; skip")
                continue
            cr3bp_reach = _max_earth_reach_ld(sys_zero, seed.state_initial, seed.period_nondim)
            fam = continue_bcr4bp_family_in_musun(seed, target_mu_sun=target, n_steps=40)
            members = fam.members
            if not members:
                print(
                    f"C={c_j}: continuation produced no members (cr3bp reach {cr3bp_reach:.2f} LD)"
                )
                continue
            last = members[-1]
            reach = _max_earth_reach_ld(
                last.orbit.system, last.orbit.state_initial, last.orbit.period_nondim
            )
            best_reach = max(best_reach, reach if np.isfinite(reach) else 0.0)
            print(
                f"C={c_j}: cr3bp reach={cr3bp_reach:.2f} LD -> BCR4BP@mu_sun "
                f"reached {last.mu_sun_value:.3e} ({len(members)} members), "
                f"final reach={reach:.2f} LD"
            )
        except Exception as exc:
            print(f"C={c_j}: EXC {type(exc).__name__}: {exc}")

    print(
        f"\n=== VERDICT: best BCR4BP reach = {best_reach:.2f} LD vs SE-L {SE_L_REACH_LD:.2f} LD ==="
    )
    print(
        "REACHABLE -> full build justified"
        if best_reach > 0.5 * SE_L_REACH_LD
        else "STAYS EM-BOUNDED -> full Phase-B build NOT justified from EM-libration family"
    )


if __name__ == "__main__":
    main()
