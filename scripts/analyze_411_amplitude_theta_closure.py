"""#411 amplitude-knob theta-closure feasibility for the cross-system cycle.

Supersedes the fixed-amplitude (n_em=41, n_se=19) feasibility map
(analyze_405_theta_closure.py). That map held both orbit amplitudes FIXED and
forced integer revolution counts to absorb the entire ~2.2 rad single-rev
theta-time-consistency gap -- which demands 41 EM-orbit + 19 SE-orbit revolutions.
This script shows that is the WRONG model: it is physically infeasible and
unnecessary.

INFEASIBLE at fixed amplitude: a manifold seeded at epsilon can only SHADOW its
orbit for ~ ln(eps_range)/ln(|lambda_u|) revolutions before exponential departure.
For the EM-L2 Lyapunov |lambda_u| ~ 1.2e3 (ln ~ 7.1/rev), so even with eps spanning
12 decades the shadow budget is only ~3-4 revolutions -- you cannot loiter 41 revs
on an orbit that unstable; the trajectory has long since left.

UNNECESSARY -- the amplitude is a continuous knob: over each Lyapunov family the
per-revolution relative-phase advance Delta-theta(C) mod 2pi sweeps the FULL [0,2pi)
(the SE term is ~38 rad before the mod, so a tiny SE-period change sweeps the whole
circle). With n_em = n_se = 1, the closure condition

    gap(c_em, c_se) + Delta-theta_em(c_em) + Delta-theta_se(c_se) == 0  (mod 2pi)

is ONE equation in TWO continuous amplitude knobs, each of which alone sweeps
[0,2pi). A 1-D solution CURVE therefore generically exists at single revolution,
inside the shadow budget. This is a NECESSARY-condition (precondition) result; it
asserts no closed orbit. Sufficiency -- whether both connection legs still converge
at low patch-dV ON that curve -- is the open question the full corrector must answer.

Run: uv run python scripts/analyze_411_amplitude_theta_closure.py
"""

from __future__ import annotations

import math

import numpy as np

from cyclerfinder.genome.cross_system_cycle import em_moon_system, se_earth_system
from cyclerfinder.genome.heteroclinic_cycle import LyapunovNode, _planar_floquet_pair

# Single-revolution theta-time-consistency gap, from the committed #411 analysis.
SINGLE_REV_GAP_RAD = 2.2
# eps spans ~12 decades (linear-seed floor ~1e-13 up to a macroscopic ~1e-1 departure).
EPS_RANGE = 1e12
TWO_PI = 2.0 * math.pi


def _shadow_budget(lam_u: float) -> float:
    """Max revolutions a manifold seeded at eps can shadow its orbit before departure."""
    return math.log(EPS_RANGE) / math.log(lam_u) if lam_u > 1.0 + 1e-12 else float("inf")


def _family_map(system, x0_guess, jacobis, period_guess, label, sgn, omega_rel):
    """Return per-amplitude (C, period, dtheta_mod2pi, |lam_u|, shadow_budget) rows."""
    rows = []
    for c in jacobis:
        try:
            node = LyapunovNode.from_libration(
                system,
                x0_guess=x0_guess,
                jacobi=float(c),
                period_guess=period_guess,
                label=label,
                ydot0_sign=sgn,
            )
        except (RuntimeError, ValueError):
            continue
        if not node.converged:
            continue
        lam_u = _planar_floquet_pair(system, node.state0, node.period)[0]
        dtheta = (omega_rel * node.period * system.t_s) % TWO_PI
        rows.append((float(c), node.period, dtheta, lam_u, _shadow_budget(lam_u)))
    return rows


def main() -> None:
    se = se_earth_system()
    em = em_moon_system()
    omega_rel = 1.0 / em.t_s - 1.0 / se.t_s  # rad/s, relative SE-EM line rate

    em_rows = _family_map(em, 1.18, np.linspace(3.108, 3.153, 12), 3.4, "EM-L2", -1.0, omega_rel)
    se_rows = _family_map(
        se, 1.009, np.linspace(3.0000, 3.0008, 17), 3.06, "SE-L2", -1.0, omega_rel
    )

    print("=== #411 amplitude-knob theta-closure feasibility ===")
    print(f"omega_rel = {omega_rel:.6e} rad/s; single-rev gap = {SINGLE_REV_GAP_RAD:.3f} rad\n")

    for name, rows in (("EM-L2", em_rows), ("SE-L2", se_rows)):
        if not rows:
            print(f"{name}: no converged members")
            continue
        dthetas = [r[2] for r in rows]
        budgets = [r[4] for r in rows]
        coverage = (max(dthetas) - min(dthetas)) / TWO_PI
        print(
            f"{name}: {len(rows)} members | dtheta mod 2pi span "
            f"[{min(dthetas):.3f}, {max(dthetas):.3f}] rad = {coverage * 100:.0f}% of 2pi | "
            f"shadow budget [{min(budgets):.1f}, {max(budgets):.1f}] rev"
        )

    # Fixed-amplitude verdict (why 41/19 is infeasible).
    em_bud = np.median([r[4] for r in em_rows]) if em_rows else float("nan")
    print(
        f"\nFixed-amplitude (41,19): n_em=41 vs EM budget ~{em_bud:.1f} rev -> "
        f"{'OK' if em_bud >= 41 else 'INFEASIBLE (>10x over budget)'}"
    )

    # Dual-knob verdict. The closure equation gap + n_em*dtheta_em + n_se*dtheta_se == 0
    # (mod 2pi) is ONE equation in TWO continuous amplitude knobs. If EITHER family's
    # dtheta(C) sweeps (essentially) the full circle WITHIN its shadow budget, then for
    # any fixed amplitude on the other side, sweeping that family crosses 0 -> a single-
    # rev (n=1) closure curve exists. One covering knob within budget is sufficient.
    def _covers(rows):
        return bool(rows) and (max(r[2] for r in rows) - min(r[2] for r in rows)) > 0.85 * TWO_PI

    em_cover = _covers(em_rows)
    em_in_budget = bool(em_rows) and min(r[4] for r in em_rows) >= 1.0  # >=1 rev shadowable
    se_cover = _covers(se_rows)
    se_in_budget = bool(se_rows) and min(r[4] for r in se_rows) >= 1.0
    em_knob = em_cover and em_in_budget
    se_knob = se_cover and se_in_budget
    print(
        f"Dual-knob (n=1): EM is a full-circle knob within budget ({'YES' if em_knob else 'NO'}); "
        f"SE is a full-circle knob within budget ({'YES' if se_knob else 'NO'})."
    )
    if em_knob or se_knob:
        print(
            "=> VERDICT: FEASIBLE at single revolution. At least one amplitude knob sweeps "
            "[0,2pi) within the shadow budget, so the theta-closure equation (1 eq, 2 "
            "unknowns) has a 1-D solution curve. Build the coupled corrector; sufficiency "
            "(low patch-dV on the curve) is the open question."
        )
    else:
        print("=> VERDICT: precondition NOT met; closure curve not guaranteed at single rev.")


if __name__ == "__main__":
    main()
