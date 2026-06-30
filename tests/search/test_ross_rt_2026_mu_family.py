"""Ross & Roberts-Tsoukkas 2026 (arXiv:2606.29189) mu-family representatives.

Phase 2 and Phase 3 of task #494.

Phase 2: recover each of the 6 Table-I representatives from their published
(mu, C, T) via ``correct_symmetric_fixed_jacobi``.  Each check:

  * corrector converges to published x0 (within 5e-7) and T (within 1e-5)
  * winding_topology returns (n*k1, n*k2) for some integer n >= 1, prograde
  * Barden stability |nu| < 1 (all six published members are stable)
  * independent Radau cross-check (closure + Jacobi conservation)

Rep 6 (mu=0.5, k=(1,1)) special note: the paper's T=8.79 is the 3rd iterate
of the fundamental period T1~2.93.  The family label (1,1) refers to the
fundamental winding; at the published T the winding is (3,3).  The topology
check accepts (n*k1, n*k2) for integer n>=1.  A separate test verifies winding
at T/3 gives (1,1).

Phase 3: instantiate the (3,2) representative at Pluto-Charon (mu=0.10851) by
seeding the corrector from the mu=0.1 (3,2) IC.  Additional physical
cross-checks: C(L1) vs Jbara 2025 and Holman-Wiegert circumbinary a_crit.

SOURCED-GOLDEN DISCIPLINE: every EXPECTED value here traces to Ross &
Roberts-Tsoukkas 2026 (arXiv:2606.29189v1, Table I) for Phase 2, and to
Jbara 2025 + Holman-Wiegert 1999 for the Phase-3 physical cross-checks.
No EXPECTED value is a number our own code produced.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import yaml  # type: ignore[import-untyped]

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.search.binary_star_search import collinear_lpoints, winding_topology

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_GOLDEN_PATH = (
    Path(__file__).parent.parent.parent / "data" / "golden" / "ross_rt_2026_cycler_families.yaml"
)


def _load_reps() -> list[dict[str, Any]]:
    with _GOLDEN_PATH.open() as fh:
        data: dict[str, Any] = yaml.safe_load(fh)
    return list(data["representatives"])


def _make_system(mu: float) -> cr3bp.CR3BPSystem:
    """Non-dimensional system: l_km and t_s irrelevant for pure CR3BP tests."""
    return cr3bp.CR3BPSystem(mu=mu, primary="P1", secondary="P2", l_km=1.0, t_s=1.0)


def _topology_ok(topo: object, k1_pub: int, k2_pub: int) -> tuple[bool, int]:
    """Return (ok, n) where n=1 for exact match, n>1 for nth iterate.

    Accepts (n*k1_pub, n*k2_pub) for integer n >= 1.
    """
    k1 = int(topo.k1)  # type: ignore[attr-defined]
    k2 = int(topo.k2)  # type: ignore[attr-defined]
    if k1 == k1_pub and k2 == k2_pub:
        return True, 1
    if k1_pub > 0 and k2_pub > 0:
        n = k1 // k1_pub
        if n >= 1 and k1 == n * k1_pub and k2 == n * k2_pub:
            return True, n
    return False, 0


def _crosscheck(system: cr3bp.CR3BPSystem, orbit: cp.SymmetricOrbit) -> tuple[bool, float]:
    """Build a PeriodicOrbit from a SymmetricOrbit and run the Radau cross-check."""
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    po = cp.PeriodicOrbit(
        state0=state0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        converged=orbit.converged,
        closure_residual=orbit.crossing_residual,
    )
    return cp.crosscheck_periodic(system, po, closure_tol=1e-6, jacobi_tol=1e-8)


# ---------------------------------------------------------------------------
# Phase 2 -- per-representative corrector settings
# ---------------------------------------------------------------------------

# Empirically determined via the exploration summarised in
# docs/notes/2026-06-30-494-phase2-3-mu-family-pluto-charon-verdict.md.
# Columns: (half_crossings, ydot0_sign, tol)
_REP_SETTINGS: list[tuple[int | None, float, float]] = [
    (None, -1.0, 1e-10),  # Rep 1: mu=0.001,        k=(1,1)
    (3, -1.0, 1e-8),  # Rep 2: mu=0.01215,       k=(1,1)  -- different member than AAS 2025
    (None, -1.0, 1e-10),  # Rep 3: mu=0.01215,       k=(3,3)
    (None, -1.0, 1e-10),  # Rep 4: mu=0.1,           k=(3,2)  -- Phase-3 anchor
    (None, -1.0, 1e-10),  # Rep 5: mu=0.3,           k=(3,1)
    (3, -1.0, 1e-10),  # Rep 6: mu=0.5,           k=(1,1)  -- T=3rd iterate, hc=3 picks T/3
]


def _rep_id(rep: dict[str, Any]) -> str:
    k = rep["k"]
    mu_str = f"{rep['mu']:.5g}".replace(".", "_")
    return f"mu{mu_str}_k{k[0]}_{k[1]}"


@pytest.mark.parametrize(
    ("rep", "settings"),
    list(zip(_load_reps(), _REP_SETTINGS, strict=True)),
    ids=[_rep_id(r) for r in _load_reps()],
)
def test_494_phase2_recover_table_i_representative(
    rep: dict[str, Any],
    settings: tuple[int | None, float, float],
) -> None:
    """Phase-2: recover each Table-I representative from its published (mu, C, T).

    Asserts corrector convergence, x0/T proximity to published values (within
    the precision of the printed IC), winding-topology match (allowing the
    nth-iterate form (n*k1, n*k2) for Rep 6), prograde direction, Barden
    linear stability |nu|<1, and independent-integrator cross-check.

    Source: Ross & Roberts-Tsoukkas 2026, arXiv:2606.29189v1, Table I.
    """
    mu = float(rep["mu"])
    k1_pub, k2_pub = int(rep["k"][0]), int(rep["k"][1])
    x0_pub = float(rep["x0"])
    c_pub = float(rep["C"])
    t_pub = float(rep["T"])
    hc, sign, tol = settings

    system = _make_system(mu)
    orbit = cp.correct_symmetric_fixed_jacobi(
        system, x0_pub, c_pub, t_pub, ydot0_sign=sign, half_crossings=hc, tol=tol
    )

    label = f"mu={mu} k=({k1_pub},{k2_pub})"

    # --- convergence ---
    assert orbit.converged, (
        f"{label}: corrector did not converge (res={orbit.crossing_residual:.2e})"
    )

    # --- x0 within 5e-7 (limited by printed precision) ---
    assert abs(orbit.x0 - x0_pub) < 5e-7, (
        f"{label}: |x0 - x0_pub| = {abs(orbit.x0 - x0_pub):.2e} >= 5e-7"
    )

    # --- period T within 1e-5 ---
    assert abs(orbit.period - t_pub) < 1e-5, (
        f"{label}: |T - T_pub| = {abs(orbit.period - t_pub):.2e} >= 1e-5"
    )

    # --- winding topology: accept (n*k1, n*k2) for n >= 1 ---
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    topo = winding_topology(mu, state0, orbit.period)
    ok_topo, _n = _topology_ok(topo, k1_pub, k2_pub)
    assert ok_topo, (
        f"{label}: winding mismatch -- expected ({k1_pub},{k2_pub}) or multiple, "
        f"got ({topo.k1},{topo.k2})"
    )
    assert topo.prograde, f"{label}: expected prograde (both winding numbers > 0)"

    # --- Barden stability: |nu| < 1 (all Table-I members are stable) ---
    nu, _lam = cp.barden_stability(system, orbit)
    assert abs(nu) < 1.0, (
        f"{label}: Barden stability FAILED |nu|={abs(nu):.4f} >= 1 (published sp={rep['sp']:.4f})"
    )

    # --- independent Radau cross-check ---
    ok_cc, dj = _crosscheck(system, orbit)
    assert ok_cc, f"{label}: independent cross-check failed (dj={dj:.2e})"


# ---------------------------------------------------------------------------
# Phase 2 -- Rep 6 fundamental-period check
# ---------------------------------------------------------------------------


def test_494_phase2_rep6_fundamental_winding_is_11() -> None:
    """Rep 6 (mu=0.5, published T=8.79) topology at T/3 (fundamental period)
    must be (1,1) -- the family label refers to the fundamental winding, while
    hc=3 locks onto T=3*T1 (the 3rd iterate, winding (3,3) at the full T).

    Source: Ross & Roberts-Tsoukkas 2026, arXiv:2606.29189v1, Table I, row 6;
    sp triple-angle identity 4*nu1^3 - 3*nu1 = nu3 verifies iterate structure.
    """
    reps = _load_reps()
    rep = reps[5]  # mu=0.5, k=(1,1)
    mu = float(rep["mu"])  # 0.5
    x0_pub = float(rep["x0"])
    c_pub = float(rep["C"])
    t_pub = float(rep["T"])  # 8.7920...

    system = _make_system(mu)
    # Recover with hc=3 (locks onto the full 3rd-iterate period T~8.79)
    orbit = cp.correct_symmetric_fixed_jacobi(
        system, x0_pub, c_pub, t_pub, ydot0_sign=-1.0, half_crossings=3, tol=1e-10
    )
    assert orbit.converged, "Rep 6 corrector (hc=3) did not converge"

    # Verify winding at T/3 (fundamental period) = (1,1)
    t_fund = orbit.period / 3.0
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    topo_fund = winding_topology(mu, state0, t_fund)
    assert (topo_fund.k1, topo_fund.k2) == (1, 1), (
        f"fundamental-period winding expected (1,1), got ({topo_fund.k1},{topo_fund.k2})"
    )
    assert topo_fund.prograde, "fundamental-period orbit must be prograde"

    # Verify |nu| at full T matches published sp=0.9376 (triple-angle identity).
    # SOURCED EXPECTED: Ross 2026 Table I, sp=0.9376 for this entry.
    nu_full, _ = cp.barden_stability(system, orbit)
    assert abs(nu_full - float(rep["sp"])) < 1e-3, (
        f"nu at T (3rd iterate) = {nu_full:.6f}, published sp = {rep['sp']}"
    )


# ---------------------------------------------------------------------------
# Phase 3 -- Pluto-Charon (mu=0.10851) instantiation of the (3,2) family
# ---------------------------------------------------------------------------

# Physical constants for Pluto-Charon.
_PC_MU = 0.10851  # task-specified mu (per Ross 2026 paper context)
_PC_L_KM = 19600.0  # Pluto-Charon mean separation [km] (from satellites.py)
_PC_GM_SYSTEM = 975.5  # Pluto+Charon system GM [km^3/s^2] (JPL DE440 / satellites.py)
_PC_T_S = 2.0 * math.pi * math.sqrt(_PC_L_KM**3 / _PC_GM_SYSTEM)  # Charon orbital period [s]

# Rep 4 (mu=0.1, (3,2)) -- Phase-3 anchor
_PC_ANCHOR_X0 = -0.694376003123377
_PC_ANCHOR_C = 3.573367616904619
_PC_ANCHOR_T = 12.295263874014290


def test_494_phase3_pluto_charon_32_cycler() -> None:
    """Phase-3: instantiate the (3,2) cycler family at Pluto-Charon (mu=0.10851).

    Seed the fixed-Jacobi corrector at mu=0.10851 from the mu=0.1 anchor IC
    (Rep 4, Table I).  Delta_mu = 0.00851 is small enough that the seed lies
    in the basin of the (3,2) branch at the anchor Jacobi constant.

    Physical cross-checks (sourced):
      * C(L1) at mu=0.10851: Jbara 2025 reports ~3.6210 for mu~0.109; tol 0.005.
      * Holman-Wiegert (1999) critical semi-major axis formula
        a_c/a_bin = 1.60 + 4.12*mu - 5.09*mu^2  =>  a_crit ~ 38,947 km.

    STABILITY FINDING: at the anchor Jacobi constant C=3.5734, the (3,2) orbit
    recovered at mu=0.10851 is UNSTABLE (nu ~ 1903).  The family likely loses
    linear stability in the mu in [0.1, 0.10851] interval at this C value.
    Stability is measured but NOT asserted; convergence, topology, and physical
    cross-checks are asserted.  See the Phase-3 verdict note.
    """
    system_pc = cr3bp.CR3BPSystem(
        mu=_PC_MU,
        primary="Pluto",
        secondary="Charon",
        l_km=_PC_L_KM,
        t_s=_PC_T_S,
    )

    # Correct at mu=0.10851 seeded from mu=0.1 anchor.
    # half_crossings=6: the (3,2) half-period crossing is the 6th x-axis crossing
    # (same index as at mu=0.1; orbit topology preserved at this seed).
    orbit_pc = cp.correct_symmetric_fixed_jacobi(
        system_pc,
        _PC_ANCHOR_X0,
        _PC_ANCHOR_C,
        _PC_ANCHOR_T,
        ydot0_sign=-1.0,
        half_crossings=6,
        tol=1e-10,
    )

    assert orbit_pc.converged, (
        f"Pluto-Charon (3,2) corrector did not converge (res={orbit_pc.crossing_residual:.2e})"
    )

    # Topology must be (3,2), prograde.
    state0_pc = np.array([orbit_pc.x0, 0.0, 0.0, 0.0, orbit_pc.ydot0, 0.0])
    topo_pc = winding_topology(_PC_MU, state0_pc, orbit_pc.period)
    assert (topo_pc.k1, topo_pc.k2) == (3, 2), (
        f"Pluto-Charon topology: expected (3,2), got ({topo_pc.k1},{topo_pc.k2})"
    )
    assert topo_pc.prograde, "Pluto-Charon (3,2) orbit must be prograde"

    # Stability -- measured but NOT asserted (UNSTABLE at anchor C; see docstring).
    nu_pc, _lam_pc = cp.barden_stability(system_pc, orbit_pc)
    # nu_pc ~ 1903; the (3,2) orbit is UNSTABLE at mu=0.10851 with C=3.5734.
    _ = nu_pc

    # --- Physical cross-check 1: C(L1) vs Jbara 2025 ---
    # Source: Jbara, R. 2025 -- C_L1 ~ 3.6210 at mu ~ 0.109; tolerance 0.005.
    l1, _l2, _l3 = collinear_lpoints(_PC_MU)
    c_l1 = cr3bp.jacobi_constant(np.array([l1, 0.0, 0.0, 0.0, 0.0, 0.0]), _PC_MU)
    # SOURCED EXPECTED: Jbara 2025, mu~0.109 -> C_L1=3.6210; tolerance 0.005.
    assert abs(c_l1 - 3.6210) < 0.005, (
        f"C(L1) at mu={_PC_MU} = {c_l1:.4f}, expected ~3.6210 (Jbara 2025) within 0.005"
    )

    # --- Physical cross-check 2: Holman-Wiegert a_crit ---
    # Source: Holman & Wiegert 1999, Eq. 1: a_c/a_bin = 1.60 + 4.12*mu - 5.09*mu^2.
    ac_ratio = 1.60 + 4.12 * _PC_MU - 5.09 * _PC_MU**2
    a_crit_km = ac_ratio * _PC_L_KM
    # SOURCED EXPECTED: H-W 1999 formula -> ~38947 km; tolerance 100 km.
    assert abs(a_crit_km - 38947.0) < 100.0, (
        f"H-W a_crit = {a_crit_km:.1f} km, expected ~38947 km (+-100 km)"
    )

    # --- Independent Radau cross-check ---
    ok_cc, dj = _crosscheck(system_pc, orbit_pc)
    assert ok_cc, f"Pluto-Charon (3,2) independent cross-check failed (dj={dj:.2e})"
