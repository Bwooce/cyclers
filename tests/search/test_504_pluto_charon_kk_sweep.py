"""Task #504: Pluto-Charon (k1,k2) family sweep test suite.

Sweeps (1,1), (2,1), (3,1), (2,2), (3,3) at Pluto-Charon
mu = cr3bp_system("Pluto","Charon").mu = 0.10876473603280369, seeding from
the nearest Ross-RT 2026 (arXiv:2606.29189) Table-I anchor.

POSITIVE CONTROL: re-finds the admitted (3,2) stable member from #494 at the
exact mu from cr3bp_system (slightly different from 0.10851 used earlier).
Seeds from the mu=0.1 Table-I anchor → C-sweep → brentq nu=0 → C ≈ 3.5795.

OTHER FAMILIES: clean-negative-aware.  If no stable window is found the test
PASSES (records a clean negative); if a stable member IS found, all validity
gates must hold.
"""

from __future__ import annotations

import numpy as np

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.binary_star_search import collinear_lpoints
from cyclerfinder.search.pluto_charon_kk_sweep import (
    PC_MU,
    SweepResult,
    make_pluto_charon_system,
    sweep_11,
    sweep_21,
    sweep_22,
    sweep_31,
    sweep_32_positive_control,
    sweep_33,
)

# ---------------------------------------------------------------------------
# Shared assertion helpers
# ---------------------------------------------------------------------------


def _c_l1_pc() -> float:
    l1, _l2, _l3 = collinear_lpoints(PC_MU)
    return float(cr3bp.jacobi_constant(np.array([l1, 0.0, 0.0, 0.0, 0.0, 0.0]), PC_MU))


def _assert_valid_stable_member(result: SweepResult, *, label: str) -> None:
    """Assert all validity gates for a found stable member."""
    assert result.stable_found, f"{label}: stable_found=False but asserted True"
    assert result.topology_ok, f"{label}: topology ({result.k1},{result.k2}) check failed"
    assert result.prograde, f"{label}: orbit must be prograde"
    assert result.reaches_secondary, f"{label}: orbit must reach the secondary realm (L1 passage)"
    assert result.nu_mid is not None
    assert abs(result.nu_mid) < 1.0, f"{label}: |nu|={abs(result.nu_mid):.4f} >= 1 (not stable)"
    assert result.crosscheck_ok, (
        f"{label}: independent Radau crosscheck failed (dJ={result.crosscheck_dj:.2e})"
    )


# ---------------------------------------------------------------------------
# Positive control: (3,2) at mu=PC_MU
# ---------------------------------------------------------------------------


def test_504_positive_control_32() -> None:
    """POSITIVE CONTROL: re-find the stable (3,2) member at PC mu=0.10876.

    Seeds from the Ross-RT 2026 Table-I (mu=0.1, (3,2)) anchor:
      x0=-0.694376003123377, C=3.573367616904619, T=12.295263874014290
    C-sweeps upward to locate the nu=0 midpoint (C ≈ 3.5795, derived).

    Sourced cross-checks:
      Jbara 2025 (arXiv:2510.13479): C_L1 ≈ 3.6210 at mu~0.109; tol 0.005.
      Holman & Wiegert 1999, Eq. 1: a_crit ≈ 38948 km; tol 150 km.

    This is a HARD assertion: if it fails, the sweep machinery is broken.
    """
    result = sweep_32_positive_control()

    assert result.stable_found, (
        "Positive control (3,2) FAILED — sweep machinery is broken. "
        f"method={result.method!r}, note={result.note!r}"
    )
    _assert_valid_stable_member(result, label="(3,2) positive control")

    # C close to derived nu=0 midpoint (~3.5795); allow ±0.002 for seed variation.
    assert result.jacobi_mid is not None
    assert abs(result.jacobi_mid - 3.5795) < 0.002, (
        f"(3,2) jacobi_mid={result.jacobi_mid:.6f} far from expected 3.5795"
    )

    # C < C_L1 (orbit must be below the Hill threshold)
    c_l1 = _c_l1_pc()
    assert result.jacobi_mid < c_l1, f"(3,2): C={result.jacobi_mid:.6f} >= C_L1={c_l1:.6f}"

    # Period sanity: Charon P ≈ 6.387 d; (3,2) orbit ≈ 12 d
    assert result.period_days is not None
    assert 10.0 < result.period_days < 20.0, (
        f"(3,2): period_days={result.period_days:.2f} d out of expected [10,20]"
    )

    # Cross-check against sourced values at PC mu
    # SOURCED: Jbara 2025, arXiv:2510.13479 — C_L1 ≈ 3.6210 at mu~0.109
    assert abs(c_l1 - 3.6210) < 0.005, f"C_L1={c_l1:.4f} vs Jbara 2025 3.6210 (tol 0.005)"
    # SOURCED: Holman & Wiegert 1999, Eq. 1
    pc = make_pluto_charon_system()
    ac_ratio = 1.60 + 4.12 * PC_MU - 5.09 * PC_MU**2
    a_crit_km = ac_ratio * pc.l_km
    assert abs(a_crit_km - 38948.0) < 150.0, (
        f"H-W a_crit={a_crit_km:.1f} km, expected ~38948 km (tol 150 km)"
    )


# ---------------------------------------------------------------------------
# (1,1) sweep
# ---------------------------------------------------------------------------


def test_504_sweep_11() -> None:
    """Sweep the (1,1) family at PC mu.

    Seeds from the Ross-RT 2026 Table-I (mu=0.001, (1,1)) anchor via
    mu-continuation, then C-sweeps for the stable window.

    Clean-negative-aware: passes either way.
    """
    result = sweep_11()

    if not result.stable_found:
        return  # clean negative

    _assert_valid_stable_member(result, label="(1,1)")
    assert result.jacobi_mid is not None
    assert result.jacobi_mid < _c_l1_pc(), (
        f"(1,1): C={result.jacobi_mid:.6f} >= C_L1={_c_l1_pc():.6f}"
    )


# ---------------------------------------------------------------------------
# (2,1) sweep
# ---------------------------------------------------------------------------


def test_504_sweep_21() -> None:
    """Sweep the (2,1) family at PC mu (grid search — not in Table-I).

    Clean-negative-aware: passes either way.
    """
    result = sweep_21()

    if not result.stable_found:
        return  # clean negative

    _assert_valid_stable_member(result, label="(2,1)")
    assert result.jacobi_mid is not None
    assert result.jacobi_mid < _c_l1_pc(), (
        f"(2,1): C={result.jacobi_mid:.6f} >= C_L1={_c_l1_pc():.6f}"
    )


# ---------------------------------------------------------------------------
# (3,1) sweep
# ---------------------------------------------------------------------------


def test_504_sweep_31() -> None:
    """Sweep the (3,1) family at PC mu.

    The mu=0.3 anchor has C=3.702 > C_L1(PC)=3.621: expected clean negative
    because the orbit cannot reach the secondary at PC mu at this energy level.
    Both direct mu-step and C-walk strategies are attempted.

    Clean-negative-aware: passes either way.
    """
    result = sweep_31()

    if not result.stable_found:
        return  # expected clean negative

    _assert_valid_stable_member(result, label="(3,1)")
    assert result.jacobi_mid is not None
    assert result.jacobi_mid < _c_l1_pc(), (
        f"(3,1): C={result.jacobi_mid:.6f} >= C_L1={_c_l1_pc():.6f}"
    )


# ---------------------------------------------------------------------------
# (2,2) sweep
# ---------------------------------------------------------------------------


def test_504_sweep_22() -> None:
    """Sweep the (2,2) family at PC mu (grid search — not in Table-I).

    Clean-negative-aware: passes either way.
    """
    result = sweep_22()

    if not result.stable_found:
        return  # clean negative

    _assert_valid_stable_member(result, label="(2,2)")
    assert result.jacobi_mid is not None
    assert result.jacobi_mid < _c_l1_pc(), (
        f"(2,2): C={result.jacobi_mid:.6f} >= C_L1={_c_l1_pc():.6f}"
    )


# ---------------------------------------------------------------------------
# (3,3) sweep
# ---------------------------------------------------------------------------


def test_504_sweep_33() -> None:
    """Sweep the (3,3) family at PC mu.

    Seeds from the Ross-RT 2026 Table-I (mu=0.012150, (3,3)) anchor via
    mu-continuation, then C-sweeps for the stable window.

    Clean-negative-aware: passes either way.
    """
    result = sweep_33()

    if not result.stable_found:
        return  # clean negative

    _assert_valid_stable_member(result, label="(3,3)")
    assert result.jacobi_mid is not None
    assert result.jacobi_mid < _c_l1_pc(), (
        f"(3,3): C={result.jacobi_mid:.6f} >= C_L1={_c_l1_pc():.6f}"
    )
