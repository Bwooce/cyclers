"""Evidence tests for task #511: Pluto-Charon (3,2) real-ephemeris lever.

Kernel-backed tests (skipped when ``plu060.bsp`` is absent -- local-only,
129 MB, never committed; mirrors the JUP365 skip convention in
``tests/nbody/test_jovian_shoot.py``): SPICE osculating-element extraction
and the ER3BP differential correction.

Kernel-free tests (always run): the structural period-incommensurability
check and the differential corrector's residual bookkeeping against a
synthetic seed, which do not need SPICE.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

pytest.importorskip("scipy")

from cyclerfinder.core.cr3bp import CR3BPSystem, cr3bp_system
from cyclerfinder.search.cr3bp_periodic import SymmetricOrbit
from cyclerfinder.verify.pluto_charon_realeph import (
    differential_correct_pc32_to_eccentricity,
)
from cyclerfinder.verify.spice_kernels import ensure_pluto_kernel

try:
    _KERNEL: str | None = ensure_pluto_kernel()
except Exception:  # plu060.bsp is local-only (129 MB, absent in CI) -> skip, don't error
    _KERNEL = None

_needs_kernel = pytest.mark.skipif(_KERNEL is None, reason="plu060.bsp not furnished (local-only)")

# Catalogue parameters (SOURCED: same as scripts/pc_v2_longspan.py, #494/#504/#505)
PC_C = 3.57951501972907
PC_X0 = -0.693198287043369
PC_T_GUESS = 11.8334625170346


def _pc32_seed() -> tuple[CR3BPSystem, SymmetricOrbit]:
    import cyclerfinder.search.cr3bp_periodic as cp

    system = cr3bp_system("Pluto", "Charon")
    orbit = cp.correct_symmetric_fixed_jacobi(
        system, PC_X0, PC_C, PC_T_GUESS, ydot0_sign=-1.0, half_crossings=6, tol=1e-10
    )
    assert orbit.converged
    return system, orbit


# --- kernel-free: structural incommensurability -----------------------------------


def test_pc32_period_is_not_commensurate_with_charon_period() -> None:
    """PC (3,2) seed period / 2*pi is far from an integer (the #441 trap precondition).

    T_nd = 11.8334625 TU means the spacecraft period is 1.883 Charon periods
    -- not an integer multiple. Under the ER3BP pulsating frame (2*pi-periodic
    in true anomaly f), exact periodicity at e>0 requires the period in f to
    be a multiple of 2*pi (docstring, ``genome/er3bp_periodic.py``; durable
    finding, #441). This test pins the ratio so the #511 verdict's structural
    claim is regression-checked independent of the (kernel-gated) numerical
    confirmation below.
    """
    _system, orbit = _pc32_seed()
    ratio = orbit.period / (2.0 * math.pi)
    nearest_int = round(ratio)
    assert abs(ratio - nearest_int) > 0.05, (
        f"period ratio {ratio:.6f} unexpectedly close to integer {nearest_int} "
        "-- the period_f trap precondition no longer holds"
    )


def test_differential_correct_e0_bridge_is_exact() -> None:
    """e=0 ER3BP bridge reproduces the CR3BP seed with near-machine residual.

    Mirrors the #441 Sec. 1 finding ("At e=0 the ER3BP IS the CR3BP, so the
    bridge is exact") for the PC (3,2) orbit specifically. Kernel-free: e=0
    needs no SPICE input.
    """
    system, orbit = _pc32_seed()
    result = differential_correct_pc32_to_eccentricity(
        0.0,
        x0_seed=orbit.x0,
        ydot0_seed=orbit.ydot0,
        period_seed=orbit.period,
        mu=system.mu,
        n_steps=1,
    )
    assert result.seed_corrector_residual < 1e-8
    assert result.seed_independent_residual < 1e-6


# --- kernel-backed: the real-eph differential correction ---------------------------


@_needs_kernel
def test_pc32_charon_real_osculating_eccentricity_is_tiny_and_stable() -> None:
    """SPICE-derived Charon osculating eccentricity: O(1e-4), epoch-stable.

    Cross-checks #506's back-of-envelope ``e < 5e-5`` (Brozovic et al. 2015
    MEAN eccentricity) against the actual SPICE osculating value at several
    epochs -- confirms both that plu060.bsp reads sanely (Charon's SMA near
    the 19600 km catalogue value) and that the real eccentricity, while a few
    times the published mean value, is still very small.
    """
    from cyclerfinder.verify.pluto_charon_realeph import charon_osculating_elements

    assert _KERNEL is not None  # guaranteed by @_needs_kernel decorator
    epochs = ["2000-01-01T12:00:00", "2026-07-01T00:00:00", "2050-01-01T00:00:00"]
    elems = [charon_osculating_elements(iso, _KERNEL) for iso in epochs]
    for el in elems:
        assert 19000.0 < el.a_km < 20200.0, (
            f"Charon SMA {el.a_km} km off the ~19600 km catalogue value"
        )
        assert 0.0 < el.eccentricity < 2e-3, (
            f"osculating e={el.eccentricity:.3e} outside the expected O(1e-4) band"
        )
        assert 6.3 < el.period_days < 6.5
    e_vals = np.array([el.eccentricity for el in elems])
    assert float(np.ptp(e_vals)) < 1e-4, "osculating e varies too much across a 50-year epoch span"


@_needs_kernel
def test_pc32_realeph_differential_correction_does_not_close() -> None:
    """The #511 headline result: PC (3,2) has NO strictly-periodic real-eph analog.

    Runs the full differential-correction pipeline (SPICE real e -> ER3BP
    e-continuation -> independent-closure gate) and pins the false-positive-
    closure signature: the symmetric corrector residual converges (it only
    tests the half-period crossing condition) while the independent
    full-orbit Radau residual stays far above the closure gate -- the #441
    period_f-trap signature, now confirmed on real SPICE-sourced input.
    """
    from cyclerfinder.verify.pluto_charon_realeph import charon_osculating_elements

    assert _KERNEL is not None  # guaranteed by @_needs_kernel decorator
    charon = charon_osculating_elements("2026-07-01T00:00:00", _KERNEL)
    system, orbit = _pc32_seed()

    result = differential_correct_pc32_to_eccentricity(
        charon.eccentricity,
        x0_seed=orbit.x0,
        ydot0_seed=orbit.ydot0,
        period_seed=orbit.period,
        mu=system.mu,
        n_steps=20,
        independent_tol=1e-8,
    )

    assert result.orbit is not None, (
        "continuation should reach the target e (corrector residual converges)"
    )
    assert result.target_corrector_residual < 1e-8, (
        f"expected the symmetric corrector to (falsely) report convergence, "
        f"got residual={result.target_corrector_residual:.3e}"
    )
    assert result.target_independent_residual > 1e-4, (
        f"expected the independent-closure gate to catch non-closure "
        f"(residual={result.target_independent_residual:.3e}); if this now passes, "
        "the real-eph analog orbit may actually exist -- re-open the #511 verdict"
    )
    assert not result.converged
