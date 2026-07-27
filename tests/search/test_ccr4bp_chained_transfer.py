"""Tests for the CCR4BP chained-transfer glue module (`#715`).

Second independent positive control for the CCR4BP heteroclinic-search
pipeline, grounded in Aryan & Fitzgerald (2024) AAS 24-103's Callisto-L1 ->
Ganymede-flyby -> Europa-L2 transit (digested at `#710`). See
``src/cyclerfinder/search/ccr4bp_chained_transfer.py``'s own module
docstring for the full account of what was and was not reproduced.

Structure, mirroring `#694`'s own ``test_ccr4bp_heteroclinic_search.py``
discipline: (1) the new generic Lyapunov-seed builder, cross-checked against
the paper's own quoted Jacobi constants for BOTH systems; (2) the torus
corrector applied to those seeds -- a genuine SECOND positive control for
the torus corrector, independent of `#690`'s own resonant-orbit-seeded JEG
torus; (3) manifold globalization succeeds mechanically on both (new,
previously-unexercised) systems; (4) the closest-approach-to-Ganymede scan,
whose live-observed result is a HONEST, DOCUMENTED negative within this
task's own search scope (a modest theta2 grid, ~90-day propagation window)
-- pinned as a regression floor, not overclaimed as "no connection exists."
"""

from __future__ import annotations

from typing import Literal

import numpy as np
import pytest

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.search.ccr4bp_chained_transfer as ct
import cyclerfinder.search.ccr4bp_manifold_globalize as mg
import cyclerfinder.search.variational_ccr4bp_torus as vt
from cyclerfinder.core.ccr4bp_callisto_ganymede import (
    L_KM as CALLISTO_L_KM,
)
from cyclerfinder.core.ccr4bp_callisto_ganymede import (
    jupiter_callisto_ganymede_default,
)
from cyclerfinder.genome.composed_moon_map import moon_config

_EUROPA_L_KM = moon_config("Europa").sma_km
_CALLISTO_PERIOD_DAYS = moon_config("Callisto").period_s / 86400.0
_EUROPA_PERIOD_DAYS = moon_config("Europa").period_s / 86400.0

# Paper's own quoted Table 1/2 Jacobi constants (#710 digest section 1/3),
# both systems' "2:1 resonant CR3BP planar Lyapunov" seed families.
_CALLISTO_JACOBI = 3.0044
_EUROPA_JACOBI = 3.0034

_TORUS_N1 = 1
_TORUS_N2 = 20
_TORUS_TR_SOLVER: Literal["exact", "lsmr"] = "exact"
_TORUS_MAX_NFEV = 600
_TORUS_GAUGE_WEIGHT = 30.0
_TORUS_RHO_WEIGHT = 100.0


# ---------------------------------------------------------------------------
# lyapunov_seed_at_jacobi: generic L1/L2 seed builder, cross-checked against
# the paper's own quoted Jacobi constants.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def callisto_system() -> ccr4bp.CCR4BPSystem:
    return jupiter_callisto_ganymede_default()


@pytest.fixture(scope="module")
def europa_system() -> ccr4bp.CCR4BPSystem:
    return ccr4bp.jupiter_europa_ganymede_default()


@pytest.fixture(scope="module")
def callisto_l1_seed(callisto_system: ccr4bp.CCR4BPSystem) -> ct.LyapunovSeed:
    return ct.lyapunov_seed_at_jacobi(
        callisto_system.mu, _CALLISTO_JACOBI, "L1", x0_offset=0.005, ydot0_sign=1.0
    )


@pytest.fixture(scope="module")
def europa_l2_seed(europa_system: ccr4bp.CCR4BPSystem) -> ct.LyapunovSeed:
    return ct.lyapunov_seed_at_jacobi(
        europa_system.mu, _EUROPA_JACOBI, "L2", x0_offset=0.002, ydot0_sign=-1.0
    )


def test_callisto_l1_seed_hits_paper_jacobi_exactly(callisto_l1_seed: ct.LyapunovSeed) -> None:
    assert callisto_l1_seed.jacobi == pytest.approx(_CALLISTO_JACOBI, abs=1e-9)
    assert callisto_l1_seed.orbit.converged
    assert callisto_l1_seed.orbit.crossing_residual < 1e-9
    # A genuine small-to-moderate-amplitude Lyapunov period, not some other
    # far family member (sanity range, not an exact target).
    assert 1.0 < callisto_l1_seed.period < 6.0


def test_europa_l2_seed_hits_paper_jacobi_exactly(europa_l2_seed: ct.LyapunovSeed) -> None:
    assert europa_l2_seed.jacobi == pytest.approx(_EUROPA_JACOBI, abs=1e-9)
    assert europa_l2_seed.orbit.converged
    assert europa_l2_seed.orbit.crossing_residual < 1e-9
    assert 1.0 < europa_l2_seed.period < 6.0


def test_lyapunov_seed_rejects_bad_l_point(callisto_system: ccr4bp.CCR4BPSystem) -> None:
    with pytest.raises(ValueError, match="L1"):
        ct.lyapunov_seed_at_jacobi(callisto_system.mu, _CALLISTO_JACOBI, "L3")


# ---------------------------------------------------------------------------
# ganymede_synodic_position cross-check against core.ccr4bp's own private
# helper (must stay byte-for-byte identical -- see module docstring).
# ---------------------------------------------------------------------------


def test_ganymede_synodic_position_matches_private_helper(
    callisto_system: ccr4bp.CCR4BPSystem, europa_system: ccr4bp.CCR4BPSystem
) -> None:
    for system in (callisto_system, europa_system):
        for t in np.linspace(0.0, 10.0, 13):
            gx, gy = ct.ganymede_synodic_position(float(t), system)
            gx_ref, gy_ref, gz_ref = ccr4bp._ganymede_position(float(t), system)
            assert gx == pytest.approx(gx_ref, abs=1e-14)
            assert gy == pytest.approx(gy_ref, abs=1e-14)
            assert gz_ref == 0.0


# ---------------------------------------------------------------------------
# Torus corrector: SECOND positive control, independent of #690's own
# resonant-orbit-seeded JEG torus -- genuine L1/L2 Lyapunov seeds at the
# paper's own exact Jacobi constants.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def callisto_l1_torus(
    callisto_system: ccr4bp.CCR4BPSystem, callisto_l1_seed: ct.LyapunovSeed
) -> vt.CCR4BPTorusVariationalResult:
    return vt.discover_ccr4bp_torus_from_resonant_orbit(
        callisto_system,
        callisto_l1_seed.state0,
        callisto_l1_seed.period,
        n1=_TORUS_N1,
        n2=_TORUS_N2,
        tr_solver=_TORUS_TR_SOLVER,
        max_nfev=_TORUS_MAX_NFEV,
        gauge_weight=_TORUS_GAUGE_WEIGHT,
        rho_weight=_TORUS_RHO_WEIGHT,
    )


@pytest.fixture(scope="module")
def europa_l2_torus(
    europa_system: ccr4bp.CCR4BPSystem, europa_l2_seed: ct.LyapunovSeed
) -> vt.CCR4BPTorusVariationalResult:
    return vt.discover_ccr4bp_torus_from_resonant_orbit(
        europa_system,
        europa_l2_seed.state0,
        europa_l2_seed.period,
        n1=_TORUS_N1,
        n2=_TORUS_N2,
        tr_solver=_TORUS_TR_SOLVER,
        max_nfev=_TORUS_MAX_NFEV,
        gauge_weight=_TORUS_GAUGE_WEIGHT,
        rho_weight=_TORUS_RHO_WEIGHT,
    )


def test_callisto_l1_torus_converges_cleanly(
    callisto_l1_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """Live-observed: residual_rms ~5.5e-5, closure ~1.6e-5 -- cleaner than
    JEG's own resonant-orbit-seeded 3:4 torus (~1.6e-4), because this is a
    genuine libration-point-family orbit at physically-realistic (small)
    amplitude, not a distant highly-eccentric resonant orbit."""
    assert callisto_l1_torus.residual_rms < 1e-3
    assert callisto_l1_torus.closure_residual < 1e-2


def test_europa_l2_torus_converges_cleanly(
    europa_l2_torus: vt.CCR4BPTorusVariationalResult,
) -> None:
    """Live-observed: residual_rms ~3.5e-5, closure ~3.1e-5."""
    assert europa_l2_torus.residual_rms < 1e-3
    assert europa_l2_torus.closure_residual < 1e-2


# ---------------------------------------------------------------------------
# Manifold globalization succeeds mechanically on both (new) systems, and
# the Ganymede closest-approach scan: an HONEST, bounded negative within
# this task's own search scope (n_theta2=12, ~90-day propagation window).
# ---------------------------------------------------------------------------

# ~90 days in each system's own nondimensional time unit (TU = base moon's
# own orbital period / 2*pi) -- see module docstring for why this must be
# computed PER SYSTEM (Callisto's and Europa's own orbital periods differ).
_SEARCH_DAYS = 90.0
_CALLISTO_TU_DAYS = _CALLISTO_PERIOD_DAYS / (2.0 * np.pi)
_EUROPA_TU_DAYS = _EUROPA_PERIOD_DAYS / (2.0 * np.pi)


@pytest.fixture(scope="module")
def callisto_unstable_tube(
    callisto_l1_torus: vt.CCR4BPTorusVariationalResult,
) -> mg.ManifoldTube:
    t_max_periods = (_SEARCH_DAYS / _CALLISTO_TU_DAYS) / callisto_l1_torus.period
    return mg.globalize_manifold_tube(
        callisto_l1_torus,
        "unstable",
        n_theta2=12,
        t_max_periods=t_max_periods,
        n_time=100,
        rtol=1e-11,
        atol=1e-11,
    )


@pytest.fixture(scope="module")
def europa_stable_tube(europa_l2_torus: vt.CCR4BPTorusVariationalResult) -> mg.ManifoldTube:
    t_max_periods = (_SEARCH_DAYS / _EUROPA_TU_DAYS) / europa_l2_torus.period
    return mg.globalize_manifold_tube(
        europa_l2_torus,
        "stable",
        n_theta2=12,
        t_max_periods=t_max_periods,
        n_time=100,
        rtol=1e-11,
        atol=1e-11,
    )


def test_callisto_unstable_tube_globalizes_without_failure(
    callisto_unstable_tube: mg.ManifoldTube,
) -> None:
    """The interior-perturber (a_gan < 1) system was never exercised by
    ``globalize_manifold_tube`` before this task -- confirms it globalizes
    mechanically (most departures valid, no all-NaN collapse)."""
    assert callisto_unstable_tube.valid.sum() >= 8  # out of 12


def test_europa_stable_tube_globalizes_without_failure(europa_stable_tube: mg.ManifoldTube) -> None:
    assert europa_stable_tube.valid.sum() >= 8  # out of 12


def test_callisto_l1_unstable_manifold_ganymede_approach_honest_negative(
    callisto_system: ccr4bp.CCR4BPSystem, callisto_unstable_tube: mg.ManifoldTube
) -> None:
    """HONEST NEGATIVE, pinned as a regression floor, not overclaimed: over a
    12-phase, ~90-day search of the Callisto-L1 unstable manifold (seeded
    from the paper's own exact C=3.0044 Lyapunov family), the closest planar
    approach to Ganymede stays far above any genuine flyby scale. The
    paper's own reported ~60-day Callisto-to-Ganymede transit is a single,
    carefully-located trajectory in a much larger search space (their own
    Poincare-map methodology, per the `#710` digest) -- this coarse scan
    does not claim to rule it out, only that it was not found at this
    search scope. Live-observed (this exact test config): dmin ~381,000 km
    at t~44.5 days (the SAME order of magnitude as the paper's own ~60-day
    figure for this leg, but nowhere near flyby-close), r reaching down to
    ~0.74 Callisto-SMA units -- real radial excursion from the L1 seed's own
    ~0.97, just not reaching all the way to Ganymede's own orbital radius
    (a_gan ~0.57)."""
    approach = ct.scan_ganymede_approach(callisto_unstable_tube, callisto_system)
    dmin_km, t_days = ct.to_physical_km_days(approach, CALLISTO_L_KM, _CALLISTO_PERIOD_DAYS)
    assert approach.n_valid_departures >= 8
    # Loose regression bound on the live-observed floor -- documents the
    # search did not find a close encounter, without asserting an exact value.
    assert dmin_km > 100_000.0, (dmin_km, t_days, approach)
    # And the excursion genuinely does move the spacecraft's radius (not a
    # frozen departure-point artifact) -- r_min is comfortably below the L1
    # seed's own ~0.97 radius, confirming the manifold IS exploring, just not
    # reaching all the way to Ganymede's orbit within this search's scope.
    assert approach.r_min < 0.9


def test_europa_l2_stable_manifold_ganymede_approach_honest_negative(
    europa_system: ccr4bp.CCR4BPSystem, europa_stable_tube: mg.ManifoldTube
) -> None:
    """Same honest-negative discipline as the Callisto-side test above,
    applied to the Europa-L2 stable manifold (backward-propagated, per the
    paper's own convention). Live-observed (this exact test config):
    dmin ~232,000 km at t~71.8 days -- a striking near-coincidence with the
    paper's own quoted ~74-day Europa-side transit TIME, but at a distance
    nowhere near flyby-close; not claimed as evidence of the same underlying
    connection (a time-scale match alone is not a connection -- see this
    project's own "it closed!" danger-signal discipline), just noted as an
    interesting, honestly-reported coincidence worth a future closer look."""
    approach = ct.scan_ganymede_approach(europa_stable_tube, europa_system)
    dmin_km, t_days = ct.to_physical_km_days(approach, _EUROPA_L_KM, _EUROPA_PERIOD_DAYS)
    assert approach.n_valid_departures >= 8
    assert dmin_km > 100_000.0, (dmin_km, t_days, approach)
