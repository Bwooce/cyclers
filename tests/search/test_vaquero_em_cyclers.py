"""Tests for the Vaquero 2013 Earth-Moon 2:1/3:1 cycler reproduction (#799).

Sourced constants are Vaquero's own printed prose values (Sec. 4.4.7,
p. 171; no digit-grade IC table exists for these families -- see the module
docstring). Heavy full-range continuation lives in
``scripts/screen_799_vaquero_em_cycler_families.py`` (results:
``data/found/799_vaquero_em_cycler_families/results.json``); these tests
pin the seed construction, one converged member per family against
Vaquero's own criteria, and a short continuation smoke.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
import cyclerfinder.search.vaquero_em_cyclers as vec


@pytest.fixture(scope="module")
def system() -> cr3bp.CR3BPSystem:
    return vec.earth_moon_system()


@pytest.fixture(scope="module")
def seed21(system: cr3bp.CR3BPSystem) -> tuple[cp.SymmetricOrbit, vec.ApoapsisSeed]:
    return vec.seed_vaquero_family(system, 2, 2.30)


@pytest.fixture(scope="module")
def seed31(system: cr3bp.CR3BPSystem) -> tuple[cp.SymmetricOrbit, vec.ApoapsisSeed]:
    return vec.seed_vaquero_family(system, 3, 2.80)


# ---------------------------------------------------------------------------
# Sourced constants (Vaquero 2013 Sec. 4.4.7, p. 171 prose), verbatim.
# ---------------------------------------------------------------------------


def test_sourced_jacobi_ranges_verbatim() -> None:
    assert vec.VAQUERO_C_RANGE_21 == (1.98, 2.66)
    assert vec.VAQUERO_C_RANGE_31 == (2.54, 3.13)


def test_sourced_endpoint_tofs_verbatim() -> None:
    assert vec.VAQUERO_TOF_DAYS_21 == {1.98: 6.39, 2.66: 4.91}
    assert vec.VAQUERO_TOF_DAYS_31 == {2.54: 4.90, 3.13: 5.04}
    assert vec.VAQUERO_TOF_CEILING_DAYS == 7.0


def test_perigee_band_is_derived_from_printed_altitudes() -> None:
    lo, hi = vec.PERIGEE_BAND_KM
    assert lo == vec.EARTH_RADIUS_KM + 180.0
    assert hi == vec.EARTH_RADIUS_KM + 35786.0


# ---------------------------------------------------------------------------
# Two-body apoapsis seed geometry (pure math).
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("p", [2, 3])
def test_tisserand_seed_semi_major_axis_is_keplerian(p: int) -> None:
    c_mid = 2.3 if p == 2 else 2.8
    g = vec.tisserand_apoapsis_seed(p, c_mid)
    assert g.a == pytest.approx(p ** (-2.0 / 3.0), rel=1e-14)
    assert 0.0 < g.e < 1.0
    assert g.r_apo == pytest.approx(g.a * (1.0 + g.e), rel=1e-14)
    # Rotating-frame apoapsis ydot is negative (sub-synodic angular rate).
    assert g.ydot_rot < 0.0


def test_tisserand_seed_apogees_span_fig444_x0_band() -> None:
    """3:1 apogee inside the Moon's orbit (cislunar), 2:1 beyond it
    (circumlunar) -- the Fig. 4.44 plotted x0 band [~0.6, ~1.2+]."""
    g21 = vec.tisserand_apoapsis_seed(2, 2.30)
    g31 = vec.tisserand_apoapsis_seed(3, 2.80)
    assert 1.0 < g21.r_apo < 1.3
    assert 0.6 < g31.r_apo < 1.0


@pytest.mark.parametrize(
    ("p", "c_range"),
    [(2, vec.VAQUERO_C_RANGE_21), (3, vec.VAQUERO_C_RANGE_31)],
)
def test_tisserand_seed_feasible_across_vaquero_range(p: int, c_range: tuple[float, float]) -> None:
    for c in c_range:
        g = vec.tisserand_apoapsis_seed(p, c)
        assert 0.0 < g.e < 1.0


def test_tisserand_seed_rejects_infeasible() -> None:
    with pytest.raises(ValueError):
        vec.tisserand_apoapsis_seed(3, 5.0)  # beyond the prograde ceiling
    with pytest.raises(ValueError):
        vec.tisserand_apoapsis_seed(1, 2.5)  # not an interior p:1 resonance


# ---------------------------------------------------------------------------
# Seed convergence + Vaquero's own criteria at one mid-range member each.
# ---------------------------------------------------------------------------


def test_seed_2to1_converges_on_family(
    system: cr3bp.CR3BPSystem, seed21: tuple[cp.SymmetricOrbit, vec.ApoapsisSeed]
) -> None:
    orbit, _geom = seed21
    assert orbit.converged
    assert orbit.crossing_residual < 1e-10
    assert orbit.jacobi == pytest.approx(2.30, abs=1e-12)
    # Full period ~ one Moon revolution (2*pi) -- the p:1 signature.
    assert orbit.period / vec.TWO_PI == pytest.approx(1.0, abs=0.05)
    assert vec.half_crossing_index(system, orbit) == vec.VAQUERO_HALF_CROSSINGS


def test_seed_3to1_converges_on_family(
    system: cr3bp.CR3BPSystem, seed31: tuple[cp.SymmetricOrbit, vec.ApoapsisSeed]
) -> None:
    orbit, _geom = seed31
    assert orbit.converged
    assert orbit.crossing_residual < 1e-10
    assert orbit.jacobi == pytest.approx(2.80, abs=1e-12)
    assert orbit.period / vec.TWO_PI == pytest.approx(1.0, abs=0.05)
    assert vec.half_crossing_index(system, orbit) == vec.VAQUERO_HALF_CROSSINGS


def test_member_report_2to1_meets_vaquero_criteria(
    system: cr3bp.CR3BPSystem, seed21: tuple[cp.SymmetricOrbit, vec.ApoapsisSeed]
) -> None:
    orbit, _geom = seed21
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    r = vec.member_report(system, state0, orbit.period, n_samples=50_000)
    assert r.perigee_in_band  # criterion 1 (LEO-GEO insertion band)
    assert r.tof_within_ceiling  # criterion 3 (<= 7 days)
    assert r.moon_min_vs_l2_dist < 1.5  # criterion 2 scale (circumlunar / L2)
    assert r.a_two_body_lu == pytest.approx(2 ** (-2.0 / 3.0), rel=0.05)  # 2:1 resonance
    assert r.stable_planar  # criterion 4: this member is linearly stable
    assert 0.6 < r.x0 < 1.3  # Fig. 4.44 plotted x0 band


def test_member_report_3to1_meets_vaquero_criteria(
    system: cr3bp.CR3BPSystem, seed31: tuple[cp.SymmetricOrbit, vec.ApoapsisSeed]
) -> None:
    orbit, _geom = seed31
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    r = vec.member_report(system, state0, orbit.period, n_samples=50_000)
    assert r.perigee_in_band  # criterion 1
    assert r.tof_within_ceiling  # criterion 3
    assert r.moon_min_vs_l1_dist < 1.5  # criterion 2 scale (cislunar / L1)
    assert r.a_two_body_lu == pytest.approx(3 ** (-2.0 / 3.0), rel=0.05)  # 3:1 resonance
    # Criterion 4, honest form: this member is NOT linearly stable but its
    # unstable mode is small (|lambda| ~ 15 per ~27-day period, vs. 2513.2
    # for the 4:3 family Vaquero herself excludes as a non-cycler).
    assert not r.stable_planar
    assert r.abs_lambda < 50.0
    assert 0.6 < r.x0 < 1.0  # Fig. 4.44 plotted x0 band, cislunar side


def test_seed_members_pass_independent_radau_crosscheck(
    system: cr3bp.CR3BPSystem,
    seed21: tuple[cp.SymmetricOrbit, vec.ApoapsisSeed],
    seed31: tuple[cp.SymmetricOrbit, vec.ApoapsisSeed],
) -> None:
    for orbit, _geom in (seed21, seed31):
        po = cp.PeriodicOrbit(
            state0=np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0]),
            period=orbit.period,
            jacobi=orbit.jacobi,
            converged=True,
            closure_residual=orbit.crossing_residual,
        )
        ok, dj = cp.crosscheck_periodic(system, po, closure_tol=1e-3, jacobi_tol=1e-8)
        assert ok, f"Radau cross-check failed (dJacobi={dj:.2e})"


# ---------------------------------------------------------------------------
# Continuation smoke (short window; the full printed-range walk is the
# screen script's job and is archived in data/found/).
# ---------------------------------------------------------------------------


def test_short_continuation_walks_both_directions(system: cr3bp.CR3BPSystem) -> None:
    rep = vec.reproduce_vaquero_family(system, 3, c_seed=2.80, c_range=(2.78, 2.82), d_jacobi=0.01)
    assert len(rep.members) >= 4
    c_vals = [m.jacobi for m in rep.members]
    assert c_vals == sorted(c_vals)
    assert min(c_vals) == pytest.approx(2.78, abs=1e-6)
    assert max(c_vals) == pytest.approx(2.82, abs=1e-6)
    for m in rep.members:
        # Gauntlet-passed members: converged, Jacobi-conserving, Radau-checked.
        assert m.crossing_residual < 1e-9
        assert m.radau_djacobi < 1e-8
        assert not math.isnan(m.nu)
    # Adjacent members vary smoothly (same family, no topology jump).
    x0s = [m.x0 for m in rep.members]
    assert all(abs(a - b) < 0.01 for a, b in itertools.pairwise(x0s))
