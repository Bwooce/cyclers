"""Tests for the Lion & Handelsman add-an-impulse refinement.

The mechanics gate is a *constructed* case (golden-test discipline): the same
fixed-time "long-way" Lambert transfer that :mod:`tests.verify.test_primer`
flags as IMPROVABLE (interior ``|p| > 1``) must, after one add-an-impulse step,
show a refined total ΔV strictly below the original. This is the qualitative
"teeth" — primer theory predicts an improvement exists and the optimiser must
realise one — with NO invented magnitude asserted.

The Aldrin application is behind ``@pytest.mark.slow`` with qualitative
assertions only (refined < original; both > 0). Its recoverable-ΔV number is
OUR computation (diagnostic, provisional) and is surfaced for the results note,
never asserted as a sourced golden.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.lambert import lambert
from cyclerfinder.verify.primer_refine import (
    AddImpulseRefinement,
    add_impulse_refine_coast,
    refine_coast_from_states,
)

MU_CANONICAL = 1.0


# ---------------------------------------------------------------------------
# Mechanics gate — the IMPROVABLE long-way Lambert transfer refines downward
# ---------------------------------------------------------------------------


def _long_way_refinement(
    r1: float,
    r2: float,
    theta_deg: float,
    tof_factor: float,
    mu: float = MU_CANONICAL,
) -> AddImpulseRefinement:
    """Add-an-impulse refine the long-way (transfer angle > 180°) Lambert coast.

    Mirrors ``tests.verify.test_primer._long_way_transfer``: a circular-orbit
    departure at angle 0 to a circular-orbit arrival at ``theta_deg`` over
    ``tof_factor x t_Hohmann``. The boundary velocities are the local circular
    velocities at each end (so the two-impulse cost matches the diagnostic's
    departure/arrival ΔVs). Such transfers are the classic Lion & Handelsman
    non-optimal case where an interior impulse lowers total ΔV.
    """
    th = np.radians(theta_deg)
    big_r1 = np.array([r1, 0.0, 0.0])
    big_r2 = np.array([r2 * np.cos(th), r2 * np.sin(th), 0.0])
    a_t = 0.5 * (r1 + r2)
    tof = tof_factor * np.pi * np.sqrt(a_t**3 / mu)
    sol = lambert(big_r1, big_r2, tof, mu=mu)[0]
    v_dep = np.asarray(sol.v1, dtype=np.float64)
    v_arr = np.asarray(sol.v2, dtype=np.float64)
    vc1 = np.sqrt(mu / r1) * np.array([0.0, 1.0, 0.0])
    vc2 = np.sqrt(mu / r2) * np.array([-np.sin(th), np.cos(th), 0.0])
    dv1 = v_dep - vc1
    dv2 = vc2 - v_arr
    return refine_coast_from_states(
        big_r1,
        v_dep,
        big_r2,
        tof,
        vc1,
        vc2,
        dv1 / np.linalg.norm(dv1),
        dv2 / np.linalg.norm(dv2),
        mu=mu,
        n_samples=400,
    )


def test_long_way_refinement_strictly_improves() -> None:
    """One add-impulse step on the IMPROVABLE 200° transfer lowers total ΔV.

    Qualitative teeth only: primer theory predicts an improvement exists on a
    coast with interior ``|p| > 1``; the optimiser must realise refined < original
    after a single add-an-impulse step. No magnitude is asserted (it is OUR
    computed value).
    """
    ref = _long_way_refinement(1.0, 4.0, theta_deg=200.0, tof_factor=1.0)
    assert ref.original_dv_kms > 0.0
    assert ref.refined_dv_kms > 0.0
    # Strictly below the original two-impulse cost.
    assert ref.refined_dv_kms < ref.original_dv_kms
    assert ref.recoverable_dv_kms > 0.0
    assert 0.0 < ref.recoverable_fraction < 1.0
    # The midcourse impulse lands in the interior of the coast.
    assert 0.0 < ref.midcourse_time_frac < 1.0


def test_deeper_long_way_refinement_improves() -> None:
    """A 250°, longer transfer (deeper |p| bulge) also refines downward."""
    ref = _long_way_refinement(1.0, 4.0, theta_deg=250.0, tof_factor=1.5)
    assert ref.refined_dv_kms < ref.original_dv_kms
    assert ref.recoverable_dv_kms > 0.0


def test_refinement_never_reports_worse_than_original() -> None:
    """Defensive clamp: refined ΔV is never above the original (zero-offset limit).

    The original single arc is always a feasible point of the refined problem
    (midcourse burn → 0), so the reported refined cost must be ≤ original even on
    a near-optimal coast where the optimiser barely moves.
    """
    # A short, near-Hohmann transfer (angle just under 180°): little to gain.
    ref = _long_way_refinement(1.0, 2.0, theta_deg=170.0, tof_factor=1.0)
    assert ref.refined_dv_kms <= ref.original_dv_kms + 1e-12
    assert ref.recoverable_dv_kms >= -1e-12


def test_add_impulse_rejects_nonpositive_duration() -> None:
    r_a = np.array([1.0, 0.0, 0.0])
    r_b = np.array([0.0, 1.0, 0.0])
    with pytest.raises(ValueError, match="duration must be positive"):
        add_impulse_refine_coast(
            r_a,
            r_b,
            0.0,
            np.array([0.0, 1.0, 0.0]),
            np.array([-1.0, 0.0, 0.0]),
            primer_peak_frac=0.5,
            primer_peak_dir=np.array([1.0, 0.0, 0.0]),
            mu=MU_CANONICAL,
        )


# ---------------------------------------------------------------------------
# Aldrin application (DIAGNOSTIC, PROVISIONAL) — qualitative, slow
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_aldrin_coast0_add_impulse_recoverable_dv() -> None:
    """Add-an-impulse refine coast 0 (E→M) of the real Aldrin maintenance schedule.

    DIAGNOSTIC / PROVISIONAL: the recoverable ΔV is OUR computed value (realised
    improvement from this optimiser), not a sourced or primer-predicted bound.
    Qualitative assertions only — refined ≤ original, both ≥ 0 — and the numbers
    are printed for the results note. The Guzman 2002 multi-rev caveats apply.

    HONEST FINDING (surfaced for the note, not a golden): on the *heliocentric*
    coast-0 two-impulse cost the realised recoverable ΔV is negligible. The
    schedule's 2.9138 km/s lives in the Earth **flyby turn-deficit** (a
    geometric quantity), not in heliocentric midcourse burns — an interior
    heliocentric impulse cannot recover a flyby turn deficit. The primer
    IMPROVABLE flag (max|p| = 1.122) is a necessary-condition violation of the
    *heliocentric* arc; the recoverable heliocentric ΔV it implies is ~0 here.
    """
    from datetime import UTC, datetime

    from cyclerfinder.core.constants import MU_SUN_KM3_S2
    from cyclerfinder.core.ephemeris import Ephemeris
    from cyclerfinder.search.maintain import optimise_aldrin_maintenance_dv

    ephem = Ephemeris("astropy")
    result = optimise_aldrin_maintenance_dv(
        ephem,
        real_window_priority_date=datetime(1985, 10, 28, tzinfo=UTC),
    )
    assert result.converged, "Aldrin maintenance solve did not converge"

    cycler = result.cycler
    # Coast 0 is the Earth->Mars leg: endpoints are encounters 0 and 1.
    enc_a = cycler.encounters[0]
    enc_b = cycler.encounters[1]
    leg0 = cycler.legs[0]
    tof_s = float(enc_b.t - enc_a.t)

    r_a = np.asarray(enc_a.r, dtype=np.float64)
    r_b = np.asarray(enc_b.r, dtype=np.float64)
    v0_depart = np.asarray(leg0.v_depart, dtype=np.float64)

    # Boundary velocities held fixed (like-for-like): pre-departure at A is the
    # spacecraft arrival velocity at the first Earth encounter (v_planet +
    # vinf_in); post-arrival at B is the spacecraft departure velocity for the
    # next leg (v_planet + vinf_out at Mars).
    v_a_before = np.asarray(enc_a.v_planet + enc_a.vinf_in, dtype=np.float64)
    v_b_after = np.asarray(enc_b.v_planet + enc_b.vinf_out, dtype=np.float64)

    # Bounding unit ΔV directions for the primer BVP (departure burn at A,
    # arrival burn at B), with a defensive fallback for a near-ballistic end.
    dv_a = v0_depart - v_a_before
    dv_b = v_b_after - np.asarray(leg0.v_arrive, dtype=np.float64)
    n_a = float(np.linalg.norm(dv_a))
    n_b = float(np.linalg.norm(dv_b))
    p0_hat = dv_a / n_a if n_a > 1e-9 else np.array([1.0, 0.0, 0.0])
    p1_hat = dv_b / n_b if n_b > 1e-9 else np.array([1.0, 0.0, 0.0])

    ref = refine_coast_from_states(
        r_a,
        v0_depart,
        r_b,
        tof_s,
        v_a_before,
        v_b_after,
        p0_hat,
        p1_hat,
        mu=MU_SUN_KM3_S2,
        n_samples=300,
    )

    assert ref.original_dv_kms >= 0.0
    assert ref.refined_dv_kms >= 0.0
    # Refined is never worse than the original (the original arc is a feasible
    # zero-offset point); recovery is >= 0. On this near-ballistic heliocentric
    # coast it is essentially nil (see the docstring finding) — we assert only
    # the direction, never a magnitude.
    assert ref.refined_dv_kms <= ref.original_dv_kms + 1e-12
    assert ref.recoverable_dv_kms >= -1e-12

    print("\nAldrin coast 0 add-an-impulse refinement (DIAGNOSTIC, PROVISIONAL):")
    print(f"  schedule maintenance dv (km/s, our value) = {result.maintenance_dv_kms:.4f}")
    print(f"  original heliocentric coast-0 dv (km/s)   = {ref.original_dv_kms:.6f}")
    print(f"  refined  heliocentric coast-0 dv (km/s)   = {ref.refined_dv_kms:.6f}")
    print(
        f"  realised recoverable (our optimiser)      = {ref.recoverable_dv_kms:.6f} km/s "
        f"({100.0 * ref.recoverable_fraction:.2f}% of the coast-0 cost)"
    )
    print(
        f"  seed primer-peak t/T = {ref.seed_peak_time_frac:.3f}; "
        f"optimised midcourse t/T = {ref.midcourse_time_frac:.3f}"
    )
