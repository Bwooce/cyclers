"""Mechanics gates + sourced-anchor probe for the one-DSM-per-leg genome.

The mechanics gates (Gate 1-3) are CONSTRUCTED: the "expected" value on every
assertion is defined BY CONSTRUCTION by a known transfer we build here (a Lambert
arc, or a hand-applied impulse), NOT by anything the DSM evaluator itself computes.
These are algorithm-mechanics tests, not golden tests -- exactly the framing the
MBH wrapper's Gate 1 used (``docs/notes/2026-06-07-mbh-wrapper.md``).

The probe (``test_dsm_644gg3_probe``) is the scientific payoff and is marked
``slow`` + wall-capped: it asks whether one-DSM-per-leg + MBH reaches the
sourced-anchor basin (V_inf E=6.44 / M=3.74) that the single-ellipse free-return
genome provably cannot (MBH Gate-3 negative). EXPECTED = the sourced anchors;
emerged values are evidence. Either outcome is reportable.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import propagate
from cyclerfinder.core.lambert import lambert
from cyclerfinder.search.dsm_leg import (
    _flyby_continuity_residual,
    dsm_chain_correct,
    dsm_chain_decision_vector,
    dsm_leg,
    evaluate_dsm_chain,
    sequence_keyed_bounds,
)

MU = MU_SUN_KM3_S2


# ---------------------------------------------------------------------------
# Gate 1 -- eta-degeneracy regression (CONSTRUCTED)
# ---------------------------------------------------------------------------


def test_eta_consistency_zero_dsm_when_v0_is_lambert_velocity() -> None:
    """If v0 IS the Lambert arc's start velocity, dV_DSM ~ 0 for ANY eta.

    Propagating along the Lambert arc for eta*tof and then Lambert-ing the rest
    recovers the SAME arc, so the interior impulse vanishes. This is the clean
    statement of the eta-degeneracy: the one-DSM leg reduces to a pure Lambert leg
    whenever the front arc already follows the Lambert solution. CONSTRUCTED: the
    Lambert arc is the reference, not a DSM-evaluator output.
    """
    r0 = np.array([1.0 * AU_KM, 0.0, 0.0])
    target = np.array([0.0, 1.52 * AU_KM, 0.0])
    tof = 200.0 * 86400.0
    v_lambert = lambert(r0, target, tof, mu=MU)[0].v1

    for eta in (1.0e-4, 0.25, 0.5, 0.75, 0.9999):
        leg = dsm_leg(r0, v_lambert, tof, eta, target, mu=MU)
        assert leg.dv_dsm_kms < 1.0e-6, f"eta={eta}: dV_DSM={leg.dv_dsm_kms}"
        # And the arrival velocity equals the pure-Lambert arrival.
        assert np.linalg.norm(leg.v_arrive - lambert(r0, target, tof, mu=MU)[0].v2) < 1.0e-6


def test_eta_to_zero_reduces_to_start_endpoint_discontinuity() -> None:
    """As eta -> 0 the DSM impulse equals the START-endpoint velocity discontinuity.

    With eta ~ 0 the DSM sits at the leg start, so the leg becomes a pure Lambert
    leg from r0 and the impulse is the mismatch between the supplied v0 and the
    Lambert departure velocity ``||v_lambert - v0||``. CONSTRUCTED: v0 is set to
    the Lambert velocity plus a KNOWN offset; dV must recover that offset.
    """
    r0 = np.array([1.0 * AU_KM, 0.0, 0.0])
    target = np.array([0.0, 1.52 * AU_KM, 0.0])
    tof = 200.0 * 86400.0
    v_lambert = lambert(r0, target, tof, mu=MU)[0].v1

    known_offset = np.array([0.5, -0.3, 0.2])
    v0 = v_lambert + known_offset
    leg = dsm_leg(r0, v0, tof, 1.0e-5, target, mu=MU)
    assert abs(leg.dv_dsm_kms - np.linalg.norm(known_offset)) < 1.0e-3


def test_eta_endpoints_are_rejected_as_singular() -> None:
    """Exact eta in {0, 1} are singular (zero-duration arc) and must raise."""
    r0 = np.array([1.0 * AU_KM, 0.0, 0.0])
    target = np.array([0.0, 1.52 * AU_KM, 0.0])
    tof = 200.0 * 86400.0
    v0 = lambert(r0, target, tof, mu=MU)[0].v1
    for bad in (0.0, 1.0, -0.1, 1.2):
        with pytest.raises(ValueError):
            dsm_leg(r0, v0, tof, bad, target, mu=MU)


# ---------------------------------------------------------------------------
# Gate 2 -- constructed broken-plane two-impulse transfer reproduced exactly
# ---------------------------------------------------------------------------


def test_constructed_broken_plane_transfer_recovered() -> None:
    """Build a known broken-plane DSM transfer by construction; recover it exactly.

    CONSTRUCTION: from an out-of-plane departure ``(r0, v0)`` propagate ballistically
    for ``eta*tof`` to a DSM point, apply a KNOWN out-of-plane impulse ``dv_known``,
    then propagate the post-impulse state forward ``(1-eta)*tof`` to a target. Feed
    ``(r0, v0, tof, eta, r_target)`` back to :func:`dsm_leg`; the Lambert back arc
    must recover the DSM position, the post-impulse departure, the arrival velocity,
    and ``||dv_known||`` -- the broken-plane nature (non-zero z impulse) is what a
    single-ellipse genome cannot represent and this primitive can.
    """
    r0 = np.array([1.0 * AU_KM, 0.0, 0.0])
    v0 = np.array([0.0, 30.0, 1.5])  # out-of-plane -> broken plane
    tof = 300.0 * 86400.0
    eta = 0.4

    r_dsm, v12 = propagate(r0, v0, eta * tof, MU)
    dv_known = np.array([0.8, -0.5, 0.3])
    v21_known = v12 + dv_known
    r_target, v22_known = propagate(r_dsm, v21_known, (1.0 - eta) * tof, MU)

    leg = dsm_leg(r0, v0, tof, eta, r_target, mu=MU)

    assert abs(leg.dv_dsm_kms - np.linalg.norm(dv_known)) < 1.0e-6
    assert np.linalg.norm(leg.r_dsm - r_dsm) < 1.0e-3  # km
    assert np.linalg.norm(leg.v_depart_post_dsm - v21_known) < 1.0e-6
    assert np.linalg.norm(leg.v_arrive - v22_known) < 1.0e-6
    assert leg.v_arrive_pre_dsm.shape == (3,)
    assert abs(leg.t_dsm_sec - eta * tof) < 1.0e-6


# ---------------------------------------------------------------------------
# Gate 3 -- determinism / audit trail (rng-free; full dV breakdown)
# ---------------------------------------------------------------------------


def test_chain_is_deterministic_and_carries_full_audit_trail() -> None:
    """A chained evaluation is rng-free and exposes the full dV breakdown.

    Two identical calls must be byte-identical (no randomness), and the result must
    carry the per-leg DSM impulses, emerged per-body V_inf, eta per leg, and DSM
    states (the audit trail this repo's culture requires).
    """
    ephem = Ephemeris(model="circular")
    sequence = ("E", "M", "E")
    kwargs = dict(
        sequence=sequence,
        t0_sec=0.0,
        vinf_out0_kms=3.0,
        alpha0=0.5,
        beta0=0.0,
        tof_days_per_leg=(200.0, 200.0),
        eta_per_leg=(0.5, 0.5),
        ephem=ephem,
    )
    r1 = evaluate_dsm_chain(**kwargs)  # type: ignore[arg-type]
    r2 = evaluate_dsm_chain(**kwargs)  # type: ignore[arg-type]

    assert r1.total_dv_kms == r2.total_dv_kms
    assert r1.dv_dsm_per_leg_kms == r2.dv_dsm_per_leg_kms

    # Full audit trail present.
    assert len(r1.dv_dsm_per_leg_kms) == 2
    assert set(r1.vinf_in_kms) == {1, 2}
    assert r1.eta_per_leg == (0.5, 0.5)
    assert len(r1.dsm_states) == 2
    # total_dv == sum of per-leg DSM impulses (flyby mission, no rendezvous arrival).
    assert abs(r1.total_dv_kms - sum(r1.dv_dsm_per_leg_kms)) < 1.0e-9
    # max_residual_kms alias matches the MBH-adapter objective.
    assert r1.max_residual_kms == r1.total_dv_kms


def test_sequence_keyed_bounds_layout_and_eta_box() -> None:
    """Takao A.1-A.3 bounds: correct arity, eta box [0,1], alpha/beta boxes."""
    sequence = ("E", "M", "E")
    t0_window = (0.0, 1.0e8)
    b = sequence_keyed_bounds(sequence=sequence, t0_window_sec=t0_window)
    # Layout: [t0, vinf0, alpha0, beta0, *tof(2), *eta(2)] -> length 8.
    assert b.lower.shape == (8,)
    assert b.upper.shape == (8,)
    # alpha in [-pi, pi], beta in [-pi/2, pi/2] (A.1).
    assert b.lower[2] == pytest.approx(-np.pi)
    assert b.upper[2] == pytest.approx(np.pi)
    assert b.lower[3] == pytest.approx(-0.5 * np.pi)
    assert b.upper[3] == pytest.approx(0.5 * np.pi)
    # eta box [0, 1] for both legs (A.1).
    assert list(b.lower[6:8]) == [0.0, 0.0]
    assert list(b.upper[6:8]) == [1.0, 1.0]
    # ToF lower for an inner E<->M leg is the 30-day A.2 floor.
    assert b.lower[4] == pytest.approx(30.0)


# ---------------------------------------------------------------------------
# THE PROBE (slow, wall-capped) -- sourced-anchor gate, 6.44Gg3
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_dsm_644gg3_probe() -> None:
    """Does one-DSM-per-leg + MBH reach the 6.44Gg3 sourced-anchor basin?

    Row ``russell-ch4-6.44Gg3`` (Russell 2004 Table 4.13) is a 3-synodic MULTI-ARC
    cycler: two generic Earth-Earth free-return arcs (g(2.087 yr) + G(4.3191 yr))
    bracketing the E->M outbound (262 d) and M->E return (262 d). SOURCED anchors:
    V_inf E=6.44, M=3.74 km/s, aphelion 1.54 AU. The single-ellipse free-return
    genome PROVABLY cannot reach this (MBH Gate-3: emerged V_inf E=3.01/M=3.06, far
    off-anchor -- ``docs/notes/2026-06-07-mbh-wrapper.md``).

    This probe drives the one-DSM-per-leg chain via the MBH wrapper. EXPECTED = the
    sourced anchors (6.44 / 3.74); the emerged V_inf is the EVIDENCE. We assert only
    that the search RUNS and produces a finite, audited result -- the scientific
    finding (basin reached or not) is reported in the note, and EITHER outcome is a
    legitimate frontier result. No catalogue writeback.
    """
    import time

    from cyclerfinder.search.dsm_leg import make_dsm_chain_step
    from cyclerfinder.search.mbh import mbh

    ephem = Ephemeris(model="circular")
    # Multi-arc topology represented as an explicit E->M->E chain with an interior
    # DSM on each leg (the primitive the single-ellipse genome lacks). The 262-day
    # sourced transits seed both legs; the departure V_inf is seeded at the sourced
    # Earth anchor 6.44 km/s so the search starts in the right energy class.
    sequence = ("E", "M", "E")
    n_legs = 2
    t0_seed = 0.0
    # alpha0 ~ +pi/2: Earth at t=0 (circular backend) sits on +x moving +y, so a
    # tangential prograde departure points the V_inf along +y (azimuth ~ pi/2). The
    # corrector + MBH refine it; this only puts the seed in the right half-plane.
    x0 = dsm_chain_decision_vector(
        t0_sec=t0_seed,
        vinf_out0_kms=6.44,
        alpha0=0.5 * np.pi,
        beta0=0.0,
        tof_days_per_leg=(262.0, 262.0),
        eta_per_leg=(0.5, 0.5),
    )
    bounds = sequence_keyed_bounds(
        sequence=sequence,
        t0_window_sec=(-2.0e7, 2.0e7),
        vinf_out0_bounds_kms=(1.0, 9.0),
    )
    step = make_dsm_chain_step(sequence=sequence, ephem=ephem, bounds=bounds, tol_kms=0.1)

    # Hop every gene that enters the dV objective: t0 (~few days), the departure
    # V_inf magnitude/azimuth/elevation, per-leg ToF (~10 d), per-leg eta (~0.1).
    abs_scale = np.full(x0.shape, np.nan)
    abs_scale[0] = 5.0 * 86400.0  # t0: a few days
    abs_scale[1] = 0.5  # vinf_out0 km/s
    abs_scale[2] = 0.2  # alpha0 rad
    abs_scale[3] = 0.1  # beta0 rad
    abs_scale[4 : 4 + n_legs] = 10.0  # tof days
    abs_scale[4 + n_legs : 4 + 2 * n_legs] = 0.1  # eta

    t_start = time.monotonic()
    result = mbh(
        step,
        x0,
        n_hops=120,
        perturbation="cauchy",
        perturbation_scale=0.0,  # use absolute scales only (frozen relative)
        perturbation_absolute_scale=[float(s) for s in abs_scale],
        rng_seed=6,
        stop_after_stall=60,
    )
    elapsed = time.monotonic() - t_start
    assert elapsed < 600.0, f"probe exceeded 10 min wall cap: {elapsed:.1f}s"

    # The probe must produce a finite, fully-audited result regardless of basin.
    assert np.isfinite(result.best_objective)
    assert result.rng_seed == 6
    assert len(result.objective_history) == result.hops_attempted

    # Report the emerged-vs-sourced evidence (printed for the note; not gated).
    info = result.best_info
    print("\n=== 6.44Gg3 one-DSM-per-leg probe ===")
    print(f"feasible={result.best_feasible}  total_dV={result.best_objective:.4f} km/s")
    print(f"hops attempted/accepted = {result.hops_attempted}/{result.hops_accepted}")
    print(f"emerged V_inf_in  (sourced E=6.44, M=3.74): {info.get('vinf_in_kms')}")
    print(f"emerged V_inf_out: {info.get('vinf_out_kms')}")
    print(f"per-leg dV_DSM km/s: {info.get('dv_dsm_per_leg_kms')}")
    print(f"eta per leg: {info.get('eta_per_leg')}")
    print(f"tof days per leg: {info.get('tof_days_per_leg')}")
    print(f"wall: {elapsed:.1f}s")


# ---------------------------------------------------------------------------
# THE FULL-SEQUENCE PROBE (slow, wall-capped) -- task #153 follow-up to #150.
#
# #150's probe used the MINIMAL E->M->E chain (the control) and recorded a clean
# negative, with the follow-up: drive the FULL multi-arc topology -- the two
# generic Earth-Earth free-return arcs the row's descriptor encodes -- as their
# own legs. This is that probe. Every seeded quantity traces to the row's
# descriptor (``data/catalogue.yaml`` ``russell-ch4-6.44Gg3``):
#
#   * the two generic-arc ToFs from ``free_return_arcs[]`` (first descriptor
#     number = TOF in YEARS, per ``docs/notes/multi-arc-classification.md``):
#     g(2.087 yr) and G(4.3191 yr) -- summing to the 6.41-yr ``period.years``;
#   * the 262-day E->M outbound transit from ``invariants.transit_times_days[0]``
#     (and ``trajectory.segments[out-em]``);
#   * the departure V_inf seeded at the sourced Earth anchor 6.44 km/s
#     (``vinf_kms_at_encounters[0]``), direction-free magnitude.
#
# SOURCED anchors (EXPECTED side): V_inf E=6.44, M=3.74 km/s. The emerged V_inf
# is EVIDENCE -- never imposed (constraint-vs-evidence separation). EITHER
# outcome is a legitimate frontier result; the verdict is reported in
# ``docs/notes/2026-06-07-dsm-full-sequence-probe.md``.
# ---------------------------------------------------------------------------

# Sourced descriptor constants (data/catalogue.yaml russell-ch4-6.44Gg3).
_YEAR_DAYS = 365.25
_G_ARC_YEARS = 2.087  # free_return_arcs[0].tof_years  g(2.087,1111.33,L)
_BIG_G_ARC_YEARS = 4.3191  # free_return_arcs[1].tof_years  G(4.3191,1194.88,L)
_TRANSIT_OUT_DAYS = 262.0  # invariants.transit_times_days[0]; segments[out-em].tof_days
_VINF_E_SOURCED = 6.44  # vinf_kms_at_encounters[0].vinf_kms
_VINF_M_SOURCED = 3.74  # vinf_kms_at_encounters[1].vinf_kms


@pytest.mark.slow
def test_dsm_644gg3_full_sequence_probe() -> None:
    """Full multi-arc sequence for 6.44Gg3: two generic E-E free-return arcs.

    The #150 minimal E->M->E chain floored at 9.4 km/s and predicted the missing
    piece was the FULL multi-leg topology -- the two generic Earth-Earth arcs
    (g 2.087 yr + G 4.3191 yr) the descriptor encodes, each unrolled into its
    Mars-bracketing pieces (E->M outbound + M->E return) so the period sums to the
    sourced 6.41 yr. Sequence ``E-M-E-M-E`` (4 legs):

      leg 1: E->M  262 d (sourced outbound transit)
      leg 2: M->E  (g_arc - 262 d) -- remainder of the 2.087-yr g arc
      leg 3: E->M  262 d (sourced outbound transit)
      leg 4: M->E  (G_arc - 262 d) -- remainder of the 4.3191-yr G arc

    Driven by the MBH wrapper (cauchy, seed 6, <=120 hops, stall 60). EXPECTED =
    the sourced anchors (6.44 / 3.74); the emerged V_inf is EVIDENCE. The probe
    asserts only that the search RUNS and produces a finite, audited result; the
    scientific verdict (basin reached or not) is in the note. No catalogue
    writeback.
    """
    import time

    from cyclerfinder.search.dsm_leg import DsmBounds, make_dsm_chain_step
    from cyclerfinder.search.mbh import mbh

    ephem = Ephemeris(model="circular")
    sequence = ("E", "M", "E", "M", "E")
    n_legs = 4

    g_arc_days = _G_ARC_YEARS * _YEAR_DAYS
    big_g_arc_days = _BIG_G_ARC_YEARS * _YEAR_DAYS
    tof_seed = (
        _TRANSIT_OUT_DAYS,
        g_arc_days - _TRANSIT_OUT_DAYS,
        _TRANSIT_OUT_DAYS,
        big_g_arc_days - _TRANSIT_OUT_DAYS,
    )
    # The two arc ToFs sum to the sourced 6.41-yr period (sanity, not a gate).
    assert abs((g_arc_days + big_g_arc_days) / _YEAR_DAYS - 6.406) < 0.01

    x0 = dsm_chain_decision_vector(
        t0_sec=0.0,
        vinf_out0_kms=_VINF_E_SOURCED,
        alpha0=0.5 * np.pi,  # tangential prograde departure (circular backend, t=0)
        beta0=0.0,
        tof_days_per_leg=tof_seed,
        eta_per_leg=(0.5,) * n_legs,
    )
    # Custom bounds: ToFs +-30% around the sourced descriptor breakdown (the long
    # loop-arc legs do not fit sequence_keyed_bounds' E<->M inner-pair window, and
    # an E->E synodic period is singular -- so the bounds are stated explicitly
    # from the descriptor here).
    lower = np.array(
        [-2.0e7, 1.0, -np.pi, -0.5 * np.pi, *[0.7 * t for t in tof_seed], *([0.0] * n_legs)],
        dtype=np.float64,
    )
    upper = np.array(
        [2.0e7, 9.0, np.pi, 0.5 * np.pi, *[1.3 * t for t in tof_seed], *([1.0] * n_legs)],
        dtype=np.float64,
    )
    bounds = DsmBounds(lower=lower, upper=upper)
    step = make_dsm_chain_step(sequence=sequence, ephem=ephem, bounds=bounds, tol_kms=0.1)

    abs_scale = np.full(x0.shape, np.nan)
    abs_scale[0] = 5.0 * 86400.0  # t0 a few days
    abs_scale[1] = 0.5  # vinf_out0 km/s
    abs_scale[2] = 0.2  # alpha0 rad
    abs_scale[3] = 0.1  # beta0 rad
    abs_scale[4 : 4 + n_legs] = 20.0  # tof days
    abs_scale[4 + n_legs : 4 + 2 * n_legs] = 0.1  # eta

    t_start = time.monotonic()
    result = mbh(
        step,
        x0,
        n_hops=120,
        perturbation="cauchy",
        perturbation_scale=0.0,
        perturbation_absolute_scale=[float(s) for s in abs_scale],
        rng_seed=6,
        stop_after_stall=60,
    )
    elapsed = time.monotonic() - t_start
    assert elapsed < 600.0, f"probe exceeded 10 min wall cap: {elapsed:.1f}s"

    # Finite, fully-audited result regardless of basin (the only gated claim).
    assert np.isfinite(result.best_objective)
    assert result.rng_seed == 6
    assert len(result.objective_history) == result.hops_attempted
    info = result.best_info
    assert len(info["dv_dsm_per_leg_kms"]) == n_legs
    assert set(info["vinf_in_kms"]) == {1, 2, 3, 4}

    print("\n=== 6.44Gg3 FULL-sequence (E-M-E-M-E, two generic arcs) probe ===")
    print(f"feasible={result.best_feasible}  total_dV={result.best_objective:.4f} km/s")
    print(f"hops attempted/accepted = {result.hops_attempted}/{result.hops_accepted}")
    print(
        f"emerged V_inf_in  (sourced E={_VINF_E_SOURCED}, M={_VINF_M_SOURCED}): "
        f"{info['vinf_in_kms']}"
    )
    print(f"per-leg dV_DSM km/s: {info['dv_dsm_per_leg_kms']}")
    print(f"eta per leg: {info['eta_per_leg']}")
    print(f"tof days per leg: {info['tof_days_per_leg']}")
    print(f"wall: {elapsed:.1f}s")


@pytest.mark.slow
def test_dsm_leg_single_rev_floors_on_multirev_loop_arc() -> None:
    """The structural reason the full-sequence probe cannot close: single-rev only.

    The two generic Earth-Earth arcs are MULTI-revolution (a 2.087-yr / 4.319-yr
    flight at a ~1.27-1.54 AU heliocentric ellipse, period ~1.4-1.9 yr, makes 1.5-3
    revolutions). :func:`dsm_leg`'s back-arc Lambert is SINGLE-rev (``max_revs=0``),
    so on a >1-period leg it is forced onto the degenerate near-radial high-energy
    branch. This test demonstrates that the single-rev branch on the G-arc M->E
    return is strictly worse (higher departure speed) than the multi-rev branches
    -- the precise mechanism behind the negative. CONSTRUCTED comparison: the
    reference is the multi-rev Lambert family, not a DSM-evaluator output.
    """
    ephem = Ephemeris(model="circular")
    big_g_arc_days = _BIG_G_ARC_YEARS * _YEAR_DAYS
    t_back = (big_g_arc_days - _TRANSIT_OUT_DAYS) * 86400.0  # G-arc M->E return ToF
    r_mars, _ = ephem.state("M", _TRANSIT_OUT_DAYS * 86400.0)
    r_earth, _ = ephem.state("E", _TRANSIT_OUT_DAYS * 86400.0 + t_back)

    single = lambert(r_mars, r_earth, t_back, mu=MU, prograde=True, max_revs=0)
    multi = lambert(r_mars, r_earth, t_back, mu=MU, prograde=True, max_revs=3)

    v1_single = min(float(np.linalg.norm(s.v1)) for s in single)
    v1_multi = min(float(np.linalg.norm(s.v1)) for s in multi)
    # The leg is genuinely multi-rev: more solution branches exist beyond single-rev.
    assert len(multi) > len(single)
    # The single-rev branch dsm_leg is locked to is strictly worse (higher |v1|)
    # than the best multi-rev branch -- the loop arc cannot be represented well by
    # the single-rev primitive.
    assert v1_multi < v1_single - 1.0


# ---------------------------------------------------------------------------
# Multi-rev mechanics gate (#157, CONSTRUCTED) -- the back-arc multi-rev fix.
#
# #153 proved the LAMBERT FAMILY (lambert(..., max_revs=3)) has a ~15.5 km/s
# best-|v1| branch vs the ~28.2 km/s single-rev branch on the 1315.5-d G-arc
# M->E return. This gate proves the same improvement is now wired into the
# dsm_leg PRIMITIVE: dsm_leg(..., max_revs=3) selects the multi-rev branch and
# its DSM impulse / departure speed drop substantially below the single-rev
# (max_revs=0) result on the very same arc. CONSTRUCTED comparison: the
# reference is the leg's own single-rev result and the multi-rev Lambert family,
# not any DSM-evaluator/optimiser output (golden-rule mechanics, not a golden
# test). Mechanics label.
# ---------------------------------------------------------------------------


def test_dsm_leg_max_revs_recovers_multirev_improvement_on_g_arc() -> None:
    """dsm_leg(max_revs=3) reproduces the #153 multi-rev |v1| improvement.

    On the 1315.5-d G-arc M->E return (the loop arc that floors the single-rev
    full sequence), a near-pure-Lambert leg (eta -> small, so the DSM sits at the
    leg start and the back arc is essentially the whole M->E transfer):

      * single-rev (max_revs=0) is forced onto the degenerate near-radial branch
        -- departure speed |v1| ~ 28 km/s (cf. #153's measured 28.18 km/s);
      * multi-rev (max_revs=3) selects a multi-revolution branch with |v1| ~ 15
        km/s (cf. #153's measured 15.51 km/s),

    and the selected DSM impulse drops strictly and substantially. This is the
    mechanism the full-sequence re-probe needs; the assertion thresholds are loose
    floors around the #153-cited numbers, not exact golden values.
    """
    ephem = Ephemeris(model="circular")
    t_back = (_BIG_G_ARC_YEARS * _YEAR_DAYS - _TRANSIT_OUT_DAYS) * 86400.0
    r_mars, v_mars = ephem.state("M", _TRANSIT_OUT_DAYS * 86400.0)
    r_earth, _ = ephem.state("E", _TRANSIT_OUT_DAYS * 86400.0 + t_back)

    # eta -> small: the front (ballistic) arc is negligible, so the back-arc
    # Lambert spans (almost) the whole G-arc M->E return -- the multi-rev regime.
    eta = 0.01
    v0 = np.asarray(v_mars) + np.array([0.0, 3.0, 0.0])  # a modest departure V_inf

    leg_single = dsm_leg(r_mars, v0, t_back, eta, np.asarray(r_earth), mu=MU, max_revs=0)
    leg_multi = dsm_leg(r_mars, v0, t_back, eta, np.asarray(r_earth), mu=MU, max_revs=3)

    v1_single = float(np.linalg.norm(leg_single.v_depart_post_dsm))
    v1_multi = float(np.linalg.norm(leg_multi.v_depart_post_dsm))

    # The single-rev path is the legacy degenerate branch (#153: ~28.2 km/s).
    assert leg_single.n_revs_chosen == 0
    assert leg_single.branch_chosen == "single"
    assert v1_single > 25.0  # near the #153-cited 28.18 km/s
    # The multi-rev path selects a genuine multi-revolution branch (#153: ~15.5).
    assert leg_multi.n_revs_chosen >= 1
    assert v1_multi < 18.0  # near the #153-cited 15.51 km/s
    # The improvement is strict and substantial (>5 km/s on |v1|, and the DSM
    # impulse itself drops well below the single-rev value).
    assert v1_multi < v1_single - 5.0
    assert leg_multi.dv_dsm_kms < leg_single.dv_dsm_kms - 3.0


def test_dsm_leg_max_revs_zero_is_bit_identical_to_default() -> None:
    """max_revs=0 (the default) is byte-identical to the historical single-rev path.

    Regression guard for the additive multi-rev change: passing max_revs=0
    explicitly, omitting it, and the converged single-rev branch must all produce
    the SAME DSM impulse, arrival velocity, and audit (n_revs=0, branch="single").
    """
    r0 = np.array([1.0 * AU_KM, 0.0, 0.0])
    target = np.array([0.0, 1.52 * AU_KM, 0.0])
    tof = 200.0 * 86400.0
    v0 = lambert(r0, target, tof, mu=MU)[0].v1 + np.array([0.4, -0.2, 0.1])

    default = dsm_leg(r0, v0, tof, 0.5, target, mu=MU)
    explicit_zero = dsm_leg(r0, v0, tof, 0.5, target, mu=MU, max_revs=0)

    assert default.dv_dsm_kms == explicit_zero.dv_dsm_kms
    assert np.array_equal(default.v_arrive, explicit_zero.v_arrive)
    assert default.n_revs_chosen == 0
    assert default.branch_chosen == "single"


def test_dsm_leg_rev_branch_explicit_selection() -> None:
    """An explicit rev_branch selector picks exactly that branch (or raises)."""
    ephem = Ephemeris(model="circular")
    t_back = (_BIG_G_ARC_YEARS * _YEAR_DAYS - _TRANSIT_OUT_DAYS) * 86400.0
    r_mars, v_mars = ephem.state("M", _TRANSIT_OUT_DAYS * 86400.0)
    r_earth, _ = ephem.state("E", _TRANSIT_OUT_DAYS * 86400.0 + t_back)
    v0 = np.asarray(v_mars) + np.array([0.0, 3.0, 0.0])

    leg = dsm_leg(r_mars, v0, t_back, 0.01, np.asarray(r_earth), mu=MU, rev_branch=(1, "low"))
    assert leg.n_revs_chosen == 1
    assert leg.branch_chosen == "low"

    # An impossible revolution count for this ToF raises LambertError.
    from cyclerfinder.core.lambert import LambertError

    with pytest.raises(LambertError):
        dsm_leg(r_mars, v0, t_back, 0.01, np.asarray(r_earth), mu=MU, rev_branch=(99, "low"))


# ---------------------------------------------------------------------------
# THE MULTI-REV RE-PROBES (#157, slow, wall-capped) -- the #153 probes re-run
# with the back-arc multi-rev branch now wired in (max_revs sized to the arc
# periods). EXPECTED = each row's OWN sourced anchors; emerged V_inf is EVIDENCE.
# Either outcome is reportable; NO catalogue writeback. Verdicts recorded in
# docs/notes/2026-06-07-dsm-multirev-probe.md.
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_dsm_644gg3_full_sequence_multirev_reprobe() -> None:
    """6.44Gg3 full E-M-E-M-E sequence re-run with back-arc multi-rev (#157).

    Identical seeding to ``test_dsm_644gg3_full_sequence_probe`` (the #153
    single-rev baseline that floored at ~30 km/s) but with ``max_revs=3`` so the
    long loop-arc legs (g 2.087 yr ~1.5 rev, G 4.319 yr ~2-3 rev) may take a
    multi-revolution branch -- the first configuration in which the multi-arc
    topology is representable. EXPECTED = the sourced anchors (E=6.44 / M=3.74);
    the emerged V_inf is EVIDENCE. Asserts only a finite, audited result; the
    verdict is in the note.
    """
    import time

    from cyclerfinder.search.dsm_leg import DsmBounds, make_dsm_chain_step
    from cyclerfinder.search.mbh import mbh

    ephem = Ephemeris(model="circular")
    sequence = ("E", "M", "E", "M", "E")
    n_legs = 4

    g_arc_days = _G_ARC_YEARS * _YEAR_DAYS
    big_g_arc_days = _BIG_G_ARC_YEARS * _YEAR_DAYS
    tof_seed = (
        _TRANSIT_OUT_DAYS,
        g_arc_days - _TRANSIT_OUT_DAYS,
        _TRANSIT_OUT_DAYS,
        big_g_arc_days - _TRANSIT_OUT_DAYS,
    )

    x0 = dsm_chain_decision_vector(
        t0_sec=0.0,
        vinf_out0_kms=_VINF_E_SOURCED,
        alpha0=0.5 * np.pi,
        beta0=0.0,
        tof_days_per_leg=tof_seed,
        eta_per_leg=(0.5,) * n_legs,
    )
    lower = np.array(
        [-2.0e7, 1.0, -np.pi, -0.5 * np.pi, *[0.7 * t for t in tof_seed], *([0.0] * n_legs)],
        dtype=np.float64,
    )
    upper = np.array(
        [2.0e7, 9.0, np.pi, 0.5 * np.pi, *[1.3 * t for t in tof_seed], *([1.0] * n_legs)],
        dtype=np.float64,
    )
    bounds = DsmBounds(lower=lower, upper=upper)
    # max_revs=3: g arc ~1.5 rev, G arc ~2-3 rev (#153 diagnosis).
    step = make_dsm_chain_step(
        sequence=sequence, ephem=ephem, bounds=bounds, tol_kms=0.1, max_revs=3
    )

    abs_scale = np.full(x0.shape, np.nan)
    abs_scale[0] = 5.0 * 86400.0
    abs_scale[1] = 0.5
    abs_scale[2] = 0.2
    abs_scale[3] = 0.1
    abs_scale[4 : 4 + n_legs] = 20.0
    abs_scale[4 + n_legs : 4 + 2 * n_legs] = 0.1

    t_start = time.monotonic()
    result = mbh(
        step,
        x0,
        n_hops=120,
        perturbation="cauchy",
        perturbation_scale=0.0,
        perturbation_absolute_scale=[float(s) for s in abs_scale],
        rng_seed=6,
        stop_after_stall=60,
    )
    elapsed = time.monotonic() - t_start
    assert elapsed < 600.0, f"probe exceeded 10 min wall cap: {elapsed:.1f}s"

    assert np.isfinite(result.best_objective)
    assert result.rng_seed == 6
    assert len(result.objective_history) == result.hops_attempted
    info = result.best_info
    assert len(info["dv_dsm_per_leg_kms"]) == n_legs
    assert set(info["vinf_in_kms"]) == {1, 2, 3, 4}

    print("\n=== 6.44Gg3 FULL-sequence MULTI-REV re-probe (max_revs=3) ===")
    print(f"feasible={result.best_feasible}  total_dV={result.best_objective:.4f} km/s")
    print(f"hops attempted/accepted = {result.hops_attempted}/{result.hops_accepted}")
    print(
        f"emerged V_inf_in  (sourced E={_VINF_E_SOURCED}, M={_VINF_M_SOURCED}): "
        f"{info['vinf_in_kms']}"
    )
    print(f"per-leg dV_DSM km/s: {info['dv_dsm_per_leg_kms']}")
    print(f"n_revs per leg: {info.get('n_revs_per_leg')}")
    print(f"branch per leg: {info.get('branch_per_leg')}")
    print(f"eta per leg: {info['eta_per_leg']}")
    print(f"tof days per leg: {info['tof_days_per_leg']}")
    print(f"wall: {elapsed:.1f}s")


# Sourced descriptor constants for russell-ch4-4.991gG2 (the S1L1 Russell
# free-return framing). data/catalogue.yaml russell-ch4-4.991gG2:
#   free_return_arcs: g(1.4612,...) + G(2.8096,...) ; transit out/in = 150/150 d;
#   vinf_kms_at_encounters: E=4.99, M=5.10 ; aphelion 1.64 AU ; period 4.27 yr.
# These are this ROW's OWN anchors -- NOT the s1l1-2syn-em-cpom 5.65/3.05 framing
# (a different idealisation of the same physical cycler; multi-arc-classification
# §7/§12). Do not mix framings.
_S1L1_G_ARC_YEARS = 1.4612  # russell-ch4-4.991gG2 free_return_arcs[0].tof_years
_S1L1_BIG_G_ARC_YEARS = 2.8096  # free_return_arcs[1].tof_years
_S1L1_TRANSIT_DAYS = 150.0  # trajectory.segments[out-em].tof_days
_S1L1_VINF_E_SOURCED = 4.99  # vinf_kms_at_encounters[0].vinf_kms
_S1L1_VINF_M_SOURCED = 5.10  # vinf_kms_at_encounters[1].vinf_kms


@pytest.mark.slow
def test_dsm_s1l1_4991gg2_two_arc_multirev_reprobe() -> None:
    """S1L1 (russell-ch4-4.991gG2) two-arc geometry, back-arc multi-rev (#157).

    The 4.991gG2 row is the Russell free-return framing of the S1L1 / "Notable"
    2-synodic Earth-Mars cycler: two generic arcs g(1.4612 yr) + G(2.8096 yr)
    bracketing the 150-day E->M outbound transit (this row's OWN sourced anchors:
    V_inf E=4.99, M=5.10 km/s, aphelion 1.64 AU -- NOT the s1l1-2syn-em-cpom
    5.65/3.05 framing; multi-arc-classification §7/§12). The arcs sum to the row's
    4.27-yr period.

    Sequence ``E-M-E-M-E`` (4 legs), the two arcs unrolled into their
    Mars-bracketing pieces, with ``max_revs=3`` so the loop arcs may take a
    multi-revolution branch. EXPECTED = the row's sourced anchors (4.99 / 5.10);
    emerged V_inf is EVIDENCE. Asserts only a finite, audited result; the verdict
    (a close-and-match would be the first multi-arc closure) is in the note. No
    catalogue writeback.
    """
    import time

    from cyclerfinder.search.dsm_leg import DsmBounds, make_dsm_chain_step
    from cyclerfinder.search.mbh import mbh

    ephem = Ephemeris(model="circular")
    sequence = ("E", "M", "E", "M", "E")
    n_legs = 4

    g_arc_days = _S1L1_G_ARC_YEARS * _YEAR_DAYS
    big_g_arc_days = _S1L1_BIG_G_ARC_YEARS * _YEAR_DAYS
    tof_seed = (
        _S1L1_TRANSIT_DAYS,
        g_arc_days - _S1L1_TRANSIT_DAYS,
        _S1L1_TRANSIT_DAYS,
        big_g_arc_days - _S1L1_TRANSIT_DAYS,
    )
    # The two arc ToFs sum to the sourced 4.27-yr period (sanity, not a gate).
    assert abs((g_arc_days + big_g_arc_days) / _YEAR_DAYS - 4.2708) < 0.01

    x0 = dsm_chain_decision_vector(
        t0_sec=0.0,
        vinf_out0_kms=_S1L1_VINF_E_SOURCED,
        alpha0=0.5 * np.pi,
        beta0=0.0,
        tof_days_per_leg=tof_seed,
        eta_per_leg=(0.5,) * n_legs,
    )
    lower = np.array(
        [-2.0e7, 1.0, -np.pi, -0.5 * np.pi, *[0.7 * t for t in tof_seed], *([0.0] * n_legs)],
        dtype=np.float64,
    )
    upper = np.array(
        [2.0e7, 9.0, np.pi, 0.5 * np.pi, *[1.3 * t for t in tof_seed], *([1.0] * n_legs)],
        dtype=np.float64,
    )
    bounds = DsmBounds(lower=lower, upper=upper)
    step = make_dsm_chain_step(
        sequence=sequence, ephem=ephem, bounds=bounds, tol_kms=0.1, max_revs=3
    )

    abs_scale = np.full(x0.shape, np.nan)
    abs_scale[0] = 5.0 * 86400.0
    abs_scale[1] = 0.5
    abs_scale[2] = 0.2
    abs_scale[3] = 0.1
    abs_scale[4 : 4 + n_legs] = 20.0
    abs_scale[4 + n_legs : 4 + 2 * n_legs] = 0.1

    t_start = time.monotonic()
    result = mbh(
        step,
        x0,
        n_hops=120,
        perturbation="cauchy",
        perturbation_scale=0.0,
        perturbation_absolute_scale=[float(s) for s in abs_scale],
        rng_seed=6,
        stop_after_stall=60,
    )
    elapsed = time.monotonic() - t_start
    assert elapsed < 600.0, f"probe exceeded 10 min wall cap: {elapsed:.1f}s"

    assert np.isfinite(result.best_objective)
    assert result.rng_seed == 6
    assert len(result.objective_history) == result.hops_attempted
    info = result.best_info
    assert len(info["dv_dsm_per_leg_kms"]) == n_legs
    assert set(info["vinf_in_kms"]) == {1, 2, 3, 4}

    print("\n=== S1L1 4.991gG2 two-arc MULTI-REV re-probe (max_revs=3) ===")
    print(f"feasible={result.best_feasible}  total_dV={result.best_objective:.4f} km/s")
    print(f"hops attempted/accepted = {result.hops_attempted}/{result.hops_accepted}")
    print(
        f"emerged V_inf_in  (sourced E={_S1L1_VINF_E_SOURCED}, M={_S1L1_VINF_M_SOURCED}): "
        f"{info['vinf_in_kms']}"
    )
    print(f"per-leg dV_DSM km/s: {info['dv_dsm_per_leg_kms']}")
    print(f"n_revs per leg: {info.get('n_revs_per_leg')}")
    print(f"branch per leg: {info.get('branch_per_leg')}")
    print(f"eta per leg: {info['eta_per_leg']}")
    print(f"tof days per leg: {info['tof_days_per_leg']}")
    print(f"wall: {elapsed:.1f}s")


# ---------------------------------------------------------------------------
# Phase 1 (#162) -- mechanics gates for the explicit flyby V_inf-continuity +
# bend-feasibility residual. CONSTRUCTED: every "expected" value is a hand-built
# feasible / infeasible flyby geometry (or a hand-computed flyby_dv), never a
# DSM-evaluator output. These prove the new per-flyby residual term is correct
# BEFORE the decisive 6.44Gg3 probe consumes it.
# ---------------------------------------------------------------------------

# Mars flyby constants (cyclerfinder.core.constants PLANETS["M"]).
_MU_MARS = 4.282837521e4  # km^3/s^2
_RP_MIN_MARS = 3396.19 + 300.0  # radius_eq_km + safe_alt_km, km


def _rotate_in_plane(vec: np.ndarray, angle_rad: float) -> np.ndarray:
    """Rotate a 3-vector by ``angle_rad`` about the z axis (in the xy plane)."""
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    return np.array(
        [c * vec[0] - s * vec[1], s * vec[0] + c * vec[1], vec[2]],
        dtype=np.float64,
    )


def test_flyby_continuity_residual_zero_on_feasible_flyby() -> None:
    """The per-flyby residual term is ~0 on a CONSTRUCTED bend-feasible flyby.

    Build a v_inf^- / v_inf^+ pair with EQUAL magnitude and a turn angle strictly
    inside ``max_bend(mu_M, rp_min, |v_inf|)``. A ballistic Mars flyby can realise
    this exactly, so ``flyby_dv == 0`` and the continuity residual term must vanish.
    Reference = the known feasible geometry (the bend cone), NOT the evaluator.
    Label: mechanics.
    """
    from cyclerfinder.core.flyby import max_bend

    vinf = 3.74  # the sourced Mars anchor magnitude (a physically relevant speed)
    vin = np.array([vinf, 0.0, 0.0])
    delta_max = max_bend(_MU_MARS, _RP_MIN_MARS, vinf)
    # A turn strictly inside the cone (half of the achievable maximum).
    vout = _rotate_in_plane(vin, 0.5 * delta_max)
    assert abs(np.linalg.norm(vout) - vinf) < 1.0e-12  # equal magnitude by construction

    r = _flyby_continuity_residual(vin, vout, "M")
    assert r < 1.0e-9, f"feasible flyby must charge ~0, got {r}"


def test_flyby_continuity_residual_charges_infeasible_turn() -> None:
    """The per-flyby residual equals the hand-computed flyby_dv on an infeasible turn.

    Build an equal-magnitude pair whose turn angle is BEYOND ``max_bend`` (a
    near-reversal at high V_inf, where the bend cone is narrow). The residual must
    equal the independently hand-computed ``flyby_dv`` (> 0). This is precisely the
    bend-feasibility term the old scalar objective could not see. Label: mechanics.
    """
    from cyclerfinder.core.flyby import flyby_dv, max_bend

    vinf = 8.0  # high V_inf -> narrow cone
    vin = np.array([vinf, 0.0, 0.0])
    delta_max = max_bend(_MU_MARS, _RP_MIN_MARS, vinf)
    # A turn well beyond the cone (160 deg, far past delta_max which is small here).
    turn = np.radians(160.0)
    assert turn > delta_max  # infeasible by construction
    vout = _rotate_in_plane(vin, turn)
    assert abs(np.linalg.norm(vout) - vinf) < 1.0e-12

    expected = flyby_dv(vin, vout, _MU_MARS, _RP_MIN_MARS)
    assert expected > 0.0  # constructed-infeasible -> strictly positive
    r = _flyby_continuity_residual(vin, vout, "M")
    assert abs(r - expected) < 1.0e-6, f"residual {r} != hand-computed flyby_dv {expected}"


def test_charge_flyby_continuity_default_off_is_bit_identical() -> None:
    """charge_flyby_continuity=False (default) is bit-identical to current code.

    The regression that guarantees no existing caller/test changes: with the flag
    off, ``evaluate_dsm_chain`` and ``dsm_chain_correct`` return exactly the current
    scalar-objective results on the E-M-E mechanics fixture, and the residual is
    still the length-1 scalar. Label: regression.
    """
    ephem = Ephemeris(model="circular")
    sequence = ("E", "M", "E")
    kwargs = dict(
        sequence=sequence,
        t0_sec=0.0,
        vinf_out0_kms=3.0,
        alpha0=0.5,
        beta0=0.0,
        tof_days_per_leg=(200.0, 200.0),
        eta_per_leg=(0.5, 0.5),
        ephem=ephem,
    )
    base = evaluate_dsm_chain(**kwargs)  # type: ignore[arg-type]
    flagged = evaluate_dsm_chain(charge_flyby_continuity=False, **kwargs)  # type: ignore[arg-type]

    # Bit-identical scalar objective and per-leg breakdown.
    assert flagged.total_dv_kms == base.total_dv_kms
    assert flagged.dv_dsm_per_leg_kms == base.dv_dsm_per_leg_kms
    assert flagged.vinf_in_kms == base.vinf_in_kms
    assert flagged.vinf_out_kms == base.vinf_out_kms
    # The default residual_vector is the length-1 scalar [total_dv].
    assert flagged.residual_vector.shape == (1,)
    assert flagged.residual_vector[0] == base.total_dv_kms

    # dsm_chain_correct default path is unchanged too (same converged decision).
    x0 = dsm_chain_decision_vector(
        t0_sec=0.0,
        vinf_out0_kms=3.0,
        alpha0=0.5,
        beta0=0.0,
        tof_days_per_leg=(200.0, 200.0),
        eta_per_leg=(0.5, 0.5),
    )
    c_base = dsm_chain_correct(x0, sequence=sequence, ephem=ephem, max_nfev=20)
    c_flag = dsm_chain_correct(
        x0, sequence=sequence, ephem=ephem, max_nfev=20, charge_flyby_continuity=False
    )
    assert c_flag.total_dv_kms == c_base.total_dv_kms
    assert c_flag.tof_days_per_leg == c_base.tof_days_per_leg
    assert c_flag.eta_per_leg == c_base.eta_per_leg
