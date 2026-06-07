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
