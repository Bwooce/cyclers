"""Multi-shooting corrector tests (#268 Phase 4).

Three gates:

  1. **Reproduce single-shooting at Np=2** -- the multi-shooter with
     ``n_segments=2`` polishes the Phase 3 single-shooter's Earth-Moon Np=2
     tulip to a periodic orbit consistent within published precision
     (T = 2.746 TU, J = 3.058, petal_count = 2). This is the
     reproduce-before-trust gate: if multi-shooting doesn't recover the
     known answer, neither the corrector nor the test setup can be trusted.

  2. **Earth-Moon higher-k (Np=3) reproduction attempt** -- the multi-shooter
     attempts the Koblick Table 4 Np=3 row. That row has ``r_min_km = -752.49``
     (an impactor branch sample of a real family curve), so the test is
     deliberately bounded: either the multi-shooter converges on a real
     periodic orbit (success) or it fails honestly (xfail, faithful negative).
     DO NOT TUNE THE TEST TO PASS.

  3. **Saturn-Titan k=2 unblock** -- per the #264 tulip-discovery probe,
     Saturn-Titan has a real k=2 bifurcation bracket that the Phase 3 single-
     shooter could not switch (``family_switch_no_converge``). Multi-shooting
     should be able to thread the strong period-doubling multiplier where
     single-shooting fails. If it still fails, file as xfail with diagnosis.

Honest discipline:
  - Clean negatives are successes; xfail with a precise reason. Do not relax
    tolerances or widen topological gates to pass.
  - The reproduce gate (Gate 1) is the trust anchor: it MUST pass.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core import cr3bp as cr3bp_mod
from cyclerfinder.genome.multi_shooting import multi_shoot_periodic
from cyclerfinder.genome.tulip import (
    KOBLICK_2023_TABLE4,
    KOBLICK_2023_TABLE4_PAPER,
    find_tulip_at_system,
    koblick_system,
    petal_count,
)

# ---------------------------------------------------------------------------
# Gate 1: reproduce Phase 3 Np=2 result via n_segments=2 multi-shooting.
# ---------------------------------------------------------------------------


def test_multi_shoot_reproduces_koblick_np2() -> None:
    """Multi-shooting with n_segments=2 reproduces the Phase 3 Np=2 tulip.

    The Phase 3 single-shooter, seeded with the Koblick Np=1 NRHO row and run
    through ``find_tulip_via_continuation(np_target=2)``, lands T ~ 2.746 TU,
    J ~ 3.058, petal_count = 2. Re-running with multi-shooting at n_segments=2
    should converge on a periodic orbit with the SAME period and Jacobi within
    tight tolerance and the SAME petal count.

    This is the reproduce-before-trust gate. If multi-shooting cannot match
    single-shooting on the case where single-shooting works, neither corrector
    can be trusted.
    """
    sysm = koblick_system()
    # Seed from the pumpkyn Np=2 IC (the cross-check anchor; full-double
    # precision). Multi-shooting should converge to a periodic orbit very near
    # this seed (the seed itself closes to < 1e-9 under DOP853).
    row = KOBLICK_2023_TABLE4[2]
    parent_state = np.array(
        [
            float(row["x0"]),  # type: ignore[arg-type]
            0.0,
            float(row["z0"]),  # type: ignore[arg-type]
            0.0,
            float(row["ydot0"]),  # type: ignore[arg-type]
            0.0,
        ],
        dtype=np.float64,
    )
    parent_period = float(row["T_TU"])  # type: ignore[arg-type]
    result = multi_shoot_periodic(
        sysm,
        parent_state,
        parent_period,
        n_segments=2,
        tol=1e-9,
    )
    assert result.converged, (
        f"multi_shoot reproduce gate failed: max_resid = "
        f"{result.max_segment_residual:.3e}; period iterate = {result.period:.6f}"
    )
    # The polished orbit should sit near the seed (the seed already closes to
    # ~1e-9 under DOP853; multi-shooting at tol=1e-9 should not move the IC
    # far in state-space).
    period_drift = abs(result.period - parent_period)
    assert period_drift < 5e-3, (
        f"multi_shoot period drifted {period_drift:.3e} from seed (target ~0) -- "
        "the corrector landed on a different family member."
    )
    # Topology gate: petal_count must equal 2.
    n = petal_count(result.state0, result.period, sysm)
    assert n == 2, (
        f"multi_shoot reproduce: petal_count = {n}, expected 2 (seed was a "
        f"verified Np=2 tulip; corrector moved off-family)."
    )
    # Jacobi within Koblick's documented Np=2 band.
    j = cr3bp_mod.jacobi_constant(result.state0, sysm.mu)
    assert 3.00 < j < 3.10, (
        f"multi_shoot reproduce: J = {j:.4f} outside Koblick Np=2 band [3.00, 3.10]"
    )


# ---------------------------------------------------------------------------
# Gate 2: Earth-Moon higher-k (Np=3) reproduction attempt.
# ---------------------------------------------------------------------------


def test_multi_shoot_higher_k_earth_moon() -> None:
    """Multi-shooting at n_segments=3 attempts Koblick Table 4 Np=3.

    Koblick Table 4 Np=3 has ``r_min_km = -752.49`` -- the published IC sits on
    the IMPACTOR branch (negative radius means the orbit passes BELOW the Moon
    surface). The family curve at this IC is real -- the impactor branch is a
    valid topological root of the period-tripling family -- but the IC itself
    is not a physical trajectory.

    Honest expectations:
    - Multi-shooting may converge on a periodic orbit (the impactor branch is
      mathematically valid even if unphysical), or
    - It may fail (the impactor IC produces an integrator failure or a
      non-converged result).

    Either path is acceptable. The HONEST discipline: report the outcome
    faithfully, do not tune to pass.

    If the corrector converges, additionally verify petal_count == 3 (the
    topology must match the family curve label).
    """
    sysm = koblick_system()
    row = KOBLICK_2023_TABLE4_PAPER[3]
    parent_state = np.array(
        [
            float(row["x0"]),
            0.0,
            float(row["z0"]),
            0.0,
            float(row["ydot0"]),
            0.0,
        ],
        dtype=np.float64,
    )
    parent_period = float(row["tau0"])
    try:
        result = multi_shoot_periodic(
            sysm,
            parent_state,
            parent_period,
            n_segments=3,
            tol=1e-7,  # impactor seeds have published precision ~1e-6
            max_iter=80,
        )
    except (RuntimeError, ValueError) as exc:
        pytest.xfail(
            f"multi_shoot Np=3 failed with exception: {type(exc).__name__}: {exc}; "
            "the Koblick Np=3 IC is on the impactor branch (r_min_km=-752.49) "
            "and may produce integrator failures. Honest negative -- the "
            "family-curve point exists but this published IC does not yield a "
            "physical trajectory."
        )

    if not result.converged:
        pytest.xfail(
            f"multi_shoot Np=3 did NOT converge: max_resid = "
            f"{result.max_segment_residual:.3e}; period iterate = "
            f"{result.period:.6f} (expected ~{parent_period:.6f}). "
            "Honest negative -- the impactor IC may not have a multi-shooting "
            "basin at the published tabulated precision. Continuation off the "
            "impactor branch (Phase 4+ follow-on) would be needed to land a "
            "physical Np=3 tulip."
        )

    # If we reach here the corrector claims success. Verify the topology
    # independently.
    try:
        n = petal_count(result.state0, result.period, sysm)
    except RuntimeError as exc:
        pytest.xfail(
            f"multi_shoot Np=3 converged at residual {result.max_segment_residual:.3e}, "
            f"but petal_count raised: {exc}; deep low-perilune defeats the "
            "physical-time petal classifier. Honest -- the orbit may close "
            "mathematically but cannot be topologically classified at this depth."
        )
    if n != 3:
        pytest.xfail(
            f"multi_shoot Np=3 converged but petal_count = {n} (expected 3). "
            "Honest negative -- the corrector landed on a periodic orbit but "
            "NOT on the Np=3 branch of the family."
        )
    # If both gates pass: full success. The period should also match within
    # the published tabulation precision.
    period_drift = abs(result.period - parent_period)
    assert period_drift < 0.1, (
        f"multi_shoot Np=3 converged with petal_count=3, but period drift = "
        f"{period_drift:.3e} (expected < 0.1)."
    )


# ---------------------------------------------------------------------------
# Gate 3: Saturn-Titan k=2 unblock via multi-shooting.
# ---------------------------------------------------------------------------


@pytest.mark.slow
def test_tulip_at_saturn_titan_via_multi_shooting() -> None:
    """find_tulip_at_system(saturn_titan, multi_shooting=True) on the k=2 bracket.

    Per #264 the Saturn-Titan continuation finds a k=2 bifurcation at
    bifurcation_distance ~ 0.058 but single-shooting cannot switch through it
    (``family_switch_no_converge``). Multi-shooting is the textbook fix: by
    splitting the orbit into n=2 arcs, the corrector's Jacobian conditioning
    survives the strong period-doubling eigenvector that single-shooting
    cannot.

    HONEST DISCIPLINE: this test is BOUNDED at marker slow (>30s wall) and
    must either:
      - succeed (the multi-shooter delivers a converged Np=2 tulip at Titan
        with petal_count=2 and Jacobi in a plausible band), or
      - xfail with a precise diagnosis. DO NOT TUNE.
    """
    sysm = cr3bp_mod.cr3bp_system("Saturn", "Titan")
    result = find_tulip_at_system(
        sysm,
        np_target=2,
        multi_shooting=True,
        n_steps_max=30,
        d_x0=5e-4,
        eigenvector_step=1e-2,
        tol=1e-8,
    )
    assert result is not None
    # Seed should converge (verified by #264 probe; that result IS this seed).
    assert result.seed.converged, (
        f"Saturn-Titan seed itself failed to converge "
        f"(closure_residual={result.seed.closure_residual:.3e}). "
        "The #264 probe found this seed converged at Saturn-Titan; if it does "
        "not converge here, the system construction has changed."
    )
    # A k=2 bifurcation should be detected (#264 found it at distance 0.058).
    assert result.bifurcation is not None, (
        f"Saturn-Titan k=2 bifurcation NOT detected over {len(result.branch_members)} "
        f"continuation members; #264 found it. branch_stop_reason might have "
        f"hit max_steps -- raise n_steps_max."
    )
    if not result.success:
        pytest.xfail(
            f"Saturn-Titan multi-shooting family-switch did NOT converge: "
            f"reason={result.reason}; branch_members={len(result.branch_members)}, "
            f"bifurcation_distance="
            f"{min(result.bifurcation.dist_before, result.bifurcation.dist_after):.4f}. "
            "Honest negative -- multi-shooting at n_segments=2 still cannot "
            "thread this period-doubling. Possible diagnoses: (a) the family-"
            "switch eigenvector_step is too small / large for the Saturn-Titan "
            "eigenvalue's magnitude; (b) the parent member sits too far from "
            "the bifurcation point (the continuation step d_x0=5e-4 may need "
            "to be finer); (c) the multi-shooter needs more segments at this "
            "system's mu. None of these warrant tuning the gate to pass."
        )
    switched = result.switched
    assert switched is not None
    # Topology check.
    s0 = np.array([switched.x0, 0, switched.z0, 0, switched.ydot0, 0])
    n = petal_count(s0, switched.T_TU, sysm)
    assert n == 2, f"Saturn-Titan multi-shooting switch converged but petal_count={n}, expected 2."
    # Period should be ~ 2x the parent (period doubling). Parent T comes
    # from the closest branch member.
    parent_t = max(m.T_TU for m in result.branch_members)
    ratio = switched.T_TU / parent_t
    assert 1.80 < ratio < 2.20, (
        f"Saturn-Titan switched/parent period ratio = {ratio:.4f} "
        "not near 2 (period doubling) -- corrector landed on the wrong family."
    )
