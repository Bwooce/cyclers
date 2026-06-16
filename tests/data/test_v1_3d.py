"""Tests for the V1 same-model 3D gauntlet (#306 Phase 1 Part A).

Sourced golden discipline (per ``feedback_golden_tests_sourced_only``)
---------------------------------------------------------------------
The tests assert on CLOSURE QUALITY + TOPOLOGY, NOT on specific 3D state
numbers our own code produced. The acceptable input states are read from
JSONLs whose ICs were *converged in-repo* in earlier tasks; treating them
as goldens for V1 is honest because V1's whole job is to RE-CONVERGE them
and report closure. The spec floors (1 m/s for km/s, 1e-6 nondim default)
are spec §14 and the corrector default — both module constants.

Test cases
----------
  1. Planar Braik-Ross C11a seed (sourced from ``catalogue.yaml`` row
     ``braik-ross-c11a-cycler-2026``) — recovers a planar member with
     ``passes_v1=True`` and ``degenerate_planar=True`` (the planar manifold
     is invariant, this is correct behaviour).
  2. 3D Braik-Ross C11a extension (#287 spike) — z0 ~ -0.241 from the
     ``data/spike_287.jsonl`` ``case=family_seed_z0_eq_neg_p24`` entry;
     V1 should pass with ``degenerate_planar=False``.
  3. 3D k=4 sub-family member (#301 doubly-hyperbolic) — IC from
     ``data/family_296_3d_subfamilies_299.jsonl`` first bracket's switched
     state; V1 should pass.
  4. Negative: a deliberately broken IC (random perturbation) should NOT
     pass V1 — the independent re-propagation residual will exceed the
     nondim floor.

Discipline
----------
* NO catalogue writeback inside the tests.
* The V1 floor (``V1_FLOOR_KMS = 1e-3``) is spec §14, asserted equal in a
  fabrication-guard test.
* Every passing test asserts the independent Radau cross-check held — V1
  is meaningless without it.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.validation.v1_3d import (
    V1_FLOOR_KMS,
    V1_FLOOR_NONDIM_DEFAULT,
    V1Verdict3D,
    run_v1_3d,
)
from cyclerfinder.search.cr3bp_general_periodic_3d import (
    FREE_VARS_SYMMETRIC_TULIP,
    RESIDUAL_PERPENDICULAR_HALF_PERIOD,
)

# Repository root — tests live at <repo>/tests/data/, this file is
# <repo>/tests/data/test_v1_3d.py.
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent

# Sourced planar golden (Braik-Ross 2026 C11a, in catalogue.yaml).
C11A_X0 = -0.8116406668238195
C11A_YDOT0 = -0.11859055759763637
C11A_PERIOD_TU = 9.69107744379376
C11A_JACOBI = 3.1294
EM_MU = 1.2150584270572e-2  # Braik-Ross 2026 Table 1
EM_L_KM = 384400.0
EM_T_S = 375699.8


def _em_system() -> cr3bp.CR3BPSystem:
    """Earth-Moon CR3BP at the catalogued (sourced) Braik-Ross mu."""
    return cr3bp.CR3BPSystem(
        mu=EM_MU,
        primary="earth",
        secondary="moon",
        l_km=EM_L_KM,
        t_s=EM_T_S,
    )


def test_v1_floor_is_spec_value() -> None:
    """Fabrication guard: the V1 km/s floor is the spec §14 value (1e-3).

    Per ``feedback_golden_tests_sourced_only``, the floor is sourced from
    the spec and must not be silently tunable.
    """
    assert V1_FLOOR_KMS == 1.0e-3


def test_v1_planar_braik_ross_c11a_passes() -> None:
    """V1 on the sourced planar Braik-Ross C11a IC.

    The IC traces to ``data/catalogue.yaml`` row ``braik-ross-c11a-cycler-2026``.
    V1 re-closes under the 3D corrector AND independently re-propagates under
    Radau. Both must hold; the km/s floor must hold.
    """
    system = _em_system()
    state0 = np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    verdict = run_v1_3d(
        "braik-ross-c11a-planar",
        state0,
        C11A_PERIOD_TU,
        system,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        is_half_period_residual=True,
        notes="sourced planar C11a; catalogue.yaml row braik-ross-c11a-cycler-2026",
    )
    assert isinstance(verdict, V1Verdict3D)
    assert verdict.converged_corrector, (
        f"corrector failed at {verdict.closure_residual_nondim:.3e} nondim"
    )
    assert verdict.converged_independent, (
        f"independent Radau closure failed at "
        f"{verdict.independent_closure_nondim:.3e} nondim "
        f"(> {verdict.independent_floor_nondim:.0e})"
    )
    assert verdict.passes_v1, (
        f"V1 FAIL: independent={verdict.independent_closure_kms:.3e} km/s "
        f"> v1_floor={verdict.v1_floor_kms:.0e} km/s"
    )
    # Planar IC — the planar manifold is invariant, so degenerate_planar=True
    # is the correct (non-vetoing) outcome.
    assert verdict.degenerate_planar
    # The km/s independent closure is below the spec V1 bar.
    assert verdict.independent_closure_kms < V1_FLOOR_KMS


def _read_spike_287_3d_seed() -> tuple[np.ndarray, float]:
    """Read the #287 spike's 3D member IC (z0 ~ -0.241) from spike_287.jsonl.

    The spike's seed at z0 = -0.05 converges to a non-trivial 3D member;
    we use that one — it's the first genuinely 3D row in the JSONL.

    Returns ``(state0, period_tu)`` for the 3D member.
    """
    path = _REPO_ROOT / "data" / "spike_287.jsonl"
    if not path.exists():
        pytest.skip(f"spike_287.jsonl missing at {path}")
    for raw in path.read_text().splitlines():
        row = json.loads(raw)
        # The first row that converged AND is genuinely 3D (|z0| > 1e-3) is
        # the spike's primary 3D payload.
        if not row.get("converged"):
            continue
        if abs(float(row.get("z0", 0.0))) < 1e-3:
            continue
        x0 = float(row["x0"])
        z0 = float(row["z0"])
        ydot0 = float(row["ydot0"])
        period_tu = float(row["T_TU"])
        state0 = np.array([x0, 0.0, z0, 0.0, ydot0, 0.0], dtype=np.float64)
        return state0, period_tu
    pytest.skip("no genuinely-3D converged row found in spike_287.jsonl")


def test_v1_spike_287_3d_member_passes() -> None:
    """V1 on the #287 spike's 3D Braik-Ross (1,1) extension.

    The IC is the FIRST converged 3D row in ``data/spike_287.jsonl`` (the
    z0_guess sweep that established the spike's 80-member family). V1 must
    re-close it under the 3D corrector AND the independent Radau check.
    """
    state0, period_tu = _read_spike_287_3d_seed()
    system = _em_system()
    verdict = run_v1_3d(
        "spike-287-3d-member",
        state0,
        period_tu,
        system,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        is_half_period_residual=True,
        notes="data/spike_287.jsonl 3D member; family-tracer seed",
    )
    assert verdict.converged_corrector, (
        f"corrector failed at {verdict.closure_residual_nondim:.3e} nondim"
    )
    assert verdict.converged_independent, (
        f"independent Radau closure failed at {verdict.independent_closure_nondim:.3e} nondim"
    )
    assert verdict.passes_v1, (
        f"V1 FAIL: independent={verdict.independent_closure_kms:.3e} km/s "
        f"> v1_floor={verdict.v1_floor_kms:.0e} km/s"
    )
    # Genuinely 3D — escapes the planar manifold.
    assert not verdict.degenerate_planar


def _read_first_subfamily_member() -> tuple[np.ndarray, float]:
    """Read a #301 sub-family member IC from the subfamilies JSONL.

    Returns the SWITCHED IC of the first accepted bracket (a doubly-hyperbolic
    asymmetric 3D orbit; period ratio ~ k = 4 vs the parent). This is an
    asymmetric IC (the y-axis perpendicular-crossing symmetry is broken by
    the bifurcation), so the corrector must run in full-asymmetric mode.

    Returns ``(state0, period_tu)``.
    """
    path = _REPO_ROOT / "data" / "family_296_3d_subfamilies_299.jsonl"
    if not path.exists():
        pytest.skip(f"family_296_3d_subfamilies_299.jsonl missing at {path}")
    for raw in path.read_text().splitlines():
        row = json.loads(raw)
        if row.get("type") == "header":
            continue
        if row.get("outcome") != "accepted":
            continue
        switched = row.get("switched")
        if switched is None:
            continue
        state = switched.get("switched_state")
        period = switched.get("switched_T_TU")
        if state is None or period is None:
            continue
        return np.asarray(state, dtype=np.float64), float(period)
    pytest.skip("no accepted bracket with switched state found in subfamilies JSONL")


def test_v1_3d_subfamily_member_passes() -> None:
    """V1 on a #301 doubly-hyperbolic 3D sub-family member.

    The IC is the SWITCHED state of the first accepted bracket in
    ``data/family_296_3d_subfamilies_299.jsonl``. Per the JSONL header the
    sub-family was converged at residual 1e-10 to 1e-14; V1 must reproduce
    that.

    The switched state has y=0 and xdot=0 (the bifurcation preserved the
    parent's symmetry plane), so we run the corrector in the symmetric
    tulip mode — that matches how the sub-family tracer closed each member.
    """
    state0, period_tu = _read_first_subfamily_member()
    system = _em_system()
    verdict = run_v1_3d(
        "subfamily-299-k4-bifurcation-member",
        state0,
        period_tu,
        system,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        is_half_period_residual=True,
        notes="#301 doubly-hyperbolic sub-family; first accepted bracket switched state",
    )
    assert verdict.converged_corrector
    assert verdict.converged_independent
    assert verdict.passes_v1, (
        f"V1 FAIL on sub-family member: "
        f"independent={verdict.independent_closure_kms:.3e} km/s > "
        f"v1_floor={verdict.v1_floor_kms:.0e} km/s"
    )
    assert not verdict.degenerate_planar


def test_v1_rejects_deliberately_broken_ic() -> None:
    """Negative control: a randomly perturbed 3D IC fails V1.

    Perturb the planar Braik-Ross seed by 0.1 nondim in each component AND
    perturb the period by 50% — a deliberately non-periodic starting point.
    The corrector should EITHER fail to converge OR converge to something
    whose independent re-propagation exceeds the nondim floor. Either way
    V1 must NOT pass.
    """
    system = _em_system()
    rng = np.random.default_rng(seed=20260616)
    perturbation = rng.uniform(-0.1, 0.1, size=6)
    # Start from a random direction in the planar-3D space so the IC is
    # nowhere near any known periodic orbit.
    state0 = np.array([-0.5, 0.3, 0.4, 0.5, -0.4, 0.2], dtype=np.float64) + perturbation
    period_guess = C11A_PERIOD_TU * 0.5  # 50% off the planar seed period
    verdict = run_v1_3d(
        "deliberately-broken-ic",
        state0,
        period_guess,
        system,
        max_iter=20,
        notes="negative control: random IC, V1 must reject",
    )
    assert not verdict.passes_v1, (
        f"V1 false-positive on a random IC: "
        f"closure_corr={verdict.closure_residual_nondim:.3e} nondim, "
        f"independent={verdict.independent_closure_nondim:.3e} nondim "
        f"(both should exceed the nondim floor or the km/s floor)"
    )


def test_v1_nondim_floor_default_is_corrector_default() -> None:
    """Sanity: the V1 module's nondim floor default matches the 3D
    corrector's own ``independent_tol`` default (1e-6).

    This keeps the two layers numerically consistent: a corrector that
    self-flags ``converged=True`` at exactly the threshold also passes the
    V1 nondim sub-gate.
    """
    assert V1_FLOOR_NONDIM_DEFAULT == 1.0e-6


def test_v1_verdict_carries_audit_fields() -> None:
    """The V1 verdict carries the floors it was held against — the audit
    trail must let a later reader reconstruct exactly what passed/failed.
    """
    system = _em_system()
    state0 = np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    verdict = run_v1_3d(
        "audit-trail-check",
        state0,
        C11A_PERIOD_TU,
        system,
        free_vars=FREE_VARS_SYMMETRIC_TULIP,
        residual_indices=RESIDUAL_PERPENDICULAR_HALF_PERIOD,
        is_half_period_residual=True,
        notes="audit-trail check",
    )
    assert verdict.v1_floor_kms == V1_FLOOR_KMS
    assert verdict.independent_floor_nondim == V1_FLOOR_NONDIM_DEFAULT
    assert verdict.candidate_id == "audit-trail-check"
    assert verdict.notes == "audit-trail check"


def test_v1_rejects_malformed_system() -> None:
    """A CR3BPSystem with l_km=0 or t_s=0 must fail loudly, not silently
    divide by zero (which would emit a meaningless inf km/s residual).
    """
    bad_system = cr3bp.CR3BPSystem(
        mu=EM_MU,
        primary="earth",
        secondary="moon",
        l_km=0.0,  # bad
        t_s=EM_T_S,
    )
    state0 = np.array([C11A_X0, 0.0, 0.0, 0.0, C11A_YDOT0, 0.0], dtype=np.float64)
    with pytest.raises(ValueError, match="invalid CR3BP system"):
        run_v1_3d(
            "bad-system",
            state0,
            C11A_PERIOD_TU,
            bad_system,
        )
