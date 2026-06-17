"""#347 Phase 1 — P1.1 reproduction-gate for the (3,2) Earth-Moon symmetric anchor.

The anchor is Braik-Ross 2026 Table 2 row C32: P = 78.613 days, sigma_d = 0.1583
day^-1 at CJ = 3.1294. The IC seed (x0 = -0.2752115, ydot0_sign = -1,
half_crossings = 6) lives in
:data:`cyclerfinder.search.reachable_representatives._CYCLER_SEEDS["C32"]` and was
set up at #262 from the Ross & Roberts-Tsoukkas 2025 AAS-25-621 family seed
region. The full pipeline
:func:`cyclerfinder.search.reachable_representatives.recover_all_cyclers_braik_ross`
recovers a Representative whose period matches Braik-Ross to ~0.0005%.

This test adds three orthogonal cross-checks against the sourced numerics:

  1. Period gate vs Braik-Ross 78.613 d, within 1%.
  2. Topology gate via the independent winding-number classifier
     :func:`cyclerfinder.search.binary_star_search.winding_topology`. Must yield
     (k1, k2) = (3, 2) with prograde windings.
  3. Floquet sigma gate via the project's own monodromy +
     :func:`cyclerfinder.search.bifurcation_detector.floquet_multipliers`,
     converted to sigma_d via Braik-Ross eq. (20):
     sigma = ln(|lambda_max nontrivial|) / T (TU^-1), sigma_d = sigma / TU_DAYS.
     Must match the sourced 0.1583 day^-1 within 5%.

The gates are SOURCED on the expected side (Braik-Ross Table 2). The computed
side is the project's own pipeline; this is a reproduction test, not a
self-consistency test.

Discipline cross-reference: ``feedback_golden_tests_sourced_only``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.bifurcation_detector import (
    floquet_multipliers,
    monodromy,
)
from cyclerfinder.search.binary_star_search import winding_topology
from cyclerfinder.search.reachable_representatives import (
    TU_DAYS,
    Representative,
    braik_ross_system,
    recover_all_cyclers_braik_ross,
)

# Braik-Ross 2026 Table 2 sourced values for C32 (3,2)-cycler at CJ = 3.1294.
BRAIK_ROSS_C32_PERIOD_DAYS = 78.613
BRAIK_ROSS_C32_SIGMA_D = 0.1583  # day^-1
BRAIK_ROSS_C32_SIGMA_TU = 0.6886  # TU^-1; redundant cross-check
BRAIK_ROSS_CJ = 3.1294


def _recover_c32_representative() -> tuple[cr3bp.CR3BPSystem, Representative]:
    system = braik_ross_system()
    reps = recover_all_cyclers_braik_ross(system)
    for r in reps:
        if r.label == "C32":
            return system, r
    raise RuntimeError("recover_all_cyclers_braik_ross returned no C32 row")


def test_c32_anchor_period_matches_braik_ross_table2() -> None:
    """Period gate: recovered C32 period in days matches Braik-Ross 78.613 within 1%."""
    _system, rep = _recover_c32_representative()
    assert rep.converged, "C32 corrector did not converge"
    assert rep.confirmed, (
        f"C32 recovery flag .confirmed=False (period {rep.period_days:.4f} d "
        f"vs sourced {BRAIK_ROSS_C32_PERIOD_DAYS} d)"
    )
    rel_err = abs(rep.period_days - BRAIK_ROSS_C32_PERIOD_DAYS) / BRAIK_ROSS_C32_PERIOD_DAYS
    assert rel_err < 0.01, (
        f"C32 period rel-error {rel_err:.4%} exceeds 1% gate "
        f"(recovered {rep.period_days:.4f} d, sourced {BRAIK_ROSS_C32_PERIOD_DAYS} d)"
    )


def test_c32_anchor_jacobi_pinned_to_cj_braik_ross() -> None:
    """Jacobi gate: corrector holds CJ = 3.1294 to better than 1e-9 (Newton-floor)."""
    _system, rep = _recover_c32_representative()
    assert abs(rep.jacobi - BRAIK_ROSS_CJ) < 1e-9, (
        f"C32 Jacobi {rep.jacobi:.12f} drifted from sourced {BRAIK_ROSS_CJ}"
    )


def test_c32_anchor_winding_topology_is_3_2_prograde() -> None:
    """Topology gate: independent winding classifier yields (k1, k2) = (3, 2) prograde."""
    system, rep = _recover_c32_representative()
    topo = winding_topology(system.mu, rep.state0, rep.period)
    assert topo.k1 == 3, f"C32 winding k1={topo.k1}, expected 3"
    assert topo.k2 == 2, f"C32 winding k2={topo.k2}, expected 2"
    assert topo.prograde, f"C32 windings not both positive: w1={topo.w1}, w2={topo.w2}"


def test_c32_anchor_floquet_sigma_matches_braik_ross_eq20() -> None:
    """Floquet sigma gate: sigma_d = ln(|lambda_max nontriv|)/T (Braik-Ross eq. 20).

    Computes the monodromy at the anchor IC and converts the largest non-trivial
    Floquet multiplier to sigma_d in day^-1. Compares against Braik-Ross's
    sourced 0.1583 day^-1.

    The "non-trivial" multipliers are all except the two closest to +1 (the
    energy + time-translation pair, which sits at +1 exactly for an autonomous
    Hamiltonian periodic orbit and is split only by integrator round-off).
    """
    system, rep = _recover_c32_representative()
    mono = monodromy(system, rep.state0, rep.period)
    eigs = floquet_multipliers(mono)
    # Exclude the two eigenvalues closest to +1 (the trivial pair).
    order = np.argsort([abs(complex(e) - 1.0) for e in eigs])
    nontriv_indices = order[2:]
    nontriv_mags = [abs(complex(eigs[i])) for i in nontriv_indices]
    lam_max = max(nontriv_mags)
    assert lam_max > 1.0, (
        f"C32 nontrivial |lambda_max|={lam_max:.6f} <= 1 — the (3,2) cycler is unstable, "
        "expected |lambda_max| >> 1"
    )
    sigma_tu = math.log(lam_max) / rep.period
    sigma_d = sigma_tu / TU_DAYS

    rel_err_d = abs(sigma_d - BRAIK_ROSS_C32_SIGMA_D) / BRAIK_ROSS_C32_SIGMA_D
    assert rel_err_d < 0.05, (
        f"C32 Floquet sigma_d {sigma_d:.4f} day^-1 deviates from "
        f"Braik-Ross Table 2 {BRAIK_ROSS_C32_SIGMA_D} day^-1 by {rel_err_d:.4%} "
        f"(>5% gate)"
    )
    rel_err_tu = abs(sigma_tu - BRAIK_ROSS_C32_SIGMA_TU) / BRAIK_ROSS_C32_SIGMA_TU
    assert rel_err_tu < 0.05, (
        f"C32 Floquet sigma_TU {sigma_tu:.4f} TU^-1 deviates from "
        f"Braik-Ross Table 2 {BRAIK_ROSS_C32_SIGMA_TU} TU^-1 by {rel_err_tu:.4%} "
        f"(>5% gate)"
    )


@pytest.mark.parametrize(
    "label,expected_period_days",
    [
        ("C11a", 42.140),
        ("C11b", 55.995),
        ("C21", 84.533),
        ("C32", 78.613),
    ],
)
def test_all_four_braik_ross_cyclers_recover_within_1pct(
    label: str, expected_period_days: float
) -> None:
    """Sanity-cross-check: all four Braik-Ross cyclers reproduce within 1% on period.

    This is a regression breadcrumb; the per-family Jacobi (C21 unrounded vs
    others at literal 3.1294) is encoded in
    :data:`cyclerfinder.search.reachable_representatives._CYCLER_SEEDS`. If a
    refactor breaks ANY of the four, P1.x downstream is suspect.
    """
    system = braik_ross_system()
    reps = recover_all_cyclers_braik_ross(system)
    by_label = {r.label: r for r in reps}
    rep = by_label[label]
    assert rep.converged, f"{label}: corrector did not converge"
    rel_err = abs(rep.period_days - expected_period_days) / expected_period_days
    assert rel_err < 0.01, (
        f"{label}: period {rep.period_days:.4f} d vs sourced "
        f"{expected_period_days} d (rel err {rel_err:.4%}, gate 1%)"
    )
