"""C_J=3.1294 method-validation gate for the Braik-Ross reachable-set scorer.

REPRODUCE-BEFORE-TRUST: before the scorer is allowed to rank our families it must
reproduce Braik & Ross 2026's (arXiv:2605.31543) published structural result at
their common energy C_J = 3.1294 (paper Table 4 / Fig. 10):

  * the (3,2)-cycler C32 is the DOMINANT family (rank 1 in strength, harmonic
    closeness, AND betweenness -- Table 4 values 0.2850 / 0.2891 / 0.5000), and
  * the 2:1 stable resonant R21-S is the persistent HARD-ACCESS family (last in
    strength and closeness, zero betweenness).

RECOVERY GAP (honest scoping; see the results note
``docs/notes/2026-06-13-reachable-scorer-ungate.md``): Braik-Ross publish NO state
vectors, only periods + sigma. We recover each member at C_J=3.1294 and confirm it
against BOTH the sourced Table-2 period AND the sourced Floquet rate sigma before
scoring (the sigma check rejects spurious same-period orbits). With the
network-independent free-(x0, t_half) corrector, NINE of the thirteen members
recover and confirm offline (LL1, LL2, DPO, R21-S, R21-U, R31-S, R31-U, R52-S,
C21). The four that do NOT recover at this off-stable common energy with the
available single-/free-period shooting correctors are the three unstable cyclers
(C11a sigma=1.05, C11b sigma=0.93, C32 sigma=0.69 -- each collapses onto a nearby
spurious lower-sigma orbit) and the 5:2 unstable resonant R52-U; they are excluded
rather than faked. The gate therefore runs on the nine source-confirmable nodes.
Crucially C32 -- the family whose dominance is the headline claim -- is itself
NOT faithfully recoverable here, so the C32-dominance gate cannot be tested. This
is a faithfully-scoped reproduction, not a forced one; parameters are NOT tuned.

All EXPECTED values trace to a published source: sourced periods + sigma to
Braik-Ross Table 2, and the C32-dominant / R21-S-hard-access ranking to Table 4 /
Fig. 10. Recovered ICs are derived quantities (never goldens).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

import cyclerfinder.search.reachable_network as rn
import cyclerfinder.search.reachable_representatives as rr


def _recover_subset() -> list[rr.Representative]:
    """Recover the network-independent source-confirmable gate subset at C_J=3.1294.

    Period AND Floquet sigma are both confirmed against Braik-Ross Table 2 (no
    JPL network call). Returns only the members that confirm.
    """
    sysm = rr.braik_ross_system()
    return [r for r in rr.recover_offline_set(sysm) if r.confirmed]


@pytest.mark.slow
def test_recover_representatives_periods_and_sigma_confirmed() -> None:
    """Each gate-subset member recovers to its Table-2 period AND sigma (offline)."""
    reps = _recover_subset()
    by = {r.label: r for r in reps}
    for r in reps:
        src_d, src_s = rr.SOURCED_TABLE2[r.label]
        print(
            f"{r.label:6s} conv={r.converged!s:5s} T={r.period_days:8.3f} d "
            f"(sourced {src_d:7.3f}, sigma {src_s}) C={r.jacobi:.6f} x0={r.state0[0]:+.5f} "
            f"confirmed={r.confirmed}"
        )
    expected = ("LL1", "LL2", "DPO", "R21-S", "R21-U", "R31-S", "R31-U", "R52-S", "C21")
    for label in expected:
        assert label in by, f"{label} did not confirm (period+sigma)"
        r = by[label]
        assert abs(r.period_days - r.sourced_period_days) < 0.5, (
            f"{label}: recovered {r.period_days:.3f} d vs sourced {r.sourced_period_days} d"
        )


def _run_network(reps: list[rr.Representative]) -> tuple[list[str], rn.Centralities]:
    """Build the budget-capped Braik-Ross network and its normalized centralities."""
    grid = rn.VoxelGrid(dx=0.02, dy=0.02, dtheta=math.radians(10.0))
    sysm = rr.braik_ross_system()
    # Common accessibility horizon T_a = paper reference T_cap = 1 sidereal month
    # = 2*pi TU ~ 27.32 d (so accessibility differences reflect transport).
    horizon = 2.0 * math.pi
    forward = [
        rn.build_reachable_set(
            sysm,
            r.state0,
            r.period,
            grid,
            rr.C_J_BRAIK_ROSS,
            n_seeds=10,
            n_fan=9,
            delta_max=math.radians(30.0),
            horizon=horizon,
        )
        for r in reps
    ]
    backward = [rn.mirror_reachable_set(f, grid) for f in forward]
    mat_nd = rn.proxy_matrix(forward, backward, grid)
    # Braik-Ross edge-retention: drop edges over the max-budget cap (Eq. 54), in
    # m/s -- this is what creates relay routing / nonzero betweenness.
    capped = rn.apply_budget_cap(mat_nd, dv_cap_ms=rn.DV_CAP_MS)
    # Normalize against the full Nf=13 representative set (Table-4 convention).
    return [r.label for r in reps], rn.normalized_centralities(capped, n_families=13)


@pytest.mark.slow
def test_stable_resonants_are_hard_access_on_subset() -> None:
    """REPRODUCED part of the gate: the stable resonants are hard-access (bottom).

    On the nine source-confirmable nodes, the 2:1 stable resonant R21-S ranks in
    the bottom three of both strength and harmonic closeness, alongside the other
    stable resonants -- consistent with the Braik-Ross "stable resonant tori
    resist invasion" hard-access reading (Table 4: R21-S/R31-S/R52-S occupy the
    bottom ranks 13/12/10). This is the qualitatively reproduced half of the
    published structural result.
    """
    reps = _recover_subset()
    labels, cent = _run_network(reps)
    i_r21s = labels.index("R21-S")
    bottom3_strength = set(np.argsort(cent.strength)[:3].tolist())
    bottom3_closeness = set(np.argsort(cent.harmonic_closeness)[:3].tolist())
    assert i_r21s in bottom3_strength, _rank_msg("strength", labels, cent.strength)
    assert i_r21s in bottom3_closeness, _rank_msg(
        "harmonic_closeness", labels, cent.harmonic_closeness
    )


@pytest.mark.slow
@pytest.mark.xfail(
    reason=(
        "NOT TESTABLE on the source-confirmable subset: C32 is NOT among the nine "
        "recoverable members. C32 (and C11a, C11b, R52-U) are unstable orbits that "
        "collapse onto spurious lower-sigma members under the available shooting "
        "correctors at the off-stable common energy; a robust Jacobi-constrained "
        "multiple-shooting corrector (or the JPL oracle, unavailable here) is needed. "
        "C32's hub/gateway/relay dominance is intrinsically a full-13-node property, "
        "and C32 itself is not faithfully recovered -- so the headline gate cannot be "
        "exercised. A faithful negative: the scorer stays GATED for our families. "
        "See docs/notes/2026-06-13-reachable-scorer-ungate.md."
    ),
    strict=False,
)
def test_validation_gate_c32_dominant() -> None:
    """METHOD-VALIDATION GATE (published claim): C32 is the dominant family.

    Encodes the Braik-Ross Table-4 / Fig-10 claim verbatim so the negative is
    recorded honestly rather than silently dropped. Expected to xfail: C32 is not
    in the recoverable subset, so ``labels.index('C32')`` raises -- the gate cannot
    even be evaluated, which IS the honest negative. Parameters are NOT tuned.
    """
    reps = _recover_subset()
    labels, cent = _run_network(reps)
    print()
    for name, vals in (
        ("strength", cent.strength),
        ("harmonic", cent.harmonic_closeness),
        ("betweenness", cent.betweenness),
    ):
        print(_rank_msg(name, labels, vals))
    i_c32 = labels.index("C32")  # raises ValueError: C32 not recovered -> xfail
    assert int(np.argmax(cent.strength)) == i_c32
    assert int(np.argmax(cent.harmonic_closeness)) == i_c32
    assert int(np.argmax(cent.betweenness)) == i_c32


def _rank_msg(name: str, labels: list[str], values: np.ndarray) -> str:
    order = np.argsort(-values)
    ranked = ", ".join(f"{labels[i]}={values[i]:.4g}" for i in order)
    return f"{name} (high->low): {ranked}"
