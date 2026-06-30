"""C_J=3.1294 method-validation gate for the Braik-Ross reachable-set scorer.

REPRODUCE-BEFORE-TRUST: before the scorer is allowed to rank our families it must
reproduce Braik & Ross 2026's (arXiv:2605.31543) published structural result at
their common energy C_J = 3.1294 (paper Table 4 / Fig. 10):

  * the (3,2)-cycler C32 is the DOMINANT family (rank 1 in strength, harmonic
    closeness, AND betweenness -- Table 4 values 0.2850 / 0.2891 / 0.5000), and
  * the 2:1 stable resonant R21-S is the persistent HARD-ACCESS family (last in
    strength and closeness, zero betweenness).

RECOVERY STATUS (#262, post-#249 4/4 cycler recovery): all four Braik-Ross
cycler members (C11a, C11b, C21, C32) are now rigorously reproduced via the
symmetric perpendicular-x-axis-crossing corrector at the correct per-family
Jacobi (literal :data:`C_J_BRAIK_ROSS` = 3.1294 for C11a/C11b/C32; the unrounded
:data:`C_J_C21` = 3.129389531088256 for C21, whose (2,1) family spans ΔC ~ 4e-12
and does not exist at the literal printed value). Combined with the eight
network-independent offline confirmations (LL1, LL2, DPO, R21-S, R21-U, R31-S,
R31-U, R52-S), the gate now runs on TWELVE source-confirmable nodes -- the
full Braik-Ross representative set minus only the 5:2 unstable resonant R52-U
(still unrecovered; excluded rather than faked).

All EXPECTED values trace to a published source: sourced periods + sigma to
Braik-Ross Table 2, and the C32-dominant / R21-S-hard-access ranking to Table 4 /
Fig. 10. Recovered ICs are derived quantities (never goldens).
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
import pytest

import cyclerfinder.search.reachable_network as rn
import cyclerfinder.search.reachable_representatives as rr

_PUBLISHED_DVMATRIX = Path("data/golden/braik_ross_2026_dvmatrix_mps.csv")


def _load_published_dvmatrix() -> tuple[list[str], np.ndarray]:
    """Braik-Ross 2026 published 13x13 proxy-dV matrix (m/s), NaN -> inf (no edge)."""
    rows = list(csv.reader(_PUBLISHED_DVMATRIX.read_text().splitlines()))
    header = rows[0][1:]
    mat = np.array(
        [[float(x) if x.strip().lower() != "nan" else math.inf for x in r[1:]] for r in rows[1:]],
        dtype=np.float64,
    )
    np.fill_diagonal(mat, 0.0)
    return header, mat


def test_centrality_scorer_reproduces_braik_ross_table4() -> None:
    """Centrality scorer reproduces Braik-Ross Table 4 from the PUBLISHED dV matrix (#497).

    Fed Braik & Ross 2026's own published proxy-dV matrix (data/golden/
    braik_ross_2026_dvmatrix_mps.csv, adopted #495), ``normalized_centralities``
    reproduces the Table-4 C32-dominant values to print precision (0.2850 / 0.2891 /
    0.5000) and C32 is the argmax of all three metrics. EXPECTED values are the
    PUBLISHED Table 4 (sourced, never circular). This ISOLATES the #249 gate failure
    (``test_validation_gate_c32_dominant``, xfail) to OUR proxy-dV fidelity -- the
    centrality math itself is correct.
    """
    header, mat = _load_published_dvmatrix()
    cent = rn.normalized_centralities(mat, n_families=13)
    ic = header.index("Cycler 32")
    assert cent.strength[ic] == pytest.approx(0.2850, abs=5e-4)
    assert cent.harmonic_closeness[ic] == pytest.approx(0.2891, abs=5e-4)
    assert cent.betweenness[ic] == pytest.approx(0.5000, abs=5e-4)
    assert int(np.argmax(cent.strength)) == ic
    assert int(np.argmax(cent.harmonic_closeness)) == ic
    assert int(np.argmax(cent.betweenness)) == ic


def _recover_subset() -> list[rr.Representative]:
    """Recover the source-confirmable gate set at C_J=3.1294 (#262 / post-#249).

    Combines the offline source-confirmable nodes (LL1, LL2, DPO, R21-S, R21-U,
    R31-S, R31-U, R52-S) with the four Braik-Ross cyclers (C11a, C11b, C21, C32)
    recovered via :func:`rr.recover_all_cyclers_braik_ross`. C21 is recovered at
    its own (unrounded) Jacobi :data:`rr.C_J_C21` (where the (2,1) family
    actually lives); the other three cyclers + the eight offline nodes are all
    at :data:`rr.C_J_BRAIK_ROSS`. Period AND (for the offline nodes) Floquet
    sigma are both confirmed against Braik-Ross Table 2; we dedupe by label so
    the C21 in the offline set is replaced by the cycler-route version.

    Returns only the members that confirm (no faked members enter the network).
    """
    sysm = rr.braik_ross_system()
    by_label: dict[str, rr.Representative] = {}
    for r in rr.recover_offline_set(sysm):
        if r.confirmed:
            by_label[r.label] = r
    for r in rr.recover_all_cyclers_braik_ross(sysm):
        if r.confirmed:
            by_label[r.label] = r  # cycler-route C21 supersedes offline C21
    # Canonical ordering for stable label/value reporting.
    order = (
        "LL1",
        "LL2",
        "DPO",
        "R21-S",
        "R21-U",
        "R31-S",
        "R31-U",
        "R52-S",
        "C11a",
        "C11b",
        "C21",
        "C32",
    )
    return [by_label[label] for label in order if label in by_label]


@pytest.mark.slow
def test_recover_representatives_periods_and_sigma_confirmed() -> None:
    """Each gate-subset member recovers to its Table-2 period (offline + cyclers).

    Offline nodes also confirm against the sourced Floquet sigma; the four
    cyclers are confirmed by period + (k1,k2) winding topology + prograde +
    Radau in :mod:`tests.search.test_cr3bp_ross_families` (#249 4/4 reproduction)
    -- here we only check the period close-to-sourced gate.
    """
    reps = _recover_subset()
    by = {r.label: r for r in reps}
    for r in reps:
        src_d, src_s = rr.SOURCED_TABLE2[r.label]
        print(
            f"{r.label:6s} conv={r.converged!s:5s} T={r.period_days:8.3f} d "
            f"(sourced {src_d:7.3f}, sigma {src_s}) C={r.jacobi:.6f} x0={r.state0[0]:+.5f} "
            f"confirmed={r.confirmed}"
        )
    expected = (
        "LL1",
        "LL2",
        "DPO",
        "R21-S",
        "R21-U",
        "R31-S",
        "R31-U",
        "R52-S",
        "C11a",
        "C11b",
        "C21",
        "C32",
    )
    for label in expected:
        assert label in by, f"{label} did not confirm (period+sigma or period)"
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
        "FAITHFUL NEGATIVE on the 12-node source-confirmable set (#262, post-#249 "
        "4/4 cycler recovery): with all four Braik-Ross cyclers (C11a, C11b, C21, "
        "C32) and the eight offline-confirmable resonant/libration/DPO nodes in "
        "play, C32 does NOT emerge as the dominant family node under our scorer. "
        "Observed ranking on this run: strength argmax = C11a; harmonic-closeness "
        "argmax = C21; betweenness argmax = R21-U; C32 has ZERO betweenness (no "
        "relay role at all). DIAGNOSED (#497, post-#495 golden adoption): the cause "
        "is OUR proxy-dV FIDELITY, not the centrality math nor a missing node. The "
        "centrality scorer reproduces Braik Table 4 EXACTLY from the published dV "
        "matrix (test_centrality_scorer_reproduces_braik_ross_table4 PASSES: "
        "0.2850/0.2891/0.5000). Our heading-fan proxy OVERESTIMATES dV (per #495: "
        "dc_refined < proxy in every pair), so recalibrating DV_CAP_MS to the Braik "
        "51 m/s reference EMPTIES our network (all centralities 0), and at 409.3 our "
        "proxy mis-ranks betweenness to R21-U. The fix is proxy-dV calibration / "
        "rebuilding from the adopted dc_refined golden -- NOT a cap value, so the "
        "test stays xfail and the scorer stays GATED for OUR proxy. Parameters are "
        "NOT tuned. See docs/notes/2026-06-30-497-c32-gate-diagnosis.md."
    ),
    strict=True,
)
def test_validation_gate_c32_dominant() -> None:
    """METHOD-VALIDATION GATE (published claim): C32 is the dominant family.

    Encodes the Braik-Ross Table-4 / Fig-10 claim verbatim so the negative is
    recorded honestly rather than silently dropped. With the 4/4 cycler set in
    play (#262) the gate is now evaluable; expected to xfail because our scorer's
    ranking does not promote C32 (see xfail reason above). Parameters are NOT
    tuned.
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
    i_c32 = labels.index("C32")
    assert int(np.argmax(cent.strength)) == i_c32
    assert int(np.argmax(cent.harmonic_closeness)) == i_c32
    assert int(np.argmax(cent.betweenness)) == i_c32


@pytest.mark.slow
def test_validation_gate_c32_undominant_faithful_negative() -> None:
    """RECORD the faithful negative: C32 is NOT the dominant node on our 12-set.

    Companion to ``test_validation_gate_c32_dominant`` (xfail): asserts the
    *observed* ranking on the post-#249 4/4-cycler 12-node source-confirmable
    set so the negative is captured as a passing test (not just an xfail), and
    will fail-loud if the ranking ever changes (e.g. if we add R52-U later or if
    a corrector change shifts the scoring). Parameters are NOT tuned -- this
    test reports what the unmodified scorer produces.

    Findings (commit-time): C32 is below the median in strength AND in harmonic
    closeness, and has zero betweenness. The hub/gateway/relay node identities
    differ from Braik-Ross Table 4 -- which means the scorer is not yet a
    reliable family-selection prioritizer on our families.
    """
    reps = _recover_subset()
    labels, cent = _run_network(reps)
    i_c32 = labels.index("C32")
    # C32 is NOT rank 1 in any metric on this set.
    assert int(np.argmax(cent.strength)) != i_c32, _rank_msg("strength", labels, cent.strength)
    assert int(np.argmax(cent.harmonic_closeness)) != i_c32, _rank_msg(
        "harmonic_closeness", labels, cent.harmonic_closeness
    )
    # And C32 has no relay role (betweenness == 0).
    assert cent.betweenness[i_c32] == 0.0, _rank_msg("betweenness", labels, cent.betweenness)


def _rank_msg(name: str, labels: list[str], values: np.ndarray) -> str:
    order = np.argsort(-values)
    ranked = ", ".join(f"{labels[i]}={values[i]:.4g}" for i in order)
    return f"{name} (high->low): {ranked}"
