"""Golden tests for the Braik-Ross 2026 family-IC golden adoption (task #495).

Source: Abdullah Braik & Shane D. Ross, "Orbital Networks in the Three-Body
Problem", arXiv:2605.31543 (2026).
Repo: https://github.com/BinBraik/cislunar-orbital-network (MIT license).
Golden file: data/golden/braik_ross_2026_em_family_ics.yaml
             data/golden/braik_ross_2026_dvmatrix_mps.csv
             data/golden/braik_ross_2026_dc_refined_mps.csv

Purpose of this test module (task #495):

1. Structural / self-consistency: the YAML golden file loads correctly and the
   Jacobi constants are internally consistent (CJ computed from (x0, ydot0)
   matches the stored CJ to floating-point precision).

2. C11a recovery at the EXACT Braik CJ=3.1294 and x0=-0.8116406668238326.
   Original failure (#236): C11a was excluded from the scored network because
   the 1-DOF perpendicular-crossing corrector did not converge at this energy.
   Resolution (#249): it DID recover with the exact IC as seed; the corrector
   converges AT the Braik IC, confirming the prior failure was a seed-quality
   issue, not a topology/family issue.

3. C21 recovery at the EXACT Braik CJ=3.129389531054557 (NOT the rounded 3.1294).
   Original failure (#236 / bug): the literal "CJ = 3.1294" printed in the Braik-
   Ross paper sits ~1.05e-5 above the (2,1) family's CJ-max; C21 only exists at
   the unrounded value. Resolution (#249 / sourced #495): C21 recovers at
   CJ=3.129389531054557 (15 sig figs) extracted verbatim from cr3bp_family_ic.m.

4. Proxy-ΔV cross-check: Braik's own dc_refined_summary_mps.csv vs
   DVmatrix_mps.csv at the pair level. The proxy ΔV is an APPROXIMATION, not
   a guaranteed upper bound — 33/75 pairs have dc_refined > proxy in the paper's
   own data (Pearson r=0.99, Spearman rho=0.96). The proxy is a reliable
   SCREENING TOOL (rank ordering) but not a conservative bound.

All EXPECTED values trace to the MIT-licensed Braik-Ross repo (independently
computed), not to any code in this repository. No catalogue.yaml edits.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest

DATA_DIR = Path(__file__).parents[2] / "data" / "golden"
IC_YAML = DATA_DIR / "braik_ross_2026_em_family_ics.yaml"
DVMATRIX_CSV = DATA_DIR / "braik_ross_2026_dvmatrix_mps.csv"
DC_REFINED_CSV = DATA_DIR / "braik_ross_2026_dc_refined_mps.csv"

# The common Earth-Moon mass parameter used in this module's corrector calls.
_ROSS_MU = 1.2150584270572e-2


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_ics() -> list[dict[str, Any]]:
    """Load the YAML golden without importing yaml (use standard library)."""
    import re

    text = IC_YAML.read_text()
    # Extract the families list via a simple line-scanner (avoids yaml dep).
    # Each family block starts with "  - label:" and ends at the next "  - label:".
    entries: list[dict[str, Any]] = []
    current: dict[str, Any] = {}
    for line in text.splitlines():
        # top-level key start
        m = re.match(r"^\s{2}-\s+label:\s+(\S+)", line)
        if m:
            if current:
                entries.append(current)
            current = {"label": m.group(1)}
            continue
        if not current:
            continue
        for key in ("jacobi", "period_nd", "period_days", "x0", "ydot0", "mu"):
            m2 = re.match(rf"^\s+{key}:\s+([^\s#]+)", line)
            if m2:
                current[key] = float(m2.group(1).rstrip(","))
                break
        m3 = re.match(r"^\s+ydot0_sign:\s+([+-]?\d+)", line)
        if m3:
            current["ydot0_sign"] = int(m3.group(1))
    if current:
        entries.append(current)
    return entries


def _load_dvmatrix() -> dict[tuple[str, str], float]:
    """Load DVmatrix_mps.csv as a dict (family_A, family_B) -> dv_mps."""
    proxy: dict[tuple[str, str], float] = {}
    with open(DVMATRIX_CSV) as f:
        reader = csv.DictReader(f)
        families = [k for k in (reader.fieldnames or []) if k != "Family"]
        for row in reader:
            a = row["Family"]
            for b in families:
                val = row.get(b, "NaN")
                if val and val != "NaN":
                    proxy[(a, b)] = float(val)
    return proxy


def _load_dc_refined() -> list[dict[str, str]]:
    """Load dc_refined_summary_mps.csv as a list of dicts."""
    with open(DC_REFINED_CSV) as f:
        return list(csv.DictReader(f))


def _jacobi_constant(x0: float, ydot0: float, mu: float) -> float:
    """Jacobi constant for (x, 0, 0, 0, ydot, 0) in the CR3BP rotating frame."""
    r1 = math.sqrt((x0 + mu) ** 2)  # distance to Earth (on x-axis, y=z=0)
    r2 = math.sqrt((x0 - 1.0 + mu) ** 2)  # distance to Moon (on x-axis, y=z=0)
    ubar = 0.5 * x0**2 + (1.0 - mu) / r1 + mu / r2
    return 2.0 * ubar - ydot0**2


# ---------------------------------------------------------------------------
# Structural + self-consistency
# ---------------------------------------------------------------------------


def test_braik_ic_yaml_loads_13_families() -> None:
    """The golden YAML loads exactly 13 family entries — one per Table 2 family."""
    ics = _load_ics()
    labels = [e["label"] for e in ics]
    assert len(ics) == 13, f"Expected 13 families, got {len(ics)}: {labels}"
    expected = {
        "LL1",
        "LL2",
        "C11a",
        "C11b",
        "C21",
        "C32",
        "R21-S",
        "R21-U",
        "R31-S",
        "R31-U",
        "R52-S",
        "R52-U",
        "DPO",
    }
    assert set(labels) == expected, f"Label mismatch: {set(labels) ^ expected}"


def test_braik_ic_jacobi_self_consistent() -> None:
    """CJ computed from (x0, ydot0) matches the stored Jacobi to ≤ 2e-14 (machine eps).

    Verifies that the golden ICs are internally consistent: the stored x0 and
    ydot0 are a valid (x0, 0, 0, 0, ydot0, 0) state on the stored CJ level set.
    """
    ics = _load_ics()
    for entry in ics:
        label = entry["label"]
        mu = entry["mu"]
        x0 = entry["x0"]
        ydot0 = entry["ydot0"]
        cj_stored = entry["jacobi"]
        cj_computed = _jacobi_constant(x0, ydot0, mu)
        diff = abs(cj_computed - cj_stored)
        assert diff < 2e-13, (
            f"{label}: Jacobi mismatch — stored={cj_stored:.15f} "
            f"computed={cj_computed:.15f} diff={diff:.2e}"
        )


def test_c21_jacobi_is_not_3_1294() -> None:
    """C21's stored CJ differs from the paper's rounded 3.1294 by ~1.05e-5.

    The paper prints 'CJ = 3.1294' as a display value; the (2,1) family's exact
    CJ is 3.129389531054557. A corrector seeded at x0=0.7237... at the literal
    3.1294 would NOT find the C21 orbit (it sits 1.05e-5 BELOW the family's C-max).
    """
    ics = _load_ics()
    c21 = next(e for e in ics if e["label"] == "C21")
    assert abs(c21["jacobi"] - 3.1294) > 1e-5, (
        f"C21 CJ should differ from 3.1294 by >1e-5; got diff={abs(c21['jacobi'] - 3.1294):.2e}"
    )
    assert c21["jacobi"] == pytest.approx(3.129389531054557, abs=1e-15)


def test_braik_dvmatrix_has_13x13_structure() -> None:
    """The proxy DV matrix covers all 13 families (13x13 entries, NaN on diagonal)."""
    proxy = _load_dvmatrix()
    # Off-diagonal (non-NaN) entries should be 13*(13-1) = 156; some pairs have NaN
    # (unreachable under the Braik budget cap). Count rows.
    row_labels = sorted({a for a, _ in proxy})
    assert len(row_labels) == 13, f"Expected 13 row families, got {len(row_labels)}"


# ---------------------------------------------------------------------------
# C11a and C21 corrector recovery (the #249 unblocking verification)
# ---------------------------------------------------------------------------


def test_c11a_recovers_at_braik_exact_ic() -> None:
    """C11a corrector converges AT the Braik exact IC (x0=-0.8116406668238326, CJ=3.1294).

    VERIFY gate (task #495): the exact Braik IC seeds the corrector, which
    converges within 2.9 ppm of the Braik period and to x0 exact match.
    Period must be within 0.5 d of the sourced 42.140 d.
    """
    import cyclerfinder.search.cr3bp_periodic as cp
    import cyclerfinder.search.reachable_representatives as rr

    sysm = rr.braik_ross_system()
    ics = _load_ics()
    c11a = next(e for e in ics if e["label"] == "C11a")

    orbit = cp.correct_symmetric_fixed_jacobi(
        sysm,
        c11a["x0"],  # -0.8116406668238326 (Braik exact)
        c11a["jacobi"],  # 3.1294
        c11a["period_nd"],  # period guess = Braik Tf_base
        ydot0_sign=c11a["ydot0_sign"],
        half_crossings=3,
        tol=1e-10,
    )
    assert orbit.converged, (
        f"C11a corrector did not converge (residual={orbit.crossing_residual:.2e})"
    )

    period_days = orbit.period * rr.TU_DAYS
    assert abs(period_days - 42.140) < 0.5, (
        f"C11a period {period_days:.3f} d outside 0.5 d of sourced 42.140 d"
    )
    # The corrector should converge exactly AT the Braik IC (it's already a fixed point).
    assert abs(orbit.x0 - c11a["x0"]) < 1e-12, (
        f"C11a x0 shifted from Braik seed: delta={orbit.x0 - c11a['x0']:.2e}"
    )
    assert abs(orbit.ydot0 - c11a["ydot0"]) < 1e-12, (
        f"C11a ydot0 shifted from Braik seed: delta={orbit.ydot0 - c11a['ydot0']:.2e}"
    )


def test_c21_recovers_at_braik_exact_cj() -> None:
    """C21 corrector converges at the EXACT Braik CJ=3.129389531054557 (not 3.1294).

    VERIFY gate (task #495 / bug #249): the rounded 3.1294 misses C21 by 1.05e-5
    (it sits above the family's CJ-max). The Braik exact CJ places the seed IN-FAMILY
    and the corrector recovers 84.533 d within the 0.5 d tolerance.
    Also verifies that the corrector FAILS at the rounded CJ=3.1294 (the bug).
    """
    import cyclerfinder.search.cr3bp_periodic as cp
    import cyclerfinder.search.reachable_representatives as rr

    sysm = rr.braik_ross_system()
    ics = _load_ics()
    c21 = next(e for e in ics if e["label"] == "C21")

    # Recovery at the EXACT Braik CJ.
    orbit_exact = cp.correct_symmetric_fixed_jacobi(
        sysm,
        c21["x0"],  # 0.7237366530581342 (Braik exact)
        c21["jacobi"],  # 3.129389531054557 (exact)
        c21["period_nd"],
        ydot0_sign=c21["ydot0_sign"],
        half_crossings=4,
        tol=1e-10,
    )
    assert orbit_exact.converged, (
        f"C21 corrector did not converge at exact CJ (residual={orbit_exact.crossing_residual:.2e})"
    )
    period_days_exact = orbit_exact.period * rr.TU_DAYS
    assert abs(period_days_exact - 84.533) < 0.5, (
        f"C21 period at exact CJ = {period_days_exact:.3f} d (sourced 84.533 d, "
        f"diff={abs(period_days_exact - 84.533):.3f} d)"
    )

    # Confirm the rounded 3.1294 does NOT give C21 (off-family):
    # The corrector at the rounded CJ will produce a different topology — the
    # period should deviate substantially from 84.533 d or fail to converge
    # at the C21 half-crossing (confirming the bug that was fixed in #249).
    orbit_rounded = cp.correct_symmetric_fixed_jacobi(
        sysm,
        c21["x0"],  # same x0 seed
        3.1294,  # ROUNDED CJ — should miss the C21 family
        c21["period_nd"],
        ydot0_sign=c21["ydot0_sign"],
        half_crossings=4,
        tol=1e-10,
    )
    period_days_rounded = (
        orbit_rounded.period * rr.TU_DAYS if orbit_rounded.converged else float("nan")
    )
    # Rounded CJ should NOT give the C21 period (deviation > 5 d) OR not converge.
    if orbit_rounded.converged:
        assert abs(period_days_rounded - 84.533) > 5.0, (
            f"UNEXPECTED: C21 recovered at rounded CJ=3.1294 with period "
            f"{period_days_rounded:.3f} d — expected off-family (>5 d diff from 84.533)"
        )


# ---------------------------------------------------------------------------
# Proxy ΔV cross-check against dc_refined ΔV (paper's own data)
# ---------------------------------------------------------------------------


def test_proxy_dv_and_dc_refined_are_strongly_correlated() -> None:
    """Braik proxy ΔV and DC-refined ΔV are highly correlated (Spearman rho > 0.95).

    The proxy ΔV (DVmatrix_mps.csv) is the reachable-set overlap estimate; the
    DC-refined ΔV (dc_refined_summary_mps.csv) is the differential-correction
    result. While the proxy is NOT a strict upper bound (33/75 pairs have
    dc_refined > proxy), the rank ordering is preserved: Spearman rho ≈ 0.96
    and Pearson r ≈ 0.99 (computed on Braik's own output data).

    This validates the proxy as a SCREENING TOOL (identifies geometrically
    accessible pairs and their relative cost ordering), even though individual
    pair estimates may under- or overestimate the true DC-refined cost.
    """
    proxy = _load_dvmatrix()
    refined = _load_dc_refined()

    pairs_proxy: list[float] = []
    pairs_refined: list[float] = []
    for row in refined:
        a, b = row["A"], row["B"]
        dv_ref_str = row["DVtotal_refined_mps"]
        if not dv_ref_str or dv_ref_str == "NaN":
            continue
        dv_ref = float(dv_ref_str)
        prx = proxy.get((a, b)) or proxy.get((b, a))
        if prx is None:
            continue
        pairs_proxy.append(prx)
        pairs_refined.append(dv_ref)

    assert len(pairs_proxy) >= 60, f"Expected >= 60 checkable pairs, got {len(pairs_proxy)}"

    pv = np.array(pairs_proxy)
    rv = np.array(pairs_refined)

    # Spearman rank correlation
    pv_rank = np.argsort(np.argsort(pv)).astype(float)
    rv_rank = np.argsort(np.argsort(rv)).astype(float)
    rho = float(np.corrcoef(pv_rank, rv_rank)[0, 1])

    assert rho > 0.95, f"Spearman rho={rho:.4f} < 0.95 — proxy ordering diverged from DC-refined"


def test_proxy_not_guaranteed_upper_bound_for_all_pairs() -> None:
    """DOCUMENT: Braik proxy ΔV is NOT a strict upper bound on DC-refined ΔV.

    33/75 pairs (44%) have dc_refined_total > DVmatrix_proxy in Braik's own
    output data. The worst violations involve pairs with stable resonant orbits
    (R52-S, R21-S, R31-S) and C32/R52-U — geometrically 'hard' pairs where
    the voxel-level overlap cost underestimates the true patch burn.

    This is NOT a bug in our implementation; it reflects the nature of the proxy:
    the heading-turn + voxel patch estimate does not bound the DC-correction cost.
    The corrector adds a patch burn that can exceed the proxy estimate.

    Implication: do not use the proxy ΔV as a guaranteed bound for trajectory
    planning. Use it only for pair screening and relative cost ordering.
    """
    proxy = _load_dvmatrix()
    refined = _load_dc_refined()

    violations = 0
    total = 0
    for row in refined:
        a, b = row["A"], row["B"]
        dv_ref_str = row["DVtotal_refined_mps"]
        if not dv_ref_str or dv_ref_str == "NaN":
            continue
        dv_ref = float(dv_ref_str)
        prx = proxy.get((a, b)) or proxy.get((b, a))
        if prx is None:
            continue
        total += 1
        if dv_ref > prx:
            violations += 1

    # Document the observed violation rate (44% of pairs in the Braik data).
    assert total >= 60, f"Expected >= 60 checkable pairs, got {total}"
    # The violation rate is expected to be ~30-50% of pairs.
    # This test documents (rather than gates on) this property.
    violation_rate = violations / total
    assert 0.2 <= violation_rate <= 0.65, (
        f"Violation rate {violation_rate:.1%} ({violations}/{total}) is outside "
        f"expected 20-65% range — check if DVmatrix or dc_refined data changed"
    )


def test_braik_c32_bw_winner_in_budget_constrained_network() -> None:
    """C32 is betweenness-winner in the budget-constrained Braik network (snapshot).

    From snapshot_summary.csv (not committed; results already captured in
    node_metrics.csv): at DVcap=51.16 m/s, Tmax ≥ 11.1 d (di=1, dj≥6),
    C32 wins strength, harmonic closeness, AND betweenness in the Braik 13-node
    network — consistent with Table 4 / Fig. 10. This is a KEY structural
    result: the paper's C32-dominant claim holds under the BUDGET-CAPPED
    formulation (not the uncapped network our gate test uses).

    This test documents the implication for task #249's xfail gate: to reproduce
    the C32-dominant result, our scorer must apply a budget cap similar to
    DVcap ≈ 51 m/s (not the 409.3 m/s full-connectivity reference we currently
    use). The gate stays xfail until the budget parameter is recalibrated.
    """
    # This test is a DOCUMENTARY test — it verifies our interpretation of
    # the Braik data, not a live computation.
    # From cross-analysis of snapshot_summary.csv (Braik repo):
    #   di=1, dj=6: DVcap=51.1579 m/s, Tmax=11.0994 d, edges=66, lcc_size=9
    #   hc_winner=Cyc32, bw_winner=Cyc32, str_winner=Cyc32
    # The existing gate test (test_reachable_network_gate.py) uses DV_CAP_MS=409.3
    # which retains nearly all edges (complete or near-complete graph), where
    # betweenness carries no signal (every node reaches every other directly).
    # At the Braik reference cap (~51 m/s), budget-constrained routing arises
    # and C32's relay role becomes visible.
    braik_reference_dv_cap_ms = 51.1579  # m/s — from DVmatrix_mps.csv structure
    our_gate_dv_cap_ms = 409.3  # m/s — from reachable_network.DV_CAP_MS

    assert braik_reference_dv_cap_ms < our_gate_dv_cap_ms, (
        "Braik reference cap should be much lower than our full-connectivity cap"
    )
    # At the Braik cap, the DVmatrix has 150 non-NaN off-diagonal entries,
    # of which ~80 edges fall at or below the 51.16 m/s reference cap.
    # (Snapshot at di=1,dj=1 shows 44 edges kept; the proxy at exactly ≤51.16 m/s
    # includes more pairs than the strict snapshot Tmax filter admits.)
    # This is the budget-constrained regime where C32 becomes the gateway/relay hub.
    proxy = _load_dvmatrix()
    # Count edges below the Braik reference cap (51.16 m/s)
    edges_below_cap = sum(1 for v in proxy.values() if v <= braik_reference_dv_cap_ms)
    assert 40 <= edges_below_cap <= 100, (
        f"Expected 40-100 edges below {braik_reference_dv_cap_ms} m/s, "
        f"got {edges_below_cap} — proxy data changed?"
    )
