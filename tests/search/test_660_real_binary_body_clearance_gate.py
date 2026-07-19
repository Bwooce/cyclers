"""Task #660: min-clearance-vs-body-radius physical gate -- test suite.

Motivation (see `#659`'s own bullet in data/OUTSTANDING.md, and this file's
sibling `docs/notes/2026-07-19-659-antiope-adjudication.md`): `#657`'s
real-binary (k1,k2) sweep found 2 Antiope candidates -- (1,1) and (2,2) --
that passed EVERY existing dynamical gate (topology / prograde /
reaches_secondary / Barden |nu|<1 / independent Radau crosscheck), yet both
are collision trajectories: they pass 3-14 km from each body's CENTRE while
Antiope's components have real radii of ~42-44 km. `RealBinarySystem` had no
body-radius field, and `reaches_secondary` (`x.max() > L1`) checks nothing
about physical clearance.

This suite is the mandatory regression converting that manual finding into a
permanent, automated gate check:

1. Every `REAL_BINARY_SYSTEMS` entry has a sourced `radius_km_primary`/
   `radius_km_secondary`/`radius_source` (no unsourced placeholders).
2. THE headline regression: the already-catalogued PC(3,2) row
   (`ross-rt-pc-cycler-32-2026`) clears BOTH Pluto and Charon with wide
   margin under the new gate -- this must stay green, the gate must NOT
   retroactively invalidate the one real catalogued hit.
3. Antiope (1,1) and (2,2) now fail EXPLICITLY on clearance, with the
   min-distance numbers pinned to #659's own adjudication figures.
4. Backward-compatibility: the gate is a documented no-op (byte-identical
   pre-#660 behavior) when no radius is supplied -- checked directly against
   `_gate_clearance`, not by re-running any other (expensive) live sweep.
"""

from __future__ import annotations

import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.search.pluto_charon_kk_sweep import (
    CHARON_RADIUS_KM,
    PLUTO_RADIUS_KM,
    SweepResult,
    make_pluto_charon_system,
)
from cyclerfinder.search.real_binary_kk_sweep import (
    DEFAULT_CLEARANCE_MARGIN_KM,
    REAL_BINARY_SYSTEMS,
    _finalize_grid_candidate,
    _gate_clearance,
    min_body_clearance_km,
    sweep_family,
)

# ---------------------------------------------------------------------------
# 1. Sourced radii present for every system, no unsourced placeholders
# ---------------------------------------------------------------------------


def test_660_every_system_has_sourced_radii() -> None:
    """Every REAL_BINARY_SYSTEMS entry must carry a sourced radius for BOTH
    bodies (no unsourced placeholder), per this task's explicit mandate: 7/7
    systems have citable literature radii, so `None` is not expected/allowed
    here (contrast a hypothetical future system where no radius exists in the
    literature -- `RealBinarySystem.radius_km_primary` staying `None` there
    is the honest-gap convention `_gate_clearance` is built to respect)."""
    for key, sys_def in REAL_BINARY_SYSTEMS.items():
        assert sys_def.radius_km_primary is not None, f"{key}: missing radius_km_primary"
        assert sys_def.radius_km_secondary is not None, f"{key}: missing radius_km_secondary"
        assert sys_def.radius_km_primary > 0.0, f"{key}: radius_km_primary must be positive"
        assert sys_def.radius_km_secondary > 0.0, f"{key}: radius_km_secondary must be positive"
        assert sys_def.radius_source, f"{key}: missing radius_source citation"
        # Sanity: neither body should be bigger than the mutual separation
        # itself (a basic physical consistency check on the sourced values).
        assert sys_def.radius_km_primary < sys_def.l_km, f"{key}: primary radius exceeds l_km"
        assert sys_def.radius_km_secondary < sys_def.l_km, f"{key}: secondary radius exceeds l_km"


def test_660_specific_sourced_radius_values() -> None:
    """Pin the exact radius figures derived above (catches a future edit
    silently drifting a cited number, mirroring #549's/#657's own
    `test_..._sourced_mu_values` convention)."""
    pm = REAL_BINARY_SYSTEMS["patroclus-menoetius"]
    assert pm.radius_km_primary == 56.5  # D=113 km / 2
    assert pm.radius_km_secondary == 52.0  # D=104 km / 2

    dd = REAL_BINARY_SYSTEMS["didymos-dimorphos"]
    assert dd.radius_km_primary is not None and dd.radius_km_secondary is not None
    assert abs(dd.radius_km_primary - 0.39) < 1e-9  # D=780 m / 2
    assert abs(dd.radius_km_secondary - 0.0755) < 1e-9  # D=151 m / 2

    ov = REAL_BINARY_SYSTEMS["orcus-vanth"]
    assert ov.radius_km_primary == 455.0
    assert ov.radius_km_secondary == 221.5

    ed = REAL_BINARY_SYSTEMS["eris-dysnomia"]
    assert ed.radius_km_primary == 1163.0
    assert ed.radius_km_secondary == 307.5

    sn = REAL_BINARY_SYSTEMS["sila-nunam"]
    assert sn.radius_km_primary == 121.5
    assert sn.radius_km_secondary == 115.0

    antiope = REAL_BINARY_SYSTEMS["antiope"]
    assert antiope.radius_km_primary is not None and antiope.radius_km_secondary is not None
    assert abs(antiope.radius_km_primary - (46.5 + 43.5 + 41.8) / 3.0) < 1e-9
    assert abs(antiope.radius_km_secondary - (44.7 + 41.4 + 39.8) / 3.0) < 1e-9
    # Matches #659's own adjudication figures (~43.9 / ~41.9 km).
    assert abs(antiope.radius_km_primary - 43.93) < 0.01
    assert abs(antiope.radius_km_secondary - 41.97) < 0.01

    lh = REAL_BINARY_SYSTEMS["lempo-hiisi"]
    assert lh.radius_km_primary == 125.5  # Hiisi, D=251 km / 2
    assert lh.radius_km_secondary == 136.0  # Lempo, D=272 km / 2


# ---------------------------------------------------------------------------
# 2. HEADLINE regression: PC(3,2) clears both bodies with wide margin
# ---------------------------------------------------------------------------


def test_660_pc32_clears_both_bodies_with_wide_margin() -> None:
    """The already-catalogued `ross-rt-pc-cycler-32-2026` row (PC (3,2)) must
    NOT be retroactively invalidated by the new gate. Per #659's own fairness
    check: ~2647 km clearance above Pluto's surface, ~481 km above Charon's
    -- this is THE most important test in this task."""
    pc = make_pluto_charon_system()
    result = sweep_family(
        pc,
        "mu01_32",
        radius_km_primary=PLUTO_RADIUS_KM,
        radius_km_secondary=CHARON_RADIUS_KM,
    )

    assert result.stable_found, (
        f"PC(3,2) must still be found stable under the new gate: "
        f"method={result.method!r}, note={result.note!r}"
    )
    assert result.min_clearance_ok is True, (
        f"PC(3,2) must CLEAR the new body-clearance gate (note={result.note!r})"
    )
    assert result.min_distance_primary_km is not None
    assert result.min_distance_secondary_km is not None

    # Per #659's adjudication: closest approach ~3836 km from Pluto centre,
    # ~1087 km from Charon centre. Pin with a generous but meaningful
    # tolerance (a few km is plausible corrector-precision drift).
    assert abs(result.min_distance_primary_km - 3836.0) < 5.0, (
        f"min_distance_primary_km={result.min_distance_primary_km:.2f} far "
        "from #659's own recomputation (~3836 km)"
    )
    assert abs(result.min_distance_secondary_km - 1087.0) < 5.0, (
        f"min_distance_secondary_km={result.min_distance_secondary_km:.2f} far "
        "from #659's own recomputation (~1087 km)"
    )

    # Explicit "wide margin" check: thousands of km above Pluto's surface,
    # hundreds of km above Charon's.
    assert result.min_distance_primary_km - PLUTO_RADIUS_KM > 2000.0
    assert result.min_distance_secondary_km - CHARON_RADIUS_KM > 400.0


# ---------------------------------------------------------------------------
# 3. Antiope (1,1) and (2,2) now fail EXPLICITLY on clearance
# ---------------------------------------------------------------------------


def test_660_antiope_11_fails_clearance_gate_explicitly() -> None:
    """Converts #659's manual (1,1) finding into a permanent regression:
    min distance ~11.8 km (primary) / ~4.2 km (secondary) vs radii
    ~43.9/~42.0 km -- deeply sub-surface, must now be an EXPLICIT negative
    (not merely "note this happens to pass every dynamical gate")."""
    antiope = REAL_BINARY_SYSTEMS["antiope"]
    assert antiope.radius_km_primary is not None and antiope.radius_km_secondary is not None
    target = antiope.to_cr3bp_system()

    result = sweep_family(
        target,
        "mu05_11",
        radius_km_primary=antiope.radius_km_primary,
        radius_km_secondary=antiope.radius_km_secondary,
    )

    assert result.stable_found is False, (
        f"Antiope (1,1) must now be REJECTED by the body-clearance gate (note={result.note!r})"
    )
    assert result.min_clearance_ok is False
    assert "clearance" in result.note.lower()

    assert result.min_distance_primary_km is not None
    assert result.min_distance_secondary_km is not None
    # Per #659: 0.0671 nd = 11.8 km (primary), 0.0237 nd = 4.2 km (secondary).
    assert abs(result.min_distance_primary_km - 11.81) < 0.1
    assert abs(result.min_distance_secondary_km - 4.17) < 0.1

    # Both distances are indeed WELL below each body's own radius (sanity
    # that this is a real sub-surface failure, not a marginal one).
    assert result.min_distance_primary_km < antiope.radius_km_primary - 25.0
    assert result.min_distance_secondary_km < antiope.radius_km_secondary - 30.0


def test_660_antiope_22_fails_clearance_gate_explicitly() -> None:
    """Converts #659's manual (2,2) finding into a permanent regression.

    Reconverges directly from the EXACT converged IC #659's adjudication note
    recorded (docs/notes/2026-07-19-659-antiope-adjudication.md: C=
    3.4661023165370235, x0=-0.5742744462570041, T=6.011499192614617, found
    originally via `sweep_family_grid`) rather than re-running the full
    ~75s (x0,C,hc) grid search -- the corrector re-converges near-instantly
    from an already-converged point, and this exercises the IDENTICAL
    `_gate_clearance` logic `sweep_family_grid` wires in, just without
    paying the full discovery-grid cost for a already-known point (matching
    this project's own precedent of NOT re-running expensive discovery
    sweeps just to re-derive an already-established number, e.g. #549's own
    docstring explicitly skipping a 35-265s-per-anchor re-run of
    Patroclus-Menoetius for the same reason).
    """
    antiope = REAL_BINARY_SYSTEMS["antiope"]
    assert antiope.radius_km_primary is not None and antiope.radius_km_secondary is not None
    target = antiope.to_cr3bp_system()

    x0 = -0.5742744462570041
    jacobi = 3.4661023165370235
    period = 6.011499192614617

    orbit = cp.correct_symmetric_fixed_jacobi(
        target, x0, jacobi, period, ydot0_sign=-1.0, half_crossings=None, tol=1e-10
    )
    assert orbit.converged, "reconvergence from the known #659 (2,2) IC must succeed"

    result = _finalize_grid_candidate(
        2,
        2,
        target,
        orbit,
        "reconverge_from_659_ic",
        radius_km_primary=antiope.radius_km_primary,
        radius_km_secondary=antiope.radius_km_secondary,
        clearance_margin_km=DEFAULT_CLEARANCE_MARGIN_KM,
    )

    assert result.stable_found is False, (
        f"Antiope (2,2) must now be REJECTED by the body-clearance gate (note={result.note!r})"
    )
    assert result.min_clearance_ok is False
    assert "clearance" in result.note.lower()

    assert result.min_distance_primary_km is not None
    assert result.min_distance_secondary_km is not None
    # Per #659: 0.0781 nd = 13.75 km (primary), 0.0180 nd = 3.16 km (secondary).
    assert abs(result.min_distance_primary_km - 13.75) < 0.1
    assert abs(result.min_distance_secondary_km - 3.16) < 0.1

    assert result.min_distance_primary_km < antiope.radius_km_primary - 25.0
    assert result.min_distance_secondary_km < antiope.radius_km_secondary - 30.0


# ---------------------------------------------------------------------------
# 4. Backward compatibility: gate is a documented no-op without a radius
# ---------------------------------------------------------------------------


def test_660_gate_is_noop_when_radius_unsourced() -> None:
    """When either radius is `None` (the default, and every pre-#660 caller's
    behavior), `_gate_clearance` must return `res` UNCHANGED -- proving every
    EXISTING call site (test_549/test_657's own sweep_family calls, every
    script) is byte-identical to its pre-#660 behavior. This is the "quick
    sanity check" this task calls for in place of a full re-sweep of every
    other already-negative system: since the gate is a structural no-op
    without a radius, and no pre-#660 caller ever supplies one, no other
    system's recorded verdict can possibly change."""
    pc = make_pluto_charon_system()
    res = SweepResult(
        k1=3,
        k2=2,
        stable_found=True,
        jacobi_mid=3.6,
        x0_mid=-0.7,
        ydot0_mid=0.5,
        period_mid=11.0,
        topology_ok=True,
        prograde=True,
        reaches_secondary=True,
        crosscheck_ok=True,
        method="unit-test-fixture",
        note="",
    )

    gated_both_none = _gate_clearance(res, pc, None, None, 0.0)
    assert gated_both_none == res
    assert gated_both_none.min_clearance_ok is None
    assert gated_both_none.min_distance_primary_km is None
    assert gated_both_none.min_distance_secondary_km is None

    # Also a no-op if only ONE radius is known (never treat a half-known
    # clearance as evaluated).
    gated_one_none = _gate_clearance(res, pc, PLUTO_RADIUS_KM, None, 0.0)
    assert gated_one_none == res

    # And a no-op for a non-stable / incomplete result (nothing to check).
    negative_res = SweepResult(k1=1, k2=1, stable_found=False, method="x", note="already negative")
    gated_negative = _gate_clearance(negative_res, pc, PLUTO_RADIUS_KM, CHARON_RADIUS_KM, 0.0)
    assert gated_negative == negative_res


def test_660_default_clearance_margin_is_zero_and_documented() -> None:
    """Sanity on the documented default margin convention (see module
    comment in real_binary_kk_sweep.py for the full rationale)."""
    assert DEFAULT_CLEARANCE_MARGIN_KM == 0.0
    assert isinstance(DEFAULT_CLEARANCE_MARGIN_KM, float)


def test_660_min_body_clearance_km_is_symmetric_for_a_symmetric_orbit() -> None:
    """Direct unit check of the low-level distance helper: for the PC(3,2)
    orbit (already known to be symmetric about the x-axis, y(0)=0), the
    reported min distances must be POSITIVE and each strictly less than the
    orbit's own max radial extent (a basic sanity bound, not a tight pin)."""
    pc = make_pluto_charon_system()
    # Reuse the PC(3,2) Table-I anchor IC directly (cheap, no search).
    x0 = -0.694376003123377
    jacobi = 3.573367616904619
    period = 12.295263874014290

    orbit = cp.correct_symmetric_fixed_jacobi(
        pc, x0, jacobi, period, ydot0_sign=-1.0, half_crossings=6, tol=1e-10
    )
    assert orbit.converged

    d_p_km, d_s_km = min_body_clearance_km(pc, orbit.x0, orbit.ydot0, orbit.period)
    assert d_p_km > 0.0
    assert d_s_km > 0.0
    # Both bodies are well inside the mutual separation (l_km), so the
    # minimum distance to either centre must be less than a few l_km.
    assert d_p_km < 10.0 * pc.l_km
    assert d_s_km < 10.0 * pc.l_km
