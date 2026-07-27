"""Tests for `#729`'s epoch(+torus-point) robustness scan
(``scripts/run_729_epoch_torus_robustness_scan.py``).

`#729` reuses `#726`'s own module
(:mod:`cyclerfinder.search.crnbp_real_ephemeris_consistency`) UNMODIFIED-in-
algorithm as a library (this test suite's own build DID fix one real bug in
that module -- an LSK ``spice.furnsh`` call with no idempotency guard that
exhausted CSPICE's ``MAXFIL`` kernel-pool limit after ~5300 repeated calls;
see the module's own ``_LSK_FURNISHED`` guard and comment. That fix does not
change any computed physics, only whether a multi-thousand-call scan can run
to completion at all). The physics/propagation pipeline itself is already
covered by `#726`'s own test suite
(``tests/search/test_crnbp_real_ephemeris_consistency.py``). This test file
covers only the NEW orchestration/aggregation logic `#729` adds: running
`#726`'s own dense-synodic-scan-plus-bisection methodology at several epochs
AND several torus points, and aggregating duty-cycle/local-minimum results
across them -- mirroring `#705`'s own test file's scope and style exactly
(``tests/scripts/test_705_epoch_robustness_scan.py``).

SPICE-gated (skipped if ``jup365.bsp``/the vendored LSK are not installed),
matching `#705`'s/`#726`'s own test modules' skip convention.
"""

from __future__ import annotations

import datetime as _dt
import math

import numpy as np
import pytest

import scripts.run_729_epoch_torus_robustness_scan as run
from cyclerfinder.verify.spice_kernels import NAIF_JUP365_LOCAL, VENDORED_LSK_PATH

_KERNELS_PRESENT = NAIF_JUP365_LOCAL.expanduser().exists() and VENDORED_LSK_PATH.exists()
_SKIP_REASON = (
    f"jup365.bsp not installed (looked at {NAIF_JUP365_LOCAL.expanduser()}); "
    "see cyclerfinder.verify.spice_kernels.ensure_jup365_kernel"
)


def test_duty_cycle_window_counting_matches_705_own_logic() -> None:
    """`#729`'s own :func:`_duty_cycle` must count distinct sub-threshold
    windows the same way `#705`'s own driver script did (contiguous runs of
    ``True`` in the boolean threshold mask, non-circular) -- this is the
    IDENTICAL test case `#705`'s own test file used, reproduced here to lock
    in that the duplicated logic (see module docstring: duplicated for
    import-independence between sibling scan scripts, not imported) behaves
    identically."""
    pos = np.array([60000.0, 800.0, 800.0, 60000.0, 60000.0, 1500.0, 60000.0])
    duty = run._duty_cycle(pos)
    below_1000 = duty["below_1000km"]
    assert below_1000["n_distinct_windows"] == 1
    assert below_1000["fraction_of_period"] == pytest.approx(2.0 / 7.0)
    below_2000 = duty["below_2000km"]
    assert below_2000["n_distinct_windows"] == 2
    assert below_2000["fraction_of_period"] == pytest.approx(3.0 / 7.0)


def test_narrow_near_miss_threshold_is_sourced_from_705_own_scale() -> None:
    """The 'narrow near-miss' gate (``NARROW_NEAR_MISS_KM``) is explicitly
    borrowed from `#705`'s own established scale (`#726` has no own
    committed near-miss reference to derive a threshold from the way `#705`
    derived one from `#704`'s committed best point) -- 5,000 km, matching
    `#704`'s own original "3 distinct sub-5000km-mismatch windows"
    characterization and comfortably above `#705`'s own worst per-epoch
    local minimum (142.82 km)."""
    assert pytest.approx(5000.0) == run.NARROW_NEAR_MISS_KM


def test_torus_point_axes_are_the_documented_symmetric_sample() -> None:
    """The headline point matches `#726`'s own headline choice
    (``theta1=theta2=0``); the 4 secondary points are the documented
    ``{0, pi} x {0, pi}`` corners plus the diagonal midpoint -- a
    deliberately systematic (not ad hoc) sample per the module docstring's
    own justification for adding this second axis beyond `#705`'s scope."""
    assert run.HEADLINE_TORUS_POINT == (0.0, 0.0)
    expected_secondary = {
        (np.pi, 0.0),
        (0.0, np.pi),
        (np.pi, np.pi),
        (np.pi / 2.0, np.pi / 2.0),
    }
    assert set(run.SECONDARY_TORUS_POINTS) == expected_secondary
    assert len(run.SECONDARY_TORUS_POINTS) == 4


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_lsk_furnish_guard_prevents_kernel_pool_exhaustion() -> None:
    """Regression test for the real bug `#729`'s own build found and fixed
    in `#726`'s module (:mod:`cyclerfinder.search.crnbp_real_ephemeris_consistency`):
    calling ``check_torus_survives_real_ephemeris`` many times must NOT grow
    SPICE's loaded-kernel count (``spice.ktotal('ALL')``) past a small,
    constant number -- before the fix, each call re-furnished the LSK
    (CSPICE does not dedupe identical-path ``furnsh`` calls), and ~5300
    calls raised ``SpiceNOMOREROOM`` (the ``MAXFIL`` limit), which is
    trivially reached by any real multi-epoch scan (`#729`'s own primary
    axis alone issues ~3850 calls).
    """
    import spiceypy as spice

    import cyclerfinder.core.ccr4bp as ccr4bp
    import cyclerfinder.core.crnbp as crnbp
    import cyclerfinder.search.crnbp_real_ephemeris_consistency as rec5
    import cyclerfinder.search.variational_ccr4bp_torus as vt
    import cyclerfinder.search.variational_crnbp_torus as vc

    system4 = ccr4bp.jupiter_europa_ganymede_default()
    s0, period, res = run._resonant_symmetric_orbit(system4.mu, 3, 4)
    assert res < 1e-10
    ccr4bp_torus = vt.discover_ccr4bp_torus_from_resonant_orbit(
        system4,
        s0,
        period,
        n1=2,
        n2=10,
        tr_solver="exact",
        max_nfev=300,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    target = crnbp.jupiter_europa_io_ganymede_default()
    torus = vc.discover_crnbp_torus_from_ccr4bp_seed(
        ccr4bp_torus,
        mu_io=0.0,
        a_io=target.perturbers[0].a,
        omega_io=target.perturbers[0].omega,
        theta_io0=target.perturbers[0].theta0,
        tr_solver="exact",
        max_nfev=300,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    # The FIRST call legitimately furnishes both the LSK and jup365.bsp (each
    # module's own lazy, module-level "furnish once" guard) -- so the count
    # after call 1 is the correct baseline, not the pre-call-0 count (which
    # may be 0 in a fresh worker process or already-elevated if another test
    # in this session furnished kernels first).
    first = rec5.check_torus_survives_real_ephemeris(
        "2000-01-01T00:00:00", torus, 0.0, 0.0, t_window_tu=torus.period * 0.02
    )
    assert first.propagation_success
    baseline = spice.ktotal("ALL")
    for _ in range(29):
        result = rec5.check_torus_survives_real_ephemeris(
            "2000-01-01T00:00:00", torus, 0.0, 0.0, t_window_tu=torus.period * 0.02
        )
        assert result.propagation_success
    after = spice.ktotal("ALL")
    assert after == baseline, (
        f"SPICE kernel count grew from {baseline} to {after} over 29 further calls -- "
        "the LSK furnish-idempotency guard regressed."
    )


@pytest.mark.skipif(not _KERNELS_PRESENT, reason=_SKIP_REASON)
def test_epoch_torus_robustness_scan_runs_end_to_end_and_is_finite() -> None:
    """End-to-end smoke test at reduced scale (2 epochs x 6 synodic points,
    2 torus points -- not `#729`'s own full-resolution production run)
    against the ACTUAL `#720`/`#723`/`#724` torus, deterministically
    reconstructed the same way `#726`'s own test module does. Does not
    assert narrowness/recurrence (that is `#729`'s own reported finding, not
    a frozen regression target) -- only that the orchestration runs end to
    end and produces finite, physically-scaled numbers, matching `#705`'s
    own e2e smoke test's spirit for the single-axis case.
    """
    t0 = 0.0
    torus = run.rebuild_724_final_torus(t0)
    assert torus.residual_rms < 1e-2

    base_dts = [
        _dt.datetime(2000, 1, 1, tzinfo=_dt.UTC),
        _dt.datetime(2083, 1, 1, tzinfo=_dt.UTC),
    ]
    synodic_period_days = 7.05  # ballpark check only, not a frozen target
    torus_points = [(0.0, 0.0), (np.pi, np.pi)]

    scans = []
    for th1, th2 in torus_points:
        for base_dt in base_dts:
            scan = run.scan_one_epoch_one_torus_point(
                torus, th1, th2, base_dt, synodic_period_days, 6, t0, label="smoke"
            )
            scans.append(scan)

    assert len(scans) == 4
    for scan in scans:
        assert math.isfinite(scan["synodic_scan_pos_gap_km_min"])
        assert math.isfinite(scan["synodic_scan_pos_gap_km_max"])
        assert math.isfinite(scan["synodic_scan_pos_gap_km_median"])
        lm = scan["local_minimum"]
        assert lm["propagation_success"]
        assert math.isfinite(lm["pos_gap_km"])
        assert math.isfinite(lm["vel_gap_km_s"])
        # Sanity scale: bounded well above Ganymede's own SMA (~1.07e6 km)
        # but not required to be small -- a wiring-bug blowup check, not a
        # narrowness assertion (matching #705's own e2e test's scale check).
        assert lm["pos_gap_km"] < 5.0e6
        assert len(scan["synodic_phase_scan"]) == 6
        assert isinstance(scan["narrow_near_miss"], bool)
