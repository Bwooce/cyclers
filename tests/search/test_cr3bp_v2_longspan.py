"""Pinning test for the #229 long-span discriminating inertial run (V2 evidence).

Re-runs the `ross-v2-longspan-v1` instrument (scripts/ross_v2_longspan.py,
which reuses the campaign's `cr3bp-inertial-rebound-ias15-v1` REBOUND/IAS15
harness verbatim) on ONE representative row -- ross-rt-em-cycler-11-2025, the
(1,1) family, the fastest of the five -- at a shortened-but-still-
DISCRIMINATING span, and pins the bounded-band verdict + Jacobi-drift bound
that back the proposed V1 -> V2-ballistic promotion (held for review, #229).

DISCRIMINATION (the whole point, per the #224 assessment note): the campaign's
5-period gate cannot fail for a residual-zero IC (seed error delta_1 ~ 1e-9
cannot amplify past the 3A band in 5 periods even for |lambda| ~ 10^2..10^3).
The span here is chosen so a hypothetical |lambda| = 2 instability WOULD have
amplified the actual measured delta_1 past 3A -- the test asserts that
in-instrument (span >= N_req = ln(3A/delta_1)/ln 2), so the BOUNDED verdict is
non-vacuous by construction, not by configuration.

MODEL SCOPE (like-for-like): the harness integrates the CR3BP idealisation in
inertial coordinates -- the row's DEFINING model. No real-ephemeris claim.

Runtime: ~2-4 s of integration (40 periods x 128 IAS15 samples) + the
corrector; marked slow alongside the project's other multi-second
integration-heavy tests. Skips cleanly when rebound is absent (project rule).
"""

from __future__ import annotations

import math

import pytest

pytest.importorskip("rebound")

import scripts.cr3bp_family_search as fs
from scripts.ross_v2_longspan import LAMBDA_HYPO, run_row

# 40 periods: comfortably above the measured N_req ~ 30 for the (1,1) row
# (delta_1 ~ 5e-9, A ~ 1.8 nd) while keeping the test a few seconds.
N_PERIODS_TEST = 40


@pytest.mark.slow
def test_ross_11_longspan_bounded_band() -> None:
    seed = next(s for s in fs.ROSS_SEEDS if s.label == "ross-(1,1)")
    res = run_row(seed, N_PERIODS_TEST)

    # The instrument must have teeth: the span actually used exceeds the span a
    # hypothetical |lambda| = LAMBDA_HYPO = 2 instability needs to blow the
    # MEASURED seed error past the 3A departure band. If this fails, the run
    # was NOT discriminating and the BOUNDED verdict would be vacuous.
    assert math.isfinite(res.n_periods_required)
    assert res.n_periods >= res.n_periods_required, (
        f"span {res.n_periods} periods < discriminating threshold "
        f"{res.n_periods_required:.1f} (delta_1={res.delta1_nd:.3e}, "
        f"band={res.departure_band_nd:.3f}, lambda_hypo={LAMBDA_HYPO})"
    )

    # Bounded-band verdict: never departs the 3A band, no exponential growth
    # between span halves (a real |lambda|=2 would give a half-ratio ~2^20 at
    # this span; the observed linear phase drift gives ~2).
    assert res.verdict == "BOUNDED"
    assert not res.diverged
    assert res.max_drift_span_nd <= res.departure_band_nd
    assert res.max_drift_span_nd < 1e-4  # 4+ orders inside the band in practice

    # Jacobi drift over the span: in-model invariant, IAS15 hygiene. The full
    # 100-period run measured ~1e-10; 1e-9 is the campaign's R2 bound and holds
    # comfortably at this shorter span.
    assert res.jacobi_drift_span <= 1e-9

    # The row's stability verdict input: Barden nu reproduces inside Ross's
    # published stable window (|nu| < 1, nu ~ 0 midpoint -- same gate as the
    # #212 reproduction tests).
    assert abs(res.nu) < 1.0
    assert abs(res.nu) < 0.2
