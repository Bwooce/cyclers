"""Evidence test for task #505: Pluto-Charon (3,2) V2-ballistic long-span run.

Pins the bounded-drift verdict that backs the proposed V1 -> V2-ballistic
promotion for ``ross-rt-pc-cycler-32-2026`` (held for user adjudication).

Architecture mirrors ``tests/search/test_cr3bp_v2_longspan.py`` (the EM V2
evidence pinning test, #229) but for the Pluto-Charon system.  The harness is
identical: REBOUND IAS15 in the CR3BP inertial frame (``cr3bp-inertial-rebound-
ias15-v1``), the row's defining model (CR3BP), NOT a real-ephemeris claim.

DISCRIMINATION (the whole point, per the #224/229 assessment pattern): the
span is chosen so a hypothetical |lambda|=2 instability WOULD have amplified
the measured delta_1 past the 3A departure band.  The test asserts this
criterion holds (span >= N_req = ln(3A/delta_1)/ln 2), so the BOUNDED verdict
is non-vacuous.

NOT marked slow: this test runs in ~0.6 s (50 periods x 128 IAS15 samples +
corrector).  Per project rule (feedback_delegation_fresh_agent_not_fork), a
V-claim's evidence test must run in the DEFAULT CI suite.  Skips cleanly when
rebound is absent.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

pytest.importorskip("rebound")

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from scripts.pc_v2_longspan import HALF_RATIO_MAX, LAMBDA_HYPO, run_pc_longspan

# 50 periods: comfortably above the measured N_req ~ 33.7 for this row
# (delta_1 ~ 3.5e-10, A ~ 1.64 nd) while keeping the test under ~1 s.
N_PERIODS_TEST = 50


def test_505_pc_32_longspan_bounded_band() -> None:
    """PC (3,2) long-span bounded-drift — V2-ballistic evidence (HELD for review).

    Positive control: corrector reproduces the admitted IC (period ~11.83 TU,
    nu ~ 0, winding (3,2)).  Then 50-period REBOUND/IAS15 inertial run confirms
    bounded rotating-frame drift — the V2-ballistic discriminating criterion.
    """
    res = run_pc_longspan(N_PERIODS_TEST)

    # Span must be discriminating: >= N_req periods under hypothetical |lambda|=2
    assert math.isfinite(res.n_periods_required), (
        f"N_req is not finite (delta1={res.delta1_nd:.3e}, band={res.departure_band_nd:.3f})"
    )
    assert res.n_periods >= res.n_periods_required, (
        f"span {res.n_periods} periods < discriminating threshold "
        f"{res.n_periods_required:.1f} (delta1={res.delta1_nd:.3e}, "
        f"band={res.departure_band_nd:.3f}, lambda_hypo={LAMBDA_HYPO})"
    )

    # Bounded-band verdict: no divergence, all samples inside 3A, no exponential
    # growth between span halves.
    assert res.verdict == "BOUNDED", (
        f"Expected BOUNDED, got {res.verdict}: "
        f"max_drift={res.max_drift_span_nd:.3e}, band={res.departure_band_nd:.3f}, "
        f"half_ratio={res.half_ratio:.2f} (limit {HALF_RATIO_MAX})"
    )
    assert not res.diverged, "orbit diverged during propagation"
    assert res.max_drift_span_nd <= res.departure_band_nd, (
        f"max drift {res.max_drift_span_nd:.3e} exceeds departure band {res.departure_band_nd:.3f}"
    )
    # Drift is orders of magnitude inside the band (not just barely inside)
    assert res.max_drift_span_nd < 1e-5, (
        f"drift {res.max_drift_span_nd:.3e} >= 1e-5 (expected deep inside band)"
    )

    # Jacobi conservation: IAS15 hygiene, in-model invariant
    assert res.jacobi_drift_span <= 1e-9, (
        f"Jacobi drift {res.jacobi_drift_span:.3e} > 1e-9 over {N_PERIODS_TEST} periods"
    )

    # Positive control: Barden nu reproduces near 0 (maximally stable, catalogue value nu~1e-6)
    assert abs(res.nu) < 1.0, f"Barden stability FAILED: |nu|={abs(res.nu):.4f} >= 1"
    assert abs(res.nu) < 1e-3, f"nu={res.nu:.3e} is not near 0 — expected maximally-stable midpoint"
