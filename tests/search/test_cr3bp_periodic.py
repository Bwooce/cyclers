"""CR3BP periodic-orbit corrector (plan 2026-06-10, Phase 2)."""

from __future__ import annotations

import numpy as np
import pytest

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp

# Sourced Arenstorf golden (μ, x0, vy0, period) — Arenstorf 1963 / Hairer et al.
# (Hairer, Nørsett, Wanner, "Solving ODEs I", p. 129, test problem B5)
MU = 0.012277471
X0, VY0, PERIOD = 0.994, -2.0015851063790825, 17.0652165601579625


def test_arenstorf_orbit_is_periodic() -> None:
    # The published Arenstorf IC already (very nearly) closes after one period —
    # the corrector confirms it and returns a tight closure residual.
    sysm = cr3bp.CR3BPSystem(mu=MU, primary="t", secondary="t", l_km=1.0, t_s=1.0)
    s0 = np.array([X0, 0.0, 0.0, 0.0, VY0, 0.0])
    res = cp.correct_periodic(sysm, s0, PERIOD)
    assert res.converged
    assert res.closure_residual < 1e-6
    assert res.period == pytest.approx(PERIOD, rel=1e-2)


def test_non_periodic_guess_does_not_converge() -> None:
    sysm = cr3bp.CR3BPSystem(mu=MU, primary="t", secondary="t", l_km=1.0, t_s=1.0)
    s0 = np.array([0.7, 0.0, 0.0, 0.0, -0.2, 0.0])
    res = cp.correct_periodic(sysm, s0, 3.0, max_iter=8)
    assert not res.converged
