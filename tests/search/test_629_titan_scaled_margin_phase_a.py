"""Task #629 Phase A: regression test for the scaled-margin (c_margin_alpha) logic.

:func:`real_binary_kk_sweep.mu_step_to_system_tracking_c_l1` gained a new
``c_margin_alpha`` parameter (#629 Phase A) that, when given, overrides
#627's absolute ``c_margin`` with one scaled to each step's own corridor
width: ``c_margin_step = alpha * (C_L1(mu_next) - 3)``. This holds the
scaled energy ``rho = (C - 3) / (C_L1(mu) - 3)`` at ``1 - alpha`` once the
walk starts clamping, rather than letting an absolute margin ratchet the
walk arbitrarily far below the shrinking corridor (see the #629 design read,
docs/notes/2026-07-18-629-design-read-titan-kk-grid.md, and the module
docstring for the full derivation).

This is a fast, deterministic self-consistency check -- not a golden test
(no published expected value exists for this parameter): a large, deliberately
aggressive alpha (0.5, well below both sourced anchors' own rho ~0.79-0.97)
must force the walk to clamp, and the LANDED member's own rho must sit at the
predicted ceiling ``1 - alpha`` to good precision, confirming the margin
actually scales with C_L1(mu) - 3 at each step rather than behaving like a
fixed absolute constant (which would land at a mu-dependent, non-constant
rho instead -- exactly #627's bug).
"""

from __future__ import annotations

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.pluto_charon_kk_sweep import _c_l1
from cyclerfinder.search.real_binary_kk_sweep import mu_step_to_system_tracking_c_l1


def test_c_margin_alpha_lands_at_predicted_ceiling_rho() -> None:
    """A short hop with a large c_margin_alpha lands at rho=1-alpha, not the anchor's own rho."""
    target = cr3bp.CR3BPSystem(mu=0.0115, primary="P1", secondary="P2", l_km=1.0, t_s=1.0)
    anchor_mu = 0.012150584270572
    anchor_x0 = -0.322477620583087
    anchor_jacobi = 3.183379082910527
    anchor_period = 19.503763587070285

    c_l1_anchor = _c_l1(anchor_mu)
    rho_anchor = (anchor_jacobi - 3.0) / (c_l1_anchor - 3.0)
    alpha = 0.5
    assert rho_anchor > 1.0 - alpha, (
        "test fixture assumption: the anchor's own rho must exceed the clamp ceiling "
        "1-alpha, otherwise this test does not actually exercise the clamp"
    )

    landed = mu_step_to_system_tracking_c_l1(
        anchor_mu,
        target,
        anchor_x0,
        anchor_jacobi,
        anchor_period,
        hc=7,
        sign=-1.0,
        n_steps=10,
        c_margin_alpha=alpha,
        tol=1e-10,
    )
    assert landed is not None, "short, well-behaved mu-decrease unexpectedly failed to converge"

    c_l1_target = _c_l1(target.mu)
    rho_landed = (landed.jacobi - 3.0) / (c_l1_target - 3.0)
    assert abs(rho_landed - (1.0 - alpha)) < 5e-3, (
        f"rho_landed={rho_landed:.5f} should track the predicted clamp ceiling "
        f"1-alpha={1.0 - alpha:.5f}, confirming the margin scales with C_L1(mu)-3"
    )


def test_c_margin_alpha_none_preserves_legacy_absolute_c_margin_behavior() -> None:
    """c_margin_alpha=None (the default) must reproduce #627's original absolute-margin result.

    Same short hop as tests/search/test_627_titan_pilot.py's own
    ``test_627_c_tracking_short_hop_is_self_consistent`` -- a pure
    non-regression check that adding the new parameter did not change
    default behavior at all.
    """
    target = cr3bp.CR3BPSystem(mu=0.0115, primary="P1", secondary="P2", l_km=1.0, t_s=1.0)
    anchor_x0 = -0.322477620583087
    anchor_jacobi = 3.183379082910527
    anchor_period = 19.503763587070285

    landed = mu_step_to_system_tracking_c_l1(
        0.012150584270572,
        target,
        anchor_x0,
        anchor_jacobi,
        anchor_period,
        hc=7,
        sign=-1.0,
        n_steps=10,
        c_margin=0.005,
        tol=1e-10,
    )
    assert landed is not None
    assert landed.converged
    assert landed.crossing_residual < 1e-8
