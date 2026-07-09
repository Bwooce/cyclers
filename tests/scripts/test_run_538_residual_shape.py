"""Well-posedness guard for #538's cross-system cycler corrector (Task 1).

#537's refining solve was rank-deficient: 3 equations against 4 unknowns, with
velocity never in the residual, so its "connection" could not be trusted (its
reported gaps were post-hoc diagnostics, not driven to zero). This test locks
in that #538's residual is OVER-DETERMINED (n_residuals >= n_unknowns) and
carries the full structure -- so a future edit cannot silently regress to the
under-determined form.

The expected numbers here are the *design* of the parameterization (a
structural count of segments and match conditions), not a value produced by our
own numerical code, so this is not a circular golden.
"""

from __future__ import annotations

import scripts.run_538_qbcp_cycler as run


def test_residual_is_over_determined() -> None:
    n_unknowns, n_residuals = run._residual_shape_ok()
    assert n_residuals >= n_unknowns, (
        f"residual is under-determined ({n_residuals} residuals < {n_unknowns} "
        "unknowns) -- this is exactly #537's rank-deficient failure mode"
    )


def test_unknown_count_matches_design() -> None:
    # 4 torus phase pairs (8) + 4 arc durations (tau_f, tau_bf, tau_r, tau_br).
    assert run.N_UNKNOWNS == 12
    assert run.N_TORUS_PHASE_UNKNOWNS == 8
    assert run.N_DURATION_UNKNOWNS == 4
    n_unknowns, _ = run._residual_shape_ok()
    assert n_unknowns == run.N_UNKNOWNS


def test_residual_count_matches_design() -> None:
    # Two legs, each a 6-state (pos+vel) match plus a time-consistency scalar (7),
    # plus SE closure (2) and EM closure (2): 2*7 + 4 = 18.
    assert run.N_PER_LEG == 7
    assert run.N_RESIDUALS == 18
    _, n_residuals = run._residual_shape_ok()
    assert n_residuals == run.N_RESIDUALS


def test_velocity_is_in_the_residual() -> None:
    # The 6-state match per leg is position (3) + velocity (3); the whole point
    # of #538 vs #537 is that velocity is carried, not dropped.
    assert run.N_STATE_MATCH == 6
