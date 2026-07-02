"""Tests for the generic pseudo-arclength continuation primitive (#524).

POSITIVE CONTROL for the machinery: the unit circle x^2 + y^2 - 1 = 0. This is
a closed-form, exactly-verifiable co-dimension-1 curve with a KNOWN fold in
every natural-parameter sense (fixing x and solving for y fails to turn
around at x=+-1, where dy/dx -> infinity; likewise for y at y=+-1) -- the
same structural obstruction the #249 Braik-Ross C11a/C11b/C21 fold and the
#496 cross-system phase-closure wall both hit. A working arclength walk
should sail through x=+-1 and y=+-1 without any special handling, staying on
the circle throughout. EXPECTED values here are the closed-form circle
identity itself, never a value this module computed and then asserted
against itself.
"""

from __future__ import annotations

import numpy as np

from cyclerfinder.search.pseudo_arclength import (
    ContinuationStopReason,
    compute_tangent,
    continue_curve,
)


def _circle_residual(z: np.ndarray) -> np.ndarray | None:
    x, y = z
    return np.array([x * x + y * y - 1.0])


def test_compute_tangent_orthogonal_to_gradient_on_circle() -> None:
    """The tangent at (1, 0) must be orthogonal to the gradient (2x, 2y) = (2, 0).

    Tolerance is set by the forward-difference Jacobian's own precision
    (h=1e-7 in ``_numerical_jacobian``), not by machine epsilon.
    """
    z = np.array([1.0, 0.0])
    tan = compute_tangent(_circle_residual, z)
    assert tan is not None
    assert abs(float(tan @ np.array([2.0, 0.0]))) < 1e-6
    assert abs(float(np.linalg.norm(tan)) - 1.0) < 1e-9


def test_continue_curve_stays_on_circle_through_a_fold() -> None:
    """Walking from (0, 1) must pass through the x=1 fold (dy/dx -> inf there)
    while the residual stays near zero at every corrected point -- the fold a
    natural-parameter (fix x, solve y) walk cannot turn.

    Direction is pinned explicitly with ``initial_tangent``: at (0, 1) the
    tangent's y-component is ~0 by construction (tangent = (+-1, 0) exactly
    on the closed-form circle), so the module's default "orient toward
    increasing z[-1]" heuristic is genuinely ambiguous there and must not be
    relied on -- this is a real, documented degenerate case, not a bug.
    """
    z0 = np.array([0.0, 1.0])
    curve = continue_curve(
        _circle_residual,
        z0,
        step_size=0.05,
        max_steps=40,
        tol=1e-12,
        initial_tangent=np.array([1.0, 0.0]),
    )
    assert curve.stop_reason == ContinuationStopReason.MAX_STEPS
    assert len(curve.points) == 41  # z0 plus 40 corrected steps

    zs = curve.z_values()
    residuals = zs[:, 0] ** 2 + zs[:, 1] ** 2 - 1.0
    assert np.max(np.abs(residuals)) < 1e-9, "every point must stay on the circle"

    # The walk moves toward increasing x from (0,1) (direction pinned above)
    # and must cross x=1 (the fold) -- i.e. some point has x close to 1 while
    # y crosses through ~0.
    max_x = float(np.max(zs[:, 0]))
    assert max_x > 0.99, f"walk did not reach the x=1 fold region (max x={max_x})"
    # Near the fold, y must pass close to 0 (the top-right quadrant transit).
    idx_near_fold = int(np.argmax(zs[:, 0]))
    assert abs(zs[idx_near_fold, 1]) < 0.2


def test_continue_curve_full_loop_returns_near_start() -> None:
    """A long enough walk traces the whole circle and returns near (0, 1)."""
    z0 = np.array([0.0, 1.0])
    n_steps = 400
    step_size = 2.0 * np.pi / n_steps  # arclength of a unit circle is 2*pi
    curve = continue_curve(
        _circle_residual,
        z0,
        step_size=step_size,
        max_steps=n_steps,
        tol=1e-12,
    )
    assert curve.stop_reason == ContinuationStopReason.MAX_STEPS
    final = curve.points[-1].z
    assert np.linalg.norm(final - z0) < 0.05, f"loop did not close: final={final}"
    assert abs(curve.points[-1].arclength_s - 2.0 * np.pi) < 1e-9


def test_continue_curve_target_reached_stops_at_crossing() -> None:
    """Walking toward increasing x stops once x crosses a target value."""
    z0 = np.array([0.0, 1.0])
    curve = continue_curve(
        _circle_residual,
        z0,
        step_size=0.05,
        max_steps=100,
        tol=1e-12,
        target_index=0,
        target_value=0.5,
    )
    assert curve.stop_reason == ContinuationStopReason.TARGET_REACHED
    last_two_x = [p.z[0] for p in curve.points[-2:]]
    assert min(last_two_x) <= 0.5 <= max(last_two_x)


def test_continue_curve_reports_residual_unavailable_at_bad_seed() -> None:
    """A seed the residual function can never evaluate fails cleanly, no crash."""

    def _always_none(_z: np.ndarray) -> np.ndarray | None:
        return None

    curve = continue_curve(_always_none, np.array([0.0, 1.0]), max_steps=5)
    assert curve.stop_reason == ContinuationStopReason.RESIDUAL_UNAVAILABLE
    assert curve.points == []


def test_compute_tangent_rejects_wrong_codimension() -> None:
    """A residual with the wrong number of components (not N-1) raises, not silently mis-solves."""

    def _bad_residual(z: np.ndarray) -> np.ndarray | None:
        # 2 equations for a 2-vector z -- codimension 0, not the codimension-1
        # case this module targets.
        return np.array([z[0] - 1.0, z[1] - 1.0])

    import pytest

    with pytest.raises(ValueError, match="N-1 components"):
        compute_tangent(_bad_residual, np.array([1.0, 1.0]))
