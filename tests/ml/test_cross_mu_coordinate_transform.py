"""Tests for `#649`'s coordinate-transform module
(``cyclerfinder.ml.cross_mu_coordinate_transform``).

Small constructed fixtures only, no physical/discovery claim pinned here --
see ``scripts/run_649_coordinate_fix_pilot.py`` for the actual evaluation
against `#624`'s protocol. These tests check the transform's own internal
self-consistency: identity at mu_target==mu_train, rho-invariance under the
transform (the entire point of the construction), and honest ``None`` on an
unrealizable draw.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.cr3bp import jacobi_constant
from cyclerfinder.ml.cross_mu_coordinate_transform import (
    hill_radius,
    jacobi_at_l1,
    rho_scaled_energy,
    transform_seed_to_target_mu,
)

_MU_TRAIN = 0.01215  # #628's TRAINING_MU
_MU_MU001 = 0.001
_MU_SUN_EARTH = 3.003e-6


def _sample_state0(rng: np.random.Generator) -> tuple[np.ndarray, float]:
    x0 = rng.uniform(0.5, 1.3)
    y0 = rng.uniform(-0.3, 0.3)
    z0 = rng.uniform(-0.1, 0.1)
    v = rng.uniform(-1.0, 1.0, 3)
    state0 = np.array([x0, y0, z0, *v])
    period = float(rng.uniform(1.0, 10.0))
    return state0, period


def test_jacobi_at_l1_is_above_three_and_rho_of_l1_is_exactly_one() -> None:
    for mu in (_MU_TRAIN, _MU_MU001, _MU_SUN_EARTH, 0.3):
        c_l1 = jacobi_at_l1(mu)
        assert c_l1 > 3.0
        # rho=(C-3)/(C_L1-3) evaluated AT C=C_L1 must be exactly 1 by
        # construction -- a self-consistency check of rho_scaled_energy
        # against jacobi_at_l1, not a hardcoded physical number.
        assert rho_scaled_energy(c_l1, mu) == pytest.approx(1.0, abs=1e-12)


def test_hill_radius_shrinks_with_mu() -> None:
    # #629's own established Hill-scaling direction: smaller mu -> smaller
    # natural length/energy scale around L1.
    r_train = hill_radius(_MU_TRAIN)
    r_mu001 = hill_radius(_MU_MU001)
    r_sun_earth = hill_radius(_MU_SUN_EARTH)
    assert r_train > r_mu001 > r_sun_earth > 0.0


def test_identity_transform_when_target_equals_train_mu() -> None:
    rng = np.random.default_rng(42)
    state0, period = _sample_state0(rng)
    result = transform_seed_to_target_mu(state0, period, _MU_TRAIN, _MU_TRAIN)
    assert result is not None
    np.testing.assert_allclose(result.state0, state0)
    assert result.period == period
    assert result.scale == pytest.approx(1.0)
    assert result.c_guess == pytest.approx(result.c_target)
    # rho must match direct computation from the SAME state/mu.
    c_direct = jacobi_constant(state0, _MU_TRAIN)
    expected_rho = rho_scaled_energy(c_direct, _MU_TRAIN)
    assert result.rho == pytest.approx(expected_rho)


@pytest.mark.parametrize("mu_target", [_MU_MU001, _MU_SUN_EARTH])
def test_transformed_seed_reproduces_the_target_rho_exactly(mu_target: float) -> None:
    """The entire point of the construction (task step 3): whenever the
    transform succeeds (returns non-None), the CONSTRUCTED state's own
    Jacobi constant at mu_target, run back through rho_scaled_energy, must
    reproduce the SAME rho the raw Earth-Moon guess had -- not merely close,
    exactly (up to floating point), because the velocity-magnitude solve is
    designed to hit c_target precisely.
    """
    rng = np.random.default_rng(7)
    n_checked = 0
    for _ in range(50):
        state0_guess, period_guess = _sample_state0(rng)
        result = transform_seed_to_target_mu(state0_guess, period_guess, _MU_TRAIN, mu_target)
        if result is None:
            continue
        n_checked += 1
        c_reconstructed = jacobi_constant(result.state0, mu_target)
        assert c_reconstructed == pytest.approx(result.c_target, abs=1e-9)
        rho_reconstructed = rho_scaled_energy(c_reconstructed, mu_target)
        assert rho_reconstructed == pytest.approx(result.rho, abs=1e-8)
    # sanity: the fixed rng/bounding-box combination above must actually
    # exercise the success path at least once, or this test would pass
    # vacuously.
    assert n_checked > 0


def test_transform_returns_none_on_an_unrealizable_draw() -> None:
    """A concrete, reproducible draw (found by scanning a fixed rng stream)
    where the Hill-scaled position at mu=0.001 sits outside the target
    system's zero-velocity surface at the rho-matched target Jacobi constant
    -- the transform must return None honestly rather than fabricate a
    complex/negative-energy velocity.
    """
    state0_guess = np.array(
        [0.53393922, 0.04894849, -0.01515395, 0.31708549, 0.06292507, -0.16643172]
    )
    period_guess = 4.168233774614348
    result = transform_seed_to_target_mu(state0_guess, period_guess, _MU_TRAIN, _MU_MU001)
    assert result is None


def test_transform_returns_none_on_zero_velocity_guess() -> None:
    """No velocity direction to preserve -- a separate honest-None edge case
    from the unrealizable-energy one above."""
    state0_guess = np.array([0.8, 0.0, 0.0, 0.0, 0.0, 0.0])
    result = transform_seed_to_target_mu(state0_guess, 5.0, _MU_TRAIN, _MU_MU001)
    assert result is None


def test_rejects_non_6vector_state0() -> None:
    with pytest.raises(ValueError):
        transform_seed_to_target_mu([1.0, 2.0, 3.0], 5.0, _MU_TRAIN, _MU_MU001)
