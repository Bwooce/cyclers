"""Pseudo-arclength continuation-in-e (ER3BP) — smooth-family regression golden.

On a SMOOTH family that does not fold below ``e_target`` (the Broucke 1969
Earth-Moon Family-7P floor seed used by #432 — it continues smoothly to
e=0.0549), the new fold-capable pseudo-arclength walker must agree with the
trusted secant ``continue_er3bp_family_in_e``: same family, no fold => both
land on the same final member. This pins the arclength machinery's correctness
against the secant on a non-folding case (the fold golden lives in a separate
slow-marked test).
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.genome.er3bp_continuation import (
    continue_er3bp_family_in_e,
    continue_er3bp_family_in_e_arclength,
)

# Broucke 1969 TR 32-1360 Table 12, Family 7P, Earth-Moon mu=0.0121550, Orbit 1.
# The e=0.0001 member IC; #432 continues this smoothly to the e=0.0549 member.
_MU_EM = 0.0121550
_BROUCKE_EM_IC = np.array([0.1520965, 0.0, 0.0, 0.0, 3.1608994, 0.0])
# Corrector ``period_f`` is the INTEGRATION span = the half-period (pi) when
# is_half_period_residual=True (the full 2*pi period is stored elsewhere); this
# matches how er3bp_discovery.continue_and_monitor drives the secant.
_PERIOD_F = np.pi
_E_TARGET = 0.0549


@pytest.mark.slow
def test_arclength_matches_secant_on_smooth_family() -> None:
    """Arclength walk reaches e_target and matches the secant final member.

    The Broucke EM floor family does not fold below 0.0549, so the
    pseudo-arclength continuator (which CAN walk through folds) must reproduce
    the secant continuator's endpoint on this smooth family.
    """
    sys_base = ER3BPSystem(
        mu=_MU_EM,
        e=0.0,
        primary_name="Earth",
        secondary_name="Moon",
    )

    secant = continue_er3bp_family_in_e(
        sys_base,
        _BROUCKE_EM_IC,
        _PERIOD_F,
        _E_TARGET,
        n_steps=40,
        is_half_period_residual=True,
        tol=1e-10,
    )
    secant_final = secant[-1]

    arclength = continue_er3bp_family_in_e_arclength(
        sys_base,
        _BROUCKE_EM_IC,
        _PERIOD_F,
        _E_TARGET,
        ds=0.005,
        max_steps=400,
        is_half_period_residual=True,
        tol=1e-10,
    )

    # Walk reached the target.
    assert len(arclength) >= 2
    final_e = arclength[-1].e
    assert final_e == pytest.approx(_E_TARGET, abs=1e-3), (
        f"arclength stopped at e={final_e}, expected to reach {_E_TARGET}"
    )

    # Same family, no fold => same final member as the secant (loose tol:
    # the two predictors take different paths/step counts to the same orbit).
    np.testing.assert_allclose(
        arclength[-1].state0,
        secant_final.state0,
        atol=1e-6,
        err_msg="arclength final IC diverged from secant on a smooth family",
    )
