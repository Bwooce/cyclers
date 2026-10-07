"""#970: Schwaniger 1963 (NASA TN D-1833 sec. III.E) retrograde Earth-Moon periodic free return
as a positive control for CR3BP correctors and propagators, including any future corrector that
is regularised at both primaries (#948 R4).

The orbit passes 177 km above the Earth (perigee 6555 km, Schwaniger's injection radius) and
about 465 km above the Moon once per period. The initial state below is the project
corrector's member (``scripts/run_970_schwaniger_control.py``, data/970_schwaniger/control.json),
so it is NOT used as an expected value. The expected values are:

* sourced: Schwaniger's own definition (perpendicular perigee crossing at r = 6555 km, p.2 and
  p.6) and his graph-read "about 2150 km" periselenum and "about 650 hours" period (p.6-7),
  checked at a loose tolerance (the printed values are "about" figures);
* structural: closure, perpendicular crossings, det M = 1, the trivial eigenvalue pair and the
  Jacobi constant, all model identities independent of this code;
* cross-code: the digest's independent script (2202.5 km, 625.5 h at mu = 0.012150), labelled
  as a second-code check, not a published golden.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.core.cr3bp import CR3BPSystem, jacobi_constant, propagate
from cyclerfinder.core.cr3bp_ks import MoonCentredCR3BP, propagate_ks
from cyclerfinder.core.cr3bp_regularized import (
    physical_to_regularized_span,
    propagate_regularized,
)
from cyclerfinder.search.cr3bp_periodic import SymmetricOrbit, barden_stability

MU = 0.012150  # the digest's mass ratio (Schwaniger states none)
L_KM = 384400.0
T_S = math.sqrt(L_KM**3 / (398600.4 / (1.0 - MU)))
SYSTEM = CR3BPSystem(mu=MU, primary="Earth", secondary="Moon", l_km=L_KM, t_s=T_S)
# Corrected member at the periselenum (x-axis, Earth side of the Moon), mu = 0.012150.
X0 = 0.9821202024620871
YDOT0 = -2.4712785397261854
T_HALF = 3.000751298161069  # half period (perigee crossing)
STATE0 = np.array([X0, 0.0, 0.0, 0.0, YDOT0, 0.0])


def _r1_km(s: np.ndarray) -> float:
    return math.hypot(s[0] + MU, s[1]) * L_KM


def _r2_km(s: np.ndarray) -> float:
    return math.hypot(s[0] - 1.0 + MU, s[1]) * L_KM


def test_perigee_crossing_is_perpendicular_at_schwaniger_radius() -> None:
    half = propagate(SYSTEM, STATE0, T_HALF).state_f
    assert abs(half[1]) < 1e-7  # on the Earth-Moon line (speed ~10.7 units x event-time error)
    assert abs(half[3]) < 1e-7  # perpendicular (Schwaniger's horizontal perigee)
    assert _r1_km(half) == pytest.approx(6555.0, abs=1e-3)  # 100 n.mi. altitude, p.2
    assert half[0] + MU > 0.0  # between the Earth and the Moon (longitude 180 deg in MEP)


def test_sourced_about_values() -> None:
    periselene = _r2_km(STATE0)
    period_h = 2.0 * T_HALF * T_S / 3600.0
    assert periselene == pytest.approx(2150.0, rel=0.05)  # "about 2150 km"
    assert period_h == pytest.approx(650.0, rel=0.05)  # "about 650 hours"
    # second-code check (digest script, same model): 2202.5 km, 625.5 h
    assert periselene == pytest.approx(2202.5, abs=0.1)
    assert period_h == pytest.approx(625.5, abs=0.05)


def test_full_period_closure_three_integrators() -> None:
    period = 2.0 * T_HALF
    plain = propagate(SYSTEM, STATE0, period).state_f
    assert np.linalg.norm(plain[:3] - STATE0[:3]) * L_KM < 1e-3  # km
    ks = propagate_ks(MoonCentredCR3BP(MU), STATE0, period)
    assert np.linalg.norm(ks.state[:3] - STATE0[:3]) * L_KM < 1e-3
    assert ks.r_min * L_KM == pytest.approx(_r2_km(STATE0), abs=1e-3)
    span = physical_to_regularized_span(SYSTEM, STATE0, (0.0, period))
    sund = propagate_regularized(
        SYSTEM, STATE0, (span[0], 1.5 * span[1]), regularization="r1r2", t_stop=period
    )
    assert np.linalg.norm(sund.state_at_s[:3, -1] - STATE0[:3]) * L_KM < 1e-3
    assert jacobi_constant(plain, MU) == pytest.approx(jacobi_constant(STATE0, MU), abs=1e-10)


def test_monodromy_identities_and_instability() -> None:
    arc = propagate(SYSTEM, STATE0, 2.0 * T_HALF, with_stm=True, stm_mode="fixed_path")
    assert arc.stm is not None
    idx = [0, 1, 3, 4]
    m4 = arc.stm[np.ix_(idx, idx)]
    assert np.linalg.det(arc.stm) == pytest.approx(1.0, abs=1e-6)
    eig = np.linalg.eigvals(m4)
    near_one = sorted(eig, key=lambda e: abs(e - 1.0))[:2]
    assert all(abs(e - 1.0) < 2e-3 for e in near_one)  # trivial pair
    b_h = float(np.trace(m4) - 2.0)
    assert b_h > 2.0  # unstable in the plane
    orb = SymmetricOrbit(
        x0=X0,
        ydot0=YDOT0,
        jacobi=jacobi_constant(STATE0, MU),
        t_half=T_HALF,
        period=2.0 * T_HALF,
        converged=True,
        crossing_residual=0.0,
        n_iter=0,
    )
    nu, _ = barden_stability(SYSTEM, orb)
    assert 2.0 * nu == pytest.approx(b_h, rel=1e-6)  # half-period (Barden) form agrees
    mz = arc.stm[np.ix_([2, 5], [2, 5])]
    assert abs(float(np.trace(mz))) > 2.0  # vertically unstable as well
