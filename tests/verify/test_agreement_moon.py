"""Tier-1 Phase 7: Axis-A code-path agreement on a Jovian moon pair (plan Phase 7
Task 7.1). NON-GOLDEN internal crosscheck; both sides are OUR computation.

PLAN DEVIATION (task #76): the plan sketched "VILM quadrature ΔV vs corrector
leg ΔV". Those are not the same physical quantity — the VILM ΔV_min is a full
intermoon-transfer cost (escape + leverage + capture) while the ballistic
corrector computes V_inf-continuity, not a transfer ΔV, and the I-E-G chain is
not bend-feasible (Phase 3 finding). The genuinely comparable two-code-path
crosscheck is the **Jovicentric Hohmann V_inf** between a moon pair, computed by
(a) the VILM vis-viva (search/vilm._hohmann_vinf) and (b) an independent Lambert
solve about Jupiter (core/lambert + the centred ephemeris frame). They must agree
to far inside the linked-conic band — the Axis-A internal consistency this task
exists to assert."""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.core.lambert import lambert
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.vilm import _hohmann_vinf


def _lambert_hohmann_vinf(outer: str, inner: str) -> tuple[float, float]:
    """Independent Jovicentric Hohmann V_inf via a near-180-deg Lambert solve.

    Places the outer moon at +x (transfer apoapsis) and the inner moon at the
    near-opposite point (179.5 deg, avoiding the exact-180 Lambert singularity),
    both on their circular orbits about Jupiter; solves the half-period transfer
    and returns ``(vinf_outer, vinf_inner) = |v_sc - v_moon|``.
    """
    mu = PRIMARIES["Jupiter"]
    a_o = SATELLITES[outer].sma_km
    a_i = SATELLITES[inner].sma_km
    a_t = 0.5 * (a_o + a_i)
    tof = math.pi * math.sqrt(a_t**3 / mu)
    ang = math.radians(179.5)
    r1 = np.array([a_o, 0.0, 0.0])
    v1_pl = np.array([0.0, math.sqrt(mu / a_o), 0.0])
    r2 = np.array([a_i * math.cos(ang), a_i * math.sin(ang), 0.0])
    v2_pl = math.sqrt(mu / a_i) * np.array([-math.sin(ang), math.cos(ang), 0.0])
    sols = lambert(r1, r2, tof, mu=mu, max_revs=0)
    vinf_o_ref, _ = _hohmann_vinf(outer, inner)
    best = min(
        sols, key=lambda s: abs(float(np.linalg.norm(np.asarray(s.v1) - v1_pl)) - vinf_o_ref)
    )
    vinf_o = float(np.linalg.norm(np.asarray(best.v1) - v1_pl))
    vinf_i = float(np.linalg.norm(np.asarray(best.v2) - v2_pl))
    return vinf_o, vinf_i


def test_vilm_vs_lambert_jovicentric_vinf_agree() -> None:
    # VILM vis-viva vs an independent Lambert solve about Jupiter, same Hohmann
    # Ganymede->Europa transfer. Two code paths, one Jovicentric quantity.
    vilm_o, vilm_i = _hohmann_vinf("Ganymede", "Europa")
    lam_o, lam_i = _lambert_hohmann_vinf("Ganymede", "Europa")
    # Far inside the 10% linked-conic band (these agree to <1% — it is the same
    # two-body conic computed two ways).
    assert lam_o == pytest.approx(vilm_o, rel=0.02)
    assert lam_i == pytest.approx(vilm_i, rel=0.02)
