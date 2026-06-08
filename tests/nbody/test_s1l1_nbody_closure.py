"""INDEPENDENT n-body cross-validation of the #164 S1L1 two-body closure (#165).

The decisive gate for whether S1L1 / ``russell-ch4-4.991gG2`` (#94) earns the
catalogue's first V3. #164 found a CLOSED two-arc geometry with the TWO-BODY
radial-crossing continuation SOLVER (``search/continuation_chain.py``): both arcs'
emerged V_inf land at the SOURCED anchors (E 4.99 / M 5.10) and both descriptor
ToFs are reached. This test confirms — or refutes — that geometry with an
**INDEPENDENT integrator**: REBOUND / IAS15 over the real DE440 planet ephemeris
(``nbody/propagator.py``), which shares NONE of the continuation solver's machinery.

What independence means here
----------------------------
The ONLY thing seeded from the continuation fit is the spacecraft *initial state*
(the #164 closed elements -> a heliocentric departure state at the real DE440
Earth, frame/time-converted via ``nbody/convert.py`` conventions). Everything
downstream — where the spacecraft actually goes, where the real DE440 Mars
actually is, the achieved Mars V_inf and miss-distance — is pure n-body / DE440
output, never fit. The sourced anchors (E 4.99 / M 5.10) are EXPECTED; the n-body
miss-distance and V_inf are EVIDENCE.

THE GATE
--------
* CONFIRMED: the independent propagation reaches the real DE440 Mars (miss within
  a few Mars SOI) with V_inf within tolerance of the 5.10 anchor at the closed
  encounter epoch -> S1L1 is independently confirmed.
* DRIFT: the independent trajectory departs the continuation geometry (real Mars
  is not where the construction puts the crossing) -> S1L1 is "two-body closeable,
  n-body-unconfirmed" — a continuation artifact, not a real n-body cycler. A VALID,
  publishable outcome, reported plainly; the tolerance is NOT loosened to force a
  pass.

THE FINDING (see docs/notes/2026-06-08-s1l1-nbody-crossval-results.md): **DRIFT.**
The #164 construction is a radial-crossing + V_inf + ToF closure that does NOT
enforce *longitude rendezvous* with the real ephemeris Mars. The spacecraft
crosses Mars's orbital RADIUS at the right time with the right speed but ~110 deg
away in longitude from where DE440 Mars actually is; the achieved Mars miss is
~2.6 AU (both arcs), confirmed identically by an independent pure two-body Kepler
propagation of the same absolute state (so it is structural, not an integrator
artifact). NO catalogue writeback; S1L1 promotion is the main session's call.

Wall budget: the @slow encounter test is bounded well under the 25 min cap (the
n-body legs are sub-second once the spacecraft starts outside Earth's SOI; the
dominant cost is the one-off continuation re-solve, ~2 s).
"""

from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.constants import (  # noqa: E402
    AU_KM,
    MU_SUN_KM3_S2,
    PLANETS,
    SECONDS_PER_DAY,
)
from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.core.kepler import coe_to_rv  # noqa: E402
from cyclerfinder.core.kepler import propagate as kepler_propagate  # noqa: E402
from cyclerfinder.nbody.forces import RailsEphemerisCache  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402
from cyclerfinder.search.continuation_chain import (  # noqa: E402
    _ChainSolve,
    _earth_vinf_vector_ramped,
    _effective_planet,
    _planet_phase,
    _RampedArc,
    continuation_chain_correct,
)
from cyclerfinder.search.free_return import _crossing  # noqa: E402
from cyclerfinder.search.free_return_chain import free_return_chain_correct  # noqa: E402

# --- S1L1 / 4.991gG2 SOURCED anchors (Russell 2004 Table 4.9, the row's OWN
# anchors). These are the EXPECTED side of the gate. ---------------------------
_4991_APHELION = 1.64
_4991_G_TOF = 1.4612
_4991_BIGG_TOF = 2.8096
_4991_VINF_E = 4.99  # km/s, sourced Earth anchor (EXPECTED)
_4991_VINF_M = 5.10  # km/s, sourced Mars anchor (EXPECTED)

# --- The independent-confirmation tolerance, grounded in the harness's golden
# gates and basic flyby physics (NOT loosened to force a pass). ---------------
#
# A genuine Mars flyby requires the spacecraft inside Mars's gravitational sphere
# of influence: r_SOI = a_Mars * (m_Mars / m_Sun)^(2/5) ~ 0.0038 AU ~ 5.8e5 km.
# We accept "an encounter occurred" only if the n-body miss-distance is within a
# few Mars SOI (a real targeted flyby threads the B-plane well inside one SOI;
# 3 SOI is a generous independent-confirmation band). The integrator itself is
# trustworthy to FAR below this scale — GOLDEN GATE 1 pins REBOUND/IAS15 to the
# analytic two-body solution to < 1 km over 120 d, GOLDEN GATE 3 pins the DE440
# planet read to < 1e-6 km, and GOLDEN GATE 2 pins Sun-only energy drift to
# < 1e-10 — so a multi-SOI miss is a real geometric gap, never integrator noise.
_MARS_SOI_AU = PLANETS["M"].sma_au * (PLANETS["M"].mu_km3_s2 / MU_SUN_KM3_S2) ** 0.4
_ENCOUNTER_MISS_TOL_AU = 3.0 * _MARS_SOI_AU  # ~0.0114 AU


def _closed_geometry() -> _ChainSolve:
    """Re-solve the #164 two-arc continuation to recover the CLOSED geometry.

    Returns the ``best_final`` ``_ChainSolve`` (closed arcs + encounter epoch).
    This is the ONLY contact with the continuation solver — it supplies the
    initial state; nothing downstream is fit.
    """
    seed = free_return_chain_correct(
        _4991_APHELION, _4991_G_TOF, _4991_BIGG_TOF, _4991_VINF_E, _4991_VINF_M
    )
    a1, e1 = seed.arcs[0].a_au, seed.arcs[0].e
    a2, e2 = seed.arcs[1].a_au, seed.arcs[1].e
    res = continuation_chain_correct(
        a1, e1, a2, e2, 0.0, _4991_G_TOF, _4991_BIGG_TOF, _4991_VINF_E, _4991_VINF_M
    )
    assert res.best_final is not None
    return res.best_final


def _departure_state(
    arc: _RampedArc, t_dep_sec: float, ephem: Ephemeris
) -> tuple[np.ndarray, np.ndarray]:
    """Build the absolute heliocentric departure state for one closed arc.

    The honest patched-conic departure: place the spacecraft at the **real DE440
    Earth position** at ``t_dep_sec`` with velocity ``v_Earth_DE440 + V_inf``,
    where the V_inf vector's (radial, tangential, normal) components come from the
    #164 closed geometry and are oriented in the REAL Earth's local encounter
    frame. The frame/time axis is the harness's (TDB-J2000-seconds, heliocentric
    J2000-ecliptic) — the same ``Ephemeris('astropy').state`` axis the propagator
    and ``nbody/convert.py`` use, so no extra rotation is applied (GOLDEN GATE 3).
    """
    a_km = arc.a_au * AU_KM
    # The continuation geometry's V_inf vector lives in a perifocal-aligned frame;
    # decompose it onto the local (r_hat, t_hat, n_hat) basis at the Earth crossing.
    ep_e = _effective_planet(
        "E", _planet_phase("E", t_dep_sec, MU_SUN_KM3_S2), 1.0, 1.0, mu=MU_SUN_KM3_S2
    )
    nu_e = _crossing(a_km, arc.e, ep_e.r_eff_km)
    r_sc, _v_sc = coe_to_rv(a_km, arc.e, nu_e, MU_SUN_KM3_S2)
    r_hat = np.asarray(r_sc) / np.linalg.norm(r_sc)
    t_hat = np.array([-r_hat[1], r_hat[0], 0.0])
    vinf_vec = _earth_vinf_vector_ramped(arc, 1.0, 1.0, mu=MU_SUN_KM3_S2)
    v_r = float(np.dot(vinf_vec, r_hat))
    v_t = float(np.dot(vinf_vec, t_hat))
    v_n = float(vinf_vec[2])
    # Real DE440 Earth state -> its local encounter frame.
    r_earth, v_earth = ephem.state("E", t_dep_sec)
    rhat_e = r_earth / np.linalg.norm(r_earth)
    that_e = v_earth - np.dot(v_earth, rhat_e) * rhat_e
    that_e /= np.linalg.norm(that_e)
    nhat_e = np.cross(rhat_e, that_e)
    vinf_abs = v_r * rhat_e + v_t * that_e + v_n * nhat_e
    return np.asarray(r_earth, dtype=np.float64), np.asarray(v_earth + vinf_abs, dtype=np.float64)


def _propagate_to_mars_crossing(
    arc: _RampedArc, ephem: Ephemeris
) -> tuple[float, float, float, bool, float]:
    """Independent n-body propagation of one arc to its Mars-crossing epoch.

    Returns ``(miss_au, vinf_kms, edrift, converged, miss_2body_au)``. The
    spacecraft is started just outside Earth's SOI along the outgoing V_inf
    asymptote (the patched-conic departure point) so IAS15 does not grind on the
    Earth point-mass singularity at the departure node. The reference miss is also
    computed with a pure two-body Kepler propagation of the SAME absolute state —
    if the two agree, the gap is geometric (the construction), not integrator-side.
    """
    bf_t0 = arc.t0_sec
    r0, v0 = _departure_state(arc, bf_t0, ephem)
    tof_em_sec = arc.geometry.tof_em_days * SECONDS_PER_DAY
    t1 = bf_t0 + tof_em_sec
    r_mars, v_mars = ephem.state("M", t1)

    # Independent integrator (REBOUND/IAS15, planets on DE440 rails).
    vinf_dir = v0 - ephem.state("E", bf_t0)[1]
    vinf_dir = vinf_dir / np.linalg.norm(vinf_dir)
    soi_e_km = PLANETS["E"].sma_au * (PLANETS["E"].mu_km3_s2 / MU_SUN_KM3_S2) ** 0.4 * AU_KM
    r0_outside = r0 + vinf_dir * soi_e_km
    cache = RailsEphemerisCache(
        ("E", "M"), ephem, bf_t0 - 5 * SECONDS_PER_DAY, t1 + 10 * SECONDS_PER_DAY
    )
    arc_out = RestrictedNBody("rebound").propagate(
        r0_outside,
        v0,
        bf_t0,
        t1,
        bodies=("E", "M"),
        ephem=ephem,
        cache=cache,
        accuracy=1e-10,
        max_wall_sec=120.0,
    )
    miss_au = float(np.linalg.norm(arc_out.r_km - r_mars) / AU_KM)
    vinf_kms = float(np.linalg.norm(arc_out.v_km_s - v_mars))

    # Independent two-body cross-check of the SAME absolute state.
    rk, _vk = kepler_propagate(r0, v0, tof_em_sec)
    miss_2body_au = float(np.linalg.norm(rk - r_mars) / AU_KM)

    return (
        miss_au,
        vinf_kms,
        float(arc_out.energy_rel_drift),
        bool(arc_out.converged),
        miss_2body_au,
    )


def test_mars_soi_tolerance_is_grounded_not_inflated() -> None:
    """The encounter tolerance is a few Mars SOI — a real geometric band, NOT a
    knob loosened to manufacture a pass. Pins the tolerance derivation."""
    assert abs(_MARS_SOI_AU - 0.00384) < 1e-4
    # 3 SOI is ~0.0115 AU ~ 1.7e6 km: generous for confirmation, still ~230x
    # tighter than the observed ~2.6 AU drift, so it cannot mask the verdict.
    assert _ENCOUNTER_MISS_TOL_AU < 0.02


@pytest.mark.slow
def test_s1l1_nbody_crossvalidation_verdict() -> None:
    """INDEPENDENT n-body cross-check of the #164 S1L1 closure — the V3 gate.

    Builds each closed arc's departure state at the real DE440 Earth, propagates
    it forward in REBOUND/IAS15 (planets on DE440 rails) to the closed Mars
    crossing epoch, and measures the achieved miss-distance and V_inf to the real
    DE440 Mars. The sourced anchor (M 5.10) is EXPECTED; the miss / V_inf are
    EVIDENCE.

    OUTCOME: DRIFT. Both arcs miss real Mars by ~2.6 AU — the construction closes
    on Mars's RADIUS + V_inf + ToF but real Mars is ~110 deg away in longitude at
    the crossing. The independent two-body propagation of the same state shows the
    identical ~2.6 AU miss, so it is structural, not integrator noise. This test
    PINS the DRIFT verdict: S1L1 is two-body closeable but n-body-UNCONFIRMED. It
    will flip to a CONFIRMED assertion only if a future construction (e.g. one that
    adds a longitude-rendezvous constraint) brings the miss inside the encounter
    tolerance — at which point the asserts below must be inverted deliberately."""
    ephem = Ephemeris("astropy")
    bf = _closed_geometry()

    results: dict[str, tuple[float, float, float, bool, float]] = {}
    for tag, arc in (("arc1_g", bf.arc1), ("arc2_G", bf.arc2)):
        results[tag] = _propagate_to_mars_crossing(arc, ephem)

    # The integrator is healthy on both legs (so the miss is a real geometric gap,
    # not an integration breakdown): converged, Sun-dominated energy drift tiny.
    for tag, (_m, _v, edrift, converged, _m2) in results.items():
        assert converged, f"{tag}: n-body leg did not converge — integrator issue, not geometry"
        assert abs(edrift) < 1e-6, f"{tag}: energy drift {edrift:.2e} above the §5.2 floor"

    # DRIFT verdict, pinned: BOTH arcs miss real DE440 Mars by FAR more than the
    # (generous) encounter tolerance — S1L1 is NOT independently confirmed.
    for tag, (miss_au, _vinf, _e, _c, miss_2body) in results.items():
        assert miss_au > _ENCOUNTER_MISS_TOL_AU, (
            f"{tag}: n-body Mars miss {miss_au:.3f} AU is INSIDE the "
            f"{_ENCOUNTER_MISS_TOL_AU:.4f} AU encounter band — this would be "
            "CONFIRMED; invert this test deliberately if so."
        )
        # Structural, not integrator-side: pure two-body propagation of the same
        # absolute state reproduces the same multi-AU miss (within 0.1 AU).
        assert miss_2body == pytest.approx(miss_au, abs=0.1), (
            f"{tag}: n-body miss {miss_au:.3f} AU and two-body miss {miss_2body:.3f} AU disagree — "
            "an integrator artifact would; they agree, so the gap is geometric."
        )

    # The observed drift is severe (~2.6 AU on both arcs) — record the magnitude.
    assert results["arc1_g"][0] > 1.0
    assert results["arc2_G"][0] > 1.0
