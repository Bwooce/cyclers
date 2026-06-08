"""CONTINUOUS-from-one-seed S1L1 maintenance-dv measurement (#169, the V4 attempt).

#167 (``tests/nbody/test_s1l1_corrected_nbody.py``) confirmed the corrected S1L1
topology by reconstructing each Mars-transit leg INDEPENDENTLY from its own App-C
v_inf node — i.e. it RE-ANCHORS v_inf at every node ("Russell's per-leg
reproduction recipe"). That verifies Russell's published cycler but leaves the
SINGLE-continuous-orbit maintenance dv UNMEASURED: re-anchoring silently supplies
whatever velocity each leg needs.

This module measures what #167 left out. It propagates ONE continuous trajectory
from the FIRST Earth-departure App-C state (leg 2 — depart 2026-12-15,
v_inf = (-2.278, 5.322, 0.574) km/s, the seed of the #167 results note) forward
through every App-C node WITHOUT re-anchoring v_inf, over the real DE440 ephemeris
on the INDEPENDENT REBOUND/IAS15 integrator. At each encounter it measures the
achieved miss + v_inf and the MAINTENANCE dv — the part a ballistic flyby cannot
supply (the |v_inf| magnitude change, plus any un-bendable shortfall when the
required turn exceeds the safe-periapsis maximum bend). A real flyby rotates v_inf
for free; only that residual is a dv. See
:func:`cyclerfinder.search.s1l1_corrected.continuous_chain`.

THE GATE — three-way (honest, do NOT loosen tolerances to manufacture CLEAN):

* CLEAN  — the continuous trajectory holds all encounters in/near the SOI band
  with a BOUNDED total maintenance dv (compared to the spec §14 V3 budget
  horizon_tcm_mps = 120 m/s and the 3.0 km/s/cycle engineering bar). A small
  bounded dv means an independent integrator + ephemeris reproduces the trajectory
  AND its maintenance dv — the spec §14 V4 definition (caveat: canonical V4 names
  GMAT; ours is REBOUND/IAS15 — a GMAT cross-check would be the canonical-V4
  finisher).
* PARTIAL — holds for K cycles / one model then needs growing correction.
* DIVERGES — the continuous orbit is NOT self-consistent (per-leg re-anchoring was
  doing real work / hiding instability) -> bounds the V3 claim.

VERDICT (see ``docs/notes/2026-06-08-s1l1-continuous-v4-results.md``): **PARTIAL.**
In Russell's OWN model (Sun-only patched-conic, flybys instantaneous at the
patch points) the continuous single-seed trajectory holds all 20 nodes (E miss
<= 18,434 km, M miss <= 3,173 km — all << SOI) at the App-C v_inf (|v_inf| matched
to <= 2.7 m/s) with a BOUNDED total maintenance dv of ~62 m/s — UNDER the 120 m/s
V3 budget (STRENGTHENS V3). BUT once Mars's finite continuous gravity is modelled
(Mars-perturbed), the legs AFTER each Mars flyby (the M->E returns) diverge
(>1e8 km, integrator non-converged) because the App-C nodes are patched-conic
states that do not account for the continuous deflection through the Mars
encounter — the per-leg re-anchoring across the Mars flyby is doing real work
(QUALIFIES V3). NO catalogue writeback here.

Wall budget: the Mars-perturbed continuous chain builds a rails cache per leg
(~25 s each, 19 legs) — ~8 min, under the 25-min @slow cap.
"""

from __future__ import annotations

import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.constants import (  # noqa: E402
    AU_KM,
    MU_SUN_KM3_S2,
    PLANETS,
)
from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402
from cyclerfinder.search.s1l1_corrected import continuous_chain  # noqa: E402

# Same SOI bands as #165/#167 — kept IDENTICAL, never loosened.
_MARS_SOI_AU = PLANETS["M"].sma_au * (PLANETS["M"].mu_km3_s2 / MU_SUN_KM3_S2) ** 0.4
_EARTH_SOI_AU = PLANETS["E"].sma_au * (PLANETS["E"].mu_km3_s2 / MU_SUN_KM3_S2) ** 0.4
_MARS_BAND_AU = 3.0 * _MARS_SOI_AU  # ~0.0116 AU
_EARTH_BAND_AU = 3.0 * _EARTH_SOI_AU  # ~0.0185 AU

# Spec §14 V3 horizon-TCM budget and the engineering plausibility bar.
_V3_TCM_BUDGET_KMS = 0.120  # spec §14 horizon_tcm_mps = 120 m/s
_PLAUSIBILITY_BAR_KMS_PER_CYCLE = 3.0  # verify/plausibility convention (per cycle)
_N_CYCLES = 7


def test_continuous_chain_mechanics() -> None:
    """Fast (no @slow): the continuous chain has the right shape and is sourced.

    20 forward nodes from leg 2 (7 Mars + 13 Earth), one continuous trajectory, no
    re-anchoring. v_inf breathes the real-eph 3.2-8.0 km/s span (NOT a coplanar
    5.10 / 3.05 anchor). This is the structural check; the integrator gate is the
    @slow test below."""
    ephem = Ephemeris("astropy")
    chain = continuous_chain(ephem)  # two-body-Sun Kepler (cheap)
    assert len(chain) == 20
    assert sum(n.is_mars for n in chain) == 7
    vinfs = [n.vinf_in_kms for n in chain]
    assert min(vinfs) < 3.3, "Mars v_inf should reach the low (~3.2) end"
    assert max(vinfs) > 8.0, "Mars v_inf should reach the high (~8.0) end"
    # The maintained departures track the App-C v_inf to <~3 m/s (the chain is
    # self-consistent in Russell's Sun-only model — the dv is the residual).
    mismatch = [abs(n.vinf_in_kms - n.vinf_appc_kms) for n in chain if n.vinf_appc_kms > 0]
    assert max(mismatch) < 5e-3, f"v_inf mag mismatch {max(mismatch):.4f} km/s too large"


@pytest.mark.slow
def test_continuous_sun_only_bounded_maintenance() -> None:
    """V4 GATE (Sun-only): continuous single-seed orbit holds with BOUNDED dv.

    Propagates ONE trajectory from the first App-C Earth departure forward through
    all 20 nodes on REBOUND/IAS15, Sun-only (Russell's own patched-conic cruise
    model: flybys are instantaneous turns at the patch points), with NO v_inf
    re-anchoring. EXPECTED = the App-C per-leg v_inf; the miss / v_inf / dv are
    EVIDENCE. VERDICT pinned: CLEAN in this model — every node lands << SOI at the
    App-C v_inf and the total maintenance dv is bounded UNDER the 120 m/s V3
    budget. If a future change pushes the dv over budget or a node out of band,
    this fails deliberately: re-derive, do NOT loosen the band/budget."""
    ephem = Ephemeris("astropy")
    prop = RestrictedNBody("rebound")
    chain = continuous_chain(ephem, propagate=prop, perturbers=())
    assert len(chain) == 20

    for n in chain:
        band = _MARS_BAND_AU if n.is_mars else _EARTH_BAND_AU
        miss_au = n.miss_km / AU_KM
        assert miss_au < band, (
            f"leg {n.arrival_leg_no} {n.body}: continuous miss {miss_au:.2e} AU "
            f"outside the 3-SOI band {band:.4f} AU"
        )

    total_dv = sum(n.dv_total_kms for n in chain)
    # CLEAN: bounded under the spec §14 V3 horizon-TCM budget.
    assert total_dv < _V3_TCM_BUDGET_KMS, (
        f"total maintenance dv {total_dv * 1000:.1f} m/s exceeds the V3 budget "
        f"{_V3_TCM_BUDGET_KMS * 1000:.0f} m/s"
    )
    # ... and far under the 3.0 km/s/cycle engineering plausibility bar.
    assert total_dv / _N_CYCLES < _PLAUSIBILITY_BAR_KMS_PER_CYCLE
    # The achieved v_inf reproduces the App-C real-eph values (breathing 3.2-8.0,
    # NOT a single coplanar anchor) — the chain is self-consistent in this model.
    assert min(n.vinf_in_kms for n in chain) < 3.3
    assert max(n.vinf_in_kms for n in chain) > 8.0


@pytest.mark.slow
def test_continuous_mars_perturbed_diverges_post_flyby() -> None:
    """V4 QUALIFIER (Mars-perturbed): the continuous orbit is NOT self-consistent.

    Turns on real DE440 Mars as a continuous perturber along the whole cruise (the
    flyby gravity Russell models as instantaneous, now acting continuously) and
    runs the SAME continuous single-seed chain. VERDICT pinned: the Mars-transit
    (G) legs still arrive near Mars (in-band), but the legs AFTER each Mars flyby
    (the M->E returns) DIVERGE — the App-C nodes are patched-conic states that do
    not account for the continuous deflection through the Mars encounter, so the
    per-leg re-anchoring across the Mars flyby is doing real work. This BOUNDS the
    V3 claim: a real continuous flight needs active maintenance through each Mars
    flyby, not the silent re-anchoring. The total continuous maintenance dv here is
    O(100 km/s) — far over budget. A DIVERGE is a valid, important result; it is
    NOT massaged toward CLEAN."""
    ephem = Ephemeris("astropy")
    prop = RestrictedNBody("rebound")
    chain = continuous_chain(ephem, propagate=prop, perturbers=("M",))
    assert len(chain) == 20

    # The Mars (G-transit) arrivals themselves stay near Mars (the seed leg + the
    # re-anchored departures arrive in-band; the deepened v_inf is the real flyby
    # gravity pulling the spacecraft in).
    mars_nodes = [n for n in chain if n.is_mars]
    assert len(mars_nodes) == 7
    for n in mars_nodes:
        assert n.miss_km / AU_KM < _MARS_BAND_AU, (
            f"Mars leg {n.arrival_leg_no}: miss {n.miss_km / AU_KM:.2e} AU out of band"
        )

    # But the continuous orbit diverges: the M->E return legs (the Earth node right
    # after each Mars node) blow up by >1e7 km — the patched-conic seed is not
    # self-consistent once the Mars flyby gravity acts continuously.
    earth_after_mars = [
        chain[i] for i in range(1, len(chain)) if chain[i].body == "E" and chain[i - 1].is_mars
    ]
    assert earth_after_mars, "expected M->E return legs in the chain"
    diverged = [n for n in earth_after_mars if n.miss_km > 1.0e7]
    assert diverged, (
        "expected the post-Mars-flyby return legs to DIVERGE in the perturbed "
        "model (the re-anchoring is doing real work); none did"
    )

    # The honest total: O(100 km/s), far over both the 120 m/s V3 budget and the
    # 3.0 km/s/cycle bar — the QUALIFIER on the V3 the sibling is writing.
    total_dv = sum(n.dv_total_kms for n in chain)
    assert total_dv > _PLAUSIBILITY_BAR_KMS_PER_CYCLE * _N_CYCLES, (
        f"perturbed continuous dv {total_dv:.1f} km/s should blow past the bar"
    )
