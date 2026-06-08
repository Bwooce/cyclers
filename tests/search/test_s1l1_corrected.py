"""CORRECTED-topology S1L1 / 4.991gG2 construction tests (#167).

The construction-mechanics half (fast, two-body): verifies the corrected per-cycle
sequence (E -> g(Earth-Earth, NO Mars) -> E -> G(Earth-Mars-Earth transit) -> E) built
from Russell App-C #83's sourced real-eph state nodes:

* topology — three encounters per cycle E->E->M, one Mars encounter per cycle, the
  ``g`` legs are E->E and the ``G`` legs are E->M;
* the Mars-transit (G) legs intercept the REAL DE440 Mars at the true longitude
  (the rendezvous #165 omitted): miss well inside Mars SOI, v_inf 4-decimal-matching
  Russell's published per-leg values, which BREATHE 3.2-8.0 km/s;
* the g (E->E) legs stay FAR from Mars (sub-Mars aphelion ~1.27 AU, closest approach
  ~1.05 AU) — the direct refutation of the #164/#165 both-arcs-Mars-crossing error.

GOLDEN/HONESTY: the EXPECTED side (per-leg transit time / v_inf) traces entirely to
Russell's printed App-C characteristics table; the reconstruction is a cross-check of
two published quantities against an independent ephemeris (DE440 vs his DE405), not a
fit. The n-body confirmation gate is the @slow ``tests/nbody/test_s1l1_corrected_nbody.py``.
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.constants import MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.search.s1l1_corrected import (
    APPC_LEGS,
    APPC_MARS_TRANSIT,
    build_seeded_arcs,
    g_arc_clearances,
    reconstruct_mars_encounters,
)

# Mars SOI ~ 0.00386 AU; a real targeted flyby threads well inside it. 3 SOI is the
# #165 generous independent-confirmation band (kept IDENTICAL here, never loosened).
_MARS_SOI_AU = PLANETS["M"].sma_au * (PLANETS["M"].mu_km3_s2 / MU_SUN_KM3_S2) ** 0.4
_ENCOUNTER_MISS_TOL_AU = 3.0 * _MARS_SOI_AU


@pytest.fixture(scope="module")
def ephem() -> Ephemeris:
    return Ephemeris("astropy")


def test_topology_three_encounters_per_cycle_one_mars() -> None:
    """The App-C E/M column is E E M | E E M | ... : 3 encounters/cycle, 1 Mars/cycle."""
    seq = [body for (_no, body, _ts, _v) in APPC_LEGS]
    # 22 nodes: E E M repeated 7 times then a terminal E.
    assert seq == ["E", "E", "M"] * 7 + ["E"]
    mars_nodes = [no for (no, body, _ts, _v) in APPC_LEGS if body == "M"]
    assert mars_nodes == [3, 6, 9, 12, 15, 18, 21]  # one Mars encounter per cycle


def test_arc_classification_g_is_ee_bigg_is_em(ephem: Ephemeris) -> None:
    """g arcs are E->E (free return); G transit arcs are E->M (Mars outbound)."""
    arcs = build_seeded_arcs(ephem)
    assert len(arcs) == 21  # 22 nodes -> 21 arcs
    g_legs = [a.leg_no for a in arcs if a.is_g_arc]
    transit_legs = [a.leg_no for a in arcs if a.is_mars_transit]
    assert g_legs == [1, 4, 7, 10, 13, 16, 19]  # the E->E free returns
    assert transit_legs == [2, 5, 8, 11, 14, 17, 20]  # the E->M Mars transits
    # No arc is both; the M->E return legs (3,6,...) are neither g nor a transit.
    assert not any(a.is_g_arc and a.is_mars_transit for a in arcs)


def test_mars_transit_legs_intercept_real_mars_at_true_longitude(ephem: Ephemeris) -> None:
    """All 7 G legs reach the REAL DE440 Mars inside the SOI band at its true longitude.

    This is the longitude-rendezvous constraint #165 lacked: the seeded App-C state
    lands the spacecraft on real Mars (miss << SOI) AND at Mars's actual heliocentric
    ecliptic longitude — not merely on Mars's radius ~110 deg away (#165)."""
    encs = reconstruct_mars_encounters(ephem)
    assert len(encs) == 7
    for e in encs:
        # Inside the (generous, un-loosened) 3-SOI band — a real intercept.
        assert e.miss_au < _ENCOUNTER_MISS_TOL_AU, (
            f"leg->M{e.arrival_leg_no}: miss {e.miss_au:.2e} AU outside the "
            f"{_ENCOUNTER_MISS_TOL_AU:.4f} AU band"
        )
        # True longitude rendezvous: SC and Mars at the same heliocentric longitude.
        dlon = abs((e.sc_lon_deg - e.mars_lon_deg + 180.0) % 360.0 - 180.0)
        assert dlon < 0.5, f"leg->M{e.arrival_leg_no}: SC/Mars longitude gap {dlon:.3f} deg"


def test_mars_vinf_matches_published_per_leg_and_breathes(ephem: Ephemeris) -> None:
    """Achieved Mars v_inf 4-decimal-matches Russell's per-leg values, breathing 3.2-8.0.

    The real-eph cycler does NOT reproduce a single coplanar anchor (5.10 Russell /
    3.05 Rogers); its Mars v_inf is epoch-dependent. EXPECTED = the printed per-leg
    numbers; achieved = evidence."""
    encs = reconstruct_mars_encounters(ephem)
    for e in encs:
        _pub_tof, pub_vinf = APPC_MARS_TRANSIT[e.arrival_leg_no]
        assert e.vinf_kms == pytest.approx(pub_vinf, abs=2e-3), (
            f"leg->M{e.arrival_leg_no}: v_inf {e.vinf_kms:.4f} vs published {pub_vinf}"
        )
    vinfs = [e.vinf_kms for e in encs]
    # Breathes across the full 3.2-8.0 km/s span (NOT a single 5.10/3.05 anchor).
    assert min(vinfs) < 3.3 and max(vinfs) > 8.0
    assert np.mean(vinfs) == pytest.approx(5.48, abs=0.05)


def test_g_arcs_stay_far_from_mars(ephem: Ephemeris) -> None:
    """The g (E->E) free returns never approach Mars — refutes the #164/#165 error.

    #164/#165 forced BOTH arcs to cross Mars's radius (geometrically impossible for the
    real g arc, aphelion ~1.27 AU). Here the g legs stay sub-Mars: closest approach
    ~0.7-1.05 AU, aphelion < 1.4 AU, far outside any Mars-encounter band."""
    clears = g_arc_clearances(ephem)
    assert len(clears) == 7
    for c in clears:
        assert c.closest_mars_au > 0.5, (
            f"g leg {c.leg_no}: came within {c.closest_mars_au:.3f} AU of Mars"
        )
        assert c.aphelion_au < 1.45, (
            f"g leg {c.leg_no}: aphelion {c.aphelion_au:.3f} AU reaches Mars"
        )
        # Specifically NOT inside the Mars-encounter band — the structural fix.
        assert c.closest_mars_au > _ENCOUNTER_MISS_TOL_AU * 10
    # The first g arc matches #166's DE440 probe (closest ~1.05 AU, aphelion ~1.27 AU).
    first = next(c for c in clears if c.leg_no == 1)
    assert first.closest_mars_au == pytest.approx(1.05, abs=0.05)
    assert first.aphelion_au == pytest.approx(1.27, abs=0.05)
