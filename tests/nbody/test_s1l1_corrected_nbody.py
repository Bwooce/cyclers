"""INDEPENDENT n-body confirmation of the CORRECTED-topology S1L1 closure (#167).

The decisive gate — the corrected #165. #165 propagated the #164 two-arc construction
(BOTH arcs modelled as Mars-crossing, no longitude target) and found a ~2.6 AU DRIFT.
#166 found the topology error (lowercase g is a pure Earth-to-Earth free return; only
uppercase G is the Mars-transit leg) AND the exact real-eph seed (Russell App-C #83).
This gate propagates the corrected construction — seeded entirely from Russell's
printed per-leg state — through the INDEPENDENT REBOUND/IAS15 integrator over the real
DE440 ephemeris, sharing NONE of any solver's machinery.

What independence means
-----------------------
The ONLY thing seeded is the per-leg departure state (Russell App-C v_inf + the real
DE440 planet velocity at the printed epoch). Everything downstream — where the
spacecraft actually goes, where DE440 Mars actually is, the achieved miss / v_inf — is
pure integrator / ephemeris output, never fit.

Two physically distinct independent checks, both reported (the honest characterization
of a patched-conic seed):

* **Sun-only IAS15** — matches Russell's own patched-conic cruise model (each ballistic
  leg is two-body-Sun; planet gravity enters only as the instantaneous flyby bend at a
  node, confirmed by the App-C dv ~1e-11). This is the apples-to-apples reproduction of
  the published seed by an INDEPENDENT integrator (GOLDEN GATE 1 pins IAS15 to the
  analytic two-body to < 1 km / 120 d). RESULT: all 7 Mars encounters land 380-3200 km
  from real DE440 Mars (<< 1 SOI) at the published per-leg v_inf to 4 decimals.

* **Mars-perturbed IAS15** — turns on real DE440 Mars as a continuous perturber along
  the whole cruise (the real flyby gravity Russell models as instantaneous, now acting
  continuously). RESULT: all 7 encounters land 6,600-40,900 km from real Mars (< 0.1
  Mars SOI), still comfortably inside the 3-SOI confirmation band.

(The earlier #165-style E+M run that started the spacecraft a full Earth-SOI radius
*displaced* along the departure asymptote produced larger 0.01-0.03 AU misses — that is
the patched-conic departure artifact, not real geometry: the App-C seed is a departure
VELOCITY off Earth, not a state already clear of Earth's well, so a continuous Earth
perturbation in the departure region bends it; the miss shrinks monotonically toward
the SOI band as the start point moves clear of Earth. It is documented in the results
note, not gated here, because it conflates the seed-handoff convention with the
geometry.)

THE GATE — three-way
--------------------
* CONFIRMED: Mars encounters land in-band at the App-C real-eph v_inf (breathing
  3.2-8.0 km/s, NOT the coplanar 5.10) AND the g arcs stay far from Mars. -> achieved.
* PARTIAL: closer than #165's 2.6 AU but not all in-band.
* DRIFT-AGAIN: still misses.

VERDICT: **CONFIRMED.** All 7 Mars encounters land inside the 3-SOI band (the same
band #165 used, NOT loosened) at the published per-leg v_inf, on an INDEPENDENT
integrator, and the g arcs stay sub-Mars. S1L1 / ``russell-ch4-4.991gG2`` is
independently confirmed on DE440 from a sourced real-eph seed. No catalogue writeback
here — the V3 promotion is the main session's call. See
``docs/notes/2026-06-08-s1l1-corrected-closure-results.md``.

Wall budget: ~2 s for all 14 propagations — far under the 25-min @slow cap.
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
from cyclerfinder.nbody.forces import RailsEphemerisCache  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402
from cyclerfinder.search.s1l1_corrected import (  # noqa: E402
    APPC_MARS_TRANSIT,
    build_seeded_arcs,
    g_arc_clearances,
)

# Same 3-Mars-SOI confirmation band as #165 — kept IDENTICAL, never loosened.
_MARS_SOI_AU = PLANETS["M"].sma_au * (PLANETS["M"].mu_km3_s2 / MU_SUN_KM3_S2) ** 0.4
_ENCOUNTER_MISS_TOL_AU = 3.0 * _MARS_SOI_AU  # ~0.0116 AU


def _mars_encounter(
    arc: object,
    ephem: Ephemeris,
    prop: RestrictedNBody,
    *,
    perturbers: tuple[str, ...],
) -> tuple[float, float, bool, float]:
    """Independent IAS15 propagation of one Mars-transit leg to its arrival epoch.

    Returns ``(miss_au, vinf_kms, converged, energy_rel_drift)`` against real DE440
    Mars. ``perturbers`` empty => Sun-only (Russell's patched-conic cruise model);
    ``("M",)`` => real Mars as a continuous perturber.
    """
    t0 = arc.t0_sec  # type: ignore[attr-defined]
    t1 = arc.t1_sec  # type: ignore[attr-defined]
    cache = None
    if perturbers:
        cache = RailsEphemerisCache(
            perturbers, ephem, t0 - 5 * SECONDS_PER_DAY, t1 + 10 * SECONDS_PER_DAY
        )
    out = prop.propagate(
        arc.r0_km,  # type: ignore[attr-defined]
        arc.v0_km_s,  # type: ignore[attr-defined]
        t0,
        t1,
        bodies=perturbers,
        ephem=ephem if perturbers else None,
        cache=cache,
        accuracy=1e-11 if not perturbers else 1e-10,
        max_wall_sec=120.0,
    )
    r_m, v_m = ephem.state("M", t1)
    miss_au = float(np.linalg.norm(out.r_km - np.asarray(r_m)) / AU_KM)
    vinf_kms = float(np.linalg.norm(out.v_km_s - np.asarray(v_m)))
    return miss_au, vinf_kms, bool(out.converged), float(out.energy_rel_drift)


@pytest.mark.slow
def test_s1l1_corrected_nbody_confirmed() -> None:
    """INDEPENDENT n-body confirmation of the corrected S1L1 closure — the V3 gate.

    Propagates all 7 Mars-transit (G) legs of the corrected construction, each seeded
    from Russell App-C #83's printed per-leg state, through REBOUND/IAS15 over the real
    DE440 ephemeris. EXPECTED = the published per-leg v_inf; the miss / v_inf are
    EVIDENCE. VERDICT pinned: CONFIRMED — all 7 encounters inside the 3-SOI band at the
    published v_inf, on BOTH the Sun-only (patched-conic) and Mars-perturbed models,
    and the g arcs stay far from Mars. This INVERTS the #165 DRIFT for the corrected
    topology. If a future change pushes any encounter out of band, this assert fails
    deliberately — re-derive, do not loosen the band."""
    ephem = Ephemeris("astropy")
    prop = RestrictedNBody("rebound")
    arcs = [a for a in build_seeded_arcs(ephem) if a.is_mars_transit]
    assert len(arcs) == 7

    sun_vinfs: list[float] = []
    for arc in arcs:
        arrival_no = arc.leg_no + 1
        _pub_tof, pub_vinf = APPC_MARS_TRANSIT[arrival_no]

        # (1) Sun-only independent integrator == Russell's patched-conic cruise model.
        miss_s, vinf_s, conv_s, edrift_s = _mars_encounter(arc, ephem, prop, perturbers=())
        assert conv_s, f"M{arrival_no}: Sun-only leg did not converge"
        assert abs(edrift_s) < 1e-6, f"M{arrival_no}: energy drift {edrift_s:.2e}"
        assert miss_s < _ENCOUNTER_MISS_TOL_AU, (
            f"M{arrival_no}: Sun-only miss {miss_s:.2e} AU outside the "
            f"{_ENCOUNTER_MISS_TOL_AU:.4f} AU band — NOT confirmed"
        )
        assert miss_s < 5e-5, f"M{arrival_no}: Sun-only miss {miss_s:.2e} AU larger than expected"
        # v_inf matches the published per-leg value to 4 decimals (real-eph, breathes).
        assert vinf_s == pytest.approx(pub_vinf, abs=2e-3), (
            f"M{arrival_no}: Sun-only v_inf {vinf_s:.4f} vs published {pub_vinf}"
        )
        sun_vinfs.append(vinf_s)

        # (2) Mars-perturbed: the real flyby gravity acting continuously. Still in-band.
        miss_m, _vinf_m, conv_m, _edm = _mars_encounter(arc, ephem, prop, perturbers=("M",))
        assert conv_m, f"M{arrival_no}: Mars-perturbed leg did not converge"
        assert miss_m < _ENCOUNTER_MISS_TOL_AU, (
            f"M{arrival_no}: Mars-perturbed miss {miss_m:.2e} AU outside the band"
        )

    # The Mars v_inf BREATHES across the full 3.2-8.0 km/s span (the real-eph
    # signature), NOT a single coplanar 5.10 / Rogers 3.05 anchor.
    assert min(sun_vinfs) < 3.3, "Mars v_inf should reach the low (~3.2) end"
    assert max(sun_vinfs) > 8.0, "Mars v_inf should reach the high (~8.0) end"
    assert np.mean(sun_vinfs) == pytest.approx(5.48, abs=0.05)


@pytest.mark.slow
def test_s1l1_corrected_g_arcs_stay_far_from_mars() -> None:
    """The corrected g (E->E) free returns never approach Mars (the structural fix).

    #164/#165 forced both arcs to cross Mars's radius; the real g arc has aphelion
    ~1.27 AU and cannot reach Mars (1.52 AU). Confirms every g leg's closest approach
    to real DE440 Mars is far outside any encounter band — the topology #165 got
    wrong."""
    ephem = Ephemeris("astropy")
    clears = g_arc_clearances(ephem)
    assert len(clears) == 7
    for c in clears:
        assert c.closest_mars_au > _ENCOUNTER_MISS_TOL_AU * 20, (
            f"g leg {c.leg_no}: closest Mars {c.closest_mars_au:.3f} AU — too near"
        )
        assert c.aphelion_au < 1.45, (
            f"g leg {c.leg_no}: aphelion {c.aphelion_au:.3f} AU reaches Mars"
        )
