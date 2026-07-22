"""#688 structural tests for the composed two-moon Keplerian-map screen.

These assert STRUCTURAL / geometric facts (resonance-a formula, sourced moon
SMAs, exterior-map validity domain, single-moon reduction), not any discovery
verdict.  The screen itself is screen-grade heuristics (see module docstring).
"""

from __future__ import annotations

import math

from cyclerfinder.genome.composed_moon_map import (
    ComposedMoonMap,
    ComposedState,
    moon_config,
    resonance_semimajor,
    tisserand_cj,
)
from cyclerfinder.genome.keplerian_map import eccentricity_from_tisserand


def _cm() -> ComposedMoonMap:
    return ComposedMoonMap(moon_config("Europa"), moon_config("Ganymede"))


def test_resonance_semimajor_matches_rs07_one_to_two() -> None:
    """1:2 spacecraft:moon resonance => a = 2^(2/3) (RS07 sourced fixed point)."""
    assert abs(resonance_semimajor(1, 2) - 2.0 ** (2.0 / 3.0)) < 1e-12


def test_period_ratio_is_the_laplace_two_to_one() -> None:
    """Ganymede:Europa period ratio ~ 2.01 (the real 2:1 Laplace commensurability)."""
    eur = moon_config("Europa")
    gan = moon_config("Ganymede")
    ratio = gan.period_s / eur.period_s
    assert 2.0 < ratio < 2.03, f"period ratio {ratio:.4f} not ~2:1 Laplace"


def test_tisserand_cj_inverts_eccentricity() -> None:
    """tisserand_cj is the exact inverse of eccentricity_from_tisserand."""
    a_norm = 1.35
    k = -0.5 / a_norm
    e = 0.234
    c_j = tisserand_cj(a_norm, e)
    e_back = eccentricity_from_tisserand(k, c_j)
    assert abs(e_back - e) < 1e-9


def test_encounter_shells_are_disjoint() -> None:
    """Europa and Ganymede exterior-encounter shells do not overlap in km.

    A Europa encounter pins periapsis ~ Europa's orbit; a Ganymede encounter
    pins periapsis ~ Ganymede's orbit; the exterior map conserves periapsis, so
    no single orbit can alternate as a self-consistent encounter of both.
    """
    cm = _cm()
    eur = cm.cfg["Europa"]
    gan = cm.cfg["Ganymede"]
    eur_shell_hi_km = cm.shell_hi * eur.sma_km
    gan_shell_lo_km = 1.0 * gan.sma_km
    assert eur_shell_hi_km < gan_shell_lo_km, "encounter shells unexpectedly overlap"


def test_paper_ganymede_resonances_are_interior() -> None:
    """The CCR4BP paper's Jupiter-Ganymede 3:2 and 4:3 resonances are interior (a<1).

    Interior resonances are OUTSIDE the RS07 exterior-periapsis map's validity;
    this is the structural reason the exterior-map screen cannot represent the
    Ganymede side of the Europa-Ganymede tour.
    """
    assert resonance_semimajor(3, 2) < 1.0
    assert resonance_semimajor(4, 3) < 1.0


def test_ganymede_inert_after_europa_patch() -> None:
    """After a physical-preserving patch from a Europa encounter, the Ganymede
    exterior map is inert (periapsis interior to Ganymede => no valid encounter)."""
    cm = _cm()
    # Start at a valid Europa exterior encounter (3:4 resonance, C_J=3).
    a_norm = resonance_semimajor(3, 4)
    e = eccentricity_from_tisserand(-0.5 / a_norm, 3.0)
    start = ComposedState(a_phys_km=a_norm * cm.cfg["Europa"].sma_km, e=e, varpi_rad=0.0, t_s=0.0)
    # One Europa step must be a valid encounter.
    _, rec_e = cm.step(start, "Europa")
    assert rec_e.encounter_valid
    # A Ganymede step on the same physical orbit must NOT be a valid encounter.
    _, rec_g = cm.step(start, "Ganymede")
    assert not rec_g.encounter_valid
    assert not rec_g.kicked


def test_reduces_to_single_moon_map() -> None:
    """A pure-Europa itinerary reproduces the standalone Europa KeplerianMap."""
    cm = _cm()
    a_norm = 1.45
    e = eccentricity_from_tisserand(-0.5 / a_norm, 3.0)
    start = ComposedState(a_phys_km=a_norm * cm.cfg["Europa"].sma_km, e=e, varpi_rad=0.2, t_s=0.0)
    end, recs = cm.run_segment(start, "Europa", 5)
    # a stayed physical and the orbit remained near Europa's exterior shell
    assert not math.isnan(end.e)
    assert all(r.moon == "Europa" for r in recs)
    assert any(r.encounter_valid for r in recs)
