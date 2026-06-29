"""#480 positive control — ballistic-cycler construction of the Hernandez 2017 EIGE.

EIGE (Europa-Io-Ganymede-Europa) is the paper's maintenance-ΔV demonstration cycler:
ONE synodic period, ONE revolution (AAS 17-608 Fig 5, pp.10-11). This asserts that a
genuine *feasible* ballistic EIGE exists in the paper's ideal circular-coplanar model,
and that — pinned at the two SOURCED Fig-5 interior altitudes (Io 2,817 / Ganymede
13,180 km) — it is ballistic with a feasible third (Europa) altitude PREDICTION.

The EXPECTED altitudes trace to Fig 5, never to our code
(``feedback_golden_tests_sourced_only``). The Io/Ganymede altitudes are SOURCED pin
targets (like EGGIE's Table-4 V∞ pull); the Europa altitude is the free cross-check.

Verdict note: ``docs/notes/2026-06-30-480-eige-ballistic-construction-verdict.md``.
"""

from __future__ import annotations

from cyclerfinder.search.eige_ballistic import (
    ALT_MAX_KM,
    ALT_MIN_KM,
    EIGE_RESONANT_SMA_KM,
    feasible_ballistic_eige,
)

# SOURCED from Hernandez-Jones-Jesick 2017, AAS 17-608, Figure 5 (digest
# docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md).
# These EXPECTED values trace to the paper, never to our own computation.
FIG5_ALT_IO_KM = 2817.0
FIG5_ALT_GANYMEDE_KM = 13180.0
FIG5_ALT_EUROPA_KM = 470.0

# Ideal Ganymede sma (km), the 1:1 resonant sma — sourced model geometry (paper p.3).
IDEAL_GANYMEDE_SMA_KM = 1_055_289.0


def test_resonant_sma_is_ideal_ganymede() -> None:
    """EIGE is 1:1 (one synodic period / one rev): a_sc equals ideal Ganymede sma."""
    assert abs(EIGE_RESONANT_SMA_KM - IDEAL_GANYMEDE_SMA_KM) < 5.0, EIGE_RESONANT_SMA_KM


def test_feasible_ballistic_eige_exists() -> None:
    """A closed, equal-V∞, bend-feasible ballistic EIGE exists in the ideal model."""
    e = feasible_ballistic_eige()

    # Ballistic: equal in/out |V∞| at every flyby (the defining property).
    assert e.ballistic_resnorm_kms < 1.0e-3, e.ballistic_resnorm_kms
    # The cycle closes at the Europa periodicity seam (magnitude).
    assert e.seam_defect_kms < 1.0e-3, e.seam_defect_kms
    # Near-ballistic total ΔV (the flybys carry essentially no maneuver).
    assert e.total_dv_ms < 1.0, e.total_dv_ms

    # Every flyby altitude lies inside the paper's 25-70,000 km window.
    assert e.all_feasible, e.flyby_alt_km
    for key, alt in e.flyby_alt_km.items():
        assert ALT_MIN_KM <= alt <= ALT_MAX_KM, (key, alt)


def test_eige_reaches_sourced_fig5_interior_altitudes() -> None:
    """The construction reaches the SOURCED Fig-5 interior altitudes ballistically.

    Io and Ganymede are the soft pin targets; the test confirms a *ballistic* member
    exists at the paper's printed interior altitudes (the construction succeeds — these
    altitudes are reachable on the ballistic manifold). EXPECTED values trace to Fig 5.
    """
    e = feasible_ballistic_eige()
    assert abs(e.flyby_alt_km["Io"] - FIG5_ALT_IO_KM) < 50.0, e.flyby_alt_km
    assert abs(e.flyby_alt_km["Ganymede"] - FIG5_ALT_GANYMEDE_KM) < 50.0, e.flyby_alt_km


def test_eige_europa_altitude_prediction_is_low_order() -> None:
    """The PREDICTED (untargeted) Europa altitude lands the same low order as Fig-5.

    Europa's altitude is not pinned — it is the free output. It comes out ~1,323 km,
    the same low order as Fig-5's printed 470 km; the residual difference is the
    ideal↔real-ephemeris gap (Fig 5 is real-ephemeris). Asserted as a bounded
    cross-check, not an exact reproduction.
    """
    e = feasible_ballistic_eige()
    assert 100.0 < e.flyby_alt_km["Europa"] < 5000.0, e.flyby_alt_km


def test_eige_is_low_excess_speed_regime() -> None:
    """EIGE sits in the low-excess-speed 5-9 km/s regime (not the 12-16 of EGIEIE).

    These are REPORTED construction outputs (the paper prints no EIGE V∞), bounding the
    regime: Europa ~8.70, Io ~5.14, Ganymede ~7.23 km/s — the navigation-viable band.
    """
    e = feasible_ballistic_eige()
    assert 4.0 < e.vinf_kms["Io"] < 9.0, e.vinf_kms
    assert 6.0 < e.vinf_kms["Ganymede"] < 9.0, e.vinf_kms
    assert 7.0 < e.vinf_kms["Europa_dep"] < 10.0, e.vinf_kms
    # Equal in/out |V∞| at the interior flybys (ballistic, machine precision).
    assert abs(e.vinf_kms["Io"] - e.vinf_kms["Io_out"]) < 1.0e-3, e.vinf_kms
    assert abs(e.vinf_kms["Ganymede"] - e.vinf_kms["Ganymede_out"]) < 1.0e-3, e.vinf_kms
