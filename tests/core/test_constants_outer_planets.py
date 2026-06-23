"""Sourced-value checks + registry invariants for the heliocentric bodies
added beyond the original V/E/M set: Mercury (``"Me"``), Jupiter (``"J"``),
Saturn (``"S"``), Uranus (``"U"``), Neptune (``"N"``).

SOURCING DISCIPLINE: every EXPECTED literal below IS the published source value
(it is transcribed from the cited table, NOT computed by our own code), so these
tests pin transcription. Citations:

* sma_au / ecc / inc_deg / lan_deg: Standish & Williams, "Approximate Positions
  of the Planets", JPL Solar System Dynamics, Table 1 (1800-2050 AD).
* mu_km3_s2: JPL DE440 (Park et al., AJ 2021). Planet-only for Mercury (no
  moons); system GM for the giant planets (see constants.py module docstring).
* radius_eq_km: IAU 2015 WGCCRE (Archinal et al., CMDA 2018, 130:22).
"""

from __future__ import annotations

import pytest

from cyclerfinder.core.constants import PLANETS, SUPPORTED_BODIES
from cyclerfinder.core.ephemeris import Ephemeris, inclined_planets
from cyclerfinder.search.sequence import Cell, enumerate_cells, tisserand_feasible

_NEW_CODES = ("Me", "J", "S", "U", "N")


# --- sourced spot checks (EXPECTED side == published source) ---------------


@pytest.mark.parametrize(
    ("code", "sma_au", "ecc"),
    [
        # Standish & Williams Table 1, a_0 / e_0 columns.
        ("Me", 0.38709927, 0.20563593),
        ("J", 5.20288700, 0.04838624),
        ("S", 9.53667594, 0.05386179),
        ("U", 19.18916464, 0.04725744),
        ("N", 30.06992276, 0.00859048),
    ],
)
def test_orbital_elements_match_standish_williams(code: str, sma_au: float, ecc: float) -> None:
    assert PLANETS[code].sma_au == sma_au
    assert PLANETS[code].ecc == ecc


@pytest.mark.parametrize(
    ("code", "mu", "radius"),
    [
        # mu: JPL DE440 (Park et al. 2021); radius: IAU 2015 WGCCRE (Archinal 2018).
        ("Me", 2.2031868551e4, 2440.53),
        ("J", 1.267127641e8, 71492.0),
        ("S", 3.79405848418e7, 60268.0),
        ("U", 5.7945564e6, 25559.0),
        ("N", 6.83652710058e6, 24764.0),
    ],
)
def test_mu_and_radius_match_sources(code: str, mu: float, radius: float) -> None:
    assert PLANETS[code].mu_km3_s2 == mu
    assert PLANETS[code].radius_eq_km == radius


@pytest.mark.parametrize(
    ("code", "inc", "lan"),
    [
        # Standish & Williams Table 1, I_0 / Omega_0 columns.
        ("Me", 7.00497902, 48.33076593),
        ("J", 1.30439695, 100.47390909),
        ("S", 2.48599187, 113.66242448),
        ("U", 0.77263783, 74.01692503),
        ("N", 1.77004347, 131.78422574),
    ],
)
def test_inclined_table_matches_standish_williams(code: str, inc: float, lan: float) -> None:
    inclined = inclined_planets()
    assert inclined[code].inc_deg == inc
    assert inclined[code].lan_deg == lan
    # The live coplanar table stays at the 0.0 default (never mutated).
    assert PLANETS[code].inc_deg == 0.0
    assert PLANETS[code].lan_deg == 0.0


# --- registry-shape invariants --------------------------------------------


@pytest.mark.parametrize("code", _NEW_CODES)
def test_registry_record_is_well_formed(code: str) -> None:
    p = PLANETS[code]
    assert p.code == code
    assert p.name  # non-empty full name
    assert p.mu_km3_s2 > 0.0
    assert p.radius_eq_km > 0.0
    assert p.sma_au > 0.0
    assert p.mean_motion_deg_day > 0.0
    assert p.safe_alt_km > 0.0
    assert 0.0 <= p.ecc < 1.0  # bound ellipse


@pytest.mark.parametrize("code", _NEW_CODES)
def test_perihelion_below_aphelion_via_ecc(code: str) -> None:
    p = PLANETS[code]
    peri = p.sma_au * (1.0 - p.ecc)
    apo = p.sma_au * (1.0 + p.ecc)
    assert 0.0 < peri < apo


def test_new_codes_are_supported_bodies() -> None:
    for code in _NEW_CODES:
        assert code in SUPPORTED_BODIES


def test_mean_motion_orders_with_distance() -> None:
    """Sanity: an outer planet orbits slower than an inner one (Kepler III)."""
    mm = {c: PLANETS[c].mean_motion_deg_day for c in ("Me", "E", "J", "N")}
    assert mm["Me"] > mm["E"] > mm["J"] > mm["N"]


# --- byte-identical V/E/M guard -------------------------------------------


def test_vem_entries_unchanged() -> None:
    assert PLANETS["V"].mu_km3_s2 == 3.24858592e5
    assert PLANETS["V"].radius_eq_km == 6051.8
    assert PLANETS["V"].safe_alt_km == 300.0
    assert PLANETS["E"].mu_km3_s2 == 3.98600435507e5
    assert PLANETS["E"].radius_eq_km == 6378.137
    # #426: Earth/Mars flyby floor corrected 300 -> 200 km (SOURCED, Russell 2004 p.165
    # r_p,min 6578.0 / 3598.5 km). Intentional; Venus stays 300 (no sourced revision).
    assert PLANETS["E"].safe_alt_km == 200.0
    assert PLANETS["M"].mu_km3_s2 == 4.282837521e4
    assert PLANETS["M"].radius_eq_km == 3396.19
    assert PLANETS["M"].safe_alt_km == 200.0
    # inc/lan still coplanar default for all three live entries.
    for code in ("V", "E", "M"):
        assert PLANETS[code].inc_deg == 0.0
        assert PLANETS[code].lan_deg == 0.0


# --- behavioural smoke (NO feasibility goldens) ---------------------------


def test_enumerate_and_tisserand_run_for_ej() -> None:
    """enumerate_cells + tisserand_feasible execute on a new body pair without
    error. Behavioural only: we assert the calls run and the predicate returns a
    bool, NOT that any particular cell is feasible."""
    cells = list(enumerate_cells(("E", "J"), l_max=3, k_max=1, n_max=0))
    assert cells  # the generator yields under these caps
    assert all(isinstance(c, Cell) for c in cells)

    coplanar_results = [tisserand_feasible(c, vinf_cap=8.0) for c in cells]
    assert all(isinstance(r, bool) for r in coplanar_results)

    ephem = Ephemeris.inclined_circular()
    inclined_results = [tisserand_feasible(c, vinf_cap=8.0, ephem=ephem) for c in cells]
    assert all(isinstance(r, bool) for r in inclined_results)


def test_enumerate_runs_for_vej_triple() -> None:
    cells = list(enumerate_cells(("V", "E", "J"), l_max=3, k_max=1, n_max=0))
    assert cells
    assert all(set(c.sequence) <= {"V", "E", "J"} for c in cells)
