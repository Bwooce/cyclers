"""Sourced-value checks + registry invariants for the dwarf planets / planetoids
added 2026-06-14 (#260) beyond the eight major planets: Pluto (``"Pl"``),
Ceres (``"Ce"``), Eris (``"Er"``), Makemake (``"Mk"``), Haumea (``"Ha"``),
Vesta (``"Ve"``), Pallas (``"Pa"``).

SOURCING DISCIPLINE: every EXPECTED literal below IS the published source value
(transcribed from the cited table / kernel, NOT computed by our own code), so
these tests pin transcription. Citations:

* mu_km3_s2 (GM):
    - Ceres/Pallas/Vesta/Pluto(system): JPL DE440 ``gm_de440.tpc`` (Park 2021).
    - Eris/Makemake/Haumea: GM = G*M, G = 6.67430e-20 (CODATA 2018); masses
      Eris 1.638e22 kg (Brown & Schaller 2007), Makemake 2.69e21 kg
      (Bamberger 2025), Haumea 3.952e21 kg (Proudfoot 2024). MASS-derived.
* sma_au / ecc / inc_deg / lan_deg:
    - Pluto: Standish & Williams Table 2a (3000 BC-3000 AD; the only Standish
      table retaining Pluto).
    - the rest: JPL SBDB osculating elements, epoch JD 2461200.5, 2026-06-14.
"""

from __future__ import annotations

import math

import pytest

from cyclerfinder.core.constants import PLANETS, SUPPORTED_BODIES, VINF_CEILING_KMS
from cyclerfinder.core.ephemeris import Ephemeris, inclined_planets
from cyclerfinder.search.sequence import Cell, enumerate_cells, tisserand_feasible

_PLANETOID_CODES = ("Pl", "Ce", "Er", "Mk", "Ha", "Ve", "Pa")

# G (CODATA 2018), km^3 kg^-1 s^-2 — the conversion the registry uses for the
# three mass-derived GMs. Mirrored here so the EXPECTED side stays "G*M from the
# published mass", not a copy of our own computed mu.
_G = 6.67430e-20


# --- sourced spot checks (EXPECTED side == published source) ---------------


@pytest.mark.parametrize(
    ("code", "sma_au", "ecc"),
    [
        # Pluto: Standish & Williams Table 2a (a0/e0). Rest: JPL SBDB osculating
        # (a/e) at epoch JD 2461200.5.
        ("Pl", 39.48686035, 0.24885238),
        ("Ce", 2.765552595034094, 0.07969229514816586),
        ("Er", 67.93394687853566, 0.4382385347971672),
        ("Mk", 45.57093317300052, 0.1588889953992523),
        ("Ha", 43.06029023650952, 0.1944430148898797),
        ("Ve", 2.361365965127599, 0.09020374382834395),
        ("Pa", 2.769559010737709, 0.2307000995648547),
    ],
)
def test_orbital_elements_match_sources(code: str, sma_au: float, ecc: float) -> None:
    assert PLANETS[code].sma_au == sma_au
    assert PLANETS[code].ecc == ecc


@pytest.mark.parametrize(
    ("code", "mu"),
    [
        # DE440 gm_de440.tpc (Ceres/Pallas/Vesta/Pluto-system).
        ("Pl", 9.755e2),
        ("Ce", 6.26288886444e1),
        ("Ve", 1.72882328792e1),
        ("Pa", 1.36658781460e1),
    ],
)
def test_de440_gm_matches_kernel(code: str, mu: float) -> None:
    assert PLANETS[code].mu_km3_s2 == mu


@pytest.mark.parametrize(
    ("code", "mass_kg"),
    [
        # GM = G*M from the published satellite-dynamics masses.
        ("Er", 1.638e22),
        ("Mk", 2.69e21),
        ("Ha", 3.952e21),
    ],
)
def test_mass_derived_gm_is_g_times_published_mass(code: str, mass_kg: float) -> None:
    assert PLANETS[code].mu_km3_s2 == pytest.approx(_G * mass_kg, rel=1e-12)


@pytest.mark.parametrize(
    ("code", "inc", "lan"),
    [
        # Pluto: Standish & Williams Table 2a (I0/Omega0). Rest: SBDB i/om.
        ("Pl", 17.14104260, 110.30167986),
        ("Ce", 10.58802780183462, 80.24862682043221),
        ("Er", 43.9258279471791, 36.00477044417249),
        ("Mk", 29.02785603743067, 79.2948338209406),
        ("Ha", 28.20847393040364, 121.7860561329425),
        ("Ve", 7.143925545058711, 103.701293265032),
        ("Pa", 34.93279321851542, 172.8866193357694),
    ],
)
def test_inclined_table_matches_sources(code: str, inc: float, lan: float) -> None:
    inclined = inclined_planets()
    assert inclined[code].inc_deg == inc
    assert inclined[code].lan_deg == lan
    # The live coplanar table stays at the 0.0 default (never mutated).
    assert PLANETS[code].inc_deg == 0.0
    assert PLANETS[code].lan_deg == 0.0


# --- registry-shape invariants --------------------------------------------


@pytest.mark.parametrize("code", _PLANETOID_CODES)
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


@pytest.mark.parametrize("code", _PLANETOID_CODES)
def test_mean_motion_is_kepler_third_law(code: str) -> None:
    """mean_motion is DERIVED from sma + MU_SUN (not hand-copied)."""
    from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, SECONDS_PER_DAY

    p = PLANETS[code]
    a_km = p.sma_au * AU_KM
    period_s = 2.0 * math.pi * math.sqrt(a_km**3 / MU_SUN_KM3_S2)
    expected = 360.0 / (period_s / SECONDS_PER_DAY)
    assert p.mean_motion_deg_day == pytest.approx(expected, rel=1e-12)


def test_codes_are_supported_bodies() -> None:
    for code in _PLANETOID_CODES:
        assert code in SUPPORTED_BODIES


def test_each_planetoid_has_a_finite_vinf_ceiling() -> None:
    """The elliptic-periodicity V_inf ceiling computes for every new body."""
    for code in _PLANETOID_CODES:
        assert VINF_CEILING_KMS[code] > 0.0 and math.isfinite(VINF_CEILING_KMS[code])


def test_no_short_code_collision_with_planets() -> None:
    """The two-letter planetoid codes do not collide with any existing code."""
    seen = list(PLANETS.keys())
    assert len(seen) == len(set(seen))  # all codes unique
    for code in _PLANETOID_CODES:
        assert code in PLANETS


# --- byte-identical major-planet guard ------------------------------------


def test_major_planet_entries_unchanged() -> None:
    """Adding the planetoids must not perturb any of the eight major planets."""
    assert PLANETS["E"].mu_km3_s2 == 3.98600435507e5
    assert PLANETS["E"].sma_au == 1.00000261
    assert PLANETS["N"].mu_km3_s2 == 6.83652710058e6
    assert PLANETS["N"].sma_au == 30.06992276
    # The original eight are still all present.
    for code in ("V", "E", "M", "Me", "J", "S", "U", "N"):
        assert code in PLANETS


# --- behavioural smoke (NO feasibility goldens) ---------------------------


def test_enumerate_and_tisserand_run_for_planetoid_pair() -> None:
    """enumerate_cells + tisserand_feasible execute for an Earth-Ceres pair
    without error. Behavioural only: we assert the calls run and the predicate
    returns a bool, NOT that any particular cell is feasible (the screen is free
    to prune these low-mass bodies)."""
    cells = list(enumerate_cells(("E", "Ce"), l_max=3, k_max=1, n_max=0))
    assert cells
    assert all(isinstance(c, Cell) for c in cells)

    coplanar_results = [tisserand_feasible(c, vinf_cap=8.0) for c in cells]
    assert all(isinstance(r, bool) for r in coplanar_results)

    ephem = Ephemeris.inclined_circular()
    inclined_results = [tisserand_feasible(c, vinf_cap=8.0, ephem=ephem) for c in cells]
    assert all(isinstance(r, bool) for r in inclined_results)
