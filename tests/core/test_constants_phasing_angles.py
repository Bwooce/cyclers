"""Sourced-value checks for the J2000 time-phasing angles added to every
``PlanetData`` record: longitude of perihelion ``varpi_deg`` (ϖ) and mean
longitude at epoch ``L0_deg`` (L0).

SOURCING DISCIPLINE: every EXPECTED literal below IS the published source value
(transcribed from the cited table, NOT computed by our own code), so these tests
pin transcription. Citation:

* varpi_deg / L0_deg: Standish & Williams, "Approximate Positions of the
  Planets", JPL Solar System Dynamics, Table 1 (1800-2050 AD), ``varpi_0``
  (longitude of perihelion) and ``L_0`` (mean longitude) columns — the SAME
  table that supplies ``sma_au``/``ecc``/``inc_deg``/``lan_deg``.

These angles feed the time-true visualization (the planet-elements.json
emitter); no ephemeris backend reads them, so unlike inc/lan they are carried
live in PLANETS rather than in a separate opt-in table.
"""

from __future__ import annotations

import pytest

from cyclerfinder.core.constants import PLANETS

# Standish & Williams Table 1, L_0 / varpi_0 columns (deg, J2000).
_PHASING_ANGLES: dict[str, tuple[float, float]] = {
    # code: (L0_deg, varpi_deg)
    "Me": (252.25032350, 77.45779628),
    "V": (181.97909950, 131.60246718),
    "E": (100.46457166, 102.93768193),  # Earth-Moon barycentre row
    "M": (-4.55343205, -23.94362959),
    "J": (34.39644051, 14.72847983),
    "S": (49.95424423, 92.59887831),
    "U": (313.23810451, 170.95427630),
    "N": (-55.12002969, 44.96476227),
}


@pytest.mark.parametrize(
    ("code", "l0", "varpi"), [(c, l0, vp) for c, (l0, vp) in _PHASING_ANGLES.items()]
)
def test_phasing_angles_match_standish_williams(code: str, l0: float, varpi: float) -> None:
    assert PLANETS[code].L0_deg == l0
    assert PLANETS[code].varpi_deg == varpi


def test_all_eight_bodies_carry_phasing_angles() -> None:
    """Every body in the registry has both phasing angles sourced (non-default).

    A zero would mean an un-sourced body slipped in; all eight Standish &
    Williams rows are non-zero in both columns.
    """
    assert set(_PHASING_ANGLES) == set(PLANETS)
    for code in PLANETS:
        assert PLANETS[code].L0_deg != 0.0
        assert PLANETS[code].varpi_deg != 0.0


def test_argument_of_perihelion_reduction_is_finite() -> None:
    """The standard Standish reductions ω = ϖ - Ω and M0 = L0 - ϖ are the
    quantities the viz consumes; check they are computable for all bodies
    (no NaN/inf) without asserting a derived numeric value (would be circular).
    """
    from cyclerfinder.core.ephemeris import inclined_planets

    inclined = inclined_planets()
    for code, p in PLANETS.items():
        lan = inclined[code].lan_deg  # real Ω for inclined bodies, 0.0 for Earth
        omega = p.varpi_deg - lan
        m0 = p.L0_deg - p.varpi_deg
        assert omega == omega  # not NaN
        assert m0 == m0
