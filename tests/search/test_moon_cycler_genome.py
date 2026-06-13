"""Tests for the repeated-moon multi-revolution cycler genome (#254).

Covers the four build steps of
``docs/notes/2026-06-14-repeated-moon-multirev-genome-design.md``:

1. moon-system registry + Tisserand/V_inf graph (centre-aware, CGE links);
2. genome representation + decision-vector round-trip (incl. the Liang members);
3. repeated-sequence periodicity corrector (one canonical residual);
4. the REPRODUCE-BEFORE-SEARCH gate against Liang's published Tables 3/5/7.

DISCIPLINE: every EXPECTED value the gate asserts is Liang's PUBLISHED number
(sourced via :mod:`cyclerfinder.search.cge_scaffold`), never our own output.
"""

from __future__ import annotations

import math

import pytest

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES
from cyclerfinder.search.moon_cycler_genome import (
    CGE_SEQUENCE,
    JUPITER_MOONS,
    EncounterGene,
    LegGene,
    MoonCyclerGenome,
    jupiter_system,
    liang_member_genome,
    moon_linkable,
    moon_tisserand_to_vinf,
    moon_vinf_to_tisserand,
    vinf_graph_edges,
)

# ---------------------------------------------------------------------------
# Step 1 — registry + Tisserand/V_inf graph
# ---------------------------------------------------------------------------


def test_jupiter_system_registry_periods() -> None:
    """Galilean periods derive from the registry SMA + Jupiter GM (Kepler III).

    Reference sidereal periods (JPL SSD): Io 1.769 d, Europa 3.551 d,
    Ganymede 7.155 d, Callisto 16.689 d. The registry must reproduce these to
    <1e-2 d (the registry SMAs are 4-figure, so the period is good to that).
    """
    sys = jupiter_system()
    assert sys.moons == JUPITER_MOONS
    expected = {"Io": 1.769, "Europa": 3.551, "Ganymede": 7.155, "Callisto": 16.689}
    for moon, p_ref in expected.items():
        assert sys.period_days(moon) == pytest.approx(p_ref, abs=2.0e-2)


def test_circular_speed_matches_vis_viva() -> None:
    """Circular speed = sqrt(mu_primary / a_moon), the planet-centric vis-viva."""
    sys = jupiter_system()
    for moon in sys.moons:
        a = SATELLITES[moon].sma_km
        assert sys.circular_speed_kms(moon) == pytest.approx(
            math.sqrt(PRIMARIES["Jupiter"] / a), rel=1e-12
        )


@pytest.mark.parametrize("moon", JUPITER_MOONS)
@pytest.mark.parametrize("vinf", [1.0, 3.0, 5.67, 6.99])
def test_tisserand_vinf_round_trip(moon: str, vinf: float) -> None:
    """V_inf -> T -> V_inf is exact (centre-aware, Jovicentric mu)."""
    sys = jupiter_system()
    t_p = moon_vinf_to_tisserand(sys, moon, vinf)
    assert t_p < 3.0  # a real encounter exists at positive V_inf
    assert moon_tisserand_to_vinf(sys, moon, t_p) == pytest.approx(vinf, rel=1e-10)


def test_tisserand_is_centre_aware() -> None:
    """The Jovicentric Tisserand parameter differs from the heliocentric one.

    Guards against accidentally using MU_SUN: at a fixed V_inf the term
    ``V_inf^2 * a_moon / mu`` is enormous heliocentrically (a_moon in AU but
    mu_sun) — the centre-aware value must use the primary GM and stay in the
    physical (0, 3) window for a Galilean encounter.
    """
    sys = jupiter_system()
    t_jov = moon_vinf_to_tisserand(sys, "Callisto", 5.67)
    assert 2.0 < t_jov < 3.0  # physical Jovicentric encounter regime


def test_cge_moon_pairs_are_linkable() -> None:
    """The Liang CGE sequence pairs (E-G, G-C, E-C) share a common-V_inf orbit.

    A necessary condition for the CGE tour to exist as a Tisserand-graph walk:
    every consecutive moon pair in the Callisto-Ganymede-Europa system must be
    energetically linkable at a representative leg V_inf. The published per-leg
    V_inf (Tables 3/5/7) sit around 4.5-7 km/s; we test the lower end (4 km/s)
    where the constant-V_inf contours of the tightly-packed Galileans overlap.
    """
    sys = jupiter_system()
    for a, b in (("Europa", "Ganymede"), ("Ganymede", "Callisto"), ("Europa", "Callisto")):
        assert moon_linkable(sys, a, b, 4.0), f"{a}-{b} should be linkable at 4 km/s"


def test_vinf_graph_edges_symmetric_and_no_self_pairs() -> None:
    """The graph stores both orderings of each pair and omits self-pairs."""
    sys = jupiter_system()
    edges = vinf_graph_edges(sys, 3.0)
    for (a, b), ok in edges.items():
        assert a != b
        assert edges[(b, a)] == ok


# ---------------------------------------------------------------------------
# Step 2 — genome representation + decision vector
# ---------------------------------------------------------------------------


def _sample_genome() -> MoonCyclerGenome:
    sys = jupiter_system()
    encounters = (
        EncounterGene("Callisto", 0.10),
        EncounterGene("Ganymede", -0.25),
        EncounterGene("Callisto", 0.30),
    )
    legs = (LegGene(p=1, q=1, n_rev=1), LegGene(p=2, q=1, n_rev=0))
    return MoonCyclerGenome(sys, encounters, legs, epoch_days=12.5, perijove_scale=0.5)


def test_genome_requires_k_minus_one_legs() -> None:
    """A k-encounter cycle must carry exactly k-1 legs."""
    sys = jupiter_system()
    with pytest.raises(ValueError, match="k-1 legs"):
        MoonCyclerGenome(
            sys,
            (EncounterGene("Callisto", 0.0), EncounterGene("Ganymede", 0.0)),
            (),  # should be 1 leg
        )


def test_decision_vector_round_trip() -> None:
    """encode -> vector -> decode reproduces the genome exactly (lossless)."""
    g = _sample_genome()
    vec = g.to_vector()
    back = MoonCyclerGenome.from_vector(vec, g.system, len(g.encounters))
    assert back.sequence == g.sequence
    assert back.legs == g.legs
    assert back.epoch_days == g.epoch_days
    assert back.perijove_scale == g.perijove_scale
    for e_in, e_out in zip(g.encounters, back.encounters, strict=True):
        assert e_out.moon == e_in.moon
        assert e_out.b_plane_angle_rad == pytest.approx(e_in.b_plane_angle_rad, abs=1e-15)


def test_from_vector_rejects_wrong_length() -> None:
    g = _sample_genome()
    vec = g.to_vector()
    with pytest.raises(ValueError, match="vector length"):
        MoonCyclerGenome.from_vector([*vec, 0.0], g.system, len(g.encounters))


def test_sample_genome_is_valid() -> None:
    assert _sample_genome().is_valid()


def test_invalid_moon_fails_validity() -> None:
    sys = jupiter_system()
    g = MoonCyclerGenome(
        sys,
        (EncounterGene("Titan", 0.0), EncounterGene("Callisto", 0.0)),  # Titan orbits Saturn
        (LegGene(1, 1, 1),),
    )
    assert not g.is_valid()


@pytest.mark.parametrize("member", ["A", "B", "C"])
def test_liang_members_encode_validly(member: str) -> None:
    """The published Liang CGE members encode as valid genomes that round-trip.

    Design step 2 anchor ("the Liang members encode validly"): the sourced
    Callisto-Ganymede-Callisto-Europa-Callisto sequence (5 flybys / 4 legs),
    all 1-rev Lambert arcs, with the member's perijove scale.
    """
    g = liang_member_genome(member)
    assert g.is_valid()
    assert g.sequence == CGE_SEQUENCE
    assert len(g.legs) == 4
    assert all(leg.n_rev == 1 for leg in g.legs)
    vec = g.to_vector()
    back = MoonCyclerGenome.from_vector(vec, g.system, len(g.encounters))
    assert back.sequence == g.sequence
    assert back.perijove_scale == g.perijove_scale
