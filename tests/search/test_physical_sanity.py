"""Tests for the physical-sanity flyby gate (task #324, 2026-06-16).

Motivating case: the Umbriel-Oberon-Umbriel (1,1) SILVER survivor from #312
(commit ``bebaeaf`` / ``d0d5898``) passed every automated guard but delivered
no usable bend at Umbriel. The agent's manual physical check rejected it; this
test suite locks the equivalent automated check in place.

EXPECTED values are computed from the patched-conic max-bend formula
``sin(delta_max/2) = mu / (mu + r_p * V_inf**2)`` with sourced GM / radius for
each body. They are NOT golden values in the sourced-publication sense
(:doc:`feedback_golden_tests_sourced_only`) — the gate is a derived consequence
of :func:`cyclerfinder.core.flyby.max_bend` (which IS sourced from BMW §6.4)
plus a judgment threshold. So the assertions here check *consistency* with
``max_bend`` plus *real-mission-shaped magnitudes* (Galileo Earth, Cassini
Venus, Aldrin Mars — admit; Umbriel-at-2.27-km/s — reject).
"""

from __future__ import annotations

from math import degrees

import pytest

from cyclerfinder.core.flyby import max_bend
from cyclerfinder.core.satellites import SATELLITES
from cyclerfinder.search.physical_sanity import (
    DEFAULT_MIN_USEFUL_BEND_DEG,
    FlybyPhysicalVerdict,
    candidate_passes_physical_gate,
    flyby_is_useful,
)

# ---------------------------------------------------------------------------
# Primary motivating case: Umbriel at the prompt's worst-case V_inf=2.27 km/s
# ---------------------------------------------------------------------------


def test_umbriel_at_2p27_rejected() -> None:
    """At V_inf=2.27 km/s the Umbriel flyby's max bend is ~2.7 deg << 5 deg → reject.

    The prompt's motivating case (the agent's manual rejection of an
    Umbriel-Oberon-Umbriel candidate where V_inf at Umbriel sat well above
    Umbriel's escape velocity ~0.55 km/s). The patched-conic max bend at
    r_p = R_Umbriel + 100 km is 2.6991 deg; the 5 deg floor rejects.
    """
    verdict = flyby_is_useful("Umbriel", 2.27, min_safe_altitude_km=100.0)
    assert isinstance(verdict, FlybyPhysicalVerdict)
    assert verdict.body == "Umbriel"
    assert verdict.vinf_kms == pytest.approx(2.27)
    assert verdict.min_safe_altitude_km == pytest.approx(100.0)
    # R_Umbriel = 584.7 km from satellites.py:165 → r_p = 684.7 km.
    assert verdict.periapsis_radius_km == pytest.approx(684.7)
    # Sanity vs the formula in core/flyby.py.
    expected_bend = degrees(max_bend(SATELLITES["Umbriel"].mu_km3_s2, 684.7, 2.27))
    assert verdict.max_bend_deg == pytest.approx(expected_bend)
    assert verdict.max_bend_deg == pytest.approx(2.699, abs=0.01)
    assert verdict.is_useful is False
    assert "below" in verdict.notes


def test_umbriel_actual_312_silver_admitted() -> None:
    """The #312 SILVER row's actual V_inf~0.92 km/s at Umbriel ADMITS (~14.7 deg).

    The prompt cited the prompt's worst-case 2.27 km/s; the *actual* JSONL row
    (data/scan_312_uranus_oberon_umbriel.jsonl, only SILVER survivor) has
    V_inf=0.9199 km/s at Umbriel, which the gate ADMITS at ~14.7 deg. This
    is the right behaviour — the gate is a floor, not a ceiling, and the
    actual #312 candidate fails for OTHER reasons (genome ceiling, not
    geometric vacuity).
    """
    verdict = flyby_is_useful("Umbriel", 0.9199, min_safe_altitude_km=100.0)
    assert verdict.is_useful is True
    assert verdict.max_bend_deg == pytest.approx(14.71, abs=0.05)


# ---------------------------------------------------------------------------
# Sourced real-mission anchors (Galileo, Cassini, Aldrin)
# ---------------------------------------------------------------------------


def test_galileo_earth_flyby_admitted() -> None:
    """Galileo's Earth flybys (V_inf~6.232 km/s) sit at ~74 deg max bend → admit.

    V_inf=6.232 km/s is the Galileo-1 Earth-encounter value (D'Amario,
    Bright & Wolf 1992, Space Sci Rev 60:23, Table 2). At Earth's sourced 200 km
    flyby floor (#426, Russell 2004 p.165) the max bend is ~75.1 deg — comfortably
    above any sensible useful-bend floor.
    """
    verdict = flyby_is_useful("E", 6.232)
    assert verdict.is_useful is True
    assert verdict.max_bend_deg == pytest.approx(75.09, abs=0.1)


def test_aldrin_mars_flyby_admitted() -> None:
    """Aldrin-cycler Mars flybys (V_inf~5.5 km/s) sit at ~32 deg max bend → admit.

    The classic Aldrin 1L1 cycler's Mars encounter V_inf is ~5.5 km/s
    (Russell & Ocampo 2005; McConaghy 2002 Table 4 lists Earth V_inf=8 km/s;
    the Mars side runs cooler). Max bend at Mars (sourced 200 km floor, #426) is ~32.8 deg.
    """
    verdict = flyby_is_useful("M", 5.5)
    assert verdict.is_useful is True
    assert verdict.max_bend_deg == pytest.approx(32.82, abs=0.1)


def test_cassini_venus_flyby_admitted() -> None:
    """Cassini's Venus flybys (V_inf~7 km/s) sit at ~61 deg max bend → admit.

    Cassini's V1/V2 Venus flybys were used as the inner-leg gravity assists
    of the VVEJGA trajectory (Peralta & Flanagan 1995). V_inf~7 km/s; max
    bend at the 300 km safe altitude is ~61 deg.
    """
    verdict = flyby_is_useful("V", 7.0)
    assert verdict.is_useful is True
    assert verdict.max_bend_deg == pytest.approx(61.42, abs=0.1)


# ---------------------------------------------------------------------------
# Sequence-level (the actual gate the discovery daemon will call)
# ---------------------------------------------------------------------------


def test_sequence_umbriel_oberon_umbriel_with_high_vinf_rejected() -> None:
    """Worst-case Umbriel-Oberon-Umbriel @ V_inf=(2.27,0.98,2.27) — reject on Umbriel.

    The simulated worst-case prompt sequence — sequence-level rejection happens
    on the first Umbriel pass.
    """
    passed, verdicts = candidate_passes_physical_gate(
        ("Umbriel", "Oberon", "Umbriel"),
        (2.27, 0.98, 2.27),
    )
    assert passed is False
    assert verdicts[0].is_useful is False  # Umbriel#1 fails
    assert verdicts[1].is_useful is True  # Oberon OK
    assert verdicts[2].is_useful is False  # Umbriel#2 fails


def test_sequence_actual_312_silver_admitted() -> None:
    """The actual #312 SILVER row admits — gate doesn't double-flag a real candidate.

    Per the JSONL: vinfs=(0.9199, 0.9604, 0.8947) for (Umbriel, Oberon, Umbriel).
    All three sit above the 5 deg useful-bend floor, so the gate admits — the
    candidate is correctly forwarded to lit-check / ML / gauntlet for OTHER
    judgments. The gate's job is to catch the GEOMETRIC vacuum case, not to
    validate the candidate.
    """
    passed, verdicts = candidate_passes_physical_gate(
        ("Umbriel", "Oberon", "Umbriel"),
        (0.9199258810725036, 0.9604309791298091, 0.8946936085078939),
    )
    assert passed is True
    assert all(v.is_useful for v in verdicts)


def test_sequence_galileo_veega_admitted() -> None:
    """A Galileo-shaped V-E-E-G-A profile admits at every leg."""
    # Approximate published V_inf magnitudes (D'Amario et al. 1992); these are
    # *real* tour V_infs, all well above the 5 deg floor at their respective bodies.
    passed, verdicts = candidate_passes_physical_gate(
        ("V", "E", "E"),
        (4.85, 6.232, 9.0),
    )
    assert passed is True
    assert all(v.is_useful for v in verdicts)


# ---------------------------------------------------------------------------
# Threshold + altitude knob behaviour
# ---------------------------------------------------------------------------


def test_default_threshold_is_5_deg() -> None:
    assert DEFAULT_MIN_USEFUL_BEND_DEG == 5.0


def test_threshold_tightenable() -> None:
    """A tighter (e.g. 30 deg) threshold can reject otherwise-OK flybys.

    Mars at V_inf=5.5 km/s gives ~32 deg max bend — admits at 5 and 30 deg
    floors, rejects at 35 deg. Confirms the threshold is genuinely a knob.
    """
    v5 = flyby_is_useful("M", 5.5, min_useful_bend_deg=5.0)
    v30 = flyby_is_useful("M", 5.5, min_useful_bend_deg=30.0)
    v35 = flyby_is_useful("M", 5.5, min_useful_bend_deg=35.0)
    assert v5.is_useful is True
    assert v30.is_useful is True
    assert v35.is_useful is False


def test_altitude_override_changes_periapsis() -> None:
    """A higher safe-altitude override widens r_p → reduces max_bend."""
    low = flyby_is_useful("Umbriel", 0.92, min_safe_altitude_km=100.0)
    high = flyby_is_useful("Umbriel", 0.92, min_safe_altitude_km=1000.0)
    assert high.periapsis_radius_km > low.periapsis_radius_km
    assert high.max_bend_deg < low.max_bend_deg


def test_per_body_altitude_override_in_sequence() -> None:
    overrides = {"Umbriel": 1000.0}
    passed_low, verdicts_low = candidate_passes_physical_gate(("Umbriel", "Oberon"), (0.92, 0.98))
    passed_high, verdicts_high = candidate_passes_physical_gate(
        ("Umbriel", "Oberon"),
        (0.92, 0.98),
        per_body_min_safe_altitude_km=overrides,
    )
    # Oberon unchanged.
    assert verdicts_low[1].max_bend_deg == pytest.approx(verdicts_high[1].max_bend_deg)
    # Umbriel widened → smaller bend.
    assert verdicts_high[0].max_bend_deg < verdicts_low[0].max_bend_deg
    # Both still admit at default floor.
    assert passed_low is True
    assert passed_high is True


# ---------------------------------------------------------------------------
# Body lookup + error handling
# ---------------------------------------------------------------------------


def test_unknown_body_raises() -> None:
    """The gate must NEVER silently admit an unknown body."""
    with pytest.raises(KeyError):
        flyby_is_useful("NotARealBody", 1.0)


def test_planet_code_lookup_works() -> None:
    """Single-letter planet codes resolve via PLANETS first."""
    v = flyby_is_useful("E", 6.232)
    assert v.body == "E"
    assert v.periapsis_radius_km == pytest.approx(6378.137 + 200.0)  # #426 sourced floor


def test_moon_name_lookup_works() -> None:
    """Full moon names resolve via SATELLITES."""
    v = flyby_is_useful("Europa", 5.0)
    assert v.body == "Europa"
    assert v.periapsis_radius_km == pytest.approx(1560.8 + 100.0)


def test_negative_vinf_rejected() -> None:
    with pytest.raises(ValueError):
        flyby_is_useful("E", -1.0)


def test_negative_threshold_rejected() -> None:
    with pytest.raises(ValueError):
        flyby_is_useful("E", 5.0, min_useful_bend_deg=-1.0)


def test_mismatched_sequence_length_rejected() -> None:
    with pytest.raises(ValueError):
        candidate_passes_physical_gate(("E", "V"), (5.0,))


# ---------------------------------------------------------------------------
# Direct consistency vs core/flyby.py::max_bend (the wrapped primitive)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "body, vinf",
    [
        ("E", 6.232),
        ("V", 7.0),
        ("M", 5.5),
        ("Umbriel", 2.27),
        ("Umbriel", 0.92),
        ("Oberon", 0.98),
        ("Europa", 5.0),
        ("Titan", 3.5),
    ],
)
def test_max_bend_deg_matches_core_flyby(body: str, vinf: float) -> None:
    """The reported max_bend_deg must equal core/flyby.py's max_bend (in deg)."""
    verdict = flyby_is_useful(body, vinf)
    if body in {"E", "V", "M"}:
        from cyclerfinder.core.constants import PLANETS

        p = PLANETS[body]
        mu, r = p.mu_km3_s2, p.radius_eq_km
    else:
        s = SATELLITES[body]
        mu, r = s.mu_km3_s2, s.radius_eq_km
    expected = degrees(max_bend(mu, r + verdict.min_safe_altitude_km, vinf))
    assert verdict.max_bend_deg == pytest.approx(expected, rel=1e-12)
