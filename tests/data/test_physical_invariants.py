"""Task 2: source-independent physical-consistency invariants.

These are golden in the strict sense — pure physics that must hold for any
row regardless of its source, so they catch transcription / unit errors
with no second source needed. The live-catalogue ratchet asserts the
current data/catalogue.yaml passes; the synthetic cases prove each check
has teeth and that the known data realities (multi-arc null a/e, the
sourced ~8% Sanchez period, the VEM beat token, high-V∞ Russell rows, the
Mars-perihelion reach) do NOT false-fail.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import yaml  # type: ignore[import-untyped]

from cyclerfinder.data.validate import validate_physical_invariants

CATALOGUE_PATH = Path(__file__).resolve().parent.parent.parent / "data" / "catalogue.yaml"


def _live_rows() -> list[dict[str, Any]]:
    return cast("list[dict[str, Any]]", yaml.safe_load(CATALOGUE_PATH.read_text()))


# ---------------------------------------------------------------------------
# Live catalogue ratchet
# ---------------------------------------------------------------------------


def test_live_catalogue_passes_physical_invariants() -> None:
    """The real catalogue must satisfy every source-independent physics check."""
    errs = validate_physical_invariants(_live_rows())
    assert errs == [], "physical-invariant violations:\n" + "\n".join(errs)


# ---------------------------------------------------------------------------
# orbit-element identities (teeth)
# ---------------------------------------------------------------------------


def test_a_must_equal_peri_apo_midpoint() -> None:
    bad = {
        "id": "x",
        "cycler_class": "single-ellipse",
        "orbit_elements": {"a_au": 1.6, "e": 0.393, "perihelion_au": 0.97, "aphelion_au": 3.0},
    }
    errs = validate_physical_invariants([bad])
    assert any("(peri+apo)/2" in m for m in errs), errs


def test_e_must_match_geometry() -> None:
    bad = {
        "id": "x",
        "cycler_class": "single-ellipse",
        "orbit_elements": {"a_au": 1.6, "e": 0.10, "perihelion_au": 0.97, "aphelion_au": 2.23},
    }
    errs = validate_physical_invariants([bad])
    assert any("apo-peri" in m for m in errs), errs


def test_eccentricity_out_of_range_flagged() -> None:
    bad = {"id": "x", "cycler_class": "single-ellipse", "orbit_elements": {"a_au": 1.6, "e": 1.4}}
    errs = validate_physical_invariants([bad])
    assert any("eccentricity out of" in m for m in errs), errs


def test_peri_a_apo_ordering_flagged() -> None:
    bad = {
        "id": "x",
        "cycler_class": "single-ellipse",
        # a below perihelion: ordering violated (and a/mid mismatch).
        "orbit_elements": {"a_au": 0.5, "e": 0.257, "perihelion_au": 0.97, "aphelion_au": 1.64},
    }
    errs = validate_physical_invariants([bad])
    assert any("peri <= a <= apo" in m for m in errs), errs


def test_clean_single_ellipse_orbit_passes() -> None:
    ok = {
        "id": "aldrin",
        "cycler_class": "single-ellipse",
        "primary": "Sun",
        "bodies": ["E", "M"],
        "orbit_elements": {"a_au": 1.6, "e": 0.393, "perihelion_au": 0.97, "aphelion_au": 2.23},
    }
    assert validate_physical_invariants([ok]) == []


# ---------------------------------------------------------------------------
# V_inf sanity (teeth + no false-fail on real high-V∞ rows)
# ---------------------------------------------------------------------------


def test_vinf_unit_error_flagged() -> None:
    """A V∞ entered in m/s (9700) lands well above the km/s ceiling."""
    bad = {"id": "x", "vinf_kms_at_encounters": [{"body": "E", "vinf_kms": 9700.0}]}
    errs = validate_physical_invariants([bad])
    assert any("unit error" in m for m in errs), errs


def test_vinf_negative_flagged() -> None:
    bad = {"id": "x", "vinf_kms_at_encounters": [{"body": "E", "vinf_kms": -1.0}]}
    errs = validate_physical_invariants([bad])
    assert any("negative" in m for m in errs), errs


def test_high_but_real_vinf_passes() -> None:
    """Real Russell-Ocampo cyclers reach ~20 km/s — must NOT false-fail."""
    ok = {"id": "russell", "vinf_kms_at_encounters": [{"body": "E", "vinf_kms": 20.3}]}
    assert validate_physical_invariants([ok]) == []


# ---------------------------------------------------------------------------
# Period identity: single-ellipse body pair, multi-arc skip, beat token
# ---------------------------------------------------------------------------


def test_single_ellipse_period_mismatch_flagged() -> None:
    bad = {
        "id": "x",
        "cycler_class": "single-ellipse",
        "bodies": ["E", "M"],
        "period": {"pair": "E-M", "k": 1, "years": 5.0},  # k=1 E-M ~ 2.135 yr
    }
    errs = validate_physical_invariants([bad])
    assert any("T_syn" in m for m in errs), errs


def test_single_ellipse_period_match_passes() -> None:
    ok = {
        "id": "x",
        "cycler_class": "single-ellipse",
        "bodies": ["E", "M"],
        "period": {"pair": "E-M", "k": 1, "years": 2.135},
    }
    assert validate_physical_invariants([ok]) == []


def test_multi_arc_approximate_period_label_not_flagged() -> None:
    """multi-arc synodic-integer k is an approximate label, not the true period.

    Mirrors sourced sanchez-net-2022-em-cycler2 (7.87 yr vs k=4 fit 8.54 yr,
    ~8% under) — must NOT be flagged as an error.
    """
    ok = {
        "id": "sanchez-like",
        "cycler_class": "multi-arc",
        "bodies": ["E", "M"],
        "period": {"pair": "E-M", "k": 4, "years": 7.87},
    }
    assert validate_physical_invariants([ok]) == []


def test_period_years_nonpositive_flagged() -> None:
    bad = {"id": "x", "period": {"pair": "E-M", "k": 1, "years": 0.0}}
    errs = validate_physical_invariants([bad])
    assert any("years must be > 0" in m for m in errs), errs


# ---------------------------------------------------------------------------
# Beat-token rows (VEM-syn) — the Forge R1 delta-4 / M8 R1 delta-3 trap
# ---------------------------------------------------------------------------


def test_vem_beat_token_rows_pass_in_every_itinerary_order() -> None:
    """All three VEM beat-token orderings validate against the canonical beat.

    The raw body order varies (V-E-M, E-M-V, M-E-V); the invariant must
    canonicalise (sort) the body set before the beat search, else
    multi_body_beat_days returns [] for the middle-reference orderings and
    the rows false-fail (or crash on a naive 'VEM'/'syn' pair split).
    """
    for bodies in (["V", "E", "M"], ["E", "M", "V"], ["M", "E", "V"]):
        row = {
            "id": "vem-" + "".join(bodies),
            "cycler_class": "multi-arc",
            "bodies": bodies,
            "period": {"pair": "VEM-syn", "k": 2, "years": 12.8},
        }
        errs = validate_physical_invariants([row])
        assert errs == [], f"VEM order {bodies} false-failed: {errs}"


def test_vem_beat_token_mismatch_flagged() -> None:
    """A wrong beat-token period IS flagged (the check has teeth)."""
    row = {
        "id": "vem-bad",
        "cycler_class": "multi-arc",
        "bodies": ["V", "E", "M"],
        "period": {"pair": "VEM-syn", "k": 2, "years": 99.0},
    }
    errs = validate_physical_invariants([row])
    assert any("beat(VEM-syn)" in m for m in errs), errs


def test_beat_token_does_not_crash_on_non_body_pair() -> None:
    """A beat token must never be split into synodic_period('VEM','syn')."""
    row = {
        "id": "vem",
        "cycler_class": "multi-arc",
        "bodies": ["V", "E", "M"],
        "period": {"pair": "VEM-syn", "k": 2, "years": 12.8},
    }
    # Must not raise.
    validate_physical_invariants([row])


# ---------------------------------------------------------------------------
# Encounter reach: Mars-perihelion rows must NOT false-fail
# ---------------------------------------------------------------------------


def test_reach_mars_near_perihelion_passes() -> None:
    """niehoff-visit1-like: apo 1.40 < Mars sma 1.524 but > Mars peri 1.381.

    Using mean sma would false-fail; the correct invariant uses Mars
    perihelion (e≈0.093) and passes.
    """
    ok = {
        "id": "niehoff-like",
        "cycler_class": "single-ellipse",
        "primary": "Sun",
        "bodies": ["E", "M"],
        "orbit_elements": {"a_au": 1.17, "e": 0.197, "perihelion_au": 0.94, "aphelion_au": 1.40},
    }
    assert validate_physical_invariants([ok]) == []


def test_reach_orbit_too_small_for_mars_flagged() -> None:
    """An orbit whose aphelion is below even Mars perihelion cannot reach Mars."""
    bad = {
        "id": "too-small",
        "cycler_class": "single-ellipse",
        "primary": "Sun",
        "bodies": ["E", "M"],
        "orbit_elements": {"a_au": 1.05, "e": 0.10, "perihelion_au": 0.945, "aphelion_au": 1.155},
    }
    errs = validate_physical_invariants([bad])
    assert any("cannot encounter M" in m for m in errs), errs


# ---------------------------------------------------------------------------
# Known data realities that must be skipped, not flagged
# ---------------------------------------------------------------------------


def test_multi_arc_null_ae_not_flagged() -> None:
    """multi-arc rows have null a/e by design — no orbit-identity check runs."""
    ok = {"id": "x", "cycler_class": "multi-arc", "orbit_elements": {"center": "Sun"}}
    assert validate_physical_invariants([ok]) == []
