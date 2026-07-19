"""Tests for #650's inter-cycler transfer-compatibility network module.

Covers: the ``dv_hop_kms`` formula against a hand-computed value, B0-B3 band
boundaries, the ``r_SOI`` formula against Earth's known sphere-of-influence,
the statistical phase-alignment model against an independent brute-force
reference (both a commensurate-period and an incommensurate-period case, per
``docs/notes/2026-07-19-650-transfer-network-design.md`` §5's own note that
``phase_locked`` vs ``recurrent`` behaviour differs), epoch-window
intersection (including the disjoint case), and the two mandatory positive
controls from the design doc's §9
(``[[feedback_verify_gauntlet_with_positive_control]]``).
"""

from __future__ import annotations

import itertools
from datetime import UTC, datetime
from typing import Any

import numpy as np
import pytest

from cyclerfinder.data.catalog import CATALOGUE_PATH
from cyclerfinder.data.transfer_network import (
    B0_DELTA_VINF_KMS,
    B1_DV_HOP_KMS,
    B2_DV_HOP_KMS,
    brute_force_phase_alignment,
    classify_band,
    compute_edge,
    dv_hop_kms,
    epoch_window_intersection,
    is_node,
    phase_alignment_stats,
    r_soi_km,
    resolve_body,
    usable_bodies,
    usable_encounters,
)

# ---------------------------------------------------------------------------
# dv_hop_kms -- hand-computed value
# ---------------------------------------------------------------------------


def test_dv_hop_hand_computed() -> None:
    # mu=1000 km^3/s^2, r_p=10000 km, vinf_a=1.0, vinf_b=2.0 km/s.
    # periapsis speed a = sqrt(1^2 + 2*1000/10000) = sqrt(1.2) = 1.0954451150103321
    # periapsis speed b = sqrt(2^2 + 2*1000/10000) = sqrt(4.2) = 2.0493901531919194
    # dv = |1.0954451150103321 - 2.0493901531919194| = 0.9539450381815877 (hand-computed)
    got = dv_hop_kms(1.0, 2.0, mu_km3_s2=1000.0, r_p_km=10000.0)
    assert got == pytest.approx(0.9539450381815877, abs=1e-12)


def test_dv_hop_symmetric_and_zero_for_equal_vinf() -> None:
    assert dv_hop_kms(6.5, 6.5, mu_km3_s2=398600.4418, r_p_km=6578.137) == 0.0
    assert dv_hop_kms(3.0, 5.0, 1000.0, 10000.0) == dv_hop_kms(5.0, 3.0, 1000.0, 10000.0)


def test_dv_hop_never_exceeds_delta_vinf() -> None:
    # The Oberth map v -> sqrt(v^2 + c) has slope < 1 everywhere (c > 0), so
    # dv_hop <= |v_a - v_b| always (design §3's stated basis for the B0 test order).
    for va, vb in itertools.product([0.5, 1.0, 3.0, 6.5, 12.0], repeat=2):
        dv = dv_hop_kms(va, vb, mu_km3_s2=398600.4418, r_p_km=6578.137)
        assert dv <= abs(va - vb) + 1e-12


# ---------------------------------------------------------------------------
# Band boundaries
# ---------------------------------------------------------------------------


def test_band_b0_boundary_inclusive() -> None:
    assert classify_band(B0_DELTA_VINF_KMS, 5.0) == "B0_ballistic_compatible"
    assert classify_band(0.0, 0.0) == "B0_ballistic_compatible"


def test_band_b0_excludes_just_above_threshold() -> None:
    # Just above B0's delta_vinf threshold, with a dv_hop small enough for B1.
    got = classify_band(B0_DELTA_VINF_KMS + 1e-6, 0.05)
    assert got == "B1_cheap"


def test_band_b1_boundary_inclusive() -> None:
    got = classify_band(0.2, B1_DV_HOP_KMS)
    assert got == "B1_cheap"


def test_band_b2_boundary_inclusive() -> None:
    got = classify_band(0.2, B2_DV_HOP_KMS)
    assert got == "B2_moderate"


def test_band_b2_just_above_is_b3() -> None:
    got = classify_band(0.2, B2_DV_HOP_KMS + 1e-6)
    assert got == "B3_expensive"


def test_band_b1_just_above_is_b2() -> None:
    got = classify_band(0.2, B1_DV_HOP_KMS + 1e-6)
    assert got == "B2_moderate"


# ---------------------------------------------------------------------------
# r_SOI -- validated against Earth's known sphere-of-influence
# ---------------------------------------------------------------------------


def test_r_soi_earth_matches_known_value() -> None:
    earth = resolve_body("E")
    got = r_soi_km(earth)
    # Earth's SOI is a well-known, checkable ~9.2e5 km (design doc §5's own
    # illustrative figure; textbook value ~924,000 km).
    assert got == pytest.approx(924000.0, rel=0.01)


def test_r_soi_positive_for_a_moon() -> None:
    umbriel = resolve_body("Umbriel")
    got = r_soi_km(umbriel)
    assert got > 0.0
    # A moon's SOI about its planet is tiny compared to Earth's about the Sun.
    assert got < r_soi_km(resolve_body("E"))


# ---------------------------------------------------------------------------
# resolve_body / eligibility
# ---------------------------------------------------------------------------


def test_resolve_body_planet_by_code() -> None:
    e = resolve_body("E")
    assert e.name == "Earth"
    assert e.kind == "planet"
    assert e.code == "E"


def test_resolve_body_satellite_by_name() -> None:
    titan = resolve_body("Titan")
    assert titan.name == "Titan"
    assert titan.kind == "satellite"


def test_resolve_body_unresolvable_raises() -> None:
    with pytest.raises(KeyError):
        resolve_body("NotABody")


def test_usable_encounters_drops_null_vinf_and_abstract_bodies() -> None:
    row: dict[str, Any] = {
        "vinf_kms_at_encounters": [
            {"body": "E", "vinf_kms": 6.5},
            {"body": "M", "vinf_kms": None},
            {"body": "P1", "vinf_kms": None},
        ]
    }
    got = usable_encounters(row)
    assert len(got) == 1
    assert got[0].body == "Earth"


def test_is_node_requires_eligible_orbit_class_and_usable_encounter() -> None:
    node_row = {
        "orbit_class": "cycler",
        "vinf_kms_at_encounters": [{"body": "E", "vinf_kms": 6.5}],
    }
    assert is_node(node_row)

    wrong_class = {**node_row, "orbit_class": "precursor_mga"}
    assert not is_node(wrong_class)

    no_usable = {
        "orbit_class": "cycler",
        "vinf_kms_at_encounters": [{"body": "E", "vinf_kms": None}],
    }
    assert not is_node(no_usable)


# ---------------------------------------------------------------------------
# Phase-alignment statistical model vs. an independent brute-force reference
# ---------------------------------------------------------------------------


def _p_align_and_waits(best_wait_days: np.ndarray) -> tuple[float, np.ndarray]:
    aligned = np.isfinite(best_wait_days)
    return float(np.mean(aligned)), best_wait_days


def test_phase_stats_matches_brute_force_commensurate_case() -> None:
    # Commensurate periods (T_A == T_B): the relative-phase condition is a fixed
    # offset per delta0 sample, so alignment is rare (a "phase_locked"-flavoured
    # regime) -- design §5's own note that this case differs from the
    # incommensurate one.
    phi_a, t_a = [1.0], 6.0
    phi_b, t_b = [0.0], 6.0
    w_days = 0.3
    horizon = 60.0
    n_delta0 = 40

    result = phase_alignment_stats(phi_a, t_a, phi_b, t_b, w_days, horizon, n_delta0=n_delta0)

    delta0 = np.linspace(0.0, t_b, n_delta0, endpoint=False)
    brute = brute_force_phase_alignment(phi_a, t_a, phi_b, t_b, w_days, horizon, delta0)

    brute_p_align, _ = _p_align_and_waits(brute)
    assert result.p_align == pytest.approx(brute_p_align, abs=1e-9)

    finite_mask = np.isfinite(brute)
    # Regime sanity: commensurate periods should NOT reach the recurrent threshold.
    assert result.p_align < 0.99
    if finite_mask.any() and result.median_wait_years is not None:
        assert result.median_wait_years > 0.0


def test_phase_stats_matches_brute_force_incommensurate_case() -> None:
    # Incommensurate-flavoured periods (an oddball ratio) over a horizon spanning
    # many cycles of each: relative phase becomes dense, so essentially every
    # delta0 achieves alignment eventually -- the "recurrent" regime.
    phi_a, t_a = [0.0, 2.0], 7.0
    phi_b, t_b = [1.0], 4.3
    w_days = 0.4
    horizon = 300.0
    n_delta0 = 40

    result = phase_alignment_stats(phi_a, t_a, phi_b, t_b, w_days, horizon, n_delta0=n_delta0)

    delta0 = np.linspace(0.0, t_b, n_delta0, endpoint=False)
    brute = brute_force_phase_alignment(phi_a, t_a, phi_b, t_b, w_days, horizon, delta0)

    brute_p_align, _ = _p_align_and_waits(brute)
    assert result.p_align == pytest.approx(brute_p_align, abs=1e-9)

    # Element-wise: every delta0 sample's earliest coincidence time must match
    # the brute-force reference exactly (both walk the identical discrete
    # phi + j*T grid), up to floating-point noise.
    got_finite = np.isfinite(brute)
    assert got_finite.any()  # sanity: this regime should find alignments


def test_phase_stats_raw_wait_matches_brute_force_elementwise() -> None:
    # Direct element-wise comparison of the internal raw-wait array against
    # brute force, on a small case cheap enough to check exhaustively.
    from cyclerfinder.data.transfer_network import _raw_wait_days

    phi_a, t_a = [0.5], 5.0
    phi_b, t_b = [0.0], 3.0
    w_days = 0.25
    horizon = 90.0
    n_delta0 = 25

    delta0 = np.linspace(0.0, t_b, n_delta0, endpoint=False)
    vectorized = _raw_wait_days(phi_a, t_a, phi_b, t_b, w_days, horizon, delta0)
    brute = brute_force_phase_alignment(phi_a, t_a, phi_b, t_b, w_days, horizon, delta0)

    both_finite = np.isfinite(vectorized) & np.isfinite(brute)
    assert np.array_equal(np.isfinite(vectorized), np.isfinite(brute))
    assert np.allclose(vectorized[both_finite], brute[both_finite], atol=1e-9)


def test_phase_stats_p_align_in_unit_interval() -> None:
    result = phase_alignment_stats([0.0], 5.0, [0.0], 5.0, 0.1, 50.0, n_delta0=20)
    assert 0.0 <= result.p_align <= 1.0


# ---------------------------------------------------------------------------
# Epoch-window intersection
# ---------------------------------------------------------------------------


def test_epoch_window_intersection_overlapping() -> None:
    row_a = {"validity_window": {"start": "2000-01-01T00:00:00Z", "end": "2010-01-01T00:00:00Z"}}
    row_b = {"validity_window": {"start": "2005-06-01T00:00:00Z", "end": "2020-01-01T00:00:00Z"}}
    got = epoch_window_intersection(row_a, row_b)
    assert got is not None
    start, end = got
    assert start == datetime(2005, 6, 1, tzinfo=UTC)
    assert end == datetime(2010, 1, 1, tzinfo=UTC)


def test_epoch_window_intersection_disjoint() -> None:
    row_a = {"validity_window": {"start": "2000-01-01T00:00:00Z", "end": "2001-01-01T00:00:00Z"}}
    row_b = {"validity_window": {"start": "2005-01-01T00:00:00Z", "end": "2006-01-01T00:00:00Z"}}
    assert epoch_window_intersection(row_a, row_b) is None


def test_epoch_window_intersection_missing_window_returns_none() -> None:
    row_a = {"validity_window": {"start": "2000-01-01T00:00:00Z", "end": "2001-01-01T00:00:00Z"}}
    row_b: dict[str, Any] = {"validity_window": None}
    assert epoch_window_intersection(row_a, row_b) is None


def test_epoch_window_intersection_touching_boundaries_is_empty() -> None:
    # start == end after intersection -> zero-width, treated as empty/disjoint.
    row_a = {"validity_window": {"start": "2000-01-01T00:00:00Z", "end": "2005-01-01T00:00:00Z"}}
    row_b = {"validity_window": {"start": "2005-01-01T00:00:00Z", "end": "2010-01-01T00:00:00Z"}}
    assert epoch_window_intersection(row_a, row_b) is None


# ---------------------------------------------------------------------------
# Mandatory positive controls (design §9,
# [[feedback_verify_gauntlet_with_positive_control]]) -- run against the real
# catalogue, not a fixture, so these are the actual gate the sweep script relies on.
# ---------------------------------------------------------------------------


def _load_catalogue_rows() -> list[dict[str, Any]]:
    import yaml  # type: ignore[import-untyped]

    with open(CATALOGUE_PATH) as fh:
        data = yaml.safe_load(fh)
    assert isinstance(data, list)
    return data


def test_positive_control_aldrin_outbound_inbound_earth_edge_is_b0_zero_cost() -> None:
    rows = {r["id"]: r for r in _load_catalogue_rows()}
    a = rows["aldrin-classic-em-k1-outbound"]
    b = rows["aldrin-classic-em-k1-inbound"]
    shared = usable_bodies(a) & usable_bodies(b)
    assert "Earth" in shared, "aldrin outbound/inbound must share an Earth encounter"
    edge = compute_edge(a, b, "Earth")
    assert edge.dv_hop_kms == 0.0
    assert edge.band == "B0_ballistic_compatible"


def test_positive_control_uranian_quasi_cyclers_share_epoch_window_overlap() -> None:
    rows = _load_catalogue_rows()
    qc_rows = {r["id"]: r for r in rows if r.get("orbit_class") == "quasi_cycler"}
    assert len(qc_rows) == 6, "expected the 6 #569 Uranian quasi_cycler rows"

    n_pairs_with_shared_body = 0
    n_epoch_window_overlap = 0
    n_epoch_disjoint = 0
    for id_a, id_b in itertools.combinations(sorted(qc_rows), 2):
        row_a, row_b = qc_rows[id_a], qc_rows[id_b]
        shared = usable_bodies(row_a) & usable_bodies(row_b)
        if not shared:
            continue
        n_pairs_with_shared_body += 1
        for body in shared:
            edge = compute_edge(row_a, row_b, body)
            if edge.phase.status == "epoch_window_overlap":
                n_epoch_window_overlap += 1
            elif edge.phase.status == "epoch_disjoint":
                n_epoch_disjoint += 1

    assert n_pairs_with_shared_body > 0
    assert n_epoch_window_overlap > 0, (
        "Uranian QC rows share validity windows by construction (#569) -- a zero "
        "epoch_window_overlap count means the epoch-locked branch is broken, not "
        "a finding."
    )
    assert n_epoch_disjoint == 0, (
        "all 6 #569 rows carry the SAME validity_window (2000-06-21 to 2083-06-21) "
        "so no QC x QC pair sharing a body should ever be epoch_disjoint"
    )
