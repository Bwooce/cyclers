"""Tests for the EGGIE resonant-conic initial-guess generator (#480 follow-up 1).

The make-or-break gate is the positive control (``feedback_verify_gauntlet_with_
positive_control``): does the single resonant conic put all three Galilean V∞ on
the sourced Table-4 targets? It does (``test_eggie_pure_conic_vinf_in_band`` — the
diagnosis spike-4 crux, reproduced deterministically). The full patched-conic tour
then holds V∞ in band but reaches only ΣΔV ≈ 1.8e2 m/s (not the paper's near-
ballistic 0.70 m/s) — the right basin, not the ballistic floor; that gap is the
Stage-3 corrector's job and is asserted/documented honestly here, NOT loosened.

All EXPECTED V∞/ToF values are SOURCED to Hernandez-Jones-Jesick 2017 (AAS 17-608)
Table 4 via docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md
(``feedback_golden_tests_sourced_only`` — never a value our own code computed).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from cyclerfinder.search.resonant_conic import (
    EGGIE_SEQUENCE,
    EGGIE_TOFS_TABLE4_DAYS,
    EGGIE_VINF_TARGET_KMS,
    MU_JUPITER_KM3_S2,
    conic_radius,
    conic_state,
    conic_vinf,
    crossing_true_anomalies,
    ecc_bounds,
    eggie_initial_guess,
    eggie_resonant_sma,
    ideal_moon_smas,
    ideal_t_syn,
    refine_eggie,
    resonant_period,
    resonant_sma,
    time_since_periapsis,
)

MU = MU_JUPITER_KM3_S2


# --- Ideal-model constants (sourced to paper p.3 / Table 1) -----------------------


def test_ideal_moon_smas_match_paper_p3() -> None:
    """Ideal a_Eur/a_Gan from the 5.2°-per-synodic factors (paper p.3)."""
    smas = ideal_moon_smas()
    # a_Io is the real Io sma (registry); Europa/Ganymede ~ real values (cross-check).
    assert smas["Io"] == pytest.approx(421800.0)
    assert smas["Europa"] == pytest.approx(667964.0, abs=200.0)  # cf. real 671100
    assert smas["Ganymede"] == pytest.approx(1055289.0, abs=300.0)  # cf. real 1070400


def test_resonant_sma_eggie_45() -> None:
    """EGGIE 4:5 resonant sma ≈ 9.094e5 km (plan/diagnosis sourced value)."""
    a = resonant_sma(4, 5, ideal_t_syn())
    assert a == pytest.approx(909420.0, abs=50.0)
    assert eggie_resonant_sma() == pytest.approx(a)
    # T_sc = 4/5 · T_syn (Eq. 1).
    assert resonant_period(a) == pytest.approx(0.8 * ideal_t_syn(), rel=1e-9)


def test_ecc_bounds_match_sourced() -> None:
    """e ∈ [0.536, 0.921] (diagnosis spike 4): r_p ≤ a_Io binds e_min, r_p ≥ R_Jup e_max."""
    a = eggie_resonant_sma()
    smas = ideal_moon_smas()
    e_min, e_max = ecc_bounds(a, a_io=smas["Io"], a_gan=smas["Ganymede"])
    assert e_min == pytest.approx(0.5362, abs=1e-3)
    assert e_max == pytest.approx(0.9214, abs=1e-3)


# --- Conic geometry helpers -------------------------------------------------------


def test_conic_radius_periapsis_apoapsis() -> None:
    a, e = 9.0e5, 0.62
    assert conic_radius(a, e, 0.0) == pytest.approx(a * (1.0 - e))
    assert conic_radius(a, e, math.pi) == pytest.approx(a * (1.0 + e))


def test_conic_state_radius_and_visviva() -> None:
    """conic_state position magnitude == r(nu); speed == vis-viva √(mu(2/r - 1/a))."""
    a, e, omega = 9.0e5, 0.62, 0.7
    for nu in (0.3, 1.5, 3.0, 4.5):
        r_vec, v_vec = conic_state(a, e, omega, nu, MU)
        r = conic_radius(a, e, nu)
        assert float(np.linalg.norm(r_vec)) == pytest.approx(r, rel=1e-12)
        vis_viva = math.sqrt(MU * (2.0 / r - 1.0 / a))
        assert float(np.linalg.norm(v_vec)) == pytest.approx(vis_viva, rel=1e-12)


def test_time_since_periapsis_apoapsis_is_half_period() -> None:
    a, e = 9.0e5, 0.62
    t = time_since_periapsis(a, e, math.pi, MU)
    assert t == pytest.approx(0.5 * resonant_period(a, MU), rel=1e-9)


def test_crossing_true_anomalies_hit_moon_radius() -> None:
    a, e = eggie_resonant_sma(), 0.62
    smas = ideal_moon_smas()
    nus = crossing_true_anomalies(a, e, smas["Ganymede"])
    assert nus is not None
    nu_out, nu_in = nus
    assert 0.0 < nu_out < math.pi < nu_in < 2.0 * math.pi
    assert conic_radius(a, e, nu_out) == pytest.approx(smas["Ganymede"], rel=1e-9)
    assert conic_radius(a, e, nu_in) == pytest.approx(smas["Ganymede"], rel=1e-9)
    # A circle below perijove is never reached.
    assert crossing_true_anomalies(a, e, 1.0e4) is None


# --- POSITIVE CONTROL #1: pure-conic V∞ on the Table-4 targets (spike 4) ----------


def test_eggie_pure_conic_vinf_in_band() -> None:
    """The make-or-break gate: ONE resonant conic puts all three V∞ on Table 4.

    Deterministic, no optimisation. Achieved at e=0.620: Europa 9.08, Ganymede 6.78,
    Io 8.35 km/s vs SOURCED Table-4 9.12/7.07/8.38 (Ganymede worst, 0.29 km/s low —
    inside the ~0.5 km/s ideal-vs-real-T_syn tolerance noted in the diagnosis).
    """
    e = 0.620
    for moon, target in EGGIE_VINF_TARGET_KMS.items():
        v = conic_vinf(e, moon)
        assert v == pytest.approx(target, abs=0.5), f"{moon} V∞ {v:.3f} vs Table-4 {target}"


def test_conic_vinf_raises_when_unreachable() -> None:
    # At low e the conic apojove (a(1+e)) falls short of Ganymede's orbit.
    with pytest.raises(ValueError):
        conic_vinf(0.10, "Ganymede")


# --- The conic initial guess ------------------------------------------------------


def test_eggie_initial_guess_structure_and_tofs() -> None:
    """The conic seed has the EGGIE topology, nodes on moon orbits, ToFs ~ Table 4."""
    g = eggie_initial_guess(0.620)
    assert g.sequence == EGGIE_SEQUENCE
    assert len(g.epochs_days) == 5
    assert g.epochs_days[0] == 0.0
    # Each node sits on its moon's circular orbit radius.
    smas = ideal_moon_smas()
    for k, moon in enumerate(g.sequence):
        assert float(np.linalg.norm(g.node_positions[k])) == pytest.approx(smas[moon], rel=1e-6)
    # ToFs come from the conic crossings and land near the SOURCED Table-4 values
    # (the conic G-G crossing carries the known +0.4 d resonance-closure artifact).
    for got, table4 in zip(g.tofs_days, EGGIE_TOFS_TABLE4_DAYS, strict=True):
        assert got == pytest.approx(table4, abs=0.6)


# --- POSITIVE CONTROL #2: refined tour holds V∞ in band; honest ΔV verdict --------


def test_eggie_refine_holds_vinf_in_band_and_beats_offbasin() -> None:
    """The tight refine stays in the resonant basin: V∞ in band, ΣΔV ≪ off-basin.

    SOURCED comparison: the diagnosis's per-leg-Lambert / free-ToF searches landed
    OFF basin — ΣΔV 385-1139 m/s WITH WRONG V∞ (Ganymede collapsing to ~2-4 km/s;
    docs/notes/2026-06-29-480-eggie-ideal-positive-control-diagnosis.md spikes 1-3).
    The conic-seeded refine instead holds all three V∞ on the Table-4 targets AND
    cuts ΣΔV below that off-basin floor.

    HONEST PARTIAL verdict: it is NOT near-ballistic (paper Table-4 ΔV = 0.70 m/s).
    Achieved ≈ 1.8e2 m/s interior ΔV + ≈ 4.7e2 m/s periodicity wrap — the patched-
    conic Lambert legs reach the right basin, not the ballistic floor. Closing to
    near-ballistic is the Stage-3 multiple-shooting corrector's job. Thresholds are
    set from what is ACTUALLY achieved and are NOT loosened to force a pass.
    """
    r = refine_eggie(n_restarts=3)  # deterministic; ~15 s
    # (1) all three V∞ on the SOURCED Table-4 targets (abs 0.5 km/s).
    assert r.vinf_in_band
    targets = [EGGIE_VINF_TARGET_KMS[m] for m in EGGIE_SEQUENCE]
    for k, (v, t) in enumerate(zip(r.vinf_kms, targets, strict=True)):
        assert v == pytest.approx(t, abs=0.5), f"node {k} V∞ {v:.3f} vs {t}"
    # (2) the two Ganymede flybys are at equal V∞ (ballistic-cycler property, paper §3).
    assert r.vinf_kms[1] == pytest.approx(r.vinf_kms[2], abs=0.5)
    # (3) ΣΔV is well below the off-basin per-leg-Lambert floor (385 m/s, diagnosis).
    assert r.sum_interior_dv_ms < 385.0
    # (4) but NOT near-ballistic — documents the PARTIAL verdict (achieved ~184 m/s).
    assert r.sum_interior_dv_ms > 50.0
