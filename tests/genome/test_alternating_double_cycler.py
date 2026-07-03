"""Tests for the alternating-double-cycler construction operator (#526).

Positive controls only (sourced-only rule): every EXPECTED value traces to
Liang, Yang, Li, Bai & Qin 2024 (JGCD, DOI 10.2514/1.G008387) printed
Table 1 / Eqs. 1-2, or to this project's own already-committed sourced
constants (CGE_SEQUENCE, the JPL-SSD satellite registry) -- never to a
value this module itself computes.
"""

from __future__ import annotations

import math

import pytest

from cyclerfinder.genome.alternating_double_cycler import (
    analyze_near_resonance,
    build_alternating_double_cycler_seed,
    liang_cge_alternating_seed,
    liang_cge_near_resonance,
    mean_motion_rad_day,
    synodic_period_days,
)
from cyclerfinder.search.moon_cycler_genome import (
    CGE_SEQUENCE,
    LegGene,
    jupiter_system,
    saturn_system,
)


def test_synodic_period_matches_liang_table1_values() -> None:
    # Liang et al. 2024 Sec. II.A (unlabeled eq. before Eq. 1) / Eq. 1 (p. 3):
    # S_C,G = 12.5232 d, S_G,E = 7.0509 d, from Table 1 mean motions
    # (Europa 1.7693, Ganymede 0.8782, Callisto 0.3765 rad/day, p. 4).
    s_cg = synodic_period_days(0.8782, 0.3765)
    s_ge = synodic_period_days(1.7693, 0.8782)
    assert s_cg == pytest.approx(12.5232, abs=1e-3)
    assert s_ge == pytest.approx(7.0509, abs=2e-3)


def test_synodic_period_rejects_co_orbital_bodies() -> None:
    with pytest.raises(ValueError, match="co-orbital"):
        synodic_period_days(1.0, 1.0)


def test_liang_cge_near_resonance_reproduces_eq1_eq2() -> None:
    """Sourced golden: Liang et al. 2024 Eqs. 1-2 (p. 3), from their own Table 1."""
    analysis = liang_cge_near_resonance()
    # Eq. 1: ratio ~ 1.7761 ~ 7/4.
    assert analysis.ratio == pytest.approx(1.7761, abs=2e-4)
    assert (analysis.p, analysis.q) == (7, 4)
    # Eq. 2: mismatch = 4*S_C,G - 7*S_G,E = 50.0928 - 49.3563 = 0.7365 d. Our
    # S_C,G/S_G,E (computed from the same printed Table 1 mean motions via
    # the exact 2*pi/dn formula) differ from the paper's printed 4-decimal
    # values in the 3rd decimal, so the derived mismatch/quasi-period carry
    # a matching few-1e-3 tolerance -- NOT rounded to match, genuinely close.
    assert analysis.mismatch_days == pytest.approx(0.7365, abs=3e-3)
    assert analysis.quasi_period_days == pytest.approx(50.0928, abs=3e-3)


def test_analyze_near_resonance_rejects_degenerate_denominator() -> None:
    with pytest.raises(ValueError, match="max_denominator"):
        analyze_near_resonance("A", "B", "C", {"A": 1.0, "B": 0.5, "C": 0.25}, max_denominator=0)


def test_registry_independent_cross_check_recovers_same_near_resonance() -> None:
    """Independent (non-golden) check: JPL-SSD registry mean motions, NOT Liang's
    printed Table 1, still recover the same qualitative 7:4 near-resonance within
    ~2% of Liang's own mismatch -- two independent data sources agreeing."""
    system = jupiter_system()
    motions = {m: mean_motion_rad_day(system, m) for m in ("Callisto", "Ganymede", "Europa")}
    analysis = analyze_near_resonance("Callisto", "Ganymede", "Europa", motions)
    assert (analysis.p, analysis.q) == (7, 4)
    liang = liang_cge_near_resonance()
    assert analysis.mismatch_days == pytest.approx(liang.mismatch_days, rel=0.05)
    assert analysis.quasi_period_days == pytest.approx(liang.quasi_period_days, rel=0.02)


def test_mismatch_fraction_is_scale_free() -> None:
    analysis = liang_cge_near_resonance()
    assert analysis.mismatch_fraction == pytest.approx(
        abs(analysis.mismatch_days) / analysis.quasi_period_days
    )
    assert 0.0 < analysis.mismatch_fraction < 1.0


def test_build_alternating_seed_reproduces_published_cgcec_sequence() -> None:
    """Structural positive control: stitching Liang's own two halves must
    reproduce their published CGCEC sequence exactly."""
    seed = liang_cge_alternating_seed()
    assert seed.sequence == CGE_SEQUENCE
    assert seed.is_valid()
    assert len(seed.legs) == len(CGE_SEQUENCE) - 1


def test_build_alternating_seed_multi_cycle_drops_junction_duplicates() -> None:
    system = jupiter_system()
    seed = build_alternating_double_cycler_seed(
        system,
        half_cycle_1=("Callisto", "Ganymede", "Callisto"),
        half_cycle_2=("Callisto", "Europa", "Callisto"),
        n_cycles=2,
    )
    # One cycle is 5 bodies (4 legs); two cycles share the junction hub, so
    # 5 + 4 = 9 bodies / 8 legs, never 10 bodies with a duplicated hub pair.
    assert seed.sequence == CGE_SEQUENCE + CGE_SEQUENCE[1:]
    assert len(seed.sequence) == 9
    assert len(seed.legs) == 8
    assert seed.is_valid()


def test_build_alternating_seed_rejects_mismatched_hub() -> None:
    system = jupiter_system()
    with pytest.raises(ValueError, match="must share the same hub"):
        build_alternating_double_cycler_seed(
            system,
            half_cycle_1=("Callisto", "Ganymede", "Callisto"),
            half_cycle_2=("Ganymede", "Europa", "Ganymede"),
        )


def test_build_alternating_seed_rejects_non_returning_half_cycle() -> None:
    system = jupiter_system()
    with pytest.raises(ValueError, match="start and end at the shared hub"):
        build_alternating_double_cycler_seed(
            system,
            half_cycle_1=("Callisto", "Ganymede"),
            half_cycle_2=("Callisto", "Europa", "Callisto"),
        )


def test_build_alternating_seed_rejects_wrong_leg_count() -> None:
    system = jupiter_system()
    with pytest.raises(ValueError, match="expected"):
        build_alternating_double_cycler_seed(
            system,
            half_cycle_1=("Callisto", "Ganymede", "Callisto"),
            half_cycle_2=("Callisto", "Europa", "Callisto"),
            legs=(LegGene(p=1, q=1, n_rev=1),),
        )


def test_build_alternating_seed_rejects_n_cycles_below_one() -> None:
    system = jupiter_system()
    with pytest.raises(ValueError, match="n_cycles"):
        build_alternating_double_cycler_seed(
            system,
            half_cycle_1=("Callisto", "Ganymede", "Callisto"),
            half_cycle_2=("Callisto", "Europa", "Callisto"),
            n_cycles=0,
        )


def test_reusability_on_saturn_chain_no_novelty_claim() -> None:
    """Body-agnostic reuse check ONLY: the function must run on a second,
    independent moon system (Saturn) and return a well-formed result. This
    asserts NOTHING about whether Enceladus-Dione-Rhea are actually
    near-resonant -- no novelty / discovery claim is made here."""
    system = saturn_system()
    motions = {m: mean_motion_rad_day(system, m) for m in ("Enceladus", "Dione", "Rhea")}
    analysis = analyze_near_resonance("Enceladus", "Dione", "Rhea", motions)
    assert math.isfinite(analysis.mismatch_days)
    assert analysis.quasi_period_days > 0.0
    assert analysis.p >= 1 and analysis.q >= 1
