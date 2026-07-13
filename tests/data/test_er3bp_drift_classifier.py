"""Positive control for the geocentric-ER3BP bounded-vs-divergent drift
classifier (#583).

Per the task's own mandate (mirroring
[[feedback_verify_gauntlet_with_positive_control]]): a NEW long-span gate
must be validated against known-good AND known-bad cases before it judges
anything new. The known-good cases here are #581 stage 2's own 11
recognizably-reproduced Gurfil-Kasdin (2002) families, propagated from their
PUBLISHED Table 3 ICs (sourced, not our own computed values) -- see
``scripts/run_581_gurfil_reproduction.py::TABLE34``.

HONEST RESULT (not force-fitted to "11/11 bounded"): propagating all 11
families to the classifier's full N_REVS_DEFAULT=50-year horizon -- ten times
further than Gurfil-Kasdin's own longest tested span (their Table 4
"practical stability" check only extends to 5 years) -- shows:

  * 8 families stay solidly bounded through the full 50 years: A, C, D, E,
    G, H, K, M (growth_ratio and trend_fraction both near 1.0/0.0 -- a flat,
    stationary r-band, not a borderline pass).
  * Family L is published by Gurfil-Kasdin THEMSELVES as a "3D DEO"
    (Departure/Escape Orbit) type -- i.e. the paper's OWN taxonomy already
    calls it a non-practically-stable, eventually-escaping family (#581's
    own ``characterize()`` already found ``terminated_5yr=True`` for L in
    the stage-2 run). Classifying L divergent is therefore a CORRECT
    positive-control result for a divergent case, not a classifier failure.
  * Families B and J (DRO / 3D DRO -- types the paper DOES call "practically
    stable") classify DIVERGENT at the 50-year horizon: both show a long,
    clean, oscillating-but-flat r-band for 25-38 years, THEN an
    unmistakable, tolerance-independent (checked at rtol 1e-9/1e-11/1e-12,
    escape fires at the identical window every time) run to the escape
    boundary over the final few windows. This is a genuine finding, not a
    classifier bug: Gurfil-Kasdin's own "practically stable" claim is only
    ever tested to 5 years (Table 4), and both B and J are still bounded at
    5 years (matches #581 stage 2's own ``rmax_km_5yr`` figures) -- the
    paper never asserted anything about year 25-50. A 50-year classifier
    probing ten times past the paper's own tested horizon finding real,
    clean, late-onset chaotic escape in 2 of 11 families is the classifier
    doing exactly its job (catching what a short-window fitness/PS check
    cannot see), not evidence it is broken. See
    ``docs/notes/2026-07-13-583-corpus-anchors-and-drift-classifier.md``
    for the full per-window trace and tolerance-independence check.

This module therefore tests the classifier against a DIAGNOSED, fully
accounted-for expectation for all 11 reproduced families (8 bounded, 3
divergent with a documented reason each) rather than silently dropping the
inconvenient two or fudging a threshold to manufacture "11/11".
"""

from __future__ import annotations

import numpy as np
import pytest

from cyclerfinder.core.er3bp_geocentric import table_interleaved_to_state
from cyclerfinder.data.validation.er3bp_drift_classifier import (
    N_REVS_DEFAULT,
    classify_bounded_drift,
    spot_check_theta0_robustness,
)
from scripts.run_581_gurfil_reproduction import TABLE34

# #581 stage 2 "recognizably reproduced" families (results table,
# docs/notes/2026-07-12-581-niching-ga-stage2-positive-control.md #4.1).
REPRODUCED_FAMILIES = ("A", "B", "C", "D", "E", "G", "H", "J", "K", "L", "M")

# Solidly bounded through the full 50-year horizon (flat r-band both by
# growth-ratio and trend-fraction, nowhere near either threshold).
BOUNDED_AT_50YR = ("A", "C", "D", "E", "G", "H", "K", "M")

# Divergent at 50 years for a DOCUMENTED reason (see module docstring):
# L is the paper's own DEO/escape type; B and J are genuine late-onset
# chaotic escapes beyond the paper's own 5-year-tested horizon.
DIVERGENT_AT_50YR = ("L", "B", "J")


def _family_state(fam: str) -> tuple[np.ndarray, float]:
    icv, theta0, _ftype, _rmin_km, _rmax_km = TABLE34[fam]
    return table_interleaved_to_state(np.array(icv)), theta0


@pytest.mark.parametrize("fam", BOUNDED_AT_50YR)
def test_known_bounded_family_classifies_bounded(fam: str) -> None:
    state, theta0 = _family_state(fam)
    v = classify_bounded_drift(state, theta0, n_revs=N_REVS_DEFAULT, candidate_id=fam)
    assert v.bounded, (
        f"family {fam} (published Gurfil-Kasdin IC) expected BOUNDED at "
        f"{N_REVS_DEFAULT}yr but classified divergent: {v.notes}"
    )
    assert not v.terminated_early
    assert v.n_windows_complete == N_REVS_DEFAULT
    # A genuinely stationary r-band, not a borderline pass.
    assert v.growth_ratio < 1.5
    assert v.trend_fraction < 0.15


@pytest.mark.parametrize("fam", DIVERGENT_AT_50YR)
def test_known_divergent_family_classifies_divergent(fam: str) -> None:
    state, theta0 = _family_state(fam)
    v = classify_bounded_drift(state, theta0, n_revs=N_REVS_DEFAULT, candidate_id=fam)
    assert not v.bounded, (
        f"family {fam} expected DIVERGENT at {N_REVS_DEFAULT}yr (see module "
        f"docstring for the sourced diagnosis) but classified bounded"
    )
    assert v.terminated_early
    assert v.termination_reason == "escape"


def test_all_11_reproduced_families_are_accounted_for() -> None:
    """No family silently dropped: bounded + divergent sets partition all 11."""
    assert set(BOUNDED_AT_50YR) | set(DIVERGENT_AT_50YR) == set(REPRODUCED_FAMILIES)
    assert set(BOUNDED_AT_50YR).isdisjoint(DIVERGENT_AT_50YR)


def test_b_and_j_escape_is_tolerance_independent() -> None:
    """B and J's late escape is real chaotic dynamics, not integrator noise.

    Escapes at the SAME window index across three tolerance settings
    spanning three orders of magnitude (rtol/atol 1e-9, 1e-11, 1e-12) rules
    out DOP853 error accumulation as the cause.
    """
    for fam in ("B", "J"):
        state, theta0 = _family_state(fam)
        windows = set()
        for tol in (1e-9, 1e-11, 1e-12):
            v = classify_bounded_drift(
                state, theta0, n_revs=45, rtol=tol, atol=tol, candidate_id=fam
            )
            assert v.terminated_early and v.termination_reason == "escape"
            windows.add(v.n_windows_complete)
        assert len(windows) == 1, (
            f"family {fam}: escape window varied across tolerances {windows} "
            f"-- would suggest a numerical artifact rather than real chaos"
        )


def test_deliberately_perturbed_ic_classifies_divergent() -> None:
    """A synthetic off-basin perturbation of a bounded DRO escapes almost
    immediately -- the literal 'known-escaping IC' the task asks for,
    independent of any published family."""
    icv_a, theta0, _ftype, _rmin, _rmax = TABLE34["A"]
    icv = list(icv_a)
    icv[3] = -0.02  # y'0 way off the Henon DRO ridge (published -0.07776)
    state = table_interleaved_to_state(np.array(icv))
    v = classify_bounded_drift(state, theta0, n_revs=15, candidate_id="A-perturbed-escape")
    assert not v.bounded
    assert v.terminated_early
    assert v.termination_reason == "escape"
    assert v.n_windows_complete <= 2  # escapes almost immediately


def test_n_revs_below_four_rejected() -> None:
    state, theta0 = _family_state("A")
    with pytest.raises(ValueError, match="n_revs must be"):
        classify_bounded_drift(state, theta0, n_revs=3)


def test_theta0_robustness_spot_check_shape() -> None:
    """The spot-check returns one bounded/divergent verdict per phase tested."""
    state, theta0 = _family_state("C")  # deepest, most robust bounded DRO
    out = spot_check_theta0_robustness(state, theta0, n_revs=15)
    assert len(out) == 3
    # C sits deep in the DRO ridge (rmin/rmax ratio ~1.14, closest of A/B/C
    # to Earth) -- expect it to stay bounded across all three phases tested.
    assert all(out.values()), f"family C flipped bounded->divergent across phases: {out}"


def test_growth_ratio_and_trend_fraction_thresholds_are_module_constants() -> None:
    """Thresholds are written into the code, not left as a builder TODO."""
    from cyclerfinder.data.validation.er3bp_drift_classifier import (
        ESCAPE_RADIUS_NORM_DEFAULT,
        RBAND_TREND_FRACTION_THRESHOLD,
        RMAX_GROWTH_RATIO_THRESHOLD,
        SAMPLES_PER_REV,
    )

    assert RMAX_GROWTH_RATIO_THRESHOLD == 3.0
    assert RBAND_TREND_FRACTION_THRESHOLD == 0.30
    assert ESCAPE_RADIUS_NORM_DEFAULT == 0.5
    assert SAMPLES_PER_REV == 500
    assert N_REVS_DEFAULT == 50
