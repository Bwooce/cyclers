"""Axis-A code-path agreement tests — the Forge, Phase 2.

Covers :mod:`cyclerfinder.verify.agreement`:

* Aldrin agreement gate — the resonance construction and the cycler's own
  Lambert-built V_inf (two independent in-house code paths) agree within
  tolerance, AND the forward Kepler re-propagation residual is small. The
  EXPECTED side is *agreement between two independent in-house paths*
  (a consistency predicate), not an assertion against a fabricated
  constant — golden-discipline compliant. The sourced Aldrin V_inf may be
  asserted additionally because it is independently sourced.
* Kepler re-propagation teeth — a corrupted (perturbed-state) cycler FAILS
  the residual gate.
* Combiner semantics — ``agreed`` requires >= 2 available paths all
  passing; an available-but-failing path vetoes.

The fast path runs without an ephemeris (paths (b) + (c) only); the
lamberthub integration test (path (a)) lives in
``test_agreement_lamberthub.py`` and is marked slow.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.model.cycler import Cycler, Encounter, Leg, orbit_elements_au
from cyclerfinder.search.construct import build_aldrin_seed
from cyclerfinder.search.resonant_construct import construct_resonant_cycler
from cyclerfinder.verify.agreement import (
    KEPLER_REPROP_TOL_KM,
    VINF_AGREEMENT_TOL_KMS,
    AgreementReport,
    crosscheck_code_paths,
)

# Independently sourced Aldrin V_inf (Russell 2004 Table 3.4 cycler
# 1.0.1.-1; Rogers 2012 Table 1). Asserted only as a *sourced* sanity
# anchor alongside the path-agreement predicate, never as the EXPECTED
# side of an internal-vs-internal agreement.
SOURCED_VINF_EARTH_KMS = 6.5
SOURCED_VINF_MARS_KMS = 9.7
SOURCED_VINF_TOL_KMS = 0.5


@pytest.fixture(scope="module")
def aldrin_seed() -> Cycler:
    """The canonical Aldrin E->M slice built by the Lambert constructor."""
    eph = Ephemeris(model="circular")
    return build_aldrin_seed(eph)


# ---------------------------------------------------------------------------
# Aldrin agreement gate (fast — no ephemeris, paths (b) + (c))
# ---------------------------------------------------------------------------


def test_aldrin_construction_paths_agree(aldrin_seed: Cycler) -> None:
    """GATE — resonance construction vs Lambert-built V_inf agree within tol.

    Two independent in-house derivations of the Aldrin V_inf:

    * the Lambert boundary-value solve carried on ``aldrin_seed``;
    * the closed-form vis-viva crossing construction
      (:func:`construct_resonant_cycler`) from the seed's own ``(a, e)``.

    They must agree to within :data:`VINF_AGREEMENT_TOL_KMS`. This is the
    Axis-A construction-vs-construction predicate — golden-discipline
    compliant (no fabricated constant on the EXPECTED side).
    """
    report = crosscheck_code_paths(aldrin_seed)
    path_b = report.construction_optimiser
    assert path_b.available is True
    assert path_b.passed is True, (
        f"construction paths disagree: max_diff={path_b.construction_max_diff_kms} km/s "
        f"resonant={path_b.resonant_vinf_kms} cycler={path_b.cycler_vinf_kms}"
    )
    assert path_b.construction_max_diff_kms < VINF_AGREEMENT_TOL_KMS
    # The two paths are genuinely independent yet match to ~1e-13 km/s.
    assert path_b.construction_max_diff_kms < 1.0e-6


def test_aldrin_kepler_reprop_small(aldrin_seed: Cycler) -> None:
    """GATE — forward Kepler re-propagation residual is sub-km on the seed."""
    report = crosscheck_code_paths(aldrin_seed)
    path_c = report.kepler_reprop
    assert path_c.available is True
    assert path_c.passed is True, (
        f"Kepler re-propagation residual too large: {path_c.max_residual_km} km"
    )
    assert path_c.max_residual_km < KEPLER_REPROP_TOL_KM
    # Empirically ~2e-4 km on this seed.
    assert path_c.max_residual_km < 1.0e-2


def test_aldrin_two_paths_agree_overall(aldrin_seed: Cycler) -> None:
    """GATE — with two available paths both passing, ``agreed`` is True."""
    report = crosscheck_code_paths(aldrin_seed)
    assert isinstance(report, AgreementReport)
    # No ephemeris ⇒ lamberthub path unavailable; (b) and (c) available.
    assert report.lamberthub.available is False
    assert report.n_paths_available == 2
    assert report.n_paths_passed == 2
    assert report.agreed is True


def test_aldrin_sourced_vinf_anchor(aldrin_seed: Cycler) -> None:
    """The cycler's V_inf matches the *sourced* Aldrin anchors (6.5 / 9.7).

    Asserted additionally to the path-agreement predicate because these
    values are independently sourced (Russell 2004 / Rogers 2012), not
    produced by our own code.
    """
    vinf = {
        enc.body: max(
            float(np.linalg.norm(enc.vinf_in)),
            float(np.linalg.norm(enc.vinf_out)),
        )
        for enc in aldrin_seed.encounters
    }
    assert abs(vinf["E"] - SOURCED_VINF_EARTH_KMS) < SOURCED_VINF_TOL_KMS
    assert abs(vinf["M"] - SOURCED_VINF_MARS_KMS) < SOURCED_VINF_TOL_KMS


# ---------------------------------------------------------------------------
# Kepler re-propagation teeth — a corrupted cycler must FAIL the residual gate
# ---------------------------------------------------------------------------


def _perturb_departure_velocity(cycler: Cycler, factor: float) -> Cycler:
    """Return a copy of ``cycler`` whose first encounter's ``vinf_out`` (and
    matching first leg ``v_depart``) is scaled by ``factor``.

    Scaling ``vinf_out`` corrupts the spacecraft departure state the
    Kepler re-propagation reconstructs, so the propagated arrival no
    longer matches the claimed encounter position.
    """
    enc0 = cycler.encounters[0]
    bad_vinf = np.asarray(enc0.vinf_out, dtype=np.float64) * factor
    new_enc0 = dataclasses.replace(enc0, vinf_out=bad_vinf, vinf_in=bad_vinf)
    new_encounters = [new_enc0, *cycler.encounters[1:]]
    leg0 = cycler.legs[0]
    bad_v_depart = np.asarray(enc0.v_planet, dtype=np.float64) + bad_vinf
    new_leg0 = dataclasses.replace(leg0, v_depart=bad_v_depart)
    new_legs = [new_leg0, *cycler.legs[1:]]
    return dataclasses.replace(cycler, encounters=new_encounters, legs=new_legs)


def test_kepler_reprop_teeth_corrupted_fails(aldrin_seed: Cycler) -> None:
    """A 1 %-perturbed departure state FAILS the Kepler residual gate.

    Teeth test: confirms the residual gate actually rejects a corrupted
    cycler rather than passing everything. The clean seed passes
    (:func:`test_aldrin_kepler_reprop_small`); the perturbed copy must
    diverge by orders of magnitude above :data:`KEPLER_REPROP_TOL_KM`.
    """
    corrupted = _perturb_departure_velocity(aldrin_seed, factor=1.01)
    report = crosscheck_code_paths(corrupted)
    path_c = report.kepler_reprop
    assert path_c.available is True
    assert path_c.passed is False, (
        f"corrupted cycler unexpectedly passed: residual={path_c.max_residual_km} km"
    )
    assert path_c.max_residual_km > KEPLER_REPROP_TOL_KM
    # Empirically diverges by millions of km for a 1 % perturbation.
    assert path_c.max_residual_km > 1.0e5


def test_corrupted_cycler_not_agreed(aldrin_seed: Cycler) -> None:
    """A corrupted cycler's overall verdict is not ``agreed`` (a failing
    available path vetoes)."""
    corrupted = _perturb_departure_velocity(aldrin_seed, factor=1.01)
    report = crosscheck_code_paths(corrupted)
    assert report.kepler_reprop.passed is False
    assert report.agreed is False


# ---------------------------------------------------------------------------
# Combiner semantics
# ---------------------------------------------------------------------------


def test_optimiser_witness_folds_in(aldrin_seed: Cycler) -> None:
    """A supplied optimiser V_inf that matches is folded into path (b)."""
    # Use the resonance construction's own output as a stand-in "optimiser"
    # witness: it is the same physical family, so it agrees. (The real
    # idealised optimiser does not converge on the E-M family — see the
    # documented scope boundary; this exercises the fold-in wiring.)
    a_au, e = orbit_elements_au(
        aldrin_seed.encounters[0].r, aldrin_seed.legs[0].v_depart, mu=MU_SUN_KM3_S2
    )
    rc = construct_resonant_cycler(a_au, e, bodies=("E", "M"))
    report = crosscheck_code_paths(aldrin_seed, optimiser_vinf_kms=dict(rc.vinf_kms))
    path_b = report.construction_optimiser
    assert path_b.optimiser_available is True
    assert path_b.optimiser_max_diff_kms is not None
    assert path_b.passed is True


def test_optimiser_witness_disagreement_vetoes(aldrin_seed: Cycler) -> None:
    """A supplied optimiser V_inf that disagrees fails path (b)."""
    bad = {"E": 99.0, "M": 99.0}
    report = crosscheck_code_paths(aldrin_seed, optimiser_vinf_kms=bad)
    path_b = report.construction_optimiser
    assert path_b.optimiser_available is True
    assert path_b.passed is False
    assert report.agreed is False


def test_single_path_not_agreed() -> None:
    """A cycler exposing only one available path is never ``agreed``.

    A 2-encounter cycler whose orbit is hyperbolic on its single leg makes
    the resonance construction unavailable (path b), leaving only the
    Kepler re-propagation (path c) — one path is insufficient for the
    >= 2-path predicate.
    """
    # Build a degenerate single-leg "cycler" with a hyperbolic departure so
    # the elliptic crossing construction is undefined (path b unavailable),
    # while the Kepler re-propagation still runs (path c available).
    mu = MU_SUN_KM3_S2
    r0 = np.array([1.496e8, 0.0, 0.0], dtype=np.float64)
    v_circ = float(np.sqrt(mu / np.linalg.norm(r0)))
    v_hyper = np.array([0.0, 2.0 * v_circ, 0.0], dtype=np.float64)  # > escape ⇒ hyperbolic
    enc0 = Encounter(
        body="E",
        t=0.0,
        r=r0,
        v_planet=np.array([0.0, v_circ, 0.0], dtype=np.float64),
        vinf_in=v_hyper - np.array([0.0, v_circ, 0.0], dtype=np.float64),
        vinf_out=v_hyper - np.array([0.0, v_circ, 0.0], dtype=np.float64),
    )
    # Propagate the hyperbolic state forward to define the arrival encounter
    # so path (c) closes exactly (residual ~ 0) and is the lone passing path.
    from cyclerfinder.core.kepler import propagate

    tof = 50.0 * 86400.0
    r1, _ = propagate(r0, v_hyper, tof, mu=mu)
    enc1 = Encounter(
        body="M",
        t=tof,
        r=r1,
        v_planet=np.zeros(3, dtype=np.float64),
        vinf_in=np.zeros(3, dtype=np.float64),
        vinf_out=np.zeros(3, dtype=np.float64),
    )
    leg = Leg(
        from_body="E",
        to_body="M",
        t_depart=0.0,
        t_arrive=tof,
        v_depart=v_hyper,
        v_arrive=np.zeros(3, dtype=np.float64),
    )
    cyc = Cycler(bodies=["E", "M"], period=tof, encounters=[enc0, enc1], legs=[leg])

    report = crosscheck_code_paths(cyc)
    assert report.construction_optimiser.available is False
    assert report.kepler_reprop.available is True
    assert report.kepler_reprop.passed is True
    assert report.n_paths_available == 1
    assert report.agreed is False
