"""M7 leg-targeting solver (#423) — golden + Mars-perturbed smoke.

GOLDEN GATE: in Sun-only mode (``bodies=()``) the n-body position-targeting Newton
:func:`~cyclerfinder.nbody.maintenance_shoot.target_leg` must reproduce the INDEPENDENT
two-body :func:`cyclerfinder.core.lambert.lambert` departure velocity for the same
boundary-value leg. Lambert is the sourced/independent cross-check (a different solver,
analytic universal-variable), so agreement validates the targeting kernel without
circularity.
"""

from __future__ import annotations

import numpy as np
import pytest

rebound = pytest.importorskip("rebound")

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS  # noqa: E402
from cyclerfinder.core.ephemeris import Ephemeris  # noqa: E402
from cyclerfinder.core.lambert import lambert  # noqa: E402
from cyclerfinder.nbody.maintenance_shoot import target_leg  # noqa: E402
from cyclerfinder.nbody.propagator import RestrictedNBody  # noqa: E402

_DAY_S = 86400.0


def test_target_leg_sun_only_reproduces_lambert() -> None:
    """Sun-only targeting Newton recovers the two-body Lambert departure velocity.

    Real Earth->Mars leg over DE440: take r0 = Earth(t0), r_target = Mars(t1) for a
    ~210 d transfer, solve the two-body Lambert for the departure velocity, then run
    target_leg in Sun-only mode seeded with a deliberately PERTURBED guess. It must
    converge back to the Lambert velocity (km/s) and land sub-km on the target."""
    ephem = Ephemeris("astropy")
    prop = RestrictedNBody("rebound")

    t0 = 27.0 * 365.25 * _DAY_S  # ~2027, an arbitrary in-ephemeris epoch
    tof = 210.0 * _DAY_S
    t1 = t0 + tof
    r0, _ = ephem.state("E", t0)
    r_target, _ = ephem.state("M", t1)
    r0 = np.asarray(r0, dtype=np.float64)
    r_target = np.asarray(r_target, dtype=np.float64)

    sols = lambert(r0, r_target, tof, mu=MU_SUN_KM3_S2)
    v_lambert = np.asarray(sols[0].v1, dtype=np.float64)  # single-rev departure velocity

    # Seed with a perturbed guess (off by ~0.3 km/s in each component) so the Newton
    # iteration has real work to do — it must not just echo the seed.
    v_guess = v_lambert + np.array([0.3, -0.3, 0.2])
    res = target_leg(
        prop,
        r0,
        t0,
        t1,
        r_target,
        v_guess,
        bodies=(),  # Sun-only two-body — the Lambert regime
        ephem=ephem,
        tol_km=1.0,
    )

    assert res.converged, f"Sun-only targeting did not converge (miss {res.miss_km:.1f} km)"
    assert res.miss_km < 1.0, f"arrival miss {res.miss_km:.3f} km exceeds 1 km"
    # Recovers the independent Lambert departure velocity (mm/s-level agreement).
    assert res.v_dep_km_s == pytest.approx(v_lambert, abs=1e-5), (
        f"target_leg v_dep {res.v_dep_km_s} != lambert {v_lambert}"
    )


@pytest.mark.slow
def test_target_leg_mars_perturbed_converges_in_band() -> None:
    """Mars-perturbed targeting converges and lands far inside the Mars 3-SOI band.

    Same leg, now with Mars as a continuous perturber. The Mars-perturbed departure
    velocity differs from the two-body Lambert value (the perturbation is real work),
    but the solve must still drive the arrival miss sub-km — i.e. the STM Jacobian is
    correct under the perturber (the REBOUND-variation gravity-gradient gotcha is
    handled in the propagator)."""
    ephem = Ephemeris("astropy")
    prop = RestrictedNBody("rebound")

    t0 = 27.0 * 365.25 * _DAY_S
    tof = 210.0 * _DAY_S
    t1 = t0 + tof
    r0, _ = ephem.state("E", t0)
    r_target, _ = ephem.state("M", t1)
    r0 = np.asarray(r0, dtype=np.float64)
    r_target = np.asarray(r_target, dtype=np.float64)

    v_lambert = np.asarray(lambert(r0, r_target, tof, mu=MU_SUN_KM3_S2)[0].v1, dtype=np.float64)

    res = target_leg(prop, r0, t0, t1, r_target, v_lambert, bodies=("M",), ephem=ephem, tol_km=1.0)

    mars_soi_au = PLANETS["M"].sma_au * (PLANETS["M"].mu_km3_s2 / MU_SUN_KM3_S2) ** 0.4
    assert res.converged, f"Mars-perturbed targeting did not converge (miss {res.miss_km:.1f} km)"
    assert res.miss_km < 1.0, f"arrival miss {res.miss_km:.3f} km exceeds 1 km"
    assert res.miss_km / AU_KM < 3.0 * mars_soi_au  # trivially true; documents the band
