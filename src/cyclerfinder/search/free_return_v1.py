"""§14 V1 mechanics for the #137 free-return (radial-crossing) genome.

The free-return corrector (:mod:`cyclerfinder.search.free_return`) closes a
single heliocentric ellipse ``(a, e)`` whose radial crossings of Earth and Mars
reproduce a sourced cycler's DERIVED V_inf (the #137 breakthrough, recorded in
``docs/notes/2026-06-07-russell12-freereturn-results.md``). That closure gates on
V_inf-continuity + phase to the corrector's 0.1 km/s floor; it is NOT, by itself,
the spec §14 V1 gate.

This module applies the **literal §14 V1 mechanics** to the closed geometry,
like-for-like on the CIRCULAR ephemeris (the same fidelity the sourced rows are
defined in — a circular-coplanar reproduction of a circular-coplanar source):

* **Path (a) — lamberthub re-solve.** Every leg of the reconstructed E->M->E
  free-return arc is re-solved with ``lamberthub`` izzo2015 + gooding1990 and
  compared to the in-house Lambert; agreement must be < ``V1_TOLERANCE_MPS``.
* **Path (c) — Kepler forward re-propagation.** Each leg's departure state is
  propagated forward with the universal-variable two-body propagator (NOT the
  Lambert that built it) and the arrival position compared to the planet's
  claimed position; residual must be < ``KEPLER_REPROP_TOL_KM``.

Both are reused verbatim from :func:`cyclerfinder.verify.agreement
.crosscheck_code_paths` (paths a and c) — never reimplemented here.

The §14 V1 text is exactly paths (a) and (c). The Axis-A combiner's path (b)
(resonance-construction-vs-cycler V_inf) is an additional Forge witness, not part
of §14 V1; the Aldrin INBOUND V1 evidence note already records a V1 pass resting
on the two §14 paths when (b) is unavailable. We therefore gate on (a) AND (c).

The genome-specific gate (the honesty teeth)
--------------------------------------------
A single free-return ellipse only forms a CLOSED, V_inf-continuous E->M->E
trajectory when its descending Earth crossing coincides with where Earth actually
is on the circular orbit. For the genuine single-ellipse free-return rows it does
(the Mars-flyby V_inf is continuous to ~0.01-0.2 km/s and the return arrives at
Earth ballistically). For the multi-arc rows whose return needs intermediate
phasing loops, the descending Earth crossing is ~180 deg from the real Earth, so
forcing a Lambert return leg lands a ~29 km/s Mars V_inf discontinuity — a broken,
non-physical trajectory whose lamberthub/Kepler self-consistency would pass
VACUOUSLY. :data:`VINF_CONTINUITY_TOL_KMS` rejects exactly that: V1 is awarded
only when the reconstructed arc is a genuinely closed, V_inf-continuous cycler.

Pure aside from the verify import; depends on core + the free-return corrector.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import pi, sqrt
from typing import Final

import numpy as np

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris
from cyclerfinder.core.kepler import coe_to_rv  # noqa: F401  (kept for callers/tests)
from cyclerfinder.core.lambert import (
    LambertConvergenceError,
    LambertGeometryError,
    lambert,
)
from cyclerfinder.model.cycler import Cycler, Encounter, Leg
from cyclerfinder.search.free_return import (
    FreeReturnClosureResult,
    _true_to_mean,
    free_return_geometry,
)
from cyclerfinder.verify.agreement import crosscheck_code_paths

DAY_S = SECONDS_PER_DAY

# V_inf-continuity ceiling at the Mars flyby for the reconstructed E->M->E arc,
# km/s. A genuine single-ellipse free return is continuous to ~0.01-0.2 km/s; a
# forced multi-arc return lands ~24 km/s of discontinuity. 0.5 km/s (the campaign
# match tolerance) cleanly separates the two and is NOT tightened per-row.
VINF_CONTINUITY_TOL_KMS: Final[float] = 0.5


@dataclass(frozen=True)
class FreeReturnV1Result:
    """§14 V1 mechanics verdict for a closed free-return geometry.

    Attributes
    ----------
    built:
        ``True`` when the E->M->E arc reconstructed and both legs Lambert-solved.
    lamberthub_passed:
        Path (a): in-house Lambert vs lamberthub agreement < ``V1_TOLERANCE_MPS``.
    lamberthub_max_diff_mps:
        Worst per-leg lamberthub disagreement, m/s.
    kepler_reprop_passed:
        Path (c): forward Kepler re-propagation residual < ``KEPLER_REPROP_TOL_KM``.
    kepler_reprop_max_residual_km:
        Worst per-leg re-propagation residual, km.
    vinf_continuity_kms:
        Mars-flyby ``||V_inf_in| - |V_inf_out||`` for the reconstructed arc, km/s.
    vinf_continuous:
        ``vinf_continuity_kms <= VINF_CONTINUITY_TOL_KMS`` — the genome honesty
        gate (rejects the vacuous multi-arc forced-return pass).
    v1_passed:
        ``built and lamberthub_passed and kepler_reprop_passed and
        vinf_continuous`` — the full §14 V1 (paths a+c) on a genuinely closed,
        V_inf-continuous reconstruction.
    detail:
        Human-readable note (the failure reason when not built / not passed).
    """

    built: bool
    lamberthub_passed: bool
    lamberthub_max_diff_mps: float
    kepler_reprop_passed: bool
    kepler_reprop_max_residual_km: float
    vinf_continuity_kms: float
    vinf_continuous: bool
    v1_passed: bool
    detail: str


def build_free_return_cycler(
    sol: FreeReturnClosureResult,
    ephem: Ephemeris,
    period_sec: float,
    *,
    bodies: tuple[str, str] = ("E", "M"),
    mu: float = MU_SUN_KM3_S2,
) -> Cycler:
    """Reconstruct the closed E->M->E free-return arc as a :class:`Cycler`.

    The three encounters lie on the converged free-return ellipse ``(a, e)``:
    Earth at the ascending crossing ``nu_E`` (epoch ``sol.t0_sec``), Mars at the
    ascending crossing ``nu_M``, and Earth again at the DESCENDING crossing
    ``2*pi - nu_E``. Encounter epochs are the Keplerian crossing times on the
    ellipse; each leg is an in-house Lambert solution between the ACTUAL circular
    planet positions at those epochs (exactly the §14 V1 construction).

    Raises :class:`~cyclerfinder.core.lambert.LambertConvergenceError` /
    :class:`~cyclerfinder.core.lambert.LambertGeometryError` if a leg's geometry
    is single-rev-Lambert-singular (recorded as a build failure by the caller).
    """
    inner, outer = bodies
    a_km = sol.a_au * AU_KM
    e = sol.e
    g = free_return_geometry(sol.a_au, e, bodies=bodies, mu=mu)
    nu_inner = g.nu[inner]
    nu_outer = g.nu[outer]
    n = sqrt(mu / a_km**3)  # rad/s

    def tof_between(nu0: float, nu1: float) -> float:
        return ((_true_to_mean(nu1, e) - _true_to_mean(nu0, e)) % (2.0 * pi)) / n

    t_e0 = sol.t0_sec
    t_m = t_e0 + tof_between(nu_inner, nu_outer)
    t_e1 = t_e0 + tof_between(nu_inner, 2.0 * pi - nu_inner)

    def st(body: str, t: float) -> tuple[np.ndarray, np.ndarray]:
        r, v = ephem.state(body, t)
        return np.asarray(r, dtype=np.float64), np.asarray(v, dtype=np.float64)

    r_e0, v_e0 = st(inner, t_e0)
    r_m, v_m = st(outer, t_m)
    r_e1, v_e1 = st(inner, t_e1)

    s1 = lambert(r_e0, r_m, t_m - t_e0, mu=mu, max_revs=0)[0]
    s2 = lambert(r_m, r_e1, t_e1 - t_m, mu=mu, max_revs=0)[0]
    v1_dep, v1_arr = np.asarray(s1.v1, np.float64), np.asarray(s1.v2, np.float64)
    v2_dep, v2_arr = np.asarray(s2.v1, np.float64), np.asarray(s2.v2, np.float64)

    encounters = [
        Encounter(inner, t_e0, r_e0, v_e0, v1_dep - v_e0, v1_dep - v_e0),
        Encounter(outer, t_m, r_m, v_m, v1_arr - v_m, v2_dep - v_m),
        Encounter(inner, t_e1, r_e1, v_e1, v2_arr - v_e1, v2_arr - v_e1),
    ]
    legs = [
        Leg(inner, outer, t_e0, t_m, v1_dep, v1_arr, 0, "single"),
        Leg(outer, inner, t_m, t_e1, v2_dep, v2_arr, 0, "single"),
    ]
    return Cycler([inner, outer, inner], period_sec, encounters, legs, sense="outbound")


def free_return_v1_mechanics(
    sol: FreeReturnClosureResult,
    ephem: Ephemeris,
    period_sec: float,
    *,
    bodies: tuple[str, str] = ("E", "M"),
    mu: float = MU_SUN_KM3_S2,
) -> FreeReturnV1Result:
    """Run §14 V1 mechanics (paths a+c) + the V_inf-continuity gate on *sol*.

    Reconstructs the closed E->M->E free-return arc (:func:`build_free_return_cycler`),
    runs :func:`crosscheck_code_paths` (reusing paths a and c verbatim), and
    awards ``v1_passed`` only when both §14 paths pass AND the Mars-flyby V_inf is
    continuous (the genome honesty gate). Never raises — a Lambert-singular leg is
    reported as ``built=False`` with the reason.
    """
    try:
        cycler = build_free_return_cycler(sol, ephem, period_sec, bodies=bodies, mu=mu)
    except (LambertConvergenceError, LambertGeometryError) as exc:
        return FreeReturnV1Result(
            built=False,
            lamberthub_passed=False,
            lamberthub_max_diff_mps=float("inf"),
            kepler_reprop_passed=False,
            kepler_reprop_max_residual_km=float("inf"),
            vinf_continuity_kms=float("inf"),
            vinf_continuous=False,
            v1_passed=False,
            detail=f"E->M->E reconstruction Lambert-singular: {type(exc).__name__}",
        )

    report = crosscheck_code_paths(cycler, ephem, mu=mu)
    a_pass = report.lamberthub.passed
    c_pass = report.kepler_reprop.passed

    mars = cycler.encounters[1]
    vinf_in = float(np.linalg.norm(mars.vinf_in))
    vinf_out = float(np.linalg.norm(mars.vinf_out))
    continuity = abs(vinf_in - vinf_out)
    continuous = continuity <= VINF_CONTINUITY_TOL_KMS

    v1 = a_pass and c_pass and continuous
    if v1:
        detail = "§14 V1 (paths a+c) pass on a closed, V_inf-continuous free-return arc"
    elif not continuous:
        detail = (
            f"reconstructed return leg breaks Mars V_inf continuity "
            f"({continuity:.2f} km/s > {VINF_CONTINUITY_TOL_KMS} km/s): the single "
            f"free-return ellipse does not close to Earth (multi-arc; needs phasing loops)"
        )
    else:
        detail = f"§14 path fail: lamberthub={a_pass} kepler_reprop={c_pass}"

    return FreeReturnV1Result(
        built=True,
        lamberthub_passed=a_pass,
        lamberthub_max_diff_mps=float(report.lamberthub.max_diff_mps),
        kepler_reprop_passed=c_pass,
        kepler_reprop_max_residual_km=float(report.kepler_reprop.max_residual_km),
        vinf_continuity_kms=continuity,
        vinf_continuous=continuous,
        v1_passed=v1,
        detail=detail,
    )


__all__ = [
    "VINF_CONTINUITY_TOL_KMS",
    "FreeReturnV1Result",
    "build_free_return_cycler",
    "free_return_v1_mechanics",
]
