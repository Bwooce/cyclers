"""True ballistic-cycler construction for the Hernandez 2017 EGGIE triple cycler.

This is the #480 ballistic construction in the paper's ideal circular-coplanar
Galilean model (Hernandez, Jones & Jesick 2017, AAS 17-608, pp.3-7; digest
``docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md``).

WHY this module exists (the corrected #480 diagnosis)
----------------------------------------------------
The earlier resonant-conic seed (:mod:`cyclerfinder.search.resonant_conic`) pins
each node's *departure* V∞ to the conic magnitude but does NOT enforce the
equal-in/out V∞ ballistic-flyby constraint, so its Ganymede #2 collapses to ~4 km/s
and the Io flyby goes sub-surface (see
``docs/notes/2026-06-29-480-eggie-forward-verify-correction.md``). Matching the
conic |V∞| magnitude is *necessary but not sufficient*: a true ballistic cycler is a
chain of independent Jovian conic legs linked by flybys that PRESERVE |V∞| (equal
in/out) at each encounter, with the cycle closing.

The construction (constraints, not post-hoc metrics)
----------------------------------------------------
Model EGGIE (Europa-Ganymede-Ganymede-Io-Europa, 5 s/c revs over 4 synodic periods)
as 4 Lambert conic legs between the circular-coplanar moon positions. Free variables
are the departure phase (Ganymede & Io inertial angle at t0; Europa angle is the
gauge, fixed at 0) and the 4 leg ToFs. The SOLVE residuals are:

* ``g1, g2, g3`` — equal-in/out |V∞| at the 3 interior flybys (Ganymede, Ganymede,
  Io): the ballistic property (``|V∞_in| = |V∞_out|`` -> zero flyby ΔV).
* ``g4`` — Europa periodicity seam (``|V∞_in(final)| = |V∞_out(depart)|``): in the
  coplanar model a flyby can bend the V∞ vector by any in-plane angle, so the binding
  closure condition is the magnitude (the per-cycle inertial rotation is absorbed).
* ``g5`` — Ganymede resonant return (``|V∞_out(G1)| = |V∞_in(G2)|``, both on the
  single G->G arc): this is what makes the two Ganymede flybys occur at *equal* V∞
  (the paper's defining ballistic-cycler property, Table 4 both 7.07 km/s).

Bend feasibility (each flyby's required turn achievable with periapsis altitude in
the paper's 25-70,000 km window) is checked via the paper's Eq 3-5 maneuver model
(:func:`cyclerfinder.nbody.jovian.flyby_maneuver_dv`).

The two regimes found (the honest verdict — "math decides")
-----------------------------------------------------------
1. **A genuine feasible ballistic EGGIE exists** (GATE A): equal-in/out |V∞| at all
   four flybys, the two Ganymede flybys equal, the cycle closed, ALL flyby altitudes
   in the 25-70,000 km window, total flyby ΔV = 0.0 m/s. It sits at lower excess
   speed than Table 4 (Europa ~9.01, Ganymede ~6.76, Io ~6.57 km/s) -> a real
   ballistic EGGIE-topology member, not the Table-4 member. See
   :data:`FEASIBLE_BALLISTIC_SEED` / :func:`feasible_ballistic_eggie`.
2. **The Table-4 V∞ levels are reproduced exactly on the ballistic manifold**
   (Europa 9.12, both Ganymede 7.07, Io 8.38, seam closed, total ΔV ~ 0) — but the
   flyby bends the coplanar geometry then forces are ~180° reversals, i.e. the
   periapsis altitudes go SUB-SURFACE (infeasible). The feasible (gentle-bend)
   Table-4 altitudes the paper prints (653-6263 km) require the real-ephemeris
   conversion (slightly inclined/eccentric moons supply the out-of-plane B-plane
   freedom that relieves the bend); Table 4's real departure date (29-Sep-2020) is
   consistent with a real-ephemeris solution. See :data:`TABLE4_VINF_SEED` /
   :func:`table4_vinf_eggie`.

Sourced anchors (never our own computed values; ``feedback_golden_tests_sourced_only``)
are the Table-4 invariants in
``docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md``:
Europa 9.12, Ganymede 7.07 (both), Io 8.38 km/s; total ΔV 0.70 m/s; flyby-altitude
window 25-70,000 km. GATE A asserts model-intrinsic *physical properties* (equal
in/out |V∞|, in-window altitudes), not arbitrary goldens.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

from cyclerfinder.core.lambert import lambert
from cyclerfinder.nbody.jovian import flyby_maneuver_dv
from cyclerfinder.search.resonant_conic import (
    MU_JUPITER_KM3_S2,
    ideal_moon_smas,
)

Vec3 = NDArray[np.float64]

MU: float = MU_JUPITER_KM3_S2
SECONDS_PER_DAY: float = 86400.0

#: The EGGIE flyby sequence (Europa-first by the paper's convention).
EGGIE_SEQUENCE: tuple[str, ...] = ("Europa", "Ganymede", "Ganymede", "Io", "Europa")

#: Ideal-model Galilean semi-major axes (km) and circular mean motions (rad/s).
_SMAS: dict[str, float] = ideal_moon_smas()
_MEAN_MOTION: dict[str, float] = {m: math.sqrt(MU / _SMAS[m] ** 3) for m in _SMAS}

#: Paper flyby-altitude window (km), AAS 17-608 p.7.
ALT_MIN_KM: float = 25.0
ALT_MAX_KM: float = 70000.0

#: Table-4 sourced V∞ targets (km/s) — digest pp.10; sourced, never code-derived.
EGGIE_VINF_TABLE4: dict[str, float] = {"Europa": 9.12, "Ganymede": 7.07, "Io": 8.38}

#: A leg plan entry: (n_revs, branch). branch is "single" for n_revs==0, else
#: "low"/"high" (the two multi-rev Lambert branches).
LegPlan = tuple[int, str]
Plan = tuple[LegPlan, LegPlan, LegPlan, LegPlan]


# --- Ideal circular-coplanar moon ephemeris --------------------------------------


def moon_state(moon: str, phase0_rad: float, t_sec: float) -> tuple[Vec3, Vec3]:
    """Circular-coplanar moon ``(r, v)`` (km, km/s) at inertial phase ``phase0`` + n·t."""
    theta = phase0_rad + _MEAN_MOTION[moon] * t_sec
    r = _SMAS[moon] * np.array([math.cos(theta), math.sin(theta), 0.0])
    r_hat = r / np.linalg.norm(r)
    v_mag = math.sqrt(MU / _SMAS[moon])
    v = v_mag * np.array([-r_hat[1], r_hat[0], 0.0])
    return r, v


# --- Patched-conic leg construction ----------------------------------------------


def build_legs(
    phi_gan: float, phi_io: float, tofs_sec: tuple[float, ...] | list[float], plan: Plan
) -> tuple[list[Vec3], list[Vec3]] | None:
    """Build the 4 Lambert legs of the EGGIE tour in the ideal model.

    Europa's departure phase is the gauge (fixed at 0); ``phi_gan`` / ``phi_io`` are
    the Ganymede / Io inertial phases at ``t0``. Returns ``(vinf_out, vinf_in)`` lists
    indexed by node (0..4), where ``vinf_out[k]`` is the V∞ departing node ``k`` and
    ``vinf_in[k]`` the V∞ arriving node ``k`` (``None`` at the unused boundary slots).
    Returns ``None`` if any leg is infeasible (non-positive ToF / no Lambert solution
    on the planned rev/branch).
    """
    t = [0.0]
    for tof in tofs_sec:
        t.append(t[-1] + float(tof))
    phases = {"Europa": 0.0, "Ganymede": phi_gan, "Io": phi_io}
    rpos: list[Vec3] = []
    rvel: list[Vec3] = []
    for k, moon in enumerate(EGGIE_SEQUENCE):
        r, v = moon_state(moon, phases[moon], t[k])
        rpos.append(r)
        rvel.append(v)
    vinf_out: list[Vec3 | None] = [None] * 5
    vinf_in: list[Vec3 | None] = [None] * 5
    for k in range(4):
        leg_tof = t[k + 1] - t[k]
        if leg_tof <= 0.0:
            return None
        n_revs, branch = plan[k]
        try:
            sols = lambert(rpos[k], rpos[k + 1], leg_tof, mu=MU, max_revs=max(n_revs, 1))
        except Exception:
            return None
        match = [s for s in sols if s.n_revs == n_revs and (n_revs == 0 or s.branch == branch)]
        if not match:
            return None
        sol = match[0]
        vinf_out[k] = np.asarray(sol.v1, dtype=np.float64) - rvel[k]
        vinf_in[k + 1] = np.asarray(sol.v2, dtype=np.float64) - rvel[k + 1]
    return vinf_out, vinf_in  # type: ignore[return-value]


def _nrm(v: Vec3) -> float:
    return float(np.linalg.norm(v))


def ballistic_residual(
    x: NDArray[np.float64], plan: Plan, target_w: float = 0.0, *, include_seam: bool = True
) -> NDArray[np.float64]:
    """Ballistic-cycler residual vector (km/s); +3 entries if ``target_w`` > 0.

    ``x = [phi_gan, phi_io, tof1, tof2, tof3, tof4]`` (ToFs in seconds). The core
    residuals are the equal-in/out |V∞| ballistic constraints at the interior flybys
    (Ganymede, Ganymede, Io), the Ganymede resonant-return (equal-7.07) constraint,
    and — when ``include_seam`` — the Europa periodicity seam (``|V∞_in(final)| =
    |V∞_out(depart)|``). Dropping the seam isolates the interior G->G->I sub-tour (the
    seam is the binding constraint at Table-4 V∞ in the 2-D model). With ``target_w`` > 0
    three weighted residuals pull the Europa/Ganymede/Io departure |V∞| toward the
    sourced Table-4 levels (used to slide along the ballistic manifold to Table-4 V∞).
    """
    n_core = 5 if include_seam else 4
    built = build_legs(x[0], x[1], list(x[2:]), plan)
    if built is None:
        return np.full(n_core + (3 if target_w else 0), 10.0)
    vo, vi = built
    res = [
        _nrm(vi[1]) - _nrm(vo[1]),  # Ganymede #1 ballistic
        _nrm(vi[2]) - _nrm(vo[2]),  # Ganymede #2 ballistic
        _nrm(vi[3]) - _nrm(vo[3]),  # Io ballistic
    ]
    if include_seam:
        res.append(_nrm(vi[4]) - _nrm(vo[0]))  # Europa periodicity seam
    res.append(_nrm(vo[1]) - _nrm(vi[2]))  # Ganymede resonant return (-> equal V∞)
    if target_w:
        res += [
            target_w * (_nrm(vo[0]) - EGGIE_VINF_TABLE4["Europa"]),
            target_w * (_nrm(vo[1]) - EGGIE_VINF_TABLE4["Ganymede"]),
            target_w * (_nrm(vo[3]) - EGGIE_VINF_TABLE4["Io"]),
        ]
    return np.asarray(res, dtype=np.float64)


def ballistic_resnorm(x: NDArray[np.float64], plan: Plan, *, include_seam: bool = True) -> float:
    """Norm of the core ballistic residuals (km/s) — the closure measure."""
    n_core = 5 if include_seam else 4
    res = ballistic_residual(x, plan, target_w=0.0, include_seam=include_seam)
    return float(np.linalg.norm(res[:n_core]))


# --- Result dataclass + evaluator ------------------------------------------------


@dataclass(frozen=True)
class BallisticEggie:
    """A constructed EGGIE tour evaluated against the ballistic-cycler constraints."""

    phi_gan_rad: float
    phi_io_rad: float
    tofs_days: tuple[float, float, float, float]
    plan: Plan
    #: Departure |V∞| at Europa, and arrival |V∞| at G1, G2, Io, Europa (km/s).
    vinf_kms: dict[str, float]
    #: Per-flyby maneuver ΔV (m/s), paper Eq 3-5 (zero iff ballistic & bend-feasible).
    flyby_dv_ms: dict[str, float]
    #: Per-flyby periapsis altitude (km).
    flyby_alt_km: dict[str, float]
    total_dv_ms: float
    ballistic_resnorm_kms: float
    ganymede_equal_resid_kms: float
    seam_defect_kms: float
    all_feasible: bool


def evaluate(
    phi_gan: float, phi_io: float, tofs_sec: tuple[float, ...] | list[float], plan: Plan
) -> BallisticEggie | None:
    """Evaluate a constructed EGGIE: V∞, flyby ΔV/altitudes, ballistic closure."""
    built = build_legs(phi_gan, phi_io, list(tofs_sec), plan)
    if built is None:
        return None
    vo, vi = built
    flybys = {
        "G1": (vi[1], vo[1], "Ganymede"),
        "G2": (vi[2], vo[2], "Ganymede"),
        "Io": (vi[3], vo[3], "Io"),
        "E": (vi[4], vo[0], "Europa"),
    }
    dv_ms: dict[str, float] = {}
    alt_km: dict[str, float] = {}
    feas: list[bool] = []
    for key, (vinf_in, vinf_out, moon) in flybys.items():
        dv, alt, ok = flyby_maneuver_dv(
            vinf_in, vinf_out, moon, alt_min=ALT_MIN_KM, alt_max=ALT_MAX_KM
        )
        dv_ms[key] = dv
        alt_km[key] = alt
        feas.append(ok)
    vinf = {
        "Europa_dep": _nrm(vo[0]),
        "Ganymede1": _nrm(vi[1]),
        "Ganymede1_out": _nrm(vo[1]),
        "Ganymede2": _nrm(vi[2]),
        "Ganymede2_out": _nrm(vo[2]),
        "Io": _nrm(vi[3]),
        "Io_out": _nrm(vo[3]),
        "Europa_arr": _nrm(vi[4]),
    }
    x = np.array([phi_gan, phi_io, *tofs_sec])
    return BallisticEggie(
        phi_gan_rad=float(phi_gan),
        phi_io_rad=float(phi_io),
        tofs_days=tuple(float(t) / SECONDS_PER_DAY for t in tofs_sec),  # type: ignore[arg-type]
        plan=plan,
        vinf_kms=vinf,
        flyby_dv_ms=dv_ms,
        flyby_alt_km=alt_km,
        total_dv_ms=float(sum(dv_ms.values())),
        ballistic_resnorm_kms=ballistic_resnorm(x, plan),
        ganymede_equal_resid_kms=abs(_nrm(vo[1]) - _nrm(vi[2])),
        seam_defect_kms=abs(_nrm(vi[4]) - _nrm(vo[0])),
        all_feasible=all(feas),
    )


def refine(
    seed_x: Sequence[float] | NDArray[np.float64],
    plan: Plan,
    *,
    target_w: float = 0.0,
    include_seam: bool = True,
    max_nfev: int = 600,
) -> BallisticEggie:
    """Differential-correct the EGGIE from ``seed_x`` (trf least-squares).

    With ``target_w`` = 0 it drives the core ballistic residuals to zero (closest
    ballistic point to the seed). With ``target_w`` > 0 it additionally pulls the
    departure V∞ toward the sourced Table-4 levels (then re-evaluates the core
    ballistic closure). ``include_seam=False`` drops the Europa periodicity seam so the
    corrector isolates the interior G->G->I sub-tour. Deterministic given the seed.
    """
    x0 = np.asarray(seed_x, dtype=np.float64)
    if target_w:
        r1 = least_squares(
            ballistic_residual,
            x0,
            args=(plan, target_w),
            kwargs={"include_seam": include_seam},
            method="trf",
            max_nfev=max_nfev,
            xtol=1e-14,
            ftol=1e-14,
        )
        x0 = np.asarray(r1.x, dtype=np.float64)
    sol = least_squares(
        ballistic_residual,
        x0,
        args=(plan, 0.0),
        kwargs={"include_seam": include_seam},
        method="trf",
        max_nfev=max_nfev,
        xtol=1e-14,
        ftol=1e-14,
    )
    result = evaluate(sol.x[0], sol.x[1], tuple(sol.x[2:]), plan)
    if result is None:  # pragma: no cover - the corrector stays in the feasible build region
        raise ValueError("EGGIE refine produced an infeasible leg build")
    return result


# --- Documented converged seeds (the two #480 regimes) ---------------------------

#: GATE A — the feasible ballistic EGGIE (all altitudes in window, total ΔV ~ 0,
#: Ganymede equal). V∞ Europa ~9.01, Ganymede ~6.76, Io ~6.57 km/s (a real
#: EGGIE-topology member at lower excess speed than Table 4). Seed produced by the
#: feasibility-aware corrector (``scripts/_eggie_feas_480.py``); refining it is
#: deterministic and lands on total ΔV = 0.0.
FEASIBLE_BALLISTIC_PLAN: Plan = ((0, "single"), (1, "high"), (1, "low"), (1, "high"))
FEASIBLE_BALLISTIC_SEED: tuple[float, ...] = (
    2.7950329174276782,
    15.558378851451593,
    139020.11775361068,
    733105.4873988518,
    640292.4998159712,
    1233640.1161620617,
)

#: GATE B — the Table-4 V∞ member ON the ballistic manifold (Europa 9.12, both
#: Ganymede 7.07, Io 8.38, seam closed, total flyby ΔV ~ 0) whose coplanar flyby
#: bends are ~180° reversals -> SUB-SURFACE altitudes (infeasible). Demonstrates the
#: binding constraint: in the strict 2D ideal model the Table-4 V∞ cycler is ballistic
#: but its flybys are not bend-feasible. Seed from ``scripts/_eggie_feas_480.py``.
TABLE4_VINF_PLAN: Plan = ((0, "single"), (1, "low"), (1, "low"), (2, "low"))
TABLE4_VINF_SEED: tuple[float, ...] = (
    0.07336039489270582,
    17.887062380206224,
    49601.001213030504,
    1513497.3921992257,
    621613.2105097923,
    1041265.4073508935,
)

#: GATE B (refined) — the interior G->G->I sub-tour reproduces Table-4 V∞ (Europa 9.12,
#: both Ganymede 7.07, Io 8.38) with the 3 interior flybys EXACTLY ballistic AND all
#: interior altitudes in the 25-70000 km window (G1 ~1419, G2 ~2233, Io ~7177 km, the
#: paper's ballpark) — only the Europa periodicity SEAM stays open (Europa arrival
#: ~9.37 vs 9.12 departure; Europa flyby sub-surface). This localises the binding
#: constraint to the seam: the interior is feasible at Table-4 V∞; full periodic closure
#: in the strict 2D model is what fails (needs the real-ephemeris/3D B-plane freedom).
INTERIOR_TABLE4_PLAN: Plan = ((0, "single"), (1, "high"), (1, "low"), (2, "low"))
INTERIOR_TABLE4_SEED: tuple[float, ...] = (
    -3.546238161079704,
    -46.580104582133984,
    135799.16277979105,
    727708.8738576599,
    589043.755322358,
    1144979.4600162422,
)


def feasible_ballistic_eggie() -> BallisticEggie:
    """Construct the feasible ballistic EGGIE (GATE A) — deterministic.

    Refines :data:`FEASIBLE_BALLISTIC_SEED` to a closed, equal-in/out-|V∞|,
    bend-feasible (all altitudes in [25, 70000] km) ballistic EGGIE with total flyby
    ΔV = 0.0 m/s and the two Ganymede flybys at equal V∞.
    """
    return refine(FEASIBLE_BALLISTIC_SEED, FEASIBLE_BALLISTIC_PLAN, target_w=0.0)


def table4_vinf_eggie() -> BallisticEggie:
    """Construct the Table-4-V∞ EGGIE on the ballistic manifold (GATE B finding).

    Refines :data:`TABLE4_VINF_SEED` with the Table-4 V∞ pull: the result is ballistic
    (equal in/out |V∞|, Ganymede equal) AT the sourced Table-4 levels but with
    sub-surface (infeasible) flyby altitudes — the documented binding constraint of the
    strict 2D ideal model.
    """
    return refine(TABLE4_VINF_SEED, TABLE4_VINF_PLAN, target_w=0.5)


def interior_table4_eggie() -> BallisticEggie:
    """Construct the interior-feasible Table-4 EGGIE sub-tour (GATE B refined finding).

    Refines :data:`INTERIOR_TABLE4_SEED` with the seam dropped (``include_seam=False``):
    the 3 interior flybys (Ganymede, Ganymede, Io) are EXACTLY ballistic at the sourced
    Table-4 V∞ (Europa entry 9.12, both Ganymede 7.07, Io 8.38) with all interior
    altitudes inside the 25-70000 km window — but the Europa periodicity seam stays open
    (Europa arrival ~9.37 km/s, Europa flyby sub-surface). Pinpoints the seam as the
    binding constraint for full periodic closure in the 2D ideal model.
    """
    return refine(INTERIOR_TABLE4_SEED, INTERIOR_TABLE4_PLAN, target_w=0.5, include_seam=False)


__all__ = [
    "ALT_MAX_KM",
    "ALT_MIN_KM",
    "EGGIE_SEQUENCE",
    "EGGIE_VINF_TABLE4",
    "FEASIBLE_BALLISTIC_PLAN",
    "FEASIBLE_BALLISTIC_SEED",
    "INTERIOR_TABLE4_PLAN",
    "INTERIOR_TABLE4_SEED",
    "TABLE4_VINF_PLAN",
    "TABLE4_VINF_SEED",
    "BallisticEggie",
    "ballistic_residual",
    "ballistic_resnorm",
    "build_legs",
    "evaluate",
    "feasible_ballistic_eggie",
    "interior_table4_eggie",
    "moon_state",
    "refine",
    "table4_vinf_eggie",
]
