"""True ballistic-cycler construction for the Hernandez 2017 EIGE triple cycler.

EIGE (Europa-Io-Ganymede-Europa) is the paper's **maintenance-ΔV demonstration**
cycler: ONE synodic period, ONE spacecraft revolution (Hernandez, Jones & Jesick
2017, AAS 17-608, Figure 5, pp.10-11). It is the structurally SIMPLEST Galilean
triple cycler — three Lambert legs over a single revolution at the 1:1 resonant
semi-major axis ``a = a_Ganymede`` — which is exactly why the paper uses it to
illustrate the ideal→real-ephemeris transition (its Fig-5 real-ephemeris instance is
ballistic for the first cycle, then grows to ~30 m/s over 10 repeat cycles).

WHY this module exists
----------------------
It is the EIGE analog of :mod:`cyclerfinder.search.eggie_ballistic` and the **positive
control** for the per-cycle maintenance-ΔV method (``nbody.jovian.chain_cycles``):
before quoting an EGGIE maintenance number we reproduce a paper-printed maintenance
example, the way [[feedback_verify_gauntlet_with_positive_control]] requires. This
module builds the *ideal-model* ballistic EIGE that seeds that real-ephemeris lane.

The construction (constraints, not post-hoc metrics)
----------------------------------------------------
Model EIGE as 3 Lambert conic legs between the circular-coplanar moon positions
(paper p.3 ideal model). The Europa departure phase is the gauge (fixed at 0); the
free variables are the Io and Ganymede inertial phases at ``t0`` plus the 3 leg ToFs.
The SOLVE residuals are the equal-in/out |V∞| ballistic-flyby constraints at the two
interior flybys (Io, Ganymede) and the Europa periodicity seam
(``|V∞_in(final)| = |V∞_out(depart)|``). Io and Ganymede each appear once, so — unlike
EGGIE's repeated Ganymede — there is no resonant-return equal-V∞ constraint. The
repeated Europa (depart + wrap) is self-consistent by the resonance: the tour period
is one ``T_sc = T_syn = 2·T_Europa``, so Europa returns to its departure phase exactly
([[feedback_constructed_tour_per_encounter_self_consistency]]).

That leaves 5 free variables and 3 ballistic residuals — a 2-DOF family. The two spare
DOF are pinned by **softly targeting the two sourced Fig-5 interior altitudes** (Io
2,817 km, Ganymede 13,180 km), exactly as :mod:`eggie_ballistic` softly targets the
sourced Table-4 V∞. These altitude targets are SOURCED inputs used to pin a unique,
reproducible member — not asserted discoveries. The **Europa flyby altitude is then a
free PREDICTION** (~1,323 km), the non-circular cross-check: it lands the same order as
Fig-5's printed 470 km (the residual difference is the ideal↔real-ephemeris gap, the
paper's own Fig-5 being real-ephemeris).

The verdict (the honest finding)
--------------------------------
A genuine **feasible ballistic EIGE exists** in the ideal model: equal-in/out |V∞| at
all three flybys, the cycle closed, all flyby altitudes inside the paper's 25-70,000 km
window, total flyby ΔV ~ 0. Its excess speeds are the **low-excess-speed 5-9 km/s
regime** (Europa ~8.70, Io ~5.14, Ganymede ~7.23 km/s) — the same navigation-viable
band as EGGIE, NOT the 12-16 km/s of the 1-synodic/2-rev EGIEIE (correcting a prior
scoping over-extrapolation). Topology: Europa departs INBOUND, Io inbound, Ganymede
outbound, Europa inbound (wrap), all within one revolution.

Sourced anchors (never our own computed values; ``feedback_golden_tests_sourced_only``):
the Fig-5 interior flyby altitudes Io 2,817 / Ganymede 13,180 / Europa 470 km and the
1-synodic/1-rev classification (AAS 17-608 pp.10-11), and the 25-70,000 km flyby-altitude
window (p.7). Digest:
``docs/notes/2026-06-26-digest-hernandez-2017-ieg-triple-cyclers-aas-17-608.md``.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

from cyclerfinder.core.lambert import lambert
from cyclerfinder.nbody.jovian import flyby_altitude_km, flyby_maneuver_dv
from cyclerfinder.search.eggie_ballistic import MU, moon_state
from cyclerfinder.search.resonant_conic import ideal_t_syn, resonant_sma

Vec3 = NDArray[np.float64]

SECONDS_PER_DAY: float = 86400.0

#: The EIGE flyby sequence (Europa-first by the paper's convention), with the trailing
#: Europa the periodicity wrap (start of the next cycle).
EIGE_SEQUENCE: tuple[str, ...] = ("Europa", "Io", "Ganymede", "Europa")

#: EIGE resonance (AAS 17-608 Fig 5: one synodic period, one revolution -> 1:1).
EIGE_N_SYN: int = 1
EIGE_N_REV: int = 1

#: Paper flyby-altitude window (km), AAS 17-608 p.7.
ALT_MIN_KM: float = 25.0
ALT_MAX_KM: float = 70000.0

#: Fig-5 sourced interior flyby altitudes (km) — AAS 17-608 pp.10-11. SOURCED, never
#: code-derived; used as the soft pin targets (Io, Ganymede) and the Europa cross-check.
FIG5_ALT_KM: dict[str, float] = {"Io": 2817.0, "Ganymede": 13180.0, "Europa": 470.0}

#: The 1:1 resonant spacecraft semi-major axis (km) = ideal Ganymede sma.
EIGE_RESONANT_SMA_KM: float = resonant_sma(EIGE_N_SYN, EIGE_N_REV, ideal_t_syn())


def _nrm(v: Vec3) -> float:
    return float(np.linalg.norm(v))


# --- Patched-conic leg construction ----------------------------------------------


def build_legs(
    phi_io: float, phi_gan: float, tofs_sec: Sequence[float]
) -> tuple[list[Vec3], list[Vec3]] | None:
    """Build the 3 Lambert legs of the EIGE tour in the ideal circular-coplanar model.

    Europa's departure phase is the gauge (fixed at 0); ``phi_io`` / ``phi_gan`` are the
    Io / Ganymede inertial phases at ``t0``. Each leg is a single-revolution (n_revs=0),
    prograde, short/long-way Lambert arc — the deterministic plan of the EIGE member
    (all three legs are sub-revolution arcs of the single 1-rev tour). Returns
    ``(vinf_out, vinf_in)`` lists indexed by node 0..3, where ``vinf_out[k]`` departs
    node ``k`` and ``vinf_in[k]`` arrives node ``k`` (``None`` at the boundary slots).
    Returns ``None`` if any leg is infeasible.
    """
    t = [0.0]
    for tof in tofs_sec:
        t.append(t[-1] + float(tof))
    phases = {"Europa": 0.0, "Io": phi_io, "Ganymede": phi_gan}
    rpos: list[Vec3] = []
    rvel: list[Vec3] = []
    for k, moon in enumerate(EIGE_SEQUENCE):
        r, v = moon_state(moon, phases[moon], t[k])
        rpos.append(r)
        rvel.append(v)
    vinf_out: list[Vec3 | None] = [None] * 4
    vinf_in: list[Vec3 | None] = [None] * 4
    for k in range(3):
        leg_tof = t[k + 1] - t[k]
        if leg_tof <= 0.0:
            return None
        try:
            sols = lambert(rpos[k], rpos[k + 1], leg_tof, mu=MU, max_revs=0, prograde=True)
        except Exception:
            return None
        match = [s for s in sols if s.n_revs == 0]
        if not match:
            return None
        sol = match[0]
        vinf_out[k] = np.asarray(sol.v1, dtype=np.float64) - rvel[k]
        vinf_in[k + 1] = np.asarray(sol.v2, dtype=np.float64) - rvel[k + 1]
    return vinf_out, vinf_in  # type: ignore[return-value]


def _smooth_alt_km(vinf_in: Vec3, vinf_out: Vec3, moon: str) -> float:
    """Patched-conic flyby altitude (km), capped finite for the soft target residual."""
    a = flyby_altitude_km(vinf_in, vinf_out, moon)
    if not math.isfinite(a) or a > 1.0e6:
        return 1.0e6
    return a


def ballistic_residual(x: NDArray[np.float64], *, alt_target_w: float = 0.0) -> NDArray[np.float64]:
    """EIGE ballistic-cycler residual vector (km/s); +2 entries if ``alt_target_w`` > 0.

    ``x = [phi_io, phi_gan, tof1, tof2, tof3]`` (ToFs in seconds). The 3 core residuals
    are the equal-in/out |V∞| ballistic constraints at Io and Ganymede and the Europa
    periodicity seam. With ``alt_target_w`` > 0 two weighted residuals softly pull the Io
    and Ganymede flyby altitudes toward the SOURCED Fig-5 values (Io 2,817 / Ganymede
    13,180 km) — the well-determining pin for the otherwise 2-DOF family.
    """
    built = build_legs(x[0], x[1], list(x[2:]))
    if built is None:
        return np.full(3 + (2 if alt_target_w else 0), 50.0)
    vo, vi = built
    res = [
        _nrm(vi[1]) - _nrm(vo[1]),  # Io ballistic
        _nrm(vi[2]) - _nrm(vo[2]),  # Ganymede ballistic
        _nrm(vi[3]) - _nrm(vo[0]),  # Europa periodicity seam
    ]
    if alt_target_w:
        a_io = _smooth_alt_km(vi[1], vo[1], "Io")
        a_gan = _smooth_alt_km(vi[2], vo[2], "Ganymede")
        res.append(alt_target_w * (a_io - FIG5_ALT_KM["Io"]) / 1000.0)
        res.append(alt_target_w * (a_gan - FIG5_ALT_KM["Ganymede"]) / 1000.0)
    return np.asarray(res, dtype=np.float64)


def ballistic_resnorm(x: NDArray[np.float64]) -> float:
    """Norm of the 3 core ballistic residuals (km/s) — the closure measure."""
    return float(np.linalg.norm(ballistic_residual(x, alt_target_w=0.0)))


# --- Result dataclass + evaluator ------------------------------------------------


@dataclass(frozen=True)
class BallisticEige:
    """A constructed EIGE tour evaluated against the ballistic-cycler constraints."""

    phi_io_rad: float
    phi_gan_rad: float
    tofs_days: tuple[float, float, float]
    #: Departure |V∞| at Europa, and arrival |V∞| at Io, Ganymede, Europa (km/s).
    vinf_kms: dict[str, float]
    #: Per-flyby maneuver ΔV (m/s), paper Eq 3-5 (zero iff ballistic & bend-feasible).
    flyby_dv_ms: dict[str, float]
    #: Per-flyby periapsis altitude (km).
    flyby_alt_km: dict[str, float]
    total_dv_ms: float
    ballistic_resnorm_kms: float
    seam_defect_kms: float
    all_feasible: bool


def evaluate(phi_io: float, phi_gan: float, tofs_sec: Sequence[float]) -> BallisticEige | None:
    """Evaluate a constructed EIGE: V∞, flyby ΔV/altitudes, ballistic closure."""
    built = build_legs(phi_io, phi_gan, list(tofs_sec))
    if built is None:
        return None
    vo, vi = built
    flybys = {
        "Io": (vi[1], vo[1], "Io"),
        "Ganymede": (vi[2], vo[2], "Ganymede"),
        "Europa": (vi[3], vo[0], "Europa"),
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
        "Io": _nrm(vi[1]),
        "Io_out": _nrm(vo[1]),
        "Ganymede": _nrm(vi[2]),
        "Ganymede_out": _nrm(vo[2]),
        "Europa_arr": _nrm(vi[3]),
    }
    x = np.array([phi_io, phi_gan, *tofs_sec])
    return BallisticEige(
        phi_io_rad=float(phi_io),
        phi_gan_rad=float(phi_gan),
        tofs_days=tuple(float(t) / SECONDS_PER_DAY for t in tofs_sec),  # type: ignore[arg-type]
        vinf_kms=vinf,
        flyby_dv_ms=dv_ms,
        flyby_alt_km=alt_km,
        total_dv_ms=float(sum(dv_ms.values())),
        ballistic_resnorm_kms=ballistic_resnorm(x),
        seam_defect_kms=abs(_nrm(vi[3]) - _nrm(vo[0])),
        all_feasible=all(feas),
    )


def refine(
    seed_x: Sequence[float] | NDArray[np.float64],
    *,
    alt_target_w: float = 0.05,
    max_nfev: int = 8000,
) -> BallisticEige:
    """Differential-correct the EIGE from ``seed_x`` (trf least-squares).

    The 3 core ballistic residuals are hard-driven to zero; with ``alt_target_w`` > 0 two
    soft residuals pull the Io and Ganymede altitudes toward the sourced Fig-5 values,
    pinning the 2-DOF family to a unique reproducible member. ToFs are bounded positive
    (and below 1.5·T_sc). Deterministic given the seed.
    """
    x0 = np.asarray(seed_x, dtype=np.float64)
    t_sc = 2.0 * math.pi * math.sqrt(EIGE_RESONANT_SMA_KM**3 / MU)
    lb = np.array([-np.inf, -np.inf, 60.0, 60.0, 60.0])
    ub = np.array([np.inf, np.inf, 1.5 * t_sc, 1.5 * t_sc, 1.5 * t_sc])
    sol = least_squares(
        ballistic_residual,
        np.clip(x0, lb, ub),
        kwargs={"alt_target_w": alt_target_w},
        method="trf",
        bounds=(lb, ub),
        max_nfev=max_nfev,
        xtol=1e-15,
        ftol=1e-15,
    )
    result = evaluate(sol.x[0], sol.x[1], tuple(sol.x[2:]))
    if result is None:  # pragma: no cover - the corrector stays in the feasible build region
        raise ValueError("EIGE refine produced an infeasible leg build")
    return result


# --- Documented converged seed (the feasible ballistic EIGE) ---------------------

#: The feasible ballistic EIGE member (topology Europa-in / Io-in / Ganymede-out),
#: pinned at the two sourced Fig-5 interior altitudes. Refining it is deterministic and
#: lands on: equal-in/out |V∞| at all 3 flybys, cycle closed, all altitudes in window,
#: total flyby ΔV ~ 0; V∞ Europa ~8.70, Io ~5.14, Ganymede ~7.23 km/s; predicted Europa
#: altitude ~1,323 km. Produced by the conic-seeded altitude-anchored corrector
#: (``scripts``/exploration); a range of resonant-conic seeds (e 0.61-0.635) all converge
#: to this same fixed point.
EIGE_BALLISTIC_SEED: tuple[float, ...] = (
    -0.19423305657554754,
    -3.8762436775947724,
    36446.99232771451,
    96567.0663540463,
    472939.87410505034,
)


def feasible_ballistic_eige() -> BallisticEige:
    """Construct the feasible ballistic EIGE — deterministic.

    Refines :data:`EIGE_BALLISTIC_SEED` to a closed, equal-in/out-|V∞|, bend-feasible
    (all altitudes in [25, 70000] km) ballistic EIGE with total flyby ΔV ~ 0, pinned at
    the sourced Fig-5 interior altitudes (Io 2,817 / Ganymede 13,180 km). The Europa
    flyby altitude (~1,323 km) is the free prediction (Fig-5 prints 470 km).
    """
    return refine(EIGE_BALLISTIC_SEED, alt_target_w=0.05)


def eige_tof_seed_days() -> tuple[float, float, float]:
    """The 3 EIGE leg ToFs (days) of the feasible ballistic member — the real-eph seed.

    The ToF seed handed to ``nbody.jovian.chain_cycles`` (with ``sequence`` = EIGE and a
    3-leg single-rev branch plan) to drive the real-ephemeris maintenance-ΔV lane.
    """
    e = feasible_ballistic_eige()
    return e.tofs_days


__all__ = [
    "ALT_MAX_KM",
    "ALT_MIN_KM",
    "EIGE_BALLISTIC_SEED",
    "EIGE_N_REV",
    "EIGE_N_SYN",
    "EIGE_RESONANT_SMA_KM",
    "EIGE_SEQUENCE",
    "FIG5_ALT_KM",
    "BallisticEige",
    "ballistic_residual",
    "ballistic_resnorm",
    "build_legs",
    "eige_tof_seed_days",
    "evaluate",
    "feasible_ballistic_eige",
    "moon_state",
]
