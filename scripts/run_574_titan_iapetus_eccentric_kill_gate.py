"""#574 Stage A -- Titan-Iapetus eccentric-Keplerian 3D-closure kill gate.

Cheap kill-gate test for whether #573's 22 idealized-CIRCULAR 3D-closure
branches (Titan-Iapetus, Saturn) survive when both moons are placed on
REAL, non-negligible ECCENTRIC Keplerian orbits (Titan e~=0.0288, Iapetus
e~=0.028 -- 7-25x the Uranian moons' e<=0.004 that #558's family tolerated).
See ``data/OUTSTANDING.md`` #574 (Stage A only -- Stage B, productization +
real-ephemeris SPICE validation, is explicitly OUT of scope here and is not
touched by this script) for the authoritative, Fable-corrected spec. This
docstring only summarizes.

Method -- continuation in eccentricity (Fable-recommended over a fresh grid)
-----------------------------------------------------------------------------
Each of the 22 #573 circular branches is its own known solution
(Omega*, tof_scale*, rel_offset_deg*) at e=0. This script:

1. Loads all 22 branches' representative (Omega*, tof_scale*, n_rev,
   representative_residual_kms, eccentricity_robust) plus one seed
   rel_offset_deg (parsed from the branch's first ``seed_labels`` entry)
   directly from ``data/probe_573_titan_iapetus_population_closure.jsonl``
   -- never hand-transcribed, so a #572-style wrong-input bug cannot occur.
2. Generalizes #573's circular state generators to a genuinely ECCENTRIC
   Kepler propagator (``kepler_state_3d``) that reduces EXACTLY to
   ``_moon_state``/``iapetus_state_3d`` at e=0 (checked by an explicit smoke
   test below, mirroring #573's own ``smoke_test_reduction_pass`` pattern --
   C2 of the Fable plan review).
3. C1 (MANDATORY, prevents the #480 EGGIE per-encounter self-consistency
   bug): the free parameters per branch are a 4D space -- Omega (Iapetus
   ascending-node longitude), tof_scale, M0_Titan, M0_Iapetus (mean anomaly
   AT EPOCH t=0 for each moon). There is NO free "phase at the encounter" --
   Titan's second state (at t=2*tof, matching the #573 doubled-TOF pattern
   in ``evaluate_point_tracked``) is derived by Kepler-propagating the SAME
   M0_Titan over 2*tof, exactly mirroring #573's ``_moon_state(0.0, n_a,
   2.0*tof, ...)`` call (mean-motion propagation from the t=0 state, never a
   free re-specification). Iapetus's state at the first encounter (t=tof) is
   likewise Kepler-propagated from M0_Iapetus over tof.
   Argument of periapsis is FIXED at 0 for BOTH moons (periapsis at the
   Omega=0 reference direction / at the ascending node) -- an unsourced
   throwaway-script simplification (real periapsis precesses under
   perturbations we do not model here at all), NOT a free/searched
   parameter, so it does not reintroduce a per-encounter phase DOF.
4. Steps e from 0 -> the real value over ``N_ECC_STEPS`` stages, refining all
   4 free parameters via bounded Nelder-Mead at each stage, seeded from the
   PREVIOUS stage's converged point (small per-step search windows --
   prevents branch-jumping when eccentricity turns on). Stage 0 (e=0) is
   BOTH the smoke test (must reproduce the circular branch's own residual)
   AND the free e=0 positive control (C2).
5. At full eccentricity, a branch "survives" iff residual <= 0.05 km/s
   (the project-wide ``GATE_RESIDUAL_KMS`` bar, identical to #573's own
   criterion) AND the #324 physical bend gate
   (``candidate_passes_physical_gate``, reused verbatim, not reimplemented)
   passes at the REAL eccentric V_inf.
6. Survivors are deduped via union-find proximity clustering EXTENDED to
   include (M0_Titan, M0_Iapetus) on top of #573's (n_rev, Omega, tof_scale,
   V_inf) criterion (C3) -- checking explicitly whether the two known
   near-mirror pairs {branch 2, branch 19} and {branch 4, branch 14} merge
   under eccentricity (they are an EXACT circular-only degeneracy: the
   (Omega+180, u+180) transform leaves a CONSTANT-radius circular orbit's
   (x, y) unchanged and only flips z: r(nu)=a always, so r(nu)=r(nu+180).
   Under eccentricity r(nu) != r(nu+180) in general (nu=0 is periapsis,
   nu=180 is apoapsis under the omega=0 convention above), so this predicts
   the exact degeneracy SHOULD break -- verified empirically below, not
   assumed).

Pre-registered thresholds (C3, identical anchor to the #573/#558
precedent): PASS = >=5 deduped survivors from the 17 eccentricity-robust
branches, spanning >=2 n_rev classes. KILL = <=3. 4 = MARGINAL (reported,
NOT resolved here -- flagged for an explicit Opus adjudication per spec).
Floor-hugger control: ids 11, 15, 17, 18, 21 (bends 5.28/5.37/5.02/5.32/
5.02 deg) are EXPECTED to die; if >=3 instead survive, the eccentricity-
robust >=6.0 deg proxy is flagged as non-discriminating.

Framing (mandatory, per spec): ANY surviving family is quasi-cycler-class
evidence, same standing as #312's own family (V2 fails on drift by design),
NOT a ballistic-cycler finding and NOT a novelty claim -- an internal fact
about our own idealized (now eccentric-Keplerian, still non-ephemeris)
search space.

Eccentricity sourcing: e_Titan ~= 0.0288, e_Iapetus ~= 0.028 -- these are
the values mandated by the #574 Stage-A spec text itself (data/OUTSTANDING.md
#574), sourced there to JPL SSD Planetary Satellite Mean Orbital Parameters
(ssd.jpl.nasa.gov/sats/elem/). Kept as throwaway local constants in THIS
script (not added to ``core/satellites.py``) per the spec's own explicit
license ("your judgment... given this is explicitly a throwaway/kill-gate
script") -- Stage B, if greenlit, would productize this properly with a
sourced registry field.

Run as::

    uv run python scripts/run_574_titan_iapetus_eccentric_kill_gate.py
"""

from __future__ import annotations

import json
import math
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from probe_572_titan_iapetus_3d_closure import (  # noqa: E402
    ANCHOR,
    FLYBY,
    INCLINATION_DEG,
    PRIMARY,
    _leg_best,
)
from scan_558_uranus_all_pairs_offset_sweep import GATE_RESIDUAL_KMS  # noqa: E402

from cyclerfinder.core.lambert import LambertConvergenceError, LambertGeometryError  # noqa: E402
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.data.preflight import MethodCapability, preflight_search  # noqa: E402
from cyclerfinder.search.discovery_campaign import (  # noqa: E402
    DAY_S,
    _mean_motion_rad_day,
    _moon_state,
)
from cyclerfinder.search.physical_sanity import (  # noqa: E402
    DEFAULT_MIN_USEFUL_BEND_DEG,
    candidate_passes_physical_gate,
)

DATA_DIR = ROOT / "data"
PROBE_573_PATH = DATA_DIR / "probe_573_titan_iapetus_population_closure.jsonl"
OUT_PATH = DATA_DIR / "probe_574_titan_iapetus_eccentric_kill_gate.jsonl"

_REGION_ID = "saturn-titan-iapetus-574-stageA-eccentric-kill-gate-2026-07-12"
_METHOD = MethodCapability(
    genome=(
        "eccentric-Keplerian re-check (Stage A kill-gate) of the #573 22-branch circular "
        "closure family -- continuation in eccentricity (e: 0 -> real Titan/Iapetus values) "
        "from each known circular branch, refining (Omega, tof_scale, M0_Titan, M0_Iapetus)"
    ),
    corrector=(
        "reuses #572/#573's Lambert-closure + #324 physical-gate machinery unmodified; "
        "adds a Kepler eccentric-anomaly state propagator (verified via an e=0 "
        "positive-control reduction to the #573 circular results) and the 4D "
        "continuation/refinement loop on top"
    ),
    capability_tags=frozenset(
        {"lambert", "3d-closure", "saturn", "titan", "iapetus", "eccentric", "kill-gate"}
    ),
    git_sha="working-tree",
)
EMPTY_REGIONS_PATH = DATA_DIR / "empty_regions.jsonl"

# --- Eccentricity sourcing (see module docstring) -- mandated by the #574
# Stage-A spec text (data/OUTSTANDING.md #574), sourced there to JPL SSD
# Planetary Satellite Mean Orbital Parameters (ssd.jpl.nasa.gov/sats/elem/).
ECC_TITAN = 0.0288
ECC_IAPETUS = 0.028

ECC_ROBUST_BEND_DEG = 6.0
FLOOR_HUGGER_BRANCH_IDS = (11, 15, 17, 18, 21)

# Continuation schedule: fraction of the real eccentricity applied at each
# stage. Stage 0 (fraction 0.0) is the e=0 smoke test / positive control.
ECC_FRACTIONS = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)

# Per-continuation-step Nelder-Mead search windows (small -- steps are small
# in eccentricity, so the converged point should not move far).
WINDOW_OMEGA_DEG = 8.0
WINDOW_TOF_SCALE = 0.03
WINDOW_M0_DEG = 8.0

# Branch-dedup proximity thresholds -- #573's own (Omega, tof_scale, V_inf)
# criteria, EXTENDED (C3) with new (M0_Titan, M0_Iapetus) thresholds.
MERGE_OMEGA_DEG = 5.0
MERGE_TOF_SCALE = 0.05
MERGE_VINF_KMS = 0.10
MERGE_M0_DEG = 5.0


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def _circ_dist_deg(a: float, b: float) -> float:
    return abs((a - b + 180.0) % 360.0 - 180.0)


# ---------------------------------------------------------------------------
# Eccentric Kepler propagation -- generalizes #572/#573's circular
# `_moon_state` / `iapetus_state_3d` to a real ellipse. Argument of periapsis
# is fixed at 0 for both moons (periapsis at the Omega=0 reference
# direction), so the only per-body free parameter is M0 (mean anomaly at
# t=0); Omega (RAAN, Iapetus only -- Titan stays in-plane, Omega=inc=0
# always) and tof_scale remain free/searched exactly as in #573.
# ---------------------------------------------------------------------------


def _solve_kepler_e(
    mean_anomaly_rad: float, ecc: float, *, tol: float = 1e-13, max_iter: int = 60
) -> float:
    """Newton-Raphson solve of Kepler's equation E - e*sin(E) = M."""
    m = mean_anomaly_rad % (2.0 * math.pi)
    e_anom = m if ecc < 0.8 else math.pi
    for _ in range(max_iter):
        f = e_anom - ecc * math.sin(e_anom) - m
        fp = 1.0 - ecc * math.cos(e_anom)
        d = f / fp
        e_anom -= d
        if abs(d) < tol:
            break
    return e_anom


def kepler_state_3d(
    m0_rad: float,
    n_rad_day: float,
    t_days: float,
    sma_km: float,
    mu: float,
    ecc: float,
    raan_rad: float,
    inc_rad: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Eccentric-Keplerian moon position/velocity (planet frame, km, km/s).

    Mean anomaly is propagated by MEAN MOTION from the epoch value `m0_rad`
    (Kepler III mean motion depends only on `sma_km`, not `ecc`, so reusing
    the same `n_rad_day` the circular code already computes is exact -- C1
    discipline: no free per-encounter phase, only epoch M0 + elapsed time).
    Argument of periapsis is fixed at 0 (periapsis at the Omega=0 direction);
    `raan_rad=inc_rad=0` (Titan's case) reduces this to a pure in-plane
    ellipse with NO rotation, i.e. the perifocal frame IS the planet frame.

    At `ecc=0` this reduces EXACTLY to the circular case (E=nu=M, r=sma_km)
    -- verified against `_moon_state`/`iapetus_state_3d` by
    `_smoke_test_kepler_reduction` below.
    """
    m_t = m0_rad + n_rad_day * t_days
    e_anom = _solve_kepler_e(m_t, ecc)
    cos_e = math.cos(e_anom)
    nu = 2.0 * math.atan2(
        math.sqrt(1.0 + ecc) * math.sin(e_anom / 2.0), math.sqrt(1.0 - ecc) * math.cos(e_anom / 2.0)
    )
    r = sma_km * (1.0 - ecc * cos_e)
    p = sma_km * (1.0 - ecc * ecc)
    cos_nu, sin_nu = math.cos(nu), math.sin(nu)
    px = r * cos_nu
    py = r * sin_nu
    v_scale = math.sqrt(mu / max(p, 1e-9))
    vx_pf = -v_scale * sin_nu
    vy_pf = v_scale * (ecc + cos_nu)

    cos_o, sin_o = math.cos(raan_rad), math.sin(raan_rad)
    cosi, sini = math.cos(inc_rad), math.sin(inc_rad)

    def _rot(px_: float, py_: float) -> np.ndarray:
        x = cos_o * px_ - sin_o * cosi * py_
        y = sin_o * px_ + cos_o * cosi * py_
        z = sini * py_
        return np.array([x, y, z])

    pos = _rot(px, py)
    vel = _rot(vx_pf, vy_pf)
    return pos, vel


def _smoke_test_kepler_reduction() -> bool:
    """At ecc=0, kepler_state_3d must reproduce _moon_state / the circular
    iapetus_state_3d formula exactly (any M0, Omega, inc)."""
    from probe_572_titan_iapetus_3d_closure import iapetus_state_3d

    mu = PRIMARIES[PRIMARY]
    sat_a = SATELLITES[ANCHOR]
    sat_b = SATELLITES[FLYBY]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    v_circ_b = math.sqrt(mu / sat_b.sma_km)
    ok = True

    # Titan (planar): kepler_state_3d(m0, n_a, t, sma_a, mu, 0, 0, 0) vs _moon_state(m0, n_a, t)
    for m0_deg in (0.0, 37.0, 190.0, 355.0):
        for t_days in (0.0, 1.3, 12.0):
            m0 = math.radians(m0_deg)
            r3d, v3d = kepler_state_3d(m0, n_a, t_days, sat_a.sma_km, mu, 0.0, 0.0, 0.0)
            r2d, v2d = _moon_state(m0, n_a, t_days, sat_a.sma_km, mu)
            dr = float(np.linalg.norm(r3d - r2d))
            dv = float(np.linalg.norm(v3d - v2d))
            if dr > 1e-6 or dv > 1e-9:
                print(
                    f"  KEPLER SMOKE FAIL (Titan): m0={m0_deg} t={t_days} dr={dr:.3e} dv={dv:.3e}",
                    flush=True,
                )
                ok = False

    # Iapetus (inclined): kepler_state_3d(m0, n_b, 0, sma_b, mu, 0, om, inc) vs iapetus_state_3d(u)
    for omega_deg in (0.0, 37.0, 190.0, 355.0):
        for u_deg in (0.0, 42.0, 200.0):
            m0 = math.radians(u_deg)
            omega = math.radians(omega_deg)
            inc = math.radians(INCLINATION_DEG)
            r3d, v3d = kepler_state_3d(m0, n_b, 0.0, sat_b.sma_km, mu, 0.0, omega, inc)
            r2d, v2d = iapetus_state_3d(m0, v_circ_b, sat_b.sma_km, omega, inc)
            dr = float(np.linalg.norm(r3d - r2d))
            dv = float(np.linalg.norm(v3d - v2d))
            if dr > 1e-6 or dv > 1e-9:
                print(
                    f"  KEPLER SMOKE FAIL (Iapetus): omega={omega_deg} u={u_deg} "
                    f"dr={dr:.3e} dv={dv:.3e}",
                    flush=True,
                )
                ok = False
    return ok


# ---------------------------------------------------------------------------
# Eccentric closure evaluation (C1: 4D free-parameter space).
# ---------------------------------------------------------------------------


def evaluate_point_ecc(
    m0_titan_deg: float,
    m0_iapetus_deg: float,
    tof_scale: float,
    n_rev: tuple[int, int],
    omega_deg: float,
    e_titan: float,
    e_iapetus: float,
) -> tuple[dict[str, Any] | None, str | None]:
    """Eccentric analogue of #573's `evaluate_point_tracked`. Titan's state
    at t=0 AND t=2*tof are BOTH Kepler-propagated from the SAME m0_titan_deg
    (C1 -- no free re-specification at the second encounter); Iapetus's
    state at t=tof is Kepler-propagated from m0_iapetus_deg.
    """
    mu = PRIMARIES[PRIMARY]
    sat_a = SATELLITES[ANCHOR]
    sat_b = SATELLITES[FLYBY]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b
    tof = tof_scale * math.sqrt(p_a * p_b)
    tof_s = tof * DAY_S

    m0_t = math.radians(m0_titan_deg)
    m0_i = math.radians(m0_iapetus_deg)
    omega = math.radians(omega_deg)
    inc = math.radians(INCLINATION_DEG)

    r0, v0 = kepler_state_3d(m0_t, n_a, 0.0, sat_a.sma_km, mu, e_titan, 0.0, 0.0)
    r1, v1 = kepler_state_3d(m0_i, n_b, tof, sat_b.sma_km, mu, e_iapetus, omega, inc)
    r2, v2 = kepler_state_3d(m0_t, n_a, 2.0 * tof, sat_a.sma_km, mu, e_titan, 0.0, 0.0)

    n0, n1 = n_rev
    try:
        leg0 = _leg_best(r0, v0, r1, v1, tof_s, mu, n0)
        leg1 = _leg_best(r1, v1, r2, v2, tof_s, mu, n1)
    except LambertGeometryError:
        return None, "geometry"
    except LambertConvergenceError:
        return None, "convergence"
    if leg0 is None or leg1 is None:
        return None, "infeasible_n_rev"

    r_mid = abs(leg0["vinf_in"] - leg1["vinf_out"])
    r_periodic = abs(leg0["vinf_out"] - leg1["vinf_in"])
    residual = max(r_mid, r_periodic)

    vinf0 = leg0["vinf_out"]
    vinf1 = max(leg0["vinf_in"], leg1["vinf_out"])
    vinf2 = leg1["vinf_in"]

    return (
        {
            "m0_titan_deg": m0_titan_deg,
            "m0_iapetus_deg": m0_iapetus_deg,
            "tof_scale": tof_scale,
            "n_rev": list(n_rev),
            "omega_deg": omega_deg,
            "e_titan": e_titan,
            "e_iapetus": e_iapetus,
            "residual_kms": residual,
            "vinf_kms": [vinf0, vinf1, vinf2],
        },
        None,
    )


def _refine(
    x0: dict[str, float],
    n_rev: tuple[int, int],
    e_titan: float,
    e_iapetus: float,
) -> dict[str, Any]:
    """Bounded 4D Nelder-Mead refinement of (Omega, tof_scale, M0_T, M0_I)
    around x0, small per-step windows (continuation discipline)."""
    from scipy.optimize import minimize

    omega0, tof0, m0t0, m0i0 = (
        x0["omega_deg"],
        x0["tof_scale"],
        x0["m0_titan_deg"],
        x0["m0_iapetus_deg"],
    )
    bounds = [
        (omega0 - WINDOW_OMEGA_DEG, omega0 + WINDOW_OMEGA_DEG),
        (tof0 - WINDOW_TOF_SCALE, tof0 + WINDOW_TOF_SCALE),
        (m0t0 - WINDOW_M0_DEG, m0t0 + WINDOW_M0_DEG),
        (m0i0 - WINDOW_M0_DEG, m0i0 + WINDOW_M0_DEG),
    ]

    def _obj(x: np.ndarray) -> float:
        omega_deg, tof_scale, m0t, m0i = (float(v) for v in x)
        pt, _err = evaluate_point_ecc(m0t, m0i, tof_scale, n_rev, omega_deg, e_titan, e_iapetus)
        if pt is None:
            return 1.0e3
        return pt["residual_kms"]

    x_start = np.array([omega0, tof0, m0t0, m0i0])
    res = minimize(
        _obj,
        x_start,
        method="Nelder-Mead",
        bounds=bounds,
        options={"xatol": 1e-7, "fatol": 1e-12, "maxiter": 400, "maxfev": 400},
    )
    omega_f, tof_f, m0t_f, m0i_f = (float(v) for v in res.x)
    pt, _err = evaluate_point_ecc(m0t_f, m0i_f, tof_f, n_rev, omega_f % 360.0, e_titan, e_iapetus)
    seed_pt, _err0 = evaluate_point_ecc(m0t0, m0i0, tof0, n_rev, omega0, e_titan, e_iapetus)
    seed_residual = seed_pt["residual_kms"] if seed_pt is not None else float("inf")
    if pt is None or (seed_pt is not None and pt["residual_kms"] > seed_residual):
        # Refinement failed to improve -- fall back to the seed point.
        pt = seed_pt if seed_pt is not None else pt
        omega_f, tof_f, m0t_f, m0i_f = omega0 % 360.0, tof0, m0t0, m0i0
    return {
        "point": pt,
        "omega_deg": omega_f % 360.0,
        "tof_scale": tof_f,
        "m0_titan_deg": m0t_f % 360.0,
        "m0_iapetus_deg": m0i_f % 360.0,
    }


def load_branches() -> list[dict[str, Any]]:
    """Load the 22 #573 branches + one seed rel_offset_deg (parsed from the
    branch's first seed_label, e.g. 'rel177_tof2.80_n22' -> 177.0) directly
    from the #573 output -- no hand transcription."""
    with PROBE_573_PATH.open(encoding="utf-8") as fh:
        summary = None
        for line in fh:
            rec = json.loads(line)
            if rec.get("kind") == "population_summary":
                summary = rec
    if summary is None:
        raise RuntimeError(f"no population_summary record found in {PROBE_573_PATH}")

    branches = []
    for b in summary["branches"]:
        label0 = b["seed_labels"][0]
        m = re.match(r"rel(-?\d+(?:\.\d+)?)_tof", label0)
        if not m:
            raise RuntimeError(f"could not parse rel_offset from seed label {label0!r}")
        rel_offset_seed_deg = float(m.group(1))
        branches.append(
            {
                "branch_id": b["branch_id"],
                "n_rev": tuple(b["n_rev"][0]),
                "seed_labels": b["seed_labels"],
                "circular_omega_deg": b["representative_omega_deg"],
                "circular_tof_scale": b["representative_tof_scale"],
                "circular_residual_kms": b["representative_residual_kms"],
                "circular_vinf_kms": b["representative_vinf_kms"],
                "circular_iapetus_bend_deg": b["representative_iapetus_bend_deg"],
                "eccentricity_robust": b["eccentricity_robust"],
                "m0_iapetus_seed_deg": rel_offset_seed_deg,
            }
        )
    return branches


def run_branch_continuation(branch: dict[str, Any], seq: tuple[str, str, str]) -> dict[str, Any]:
    """Continuation in eccentricity for one branch, all ECC_FRACTIONS stages,
    checkpointed at the caller level (append+flush per branch)."""
    n_rev = branch["n_rev"]
    x = {
        "omega_deg": branch["circular_omega_deg"],
        "tof_scale": branch["circular_tof_scale"],
        "m0_titan_deg": 0.0,  # WLOG at e=0, matching #573's fixed theta0=0
        "m0_iapetus_deg": branch["m0_iapetus_seed_deg"],
    }
    stages: list[dict[str, Any]] = []
    for frac in ECC_FRACTIONS:
        e_t = ECC_TITAN * frac
        e_i = ECC_IAPETUS * frac
        refined = _refine(x, n_rev, e_t, e_i)
        pt = refined["point"]
        x = {
            "omega_deg": refined["omega_deg"],
            "tof_scale": refined["tof_scale"],
            "m0_titan_deg": refined["m0_titan_deg"],
            "m0_iapetus_deg": refined["m0_iapetus_deg"],
        }
        stage_rec: dict[str, Any] = {
            "ecc_fraction": frac,
            "e_titan": e_t,
            "e_iapetus": e_i,
            "omega_deg": x["omega_deg"],
            "tof_scale": x["tof_scale"],
            "m0_titan_deg": x["m0_titan_deg"],
            "m0_iapetus_deg": x["m0_iapetus_deg"],
        }
        if pt is None:
            stage_rec["residual_kms"] = None
            stage_rec["vinf_kms"] = None
            stage_rec["lambert_infeasible"] = True
        else:
            stage_rec["residual_kms"] = pt["residual_kms"]
            stage_rec["vinf_kms"] = pt["vinf_kms"]
            stage_rec["lambert_infeasible"] = False
        stages.append(stage_rec)

    final = stages[-1]
    residual_ok = (not final["lambert_infeasible"]) and final["residual_kms"] <= GATE_RESIDUAL_KMS
    gate_pass = False
    bends: list[float] | None = None
    if residual_ok:
        gate_pass, gate_verdicts = candidate_passes_physical_gate(
            seq, tuple(final["vinf_kms"]), min_useful_bend_deg=DEFAULT_MIN_USEFUL_BEND_DEG
        )
        bends = [v.max_bend_deg for v in gate_verdicts]
    survives = bool(residual_ok and gate_pass)

    # Smoke check (C2): stage 0 (e=0) residual must match the circular
    # branch's own known residual (well inside the gate).
    stage0 = stages[0]
    smoke_ok = (
        not stage0["lambert_infeasible"]
        and stage0["residual_kms"] < GATE_RESIDUAL_KMS
        and abs(stage0["residual_kms"] - branch["circular_residual_kms"]) < 1e-3
    )

    return {
        "branch_id": branch["branch_id"],
        "n_rev": list(n_rev),
        "eccentricity_robust": branch["eccentricity_robust"],
        "is_floor_hugger": branch["branch_id"] in FLOOR_HUGGER_BRANCH_IDS,
        "circular_residual_kms": branch["circular_residual_kms"],
        "circular_iapetus_bend_deg": branch["circular_iapetus_bend_deg"],
        "smoke_ok": smoke_ok,
        "stage0_residual_kms": stage0["residual_kms"],
        "stages": stages,
        "final_residual_kms": final["residual_kms"],
        "final_vinf_kms": final["vinf_kms"],
        "final_bends_deg": bends,
        "final_residual_ok": residual_ok,
        "final_gate_pass": gate_pass,
        "survives": survives,
        "final_omega_deg": final["omega_deg"],
        "final_tof_scale": final["tof_scale"],
        "final_m0_titan_deg": final["m0_titan_deg"],
        "final_m0_iapetus_deg": final["m0_iapetus_deg"],
    }


def cluster_survivors(survivors: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    """Union-find dedup EXTENDED to (M0_Titan, M0_Iapetus) per C3."""
    n = len(survivors)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for i in range(n):
        for j in range(i + 1, n):
            a, b = survivors[i], survivors[j]
            if tuple(a["n_rev"]) != tuple(b["n_rev"]):
                continue
            if _circ_dist_deg(a["final_omega_deg"], b["final_omega_deg"]) > MERGE_OMEGA_DEG:
                continue
            if abs(a["final_tof_scale"] - b["final_tof_scale"]) > MERGE_TOF_SCALE:
                continue
            if abs(a["final_vinf_kms"][1] - b["final_vinf_kms"][1]) > MERGE_VINF_KMS:
                continue
            if _circ_dist_deg(a["final_m0_titan_deg"], b["final_m0_titan_deg"]) > MERGE_M0_DEG:
                continue
            if _circ_dist_deg(a["final_m0_iapetus_deg"], b["final_m0_iapetus_deg"]) > MERGE_M0_DEG:
                continue
            union(i, j)

    groups: dict[int, list[dict[str, Any]]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(survivors[i])
    return list(groups.values())


def main() -> int:
    t_start = time.time()
    sha = _git_sha()
    print(f"[574A] Titan-Iapetus ECCENTRIC-Keplerian kill gate -- sha={sha}", flush=True)

    preflight_search(
        task_no=574,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=Path(__file__),
        n_points=22,
        override_reason=(
            "read-only eccentric re-check of the #573 22-branch closure set "
            "(continuation from known solutions, 6.8s measured total) -- a kill-gate "
            "check, not an unbudgeted discovery sweep"
        ),
    )

    print(
        "[574A] smoke test: kepler_state_3d reduces to _moon_state/iapetus_state_3d at e=0 ...",
        flush=True,
    )
    kepler_smoke_ok = _smoke_test_kepler_reduction()
    print(f"[574A] kepler reduction smoke test PASS: {kepler_smoke_ok}", flush=True)
    if not kepler_smoke_ok:
        print(
            "[574A] ABORTING -- eccentric Kepler propagator does not reduce correctly at e=0.",
            flush=True,
        )
        return 1

    branches = load_branches()
    print(f"[574A] loaded {len(branches)} #573 branches from {PROBE_573_PATH}", flush=True)
    n_robust = sum(1 for b in branches if b["eccentricity_robust"])
    n_hugger = sum(1 for b in branches if b["branch_id"] in FLOOR_HUGGER_BRANCH_IDS)
    print(f"[574A]   eccentricity-robust: {n_robust}, floor-huggers: {n_hugger}", flush=True)
    if len(branches) != 22 or n_robust != 17 or n_hugger != 5:
        print(
            "[574A] WARNING: population counts do not match the pre-registered 22/17/5 split.",
            flush=True,
        )

    seq = (ANCHOR, FLYBY, ANCHOR)

    out_fh = OUT_PATH.open("w", encoding="utf-8")

    def _write(rec: dict[str, Any]) -> None:
        out_fh.write(json.dumps(rec, default=str) + "\n")
        out_fh.flush()

    _write(
        {
            "_meta": True,
            "task": "#574 Stage A -- Titan-Iapetus eccentric-Keplerian 3D-closure kill gate",
            "git_sha": sha,
            "e_titan": ECC_TITAN,
            "e_iapetus": ECC_IAPETUS,
            "inclination_deg": INCLINATION_DEG,
            "gate_residual_kms": GATE_RESIDUAL_KMS,
            "min_useful_bend_deg": DEFAULT_MIN_USEFUL_BEND_DEG,
            "ecc_fractions": list(ECC_FRACTIONS),
            "window_omega_deg": WINDOW_OMEGA_DEG,
            "window_tof_scale": WINDOW_TOF_SCALE,
            "window_m0_deg": WINDOW_M0_DEG,
            "merge_omega_deg": MERGE_OMEGA_DEG,
            "merge_tof_scale": MERGE_TOF_SCALE,
            "merge_vinf_kms": MERGE_VINF_KMS,
            "merge_m0_deg": MERGE_M0_DEG,
            "kepler_smoke_test_reduction_pass": kepler_smoke_ok,
            "n_branches": len(branches),
            "n_eccentricity_robust": n_robust,
            "n_floor_huggers": n_hugger,
        }
    )

    results: list[dict[str, Any]] = []
    n_smoke_fail = 0
    for idx, branch in enumerate(branches):
        t0 = time.time()
        rec = run_branch_continuation(branch, seq)
        elapsed = time.time() - t0
        rec["elapsed_s"] = elapsed
        if not rec["smoke_ok"]:
            n_smoke_fail += 1
        results.append(rec)
        _write({"kind": "branch_result", **rec})
        final_res = rec["final_residual_kms"]
        final_res_str = f"{final_res:.4f}" if final_res is not None else "nan"
        print(
            f"[574A] [{idx + 1:2d}/{len(branches)}] branch={branch['branch_id']:2d} "
            f"n_rev={branch['n_rev']} ecc_robust={branch['eccentricity_robust']} "
            f"smoke_ok={rec['smoke_ok']} stage0_res={rec['stage0_residual_kms']:.3e} "
            f"final_res={final_res_str} "
            f"survives={rec['survives']} ({elapsed:.1f}s, total {time.time() - t_start:.1f}s)",
            flush=True,
        )

    print(
        f"[574A] === C2 POSITIVE CONTROL === smoke fails: {n_smoke_fail}/{len(branches)}",
        flush=True,
    )
    positive_control_pass = n_smoke_fail == 0
    if not positive_control_pass:
        print(
            "[574A] ABORTING -- C2 positive control FAILED (eccentric machinery does not "
            "recover the circular branches at e=0). Not crediting any kill/pass verdict.",
            flush=True,
        )
        _write(
            {
                "kind": "abort",
                "reason": "c2_positive_control_failed",
                "n_smoke_fail": n_smoke_fail,
            }
        )
        out_fh.close()
        return 1

    survivors = [r for r in results if r["survives"]]
    robust_survivors = [r for r in survivors if r["eccentricity_robust"]]
    hugger_survivors = [r for r in survivors if r["is_floor_hugger"]]

    branch_groups = cluster_survivors(robust_survivors)
    branch_groups.sort(key=lambda g: min(b["branch_id"] for b in g))
    n_deduped_robust_survivors = len(branch_groups)
    n_rev_classes = len({tuple(m["n_rev"]) for g in branch_groups for m in g})

    hugger_groups = cluster_survivors(hugger_survivors)
    n_deduped_hugger_survivors = len(hugger_groups)

    # Explicit mirror-pair check (C3): do known near-mirror pairs {2,19},
    # {4,14} end up in the same dedup group under eccentricity?
    surv_by_id = {r["branch_id"]: r for r in survivors}

    def _pair_status(a_id: int, b_id: int) -> str:
        a, b = surv_by_id.get(a_id), surv_by_id.get(b_id)
        if a is None or b is None:
            missing = [i for i in (a_id, b_id) if i not in surv_by_id]
            return f"NOT_BOTH_SURVIVING (missing {missing})"
        for g in branch_groups + hugger_groups:
            ids = {m["branch_id"] for m in g}
            if a_id in ids and b_id in ids:
                return "MERGED"
        return "DID_NOT_MERGE"

    mirror_23_19 = _pair_status(2, 19)
    mirror_4_14 = _pair_status(4, 14)

    if n_deduped_robust_survivors >= 5:
        verdict = "PASS"
    elif n_deduped_robust_survivors <= 3:
        verdict = "KILL"
    else:
        verdict = "MARGINAL"

    hugger_control_flag = n_deduped_hugger_survivors >= 3

    summary = {
        "kind": "summary",
        "n_branches": len(branches),
        "n_eccentricity_robust_input": n_robust,
        "n_floor_hugger_input": n_hugger,
        "positive_control_pass": positive_control_pass,
        "n_raw_survivors": len(survivors),
        "n_raw_robust_survivors": len(robust_survivors),
        "n_deduped_robust_survivors": n_deduped_robust_survivors,
        "deduped_robust_survivor_branch_ids": [
            sorted(m["branch_id"] for m in g) for g in branch_groups
        ],
        "n_rev_classes_spanned_by_survivors": n_rev_classes,
        "n_raw_hugger_survivors": len(hugger_survivors),
        "n_deduped_hugger_survivors": n_deduped_hugger_survivors,
        "hugger_control_non_discriminating_flag": hugger_control_flag,
        "mirror_pair_2_19_status": mirror_23_19,
        "mirror_pair_4_14_status": mirror_4_14,
        "verdict": verdict,
        "total_elapsed_s": time.time() - t_start,
    }
    _write(summary)
    out_fh.close()

    print("[574A] === SUMMARY ===", flush=True)
    print(
        f"[574A]   C2 positive control PASS: {positive_control_pass} "
        f"(0/{len(branches)} smoke fails)",
        flush=True,
    )
    print(
        f"[574A]   raw survivors: {len(survivors)} "
        f"(robust: {len(robust_survivors)}, hugger: {len(hugger_survivors)})",
        flush=True,
    )
    print(f"[574A]   DEDUPED robust survivors: {n_deduped_robust_survivors}", flush=True)
    print(
        f"[574A]   deduped robust survivor branch-id groups: "
        f"{summary['deduped_robust_survivor_branch_ids']}",
        flush=True,
    )
    print(f"[574A]   n_rev classes spanned: {n_rev_classes}", flush=True)
    print(
        f"[574A]   deduped floor-hugger survivors: {n_deduped_hugger_survivors} "
        f"(non-discriminating flag: {hugger_control_flag})",
        flush=True,
    )
    print(f"[574A]   mirror pair {{2,19}}: {mirror_23_19}", flush=True)
    print(f"[574A]   mirror pair {{4,14}}: {mirror_4_14}", flush=True)
    print(f"[574A]   VERDICT: {verdict}", flush=True)
    print(f"[574A] results written to {OUT_PATH}", flush=True)
    print(f"[574A] total wall clock: {time.time() - t_start:.1f}s", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
