"""#572 -- Titan-Iapetus 3D-closure probe (throwaway, one-off; NOT #552).

Decision gate for whether the #552 3D/inclined-releg genome extension is
worth scoping (see ``data/OUTSTANDING.md`` #572, which is authoritative --
this docstring only summarizes). #571's coplanar Titan-Iapetus sweep found
187 gate-passing candidates, all under the coplanar idealization (Iapetus's
real ~15.5 deg inclination to Titan's orbital plane is ignored by the
`_moon_state()` circular-coplanar generator every #558-lineage script uses).
An Opus + Fable bend-gate re-check found 149/187 survive a *necessary-not-
sufficient* max-bend check under the real inclination correction, but that
is not a closure proof. This script attempts an actual 3D Lambert closure
for the two independently Fable-verified best Titan-anchored candidates
(see the #572 OUTSTANDING.md entry's "Fable correction" for why the
original best-by-residual exemplar, rel_offset=18/tof_scale=1.15, is WRONG
to use -- it's dead under its own corrected-bend gate).

Method
------
Titan is kept in the primary's xy reference plane (z=0), matching how
#571's own sweep (and every #558-lineage script) already treats it -- this
project's stated approach is "Titan stays equatorial, Iapetus goes 3D", not
a full two-inclined-body treatment.

Iapetus is placed on a REAL 3D circular orbit inclined by ``INCLINATION_DEG``
(15.5 deg, the documented conservative estimate cited in the #572 entry;
``core/satellites.py`` carries no inclination field) to Titan's plane, with
an explicit ascending-node longitude ``Omega`` (measured from Titan's fixed
t=0 position, exploiting the same rotational symmetry the coplanar sweep
itself relies on -- WLOG Titan's t=0 phase is 0). Iapetus's in-plane
argument of latitude at t=0 is set from the candidate's own ``rel_offset_deg``
(preserving the coplanar seed's phase pattern); ``Omega`` (node alignment) is
the genuinely NEW free variable this probe searches over, per the #572
Fable correction that a negative result without this search is
formulation-conditional, not a real closure failure.

Standard orbital-plane rotation (R3(Omega) . R1(inc) applied to a circular
in-plane state parameterized by argument of latitude ``u``):

    x = a*(cosO*cosu - sinO*sinu*cosi)
    y = a*(sinO*cosu + cosO*sinu*cosi)
    z = a*sinu*sini
    vx = -n*a*(cosO*sinu + sinO*cosu*cosi)
    vy =  n*a*(-sinO*sinu + cosO*cosu*cosi)
    vz =  n*a*cosu*sini

At ``inc=0`` this reduces EXACTLY to ``discovery_campaign._moon_state``'s
coplanar formula for any Omega (Omega+u plays the role of theta) -- checked
below as an inline smoke test before the real search, so a bug in the
rotation algebra cannot silently masquerade as "no 3D closure".

``core/lambert.py`` is fed genuinely 3D ``r1, r2`` (already confirmed 3D-
capable by two independent reviews per the OUTSTANDING.md entry: 3D dot
product for the transfer angle, ``cross_z``-based branch selection is
exactly the z-component of the true 3D ``r1 x r2``, 3D f/g velocity
construction) -- no coplanar assumption is introduced in the solve itself,
only in the (now-replaced) upstream state generator.

Residual and gate are IDENTICAL in form to #558/#571's own
``residual_at_point`` / ``candidate_passes_physical_gate`` (the latter is
imported and reused verbatim, per task scope -- do not reimplement it).

Decision criterion ("near the coplanar residual"): the project-wide
GATE_RESIDUAL_KMS = 0.05 km/s bar (the same bar #558/#571's own gate uses
everywhere) is used as "near" -- both coplanar seeds (3.3e-3, 6.2e-3 km/s)
sit well inside it, so a 3D solution clearing this SAME bar is a fair
apples-to-apples closure claim, not a loosened one. The realized minimum
residual is reported in full either way so the reader can judge "near" for
themselves.

Run as::

    uv run python scripts/probe_572_titan_iapetus_3d_closure.py
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from scan_558_uranus_all_pairs_offset_sweep import GATE_RESIDUAL_KMS  # noqa: E402

from cyclerfinder.core.lambert import (  # noqa: E402
    LambertConvergenceError,
    LambertGeometryError,
    lambert,
)
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.discovery_campaign import (  # noqa: E402
    DAY_S,
    _mean_motion_rad_day,
    _moon_state,
)
from cyclerfinder.search.physical_sanity import (  # noqa: E402
    DEFAULT_MIN_USEFUL_BEND_DEG,
    candidate_passes_physical_gate,
)

PRIMARY = "Saturn"
ANCHOR = "Titan"
FLYBY = "Iapetus"
INCLINATION_DEG = 15.5  # documented conservative estimate; core/satellites.py has no field

DATA_DIR = ROOT / "data"

# The 2 Fable-verified Titan-anchored candidates from
# data/scan_571_saturn_titan_iapetus.jsonl (confirmed by direct grep against
# the raw file before writing this script -- see task transcript). DO NOT
# substitute the "best-by-residual" record (rel_offset=18, tof_scale=1.15):
# that one is dead under its own corrected-bend gate (see OUTSTANDING.md
# #572 Fable correction).
CANDIDATES: tuple[dict[str, Any], ...] = (
    {
        "label": "cand1_rel255_tof1.80_n11",
        "rel_offset_deg": 255.00000000000003,
        "tof_scale": 1.8,
        "n_rev": (1, 1),
        "coplanar_residual_kms": 0.003346388794962518,
        "coplanar_iapetus_vinf_kms": 0.9815144709583703,
        "coplanar_bend_deg": 8.86,
    },
    {
        "label": "cand2_rel89_tof0.70_n00",
        "rel_offset_deg": 89.00000000000001,
        "tof_scale": 0.7,
        "n_rev": (0, 0),
        "coplanar_residual_kms": 0.006203942031542287,
        "coplanar_iapetus_vinf_kms": 0.9858131311412763,
        "coplanar_bend_deg": 8.83,
    },
)


def _git_sha() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True, cwd=ROOT
        ).strip()
    except Exception:
        return "unknown"


def iapetus_state_3d(
    u_rad: float, v_circ_km_s: float, sma_km: float, omega_rad: float, inc_rad: float
) -> tuple[np.ndarray, np.ndarray]:
    """Iapetus circular state on a plane inclined ``inc_rad`` to the primary's
    xy plane, ascending node at longitude ``omega_rad`` (measured from the
    same fixed axis Titan's t=0 phase uses), argument of latitude ``u_rad``.

    ``v_circ_km_s`` is the circular orbital SPEED (``sqrt(mu/sma_km)``, km/s
    -- matching ``_moon_state``'s own ``v_circ``, NOT the rad/day mean
    motion; ``du/dt`` in rad/s is ``v_circ_km_s / sma_km``, so scaling the
    velocity components by ``v_circ_km_s`` directly keeps everything in
    km/s without a day/second unit-mismatch).

    Reduces EXACTLY to ``_moon_state``'s coplanar formula at ``inc_rad=0``
    (verified in ``_smoke_test_reduction`` below).
    """
    cos_o, sin_o = math.cos(omega_rad), math.sin(omega_rad)
    cosi, sini = math.cos(inc_rad), math.sin(inc_rad)
    cosu, sinu = math.cos(u_rad), math.sin(u_rad)
    x = sma_km * (cos_o * cosu - sin_o * sinu * cosi)
    y = sma_km * (sin_o * cosu + cos_o * sinu * cosi)
    z = sma_km * sinu * sini
    vx = -v_circ_km_s * (cos_o * sinu + sin_o * cosu * cosi)
    vy = v_circ_km_s * (-sin_o * sinu + cos_o * cosu * cosi)
    vz = v_circ_km_s * cosu * sini
    pos = np.array([x, y, z])
    vel = np.array([vx, vy, vz])
    return pos, vel


def _smoke_test_reduction() -> bool:
    """At inc=0, iapetus_state_3d must reproduce _moon_state exactly (any Omega)."""
    mu = PRIMARIES[PRIMARY]
    sat_b = SATELLITES[FLYBY]
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    v_circ = math.sqrt(mu / sat_b.sma_km)
    ok = True
    for omega_deg in (0.0, 37.0, 190.0, 355.0):
        for u_deg in (0.0, 42.0, 200.0):
            u = math.radians(u_deg)
            omega = math.radians(omega_deg)
            r3d, v3d = iapetus_state_3d(u, v_circ, sat_b.sma_km, omega, 0.0)
            r2d, v2d = _moon_state(u + omega, n_b, 0.0, sat_b.sma_km, mu)
            dr = float(np.linalg.norm(r3d - r2d))
            dv = float(np.linalg.norm(v3d - v2d))
            if dr > 1e-6 or dv > 1e-9:
                print(
                    f"  SMOKE TEST FAIL: omega={omega_deg} u={u_deg} dr={dr:.3e} dv={dv:.3e}",
                    flush=True,
                )
                ok = False
    return ok


def _leg_best(
    r_a: np.ndarray,
    v_a: np.ndarray,
    r_b: np.ndarray,
    v_b: np.ndarray,
    tof_s: float,
    mu: float,
    n_rev: int,
    *,
    retry_perturb_rad: float = 0.0,
) -> dict[str, Any] | None:
    """Solve one leg for the exact requested n_rev; return vinf_out/vinf_in/v1/v2.

    Returns None if n_rev infeasible (Lambert did not admit that revolution
    count). Raises LambertGeometryError to the caller (handled by the node-
    alignment sweep loop, which perturbs Omega and retries -- see #572's
    "do not count LambertGeometryError as no-closure" note).
    """
    sols = lambert(r_a, r_b, tof_s, mu=mu, max_revs=max(0, n_rev))
    cands = [s for s in sols if s.n_revs == n_rev]
    if not cands:
        return None
    best = min(cands, key=lambda s: float(np.linalg.norm(s.v1 - v_a)))
    return {
        "vinf_out": float(np.linalg.norm(best.v1 - v_a)),
        "vinf_in": float(np.linalg.norm(best.v2 - v_b)),
    }


def evaluate_point(
    rel_offset_deg: float,
    tof_scale: float,
    n_rev: tuple[int, int],
    omega_deg: float,
    inc_deg: float,
) -> dict[str, Any] | None:
    """One (rel_offset, tof_scale, n_rev, Omega, inc) 3D evaluation.

    Returns None on LambertGeometryError/ConvergenceError at this exact
    Omega (caller retries at a nearby Omega -- these are solver-domain
    artifacts, not "no closure" per #572's explicit instruction).
    """
    mu = PRIMARIES[PRIMARY]
    sat_a = SATELLITES[ANCHOR]
    sat_b = SATELLITES[FLYBY]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    v_circ_b = math.sqrt(mu / sat_b.sma_km)
    p_a = 2.0 * math.pi / n_a
    p_b = 2.0 * math.pi / n_b
    tof = tof_scale * math.sqrt(p_a * p_b)
    tof_s = tof * DAY_S

    omega = math.radians(omega_deg)
    inc = math.radians(inc_deg)
    u0 = math.radians(rel_offset_deg)  # Iapetus argument-of-latitude at t=0

    r0, v0 = _moon_state(0.0, n_a, 0.0, sat_a.sma_km, mu)  # Titan @ t=0, planar
    r1, v1 = iapetus_state_3d(u0 + n_b * tof, v_circ_b, sat_b.sma_km, omega, inc)  # Iapetus @ tof
    r2, v2 = _moon_state(0.0, n_a, 2.0 * tof, sat_a.sma_km, mu)  # Titan @ 2*tof, planar

    n0, n1 = n_rev
    try:
        leg0 = _leg_best(r0, v0, r1, v1, tof_s, mu, n0)
        leg1 = _leg_best(r1, v1, r2, v2, tof_s, mu, n1)
    except (LambertGeometryError, LambertConvergenceError):
        return None
    if leg0 is None or leg1 is None:
        return None

    r_mid = abs(leg0["vinf_in"] - leg1["vinf_out"])
    r_periodic = abs(leg0["vinf_out"] - leg1["vinf_in"])
    residual = max(r_mid, r_periodic)

    vinf0 = leg0["vinf_out"]  # Titan, departure
    vinf1 = max(leg0["vinf_in"], leg1["vinf_out"])  # Iapetus, in/out
    vinf2 = leg1["vinf_in"]  # Titan, arrival

    return {
        "rel_offset_deg": rel_offset_deg,
        "tof_scale": tof_scale,
        "n_rev": list(n_rev),
        "omega_deg": omega_deg,
        "inc_deg": inc_deg,
        "residual_kms": residual,
        "vinf_kms": [vinf0, vinf1, vinf2],
    }


def sweep_node_alignment(cand: dict[str, Any], *, n_omega: int = 3600) -> dict[str, Any]:
    """Fine grid-search of Omega in [0, 360) deg at fixed (rel_offset, tof_scale,
    n_rev), handling LambertGeometryError via a small retry perturbation, then
    enumerate ALL local-minimum BASINS in the residual-vs-Omega landscape and
    locally refine EACH ONE (bounded Nelder-Mead, window restricted around the
    basin) rather than only the single globally-lowest-residual point.

    This matters: a first pass (see task transcript) found the GLOBAL
    residual minimum sits on a physically different, high-V_inf Lambert
    branch (~8.6 km/s vs the coplanar seed's ~1.0-1.4 km/s) that fails the
    #324 bend gate outright -- an unconstrained single-seed Nelder-Mead
    started from that global minimum wanders further into that same wrong
    branch and would silently miss a DIFFERENT, lower-residual-than-gate,
    gate-PASSING basin at a different Omega. Per #572's own node-alignment
    discipline (a negative result without searching node alignment is
    formulation-conditional), every local minimum is checked, not just the
    deepest one -- 1 basin surviving is enough for closure_found=True.
    """
    from scipy.optimize import minimize

    rel_offset_deg = cand["rel_offset_deg"]
    tof_scale = cand["tof_scale"]
    n_rev = tuple(cand["n_rev"])

    grid: list[tuple[float, dict[str, Any] | None]] = []
    n_geometry_errors = 0
    n_geometry_retried_ok = 0
    for i in range(n_omega):
        omega_deg = 360.0 * i / n_omega
        pt = evaluate_point(rel_offset_deg, tof_scale, n_rev, omega_deg, INCLINATION_DEG)
        if pt is None:
            n_geometry_errors += 1
            # Retry at small perturbations near the singular point (#572
            # note 2: LambertGeometryError is a solver-domain artifact near
            # the 180-deg-transfer singularity, not evidence of no closure).
            for eps in (0.02, -0.02, 0.05, -0.05, 0.15, -0.15):
                pt = evaluate_point(
                    rel_offset_deg, tof_scale, n_rev, omega_deg + eps, INCLINATION_DEG
                )
                if pt is not None:
                    n_geometry_retried_ok += 1
                    break
        grid.append((omega_deg, pt))

    n_feasible = sum(1 for _o, p in grid if p is not None)

    # Local-minimum (basin) detection on the feasible subset of the periodic grid.
    basins_grid: list[dict[str, Any]] = []
    for i in range(n_omega):
        _o, p = grid[i]
        if p is None:
            continue
        _op, pp = grid[i - 1]
        _on, pn = grid[(i + 1) % n_omega]
        if pp is None or pn is None:
            continue
        if p["residual_kms"] <= pp["residual_kms"] and p["residual_kms"] <= pn["residual_kms"]:
            basins_grid.append(p)
    basins_grid.sort(key=lambda r: r["residual_kms"])

    # Refine each basin (bounded Nelder-Mead, window restricted to +-15 deg in
    # Omega and +-0.1 in tof_scale around that basin's own grid point --
    # prevents the optimizer from wandering into a DIFFERENT branch, which is
    # exactly what an unbounded single-seed search did on the first pass).
    refined_basins: list[dict[str, Any]] = []
    for seed in basins_grid[:12]:  # cap: enough to cover every distinct branch, bounded cost
        omega_lo, omega_hi = seed["omega_deg"] - 15.0, seed["omega_deg"] + 15.0
        tof_lo, tof_hi = tof_scale - 0.1, tof_scale + 0.1

        def _obj(x: np.ndarray, _n_rev: tuple[int, int] = n_rev) -> float:
            omega_deg_x, tof_scale_x = float(x[0]), float(x[1])
            pt = evaluate_point(rel_offset_deg, tof_scale_x, _n_rev, omega_deg_x, INCLINATION_DEG)
            if pt is None:
                return 1.0e3
            return pt["residual_kms"]

        x0 = np.array([seed["omega_deg"], tof_scale])
        res = minimize(
            _obj,
            x0,
            method="Nelder-Mead",
            bounds=[(omega_lo, omega_hi), (tof_lo, tof_hi)],
            options={"xatol": 1e-5, "fatol": 1e-9, "maxiter": 200, "maxfev": 200},
        )
        refined_pt = evaluate_point(
            rel_offset_deg, float(res.x[1]), n_rev, float(res.x[0]) % 360.0, INCLINATION_DEG
        )
        if refined_pt is not None and refined_pt["residual_kms"] <= seed["residual_kms"]:
            refined_basins.append(refined_pt)
        else:
            refined_basins.append(seed)

    refined_basins.sort(key=lambda r: r["residual_kms"])

    return {
        "label": cand["label"],
        "rel_offset_deg": rel_offset_deg,
        "tof_scale_seed": tof_scale,
        "n_rev": list(n_rev),
        "coplanar_residual_kms": cand["coplanar_residual_kms"],
        "coplanar_iapetus_vinf_kms": cand["coplanar_iapetus_vinf_kms"],
        "coplanar_bend_deg": cand["coplanar_bend_deg"],
        "inclination_deg": INCLINATION_DEG,
        "n_omega_grid": n_omega,
        "n_geometry_errors_encountered": n_geometry_errors,
        "n_geometry_errors_resolved_by_retry": n_geometry_retried_ok,
        "n_feasible_omega_points": n_feasible,
        "n_basins_found": len(basins_grid),
        "basins": refined_basins,
        "best_overall": refined_basins[0] if refined_basins else None,
    }


def main() -> int:
    sha = _git_sha()
    print(f"[572] Titan-Iapetus 3D-closure probe -- sha={sha}", flush=True)

    print("[572] smoke test: iapetus_state_3d reduces to _moon_state at inc=0 ...", flush=True)
    smoke_ok = _smoke_test_reduction()
    print(f"[572] smoke test PASS: {smoke_ok}", flush=True)
    if not smoke_ok:
        print("[572] ABORTING -- 3D state generator does not reduce correctly.", flush=True)
        return 1

    out_records: list[dict[str, Any]] = [
        {
            "_meta": True,
            "task": "#572 Titan-Iapetus 3D-closure probe",
            "git_sha": sha,
            "inclination_deg": INCLINATION_DEG,
            "gate_residual_kms": GATE_RESIDUAL_KMS,
            "min_useful_bend_deg": DEFAULT_MIN_USEFUL_BEND_DEG,
            "smoke_test_reduction_pass": smoke_ok,
        }
    ]

    verdicts: list[dict[str, Any]] = []
    seq = (ANCHOR, FLYBY, ANCHOR)
    for cand in CANDIDATES:
        print(f"[572] --- {cand['label']} ---", flush=True)
        t0 = time.time()
        sweep = sweep_node_alignment(cand, n_omega=3600)
        elapsed = time.time() - t0
        sweep["elapsed_s"] = elapsed
        print(
            f"[572]   n_feasible_omega={sweep['n_feasible_omega_points']}/{sweep['n_omega_grid']}  "
            f"n_basins={sweep['n_basins_found']}  "
            f"geom_errors={sweep['n_geometry_errors_encountered']} "
            f"(resolved_by_retry={sweep['n_geometry_errors_resolved_by_retry']})  "
            f"({elapsed:.1f}s)",
            flush=True,
        )

        # Evaluate EVERY refined basin against BOTH gates (residual-near-
        # coplanar AND #324 physical bend) -- per the diagnostic run in the
        # task transcript, the single globally-lowest-residual basin can be a
        # physically different (high-V_inf, low-bend) branch that FAILS the
        # bend gate while a higher-(but-still-under-gate)-residual basin at a
        # different Omega PASSES both gates. A closure is "found" if ANY
        # basin clears both, not only the deepest one.
        basin_evals: list[dict[str, Any]] = []
        for b in sweep["basins"]:
            residual_near = b["residual_kms"] < GATE_RESIDUAL_KMS
            gate_pass, gate_verdicts = candidate_passes_physical_gate(
                seq, tuple(b["vinf_kms"]), min_useful_bend_deg=DEFAULT_MIN_USEFUL_BEND_DEG
            )
            bends = [v.max_bend_deg for v in gate_verdicts]
            closure = bool(residual_near and gate_pass)
            basin_evals.append(
                {
                    "omega_deg": b["omega_deg"],
                    "tof_scale": b["tof_scale"],
                    "residual_kms": b["residual_kms"],
                    "vinf_kms": b["vinf_kms"],
                    "residual_near_coplanar": residual_near,
                    "physical_gate_pass": gate_pass,
                    "max_bend_deg_per_encounter": bends,
                    "closure": closure,
                }
            )
            vinf_str = [f"{v:.3f}" for v in b["vinf_kms"]]
            bends_str = [f"{x:.2f}" for x in bends]
            print(
                f"[572]     basin omega={b['omega_deg']:7.3f} tof_scale={b['tof_scale']:.4f} "
                f"residual={b['residual_kms']:.6f}  vinf={vinf_str}  "
                f"bends={bends_str}  "
                f"near_gate={residual_near} phys_gate={gate_pass} CLOSURE={closure}",
                flush=True,
            )

        closing_basins = [be for be in basin_evals if be["closure"]]
        closure = len(closing_basins) > 0
        if not sweep["basins"]:
            print("[572]   NO feasible 3D point found at ANY node alignment tried.", flush=True)
            verdict = {
                "label": cand["label"],
                "closure_found": False,
                "reason": "no_feasible_lambert_point_at_any_omega",
            }
        else:
            best_closing = (
                min(closing_basins, key=lambda b: b["residual_kms"]) if closing_basins else None
            )
            print(f"[572]   >>> CLOSURE VERDICT for {cand['label']}: {closure} <<<", flush=True)
            verdict = {
                "label": cand["label"],
                "closure_found": closure,
                "n_basins_evaluated": len(basin_evals),
                "n_closing_basins": len(closing_basins),
                "best_closing_basin": best_closing,
                "all_basins": basin_evals,
            }
        verdicts.append(verdict)
        out_records.append({"kind": "candidate_result", **sweep})
        out_records.append({"kind": "candidate_verdict", **verdict})

    out_path = DATA_DIR / "probe_572_titan_iapetus_3d_closure.jsonl"
    with out_path.open("w", encoding="utf-8") as fh:
        for rec in out_records:
            fh.write(json.dumps(rec, default=str) + "\n")
    print(f"[572] results written to {out_path}", flush=True)

    print("[572] === SUMMARY ===", flush=True)
    for v in verdicts:
        print(f"[572]   {v['label']}: closure_found={v['closure_found']}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
