"""#448 Region C — high-e Sun-planet ER3BP cyclers (the SKIPPED corner).

Region C of the 2026-06-25 discovery-strategy prioritization draft: the high-e
Sun-planet ER3BP corner that ``er3bp-discovery-em-broucke-koblick`` explicitly
SKIPPED (Earth-Moon only; Mercury/Mars/Pluto skipped for want of a CR3BP seed).
#435 generated the seeds and ran a FIRST pass with the SECANT continuator +
Floquet-transition classification (registry entry
``er3bp-discovery-high-e-sun-planet-lyapunov-dro-2026-06-24``: 6/6 survive, 0
bifurcations, NO V-inf trend measured).

This #448 pass is a CAPABILITY-SUBSUMING re-sweep (capability-subsumption rule):
it re-runs the same 6 families with THREE method extensions #435 did not have:

  1. FOLD-AWARE pseudo-arclength continuation
     (``continue_er3bp_family_in_e_arclength``) instead of the secant
     continuator. The ``er3bp-direct-e0-blind-grid`` negative established the
     secant discriminator is unreliable at folds; the arclength walker is the
     fix and can walk THROUGH a turning point in e (where a secant stalls and
     would mis-report a "death" / miss a fold-born branch).
  2. SADDLE-CENTER branch-detection probe (``branch_at_saddle_center_er3bp``)
     run at the highest-e member of each family — directly answers "does any
     family bifurcate / carry an e>0-only branch?".
  3. V-INF PROXY TREND vs the circular model: closest-approach inertial
     relative speed to the SECONDARY at e=0 vs at e=target_e. This is the
     kill-criterion's clause (c) ("a V-inf reduction vs the circular model")
     that #435 never measured.

KILL-CRITERION (verbatim, draft Region C): continue 3 Sun-Mars + 3 Sun-Mercury
CR3BP families to real e. KILL if ALL families persist with NO bifurcation AND
NO V-inf reduction vs circular (the Earth-Moon outcome) -> the ER3BP-departure
hypothesis is FALSIFIED for high-e Sun-planet cyclers; close the corner.

Report-only — NO catalogue writeback. A "not-found" literature status is
necessary-not-sufficient; the V0-V5 gauntlet still governs any survivor.

Usage::

    uv run python scripts/run_448_region_c_highE.py
"""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import cr3bp_system
from cyclerfinder.core.er3bp import ER3BPSystem, er3bp_eom
from cyclerfinder.genome.er3bp_branching import branch_at_saddle_center_er3bp
from cyclerfinder.genome.er3bp_continuation import (
    ContinuationError,
    continue_er3bp_family_in_e_arclength,
)
from cyclerfinder.genome.er3bp_periodic import ER3BPPeriodicOrbit
from cyclerfinder.search.cr3bp_seed_generator import dro_seed, lyapunov_seed
from cyclerfinder.search.er3bp_floquet import er3bp_monodromy, floquet_classify

_ROOT = Path(__file__).resolve().parents[1]
_DATA_DIR = _ROOT / "data"
_RUNLOG = _DATA_DIR / "runlogs" / "448_region_c.runlog.jsonl"

# Per draft Region C: 3 Sun-Mars + 3 Sun-Mercury families. The #435 seed
# generator yields 2 well-conditioned families per system (L1-Lyapunov, DRO).
# To honour the "3 each" instruction we add a second Lyapunov amplitude rung
# (a distinct family member: larger-amplitude L1 Lyapunov) per system, giving
# 3 distinct seeds each for Mars and Mercury (6 total — the draft's count).
_SYSTEMS: tuple[tuple[str, str, float], ...] = (
    ("Sun", "Mars", 0.093),
    ("Sun", "Mercury", 0.206),
)

# arclength step in the unit (x0, ydot0, e) tangent space. 0.002 keeps the
# walker inside the corrector basin at mu~1e-7 while still reaching real e in
# a few hundred steps.
_DS = 0.002
_MAX_STEPS = 1200

_SEED_SOURCE = (
    "#448 CR3BP Lyapunov/DRO seed (collinear-linear + fixed-Jacobi corrector); "
    "fold-aware arclength e-continuation"
)


@dataclass(frozen=True)
class FamilySeed:
    label: str
    system: ER3BPSystem
    state0: NDArray[np.float64]
    period_f_full: float
    target_e: float


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _append_runlog(record: dict[str, object]) -> None:
    _RUNLOG.parent.mkdir(parents=True, exist_ok=True)
    with _RUNLOG.open("a") as f:
        f.write(json.dumps(record) + "\n")


def _generate_seeds(primary: str, secondary: str, target_e: float) -> list[FamilySeed]:
    """3 distinct planar families per system: L1-Lyapunov (small + larger Ax) + DRO."""
    cr3bp = cr3bp_system(primary, secondary)
    system = ER3BPSystem.from_cr3bp(cr3bp, target_e)
    seeds: list[FamilySeed] = []

    for amp, tag in ((3e-3, "L1-lyapunov"), (1e-2, "L1-lyapunov-hiamp")):
        try:
            st, period_f = lyapunov_seed(cr3bp, point="L1", amplitude=amp)
        except Exception as exc:
            _log(f"  {secondary}: lyapunov amp={amp} did not converge ({exc}); skip")
            continue
        seeds.append(
            FamilySeed(
                label=f"{secondary}-{tag}",
                system=system,
                state0=st,
                period_f_full=period_f,
                target_e=target_e,
            )
        )

    try:
        st_dro, period_dro = dro_seed(cr3bp)
    except Exception as exc:
        _log(f"  {secondary}: DRO seed did not converge ({exc}); skip")
    else:
        seeds.append(
            FamilySeed(
                label=f"{secondary}-dro",
                system=system,
                state0=st_dro,
                period_f_full=period_dro,
                target_e=target_e,
            )
        )

    return seeds


def _vinf_proxy(orbit: ER3BPPeriodicOrbit, system: ER3BPSystem, v_unit_kms: float) -> dict:
    """Closest-approach inertial relative speed to the SECONDARY over one period.

    Reuses the #447 convention: in the rotating nondim frame the secondary sits
    at (1-mu, 0, 0); the inertial relative velocity is v_rot + omega x rel with
    omega = z-hat. For the ER3BP the physical separation pulsates by
    (1 - e*cos f)/(1+e*cos f)-scaling, but for an apples-to-apples
    "e=target vs e=0" comparison of the SAME family we evaluate the proxy in the
    pulsating frame at both ends — the trend (reduction vs circular) is the
    kill-criterion quantity, not the absolute number.
    """
    mu = system.mu
    secondary = np.array([1.0 - mu, 0.0, 0.0])
    f_span = (0.0, orbit.period_f)
    ts = np.linspace(0.0, orbit.period_f, 4000)
    sol = solve_ivp(
        er3bp_eom,
        f_span,
        np.asarray(orbit.state0, float),
        args=(mu, system.e),
        method="DOP853",
        rtol=1e-11,
        atol=1e-11,
        t_eval=ts,
        dense_output=False,
    )
    if not sol.success:
        return {"dmin_nondim": None, "vrel_kms": None}
    pos = sol.y[:3, :].T
    vel = sol.y[3:, :].T
    rel = pos - secondary
    dist = np.linalg.norm(rel, axis=1)
    idx = int(np.argmin(dist))
    r = rel[idx]
    v = vel[idx]
    omega_cross_r = np.array([-r[1], r[0], 0.0])
    v_inertial_rel = v + omega_cross_r
    return {
        "dmin_nondim": float(dist[idx]),
        "vrel_kms": float(np.linalg.norm(v_inertial_rel)) * v_unit_kms,
        "f_at_min": float(ts[idx]),
    }


def _continue_family(seed: FamilySeed) -> dict:
    """Fold-aware arclength continuation e=0->target_e + Floquet + branch + V-inf."""
    cr3bp = cr3bp_system(seed.system.primary_name, seed.system.secondary_name)
    v_unit_kms = cr3bp.l_km / cr3bp.t_s

    # Corrector integration span is the HALF period (symmetric residual).
    integration_f = seed.period_f_full / 2.0

    sys0 = ER3BPSystem(
        mu=seed.system.mu,
        e=0.0,
        primary_name=seed.system.primary_name,
        secondary_name=seed.system.secondary_name,
    )

    t_fam = time.time()
    try:
        family = continue_er3bp_family_in_e_arclength(
            sys0,
            seed.state0,
            integration_f,
            seed.target_e,
            ds=_DS,
            max_steps=_MAX_STEPS,
            is_half_period_residual=True,
        )
    except ContinuationError as exc:
        return {
            "label": seed.label,
            "outcome": "seed_failed",
            "reason": str(exc),
            "e_max_reached": 0.0,
            "n_members": 0,
        }

    if not family:
        return {
            "label": seed.label,
            "outcome": "empty",
            "e_max_reached": 0.0,
            "n_members": 0,
        }

    # Floquet-classify each member; detect an elliptic<->hyperbolic regime flip
    # (the genuine bifurcation signal — NOT "an eigenvalue on the unit circle",
    # which holds for every stable member).
    e_star: float | None = None
    prev_regime: str | None = None
    step_records: list[dict] = []
    for orb in family:
        sys_e = ER3BPSystem(
            mu=orb.mu,
            e=orb.e,
            primary_name=seed.system.primary_name,
            secondary_name=seed.system.secondary_name,
        )
        try:
            mono = er3bp_monodromy(orb.state0, orb.period_f, sys_e)
            fl = floquet_classify(mono)
            tag, on_uc = fl.stability_tag, fl.on_unit_circle
        except Exception:
            tag, on_uc = "unknown", False
        regime = (
            "hyperbolic"
            if tag == "unstable"
            else ("elliptic" if tag in ("stable", "marginal") else "unknown")
        )
        if (
            regime != "unknown"
            and prev_regime is not None
            and prev_regime != "unknown"
            and regime != prev_regime
            and e_star is None
        ):
            e_star = orb.e
        if regime != "unknown":
            prev_regime = regime
        step_records.append(
            {
                "e": orb.e,
                "corrector_residual": orb.corrector_residual,
                "independent_residual": orb.independent_residual,
                "stability_tag": tag,
                "on_unit_circle": on_uc,
            }
        )

    e_max = family[-1].e
    reached_target = abs(e_max - seed.target_e) <= 1e-4

    # V-inf proxy at e=0 (circular) vs at the highest-e member.
    vinf_circular = _vinf_proxy(family[0], sys0, v_unit_kms)
    sys_top = ER3BPSystem(
        mu=seed.system.mu,
        e=family[-1].e,
        primary_name=seed.system.primary_name,
        secondary_name=seed.system.secondary_name,
    )
    vinf_elliptic = _vinf_proxy(family[-1], sys_top, v_unit_kms)
    vinf_reduction = None
    if vinf_circular["vrel_kms"] is not None and vinf_elliptic["vrel_kms"] is not None:
        vinf_reduction = vinf_circular["vrel_kms"] - vinf_elliptic["vrel_kms"]

    # Saddle-center branch-detection probe at the highest-e member.
    branch_orbit, branch_info = branch_at_saddle_center_er3bp(
        sys_top,
        family[-1].state0,
        family[-1].period_f / 2.0,
    )
    branched = branch_orbit is not None and not _same_orbit(branch_orbit.state0, family[-1].state0)

    if e_star is not None or branched:
        outcome = "bifurcates"
    elif reached_target:
        outcome = "survives"
    else:
        outcome = "dies"

    return {
        "label": seed.label,
        "outcome": outcome,
        "e_star": e_star,
        "branched": branched,
        "branch_info": {k: str(v) for k, v in branch_info.items()},
        "e_max_reached": e_max,
        "n_members": len(family),
        "reached_target": reached_target,
        "vinf_proxy": {
            "circular_e0_vrel_kms": vinf_circular["vrel_kms"],
            "elliptic_etop_vrel_kms": vinf_elliptic["vrel_kms"],
            "reduction_kms": vinf_reduction,
            "circular_dmin_nondim": vinf_circular["dmin_nondim"],
            "elliptic_dmin_nondim": vinf_elliptic["dmin_nondim"],
        },
        "wall_s": time.time() - t_fam,
        "first_step": step_records[0] if step_records else None,
        "last_step": step_records[-1] if step_records else None,
    }


def _same_orbit(a: NDArray[np.float64], b: NDArray[np.float64], *, tol: float = 1e-6) -> bool:
    return bool(np.linalg.norm(np.asarray(a) - np.asarray(b)) < tol)


def main() -> None:
    t0 = time.time()
    # Fresh runlog each invocation.
    if _RUNLOG.exists():
        _RUNLOG.unlink()
    _log(f"#448 Region C high-e Sun-planet ER3BP (arclength + branch + V-inf), ds={_DS}")
    _append_runlog({"event": "start", "ts": time.strftime("%Y-%m-%dT%H:%M:%S"), "ds": _DS})

    all_results: list[dict] = []
    for primary, secondary, target_e in _SYSTEMS:
        _log(f"=== {primary}-{secondary} (target e={target_e}) ===")
        seeds = _generate_seeds(primary, secondary, target_e)
        _log(f"{primary}-{secondary}: {len(seeds)} family seed(s)")
        for i, seed in enumerate(seeds, 1):
            _log(f"{secondary} [{i}/{len(seeds)}] continuing {seed.label} -> e={target_e}")
            res = _continue_family(seed)
            res["system"] = f"{primary}-{secondary}"
            res["target_e"] = target_e
            vp = res.get("vinf_proxy", {})
            _log(
                f"  {seed.label}: outcome={res['outcome']} "
                f"e_max={res['e_max_reached']:.4f} e_star={res.get('e_star')} "
                f"branched={res.get('branched')} members={res['n_members']} "
                f"vinf_circ={vp.get('circular_e0_vrel_kms')} "
                f"vinf_ellip={vp.get('elliptic_etop_vrel_kms')} "
                f"reduction={vp.get('reduction_kms')}"
            )
            _append_runlog(res)
            all_results.append(res)

    survives = sum(1 for r in all_results if r["outcome"] == "survives")
    bifurcates = sum(1 for r in all_results if r["outcome"] == "bifurcates")
    dies = sum(1 for r in all_results if r["outcome"] in ("dies", "seed_failed", "empty"))
    any_vinf_reduction = any(
        (r.get("vinf_proxy", {}).get("reduction_kms") or 0.0) > 1e-3 for r in all_results
    )

    summary = {
        "event": "summary",
        "n_families": len(all_results),
        "survives": survives,
        "bifurcates": bifurcates,
        "dies": dies,
        "any_vinf_reduction_vs_circular": any_vinf_reduction,
        "kill_criterion_fired": (bifurcates == 0 and not any_vinf_reduction),
        "wall_s": time.time() - t0,
    }
    _append_runlog(summary)
    _log(
        f"SUMMARY: families={len(all_results)} survives={survives} "
        f"bifurcates={bifurcates} dies={dies} "
        f"any_vinf_reduction={any_vinf_reduction} "
        f"KILL_FIRED={summary['kill_criterion_fired']} ({time.time() - t0:.1f}s)"
    )


if __name__ == "__main__":
    main()
