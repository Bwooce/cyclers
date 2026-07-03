"""#535 sensitivity pre-check: does the x0=0.9920 candidate survive real Earth
eccentricity (e=0.0167, NASA fact sheet), via the existing ER3BP core
(core/er3bp.py), before committing to a full V0-V5 gauntlet pass.

Positive control first: reproduce the CR3BP result (e=0) through the ER3BP
code path (should be identical to run_535's own CR3BP propagation) before
trusting the e=0.0167 perturbed result.
"""

from __future__ import annotations

import pathlib
import sys

_REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))

import numpy as np  # noqa: E402

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
from cyclerfinder.core.er3bp import ER3BPSystem, propagate_er3bp  # noqa: E402
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import PreflightBlockedError, preflight_search  # noqa: E402
from cyclerfinder.search.hill_sphere_return_detector import (  # noqa: E402
    find_admission_windows,
    find_returns,
)

_METHOD = MethodCapability(
    genome="ER3BP re-propagation of a single already-known #535 candidate IC (core/er3bp.py) "
    "at e=0 (positive control) and e=0.0167 (real Earth eccentricity)",
    corrector="none (direct propagation + admission-window classification, no correction)",
    capability_tags=frozenset(
        {"er3bp", "ballistic", "coplanar", "single-arc", "sensitivity-check"}
    ),
    git_sha="working-tree",
)

CONTROL_X0 = 0.9920
CONTROL_XDOT0 = 0.0
CONTROL_YDOT0 = 0.13033911  # matches C=2.9839437412 at e=0 (verified in run_535)
YEARS = 50
MIN_SEPARATION_YEARS = 1.0
WINDOW_LO_YEARS, WINDOW_HI_YEARS = 10.0, 15.0
N_RETURNS_LO, N_RETURNS_HI = 3, 15
GEOMETRY_FACTOR = 3.0
EARTH_ECCENTRICITY = 0.0167  # NASA Earth fact sheet (standard sourced value)


def time_from_true_anomaly(f: np.ndarray, e: float) -> np.ndarray:
    f = np.asarray(f, dtype=float)
    n_rev = np.floor(f / (2 * np.pi))
    f_local = f - n_rev * 2 * np.pi
    e_anom = 2 * np.arctan2(
        np.sqrt(1 - e) * np.sin(f_local / 2), np.sqrt(1 + e) * np.cos(f_local / 2)
    )
    e_anom = np.mod(e_anom, 2 * np.pi)
    m_anom = e_anom - e * np.sin(e_anom)
    return n_rev * 2 * np.pi + m_anom


def evaluate(e: float, mu: float, r_hill: float, label: str) -> tuple[list, list]:
    sys_ = ER3BPSystem(mu=mu, e=e, primary_name="Sun", secondary_name="Earth")
    state0 = np.array([CONTROL_X0, 0.0, 0.0, CONTROL_XDOT0, CONTROL_YDOT0, 0.0])
    f_max = YEARS * 2.0 * np.pi
    f_eval, y, _ = propagate_er3bp(
        state0,
        (0.0, f_max),
        sys_,
        rtol=1e-11,
        atol=1e-11,
        with_stm=False,
        method="DOP853",
    )
    t = time_from_true_anomaly(f_eval, e)
    order = np.argsort(t)
    t_sorted = t[order]
    pos_rel_earth = np.stack([y[0, order] - (1.0 - mu), y[1, order]], axis=1)
    returns = find_returns(
        t_sorted, pos_rel_earth, r_hill=r_hill, min_separation=MIN_SEPARATION_YEARS * 2.0 * np.pi
    )
    windows = find_admission_windows(
        returns,
        float(t_sorted[0]),
        float(t_sorted[-1]),
        window_lo=WINDOW_LO_YEARS * 2.0 * np.pi,
        window_hi=WINDOW_HI_YEARS * 2.0 * np.pi,
        n_returns_lo=N_RETURNS_LO,
        n_returns_hi=N_RETURNS_HI,
        geometry_factor=GEOMETRY_FACTOR,
    )
    admissible = [w for w in windows if w.geometry_ok]
    print(f"[{label}] e={e}: {len(returns)} returns total, {len(admissible)} admissible window(s)")
    for r in returns:
        print(
            f"    return at t={r.t_closest / (2 * np.pi):.3f} yr, "
            f"closest={r.closest_distance / r_hill:.3f}x Hill"
        )
    return returns, admissible


def main() -> None:
    preflight_search(
        task_no=535,
        region_id="sun-earth-transient-quasi-cycler-er3bp-sensitivity-check",
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        # Exactly 2 propagations of one already-known IC (e=0 control, e=0.0167 check) --
        # not a search/sweep, but the AST gate (tests/scripts/test_scripts_call_preflight.py)
        # requires every new script call preflight_search rather than a silent exemption.
        n_points=2,
        timing_pilot_seconds_per_point=5.0,
    )

    system = cr3bp.cr3bp_system("Sun", "Earth")
    mu = float(system.mu)
    r_hill = (mu / 3.0) ** (1.0 / 3.0)

    print("=== Positive control: ER3BP with e=0 must match CR3BP result ===")
    _, admissible0 = evaluate(0.0, mu, r_hill, "control-e0")
    if not admissible0:
        raise RuntimeError(
            "POSITIVE CONTROL FAILED: ER3BP e=0 path does not reproduce the known "
            "CR3BP admissible result. Do not trust the e=0.0167 result below."
        )
    print()
    print("=== Sensitivity check: real Earth eccentricity e=0.0167 (NASA fact sheet) ===")
    _, admissible_e = evaluate(EARTH_ECCENTRICITY, mu, r_hill, "real-e")

    print()
    print("=== VERDICT ===")
    if admissible_e:
        print("SURVIVES: admission holds under real Earth eccentricity perturbation.")
    else:
        print(
            "COLLAPSES: admission does NOT survive real Earth eccentricity perturbation "
            "-- the candidate is specific to the idealized circular-restricted model."
        )


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
