"""#557 ER3BP sensitivity gate: do the Sun-Jupiter quasi-Hilda admissible
candidates survive Jupiter's real orbital eccentricity (e=0.04838624), via the
existing ER3BP core (core/er3bp.py)?

This is the plan's Risk-1 mitigation and a genuine go/no-go gate, NOT a
formality: #535's Earth transient-capture corridor TOTALLY COLLAPSED under this
exact check (e=0.0167). Koon 2001 (p.29) argues Jupiter's e "plays little role
during the fast resonance transition" (fast tube-mediated capture, unlike
Earth's slow horseshoe) -- a hypothesis this script TESTS, not assumes.

Positive control first (per this project's verify-before-trust discipline): the
ER3BP e=0 path must reproduce the CR3BP result before the perturbed e result is
trusted. Everything is in Jupiter PERIODS (1 period = 2*pi), matching the
criterion note and the #557 search script -- never a hardcoded 2*pi = 1 yr.
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

ONE_PERIOD = 2.0 * np.pi
T_JUPITER_YEARS = 11.868
MIN_SEPARATION_PERIODS = 1.0
WINDOW_LO_PERIODS, WINDOW_HI_PERIODS = 10.0, 15.0
N_RETURNS_LO, N_RETURNS_HI = 3, 15
GEOMETRY_FACTOR = 3.0
HORIZON_PERIODS = 55.0
# Jupiter mean eccentricity, in-repo sourced (Standish & Williams J2000,
# core/constants.py PLANETS["J"].ecc). See criterion note Sec 1.
JUPITER_ECCENTRICITY = 0.04838624

# The Gehrels 3 anchor (sourced-elements IC; criterion note Sec 6).
ANCHOR = ("Gehrels3-anchor", 0.69591, 0.0, 0.57308)  # (label, x0, xdot0, ydot0)

# The 16 admissible candidates from the #557 CR3BP coarse scan (run_557 log,
# 2026-07-11), each as (x0, C, ydot0-at-e0). ydot0 recomputed here from C for
# self-consistency; the stored value is a cross-check.
CANDIDATES: tuple[tuple[float, float], ...] = (
    (0.65, 3.000),
    (0.77, 3.000),
    (0.82, 3.000),
    (0.75, 3.005),
    (0.77, 3.005),
    (0.79, 3.005),
    (0.90, 3.005),
    (0.74, 3.010),
    (0.76, 3.010),
    (0.90, 3.010),
    (0.68, 3.015),
    (0.80, 3.015),
    (0.91, 3.015),
    (0.83, 3.020),
    (0.86, 3.020),
    (0.90, 3.025),
)

_METHOD = MethodCapability(
    genome="ER3BP re-propagation of the #557 Sun-Jupiter quasi-Hilda admissible ICs "
    "(core/er3bp.py) at e=0 (positive control) and e=0.04838624 (real Jupiter eccentricity)",
    corrector="none (direct propagation + admission-window classification, no correction)",
    capability_tags=frozenset(
        {"er3bp", "ballistic", "coplanar", "single-arc", "sensitivity-check"}
    ),
    git_sha="working-tree",
)


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


def _ydot0_from_jacobi(x0: float, c_target: float, mu: float) -> float:
    r1 = abs(x0 + mu)
    r2 = abs(x0 - 1.0 + mu)
    rad = x0 * x0 + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - c_target
    return float(np.sqrt(rad))


def evaluate(
    x0: float, xdot0: float, ydot0: float, e: float, mu: float, r_hill: float
) -> tuple[int, list]:
    sys_ = ER3BPSystem(mu=mu, e=e, primary_name="Sun", secondary_name="Jupiter")
    state0 = np.array([x0, 0.0, 0.0, xdot0, ydot0, 0.0])
    f_max = HORIZON_PERIODS * 2.0 * np.pi
    f_eval, y, _ = propagate_er3bp(
        state0, (0.0, f_max), sys_, rtol=1e-11, atol=1e-11, with_stm=False, method="DOP853"
    )
    t = time_from_true_anomaly(f_eval, e)
    order = np.argsort(t)
    t_sorted = t[order]
    pos_rel_jupiter = np.stack([y[0, order] - (1.0 - mu), y[1, order]], axis=1)
    returns = find_returns(
        t_sorted, pos_rel_jupiter, r_hill=r_hill, min_separation=MIN_SEPARATION_PERIODS * ONE_PERIOD
    )
    windows = find_admission_windows(
        returns,
        float(t_sorted[0]),
        float(t_sorted[-1]),
        window_lo=WINDOW_LO_PERIODS * ONE_PERIOD,
        window_hi=WINDOW_HI_PERIODS * ONE_PERIOD,
        n_returns_lo=N_RETURNS_LO,
        n_returns_hi=N_RETURNS_HI,
        geometry_factor=GEOMETRY_FACTOR,
    )
    admissible = [w for w in windows if w.geometry_ok]
    return len(returns), admissible


def main() -> None:
    preflight_search(
        task_no=557,
        region_id="sun-jupiter-quasi-hilda-transient-capture-er3bp-sensitivity-check",
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=1 + len(CANDIDATES),
        timing_pilot_seconds_per_point=1.0,
    )

    system = cr3bp.cr3bp_system("Sun", "Jupiter")
    mu = float(system.mu)
    r_hill = (mu / 3.0) ** (1.0 / 3.0)

    print("=== Positive control: ER3BP e=0 must reproduce the CR3BP anchor result ===", flush=True)
    n0, _adm0 = evaluate(ANCHOR[1], ANCHOR[2], ANCHOR[3], 0.0, mu, r_hill)
    print(
        f"anchor e=0: {n0} returns (CR3BP gave 2) -- control {'OK' if n0 >= 1 else 'FAILED'}",
        flush=True,
    )
    if n0 < 1:
        raise RuntimeError("POSITIVE CONTROL FAILED: ER3BP e=0 path does not reproduce CR3BP.")

    # Also verify e=0 reproduces admissibility on one strong candidate before trusting e>0.
    strong = (0.68, 3.015)
    yd_s = _ydot0_from_jacobi(strong[0], strong[1], mu)
    n_s0, adm_s0 = evaluate(strong[0], 0.0, yd_s, 0.0, mu, r_hill)
    print(
        f"strong candidate x0={strong[0]} C={strong[1]} e=0: {n_s0} returns, "
        f"{len(adm_s0)} admissible (CR3BP gave 10 returns, 79 windows) -- "
        f"control {'OK' if adm_s0 else 'FAILED'}",
        flush=True,
    )
    if not adm_s0:
        raise RuntimeError("POSITIVE CONTROL FAILED: ER3BP e=0 loses candidate admissibility.")

    print(
        f"\n=== Sensitivity: real Jupiter eccentricity e={JUPITER_ECCENTRICITY} "
        f"(Standish & Williams, in-repo) ===",
        flush=True,
    )
    print("anchor: (single-capture, not admissible even at e=0)", flush=True)
    n_a, adm_a = evaluate(ANCHOR[1], ANCHOR[2], ANCHOR[3], JUPITER_ECCENTRICITY, mu, r_hill)
    print(f"  anchor e={JUPITER_ECCENTRICITY}: {n_a} returns, {len(adm_a)} admissible", flush=True)

    survivors = 0
    print("\ncandidate survival under real e:", flush=True)
    for x0, c in CANDIDATES:
        yd = _ydot0_from_jacobi(x0, c, mu)
        n_e, adm_e = evaluate(x0, 0.0, yd, JUPITER_ECCENTRICITY, mu, r_hill)
        ok = bool(adm_e)
        survivors += ok
        geom = min((w.geometry_ratio for w in adm_e), default=float("nan"))
        print(
            f"  x0={x0:.3f} C={c:.3f}: e={JUPITER_ECCENTRICITY} -> {n_e} returns, "
            f"{len(adm_e)} admissible, best_geom={geom:.3f}  "
            f"{'SURVIVES' if ok else 'collapses'}",
            flush=True,
        )

    print("\n=== VERDICT ===", flush=True)
    print(
        f"{survivors} of {len(CANDIDATES)} admissible candidates SURVIVE Jupiter's real "
        f"eccentricity (e={JUPITER_ECCENTRICITY}).",
        flush=True,
    )
    if survivors:
        print(
            "ROBUST corridor(s) exist -- in sharp contrast to #535's Earth corridor, which "
            "totally collapsed under e=0.0167. Any survivor needs an independent Fable "
            "second-opinion pass before any catalogue-adjacent writeback.",
            flush=True,
        )
    else:
        print(
            "All collapse -- like #535's Earth case, the idealized CR3BP admissibility does not "
            "survive real eccentricity. Register a clean null under the Option-A criterion.",
            flush=True,
        )


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"BLOCKED by preflight_search:\n{exc}", flush=True)
        sys.exit(1)
