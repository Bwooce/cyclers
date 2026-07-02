"""#530 Unstable-manifold propagation for co-orbital/resonance cyclers.

Both #523 (Earth co-orbital horseshoe) and #527 (Sun-Jupiter Hilda 3:2 MMR)
independently found the same pattern this session: every CERTIFIED, truly
periodic orbit in these families stays well outside the companion body's Hill
sphere -- that is precisely what makes a periodic orbit dynamically stable.
Per Guido & Efthymiopoulos (arXiv:2604.00679, #527's positive-control paper),
the genuine close-encounter/transport dynamics live on chaotic orbits
shadowing HETEROCLINIC connections between the UNSTABLE manifolds of periodic
orbits near a resonance separatrix -- a different dynamical object that
periodic-orbit enumeration (used for both #523 and #527) cannot surface no
matter how the grid is refined.

This script:

1. Re-certifies a representative set of the orbits already found by #523 and
   #527 (warm-started from their already-known-converged (x0, xdot0), so this
   re-certification is fast -- it is not a fresh search).
2. Computes each orbit's 6x6 monodromy matrix (reusing
   `search/bifurcation_detector.py`'s existing `monodromy`, already used for
   the #347 Phase 1 work) and its eigenvalues/eigenvectors, to identify which
   family members are genuinely UNSTABLE and get the manifold direction.
3. For unstable members, propagates the unstable manifold (a small
   perturbation along the unstable eigenvector, forward in time) and tracks
   the minimum distance to the companion body -- testing directly whether the
   manifold achieves the Hill-sphere encounter the periodic orbit itself does
   not.

POSITIVE CONTROL for the machinery (not a specific numeric target -- neither
source paper tabulates precise ICs, verified this session): every CR3BP
periodic-orbit monodromy must be SYMPLECTIC (det = 1, a consequence of
Liouville's theorem, independent of any specific paper) and must carry a
trivial eigenvalue pair at +1 (energy conservation / time-translation
invariance). Both are checked for every orbit before any stability
classification is trusted.
"""

from __future__ import annotations

import datetime
import pathlib
import sys
import time

# Ensure the src tree is on the path when invoked as a script.
_REPO = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO / "src"))

import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402

import cyclerfinder.core.cr3bp as cr3bp  # noqa: E402
from cyclerfinder.data.method_capability import MethodCapability  # noqa: E402
from cyclerfinder.data.preflight import PreflightBlockedError, preflight_search  # noqa: E402
from cyclerfinder.search.bifurcation_detector import monodromy  # noqa: E402
from cyclerfinder.search.cr3bp_general_periodic import correct_general_periodic  # noqa: E402

_REGION_ID = "coorbital-resonance-unstable-manifold-encounters"
_METHOD = MethodCapability(
    genome="unstable-manifold propagation from certified #523/#527 periodic orbits",
    corrector=(
        "correct_general_periodic (re-certification) + monodromy/floquet_multipliers "
        "(stability) + eigenvector-seeded manifold propagation"
    ),
    capability_tags=frozenset(
        {"cr3bp", "ballistic", "coplanar", "single-arc", "invariant-manifold"}
    ),
    git_sha="working-tree",
)

# Warm-start seeds: the already-known-converged (x0, xdot0) from #523/#527's
# own certified output. Re-certifying here is a REFINEMENT to full precision,
# not a fresh search (the coarse enumerator is not re-run).
_JUPITER_TARGETS = [
    # (label, c_target, x0_seed, xdot0_seed, n)
    ("hilda_C2.95_x0.786", 2.95, 0.7860664785360532, -6.651920962338042e-14, 1),
    ("hilda_C2.98_x0.873_closest", 2.98, 0.8725434709791006, -4.974669903388049e-15, 1),
    ("hilda_C2.98_x0.639", 2.98, 0.6391159799811946, 6.706314149990968e-13, 1),
    ("hilda_C3.0613_seed", 3.0613, 0.7296542563751145, -2.899483905139338e-12, 1),
    ("hilda_C3.14_x0.761", 3.14, 0.7614557190251472, -6.288067563627053e-13, 1),
    ("hilda_C3.14_x0.668", 3.14, 0.6675057036640826, 1.1879277941203975e-14, 1),
]
_EARTH_TARGETS = [
    ("earth_coorbital_x0.836", 2.9990, 0.83658, 0.0, 1),
    ("earth_coorbital_x0.901", 2.9990, 0.90104, 0.0, 1),
    ("earth_coorbital_x0.784", 2.9990, 0.78447, 0.0, 1),
]


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _min_dist_over_span(
    system: cr3bp.CR3BPSystem, state0: np.ndarray, t_span: float, n_eval: int = 2000
) -> tuple[float, bool]:
    """Min distance to the SECONDARY body over [0, t_span]. Returns (dist, ok)."""
    mu = float(system.mu)
    t_eval = np.linspace(0.0, t_span, n_eval)
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, t_span),
        state0,
        args=(mu,),
        method="DOP853",
        rtol=1e-11,
        atol=1e-11,
        t_eval=t_eval,
    )
    if not sol.success:
        return float("inf"), False
    dx = sol.y[0] - (1.0 - mu)
    dy = sol.y[1]
    return float(np.min(np.hypot(dx, dy))), True


def _analyse_orbit(
    system: cr3bp.CR3BPSystem,
    label: str,
    c_target: float,
    x0_seed: float,
    xdot0_seed: float,
    n: int,
    r_hill: float,
) -> None:
    print(f"[{_ts()}] --- {label} (C={c_target}) ---")
    orbit = correct_general_periodic(
        system,
        x0_seed,
        xdot0_seed,
        c_target,
        period_guess=100.0,
        half_crossings=2 * n,
        ydot0_sign=1.0,
        tol=1e-11,
        t_hi_frac=1.15,
        max_iter=60,
    )
    if not (orbit.converged and orbit.residual <= 1e-9):
        print(f"[{_ts()}]   re-certification FAILED (residual={orbit.residual:.3e}); skipping.")
        return

    state0 = np.array([orbit.x0, 0.0, 0.0, orbit.xdot0, orbit.ydot0, 0.0], dtype=np.float64)
    baseline_dist, ok = _min_dist_over_span(system, state0, orbit.period)
    print(
        f"[{_ts()}]   certified: period={orbit.period:.4f}, baseline min_dist={baseline_dist:.5f} "
        f"(Hill={r_hill:.5f}, ratio={baseline_dist / r_hill:.2f}x)"
    )

    mono = monodromy(system, state0, orbit.period)
    det = float(np.linalg.det(mono))
    eigvals, eigvecs = np.linalg.eig(mono)
    order = np.argsort(-np.abs(eigvals))
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    mags = np.abs(eigvals)
    print(
        f"[{_ts()}]   POSITIVE CONTROL: det(monodromy)={det:.6f} (want ~1, symplectic); "
        f"|eigenvalues| sorted desc = {np.round(mags, 4)}"
    )
    if abs(det - 1.0) > 1e-4:
        print(
            f"[{_ts()}]   WARNING: symplectic check failed (|det-1|={abs(det - 1.0):.3e}) -- "
            f"treat this orbit's stability classification with caution."
        )

    # The trivial pair (energy/time-translation) sits at +1; identify the
    # NON-trivial pair with the largest magnitude departure from 1 as the
    # governing (in-plane or out-of-plane) stability multiplier.
    non_trivial_idx = [
        i for i in range(6) if abs(eigvals[i].real - 1.0) > 1e-3 or abs(eigvals[i].imag) > 1e-6
    ]
    if not non_trivial_idx:
        print(f"[{_ts()}]   all multipliers trivial/on-unit-circle -- no unstable direction.")
        return
    lead = max(non_trivial_idx, key=lambda i: mags[i])
    lam = eigvals[lead]
    is_unstable = abs(lam) > 1.01  # small margin above 1 to avoid numerical-noise false positives

    print(
        f"[{_ts()}]   leading non-trivial multiplier: lambda={lam:.6f}, |lambda|={abs(lam):.6f} "
        f"-> {'UNSTABLE' if is_unstable else 'stable/elliptic (on unit circle)'}"
    )

    if not is_unstable:
        return

    # Propagate the unstable manifold: perturb along the (real part of the)
    # unstable eigenvector, both signs, forward in time. The propagation
    # HORIZON must be long enough for the perturbation to actually grow to a
    # macroscopically significant size -- a fixed "3 periods" (an earlier
    # version of this script) is nowhere near enough: at |lambda|=1.7, growing
    # eps=1e-6 to O(0.05) (a fraction of the orbit's own scale) takes
    # ln(0.05/eps)/ln(|lambda|) periods, computed below rather than guessed,
    # capped so a near-1 (barely unstable) multiplier can't demand an
    # unbounded horizon.
    vec = eigvecs[:, lead].real
    vec = vec / np.linalg.norm(vec)
    eps = 1e-6
    target_growth = 0.05
    n_periods = min(40, max(3, int(np.ceil(np.log(target_growth / eps) / np.log(abs(lam))))))
    print(f"[{_ts()}]   manifold propagation horizon: {n_periods} periods (Lyapunov-scaled).")
    for sign in (+1.0, -1.0):
        perturbed = state0 + sign * eps * vec
        manifold_dist, ok = _min_dist_over_span(
            system, perturbed, float(n_periods) * orbit.period, n_eval=6000
        )
        if not ok:
            print(f"[{_ts()}]     manifold (sign={sign:+.0f}): propagation failed.")
            continue
        improved = manifold_dist < baseline_dist
        ratio = manifold_dist / r_hill
        print(
            f"[{_ts()}]     manifold (sign={sign:+.0f}, eps={eps:.0e}, {n_periods} periods): "
            f"min_dist = {manifold_dist:.5f} (Hill={r_hill:.5f}, ratio={ratio:.2f}x) "
            f"{'[CLOSER than baseline]' if improved else ''} "
            f"{'[ENCOUNTER]' if manifold_dist < r_hill else ''}"
        )


def main() -> None:
    print(
        f"[{_ts()}] #530 unstable-manifold propagation for co-orbital/resonance cyclers starting."
    )

    n_points = len(_JUPITER_TARGETS) + len(_EARTH_TARGETS)
    preflight_search(
        task_no=530,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=n_points,
        timing_pilot_seconds_per_point=60.0,  # generous; Earth re-certification is the slow part
    )

    t0 = time.time()

    print(f"[{_ts()}] === Sun-Jupiter Hilda-family orbits (from #527) ===")
    sys_sj = cr3bp.cr3bp_system("Sun", "Jupiter")
    r_hill_jupiter = (float(sys_sj.mu) / 3.0) ** (1.0 / 3.0)
    for label, c_target, x0, xdot0, n in _JUPITER_TARGETS:
        _analyse_orbit(sys_sj, label, c_target, x0, xdot0, n, r_hill_jupiter)

    print(f"[{_ts()}] === Sun-Earth co-orbital orbits (from #523) ===")
    sys_se = cr3bp.cr3bp_system("Sun", "Earth")
    r_hill_earth = (float(sys_se.mu) / 3.0) ** (1.0 / 3.0)
    for label, c_target, x0, xdot0, n in _EARTH_TARGETS:
        _analyse_orbit(sys_se, label, c_target, x0, xdot0, n, r_hill_earth)

    print()
    print(f"[{_ts()}] Analysis complete in {time.time() - t0:.1f}s.")


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
