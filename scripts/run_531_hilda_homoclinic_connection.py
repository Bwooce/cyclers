"""#531 Homoclinic-connection test for #530's unstable Hilda orbit.

#530 found exactly one genuinely unstable orbit among 9 re-certified #523/#527
periodic orbits (Sun-Jupiter Hilda, C=3.14, x0~=0.7615, lambda~=1.718), and
showed its manifold's LOCAL/LINEAR-REGIME growth (a single perturbed
trajectory, Lyapunov-scaled propagation horizon) does not bring it within
Earth's -- sorry, Jupiter's -- Hill sphere. This script tests the next,
cheaper step #530's writeup flagged: does this orbit have a genuine
HOMOCLINIC CONNECTION (its own unstable manifold Wu meeting its own stable
manifold Ws transversally on the {y=0} section, a real invariant-manifold
object, not just a short local propagation), and does the CONNECTED path
(the two manifold legs joined at the matched crossing) pass through the Hill
sphere anywhere along its length?

Reuses the #314 heteroclinic_cycle framework UNCHANGED (correct_connection,
assemble_cycle, crosscheck_cycle, _planar_floquet_pair, _seed_on_manifold) --
this script only supplies a new LyapunovNode built directly from #530's
already-computed orbit (bypassing the libration-specific from_libration
constructor, which LyapunovNode does not require) and explores the
(k_u, k_s, branch_u, branch_s) manifold-branch/crossing-index space for THIS
orbit's geometry, since the framework's defaults are tuned for the W-Z
Oterma L1<->L2 case, not a Hilda resonant orbit's homoclinic tangle.

POSITIVE CONTROL: re-runs the framework's own existing golden test live
(test_connection_l1_to_l2_converges, the W-Z Sun-Jupiter-Oterma L1->L2
heteroclinic connection, arXiv:math/0201278) before trusting anything from
the exploration below -- confirms the machinery works TODAY, not just that
it worked when #314 was built.
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
from cyclerfinder.genome.heteroclinic_cycle import (  # noqa: E402
    HeteroclinicConnection,
    LyapunovNode,
    _planar_floquet_pair,
    _seed_on_manifold,
    assemble_cycle,
    correct_connection,
    crosscheck_cycle,
)
from cyclerfinder.search.cr3bp_general_periodic import correct_general_periodic  # noqa: E402

_REGION_ID = "hilda-c3.14-homoclinic-connection-hill-encounter"
_METHOD = MethodCapability(
    genome="#314 heteroclinic_cycle LyapunovNode built from a #530 unstable Hilda orbit",
    corrector="correct_connection / assemble_cycle / crosscheck_cycle (#314, unchanged)",
    capability_tags=frozenset(
        {"cr3bp", "ballistic", "coplanar", "single-arc", "invariant-manifold", "heteroclinic"}
    ),
    git_sha="working-tree",
)

# The #530 unstable Hilda orbit's seed (re-certified fresh here for
# full-precision state0/period/jacobi -- not read from a log).
HILDA_C = 3.14
HILDA_X0_SEED = 0.7614557190251472
HILDA_XDOT0_SEED = -6.288067563627053e-13

# Known-good (k_u, k_s, branch_u, branch_s), found by the first full 36-combo
# exploration pass this session (384s) and independently cross-checked
# (DOP853 vs Radau, residual 3.582e-08 both): (1, 1, -1, +1). Set to None to
# force a fresh full exploration instead of reusing this result.
KNOWN_GOOD_COMBO: tuple[int, int, int, int] | None = (1, 1, -1, +1)


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


def _positive_control() -> None:
    """Live re-run of the #314 framework's own W-Z Oterma golden (not cached)."""
    print(
        f"[{_ts()}] POSITIVE CONTROL: re-running the #314 framework's W-Z Oterma "
        f"L1->L2 heteroclinic-connection golden live ..."
    )
    wz_mu = 0.0009537
    wz_c_ours = 3.03 - wz_mu * (1.0 - wz_mu)
    system = cr3bp.CR3BPSystem(
        mu=wz_mu, primary="sun", secondary="jupiter", l_km=778.57e6, t_s=5.957e8
    )
    l1 = LyapunovNode.from_libration(
        system, x0_guess=0.9208034913207400196, jacobi=wz_c_ours, period_guess=3.0, label="L1"
    )
    l2 = LyapunovNode.from_libration(
        system,
        x0_guess=1.081929486841799903,
        jacobi=wz_c_ours,
        period_guess=3.0,
        label="L2",
        ydot0_sign=-1.0,
    )
    conn = correct_connection(system, l1, l2, tol=1e-7)
    if not conn.converged:
        raise RuntimeError(
            f"POSITIVE CONTROL FAILED: the #314 framework's own W-Z Oterma L1->L2 "
            f"golden did not converge live (residual={conn.residual:.3e}). Do not "
            f"trust anything below -- the machinery itself is not working today."
        )
    print(
        f"[{_ts()}] POSITIVE CONTROL PASSED: W-Z L1->L2 connection converged, "
        f"residual={conn.residual:.3e}, crossing x={conn.crossing_xv[0]:.5f} "
        f"(W-Z Part I: ~0.95792)."
    )


def _build_hilda_node(system: cr3bp.CR3BPSystem) -> LyapunovNode:
    orbit = correct_general_periodic(
        system,
        HILDA_X0_SEED,
        HILDA_XDOT0_SEED,
        HILDA_C,
        period_guess=100.0,
        half_crossings=2,
        ydot0_sign=1.0,
        tol=1e-11,
        t_hi_frac=1.15,
        max_iter=60,
    )
    if not (orbit.converged and orbit.residual <= 1e-9):
        raise RuntimeError(
            f"Hilda C=3.14 orbit re-certification failed: residual={orbit.residual:.3e}"
        )
    state0 = np.array([orbit.x0, 0.0, 0.0, orbit.xdot0, orbit.ydot0, 0.0], dtype=np.float64)
    lam_u, v_u, lam_s, v_s = _planar_floquet_pair(system, state0, orbit.period)
    print(
        f"[{_ts()}] Hilda node built: period={orbit.period:.4f}, jacobi={HILDA_C}, "
        f"lambda_u={lam_u:.4f}, lambda_s={lam_s:.4f}"
    )
    return LyapunovNode(
        label="hilda_C3.14",
        state0=state0,
        period=orbit.period,
        jacobi=HILDA_C,
        unstable_eigvec=v_u,
        stable_eigvec=v_s,
        converged=True,
    )


def _min_dist_over_span(
    system: cr3bp.CR3BPSystem, state0: np.ndarray, t_span: float, n_eval: int = 3000
) -> float:
    """Min distance to the secondary over [0, t_span] (t_span may be negative).

    t_eval must be monotonic in the SAME direction as the integration span
    (0.0 -> t_span), not artificially re-sorted ascending -- scipy requires
    t_eval to match the span direction (an earlier version got this backwards
    for negative t_span and crashed with 'not properly sorted').
    """
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
        return float("inf")
    dx = sol.y[0] - (1.0 - mu)
    dy = sol.y[1]
    return float(np.min(np.hypot(dx, dy)))


def main() -> None:
    print(f"[{_ts()}] #531 homoclinic-connection test for the #530 unstable Hilda orbit starting.")

    system = cr3bp.cr3bp_system("Sun", "Jupiter")
    r_hill = (float(system.mu) / 3.0) ** (1.0 / 3.0)

    # n_points: the (k_u, k_s, branch) exploration grid below.
    k_values = (1, 2, 3)
    branch_values = (+1, -1)
    n_points = len(k_values) * len(k_values) * len(branch_values) * len(branch_values)

    preflight_search(
        task_no=531,
        region_id=_REGION_ID,
        method=_METHOD,
        script_path=pathlib.Path(__file__),
        n_points=n_points,
        # Measured this session: correct_connection's internal scan+Newton at
        # reduced exploration settings (scan_n=8, max_iter=15) is a few seconds
        # per (k_u, k_s, branch) combination.
        timing_pilot_seconds_per_point=5.0,
    )

    _positive_control()

    node = _build_hilda_node(system)

    if KNOWN_GOOD_COMBO is not None:
        k_u, k_s, branch_u, branch_s = KNOWN_GOOD_COMBO
        print(
            f"[{_ts()}] Reusing known-good combo from a prior full exploration this "
            f"session (skipping the {n_points}-combination search): k_u={k_u} k_s={k_s} "
            f"branch_u={branch_u:+d} branch_s={branch_s:+d}."
        )
    else:
        print(
            f"[{_ts()}] Exploring (k_u, k_s, branch_u, branch_s) for the homoclinic "
            f"self-connection ({n_points} combinations, reduced-budget pass) ..."
        )
        t0 = time.time()
        best: tuple[float, HeteroclinicConnection, tuple[int, int, int, int]] | None = None
        n_done = 0
        for k_u in k_values:
            for k_s in k_values:
                for branch_u in branch_values:
                    for branch_s in branch_values:
                        n_done += 1
                        conn = correct_connection(
                            system,
                            node,
                            node,
                            k_u=k_u,
                            k_s=k_s,
                            branch_u=branch_u,
                            branch_s=branch_s,
                            scan_n=8,
                            max_iter=15,
                            tol=1e-5,
                        )
                        print(
                            f"[{_ts()}]   [{n_done}/{n_points}] k_u={k_u} k_s={k_s} "
                            f"branch_u={branch_u:+d} branch_s={branch_s:+d}: "
                            f"residual={conn.residual:.3e}"
                        )
                        if best is None or conn.residual < best[0]:
                            best = (conn.residual, conn, (k_u, k_s, branch_u, branch_s))
        print(f"[{_ts()}] Exploration complete in {time.time() - t0:.1f}s.")

        if best is None or not np.isfinite(best[0]):
            print(f"[{_ts()}] No (k_u, k_s, branch) combination reached the section at all.")
            return

        best_residual, _best_conn, (k_u, k_s, branch_u, branch_s) = best
        print(
            f"[{_ts()}] Best exploration candidate: k_u={k_u} k_s={k_s} branch_u={branch_u:+d} "
            f"branch_s={branch_s:+d}, residual={best_residual:.3e}"
        )

    print(f"[{_ts()}] Refining the best candidate to full precision (tol=1e-7) ...")
    refined = correct_connection(
        system,
        node,
        node,
        k_u=k_u,
        k_s=k_s,
        branch_u=branch_u,
        branch_s=branch_s,
        scan_n=20,
        max_iter=60,
        tol=1e-7,
    )
    print(
        f"[{_ts()}]   refined: converged={refined.converged}, residual={refined.residual:.3e}, "
        f"crossing=({refined.crossing_xv[0]:.5f}, {refined.crossing_xv[1]:.5f})"
    )

    if not refined.converged:
        print(
            f"[{_ts()}] No genuine homoclinic connection certified for this orbit at this "
            f"(k_u, k_s, branch) choice. Not registering an encounter claim."
        )
        return

    cycle = assemble_cycle(
        system,
        [node],
        tol=1e-7,
        per_leg_kwargs=[{"k_u": k_u, "k_s": k_s, "branch_u": branch_u, "branch_s": branch_s}],
    )
    print(
        f"[{_ts()}] assemble_cycle: closed={cycle.closed}, "
        f"max_leg_residual={cycle.max_leg_residual:.3e}"
    )

    checked = crosscheck_cycle(system, [node], cycle)
    print(
        f"[{_ts()}] crosscheck_cycle (independent Radau): independent_residual="
        f"{checked.independent_residual:.3e}"
    )

    # Track Hill-sphere distance along BOTH manifold legs of the certified
    # connection -- the actual test #530 could not do (it only propagated a
    # short local segment, not a full manifold leg to a certified crossing).
    max_time = 8.0 * node.period
    seed_u = _seed_on_manifold(
        system, node, tau=refined.tau_u, direction="unstable", branch=branch_u, epsilon=1e-6
    )
    seed_s = _seed_on_manifold(
        system, node, tau=refined.tau_s, direction="stable", branch=branch_s, epsilon=1e-6
    )
    dist_unstable_leg = _min_dist_over_span(system, seed_u, max_time)
    dist_stable_leg = _min_dist_over_span(system, seed_s, -max_time)
    print(
        f"[{_ts()}] Unstable leg (forward, {max_time:.1f} nondim): min_dist_to_jupiter="
        f"{dist_unstable_leg:.5f} (Hill={r_hill:.5f}, ratio={dist_unstable_leg / r_hill:.2f}x) "
        f"{'[ENCOUNTER]' if dist_unstable_leg < r_hill else ''}"
    )
    print(
        f"[{_ts()}] Stable leg (backward, {max_time:.1f} nondim): min_dist_to_jupiter="
        f"{dist_stable_leg:.5f} (Hill={r_hill:.5f}, ratio={dist_stable_leg / r_hill:.2f}x) "
        f"{'[ENCOUNTER]' if dist_stable_leg < r_hill else ''}"
    )

    print()
    print(f"[{_ts()}] #531 complete.")


if __name__ == "__main__":
    try:
        main()
    except PreflightBlockedError as exc:
        print(f"[{_ts()}] BLOCKED by preflight_search:\n{exc}")
        sys.exit(1)
