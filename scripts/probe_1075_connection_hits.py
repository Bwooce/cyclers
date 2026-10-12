"""#1075 probe: do the stored manifold-connection known hits still converge today?

Cases (the arguments are copied from the tests):
  jovian  - tests/search/test_jovian_resonant_connections.py _KNOWN_HIT_T3 via the
            scan_n=8 grid seed (test_table3_gate_honest_result), stored crossing
            x=-1.3068706902782006, xdot=0.02869075783271665.
  tertiary - tests/search/test_saturn_titan_resonant_connections.py _TERTIARY with its
            stored tau, stored crossing x=-1.12647328, xdot=0.13211509.
For each case it prints the run with the test's own settings, then variants
(fd_step, max_iter) so a robustness change can be judged by identity with the stored point.

Usage: uv run python scripts/probe_1075_connection_hits.py {jovian,tertiary} [--variants]
"""

from __future__ import annotations

import argparse
import time

import numpy as np

import cyclerfinder.search.jovian_resonant_connections as jrc
import cyclerfinder.search.saturn_titan_resonant_connections as stc
import cyclerfinder.search.saturn_titan_resonant_families as stf
from cyclerfinder.genome.heteroclinic_cycle import correct_connection


def _show(tag: str, conn: object, t0: float, ref: tuple[float, float]) -> None:
    xv = np.asarray(conn.crossing_xv)  # type: ignore[attr-defined]
    d = float(np.hypot(xv[0] - ref[0], xv[1] - ref[1]))
    print(
        f"{tag}: conv={conn.converged} res={conn.residual:.3e} it={conn.n_iter} "  # type: ignore[attr-defined]
        f"tau=({conn.tau_u!r},{conn.tau_s!r}) xv=({xv[0]!r},{xv[1]!r}) "  # type: ignore[attr-defined]
        f"dist_to_stored={d:.3e} notes={conn.notes!r} [{time.time() - t0:.0f}s]",  # type: ignore[attr-defined]
        flush=True,
    )


def jovian(variants: bool) -> None:
    system = jrc.jupiter_europa_system()
    _s, n34 = jrc.build_3_4_lo_node(system)
    _s, n56 = jrc.build_5_6_lo_node(system)
    ref = (-1.3068706902782006, 0.02869075783271665)
    base = dict(
        epsilon=jrc.ANDERSON_LO_EPSILON,
        branch_u=+1,
        branch_s=+1,
        k_u=4,
        k_s=4,
        ydot_sign_u=jrc.SECTION_YDOT_SIGN,
        ydot_sign_s=jrc.SECTION_YDOT_SIGN,
        x_sign_u=jrc.SECTION_X_SIGN,
        x_sign_s=jrc.SECTION_X_SIGN,
        max_time_factor=3.0,
        scan_n=8,
        tol=1e-7,
        max_iter=20,
    )
    t0 = time.time()
    c = correct_connection(system, n34, n56, **base)  # type: ignore[arg-type]
    _show("jovian test-settings", c, t0, ref)
    if variants:
        for kw in ({"max_iter": 60}, {"fd_step": 1e-7}, {"fd_step": 1e-8}, {"scan_n": 12}):
            t0 = time.time()
            c = correct_connection(system, n34, n56, **{**base, **kw})  # type: ignore[arg-type]
            _show(f"jovian {kw}", c, t0, ref)


def tertiary(variants: bool) -> None:
    system = stf.saturn_titan_system()
    _s, node, _r = stc.build_34_node(system)
    ref = (-1.12647328, 0.13211509)
    base = dict(
        k_u=4,
        k_s=4,
        epsilon=stc.EPSILON,
        branch_u=+1,
        branch_s=-1,
        tau_u0=3.570427157637045,
        tau_s0=23.179487310480308,
        ydot_sign_u=stc.SECTION_YDOT_SIGN,
        ydot_sign_s=stc.SECTION_YDOT_SIGN,
        x_sign_u=None,
        x_sign_s=None,
        max_time_factor=3.0,
        tol=1e-9,
        max_iter=60,
        fd_step=1e-7,
    )
    t0 = time.time()
    c = correct_connection(system, node, node, **base)  # type: ignore[arg-type]
    _show("tertiary test-settings", c, t0, ref)
    if variants:
        for kw in ({"fd_step": 1e-8}, {"fd_step": 1e-6}, {"max_iter": 200}, {"tol": 1e-8}):
            t0 = time.time()
            c = correct_connection(system, node, node, **{**base, **kw})  # type: ignore[arg-type]
            _show(f"tertiary {kw}", c, t0, ref)


def _floor(system: object, a: object, b: object, tau_u: float, tau_s: float, kw: dict) -> None:  # type: ignore[type-arg]
    """Sample R(tau) on a 13-point line through (tau_u, tau_s) along each phase at steps
    1e-11, 1e-10, 1e-9; fit a straight line per component; print the scatter (std of the
    fit residuals). The scatter is the noise of the section-gap map (integration error
    amplified by the manifold's instability); Newton cannot drive ||R|| reliably below it."""
    from cyclerfinder.genome.heteroclinic_cycle import _connection_residual

    for axis in ("u", "s"):
        for h in (1e-11, 1e-10, 1e-9):
            ks = np.arange(-6, 7)
            rs = []
            for k in ks:
                du, ds = (k * h, 0.0) if axis == "u" else (0.0, k * h)
                r, _ = _connection_residual(system, a, b, tau_u=tau_u + du, tau_s=tau_s + ds, **kw)  # type: ignore[arg-type]
                assert r is not None
                rs.append(np.asarray(r))
            arr = np.array(rs)
            scatter = []
            for c in range(2):
                p = np.polyfit(ks * h, arr[:, c], 1)
                scatter.append(float(np.std(arr[:, c] - np.polyval(p, ks * h))))
            print(
                f"axis={axis} h={h:g} slope~{np.polyfit(ks * h, arr[:, 0], 1)[0]:.3e} "
                f"scatter(x,xdot)=({scatter[0]:.3e},{scatter[1]:.3e}) "
                f"|R(0)|={np.linalg.norm(arr[6]):.3e}",
                flush=True,
            )


def tertiary_floor() -> None:
    system = stf.saturn_titan_system()
    _s, node, _r = stc.build_34_node(system)
    kw = dict(
        k_u=4,
        k_s=4,
        epsilon=stc.EPSILON,
        branch_u=+1,
        branch_s=-1,
        max_time=3.0 * node.period,
        ydot_sign_u=stc.SECTION_YDOT_SIGN,
        ydot_sign_s=stc.SECTION_YDOT_SIGN,
        x_sign_u=None,
        x_sign_s=None,
    )
    # today's converged phases (probe log, test settings)
    _floor(system, node, node, 3.5704271625549313, 23.17948731233368, kw)


def symmetric_floor() -> None:
    """Noise floor at the #782 near-perpendicular (k_u = k_s = 4, branch -1/-1) hit behind
    test 9 (find_symmetric_chain_seed), at the phases today's Newton stalls at (saturn_scans)."""
    system = stf.saturn_titan_system()
    _s, node, _r = stc.build_34_node(system)
    kw = dict(
        k_u=4,
        k_s=4,
        epsilon=stc.EPSILON,
        branch_u=-1,
        branch_s=-1,
        max_time=3.0 * node.period,
        ydot_sign_u=stc.SECTION_YDOT_SIGN,
        ydot_sign_s=stc.SECTION_YDOT_SIGN,
        x_sign_u=None,
        x_sign_s=None,
    )
    _floor(system, node, node, 3.188855090798739, 22.951944857171487, kw)


def jovian_floor() -> None:
    system = jrc.jupiter_europa_system()
    _s, n34 = jrc.build_3_4_lo_node(system)
    _s, n56 = jrc.build_5_6_lo_node(system)
    kw = dict(
        k_u=4,
        k_s=4,
        epsilon=jrc.ANDERSON_LO_EPSILON,
        branch_u=+1,
        branch_s=+1,
        max_time=3.0 * max(n34.period, n56.period),
        ydot_sign_u=jrc.SECTION_YDOT_SIGN,
        ydot_sign_s=jrc.SECTION_YDOT_SIGN,
        x_sign_u=jrc.SECTION_X_SIGN,
        x_sign_s=jrc.SECTION_X_SIGN,
    )
    # the re-found phases (jovian_legs, probe log)
    _floor(system, n34, n56, 14.368686589812885, 7.2781296173548675, kw)


def jovian_legs(n: int) -> None:
    """Re-find the Jovian Table-3 hit's phases one leg at a time. The Wu(3:4-LO) k_u=4 crossing
    depends on tau_u only and the Ws(5:6-LO) k_s=4 crossing on tau_s only, so a 1-D scan of each
    leg toward the STORED crossing gives a seed; refine each leg by a local re-scan; then Newton
    (test settings) from the refined phases. Identity = distance of the result to the stored
    crossing."""
    from cyclerfinder.genome.heteroclinic_cycle import _section_crossing, _seed_on_manifold

    system = jrc.jupiter_europa_system()
    _s, n34 = jrc.build_3_4_lo_node(system)
    _s, n56 = jrc.build_5_6_lo_node(system)
    ref = np.array([-1.3068706902782006, 0.02869075783271665])
    max_time = 3.0 * max(n34.period, n56.period)

    def leg(node: object, tau: float, direction: str) -> np.ndarray | None:
        seed = _seed_on_manifold(
            system,
            node,
            tau=tau,
            direction=direction,
            branch=+1,  # type: ignore[arg-type]
            epsilon=jrc.ANDERSON_LO_EPSILON,
        )
        return _section_crossing(
            system,
            seed,
            direction=direction,
            k=4,
            max_time=max_time,
            ydot_sign=jrc.SECTION_YDOT_SIGN,
            x_sign=jrc.SECTION_X_SIGN,
        )

    best: dict[str, float] = {}
    for name, node, direction in (("u", n34, "unstable"), ("s", n56, "stable")):
        period = float(node.period)  # type: ignore[attr-defined]
        lo, hi, m = 0.0, period, n
        t_best = float("nan")
        for level in range(4):
            taus = [lo + (hi - lo) * (i + 0.5) / m for i in range(m)]
            dists = []
            for t in taus:
                p = leg(node, t % period, direction)
                dists.append(float("inf") if p is None else float(np.linalg.norm(p - ref)))
            i_best = int(np.argmin(dists))
            t_best = taus[i_best]
            width = (hi - lo) / m
            print(
                f"leg {name} level {level}: best tau={t_best!r} dist={dists[i_best]:.3e} "
                f"(reached {sum(np.isfinite(dists))}/{m})",
                flush=True,
            )
            lo, hi, m = t_best - 2 * width, t_best + 2 * width, 40
        best[name] = t_best % period
    t0 = time.time()
    conn = correct_connection(
        system,
        n34,
        n56,
        epsilon=jrc.ANDERSON_LO_EPSILON,
        branch_u=+1,
        branch_s=+1,
        k_u=4,
        k_s=4,
        tau_u0=best["u"],
        tau_s0=best["s"],
        ydot_sign_u=jrc.SECTION_YDOT_SIGN,
        ydot_sign_s=jrc.SECTION_YDOT_SIGN,
        x_sign_u=jrc.SECTION_X_SIGN,
        x_sign_s=jrc.SECTION_X_SIGN,
        max_time_factor=3.0,
        tol=1e-7,
        max_iter=20,
    )
    _show("newton from refined leg phases", conn, t0, (float(ref[0]), float(ref[1])))


def saturn_scans() -> None:
    """Per-combo correct_connection runs behind tests 8 and 9 (find_homoclinic with the tests'
    settings), printing residual and crossing for each combination, converged or not."""
    system = stf.saturn_titan_system()
    _s, node, _r = stc.build_34_node(system)
    common = dict(
        epsilon=stc.EPSILON,
        ydot_sign_u=stc.SECTION_YDOT_SIGN,
        ydot_sign_s=stc.SECTION_YDOT_SIGN,
        x_sign_u=None,
        x_sign_s=None,
        max_time_factor=3.0,
    )
    runs = [
        ("test8", +1, ku, ks, dict(scan_n=8, max_iter=20, tol=1e-7, fd_step=1e-6))
        for ku in (4, 5)
        for ks in (4, 5)
    ] + [
        ("test9", -1, ku, ks, dict(scan_n=12, max_iter=40, tol=1e-9, fd_step=1e-6))
        for ku in (4, 5)
        for ks in (4, 5)
    ]
    for tag, br, ku, ks, kw in runs:
        t0 = time.time()
        c = correct_connection(
            system,
            node,
            node,
            k_u=ku,
            k_s=ks,
            branch_u=br,
            branch_s=br,
            **common,
            **kw,  # type: ignore[arg-type]
        )
        _show(f"{tag} branch={br:+d} k=({ku},{ks})", c, t0, (float("nan"), float("nan")))


def jovian_grid(top: int) -> None:
    """Rank the scan_n=8 grid cells by section gap; Newton (test settings) from the best ``top``."""
    from cyclerfinder.genome.heteroclinic_cycle import _connection_residual

    system = jrc.jupiter_europa_system()
    _s, n34 = jrc.build_3_4_lo_node(system)
    _s, n56 = jrc.build_5_6_lo_node(system)
    ref = (-1.3068706902782006, 0.02869075783271665)
    max_time = 3.0 * max(n34.period, n56.period)
    kw = dict(
        epsilon=jrc.ANDERSON_LO_EPSILON,
        branch_u=+1,
        branch_s=+1,
        k_u=4,
        k_s=4,
        ydot_sign_u=jrc.SECTION_YDOT_SIGN,
        ydot_sign_s=jrc.SECTION_YDOT_SIGN,
        x_sign_u=jrc.SECTION_X_SIGN,
        x_sign_s=jrc.SECTION_X_SIGN,
    )
    n = 8
    cells = []
    for i in range(n):
        for j in range(n):
            tu = n34.period * (i + 0.5) / n
            ts = n56.period * (j + 0.5) / n
            r, _xv = _connection_residual(
                system,
                n34,
                n56,
                tau_u=tu,
                tau_s=ts,
                max_time=max_time,
                **kw,  # type: ignore[arg-type]
            )
            if r is not None:
                cells.append((float(np.linalg.norm(r)), i, j, tu, ts))
    cells.sort()
    print(f"{len(cells)}/64 cells reach the section; best 8:")
    for c in cells[:8]:
        print(f"  gap={c[0]!r} cell=({c[1]},{c[2]}) tau=({c[3]!r},{c[4]!r})", flush=True)
    for c in cells[:top]:
        t0 = time.time()
        conn = correct_connection(
            system,
            n34,
            n56,
            tau_u0=c[3],
            tau_s0=c[4],
            max_time_factor=3.0,
            tol=1e-7,
            max_iter=20,
            **kw,  # type: ignore[arg-type]
        )
        _show(f"newton from cell ({c[1]},{c[2]})", conn, t0, ref)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "case",
        choices=[
            "jovian",
            "tertiary",
            "jovian_grid",
            "tertiary_floor",
            "jovian_legs",
            "jovian_floor",
            "saturn_scans",
            "symmetric_floor",
        ],
    )
    ap.add_argument("--variants", action="store_true")
    a = ap.parse_args()
    cases = {
        "jovian": lambda: jovian(a.variants),
        "tertiary": lambda: tertiary(a.variants),
        "jovian_grid": lambda: jovian_grid(4),
        "jovian_legs": lambda: jovian_legs(160),
        "saturn_scans": saturn_scans,
        "jovian_floor": jovian_floor,
        "tertiary_floor": tertiary_floor,
        "symmetric_floor": symmetric_floor,
    }
    cases[a.case]()
