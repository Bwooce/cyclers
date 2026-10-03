"""#882 -- stroboscopic-map invariant circles and connections in the planar CCR4BP.

Driver for the gate measurements and runs R1-R3 of #882, built on
:mod:`cyclerfinder.search.ccr4bp_strob_connection`. Each ``--stage`` is
independent, finishes in a few minutes, prints timestamped progress and
checkpoints to ``data/found/882_ccr4bp_strob_connection/<stage>.json`` (plus
``.npz`` arrays where useful). Stages that depend on an earlier stage read its
checkpoint.

Stages
------
gates  measured numbers behind the G1-G5 gate tests (the tests assert bounds;
       this records the values).
g5     G5 as specified: the #694 Europa 3:4 pseudospectral torus (built exactly
       as scripts/screen_694_ccr4bp_heteroclinic_search.py builds it) as a
       seed for the stroboscopic invariant-circle corrector.
r1     unperturbed (mu_gan = 0) homoclinic connection of the unstable 3:4
       resonant orbit: cloud, coarse search, refine, verify.
r2     continuation of the R1 connection in mu_gan up to physical Ganymede mass.
r3a    Uranus-Umbriel-Titania, 1:2 exterior Umbriel resonance (crosses Titania's orbit).
r3b    Uranus-Umbriel-Titania, 2:3 and 3:4 exterior unstable members.

No catalogue writeback.

Run:  OMP_NUM_THREADS=1 uv run python scripts/screen_882_ccr4bp_strob_connection.py --stage r1
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

import cyclerfinder.core.ccr4bp as ccr4bp  # noqa: E402
import cyclerfinder.search.ccr4bp_strob_connection as sc  # noqa: E402
from cyclerfinder.core.constants import PLANETS  # noqa: E402
from cyclerfinder.core.satellites import SATELLITES  # noqa: E402

OUT_DIR = ROOT / "data" / "found" / "882_ccr4bp_strob_connection"
T_START = time.time()


def log(msg: str) -> None:
    stamp = time.strftime("%Y-%m-%dT%H:%M:%S")
    print(f"[{stamp} +{time.time() - T_START:7.1f}s] {msg}", flush=True)


def save(name: str, payload: dict[str, Any], arrays: dict[str, np.ndarray] | None = None) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / f"{name}.json").write_text(json.dumps(payload, indent=2, default=_jsonable))
    if arrays:
        np.savez(OUT_DIR / f"{name}.npz", **arrays)
    log(f"wrote {OUT_DIR / name}.json")


def _jsonable(o: Any) -> Any:
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    if dataclasses.is_dataclass(o) and not isinstance(o, type):
        return dataclasses.asdict(o)
    return str(o)


# ---------------------------------------------------------------------------
# Systems, orbits, radii.
# ---------------------------------------------------------------------------


def jeg_radii() -> sc.CollisionRadii:
    """Jupiter (IAU 2015 equatorial), Europa and Ganymede mean radii, in Europa-SMA units."""
    lu = SATELLITES["Europa"].sma_km
    return sc.CollisionRadii(
        planet=PLANETS["J"].radius_eq_km / lu,
        base_moon=SATELLITES["Europa"].radius_eq_km / lu,
        perturber=SATELLITES["Ganymede"].radius_eq_km / lu,
    )


def uut_radii() -> sc.CollisionRadii:
    lu = SATELLITES["Umbriel"].sma_km
    return sc.CollisionRadii(
        planet=PLANETS["U"].radius_eq_km / lu,
        base_moon=SATELLITES["Umbriel"].radius_eq_km / lu,
        perturber=SATELLITES["Titania"].radius_eq_km / lu,
    )


def resonant_orbit(
    mu: float, p: int, q: int, e: float, apse: str = "peri", side: int = 1
) -> tuple[np.ndarray, float, float]:
    """Symmetric p:q orbit from a Keplerian apse guess (see the module helper)."""
    a = (q / p) ** (2.0 / 3.0)
    r = a * (1.0 - e) if apse == "peri" else a * (1.0 + e)
    vin = math.sqrt((1.0 - mu) * (2.0 / r - 1.0 / a))
    x0 = side * r - mu
    vy0 = side * vin - x0
    return sc.symmetric_periodic_orbit(mu, x0, vy0, math.pi * q)


def floquet(mu: float, s4: np.ndarray, period: float) -> np.ndarray:
    """Planar monodromy eigenvalues, integrated independently with cr3bp_stm_eom."""
    from scipy.integrate import solve_ivp

    import cyclerfinder.core.cr3bp as cr3bp

    y0 = np.concatenate([sc.to_state6(s4), np.eye(6).reshape(-1)])
    sol = solve_ivp(
        cr3bp.cr3bp_stm_eom, (0.0, period), y0, args=(mu,), method="DOP853", rtol=1e-13, atol=1e-13
    )
    mono = sol.y[6:, -1].reshape(6, 6)[np.ix_(sc._PLANAR, sc._PLANAR)]
    return np.linalg.eigvals(mono)


def jeg_34_circle(n_nodes: int = 201) -> tuple[sc.InvariantCircle, np.ndarray, float]:
    phys = ccr4bp.jupiter_europa_ganymede_default()
    sys0 = dataclasses.replace(phys, mu_gan=0.0)
    s4, period, res = resonant_orbit(phys.mu, 3, 4, 0.1, "peri", 1)
    assert res < 1e-11, res
    return sc.seed_circle_from_periodic_orbit(sys0, s4, period, n_nodes), s4, period


def circle_summary(c: sc.InvariantCircle) -> dict[str, Any]:
    return {
        "n_nodes": c.n_nodes,
        "t0": c.t0,
        "rho": c.rho,
        "rho_over_2pi": c.rho / (2 * math.pi),
        "residual": c.residual,
        "converged": c.converged,
        "residual_history": list(c.residual_history),
        "fourier_tail": c.fourier_tail(),
        "mu_gan": c.system.mu_gan,
    }


def bundle_summary(b: sc.HyperbolicBundles) -> dict[str, Any]:
    return {
        "lam_u": b.lam_u,
        "lam_s": b.lam_s,
        "lam_product": b.lam_u * b.lam_s,
        "residual_u_offgrid": b.residual_u,
        "residual_s_offgrid": b.residual_s,
        "tail_u": b.tail_u,
        "tail_s": b.tail_s,
    }


def verification_summary(v: sc.ConnectionVerification) -> dict[str, Any]:
    return {
        "genuine": v.genuine,
        "reasons": list(v.reasons),
        "phase_consistent": v.phase_consistent,
        "junction_residual": v.junction_residual,
        "dist_to_circle_per_period": v.dist_to_circle.tolist(),
        "max_excursion": v.max_excursion,
        "final_distance": v.final_distance,
        "offset_size": v.offset_size,
        "integrator_disagreement": v.integrator_disagreement,
        "min_dist_planet_moon_perturber": v.min_dist.tolist(),
        "jacobi_drift": v.jacobi_drift,
        "jacobi_start": v.jacobi_start,
    }


def connection_summary(c: sc.Connection) -> dict[str, Any]:
    d = dataclasses.asdict(c)
    d["singular_values"] = c.singular_values.tolist()
    d["junction_state"] = c.junction_state.tolist()
    return d


def print_verification(v: sc.ConnectionVerification) -> None:
    log(
        f"  verify: genuine={v.genuine} junction={v.junction_residual:.2e} "
        f"max_exc={v.max_excursion:.3e} final={v.final_distance:.3e} "
        f"dop853-vs-radau={v.integrator_disagreement:.2e} "
        f"min dist (planet, moon, perturber)={np.array2string(v.min_dist, precision=4)} "
        f"jacobi drift={v.jacobi_drift:.2e}"
    )
    for k, d in enumerate(v.dist_to_circle):
        print(f"      period {k:3d}: distance to circle {d:.4e}", flush=True)
    for r in v.reasons:
        log(f"  reason: {r}")


# ---------------------------------------------------------------------------
# Stages.
# ---------------------------------------------------------------------------


def stage_gates(_: argparse.Namespace) -> None:
    phys = ccr4bp.jupiter_europa_ganymede_default()
    c0, s4, period = jeg_34_circle()
    out: dict[str, Any] = {
        "orbit_state": s4,
        "orbit_period": period,
        "forcing_period": phys.ganymede_synodic_period,
    }
    # G1
    x = c0.nodes[17] + np.array([1e-3, -2e-3, 5e-4, 1e-3])
    fx, _ = sc.strob_map(phys, x, n=1, t0=0.0)
    back, _ = sc.strob_map(phys, fx, n=-1, t0=phys.ganymede_synodic_period)
    out["g1_inverse_error"] = float(np.max(np.abs(back - x)))
    a, _ = sc.strob_map(phys, x, n=1, t0=2.3)
    b, _ = sc.strob_map(phys, x, n=1, t0=2.3 + phys.ganymede_synodic_period)
    out["g1_t0_vs_t0_plus_P"] = float(np.max(np.abs(a - b)))
    _, phi = sc.strob_map(phys, c0.nodes[40], n=1, t0=1.0, with_stm=True)
    assert phi is not None
    h = 1e-6
    fd = np.empty((4, 4))
    for k in range(4):
        dx = np.zeros(4)
        dx[k] = h
        fp, _ = sc.strob_map(phys, c0.nodes[40] + dx, n=1, t0=1.0)
        fm, _ = sc.strob_map(phys, c0.nodes[40] - dx, n=1, t0=1.0)
        fd[:, k] = (fp - fm) / (2 * h)
    out["g1_stm_vs_fd_rel"] = float(np.max(np.abs(fd - phi)) / np.max(np.abs(phi)))
    _, phi2 = sc.strob_map(phys, c0.nodes[5], n=2, t0=0.7, with_stm=True)
    assert phi2 is not None
    t = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, -1, 1, 0], [1, 0, 0, 1]], dtype=float)
    m = t @ phi2 @ np.linalg.inv(t)
    j = np.block([[np.zeros((2, 2)), np.eye(2)], [-np.eye(2), np.zeros((2, 2))]])
    out["g1_det_minus_1"] = float(np.linalg.det(phi2) - 1.0)
    out["g1_symplectic_error_canonical"] = float(np.max(np.abs(m.T @ j @ m - j)))
    out["g1_stm_max_entry"] = float(np.max(np.abs(phi2)))
    log(f"G1 {json.dumps({k: v for k, v in out.items() if k.startswith('g1')})}")
    # G2
    res0, _ = sc.circle_residual(c0)
    corr = sc.correct_invariant_circle(c0.system, c0.nodes, c0.rho, tol=1e-9)
    import cyclerfinder.core.cr3bp as cr3bp

    jac = [cr3bp.jacobi_constant(sc.to_state6(u), phys.mu) for u in c0.nodes]
    out["g2"] = {
        "n_nodes": c0.n_nodes,
        "seed_residual": res0,
        "after_corrector_residual": corr.residual,
        "corrector_node_change": float(np.max(np.abs(corr.nodes - c0.nodes))),
        "jacobi_spread": max(jac) - min(jac),
        "fourier_tail": c0.fourier_tail(),
        "rho": c0.rho,
    }
    log(f"G2 {out['g2']}")
    # G3
    b0 = sc.hyperbolic_bundles(c0)
    ev = floquet(phys.mu, s4, period)
    lam_f = float(np.max(np.abs(ev)))
    expected = lam_f ** (phys.ganymede_synodic_period / period)
    out["g3"] = {
        **bundle_summary(b0),
        "floquet_multipliers": [complex(e) for e in ev],
        "floquet_lam": lam_f,
        "expected_lam_u": expected,
        "rel_error": b0.lam_u / expected - 1.0,
    }
    log(f"G3 {out['g3']}")
    # G4
    theta, s, n = 0.9, 1.3, 3
    mism = []
    for eps in (1e-4, 5e-5, 2.5e-5):
        pa, _ = sc.manifold_point(c0, b0, "unstable", theta, s, n + 1, eps=eps)
        pb, _ = sc.manifold_point(c0, b0, "unstable", theta + c0.rho, b0.lam_u * s, n, eps=eps)
        mism.append(float(np.linalg.norm(pa - pb)))
    out["g4"] = {"mismatch": mism, "ratios": [mism[0] / mism[1], mism[1] / mism[2]]}
    log(f"G4 {out['g4']}")
    # G5 (hyperbolic circle at physical mass)
    steps = sc.continue_circle_in_mass(c0, phys.mu_gan, verbose=False)
    cp = steps[-1]
    orb = sc.strob_iterates(phys, cp.nodes, n=3, t0=cp.t0)
    target = cp.state(cp.thetas + 3.0 * cp.rho)
    bp = sc.hyperbolic_bundles(cp)
    out["g5_hyperbolic"] = {
        "steps": [circle_summary(c) for c in steps],
        "closure_3_periods": float(np.max(np.abs(orb.states[3] - target))),
        **bundle_summary(bp),
    }
    log(f"G5 (hyperbolic 3:4) {out['g5_hyperbolic']['closure_3_periods']:.2e} lam_u {bp.lam_u}")
    save("gates", out, {"circle_phys_nodes": cp.nodes, "v_u": bp.v_u, "v_s": bp.v_s})


def stage_g5(args: argparse.Namespace) -> None:
    sys.path.insert(0, str(ROOT / "scripts"))
    from screen_694_ccr4bp_heteroclinic_search import _resonant_symmetric_orbit

    import cyclerfinder.search.variational_ccr4bp_torus as vt

    phys = ccr4bp.jupiter_europa_ganymede_default()
    s0, period, res = _resonant_symmetric_orbit(phys.mu, 3, 4)
    log(f"#694 seed orbit: perp residual {res:.2e}, period {period:.6f}")
    ev = floquet(phys.mu, sc.to_state4(s0), period)
    log(f"#694 seed orbit Floquet multipliers (mu_gan=0): {ev}")
    torus = vt.discover_ccr4bp_torus_from_resonant_orbit(
        phys,
        s0,
        period,
        n1=1,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    log(
        f"#694 torus: residual_rms {torus.residual_rms:.3e}, "
        f"closure {torus.closure_residual:.3e}, rho_strob {torus.rho_strob:.6f}"
    )
    rec: dict[str, Any] = {
        "seed_orbit_floquet": [complex(e) for e in ev],
        "torus_residual_rms": torus.residual_rms,
        "torus_closure_residual": torus.closure_residual,
        "torus_rho_strob": torus.rho_strob,
        "runs": [],
    }
    radii = jeg_radii()
    for n_nodes in args.nodes:
        seed = sc.seed_circle_from_pseudospectral(torus, n_nodes, theta1=0.0)
        res0, _ = sc.circle_residual(seed)
        orb = sc.strob_iterates(phys, seed.nodes, n=1, t0=seed.t0, radii=radii)
        log(
            f"N={n_nodes}: seed residual {res0:.3e}; "
            f"min dist to Ganymede over one period {orb.min_dist[:, 2].min():.4e}"
        )
        c = sc.correct_invariant_circle(
            phys, seed.nodes, seed.rho, t0=seed.t0, tol=1e-10, max_iter=15, verbose=True
        )
        amps = c.fourier_amplitudes()
        entry = {
            "n_nodes": n_nodes,
            "seed_residual": res0,
            "seed_min_dist_ganymede": float(orb.min_dist[:, 2].min()),
            **circle_summary(c),
            "fourier_amplitudes": amps.tolist(),
            "seed_fourier_amplitudes": seed.fourier_amplitudes().tolist(),
        }
        if c.converged:
            o3 = sc.strob_iterates(phys, c.nodes, n=3, t0=c.t0)
            entry["closure_3_periods"] = float(
                np.max(np.abs(o3.states[3] - c.state(c.thetas + 3 * c.rho)))
            )
            try:
                b = sc.hyperbolic_bundles(c)
                entry["bundles"] = bundle_summary(b)
            except ValueError as exc:
                entry["bundles"] = f"not hyperbolic: {exc}"
        log(f"N={n_nodes}: {circle_summary(c)}")
        rec["runs"].append(entry)
        save("g5", rec)


def _r1_setup(eps: float) -> tuple[sc.InvariantCircle, sc.HyperbolicBundles, np.ndarray, float]:
    c0, s4, period = jeg_34_circle()
    b0 = sc.hyperbolic_bundles(c0)
    log(f"circle N={c0.n_nodes} rho={c0.rho:.6f}; lam_u={b0.lam_u:.6f} lam_s={b0.lam_s:.6f}")
    return c0, b0, s4, period


def stage_r1(args: argparse.Namespace) -> None:
    import cyclerfinder.core.cr3bp as cr3bp

    eps = args.eps
    c0, b0, s4, period = _r1_setup(eps)
    radii = jeg_radii()
    c_orbit = cr3bp.jacobi_constant(sc.to_state6(s4), c0.system.mu)
    clouds: dict[tuple[str, float], sc.ManifoldCloud] = {}
    for branch in ("unstable", "stable"):
        for sign in (1.0, -1.0):
            cl = sc.manifold_cloud(
                c0,
                b0,
                branch,
                n_theta=args.n_theta,
                n_s=args.n_s,
                n_max=args.n_max,
                eps=eps,
                sign=sign,
                radii=radii,
            )
            clouds[(branch, sign)] = cl
            log(
                f"cloud {branch} sign {sign:+.0f}: {cl.thetas.size} trajectories x {args.n_max} "
                f"periods, dropped {cl.n_dropped}"
            )
    cands: list[sc.ConnectionCandidate] = []
    for su in (1.0, -1.0):
        for ss in (1.0, -1.0):
            cc = sc.coarse_intersections(
                clouds[("unstable", su)],
                clouds[("stable", ss)],
                n_best=args.n_cand,
                min_excursion=args.min_exc,
                circle=c0,
            )
            log(
                f"coarse (u{su:+.0f}, s{ss:+.0f}): best distances "
                f"{[round(c.distance, 5) for c in cc[:5]]} "
                f"(n_u, n_s) {[(c.n_u, c.n_s) for c in cc[:5]]}"
            )
            cands.extend(cc)
    cands.sort(key=lambda c: c.distance)
    results: list[dict[str, Any]] = []
    verified = 0
    for cand in cands[: args.n_refine]:
        con = sc.refine_connection(c0, b0, c0, b0, cand, eps=eps)
        log(
            f"refine {cand.n_u},{cand.n_s} signs ({cand.sign_u:+.0f},{cand.sign_s:+.0f}) "
            f"coarse {cand.distance:.3e} -> "
            f"residual {con.residual:.3e} in {con.n_iter} it; "
            f"sv {np.array2string(con.singular_values, precision=3)}"
        )
        entry: dict[str, Any] = {
            "candidate": dataclasses.asdict(cand),
            "connection": connection_summary(con),
        }
        if con.converged and verified < args.n_verify:
            v = sc.verify_connection(c0, b0, c0, b0, con)
            print_verification(v)
            entry["verification"] = verification_summary(v)
            entry["jacobi_departure_minus_orbit"] = v.jacobi_start - c_orbit
            log(
                f"  Jacobi: orbit {c_orbit:.12f}, departure - orbit "
                f"{v.jacobi_start - c_orbit:.3e}, drift {v.jacobi_drift:.3e}"
            )
            verified += 1
        results.append(entry)
        save(
            f"r1_eps{eps:g}",
            {
                "eps": eps,
                "circle": circle_summary(c0),
                "bundles": bundle_summary(b0),
                "orbit_state": s4,
                "orbit_period": period,
                "jacobi_orbit": c_orbit,
                "clouds": {
                    f"{k[0]}{k[1]:+.0f}": {"dropped": v.n_dropped, "n": int(v.thetas.size)}
                    for k, v in clouds.items()
                },
                "results": results,
            },
        )


def _load_r1(eps: float) -> dict[str, Any]:
    return json.loads((OUT_DIR / f"r1_eps{eps:g}.json").read_text())


def stage_r2(args: argparse.Namespace) -> None:
    """Continue R1's verified connection in mu_gan."""
    eps = args.eps
    r1 = _load_r1(eps)
    good = [r for r in r1["results"] if r.get("verification", {}).get("genuine")]
    if not good:
        log("no genuine R1 connection to continue")
        return
    base = good[0]["connection"]
    phys = ccr4bp.jupiter_europa_ganymede_default()
    c0, b0, _, _ = _r1_setup(eps)
    # Unperturbed family: time-shifted copies of the homoclinic, by seeding the
    # unperturbed refine at shifted theta (rank-deficient Newton slides onto the curve).
    shifts = np.linspace(0.0, 2 * math.pi, args.n_family, endpoint=False)
    family: list[sc.Connection] = []
    for dth in shifts:
        cand = sc.ConnectionCandidate(
            theta_u=base["theta_u"] + dth,
            s_u=base["s_u"],
            theta_s=base["theta_s"] + dth,
            s_s=base["s_s"],
            n_u=base["n_u"],
            n_s=base["n_s"],
            sign_u=base["sign_u"],
            sign_s=base["sign_s"],
            distance=0.0,
        )
        con = sc.refine_connection(c0, b0, c0, b0, cand, eps=eps)
        family.append(con)
    ok = [f for f in family if f.converged]
    log(f"unperturbed family: {len(ok)}/{len(family)} members converged")
    fracs = [float(f) for f in args.fractions]
    # first perturbed step: Jacobi-mismatch scan along the family to bracket Melnikov zeros
    sys1 = dataclasses.replace(phys, mu_gan=fracs[0] * phys.mu_gan)
    steps = sc.continue_circle_in_mass(c0, sys1.mu_gan, fractions=(1.0,))
    c1 = steps[-1]
    b1 = sc.hyperbolic_bundles(c1, ref_u=b0.v_u, ref_s=b0.v_s)
    log(f"mu_gan={sys1.mu_gan:.3e}: circle {c1.residual:.2e}, lam_u {b1.lam_u:.6f}")
    scan = []
    for f in ok:
        xu, _ = sc.manifold_point(
            c1, b1, "unstable", f.theta_u, f.s_u, f.n_u, eps=eps, sign=f.sign_u
        )
        xs, _ = sc.manifold_point(c1, b1, "stable", f.theta_s, f.s_s, f.n_s, eps=eps, sign=f.sign_s)
        dc = float(sc._jacobi4(xu, phys.mu) - sc._jacobi4(xs, phys.mu))
        scan.append((f, dc, float(np.linalg.norm(xu - xs))))
    for f, dc, gap in scan:
        print(
            f"      theta_u {f.theta_u:.4f} s_u {f.s_u:.4f}: dC {dc:+.3e}  gap {gap:.3e}",
            flush=True,
        )
    seeds = []
    for i in range(len(scan)):
        a, b = scan[i], scan[(i + 1) % len(scan)]
        if a[1] * b[1] < 0:
            seeds.append(a[0] if abs(a[1]) < abs(b[1]) else b[0])
    log(f"{len(seeds)} sign changes of the Jacobi mismatch along the family")
    rec: dict[str, Any] = {
        "eps": eps,
        "fractions": fracs,
        "family_size": len(ok),
        "scan": [{"theta_u": f.theta_u, "s_u": f.s_u, "dC": dc, "gap": g} for f, dc, g in scan],
        "branches": [],
    }
    for bi, seed in enumerate(seeds):
        branch: list[dict[str, Any]] = []
        cur_conn = seed
        cur_circle, cur_b = c0, b0
        for frac in fracs:
            mg = frac * phys.mu_gan
            st = sc.continue_circle_in_mass(cur_circle, mg, fractions=(1.0,))
            circ = st[-1]
            if not circ.converged:
                log(f"branch {bi}: circle failed at mu_gan frac {frac}: {circ.residual_history}")
                branch.append({"frac": frac, "circle": circle_summary(circ), "failed": "circle"})
                break
            bun = sc.hyperbolic_bundles(circ, ref_u=cur_b.v_u, ref_s=cur_b.v_s)
            cand = sc.ConnectionCandidate(
                theta_u=cur_conn.theta_u,
                s_u=cur_conn.s_u,
                theta_s=cur_conn.theta_s,
                s_s=cur_conn.s_s,
                n_u=cur_conn.n_u,
                n_s=cur_conn.n_s,
                sign_u=cur_conn.sign_u,
                sign_s=cur_conn.sign_s,
                distance=0.0,
            )
            con = sc.refine_connection(circ, bun, circ, bun, cand, eps=eps)
            e: dict[str, Any] = {
                "frac": frac,
                "mu_gan": mg,
                "circle": circle_summary(circ),
                "bundles": bundle_summary(bun),
                "connection": connection_summary(con),
            }
            log(
                f"branch {bi} frac {frac:g}: circle {circ.residual:.1e} lam_u {bun.lam_u:.5f}; "
                f"connection residual "
                f"{con.residual:.2e} sv {np.array2string(con.singular_values, precision=3)} "
                f"theta_u {con.theta_u:.5f}"
            )
            if not con.converged:
                e["failed"] = "connection"
                branch.append(e)
                break
            if frac in (fracs[0], fracs[-1]) or args.verify_all:
                v = sc.verify_connection(circ, bun, circ, bun, con)
                print_verification(v)
                e["verification"] = verification_summary(v)
            branch.append(e)
            cur_conn, cur_circle, cur_b = con, circ, bun
            rec["branches"] = [*rec["branches"][:bi], branch]
            save(f"r2_eps{eps:g}", rec)
        rec["branches"] = [*rec["branches"][:bi], branch]
        save(f"r2_eps{eps:g}", rec)


def _r3_system() -> ccr4bp.CCR4BPSystem:
    import cyclerfinder.core.ccr4bp_umbriel_titania as ut

    return ut.uranus_umbriel_titania_default()


def stage_r3a(args: argparse.Namespace) -> None:
    sys.path.insert(0, str(ROOT / "scripts"))
    from screen_701_ccr4bp_umbriel_titania_search import _resonant_symmetric_orbit

    uut = _r3_system()
    s0, period, res = _resonant_symmetric_orbit(uut.mu, 1, 2)
    ev = floquet(uut.mu, sc.to_state4(s0), period)
    log(f"Umbriel 1:2 orbit: residual {res:.2e}, period {period:.5f}, Floquet {ev}")
    sys0 = dataclasses.replace(uut, mu_gan=0.0)
    rec: dict[str, Any] = {"floquet": [complex(e) for e in ev], "runs": []}
    for n_nodes in args.nodes:
        c0 = sc.seed_circle_from_periodic_orbit(sys0, sc.to_state4(s0), period, n_nodes)
        r = np.hypot(c0.nodes[:, 0] + uut.mu, c0.nodes[:, 1])
        log(
            f"N={n_nodes}: orbit radius range [{r.min():.4f}, {r.max():.4f}], "
            f"Titania at {uut.a_gan:.4f}, rho {c0.rho:.5f}"
        )
        steps = sc.continue_circle_in_mass(
            c0, uut.mu_gan, fractions=tuple(args.fractions), verbose=True
        )
        for st in steps:
            log(f"   mu_gan {st.system.mu_gan:.3e}: {circle_summary(st)}")
        last = steps[-1]
        rec["runs"].append(
            {
                "n_nodes": n_nodes,
                "radius_range": [float(r.min()), float(r.max())],
                "steps": [circle_summary(s) for s in steps],
                "last_fourier_amplitudes": last.fourier_amplitudes().tolist(),
            }
        )
        save("r3a", rec)


def stage_r3b(args: argparse.Namespace) -> None:
    uut = _r3_system()
    radii = uut_radii()
    sys0 = dataclasses.replace(uut, mu_gan=0.0)
    rec: dict[str, Any] = {"orbits": []}
    for p, q in ((2, 3), (3, 4)):
        for e in args.ecc:
            for apse in ("peri", "apo"):
                s4, period, res = resonant_orbit(uut.mu, p, q, e, apse, 1)
                if res > 1e-10:
                    continue
                ev = floquet(uut.mu, s4, period)
                c0 = sc.seed_circle_from_periodic_orbit(sys0, s4, period, args.nodes[0])
                r = np.hypot(c0.nodes[:, 0] + uut.mu, c0.nodes[:, 1])
                dmoon = float(np.min(np.hypot(c0.nodes[:, 0] - 1 + uut.mu, c0.nodes[:, 1])))
                lam = float(np.max(np.abs(ev)))
                hyper = bool(
                    np.any(np.abs(np.abs(ev) - 1.0) > 1e-4) and np.all(np.abs(ev.imag) < 1e-9)
                )
                entry: dict[str, Any] = {
                    "p": p,
                    "q": q,
                    "e_seed": e,
                    "apse": apse,
                    "period": period,
                    "floquet_max": lam,
                    "hyperbolic": hyper,
                    "r_range": [float(r.min()), float(r.max())],
                    "min_dist_umbriel_nodes": dmoon,
                }
                log(
                    f"{p}:{q} e={e} {apse}: T={period:.4f} lam={lam:.4f} hyperbolic={hyper} "
                    f"r=[{r.min():.3f},{r.max():.3f}] dUmb={dmoon:.3f}"
                )
                rec["orbits"].append(entry)
                if not hyper or r.max() > uut.a_gan - 0.1:
                    continue
                steps = sc.continue_circle_in_mass(c0, uut.mu_gan, fractions=tuple(args.fractions))
                last = steps[-1]
                entry["circle_steps"] = [circle_summary(s) for s in steps]
                log(f"   circle at physical Titania mass: {circle_summary(last)}")
                if not last.converged:
                    save("r3b", rec)
                    continue
                b = sc.hyperbolic_bundles(last)
                entry["bundles"] = bundle_summary(b)
                log(f"   bundles: {bundle_summary(b)}")
                save("r3b", rec)
                if args.homoclinic and lam > 1.5:
                    entry["homoclinic"] = _homoclinic_search(last, b, radii, args)
                    save("r3b", rec)
    save("r3b", rec)


def _homoclinic_search(
    c: sc.InvariantCircle,
    b: sc.HyperbolicBundles,
    radii: sc.CollisionRadii,
    args: argparse.Namespace,
) -> list[dict[str, Any]]:
    eps = args.eps
    clouds = {}
    for branch in ("unstable", "stable"):
        for sign in (1.0, -1.0):
            clouds[(branch, sign)] = sc.manifold_cloud(
                c,
                b,
                branch,
                n_theta=args.n_theta,
                n_s=args.n_s,
                n_max=args.n_max,
                eps=eps,
                sign=sign,
                radii=radii,
            )
            log(f"   cloud {branch}{sign:+.0f}: dropped {clouds[(branch, sign)].n_dropped}")
    cands: list[sc.ConnectionCandidate] = []
    for su in (1.0, -1.0):
        for ss in (1.0, -1.0):
            cands.extend(
                sc.coarse_intersections(
                    clouds[("unstable", su)],
                    clouds[("stable", ss)],
                    n_best=args.n_cand,
                    min_excursion=args.min_exc,
                    circle=c,
                )
            )
    cands.sort(key=lambda x: x.distance)
    out = []
    verified = 0
    for cand in cands[: args.n_refine]:
        con = sc.refine_connection(c, b, c, b, cand, eps=eps)
        log(
            f"   refine ({cand.n_u},{cand.n_s}) coarse {cand.distance:.3e} -> {con.residual:.3e} "
            f"sv {np.array2string(con.singular_values, precision=3)}"
        )
        e: dict[str, Any] = {
            "candidate": dataclasses.asdict(cand),
            "connection": connection_summary(con),
        }
        if con.converged and verified < args.n_verify:
            v = sc.verify_connection(c, b, c, b, con)
            print_verification(v)
            e["verification"] = verification_summary(v)
            verified += 1
        out.append(e)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--stage", required=True, choices=["gates", "g5", "r1", "r2", "r3a", "r3b"])
    ap.add_argument("--eps", type=float, default=1e-4)
    ap.add_argument("--n-theta", type=int, default=96)
    ap.add_argument("--n-s", type=int, default=6)
    ap.add_argument("--n-max", type=int, default=14)
    ap.add_argument("--n-cand", type=int, default=8)
    ap.add_argument("--n-refine", type=int, default=8)
    ap.add_argument("--n-verify", type=int, default=2)
    ap.add_argument("--min-exc", type=float, default=0.05)
    ap.add_argument("--n-family", type=int, default=24)
    ap.add_argument("--verify-all", action="store_true")
    ap.add_argument("--homoclinic", action="store_true")
    ap.add_argument("--nodes", type=int, nargs="+", default=[201])
    ap.add_argument("--ecc", type=float, nargs="+", default=[0.05, 0.1, 0.15])
    ap.add_argument(
        "--fractions", type=float, nargs="+", default=[0.001, 0.01, 0.03, 0.1, 0.25, 0.5, 0.75, 1.0]
    )
    args = ap.parse_args()
    log(f"stage {args.stage} start")
    {
        "gates": stage_gates,
        "g5": stage_g5,
        "r1": stage_r1,
        "r2": stage_r2,
        "r3a": stage_r3a,
        "r3b": stage_r3b,
    }[args.stage](args)
    log(f"stage {args.stage} done")


if __name__ == "__main__":
    main()
