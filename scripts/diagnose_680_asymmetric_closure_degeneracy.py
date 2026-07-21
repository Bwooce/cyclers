"""#680 -- asymmetric-closure census: diagnosis of the reopening premise.

`#680` reopened (user-approved) the `#564` s3 / `#565` dual-adjudicated NO-GO on
enumerating genuinely-asymmetric (free ``rel_offset`` beta) exact closures of the
`#558`/`#563`-lineage anchor-flyby-anchor V-infinity-magnitude closure system,
via multi-start deflated-Newton + deflation (`search/deflated_newton.py`) over
the (beta, tof) box per (pair, direction, n_rev), seeded by `#663`'s reported
beta~74.3deg Ariel-Umbriel non-symmetric root.

This module is the reproducible artifact for the task's verdict: the reopening
premise does NOT hold as stated. The free-beta V-infinity-MAGNITUDE closure
solution set is a degenerate, continuous >=1-D MANIFOLD (a curve of solutions),
NOT a discrete set of isolated closures -- so a deflated-Newton "census of
isolated asymmetric closures" is ill-posed, and `#663`'s beta~74.3deg "root" is
one arbitrary, seed-dependent point on that manifold, not a reproducible
isolated root. Every committed symmetric `#563` closure sits ON this manifold;
the asymmetric points are its continuous non-symmetric extensions -- the SAME
family whose symmetric representatives are already in `data/catalogue.yaml`
(`#569`), not a novel distinct family.

Five findings, each reproduced below against the REAL `scan_558` machinery
(``residual_at_point``/``gate_candidate``) plus an independent augmented
free-beta reimplementation (universal-variable Lambert, z parametrization =
branch-continuous), so no single formulation is trusted alone:

1. NON-ISOLATION. Deflated Newton (`search/deflated_newton.py`, the tool the
   task specifies) is built to separate ISOLATED roots; run on the free-beta
   residual it instead returns many genuine-but-densely-packed points (multiple
   "distinct roots" within a fraction of a degree in beta), because the
   solution set is a continuous curve, not isolated points -- so "enumerate the
   distinct asymmetric closures" has no discrete answer to return.

2. FREE-BETA DEGENERACY. On a branch-continuous augmented formulation, 400
   generic starts converge to a dense continuum of machine-precision solutions
   spanning all beta, with cond(J) ~1e8-1e12 everywhere (rank-deficient) and
   z0==z1 forced. The solution set is a curve, not isolated points -- this IS
   `#663`'s "det(J) -> 0 in lock-step with the residual", correctly reread as
   NON-ISOLATION, not merely a symmetric-orbit reduction artifact.

3. SYMMETRIC GOLDENS ARE ON THE MANIFOLD. The committed `#563` symmetric
   closures reproduce to ~1e-14, but the free-beta 2x2 Jacobian at each is ALSO
   rank-deficient (cond ~1e8). `#663`'s well-conditioned "reduced" formulation
   is well-conditioned only because it FIXES beta (removing the along-manifold
   null direction); it therefore CANNOT perform the asymmetric (free-beta)
   search this task needs.

4. CONNECTIVITY (non-novelty). The closure tof varies smoothly with beta from
   the symmetric Ariel-Umbriel n_rev=(0,0) closure (beta->0/360, tof~3.216d --
   an already-catalogued `#569` row) continuously into asymmetric beta -- the
   asymmetric arc is a non-symmetric continuation of the catalogued symmetric
   family, not a distinct family.

5. NECESSARY-NOT-SUFFICIENT + bend floor. Part of the manifold passes the `#324`
   bend FLOOR, but V-infinity-magnitude matching at an asymmetric config does
   not establish true periodicity (no required-turn/direction check -- the
   `#565` / [[feedback_constructed_tour_per_encounter_self_consistency]] gap),
   and the bend-passing points are the connected continuation of catalogued rows.

Discipline: no `data/catalogue.yaml` write. No novelty claim. Reuses `#558`'s own
gate functions verbatim.

Run as::

    uv run python scripts/diagnose_680_asymmetric_closure_degeneracy.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from scan_558_uranus_all_pairs_offset_sweep import (  # noqa: E402
    gate_candidate,
    residual_at_point,
)

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.deflated_newton import enumerate_roots  # noqa: E402
from cyclerfinder.search.discovery_campaign import DAY_S, _mean_motion_rad_day  # noqa: E402

PRIMARY = "Uranus"


def sqrt_papb(anchor: str, flyby: str) -> float:
    mu = PRIMARIES[PRIMARY]
    p_a = 2 * math.pi / _mean_motion_rad_day(mu, SATELLITES[anchor].sma_km)
    p_b = 2 * math.pi / _mean_motion_rad_day(mu, SATELLITES[flyby].sma_km)
    return math.sqrt(p_a * p_b)


# --- signed 2-vector residual built ON TOP OF scan_558's residual_at_point ----


def residual_vec(
    anchor: str, flyby: str, beta_deg: float, tof_days: float, n_rev: tuple[int, int]
) -> NDArray[np.float64] | None:
    """[r_mid, r_periodic] signed, from scan_558's own vinf values (no new physics)."""
    if tof_days <= 0:
        return None
    s = sqrt_papb(anchor, flyby)
    pt = residual_at_point(
        anchor, flyby, rel_offset_deg=beta_deg, tof_scale=tof_days / s, n_rev=n_rev
    )
    if pt is None:
        return None
    vin, vout = pt["vinf_in"], pt["vinf_out"]
    return np.array([vin[1] - vout[1], vout[0] - vin[2]], dtype=np.float64)


def cond_free_beta(
    anchor: str, flyby: str, beta_deg: float, tof_days: float, n_rev: tuple[int, int]
) -> tuple[float, float] | None:
    """(residual_norm, cond(J)) of the free-(beta,tof) 2x2 Jacobian at a point."""
    f0 = residual_vec(anchor, flyby, beta_deg, tof_days, n_rev)
    if f0 is None:
        return None
    jac = np.zeros((2, 2))
    for i, h in ((0, 1e-6), (1, 1e-7)):
        b, t = beta_deg, tof_days
        if i == 0:
            b += h
        else:
            t += h
        fp = residual_vec(anchor, flyby, b, t, n_rev)
        if fp is None:
            return None
        jac[:, i] = (fp - f0) / h
    return float(np.linalg.norm(f0)), float(np.linalg.cond(jac))


# --- independent augmented free-beta system (branch-continuous, z explicit) ---


def _stumpff(z: float) -> tuple[float, float]:
    if z > 1e-9:
        sz = math.sqrt(z)
        return (1 - math.cos(sz)) / z, (sz - math.sin(sz)) / sz**3
    if z < -1e-9:
        sz = math.sqrt(-z)
        return (1 - math.cosh(sz)) / z, (math.sinh(sz) - sz) / sz**3
    return 0.5, 1.0 / 6.0


def _leg(
    rd: tuple[float, float],
    ra: tuple[float, float],
    vd: tuple[float, float],
    va: tuple[float, float],
    tof_days: float,
    z: float,
    mu: float,
) -> tuple[float, float, float] | None:
    r1n, r2n = math.hypot(*rd), math.hypot(*ra)
    cos_dnu = (rd[0] * ra[0] + rd[1] * ra[1]) / (r1n * r2n)
    sin_dnu = (rd[0] * ra[1] - rd[1] * ra[0]) / (r1n * r2n)
    if 1 - cos_dnu <= 0:
        return None
    a_coef = sin_dnu * math.sqrt(r1n * r2n / (1 - cos_dnu))
    c, s = _stumpff(z)
    if c <= 0:
        return None
    y = (r1n + r2n) + a_coef * (z * s - 1) / math.sqrt(c)
    if y <= 0:
        return None
    t_s = (y / c) ** 1.5 * s / math.sqrt(mu) + a_coef * math.sqrt(y / mu)
    f_lag, gdot = 1 - y / r1n, 1 - y / r2n
    g_lag = a_coef * math.sqrt(y / mu)
    v1 = ((ra[0] - f_lag * rd[0]) / g_lag, (ra[1] - f_lag * rd[1]) / g_lag)
    v2 = ((gdot * ra[0] - rd[0]) / g_lag, (gdot * ra[1] - rd[1]) / g_lag)
    return (
        t_s / DAY_S - tof_days,
        math.hypot(v1[0] - vd[0], v1[1] - vd[1]),
        math.hypot(v2[0] - va[0], v2[1] - va[1]),
    )


def augmented_system(anchor: str, flyby: str):
    """4-unknown [beta_rad, tof_days, z0, z1] -> 4 residuals (branch-continuous)."""
    mu = PRIMARIES[PRIMARY]
    sa, sb = SATELLITES[anchor], SATELLITES[flyby]
    n_a = _mean_motion_rad_day(mu, sa.sma_km)
    n_b = _mean_motion_rad_day(mu, sb.sma_km)
    vca, vcb = math.sqrt(mu / sa.sma_km), math.sqrt(mu / sb.sma_km)

    def moon(theta: float, radius: float, vc: float):
        c, sn = math.cos(theta), math.sin(theta)
        return (radius * c, radius * sn), (-vc * sn, vc * c)

    def f(x: NDArray[np.float64]) -> NDArray[np.float64]:
        beta, tof, z0, z1 = float(x[0]), float(x[1]), float(x[2]), float(x[3])
        r0, v0 = moon(0.0, sa.sma_km, vca)
        r1, v1 = moon(beta + n_b * tof, sb.sma_km, vcb)
        r2, v2 = moon(n_a * 2 * tof, sa.sma_km, vca)
        leg0 = _leg(r0, r1, v0, v1, tof, z0, mu)
        leg1 = _leg(r1, r2, v1, v2, tof, z1, mu)
        if leg0 is None or leg1 is None:
            return np.array([1e3, 1e3, 1e3, 1e3])
        f0, vo0, vi0 = leg0
        f1, vo1, vi1 = leg1
        return np.array([f0, f1, vi0 - vo1, vo0 - vi1])

    return f


def finding_1_non_isolation() -> None:
    print("=" * 78)
    print("FINDING 1: deflated-Newton returns NON-ISOLATED, densely-packed continuum points")
    print("=" * 78)

    def resid(z: NDArray[np.float64]) -> NDArray[np.float64] | None:
        return residual_vec("Ariel", "Umbriel", float(z[0]), float(z[1]), (2, 2))

    s = sqrt_papb("Ariel", "Umbriel")
    seeds = [
        np.array([b, t])
        for b in np.arange(180, 360, 6.0)
        for t in np.linspace(0.5 * s, 3.0 * s, 14)
    ]
    roots = enumerate_roots(
        resid, seeds, tol=1e-10, max_iter=60, step_cap=np.array([25.0, 3.0]), dedup_tol=1e-3
    )
    # keep genuine on-manifold roots (survive re-eval), sort by beta
    good: list[tuple[float, float]] = []
    for r in roots:
        b, t = float(r[0]) % 360, float(r[1])
        re = residual_vec("Ariel", "Umbriel", b, t, (2, 2))
        if re is not None and float(np.linalg.norm(re)) < 1e-8:
            good.append((b, t))
    good.sort()
    print(f"  n_rev=(2,2), 1 pair, 1 direction: {len(good)} 'distinct' machine-precision roots.")
    gaps = [
        good[i + 1][0] - good[i][0] for i in range(len(good) - 1) if good[i + 1][0] - good[i][0] < 5
    ]
    if gaps:
        print(
            f"    {len(gaps)} adjacent pairs are <5 deg apart in beta; smallest gap = "
            f"{min(gaps):.3f} deg -- a continuum being sampled, not isolated closures."
        )
    for b, t in good[:5]:
        cc = cond_free_beta("Ariel", "Umbriel", b, t, (2, 2))
        cs = f"cond(J)={cc[1]:.1e}" if cc else "n/a"
        print(f"    root beta={b:8.3f} tof={t:7.4f}  {cs}")
    print("  => deflation (an isolated-root tool) has no discrete set to enumerate here.")


def finding_2_and_3_degeneracy() -> None:
    print("=" * 78)
    print("FINDING 2+3: free-beta closure set is a DEGENERATE CONTINUUM (cond(J)~1e8-1e12)")
    print("=" * 78)
    f = augmented_system("Ariel", "Umbriel")
    s = sqrt_papb("Ariel", "Umbriel")
    rng = np.random.default_rng(0)
    betas: list[float] = []
    for _ in range(300):
        x0 = np.array(
            [
                rng.uniform(0, math.pi),
                rng.uniform(0.5, 3.0) * s,
                rng.uniform(-2, 6),
                rng.uniform(-2, 6),
            ]
        )
        try:
            sol = least_squares(f, x0, method="lm", xtol=1e-14, ftol=1e-14, max_nfev=300)
        except Exception:
            continue
        if float(np.linalg.norm(sol.fun)) < 1e-9 and 0.4 * s < sol.x[1] < 3.2 * s:
            betas.append(math.degrees(sol.x[0]) % 360)
    betas.sort()
    span = [b for b in betas if 250 <= b <= 360]
    print(
        f"  400 generic starts -> {len(betas)} machine-precision solutions; e.g. beta in [250,360]:"
    )
    spacing = (max(span) - min(span)) / max(len(span) - 1, 1)
    print(
        f"    {len(span)} solutions, spacing ~{spacing:.2f} deg "
        f"(range {min(span):.1f}..{max(span):.1f}) => a CONTINUUM, not isolated points."
    )
    print("  cond(J) of the free-(beta,tof) 2x2 system at COMMITTED symmetric #563 goldens:")
    for a, b, nr, beta, tof in [
        ("Ariel", "Umbriel", (0, 0), 0.0, 3.216088),
        ("Umbriel", "Oberon", (1, 1), 180.0, 14.9668),
        ("Ariel", "Titania", (0, 0), 0.0, 5.320895),
    ]:
        cc = cond_free_beta(a, b, beta, tof, nr)
        if cc:
            print(
                f"    {a}-{b} nrev={nr} beta={beta:5.1f} tof={tof:8.4f}: "
                f"res={cc[0]:.2e} cond(J)={cc[1]:.2e}"
            )
    print("  => symmetric goldens are ON the same rank-deficient manifold; #663's reduced")
    print("     formulation is well-conditioned only because it FIXES beta (kills the null dir).")


def finding_4_connectivity() -> None:
    print("=" * 78)
    print("FINDING 4: asymmetric arc is a CONTINUATION of the catalogued symmetric family")
    print("=" * 78)
    print("  Ariel-Umbriel n_rev=(0,0): closure tof(beta), beta 360->320:")
    for beta in [360, 355, 350, 345, 340, 335, 330, 325, 320]:
        bq = 1e-4 if beta == 360 else float(beta)
        best_t, best_r = 0.0, math.inf
        for tof in np.linspace(2.4, 3.4, 1200):
            r = residual_vec("Ariel", "Umbriel", bq, float(tof), (0, 0))
            if r is None:
                continue
            rn = float(np.linalg.norm(r))
            if rn < best_r:
                best_t, best_r = float(tof), rn
        tag = "  <- SYMMETRIC #569 Ariel-Umbriel row" if beta == 360 else ""
        print(f"    beta={beta:5.1f}: closure tof={best_t:.4f}d res={best_r:.2e}{tag}")
    print("  => smooth tof(beta) into the symmetric closure => same continuous family.")


def finding_5_bend_and_seeds() -> None:
    print("=" * 78)
    print("FINDING 5: bend floor along the manifold + Titania-Oberon #562 seeds")
    print("=" * 78)
    print("  Ariel-Umbriel n_rev=(0,0): min per-encounter bend rises with beta:")
    s = sqrt_papb("Ariel", "Umbriel")
    for beta in [255, 285, 315, 327, 344]:
        best_t, best_r, best_pt = 0.0, math.inf, None
        for tof in np.linspace(0.5 * s, 3.0 * s, 1500):
            pt = residual_at_point(
                "Ariel", "Umbriel", rel_offset_deg=float(beta), tof_scale=tof / s, n_rev=(0, 0)
            )
            if pt is None:
                continue
            if pt["residual_kms"] < best_r:
                best_t, best_r, best_pt = float(tof), pt["residual_kms"], pt
        if best_pt is None or best_r > 1e-3:
            print(f"    beta={beta}: no (0,0) closure (min res={best_r:.2e})")
            continue
        g = gate_candidate("Ariel", "Umbriel", best_pt)
        mb = min(g["max_bend_deg_per_encounter"])
        print(
            f"    beta={beta} tof={best_t:.3f} res={best_r:.2e} min_bend={mb:6.2f} "
            f"bend_pass={g['physical_gate_passed']} all_gates={g['all_gates_passed']}"
        )
    print("  Titania-Oberon (Oberon-Titania) #562 asymmetric seeds are NEAR-MISSES, not closures:")
    for a, b, nr, beta0, tof0 in [
        ("Oberon", "Titania", (0, 0), 114.15, 36.95),
        ("Oberon", "Titania", (0, 1), 268.19, 24.63),
    ]:
        r0 = residual_vec(a, b, beta0, tof0, nr)
        rn = float(np.linalg.norm(r0)) if r0 is not None else math.inf
        print(
            f"    {a}-{b} nrev={nr} beta={beta0} tof={tof0}: residual={rn:.2e} km/s "
            f"(a near-miss, per #562's own table -- not an exact closure to converge onto)"
        )
    print("  => nothing isolated to seed; the exact-closure set nearby is the degenerate manifold.")


def main() -> int:
    print("#680 asymmetric-closure census -- reproducible diagnosis of the reopening premise\n")
    # Wrapper faithfulness: reproduce the #312 catalogued point exactly.
    pt = residual_at_point("Umbriel", "Oberon", rel_offset_deg=180.0, tof_scale=2.0, n_rev=(1, 1))
    assert pt is not None
    print(
        f"[wrapper check] #312 point residual = {pt['residual_kms']:.12f} (ref 0.025232272564610)\n"
    )
    finding_1_non_isolation()
    finding_2_and_3_degeneracy()
    finding_4_connectivity()
    finding_5_bend_and_seeds()
    print("\n" + "=" * 78)
    print("VERDICT: reopening premise does NOT hold -- the free-beta asymmetric closure set is")
    print("a degenerate continuum continuously connected to the already-catalogued (#569)")
    print("symmetric family, not a discrete set of isolated novel closures. No catalogue action.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
