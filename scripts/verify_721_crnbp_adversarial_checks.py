"""#721: adversarial, independent verification of #717's N=5 CRNBP EOM and #720's torus.

This script deliberately does NOT reuse `#717`'s own derivation/tests or `#720`'s
own pseudospectral residual as evidence. It provides:

``eom``
    (1) A from-raw-Newton re-derivation check of the zero-coupling claim:
    build the EXACT relative acceleration of a massless particle w.r.t. the
    M1-M2 barycenter directly from Newton's law (spacecraft acceleration minus
    the barycenter's own acceleration ``(M1*a1 + M2*a2)/(M1+M2)`` -- a route
    on which extra-body-to-extra-body coupling terms NEVER appear), and compare
    it, at random configurations, against (a) Negri & Prado's Eq. 8 bracket
    form WITH the inner-sum coupling term and (b) the same with the
    extra-extra coupling omitted (superposition). All three must agree to
    machine precision, for 2 AND 3 extra bodies.
    (2) A faithful-Eq.-9/10/11 rotating-frame acceleration (WITH the ``-mu2``
    barycentric shift of the perturber positions and the full k=1,2 indirect
    terms) compared against ``core.crnbp.crnbp_eom``, to QUANTIFY the module's
    inherited idealisations (unshifted perturber positions, Jupiter-only
    indirect term) rather than assume them small.
    (3) The Laplace-configuration check: the physical libration center is
    ``theta_io + 2*theta_gan = 180 deg`` (lambda_Io - 3*lambda_Eu + 2*lambda_Gan
    = 180 deg, librating with ~0.064 deg amplitude -- Sinclair 1975 / Paita et
    al. 2018); report what the in-repo default actually uses.

``build-seed`` / ``continue-mu-io``
    Rebuild `#720`'s own pipeline (3:4 resonant orbit -> #690 CCR4BP torus at
    n1=2 -> mu_io=0 CRNBP seed -> mu_Io continuation to the physical value),
    checkpointing each stage to ``.npz`` so no single run exceeds a shell
    timeout. ``--theta-io0`` allows re-running the continuation at the
    PHYSICAL Laplace phase (pi) instead of the repo default (0).

``flow``
    The genuinely independent convergence check (`ghost_guard`-style): sample
    random torus angles, propagate the torus state through the FULL
    coupling-included ``crnbp_eom`` with BOTH DOP853 and Radau at tight
    tolerances over multiple horizons (up to a full forcing period -- 50x the
    corrector's own 0.02-period closure check), and measure the distance from
    the flowed state to the torus's own predicted point
    ``u(theta1 + omega1*t, theta2 + omega2*t)``, plus the Radau-vs-DOP853
    endpoint delta (integrator-artifact bound). Also an independent
    invariance-PDE residual via CENTRAL FINITE DIFFERENCES of the evaluated
    torus surface (no reuse of the corrector's basis-derivative machinery).
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.ccr4bp as ccr4bp
import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.crnbp as crnbp
import cyclerfinder.search.variational_ccr4bp_torus as vt
import cyclerfinder.search.variational_crnbp_torus as vc
from cyclerfinder.genome.composed_moon_map import resonance_semimajor

Vec = NDArray[np.float64]


# ---------------------------------------------------------------------------
# Part 1: from-raw-Newton re-derivation of the zero-coupling claim.
# ---------------------------------------------------------------------------


def _accel_on(target: Vec, sources: list[tuple[float, Vec]]) -> Vec:
    """Newtonian acceleration on a body at ``target`` from point masses (G=1)."""
    acc = np.zeros(3)
    for gm, pos in sources:
        d = pos - target
        acc += gm * d / float(np.linalg.norm(d)) ** 3
    return acc


def _raw_newton_relative_accel(
    m1: float, m2: float, p1: Vec, p2: Vec, extras: list[tuple[float, Vec]], p_sc: Vec
) -> Vec:
    """EXACT relative acceleration of the massless particle w.r.t. the M1-M2
    barycenter, from raw Newton only.

    ``a_rel = a_sc - (M1*a1 + M2*a2)/(M1+M2)`` where ``a1``/``a2`` are the
    primaries' own Newtonian accelerations. On THIS route the mutual
    extra-extra coupling NEVER appears (each primary's acceleration is linear
    in each extra body's mass) -- so agreement with the Eq. 8 bracket form is
    an independent proof that Eq. 8's coupling term nets to exactly zero.
    """
    all_for_sc = [(m1, p1), (m2, p2), *extras]
    a_sc = _accel_on(p_sc, all_for_sc)
    a_1 = _accel_on(p1, [(m2, p2), *extras])
    a_2 = _accel_on(p2, [(m1, p1), *extras])
    a_bary = (m1 * a_1 + m2 * a_2) / (m1 + m2)
    return a_sc - a_bary


def _eq8_relative_accel(
    m1: float,
    m2: float,
    p1: Vec,
    p2: Vec,
    extras: list[tuple[float, Vec]],
    p_sc: Vec,
    *,
    include_coupling: bool,
) -> Vec:
    """Negri & Prado (2022) Eq. 8, transcribed directly from the source PDF.

    ``include_coupling=False`` keeps the k=1,2 (primary/indirect) part of the
    inner sum but drops the extra-extra (k>=3) coupling terms -- the
    "superposition" form `#717` claims is mathematically identical.
    """
    acc = np.zeros(3)
    for gm, pos in ((m1, p1), (m2, p2)):
        d = pos - p_sc
        acc += gm * d / float(np.linalg.norm(d)) ** 3
    for j, (gm_j, pos_j) in enumerate(extras):
        d = pos_j - p_sc
        bracket = d / float(np.linalg.norm(d)) ** 3
        inner = np.zeros(3)
        for gm_k, pos_k in ((m1, p1), (m2, p2)):
            dk = pos_k - pos_j
            inner += (gm_k / (m1 + m2)) * dk / float(np.linalg.norm(dk)) ** 3
        if include_coupling:
            for k, (gm_k, pos_k) in enumerate(extras):
                if k == j:
                    continue
                dk = pos_k - pos_j
                inner += (gm_k / (m1 + m2)) * dk / float(np.linalg.norm(dk)) ** 3
        acc += gm_j * (bracket + inner)
    return acc


def check_zero_coupling_from_raw_newton(n_extras: int, n_trials: int, seed: int) -> float:
    """Max deviation between raw-Newton, Eq.-8-with-coupling and superposition."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(n_trials):
        m1 = 1.0
        m2 = float(rng.uniform(1e-5, 0.3))
        m_ex = [float(rng.uniform(1e-6, 1e-2)) for _ in range(n_extras)]
        p1 = rng.normal(size=3)
        p2 = rng.normal(size=3) * 2.0
        extras = [(m, rng.normal(size=3) * 3.0) for m in m_ex]
        p_sc = rng.normal(size=3) * 1.5
        a_raw = _raw_newton_relative_accel(m1, m2, p1, p2, extras, p_sc)
        a_eq8 = _eq8_relative_accel(m1, m2, p1, p2, extras, p_sc, include_coupling=True)
        a_sup = _eq8_relative_accel(m1, m2, p1, p2, extras, p_sc, include_coupling=False)
        scale = max(1.0, float(np.max(np.abs(a_raw))))
        worst = max(
            worst,
            float(np.max(np.abs(a_eq8 - a_raw))) / scale,
            float(np.max(np.abs(a_sup - a_raw))) / scale,
            float(np.max(np.abs(a_sup - a_eq8))) / scale,
        )
    return worst


def faithful_eq11_accel(state6: Vec, t: float, system: crnbp.CRNBPSystem) -> Vec:
    """FAITHFUL Negri-Prado Eq. 9+10 rotating-frame acceleration (independent
    transcription): perturbers at their BARYCENTRIC positions (``-mu2`` shift
    included), full k=1,2 indirect terms with the correct ``mu1``/``mu2``
    coefficients, plus the (provably net-zero) extra-extra coupling.

    Used to QUANTIFY ``crnbp_eom``'s inherited idealisations, not to assert
    exact agreement.
    """
    mu2 = system.mu
    mu1 = 1.0 - mu2
    x, y, z, vx, vy, _vz = (float(v) for v in state6)
    rho = np.array([x, y, z])
    p1 = np.array([-mu2, 0.0, 0.0])
    p2 = np.array([1.0 - mu2, 0.0, 0.0])
    extras: list[tuple[float, Vec]] = []
    for p in system.perturbers:
        gx, gy = p.position(t)
        extras.append((p.mu, np.array([gx - mu2, gy, 0.0])))
    grav = _eq8_relative_accel(mu1, mu2, p1, p2, extras, rho, include_coupling=True)
    # Rotating-frame terms (Omega = z-hat, nondim): -2 Omega x v - Omega x (Omega x rho)
    ax = 2.0 * vy + x + grav[0]
    ay = -2.0 * vx + y + grav[1]
    az = grav[2]
    return np.array([ax, ay, az])


def run_eom_checks() -> None:
    print("== #721 independent EOM checks ==")
    for n_extras in (2, 3):
        worst = check_zero_coupling_from_raw_newton(n_extras, n_trials=200, seed=42 + n_extras)
        print(
            f"raw-Newton vs Eq.8(+coupling) vs superposition, {n_extras} extras, "
            f"200 random configs: max rel deviation = {worst:.3e}"
        )
        assert worst < 1e-12, "zero-coupling claim FAILED the raw-Newton check"

    sys5 = crnbp.jupiter_europa_io_ganymede_default()
    rng = np.random.default_rng(7)
    worst_ideal = 0.0
    io_scale = 0.0
    for _ in range(50):
        s = np.array(
            [
                rng.uniform(0.4, 1.3),
                rng.uniform(-0.6, 0.6),
                rng.uniform(-0.05, 0.05),
                rng.uniform(-0.3, 0.3),
                rng.uniform(-0.3, 0.3),
                rng.uniform(-0.05, 0.05),
            ]
        )
        t = float(rng.uniform(0.0, 25.0))
        a_code = crnbp.crnbp_eom(t, s, sys5)[3:6]
        a_faithful = faithful_eq11_accel(s, t, sys5)
        worst_ideal = max(worst_ideal, float(np.max(np.abs(a_code - a_faithful))))
        # Io's own direct+indirect forcing magnitude, for scale.
        io = sys5.perturbers[0]
        gx, gy = io.position(t)
        d = math.hypot(s[0] - gx, s[1] - gy)
        io_scale = max(io_scale, io.mu / d**2 + io.mu / io.a**2)
    print(
        f"crnbp_eom vs faithful Eq.9/10 (barycentric shift + full indirect): "
        f"max |delta a| = {worst_ideal:.3e} (Io forcing scale ~{io_scale:.3e})"
    )

    io, gan = sys5.perturbers
    arg_deg = math.degrees(io.theta0 + 2.0 * gan.theta0) % 360.0
    print(
        f"Laplace argument theta_io0 + 2*theta_gan0 of the repo default = {arg_deg:.1f} deg "
        f"(physical libration center = 180.0 deg, amplitude ~0.064 deg)"
    )
    dt_arg = abs(io.omega + 2.0 * gan.omega)
    print(f"d/dt of the Laplace argument under the projected rates = {dt_arg:.3e} (0 = locked)")


# ---------------------------------------------------------------------------
# Part 2: rebuild #720's pipeline with stage checkpoints.
# ---------------------------------------------------------------------------


def _resonant_symmetric_orbit(
    mu: float, p_sc: int, q_moon: int, *, max_iter: int = 60, tol: float = 1e-12, cap: float = 0.05
) -> tuple[Vec, float, float]:
    """Symmetric p:q resonant CR3BP periodic orbit (same construction as
    `#690`'s / `#720`'s own test scaffolding, re-typed here)."""
    a = resonance_semimajor(p_sc, q_moon)
    period = 2.0 * np.pi * q_moon
    th = 0.5 * period
    x0 = a - mu
    vy0 = float(np.sqrt((1.0 - mu) / a)) - x0
    res = np.inf
    for k in range(max_iter):
        s0 = np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0])
        y42 = np.concatenate([s0, np.eye(6).reshape(36)])
        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom, (0.0, th), y42, args=(mu,), method="DOP853", rtol=1e-12, atol=1e-12
        )
        sf = sol.y[:, -1]
        phi = sf[6:].reshape(6, 6)
        g = np.array([sf[1], sf[3]])
        res = float(np.linalg.norm(g))
        if res < tol:
            break
        jac = np.array([[phi[1, 0], phi[1, 4]], [phi[3, 0], phi[3, 4]]])
        dz = np.linalg.solve(jac, -g)
        dz = dz * (0.3 if k < 8 else 1.0)
        norm = float(np.linalg.norm(dz))
        if norm > cap:
            dz = dz / norm * cap
        x0 += dz[0]
        vy0 += dz[1]
    return np.array([x0, 0.0, 0.0, 0.0, vy0, 0.0]), period, res


@dataclass(frozen=True)
class TorusRecord:
    """Minimal serializable snapshot of a CRNBP torus result."""

    mu: float
    perturbers: tuple[tuple[float, float, float, float], ...]  # (mu, a, omega, theta0)
    coeffs: Vec
    omega1: float
    omega2: float
    n1: int
    n2: int
    m1: int
    m2: int
    period_multiple: int
    transverse_amplitude: float
    residual_rms: float
    closure_residual: float


def _record_from_result(r: vc.CRNBPTorusVariationalResult) -> TorusRecord:
    perts = tuple((p.mu, p.a, p.omega, p.theta0) for p in r.system.perturbers)
    return TorusRecord(
        mu=r.system.mu,
        perturbers=perts,
        coeffs=r.coeffs,
        omega1=r.omega1,
        omega2=r.omega2,
        n1=r.n1,
        n2=r.n2,
        m1=r.m1,
        m2=r.m2,
        period_multiple=r.period_multiple,
        transverse_amplitude=r.transverse_amplitude,
        residual_rms=r.residual_rms,
        closure_residual=r.closure_residual,
    )


def _save_record(path: Path, rec: TorusRecord) -> None:
    np.savez(
        path,
        mu=rec.mu,
        perturbers=np.array(rec.perturbers),
        coeffs=rec.coeffs,
        omega1=rec.omega1,
        omega2=rec.omega2,
        n1=rec.n1,
        n2=rec.n2,
        m1=rec.m1,
        m2=rec.m2,
        period_multiple=rec.period_multiple,
        transverse_amplitude=rec.transverse_amplitude,
        residual_rms=rec.residual_rms,
        closure_residual=rec.closure_residual,
    )


def _load_record(path: Path) -> TorusRecord:
    z = np.load(path)
    perts = tuple(
        (float(row[0]), float(row[1]), float(row[2]), float(row[3])) for row in z["perturbers"]
    )
    return TorusRecord(
        mu=float(z["mu"]),
        perturbers=perts,
        coeffs=z["coeffs"],
        omega1=float(z["omega1"]),
        omega2=float(z["omega2"]),
        n1=int(z["n1"]),
        n2=int(z["n2"]),
        m1=int(z["m1"]),
        m2=int(z["m2"]),
        period_multiple=int(z["period_multiple"]),
        transverse_amplitude=float(z["transverse_amplitude"]),
        residual_rms=float(z["residual_rms"]),
        closure_residual=float(z["closure_residual"]),
    )


def _result_from_record(rec: TorusRecord) -> vc.CRNBPTorusVariationalResult:
    perts = tuple(
        crnbp.CRNBPPerturber(mu=m, a=a, omega=om, theta0=th) for (m, a, om, th) in rec.perturbers
    )
    system = crnbp.CRNBPSystem(mu=rec.mu, perturbers=perts)
    return vc.CRNBPTorusVariationalResult(
        system=system,
        coeffs=rec.coeffs,
        omega1=rec.omega1,
        omega2=rec.omega2,
        rotation_number=rec.omega2 / rec.omega1,
        rho_strob=rec.omega2 * (2.0 * np.pi / rec.omega1),
        period=2.0 * np.pi / rec.omega1,
        n1=rec.n1,
        n2=rec.n2,
        m1=rec.m1,
        m2=rec.m2,
        period_multiple=rec.period_multiple,
        transverse_amplitude=rec.transverse_amplitude,
        residual_rms=rec.residual_rms,
        closure_residual=rec.closure_residual,
        converged=True,
        n_iter=0,
    )


def run_build_seed(workdir: Path, theta_io0: float) -> None:
    workdir.mkdir(parents=True, exist_ok=True)
    sys4 = ccr4bp.jupiter_europa_ganymede_default()
    s0, period, res = _resonant_symmetric_orbit(sys4.mu, 3, 4)
    print(f"3:4 resonant orbit: perp residual {res:.2e}")
    assert res < 1e-10
    phys = vt.discover_ccr4bp_torus_from_resonant_orbit(
        sys4,
        s0,
        period,
        n1=2,
        n2=20,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    print(
        f"CCR4BP n1=2 torus: residual_rms={phys.residual_rms:.6e} "
        f"closure={phys.closure_residual:.6e} rot={phys.rotation_number:.9f}"
    )
    target = crnbp.jupiter_europa_io_ganymede_default()
    io = target.perturbers[0]
    seed = vc.discover_crnbp_torus_from_ccr4bp_seed(
        phys,
        mu_io=0.0,
        a_io=io.a,
        omega_io=io.omega,
        theta_io0=theta_io0,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    print(
        f"CRNBP mu_io=0 seed (theta_io0={theta_io0:.6f}): "
        f"residual_rms={seed.residual_rms:.6e} closure={seed.closure_residual:.6e}"
    )
    _save_record(workdir / "seed.npz", _record_from_result(seed))
    print(f"saved {workdir / 'seed.npz'}")


def run_continue_mu_io(workdir: Path, theta_io0: float, out_name: str) -> None:
    seed = _result_from_record(_load_record(workdir / "seed.npz"))
    target0 = crnbp.jupiter_europa_io_ganymede_default()
    io_t = target0.perturbers[0]
    gan_t = target0.perturbers[1]
    # Re-anchor Io's epoch phase (the seed itself was built at this theta_io0
    # already; at mu_io=0 the phase is inert, so re-stamping is safe).
    io_t = crnbp.CRNBPPerturber(mu=io_t.mu, a=io_t.a, omega=io_t.omega, theta0=theta_io0)
    target = crnbp.CRNBPSystem(mu=target0.mu, perturbers=(io_t, gan_t))
    steps = vc.continue_crnbp_torus_mu_io(
        seed,
        target,
        n_steps=8,
        tr_solver="exact",
        max_nfev=600,
        gauge_weight=30.0,
        rho_weight=100.0,
    )
    print(f"continuation (theta_io0={theta_io0:.6f}): {len(steps)} steps")
    for st in steps:
        print(
            f"  mu_io={st.system.perturbers[0].mu:.6e} residual_rms={st.residual_rms:.6e} "
            f"closure={st.closure_residual:.6e} rot={st.rotation_number:.9f}"
        )
    final = steps[-1]
    assert final.system.perturbers[0].mu == io_t.mu, "continuation did not reach physical mu_io"
    _save_record(workdir / out_name, _record_from_result(final))
    print(f"saved {workdir / out_name}")


# ---------------------------------------------------------------------------
# Part 3: independent flow-invariance check (ghost_guard-style, Radau+DOP853).
# ---------------------------------------------------------------------------


def _fd_invariance_residual(result: vc.CRNBPTorusVariationalResult, n_samples: int) -> float:
    """Invariance-PDE residual via CENTRAL FINITE DIFFERENCES of the evaluated
    torus surface at random angles (independent of the corrector's own basis
    derivative matrices and its collocation grid), against the FULL
    coupling-included scalar ``crnbp_eom``. Returns the RMS over samples."""
    rng = np.random.default_rng(2026)
    h = 1e-5
    sq = 0.0
    for _ in range(n_samples):
        th1 = float(rng.uniform(0, 2 * np.pi))
        th2 = float(rng.uniform(0, 2 * np.pi))
        u = vc.evaluate_torus_state(result, th1, th2)
        du1 = (
            vc.evaluate_torus_state(result, th1 + h, th2)
            - vc.evaluate_torus_state(result, th1 - h, th2)
        ) / (2 * h)
        du2 = (
            vc.evaluate_torus_state(result, th1, th2 + h)
            - vc.evaluate_torus_state(result, th1, th2 - h)
        ) / (2 * h)
        s6 = np.array([u[0], u[1], 0.0, u[2], u[3], 0.0])
        t0 = th1 / result.omega1
        f6 = crnbp.crnbp_eom(t0, s6, result.system)
        f4 = np.array([f6[0], f6[1], f6[3], f6[4]])
        r = result.omega1 * du1 + result.omega2 * du2 - f4
        sq += float(np.sum(r * r)) / 4.0
    return math.sqrt(sq / n_samples)


def run_flow_check(workdir: Path, torus_file: str, n_samples: int) -> None:
    result = _result_from_record(_load_record(workdir / torus_file))
    print(
        f"== flow check on {torus_file}: mu_io={result.system.perturbers[0].mu:.6e} "
        f"theta_io0={result.system.perturbers[0].theta0:.6f} "
        f"stored residual_rms={result.residual_rms:.6e} =="
    )
    fd_rms = _fd_invariance_residual(result, 200)
    print(f"independent FD invariance residual (200 random pts): rms={fd_rms:.6e}")

    period = result.period
    horizons = [0.05 * period, 0.25 * period, 1.0 * period]
    rng = np.random.default_rng(0xBEEF)
    rows: list[tuple[float, float, list[float], list[float], float]] = []
    for _ in range(n_samples):
        th1 = float(rng.uniform(0, 2 * np.pi))
        th2 = float(rng.uniform(0, 2 * np.pi))
        u0 = vc.evaluate_torus_state(result, th1, th2)
        s6 = np.array([u0[0], u0[1], 0.0, u0[2], u0[3], 0.0])
        t0 = th1 / result.omega1
        errs: dict[str, list[float]] = {}
        finals: dict[str, Vec] = {}
        for method in ("DOP853", "Radau"):
            sol = solve_ivp(
                crnbp.crnbp_eom,
                (t0, t0 + horizons[-1]),
                s6,
                args=(result.system,),
                method=method,
                t_eval=[t0 + h for h in horizons],
                rtol=1e-12,
                atol=1e-12,
            )
            assert sol.success, f"{method} propagation failed"
            errs[method] = []
            for col, h in zip(sol.y.T, horizons, strict=True):
                target = vc.evaluate_torus_state(
                    result, th1 + result.omega1 * h, th2 + result.omega2 * h
                )
                planar = np.array([col[0], col[1], col[3], col[4]])
                errs[method].append(float(np.linalg.norm(planar - target)))
            finals[method] = np.asarray(sol.y[:, -1])
        integ_delta = float(np.linalg.norm(finals["DOP853"] - finals["Radau"]))
        rows.append((th1, th2, errs["DOP853"], errs["Radau"], integ_delta))

    for label, idx in (("0.05P", 0), ("0.25P", 1), ("1.00P", 2)):
        d_dop = np.array([r[2][idx] for r in rows])
        d_rad = np.array([r[3][idx] for r in rows])
        print(
            f"horizon {label}: off-torus |err| DOP853 median={np.median(d_dop):.3e} "
            f"max={np.max(d_dop):.3e} | Radau median={np.median(d_rad):.3e} "
            f"max={np.max(d_rad):.3e}"
        )
    deltas = np.array([r[4] for r in rows])
    print(
        f"Radau-vs-DOP853 endpoint delta at 1.00P: median={np.median(deltas):.3e} "
        f"max={np.max(deltas):.3e} (nondim; 1 DU = 671100 km)"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("eom")
    p_seed = sub.add_parser("build-seed")
    p_seed.add_argument("--workdir", type=Path, required=True)
    p_seed.add_argument("--theta-io0", type=float, default=0.0)
    p_cont = sub.add_parser("continue-mu-io")
    p_cont.add_argument("--workdir", type=Path, required=True)
    p_cont.add_argument("--theta-io0", type=float, default=0.0)
    p_cont.add_argument("--out-name", type=str, default="final.npz")
    p_flow = sub.add_parser("flow")
    p_flow.add_argument("--workdir", type=Path, required=True)
    p_flow.add_argument("--torus-file", type=str, default="final.npz")
    p_flow.add_argument("--n-samples", type=int, default=8)
    args = ap.parse_args()
    if args.cmd == "eom":
        run_eom_checks()
    elif args.cmd == "build-seed":
        run_build_seed(args.workdir, args.theta_io0)
    elif args.cmd == "continue-mu-io":
        run_continue_mu_io(args.workdir, args.theta_io0, args.out_name)
    elif args.cmd == "flow":
        run_flow_check(args.workdir, args.torus_file, args.n_samples)


if __name__ == "__main__":
    main()
