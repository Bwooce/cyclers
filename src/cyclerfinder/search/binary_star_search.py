"""Direct fixed-mu binary-star cycler search (#255).

Roberts-Tsoukkas & Ross 2026 (journal: "Stable Prograde Earth-Moon Multi-Orbiter
Cyclers") depict -- in FIGURES ONLY (Fig. 3), no printed numbers -- stable
prograde cyclers persisting into the binary-star mass-parameter range: a (1,3)
exterior cycler at mu=0.1, a (3,1) cycler at mu=0.3, and a (1,1) equal-mass
cycler at mu=0.5.

The #252 pseudo-arclength *continuation* of the held Earth-Moon members FAILED to
reach these: continued to the target mu as genuine periodic orbits, but
branch-switched at folds and lost the cycler topology (the depicted orbit lives
on a different (x0, C) branch not reached from the EM seed). This module takes
the complementary route the #252 note recommends: **figure-read the approximate
x0 directly off Fig. 3 and seed the fixed-mu corrector in the cycler basin**, so
the corrector converges onto the cycler branch rather than the one-sided
librational branch the continuation fell into.

DISCOVERY DISCIPLINE: any orbit found here is OUR OWN computed
``(mu, C, T, IC, nu)`` -- the paper prints none for the binary-star members. A
match is a **discovery candidate** requiring the literature-novelty check (#261)
and the full V0-V5 gauntlet, NOT a sourced row. This module emits review-gated
evidence only; it never writes the catalogue.

Topology is classified by the **winding number around each primary** (k1 = signed
revolutions about P1 at ``(-mu, 0)``, k2 = about P2 at ``(1-mu, 0)``), which
reproduces the published Earth-Moon (3,1) and (1,1) labels exactly and is
prograde (positive) -- see :func:`tests.search.test_binary_star_search`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp


def collinear_lpoints(mu: float) -> tuple[float, float, float]:
    """Return the collinear libration points ``(L1, L2, L3)`` on the x-axis.

    Roots of ``dOmega/dx = 0`` at ``y = 0``: ``L1`` between the primaries,
    ``L2`` beyond the secondary, ``L3`` beyond the primary.
    """

    def dom(x: float) -> float:
        r1 = x + mu
        r2 = x - (1.0 - mu)
        return x - (1.0 - mu) * r1 / abs(r1) ** 3 - mu * r2 / abs(r2) ** 3

    l1 = brentq(dom, -mu + 1e-6, 1.0 - mu - 1e-6)
    l2 = brentq(dom, 1.0 - mu + 1e-6, 2.5)
    l3 = brentq(dom, -2.5, -mu - 1e-6)
    return float(l1), float(l2), float(l3)


@dataclass(frozen=True)
class Topology:
    """Winding-number topology of one full-period planar orbit."""

    k1: int  # rounded |winding| about the primary P1 at (-mu, 0)
    k2: int  # rounded |winding| about the secondary P2 at (1-mu, 0)
    w1: float  # raw signed winding about P1 (revolutions)
    w2: float  # raw signed winding about P2
    prograde: bool  # both windings positive (paper's prograde cyclers)
    x_min: float
    x_max: float
    reaches_secondary: bool  # x-extent passes the L1 neck into the secondary realm


def winding_topology(
    mu: float,
    state0: NDArray[np.float64],
    period: float,
    *,
    rtol: float = 1e-11,
    atol: float = 1e-11,
    n_samples: int = 4000,
) -> Topology:
    """Classify a planar CR3BP periodic orbit by winding number about each primary."""
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        np.asarray(state0, float),
        args=(mu,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        max_step=period / n_samples,
    )
    x, y = sol.y[0], sol.y[1]

    def wind(px: float, py: float) -> float:
        th = np.unwrap(np.arctan2(y - py, x - px))
        return float((th[-1] - th[0]) / (2.0 * np.pi))

    w1 = wind(-mu, 0.0)
    w2 = wind(1.0 - mu, 0.0)
    l1, _l2, _l3 = collinear_lpoints(mu)
    return Topology(
        k1=round(abs(w1)),
        k2=round(abs(w2)),
        w1=w1,
        w2=w2,
        prograde=(w1 > 0.0 and w2 > 0.0),
        x_min=float(x.min()),
        x_max=float(x.max()),
        reaches_secondary=bool(x.max() > l1),
    )


@dataclass(frozen=True)
class BinaryStarCandidate:
    mu: float
    x0: float
    ydot0: float
    jacobi: float
    period: float
    nu: float
    abs_lambda: float
    k1: int
    k2: int
    prograde: bool
    x_min: float
    x_max: float
    crossing_residual: float
    radau_djacobi: float
    stable: bool


def _system(mu: float) -> cr3bp.CR3BPSystem:
    # l_km / t_s are the Earth-Moon scale; they do not enter the nondimensional
    # corrector or the winding classifier (only Jacobi/period are reported in TU).
    return cr3bp.CR3BPSystem(
        mu=float(mu), primary="Earth", secondary="Moon", l_km=384400.0, t_s=375699.8
    )


def figure_seeded_search(
    mu: float,
    target_k1: int,
    target_k2: int,
    *,
    x0_lo: float,
    x0_hi: float,
    c_lo: float,
    c_hi: float,
    n_x0: int = 24,
    n_c: int = 16,
    half_crossings_set: tuple[int, ...] = (1, 2, 3, 4, 5),
    ydot0_signs: tuple[float, ...] = (-1.0, 1.0),
    period_guess: float = 14.0,
    corrector_tol: float = 1e-12,
    radau_closure_tol: float = 1e-3,
    radau_jacobi_tol: float = 1e-8,
    require_stable: bool = True,
    require_prograde: bool = True,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> list[BinaryStarCandidate]:
    """Seed the fixed-mu symmetric corrector across the figure-read ``(x0, C)``
    box and keep converged orbits whose winding topology matches the depicted
    ``(target_k1, target_k2)`` cycler.

    The cheap winding/topology screen runs first; the expensive Barden monodromy
    and the independent-Radau Jacobi cross-check run only on topology survivors.
    Returns deduplicated candidates (by rounded ``(x0, C, T)``).
    """
    system = _system(mu)
    x0_grid = np.linspace(x0_lo, x0_hi, n_x0)
    c_grid = np.linspace(c_lo, c_hi, n_c)
    seen: set[tuple[float, float, float]] = set()
    out: list[BinaryStarCandidate] = []

    for hc in half_crossings_set:
        for sign in ydot0_signs:
            for c in c_grid:
                for x0g in x0_grid:
                    try:
                        o = cp.correct_symmetric_fixed_jacobi(
                            system,
                            float(x0g),
                            float(c),
                            period_guess,
                            ydot0_sign=sign,
                            half_crossings=hc,
                            tol=corrector_tol,
                            max_iter=40,
                            rtol=rtol,
                            atol=atol,
                        )
                    except (ValueError, RuntimeError):
                        continue
                    if not o.converged:
                        continue
                    key = (round(o.x0, 4), round(o.jacobi, 4), round(o.period, 3))
                    if key in seen:
                        continue
                    seen.add(key)

                    state0 = np.array([o.x0, 0.0, 0.0, 0.0, o.ydot0, 0.0])
                    topo = winding_topology(mu, state0, o.period)
                    if topo.k1 != target_k1 or topo.k2 != target_k2:
                        continue
                    if require_prograde and not topo.prograde:
                        continue

                    try:
                        nu, lam = cp.barden_stability(system, o, rtol=rtol, atol=atol)
                    except Exception:
                        continue
                    stable = abs(float(nu)) < 1.0
                    if require_stable and not stable:
                        continue

                    po = cp.PeriodicOrbit(
                        state0=state0,
                        period=o.period,
                        jacobi=o.jacobi,
                        converged=True,
                        closure_residual=o.crossing_residual,
                    )
                    _ok, radau_dj = cp.crosscheck_periodic(
                        system, po, closure_tol=radau_closure_tol, jacobi_tol=radau_jacobi_tol
                    )
                    out.append(
                        BinaryStarCandidate(
                            mu=float(mu),
                            x0=float(o.x0),
                            ydot0=float(o.ydot0),
                            jacobi=float(o.jacobi),
                            period=float(o.period),
                            nu=float(nu),
                            abs_lambda=float(abs(lam)),
                            k1=topo.k1,
                            k2=topo.k2,
                            prograde=topo.prograde,
                            x_min=topo.x_min,
                            x_max=topo.x_max,
                            crossing_residual=float(o.crossing_residual),
                            radau_djacobi=float(radau_dj),
                            stable=stable,
                        )
                    )
    return out
