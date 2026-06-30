"""Task #504: Pluto-Charon (k1,k2)-cycler family sweep utilities.

Sweeps (k1,k2) families at Pluto-Charon mu=cr3bp_system("Pluto","Charon").mu,
seeding from the nearest Ross-RT 2026 Table-I anchor and locating the |nu|<1
stable window (brentq on the Barden nu(C)=0 midpoint).

Public API
----------
make_pluto_charon_system()     -> CR3BPSystem (PC physical scales)
mu_step_to_orbit(...)          -> SymmetricOrbit | None   (mu-continuation)
c_sweep_find_nu_zero(...)      -> SymmetricOrbit | None   (C-sweep + brentq)
sweep_32_positive_control()    -> SweepResult
sweep_11()                     -> SweepResult
sweep_21()                     -> SweepResult
sweep_22()                     -> SweepResult
sweep_31()                     -> SweepResult
sweep_33()                     -> SweepResult

Ross-RT 2026 Table-I anchors:
  data/golden/ross_rt_2026_cycler_families.yaml
"""

from __future__ import annotations

import signal
from collections.abc import Callable
from dataclasses import dataclass
from typing import TypeVar

import numpy as np
from scipy.optimize import brentq

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.core.cr3bp import cr3bp_system
from cyclerfinder.search.binary_star_search import collinear_lpoints, winding_topology

# ---------------------------------------------------------------------------
# Per-call SIGALRM timeout helper
# ---------------------------------------------------------------------------

_T = TypeVar("_T")


class _CallTimeoutError(Exception):
    """Raised when a per-call SIGALRM fires."""


def _run_with_timeout(fn: Callable[[], _T], seconds: int = 3) -> _T | None:
    """Run ``fn()`` with a SIGALRM wall-clock limit.

    Returns ``None`` if the call exceeds ``seconds``.  Unix-only (SIGALRM).
    """

    def _handler(signum: int, frame: object) -> None:
        raise _CallTimeoutError()

    old = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(seconds)
    try:
        return fn()
    except _CallTimeoutError:
        return None
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


# ---------------------------------------------------------------------------
# Pluto-Charon constants
# ---------------------------------------------------------------------------

#: Mass-parameter from cr3bp_system("Pluto","Charon")  (GM_Charon / GM_system)
PC_MU: float = cr3bp_system("Pluto", "Charon").mu  # 0.10876473603280369

_PC_SYS_CACHED: cr3bp.CR3BPSystem | None = None


def make_pluto_charon_system() -> cr3bp.CR3BPSystem:
    """Return the CR3BPSystem for Pluto-Charon with physical scales."""
    global _PC_SYS_CACHED
    if _PC_SYS_CACHED is None:
        _PC_SYS_CACHED = cr3bp_system("Pluto", "Charon")
    return _PC_SYS_CACHED


def _c_l1(mu: float) -> float:
    """Jacobi constant at L1 for mass-ratio mu."""
    l1, _l2, _l3 = collinear_lpoints(mu)
    return float(cr3bp.jacobi_constant(np.array([l1, 0.0, 0.0, 0.0, 0.0, 0.0]), mu))


# ---------------------------------------------------------------------------
# Non-dimensional helper (for mu-stepping)
# ---------------------------------------------------------------------------


def _nd_system(mu: float) -> cr3bp.CR3BPSystem:
    """Non-dimensional CR3BP system (scales don't enter the corrector)."""
    return cr3bp.CR3BPSystem(mu=float(mu), primary="P1", secondary="P2", l_km=1.0, t_s=1.0)


# ---------------------------------------------------------------------------
# Mu-continuation
# ---------------------------------------------------------------------------


def mu_step_to_orbit(
    anchor_mu: float,
    target_mu: float,
    anchor_x0: float,
    anchor_jacobi: float,
    anchor_period: float,
    *,
    hc: int | None = None,
    sign: float = -1.0,
    n_steps: int = 40,
    tol: float = 1e-10,
    rtol: float = 1e-12,
    atol: float = 1e-12,
) -> cp.SymmetricOrbit | None:
    """Step mu from anchor_mu to target_mu in n_steps, returning the orbit at target_mu.

    Holds the Jacobi constant fixed (anchor_jacobi) and re-corrects x0/T at each
    mu step.  Returns None if convergence fails at any step.
    """
    mus = np.linspace(anchor_mu, target_mu, n_steps + 1)[1:]  # skip anchor itself
    x0_cur = anchor_x0
    t_cur = anchor_period
    jacobi = anchor_jacobi

    # Establish orbit at anchor_mu first
    sys_anc = _nd_system(anchor_mu)
    try:
        o = cp.correct_symmetric_fixed_jacobi(
            sys_anc,
            x0_cur,
            jacobi,
            t_cur,
            ydot0_sign=sign,
            half_crossings=hc,
            tol=tol,
            rtol=rtol,
            atol=atol,
        )
    except ValueError:
        return None
    if not o.converged:
        return None
    x0_cur, t_cur = o.x0, o.period

    for mu_next in mus:
        sys_next = _nd_system(mu_next)
        try:
            o = cp.correct_symmetric_fixed_jacobi(
                sys_next,
                x0_cur,
                jacobi,
                t_cur,
                ydot0_sign=sign,
                half_crossings=hc,
                tol=tol,
                rtol=rtol,
                atol=atol,
            )
        except ValueError:
            return None
        if not o.converged:
            return None
        x0_cur, t_cur = o.x0, o.period

    # Final correction in the PC system
    sys_pc = make_pluto_charon_system()
    try:
        o_final = cp.correct_symmetric_fixed_jacobi(
            sys_pc,
            x0_cur,
            jacobi,
            t_cur,
            ydot0_sign=sign,
            half_crossings=hc,
            tol=tol,
            rtol=rtol,
            atol=atol,
        )
    except ValueError:
        return None
    if not o_final.converged:
        return None
    return o_final


# ---------------------------------------------------------------------------
# C-sweep + brentq nu=0 finder
# ---------------------------------------------------------------------------


def c_sweep_find_nu_zero(
    system: cr3bp.CR3BPSystem,
    x0_start: float,
    jacobi_start: float,
    period_start: float,
    *,
    hc: int | None,
    sign: float = -1.0,
    c_lo: float,
    c_hi: float,
    n_coarse: int = 60,
    tol: float = 1e-10,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    nu_tol: float = 1e-10,
) -> cp.SymmetricOrbit | None:
    """C-sweep the (k1,k2) family branch and locate the nu=0 stable midpoint.

    Starts from (x0_start, jacobi_start, period_start), walks to c_lo if needed,
    then sweeps C from c_lo to c_hi in n_coarse steps computing Barden nu at
    each step.  Exits early at the first nu sign-change bracket and runs brentq
    to locate the nu=0 midpoint.

    Returns the SymmetricOrbit at the nu=0 midpoint, or None if no sign change
    is found (= no stable window on this branch in [c_lo, c_hi]).
    """
    # Walk from jacobi_start to c_lo if they differ
    x0_cur, t_cur = x0_start, period_start
    if jacobi_start != c_lo:
        for c_walk in np.linspace(jacobi_start, c_lo, 20)[1:]:
            try:
                o = cp.correct_symmetric_fixed_jacobi(
                    system,
                    x0_cur,
                    c_walk,
                    t_cur,
                    ydot0_sign=sign,
                    half_crossings=hc,
                    tol=tol,
                    rtol=rtol,
                    atol=atol,
                )
            except ValueError:
                break
            if o.converged:
                x0_cur, t_cur = o.x0, o.period

    # Coarse sweep: stop early on first sign change.
    c_grid = np.linspace(c_lo, c_hi, n_coarse)
    nu_prev: float | None = None
    orbit_prev: cp.SymmetricOrbit | None = None
    bracket: tuple[float, float, float, float] | None = None
    x0_sweep, t_sweep = x0_cur, t_cur

    for i, c_val in enumerate(c_grid):
        try:
            o = cp.correct_symmetric_fixed_jacobi(
                system,
                x0_sweep,
                c_val,
                t_sweep,
                ydot0_sign=sign,
                half_crossings=hc,
                tol=tol,
                rtol=rtol,
                atol=atol,
            )
        except ValueError:
            nu_prev = None
            orbit_prev = None
            continue
        if not o.converged:
            nu_prev = None
            orbit_prev = None
            continue
        nu, _ = cp.barden_stability(system, o, rtol=rtol, atol=atol)
        if nu_prev is not None and nu_prev * nu < 0.0:
            bracket = (c_grid[i - 1], nu_prev, c_val, nu)
            break  # early exit: first sign change found
        nu_prev = nu
        orbit_prev = o
        x0_sweep, t_sweep = o.x0, o.period

    if bracket is None:
        return None

    c_a, _nu_a, c_b, _nu_b = bracket

    # Brentq seed: orbit at c_a (last converged before sign change)
    assert orbit_prev is not None
    x0_brent = orbit_prev.x0
    t_brent = orbit_prev.period

    def _nu_at(c_val_inner: float) -> float:
        try:
            o = cp.correct_symmetric_fixed_jacobi(
                system,
                x0_brent,
                c_val_inner,
                t_brent,
                ydot0_sign=sign,
                half_crossings=hc,
                tol=1e-11,
                rtol=1e-13,
                atol=1e-13,
            )
        except ValueError:
            return float("nan")
        if not o.converged:
            return float("nan")
        nu_inner, _ = cp.barden_stability(system, o, rtol=1e-13, atol=1e-13)
        return float(nu_inner)

    try:
        c_mid = brentq(_nu_at, c_a, c_b, xtol=nu_tol, rtol=nu_tol, maxiter=60)
    except ValueError:
        return None

    try:
        o_mid = cp.correct_symmetric_fixed_jacobi(
            system,
            x0_brent,
            c_mid,
            t_brent,
            ydot0_sign=sign,
            half_crossings=hc,
            tol=1e-11,
            rtol=1e-13,
            atol=1e-13,
        )
    except ValueError:
        return None
    if not o_mid.converged:
        return None
    return o_mid


# ---------------------------------------------------------------------------
# Per-family sweep result
# ---------------------------------------------------------------------------


@dataclass
class SweepResult:
    """Result of sweeping one (k1,k2) family at Pluto-Charon."""

    k1: int
    k2: int
    stable_found: bool
    #
    jacobi_mid: float | None = None  # Jacobi constant at nu=0 midpoint
    x0_mid: float | None = None  # x0 at nu=0 midpoint
    ydot0_mid: float | None = None  # ydot0 at nu=0 midpoint
    period_mid: float | None = None  # period (TU) at nu=0 midpoint
    period_days: float | None = None  # period (days)
    nu_mid: float | None = None  # Barden nu at midpoint (should be ~0)
    topology_ok: bool = False  # winding matches (k1,k2)
    prograde: bool = False
    reaches_secondary: bool = False
    crosscheck_ok: bool = False
    crosscheck_dj: float | None = None
    #
    method: str = ""  # how the seed was obtained
    note: str = ""  # clean-negative reason or other note

    @property
    def c_mid(self) -> float | None:
        """Alias for jacobi_mid (backwards-compat)."""
        return self.jacobi_mid


# ---------------------------------------------------------------------------
# (3,2) positive control
# ---------------------------------------------------------------------------

# Anchor from Ross-RT 2026 Table I, row 4 (mu=0.1, (3,2))
_PC_ANCHOR_X0 = -0.694376003123377
_PC_ANCHOR_C = 3.573367616904619
_PC_ANCHOR_T = 12.295263874014290
_PC_ANCHOR_HC = 6


def _crosscheck(system: cr3bp.CR3BPSystem, orbit: cp.SymmetricOrbit) -> tuple[bool, float]:
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    po = cp.PeriodicOrbit(
        state0=state0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        converged=orbit.converged,
        closure_residual=orbit.crossing_residual,
    )
    return cp.crosscheck_periodic(system, po, closure_tol=1e-6, jacobi_tol=1e-8)


def _build_result(
    k1: int,
    k2: int,
    system: cr3bp.CR3BPSystem,
    orbit: cp.SymmetricOrbit,
    method: str,
) -> SweepResult:
    """Build a SweepResult from a converged SymmetricOrbit."""
    nu, _ = cp.barden_stability(system, orbit, rtol=1e-13, atol=1e-13)
    stable = abs(nu) < 1.0

    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    topo = winding_topology(system.mu, state0, orbit.period)
    topology_ok = topo.k1 == k1 and topo.k2 == k2

    ok_cc, dj = _crosscheck(system, orbit)

    period_days = orbit.period * system.t_s / 86400.0

    return SweepResult(
        k1=k1,
        k2=k2,
        stable_found=stable,
        jacobi_mid=float(orbit.jacobi),
        x0_mid=float(orbit.x0),
        ydot0_mid=float(orbit.ydot0),
        period_mid=float(orbit.period),
        period_days=float(period_days),
        nu_mid=float(nu),
        topology_ok=topology_ok,
        prograde=topo.prograde,
        reaches_secondary=topo.reaches_secondary,
        crosscheck_ok=ok_cc,
        crosscheck_dj=float(dj),
        method=method,
        note="",
    )


def sweep_32_positive_control() -> SweepResult:
    """Re-find the admitted (3,2) stable member at PC mu. Positive control."""
    sys_pc = make_pluto_charon_system()
    c_l1 = _c_l1(PC_MU)

    # Correct at anchor C (establishes the branch)
    orbit_anc = cp.correct_symmetric_fixed_jacobi(
        sys_pc,
        _PC_ANCHOR_X0,
        _PC_ANCHOR_C,
        _PC_ANCHOR_T,
        ydot0_sign=-1.0,
        half_crossings=_PC_ANCHOR_HC,
        tol=1e-10,
    )
    if not orbit_anc.converged:
        return SweepResult(
            k1=3,
            k2=2,
            stable_found=False,
            method="direct_seed",
            note="anchor correction failed",
        )

    # C-sweep upward from anchor to just below C_L1
    orbit_stable = c_sweep_find_nu_zero(
        sys_pc,
        orbit_anc.x0,
        orbit_anc.jacobi,
        orbit_anc.period,
        hc=_PC_ANCHOR_HC,
        sign=-1.0,
        c_lo=_PC_ANCHOR_C,
        c_hi=c_l1 - 0.002,
        n_coarse=60,
    )
    if orbit_stable is None:
        return SweepResult(
            k1=3,
            k2=2,
            stable_found=False,
            method="c_sweep",
            note="no stable window found",
        )

    return _build_result(3, 2, sys_pc, orbit_stable, "c_sweep_from_mu01_anchor")


# ---------------------------------------------------------------------------
# (1,1) sweep  — two Table-I anchors at mu=0.001 and mu=0.012150
# ---------------------------------------------------------------------------


def sweep_11() -> SweepResult:
    """Sweep the (1,1) family at Pluto-Charon mu.

    Tries mu-continuation from the nearest Table-I anchor (mu=0.001, k=(1,1))
    then C-sweeps for the stable window.  Falls back to the mu=0.012150 anchor.
    """
    sys_pc = make_pluto_charon_system()
    c_l1 = _c_l1(PC_MU)

    # Anchor 1: mu=0.001
    anc_mu = 0.001
    anc_x0 = -0.647047499999966
    anc_jacobi = 3.031605708907296
    anc_period = 14.774502790974823
    orbit_seed = mu_step_to_orbit(
        anc_mu,
        PC_MU,
        anc_x0,
        anc_jacobi,
        anc_period,
        hc=None,
        sign=-1.0,
        n_steps=40,
    )

    # Anchor 2 fallback: mu=0.012150
    if orbit_seed is None:
        anc_mu = 0.012150584270572
        anc_x0 = -0.768217354461248
        anc_jacobi = 3.151175879917331
        anc_period = 10.291893641936499
        orbit_seed = mu_step_to_orbit(
            anc_mu,
            PC_MU,
            anc_x0,
            anc_jacobi,
            anc_period,
            hc=None,
            sign=-1.0,
            n_steps=40,
        )

    if orbit_seed is None:
        return SweepResult(
            k1=1,
            k2=1,
            stable_found=False,
            method="mu_step",
            note="mu-continuation failed at both anchors",
        )

    # C-sweep around the landed jacobi
    c_lo = max(2.9, orbit_seed.jacobi - 0.3)
    c_hi = min(c_l1 - 0.002, orbit_seed.jacobi + 0.3)
    orbit_stable = c_sweep_find_nu_zero(
        sys_pc,
        orbit_seed.x0,
        orbit_seed.jacobi,
        orbit_seed.period,
        hc=None,
        sign=-1.0,
        c_lo=c_lo,
        c_hi=c_hi,
        n_coarse=60,
    )
    method = f"mu_step_from_mu{anc_mu:.3f}_then_c_sweep"
    if orbit_stable is None:
        return SweepResult(
            k1=1,
            k2=1,
            stable_found=False,
            method=method,
            note="no stable window in C-sweep range",
        )
    return _build_result(1, 1, sys_pc, orbit_stable, method)


# ---------------------------------------------------------------------------
# (3,1) sweep — Table-I anchor at mu=0.3
# ---------------------------------------------------------------------------


def sweep_31() -> SweepResult:
    """Sweep the (3,1) family at Pluto-Charon mu.

    The mu=0.3 anchor has C=3.702 > C_L1(PC)=3.621: above the Hill threshold
    at the target mu.  Mu-step from 0.3 to PC_MU at fixed C=3.702; the corrector
    will fail near where C crosses C_L1(mu_step).  Also tries C-walking the (3,1)
    family at mu=0.3 to a lower C before stepping mu.
    """
    sys_pc = make_pluto_charon_system()
    c_l1_pc = _c_l1(PC_MU)

    anc_mu = 0.3
    anc_x0 = -0.804725783387797
    anc_jacobi = 3.701958166478617
    anc_period = 9.094576400494693

    # Strategy A: mu-step from mu=0.3 at anchor jacobi
    orbit_seed = mu_step_to_orbit(
        anc_mu,
        PC_MU,
        anc_x0,
        anc_jacobi,
        anc_period,
        hc=None,
        sign=-1.0,
        n_steps=40,
    )

    if orbit_seed is not None and orbit_seed.jacobi < c_l1_pc - 0.001:
        c_lo = max(3.0, orbit_seed.jacobi - 0.3)
        c_hi = min(c_l1_pc - 0.002, orbit_seed.jacobi + 0.3)
        orbit_stable = c_sweep_find_nu_zero(
            sys_pc,
            orbit_seed.x0,
            orbit_seed.jacobi,
            orbit_seed.period,
            hc=None,
            sign=-1.0,
            c_lo=c_lo,
            c_hi=c_hi,
            n_coarse=60,
        )
        method = "mu_step_from_mu03_anchor"
        if orbit_stable is None:
            return SweepResult(
                k1=3,
                k2=1,
                stable_found=False,
                method=method,
                note="no stable window in C-sweep",
            )
        res = _build_result(3, 1, sys_pc, orbit_stable, method)
        if not res.topology_ok:
            return SweepResult(
                k1=3,
                k2=1,
                stable_found=False,
                method=method,
                note=(
                    "Strategy A found stable orbit but wrong topology "
                    f"(reaches_secondary={res.reaches_secondary}); clean negative"
                ),
            )
        return res

    # Strategy B: C-walk at mu=0.3 to a lower C, then mu-step
    c_l1_03 = _c_l1(0.3)
    c_target = min(c_l1_pc - 0.05, c_l1_03 - 0.01)
    if c_target < anc_jacobi:
        sys_03 = _nd_system(0.3)
        try:
            orbit_03 = cp.correct_symmetric_fixed_jacobi(
                sys_03, anc_x0, anc_jacobi, anc_period, ydot0_sign=-1.0, tol=1e-10
            )
        except ValueError:
            orbit_03 = None
        if orbit_03 is not None and orbit_03.converged:
            x0_walk, t_walk = orbit_03.x0, orbit_03.period
            for c_walk in np.linspace(anc_jacobi, c_target, 20)[1:]:
                try:
                    o = cp.correct_symmetric_fixed_jacobi(
                        sys_03, x0_walk, c_walk, t_walk, ydot0_sign=-1.0, tol=1e-10
                    )
                except ValueError:
                    break
                if not o.converged:
                    break
                x0_walk, t_walk = o.x0, o.period
            orbit_seed2 = mu_step_to_orbit(
                0.3, PC_MU, x0_walk, c_target, t_walk, hc=None, sign=-1.0, n_steps=40
            )
            if orbit_seed2 is not None and orbit_seed2.jacobi < c_l1_pc - 0.001:
                c_lo = max(3.0, orbit_seed2.jacobi - 0.3)
                c_hi = min(c_l1_pc - 0.002, orbit_seed2.jacobi + 0.3)
                orbit_stable = c_sweep_find_nu_zero(
                    sys_pc,
                    orbit_seed2.x0,
                    orbit_seed2.jacobi,
                    orbit_seed2.period,
                    hc=None,
                    sign=-1.0,
                    c_lo=c_lo,
                    c_hi=c_hi,
                    n_coarse=60,
                )
                method = "c_walk_at_mu03_then_mu_step"
                if orbit_stable is None:
                    return SweepResult(
                        k1=3,
                        k2=1,
                        stable_found=False,
                        method=method,
                        note="no stable window",
                    )
                res = _build_result(3, 1, sys_pc, orbit_stable, method)
                if not res.topology_ok:
                    return SweepResult(
                        k1=3,
                        k2=1,
                        stable_found=False,
                        method=method,
                        note=(
                            "Strategy B found a stable orbit but wrong topology "
                            f"(reaches_secondary={res.reaches_secondary}); clean negative"
                        ),
                    )
                return res

    return SweepResult(
        k1=3,
        k2=1,
        stable_found=False,
        method="mu_step_from_mu03",
        note=(
            f"anchor jacobi={anc_jacobi:.4f} > C_L1(PC)={c_l1_pc:.4f}; "
            "mu-continuation could not reach PC mu"
        ),
    )


# ---------------------------------------------------------------------------
# (3,3) sweep — Table-I anchor at mu=0.012150
# ---------------------------------------------------------------------------


def sweep_33() -> SweepResult:
    """Sweep the (3,3) family at Pluto-Charon mu.

    Mu-step from the mu=0.012150 anchor (k=(3,3)) to PC_MU in 40 steps,
    then C-sweep for the stable window.
    """
    sys_pc = make_pluto_charon_system()
    c_l1 = _c_l1(PC_MU)

    anc_mu = 0.012150584270572
    anc_x0 = -0.322477620583087
    anc_jacobi = 3.183379082910527
    anc_period = 19.503763587070285

    orbit_seed = mu_step_to_orbit(
        anc_mu,
        PC_MU,
        anc_x0,
        anc_jacobi,
        anc_period,
        hc=None,
        sign=-1.0,
        n_steps=40,
    )
    if orbit_seed is None:
        return SweepResult(
            k1=3,
            k2=3,
            stable_found=False,
            method="mu_step",
            note="mu-continuation failed",
        )

    c_lo = max(3.0, orbit_seed.jacobi - 0.3)
    c_hi = min(c_l1 - 0.002, orbit_seed.jacobi + 0.3)
    orbit_stable = c_sweep_find_nu_zero(
        sys_pc,
        orbit_seed.x0,
        orbit_seed.jacobi,
        orbit_seed.period,
        hc=None,
        sign=-1.0,
        c_lo=c_lo,
        c_hi=c_hi,
        n_coarse=60,
    )
    method = f"mu_step_from_mu{anc_mu:.5f}_then_c_sweep"
    if orbit_stable is None:
        return SweepResult(
            k1=3,
            k2=3,
            stable_found=False,
            method=method,
            note="no stable window",
        )
    return _build_result(3, 3, sys_pc, orbit_stable, method)


# ---------------------------------------------------------------------------
# (2,1) sweep — grid search (not in Table-I)
# ---------------------------------------------------------------------------


def _grid_seed_search(
    sys_pc: cr3bp.CR3BPSystem,
    k1_target: int,
    k2_target: int,
    x0_grid: np.ndarray,
    c_grid: np.ndarray,
    hc_list: tuple[int, ...],
    period_guess: float,
    per_call_timeout: int = 4,
) -> cp.SymmetricOrbit | None:
    """Brute-force (x0, C, hc) grid search for a (k1_target, k2_target) orbit.

    Each corrector call is bounded to ``per_call_timeout`` seconds via SIGALRM.
    Returns the first converged orbit with the correct topology, or None.
    """
    for hc_try in hc_list:
        for c_val in c_grid:
            for x0_try in x0_grid:
                x0_v = float(x0_try)
                c_v = float(c_val)
                hc_v = int(hc_try)

                def _fn(
                    _x0: float = x0_v,
                    _c: float = c_v,
                    _hc: int = hc_v,
                ) -> cp.SymmetricOrbit:
                    return cp.correct_symmetric_fixed_jacobi(
                        sys_pc,
                        _x0,
                        _c,
                        period_guess,
                        ydot0_sign=-1.0,
                        half_crossings=_hc,
                        tol=1e-10,
                    )

                try:
                    o = _run_with_timeout(_fn, seconds=per_call_timeout)
                except (ValueError, RuntimeError):
                    continue
                if o is None or not o.converged:
                    continue
                state0 = np.array([o.x0, 0.0, 0.0, 0.0, o.ydot0, 0.0])
                topo = winding_topology(sys_pc.mu, state0, o.period)
                if topo.k1 == k1_target and topo.k2 == k2_target and topo.prograde:
                    return o
    return None


def sweep_21() -> SweepResult:
    """Sweep the (2,1) family at Pluto-Charon mu via grid search.

    (2,1) is not in the Ross-RT Table-I; seed by brute-force (x0, C) grid
    with multiple half-crossings values, classify winding topology, then
    C-sweep any converged (2,1) orbit for the stable window.

    Each corrector call is bounded by a SIGALRM timeout to avoid runaway
    integrations near P2.
    """
    sys_pc = make_pluto_charon_system()
    c_l1 = _c_l1(PC_MU)

    # Coarser grid: 8 x0 x 8 C x 2 hc = 128 trials max; 4s each <= 512s worst case.
    x0_grid = np.linspace(-0.80, -0.35, 8)
    c_grid = np.linspace(3.10, c_l1 - 0.01, 8)
    hc_list = (3, 4)

    seed_orbit = _grid_seed_search(sys_pc, 2, 1, x0_grid, c_grid, hc_list, 14.0)

    if seed_orbit is None:
        return SweepResult(
            k1=2,
            k2=1,
            stable_found=False,
            method="grid_search_8x8x2",
            note="no (2,1) orbit found in grid (x0∈[-0.80,-0.35], C∈[3.10,C_L1-0.01], hc∈(3,4))",
        )

    c_lo = max(3.0, seed_orbit.jacobi - 0.3)
    c_hi = min(c_l1 - 0.002, seed_orbit.jacobi + 0.3)
    orbit_stable = c_sweep_find_nu_zero(
        sys_pc,
        seed_orbit.x0,
        seed_orbit.jacobi,
        seed_orbit.period,
        hc=None,
        sign=-1.0,
        c_lo=c_lo,
        c_hi=c_hi,
        n_coarse=60,
    )
    if orbit_stable is None:
        return SweepResult(
            k1=2,
            k2=1,
            stable_found=False,
            method="grid_then_c_sweep",
            note="no stable window",
        )
    return _build_result(2, 1, sys_pc, orbit_stable, "grid_seed_then_c_sweep")


# ---------------------------------------------------------------------------
# (2,2) sweep — grid search (not in Table-I)
# ---------------------------------------------------------------------------


def sweep_22() -> SweepResult:
    """Sweep the (2,2) family at Pluto-Charon mu via grid search.

    (2,2) is not in the Ross-RT Table-I; grid-search then C-sweep.

    Each corrector call is bounded by a SIGALRM timeout to avoid runaway
    integrations near P2.
    """
    sys_pc = make_pluto_charon_system()
    c_l1 = _c_l1(PC_MU)

    # Coarser grid: 8 x0 x 8 C x 2 hc = 128 trials max; 4s each <= 512s worst case.
    x0_grid = np.linspace(-0.70, -0.20, 8)
    c_grid = np.linspace(3.05, c_l1 - 0.01, 8)
    hc_list = (4, 5)

    seed_orbit = _grid_seed_search(sys_pc, 2, 2, x0_grid, c_grid, hc_list, 18.0)

    if seed_orbit is None:
        return SweepResult(
            k1=2,
            k2=2,
            stable_found=False,
            method="grid_search_8x8x2",
            note="no (2,2) orbit found in grid (x0∈[-0.70,-0.20], C∈[3.05,C_L1-0.01], hc∈(4,5))",
        )

    c_lo = max(3.0, seed_orbit.jacobi - 0.3)
    c_hi = min(c_l1 - 0.002, seed_orbit.jacobi + 0.3)
    orbit_stable = c_sweep_find_nu_zero(
        sys_pc,
        seed_orbit.x0,
        seed_orbit.jacobi,
        seed_orbit.period,
        hc=None,
        sign=-1.0,
        c_lo=c_lo,
        c_hi=c_hi,
        n_coarse=60,
    )
    if orbit_stable is None:
        return SweepResult(
            k1=2,
            k2=2,
            stable_found=False,
            method="grid_then_c_sweep",
            note="no stable window",
        )
    return _build_result(2, 2, sys_pc, orbit_stable, "grid_seed_then_c_sweep")
