"""Poincare-section return-map backends for the DA/HOTM enumeration lane (#450).

The lane reduces the planar CR3BP periodic-orbit search to fixed points of the
n-th iterate of a Poincare-section return map. On the section

    Sigma = {y = 0, ydot >= 0, 0 < x < 1 - mu}

the Jacobi integral (one DOF) and the section constraint y=0 (one DOF) leave a
2-D section state ``(x, xdot)``; ``ydot`` is recovered from the Jacobi constant on
the ``ydot >= 0`` branch (paper convention; the Png' family has ydot0 > 0). This is
exactly the reduction the design draft §2 specifies and the same algebra
:func:`cyclerfinder.search.cr3bp_periodic.ydot0_from_jacobi` already implements.

A SINGLE REVOLUTION is the first return to ``Sigma`` *with the same ydot sign*
(``ydot > 0``) -- i.e. the next y=0 crossing whose ydot matches the section
orientation. With this definition a period-n Png' orbit is a fixed point of the
n-th iterate ``P^n`` (verified for P5g': its 5 ydot>0 crossings end exactly at the
published period), matching the paper's "compose one single-rev map n times" /
"revolution count n" parameterization (design draft §2).

Two backends share the :class:`SectionMap` interface (the swappable seam):

* :class:`SamplingSectionMap` -- brute-force float-propagator realization (this
  module), the validation ORACLE. Mirrors ``search/reachable_impulsive.py``: the
  method's geometry without the paper's DA speed.
* :class:`DASectionMap` -- pure-Python truncated Taylor-map (#450 Task 8, added to
  the SAME interface), validated against the sampling oracle. NO MOSEK / DACEyPy.

Pure: math / numpy / scipy + ``cyclerfinder.core.cr3bp`` /
``cyclerfinder.search.cr3bp_periodic``.
"""

from __future__ import annotations

import abc
import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp


@dataclass(frozen=True)
class SectionPoint:
    """A point on the Poincare section ``Sigma``: ``(x, xdot)`` at ``y = 0``."""

    x: float
    xdot: float

    def as_array(self) -> NDArray[np.float64]:
        return np.array([self.x, self.xdot], dtype=np.float64)


@dataclass(frozen=True)
class SectionReturn:
    """The image of a section point under (a power of) the return map.

    ``point`` is the returned section state; ``t`` is the elapsed nondimensional
    time to reach it (the sum over the chained single revolutions for ``compose``).
    """

    point: SectionPoint
    t: float


class SectionMap(abc.ABC):
    """Abstract Poincare-section return map at a fixed Jacobi constant.

    Concrete backends implement :meth:`single_rev`; :meth:`compose` and
    :meth:`residual` are provided generically so the enumerator/driver are
    backend-agnostic.
    """

    def __init__(self, system: cr3bp.CR3BPSystem, c_target: float, *, ydot_sign: float = 1.0):
        self.system = system
        self.c_target = float(c_target)
        self.ydot_sign = float(ydot_sign)

    def lift(self, s: SectionPoint) -> NDArray[np.float64]:
        """Lift a section point ``(x, xdot)`` to the full state on ``Sigma``.

        ``ydot`` is recovered from the Jacobi constant on the ``ydot_sign`` branch
        (default ``+1`` for the ``ydot >= 0`` section). Raises ``ValueError`` if
        the requested Jacobi constant is infeasible at ``(x, xdot)`` (negative
        radicand).
        """
        mu = float(self.system.mu)
        rad = cp._ubar_x_at_axis(s.x, mu) - self.c_target - s.xdot * s.xdot
        if rad < 0.0:
            raise ValueError(
                f"lift: negative Jacobi radicand {rad:.3e} at x={s.x:.6f}, "
                f"xdot={s.xdot:.6f}, C={self.c_target:.6f}"
            )
        ydot = self.ydot_sign * math.sqrt(rad)
        return np.array([s.x, 0.0, 0.0, s.xdot, ydot, 0.0], dtype=np.float64)

    @abc.abstractmethod
    def single_rev(self, s: SectionPoint) -> SectionReturn:
        """First return to ``Sigma`` (same ydot sign) from ``s``.

        Raises ``ValueError`` if ``s`` is infeasible (lift fails) or no such
        return exists within the backend's horizon.
        """

    def compose(self, s: SectionPoint, n: int) -> SectionReturn:
        """The n-th iterate ``P^n(s)`` (n chained single revolutions)."""
        if n < 1:
            raise ValueError(f"compose: n must be >= 1, got {n}")
        cur = s
        total_t = 0.0
        for _ in range(n):
            step = self.single_rev(cur)
            cur = step.point
            total_t += step.t
        return SectionReturn(point=cur, t=total_t)

    def residual(self, s: SectionPoint, n: int) -> float:
        """Section-map fixed-point residual ``||P^n(s) - s||``.

        Returns ``+inf`` for an infeasible point or a missing return (so a grid
        sweep can treat it as "not a candidate" without crashing).
        """
        try:
            img = self.compose(s, n)
        except (ValueError, RuntimeError):
            return float("inf")
        d = img.point.as_array() - s.as_array()
        return float(np.linalg.norm(d))


class SamplingSectionMap(SectionMap):
    """Brute-force float-propagator realization of the section map (the oracle).

    Integrates the lifted IC with ``core.cr3bp`` and detects the first y=0
    crossing whose ``ydot`` matches the section orientation. No differential
    algebra -- the method's geometry without its speed (design draft §0, the
    ``reachable_impulsive.py`` precedent).
    """

    def __init__(
        self,
        system: cr3bp.CR3BPSystem,
        c_target: float,
        *,
        ydot_sign: float = 1.0,
        rtol: float = 1e-12,
        atol: float = 1e-12,
        t_max: float = 8.0,
    ):
        super().__init__(system, c_target, ydot_sign=ydot_sign)
        self.rtol = float(rtol)
        self.atol = float(atol)
        self.t_max = float(t_max)

    def single_rev(self, s: SectionPoint) -> SectionReturn:
        state0 = self.lift(s)  # raises ValueError if infeasible
        mu = float(self.system.mu)

        def _y_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
            return float(y[1])

        _y_event.direction = 0.0  # type: ignore[attr-defined]

        sol = solve_ivp(
            cr3bp.cr3bp_eom,
            (0.0, self.t_max),
            state0,
            args=(mu,),  # type: ignore[call-overload]
            method="DOP853",
            rtol=self.rtol,
            atol=self.atol,
            events=_y_event,
        )
        if not sol.success:
            raise RuntimeError(f"single_rev: propagation failed: {sol.message}")
        # First y=0 crossing (t>0) whose ydot has the section orientation.
        sign = math.copysign(1.0, self.ydot_sign)
        for t_c, yf in zip(sol.t_events[0], sol.y_events[0], strict=True):
            if t_c <= 1e-9 * self.t_max:
                continue
            if math.copysign(1.0, float(yf[4])) == sign:
                return SectionReturn(
                    point=SectionPoint(x=float(yf[0]), xdot=float(yf[3])), t=float(t_c)
                )
        raise RuntimeError(
            f"single_rev: no ydot{'>' if sign > 0 else '<'}0 y=0 return within t_max={self.t_max}"
        )


# ---------------------------------------------------------------------------
# Pure-Python truncated Taylor-map backend (#450 Task 8; USER DECISION 2026-06-25
# option (b): NO MOSEK / DACEyPy / dace / pyaudi). A 2-variable truncated
# polynomial in the section offset (dx, dxdot) about a reference; the single-rev
# map is fitted by finite differences of the float propagator, composed to P^n by
# truncated polynomial composition, and its fixed point found by an iterated
# (re-expanded) trust-region Newton -- then refined (#805) by an EXACT-derivative
# multiple-shooting Newton on the section chain (STM/variational propagation).
# ---------------------------------------------------------------------------

# A monomial dict maps an exponent pair (a, b) -> coefficient for dx^a dxdot^b.
_Poly = dict[tuple[int, int], float]


def _exponents(order: int) -> list[tuple[int, int]]:
    """All exponent pairs (a, b) with a + b <= ``order``."""
    return [(a, total - a) for total in range(order + 1) for a in range(total + 1)]


def _poly_mul(p: _Poly, q: _Poly, order: int) -> _Poly:
    """Multiply two scalar polynomials, truncating at total degree ``order``."""
    out: _Poly = {}
    for (a1, b1), c1 in p.items():
        for (a2, b2), c2 in q.items():
            if a1 + a2 + b1 + b2 <= order:
                key = (a1 + a2, b1 + b2)
                out[key] = out.get(key, 0.0) + c1 * c2
    return out


def _poly_subst(p: _Poly, sx: _Poly, sy: _Poly, order: int) -> _Poly:
    """Substitute (dx, dxdot) -> (sx, sy) into ``p``, truncated at ``order``."""
    max_a = max((a for a, _ in p), default=0)
    max_b = max((b for _, b in p), default=0)
    pow_x: dict[int, _Poly] = {0: {(0, 0): 1.0}}
    for i in range(1, max_a + 1):
        pow_x[i] = _poly_mul(pow_x[i - 1], sx, order)
    pow_y: dict[int, _Poly] = {0: {(0, 0): 1.0}}
    for i in range(1, max_b + 1):
        pow_y[i] = _poly_mul(pow_y[i - 1], sy, order)
    out: _Poly = {}
    for (a, b), c in p.items():
        term = _poly_mul(pow_x[a], pow_y[b], order)
        for e, v in term.items():
            out[e] = out.get(e, 0.0) + c * v
    return out


def _poly_eval(p: _Poly, dx: float, dxd: float) -> float:
    return float(sum(c * dx**a * dxd**b for (a, b), c in p.items()))


@dataclass(frozen=True)
class TaylorMap2:
    """A 2-D truncated Taylor map: section offset (dx, dxdot) -> output offset.

    ``px`` / ``pxd`` are the polynomials for the x and xdot output offsets
    (relative to the reference about which the map was expanded), of total degree
    <= ``order``.
    """

    px: _Poly
    pxd: _Poly
    order: int

    def evaluate(self, dx: float, dxd: float) -> tuple[float, float]:
        """Output offset (dx_out, dxdot_out) at input offset ``(dx, dxd)``."""
        return _poly_eval(self.px, dx, dxd), _poly_eval(self.pxd, dx, dxd)

    def compose_self(self, n: int) -> TaylorMap2:
        """The n-fold composition ``map o map o ... o map`` (n times)."""
        if n < 1:
            raise ValueError(f"compose_self: n must be >= 1, got {n}")
        cx, cy = self.px, self.pxd
        for _ in range(n - 1):
            cx, cy = (
                _poly_subst(cx, self.px, self.pxd, self.order),
                _poly_subst(cy, self.px, self.pxd, self.order),
            )
        return TaylorMap2(px=cx, pxd=cy, order=self.order)


class DASectionMap(SamplingSectionMap):
    """Pure-Python truncated Taylor-map section-map backend (the deliverable).

    Inherits the float ``single_rev`` from :class:`SamplingSectionMap` (so
    single_rev / compose are bit-for-bit the same geometry -- the swappable-seam
    parity the design requires). The Taylor layer adds:

    * :meth:`taylor_single_rev` -- fit the single-rev map to a truncated
      polynomial about a reference by finite differences (deliberately so --
      see its docstring for the #805 finding);
    * :meth:`taylor_fixed_point` -- compose to ``P^n``, find its fixed point by
      an iterated (re-expanded) trust-region Newton, then refine through the
      exact-derivative endgame;
    * :meth:`single_rev_stm` -- single revolution plus the EXACT 2x2 section-map
      Jacobian assembled from the state+STM variational propagation (#805);
    * :meth:`section_chain_newton` -- exact-derivative multiple-shooting Newton
      on the n-rev section chain (the #805 endgame).

    No differential-algebra library, no MOSEK. The FD-coefficient accuracy floors
    the POLYNOMIAL stage's fixed-point distance for strongly-unstable multi-rev
    orbits, and that floor is MACHINE-DEPENDENT (~3e-5 Linux/OpenBLAS build
    machine, ~2.8e-4 macOS/Accelerate; #450 decision note + #804 note):
    compose_self evaluates the [-h, h]-fitted single-rev polynomial at the
    orbit's other chain points far outside the fit domain, so tiny
    BLAS/integrator coefficient differences shift where the re-expansion
    iteration lands. Since #805 the chain-Newton endgame converges from that
    floor to the TRUE float-map fixed point (~1e-11 chain residual; P5g' lands
    ~4e-9 from the published IC on the macOS machine), eliminating the
    machine-dependence of the landing whenever the endgame accepts. The
    corrector (Task 5) still performs the certified closure.
    """

    def single_rev_stm(self, s: SectionPoint) -> tuple[SectionReturn, NDArray[np.float64]]:
        """First return to ``Sigma`` from ``s`` plus the EXACT section-map Jacobian.

        Integrates the state + 6x6 STM variational equations
        (:func:`cyclerfinder.core.cr3bp.cr3bp_stm_eom`) to the same first
        same-ydot-sign y=0 crossing as :meth:`single_rev` and assembles the 2x2
        section-map Jacobian ``d(x', xdot') / d(x, xdot)`` from three exact
        pieces (#805 -- no finite differencing anywhere):

        * the LIFT derivative ``d(state0)/d(x, xdot)`` -- ``ydot0`` is a function
          of ``(x, xdot)`` through the Jacobi constant, so
          ``d ydot0/dx = -(dUbar/dx)/ydot0`` and ``d ydot0/dxdot = -xdot/ydot0``;
        * the flow STM ``Phi(t_c)`` from the variational propagation;
        * the SECTION projection -- the crossing time varies with the initial
          condition to hold ``y(t_c) = 0``, so
          ``delta t = -(Phi L ds)_y / ydot_f`` and the full-state variation is
          ``Phi L ds + f(X_f) delta t`` (rows x and xdot are kept).

        Raises ``ValueError`` if ``s`` is infeasible or section-tangent
        (``ydot0 ~ 0``: the sqrt-branch lift derivative is singular), and
        ``RuntimeError`` if the propagation fails or no same-sign return exists
        within ``t_max`` (same contract as :meth:`single_rev`).
        """
        state0 = self.lift(s)  # raises ValueError if infeasible
        mu = float(self.system.mu)
        ydot0 = float(state0[4])
        if abs(ydot0) < 1e-12:
            raise ValueError(
                f"single_rev_stm: section tangency (ydot0={ydot0:.3e}) at "
                f"x={s.x:.6f}, xdot={s.xdot:.6f} -- lift derivative singular"
            )
        aug0 = np.concatenate([state0, np.eye(6).reshape(36)])

        def _y_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
            return float(y[1])

        _y_event.direction = 0.0  # type: ignore[attr-defined]

        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom,
            (0.0, self.t_max),
            aug0,
            args=(mu,),  # type: ignore[call-overload]
            method="DOP853",
            rtol=self.rtol,
            atol=self.atol,
            events=_y_event,
        )
        if not sol.success:
            raise RuntimeError(f"single_rev_stm: propagation failed: {sol.message}")
        sign = math.copysign(1.0, self.ydot_sign)
        for t_c, yf in zip(sol.t_events[0], sol.y_events[0], strict=True):
            if t_c <= 1e-9 * self.t_max:
                continue
            if math.copysign(1.0, float(yf[4])) == sign:
                break
        else:
            raise RuntimeError(
                f"single_rev_stm: no ydot{'>' if sign > 0 else '<'}0 y=0 return "
                f"within t_max={self.t_max}"
            )
        state_f = np.asarray(yf[:6], dtype=np.float64)
        phi = np.asarray(yf[6:], dtype=np.float64).reshape(6, 6)
        # Lift derivative L = d(state0)/d(x, xdot) (6x2). _ubar_x_at_axis is
        # -2*Ubar on the axis, so d(rad)/dx = -2*dUbar/dx and
        # d ydot0/dx = (d rad/dx)/(2 ydot0) = -(dUbar/dx)/ydot0.
        lmat = np.zeros((6, 2), dtype=np.float64)
        lmat[0, 0] = 1.0
        lmat[3, 1] = 1.0
        lmat[4, 0] = -cp._ubar_grad_x_at_axis(s.x, mu) / ydot0
        lmat[4, 1] = -s.xdot / ydot0
        m = phi @ lmat  # free (fixed-time) final-state variation, 6x2
        f = cr3bp.cr3bp_eom(float(t_c), state_f, mu)
        fy = float(f[1])  # = ydot_f, nonzero by the same-sign crossing selection
        jac = np.empty((2, 2), dtype=np.float64)
        jac[0, :] = m[0, :] - (float(f[0]) / fy) * m[1, :]
        jac[1, :] = m[3, :] - (float(f[3]) / fy) * m[1, :]
        ret = SectionReturn(
            point=SectionPoint(x=float(state_f[0]), xdot=float(state_f[3])), t=float(t_c)
        )
        return ret, jac

    def taylor_single_rev(
        self, s_ref: SectionPoint, *, order: int, h: float, samples: int
    ) -> TaylorMap2:
        """Fit the single-rev map about ``s_ref`` to a degree-``order`` polynomial.

        Samples a ``samples x samples`` tensor grid of offsets in ``[-h, h]^2``,
        evaluates the float single-rev image at each feasible sample, and
        least-squares fits the output OFFSET (image minus reference) to the
        monomial basis. Raises ``ValueError`` if too few samples survive.

        DELIBERATELY finite-difference, including the affine block (#805
        finding, docs/notes/2026-08-09-805): substituting the EXACT constant /
        first-order coefficients from :meth:`single_rev_stm` makes the iterated
        fixed-point descent WORSE (P5g': 2.8e-4 -> ~1e-3, wandering then
        divergence), because the least-squares affine block is a [-h, h]^2
        domain-AVERAGED linearization (it differs from the true point Jacobian
        by 6-25% under this map's violent curvature) and the out-of-domain
        compose_self extrapolation that drives the descent depends on exactly
        that averaging. The exact-derivative machinery enters through the
        :meth:`section_chain_newton` endgame instead.
        """
        exps = _exponents(order)
        offs = np.linspace(-h, h, samples)
        rows: list[tuple[float, float]] = []
        outs: list[tuple[float, float]] = []
        for a in offs:
            for c in offs:
                try:
                    img = self.single_rev(
                        SectionPoint(x=s_ref.x + float(a), xdot=s_ref.xdot + float(c))
                    )
                except (ValueError, RuntimeError):
                    continue
                rows.append((float(a), float(c)))
                outs.append((img.point.x - s_ref.x, img.point.xdot - s_ref.xdot))
        if len(rows) < len(exps) + 2:
            raise ValueError(
                f"taylor_single_rev: only {len(rows)} feasible samples for "
                f"{len(exps)} coefficients at s_ref={s_ref}, h={h}"
            )
        offs_arr = np.array(rows, dtype=np.float64)
        out_arr = np.array(outs, dtype=np.float64)
        design = np.column_stack([(offs_arr[:, 0] ** a) * (offs_arr[:, 1] ** b) for a, b in exps])
        cx, *_ = np.linalg.lstsq(design, out_arr[:, 0], rcond=None)
        cxd, *_ = np.linalg.lstsq(design, out_arr[:, 1], rcond=None)
        px = {e: float(cx[i]) for i, e in enumerate(exps)}
        pxd = {e: float(cxd[i]) for i, e in enumerate(exps)}
        return TaylorMap2(px=px, pxd=pxd, order=order)

    def section_chain_newton(
        self,
        s: SectionPoint,
        n: int,
        *,
        max_iter: int = 20,
        accept_tol: float = 1e-9,
    ) -> SectionPoint | None:
        """Exact-derivative multiple-shooting Newton on the n-rev section chain.

        Solves ``P(s_i) = s_{i+1 mod n}`` for all n chain nodes simultaneously,
        with each 2x2 block Jacobian EXACT from :meth:`single_rev_stm` (#805).
        Multiple shooting distributes the orbit's instability across the
        per-rev blocks (|J_rev| ~ 5e2 for P5g') instead of compounding it into
        the ``P^n`` composition (|J^n| ~ 3.6e3, whose direct Newton diverges --
        measured, #805 note), so from a Taylor-stage landing whose float chain
        exists (residual finite, ~3e-5..3e-4 from the fixed point) it converges
        to the TRUE float-map fixed point at ~1e-11 chain residual. Globalized
        by backtracking line search on the chain-residual norm.

        Returns the refined first node, or ``None`` if the initial chain does
        not exist (off-family point: some rev has no same-sign return -- the
        cheap fail-fast for enumerator sweeps), the Newton system is singular,
        or the best chain residual achieved is above ``accept_tol``. Callers
        fall back to the unrefined point on ``None``, so this endgame can only
        improve a landing, never worsen it.
        """
        if n < 1:
            raise ValueError(f"section_chain_newton: n must be >= 1, got {n}")
        # Initial chain from s: fails fast for off-family points.
        nodes = [s]
        try:
            for _ in range(n - 1):
                nodes.append(self.single_rev(nodes[-1]).point)
        except (ValueError, RuntimeError):
            return None

        def _chain_residual(pts: list[SectionPoint]) -> NDArray[np.float64] | None:
            r = np.zeros(2 * n, dtype=np.float64)
            for i in range(n):
                try:
                    img = self.single_rev(pts[i]).point
                except (ValueError, RuntimeError):
                    return None
                j = (i + 1) % n
                r[2 * i] = img.x - pts[j].x
                r[2 * i + 1] = img.xdot - pts[j].xdot
            return r

        best: tuple[float, list[SectionPoint]] | None = None
        for _ in range(max_iter):
            rvec = np.zeros(2 * n, dtype=np.float64)
            jmat = np.zeros((2 * n, 2 * n), dtype=np.float64)
            try:
                for i in range(n):
                    ret, jac = self.single_rev_stm(nodes[i])
                    j = (i + 1) % n
                    rvec[2 * i] = ret.point.x - nodes[j].x
                    rvec[2 * i + 1] = ret.point.xdot - nodes[j].xdot
                    jmat[2 * i : 2 * i + 2, 2 * i : 2 * i + 2] = jac
                    jmat[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] -= np.eye(2)
            except (ValueError, RuntimeError):
                break
            rnorm = float(np.linalg.norm(rvec))
            if best is None or rnorm < best[0]:
                best = (rnorm, nodes)
            if rnorm < 1e-12:
                break
            try:
                step = np.linalg.solve(jmat, -rvec)
            except np.linalg.LinAlgError:
                break
            # Backtracking on the chain-residual norm; a trial whose chain walls
            # off is rejected the same as a non-decreasing one.
            lam = 1.0
            accepted = False
            for _ in range(8):
                trial = [
                    SectionPoint(
                        x=nodes[i].x + lam * float(step[2 * i]),
                        xdot=nodes[i].xdot + lam * float(step[2 * i + 1]),
                    )
                    for i in range(n)
                ]
                r_trial = _chain_residual(trial)
                if r_trial is not None and float(np.linalg.norm(r_trial)) < rnorm:
                    accepted = True
                    break
                lam *= 0.5
            if not accepted:
                break
            nodes = trial
        if best is not None and best[0] <= accept_tol:
            return best[1][0]
        return None

    @staticmethod
    def _poly_fixed_point(pn: TaylorMap2, trust: float) -> NDArray[np.float64] | None:
        """Solve ``P^n(d) = d`` for the offset ``d`` by trust-region Newton.

        ``P^n`` is the offset-out map, so the fixed-point condition is
        ``px(d) - dx = 0, pxd(d) - dxd = 0``. Returns ``None`` if the iterate
        leaves the trust region (a spurious out-of-domain polynomial root).
        """
        d = np.zeros(2, dtype=np.float64)

        def _resid(v: NDArray[np.float64]) -> NDArray[np.float64]:
            ex, ey = pn.evaluate(float(v[0]), float(v[1]))
            return np.array([ex - v[0], ey - v[1]], dtype=np.float64)

        for _ in range(80):
            f = _resid(d)
            if float(np.linalg.norm(f)) < 1e-15:
                break
            eps = 1e-8
            jac = np.zeros((2, 2))
            for k in range(2):
                dp = d.copy()
                dp[k] += eps
                jac[:, k] = (_resid(dp) - f) / eps
            try:
                step = np.linalg.solve(jac, -f)
            except np.linalg.LinAlgError:
                break
            ns = float(np.linalg.norm(step))
            if ns > trust:
                step = step * (trust / ns)
            d = d + step
            if float(np.linalg.norm(d)) > 4.0 * trust:
                return None
        return d

    def taylor_fixed_point(
        self,
        s_ref: SectionPoint,
        *,
        n: int,
        order: int,
        h: float,
        samples: int,
        max_iter: int = 30,
        refine: bool = True,
    ) -> SectionPoint:
        """Iterated Taylor-map fixed point of ``P^n`` from a coarse reference.

        Each pass re-expands the single-rev map about the current iterate, composes
        to ``P^n``, solves the polynomial fixed point in a trust region, and moves
        the reference there. If a pass fails (too few feasible samples or an
        out-of-domain root) the iteration stops at the best iterate so far.

        With ``refine=True`` (default, #805) the FD-floor landing is then handed
        to the exact-derivative :meth:`section_chain_newton` endgame: when the
        landing's float chain exists (i.e. the descent genuinely reached a
        family), the endgame converges to the TRUE float-map fixed point
        (~1e-11 chain residual), replacing the machine-dependent FD floor
        (~3e-5 Linux .. ~3e-4 macOS, #804 note) with an integrator-accuracy
        landing. If the endgame declines (off-family landing, no chain), the
        FD-floor landing is returned unchanged -- the corrector still finishes
        from there (Task 5).
        """
        cur = s_ref
        prev_step = float("inf")
        stagnant = 0
        for _ in range(max_iter):
            try:
                tmap = self.taylor_single_rev(cur, order=order, h=h, samples=samples)
            except ValueError:
                break
            pn = tmap.compose_self(n)
            d = self._poly_fixed_point(pn, trust=h * 3.0)
            if d is None:
                break
            cur = SectionPoint(x=cur.x + float(d[0]), xdot=cur.xdot + float(d[1]))
            step = math.hypot(float(d[0]), float(d[1]))
            if step < 1e-11:
                break
            # Stagnation break: if the step stops shrinking, further passes only
            # burn propagations at the FD floor. Bail after a couple of stalls.
            if step >= 0.9 * prev_step:
                stagnant += 1
                if stagnant >= 2:
                    break
            else:
                stagnant = 0
            prev_step = step
        if refine:
            refined = self.section_chain_newton(cur, n)
            if refined is not None:
                return refined
        return cur


__all__ = [
    "DASectionMap",
    "SamplingSectionMap",
    "SectionMap",
    "SectionPoint",
    "SectionReturn",
    "TaylorMap2",
]
