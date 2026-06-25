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
# (re-expanded) trust-region Newton.
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
      polynomial about a reference by finite differences;
    * :meth:`taylor_fixed_point` -- compose to ``P^n`` and find its fixed point by
      an iterated (re-expanded) trust-region Newton, returning the section point.

    No differential-algebra library, no MOSEK. The FD-coefficient accuracy floors
    the achievable fixed-point distance for strongly-unstable multi-rev orbits
    (~3e-5 for P5g'); the corrector finishes the closure (Task 5).
    """

    def taylor_single_rev(
        self, s_ref: SectionPoint, *, order: int, h: float, samples: int
    ) -> TaylorMap2:
        """Fit the single-rev map about ``s_ref`` to a degree-``order`` polynomial.

        Samples a ``samples x samples`` tensor grid of offsets in ``[-h, h]^2``,
        evaluates the float single-rev image at each feasible sample, and
        least-squares fits the output OFFSET (image minus reference) to the
        monomial basis. Raises ``ValueError`` if too few samples survive.
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
    ) -> SectionPoint:
        """Iterated Taylor-map fixed point of ``P^n`` from a coarse reference.

        Each pass re-expands the single-rev map about the current iterate, composes
        to ``P^n``, solves the polynomial fixed point in a trust region, and moves
        the reference there. Returns the converged section point. If a pass fails
        (too few feasible samples or an out-of-domain root) the best iterate so far
        is returned -- the corrector still finishes from there (Task 5).
        """
        cur = s_ref
        prev_step = float("inf")
        stagnant = 0
        for _ in range(max_iter):
            try:
                tmap = self.taylor_single_rev(cur, order=order, h=h, samples=samples)
            except ValueError:
                return cur
            pn = tmap.compose_self(n)
            d = self._poly_fixed_point(pn, trust=h * 3.0)
            if d is None:
                return cur
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
        return cur


__all__ = [
    "DASectionMap",
    "SamplingSectionMap",
    "SectionMap",
    "SectionPoint",
    "SectionReturn",
    "TaylorMap2",
]
