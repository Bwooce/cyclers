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


__all__ = [
    "SamplingSectionMap",
    "SectionMap",
    "SectionPoint",
    "SectionReturn",
]
