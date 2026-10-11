"""#1024: a barycentric circular-model option for binaries in the two-working-body generator.

The default ideal model (:class:`two_working_body.CircularSystem`) puts every body on a circle
about the origin whose radius follows from its period by Kepler III with the system GM ``mu``. For
a binary with a heavy secondary (Pluto-Charon, mass fraction 0.109) that places the secondary at
the radius of the RELATIVE orbit (19,596 km for Charon) at the speed sqrt(mu / a), about 12 %
above its barycentric radius and speed.

This module adds the option (docs/notes/2026-10-11-1024-barycentric-binary-model.md sec. 2):

- origin at the barycentre, spacecraft legs two-body about it with the system GM (unchanged ``mu``);
- the secondary on a circle of radius ``(1 - mu_S) a_rel`` at its period, moving at the KINEMATIC
  speed ``2 pi r / P`` (it is not a test particle of the barycentre monopole);
- the primary on the opposite circle of radius ``mu_S a_rel``, used only for distances (the surface
  screen), not as a gravitating point of the legs;
- every other body unchanged.

The option is opt-in: :func:`barycentric_binary` builds a :class:`BarycentricCircularSystem` from an
existing system; nothing else changes. Code that assumes a body is Keplerian about ``mu`` checks
:meth:`BarycentricCircularSystem.is_keplerian`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np
from scipy.optimize import minimize_scalar

from cyclerfinder.search.two_working_body import (
    CircularSystem,
    Cycle,
    Flyby,
    LambertLeg,
    Vec,
    _blocks,
    eval_lambert_legs,
    fixed_duration_s,
    kepler_step,
)


@dataclass
class BarycentricCircularSystem(CircularSystem):
    """A :class:`CircularSystem` with kinematic-speed bodies and a moving primary.

    ``kinematic`` names bodies whose speed is ``2 pi a / P`` (their circle is not a Keplerian orbit
    about ``mu``). ``primary_orbit`` is ``(radius_km, period_s, theta0_rad)`` of the primary's own
    circle about the barycentre, or ``None``.
    """

    kinematic: frozenset[str] = frozenset()
    primary_orbit: tuple[float, float, float] | None = None

    def state(self, code: str, t_s: float) -> tuple[Vec, Vec]:
        if code not in self.kinematic:
            return super().state(code, t_s)
        a, per, th0 = self.bodies[code]
        th = th0 + 2.0 * math.pi * t_s / per
        v = 2.0 * math.pi * a / per
        return (
            np.array([a * math.cos(th), a * math.sin(th), 0.0]),
            np.array([-v * math.sin(th), v * math.cos(th), 0.0]),
        )

    def is_keplerian(self, code: str) -> bool:
        """True when the body's circle is a Keplerian orbit about ``mu``."""
        return code not in self.kinematic

    def primary_state(self, t_s: float) -> tuple[Vec, Vec]:
        """Position and velocity of the moving primary (km, km/s)."""
        if self.primary_orbit is None:
            return np.zeros(3), np.zeros(3)
        a, per, th0 = self.primary_orbit
        th = th0 + 2.0 * math.pi * t_s / per
        v = 2.0 * math.pi * a / per
        return (
            np.array([a * math.cos(th), a * math.sin(th), 0.0]),
            np.array([-v * math.sin(th), v * math.cos(th), 0.0]),
        )


def barycentric_binary(
    system: CircularSystem, secondary: str, mass_fraction: float
) -> BarycentricCircularSystem:
    """The barycentric option of ``system`` for the binary (primary, ``secondary``).

    ``mass_fraction`` is ``GM_S / (GM_P + GM_S)``; the secondary's existing radius is the relative
    orbit's ``a_rel``. ``mass_fraction = 0`` reproduces the default states (to round-off: the
    secondary's speed is then ``2 pi a / P``, equal to ``sqrt(mu / a)`` for a Keplerian circle).
    """
    if not 0.0 <= mass_fraction < 1.0:
        raise ValueError(f"mass_fraction must be in [0, 1): {mass_fraction}")
    a_rel, per, th0 = system.bodies[secondary]
    bodies = dict(system.bodies)
    bodies[secondary] = ((1.0 - mass_fraction) * a_rel, per, th0)
    return BarycentricCircularSystem(
        mu=system.mu,
        bodies=bodies,
        massless=system.massless,
        flyby_overrides=dict(system.flyby_overrides),
        kinematic=frozenset({secondary}),
        primary_orbit=(mass_fraction * a_rel, per, th0 + math.pi),
    )


@dataclass
class _Arc:
    t0: float
    r0: Vec
    v0: Vec
    dt: float


@dataclass
class PrimaryApproach:
    """Closest approach to the moving primary over a cycle."""

    min_km: float
    t_days: float
    arcs: int = field(default=0)


def _cycle_arcs(
    system: CircularSystem, cycle: Cycle, x: np.ndarray, flybys: list[Flyby]
) -> list[_Arc]:
    legs = eval_lambert_legs(system, cycle, x)
    if legs is None:
        raise ValueError("Lambert legs do not evaluate at these dates")
    arcs = []
    for ev, li in zip(legs, cycle.lambert_index, strict=True):
        leg = cycle.legs[li]
        assert isinstance(leg, LambertLeg)
        r0, w0 = system.state(leg.frm, ev.t_dep)
        arcs.append(_Arc(ev.t_dep, r0, w0 + ev.vinf_dep, ev.t_arr - ev.t_dep))
    for blk in _blocks(system, cycle, legs):
        for leg_f, t0 in blk.fixed:
            fb = next(f for f in flybys if f.body == blk.body and abs(f.t_s - t0) < 1.0)
            r0, w0 = system.state(blk.body, t0)
            arcs.append(_Arc(t0, r0, w0 + fb.vinf_out, fixed_duration_s(system, leg_f)))
    return arcs


def primary_min_distance_km(
    system: BarycentricCircularSystem,
    cycle: Cycle,
    x: np.ndarray,
    flybys: list[Flyby],
    n_samples: int = 400,
) -> PrimaryApproach:
    """Closest approach of the spacecraft to the MOVING primary over every leg of the cycle
    (Lambert and fixed legs), sampled and then refined by bounded minimisation."""
    best = (math.inf, 0.0)
    arcs = _cycle_arcs(system, cycle, x, flybys)
    for arc in arcs:

        def dist(t: float, arc: _Arc = arc) -> float:
            r = kepler_step(arc.r0, arc.v0, float(t), system.mu)[0] if t > 0.0 else arc.r0
            return float(np.linalg.norm(r - system.primary_state(arc.t0 + t)[0]))

        ts = np.linspace(0.0, arc.dt, n_samples)
        ds = np.array([dist(float(t)) for t in ts])
        i = int(np.argmin(ds))
        lo, hi = float(ts[max(i - 1, 0)]), float(ts[min(i + 1, n_samples - 1)])
        res = minimize_scalar(dist, bounds=(lo, hi), method="bounded", options={"xatol": 1.0})
        d_min, t_min = (float(res.fun), float(res.x)) if res.fun < ds[i] else (ds[i], float(ts[i]))
        if d_min < best[0]:
            best = (d_min, (arc.t0 + t_min) / 86400.0)
    return PrimaryApproach(best[0], best[1], len(arcs))


def with_mass_fraction(
    system: BarycentricCircularSystem, secondary: str, a_rel: float, mass_fraction: float
) -> BarycentricCircularSystem:
    """The same option at another mass fraction (a continuation step), from the relative radius."""
    _, per, th0 = system.bodies[secondary]
    bodies = dict(system.bodies)
    bodies[secondary] = ((1.0 - mass_fraction) * a_rel, per, th0)
    return replace(
        system,
        bodies=bodies,
        primary_orbit=(mass_fraction * a_rel, per, th0 + math.pi),
        _fb={},
    )
