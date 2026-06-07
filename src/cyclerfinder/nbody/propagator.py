"""The Propagator interface + the REBOUND backend (design §1, §3).

A single ``Propagator`` protocol over the **TDB-J2000-seconds axis +
heliocentric J2000-ecliptic frame** that ``core/ephemeris.py`` defines, so the
integrator backend (REBOUND now, a deferred Tudat slot per design Q1) is a swap,
not a rewrite — mirroring how ``core/ephemeris.py`` hides backends behind one
``state()``.

The baseline backend is :class:`RestrictedNBody` on REBOUND / IAS15. The
spacecraft is a massless test particle in a heliocentric simulation; the Sun is
the central mass and the (optional) planet perturbers act through the rails
third-body force callback (:mod:`cyclerfinder.nbody.forces`), read on rails from
the shared DE440 BSP. Units are fixed to km, km/s, s with ``mu`` from
``core/constants`` (REBOUND is unit-agnostic; we pin G via the Sun mass so that
``G * M_sun = MU_SUN_KM3_S2``).
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.nbody.forces import ingest_planet_state

if TYPE_CHECKING:
    from cyclerfinder.core.ephemeris import Ephemeris

Vec3 = NDArray[np.float64]

# REBOUND is unit-agnostic; we work in km, km/s, s. Pin G so that the central
# Sun particle of mass 1.0 reproduces MU_SUN_KM3_S2 exactly (G * 1.0 == mu_sun).
_G_KM3_KG_S2 = MU_SUN_KM3_S2  # with M_sun := 1.0 mass unit


@dataclass(frozen=True)
class NBodyArc:
    """Frozen result of one propagation (design §5 result-dataclass discipline)."""

    r_km: Vec3
    v_km_s: Vec3
    t1_sec: float
    energy_rel_drift: float
    anchor_err_km: float
    integrator_accuracy: float
    bodies: tuple[str, ...]
    converged: bool


class Propagator(Protocol):
    """One method: propagate a state over the TDB-J2000 axis (design §1)."""

    def propagate(
        self,
        r0_km: Vec3,
        v0_km_s: Vec3,
        t0_sec: float,
        t1_sec: float,
        *,
        bodies: Sequence[str],
        accuracy: float,
    ) -> NBodyArc: ...


class RestrictedNBody:
    """Restricted-n-body propagator (REBOUND / IAS15 baseline; design §1, Q1)."""

    _SUPPORTED = ("rebound",)

    def __init__(self, backend: str = "rebound") -> None:
        if backend not in self._SUPPORTED:
            raise ValueError(
                f"unknown propagator backend {backend!r}; supported: "
                f"{self._SUPPORTED}. Tudat is a deferred slot (design Q1): wired "
                "only on demonstrated finite-diff STM need, with its own extra."
            )
        self.backend = backend

    def propagate(
        self,
        r0_km: Vec3,
        v0_km_s: Vec3,
        t0_sec: float,
        t1_sec: float,
        *,
        bodies: Sequence[str] = (),
        accuracy: float = 1e-10,
        ephem: Ephemeris | None = None,
    ) -> NBodyArc:
        """Propagate ``(r0, v0)`` from ``t0_sec`` to ``t1_sec`` in restricted n-body.

        ``bodies`` are perturber body codes (read on rails via ``ephem``); empty
        means Sun-only two-body (GOLDEN GATE 1). ``accuracy`` is the IAS15
        ``epsilon`` tolerance. Returns a frozen :class:`NBodyArc` with the final
        state and the relative-energy-drift diagnostic.
        """
        import rebound

        bodies = tuple(bodies)
        if bodies and ephem is None:
            raise ValueError("perturber bodies requested but no ephem supplied")

        r0 = np.asarray(r0_km, dtype=np.float64)
        v0 = np.asarray(v0_km_s, dtype=np.float64)

        sim = rebound.Simulation()
        sim.G = _G_KM3_KG_S2
        sim.integrator = "ias15"
        # REBOUND 5.0: the IAS15 accuracy parameter lives on the integrator
        # configuration object (sim.integrator.epsilon), not the removed
        # sim.ri_ias15.epsilon accessor of the 3.x/4.x line.
        sim.integrator.epsilon = accuracy
        # Central Sun (mass 1.0 => G*M = mu_sun) at the heliocentric origin.
        sim.add(m=1.0, x=0.0, y=0.0, z=0.0, vx=0.0, vy=0.0, vz=0.0)
        # Massless spacecraft test particle.
        sim.add(
            m=0.0,
            x=float(r0[0]),
            y=float(r0[1]),
            z=float(r0[2]),
            vx=float(v0[0]),
            vy=float(v0[1]),
            vz=float(v0[2]),
        )
        sim.t = float(t0_sec)

        # Rails third-body perturbation as an additional (velocity-independent)
        # force. The Sun central term is handled by REBOUND's own gravity; this
        # callback adds only the planet perturbations (direct + indirect).
        if bodies:
            _install_rails_forces(sim, bodies, ephem)

        e0 = float(sim.energy())
        sim.integrate(float(t1_sec))
        e1 = float(sim.energy())
        energy_rel_drift = (e1 - e0) / abs(e0) if e0 != 0.0 else 0.0

        sc = sim.particles[1]
        r_km = np.array([sc.x, sc.y, sc.z], dtype=np.float64)
        v_km_s = np.array([sc.vx, sc.vy, sc.vz], dtype=np.float64)

        return NBodyArc(
            r_km=r_km,
            v_km_s=v_km_s,
            t1_sec=float(t1_sec),
            energy_rel_drift=float(energy_rel_drift),
            anchor_err_km=0.0,
            integrator_accuracy=float(accuracy),
            bodies=bodies,
            converged=True,
        )


def _install_rails_forces(sim: object, bodies: tuple[str, ...], ephem: Ephemeris | None) -> None:
    """Attach the rails third-body perturbation to a REBOUND simulation.

    The callback reads each perturber's heliocentric state at the integrator's
    *current* sub-step time (``sim.t``) — this is the per-force-eval ephemeris
    read that is the harness's cost (design Risk 3). Sun-frame indirect term
    included so the heliocentric (non-inertial) frame is consistent.
    """
    from cyclerfinder.core.constants import PLANETS

    assert ephem is not None

    def additional_forces(reb_sim_pointer: object) -> None:
        reb_sim = reb_sim_pointer.contents  # type: ignore[attr-defined]
        t_sec = float(reb_sim.t)
        sc = reb_sim.particles[1]
        r = np.array([sc.x, sc.y, sc.z], dtype=np.float64)
        ax = ay = az = 0.0
        for body in bodies:
            mu_p = PLANETS[body].mu_km3_s2
            r_p, _ = ingest_planet_state(body, t_sec, ephem)
            d = r_p - r
            d3 = float(np.linalg.norm(d)) ** 3
            rp3 = float(np.linalg.norm(r_p)) ** 3
            acc = mu_p * (d / d3 - r_p / rp3)
            ax += float(acc[0])
            ay += float(acc[1])
            az += float(acc[2])
        sc.ax += ax
        sc.ay += ay
        sc.az += az

    # The setter wraps `func` in a CFUNCTYPE and stores it in the simulation's
    # own `_afp` slot, so the callback is kept alive for us (REBOUND 5.0); we do
    # NOT attach our own attribute — Simulation.__setattr__ rejects unknown names.
    sim.additional_forces = additional_forces  # type: ignore[attr-defined]
    sim.force_is_velocity_dependent = 0  # type: ignore[attr-defined]


__all__ = ["NBodyArc", "Propagator", "RestrictedNBody"]
