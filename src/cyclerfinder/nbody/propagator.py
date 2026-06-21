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
from typing import TYPE_CHECKING, Any, Protocol

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import MU_SUN_KM3_S2
from cyclerfinder.nbody.forces import RailsEphemerisCache

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
    """Reserved closure-error slot — ALWAYS 0.0 from :class:`RestrictedNBody`.

    The propagator seeds the spacecraft at the *exact* ``(r0, v0)`` it is handed
    and performs no ephemeris-anchoring step, so it has no anchor/closure error of
    its own to report; this field is reserved for a future backend that does anchor
    against an ephemeris. **Callers that need a closure error must compute it
    themselves** from the returned terminal state — :mod:`cyclerfinder.nbody.rung`
    does exactly this (``terminal_closure_km`` = ``|r_wrap - r_seed|``). Do NOT
    read this field as a meaningful closure metric for the REBOUND backend.
    """
    integrator_accuracy: float
    bodies: tuple[str, ...]
    converged: bool
    stm: NDArray[np.float64] | None = None
    """6x6 state-transition matrix ``Phi[:, k] = d(final_state) / d(init_state_k)``.

    ``None`` unless ``propagate(with_stm=True)`` (default), so existing callers are
    unchanged. The state ordering is ``[x, y, z, vx, vy, vz]`` (km, km/s) — the same
    convention as :class:`cyclerfinder.core.cr3bp.CR3BPArc.stm`. Co-integrated via
    REBOUND's native first-order variational particles in the SAME adaptive
    integration as the state path (the variational ODE is linear given that path),
    so it costs one extra integration, not six finite-difference re-propagations.
    If the integration budgets out / diverges (``converged=False``), ``stm`` carries
    the best-effort variational state at the last finite step rather than crashing.
    """


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
        with_stm: bool = ...,
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
        cache: RailsEphemerisCache | None = None,
        max_steps: int = 2_000_000,
        max_wall_sec: float = 90.0,
        with_stm: bool = False,
    ) -> NBodyArc:
        """Propagate ``(r0, v0)`` from ``t0_sec`` to ``t1_sec`` in restricted n-body.

        ``bodies`` are perturber body codes (read on rails via ``ephem``); empty
        means Sun-only two-body (GOLDEN GATE 1). ``accuracy`` is the IAS15
        ``epsilon`` tolerance. Returns a frozen :class:`NBodyArc` with the final
        state and the relative-energy-drift diagnostic.

        ``cache`` is an optional pre-built :class:`RailsEphemerisCache` spanning
        ``[t0, t1]`` (or wider). Building the spline cache is the dominant per-call
        cost (~0.5 s for a multi-month arc — the planet-state sampling); a
        multiple-shooting solve evaluates many leg propagations over the SAME
        itinerary window, so the shooter builds ONE cache spanning the whole
        itinerary and passes it to every leg, turning the per-leg cost into the
        bare integration. When ``cache is None`` a per-call cache is built (the
        rung / golden-gate path, where each call is independent).

        Divergence is a first-class outcome (design §3, mirror ``correct.py``):
        a pathological seed (e.g. one that grazes a perturber's softened core)
        forces IAS15 into ever-smaller steps. Rather than spin forever, the
        integration is budgeted by ``max_steps`` and ``max_wall_sec``; hitting
        either returns ``converged=False`` with the last finite state — the honest
        DIVERGENT signal the SILVER rung grades as an ARTIFACT.

        When ``with_stm`` is True the returned :class:`NBodyArc` carries the 6x6
        per-leg state-transition matrix ``Phi[:, k] = d(final_state)/d(init_k)``
        (#388), co-integrated via REBOUND's native first-order variational
        particles on the same adaptive step. When False the call is byte-identical
        to the prior behaviour (no variations are added) and ``stm`` stays ``None``.
        On a budget-out / divergence the best-effort variational state at the last
        finite step is returned (``converged=False``), never a crash.
        """
        import time

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
            assert ephem is not None
            if cache is None:
                cache = RailsEphemerisCache(bodies, ephem, t0_sec, t1_sec)
            _install_rails_forces(sim, bodies, cache)

        # First-order variational particles for the STM (#388). The spacecraft is
        # particle index 1 (Sun is 0). Six variations seed the six unit basis
        # perturbations of the spacecraft initial state (e_0..e_2 position,
        # e_3..e_5 velocity); each co-integrates with the state on the SAME
        # adaptive IAS15 step, and its final spacecraft state is column k of Phi.
        i_sc = 1
        variations: list[Any] = []
        if with_stm:
            for k in range(6):
                var = sim.add_variation(order=1)
                vp = var.particles[i_sc]
                vp.x = 1.0 if k == 0 else 0.0
                vp.y = 1.0 if k == 1 else 0.0
                vp.z = 1.0 if k == 2 else 0.0
                vp.vx = 1.0 if k == 3 else 0.0
                vp.vy = 1.0 if k == 4 else 0.0
                vp.vz = 1.0 if k == 5 else 0.0
                variations.append(var)

        e0 = float(sim.energy())
        t_target = float(t1_sec)
        # Walk toward the target in coarse time chunks so the step / wall budgets
        # can be checked between chunks (REBOUND has no native step cap). With no
        # perturbers (two-body) this completes in one chunk — the golden gates are
        # unaffected. With perturbers the chunk count scales with the arc length
        # (~1 chunk per 30 days, clamped to [4, 200]): a fixed 200-chunk cadence
        # adds ~10x Python-call overhead on short multi-month legs for an identical
        # result (the chunking is only a budget-check cadence — IAS15 adapts its own
        # internal substeps), which is the difference between a tractable and an
        # 8-min-timeout multiple-shooting Jacobian.
        if not bodies:
            n_chunks = 1
        else:
            from cyclerfinder.core.constants import SECONDS_PER_DAY as _SPD

            arc_days = abs(t_target - float(sim.t)) / _SPD
            n_chunks = int(min(200, max(4, arc_days / 30.0)))
        chunk = (t_target - float(sim.t)) / n_chunks if n_chunks else 0.0
        wall_start = time.monotonic()
        converged = True
        target = float(sim.t)
        try:
            for _ in range(n_chunks):
                target = t_target if _ == n_chunks - 1 else float(sim.t) + chunk
                sim.integrate(target)
                sc_p = sim.particles[1]
                if not (np.isfinite(sc_p.x) and np.isfinite(sc_p.y) and np.isfinite(sc_p.z)):
                    converged = False
                    break
                # REBOUND 5.0: the cumulative integration-step count is the
                # `steps_done` attribute; `steps` is a (callable) method, so
                # comparing it to an int raises and would be silently swallowed by
                # the divergence handler — always flagging converged=False.
                if sim.steps_done > max_steps or (time.monotonic() - wall_start) > max_wall_sec:
                    converged = False
                    break
        except Exception:
            # REBOUND raises on a genuine integration breakdown; treat as divergent.
            converged = False

        e1 = float(sim.energy())
        energy_rel_drift = (e1 - e0) / abs(e0) if e0 != 0.0 and np.isfinite(e1) else 0.0

        sc = sim.particles[1]
        r_km = np.array([sc.x, sc.y, sc.z], dtype=np.float64)
        v_km_s = np.array([sc.vx, sc.vy, sc.vz], dtype=np.float64)

        stm: NDArray[np.float64] | None = None
        if with_stm:
            phi = np.empty((6, 6), dtype=np.float64)
            for k, var in enumerate(variations):
                vp = var.particles[i_sc]
                phi[:, k] = [vp.x, vp.y, vp.z, vp.vx, vp.vy, vp.vz]
            stm = phi

        return NBodyArc(
            r_km=r_km,
            v_km_s=v_km_s,
            t1_sec=float(sim.t),
            energy_rel_drift=float(energy_rel_drift),
            # Reserved: RestrictedNBody seeds at the exact (r0, v0) and never
            # anchors against an ephemeris, so it has no anchor error of its own.
            # Callers compute closure themselves (see NBodyArc.anchor_err_km doc;
            # rung.py derives terminal_closure_km from the returned state).
            anchor_err_km=0.0,
            integrator_accuracy=float(accuracy),
            bodies=bodies,
            converged=converged,
            stm=stm,
        )


def _install_rails_forces(sim: object, bodies: tuple[str, ...], cache: RailsEphemerisCache) -> None:
    """Attach the rails third-body perturbation to a REBOUND simulation.

    The callback reads each perturber's heliocentric position at the integrator's
    *current* sub-step time (``sim.t``) from the spline rails cache (built once
    from the shared DE440 reader) — the per-force-eval read that is the harness's
    cost (design Risk 3), made tractable by interpolation. Sun-frame indirect term
    included so the heliocentric (non-inertial) frame is consistent.

    **Variational (STM) coupling.** REBOUND's native first-order variational
    particles co-integrate only the *in-sim* gravity (the central Sun); the rails
    perturbation lives in this Python ``additional_forces`` callback, which REBOUND
    does NOT differentiate. So when variational particles are present (``with_stm``)
    this callback must also apply the linearised perturber force to them, otherwise
    the STM is correct on cruise legs but qualitatively wrong near a perturber. The
    linearisation is the gravity-gradient tensor of the *direct* term,
    ``G_p = mu_p * (3 d d^T / |d|^5 - I/|d|^3)`` with ``d = r_p - r_sc`` (the
    indirect ``-r_p/|r_p|^3`` term is spacecraft-independent, so it drops out), and
    in the softened core (``|d| <= d_safe``, ``d_safe`` constant) it degenerates to
    ``-mu_p/d_safe^3 * I`` — matching the gradient of the clamped direct force so
    the analytic STM agrees with a finite-difference of THIS callback. The summed
    tensor ``G`` is applied as ``delta_a = G @ delta_r`` to each variational
    particle's spacecraft component. Built only when ``N_var_config > 0`` so the
    non-STM path keeps its exact prior cost.
    """
    from cyclerfinder.core.constants import PLANETS

    def additional_forces(reb_sim_pointer: object) -> None:
        reb_sim = reb_sim_pointer.contents  # type: ignore[attr-defined]
        t_sec = float(reb_sim.t)
        sc = reb_sim.particles[1]
        r = np.array([sc.x, sc.y, sc.z], dtype=np.float64)
        n_var = int(reb_sim.N_var_config)
        grad = np.zeros((3, 3), dtype=np.float64) if n_var else None
        ax = ay = az = 0.0
        for body in bodies:
            pdata = PLANETS[body]
            mu_p = pdata.mu_km3_s2
            r_p = cache.position(body, t_sec)
            d = r_p - r
            d_norm = float(np.linalg.norm(d))
            # Soften the point-mass singularity at the safe-flyby periapsis: inside
            # the SOI a real flyby is a B-plane-targeted patched-conic turn, NOT a
            # heliocentric point-mass integration (design §2: flyby bending lives
            # at the patch points). Clamping |d| to the safe periapsis keeps the
            # heliocentric arc finite; a trajectory that would dive inside it is a
            # divergent-seed signal the rung surfaces, not a NaN crash.
            d_safe = pdata.radius_eq_km + pdata.safe_alt_km
            d_eff = max(d_norm, d_safe)
            d3 = d_eff**3
            rp3 = float(np.linalg.norm(r_p)) ** 3
            acc = mu_p * (d / d3 - r_p / rp3)
            ax += float(acc[0])
            ay += float(acc[1])
            az += float(acc[2])
            if grad is not None:
                # Gravity-gradient tensor of the direct term wrt the spacecraft
                # position, consistent with the (possibly softened) force above.
                if d_norm > d_safe:
                    grad += mu_p * (3.0 * np.outer(d, d) / (d_norm**5) - np.eye(3) / (d_norm**3))
                else:
                    grad += -mu_p / (d_safe**3) * np.eye(3)
        sc.ax += ax
        sc.ay += ay
        sc.az += az
        if grad is not None:
            pv = reb_sim.particles_var
            for i in range(n_var):
                idx = int(reb_sim.var_config[i].index)
                vp = pv[idx + 1]  # spacecraft (i_sc=1) variational component of config i
                dr = np.array([vp.x, vp.y, vp.z], dtype=np.float64)
                da = grad @ dr
                vp.ax += float(da[0])
                vp.ay += float(da[1])
                vp.az += float(da[2])

    # The setter wraps `func` in a CFUNCTYPE and stores it in the simulation's
    # own `_afp` slot, so the callback is kept alive for us (REBOUND 5.0); we do
    # NOT attach our own attribute — Simulation.__setattr__ rejects unknown names.
    sim.additional_forces = additional_forces  # type: ignore[attr-defined]
    sim.force_is_velocity_dependent = 0  # type: ignore[attr-defined]


__all__ = ["NBodyArc", "Propagator", "RestrictedNBody"]
