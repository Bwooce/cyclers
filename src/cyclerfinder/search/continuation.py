"""Circular -> ephemeris CONTINUATION driver (Russell 2004 dissertation §5.4).

This is the #143 rank-4 build: the *continuation driver* that walks a cycler
closed on the circular-coplanar planet model out to the true ephemeris, the piece
the deep-dive (``docs/notes/2026-06-07-russell-2004-continuation-deepdive.md``)
identifies as genuinely missing from the tree.

What is continued
-----------------
Per Russell §5.4.1 / §5.4.8 (deep-dive §1): **the solar-system MODEL fidelity, NOT
the trajectory genome.** The cycler's free variables stay the optimisation
unknowns; the homotopy parameter is the *planet model* that the corrector's
position/phase constraints reference. The schedule (Fig.5.3, p.169) is, verbatim:

0. Circular-coplanar model (J2000 mean ``a, Ω, ω, nu`` for Earth and Mars; ``e=i=0``).
1. Walk towards the J2000 mean **eccentricities** in ``nstep`` equal steps.
2. Walk towards the J2000 mean **inclinations** in ``nstep`` equal steps.
3. One final step to the accurate ephemeris.

The large elements ``a, Ω, ω, nu`` are **frozen** at their J2000 mean values
throughout (deep-dive failure mode 5, p.169: "otherwise the model changes too
rapidly"); only the two small perturbations ``e`` and ``i`` are ramped, by factors
``λ_e, λ_i ∈ [0, 1]``. Each step re-converges the corrector seeded from the
previous step's solution (classic imbedding, §5.4.1).

The ``nstep`` ladder
--------------------
``nstep = 3^(steploop-1)`` -> ``{1, 3, 9, 27, 81, 243}`` (Fig.5.4, p.170). There is
**no universal step size** (deep-dive failure mode 6, p.171): the best ``nstep``
depends on the parent and launch window, so the driver runs several rungs and keeps
the lowest-residual result. We skip 243 by default (wall-time; the deep-dive notes
the full ladder buys ~50% improvement at large cost and ``nstep=5`` is "an excellent
compromise"). Skipped rungs are recorded in the audit trail, never silently capped.

Conditioning near e=i=0 (deep-dive §3 / failure mode 4)
-------------------------------------------------------
Russell switches the *analytic Jacobian* to a non-singular ``β`` element set because
the classic-element partials blow up at ``e→0, i→0``. Our inner solve is
:func:`cyclerfinder.search.free_return.free_return_correct`, whose genome is
``(a, e, t0)`` with the spacecraft ``e`` bounded away from zero
(``0 < e < 0.95``) and whose residual is built from heliocentric *longitudes* and a
Mars-reach margin — it never forms a classic ``(ω, nu)`` partial at the circular
endpoint. The singularity Russell guards against is in the *planet* model elements,
which here are supplied by the backend as Cartesian states, not differentiated. So
the ``(a, e, t0)`` genome is already non-singular for this homotopy and needs no
local reparameterisation; this is documented rather than re-derived.

The intermediate models (this module, no core edit)
---------------------------------------------------
:class:`_RampedElementsBackend` is a thin backend (conforms to the
``Ephemeris.state(body, t_sec)`` duck type) that places each planet on a *Keplerian
ellipse* with frozen J2000 mean ``a, Ω, ω`` and ramped ``e = λ_e·e_J2000``,
``i = λ_i·i_J2000``. It reuses :func:`cyclerfinder.core.kepler.coe_to_rv` for the
in-plane state and rotates by ``(i, Ω)``. The phase is anchored so that at
``λ_e = λ_i = 0`` and at ``t_sec = 0`` the body sits at ``θ = 0`` — i.e. the
``λ_e = λ_i = 0`` backend reproduces ``Ephemeris('circular')`` bit-identically (the
fast mechanics gate). The final step swaps to ``Ephemeris('astropy')`` (DE440).

NOTE: ``Ephemeris('circular')`` ignores ``Ω, ω`` and starts every body at
``θ = 0``; the J2000 mean ``Ω, ω, nu`` are therefore NOT used at ``λ = 0`` (the seed
is the circular model the seed cycler was closed on). They become meaningful only as
the true ephemeris is approached and are carried for completeness / future use, but
to preserve the bit-identical ``λ=0`` gate the ramped ellipse phase is built on the
same mean-motion clock as the circular backend (see ``_RampedElementsBackend``).

Pure: depends only on core/constants, core/kepler, core/ephemeris, search/free_return.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import cos, pi, radians, sin

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS, SECONDS_PER_DAY
from cyclerfinder.core.ephemeris import Ephemeris, Vec3
from cyclerfinder.core.kepler import coe_to_rv
from cyclerfinder.search.free_return import FreeReturnClosureResult, free_return_correct

DAY_S = SECONDS_PER_DAY

# Russell Table 5.4 (p.169): J2000 MEAN eccentricities and inclinations — the
# ramp ENDPOINTS. Sourced (Russell 2004 dissertation Table 5.4); used as the
# λ=1 targets. Earth defines the ecliptic so its J2000 mean inclination is 0.
_J2000_MEAN_E: dict[str, float] = {
    "E": 1.67086171540e-2,
    "M": 9.34006199474e-2,
}
_J2000_MEAN_I_DEG: dict[str, float] = {
    "E": 0.0,
    "M": 1.84972647778e0,
}
# Frozen large angles (Table 5.4): node Ω and argument of periapsis ω.
_J2000_MEAN_LAN_DEG: dict[str, float] = {"E": 0.0, "M": 4.95655237028e1}
_J2000_MEAN_ARGP_DEG: dict[str, float] = {"E": 1.02937348083e2, "M": 2.86494710278e2}
# J2000 true anomaly nu (Table 5.4). The frozen J2000 mean LONGITUDE the ramped
# model must carry so its planet phase tracks the true ephemeris (otherwise the
# final mean-element -> DE440 step is seeded ~100° off in absolute longitude and
# the corrector diverges). At e=i=0, Ω+ω+nu = the mean longitude; we phase-anchor
# on Ω+ω+nu directly (Table 5.4) so the intermediate model's longitudes land on
# DE440's at J2000 (Earth: 100.47° vs DE440 100.38°; Mars: 355.43° vs 359.45°).
_J2000_MEAN_NU_DEG: dict[str, float] = {"E": -2.47089957222e0, "M": 1.93730406472e1}


def _j2000_mean_longitude_rad(body: str) -> float:
    """Frozen J2000 mean longitude ``Ω + ω + nu`` (rad) from Russell Table 5.4."""
    deg = (
        _J2000_MEAN_LAN_DEG.get(body, 0.0)
        + _J2000_MEAN_ARGP_DEG.get(body, 0.0)
        + _J2000_MEAN_NU_DEG.get(body, 0.0)
    )
    return radians(deg % 360.0)


# The nstep ladder nstep = 3^(steploop-1) (Fig.5.4 step 7.1, p.170). 243 skipped
# by default for wall-time (deep-dive §1: full ladder buys ~50% at large cost).
LADDER: tuple[int, ...] = (1, 3, 9, 27, 81)


class _RampedElementsBackend:
    """Planet states on a Keplerian ellipse with ramped ``e`` and ``i``.

    Frozen J2000 mean ``a`` (from :data:`PLANETS`) and the circular mean-motion
    clock; ``e = lam_e · e_J2000`` and ``i = lam_i · i_J2000`` (Russell Table 5.4
    endpoints). At ``lam_e = lam_i = 0`` the state is byte-identical to
    :class:`cyclerfinder.core.ephemeris._CircularBackend` (the bit-identical λ=0
    gate): the in-plane circular state with the body at ``θ = 0`` at ``t = 0``.

    Eccentricity is introduced by treating the circular mean anomaly ``M = n·t``
    as the ellipse's mean anomaly (same a, same period, same θ=0 reference), so the
    transition from circle to ellipse is continuous in ``lam_e`` (a small e bends
    the circle into a low-eccentricity ellipse with periapsis along ``+x``).
    Inclination tilts the orbit plane about the ``+x`` node line by
    ``R_z(Ω) R_x(i)`` (the same convention as the inclined-circular backend), so
    ``lam_i`` lifts the plane continuously off the ecliptic.
    """

    def __init__(self, lam_e: float, lam_i: float, lam_p: float = 1.0) -> None:
        self.lam_e = float(lam_e)
        self.lam_i = float(lam_i)
        # λ_p ramps the per-body J2000 mean-LONGITUDE phase from 0 (the circular
        # seed's θ=0-at-J2000 frame) to 1 (aligned with the true ephemeris). It is
        # a SEPARATE leading ramp, NOT coupled to e: the phase offset is a per-body
        # frame/epoch alignment (different bodies need different offsets, so a
        # single t0 shift cannot do it) and ~100° at Earth, far too large to land
        # in one e-step. Ramping it first and gradually (its own nstep sub-steps)
        # keeps each step a small perturbation the corrector can track, then the
        # e- and i-ramps proceed from the phase-aligned model.
        self.lam_p = float(lam_p)

    def _elements(self, body: str) -> tuple[float, float, float]:
        e = self.lam_e * _J2000_MEAN_E.get(body, 0.0)
        i_rad = radians(self.lam_i * _J2000_MEAN_I_DEG.get(body, 0.0))
        lan_rad = radians(_J2000_MEAN_LAN_DEG.get(body, 0.0))
        return e, i_rad, lan_rad

    @staticmethod
    def _mean_to_true(mean_anom: float, e: float) -> float:
        """Mean -> true anomaly (Kepler's equation, Newton). e small here."""
        m = mean_anom % (2.0 * pi)
        if e == 0.0:
            return m
        ecc = m if e < 0.8 else pi
        for _ in range(60):
            f = ecc - e * sin(ecc) - m
            fp = 1.0 - e * cos(ecc)
            d = f / fp
            ecc -= d
            if abs(d) < 1.0e-14:
                break
        beta = e / (1.0 + (1.0 - e * e) ** 0.5)
        return float(ecc + 2.0 * np.arctan2(beta * sin(ecc), 1.0 - beta * cos(ecc)))

    def state(self, body: str, t_sec: float) -> tuple[Vec3, Vec3]:
        planet = PLANETS[body]
        a_km = planet.sma_au * AU_KM
        n_rad_s = planet.mean_motion_deg_day * (pi / 180.0) / SECONDS_PER_DAY
        e, i_rad, lan_rad = self._elements(body)
        # Phase-anchor: add the frozen J2000 mean longitude ramped by λ_p so the
        # ramped model's planet phase tracks the true ephemeris. At λ_p=0 the
        # offset is zero -> bit-identical to Ephemeris('circular'); at λ_p=1 the
        # body sits at its J2000 mean longitude at t=0, aligned with DE440 (so the
        # final ephemeris step is a SMALL perturbation, not a ~100° phase jump).
        mean_anom = n_rad_s * t_sec + self.lam_p * _j2000_mean_longitude_rad(body)
        if e == 0.0:
            # Bit-identical circular path (no Kepler solve, no perifocal scaling).
            theta = mean_anom
            cos_t, sin_t = cos(theta), sin(theta)
            speed = a_km * n_rad_s
            r_pf = np.array([a_km * cos_t, a_km * sin_t, 0.0], dtype=np.float64)
            v_pf = np.array([-speed * sin_t, speed * cos_t, 0.0], dtype=np.float64)
        else:
            nu = float(self._mean_to_true(mean_anom, e))
            r_pf, v_pf = coe_to_rv(a_km, e, nu, MU_SUN_KM3_S2)
        if i_rad == 0.0:
            # No tilt -> the node is undefined / irrelevant; the in-plane state IS
            # the ecliptic state. This is the λ_i=0 branch (the whole e-ramp phase)
            # and keeps it bit-identical to Ephemeris('circular') for any frozen Ω.
            return r_pf, v_pf
        rot = self._tilt(i_rad, lan_rad)
        return (
            np.asarray(rot @ r_pf, dtype=np.float64),
            np.asarray(rot @ v_pf, dtype=np.float64),
        )

    @staticmethod
    def _tilt(i_rad: float, lan_rad: float) -> NDArray[np.float64]:
        """Tilt the orbit plane about the node line WITHOUT rotating the in-plane
        phase: ``R_z(Ω) R_x(i) R_z(-Ω)``.

        Unlike the inclined-circular backend's ``R_z(Ω) R_x(i)`` (which places the
        body ON the ascending node at ``t=0``), this conjugated form is the
        IDENTITY at ``i=0`` for any ``Ω`` — so as ``λ_i`` ramps off zero the plane
        lifts continuously off the ecliptic with NO discontinuous longitude jump
        from activating a frozen non-zero node (Mars Ω=49.6°). The tilt axis is the
        line at ecliptic longitude ``Ω``; the heliocentric longitude of the body is
        preserved at ``i=0``, which is exactly what the continuation needs (ramp i,
        not phase).
        """
        ci, si = cos(i_rad), sin(i_rad)
        cl, sl = cos(lan_rad), sin(lan_rad)
        rx = np.array([[1.0, 0.0, 0.0], [0.0, ci, si], [0.0, -si, ci]], dtype=np.float64)
        rz = np.array([[cl, -sl, 0.0], [sl, cl, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)
        rz_neg = np.array([[cl, sl, 0.0], [-sl, cl, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)
        return rz @ rx @ rz_neg


def ramped_ephemeris(lam_e: float, lam_i: float, lam_p: float = 0.0) -> Ephemeris:
    """Build an :class:`Ephemeris`-compatible provider at ramp factors.

    Returns a real :class:`Ephemeris` whose backend is swapped for a
    :class:`_RampedElementsBackend`, so it is a drop-in for any corrector that
    consumes ``ephem.state(body, t_sec)``. NO core file is edited: the backend is a
    private helper in this module and is injected post-construction.

    ``(0.0, 0.0, 0.0)`` reproduces ``Ephemeris('circular')`` bit-identically.
    ``lam_p`` ramps the J2000 mean-longitude phase anchor (default 0).
    """
    ephem = Ephemeris("circular")
    ephem._backend = _RampedElementsBackend(
        lam_e, lam_i, lam_p
    )  # deliberate backend injection (no core edit)
    ephem._model = f"ramped(p={lam_p:.6g},e={lam_e:.6g},i={lam_i:.6g})"
    return ephem


# ---------------------------------------------------------------------------
# Audit trail
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ContinuationStep:
    """One homotopy step's outcome (the per-step audit record)."""

    phase: str  # "p-ramp" | "e-ramp" | "i-ramp" | "ephemeris"
    lam_e: float
    lam_i: float
    lam_p: float
    model: str
    max_residual_kms: float
    a_au: float
    e: float
    t0_sec: float
    converged: bool


@dataclass(frozen=True)
class ContinuationRung:
    """One ``nstep`` rung: its step trail and final result."""

    nstep: int
    steps: tuple[ContinuationStep, ...]
    final: FreeReturnClosureResult
    completed: bool  # reached the ephemeris step without diverging


@dataclass(frozen=True)
class ContinuationResult:
    """Continuation driver outcome (deep-dive §5 ``optimise_cell_ephemeris``).

    Attributes
    ----------
    best_final:
        The lowest-residual final (true-ephemeris) :class:`FreeReturnClosureResult`
        across all completed rungs — the headline.
    winning_nstep:
        Which ladder rung produced ``best_final``.
    rungs:
        Per-rung audit trail (every step's residual + which rung won).
    ladder:
        The ``nstep`` values actually attempted.
    skipped:
        Ladder values deliberately NOT run (e.g. 243), recorded — never silent.
    """

    best_final: FreeReturnClosureResult
    winning_nstep: int
    rungs: tuple[ContinuationRung, ...]
    ladder: tuple[int, ...]
    skipped: tuple[int, ...] = field(default_factory=tuple)

    @property
    def ballistic_within_tol(self) -> bool:
        """``best_final`` converged at the true-ephemeris step (residual < tol)."""
        return self.best_final.converged


def _ramp_schedule(nstep: int) -> list[tuple[str, float, float, float]]:
    """Build the (phase, λ_e, λ_i, λ_p) sequence for one rung.

    Russell's schedule is e-ramp then i-ramp (Fig.5.3). We PREPEND a phase ramp
    (``λ_p: 0 -> 1``) because our circular seed lives in a θ=0-at-J2000 frame while
    the true ephemeris does not, and that per-body ~100° frame offset cannot be
    absorbed by the corrector's single t0 (different bodies need different shifts).
    The phase ramp is the epoch/frame-alignment that Russell folds into his
    Fig.5.4 seeding step ("propagate the simple model until the parent's beginning
    phase angle is achieved"); here it is an explicit homotopy leg so each step
    stays small. Then ``nstep`` e-steps (``λ_e: 0 -> 1`` at λ_p=1, λ_i=0) and
    ``nstep`` i-steps (``λ_i: 0 -> 1`` at λ_p=λ_e=1). The final true-ephemeris step
    is appended by the driver.
    """
    sched: list[tuple[str, float, float, float]] = []
    for k in range(1, nstep + 1):
        sched.append(("p-ramp", 0.0, 0.0, k / nstep))
    for k in range(1, nstep + 1):
        sched.append(("e-ramp", k / nstep, 0.0, 1.0))
    for k in range(1, nstep + 1):
        sched.append(("i-ramp", 1.0, k / nstep, 1.0))
    return sched


def continuation_correct(
    t0_seed_sec: float,
    a_seed_au: float,
    e_seed: float,
    period_sec: float,
    *,
    bodies: tuple[str, str] = ("E", "M"),
    mu: float = MU_SUN_KM3_S2,
    tol_kms: float = 0.1,
    ladder: tuple[int, ...] = LADDER,
    final_ephemeris: Ephemeris | None = None,
) -> ContinuationResult:
    """Walk a circular-coplanar free-return closure out to the true ephemeris.

    Implements the deep-dive §5 ``optimise_cell_ephemeris`` sketch: seed from the
    circular-coplanar ``(a, e, t0)`` closure, then for each ``nstep`` in ``ladder``
    ramp ``e`` (then ``i``) in ``nstep`` equal steps and finally step to
    ``final_ephemeris`` (default ``Ephemeris('astropy')`` = DE440). Each step
    re-converges :func:`free_return_correct` seeded from the previous step's
    ``(a, e, t0)``. Keeps the lowest-residual final result across rungs.

    The seed ``(t0, a, e)`` is the SOURCED / circular-closed input (a constraint);
    the per-body V_inf and ToFs on ``best_final`` EMERGE and are evidence (the
    golden-rule separation inherited from :mod:`free_return`).

    Parameters
    ----------
    final_ephemeris:
        The λ=1 true-ephemeris target. Defaults to a fresh
        ``Ephemeris('astropy')``. Injectable for tests (e.g. a ramped backend at
        ``(1, 1)`` to isolate the ramp from the DE440 swap).
    ladder:
        ``nstep`` rungs to try. The default skips 243 (recorded in
        ``ContinuationResult.skipped``).
    """
    skipped = tuple(n for n in (1, 3, 9, 27, 81, 243) if n not in ladder)
    if final_ephemeris is None:
        final_ephemeris = Ephemeris("astropy")

    def _correct(seed: FreeReturnClosureResult | None, ephem: Ephemeris) -> FreeReturnClosureResult:
        if seed is None:
            t0, a, e = t0_seed_sec, a_seed_au, e_seed
        else:
            t0, a, e = seed.t0_sec, seed.a_au, seed.e
        return free_return_correct(
            t0_seed_sec=t0,
            a_seed_au=a,
            e_seed=e,
            period_sec=period_sec,
            ephem=ephem,
            bodies=bodies,
            mu=mu,
            tol_kms=tol_kms,
        )

    rungs: list[ContinuationRung] = []
    for nstep in ladder:
        steps: list[ContinuationStep] = []
        current: FreeReturnClosureResult | None = None
        completed = True
        for phase, lam_e, lam_i, lam_p in _ramp_schedule(nstep):
            ephem = ramped_ephemeris(lam_e, lam_i, lam_p)
            current = _correct(current, ephem)
            steps.append(
                ContinuationStep(
                    phase=phase,
                    lam_e=lam_e,
                    lam_i=lam_i,
                    lam_p=lam_p,
                    model=ephem.model,
                    max_residual_kms=current.max_residual_kms,
                    a_au=current.a_au,
                    e=current.e,
                    t0_sec=current.t0_sec,
                    converged=current.converged,
                )
            )
            if not np.isfinite(current.max_residual_kms):
                completed = False
                break
        if completed:
            final = _correct(current, final_ephemeris)
            steps.append(
                ContinuationStep(
                    phase="ephemeris",
                    lam_e=1.0,
                    lam_i=1.0,
                    lam_p=1.0,
                    model=final_ephemeris.model,
                    max_residual_kms=final.max_residual_kms,
                    a_au=final.a_au,
                    e=final.e,
                    t0_sec=final.t0_sec,
                    converged=final.converged,
                )
            )
        else:
            # Diverged mid-ramp: the last finite-or-not step is the rung's "final".
            final = current if current is not None else _correct(None, final_ephemeris)
        rungs.append(
            ContinuationRung(nstep=nstep, steps=tuple(steps), final=final, completed=completed)
        )

    # Keep the lowest-residual final across COMPLETED rungs; fall back to the best
    # of all rungs if none completed (so a result object always exists).
    completed_rungs = [r for r in rungs if r.completed]
    pool = completed_rungs if completed_rungs else rungs
    best = min(pool, key=lambda r: r.final.max_residual_kms)
    return ContinuationResult(
        best_final=best.final,
        winning_nstep=best.nstep,
        rungs=tuple(rungs),
        ladder=tuple(ladder),
        skipped=skipped,
    )
