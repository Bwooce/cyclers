"""V3 independent-integrator 3D gauntlet for CR3BP periodic orbits (#306).

Spec reference
--------------
* §14 V3 — independent-cross-check at scale: the candidate's defining
  dynamical signature must survive re-evaluation under an INDEPENDENT
  integrator ARCHITECTURE.

Why this module exists
----------------------
The headline :func:`cyclerfinder.data.validation.v3_3d.run_v3_3d` is
**moontour-specific** (it re-solves Lambert legs in the planet frame) and
does NOT apply to a full-3D CR3BP *periodic orbit* defined by a 6-vector IC
``(x,y,z,ẋ,ẏ,ż)`` + period ``T``. The in-file
:func:`v3_3d.run_v3_periodic_regression` is INTRA-model — it re-runs the same
scipy DOP853 propagator and asserts the same closure, which is a regression
guard, not an independent-architecture cross-check.

This module supplies the genuine independent-architecture V3 for a periodic
orbit: re-propagate the IC in the CR3BP **rotating frame** under **REBOUND
IAS15** (high-order adaptive Gauss-Radau), a fundamentally different family
from the project's scipy DOP853 / Radau, and assert

* (i) the one-period 6D closure ``||X(T) - X(0)||`` that V1 asserts, AND
* (ii) the ≥3-cycle bounded-drift signature that V2 asserts — by agreeing
  with the V2 (DOP853) per-cycle drift series to << the V2 floor.

Agreement ⇒ the periodicity is a real property of the CR3BP model, not a
scipy-DOP853/Radau artefact.

Integrator architecture (the load-bearing risk + its mitigation)
----------------------------------------------------------------
REBOUND integrates inertial N-body. The CR3BP rotating-frame EOM has Coriolis
(``±2 v``) + centrifugal terms that are **velocity-dependent** and must live
in an ``additional_forces`` Python callback (with
``force_is_velocity_dependent = 1``). Per the project memory *REBOUND
variation + custom force gotcha*, REBOUND's native variational particles do
NOT differentiate a Python callback — but V3 here is STATE-ONLY (closure +
drift), so the variational hazard is sidestepped entirely.

The residual risk — IAS15 step/tolerance behaviour under a velocity-dependent
callback — is gated by the MANDATORY parity test
(:func:`tests.data.test_v3_3d_periodic`): IAS15 must match DOP853 to a tight
tolerance over one period on a known CR3BP orbit BEFORE any verdict is
trusted. The parity test passes for the C21 candidate at ~1e-11 nondim, so
the preferred **mode (A) rotating-frame callback** is used. A documented
**mode (B) inertial two-primary fallback** exists for the contingency that
(A) ever fails parity on a future system; it is selected via ``mode`` /
``prefer_rotating_callback``.

When REBOUND is not importable, the helper falls back to scipy **LSODA**
(multistep BDF — a distinct family from the project's single-step DOP853) and
records the honest label, mirroring
:func:`v3_3d._ias15_propagate_planet_frame`.

Discipline
----------
* NO catalogue writeback. A V3 PASS is necessary-not-sufficient — V4 + V5 +
  the literature miss still gate.
* The V3 verdict is whatever the math says. Specific drift numbers are NOT
  golden; integrator AGREEMENT is.
* The parity test GATES everything: if IAS15 cannot match DOP853 on a known
  orbit, the V3 verdict is meaningless and that must be surfaced honestly.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.data.validation.v2_3d import V2Verdict3D

V3_PERIODIC_CLOSURE_FLOOR_NONDIM: Final[float] = 1.0e-7
"""One-period IAS15 closure floor in nondim units.

Matches the :func:`v3_3d.run_v3_periodic_regression` placeholder
(~38 m in Earth-Moon units). This is an AGREEMENT / closure bar for the
independent integrator, NOT the spec §14 V1 1 m/s floor — it asserts IAS15
reproduces the same one-period closure the corrector-clean IC has under
DOP853 (typically ~1e-12 nondim for a converged member)."""

V3_PERIODIC_DRIFT_AGREEMENT_FLOOR_KMS: Final[float] = 1.0
"""IAS15-vs-DOP853 per-cycle drift agreement floor in km.

The two integrator architectures must agree on the per-cycle bounded-drift
series to within 1 km — << the 50,000 km V2 same-model drift floor. If they
agree this tightly, the V2 bounded-drift signature is integrator-architecture
independent (a real model property, not a DOP853 artefact)."""

V3_PERIODIC_N_CYCLES_MIN: Final[int] = 3
"""Spec §14 V2-ballistic minimum carried into V3: ``>= 3`` cycles to even
report a drift-agreement verdict (V3 compares against V2's ≥3-cycle series)."""

#: Internal helper modes (mode-A preferred, mode-B documented fallback).
_MODE_ROTATING_CALLBACK: Final[str] = "rotating_callback"
_MODE_INERTIAL_TWOBODY: Final[str] = "inertial_twobody"


@dataclass(frozen=True)
class V3PeriodicVerdict3D:
    """Frozen V3 verdict for a 3D CR3BP periodic-orbit candidate.

    Attributes
    ----------
    candidate_id:
        Identifier carried for the audit trail.
    integrator:
        Human-readable integrator-architecture label, e.g.
        ``"REBOUND IAS15 (rotating callback)"``,
        ``"REBOUND IAS15 (inertial two-body)"``, or
        ``"scipy LSODA fallback (rotating callback)"``. The fallback is
        honest, not silent.
    closure_residual_nondim_ias15:
        One-period closure ``||X(T) - X(0)||`` under IAS15, nondim.
    closure_residual_kms_ias15:
        The same converted to km/s via the system velocity unit
        (``l_km / t_s``). NOTE: this is a 6D state-norm scaled by the
        velocity unit, reported for human scale only; the gate is on the
        nondim closure.
    per_cycle_drift_kms_ias15:
        Cumulative position drift at each cycle boundary under IAS15, km
        (``||X_pos((k+1)T) - X_pos(0)||``), mirroring
        :attr:`V2Verdict3D.per_cycle_drift_kms`.
    per_cycle_drift_kms_dop853:
        The V2 (DOP853) per-cycle drift series, sliced to the cycles V3
        propagated, for direct comparison.
    drift_agreement_kms:
        ``max_k |ias15[k] - dop853[k]|`` — the headline integrator-agreement
        number.
    closure_floor_nondim:
        The nondim closure bar this verdict was held against.
    drift_agreement_floor_kms:
        The km drift-agreement bar this verdict was held against.
    n_cycles_propagated:
        Cycles IAS15 actually completed.
    converged_at_each_return:
        Whether IAS15 completed every requested cycle.
    passes_v3:
        ``converged_at_each_return AND n_cycles_propagated >= 3 AND
        closure_residual_nondim_ias15 <= closure_floor_nondim AND
        drift_agreement_kms <= drift_agreement_floor_kms``. The headline.
    degenerate_planar:
        Carried from the candidate (the IC's planar degeneracy flag, if the
        caller supplies it). Diagnostic; does not veto.
    notes:
        Free-form audit string.
    """

    candidate_id: str
    integrator: str
    closure_residual_nondim_ias15: float
    closure_residual_kms_ias15: float
    per_cycle_drift_kms_ias15: tuple[float, ...]
    per_cycle_drift_kms_dop853: tuple[float, ...]
    drift_agreement_kms: float
    closure_floor_nondim: float
    drift_agreement_floor_kms: float
    n_cycles_propagated: int
    converged_at_each_return: bool
    passes_v3: bool
    degenerate_planar: bool
    notes: str = ""


def _ias15_propagate_cr3bp_rotating(
    state0_nondim: NDArray[np.float64],
    t_nondim: float,
    mu: float,
    *,
    mode: str = _MODE_ROTATING_CALLBACK,
) -> tuple[NDArray[np.float64], str]:
    """Propagate a 6D CR3BP rotating-frame state for ``t_nondim`` under REBOUND IAS15.

    Returns ``(state_f_nondim, integrator_label)``.

    Mode (A) ``"rotating_callback"`` (preferred, truly independent): a single
    massless particle, no in-sim gravity (``sim.G = 0``), and an
    ``additional_forces`` callback that applies the full CR3BP rotating-frame
    acceleration

        a = (x + 2 ẏ - (1-µ)(x+µ)/r1³ - µ(x-1+µ)/r2³,
             y - 2 ẋ - (1-µ)y/r1³    - µ y/r2³,
                       - (1-µ)z/r1³    - µ z/r2³)

    identical term-for-term to :func:`cyclerfinder.core.cr3bp.cr3bp_eom`. The
    Coriolis terms are velocity-dependent, so
    ``force_is_velocity_dependent = 1`` is set. IAS15's Gauss-Radau stepper is
    the independent architecture vs scipy DOP853.

    Mode (B) ``"inertial_twobody"`` (documented fallback): realise the
    rotating IC into the CR3BP's own nondim inertial frame (synodic→inertial
    rotation at ω = 1), integrate REBOUND native two-primary gravity (Earth at
    ``-µ`` mass ``1-µ``, Moon at ``1-µ`` mass ``µ``, on the CR3BP circular
    orbit, ``sim.G = 1``), then rotate the terminal state back to the rotating
    frame. No Python force callback (cleaner numerically) but it re-imposes the
    CR3BP circular-primary assumption, so it is still same-model.

    When ``rebound`` is unimportable, falls back to scipy LSODA (multistep BDF,
    a distinct family from DOP853) integrating :func:`cr3bp.cr3bp_eom` directly,
    with an honest ``"scipy LSODA fallback (...)"`` label.
    """
    state0 = np.asarray(state0_nondim, dtype=np.float64)
    if state0.shape != (6,):
        raise ValueError(f"state0_nondim must be 6D; got shape {state0.shape}")
    if mode not in (_MODE_ROTATING_CALLBACK, _MODE_INERTIAL_TWOBODY):
        raise ValueError(
            f"mode must be {_MODE_ROTATING_CALLBACK!r} or {_MODE_INERTIAL_TWOBODY!r}; got {mode!r}"
        )

    try:
        import rebound
    except ImportError:
        return _lsoda_fallback_cr3bp_rotating(state0, t_nondim, mu, mode=mode)

    if mode == _MODE_ROTATING_CALLBACK:
        return _ias15_rotating_callback(rebound, state0, t_nondim, mu)
    return _ias15_inertial_twobody(rebound, state0, t_nondim, mu)


def _ias15_rotating_callback(
    rebound: object,
    state0: NDArray[np.float64],
    t_nondim: float,
    mu: float,
) -> tuple[NDArray[np.float64], str]:
    """Mode (A): IAS15 with the CR3BP rotating-frame acceleration in a callback."""
    sim = rebound.Simulation()  # type: ignore[attr-defined]
    sim.integrator = "ias15"
    # No in-sim gravity — the whole acceleration is supplied by the callback.
    sim.G = 0.0
    sim.add(
        m=0.0,
        x=float(state0[0]),
        y=float(state0[1]),
        z=float(state0[2]),
        vx=float(state0[3]),
        vy=float(state0[4]),
        vz=float(state0[5]),
    )
    om1 = 1.0 - mu

    def additional_forces(reb_sim_pointer: object) -> None:
        reb_sim = reb_sim_pointer.contents  # type: ignore[attr-defined]
        p = reb_sim.particles[0]
        x, y, z = p.x, p.y, p.z
        vx, vy = p.vx, p.vy  # z-accel has no Coriolis term, so vz is unused
        r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
        r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
        r1c = r1**3
        r2c = r2**3
        p.ax += x + 2.0 * vy - om1 * (x + mu) / r1c - mu * (x - 1.0 + mu) / r2c
        p.ay += y - 2.0 * vx - om1 * y / r1c - mu * y / r2c
        p.az += -om1 * z / r1c - mu * z / r2c

    # REBOUND 5.0 wraps the callback in a CFUNCTYPE in its own slot (kept
    # alive for us). The Coriolis term is velocity-dependent, so this MUST be
    # flagged or IAS15's predictor will treat the force as velocity-free.
    sim.additional_forces = additional_forces
    sim.force_is_velocity_dependent = 1
    sim.integrate(float(t_nondim))
    p = sim.particles[0]
    state_f = np.array([p.x, p.y, p.z, p.vx, p.vy, p.vz], dtype=np.float64)
    return state_f, "REBOUND IAS15 (rotating callback)"


def _rotating_to_inertial(state6: NDArray[np.float64], theta: float) -> NDArray[np.float64]:
    """Synodic (rotating, ω=1) → nondim inertial at rotation angle ``theta``.

    Position: ``r_in = R(theta) r_rot``. Velocity: ``v_in = R(theta)(v_rot +
    omega_cross_r)`` with omega = z-hat (rotating frame spins at unit rate about z).
    """
    c, s = math.cos(theta), math.sin(theta)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)
    r = state6[:3]
    v = state6[3:]
    omega_cross_r = np.array([-r[1], r[0], 0.0], dtype=np.float64)  # z-hat cross r
    r_in = rot @ r
    v_in = rot @ (v + omega_cross_r)
    return np.concatenate([r_in, v_in])


def _inertial_to_rotating(state6: NDArray[np.float64], theta: float) -> NDArray[np.float64]:
    """Nondim inertial → synodic (rotating, ω=1) at rotation angle ``theta``.

    Inverse of :func:`_rotating_to_inertial`: ``r_rot = R(-theta) r_in``;
    ``v_rot = R(-theta) v_in - omega_cross_r``.
    """
    c, s = math.cos(theta), math.sin(theta)
    rot_inv = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)
    r_in = state6[:3]
    v_in = state6[3:]
    r_rot = rot_inv @ r_in
    v_body = rot_inv @ v_in
    omega_cross_r = np.array([-r_rot[1], r_rot[0], 0.0], dtype=np.float64)
    v_rot = v_body - omega_cross_r
    return np.concatenate([r_rot, v_rot])


def _ias15_inertial_twobody(
    rebound: object,
    state0: NDArray[np.float64],
    t_nondim: float,
    mu: float,
) -> tuple[NDArray[np.float64], str]:
    """Mode (B) fallback: realise the rotating IC into the nondim inertial frame.

    Two massive primaries (Earth mass ``1-µ`` at ``-µ``, Moon mass ``µ`` at
    ``1-µ``) on the CR3BP circular orbit, ``sim.G = 1`` (native REBOUND
    gravity — no Python callback). Integrate inertial IAS15, rotate the
    terminal state back to the rotating frame.

    NOTE: this re-imposes the CR3BP circular-primary geometry exactly, so it
    is the same model as mode (A) — it is the numerically-cleaner fallback
    used only if the velocity-dependent callback of mode (A) ever fails the
    parity gate. The primaries themselves are integrated (a 3-body inertial
    simulation), so this is a faithful inertial realisation of the CR3BP.
    """
    sim = rebound.Simulation()  # type: ignore[attr-defined]
    sim.integrator = "ias15"
    sim.G = 1.0
    om1 = 1.0 - mu
    # Primaries on the circular orbit: at theta=0 they lie on the inertial
    # x-axis (Earth at -mu, Moon at 1-mu), moving with the rotating frame
    # (omega = z-hat): inertial velocity v = omega cross r = (-y, x, 0) = (0, x, 0).
    sim.add(m=om1, x=-mu, y=0.0, z=0.0, vx=0.0, vy=-mu, vz=0.0)
    sim.add(m=mu, x=om1, y=0.0, z=0.0, vx=0.0, vy=om1, vz=0.0)
    # Spacecraft (massless) — rotating IC realised into inertial at theta=0.
    sc0 = _rotating_to_inertial(state0, 0.0)
    sim.add(
        m=0.0,
        x=float(sc0[0]),
        y=float(sc0[1]),
        z=float(sc0[2]),
        vx=float(sc0[3]),
        vy=float(sc0[4]),
        vz=float(sc0[5]),
    )
    sim.move_to_com()
    sim.integrate(float(t_nondim))
    sc = sim.particles[2]
    sc_inertial = np.array([sc.x, sc.y, sc.z, sc.vx, sc.vy, sc.vz], dtype=np.float64)
    # The COM shift is constant in the inertial frame; undo it for the
    # spacecraft by re-referencing to the instantaneous primary barycentre.
    # move_to_com shifts everything by a constant velocity*t + offset; the
    # rotating-frame transform is defined about the (fixed) system barycentre
    # at the origin, so we transform the COM-frame state directly — the COM is
    # at the origin after move_to_com, matching the rotating-frame origin.
    state_f = _inertial_to_rotating(sc_inertial, float(t_nondim))
    return state_f, "REBOUND IAS15 (inertial two-body)"


def _lsoda_fallback_cr3bp_rotating(
    state0: NDArray[np.float64],
    t_nondim: float,
    mu: float,
    *,
    mode: str,
) -> tuple[NDArray[np.float64], str]:
    """scipy LSODA fallback when rebound is unimportable.

    Integrates :func:`cr3bp.cr3bp_eom` directly under LSODA (multistep BDF —
    a distinct integrator family from the project's single-step DOP853). The
    ``mode`` is carried into the label for the audit trail even though the
    fallback always integrates the rotating-frame EOM directly.
    """
    from scipy.integrate import solve_ivp

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, float(t_nondim)),
        state0,
        args=(mu,),
        method="LSODA",
        rtol=1e-12,
        atol=1e-12,
    )
    if not sol.success:
        raise RuntimeError(f"LSODA fallback failed: {sol.message}") from None
    return sol.y[:, -1].copy(), f"scipy LSODA fallback ({mode})"


def run_v3_3d_periodic(
    candidate_id: str,
    state0: NDArray[np.float64],
    period_nondim: float,
    system: cr3bp.CR3BPSystem,
    *,
    v2_verdict: V2Verdict3D,
    n_cycles: int = V3_PERIODIC_N_CYCLES_MIN,
    closure_floor_nondim: float = V3_PERIODIC_CLOSURE_FLOOR_NONDIM,
    drift_agreement_floor_kms: float = V3_PERIODIC_DRIFT_AGREEMENT_FLOOR_KMS,
    prefer_rotating_callback: bool = True,
    degenerate_planar: bool = False,
    notes: str = "",
) -> V3PeriodicVerdict3D:
    """Run V3 for a CR3BP periodic orbit: REBOUND IAS15 independent cross-check.

    Pipeline:
      1. Propagate the IC one period under IAS15 in the rotating frame (mode A
         callback by default; mode B inertial two-body if
         ``prefer_rotating_callback=False``). Record the one-period closure
         ``||X(T) - X(0)||``.
      2. Propagate the IC ``n_cycles`` consecutive periods under IAS15,
         recording the cumulative position drift at each cycle boundary (the
         IAS15 analogue of V2's DOP853 drift series).
      3. Slice the V2 (DOP853) drift series to the cycles propagated and
         compute ``drift_agreement_kms = max_k |ias15[k] - dop853[k]|``.
      4. Verdict: PASS iff every cycle completed AND ``n >= 3`` AND the
         one-period IAS15 closure is below ``closure_floor_nondim`` AND the
         drift agreement is below ``drift_agreement_floor_kms``.

    Parameters
    ----------
    candidate_id:
        Identifier carried into the verdict.
    state0:
        6-vector IC ``(x, y, z, ẋ, ẏ, ż)``, nondim CR3BP rotating frame.
    period_nondim:
        Full nondim period (the V1/corrector-refined period).
    system:
        CR3BP system. Supplies ``mu`` (EOM) and the km / km-s conversions.
    v2_verdict:
        The :class:`V2Verdict3D` from :func:`run_v2_3d` on this candidate at
        the SAME (or larger) ``n_cycles``. V3 reads
        ``per_cycle_drift_kms`` as the DOP853 series to agree with.
    n_cycles:
        Cycles to propagate under IAS15. Must be ``>= 3`` (spec §14).
    closure_floor_nondim:
        Bar for the one-period IAS15 closure. Default 1e-7 nondim.
    drift_agreement_floor_kms:
        Bar for the IAS15-vs-DOP853 drift agreement. Default 1 km.
    prefer_rotating_callback:
        ``True`` → mode (A) rotating-frame callback (preferred, truly
        independent). ``False`` → mode (B) inertial two-body fallback.
    degenerate_planar:
        Carried from the candidate (e.g. from a V1 verdict). Diagnostic.
    notes:
        Free-form audit note.

    Returns
    -------
    V3PeriodicVerdict3D
        ``passes_v3`` is the headline.

    Notes
    -----
    A V3 PASS does NOT admit to the catalogue. V4 + V5 + the literature miss
    still gate. V3 asserts the V1/V2 verdict is integrator-architecture
    independent — nothing more.
    """
    state0_arr = np.asarray(state0, dtype=np.float64)
    if state0_arr.shape != (6,):
        raise ValueError(f"state0 must be 6D; got shape {state0_arr.shape}")
    if period_nondim <= 0.0:
        raise ValueError(f"period_nondim must be > 0; got {period_nondim}")
    if n_cycles < V3_PERIODIC_N_CYCLES_MIN:
        raise ValueError(
            f"V3 requires n_cycles >= {V3_PERIODIC_N_CYCLES_MIN} (spec §14); got {n_cycles}"
        )
    if closure_floor_nondim <= 0.0:
        raise ValueError(f"closure_floor_nondim must be > 0; got {closure_floor_nondim}")
    if drift_agreement_floor_kms <= 0.0:
        raise ValueError(f"drift_agreement_floor_kms must be > 0; got {drift_agreement_floor_kms}")
    if system.l_km <= 0.0 or system.t_s <= 0.0:
        raise ValueError(
            f"invalid CR3BP system for V3 km conversion: l_km={system.l_km} t_s={system.t_s}"
        )
    if len(v2_verdict.per_cycle_drift_kms) < n_cycles:
        raise ValueError(
            f"v2_verdict has only {len(v2_verdict.per_cycle_drift_kms)} cycles "
            f"but V3 wants {n_cycles}"
        )

    mu = float(system.mu)
    l_km = float(system.l_km)
    v_unit_kms = l_km / float(system.t_s)
    mode = _MODE_ROTATING_CALLBACK if prefer_rotating_callback else _MODE_INERTIAL_TWOBODY

    # --- (1) one-period closure under IAS15 ---
    state_one_period, integrator_label = _ias15_propagate_cr3bp_rotating(
        state0_arr, float(period_nondim), mu, mode=mode
    )
    closure_nondim = float(np.linalg.norm(state_one_period - state0_arr))
    closure_kms = closure_nondim * v_unit_kms

    # --- (2) multi-cycle drift under IAS15 (sequential, cycle-by-cycle) ---
    drifts_kms_ias15: list[float] = []
    current = state0_arr.copy()
    converged = True
    n_done = 0
    for _ in range(n_cycles):
        try:
            current, _label = _ias15_propagate_cr3bp_rotating(
                current, float(period_nondim), mu, mode=mode
            )
        except RuntimeError:
            converged = False
            break
        n_done += 1
        pos_delta_nondim = float(np.linalg.norm(current[:3] - state0_arr[:3]))
        drifts_kms_ias15.append(pos_delta_nondim * l_km)

    # --- (3) agreement vs the V2 (DOP853) series, sliced to cycles done ---
    dop853_series = tuple(float(d) for d in v2_verdict.per_cycle_drift_kms[:n_done])
    if n_done == 0:
        drift_agreement = float("inf")
    else:
        drift_agreement = max(
            abs(ias15 - dop853)
            for ias15, dop853 in zip(drifts_kms_ias15, dop853_series, strict=True)
        )

    converged_at_each_return = bool(converged and n_done == n_cycles)
    passes_v3 = bool(
        converged_at_each_return
        and n_done >= V3_PERIODIC_N_CYCLES_MIN
        and math.isfinite(closure_nondim)
        and closure_nondim <= closure_floor_nondim
        and math.isfinite(drift_agreement)
        and drift_agreement <= drift_agreement_floor_kms
    )

    return V3PeriodicVerdict3D(
        candidate_id=candidate_id,
        integrator=integrator_label,
        closure_residual_nondim_ias15=closure_nondim,
        closure_residual_kms_ias15=closure_kms,
        per_cycle_drift_kms_ias15=tuple(drifts_kms_ias15),
        per_cycle_drift_kms_dop853=dop853_series,
        drift_agreement_kms=float(drift_agreement),
        closure_floor_nondim=float(closure_floor_nondim),
        drift_agreement_floor_kms=float(drift_agreement_floor_kms),
        n_cycles_propagated=int(n_done),
        converged_at_each_return=converged_at_each_return,
        passes_v3=passes_v3,
        degenerate_planar=bool(degenerate_planar),
        notes=notes,
    )


__all__ = [
    "V3_PERIODIC_CLOSURE_FLOOR_NONDIM",
    "V3_PERIODIC_DRIFT_AGREEMENT_FLOOR_KMS",
    "V3_PERIODIC_N_CYCLES_MIN",
    "V3PeriodicVerdict3D",
    "run_v3_3d_periodic",
]
