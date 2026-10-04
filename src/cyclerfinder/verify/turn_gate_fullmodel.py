"""Demanded-turn gate applied to full-model passes in a rotating frame (#937).

The demanded-turn gate of :mod:`cyclerfinder.verify.turn_gate` compares the
incoming and outgoing V-infinity VECTORS of a flyby. For a patched-conic chain
those vectors are given. For an orbit of the restricted problem (CR3BP,
bicircular) they must be measured, and the measurement has one trap: the
state is in a ROTATING frame, so two body-relative velocities taken at the
entry and the exit of the sphere of influence, a day or more apart, are
expressed in axes that have themselves turned by ``n (t_out - t_in)`` in the
meantime (13.2 degrees per day in the Earth-Moon system). Differencing their
rotating-frame components adds that rotation to the turn.

That is what the ``#899`` review did (``docs/notes/2026-10-05-899-second-
species-scan-adversarial-review.md`` section 5 item 6): its "32.1 degrees
demanded against 19.1 available" for Casoliva's published 7-3b/c orbit is the
rotating-axes angle (32.12 reproduced here) against the turn of the
osculating hyperbola at the pass radius (19.11 reproduced); in inertial axes
the pass turns by 18.38 degrees, against 70.0 available at the 100 km floor.
All six of the review's turn numbers (orbits 3 and 7 and 7-3b/c) are
reproduced to 0.1 degree by the rotating-axes angle; see the ``#937`` note.

What this module measures (per pass of the secondary, model units then km/s)
---------------------------------------------------------------------------
* periapsis ``r_p`` (``rho . w = 0`` event), the Moon-relative inertial speed
  ``v_p`` there, and the OSCULATING two-body hyperbola at periapsis:
  ``e = r_p v_p^2 / mu - 1``, ``v_inf = sqrt(v_p^2 - 2 mu / r_p)``, turn
  ``2 asin(1/e)`` (Breakwell & Perko 1974 and Guillaume 1975 digests: exact for
  the two-body hyperbola);
* the window: entry and exit of the sphere ``|rho| = r_w`` around the
  periapsis (default ``r_w`` = Laplace sphere of influence
  ``(mu / (1 - mu))^(2/5)``), the secondary-relative inertial velocity vectors
  ``w_in`` and ``w_out`` there (rotating components turned into inertial axes
  by ``R(t)``), and the demanded turn ``angle(w_in, w_out)``;
* the decomposition of that turn: the rotation rate of ``w`` is
  ``(w x dw/dt) . h / |w|^2``; ``dw/dt`` (inertial) is split into the
  secondary's own pull ``-mu rho / |rho|^3`` and the rest (the primary's tide,
  the Sun in the bicircular model), and both are integrated across the window;
* the rotating-axes angle the ``#899`` review used, as a diagnostic only.

Verdict per pass (``status``)
-----------------------------
* ``"fail"``: the periapsis is below ``R + floor`` (inside the body when
  below ``R``). A full-model orbit is flyable at a pass only if it clears the
  floor; a published point-mass orbit that does not clear it fails here, which
  is correct (it is a mathematical orbit, not a flyable one).
* ``"indeterminate"`` (#906: never a rejection):
  - the pass is not complete inside the propagated arc;
  - the osculating orbit at periapsis is not a hyperbola (``e <= 1``: a slow,
    temporarily captured pass, the Henon / Brjuno low-energy regime, where no
    V-infinity and no patched-conic turn is defined);
  - a radial pass (zero osculating angular momentum, ``e = 1`` collision arcs);
  - the window turn exceeds the floor bend although the periapsis clears the
    floor: a two-body hyperbola above the floor cannot turn that much, so the
    excess is the non-secondary part, i.e. the patched conic does not describe
    this pass.
* ``"pass"``: otherwise, i.e. the periapsis clears the floor and the gate of
  :mod:`turn_gate` on ``(w_in, w_out)`` finds the turn within the floor bend.

The gate's own patched-conic verdict (``gate.turn_feasible``) is reported
unchanged for every pass.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.verify.turn_gate import Encounter, EncounterTurn, evaluate_encounter

Vec = NDArray[np.float64]
Rhs = Callable[[float, Vec], Vec]

_Z = np.array([0.0, 0.0, 1.0])

#: Below this periapsis (model length units) a pass is also propagated with the
#: KS-regularised integrator of :mod:`cyclerfinder.core.cr3bp_ks` (CR3BP only) and
#: the two window-exit states are compared (``ks_exit_mismatch``).
KS_CHECK_RP: float = 1.0e-3


@dataclass(frozen=True)
class RotatingModel:
    """A model in a frame rotating uniformly at rate 1 about +z, with the flyby
    body (mass parameter ``mu``) fixed at ``(1 - mu, 0, 0)``.

    ``rhs(t, state6)`` returns the rotating-frame derivative. ``l_km`` and
    ``t_s`` convert to km and seconds. ``radius_km`` and ``alt_floor_km`` are
    the flyby body's radius and altitude floor. ``ks_mu`` is set for the plain
    CR3BP, enabling the KS cross-check of deep passes.
    """

    name: str
    mu: float
    rhs: Rhs
    l_km: float
    t_s: float
    body: str
    radius_km: float
    alt_floor_km: float
    ks_mu: float | None = None

    @property
    def v_kms(self) -> float:
        return self.l_km / self.t_s

    @property
    def gm_km3_s2(self) -> float:
        """The body's GM implied by the model's own ``mu`` and scales."""
        return self.mu * self.l_km**3 / self.t_s**2

    @property
    def soi(self) -> float:
        """Laplace sphere of influence ``(mu / (1 - mu))^(2/5)`` (model units)."""
        return float((self.mu / (1.0 - self.mu)) ** 0.4)

    @property
    def tidal_turn_scale_speed(self) -> float:
        """``n r_soi`` in model units (n = 1): the speed scale of the tidal turn
        ``(n r_soi / v)^2`` across the sphere of influence."""
        return self.soi


def cr3bp_model(
    mu: float,
    *,
    l_km: float,
    t_s: float,
    body: str,
    radius_km: float,
    alt_floor_km: float,
) -> RotatingModel:
    """The planar/spatial CR3BP of :func:`cyclerfinder.core.cr3bp.cr3bp_eom`."""
    from cyclerfinder.core.cr3bp import cr3bp_eom

    def rhs(t: float, s: Vec) -> Vec:
        return cr3bp_eom(t, s, mu)

    return RotatingModel("cr3bp", mu, rhs, l_km, t_s, body, radius_km, alt_floor_km, ks_mu=mu)


def bcr4bp_model(
    system: object,
    *,
    l_km: float,
    t_s: float,
    body: str,
    radius_km: float,
    alt_floor_km: float,
) -> RotatingModel:
    """The bicircular model of :func:`cyclerfinder.core.bcr4bp.bcr4bp_eom` (its
    synodic frame rotates uniformly with the Earth-Moon line, the Moon fixed at
    ``1 - mu``). The quasi-bicircular model is NOT supported: its frame does not
    rotate uniformly, so ``R(t)`` below would be wrong there."""
    from cyclerfinder.core.bcr4bp import BCR4BPSystem, bcr4bp_eom

    assert isinstance(system, BCR4BPSystem)

    def rhs(t: float, s: Vec) -> Vec:
        return bcr4bp_eom(t, s, system)

    return RotatingModel("bcr4bp", float(system.mu), rhs, l_km, t_s, body, radius_km, alt_floor_km)


def _rot(t: float) -> Vec:
    c, s = math.cos(t), math.sin(t)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def body_relative_inertial(model: RotatingModel, t: float, state: Vec) -> tuple[Vec, Vec]:
    """Body-relative position and inertial velocity, both in INERTIAL axes
    (axes aligned with the rotating ones at ``t = 0``)."""
    rho = np.asarray(state[:3], dtype=np.float64) - np.array([1.0 - model.mu, 0.0, 0.0])
    w = np.asarray(state[3:6], dtype=np.float64) + np.cross(_Z, rho)
    r = _rot(t)
    return r @ rho, r @ w


def rotating_axes_velocity(model: RotatingModel, state: Vec) -> Vec:
    """Body-relative inertial velocity in the ROTATING axes of the same instant.
    Differencing two of these taken at different times is the ``#899`` error."""
    rho = np.asarray(state[:3], dtype=np.float64) - np.array([1.0 - model.mu, 0.0, 0.0])
    return np.asarray(state[3:6], dtype=np.float64) + np.cross(_Z, rho)


def _accel_split(model: RotatingModel, t: float, state: Vec) -> tuple[Vec, Vec, Vec]:
    """(w, body part, rest) of the inertial derivative of ``w``, rotating axes."""
    d = model.rhs(t, state)
    rho = state[:3] - np.array([1.0 - model.mu, 0.0, 0.0])
    v = state[3:6]
    w = v + np.cross(_Z, rho)
    total = d[3:6] + 2.0 * np.cross(_Z, v) + np.cross(_Z, np.cross(_Z, rho))
    body = -model.mu * rho / float(np.linalg.norm(rho)) ** 3
    return w, body, total - body


def _angle(a: Vec, b: Vec) -> float:
    c = float(a @ b) / (float(np.linalg.norm(a)) * float(np.linalg.norm(b)))
    return math.acos(max(-1.0, min(1.0, c)))


@dataclass(frozen=True)
class PassTurn:
    """One pass of the flyby body (angles in degrees, speeds in km/s)."""

    t_periapsis: float
    rp_km: float
    alt_km: float
    window_radius_km: float
    t_in: float
    t_out: float
    window_days: float
    complete: bool
    n_periapses: int
    vp_kms: float
    e_osc: float
    vinf_osc_kms: float
    turn_osc_deg: float
    v_window_in_kms: float
    v_window_out_kms: float
    demanded_turn_deg: float
    body_part_deg: float
    rest_part_deg: float
    rotating_axes_turn_deg: float
    hyperbola_turn_at_rp_deg: float
    tidal_turn_scale_deg: float
    ks_exit_mismatch: float
    gate: EncounterTurn | None
    status: str
    reason: str

    def as_dict(self) -> dict[str, object]:
        d = {k: getattr(self, k) for k in self.__dataclass_fields__ if k != "gate"}
        d["gate"] = None if self.gate is None else self.gate.as_dict()
        return d


@dataclass(frozen=True)
class FullModelTurnReport:
    """All passes of one trajectory inside the window radius."""

    passes: tuple[PassTurn, ...]
    min_distance_km: float
    notes: tuple[str, ...] = field(default=())

    @property
    def status(self) -> str:
        """``"no_encounter"`` when no pass enters the window; else ``"fail"`` if any
        pass fails, ``"indeterminate"`` if any is indeterminate, else ``"pass"``."""
        if not self.passes:
            return "no_encounter"
        st = {p.status for p in self.passes}
        if "fail" in st:
            return "fail"
        if "indeterminate" in st:
            return "indeterminate"
        return "pass"

    @property
    def worst_ratio(self) -> float:
        return max((p.gate.ratio for p in self.passes if p.gate is not None), default=0.0)


def full_model_pass_turns(
    model: RotatingModel,
    state0: Vec,
    period: float,
    *,
    periodic: bool = True,
    window_radius: float | None = None,
    rtol: float = 1e-12,
    atol: float = 1e-14,
    decompose: bool = True,
) -> FullModelTurnReport:
    """Measure every pass of the flyby body on a trajectory and gate it.

    ``periodic=True``: ``state0`` is on a periodic orbit of period ``period``;
    the orbit is propagated over two periods and the passes whose periapsis
    falls in ``[T/2, 3T/2)`` are kept, so a pass straddling ``t = 0`` is seen
    whole, once. ``periodic=False``: the arc ``[0, period]`` only; a pass whose
    window is cut by an end of the arc is ``"indeterminate"``.
    """
    s0 = np.asarray(state0, dtype=np.float64).reshape(-1)
    if s0.size == 4:
        s0 = np.array([s0[0], s0[1], 0.0, s0[2], s0[3], 0.0])
    r_w = model.soi if window_radius is None else float(window_radius)
    c = np.array([1.0 - model.mu, 0.0, 0.0])
    # periodic: one period each way from state0, passes with periapsis in [-T/2, T/2) (each
    # pass once, seen whole, with the error growth of an unstable orbit limited to one period)
    spans = [(0.0, period), (0.0, -period)] if periodic else [(0.0, period)]
    t_lo, t_hi = (-period, period) if periodic else (0.0, period)

    def ev_peri(t: float, s: Vec) -> float:
        rho = s[:3] - c
        return float(rho @ s[3:6])

    def ev_window(t: float, s: Vec) -> float:
        return float(np.linalg.norm(s[:3] - c)) - r_w

    dense_parts: list[Callable[[float], Vec]] = []
    t_peri_l: list[float] = []
    t_win_l: list[float] = []
    min_dist = math.inf
    for span in spans:
        sol = solve_ivp(
            lambda t, s: model.rhs(t, s),
            span,
            s0,
            method="DOP853",
            rtol=rtol,
            atol=atol,
            dense_output=True,
            events=[ev_peri, ev_window],
        )
        if sol.status != 0 or sol.sol is None or sol.t_events is None:
            raise RuntimeError(f"full_model_pass_turns: propagation failed ({sol.message})")
        sol_dense = sol.sol
        dense_parts.append(sol_dense)
        t_win_l.extend(float(t) for t in sol.t_events[1])
        for t in sol.t_events[0]:
            st = sol_dense(t)
            # keep minima of |rho| only: d2(|rho|^2 / 2)/dt2 = |v|^2 + rho . a > 0
            acc = model.rhs(float(t), st)[3:6]
            if float(st[3:6] @ st[3:6] + (st[:3] - c) @ acc) > 0.0:
                t_peri_l.append(float(t))
        in_first_period = (np.abs(sol.t) <= period) if periodic else np.ones(sol.t.size, bool)
        rho_steps = np.linalg.norm(np.asarray(sol.y)[:3, in_first_period] - c[:, None], axis=0)
        if rho_steps.size:
            min_dist = min(min_dist, float(rho_steps.min()))

    def dense(t: float) -> Vec:
        part = dense_parts[0] if (t >= 0.0 or len(dense_parts) == 1) else dense_parts[1]
        return np.asarray(part(t), dtype=np.float64)

    t_peri = sorted({round(t, 12) for t in t_peri_l})
    t_win = np.array(sorted({round(t, 12) for t in t_win_l}))
    for tp in t_peri:
        if periodic and -0.5 * period <= tp < 0.5 * period:
            min_dist = min(min_dist, float(np.linalg.norm(dense(tp)[:3] - c)))

    # one pass per window: a window holding several periapses is a captured pass, judged at
    # its deepest periapsis and counted
    groups: dict[tuple[float, float, bool], tuple[float, float, int]] = {}
    for tp in t_peri:
        if periodic and not (-0.5 * period <= tp < 0.5 * period):
            continue
        rp = float(np.linalg.norm(dense(tp)[:3] - c))
        if rp >= r_w:
            continue
        before = t_win[t_win < tp]
        after = t_win[t_win > tp]
        complete = before.size > 0 and after.size > 0
        t_in = float(before[-1]) if before.size else t_lo
        t_out = float(after[0]) if after.size else t_hi
        key = (t_in, t_out, complete)
        if key in groups:
            tp0, rp0, n0 = groups[key]
            groups[key] = (tp, rp, n0 + 1) if rp < rp0 else (tp0, rp0, n0 + 1)
        else:
            groups[key] = (tp, rp, 1)
    passes = [
        _measure_pass(model, dense, tp, rp, t_in, t_out, r_w, complete, n, rtol, atol, decompose)
        for (t_in, t_out, complete), (tp, rp, n) in sorted(groups.items(), key=lambda kv: kv[1][0])
    ]
    min_km = min_dist * model.l_km
    return FullModelTurnReport(tuple(passes), min_km)


def _measure_pass(
    model: RotatingModel,
    dense: Callable[[float], Vec],
    tp: float,
    rp: float,
    t_in: float,
    t_out: float,
    r_w: float,
    complete: bool,
    n_periapses: int,
    rtol: float,
    atol: float,
    decompose: bool,
) -> PassTurn:
    mu = model.mu
    vk = model.v_kms
    sp = dense(tp)
    s_in = dense(t_in)
    s_out = dense(t_out)
    rho_p, w_p = body_relative_inertial(model, tp, sp)
    vp = float(np.linalg.norm(w_p))
    h = np.cross(rho_p, w_p)
    h_norm = float(np.linalg.norm(h))
    e_osc = rp * vp * vp / mu - 1.0
    energy = 0.5 * vp * vp - mu / rp
    vinf_osc = math.sqrt(2.0 * energy) if energy > 0.0 else math.nan
    turn_osc = math.degrees(2.0 * math.asin(1.0 / e_osc)) if e_osc > 1.0 else math.nan

    _, w_in = body_relative_inertial(model, t_in, s_in)
    _, w_out = body_relative_inertial(model, t_out, s_out)
    demanded = math.degrees(_angle(w_in, w_out))
    rot_axes = math.degrees(
        _angle(rotating_axes_velocity(model, s_in), rotating_axes_velocity(model, s_out))
    )
    v_in = float(np.linalg.norm(w_in))
    v_out = float(np.linalg.norm(w_out))
    v_min = min(v_in, v_out)
    e_rp = 1.0 + rp * v_min * v_min / mu
    hyp_rp = math.degrees(2.0 * math.asin(1.0 / e_rp))
    tidal = math.degrees((model.tidal_turn_scale_speed / v_min) ** 2) if v_min > 0 else math.inf

    body_part = rest_part = math.nan
    if decompose and complete:
        n_hat_inertial = h / h_norm if h_norm > 0 else _Z

        def aug(t: float, y: Vec) -> Vec:
            s = y[:6]
            w, a_b, a_r = _accel_split(model, t, s)
            n_rot = _rot(t).T @ n_hat_inertial
            w2 = float(w @ w)
            return np.concatenate(
                [
                    model.rhs(t, s),
                    [float(np.cross(w, a_b) @ n_rot) / w2, float(np.cross(w, a_r) @ n_rot) / w2],
                ]
            )

        a = solve_ivp(
            aug,
            (t_in, t_out),
            np.concatenate([s_in, [0.0, 0.0]]),
            method="DOP853",
            rtol=rtol,
            atol=atol,
        )
        body_part = math.degrees(float(a.y[6, -1]))
        rest_part = math.degrees(float(a.y[7, -1]))

    ks_mismatch = math.nan
    if model.ks_mu is not None and complete and rp < KS_CHECK_RP:
        from cyclerfinder.core.cr3bp_ks import MoonCentredCR3BP, propagate_ks

        arc = propagate_ks(MoonCentredCR3BP(model.ks_mu), s_in, t_out, t0=t_in, rtol=1e-13)
        ks_mismatch = float(np.max(np.abs(arc.state - s_out)))

    enc = Encounter(
        body=model.body,
        mu_km3_s2=model.gm_km3_s2,
        radius_km=model.radius_km,
        alt_floor_km=model.alt_floor_km,
        vinf_in=w_in * vk,
        vinf_out=w_out * vk,
        label=f"t={tp:.6g}",
        tidal_speed_kms=model.tidal_turn_scale_speed * vk,
    )
    gate = evaluate_encounter(enc) if v_min > 0 else None

    rp_km = rp * model.l_km
    floor_km = model.radius_km + model.alt_floor_km
    if rp_km < floor_km:
        status = "fail"
        where = "inside the body" if rp_km < model.radius_km else "below the altitude floor"
        reason = f"periapsis {rp_km - model.radius_km:.1f} km altitude: {where}"
    elif not complete:
        status, reason = "indeterminate", "pass not complete inside the propagated arc"
    elif n_periapses > 1:
        status = "indeterminate"
        reason = f"{n_periapses} periapses inside one window: a temporarily captured pass"
    elif h_norm <= 1e-12 * rp * vp:
        status, reason = "indeterminate", "radial pass: osculating angular momentum zero"
    elif e_osc <= 1.0:
        status = "indeterminate"
        reason = f"osculating orbit at periapsis is not a hyperbola (e = {e_osc:.4f}): slow pass"
    elif gate is None or not gate.turn_feasible:
        status = "indeterminate"
        reason = (
            "window turn exceeds the floor bend although the periapsis clears the floor: "
            "the non-body part of the turn is outside the patched-conic model"
        )
    else:
        status, reason = "pass", ""
    return PassTurn(
        t_periapsis=tp,
        rp_km=rp_km,
        alt_km=rp_km - model.radius_km,
        window_radius_km=r_w * model.l_km,
        t_in=t_in,
        t_out=t_out,
        window_days=(t_out - t_in) * model.t_s / 86400.0,
        complete=complete,
        n_periapses=n_periapses,
        vp_kms=vp * vk,
        e_osc=e_osc,
        vinf_osc_kms=vinf_osc * vk,
        turn_osc_deg=turn_osc,
        v_window_in_kms=v_in * vk,
        v_window_out_kms=v_out * vk,
        demanded_turn_deg=demanded,
        body_part_deg=body_part,
        rest_part_deg=rest_part,
        rotating_axes_turn_deg=rot_axes,
        hyperbola_turn_at_rp_deg=hyp_rp,
        tidal_turn_scale_deg=tidal,
        ks_exit_mismatch=ks_mismatch,
        gate=gate,
        status=status,
        reason=reason,
    )
