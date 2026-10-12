"""Varin's closed families i_k of symmetric periodic orbits in the planar CR3BP (#975).

Source: Bruno & Varin 2009, "Closed Families of Periodic Solutions of a Restricted Three-Body
Problem", Solar System Research 43(3):253-276 (the journal form of Varin's KIAM preprint 16/2008);
printed values in data/sources/varin-bruno-2009-closed-families-tables.yaml. Pre-registration:
docs/notes/2026-10-11-975-varin-cascade-satellite-mass-ratios.md.

Varin's frame: synodic, origin at P1 (mass 1 - mu), P2 at (1, 0), canonical momenta
y1 = dx1/dt - x2, y2 = dx2/dt + x1. The project's frame is barycentric (P1 at (-mu, 0)), so
x = x1 - mu, y = x2, u = y1 + x2, v = y2 - x1 (rotating-frame velocities are unchanged by the
shift of origin). With H = (u^2 + v^2)/2 - (x1^2 + x2^2)/2 - (1-mu)/r1 - mu/r2 + mu x1, the
barycentric Jacobi constant is C_J = -2H + mu^2 = C_printed - mu (1 - mu) (Varin prints -2H + mu).

A closed family i_k branches off family i at a saddle point mu'_k of the characteristic surface and
shrinks to one orbit at an elliptic point mu''_k. Both are solutions of the 4 x 4 system
x2(T/2) = 0, y1(T/2) = 0, m_41 = 0, m_14 = 0 in (x1(0), y2(0), T, mu) (Bruno & Varin 2009
eq. (2)-(3)); m_ij are elements of the 4 x 4 monodromy matrix in Varin's coordinates.

Everything is in the CR3BP's nondimensional units at the system's mu; the caller converts lengths to
km with the system's l* (`CR3BPSystem.l_km`).
"""

from __future__ import annotations

import math
import time
from collections.abc import Callable
from dataclasses import dataclass, replace

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core import cr3bp

Arr = NDArray[np.float64]

#: Varin coordinates z = (x1, x2, y1, y2) -> project rotating coordinates w = (x, y, u, v)
#: (without the constant shift x = x1 - mu): w = A z.
_A = np.array(
    [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 1.0, 1.0, 0.0], [-1.0, 0.0, 0.0, 1.0]]
)
_A_INV = np.linalg.inv(_A)
_PLANAR = [0, 1, 3, 4]
RTOL = ATOL = 1e-12


def system_at(mu: float, base: cr3bp.CR3BPSystem | None = None) -> cr3bp.CR3BPSystem:
    """A CR3BP system with mass ratio ``mu`` (the base system's l* and t*, or unit scales)."""
    if base is None:
        return cr3bp.CR3BPSystem(mu=mu, primary="P1", secondary="P2", l_km=1.0, t_s=1.0)
    return replace(base, mu=mu)


def crossing_from_varin(a: float, e: float, mu: float) -> tuple[float, float]:
    """Varin's astronomical coordinates (a~, e~) of a symmetric crossing -> (x1(0), y2(0))."""
    x1 = a * abs(2.0 - abs(e))
    y2 = math.copysign(1.0, a * e) * math.sqrt(abs(e) * (1.0 - mu) / abs(x1))
    return x1, y2


def state_from_crossing(x1: float, y2: float, mu: float) -> Arr:
    """Barycentric 6-state of the symmetric crossing x2 = y1 = 0 at (x1, y2)."""
    return np.array([x1 - mu, 0.0, 0.0, 0.0, y2 - x1, 0.0])


def jacobi_from_printed(c_printed: float, mu: float) -> float:
    """Varin's printed C (-2H + mu) -> the barycentric Jacobi constant."""
    return c_printed - mu * (1.0 - mu)


def _prop(state: Arr, t: float, mu: float, with_stm: bool) -> cr3bp.CR3BPArc:
    return cr3bp.propagate(
        system_at(mu), state, t, with_stm=with_stm, rtol=RTOL, atol=ATOL, stm_mode="fixed_path"
    )


def half_residuals(x1: float, y2: float, t_full: float, mu: float) -> tuple[float, float]:
    """(x2, y1) at T/2 in Varin's coordinates; both vanish for a symmetric periodic orbit."""
    sf = _prop(state_from_crossing(x1, y2, mu), 0.5 * t_full, mu, False).state_f
    return float(sf[1]), float(sf[3] - sf[1])


def monodromy_varin(x1: float, y2: float, t_full: float, mu: float) -> tuple[Arr, Arr]:
    """(4 x 4 monodromy in Varin's coordinates, 2 x 2 vertical block) over the full period."""
    arc = _prop(state_from_crossing(x1, y2, mu), t_full, mu, True)
    assert arc.stm is not None
    m6 = np.asarray(arc.stm)
    m4 = m6[np.ix_(_PLANAR, _PLANAR)]
    mz = m6[np.ix_([2, 5], [2, 5])]
    return _A_INV @ m4 @ _A, mz


def traces(x1: float, y2: float, t_full: float, mu: float) -> tuple[float, float]:
    """Varin's Tr = trace(M) - 2 (plane) and Tr_v (vertical)."""
    m, mz = monodromy_varin(x1, y2, t_full, mu)
    return float(np.trace(m) - 2.0), float(np.trace(mz))


def modified_trace(tr: float) -> float:
    """Varin's Tr~: Tr if |Tr| <= 2, else (1 + log2|Tr|) sgn Tr."""
    return tr if abs(tr) <= 2.0 else math.copysign(1.0 + math.log2(abs(tr)), tr)


def crossing_near_half(x1: float, y2: float, t_full: float, mu: float) -> tuple[float, float]:
    """Integrate the crossing state and return (time, |y1|) at the x2 = 0 crossing nearest T/2."""

    def ev(t: float, s: Arr, mu_: float) -> float:
        return float(s[1])

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, 0.75 * t_full),
        state_from_crossing(x1, y2, mu),
        args=(mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=RTOL,
        atol=ATOL,
        events=ev,
    )
    ts = sol.t_events[0]
    ys = sol.y_events[0]
    i = int(np.argmin(np.abs(ts - 0.5 * t_full)))
    return float(ts[i]), float(abs(ys[i][3] - ys[i][1]))


def critical_residuals(z: Arr) -> Arr:
    """Varin's 4 x 4 system at z = (x1(0), y2(0), T, mu)."""
    x1, y2, t_full, mu = (float(v) for v in z)
    r1, r2 = half_residuals(x1, y2, t_full, mu)
    m, _ = monodromy_varin(x1, y2, t_full, mu)
    return np.array([r1, r2, m[3, 0], m[0, 3]])


def critical_residuals_half(z: Arr) -> Arr:
    """The same critical points from the half-period STM only: (x2, y1, D_y, D_x), with D_y, D_x the
    determinants of the half-residual Jacobian's (x1, T) and (y2, T) columns (extrema of y2(0) and
    of x1(0) along the family at fixed mu; Bruno & Varin's m_41 = 0 and m_14 = 0)."""
    x1, y2, t_full, mu = (float(v) for v in z)
    r, jac = half_residuals_jac(x1, y2, t_full, mu)
    return np.array(
        [r[0], r[1], float(np.linalg.det(jac[:, [0, 2]])), float(np.linalg.det(jac[:, [1, 2]]))]
    )


@dataclass(frozen=True)
class CriticalPoint:
    x1: float
    y2: float
    t_full: float
    mu: float
    converged: bool
    residual: float
    iterations: int
    tr: float
    tr_v: float


def solve_critical_point(
    z0: Arr,
    *,
    tol: float = 1e-11,
    max_it: int = 30,
    rel_step: float = 1e-7,
    half: bool = False,
) -> CriticalPoint:
    """Newton on Varin's 4 x 4 system with a central finite-difference Jacobian (steps relative to
    each unknown's magnitude, so the step scales with mu at any mass ratio). ``half`` uses
    :func:`critical_residuals_half` (half-period STM; determinant rows held to 100 x tol)."""
    fun = critical_residuals_half if half else critical_residuals
    tols = np.array([tol, tol, DET_TOL_FACTOR * tol, DET_TOL_FACTOR * tol]) if half else tol
    z = np.asarray(z0, dtype=float).copy()
    f = fun(z)
    it = 0
    while it < max_it and not bool(np.all(np.abs(f) < tols)):
        it += 1
        jac = np.empty((4, 4))
        for j in range(4):
            h = rel_step * max(abs(z[j]), 1e-300)
            zp, zm = z.copy(), z.copy()
            zp[j] += h
            zm[j] -= h
            jac[:, j] = (fun(zp) - fun(zm)) / (2.0 * h)
        dz = np.linalg.solve(jac, -f)
        # keep mu positive and the step inside a trust region of 10% of each unknown
        scale = float(np.max(np.abs(dz) / (0.1 * np.abs(z))))
        if scale > 1.0:
            dz /= scale
        z = z + dz
        f = fun(z)
    res = float(np.max(np.abs(f)))
    conv = bool(np.all(np.abs(f) < tols))
    tr, tr_v = traces(z[0], z[1], z[2], z[3])
    return CriticalPoint(
        x1=float(z[0]),
        y2=float(z[1]),
        t_full=float(z[2]),
        mu=float(z[3]),
        converged=conv,
        residual=res,
        iterations=it,
        tr=tr,
        tr_v=tr_v,
    )


# ---------------------------------------------------------------------------
# Continuation: the y2-extremum curve (x2 = y1 = m_41 = 0) in mu, and a closed family at fixed mu
# ---------------------------------------------------------------------------

_DR_DW = np.array([[0.0, 1.0, 0.0, 0.0], [0.0, -1.0, 1.0, 0.0]])  # (x2, y1) from (x, y, u, v)
_DW0 = np.array([[1.0, 0.0], [0.0, 0.0], [0.0, 0.0], [-1.0, 1.0]])  # d w0 / d (x1, y2)


def half_residuals_jac(x1: float, y2: float, t_full: float, mu: float) -> tuple[Arr, Arr]:
    """(x2, y1) at T/2 and their 2 x 3 Jacobian in (x1, y2, T), from the half-period STM."""
    arc = _prop(state_from_crossing(x1, y2, mu), 0.5 * t_full, mu, True)
    assert arc.stm is not None
    sf = np.asarray(arc.state_f)
    phi4 = np.asarray(arc.stm)[np.ix_(_PLANAR, _PLANAR)]
    r = np.array([sf[1], sf[3] - sf[1]])
    f = cr3bp.cr3bp_eom(0.5 * t_full, sf, mu)[[0, 1, 3, 4]]
    jac = np.column_stack([_DR_DW @ phi4 @ _DW0, 0.5 * (_DR_DW @ f)])
    return r, jac


def extremum_residuals(w: Arr) -> Arr:
    """(x2(T/2), y1(T/2), D_y) at w = (x1(0), y2(0), T, log mu).

    D_y = det of the half-residual Jacobian's (x1, T) columns: it vanishes where y2(0) is extremal
    along the family at fixed mu, the condition Bruno & Varin write as m_41 = 0 (checked on the
    eight Table 1 orbits: D_y ~ 1e-7 there against O(1) nearby). It needs only the half-period
    STM, not the monodromy."""
    x1, y2, t_full, lmu = (float(v) for v in w)
    r, jac = half_residuals_jac(x1, y2, t_full, math.exp(lmu))
    return np.array([r[0], r[1], float(np.linalg.det(jac[:, [0, 2]]))])


def correct_extremum_at_mu(
    w_guess: Arr, mu: float, *, tol: float = 1e-10, max_it: int = 20, rel_step: float = 1e-7
) -> Arr | None:
    """Newton on x2 = y1 = m_41 = 0 in (x1, y2, T) at fixed mu; returns w = (x1, y2, T, log mu)."""
    lmu = math.log(mu)

    def fun(v: Arr) -> Arr:
        return extremum_residuals(np.array([v[0], v[1], v[2], lmu]))

    v = np.asarray(w_guess[:3], dtype=float).copy()
    for _ in range(max_it):
        f = fun(v)
        if _ext_ok(f, tol):
            return np.array([v[0], v[1], v[2], lmu])
        steps = rel_step * np.abs(v)
        v = v + np.linalg.solve(_fd_jac(fun, v, steps), -f)
    return None


def _fd_jac(fun: object, w: Arr, steps: Arr) -> Arr:
    f0 = fun(w)  # type: ignore[operator]
    cols = []
    for j in range(len(w)):
        wp = w.copy()
        wp[j] += steps[j]
        cols.append((fun(wp) - f0) / steps[j])  # type: ignore[operator]
    return np.column_stack(cols)


def _tangent(jac: Arr, prev: Arr | None) -> Arr:
    _, _, vt = np.linalg.svd(jac)
    t = vt[-1]
    if prev is not None and float(t @ prev) < 0.0:
        t = -t
    return np.asarray(t / np.linalg.norm(t))


@dataclass
class CurvePoint:
    w: list[float]
    mu: float
    residual: float


_SCALE = np.array([1.0, 1.0, 2.0 * math.pi, 1.0])  # q = (x1, y2, T/2pi, log mu) -> w = q * _SCALE
#: tolerance on D_y relative to the orbit tolerance: D_y is a determinant of O(1-10) Jacobian
#: entries, so its integration noise is about 100 times the state's
DET_TOL_FACTOR = 100.0

#: Tool amendment A2 (#975 note sec. 11.4): a predictor-corrector step is accepted only if the
#: corrected point lies within MAX_DISP_FRAC |step| of the predictor and the tangent turns by at
#: most MAX_TURN_DEG; otherwise the step is halved. Without it a 0.02 step at a hairpin of the i_1
#: oval at the Earth-Moon mu landed 1.02 away on another family (the 3-fold cover of a T~ 1 orbit).
MAX_DISP_FRAC = 0.5
MAX_TURN_DEG = 10.0
_COS_MAX_TURN = math.cos(math.radians(MAX_TURN_DEG))


def _step_ok(corrected: Arr, predicted: Arr, step: float, t_old: Arr, t_new: Arr) -> str | None:
    """None if the step passes the A2 acceptance test, else the reason ("displacement", "turn")."""
    if float(np.linalg.norm(corrected - predicted)) > MAX_DISP_FRAC * abs(step):
        return "displacement"
    if float(t_old @ t_new) < _COS_MAX_TURN:
        return "turn"
    return None


def _ext_ok(f: Arr, tol: float) -> bool:
    return bool(max(abs(f[0]), abs(f[1])) < tol and abs(f[2]) < DET_TOL_FACTOR * tol)


def _ext_q(q: Arr) -> Arr:
    return extremum_residuals(q * _SCALE)


def follow_extremum_curve(
    w0: Arr,
    *,
    direction: Arr,
    mu_stop: float,
    ds: float = 2e-2,
    ds_min: float = 1e-5,
    ds_max: float = 2e-1,
    max_steps: int = 400,
    tol: float = 1e-10,
    rel_step: float = 1e-7,
    on_point: Callable[[CurvePoint], None] | None = None,
    deadline: float | None = None,
    guard: bool = True,
    stats: dict[str, int] | None = None,
) -> tuple[list[CurvePoint], str]:
    """Pseudo-arclength continuation of x2 = y1 = m_41 = 0 in w = (x1, y2, T, log mu).

    The arclength is measured in q = (x1, y2, T/2pi, log mu); the corrector is a chord Newton with
    the Jacobian of the last converged point (a fresh Jacobian on a retry). Stops when mu crosses
    ``mu_stop`` (the last point is then corrected at mu = mu_stop), at a fold in mu (the mu
    component of the tangent changes sign; returned as "fold"), or after max_steps. ``on_point`` is
    called with every new point (for checkpoints). ``deadline`` (time.time()) stops it with
    "time budget". ``guard`` applies the A2 step acceptance (:func:`_step_ok`); ``stats`` counts the
    rejections by reason."""
    q = np.asarray(w0, dtype=float) / _SCALE

    def steps_at(v: Arr) -> Arr:
        return np.array(
            [rel_step * abs(v[0]), rel_step * abs(v[1]), rel_step * abs(v[2]), rel_step]
        )

    jac = _fd_jac(_ext_q, q, steps_at(q))
    t = _tangent(jac, np.asarray(direction, dtype=float) / _SCALE)
    pts = [CurvePoint((q * _SCALE).tolist(), math.exp(q[3]), float(np.max(np.abs(_ext_q(q)))))]
    lstop = math.log(mu_stop)
    side = math.copysign(1.0, lstop - q[3])
    for _ in range(max_steps):
        if deadline is not None and time.time() > deadline:
            return pts, "time budget"
        ok = False
        qc = q
        fresh = False
        jac_new: Arr | None = None
        while ds >= ds_min:
            qp = q + ds * t
            qc = qp.copy()
            ja = np.vstack([jac, t])
            for _ in range(10):
                f = _ext_q(qc)
                if _ext_ok(f, tol):
                    ok = True
                    break
                g = np.concatenate([f, [float(t @ (qc - qp))]])
                qc = qc + np.linalg.solve(ja, -g)
            if ok and guard:
                jac_new = _fd_jac(_ext_q, qc, steps_at(qc))
                why_not = _step_ok(qc, qp, ds, t, _tangent(jac_new, t))
                if why_not is not None:
                    ok = False
                    if stats is not None:
                        stats[why_not] = stats.get(why_not, 0) + 1
                    ds *= 0.5
                    continue
            if ok:
                break
            if not fresh:
                jac = _fd_jac(_ext_q, q, steps_at(q))
                fresh = True
            else:
                ds *= 0.5
        if not ok:
            return pts, "corrector failed"
        jac = jac_new if jac_new is not None else _fd_jac(_ext_q, qc, steps_at(qc))
        t_new = _tangent(jac, t)
        if t_new[3] * t[3] < 0.0 and t[3] * side > 0.0:
            pts.append(CurvePoint((qc * _SCALE).tolist(), math.exp(qc[3]), 0.0))
            if on_point is not None:
                on_point(pts[-1])
            return pts, "fold"
        if (qc[3] - lstop) * side >= 0.0:
            qi = q + (lstop - q[3]) / (qc[3] - q[3]) * (qc - q)
            wf = correct_extremum_at_mu(qi * _SCALE, mu_stop)
            if wf is None:
                return pts, "correction at mu_stop failed"
            pts.append(CurvePoint(wf.tolist(), mu_stop, 0.0))
            if on_point is not None:
                on_point(pts[-1])
            return pts, "reached mu_stop"
        q, t = qc, t_new
        pts.append(CurvePoint((q * _SCALE).tolist(), math.exp(q[3]), 0.0))
        if on_point is not None:
            on_point(pts[-1])
        ds = min(ds * 1.5, ds_max)
    return pts, "max steps"


@dataclass
class LoopMember:
    x1: float
    y2: float
    t_full: float
    residual: float


def _loop_eval(p: Arr, mu: float) -> tuple[Arr, Arr]:
    """Half residuals and their 2 x 3 Jacobian in p = (x1, y2, T/2pi)."""
    r, jac = half_residuals_jac(float(p[0]), float(p[1]), 2.0 * math.pi * float(p[2]), mu)
    jac = jac.copy()
    jac[:, 2] *= 2.0 * math.pi
    return r, jac


def _loop_step(p: Arr, t: Arr, step: float, mu: float, tol: float) -> Arr | None:
    """One pseudo-arclength predictor-corrector step; None if the corrector does not converge."""
    pp = p + step * t
    pc = pp.copy()
    for _ in range(10):
        if pc[2] <= 0.0:
            return None
        try:
            r, jac = _loop_eval(pc, mu)
        except (RuntimeError, ValueError):
            return None
        if float(np.max(np.abs(r))) < tol:
            return pc
        g = np.concatenate([r, [float(t @ (pc - pp))]])
        pc = pc + np.linalg.solve(np.vstack([jac, t]), -g)
    return None


def _accepted_step(
    p: Arr, t: Arr, step: float, mu: float, tol: float, guard: bool, stats: dict[str, int] | None
) -> tuple[Arr, Arr] | None:
    """One loop step and the tangent there (oriented along ``t``), or None if the corrector fails
    or (with ``guard``) the A2 acceptance test rejects the step."""
    pc = _loop_step(p, t, step, mu, tol)
    if pc is None:
        if stats is not None:
            stats["corrector"] = stats.get("corrector", 0) + 1
        return None
    tn = _tangent(_loop_eval(pc, mu)[1], t)
    if guard:
        why_not = _step_ok(pc, p + step * t, step, t, tn)
        if why_not is not None:
            if stats is not None:
                stats[why_not] = stats.get(why_not, 0) + 1
            return None
    return pc, tn


def trace_loop(
    x1: float,
    y2: float,
    t_full: float,
    mu: float,
    *,
    direction: Arr,
    ds: float = 5e-3,
    ds_min: float = 1e-6,
    ds_max: float = 2e-2,
    min_members: int = 20,
    max_members: int = 3000,
    tol: float = 1e-10,
    on_member: Callable[[LoopMember], None] | None = None,
    origin: Arr | None = None,
    n_done: int = 0,
    deadline: float | None = None,
    guard: bool = True,
    start_tangent: Arr | None = None,
    stats: dict[str, int] | None = None,
) -> tuple[list[LoopMember], str, float]:
    """Pseudo-arclength continuation of the symmetric family through (x1, y2, T) at fixed mu, in
    p = (x1(0), y2(0), T/2pi), with the analytic half-period Jacobian.

    Returns (members, status, closure): status "closed" when the trace comes back to its start after
    at least ``min_members`` members (closure = the distance of the corrected final member from the
    start in p), "corrector failed", "time budget" or "max members" (closure = nan). To resume a
    checkpointed trace, pass the last member as the start, the first member as ``origin`` (in p) and
    the number of members already traced as ``n_done``.

    ``guard`` applies the A2 step acceptance (#975 note sec. 11.4) to every step, the closing
    landing included. A closure also needs the closing tangent to point the same way as the start
    tangent (dot > 0), so a trace that doubles back is not "closed"; ``start_tangent`` gives it on a
    resume (default: the tangent at the start of this call when it is not a resume, else no check).
    ``stats`` counts rejected steps by reason."""
    p = np.array([x1, y2, t_full / (2.0 * math.pi)])
    p0 = p.copy() if origin is None else np.asarray(origin, dtype=float)
    _, jac = _loop_eval(p, mu)
    t = _tangent(jac, np.asarray(direction, dtype=float))
    if start_tangent is not None:
        t_start: Arr | None = np.asarray(start_tangent, dtype=float)
    else:
        t_start = t.copy() if origin is None else None
    members = [LoopMember(x1, y2, t_full, float(np.max(np.abs(_loop_eval(p, mu)[0]))))]
    for n in range(n_done, max_members):
        if deadline is not None and time.time() > deadline:
            return members, "time budget", float("nan")
        if n + 1 >= min_members and float(np.linalg.norm(p - p0)) < 2.0 * ds:
            # try to land on the origin: repeated projected steps toward it, each corrected
            pc, tc = p, t
            for _ in range(6):
                proj = float(tc @ (p0 - pc))
                if abs(proj) > 2.0 * ds:
                    break
                acc = _accepted_step(pc, tc, proj, mu, tol, guard, stats)
                if acc is None:
                    break
                pc, tc = acc
                if float(np.linalg.norm(pc - p0)) < 1e-9:
                    break
            dist = float(np.linalg.norm(pc - p0))
            forward = t_start is None or float(tc @ t_start) > 0.0
            if dist < 1e-6 and forward:
                m = LoopMember(float(pc[0]), float(pc[1]), 2.0 * math.pi * float(pc[2]), 0.0)
                members.append(m)
                if on_member is not None:
                    on_member(m)
                return members, "closed", dist
            # passed near the origin without landing on it: keep tracing
        step = ds
        acc = None
        while abs(step) >= ds_min:
            acc = _accepted_step(p, t, step, mu, tol, guard, stats)
            if acc is not None:
                break
            step *= 0.5
            ds = step
        if acc is None:
            return members, "corrector failed", float("nan")
        p, t = acc
        m = LoopMember(float(p[0]), float(p[1]), 2.0 * math.pi * float(p[2]), 0.0)
        members.append(m)
        if on_member is not None:
            on_member(m)
        ds = min(ds * 1.3, ds_max)
    return members, "max members", float("nan")


@dataclass(frozen=True)
class MemberGeometry:
    """Distances of a symmetric member over one period, in km at the system's l*."""

    min_r1_km: float
    max_r1_km: float
    min_r2_km: float
    jacobi: float


def member_geometry(
    x1: float, y2: float, t_full: float, system: cr3bp.CR3BPSystem, n: int = 20001
) -> MemberGeometry:
    """Dense DOP853 re-propagation of the member (at ``system.mu``) and its extreme distances from
    P1 and P2 in km (``system.l_km``)."""
    mu = system.mu
    s0 = state_from_crossing(x1, y2, mu)
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, t_full),
        s0,
        args=(mu,),
        method="DOP853",
        rtol=RTOL,
        atol=ATOL,
        dense_output=True,
    )
    assert sol.sol is not None
    z = sol.sol(np.linspace(0.0, t_full, n))
    r1 = np.hypot(z[0] + mu, z[1])
    r2 = np.hypot(z[0] - 1.0 + mu, z[1])
    return MemberGeometry(
        min_r1_km=float(np.min(r1)) * system.l_km,
        max_r1_km=float(np.max(r1)) * system.l_km,
        min_r2_km=float(np.min(r2)) * system.l_km,
        jacobi=float(cr3bp.jacobi_constant(s0, mu)),
    )


def distance_to_loop(
    point: Arr, members: list[Arr], mu: float, *, tol: float = 1e-10, max_refine: int = 8
) -> float:
    """Distance (in p = (x1(0), y2(0), T/2pi)) from ``point`` to the traced family through
    ``members``: from the nearest member, repeated projected pseudo-arclength steps toward the
    point, each corrected onto the family. Small (< 1e-6) when the point lies on that family."""
    pts = np.asarray(members, dtype=float)
    target = np.asarray(point, dtype=float)
    i = int(np.argmin(np.linalg.norm(pts - target, axis=1)))
    p = pts[i]
    t = _tangent(_loop_eval(p, mu)[1], None)
    best = float(np.linalg.norm(p - target))
    for _ in range(max_refine):
        nxt = _loop_step(p, t, float(t @ (target - p)), mu, tol)
        if nxt is None:
            break
        p = nxt
        t = _tangent(_loop_eval(p, mu)[1], t)
        best = min(best, float(np.linalg.norm(p - target)))
        if best < 1e-10:
            break
    return best
