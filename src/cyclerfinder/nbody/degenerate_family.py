"""Seed correction along a degenerate patched-conic family (#1055, Lyapunov-Schmidt reduction).

A full-revolution (resonant) leg leaves a one-parameter family of patched-conic chains. In a
continuous-gravity shooting problem scaled by the moon-mass factor sigma, the family shows as a
singular value of the column-scaled Jacobian proportional to sigma, and the residual's projection
on its left singular vector is also proportional to sigma. The Newton move along the family then
stays O(1) as sigma -> 0, and a damped solver started from the seed stalls. The correction moves the
seed ALONG the family, at fixed sigma, to the member where that projection vanishes.

Formulation-agnostic: a residual function ``rj(z) -> (r, J)`` at fixed sigma (square system).
Variables are column-scaled with fixed scales ``d`` (the Jacobian column norms at the base point).
"""

from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from numpy.typing import NDArray

Arr = NDArray[np.float64]
ResJac = Callable[[Arr], tuple[Arr, Arr]]


def column_scales(j: Arr) -> Arr:
    d = np.linalg.norm(j, axis=0)
    d[d == 0.0] = 1.0
    return np.asarray(d, dtype=np.float64)


def scaled_svd(j: Arr, d: Arr) -> tuple[Arr, Arr, Arr]:
    """(U, s, Vt) of J diag(1/d), singular values in descending order."""
    u, s, vt = np.linalg.svd(j / d)
    return u, s, vt


def detect_degenerate(
    j_hi: Arr,
    j_lo: Arr,
    sigma_ratio: float,
    min_drop: float = 2.5,
    max_other_change: float = 0.25,
    q_max: int = 3,
) -> int:
    """Number of trailing singular values that fall with sigma.

    ``j_hi`` and ``j_lo`` are the Jacobians at sigma and sigma / ``sigma_ratio`` (each scaled by
    its own column norms). A trailing singular value is degenerate if it drops by more than
    ``min_drop`` (proportional to sigma would be ``sigma_ratio``), and the first non-degenerate one
    above them changes by less than ``max_other_change`` (relative).
    """
    s_hi = scaled_svd(j_hi, column_scales(j_hi))[1]
    s_lo = scaled_svd(j_lo, column_scales(j_lo))[1]
    q = 0
    for i in range(1, q_max + 1):
        if s_hi[-i] / s_lo[-i] > min_drop:
            q = i
        else:
            break
    if q == 0:
        return 0
    k = len(s_hi) - q - 1
    if abs(s_lo[k] / s_hi[k] - 1.0) > max_other_change:
        return 0
    return q


@dataclass
class BorderedResult:
    z: Arr
    converged: bool
    iterations: int
    complement_norm: float
    g: Arr
    paused: bool = False


def bordered_solve(
    rj: ResJac,
    z_pred: Arr,
    d: Arr,
    u_null: Arr,
    v_null: Arr,
    z_ref: Arr,
    alpha: Arr,
    tol: float,
    max_iter: int = 25,
    deadline: float | None = None,
) -> BorderedResult:
    """Gauss-Newton on the complement rows (I - U U^T) r = 0 and the q arclength rows
    V^T (d * (z - z_ref)) = alpha. ``u_null`` (m x q) and ``v_null`` (n x q, scaled coordinates,
    orthonormal columns) are held fixed. The arclength rows are linear, so every iterate is
    projected onto them exactly; the line search acts on the complement residual alone. Converged
    when the complement residual (max norm) is below ``tol``."""

    def onto(zz: Arr) -> Arr:
        e = v_null.T @ (d * (zz - z_ref)) - alpha
        return np.asarray(zz - (v_null @ e) / d, dtype=np.float64)

    z = onto(np.asarray(z_pred, dtype=np.float64))
    proj = np.eye(len(u_null)) - u_null @ u_null.T
    it = 0
    r, j = rj(z)
    for it in range(1, max_iter + 1):
        comp = proj @ r
        if float(np.max(np.abs(comp))) < tol:
            return BorderedResult(z, True, it - 1, float(np.max(np.abs(comp))), u_null.T @ r)
        a = np.vstack([proj @ (j / d), v_null.T])
        f = np.concatenate([comp, np.zeros(v_null.shape[1])])
        step = np.linalg.lstsq(a, -f, rcond=None)[0] / d
        f0 = float(np.linalg.norm(comp))
        for lam in [0.5**i for i in range(10)]:
            zt = onto(z + lam * step)
            try:
                rt, jt = rj(zt)
            except RuntimeError:
                continue
            ct = proj @ rt
            if np.all(np.isfinite(ct)) and float(np.linalg.norm(ct)) < f0:
                z, r, j = zt, rt, jt
                break
        else:
            break
        if deadline is not None and time.monotonic() > deadline:  # after at least one step
            cn = float(np.max(np.abs(proj @ r)))
            if cn >= tol:
                return BorderedResult(z, False, it, cn, u_null.T @ r, paused=True)
    comp = proj @ r
    ok = float(np.max(np.abs(comp))) < tol
    return BorderedResult(z, ok, it, float(np.max(np.abs(comp))), u_null.T @ r)


@dataclass
class FamilyFrame:
    """The linearisation at the BASE: the seed with its complement rows converged (#1055
    amendment 1). ``d`` column scales, ``u_n`` / ``v_n`` the q left / right null vectors (scaled),
    ``big_l`` the Newton-step length along them, ``base`` the converged base point."""

    base: BorderedResult
    d: Arr
    u_n: Arr
    v_n: Arr
    big_l: float
    converged: bool


def family_frame(rj: ResJac, z0: Arr, q: int, tol: float) -> FamilyFrame:
    """Converge the complement at the seed with the seed's null vectors, then re-linearise at the
    converged point and converge again with the base's own vectors (arclength rows about it)."""
    _, j0 = rj(z0)
    d = column_scales(j0)
    u, s, vt = scaled_svd(j0, d)
    first = bordered_solve(rj, z0, d, u[:, -q:], vt[-q:].T, z0, np.zeros(q), tol)
    if not first.converged:
        return FamilyFrame(first, d, u[:, -q:], vt[-q:].T, 0.0, False)
    rb_, jb = rj(first.z)
    d = column_scales(jb)
    u, s, vt = scaled_svd(jb, d)
    u_n, v_n = u[:, -q:], vt[-q:].T
    base = bordered_solve(rj, first.z, d, u_n, v_n, first.z, np.zeros(q), tol)
    big_l = float(np.linalg.norm((u_n.T @ rb_) / s[-q:]))
    return FamilyFrame(base, d, u_n, v_n, big_l, base.converged)


@dataclass
class FamilyResult:
    z_star: Arr | None
    alpha_star: float | None
    status: str
    scan: list[tuple[float, float, float]] = field(default_factory=list)  # (alpha, g, |comp|)
    newton_length: float = 0.0


def family_root_1d(
    rj: ResJac,
    z0: Arr,
    tol: float,
    max_lengths: float = 3.0,
    h_frac: float = 0.01,
    log: Callable[[str], None] | None = None,
    g_floor: float = 0.0,
    state: dict[str, Any] | None = None,
    save: Callable[[dict[str, Any]], None] | None = None,
    deadline: float | None = None,
) -> FamilyResult:
    """q = 1: converge the complement at the base, scan alpha both ways along the fixed null
    direction for a sign change of g = u . r (the base's u), then a secant to
    |g| < max(1e-3 |g(0)|, g_floor).

    The arclength alpha is in scaled units along v (the base's right null vector). The scan length
    is ``max_lengths`` times the Newton-step length L = |u . r| / s_min at the base; the first step
    is ``h_frac`` L, doubled after a corrector converging in 4 iterations or fewer, halved on
    failure (down to 1e-4 L).

    Resumable: ``state`` (a JSON-serialisable dict, mutated in place) holds the phase, the frame,
    the scan position and the bracket; ``save`` is called after every step; past ``deadline``
    (time.monotonic) the function returns with status "paused" and a later call with the same
    ``state`` continues."""
    say = log or (lambda _m: None)
    st: dict[str, Any] = state if state is not None else {}
    put = save or (lambda _s: None)
    res = FamilyResult(None, None, "paused")

    progressed = False

    def late() -> bool:
        return progressed and deadline is not None and time.monotonic() > deadline

    def arr(key: str) -> Arr:
        return np.asarray(st[key], dtype=np.float64)

    if "phase" not in st:
        _, j0 = rj(z0)
        d = column_scales(j0)
        u, _s, vt = scaled_svd(j0, d)
        st.update(phase="frame1", d=d.tolist(), u=u[:, -1].tolist(), v=vt[-1].tolist())
        st.update(z=np.asarray(z0, dtype=np.float64).tolist(), z_seed=list(map(float, z0)))
        put(st)
    while True:
        ph = st["phase"]
        d = arr("d")
        u_n, v_n = arr("u")[:, None], arr("v")[:, None]
        if ph in ("frame1", "frame2"):
            zr = arr("z") if ph == "frame2" else arr("z_seed")
            br = bordered_solve(rj, arr("z"), d, u_n, v_n, zr, np.zeros(1), tol, deadline=deadline)
            st["z"] = br.z.tolist()
            if br.paused:
                put(st)
                return res
            if not br.converged:
                st["phase"] = "done"
                st["status"] = "base complement did not converge"
                put(st)
                return FamilyResult(None, None, st["status"])
            if ph == "frame1":
                rb_, jb = rj(br.z)
                d = column_scales(jb)
                u, s_, vt = scaled_svd(jb, d)
                st.update(d=d.tolist(), u=u[:, -1].tolist(), v=vt[-1].tolist())
                st["big_l"] = abs(float(u[:, -1] @ rb_)) / float(s_[-1])
                st["phase"] = "frame2"
            else:
                g0 = float(br.g[0])
                st.update(base=br.z.tolist(), g0=g0, scan=[[0.0, g0, br.complement_norm]])
                say(f"family: L {st['big_l']:.3e} g(0) {g0:.3e} comp {br.complement_norm:.2e}")
                if abs(g0) <= tol:
                    st.update(phase="done", status="root at base", z_star=br.z.tolist(), a_star=0.0)
                else:
                    st.update(phase="scan", sgn=1.0, a_prev=0.0, z_prev=br.z.tolist(), g_prev=g0)
                    st.update(dzda=(v_n[:, 0] / d).tolist(), h=h_frac * st["big_l"])
            put(st)
        elif ph == "scan":
            big_l = float(st["big_l"])
            sgn, a_prev, h = float(st["sgn"]), float(st["a_prev"]), float(st["h"])
            if abs(a_prev) >= max_lengths * big_l or h < 1e-4 * big_l:
                if sgn > 0:
                    st.update(sgn=-1.0, a_prev=0.0, z_prev=st["base"], g_prev=st["g0"])
                    st.update(dzda=(v_n[:, 0] / d).tolist(), h=h_frac * big_l)
                else:
                    st.update(phase="done", status="no sign change of g within the scanned range")
                put(st)
                continue
            if late():
                return res
            a = a_prev + sgn * h
            z_prev = arr("z_prev")
            pred = z_prev + arr("dzda") * (a - a_prev)
            br = bordered_solve(rj, pred, d, u_n, v_n, arr("base"), np.array([a]), tol)
            if not br.converged:
                st["h"] = 0.5 * h
                put(st)
                progressed = True
                continue
            g = float(br.g[0])
            st["scan"].append([a, g, br.complement_norm])
            say(f"family: alpha {a:.4e} g {g:.4e} comp {br.complement_norm:.2e} it {br.iterations}")
            if np.sign(g) != np.sign(float(st["g_prev"])):
                st.update(phase="secant", lo=[a_prev, st["z_prev"], st["g_prev"]])
                st["hi"] = [a, br.z.tolist(), g]
            else:
                st["dzda"] = ((br.z - z_prev) / (a - a_prev)).tolist()
                st.update(a_prev=a, z_prev=br.z.tolist(), g_prev=g)
                if br.iterations <= 4:
                    st["h"] = 2.0 * h
            put(st)
            progressed = True
        elif ph == "secant":
            if late():
                return res
            a_lo, z_lo, g_lo = st["lo"]
            a_hi, z_hi, g_hi = st["hi"]
            z_lo, z_hi = np.asarray(z_lo), np.asarray(z_hi)
            a = a_hi - g_hi * (a_hi - a_lo) / (g_hi - g_lo)
            w = (a - a_lo) / (a_hi - a_lo)
            br = bordered_solve(
                rj, z_lo + w * (z_hi - z_lo), d, u_n, v_n, arr("base"), np.array([a]), tol
            )
            st["n_secant"] = int(st.get("n_secant", 0)) + 1
            if not br.converged:
                st.update(phase="done", status="secant corrector failed")
            else:
                g = float(br.g[0])
                st["scan"].append([a, g, br.complement_norm])
                say(f"family: secant alpha {a:.4e} g {g:.4e}")
                if abs(g) < max(1e-3 * abs(float(st["g0"])), g_floor):
                    st.update(phase="done", status="root", z_star=br.z.tolist(), a_star=a)
                elif st["n_secant"] >= 40:
                    st.update(phase="done", status="secant did not converge")
                elif np.sign(g) == np.sign(g_lo):
                    st["lo"] = [a, br.z.tolist(), g]
                else:
                    st["hi"] = [a, br.z.tolist(), g]
            put(st)
            progressed = True
        else:  # done
            out = FamilyResult(None, None, st["status"], newton_length=float(st.get("big_l", 0.0)))
            out.scan = [(float(a), float(g), float(c)) for a, g, c in st.get("scan", [])]
            if "z_star" in st:
                out.z_star = arr("z_star")
                out.alpha_star = float(st["a_star"])
            return out


def family_root_nd(
    rj: ResJac,
    z0: Arr,
    q: int,
    tol: float,
    max_lengths: float = 3.0,
    max_iter: int = 30,
    log: Callable[[str], None] | None = None,
    g_floor: float = 0.0,
) -> FamilyResult:
    """q >= 1: damped Newton on the reduced equations g(alpha) = U^T r = 0, alpha in R^q (scaled
    arclength along the base's q right null vectors). g's Jacobian by forward differences (step
    1 % of L, the Newton-step length of the base); steps limited to L; stop at |g| < 1e-3 |g(0)|,
    or when |alpha| exceeds ``max_lengths`` L."""
    say = log or (lambda _m: None)
    fr = family_frame(rj, z0, q, tol)
    d, u_n, v_n, big_l, cur = fr.d, fr.u_n, fr.v_n, fr.big_l, fr.base
    res = FamilyResult(None, None, "base complement did not converge", newton_length=big_l)
    alpha = np.zeros(q)
    if not fr.converged:
        return res
    z0 = cur.z  # arclength is measured from the base
    g0n = float(np.linalg.norm(cur.g))
    say(f"family nd: q {q} L {big_l:.3e} |g(0)| {g0n:.3e}")
    if g0n <= tol:
        res.z_star, res.alpha_star, res.status = cur.z, 0.0, "root at base"
        return res
    for _ in range(max_iter):
        gn = float(np.linalg.norm(cur.g))
        res.scan.append((float(np.linalg.norm(alpha)), gn, cur.complement_norm))
        if gn < max(1e-3 * g0n, g_floor):
            res.z_star, res.alpha_star, res.status = cur.z, float(np.linalg.norm(alpha)), "root"
            return res
        if float(np.linalg.norm(alpha)) > max_lengths * big_l:
            res.status = "no root within the scanned range"
            return res
        h = 0.01 * big_l
        jg = np.zeros((q, q))
        for i in range(q):
            da = np.zeros(q)
            da[i] = h
            b = bordered_solve(rj, cur.z + (v_n[:, i] * h) / d, d, u_n, v_n, z0, alpha + da, tol)
            if not b.converged:
                res.status = "finite-difference corrector failed"
                return res
            jg[:, i] = (b.g - cur.g) / h
        step = np.linalg.lstsq(jg, -cur.g, rcond=None)[0]
        if float(np.linalg.norm(step)) > big_l:
            step *= big_l / float(np.linalg.norm(step))
        for lam in [0.5**i for i in range(8)]:
            a_new = alpha + lam * step
            b = bordered_solve(rj, cur.z + (v_n @ (lam * step)) / d, d, u_n, v_n, z0, a_new, tol)
            if b.converged and float(np.linalg.norm(b.g)) < gn:
                alpha, cur = a_new, b
                say(f"family nd: |alpha| {np.linalg.norm(alpha):.3e} |g| {np.linalg.norm(b.g):.3e}")
                break
        else:
            res.status = "line search failed"
            return res
    res.status = "Newton did not converge"
    return res
