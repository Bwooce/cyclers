"""#668 -- validated (rigorous, interval-enclosed) Taylor ODE integrator.

This is the single load-bearing primitive that ``#653``'s Wilczak-Zgliczynski
(W-Z) proof-machinery scoping identified as *the* thing this codebase does not
have: a rigorous ODE integrator that returns a mathematically guaranteed
*enclosure* of the true solution, not a floating-point approximation with an
error estimate.  ``#403``'s Oterma golden
(``data/golden/wz_oterma_heteroclinic.yaml``) reproduced W-Z's published
crossing coordinates in ordinary floating point; W-Z's *proof* rests on rigorous
interval enclosures of Poincare section maps produced by the CAPD C++ library's
C0/C1-Lohner algorithms.  There is no maintained Python binding to CAPD (the
PyPI ``capd`` name is an empty placeholder), and this ``uv``-managed environment
has no Julia / ``TaylorIntegration.jl`` route either (``#662``).  The only
rigorous substrate available is ``mpmath.iv`` (directed-rounding interval
arithmetic), already vetted in-repo by the ``#610``/``#625`` closed-form
bend-gate certificates -- so this module builds the validated integrator on it
from scratch.

Method -- Lohner's C0 algorithm / the "high-order enclosure" (interval Taylor)
method (R. J. Lohner, *Computation of guaranteed enclosures for the solutions of
ODEs*, 1992; N. S. Nedialkov, K. R. Jackson & G. F. Corliss, *Validated
solutions of IVPs for ODEs*, Appl. Math. Comput. 105 (1999) 21-68; R. E. Moore,
*Interval Analysis*, 1966).  One validated step from ``t0`` to ``t0 + h``:

  1. **A-priori (rough) enclosure.**  Find a box ``[W]`` and *verify* the
     Picard/Banach inclusion  ``[y0] + [0, h] * f([W])  subset of  [W]``.  By the
     Schauder fixed-point theorem this inclusion is a rigorous certificate that
     the true solution ``y(t)`` stays in ``[W]`` for every ``t in [t0, t0+h]``.
     Not "probably" -- the inclusion is checked with interval endpoints, so if it
     holds it holds for the exact flow.

  2. **Taylor coefficients.**  Bootstrap the solution's Taylor coefficients
     ``y_0, y_1, ..., y_p`` order by order via the standard autonomous-ODE
     recurrence  ``y_{k+1} = F_k(y_0..y_k) / (k+1)``  where ``F_k`` is the k-th
     Taylor coefficient of ``f(y(.))`` -- computed here in interval Taylor
     arithmetic (Cauchy products and a real-power recurrence), evaluated at the
     *point* interval ``[y0]``.

  3. **Rigorous remainder.**  Re-run the same recurrence one order further with
     the *rough enclosure* ``[W]`` as the 0-th coefficient, and take its order
     ``p+1`` coefficient ``R``.  Then

         [y1]  =  sum_{k=0}^{p} y_k([y0]) * h^k   +   R([W]) * h^{p+1}

     is a rigorous enclosure of ``y(t0+h)``: the Lagrange remainder's
     ``(p+1)``-th derivative is evaluated at some ``xi in [t0, t0+h]`` whose state
     is provably in ``[W]``, so ``R([W])`` bounds it for the exact solution.

Honest scope -- what this is and is NOT (carried into the driver's reporting too):

  * This is a genuine, interval-*rigorous* C0 one-step + multi-step integrator.
    Its enclosures are guaranteed, validated against closed-form goldens
    (exp, harmonic oscillator) with a positive control proving the enclosure
    both *contains* the true value and *excludes* nearby wrong values (it is not
    a vacuous ``[-inf, inf]`` bound).
  * It is NOT yet the full W-Z proof.  The *naive C0* multi-step ``integrate``
    keeps its enclosure axis-aligned in the fixed frame, so over long
    integrations near the ~2e4x/period L1/L2 saddle the enclosure suffers the
    classical *wrapping effect* and eventually blows up (quantified honestly by
    the ``#668`` driver).

  * ``#669`` adds the standard wrapping-effect fix -- **Lohner's C1 algorithm
    with QR-coordinate reframing** (``integrate_c1_qr`` below): the enclosure is
    carried as ``yhat + A @ [r]`` where ``A`` is a point orthogonal frame
    re-chosen each step (QR of the propagated variational/Jacobian flow) and
    ``[r]`` is a box in that *rotating* frame, so spurious volume no longer
    accumulates under the flow's rotation/shear.  This is genuine and effective:
    on the textbook wrapping problem (rigid rotation of a wide box, ~3 turns) the
    naive ``integrate`` blows the enclosure up by ~5e8x while ``integrate_c1_qr``
    holds it essentially flat (see ``#669``'s tests / driver).  QR is the
    ingredient the future wide-h-set covering-relations stage *requires*.

    Honest caveat, though: on the ``#668`` point-trajectory Oterma L1->L2 arc, QR
    does NOT extend the reach -- because that arc is not wrapping-limited.  Its
    true flow stretching is modest (``||STM||_2 ~ 3e2`` by ``t ~ 0.45``) and does
    not rotate strongly against the fixed frame, so wrapping is minor and QR ~=
    naive C0 there.  The real horizon is a *physical* close flyby of the secondary
    (Jupiter perijove ``r2 ~ 0.0045`` at ``t ~ 0.465``), where the near-collision
    field's high derivatives and the enclosure width surpassing the miss distance
    stop the a-priori box -- a bottleneck QR cannot and should not fix (that needs
    close-approach regularisation).  ``integrate_c1_qr`` still does NOT provide
    rigorous Poincare-section-map derivatives or the covering-relations / h-set
    machinery -- those remain future stages.

All arithmetic is ``mpmath.iv`` at ``mp.iv.dps`` precision set by the caller
(the QR frame choice and approximate inverse are done in ``float`` -- they only
*choose* the coordinate system; rigor is re-established by enclosing every
propagation step, and the frame inverse, in interval arithmetic).
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, TypeAlias

import numpy as np

try:
    import mpmath as mp  # noqa: F401 -- availability probe (HAVE_MPMATH)

    HAVE_MPMATH = True
except ImportError:  # pragma: no cover - exercised by the skip-clean test
    HAVE_MPMATH = False


# --------------------------------------------------------------------------- #
# Interval Taylor-series arithmetic.                                           #
#                                                                             #
# A series is a plain Python list of mpmath.iv intervals [c0, c1, ..., cp]     #
# representing  sum_k c_k * (t - t0)^k  truncated at order p.  Every operation #
# preserves the invariant that coefficient k of the result depends only on     #
# coefficients <= k of the inputs -- this is what lets the ODE recurrence      #
# bootstrap order by order.                                                    #
# --------------------------------------------------------------------------- #
Series: TypeAlias = list[Any]  # list of mpmath.iv intervals [c0, c1, ..., cp]


def ts_const(iv: Any, value: Any, order: int) -> Series:
    """Constant series: coefficient 0 is ``value`` (an iv interval), rest zero."""
    zero = iv.mpf(0)
    return [value if k == 0 else zero for k in range(order + 1)]


def ts_add(iv: Any, a: Series, b: Series) -> Series:
    return [a[k] + b[k] for k in range(len(a))]


def ts_sub(iv: Any, a: Series, b: Series) -> Series:
    return [a[k] - b[k] for k in range(len(a))]


def ts_scale(iv: Any, a: Series, s: Any) -> Series:
    """Multiply a whole series by an interval (or real) scalar ``s``."""
    return [a[k] * s for k in range(len(a))]


def ts_add_scalar(iv: Any, a: Series, s: Any) -> Series:
    """Add an interval/real scalar ``s`` to the 0-th coefficient only."""
    out = list(a)
    out[0] = out[0] + iv.mpf(s) if not isinstance(s, type(a[0])) else out[0] + s
    return out


def ts_mul(iv: Any, a: Series, b: Series) -> Series:
    """Truncated Cauchy product:  c_k = sum_{j=0}^{k} a_j b_{k-j}."""
    n = len(a)
    zero = iv.mpf(0)
    out = [zero for _ in range(n)]
    for k in range(n):
        acc = zero
        for j in range(k + 1):
            acc = acc + a[j] * b[k - j]
        out[k] = acc
    return out


def ts_pow(iv: Any, a: Series, alpha: float) -> Series:
    """Series real power  u = a**alpha  via the standard J.C.P. Miller recurrence.

    Requires ``0`` not in ``a[0]`` (the base's constant term).  For the CR3BP
    ``s**(-3/2)`` distance-cubed terms this always holds (a squared distance is
    strictly positive off-collision).

        u_0 = a_0**alpha
        u_k = (1 / (k * a_0)) * sum_{j=0}^{k-1} (alpha*(k - j) - j) * a_{k-j} * u_j
    """
    n = len(a)
    a0 = a[0]
    # a0**alpha for an interval base via exp(alpha*log(a0)); a0 > 0 required.
    u0 = iv.exp(iv.mpf(alpha) * iv.log(a0))
    zero = iv.mpf(0)
    u = [zero for _ in range(n)]
    u[0] = u0
    for k in range(1, n):
        acc = zero
        for j in range(k):
            coeff = iv.mpf(alpha) * (k - j) - iv.mpf(j)
            acc = acc + coeff * a[k - j] * u[j]
        u[k] = acc / (iv.mpf(k) * a0)
    return u


# --------------------------------------------------------------------------- #
# Planar CR3BP vector field, as an interval-Taylor "jet".                      #
# --------------------------------------------------------------------------- #
def cr3bp_planar_jet(iv: Any, state: list[Series], mu: float) -> list[Series]:
    """Return the derivative series [X', Y', VX', VY'] from state series.

    ``state`` is ``[X, Y, VX, VY]`` (each a truncated interval Taylor series).
    Mirrors :func:`cyclerfinder.core.cr3bp.cr3bp_eom` restricted to the plane
    (z = vz = 0):

        x'  = vx
        y'  = vy
        vx' = x + 2 vy - (1-mu)(x+mu) s1^{-3/2} - mu (x-1+mu) s2^{-3/2}
        vy' = y - 2 vx - (1-mu) y     s1^{-3/2} - mu     y   s2^{-3/2}

    with  s1 = (x+mu)^2 + y^2,  s2 = (x-1+mu)^2 + y^2.
    """
    X, Y, VX, VY = state
    om1 = 1.0 - mu

    xm1 = ts_add_scalar(iv, X, iv.mpf(mu))  # x + mu
    xm2 = ts_add_scalar(iv, X, iv.mpf(mu - 1.0))  # x - 1 + mu
    y2 = ts_mul(iv, Y, Y)
    s1 = ts_add(iv, ts_mul(iv, xm1, xm1), y2)
    s2 = ts_add(iv, ts_mul(iv, xm2, xm2), y2)
    r1i3 = ts_pow(iv, s1, -1.5)
    r2i3 = ts_pow(iv, s2, -1.5)

    # gravitational terms
    g1x = ts_scale(iv, ts_mul(iv, xm1, r1i3), iv.mpf(om1))
    g2x = ts_scale(iv, ts_mul(iv, xm2, r2i3), iv.mpf(mu))
    g1y = ts_scale(iv, ts_mul(iv, Y, r1i3), iv.mpf(om1))
    g2y = ts_scale(iv, ts_mul(iv, Y, r2i3), iv.mpf(mu))

    ax = ts_sub(iv, ts_sub(iv, ts_add(iv, X, ts_scale(iv, VY, iv.mpf(2.0))), g1x), g2x)
    ay = ts_sub(iv, ts_sub(iv, ts_sub(iv, Y, ts_scale(iv, VX, iv.mpf(2.0))), g1y), g2y)
    return [VX, VY, ax, ay]


def cr3bp_planar_jacobi(iv: Any, state6: list[Any], mu: float) -> Any:
    """Interval enclosure of the planar Jacobi constant C over a state box.

    ``state6`` is ``[x, y, vx, vy]`` intervals.  Closed-form, no integration:
    ``C = (x^2+y^2) + 2(1-mu)/r1 + 2 mu/r2 - (vx^2+vy^2)``.
    """
    x, y, vx, vy = state6
    om1 = iv.mpf(1.0 - mu)
    r1 = iv.sqrt((x + iv.mpf(mu)) ** 2 + y**2)
    r2 = iv.sqrt((x + iv.mpf(mu - 1.0)) ** 2 + y**2)
    return (x**2 + y**2) + 2 * om1 / r1 + 2 * iv.mpf(mu) / r2 - (vx**2 + vy**2)


# --------------------------------------------------------------------------- #
# The validated integrator.                                                    #
# --------------------------------------------------------------------------- #
def _f_over_box(iv: Any, jet: Any, box: list[Any], mu: float) -> list[Any]:
    """Order-0 interval evaluation f([box]): coeff 0 of the jet on constant series."""
    order0 = [ts_const(iv, box[i], 0) for i in range(len(box))]
    deriv = jet(iv, order0, mu)
    return [d[0] for d in deriv]


def _subset(iv: Any, inner: Any, outer: Any) -> bool:
    """Rigorous interval inclusion  inner subset of outer  (endpoint compare)."""
    return bool(outer.a <= inner.a) and bool(inner.b <= outer.b)


def apriori_enclosure(
    iv: Any,
    jet: Any,
    y0: list[Any],
    h: float,
    mu: float,
    *,
    max_iter: int = 60,
    inflate: float = 1.2,
) -> list[Any] | None:
    """Rough a-priori enclosure ``[W]`` with the Picard inclusion verified.

    Returns ``[W]`` (a list of intervals) such that the true solution provably
    stays inside it over ``[0, h]``, or ``None`` if the step ``h`` is too large
    to validate within ``max_iter`` inflations.
    """
    dim = len(y0)
    hint = iv.mpf([0.0, h])
    try:
        fy0 = _f_over_box(iv, jet, y0, mu)
    except (ValueError, ArithmeticError):  # vector field undefined at IC
        return None
    W = [y0[i] + hint * fy0[i] for i in range(dim)]
    for _ in range(max_iter):
        # inflate W around each midpoint
        Wi = []
        for w in W:
            m = w.mid
            r = (w.b - w.a) / 2
            pad = r * iv.mpf(inflate) + iv.mpf(1e-30)
            Wi.append(iv.mpf([m.a - pad.b, m.b + pad.b]))
        # A widened box can cross a collision manifold (s1/s2 <= 0), making the
        # r^{-3/2} terms undefined -- that is a rigorous "cannot validate this
        # step" signal (the enclosure has grown too large), reported as None.
        try:
            fW = _f_over_box(iv, jet, Wi, mu)
        except (ValueError, ArithmeticError):
            return None
        Wnew = [y0[i] + hint * fW[i] for i in range(dim)]
        if all(_subset(iv, Wnew[i], Wi[i]) for i in range(dim)):
            return Wi
        W = Wnew
    return None


def solution_series(iv: Any, jet: Any, y0: list[Any], order: int, mu: float) -> list[Series]:
    """Bootstrap solution Taylor coefficients to ``order`` from IC ``y0``.

    ``y0`` may be a point box (``[y0]``) for the polynomial part, or the rough
    enclosure (``[W]``) for the remainder term -- the recurrence is identical.
    """
    dim = len(y0)
    state = [ts_const(iv, y0[i], order) for i in range(dim)]
    for k in range(order):
        deriv = jet(iv, state, mu)
        for i in range(dim):
            state[i][k + 1] = deriv[i][k] / iv.mpf(k + 1)
    return state


def validated_step(
    iv: Any,
    jet: Any,
    y0: list[Any],
    h: float,
    mu: float,
    *,
    order: int = 12,
) -> list[Any] | None:
    """One rigorous step: enclosure of ``y(t0 + h)`` given ``y(t0) in y0``.

    Returns the enclosure list, or ``None`` if the a-priori enclosure failed
    (caller should shrink ``h``).
    """
    W = apriori_enclosure(iv, jet, y0, h, mu)
    if W is None:
        return None
    dim = len(y0)
    poly = solution_series(iv, jet, y0, order, mu)
    rem = solution_series(iv, jet, W, order + 1, mu)
    hpow = [iv.mpf(1)]
    hi = iv.mpf(h)
    for _ in range(order + 1):
        hpow.append(hpow[-1] * hi)
    out = []
    for i in range(dim):
        acc = iv.mpf(0)
        for k in range(order + 1):
            acc = acc + poly[i][k] * hpow[k]
        acc = acc + rem[i][order + 1] * hpow[order + 1]
        out.append(acc)
    return out


def integrate(
    iv: Any,
    jet: Any,
    y0: list[Any],
    t_final: float,
    mu: float,
    *,
    n_steps: int,
    order: int = 12,
) -> dict[str, Any]:
    """Multi-step validated integration over ``[0, t_final]`` in ``n_steps``.

    Returns a dict with the final enclosure, per-step max half-width (the
    wrapping-effect growth curve), and whether the whole run stayed validated.
    """
    h = t_final / n_steps
    y = list(y0)
    widths = []
    ok = True
    for _ in range(n_steps):
        nxt = validated_step(iv, jet, y, h, mu, order=order)
        if nxt is None:
            ok = False
            break
        y = nxt
        # rigorous upper bound on each component's half-width via .delta (width)
        widths.append(max(float(c.delta.b) for c in y) / 2.0)
    return {
        "final": y,
        "validated": ok,
        "n_completed": len(widths),
        "max_halfwidths": widths,
        "final_max_halfwidth": widths[-1] if widths else None,
        "step_h": h,
        "order": order,
    }


# --------------------------------------------------------------------------- #
# C1 variational jet: augmented state [y (n) ; V (n*n, row-major)] carrying the #
# state-transition matrix V = d phi_t / d y0 alongside the state.  The V block  #
# obeys  V' = Df(y) . V,  V(0) = I,  so V(h) = D phi_h(y0) is the flow          #
# Jacobian -- the object the Lohner QR reframing needs.  Same jet primitives.   #
# --------------------------------------------------------------------------- #
def cr3bp_planar_variational_jet(iv: Any, aug: list[Series], mu: float) -> list[Series]:
    """Augmented derivative series for the planar CR3BP state + 4x4 STM.

    ``aug`` is ``[X, Y, VX, VY, V00, V01, ..., V33]`` (4 state + 16 STM series,
    the STM stored row-major).  Returns the 20 derivative series.  The state
    block matches :func:`cr3bp_planar_jet`; the STM block is ``V' = Df(y) V`` with
    the pseudo-potential second derivatives ``Uxx, Uxy, Uyy`` mirroring
    :func:`cyclerfinder.core.cr3bp.cr3bp_stm_eom` restricted to the plane.
    """
    order = len(aug[0]) - 1
    X, Y, VX, VY = aug[0], aug[1], aug[2], aug[3]
    V = aug[4:20]  # 16 series, row-major 4x4
    om1 = 1.0 - mu

    xm1 = ts_add_scalar(iv, X, iv.mpf(mu))  # x + mu
    xm2 = ts_add_scalar(iv, X, iv.mpf(mu - 1.0))  # x - 1 + mu
    y2 = ts_mul(iv, Y, Y)
    s1 = ts_add(iv, ts_mul(iv, xm1, xm1), y2)
    s2 = ts_add(iv, ts_mul(iv, xm2, xm2), y2)
    r1i3 = ts_pow(iv, s1, -1.5)
    r2i3 = ts_pow(iv, s2, -1.5)
    r1i5 = ts_pow(iv, s1, -2.5)
    r2i5 = ts_pow(iv, s2, -2.5)

    # state derivatives (identical to cr3bp_planar_jet)
    g1x = ts_scale(iv, ts_mul(iv, xm1, r1i3), iv.mpf(om1))
    g2x = ts_scale(iv, ts_mul(iv, xm2, r2i3), iv.mpf(mu))
    g1y = ts_scale(iv, ts_mul(iv, Y, r1i3), iv.mpf(om1))
    g2y = ts_scale(iv, ts_mul(iv, Y, r2i3), iv.mpf(mu))
    ax = ts_sub(iv, ts_sub(iv, ts_add(iv, X, ts_scale(iv, VY, iv.mpf(2.0))), g1x), g2x)
    ay = ts_sub(iv, ts_sub(iv, ts_sub(iv, Y, ts_scale(iv, VX, iv.mpf(2.0))), g1y), g2y)

    # pseudo-potential second derivatives (series), 3*om1 / 3*mu prefactors
    base = ts_sub(
        iv,
        ts_sub(iv, ts_const(iv, iv.mpf(1.0), order), ts_scale(iv, r1i3, iv.mpf(om1))),
        ts_scale(iv, r2i3, iv.mpf(mu)),
    )
    uxx = ts_add(
        iv,
        base,
        ts_add(
            iv,
            ts_scale(iv, ts_mul(iv, ts_mul(iv, xm1, xm1), r1i5), iv.mpf(3.0 * om1)),
            ts_scale(iv, ts_mul(iv, ts_mul(iv, xm2, xm2), r2i5), iv.mpf(3.0 * mu)),
        ),
    )
    uyy = ts_add(
        iv,
        base,
        ts_add(
            iv,
            ts_scale(iv, ts_mul(iv, y2, r1i5), iv.mpf(3.0 * om1)),
            ts_scale(iv, ts_mul(iv, y2, r2i5), iv.mpf(3.0 * mu)),
        ),
    )
    uxy = ts_add(
        iv,
        ts_scale(iv, ts_mul(iv, ts_mul(iv, xm1, Y), r1i5), iv.mpf(3.0 * om1)),
        ts_scale(iv, ts_mul(iv, ts_mul(iv, xm2, Y), r2i5), iv.mpf(3.0 * mu)),
    )

    two = iv.mpf(2.0)
    mtwo = iv.mpf(-2.0)
    dV: list[Series] = [[] for _ in range(16)]
    for c in range(4):  # column of V
        v0, v1, v2, v3 = V[c], V[4 + c], V[8 + c], V[12 + c]
        dV[c] = v2  # row 0 of Df = [0,0,1,0]
        dV[4 + c] = v3  # row 1 of Df = [0,0,0,1]
        # row 2 of Df = [Uxx, Uxy, 0, 2]
        dV[8 + c] = ts_add(
            iv, ts_add(iv, ts_mul(iv, uxx, v0), ts_mul(iv, uxy, v1)), ts_scale(iv, v3, two)
        )
        # row 3 of Df = [Uxy, Uyy, -2, 0]
        dV[12 + c] = ts_add(
            iv, ts_add(iv, ts_mul(iv, uxy, v0), ts_mul(iv, uyy, v1)), ts_scale(iv, v2, mtwo)
        )
    return [VX, VY, ax, ay, *dV]


# --------------------------------------------------------------------------- #
# Rigorous interval linear algebra (small dense matrices as list-of-lists).     #
# --------------------------------------------------------------------------- #
IMat: TypeAlias = list[list[Any]]  # interval matrix


def _mid(x: Any) -> Any:
    """Midpoint of an interval as a *thin* interval (tiny rounding width).

    Kept as an interval (not cast to a bare float/mpf): the ``yhat + A[r]``
    representation stays rigorous with a thin-interval centre, whose ~1e-40 width
    simply folds into the local remainder -- no precision is lost, and no fragile
    float cast of a non-degenerate interval is attempted.
    """
    return (x.a + x.b) / 2


def _midf(x: Any) -> float:
    """Float midpoint of an interval, for the (heuristic-only) numpy frame ops."""
    return 0.5 * (float(x.a) + float(x.b))


def _mag(x: Any) -> Any:
    """Rigorous sup |x| of an interval, as an mpf real."""
    lo = -x.a if x.a < 0 else x.a
    hi = -x.b if x.b < 0 else x.b
    return hi if hi > lo else lo


def _imatmul(iv: Any, a: IMat, b: IMat) -> IMat:
    n, k, m = len(a), len(b), len(b[0])
    out: IMat = [[iv.mpf(0) for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for j in range(m):
            acc = iv.mpf(0)
            for t in range(k):
                acc = acc + a[i][t] * b[t][j]
            out[i][j] = acc
    return out


def _imatvec(iv: Any, a: IMat, v: list[Any]) -> list[Any]:
    n, k = len(a), len(v)
    out: list[Any] = []
    for i in range(n):
        acc = iv.mpf(0)
        for t in range(k):
            acc = acc + a[i][t] * v[t]
        out.append(acc)
    return out


def rigorous_inverse(iv: Any, a_pt: IMat) -> IMat | None:
    """Rigorous interval enclosure of the inverse of a *point* matrix ``a_pt``.

    ``a_pt`` is a thin (point) interval matrix.  Computes an approximate inverse
    ``B`` in float, then encloses the true inverse via the Neumann bound: with
    ``E = I - B A`` and ``q = ||E||_inf < 1`` (verified rigorously),

        A^{-1} = (I - E)^{-1} B = B + E B + sum_{k>=2} E^k B,

    and each entry of the tail is bounded by ``bmax * q^2 / (1 - q)`` where
    ``bmax = max |B_ij|``.  Returns ``None`` if ``q < 1`` cannot be certified
    (the frame is too ill-conditioned to invert rigorously -- caller fails).
    """
    n = len(a_pt)
    a_float = np.array([[_midf(a_pt[i][j]) for j in range(n)] for i in range(n)])
    try:
        b_float = np.linalg.inv(a_float)
    except np.linalg.LinAlgError:  # pragma: no cover - singular frame
        return None
    biv: IMat = [[iv.mpf(float(b_float[i][j])) for j in range(n)] for i in range(n)]
    ba = _imatmul(iv, biv, a_pt)
    e: IMat = [
        [(iv.mpf(1) if i == j else iv.mpf(0)) - ba[i][j] for j in range(n)] for i in range(n)
    ]
    # rigorous ||E||_inf upper bound
    q = iv.mpf(0)
    for i in range(n):
        rs = iv.mpf(0)
        for j in range(n):
            rs = rs + iv.mpf([0, _mag(e[i][j])])
        if rs.b > q.b:
            q = rs
    if not bool(q.b < 1):
        return None
    bmax = iv.mpf(0)
    for i in range(n):
        for j in range(n):
            m = _mag(biv[i][j])
            if m > bmax.a:
                bmax = iv.mpf([m, m])
    qb = iv.mpf([q.b, q.b])
    delta = (bmax * qb * qb / (iv.mpf(1) - qb)).b  # rigorous per-entry tail bound
    eb = _imatmul(iv, e, biv)
    return [[biv[i][j] + eb[i][j] + iv.mpf([-delta, delta]) for j in range(n)] for i in range(n)]


# --------------------------------------------------------------------------- #
# Lohner C1 algorithm with QR-coordinate reframing (the wrapping-effect fix).   #
# --------------------------------------------------------------------------- #
def _qr_step(
    iv: Any,
    jet: Callable[..., Any],
    var_jet: Callable[..., Any],
    yhat: list[Any],
    a_frame: IMat,
    rbox: list[Any],
    h: float,
    mu: float,
    order: int,
) -> tuple[list[Any], IMat, list[Any]] | None:
    """One Lohner-C1/QR step of the enclosure ``yhat + a_frame @ [rbox]``.

    Returns the updated ``(yhat, a_frame, rbox)`` (``yhat`` a point vector,
    ``a_frame`` a point orthogonal matrix, ``rbox`` a box in the rotating frame),
    or ``None`` if any sub-step cannot be validated (caller shrinks ``h``).
    """
    n = len(yhat)
    # (1) full interval-box hull of the current set  yhat (+) A[rbox]
    dev = _imatvec(iv, a_frame, rbox)
    ybox = [yhat[i] + dev[i] for i in range(n)]

    # (2) validated C0 flow of the CENTER point -> new center + thin remainder
    cflow = validated_step(iv, jet, list(yhat), h, mu, order=order)
    if cflow is None:
        return None
    yhat_new = [_mid(cflow[i]) for i in range(n)]  # thin-interval centre
    e_c = [cflow[i] - yhat_new[i] for i in range(n)]  # thin, centred ~0

    # (3) rigorous enclosure of D phi_h over the whole IC box, via the augmented
    #     (state + STM) validated flow started from  ybox (state) (+) I (STM).
    aug0 = list(ybox)
    for r in range(n):
        for cc in range(n):
            aug0.append(iv.mpf(1) if r == cc else iv.mpf(0))
    augf = validated_step(iv, var_jet, aug0, h, mu, order=order)
    if augf is None:
        return None
    jflat = augf[n:]
    jmat: IMat = [[jflat[r * n + cc] for cc in range(n)] for r in range(n)]

    # (4) M = J . A   (interval matrix)
    mmat = _imatmul(iv, jmat, a_frame)

    # (5) pick the new orthogonal frame: QR of mid(M)  (float, heuristic only)
    mid_m = np.array([[_midf(mmat[i][j]) for j in range(n)] for i in range(n)])
    q_float, _r_float = np.linalg.qr(mid_m)
    a_new: IMat = [[iv.mpf(float(q_float[i][j])) for j in range(n)] for i in range(n)]

    # (6) rigorous inverse of the new frame
    a_new_inv = rigorous_inverse(iv, a_new)
    if a_new_inv is None:
        return None

    # (7) re-express the deviation set in the new frame:
    #     rbox_new = (A_new^{-1} M) rbox  +  A_new^{-1} e_c
    ainv_m = _imatmul(iv, a_new_inv, mmat)
    t1 = _imatvec(iv, ainv_m, rbox)
    t2 = _imatvec(iv, a_new_inv, e_c)
    rbox_new = [t1[i] + t2[i] for i in range(n)]
    return yhat_new, a_new, rbox_new


def enclosure_box(iv: Any, yhat: list[Any], a_frame: IMat, rbox: list[Any]) -> list[Any]:
    """Interval-box hull of the parallelepiped enclosure ``yhat + a_frame @ rbox``."""
    dev = _imatvec(iv, a_frame, rbox)
    return [yhat[i] + dev[i] for i in range(len(yhat))]


def integrate_c1_qr(
    iv: Any,
    jet: Callable[..., Any],
    var_jet: Callable[..., Any],
    y0: list[Any],
    t_final: float,
    mu: float,
    *,
    n_steps: int,
    order: int = 10,
) -> dict[str, Any]:
    """Multi-step Lohner-C1/QR validated integration over ``[0, t_final]``.

    ``jet`` is the plain state jet (for the point-centre C0 flow); ``var_jet`` is
    the augmented state+STM jet (for the flow Jacobian).  ``y0`` is a list of
    (point) intervals.  Returns the same result shape as :func:`integrate` plus
    the running frame, so the enclosure box at any step is
    ``enclosure_box(iv, yhat, a_frame, rbox)``.
    """
    n = len(y0)
    yhat = [_mid(y0[i]) for i in range(n)]
    a_frame: IMat = [[iv.mpf(1) if i == j else iv.mpf(0) for j in range(n)] for i in range(n)]
    rbox = [y0[i] - yhat[i] for i in range(n)]  # thin (0 if y0 is a point)
    h = t_final / n_steps
    widths: list[float] = []
    ok = True
    for _ in range(n_steps):
        nxt = _qr_step(iv, jet, var_jet, yhat, a_frame, rbox, h, mu, order)
        if nxt is None:
            ok = False
            break
        yhat, a_frame, rbox = nxt
        box = enclosure_box(iv, yhat, a_frame, rbox)
        widths.append(max(float(c.delta.b) for c in box) / 2.0)
    return {
        "final": enclosure_box(iv, yhat, a_frame, rbox),
        "final_yhat": yhat,
        "final_frame": a_frame,
        "final_rbox": rbox,
        "validated": ok,
        "n_completed": len(widths),
        "max_halfwidths": widths,
        "final_max_halfwidth": widths[-1] if widths else None,
        "step_h": h,
        "order": order,
    }


__all__ = [
    "HAVE_MPMATH",
    "apriori_enclosure",
    "cr3bp_planar_jacobi",
    "cr3bp_planar_jet",
    "cr3bp_planar_variational_jet",
    "enclosure_box",
    "integrate",
    "integrate_c1_qr",
    "rigorous_inverse",
    "solution_series",
    "ts_add",
    "ts_add_scalar",
    "ts_const",
    "ts_mul",
    "ts_pow",
    "ts_scale",
    "ts_sub",
    "validated_step",
]
