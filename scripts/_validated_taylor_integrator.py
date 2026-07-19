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
  * It is NOT yet the full W-Z proof.  It has no QR-coordinate ("Lohner
    reframing") wrapping-effect control, so over long integrations near the
    ~2e4x/period L1/L2 saddle the enclosure width grows and eventually blows up
    (this is expected and is quantified honestly by the driver).  It has no
    C1-variational (STM) enclosure, no Poincare-section-map rigor, and no
    covering-relations / h-set machinery.  Those are the remaining stages of a
    full W-Z reproduction and are explicitly future work.  A short, low-stretch
    arc is what this dispatch rigorously certifies -- clearly a *fraction* of the
    published proof, not the whole thing.

All arithmetic is ``mpmath.iv`` at ``mp.iv.dps`` precision set by the caller.
"""

from __future__ import annotations

from typing import Any, TypeAlias

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


__all__ = [
    "HAVE_MPMATH",
    "apriori_enclosure",
    "cr3bp_planar_jacobi",
    "cr3bp_planar_jet",
    "integrate",
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
