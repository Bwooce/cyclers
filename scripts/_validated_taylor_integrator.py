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
# #670 -- Levi-Civita close-approach REGULARIZATION around the secondary.       #
#                                                                             #
# The plain CR3BP jet's  mu/r2  pseudo-potential term is genuinely singular as  #
# the trajectory approaches the secondary (Jupiter): #669 proved the rigorous   #
# Oterma-arc enclosure's real horizon is a physical close flyby (perijove r2 ~  #
# 0.004 at t ~ 0.466), where the a-priori Picard box's width overtakes the      #
# shrinking miss distance -- a mathematical singularity in the equations of     #
# motion themselves, which no coordinate REFRAMING (QR/Lohner-C1) can remove.   #
#                                                                             #
# Levi-Civita regularization (Levi-Civita 1920; Szebehely, *Theory of Orbits*   #
# 1967, Ch.3; in-corpus digest 2026-06-19-digest-szebehely-1967.md) removes the #
# 1/r singularity by construction.  Put the secondary at the origin via         #
# xi = x-(1-mu), eta = y, so r2 = sqrt(xi^2 + eta^2).  Introduce the complex    #
# square map  z = w^2  (z = xi + i*eta, w = u + i*v):                           #
#                                                                             #
#     xi = u^2 - v^2,   eta = 2 u v,   r2 = u^2 + v^2 = |w|^2                    #
#                                                                             #
# together with the fictitious-time reparametrization  dt = r2 * dtau  (so      #
# d/dt = (1/r2) d/dtau).  Working canonically (Poincare time transformation on  #
# the extended-phase-space Hamiltonian  Gamma = r2 * (H - h)  with  h = -C/2    #
# the fixed Jacobi energy), the singular  -r2 * mu/r2 = -mu  term becomes a     #
# CONSTANT -- the quadratic map's Jacobian (|dz/dw|^2 = 4 r2) vanishes at        #
# exactly the rate that cancels the 1/r2 blow-up -- so the regularized field is  #
# analytic straight through r2 = 0.  A genuine close approach becomes an         #
# ordinary, smoothly-traversable arc.  The regularized Hamiltonian is           #
#                                                                             #
#   Gamma = (Pu^2 + Pv^2)/8 + r2 (v Pu - u Pv)/2 - (1-mu)(v Pu + u Pv)/2         #
#             - (1-mu) r2 / r1 - mu - r2 h                                       #
#                                                                             #
# with  Pu, Pv  the momenta conjugate to  u, v  and  r1 = sqrt((u^2-v^2+1)^2 +   #
# (2uv)^2)  (distance to the PRIMARY, which stays ~1 near Jupiter, hence still   #
# analytic).  Hamilton's equations d/dtau = J grad Gamma give the state jet      #
# below.  A 5th quadrature coordinate  T' = r2  carries the elapsed PHYSICAL     #
# time t = integral r2 dtau, so the rigorous enclosure yields rigorous bounds    #
# on the physical time reached.  The same interval-Taylor primitives are reused. #
# --------------------------------------------------------------------------- #
def lc_secondary_from_physical(iv: Any, x: Any, y: Any, vx: Any, vy: Any, mu: float) -> list[Any]:
    """Physical planar state ``[x,y,vx,vy]`` -> regularized ``[u,v,Pu,Pv,T=0]``.

    Secondary-centred ``xi = x-(1-mu)``, ``eta = y``; ``w = sqrt(z)`` with
    ``z = xi + i eta`` (principal branch, real-axis sign carried by ``eta``);
    canonical momenta ``Pu = 2u pxi + 2v peta``, ``Pv = -2v pxi + 2u peta`` from
    the secondary-frame momenta ``pxi = vx - eta``, ``peta = vy + xi + (1-mu)``.
    ``T`` (elapsed physical time) starts at 0.
    """
    om1 = iv.mpf(1.0 - mu)
    xi = x - om1
    eta = y
    r2 = iv.sqrt(xi**2 + eta**2)
    u = iv.sqrt((r2 + xi) / 2)
    v_mag = iv.sqrt((r2 - xi) / 2)
    v = v_mag if bool(eta.b >= 0) or bool(eta.a >= 0) else -v_mag
    pxi = vx - eta
    peta = vy + xi + om1
    pu = 2 * u * pxi + 2 * v * peta
    pv = -2 * v * pxi + 2 * u * peta
    return [u, v, pu, pv, iv.mpf(0.0)]


def lc_secondary_to_physical(iv: Any, state5: list[Any], mu: float) -> list[Any]:
    """Regularized ``[u,v,Pu,Pv,(T)]`` -> physical planar state ``[x,y,vx,vy]``.

    Inverse of :func:`lc_secondary_from_physical` (``T`` ignored):
    ``xi = u^2-v^2``, ``eta = 2uv``, ``r2 = u^2+v^2``; secondary-frame momenta
    ``pxi = (u Pu - v Pv)/(2 r2)``, ``peta = (v Pu + u Pv)/(2 r2)``; then
    ``vx = pxi + eta``, ``vy = peta - xi - (1-mu)``, ``x = xi + (1-mu)``.
    """
    om1 = iv.mpf(1.0 - mu)
    u, v, pu, pv = state5[0], state5[1], state5[2], state5[3]
    xi = u * u - v * v
    eta = 2 * u * v
    r2 = u * u + v * v
    pxi = (u * pu - v * pv) / (2 * r2)
    peta = (v * pu + u * pv) / (2 * r2)
    vx = pxi + eta
    vy = peta - xi - om1
    return [xi + om1, eta, vx, vy]


def make_cr3bp_lc_secondary_jet(h: float) -> Callable[..., list[Series]]:
    """Build the Levi-Civita regularized planar-CR3BP jet for fixed energy ``h``.

    Returns a ``jet(iv, state, mu)`` compatible with the ``ts_*`` machinery, whose
    5-component state is ``[u, v, Pu, Pv, T]`` in fictitious time ``tau``.  ``h``
    is the fixed Jacobi energy ``-C/2`` (constant along the arc).  The field is
    ANALYTIC through the secondary close approach (no ``1/r2``); the only residual
    denominators are powers of ``r1 = dist-to-primary``, which stays ``~1``.
    """

    def jet(iv: Any, state: list[Series], mu: float) -> list[Series]:
        u, v, pu, pv, _t = state
        om1 = iv.mpf(1.0 - mu)
        two_h = iv.mpf(2.0 * h)
        half = iv.mpf(0.5)
        quarter = iv.mpf(0.25)

        r2 = ts_add(iv, ts_mul(iv, u, u), ts_mul(iv, v, v))  # u^2 + v^2
        amom = ts_sub(iv, ts_mul(iv, v, pu), ts_mul(iv, u, pv))  # v Pu - u Pv
        # P = u^2 - v^2 + 1,  Q = 2 u v  ->  r1^2 = P^2 + Q^2
        pp = ts_add_scalar(iv, ts_sub(iv, ts_mul(iv, u, u), ts_mul(iv, v, v)), iv.mpf(1.0))
        qq = ts_scale(iv, ts_mul(iv, u, v), iv.mpf(2.0))
        r1sq = ts_add(iv, ts_mul(iv, pp, pp), ts_mul(iv, qq, qq))
        r1i = ts_pow(iv, r1sq, -0.5)  # 1 / r1
        r1i3 = ts_pow(iv, r1sq, -1.5)  # 1 / r1^3

        # u' = Pu/4 + (r2 - om1) v / 2
        up = ts_add(
            iv,
            ts_scale(iv, pu, quarter),
            ts_scale(iv, ts_mul(iv, ts_add_scalar(iv, r2, -om1), v), half),
        )
        # v' = Pv/4 - (r2 + om1) u / 2
        vp = ts_sub(
            iv,
            ts_scale(iv, pv, quarter),
            ts_scale(iv, ts_mul(iv, ts_add_scalar(iv, r2, om1), u), half),
        )
        # Pu' = -u A + (r2 + om1) Pv/2 + 2 om1 u / r1 - om1 r2 (2 u P + 2 v Q)/r1^3 + 2 h u
        inner1 = ts_scale(iv, ts_add(iv, ts_mul(iv, u, pp), ts_mul(iv, v, qq)), iv.mpf(2.0))
        pup = ts_scale(iv, ts_mul(iv, u, amom), iv.mpf(-1.0))
        pup = ts_add(iv, pup, ts_scale(iv, ts_mul(iv, ts_add_scalar(iv, r2, om1), pv), half))
        pup = ts_add(iv, pup, ts_scale(iv, ts_mul(iv, u, r1i), 2 * om1))
        pup = ts_sub(iv, pup, ts_scale(iv, ts_mul(iv, ts_mul(iv, r2, inner1), r1i3), om1))
        pup = ts_add(iv, pup, ts_scale(iv, u, two_h))
        # Pv' = -v A + (om1 - r2) Pu/2 + 2 om1 v / r1 - om1 r2 (-2 v P + 2 u Q)/r1^3 + 2 h v
        inner2 = ts_scale(iv, ts_sub(iv, ts_mul(iv, u, qq), ts_mul(iv, v, pp)), iv.mpf(2.0))
        om1_m_r2 = ts_add_scalar(iv, ts_scale(iv, r2, iv.mpf(-1.0)), om1)
        pvp = ts_scale(iv, ts_mul(iv, v, amom), iv.mpf(-1.0))
        pvp = ts_add(iv, pvp, ts_scale(iv, ts_mul(iv, om1_m_r2, pu), half))
        pvp = ts_add(iv, pvp, ts_scale(iv, ts_mul(iv, v, r1i), 2 * om1))
        pvp = ts_sub(iv, pvp, ts_scale(iv, ts_mul(iv, ts_mul(iv, r2, inner2), r1i3), om1))
        pvp = ts_add(iv, pvp, ts_scale(iv, v, two_h))
        # T' = r2  (dt = r2 dtau) -- elapsed physical time as a quadrature
        return [up, vp, pup, pvp, r2]

    return jet


def make_cr3bp_lc_secondary_variational_jet(h: float) -> Callable[..., list[Series]]:
    """Build the C1 (state + STM) Levi-Civita regularized jet for fixed energy ``h``.

    This is ``#671``'s Stage-5 addition: the REGULARIZED-coordinate analogue of
    :func:`cr3bp_planar_variational_jet` (which carries the STM for the *physical*
    CR3BP).  ``#670`` regularized the secondary close approach but integrated it
    with the plain (non-variational) jet, so ``#669``'s QR/Lohner-C1 wrapping fix
    could not be applied in the regularized frame.  This jet supplies the missing
    piece: the variational flow ``V' = Dg(w) V`` of the Levi-Civita field
    ``g(w)``, ``w = [u, v, Pu, Pv, T]``, so :func:`integrate_c1_qr` can reframe the
    enclosure in the regularized coordinates' own local stretching directions.

    The augmented state is ``[u, v, Pu, Pv, T, V00, V01, ..., V44]`` -- the 5-state
    plus a row-major 5x5 STM (25 series).  ``Dg`` is the analytic Jacobian of the
    exact same field :func:`make_cr3bp_lc_secondary_jet` produces (re-derived and
    validated against finite differences of that field, and against the closed-form
    harmonic STM of the radial-collision anchor).  The Jacobian is worked out from

        u'  = Pu/4 + (r2 - w1) v/2,                w1 := 1-mu
        v'  = Pv/4 - (r2 + w1) u/2,
        Pu' = -u A + (r2 + w1) Pv/2 + 2 w1 u/r1 - w1 r2 I1/r1^3 + 2h u,
        Pv' = -v A + (w1 - r2) Pu/2 + 2 w1 v/r1 - w1 r2 I2/r1^3 + 2h v,
        T'  = r2,

    with ``A = v Pu - u Pv``, ``P = u^2-v^2+1``, ``Q = 2uv``, ``r1^2 = P^2+Q^2``,
    ``I1 = 2(uP+vQ)``, ``I2 = 2(uQ-vP)``, and the key identities ``d(r1^2)/du =
    2 I1``, ``d(r1^2)/dv = 2 I2``.  The resulting 4x4 Jacobian block is (Hamiltonian
    symmetry: the (Pu',Pv') cross terms coincide, and ``du'/dv = -dPv'/dPu`` etc.,
    all verified below by construction)::

        Df[u']  = [ uv,              (r2-w1)/2 + v^2,  1/4,               0            ]
        Df[v']  = [ -(r2+w1)/2 - u^2, -uv,             0,                 1/4          ]
        Df[Pu'] = [ J20,             J21,              -uv,               u^2+(r2+w1)/2 ]
        Df[Pv'] = [ J21,             J31,              (w1-r2)/2 - v^2,   uv           ]
        Df[T']  = [ 2u,              2v,               0,                 0            ]

    (the 5th, T, column of ``Df`` is identically zero -- nothing depends on ``T``),
    with
        J20 = -A + 2u Pv + 2 w1/r1 - 4 w1 u I1/r1^3 - w1 r2 (2P+4r2)/r1^3
                 + 3 w1 r2 I1^2/r1^5 + 2h,
        J21 = -u Pu + v Pv - 2 w1 u I2/r1^3 - 2 w1 v I1/r1^3 - 2 w1 r2 Q/r1^3
                 + 3 w1 r2 I1 I2/r1^5,
        J31 = -A - 2v Pu + 2 w1/r1 - 4 w1 v I2/r1^3 - w1 r2 (4r2-2P)/r1^3
                 + 3 w1 r2 I2^2/r1^5 + 2h.
    """

    def jet(iv: Any, aug: list[Series], mu: float) -> list[Series]:
        u, v, pu, pv, _t = aug[0], aug[1], aug[2], aug[3], aug[4]
        V = aug[5:30]  # 25 series, row-major 5x5
        om1 = iv.mpf(1.0 - mu)
        two_h = iv.mpf(2.0 * h)
        half = iv.mpf(0.5)
        quarter = iv.mpf(0.25)

        # --- shared field quantities (identical to make_cr3bp_lc_secondary_jet) ---
        r2 = ts_add(iv, ts_mul(iv, u, u), ts_mul(iv, v, v))  # u^2 + v^2
        amom = ts_sub(iv, ts_mul(iv, v, pu), ts_mul(iv, u, pv))  # A = v Pu - u Pv
        pp = ts_add_scalar(iv, ts_sub(iv, ts_mul(iv, u, u), ts_mul(iv, v, v)), iv.mpf(1.0))  # P
        qq = ts_scale(iv, ts_mul(iv, u, v), iv.mpf(2.0))  # Q = 2uv
        r1sq = ts_add(iv, ts_mul(iv, pp, pp), ts_mul(iv, qq, qq))
        r1i = ts_pow(iv, r1sq, -0.5)  # 1/r1
        r1i3 = ts_pow(iv, r1sq, -1.5)  # 1/r1^3
        r1i5 = ts_pow(iv, r1sq, -2.5)  # 1/r1^5 (STM only)
        inner1 = ts_scale(iv, ts_add(iv, ts_mul(iv, u, pp), ts_mul(iv, v, qq)), iv.mpf(2.0))  # I1
        inner2 = ts_scale(iv, ts_sub(iv, ts_mul(iv, u, qq), ts_mul(iv, v, pp)), iv.mpf(2.0))  # I2

        # --- state derivatives (identical to make_cr3bp_lc_secondary_jet) ---------
        up = ts_add(
            iv,
            ts_scale(iv, pu, quarter),
            ts_scale(iv, ts_mul(iv, ts_add_scalar(iv, r2, -om1), v), half),
        )
        vp = ts_sub(
            iv,
            ts_scale(iv, pv, quarter),
            ts_scale(iv, ts_mul(iv, ts_add_scalar(iv, r2, om1), u), half),
        )
        i1_pot = ts_scale(iv, ts_add(iv, ts_mul(iv, u, pp), ts_mul(iv, v, qq)), iv.mpf(2.0))
        pup = ts_scale(iv, ts_mul(iv, u, amom), iv.mpf(-1.0))
        pup = ts_add(iv, pup, ts_scale(iv, ts_mul(iv, ts_add_scalar(iv, r2, om1), pv), half))
        pup = ts_add(iv, pup, ts_scale(iv, ts_mul(iv, u, r1i), 2 * om1))
        pup = ts_sub(iv, pup, ts_scale(iv, ts_mul(iv, ts_mul(iv, r2, i1_pot), r1i3), om1))
        pup = ts_add(iv, pup, ts_scale(iv, u, two_h))
        i2_pot = ts_scale(iv, ts_sub(iv, ts_mul(iv, u, qq), ts_mul(iv, v, pp)), iv.mpf(2.0))
        om1_m_r2 = ts_add_scalar(iv, ts_scale(iv, r2, iv.mpf(-1.0)), om1)
        pvp = ts_scale(iv, ts_mul(iv, v, amom), iv.mpf(-1.0))
        pvp = ts_add(iv, pvp, ts_scale(iv, ts_mul(iv, om1_m_r2, pu), half))
        pvp = ts_add(iv, pvp, ts_scale(iv, ts_mul(iv, v, r1i), 2 * om1))
        pvp = ts_sub(iv, pvp, ts_scale(iv, ts_mul(iv, ts_mul(iv, r2, i2_pot), r1i3), om1))
        pvp = ts_add(iv, pvp, ts_scale(iv, v, two_h))
        tp = r2

        # --- Jacobian Df (4x4 on [u,v,Pu,Pv]; row T and column T handled apart) ---
        uv = ts_mul(iv, u, v)
        u2 = ts_mul(iv, u, u)
        v2 = ts_mul(iv, v, v)
        neg_uv = ts_scale(iv, uv, iv.mpf(-1.0))
        # row u'
        d_uu = uv
        d_uv = ts_add(iv, ts_scale(iv, ts_add_scalar(iv, r2, -om1), half), v2)  # (r2-w1)/2 + v^2
        # row v'
        d_vu = ts_sub(iv, ts_scale(iv, ts_add_scalar(iv, r2, om1), -half), u2)  # -(r2+w1)/2 - u^2
        d_vv = neg_uv
        # row Pu'  (J20, J21, -uv, u^2 + (r2+w1)/2)
        j20 = ts_scale(iv, amom, iv.mpf(-1.0))
        j20 = ts_add(iv, j20, ts_scale(iv, ts_mul(iv, u, pv), iv.mpf(2.0)))
        j20 = ts_add(iv, j20, ts_scale(iv, r1i, 2 * om1))
        j20 = ts_sub(iv, j20, ts_scale(iv, ts_mul(iv, ts_mul(iv, u, inner1), r1i3), 4 * om1))
        twoP_4r2 = ts_add(iv, ts_scale(iv, pp, iv.mpf(2.0)), ts_scale(iv, r2, iv.mpf(4.0)))
        j20 = ts_sub(iv, j20, ts_scale(iv, ts_mul(iv, ts_mul(iv, r2, twoP_4r2), r1i3), om1))
        j20 = ts_add(
            iv,
            j20,
            ts_scale(iv, ts_mul(iv, ts_mul(iv, r2, ts_mul(iv, inner1, inner1)), r1i5), 3 * om1),
        )
        j20 = ts_add_scalar(iv, j20, two_h)
        # J21 = dPu'/dv = dPv'/du  (Hamiltonian cross term)
        j21 = ts_scale(iv, ts_mul(iv, u, pu), iv.mpf(-1.0))
        j21 = ts_add(iv, j21, ts_mul(iv, v, pv))
        j21 = ts_sub(iv, j21, ts_scale(iv, ts_mul(iv, ts_mul(iv, u, inner2), r1i3), 2 * om1))
        j21 = ts_sub(iv, j21, ts_scale(iv, ts_mul(iv, ts_mul(iv, v, inner1), r1i3), 2 * om1))
        j21 = ts_sub(iv, j21, ts_scale(iv, ts_mul(iv, ts_mul(iv, r2, qq), r1i3), 2 * om1))
        j21 = ts_add(
            iv,
            j21,
            ts_scale(iv, ts_mul(iv, ts_mul(iv, r2, ts_mul(iv, inner1, inner2)), r1i5), 3 * om1),
        )
        d_2Pu = neg_uv
        d_2Pv = ts_add(iv, u2, ts_scale(iv, ts_add_scalar(iv, r2, om1), half))  # u^2+(r2+w1)/2
        # row Pv'  (J21, J31, (w1-r2)/2 - v^2, uv)
        j31 = ts_scale(iv, amom, iv.mpf(-1.0))
        j31 = ts_sub(iv, j31, ts_scale(iv, ts_mul(iv, v, pu), iv.mpf(2.0)))
        j31 = ts_add(iv, j31, ts_scale(iv, r1i, 2 * om1))
        j31 = ts_sub(iv, j31, ts_scale(iv, ts_mul(iv, ts_mul(iv, v, inner2), r1i3), 4 * om1))
        fourr2_2P = ts_sub(iv, ts_scale(iv, r2, iv.mpf(4.0)), ts_scale(iv, pp, iv.mpf(2.0)))
        j31 = ts_sub(iv, j31, ts_scale(iv, ts_mul(iv, ts_mul(iv, r2, fourr2_2P), r1i3), om1))
        j31 = ts_add(
            iv,
            j31,
            ts_scale(iv, ts_mul(iv, ts_mul(iv, r2, ts_mul(iv, inner2, inner2)), r1i5), 3 * om1),
        )
        j31 = ts_add_scalar(iv, j31, two_h)  # + 2h  (from the +2h v term of Pv')
        d_3Pu = ts_sub(iv, ts_scale(iv, om1_m_r2, half), v2)  # (w1-r2)/2 - v^2
        d_3Pv = uv

        # Df rows as [c_u, c_v, c_Pu, c_Pv] series lists; column T is zero.
        zero = ts_const(iv, iv.mpf(0), len(u) - 1)
        quarter_s = ts_const(iv, quarter, len(u) - 1)
        two_u = ts_scale(iv, u, iv.mpf(2.0))
        two_v = ts_scale(iv, v, iv.mpf(2.0))
        df = [
            [d_uu, d_uv, quarter_s, zero, zero],
            [d_vu, d_vv, zero, quarter_s, zero],
            [j20, j21, d_2Pu, d_2Pv, zero],
            [j21, j31, d_3Pu, d_3Pv, zero],
            [two_u, two_v, zero, zero, zero],
        ]

        # --- STM block  dV = Df . V  (5x5, row-major) ---
        # Precompute the non-zero (row, col) index pattern of Df (``zero`` entries
        # -- the T column, and the structural zeros of rows u'/v'/T' -- are skipped;
        # identity test is valid since every zero slot reuses the same object).
        nz = [[kk for kk in range(5) if df[rrow][kk] is not zero] for rrow in range(5)]
        dV: list[Series] = [zero for _ in range(25)]
        for c in range(5):  # column of V
            col = [V[0 * 5 + c], V[1 * 5 + c], V[2 * 5 + c], V[3 * 5 + c], V[4 * 5 + c]]
            for rrow in range(5):
                ks = nz[rrow]
                acc = ts_mul(iv, df[rrow][ks[0]], col[ks[0]])
                for kk in ks[1:]:
                    acc = ts_add(iv, acc, ts_mul(iv, df[rrow][kk], col[kk]))
                dV[rrow * 5 + c] = acc
        return [up, vp, pup, pvp, tp, *dV]

    return jet


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
    bmat: IMat,
) -> tuple[list[Any], IMat, list[Any], IMat] | None:
    """One Lohner-C1/QR step of the enclosure ``yhat + a_frame @ [rbox]``.

    Returns the updated ``(yhat, a_frame, rbox, bmat)`` (``yhat`` a point vector,
    ``a_frame`` a point orthogonal matrix, ``rbox`` a box in the rotating frame),
    or ``None`` if any sub-step cannot be validated (caller shrinks ``h``).

    ``bmat`` is the accumulated frame-relative state-transition matrix: if the running
    frame product satisfies ``Phi_{0->k} = a_frame @ bmat`` on entry, then it does
    again on exit with the returned ``bmat``.  This is the standard Lohner-C1 route to
    a rigorous *and* tight accumulated flow Jacobian: the per-step factor
    ``a_new^{-1} M`` (already formed for the deviation update) is the well-conditioned
    R-like factor of the QR, so the interval product of these factors stays tight even
    as the true STM magnitude grows exponentially -- the same reframing that keeps the
    enclosure box from wrapping keeps ``bmat`` from wrapping.  (Pass an identity
    ``bmat`` at step 0; ignore the returned value if the accumulated STM is not needed.)
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
    # (8) accumulated frame-relative STM:  Phi = a_new @ bmat_new,  bmat_new =
    #     (a_new^{-1} M) @ bmat  (same well-conditioned factor as the deviation map)
    bmat_new = _imatmul(iv, ainv_m, bmat)
    return yhat_new, a_new, rbox_new, bmat_new


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
    bmat: IMat = [[iv.mpf(1) if i == j else iv.mpf(0) for j in range(n)] for i in range(n)]
    h = t_final / n_steps
    widths: list[float] = []
    ok = True
    for _ in range(n_steps):
        nxt = _qr_step(iv, jet, var_jet, yhat, a_frame, rbox, h, mu, order, bmat)
        if nxt is None:
            ok = False
            break
        yhat, a_frame, rbox, bmat = nxt
        box = enclosure_box(iv, yhat, a_frame, rbox)
        widths.append(max(float(c.delta.b) for c in box) / 2.0)
    return {
        "final": enclosure_box(iv, yhat, a_frame, rbox),
        "final_yhat": yhat,
        "final_frame": a_frame,
        "final_rbox": rbox,
        "final_stm": _imatmul(iv, a_frame, bmat),
        "validated": ok,
        "n_completed": len(widths),
        "max_halfwidths": widths,
        "final_max_halfwidth": widths[-1] if widths else None,
        "step_h": h,
        "order": order,
    }


# --------------------------------------------------------------------------- #
# #672 -- Rigorous Poincare SECTION-CROSSING MAP with its own rigorous Jacobian.#
#                                                                             #
# The W-Z existence proof is built on covering relations between h-sets placed  #
# along the Oterma orbit, verified through the Poincare SECTION MAP P (not the   #
# raw flow).  Their section is Theta = {y = 0}, parameterized by (x, xdot); in   #
# the Levi-Civita regularized coordinates carried by #670/#671 the physical      #
# ordinate is  y = eta = 2 u v, so the section condition is  sigma(w) = 2 u v = 0#
# (a coordinate section in physical space; the same code handles ANY smooth      #
# scalar section via the (sigma_val, sigma_grad) descriptor below).             #
#                                                                             #
# This builds the necessary PRECURSOR to h-sets: given a rigorously-enclosed     #
# initial condition, (1) rigorously ISOLATE the fictitious-time interval [tau*]  #
# of the first section crossing -- an interval-Newton root of sigma(state(s))=0  #
# on the crossing step's rigorous Taylor MODEL (poly coeffs from the IC box +    #
# an order-(p+1) Lagrange-remainder coefficient bounded over the whole-step      #
# a-priori enclosure), NOT a floating-point event trigger; the same step         #
# certifies transversality (0 not in dsigma/ds over the step); (2) evaluate the  #
# enclosed state AT [tau*]; and (3) form the rigorous section-map Jacobian        #
#                                                                             #
#   DP = Dphi_tau*  -  f(P) . (Dsigma(P) . Dphi_tau*) / (Dsigma(P) . f(P))       #
#                                                                             #
# (the flow's variational matrix Dphi_tau* PLUS the implicit-function correction #
# for the crossing time's dependence on the IC -- derived by differentiating the #
# crossing condition sigma(phi_{tau*(w)}(w)) = 0, chain rule + implicit function #
# theorem; standard Poincare-map-derivative material, e.g. Simo's / Guckenheimer #
# & Holmes; get it exactly, do not approximate).  Dphi_tau* is the accumulated   #
# STM Phi_{0->tau*} obtained tightly from the QR stepper's ``bmat`` (Phi=A@B)     #
# composed with the crossing step's partial-step STM Taylor model.              #
#                                                                             #
# Validated to full interval rigor against a closed-form positive control (the   #
# 2D harmonic oscillator, section {y1=0}: crossing time, post-map state, AND the #
# full DP are hand-computable and independently verified).  Building the actual  #
# h-sets and proving covering relations remains the NEXT stage, not this one.    #
# --------------------------------------------------------------------------- #
def eval_series_horner(iv: Any, series: Series, s: Any) -> Any:
    """Rigorously enclose ``sum_k series[k] s^k`` (Horner), ``s`` an iv interval."""
    acc = series[-1]
    for k in range(len(series) - 2, -1, -1):
        acc = acc * s + series[k]
    return acc


def deriv_series(iv: Any, series: Series) -> Series:
    """Termwise derivative series:  d/ds sum_k c_k s^k = sum_k k c_k s^{k-1}."""
    return [series[k] * iv.mpf(k) for k in range(1, len(series))]


def flow_taylor_model(
    iv: Any,
    jet: Any,
    ybox: list[Any],
    h: float,
    mu: float,
    order: int,
) -> list[Series] | None:
    """Rigorous Taylor MODEL of the flow over ``[0, h]`` from the IC box ``ybox``.

    Returns one :data:`Series` per state component, length ``order + 2``:
    coefficients ``0..order`` are the interval Taylor coefficients of the solution
    started anywhere in ``ybox`` (``solution_series(ybox)``); coefficient
    ``order + 1`` is the Lagrange-remainder coefficient bounded over the whole-step
    a-priori enclosure ``W`` (``solution_series(W)[order+1]``).  Then for every
    ``ic in ybox`` and every ``s in [0, h]``,

        y_i(t0 + s; ic)  in  sum_{k=0}^{order} c[i][k] s^k  +  c[i][order+1] s^{order+1}

    i.e. :func:`eval_series_horner` of the returned model at any ``s in [0, h]``
    rigorously encloses the true state -- the object the crossing isolation and the
    partial-step STM both evaluate.  ``None`` if the a-priori step cannot validate.
    """
    W = apriori_enclosure(iv, jet, ybox, h, mu)
    if W is None:
        return None
    poly = solution_series(iv, jet, ybox, order, mu)
    rem = solution_series(iv, jet, W, order + 1, mu)
    return [[*poly[i], rem[i][order + 1]] for i in range(len(ybox))]


def isolate_section_crossing(
    iv: Any,
    state_tm: list[Series],
    sigma_val: Callable[..., Any],
    sigma_grad: Callable[..., Any],
    h: float,
    *,
    iters: int = 60,
) -> Any | None:
    """Interval-Newton isolation of the section crossing on a step's Taylor model.

    ``state_tm`` is a :func:`flow_taylor_model` (list of state Series over ``[0,h]``).
    ``sigma_val(iv, statevec) -> interval`` evaluates the scalar section function on
    a state VECTOR; ``sigma_grad(iv, statevec) -> list`` its gradient row.  Returns a
    TIGHT interval ``[s*]`` (subset of ``[0, h]``) proven to contain a unique,
    transversal crossing ``sigma(state(s*)) = 0``, or ``None`` if the step does not
    rigorously bracket exactly one transversal crossing (no strict sign change over
    ``[0, h]``, or ``0`` in ``dsigma/ds`` over the step -- i.e. not certified
    monotone / not transversal).

    ``dsigma/ds`` is formed rigorously as ``grad sigma(state(s)) . state'(s)`` with
    ``state'`` the termwise derivative of the Taylor model -- so ``0 not in`` it over
    ``[0, h]`` is a genuine transversality certificate (``Dsigma . f != 0``), exactly
    the non-degeneracy the section-map Jacobian formula requires.
    """
    dstate = [deriv_series(iv, comp) for comp in state_tm]

    def sig_at(s: Any) -> Any:
        st = [eval_series_horner(iv, c, s) for c in state_tm]
        return sigma_val(iv, st)

    def dsig_at(s: Any) -> Any:
        st = [eval_series_horner(iv, c, s) for c in state_tm]
        dst = [eval_series_horner(iv, c, s) for c in dstate]
        g = sigma_grad(iv, st)
        acc = iv.mpf(0)
        for k in range(len(st)):
            acc = acc + g[k] * dst[k]
        return acc

    s0 = sig_at(iv.mpf([0.0, 0.0]))
    sh = sig_at(iv.mpf([h, h]))
    up = bool(s0.b < 0) and bool(sh.a > 0)
    down = bool(s0.a > 0) and bool(sh.b < 0)
    if not (up or down):
        return None  # no rigorous strict sign change across the step
    dfull = dsig_at(iv.mpf([0.0, h]))
    if not (bool(dfull.a > 0) or bool(dfull.b < 0)):
        return None  # 0 in dsigma/ds over the step: not certified monotone/transversal

    bracket = iv.mpf([0.0, h])
    for _ in range(iters):
        sc = (bracket.a + bracket.b) / 2
        sc_iv = iv.mpf([sc, sc])
        num = sig_at(sc_iv)
        den = dsig_at(bracket)  # 0 not in den (certified above)
        nn = sc_iv - num / den
        lo = nn.a if bool(nn.a > bracket.a) else bracket.a
        hi = nn.b if bool(nn.b < bracket.b) else bracket.b
        if bool(lo > hi):
            return None  # Newton image disjoint -> no root (contradiction-safe)
        width_old = bracket.b - bracket.a
        bracket = iv.mpf([lo, hi])
        if not bool(bracket.b - bracket.a < width_old):
            break  # converged: interval no longer shrinking
    return bracket


def section_map_jacobian(
    iv: Any,
    stm: IMat,
    fvec: list[Any],
    dsig: list[Any],
) -> tuple[IMat, Any, list[Any]]:
    """Rigorous Poincare section-map Jacobian from the flow Jacobian at the crossing.

    ``stm`` is ``Dphi_tau* = Phi_{0->tau*}`` (the accumulated flow STM to the
    crossing), ``fvec = f(P)`` the vector field at the crossing state, ``dsig =
    Dsigma(P)`` the section gradient there.  Returns ``(DP, Dsigma.f, Dsigma.Phi)``
    with

        DP = Phi  -  f (Dsigma . Phi) / (Dsigma . f)

    the standard implicit-function correction for the crossing time's dependence on
    the initial condition (differentiate ``sigma(phi_{tau*(w)}(w)) = 0``:
    ``Dsigma (f Dtau* + Phi) = 0`` gives ``Dtau* = -(Dsigma Phi)/(Dsigma f)``, and
    ``DP = f Dtau* + Phi``).  The caller must have certified ``Dsigma . f != 0``
    (transversality); it is returned so the caller can assert non-vacuous division.
    """
    n = len(fvec)
    dsf = iv.mpf(0)
    for k in range(n):
        dsf = dsf + dsig[k] * fvec[k]
    dsphi: list[Any] = []
    for j in range(n):
        acc = iv.mpf(0)
        for k in range(n):
            acc = acc + dsig[k] * stm[k][j]
        dsphi.append(acc)
    dp: IMat = [[stm[i][j] - fvec[i] * dsphi[j] / dsf for j in range(n)] for i in range(n)]
    return dp, dsf, dsphi


def rigorous_section_map(
    iv: Any,
    jet: Callable[..., Any],
    var_jet: Callable[..., Any],
    y0: list[Any],
    mu: float,
    *,
    sigma_val: Callable[..., Any],
    sigma_grad: Callable[..., Any],
    tau_max: float,
    n_steps: int,
    order: int = 10,
    tau_min: float = 0.0,
) -> dict[str, Any]:
    """Rigorously isolate the first section crossing of the flow from ``y0`` and
    return the enclosed crossing state, crossing time, and section-map Jacobian.

    Integrates the enclosure with the Lohner-C1/QR stepper (so both the enclosure
    box and the accumulated flow Jacobian ``Phi = A @ bmat`` stay tight through
    close approaches / stretching), monitoring the section value ``sigma`` on each
    step endpoint.  On the first step whose endpoints straddle ``sigma = 0`` (with
    ``tau >= tau_min``, so the trivial start-on-section crossing at ``tau = 0`` is
    skipped), it rebuilds the step's rigorous Taylor model from the (tight) step-start
    enclosure box, isolates ``[s*]`` by interval-Newton, evaluates the crossing state
    and the partial-step STM at ``[s*]``, composes the full ``Phi_{0->tau*}``, and
    forms ``DP``.  Returns a dict; ``found`` is ``True`` only if a crossing was
    rigorously isolated (with a transversality certificate).
    """
    n = len(y0)
    yhat = [_mid(y0[i]) for i in range(n)]
    a_frame: IMat = [[iv.mpf(1) if i == j else iv.mpf(0) for j in range(n)] for i in range(n)]
    rbox = [y0[i] - yhat[i] for i in range(n)]
    bmat: IMat = [[iv.mpf(1) if i == j else iv.mpf(0) for j in range(n)] for i in range(n)]
    h = tau_max / n_steps

    def sigma_box() -> Any:
        return sigma_val(iv, enclosure_box(iv, yhat, a_frame, rbox))

    prev_sign = 0  # -1 / +1 once sigma is rigorously nonzero
    result: dict[str, Any] = {
        "found": False,
        "validated_to": 0.0,
        "n_steps_taken": 0,
        "reason": "no crossing isolated in [tau_min, tau_max]",
    }
    for step in range(n_steps):
        tau0 = step * h
        # snapshot the step-START enclosure (needed if THIS step brackets the crossing)
        start = (list(yhat), [row[:] for row in a_frame], [row[:] for row in bmat], list(rbox))
        nxt = _qr_step(iv, jet, var_jet, yhat, a_frame, rbox, h, mu, order, bmat)
        if nxt is None:
            result["reason"] = f"enclosure wall at tau~{tau0:.4f} (step {step}) before any crossing"
            result["validated_to"] = tau0
            result["n_steps_taken"] = step
            return result
        yhat, a_frame, rbox, bmat = nxt
        result["validated_to"] = tau0 + h
        result["n_steps_taken"] = step + 1

        sig = sigma_box()
        cur_sign = 1 if bool(sig.a > 0) else (-1 if bool(sig.b < 0) else 0)
        if tau0 + h < tau_min or cur_sign == 0:
            if cur_sign != 0:
                prev_sign = cur_sign
            continue
        if prev_sign != 0 and cur_sign != prev_sign:
            # crossing is inside THIS step: [tau0, tau0+h]
            yhat_s, a_s, b_s, rbox_s = start
            ybox_start = enclosure_box(iv, yhat_s, a_s, rbox_s)
            phi_start = _imatmul(iv, a_s, b_s)  # Phi_{0->tau0}
            # augmented Taylor model over the step: state (n) + STM (n*n), V(0)=I
            aug_ic = list(ybox_start)
            for r in range(n):
                for c in range(n):
                    aug_ic.append(iv.mpf(1) if r == c else iv.mpf(0))
            aug_tm = flow_taylor_model(iv, var_jet, aug_ic, h, mu, order)
            if aug_tm is None:
                result["reason"] = f"crossing-step Taylor model failed at tau~{tau0:.4f}"
                return result
            state_tm = aug_tm[:n]
            s_star = isolate_section_crossing(
                iv, state_tm, sigma_val, sigma_grad, h, iters=order + 40
            )
            if s_star is None:
                # endpoints flipped but this step's model did not certify a single
                # transversal crossing (grazing / multiple / boundary): report honestly
                result["reason"] = (
                    f"endpoints straddle sigma=0 near tau~{tau0:.4f} but interval-Newton "
                    "could not certify a unique transversal crossing on the step model"
                )
                return result
            # crossing state (regularized) and partial-step STM at [s*]
            cross_state = [eval_series_horner(iv, c, s_star) for c in state_tm]
            phi_step = [
                [eval_series_horner(iv, aug_tm[n + r * n + c], s_star) for c in range(n)]
                for r in range(n)
            ]
            phi_full = _imatmul(iv, phi_step, phi_start)  # Phi_{0->tau*}
            fvec = _f_over_box(iv, jet, cross_state, mu)
            dsig = sigma_grad(iv, cross_state)
            dp, dsf, _dsphi = section_map_jacobian(iv, phi_full, fvec, dsig)
            tau_star = iv.mpf([tau0, tau0]) + s_star
            result.update(
                {
                    "found": True,
                    "reason": "rigorously isolated a unique transversal section crossing",
                    "crossing_step": step,
                    "tau_star": tau_star,
                    "s_star_local": s_star,
                    "crossing_state": cross_state,
                    "stm_full": phi_full,
                    "section_jacobian": dp,
                    "dsigma_dot_f": dsf,
                    "vector_field_at_crossing": fvec,
                    "transversal": bool(dsf.a > 0) or bool(dsf.b < 0),
                }
            )
            return result
        prev_sign = cur_sign
    return result


# --------------------------------------------------------------------------- #
# #673 -- h-sets and the RIGOROUS 2D COVERING RELATION (the existence-proof step).#
#                                                                             #
# W-Z prove existence of a connecting orbit with a chain of *covering relations* #
# between *h-sets* placed along the orbit, verified through the Poincare section  #
# map P.  This is the first stage that attempts that actual proof step.  We work  #
# in the section's own 2D coordinates (physical (x, xdot) on {y=0}), where every  #
# h-set has one unstable and one stable direction (u=s=1) -- the setting of every #
# W-Z Oterma h-set.                                                              #
#                                                                             #
# h-set (Zgliczynski-Gidea 2004, "Covering relations for multidimensional        #
# dynamical systems", J. Diff. Eq. 202; and W-Z Part I): an h-set N is a box with #
# a distinguished product structure.  Here N = { c + a*uvec + b*svec : a,b in     #
# [-1,1] }, with ``uvec`` the UNSTABLE (exit) direction and ``svec`` the STABLE    #
# (entry) direction.  Its coordinate chart c_N maps N onto the unit square        #
# [-1,1]^2 = B_u x B_s (a = unstable coord, b = stable coord).  The distinguished  #
# faces are                                                                       #
#     N^-  (exit set / unstable boundary)  = { a = +-1 }   (left/right edges),     #
#     N^+  (entry set / stable  boundary)  = { b = +-1 }   (top/bottom edges).     #
#                                                                             #
# Covering relation  N ==f==> M  (both h-sets with u=s=1, f = section map): the    #
# image of N must "stretch all the way across" M in the unstable direction while   #
# staying "pinched" inside M in the stable direction, in a topological (degree)    #
# sense -- NOT merely "the boxes overlap".  For 2D box h-sets the ZG topological    #
# definition is implied by these two rigorously-checkable conditions, expressed    #
# in M's OWN coordinates f_c = c_M o f o c_N^{-1} (Zgliczynski-Gidea Thm 4/16; the  #
# planar "correctly aligned" case used throughout the W-Z CAP verifications):      #
#                                                                             #
#   (S) STABLE containment: pi_s f_c(N) subset (-1, 1) -- the WHOLE image of N is  #
#       strictly inside M in the stable coordinate (pinched, no exit top/bottom);  #
#   (U) UNSTABLE exit, opposite sides: the two unstable edges map strictly beyond  #
#       opposite unstable ends of M --                                             #
#         (U+)  pi_u f_c(N^-_left) < -1  AND  pi_u f_c(N^-_right) > +1,   OR       #
#         (U-)  pi_u f_c(N^-_left) > +1  AND  pi_u f_c(N^-_right) < -1.            #
#                                                                             #
#   (S)+(U) => N f-covers M: the image is connected, its stable part lies in M's   #
#   stable strip, and its unstable extent runs strictly past both unstable ends    #
#   of M on opposite sides, so by the intermediate-value / degree argument the      #
#   image crosses M with nonzero degree.  This is a TOPOLOGICAL statement, not a    #
#   size comparison -- interval OVERLAP alone proves nothing and is deliberately    #
#   NOT what is checked here.                                                       #
#                                                                             #
# The checker takes the map images as already-computed rigorous ENCLOSURES (boxes  #
# in the ambient (x, xdot) section coords) of  f(N^-_left), f(N^-_right), f(N):     #
# an interval OVER-approximation is on the SOUND side of every condition -- for     #
# (S) a superset still inside (-1,1) proves the true image is; for (U) a superset   #
# still strictly beyond an end proves the true edge is.  So a verdict of "covers"   #
# from over-approximated images is a genuine rigorous certificate.                  #
# --------------------------------------------------------------------------- #
def hset_chart_inverse(iv: Any, uvec: list[Any], svec: list[Any]) -> IMat | None:
    """Rigorous inverse of an h-set's 2x2 frame ``B = [uvec | svec]`` (columns).

    ``c_N^{-1}(a,b) = center + B (a,b)``, so ``c_N(p) = B^{-1}(p - center)``; this
    returns the rigorous interval enclosure of ``B^{-1}`` (via :func:`rigorous_inverse`).
    ``uvec``/``svec`` are the unstable/stable direction column vectors (thin
    intervals or reals).  ``None`` if the frame is degenerate (collinear directions).
    """
    b_mat: IMat = [[uvec[0], svec[0]], [uvec[1], svec[1]]]
    return rigorous_inverse(iv, b_mat)


def _to_local_box(iv: Any, binv: IMat, center: list[Any], image_box: list[Any]) -> tuple[Any, Any]:
    """Map an ambient interval box into an h-set's local (unstable, stable) coords.

    ``image_box`` is ``[x-interval, xdot-interval]`` (an enclosure in ambient
    section coordinates); returns ``(a_interval, b_interval)`` enclosing
    ``B^{-1}(image_box - center)`` -- a (possibly inflated, always sound) box in the
    target h-set's local coordinates.
    """
    d = [image_box[0] - center[0], image_box[1] - center[1]]
    a = binv[0][0] * d[0] + binv[0][1] * d[1]
    b = binv[1][0] * d[0] + binv[1][1] * d[1]
    return a, b


def covering_relation_2d_local(
    iv: Any,
    *,
    au_left: Any,
    au_right: Any,
    as_box: Any,
) -> dict[str, Any]:
    """Pure topological core of the 2D covering check, in M's OWN local coordinates.

    Inputs are rigorous interval enclosures, expressed in the target h-set ``M``'s
    unstable/stable coordinates:

      * ``au_left``  -- M-UNSTABLE coordinate of ``f(N^-_left)``  (image of N's left edge)
      * ``au_right`` -- M-UNSTABLE coordinate of ``f(N^-_right)`` (image of N's right edge)
      * ``as_box``   -- M-STABLE   coordinate of ``f(N)``         (image of the whole set)

    Verdict ``covers = (S) and (U)`` with

      (S) STABLE containment:  ``as_box subset (-1, 1)``  (image pinched inside M), and
      (U) UNSTABLE exit:  the two edges map strictly beyond OPPOSITE unstable ends of M
          -- ``(au_left < -1 and au_right > +1)`` or the mirror.

    This is the ZG "correctly aligned" sufficient condition for ``N ==f==> M`` -- a
    topological (degree) statement.  Over-approximated inputs stay on the sound side of
    every strict inequality, so ``covers=True`` is a genuine covering CERTIFICATE.  Keeping
    the check in M's coordinates (rather than transforming an ambient box hull) is what
    makes it valid when M's stable/unstable directions are NOT axis-aligned: the caller
    must supply ``as_box`` computed WITHOUT decorrelating the (generically sheared) image.
    """
    cond_s = bool(as_box.a > -1) and bool(as_box.b < 1)
    orient_pos = bool(au_left.b < -1) and bool(au_right.a > 1)
    orient_neg = bool(au_left.a > 1) and bool(au_right.b < -1)
    cond_u = orient_pos or orient_neg
    return {
        "covers": bool(cond_s and cond_u),
        "stable_containment": cond_s,
        "unstable_exit": cond_u,
        "orientation": ("preserving" if orient_pos else ("reversing" if orient_neg else None)),
        "stable_coord_image": [float(as_box.a), float(as_box.b)],
        "unstable_coord_left": [float(au_left.a), float(au_left.b)],
        "unstable_coord_right": [float(au_right.a), float(au_right.b)],
    }


def covering_relation_2d(
    iv: Any,
    *,
    m_center: list[Any],
    m_uvec: list[Any],
    m_svec: list[Any],
    image_left: list[Any],
    image_right: list[Any],
    image_box: list[Any],
) -> dict[str, Any]:
    """2D covering check from AMBIENT-coordinate image box enclosures.

    Convenience wrapper over :func:`covering_relation_2d_local`: transforms the ambient
    ``[x-iv, xdot-iv]`` image enclosures of ``f(N^-_left)``, ``f(N^-_right)``, ``f(N)``
    into ``M``'s local coordinates (via the rigorous frame inverse ``B_M^{-1}``) and
    applies the covering conditions.

    IMPORTANT rigor caveat: transforming an ambient interval BOX into M's coordinates is
    EXACT only when M's frame is axis-aligned with the ambient section coordinates; for a
    rotated M the box hull DECORRELATES the (generically thin, sheared) image and inflates
    its stable coordinate, so a genuine covering can fail to verify through this wrapper.
    When M's directions are rotated, supply M-local enclosures that preserve the image's
    shear (a parallelepiped, not a box hull) and call :func:`covering_relation_2d_local`
    directly.  A ``covers=True`` verdict from this wrapper is always sound; a ``False`` may
    be an artifact of the box-hull decorrelation, not a true non-covering.
    """
    binv = hset_chart_inverse(iv, m_uvec, m_svec)
    if binv is None:
        return {"covers": False, "reason": "degenerate target h-set frame M"}
    au_left, _ = _to_local_box(iv, binv, m_center, image_left)
    au_right, _ = _to_local_box(iv, binv, m_center, image_right)
    _, as_box = _to_local_box(iv, binv, m_center, image_box)
    return covering_relation_2d_local(iv, au_left=au_left, au_right=au_right, as_box=as_box)


# --------------------------------------------------------------------------- #
# #674 -- CORRELATION-PRESERVING (mean-value) section-map image of a sub-box.   #
#                                                                             #
# #673's covering attempt hit a representation wall: propagating each face of   #
# an h-set (whole box / left edge / right edge) INDEPENDENTLY through the        #
# section map and box-hulling each image re-introduces the wrapping effect --    #
# the axis-aligned hull of the thin, sheared true image is ~9x too wide, and     #
# that slop masks the genuine hyperbolic edge separation.  The fix is the same   #
# discipline #669's QR reframing applied to the FLOW enclosure, now applied to   #
# how a single section-map IMAGE is described: carry the image in a correlation- #
# preserving (parallelepiped) form built from ONE rigorous linearization of the  #
# map rather than re-boxing per face.                                            #
#                                                                             #
# Rigorous mean-value enclosure.  Run the section map ONCE on the whole h-set    #
# box (giving [DP], a rigorous enclosure of the section-map Jacobian over the     #
# whole IC box W0) and ONCE on the box CENTER (giving a TIGHT image P(w_hat)).    #
# Then, by the interval mean-value theorem, for every IC ``w0`` in ``W0`` (in     #
# particular every point of any face / sub-box):                                 #
#                                                                             #
#     P(w0)  in  P(w_hat)  +  [DP] . (w0 - w_hat).                                #
#                                                                             #
# The face's parameter structure enters only through the small, tight offset     #
# ``w0 - w_hat`` -- so the large per-step stretching that wraps a direct box      #
# propagation is applied ONCE (as the tight [DP]) instead of accumulating, and    #
# the image width collapses toward the true linear stretch.  The result is a      #
# genuine over-enclosure of P(face) (mean value is exact to first order; [DP]'s   #
# range over W0 rigorously bounds the second-order variation), so feeding it to   #
# the covering checker keeps every verdict sound.                                 #
# --------------------------------------------------------------------------- #
def section_map_meanvalue_image(
    iv: Any,
    center_image: list[Any],
    jac_box: IMat,
    offset: list[Any],
) -> list[Any]:
    """Correlation-preserving (mean-value) enclosure of a sub-box's section-map image.

    ``center_image`` is a TIGHT enclosure of ``P(w_hat)`` (the section-map image of the
    IC-box CENTER, from a center-point run); ``jac_box`` is a rigorous enclosure of the
    section-map Jacobian ``DP`` over the WHOLE IC box (from a whole-box run's
    ``section_jacobian``); ``offset`` is the regularized IC offset ``w0_face - w_hat`` of
    the face/sub-box whose image is wanted (a thin interval vector).  Returns the rigorous
    interval enclosure

        P(face)  subset  center_image + jac_box @ offset

    valid for every IC in the (convex) sub-box, by the interval mean-value theorem
    (``jac_box`` bounds ``DP`` over the convex hull of ``{w_hat} U face subset W0``).
    Unlike an independent box propagation of the face, the map's stretching is applied
    once through the tight ``jac_box`` rather than accumulated with per-step wrapping, so
    the enclosure preserves the image's true (thin, sheared) shape.
    """
    n = len(center_image)
    out: list[Any] = []
    for i in range(n):
        acc = center_image[i]
        for k in range(n):
            acc = acc + jac_box[i][k] * offset[k]
        out.append(acc)
    return out


# --------------------------------------------------------------------------- #
# #675 -- COMPOSED (multi-return) mean-value section-map image.                 #
#                                                                             #
# #674 certified a SINGLE-return covering only at a near-vacuous h-set size    #
# (ru=1e-8); at the more meaningful ru=1e-6 the single-return covering fails    #
# (ratio 0.28) because [DP] genuinely varies ~0.9 across the box.  A genuine    #
# W-Z proof instead COMPOSES many section returns: the map's hyperbolic stretch #
# compounds multiplicatively (P^N stretch ~ lambda^N), and the composed         #
# mean-value image is built by the CHAIN RULE with the SAME one-linearization   #
# discipline as #674.                                                           #
#                                                                             #
# Rigor.  For the N-fold composed section map P^N = P o ... o P and any IC w0    #
# in the (convex) h-set W0, the interval mean-value theorem gives                #
#                                                                             #
#     P^N(w0)  in  P^N(w_hat)  +  [D(P^N) over W0] . (w0 - w_hat),               #
#                                                                             #
# and by the chain rule  D(P^N)(w) = DP(P^{N-1}(w)) ... DP(P(w)) DP(w).  If      #
# ``jac_boxes = [DP_1, ..., DP_N]`` are rigorous section-Jacobian enclosures     #
# with DP_1 taken over W0, DP_2 over an enclosure of P(W0), ..., DP_N over an     #
# enclosure of P^{N-1}(W0) (each leg's whole-box run), then for every w in W0     #
# the k-th factor DP(P^{k-1}(w)) lies in DP_k, so the interval matrix product     #
# ``DP_N ... DP_1`` (inclusion-monotone) rigorously encloses D(P^N) over W0.      #
# ``compose_section_jacobians`` forms that product; feeding it (with the tight    #
# centre image P^N(w_hat) from the centre chain) to ``section_map_meanvalue_image #
# `` yields a sound enclosure of P^N(face) at the ORIGINAL h-set's face offset.   #
# --------------------------------------------------------------------------- #
def compose_section_jacobians(iv: Any, jac_boxes: list[IMat]) -> IMat:
    """Chain-rule composition ``DP_N @ DP_{N-1} @ ... @ DP_1`` of per-leg section
    Jacobian enclosures.

    ``jac_boxes[k]`` is a rigorous interval enclosure of the section-map Jacobian over
    the (k+1)-th leg's IC box (leg 1 over the h-set ``W0``, leg 2 over an enclosure of
    ``P(W0)``, ...).  Returns a rigorous enclosure of the composed map ``P^N``'s
    Jacobian over ``W0`` by the chain rule -- the matrix product is inclusion-monotone,
    so composing the leg enclosures is sound.  ``len(jac_boxes) == 1`` reproduces the
    single-return Jacobian unchanged.  Pair with :func:`section_map_meanvalue_image`
    (centre image = the composed centre chain's final crossing state) for the composed
    correlation-preserving image of a sub-box.
    """
    if not jac_boxes:
        raise ValueError("compose_section_jacobians needs at least one leg Jacobian")
    composed = [row[:] for row in jac_boxes[0]]
    for k in range(1, len(jac_boxes)):
        composed = _imatmul(iv, jac_boxes[k], composed)
    return composed


__all__ = [
    "HAVE_MPMATH",
    "apriori_enclosure",
    "compose_section_jacobians",
    "covering_relation_2d",
    "covering_relation_2d_local",
    "cr3bp_planar_jacobi",
    "cr3bp_planar_jet",
    "cr3bp_planar_variational_jet",
    "deriv_series",
    "enclosure_box",
    "eval_series_horner",
    "flow_taylor_model",
    "hset_chart_inverse",
    "integrate",
    "integrate_c1_qr",
    "isolate_section_crossing",
    "lc_secondary_from_physical",
    "lc_secondary_to_physical",
    "make_cr3bp_lc_secondary_jet",
    "make_cr3bp_lc_secondary_variational_jet",
    "rigorous_inverse",
    "rigorous_section_map",
    "section_map_jacobian",
    "section_map_meanvalue_image",
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
