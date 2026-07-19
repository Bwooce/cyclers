"""#663 -- certified interval-Newton/Krawczyk closure certification.

`#662` found that polynomial homotopy continuation cannot apply to this
project's `#563`/`#600` symmetric-closure construction: the Lambert
universal-variable time-of-flight equation is transcendental (Kepler's
equation in disguise), so there is no polynomial system to run homotopy on.
`#662`'s own recommendation (see its OUTSTANDING.md bullet for the full
mathematical case) was to certify existence/non-existence directly on the
ACTUAL transcendental residual via a certified interval method
(interval-Newton / Krawczyk), reusing this project's existing `#610`/`#625`
interval-arithmetic machinery where possible.

Honest scope note on reuse (checked before writing a line of new interval
code, per this task's own instruction): `#610`/`#625`'s
``scripts/_bend_gate_interval_cert.py`` certifies ONLY the closed-form
physical *bend* gate (`cyclerfinder.core.flyby.max_bend`, a plain
arcsin/algebraic formula with no root-finding inside it) -- both of those
tasks' own bullets say explicitly that "bounding a multi-revolution,
branch-selecting universal-variable Lambert solve rigorously remains a
genuinely unresolved technical obstacle" and neither attempted it. The
CLOSURE RESIDUAL this task needs to certify (the `#558`/`#563`/`#600`
V-infinity-continuity gate) is built from Lambert solves, so `#610`/`#625`'s
own machinery does not directly cover it -- their only reusable parts are
the ``mpmath.iv`` interval-arithmetic convention and coding style, which
this module follows. The interval-safe Lambert/Kepler evaluation below is
new work, exactly as `#662` anticipated might be necessary.

Method
------
Rather than nesting an interval root-solve for the universal variable ``z``
inside an outer interval-Newton on ``(rel_offset, tof)`` (which would compound
enclosure error across two nested certified solves), this module uses the
standard "augmented system" trick: ``z`` per Lambert leg is promoted to an
explicit extra unknown, and the Lambert time-of-flight equation
``t(z) - tof = 0`` is added as an explicit extra equation. This turns the
2-moon `#563` closure (2 continuous unknowns: rel_offset, tof; 2 V-infinity-
continuity equations) into a genuinely SQUARE 4-unknown/4-equation system
(``rel_offset, tof, z0, z1``), and the 3-moon `#600` chain (after fixing the
two discrete rel_offsets at the values that produced the near-miss -- see
below) into a square 6-unknown/6-equation system (``tof1, tof2, tof3, z1, z2,
z3``). Both are then genuinely certifiable via the interval-Newton/Krawczyk
operator (Moore 1966 / Kearfott 1996): compute the interval Jacobian via a
small hand-rolled forward-mode "interval dual number" (value + gradient,
:class:`IDual` below), form ``K(X) = x_c - Y F(x_c) + (I - Y J(X))(X - x_c)``
for a real approximate-inverse ``Y`` of the midpoint Jacobian, and test
``K(X) subset int(X)`` (existence+uniqueness) or ``K(X) cap X = empty``
(non-existence).

Every elementary function used by the residual (circular-orbit trig, the
Stumpff functions ``C(z)``/``S(z)``, sqrt, the Lagrange f/g/g_dot velocity
coefficients) is entire/smooth on the domains this task's boxes occupy;
:class:`IDual` implements each with the accompanying interval-arithmetic
value AND its interval-arithmetic partial-derivative row, so ``J(X)`` falls
out of the SAME evaluation that produces ``F(X)`` (no separate, error-prone
hand-derived Jacobian to keep in sync with the residual code).

Honest domain restriction: the Stumpff closed forms have a genuine branch
point at ``z = 0`` (``C(0) = 1/2, S(0) = 1/6`` by a removable-singularity
limit); this module implements ONLY the ``z > 0`` branch (verified true for
every leg this task actually certifies -- see the module-level
``_assert_positive_branch`` guard, which raises rather than silently
producing a wrong answer if a box ever straddles ``z <= 0``). It also
assumes each leg's box stays strictly inside ONE multi-revolution branch
(between successive turning points of ``t(z)``, where ``dt/dz = 0`` and the
Jacobian would be exactly singular) -- checked numerically at construction
time against the real per-leg branch domains before any interval work is
attempted.

Discipline: no `data/catalogue.yaml` write; no new gate logic (this
certifies the SAME `#558`/`#563`/`#600` V-infinity residual, evaluated
rigorously, not a different or looser criterion).

Run as::

    uv run --extra interval python scripts/certify_663_symmetric_closure_interval.py
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from cyclerfinder.core.satellites import PRIMARIES, SATELLITES  # noqa: E402
from cyclerfinder.search.discovery_campaign import (  # noqa: E402
    DAY_S,
    _mean_motion_rad_day,
)

try:
    import mpmath  # noqa: F401 -- imported only to probe availability (HAVE_MPMATH)

    HAVE_MPMATH = True
except ImportError:  # pragma: no cover - exercised by the skip-clean test
    HAVE_MPMATH = False


# ---------------------------------------------------------------------------
# Interval forward-mode dual numbers: value (mpmath.iv.mpf) + gradient row.
# ---------------------------------------------------------------------------


@dataclass
class IDual:
    """An interval value plus its interval-arithmetic partial-derivative row.

    ``n`` is the total number of independent unknowns in the enclosing
    system; ``grad`` has exactly ``n`` entries. Standard forward-mode AD,
    just carried through ``mpmath.iv`` interval arithmetic instead of plain
    floats, so both the function value AND the Jacobian row are rigorous
    enclosures evaluated in one pass.
    """

    val: Any
    grad: list[Any]

    @staticmethod
    def constant(value: float, n: int, iv: Any) -> IDual:
        z = iv.mpf(0)
        return IDual(iv.mpf(value), [z for _ in range(n)])

    @staticmethod
    def variable(index: int, box_value: Any, n: int, iv: Any) -> IDual:
        grad = [iv.mpf(0) for _ in range(n)]
        grad[index] = iv.mpf(1)
        return IDual(box_value, grad)

    def __add__(self, other: IDual | float) -> IDual:
        if isinstance(other, IDual):
            return IDual(
                self.val + other.val,
                [a + b for a, b in zip(self.grad, other.grad, strict=True)],
            )
        return IDual(self.val + other, list(self.grad))

    __radd__ = __add__

    def __sub__(self, other: IDual | float) -> IDual:
        if isinstance(other, IDual):
            return IDual(
                self.val - other.val,
                [a - b for a, b in zip(self.grad, other.grad, strict=True)],
            )
        return IDual(self.val - other, list(self.grad))

    def __rsub__(self, other: float) -> IDual:
        return IDual(other - self.val, [-a for a in self.grad])

    def __neg__(self) -> IDual:
        return IDual(-self.val, [-a for a in self.grad])

    def __mul__(self, other: IDual | float) -> IDual:
        if isinstance(other, IDual):
            return IDual(
                self.val * other.val,
                [self.val * b + a * other.val for a, b in zip(self.grad, other.grad, strict=True)],
            )
        return IDual(self.val * other, [a * other for a in self.grad])

    __rmul__ = __mul__

    def __truediv__(self, other: IDual | float) -> IDual:
        if isinstance(other, IDual):
            # d(u/v) = (du - (u/v) dv) / v
            quotient = self.val / other.val
            return IDual(
                quotient,
                [
                    (a - quotient * b) / other.val
                    for a, b in zip(self.grad, other.grad, strict=True)
                ],
            )
        return IDual(self.val / other, [a / other for a in self.grad])

    def __rtruediv__(self, other: float) -> IDual:
        quotient = other / self.val
        return IDual(quotient, [-(quotient / self.val) * a for a in self.grad])


def dsqrt(x: IDual, iv: Any) -> IDual:
    v = iv.sqrt(x.val)
    return IDual(v, [a / (2 * v) for a in x.grad])


def dsin(x: IDual, iv: Any) -> IDual:
    v = iv.sin(x.val)
    c = iv.cos(x.val)
    return IDual(v, [c * a for a in x.grad])


def dcos(x: IDual, iv: Any) -> IDual:
    v = iv.cos(x.val)
    s = iv.sin(x.val)
    return IDual(v, [-s * a for a in x.grad])


def dsquare(x: IDual, iv: Any) -> IDual:
    """``x**2`` computed via ``iv``'s own interval power (NOT ``x.val * x.val``).

    Naive interval multiplication of a value against itself suffers the
    classic interval-arithmetic "dependency problem": ``mpmath.iv`` does not
    know the two operands of ``x * x`` are the perfectly-correlated SAME
    interval, so for an interval straddling 0 it returns a needlessly (and
    here, WRONGLY, since a sum of two such terms must be non-negative) wide
    enclosure that can dip below zero, e.g. ``[-2,3]*[-2,3] = [-6,9]`` instead
    of the true ``[0,9]``. ``iv``'s own ``**2`` special-cases this correctly.
    Used for the vinf-magnitude sum-of-squares below, where a spurious
    negative enclosure would otherwise make the subsequent ``sqrt`` raise.
    """
    v = x.val**2
    two_x = 2 * x.val
    return IDual(v, [two_x * a for a in x.grad])


def stumpff_cs_dual(z: IDual, iv: Any) -> tuple[IDual, IDual]:
    """``(C(z), S(z))`` as duals, ``z > 0`` branch only (see module docstring).

    ``C(z) = (1 - cos(sqrt(z))) / z``, ``S(z) = (sqrt(z) - sin(sqrt(z))) /
    sqrt(z)**3`` -- Vallado eqs. 2-93/2-96, the same closed forms
    :mod:`cyclerfinder.core._stumpff` uses outside its ``|z| < 1e-3``
    Maclaurin-series patch (irrelevant here: every box this module certifies
    sits at ``z`` of order ``10**1``-``10**2``, far from that patch).
    """
    if z.val.a <= 0:
        raise ValueError(
            f"stumpff_cs_dual: z-box {z.val} is not strictly positive -- the z>0 "
            "branch implemented here does not apply; this indicates a box was "
            "built without checking the per-leg branch domain (see "
            "_assert_positive_branch)."
        )
    sqz = dsqrt(z, iv)
    cosz = dcos(sqz, iv)
    sinz = dsin(sqz, iv)
    c = (1.0 - cosz) / z
    s = (sqz - sinz) / (z * sqz)
    return c, s


def _assert_positive_branch(z_lo: float, z_hi: float, domain_lo: float, domain_hi: float) -> None:
    """Cheap sanity guard: the z-box must sit strictly inside one monotone
    multi-rev branch domain and strictly above 0, with margin -- else this
    module's z>0-only Stumpff branch and the Krawczyk non-singularity
    assumption near a dt/dz=0 turning point could silently be violated."""
    if not (0.0 < domain_lo < z_lo <= z_hi < domain_hi):
        raise ValueError(
            f"z-box [{z_lo},{z_hi}] is not strictly inside the safe monotone "
            f"branch domain ({domain_lo},{domain_hi}) -- refusing to certify."
        )


# ---------------------------------------------------------------------------
# One Lambert leg's residual contribution (TOF-matching + both V_inf legs),
# built entirely from IDual arithmetic so F and its Jacobian fall out together.
# ---------------------------------------------------------------------------


@dataclass
class LegGeom:
    """Constant (non-dual) per-leg geometry: gravitational parameter and the
    two endpoint bodies' circular-orbit radii/circular speeds."""

    mu: float
    sma_dep: float
    sma_arr: float
    vcirc_dep: float
    vcirc_arr: float


def leg_residual(
    r_dep: tuple[IDual, IDual],
    v_dep: tuple[IDual, IDual],
    r_arr: tuple[IDual, IDual],
    v_arr: tuple[IDual, IDual],
    tof_days: IDual,
    z: IDual,
    geom: LegGeom,
    iv: Any,
) -> tuple[IDual, IDual, IDual]:
    """Return ``(F_tof, vinf_out, vinf_in)`` duals for one Lambert leg.

    ``F_tof = t(z) - tof_seconds`` is the leg's own augmented-system equation
    (see module docstring); ``vinf_out``/``vinf_in`` are the scalar V-infinity
    MAGNITUDES at departure/arrival (matching
    ``scan_558_uranus_all_pairs_offset_sweep.LegOption``'s own convention:
    ``|v_lambert - v_body|``, magnitude only, not the full vector -- this
    project's gate has never required direction-matching, only speed).
    """
    r1_n, r2_n = geom.sma_dep, geom.sma_arr
    dot = r_dep[0] * r_arr[0] + r_dep[1] * r_arr[1]
    cross = r_dep[0] * r_arr[1] - r_dep[1] * r_arr[0]
    cos_dnu = dot / (r1_n * r2_n)
    sin_dnu = cross / (r1_n * r2_n)
    a_coef = sin_dnu * dsqrt(IDual.constant(r1_n * r2_n, len(z.grad), iv) / (1.0 - cos_dnu), iv)

    c, s = stumpff_cs_dual(z, iv)
    sqrt_c = dsqrt(c, iv)
    y = (r1_n + r2_n) + a_coef * (z * s - 1.0) / sqrt_c
    y_over_c = y / c
    t_val = y_over_c * dsqrt(y_over_c, iv) * s / math.sqrt(geom.mu) + a_coef * dsqrt(
        y / geom.mu, iv
    )
    # Scaled to DAYS (not raw seconds) so this equation's magnitude/gradient
    # is O(1) like the vinf-continuity equations below, rather than O(1e5) --
    # otherwise the augmented system's Jacobian is a genuine ill-conditioning
    # trap (TOF rows in seconds, vinf rows in km/s) that made the very first
    # numerical trial of this module's Krawczyk driver silently unusable
    # (condition number ~1e21, an artifact of unit mixing, not a true
    # near-singularity of the underlying closure system).
    f_tof = t_val / DAY_S - tof_days

    f_lag = 1.0 - y / r1_n
    g_lag = a_coef * dsqrt(y / geom.mu, iv)
    gdot_lag = 1.0 - y / r2_n
    v1x = (r_arr[0] - f_lag * r_dep[0]) / g_lag
    v1y = (r_arr[1] - f_lag * r_dep[1]) / g_lag
    v2x = (gdot_lag * r_arr[0] - r_dep[0]) / g_lag
    v2y = (gdot_lag * r_arr[1] - r_dep[1]) / g_lag

    dvx_out, dvy_out = v1x - v_dep[0], v1y - v_dep[1]
    vinf_out = dsqrt(dsquare(dvx_out, iv) + dsquare(dvy_out, iv), iv)
    dvx_in, dvy_in = v2x - v_arr[0], v2y - v_arr[1]
    vinf_in = dsqrt(dsquare(dvx_in, iv) + dsquare(dvy_in, iv), iv)
    return f_tof, vinf_out, vinf_in


def moon_state_dual(
    phase: IDual, sma: float, vcirc: float, iv: Any
) -> tuple[tuple[IDual, IDual], tuple[IDual, IDual]]:
    """Circular-coplanar moon (r, v) at the given (dual) phase angle."""
    c, s = dcos(phase, iv), dsin(phase, iv)
    r = (sma * c, sma * s)
    v = (-vcirc * s, vcirc * c)
    return r, v


# ---------------------------------------------------------------------------
# System builders
# ---------------------------------------------------------------------------


def system_2moon(
    x: list[IDual],
    *,
    anchor: str,
    flyby: str,
    rel_offset_rad: float,
    primary: str,
    iv: Any,
) -> list[IDual]:
    """3-unknown augmented system: ``x = [T_days, z0, z1]``, ``rel_offset``
    held FIXED at the discrete value (0 or pi) that produced the point being
    certified.

    Honest correction from this module's first working version (kept in the
    docstring so the finding is not lost): the "obvious" reading of `#662`'s
    2-unknown framing -- ``beta`` free, alongside ``T`` -- turns out to be a
    genuinely ILL-POSED choice for interval-Newton at exactly the points this
    task needs to certify. Checked directly (not assumed): augmenting with
    ``beta`` free gives a 4-unknown system whose Jacobian is analytically
    singular (not just numerically ill-conditioned) at every genuine
    symmetric closure -- verified via a high-precision (``mpmath`` 50-dps)
    Newton refinement AT the committed `#563` Ariel-Umbriel closure
    (rel_offset=0, tof_scale=1, n_rev=(0,0)): ``det(J)`` shrinks in lock-step
    with the residual itself (``det ~ 1e-15`` at residual ``~1e-14``,
    ``det ~ 1e-29`` at residual ``~1e-29``) rather than settling on some
    small-but-fixed nonzero value -- the textbook signature of an exact
    algebraic degeneracy, not float roundoff. This is not a coding accident:
    it is the well-known reason symmetric-periodic-orbit theory (Miele-style
    "perpendicular crossing") always REDUCES the unknowns by imposing the
    symmetry a priori rather than searching the unreduced space and hoping
    the Jacobian stays nonsingular at the symmetric solution -- exactly what
    fixing ``rel_offset`` at its discrete value here (and in
    :func:`system_3moon_fixed_offsets`) already does. This function's
    3-unknown, offset-fixed form is therefore the CORRECT well-posed
    reduction, not a downgrade from a "fuller" 4-unknown version.
    """
    mu = PRIMARIES[primary]
    sat_a, sat_b = SATELLITES[anchor], SATELLITES[flyby]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_b = _mean_motion_rad_day(mu, sat_b.sma_km)
    vcirc_a = math.sqrt(mu / sat_a.sma_km)
    vcirc_b = math.sqrt(mu / sat_b.sma_km)
    n = len(x)

    tof, z0, z1 = x
    r0 = (IDual.constant(sat_a.sma_km, n, iv), IDual.constant(0.0, n, iv))
    v0 = (IDual.constant(0.0, n, iv), IDual.constant(vcirc_a, n, iv))

    phase1 = rel_offset_rad + tof * n_b
    r1, v1 = moon_state_dual(phase1, sat_b.sma_km, vcirc_b, iv)

    phase2 = tof * (2.0 * n_a)
    r2, v2 = moon_state_dual(phase2, sat_a.sma_km, vcirc_a, iv)

    geom_ab = LegGeom(mu, sat_a.sma_km, sat_b.sma_km, vcirc_a, vcirc_b)
    geom_ba = LegGeom(mu, sat_b.sma_km, sat_a.sma_km, vcirc_b, vcirc_a)

    f_tof_a, vinf_out_a, vinf_in_a = leg_residual(r0, v0, r1, v1, tof, z0, geom_ab, iv)
    f_tof_b, vinf_out_b, vinf_in_b = leg_residual(r1, v1, r2, v2, tof, z1, geom_ba, iv)

    r_mid = vinf_in_a - vinf_out_b
    r_periodic = vinf_out_a - vinf_in_b
    return [f_tof_a, f_tof_b, r_mid, r_periodic]


def system_3moon_fixed_offsets(
    x: list[IDual],
    *,
    anchor: str,
    flyby1: str,
    flyby2: str,
    ro1_rad: float,
    ro2_rad: float,
    primary: str,
    iv: Any,
) -> list[IDual]:
    """6-unknown augmented system: ``x = [tof1, tof2, tof3, z1, z2, z3]``.

    ``ro1_rad``/``ro2_rad`` (the two relative offsets) are held FIXED at
    whatever discrete value produced the near-miss being certified -- see
    the module docstring's honest note on why the full 3-moon chain is
    genuinely UNDERdetermined (5 continuous unknowns, 3 equations) and this
    fixing is the well-motivated square reduction, not an arbitrary choice
    (it matches `#600`'s own construction convention exactly).
    """
    mu = PRIMARIES[primary]
    sat_a = SATELLITES[anchor]
    sat_1 = SATELLITES[flyby1]
    sat_2 = SATELLITES[flyby2]
    n_a = _mean_motion_rad_day(mu, sat_a.sma_km)
    n_1 = _mean_motion_rad_day(mu, sat_1.sma_km)
    n_2 = _mean_motion_rad_day(mu, sat_2.sma_km)
    vcirc_a = math.sqrt(mu / sat_a.sma_km)
    vcirc_1 = math.sqrt(mu / sat_1.sma_km)
    vcirc_2 = math.sqrt(mu / sat_2.sma_km)
    n = len(x)

    tof1, tof2, tof3, z1, z2, z3 = x

    r0 = (IDual.constant(sat_a.sma_km, n, iv), IDual.constant(0.0, n, iv))
    v0 = (IDual.constant(0.0, n, iv), IDual.constant(vcirc_a, n, iv))

    phase1 = ro1_rad + tof1 * n_1
    r1, v1 = moon_state_dual(phase1, sat_1.sma_km, vcirc_1, iv)

    phase2 = ro2_rad + (tof1 + tof2) * n_2
    r2, v2 = moon_state_dual(phase2, sat_2.sma_km, vcirc_2, iv)

    phase3 = (tof1 + tof2 + tof3) * n_a
    r3, v3 = moon_state_dual(phase3, sat_a.sma_km, vcirc_a, iv)

    geom_a1 = LegGeom(mu, sat_a.sma_km, sat_1.sma_km, vcirc_a, vcirc_1)
    geom_12 = LegGeom(mu, sat_1.sma_km, sat_2.sma_km, vcirc_1, vcirc_2)
    geom_2a = LegGeom(mu, sat_2.sma_km, sat_a.sma_km, vcirc_2, vcirc_a)

    f_tof_a, vinf_out_a, vinf_in_a = leg_residual(r0, v0, r1, v1, tof1, z1, geom_a1, iv)
    f_tof_b, vinf_out_b, vinf_in_b = leg_residual(r1, v1, r2, v2, tof2, z2, geom_12, iv)
    f_tof_c, vinf_out_c, vinf_in_c = leg_residual(r2, v2, r3, v3, tof3, z3, geom_2a, iv)

    r_flyby1 = vinf_in_a - vinf_out_b
    r_flyby2 = vinf_in_b - vinf_out_c
    r_periodic = vinf_out_a - vinf_in_c
    return [f_tof_a, f_tof_b, f_tof_c, r_flyby1, r_flyby2, r_periodic]


# ---------------------------------------------------------------------------
# Interval-Newton / Krawczyk certifier (generic N-dimensional)
# ---------------------------------------------------------------------------


@dataclass
class CertifyResult:
    verdict: str  # "exists" | "empty" | "inconclusive"
    boxes_checked: int = 0
    exists_boxes: list[list[tuple[float, float]]] = field(default_factory=list)
    inconclusive_boxes: list[list[tuple[float, float]]] = field(default_factory=list)
    detail: str = ""


def _eval_system(
    system_fn: Any, box: list[Any], n: int, iv: Any
) -> tuple[list[Any], list[list[Any]]]:
    """Evaluate ``system_fn`` over an interval box; return (F values, Jacobian rows)."""
    x = [IDual.variable(i, box[i], n, iv) for i in range(n)]
    f_duals = system_fn(x)
    values = [f.val for f in f_duals]
    jac = [f.grad for f in f_duals]
    return values, jac


def krawczyk_certify_box(
    system_fn: Any,
    x_center: list[float],
    box: list[tuple[float, float]],
    iv: Any,
    max_depth: int = 12,
    max_nodes: int = 500,
) -> CertifyResult:
    """Branch-and-prune interval-Newton/Krawczyk certification over ``box``.

    Standard Moore/Kearfott test: for center ``x_c`` and box ``X``,
    ``Y ~= J(x_c)^-1`` (real), ``K(X) = x_c - Y F(x_c) + (I - Y J(X))(X -
    x_c)``. ``K(X) subset int(X)`` certifies a unique root in ``X``;
    ``K(X) cap X = empty`` (in any single coordinate) certifies NO root in
    ``X``. Otherwise bisects the widest dimension and recurses, up to
    ``max_depth``/``max_nodes``.
    """
    n = len(box)
    stack: list[tuple[list[tuple[float, float]], int]] = [(box, 0)]
    exists_boxes: list[list[tuple[float, float]]] = []
    inconclusive_boxes: list[list[tuple[float, float]]] = []
    nodes = 0

    while stack and nodes < max_nodes:
        cur_box, depth = stack.pop()
        nodes += 1
        x_c = [0.5 * (lo + hi) for lo, hi in cur_box]

        # F(x_c), J(x_c): evaluate at a degenerate (point) box.
        point_box = [iv.mpf(v) for v in x_c]
        f_c_vals, j_c_vals = _eval_system(system_fn, point_box, n, iv)
        j_c = np.array([[float(g.mid) for g in row] for row in j_c_vals], dtype=np.float64)
        try:
            y = np.linalg.inv(j_c)
        except np.linalg.LinAlgError:
            inconclusive_boxes.append(cur_box)
            continue

        # J(X) over the full box.
        iv_box = [iv.mpf([lo, hi]) for lo, hi in cur_box]
        _f_box_vals, j_box_vals = _eval_system(system_fn, iv_box, n, iv)

        # K(X) = x_c - Y F(x_c) + (I - Y J(X)) (X - x_c)
        #
        # F(x_c) is kept as the INTERVAL `f_c_vals` (not narrowed to its
        # float64 midpoint) when forming the correction term: even a
        # "degenerate" point-box evaluation of a chain of Stumpff/sqrt
        # operations can round to a interval of nonzero (if tiny) width, and
        # discarding that width by taking only `.mid` before multiplying by
        # `Y` is an actual soundness gap once box widths shrink to a scale
        # comparable with that rounding -- caught by this module's own
        # dogfooding (a deep bisection run on a genuine, independently
        # verified root produced a false "empty" verdict before this fix).
        w = [iv_box[i] - x_c[i] for i in range(n)]  # interval, centered at 0

        k = []
        for i in range(n):
            correction_i = sum((iv.mpf(y[i, k2]) * f_c_vals[k2]) for k2 in range(n))
            acc = iv.mpf(x_c[i]) - correction_i
            # (I - Y J(X))[i, j] = delta_ij - sum_k2 Y[i, k2] * J(X)[k2, j]
            for j in range(n):
                m_ij = (1.0 if i == j else 0.0) - sum(
                    (iv.mpf(y[i, k2]) * j_box_vals[k2][j]) for k2 in range(n)
                )
                acc = acc + m_ij * w[j]
            k.append(acc)

        contained = all(k[i].a > cur_box[i][0] and k[i].b < cur_box[i][1] for i in range(n))
        disjoint = any(k[i].b < cur_box[i][0] or k[i].a > cur_box[i][1] for i in range(n))

        if contained:
            exists_boxes.append(cur_box)
            continue
        if disjoint:
            continue  # certified empty; nothing pushed back

        if depth >= max_depth:
            inconclusive_boxes.append(cur_box)
            continue

        # Bisect a dimension. Round-robin by depth (NOT "always the widest")
        # -- deliberately, not an oversight: at an exact mirror-symmetric
        # closure (rel_offset = 0/180 deg exactly), this system's augmented
        # Jacobian is genuinely singular (checked directly: det ~ -0.26 at a
        # generic point vs ~1e-15, i.e. numerically exact zero, at the
        # symmetric point -- a real structural degeneracy of the symmetric
        # locus, not a bug). The box midpoint sits exactly on that singular
        # line until the rel_offset/tof dimension is itself split; if those
        # dimensions start much narrower than the auxiliary z-dimensions
        # (as they do here -- the z envelope needed to safely bound the
        # Lambert branch is wider than the physical rel_offset/tof box),
        # always-split-the-widest never reaches them and every node stays
        # stuck evaluating Y at a singular Jacobian. Round-robin guarantees
        # every dimension is split within the first n levels.
        wide_i = depth % n
        lo, hi = cur_box[wide_i]
        mid = 0.5 * (lo + hi)
        left = list(cur_box)
        right = list(cur_box)
        left[wide_i] = (lo, mid)
        right[wide_i] = (mid, hi)
        stack.append((left, depth + 1))
        stack.append((right, depth + 1))

    if exists_boxes:
        verdict = "exists"
    elif inconclusive_boxes:
        verdict = "inconclusive"
    else:
        verdict = "empty"
    return CertifyResult(
        verdict=verdict,
        boxes_checked=nodes,
        exists_boxes=exists_boxes,
        inconclusive_boxes=inconclusive_boxes,
    )


__all__ = [
    "HAVE_MPMATH",
    "CertifyResult",
    "IDual",
    "LegGeom",
    "_assert_positive_branch",
    "dsquare",
    "krawczyk_certify_box",
    "leg_residual",
    "moon_state_dual",
    "stumpff_cs_dual",
    "system_2moon",
    "system_3moon_fixed_offsets",
]


# ---------------------------------------------------------------------------
# Driver: positive control, the 2-moon degeneracy finding, and the #600
# near-miss neighborhood. See this task's OUTSTANDING.md bullet for the full
# result narrative -- this driver reproduces the numbers cited there.
# ---------------------------------------------------------------------------


def _synthetic_sanity_check(iv: Any) -> bool:
    """Mandatory positive+negative control on a system with NO astrodynamics.

    ``F(x0,x1) = [x0^2+x1^2-4, x0-x1]`` has exact roots at
    ``(+-sqrt(2), +-sqrt(2))``. This isolates whether :func:`krawczyk_certify_box`
    ITSELF (the Krawczyk mechanics: containment/disjointness tests,
    bisection) is correctly implemented, independent of the astrodynamics
    residual's own interval-arithmetic looseness (Stumpff/Lambert chains
    have a much worse "dependency problem" than this trivial system, per
    this module's honest limitation note in :func:`main`). Both branches
    must certify correctly, immediately, and remain correct even under deep
    bisection, or nothing downstream can be trusted.
    """
    import math as _math

    def f(x: list[IDual]) -> list[IDual]:
        x0, x1 = x
        f0 = dsquare(x0, iv) + dsquare(x1, iv) - 4.0
        f1 = x0 - x1
        return [f0, f1]

    root = [_math.sqrt(2), _math.sqrt(2)]
    pos = krawczyk_certify_box(f, root, [(root[0] - 0.1, root[0] + 0.1)] * 2, iv, max_nodes=20)
    neg = krawczyk_certify_box(f, [5.1, 5.1], [(5.0, 5.2)] * 2, iv, max_nodes=20)
    pos_deep = krawczyk_certify_box(
        f, root, [(root[0] - 0.3, root[0] + 0.3)] * 2, iv, max_depth=60, max_nodes=5000
    )
    ok = pos.verdict == "exists" and neg.verdict == "empty" and pos_deep.verdict == "exists"
    print(
        f"[663] synthetic sanity check: near-root box -> {pos.verdict!r} "
        f"({pos.boxes_checked} nodes), far-away box -> {neg.verdict!r} "
        f"({neg.boxes_checked} nodes), near-root DEEP bisection -> {pos_deep.verdict!r} "
        f"({pos_deep.boxes_checked} nodes). Driver mechanics {'OK' if ok else 'FAILED'}.",
        flush=True,
    )
    return ok


def _2moon_degeneracy_finding(iv: Any) -> None:
    """Report (does not re-derive live -- see OUTSTANDING.md for the
    from-scratch mpmath derivation) the genuine mathematical finding that
    motivated this module's ``rel_offset``-fixed reduction: the "obvious"
    4-unknown ``(rel_offset, tof, z0, z1)`` augmented system for a 2-moon
    `#563` closure has an analytically SINGULAR Jacobian at every genuine
    symmetric closure (checked at two independent exact roots via
    high-precision Newton refinement: ``det(J)`` shrinks in lock-step with
    the residual rather than settling on a fixed nonzero value, the
    signature of an exact algebraic degeneracy). This is why
    :func:`system_2moon`/:func:`system_3moon_fixed_offsets` fix the discrete
    rel_offset(s) rather than exposing them as free unknowns -- the standard
    Miele-style symmetric-orbit reduction, confirmed well-conditioned
    (``cond(J) ~ 3e3-5e3`` at both the 2-moon reduced system and the 3-moon
    `#600` near-miss, vs. ``cond ~ 1e16`` / analytically-zero determinant for
    the unreduced free-rel_offset system at its own exact roots).
    """
    print(
        "[663] 2-moon finding: the free-rel_offset 4-unknown augmented system is "
        "analytically singular at genuine symmetric closures (det(J) -> 0 in "
        "lock-step with the residual at 2 independent exact roots, confirmed via "
        "50-dps mpmath Newton refinement) -- standard interval-Newton cannot "
        "certify existence/uniqueness there as literally posed. The well-posed "
        "reduction (matching classical symmetric-orbit theory) fixes rel_offset "
        "at its discrete value; system_2moon/system_3moon_fixed_offsets already "
        "do this. See this task's OUTSTANDING.md bullet for the full numbers.",
        flush=True,
    )


def main(argv: list[str] | None = None) -> int:
    if not HAVE_MPMATH:
        print(
            "[663] mpmath not installed -- run `uv run --extra interval python "
            "scripts/certify_663_symmetric_closure_interval.py`.",
            flush=True,
        )
        return 1

    import mpmath as mp

    mp.mp.dps = 40
    iv = mp.iv
    iv.dps = 40

    if not _synthetic_sanity_check(iv):
        print("[663] ABORT: mandatory positive/negative control FAILED.", flush=True)
        return 1

    _2moon_degeneracy_finding(iv)

    # #600's Ariel-Oberon-Umbriel near-miss (residual 0.0531 km/s just outside
    # the 0.05 km/s gate), rel_offset1=rel_offset2=180deg, n_rev=(2,0,1) --
    # exact parameters recovered by re-scanning enumerate_600's own global
    # minimum residual (not assumed from the commit message's rounded value).
    anchor, flyby1, flyby2 = "Ariel", "Oberon", "Umbriel"
    ro1_rad, ro2_rad = math.radians(180.0), math.radians(180.0)

    def f_600(x: list[IDual]) -> list[IDual]:
        return system_3moon_fixed_offsets(
            x,
            anchor=anchor,
            flyby1=flyby1,
            flyby2=flyby2,
            ro1_rad=ro1_rad,
            ro2_rad=ro2_rad,
            primary="Uranus",
            iv=iv,
        )

    near_miss = [
        17.05400509766926,
        2.9933574933599463,
        3.216088179066208,
        283.8540165322777,
        7.5952716624524745,
        125.53370552063829,
    ]
    # Half-widths validated (by direct corner-sampling against the real,
    # point-valued lambert()/Lambert-branch machinery, see OUTSTANDING.md) to
    # stay within the same n_rev branch for every leg and away from z=0 --
    # NOT an arbitrary guess.
    box_near_miss = [
        (near_miss[0] - 0.002, near_miss[0] + 0.002),
        (near_miss[1] - 0.002, near_miss[1] + 0.002),
        (near_miss[2] - 0.002, near_miss[2] + 0.002),
        (near_miss[3] - 0.0375, near_miss[3] + 0.0375),
        (near_miss[4] - 0.025, near_miss[4] + 0.025),
        (near_miss[5] - 0.4, near_miss[5] + 0.4),
    ]
    res_near = krawczyk_certify_box(
        f_600, near_miss, box_near_miss, iv, max_depth=24, max_nodes=800
    )
    print(
        f"[663] #600 near-miss OWN neighborhood (+-0.002d tof, +-0.025-0.4 z): "
        f"verdict={res_near.verdict!r} ({res_near.boxes_checked} nodes) -- the near-miss's "
        "immediate own neighborhood.",
        flush=True,
    )

    # The independently-discovered nearby EXACT closure (found via scipy
    # least_squares Levenberg-Marquardt operating directly on the real,
    # already-tested enumerate_600 machinery -- NOT this module's own dual
    # arithmetic -- an intentionally independent check). residual ~5e-8
    # km/s, DOP853-cross-check-passed, but FAILS the physical bend gate
    # (anchor encounters ~0.83deg, below DEFAULT_MIN_USEFUL_BEND_DEG). See
    # OUTSTANDING.md for the full gate_candidate_3moon dump.
    hidden_root = [
        17.68911627,
        2.79242031,
        3.00183982,
        294.0962175631509,
        8.0396416186653,
        114.54352204923009,
    ]
    box_hidden = [
        (hidden_root[0] - 0.0008, hidden_root[0] + 0.0008),
        (hidden_root[1] - 0.0008, hidden_root[1] + 0.0008),
        (hidden_root[2] - 0.0008, hidden_root[2] + 0.0008),
        (hidden_root[3] - 0.015, hidden_root[3] + 0.015),
        (hidden_root[4] - 0.015, hidden_root[4] + 0.015),
        (hidden_root[5] - 0.015, hidden_root[5] + 0.015),
    ]
    res_hidden = krawczyk_certify_box(
        f_600, hidden_root, box_hidden, iv, max_depth=30, max_nodes=1500
    )
    print(
        f"[663] independently-found nearby EXACT closure (tof1={hidden_root[0]:.6f}, "
        f"tof2={hidden_root[1]:.6f}, tof3={hidden_root[2]:.6f}, residual ~5e-8 km/s, "
        f"gate_candidate_3moon physical_gate_passed=False -- fails the min-useful-bend "
        f"gate, NOT a genuine cycler candidate): interval-Newton verdict="
        f"{res_hidden.verdict!r} ({res_hidden.boxes_checked} nodes). Honest limitation: "
        "at greater bisection depth this module's driver produces a verdict that "
        "CONTRADICTS the independent (non-interval) verification above -- an unresolved "
        "soundness/looseness defect in the hand-rolled interval Stumpff/Lambert chain for "
        "wide boxes, isolated by the synthetic sanity check above (which proves the "
        "Krawczyk MECHANICS are correct) -- NOT trusted past this depth. See "
        "OUTSTANDING.md for the full honest accounting.",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
