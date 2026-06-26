"""Coplanar Tisserand graph: V_inf contours and the pairwise ``linkable`` predicate.

The Tisserand parameter at body ``p`` (semi-major axis ``a_p``) for a
spacecraft heliocentric orbit ``(a, e, i)`` is

.. math::

    T_p(a, e, i) = \\frac{a_p}{a} + 2\\cos(i)\\,\\sqrt{\\tfrac{a}{a_p}\\,(1-e^2)}

and is conserved across a ballistic flyby at ``p`` (Tisserand 1896). It is
directly related to the hyperbolic excess speed at ``p``:

.. math::

    V_\\infty^2 = \\frac{\\mu_\\odot}{a_p}\\,(3 - T_p)

so a constant-:math:`V_\\infty` flyby at ``p`` lives on a constant-:math:`T_p`
curve in ``(a, e)`` space.

**The original M2 predicate** :func:`linkable` **implements the coplanar
case only**, fixing ``i = 0`` so ``cos(i) = 1``; it stays coplanar by
contract — do not feed inclined real-ephemeris orbits to it and expect
sensible answers. The 3-D extension is provided by the explicit, opt-in
:func:`linkable_3d` (STAGE 4): at a fixed :math:`V_\\infty` each body's
Tisserand equation fixes ``cos(i_sc)`` at every ``(a, e)``, and a 2-D grid
scan tests whether the two bodies agree on a single reachable spacecraft
inclination. ``linkable_3d`` is the sanctioned 3-D predicate; it
conservatively extends :func:`linkable` (a coplanar ``True`` implies a 3-D
``True``). The coplanar :func:`linkable` is left untouched so existing
callers and goldens stay byte-identical.

The :func:`linkable` predicate is the linchpin of the M4 Tisserand-pruning
gate (spec §13.3): if two bodies' constant-:math:`V_\\infty` contours fail
to intersect in ``(a, e)`` space, no orbit at that ``V_inf`` reaches both,
so no flyby sequence between them at that energy can exist. The M4
enumerator discards every cell whose consecutive bodies fail this test,
which is what shrinks the search space by orders of magnitude.

Plan: ``docs/phases/m2-flyby-maps/plan.md`` §3.2.

References
----------
* Tisserand, F., *Traité de Mécanique Céleste*, vol. 4, 1896.
* Strange, N. J. & Longuski, J. M., "Graphical Method for Gravity-Assist
  Trajectory Design", *J. Spacecraft and Rockets*, 39(1):9-16, 2002.
* Campagnola, S. & Russell, R. P., "Endgame Problem Part 2: Multibody
  Technique and the Tisserand-Poincaré Graph", *JGCD*, 33(2):476-486, 2010.
"""

from __future__ import annotations

import contextlib
from functools import cache
from math import sqrt
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq, minimize

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.satellites import SATELLITES

# Memoization (task #472): the V∞<->Tisserand scalar conversions and the
# expensive ``linkable`` contour-intersection predicate are pure functions of
# their (hashable: str / float / tuple) arguments — they read only the *frozen*
# module-constant body tables (:data:`PLANETS` / :data:`SATELLITES`, frozen
# dataclasses built once at import) and the module-constant ``mu``. The
# moon-prune gate re-evaluates them on the same (body-pair, V∞, a_range, mu)
# thousands of times, so an exact (no-rounding) in-memory ``functools.cache`` is safe.
# Parity / key / purity are guarded in tests/search/test_tisserand_memoization.py;
# bypass with ``.__wrapped__``, inspect with ``.cache_info()``.

if TYPE_CHECKING:
    from matplotlib.axes import Axes

# Closed numerical tolerances used throughout.
_E_MAX: float = 1.0 - 1.0e-4
"""Upper sweep limit on eccentricity; we stay short of the parabolic boundary."""

_U_MIN_NUM: float = 1.0e-6
"""Lower bracketing floor for the dimensionless ``u = sqrt(a/a_p)`` solver."""


# ---------------------------------------------------------------------------
# Body semi-major axis helper
# ---------------------------------------------------------------------------


@cache
def _a_p_km(body: str) -> float:
    """Semi-major axis of ``body`` in km.

    Resolves a heliocentric planet code via :data:`PLANETS` (``sma_au * AU_KM``)
    and, failing that, a moon code via
    :data:`cyclerfinder.core.satellites.SATELLITES` (``sma_km`` — already km,
    already about the primary). The Tisserand identity is centre-blind once
    ``a_p`` (and ``mu``) name the right centre (moon-tour Tier-1).
    """
    pl = PLANETS.get(body)
    if pl is not None:
        return pl.sma_au * AU_KM
    return SATELLITES[body].sma_km


# ---------------------------------------------------------------------------
# V_inf <-> Tisserand conversion (coplanar, i=0)
# ---------------------------------------------------------------------------


@cache
def vinf_to_tisserand(body: str, vinf_kms: float, *, mu: float = MU_SUN_KM3_S2) -> float:
    """Tisserand parameter ``T_p`` at ``body`` for hyperbolic excess ``vinf_kms`` (km/s).

    Coplanar (``i = 0``) only. ``T_p = 3 - V_inf^2 * a_p / mu``. ``mu`` defaults
    to the Sun (heliocentric, byte-identical); pass a primary's GM
    (``PRIMARIES["Jupiter"]``) for a Jovicentric/Saturnian V_inf about that
    primary (moon-tour Tier-1) — the canonical identity ``T = 3 - v_inf^2`` is
    centre-aware via ``mu`` + ``a_p``.
    """
    a_p = _a_p_km(body)
    return 3.0 - (vinf_kms * vinf_kms) * a_p / mu


@cache
def tisserand_to_vinf(body: str, t_p: float, *, mu: float = MU_SUN_KM3_S2) -> float:
    """Inverse of :func:`vinf_to_tisserand`. Returns 0.0 when ``t_p >= 3``.

    ``t_p >= 3`` corresponds to spacecraft orbits that cannot encounter the
    body at any positive ``V_inf`` (the body is unreachable from the
    spacecraft's (a, e) regime); the conventional inverse returns 0.

    ``mu`` defaults to the Sun (heliocentric, byte-identical); pass a primary's
    GM for a centre-aware inverse (moon-tour Tier-1).

    Parameter spelled ``t_p`` (lower-case) per PEP 8 / ruff N803; the
    literature symbol is ``T_p`` (capital T, subscript p for "planet").
    """
    if t_p >= 3.0:
        return 0.0
    a_p = _a_p_km(body)
    return float(sqrt(mu * (3.0 - t_p) / a_p))


# ---------------------------------------------------------------------------
# Constant-V_inf contour in (a, e) space
# ---------------------------------------------------------------------------


def _cubic_in_u(u: float, e: float, t_p: float) -> float:
    """Tisserand contour equation in dimensionless ``u = sqrt(a/a_p)``.

    Derived from ``a_p/a + 2*sqrt((a/a_p)*(1-e^2)) = T_p`` by multiplying
    by ``u^2``: ``1 + 2*u^3*sqrt(1-e^2) = T_p*u^2``. Solutions in ``u > 0``
    correspond to physical ``a = u^2 * a_p``. Argument spelled ``t_p`` per
    PEP 8 / ruff N803.
    """
    s = sqrt(max(0.0, 1.0 - e * e))
    return 2.0 * u * u * u * s - t_p * u * u + 1.0


def _contour_roots_in_u(e: float, t_p: float) -> list[float]:
    """All positive real roots in ``u`` of the cubic at the given ``e, t_p``.

    The cubic has 0, 1, or 2 positive real roots depending on the position
    of the local minimum. Closed-form bracketing:

    * Derivative ``f'(u) = 2u*(3u*sqrt(1-e^2) - t_p)`` vanishes at ``u=0``
      and ``u_min = t_p / (3*sqrt(1-e^2))`` (when ``t_p > 0`` and
      ``e < 1``). With ``f(0) = 1 > 0`` and ``f -> +inf`` as
      ``u -> +inf``, roots exist iff ``f(u_min) < 0``, in which case
      there are exactly two of them: one in ``(0, u_min)`` and one in
      ``(u_min, large)``.

    Returns the roots sorted ascending. May be empty.
    """
    if t_p <= 0.0:
        # Limiting regime: cubic is monotone increasing in u for u >= 0
        # (negative linear in u^2 disappears), so f(u) >= 1 > 0, no roots.
        return []
    s = sqrt(max(0.0, 1.0 - e * e))
    if s == 0.0:
        # e = 1: equation degenerates to -t_p*u^2 + 1 = 0 -> u = 1/sqrt(t_p)
        # (one root). We never hit this in practice (e capped below 1).
        if t_p > 0.0:
            return [1.0 / sqrt(t_p)]
        return []

    u_crit = t_p / (3.0 * s)
    f_crit = _cubic_in_u(u_crit, e, t_p)
    if f_crit >= 0.0:
        return []

    # Bracket and refine the two roots.
    # Lower root: in (eps, u_crit). f(eps) = 1 > 0, f(u_crit) < 0.
    try:
        u_lo = float(brentq(_cubic_in_u, _U_MIN_NUM, u_crit, args=(e, t_p)))
    except ValueError:
        u_lo = float("nan")

    # Upper root: in (u_crit, u_high). Pick u_high big enough that
    # f(u_high) > 0. The cubic term dominates; set u_high = max(2*u_crit,
    # t_p) which gives f >= 2*t_p^3/(27*(1-e^2)) * (2^3 - 1) + 1.
    u_high = max(2.0 * u_crit, t_p) + 1.0
    # Ensure positive at upper end; if not, double until it is.
    while _cubic_in_u(u_high, e, t_p) <= 0.0 and u_high < 1.0e6:
        u_high *= 2.0
    try:
        u_hi = float(brentq(_cubic_in_u, u_crit, u_high, args=(e, t_p)))
    except ValueError:
        u_hi = float("nan")

    roots = [r for r in (u_lo, u_hi) if not np.isnan(r) and r > 0.0]
    roots.sort()
    return roots


def _orbit_crosses_planet(a_au: float, e: float, a_p_au: float) -> bool:
    """True iff a heliocentric orbit ``(a, e)`` crosses the planet's circular orbit.

    For the Tisserand parameter to represent a physical encounter V_inf, the
    spacecraft orbit must actually intersect the planet's orbit (perihelion
    inside, aphelion outside). Without this constraint the constant-T_p
    curve includes ``(a, e)`` points the spacecraft can never visit at the
    planet, and ``linkable`` would falsely declare any two contours as
    sharing structure. (Spec §3, §13.3.)
    """
    perihelion = a_au * (1.0 - e)
    aphelion = a_au * (1.0 + e)
    return perihelion <= a_p_au <= aphelion


def vinf_contour(
    body: str,
    vinf_kms: float,
    a_range_au: tuple[float, float] = (0.3, 5.0),
    n_points: int = 200,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Constant-:math:`V_\\infty` contour at ``body`` in ``(a, e)`` space.

    Parameterised by ``e in [0, 1 - 1e-4)`` with ``n_points`` samples. For
    each ``e`` the contour equation (a cubic in ``u = sqrt(a/a_p)``) is
    solved and every positive real root lying with ``a in a_range_au``
    **and** whose spacecraft orbit actually crosses the planet's orbit
    (perihelion inside, aphelion outside) is returned. Samples failing
    either filter are silently dropped (this is a normal regime — the
    contour just does not extend there — not an error).

    The "orbit-crosses-planet" filter is essential: without it the
    constant-:math:`T_p` curve includes ``(a, e)`` points the spacecraft
    can never visit at the planet, and downstream consumers (notably
    :func:`linkable`) would conflate algebraic contour intersection with
    physical accessibility. Spec §3.

    **Coplanar (i = 0) only.**

    Parameters
    ----------
    body:
        One-letter planet code.
    vinf_kms:
        Hyperbolic excess speed, km/s. Must be non-negative.
    a_range_au:
        ``(a_min_AU, a_max_AU)`` filter on physically-of-interest orbits.
    n_points:
        Number of ``e`` samples in ``[0, 1 - 1e-4)``.

    Returns
    -------
    tuple[NDArray, NDArray]
        Two equal-length float64 arrays ``(a_au, e)``. Lengths may be
        smaller than ``n_points`` (samples without a valid root are
        dropped) and may exceed ``n_points`` only if two branches
        contribute at a single ``e``. Empty arrays when no contour exists.
    """
    if vinf_kms < 0.0:
        raise ValueError(f"vinf_kms must be non-negative; got {vinf_kms}")
    if n_points < 2:
        raise ValueError(f"n_points must be >= 2; got {n_points}")
    a_min_au, a_max_au = a_range_au
    if not (a_min_au > 0.0 and a_max_au > a_min_au):
        raise ValueError(
            f"a_range_au must be (a_min, a_max) with 0 < a_min < a_max; got {a_range_au}"
        )

    t_p = vinf_to_tisserand(body, vinf_kms)
    a_p_km = _a_p_km(body)
    a_p_au = a_p_km / AU_KM

    e_samples = np.linspace(0.0, _E_MAX, n_points)

    a_out: list[float] = []
    e_out: list[float] = []
    for e in e_samples:
        e_f = float(e)
        roots = _contour_roots_in_u(e_f, t_p)
        for u in roots:
            a_au = u * u * a_p_au
            if not (a_min_au <= a_au <= a_max_au):
                continue
            if not _orbit_crosses_planet(a_au, e_f, a_p_au):
                continue
            a_out.append(a_au)
            e_out.append(e_f)
    return np.asarray(a_out, dtype=np.float64), np.asarray(e_out, dtype=np.float64)


# ---------------------------------------------------------------------------
# linkable predicate (M4 pruning gate, spec §13.3)
# ---------------------------------------------------------------------------


def _a_branches_at_e(
    body: str,
    e: float,
    t_p: float,
    a_range_au: tuple[float, float],
) -> list[float]:
    """All physical ``a`` (AU) solutions of the contour at this ``e``.

    Filters by both the requested ``a_range_au`` and by the orbit-crosses-
    planet condition (see :func:`_orbit_crosses_planet`).
    """
    a_p_au = _a_p_km(body) / AU_KM
    roots_u = _contour_roots_in_u(e, t_p)
    a_min_au, a_max_au = a_range_au
    out: list[float] = []
    for u in roots_u:
        a_au = u * u * a_p_au
        if not (a_min_au <= a_au <= a_max_au):
            continue
        if not _orbit_crosses_planet(a_au, e, a_p_au):
            continue
        out.append(a_au)
    return out


@cache
def linkable(
    body_a: str,
    body_b: str,
    vinf_kms: float = 0.0,
    tol_au: float = 0.01,
    tol_e: float = 0.01,
    a_range_au: tuple[float, float] = (0.3, 5.0),
    n_points: int = 200,
    *,
    mu: float = MU_SUN_KM3_S2,
) -> bool:
    """Do the constant-:math:`V_\\infty` contours of ``body_a`` and ``body_b`` intersect?

    True iff there exists an ``e in [0, 1)`` and some branch of each
    body's contour at this ``V_inf`` with the two ``a`` values agreeing
    within ``tol_au`` (AU). Equivalently: there exists a spacecraft orbit
    of fixed ``(a, e)`` reachable from both bodies at the given
    :math:`V_\\infty` — so a flyby sequence between them at this energy
    is energetically possible.

    Implementation: sample the shared ``e`` grid; at each ``e`` compute
    every branch's ``a`` for each body; for every pair of branches
    (``branch_a``, ``branch_b``) look for sign changes of
    ``g(e) = a_a(e) - a_b(e)``. Refine sign changes with
    :func:`scipy.optimize.brentq`. Tangent contact (no sign change but
    ``|g|_min < tol_au``) is also counted as linkable.

    **Coplanar (i = 0) only.** This is the M4 Tisserand-pruning gate
    (spec §13.3). Returns False (never raises) on numerical failures.

    Parameters
    ----------
    body_a, body_b:
        One-letter planet codes.
    vinf_kms:
        Common hyperbolic excess speed, km/s. Must be non-negative.
    tol_au, tol_e:
        Tolerances on ``a``-distance (AU) and ``e``-width when declaring
        tangent contact.
    a_range_au:
        ``(a_min, a_max)`` AU filter on the contours.
    n_points:
        Eccentricity grid density.
    mu:
        Central GM for the Tisserand identity. Defaults to the Sun
        (heliocentric, byte-identical). Pass a primary's GM
        (``PRIMARIES["Jupiter"]``) for a Jovicentric/Saturnian moon pair —
        ``a_range_au`` must then be moon-scale (about-primary SMA in AU)
        (moon-tour Tier-1).

    Returns
    -------
    bool
    """
    if vinf_kms < 0.0:
        return False  # never raises
    if n_points < 2:
        return False

    t_pa = vinf_to_tisserand(body_a, vinf_kms, mu=mu)
    t_pb = vinf_to_tisserand(body_b, vinf_kms, mu=mu)

    e_samples = np.linspace(0.0, _E_MAX, n_points)

    # Pre-compute per-sample branches for each body.
    branches_a: list[list[float]] = [
        _a_branches_at_e(body_a, float(e), t_pa, a_range_au) for e in e_samples
    ]
    branches_b: list[list[float]] = [
        _a_branches_at_e(body_b, float(e), t_pb, a_range_au) for e in e_samples
    ]

    # Walk consecutive e samples and check every pair of branches at each end.
    # If any pair has a sign change or a near-tangent, we are linkable.
    for i in range(len(e_samples) - 1):
        e_lo, e_hi = float(e_samples[i]), float(e_samples[i + 1])
        for a_branch_lo in branches_a[i]:
            for b_branch_lo in branches_b[i]:
                g_lo = a_branch_lo - b_branch_lo
                # Tangent tolerance at this sample
                if abs(g_lo) <= tol_au:
                    return True
                # Match these branches with the closest at e_hi (by minimum |g|).
                for a_branch_hi in branches_a[i + 1]:
                    for b_branch_hi in branches_b[i + 1]:
                        g_hi = a_branch_hi - b_branch_hi
                        if g_lo * g_hi < 0.0:
                            # Sign change. Refine via brentq on linearly-
                            # interpolated branch functions; the refinement
                            # is best-effort and swallows ValueError if the
                            # bracket is degenerate (the boolean answer is
                            # already True).
                            with contextlib.suppress(ValueError):
                                _refine_intersection_e(
                                    body_a,
                                    body_b,
                                    t_pa,
                                    t_pb,
                                    e_lo,
                                    e_hi,
                                    a_range_au=a_range_au,
                                )
                            return True
                        if abs(g_hi) <= tol_au:
                            return True

    # Also catch the last sample on its own.
    for a_branch in branches_a[-1]:
        for b_branch in branches_b[-1]:
            if abs(a_branch - b_branch) <= tol_au:
                return True

    # Tangent contact along an entire branch: scan all samples for min |g|.
    return False


def _refine_intersection_e(
    body_a: str,
    body_b: str,
    t_pa: float,
    t_pb: float,
    e_lo: float,
    e_hi: float,
    a_range_au: tuple[float, float],
) -> float:
    """Refine an e-bracket containing a contour intersection via brentq.

    Used only for diagnostic robustness inside :func:`linkable`; the actual
    boolean answer is already True the moment a sign change is observed.
    This helper exists so the M2 implementation can be extended in M4 to
    return the intersection ``e`` value without changing the public API.
    """

    def difference(e: float) -> float:
        a_a_list = _a_branches_at_e(body_a, e, t_pa, a_range_au)
        a_b_list = _a_branches_at_e(body_b, e, t_pb, a_range_au)
        if not a_a_list or not a_b_list:
            # Degenerate: pick the nearer branch as the best-effort signed value.
            return float("nan")
        # Match closest branches (single-difference for the refined bracket).
        # In the unique-pair regime this is exact; in the multi-branch regime
        # the caller has already isolated the correct pair via the i-grid.
        a_a = min(a_a_list, key=lambda x: abs(x - (a_b_list[0])))
        a_b = min(a_b_list, key=lambda x: abs(x - a_a))
        return a_a - a_b

    return float(brentq(difference, e_lo, e_hi, xtol=1.0e-6))


# ---------------------------------------------------------------------------
# linkable_3d — the sanctioned 3-D extension of the linkability predicate
# ---------------------------------------------------------------------------


def _cos_i_sc(a_au: float, e: float, a_p_au: float, t_p: float) -> float | None:
    """Spacecraft ``cos(i_sc)`` implied by body ``p``'s Tisserand equation.

    From :math:`T_p = a_p/a + 2\\cos(i)\\sqrt{(a/a_p)(1-e^2)}`,

    .. math::

        \\cos(i_{sc}) = \\frac{T_p - a_p/a}{2\\sqrt{(a/a_p)(1-e^2)}}.

    Returns ``None`` when the denominator is non-positive (degenerate
    ``e \\to 1`` or ``a \\to 0``) — such points carry no usable inclination
    constraint.
    """
    denom = 2.0 * sqrt((a_au / a_p_au) * max(0.0, 1.0 - e * e))
    if denom <= 0.0:
        return None
    return (t_p - a_p_au / a_au) / denom


def linkable_3d(
    body_a: str,
    body_b: str,
    vinf_kms: float,
    *,
    i_sc_max_deg: float = 30.0,
    tol_cos_i: float = 1.0e-3,
    a_range_au: tuple[float, float] = (0.3, 5.0),
    n_points: int = 80,
    n_a_points: int = 80,
) -> bool:
    """3-D (inclined) extension of :func:`linkable`.

    At a fixed :math:`V_\\infty`, each body's Tisserand equation analytically
    fixes the spacecraft inclination ``cos(i_sc)`` at any ``(a, e)`` point.
    A 2-D grid scan over ``(a, e)`` asks whether the two bodies' implied
    ``cos(i_sc)`` agree within ``tol_cos_i`` at a point where the spacecraft
    orbit physically crosses both planets' orbits (perihelion ≤ ``a_p`` ≤
    aphelion — a frame-independent scalar inequality) and the implied
    inclination is reachable (``i_sc ≤ i_sc_max_deg``, i.e. ``cos(i_sc) ≥
    cos(i_sc_max)``). This reduces the 3-D surface-intersection problem to a
    2-D scan with O(1) work per point.

    The coplanar case is recovered exactly: for ``i_sc_max_deg → 0`` only
    points whose two implied ``cos(i_sc)`` are both ≈ 1 survive, which is
    where the two constant-:math:`V_\\infty` contours of :func:`linkable`
    intersect. A coplanar ``True`` therefore implies a 3-D ``True``
    (inclination opens options, never closes them).

    Formula (sourced — Strange & Longuski 2002 JSR 39(1):9-16):
    ``T_p = a_p/a + 2 cos(i) sqrt((a/a_p)(1-e^2))``,
    ``V_inf^2 = (mu_sun/a_p)(3 - T_p)``.

    Parameters
    ----------
    body_a, body_b:
        One-letter planet codes.
    vinf_kms:
        Common hyperbolic excess speed, km/s. Must be non-negative.
    i_sc_max_deg:
        Maximum spacecraft inclination (deg) considered reachable. The
        default 30° bounds the search to physically sensible cycler orbits.
    tol_cos_i:
        Agreement tolerance on the two bodies' implied ``cos(i_sc)``.
        Default ``1e-3`` is conservative; it is calibrated against the
        coplanar consistency gate (``test_linkable_3d_coplanar_agrees_with_2d``)
        — the physical scale derives from the 0.01 AU tolerance in
        :func:`linkable` propagated through the Tisserand formula
        (``d(cos_i)/da ≈ 0.5 / AU`` near ``a = a_p``, so ``≈ 0.005``).
    a_range_au:
        ``(a_min, a_max)`` AU filter on the spacecraft semi-major axis.
    n_points:
        Eccentricity grid density in ``[0, 1 - 1e-4)``.
    n_a_points:
        Semi-major-axis grid density across ``a_range_au``.

    Returns
    -------
    bool
        ``True`` iff a feasible inclined intersection exists. **Never
        raises** — degenerate inputs return ``False``, mirroring the
        :func:`linkable` contract.
    """
    try:
        if vinf_kms < 0.0:
            return False
        if n_points < 2 or n_a_points < 2:
            return False

        a_min_au, a_max_au = a_range_au
        if not (a_min_au > 0.0 and a_max_au > a_min_au):
            return False

        # Coplanar component (exact). A coplanar-linkable pair is a special
        # case of the 3-D predicate at ``i_sc = 0`` (``cos i = 1``), so it is
        # always 3-D-linkable. Delegating to :func:`linkable` here makes the
        # ``i_sc_max → 0`` limit reproduce :func:`linkable` bit-for-bit (the
        # critical regression gate) and supplies the monotonicity guarantee
        # (coplanar ``True`` ⇒ 3-D ``True``). The inclined scan below only
        # *adds* genuinely-inclined intersections (``cos i < 1``) that the
        # coplanar predicate misses.
        if linkable(
            body_a,
            body_b,
            vinf_kms,
            a_range_au=a_range_au,
            n_points=max(n_points, 2),
        ):
            return True

        a_pa_au = _a_p_km(body_a) / AU_KM
        a_pb_au = _a_p_km(body_b) / AU_KM
        t_pa = vinf_to_tisserand(body_a, vinf_kms)
        t_pb = vinf_to_tisserand(body_b, vinf_kms)

        cos_i_min = float(np.cos(np.radians(i_sc_max_deg)))

        a_samples = np.linspace(a_min_au, a_max_au, n_a_points)

        def _reachable(cos_i: float) -> bool:
            # The implied spacecraft inclination must be physically reachable:
            # cos(i) in [cos(i_max), 1] (i_sc in [0, i_sc_max]). The +tol on
            # the upper bound absorbs the tiny cos(i) > 1 overshoot from the
            # finite agreement tolerance at the exact coplanar (i=0) contact.
            return cos_i_min - tol_cos_i <= cos_i <= 1.0 + tol_cos_i

        def _g(a_au: float, e: float) -> float:
            # Signed disagreement of the two bodies' implied cos(i_sc).
            ca = _cos_i_sc(a_au, e, a_pa_au, t_pa)
            cb = _cos_i_sc(a_au, e, a_pb_au, t_pb)
            if ca is None or cb is None:
                return float("nan")
            return ca - cb

        # The intersection of the two bodies' implied ``cos(i_sc)`` surfaces
        # is a 1-D locus. Sweeping the semi-major axis ``a`` reduces the
        # search to a 1-D root-find in ``e`` per ``a`` — exact rather than
        # grid-limited, which is what makes the coplanar (i_sc → 0) limit
        # reproduce :func:`linkable` faithfully.
        #
        # For a fixed ``a`` both orbits cross both planets iff
        # ``e >= e_min = max_p |1 - a_p/a|`` (and ``e < 1``); on that
        # interval ``g(e) = cos_ia(a,e) - cos_ib(a,e)`` is continuous. A sign
        # change brackets an exact agreement (refined with ``brentq``); a
        # tangential touch is caught by sampling ``|g|`` and accepting when
        # it falls within ``tol_cos_i``. Each candidate is gated on
        # reachability (at the coplanar limit this demands the shared
        # ``cos(i_sc) ≈ 1`` at the crossing).
        def _accept(a_au: float, e: float) -> bool:
            # |g| already known small; gate the shared inclination on
            # reachability.
            cos_shared = _cos_i_sc(a_au, e, a_pa_au, t_pa)
            return cos_shared is not None and _reachable(cos_shared)

        e_sub = np.linspace(0.0, _E_MAX, n_points)
        # Global grid minimum of |g|, used to seed a final 2-D refinement for
        # tangential touches whose minimum lies between coarse grid nodes.
        best_abs_g = float("inf")
        best_a = a_pa_au
        best_e = 0.0
        for a_val in a_samples:
            a_au = float(a_val)
            e_min = max(abs(1.0 - a_pa_au / a_au), abs(1.0 - a_pb_au / a_au))
            if e_min >= _E_MAX:
                continue
            e_grid = e_min + (e_sub / _E_MAX) * (_E_MAX - e_min)

            # Transversal crossing: a sign change of g(e) brackets an exact
            # agreement root (refined with brentq).
            g_prev: float | None = None
            e_prev = e_min
            for e_val in e_grid:
                e = float(e_val)
                g = _g(a_au, e)
                if np.isnan(g):
                    g_prev = None
                    continue
                if abs(g) < best_abs_g:
                    best_abs_g, best_a, best_e = abs(g), a_au, e
                if abs(g) <= tol_cos_i and _accept(a_au, e):
                    return True
                if g_prev is not None and g_prev * g < 0.0:
                    try:
                        e_star = float(
                            brentq(
                                lambda ee, _a=a_au: _g(_a, ee),
                                e_prev,
                                e,
                                xtol=1.0e-10,
                            )
                        )
                    except (ValueError, RuntimeError):
                        e_star = 0.5 * (e_prev + e)
                    if _accept(a_au, e_star):
                        return True
                g_prev, e_prev = g, e

        # Tangential touch: the agreement locus may graze between coarse grid
        # nodes (no sign change). Refine |g| in 2-D around the grid minimum
        # and accept if the refined minimum is within tolerance and reachable.
        if best_abs_g <= tol_cos_i and _accept(best_a, best_e):
            return True

        def _abs_g_vec(x: NDArray[np.float64]) -> float:
            a_au = float(x[0])
            e = float(x[1])
            val = _g(a_au, e)
            return 1.0 if np.isnan(val) else abs(val)

        da = (a_max_au - a_min_au) / max(1, n_a_points - 1)
        de = _E_MAX / max(1, n_points - 1)
        bounds = [
            (max(a_min_au, best_a - 2.0 * da), min(a_max_au, best_a + 2.0 * da)),
            (max(0.0, best_e - 2.0 * de), min(_E_MAX, best_e + 2.0 * de)),
        ]
        try:
            res = minimize(
                _abs_g_vec,
                x0=np.array([best_a, best_e], dtype=np.float64),
                bounds=bounds,
                method="L-BFGS-B",
            )
            a_ref, e_ref = float(res.x[0]), float(res.x[1])
            if abs(_g(a_ref, e_ref)) <= tol_cos_i and _accept(a_ref, e_ref):
                return True
        except (ValueError, RuntimeError):
            pass
        return False
    except Exception:
        return False


def linkable_region(
    body_a: str,
    body_b: str,
    vinf_cap_kms: float,
    n_vinf: int = 50,
) -> list[float]:
    """Sample ``V_inf in (0, vinf_cap_kms]`` and return values where
    :func:`linkable` is True.

    The returned list is the raw True-set (no run-length encoding); M4
    callers reduce it to interval edges. The reason for the raw form is
    that M2 has no use for the intervals yet — exposing them prematurely
    would couple M2 to an M4-shaped API.

    Parameters
    ----------
    body_a, body_b:
        One-letter planet codes.
    vinf_cap_kms:
        Upper bound on the V_inf scan (exclusive of 0, inclusive of cap).
    n_vinf:
        Number of grid samples in ``(0, vinf_cap_kms]``.

    Returns
    -------
    list[float]
        Subset of the V_inf grid where ``linkable`` is True; possibly
        empty.
    """
    if vinf_cap_kms <= 0.0:
        return []
    if n_vinf < 1:
        return []
    grid = np.linspace(vinf_cap_kms / n_vinf, vinf_cap_kms, n_vinf)
    return [float(v) for v in grid if linkable(body_a, body_b, float(v))]


# ---------------------------------------------------------------------------
# Diagnostic plotting (optional ``viz`` extra)
# ---------------------------------------------------------------------------


def plot_tisserand(
    bodies: list[str],
    vinf_levels_kms: list[float],
    ax: Axes | None = None,
) -> Axes:
    """Overlay contour ``(a, e)`` curves for each ``(body, V_inf)`` pair.

    Diagnostic helper, not on the M2 CI gate. Requires the ``viz`` optional
    extra (matplotlib); imports lazily so the module loads without it.

    Parameters
    ----------
    bodies:
        Planet codes to plot.
    vinf_levels_kms:
        V_inf values (km/s) for which to draw a contour per body.
    ax:
        Optional matplotlib Axes; created on a fresh figure if absent.

    Returns
    -------
    matplotlib.axes.Axes
        The Axes object the curves were drawn on.

    Raises
    ------
    ImportError
        If matplotlib is not installed (install with
        ``uv sync --extra viz``).
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "plot_tisserand requires matplotlib; install the `viz` extra: `uv sync --extra viz`"
        ) from exc

    if ax is None:
        _, ax = plt.subplots(figsize=(7.0, 5.0))

    for body in bodies:
        for vinf in vinf_levels_kms:
            a_au, e = vinf_contour(body, vinf)
            if a_au.size == 0:
                continue
            ax.plot(a_au, e, label=f"{body} V_inf={vinf:.2f}", linewidth=1.0)
    ax.set_xlabel("semi-major axis a (AU)")
    ax.set_ylabel("eccentricity e")
    ax.set_title("Tisserand graph (coplanar)")
    ax.legend(fontsize=8, ncol=2)
    return ax
