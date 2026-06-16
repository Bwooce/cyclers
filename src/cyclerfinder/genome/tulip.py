"""Tulip-orbit genome (#266 Phase 2): sourced ICs + reproduce gate + petal classifier.

A *tulip orbit* is an Np-petal generalisation of a halo / NRHO in the CR3BP:
each "petal" is one perilune passage per period, so an Np=1 tulip is a vanilla
halo / NRHO and Np=2, 3, ... are period-multiplying bifurcations of the parent
family. The geometry resembles a Tulip viewed end-on -- each petal closes near
the secondary while the apolune loop opens out into a near-rectilinear excursion
away from the secondary.

This module is the Track-A genome surface for tulip orbits:

  * ``KOBLICK_2023_TABLE4`` -- sourced IC dictionary keyed by petal count Np.
    Each entry pins down ``(x0, y0=0, z0, xdot0=0, ydot0, zdot0=0, T_TU, ...)``
    in Koblick's CR3BP normalisation. **Provenance is mandatory**: every entry
    carries a citation comment naming the upstream source. Today the table is
    seeded with the Np=2 IC published in Coorbital's MATLAB ``pumpkyn`` package
    (``getTulip.m``), which provides the **independent cross-check anchor**
    required by the orbit-closure discipline; the Koblick 2023 AMOSTECH Table 4
    rows for the other petal counts are placeholders to be filled when the
    paper is digitised. The Np=1 parent is the widely-published Earth-Moon L2
    Southern 9:2 NRHO (NASA Gateway baseline; same seed as
    :func:`tests.search.test_reachable_impulsive._recover_92_nrho`).
  * :func:`koblick_system` -- the CR3BPSystem in Koblick's normalisation
    (mu = 1.215058560962404e-2, LU = 389703 km, TU = 382981 s). The standard
    :func:`cyclerfinder.core.cr3bp.cr3bp_system` constructs a slightly different
    Earth-Moon system from the JPL system GM; the tulip ICs are precisely
    referenced to Koblick's values, so we ship the upstream constants directly.
  * :func:`reproduce_tulip` -- the reproduce-before-trust gate. Propagates the
    sourced IC with the Sundman-regularised propagator (justified because
    tulip orbits graze the secondary at every petal) and reports closure,
    Jacobi, observed petal count, and Floquet multipliers.
  * :func:`petal_count` -- the topological classifier. Counts perilune passages
    over one period; orthogonal to the (k1, k2) winding-number classifier in
    :func:`cyclerfinder.search.binary_star_search.winding_topology`.

NOT in this module: family-switching corrector across bifurcations, L2-NRHO
continuation. Those land in Phase 3 (#266 follow-up).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from cyclerfinder.search.bifurcation_detector import BifurcationPoint
    from cyclerfinder.search.nrho_continuation import SymmetricNRHO
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.core.cr3bp_regularized as creg

# ---------------------------------------------------------------------------
# Koblick / pumpkyn CR3BP constants.
# ---------------------------------------------------------------------------

# Earth-Moon mass parameter mu = m_Moon / (m_Earth + m_Moon).
#
# Source: JPL DE405 planetary constants -- the canonical Koblick 2023 AMOSTECH
# value, also exactly matching Coorbital's pumpkyn `getTulip.m` constants.
# Numerically distinct from the project's default cr3bp_system("Earth", "Moon")
# value (~1.21505860e-2 from gm_de440 system GM) by ~1e-9 -- below the
# tolerance any IC in this module is published to, but kept explicit because
# every tulip IC below was *measured* against this mu and would not close
# perfectly against the project default.
KOBLICK_2023_MU: float = 1.215058560962404e-2

# Earth-Moon characteristic length (Koblick / pumpkyn).  In Koblick's
# normalisation LU is fixed at 389703 km (the JPL DE405 lunar mean-orbit radius
# proxy). The project default cr3bp_system("Earth", "Moon") uses
# SATELLITES["Moon"].sma_km = 384400 km (the more widely cited value); pinning
# the upstream LU here keeps the published tulip ICs dimensionally consistent
# without forcing every caller to re-scale.
KOBLICK_2023_LU_KM: float = 389703.0

# Earth-Moon characteristic time TU = sqrt(LU**3 / G(M_E + M_M)).
# Source: Koblick 2023 -- consistent with LU_KM above and the JPL DE405 GM.
KOBLICK_2023_TU_S: float = 382981.0


def koblick_system() -> cr3bp.CR3BPSystem:
    """Return the CR3BP system in Koblick's normalisation.

    Every IC in :data:`KOBLICK_2023_TABLE4` is referenced to this system; do not
    substitute the project's default Earth-Moon system when reproducing those
    ICs -- the mu differs at the 1e-9 level and the LU/TU differ by ~1.4%.
    """
    return cr3bp.CR3BPSystem(
        mu=KOBLICK_2023_MU,
        primary="Earth",
        secondary="Moon",
        l_km=KOBLICK_2023_LU_KM,
        t_s=KOBLICK_2023_TU_S,
    )


# ---------------------------------------------------------------------------
# Sourced IC table.
# ---------------------------------------------------------------------------

# Tulip-orbit ICs keyed by petal count Np.
#
# Provenance per row:
#   Np = 1 -- NASA Gateway baseline Earth-Moon L2 Southern 9:2 NRHO
#             (see Lee 2019, NASA / e.g. the widely cited seed
#             (x, z, ydot) = (1.0213, -0.1824, -0.1031), period ~ 1.5111 TU
#             also used by tests/search/test_reachable_impulsive.py
#             `_recover_92_nrho` for Zhou error-index validation).
#   Np = 2 -- Coorbital / pumpkyn MATLAB package, ``getTulip.m`` initial
#             conditions for the Np = 2 "butterfly" tulip
#             ``(x0, z0, ydot0, tau0) = (1.01173097566369,
#             0.173908177213147, -0.0799007564597011,
#             2.74865559389775)``. ``tau0`` is the FULL nondimensional period
#             in Koblick's TU. The pumpkyn project distributes these as
#             public-domain MATLAB code; this row is the **independent cross-
#             check anchor** required by the orbit-closure discipline.
#
# Each entry exposes:
#   x0, y0, z0, xdot0, ydot0, zdot0 -- CR3BP nondim initial state
#   T_TU                            -- full nondim period
#   n_petals                        -- the topological petal count (sourced)
#   r_min_km, r_max_km              -- optional perilune / apolune annotations
#                                      (in km); ``None`` when not published.
#
# The Koblick 2023 AMOSTECH Table 4 rows for Np = 3 ... 15 are deliberately
# absent here: the AMOSTECH paper was not accessible during Phase 2 build, and
# the discovery discipline forbids manufacturing ICs (see the project memory
# entries `feedback_orbit_closure_discipline` and
# `feedback_published_rounded_values_are_display`). When the table is digitised
# add rows here verbatim and verify each via the existing reproduce gate.
KOBLICK_2023_TABLE4: dict[int, dict[str, float | int | None]] = {
    1: {
        "x0": 1.0213,
        "y0": 0.0,
        "z0": -0.1824,
        "xdot0": 0.0,
        "ydot0": -0.1031,
        "zdot0": 0.0,
        "T_TU": 1.5111,
        "n_petals": 1,
        "r_min_km": 3000.0,  # NRHO perilune ~3000 km (Lee 2019)
        "r_max_km": None,
        # Sourced as "NASA Gateway baseline 9:2 NRHO seed" -- the seed is sub-
        # display precision and the test gate re-corrects it to a true periodic
        # orbit before drawing topology / Jacobi conclusions.
    },
    2: {
        "x0": 1.01173097566369,
        "y0": 0.0,
        "z0": 0.173908177213147,
        "xdot0": 0.0,
        "ydot0": -0.0799007564597011,
        "zdot0": 0.0,
        "T_TU": 2.74865559389775,
        "n_petals": 2,
        "r_min_km": None,
        "r_max_km": None,
        # Sourced via Coorbital / pumpkyn ``getTulip.m`` (full-double seed).
        # The Koblick 2023 AMOSTECH Table 4 row at fixed x0=1.023731 is a DIFFERENT
        # member of the same Np=2 family (Koblick z0=0.174305, ydot0=-0.082095,
        # tau0=2.756426 -> T=5.512852, nu=1.143759, r_min=252 km). See the variant
        # below (`KOBLICK_2023_TABLE4_PAPER`) — both are kept because they index
        # different points on the SAME Np=2 family curve; the pumpkyn seed sits
        # near a bifurcation (nu~1.0), the paper row at nu=1.144.
    },
}


# Koblick 2023 AMOSTECH Table 4 (page 6): the 15-family IC table at the FIXED
# crossing x0 = 1.023731 (so a one-parameter family indexed by Np). All rows are
# planar perpendicular x-z-plane crossings (y0 = xdot0 = zdot0 = 0); the
# published columns are (x0, z0, ydot0, tau0). #266 Phase 3 verified that
# ``tau0`` IS THE FULL nondim period (the IC closes to machine precision at
# t=tau0 under DOP853 at rtol=1e-12) -- NOT the half-period as initially
# misread. This is consistent with Koblick's Fig 6 captions, which label the
# horizontal axis as "period (TU)" with values in the 1.5-5.8 TU range matching
# the tau0 column directly. Each member has ONE half-period perpendicular
# x-z-plane re-crossing at t = tau0/2 (the orbit's symmetry partner of the IC).
#
# Source: Koblick (2023) "Novel Tulip-Shaped Three-body Orbits for Cislunar Space
# Domain Awareness Missions", AMOSTECH 2023, Poster, Table 4 (p.6). Sourced
# precision: 5-6 sig fig (x), 5-6 (z, ydot), 7 (tau0). The nu column (stability
# index from Koblick Eqn 6) and r_min/r_max (km) carry over from the table.
#
# IMPORTANT: r_min in the paper appears to be a SIGNED quantity (radius vs Moon
# centre, lunar radius = 1737 km), so Np>=3 with r_min < 0 km are LUNAR
# IMPACTORS at this particular IC sample - they are valid bifurcation roots of
# the family curve but unphysical as trajectories. The reproduce gate should not
# call them physical orbits without continuation off the impactor branch.
KOBLICK_2023_TABLE4_PAPER: dict[int, dict[str, float | int]] = {
    1: dict(
        x0=1.023731,
        z0=0.183250,
        ydot0=-0.106950,
        tau0=1.533637,
        n_petals=1,
        nu=1.369020,
        r_min_km=1829.78,
        r_max_km=71031.98,
    ),
    2: dict(
        x0=1.023731,
        z0=0.174305,
        ydot0=-0.082095,
        tau0=2.756426,
        n_petals=2,
        nu=1.143759,
        r_min_km=252.10,
        r_max_km=67614.28,
    ),
    3: dict(
        x0=1.023731,
        z0=0.159022,
        ydot0=-0.049901,
        tau0=3.588824,
        n_petals=3,
        nu=1.000000,
        r_min_km=-752.49,
        r_max_km=61792.20,
    ),
    4: dict(
        x0=1.023731,
        z0=0.138427,
        ydot0=-0.016770,
        tau0=4.050042,
        n_petals=4,
        nu=1.006718,
        r_min_km=-985.88,
        r_max_km=53991.24,
    ),
    5: dict(
        x0=1.023731,
        z0=0.122012,
        ydot0=+0.006413,
        tau0=4.380388,
        n_petals=5,
        nu=1.000000,
        r_min_km=-868.06,
        r_max_km=47824.86,
    ),
    6: dict(
        x0=1.023731,
        z0=0.108984,
        ydot0=+0.026199,
        tau0=4.628880,
        n_petals=6,
        nu=1.000070,
        r_min_km=-652.02,
        r_max_km=42976.89,
    ),
    7: dict(
        x0=1.023731,
        z0=0.098423,
        ydot0=+0.044549,
        tau0=4.827278,
        n_petals=7,
        nu=1.000000,
        r_min_km=-399.67,
        r_max_km=39087.92,
    ),
    8: dict(
        x0=1.023731,
        z0=0.089585,
        ydot0=+0.062522,
        tau0=4.991642,
        n_petals=8,
        nu=1.000000,
        r_min_km=-125.78,
        r_max_km=35870.64,
    ),
    9: dict(
        x0=1.023731,
        z0=0.081981,
        ydot0=+0.080721,
        tau0=5.131974,
        n_petals=9,
        nu=1.000000,
        r_min_km=167.60,
        r_max_km=33137.37,
    ),
    10: dict(
        x0=1.023731,
        z0=0.075272,
        ydot0=+0.099631,
        tau0=5.254942,
        n_petals=10,
        nu=1.000000,
        r_min_km=485.03,
        r_max_km=30759.16,
    ),
    11: dict(
        x0=1.023731,
        z0=0.069206,
        ydot0=+0.119765,
        tau0=5.365409,
        n_petals=11,
        nu=1.000000,
        r_min_km=836.20,
        r_max_km=28641.99,
    ),
    12: dict(
        x0=1.023731,
        z0=0.063570,
        ydot0=+0.141824,
        tau0=5.467385,
        n_petals=12,
        nu=1.000000,
        r_min_km=1238.55,
        r_max_km=26710.27,
    ),
    13: dict(
        x0=1.023731,
        z0=0.058150,
        ydot0=+0.167010,
        tau0=5.564950,
        n_petals=13,
        nu=1.000000,
        r_min_km=1724.83,
        r_max_km=24891.06,
    ),
    14: dict(
        x0=1.023731,
        z0=0.052617,
        ydot0=+0.198057,
        tau0=5.664248,
        n_petals=14,
        nu=1.000000,
        r_min_km=2376.01,
        r_max_km=23082.05,
    ),
    15: dict(
        x0=1.023731,
        z0=0.045796,
        ydot0=+0.247118,
        tau0=5.787296,
        n_petals=15,
        nu=1.000000,
        r_min_km=3566.42,
        r_max_km=20935.32,
    ),
}


# ---------------------------------------------------------------------------
# Reproduce-before-trust result type.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ReproductionResult:
    """Outcome of :func:`reproduce_tulip`.

    Attributes
    ----------
    np_target :
        The petal count requested (the key into :data:`KOBLICK_2023_TABLE4`).
    state0 :
        The sourced IC the orbit was integrated from.
    period :
        The sourced period T_TU (NOT re-fitted; this is the trust gate, not a
        corrector).
    closed :
        True iff ``|state(T) - state(0)|`` is below the closure tolerance.
    closure_residual :
        Raw L2 closure residual at one period.
    jacobi_constant :
        ``cr3bp.jacobi_constant(state0, mu)`` -- the Jacobi value at the IC.
    n_petals_observed :
        Integer petal count derived from :func:`petal_count`. A matching
        ``n_petals_observed == np_target`` is the topological half of the gate.
    monodromy_eigs :
        The six Floquet multipliers (sorted by magnitude descending). Computed
        from the STM at one full period via the standard variational EOM.
    """

    np_target: int
    state0: NDArray[np.float64]
    period: float
    closed: bool
    closure_residual: float
    jacobi_constant: float
    n_petals_observed: int
    monodromy_eigs: NDArray[np.complex128]


# ---------------------------------------------------------------------------
# Petal-count topological classifier.
# ---------------------------------------------------------------------------


def petal_count(
    state0: NDArray[np.float64],
    period: float,
    system: cr3bp.CR3BPSystem,
    *,
    rtol: float = 1e-11,
    atol: float = 1e-11,
) -> int:
    """Count secondary-relative perilune passages over one period.

    A petal is one local minimum of ``r2(t)`` (the spacecraft's distance to the
    secondary) over the period ``[0, T]``. Perilune passages are identified by
    sign changes of the secondary-relative radial velocity ``dr2/dt`` -- a
    zero-crossing of ``r * vr_to_secondary`` from negative to positive marks a
    minimum.

    The classifier runs in **physical time** (not regularised time): perilune
    counts are unambiguous in physical time and ``solve_ivp`` event detection
    in physical time is the simplest correct implementation. For a low-perilune
    orbit the integrator may take many steps near each petal, but the event
    detection is cheap relative to the petal cost itself.

    Parameters
    ----------
    state0, period, system :
        IC, full nondim period, and CR3BP system the IC is referenced to.

    Returns
    -------
    int :
        The integer petal count over ``[0, T]``. The endpoint ``t = T`` is
        excluded from the count: if the orbit's periodicity puts a perilune
        exactly at t = T (which it does for a sourced IC at the family node),
        it counts the t = 0 perilune (which equals t = T modulo periodicity)
        rather than double-counting both.

    Raises
    ------
    RuntimeError
        If the underlying integrator fails. The classifier integrates in
        physical time without regularisation; deep low perilune may stress the
        solver. In that case fall back to the regularised propagator and
        manually extract the perilune count from the dense ``state_at_s`` grid.
    """
    s0 = np.asarray(state0, dtype=np.float64)
    mu = system.mu
    moon_x = 1.0 - mu

    # Event: dr2/dt = 0, direction +1 (rising), i.e. a local minimum of r2.
    # Equivalent (and slightly cheaper) form: the sign of the secondary-relative
    # radial velocity ``rx * vx + ry * vy + rz * vz``.
    def _peri_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        rx = y[0] - moon_x
        ry = y[1]
        rz = y[2]
        vx, vy, vz = y[3], y[4], y[5]
        return float(rx * vx + ry * vy + rz * vz)

    _peri_event.terminal = False  # type: ignore[attr-defined]
    _peri_event.direction = 1.0  # type: ignore[attr-defined]

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        s0,
        args=(mu,),  # type: ignore[call-overload]
        method="DOP853",
        rtol=rtol,
        atol=atol,
        events=_peri_event,
    )
    if not sol.success:
        raise RuntimeError(f"petal_count: integrator failed at t={sol.t[-1]}: {sol.message}")
    t_events = sol.t_events[0] if sol.t_events is not None else np.array([])

    # Exclude the endpoint t = T (within solver tolerance): if the orbit is
    # truly periodic with a perilune at t = 0, then t = T is the same petal
    # modulo periodicity, not an additional one.
    eps = max(1e-9, 1e-9 * abs(period))
    n = int(np.sum((t_events > eps) & (t_events < period - eps)))

    # If t = 0 is itself a perilune (the orbit IC sits at perilune, which is
    # standard for the NRHO Np=1 seed -- not the case for the pumpkyn Np=2
    # tulip IC, which sits at apolune), count it as the first petal explicitly.
    # Heuristic: |dr2/dt(0)| ~ 0 AND d2(r2)/dt2(0) > 0 (true local min, not max).
    rx0, ry0, rz0 = s0[0] - moon_x, s0[1], s0[2]
    vx0, vy0, vz0 = s0[3], s0[4], s0[5]
    r_dot_v0 = rx0 * vx0 + ry0 * vy0 + rz0 * vz0
    r0 = float(np.sqrt(rx0 * rx0 + ry0 * ry0 + rz0 * rz0))
    if r0 > 0 and abs(r_dot_v0) / r0 < 1e-4:
        # Sign of d/dt(dr2/dt) = d/dt(r dot v) at t=0:
        #   d/dt(r dot v) = |v|^2 + r dot a
        # so a true minimum has |v|^2 + r dot a > 0.
        ax0, ay0, az0 = cr3bp.cr3bp_eom(0.0, s0, mu)[3:6]
        v2 = vx0 * vx0 + vy0 * vy0 + vz0 * vz0
        r_dot_a = rx0 * ax0 + ry0 * ay0 + rz0 * az0
        if v2 + r_dot_a > 0.0:
            n += 1
    return n


# ---------------------------------------------------------------------------
# 3D-topology gate (#322 bug-fix).
# ---------------------------------------------------------------------------

# Out-of-plane amplitude floor for the "genuine 3D tulip" gate.
#
# Bug context: ``petal_count`` counts local minima of the spacecraft-to-secondary
# distance in 3D, but a planar Np-petal orbit (z(t) ≡ 0) still produces exactly
# Np in-plane perilune minima — the classifier alone CANNOT distinguish a
# planar Np-petal orbit from a genuine 3D Np-tulip. At very small mu (e.g. the
# Mars-Phobos system, mu ~ 1.65e-8), the symmetric corrector drives the seed's
# z0 toward zero and the orbit collapses to planar; petal_count then fires a
# FALSE POSITIVE on the tulip topology gate. See #313 negative + #322 fix.
#
# Sourced threshold rationale:
# ``KOBLICK_2023_TABLE4_PAPER`` publishes z0 values for Np=1..15 in the Koblick
# Earth-Moon normalisation. The SMALLEST z0 in that table is the Np=15 row at
# z0 = 0.045796 nondim; the LARGEST is the Np=1 row at z0 = 0.183250 nondim. We
# use a conservative floor of ``5e-3`` nondim (~5% of the smallest Koblick z0
# at Np=15, and ~3% of the canonical Np=2 z0=0.174305). This is a Phase 1
# floor: anything below 5e-3 nondim is definitely planar-collapse rather than
# a 3D tulip on the published family. The number is a conservative cliff
# (sourced family minimum is 0.046 -- two decades above the floor), not a
# precision threshold. Phase 2 (#322 follow-up) may refine this against
# Koblick's Fig.~7 out-of-plane envelopes if/when they are digitised.
#
# Concretely at Earth-Moon scale, 5e-3 nondim ≈ 1.95e3 km (LU=389703 km, so
# ~1948 km out-of-plane). At Mars-Phobos (LU = 9375 km) the same floor is
# 46.9 km; at Pluto-Charon (LU = 19600 km) it is 98 km.
TULIP_Z_AMPLITUDE_FLOOR_NONDIM: float = 5e-3


def out_of_plane_amplitude(
    state0: NDArray[np.float64],
    period: float,
    system: cr3bp.CR3BPSystem,
    *,
    rtol: float = 1e-11,
    atol: float = 1e-11,
    n_samples: int = 401,
) -> float:
    """Compute ``max |z(t)|`` over one period of a CR3BP orbit.

    Used by :func:`is_three_dimensional` and the ``find_tulip_at_system``
    topology gate to distinguish a genuine 3D tulip from a planar Np-petal
    orbit that happens to share the in-plane petal count.

    The integrator runs in physical time with ``DOP853`` and a dense ``t_eval``
    grid (the cost is one cheap forward propagation per gate call). For
    near-collision low-perilune orbits the regularised propagator might be
    safer; in practice the Koblick / Mars-moon ICs the gate has to discriminate
    do not hit collisions during this evaluation.

    Parameters
    ----------
    state0, period, system :
        IC, full nondim period, and CR3BP system the IC is referenced to.
    rtol, atol :
        Integrator tolerances. Loosened from the corrector's ``1e-12`` because
        we want ``max|z|`` to a few significant figures only.
    n_samples :
        Number of evenly-spaced sample points on ``[0, T]``. Default 401 is
        cheap and captures any plausible z-extremum even for the Np=15 row
        (one out-of-plane lobe per ~6.3 TU / 15 petals ~ 0.4 TU spacing).

    Returns
    -------
    float :
        The maximum of ``|z(t)|`` over the dense sample grid.

    Raises
    ------
    RuntimeError
        If the integrator fails.
    """
    from scipy.integrate import solve_ivp

    s0 = np.asarray(state0, dtype=np.float64)
    t_eval = np.linspace(0.0, period, int(n_samples))
    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        s0,
        args=(system.mu,),
        method="DOP853",
        rtol=rtol,
        atol=atol,
        t_eval=t_eval,
    )
    if not sol.success:
        raise RuntimeError(
            f"out_of_plane_amplitude: integrator failed at t={sol.t[-1]}: {sol.message}"
        )
    return float(np.max(np.abs(sol.y[2])))


def is_three_dimensional(
    state0: NDArray[np.float64],
    period: float,
    system: cr3bp.CR3BPSystem,
    *,
    z_floor: float = TULIP_Z_AMPLITUDE_FLOOR_NONDIM,
    rtol: float = 1e-11,
    atol: float = 1e-11,
) -> tuple[bool, float]:
    """Return ``(is_3d, max_abs_z)`` for the 3D-topology gate.

    An orbit qualifies as 3D iff BOTH ``|z0| >= z_floor`` AND
    ``max|z(t)| >= z_floor``. Either-or would be insufficient: a sourced
    Koblick IC can have z0 near zero at the perpendicular-crossing IC instant
    while still having a non-trivial out-of-plane lobe (so we check max|z|);
    and a degenerate IC that crosses z=0 only at sample times could pass a
    pure max|z| check spuriously (so we also pin z0).

    See :data:`TULIP_Z_AMPLITUDE_FLOOR_NONDIM` for the floor's sourced
    justification and units.
    """
    z0_ok = abs(float(state0[2])) >= z_floor
    if not z0_ok:
        # Short-circuit: the IC itself is planar, no need to integrate.
        return False, abs(float(state0[2]))
    max_abs_z = out_of_plane_amplitude(state0, period, system, rtol=rtol, atol=atol)
    return (max_abs_z >= z_floor), max_abs_z


TopologyVerdict = Literal[
    "3D tulip",
    "planar Np-petal collapse",
    "petal count mismatch",
    "seed no converge",
    "unknown",
]


# ---------------------------------------------------------------------------
# Reproduce-before-trust gate.
# ---------------------------------------------------------------------------


def _state0_from_row(row: dict[str, float | int | None]) -> NDArray[np.float64]:
    return np.array(
        [
            float(row["x0"]),  # type: ignore[arg-type]
            float(row["y0"]),  # type: ignore[arg-type]
            float(row["z0"]),  # type: ignore[arg-type]
            float(row["xdot0"]),  # type: ignore[arg-type]
            float(row["ydot0"]),  # type: ignore[arg-type]
            float(row["zdot0"]),  # type: ignore[arg-type]
        ],
        dtype=np.float64,
    )


def _monodromy_at_period(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    rtol: float,
    atol: float,
) -> NDArray[np.float64]:
    """Integrate the variational EOM for one full period; return the STM."""
    arc = cr3bp.propagate(system, state0, period, with_stm=True, rtol=rtol, atol=atol)
    assert arc.stm is not None
    return arc.stm


def reproduce_tulip(
    np_target: int,
    *,
    n_periods: int = 1,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    closure_tol: float = 1e-4,
) -> ReproductionResult:
    """Propagate a sourced tulip IC and return the reproduce-gate evidence.

    Parameters
    ----------
    np_target :
        Petal count -- key into :data:`KOBLICK_2023_TABLE4`. ``1`` is the
        parent NRHO, ``2`` is the butterfly (the primary reproduce target).
    n_periods :
        Number of full periods to propagate. ``1`` is the gate; ``>= 2`` is
        useful for visual stability checks.
    rtol, atol :
        Integrator tolerances forwarded to the propagators.
    closure_tol :
        L2 residual under which the orbit is declared "closed". The default
        ``1e-4`` matches the **published precision** of the sourced ICs --
        Koblick's 5-6 sig fig and pumpkyn's full-double-but-itself-corrected
        values both leave residual floors above machine precision. **Do not
        tighten this without first re-correcting the IC**; tightening to
        ``1e-12`` and watching the test fail is a circular validation
        (see the project memory ``feedback_orbit_closure_discipline``).

    Notes
    -----
    The state is propagated with the Sundman-regularised propagator (#266
    Phase 1) because tulip orbits graze the secondary at every petal; the
    standard propagator can take many more RHS evaluations at the same
    tolerance. Closure residual is computed at physical time = ``period``
    via the regularised propagator's ``t_stop`` terminal-time event.

    Raises
    ------
    KeyError
        If ``np_target`` is not present in :data:`KOBLICK_2023_TABLE4`.
    """
    if np_target not in KOBLICK_2023_TABLE4:
        raise KeyError(
            f"reproduce_tulip: np_target={np_target} is not in KOBLICK_2023_TABLE4; "
            f"available: {sorted(KOBLICK_2023_TABLE4.keys())}"
        )
    if n_periods < 1:
        raise ValueError(f"reproduce_tulip: n_periods must be >= 1, got {n_periods}")
    row = KOBLICK_2023_TABLE4[np_target]
    state0 = _state0_from_row(row)
    period = float(row["T_TU"])  # type: ignore[arg-type]
    system = koblick_system()

    # Propagate over n_periods * T via the regularised propagator with an
    # event-locked terminal physical time. The regularised propagator handles
    # the low-perilune grazes that would otherwise burn integrator steps.
    t_final = n_periods * period
    s_span = creg.physical_to_regularized_span(system, state0, (0.0, t_final))
    s_span = (s_span[0], s_span[1] * 1.5)
    arc = creg.propagate_regularized(
        system,
        state0,
        s_span,
        rtol=rtol,
        atol=atol,
        regularization="r2",  # Moon-only regularisation is cheapest at lunar perilune
        t_stop=t_final,
    )
    state_f = arc.state_at_s[:, -1]
    closure_residual = float(np.linalg.norm(state_f - state0))
    closed = closure_residual < closure_tol

    jacobi = cr3bp.jacobi_constant(state0, system.mu)

    n_petals_observed = petal_count(state0, period, system, rtol=rtol, atol=atol)

    stm = _monodromy_at_period(system, state0, period, rtol=rtol, atol=atol)
    eigs = np.linalg.eigvals(stm)
    # Sort by magnitude descending; complex.
    order = np.argsort(-np.abs(eigs))
    eigs_sorted = np.asarray(eigs[order], dtype=np.complex128)

    return ReproductionResult(
        np_target=int(np_target),
        state0=state0,
        period=period,
        closed=closed,
        closure_residual=closure_residual,
        jacobi_constant=jacobi,
        n_petals_observed=int(n_petals_observed),
        monodromy_eigs=eigs_sorted,
    )


# ---------------------------------------------------------------------------
# End-to-end Phase 3 reproduce gate: continuation + family-switching.
# ---------------------------------------------------------------------------


def find_tulip_via_continuation(
    np_target: int = 2,
    *,
    system: cr3bp.CR3BPSystem | None = None,
    d_x0: float = 2e-3,
    n_steps_max: int = 60,
    perilune_floor_km: float | None = None,
    eigenvector_step: float = 1e-3,
    tol: float = 1e-10,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    multi_shooting: bool = False,
    multi_shoot_segments: int | None = None,
    seed_row: dict[str, float | int] | None = None,
) -> FindTulipResult:
    """End-to-end Phase 3 reproduce gate: continuation + family-switching.

    1. From the Koblick 2023 AMOSTECH Table 4 Np=1 (paper-row) seed at
       ``x0=1.023731``, correct it via the symmetric corrector.
    2. Continue the family in ``x0`` (direction=-1: decreasing x0) until a
       period-doubling (k=2) bifurcation bracket is found OR the perilune floor
       (if any) is reached.
    3. At the first k=2 bracket, call
       :func:`cyclerfinder.genome.family_switch.switch_family` with ``k=2``.
    4. Independently verify ``petal_count == 2`` on the switched orbit.

    Parameters
    ----------
    np_target :
        The petal count to land on. Currently only ``np_target=2`` is wired in.
    system :
        CR3BP system. Defaults to the Koblick / pumpkyn normalisation.
    d_x0 :
        Continuation step in x0.
    n_steps_max :
        Step budget for the continuation.
    perilune_floor_km :
        If given, stop the continuation when perilune drops below this floor.
    eigenvector_step :
        Forwarded to :func:`switch_family`.
    tol, rtol, atol :
        Forwarded to the corrector and propagators.
    multi_shooting :
        If True, the family-switch step uses the Phase 4 multi-shooter
        instead of single-shooting. Default ``False`` keeps Phase 3 behaviour.
    multi_shoot_segments :
        Number of segments for the multi-shooter. ``None`` selects
        ``max(np_target, 2)`` inside :func:`switch_family`.
    seed_row :
        Explicit seed row (``x0, z0, ydot0, tau0`` keys). When ``None`` (default)
        the Koblick AMOSTECH Table 4 Np=1 row is used. Pass a custom seed when
        running at a non-Earth-Moon system where the Koblick IC may not be a
        periodic orbit -- callers obtain a system-appropriate seed via
        :func:`find_tulip_at_system`.

    Returns
    -------
    FindTulipResult :
        Carries the seed orbit, the continuation branch, the first matched
        bifurcation, the switched member (or ``None``), and a high-level
        ``success`` flag.

    Notes
    -----
    This function does NOT modify the catalogue or write to disk. It is a
    diagnostic / reproduce-gate routine for the genome layer.
    """
    # Local imports to avoid a circular dependency at module-load time
    # (family_switch imports genome.tulip).
    from cyclerfinder.genome.family_switch import switch_family
    from cyclerfinder.search.nrho_continuation import (
        continue_nrho_family,
        correct_symmetric_nrho,
    )

    if np_target != 2:
        raise NotImplementedError(
            f"find_tulip_via_continuation: only np_target=2 wired in Phase 3, got {np_target}"
        )
    if system is None:
        system = koblick_system()
    if seed_row is None:
        seed_row = KOBLICK_2023_TABLE4_PAPER[1]
    # tau0 is the FULL period for the NRHO family rows (verified by closure of
    # the IC at t=tau0 to machine precision; the half-period reading of the
    # Phase 2 docstring was incorrect).
    t_seed = float(seed_row["tau0"])
    seed = correct_symmetric_nrho(
        system,
        float(seed_row["x0"]),
        float(seed_row["z0"]),
        float(seed_row["ydot0"]),
        t_seed,
        tol=1e-11,
        rtol=rtol,
        atol=atol,
    )
    if not seed.converged:
        return FindTulipResult(
            seed=seed,
            branch_members=[],
            bifurcation=None,
            switched=None,
            success=False,
            reason="seed_no_converge",
        )

    # The period-doubling bifurcation on this family is a TANGENT bifurcation:
    # the real hyperbolic pair (-2.30, -0.43) collide AT -1 as x0 decreases,
    # then turn into a complex pair on the unit circle. The "distance to -1"
    # signal has a TANGENT MINIMUM at the bifurcation (rather than a crossing
    # through zero), so a coarse tolerance like bif_tol=0.1 brackets it.
    branch = continue_nrho_family(
        seed,
        system,
        label="koblick_l2_southern_nrho",
        direction=-1,
        d_x0=d_x0,
        n_steps_max=n_steps_max,
        perilune_floor_km=perilune_floor_km,
        tol=1e-10,
        rtol=rtol,
        atol=atol,
        bif_k_max=4,
        bif_tol=1e-1,
        stop_on_first_bifurcation=False,
    )
    k2_bifs = [b for b in branch.bifurcations if b.k == 2]
    if not k2_bifs:
        return FindTulipResult(
            seed=seed,
            branch_members=list(branch.members),
            bifurcation=None,
            switched=None,
            success=False,
            reason=f"no_k2_bifurcation:{branch.stop_reason.value}:{branch.n_steps}steps",
        )
    # Pick the bracket whose MIDPOINT member has the smallest |lam + 1|: this is
    # the closest-to-bifurcation member on the family, the best parent for the
    # family-switching corrector.
    from cyclerfinder.search.bifurcation_detector import floquet_multipliers

    best_idx = -1
    best_dist = float("inf")
    for i, m in enumerate(branch.members):
        if m.monodromy is None:
            continue
        eigs = floquet_multipliers(m.monodromy)
        d = min(abs(complex(e) + 1.0) for e in eigs)
        if d < best_dist:
            best_dist = d
            best_idx = i
    if best_idx < 0:
        return FindTulipResult(
            seed=seed,
            branch_members=list(branch.members),
            bifurcation=k2_bifs[0],
            switched=None,
            success=False,
            reason="no_monodromy_on_any_member",
        )
    parent = branch.members[best_idx]
    bif = k2_bifs[0]
    switched = switch_family(
        parent,
        bif,
        system,
        k=2,
        eigenvector_step=eigenvector_step,
        tol=tol,
        rtol=rtol,
        atol=atol,
        multi_shooting=multi_shooting,
        multi_shoot_segments=multi_shoot_segments,
    )
    return FindTulipResult(
        seed=seed,
        branch_members=list(branch.members),
        bifurcation=bif,
        switched=switched,
        success=switched is not None,
        reason=("ok" if switched is not None else "family_switch_no_converge"),
    )


@dataclass(frozen=True)
class FindTulipResult:
    """Outcome of :func:`find_tulip_via_continuation`.

    Attributes
    ----------
    topology_verdict :
        3D-topology classification of the returned orbit (#322). One of:

          * ``"3D tulip"`` -- genuine 3D Np-tulip (passes both petal_count and
            out-of-plane amplitude gates).
          * ``"planar Np-petal collapse"`` -- petal_count matches but the
            orbit's z(t) collapsed below
            :data:`TULIP_Z_AMPLITUDE_FLOOR_NONDIM`; NOT a real tulip.
          * ``"petal count mismatch"`` -- corrected orbit has the wrong petal
            count (e.g. seed at Np=1 family, target Np=2 -- the canonical
            Tier B fallback signal).
          * ``"seed no converge"`` -- the corrector failed before topology
            could be classified.
          * ``"unknown"`` -- topology not classified (default for paths that
            pre-date the gate, kept for backward compat).
    max_abs_z :
        ``max |z(t)|`` over one period of the returned orbit (nondim), or
        ``None`` if topology was not classified. Carries the 3D-topology
        evidence alongside the verdict.
    """

    seed: SymmetricNRHO
    branch_members: list[SymmetricNRHO]
    bifurcation: BifurcationPoint | None
    switched: SymmetricNRHO | None
    success: bool
    reason: str
    topology_verdict: TopologyVerdict = "unknown"
    max_abs_z: float | None = None


def find_tulip_at_system(
    system: cr3bp.CR3BPSystem,
    *,
    np_target: int = 2,
    multi_shooting: bool = True,
    multi_shoot_segments: int | None = None,
    seed_row: dict[str, float | int] | None = None,
    d_x0: float = 5e-4,
    n_steps_max: int = 60,
    perilune_floor_km: float | None = None,
    eigenvector_step: float = 1e-2,
    tol: float = 1e-9,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    try_direct_seed: bool = True,
) -> FindTulipResult | None:
    """Cross-system tulip-orbit finder (#268 Phase 4 + #280 Phase 5).

    Find a Np=``np_target`` tulip at an arbitrary CR3BP system. The routine
    walks two escalation tiers in order:

    **Tier A (Phase 5 #280, ``try_direct_seed=True``).**  Correct the seed at
    the target system's mu and INDEPENDENTLY verify the topology via
    :func:`petal_count`. The symmetric corrector's basin of attraction shifts
    with mu, so a Koblick Earth-Moon NRHO seed CAN converge to a
    *system-native* Np=``np_target`` orbit on a different (typically planar)
    family at the target mu (verified empirically at Jupiter-Ganymede,
    Saturn-Titan, and Pluto-Charon -- the Koblick Np=1 seed lands Np=2 there).
    When the corrected seed's petal count already matches ``np_target``, return
    it directly: it IS the answer. No continuation, no family-switching needed.

    **Tier B (Phase 3/4 fallback).**  When the direct-seed petal count does NOT
    match (the canonical Earth-Moon path: Koblick seed -> Np=1 NRHO -> continue
    in x0 -> k=2 bifurcation -> family-switch -> Np=2), fall through to
    :func:`find_tulip_via_continuation`. The Phase 4 multi-shooter
    handles strong period-doublings where single-shooting fails.

    Parameters
    ----------
    system :
        Any CR3BP system (built via
        :func:`cyclerfinder.core.cr3bp.cr3bp_system` or the system-specific
        :func:`koblick_system`).
    np_target :
        Petal count to land on. Currently only ``np_target=2`` is supported by
        Tier B's continuation routine; Tier A in principle accepts any
        ``np_target`` (whatever the system-native seed-correction produces).
    multi_shooting :
        Tier B only. If True (default), the family-switch step uses
        multi-shooting. Ignored in Tier A (no family-switch).
    multi_shoot_segments :
        Tier B only. Number of multi-shooting segments. ``None`` picks
        ``max(np_target, 2)``.
    seed_row :
        Explicit IC seed for the NRHO -- a dict with keys ``x0, z0, ydot0,
        tau0``. When ``None`` (default) the Koblick AMOSTECH Table 4 Np=1
        paper-row is used. Used by BOTH tiers.
    d_x0, n_steps_max, perilune_floor_km :
        Tier B continuation parameters forwarded to
        :func:`continue_nrho_family`.
    eigenvector_step :
        Tier B family-switch parameter.
    tol, rtol, atol :
        Forwarded to the corrector and propagators in both tiers.
    try_direct_seed :
        When True (default), Tier A runs first. Set False to force the Tier B
        Phase 3/4 path (useful for testing the bifurcation pipeline at systems
        where Tier A also works).

    Returns
    -------
    FindTulipResult | None :
        ``None`` only when the seed orbit fails to converge at the given
        system (the corrector cannot find ANY periodic orbit near the seed).
        Otherwise a :class:`FindTulipResult` carrying the seed, branch,
        bifurcation, switched member (or ``None``), and ``success`` flag.
        When Tier A succeeds, the seed itself IS the discovered orbit:
        ``switched is seed``, ``branch_members=[]``, ``bifurcation=None``,
        and ``reason="direct_seed_match"``.

    Notes
    -----
    No catalogue writeback. This is a discovery / cross-system reproduction
    routine; the caller (e.g. :mod:`scripts.tulip_discovery_probe`) decides
    whether to log the outcome.

    Independent cross-check (orbit-closure discipline). Tier A's "the seed is
    the answer" claim is gated on petal_count, which is INDEPENDENT of the
    corrector's residual: the corrector only ensures perpendicular-crossing
    closure, while petal_count is a topological classifier (counts perilune
    passes). The two agreeing on np_target IS the cross-check; this is the
    same petal_count gate that ``switch_family(..., verify_petal_count=True)``
    uses in Tier B.
    """
    if try_direct_seed:
        from cyclerfinder.search.nrho_continuation import correct_symmetric_nrho

        active_seed_row: dict[str, float | int] = (
            seed_row if seed_row is not None else KOBLICK_2023_TABLE4_PAPER[1]
        )
        t_seed = float(active_seed_row["tau0"])
        direct_seed = correct_symmetric_nrho(
            system,
            float(active_seed_row["x0"]),
            float(active_seed_row["z0"]),
            float(active_seed_row["ydot0"]),
            t_seed,
            tol=1e-11,
            rtol=rtol,
            atol=atol,
        )
        if direct_seed.converged:
            state0 = np.array(
                [direct_seed.x0, 0.0, direct_seed.z0, 0.0, direct_seed.ydot0, 0.0],
                dtype=np.float64,
            )
            try:
                n_direct = petal_count(state0, direct_seed.T_TU, system, rtol=rtol, atol=atol)
            except RuntimeError:
                n_direct = -1
            if n_direct == np_target:
                # #322 fix: petal_count alone is NOT sufficient -- a planar
                # Np-petal orbit (z(t) ≡ 0) shares the in-plane petal count
                # with a genuine 3D tulip. At very small mu the symmetric
                # corrector collapses z0 -> 0; we must independently verify
                # 3D topology via out-of-plane amplitude.
                try:
                    is_3d, max_abs_z = is_three_dimensional(
                        state0, direct_seed.T_TU, system, rtol=rtol, atol=atol
                    )
                except RuntimeError:
                    is_3d, max_abs_z = False, float("nan")
                if is_3d:
                    # Tier A hit: the seed IS already the target Np tulip at this mu.
                    return FindTulipResult(
                        seed=direct_seed,
                        branch_members=[],
                        bifurcation=None,
                        switched=direct_seed,
                        success=True,
                        reason="direct_seed_match",
                        topology_verdict="3D tulip",
                        max_abs_z=max_abs_z,
                    )
                # Petal count matched but the orbit collapsed to planar -- this
                # is the #322 false-positive pattern. Refuse to claim success;
                # fall through to Tier B (which uses bifurcation tracking +
                # family-switching, immune to the z0-collapse mode).
                # Carry the diagnostic verdict so the caller can log it.
                tier_a_collapse_result: FindTulipResult | None = FindTulipResult(
                    seed=direct_seed,
                    branch_members=[],
                    bifurcation=None,
                    switched=None,
                    success=False,
                    reason="direct_seed_planar_collapse",
                    topology_verdict="planar Np-petal collapse",
                    max_abs_z=max_abs_z,
                )
            else:
                tier_a_collapse_result = None
        else:
            tier_a_collapse_result = None
    else:
        tier_a_collapse_result = None

    # Tier B fallback: full Phase 3/4 pipeline (continuation + family-switch).
    try:
        tier_b = find_tulip_via_continuation(
            np_target=np_target,
            system=system,
            d_x0=d_x0,
            n_steps_max=n_steps_max,
            perilune_floor_km=perilune_floor_km,
            eigenvector_step=eigenvector_step,
            tol=tol,
            rtol=rtol,
            atol=atol,
            multi_shooting=multi_shooting,
            multi_shoot_segments=multi_shoot_segments,
            seed_row=seed_row,
        )
    except NotImplementedError:
        # Tier B currently only handles np_target=2. For higher Np targets
        # where Tier A flagged a planar collapse, the caller still benefits
        # from the diagnostic; surface the Tier A verdict.
        if tier_a_collapse_result is not None:
            return tier_a_collapse_result
        raise

    # If Tier B succeeded, classify the topology of the switched member too.
    if tier_b.success and tier_b.switched is not None:
        sw = tier_b.switched
        sw_state0 = np.array([sw.x0, 0.0, sw.z0, 0.0, sw.ydot0, 0.0], dtype=np.float64)
        try:
            is_3d_b, max_abs_z_b = is_three_dimensional(
                sw_state0, sw.T_TU, system, rtol=rtol, atol=atol
            )
        except RuntimeError:
            is_3d_b, max_abs_z_b = False, float("nan")
        verdict_b: TopologyVerdict = "3D tulip" if is_3d_b else "planar Np-petal collapse"
        # Re-emit with topology fields populated.
        return FindTulipResult(
            seed=tier_b.seed,
            branch_members=tier_b.branch_members,
            bifurcation=tier_b.bifurcation,
            switched=tier_b.switched if is_3d_b else None,
            success=is_3d_b,
            reason=(tier_b.reason if is_3d_b else "tier_b_planar_collapse"),
            topology_verdict=verdict_b,
            max_abs_z=max_abs_z_b,
        )
    # Tier B did not produce a switched orbit. Prefer the Tier A diagnostic
    # if we have one (it's more specific than "no_k2_bifurcation").
    if tier_a_collapse_result is not None:
        return tier_a_collapse_result
    return tier_b
