"""Seedless spectral (harmonic-balance) periodic-orbit discovery for the QBCP (#611).

Follow-up to `#606`'s `variational_periodic_orbit.py`, which built the CR3BP-
specific version of this idea and used it to cross `#556`'s L1 quasi-halo
wall. This module asks whether the SAME harmonic-balance philosophy attacks
`#538`/`#544`'s much harder wall: the QBCP (Quasi-Bicircular Restricted
Four-Body Problem) Earth-Moon L1/L2 region, which `#544` root-caused as
"violently unstable" -- a frozen-time linearization rate of ~2-3 over one
Sun-synodic period `T_s ~ 6.79` implies a one-period amplification of
``exp(rate * T_s) ~ 1e6-1e8``, which is why every SHOOTING-based corrector in
this codebase (single-period GMOS, and even naive forward propagation of the
published golden ICs themselves) blows up here.

Does the CR3BP module's approach generalize by simply swapping in the QBCP
EOM residual? NO -- investigated and found genuinely different, not just
assumed:

1. **Non-autonomous with a FIXED, known period, not a free unknown.** The
   CR3BP module's Fourier series is built over an arbitrary loop phase
   ``theta = 2*pi*t/period`` with ``period`` itself a free unknown (``logT``
   in the free-variable vector) -- exploiting the CR3BP's autonomous
   time-shift symmetry (any phase of a periodic orbit is an equally valid
   representative, so a phase GAUGE must be pinned to remove that
   degeneracy). The QBCP's ``qbcp_eom`` depends EXPLICITLY on absolute time
   ``t`` through ``evaluate_alphas`` (the Sun's synodic motion) with a fixed
   epoch convention (``theta_sun0 = 0``) -- shifting a candidate trajectory in
   time does NOT produce another valid solution unless the epoch is shifted
   with it, so there is no continuous time-shift symmetry to gauge-fix.
   Instead, only a discrete family of period choices is physically
   meaningful: the orbit must close after exactly one Sun-synodic period
   ``T_s`` (or an integer multiple, for subharmonic/resonant substitutes) --
   this is precisely what the published POL1/POL2 "dynamical substitutes"
   are (Rosales & Jorba 2023). So the period becomes a FIXED input
   (``period_multiple * T_s``), not a solved-for unknown, and the CR3BP
   module's phase-gauge machinery (fixing ``sin_x1 = 0``) has no analogue
   here and is dropped entirely -- not carried over.
2. **A first-order 6-state canonical system, not a second-order 3-position
   system.** The CR3BP module Fourier-expands ONLY (x, y, z) and derives
   velocity/acceleration analytically via theta-derivatives of that same
   series, because the CR3BP rotating-frame EOM is exactly second-order in
   position. The QBCP's canonical (PM) representation is a genuinely
   first-order system in SIX state variables (x, y, z, px, py, pz) --
   ``px``/``py``/``pz`` are conjugate momenta linearly mixed with position
   through the time-varying ``alpha_1``/``alpha_2``/``alpha_3`` coefficients
   (see ``core.qbcp.qbcp_eom``), not simply the position series'
   derivative. Reusing the CR3BP residual code unchanged is not possible;
   this module independently Fourier-expands all SIX canonical state
   components and matches FIRST time-derivatives against ``qbcp_eom``'s six
   right-hand sides at each collocation point.

Net: this is a genuine, bounded ADAPTATION (same harmonic-balance/collocation
principle, same ``scipy.optimize.least_squares`` machinery, no new
dependency), not a drop-in EOM swap. It is also, in one respect, SIMPLER than
the CR3BP case: fixing the period removes both the ``logT`` free variable and
the phase-gauge convention, and empirically (see
``tests/search/test_variational_periodic_orbit_qbcp.py``) the QBCP has no
analogue of the CR3BP's "trivial zero-amplitude equilibrium" attractor for
the optimizer to collapse onto (the Sun's forcing is never zero, so there is
no autonomous fixed point at all) -- no amplitude anchor is required.

Positive control and honest scope boundary
-------------------------------------------
The positive control target is the published POL1/POL2 golden (Rosales &
Jorba 2023 Table 4, already used as the sourced reference in
``tests/core/test_qbcp.py`` and independently reconstructed from scratch via
a 12-segment multiple-shooting corrector in
``scripts/analyze_593_qbcp_l1_substitute_reconciliation.py``, converging to
periodicity residual ~1e-11 and landing ``dist-to-POL1 = 1.81e-2`` under the
current, `#592`-fixed ``alpha_6`` scaling -- a gap independently attributed in
`#544`'s own investigation to a Gimeno-2018-vs-Rosales-2023 model-instance
Fourier-refit difference, NOT a corrector defect). This module, cold-started
with NO continuation bootstrap at all (no CR3BP L1 fixed point, no BCR4BP
mu_sun ramp -- just a rough center guess near the collinear point and the
KNOWN period ``T_s``), converges directly to the same object, landing at an
almost identical ``dist-to-POL1`` (see the test file for the exact number)
-- an independent cross-check of the multi-shooting result AND a genuine
wall-crossing: `#544`'s whole multi-stage bootstrap chain existed specifically
because a single-shot/single-period method could not survive this region's
instability, and this method reaches the same answer from a plain cold start
in a fraction of a second.

**What this module does NOT do (explicit scope boundary, not overclaimed):**
`#538`/`#544`'s actual named target is the QBCP EM-L1/L2 invariant 2-TORUS
(``genome.qbcp_torus.correct_qbcp_torus``) -- a genuinely quasi-periodic,
two-angle family (longitudinal + transverse Fourier modes, a free rotation
number, a stroboscopic-map formulation) that has the periodic orbit found
here as its zero-amplitude CENTER, not as itself. Generalizing this harmonic-
balance approach to a genuine 2D quasi-periodic torus (an extra angle
variable, an unknown rotation number entering the residual, matching
invariance under the stroboscopic map rather than a plain time-derivative
residual) is a materially larger capability build, out of this task's
bounded scope, and was NOT attempted here. This module crosses the
instability wall for the periodic-orbit case that anchors the torus family;
it does not solve the torus corrector itself.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import cyclerfinder.core.qbcp as qbcp
from cyclerfinder.search.outcome_log import log_outcome

_N_STATE = 6  # canonical PM state: x, y, z, px, py, pz


@dataclass(frozen=True)
class QBCPVariationalOrbitResult:
    """A QBCP periodic-orbit candidate found by seedless spectral collocation.

    ``state0_pm`` is the reconstructed canonical (PM) 6-state at ``t=0``;
    ``period`` is FIXED at the caller's chosen ``period_multiple * T_s`` (not
    solved for -- see module docstring). ``residual_rms`` is the RMS of the
    harmonic-balance residual (the six first-order QBCP EOM right-hand
    sides) at the collocation points. ``closure_residual`` is an INDEPENDENT
    check: ``state0_pm`` propagated for ``period`` through the TRUE
    (non-truncated) nonlinear ``qbcp_eom`` via ``solve_ivp``, compared
    against ``state0_pm`` itself -- not circular with ``residual_rms``.
    """

    state0_pm: NDArray[np.float64]
    period: float
    converged: bool
    residual_rms: float
    closure_residual: float
    n_harmonics: int
    n_collocation: int
    n_restarts_tried: int
    cost: float
    raw_coeffs: NDArray[np.float64]


def _n_free(n_harmonics: int) -> int:
    return _N_STATE * (1 + 2 * n_harmonics)


def _unpack(
    z: NDArray[np.float64], n_harmonics: int
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Unpack the free-variable vector into (dc[6], cos[6,n], sin[6,n]).

    Row order is (x, y, z, px, py, pz), matching ``qbcp.qbcp_eom``'s state
    ordering. No phase gauge and no amplitude anchor: the period is fixed
    (see module docstring), so there is no continuous time-shift degeneracy
    to remove, and empirically there is no trivial zero-amplitude attractor
    to guard against (the Sun's forcing never vanishes).
    """
    n = n_harmonics
    dc = np.zeros(_N_STATE, dtype=np.float64)
    cosc = np.zeros((_N_STATE, n), dtype=np.float64)
    sinc = np.zeros((_N_STATE, n), dtype=np.float64)
    idx = 0
    for v in range(_N_STATE):
        dc[v] = z[idx]
        idx += 1
        if n > 0:
            cosc[v] = z[idx : idx + n]
            idx += n
            sinc[v] = z[idx : idx + n]
            idx += n
    return dc, cosc, sinc


def _eval_series(
    theta: NDArray[np.float64], dc: float, ccos: NDArray[np.float64], csin: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Evaluate a truncated Fourier series and its 1st theta-derivative.

    Only a first derivative is needed (unlike the CR3BP module's second
    derivative): the QBCP canonical system is first-order in all six state
    components.
    """
    n = ccos.shape[0]
    if n == 0:
        return np.full_like(theta, dc), np.zeros_like(theta)
    k = np.arange(1, n + 1)
    ang = np.outer(theta, k)
    cos_ang, sin_ang = np.cos(ang), np.sin(ang)
    f = dc + cos_ang @ ccos + sin_ang @ csin
    fp = cos_ang @ (k * csin) - sin_ang @ (k * ccos)
    return f, fp


def _harmonic_balance_residual(
    z: NDArray[np.float64],
    t_grid: NDArray[np.float64],
    omega_loop: float,
    n_harmonics: int,
    system: qbcp.QBCPSystem,
) -> NDArray[np.float64]:
    """QBCP EOM residual of the Fourier-parameterized 6-state loop at collocation points."""
    dc, cosc, sinc = _unpack(z, n_harmonics)
    theta = omega_loop * t_grid
    n_pts = len(t_grid)
    states = np.empty((n_pts, _N_STATE), dtype=np.float64)
    derivs = np.empty((n_pts, _N_STATE), dtype=np.float64)
    for v in range(_N_STATE):
        f, fp = _eval_series(theta, dc[v], cosc[v], sinc[v])
        states[:, v] = f
        derivs[:, v] = omega_loop * fp
    res = np.empty((n_pts, _N_STATE), dtype=np.float64)
    for i in range(n_pts):
        rhs = qbcp.qbcp_eom(float(t_grid[i]), states[i], system)
        res[i] = derivs[i] - rhs
    return res.flatten()


def _reconstruct_state0(z: NDArray[np.float64], n_harmonics: int) -> NDArray[np.float64]:
    dc, cosc, _sinc = _unpack(z, n_harmonics)
    state0 = dc + cosc.sum(axis=1)  # cos(0) = 1 for every harmonic
    return np.asarray(state0, dtype=np.float64)


def discover_qbcp_periodic_orbit(
    system: qbcp.QBCPSystem,
    *,
    n_harmonics: int = 32,
    n_collocation: int | None = None,
    period_multiple: int = 1,
    center_guess: tuple[float, float, float, float, float, float] = (0.85, 0.0, 0.0, 0.0, 0.8, 0.0),
    coefficient_noise: float = 0.01,
    n_restarts: int = 3,
    rng: np.random.Generator | None = None,
    tol: float = 1e-6,
    max_nfev: int = 1500,
    warm_start: NDArray[np.float64] | None = None,
) -> QBCPVariationalOrbitResult:
    """Seedless QBCP periodic-orbit discovery via 6-state Fourier harmonic balance.

    Represents a candidate closed loop as a degree-``n_harmonics`` truncated
    real Fourier series in EACH of the six canonical (PM) state components
    (x, y, z, px, py, pz) over a FIXED period ``period_multiple * T_s`` (``T_s``
    = ``system.sun_period_tu``, the Sun's synodic period -- the only period at
    which a genuine QBCP periodic orbit can exist, since the model is
    explicitly time-periodic with this period, not autonomous). Minimizes the
    L2 residual of the six first-order QBCP equations of motion
    (``qbcp.qbcp_eom``) at ``n_collocation`` points spread uniformly over the
    period, via ``scipy.optimize.least_squares`` (Levenberg-Marquardt).

    NO integration, no continuation bootstrap (no CR3BP-L1-to-BCR4BP-to-QBCP
    chain), and no near-orbit initial guess are required: the caller supplies
    only a rough center state near the target collinear point. This is the
    QBCP analogue of `#606`'s CR3BP seedless discovery, and directly attacks
    the "violently unstable" EM-L1/L2 region `#544` root-caused as the reason
    every shooting-based QBCP corrector needs an expensive multi-stage
    bootstrap (or fails outright).

    Multi-start: ``n_restarts`` independent random cold starts (random
    higher-harmonic coefficients, ``center_guess`` jittered) are tried; the
    FIRST one achieving ``residual_rms < tol`` is returned. If none converge,
    the best (lowest-cost) attempt is returned with ``converged=False``.

    ``warm_start``, if given, is tried FIRST as the exact initial coefficient
    vector (e.g. a previous call's ``raw_coeffs``, for a caller-driven
    continuation using only this module's own solves). Must match
    ``n_harmonics``.

    Why the default ``n_harmonics=32`` is much higher than the CR3BP sibling
    module's default of 8 (empirically measured, not a guess): the
    ``residual_rms < tol`` criterion alone is NOT a reliable convergence
    indicator here at low harmonic counts, because it only checks the
    collocation-point residual, not what happens BETWEEN points once the
    truncated series is propagated through the true nonlinear (and violently
    unstable) flow. At ``n_harmonics=8``, a fit can satisfy
    ``residual_rms < 1e-6`` yet still carry a ``closure_residual`` of
    ``~0.6`` (an O(1) fraction of the state itself -- not remotely periodic);
    at 16, ``closure_residual`` is still ``~0.26``; only from ``n_harmonics
    >= 32`` does ``closure_residual`` collapse to ``1e-6``-``1e-7`` (see
    ``tests/search/test_variational_periodic_orbit_qbcp.py``). This is the
    harmonic-balance signature of the same ~1e6-1e8 per-period amplification
    #544 root-caused: with too few Fourier degrees of freedom, the truncated
    series can zero the residual AT the sampled points while still being
    wrong enough BETWEEN them that the violent instability blows the gap up
    to O(1) over one period. **Callers must not treat ``converged=True``
    alone as proof of a genuine periodic orbit for this system -- always
    check ``closure_residual`` too.**

    Why ``method="trf"`` here, unlike the CR3BP sibling's ``method="lm"``
    (found empirically, not a style choice): at ``n_harmonics=32`` there are
    390 free variables, and ``scipy.optimize.least_squares``'s ``max_nfev``
    is NOT a reliable wall-clock bound for ``method="lm"`` (MINPACK's
    ``lmdif``) on this problem -- measured directly: requesting
    ``max_nfev=100`` cost 38,710 actual residual calls (~390x over),
    because MINPACK's own internal step-search/Jacobian-refresh loop is not
    gated by the reported ``nfev`` the way the docs imply, and this QBCP
    problem's ill-conditioned, violently-unstable Jacobian triggers that
    loop heavily for an unlucky random restart. ``method="trf"`` (Trust
    Region Reflective, pure-Python) still overshoots the requested
    ``max_nfev`` somewhat (measured ~4x-23x over, not exact), but stays
    bounded to tens of seconds per attempt in practice rather than the
    multi-minute-plus runs ``"lm"`` produced for the same unlucky seeds --
    this is why ``max_nfev`` defaults to a modest 1500 rather than a large
    number: most random restarts converge in a few seconds, and the rare
    unlucky one is capped rather than left open-ended. This is a genuine,
    problem-specific difference from the CR3BP module, not a copy-paste
    inconsistency to "fix" by matching it.
    """
    if n_harmonics < 1:
        raise ValueError(
            f"discover_qbcp_periodic_orbit: n_harmonics must be >= 1, got {n_harmonics}"
        )
    if period_multiple < 1:
        raise ValueError(
            f"discover_qbcp_periodic_orbit: period_multiple must be >= 1, got {period_multiple}"
        )
    period = period_multiple * system.sun_period_tu
    omega_loop = 2.0 * np.pi / period
    n_coll = n_collocation if n_collocation is not None else 6 * n_harmonics
    t_grid = np.linspace(0.0, period, n_coll, endpoint=False)
    n_free = _n_free(n_harmonics)
    gen = rng if rng is not None else np.random.default_rng()

    if warm_start is not None and warm_start.shape != (n_free,):
        raise ValueError(
            f"discover_qbcp_periodic_orbit: warm_start shape {warm_start.shape} does not "
            f"match the expected free-variable count {(n_free,)} for n_harmonics={n_harmonics}."
        )

    best_z: NDArray[np.float64] | None = None
    best_cost = float("inf")
    n_random_restarts = max(1, n_restarts)
    n_attempts = n_random_restarts + (1 if warm_start is not None else 0)
    attempt = 0
    for attempt in range(n_attempts):
        if warm_start is not None and attempt == 0:
            z0 = warm_start.copy()
        else:
            z0 = np.zeros(n_free)
            for v in range(_N_STATE):
                idx = v * (1 + 2 * n_harmonics)
                spread = 0.03 if v == 0 else 0.02
                z0[idx] = center_guess[v] + gen.normal(scale=spread)
                if n_harmonics > 0:
                    z0[idx + 1 : idx + 1 + 2 * n_harmonics] = gen.normal(
                        scale=coefficient_noise, size=2 * n_harmonics
                    )
        sol = least_squares(
            _harmonic_balance_residual,
            z0,
            args=(t_grid, omega_loop, n_harmonics, system),
            method="trf",
            xtol=1e-15,
            ftol=1e-15,
            gtol=1e-15,
            max_nfev=max_nfev,
        )
        if sol.cost < best_cost:
            best_cost = sol.cost
            best_z = sol.x.copy()
        rms = float(np.sqrt(2.0 * sol.cost / (_N_STATE * n_coll)))
        if rms < tol:
            best_z = sol.x.copy()
            best_cost = sol.cost
            break

    assert best_z is not None
    state0_pm = _reconstruct_state0(best_z, n_harmonics)
    residual_rms = float(np.sqrt(2.0 * best_cost / (_N_STATE * n_coll)))
    sol_prop = solve_ivp(
        lambda t, y: qbcp.qbcp_eom(t, y, system),
        (0.0, period),
        state0_pm,
        method="DOP853",
        rtol=1e-13,
        atol=1e-13,
    )
    closure_residual = float(np.linalg.norm(sol_prop.y[:, -1] - state0_pm))
    converged = residual_rms < tol
    log_outcome(
        solver="variational_periodic_orbit_qbcp.discover_qbcp_periodic_orbit",
        inputs={
            "mu": float(system.mu),
            "n_harmonics": int(n_harmonics),
            "n_collocation": int(n_coll),
            "period_multiple": int(period_multiple),
            "center_guess": list(center_guess),
        },
        outcome={
            "converged": bool(converged),
            "residual_rms": residual_rms,
            "closure_residual": closure_residual,
            "period": float(period),
        },
        meta={"model": "qbcp"},
    )
    return QBCPVariationalOrbitResult(
        state0_pm=state0_pm,
        period=period,
        converged=converged,
        residual_rms=residual_rms,
        closure_residual=closure_residual,
        n_harmonics=n_harmonics,
        n_collocation=n_coll,
        n_restarts_tried=attempt + 1,
        cost=best_cost,
        raw_coeffs=best_z.copy(),
    )
