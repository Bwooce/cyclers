"""Seedless spectral (harmonic-balance) periodic-orbit discovery for the CR3BP (#606).

Motivation
----------
Every periodic-orbit corrector already in this codebase (``cr3bp_periodic.py``,
``pseudo_arclength.py``, ``cr3bp_continuation.py``, ``deflated_newton.py``) is a
Newton/shooting method: it starts from a state already close to a genuine
periodic orbit and drives a closure residual to zero. That only ever converges
to whichever family's basin of attraction contains the seed -- documented
repeatedly in this project's own history as a structural wall (#388/S1L1's
continuation relaxing to an off-anchor basin; #538/#544's QBCP torus corrector
that cannot converge the violently-unstable EM-L1/L2 region; #556's
``richardson_halo_seed``/``correct_general_periodic_3d`` combination, which
this module's own pilot re-confirms fails for essentially the entire EM L1
halo family between its bifurcation (C=3.1745) and C~3.146 -- see
``scripts/run_606_variational_pilot.py``).

This module builds a genuinely different, SEEDLESS route: represent a
candidate closed loop as a truncated real Fourier series in time over one
candidate period, and require the loop to satisfy the CR3BP equations of
motion (not a shooting/closure condition) at a set of collocation points
spread over the period. No integration is involved and no initial guess needs
to already be near a genuine trajectory arc -- only a rough, generic starting
Fourier series (small-amplitude loop, offset location, mismatched period) is
required, exactly the "cold start" a basin-restricted shooting method cannot
tolerate.

Mathematical formulation and precedent
---------------------------------------
The classical precedent is the figure-eight three-body orbit: Moore, C.
(1993), "Braids in classical dynamics", Phys. Rev. Lett. 70, 3675-3679,
discretized the three-body action functional via a truncated Fourier series
and minimized it by gradient descent from a symmetry-constrained topological
guess, then refined the result with a shooting method. Chenciner, A. &
Montgomery, R. (2000), "A remarkable periodic solution of the three-body
problem in the case of equal masses", Annals of Mathematics 152, 881-901,
proved existence via direct minimization of the action functional over a
symmetry-constrained loop space, building on Moore's numerics.

This module implements the mathematically-equivalent HARMONIC-BALANCE /
spectral-collocation variant rather than literal action-gradient-descent:
minimize (via ``scipy.optimize.least_squares``, an off-the-shelf nonlinear
least-squares solver -- no new dependency) the L2 residual of the CR3BP
EQUATIONS OF MOTION evaluated at collocation points along the Fourier-
parameterized loop, rather than the action integral itself. The two are
stationary at the same point: the CR3BP rotating-frame Lagrangian is
``L = 1/2(v + omega x r)^2 + Ubar(r)`` (kinetic term including the Coriolis
cross term plus the effective potential; standard form, e.g. Szebehely 1967),
and its Euler-Lagrange equations ARE the CR3BP equations of motion
(``cyclerfinder.core.cr3bp.cr3bp_eom``) -- so a trajectory makes the residual
of the EOM vanish everywhere IFF it is a stationary point of the action.
Driving the EOM residual to zero at a discrete set of points (harmonic
balance / weighted-residual collocation, a textbook technique for
approximating periodic solutions of nonlinear ODEs by truncated Fourier
series) is far more tractable with ``scipy.optimize.least_squares`` --
robust trust-region/Levenberg-Marquardt solvers with off-the-shelf Jacobians
-- than a hand-rolled action-gradient-descent loop, and avoids adding a new
autodiff dependency (this project has no ``jax``/``torch`` dependency
currently and the CR3BP EOM is smooth and cheap to evaluate at collocation
points, so finite-difference Jacobians via ``least_squares`` are adequate).
Moore's own two-stage pipeline (spectral minimization for the shape, then a
shooting method for high-precision refinement) is mirrored here: this
module's raw solution is already very tightly converged (see
``tests/search/test_variational_periodic_orbit.py``), but callers wanting
machine-precision closure can still polish it with
:func:`cyclerfinder.search.cr3bp_periodic.correct_periodic`.

Gauge fixing (breaking trivial/continuous symmetries)
------------------------------------------------------
A truncated Fourier loop has two symmetries the raw EOM residual does not
break: (1) an overall time-shift (any phase of the loop is an equally valid
parameterization) and (2) the L-point EQUILIBRIUM is itself an exact
zero-residual "periodic orbit" of zero amplitude, a trivial attractor the
optimizer can collapse onto. Both are fixed by construction, not learned:

* Phase gauge: the sine coefficient of x's first harmonic is fixed at 0. This
  pins the time origin to an extremum of x(t) (empirically, for every genuine
  symmetric CR3BP periodic orbit tested here, this coincides with the
  perpendicular-crossing convention used elsewhere in this codebase --
  ``y(0) = 0``, ``xdot(0) = 0`` -- to within the solver's own residual
  tolerance; not assumed, just observed).
* Amplitude anchor(s): the caller fixes one or more leading Fourier
  coefficients (``anchor_x1``/``anchor_y1``/``anchor_z1``) to a nonzero
  value, analogous to how
  :func:`cyclerfinder.search.cr3bp_seed_generator.lyapunov_seed` fixes an
  amplitude ``Ax`` or
  :func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`
  fixes a Jacobi constant -- a legitimate free-variable-elimination choice
  ("how big an orbit do you want"), not information leaked from a target
  orbit's own converged state. At least one anchor must be given (checked)
  so the trivial equilibrium branch cannot be reached.

For a genuinely 3D (halo-topology) orbit, anchoring z ALONE is not enough:
this module's own pilot found that leaving the in-plane (y) amplitude
completely free lets the optimizer slide onto the ALREADY-KNOWN, easier,
nearly-decoupled small-amplitude vertical-Lyapunov/"tulip" branch (Koon, Lo,
Marsden & Ross 2011 Ch. 2.5-2.7's linearized out-of-plane mode -- a genuine,
distinct family, not the sought halo) instead of the coupled halo family.
Reproducing a genuine halo requires anchoring BOTH ``anchor_y1`` and
``anchor_z1`` together, in a ratio Fourier-fit from a nearby already-built
family member -- but a CONSTANT ratio is only approximately right (the true
ratio drifts with amplitude, per ``scripts/run_606_variational_pilot.py``'s
own measurements), so a fixed-ratio dual anchor plateaus at a modest
(non-machine-precision) residual once the anchor pair drifts from where the
ratio was calibrated. The robust recipe (validated in the pilot script) is
instead: anchor Z ALONE, but supply ``warm_start`` -- a full coefficient
vector Fourier-fit from ONE already-known family member (which already
carries the correct nonlinear y-z coupling) -- and then walk the amplitude
anchor in small steps, re-using each step's ``raw_coeffs`` as the next
``warm_start``. This is a genuine natural-parameter continuation done
ENTIRELY with this module's own solves (never the existing correctors after
the single bootstrap point), and in the pilot it crosses the ENTIRE #556
near-bifurcation wall to near-machine precision. This is itself an honest
finding: a seedless method still has its OWN family-selection bias (which
branch is "cheapest" to reach from a generic cold start, and how well a
fixed gauge tracks a curving family); it just has a DIFFERENT bias than
shooting/continuation, not none at all.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import least_squares

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.search.outcome_log import log_outcome

_COORDS = ("x", "y", "z")


@dataclass(frozen=True)
class VariationalOrbitResult:
    """A CR3BP periodic-orbit candidate found by seedless spectral collocation.

    ``state0``/``period`` is the reconstructed 6-state (at the Fourier
    series' ``theta=0`` phase) and full period; ``jacobi`` is evaluated at
    ``state0``. ``residual_rms`` is the RMS of the harmonic-balance EOM
    residual at the collocation points (the quantity actually minimized);
    ``closure_residual`` is an INDEPENDENT check -- ``state0`` numerically
    propagated for ``period`` via the true (non-truncated) nonlinear EOM,
    compared against ``state0`` itself -- so a tiny ``closure_residual`` is
    not circular with a tiny ``residual_rms`` (one is a property of the
    truncated series, the other of the actual flow).
    """

    state0: NDArray[np.float64]
    period: float
    jacobi: float
    converged: bool
    residual_rms: float
    closure_residual: float
    n_harmonics: int
    n_collocation: int
    n_restarts_tried: int
    cost: float
    raw_coeffs: NDArray[np.float64]
    """The full free-variable vector at the returned solution. Reusable as
    ``warm_start`` in a subsequent :func:`discover_periodic_orbit` call with
    the SAME ``n_harmonics`` and the same anchors set to non-``None`` (their
    VALUES may differ) -- e.g. to continue a family through a caller-driven
    amplitude schedule using this method's own solves throughout, never the
    existing shooting correctors. See ``scripts/run_606_variational_pilot.py``.
    """


def _unpack(
    z: NDArray[np.float64],
    n_harmonics: int,
    anchor_x1: float | None,
    anchor_y1: float | None,
    anchor_z1: float | None,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], float]:
    """Unpack the free-variable vector into (dc[3], cos[3,N], sin[3,N], period).

    ``cos``/``sin`` rows are indexed 0=x, 1=y, 2=z; columns are harmonics
    k=1..N (0-based column k-1). ``sin[0,0]`` (x's first sine coefficient) is
    ALWAYS fixed at 0 (phase gauge); any of ``anchor_x1``/``anchor_y1``/
    ``anchor_z1`` supplied fixes ``cos[0,0]``/``sin[1,0]``/``cos[2,0]``
    respectively (the amplitude-anchor convention: x, z use the cosine
    fundamental, y uses the sine fundamental, matching a symmetric orbit's
    natural even/odd split about the ``theta=0`` phase).
    """
    n = n_harmonics
    i = 0
    dc = np.array([z[0], z[1], z[2]], dtype=np.float64)
    i = 3
    cosc = np.zeros((3, n))
    sinc = np.zeros((3, n))
    if anchor_x1 is None:
        cosc[0, 0] = z[i]
        i += 1
    else:
        cosc[0, 0] = anchor_x1
    # sinc[0, 0] stays 0 (phase gauge -- never a free variable, never anchored).
    if n > 1:
        cosc[0, 1:] = z[i : i + n - 1]
        i += n - 1
        sinc[0, 1:] = z[i : i + n - 1]
        i += n - 1
    if anchor_y1 is None:
        sinc[1, 0] = z[i]
        i += 1
    else:
        sinc[1, 0] = anchor_y1
    cosc[1, 0] = z[i]
    i += 1
    if n > 1:
        cosc[1, 1:] = z[i : i + n - 1]
        i += n - 1
        sinc[1, 1:] = z[i : i + n - 1]
        i += n - 1
    if anchor_z1 is None:
        cosc[2, 0] = z[i]
        i += 1
    else:
        cosc[2, 0] = anchor_z1
    sinc[2, 0] = z[i]
    i += 1
    if n > 1:
        cosc[2, 1:] = z[i : i + n - 1]
        i += n - 1
        sinc[2, 1:] = z[i : i + n - 1]
        i += n - 1
    period = float(np.exp(z[i]))
    return dc, cosc, sinc, period


def _n_free_actual(
    n_harmonics: int, anchor_x1: float | None, anchor_y1: float | None, anchor_z1: float | None
) -> int:
    n = n_harmonics
    total = 3  # dc
    total += 0 if anchor_x1 is not None else 1
    total += 2 * (n - 1)  # x-cos[1:], x-sin[1:]
    total += 0 if anchor_y1 is not None else 1
    total += 1 + 2 * (n - 1)  # y-cos[0] always free, plus [1:] both
    total += 0 if anchor_z1 is not None else 1
    total += 1 + 2 * (n - 1)  # z-sin[0] always free, plus [1:] both
    total += 1  # logT
    return total


def _eval_series(
    theta: NDArray[np.float64], dc: float, ccos: NDArray[np.float64], csin: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Evaluate a truncated Fourier series and its 1st/2nd theta-derivatives."""
    n = ccos.shape[0]
    k = np.arange(1, n + 1)
    ang = np.outer(theta, k)
    cos_ang, sin_ang = np.cos(ang), np.sin(ang)
    f = dc + cos_ang @ ccos + sin_ang @ csin
    fp = cos_ang @ (k * csin) - sin_ang @ (k * ccos)
    fpp = -(cos_ang @ (k * k * ccos) + sin_ang @ (k * k * csin))
    return f, fp, fpp


def _harmonic_balance_residual(
    z: NDArray[np.float64],
    theta: NDArray[np.float64],
    mu: float,
    n_harmonics: int,
    anchor_x1: float | None,
    anchor_y1: float | None,
    anchor_z1: float | None,
) -> NDArray[np.float64]:
    """CR3BP EOM residual of the Fourier-parameterized loop at collocation points."""
    dc, cosc, sinc, period = _unpack(z, n_harmonics, anchor_x1, anchor_y1, anchor_z1)
    x, xp, xpp = _eval_series(theta, dc[0], cosc[0], sinc[0])
    y, yp, ypp = _eval_series(theta, dc[1], cosc[1], sinc[1])
    zc, zp, zpp = _eval_series(theta, dc[2], cosc[2], sinc[2])
    w = 2.0 * np.pi / period
    xdot, ydot, _zdot = w * xp, w * yp, w * zp
    xddot, yddot, zddot = w * w * xpp, w * w * ypp, w * w * zpp
    r1 = np.sqrt((x + mu) ** 2 + y * y + zc * zc)
    r2 = np.sqrt((x - 1.0 + mu) ** 2 + y * y + zc * zc)
    r1c, r2c = r1**3, r2**3
    ubar_x = x - (1.0 - mu) * (x + mu) / r1c - mu * (x - 1.0 + mu) / r2c
    ubar_y = y - (1.0 - mu) * y / r1c - mu * y / r2c
    ubar_z = -(1.0 - mu) * zc / r1c - mu * zc / r2c
    fx = xddot - 2.0 * ydot - ubar_x
    fy = yddot + 2.0 * xdot - ubar_y
    fz = zddot - ubar_z
    return np.concatenate([fx, fy, fz])


def _reconstruct_state0(
    z: NDArray[np.float64],
    n_harmonics: int,
    anchor_x1: float | None,
    anchor_y1: float | None,
    anchor_z1: float | None,
) -> tuple[NDArray[np.float64], float]:
    dc, cosc, sinc, period = _unpack(z, n_harmonics, anchor_x1, anchor_y1, anchor_z1)
    k = np.arange(1, n_harmonics + 1)
    pos0 = dc + cosc.sum(axis=1)  # cos(0)=1, sin(0)=0
    w = 2.0 * np.pi / period
    vel0 = w * (k * sinc).sum(axis=1)
    state0 = np.array([pos0[0], pos0[1], pos0[2], vel0[0], vel0[1], vel0[2]], dtype=np.float64)
    return state0, period


def _residual_fixed_jacobi(
    z: NDArray[np.float64],
    theta: NDArray[np.float64],
    mu: float,
    n_harmonics: int,
    jacobi_target: float,
    jacobi_weight: float,
) -> NDArray[np.float64]:
    """Harmonic-balance EOM residual PLUS one extra row pinning the Jacobi
    constant, with ALL THREE amplitude anchors free (see
    :func:`discover_periodic_orbit_fixed_jacobi`).
    """
    base = _harmonic_balance_residual(z, theta, mu, n_harmonics, None, None, None)
    state0, _period = _reconstruct_state0(z, n_harmonics, None, None, None)
    c = cr3bp.jacobi_constant(state0, mu)
    jac_row = np.array([jacobi_weight * (c - jacobi_target)])
    return np.concatenate([base, jac_row])


def _fixed_jacobi_continuation_solve(
    z_start: NDArray[np.float64],
    theta: NDArray[np.float64],
    mu: float,
    n_harmonics: int,
    jacobi_target: float,
    jacobi_weight: float,
    n_continuation_steps: int,
    max_nfev: int,
) -> tuple[NDArray[np.float64], float]:
    """Natural-parameter continuation IN THE TARGET JACOBI CONSTANT, from
    whatever Jacobi constant the cold-start ``z_start`` happens to reconstruct
    to, in ``n_continuation_steps`` steps up to ``jacobi_target``.

    A single ``least_squares`` solve jumping straight from a random cold
    start to the final ``jacobi_target`` empirically stalls at a
    compromise point (small-but-nonzero EOM residual AND small-but-nonzero
    Jacobi error simultaneously, neither actually zero) far more often than
    a solve that only has to move a SMALL step in target Jacobi each time,
    warm-started from the previous step's converged coefficients -- the
    same warm-started-natural-parameter-continuation idea this module's own
    ``warm_start`` docstring already uses for amplitude, applied here to the
    Jacobi constant instead (mirrors this project's existing corridor/
    continuation conventions, e.g. #629's rho-corridor). Returns the final
    ``(raw_coeffs, cost)``.
    """
    state0_0, _period0 = _reconstruct_state0(z_start, n_harmonics, None, None, None)
    c0 = cr3bp.jacobi_constant(state0_0, mu)
    schedule = np.linspace(c0, jacobi_target, max(2, n_continuation_steps))[1:]
    z = z_start.copy()
    sol = None
    for c_step in schedule:
        sol = least_squares(
            _residual_fixed_jacobi,
            z,
            args=(theta, mu, n_harmonics, float(c_step), jacobi_weight),
            method="lm",
            xtol=1e-15,
            ftol=1e-15,
            gtol=1e-15,
            max_nfev=max_nfev,
        )
        z = sol.x
    assert sol is not None
    return z, float(sol.cost)


def discover_periodic_orbit_fixed_jacobi(
    system: cr3bp.CR3BPSystem,
    jacobi_target: float,
    *,
    n_harmonics: int = 16,
    n_collocation: int | None = None,
    z0: NDArray[np.float64] | None = None,
    center_guess: tuple[float, float, float] = (0.8, 0.0, 0.0),
    period_guess: float = 3.0,
    coefficient_noise: float = 0.05,
    jacobi_weight: float = 50.0,
    n_continuation_steps: int = 8,
    rng: np.random.Generator | None = None,
    tol: float = 1e-8,
    jacobi_tol: float = 1e-6,
    max_nfev: int = 15000,
) -> VariationalOrbitResult:
    """Single-attempt fixed-Jacobi variant of :func:`discover_periodic_orbit` (#648).

    Unlike :func:`discover_periodic_orbit`, no amplitude anchor is fixed by
    the caller -- ALL THREE (``anchor_x1``/``anchor_y1``/``anchor_z1``) are
    free variables. In their place, one extra residual row pins the Jacobi
    constant of the reconstructed ``state0`` to ``jacobi_target`` (weighted
    by ``jacobi_weight``), reached via ``n_continuation_steps`` of natural-
    parameter continuation in the target Jacobi constant itself (see
    :func:`_fixed_jacobi_continuation_solve` -- empirically much more
    reliable than a single solve jumping straight to ``jacobi_target`` from a
    cold start). This is the natural fixed-energy analogue of
    ``cr3bp_periodic.correct_symmetric_fixed_jacobi`` (which fixes ``C`` by
    releasing ``x0`` instead of an amplitude anchor), needed so a caller can
    ask "what periodic orbits exist at THIS Jacobi constant" without already
    knowing which family's amplitude to anchor.

    The trivial zero-amplitude libration-point equilibrium remains excluded:
    its OWN Jacobi constant is pinned at ``C_L(mu)`` for whichever L-point it
    sits at, so as long as ``jacobi_target`` is not exactly one of those five
    values (never true for a generic target), the Jacobi-fixing row is
    already unsatisfied at zero amplitude and the trivial branch is not a
    root of the augmented system.

    Single attempt only (no internal multistart, unlike
    :func:`discover_periodic_orbit`) -- callers wanting multistart/deflated
    enumeration should drive repeated calls themselves (see
    :mod:`cyclerfinder.search.deflated_variational_periodic_orbit`).
    ``converged`` requires BOTH the harmonic-balance EOM residual RMS below
    ``tol`` AND the Jacobi-constant error below ``jacobi_tol`` -- a candidate
    passing this is STILL only a Fourier-residual claim, not yet a certified
    trajectory (see the Radau cross-check discipline in
    :mod:`cyclerfinder.search.deflated_variational_periodic_orbit`).
    """
    if n_harmonics < 1:
        raise ValueError(
            f"discover_periodic_orbit_fixed_jacobi: n_harmonics must be >= 1, got {n_harmonics}"
        )
    mu = system.mu
    n_coll = n_collocation if n_collocation is not None else 6 * n_harmonics
    theta = np.linspace(0.0, 2.0 * np.pi, n_coll, endpoint=False)
    n_free = _n_free_actual(n_harmonics, None, None, None)
    gen = rng if rng is not None else np.random.default_rng()

    if z0 is not None:
        if z0.shape != (n_free,):
            raise ValueError(
                f"discover_periodic_orbit_fixed_jacobi: z0 shape {z0.shape} does not match "
                f"the expected free-variable count {(n_free,)} for n_harmonics={n_harmonics}."
            )
        z_start = z0.copy()
    else:
        z_start = np.zeros(n_free)
        z_start[0] = center_guess[0] + gen.normal(scale=0.05)
        z_start[1] = center_guess[1] + gen.normal(scale=0.02)
        z_start[2] = center_guess[2] + gen.normal(scale=0.02)
        z_start[3:-1] = gen.normal(scale=coefficient_noise, size=n_free - 4)
        z_start[-1] = np.log(period_guess * float(gen.uniform(0.5, 1.8)))

    z_final, cost = _fixed_jacobi_continuation_solve(
        z_start,
        theta,
        mu,
        n_harmonics,
        jacobi_target,
        jacobi_weight,
        n_continuation_steps,
        max_nfev,
    )
    state0, period = _reconstruct_state0(z_final, n_harmonics, None, None, None)
    base_res = _harmonic_balance_residual(z_final, theta, mu, n_harmonics, None, None, None)
    residual_rms = float(np.sqrt(np.mean(base_res**2)))
    jacobi = cr3bp.jacobi_constant(state0, mu)
    jacobi_err = abs(jacobi - jacobi_target)
    arc = cr3bp.propagate(system, state0, period, rtol=1e-13, atol=1e-13)
    closure_residual = float(np.linalg.norm(arc.state_f - state0))
    converged = residual_rms < tol and jacobi_err < jacobi_tol
    log_outcome(
        solver="variational_periodic_orbit.discover_periodic_orbit_fixed_jacobi",
        inputs={
            "mu": float(mu),
            "n_harmonics": int(n_harmonics),
            "n_collocation": int(n_coll),
            "jacobi_target": float(jacobi_target),
            "jacobi_weight": float(jacobi_weight),
            "center_guess": list(center_guess),
            "period_guess": float(period_guess),
        },
        outcome={
            "converged": bool(converged),
            "residual_rms": residual_rms,
            "jacobi_error": float(jacobi_err),
            "closure_residual": closure_residual,
            "period": float(period),
            "jacobi": float(jacobi),
        },
        meta={"primary": system.primary, "secondary": system.secondary},
    )
    return VariationalOrbitResult(
        state0=state0,
        period=period,
        jacobi=jacobi,
        converged=converged,
        residual_rms=residual_rms,
        closure_residual=closure_residual,
        n_harmonics=n_harmonics,
        n_collocation=n_coll,
        n_restarts_tried=1,
        cost=cost,
        raw_coeffs=z_final.copy(),
    )


def discover_periodic_orbit(
    system: cr3bp.CR3BPSystem,
    *,
    n_harmonics: int = 16,
    n_collocation: int | None = None,
    anchor_x1: float | None = None,
    anchor_y1: float | None = None,
    anchor_z1: float | None = None,
    center_guess: tuple[float, float, float] = (0.8, 0.0, 0.0),
    period_guess: float = 3.0,
    coefficient_noise: float = 0.01,
    n_restarts: int = 8,
    rng: np.random.Generator | None = None,
    tol: float = 1e-8,
    max_nfev: int = 30000,
    warm_start: NDArray[np.float64] | None = None,
) -> VariationalOrbitResult:
    """Seedless CR3BP periodic-orbit discovery via Fourier harmonic-balance collocation.

    Represents a candidate closed loop as a degree-``n_harmonics`` truncated
    real Fourier series in each of (x, y, z) over one period, and minimizes
    the L2 residual of the CR3BP equations of motion (see module docstring
    for the action-functional/harmonic-balance equivalence and its
    citation) at ``n_collocation`` points spread uniformly over the period,
    via ``scipy.optimize.least_squares`` (Levenberg-Marquardt).

    NO integration and no near-orbit initial guess are required: the caller
    supplies only a rough center location, a rough period, and at least one
    nonzero amplitude anchor (breaking the trivial time-phase and
    zero-amplitude-equilibrium degeneracies -- see module docstring). This is
    the key structural difference from every shooting/continuation corrector
    already in this codebase, which requires a seed already inside the
    target family's basin of attraction.

    Multi-start: ``n_restarts`` independent random cold starts (random
    higher-harmonic coefficients, ``center_guess``/``period_guess`` jittered)
    are tried; the FIRST one achieving ``residual_rms < tol`` with a period
    inside ``(0.05, 50)`` (a divergence guard, not a physical constraint) is
    returned. If none converge, the best (lowest-cost) attempt is returned
    with ``converged=False``.

    ``warm_start``, if given, is tried FIRST (before any random restart) as
    the exact initial coefficient vector -- e.g. a previous call's
    ``raw_coeffs`` field, at a nearby anchor value, for a caller-driven
    continuation entirely within this seedless method (never touching the
    existing shooting correctors). Must match ``n_harmonics`` and which of
    ``anchor_x1``/``anchor_y1``/``anchor_z1`` are ``None`` vs. given (their
    VALUES may differ from the call that produced it).

    Raises
    ------
    ValueError
        If none of ``anchor_x1``/``anchor_y1``/``anchor_z1`` is given (the
        trivial equilibrium branch would otherwise be reachable), or if
        ``n_harmonics < 1``.
    """
    if n_harmonics < 1:
        raise ValueError(f"discover_periodic_orbit: n_harmonics must be >= 1, got {n_harmonics}")
    if anchor_x1 is None and anchor_y1 is None and anchor_z1 is None:
        raise ValueError(
            "discover_periodic_orbit: at least one of anchor_x1/anchor_y1/anchor_z1 must be "
            "given (a nonzero amplitude anchor) -- otherwise the optimizer can trivially "
            "collapse to the exact zero-amplitude libration-point equilibrium."
        )
    mu = system.mu
    n_coll = n_collocation if n_collocation is not None else 6 * n_harmonics
    theta = np.linspace(0.0, 2.0 * np.pi, n_coll, endpoint=False)
    n_free = _n_free_actual(n_harmonics, anchor_x1, anchor_y1, anchor_z1)
    gen = rng if rng is not None else np.random.default_rng()

    if warm_start is not None and warm_start.shape != (n_free,):
        raise ValueError(
            f"discover_periodic_orbit: warm_start shape {warm_start.shape} does not match "
            f"the expected free-variable count {(n_free,)} for n_harmonics={n_harmonics} and "
            "this call's anchor pattern."
        )

    best_z: NDArray[np.float64] | None = None
    best_cost = float("inf")
    n_random_restarts = max(1, n_restarts)
    n_attempts = n_random_restarts + (1 if warm_start is not None else 0)
    for attempt in range(n_attempts):
        if warm_start is not None and attempt == 0:
            z0 = warm_start.copy()
        else:
            z0 = np.zeros(n_free)
            z0[0] = center_guess[0] + gen.normal(scale=0.05)
            z0[1] = center_guess[1] + gen.normal(scale=0.02)
            z0[2] = center_guess[2] + gen.normal(scale=0.02)
            z0[3:-1] = gen.normal(scale=coefficient_noise, size=n_free - 4)
            z0[-1] = np.log(period_guess * float(gen.uniform(0.5, 1.8)))
        sol = least_squares(
            _harmonic_balance_residual,
            z0,
            args=(theta, mu, n_harmonics, anchor_x1, anchor_y1, anchor_z1),
            method="lm",
            xtol=1e-15,
            ftol=1e-15,
            gtol=1e-15,
            max_nfev=max_nfev,
        )
        _state0_try, period_try = _reconstruct_state0(
            sol.x, n_harmonics, anchor_x1, anchor_y1, anchor_z1
        )
        if not (0.05 < period_try < 50.0):
            continue
        if sol.cost < best_cost:
            best_cost = sol.cost
            best_z = sol.x.copy()
        rms = float(np.sqrt(2.0 * sol.cost / (3 * n_coll)))
        if rms < tol:
            best_z = sol.x.copy()
            best_cost = sol.cost
            break

    assert best_z is not None
    state0, period = _reconstruct_state0(best_z, n_harmonics, anchor_x1, anchor_y1, anchor_z1)
    residual_rms = float(np.sqrt(2.0 * best_cost / (3 * n_coll)))
    jacobi = cr3bp.jacobi_constant(state0, mu)
    arc = cr3bp.propagate(system, state0, period, rtol=1e-13, atol=1e-13)
    closure_residual = float(np.linalg.norm(arc.state_f - state0))
    converged = residual_rms < tol
    log_outcome(
        solver="variational_periodic_orbit.discover_periodic_orbit",
        inputs={
            "mu": float(mu),
            "n_harmonics": int(n_harmonics),
            "n_collocation": int(n_coll),
            "anchor_x1": anchor_x1,
            "anchor_y1": anchor_y1,
            "anchor_z1": anchor_z1,
            "center_guess": list(center_guess),
            "period_guess": float(period_guess),
        },
        outcome={
            "converged": bool(converged),
            "residual_rms": residual_rms,
            "closure_residual": closure_residual,
            "period": float(period),
            "jacobi": float(jacobi),
        },
        meta={"primary": system.primary, "secondary": system.secondary},
    )
    return VariationalOrbitResult(
        state0=state0,
        period=period,
        jacobi=jacobi,
        converged=converged,
        residual_rms=residual_rms,
        closure_residual=closure_residual,
        n_harmonics=n_harmonics,
        n_collocation=n_coll,
        n_restarts_tried=attempt + 1,
        cost=best_cost,
        raw_coeffs=best_z.copy(),
    )
