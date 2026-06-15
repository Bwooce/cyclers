"""Sundman time-regularized CR3BP propagator (#266 Phase 1, tulip-orbit builder).

The standard CR3BP propagator (:mod:`cyclerfinder.core.cr3bp`) loses precision and
runs out of integrator steps as a trajectory grazes a primary (r -> 0): the
acceleration goes as 1/r^2 and the step size collapses faster than DOP853's
error estimator can adapt. Sundman's classical regularisation trades physical
time ``t`` for a *regularised* time ``s`` via

    dt/ds = r1 * r2     (or r1 only, or r2 only)

so that the step size in ``s`` expands as the spacecraft approaches a primary,
keeping the integrator on a benign quasi-analytic curve. The augmented state
``(x, y, z, vx, vy, vz, t)`` is integrated in ``s`` and physical time is
recovered from the last component.

This module COMPOSES with :mod:`cyclerfinder.core.cr3bp` (re-uses ``_r1_r2``
implicitly via the same algebra and ``CR3BPSystem``); it does NOT replace
``cr3bp_eom`` / ``propagate``, which remain the hot path for every test in the
tree. The regularised propagator is opt-in for low-perilune work — tulip orbits
(#266) and any future close-encounter studies.

Reference for the technique: Szebehely, "Theory of Orbits" (1967), §3.7
"Regularization of motion near singularities"; the modern multi-primary form
``dt/ds = r1 * r2`` is the Lemaitre-type two-body regularisation generalised
to the CR3BP (Stiefel & Scheifele, "Linear and Regular Celestial Mechanics",
1971, Ch. III). Pure: math / numpy / scipy + cyclerfinder.core.cr3bp.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

from cyclerfinder.core.cr3bp import CR3BPSystem, cr3bp_eom

Regularization = Literal["r1r2", "r1", "r2"]
_VALID_REGULARIZATIONS: tuple[Regularization, ...] = ("r1r2", "r1", "r2")


def _dt_ds(x: float, y: float, z: float, mu: float, regularization: Regularization) -> float:
    """Return dt/ds for the requested regularisation choice."""
    r1 = math.sqrt((x + mu) ** 2 + y * y + z * z)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y + z * z)
    if regularization == "r1r2":
        return r1 * r2
    if regularization == "r1":
        return r1
    if regularization == "r2":
        return r2
    # Defensive: validated at API entry, but mypy needs the explicit branch.
    raise ValueError(f"sundman: unknown regularization {regularization!r}")


def sundman_rhs(
    s: float,
    state_aug: NDArray[np.float64],
    mu: float,
    regularization: Regularization = "r1r2",
) -> NDArray[np.float64]:
    """Regularised RHS for the augmented state ``(x, y, z, vx, vy, vz, t)``.

    Multiplies every component of the standard CR3BP RHS (in physical time) by
    ``dt/ds`` so that integrating in ``s`` produces the same trajectory but with
    a step size that expands near a primary. The trailing component's RHS is
    ``dt/ds`` itself, so physical time is recovered as the integral of dt/ds
    along the regularised arc.

    ``regularization`` selects the multiplier:
      - ``"r1r2"`` (default): both-primary regularisation, dt/ds = r1*r2;
      - ``"r1"``: Earth-only, dt/ds = r1;
      - ``"r2"``: Moon-only, dt/ds = r2.

    For tulip orbits (low perilune) ``"r1r2"`` or ``"r2"`` are appropriate;
    ``"r1r2"`` is the safe default because it also handles any Earth-side
    close approaches that arise in continuation.
    """
    x, y, z = float(state_aug[0]), float(state_aug[1]), float(state_aug[2])
    dtds = _dt_ds(x, y, z, mu, regularization)
    f6 = cr3bp_eom(0.0, state_aug[:6], mu)
    out = np.empty(7, dtype=np.float64)
    out[:6] = f6 * dtds
    out[6] = dtds
    return out


@dataclass(frozen=True)
class RegularizedArc:
    """Result of a Sundman-regularized propagation.

    Attributes
    ----------
    s : NDArray[np.float64]
        The regularised-time grid returned by ``solve_ivp`` (shape ``(N,)``).
    state_at_s : NDArray[np.float64]
        The 6-component CR3BP state at each ``s`` (shape ``(6, N)``). The
        augmented time component is dropped here and exposed via ``t_at_s``.
    t_at_s : NDArray[np.float64]
        The physical time at each ``s`` (shape ``(N,)``).
    solver_success : bool
        Mirrors ``solve_ivp(...).success``.
    nfev : int
        RHS evaluations the integrator performed (the headline integrator-cost
        metric — used to demonstrate the regularised win at low perilune).
    """

    s: NDArray[np.float64]
    state_at_s: NDArray[np.float64]
    t_at_s: NDArray[np.float64]
    solver_success: bool
    nfev: int


def _validate_regularization(regularization: str) -> Regularization:
    if regularization not in _VALID_REGULARIZATIONS:
        raise ValueError(
            f"propagate_regularized: unknown regularization {regularization!r}; "
            f"expected one of {_VALID_REGULARIZATIONS}"
        )
    # Narrow the type after the membership check (Literal-cast).
    return regularization  # mypy treats this as Regularization via narrowing


def propagate_regularized(
    system: CR3BPSystem,
    state0: NDArray[np.float64],
    s_span: tuple[float, float],
    *,
    mu: float | None = None,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    max_step: float = np.inf,
    dense_output: bool = False,
    regularization: str = "r1r2",
    t_stop: float | None = None,
) -> RegularizedArc:
    """Propagate a CR3BP state in regularised time ``s`` over ``s_span``.

    Augments ``state0`` (6-vector) with ``t = 0`` and integrates the
    :func:`sundman_rhs` augmented system with DOP853 (matches the convention in
    :func:`cyclerfinder.core.cr3bp.propagate`).

    Parameters
    ----------
    system :
        CR3BP system carrying the mass ratio. If ``mu`` is given it overrides
        ``system.mu`` (useful for explicit unit testing against a known mu).
    state0 :
        6-vector ``(x, y, z, vx, vy, vz)`` in CR3BP nondimensional units.
    s_span :
        Half-open span ``(s0, sf)`` in regularised time. The caller is
        responsible for choosing ``sf`` large enough to cover the physical time
        range of interest — use :func:`physical_to_regularized_span` to map a
        physical-time horizon to ``s``.
    rtol, atol, max_step :
        Forwarded to ``scipy.integrate.solve_ivp``.
    dense_output :
        Forwarded; when ``True`` the underlying solver builds a dense
        interpolant (not exposed here — the public surface returns the discrete
        ``s`` / ``state_at_s`` / ``t_at_s`` arrays; dense_output is honoured
        only to control solver behaviour, not to leak the interpolant).
    regularization :
        ``"r1r2"`` (default), ``"r1"``, or ``"r2"``. See :func:`sundman_rhs`.
    t_stop :
        Optional terminal physical time. When supplied, a terminal solve_ivp
        event ``t(s) = t_stop`` is installed; the integration stops at the
        precise ``s`` where physical time reaches ``t_stop`` (event-resolved to
        the solver's tolerance, not to the discrete-step grid). ``s_span[1]``
        is treated as an upper bound the integrator may stop short of. Useful
        for landing exactly at a target physical time — the standard
        :func:`cyclerfinder.core.cr3bp.propagate` call returns the state at a
        physical time, so this option makes round-tripping the two propagators
        natural.

    Raises
    ------
    ValueError
        On invalid ``regularization`` key, non-finite tolerances, or
        non-positive ``rtol``/``atol``. Mirrors the validation in
        :func:`cyclerfinder.core.cr3bp.propagate`'s scipy underlay.
    RuntimeError
        If the integrator reports failure. The regularised system rarely fails
        — that is the point — but a genuine collision (state passes *through*
        a primary in finite ``s``) will still trip it.
    """
    reg = _validate_regularization(regularization)
    if not (math.isfinite(rtol) and rtol > 0.0):
        raise ValueError(f"propagate_regularized: rtol must be positive finite, got {rtol}")
    if not (math.isfinite(atol) and atol > 0.0):
        raise ValueError(f"propagate_regularized: atol must be positive finite, got {atol}")
    if not (math.isfinite(max_step) or max_step == np.inf) or max_step <= 0.0:
        raise ValueError(
            f"propagate_regularized: max_step must be positive (or +inf), got {max_step}"
        )
    state0_arr = np.asarray(state0, dtype=np.float64)
    if state0_arr.shape != (6,):
        raise ValueError(
            f"propagate_regularized: state0 must be a 6-vector, got shape {state0_arr.shape}"
        )
    mu_used = float(system.mu if mu is None else mu)
    s0, sf = float(s_span[0]), float(s_span[1])
    if sf == s0:
        raise ValueError(f"propagate_regularized: empty s_span ({s0}, {sf})")

    y0 = np.empty(7, dtype=np.float64)
    y0[:6] = state0_arr
    y0[6] = 0.0

    def _rhs(s: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        return sundman_rhs(s, y, mu_used, reg)

    if t_stop is not None and not math.isfinite(t_stop):
        raise ValueError(f"propagate_regularized: t_stop must be finite, got {t_stop}")
    if t_stop is None:
        sol = solve_ivp(
            _rhs,
            (s0, sf),
            y0,
            method="DOP853",
            rtol=rtol,
            atol=atol,
            max_step=max_step,
            dense_output=dense_output,
        )
    else:
        t_stop_val = float(t_stop)

        def _t_event(s: float, y: NDArray[np.float64]) -> float:
            return float(y[6]) - t_stop_val

        _t_event.terminal = True  # type: ignore[attr-defined]
        _t_event.direction = 1.0  # type: ignore[attr-defined]
        sol = solve_ivp(
            _rhs,
            (s0, sf),
            y0,
            method="DOP853",
            rtol=rtol,
            atol=atol,
            max_step=max_step,
            dense_output=dense_output,
            events=_t_event,
        )
    if not sol.success:
        raise RuntimeError(
            f"propagate_regularized: integrator failed at s={sol.t[-1]}: {sol.message}"
        )
    return RegularizedArc(
        s=np.asarray(sol.t, dtype=np.float64),
        state_at_s=np.asarray(sol.y[:6, :], dtype=np.float64),
        t_at_s=np.asarray(sol.y[6, :], dtype=np.float64),
        solver_success=bool(sol.success),
        nfev=int(sol.nfev),
    )


def physical_to_regularized_span(
    system: CR3BPSystem,
    state0: NDArray[np.float64],
    t_span: tuple[float, float],
    *,
    regularization: str = "r1r2",
    initial_ds: float | None = None,
    max_ds: float | None = None,
) -> tuple[float, float]:
    """Estimate an ``s_span`` that brackets a requested physical-time range.

    Strategy: integrate the augmented system forward in ``s`` with a generous
    ``max_step`` and an event that triggers the first time the physical time
    reaches ``t_span[1]``. The ``s`` at the event is returned as the upper
    bound. No optimisation theatre: a single coarse pass.

    ``initial_ds`` / ``max_ds`` are optional scale hints for the solver. When
    not provided, sensible defaults are chosen from the typical magnitude of
    ``dt/ds`` at ``state0`` — a few orbital revolutions' worth of ``s`` at most.

    Returns
    -------
    (s0, sf) :
        ``s0`` mirrors ``t_span[0]`` interpreted as the starting ``s`` (we
        always start the integrator at ``s = 0`` regardless of ``t_span[0]``,
        because regularised time and physical time both start at zero in the
        augmented system; callers needing a shifted ``t0`` should pre-translate
        externally).

    Raises
    ------
    ValueError
        On invalid inputs (empty ``t_span``, negative span, bad regularisation).
    RuntimeError
        If the search fails to reach ``t_span[1]`` within an internal cap of
        ``s = 1e6``. This indicates either a malformed call or a trajectory
        that escapes to infinity in ``s`` before reaching the requested ``t``;
        callers should narrow the time horizon and retry.
    """
    reg = _validate_regularization(regularization)
    t0, tf = float(t_span[0]), float(t_span[1])
    if tf <= t0:
        raise ValueError(f"physical_to_regularized_span: require tf > t0, got ({t0}, {tf})")
    state0_arr = np.asarray(state0, dtype=np.float64)
    if state0_arr.shape != (6,):
        raise ValueError(
            f"physical_to_regularized_span: state0 must be a 6-vector, got shape {state0_arr.shape}"
        )
    mu_used = float(system.mu)
    duration = tf - t0

    # Cheap scale: dt/ds at the seed state. If r1 ~ r2 ~ O(1) this is O(1) and
    # ds ~ dt; near a primary it may be tiny, so we cap below.
    x0, y0, z0 = float(state0_arr[0]), float(state0_arr[1]), float(state0_arr[2])
    dtds0 = max(_dt_ds(x0, y0, z0, mu_used, reg), 1e-6)
    s_cap = max(1.0, duration / dtds0) * 100.0  # generous upper bracket
    s_cap = min(s_cap, 1.0e6)

    y_aug = np.empty(7, dtype=np.float64)
    y_aug[:6] = state0_arr
    y_aug[6] = t0

    def _rhs(s: float, y: NDArray[np.float64]) -> NDArray[np.float64]:
        return sundman_rhs(s, y, mu_used, reg)

    def _time_event(s: float, y: NDArray[np.float64]) -> float:
        return float(y[6]) - tf

    _time_event.terminal = True  # type: ignore[attr-defined]
    _time_event.direction = 1.0  # type: ignore[attr-defined]

    # initial_step / max_step: take caller hints if supplied; otherwise let
    # DOP853 choose. The default values are deliberately conservative.
    if initial_ds is not None and (initial_ds <= 0.0 or not math.isfinite(initial_ds)):
        raise ValueError(
            f"physical_to_regularized_span: initial_ds must be positive finite, got {initial_ds}"
        )
    if max_ds is not None and (max_ds <= 0.0 or not math.isfinite(max_ds)):
        raise ValueError(
            f"physical_to_regularized_span: max_ds must be positive finite, got {max_ds}"
        )

    sol = solve_ivp(
        _rhs,
        (0.0, s_cap),
        y_aug,
        method="DOP853",
        rtol=1e-9,
        atol=1e-9,
        events=_time_event,
        first_step=initial_ds,
        max_step=max_ds if max_ds is not None else np.inf,
    )
    if not sol.success:
        raise RuntimeError(
            f"physical_to_regularized_span: integrator failed at s={sol.t[-1]}: {sol.message}"
        )
    t_events = sol.t_events[0] if sol.t_events is not None else np.array([])
    if len(t_events) == 0:
        raise RuntimeError(
            f"physical_to_regularized_span: physical time did not reach {tf} within s_cap={s_cap}; "
            "increase the cap or narrow t_span"
        )
    return (0.0, float(t_events[0]))


def extract_perilune_distance(arc: RegularizedArc, system: CR3BPSystem) -> float:
    """Minimum r2 (secondary-relative distance) along a regularised arc.

    Tulip-orbit characterisation cares about the closest approach to the
    secondary (the Moon in Earth-Moon); this is the discrete minimum over the
    ``s``-grid stored in ``arc``. The regularised grid is naturally dense near
    perilune (because ``dt/ds`` shrinks there), so the discrete minimum is a
    good proxy for the true minimum without any interpolation.

    The returned value is in CR3BP nondimensional units of length; multiply by
    ``system.l_km`` to convert to kilometres.
    """
    states = arc.state_at_s
    if states.shape[0] != 6 or states.shape[1] == 0:
        raise ValueError(
            f"extract_perilune_distance: arc.state_at_s has shape {states.shape}; "
            "expected (6, N) with N >= 1"
        )
    x = states[0, :]
    y = states[1, :]
    z = states[2, :]
    # Secondary at (1 - mu, 0, 0).
    r2 = np.sqrt((x - (1.0 - system.mu)) ** 2 + y * y + z * z)
    return float(np.min(r2))
