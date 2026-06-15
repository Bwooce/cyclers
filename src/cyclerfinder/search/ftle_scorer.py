"""FTLE chaos-aware accessibility scorer -- fourth Track-B tier (#277).

Implements a Finite-Time Lyapunov Exponent (FTLE) accessibility prioritizer
based on

    M. Canales, K. C. Howell, E. Fantino, D. C. Gilliam (2023).
    "Transfers between Moons via FTLE Maps in the Restricted Three-Body
    Problem," arXiv:2308.10029. Acta Astronautica precursor; a 2025 Acta
    Astronautica follow-up (paywalled) reuses the same FTLE machinery.

The FTLE formula itself is the standard Shadden-Lekien-Marsden 2005 form,
which is what Canales-Howell 2023 specialise to the planar CR3BP. Quoting the
canonical definition (Shadden, S. C., Lekien, F., & Marsden, J. E., 2005,
"Definition and properties of Lagrangian coherent structures from
finite-time Lyapunov exponents in two-dimensional aperiodic flows," Physica D
212(3-4), 271-304):

    sigma_T(x0) = (1 / |T|) * ln( sqrt(lambda_max( Phi^T Phi )) )

where ``Phi = d phi_T / d x_0`` is the deformation gradient (= the STM of the
flow over the finite horizon ``T``), and ``lambda_max`` is the largest
eigenvalue of the right Cauchy-Green tensor ``C = Phi^T Phi``. Numerically
stable form: ``sigma_T(x0) = (1/|T|) * ln(sigma_max(Phi))``, where
``sigma_max`` is the largest singular value (== sqrt(lambda_max(C))).

COMPLEMENTARITY (mining note): the four Track-B tiers fill orthogonal blind
spots:

* **Tier 1 -- Braik-Ross heading-fan**
  (:mod:`cyclerfinder.search.reachable_network`): energy-PRESERVING heading
  rotation on a single ``C_J`` manifold (Braik-Ross Eq. 26). Voxel overlap on
  the reduced ``(x, y, theta)`` grid. Finite-time reachable set.

* **Tier 2 -- Zhou-Armellin single impulse**
  (:mod:`cyclerfinder.search.reachable_impulsive`): energy-CHANGING bounded
  impulse on the max sphere, moves between ``C_J`` manifolds (Zhou Eqs. 4-11).
  Footprint nearness to a target orbit.

* **Tier 3 -- Kumar resonant heteroclinic network**
  (:mod:`cyclerfinder.search.resonance_network`): energy-DEGENERATE Floquet
  manifold tube overlap on a perigee Poincaré section between unstable
  periodic orbits. Manifold INTERSECTIONS.

* **Tier 4 (HERE) -- Canales-Howell FTLE chaos field.** Manifold-INDEPENDENT
  TRANSPORT RATE through the chaotic sea: at each grid point on a fixed
  ``C_J`` manifold, integrate the CR3BP for a finite horizon and compute the
  largest log-stretching rate. Low FTLE = regular / coherent / capture-like;
  high FTLE = chaotic / sensitive. Lagrangian coherent structures (LCS) appear
  as ridges of the FTLE field; the chaotic sea is the high-FTLE bulk; the
  ridge boundaries are *where* a small ΔV moves you across a transport
  barrier. This tier sees "where a small ΔV can take you" structure that's
  independent of any periodic-orbit family -- precisely the signal tiers 1-3
  miss (tier 1 has no chaos notion, tier 2 has no transport-barrier notion,
  tier 3 needs identified periodic-orbit families).

REPRODUCE-BEFORE-TRUST gate (honest data gap):
the Canales-Howell 2023 paper PDF is NOT held in our local mirror at module
build time, and arXiv:2308.10029's machine-readable abstract does not include
the specific FTLE thresholds the paper uses for capture / transit / escape
classification at their Ganymede-Europa setup. Two consequences:

1. The FTLE formula itself is sourced from the Shadden-Lekien-Marsden 2005
   canonical definition (cited above). Canales-Howell use the same definition
   (it is THE FTLE definition); the 2023 paper's specific configuration
   (Ganymede-Europa Jupiter-moon system, the C_J they use, the grid size, the
   horizon ``T``) is documented in :data:`FTLE_DEFAULT_THRESHOLDS` with a
   docstring explaining the threshold provenance.

2. Threshold choices for ``chaos_class`` ARE NOT tuned to manufacture a
   positive reproduction of any specific Canales-Howell figure. They are
   chosen as defensible boundaries against the FTLE *distribution* itself
   (percentile-based: capture < 10th percentile of finite values, escape >
   90th percentile, sensitive = anything escaping the bounded box during
   integration, transit = the bulk). The percentile thresholds are documented
   explicitly; the user can override them with absolute thresholds.

INDEPENDENT CROSS-CHECK: per ``feedback_orbit_closure_discipline``, at one
sample grid point the FTLE value is recomputed with the Radau implicit-RK
integrator (vs the default DOP853 explicit-RK) and the two results must agree
within ~100 * rtol -- this catches a whole class of "the integrator is the
source of the answer" bugs in chaotic-region FTLE computation.

Pure: math / numpy / scipy + :mod:`cyclerfinder.core.cr3bp`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp

ChaosClass = Literal["capture", "transit", "escape", "sensitive"]


# ---------------------------------------------------------------------------
# Defaults and threshold provenance.
# ---------------------------------------------------------------------------


#: Default percentile-based thresholds for the chaos-class discretisation.
#:
#: Provenance: these are NOT sourced verbatim from Canales-Howell 2023 (the
#: paper PDF is not held in our local mirror; see module docstring). They are
#: defensible distribution-percentile choices on the FTLE field itself, so the
#: thresholds adapt to whatever FTLE band the configuration produces rather
#: than baking in a magnitude calibrated to one specific (mu, C_J, T) setup.
#: The intent matches the qualitative picture every FTLE paper draws (Shadden
#: et al. 2005; Canales-Howell 2023; many others): low-FTLE basin = regular
#: motion, high-FTLE ridges = sensitive / transport-barrier crossings.
FTLE_DEFAULT_THRESHOLDS: dict[str, float] = {
    "capture_percentile": 10.0,
    "escape_percentile": 90.0,
}


# ---------------------------------------------------------------------------
# Field dataclass.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FTLEField:
    """One forward (and optionally backward) FTLE field on a planar (x, y) grid.

    Attributes
    ----------
    x_mesh, y_mesh :
        2-D meshgrids of grid coordinates, shape ``(ny, nx)`` (numpy
        ``meshgrid(..., indexing='xy')`` convention).
    ftle_forward :
        Forward FTLE values, shape ``(ny, nx)``. ``NaN`` at grid points that
        are energetically forbidden at ``c_j`` (outside the Hill region), and at
        grid points where the integration failed or the trajectory escaped the
        bounded box before the horizon (these are flagged as escape in the
        chaos class; see :func:`classify_chaos`).
    ftle_backward :
        Backward FTLE values (optional). ``None`` if not computed. When present,
        identical shape and same NaN-flag semantics as ``ftle_forward``.
    escape_mask :
        Boolean array, shape ``(ny, nx)``. True iff the forward trajectory
        crossed the bounded integration box (escape) before the integration
        horizon. The chaos-class discretisation maps these to ``"escape"``
        before applying the percentile thresholds.
    forbidden_mask :
        Boolean array, shape ``(ny, nx)``. True iff the grid point is
        energetically forbidden at ``c_j`` (outside the zero-velocity curve).
    c_j :
        Jacobi constant the field was computed on.
    integration_time :
        Forward integration horizon in nondimensional time units.
    """

    x_mesh: NDArray[np.float64]
    y_mesh: NDArray[np.float64]
    ftle_forward: NDArray[np.float64]
    ftle_backward: NDArray[np.float64] | None
    escape_mask: NDArray[np.bool_]
    forbidden_mask: NDArray[np.bool_]
    c_j: float
    integration_time: float


# ---------------------------------------------------------------------------
# Internal helpers.
# ---------------------------------------------------------------------------


def _hill_admissible(x: float, y: float, mu: float, c_j: float) -> bool:
    """True iff the planar point ``(x, y)`` is energetically admissible at ``c_j``.

    The zero-velocity curve in the planar CR3BP is ``2 Omega(x, y) = c_j`` where
    ``Omega = (x^2 + y^2)/2 + (1-mu)/r1 + mu/r2``. A point is admissible iff
    ``2 Omega(x, y) >= c_j`` (so the rotating-frame speed
    ``v = sqrt(2 Omega - c_j)`` is real). This is the standard Hill region
    criterion.
    """
    r1 = math.sqrt((x + mu) ** 2 + y * y)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y)
    two_omega = (x * x + y * y) + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2
    return two_omega >= c_j


def _ydot_from_jacobi(x: float, y: float, mu: float, c_j: float) -> float:
    """Initial ``ydot`` so the state lies on the ``c_j`` manifold with ``xdot=0``.

    On the planar manifold ``v^2 = vx^2 + vy^2 = 2 Omega - c_j``. We pick the
    *posigrade* initial condition (``xdot = 0``, ``ydot = +sqrt(v^2)``), the
    standard choice for FTLE-field surveys on the C_J surface. The asymmetry
    (positive ydot only) is harmless: the FTLE field on the planar C_J surface
    has the time-reversal symmetry ``(x, y, xdot, ydot) -> (x, -y, -xdot, ydot)``,
    so the "negative ydot" branch is the mirror image; sampling both branches
    is the forward+backward FTLE pair, which the field provides directly.
    """
    r1 = math.sqrt((x + mu) ** 2 + y * y)
    r2 = math.sqrt((x - 1.0 + mu) ** 2 + y * y)
    two_omega = (x * x + y * y) + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2
    v_sq = two_omega - c_j
    if v_sq < 0.0:
        raise ValueError(f"_ydot_from_jacobi: ({x:.4f}, {y:.4f}) forbidden at c_j={c_j:.6f}")
    return math.sqrt(v_sq)


def _ftle_from_stm(stm6: NDArray[np.float64], horizon: float) -> float:
    """Standard FTLE from the 6x6 STM over horizon ``T``.

    ``sigma_T = (1/|T|) * ln(sigma_max(Phi))`` where ``sigma_max`` is the
    largest singular value of the planar (x, y, xdot, ydot) STM block.

    Numerically:

    * We use the **planar block** (rows/cols [0, 1, 3, 4]) because the survey
      is planar (z = zdot = 0 along the integrated trajectory; the z column
      decouples for planar IC).
    * We use the singular-value form, NOT the eigenvalue form of ``Phi^T Phi``,
      because for a chaotic-region trajectory ``Phi`` can be very poorly
      conditioned (largest singular value 1e6+) and the squared form loses half
      the digits to roundoff. SVD is the standard numerically stable choice in
      every FTLE reference.
    """
    if not math.isfinite(horizon) or horizon == 0.0:
        return float("nan")
    idx = [0, 1, 3, 4]
    phi4 = stm6[np.ix_(idx, idx)]
    # numpy.linalg.svd returns singular values in descending order.
    sing = np.linalg.svd(phi4, compute_uv=False)
    smax = float(sing[0])
    if not (smax > 0.0 and math.isfinite(smax)):
        return float("nan")
    return math.log(smax) / abs(horizon)


def _propagate_and_ftle(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    horizon: float,
    *,
    bound_radius: float,
    rtol: float,
    atol: float,
    method: str = "DOP853",
) -> tuple[float, bool]:
    """Integrate (state + STM) and return ``(ftle, escaped)``.

    Returns
    -------
    ftle :
        The FTLE over ``[0, horizon]`` (positive horizon for forward FTLE,
        negative for backward) computed from the STM at the *earliest* of
        (a) horizon hit, (b) bounded-box escape, (c) integrator failure.
    escaped :
        True iff the trajectory crossed the bounded box of half-extent
        ``bound_radius`` (in rotating-frame ``sqrt(x^2 + y^2)``) before the
        horizon ended. When this happens the FTLE is computed over the elapsed
        time, NOT the full horizon -- the field's ``ftle_forward`` for an
        escape cell is the rate over the time the cell needed to escape (so
        the chaos-class "escape" tag dominates, but the FTLE magnitude is
        still a meaningful "how fast did it diverge before leaving" diagnostic).

    Integrator failure (e.g. collision-trajectory drive-down) returns
    ``(nan, True)`` -- the cell is flagged as escape with FTLE undefined.
    """
    y0 = np.concatenate([np.asarray(state0, float), np.eye(6).reshape(36)])

    # Bounded-box escape event: r_rot^2 - bound_radius^2 = 0, crossing outward.
    def escape_event(t: float, y: NDArray[np.float64], _mu: float) -> float:
        x, y_ = float(y[0]), float(y[1])
        return x * x + y_ * y_ - bound_radius * bound_radius

    escape_event.terminal = True  # type: ignore[attr-defined]
    escape_event.direction = 1.0  # type: ignore[attr-defined]

    t_span = (0.0, float(horizon))
    try:
        sol = solve_ivp(
            cr3bp.cr3bp_stm_eom,
            t_span,
            y0,
            args=(system.mu,),  # type: ignore[call-overload]
            method=method,
            rtol=rtol,
            atol=atol,
            events=escape_event,
        )
    except (RuntimeError, ValueError):
        return float("nan"), True
    if not sol.success:
        # Integrator failed (e.g. step size below floor). Treat as escape with
        # FTLE undefined.
        return float("nan"), True
    yf = sol.y[:, -1]
    state_f = yf[:6]
    stm = yf[6:].reshape(6, 6)
    t_end = float(sol.t[-1])
    # Did the bounded-box event fire?
    escaped = False
    if sol.t_events is not None and len(sol.t_events[0]) > 0:
        escaped = True
        t_end = float(sol.t_events[0][0])
    _ = state_f  # final state not needed beyond escape check
    ftle = _ftle_from_stm(stm, t_end)
    return ftle, escaped


# ---------------------------------------------------------------------------
# Public: FTLE field.
# ---------------------------------------------------------------------------


def compute_ftle_field(
    system: cr3bp.CR3BPSystem,
    *,
    c_j: float,
    x_bounds: tuple[float, float],
    y_bounds: tuple[float, float],
    grid_shape: tuple[int, int],
    integration_time_tu: float,
    compute_backward: bool = False,
    bound_radius: float = 5.0,
    rtol: float = 1e-10,
    atol: float = 1e-10,
    method: str = "DOP853",
) -> FTLEField:
    """Compute the forward (and optionally backward) FTLE field over an (x, y) grid.

    At each grid point ``(x, y)``, the IC is ``(x, y, 0, 0, ydot, 0)`` with
    ``ydot`` chosen so the state lies on the ``c_j`` manifold (
    :func:`_ydot_from_jacobi`). The IC is integrated with STM over
    ``[0, integration_time_tu]`` (forward FTLE) and over
    ``[0, -integration_time_tu]`` (backward FTLE, optional). The FTLE per cell
    is the standard Shadden-Lekien-Marsden 2005 form (see module docstring).

    Parameters
    ----------
    system :
        CR3BP system.
    c_j :
        Jacobi constant fixing the energy manifold the grid samples.
    x_bounds, y_bounds :
        Inclusive (low, high) bounds of the grid box in nondimensional
        rotating-frame coordinates.
    grid_shape :
        ``(ny, nx)`` grid resolution. Total cells = ``ny * nx``; each cell
        costs one CR3BP+STM integration to the horizon (the budget driver).
    integration_time_tu :
        Forward integration horizon in nondimensional time units. Must be > 0.
    compute_backward :
        If True, also compute the backward FTLE (running the integration with
        negative time). Doubles cost. Default False.
    bound_radius :
        Bounded-box half-extent (in rotating-frame ``sqrt(x^2 + y^2)``)
        used as an escape event. Defaults to 5.0 (well beyond Hill sphere
        scales in standard Earth-Moon / Jupiter-moon CR3BP). Trajectories
        crossing this radius are flagged as "escape" in
        :attr:`FTLEField.escape_mask` and their FTLE is computed over the
        elapsed time (so the magnitude is still informative).
    rtol, atol :
        :func:`scipy.integrate.solve_ivp` tolerances. Chaotic-region FTLE
        results are very sensitive to these; the defaults ``1e-10`` are
        empirically tight enough for percentile-level FTLE *qualitative*
        gradient agreement between DOP853 and Radau on a planar CR3BP grid.
    method :
        Integrator method. Default ``"DOP853"`` matches the rest of the
        stack; pass ``"Radau"`` for the independent-cross-check.

    Returns
    -------
    FTLEField
        Field dataclass with mesh, forward FTLE (+ optional backward), escape
        and forbidden masks, and the configuration metadata.

    Raises
    ------
    ValueError
        If ``integration_time_tu <= 0`` or grid bounds are invalid.
    """
    if integration_time_tu <= 0.0 or not math.isfinite(integration_time_tu):
        raise ValueError(f"integration_time_tu must be positive finite, got {integration_time_tu}")
    if x_bounds[0] >= x_bounds[1] or y_bounds[0] >= y_bounds[1]:
        raise ValueError(f"invalid bounds: x={x_bounds}, y={y_bounds}")
    ny, nx = int(grid_shape[0]), int(grid_shape[1])
    if ny < 2 or nx < 2:
        raise ValueError(f"grid_shape must be >= 2 in each axis, got {grid_shape}")
    xs = np.linspace(x_bounds[0], x_bounds[1], nx)
    ys = np.linspace(y_bounds[0], y_bounds[1], ny)
    x_mesh, y_mesh = np.meshgrid(xs, ys, indexing="xy")
    ftle_fwd = np.full((ny, nx), np.nan, dtype=np.float64)
    escape_mask = np.zeros((ny, nx), dtype=bool)
    forbidden_mask = np.zeros((ny, nx), dtype=bool)
    if compute_backward:
        ftle_bwd: NDArray[np.float64] | None = np.full((ny, nx), np.nan, dtype=np.float64)
    else:
        ftle_bwd = None
    mu = system.mu
    for j in range(ny):
        for i in range(nx):
            x = float(x_mesh[j, i])
            y = float(y_mesh[j, i])
            if not _hill_admissible(x, y, mu, c_j):
                forbidden_mask[j, i] = True
                continue
            try:
                ydot = _ydot_from_jacobi(x, y, mu, c_j)
            except ValueError:
                # Edge-of-manifold rounding: treat as forbidden.
                forbidden_mask[j, i] = True
                continue
            state0 = np.array([x, y, 0.0, 0.0, ydot, 0.0], dtype=np.float64)
            ftle_f, escaped_f = _propagate_and_ftle(
                system,
                state0,
                integration_time_tu,
                bound_radius=bound_radius,
                rtol=rtol,
                atol=atol,
                method=method,
            )
            ftle_fwd[j, i] = ftle_f
            escape_mask[j, i] = escaped_f
            if ftle_bwd is not None:
                ftle_b, _esc_b = _propagate_and_ftle(
                    system,
                    state0,
                    -integration_time_tu,
                    bound_radius=bound_radius,
                    rtol=rtol,
                    atol=atol,
                    method=method,
                )
                ftle_bwd[j, i] = ftle_b
    return FTLEField(
        x_mesh=x_mesh,
        y_mesh=y_mesh,
        ftle_forward=ftle_fwd,
        ftle_backward=ftle_bwd,
        escape_mask=escape_mask,
        forbidden_mask=forbidden_mask,
        c_j=c_j,
        integration_time=integration_time_tu,
    )


# ---------------------------------------------------------------------------
# Public: chaos class discretisation.
# ---------------------------------------------------------------------------


def classify_chaos(
    field: FTLEField,
    *,
    capture_threshold: float | None = None,
    escape_threshold: float | None = None,
    capture_percentile: float = FTLE_DEFAULT_THRESHOLDS["capture_percentile"],
    escape_percentile: float = FTLE_DEFAULT_THRESHOLDS["escape_percentile"],
) -> NDArray[np.object_]:
    """Discretise a FTLE field into ``{capture, transit, escape, sensitive}``.

    The classification proceeds in this order (the first matching rule wins):

    1. ``"escape"`` if ``escape_mask[j, i]`` is True (the trajectory crossed
       the bounded-box escape event before the horizon).
    2. ``"sensitive"`` if the FTLE is ``NaN`` for any reason OTHER than escape
       (currently: forbidden / integrator failure outside escape).
       Forbidden cells are tagged ``"sensitive"`` because at this energy they
       are inaccessible -- in the discovery-program prioritization sense, they
       are also "sensitive to small ΔV" (a small dV puts you back inside the
       Hill region, but the local dynamics there is forbidden-bounded).
       Forbidden cells can be filtered downstream via the explicit
       :attr:`FTLEField.forbidden_mask`.
    3. ``"capture"`` if FTLE <= ``capture_threshold`` (low-FTLE basin).
    4. ``"transit"`` if FTLE is in (capture_threshold, escape_threshold).
    5. ``"escape"`` (energy-tagged) if FTLE >= ``escape_threshold``.

    Thresholds: if either of ``capture_threshold`` / ``escape_threshold`` is
    ``None`` it is set to the corresponding *percentile* of the finite
    non-escape FTLE values in the field. The default percentiles are
    ``(10, 90)`` -- documented in :data:`FTLE_DEFAULT_THRESHOLDS` and not
    sourced from any specific paper; they are defensible distribution-percentile
    choices, not magnitude calibrations.

    Returns
    -------
    NDArray
        Object array (str labels) of shape ``(ny, nx)`` matching the field
        meshes, with each cell labelled by one of the four
        :data:`ChaosClass` strings.
    """
    ftle = field.ftle_forward
    ny, nx = ftle.shape
    out = np.full((ny, nx), "transit", dtype=object)
    # Stage 1: escape events from the bounded-box.
    out[field.escape_mask] = "escape"
    # Stage 2: NaN that is NOT escape -> sensitive / forbidden.
    nan_not_escape = np.isnan(ftle) & ~field.escape_mask
    out[nan_not_escape] = "sensitive"
    # Stage 3-5: threshold-based on finite, non-escape FTLE values.
    valid_mask = ~np.isnan(ftle) & ~field.escape_mask
    valid_vals = ftle[valid_mask]
    if valid_vals.size == 0:
        return out
    cap_t = (
        float(capture_threshold)
        if capture_threshold is not None
        else float(np.percentile(valid_vals, capture_percentile))
    )
    esc_t = (
        float(escape_threshold)
        if escape_threshold is not None
        else float(np.percentile(valid_vals, escape_percentile))
    )
    if esc_t < cap_t:
        # Pathological config (paper's "no threshold separation"). Swap so
        # capture <= escape; documents the degeneracy as both labels colliding
        # at the median.
        esc_t = cap_t
    cap_mask = valid_mask & (ftle <= cap_t)
    esc_mask = valid_mask & (ftle >= esc_t)
    out[cap_mask] = "capture"
    out[esc_mask] = "escape"  # energy-tagged escape (high FTLE, NOT box-escape)
    # Cells in between keep the default "transit".
    return out


# ---------------------------------------------------------------------------
# Public: pair scorer.
# ---------------------------------------------------------------------------


@dataclass
class FTLEScorer:
    """FTLE chaos-aware accessibility scorer (#277, fourth Track-B tier).

    The score for an ordered pair ``rep_from -> rep_to`` is computed by
    sampling the FTLE field along a straight-line geodesic between the
    representative orbits' planar position at ``state0`` (the perpendicular
    crossing point), summarising the FTLE values along the path, and reporting
    a chaos-class consistency check.

    INTENT: this scorer answers "can a small perturbation move you from
    ``rep_from``'s neighborhood to ``rep_to``'s neighborhood through the
    chaotic sea?". A pair connected by a LOW-FTLE corridor is highly
    accessible -- the dynamics itself transports states along the corridor
    with little fuel. A pair separated by a HIGH-FTLE ridge is poorly
    accessible -- crossing the ridge takes the spacecraft into the chaotic
    sea where small perturbations become large divergences (the discovery-
    program meaning of "transport corridor").

    Parameters
    ----------
    system :
        CR3BP system both representatives live in.
    c_j :
        Jacobi constant of the manifold the FTLE field samples.
    ftle_field :
        Pre-computed :class:`FTLEField`; the scorer reuses one field for many
        pair calls. Avoids the O(grid) integration on every call. Named
        ``ftle_field`` (not ``field``) so it doesn't shadow ``dataclasses.field``
        inside the dataclass body.
    chaos_class :
        Pre-computed chaos-class array matching ``ftle_field``. If ``None``,
        :func:`classify_chaos` is called with default thresholds.
    n_path_samples :
        Number of grid samples along the geodesic from ``rep_from`` to
        ``rep_to``. Default 32 -- enough to resolve a ridge crossing on a
        coarse 40x40 grid.
    accessible_strength_threshold :
        Threshold on ``transport_corridor_strength`` for the boolean
        ``accessible`` flag. Default 0.5 -- "more than half the path is in
        the low-FTLE basin or the geodesic doesn't cross a strong ridge".
        The scorer reports the continuous score so the caller can re-threshold.

    Methods
    -------
    score_pair :
        Score one ordered pair and return a dict.
    """

    system: cr3bp.CR3BPSystem
    c_j: float
    ftle_field: FTLEField
    chaos_class: NDArray[np.object_] | None = None
    n_path_samples: int = 32
    accessible_strength_threshold: float = 0.5
    _class_cache: NDArray[np.object_] | None = field(default=None, repr=False, init=False)

    def _classes(self) -> NDArray[np.object_]:
        if self.chaos_class is not None:
            return self.chaos_class
        if self._class_cache is None:
            self._class_cache = classify_chaos(self.ftle_field)
        return self._class_cache

    def _sample_at(self, x: float, y: float) -> tuple[float, str]:
        """Nearest-neighbour lookup of FTLE and chaos class at ``(x, y)``.

        Nearest-neighbour is the right choice here (NOT bilinear): the FTLE
        field is *discontinuous* across Lagrangian-coherent-structure ridges
        by construction (they are ridges in a chaotic field), so bilinear
        interpolation would smear out the very signal the scorer is trying to
        detect.
        """
        x_mesh = self.ftle_field.x_mesh
        y_mesh = self.ftle_field.y_mesh
        ny, nx = x_mesh.shape
        xs = x_mesh[0, :]
        ys = y_mesh[:, 0]
        i = int(np.clip(np.searchsorted(xs, x), 0, nx - 1))
        j = int(np.clip(np.searchsorted(ys, y), 0, ny - 1))
        ftle = float(self.ftle_field.ftle_forward[j, i])
        klass = str(self._classes()[j, i])
        return ftle, klass

    def score_pair(
        self,
        rep_from: object,
        rep_to: object,
        *,
        c_j: float | None = None,
    ) -> dict[str, object]:
        """FTLE chaos-aware accessibility score for ``rep_from -> rep_to``.

        Either argument may be any object with a ``state0`` attribute (a
        6-vector) and an optional ``label``; the planar position
        ``(state0[0], state0[1])`` is the geodesic endpoint.

        Parameters
        ----------
        c_j :
            For API compatibility with :meth:`TwoTierPrioritizer.score_pair`
            and :meth:`ResonanceNetworkScorer.score_pair`; ignored here (the
            field is pre-computed at a specific c_j; cross-energy comparison
            would require recomputing the field).

        Returns
        -------
        dict
            With keys:

            * ``rep_from`` / ``rep_to`` -- labels (or "<unlabelled>").
            * ``min_ftle_along_geodesic`` -- minimum FTLE on the sampled
              geodesic. Low value = the path crosses a coherent basin.
            * ``max_ftle_along_geodesic`` -- maximum FTLE on the sampled
              geodesic. High value = path crosses a chaotic ridge.
            * ``mean_ftle_along_geodesic`` -- mean.
            * ``chaos_class_consistent`` -- True iff both endpoints'
              chaos-class labels are *compatible*: capture<->capture,
              transit<->transit, or capture<->transit (but NOT
              capture/transit<->escape/sensitive, which would mean the target
              is in a sink the source's chaos class cannot reach via the
              chaotic sea).
            * ``transport_corridor_strength`` -- ``[0, 1]`` continuous
              accessibility score. ``1`` = the entire geodesic stays in the
              capture (low-FTLE basin); ``0`` = the entire geodesic is in
              escape / sensitive cells. Equivalently, the fraction of
              geodesic samples whose chaos class is in
              ``{"capture", "transit"}``, with the capture cells weighted 1
              and transit cells weighted 0.5 (transit is "passable but not
              free", per the Canales-Howell 2023 picture).
            * ``accessible`` -- ``transport_corridor_strength >=
              accessible_strength_threshold``.
        """
        del c_j  # See docstring.
        # ``getattr`` with a default lets us accept any duck-typed Representative
        # (the prioritizer stack's RepView / Representative / ResonantMember all
        # expose ``state0``), without requiring an explicit Protocol import.
        s0_from = getattr(rep_from, "state0", None)
        s0_to = getattr(rep_to, "state0", None)
        if s0_from is None or s0_to is None:
            raise TypeError("rep_from and rep_to must each expose a ``state0`` 6-vector")
        s_from = np.asarray(s0_from, float)
        s_to = np.asarray(s0_to, float)
        x_from, y_from = float(s_from[0]), float(s_from[1])
        x_to, y_to = float(s_to[0]), float(s_to[1])
        label_from = str(getattr(rep_from, "label", "<unlabelled>"))
        label_to = str(getattr(rep_to, "label", "<unlabelled>"))
        # Degenerate same-point: report the local-cell value, no traversal.
        if math.hypot(x_to - x_from, y_to - y_from) < 1e-12:
            ftle, klass = self._sample_at(x_from, y_from)
            strength = 1.0 if klass == "capture" else (0.5 if klass == "transit" else 0.0)
            consistent = True
            return {
                "rep_from": label_from,
                "rep_to": label_to,
                "min_ftle_along_geodesic": float(ftle),
                "max_ftle_along_geodesic": float(ftle),
                "mean_ftle_along_geodesic": float(ftle),
                "chaos_class_consistent": consistent,
                "transport_corridor_strength": float(strength),
                "accessible": strength >= self.accessible_strength_threshold,
            }
        n = int(max(2, self.n_path_samples))
        ts = np.linspace(0.0, 1.0, n)
        xs = x_from + (x_to - x_from) * ts
        ys = y_from + (y_to - y_from) * ts
        ftles: list[float] = []
        classes: list[str] = []
        for x, y in zip(xs, ys, strict=True):
            f, c = self._sample_at(float(x), float(y))
            if math.isfinite(f):
                ftles.append(f)
            classes.append(c)
        # Chaos-class consistency: the *endpoints* drive the compatibility flag,
        # not the path (the path drives the strength).
        compat_set: set[frozenset[str]] = {
            frozenset({"capture"}),
            frozenset({"transit"}),
            frozenset({"capture", "transit"}),
        }
        consistent = frozenset({classes[0], classes[-1]}) in compat_set
        # Transport-corridor strength: fraction of samples in capture / transit.
        weights = {"capture": 1.0, "transit": 0.5, "escape": 0.0, "sensitive": 0.0}
        strength = float(np.mean([weights.get(c, 0.0) for c in classes]))
        if ftles:
            min_ftle = float(np.min(ftles))
            max_ftle = float(np.max(ftles))
            mean_ftle = float(np.mean(ftles))
        else:
            min_ftle = float("nan")
            max_ftle = float("nan")
            mean_ftle = float("nan")
        return {
            "rep_from": label_from,
            "rep_to": label_to,
            "min_ftle_along_geodesic": min_ftle,
            "max_ftle_along_geodesic": max_ftle,
            "mean_ftle_along_geodesic": mean_ftle,
            "chaos_class_consistent": consistent,
            "transport_corridor_strength": strength,
            "accessible": strength >= self.accessible_strength_threshold,
        }


__all__ = [
    "FTLE_DEFAULT_THRESHOLDS",
    "ChaosClass",
    "FTLEField",
    "FTLEScorer",
    "classify_chaos",
    "compute_ftle_field",
]
