r"""Flyby-continuity-wired FBS-analytic ΔV-chain optimiser (#244, Stage 2).

This extends the #243 match-point NLP (:mod:`cyclerfinder.search.fbs_optimize`)
with the patched-conic flyby-continuity constraints (Ellison Eqs. 3-4) at every
interior body — the #243 caveat-2 wiring — and supplies their ANALYTIC gradients
via :func:`cyclerfinder.core.fbs_match_point.flyby_coupling_block`.

What changes versus #243's :func:`optimize_chain_fbs`
----------------------------------------------------
#243's NLP shares ONE heliocentric boundary velocity ``v_j`` at every interior
body, so heliocentric velocity is continuous by construction and a real turning
flyby cannot be represented (v∞ in == v∞ out trivially). Wiring a patched-conic
flyby means the interior body carries SEPARATE arrival/departure heliocentric
velocities ``v_arr_j`` and ``v_dep_j`` (the spacecraft turns its v∞ vector across
the flyby), coupled by the two Ellison continuity constraints:

* ``c_vinf_j  = ‖v_dep_j - v_planet_j‖ - ‖v_arr_j - v_planet_j‖ = 0`` (Eq. 3,
  equality — an unpowered flyby conserves v∞ magnitude);
* ``c_alt_j   = r_periapse(v∞) - r_safe >= 0`` (Eq. 4, inequality — the bend is
  ballistically realizable above the body's safe altitude).

Decision vector
---------------
``x = [ Δv_0..Δv_{M-1} (3M) | v_dep_0 (3) | (v_arr_j, v_dep_j) for j=1..M-1
        (6(M-1)) | v_arr_M (3) ]``

i.e. the departure node 0 has a departure velocity only, the arrival node M an
arrival velocity only, and each interior node ``j`` (0<j<M) BOTH. Leg ``i`` uses
the departure velocity at node ``i`` as its left-boundary ``v0`` and the arrival
velocity at node ``i+1`` as its right-boundary ``vf``.

Constraint stack handed to SLSQP (all row-scaled like #243):
``[ chain match-point defect (6M) ; c_vinf (M-1) ; c_alt (M-1) ]`` with the
match-point block driven to 0, the v∞-continuity block to 0, and the altitude
block kept ``>= 0``.

The analytic Jacobian blends the FBS per-leg STM columns (match-point rows; the
#226 :func:`match_point_defect_jacobian`) with the
:func:`flyby_coupling_block` gradients (continuity / altitude rows). The only
available correctness check is FD-vs-analytic on the whole stacked constraint
(Ellison publishes no numeric gradient) — see ``test_fbs_optimize_flyby.py``.

NO catalogue writeback. Method evaluation only. Additive: the #243 module is
untouched.

Source: D. H. Ellison et al., "Analytic Gradient Computation for Bounded-Impulse
Trajectory Models Using Two-Sided Shooting," JGCD 41(7), 2018,
doi:10.2514/1.G003077 (Eqs. 2-4, 31-32, 42, A1-A6).
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import NonlinearConstraint, minimize

from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2
from cyclerfinder.core.fbs_match_point import (
    FbsLeg,
    chain_defect,
    chain_defect_jacobian,
    flyby_coupling_block,
)
from cyclerfinder.core.kepler import KeplerError
from cyclerfinder.search.fbs_optimize import ChainLegSpec

Vec3 = NDArray[np.float64]

_DEFECT_PENALTY: float = 1.0e3


@dataclass(frozen=True)
class FlybyChainSpec:
    """Fixed (non-decision) geometry of a flyby chain.

    Attributes
    ----------
    legs:
        ``M`` :class:`ChainLegSpec` (boundary positions, ToF, alpha). FIXED.
    bodies:
        The ``M+1`` body codes at each node (``len == M+1``). The interior
        bodies (indices ``1..M-1``) are the flyby bodies whose continuity /
        altitude constraints are wired.
    v_planets:
        The ``M+1`` heliocentric planet velocities at each node's epoch (km/s),
        FIXED. Used to form ``v∞ = v_sc - v_planet`` at each interior flyby.
    """

    legs: tuple[ChainLegSpec, ...]
    bodies: tuple[str, ...]
    v_planets: tuple[Vec3, ...]

    def __post_init__(self) -> None:
        m = len(self.legs)
        if m < 1:
            raise ValueError("chain must contain at least one leg")
        if len(self.bodies) != m + 1:
            raise ValueError(f"bodies must have {m + 1} entries, got {len(self.bodies)}")
        if len(self.v_planets) != m + 1:
            raise ValueError(f"v_planets must have {m + 1} entries, got {len(self.v_planets)}")


@dataclass(frozen=True)
class FbsFlybyOptimizeResult:
    """Outcome of one flyby-wired FBS-gradient (or FD) chain ΔV optimisation."""

    total_dv_kms: float
    dv_per_leg_kms: tuple[float, ...]
    dvs: tuple[Vec3, ...]
    v_dep_per_node: tuple[Vec3, ...]
    v_arr_per_node: tuple[Vec3, ...]
    vinf_per_node_kms: tuple[float, ...]
    max_defect: float
    max_cvinf: float
    min_calt_km: float
    feasible: bool
    flyby_feasible: bool
    success: bool
    nfev: int
    njev: int
    nit: int
    constr_nfev: int
    constr_njev: int
    wall_s: float
    used_analytic_jac: bool


# ---------------------------------------------------------------------------
# Decision-vector layout helpers
# ---------------------------------------------------------------------------


def _layout(m: int) -> tuple[int, int]:
    """Return ``(n_dv, n_total)`` for an ``M``-leg flyby chain.

    ``n_dv = 3M`` impulse components; velocity block is
    ``3 (dep node 0) + 6 (M-1 interior arr+dep) + 3 (arr node M) = 6M`` components,
    so ``n_total = 3M + 6M = 9M``.
    """
    n_dv = 3 * m
    n_total = 9 * m
    return n_dv, n_total


def _vel_offset(m: int, node: int, *, arr: bool) -> int:
    """Column offset (within the velocity block) of node ``node``'s arr/dep slot.

    Velocity block ordering: ``[v_dep_0 | v_arr_1 v_dep_1 | ... | v_arr_M]``.
    Node 0 has dep only; node M has arr only; interior nodes have ``arr`` then
    ``dep``. Returns the absolute column index into the full decision vector.
    """
    n_dv = 3 * m
    if node == 0:
        if arr:
            raise ValueError("departure node 0 has no arrival slot")
        return n_dv
    # node 0 occupies 3; each interior node occupies 6 (arr, dep)
    base = n_dv + 3 + 6 * (node - 1)
    if node == m:
        if not arr:
            raise ValueError("arrival node M has no departure slot")
        return base
    return base if arr else base + 3


def _pack_x0(
    spec: FlybyChainSpec,
    dvs: list[Vec3] | tuple[Vec3, ...],
    varr: list[Vec3] | tuple[Vec3, ...],
    vdep: list[Vec3] | tuple[Vec3, ...],
) -> NDArray[np.float64]:
    """Pack seeds into the flat decision vector.

    ``varr`` carries the arrival velocities at nodes ``1..M`` (length ``M``);
    ``vdep`` carries the departure velocities at nodes ``0..M-1`` (length ``M``).
    """
    m = len(spec.legs)
    if len(dvs) != m:
        raise ValueError(f"dvs must have {m} entries, got {len(dvs)}")
    if len(varr) != m:
        raise ValueError(f"varr must have {m} entries (nodes 1..M), got {len(varr)}")
    if len(vdep) != m:
        raise ValueError(f"vdep must have {m} entries (nodes 0..M-1), got {len(vdep)}")
    _, n_total = _layout(m)
    x = np.zeros(n_total, dtype=np.float64)
    for i in range(m):
        x[3 * i : 3 * i + 3] = np.asarray(dvs[i], dtype=np.float64)
    # departure node 0
    x[_vel_offset(m, 0, arr=False) : _vel_offset(m, 0, arr=False) + 3] = np.asarray(
        vdep[0], dtype=np.float64
    )
    # interior nodes 1..M-1: arr then dep
    for node in range(1, m):
        a = _vel_offset(m, node, arr=True)
        d = _vel_offset(m, node, arr=False)
        x[a : a + 3] = np.asarray(varr[node - 1], dtype=np.float64)
        x[d : d + 3] = np.asarray(vdep[node], dtype=np.float64)
    # arrival node M
    a = _vel_offset(m, m, arr=True)
    x[a : a + 3] = np.asarray(varr[m - 1], dtype=np.float64)
    return x


def _split(
    spec: FlybyChainSpec, x: NDArray[np.float64]
) -> tuple[list[Vec3], list[Vec3], list[Vec3]]:
    """Inverse of :func:`_pack_x0`: ``(dvs, v_dep_per_node, v_arr_per_node)``.

    ``v_dep_per_node`` has length ``M`` (nodes 0..M-1); ``v_arr_per_node`` has
    length ``M`` (nodes 1..M). Leg ``i`` uses ``v_dep_per_node[i]`` as ``v0`` and
    ``v_arr_per_node[i]`` as ``vf``.
    """
    m = len(spec.legs)
    dvs = [np.asarray(x[3 * i : 3 * i + 3], dtype=np.float64) for i in range(m)]
    v_dep: list[Vec3] = []
    for i in range(m):
        off = _vel_offset(m, i, arr=False)
        v_dep.append(np.asarray(x[off : off + 3], dtype=np.float64))
    v_arr: list[Vec3] = []
    for node in range(1, m + 1):
        off = _vel_offset(m, node, arr=True)
        v_arr.append(np.asarray(x[off : off + 3], dtype=np.float64))
    return dvs, v_dep, v_arr


def _fbs_legs(
    spec: FlybyChainSpec, v_dep: list[Vec3], v_arr: list[Vec3], mu: float
) -> tuple[FbsLeg, ...]:
    """Build FBS legs: leg i = (r0, v_dep[i]) -> (rf, v_arr[i])."""
    return tuple(
        FbsLeg(
            r0=leg.r0,
            v0=np.asarray(v_dep[i], dtype=np.float64),
            rf=leg.rf,
            vf=np.asarray(v_arr[i], dtype=np.float64),
            tof_s=leg.tof_s,
            alpha=leg.alpha,
            mu=mu,
        )
        for i, leg in enumerate(spec.legs)
    )


def _defect_scale(spec: FlybyChainSpec, mu: float) -> NDArray[np.float64]:
    """Per-row non-dimensionalisation of the match-point defect (pos/AU, vel/v_circ)."""
    rows: list[float] = []
    for leg in spec.legs:
        v_scale = float(np.sqrt(mu / float(np.linalg.norm(leg.r0))))
        rows.extend([AU_KM, AU_KM, AU_KM, v_scale, v_scale, v_scale])
    return np.asarray(rows, dtype=np.float64)


def _flyby_scale(
    spec: FlybyChainSpec, mu: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Row scales for the (c_vinf, c_alt) blocks: vel by v_circ, altitude by AU."""
    m = len(spec.legs)
    n_int = m - 1
    vinf_scale = np.empty(n_int, dtype=np.float64)
    for node in range(1, m):
        r0 = spec.legs[node - 1].rf
        vinf_scale[node - 1] = float(np.sqrt(mu / float(np.linalg.norm(r0))))
    alt_scale = np.full(n_int, AU_KM, dtype=np.float64)
    return vinf_scale, alt_scale


# ---------------------------------------------------------------------------
# Constraint vector + analytic Jacobian
# ---------------------------------------------------------------------------


def _flyby_constraint_vector(
    spec: FlybyChainSpec, x: NDArray[np.float64], *, mu: float = MU_SUN_KM3_S2
) -> NDArray[np.float64]:
    """Stacked, row-scaled constraint ``[defect (6M) ; c_vinf (M-1) ; c_alt (M-1)]``.

    The match-point and v∞-continuity blocks are equalities (driven to 0); the
    altitude block is the inequality slack (kept ``>= 0``). Row-scaled so SLSQP
    sees a well-conditioned constraint and a single ``feas_tol`` is meaningful.
    """
    m = len(spec.legs)
    n_int = m - 1
    dvs, v_dep, v_arr = _split(spec, x)
    scale = _defect_scale(spec, mu)
    vinf_scale, alt_scale = _flyby_scale(spec, mu)

    fbs_legs = _fbs_legs(spec, v_dep, v_arr, mu)
    try:
        defect = chain_defect(fbs_legs, tuple(dvs)) / scale
    except KeplerError:
        defect = np.full(6 * m, _DEFECT_PENALTY, dtype=np.float64)

    c_vinf = np.zeros(n_int, dtype=np.float64)
    c_alt = np.zeros(n_int, dtype=np.float64)
    for node in range(1, m):
        v_in = v_arr[node - 1]  # arrival at this interior node
        v_out = v_dep[node]  # departure from this interior node
        v_pl = np.asarray(spec.v_planets[node], dtype=np.float64)
        try:
            blk = flyby_coupling_block(v_in, v_out, v_pl, spec.bodies[node])
            c_vinf[node - 1] = blk.c_vinf_kms / vinf_scale[node - 1]
            c_alt[node - 1] = blk.c_altitude_km / alt_scale[node - 1]
        except (ValueError, KeplerError):
            c_vinf[node - 1] = _DEFECT_PENALTY
            c_alt[node - 1] = -_DEFECT_PENALTY
    return np.concatenate([defect, c_vinf, c_alt])


def _flyby_constraint_jac_analytic(
    spec: FlybyChainSpec, x: NDArray[np.float64], *, mu: float = MU_SUN_KM3_S2
) -> NDArray[np.float64]:
    """Analytic Jacobian of :func:`_flyby_constraint_vector` (row-scaled).

    Match-point rows: the #226 ``chain_defect_jacobian`` expressed in this
    layout's columns (each leg's ``v0`` -> the node's departure slot, ``vf`` ->
    the node's arrival slot). Flyby rows: the ``flyby_coupling_block`` analytic
    gradients placed in the interior node's arr / dep columns.
    """
    m = len(spec.legs)
    n_int = m - 1
    _, n_total = _layout(m)
    n_rows = 6 * m + 2 * n_int
    dvs, v_dep, v_arr = _split(spec, x)
    scale = _defect_scale(spec, mu)
    vinf_scale, alt_scale = _flyby_scale(spec, mu)

    jac = np.zeros((n_rows, n_total), dtype=np.float64)

    # --- match-point defect rows (the #226 per-leg STM columns) ---
    fbs_legs = _fbs_legs(spec, v_dep, v_arr, mu)
    try:
        # chain_defect_jacobian returns (6M) x (3M + 3(M+1)) in the SHARED-v layout
        # [Δv | v_0..v_M]; we re-place its v columns into this layout's dep/arr slots.
        from cyclerfinder.core.fbs_match_point import match_point_defect_jacobian

        for i, (leg, dv) in enumerate(zip(fbs_legs, dvs, strict=True)):
            per_leg = match_point_defect_jacobian(leg, dv)  # 6x9 [Δv | v0 | vf]
            row = 6 * i
            # ∂c_i/∂Δv_i
            jac[row : row + 6, 3 * i : 3 * i + 3] = per_leg[:, 0:3]
            # ∂c_i/∂v0 = leg i departure velocity at node i
            d = _vel_offset(m, i, arr=False)
            jac[row : row + 6, d : d + 3] = per_leg[:, 3:6]
            # ∂c_i/∂vf = leg i arrival velocity at node i+1
            a = _vel_offset(m, i + 1, arr=True)
            jac[row : row + 6, a : a + 3] = per_leg[:, 6:9]
        _ = chain_defect_jacobian  # documented provenance of the block structure
    except KeplerError:
        return np.zeros((n_rows, n_total), dtype=np.float64)
    # row-scale the match-point block
    jac[: 6 * m, :] /= scale[:, None]

    # --- flyby continuity / altitude rows ---
    for node in range(1, m):
        v_in = v_arr[node - 1]
        v_out = v_dep[node]
        v_pl = np.asarray(spec.v_planets[node], dtype=np.float64)
        a_col = _vel_offset(m, node, arr=True)
        d_col = _vel_offset(m, node, arr=False)
        try:
            blk = flyby_coupling_block(v_in, v_out, v_pl, spec.bodies[node])
        except (ValueError, KeplerError):
            continue
        r_vinf = 6 * m + (node - 1)
        r_alt = 6 * m + n_int + (node - 1)
        sv = vinf_scale[node - 1]
        sa = alt_scale[node - 1]
        jac[r_vinf, a_col : a_col + 3] = blk.d_cvinf_d_v_arr / sv
        jac[r_vinf, d_col : d_col + 3] = blk.d_cvinf_d_v_dep / sv
        jac[r_alt, a_col : a_col + 3] = blk.d_calt_d_v_arr / sa
        jac[r_alt, d_col : d_col + 3] = blk.d_calt_d_v_dep / sa
    return jac


def _flyby_constraint_jac_fd(
    spec: FlybyChainSpec,
    x: NDArray[np.float64],
    *,
    analytic: bool = True,
    mu: float = MU_SUN_KM3_S2,
    h: float = 1.0e-5,
) -> NDArray[np.float64]:
    """Return the constraint Jacobian — analytic (``analytic=True``) or central FD.

    A single entry point so the test can request either lane and compare. The FD
    lane central-differences :func:`_flyby_constraint_vector` column-by-column.
    """
    if analytic:
        return _flyby_constraint_jac_analytic(spec, x, mu=mu)
    n = x.size
    c0 = _flyby_constraint_vector(spec, x, mu=mu)
    jac = np.zeros((c0.size, n), dtype=np.float64)
    for k in range(n):
        xp = x.copy()
        xm = x.copy()
        step = h * max(abs(float(x[k])), 1.0)
        xp[k] += step
        xm[k] -= step
        jac[:, k] = (
            _flyby_constraint_vector(spec, xp, mu=mu) - _flyby_constraint_vector(spec, xm, mu=mu)
        ) / (2.0 * step)
    return jac


# ---------------------------------------------------------------------------
# The optimiser
# ---------------------------------------------------------------------------


def optimize_chain_fbs_flyby(
    spec: FlybyChainSpec,
    dv0_per_leg: list[Vec3] | tuple[Vec3, ...],
    varr0: list[Vec3] | tuple[Vec3, ...],
    vdep0: list[Vec3] | tuple[Vec3, ...],
    *,
    mu: float = MU_SUN_KM3_S2,
    rendezvous_vplanet: Vec3 | None = None,
    use_analytic_jac: bool = True,
    feas_tol: float = 1.0e-8,
    flyby_tol: float = 1.0e-6,
    maxiter: int = 200,
    ftol: float = 1.0e-9,
) -> FbsFlybyOptimizeResult:
    r"""Minimise chain ΔV with match-point + patched-conic flyby-continuity constraints.

    The objective is ``Σ‖Δv_i‖`` (plus the terminal arrival v∞ when
    ``rendezvous_vplanet``). The equality constraints are the stacked match-point
    defect (driven to 0) AND the per-flyby v∞-magnitude continuity (driven to 0);
    the inequality constraints are the per-flyby periapsis-altitude feasibilities
    (kept ``>= 0``). The constraint Jacobian is the ANALYTIC blend
    (:func:`_flyby_constraint_jac_analytic`) when ``use_analytic_jac`` (the FBS
    lane) or SLSQP finite-difference when not (the FD baseline) — the ONLY
    difference between the lanes is the gradient source.

    ``varr0`` seeds arrival velocities at nodes ``1..M`` (length ``M``); ``vdep0``
    seeds departures at nodes ``0..M-1`` (length ``M``).
    """
    m = len(spec.legs)
    n_int = m - 1
    _, n_total = _layout(m)
    counts = {"obj": 0, "grad": 0, "con": 0, "cjac": 0}

    def _obj(x: NDArray[np.float64]) -> float:
        counts["obj"] += 1
        dvs, _, v_arr = _split(spec, x)
        total = float(sum(float(np.linalg.norm(dv)) for dv in dvs))
        if rendezvous_vplanet is not None:
            vinf_arr = v_arr[-1] - np.asarray(rendezvous_vplanet, dtype=np.float64)
            total += float(np.linalg.norm(vinf_arr))
        return total

    def _obj_grad(x: NDArray[np.float64]) -> NDArray[np.float64]:
        counts["grad"] += 1
        dvs, _, v_arr = _split(spec, x)
        g = np.zeros(n_total, dtype=np.float64)
        for i, dv in enumerate(dvs):
            nrm = float(np.linalg.norm(dv))
            if nrm > 0.0:
                g[3 * i : 3 * i + 3] = dv / nrm
        if rendezvous_vplanet is not None:
            vinf_arr = v_arr[-1] - np.asarray(rendezvous_vplanet, dtype=np.float64)
            nrm = float(np.linalg.norm(vinf_arr))
            if nrm > 0.0:
                a = _vel_offset(m, m, arr=True)
                g[a : a + 3] = vinf_arr / nrm
        return g

    # --- equality block: match-point defect + v∞ continuity (rows 0 .. 6M+n_int) ---
    def _eq(x: NDArray[np.float64]) -> NDArray[np.float64]:
        counts["con"] += 1
        c = _flyby_constraint_vector(spec, x, mu=mu)
        return c[: 6 * m + n_int]

    def _eq_jac(x: NDArray[np.float64]) -> NDArray[np.float64]:
        counts["cjac"] += 1
        j = _flyby_constraint_jac_analytic(spec, x, mu=mu)
        return j[: 6 * m + n_int, :]

    # --- inequality block: altitude (rows 6M+n_int .. end), >= 0 ---
    def _ineq(x: NDArray[np.float64]) -> NDArray[np.float64]:
        c = _flyby_constraint_vector(spec, x, mu=mu)
        return c[6 * m + n_int :]

    def _ineq_jac(x: NDArray[np.float64]) -> NDArray[np.float64]:
        j = _flyby_constraint_jac_analytic(spec, x, mu=mu)
        return j[6 * m + n_int :, :]

    eq_jac: Any = _eq_jac if use_analytic_jac else "2-point"
    constraints: list[Any] = [NonlinearConstraint(_eq, 0.0, 0.0, jac=eq_jac)]
    if n_int > 0:
        ineq_jac: Any = _ineq_jac if use_analytic_jac else "2-point"
        constraints.append(NonlinearConstraint(_ineq, 0.0, np.inf, jac=ineq_jac))

    x0 = _pack_x0(spec, list(dv0_per_leg), list(varr0), list(vdep0))

    t_start = time.perf_counter()
    res: Any = minimize(
        _obj,
        x0,
        jac=_obj_grad,
        method="SLSQP",
        constraints=constraints,
        options={"maxiter": maxiter, "ftol": ftol},
    )
    wall = time.perf_counter() - t_start

    x_sol = np.asarray(res.x, dtype=np.float64)
    dvs_sol, v_dep_sol, v_arr_sol = _split(spec, x_sol)
    c_sol = _flyby_constraint_vector(spec, x_sol, mu=mu)
    defect = c_sol[: 6 * m]
    c_vinf = c_sol[6 * m : 6 * m + n_int]
    c_alt = c_sol[6 * m + n_int :]
    max_defect = float(np.max(np.abs(defect)))
    max_cvinf = float(np.max(np.abs(c_vinf))) if n_int else 0.0
    min_calt = float(np.min(c_alt)) if n_int else float("inf")

    # emerged v∞ at each node (arrival side for nodes 1..M, departure for node 0)
    vinf_nodes: list[float] = []
    vdep0_node = v_dep_sol[0] - np.asarray(spec.v_planets[0], dtype=np.float64)
    vinf_nodes.append(float(np.linalg.norm(vdep0_node)))
    for node in range(1, m + 1):
        v_in = v_arr_sol[node - 1] - np.asarray(spec.v_planets[node], dtype=np.float64)
        vinf_nodes.append(float(np.linalg.norm(v_in)))

    total_dv = _obj(x_sol)
    dv_per_leg = tuple(float(np.linalg.norm(dv)) for dv in dvs_sol)
    feasible = bool(max_defect < feas_tol and max_cvinf < flyby_tol)
    flyby_feasible = bool(max_cvinf < flyby_tol and (min_calt >= -flyby_tol * AU_KM))

    return FbsFlybyOptimizeResult(
        total_dv_kms=total_dv,
        dv_per_leg_kms=dv_per_leg,
        dvs=tuple(dvs_sol),
        v_dep_per_node=tuple(v_dep_sol),
        v_arr_per_node=tuple(v_arr_sol),
        vinf_per_node_kms=tuple(vinf_nodes),
        max_defect=max_defect,
        max_cvinf=max_cvinf,
        min_calt_km=min_calt,
        feasible=feasible,
        flyby_feasible=flyby_feasible,
        success=bool(res.success),
        nfev=int(res.get("nfev", 0)),
        njev=int(res.get("njev", 0)),
        nit=int(res.get("nit", 0)),
        constr_nfev=int(counts["con"]),
        constr_njev=int(counts["cjac"]),
        wall_s=float(wall),
        used_analytic_jac=use_analytic_jac,
    )


__all__ = [
    "FbsFlybyOptimizeResult",
    "FlybyChainSpec",
    "optimize_chain_fbs_flyby",
]
