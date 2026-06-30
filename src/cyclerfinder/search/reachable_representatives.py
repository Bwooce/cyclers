"""Recover the Braik-Ross 2026 common-energy representative orbits at C_J=3.1294.

The paper (Table 2) lists thirteen planar Earth-Moon periodic orbits, one per
family, all re-selected at a single common Jacobi constant ``C_J = 3.1294``, but
publishes ONLY period-in-days and the Floquet rate sigma -- no state vectors. To
score the accessibility network we must first *recover the member* at that
energy, then confirm it against the sourced period before trusting it.

Two recovery routes (mining note Q3 / proposed task 2):

* Lyapunov (LL1, LL2), distant prograde (DPO), and the resonant families are
  standard JPL 3-Body Periodic-Orbit families. We pull the member nearest
  ``C_J = 3.1294`` from the JPL oracle
  (:func:`cyclerfinder.verify.jpl_periodic_orbits.query`) as a seed, then
  re-correct it under OUR mass ratio with the fixed-Jacobi symmetric corrector so
  the recovered orbit sits exactly on ``C_J = 3.1294`` (the JPL mu differs from
  ours by ~1e-7; see ``verify.jpl_periodic_orbits`` CONVENTION RECONCILIATION).
* The four cyclers (C11a, C11b, C21, C32) are NOT in the JPL DB. We recover them
  with our own fixed-Jacobi symmetric corrector from the Ross & Roberts-Tsoukkas
  2025 (AAS 25-621) family seed regions (x0 region + half-crossing index), driven
  to the sourced ``C_J = 3.1294`` period.

SOURCED-CONFIRMATION DISCIPLINE: every recovered member is only *trusted* once its
period matches the Braik-Ross Table-2 sourced period (in days) within tolerance.
A member whose period does not match is reported as unconfirmed and must NOT enter
the scored network -- selecting members by our own criteria instead of a
published anchor would make the validation gate circular.

This module performs a network call (JPL) for the non-cycler families; the
cyclers are recovered offline. Pure-corrector helpers are network-free.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import cyclerfinder.core.cr3bp as cr3bp
import cyclerfinder.search.cr3bp_periodic as cp
from cyclerfinder.verify.jpl_periodic_orbits import query

#: Braik-Ross common Jacobi constant (Table 2 header / Sec. 2).
C_J_BRAIK_ROSS = 3.1294

#: Exact Jacobi constant for the (2,1) C21 cycler family at the Braik-Ross
#: common energy. Braik-Ross 2026 prints "C_J = 3.1294" throughout but the (2,1)
#: family has full Jacobi extent ΔC ~ 4e-12 — essentially a single point in C at
#: this value. The literal 3.1294 sits ~1e-5 above the family's C-max, so C21 is
#: only recoverable at this unrounded Jacobi (see #249 / task #495).
#:
#: PRIMARY SOURCE (task #495): Braik & Ross 2026 repo (MIT license),
#: ``src/cr3bp_family_ic.m``, independently computed via correct_po_to_cj313_v2.m.
#: Value: 3.129389531054557 (15 sig figs), stored verbatim in
#: ``data/golden/braik_ross_2026_em_family_ics.yaml``.
#:
#: CROSS-SOURCE: Ross & Roberts-Tsoukkas 2025 (AAS 25-621, Table 4) gives
#: 3.129389531088256 — differs at the 12th decimal place (~3e-11); both values
#: produce a confirmed period recovery. The Braik-repo value is used as the
#: primary because it is independently computed against the same ICs used in
#: the accessibility-network analysis.
#:
#: Standing rule ``feedback_published_rounded_values_are_display``: a printed
#: C/T/V_inf is a display value, not a literal — when a search returns no
#: topology-correct matches, suspect the printed value is rounded.
C_J_C21 = 3.129389531054557  # Braik & Ross 2026 repo (MIT) — cr3bp_family_ic.m
C_J_C21_AAS = 3.129389531088256  # cross-source: Ross-RT 2025 AAS-25-621 Table 4

#: Ross & Roberts-Tsoukkas 2025 Earth-Moon mass ratio (AAS 25-621, p. 3). Used
#: for the cycler recovery so the corrector matches the family-defining paper.
ROSS_MU = 1.2150584270572e-2

#: 1 nondimensional time unit in days (T_EM = 27.321661 d, TU = T_EM / 2pi).
TU_DAYS = 27.321661 / (2.0 * math.pi)


def braik_ross_system() -> cr3bp.CR3BPSystem:
    """Earth-Moon CR3BP system in the Braik-Ross / Ross-RT nondimensional scales.

    Uses the AAS 25-621 mass ratio (identical to Braik-Ross Table 1 to all
    printed digits) and the standard a_M = 384400 km length / T_EM/2pi time
    scales, so periods convert to days via :data:`TU_DAYS`.
    """
    return cr3bp.CR3BPSystem(
        mu=ROSS_MU,
        primary="Earth",
        secondary="Moon",
        l_km=384400.0,
        t_s=TU_DAYS * 86400.0,
    )


@dataclass(frozen=True)
class Representative:
    """One recovered common-energy representative orbit.

    ``state0`` is the planar IC ``(x0, 0, 0, 0, ydot0, 0)``; ``period`` is
    nondimensional (multiply by :data:`TU_DAYS` for days). ``sourced_period_days``
    is the Braik-Ross Table-2 value; ``confirmed`` is True iff the recovered
    period matches it within ``tol_days``.
    """

    label: str
    state0: NDArray[np.float64]
    period: float
    jacobi: float
    sourced_period_days: float
    period_days: float
    confirmed: bool
    converged: bool


# Braik-Ross Table 2 sourced periods (days) at C_J = 3.1294 for the families we
# can source-confirm. (Resonant R31/R52 and the U-branches carry no period in the
# data available to us, so they are intentionally excluded from the gate set --
# including them would assert un-cross-checkable members; see module docstring.)
SOURCED_PERIODS_DAYS: dict[str, float] = {
    "LL1": 12.811,
    "LL2": 15.117,
    "DPO": 11.184,
    "R21-S": 26.500,
    "C11a": 42.140,
    "C11b": 55.995,
    "C21": 84.533,
    "C32": 78.613,
}

# Ross & Roberts-Tsoukkas 2025 cycler family seed regions. Each entry carries:
#   x0           — perpendicular-x-axis IC seed (x at the half-period crossing).
#   ydot0_sign   — sign of ydot0 (the symmetric corrector's velocity branch).
#   half_crossings — 1-based index of the perpendicular x-axis crossing taken as
#                    the half-period (a (k1,k2) cycler crosses the x-axis many
#                    times per period; this selects the right branch).
#   c_override   — per-family Jacobi for recovery, or ``None`` to use
#                    :data:`C_J_BRAIK_ROSS`. C21's (2,1) family has ΔC ~ 4e-12
#                    (essentially a single point in C); the literal 3.1294 sits
#                    ~1e-5 above the family's C-max so we must use the unrounded
#                    :data:`C_J_C21` from the Braik-Ross 2026 repo (MIT) to
#                    land in-family.
#
# Source for the x0 region / half-crossing index (task #495, post-#249):
# Braik & Ross 2026 repo (MIT), src/cr3bp_family_ic.m, independently computed
# via correct_po_to_cj313_v2.m; stored in data/golden/braik_ross_2026_em_family_ics.yaml.
# C11a / C21 exact ICs adopted from the MIT-licensed repo (task #495).
#
# RECOVERY STATUS at the per-family Jacobi (post-#249, 4/4, sourced #495):
#   * C11a (1,1): 42.1405 d  — corrector converges AT Braik IC (x0 match exact).
#   * C11b (1,1): 55.9590 d  — recovered at C_J_BRAIK_ROSS (0.06%).
#   * C21  (2,1): 84.5331 d  — recovered at C_J_C21=3.129389531054557 (0.4 ppm).
#   * C32  (3,2): 78.6126 d  — recovered at C_J_BRAIK_ROSS (0.0005%).
_CYCLER_SEEDS: dict[str, dict[str, float | int | None]] = {
    "C11a": {
        # Exact x0 from Braik & Ross 2026 repo cr3bp_family_ic.m (#495)
        "x0": -0.8116406668238326,
        "ydot0_sign": -1.0,
        "half_crossings": 3,
        "c_override": None,
    },
    "C11b": {
        "x0": -0.7684981,
        "ydot0_sign": -1.0,
        "half_crossings": 6,
        "c_override": None,
    },
    "C21": {
        # Exact x0 from Braik & Ross 2026 repo cr3bp_family_ic.m (#495)
        "x0": 7.237366530581342e-01,
        "ydot0_sign": +1.0,
        "half_crossings": 4,
        "c_override": C_J_C21,  # 3.129389531054557 — Braik-repo sourced
    },
    "C32": {
        "x0": -0.2752115,
        "ydot0_sign": -1.0,
        "half_crossings": 6,
        "c_override": None,
    },
}


def _jpl_seed_near_cj(
    family: str,
    *,
    libr: int | None = None,
    branch: str | None = None,
    cj: float = C_J_BRAIK_ROSS,
    stable: bool | None = None,
) -> NDArray[np.float64]:
    """JPL member nearest ``cj`` (optionally filtered to the stable branch).

    ``stable=True`` keeps members with stability index near 1 (|stab| < ~1+eps,
    bounded/stable in JPL's reduced index convention); ``stable=False`` keeps the
    strongly unstable members; ``None`` takes the global nearest-in-C member.
    """
    _constants, orbits = query("earth-moon", family, libr=libr, branch=branch)
    js = np.array([o.jacobi for o in orbits])
    stab = np.array([abs(o.stability) for o in orbits])
    mask = np.ones(len(orbits), dtype=bool)
    if stable is True:
        mask = stab <= 1.0 + 1e-6
    elif stable is False:
        mask = stab > 1.0 + 1e-6
    if not mask.any():
        mask = np.ones(len(orbits), dtype=bool)
    cand = np.where(mask)[0]
    i = int(cand[np.argmin(np.abs(js[cand] - cj))])
    return np.asarray(orbits[i].state0, dtype=np.float64)


def recover_from_seed(
    system: cr3bp.CR3BPSystem,
    label: str,
    x0_seed: float,
    period_guess_days: float,
    *,
    ydot0_sign: float,
    half_crossings: int,
    tol_days: float = 0.5,
    corrector_tol: float = 1e-10,
    jacobi: float | None = None,
) -> Representative:
    """Correct a symmetric member to a fixed Jacobi and confirm its period.

    Wraps :func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`
    (Jacobi held at ``jacobi`` if given, else :data:`C_J_BRAIK_ROSS`) and compares
    the recovered period to the Braik-Ross sourced period. ``period_guess_days``
    seeds the corrector's period; the sourced period is the confirmation target.

    Per-family Jacobi (#262 / post-#249): the (2,1) C21 family has full Jacobi
    extent ~ 4e-12 and only exists at the unrounded :data:`C_J_C21`, NOT at the
    literal :data:`C_J_BRAIK_ROSS` printed in the Braik-Ross paper -- pass
    ``jacobi=C_J_C21`` for that recovery. C11a/C11b/C32 use the literal value.
    """
    sourced = SOURCED_PERIODS_DAYS[label]
    period_guess = period_guess_days / TU_DAYS
    c_target = C_J_BRAIK_ROSS if jacobi is None else float(jacobi)
    orbit = cp.correct_symmetric_fixed_jacobi(
        system,
        x0_seed,
        c_target,
        period_guess,
        ydot0_sign=ydot0_sign,
        half_crossings=half_crossings,
        tol=corrector_tol,
    )
    period_days = orbit.period * TU_DAYS
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    confirmed = orbit.converged and abs(period_days - sourced) <= tol_days
    return Representative(
        label=label,
        state0=state0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        sourced_period_days=sourced,
        period_days=period_days,
        confirmed=confirmed,
        converged=orbit.converged,
    )


def recover_all_cyclers_braik_ross(
    system: cr3bp.CR3BPSystem,
    *,
    tol_days: float = 0.5,
    corrector_tol: float = 1e-10,
) -> list[Representative]:
    """Recover all four Braik-Ross cycler members (C11a, C11b, C21, C32) (#262).

    Uses the seeds + per-family Jacobi pinned in :data:`_CYCLER_SEEDS` to recover
    each member with the symmetric perpendicular-x-axis-crossing corrector and
    confirms the period against the Braik-Ross Table-2 sourced value. C21 uses
    its unrounded Jacobi :data:`C_J_C21`; the other three use :data:`C_J_BRAIK_ROSS`.

    Returns one :class:`Representative` per cycler in canonical order
    (C11a, C11b, C21, C32) regardless of confirmation status — the caller filters
    on ``.confirmed`` before scoring.
    """
    out: list[Representative] = []
    for label in ("C11a", "C11b", "C21", "C32"):
        seed = _CYCLER_SEEDS[label]
        x0_raw = seed["x0"]
        sign_raw = seed["ydot0_sign"]
        hc_raw = seed["half_crossings"]
        c_override = seed["c_override"]
        assert x0_raw is not None and sign_raw is not None and hc_raw is not None
        out.append(
            recover_from_seed(
                system,
                label,
                float(x0_raw),
                SOURCED_PERIODS_DAYS[label],
                ydot0_sign=float(sign_raw),
                half_crossings=int(hc_raw),
                tol_days=tol_days,
                corrector_tol=corrector_tol,
                jacobi=None if c_override is None else float(c_override),
            )
        )
    return out


def recover_jpl_family(
    system: cr3bp.CR3BPSystem,
    label: str,
    family: str,
    *,
    libr: int | None = None,
    branch: str | None = None,
    half_crossings: int = 1,
    stable: bool | None = None,
    tol_days: float = 0.5,
) -> Representative:
    """Recover a JPL-DB family member at ``C_J=3.1294`` (seed from JPL, re-correct).

    Pulls the member nearest the common energy from the JPL oracle, then
    re-corrects it under our mass ratio with the fixed-Jacobi symmetric corrector
    so it sits exactly on ``C_J = 3.1294``; confirms against the sourced period.

    NETWORK ROUTE: requires the JPL oracle. For a network-independent recovery use
    :func:`recover_free_period` with an offline seed from :data:`OFFLINE_SEEDS`.
    """
    seed = _jpl_seed_near_cj(family, libr=libr, branch=branch, stable=stable)
    x0_seed = float(seed[0])
    ydot0_sign = math.copysign(1.0, float(seed[4]))
    return recover_from_seed(
        system,
        label,
        x0_seed,
        SOURCED_PERIODS_DAYS[label],
        ydot0_sign=ydot0_sign,
        half_crossings=half_crossings,
        tol_days=tol_days,
    )


# ---------------------------------------------------------------------------
# Network-INDEPENDENT (offline) member recovery (#247).
#
# The JPL oracle is the intended seed source for the Lyapunov / resonant / DPO
# families, but it requires a live network call. The offline route below recovers
# those members from analytic seeds (collinear-point linear Lyapunov amplitude,
# resonant x0 region) using the free-(x0, t_half) perpendicular-crossing corrector
# (:func:`correct_symmetric_free_period`), which -- unlike the 1-DOF fixed-crossing
# corrector -- frees the half-period time and so recovers the member of a *target*
# period region rather than collapsing onto the nearest fixed-crossing branch.
#
# SOURCED-CONFIRMATION DISCIPLINE is unchanged: the seed (x0 region + velocity
# sign + target half-period) is the only family input; the recovered period is a
# PREDICTION confirmed against the Braik-Ross Table-2 sourced period before the
# member is admitted. Where a member does not confirm it is reported unconfirmed
# and excluded rather than faked.
# ---------------------------------------------------------------------------

#: Full Braik-Ross Table 2 (arXiv:2605.31543, p. 11) sourced periods (days) AND
#: Floquet instability rates sigma (TU^-1) at C_J = 3.1294, for ALL thirteen
#: representatives. sigma = 0 for the stable resonants (|lambda_max| = 1).
SOURCED_TABLE2: dict[str, tuple[float, float]] = {
    "LL1": (12.811, 2.4884),
    "LL2": (15.117, 1.9797),
    "C11a": (42.140, 1.0482),
    "C11b": (55.995, 0.9255),
    "C21": (84.533, 0.1358),
    "C32": (78.613, 0.6886),
    "R21-S": (26.500, 0.0),
    "R21-U": (31.039, 0.8397),
    "R31-S": (27.252, 0.0),
    "R31-U": (28.066, 0.40124),
    "R52-S": (54.802, 0.0),
    "R52-U": (56.436, 0.36547),
    "DPO": (11.184, 1.5886),
}

#: Offline recovery seeds: ``label -> (x0_seed, ydot0_sign, jacobi)``. The target
#: half-period is taken from the sourced period; ``jacobi`` is the common energy
#: 3.1294 except for C21, which is the Ross & Roberts-Tsoukkas 2025 (2,1) STABLE
#: sourced member at its own C = 3.129389531 (Braik-Ross's 3.1294 is the rounded
#: value of this same member). The seed x0 region is sourced from the family
#: geometry (collinear-point Lyapunov amplitude, resonant x0 band, AAS-25-621
#: cycler seed); the period AND the Floquet sigma are PREDICTIONS confirmed against
#: Table 2. Only families that CONFIRM (period within tolerance AND sigma within
#: tolerance) are admitted by :func:`recover_offline_set` -- no faked members.
#:
#: RECOVERED (period + sigma both confirmed offline, #247): the nine below.
#: NOT RECOVERED at this off-stable common energy with the available single-/
#: free-period shooting correctors (excluded, not faked): the three unstable
#: cyclers C11a (sigma 1.05), C11b (sigma 0.93), C32 (sigma 0.69) -- which all
#: collapse onto nearby spurious lower-sigma orbits -- and the 5:2 unstable
#: resonant R52-U (sigma 0.37). A robust Jacobi-constrained multiple-shooting
#: corrector (or the JPL oracle, unavailable in this environment) would be needed.
#: Alias: Braik-repo sourced CJ for C21 (task #495); equals C_J_C21.
ROSS_C21_JACOBI = C_J_C21  # 3.129389531054557 — Braik & Ross 2026 repo (MIT)
#: Cross-source value from AAS-25-621 Table 3 for reference only.
ROSS_C21_JACOBI_AAS = C_J_C21_AAS  # 3.129389531088256 — Ross-RT 2025 AAS-25-621
OFFLINE_SEEDS: dict[str, tuple[float, float, float]] = {
    # Exact x0 seeds from Braik & Ross 2026 repo cr3bp_family_ic.m (#495);
    # others use approximate x0 from literature / manual estimation.
    "LL1": (0.8115256290557147, 1.0, C_J_BRAIK_ROSS),  # Braik exact x0
    "LL2": (1.100554841329441, 1.0, C_J_BRAIK_ROSS),  # Braik exact x0
    "DPO": (1.060820806800189, 1.0, C_J_BRAIK_ROSS),  # Braik exact x0
    "R21-S": (0.4485378415692721, 1.0, C_J_BRAIK_ROSS),  # Braik exact x0
    "R21-U": (-0.812, -1.0, C_J_BRAIK_ROSS),  # 2:1 unstable resonant (no exact x0 adopted)
    "R31-S": (0.3568, -1.0, C_J_BRAIK_ROSS),  # 3:1 stable resonant
    "R31-U": (0.138, 1.0, C_J_BRAIK_ROSS),  # 3:1 unstable resonant
    "R52-S": (0.2278881086717652, 1.0, C_J_BRAIK_ROSS),  # Braik exact x0
    "C21": (7.237366530581342e-01, 1.0, ROSS_C21_JACOBI),  # Braik exact x0 + CJ (#495)
}


def lagrange_collinear_x(mu: float, point: str) -> float:
    """x-coordinate of the collinear libration point ``point`` (``L1`` or ``L2``).

    Root of ``dUbar/dx = 0`` on the x-axis. ``L1`` lies between the primaries,
    ``L2`` beyond the secondary. Offline (no network); used to place the linear
    Lyapunov seed.
    """
    from scipy.optimize import brentq

    def f(x: float) -> float:
        return cp._ubar_grad_x_at_axis(x, mu)

    if point == "L1":
        return float(brentq(f, -mu + 0.05, 1.0 - mu - 1e-3))
    if point == "L2":
        return float(brentq(f, 1.0 - mu + 1e-3, 1.7))
    raise ValueError(f"point must be 'L1' or 'L2', got {point!r}")


def correct_symmetric_free_period(
    system: cr3bp.CR3BPSystem,
    x0_guess: float,
    jacobi: float,
    t_half_guess: float,
    *,
    ydot0_sign: float = 1.0,
    tol: float = 1e-10,
    max_iter: int = 50,
    rtol: float = 1e-12,
    atol: float = 1e-12,
    x0_max_step: float = 0.1,
    t_half_step_frac: float = 0.2,
    x0_bounds: tuple[float, float] = (-2.0, 2.0),
) -> cp.SymmetricOrbit:
    """Free-(x0, t_half) fixed-Jacobi perpendicular-crossing corrector (#247).

    A more general member-recovery corrector than
    :func:`cyclerfinder.search.cr3bp_periodic.correct_symmetric_fixed_jacobi`:
    instead of snapping to a fixed *event crossing index* (which makes the period
    an output of whichever branch the Newton iterate slides onto), this corrector
    treats the half-period time ``t_half`` as a *free variable* and drives the two
    perpendicular-crossing residuals to zero simultaneously,

        F1: y(t_half)    = 0,
        F2: xdot(t_half) = 0,

    with ``ydot0 = ydot0_from_jacobi(x0)`` holding the Jacobi constant fixed (Ross
    Eq. 9). The 2x2 Newton system uses the half-period STM (columns 0 and 4, the
    latter coupled through ``dydot0/dx0 = -(dUbar/dx)/ydot0``) and the EOM time
    derivative at ``t_half``. Because ``t_half`` is continuous, seeding it near a
    *target* half-period recovers the member of that period region rather than the
    nearest fixed-index attractor -- this separates distinct same-energy members
    (e.g. a resonant stable/unstable pair) that the 1-DOF corrector collapses.

    Returns a :class:`~cyclerfinder.search.cr3bp_periodic.SymmetricOrbit` with
    ``period = 2 * t_half``; ``converged`` iff ``sqrt(F1^2 + F2^2) < tol``.
    """
    mu = system.mu
    x0 = float(x0_guess)
    t_half = abs(float(t_half_guess))
    lo, hi = x0_bounds
    res = float("inf")
    n_iter = 0
    for n_iter in range(1, max_iter + 1):  # noqa: B007 -- returned as iteration count
        ydot0 = cp.ydot0_from_jacobi(x0, jacobi, mu, sign=ydot0_sign)
        state0 = np.array([x0, 0.0, 0.0, 0.0, ydot0, 0.0])
        arc = cr3bp.propagate(system, state0, t_half, with_stm=True, rtol=rtol, atol=atol)
        yf = arc.state_f
        assert arc.stm is not None
        stm = arc.stm
        f1 = float(yf[1])  # y at t_half
        f2 = float(yf[3])  # xdot at t_half
        res = math.hypot(f1, f2)
        if res < tol:
            break
        fdot = cr3bp.cr3bp_eom(t_half, yf, mu)
        dydot0_dx0 = -cp._ubar_grad_x_at_axis(x0, mu) / ydot0
        ds_dx0_y = float(stm[1, 0]) + float(stm[1, 4]) * dydot0_dx0
        ds_dx0_xdot = float(stm[3, 0]) + float(stm[3, 4]) * dydot0_dx0
        jac = np.array(
            [
                [ds_dx0_y, float(fdot[1])],  # dy/dx0,    dy/dt_half
                [ds_dx0_xdot, float(fdot[3])],  # dxdot/dx0, dxdot/dt_half
            ]
        )
        try:
            step = np.linalg.solve(jac, np.array([-f1, -f2]))
        except np.linalg.LinAlgError:
            break
        dx0 = float(step[0])
        dth = float(step[1])
        if abs(dx0) > x0_max_step:
            dx0 = math.copysign(x0_max_step, dx0)
        if abs(dth) > t_half_step_frac * t_half:
            dth = math.copysign(t_half_step_frac * t_half, dth)
        x0 = min(max(x0 + dx0, lo), hi)
        t_half = t_half + dth
        if t_half <= 0.0:
            t_half = 0.5 * abs(float(t_half_guess))
    ydot0_final = cp.ydot0_from_jacobi(x0, jacobi, mu, sign=ydot0_sign)
    period = 2.0 * t_half
    converged = res < tol
    return cp.SymmetricOrbit(
        x0=x0,
        ydot0=ydot0_final,
        jacobi=cr3bp.jacobi_constant(np.array([x0, 0.0, 0.0, 0.0, ydot0_final, 0.0]), mu),
        t_half=t_half,
        period=period,
        converged=converged,
        crossing_residual=res,
        n_iter=n_iter,
    )


def recover_free_period(
    system: cr3bp.CR3BPSystem,
    label: str,
    x0_seed: float,
    *,
    ydot0_sign: float,
    jacobi: float = C_J_BRAIK_ROSS,
    tol_days: float = 0.5,
    sigma_tol: float = 0.15,
    corrector_tol: float = 1e-10,
) -> Representative:
    """Recover a member at ``jacobi`` offline via the free-period corrector.

    Seeds the corrector at the sourced half-period and the supplied ``x0`` region;
    the recovered period is confirmed against the Braik-Ross sourced period AND the
    Floquet rate ``sigma = ln(lambda_max)/T`` is checked against the sourced sigma
    (stable resonants must give sigma ~ 0). ``confirmed`` is True only if BOTH the
    period (within ``tol_days``) and sigma (within ``sigma_tol``) match -- the
    sigma check is what rejects spurious nearby orbits that merely share a period.
    """
    sourced_days, sourced_sigma = SOURCED_TABLE2[label]
    t_half_guess = 0.5 * sourced_days / TU_DAYS
    orbit = correct_symmetric_free_period(
        system,
        x0_seed,
        jacobi,
        t_half_guess,
        ydot0_sign=ydot0_sign,
        tol=corrector_tol,
    )
    period_days = orbit.period * TU_DAYS
    _nu, lam = cp.barden_stability(system, orbit)
    sigma = math.log(abs(lam)) / orbit.period if abs(lam) > 1.0 + 1e-9 else 0.0
    state0 = np.array([orbit.x0, 0.0, 0.0, 0.0, orbit.ydot0, 0.0])
    confirmed = (
        orbit.converged
        and abs(period_days - sourced_days) <= tol_days
        and abs(sigma - sourced_sigma) <= sigma_tol
    )
    return Representative(
        label=label,
        state0=state0,
        period=orbit.period,
        jacobi=orbit.jacobi,
        sourced_period_days=sourced_days,
        period_days=period_days,
        confirmed=confirmed,
        converged=orbit.converged,
    )


def recover_offline_set(system: cr3bp.CR3BPSystem) -> list[Representative]:
    """Recover every member that has an :data:`OFFLINE_SEEDS` entry, offline.

    Returns one :class:`Representative` per seeded family (confirmed or not). The
    caller filters on ``.confirmed`` before scoring -- no faked members.
    """
    out: list[Representative] = []
    for label, (x0, sign, cj) in OFFLINE_SEEDS.items():
        out.append(recover_free_period(system, label, x0, ydot0_sign=sign, jacobi=cj))
    return out
