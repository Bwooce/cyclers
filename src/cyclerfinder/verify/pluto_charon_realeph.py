"""Pluto-Charon real-ephemeris differential corrector (#511).

This is the PROPER real-eph lever for the catalogue row
``ross-rt-pc-cycler-32-2026`` (Pluto-Charon (3,2) CR3BP cycler, V2-ballistic
since #505), as scoped by #506: "differentially correcting the CR3BP IC in
the real-eph model to find the real-eph analog of the periodic orbit" -- NOT
a naive propagation test (#506 already rejected that: it only measures a
~756 m model-mismatch oscillation, not stability).

Two pieces:

1. :func:`charon_osculating_elements` reads Charon's REAL barycentric state
   relative to Pluto from the furnished ``plu060.bsp`` SPICE kernel (#510)
   and computes classical two-body osculating orbital elements. This
   replaces #506's back-of-envelope ``e < 5e-5`` (Brozovic et al. 2015 MEAN
   eccentricity) with the actual OSCULATING eccentricity the real ephemeris
   exhibits at a chosen epoch (short-period terms from Nix/Hydra/Kerberos/
   Styx + solar perturbation included) -- ~2.3e-4, about 5x the mean value.

2. :func:`differential_correct_pc32_to_eccentricity` re-targets the CR3BP
   (3,2) periodic orbit into the Elliptic Restricted 3-Body Problem (ER3BP,
   ``core/er3bp.py`` / ``genome/er3bp_periodic.py``) at that real eccentricity
   via the existing, independently-validated (#293) e-continuation corrector
   (``genome/er3bp_continuation.py``). The ER3BP pulsating-frame model IS the
   real (non-circular) Pluto-Charon two-body dynamics -- Nix/Hydra/Kerberos/
   Styx contribute <2e-6 of the system mass fraction (#506 Gate-b analysis),
   so a Keplerian ellipse for Charon's motion (which the real SPICE state
   confirms to ~2e-4 eccentricity) is the physically complete "real dynamics"
   model at this row's precision.

Both pieces are gated on ``verify.spice_kernels.ensure_pluto_kernel()`` --
the plu060.bsp kernel is local-only (129 MB) and never committed; callers
(scripts, tests) must catch the ``RuntimeError`` and skip cleanly, mirroring
the ``ensure_jup365_kernel`` convention.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.er3bp import ER3BPSystem
from cyclerfinder.genome.er3bp_continuation import ContinuationError, continue_er3bp_family_in_e
from cyclerfinder.genome.er3bp_periodic import ER3BPPeriodicOrbit, correct_er3bp_periodic

#: Pluto+Charon system GM (satellites.py PRIMARIES["Pluto"] / core/cr3bp.py
#: convention): the JPL DE440/PLU060 SYSTEM GM, already includes Charon.
GM_PLUTO_CHARON_SYSTEM_KM3_S2 = 975.5

_FURNISHED: set[str] = set()


def _furnish(kernel_path: str) -> None:
    import spiceypy

    if kernel_path not in _FURNISHED:
        spiceypy.furnsh(kernel_path)
        _FURNISHED.add(kernel_path)


@dataclass(frozen=True)
class CharonOsculatingElements:
    """Charon's real osculating orbit about Pluto at one epoch (SPICE-derived)."""

    epoch_iso: str
    a_km: float
    eccentricity: float
    period_days: float
    r_km: float


def charon_osculating_elements(epoch_iso: str, kernel_path: str) -> CharonOsculatingElements:
    """Charon's real (SPICE) osculating orbital elements relative to Pluto.

    Reads the state of Charon (NAIF 901) relative to Pluto (NAIF 999) from the
    furnished ``plu060.bsp`` kernel in the J2000 frame at ``epoch_iso`` (TDB
    calendar string), then computes classical two-body osculating elements
    (vis-viva energy + Laplace-Runge-Lenz eccentricity vector) using
    ``GM_PLUTO_CHARON_SYSTEM_KM3_S2``. This is the REAL (not circular-model)
    Charon orbit -- the input the ER3BP differential corrector needs.
    """
    import spiceypy
    from astropy.time import Time

    _furnish(kernel_path)

    et = float((Time(epoch_iso, scale="tdb") - Time(2451545.0, format="jd", scale="tdb")).sec)
    state, _lt = spiceypy.spkezr("901", et, "J2000", "NONE", "999")
    r: NDArray[np.float64] = np.asarray(state[:3], dtype=np.float64)
    v: NDArray[np.float64] = np.asarray(state[3:], dtype=np.float64)
    mu = GM_PLUTO_CHARON_SYSTEM_KM3_S2

    rn = float(np.linalg.norm(r))
    vn = float(np.linalg.norm(v))
    energy = vn * vn / 2.0 - mu / rn
    a = -mu / (2.0 * energy)
    h = np.cross(r, v)
    e_vec = np.cross(v, h) / mu - r / rn
    e = float(np.linalg.norm(e_vec))
    period_days = 2.0 * math.pi * math.sqrt(a**3 / mu) / 86400.0

    return CharonOsculatingElements(
        epoch_iso=epoch_iso, a_km=float(a), eccentricity=e, period_days=period_days, r_km=rn
    )


@dataclass(frozen=True)
class PC32RealEphResult:
    """Differential-correction verdict for the PC (3,2) row at real eccentricity."""

    e_real: float
    period_ratio_sc_to_charon: float
    """T_(3,2) / T_Charon in the CR3BP e=0 seed -- commensurability check."""
    seed_corrector_residual: float
    seed_independent_residual: float
    target_corrector_residual: float
    target_independent_residual: float
    """L2 norm of X(period_f) - X(0) from an INDEPENDENT Radau re-propagation
    at the target eccentricity -- the #441 period_f-trap gate. A converged
    ``target_corrector_residual`` with a large ``target_independent_residual``
    is the false-positive-closure signature: the symmetric half-period
    condition is satisfied but the orbit does not actually close."""
    converged: bool
    """True iff BOTH residuals are below tolerance at the target eccentricity."""
    orbit: ER3BPPeriodicOrbit | None


def differential_correct_pc32_to_eccentricity(
    e_real: float,
    *,
    x0_seed: float,
    ydot0_seed: float,
    period_seed: float,
    mu: float,
    n_steps: int = 20,
    tol: float = 1e-9,
    independent_tol: float = 1e-8,
) -> PC32RealEphResult:
    """Differentially correct the PC (3,2) CR3BP orbit into ER3BP at ``e_real``.

    Bridges the CR3BP seed into the ER3BP pulsating-frame corrector at e=0
    (exact -- #441 Sec. 1), then walks the secant e-continuator
    (``continue_er3bp_family_in_e``) from e=0 to ``e_real`` in ``n_steps``
    steps. Reports BOTH the corrector residual (which only tests the
    symmetric half-period crossing condition) and the independent Radau
    full-orbit-closure residual (the #441 period_f-trap gate) at the target
    eccentricity, so a false-positive "converged" result is visible, not
    hidden.
    """
    seed_state = np.array([x0_seed, 0.0, 0.0, 0.0, ydot0_seed, 0.0], dtype=np.float64)
    half_period = period_seed / 2.0
    period_ratio = period_seed / (2.0 * math.pi)

    sys0 = ER3BPSystem(mu=mu, e=0.0, primary_name="Pluto", secondary_name="Charon")
    seed_orbit = correct_er3bp_periodic(
        sys0, seed_state, half_period, is_half_period_residual=True, tol=tol
    )

    try:
        history = continue_er3bp_family_in_e(
            sys0,
            seed_state,
            half_period,
            e_real,
            n_steps,
            is_half_period_residual=True,
            tol=tol,
        )
        _last = history[-1]
        target_corrector_residual = _last.corrector_residual
        target_independent_residual = _last.independent_residual
        target_orbit: ER3BPPeriodicOrbit | None = _last
    except ContinuationError:
        target_orbit = None
        target_corrector_residual = float("nan")
        target_independent_residual = float("nan")

    converged = (
        target_orbit is not None
        and target_corrector_residual < tol
        and target_independent_residual < independent_tol
    )

    return PC32RealEphResult(
        e_real=e_real,
        period_ratio_sc_to_charon=period_ratio,
        seed_corrector_residual=seed_orbit.corrector_residual,
        seed_independent_residual=seed_orbit.independent_residual,
        target_corrector_residual=target_corrector_residual,
        target_independent_residual=target_independent_residual,
        converged=converged,
        orbit=target_orbit,
    )
