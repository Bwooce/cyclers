"""#391 — Amplitude-vs-Hill-fraction pre-screen for CR3BP cycler candidates.

A cheap, pre-V4 structural predictor of real-ephemeris survivability for an
Earth-Moon CR3BP periodic orbit.

Motivation (#389 V4 HALT)
-------------------------
The #347 Floquet Phase 2 framework discovered ``branch_C32_b0``, a planar (3,3)
Earth-Moon CR3BP periodic orbit that is *spectrally* almost perfectly stable
(max |Floquet| = 1.000000000000617). It cleared V1/V2/V3 — the idealized
CR3BP / 2-body gates — by 6-9 orders of magnitude, then **failed V4 by 4-5
orders** (real DE440 solar tide drives it to a ~10⁹ km escape, 0/100 launch
epochs 2000-2099). See ``docs/notes/2026-06-18-389-branch_C32_b0-admission.md``.

The structural cause, isolated in ``data/branch_c32_b0_v4_verdict.jsonl``
(``structural_diagnostic`` row + Sun-only control), is amplitude: the orbit's
farthest excursion reaches **0.77 of the Earth-Sun Hill radius**, where the
solar tidal acceleration is ~30% of Earth's gravity. The CR3BP — and V1-V3,
which live in or near it — ignore the Sun entirely, so spectral stability in the
autonomous CR3BP says *nothing* about survival once the real solar tide is
switched on at near-Hill amplitude.

The lesson: an orbit whose amplitude is a large fraction of the Earth-Sun Hill
radius is V4-doomed regardless of its CR3BP Floquet character. A one-propagation
amplitude/Hill-fraction screen flags such families *before* anyone spends a
multi-tier V0-V5 gauntlet on them.

Threshold rationale
-------------------
The solar tidal acceleration on a body at Earth-distance ``r`` scales as
``a_tide ≈ 2 G M_sun r / a_ES³`` (linearized, leading order), i.e. **linearly
in amplitude**, and equals Earth's own gravity at the Hill radius by definition.
So the tide-to-Earth-gravity ratio rises roughly as ``(r / r_Hill)³`` near the
Hill radius but is well approximated by the Hill fraction itself as the
order-of-magnitude knob.

* ``branch_C32_b0`` at Hill-fraction **0.77** → tide 30% of Earth gravity →
  catastrophic V4 failure (escape).
* We want the gate to flag *well below* 0.77, with margin, because the failure
  is steep and the V4 cost is high.

Classification bands:

* ``PASS`` (< 0.3): solar tide < ~10% of Earth gravity at apoapsis; deep inside
  the Hill sphere. A genuine V4-survival *candidate* (necessary, not sufficient
  — V4 must still be run).
* ``MARGINAL`` (0.3-0.5): non-negligible solar tide; V4 outcome uncertain. Worth
  a gauntlet but flag the risk.
* ``V4_DOOMED`` (> 0.5): solar tide is a large fraction of Earth gravity; expect
  a V4 escape failure like ``branch_C32_b0``. Tag ``cr3bp-only`` up front; do
  not spend a gauntlet without an amplitude-reducing re-scope first.

Sourcing discipline
-------------------
Every physical constant traces to :mod:`cyclerfinder.core.constants` or
:mod:`cyclerfinder.core.satellites`:

* ``MU_SUN_KM3_S2`` — JPL DE440 / IAU 2015 solar GM.
* ``AU_KM`` — IAU 2012 Resolution B2 (exact).
* ``PLANETS["E"].sma_au`` — Earth heliocentric SMA (Standish & Williams).
* ``PRIMARIES["Earth"]`` — G(M_Earth + M_Moon); the correct lumped mass for the
  Earth-Moon barycenter's Hill sphere about the Sun.
* characteristic length from the :class:`~cyclerfinder.core.cr3bp.CR3BPSystem`
  (``l_km`` = the secondary's SMA about the primary, 384400 km for Earth-Moon).

No hardcoded magic numbers.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import solve_ivp

import cyclerfinder.core.cr3bp as cr3bp
from cyclerfinder.core.constants import AU_KM, MU_SUN_KM3_S2, PLANETS
from cyclerfinder.core.satellites import PRIMARIES

#: Hill-fraction band edges (orbit max-amplitude / Earth-Sun Hill radius).
#: See module docstring for the #389 ``branch_C32_b0`` (0.77 → escape) rationale.
PASS_HILL_FRACTION: float = 0.3
DOOMED_HILL_FRACTION: float = 0.5

#: Classification labels.
CLASS_PASS = "PASS"
CLASS_MARGINAL = "MARGINAL"
CLASS_V4_DOOMED = "V4_DOOMED"


def earth_sun_hill_radius_km() -> float:
    """Earth-Sun (Earth-Moon barycenter) Hill radius, km, from sourced constants.

    ``r_H = a_ES * (G(M_E+M_M) / (3 G M_sun))^(1/3)``.

    Uses the lumped Earth+Moon GM (``PRIMARIES["Earth"]``) because the body whose
    Hill sphere matters for an Earth-Moon CR3BP orbit is the Earth-Moon
    barycenter orbiting the Sun, not the Earth alone. The numeric value
    (~1.50e6 km) matches the #389 V4 verdict's ``earth_sun_hill_radius_km`` to
    within the difference between the sourced Earth SMA (1.0000026 au) and the
    exact 1 au the verdict used.
    """
    a_es_km = float(PLANETS["E"].sma_au) * AU_KM
    gm_emb = float(PRIMARIES["Earth"])  # G(M_Earth + M_Moon)
    return float(a_es_km * (gm_emb / (3.0 * MU_SUN_KM3_S2)) ** (1.0 / 3.0))


def _classify(hill_fraction: float) -> str:
    if hill_fraction < PASS_HILL_FRACTION:
        return CLASS_PASS
    if hill_fraction < DOOMED_HILL_FRACTION:
        return CLASS_MARGINAL
    return CLASS_V4_DOOMED


@dataclass(frozen=True)
class HillScreenResult:
    """Outcome of the amplitude-vs-Hill-fraction pre-screen for one orbit.

    Attributes
    ----------
    max_amplitude_km :
        The orbit's farthest excursion from the primary (Earth) over one full
        period, in physical km (mapped via the CR3BP characteristic length).
    earth_sun_hill_radius_km :
        Earth-Sun (Earth-Moon barycenter) Hill radius, km (sourced constants).
    hill_fraction :
        ``max_amplitude_km / earth_sun_hill_radius_km``.
    solar_tide_to_earth_gravity_ratio :
        Linearized solar tidal acceleration at ``max_amplitude_km`` divided by
        Earth's gravity at that distance — the physical quantity the Hill
        fraction is a proxy for. ~0.30 at ``branch_C32_b0``'s 0.77 fraction.
    classification :
        One of ``PASS`` / ``MARGINAL`` / ``V4_DOOMED`` (see module docstring).
    n_samples :
        Number of dense-output samples used to find the max excursion.
    """

    max_amplitude_km: float
    earth_sun_hill_radius_km: float
    hill_fraction: float
    solar_tide_to_earth_gravity_ratio: float
    classification: str
    n_samples: int


def screen_orbit(
    system: cr3bp.CR3BPSystem,
    state0: NDArray[np.float64],
    period: float,
    *,
    n_samples: int = 2000,
    rtol: float = 1e-11,
    atol: float = 1e-11,
) -> HillScreenResult:
    """Pre-screen a CR3BP periodic orbit for real-ephemeris (V4) survivability.

    Propagates ``state0`` over one full ``period`` with DOP853 dense output,
    finds the orbit's farthest excursion from the **primary** (Earth, located at
    ``(-mu, 0, 0)`` in the rotating frame), maps it to physical km via
    ``system.l_km``, and classifies the resulting Earth-Sun Hill fraction.

    Using the true propagated max excursion (not just the IC distance) makes the
    screen conservative: it sees the orbit's apoapsis even when the IC is near
    periapsis. For ``branch_C32_b0`` the max (1.157e6 km) is ~0.4% above the
    IC-only norm the #389 V4 ``structural_diagnostic`` recorded (1.152e6 km).

    Parameters
    ----------
    system :
        CR3BP system (only ``mu`` and ``l_km`` are read).
    state0 :
        6-vector rotating-frame nondimensional IC.
    period :
        Full nondimensional period.
    n_samples :
        Dense-output sample count over ``[0, period]`` for the max search.
    rtol, atol :
        DOP853 tolerances for the propagation.

    Returns
    -------
    HillScreenResult

    Raises
    ------
    ValueError
        If ``state0`` is not length 6 or ``period`` is non-positive / non-finite.
    RuntimeError
        If the propagation fails (e.g. a collision trajectory).
    """
    state0 = np.asarray(state0, dtype=np.float64).reshape(-1)
    if state0.shape != (6,):
        raise ValueError(f"state0 must have shape (6,); got {state0.shape}")
    if not (period > 0.0) or not np.isfinite(period):
        raise ValueError(f"period must be > 0 finite; got {period}")
    if n_samples < 2:
        raise ValueError(f"n_samples must be >= 2; got {n_samples}")

    sol = solve_ivp(
        cr3bp.cr3bp_eom,
        (0.0, period),
        state0,
        args=(system.mu,),
        rtol=rtol,
        atol=atol,
        method="DOP853",
        dense_output=True,
    )
    if not sol.success or sol.sol is None:
        raise RuntimeError(f"hill_screen propagation failed: {sol.message}")

    ts = np.linspace(0.0, period, n_samples)
    ys = sol.sol(ts)  # shape (6, n_samples)
    earth = np.array([-system.mu, 0.0, 0.0], dtype=np.float64)
    dist_from_earth = np.linalg.norm(ys[:3, :].T - earth, axis=1)
    max_amp_nondim = float(np.max(dist_from_earth))
    max_amp_km = max_amp_nondim * float(system.l_km)

    r_hill_km = earth_sun_hill_radius_km()
    hill_fraction = max_amp_km / r_hill_km

    # Linearized solar tide vs Earth gravity at the apoapsis distance.
    gm_earth = float(PRIMARIES["Earth"])
    a_es_km = float(PLANETS["E"].sma_au) * AU_KM
    a_earth_g = gm_earth / max_amp_km**2
    a_sun_tide = 2.0 * MU_SUN_KM3_S2 * max_amp_km / a_es_km**3
    tide_ratio = a_sun_tide / a_earth_g

    return HillScreenResult(
        max_amplitude_km=max_amp_km,
        earth_sun_hill_radius_km=r_hill_km,
        hill_fraction=hill_fraction,
        solar_tide_to_earth_gravity_ratio=tide_ratio,
        classification=_classify(hill_fraction),
        n_samples=n_samples,
    )
