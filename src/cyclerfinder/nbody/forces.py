"""Rails third-body force model for the restricted n-body harness (design §2).

**Restricted n-body, heliocentric, ephemeris-driven point-mass perturbers.** The
spacecraft is a massless test particle; the planets are *on rails* from the shared
DE440 BSP (no back-reaction — exact for a massless particle). At a query time the
acceleration on the spacecraft at heliocentric position ``r`` is

    a = -mu_sun * r_hat / r^2
        + sum_p mu_p * ( (r_p - r) / |r_p - r|^3  -  r_p / |r_p|^3 )

i.e. the central Sun term plus the standard third-body **indirect** perturbation
(the second bracket term corrects for the Sun-centred — non-inertial —
heliocentric frame). ``mu`` values come from ``core/constants`` (``MU_SUN_KM3_S2``,
``PLANETS[code].mu_km3_s2``). Perturber states are read through
:func:`ingest_planet_state`, which is the §0-conversion boundary (here the
identity over ``Ephemeris('astropy').state`` — the rung's independence is the
*integrator*, not the reader; design §4 / §5.4 anchor).

The force is **velocity-independent** (no drag / SRP / GR, design §2), which IAS15
handles well; the energy-conservation gate (Task A.4/A.5) checks this rather than
assuming it.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from cyclerfinder.core.constants import MU_SUN_KM3_S2, PLANETS

if TYPE_CHECKING:
    from cyclerfinder.core.ephemeris import Ephemeris

Vec3 = NDArray[np.float64]


def ingest_planet_state(body: str, t_sec: float, ephem: Ephemeris) -> tuple[Vec3, Vec3]:
    """Read a perturber's heliocentric J2000-ecliptic state (the §0 boundary).

    For the shared-DE440 path this is the identity over
    ``Ephemeris('astropy').state(body, t_sec)``: astropy already returns
    heliocentric J2000-ecliptic km / km/s (``ephemeris.py:269-279``), the exact
    frame and axis the harness uses, so no extra rotation is applied. This is
    *why* the GOLDEN GATE 3 anchor (Task A.5) passes to numerical precision — the
    rung's independence is the integrator (REBOUND), not the reader (design §4).
    """
    r, v = ephem.state(body, t_sec)
    return (
        np.asarray(r, dtype=np.float64),
        np.asarray(v, dtype=np.float64),
    )


def rails_acceleration(
    r_km: Vec3,
    t_sec: float,
    bodies: Sequence[str],
    ephem: Ephemeris | None,
) -> Vec3:
    """Heliocentric acceleration on the massless spacecraft (km/s^2).

    Central Sun inverse-square term plus the third-body indirect perturbation from
    each body in ``bodies`` (read on rails via :func:`ingest_planet_state`). With
    no perturbers (``bodies=()``) this reduces to the pure two-body central force —
    the basis of GOLDEN GATE 1 (Sun-only ≡ ``core/kepler.propagate``).
    """
    r = np.asarray(r_km, dtype=np.float64)
    r_norm = float(np.linalg.norm(r))
    a = -MU_SUN_KM3_S2 * r / (r_norm**3)

    for body in bodies:
        if ephem is None:
            raise ValueError(
                "rails_acceleration: bodies were requested but no ephem was "
                "supplied to read their on-rails states"
            )
        mu_p = PLANETS[body].mu_km3_s2
        r_p, _ = ingest_planet_state(body, t_sec, ephem)
        d = r_p - r
        d_norm = float(np.linalg.norm(d))
        rp_norm = float(np.linalg.norm(r_p))
        # Direct attraction toward the perturber + indirect (Sun-frame) term.
        a = a + mu_p * (d / (d_norm**3) - r_p / (rp_norm**3))

    return np.asarray(a, dtype=np.float64)


__all__ = ["ingest_planet_state", "rails_acceleration"]
