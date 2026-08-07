"""CCR4BP system constructor for Saturn-Titan-Hyperion (`#716`).

Targeting the Saturn-Titan-Hyperion 4:3 mean-motion resonance
(P_Hyperion / P_Titan ≈ 1.333). Titan plays Europa's structural role (base moon, mu),
and Hyperion plays Ganymede's structural role (outer perturber, mu_gan).
"""

from __future__ import annotations

import math

from cyclerfinder.core.ccr4bp import CCR4BPSystem, two_body_synodic_rate
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

# Physical length unit for THIS system: Titan's semi-major axis about Saturn (km)
L_KM = SATELLITES["Titan"].sma_km


def saturn_titan_hyperion_default() -> CCR4BPSystem:
    """CCR4BP parameters for Saturn-Titan-Hyperion from the JPL SSD registry (`#716`).

    * mu        = GM_Titan / (GM_Saturn_sys + GM_Titan)
    * mu_gan    = GM_Hyperion / (GM_Saturn_sys + GM_Titan)
    * a_gan     = SMA_Hyperion / SMA_Titan                     (~1.212)
    * omega_gan = two_body_synodic_rate
    """
    gm_s = PRIMARIES["Saturn"]
    gm_titan = SATELLITES["Titan"].mu_km3_s2
    gm_hyperion = SATELLITES["Hyperion"].mu_km3_s2
    sma_titan = SATELLITES["Titan"].sma_km
    sma_hyperion = SATELLITES["Hyperion"].sma_km

    denom = gm_s + gm_titan
    mu = gm_titan / denom
    mu_gan = gm_hyperion / denom
    a_gan = sma_hyperion / sma_titan
    omega_gan = two_body_synodic_rate(mu, mu_gan, a_gan)

    return CCR4BPSystem(mu=mu, mu_gan=mu_gan, a_gan=a_gan, omega_gan=omega_gan, theta_gan0=0.0)


def v_unit_km_s() -> float:
    """Physical velocity unit (km/s) for Saturn-Titan system: ``L_KM * n_titan``."""
    gm_s = PRIMARIES["Saturn"]
    gm_titan = SATELLITES["Titan"].mu_km3_s2
    n_titan = math.sqrt((gm_s + gm_titan) / L_KM**3)
    return L_KM * n_titan
