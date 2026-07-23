"""Jupiter-Io-Ganymede :class:`~cyclerfinder.core.ccr4bp.CCR4BPSystem` constructor (#696).

`#693`'s cross-solar-system moon-pair screening pass
(``docs/notes/2026-07-23-693-ccr4bp-moonpair-screening.md``) ranked Jupiter
Io-Ganymede as the second-strongest novelty-cleared CCR4BP discovery candidate
after Io-Europa (`#695`, a separate, parallel task -- do NOT touch that pair
here): ``mu_pert=7.80e-5`` is essentially IDENTICAL forcing strength to the
already-validated Jupiter-Europa-Ganymede (JEG, `#689`-`#694`) system's own
``mu_gan`` perturber term (since both involve Ganymede/Jupiter), a 4.06 period
ratio (a real 4:2:1 Laplace-chain consequence, not exact), both moons'
eccentricities/inclinations well under `#693`'s own tractability bar, and NOT
FOUND in `#693`'s own targeted literature checks.

This module is `#689`'s ``jupiter_europa_ganymede_default()`` construction
pattern applied to Io-Ganymede instead of Europa-Ganymede -- Io plays Europa's
STRUCTURAL role (the base/forced body the torus is seeded from), Ganymede
plays the SAME outer-perturber role it already plays in the JEG system (same
body, same mass, now perturbing Io's orbit instead of Europa's). Deliberately
a SEPARATE module (not appended to ``core/ccr4bp.py`` itself) so this task and
the parallel `#695` (Io-Europa) task cannot collide on the same file while
both add a new system constructor concurrently -- `#689`'s own
:class:`~cyclerfinder.core.ccr4bp.CCR4BPSystem` and
:func:`~cyclerfinder.core.ccr4bp.two_body_synodic_rate` are reused UNMODIFIED,
imported from ``core.ccr4bp``.

Sourced constants (JPL SSD ``core.satellites`` registry -- cross-checked
against `#693`'s own screening table and `#689`'s own registry access pattern;
see ``tests/core/test_ccr4bp_io_ganymede.py`` for the numeric cross-check
against both `#693`'s table and ``genome.composed_moon_map.moon_config``):

    mu        = GM_Io / (GM_Jupiter_sys + GM_Io)        ~ 4.7042e-5
    mu_gan    = GM_Ganymede / (GM_Jupiter_sys + GM_Io)   ~ 7.8046e-5
    a_gan     = SMA_Ganymede / SMA_Io                    ~ 2.5377
    omega_gan = two_body_synodic_rate(mu, mu_gan, a_gan)  ~ -0.7526

The physical Io:Ganymede mean-motion ratio (``n_Ganymede/n_Io``) implied by
the two-body synodic rate is ``1 + omega_gan ~ 0.2474``, i.e. a period ratio
``T_Ganymede/T_Io ~ 4.043`` from the pure two-body SMA-Kepler relation used
here -- close to but not identical to `#693`'s own externally-sourced
JPL-fitted value of 4.059 (a ~0.4% gap of the SAME kind `#689`'s own
docstring documents for JEG's 2.014 (physical/JPL-fitted) vs 2.030
(two-body-Kepler-on-SMA) Europa:Ganymede ratio -- the moons' own masses and
real (non-two-body) orbital fits account for the residual, not a modeling
error here). This model does not require exact commensurability (`#690`'s own
established precedent, reused verbatim); the CCR4BP is time-periodic at
Ganymede's synodic rate regardless.
"""

from __future__ import annotations

from cyclerfinder.core.ccr4bp import CCR4BPSystem, two_body_synodic_rate
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES


def jupiter_io_ganymede_default() -> CCR4BPSystem:
    """CCR4BP parameters for Jupiter-Io-Ganymede from the JPL SSD registry.

    Direct analogue of :func:`cyclerfinder.core.ccr4bp.jupiter_europa_ganymede_default`,
    with Io in Europa's structural role (the base/forced moon) and Ganymede
    unchanged as the outer perturber. Reuses the SAME ``core.satellites``
    registry and the SAME two-body synodic-rate formula
    (:func:`cyclerfinder.core.ccr4bp.two_body_synodic_rate`) so the two
    constructors cannot drift apart -- cross-checked in
    ``tests/core/test_ccr4bp_io_ganymede.py``.
    """
    gm_j = PRIMARIES["Jupiter"]
    gm_io = SATELLITES["Io"].mu_km3_s2
    gm_g = SATELLITES["Ganymede"].mu_km3_s2
    sma_io = SATELLITES["Io"].sma_km
    sma_g = SATELLITES["Ganymede"].sma_km
    denom = gm_j + gm_io  # mass unit G(m1+m2), m2 = Io here (Europa's role in #689)
    mu = gm_io / denom
    mu_gan = gm_g / denom
    a_gan = sma_g / sma_io
    omega_gan = two_body_synodic_rate(mu, mu_gan, a_gan)
    return CCR4BPSystem(mu=mu, mu_gan=mu_gan, a_gan=a_gan, omega_gan=omega_gan, theta_gan0=0.0)
