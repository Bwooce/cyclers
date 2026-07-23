"""CCR4BP system constructor for Jupiter-Io-Europa (`#695`).

`#693`'s cross-solar-system moon-pair screening pass
(``docs/notes/2026-07-23-693-ccr4bp-moonpair-screening.md``) ranked Jupiter
Io-Europa as the TOP novelty-cleared candidate for the now-proven
`#689`-`#694` CCR4BP discovery pipeline: an essentially exact 2.0000 physical
period ratio (the tightest resonance in the solar system), a perturber mass
ratio (``mu_gan`` here, ``2.53e-5``) comparable in order to the already-
validated Jupiter-Europa-Ganymede (JEG) system's own base leg, both
eccentricities/inclinations well under JEG's own reference bar, and no hit in
three targeted literature checks.

This module supplies ONLY the new ``CCR4BPSystem`` parameter set -- it adds
NO new EOM/STM/corrector code. :class:`~cyclerfinder.core.ccr4bp.CCR4BPSystem`
is already fully system-agnostic (`#689`'s own design), so this is a direct
structural analogue of :func:`cyclerfinder.core.ccr4bp.jupiter_europa_ganymede_default`
with the moon roles shifted one rung down the Laplace chain:

  * Io plays Europa's STRUCTURAL role: the base/forced moon whose CR3BP
    resonant periodic orbit seeds the torus (``mu`` here == Io's mass ratio).
  * Europa plays Ganymede's STRUCTURAL role: the outer periodic perturber
    (``mu_gan`` here == Europa's nondimensional mass in the Io-Europa mass
    unit; ``a_gan`` == Europa's SMA in Io-SMA units; ``omega_gan`` == Europa's
    synodic rate in the Jupiter-Io synodic frame).

Sourcing (mirrors `#689`'s own discipline exactly)
---------------------------------------------------
Reuses the SAME ``core.satellites`` JPL SSD registry `#689` and `#688`'s
``genome/composed_moon_map.py`` both already draw from -- no independently
re-sourced constant set:

  * ``mu``     = GM_Io / (GM_Jupiter_sys + GM_Io)
  * ``mu_gan`` = GM_Europa / (GM_Jupiter_sys + GM_Io)
  * ``a_gan``  = SMA_Europa / SMA_Io                       (~1.591)
  * ``omega_gan`` = :func:`cyclerfinder.core.ccr4bp.two_body_synodic_rate`
    (the SAME physically-correct two-body-mean-motion formula `#689` built and
    proved makes the ``mu -> 0`` Jupiter-perturber reduction exact -- reused
    UNMODIFIED, not re-derived).

Cross-checked against `#693`'s own screening table (``mu_base=4.704e-5``,
``mu_pert=2.528e-5``, period ratio 2.0000) in this module's own test file,
``tests/core/test_ccr4bp_io_europa.py``.

Scope / discipline
-------------------
Pure parameter-construction module. Does NOT modify
:mod:`cyclerfinder.core.ccr4bp`, :mod:`cyclerfinder.search.variational_ccr4bp_torus`,
:mod:`cyclerfinder.search.ccr4bp_whisker`,
:mod:`cyclerfinder.search.ccr4bp_manifold_globalize`, or
:mod:`cyclerfinder.search.ccr4bp_heteroclinic_search` -- every later pipeline
stage consumes the returned :class:`~cyclerfinder.core.ccr4bp.CCR4BPSystem`
through those modules' existing, unmodified public API.
"""

from __future__ import annotations

from cyclerfinder.core.ccr4bp import CCR4BPSystem, two_body_synodic_rate
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

# Physical length unit for THIS system: Io's semi-major axis about Jupiter
# (km). NOTE this is DIFFERENT from the JEG system's own length unit (Europa's
# SMA, 671,100 km) -- any physical-unit (km) conversion reused from `#694`'s
# ``search/ccr4bp_heteroclinic_search.py`` (whose ``_L_KM``/``_v_unit_km_s``
# are hardcoded to Europa's SMA) must be rescaled by ``L_KM / 671_100.0``
# before being trusted for THIS system -- see this task's own driver script
# for the explicit correction (`#695`'s own finding: those two helpers are not
# actually system-agnostic despite every other `#694` function being so).
L_KM = SATELLITES["Io"].sma_km


def jupiter_io_europa_default() -> CCR4BPSystem:
    """CCR4BP parameters for Jupiter-Io-Europa from the JPL SSD registry.

    Reuses the SAME sourced constants and formulae as
    :func:`cyclerfinder.core.ccr4bp.jupiter_europa_ganymede_default`, with Io
    playing Europa's base-moon role and Europa playing Ganymede's perturber
    role (see module docstring):

      * mu     = GM_Io     / (GM_Jupiter_sys + GM_Io)
      * mu_gan = GM_Europa / (GM_Jupiter_sys + GM_Io)
      * a_gan  = SMA_Europa / SMA_Io                          (~1.591)
      * omega_gan = two-body synodic rate (`#689`'s own
        :func:`~cyclerfinder.core.ccr4bp.two_body_synodic_rate`, reused
        unmodified)

    The physical Io:Europa mean-motion ratio is ``n_io/n_europa ~ 2.007`` --
    the tightest (closest-to-exact) commensurability of any moon pair
    surveyed in `#693`'s screening pass, though still not imposed as an exact
    2:1 idealisation here (faithful to the sourced physical constants, per
    the same precedent `#689` set for Europa-Ganymede's own 2.014 ratio).
    """
    gm_j = PRIMARIES["Jupiter"]
    gm_io = SATELLITES["Io"].mu_km3_s2
    gm_eu = SATELLITES["Europa"].mu_km3_s2
    sma_io = SATELLITES["Io"].sma_km
    sma_eu = SATELLITES["Europa"].sma_km
    denom = gm_j + gm_io  # mass unit G(m1+m2), m2 = Io (the base moon here)
    mu = gm_io / denom
    mu_gan = gm_eu / denom
    a_gan = sma_eu / sma_io
    omega_gan = two_body_synodic_rate(mu, mu_gan, a_gan)
    return CCR4BPSystem(mu=mu, mu_gan=mu_gan, a_gan=a_gan, omega_gan=omega_gan, theta_gan0=0.0)


def v_unit_km_s() -> float:
    """Physical velocity unit (km/s) for THIS system: ``L_KM * n_io``.

    ``n_io = sqrt((GM_Jupiter_sys + GM_Io) / L_KM**3)`` -- the Jupiter-Io mean
    motion the CCR4BP nondimensionalises time by (frame rate == 1 == n_io in
    these units). Provided so downstream physical-unit reporting for THIS
    system does not have to reuse `#694`'s Europa-hardcoded
    ``ccr4bp_heteroclinic_search._v_unit_km_s`` (see module docstring / ``L_KM``
    NOTE)."""
    import math

    gm_j = PRIMARIES["Jupiter"]
    gm_io = SATELLITES["Io"].mu_km3_s2
    n_io = math.sqrt((gm_j + gm_io) / L_KM**3)
    return L_KM * n_io
