"""CCR4BP system constructor for Uranus-Umbriel-Titania (`#701`).

`#693`'s cross-solar-system moon-pair screening pass
(``docs/notes/2026-07-23-693-ccr4bp-moonpair-screening.md``), independently
re-confirmed by `#699`'s deep novelty check
(``docs/notes/2026-07-24-699-umbriel-titania-deep-litcheck.md``), ranked
Uranus Umbriel-Titania as the strongest non-Jovian CCR4BP discovery candidate
for the now-proven `#689`-`#694` pipeline: ``mu_pert=3.916e-5`` (only ~2x below
the already-validated Jupiter-Europa-Ganymede (JEG) system's own ``mu_gan``
perturber term -- the best-conditioned non-Jovian candidate surveyed), both
eccentricities under the tractability bar (Umbriel 0.004, Titania 0.002),
``Delta i ~ 0`` (both moons within a fraction of a degree of Uranus's own
equatorial/Laplace plane), and a real 2.101 near-2:1 period ratio (comparable
looseness to JEG's own 2.030) to anchor the base resonant orbit. `#699`
cleared this pair on ten independent literature-search queries plus a full-text
check of the closest adjacent paper (Kumar arXiv:2509.03655, whose CCR4BP
secondary-resonance work covers Uranus-Oberon/Titania-Oberon, NOT
Umbriel-Titania) -- no disqualifying prior CCR4BP/torus/heteroclinic work
found.

This module supplies ONLY the new ``CCR4BPSystem`` parameter set -- it adds NO
new EOM/STM/corrector code. :class:`~cyclerfinder.core.ccr4bp.CCR4BPSystem` is
already fully system-agnostic (`#689`'s own design), so this is a direct
structural analogue of
:func:`cyclerfinder.core.ccr4bp.jupiter_europa_ganymede_default` (and its
Galilean siblings `#695`'s ``core.ccr4bp_io_europa`` /`#696`'s
``core.ccr4bp_io_ganymede``) with the primary changed to Uranus and the moon
roles:

  * Umbriel plays Europa's STRUCTURAL role: the base/forced moon whose CR3BP
    resonant periodic orbit seeds the torus (``mu`` here == Umbriel's mass
    ratio).
  * Titania plays Ganymede's STRUCTURAL role: the outer periodic perturber
    (``mu_gan`` here == Titania's nondimensional mass in the Uranus-Umbriel
    mass unit; ``a_gan`` == Titania's SMA in Umbriel-SMA units; ``omega_gan``
    == Titania's synodic rate in the Uranus-Umbriel synodic frame).

Sourcing (mirrors `#689`'s own discipline exactly)
---------------------------------------------------
Reuses the SAME ``core.satellites`` JPL SSD registry `#689`/`#695`/`#696` all
already draw from -- no independently re-sourced constant set (cross-checked
in this module's own test file against BOTH ``core.satellites`` directly and
`#693`'s own independently-sourced screening-pass table):

  * ``mu``     = GM_Umbriel / (GM_Uranus_sys + GM_Umbriel)
  * ``mu_gan`` = GM_Titania / (GM_Uranus_sys + GM_Umbriel)
  * ``a_gan``  = SMA_Titania / SMA_Umbriel                  (~1.640)
  * ``omega_gan`` = :func:`cyclerfinder.core.ccr4bp.two_body_synodic_rate`
    (the SAME physically-correct two-body-mean-motion formula `#689` built and
    proved makes the ``mu -> 0`` outer-perturber reduction exact -- reused
    UNMODIFIED, not re-derived).

``GM_Uranus_sys`` (``core.satellites.PRIMARIES["Uranus"]``,
5.7945564e6 km^3/s^2) is the JPL DE440 *system* GM (Uranus + all five regular
moons); the residual folding of the other four moons' mass into the Uranus
term is the same small approximation ``jupiter_europa_ganymede_default``
already accepts for the Jovian system (documented there), scaled down here
since the Uranian regular satellites are collectively far less massive
relative to their primary than the Galilean moons are to Jupiter.

``L_KM`` / ``v_unit_km_s`` (this system's own physical-unit conversions)
--------------------------------------------------------------------------
`#695` discovered that `#694`'s ``search/ccr4bp_heteroclinic_search.py``
hardcodes its km-denominated output fields to EUROPA's SMA/GM (module-level
``_L_KM = 671_100.0`` and a ``_v_unit_km_s`` that recomputes
``GM_Jupiter+GM_Europa`` regardless of the system actually passed in) -- a
real, documented, NOT-YET-FIXED generality gap in that module (left unmodified
per this task's own reuse-only scope). ``#701``'s own driver script must
independently recompute every km-denominated quantity using Umbriel's OWN
SMA/GM via this module's ``L_KM`` / :func:`v_unit_km_s`, exactly as `#695`'s
and `#696`'s own driver scripts did for Io's SMA/GM.

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

# Physical length unit for THIS system: Umbriel's semi-major axis about Uranus
# (km). DIFFERENT from both the JEG system's own length unit (Europa's SMA,
# 671,100 km) and the Io-Europa/Io-Ganymede systems' (Io's SMA, 421,800 km) --
# any physical-unit (km) conversion reused from #694's
# ``search/ccr4bp_heteroclinic_search.py`` (whose ``_L_KM``/``_v_unit_km_s``
# are hardcoded to Europa's SMA) must be rescaled by ``L_KM / 671_100.0``
# before being trusted for THIS system -- see this task's own driver script
# for the explicit correction (the same #695/#696 finding, re-applied here).
L_KM = SATELLITES["Umbriel"].sma_km


def uranus_umbriel_titania_default() -> CCR4BPSystem:
    """CCR4BP parameters for Uranus-Umbriel-Titania from the JPL SSD registry.

    Reuses the SAME sourced constants and formulae as
    :func:`cyclerfinder.core.ccr4bp.jupiter_europa_ganymede_default`, with
    Umbriel playing Europa's base-moon role and Titania playing Ganymede's
    perturber role (see module docstring):

      * mu     = GM_Umbriel / (GM_Uranus_sys + GM_Umbriel)
      * mu_gan = GM_Titania / (GM_Uranus_sys + GM_Umbriel)
      * a_gan  = SMA_Titania / SMA_Umbriel                     (~1.640)
      * omega_gan = two-body synodic rate (`#689`'s own
        :func:`~cyclerfinder.core.ccr4bp.two_body_synodic_rate`, reused
        unmodified)

    The physical Umbriel:Titania mean-motion ratio is ``n_umbriel/n_titania ~
    2.10`` -- near, but not exactly, a 2:1 commensurability (no Uranian moon
    pair is presently exactly locked), not imposed as an exact 2:1
    idealisation here (faithful to the sourced physical constants, per the
    same precedent `#689` set for Europa-Ganymede's own 2.014-2.030 ratio and
    `#695`/`#696` set for their own Galilean pairs).
    """
    gm_u = PRIMARIES["Uranus"]
    gm_umbriel = SATELLITES["Umbriel"].mu_km3_s2
    gm_titania = SATELLITES["Titania"].mu_km3_s2
    sma_umbriel = SATELLITES["Umbriel"].sma_km
    sma_titania = SATELLITES["Titania"].sma_km
    denom = gm_u + gm_umbriel  # mass unit G(m1+m2), m2 = Umbriel (the base moon here)
    mu = gm_umbriel / denom
    mu_gan = gm_titania / denom
    a_gan = sma_titania / sma_umbriel
    omega_gan = two_body_synodic_rate(mu, mu_gan, a_gan)
    return CCR4BPSystem(mu=mu, mu_gan=mu_gan, a_gan=a_gan, omega_gan=omega_gan, theta_gan0=0.0)


def v_unit_km_s() -> float:
    """Physical velocity unit (km/s) for THIS system: ``L_KM * n_umbriel``.

    ``n_umbriel = sqrt((GM_Uranus_sys + GM_Umbriel) / L_KM**3)`` -- the
    Uranus-Umbriel mean motion the CCR4BP nondimensionalises time by (frame
    rate == 1 == n_umbriel in these units). Provided so downstream
    physical-unit reporting for THIS system does not have to reuse `#694`'s
    Europa-hardcoded ``ccr4bp_heteroclinic_search._v_unit_km_s`` (see module
    docstring / ``L_KM`` NOTE)."""
    import math

    gm_u = PRIMARIES["Uranus"]
    gm_umbriel = SATELLITES["Umbriel"].mu_km3_s2
    n_umbriel = math.sqrt((gm_u + gm_umbriel) / L_KM**3)
    return L_KM * n_umbriel
