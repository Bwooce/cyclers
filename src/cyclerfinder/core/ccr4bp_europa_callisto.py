"""CCR4BP system constructor for Jupiter Europa-Callisto (`#703`).

`#693`'s cross-solar-system moon-pair screening pass
(``docs/notes/2026-07-23-693-ccr4bp-moonpair-screening.md``), given its own
dedicated deep literature check by `#700`
(``docs/notes/2026-07-24-700-europa-callisto-deep-litcheck.md``), ranked
Jupiter Europa-Callisto as the last remaining `#693` CCR4BP discovery
candidate for the now-proven `#689`-`#694` pipeline (`#702`'s ``ghost_guard``
fix included): ``mu_pert=5.667e-5`` (comparable order to the already-validated
Jupiter-Europa-Ganymede (JEG) system's own ``mu_gan``), novelty-CLEAR (no
direct Europa-Callisto CCR4BP/PCCFBP work found after `#700`'s independent
``pdftotext``-verified check of both in-corpus Kumar-group papers plus five
further live-search queries), and a comfortably collision-clear base orbit
findable with WIDER safety margins than Io-Ganymede's own build needed
(`#700`'s own quick numeric check, re-derived properly by this task's own
test module rather than just cited).

**Real, flagged scientific-motivation caveat (`#700`'s own honest finding,
carried forward here, not smoothed over)**: unlike every other CCR4BP
candidate this project has built (JEG's literature-sourced Europa 3:4, Io-
Europa's exact 2:1, Io-Ganymede's 4.06 near-4:1 Laplace-chain ratio, Umbriel-
Titania's 2.10 near-2:1), Europa:Callisto's real period ratio is ~4.734/4.73
with NO clean low-integer commensurability (best convergents 4/1 off by 0.70,
5/1 off by 0.30, next tier 14/3 or 19/4) -- every seed orbit tested converges
to essentially the SAME near-degenerate near-circular family at whatever
radius Kepler's third law assigns it, rather than a physically-motivated
resonant orbit. This is a genuine difference in scientific motivation from
every prior CCR4BP build in this project's arc, not a build blocker.

This module supplies ONLY the new ``CCR4BPSystem`` parameter set -- it adds NO
new EOM/STM/corrector code. :class:`~cyclerfinder.core.ccr4bp.CCR4BPSystem` is
already fully system-agnostic (`#689`'s own design), so this is a direct
structural analogue of
:func:`cyclerfinder.core.ccr4bp.jupiter_europa_ganymede_default` (and its
Galilean siblings `#695`'s ``core.ccr4bp_io_europa`` / `#696`'s
``core.ccr4bp_io_ganymede`` / Uranian sibling `#701`'s
``core.ccr4bp_umbriel_titania``) with the moon roles:

  * Europa plays the SAME base/forced-body structural role it already plays
    in the JEG system (``mu`` here == Europa's mass ratio, identical to
    :func:`cyclerfinder.core.ccr4bp.jupiter_europa_ganymede_default`'s own
    ``mu``).
  * Callisto plays Ganymede's outer-perturber structural role (``mu_gan``
    here == Callisto's nondimensional mass in the Jupiter-Europa mass unit;
    ``a_gan`` == Callisto's SMA in Europa-SMA units; ``omega_gan`` ==
    Callisto's synodic rate in the Jupiter-Europa synodic frame).

Sourcing (mirrors `#689`'s own discipline exactly)
---------------------------------------------------
Reuses the SAME ``core.satellites`` JPL SSD registry `#689`/`#695`/`#696`/
`#701` all already draw from -- no independently re-sourced constant set
(cross-checked in this module's own test file against BOTH
``core.satellites`` directly and `#693`'s/`#700`'s own independently-derived
screening numbers, ``mu_base=2.528e-5``, ``mu_pert=5.667e-5``,
``a_pert=2.8054``, period ratio ``4.734``):

  * ``mu``     = GM_Europa / (GM_Jupiter_sys + GM_Europa)     (== JEG's own mu)
  * ``mu_gan`` = GM_Callisto / (GM_Jupiter_sys + GM_Europa)
  * ``a_gan``  = SMA_Callisto / SMA_Europa                    (~2.8054)
  * ``omega_gan`` = :func:`cyclerfinder.core.ccr4bp.two_body_synodic_rate`
    (the SAME physically-correct two-body-mean-motion formula `#689` built and
    proved makes the ``mu -> 0`` outer-perturber reduction exact -- reused
    UNMODIFIED, not re-derived).

``GM_Jupiter_sys`` (``core.satellites.PRIMARIES["Jupiter"]``, 1.26686534e8
km^3/s^2) is the JPL DE440 *system* GM, the SAME approximation `#689`'s own
``jupiter_europa_ganymede_default`` already accepts and documents.

``L_KM`` / ``v_unit_km_s`` (this system's own physical-unit conversions) --
a coincidence, not a fix, verified not assumed
------------------------------------------------------------------------------
`#695` discovered that `#694`'s ``search/ccr4bp_heteroclinic_search.py``
hardcodes its km-denominated output fields to EUROPA's SMA/GM (module-level
``_L_KM = 671_100.0`` and a ``_v_unit_km_s`` that recomputes
``GM_Jupiter+GM_Europa`` regardless of the system actually passed in) -- a
real, documented, NOT-YET-FIXED generality gap in that module (left unmodified
per this task's own reuse-only scope), requiring every prior non-JEG build
(`#695` Io-Europa, `#696` Io-Ganymede, `#701` Umbriel-Titania) to independently
rescale every km-denominated field by ``L_KM / 671_100.0``. **For THIS system
the base moon IS Europa itself** (unlike every other non-JEG build), so
``L_KM`` here equals ``671_100.0`` and the correction factor equals ``1.0`` --
this module's own test file explicitly verifies this equality rather than
merely asserting it, so `#703`'s own driver script can report both
``module_native_*`` and ``corrected_*`` fields (for consistency with the other
three discovery-attempt scripts) while confirming they agree exactly, not by
construction-only reasoning.

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

import math

from cyclerfinder.core.ccr4bp import CCR4BPSystem, two_body_synodic_rate
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

# Physical length unit for THIS system: Europa's semi-major axis about
# Jupiter (km). Identical to the JEG system's own length unit (also Europa's
# SMA) -- unlike the Io-based and Umbriel-based non-JEG systems this project
# has already built, whose base moon is NOT Europa. See module docstring.
L_KM = SATELLITES["Europa"].sma_km


def jupiter_europa_callisto_default() -> CCR4BPSystem:
    """CCR4BP parameters for Jupiter-Europa-Callisto from the JPL SSD registry.

    Reuses the SAME sourced constants and formulae as
    :func:`cyclerfinder.core.ccr4bp.jupiter_europa_ganymede_default`, with
    Europa UNCHANGED as the base/forced moon (same role, same mu) and
    Callisto playing Ganymede's perturber role (see module docstring):

      * mu     = GM_Europa / (GM_Jupiter_sys + GM_Europa)     (== JEG's mu)
      * mu_gan = GM_Callisto / (GM_Jupiter_sys + GM_Europa)
      * a_gan  = SMA_Callisto / SMA_Europa                    (~2.8054)
      * omega_gan = two-body synodic rate (`#689`'s own
        :func:`~cyclerfinder.core.ccr4bp.two_body_synodic_rate`, reused
        unmodified)

    The physical Europa:Callisto mean-motion ratio is ~4.73 -- the loosest,
    least commensurate period ratio of any CCR4BP pair this project has built
    (contrast JEG's 2.03, Io-Europa's exact 2.00, Io-Ganymede's 4.06,
    Umbriel-Titania's 2.10); not imposed as any low-integer idealisation here
    (faithful to the sourced physical constants, per the same precedent every
    prior system constructor in this project sets).
    """
    gm_j = PRIMARIES["Jupiter"]
    gm_europa = SATELLITES["Europa"].mu_km3_s2
    gm_callisto = SATELLITES["Callisto"].mu_km3_s2
    sma_europa = SATELLITES["Europa"].sma_km
    sma_callisto = SATELLITES["Callisto"].sma_km
    denom = gm_j + gm_europa  # mass unit G(m1+m2), m2 = Europa (same base role as JEG)
    mu = gm_europa / denom
    mu_gan = gm_callisto / denom
    a_gan = sma_callisto / sma_europa
    omega_gan = two_body_synodic_rate(mu, mu_gan, a_gan)
    return CCR4BPSystem(mu=mu, mu_gan=mu_gan, a_gan=a_gan, omega_gan=omega_gan, theta_gan0=0.0)


def v_unit_km_s() -> float:
    """Physical velocity unit (km/s) for THIS system: ``L_KM * n_europa``.

    ``n_europa = sqrt((GM_Jupiter_sys + GM_Europa) / L_KM**3)`` -- the
    Jupiter-Europa mean motion the CCR4BP nondimensionalises time by (frame
    rate == 1 == n_europa in these units). Provided so downstream
    physical-unit reporting for THIS system does not have to reuse `#694`'s
    Europa-hardcoded ``ccr4bp_heteroclinic_search._v_unit_km_s`` even though,
    for this particular system, the two are expected to agree exactly (see
    module docstring / ``L_KM`` NOTE) -- kept as an independent
    recomputation, not a re-export, so the agreement is a VERIFIED fact, not
    an assumed one."""
    gm_j = PRIMARIES["Jupiter"]
    gm_europa = SATELLITES["Europa"].mu_km3_s2
    n_europa = math.sqrt((gm_j + gm_europa) / L_KM**3)
    return L_KM * n_europa
