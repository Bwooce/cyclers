"""Jupiter-Callisto-Ganymede :class:`~cyclerfinder.core.ccr4bp.CCR4BPSystem` constructor (`#715`).

Second-positive-control build for the CCR4BP manifold/heteroclinic-search
pipeline (`#689`-`#694`), grounded in Aryan & Fitzgerald (2024), "Four Body
Invariant Structures And Chaos Analysis For Jovian Multi-Moon Ballistic
Transfers," AAS 24-103 (digested in full at
``docs/notes/2026-07-26-710-digest-aryan-fitzgerald-2024-jovian-pccfbp.md``).

That paper's Table 2 system is Jupiter-Callisto-Ganymede: Callisto plays the
"second primary"/base-moon structural role (their ``mu1``), Ganymede is the
periodic PERTURBER (their ``mu2``) -- i.e. the SAME base+perturber structure
this project's own ``core/ccr4bp.py`` already implements, just with Callisto
(not Europa/Io) as the base moon.

**The key structural wrinkle this module deliberately confronts, not
elides**: in every CCR4BP system this project has built so far (Jupiter-
Europa-Ganymede, Io-Europa, Io-Ganymede, Umbriel-Titania), the perturber is
OUTER to the base moon (``a_gan > 1``). Here Ganymede orbits INSIDE
Callisto's own orbit (``a_gan = SMA_Ganymede / SMA_Callisto ~ 0.569 < 1``) --
the paper's own base+perturber pairing puts the FARTHER moon in the
frame-defining "second primary" role and the CLOSER moon as the perturber,
the opposite radial arrangement from this project's prior builds. Nothing in
``core/ccr4bp.py``'s ``CCR4BPSystem``/``ccr4bp_eom``/``two_body_synodic_rate``
requires ``a_gan > 1`` -- the direct+indirect Ganymede acceleration formula
and the two-body synodic-rate formula are both algebraically indifferent to
whether the perturber is inside or outside the base moon's orbit, ONLY the
docstring comment `"(> 1)"` on :attr:`CCR4BPSystem.a_gan` implicitly assumed
the previously-built cases -- this module is the first to exercise
``a_gan < 1`` and ``omega_gan > 0`` (the perturber is now FASTER than the
base moon's frame rate, so its synodic angle ADVANCES rather than
regresses), and ``tests/core/test_ccr4bp_callisto_ganymede.py`` verifies the
``mu_gan -> 0`` structural reduction and the sign of ``omega_gan`` explicitly
rather than assuming the existing sibling-module test pattern transfers
unchanged.

This module supplies ONLY the new ``CCR4BPSystem`` parameter set -- no new
EOM/STM/corrector code, mirroring `#703`'s ``core/ccr4bp_europa_callisto.py``
and `#696`'s ``core/ccr4bp_io_ganymede.py`` reuse-only discipline exactly.

Sourcing (mirrors `#689`'s own discipline)
-------------------------------------------
Reuses the SAME ``core.satellites`` JPL SSD registry every prior CCR4BP
system constructor draws from:

    mu        = GM_Callisto / (GM_Jupiter_sys + GM_Callisto)   ~ 5.6667e-5
    mu_gan    = GM_Ganymede / (GM_Jupiter_sys + GM_Callisto)    ~ 7.8045e-5
    a_gan     = SMA_Ganymede / SMA_Callisto                     ~ 0.5685
    omega_gan = two_body_synodic_rate(mu, mu_gan, a_gan)         > 0

Cross-checked in this module's own test file against the paper's own quoted
Table 2 values (``mu1=5.6623e-05``, ``mu2=7.7890e-05``) -- this project's
independently-sourced JPL DE440 GM/SMA registry reproduces both to within
~0.1%/0.02% relative, the same kind of small cross-source residual `#689`'s
own docstring documents for the Europa-Ganymede pair against Kumar 2021's
quoted values.

``L_KM`` / ``v_unit_km_s`` -- this system's own physical-unit conversions,
following `#703`'s precedent that `#694`'s ``search/ccr4bp_heteroclinic_
search.py`` hardcodes Europa's SMA/GM into its km-denominated output fields
(``_L_KM``, ``_v_unit_km_s``) -- a documented, NOT-fixed generality gap in
that module, left unmodified here (this task's own new glue module,
``search/ccr4bp_chained_transfer.py``, does its OWN physical-unit conversion
using this module's ``L_KM``/``v_unit_km_s``, not `#694`'s hardcoded ones).
"""

from __future__ import annotations

import math

from cyclerfinder.core.ccr4bp import CCR4BPSystem, two_body_synodic_rate
from cyclerfinder.core.satellites import PRIMARIES, SATELLITES

# Physical length unit for THIS system: Callisto's semi-major axis about
# Jupiter (km) -- Callisto plays the base/frame-defining "second primary"
# role here (mirrors #703's L_KM = Europa's SMA for the Europa-Callisto
# system, whose base moon is likewise Europa there).
L_KM = SATELLITES["Callisto"].sma_km


def jupiter_callisto_ganymede_default() -> CCR4BPSystem:
    """CCR4BP parameters for Jupiter-Callisto-Ganymede from the JPL SSD registry.

    Callisto is the base/forced moon (frame-defining "second primary",
    Aryan & Fitzgerald's ``mu1``); Ganymede is the periodic perturber
    (their ``mu2``), here INTERIOR to the base moon's orbit
    (``a_gan ~ 0.569 < 1``, unlike every other CCR4BP system this project has
    built -- see module docstring). Reuses `#689`'s own
    :func:`cyclerfinder.core.ccr4bp.two_body_synodic_rate` UNMODIFIED; that
    formula is algebraically indifferent to whether ``a_gan`` is above or
    below 1, so ``omega_gan`` comes out positive here (Ganymede, the faster
    inner body, ADVANCES in the Callisto-synodic frame, the opposite sign
    convention from the outer-perturber systems already built).
    """
    gm_j = PRIMARIES["Jupiter"]
    gm_callisto = SATELLITES["Callisto"].mu_km3_s2
    gm_ganymede = SATELLITES["Ganymede"].mu_km3_s2
    sma_callisto = SATELLITES["Callisto"].sma_km
    sma_ganymede = SATELLITES["Ganymede"].sma_km
    denom = gm_j + gm_callisto  # mass unit G(m1+m2), m2 = Callisto (base role)
    mu = gm_callisto / denom
    mu_gan = gm_ganymede / denom
    a_gan = sma_ganymede / sma_callisto
    omega_gan = two_body_synodic_rate(mu, mu_gan, a_gan)
    return CCR4BPSystem(mu=mu, mu_gan=mu_gan, a_gan=a_gan, omega_gan=omega_gan, theta_gan0=0.0)


def v_unit_km_s() -> float:
    """Physical velocity unit (km/s) for THIS system: ``L_KM * n_callisto``.

    ``n_callisto = sqrt((GM_Jupiter_sys + GM_Callisto) / L_KM**3)`` -- the
    Jupiter-Callisto mean motion the CCR4BP nondimensionalises time by (base
    frame rate == 1 == n_callisto in these units). Independent
    recomputation, not a re-export of `#694`'s Europa-hardcoded
    ``ccr4bp_heteroclinic_search._v_unit_km_s`` (mirrors `#703`'s own
    precedent and its stated rationale)."""
    gm_j = PRIMARIES["Jupiter"]
    gm_callisto = SATELLITES["Callisto"].mu_km3_s2
    n_callisto = math.sqrt((gm_j + gm_callisto) / L_KM**3)
    return L_KM * n_callisto
