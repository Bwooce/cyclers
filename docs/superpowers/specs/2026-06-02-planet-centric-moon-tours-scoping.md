# Planet-centric (moon-tour) cyclers — scoping note

**Status:** scoped, not started. Deliberately out of scope for #75 (body-agnostic
compute machinery), which generalised the *heliocentric* V/E/M machinery only.

**Why this is its own task:** #75 removed the hardcoded body assumptions inside a
fixed *central body* (the Sun): the `PLANETS` registry, the astropy body-name
map, `synodic_omega`, the Aldrin E-M-E optimiser, and the score taxi-body are now
parameterised. What it did **not** touch is the assumption that there *is one*
central body, that it is the Sun, and that every encounter body orbits it on a
heliocentric ephemeris. Moon tours (e.g. a Galilean-moon cycler around Jupiter, or
an Earth-Moon cycler) break that assumption: the bodies orbit a *planet*, not the
Sun. That is a structural change to the propagation/ephemeris/frame layer, not a
registry edit.

## What changes

1. **Central-body abstraction.** Today `MU_SUN_KM3_S2` is threaded through
   propagation, Lambert, frames, and `orbit_elements_au`. A moon tour propagates
   in a *planet-centric* two-body field with `mu = planet.GM`. The central body
   (and its μ) must become a parameter of the propagation context rather than a
   module constant. Lambert and the rotating-frame transform already take μ as an
   argument in most places; the audit is to find the spots that still reach for
   the Sun constant directly.

2. **Satellite ephemeris source.** `Ephemeris("astropy")` returns *heliocentric*
   planet states. Moon positions relative to their primary are a different data
   source:
   - Earth's Moon is available in astropy's built-in `solar_system_ephemeris`.
   - Jovian / Saturnian satellites are **not** in the built-in kernels; they need
     SPICE satellite kernels (e.g. `jup365.bsp`) read via `jplephem` or
     `spiceypy`. This is a new dependency + data-acquisition step, and the kernels
     are large — decide whether to vendor, download-on-demand, or gate behind an
     extra.
   - A new `Ephemeris` mode (e.g. `Ephemeris("astropy", center="J")`) or a
     distinct provider is the likely shape; the provider must report states in the
     planet-centric frame consistently with the chosen μ.

3. **Units.** `AU_KM` and `*_au` element helpers assume heliocentric scales.
   Moon orbits are in 10^5–10^6 km, not AU. The element/reporting helpers should
   either become unit-agnostic or grow a planet-centric variant; golden anchors
   for moon cyclers will be in km / planet radii, not AU.

4. **Registry.** A `SATELLITES` registry analogous to `PLANETS`, keyed per
   primary (each moon needs: primary, GM, radius, mean motion about the primary,
   safe flyby altitude). `synodic_omega` and the optimiser already accept any body
   present in the registry, so once satellites are registered with the same
   `PlanetData` shape, the optimiser generalises with no further change — provided
   the propagation context carries the right central μ.

## What does NOT change

- The optimiser (`optimise_maintenance_dv`), the turn-deficit geometry
  (`idealized_flyby_turn_deficit`), the scoring taxi-cost, and the
  cycler/encounter/leg model are already body- and length-agnostic after #75. They
  operate on whatever bodies the ephemeris + registry provide; they do not assume
  the Sun. The only coupling left is the central-μ plumbing in items 1–3.

## Open questions (resolve at brainstorming time)

- One central body per cycler, or allow a heliocentric leg that hands off to a
  planet-centric capture (a true "interplanetary + moon tour" chain)? The latter
  is a much larger multi-patched-conic change; recommend starting single-central-body.
- Satellite ephemeris: vendor SPICE kernels vs download-on-demand vs an analytic
  circular-coplanar surrogate (mirrors the existing `Ephemeris("circular")`) for a
  first pass that needs no new data dependency.
- Which target system first? Earth-Moon (data already in astropy, smallest step)
  vs Galilean (richer literature, needs SPICE).

## Recommended first slice

Earth-Moon, analytic-circular ephemeris surrogate first (no new data dependency),
to exercise the central-body μ plumbing end-to-end; then add the astropy Moon
provider; defer SPICE-kernel satellite systems to a follow-up.

## Asteroids — explicitly separate

Asteroid targets are a *different* axis (real eccentric/inclined heliocentric
ephemeris; negligible GM so they are flyby *targets*, not bending bodies). Not
covered here.
