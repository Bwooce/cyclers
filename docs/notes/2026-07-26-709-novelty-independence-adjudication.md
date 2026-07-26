# #709 — Independence adjudication: #312 family vs #708 CCR4BP torus-homoclinic

Adversarial check (Fable, 2026-07-26): are the catalogue's two novel-status findings —
`umbriel-titania-1-1-uranian-quasi-cycler-2026` (#312/#563/#569 family representative) and
`umbriel-1-2-torus-homoclinic-uranus-2026` (#689→#708 CCR4BP chain) — two distinct discoveries,
or one structure double-counted at two fidelities?

**VERDICT: GENUINELY INDEPENDENT.** Two distinct mathematical objects with disjoint phase-space
geometry, different base orbits, different closure mechanisms, different equations of motion,
and zero initial-condition/seed lineage. Safe to count as two.

## Numeric comparison (all values from primary artifacts, not task prose)

Sources: `data/gauntlet_566_five_representatives.jsonl` (candidate
`enum563-line26-umbriel-titania-umbriel`), `data/found/701_ccr4bp_umbriel_titania_search/result.json`,
`data/found/705_ccr4bp_epoch_robustness_scan/result.json`, both full catalogue rows,
`core/cr3bp.py::cr3bp_eom`, `core/ccr4bp.py::ccr4bp_eom`, `core/ccr4bp_umbriel_titania.py`,
`scripts/screen_701_ccr4bp_umbriel_titania_search.py`.

| Axis | #312 rep (UTU symmetric closure) | #708 (torus homoclinic) |
|---|---|---|
| Base structure | Two Lambert arcs Umbriel→Titania→Umbriel, tof 3.9545 d each, n=1, n_rev=(0,0), rel_offset=0° | 2-D invariant torus around spacecraft:Umbriel=1:2-EXTERIOR resonant periodic orbit (seeded from Keplerian a=2^(2/3), corrected in Uranus-Umbriel CR3BP) |
| Radial extent | Touches r=265,986 km (Umbriel flyby) and r=436,298 km (Titania flyby) | r ∈ [392,861, 451,644] km; never within 126,875 km of Umbriel's orbit radius |
| Moon encounters | REQUIRED: flybys at both moons, V∞ 1.2296 / 1.0058 km/s, ~6° gravity-assist turn ⇒ periapsis at ~1,000 km scale (Umbriel GM 85.1 km³/s², SOI ≈ 3,102 km) | NONE: min Umbriel distance 181,244 km = 58× SOI; min Titania distance 534,639–828,573 km across all 10 #705 epochs (beyond Titania's own orbital radius) |
| Closure mechanism | Gravity-assist turns at each flyby + symmetric perpendicular-crossing (Miele) condition; residual 4.44e-15 km/s | Stable/unstable manifold asymptotics of the SAME torus; one-shot excursion t_u=12.534 d / t_s=12.593 d; corrected gap 5.66e-10 km / 5.11e-14 km/s |
| Equations of motion | Patched-conic/Lambert with per-leg Uranus-moon CR3BP V∞ matching (row's own data_gap note); Titania is a flyby TARGET, not a force on the Umbriel-anchored legs | CCR4BP: CR3BP + time-periodic Titania direct+indirect acceleration (`ccr4bp.py::_ganymede_acceleration`, verified in code); reduces to CR3BP only at mu_gan=0. Titania is a continuous FORCE, never a target |
| Spacecraft period | Cycle repeats each U-T synodic period 7.9089 d (two legs) | Base orbit sidereal period 8.2884 d (2× Umbriel's 4.1442 d); torus stroboscopic period 11.9911 TU = 7.9090 d = the FORCING period |
| Seed/IC lineage | #558 grid sweep → #563 direct enumeration (Lambert machinery) | `_resonant_symmetric_orbit(mu, 1, 2)` cold Keplerian seed; grep of the whole #701 pipeline finds ZERO references to #563/#312/enum563/gauntlet_566 data |

## Sharpest double-count hypotheses, tested

1. **Same base resonance?** No. #312's rep is not on a spacecraft:Umbriel resonant orbit at all —
   it is a flyby-re-turned transfer-ellipse chain repeating at the synodic period (7.909 d ≠
   2×T_Umbriel = 8.288 d). #708's base is a 1:2-exterior Umbriel-resonant orbit bracketing
   Titania's SMA. Different objects even before Titania forcing is switched on.
2. **Continuation/perturbation identity?** Structurally impossible. #312's object intersects the
   Umbriel close-encounter set twice per cycle (flyby periapsis ~10³ km); #708's object never
   comes within 58× Umbriel's SOI. A smooth CR3BP→CCR4BP continuation cannot delete two flybys
   per period — that is a topology change, not a fidelity refinement. Confirmed by the #566
   V4-strict gauntlet itself: the #312 rep propagated in a FULL n-body model *including Titania
   forcing* remains an encounter tour (bounded drift 289,483 km), i.e., the high-fidelity image
   of #312's rep already exists and is not #708's object. The two structures coexist as distinct
   objects in the same superset model.
3. **Is #708 implicitly contained in the #563 30-member family?** No. The #563 enumeration space
   is by construction {12 distinct-moon-pair flyby directions} × {n ≤ n_max} × {n_rev ∈ {0..3}²} ×
   {0°,180°} Lambert symmetric closures with mandatory per-encounter bend gates — every member has
   two real flybys. An encounter-free torus-homoclinic is not representable in that model (patched
   conic has no invariant tori, no manifolds, no Titania forcing term), and #563's own write-up
   claims exhaustiveness only "FOR THE SYMMETRIC-CLOSURE CLASS specifically."

## The one genuine numeric coincidence

Both rows carry 7.909 d: #312's tour repeat period and #708's torus stroboscopic period are both
the Umbriel-Titania synodic period. This is SYSTEM-level, not object-level: the CCR4BP stroboscopic
map period is the forcing period by definition, and #563's commensurate closures are built at
tof = n·T_syn/2 by construction. Every synodic-locked structure in this three-body pair shares this
clock; it is not evidence of shared identity.

## Honest caveats

- No numerical homotopy/continuation experiment was run between the two objects; the
  impossibility argument is geometric (58× SOI encounter-distance gap, 126,875 km radial-range
  gap, opposite closure mechanisms), which I judge decisive without it.
- #701's result.json carries a documented `physical_unit_caveat` (the reused heteroclinic-search
  module hardcodes Europa's L_KM for km conversions; the #701 driver recomputed km quantities in
  Umbriel units). I spot-checked TU/L consistency (torus period TU→days reproduces the row's
  synodic period to 1e-15) but did not re-verify every km-denominated field.
- The "1:2-exterior" label is nominal: result.json's own note says osculating a~1.698 at t=0,
  extent 1.477–1.698 (near-2:1, not exact — the pair's real period ratio is 2.101). This does not
  affect the verdict; every value of that osculating range stays far above Umbriel's orbit.
- This adjudication addresses INDEPENDENCE only, not either row's absolute novelty vs the
  literature (that is #349/#699/#706's ground, not re-litigated here).
