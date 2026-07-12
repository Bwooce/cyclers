# Digest: Gurfil & Kasdin 2002 — Niching GA characterization of geocentric ER3BP orbits

Single-paper digest. Read 24/24 pages of the user-supplied PDF on 2026-07-12 AET.
Filed to the private corpus as
`gurfil-kasdin-2002-niching-genetic-algorithms-geocentric-orbits-3d-er3bp-cmame-191-doi-10.1016-S0045-7825(02)00481-4.pdf`.

## 1. Header

- **Title (verbatim)**: *Niching genetic algorithms-based characterization of
  geocentric orbits in the 3D elliptic restricted three-body problem*
- **Authors**: Pini Gurfil, N. Jeremy Kasdin — Dept. of Mechanical and
  Aerospace Engineering, Princeton University
- **Venue**: *Computer Methods in Applied Mechanics and Engineering* **191**
  (2002) 5683-5706
- **DOI**: 10.1016/S0045-7825(02)00481-4 (from the printed PII
  `S0045-7825(02)00481-4`, p. 5683 — directly sourced from the page, not
  looked up)
- **Received / revised**: 5 February 2002 / 27 June 2002
- **Sponsor**: Jet Propulsion Laboratory / NASA (Acknowledgements, p. 5706)

Note: this is a **different Gurfil work** from the two already in the corpus
(`gurfil-2007-modern-astrodynamics` digest and the `willis-2008-book-review`
of that same book) — those are the *Modern Astrodynamics* textbook (edited by
Gurfil); this is an original 2002 CMAME research paper by Gurfil & Kasdin. No
overlap; confirmed via `grep -ril gurfil` before filing.

## 2. What the paper actually is

A **Sun-Earth Elliptic Restricted Three-Body Problem (ER3BP)** study (Earth's
orbital eccentricity `e = 0.0167` is retained, unlike the CR3BP) that uses a
**Deterministic Crowding niching genetic algorithm** to search for diverse
families of **geocentric** (Earth-centered, not libration-point-centered)
orbits satisfying practical deep-space-observatory mission constraints
(bounded distance from Earth, no Earth collision, 1-5 year practical
stability). This is NOT a libration-point paper — the origin of the rotating
coordinate frame is Earth's center, not a Lagrange point (Fig. 1, p. 5685),
motivated by NGST/TPF/Starlight/SIRTF-class observatory missions (p.
5683-5684) that need to stay a bounded, non-huge distance from Earth while
avoiding Earth's thermal/radiation/occultation environment.

**Method, not just result:** the actual research contribution is the niching
GA search *methodology* for globally characterizing multiple co-existing
orbit families from one optimization run (deterministic-crowding tournament
selection, Eq. pseudocode p. 5687) — the 14 found families (§5, below) are a
worked demonstration of the method on one specific Sun-Earth mission-
constraint formulation, not claimed as exhaustive.

## 3. Model (p. 5685-5687, §2)

Earth-centered rotating-pulsating coordinate frame (x radially outward, y
along Earth's velocity, z normal to ecliptic — **note this is a different
origin convention from the barycentric-rotating-pulsating frame used in most
other ER3BP literature**, chosen specifically "for design of Earth-centric
orbits," p. 5686). Equations of motion (Eq. 1, then normalized Eq. 9-11) use
true anomaly θ as the independent variable (standard ER3BP practice — makes
the equations autonomous in θ, Eq. 3-4 give the Kepler time↔anomaly map).
Full state `x = [x,y,z,x',y',z']`. `μ = 3.0034495182e-6` (normalized Sun-
Earth mass parameter — **cross-check**: this is essentially identical to
Richardson 1980's Table I value `μ = 3.04036e-6` for the *same* Sun-Earth
system with the Moon's mass folded in — the two papers' μ values agree to
~1% despite 22 years and different sourcing, a useful independent
cross-validation point between the two just-filed papers).

## 4. GA search setup (p. 5687-5690, §3-4)

Objective function `x0* = argmax 1/[(rmax-rmin)²+1]` over a hyper-rectangular
initial-condition search domain (Eq. 15-16), subject to `rmin > R_Earth`
(collision avoidance, Eq. 17). **12 distinct optimization sets** (Table 2, p.
5689) vary the search-domain bounds and initial true anomaly θ0 ∈ {0, π},
yielding **14 total orbit families A-N** (3 from set 1, one each from sets
2-12). GA constants (Table 1): population 200, 400 generations, crossover
p=0.999, mutation p=0.001. Stability classified via a **1-year GA-optimized
run then extended to 5 years**: "Practically Stable" (PS) if distance from
Earth never exceeds 0.1 AU (~15M km) over 5 years, else "Practically
Unstable" (PUS) — a pragmatic, non-Lyapunov stability notion specific to this
mission-constraint framing, not a dynamical-systems stability index.

## 5. The 14 families (§5, p. 5690-5703) — sourced numeric table

**Table 3** (p. 5690, initial conditions) and **Table 4** (p. 5702,
quantitative features) give full sourced numerics for a representative orbit
of every family — reproduced here directly from the typeset tables (Claude
vision, not OCR/digitized):

| Family | Type | r_min (km) | r_max (km) | r0 (km) | v0 (km/s) | Stability | Notes |
|---|---|---|---|---|---|---|---|
| A | DRO | 5,769,577 | 11,740,892 | 5,825,955 | 2.3557 | PS | simple quasi-elliptic DRO, ẏ0≈-2x0, 2 xz-crossings/yr |
| B | DRO | 2,191,130 | 3,270,091 | 2,218,402 | 1.0718 | PS | 4 xz-crossings/yr |
| C | DRO | 986,996 | 1,126,563 | 1,017,578 | 0.8792 | PS | quasi-circular, ~1M km, 8 xz-crossings/yr |
| D | DPO | 317,151 | 1,002,197 | 1,002,197 | 0.3252 | PS | 10 xz-crossings/yr |
| E | DPO | 311,893 | 1,002,197 | 1,002,197 | 0.3145 | PS | 8 xz-crossings/yr, resembles halo/quasi-halo |
| F | DRO | 224,900 | 1,002,197 | 1,000,000 | 0.6865 | PS | along+cross-track librating ellipse, 12 xz-cr/yr, closest-approach DRO (224,900 km) |
| G | ERO | 11,532 | 1,084,448 | 1,000,000 | 0 | PS | quasi-periodic Earth flybys, min approach 11,532 km (~5154 km altitude), 16 xz-cr/yr |
| H | ERO | 15,883 | 750,000 | 750,000 | 0 | PS | **periodic** Earth approach every 25 days, 16 xz-crossings/yr |
| I | DEO | 2,328,213 | 3,962,617 | 3,177,032 | 0.8503 | PUS | bounded ≥1 yr, escapes into heliocentric at 1.233 yr |
| J | 3D DRO | 6,892,060 | 8,510,975 | 7,527,807 | 2.0339 | PS | ẏ0≈-2x0 (3D analog), 2 xz/2xy/2yz crossings/yr, max out-of-ecliptic 0.04 AU |
| K | 3D DRO | 854,531 | 939,889 | 883,956 | 0.8838 | PS | quasi-circular xy-projection, out-of-ecliptic ~145,000 km, <10% rmin/rmax spread |
| L | 3D DEO | 755,223 | 1,662,761 | 1,638,694 | 0.25619 | PUS | PS for 2.686 yr, then escapes heliocentric |
| M | 3D ERO | 17,282 | 1,035,854 | 1,000,221 | 0 | PS | closest approach 17,282 km (10,904 km altitude), Hill-stable, 16 xz/17 xy/17 yz crossings/yr |
| N | 3D DEO | 409,758 | 1,521,047 | 1,521,047 | 0 | PUS | PS for 3.58 yr, then escapes heliocentric |

Family taxonomy (p. 5691): **DRO** = distant retrograde orbit (`x0>0, y'0<0`),
**DPO** = distant prograde orbit (`x0>0, y'0>0`), **ERO** = "Earth-return
orbit" — `v0=0` initial condition, characterized by a close Earth flyby that
increases the spacecraft's inertial velocity (used as a DRO/DPO transfer
trajectory), **DEO** = "delayed escape orbit" — bounded near Earth for
≥1 year, then escapes to heliocentric motion.

## 6. Disturbance robustness (§6, p. 5703-5705)

Lunar third-body gravity (Eq. 19-20, idealized circular 384,320 km orbit,
inclination neglected) and SRP (Eq. 21-23, flat-plate `κ=5e-8 m/s²`
disturbance model) were added to Families A-D and re-integrated. Fig. 18
shows the disturbed vs. nominal trajectories are visually near-identical over
5 years — **the paper's own conclusion (p. 5705) is that these orbits do NOT
drift considerably under realistic perturbations and require no
station-keeping**, in explicit contrast to the well-known high sensitivity of
Sun-Earth L1/L2 halo orbits to the same disturbances. This robustness claim
is qualitative (visual figure comparison only, no quantitative drift number
given) — flagged as **not independently reproducible from this paper alone**
without re-running the actual disturbed integration.

## 7. Relevance to the cyclerfinder codebase (context for the Fable review)

- The codebase already has a substantial, actively-developed ER3BP capability
  axis (`core/er3bp.py`, `core/er3bp_paper_frame.py`,
  `search/er3bp_{isolated_seeds,direct_seeding,discovery,floquet,periodic}.py`,
  `genome/er3bp_{branching,continuation,periodic}.py` — per the project's own
  #286 Track-A axis inventory). **Checked directly (2026-07-12): NEITHER
  existing frame matches Gurfil-Kasdin's.** `core/er3bp.py` uses the Szebehely
  Ch. 10 pulsating/Nechvile frame (independent variable = true anomaly, both
  primaries fixed) — this actually matches Gurfil-Kasdin's *independent
  variable* choice (both use true anomaly θ, Eq. 3-4 vs. their own Nechvile
  derivation) but NOT its *origin* (Szebehely/Nechvile is barycentric,
  Gurfil-Kasdin's Eq. 9-11 is Earth-centered). `core/er3bp_paper_frame.py`
  matches Antoniadou & Libert (2018)'s NON-pulsating barycentric frame
  (independent variable = time, not true anomaly) — different on both axes.
  So Gurfil-Kasdin's specific Earth-centered pulsating formulation is not
  currently implemented anywhere in the codebase; it would need its own
  frame-transform module (a geocentric offset applied to the existing
  Nechvile/true-anomaly machinery) if ever adopted, not a drop-in reuse of
  either existing module.
- **Families G and H are the most cyclerfinder-relevant result in the whole
  paper**: Family H is a genuinely *periodic* (not just quasi-periodic)
  heliocentric orbit with a periodic close Earth flyby every 25 days, PS for
  5 years, remaining bounded the entire time. This is single-body (Earth
  only) rather than a multi-body cycler in this project's usual sense, but it
  is structurally the "Earth half" of a repeated-encounter free-return
  trajectory in a *fully eccentric, single-primary-plus-Sun* dynamical model
  — a candidate reproduction/comparison target for the project's ER3BP
  periodic-orbit search machinery (`search/er3bp_periodic.py`,
  `genome/er3bp_periodic.py`) if that machinery's own positive controls don't
  already cover this specific "periodic Earth-return ER3BP orbit" case.
- The μ cross-check with Richardson 1980 (§3 above) is a small, free,
  independent validation of the project's own `PRIMARIES["Sun"]`/Earth mass
  constants if useful.

## 8. References cited (p. 5706)

Key ones not already flagged above: Szebehely (1967, already in corpus, ref
[6]); Howell, Barden, Lo (1997), "Application of dynamical systems theory to
trajectory design for a libration point mission," J. Astronaut. Sci. 45,
161-178 (ref [8], NOT yet checked against the corpus); Ocampo & Rosborough
(1993), "Transfer trajectories for distant retrograde orbits around the
Moon," AAS Paper 93-180 (ref [10] — DRO terminology precedent, relevant given
the project's own DRO-adjacent work); Henon (1969), "Numerical exploration of
the restricted problem. V. Hill's case," Astron. Astrophys. 1, 223-238 (ref
[9] — the `ẏ0 ≈ -2x0` simple-DRO condition cited in this paper as originally
Henon's, in Hill's simplified problem, not the ER3BP itself). None of these
were checked against the corpus/`KNOWN_CORPUS` as part of this digest pass —
flagging as candidates if the Fable review or a follow-up decides they're
load-bearing.

End of digest.
