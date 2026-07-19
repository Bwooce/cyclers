# #659 — Independent adjudication of the #657 Antiope (1,1)/(2,2) candidates

Role: rigorous numerical-methods / physical-realizability judgment. Read-only
(no catalogue or source writes). Parallel to a separate "Fable" pass.

## Verdict (short)

**Do NOT admit.** The two candidates are genuine, reproducible, linearly-stable,
correctly-classified periodic orbits of the *idealized point-mass* CR3BP at
Antiope's mass ratio — the machinery is sound, no ghost minimum. But they are
**not physically realizable cyclers of the real Antiope system**: both orbits
pass ~30–38 km *below the surface* of both asteroids. Antiope is a near-contact,
doubly-synchronous, near-equal-mass binary, and the point-mass CR3BP is an
invalid model for close orbits in exactly this system class. Novelty is
independently undermined: mu≈0.5 is the century-old **Copenhagen problem**, and
there is recent, directly on-topic periodic-orbit literature for close
near-equal non-spherical binary asteroids that the #657 keyword check missed.

## 1. Reproduction (bit-for-bit, current committed code)

- (1,1) via `sweep_family(antiope,"mu05_11")`: C=3.4859955526687223,
  x0=-0.5632517806525419, T=2.946983620993779, nu=5.09e-10, topo (1,1),
  prograde, xcheck dj=8.88e-15. **Matches #657 exactly.**
- (2,2) via `sweep_family_grid(...)`: C=3.4661023165370235,
  x0=-0.5742744462570041, T=6.011499192614617, nu=-7.19e-10, topo (2,2),
  prograde, xcheck dj=1.72e-12. **Matches #657 exactly.**

## 2. Stability / periodicity are genuine (not a ghost minimum)

- Independent Radau re-propagation closes: |s_f - s_0| = 1.8e-10 (both).
- Full 6×6 monodromy eigenvalues all on the unit circle: (1,1) |eig| in
  [0.99963, 1.00037]; (2,2) all = 1.0. Barden nu ≈ 0. Linearly stable.
- 50-period propagation of an off-orbit y+1e-6 perturbation stays bounded
  (max r ≈ |x0|, no escape). The corrector, gates, and stability call are
  all behaving correctly; #620-style ghost minima are not present here.

## 3. DECISIVE: the orbits are deeply subsurface (both bodies pierced)

Antiope components: D_A≈87.8 km, D_B≈83.8 km → radii ≈ 43.9 / 41.9 km;
separation a = 176 km (l_km). So each body radius ≈ **0.24–0.25 nd** — the
bodies fill roughly a quarter of the separation each (near-contact).

| orbit | min dist to PRIMARY centre | to SECONDARY centre |
|-------|----------------------------|---------------------|
| (1,1) | 0.0671 nd = 11.8 km        | 0.0237 nd = 4.2 km  |
| (2,2) | 0.0781 nd = 13.8 km        | 0.0180 nd = 3.2 km  |

Closest approach is **~30 km (primary) and ~38 km (secondary) below each
surface**. The trajectories thread through the interiors of both asteroids —
collision trajectories, physically impossible around the real Antiope.

Contrast the admitted precedent **PC(3,2)** (`ross-rt-pc-cycler-32-2026`),
recomputed here: closest approach 3836 km from Pluto centre (radius 1188 km →
**2647 km above surface**) and 1087 km from Charon centre (radius 606 km →
**481 km above surface**). PC's r/separation ≈ 0.03–0.06 (near-point-mass);
its cycler is a genuine *external* flyby of both bodies. Antiope's r/separation
≈ 0.24–0.25 breaks the point-mass idealization. The "same admission pattern as
PC(3,2)" framing does not hold: PC(3,2) is physically admissible; Antiope's are
not.

## 4. `reaches_secondary` gate is uninformative here

`reaches_secondary = (x.max() > L1)`. At mu=0.4961, L1 = 0.0055 (≈ barycentre),
so the check is "does the orbit reach positive x" — nearly trivially true for
any symmetric orbit and carrying almost no information at mu≈0.5. It does NOT
verify body clearance or a genuine external close-encounter geometry
(cf. [[feedback_constructed_tour_per_encounter_self_consistency]] — a
V∞/geometry-only pass is not sufficient; here it is weaker still). A physical
minimum-clearance-vs-body-radius gate would have flagged these immediately.

## 5. (1,1) vs (2,2): distinct members, same family

Period ratio T22/T11 = 2.0399 (not exactly 2), different C and x0, windings
(1,1) vs (2,2). So NOT one orbit doubled — genuinely two distinct members, but
of the *same* abstract Ross-RT mu=0.5 / Copenhagen family continued to Antiope's
mu. Not two independent discoveries.

## 6. Novelty — undermined on two independent grounds

1. **Copenhagen problem.** mu=0.4961 is 0.8% off the equal-mass CR3BP, the most
   exhaustively studied mass ratio in celestial mechanics (Strömgren's
   Copenhagen school 1913–39; Hénon 1965; Papadakis 1996) — 22 classically named
   periodic-orbit families (a,b,c,f,g,…,z). A symmetric perpendicular-crossing
   periodic orbit at mu≈0.5 is the classical inventory, not a new species.
2. **Directly on-topic recent literature the #657 keyword search missed:**
   - Bakker & Freeman, *Relative Equilibria and Periodic Orbits in a Binary
     Asteroid Model* (arXiv:2306.00273) — explicitly the equal-mass case, 10
     one-parameter families of periodic orbits.
   - *Planar dynamics of non-spherical close tidally locked binaries and the
     restricted three body problem*, Nonlinear Dynamics (2025),
     10.1007/s11071-025-11014-5 — "almost equal masses… two non-spherical
     asteroids tidally locked close to each other," CR3BP **+ ellipsoidal
     geometry**. Literally Antiope's system class, and its whole point is that
     close near-equal binaries require the non-spherical extended-body model —
     independent confirmation of the §3 subsurface objection.
   - Even Aljbaae et al. (2020), the paper #657 already cited, maps stable
     direct/retrograde and internal-resonant orbits in this system — the
     dynamical environment is more mapped than the #657 writeup credited.

   #657's "not-found on both backends" is a keyword-coverage artifact (searched
   "cycler"/RRT genome terms), not evidence of novelty —
   necessary-not-sufficient, as the rule says.

## 7. V-tier

Existing evidence is V0 (single idealized CR3BP model, no real-ephemeris, no
perturbation model). But this is not a "needs more validation to climb" case:
the candidates fail a **physical-admissibility precondition beneath V0** — the
trajectory intersects the physical bodies. The correct model for the cycler
question here is the non-spherical extended-body binary-asteroid RTBP (the 2025
Nonlinear Dynamics approach); in that model these specific point-mass orbits do
not exist (they are collisions). No modest additional work lifts them.

## Recommended next step

Close the Antiope line as a **clean negative for a physically-realizable
cycler** (not an empty-region stamp of the abstract family — the point-mass
orbits genuinely exist; they are just physically inadmissible). Record in
`empty_regions.jsonl`/OUTSTANDING that Antiope's near-contact geometry
(r_body/sep ≈ 0.25) disqualifies point-mass CR3BP cyclers, with the two
missed references. Optional, low priority: add a physical minimum-clearance
(vs body radius) check to `_build_result`/the gate set so future near-contact
binaries self-flag. Do NOT write `data/catalogue.yaml`.
