# #659 — Fable adjudication of #657's two Antiope gate-passing candidates

Date: 2026-07-19. Scope: strategic/literature/design-risk lens (parallel to an independent
Opus numerical pass). Read-only on catalogue/source; numbers below computed from the
committed code + the #657-logged (C, x0, T) values only.

## Verdict: NO-GO for catalogue admission as Antiope orbits. Both candidates are
## physically impossible at the real (90) Antiope — they pass ~30-40 km BELOW the
## surfaces of both components.

### The decisive check (surface clearance — a gate the machinery does not have)

Reconstructing both orbits from their logged (C, x0) via the Jacobi relation and
integrating one period (DOP853, rtol=atol=1e-12; closure |dr| = 1.3e-5 / 2.6e-6 nd —
the orbits themselves are genuine, exactly as #657 and the coordinator found):

| candidate | min dist to A center | min dist to B center | component radii |
|---|---|---|---|
| (1,1) | 0.06712 nd = **11.8 km** | 0.02370 nd = **4.2 km** | ~43-45 km each |
| (2,2) | 0.07814 nd = **13.8 km** | 0.01796 nd = **3.2 km** | ~43-45 km each |

Antiope's components are enormous relative to their 176-km separation: Descamps et al.
(2007, Icarus 187, 482) Roche-ellipsoid semi-axes 46.5x43.5x41.8 km and 44.7x41.4x39.8 km
(volume-equivalent radius ~44 km from the Aljbaae/Bartczak masses at rho~1.25-1.28 g/cc);
R/d ~ 0.25 in CR3BP units — vs Pluto 0.061 and Charon 0.031 for the PC(3,2) precedent.
Both candidates' closest approaches are at ~7-30% of the body radius — inside solid rock
by a factor of 3-14, not a marginal graze. Even ignoring collision, point-mass gravity is
meaningless at 3-14 km from the center of a 44-km triaxial ellipsoid; the field there is
dominated by the extended non-spherical mass distribution (exactly why the published
Antiope/doubly-synchronous dynamics literature uses polyhedra/ellipsoids, below).

Why every gate still passed: the gate profile (winding topology, prograde, x_max > L1,
Barden |nu|<1, independent Radau crosscheck) is purely dynamical. `RealBinarySystem`
carries no component radii and no clearance gate exists anywhere in
`real_binary_kk_sweep.py` / `pluto_charon_kk_sweep.py` / `binary_star_search.py`. This
never mattered before because Pluto-Charon's bodies are small in nd units. Same-family
failure mode as the #480 EGGIE lesson ("analytic match without per-encounter physical
self-consistency") — see `[[feedback_constructed_tour_per_encounter_self_consistency]]`.

**Fairness check on the catalogued precedent**: PC(3,2) (`ross-rt-pc-cycler-32-2026`)
run through the identical test clears both surfaces — min r1 = 3836 km (Pluto R=1188,
x3.2), min r2 = 1087 km (Charon R=606, x1.8). The existing row is safe; a retroactive
clearance gate would not disturb it (though x1.8 at Charon is worth recording).

### Literature (broader than #657's cycler-keyword check)

The narrow "not-found" is real but misleading — the surrounding field is well-populated,
and the gap exists for a physical reason:

- **Shang, Wu & Cui (2015), Ap&SS 355, 69** — "Periodic orbits in the doubly synchronous
  binary asteroid systems and their applications in space missions": systematic
  grid+differential-correction search of periodic-orbit families in doubly-synchronous
  binaries (809 Lundia, 3169 Ostro; 30 and 28 families), modeling both bodies as Roche
  triaxial ellipsoids (Ivory's theorem). This is exactly the Antiope object class
  (Antiope is the archetype doubly-synchronous system), done with the physically
  required non-point-mass gravity. Closest-prior-art for any future Antiope-class claim.
- **Bellerose & Scheeres (2008)** RF3BP applied to 1999 KW4, and successors (Shi et al.
  2018 CeMDA 130:32 polyhedral KW4 periodic orbits; CeMDA 2024 circular restricted FULL
  3BP with rigid-body spacecraft) — the field standard for binary-asteroid periodic
  orbits is the restricted FULL three-body problem, not point-mass CR3BP.
- **Santos, Sousa-Silva, Terra et al. (2023), arXiv:2307.09657 / P&SS** — spacecraft
  dynamics near an equal-mass binary (2017 YE5) via point-mass CR3BP + Poincare
  sections; prograde+retrograde periodic orbits around each primary. Point-mass CR3BP
  at mu~0.5 motivated by a binary asteroid IS published territory — for a system whose
  components are small relative to separation, unlike Antiope.
- **Aljbaae et al. (2020), MNRAS 496, 1645** (the mu source): polyhedral-gravity
  test-particle stability zones around Antiope specifically — again non-point-mass,
  again no periodic-orbit families. Consistent with: nobody publishes point-mass CR3BP
  periodic orbits grazing Antiope's components because the model is invalid there.
- Abstract-mu context: mu=0.5 is the Copenhagen problem — a century of periodic-orbit
  cataloguing (Stromgren onward). Any "abstract near-mu=0.5 family" claim faces that
  literature, but that claim isn't being made here: Ross-RT 2026 own the (1,1) family
  and our candidate is a 0.8% mu-continuation of the already-catalogued
  `ross-rt-mu05-cycler-11-2026` anchor.

### Distinctness of (1,1) vs (2,2)

Related, almost certainly parent-and-child, not two independent families:
T(2,2)/T(1,1) = 2.0399, delta-C = 0.0199, delta-x0 = 0.011. Not a literal re-traversal
(that would be identical C/x0 at exactly 2T), but the near-2:1 period at nearby (C, x0)
is the signature of a period-doubled branch of the (1,1) family (doubling bifurcations
emanate where the family's nu crosses the stability boundary; a doubled (1,1) has (2,2)
winding by construction). Note this genome already has an iterate-labeling subtlety on
record: the published mu=0.5 (1,1) T is the 3rd iterate of the fundamental. For novelty
accounting these should be treated as ONE family line. Moot for admission given the
surface-collision finding, but relevant if the abstract-mu question is ever pursued.

### Base-rate read (task item 3)

Even before the clearance check, the prior was weak. PC(3,2) — the direct precedent this
candidate's admission pattern mirrors — was catalogued as `our_status:
known-class-member` of the published Ross-RT family, NOT novel. The Antiope (1,1) is a
0.8% mu-continuation of a catalogued anchor row — an even more incremental instantiation
than PC(3,2)'s 8.8% continuation. Best case it was ever going to be a second
known-class-member row, not a #312-class novel family. Against the project's honest
record (one confirmed novel family in the entire program; Galilean 0/36, Titan-Iapetus,
Saturn mid-moons, #549's 32/32, #656's 0/9 all resolving negative or known-class), a
single lightly-searched not-found on an obscure system was never sufficient — and the
targeted search above shows the adjacent field is in fact well-populated.

### Recommended next steps (ranked by value/cost)

1. **Add a surface-clearance gate to the binary-genome machinery** (HIGH value, LOW cost).
   Add sourced component radii to `RealBinarySystem` (+ Pluto/Charon), assert
   min-approach > radius (with a margin) as a 6th gate, and record min-approach distances
   in `SweepResult`. Re-run #657's Antiope cell to convert these 2 candidates into
   documented physical-infeasibility negatives; re-confirm PC(3,2)'s x3.2/x1.8 clearance
   in a test. This also future-proofs every subsequent real-binary sweep — all remaining
   plausible targets (small-body binaries) have large R/d.
2. **Stamp Antiope appropriately in `empty_regions.jsonl`** (LOW cost) — after (1), as
   "point-mass solutions exist but violate physical-body clearance; no physically
   admissible (k1,k2) cycler under the point-mass model," method-conditional per
   `[[project_negative_results_registry]]`. Do NOT stamp it as a plain clean negative —
   the honest statement is model-invalidity, not emptiness.
3. **Only if Antiope-class systems stay strategically interesting**: the physically
   meaningful follow-on is an ellipsoid/RF3BP-gravity variant of the genome
   (Shang-2015-style Ivory potential). MEDIUM-HIGH cost (new dynamics + new gates), and
   the novelty bar is now Shang 2015 + the KW4/RF3BP line, not a blank page. Park unless
   a forcing function appears.
4. **Do not** spend a V0-V2 gauntlet or real-ephemeris run on these two candidates as
   Antiope orbits — validation effort on a physically impossible orbit is waste. The
   abstract mu=0.4961 solutions need no row either: they are ordinary continuation
   members 0.8% from the existing mu=0.5 anchor.

No catalogue writeback recommended, none performed. `data/catalogue.yaml` untouched.
