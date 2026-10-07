# #998: Pluto-Charon one-working-node cyclers with passive small moons (2026-10-07)

Status: enumeration, gauntlet and real-ephemeris check DONE (2026-10-08). The literature step is
DEFERRED (the gate is being reworked under #972); no literature_check module was imported or run.
Nothing here is called novel and no catalogue row is written.

Driver: `scripts/run_998_enumerate.py`. Data: `data/998_pluto_smallmoons/`.

## 1. Model

Circular, coplanar, ideal model. Primary mass: Pluto+Charon system GM 975.5 km^3/s^2 (registry
`PRIMARIES["Pluto"]`). Each body rides a circle whose radius follows from its period by Kepler III
(the Russell-Strange convention). Charon is the only massive flyby body (GM 106.1, R 606 km, floor
100 km, registry). The small moon is a massless target.

Periods are fitted to `plu060.bsp` (2030-01-01 plus 3 yr, 3000 samples, J2000, relative to the
Pluto-system barycentre NAIF 9), not taken from the registry:

| Body | Fitted period (d) | Ratio to Charon (fit) | Kepler a (km) | Registry a (km) | Registry ratio |
|---|---|---|---|---|---|
| Charon | 6.38722 | 1 | 19596 | 19600 | 1.000 |
| Styx | 20.16195 | 3.1566 | 42168 | not in registry | - |
| Nix | 24.85472 | 3.8913 | 48481 | 49300 | 3.989 |
| Kerberos | 32.16798 | 5.0363 | 57577 | not in registry | - |
| Hydra | 38.20202 | 5.9810 | 64569 | 65200 | 6.067 |

The brief's ratios (3.16, 3.89, 5.04, 5.98) agree with the fit. The registry Nix and Hydra values
are 1.7 and 1.0 percent off in a (papercut filed). The registry is not edited.

The four small moons are coplanar with Charon's orbit: the angle between each moon's mean orbit
normal and Charon's is 0.006 (Styx), 0.010 (Nix), 0.348 (Kerberos), 0.260 (Hydra) degrees (same
fit). The coplanar model is therefore a good fit to the geometry.

Known model offset: Charon is at 19596 km in the ideal model (Kepler a of the relative orbit) but
at a mean 17464 km from the barycentre in the kernel. Its speed in the ideal model is therefore
about 12 percent above its barycentric speed. Wobble of Pluto is ignored. The real-ephemeris step
(sec. 6) is where this shows.

## 2. Controls

Two controls. Neither is a Charon-flyby control; the #320 result cannot be recalled by the #998
generator (reasons below). Lead ruling 2026-10-07: accepted. The #320 re-run is a self-regression of the
old pipeline; VenMar#45 is the only control of the generator path; no published Charon-flyby cycler
exists to serve as a positive control. Every #998 result is therefore generator-level and
real-ephemeris-decided: the ideal model is a SEED generator only, and the plu060 step (sec. 6) is the
arbiter for every gate-passing candidate.

1. #320 self-regression. `run_998_enumerate.py control320` re-runs #320's own Pluto sweep through
   its own code (`scan_320_epoch_aware_moon_systems._per_system_sweep`; the module-level import of
   `literature_check` is satisfied by an inert stand-in, so the anchor-overlap field is not
   reproduced). 90 cells evaluated. All 51 committed near-miss rows are reproduced: residual and
   per-encounter V_inf agree to 4.3e-14 km/s, the physical-gate flag and SILVER flag agree on every
   row, no extra rows. The two SILVERs recur: Hydra-Nix-Hydra (1,1) residual 0.00137 km/s and
   Nix-Hydra-Nix (1,1) 0.00070 km/s. Data: `control_320_pluto_rerun.jsonl`,
   `control_320_pluto_rerun.compare.json`.
   Why this is a self-regression only and not a recall by the #998 pipeline:
   - neither SILVER contains Charon; Hydra and Nix are both flyby bodies with GM about 0.002 at
     V_inf of 15-30 m/s;
   - #320 closes V_inf magnitudes plus an anchor wrap, while the generator needs an exact periodic
     closure over k synodic periods;
   - the two legs sum to 94.3 d, and the registry Nix-Hydra synodic period is 74.4 d (ratio 1.267,
     not an integer).
2. Generator code path. `run_998_enumerate.py controlgen` solves VenMar#45 (R-S 2007, Venus flyby,
   Mars massless, key `k2|LV>M/0s|LM>V/0s`) through the same `solve_structure` and `assess`. Exact
   zero (residual 7.1e-15 km/s), gate pass, V_inf Venus 8.2196 and Mars 12.9648 km/s; published
   8.22 and 12.96 (R-S 2007 Table 3), stored #942 values identical to 1e-9. Data:
   `control_generator_venmar45.json`.

## 3. Enumeration

Generator: `search/two_working_body*.py` (the #942/#943 generator), driver `scripts/run_998_enumerate.py`
(wraps `scripts/run_942_enumerate.py`; cells `ps pn pk ph`, Charon body A with returns, the small
moon body B massless). Cycle template [Charon-block, Charon->moon, moon (pass-through), moon->Charon],
one visit per cycle, period k synodic periods of the pair.

Method scope (a zero count holds only under these settings): Charon block of up to 2 returns
(full-rev n:m, half-rev n-pi, same-body generic Lambert returns with 1 or 2 revolutions), transfer
legs with 0 or 1 revolution, k = 1..6 (period up to 56 d, 8.8 Charon periods), seeds n_phase 24,
n_split 10, n_refine 60, exact zero when max residual < 1e-8 km/s. Seed sensitivity: 121 structures
re-solved at the 36/12/40 production seeds of the Jovian cells gave identical zero counts and gate
statuses (k = 1-3 only). Gate: #888/#937 demanded-turn gate at the registry floor (Charon 100 km),
tri-state with the tidal band active (`tidal_speed_kms` is set for Charon), near-180 rule, re-propagation
miss < 1 km.

| Moon | Period ratio to Charon | Structures searched | Structure errors | Exact zeros | pass / indeterminate / fail / no-directions | Merged cyclers | Gate-passing | Strong |
|---|---|---|---|---|---|---|---|---|
| Styx | 3.157 | 7104 | 0 | 11414 | 836 / 6633 / 1664 / 2281 | 2381 | 183 | 108 |
| Nix | 3.891 | 6908 | 0 | 6503 | 637 / 3461 / 1015 / 1390 | 1378 | 149 | 46 |
| Kerberos | 5.036 | 6532 | 0 | 3111 | 218 / 1607 / 618 / 668 | 759 | 62 | 41 |
| Hydra | 5.981 | 6287 | 0 | 1837 | 208 / 926 / 306 / 397 | 493 | 63 | 29 |

"Merged cyclers": zeros with the same cyclic sequence of massive flybys (body, V_inf to 1 m/s, turn
to 0.1 deg), mirror twins and split labels merged. "Gate-passing": a merged group with a member whose
gate status is pass and whose independent re-propagation miss is under 1 km. Zero assessment errors: 0
in the final data. (A first Kerberos run had 3111 assessment errors, KeyError in
`sphere_of_influence_km`, because Styx and Kerberos are not in the satellite registry; the run was
discarded and redone after the driver added them to the in-process SATELLITES dict. Nix and Hydra use
their registry entries, whose sma is 1.0-1.7 percent off; this touches only the SOI diagnostic.)

Gate-passing by k (strong in brackets):

| Moon | k=2 | k=3 | k=4 | k=5 | k=6 |
|---|---|---|---|---|---|
| Styx | 1 (0) | 8 (6) | 40 (25) | 44 (32) | 90 (45) |
| Nix | 0 | 7 (4) | 17 (12) | 50 (22) | 75 (8) |
| Kerberos | 0 | 2 (1) | 6 (0) | 30 (16) | 24 (24) |
| Hydra | 0 | 1 (1) | 8 (4) | 16 (7) | 38 (17) |

No gate-passing cycler has k = 1; at k = 2 only one (Styx).

The two filters beyond the project gate, added in this task's gauntlet (`scripts/gauntlet_998.py`):
- Pluto surface: r_min of the whole cycle (recomputed with `leg_extent`) must exceed Pluto's radius
  1188.3 km. This removes 23 / 53 / 9 / 16 of the gate-passing cyclers (Styx / Nix / Kerberos /
  Hydra). Pluto's atmosphere and the neglected 2100 km barycentric wobble of Pluto are not modelled;
  r_min below about 4000 km (under 3.4 Pluto radii) is not credible in this patched-conic model.
  Strong cyclers with r_min above 4000 km: 29 / 23 / 24 / 18.
- Charon patched-conic validity: the periapsis radius the demanded turn requires (required altitude
  plus 606 km) must lie inside Charon's Laplace sphere of influence (8070 km). A small demanded turn
  at V_inf of 0.15 km/s needs a periapsis of tens of thousands of km, outside the SOI, where Charon
  cannot deliver it as a two-body hyperbola; the project gate passes these because the ratio is
  small. This removes 69 / 72 / 20 / 28 more. The 2-body limit is my choice of filter, not a project
  rule; both counts are reported.

"Strong" = gate-passing, Pluto-surface clear and every Charon periapsis inside the SOI:
Styx 108, Nix 46, Kerberos 41, Hydra 29 (224 in all).

Strong set, ideal model (ranges over the set):

| Moon | Period (d) | V_inf Charon (km/s) | V_inf moon (km/s) | Max demanded turn (deg) | Demanded/available at floor | r_min (km) | r_max (km) |
|---|---|---|---|---|---|---|---|
| Styx | 28.0-56.1 | 0.162-0.305 | 0.079-0.198 | 42.2-67.1 | 0.479-0.743 | 1256-14713 | 42364-135316 |
| Nix | 25.8-51.6 | 0.173-0.342 | 0.070-0.191 | 35.6-67.3 | 0.359-0.812 | 1867-15251 | 49960-103274 |
| Kerberos | 23.9-47.8 | 0.160-0.274 | 0.077-0.146 | 47.2-65.2 | 0.510-0.757 | 1802-16067 | 63203-115095 |
| Hydra | 23.0-46.0 | 0.242-0.359 | 0.095-0.161 | 36.7-66.8 | 0.538-0.767 | 1658-15251 | 65690-118180 |

(Over all gate-passing cyclers the lowest demanded/available ratio is 0.010 for Kerberos, a 0.8 degree
turn that needs a 196,000 km periapsis; those are the cases the SOI filter removes.) Every strong
cycler also passes at H&M's 1.1-radius floor. Almost all strong cyclers (105 / 42 / 40 / 27) have three
Charon flybys per cycle; 3 / 4 / 1 / 2 have two. The moon is met once per cycle, with no turn.

## 4. Independent re-fly, SOI self-consistency

Every gate-passing cycler (457) was re-flown leg by leg with DOP853 (rtol 1e-13, `scripts/gauntlet_942.py`
cross_check, imported with the literature module stubbed out): the integrated gate status agrees with
the generator's on all 457. For the strong set the largest arrival miss at any body is under 0.005 km,
the largest V_inf vector error 3.2e-9 km/s, and the miss as a fraction of the smaller of the Charon
and moon SOI at most 2.2e-6 (Charon SOI 8070 km; the four small-moon SOIs are about 150-300 km). Over all 457 the largest miss is 3.4 km (Kerberos, none in the strong set) and the
largest SOI fraction 1.4e-2. Data: `data/998_pluto_smallmoons/<cell>_gauntlet.json`. A small moon is
met at its centre by construction; whether a body 5-18 km across can be met to the km level is not
addressed.

## 5. Literal-collision checks

- #320 Pluto rows (`data/scan_320_epoch_aware_pluto.jsonl`, rows containing Charon and the moon). Those
  rows are two-leg near-misses (residual 0.0026-0.044 km/s, not closed) and none passes its physical
  gate, so a literal match needs both V_inf values within 0.02 km/s and the flight time within 1
  percent of the period. Result: 0 literal matches for all four moons. A looser NEAR tag (Charon V_inf
  within 0.05 km/s, time within 2 percent) fires for 454 pairings on Nix only, because the 51.0 d total
  of the #320 Nix rows coincides with the k = 6 period (51.6 d); that tag carries no weight.
- Catalogue Pluto rows (`ross-rt-pc-cycler-32-2026` family and the two computed Pluto-Charon
  (3,2) and (5,1) CR3BP rows): CR3BP periodic orbits of a spacecraft in the Pluto-Charon rotating frame.
  They have no V_inf at encounters and no small-moon target, so there is nothing to compare numerically
  against a patched-conic cycler. Not comparable.
- Howett et al. 2021 (Persephone, PSJ 2:75; digest `2026-10-03-digest-howett-2021-persephone.md`): four
  periodic orbits in the Pluto-Charon restricted problem, shown only as figures (Figs. 13-14); the text
  gives no initial conditions, periods or Jacobi constants and no patched-conic or Lambert
  construction, and says encounters with the minor satellites are below 300 m/s. What can be compared:
  the class (repeated close encounters in the Pluto system) and the V_inf scale (this task's V_inf at
  the moons, 0.07-0.2 km/s, is below 300 m/s). What cannot: any orbit, because there is nothing to
  reproduce. This is prior art at the level of a class, not of an object.
- No comparison with Russell-Strange: they did not run Pluto (Jovian and Saturnian sets only).

## 6. Real-ephemeris check (plu060.bsp)

Script `scripts/realeph_998.py`. Top three strong cyclers per moon by lowest demanded/available ratio
(12 in all). The ideal cycle is placed on the kernel states of Charon and the moon (NAIF 901-905
relative to barycentre 9, in Charon's mean orbit plane over 2030-2034) with the central mass
975.5 km^3/s^2 at the barycentre. Epochs: the first conjunction (moon and Charon at the same in-plane
angle) after 2030-01-01 and each +400 d, five epochs (2030-01 to 2034-05). The ideal start dates are
corrected on the real states (`correct_dates`, least squares on the H&M residual), for one cycle (n1)
and for three chained cycles (n3). Check of the check: the same residual evaluated on the ideal
system is 1e-15 km/s, so the script reproduces the ideal cycle before the kernel is introduced.

| Moon | Cycle key | Period (d) | V_inf Charon / moon (km/s) | Max turn (deg) | Ideal ratio | r_min / r_max (km) | Charon rp needed (km) | Defect on ideal dates (m/s) | Residual after refit, 1 cycle (m/s) | After refit, 3 cycles (m/s) | Real-eph ratio, 1 cycle |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Styx | `k4|RCharon/1:2|RCharon/1:2|LCharon>Styx/1h|LStyx>Charon/1h` | 37.40 | 0.253 / 0.112 | 42.6 | 0.479 | 1769 / 42911 | 2897 | 11.8 | 0.131 | 0.492 | 0.563 |
| Styx | `k3|RCharon/1:2|RCharon/2:3|LCharon>Styx/0s|LStyx>Charon/0s` | 28.05 | 0.162 / 0.079 | 56.3 | 0.482 | 3584 / 46012 | 4525 | 7.9 | 0.052 | 0.621 | 0.558 |
| Styx | `k5|RCharon/2:3|RCharon/1:2|LCharon>Styx/1l|LStyx>Charon/1h` | 46.74 | 0.186 / 0.086 | 53.1 | 0.488 | 3366 / 45182 | 3793 | 4.3 | 0.082 | 0.669 | 0.568 |
| Nix | `k4|LCharon>Charon/1l|LCharon>Nix/0s|LNix>Charon/0s` | 34.39 | 0.178 / 0.130 | 40.0 | 0.359 | 2622 / 99663 | 6473 | 31.1 | 0.012 | 0.285 | 0.357 |
| Nix | `k4|RCharon/1:2|RCharon/1:2|LCharon>Nix/1l|LNix>Charon/0s` | 34.39 | 0.173 / 0.070 | 54.3 | 0.481 | 3239 / 49960 | 4217 | 12.3 | 0.019 | 0.686 | 0.567 |
| Nix | `k3|RCharon/1:2|RCharon/1:2|LCharon>Nix/0s|LNix>Charon/0s` | 25.79 | 0.332 / 0.169 | 35.6 | 0.505 | 2218 / 66296 | 2190 | 10.2 | 0.058 | 0.068 | 0.578 |
| Kerberos | `k6|RCharon/1:2|RCharon/1:2|LCharon>Kerberos/1l|LKerberos>Charon/0s` | 47.82 | 0.163 / 0.082 | 59.3 | 0.510 | 3164 / 68138 | 4064 | 5.4 | 0.031 | 0.556 | 0.593 |
| Kerberos | `k6|RCharon/2:1|RCharon/2:1|LCharon>Kerberos/0s|LKerberos>Charon/0s` | 47.82 | 0.175 / 0.103 | 59.9 | 0.533 | 16067 / 83628 | 3466 | 14.0 | 0.044 | 0.340 | 0.541 |
| Kerberos | `k6|RCharon/2:1|RCharon/3:2|LCharon>Kerberos/0s|LKerberos>Charon/0s` | 47.82 | 0.161 / 0.077 | 60.0 | 0.512 | 15602 / 65719 | 4102 | 7.3 | 0.017 | 0.376 | 0.545 |
| Hydra | `k4|RCharon/1:2|RCharon/1:2|LCharon>Hydra/0s|LHydra>Charon/0s` | 30.68 | 0.273 / 0.105 | 45.2 | 0.538 | 2052 / 71853 | 2286 | 6.3 | 0.438 | 0.760 | 0.618 |
| Hydra | `k6|RCharon/1:2|RCharon/1:2|LCharon>Hydra/1h|LHydra>Charon/0s` | 46.02 | 0.271 / 0.106 | 45.6 | 0.540 | 2043 / 73425 | 2285 | 5.2 | 0.417 | 0.861 | 0.619 |
| Hydra | `k5|RCharon/2:3|RCharon/2:3|LCharon>Hydra/0s|LHydra>Charon/0s` | 38.35 | 0.344 / 0.124 | 36.7 | 0.541 | 5419 / 65690 | 1943 | 32.0 | 0.518 | 0.511 | 0.608 |

Read it as follows. The ideal start dates leave a V_inf defect of 4-32 m/s on the real states
(2-17 percent of V_inf at Charon). Re-fitting the dates removes all but 0.01-0.52 m/s for one cycle
and 0.07-0.86 m/s for three chained cycles (worst over the five epochs). The demanded-turn gate
passes at every epoch for all 12 candidates, with the demanded/available ratio about 17 percent higher
than in the ideal model (0.54-0.62 against 0.48-0.54; Nix k4 unchanged at 0.357). The re-fit date
shifts are 1-12 h (Hydra k5: 59 h at one epoch). The residual is not zero, so these are not exact
real-ephemeris cyclers: they close to under 1 m/s (0.5 percent of V_inf) in a patched-conic ephemeris
model, which is the level at which Russell and Strange call a cycler ballistic before a full n-body
optimisation. `converged` at the 1e-6 km/s tolerance was false for all of them.

Limits of this check, stated plainly (the first is a modelling limit of the one-centre patched conic
for a binary with mass ratio 0.11: Charon's ideal radius from the system GM is 19,596 km, its
barycentric radius in the kernel 17,464 km, about 12 percent in speed; a barycentric ideal model for
binaries is registered as its own task, and src is not patched here):
- The model is the barycentre-centred patched conic. Pluto's own 2100 km wobble and the 12 percent
  offset between Charon's ideal radius (19,596 km) and its barycentric radius (17,464 km) enter only
  through the kernel states of the bodies; the spacecraft still feels one point mass at the barycentre.
  The strong cyclers with r_min near 1.5-3 Pluto radii (several of the 12 above) would see the real
  Pluto-Charon field, not this one. The trajectories were not integrated in an n-body model.
- Closure is the H&M magnitude/vector residual, not a periodic-orbit solve; a three-cycle chain is the
  longest tested.
- Only the 12 best-ratio strong cyclers were checked, not the 224.

### 6.1 Real-ephemeris check of every gate-passing candidate (lead ruling)

The same script, run on all 457 gate-passing cyclers (`--pool all`, one cycle, five epochs;
`<cell>_realeph_all_n1.json`). Strong set (224): the demanded-turn gate passes at all five epochs for
219 (Styx 106 of 108, Nix 44 of 46, Kerberos 40 of 41, Hydra 29 of 29), and 217 of those also refit to
under 1 m/s at every epoch (the worst strong residual is 17 m/s, Hydra, in 2 cycles that close
above 1 m/s but below 10). The 233 gate-passing cyclers outside the strong set refit as well
(residual under 1 m/s at every epoch for 212, gate pass at every epoch for 209), so the real-ephemeris
step does not separate them; the strong filters rest on the patched-conic argument in sec. 3, not on
this step. Not done: the three-cycle chain was run only for the 12 best-ratio cyclers.

## 7. What was verified and what was assumed

Verified in this task:
- the #320 Pluto sweep re-runs to 4e-14 km/s (self-regression only);
- VenMar#45 comes out of the generator at the published V_inf;
- the kernel periods and coplanarity (sec. 1);
- every gate-passing cycler re-flown with DOP853;
- the real-ephemeris residuals above.

Assumed or not done:
- the circular coplanar ideal model, with the Charon radius set by Kepler III from the system GM;
- GMs of Styx (0.0007) and Kerberos (0.0011) are approximate and enter only the SOI diagnostic;
- the strong-set filters (Pluto surface, Charon SOI) are mine;
- no screen against crossings of the other small moons' orbits or against Pluto's atmosphere;
- no n-body integration; no 3D (inclined) extension;
- the literature step (deferred to #972): no cycler here is called novel, and "gate-passing" does not
  mean "unpublished".

## 8. Result

This is not a zero. In the ideal circular-coplanar model, with Charon as the only massive node and
Charon returns up to 2 per cycle, the enumeration finds 457 gate-passing one-working-node cyclers
(183 Styx, 149 Nix, 62 Kerberos, 63 Hydra), 224 of them after the two extra filters. Their V_inf at
Charon is 0.16-0.36 km/s, their periods 23-56 d (k = 3-6), and the 12 best close to under 1 m/s on
plu060 states at five epochs. They are candidates pending the literature gate (#972) and adjudication.
They are not catalogue rows and not novel claims. The #320 Hydra-Nix result is not recalled by this
generator (sec. 2), so the new run has no Pluto-specific positive control; its controls are VenMar#45
and the #320 self-regression.

Files: `scripts/run_998_enumerate.py`, `scripts/gauntlet_998.py`, `scripts/realeph_998.py`;
`data/998_pluto_smallmoons/` (control files, `<cell>/structures.jsonl` and `settings.json`,
`<cell>_gauntlet.json`, `<cell>_realeph_n1.json`, `<cell>_realeph_n3.json`). The per-zero records
(`<cell>/zeros.jsonl`, 34 MB) are not committed; they regenerate by re-running the enumeration with the
settings in `settings.json`.
