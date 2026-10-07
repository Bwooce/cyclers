# #997: Earth-Moon both-primary lineage (Newton alpha/beta types, Hoelker-Winston Fig. 90, Schwaniger perigee continuation)

Task `#997` (earthmoon-opus, opened by the owner 2026-10-07 ~21:50 AEDT; the Earth-Moon gate is
open for this task). Source of the task: `docs/notes/2026-10-07-971-fable-corpus-review-2.md`
sec. 5. Control: `#970` (`docs/notes/2026-10-07-970-schwaniger-control.md`).
The literature step is DEFERRED to `#972`. Nothing in this note is called novel. A clean "all
known-class" is a result.

## 0. Pre-registration (written and committed BEFORE any #997 continuation was run)

### 0.1 Model and units

- Planar CR3BP, the project's frame (Earth at (-mu, 0), Moon at (1 - mu, 0)), registry mass ratio
  mu = `cr3bp_system("Earth", "Moon").mu` = 0.01215058439469525, L = 384,400 km, TU =
  375,190.26 s. Seeds printed at another mu (Newton 1/82.45, Hoelker-Winston 1/80) are first
  continued in mu to the registry value at fixed topology, with the same stopping rules.
- Jacobi constant C in the project convention (no mu(1 - mu) term); the table also gives the value
  with the term.
- Floors (the registry): Earth 200 km altitude (`PLANETS["E"].safe_alt_km`, R_E = 6378.137 km,
  r1 >= 6578.137 km); Moon 100 km altitude (`SATELLITES["Moon"].safe_alt_km`, R_M = 1737.4 km,
  r2 >= 1837.4 km). The minimum distances are taken over the WHOLE period (event-refined extrema
  of r1 and r2), not only at the symmetric crossings.

### 0.2 Families

| id | family | seed | topology fixed by |
|---|---|---|---|
| F1 | Schwaniger retrograde cislunar free return | the `#970` member (periselene start, Earth side of the Moon) | target = 2nd y = 0 crossing (perigee) |
| F2 | Newton type 1/2, retrograde about E | Newton 1959 Table 1 rows 1-5 (x0 = 1.00 ... 1.20 at mu = 1/82.45) | target = 3rd y = 0 crossing (Newton's own rule) |
| F3 | Newton type 1/2, direct about E (mixed, zero-velocity and direct subclasses are one curve) | Table 1 rows 6-11 | 3rd crossing |
| F4 | Newton's other alpha/beta types: generating Kepler ellipse with inertial period 2 pi alpha/beta and apocentre at the Moon's distance, direct and retrograde | mu = 0 rotating-Kepler start at x0 = 1 - mu + rho, then corrected at the registry mu | the crossing nearest pi alpha on the seed |
| F5 | Hoelker-Winston 1968 Fig. 90 lemniscate completion (and the Fig. 89 neighbour) | x0 = 1.4875, ydot0 near -0.84835 at mu = 1/80; corrected for xdot = 0 at each crossing index up to 10 | the index found |

F4's alpha/beta set is fixed now, before running. alpha/beta must exceed 2^(-3/2) = 0.354 so that
the generating ellipse reaches the Moon's distance (a > 1/2). The pre-registered set is {2/5,
3/8, 4/11, 3/7, 2/3, 3/4} (perigee of the mu = 0 ellipse: 0.086, 0.040, 0.019, 0.136, 0.526,
0.650 L). This covers the three with perigees below GEO distance (2/5, 3/8, 4/11) and three
low-beta controls. Each is tried in both senses. Nothing outside this set is run under `#997`.

### 0.3 Corrector and continuation

- Symmetric perpendicular-crossing corrector. Unknowns are (x0, ydot0) at a perpendicular
  x-axis crossing. The residual is xdot at the target y = 0 crossing (index fixed per family, as in
  the table), with the crossing time projected out in the usual way
  (`correct_symmetric_fixed_jacobi`'s identity). The STM comes from the variational equations
  (`core.cr3bp.cr3bp_stm_eom`), DOP853, rtol = atol = 1e-12.
- Continuation: pseudo-arclength along the solution curve G(x0, ydot0) = 0 in the (x0, ydot0)
  plane, with the tangent the null vector of [dG/dx0, dG/dydot0]. This passes folds in x0 and in
  C. Step control: initial ds = 2e-3; halve on a Newton failure (more than 8 iterations, or
  |G| > 1e-11 after them); grow x1.5 after 3 easy steps (at most 4 iterations); ds_max = 2e-2,
  ds_min = 1e-6. Each member records C, T, perigee and periselene distances and altitudes, the
  apogee, b_h = tr(M4) - 2 and b_v = tr(Mz) from the full-period monodromy (fixed_path STM), and
  the crossing count.
- Both directions from every seed. Checkpoint: one JSONL per family under `data/997_lineage/`,
  appended and flushed per member; a re-run resumes from the last member.

### 0.4 Stopping rules (per direction)

1. **Impact at the floors.** The direction stops at the first member whose minimum r1 over the
   period is below 6578.137 km, or whose minimum r2 is below 1837.4 km. That member is recorded
   and flagged. Exception, labelled separately: F1's start member is itself 23 km below the Earth
   floor (177 km altitude). Its down-perigee direction is run on to the Earth's SURFACE (r1 =
   6378.137 km), and those members are flagged "below floor". They are for `#948` interest only
   and are never cycler-class candidates.
2. **Period doubling.** b_h crosses -2: the direction stops at the crossing. A b_h crossing +2 (a
   fold or tangent bifurcation) is recorded and the continuation goes on through it.
3. **Loss of symmetry or topology.** The direction stops if the number of y = 0 crossings
   before the target crossing changes, or if the target crossing changes type (for example, the
   perigee crossing stops being the global Earth minimum, or a lunar loop appears or vanishes:
   the M2 rotation count changes). The corrector enforces mirror symmetry, so "loss of symmetry"
   shows up as this kind of topology change, or as a corrector failure at ds_min.
4. **Escape or runaway.** max r1 > 5 L (Restrepo-Russell's own r1max), T > 60 TU, or 400 members
   per direction.
5. **Step failure.** ds < ds_min.

### 0.5 Cycler-class candidate test

A member is a **cycler-class candidate** if it is periodic and ballistic, and EVERY period it has
at least one Earth pass and at least one Moon pass, with:
- every pass clearing both floors (min r1 >= 6578.137 km, min r2 >= 1837.4 km);
- an Earth pass at perigee r1 <= 42,164 km (GEO radius, the catalogue's Vaquero LEO-GEO band
  ceiling, catalogue line 56634; the band's floor is replaced by the registry 200 km floor);
- a Moon pass inside the lunar Hill radius, (mu/3)^(1/3) L = 61,600 km (Hoelker-Winston digest test).

Families are reported whole. Candidates are the members that pass this test. Nothing here is a
novelty claim.

### 0.6 Known-class gate and database exclusions

- **Restrepo & Russell 2018** (digest `docs/notes/2026-06-17-digest-restrepo-russell-2018.md`; read
  p.3, 5, 7-8): planar AXISYMMETRIC periodic orbits only. Excluded: asymmetric orbits, orbits
  centred on L3/L4/L5, and 3-D orbits (p.3). The global search covers -5 x_L1 < x0 < 5 x_L1 in the
  MOON-centred frame (x_L1 = the Moon-L1 distance, 0.1509 L, so |x0 - x_Moon| < 0.754 L), stops
  ydot0 once a trajectory reaches r1 = 5 L, and terminates at the N-th crossing for N <= 10 (p.7-8).
  The resonance families go up to q_max = 10 for Earth-Moon (Table 1). The paper states no
  exclusion for Earth or Moon impact or for close passes; whether their Fortran pipeline kept
  near-collision orbits is unknown. A family of this note is IN their search domain if one
  perpendicular crossing satisfies |x0 - x_Moon| < 0.754 L with N <= 10 crossings to T/2 and max
  r1 <= 5. Membership itself can only be settled against the Earth-Moon data files. `#1019` lists
  only the Jupiter-Ganymede and Saturn-Titan folders, so the Earth-Moon files are NOT held. Until
  they are, the label is "in RR domain, membership unchecked" or "outside RR domain (reason)".
- **Franz & Russell 2022** (`docs/notes/2026-07-28-747-franz-russell-casoliva-crosscheck.md`):
  "Orbits that escape the vicinity of the Moon (defined as ever being more than 4 Hill's units
  from the Moon, approximately 350,000 km ...) or impact the Moon ... are not considered". Any orbit
  with an Earth pass below about 34,000 km from Earth's centre is more than 350,000 km from the
  Moon at that pass. **So every cycler-class candidate in the 0.5 sense is excluded from
  Franz-Russell by construction**, and absence there means nothing.
- **Catalogue Earth-Moon rows** (Casoliva, Vaquero, Ross/Roberts-Tsoukkas, Braik-Ross, Arenstorf,
  Genova-Aldrin, Wittal): compared by (C, T, perigee, periselene, crossing topology). Hits are
  reported as known-class.
- **Newton 1959, Hoelker-Winston 1968, Schwaniger 1963** are themselves the published sources: a
  member continued from their seed is known-class (`known-class-member`) by construction.

### 0.7 Amendments (each recorded before the runs it affects)

- 2026-10-07 22:06 AEDT, after the first F1 run, before the F2-F5 runs. New stop rule: **the start
  crossing and the T/2 crossing swap.** At that point the curve in (x0, ydot0) passes through an
  orbit of period T/2 and then retraces the same orbits from the other symmetric crossing. This is
  a period-halving end point (b_h -> +2 for this family, -2 for the half-period family). The first
  F1 run went past it and retraced 300 members that were already computed; I checked this by
  matching (C, T, perigee) and by showing that member 242's T/2 crossing is member 15's start.
  That run is kept in scratch only.
- 2026-10-07 22:14 AEDT, after the first F4 seed pass, before any F4 continuation or rho < 0.03 run:
  - The seed rho set for F4 is {0.005, 0.01, 0.02, 0.03, 0.06, 0.1, 0.15, 0.2}. The script's
    first pass had used only 0.03-0.2, which cannot reach 4/11 and most of 3/8: they need
    r_a < 2a, so rho < 0.019 (4/11) or rho < 0.040 (3/8).
  - A seed that closes with |t_half - pi alpha| > 0.15 pi alpha is dropped as "not the
    alpha/beta type". In the first pass such seeds were lunar captures with many loops, for
    example 2/5-d at rho = 0.03 closed with t_half = 3.1, not 6.3.

## 1. Results

Runs 2026-10-07 22:01-23:20 AEDT. Script `scripts/run_997_lineage.py` (seeds + continuation),
`scripts/run_997_gate.py` (known-class gate and summary). Data in `data/997_lineage/`: one JSONL
per seed and direction (`<seed>_p` / `_m`, every member with C, T, x0, ydot0, perigee/periselene/
apogee, b_h, b_v, det M4, closure, crossing count, winding numbers about each primary, floor and
candidate flags, Restrepo-Russell and Franz-Russell domain flags), `seeds.json`, `gate.json`,
`summary.json`, and `F5_mu_1_80_closures.json`.

### 1.1 Verdict

- **All cycler-class candidates are known-class.** There are three candidate stretches, on F1, F3
  and F4 3/7-retrograde. All three are continuations of published seeds (Schwaniger 1963, Newton
  1959). The 3/7-retrograde stretch also contains the catalogued Casoliva et al. 2010 7-3a orbit.
  F3 contains the catalogued Vaquero 2:1 rows. Nothing is called novel; the literature step stays
  with `#972`.
- **No new catalogue row is proposed.** Members are known-class-member by construction (sec. 0.6).

| family | distinct members | C range | T range (d) | cycler-class candidates | catalogue rows on the family |
|---|---|---|---|---|---|
| F1 Schwaniger retrograde cislunar | 46 | -0.612 ... 1.094 | 23.1-26.1 | 19: perigee 507-34,561 km alt, periselene 475-1,332 km alt, C 0.582-1.077, T 24.5-26.0 d | none (Schwaniger's row is the `#970` draft) |
| F2 Newton 1/2 retrograde | 146 | -0.489 ... 1.218 | 27.3-34.3 | 0 (low-perigee members have periselene 80,000-91,000 km alt, outside the Hill radius; low-periselene members have perigee > 44,000 km alt) | casoliva-2-1b (on segment B, exact crossing match) |
| F3 Newton 1/2 direct (mixed, zero-velocity, direct) | 241 | 1.928 ... 2.670 | 23.5-27.0 | 10: perigee 27,800-35,730 km alt, periselene 105-57,500 km alt, C 2.240-2.441 | vaquero-21-c198, -c246, -c247, -c266 |
| F4 3/7 retrograde | 38 | 0.900 ... 1.355 | 81.4-82.3 | 38 (all but the floor member): perigee 2,087-28,880 km alt, periselene 14,335-29,325 km alt | casoliva-7-3a (exact crossing match) |
| F4 2/3 retrograde | 76 | -0.637 ... -0.035 | 54.9-61.8 | 0 (perigee > 100,000 km alt) | none |
| F4 3/4 retrograde | 56 | -0.682 ... 0.295 | 78.2-88.5 | 0 (perigee > 44,000 km alt) | none |
| F4 2/5, 3/8, 4/11 and the direct senses | - | - | - | no seed of the type (sec. 1.4) | - |
| F5 Hoelker-Winston Fig. 89/90 | 0 at registry mu | - | - | 0 (sec. 1.5) | - |

### 1.2 Known-class gate

- **Catalogue (all planar Earth-Moon rows with a state, `data/997_lineage/gate.json`).** Each row's
  state was propagated over its period at its own mu. Its perpendicular crossings were compared
  with the family curves, interpolated in C between neighbouring members.
  - Exact matches, where the crossing state agrees to < 1e-3 and T to < 1e-3 relative:
    - casoliva-7-3a-em-cycler-2010 on F4 3/7-r: perpendicular crossing (1.05875, -1.51230)
      against the interpolated (1.05874, -1.51231). Perigee 19,267 against 19,268 km alt;
      periselene 25,515 against 25,515 km alt. b_h = -4.97 in my convention; the row's own
      Casoliva k convention differs.
    - casoliva-2-1b on F2.
    - vaquero-21-c198, -c246, -c247 and -c266 on F3.
  - No match: casoliva-7-3b and -7-3c (not x-axis symmetric, so outside this symmetric lineage);
    Casoliva 1-2c/d/e, 2-1a and 3-2c; Vaquero 3:1 c254 and c313; Ross/Roberts-Tsoukkas;
    Braik-Ross. These are other classes; none of them is a member of these families.
- **Restrepo & Russell 2018.** Every computed member has a perpendicular crossing within 5 x_L1
  of the Moon, N <= 10 to T/2 and max r1 < 5. So every member is **inside the RR global-search
  domain**. Membership is **unchecked**: the Earth-Moon data files are not held (`#1019` covers
  JG/ST only; question sent to the lead). Their paper states no impact filter. F1's 177 km
  perigee member and F2-F3's floor members could have been kept or dropped by their pipeline;
  unknown.
- **Franz & Russell 2022.** Every candidate has an Earth pass far beyond 350,000 km from the Moon,
  so it is **excluded by construction**. Absence from that database means nothing here.
- **Published seeds.** F1 is Schwaniger's own family (he published only the 6555 km member). F2
  and F3 are Newton's type 1/2. F4 3/7-r is Newton's general alpha/beta type ("T_P = 2 pi
  alpha/beta", described but not computed by Newton). It is also Casoliva's 7:3 class and
  Arenstorf's m/k = 3/(-7) second-kind continuation class.

### 1.3 Positive controls and checks

- All 11 Newton 1959 Table 1 rows re-close at Newton's mu = 1/82.45 with my corrector. The
  corrected ydot0 agrees with the printed value to 1e-4 to 5e-4 (`seeds.json`,
  `newton_mu_ydot0`), and the periods to 1e-3. The rows were then continued in mu (4 steps) to
  the registry value at fixed x0.
- Hoelker-Winston Fig. 89 (mu = 1/80): ydot0 = -0.8483581, P/2 = 9.97328, which reproduces the
  digest. Fig. 94 (libration): ydot0 = -0.8478460, P/2 = 4.19518, also reproduced.
- The F1 seed is the `#970` member. F1 continuation passes back through the seed's numbers.
- Every member records det M4 and the full-period closure of the half-period-corrected state.
  Over all F1-F4 members: |det M4 - 1| <= 2.6e-7 (2.7e-8 on the candidates), and full-period
  closure <= 2.7e-9 (nondimensional).

### 1.4 F4 seeds (alpha/beta types)

- Seeds were built from the mu = 0 rotating-Kepler ellipse with apocentre r_a = 1 + rho beyond
  the Moon, then corrected at the registry mu with the crossing index fixed on the mu = 0 path
  (sec. 0.2, 0.7).
- Closed as the type: 3/7-r (rho = 0.06), 2/3-r (rho = 0.02 ... 0.2; three continuation
  segments separated by period doubling) and 3/4-r (rho = 0.1).
- 2/5-d (rho = 0.06) closed onto an orbit whose T/2 crossing equals its start: the F3 type-1/2
  orbit traversed twice (C = 2.6385; F3-9 has it at k = 14-15). It is a duplicate, not a 2/5
  member.
- **4/11 and most of 3/8 cannot exist in this construction.** The generating ellipse with apocentre
  beyond the Moon has perigee 2a - r_a < R_E for rho > 0. So "4/11 with a LEO perigee" needs the
  apocentre inside the Moon, and no 4/11 seed exists. 3/8 and 2/5 with low rho failed to close or
  closed as lunar captures (half period far from pi alpha; dropped by the sec. 0.7 rule).
- The direct senses of 3/7, 2/3 and 3/4 did not close as the type.
- Not tried (outside the pre-registration): starts with the apocentre short of the Moon (rho < 0),
  and other alpha/beta.

### 1.5 F5 (Hoelker-Winston Fig. 90 lemniscate)

- At mu = 1/80 the Fig. 89/90 initial conditions close at crossing indices n = 2-6
  (`F5_mu_1_80_closures.json`).
- The n = 3 closure from Fig. 90 (ydot0 = -0.8483497, P/2 = 14.951, C = 2.85963) has winding
  numbers -1 about the Earth and +1 about the Moon. That is the lemniscate (figure-eight)
  topology the paper predicts.
- It is numerically marginal: b_h about 1.5e9, and the half-period-corrected state misses
  full-period closure by 3e-3.
- It passes 200 km from the Moon's CENTRE, inside the Moon, and its nearest Earth approach is
  247,000 km.
- Every F5 closure passes 49-200 km from the Moon's centre, which is an impact in the physical
  sense. Every one has its Earth pass beyond 220,000 km. So **no F5 orbit can be a cycler-class
  candidate at any nearby mu**.
- The fixed-x0 mu homotopy to the registry mu failed for n >= 2.
- One n = 4 run "succeeded" but changed orbit (t_half went from 24.6 to 10.3; files
  `F5-fig90-n4_*`, labelled in `summary.json` as not an F5 member).
- The Fig. 94 libration orbit, continued to the registry mu, passes inside the Moon. Its
  continuation stops at once (rule 1).

### 1.6 Stops, per family

- F1: up-perigee, the start and T/2 crossings swap at C = -0.612. The orbit approaches a
  period-T/2 orbit (b_h -> 2.00) at perigee 314,000 km alt. Down-perigee, impact with the Earth's
  surface between 177 km and -148 km alt (C 1.0854-1.0939).
- F2: segment A (x0 0.993-1.002, Newton row 1): period doubling at C = -0.488 and lunar-floor
  impact at C = 0.126. Segment B (rows 2-5): period doubling at C = 0.041-0.058 and Earth-surface
  impact at C = 1.218. The flip-unstable stretch between A and B was not continued (rule 2). In
  the plane, segment B is linearly stable (b_h 1.0-1.5); vertically it is just unstable (b_v about 2.04-2.4).
- F3: several period-doubling stops (C = 2.466-2.478 and 2.667); lunar-floor impact at C = 2.227-2.230 (near-Moon branch); Earth-floor impact at C = 1.928.
- F4 3/7-r: period doubling at C = 0.900; Earth-floor impact at C = 1.355 (perigee 189 km alt).
  Flip-unstable throughout (b_h -1.8 to -30.5). Vertically stable for C below about 1.065 (|b_v| < 2).
- F4 2/3-r and 3/4-r: period doubling at both ends of every segment.

### 1.7 Member tables (sampled; every member is in the JSONL)

Altitudes over Earth R = 6378.137 km and Moon R = 1737.4 km; C in the project convention (add
mu(1 - mu) = 0.0120030 for the other convention); b_h = tr(M4) - 2, b_v = tr(Mz) (stable iff in
[-2, 2]); "cand." = the sec. 0.5 test. x0 is the start crossing (Earth-Moon line).

#### F1 (46 distinct members)

| C | T (TU) | T (d) | x0 | perigee alt (km) | periselene alt (km) | b_h | b_v | cand. | note |
|---|---|---|---|---|---|---|---|---|---|
| -0.61230 | 5.5644 | 24.16 | 0.82262 | 314250 | 61777 | 2 | 1.03 |  |  |
| -0.61032 | 5.5623 | 24.15 | 0.80271 | 306854 | 54907 | 2.02 | 1.01 |  | start and T/2 crossings swap (period-halving end point; retrace beyond) |
| -0.58519 | 5.5367 | 24.04 | 0.88258 | 284061 | 38728 | 2.29 | 0.686 |  |  |
| -0.37927 | 5.3835 | 23.38 | 0.94881 | 206861 | 13271 | 8.33 | -5.3 |  |  |
| -0.12833 | 5.3275 | 23.13 | 0.96745 | 144217 | 6104 | 31.6 | -24 |  |  |
| 0.08792 | 5.3689 | 23.31 | 0.97392 | 101704 | 3618 | 73.6 | -50.5 |  |  |
| 0.29994 | 5.4615 | 23.72 | 0.97731 | 68260 | 2314 | 139 | -82.5 |  |  |
| 0.52318 | 5.5929 | 24.29 | 0.97945 | 40677 | 1493 | 231 | -118 |  |  |
| 0.73942 | 5.7402 | 24.93 | 0.98078 | 20514 | 980 | 340 | -149 | yes |  |
| 0.89557 | 5.8549 | 25.42 | 0.98148 | 9593 | 713 | 424 | -167 | yes |  |
| 0.98824 | 5.9256 | 25.73 | 0.98181 | 4478 | 583 | 474 | -175 | yes |  |
| 1.04799 | 5.9720 | 25.93 | 0.98201 | 1705 | 509 | 506 | -180 | yes |  |
| 1.08542 | 6.0015 | 26.06 | 0.98212 | 177 | 465 | 525 | -182 |  |  |
| 1.09392 | 6.0082 | 26.09 | 0.98214 | -148 | 456 | 530 | -183 |  | impact at the physical surface |

#### F2 (146 distinct members)

| C | T (TU) | T (d) | x0 | perigee alt (km) | periselene alt (km) | b_h | b_v | cand. | note |
|---|---|---|---|---|---|---|---|---|---|
| -0.48906 | 7.8695 | 34.17 | 1.00112 | 174164 | 3365 | 2.26 | 42.9 |  |  |
| -0.48823 | 7.8543 | 34.11 | 1.00158 | 175908 | 3541 | -2.48 | 41.5 |  | period doubling (b_h crosses -2) |
| -0.48164 | 7.8929 | 34.27 | 0.99988 | 167302 | 2887 | 21.2 | 46.8 |  |  |
| -0.42256 | 7.8601 | 34.13 | 0.99753 | 145195 | 1986 | 97.7 | 55.2 |  |  |
| -0.07780 | 7.4843 | 32.50 | 0.99358 | 72326 | 467 | 472 | 150 |  |  |
| 0.04113 | 6.7037 | 29.11 | 1.02908 | 116721 | 14113 | -2.56 | 5.95 |  | period doubling (b_h crosses -2) |
| 0.05810 | 6.6816 | 29.01 | 1.03074 | 114356 | 14750 | -2.18 | 5.54 |  | period doubling (b_h crosses -2) |
| 0.12612 | 7.2581 | 31.52 | 0.99261 | 44771 | 93 | 638 | 246 |  | impact at floor |
| 0.15581 | 6.5699 | 28.53 | 1.04240 | 100652 | 19234 | -0.501 | 3.81 |  |  |
| 0.35857 | 6.4206 | 27.88 | 1.07981 | 72355 | 33613 | 1.01 | 2.42 |  |  |
| 0.45014 | 6.3829 | 27.72 | 1.10103 | 60306 | 41769 | 1.25 | 2.23 |  |  |
| 0.58126 | 6.3482 | 27.57 | 1.13193 | 44681 | 53647 | 1.41 | 2.12 |  |  |
| 0.81832 | 6.3136 | 27.42 | 1.18053 | 22161 | 72330 | 1.51 | 2.06 |  |  |
| 0.95110 | 6.3013 | 27.36 | 1.20168 | 12668 | 80461 | 1.53 | 2.05 |  |  |
| 1.11365 | 6.2889 | 27.31 | 1.22132 | 3819 | 88010 | 1.52 | 2.04 |  |  |
| 1.21722 | 6.2818 | 27.28 | 1.23031 | -332 | 91463 | 1.51 | 2.04 |  | impact at the physical surface |
| 1.21802 | 6.2817 | 27.28 | 1.23036 | -360 | 91485 | 1.51 | 2.04 |  | impact at the physical surface |

#### F3 (241 distinct members)

| C | T (TU) | T (d) | x0 | perigee alt (km) | periselene alt (km) | b_h | b_v | cand. | note |
|---|---|---|---|---|---|---|---|---|---|
| 1.92811 | 6.2191 | 27.01 | 1.21982 | 10 | 87432 | 1.06 | 2.06 |  | impact at floor |
| 2.13527 | 6.1801 | 26.84 | 1.19164 | 9893 | 76599 | 0.59 | 2.1 |  |  |
| 2.22742 | 5.5010 | 23.89 | 0.99256 | 32123 | 75 | 440 | 96.8 |  | impact at floor |
| 2.22975 | 5.4999 | 23.88 | 0.99258 | 32291 | 80 | 437 | 96.3 |  | impact at floor |
| 2.31035 | 5.4656 | 23.73 | 0.99316 | 38308 | 303 | 355 | 78.5 |  |  |
| 2.40938 | 5.4331 | 23.59 | 0.99426 | 46314 | 728 | 255 | 58 |  |  |
| 2.46595 | 6.0228 | 26.15 | 1.11532 | 37442 | 47261 | -1.95 | 2.46 |  | period doubling (b_h crosses -2) |
| 2.47786 | 6.0118 | 26.11 | 1.11159 | 38780 | 45830 | -2.15 | 2.5 |  | period doubling (b_h crosses -2) |
| 2.48689 | 5.4187 | 23.53 | 0.99572 | 53032 | 1289 | 179 | 43.1 |  |  |
| 2.53356 | 5.4174 | 23.52 | 0.99711 | 57238 | 1821 | 135 | 34.7 |  |  |
| 2.58078 | 5.8781 | 25.53 | 1.07481 | 51775 | 31691 | -4.59 | 3.18 |  |  |
| 2.59843 | 5.4321 | 23.59 | 1.00051 | 63162 | 3128 | 75.2 | 23.4 |  |  |
| 2.63551 | 5.4596 | 23.71 | 1.00453 | 66384 | 4676 | 41.5 | 16.9 |  |  |
| 2.66272 | 5.6479 | 24.53 | 1.03146 | 65346 | 15026 | -4.07 | 6.08 |  |  |
| 2.66661 | 5.6208 | 24.41 | 1.02726 | 66337 | 13414 | -2.56 | 6.73 |  | period doubling (b_h crosses -2) |
| 2.66668 | 5.6202 | 24.41 | 1.02718 | 66357 | 13379 | -2.52 | 6.75 |  | period doubling (b_h crosses -2) |
| 2.66755 | 5.6121 | 24.37 | 1.02594 | 66624 | 12903 | -1.93 | 6.98 |  | period doubling (b_h crosses -2) |
| 2.66953 | 5.5767 | 24.22 | 1.02072 | 67577 | 10898 | 1.59 | 8.12 |  |  |

#### F4 3/7-r (38 distinct members)

| C | T (TU) | T (d) | x0 | perigee alt (km) | periselene alt (km) | b_h | b_v | cand. | note |
|---|---|---|---|---|---|---|---|---|---|
| 0.89966 | 18.9499 | 82.29 | 1.02966 | 28880 | 14335 | -1.8 | -0.447 | yes | period doubling (b_h crosses -2) |
| 0.92277 | 18.9217 | 82.17 | 1.03578 | 27079 | 16687 | -2.52 | -0.565 | yes |  |
| 0.94157 | 18.9035 | 82.09 | 1.04045 | 25575 | 18481 | -2.99 | -0.647 | yes |  |
| 0.95567 | 18.8917 | 82.04 | 1.04383 | 24442 | 19782 | -3.32 | -0.718 | yes |  |
| 0.96578 | 18.8840 | 82.00 | 1.04621 | 23632 | 20696 | -3.55 | -0.778 | yes |  |
| 0.97284 | 18.8790 | 81.98 | 1.04785 | 23068 | 21327 | -3.72 | -0.825 | yes |  |
| 0.98016 | 18.8740 | 81.96 | 1.04953 | 22487 | 21973 | -3.89 | -0.879 | yes |  |
| 0.99162 | 18.8667 | 81.93 | 1.05213 | 21583 | 22973 | -4.17 | -0.974 | yes |  |
| 1.00992 | 18.8559 | 81.88 | 1.05620 | 20159 | 24538 | -4.64 | -1.16 | yes |  |
| 1.03982 | 18.8403 | 81.81 | 1.06265 | 17896 | 27017 | -5.51 | -1.56 | yes |  |
| 1.09008 | 18.8182 | 81.72 | 1.07293 | 14296 | 28266 | -7.35 | -2.55 | yes |  |
| 1.17678 | 18.7874 | 81.58 | 1.08894 | 8763 | 23479 | -12.2 | -5.53 | yes |  |
| 1.30792 | 18.7539 | 81.44 | 1.10876 | 2087 | 17957 | -24.8 | -13.6 | yes |  |
| 1.35473 | 18.7460 | 81.40 | 1.11456 | 189 | 16488 | -30.5 | -17.2 |  | impact at floor |

#### F4 2/3-r (76 distinct members)

| C | T (TU) | T (d) | x0 | perigee alt (km) | periselene alt (km) | b_h | b_v | cand. | note |
|---|---|---|---|---|---|---|---|---|---|
| -0.63672 | 14.2350 | 61.82 | 1.00785 | 220985 | 5951 | 1.11 | 0.274 |  |  |
| -0.63615 | 14.2143 | 61.73 | 1.00848 | 222296 | 6194 | -1.26 | 0.909 |  | period doubling (b_h crosses -2) |
| -0.63576 | 14.2056 | 61.69 | 1.00873 | 222756 | 6291 | -2.17 | 1.13 |  | period doubling (b_h crosses -2) |
| -0.63227 | 14.1569 | 61.48 | 1.01005 | 224668 | 6795 | -6.44 | 2.03 |  |  |
| -0.52643 | 13.5372 | 58.78 | 1.02498 | 216086 | 12534 | -11.5 | 2.53 |  |  |
| -0.44418 | 13.2042 | 57.34 | 1.03891 | 200508 | 17890 | -6.08 | 1.67 |  |  |
| -0.41182 | 13.0955 | 56.87 | 1.04645 | 193508 | 20788 | -4.45 | 1.5 |  |  |
| -0.39025 | 13.0304 | 56.58 | 1.05238 | 188590 | 23067 | -3.56 | 1.44 |  |  |
| -0.36283 | 12.9564 | 56.26 | 1.06111 | 182065 | 26423 | -2.66 | 1.4 |  | period doubling (b_h crosses -2) |
| -0.34880 | 12.9225 | 56.12 | 1.06614 | 178610 | 28356 | -2.29 | 1.4 |  |  |
| -0.32488 | 12.8708 | 55.89 | 1.07559 | 172574 | 31992 | -1.79 | 1.41 |  | period doubling (b_h crosses -2) |
| -0.23842 | 12.7437 | 55.34 | 1.11778 | 149929 | 48206 | -1.05 | 1.51 |  |  |
| -0.15159 | 12.6828 | 55.07 | 1.16580 | 127991 | 58633 | -1.24 | 1.55 |  |  |
| -0.11458 | 12.6683 | 55.01 | 1.18600 | 119313 | 58247 | -1.46 | 1.55 |  |  |
| -0.08478 | 12.6596 | 54.97 | 1.20181 | 112656 | 57989 | -1.66 | 1.55 |  |  |
| -0.03495 | 12.6493 | 54.93 | 1.22716 | 102178 | 57703 | -2.03 | 1.54 |  | period doubling (b_h crosses -2) |

#### F4 3/4-r (56 distinct members)

| C | T (TU) | T (d) | x0 | perigee alt (km) | periselene alt (km) | b_h | b_v | cand. | note |
|---|---|---|---|---|---|---|---|---|---|
| -0.68198 | 20.3779 | 88.49 | 1.01119 | 236669 | 7234 | -0.188 | -3.32 |  | period doubling (b_h crosses -2) |
| -0.63242 | 19.8351 | 86.13 | 1.02414 | 239031 | 12212 | -28.1 | -2.41 |  |  |
| -0.54272 | 19.2315 | 83.51 | 1.04581 | 221754 | 20542 | -19.4 | -1.73 |  |  |
| -0.48867 | 18.9644 | 82.35 | 1.06560 | 207591 | 27069 | -16.3 | -2.19 |  |  |
| -0.45834 | 18.8439 | 81.83 | 1.07963 | 198823 | 25603 | -16.5 | -3.07 |  |  |
| -0.44225 | 18.7880 | 81.59 | 1.08785 | 194023 | 24692 | -17.3 | -3.81 |  |  |
| -0.42588 | 18.7361 | 81.36 | 1.09667 | 189079 | 23702 | -18.6 | -4.8 |  |  |
| -0.39332 | 18.6461 | 80.97 | 1.11518 | 179215 | 21655 | -23.2 | -7.67 |  |  |
| -0.32839 | 18.5050 | 80.36 | 1.15355 | 160165 | 17770 | -42.3 | -18.3 |  |  |
| -0.18372 | 18.2877 | 79.41 | 1.23376 | 123305 | 11629 | -164 | -78.4 |  |  |
| 0.02625 | 18.0948 | 78.58 | 1.32951 | 82237 | 7074 | -553 | -237 |  |  |
| 0.29524 | 17.9972 | 78.15 | 1.42351 | 44626 | 4297 | 223 | -90 |  | period doubling (b_h crosses -2) |


## 2. Verified vs assumed

- Verified: every listed member is a corrected symmetric periodic orbit. It has |xdot(T/2)| <
  1e-11, or < 1e-7 at the noise floor of the very sensitive F5 orbits, and det M4 = 1. The
  catalogue matches are crossing-state matches, not C-only matches. Newton's 11 rows and
  Hoelker-Winston's Figs. 89 and 94 reproduce.
- Assumed or not checked:
  - Restrepo-Russell membership (data not held).
  - Behaviour inside the period-doubling gaps.
  - Members between sampled steps; the steps are adaptive, and every candidate interval is
    bounded by sampled members.
  - The planar model only. The vertical index is computed but no 3-D family is followed.

## 3. Hand-offs

- `#972` (literature): the three candidate stretches above are known-class by construction; no
  novelty check needed unless the owner wants the F1 interior (perigee 200-35,000 km) checked
  against RR's Earth-Moon files.
- `#1019` / lead: whether the RR Earth-Moon folder is fetched.
- `#948`: F1's down-perigee end (Earth-surface impact at C = 1.09) and F2/F3's floor stops are
  natural entry points for a both-primary regularised continuation to collision.
