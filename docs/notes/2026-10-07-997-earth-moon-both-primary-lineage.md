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

## 1. Results

(filled in after the runs; see sections below)
