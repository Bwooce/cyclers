# #1000: Earth-Moon complement search (asymmetric, Earth-circulating periodic orbits that the Franz-Russell and Restrepo-Russell databases exclude by construction)

Task `#1000` (earthmoon-opus, 2026-10-08). The owner opened it in the ruling of 2026-10-07 ~21:50
AEDT, and the Earth-Moon gate is open for this task. Source:
`docs/notes/2026-10-05-938-fable-corpus-novel-paths-review.md` sec. 2 (line 155) and sec. 4.
**Nothing here is called novel. The literature step is deferred to `#972`.** A clean "only known
classes found" is a result.

## 0. Pre-registration (committed before any `#1000` computation)

### 0.1 Why this region is a complement

- **Restrepo & Russell 2018** cover only planar AXISYMMETRIC orbits (p.3; `#997` sec. 0.6). Their
  Earth-Moon files hold no record with J < 2.119 (`#997` sec. 1.8).
- **Franz & Russell 2022** discard any orbit that is ever more than about 350,000 km from the Moon.
  Their planar grid is built from perpendicular-crossing starts, so it too holds symmetric orbits
  only (`docs/notes/2026-07-28-747-franz-russell-casoliva-crosscheck.md`).
- So an ASYMMETRIC Earth-Moon periodic orbit with a low perigee and a lunar pass is in neither
  database by construction. Its absence from them says nothing about novelty.
- The catalogue already holds such orbits: casoliva-7-3b and casoliva-7-3c (asymmetric, C =
  1.0687/1.0672, T = 6 pi, perigee about 9,300 km alt, periselene about 11,450 km alt; `#997`
  `gate.json`). They are this search's positive controls.

### 0.2 Model, floors, units

- Planar CR3BP at the registry mu = `cr3bp_system("Earth", "Moon").mu` = 0.01215058439469525,
  L = 384,400 km, TU = 375,190.26 s.
- **Floors:** the registry values, Earth 200 km (`PLANETS["E"].safe_alt_km`) and **Moon 100 km**
  (`SATELLITES["Moon"].safe_alt_km` = 100.0, checked 2026-10-08; the lead's brief said 50 km). I
  use the registry's 100 km. I also report which candidates would pass at 50 km, for the
  owner's information. The 50 km column is not a test.
- Minimum distances are taken over the WHOLE period from event-refined extrema of r1 and r2.

### 0.3 Search space: a perigee return map at fixed C

- **Section.** Earth perigee passages (dr1/dt = 0 with d2r1/dt2 > 0). A point on the section at
  fixed C is (r_p, omega): perigee radius, and argument of perigee measured from the +x axis in the
  rotating frame. The velocity is perpendicular to the Earth radius vector, its magnitude comes
  from C, and its sense is prograde or retrograde about the Earth (two sheets).
- **Map.** P: (r_p, omega) -> the next perigee. A periodic orbit with k perigees per period is a
  fixed point of P^k.
  - Symmetric orbits appear as fixed points too. The method does not assume symmetry, so it
    finds both kinds.
- **Grid** (fixed before running):
  - C from 0.5 to 3.1 in steps of 0.1 (27 levels). 3.1 is below C_L1 = 3.188; above that, Earth
    and Moon cannot exchange.
  - r_p: 8 values, log-spaced from 6,578 km (the Earth floor) to 42,164 km (GEO radius).
  - omega: 36 values, step 10 deg.
  - Both senses.
  - k = 1 ... 7. One integration per seed records the first 7 perigee returns, so all k come from
    the same run. k = 7 is needed for the 7-3 controls.
  - Total: 27 x 8 x 36 x 2 = 15,552 seeds.
- **Seed filter.** A seed is a fixed-point candidate for P^k if three things hold:
  - the scaled residual |P^k(z) - z|, with r_p in units of 6,578 km and omega in radians, is a
    local minimum over its 8 grid neighbours (omega periodic);
  - that residual is below 0.3;
  - the k-return arc passes the Moon inside the lunar Hill radius, 61,600 km.
- **Refinement.**
  - Newton on P^k(z) - z = 0 in (r_p, omega) at fixed C, with a central finite-difference
    Jacobian (h = 1e-7, scaled); at most 15 iterations.
  - Converged: |P^k(z) - z| < 1e-9 scaled, and the full-period state closure is < 1e-8.
  - Then the orbit is reduced to its minimal period: k is replaced by the smallest divisor k'
    for which P^k'(z) = z.
- **Continuation** of each NEW asymmetric orbit (not of symmetric ones; those are `#997`'s and
  RR's), in C, both directions:
  - steps of 0.01, with Newton seeded from the previous member; at most 40 steps per direction.
  - Stop on any of `#997` sec. 0.4's rules: impact at either floor; a stability index crossing
    -2 (period doubling); loss of convergence after 3 halvings of the step; a topology change
    (the k' perigee count, or the winding numbers about either primary, change).

### 0.4 Classification of each converged orbit

- **Symmetric or asymmetric.** Symmetric means the orbit has a perpendicular x-axis crossing
  (y = 0, |xdot| < 1e-7) within one period; then it is mirror-symmetric (two such crossings per
  period). Asymmetric means it has none.
- **Recorded for each orbit:**
  - C and T (TU and days), and k';
  - perigee and periselene distances and altitudes (whole-period minima), and the apogee;
  - the maximum Moon distance;
  - the planar monodromy eigenvalues and b = lambda + 1/lambda for the nontrivial pair, plus the
    vertical index k_perp = tr(Mz);
  - the winding numbers about the Earth and the Moon.
- **Cycler-class candidate.** The test is `#997` sec. 0.5: every period, at least one Earth pass
  and one Moon pass, with:
  - the perigee at or above the Earth floor and at or below the GEO radius;
  - the periselene at or above the Moon floor and inside the Hill radius.
- **Dedupe.** Two orbits are the same if their C and T match to 1e-8 and one's perigee set
  contains the other's start to 1e-7. A mirror pair (omega -> -omega with time reversed) is
  counted as ONE family with two branches.

### 0.5 Known-class gates (applied to every converged orbit)

1. **Catalogue rows.** Each planar Earth-Moon row's state is propagated, and its perigee section
   points are matched to 1e-6 scaled, with T to 1e-6 relative.
   - The controls casoliva-7-3b and 7-3c must be recovered.
   - The symmetric rows found incidentally are labelled with their row id.
2. **`#997` families.** Symmetric orbits are matched against `data/997_lineage/*.jsonl` by C,
   interpolated, and by crossing state, as in `#997` sec. 1.2.
3. **Restrepo-Russell 2018.**
   - Symmetric orbits: the RR cross-match of `#997` sec. 1.8 (same rule and tolerances).
   - Asymmetric orbits are excluded by RR's construction. That is recorded; it is not a lookup.
4. **Franz-Russell 2022.** An orbit is excluded when its maximum Moon distance is above
   350,000 km. FR is also symmetric-only.
5. **Asymmetric families and pitchforks.** An asymmetric family is usually born at a symmetric
   family's pitchfork, where b crosses +2. For each new asymmetric family, the continuation
   end-points are checked for a symmetric orbit at the same C and T (within 1e-4). Any parent
   found is named. A parent in a published family makes the asymmetric branch a known-class
   candidate for `#972` to adjudicate.

### 0.6 Positive controls and acceptance

- **Required: casoliva-7-3b and 7-3c** (asymmetric, k = 7, C about 1.07) are recovered by the grid
  plus Newton, without seeding from the rows.
  - The grid level nearest to them is C = 1.1. Two extra C levels, 1.0687 and 1.0672 (the rows'
    own C values), are added to the scan for the control only.
  - **If they are not recovered, the search formulation is rejected.** I stop and report. No
    "nothing found" will be stated from a formulation that misses its own controls.
- **Also expected:**
  - symmetric members of `#997` F1 (k = 1, C 0.58-1.08);
  - F3 (k = 2 or 1, C 1.93-2.67);
  - F4 3/7-r (k = 7, C 0.90-1.35);
  - Vaquero 2:1 rows (C 1.98-2.66).
  These are reported as the recovery rate on the levels where they exist. A low rate is reported
  as a coverage limit of the grid; it is not an error.

### 0.7 Expected outcome (stated before running)

- Most converged orbits will be symmetric and known-class (`#997` families, Vaquero, Casoliva,
  RR).
- Asymmetric orbits should be rarer. They should appear near the pitchforks of symmetric
  families, and the Casoliva 7-3b/c class is one such.
- My prior: about 0-5 new asymmetric families at this grid resolution, most with a published
  symmetric parent.
- Any asymmetric cycler-class candidate is outside RR and FR by construction. That is a reason
  to send it to `#972`, not evidence of novelty.

### 0.8 Compute plan

- Script `scripts/run_1000_complement.py`, calling `preflight_search()` with task 1000, subcommands
  `scan`, `refine`, `continue` and `gate`. Data in `data/1000_complement/`.
- **Checkpoints.** One JSONL per C level and sheet, appended and flushed per seed. A re-run skips
  seeds that are done.
- **Workers.** At most 2 (`multiprocessing` Pool(2)); CI also uses the machine.
- **Every call runs in the FOREGROUND, under 8 minutes** (`--max-seconds 420`), and resumes from
  its checkpoint. No backgrounding (papercut 2026-10-08).
- **Estimate.** About 0.5 s per seed for 7 returns, so about 65 minutes of wall time at 2 workers
  over about 12 foreground calls for the scan. Refinement and continuation cost about 1-3 s per
  Newton iteration.
- Before the scan, a timing pilot of 50 seeds. If the measured cost is more than 2x the estimate,
  I report it and either hand the lead one launch command or reduce the grid. A reduction is
  recorded as an amendment BEFORE it is run.

### 0.9 Amendments (each recorded before the runs it affects)

- 2026-10-08, before the pilot: the seed residual is |P^k(z) - z| with z = (ln r_p, omega), not r_p
  in units of 6,578 km. The r_p grid is log-spaced (ratio 1.31 per step, ln step 0.27), so the
  logarithm makes the 0.3 threshold about one grid cell in both coordinates. Newton refinement
  uses the same z.

- 2026-10-08 11:45 AEDT, after the scan and before any Newton run. The scan finished: 16,704 seeds
  at 0.02-0.13 s per seed, and 329 seeds pass the candidate filter. The tolerances are set to
  what a 1e-12 integration of orbits with lambda up to about 1e3 can deliver:
  - Newton converged: |P^k(z) - z| < 1e-8 (was 1e-9);
  - full-period closure < 1e-7 (was 1e-8);
  - symmetric test: |xdot| < 1e-6 at a y = 0 crossing (was 1e-7).
  Also: the k-minimality test uses |P^d(z) - z| < 1e-6, and the control-level candidates are
  refined first.

## 1. Results

(after the runs)
