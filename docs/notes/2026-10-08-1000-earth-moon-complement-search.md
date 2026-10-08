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

## 1. Results: the formulation FAILS its positive control and is rejected (sec. 0.6)

Runs 2026-10-08 11:20-11:47 AEDT.
- The scan is complete: 16,704 seeds, 0.02-0.13 s per seed, Pool(2), foreground calls of at most
  7 minutes. 329 seeds passed the candidate filter.
- Newton ran on 96 of them, control levels first. **None converged.** 76 reached the iteration
  limit, and 20 lost the k-th return.

**Why: a check on the control itself** (scratch scripts; numbers from the catalogue row
casoliva-7-3b, which is asymmetric, C = 1.068655, T = 18.8496):
- **The map is right.** The row's orbit has exactly 7 perigees per period (radii 15,710-33,919 km,
  all inside the grid's r_p range). At its 3rd perigee, z = (ln r_p, omega) = (-3.19736, 2.41093
  rad), retrograde sheet, the residual is |P^7(z) - z| = 4.5e-8. Newton started there converges
  in 2 iterations to 8.5e-10.
- **The basin is far smaller than the grid.** The residual grows about 66 times faster than a
  perturbation in ln r_p:
  - perturbation 0.001 gives residual 0.066;
  - perturbation 0.005 gives 0.22, already nonlinear;
  - perturbation 0.05 gives 0.79.
  The monodromy has lambda about 57 for 7-3b. Newton from a 0.05 offset fails. The grid spacing
  is 0.27 in ln r_p and 0.17 rad in omega.
- **The filter cannot see the control.** At the cells next to the control fixed point (C =
  1.0687, retrograde sheet, k = 7), the scan residuals are 0.47-3.1, all above the 0.3 threshold.
  Over the whole level, only 14 cells are below 0.5.

**Conclusion.**
- A coarse return-map grid with a residual-threshold filter cannot find strongly unstable fixed
  points. The Earth-Moon cycler-class orbits are all strongly unstable (lambda 10-1e3). Their
  linear basin in the section is about 1e-3 wide, 100 times finer than any affordable global grid.
- Under the pre-registration, this formulation is rejected. **No result of the form "nothing
  found" is stated from it.** The 96 non-convergences are a property of the method, not of the
  region.
- Data kept: the scan JSONL (`data/1000_complement/scan/`, 6.9 MB, not committed; kept locally
  for re-use) and `data/1000_complement/refined.jsonl` (the 96 failed refinements).

**Proposed redesign** (for the lead's decision; nothing has run):
1. **Multiple shooting from near-Keplerian resonant skeletons.** This is the construction Casoliva
   et al. used for the 7-3 class.
   - For each p:q Earth resonance (a = (q/p)^(2/3) about the Earth, with the apocentre or
     pericentre chosen to meet the Moon), take the Kepler ellipse at an ARBITRARY orientation
     and phase (that is what lets asymmetric orbits appear).
   - Correct it by multiple shooting: N arcs, one per perigee, with continuity, fixed C and
     periodicity constraints, and the phase condition removed.
   - Multiple shooting spreads the instability over N arcs, so the basin per arc is lambda^(1/N)
     wider.
   - Positive control first: 7-3b/c from p:q = 7:3 skeletons.
2. **Symmetry breaking at pitchforks.** Asymmetric families branch from symmetric ones where b
   crosses +2.
   - Locate b = +2 points along the `#997` families and along RR's symmetric records (RR stores
     b_h).
   - Branch along the antisymmetric eigenvector.
   - This finds only asymmetric families connected to known symmetric parents, so it needs no
     separate known-class gate, but it is not a complement search.
3. A local fine grid (ln r_p and omega steps of about 1e-3) around the scan's best cells. This is
   too expensive globally: about 1e5 times the current grid.

My recommendation is (1) with the 7-3b/c control, then (2) as the pitchfork-parent gate. Each
needs its own pre-registration amendment before it runs.

**Archive of the rejected scan.** The rejected scan was copied to
`~/dev/references/cyclers-runs/1000/scan/` (outside the repo), with `cp -R`. I confirmed the copy
with `diff -rq`: identical, 58 files, 16,704 records.
- md5 of the concatenated JSONL (`cat *.jsonl | md5`, in name order): a70eff2a8c5541e0f1768ff39fd27dfc.
- md5 of the per-file md5 list: f230a4c921a94f25c6b304c5d340f8e2.
- The working-tree copy was then deleted. `refined.jsonl` stays committed.

## 2. Redesign: pre-registration amendment A (lead ruling 2026-10-08; committed before any run of it)

### 2.1 Method (1): multiple shooting from near-Keplerian p:q skeletons

- **Skeleton.** An Earth-centred two-body ellipse (GM = 1 - mu) in the INERTIAL frame, either
  sense, with period T_P = 2 pi q / p. So p revolutions take q lunar sidereal months (2 pi q TU),
  and a = (q/p)^(2/3) (1 - mu)^(1/3).
  - It is given by its perigee radius r_p and the rotating-frame angle omega of the perigee at
    t = 0. omega carries both the orientation and the phase relative to the Moon.
  - It must cross the lunar orbit: r_p < 1 < r_a = 2a - r_p.
- **Seed grid** (fixed now):
  - **p:q in {2:1, 1:2, 3:2, 5:2, 1:3, 2:3, 4:3, 5:3, 7:3, 8:3}.** These are all ratios with
    q <= 3 (T <= 6 pi) whose ellipse can cross the lunar orbit with r_p at or above the Earth
    floor. 3:1 and 7:2 cannot (2a < 1), and 4:1 cannot either.
  - Sense: inertial prograde or retrograde.
  - r_p: 8 values, log-spaced from 6,578 to 42,164 km, keeping only those with r_a > 1.
  - omega: 36 values, step 10 deg.
  - At most 10 x 2 x 8 x 36 = 5,760 seeds.
- **Multiple shooting.**
  - N = p arcs, the nodes at the skeleton's p perigees (one per revolution), each node mapped to
    the rotating frame.
  - Unknowns: the N planar node states and the N arc durations (5N).
  - Constraints: continuity at every node, with the last arc closing on node 0 (4N); node 0 at
    a perigee, (r1 . v) = 0 (1); and, in the control runs only, C = C_target (1).
  - Newton step: the minimum-norm least-squares solution. The STM of each arc comes from the
    4 x 4 variational equations (DOP853, 1e-12).
  - At most 25 iterations; the step is limited to 0.1 in the node-state norm.
  - Converged when the continuity residual is < 1e-10 (nondimensional). Then:
    - a full single-shooting closure from node 0 over T = sum of the durations must be < 1e-6.
      This is looser than sec. 0.9 because lambda^(1/1) applies over the whole period;
    - the k-minimal reduction is applied;
    - the orbit is classified as in sec. 0.4.
- **Control FIRST, hard stop.** Only the 7:3 seeds are run first: C fixed to the row values,
  1.068655371747616 (7-3b) and 1.0671969118897233 (7-3c), both senses (1,152 solves).
  - **Pass:** a converged ASYMMETRIC orbit matches each row. Its perigee set contains the row's
    perigees to 1e-6 (ln r_p, omega) and T matches to 1e-6 relative.
  - **If either row is not recovered, method (1) is rejected and I stop.**
  - Then the full grid runs with C free (no C constraint). Every converged orbit is classified,
    deduped (sec. 0.4) and gated (sec. 0.5).
  - Each new asymmetric cycler-class family is continued in C under sec. 0.3's rules.

### 2.2 Method (2): symmetry-breaking at b = +2 pitchforks (also the parent gate)

- **Parents.**
  - (a) Every adjacent member pair along the `#997` family checkpoints
    (`data/997_lineage/*.jsonl`) where the in-plane index b_h crosses +2.
  - (b) Restrepo-Russell Earth-Moon records with |b_h - 2| < 0.02 whose orbit has a perigee at
    or below the GEO radius and a lunar pass inside the Hill radius (checked by propagation).
    At most 200 of these, chosen by smallest |b_h - 2|.
- **Branching.**
  - At the bisected pitchfork member (b_h = 2 to 1e-6), take the monodromy eigenvector of the
    double eigenvalue +1 that is ANTISYMMETRIC under the x-axis reflection.
  - Perturb along it by +-1e-4. Solve with the multiple shooting of 2.1 at C = C_parent +- 1e-5
    (both sides, both signs).
  - Accept an asymmetric solution that does not collapse back onto the parent: its distance in
    the section is > 1e-5.
  - Continue each branch in C with sec. 0.3's rules.
- **Control FIRST, hard stop.** The 7-3b/c family recovered by (1) is continued in C until its
  b_h reaches +2 or it meets a symmetric orbit. That is its pitchfork parent.
  - Method (2), applied to that parent, must regenerate the 7-3b/c branch.
  - **If (1) finds no pitchfork along 7-3b/c within its continuation range, the control for (2)
    is "not applicable".** That is recorded; then (2) runs as a parent gate only, and makes no
    complement claim.
  - **If the parent exists and (2) fails to regenerate the branch, method (2) is rejected.**

### 2.3 Compute

- Script `scripts/run_1000_shooting.py` (preflight task 1000), Pool(2).
- Foreground calls under 8 minutes (`--max-seconds 420`), with JSONL checkpoints under
  `data/1000_complement/shooting/`.
- A timing pilot of 20 control seeds goes first.

- **Amendment A2** (2026-10-08 11:58 AEDT). Made after the 20-seed timing pilot and before any
  control run.
  - What the pilot showed: 17 of 20 solves impacted, and 3 hit the iteration limit. The skeleton
    transform is right: with mu = 1e-9 one arc closes on the next node to 6e-7. But with the real
    mu, the lunar perturbation shifts the revolution time by a few per cent. At a PERIGEE node,
    that becomes an O(1) velocity mismatch, because the angular rate there is about 100 rad/TU.
  - Change: the multiple-shooting nodes move to the skeleton's p APOCENTRES, where the dynamics
    are slow. The node-0 phase condition (r1 . v = 0) is unchanged, and now picks an apogee.
    Classification and the k-minimal test start from the orbit's first perigee after node 0.


## 3. Method (1) control: PASSED (2026-10-08 12:10 AEDT)

- **Run.** 7:3 seeds at the two row C values, 1,152 solves at about 0.5 s per seed with Pool(2),
  in two foreground calls. Output: `data/1000_complement/shooting/control.jsonl`.
  - Status: 176 converged, 950 impacted (the arcs pass inside the Earth or the Moon during the
    iterations), 26 reached the iteration limit.
- **Converged orbits.**
  - Symmetric, T = 18.827-18.828 (k = 7): 126. Their perigee and periselene are about 15,800
    and 29,300 km alt, which is the `#997` F4 3/7-r family at these C. There are also 3 other
    symmetric k = 7 orbits (T = 18.446 and 18.961).
  - Asymmetric: 47. Every one is either the catalogue row or its mirror image under
    (x, y, t) -> (x, -y, -t), matched by perigee section point to 1e-6 and T to 1e-6:
    - casoliva-7-3b-em-cycler-2010 (C = 1.068655): **22 seeds converge onto the row**, and 1
      onto its mirror. Best section distance 2.2e-11.
    - casoliva-7-3c-em-cycler-2010 (C = 1.067197): **1 seed onto the row** (distance 8.3e-10),
      and 23 onto its mirror.
- **Both positive controls are recovered from skeletons alone**, without seeding from the rows.
  Method (1) is accepted and proceeds to the full grid (sec. 2.1) with C free.

## 4. Method (1) grid, continuations, and the method (2) control: PASSED (2026-10-08 13:40 AEDT)

### 4.1 Method (1) full grid (C free)

- **Run.** 5,760 seeds in `data/1000_complement/shooting/grid.jsonl`.
  - Status: 1,505 converged, 3,101 impacted during the iterations, 650 reached the iteration
    limit, and 504 skipped because their skeleton does not cross the lunar orbit.
  - Of the 1,505 converged orbits, 213 are cycler-class candidates: 191 symmetric and 22
    asymmetric.
  - The asymmetric candidates form four families (p:q, inertial sense, winding about the Earth,
    winding about the Moon):

| family | grid members | C | T (TU) | perigee alt (km) | periselene alt (km) | b (in-plane, from arc STMs; sec. 6) | max Moon distance (km) |
|---|---|---|---|---|---|---|---|
| A: 5:2 prograde, wE 3, wM 0 | 3 | 2.478-2.592 | 12.71-12.73 (about 4 pi) | 3,164-8,748 | 13,153-16,348 | -110 to -121 | about 790,000 |
| B: 5:3 retrograde, wE -8, wM -2 | 4 | 0.728-0.864 | 18.89-18.90 (about 6 pi) | 3,385-10,789 | 18,881-21,188 | +115 to +130 | about 950,000 |
| C: 7:3 retrograde, wE -10, wM -1 | 13 | 0.882-1.163 | 18.80-19.12 | 4,229-26,908 | 5,673-12,520 | +2.0 to +80 | about 800,000 |
| D: 7:3 prograde, wE 4, wM -1 | 2 | 2.422-2.547 | 18.66-18.82 | 365-5,754 | 10,768-13,774 | +512 to +718 | about 810,000 |

- **Continuations** (`data/1000_complement/continuation/`, `scripts/run_1000_continue.py`).
  - Fixed-C multiple shooting, with the sec. 0.3 rules.
  - Two corrector aids, neither a rule change:
    - a secant predictor, with a fall-back to the previous nodes;
    - each new member keeps the perigee nearest to the previous member's.
  - b is taken from the product of the per-arc STMs. On these orbits, single shooting over the
    whole period drifts off the orbit and gives a wrong b: one member gave -508 by single
    shooting against +115 from the arc STMs.

| family | continued C range | distinct members (candidates) | stop (each direction) |
|---|---|---|---|
| A | 2.4109-2.5804 | 81 (81) | max steps / max steps |
| B | 0.6703-0.9418 | 77 (76) | max steps / Earth floor at C 0.9418 |
| C | 0.9239-1.1188 | 43 (43) | loss of convergence at 0.924 / max steps |
| D | 2.4136-2.5601 | 33 (32) | Earth floor at 2.4136 / loss of convergence at 2.5601 |

Counts are de-duplicated: the seed appears in both direction files, and stop rows repeat a
member. A's grid member at C = 2.5916 lies just beyond A's continued range (which ends at
2.5804). It has the same topology and a b on trend (-121 against -120 at 2.580). It is
CONSISTENT with family A, but that is not shown by continuation.

- Family C passes through casoliva-7-3b. Interpolated at the row's C = 1.068655, the
  continuation members (C = 1.0658 and 1.0700) give T = 18.84964 against the row's 18.8496,
  and perigee 9,333 km alt against the row's 9,332. **Family C is the Casoliva 7-3b/c class.**
  It is known-class.
- No continuation of A, B or D reached b = +2 or met a symmetric orbit. Within the computed
  ranges, **no pitchfork parent is found for A, B or D** (sec. 0.5 item 5).

### 4.2 Method (2) control (see sec. 6, amendment A3: the pre-registered variant FAILED; the replacement passed)

- **Approach** (`scripts/run_1000_pitchfork.py control`).
  - Family C's lowest-C grid member (C = 0.881505, T = 19.12309, b = 2.023) was continued
    downward. Its minimum |xdot| at a y = 0 crossing fell from 9.5e-3 to 4e-10: the asymmetric
    branch closes onto a SYMMETRIC orbit.
  - Bisection on that symmetric family gives the **pitchfork parent: C = 0.881428101,
    T = 19.123728, b = +2.000000** (`control_approach.json`; states in
    `control_parent_and_seed_states.json`).
- **Branching.**
  - The antisymmetric-eigenvector perturbation with C fixed failed: all 36 solves were
    ill-conditioned at the bifurcation. It was replaced by the standard symmetry-breaking
    parameter: node 0 is placed at the parent's perpendicular crossing and solved with
    y(node 0) = 0 and xdot(node 0) = eta, with C free.
  - For eta = +-1e-4 ... +-3e-2, every solve converges to an asymmetric orbit
    (`control_branch.json`). C rises roughly as eta^2, from 0.881428 to 0.882275. The +eta
    and -eta solutions are mirror images: same C and T to 1e-12.
- **Match.** At eta = +-0.01 the branch is at C = 0.8815145, T = 19.1230159. Interpolated to the
  grid member's C = 0.881505, it gives T = 19.123094 against the grid member's 19.12309.
- **So the eta variant, applied at the parent, regenerates the 7-3b/c branch.** With 4.1, the
  chain is: pitchfork parent (C 0.8814) -> branch -> grid members (C 0.88-0.92) -> continuation
  (0.92-1.12) -> casoliva-7-3b (C 1.0687).
- **Status.** The branching step AS PRE-REGISTERED (eigenvector perturbation at fixed C +- dC)
  failed this control, 36 of 36. Under the hard-stop rule that rejects it. The eta variant was
  introduced after that failure and checked on the SAME control, so the control is not an
  independent test of it (sec. 6, A3). Results from the eta variant are labelled as such and
  await the lead's ruling.
- The parent is a symmetric 7:3 retrograde family at C = 0.8814 and T = 19.124.
  - It is not `#997` F4 3/7-r: that family stops at C 0.8997 with T 18.95.
  - It lies below RR's J range (2.119), so RR cannot contain it.
  - Its identity against the catalogue's symmetric casoliva-7-3a (C 1.0216, T = 6 pi) is
    checked next, by continuing the parent family (sec. 5).

## 5. Method (2) runs, gates and verdict (2026-10-08 14:00 AEDT)

### 5.1 Method (2) on the `#997` and RR parents

- **`#997` b = +2 crossings** (`data/1000_complement/pitchfork/parents997.jsonl`; eta variant,
  pending the ruling in sec. 6): six member pairs.
  - Folds (C extremal; no asymmetric branch, as expected): F2-1 near C -0.489, and F3-6/F3-7
    near C 2.6695. The JSONL's automatic guess marks F3-6/F3-7 "pitchfork?". The fold call
    rests on C being at its maximum there (`#997` F3's C range ends at 2.6695) and on no
    asymmetric branch appearing.
  - Pitchforks that give an asymmetric branch:
    - 2:3 retrograde at C = -0.620969, T 14.325. The branch has perigee about 206,000 km alt:
      **not cycler-class**.
    - 3:4 retrograde at C = 0.274826, T 17.999. The branch has perigee about 47,000 km alt,
      above the GEO radius: **not cycler-class**.
- **RR records:** 10,223 have |b_h - 2| < 0.02 (mostly near-Moon families). **None** has a perigee
  at or below the GEO radius and a lunar pass inside the Hill radius (checked by propagation at
  RR's own mu). So no RR parent is in cycler geometry.
- **Family C's parent** (sec. 4.2) is a symmetric 7:3 retrograde family at C = 0.8814,
  T = 19.124.
  - Continued upward to C = 0.939, its T rises to 19.157, moving away from casoliva-7-3a
    (C 1.0216, T 18.8496), which lies on `#997` F4 3/7-r. So it is a separate symmetric family.
  - It is below RR's J range.
  - File: `pitchfork/parent_family_toward_7-3a.json`.

### 5.2 Gates (`scripts/run_1000_gate.py`, `data/1000_complement/gate_summary.json`)

All 213 method (1) cycler-class candidates are excluded from Franz-Russell: their maximum Moon
distance is 750,000-1,950,000 km.

| family (p:q, sense, wE, wM) | sym | candidates | C | perigee alt (km) | periselene alt (km) | `#997` | RR | label |
|---|---|---|---|---|---|---|---|---|
| 2:1 pro, 1, -1 | S | 4 | 2.382-2.443 | 28,777-34,958 | 49,910-56,478 | F3 | 1 of 4 matched | known-class (Newton 1/2 direct = Vaquero 2:1) |
| 7:3 ret, -10, -1 | S | 70 | 0.864-1.341 | 707-30,987 | 6,902-29,312 | F4 3/7-r (C >= 0.90) | below J range | known-class (Newton's alpha/beta, contains casoliva-7-3a); the C < 0.90 part is family C's parent family (5.1) |
| 7:3 ret, -10, -1 | A | 13 (+46 continued) | 0.882-1.163 | 4,229-26,908 | 5,673-12,520 | - | excluded (asym) | **known-class: Casoliva 7-3b/c** (passes through the row) |
| 5:2 pro, 3, 0 | A | 3 (+84) | 2.411-2.580 | 495-8,748 | 11,772-16,348 | - | excluded (asym) | asymmetric, no catalogue row |
| 5:3 ret, -8, -2 | A | 4 (+78) | 0.670-0.942 | 200-14,548 | 17,876-22,444 | - | excluded (asym) | asymmetric, no catalogue row |
| 7:3 pro, 4, -1 | A | 2 (+34) | 2.414-2.560 | 267-6,373 | 10,716-14,064 | - | excluded (asym) | asymmetric, no catalogue row |
| 7:3 pro, 4, -1 | S | 50 | 2.432-2.468 | 11,007-12,865 | 5,900-7,515 | - | 0 (in J range; RR sparse here) | symmetric, unmatched |
| 5:2 pro, 3, -1 | S | 11 | 2.219-2.242 | 1,503-2,803 | 8,192-15,448 | - | 0 (in J range) | symmetric, unmatched |
| 5:2 ret, -7, 0 / -7, -1 | S | 26 + 12 | 1.009-1.469 | 605-32,712 | 8,854-39,527 | - | below J range | symmetric, unmatched |
| 2:3 ret, -5, -2 | S | 10 | -0.128-0.193 | 10,391-35,077 | 22,104-51,262 | - | below J range | symmetric, unmatched |
| 4:3 ret, -7, -2 | S | 4 | 0.313-0.351 | 30,814-34,447 | 57,253-59,209 | - | below J range | symmetric, unmatched |
| 7:3 ret, -10, 0; 7:3 ret, -9, 0; 8:3 ret, -11, 0 | S | 1; 1; 2 | 0.84-1.28 | 14,133-24,023 | 7,213-29,594 | - | below J range | symmetric, unmatched |

### 5.3 Verdict (nothing called novel; literature with `#972`)

1. **The positive controls pass for both methods.** casoliva-7-3b and 7-3c are recovered from
   skeletons. Their family, continued, passes through the row. Its pitchfork parent regenerates
   it.
2. **Asymmetric cycler-class families outside RR and FR by construction:** A (5:2 prograde),
   B (5:3 retrograde) and D (7:3 prograde). No catalogue row lies on them, and no symmetric
   pitchfork parent was reached within their continuation ranges.
   - They are near-Keplerian p:q resonant Earth-Moon cyclers. That CLASS is Casoliva et al.
     2010's Class 1, which includes the asymmetric 7-3b/c. Casoliva tabulates only p-q = 1-2,
     2-1, 3-2 and 7-3. So under spec sec. 16.4 these are at most computed members of a
     published class, at p:q values or senses the source did not tabulate.
   - That is my reading, not an adjudication. **`#972` decides.** All three are strongly
     unstable (|b| about 100-800). B and D reach the Earth floor at one end (C 0.9418 and
     2.4136). A stopped at the step cap in both directions; its lowest perigee is 495 km alt.
3. **Family C is the Casoliva 7-3b/c class:** known-class (the rows are on it).
4. **Symmetric candidates.**
   - F3 (Vaquero 2:1) and F4 3/7-r (Casoliva 7-3a) are known.
   - Nine other symmetric resonant families are unmatched by `#997` and RR: most lie below RR's
     J range (2.119), and RR is sparse below J 2.5. Symmetric, so in RR's class but not in its
     files. Also for `#972`.
5. **Coverage limits** (so that "not found" is not over-read):
   - The p:q set excludes 3:1, 7:2 and 4:1, because their Kepler skeleton cannot cross the lunar
     orbit. Yet a catalogued 3:1 cycler exists (vaquero-31-c254, perigee 792 km alt), held there
     by lunar perturbation. **The seed rule misses perturbation-enabled resonances.**
   - 54% of the grid solves impacted and 11% hit the iteration limit.
   - Continuations stopped mostly at the 40-step cap.
   - The grid is one 10 deg orientation sampling at a single (free) C per seed.

## 6. Checks after review, and amendment A3 (2026-10-08 14:10 AEDT)

- **A, B and D are asymmetric** (`scripts/run_1000_asym_check.py`,
  `data/1000_complement/asym_check.json`).
  - Every asymmetric grid member, and every 5th continuation member, was re-solved by multiple
    shooting. The minimum |xdot| at any y = 0 crossing was then taken ARC BY ARC from the
    nodes, not from one full-period integration (which drifts).

| family | min over members of min |xdot| at y = 0 |
|---|---|
| A | 0.0917-0.0919 |
| B | 0.150-0.174 |
| D | 0.012-0.038 |
| C (Casoliva 7-3b/c) | 0.0095 at its lowest grid member, next to the pitchfork; 0.08-0.26 elsewhere |

  - All values are 1e4 or more times the 1e-6 symmetry threshold, so the asymmetric labels
    stand.
  - One grid member of A (C 2.478) did not re-solve and is not counted.
  - The b values in sec. 4.1 come from these arc-STM products.
- **Amendment A3 (method (2) branching).**
  - The pre-registered step perturbed along the monodromy eigenvector at fixed C +- dC. It
    failed the control 36 of 36 times: every solve was ill-conditioned at the bifurcation, a
    known property of fixed-C shooting at a pitchfork. By the sec. 2.2 hard stop, **that
    variant is rejected.**
  - Replacement, introduced after the failure (13:30 AEDT): the symmetry-breaking parameter.
    Node 0 is placed at the parent's perpendicular crossing, with y(node 0) = 0 and xdot(node 0)
    = eta, and C free.
  - On the same control it regenerates the 7-3b/c branch (sec. 4.2). Because the control was
    used to develop it, that is not an independent validation.
  - Everything produced by it is labelled "eta variant": the control branch, `parents997.jsonl`,
    and sec. 5.1. **Ruling requested from the lead.** None of it changes the cycler-class
    verdict: A, B and D come from method (1), and method (2) found nothing cycler-class.
- **Evidence-only flag (no catalogue edit).** Casoliva 2010 Table 3 prints 7-3b and 7-3c with the
  SAME C (1.0687623900) and the SAME k (57.3519357), which suggests they are mirror images of
  each other. The catalogue's casoliva-7-3c row sits at C = 1.067197 (`#801`'s reproduction),
  and here its mirror is recovered 23 times against 1 for the row orientation.
- **Correction to commit e450d69d.** That commit's message describes these note edits, but they
  failed to apply (a text-anchor mismatch). It committed only the check script and its output.
  The edits are applied in the next commit.
