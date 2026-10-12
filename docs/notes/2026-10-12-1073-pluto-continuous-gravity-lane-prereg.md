# #1073: Pluto continuous-gravity lane, pre-registration

Status: PRE-REGISTRATION. Written and committed before any lane code or lane run. The build waits for the lead's
"go" and for the lead's ruling on the moon-node question of sec. 4.3. No catalogue writes. Nothing is called
novel.

Origin: the owner's decision on #1036, 2026-10-12, verbatim "(c) Hold for continuous gravity". No Pluto-Charon
rows are written until a continuous-gravity check answers the deep-dip caveat of
`docs/notes/2026-10-12-1036-owner-summary.md` sec. 4. Every strong member spends 5 to 50 percent of its Lambert-leg
time inside 2 a_C of the barycentre (r_min / a_C 0.018 to 0.83), where the one-centre patched conic is untested.

## 1. The question, and a finding that comes before it

The question: does each #1065 strong_1065 member (171) exist as an orbit in continuous Pluto-Charon gravity, and
does it stay the same object as its patched-conic member all the way to the full masses?

Before writing this note I measured one property of the 171 that decides how the question can be posed.
`scripts/analyse_1073_crank_planarity.py` produces `data/1073_pluto_cg/crank_planarity.json`. It is read-only on
the committed #1065 gauntlets and does no search and no solve.

- **166 of the 171 members are SPATIAL.** Each has one or two resonant Charon returns (`RCharon/n:m`), and the
  generator chooses the crank of each such return by its minimax-turn convention (`optimise_block`). That choice
  puts the V_inf out of the plane. The largest |v_z| / |V_inf| per member has a median of 0.79 (10th percentile
  0.56, maximum 0.92).
- **The other 5 are PLANAR.** They have no resonant leg, only Lambert legs. They are listed in sec. 8.1.
- **Planar-crank siblings.** Every resonant crank can instead be set to one of its two in-plane values. That
  changes no date and no V_inf magnitude, only the flyby turns. Of the 166 spatial members, none has a planar
  sibling that passes the #1065 strong filters: the gate, the Charon periapsis inside the SOI, the moving-Pluto
  surface and the unscheduled-pass screen. The best in-plane sibling's worst turn ratio has a median of 1.36. The
  minimax convention went out of the plane because the in-plane turns are too large.

The cranks matter because of the condition count in sec. 4.2. In the circular coplanar continuous model, with
the small moon as a pinned pass-through point (the definition the members were found under), a spatial member
has one more independent condition than unknowns. So a spatial member generically has no exact continuous-gravity
counterpart. This is a hypothesis, H1. It is tested first, in sec. 8.2, and it can be refuted. If it holds, the
deep-dip question can only be answered for the 5 planar members under the strict definition. The lead then has
to rule on the definition of the moon visit (sec. 4.3) before the other 166 can be judged at all.

## 2. Model

**Rung (a): the circular restricted three-body problem (CR3BP), spatial, Pluto and Charon massive.** This is the
continuous analogue of the #1024 barycentric patched conic, in the same way that #968/#1025 made the R-S ideal
model continuous. The plu060 ephemeris is a later rung (b) and is not part of #1073. Why rung (a) first:
- The owner's question is about the field inside the binary. The CR3BP isolates that field. In plu060 the
  differences from the circular model (small-moon e and i of order 1e-3 to 1e-2, epoch alignment) would be
  confounded with it.
- The #1065 members are zeros of the circular model. Their continuous counterparts must be sought in the matching
  continuous model first.
- The Jovian lane validated its rung (a) first. Its rung (b) has not closed a published control yet (#968 9.7,
  #1039, #1055), so a plu060 rung here would have no validated template.

Frame and units (production scale, as the project rule requires): barycentric inertial, x-y plane = Charon's
orbit plane, km, s, km/s, km^3/s^2.
- GM_total = 975.5 (`PRIMARIES["Pluto"]`, the #998 cell `mu`), GM_Charon = 106.1 (registry), GM_Pluto = 869.4
  (the difference); mass fraction mu_C = 0.108765.
- Charon's period P_C = 6.38722 d (the plu060-fitted value of #998). a_rel = (GM_total (P_C / 2 pi)^2)^(1/3) =
  19,596 km. Charon is on its 17,464.5 km barycentric circle and Pluto on the opposite 2,131.5 km circle, the
  #1024 geometry. With these values n^2 a_rel^3 = GM_total holds exactly, so at full mass the rails are the exact
  two-body solution of the binary and the model is the CR3BP, not an approximation of it.
- The small moon (Styx, Nix, Kerberos or Hydra) is massless, on its #998 circle at its fitted period: a kinematic
  waypoint, as in #998/#1065 and as Callisto in #968 rung (a). Its real GM (at most 0.003 km^3/s^2) is left out,
  because the patched conic passes through its centre.
- Radii: Charon 606 km (registry), floor 100 km (registry `safe_alt_km`). Pluto 1,188.3 km (Nimmo et al. 2017,
  the value in `scripts/gauntlet_998.py` and `run_1024_barycentric.py`). The #1036 summary's "1,194.3 km" is the
  smallest Pluto distance over the 171 (pnb 10), not a radius. The Pluto atmosphere is not modelled. A flag is
  raised below 1,288.3 km, the 100 km design floor of Stern et al. 2020 (corpus-sourced, #429).
- No J2, no solar tide, no small-moon gravity.

## 3. The sigma homotopy

The field at homotopy parameter sigma in [0, 1] is three point masses on rails that co-rotate with Charon:
- (1 - sigma) GM_total at the barycentre;
- sigma GM_Pluto on Pluto's circle;
- sigma GM_Charon on Charon's circle.

Properties, each pinned by a unit test at production scale before any run:
- At sigma = 0 the field is exactly the #1024 leg field (the 975.5 monopole at the barycentre). Charon's term
  vanishes, and with the periapsis offsets scaled by sigma (sec. 4.1) each Charon flyby tends to the patched-conic
  turn. This is the Jovian lane's limit (#968, #1025).
- At sigma = 1 the field is the CR3BP of sec. 2.
- The barycentre is the centre of mass at every sigma, and all three masses are fixed in Charon's rotating frame.
  So the field is time-independent in that frame at every sigma, and a Jacobi integral exists at every sigma.
  This is what makes the count of sec. 4.2 hold from the first step.
- The Charon point mass is softened at sigma x 606 km (the #1025 rule). The barycentre point mass gets a 1 km
  numerical clamp only, and the closest approach to the barycentre is recorded. Pluto's point mass is softened at
  1,188.3 km, unscaled, because Pluto is not a flyby node and distances to it do not scale with sigma.

Alternatives considered:
- Ramping Charon's mass fraction with Charon's radius following it (the #1024 ramp). Rejected: at sigma = 0 Charon
  would sit at 19,596 km, the old #998 model, not the #1065 zero.
- Direct Newton at sigma = 1 from the patched-conic seed. Kept as a recorded check only (#1025 sec. 3 item 6). A
  converged orbit found that way is not identified with the member.

## 4. Formulation

### 4.1 Closed chain (the #1025 / #968 corrector A method)

- Nodes in time order, over one period T = k S, where S is the Charon-moon synodic period.
- **Charon nodes:** Charon-relative position and velocity offsets plus the epoch (7 unknowns). Periapsis gauge
  (r - r_C) . (v - v_C) = 0 (1 row, weight 1e4). The seed comes from `periapsis_node` on the member's V_inf in and
  out (unclamped, #1025 deviation 3), with the offsets scaled by sigma.
- **Moon node (strict, M1):** the position is pinned to the moon at the node epoch. Unknowns: velocity and epoch
  (4).
- **Legs:** forward and backward half-arcs from consecutive nodes, matched at mid-time (6 rows, velocity rows
  weighted 1e3). A resonant leg is one arc across its full revolutions. The last leg closes on the first node
  rotated by Charon's advance over T (the Q wrap; Pluto and the moon return with the same rotation because T is a
  whole number of synodic periods).
- **Integration:** DOP853 with the analytic state transition matrix (STM), rtol 1e-13 and atol 1e-12, as #1025.
  The Jacobian is the STM plus the epoch columns. It is checked once against finite differences per member before
  use.
- **Floors:** 1e-3 km, 1e-6 km/s, gauge 1e-9. The #968 noise-floor rule applies to intermediate points only
  (1e-2 km, 1e-6 km/s). The sigma = 1 end point needs the floors.

### 4.2 The condition count, and hypothesis H1

With N_C Charon nodes and one pinned moon node:

| | unknowns | rows |
|---|---|---|
| total | 7 N_C + 4 | 6 (N_C + 1) + N_C = 7 N_C + 6 |
| Jacobi redundancy | | -1: the sum around the chain of J(forward end) - J(backward end) is zero for any unknowns, because J is constant on every arc and the Q wrap preserves J |
| effective | 7 N_C + 4 | 7 N_C + 5 |

The system is over-determined by one. Split at a planar orbit (z = v_z = 0 everywhere):
- **In-plane block:** 5 N_C + 3 unknowns, 5 N_C + 4 rows, minus the Jacobi redundancy, which is purely in-plane
  at a planar orbit. Square.
- **Out-of-plane block:** 2 N_C + 1 unknowns (z and v_z per Charon node, v_z at the moon node, whose z is pinned
  to 0), 2 (N_C + 1) rows. Over by one, and homogeneous at first order: z = 0 always solves it.

Consequences:
- A planar member is a square problem. Its continuation is the #1025 procedure unchanged.
- **H1: a spatial member has no exact solution in the strict model at any sigma > 0, generically.** Its
  out-of-plane block can only be satisfied on a codimension-one set.
- Cross-checks on earlier lane runs: #968 corrector A counted 18 unknowns, 20 rows and one redundancy, which is
  also over by one. It converged because GanCal#5 is planar. #1025 is square because two massive moons on rails
  of different periods have no Jacobi integral. That is why the square #1055 reduction found a root for gc4-19.
- Physical reading: an inclined resonant orbit returns to the plane at Charon only in the Kepler limit. In the
  binary field its line of nodes regresses. The Charon flyby offsets can absorb all but one of the resulting
  out-of-plane conditions.

**Test of H1** (the first spatial run after go; sec. 8.2). Take two spatial members among the 12 best. At sigma0
= 0.02, and again at 0.005, run the strict corrector to its stall. H1 is REFUTED if either converges to the
floors. H1 is SUPPORTED if both stall at a residual floor that does not fall with more iterations (40 Gauss-Newton
iterations, then 200 of trf), and whose largest component lies on the out-of-plane rows. The column-scaled
Jacobian must also show 7 N_C + 4 nonzero singular values against 7 N_C + 6 rows, with the left null space
(dimension 2) spanned by the Jacobi direction and an out-of-plane direction. The residual floor and its split
into in-plane and out-of-plane parts are recorded either way.

### 4.3 The moon visit: a lead decision before the 166 spatial members are run

- **M1, strict pass-through.** This is the definition the members were found under (#998, #1065: the craft
  passes through the massless moon's centre). Under H1 only planar members can exist. The 166 spatial members
  then end as "no continuous counterpart under the strict definition". None has a planar sibling that passes the
  strong filters (sec. 1).
- **M2, the moon visit as a close approach.** The moon node becomes a closest approach. The offset from the moon
  lies in the plane normal to the relative velocity, its in-plane component stays 0, and its out-of-plane
  component b_z is a new unknown. That makes the system square (7 N_C + 5 unknowns). Accept only if |b_z| at
  sigma = 1 is at most d_max. **d_max must be stated before any M2 run.** Candidates are the moon's Hill radius
  (about 150-670 km for Styx to Hydra) or a fixed encounter distance named by the owner. M2 changes the object
  from "passes through the moon" to "passes within d_max of the moon each cycle", so it is the owner's call
  through the lead.
- **M2 as a diagnostic under M1.** The converged b_z(sigma) of a spatial member says how far the member is from
  existing under M1. A value that grows with sigma means the miss is set by the field, not by the seed. It can be
  reported under either ruling.

Under M2 the system is square, and the #1055 degenerate-family method applies as written (sec. 5). Under M1 it
does not port: see sec. 5.

## 5. The degenerate full-revolution family (#1055)

- Each resonant leg (`RCharon/n:m`) leaves its crank free in the patched conic. The structure is then a q-parameter
  family, q = the number of resonant legs (2 for most of the 166; `n_resonant_legs` in the planarity file).
  Continuous gravity selects members by a solvability condition at first order in sigma (#1055 sec. 1).
- The minimax crank is a convention of the generator. Nothing ties it to a root of the solvability condition.
  GanCal#1 showed that a generator's member can be a non-root (#1055 sec. 12.7). In a coplanar model, the planar
  cranks are roots by mirror symmetry (z to -z sends each crank to its mirror). That is why gc-1's forcing
  vanished in the ideal model (#1055 sec. 1, item 3).
- **Under M1:** the reduced system has q unknowns (the cranks) against q + 1 conditions, the extra one being H1's
  row. `degenerate_family.py` solves the square bordered system, so it does not port as written. No family stage
  is run under M1. The planar members have no resonant leg: the detection of #1055 sec. 2 is run on them and q = 0
  is expected. A Charon-to-Charon Lambert leg is not degenerate unless its flight time is a whole number of Charon
  periods. Any q > 0 found there is reported, and that member stops for the lead.
- **Under M2:** the #1055 method, imported from `src/cyclerfinder/nbody/degenerate_family.py` and not edited.
  Detection at sigma0 and sigma0 / 4, the family scan of g(alpha) both ways up to 3 Newton-step lengths, the
  root, the hand-off. Then amendment 10 (member identity, sec. 6.2) decides whether the closure is the generator's
  member or a sibling with other cranks. If a root's cranks differ from the minimax cranks, the member is reported
  as a SIBLING (sec. 6.3). A sibling's own patched-conic gate is evaluated at its cranks and recorded.

## 6. Per-member procedure and acceptance (fixed now)

### 6.1 Steps

1. Build: the member's dates and minimax flybys (`cycle_flybys`), from the gauntlet record addressed by cell and
   index.
2. Start at sigma0 = 0.02. If that fails, try 0.005, then 0.00125, each from the scaled seed (Charon's effective
   mass fraction at 0.02 is about 2e-3, against 2e-6 for Ganymede in #1025, so a smaller start may be needed). If
   none converges: NO START. Under M1 for a spatial member this is H1's predicted outcome, and it is reported with
   the residual split of sec. 4.2.
3. Identity ladder going down from the start: halve sigma until 0.00125, each point from the previous one with the
   offsets scaled (sec. 6.2).
4. Natural continuation up to sigma = 1: the #1025 procedure (secant in log sigma, factor 1.2, halving to 1e-5,
   monitor-variable continuation at a stop, the two-solution fold check, direct Newton at sigma = 1 recorded
   only).
5. Impact ends a member at that sigma:
   - a Charon periapsis at or inside sigma x 606 km;
   - any point of a converged orbit within 1,188.3 km of Pluto's centre, on a 0.002-d scan of the arcs plus local
     minimisation. Pluto's radius is not scaled, because a Pluto pass is not a node.

   Flags only: the Charon floor sigma x (606 + 100) km, and Pluto below 1,288.3 km.
6. At sigma = 1, three end checks:
   - IAS15 re-fly of every half-arc (REBOUND, epsilon 1e-11, the same rails). PASS below 1e-2 km and 1e-7 km/s.
   - **Live-body re-fly:** REBOUND with Pluto and Charon as live massive particles started on their circular
     states, not rails, and the craft as a test particle, re-flying every half-arc. Same tolerance. This checks
     the rails.
   - **The #1025 encounter scan:** off-node distance minima to Charon inside its 7,192 km SOI, and to Pluto,
     recorded.

### 6.2 Identity (amendment 10, adapted to Pluto)

The patched-conic reference is the generator's member: dates, minimax cranks, and the Charon r_p and turn from
`cycle_flybys`. Tolerances are scaled from the Jovian ones by the ratio of the speeds and time scales:
- **Node V_inf:** within 0.4 percent of the patched-conic value. Amendment 10's 0.01 km/s is 0.4 percent of a gc
  member's V_inf. Here that is 0.7-1.4 m/s at Charon and 0.3-0.7 m/s at the moon.
- **r_p / sigma:** within 5 percent (dimensionless, unchanged).
- **Turn:** within 0.2 deg or 5 percent, whichever is larger (unchanged).
- **Node epoch:** within 2.5 percent of the member's median Charon flyby time scale r_p / V_inf. That is about
  400 s for a 4,000 km, 0.25 km/s flyby. The same rule gives 60 s for a gc flyby (2,500 s), which is amendment
  10's value.

Fit rule, fixed now (#968 sec. 8.1):
- The verdict quantity is the first-order (Richardson) limit from the two smallest converged sigma, s and s / 2:
  q* = 2 q(s / 2) - q(s).
- PASS needs all three:
  - q* within tolerance for every node quantity;
  - the last computed point within twice the tolerance;
  - for the V_inf shifts, the ratio |dq(s)| / |dq(s / 2)| in [1.4, 2.8] (first-order behaviour), unless both
    shifts are below 0.1 of the tolerance.
- The smallest sigma must be 0.0025 or less.
- FAIL on V_inf or epoch: IDENTITY FAIL, and the run stops.
- FAIL on r_p or turn only, with V_inf and epoch passing: the closure is a different member of the same structure
  (amendment 10). It continues and is labelled SIBLING.

### 6.3 Outcomes

| Outcome | Meaning |
|---|---|
| EXISTS | identity PASS, sigma = 1 reached, no impact, IAS15 and the live-body re-fly pass |
| EXISTS-SIBLING | as EXISTS, but on a family member with other cranks or periapses (amendment 10 FAIL on r_p or turn); not the generator's object |
| FOLDS at sigma_f | the #1025 fold checks are supported |
| IMPACT at sigma_i | Charon or Pluto, with which one named |
| NO START | no convergence at any start sigma; under M1 for a spatial member this is H1's outcome |
| IDENTITY FAIL | |
| NUMERICAL STOP (reason) | |

Only EXISTS means "the member exists in the continuous ideal model". EXISTS-SIBLING is reported with the sibling's
own patched-conic gate at its cranks.

## 7. Controls (before any member is judged)

### 7.1 Positive control PC1: published periodic orbits at the Pluto mass ratio

- Source: Bruno and Varin 2009, Sol. Syst. Res. 43(1):2-25, Table 4 (mu = 0.1, 28 family-h orbits), transcribed in
  `data/sources/bruno-varin-2009-family-h-small-mu-tables.yaml`. The digest is
  `docs/notes/2026-10-07-digest-bruno-varin-2009-family-h-small-mu-ssr.md`.
- Family h passes P2 (the Charon analogue) at T / 2, so these are published deep-secondary-pass periodic orbits
  at mu = 0.1. That is within 9 percent of Pluto-Charon's 0.1088.
- Conventions, from the digest: the frame origin is at P1 with P2 at (+1, 0), and y1, y2 are canonical momenta.
  The printed period is T~ = (1 - mu) T / (2 pi). The printed C is -2H + mu.
- Rows, by a rule fixed now, from the 12 rows the digest integrated (1-4, 12-14, 20, 21, 24-26):
  - (a) the row with x1(T/2) closest to P2, among those whose digest check is within 5e-5 (this excludes the
    y2 = 6.2758 row the digest flags);
  - (b) the row with the largest T~;
  - (c) row 1.
- The orbit is run in the lane's code at Pluto scale, with mu = 0.1 exactly (the paper's), GM_total = 975.5 and
  a_rel = 19,596 km. It is posed as a closed chain of the same corrector, with nodes at t = 0 and T / 2 and the Q
  wrap by n T.
- PASS needs every one of these:
  - convergence to the floors from the printed state;
  - y2(0) within 2e-5 of the printed value at the printed x1(0);
  - (1 - mu) T / (2 pi) within 3e-5 of T~;
  - x1(T/2) and y2(T/2) within 5e-5;
  - -2H + mu within 5e-5 of C;
  - the IAS15 re-fly passing;
  - J conserved along the orbit to 1e-10 relative.

  The nondimensional tolerances are the digest's own integration levels for 5-decimal printing. 5e-5 is 1 km and
  1.1e-5 km/s at this scale.

### 7.2 Positive control PC2: a published one-working-body cycler with a known continuous-gravity outcome

- GanCal#5 (Russell and Strange 2009). It is the same architecture as the Pluto members: one working body
  (Ganymede) and a massless pass-through target (Callisto).
- It runs through the Pluto lane code in the Jupiter-Ganymede binary: `barycentric_binary` of the R-S gc cell
  with Ganymede as the secondary, mass fraction 7.80e-5. Its seed is its own barycentric patched-conic zero (the
  #1024 on-path control route).
- PASS needs both:
  - identity by the sec. 6.2 rule against that zero;
  - sigma = 1 within 10 mu v_G = 8.5e-3 km/s of #968's full-mass orbit (V_inf G 3.0403 km/s, Callisto speed
    3.2513 km/s), with r_p within 1 percent of 3,210.5 km. 10 mu v_G is the #1024 tolerance for the O(mu) frame
    difference: #968 is Jupiter-centred, with Callisto's circle about Jupiter.
- What it does and does not test:
  - It tests the homotopy, the nodes, the moon pin, the wrap, the identity ladder and the continuation, against
    an independent code path (#968 corrector A) on a published object.
  - It is planar with no resonant leg, so it cannot test the spatial question.
  - PC1 tests the field at the Pluto mass ratio.
  - No published spatial periodic orbit at mu near 0.1 is in the corpus. The spatial code is covered by unit
    tests: the 3-D analytic STM against finite differences, and the Jacobi integral conserved in 3-D at sigma = 1
    and at sigma = 0.5. It is also covered by the H1 test itself.

### 7.3 Negative controls

- **N1, the gate.** The #1065 gate-failing zero (from the archived zeros, `~/dev/references/cyclers-runs/1065/`)
  of the pnb cell with the smallest worst ratio above 1.
  - Expected: no EXISTS. In the limit its demanded Charon periapsis is below the surface, so it should end as
    IMPACT at or near the start.
  - An EXISTS fails the control and stops the lane.
- **N2, identity discrimination.** phb 26 against its sibling zero from #1065 sec. 10.2: same key, the same V_inf
  to 1e-15 km/s, dates 2.46 d apart, gate 1.287.
  - The sec. 6.2 identity check, run on phb 26's ladder against the sibling's patched-conic member, must FAIL on
    the epochs.
  - Cheap, and it guards against the key + V_inf matching artefact of #1065.
- **N3 (M2 only), the #1055 mechanism control.** Starts from the family 0.5 L either side of a root must fail
  (#1055 sec. 3, item 3).

A failed control stops the work. It is reported with a diagnosis, and no member is judged.

## 8. Order of work after go

1. Code with tests first (TDD), in `src/cyclerfinder/nbody/pluto_cg.py` (field, STM, chain residual, Jacobian)
   and a driver `scripts/run_1073_pluto_cg.py` (preflight, task 1073). The tests are at production scale (km,
   km^3/s^2) and pin:
   - the sigma = 0 field equal to the 975.5 monopole;
   - the sigma = 1 field equal to the CR3BP;
   - the analytic STM against finite differences (3-D);
   - the Jacobi integral conserved at sigma = 0, 0.5 and 1;
   - the Q wrap preserving J;
   - the rails equal to `barycentric_binary`'s states.
2. PC1, PC2, N1, N2.
3. **The 5 planar members (sec. 8.1), strict M1.** These are the clean test of the deep-dip question, with no
   degeneracy and no H1 issue. pnb 10 is one of the 12 best.
4. **The H1 test (sec. 4.2)** on two spatial members among the 12 best (sec. 8.2).
5. The rest under the lead's ruling on sec. 4.3.

### 8.1 The 5 planar members (strong_1065, no resonant leg)

| Cell, index | Key | Worst ratio | Smallest Pluto distance (km) |
|---|---|---|---|
| psb 130 | k6, LCharon>Charon/2h, LCharon>Styx/1l, LStyx>Charon/0s | 0.253 | 1,383 |
| pnb 10 | k4, LCharon>Charon/1l, LCharon>Nix/0s, LNix>Charon/0s (one of the 12 best) | 0.357 | 1,194 |
| pnb 14 | k4, LCharon>Charon/2l, LCharon>Nix/1l, LNix>Charon/0s | 0.524 | 1,605 |
| pnb 34 | k5, LCharon>Charon/2h, LCharon>Nix/0s, LNix>Charon/0s | 0.275 | 2,756 |
| phb 31 | k6, LCharon>Charon/2l, LCharon>Hydra/0s, LHydra>Charon/0s | 0.257 | 1,575 |

All five pass within about 1,400 km of Pluto's surface. pnb 10 passes 6 km above it in the patched conic.
Continuous Pluto gravity acts hardest exactly there, so IMPACT is a live outcome.

### 8.2 The 12 best, and which are run when

- The 12 best of #998 sec. 6 are identified in #1065 by DATES, not by key plus V_inf. Each #1024 continued record
  (`data/1024_barycentric/<cell>_continue.json`, by #998 gauntlet index) is matched to the #1065 merged group that
  contains a zero equal to it under `_same` (dates modulo the phase period), using the archived zeros. A
  representative's key can differ from the member's: a key lookup found only 7 of the 12, and key plus V_inf is
  what produced the two #1065 artefacts.
- #1036 reports that 9 of the 12 are among the 171. By key, 11 of the 12 have two resonant legs. The exception is
  Nix k4 with the Charon-to-Charon Lambert leg, which is pnb 10, planar.
- For the H1 test: the two spatial best members with the lowest worst ratio among those found in the 171.
- The remaining best members run after the lead's sec. 4.3 ruling. Under M1 they are expected to end as NO START.

## 9. Cost

- Measured 2026-10-12 13:53 AEDT at load average 35: one DOP853 state-plus-STM propagation of 34.4 d in the
  Pluto-Charon CR3BP (rtol 1e-13; it passed 17 km from Pluto) took 0.86 s and 24,350 right-hand-side evaluations.
  A residual plus Jacobian covers one period in half-arcs, so 0.3-1.5 s, the same order as the Jovian 0.6-0.7 s.
- Per member:
  - start and identity ladder: about 40 evaluations;
  - continuation: 20-40 steps of 3-6 evaluations;
  - fold checks: about 50;
  - end checks: about 2 min.

  Total about 250-500 evaluations, 2-12 min. Under M2, the #1055 family stage adds 50-150 evaluations.
- Controls: under 10 min in all.
- The 5 planar members: under an hour.
- The H1 test: about 10 min.
- All 171 under M2: about 6-35 h serial, at the current load.
- A pilot at the actual load sets the shard size before any launch. Any call longer than about 10 minutes goes to
  the lead as an exact command, with the expected duration and the output path. Runs are checkpointed per member
  and per sigma point, with per-stage timestamps and running counts. Output goes to `data/1073_pluto_cg/`, live
  logs to `data/1073_pluto_cg/live/`.

## 10. Expectations (fixed now)

| Item | Expectation |
|---|---|
| PC1 passes | 0.9 |
| PC2 passes | 0.8 |
| N1 and N2 behave as required | 0.95 |
| H1 supported by the test | 0.85 |
| Each planar member EXISTS | about 0.3 each. Charon is 11 percent of the mass, V_inf at Charon is 0.17-0.36 km/s, and every one of the five dives within 1,600 km of Pluto's centre (pnb 34: 2,756 km), so folds and Pluto impacts are likely |
| Each spatial member under M2 | EXISTS or EXISTS-SIBLING about 0.3 each, with sibling likelier than exact identity. No probability is given for the share within d_max until d_max is fixed |

## 11. What would change the owner's picture

- **EXISTS:** the member carries continuous-gravity evidence in the ideal model, as gc-1 does under #1034. The
  deep-dip caveat is answered for that member in the circular model. plu060 continuous gravity stays a separate
  rung.
- **EXISTS-SIBLING:** the patched-conic member as written does not continue. A neighbour with other cranks or
  periapses does. A row could only describe the sibling, with its own gate result. That is the owner's call.
- **FOLDS, IMPACT, NO START:** no continuous counterpart on this branch. Existence elsewhere is untested. The member
  is lost for option (a)/(b) rows.
- **If H1 holds and the lead keeps M1:** at most the 5 planar members can carry continuous-gravity evidence. 166 of
  the 171 have no exact counterpart, by a structural count, not by a numerical failure. That is a finding about
  the whole set: the generator's minimax crank produces orbits that need the patched conic's decoupling to
  exist.
- **If H1 is refuted:** the spatial members are judged one by one under M1, and sec. 4.3 is moot.

## 12. Out of scope

- The 48 half-rev members: they have no barycentric zero.
- The 17 #1065 seed-grid misses: they can be added later by the same procedure, from `data/1024_barycentric/`.
- The plu060 rung, Pluto's atmosphere, J2, and small-moon gravity.

Code reuse: `degenerate_family.py`, `fold_tools.py` and the #1025 driver's functions are imported, not edited (the
Jovian lane files are owned elsewhere).
