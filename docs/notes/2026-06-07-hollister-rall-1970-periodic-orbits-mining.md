# Rall (Hollister adv.) 1970 — Free-fall periodic orbits Earth-Mars (method mining)

Mined 2026-06-07 (Task #142). **NEVER MINED before.** The foundational
periodic-orbit search methodology. Directly informs our N-arc corrector
(`src/cyclerfinder/search/correct.py`), our family-continuation strategy, and the
**S1L1 multi-arc blocker** (see project memory `project_s1l1_realeph_closure_blocker`).

**Source (cite exactly, no file path):**
Rall, C. S., "Free-Fall Periodic Orbits Connecting Earth and Mars," Sc.D. thesis,
MIT Dept. of Aeronautics and Astronautics, report **TE-34** (also NASA CR /
Measurement Systems Laboratory, NASA Grant NGR 22-009-010), October 1969 (pub.
1970). Thesis advisor W. M. Hollister. 225 pp.

> Scanned typescript, clean. Vision read of Ch.1 front matter + Ch.2 (Numerical
> Techniques, pp.15-29), Ch.3.4-3.5 (symmetric/direct returns, pp.47-61), Ch.5
> results (pp.111-118) + Ch.6 summary/conclusions (pp.119-125). All numerals
> below read unambiguously.

---

## 1. The algorithm in 3 lines

Pose the periodic orbit as **N independent encounter dates `t`**; Lambert-solve
each leg between the resulting planet positions; the residual is the vector of
**hyperbolic-excess-SPEED differences `v(t)` at every flyby** (departure minus
arrival |V∞|); drive `v(t)=0` by **steepest-descent then Newton-Raphson** on a
numerically-differenced Jacobian `dv/dt`; close the cycle with `t_{N+1}=t_1+T_cycle`;
*after* convergence check flyby passing-distance (does the vehicle miss the
planet). Build the periodic orbit by **combining two reciprocal Earth-Mars-Earth
round trips with two symmetric series of "direct returns" at Earth**, enumerated
combinatorially. Mapping verdict: **this IS our N-arc magnitude-continuity
corrector, 55 years early** — same free vars, same residual, same Newton step,
same periodicity closure; we even share the "check bend/passing-distance
post-hoc" design.

---

## 2. The algorithm stated precisely

### 2.1 The basic numerical problem (Ch.2.1, pp.15-19)
> "One must first choose initial guesses for the N independent dates which
> determine the periodic orbit. These N dates then determine uniquely the
> positions of the encountered planets ... The trajectories between the encounter
> points are uniquely determined through the solution of Lambert's problem. The
> hyperbolic excess velocities are consequently determined ... One then has a
> difference in hyperbolic excess speed before and after each planetary encounter;
> this results in N speed differences. The remaining problem is to change the N
> dates so that the N speed differences can be made to approach zero." (p.15)

> "Only after the N speed differences are made very close to zero are the
> conditions during each flyby checked to see if the vehicle always misses the
> encountered planets." (p.15) — **bend/passing-distance is a POST-HOC check, not
> in the residual. Exactly our `residual_mode="magnitude"` default.**

Residual + Jacobian (Eqs.2-1..2-4, pp.16-17):
- `v(t) = 0` (Eq.2-1) — N speed-difference residuals, function of N dates.
- Two iteration schemes (verbatim, both "used in the computer program of
  Menning"):
  - **steepest descent** `t_new = t − ε [dv/dt]^T v(t)` (Eq.2-2)
  - **Newton-Raphson** `t_new = t − [dv/dt]^{-1} v(t)` (Eq.2-3)
- `ε` chosen per step so the merit `v^T v` decreases; step halved if needed.
- **Sparsity exploited** (Eq.2-4): `v_i = v_out,i(t_i,t_{i+1}) − v_in,i(t_{i-1},t_i)`
  — each speed-difference depends on only THREE dates (tridiagonal Jacobian), so
  one-sided differencing needs only `2N` extra trajectories, vs `N²` in the
  general case (p.17). **This is the structural reason a date-parameterised
  multi-flyby corrector is cheap — a sparsity insight we should exploit if we
  ever add an analytic/banded Jacobian to `correct.py`.**

Program loop (verbatim summary, p.18):
```
1. Read starting encounter dates.
2. Lambert-solve N or N+1 trajectories; form v = v_out − v_in at each encounter.
3. If not first pass: check v^T v reduced; if not, halve step, redo (2).
4. If Σ|v_i| < convergence requirement, skip to (8).
5. Form [dv/dt] by numerical differencing (2N more Lambert solves).
6. Newton-Raphson step if Σ|v_i| < preset value, else steepest-descent step.
7. Repeat (2)-(6) until convergence at (4).
8. Compute turn angle + passing distance at each encounter; print.
```
Lambert subroutine (p.18-19): planet positions+velocities → space triangle +
ToF → iterate semimajor axis → heliocentric velocities at ends → `v_out`,`v_in`
by differencing planet vs trajectory velocities.

### 2.2 Closing the cycle — "ends" of a periodic orbit (Ch.2.2, pp.19-21)
A real periodic orbit has no ends; the program enforces periodicity with
`t_{N+1} = t_1 + T_cycle` (Eq.2-5). Two ways to supply the Nth speed difference:
- **"A" modification** (Menning's): add an extra flyby/leg so N+2 encounters, and
  require `t_{N+1}=t_1+T_cycle`, `t_{N+2}=t_2+T_cycle` (Eq.2-6) — the (N+1)th leg
  equals the first leg in time.
- **"B" modification**: difference the departure |V∞| at `t_1` against the arrival
  |V∞| at `t_{N+1}` directly (one extra leg, N+1 encounters, N-1 flybys).
Both identical if the model is exactly periodic; differ for realistic models.
**Our `correct.py` pins one leg by the sourced period — this is the "A"/"B"
periodicity closure restated.**

### 2.3 Homotopy continuation circular→eccentric/inclined (Ch.2.4, pp.26-28)
THE key continuation idea, verbatim:
> "A method of dealing with this convergence difficulty is to increment the
> eccentricities and mutual inclinations of the planets' orbits in small steps
> while going from the circular coplanar case to the eccentric inclined case. One
> should start with the encounter dates for the circular coplanar case and use
> these dates ... with very small values for the eccentricities and mutual
> inclinations. Then, with convergence in this case, one should use these new
> numbers ... with slightly larger values ..." (pp.26-27)
> "Three or four proportional steps in eccentricities and inclination were found
> adequate for the worst cases." (p.27)
General principle (p.27): "one should not try to solve a numerical problem with
an answer 'too different' from the initial guess." Caveat: continuation does NOT
guarantee a solution exists — scheme M5-3 converged at e≈0.4×actual but had **no
solution above some e/inclination** (pp.27-28). **This is precisely the
seed-continuation ladder we want for taking circular-coplanar Russell rows onto
the real ephemeris — increment e,i in 3-4 steps, re-seed each from the previous
converged dates.** Maps to a strategy we have only partially built
(`search/seed_ladder.py`).

### 2.4 Solar-system model ladder (Ch.2.3, pp.22-24)
Three models, used as a fidelity ladder for convergence:
- **I. Circular coplanar** (e=i=0). I.A periodicity-forced `P_Mars=32/17 yr`,
  `a_Mars=(32/17)^{2/3} AU`; `P_Venus=8/13 yr`, `a_Venus=(8/13)^{2/3} AU`.
  I.B = correct a,P (includes alignment drift). (p.22)
- **II. Eccentric inclined, exactly periodic**: correct e,i but a,P forced to I.A
  so positions repeat **exactly after 32 years**. (p.23)
- **III. Eccentric inclined, constant mean elements** from a published series —
  the most realistic; sub-variants III.A ("A" mod) / III.B ("B" mod). (p.23)
Patched conic (Ch.1.2.1) used throughout, all models (p.24).

### 2.5 The building blocks — direct returns (Ch.3)
A "direct return" = a free-fall leg that returns the vehicle to the SAME planet
after the planet completes some revolutions. Taxonomy (List of Symbols + Ch.3):
- **FR** full-revolution return (integer planetary periods; in-plane).
- **HR** half-revolution return (occurs perpendicular to the orbit plane — the
  out-of-plane Z direction).
- **SiSR / LiSR** symmetric returns Shorter/Longer than `(i+0.45)` revolutions
  (linearised), `i` = revolutions of the Sun. Linear-case solution to
  `4(1−cos2πt) = 3πt sin2πt` (Eq.3-6) gives Table 3-1 (see §4).
- **Turn-angle selection** (Ch.3.2-3.3): for a series of returns, choose the V∞
  tip on the return locus to **minimize the maximum turn angle** across the
  series — "one test and one possible adjustment" (p.49): if the minimum-turn
  onto the FR locus `T_I` already ≥ the inter-FR angle `A_INT`, done; else move
  the contact point to equalise them.

### 2.6 Assembling the periodic orbit — the multi-arc construction (Ch.3.5, Ch.5, Ch.6)
THE construction, verbatim (p.121):
> "The approach used to obtain periodic orbits joining Earth and Mars was to
> combine two Earth-Mars-Earth round trips of Ross with two separate series of
> direct returns at Earth in a 'symmetric' manner. Use of the Earth-Mars-Earth
> round trips as segments of the periodic orbit schemes tried **avoided the
> difficulty of making direct returns at Mars**."
Why no Mars direct returns (p.121):
> "Obtaining periodic orbits to Mars was ... a more difficult problem ... because
> the small mass of Mars makes the necessary flyby maneuvers at Mars impossible,
> because the calculated trajectories intersect the surface of the planet."
Reciprocity (p.122): the two round-trip segments are "reciprocal" — encounter
dates are negatives of each other relative to opposition; centred on different
oppositions. Combinatorial enumeration: build the *list of all combinations* of
{HR, FR, S1SR, L1SR} up to 6-at-a-time / ≤3 synodic periods (Appendix C/D),
ordered by elapsed time, then prune (p.61, p.120, p.124).

---

## 3. Maps to our X / does not map

| Rall construct | Our code / concept | Verdict |
|---|---|---|
| Free vars = N encounter dates `t`; residual = per-flyby |V∞| speed-difference `v(t)=0` | `search/correct.py` free vars `[t0, leg ToFs]`, residual = flyby V∞-magnitude continuity | **MAPS — same problem statement.** Our ToF parameterisation is equivalent to his date parameterisation. |
| Steepest-descent → Newton-Raphson on `dv/dt` (Eqs.2-2,2-3) | `scipy.least_squares` (trust-region / LM) | **MAPS in spirit** — both gradient solvers driving `v^T v→0`; we use a library solver, he hand-rolled SD+NR. |
| Tridiagonal Jacobian (each `v_i` depends on 3 dates, Eq.2-4) | (we use dense numerical Jacobian) | **DOES NOT MAP — banded-Jacobian win available** for `correct.py`. |
| Periodicity `t_{N+1}=t_1+T_cycle` ("A"/"B" mods) | `correct.py` pins one leg by sourced period | **MAPS exactly.** |
| Post-hoc passing-distance/turn check after `v(t)=0` | `correct.py` bend feasibility checked post-hoc in `magnitude` mode | **MAPS exactly** — same design choice. |
| Homotopy continuation in (e, i), 3-4 steps, re-seed each step | `search/seed_ladder.py` (partial) | **PARTIAL MAP — directly actionable.** His circular→eccentric ladder is the recipe for putting circular-coplanar Russell rows on the real ephemeris. |
| **Multi-arc construction: 2 reciprocal E-M-E round trips + 2 series of Earth direct returns; NO Mars direct returns (Mars flyby hits surface)** | our S1L1 multi-arc finding (memory `project_s1l1_realeph_closure_blocker`) | **MAPS — and CONFIRMS the S1L1 blocker.** Rall states outright that single Mars-return arcs are infeasible (trajectory intersects Mars) — a periodic E-M orbit is intrinsically multi-arc (round trips stitched by Earth-side returns). This is independent 1969 corroboration that S1L1 never closing as a single ellipse is a *modelling-truth*, not an infra bug. |
| Turn-angle min-max selection over a return series (Ch.3.2-3.3) | `core/flyby.max_bend` / `flyby_bend_slacks` (per-flyby only) | **PARTIAL** — we check each bend; he *optimises the worst bend across a series* by sliding the V∞ contact point. A series-level min-max-bend objective we do not have. |
| Solar-system model ladder (circular→ecc-incl-periodic→ecc-incl-real) | our circular-coplanar vs analytic-ephemeris model_assumption split | **MAPS** — his three models ≈ our model_assumption tiers. |

---

## 4. Candidate test anchors (tabulated → golden-eligible)

Two clean tabulations whose EXPECTED side is source-traced (analytic/linear and
computed, but PUBLISHED — usable as candidate reference anchors, with the
standard "confirm before promoting to a hard golden" caveat since they are
patched-conic, circular-coplanar values).

**Anchor R1 — Table 3-1 (p.52): linear-case symmetric-return times + departure
angles.** Solutions of `4(1−cos2πt)=3πt sin2πt`. First 10 (time in planetary
periods | SiSR departure angle deg | LiSR departure angle deg):
1.4067 | 98.578 | −81.422; 2.4453 | 94.960 | −85.040; 3.4612 | 93.509 | −86.491;
4.4699 | 92.718 | −87.282; 5.4754 | 92.219 | −87.781; 6.4792 | 91.876 | −88.124;
7.4820 | 91.625 | −88.375; 8.4841 | 91.433 | −88.567; 9.4858 | 91.282 | −88.718;
10.4871 | 91.159 | −88.841.
> Closed-form-traceable (root of Eq.3-6) → the cleanest golden candidate in the
> paper; arrival angle = −departure angle (symmetry). EMOS = 29.77 km/s.

**Anchor R2 — Table 5-1 (pp.114-115): per-orbit statistics for M4-1, M5-1, M5-2.**
Circular-coplanar (CC) headline values per periodic orbit (V∞ in EMOS = ×29.77
km/s; passing distance in planet radii; turn angle deg):
- Earth encounters next to SHORT M-transfers — V∞_CC: M4-1 0.257, M5-1 0.249,
  M5-2 0.245; passing dist 1.54 / 1.78 / (1.42,2.06); turn 48.3° / 46.0° /
  (54.1,42.7)°.
- Mars encounters — V∞_CC: 0.314 / 0.316 / 0.314; passing dist (Mars radii) 3.77 /
  7.30 / 4.79; turn 4.3° / 2.3° / 3.4°. (Note tiny Mars turn angles — the surface
  -grazing problem.)
- Earth encounters next to LONG M-transfers — V∞_CC: 0.181 / 0.211 / 0.183;
  passing dist 1.30 / 1.55 / 1.37; turn 77.4° / 60.7° / 74.8°.
Eccentric-inclined adds avg/highest/lowest spreads + "change in encounter date
from CC" (M4-1 Earth-short: avg −0.8 d, RMS 17.7 d; Mars: RMS 28.2 d).
> CC V∞ at Mars ≈ 0.314 EMOS ≈ **9.35 km/s** — consistent with our memory that
> Mars-side V∞ is high (Russell rows ~20 km/s are a different, higher-energy
> family; Rall's are low-energy ~9 km/s at Mars). Earth V∞ 0.18-0.26 EMOS ≈
> **5.4-7.7 km/s** — matches the Jones-note S1L1 ">7 km/s max" cross-check.

**Conclusions hard numbers (p.125):** shortest period found = **4 synodic
periods** (M4 family); avg Earth speeds 0.260 & 0.181 EMOS; avg Mars speed 0.324
EMOS; avg min geocentric altitudes 1.40 & 3.60 Earth radii; Mars 3.04 Mars radii.
Abstract: most efficient orbit needs **4 E-M synodic periods (≈8.33 yr avg per 2
round trips)**; others need 5+. Circular-coplanar: ALL found orbits have exactly
**3 full-revolution returns (3FR) at Earth** between the two long E-M legs (p.111).

---

## 5. Single most implementable finding (this paper)

**The circular→eccentric/inclined homotopy continuation ladder (§2.3), re-seeding
each step from the previous converged dates.** This is the precise, field-tested
procedure for migrating our circular-coplanar Russell/Rall-class seeds onto the
real ephemeris — the exact gap flagged in the S1L1 blocker memory. Rall reports
3-4 proportional steps in (e,i) suffice for the hard cases, and warns the ladder
*surfaces non-existence* (M5-3 has no eccentric solution) rather than failing
silently — which is exactly the diagnostic behaviour we want from
`seed_ladder.py`. Secondary win: the **tridiagonal Jacobian** (Eq.2-4) — each
flyby residual depends only on 3 adjacent dates — would make `correct.py`'s
Jacobian banded (2N solves not N²).

---

## 6. v4.2 backfill checks

- **center**: all dynamics heliocentric (Sun-centred); no catalogue center
  ambiguity introduced.
- **tof_days_bounds**: the symmetric returns vary continuously in length (Table
  3-1 gives discrete linear-case anchors at 1.41, 2.45, 3.46, … planetary
  periods); short E-M transfer ToFs are NOT tabulated numerically here (only the
  return-leg durations). No new bounds to backfill.
- **source_ephemeris**: Rall uses **patched-conic on a periodic/constant-mean-
  element model** (his Models I/II/III; mean elements from a "published series",
  refs 7/8 — not a modern DE ephemeris). Any catalogue row tracing to Rall must
  carry `model_assumption: circular-coplanar` (Model I) or `analytic-ephemeris`
  (Models II/III, mean-element), `source_ephemeris: patched-conic mean-element
  (Rall TE-34, 1969)` — NOT a numerically-integrated ephemeris.

No existing catalogue rows are known to cite Rall directly; if M4-1/M5-1/M5-2 are
ever added they are **multi-arc, circular-coplanar (or mean-element), 4-6 synodic
period** Earth-Mars cyclers.

---

## 7. Honest "not extractable" list

- The full numeric encounter-date listings for M4-1/M5-1/M5-2 are in Appendices
  E/F (pp.179-204), NOT read in this pass (results focus was Table 5-1 + Ch.6).
  If golden encounter-date sequences are wanted, those appendices must be mined.
- Orbital elements (a,e,i) of the spacecraft arcs are not tabulated in the main
  text (the per-orbit signature is V∞/passing-distance/turn-angle, like Jones).
- Table 3-1 values are LINEAR-case limits (V∞→0); the finite-V∞ symmetric-return
  times/angles are given only as plots (Figs.3-3..3-7), not tabulated.
- "About 100 schemes attempted, 18 missed/grazed all planets, several converged"
  (p.123) — no exact success count.
